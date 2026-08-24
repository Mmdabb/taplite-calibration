"""Mode-aware workflow orchestration and auditable run manifests."""

from __future__ import annotations

import hashlib
import json
import logging
import platform
import re
import traceback
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

from . import __version__
from .auto_calibration import run_auto_calibration
from .config import ProjectConfig, load_config
from .odme import run_odme
from .paths import portable_path
from .prepare import ensure_prepared_inputs


LOGGER = logging.getLogger("taplite_calibration.pipeline")
RUN_ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,79}$")


@dataclass(frozen=True)
class RunResult:
    run_dir: Path
    manifest_path: Path
    manifest: Dict[str, object]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def _default_run_id(mode: str) -> str:
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return "{}-{}".format(timestamp, mode.replace("auto-calibration", "auto"))


def _write_manifest(path: Path, manifest: Dict[str, object]) -> None:
    path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def run_loaded_config(config: ProjectConfig, run_id: Optional[str] = None) -> RunResult:
    run_id = run_id or _default_run_id(config.mode)
    if not RUN_ID_PATTERN.fullmatch(run_id):
        raise ValueError(
            "run_id must start with a letter or digit and contain at most 80 "
            "letters, digits, dots, underscores, or hyphens"
        )
    run_dir = config.output_dir / "runs" / run_id
    run_dir.mkdir(parents=True, exist_ok=False)
    log_dir = run_dir / "logs"
    log_dir.mkdir()
    log_path = log_dir / "run.log"
    file_handler = logging.FileHandler(log_path, encoding="utf-8")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(
        logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s")
    )
    package_logger = logging.getLogger("taplite_calibration")
    if package_logger.level == logging.NOTSET:
        package_logger.setLevel(logging.INFO)
    package_logger.addHandler(file_handler)
    manifest_path = run_dir / "run_manifest.json"
    if config.mode == "both":
        stage_order = ["prepare", "odme", "auto-calibration"]
    else:
        stage_order = ["prepare", config.mode]
    manifest: Dict[str, object] = {
        "schema_version": 1,
        "status": "running",
        "package": "taplite-calibration",
        "package_version": __version__,
        "run_id": run_id,
        "mode": config.mode,
        "stage_order": stage_order,
        "processors": config.processors,
        "project_config": portable_path(config.config_path, config.project_dir),
        "project_config_sha256": _sha256(config.config_path),
        "runtime": {
            "python": platform.python_version(),
            "implementation": platform.python_implementation(),
            "platform": platform.platform(),
        },
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "stages": [],
        "artifacts": {},
        "log": portable_path(log_path, config.project_dir),
    }
    _write_manifest(manifest_path, manifest)
    artifact_dir = run_dir / "artifacts"
    artifact_dir.mkdir()
    stages: List[Dict[str, object]] = []

    try:
        preparation = ensure_prepared_inputs(config, run_dir)
        config = preparation.config
        stages.append(preparation.manifest)
        manifest["artifacts"]["preparation_manifest"] = portable_path(  # type: ignore
            run_dir / "00-prepare" / "manifest.json", config.project_dir
        )
        scenario_for_auto = (
            config.auto_calibration.scenario_root
            if config.auto_calibration is not None
            else None
        )
        if config.mode in ("odme", "both"):
            assert config.odme is not None
            stage_root = run_dir / "01-odme"
            result_root = stage_root / "results"
            scenario_root = stage_root / "adjusted-scenarios"
            result = run_odme(
                config.odme,
                result_root,
                scenario_root,
                config.project_dir,
            )
            scenario_for_auto = scenario_root
            stages.append(result)
            manifest["artifacts"]["od_factor_dictionary"] = portable_path(  # type: ignore
                result_root / "od_factor_dictionary.npy", config.project_dir
            )
            manifest["artifacts"]["odme_screen_comparison"] = portable_path(  # type: ignore
                result_root / "screen_joint_daily_fixed_policy.csv", config.project_dir
            )

        if config.mode in ("auto-calibration", "both"):
            assert config.auto_calibration is not None
            assert scenario_for_auto is not None
            stage_number = "02" if config.mode == "both" else "01"
            stage_root = run_dir / "{}-auto-calibration".format(stage_number)
            calibration_root = stage_root / "period-runs"
            finalized_root = stage_root / "finalized"
            dictionary_output = artifact_dir / "calibrated_qvdf_node_pair_dict.npy"
            result = run_auto_calibration(
                config.auto_calibration,
                scenario_for_auto,
                calibration_root,
                finalized_root,
                dictionary_output,
                config.processors,
                config.project_dir,
            )
            stages.append(result)
            manifest["artifacts"]["calibrated_qvdf_dictionary"] = portable_path(  # type: ignore
                dictionary_output, config.project_dir
            )
            for period in ("am", "md", "pm"):
                manifest["artifacts"]["{}_calibrated_link".format(period)] = portable_path(  # type: ignore
                    finalized_root / "assignment" / period / "link_calibrated.csv",
                    config.project_dir,
                )

        artifact_index = {
            "status": "complete",
            "run_id": run_id,
            "mode": config.mode,
            "artifacts": manifest["artifacts"],
        }
        (artifact_dir / "artifact_index.json").write_text(
            json.dumps(artifact_index, indent=2) + "\n", encoding="utf-8"
        )
        manifest["stages"] = stages
        manifest["status"] = "complete"
        manifest["completed_utc"] = datetime.now(timezone.utc).isoformat()
        manifest["artifact_index"] = portable_path(
            artifact_dir / "artifact_index.json", config.project_dir
        )
        _write_manifest(manifest_path, manifest)
        LOGGER.info("Run %s completed", run_id)
        file_handler.flush()
        package_logger.removeHandler(file_handler)
        file_handler.close()
        return RunResult(run_dir, manifest_path, manifest)
    except BaseException as error:
        manifest["stages"] = stages
        manifest["status"] = "failed"
        manifest["failed_utc"] = datetime.now(timezone.utc).isoformat()
        manifest["failure"] = {
            "type": type(error).__name__,
            "message": str(error),
            "traceback_file": portable_path(run_dir / "failure_traceback.txt", config.project_dir),
        }
        (run_dir / "failure_traceback.txt").write_text(
            traceback.format_exc(), encoding="utf-8"
        )
        _write_manifest(manifest_path, manifest)
        file_handler.flush()
        package_logger.removeHandler(file_handler)
        file_handler.close()
        raise


def run_project(
    project_dir: Path,
    config_path: Optional[Path] = None,
    mode: Optional[str] = None,
    processors: Optional[int] = None,
    run_id: Optional[str] = None,
) -> RunResult:
    """Public Python API for ODME, auto calibration, or the ordered pipeline."""
    config = load_config(
        project_dir,
        config_path=config_path,
        mode_override=mode,
        processors_override=processors,
    )
    return run_loaded_config(config, run_id=run_id)
