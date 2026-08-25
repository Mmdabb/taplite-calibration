"""TAPLite refined auto-calibration staging and postprocessing."""

from __future__ import annotations

import csv
import hashlib
import json
import logging
import math
import os
import shutil
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Set, Tuple

import numpy as np

from .config import AutoCalibrationConfig, CALIBRATION_PERIODS
from . import native as taplite_native
from .paths import portable_path


LOGGER = logging.getLogger("taplite_calibration.auto_calibration")
MUTABLE_INPUTS = {
    "link.csv",
    "settings.csv",
    "auto_calibration_settings.csv",
    "departure_profiles.csv",
    "volume_constraints.csv",
}
EXCLUDED_OUTPUTS = {
    "link_performance.csv",
    "route_assignment.csv",
    "agent.csv",
    "summary_log_file.txt",
    "final_summary.csv",
    "route_columns.bin",
    "destination_accessibility.csv",
    "origin_accessibility.csv",
    "od_performance.csv",
    "system_performance.csv",
    "inaccessible_od.csv",
    "dtalite_run.log",
    "TAP_log.csv",
    "RUN_FAILURE.json",
    "RUN_SUMMARY.md",
    "RUN_CARD.md",
}
PARAMETER_COLUMNS = {
    "vdf_plf": "final_plf",
    "vdf_cd": "final_qcd",
    "vdf_cp": "final_qcp",
    "vdf_n": "qn",
    "vdf_s": "qs",
    "vdf_alpha": "alpha",
    "vdf_beta": "beta",
    "mode1_plf": "final_plf",
    "mode1_qcd": "final_qcd",
    "mode1_qcp": "final_qcp",
}
DICTIONARY_PARAMETERS = {
    "plf": "vdf_plf",
    "qdf": "vdf_qdf",
    "n": "vdf_n",
    "s": "vdf_s",
    "cp": "vdf_cp",
    "cd": "vdf_cd",
    "alpha": "vdf_alpha",
    "beta": "vdf_beta",
}
PERIOD_SEQUENCE = {"am": 1, "md": 2, "pm": 3}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def _finite_positive(value: object) -> bool:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return False
    return bool(np.isfinite(number) and number > 0.0)


def _text(value: object) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    return str(value)


def stage_scenario(source: Path, destination: Path) -> Dict[str, int]:
    """Stage immutable inputs by hard link where supported, copying mutable CSVs."""
    if destination.exists():
        raise FileExistsError(destination)
    destination.mkdir(parents=True)
    copied = 0
    linked = 0
    for path in source.iterdir():
        if (
            not path.is_file()
            or path.name in EXCLUDED_OUTPUTS
            or path.name.startswith("auto_calibration_")
            and path.name not in MUTABLE_INPUTS
        ):
            continue
        target = destination / path.name
        if path.name in MUTABLE_INPUTS:
            shutil.copy2(path, target)
            copied += 1
        else:
            try:
                os.link(path, target)
                linked += 1
            except OSError:
                shutil.copy2(path, target)
                copied += 1
    return {"copied_inputs": copied, "hard_linked_inputs": linked}


def validate_qvdf_profile(link_path: Path) -> Dict[str, object]:
    """Enforce the production QVDF/profile contract before a costly solve."""
    with link_path.open("r", newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        fields = set(reader.fieldnames or ())
        required = {"vdf_type", "qvdf_profile_mode"}
        missing = sorted(required - fields)
        if missing:
            raise ValueError("{} is missing columns: {}".format(link_path, missing))
        rows = 0
        bad_vdf = 0
        bad_profile = 0
        for row in reader:
            rows += 1
            try:
                bad_vdf += int(float(row["vdf_type"]) != 2.0)
                bad_profile += int(float(row["qvdf_profile_mode"]) != 1.0)
            except (TypeError, ValueError):
                bad_vdf += 1
                bad_profile += 1
    if bad_vdf or bad_profile:
        raise ValueError(
            "{} violates the required vdf_type=2/qvdf_profile_mode=1 contract "
            "(bad vdf_type={}, bad profile_mode={})".format(
                link_path, bad_vdf, bad_profile
            )
        )
    return {"links": rows, "vdf_type": 2, "qvdf_profile_mode": 1}


def _run_native_period(
    config: AutoCalibrationConfig,
    source: Path,
    destination: Path,
    period: str,
    processors: int,
    project_dir: Path,
) -> Dict[str, object]:
    staging = stage_scenario(source, destination)
    settings_source = config.calibration_settings_csv or (
        source / "auto_calibration_settings.csv"
    )
    shutil.copy2(settings_source, destination / "auto_calibration_settings.csv")
    profile_contract = validate_qvdf_profile(destination / "link.csv")
    started = time.time()
    result = taplite_native.auto_calibrate(
        destination,
        in_place=True,
        isolated=True,
        timeout=config.timeout_seconds,
        settings_overrides={
            "number_of_processors": processors,
            "auto_calibration": 1,
            "column_output": 2,
            "column_file_output": 0,
            "route_output": 0,
            "vehicle_output": 0,
            "link_output": 1,
            "accessibility_output": 0,
        },
    )
    required_outputs = (
        "link_performance.csv",
        "auto_calibration_history.csv",
        "auto_calibration_link_audit.csv",
        "auto_calibration_summary.json",
    )
    missing = [name for name in required_outputs if not (destination / name).is_file()]
    if result.returncode != 0 or missing:
        raise RuntimeError(
            "{} calibration failed: rc={}, missing={}".format(
                period.upper(), result.returncode, missing
            )
        )
    manifest: Dict[str, object] = {
        "status": "complete",
        "period": period.upper(),
        "backend": "native",
        "source": portable_path(source, project_dir),
        "destination": portable_path(destination, project_dir),
        "processors": processors,
        "kernel": "bundled taplite_calibration._native",
        "profile_contract": profile_contract,
        "staging": staging,
        "elapsed_seconds": time.time() - started,
        "intermediate_output_policy": (
            "column_output=2 in memory; column_file_output=0; final link output only"
        ),
        "summary": result.summary(),
    }
    (destination / "calibration_run_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest

def _smoke_audit(link_path: Path, destination: Path, period: str) -> None:
    """Create deterministic native-contract outputs for CI/orchestration tests."""
    with link_path.open("r", newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        rows = list(reader)
    required = {
        "link_id",
        "from_node_id",
        "to_node_id",
        "vdf_plf",
        "vdf_cd",
        "vdf_cp",
        "vdf_n",
        "vdf_s",
        "vdf_alpha",
        "vdf_beta",
    }
    missing = sorted(required - set(reader.fieldnames or ()))
    if missing:
        raise ValueError("smoke link.csv is missing columns: {}".format(missing))
    audit_columns = [
        "link_id",
        "from_node_id",
        "to_node_id",
        "modeled_avg_speed_mph",
        "final_plf",
        "final_qcd",
        "final_qcp",
        "qn",
        "qs",
        "alpha",
        "beta",
    ]
    with (destination / "auto_calibration_link_audit.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=audit_columns, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            speed = row.get("speed_mph") or row.get("free_speed")
            if not speed:
                speed = 0.5 * (
                    float(row.get("qvdf_start_speed_mph", 40))
                    + float(row.get("qvdf_end_speed_mph", 40))
                )
            writer.writerow(
                {
                    "link_id": row["link_id"],
                    "from_node_id": row["from_node_id"],
                    "to_node_id": row["to_node_id"],
                    "modeled_avg_speed_mph": speed,
                    "final_plf": row["vdf_plf"],
                    "final_qcd": row["vdf_cd"],
                    "final_qcp": row["vdf_cp"],
                    "qn": row["vdf_n"],
                    "qs": row["vdf_s"],
                    "alpha": row["vdf_alpha"],
                    "beta": row["vdf_beta"],
                }
            )
    (destination / "link_performance.csv").write_text(
        "link_id,volume,speed_mph\n{}\n".format(
            "\n".join(
                "{},0,{}".format(row["link_id"], row.get("speed_mph", 40))
                for row in rows
            )
        ),
        encoding="utf-8",
    )
    (destination / "auto_calibration_history.csv").write_text(
        "outer_iteration,objective\n0,0\n", encoding="utf-8"
    )
    (destination / "auto_calibration_volume_constraint_audit.csv").write_text(
        "link_id,status\n", encoding="utf-8"
    )
    (destination / "auto_calibration_summary.json").write_text(
        json.dumps(
            {
                "status": "complete",
                "backend": "smoke",
                "period": period.upper(),
                "guardrails_pass": True,
                "final_objective": 0.0,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def _run_smoke_period(
    config: AutoCalibrationConfig,
    source: Path,
    destination: Path,
    period: str,
    processors: int,
    project_dir: Path,
) -> Dict[str, object]:
    staging = stage_scenario(source, destination)
    settings_source = config.calibration_settings_csv or (
        source / "auto_calibration_settings.csv"
    )
    shutil.copy2(settings_source, destination / "auto_calibration_settings.csv")
    contract = validate_qvdf_profile(destination / "link.csv")
    _smoke_audit(destination / "link.csv", destination, period)
    manifest: Dict[str, object] = {
        "status": "complete",
        "period": period.upper(),
        "backend": "smoke",
        "source": portable_path(source, project_dir),
        "destination": portable_path(destination, project_dir),
        "processors": processors,
        "profile_contract": contract,
        "staging": staging,
        "intermediate_output_policy": "simulated native contract; no route files",
    }
    (destination / "calibration_run_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def _read_audit(path: Path) -> Dict[Tuple[int, int], Dict[str, str]]:
    required = {
        "link_id",
        "from_node_id",
        "to_node_id",
        "modeled_avg_speed_mph",
        *PARAMETER_COLUMNS.values(),
    }
    rows: Dict[Tuple[int, int], Dict[str, str]] = {}
    with path.open("r", newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        missing = sorted(required - set(reader.fieldnames or ()))
        if missing:
            raise ValueError("{} is missing columns: {}".format(path, missing))
        for row_number, row in enumerate(reader, start=2):
            pair = (int(row["from_node_id"]), int(row["to_node_id"]))
            if pair in rows:
                raise ValueError("{}:{}: duplicate node pair".format(path, row_number))
            for name in {"modeled_avg_speed_mph", *PARAMETER_COLUMNS.values()}:
                if not _finite_positive(row[name]):
                    raise ValueError("{}:{}: invalid {}".format(path, row_number, name))
            rows[pair] = row
    return rows


def _anchors(
    pair: Tuple[int, int],
    audits: Mapping[str, Mapping[Tuple[int, int], Mapping[str, str]]],
) -> Dict[str, Tuple[float, float]]:
    def nearest(preferred: str) -> float:
        order = {
            "am": ("am", "md", "pm"),
            "md": ("md", "am", "pm"),
            "pm": ("pm", "md", "am"),
        }[preferred]
        for period in order:
            row = audits[period].get(pair)
            if row is not None:
                return float(row["modeled_avg_speed_mph"])
        raise KeyError(pair)

    am = nearest("am")
    md = nearest("md")
    pm = nearest("pm")
    am_md = 0.5 * (am + md)
    md_pm = 0.5 * (md + pm)
    return {"am": (am, am_md), "md": (am_md, md_pm), "pm": (md_pm, pm)}


def _write_calibrated_link(
    source: Path,
    destination: Path,
    period: str,
    audits: Mapping[str, Mapping[Tuple[int, int], Mapping[str, str]]],
) -> Dict[str, int]:
    audit = audits[period]
    destination.parent.mkdir(parents=True, exist_ok=True)
    matched: Set[Tuple[int, int]] = set()
    with source.open("r", newline="", encoding="utf-8-sig") as src, destination.open(
        "w", newline="", encoding="utf-8"
    ) as dst:
        reader = csv.reader(src)
        writer = csv.writer(dst, lineterminator="\n")
        header = next(reader)
        required = {
            "from_node_id",
            "to_node_id",
            "qvdf_start_speed_mph",
            "qvdf_end_speed_mph",
            *PARAMETER_COLUMNS,
        }
        missing = sorted(required - set(header))
        if missing:
            raise ValueError("{} is missing columns: {}".format(source, missing))
        index = {name: header.index(name) for name in required}
        writer.writerow(header)
        row_count = 0
        for row_number, row in enumerate(reader, start=2):
            row_count += 1
            pair = (
                int(float(row[index["from_node_id"]])),
                int(float(row[index["to_node_id"]])),
            )
            accepted = audit.get(pair)
            if accepted is None:
                raise ValueError("{}:{}: no audit row for {}".format(source, row_number, pair))
            matched.add(pair)
            for target, source_name in PARAMETER_COLUMNS.items():
                row[index[target]] = accepted[source_name]
            start, end = _anchors(pair, audits)[period]
            row[index["qvdf_start_speed_mph"]] = "{:.12g}".format(start)
            row[index["qvdf_end_speed_mph"]] = "{:.12g}".format(end)
            writer.writerow(row)
    if set(audit) - matched:
        raise ValueError("{} contains audit node pairs absent from link.csv".format(source))
    return {"network_links": row_count, "matched_audit_links": len(matched)}


def finalize_periods(
    calibration_root: Path, output: Path, project_dir: Path
) -> Dict[str, object]:
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    audits = {
        period: _read_audit(
            calibration_root / period / "auto_calibration_link_audit.csv"
        )
        for period in CALIBRATION_PERIODS
    }
    period_manifest: List[Dict[str, object]] = []
    for period in CALIBRATION_PERIODS:
        source_dir = calibration_root / period
        destination_dir = output / "assignment" / period
        destination = destination_dir / "link_calibrated.csv"
        stats = _write_calibrated_link(source_dir / "link.csv", destination, period, audits)
        for name in (
            "auto_calibration_summary.json",
            "auto_calibration_history.csv",
            "auto_calibration_link_audit.csv",
            "auto_calibration_oracle_audit.csv",
            "auto_calibration_volume_constraint_audit.csv",
            "calibration_run_manifest.json",
        ):
            source = source_dir / name
            if source.is_file():
                shutil.copy2(source, destination_dir / name)
        period_manifest.append(
            {
                "period": period.upper(),
                "link_calibrated": portable_path(destination, project_dir),
                "link_calibrated_sha256": _sha256(destination),
                **stats,
            }
        )
    manifest: Dict[str, object] = {
        "status": "complete",
        "stage": "auto-calibration-finalization",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "parameter_merge": PARAMETER_COLUMNS,
        "anchor_speed_policy": {
            "AM_start": "AM final period-average speed",
            "AM_end_MD_start": "mean of AM and MD final period-average speeds",
            "MD_end_PM_start": "mean of MD and PM final period-average speeds",
            "PM_end": "PM final period-average speed",
        },
        "periods": period_manifest,
    }
    (output / "finalization_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def build_qvdf_dictionary(
    finalized_root: Path,
    output: Path,
    fallback_path: Optional[Path],
    project_dir: Path,
) -> Dict[str, object]:
    fallback: Dict[Tuple[int, int], Dict[str, Any]] = {}
    if fallback_path is not None:
        loaded = np.load(fallback_path, allow_pickle=True).item()
        if not isinstance(loaded, dict):
            raise TypeError("fallback QVDF file does not contain a dictionary")
        fallback = loaded
    lookup: Dict[Tuple[int, int], Dict[str, Any]] = {
        (int(pair[0]), int(pair[1])): dict(record)
        for pair, record in fallback.items()
    }
    period_pairs: Dict[str, Set[Tuple[int, int]]] = {}
    required = {"link_id", "from_node_id", "to_node_id", *DICTIONARY_PARAMETERS.values()}
    for period, sequence in PERIOD_SEQUENCE.items():
        link_path = finalized_root / "assignment" / period / "link_calibrated.csv"
        seen: Set[Tuple[int, int]] = set()
        with link_path.open("r", encoding="utf-8-sig", newline="") as stream:
            reader = csv.DictReader(stream)
            missing = sorted(required - set(reader.fieldnames or ()))
            if missing:
                raise ValueError("{} is missing columns: {}".format(link_path, missing))
            for row_number, row in enumerate(reader, start=2):
                pair = (int(float(row["from_node_id"])), int(float(row["to_node_id"])))
                if pair in seen:
                    raise ValueError("{}:{}: duplicate node pair".format(link_path, row_number))
                seen.add(pair)
                record = lookup.setdefault(
                    pair,
                    {
                        "data_type": "node_pair",
                        "from_node_id": pair[0],
                        "to_node_id": pair[1],
                    },
                )
                record.update(
                    {
                        "data_type": "node_pair",
                        "link_id": _text(row.get("link_id", "")),
                        "tmc_corridor_name": _text(row.get("corridor", "")),
                        "from_node_id": pair[0],
                        "to_node_id": pair[1],
                        "vdf_code": _text(row.get("link_type", "")),
                    }
                )
                for parameter, column in DICTIONARY_PARAMETERS.items():
                    value = row[column]
                    if not _finite_positive(value):
                        raise ValueError(
                            "{}:{}: invalid {}={!r}".format(
                                link_path, row_number, column, value
                            )
                        )
                    record["QVDF_{}{}".format(parameter, sequence)] = float(value)
        period_pairs[period] = seen

    expected_keys = {
        "QVDF_{}{}".format(parameter, sequence)
        for sequence in PERIOD_SEQUENCE.values()
        for parameter in DICTIONARY_PARAMETERS
    }
    missing_cells: List[Tuple[Tuple[int, int], str]] = []
    for pair, record in lookup.items():
        for sequence in PERIOD_SEQUENCE.values():
            for parameter in DICTIONARY_PARAMETERS:
                key = "QVDF_{}{}".format(parameter, sequence)
                if _finite_positive(record.get(key)):
                    continue
                candidates = [
                    (abs(other - sequence), other)
                    for other in PERIOD_SEQUENCE.values()
                    if _finite_positive(record.get("QVDF_{}{}".format(parameter, other)))
                ]
                if candidates:
                    _, nearest = min(candidates)
                    record[key] = float(record["QVDF_{}{}".format(parameter, nearest)])
                else:
                    missing_cells.append((pair, key))
        if expected_keys - set(record):
            continue
    if missing_cells:
        preview = ", ".join("{}:{}".format(pair, key) for pair, key in missing_cells[:5])
        raise ValueError(
            "{} QVDF cells are unavailable; provide a fallback dictionary. {}".format(
                len(missing_cells), preview
            )
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    np.save(output, lookup, allow_pickle=True)
    reloaded = np.load(output, allow_pickle=True).item()
    if set(reloaded) != set(lookup):
        raise RuntimeError("saved QVDF dictionary failed round-trip validation")
    manifest: Dict[str, object] = {
        "status": "PASS",
        "schema": "mode15_link_qvdf_node_pair_dict-compatible",
        "output_dictionary": portable_path(output, project_dir),
        "output_sha256": _sha256(output),
        "dictionary_node_pairs": len(reloaded),
        "parameters": list(DICTIONARY_PARAMETERS),
        "period_sequence": PERIOD_SEQUENCE,
        "fallback_dictionary": (
            portable_path(fallback_path, project_dir) if fallback_path is not None else None
        ),
    }
    manifest_path = output.with_name(output.stem + "_manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def run_auto_calibration(
    config: AutoCalibrationConfig,
    source_root: Path,
    calibration_root: Path,
    finalized_root: Path,
    dictionary_output: Path,
    processors: int,
    project_dir: Path,
) -> Dict[str, object]:
    started = time.time()
    calibration_root.mkdir(parents=True, exist_ok=False)
    runner = _run_native_period if config.backend == "native" else _run_smoke_period
    periods: List[Dict[str, object]] = []
    # Period solves are intentionally sequential. Each native solve owns the full
    # processor budget, avoiding the memory failures seen when regional solves overlap.
    for period in config.periods:
        LOGGER.info(
            "Starting %s auto calibration (%s backend, %d processors)",
            period.upper(),
            config.backend,
            processors,
        )
        periods.append(
            runner(
                config,
                source_root / period,
                calibration_root / period,
                period,
                processors,
                project_dir,
            )
        )
    finalization = finalize_periods(calibration_root, finalized_root, project_dir)
    dictionary = build_qvdf_dictionary(
        finalized_root,
        dictionary_output,
        config.fallback_qvdf_dictionary,
        project_dir,
    )
    manifest: Dict[str, object] = {
        "status": "complete",
        "stage": "auto-calibration",
        "backend": config.backend,
        "period_execution": "sequential, with native OpenMP inside each solve",
        "processors_per_solve": processors,
        "source_scenario_root": portable_path(source_root, project_dir),
        "periods": periods,
        "finalization": finalization,
        "qvdf_dictionary": dictionary,
        "elapsed_seconds": time.time() - started,
    }
    (calibration_root / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest
