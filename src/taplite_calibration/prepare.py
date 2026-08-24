"""Prepared-input QA gate and recoverable preprocessing orchestration."""

from __future__ import annotations

import json
import logging
import os
import re
import shutil
import time
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd

from .auto_targets import prepare_auto_targets
from .config import (
    AutoCalibrationConfig,
    OdmeConfig,
    PERIODS,
    PrepareAutoCalibrationConfig,
    ProjectConfig,
)
from .errors import ConfigurationError
from .paths import portable_path
from .policy_builder import build_period_policies
from .qa import QaReport, qa_auto_calibration, qa_odme


LOGGER = logging.getLogger("taplite_calibration.prepare")
SPACE = re.compile(r"\s+")
SCENARIO_OUTPUTS = {
    "route_columns.bin",
    "link_performance.csv",
    "route_assignment.csv",
    "agent.csv",
    "summary_log_file.txt",
    "final_summary.csv",
    "destination_accessibility.csv",
    "origin_accessibility.csv",
    "od_performance.csv",
    "system_performance.csv",
    "inaccessible_od.csv",
    "dtalite_run.log",
    "TAP_log.csv",
}


@dataclass(frozen=True)
class PreparationResult:
    config: ProjectConfig
    manifest: Dict[str, object]


def selected_targets(config: ProjectConfig) -> Tuple[str, ...]:
    if config.mode == "odme":
        return ("odme",)
    if config.mode == "auto-calibration":
        return ("auto-calibration",)
    return ("odme", "auto-calibration")


def run_qa(config: ProjectConfig) -> Dict[str, QaReport]:
    reports: Dict[str, QaReport] = {}
    for target in selected_targets(config):
        if target == "odme":
            if config.odme is None:
                raise ConfigurationError("ODME configuration is unavailable")
            reports[target] = qa_odme(config.odme)
        else:
            if config.auto_calibration is None:
                raise ConfigurationError("auto-calibration configuration is unavailable")
            reports[target] = qa_auto_calibration(config.auto_calibration)
    return reports


def _portable_report(report: QaReport, project_dir: Path) -> Dict[str, object]:
    data = report.to_dict()
    for issue in data["issues"]:  # type: ignore
        value = issue.get("path")
        if value:
            issue["path"] = portable_path(Path(value), project_dir)
    return data


def _issue_codes(report: QaReport) -> Set[str]:
    return {issue.code for issue in report.issues if issue.severity == "error"}


def _column_by_normalized_name(columns: Sequence[str], wanted: str) -> Optional[str]:
    normalized = SPACE.sub(" ", wanted).strip().lower()
    matches = [column for column in columns if SPACE.sub(" ", str(column)).strip().lower() == normalized]
    return matches[0] if len(matches) == 1 else None


def normalize_daily_screens(
    source: Path,
    destination: Path,
    screen_id_column: Optional[str],
    observed_column: Optional[str],
) -> Dict[str, object]:
    if not source.is_file():
        raise FileNotFoundError(source)
    frame = pd.read_csv(source, low_memory=False)
    if screen_id_column is None:
        candidates = [
            column
            for column in frame.columns
            if SPACE.sub(" ", str(column)).strip().lower()
            in {"screen_id", "screen id", "screen", "screen code", "screen_code"}
        ]
        if len(candidates) != 1:
            raise ValueError(
                "cannot infer the screen ID column; configure [prepare.odme].screen_id_column"
            )
        screen_id_column = candidates[0]
    elif screen_id_column not in frame:
        normalized = _column_by_normalized_name(frame.columns, screen_id_column)
        if normalized is None:
            raise ValueError("screen ID column {!r} is absent".format(screen_id_column))
        screen_id_column = normalized
    if observed_column is None:
        excluded = {screen_id_column}
        candidates = [
            column
            for column in frame.columns
            if column not in excluded
            and ("obs" in str(column).lower() or "count" in str(column).lower())
        ]
        positive_candidates = []
        for column in candidates:
            values = pd.to_numeric(frame[column], errors="coerce")
            if values.notna().all() and values.gt(0).all():
                positive_candidates.append(column)
        if len(positive_candidates) != 1:
            raise ValueError(
                "cannot infer one observed-volume column; configure "
                "[prepare.odme].observed_column"
            )
        observed_column = positive_candidates[0]
    elif observed_column not in frame:
        normalized = _column_by_normalized_name(frame.columns, observed_column)
        if normalized is None:
            raise ValueError("observed column {!r} is absent".format(observed_column))
        observed_column = normalized
    result = pd.DataFrame(
        {
            "screen_id": pd.to_numeric(frame[screen_id_column], errors="coerce"),
            "observed_daily_vehicle_volume": pd.to_numeric(
                frame[observed_column], errors="coerce"
            ),
        }
    ).dropna()
    result["screen_id"] = result["screen_id"].astype(np.int32)
    if result.empty or result["screen_id"].duplicated().any():
        raise ValueError("normalized screen IDs must be nonempty and unique")
    if result["observed_daily_vehicle_volume"].le(0).any():
        raise ValueError("normalized daily screen volumes must be positive")
    result = result.sort_values("screen_id", kind="mergesort")
    destination.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(destination, index=False)
    return {
        "source_rows": int(len(frame)),
        "prepared_screens": int(len(result)),
        "source_screen_id_column": str(screen_id_column),
        "source_observed_column": str(observed_column),
    }


def _stage_scenario_root(source_root: Path, destination_root: Path) -> Dict[str, object]:
    destination_root.mkdir(parents=True, exist_ok=False)
    stats: Dict[str, object] = {"periods": {}}
    for period in PERIODS:
        source = source_root / period
        if not source.is_dir():
            raise FileNotFoundError(source)
        destination = destination_root / period
        destination.mkdir()
        copied = 0
        linked = 0
        for path in source.iterdir():
            if not path.is_file() or path.name in SCENARIO_OUTPUTS:
                continue
            target = destination / path.name
            if path.name in {"link.csv", "settings.csv"}:
                shutil.copy2(path, target)
                copied += 1
            else:
                try:
                    os.link(path, target)
                    linked += 1
                except OSError:
                    shutil.copy2(path, target)
                    copied += 1
        stats["periods"][period] = {"copied": copied, "hard_linked": linked}  # type: ignore
    return stats


def _stage_policy_root(source_root: Path, destination_root: Path) -> Dict[str, object]:
    destination_root.mkdir(parents=True, exist_ok=False)
    counts: Dict[str, int] = {}
    for period in PERIODS:
        source = source_root / period
        destination = destination_root / period
        destination.mkdir()
        files = list(source.glob("od_screen_policy_*.npz"))
        if not files:
            raise FileNotFoundError(
                "no od_screen_policy_*.npz files in {}".format(source)
            )
        for path in files:
            target = destination / path.name
            try:
                os.link(path, target)
            except OSError:
                shutil.copy2(path, target)
        counts[period] = len(files)
    return {"policy_files_by_period": counts}


def _can_use_existing_daily(report: QaReport) -> bool:
    codes = _issue_codes(report)
    return not bool(
        codes
        & {
            "missing_daily_targets",
            "daily_target_schema",
            "daily_target_screen_ids",
            "daily_target_values",
            "daily_target_read",
        }
    )


def _can_use_existing_policies(report: QaReport) -> bool:
    codes = _issue_codes(report)
    return not bool(
        codes & {"missing_policy", "policy_schema", "screen_axis_mismatch"}
    )


def _can_use_existing_scenarios(report: QaReport) -> bool:
    codes = _issue_codes(report)
    return not bool(
        codes
        & {
            "missing_period_scenario",
            "missing_mode_type",
            "mode_type_read",
            "mode_type_schema",
            "missing_demand",
        }
    )


def _prepare_odme(
    config: ProjectConfig,
    report: QaReport,
    stage_root: Path,
) -> Tuple[OdmeConfig, Dict[str, object]]:
    assert config.odme is not None
    original = config.odme
    prepared = config.prepare.odme
    audit: Dict[str, object] = {}

    daily_target = original.daily_target_csv
    if not _can_use_existing_daily(report):
        if prepared.daily_screen_source_csv is None:
            raise ConfigurationError(
                "daily screen QA failed and [prepare.odme].daily_screen_source_csv "
                "is not configured"
            )
        daily_target = stage_root / "daily_screens.csv"
        audit["daily_screens"] = normalize_daily_screens(
            prepared.daily_screen_source_csv,
            daily_target,
            prepared.screen_id_column,
            prepared.observed_column,
        )

    scenario_root = original.scenario_root
    if not _can_use_existing_scenarios(report):
        source = prepared.source_scenario_root or prepared.route_run_root
        if source is None:
            raise ConfigurationError(
                "ODME scenario QA failed and neither [prepare.odme].source_scenario_root "
                "nor route_run_root is configured"
            )
        scenario_root = stage_root / "source-scenarios"
        audit["scenarios"] = _stage_scenario_root(source, scenario_root)

    policy_root = original.policy_root
    if not _can_use_existing_policies(report):
        policy_root = stage_root / "policies"
        if prepared.policy_source_root is not None:
            audit["policies"] = _stage_policy_root(
                prepared.policy_source_root, policy_root
            )
        elif prepared.route_run_root is not None:
            target_frame = pd.read_csv(daily_target)
            screens = pd.to_numeric(
                target_frame["screen_id"], errors="raise"
            ).astype(np.int32).sort_values().to_numpy()
            policy_root.mkdir(parents=True, exist_ok=False)
            period_results = []
            for period in PERIODS:
                period_results.append(
                    build_period_policies(
                        prepared.route_run_root / period,
                        policy_root / period,
                        screens,
                        period,
                    )
                )
            audit["policies"] = {
                "strategy": "streamed from DTAC-v2 route pools",
                "path_screen_incidence": "integer counted-link multiplicity",
                "periods": period_results,
            }
        else:
            raise ConfigurationError(
                "OD policy QA failed and neither [prepare.odme].policy_source_root "
                "nor route_run_root is configured"
            )
    return (
        replace(
            original,
            daily_target_csv=daily_target,
            daily_target_column="observed_daily_vehicle_volume",
            scenario_root=scenario_root,
            policy_root=policy_root,
        ),
        audit,
    )


def _auto_prepare_capability(
    prepare: PrepareAutoCalibrationConfig,
    source_scenario_root: Optional[Path],
) -> List[str]:
    missing: List[str] = []
    if source_scenario_root is None or not source_scenario_root.is_dir():
        missing.append("source_scenario_root")
    if prepare.coverage_root is None:
        if prepare.cbi_actual_root is None:
            missing.append("cbi_actual_root or coverage_root")
        if prepare.canonical_mapping_csv is None:
            missing.append("canonical_mapping_csv or coverage_root")
    elif not prepare.coverage_root.is_dir():
        missing.append("coverage_root")
    if prepare.departure_profile_csv is not None and not prepare.departure_profile_csv.is_file():
        missing.append("departure_profile_csv")
    if prepare.calibration_settings_csv is not None and not prepare.calibration_settings_csv.is_file():
        missing.append("calibration_settings_csv")
    return missing


def ensure_prepared_inputs(
    config: ProjectConfig,
    run_dir: Path,
) -> PreparationResult:
    started = time.time()
    stage_root = run_dir / "00-prepare"
    stage_root.mkdir(parents=True, exist_ok=False)
    before = run_qa(config)
    failed = [target for target, report in before.items() if not report.passed]
    warnings: List[str] = []
    for target in failed:
        warning = (
            "{} prepared-input QA did not pass; checking whether configured raw "
            "inputs can be preprocessed".format(target)
        )
        warnings.append(warning)
        LOGGER.warning(warning)

    updated = config
    preparation: Dict[str, object] = {}
    try:
        if "odme" in failed:
            warning = (
                "ODME inputs are not prepared; attempting run-local preprocessing"
            )
            warnings.append(warning)
            LOGGER.warning(warning)
            assert updated.odme is not None
            odme_config, audit = _prepare_odme(
                updated, before["odme"], stage_root / "odme"
            )
            updated = replace(updated, odme=odme_config)
            preparation["odme"] = audit

        if "auto-calibration" in failed:
            warning = (
                "auto-calibration inputs are not prepared; attempting run-local "
                "preprocessing"
            )
            warnings.append(warning)
            LOGGER.warning(warning)
            assert updated.auto_calibration is not None
            prepare_auto = updated.prepare.auto_calibration
            if (
                prepare_auto.calibration_settings_csv is None
                and updated.auto_calibration.calibration_settings_csv is not None
            ):
                prepare_auto = replace(
                    prepare_auto,
                    calibration_settings_csv=(
                        updated.auto_calibration.calibration_settings_csv
                    ),
                )
            # A shared scenario prepared for ODME is preferred so both stages
            # consume the same network/demand baseline. Otherwise use the
            # explicit raw source.
            scenario_source = (
                updated.odme.scenario_root
                if updated.odme is not None and updated.odme.scenario_root.is_dir()
                else prepare_auto.source_scenario_root
                or updated.auto_calibration.scenario_root
            )
            missing = _auto_prepare_capability(prepare_auto, scenario_source)
            if missing:
                raise ConfigurationError(
                    "auto-calibration QA failed and preprocessing cannot start; "
                    "missing: {}".format(", ".join(missing))
                )
            destination = stage_root / "prepared-scenarios"
            preparation["auto-calibration"] = prepare_auto_targets(
                prepare_auto,
                scenario_source,
                destination,
                stage_root / "auto-target-audits",
                updated.processors,
            )
            auto_config = replace(updated.auto_calibration, scenario_root=destination)
            odme_config = (
                replace(updated.odme, scenario_root=destination)
                if updated.odme is not None
                else None
            )
            updated = replace(
                updated,
                auto_calibration=auto_config,
                odme=odme_config,
            )
    except BaseException as error:
        manifest = {
            "status": "FAIL",
            "stage": "prepare",
            "action": "repair_failed",
            "warnings": warnings,
            "qa_before": {
                target: _portable_report(report, config.project_dir)
                for target, report in before.items()
            },
            "preparation": preparation,
            "failure": {
                "type": type(error).__name__,
                "message": str(error),
            },
            "elapsed_seconds": time.time() - started,
        }
        (stage_root / "manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
        )
        raise

    after = run_qa(updated)
    remaining = [target for target, report in after.items() if not report.passed]
    manifest: Dict[str, object] = {
        "status": "PASS" if not remaining else "FAIL",
        "stage": "prepare",
        "action": "prepared" if failed else "reused_prepared_inputs",
        "warnings": warnings,
        "qa_before": {
            target: _portable_report(report, config.project_dir)
            for target, report in before.items()
        },
        "preparation": preparation,
        "qa_after": {
            target: _portable_report(report, config.project_dir)
            for target, report in after.items()
        },
        "elapsed_seconds": time.time() - started,
    }
    (stage_root / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    if remaining:
        details = []
        for target in remaining:
            details.extend(
                "{}:{}".format(issue.code, issue.message)
                for issue in after[target].issues
                if issue.severity == "error"
            )
        raise ConfigurationError(
            "preprocessing completed but prepared-input QA still failed: {}".format(
                "; ".join(details[:10])
            )
        )
    if failed:
        LOGGER.warning(
            "Prepared-input QA initially failed, preprocessing completed, and the "
            "post-preparation QA now passes"
        )
    return PreparationResult(updated, manifest)
