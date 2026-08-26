"""Structured QA for prepared ODME and auto-calibration inputs."""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .config import AutoCalibrationConfig, OdmeConfig, PERIODS


AUTO_LINK_COLUMNS = {
    "link_id",
    "from_node_id",
    "to_node_id",
    "vdf_type",
    "qvdf_profile_mode",
    "calibration_observation_class",
    "facility_class",
    "target_tmc",
    "observed_p_hr",
    "observed_vt2_mph",
    "observed_avg_speed_mph",
    "s3_volume",
    "cube_vehicle_volume",
    "observation_quality",
    "mode1_plf",
    "mode1_qcd",
    "mode1_qcp",
    "cutoff_speed",
}


@dataclass(frozen=True)
class QaIssue:
    severity: str
    code: str
    message: str
    path: Optional[str] = None

    def to_dict(self) -> Dict[str, object]:
        result: Dict[str, object] = {
            "severity": self.severity,
            "code": self.code,
            "message": self.message,
        }
        if self.path is not None:
            result["path"] = self.path
        return result


@dataclass
class QaReport:
    target: str
    checks: Dict[str, object] = field(default_factory=dict)
    issues: List[QaIssue] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        return not any(issue.severity == "error" for issue in self.issues)

    def error(self, code: str, message: str, path: Optional[Path] = None) -> None:
        self.issues.append(
            QaIssue("error", code, message, str(path) if path is not None else None)
        )

    def warning(self, code: str, message: str, path: Optional[Path] = None) -> None:
        self.issues.append(
            QaIssue("warning", code, message, str(path) if path is not None else None)
        )

    def to_dict(self) -> Dict[str, object]:
        return {
            "target": self.target,
            "status": "PASS" if self.passed else "FAIL",
            "checks": self.checks,
            "issues": [issue.to_dict() for issue in self.issues],
        }


def qa_odme(config: OdmeConfig) -> QaReport:
    report = QaReport("odme")
    target = config.daily_target_csv
    screens: Optional[np.ndarray] = None
    if not target.is_file():
        report.error("missing_daily_targets", "prepared daily screen target CSV is missing", target)
    else:
        try:
            frame = pd.read_csv(target)
            required = {"screen_id", config.daily_target_column}
            missing = sorted(required - set(frame.columns))
            if missing:
                report.error(
                    "daily_target_schema",
                    "daily target CSV is missing columns: {}".format(missing),
                    target,
                )
            else:
                screen_series = pd.to_numeric(frame["screen_id"], errors="coerce")
                values = pd.to_numeric(frame[config.daily_target_column], errors="coerce")
                if screen_series.isna().any() or screen_series.duplicated().any():
                    report.error(
                        "daily_target_screen_ids",
                        "screen_id values must be numeric and unique",
                        target,
                    )
                if values.isna().any() or values.le(0).any():
                    report.error(
                        "daily_target_values",
                        "daily screen targets must be finite and positive",
                        target,
                    )
                if not report.issues:
                    screens = screen_series.astype(np.int32).sort_values().to_numpy()
                    report.checks["daily_screen_count"] = int(len(screens))
        except Exception as error:
            report.error("daily_target_read", str(error), target)

    common_zones: Optional[np.ndarray] = None
    common_screens: Optional[np.ndarray] = None
    policy_count = 0
    supported_cells = 0
    for period in PERIODS:
        scenario = config.scenario_root / period
        mode_path = scenario / "mode_type.csv"
        if not scenario.is_dir():
            report.error(
                "missing_period_scenario",
                "prepared {} scenario directory is missing".format(period.upper()),
                scenario,
            )
            continue
        if not mode_path.is_file():
            report.error("missing_mode_type", "mode_type.csv is missing", mode_path)
            continue
        try:
            modes = pd.read_csv(mode_path, low_memory=False)
        except Exception as error:
            report.error("mode_type_read", str(error), mode_path)
            continue
        if not {"mode_type", "demand_file"}.issubset(modes.columns):
            report.error(
                "mode_type_schema",
                "mode_type.csv requires mode_type and demand_file",
                mode_path,
            )
            continue
        for row in modes.to_dict("records"):
            mode = str(row["mode_type"])
            demand = scenario / str(row["demand_file"])
            if not demand.is_file() and not demand.with_suffix(".bin").is_file():
                report.error(
                    "missing_demand",
                    "demand input is missing for {}:{}".format(period, mode),
                    demand,
                )
            policy = config.policy_root / period / "od_screen_policy_{}.npz".format(mode)
            if not policy.is_file():
                report.error(
                    "missing_policy",
                    "OD-to-screen policy is missing for {}:{}".format(period, mode),
                    policy,
                )
                continue
            try:
                with np.load(policy, allow_pickle=False) as data:
                    required_arrays = {
                        "origin",
                        "destination",
                        "q0",
                        "screen_ids",
                        "zone_external",
                        "data",
                        "indices",
                        "indptr",
                    }
                    missing_arrays = sorted(required_arrays - set(data.files))
                    if missing_arrays:
                        raise ValueError("missing arrays: {}".format(missing_arrays))
                    q0 = data["q0"]
                    policy_screens = data["screen_ids"].astype(np.int32)
                    zones = data["zone_external"].astype(np.int32)
                    if q0.size == 0 or np.any(~np.isfinite(q0)) or np.any(q0 <= 0):
                        raise ValueError("q0 must contain positive finite supported cells")
                    if data["origin"].size != q0.size or data["destination"].size != q0.size:
                        raise ValueError("origin/destination/q0 lengths differ")
                    if data["indptr"].size != q0.size + 1:
                        raise ValueError("CSR indptr length must equal q0 length + 1")
                    if common_screens is None:
                        common_screens = policy_screens
                        common_zones = zones
                    elif not np.array_equal(common_screens, policy_screens) or not np.array_equal(
                        common_zones, zones
                    ):
                        raise ValueError("policy screen or zone axes differ")
                    policy_count += 1
                    supported_cells += int(q0.size)
            except Exception as error:
                report.error("policy_schema", str(error), policy)
    if screens is not None and common_screens is not None:
        if not np.array_equal(np.sort(common_screens), screens):
            report.error(
                "screen_axis_mismatch",
                "daily target and OD policy screen axes differ",
                target,
            )
    report.checks["policy_files"] = policy_count
    report.checks["supported_od_cells"] = supported_cells
    return report


def _positive(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    return values.notna() & np.isfinite(values) & values.gt(0)


def qa_auto_calibration(
    config: AutoCalibrationConfig, scenario_root: Optional[Path] = None
) -> QaReport:
    report = QaReport("auto-calibration")
    root = scenario_root or config.scenario_root
    total_links = 0
    target_links = 0
    episode_links = 0
    links_without_cube_volume = 0
    for period in config.periods:
        scenario = root / period
        link_path = scenario / "link.csv"
        if not scenario.is_dir():
            report.error(
                "missing_period_scenario",
                "prepared {} scenario directory is missing".format(period.upper()),
                scenario,
            )
            continue
        required_files = [scenario / "settings.csv", scenario / "departure_profiles.csv"]
        settings_path = config.calibration_settings_csv or (
            scenario / "auto_calibration_settings.csv"
        )
        required_files.append(settings_path)
        for path in required_files:
            if not path.is_file():
                report.error("missing_auto_input", "required prepared file is missing", path)
        if not link_path.is_file():
            report.error("missing_link", "prepared link.csv is missing", link_path)
            continue
        try:
            frame = pd.read_csv(link_path, low_memory=False)
        except Exception as error:
            report.error("link_read", str(error), link_path)
            continue
        missing = sorted(AUTO_LINK_COLUMNS - set(frame.columns))
        if missing:
            report.error(
                "auto_target_schema",
                "link.csv is missing prepared target columns: {}".format(missing),
                link_path,
            )
            continue
        total_links += len(frame)
        vdf = pd.to_numeric(frame["vdf_type"], errors="coerce")
        profile = pd.to_numeric(frame["qvdf_profile_mode"], errors="coerce")
        if not vdf.eq(2).all() or not profile.eq(1).all():
            report.error(
                "qvdf_contract",
                "every link must use vdf_type=2 and qvdf_profile_mode=1",
                link_path,
            )
        observation_class = frame["calibration_observation_class"].astype(str).str.upper()
        invalid_class = ~observation_class.isin(["E", "N", "U"])
        if invalid_class.any():
            report.error(
                "observation_class",
                "calibration_observation_class must be E, N, or U",
                link_path,
            )
        calibrated = observation_class.isin(["E", "N"])
        episode = observation_class.eq("E")
        if calibrated.any():
            for column in (
                "observed_avg_speed_mph",
                "s3_volume",
                "observation_quality",
                "mode1_plf",
                "mode1_qcd",
                "mode1_qcp",
                "cutoff_speed",
            ):
                if not _positive(frame.loc[calibrated, column]).all():
                    report.error(
                        "invalid_link_target",
                        "{} must be positive for every E/N link".format(column),
                        link_path,
                    )
            cube = pd.to_numeric(
                frame.loc[calibrated, "cube_vehicle_volume"], errors="coerce"
            )
            missing_cube = ~(cube.notna() & np.isfinite(cube) & cube.gt(0))
            if missing_cube.any():
                count = int(missing_cube.sum())
                links_without_cube_volume += count
                report.warning(
                    "missing_cube_volume",
                    (
                        "{} E/N links have no positive CUBE volume; the controller "
                        "will use their required S3 target without a two-source envelope"
                    ).format(count),
                    link_path,
                )
        if episode.any():
            for column in ("observed_p_hr", "observed_vt2_mph"):
                if not _positive(frame.loc[episode, column]).all():
                    report.error(
                        "invalid_episode_target",
                        "{} must be positive for every E link".format(column),
                        link_path,
                    )
        target_links += int(calibrated.sum())
        episode_links += int(episode.sum())
        if not calibrated.any():
            report.warning(
                "no_calibrated_links",
                "{} has no E/N calibration targets".format(period.upper()),
                link_path,
            )
    report.checks.update(
        {
            "network_links_across_periods": int(total_links),
            "calibrated_link_period_rows": int(target_links),
            "episode_link_period_rows": int(episode_links),
            "calibrated_rows_without_positive_cube_volume": int(
                links_without_cube_volume
            ),
        }
    )
    return report
