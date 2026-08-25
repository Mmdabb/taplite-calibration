"""TOML configuration loading and mode-aware validation."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple

try:  # Python 3.11+
    import tomllib  # type: ignore
except ModuleNotFoundError:  # pragma: no cover - exercised on Python 3.9/3.10
    import tomli as tomllib  # type: ignore

from .errors import ConfigurationError
from .paths import (
    ensure_file,
    optional_project_path,
    resolve_project_path,
)


MODES = ("odme", "auto-calibration", "both")
PERIODS = ("nt", "am", "md", "pm")
CALIBRATION_PERIODS = ("am", "md", "pm")
MAX_PROCESSORS = 20


@dataclass(frozen=True)
class OdmeConfig:
    policy_root: Path
    scenario_root: Path
    daily_target_csv: Path
    daily_target_column: str
    factor_cap: float
    kl_weight: float
    tv_budget: Optional[float]
    max_iterations: int
    initial_step: float
    maximum_step: float
    minimum_step: float
    minimum_improvement: float
    max_line_search: int
    projection_iterations: int
    projection_tolerance: float


@dataclass(frozen=True)
class AutoCalibrationConfig:
    scenario_root: Path
    backend: str
    periods: Tuple[str, ...]
    timeout_seconds: int
    calibration_settings_csv: Optional[Path]
    fallback_qvdf_dictionary: Optional[Path]


@dataclass(frozen=True)
class PrepareOdmeConfig:
    source_scenario_root: Optional[Path]
    daily_screen_source_csv: Optional[Path]
    screen_id_column: Optional[str]
    observed_column: Optional[str]
    policy_source_root: Optional[Path]
    route_run_root: Optional[Path]


@dataclass(frozen=True)
class PrepareAutoCalibrationConfig:
    source_scenario_root: Optional[Path]
    coverage_root: Optional[Path]
    cbi_actual_root: Optional[Path]
    cbi_virtual_root: Optional[Path]
    canonical_mapping_csv: Optional[Path]
    facility_mapping_csv: Optional[Path]
    departure_profile_csv: Optional[Path]
    calibration_settings_csv: Optional[Path]
    converted_network_root: Optional[Path]
    episode_period_policy: str


@dataclass(frozen=True)
class PrepareConfig:
    odme: PrepareOdmeConfig
    auto_calibration: PrepareAutoCalibrationConfig


@dataclass(frozen=True)
class ProjectConfig:
    project_dir: Path
    config_path: Path
    input_dir: Path
    output_dir: Path
    mode: str
    processors: int
    prepare: PrepareConfig
    odme: Optional[OdmeConfig]
    auto_calibration: Optional[AutoCalibrationConfig]


def _table(data: Mapping[str, Any], name: str, required: bool = True) -> Mapping[str, Any]:
    value = data.get(name)
    if value is None and not required:
        return {}
    if not isinstance(value, Mapping):
        raise ConfigurationError("[{}] must be a TOML table".format(name))
    return value


def _required(table: Mapping[str, Any], key: str, section: str) -> Any:
    if key not in table:
        raise ConfigurationError("[{}].{} is required".format(section, key))
    return table[key]


def _number(
    table: Mapping[str, Any], key: str, default: float, *, minimum: Optional[float] = None
) -> float:
    value = table.get(key, default)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError("{} must be numeric".format(key))
    result = float(value)
    if minimum is not None and result < minimum:
        raise ConfigurationError("{} must be at least {}".format(key, minimum))
    return result


def _integer(
    table: Mapping[str, Any], key: str, default: int, *, minimum: int = 1
) -> int:
    value = table.get(key, default)
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ConfigurationError("{} must be an integer >= {}".format(key, minimum))
    return value


def load_config(
    project_dir: Path,
    config_path: Optional[Path] = None,
    mode_override: Optional[str] = None,
    processors_override: Optional[int] = None,
) -> ProjectConfig:
    project_dir = project_dir.resolve()
    if config_path is None:
        config_path = project_dir / "calibration.toml"
    elif not config_path.is_absolute():
        config_path = (project_dir / config_path).resolve()
    else:
        config_path = config_path.resolve()
    ensure_file(config_path, "configuration file")
    with config_path.open("rb") as stream:
        data = tomllib.load(stream)
    if data.get("schema_version") != 1:
        raise ConfigurationError("schema_version must be 1")

    project = _table(data, "project")
    pipeline = _table(data, "pipeline")
    input_dir = resolve_project_path(
        project_dir, project.get("input_dir", "inputs"), "[project].input_dir"
    )
    output_dir = resolve_project_path(
        project_dir, project.get("output_dir", "outputs"), "[project].output_dir"
    )
    mode = mode_override or str(pipeline.get("mode", "both"))
    if mode not in MODES:
        raise ConfigurationError("mode must be one of {}".format(", ".join(MODES)))
    processors = (
        processors_override
        if processors_override is not None
        else _integer(pipeline, "processors", MAX_PROCESSORS)
    )
    if not 1 <= processors <= MAX_PROCESSORS:
        raise ConfigurationError("processors must be between 1 and {}".format(MAX_PROCESSORS))

    prepare_table = _table(data, "prepare", required=False)
    prepare_odme_table = prepare_table.get("odme", {})
    if not isinstance(prepare_odme_table, Mapping):
        raise ConfigurationError("[prepare.odme] must be a TOML table")
    prepare_auto_table = prepare_table.get("auto_calibration", {})
    if not isinstance(prepare_auto_table, Mapping):
        raise ConfigurationError("[prepare.auto_calibration] must be a TOML table")

    prepare_odme = PrepareOdmeConfig(
        source_scenario_root=optional_project_path(
            project_dir,
            prepare_odme_table.get("source_scenario_root"),
            "[prepare.odme].source_scenario_root",
        ),
        daily_screen_source_csv=optional_project_path(
            project_dir,
            prepare_odme_table.get("daily_screen_source_csv"),
            "[prepare.odme].daily_screen_source_csv",
        ),
        screen_id_column=(
            str(prepare_odme_table["screen_id_column"])
            if "screen_id_column" in prepare_odme_table
            else None
        ),
        observed_column=(
            str(prepare_odme_table["observed_column"])
            if "observed_column" in prepare_odme_table
            else None
        ),
        policy_source_root=optional_project_path(
            project_dir,
            prepare_odme_table.get("policy_source_root"),
            "[prepare.odme].policy_source_root",
        ),
        route_run_root=optional_project_path(
            project_dir,
            prepare_odme_table.get("route_run_root"),
            "[prepare.odme].route_run_root",
        ),
    )
    episode_policy = str(
        prepare_auto_table.get("episode_period_policy", "split_intersection")
    )
    if episode_policy not in ("assigned_period", "split_intersection"):
        raise ConfigurationError(
            "[prepare.auto_calibration].episode_period_policy must be "
            "assigned_period or split_intersection"
        )
    prepare_auto = PrepareAutoCalibrationConfig(
        source_scenario_root=optional_project_path(
            project_dir,
            prepare_auto_table.get("source_scenario_root"),
            "[prepare.auto_calibration].source_scenario_root",
        ),
        coverage_root=optional_project_path(
            project_dir,
            prepare_auto_table.get("coverage_root"),
            "[prepare.auto_calibration].coverage_root",
        ),
        cbi_actual_root=optional_project_path(
            project_dir,
            prepare_auto_table.get("cbi_actual_root"),
            "[prepare.auto_calibration].cbi_actual_root",
        ),
        cbi_virtual_root=optional_project_path(
            project_dir,
            prepare_auto_table.get("cbi_virtual_root"),
            "[prepare.auto_calibration].cbi_virtual_root",
        ),
        canonical_mapping_csv=optional_project_path(
            project_dir,
            prepare_auto_table.get("canonical_mapping_csv"),
            "[prepare.auto_calibration].canonical_mapping_csv",
        ),
        facility_mapping_csv=optional_project_path(
            project_dir,
            prepare_auto_table.get("facility_mapping_csv"),
            "[prepare.auto_calibration].facility_mapping_csv",
        ),
        departure_profile_csv=optional_project_path(
            project_dir,
            prepare_auto_table.get("departure_profile_csv"),
            "[prepare.auto_calibration].departure_profile_csv",
        ),
        calibration_settings_csv=optional_project_path(
            project_dir,
            prepare_auto_table.get("calibration_settings_csv"),
            "[prepare.auto_calibration].calibration_settings_csv",
        ),
        converted_network_root=optional_project_path(
            project_dir,
            prepare_auto_table.get("converted_network_root"),
            "[prepare.auto_calibration].converted_network_root",
        ),
        episode_period_policy=episode_policy,
    )
    prepare_config = PrepareConfig(
        odme=prepare_odme,
        auto_calibration=prepare_auto,
    )

    needs_odme = mode in ("odme", "both")
    needs_auto = mode in ("auto-calibration", "both")

    odme_config: Optional[OdmeConfig] = None
    if needs_odme:
        section = _table(data, "odme")
        factor_cap = _number(section, "factor_cap", 0.20, minimum=1e-12)
        if factor_cap > 1.0:
            raise ConfigurationError("[odme].factor_cap cannot exceed 1.0")
        tv_raw = section.get("tv_budget")
        tv_budget = None if tv_raw is None else float(tv_raw)
        if tv_budget is not None and not 0.0 <= tv_budget <= 1.0:
            raise ConfigurationError("[odme].tv_budget must be between 0 and 1")
        odme_config = OdmeConfig(
            policy_root=resolve_project_path(
                project_dir, _required(section, "policy_root", "odme"), "[odme].policy_root"
            ),
            scenario_root=resolve_project_path(
                project_dir, _required(section, "scenario_root", "odme"), "[odme].scenario_root"
            ),
            daily_target_csv=resolve_project_path(
                project_dir,
                _required(section, "daily_target_csv", "odme"),
                "[odme].daily_target_csv",
            ),
            daily_target_column=str(
                section.get("daily_target_column", "observed_daily_vehicle_volume")
            ),
            factor_cap=factor_cap,
            kl_weight=_number(section, "kl_weight", 0.02, minimum=0.0),
            tv_budget=tv_budget,
            max_iterations=_integer(section, "max_iterations", 12),
            initial_step=_number(section, "initial_step", 0.15, minimum=1e-12),
            maximum_step=_number(section, "maximum_step", 0.35, minimum=1e-12),
            minimum_step=_number(section, "minimum_step", 1e-4, minimum=1e-15),
            minimum_improvement=_number(
                section, "minimum_improvement", 1e-11, minimum=0.0
            ),
            max_line_search=_integer(section, "max_line_search", 8),
            projection_iterations=_integer(section, "projection_iterations", 10000),
            projection_tolerance=_number(
                section, "projection_tolerance", 2e-5, minimum=1e-15
            ),
        )

    auto_config: Optional[AutoCalibrationConfig] = None
    if needs_auto:
        section = _table(data, "auto_calibration")
        raw_periods = section.get("periods", list(CALIBRATION_PERIODS))
        if not isinstance(raw_periods, Sequence) or isinstance(raw_periods, str):
            raise ConfigurationError("[auto_calibration].periods must be an array")
        periods = tuple(str(value).lower() for value in raw_periods)
        if not periods or any(value not in CALIBRATION_PERIODS for value in periods):
            raise ConfigurationError("auto-calibration periods must be AM, MD, and/or PM")
        if len(set(periods)) != len(periods):
            raise ConfigurationError("auto-calibration periods cannot contain duplicates")
        if set(periods) != set(CALIBRATION_PERIODS):
            raise ConfigurationError(
                "this version calibrates AM, MD, and PM together so final anchors and "
                "the three-period QVDF dictionary remain consistent"
            )
        backend = str(section.get("backend", "native")).lower()
        if backend not in ("native", "smoke"):
            raise ConfigurationError("[auto_calibration].backend must be native or smoke")
        auto_config = AutoCalibrationConfig(
            scenario_root=resolve_project_path(
                project_dir,
                _required(section, "scenario_root", "auto_calibration"),
                "[auto_calibration].scenario_root",
            ),
            backend=backend,
            periods=periods,
            timeout_seconds=_integer(section, "timeout_seconds", 86400),
            calibration_settings_csv=optional_project_path(
                project_dir,
                section.get("calibration_settings_csv"),
                "[auto_calibration].calibration_settings_csv",
            ),
            fallback_qvdf_dictionary=optional_project_path(
                project_dir,
                section.get("fallback_qvdf_dictionary"),
                "[auto_calibration].fallback_qvdf_dictionary",
            ),
        )
    return ProjectConfig(
        project_dir=project_dir,
        config_path=config_path,
        input_dir=input_dir,
        output_dir=output_dir,
        mode=mode,
        processors=processors,
        prepare=prepare_config,
        odme=odme_config,
        auto_calibration=auto_config,
    )
