from pathlib import Path

import numpy as np

from taplite_calibration import assign, auto_calibrate, native_status
from taplite_calibration.auto_calibration import run_auto_calibration
from taplite_calibration.config import AutoCalibrationConfig


def test_bundled_native_kernel_auto_calibrates(tmp_path: Path) -> None:
    source = Path(__file__).parents[1] / "examples" / "native-kernel"
    scenario = tmp_path / "scenario"
    scenario.mkdir()
    for path in source.glob("*.csv"):
        (scenario / path.name).write_bytes(path.read_bytes())

    status = native_status(1)
    assert status["probe_team_size"] >= 1

    result = auto_calibrate(
        scenario,
        timeout=120,
        settings_overrides={"number_of_processors": 1},
    )
    assert result.returncode == 0, result.log[-4000:]
    assert len(result.links) == 4
    assert (scenario / "auto_calibration_history.csv").is_file()
    assert (scenario / "auto_calibration_link_audit.csv").is_file()
    assert (scenario / "auto_calibration_summary.json").is_file()


def test_bundled_native_kernel_runs_ordinary_assignment(tmp_path: Path) -> None:
    source = Path(__file__).parents[1] / "examples" / "native-kernel"
    scenario = tmp_path / "assignment"
    scenario.mkdir()
    for path in source.glob("*.csv"):
        (scenario / path.name).write_bytes(path.read_bytes())

    result = assign(
        scenario,
        timeout=120,
        settings_overrides={"number_of_processors": 1, "auto_calibration": 0},
    )
    assert result.returncode == 0, result.log[-4000:]
    assert len(result.links) == 4
    assert not (scenario / "auto_calibration_history.csv").exists()


def test_three_period_native_orchestration_builds_dictionary(tmp_path: Path) -> None:
    source = Path(__file__).parents[1] / "examples" / "native-kernel"
    scenario_root = tmp_path / "inputs" / "scenarios"
    for period in ("am", "md", "pm"):
        scenario = scenario_root / period
        scenario.mkdir(parents=True)
        for path in source.glob("*.csv"):
            (scenario / path.name).write_bytes(path.read_bytes())

    config = AutoCalibrationConfig(
        scenario_root=scenario_root,
        backend="native",
        periods=("am", "md", "pm"),
        timeout_seconds=120,
        calibration_settings_csv=None,
        fallback_qvdf_dictionary=None,
    )
    dictionary = tmp_path / "outputs" / "calibrated_qvdf.npy"
    manifest = run_auto_calibration(
        config,
        scenario_root,
        tmp_path / "outputs" / "period-runs",
        tmp_path / "outputs" / "finalized",
        dictionary,
        processors=1,
        project_dir=tmp_path,
    )
    assert manifest["status"] == "complete"
    assert manifest["qvdf_dictionary"]["dictionary_node_pairs"] == 4
    assert dictionary.is_file()
    calibrated = np.load(dictionary, allow_pickle=True).item()
    assert len(calibrated) == 4
    assert all(
        "QVDF_{}{}".format(parameter, sequence) in values
        for values in calibrated.values()
        for parameter in ("plf", "qdf", "n", "s", "cp", "cd", "alpha", "beta")
        for sequence in (1, 2, 3)
    )
