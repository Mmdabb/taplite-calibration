from __future__ import annotations

from pathlib import Path

from taplite_calibration.scaffold import initialize_project


def test_scaffold_includes_prepared_and_raw_input_layout(tmp_path: Path) -> None:
    created = initialize_project(tmp_path / "project")
    project = tmp_path / "project"

    assert created
    assert (project / "inputs" / "scenarios" / "am").is_dir()
    assert (project / "inputs" / "raw" / "scenarios" / "am").is_dir()
    assert (project / "inputs" / "raw" / "odme-policies" / "nt").is_dir()
    assert (project / "inputs" / "raw" / "route-runs" / "pm").is_dir()
    assert (project / "inputs" / "raw" / "tmc-observation-coverage").is_dir()
    settings = (
        project
        / "inputs"
        / "auto-calibration"
        / "refined_auto_calibration_settings.csv"
    )
    assert settings.is_file()
    assert "calibration_fit_mode,refined_fixed_point" in settings.read_text(
        encoding="utf-8"
    )
    config = (project / "calibration.toml").read_text(encoding="utf-8")
    assert '[prepare.odme]' in config
    assert '[prepare.auto_calibration]' in config
