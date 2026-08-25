"""Create a conventional, empty calibration project without copying datasets."""

from __future__ import annotations

import importlib.resources
from pathlib import Path
from typing import List


INPUT_README = """# Inputs

Large network, demand, policy, and observation files belong here and are never
included in the Python distribution.

- `scenarios/{nt,am,md,pm}/`: TAPLite scenarios. Each period needs
  `mode_type.csv`; AM/MD/PM also need a profile-consistent `link.csv` and
  `auto_calibration_settings.csv`.
- `odme/policies/{nt,am,md,pm}/`: `od_screen_policy_<mode>.npz` stores.
- `observations/daily_screens.csv`: `screen_id` and observed daily volume.
- `auto-calibration/`: packaged refined settings template and optional fallback
  QVDF dictionary.
- `raw/`: optional producer outputs used only when prepared-input QA fails.
  Keep clean source scenarios, CBI coverage snapshots, original screen counts,
  and either prebuilt policies or DTAC-v2 route runs under this directory.
"""

OUTPUT_README = """# Outputs

Every execution creates `runs/<run-id>/` with a top-level `run_manifest.json`,
numbered stage directories, an `artifacts/` index, and portable dictionaries.
This directory should not be committed to the package repository.
"""


def initialize_project(project_dir: Path) -> List[Path]:
    project_dir = project_dir.resolve()
    project_dir.mkdir(parents=True, exist_ok=True)
    config_path = project_dir / "calibration.toml"
    if config_path.exists():
        raise FileExistsError(config_path)
    resource = importlib.resources.files("taplite_calibration.resources").joinpath(
        "default_config.toml"
    )
    config_path.write_text(resource.read_text(encoding="utf-8"), encoding="utf-8")
    created = [config_path]
    directories = [
        "inputs/scenarios/nt",
        "inputs/scenarios/am",
        "inputs/scenarios/md",
        "inputs/scenarios/pm",
        "inputs/odme/policies/nt",
        "inputs/odme/policies/am",
        "inputs/odme/policies/md",
        "inputs/odme/policies/pm",
        "inputs/observations",
        "inputs/auto-calibration",
        "inputs/raw/scenarios/nt",
        "inputs/raw/scenarios/am",
        "inputs/raw/scenarios/md",
        "inputs/raw/scenarios/pm",
        "inputs/raw/odme-policies/nt",
        "inputs/raw/odme-policies/am",
        "inputs/raw/odme-policies/md",
        "inputs/raw/odme-policies/pm",
        "inputs/raw/route-runs/nt",
        "inputs/raw/route-runs/am",
        "inputs/raw/route-runs/md",
        "inputs/raw/route-runs/pm",
        "inputs/raw/tmc-observation-coverage",
        "outputs/runs",
    ]
    for value in directories:
        path = project_dir / value
        path.mkdir(parents=True, exist_ok=True)
        created.append(path)
    settings_resource = importlib.resources.files(
        "taplite_calibration.resources"
    ).joinpath("refined_auto_calibration_settings.csv")
    settings_path = (
        project_dir
        / "inputs"
        / "auto-calibration"
        / "refined_auto_calibration_settings.csv"
    )
    settings_path.write_text(
        settings_resource.read_text(encoding="utf-8"), encoding="utf-8"
    )
    created.append(settings_path)
    (project_dir / "inputs" / "README.md").write_text(INPUT_README, encoding="utf-8")
    (project_dir / "outputs" / "README.md").write_text(OUTPUT_README, encoding="utf-8")
    created.extend([project_dir / "inputs" / "README.md", project_dir / "outputs" / "README.md"])
    return created
