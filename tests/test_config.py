from pathlib import Path

import pytest

from taplite_calibration.config import load_config
from taplite_calibration.errors import ConfigurationError


def test_absolute_config_wires_are_rejected(tmp_path: Path) -> None:
    absolute = (tmp_path / "inputs").resolve()
    (tmp_path / "calibration.toml").write_text(
        """schema_version = 1
[project]
input_dir = "{}"
output_dir = "outputs"
[pipeline]
mode = "odme"
processors = 1
[odme]
policy_root = "inputs/policies"
scenario_root = "inputs/scenarios"
daily_target_csv = "inputs/screens.csv"
factor_cap = 0.2
""".format(str(absolute).replace("\\", "\\\\")),
        encoding="utf-8",
    )
    with pytest.raises(ConfigurationError, match="must be relative"):
        load_config(tmp_path)


def test_prepare_is_not_a_public_mode(tmp_path: Path) -> None:
    (tmp_path / "calibration.toml").write_text(
        """schema_version = 1
[project]
input_dir = "inputs"
output_dir = "outputs"
[pipeline]
mode = "prepare"
processors = 1
""",
        encoding="utf-8",
    )

    with pytest.raises(ConfigurationError, match="mode must be one of"):
        load_config(tmp_path)
