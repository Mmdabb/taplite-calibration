from __future__ import annotations

import csv
import shutil
from pathlib import Path

import numpy as np

from taplite_calibration.pipeline import run_project


def test_checked_in_quickstart_runs_end_to_end(tmp_path: Path) -> None:
    source = Path(__file__).resolve().parents[1] / "examples" / "quickstart"
    project = tmp_path / "quickstart"
    shutil.copytree(
        source,
        project,
        ignore=shutil.ignore_patterns("outputs", "__pycache__"),
    )

    result = run_project(project, run_id="checked-in-example")

    assert result.manifest["status"] == "complete"
    assert result.manifest["stage_order"] == [
        "prepare",
        "odme",
        "auto-calibration",
    ]
    assert result.manifest["stages"][0]["status"] == "PASS"
    comparison = (
        result.run_dir
        / "01-odme"
        / "results"
        / "screen_joint_daily_fixed_policy.csv"
    )
    with comparison.open("r", newline="", encoding="utf-8-sig") as stream:
        row = next(csv.DictReader(stream))
    observed = float(row["observed_daily_vehicle_volume"])
    before = float(row["baseline_daily_fixed_policy_vehicle_volume"])
    after = float(row["joint_odme_daily_fixed_policy_vehicle_volume"])
    assert abs(after - observed) < abs(before - observed)

    factor_path = (
        result.run_dir / "01-odme" / "results" / "od_factor_dictionary.npy"
    )
    factor_dictionary = np.load(factor_path, allow_pickle=True).item()
    assert set(factor_dictionary["periods"]) == {"nt", "am", "md", "pm"}
    factor_count = 0
    for period_entry in factor_dictionary["periods"].values():
        for shard_entry in period_entry["modes"].values():
            shard = np.load(factor_path.parent / shard_entry["file"])
            assert set(shard.files) == {
                "origin_zone_id",
                "destination_zone_id",
                "factor",
                "original_volume",
                "adjusted_volume",
            }
            assert len(shard["factor"]) == shard_entry["positive_od_cells"]
            assert np.allclose(
                shard["factor"],
                shard["adjusted_volume"] / shard["original_volume"],
                rtol=1e-6,
            )
            factor_count += len(shard["factor"])
    assert factor_count > 0

    dictionary = np.load(
        result.run_dir / "artifacts" / "calibrated_qvdf_node_pair_dict.npy",
        allow_pickle=True,
    ).item()
    assert set(dictionary) == {(100, 101), (101, 102)}
    assert all(
        "QVDF_{}{}".format(parameter, sequence) in values
        for values in dictionary.values()
        for parameter in ("plf", "qdf", "n", "s", "cp", "cd", "alpha", "beta")
        for sequence in (1, 2, 3)
    )
