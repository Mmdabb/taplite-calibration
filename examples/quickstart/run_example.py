"""Run and verify the complete synthetic quickstart."""

from __future__ import annotations

import csv
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from taplite_calibration import run_project


PROJECT = Path(__file__).resolve().parent


def main() -> None:
    run_id = "quickstart-{}".format(
        datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    )
    result = run_project(PROJECT, run_id=run_id)
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
    if abs(after - observed) >= abs(before - observed):
        raise RuntimeError("quickstart ODME did not improve the screen error")
    dictionary_path = (
        result.run_dir / "artifacts" / "calibrated_qvdf_node_pair_dict.npy"
    )
    dictionary = np.load(dictionary_path, allow_pickle=True).item()
    if set(dictionary) != {(100, 101), (101, 102)}:
        raise RuntimeError("quickstart QVDF dictionary has unexpected node pairs")
    print(
        json.dumps(
            {
                "status": "complete",
                "run_directory": str(result.run_dir),
                "automatic_qa": result.manifest["stages"][0]["status"],
                "screen_observed": observed,
                "screen_before_odme": before,
                "screen_after_odme": after,
                "qvdf_dictionary": str(dictionary_path),
                "manifest": str(result.manifest_path),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
