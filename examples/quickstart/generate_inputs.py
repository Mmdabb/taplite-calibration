"""Regenerate the tiny, public quickstart input dataset."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parent
PERIODS = ("nt", "am", "md", "pm")


def write_link(path: Path, speed: float) -> None:
    columns = [
        "link_id",
        "from_node_id",
        "to_node_id",
        "vdf_type",
        "qvdf_profile_mode",
        "qvdf_start_speed_mph",
        "qvdf_end_speed_mph",
        "speed_mph",
        "vdf_plf",
        "vdf_qdf",
        "vdf_n",
        "vdf_s",
        "vdf_cp",
        "vdf_cd",
        "vdf_alpha",
        "vdf_beta",
        "mode1_plf",
        "mode1_qcd",
        "mode1_qcp",
        "link_type",
        "corridor",
        "calibration_observation_class",
        "calibration_exclusion_reason",
        "facility_class",
        "target_tmc",
        "observed_p_hr",
        "observed_vt2_mph",
        "observed_avg_speed_mph",
        "s3_volume",
        "cube_vehicle_volume",
        "observation_quality",
        "observation_source",
        "virtual_treatment",
        "virtual_confidence",
        "cutoff_speed",
        "I4AMVOL",
        "I4MDVOL",
        "I4PMVOL",
    ]
    rows = []
    for link_id, from_node, to_node in ((1, 100, 101), (2, 101, 102)):
        episode = link_id == 1
        rows.append(
            {
                "link_id": link_id,
                "from_node_id": from_node,
                "to_node_id": to_node,
                "vdf_type": 2,
                "qvdf_profile_mode": 1,
                "qvdf_start_speed_mph": speed,
                "qvdf_end_speed_mph": speed,
                "speed_mph": speed,
                "vdf_plf": 1.0,
                "vdf_qdf": 1.0,
                "vdf_n": 1.2,
                "vdf_s": 4.0,
                "vdf_cp": 0.4,
                "vdf_cd": 1.0,
                "vdf_alpha": 0.15,
                "vdf_beta": 4.0,
                "mode1_plf": 1.0,
                "mode1_qcd": 1.0,
                "mode1_qcp": 0.4,
                "link_type": 2,
                "corridor": "synthetic",
                "calibration_observation_class": "E" if episode else "N",
                "calibration_exclusion_reason": "",
                "facility_class": "gp",
                "target_tmc": "SYNTHETIC-{}".format(link_id),
                "observed_p_hr": 1.0 if episode else "",
                "observed_vt2_mph": max(speed - 12.0, 5.0) if episode else "",
                "observed_avg_speed_mph": speed,
                "s3_volume": 100.0,
                "cube_vehicle_volume": 100.0,
                "observation_quality": 1.0,
                "observation_source": "synthetic",
                "virtual_treatment": "",
                "virtual_confidence": "",
                "cutoff_speed": 35.0,
                "I4AMVOL": 100.0,
                "I4MDVOL": 100.0,
                "I4PMVOL": 100.0,
            }
        )
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    zones = np.array([1, 2, 3], dtype=np.int32)
    origin = np.array([0, 0, 1, 1, 2, 2], dtype=np.int32)
    destination = np.array([1, 2, 0, 2, 0, 1], dtype=np.int32)
    q0 = np.full(6, 10.0, dtype=np.float64)
    incidence = sparse.csr_matrix(
        np.array([[1.0], [0.0], [0.0], [1.0], [1.0], [0.0]])
    )
    demand = pd.DataFrame(
        {
            "o_zone_id": zones[origin],
            "d_zone_id": zones[destination],
            "volume": q0,
        }
    )
    for index, period in enumerate(PERIODS):
        scenario = ROOT / "inputs" / "scenarios" / period
        policy = ROOT / "inputs" / "odme" / "policies" / period
        scenario.mkdir(parents=True, exist_ok=True)
        policy.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(
            [{"mode_type": "sov", "demand_file": "demand.csv", "pce": 1.0}]
        ).to_csv(scenario / "mode_type.csv", index=False)
        demand.to_csv(scenario / "demand.csv", index=False)
        write_link(scenario / "link.csv", 42.0 + index)
        (scenario / "auto_calibration_settings.csv").write_text(
            "key,value\nworkers,2\nauto_calibration,1\n", encoding="utf-8"
        )
        (scenario / "settings.csv").write_text(
            "number_of_processors\n2\n", encoding="utf-8"
        )
        (scenario / "departure_profiles.csv").write_text(
            "time,ratio\n0,1\n", encoding="utf-8"
        )
        np.savez_compressed(
            policy / "od_screen_policy_sov.npz",
            origin=origin,
            destination=destination,
            q0=q0,
            screen_ids=np.array([1], dtype=np.int32),
            zone_external=zones,
            data=incidence.data,
            indices=incidence.indices,
            indptr=incidence.indptr,
        )
    observations = ROOT / "inputs" / "observations"
    observations.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [{"screen_id": 1, "observed_daily_vehicle_volume": 132.0}]
    ).to_csv(observations / "daily_screens.csv", index=False)
    print("Regenerated quickstart inputs under {}".format(ROOT / "inputs"))


if __name__ == "__main__":
    main()
