from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from taplite_calibration.errors import ConfigurationError
from taplite_calibration.pipeline import run_project


PERIODS = ("nt", "am", "md", "pm")


def _write_link(path: Path, speed: float) -> None:
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
                "corridor": "smoke",
                "calibration_observation_class": "N",
                "calibration_exclusion_reason": "",
                "facility_class": "gp",
                "target_tmc": "TMC-{}".format(link_id),
                "observed_p_hr": "",
                "observed_vt2_mph": "",
                "observed_avg_speed_mph": speed,
                "s3_volume": 100.0,
                "cube_vehicle_volume": 100.0,
                "observation_quality": 1.0,
                "observation_source": "actual",
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


def _make_project(root: Path) -> None:
    zones = np.array([1, 2, 3], dtype=np.int32)
    origin = np.array([0, 0, 1, 1, 2, 2], dtype=np.int32)
    destination = np.array([1, 2, 0, 2, 0, 1], dtype=np.int32)
    q0 = np.full(6, 10.0, dtype=np.float64)
    # The counted routes 1->2, 2->3, and 3->1 form a zero-marginal cycle.
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
        scenario = root / "inputs" / "scenarios" / period
        policy = root / "inputs" / "odme" / "policies" / period
        scenario.mkdir(parents=True)
        policy.mkdir(parents=True)
        pd.DataFrame(
            [{"mode_type": "sov", "demand_file": "demand.csv", "pce": 1.0}]
        ).to_csv(scenario / "mode_type.csv", index=False)
        demand.to_csv(scenario / "demand.csv", index=False)
        _write_link(scenario / "link.csv", 42.0 + index)
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
    observations = root / "inputs" / "observations"
    observations.mkdir(parents=True)
    pd.DataFrame(
        [{"screen_id": 1, "observed_daily_vehicle_volume": 132.0}]
    ).to_csv(observations / "daily_screens.csv", index=False)
    (root / "calibration.toml").write_text(
        """schema_version = 1

[project]
input_dir = "inputs"
output_dir = "outputs"

[pipeline]
mode = "both"
processors = 2

[odme]
policy_root = "inputs/odme/policies"
scenario_root = "inputs/scenarios"
daily_target_csv = "inputs/observations/daily_screens.csv"
daily_target_column = "observed_daily_vehicle_volume"
factor_cap = 0.20
kl_weight = 0.001
max_iterations = 5
projection_iterations = 1000
projection_tolerance = 0.000001

[auto_calibration]
scenario_root = "inputs/scenarios"
backend = "smoke"
periods = ["am", "md", "pm"]
timeout_seconds = 60
""",
        encoding="utf-8",
    )


def test_combined_smoke_runs_odme_then_auto_calibration(tmp_path: Path) -> None:
    _make_project(tmp_path)
    result = run_project(tmp_path, run_id="smoke-both")

    assert result.manifest["status"] == "complete"
    assert result.manifest["stage_order"] == [
        "prepare",
        "odme",
        "auto-calibration",
    ]
    assert result.manifest["stages"][0]["action"] == "reused_prepared_inputs"
    run_dir = tmp_path / "outputs" / "runs" / "smoke-both"
    screen = pd.read_csv(
        run_dir / "01-odme" / "results" / "screen_joint_daily_fixed_policy.csv"
    ).iloc[0]
    before = abs(
        screen["baseline_daily_fixed_policy_vehicle_volume"]
        - screen["observed_daily_vehicle_volume"]
    )
    after = abs(
        screen["joint_odme_daily_fixed_policy_vehicle_volume"]
        - screen["observed_daily_vehicle_volume"]
    )
    assert after < before

    factor_dictionary = np.load(
        run_dir / "01-odme" / "results" / "od_factor_dictionary.npy",
        allow_pickle=True,
    ).item()
    assert set(factor_dictionary["periods"]) == set(PERIODS)
    for period in PERIODS:
        adjusted = pd.read_csv(
            run_dir / "01-odme" / "adjusted-scenarios" / period / "demand.csv"
        )
        assert np.isclose(adjusted["volume"].sum(), 60.0)
        by_origin = adjusted.groupby("o_zone_id")["volume"].sum().sort_index()
        by_destination = adjusted.groupby("d_zone_id")["volume"].sum().sort_index()
        assert np.allclose(by_origin.to_numpy(), [20.0, 20.0, 20.0], atol=1e-4)
        assert np.allclose(by_destination.to_numpy(), [20.0, 20.0, 20.0], atol=1e-4)

    qvdf_path = run_dir / "artifacts" / "calibrated_qvdf_node_pair_dict.npy"
    qvdf = np.load(qvdf_path, allow_pickle=True).item()
    assert set(qvdf) == {(100, 101), (101, 102)}
    assert all(
        "QVDF_{}{}".format(parameter, sequence) in qvdf[(100, 101)]
        for parameter in ("plf", "qdf", "n", "s", "cp", "cd", "alpha", "beta")
        for sequence in (1, 2, 3)
    )
    assert (
        run_dir
        / "02-auto-calibration"
        / "finalized"
        / "assignment"
        / "am"
        / "link_calibrated.csv"
    ).is_file()
    manifest_text = result.manifest_path.read_text(encoding="utf-8")
    assert str(tmp_path.resolve()) not in manifest_text


def test_failed_target_qa_is_repaired_before_odme(tmp_path: Path) -> None:
    _make_project(tmp_path)
    prepared_target = tmp_path / "inputs" / "observations" / "daily_screens.csv"
    prepared_target.unlink()
    raw = tmp_path / "raw"
    raw.mkdir()
    pd.DataFrame([{"Screen Code": 1, "Observed 2025": 132.0}]).to_csv(
        raw / "screen-counts.csv", index=False
    )
    config_path = tmp_path / "calibration.toml"
    text = config_path.read_text(encoding="utf-8")
    text = text.replace('mode = "both"', 'mode = "odme"')
    text += """

[prepare.odme]
daily_screen_source_csv = "raw/screen-counts.csv"
screen_id_column = "Screen Code"
observed_column = "Observed 2025"
"""
    config_path.write_text(text, encoding="utf-8")

    result = run_project(tmp_path, run_id="repair-target")
    assert result.manifest["status"] == "complete"
    preparation = result.manifest["stages"][0]
    assert preparation["action"] == "prepared"
    assert any("not prepared" in warning for warning in preparation["warnings"])
    assert preparation["qa_before"]["odme"]["status"] == "FAIL"
    assert preparation["qa_after"]["odme"]["status"] == "PASS"
    assert (
        tmp_path
        / "outputs"
        / "runs"
        / "repair-target"
        / "00-prepare"
        / "odme"
        / "daily_screens.csv"
    ).is_file()


def test_selected_mode_gets_its_automatic_prepare_gate(tmp_path: Path) -> None:
    _make_project(tmp_path)
    for period in ("am", "md", "pm"):
        link_path = tmp_path / "inputs" / "scenarios" / period / "link.csv"
        frame = pd.read_csv(link_path)
        frame = frame.drop(columns=["calibration_observation_class"])
        frame.to_csv(link_path, index=False)

    result = run_project(tmp_path, mode="odme", run_id="odme-only")

    assert result.manifest["status"] == "complete"
    assert result.manifest["stage_order"] == ["prepare", "odme"]
    assert len(result.manifest["stages"]) == 2
    assert result.manifest["stages"][0]["action"] == "reused_prepared_inputs"
    assert (result.run_dir / "01-odme").is_dir()
    assert not (result.run_dir / "01-auto-calibration").exists()


def test_failed_auto_target_qa_builds_targets_then_runs(tmp_path: Path) -> None:
    _make_project(tmp_path)
    for period in ("am", "md", "pm"):
        path = tmp_path / "inputs" / "scenarios" / period / "link.csv"
        frame = pd.read_csv(path)
        frame = frame.drop(columns=["calibration_observation_class"])
        frame.to_csv(path, index=False)

    corridor = tmp_path / "raw" / "cbi" / "corridors" / "C1"
    profile_dir = corridor / "03-profiles"
    episode_dir = corridor / "04-episode-detection"
    profile_dir.mkdir(parents=True)
    episode_dir.mkdir(parents=True)
    pd.DataFrame(
        [
            {
                "tmc_code": "TMC-L1",
                "corridor": "C1",
                "t_min": minute,
                "avg_weekday_speed_mph": speed,
                "avg_weekday_flow_veh_per_hr_lane": 500.0,
                "lanes": 2,
                "n_days": 10,
                "speed_at_capacity_mph": 35.0,
            }
            for minute, speed in ((360, 40.0), (540, 45.0), (900, 38.0))
        ]
    ).to_csv(profile_dir / "average_weekday_profile.csv", index=False)
    pd.DataFrame(
        [
            {
                "tmc_code": "TMC-L1",
                "corridor": "C1",
                "period": "AM",
                "P_hr": 1.0,
                "t0_hour": 7.0,
                "t2_hour": 7.5,
                "t3_hour": 8.0,
                "min_speed_mph": 30.0,
                "qdf": 1.0,
                "plf": 1.0,
                "episode_id": "E1",
            }
        ]
    ).to_csv(
        episode_dir / "average_weekday_episode_candidates.csv", index=False
    )
    pd.DataFrame(
        [
            {
                "tmc": "TMC-L1",
                "link_id": 1,
                "from_node_id": 100,
                "to_node_id": 101,
                "distance_to_tmc_ft": 0.0,
                "road_order": 1,
                "facility_class": "gp",
            }
        ]
    ).to_csv(tmp_path / "raw" / "canonical_node_pair_tmc.csv", index=False)

    config_path = tmp_path / "calibration.toml"
    text = config_path.read_text(encoding="utf-8")
    text = text.replace('mode = "both"', 'mode = "auto-calibration"')
    text += """

[prepare.auto_calibration]
source_scenario_root = "inputs/scenarios"
cbi_actual_root = "raw/cbi"
canonical_mapping_csv = "raw/canonical_node_pair_tmc.csv"
episode_period_policy = "split_intersection"
"""
    config_path.write_text(text, encoding="utf-8")

    result = run_project(tmp_path, run_id="repair-auto-targets")
    preparation = result.manifest["stages"][0]
    assert preparation["action"] == "prepared"
    assert preparation["qa_before"]["auto-calibration"]["status"] == "FAIL"
    assert preparation["qa_after"]["auto-calibration"]["status"] == "PASS"
    prepared_link = pd.read_csv(
        tmp_path
        / "outputs"
        / "runs"
        / "repair-auto-targets"
        / "00-prepare"
        / "prepared-scenarios"
        / "am"
        / "link.csv"
    )
    selected = prepared_link.loc[prepared_link["link_id"].eq(1)].iloc[0]
    assert selected["calibration_observation_class"] == "E"
    assert selected["facility_class"] == "gp"
    assert selected["observed_p_hr"] == 1.0


def test_unprepared_and_unrepairable_input_exits_with_error(tmp_path: Path) -> None:
    _make_project(tmp_path)
    (tmp_path / "inputs" / "observations" / "daily_screens.csv").unlink()
    config_path = tmp_path / "calibration.toml"
    config_path.write_text(
        config_path.read_text(encoding="utf-8").replace(
            'mode = "both"', 'mode = "odme"'
        ),
        encoding="utf-8",
    )

    with pytest.raises(ConfigurationError, match="daily_screen_source_csv"):
        run_project(tmp_path, run_id="unrepairable")
    manifest_path = (
        tmp_path
        / "outputs"
        / "runs"
        / "unrepairable"
        / "run_manifest.json"
    )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["status"] == "failed"
    prepare_manifest = json.loads(
        (
            tmp_path
            / "outputs"
            / "runs"
            / "unrepairable"
            / "00-prepare"
            / "manifest.json"
        ).read_text(encoding="utf-8")
    )
    assert prepare_manifest["status"] == "FAIL"
    assert prepare_manifest["action"] == "repair_failed"
    assert any("not prepared" in warning for warning in prepare_manifest["warnings"])
