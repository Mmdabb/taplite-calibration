"""Prepare link-level refined auto-calibration targets from CBI producer output."""

from __future__ import annotations

import logging
import math
import os
import shutil
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .config import CALIBRATION_PERIODS, PrepareAutoCalibrationConfig


LOGGER = logging.getLogger("taplite_calibration.prepare.auto_targets")
PERIOD_MINUTES = {"AM": (360, 540), "MD": (540, 900), "PM": (900, 1140)}
TARGET_COLUMNS = [
    "calibration_observation_class",
    "calibration_exclusion_reason",
    "facility_class",
    "target_tmc",
    "corridor",
    "observed_p_hr",
    "observed_vt2_mph",
    "observed_avg_speed_mph",
    "s3_volume",
    "cube_vehicle_volume",
    "observation_quality",
    "observation_source",
    "virtual_treatment",
    "virtual_confidence",
    "mode1_plf",
    "mode1_qcd",
    "mode1_qcp",
    "cutoff_speed",
]


def resolve_sources(
    config: PrepareAutoCalibrationConfig,
) -> Tuple[List[Tuple[Path, str]], Path, Optional[pd.DataFrame], Optional[pd.DataFrame]]:
    facility_frame: Optional[pd.DataFrame] = None
    metadata: Optional[pd.DataFrame] = None
    if config.coverage_root is not None:
        coverage = config.coverage_root
        actual = config.cbi_actual_root or coverage / "cbi" / "actual"
        virtual = config.cbi_virtual_root or coverage / "cbi" / "virtual"
        runs = [(actual, "actual")]
        if virtual.is_dir():
            runs.append((virtual, "virtual"))
        mapping = config.canonical_mapping_csv or (
            coverage
            / "congestion-resources"
            / "direct-mapping"
            / "canonical_node_pair_tmc.csv"
        )
        if config.facility_mapping_csv is None:
            facility_frame, metadata = _coverage_facility_mapping(coverage, mapping)
        return runs, mapping, facility_frame, metadata
    if config.cbi_actual_root is None or config.canonical_mapping_csv is None:
        raise ValueError(
            "auto target preparation requires coverage_root or both "
            "cbi_actual_root and canonical_mapping_csv"
        )
    runs = [(config.cbi_actual_root, "actual")]
    if config.cbi_virtual_root is not None:
        runs.append((config.cbi_virtual_root, "virtual"))
    return runs, config.canonical_mapping_csv, None, None


def _coverage_facility_mapping(
    coverage_root: Path, canonical_path: Path
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    canonical = pd.read_csv(
        canonical_path,
        usecols=["tmc", "link_id", "from_node_id", "to_node_id"],
        low_memory=False,
    )
    parts: List[pd.DataFrame] = []
    sources = [
        (
            coverage_root / "actual" / "gp-canonical" / "canonical_gp_actual.csv",
            "gp",
            "tmc",
        ),
        (
            coverage_root
            / "actual"
            / "managed-canonical"
            / "canonical_managed_actual.csv",
            "managed",
            "tmc",
        ),
        (
            coverage_root
            / "actual"
            / "managed-canonical"
            / "supplemental_managed_actual.csv",
            "managed",
            "source_tmc_primary",
        ),
    ]
    for path, facility_class, tmc_column in sources:
        if not path.is_file():
            continue
        frame = pd.read_csv(
            path,
            usecols=[tmc_column, "link_id", "from_node_id", "to_node_id"],
            low_memory=False,
        ).rename(columns={tmc_column: "tmc"})
        frame["facility_class"] = facility_class
        parts.append(frame)
    virtual_mapping_path = coverage_root / "virtual" / "virtual_tmc_to_link.csv"
    virtual_treatment_path = coverage_root / "virtual" / "virtual_link_treatments.csv"
    metadata = pd.DataFrame(
        columns=["tmc_code", "observation_source", "virtual_treatment", "virtual_confidence"]
    )
    if virtual_mapping_path.is_file() and virtual_treatment_path.is_file():
        virtual_mapping = pd.read_csv(virtual_mapping_path, low_memory=False)
        treatments = pd.read_csv(virtual_treatment_path, low_memory=False).rename(
            columns={
                "virtual_tmc": "tmc",
                "treatment": "virtual_treatment",
                "confidence": "virtual_confidence",
            }
        )
        keys = ["tmc", "link_id", "from_node_id", "to_node_id"]
        virtual_mapping = virtual_mapping.merge(
            treatments[
                keys + ["facility_class", "virtual_treatment", "virtual_confidence"]
            ],
            on=keys,
            how="left",
            validate="one_to_one",
        )
        parts.append(virtual_mapping[keys + ["facility_class"]])
        metadata = virtual_mapping[
            ["tmc", "virtual_treatment", "virtual_confidence"]
        ].rename(columns={"tmc": "tmc_code"})
        metadata["observation_source"] = "virtual"
    if not parts:
        raise FileNotFoundError(
            "coverage snapshot contains no GP/managed facility classification files"
        )
    facility = pd.concat(parts, ignore_index=True)
    key = ["tmc", "link_id", "from_node_id", "to_node_id"]
    for frame in (canonical, facility):
        frame["tmc"] = frame["tmc"].astype(str).str.strip()
        for column in key[1:]:
            frame[column] = pd.to_numeric(frame[column], errors="raise").astype(np.int64)
    facility = facility.drop_duplicates(key)
    classified = canonical.merge(facility, on=key, how="left", validate="one_to_one")
    if classified["facility_class"].isna().any():
        raise ValueError(
            "coverage sources do not classify {} canonical mappings".format(
                int(classified["facility_class"].isna().sum())
            )
        )
    return classified, metadata.drop_duplicates(["tmc_code", "observation_source"])


def _read_profiles(cbi_runs: Sequence[Tuple[Path, str]]) -> pd.DataFrame:
    wanted = {
        "tmc_code",
        "corridor",
        "t_min",
        "avg_weekday_speed_mph",
        "avg_weekday_flow_veh_per_hr_lane",
        "lanes",
        "n_days",
        "speed_at_capacity_mph",
    }
    frames: List[pd.DataFrame] = []
    for root, source in cbi_runs:
        for path in sorted(
            (root / "corridors").glob("*/03-profiles/average_weekday_profile.csv")
        ):
            header = set(pd.read_csv(path, nrows=0).columns)
            missing = sorted(wanted - header)
            if missing:
                raise ValueError("{} is missing columns: {}".format(path, missing))
            frame = pd.read_csv(path, usecols=sorted(wanted), low_memory=False)
            frame["observation_source"] = source
            frames.append(frame)
    if not frames:
        raise FileNotFoundError("no average-weekday CBI profiles were found")
    profiles = pd.concat(frames, ignore_index=True)
    profiles["tmc_code"] = profiles["tmc_code"].astype(str).str.strip()
    profiles["corridor"] = profiles["corridor"].astype(str).str.strip()
    for column in wanted - {"tmc_code", "corridor"}:
        profiles[column] = pd.to_numeric(profiles[column], errors="coerce")
    return profiles


def _read_episodes(
    cbi_runs: Sequence[Tuple[Path, str]], policy: str
) -> pd.DataFrame:
    wanted = {
        "tmc_code",
        "corridor",
        "period",
        "P_hr",
        "t0_hour",
        "t2_hour",
        "t3_hour",
        "min_speed_mph",
        "qdf",
        "plf",
        "episode_id",
    }
    frames: List[pd.DataFrame] = []
    for root, source in cbi_runs:
        for path in sorted(
            (root / "corridors").glob(
                "*/04-episode-detection/average_weekday_episode_candidates.csv"
            )
        ):
            header = set(pd.read_csv(path, nrows=0).columns)
            missing = sorted(wanted - header)
            if missing:
                raise ValueError("{} is missing columns: {}".format(path, missing))
            frame = pd.read_csv(path, usecols=sorted(wanted), low_memory=False)
            frame["observation_source"] = source
            frames.append(frame)
    if not frames:
        return pd.DataFrame(columns=list(wanted) + ["observation_source"])
    episodes = pd.concat(frames, ignore_index=True)
    episodes["tmc_code"] = episodes["tmc_code"].astype(str).str.strip()
    episodes["corridor"] = episodes["corridor"].astype(str).str.strip()
    episodes["period"] = episodes["period"].astype(str).str.strip().str.upper()
    numeric = ["P_hr", "t0_hour", "t2_hour", "t3_hour", "min_speed_mph", "qdf", "plf"]
    for column in numeric:
        episodes[column] = pd.to_numeric(episodes[column], errors="coerce")
    episodes = episodes[episodes["P_hr"].gt(0) & episodes["min_speed_mph"].gt(0)].copy()
    episodes["source_period"] = episodes["period"]
    episodes["raw_p_hr"] = episodes["P_hr"]
    if policy == "split_intersection":
        split_rows: List[Dict[str, object]] = []
        for row in episodes.to_dict("records"):
            start = float(row["t0_hour"])
            end = float(row["t3_hour"])
            if not (math.isfinite(start) and math.isfinite(end)):
                continue
            if end < start:
                end += 24.0
            for period, (start_minute, end_minute) in PERIOD_MINUTES.items():
                clip_start = max(start, start_minute / 60.0)
                clip_end = min(end, end_minute / 60.0)
                if clip_end <= clip_start:
                    continue
                split = dict(row)
                split.update(
                    {
                        "period": period,
                        "P_hr": clip_end - clip_start,
                        "clip_start_hour": clip_start,
                        "clip_end_hour": clip_end,
                    }
                )
                split_rows.append(split)
        episodes = pd.DataFrame(split_rows)
    else:
        episodes = episodes[episodes["period"].isin(PERIOD_MINUTES)].copy()
        episodes["clip_start_hour"] = np.nan
        episodes["clip_end_hour"] = np.nan
    if episodes.empty:
        return episodes
    episodes = episodes.sort_values(
        [
            "observation_source",
            "corridor",
            "tmc_code",
            "period",
            "min_speed_mph",
            "P_hr",
            "episode_id",
        ],
        ascending=[True, True, True, True, True, False, True],
        kind="mergesort",
    ).drop_duplicates(
        ["observation_source", "corridor", "tmc_code", "period"], keep="first"
    )
    return episodes


def build_tmc_targets(
    cbi_runs: Sequence[Tuple[Path, str]],
    episode_policy: str,
    observation_metadata: Optional[pd.DataFrame],
) -> pd.DataFrame:
    profiles = _read_profiles(cbi_runs)
    episodes = _read_episodes(cbi_runs, episode_policy)
    episode_index = {
        (row.observation_source, row.corridor, row.tmc_code, row.period): row
        for row in episodes.itertuples(index=False)
    }
    actual_days = pd.to_numeric(
        profiles.loc[profiles["observation_source"].eq("actual"), "n_days"],
        errors="coerce",
    )
    maximum_days = float(actual_days.max()) if not actual_days.empty else math.nan
    if not math.isfinite(maximum_days) or maximum_days <= 0:
        maximum_days = float(pd.to_numeric(profiles["n_days"], errors="coerce").max())
    if not math.isfinite(maximum_days) or maximum_days <= 0:
        maximum_days = 1.0
    rows: List[Dict[str, object]] = []
    for (source, corridor, tmc), group in profiles.groupby(
        ["observation_source", "corridor", "tmc_code"], sort=False
    ):
        for period, (start_minute, end_minute) in PERIOD_MINUTES.items():
            panel = group[group["t_min"].ge(start_minute) & group["t_min"].lt(end_minute)]
            speed = pd.to_numeric(panel["avg_weekday_speed_mph"], errors="coerce")
            if panel.empty or not speed.notna().any():
                continue
            flow = pd.to_numeric(
                panel["avg_weekday_flow_veh_per_hr_lane"], errors="coerce"
            )
            lanes = pd.to_numeric(panel["lanes"], errors="coerce").fillna(1).clip(lower=1)
            days = float(pd.to_numeric(panel["n_days"], errors="coerce").median())
            quality = (
                float(np.clip(days / maximum_days, 0, 1))
                if math.isfinite(days)
                else 0.5
            )
            episode = episode_index.get((source, corridor, tmc, period))
            vt2 = math.nan
            if episode is not None:
                vt2 = float(episode.min_speed_mph)
                if episode_policy == "split_intersection":
                    overlap = pd.to_numeric(
                        group.loc[
                            group["t_min"].ge(float(episode.clip_start_hour) * 60)
                            & group["t_min"].lt(float(episode.clip_end_hour) * 60),
                            "avg_weekday_speed_mph",
                        ],
                        errors="coerce",
                    )
                    if overlap.notna().any():
                        vt2 = float(overlap.min())
            plf = 1.0
            if episode is not None:
                supplied = float(episode.plf)
                if math.isfinite(supplied) and supplied > 0:
                    plf = supplied
                else:
                    qdf = float(episode.qdf)
                    hours = (end_minute - start_minute) / 60.0
                    if math.isfinite(qdf) and qdf > 0:
                        plf = 1.0 / (qdf * hours)
            cutoff = float(
                pd.to_numeric(panel["speed_at_capacity_mph"], errors="coerce").median()
            )
            rows.append(
                {
                    "corridor": corridor,
                    "tmc_code": tmc,
                    "observation_source": source,
                    "period": period,
                    "calibration_observation_class": "E" if episode is not None else "N",
                    "observed_p_hr": float(episode.P_hr) if episode is not None else np.nan,
                    "observed_vt2_mph": vt2,
                    "observed_avg_speed_mph": float(speed.mean()),
                    "s3_volume": float(np.nansum(flow.to_numpy() * lanes.to_numpy() * 0.25)),
                    "observation_quality": quality,
                    "mode1_plf_target": plf,
                    "cutoff_speed": cutoff,
                    "profile_n_days": days,
                    "profile_time_bins": int(speed.notna().sum()),
                }
            )
    targets = pd.DataFrame(rows)
    if targets.empty:
        raise ValueError("CBI profiles produced no AM/MD/PM targets")
    if observation_metadata is not None and not observation_metadata.empty:
        targets = targets.merge(
            observation_metadata,
            on=["tmc_code", "observation_source"],
            how="left",
            validate="many_to_one",
        )
    else:
        targets["virtual_treatment"] = ""
        targets["virtual_confidence"] = ""
    targets["virtual_treatment"] = targets["virtual_treatment"].fillna("")
    targets["virtual_confidence"] = targets["virtual_confidence"].fillna("")
    confidence = targets["virtual_confidence"].astype(str).str.lower().map(
        {"high": 1.0, "medium": 0.65, "low": 0.35}
    )
    virtual = targets["observation_source"].eq("virtual")
    if virtual.any() and confidence.loc[virtual].isna().any():
        raise ValueError("virtual targets require high, medium, or low confidence")
    targets.loc[virtual, "observation_quality"] = confidence.loc[virtual]
    return targets.sort_values(
        ["period", "observation_source", "corridor", "tmc_code"],
        kind="mergesort",
    )


def build_link_targets(
    targets: pd.DataFrame,
    mapping_path: Path,
    facility_path: Optional[Path],
    facility_frame: Optional[pd.DataFrame],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    mapping = pd.read_csv(mapping_path, low_memory=False)
    required = {
        "tmc",
        "link_id",
        "from_node_id",
        "to_node_id",
        "distance_to_tmc_ft",
        "road_order",
    }
    missing = sorted(required - set(mapping.columns))
    if missing:
        raise ValueError("{} is missing columns: {}".format(mapping_path, missing))
    mapping["tmc"] = mapping["tmc"].astype(str).str.strip()
    for column in ("link_id", "from_node_id", "to_node_id"):
        mapping[column] = pd.to_numeric(mapping[column], errors="raise").astype(
            np.int64
        )
    if "facility_class" not in mapping:
        if facility_frame is None:
            if facility_path is None:
                raise ValueError("facility_class is missing and no facility mapping was supplied")
            facility_frame = pd.read_csv(facility_path, low_memory=False)
        keys = ["tmc", "link_id", "from_node_id", "to_node_id"]
        facility = facility_frame[keys + ["facility_class"]].copy()
        facility["tmc"] = facility["tmc"].astype(str).str.strip()
        for column in ("link_id", "from_node_id", "to_node_id"):
            facility[column] = pd.to_numeric(
                facility[column], errors="raise"
            ).astype(np.int64)
        mapping = mapping.merge(
            facility.drop_duplicates(keys),
            on=keys,
            how="left",
            validate="one_to_one",
        )
    if mapping["facility_class"].isna().any():
        raise ValueError("facility mapping does not classify every canonical mapping")
    mapping["facility_class"] = mapping["facility_class"].astype(str).str.lower().str.strip()
    for column in ("distance_to_tmc_ft", "road_order"):
        mapping[column] = pd.to_numeric(mapping[column], errors="coerce")
    joined = mapping.merge(
        targets,
        left_on="tmc",
        right_on="tmc_code",
        how="inner",
        validate="many_to_many",
    )
    joined["_managed_rank"] = joined["facility_class"].ne("gp").astype(int)
    joined = joined.sort_values(
        ["period", "link_id", "_managed_rank", "distance_to_tmc_ft", "road_order", "tmc"],
        kind="mergesort",
        na_position="last",
    )
    joined["candidate_tmc_count"] = joined.groupby(
        ["period", "link_id"], sort=False
    )["tmc"].transform("size")
    joined["selected_for_link"] = ~joined.duplicated(["period", "link_id"])
    selected = joined[joined["selected_for_link"]].copy()
    selected["target_tmc"] = selected["tmc"]
    selected["calibration_exclusion_reason"] = np.where(
        selected["facility_class"].eq("gp"), "", "managed_facility"
    )
    selected.loc[
        selected["facility_class"].ne("gp"), "calibration_observation_class"
    ] = "U"
    return joined.drop(columns="_managed_rank"), selected.drop(columns="_managed_rank")


def _stage_scenario(source: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=False)
    excluded = {
        "link_performance.csv",
        "route_assignment.csv",
        "agent.csv",
        "summary_log_file.txt",
        "route_columns.bin",
        "TAP_log.csv",
    }
    for path in source.iterdir():
        if not path.is_file() or path.name in excluded or path.name.startswith("auto_calibration_"):
            continue
        target = destination / path.name
        if path.name in {"link.csv", "settings.csv", "departure_profiles.csv"}:
            shutil.copy2(path, target)
        else:
            try:
                os.link(path, target)
            except OSError:
                shutil.copy2(path, target)


def enrich_scenario(
    source_root: Path,
    destination_root: Path,
    selected: pd.DataFrame,
    config: PrepareAutoCalibrationConfig,
    processors: int,
) -> List[Dict[str, object]]:
    summaries: List[Dict[str, object]] = []
    destination_root.mkdir(parents=True, exist_ok=False)
    if (source_root / "nt").is_dir():
        _stage_scenario(source_root / "nt", destination_root / "nt")
    for period_lower in CALIBRATION_PERIODS:
        period = period_lower.upper()
        source = source_root / period_lower
        destination = destination_root / period_lower
        _stage_scenario(source, destination)
        if config.converted_network_root is not None:
            for filename in ("node.csv", "link.csv"):
                replacement = config.converted_network_root / period_lower / filename
                if not replacement.is_file():
                    raise FileNotFoundError(replacement)
                shutil.copy2(replacement, destination / filename)
        link_path = destination / "link.csv"
        links = pd.read_csv(link_path, low_memory=False)
        stale = [column for column in TARGET_COLUMNS if column in links]
        if stale:
            links = links.drop(columns=stale)
        period_targets = selected[selected["period"].eq(period)].copy()
        keep = [
            "link_id",
            "calibration_observation_class",
            "calibration_exclusion_reason",
            "facility_class",
            "target_tmc",
            "corridor",
            "observed_p_hr",
            "observed_vt2_mph",
            "observed_avg_speed_mph",
            "s3_volume",
            "observation_quality",
            "observation_source",
            "virtual_treatment",
            "virtual_confidence",
            "mode1_plf_target",
            "cutoff_speed",
        ]
        links = links.merge(
            period_targets[keep], on="link_id", how="left", validate="one_to_one"
        ).copy()
        links["facility_class"] = links["facility_class"].fillna("")
        links["calibration_observation_class"] = links[
            "calibration_observation_class"
        ].fillna("U")
        links["calibration_exclusion_reason"] = links[
            "calibration_exclusion_reason"
        ].fillna("unmapped_observation")
        for column in (
            "target_tmc",
            "corridor",
            "observation_source",
            "virtual_treatment",
            "virtual_confidence",
        ):
            links[column] = links[column].fillna("")
        links["mode1_plf"] = pd.to_numeric(links["vdf_plf"], errors="coerce")
        target_plf = pd.to_numeric(links.pop("mode1_plf_target"), errors="coerce")
        links.loc[target_plf.notna(), "mode1_plf"] = target_plf[target_plf.notna()]
        links["mode1_qcd"] = pd.to_numeric(links["vdf_cd"], errors="coerce")
        links["mode1_qcp"] = pd.to_numeric(links["vdf_cp"], errors="coerce")
        cube_column = "I4{}VOL".format(period)
        if cube_column not in links:
            raise ValueError("{} is missing {}".format(link_path, cube_column))
        links["cube_vehicle_volume"] = pd.to_numeric(
            links[cube_column], errors="coerce"
        )
        links["vdf_type"] = 2
        links["qvdf_profile_mode"] = 1
        links.to_csv(link_path, index=False)

        departure = config.departure_profile_csv
        if departure is None:
            departure = source / "departure_profiles.csv"
        if not departure.is_file():
            raise FileNotFoundError(departure)
        shutil.copy2(departure, destination / "departure_profiles.csv")
        calibration_settings = config.calibration_settings_csv
        if calibration_settings is None:
            calibration_settings = source / "auto_calibration_settings.csv"
        if not calibration_settings.is_file():
            raise FileNotFoundError(calibration_settings)
        settings_frame = pd.read_csv(calibration_settings, low_memory=False)
        if set(settings_frame.columns) != {"key", "value"}:
            raise ValueError("auto calibration settings must have key,value columns")
        # Settings values may mix numbers and empty strings. Convert to object
        # before updates so pandas 3's strict dtype assignment does not reject
        # valid string-valued settings.
        settings_frame["value"] = settings_frame["value"].astype(object)
        workers = settings_frame["key"].astype(str).str.strip().eq("workers")
        if workers.sum() != 1:
            raise ValueError("auto calibration settings require exactly one workers row")
        settings_frame.loc[workers, "value"] = str(processors)
        constraint = settings_frame["key"].astype(str).str.strip().eq(
            "volume_constraint_file"
        )
        if constraint.any():
            settings_frame.loc[constraint, "value"] = ""
        settings_frame.to_csv(destination / "auto_calibration_settings.csv", index=False)
        settings_path = destination / "settings.csv"
        settings = pd.read_csv(settings_path, low_memory=False)
        if len(settings) != 1:
            raise ValueError("settings.csv must contain exactly one row")
        overrides = {
            "number_of_iterations": 15,
            "number_of_processors": processors,
            "route_output": 0,
            "vehicle_output": 0,
            "auto_calibration": 1,
            "column_output": 2,
            "column_file_output": 0,
            "link_output": 1,
            "accessibility_output": 0,
            "convergence_gap_pct": 0.35,
            "auto_calibration_config": "auto_calibration_settings.csv",
            "departure_profile_file": "departure_profiles.csv",
        }
        for key, value in overrides.items():
            settings.loc[0, key] = value
        settings.to_csv(settings_path, index=False)
        counts = links["calibration_observation_class"].value_counts().to_dict()
        summaries.append(
            {
                "period": period,
                "links": int(len(links)),
                "episode_links": int(counts.get("E", 0)),
                "no_episode_links": int(counts.get("N", 0)),
                "excluded_links": int(counts.get("U", 0)),
            }
        )
    return summaries


def prepare_auto_targets(
    config: PrepareAutoCalibrationConfig,
    source_scenario_root: Path,
    destination_root: Path,
    audit_root: Path,
    processors: int,
) -> Dict[str, object]:
    cbi_runs, mapping, facility_frame, metadata = resolve_sources(config)
    for root, _source in cbi_runs:
        if not root.is_dir():
            raise FileNotFoundError(root)
    if not mapping.is_file():
        raise FileNotFoundError(mapping)
    if config.facility_mapping_csv is not None and not config.facility_mapping_csv.is_file():
        raise FileNotFoundError(config.facility_mapping_csv)
    targets = build_tmc_targets(cbi_runs, config.episode_period_policy, metadata)
    candidates, selected = build_link_targets(
        targets,
        mapping,
        config.facility_mapping_csv,
        facility_frame,
    )
    audit_root.mkdir(parents=True, exist_ok=True)
    targets.to_csv(audit_root / "tmc_period_targets.csv", index=False)
    candidates.to_csv(audit_root / "link_target_mapping_audit.csv", index=False)
    selected.to_csv(audit_root / "selected_link_period_targets.csv", index=False)
    summaries = enrich_scenario(
        source_scenario_root, destination_root, selected, config, processors
    )
    return {
        "status": "complete",
        "episode_period_policy": config.episode_period_policy,
        "observation_sources": [source for _path, source in cbi_runs],
        "tmc_period_targets": int(len(targets)),
        "selected_link_period_targets": int(len(selected)),
        "managed_lane_policy": "excluded and encoded as class U",
        "episode_source": "candidate episodes",
        "average_speed_policy": "full assignment-period arithmetic mean",
        "volume_policy": "S3 plus CUBE envelope; no sparse link counts",
        "periods": summaries,
    }
