"""Joint four-period ODME against observed daily screen totals.

The implementation preserves each period/mode's original positive OD support,
productions, attractions, and total demand. NT, AM, MD, and PM contribute to a
single daily screen objective; demand is never transferred between periods.
"""

from __future__ import annotations

import csv
import json
import logging
import math
import os
import shutil
import struct
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import sparse

from .config import OdmeConfig, PERIODS
from .paths import portable_path


LOGGER = logging.getLogger("taplite_calibration.odme")
DTAB_HEADER = struct.Struct("<4siq")
DTAB_RECORD_DTYPE = np.dtype(
    [("o_zone_id", "<i4"), ("d_zone_id", "<i4"), ("volume", "<f8")],
    align=False,
)


@dataclass
class JointPolicy:
    period: str
    name: str
    origin: np.ndarray
    destination: np.ndarray
    q0: np.ndarray
    screens: np.ndarray
    zones: np.ndarray
    matrix: sparse.csr_matrix
    production: np.ndarray
    attraction: np.ndarray
    lower: np.ndarray
    upper: np.ndarray

    @property
    def key(self) -> str:
        return "{}:{}".format(self.period, self.name)


def load_policy(path: Path) -> Tuple[
    np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, sparse.csr_matrix
]:
    with np.load(path, allow_pickle=False) as data:
        origin = data["origin"].astype(np.int32, copy=False)
        destination = data["destination"].astype(np.int32, copy=False)
        q0 = data["q0"].astype(np.float64)
        screens = data["screen_ids"].astype(np.int32, copy=False)
        zones = data["zone_external"].astype(np.int32, copy=False)
        matrix = sparse.csr_matrix(
            (
                data["data"].astype(np.float64, copy=False),
                data["indices"].astype(np.int32, copy=False),
                data["indptr"].astype(np.int64, copy=False),
            ),
            shape=(q0.size, screens.size),
        )
    return origin, destination, q0, screens, zones, matrix


def marginal_error(
    q: np.ndarray,
    origin: np.ndarray,
    destination: np.ndarray,
    production: np.ndarray,
    attraction: np.ndarray,
) -> Tuple[float, float, float, float]:
    row = np.bincount(origin, weights=q, minlength=production.size)
    column = np.bincount(destination, weights=q, minlength=attraction.size)
    row_abs = float(np.max(np.abs(row - production)))
    col_abs = float(np.max(np.abs(column - attraction)))
    row_rel = float(np.max(np.abs(row - production) / np.maximum(production, 1.0)))
    col_rel = float(np.max(np.abs(column - attraction) / np.maximum(attraction, 1.0)))
    return row_abs, col_abs, row_rel, col_rel


def project_box_marginals(
    proposal: np.ndarray,
    policy: JointPolicy,
    max_iterations: int,
    tolerance: float,
) -> Tuple[np.ndarray, int, float]:
    """KL-project a sparse OD vector onto box and marginal constraints."""
    q = np.clip(np.asarray(proposal, dtype=np.float64), policy.lower, policy.upper)
    n_zones = policy.production.size
    active_row = policy.production > 0
    active_col = policy.attraction > 0
    error = math.inf
    for iteration in range(1, max_iterations + 1):
        row = np.bincount(policy.origin, weights=q, minlength=n_zones)
        if np.any(row[active_row] <= 0):
            raise ValueError("bounded projection lost a positive production row")
        row_scale = np.ones(n_zones, dtype=np.float64)
        row_scale[active_row] = policy.production[active_row] / row[active_row]
        q = np.clip(q * row_scale[policy.origin], policy.lower, policy.upper)

        column = np.bincount(policy.destination, weights=q, minlength=n_zones)
        if np.any(column[active_col] <= 0):
            raise ValueError("bounded projection lost a positive attraction column")
        col_scale = np.ones(n_zones, dtype=np.float64)
        col_scale[active_col] = policy.attraction[active_col] / column[active_col]
        q = np.clip(q * col_scale[policy.destination], policy.lower, policy.upper)

        if iteration == 1 or iteration % 5 == 0:
            _, _, row_rel, col_rel = marginal_error(
                q,
                policy.origin,
                policy.destination,
                policy.production,
                policy.attraction,
            )
            error = max(row_rel, col_rel)
            if error <= tolerance:
                break
    if error > tolerance:
        raise RuntimeError(
            "bounded marginal projection did not converge: relative error={:.3e}".format(
                error
            )
        )
    return q, iteration, error


def period_total_variation(
    policies: Sequence[JointPolicy], adjusted: Mapping[str, np.ndarray]
) -> float:
    demand = float(sum(policy.q0.sum() for policy in policies))
    if demand <= 0:
        return 0.0
    moved = 0.5 * sum(
        float(np.abs(adjusted[policy.key] - policy.q0).sum()) for policy in policies
    )
    return moved / demand


def enforce_period_tv_budget(
    policies: Sequence[JointPolicy],
    adjusted: Mapping[str, np.ndarray],
    budget: Optional[float],
) -> Tuple[Dict[str, np.ndarray], float, float]:
    tv = period_total_variation(policies, adjusted)
    if budget is None or tv <= budget + 1e-15:
        return dict(adjusted), tv, 1.0
    scale = budget / tv if tv > 0 else 1.0
    result = {
        policy.key: policy.q0 + scale * (adjusted[policy.key] - policy.q0)
        for policy in policies
    }
    return result, period_total_variation(policies, result), scale


def screen_totals(
    policies: Sequence[JointPolicy], adjusted: Mapping[str, np.ndarray]
) -> np.ndarray:
    result = np.zeros(policies[0].matrix.shape[1], dtype=np.float64)
    for policy in policies:
        result += np.asarray(policy.matrix.T @ adjusted[policy.key]).ravel()
    return result


def objective(
    policies: Sequence[JointPolicy],
    adjusted: Mapping[str, np.ndarray],
    screen: np.ndarray,
    target: np.ndarray,
    kl_weight: float,
    total_demand: float,
) -> Tuple[float, float, float]:
    relative = (screen - target) / target
    screen_loss = 0.5 * float(np.mean(relative ** 2))
    kl = 0.0
    for policy in policies:
        q = adjusted[policy.key]
        kl += float(
            np.sum(q * np.log(np.maximum(q, 1e-300) / policy.q0) - q + policy.q0)
        )
    normalized_kl = kl / total_demand
    return screen_loss + kl_weight * normalized_kl, screen_loss, normalized_kl


def _copy_scenario(
    source: Path, destination: Path, mutable_names: Sequence[str]
) -> None:
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
        if path.is_file() and path.name not in excluded:
            target = destination / path.name
            if path.name in mutable_names:
                shutil.copy2(path, target)
            else:
                try:
                    os.link(path, target)
                except OSError:
                    shutil.copy2(path, target)


def _dense_replacement(policy: JointPolicy, q: np.ndarray) -> np.ndarray:
    dense = np.full(policy.zones.size * policy.zones.size, np.nan, dtype=np.float32)
    key = policy.origin.astype(np.int64) * policy.zones.size + policy.destination
    dense[key] = q.astype(np.float32)
    return dense


def _external_to_internal(values: np.ndarray, zones: np.ndarray) -> np.ndarray:
    lookup = {int(zone): index for index, zone in enumerate(zones)}
    return np.array([lookup.get(int(value), -1) for value in values], dtype=np.int32)


def _write_adjusted_demands(
    source: Path,
    destination: Path,
    policies: Sequence[JointPolicy],
    adjusted: Mapping[str, np.ndarray],
    modes: pd.DataFrame,
) -> List[Dict[str, object]]:
    by_name = {policy.name: policy for policy in policies}
    audit: List[Dict[str, object]] = []
    for row in modes.to_dict("records"):
        name = str(row["mode_type"])
        demand_file = str(row["demand_file"])
        policy = by_name[name]
        q = adjusted[policy.key]
        dense = _dense_replacement(policy, q)
        csv_source = source / demand_file
        binary_source = csv_source.with_suffix(".bin")
        if csv_source.is_file():
            frame = pd.read_csv(csv_source, low_memory=False)
            o_external = pd.to_numeric(frame["o_zone_id"], errors="raise").to_numpy(int)
            d_external = pd.to_numeric(frame["d_zone_id"], errors="raise").to_numpy(int)
            o_internal = _external_to_internal(o_external, policy.zones)
            d_internal = _external_to_internal(d_external, policy.zones)
            use = (o_internal >= 0) & (d_internal >= 0) & (o_internal != d_internal)
            row_key = o_internal[use].astype(np.int64) * policy.zones.size + d_internal[use]
            replacement = dense[row_key]
            has = np.isfinite(replacement)
            selected = np.flatnonzero(use)[has]
            # pandas 3 may expose a read-only NumPy view. ODME intentionally
            # mutates this independent output buffer before writing the staged
            # demand file, so request an explicit writable copy.
            volume = (
                pd.to_numeric(frame["volume"], errors="coerce")
                .fillna(0)
                .to_numpy(dtype=float, copy=True)
            )
            volume[selected] = replacement[has]
            frame["volume"] = volume
            frame.to_csv(destination / demand_file, index=False)
            replaced = int(has.sum())
            storage = "csv"
        elif binary_source.is_file():
            with binary_source.open("rb") as stream:
                raw = stream.read(DTAB_HEADER.size)
            if len(raw) != DTAB_HEADER.size:
                raise ValueError("truncated DTAB header: {}".format(binary_source))
            magic, version, count = DTAB_HEADER.unpack(raw)
            if magic != b"DTAB" or version != 1 or count < 0:
                raise ValueError("unsupported DTAB file: {}".format(binary_source))
            records = np.memmap(
                binary_source,
                dtype=DTAB_RECORD_DTYPE,
                mode="r",
                offset=DTAB_HEADER.size,
                shape=(count,),
            )
            output_records = np.array(records, copy=True)
            o_internal = _external_to_internal(output_records["o_zone_id"], policy.zones)
            d_internal = _external_to_internal(output_records["d_zone_id"], policy.zones)
            use = (o_internal >= 0) & (d_internal >= 0) & (o_internal != d_internal)
            row_key = o_internal[use].astype(np.int64) * policy.zones.size + d_internal[use]
            replacement = dense[row_key]
            has = np.isfinite(replacement)
            selected = np.flatnonzero(use)[has]
            output_records["volume"][selected] = replacement[has]
            destination_binary = destination / binary_source.name
            with destination_binary.open("wb") as stream:
                stream.write(DTAB_HEADER.pack(b"DTAB", 1, len(output_records)))
                output_records.tofile(stream)
            replaced = int(has.sum())
            storage = "dtab"
        else:
            raise FileNotFoundError("missing demand input for {}".format(policy.key))
        if replaced != q.size:
            raise RuntimeError(
                "{}: replaced {:,} records but policy has {:,} cells".format(
                    policy.key, replaced, q.size
                )
            )
        audit.append(
            {
                "period": policy.period,
                "mode": name,
                "storage": storage,
                "positive_interzonal_records_replaced": replaced,
            }
        )
    return audit


def _write_factor_shard(
    output: Path, policy: JointPolicy, q: np.ndarray
) -> Dict[str, object]:
    destination = output / "factor_shards"
    destination.mkdir(parents=True, exist_ok=True)
    filename = "{}_{}_od_factors.npz".format(policy.period, policy.name)
    shard = destination / filename
    np.savez_compressed(
        shard,
        origin_zone_id=policy.zones[policy.origin].astype(np.int32, copy=False),
        destination_zone_id=policy.zones[policy.destination].astype(np.int32, copy=False),
        factor=(q / policy.q0).astype(np.float32),
        original_volume=policy.q0.astype(np.float32),
        adjusted_volume=q.astype(np.float32),
    )
    return {
        "file": "factor_shards/{}".format(filename),
        "positive_od_cells": int(q.size),
        "default_factor_for_absent_or_intrazonal_od": 1.0,
    }


def run_odme(
    config: OdmeConfig,
    output: Path,
    scenario_output: Path,
    project_dir: Path,
) -> Dict[str, object]:
    started = time.time()
    output.mkdir(parents=True, exist_ok=False)
    policies: List[JointPolicy] = []
    policies_by_period: Dict[str, List[JointPolicy]] = {}
    modes_by_period: Dict[str, pd.DataFrame] = {}
    common_screens: Optional[np.ndarray] = None
    common_zones: Optional[np.ndarray] = None

    for period in PERIODS:
        source = config.scenario_root / period
        modes = pd.read_csv(source / "mode_type.csv", low_memory=False)
        modes_by_period[period] = modes
        period_policies: List[JointPolicy] = []
        for row in modes.to_dict("records"):
            name = str(row["mode_type"])
            path = config.policy_root / period / "od_screen_policy_{}.npz".format(name)
            origin, destination, q0, screens, zones, matrix = load_policy(path)
            if np.any(q0 <= 0):
                raise ValueError("{} contains non-positive supported OD cells".format(path))
            if common_screens is None:
                common_screens = screens
                common_zones = zones
            elif not np.array_equal(screens, common_screens) or not np.array_equal(
                zones, common_zones
            ):
                raise ValueError("all periods and modes must share screen and zone axes")
            production = np.bincount(origin, weights=q0, minlength=zones.size)
            attraction = np.bincount(destination, weights=q0, minlength=zones.size)
            policy = JointPolicy(
                period=period,
                name=name,
                origin=origin,
                destination=destination,
                q0=q0,
                screens=screens,
                zones=zones,
                matrix=matrix,
                production=production,
                attraction=attraction,
                lower=(1.0 - config.factor_cap) * q0,
                upper=(1.0 + config.factor_cap) * q0,
            )
            policies.append(policy)
            period_policies.append(policy)
        policies_by_period[period] = period_policies
        LOGGER.info(
            "Loaded %s: %d modes, %d supported OD cells",
            period.upper(),
            len(period_policies),
            sum(policy.q0.size for policy in period_policies),
        )

    assert common_screens is not None and common_zones is not None
    target_frame = pd.read_csv(config.daily_target_csv)
    if "screen_id" not in target_frame or config.daily_target_column not in target_frame:
        raise ValueError("daily target CSV is missing screen_id or target column")
    target_frame["screen_id"] = pd.to_numeric(target_frame["screen_id"], errors="raise").astype(int)
    target_frame = target_frame.set_index("screen_id")
    missing_screens = sorted(set(common_screens.tolist()) - set(target_frame.index.tolist()))
    if missing_screens:
        raise ValueError("daily target CSV is missing screens: {}".format(missing_screens))
    target = pd.to_numeric(
        target_frame.loc[common_screens, config.daily_target_column], errors="raise"
    ).to_numpy(float)
    if np.any(target <= 0):
        raise ValueError("daily screen targets must be positive")

    adjusted = {policy.key: policy.q0.copy() for policy in policies}
    baseline_by_period = {
        period: screen_totals(period_policies, adjusted)
        for period, period_policies in policies_by_period.items()
    }
    baseline = sum(baseline_by_period.values(), np.zeros_like(target))
    total_demand = float(sum(policy.q0.sum() for policy in policies))
    screen = screen_totals(policies, adjusted)
    loss, screen_loss, normalized_kl = objective(
        policies, adjusted, screen, target, config.kl_weight, total_demand
    )
    history: List[Dict[str, object]] = []
    projection_rows: List[Dict[str, object]] = []

    def history_row(iteration: int, accepted_step: float) -> Dict[str, object]:
        row: Dict[str, object] = {
            "iteration": iteration,
            "objective": loss,
            "relative_screen_loss": screen_loss,
            "normalized_kl": normalized_kl,
            "screen_wape_pct": float(100 * np.abs(screen - target).sum() / target.sum()),
            "screen_max_abs_pct": float(100 * np.max(np.abs((screen - target) / target))),
            "accepted_step": accepted_step,
        }
        for period, period_policies in policies_by_period.items():
            row["{}_tv_share".format(period)] = period_total_variation(
                period_policies, adjusted
            )
        return row

    history.append(history_row(0, 0.0))
    step = config.initial_step
    for iteration in range(1, config.max_iterations + 1):
        screen_gradient = (screen - target) / (target * target * target.size)
        gradients: Dict[str, np.ndarray] = {}
        gradient_scale = 0.0
        for policy in policies:
            q = adjusted[policy.key]
            gradient = np.asarray(policy.matrix @ screen_gradient).ravel()
            if config.kl_weight > 0:
                gradient += (
                    config.kl_weight
                    / total_demand
                    * np.log(np.maximum(q, 1e-300) / policy.q0)
                )
            gradients[policy.key] = gradient
            gradient_scale = max(gradient_scale, float(np.max(np.abs(gradient))))
        if gradient_scale <= 1e-18:
            LOGGER.info("ODME gradient is numerically zero; stopping")
            break

        accepted = False
        trial_step = step
        for _line_search in range(config.max_line_search):
            candidate: Dict[str, np.ndarray] = {}
            candidate_projection: List[Dict[str, object]] = []
            for policy in policies:
                proposal = adjusted[policy.key] * np.exp(
                    np.clip(
                        -trial_step * gradients[policy.key] / gradient_scale,
                        -0.7,
                        0.7,
                    )
                )
                q, sweeps, error = project_box_marginals(
                    proposal,
                    policy,
                    config.projection_iterations,
                    config.projection_tolerance,
                )
                candidate[policy.key] = q
                candidate_projection.append(
                    {
                        "stage": "mirror_iteration_{}".format(iteration),
                        "period": policy.period,
                        "mode": policy.name,
                        "iterations": sweeps,
                        "relative_marginal_error": error,
                    }
                )
            for period, period_policies in policies_by_period.items():
                subset = {policy.key: candidate[policy.key] for policy in period_policies}
                constrained, tv, scale = enforce_period_tv_budget(
                    period_policies, subset, config.tv_budget
                )
                candidate.update(constrained)
                candidate_projection.append(
                    {
                        "stage": "mirror_iteration_{}_tv".format(iteration),
                        "period": period,
                        "mode": "ALL",
                        "iterations": 1,
                        "relative_marginal_error": 0.0,
                        "tv_share": tv,
                        "tv_retraction_scale": scale,
                    }
                )
            candidate_screen = screen_totals(policies, candidate)
            candidate_values = objective(
                policies,
                candidate,
                candidate_screen,
                target,
                config.kl_weight,
                total_demand,
            )
            if candidate_values[0] < loss - config.minimum_improvement:
                adjusted = candidate
                screen = candidate_screen
                loss, screen_loss, normalized_kl = candidate_values
                projection_rows.extend(candidate_projection)
                accepted = True
                break
            trial_step *= 0.5
            if trial_step < config.minimum_step:
                break
        if not accepted:
            LOGGER.info("ODME iteration %d found no improving feasible step", iteration)
            break
        step = min(config.maximum_step, trial_step * 1.35)
        row = history_row(iteration, trial_step)
        history.append(row)
        LOGGER.info(
            "ODME iteration %d: WAPE %.4f%%, objective %.8g",
            iteration,
            row["screen_wape_pct"],
            loss,
        )

    pd.DataFrame(history).to_csv(output / "optimization_history.csv", index=False)
    pd.DataFrame(projection_rows).to_csv(output / "projection_audit.csv", index=False)

    factor_index: Dict[str, object] = {
        "schema_version": 1,
        "factor_cap": config.factor_cap,
        "tv_budget_per_period": config.tv_budget,
        "default_factor": 1.0,
        "path_base": "directory containing od_factor_dictionary.npy",
        "coverage": "original positive interzonal OD support",
        "periods": {},
    }
    mode_rows: List[Dict[str, object]] = []
    for policy in policies:
        q = adjusted[policy.key]
        factor = q / policy.q0
        row_abs, col_abs, row_rel, col_rel = marginal_error(
            q,
            policy.origin,
            policy.destination,
            policy.production,
            policy.attraction,
        )
        mode_rows.append(
            {
                "period": policy.period,
                "mode": policy.name,
                "vehicle_demand": float(policy.q0.sum()),
                "total_variation_share": float(
                    0.5 * np.abs(q - policy.q0).sum() / policy.q0.sum()
                ),
                "kl_divergence": float(
                    np.sum(q * np.log(np.maximum(q, 1e-300) / policy.q0) - q + policy.q0)
                ),
                "max_production_absolute_error": row_abs,
                "max_attraction_absolute_error": col_abs,
                "max_production_relative_error": row_rel,
                "max_attraction_relative_error": col_rel,
                "factor_min": float(factor.min()),
                "factor_p50": float(np.quantile(factor, 0.50)),
                "factor_max": float(factor.max()),
            }
        )
        period_entry = factor_index["periods"].setdefault(policy.period, {"modes": {}})  # type: ignore
        period_entry["modes"][policy.name] = _write_factor_shard(output, policy, q)
        adjusted_dir = output / "adjusted_q" / policy.period
        adjusted_dir.mkdir(parents=True, exist_ok=True)
        np.save(adjusted_dir / "adjusted_q_{}.npy".format(policy.name), q.astype(np.float32))
    pd.DataFrame(mode_rows).to_csv(output / "od_adjustment_by_period_mode.csv", index=False)
    np.save(output / "od_factor_dictionary.npy", factor_index, allow_pickle=True)

    period_rows: List[pd.DataFrame] = []
    for period, period_policies in policies_by_period.items():
        period_rows.append(
            pd.DataFrame(
                {
                    "screen_id": common_screens,
                    "period": period.upper(),
                    "baseline_fixed_policy_vehicle_volume": baseline_by_period[period],
                    "joint_odme_fixed_policy_vehicle_volume": screen_totals(
                        period_policies, adjusted
                    ),
                }
            )
        )
    pd.concat(period_rows, ignore_index=True).to_csv(
        output / "screen_period_contributions_no_period_targets.csv", index=False
    )

    comparison = pd.DataFrame(
        {
            "screen_id": common_screens,
            "observed_daily_vehicle_volume": target,
            "baseline_daily_fixed_policy_vehicle_volume": baseline,
            "joint_odme_daily_fixed_policy_vehicle_volume": screen,
        }
    )
    for prefix in ("baseline_daily_fixed_policy", "joint_odme_daily_fixed_policy"):
        comparison["{}_error".format(prefix)] = (
            comparison["{}_vehicle_volume".format(prefix)]
            - comparison["observed_daily_vehicle_volume"]
        )
        comparison["{}_pct_error".format(prefix)] = (
            100
            * comparison["{}_error".format(prefix)]
            / comparison["observed_daily_vehicle_volume"]
        )
    comparison.to_csv(output / "screen_joint_daily_fixed_policy.csv", index=False)

    scenario_output.mkdir(parents=True, exist_ok=False)
    demand_audit: List[Dict[str, object]] = []
    for period, period_policies in policies_by_period.items():
        source = config.scenario_root / period
        destination = scenario_output / period
        demand_names: List[str] = []
        for demand_file in modes_by_period[period]["demand_file"].astype(str):
            demand_names.extend(
                [demand_file, str(Path(demand_file).with_suffix(".bin"))]
            )
        _copy_scenario(source, destination, demand_names)
        demand_audit.extend(
            _write_adjusted_demands(
                source,
                destination,
                period_policies,
                adjusted,
                modes_by_period[period],
            )
        )
    pd.DataFrame(demand_audit).to_csv(output / "demand_write_audit.csv", index=False)

    manifest: Dict[str, object] = {
        "status": "complete",
        "stage": "odme",
        "target_semantics": "observed daily screens only; no CUBE or Vs3 period targets",
        "periods": list(PERIODS),
        "factor_cap": config.factor_cap,
        "factor_range": [1.0 - config.factor_cap, 1.0 + config.factor_cap],
        "kl_weight": config.kl_weight,
        "tv_budget_per_period": config.tv_budget,
        "screen_count": int(target.size),
        "fixed_policy_wape_before_pct": float(
            100 * np.abs(baseline - target).sum() / target.sum()
        ),
        "fixed_policy_wape_after_pct": float(
            100 * np.abs(screen - target).sum() / target.sum()
        ),
        "optimization_iterations_accepted": len(history) - 1,
        "total_interzonal_vehicle_demand": total_demand,
        "constraints": {
            "period_mode_productions_preserved": True,
            "period_mode_attractions_preserved": True,
            "period_mode_total_demand_preserved": True,
            "cross_period_demand_transfer": False,
            "positive_support_preserved": True,
        },
        "od_factor_dictionary": portable_path(
            output / "od_factor_dictionary.npy", project_dir
        ),
        "adjusted_scenario_root": portable_path(scenario_output, project_dir),
        "elapsed_seconds": time.time() - started,
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest
