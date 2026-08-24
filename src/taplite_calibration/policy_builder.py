"""Streaming DTAC-v2 route-pool to OD/screen policy preparation."""

from __future__ import annotations

import logging
import struct
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO, Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd


LOGGER = logging.getLogger("taplite_calibration.prepare.policy")
DTAC_MAGIC = 0x43415444
DTAB_HEADER = struct.Struct("<4siq")
DTAB_RECORD_DTYPE = np.dtype(
    [("o_zone_id", "<i4"), ("d_zone_id", "<i4"), ("volume", "<f8")],
    align=False,
)


def _read_array(stream: BinaryIO, dtype: str, count: int) -> np.ndarray:
    if count < 0:
        raise ValueError("negative DTAC array length {}".format(count))
    array = np.fromfile(stream, dtype=np.dtype(dtype), count=count)
    if array.size != count:
        raise ValueError(
            "truncated DTAC array: wanted {}, read {}".format(count, array.size)
        )
    return array


@dataclass(frozen=True)
class DtacHeader:
    version: int
    n_modes: int
    n_zones: int
    fingerprint: Tuple[float, float, float, float]


def read_dtac_header(stream: BinaryIO) -> DtacHeader:
    raw = stream.read(16)
    if len(raw) != 16:
        raise ValueError("truncated DTAC header")
    magic, version, n_modes, n_zones = struct.unpack("<4i", raw)
    if magic != DTAC_MAGIC or version != 2:
        raise ValueError(
            "expected DTAC v2, got magic=0x{:08X}, version={}".format(
                magic, version
            )
        )
    fingerprint_raw = stream.read(16)
    if len(fingerprint_raw) != 16:
        raise ValueError("truncated DTAC fingerprint")
    fingerprint = struct.unpack("<4f", fingerprint_raw)
    return DtacHeader(version, n_modes, n_zones, fingerprint)


def read_dtac_block(stream: BinaryIO) -> Tuple[np.ndarray, ...]:
    raw = stream.read(4)
    if len(raw) != 4:
        raise ValueError("truncated DTAC origin block")
    (n_dest,) = struct.unpack("<i", raw)
    if n_dest < 0:
        raise ValueError("negative DTAC destination count")
    if n_dest == 0:
        empty_i = np.empty(0, dtype="<i4")
        empty_f = np.empty(0, dtype="<f4")
        return (
            empty_i,
            np.zeros(1, dtype="<i4"),
            empty_f,
            np.zeros(1, dtype="<i4"),
            empty_i,
        )
    destinations = _read_array(stream, "<i4", n_dest)
    path_offsets = _read_array(stream, "<i4", n_dest + 1)
    theta = _read_array(stream, "<f4", int(path_offsets[-1]))
    link_offsets = _read_array(stream, "<i4", theta.size + 1)
    links = _read_array(stream, "<i4", int(link_offsets[-1]))
    return destinations, path_offsets, theta, link_offsets, links


def zone_sequence(node_path: Path) -> np.ndarray:
    nodes = pd.read_csv(node_path, usecols=["zone_id"], low_memory=False)
    zone = pd.to_numeric(nodes["zone_id"], errors="coerce").fillna(0).astype(np.int64)
    values = zone[zone.ge(1)].to_numpy(dtype=np.int32)
    if len(np.unique(values)) != len(values):
        raise ValueError("node.csv has duplicate positive zone_id values")
    return values


def read_modes(run_dir: Path) -> List[Dict[str, object]]:
    modes = pd.read_csv(run_dir / "mode_type.csv", low_memory=False)
    if not {"mode_type", "demand_file"}.issubset(modes.columns):
        raise ValueError("mode_type.csv requires mode_type and demand_file")
    return [
        {
            "name": str(row["mode_type"]).strip(),
            "demand_file": str(row["demand_file"]).strip(),
        }
        for row in modes.to_dict("records")
    ]


class DemandIndex:
    """Compact, sorted vehicle-demand lookup for one mode."""

    def __init__(self, path: Path):
        if path.is_file():
            frame = pd.read_csv(
                path,
                usecols=["o_zone_id", "d_zone_id", "volume"],
                dtype={
                    "o_zone_id": np.int32,
                    "d_zone_id": np.int32,
                    "volume": np.float32,
                },
                low_memory=False,
            )
            frame = frame.loc[
                frame["volume"].gt(0) & frame["o_zone_id"].ne(frame["d_zone_id"])
            ].sort_values(["o_zone_id", "d_zone_id"], kind="mergesort")
            if frame.duplicated(["o_zone_id", "d_zone_id"]).any():
                frame = frame.groupby(
                    ["o_zone_id", "d_zone_id"], as_index=False
                )["volume"].sum()
            self.origin = frame["o_zone_id"].to_numpy(dtype=np.int32, copy=False)
            self.destination = frame["d_zone_id"].to_numpy(dtype=np.int32, copy=False)
            self.volume = frame["volume"].to_numpy(dtype=np.float64, copy=False)
        else:
            binary = path.with_suffix(".bin")
            if not binary.is_file():
                raise FileNotFoundError(
                    "missing demand CSV and DTAB binary: {}".format(path)
                )
            with binary.open("rb") as stream:
                raw = stream.read(DTAB_HEADER.size)
            if len(raw) != DTAB_HEADER.size:
                raise ValueError("truncated DTAB header: {}".format(binary))
            magic, version, count = DTAB_HEADER.unpack(raw)
            if magic != b"DTAB" or version != 1 or count < 0:
                raise ValueError("unsupported DTAB demand file: {}".format(binary))
            expected = DTAB_HEADER.size + count * DTAB_RECORD_DTYPE.itemsize
            if binary.stat().st_size != expected:
                raise ValueError("DTAB size mismatch: {}".format(binary))
            records = np.memmap(
                binary,
                dtype=DTAB_RECORD_DTYPE,
                mode="r",
                offset=DTAB_HEADER.size,
                shape=(count,),
            )
            mask = (records["volume"] > 0) & (
                records["o_zone_id"] != records["d_zone_id"]
            )
            self.origin = np.asarray(records["o_zone_id"][mask], dtype=np.int32)
            self.destination = np.asarray(records["d_zone_id"][mask], dtype=np.int32)
            self.volume = np.asarray(records["volume"][mask], dtype=np.float64)
            order = np.lexsort((self.destination, self.origin))
            self.origin = self.origin[order]
            self.destination = self.destination[order]
            self.volume = self.volume[order]
            if self.origin.size > 1:
                duplicate = (self.origin[1:] == self.origin[:-1]) & (
                    self.destination[1:] == self.destination[:-1]
                )
                if np.any(duplicate):
                    starts = np.r_[0, np.flatnonzero(~duplicate) + 1]
                    self.volume = np.add.reduceat(self.volume, starts)
                    self.origin = self.origin[starts]
                    self.destination = self.destination[starts]
        origins, first, counts = np.unique(
            self.origin, return_index=True, return_counts=True
        )
        self.slices = {
            int(origin): (int(index), int(index + count))
            for origin, index, count in zip(origins, first, counts)
        }

    def match(self, origin: int, destinations: np.ndarray) -> np.ndarray:
        result = np.zeros(destinations.size, dtype=np.float64)
        limits = self.slices.get(int(origin))
        if limits is None or destinations.size == 0:
            return result
        start, end = limits
        available = self.destination[start:end]
        volume = self.volume[start:end]
        positions = np.searchsorted(available, destinations)
        valid = positions < available.size
        valid_index = np.flatnonzero(valid)
        exact = available[positions[valid_index]] == destinations[valid_index]
        chosen = valid_index[exact]
        result[chosen] = volume[positions[chosen]]
        return result


def _screen_lookup(link_path: Path, screens: np.ndarray) -> np.ndarray:
    links = pd.read_csv(link_path, usecols=["link_id", "SCREEN"], low_memory=False)
    link_ids = pd.to_numeric(links["link_id"], errors="raise").astype(np.int64)
    raw_screen = pd.to_numeric(links["SCREEN"], errors="coerce").fillna(0).astype(int)
    maximum = int(link_ids.max())
    lookup = np.full(maximum + 1, -1, dtype=np.int32)
    screen_index = {int(screen): index for index, screen in enumerate(screens)}
    lookup[link_ids.to_numpy()] = raw_screen.map(screen_index).fillna(-1).to_numpy(
        dtype=np.int32
    )
    return lookup


def _path_screen_pairs(
    links: np.ndarray,
    link_offsets: np.ndarray,
    screen_lookup: np.ndarray,
    screen_count: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    if links.size == 0:
        empty = np.empty(0, dtype=np.int64)
        return empty, empty.astype(np.int32), empty.astype(np.int16)
    if int(links.max(initial=0)) >= screen_lookup.size or int(links.min(initial=0)) < 0:
        raise ValueError("DTAC contains a link_id absent from link.csv")
    membership = screen_lookup[links]
    positions = np.flatnonzero(membership >= 0)
    if positions.size == 0:
        empty = np.empty(0, dtype=np.int64)
        return empty, empty.astype(np.int32), empty.astype(np.int16)
    path_index = np.searchsorted(link_offsets[1:], positions, side="right")
    codes = path_index.astype(np.int64) * screen_count + membership[positions]
    unique, counts = np.unique(codes, return_counts=True)
    return (
        unique // screen_count,
        (unique % screen_count).astype(np.int32),
        counts.astype(np.int16),
    )


def build_period_policies(
    run_dir: Path,
    output: Path,
    screens: np.ndarray,
    period: str,
) -> Dict[str, object]:
    output.mkdir(parents=True, exist_ok=False)
    route_path = run_dir / "route_columns.bin"
    if not route_path.is_file():
        raise FileNotFoundError(route_path)
    zones = zone_sequence(run_dir / "node.csv")
    modes = read_modes(run_dir)
    screen_lookup = _screen_lookup(run_dir / "link.csv", screens)
    mode_stats: List[Dict[str, object]] = []
    with route_path.open("rb") as stream:
        header = read_dtac_header(stream)
        if header.n_modes != len(modes) or header.n_zones != len(zones):
            raise ValueError("DTAC header does not match mode or zone inputs")
        for mode in modes:
            demand = DemandIndex(run_dir / str(mode["demand_file"]))
            policy_origin: List[np.ndarray] = []
            policy_destination: List[np.ndarray] = []
            policy_q: List[np.ndarray] = []
            policy_indptr: List[np.ndarray] = [np.array([0], dtype=np.int64)]
            policy_indices: List[np.ndarray] = []
            policy_data: List[np.ndarray] = []
            policy_nnz = 0
            external_to_internal = {
                int(zone): index for index, zone in enumerate(zones)
            }
            theta_error = 0.0
            for origin_index, external_origin in enumerate(zones):
                destinations, path_offsets, theta, link_offsets, path_links = (
                    read_dtac_block(stream)
                )
                if destinations.size == 0:
                    continue
                q = demand.match(int(external_origin), destinations)
                path_od = np.repeat(
                    np.arange(destinations.size, dtype=np.int32),
                    np.diff(path_offsets),
                )
                theta_sum = np.add.reduceat(theta.astype(np.float64), path_offsets[:-1])
                theta_error = max(
                    theta_error,
                    float(np.max(np.abs(theta_sum - 1.0), initial=0.0)),
                )
                path_index, screen_index, count = _path_screen_pairs(
                    path_links, link_offsets, screen_lookup, len(screens)
                )
                matrix = np.zeros((destinations.size, len(screens)), dtype=np.float64)
                if path_index.size:
                    pair = path_od[path_index].astype(np.int64) * len(screens) + screen_index
                    flat = np.bincount(
                        pair,
                        weights=theta[path_index].astype(np.float64) * count,
                        minlength=destinations.size * len(screens),
                    )
                    matrix = flat.reshape(destinations.size, len(screens))
                positive = np.flatnonzero(q > 0)
                if positive.size == 0:
                    continue
                selected = matrix[positive]
                nz_row, nz_screen = np.nonzero(selected)
                counts = np.bincount(nz_row, minlength=positive.size)
                internal_destination = np.array(
                    [external_to_internal.get(int(value), -1) for value in destinations[positive]],
                    dtype=np.int32,
                )
                if np.any(internal_destination < 0):
                    raise ValueError("DTAC destination is absent from node.csv zone map")
                policy_origin.append(
                    np.full(positive.size, origin_index, dtype=np.int32)
                )
                policy_destination.append(internal_destination)
                policy_q.append(q[positive].astype(np.float64))
                policy_indices.append(nz_screen.astype(np.int32))
                policy_data.append(selected[nz_row, nz_screen].astype(np.float64))
                policy_indptr.append(
                    np.cumsum(counts, dtype=np.int64) + policy_nnz
                )
                policy_nnz += int(nz_row.size)
                if (origin_index + 1) % 500 == 0:
                    LOGGER.info(
                        "%s %s policy: origin %d/%d",
                        period.upper(),
                        mode["name"],
                        origin_index + 1,
                        len(zones),
                    )
            q0 = np.concatenate(policy_q) if policy_q else np.empty(0, np.float64)
            path = output / "od_screen_policy_{}.npz".format(mode["name"])
            np.savez(
                path,
                origin=(
                    np.concatenate(policy_origin)
                    if policy_origin
                    else np.empty(0, np.int32)
                ),
                destination=(
                    np.concatenate(policy_destination)
                    if policy_destination
                    else np.empty(0, np.int32)
                ),
                q0=q0,
                indptr=np.concatenate(policy_indptr),
                indices=(
                    np.concatenate(policy_indices)
                    if policy_indices
                    else np.empty(0, np.int32)
                ),
                data=(
                    np.concatenate(policy_data)
                    if policy_data
                    else np.empty(0, np.float64)
                ),
                zone_external=zones,
                screen_ids=screens,
            )
            mode_stats.append(
                {
                    "mode": mode["name"],
                    "supported_od_cells": int(q0.size),
                    "policy_nonzeros": int(policy_nnz),
                    "theta_sum_max_abs_error": theta_error,
                    "file": path.name,
                }
            )
        if stream.read(1):
            raise ValueError("unexpected trailing bytes in route_columns.bin")
    return {
        "period": period.upper(),
        "screen_count": int(len(screens)),
        "zone_count": int(len(zones)),
        "modes": mode_stats,
    }
