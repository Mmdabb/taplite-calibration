"""Python interface to the bundled TAPLite C++ kernel.

The solver is compiled into :mod:`taplite_calibration._native`; no TAPLite
executable, adjacent checkout, or external ``pytaplite`` package is required.
Assignment calls use a short-lived Python worker by default because TAPLite's
large global route store is intentionally process-scoped.
"""

from __future__ import annotations

import csv
import os
import shutil
import subprocess
import sys
import tempfile
import threading
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Union

from .errors import CapabilityError


PathLike = Union[str, os.PathLike]
_DIRECT_CALL_LOCK = threading.Lock()


def _load_extension():
    try:
        from . import _native
    except ImportError as error:  # pragma: no cover - exercised by sdist misuse
        raise CapabilityError(
            "the bundled TAPLite native extension is unavailable; install a "
            "platform wheel or build this package with a C++17 compiler"
        ) from error
    return _native


def _read_links(run_dir: Path) -> List[Dict[str, str]]:
    path = run_dir / "link_performance.csv"
    if not path.is_file():
        return []
    with path.open("r", newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


@dataclass(frozen=True)
class NativeResult:
    """Outcome and final link table from one native kernel invocation."""

    run_dir: Path
    returncode: int
    log: str
    links: List[Dict[str, str]]
    isolated: bool

    @staticmethod
    def _number(row: Mapping[str, str], *keys: str) -> float:
        for key in keys:
            value = row.get(key)
            if value not in (None, ""):
                try:
                    return float(value)
                except (TypeError, ValueError):
                    pass
        return 0.0

    def summary(self) -> Dict[str, object]:
        loaded = [
            row
            for row in self.links
            if self._number(row, "volume", "vehicle_volume") > 0.0
        ]
        return {
            "links": len(self.links),
            "loaded_links": len(loaded),
            "total_VMT": round(sum(self._number(row, "VMT") for row in self.links), 1),
            "total_VHT": round(sum(self._number(row, "VHT") for row in self.links), 1),
            "returncode": self.returncode,
            "isolated_worker": self.isolated,
        }


def _copy_scenario(source: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    for source_file in source.iterdir():
        if source_file.is_file():
            shutil.copy2(source_file, destination / source_file.name)


def _prepare_run_directory(
    scenario: PathLike, in_place: bool, work_dir: Optional[PathLike], prefix: str
) -> Path:
    source = Path(scenario).expanduser().resolve()
    if not source.is_dir():
        raise FileNotFoundError("scenario folder not found: {}".format(source))
    if in_place:
        return source
    destination = (
        Path(work_dir).expanduser().resolve()
        if work_dir is not None
        else Path(tempfile.mkdtemp(prefix=prefix)).resolve()
    )
    if destination != source:
        _copy_scenario(source, destination)
    return destination


def _apply_settings_overrides(run_dir: Path, overrides: Mapping[str, object]) -> None:
    path = run_dir / "settings.csv"
    if not path.is_file():
        raise FileNotFoundError("TAPLite settings file not found: {}".format(path))
    with path.open("r", newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 1:
        raise ValueError("{} must contain exactly one settings row".format(path))
    row: Dict[str, object] = dict(rows[0])
    row.update(overrides)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(row), lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)


def _run_direct(kind: str, run_dir: Path) -> tuple[int, str]:
    extension = _load_extension()
    function = (
        extension.run_in_dir
        if kind == "assignment"
        else extension.run_auto_calibration_in_dir
    )
    previous = Path.cwd()
    with _DIRECT_CALL_LOCK:
        try:
            os.chdir(run_dir)
            return int(function(str(run_dir))), "in-process bundled native extension"
        finally:
            os.chdir(previous)


def _run_isolated(kind: str, run_dir: Path, timeout: int) -> tuple[int, str]:
    process = subprocess.run(
        [
            sys.executable,
            "-m",
            "taplite_calibration._native_worker",
            kind,
            str(run_dir),
        ],
        cwd=str(run_dir),
        timeout=timeout,
        capture_output=True,
        text=True,
        check=False,
    )
    log = (process.stdout or "") + (process.stderr or "")
    return int(process.returncode), log


def _run(
    kind: str,
    scenario: PathLike,
    *,
    in_place: bool,
    work_dir: Optional[PathLike],
    timeout: int,
    isolated: bool,
    settings_overrides: Optional[Mapping[str, object]],
) -> NativeResult:
    run_dir = _prepare_run_directory(
        scenario,
        in_place,
        work_dir,
        "taplite_calibration_" if kind == "auto-calibration" else "taplite_assignment_",
    )
    overrides = dict(settings_overrides or {})
    if kind == "auto-calibration":
        overrides["auto_calibration"] = 1
        overrides["column_output"] = 2
    if overrides:
        _apply_settings_overrides(run_dir, overrides)
    if timeout < 1:
        raise ValueError("timeout must be at least one second")
    returncode, log = (
        _run_isolated(kind, run_dir, timeout)
        if isolated
        else _run_direct(kind, run_dir)
    )
    return NativeResult(
        run_dir=run_dir,
        returncode=returncode,
        log=log,
        links=_read_links(run_dir),
        isolated=isolated,
    )


def assign(
    scenario: PathLike,
    *,
    in_place: bool = True,
    work_dir: Optional[PathLike] = None,
    timeout: int = 3600,
    isolated: bool = True,
    settings_overrides: Optional[Mapping[str, object]] = None,
) -> NativeResult:
    """Run ordinary TAPLite user-equilibrium assignment with the bundled kernel."""

    return _run(
        "assignment",
        scenario,
        in_place=in_place,
        work_dir=work_dir,
        timeout=timeout,
        isolated=isolated,
        settings_overrides=settings_overrides,
    )


def auto_calibrate(
    scenario: PathLike,
    *,
    in_place: bool = True,
    work_dir: Optional[PathLike] = None,
    timeout: int = 86400,
    isolated: bool = True,
    settings_overrides: Optional[Mapping[str, object]] = None,
) -> NativeResult:
    """Run the dedicated equilibrium-coupled refined QVDF calibration API.

    Network and demand are initialized once inside the worker. Every inner
    equilibrium keeps route columns and link-arrival state in memory with
    ``column_output=2``; only accepted final outputs are written.
    """

    return _run(
        "auto-calibration",
        scenario,
        in_place=in_place,
        work_dir=work_dir,
        timeout=timeout,
        isolated=isolated,
        settings_overrides=settings_overrides,
    )


def native_status(requested_threads: int = 0) -> Dict[str, object]:
    """Return build/runtime diagnostics for the bundled OpenMP kernel."""

    if requested_threads < 0:
        raise ValueError("requested_threads must be nonnegative")
    return dict(_load_extension().openmp_status(requested_threads))

