"""Private one-run worker for the process-scoped TAPLite native kernel."""

from __future__ import annotations

import sys
from pathlib import Path


def main() -> int:
    if len(sys.argv) != 3 or sys.argv[1] not in ("assignment", "auto-calibration"):
        sys.stderr.write(
            "usage: python -m taplite_calibration._native_worker "
            "{assignment|auto-calibration} SCENARIO\n"
        )
        return 64
    from . import _native

    run_dir = str(Path(sys.argv[2]).resolve())
    function = (
        _native.run_in_dir
        if sys.argv[1] == "assignment"
        else _native.run_auto_calibration_in_dir
    )
    return int(function(run_dir))


if __name__ == "__main__":
    raise SystemExit(main())

