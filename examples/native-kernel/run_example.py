"""Run the bundled native TAPLite auto-calibration API on a tiny network."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

from taplite_calibration import auto_calibrate, native_status


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("native-example-output"))
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError("output already exists: {}".format(output))

    source = Path(__file__).resolve().parent
    output.mkdir(parents=True)
    for path in source.glob("*.csv"):
        shutil.copy2(path, output / path.name)

    print("Native status:", native_status(2))
    result = auto_calibrate(
        output,
        settings_overrides={"number_of_processors": 2},
    )
    print("Result:", result.summary())
    print("Outputs:", result.run_dir)
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())

