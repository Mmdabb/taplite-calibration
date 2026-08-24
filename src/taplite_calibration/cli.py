"""Command-line interface; also available as ``python -m taplite_calibration``."""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import List, Optional

from . import __version__
from .config import MAX_PROCESSORS, MODES, load_config
from .errors import CalibrationError
from .pipeline import run_loaded_config
from .prepare import run_qa
from .scaffold import initialize_project


def _common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--project",
        type=Path,
        default=Path("."),
        help="Project directory containing calibration.toml (default: current directory).",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("calibration.toml"),
        help="Configuration path relative to --project.",
    )
    parser.add_argument(
        "--mode",
        choices=MODES,
        help="Override [pipeline].mode for this command.",
    )
    parser.add_argument(
        "--processors",
        type=int,
        choices=range(1, MAX_PROCESSORS + 1),
        metavar="1..20",
        help="Override the processor count; the package never permits more than 20.",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="taplite-calibration",
        description=(
            "Prepare inputs and run screen ODME, TAPLite QVDF auto calibration, "
            "or both in order."
        ),
    )
    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument(
        "--log-level",
        choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        default="INFO",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    init_parser = subparsers.add_parser("init", help="Create an empty portable project layout.")
    init_parser.add_argument("project", type=Path)

    validate_parser = subparsers.add_parser("validate", help="Validate config and selected inputs.")
    _common(validate_parser)

    run_parser = subparsers.add_parser("run", help="Execute the selected workflow.")
    _common(run_parser)
    run_parser.add_argument(
        "--run-id",
        help="Portable output directory name; defaults to a UTC timestamp plus mode.",
    )
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    try:
        if args.command == "init":
            created = initialize_project(args.project)
            print(
                json.dumps(
                    {
                        "status": "complete",
                        "project": str(args.project.resolve()),
                        "created_entries": len(created),
                        "next": "edit calibration.toml, add external inputs, then run validate",
                    },
                    indent=2,
                )
            )
            return 0
        config = load_config(
            args.project,
            config_path=args.config,
            mode_override=args.mode,
            processors_override=args.processors,
        )
        if args.command == "validate":
            reports = run_qa(config)
            passed = all(report.passed for report in reports.values())
            print(
                json.dumps(
                    {
                        "status": "PASS" if passed else "FAIL",
                        "mode": config.mode,
                        "processors": config.processors,
                        "stage_order": (
                            ["prepare", "odme", "auto-calibration"]
                            if config.mode == "both"
                            else ["prepare", config.mode]
                        ),
                        "qa": {
                            target: report.to_dict()
                            for target, report in reports.items()
                        },
                        "next": (
                            "inputs are prepared"
                            if passed
                            else "run will attempt configured preprocessing; otherwise it will exit"
                        ),
                    },
                    indent=2,
                )
            )
            return 0 if passed else 1
        result = run_loaded_config(config, run_id=args.run_id)
        print(
            json.dumps(
                {
                    "status": "complete",
                    "run_directory": str(result.run_dir),
                    "manifest": str(result.manifest_path),
                    "artifacts": result.manifest["artifacts"],
                },
                indent=2,
            )
        )
        return 0
    except (
        CalibrationError,
        FileNotFoundError,
        FileExistsError,
        RuntimeError,
        ValueError,
    ) as error:
        logging.getLogger("taplite_calibration").error("%s", error)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
