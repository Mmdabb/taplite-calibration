"""Portable path helpers.

Configuration files store only relative paths. Runtime code resolves them once
against the project directory and writes relative paths back to manifests.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

from .errors import ConfigurationError


PathLike = Union[str, Path]


def resolve_project_path(project_dir: Path, value: PathLike, field: str) -> Path:
    raw = Path(value)
    if raw.is_absolute():
        raise ConfigurationError(
            "{} must be relative to the project directory, got {!r}".format(
                field, str(value)
            )
        )
    return (project_dir / raw).resolve()


def optional_project_path(
    project_dir: Path, value: Optional[PathLike], field: str
) -> Optional[Path]:
    if value is None or str(value).strip() == "":
        return None
    return resolve_project_path(project_dir, value, field)


def portable_path(path: Path, project_dir: Path) -> str:
    """Return a POSIX-style path relative to the project where possible."""
    try:
        return path.resolve().relative_to(project_dir.resolve()).as_posix()
    except ValueError:
        # A user may intentionally wire an external dataset through ``..``.
        # relpath still keeps the persisted manifest portable on the same layout.
        import os

        return Path(os.path.relpath(path.resolve(), project_dir.resolve())).as_posix()


def ensure_file(path: Path, field: str) -> None:
    if not path.is_file():
        raise ConfigurationError("{} does not exist: {}".format(field, path))


def ensure_directory(path: Path, field: str) -> None:
    if not path.is_dir():
        raise ConfigurationError("{} does not exist: {}".format(field, path))
