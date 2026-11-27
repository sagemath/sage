#!/usr/bin/env python3
"""Synchronize shared BibTeX files into the Meson documentation build tree."""

import argparse
from pathlib import Path
from tempfile import NamedTemporaryFile


def copy_bibliographies(source: Path, destination: Path, stamp: Path) -> None:
    """Copy changed files, remove obsolete copies, and stamp actual changes.

    The destination is the build tree's dedicated bibliography directory.
    Its ``*.bib`` files are owned by this helper; other files are untouched.
    Keeping an unchanged stamp lets Ninja avoid rebuilding the documentation.
    """
    if source.resolve() == destination.resolve():
        raise ValueError("source and destination must be different directories")
    if not source.is_dir():
        raise FileNotFoundError(source)
    contents = {path.name: path.read_bytes() for path in sorted(source.glob("*.bib"))}
    destination.mkdir(parents=True, exist_ok=True)
    changed = False
    for name, data in contents.items():
        target = destination / name
        if target.is_file() and target.read_bytes() == data:
            continue
        with NamedTemporaryFile(dir=destination, delete=False) as temporary:
            pending = Path(temporary.name)
            temporary.write(data)
        try:
            pending.replace(target)
        finally:
            pending.unlink(missing_ok=True)
        changed = True
    for target in destination.glob("*.bib"):
        if target.name not in contents:
            target.unlink()
            changed = True
    if changed or not stamp.exists():
        stamp.touch()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("stamp", type=Path)
    args = parser.parse_args()
    copy_bibliographies(args.source, args.destination, args.stamp)
