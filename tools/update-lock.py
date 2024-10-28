#!/usr/bin/env python3
# pyright: strict

import argparse
import subprocess
from pathlib import Path

import toml


def compile_requirements(pyproject_file: Path, groups: dict[str, list[str]]) -> None:
    requirements_dir = pyproject_file.parent / "requirements"
    if not requirements_dir.exists():
        requirements_dir.mkdir()
    subprocess.run(
        [
            "pip-compile",
            "-v",
            "--upgrade",
            "--no-build-isolation",
            "--output-file",
            str(requirements_dir / "requirements.txt"),
            pyproject_file,
        ],
        check=True,
    )
    for group in groups:
        requirements_file = requirements_dir / f"{group}-requirements.txt"
        requirements_in_file = requirements_dir / f"{group}.in"
        with requirements_in_file.open("w") as f:
            f.write("\n".join(groups[group]))
        subprocess.run(
            [
                "pip-compile",
                "-v",
                "--upgrade",
                "--no-build-isolation",
                "--output-file",
                str(requirements_file),
                str(requirements_in_file),
            ],
            check=True,
        )


parser = argparse.ArgumentParser(
    description="Update lock files under requirements using pip-tools compile."
)
parser.add_argument(
    "sourcedir", help="Source directory", nargs="?", default=".", type=Path
)
options = parser.parse_args()

pyproject_file = Path(options.sourcedir) / "pyproject.toml"
with pyproject_file.open("r") as f:
    pyproject = toml.load(f)

groups = pyproject.get("dependency-groups", {})
compile_requirements(pyproject_file, groups)
