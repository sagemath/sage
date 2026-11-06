#!/usr/bin/env python3
# See README.md for more details

import argparse
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from environment_dependencies import PLATFORMS, PYTHONS, get_dependencies


def write_env_file(env_file: Path, dependencies: list[str]) -> None:
    env_file.write_text(
        """name: sage
channels:
  - conda-forge
  - nodefaults
dependencies:
"""
        + "".join(f"  - {req}" + "\n" for req in dependencies)
    )
    print(f"Conda environment file written to {env_file}")


def update_conda(source_dir: Path, systems: list[str] | None) -> None:
    pyproject_toml = source_dir / "pyproject.toml"
    if not pyproject_toml.exists():
        print(f"pyproject.toml not found in {pyproject_toml}")
        return

    def process_platform_python(platform_key, platform_value, python):
        dependencies = get_dependencies(pyproject_toml, python, platform_key)
        for tag in [""]:
            # Pin Python version
            pinned_dependencies = {
                f"python={python}" if dep == "python" else dep
                for dep in dependencies
            }
            pinned_dependencies = sorted(pinned_dependencies)

            env_file = source_dir / f"environment{tag}-{python}.yml"
            write_env_file(env_file, pinned_dependencies)
            lock_file = source_dir / f"environment{tag}-{python}-{platform_value}"
            lock_file_gen = (
                source_dir / f"environment{tag}-{python}-{platform_value}.yml"
            )
            print(
                f"Updating lock file for {env_file} at {lock_file_gen}", flush=True
            )
            subprocess.run(
                [
                    # conda-lock requires virtualenv <21, unlike Sage's
                    # development environment. Run it in its own environment.
                    "uv", "tool", "run", "--python", sys.executable,
                    "--from", "conda-lock", "conda-lock",
                    "--mamba",
                    "--channel",
                    "conda-forge",
                    "--kind",
                    "env",
                    "--platform",
                    platform_key,
                    "--file",
                    str(env_file),
                    "--lockfile",
                    str(lock_file),
                    "--filename-template",
                    str(lock_file),
                ],
                check=True,
            )

            # Add conda env name to lock file at beginning
            with open(lock_file_gen, "r+") as f:
                content = f.read()
                f.seek(0, 0)
                f.write(
                    f"name: sage{tag or '-dev'}\n{content}"
                )

    with ThreadPoolExecutor(max_workers=3) as executor:
        futures = [
            executor.submit(process_platform_python, platform_key, platform_value, python)
            for platform_key, platform_value in PLATFORMS.items()
            for python in PYTHONS
            if not (systems and platform_key not in systems)
        ]
        for future in futures:
            future.result()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "sourcedir", help="Source directory", nargs="?", default=".", type=Path
    )
    parser.add_argument(
        "-s",
        "--systems",
        help="Operating systems to build for; default is all",
        nargs="+",
        type=str,
        choices=PLATFORMS.keys(),
    )
    options = parser.parse_args()
    update_conda(options.sourcedir, options.systems)


if __name__ == "__main__":
    main()
