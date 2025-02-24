#!/usr/bin/env python3

from pathlib import Path
import subprocess

script_dir = Path(__file__).resolve().parent
root_dir = script_dir / ".." / ".."

subprocess.run([str(root_dir / "bootstrap-conda")])

platforms = {
    "linux-64": "linux",
    "linux-aarch64": "linux-aarch64",
    "osx-64": "macos-x86_64",
    "osx-arm64": "macos",
    # "win-64": "win",
}
pythons = ["3.9", "3.10", "3.11"]
tags = ["", "-dev"]

for platform_key, platform_value in platforms.items():
    for python in pythons:
        for tag in tags:
            env_file = root_dir / f"environment{tag}-{python}.yml"
            lock_file = root_dir / f"environment{tag}-{python}-{platform_value}"
            lock_file_gen = root_dir / f"environment{tag}-{python}-{platform_value}.yml"

            if not env_file.exists():
                continue

            print(f"Updating lock file for {env_file} at {lock_file_gen}", flush=True)
            subprocess.run(
                [
                    "conda-lock",
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
                f.write(f"name: sage{tag}\n{content}")
