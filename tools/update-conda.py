#!/usr/bin/env python3
# Requirements: pip install https://github.com/conda/grayskull conda-lock
# Usage: python tools/update-conda.py .

import argparse
import subprocess
from pathlib import Path

import toml as tomllib
from grayskull.config import Configuration
from grayskull.strategy.py_base import merge_setup_toml_metadata
from grayskull.strategy.py_toml import get_all_toml_info
from grayskull.strategy.pypi import extract_requirements, normalize_requirements_list
from packaging.requirements import Requirement

# Get source directory from command line arguments
parser = argparse.ArgumentParser()
parser.add_argument(
    "sourcedir", help="Source directory", nargs="?", default=".", type=Path
)
options = parser.parse_args()

platforms = {
    "linux-64": "linux",
    "linux-aarch64": "linux-aarch64",
    "osx-64": "macos-x86_64",
    "osx-arm64": "macos",
    # "win-64": "win",
}
pythons = ["3.9", "3.10", "3.11"]
tags = [""]


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


def filter_requirements(dependencies: set[str], python: str) -> set[str]:
    def filter_dep(dep: str):
        req = Requirement(dep)
        env = {"python_version": python}
        if not req.marker or req.marker.evaluate(env):
            # Serialize the requirement without the marker
            req.marker = None
            return str(req)
        return None

    return set(filter(None, map(filter_dep, dependencies)))


def update_conda(source_dir: Path) -> None:
    pyproject_toml = source_dir / "pyproject.toml"
    if not pyproject_toml.exists():
        print(f"pyproject.toml not found in {pyproject_toml}")
        return

    for platform_key, platform_value in platforms.items():
        for python in pythons:
            dependencies = get_dependencies(pyproject_toml, python)
            for tag in tags:
                # Pin Python version
                pinned_dependencies = {
                    f"python={python}" if dep == "python" else dep
                    for dep in dependencies
                }

                dev_dependencies = get_dev_dependencies(pyproject_toml)
                print(f"Adding dev dependencies: {dev_dependencies}")
                pinned_dependencies = pinned_dependencies.union(dev_dependencies)

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
                    f.write(f"name: sage{tag or '-dev'}\n{content}")


def get_dependencies(pyproject_toml: Path, python: str) -> list[str]:
    grayskull_config = Configuration("sagemath")
    pyproject_metadata = merge_setup_toml_metadata(
        {}, get_all_toml_info(pyproject_toml)
    )
    requirements = extract_requirements(pyproject_metadata, grayskull_config, {})
    all_requirements = (
        requirements.get("build", [])
        + requirements.get("host", [])
        + requirements.get("run", [])
    )

    # Specify concrete package for some virtual packages
    all_requirements.remove("{{ blas }}")
    all_requirements.append("blas=2.*=openblas")
    all_requirements.remove("{{ compiler('c') }}")
    all_requirements.append("c-compiler")
    all_requirements.remove("{{ compiler('cxx') }}")
    all_requirements.append("cxx-compiler")

    # Correct pypi name for some packages
    python_requirements = set(pyproject_metadata.get("install_requires", []))
    # Specify concrete packages for some packages not yet in grayskull
    python_requirements.remove("pkg:generic/tachyon")
    python_requirements.add("tachyon")
    python_requirements.remove("pkg:generic/sagemath-elliptic-curves")
    python_requirements.add("sagemath-db-elliptic-curves")
    python_requirements.remove("pkg:generic/sagemath-polytopes-db")
    python_requirements.add("sagemath-db-polytopes")
    python_requirements.discard("pkg:generic/sagemath-graphs")
    python_requirements.add("sagemath-db-graphs")
    python_requirements.remove("memory_allocator")
    python_requirements.add("memory-allocator")
    # Following can be removed once https://github.com/regro/cf-scripts/pull/2176 is used in grayskull
    python_requirements = {
        req.replace("lrcalc", "python-lrcalc") for req in python_requirements
    }
    python_requirements = filter_requirements(python_requirements, python)
    all_requirements += normalize_requirements_list(
        python_requirements, grayskull_config
    )
    all_requirements.remove("<{ pin_compatible('numpy') }}")
    all_requirements.remove("memory_allocator")
    # Needed to run configure/bootstrap, can be deleted once we fully migrated to meson
    all_requirements.append("autoconf")
    all_requirements.append("automake")
    all_requirements.append("m4")
    return all_requirements


def get_dev_dependencies(pyproject_toml: Path) -> list[str]:
    pyproject = tomllib.load(pyproject_toml)
    dependency_groups = pyproject.get("dependency-groups", {})
    dev_dependencies = dependency_groups.get("test", []) + dependency_groups.get(
        "docs", []
    )
    return dev_dependencies


update_conda(options.sourcedir)
