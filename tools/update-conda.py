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

# Get source directory from command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("sourcedir", help="Source directory")
options = parser.parse_args()

platforms = {
    "linux-64": "linux",
    "linux-aarch64": "linux-aarch64",
    "osx-64": "macos",
    "osx-arm64": "macos-arm64",
    #"win-64": "win",
}
pythons = ["3.9", "3.10", "3.11"]
tags = ["", "-dev"]
sources = ["src"]


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


def update_conda(source_dir: Path) -> None:
    pyproject_toml = source_dir / "pyproject.toml"
    if not pyproject_toml.exists():
        print(f"pyproject.toml not found in {pyproject_toml}")
        return

    dependencies = get_dependencies(pyproject_toml)

    for platform_key, platform_value in platforms.items():
        for python in pythons:
            for tag in tags:
                # Pin Python version
                pinned_dependencies = {
                    f"python={python}" if dep == "python" else dep
                    for dep in dependencies
                }
                if tag == "-dev":
                    dev_dependencies = get_dev_dependencies(pyproject_toml)
                    print(f"Adding dev dependencies: {dev_dependencies}")
                    pinned_dependencies = pinned_dependencies.union(dev_dependencies)
                pinned_dependencies = sorted(pinned_dependencies)

                for src in sources:
                    env_file = source_dir / src / f"environment{tag}-{python}.yml"
                    write_env_file(env_file, pinned_dependencies)
                    lock_file = (
                        source_dir / src / f"environment{tag}-{python}-{platform_value}"
                    )
                    print(
                        f"Updating lock file for {env_file} at {lock_file}", flush=True
                    )
                    subprocess.run(
                        [
                            "conda-lock",
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
                        ]
                    )


def get_dependencies(pyproject_toml: Path) -> list[str]:
    graystull_config = Configuration("sagemath")
    pyproject_metadata = merge_setup_toml_metadata(
        {}, get_all_toml_info(pyproject_toml)
    )
    requirements = extract_requirements(pyproject_metadata, graystull_config, {})
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
    python_requirements = pyproject_metadata.get("install_requires", [])
    # Specify concrete packages for some packages not yet in grayskull
    # TODO: It seems to be a bug that these external.dependencies are added to install_requires by grayskull
    python_requirements.remove("pkg:generic/tachyon")
    python_requirements.append("tachyon")
    python_requirements.remove("pkg:generic/sagemath-elliptic-curves")
    python_requirements.append("sagemath-db-elliptic-curves")
    python_requirements.remove("pkg:generic/sagemath-polytopes-db")
    python_requirements.append("sagemath-db-polytopes")
    python_requirements.remove("pkg:generic/sagemath-graphs")
    python_requirements.append("sagemath-db-graphs")
    # Following can be removed once https://github.com/regro/cf-scripts/pull/2176 is used in grayskull
    python_requirements = [
        req.replace("lrcalc", "python-lrcalc") for req in python_requirements
    ]
    all_requirements += normalize_requirements_list(
        python_requirements, graystull_config
    )
    all_requirements.remove("<{ pin_compatible('numpy') }}")

    # Add version constraints for some packages (not yet supported by grayskull/PEP 725)
    all_requirements.remove("c-compiler")
    all_requirements.append("c-compiler <=1.6")
    all_requirements.remove("cxx-compiler")
    all_requirements.append("cxx-compiler <=1.6")
    return all_requirements


def get_dev_dependencies(pyproject_toml: Path) -> list[str]:
    pyproject = tomllib.load(pyproject_toml)
    dependency_groups = pyproject["dependency-groups"]
    dev_dependencies = dependency_groups.get("test", []) + dependency_groups.get(
        "docs", []
    )
    return dev_dependencies


update_conda(Path(options.sourcedir))
