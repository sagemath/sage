#!/usr/bin/env python3
# See README.md for more details

import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import toml as tomllib
from grayskull.config import Configuration
from grayskull.strategy.py_base import merge_setup_toml_metadata
from grayskull.strategy.py_toml import get_all_toml_info
from grayskull.strategy.pypi import extract_requirements, normalize_requirements_list
from packaging.requirements import Requirement

platforms = {
    "linux-64": "linux",
    "linux-aarch64": "linux-aarch64",
    "osx-64": "macos-x86_64",
    "osx-arm64": "macos",
    "win-64": "win",
}

# Get source directory from command line arguments
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
    choices=platforms.keys(),
)
options = parser.parse_args()
pythons = ["3.11", "3.12", "3.13"]
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


def filter_requirements(dependencies: set[str], python: str, platform: str) -> set[str]:
    sys_platform = {
        "linux-64": "linux",
        "linux-aarch64": "linux",
        "osx-64": "darwin",
        "osx-arm64": "darwin",
        "win-64": "win32",
    }[platform]
    platform_machine = {
        "linux-64": "x86_64",
        "linux-aarch64": "aarch64",
        "osx-64": "x86_64",
        "osx-arm64": "arm64",
        "win-64": "x86_64",
    }[platform]
    env = {
        "python_version": python,
        "sys_platform": sys_platform,
        "platform_machine": platform_machine,
    }

    def filter_dep(dep: str):
        req = Requirement(dep)
        if not req.marker or req.marker.evaluate(env):
            # Serialize the requirement without the marker
            req.marker = None
            return str(req)
        return None

    return set(filter(None, map(filter_dep, dependencies)))


def update_conda(source_dir: Path, systems: list[str] | None) -> None:
    pyproject_toml = source_dir / "pyproject.toml"
    if not pyproject_toml.exists():
        print(f"pyproject.toml not found in {pyproject_toml}")
        return

    def process_platform_python(platform_key, platform_value, python):
        dependencies = get_dependencies(pyproject_toml, python, platform_key)
        for tag in tags:
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

    with ThreadPoolExecutor(max_workers=3) as executor:
        futures = [
            executor.submit(process_platform_python, platform_key, platform_value, python)
            for platform_key, platform_value in platforms.items()
            for python in pythons
            if not (systems and platform_key not in systems)
        ]
        for future in futures:
            future.result()

def get_dependencies(pyproject_toml: Path, python: str, platform: str) -> set[str]:
    grayskull_config = Configuration("sagemath")
    pyproject = tomllib.load(pyproject_toml)
    pyproject_metadata = merge_setup_toml_metadata(
        {}, get_all_toml_info(pyproject_toml)
    )
    requirements = extract_requirements(pyproject_metadata, grayskull_config, {})
    all_requirements: set[str] = set(
        requirements.get("build", [])
        + requirements.get("host", [])
        + requirements.get("run", [])
        + pyproject_metadata.get("install_requires", [])
        + get_dev_dependencies(pyproject)
        + get_optional_dependencies(pyproject)
    )

    # Fix requirements that are not available on conda
    all_requirements = {
        # Following can be removed once https://github.com/regro/cf-scripts/pull/2176 is used in grayskull
        req.replace("lrcalc", "python-lrcalc")
        .replace("symengine", "python-symengine")
        .replace("memory_allocator", "memory-allocator")
        .replace("pkg:generic/r-lattice", "r-lattice")
        .replace("pkg:generic/latexmk", "latexmk")
        .replace("pkg:generic/sagemath-elliptic-curves", "sagemath-db-elliptic-curves")
        .replace("pkg:generic/sagemath-graphs", "sagemath-db-graphs")
        .replace("pkg:generic/sagemath-polytopes-db", "sagemath-db-polytopes")
        .replace("pkg:generic/tachyon", "tachyon")
        .replace("brial", "libbrial") # on Conda, 'brial' refers to the Python package
        for req in all_requirements
    }
    # Exclude requirements not available on conda (for a given platform)
    exclude_packages: set[str] = {
        "p_group_cohomology",
        "sage_numerical_backends_coin",
        "sagemath_giac",
        "pynormaliz",  # due to https://github.com/sagemath/sage/issues/40214
        "latte-integrale",  # due to https://github.com/sagemath/sage/issues/40216
    }
    if platform in ("linux-aarch64", "osx-arm64"):
        exclude_packages |= {
            "4ti2",
            "latte-integrale",
            "lrslib",
        }
    elif platform == "win-64":
        exclude_packages |= {
            "4ti2",
            "bc",
            "libbrial",
            "bliss",
            "cddlib",
            "cliquer",
            "ecl",
            "eclib",
            "ecm",
            "fflas-ffpack",
            "fplll",
            "gap-defaults",
            "gengetopt",
            "gfan",
            "giac",
            "givaro",
            "iml",
            "latte-integrale",
            "lcalc",
            "libatomic_ops",
            "libbraiding",
            "libhomfly",
            "linbox",
            "lrcalc",
            "lrslib",
            "m4",
            "m4rie",
            "maxima",
            "mpfi",
            "ncurses",
            "ntl",
            "palp",
            "patch",
            "ppl",
            "primecount",
            "pynormaliz",
            "python-lrcalc",
            "readline",
            "rpy2",
            "rw",
            "singular",
            "sirocco",
            "sympow",
            "tachyon",
            "tar",
            "texinfo",
        }
    all_requirements = {
        req
        for req in all_requirements
        if not any(
            req == package or req.startswith(package + " ")
            for package in exclude_packages
        )
    }

    # Remove virtual packages to not confuse 'filter_requirements'
    all_requirements.remove("{{ blas }}")
    all_requirements.remove("{{ compiler('c') }}")
    all_requirements.remove("{{ compiler('cxx') }}")
    all_requirements.discard("<{ pin_compatible('numpy') }}")
    # For some reason, grayskull mishandles the fortran compiler sometimes
    # so handle both cases
    for item in ["{{ compiler('fortran') }}", "{{ compiler'fortran' }}"]:
        try:
            all_requirements.remove(item)
        except (ValueError, KeyError):
            pass
    for with_comment in {req for req in all_requirements if "#" in req}:
        all_requirements.discard(with_comment)

    all_requirements = filter_requirements(all_requirements, python, platform)
    all_requirements = set(
        normalize_requirements_list(list(all_requirements), grayskull_config)
    )
    # Specify concrete package for some virtual packages
    all_requirements.add("blas=2.*=openblas")
    all_requirements.add("fortran-compiler")
    if platform == "win-64":
        all_requirements.add("vs2022_win-64")
        # For mingw:
        # all_requirements.add("gcc_win-64 >= 14.2.0")
        # all_requirements.add("gxx_win-64")
    else:
        all_requirements.add("c-compiler")
        all_requirements.add("cxx-compiler")

    # Add additional dependencies based on platform
    if platform == "win-64":
        # Flint needs pthread.h
        all_requirements.add("winpthreads-devel")
        # Workaround for https://github.com/conda-forge/libpng-feedstock/issues/47
        all_requirements.add("zlib")
    if platform != "win-64":
        # Needed to run configure/bootstrap, can be deleted once we fully migrated to meson
        all_requirements.add("autoconf")
        all_requirements.add("automake")
        all_requirements.add("m4")
    # Needed to fix a bug on Macos with broken pkg-config
    all_requirements.add("expat")

    # Packages with version constraints
    # https://github.com/sagemath/sage/pull/40679
    if platform != "win-64":
        all_requirements.remove("maxima")
        all_requirements.add("maxima < 5.48.0")

    return all_requirements


def get_dev_dependencies(pyproject: dict) -> list[str]:
    dependency_groups = pyproject.get("dependency-groups", {})
    dev_dependencies = (
        dependency_groups.get("test", [])
        + dependency_groups.get("docs", [])
        + dependency_groups.get("lint", [])
        + dependency_groups.get("dev", [])
    )
    # Remove dependencies that are not available on conda
    dev_dependencies.remove("relint")
    return dev_dependencies


def get_optional_dependencies(pyproject: dict) -> list[str]:
    optional_dependencies = []
    optional_groups = pyproject.get("project", {}).get("optional-dependencies", {})
    for _, dependencies in optional_groups.items():
        optional_dependencies.extend(dependencies)
    # print(f"Optional dependencies: {optional_dependencies}")  # Uncommented for debugging
    return optional_dependencies


update_conda(options.sourcedir, options.systems)
