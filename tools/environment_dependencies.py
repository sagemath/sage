"""Dependency mappings shared by the Conda and Pixi development environments."""

from __future__ import annotations

import tomllib
from typing import TYPE_CHECKING

from grayskull.config import Configuration
from grayskull.strategy.py_base import merge_setup_toml_metadata
from grayskull.strategy.py_toml import get_all_toml_info
from grayskull.strategy.pypi import extract_requirements, normalize_requirements_list
from packaging.requirements import Requirement

if TYPE_CHECKING:
    from pathlib import Path

PLATFORMS = {
    "linux-64": "linux",
    "linux-aarch64": "linux-aarch64",
    "osx-64": "macos-x86_64",
    "osx-arm64": "macos",
    "win-64": "win",
}

PYTHONS = ["3.13", "3.14"]


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


def get_dependencies(pyproject_toml: Path, python: str, platform: str) -> set[str]:
    grayskull_config = Configuration("sagemath")
    with open(pyproject_toml, "rb") as f:
        pyproject = tomllib.load(f)
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
        .replace("memory_allocator", "memory-allocator")
        .replace("pkg:generic/r-lattice", "r-lattice")
        .replace("pkg:generic/latexmk", "latexmk")
        .replace("pkg:generic/sagemath-elliptic-curves", "sagemath-db-elliptic-curves")
        .replace("pkg:generic/sagemath-graphs", "sagemath-db-graphs")
        .replace("pkg:generic/sagemath-polytopes-db", "sagemath-db-polytopes")
        .replace("pkg:generic/tachyon", "tachyon")
        .replace("pkg:generic/highs", "highs")
        .replace("pkg:generic/libatomic_ops", "libatomic_ops")
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
        "cibuildwheel",  # fails pip check since it claims to require the PyPI package patchelf which is not available on conda-forge yet (however, it just needs a patchelf which is installed from conda-forge already.)
    }
    if platform in ("linux-aarch64", "osx-arm64"):
        exclude_packages |= {
            "4ti2",
            "latte-integrale",
            "lrslib",
        }
    elif platform == "win-64":
        exclude_packages |= {
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
    if platform in ("osx-64", "osx-arm64"):
        all_requirements.add("libblas=*=*_newaccelerate")
        # Older ld64 can strip libSystem from executables on recent macOS.
        # https://github.com/conda-forge/cctools-and-ld64-feedstock/pull/90
        all_requirements.add("ld64 >=956.6")
        # GFortran 12 omits the stdin argument when piping to the current
        # cctools assembler, causing inline Fortran/f2py compilation to fail.
        all_requirements.add("gfortran >=14.4")
    else:
        all_requirements.add("openblas")
        all_requirements.add("libblas=*=*_openblas")
    # Issue #41555: meson uses pkg-config to search for cblas whose .pc files are in blas-devel
    all_requirements.add("blas-devel")
    all_requirements.add("fortran-compiler")
    # Meson uses CMake for dependency detection; use the environment's copy.
    all_requirements.add("cmake")

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
        all_requirements.add("maxima ==5.49.0")
        all_requirements.remove("singular")
        all_requirements.add("singular ==4.4.1.p5")

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
    # Remove cibuildwheel because its dependency bashlex is not available on linux-aarch64
    dev_dependencies.remove("cibuildwheel")
    return dev_dependencies


def get_optional_dependencies(pyproject: dict) -> list[str]:
    optional_dependencies = []
    optional_groups = pyproject.get("project", {}).get("optional-dependencies", {})
    for dependencies in optional_groups.values():
        optional_dependencies.extend(dependencies)
    # print(f"Optional dependencies: {optional_dependencies}")  # Uncommented for debugging
    return optional_dependencies
