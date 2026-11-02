"""Checks shared by the Homebrew dependency bottle workflow phases."""

import importlib
from importlib import metadata
from importlib.resources import files
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tomllib


def requirements():
    """Use the candidate's locked versions for the two pytest entry points."""
    with open("uv.lock", "rb") as stream:
        packages = tomllib.load(stream)["package"]
    for name in ("pytest", "pytest-xdist"):
        versions = {p["version"] for p in packages if p["name"] == name}
        if len(versions) != 1:
            raise RuntimeError(f"Expected one locked version of {name}: {versions}")
        print(f"{name}=={versions.pop()}")


def data():
    """Expose bottled databases also to checkouts predating Sage PR #42878."""
    root = Path(os.environ["SAGE_DATA_PATH"])
    root.mkdir()
    directories = {
        "cremona": files("sage_data_elliptic_curves.data") / "cremona",
        "ellcurves": files("sage_data_elliptic_curves.data") / "ellcurves",
        "graphs": files("sage_data_graphs.data"),
        "reflexive_polytopes": files("sage_data_polytopes") / "data",
    }
    for name, path in directories.items():
        if not path.is_dir():
            raise RuntimeError(f"Missing bottled data: {path}")
        (root / name).symlink_to(path, target_is_directory=True)
    print(json.dumps({name: str(path) for name, path in directories.items()}, indent=2))


def provenance():
    """Reject an installed Sage bottle accidentally taking precedence."""
    source = Path(os.environ["CI_SOURCE"]).resolve()
    build = Path(os.environ["CI_RUN_DIR"], "build").resolve()
    result = {
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "python": sys.executable,
        "modules": {},
        "dependencies": {},
    }
    # Check both pure Python and compiled modules.
    for name, expected in (("sage.all", source), ("sage.rings.integer", build)):
        path = Path(importlib.import_module(name).__file__).resolve()
        if not path.is_relative_to(expected):
            raise RuntimeError(f"{name} came from {path}, expected it under {expected}")
        result["modules"][name] = str(path)
    dist = metadata.distribution("sagemath")
    direct_url = json.loads(dist.read_text("direct_url.json"))
    if direct_url != {"dir_info": {"editable": True}, "url": source.as_uri()}:
        raise RuntimeError(f"Not the candidate's editable installation: {direct_url}")
    for name, formula in (("sympy", "sagemath-deps"), ("sphinx", "sagemath-docbuild-deps")):
        prefix = Path(subprocess.check_output(
            ["brew", "--prefix", f"dimpase/tap/{formula}"], text=True).strip()).resolve()
        dist = metadata.distribution(name)
        path = Path(dist.locate_file("")).resolve()
        if not path.is_relative_to(prefix):
            raise RuntimeError(f"{name} came from {path}, expected the {formula} bottle")
        result["dependencies"][name] = {"version": dist.version, "path": str(path)}
    print(json.dumps(result, indent=2))


def doctest():
    paths = shlex.split(os.environ.get("CI_DOCTEST_PATHS", ""))
    if not paths:
        if os.environ.get("CI_SUITE", "smoke") == "all":
            paths = ["src/sage"]
        else:
            paths = [
                "src/sage/arith/misc.py",
                "src/sage/rings/integer.pyx",
                "src/sage/matrix/matrix_double_dense.pyx",
                "src/sage/libs/pari",
                "src/sage/libs/eclib",
                "src/sage/libs/gap",
                "src/sage/interfaces/maxima_lib.py",
            ]
    root = Path("src/sage").resolve()
    for value in paths:
        path = Path(value).resolve()
        if not path.is_relative_to(root) or not path.exists():
            raise ValueError(f"Expected an existing path under src/sage: {value}")
    arguments = ["-p", os.environ.get("CI_JOBS", "4"),
                 "--format", "github", "--random-seed=0", "--timeout=900",
                 f"--stats-path={os.environ['CI_RESULTS']}/doctest-stats.json", *paths]
    from sage.doctest.__main__ import main
    print("Running Sage doctests: " + shlex.join(arguments), flush=True)
    sys.argv = ["sage.doctest", *arguments]
    raise SystemExit(main())


if __name__ == "__main__":
    {"requirements": requirements, "data": data,
     "provenance": provenance, "doctest": doctest}[sys.argv[1]]()
