#!/usr/bin/env python3
"""Build Sage using the dependencies installed in the active Pixi environment."""

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path


def cmake_executable(prefix):
    """Use CMake installed by the lockfile, including its Windows layout."""
    return prefix / ("Library/bin/cmake.exe" if sys.platform == "win32" else "bin/cmake")


def native_file_contents(prefix):
    """Keep CMake dependency discovery inside the active native environment."""
    roots = [prefix]
    if sys.platform == "win32":
        roots.append(prefix / "Library")

    def quoted(path):
        return "'" + path.as_posix().replace("'", "\\'") + "'"

    paths = "[" + ", ".join(map(quoted, roots)) + "]"
    return (
        "[binaries]\n"
        f"cmake = {quoted(cmake_executable(prefix))}\n"
        "[built-in options]\n"
        f"cmake_prefix_path = {paths}\n"
        "[cmake]\n"
        f"CMAKE_FIND_ROOT_PATH = {paths}\n"
        "CMAKE_FIND_ROOT_PATH_MODE_PACKAGE = 'ONLY'\n"
        "CMAKE_FIND_ROOT_PATH_MODE_INCLUDE = 'ONLY'\n"
        "CMAKE_FIND_ROOT_PATH_MODE_LIBRARY = 'ONLY'\n"
        # These control dependency lookup, not the compiler's system SDK.
        "CMAKE_FIND_ROOT_PATH_MODE_PROGRAM = 'NEVER'\n"
        "CMAKE_FIND_USE_CMAKE_ENVIRONMENT_PATH = 'FALSE'\n"
        "CMAKE_FIND_USE_PACKAGE_REGISTRY = 'FALSE'\n"
        "CMAKE_FIND_USE_SYSTEM_PACKAGE_REGISTRY = 'FALSE'\n"
    )


def environment_fingerprint(prefix):
    """Identify the installed native dependencies, independent of other environments."""
    packages = []
    for record in (prefix / "conda-meta").glob("*.json"):
        metadata = json.loads(record.read_text())
        packages.append(tuple(metadata[key] for key in ("name", "version", "build")))
    if not packages:
        raise RuntimeError(f"No Conda package records found in {prefix}")
    data = json.dumps([sys.version, sorted(packages), native_file_contents(prefix)]).encode()
    return hashlib.sha256(data).hexdigest()


def prepare_build_directory(prefix, clean=False):
    """Discard cached compiler/dependency results when this environment changes."""
    fingerprint = environment_fingerprint(prefix)
    build_dir = prefix / "sage-build"
    state = prefix / "sage-build.environment"
    previous = state.read_text().strip() if state.is_file() else None
    if build_dir.is_symlink():
        raise RuntimeError(f"Refusing to use a symbolic link as build directory: {build_dir}")
    if build_dir.exists() and (clean or previous != fingerprint):
        shutil.rmtree(build_dir)
    return state, fingerprint


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wheel", action="store_true",
                        help="install a regular wheel instead of an editable installation")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="show detailed pip and compiler output")
    parser.add_argument("--jobs", type=int, help="maximum number of parallel build jobs")
    parser.add_argument("--clean", action="store_true",
                        help="discard this environment's Sage build directory before building")
    parser.add_argument("pip_args", nargs=argparse.REMAINDER,
                        help="additional pip arguments, following --")
    args = parser.parse_args()
    if args.jobs is not None and args.jobs < 1:
        parser.error("--jobs must be positive")
    prefix = Path(sys.prefix).resolve()
    if not os.environ.get("PIXI_ENVIRONMENT_NAME") or prefix != Path(
            os.environ.get("CONDA_PREFIX", "")).resolve():
        parser.error("run this command with 'pixi run build'")
    if not cmake_executable(prefix).is_file():
        parser.error("CMake is missing from the Pixi environment; run 'pixi install'")
    native_file = prefix / "sage-build-native.ini"

    # Editable Meson installations rebuild on import. Their build directory must
    # survive pip and must not be shared with another Python/native environment.
    command = [sys.executable, "-m", "pip", "install", "--no-build-isolation",
               "--no-deps", f"--config-settings=build-dir={prefix / 'sage-build'}",
               "--config-settings=setup-args=--wrap-mode=nofallback",
               f"--config-settings=setup-args=--native-file={native_file}"]
    if not args.wheel:
        command.append("--editable")
    command.append(".")
    if args.verbose:
        command.append("-v")
    if args.jobs is not None:
        command.append(f"--config-settings=compile-args=-j{args.jobs}")
    command.extend(args.pip_args[1:] if args.pip_args[:1] == ["--"] else args.pip_args)

    environment = os.environ.copy()
    if sys.platform == "win32":
        # MSVC needs the native libraries and headers provided by conda-forge.
        for variable, subdirectory in (("LIB", "lib"), ("INCLUDE", "include")):
            directory = str(prefix / "Library" / subdirectory)
            entries = environment.get(variable, "").split(os.pathsep)
            if directory not in entries:
                environment[variable] = os.pathsep.join(filter(None, [*entries, directory]))
    state, fingerprint = prepare_build_directory(prefix, clean=args.clean)
    contents = native_file_contents(prefix)
    if not native_file.is_file() or native_file.read_text() != contents:
        native_file.write_text(contents)
    # This records the dependency identity, not build success. Meson/ninja can
    # resume a failed build as long as its native dependencies have not changed.
    state.write_text(fingerprint + "\n")
    subprocess.run(command, cwd=Path(__file__).resolve().parents[1],
                   env=environment, check=True)


if __name__ == "__main__":
    main()
