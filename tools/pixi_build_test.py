"""Regression tests for keeping native build state tied to a Pixi environment."""

import importlib.util
import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest

spec = importlib.util.spec_from_file_location("pixi_build", Path(__file__).with_name("pixi-build.py"))
pixi_build = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pixi_build)


@pytest.fixture
def environment(tmp_path, monkeypatch):
    prefix = tmp_path / "env"
    records = prefix / "conda-meta"
    records.mkdir(parents=True)
    record = records / "native-library.json"
    record.write_text(json.dumps({"name": "native-library", "version": "1.0", "build": "0"}))
    monkeypatch.setattr(pixi_build.sys, "prefix", str(prefix))
    monkeypatch.setattr(pixi_build.sys, "argv", ["pixi-build.py"])
    monkeypatch.setenv("CONDA_PREFIX", str(prefix))
    monkeypatch.setenv("PIXI_ENVIRONMENT_NAME", "test")
    cmake = pixi_build.cmake_executable(prefix)
    cmake.parent.mkdir(parents=True)
    cmake.touch()
    return prefix, record


def test_native_package_update_discards_old_build(environment, monkeypatch):
    prefix, record = environment
    build = prefix / "sage-build"
    state = prefix / "sage-build.environment"
    old_object = build / "old-object.o"
    should_preserve = False

    def install(*args, **kwargs):
        assert old_object.exists() == should_preserve
        build.mkdir(exist_ok=True)
        old_object.write_text("compiled against the current environment")

    monkeypatch.setattr(pixi_build.subprocess, "run", install)
    pixi_build.main()
    original_state = state.read_text()

    should_preserve = True
    pixi_build.main()
    assert state.read_text() == original_state

    # A new build of the same library version can change its ABI as well.
    record.write_text(json.dumps({"name": "native-library", "version": "1.0", "build": "1"}))
    should_preserve = False
    pixi_build.main()
    assert state.read_text() != original_state
    assert record.is_file()


def test_failed_install_resumes_with_unchanged_dependencies(environment, monkeypatch):
    prefix, _ = environment
    build = prefix / "sage-build"
    state = prefix / "sage-build.environment"

    def successful_install(*args, **kwargs):
        build.mkdir(exist_ok=True)

    monkeypatch.setattr(pixi_build.subprocess, "run", successful_install)
    pixi_build.main()
    assert state.is_file()

    def failing_install(command, **kwargs):
        (build / "partial-object.o").touch()
        raise subprocess.CalledProcessError(1, command)

    monkeypatch.setattr(pixi_build.subprocess, "run", failing_install)
    with pytest.raises(subprocess.CalledProcessError):
        pixi_build.main()
    assert state.is_file()
    assert (build / "partial-object.o").is_file()

    def retry_install(*args, **kwargs):
        assert (build / "partial-object.o").is_file()

    monkeypatch.setattr(pixi_build.subprocess, "run", retry_install)
    pixi_build.main()
    assert state.is_file()
    assert (build / "partial-object.o").is_file()


def test_discovery_policy_change_discards_cached_dependencies(environment, monkeypatch):
    prefix, _ = environment
    build = prefix / "sage-build"
    build.mkdir()
    (build / "old-object.o").touch()
    state = prefix / "sage-build.environment"
    state.write_text(pixi_build.environment_fingerprint(prefix) + "\n")
    original = pixi_build.native_file_contents(prefix)
    monkeypatch.setattr(pixi_build, "native_file_contents", lambda prefix: original + "# new policy\n")
    pixi_build.prepare_build_directory(prefix)
    assert not build.exists()


def test_clean_only_removes_the_owned_build_directory(environment, monkeypatch):
    prefix, record = environment
    build = prefix / "sage-build"
    state = prefix / "sage-build.environment"
    build.mkdir()
    (build / "old-object.o").touch()
    state.write_text(pixi_build.environment_fingerprint(prefix) + "\n")
    neighbour = prefix / "keep-this-file"
    neighbour.touch()
    monkeypatch.setattr(pixi_build.sys, "argv", ["pixi-build.py", "--clean", "--wheel", "--verbose"])

    def install(command, **kwargs):
        assert not build.exists()
        assert "--editable" not in command
        assert "-v" in command
        native = prefix / "sage-build-native.ini"
        assert f"--config-settings=setup-args=--native-file={native}" in command
        assert native.read_text() == pixi_build.native_file_contents(prefix)
        build.mkdir()

    monkeypatch.setattr(pixi_build.subprocess, "run", install)
    pixi_build.main()
    assert state.is_file()
    assert record.is_file()
    assert neighbour.is_file()


def test_refuses_symlink_build_directory(environment, tmp_path):
    prefix, _ = environment
    outside = tmp_path / "outside"
    outside.mkdir()
    sentinel = outside / "keep-this-file"
    sentinel.touch()
    try:
        (prefix / "sage-build").symlink_to(outside, target_is_directory=True)
    except OSError:
        pytest.skip("creating symbolic links is not permitted")
    with pytest.raises(RuntimeError, match="symbolic link"):
        pixi_build.prepare_build_directory(prefix, clean=True)
    assert sentinel.is_file()


def test_cmake_discovery_ignores_host_packages(environment, tmp_path):
    meson = shutil.which("meson")
    cmake = shutil.which("cmake")
    if not meson or not cmake:
        pytest.skip("Meson and CMake are required for the discovery regression test")
    prefix, _ = environment
    host = tmp_path / "host"
    for directory, name in ((prefix, "LockedPackage"), (host, "HostPackage")):
        package = directory / "lib" / "cmake" / name
        package.mkdir(parents=True)
        (package / f"{name}Config.cmake").write_text(
            f"add_library({name}::pkg INTERFACE IMPORTED)\n")
    source = tmp_path / "source"
    source.mkdir()
    (source / "meson.build").write_text("""\
project('discovery-test')
dependency('LockedPackage', method: 'cmake', modules: ['LockedPackage::pkg'])
host = dependency('HostPackage', method: 'cmake', modules: ['HostPackage::pkg'], required: false)
message('host found: ' + host.found().to_string())
""")
    native = tmp_path / "native.ini"
    native.write_text(pixi_build.native_file_contents(prefix))
    # The temporary prefix has only test packages. Use the real CMake executable
    # while retaining the generated dependency discovery policy.
    tool = tmp_path / "tools.ini"
    tool.write_text("[binaries]\ncmake = '" + Path(cmake).as_posix().replace("'", "\\'") + "'\n")
    build_environment = os.environ.copy()
    build_environment["CMAKE_PREFIX_PATH"] = os.pathsep.join((str(prefix), str(host)))
    build_environment["HostPackage_DIR"] = str(host / "lib" / "cmake" / "HostPackage")
    for restricted in (False, True):
        command = [meson, "setup", str(tmp_path / f"build-{restricted}"), str(source)]
        if restricted:
            command.extend(("--native-file", str(native)))
        command.extend(("--native-file", str(tool)))
        result = subprocess.run(command, env=build_environment, capture_output=True, text=True,
                                check=False)
        assert result.returncode == 0, result.stdout + result.stderr
        assert f"host found: {str(not restricted).lower()}" in result.stdout
