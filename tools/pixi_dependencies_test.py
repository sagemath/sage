"""Checks for generated Pixi environments and platform dependency constraints."""

import importlib.util
import sys
import tomllib
from pathlib import Path

import pytest

TOOLS = Path(__file__).resolve().parent
sys.path.insert(0, str(TOOLS))
spec = importlib.util.spec_from_file_location("update_pixi", TOOLS / "update-pixi.py")
pixi = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pixi)


@pytest.mark.parametrize(("requirement", "expected"), [
    ("numpy >=2.3, <3", ("numpy", {"version": ">=2.3,<3"})),
    ("singular ==4.4.1.p5", ("singular", {"version": "==4.4.1.p5"})),
    ("libblas=*=*_newaccelerate", ("libblas", {"version": "*", "build": "*_newaccelerate"})),
    ("numpy=2.3=py313_0", ("numpy", {"version": "==2.3", "build": "py313_0"})),
    ("numpy=2.3", ("numpy", {"version": "2.3.*"})),
    ("python-lrcalc~=2.1", ("python-lrcalc", {"version": "~=2.1"})),
])
def test_version_and_build_constraints(requirement, expected):
    assert pixi.parse_requirement(requirement) == expected


def test_combines_constraints_and_keeps_build():
    assert pixi.merge_requirements({"numpy", "numpy >=2", "numpy <3"}) == {
        "numpy": {"version": "<3,>=2"}}
    assert pixi.merge_requirements({"libblas >=3.9", "libblas=*=*_openblas"}) == {
        "libblas": {"version": ">=3.9", "build": "*_openblas"}}
    with pytest.raises(ValueError, match="Incompatible build"):
        pixi.merge_requirements({"libblas=*=*_openblas", "libblas=*=*_newaccelerate"})
    with pytest.raises(ValueError, match="Unsupported dependency"):
        pixi.parse_requirement("conda-forge::numpy[build=py313_0]")


def environment_dependencies(manifest, python, platform):
    """Combine the actual TOML tables the way Pixi selects a feature/target."""
    dependencies = dict(manifest.get("dependencies", {}))
    dependencies.update(manifest.get("target", {}).get(platform, {}).get("dependencies", {}))
    feature = manifest.get("feature", {}).get("py" + python.replace(".", ""), {})
    dependencies.update(feature.get("dependencies", {}))
    dependencies.update(feature.get("target", {}).get(platform, {}).get("dependencies", {}))
    return {name: {"version": value} if isinstance(value, str) else value
            for name, value in dependencies.items()}


def test_generated_environments_preserve_all_dependencies():
    source = TOOLS.parent / "pyproject.toml"
    manifest = tomllib.loads(pixi.render_dependencies(pixi.dependency_tables(source)))
    for python in pixi.PYTHONS:
        for platform in pixi.PLATFORMS:
            resolved = environment_dependencies(manifest, python, platform)
            expected = pixi.merge_requirements(pixi.get_dependencies(source, python, platform))
            assert resolved == expected, (python, platform)
            if platform == "win-64":
                assert "singular" not in resolved
            else:
                assert resolved["singular"]["version"] == "==4.4.1.p5"
            blas = "*_newaccelerate" if platform.startswith("osx") else "*_openblas"
            assert resolved["libblas"]["build"] == blas
            assert ("c-compiler" in resolved) == (platform != "win-64")
            assert ("vs2022_win-64" in resolved) == (platform == "win-64")
            assert ("4ti2" in resolved) == (platform in {"linux-64", "osx-64", "win-64"})
            assert "python" not in resolved


def test_python_dependent_requirements_are_kept_separate(monkeypatch):
    def dependencies(source, python, platform):
        result = {"python", "shared", f"numpy=={1 if python == '3.13' else 2}"}
        if platform == "win-64":
            result.add("windows-only")
        if python == "3.14" and platform == "linux-64":
            result.add("linux-py314")
        return result

    monkeypatch.setattr(pixi, "get_dependencies", dependencies)
    manifest = tomllib.loads(pixi.render_dependencies(pixi.dependency_tables(Path("unused"))))
    for python in pixi.PYTHONS:
        for platform in pixi.PLATFORMS:
            actual = environment_dependencies(manifest, python, platform)
            assert actual["numpy"] == {"version": "==1" if python == "3.13" else "==2"}
            assert ("windows-only" in actual) == (platform == "win-64")
            assert ("linux-py314" in actual) == (python == "3.14" and platform == "linux-64")


def test_check_detects_drift_without_writing_and_preserves_prefix(tmp_path, monkeypatch):
    manifest = tmp_path / "pixi.toml"
    prefix = '[workspace]\nname = "sage-dev"\n\n[tasks]\nbuild = "keep this task"\n\n'
    original = prefix + pixi.MARKER + '\n\n[dependencies]\nnumpy = "<2"\n'
    manifest.write_text(original)
    monkeypatch.setattr(pixi, "dependency_tables", lambda source: {
        "dependencies": {"numpy": {"version": ">=2.3"}}})
    assert not pixi.update_pixi(tmp_path, check=True)
    assert manifest.read_text() == original
    assert pixi.update_pixi(tmp_path)
    assert manifest.read_text().startswith(prefix + pixi.MARKER)
    assert tomllib.loads(manifest.read_text())["dependencies"]["numpy"] == ">=2.3"
    assert pixi.update_pixi(tmp_path, check=True)
