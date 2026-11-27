from importlib.resources import files
from pathlib import Path

import pytest

STANDARD_DATABASES = [
    ("DatabaseGraphs", (), "graph_databases", "data", "graphs.db"),
    ("DatabaseCremona", ("cremona_mini",), "elliptic_curve_databases", "data/cremona", "cremona_mini.db"),
    ("DatabaseEllcurves", (), "elliptic_curve_databases", "data/ellcurves", "rank0"),
    ("DatabaseReflexivePolytopes", (), "polytope_databases", "data", "Full3d"),
]


@pytest.fixture(autouse=True)
def fresh_features(monkeypatch):
    import sage.features

    # Do not retain feature results or temporary search paths between tests.
    monkeypatch.setattr(sage.features, "_trivial_unique_representation_cache", {})


@pytest.mark.parametrize(("factory", "args", "package", "directory", "filename"), STANDARD_DATABASES)
def test_presence_only_imports_package(monkeypatch, factory, args, package,
                                      directory, filename):
    from sage.features import databases

    def unexpected_lookup(*args):
        pytest.fail("checking package presence must not search or extract data")

    monkeypatch.setattr(databases, "sage_data_paths", unexpected_lookup)
    monkeypatch.setattr(databases, "_packaged_data_path", unexpected_lookup)
    feature = getattr(databases, factory)(*args)
    assert feature.is_present()


@pytest.mark.parametrize(("factory", "args", "package", "directory", "filename"), STANDARD_DATABASES)
def test_legacy_files_cannot_replace_missing_package(monkeypatch, tmp_path, factory,
                                                  args, package, directory, filename):
    import importlib

    from sage.features import databases

    (tmp_path / filename).touch()
    monkeypatch.setattr(databases, "sage_data_paths", lambda name: {str(tmp_path)})
    import_module = importlib.import_module

    def missing_package(name, *args, **kwds):
        if name == package:
            raise ModuleNotFoundError(f"No module named '{name}'")
        return import_module(name, *args, **kwds)

    monkeypatch.setattr(importlib, "import_module", missing_package)
    feature = getattr(databases, factory)(*args)
    result = feature.is_present()
    assert not result
    assert package in result.reason


@pytest.mark.parametrize(("factory", "args", "package", "directory", "filename"), STANDARD_DATABASES)
def test_packaged_data_ignores_legacy_directories(monkeypatch, tmp_path, factory,
                                               args, package, directory, filename):
    import sage.env
    from sage.features import databases

    (tmp_path / filename).touch()
    monkeypatch.setattr(databases, "sage_data_paths", lambda name: {str(tmp_path)})
    for name in ("CREMONA_MINI_DATA_DIR", "CREMONA_LARGE_DATA_DIR", "ELLCURVE_DATA_DIR",
                 "GRAPHS_DATA_DIR", "POLYTOPE_DATA_DIR"):
        monkeypatch.setenv(name, str(tmp_path))
        monkeypatch.setattr(sage.env, name, str(tmp_path), raising=False)
    feature = getattr(databases, factory)(*args)
    path = Path(feature.absolute_filename())
    assert feature.filename == filename
    assert path.exists()
    assert path == files(package).joinpath(directory, filename)


@pytest.mark.parametrize("explicit_override", [False, True])
@pytest.mark.parametrize(
    ("factory", "args", "variable", "filename"),
    [
        ("DatabaseCremona", (), "CREMONA_LARGE_DATA_DIR", "cremona.db"),
        ("DatabaseReflexivePolytopes", ("polytopes_db_4d",), "POLYTOPE_DATA_DIR", "Hodge4d"),
    ],
)
def test_optional_database_directories(monkeypatch, tmp_path, factory, args,
                                      variable, filename, explicit_override):
    import sage.env
    from sage.features import databases

    path = tmp_path / filename
    if filename == "Hodge4d":
        path.mkdir()
    else:
        path.touch()
    monkeypatch.setattr(sage.env, variable, str(tmp_path) if explicit_override else None)
    monkeypatch.setattr(databases, "sage_data_paths",
                        lambda name: set() if explicit_override else {str(tmp_path)})
    feature = getattr(databases, factory)(*args)
    assert feature.is_present()
    assert Path(feature.absolute_filename()) == path


def test_extracted_data_directory_stays_available(monkeypatch, tmp_path):
    import zipfile

    from sage.features import databases

    archive = tmp_path / "resources.zip"
    with zipfile.ZipFile(archive, "w") as stream:
        stream.writestr("data/graphs.db", b"database")
        stream.writestr("data/isgci_sage.xml", b"sibling resource")
    with zipfile.ZipFile(archive) as stream:
        resource = zipfile.Path(stream, "data/")
        monkeypatch.setattr(databases, "files", lambda package: resource)
        path = databases._packaged_data_path("test_zipped_data", "")
    # The whole directory must outlive the Traversable and its archive.
    assert Path(path, "graphs.db").read_bytes() == b"database"
    assert Path(path, "isgci_sage.xml").read_bytes() == b"sibling resource"
