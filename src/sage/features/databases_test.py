from pathlib import Path

import pytest


@pytest.mark.parametrize(
    ("name", "package", "subdirectory", "filename"),
    [
        ("graphs", "sage_data_graphs.data", "", "graphs.db"),
        ("cremona", "sage_data_elliptic_curves.data", "cremona", "cremona_mini.db"),
        ("ellcurves", "sage_data_elliptic_curves.data", "ellcurves", "rank0"),
        ("reflexive_polytopes", "sage_data_polytopes", "data", "Full3d"),
    ],
)
def test_packaged_data_without_legacy_install(monkeypatch, name, package,
                                             subdirectory, filename):
    from sage.features import StaticFile, databases

    monkeypatch.setattr(databases, "sage_data_paths", lambda name: set())
    paths = databases._data_search_path(name, None, package, subdirectory)
    feature = StaticFile(name=f"packaged_{name}", filename=filename, search_path=tuple(paths))
    assert Path(feature.absolute_filename()).exists()


def test_legacy_data_precedes_package(monkeypatch, tmp_path):
    from sage.features import databases

    monkeypatch.setattr(databases, "sage_data_paths", lambda name: {str(tmp_path)})
    paths = databases._data_search_path("graphs", None, "sage_data_graphs.data")
    assert paths[0] == str(tmp_path)
    assert Path(paths[1], "graphs.db").is_file()


def test_explicit_data_directory(monkeypatch, tmp_path):
    from sage.features import databases

    def unexpected_lookup(*args):
        pytest.fail("an explicit data directory must bypass packaged resources")

    monkeypatch.setattr(databases, "_packaged_data_path", unexpected_lookup)
    assert databases._data_search_path("graphs", str(tmp_path), "sage_data_graphs.data") == str(tmp_path)


def test_missing_data_package_keeps_legacy_paths(monkeypatch, tmp_path):
    from sage.features import databases

    monkeypatch.setattr(databases, "sage_data_paths", lambda name: {str(tmp_path)})
    paths = databases._data_search_path("graphs", None, "nonexistent_sage_data_package")
    assert paths == [str(tmp_path)]


def test_extracted_data_directory_stays_available(monkeypatch, tmp_path):
    import zipfile

    from sage.features import databases

    archive = tmp_path / "resources.zip"
    with zipfile.ZipFile(archive, "w") as stream:
        stream.writestr("data/graphs.db", b"database")
    with zipfile.ZipFile(archive) as stream:
        resource = zipfile.Path(stream, "data/")
        monkeypatch.setattr(databases, "files", lambda package: resource)
        path = databases._packaged_data_path("test_zipped_data", "")
    # The returned directory must outlive the Traversable and its archive.
    assert Path(path, "graphs.db").read_bytes() == b"database"
