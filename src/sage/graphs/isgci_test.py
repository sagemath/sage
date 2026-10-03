import io
import zipfile
from types import SimpleNamespace

import pytest


def test_download_uses_user_data_directory(monkeypatch, tmp_path):
    import sage.env
    from sage.graphs import isgci

    archive = io.BytesIO()
    contents = {
        isgci._XML_FILE: b"updated XML",
        isgci._SMALLGRAPHS_FILE: b"updated small graphs",
    }
    with zipfile.ZipFile(archive, "w") as stream:
        for name, content in contents.items():
            stream.writestr(name, content)
    monkeypatch.setattr(sage.env, "DOT_SAGE", str(tmp_path))
    monkeypatch.setattr(isgci, "urlopen", lambda *args, **kwds: io.BytesIO(archive.getvalue()))
    # Updating must not ask for a writable installed database directory.
    def unexpected_lookup():
        pytest.fail("download must not modify installed package data")
    monkeypatch.setattr(isgci, "DatabaseGraphs", unexpected_lookup)

    isgci.graph_classes._download_db()
    for name, content in contents.items():
        assert (tmp_path / "db" / "graphs" / name).read_bytes() == content


@pytest.mark.parametrize("explicit_override", [False, True])
def test_parse_respects_user_updates_and_explicit_override(monkeypatch, tmp_path,
                                                         explicit_override):
    import xml.etree.ElementTree as ET

    import sage.env
    from sage.graphs import isgci

    installed = tmp_path / "installed"
    updated = tmp_path / "db" / "graphs"
    updated.mkdir(parents=True)
    for name in (isgci._XML_FILE, isgci._SMALLGRAPHS_FILE):
        (updated / name).touch()
    monkeypatch.setattr(sage.env, "DOT_SAGE", str(tmp_path))
    monkeypatch.setattr(sage.env, "GRAPHS_DATA_DIR", str(installed) if explicit_override else None)
    monkeypatch.setattr(isgci, "DatabaseGraphs", lambda: SimpleNamespace(
        absolute_filename=lambda: str(installed / "graphs.db")))
    selected = []

    def capture_input(*, file):
        selected.append(file)
        raise RuntimeError("captured XML path")

    monkeypatch.setattr(ET, "ElementTree", capture_input)
    with pytest.raises(RuntimeError, match="captured XML path"):
        isgci.graph_classes._parse_db()
    expected = installed if explicit_override else updated
    assert selected == [str(expected / isgci._XML_FILE)]
