from pathlib import Path

import pytest


@pytest.fixture
def graph_classes():
    """Restore the singleton's cached data after each test."""
    from sage.graphs.isgci import graph_classes

    methods = (graph_classes.classes, graph_classes.inclusions,
               graph_classes.smallgraphs, graph_classes.inclusion_digraph)
    saved = [(method, method() if method.is_in_cache() else None)
             for method in methods]
    for method in methods:
        method.clear_cache()
    try:
        yield graph_classes
    finally:
        for method, value in saved:
            method.clear_cache()
            if value is not None:
                method.set_cache(value)


@pytest.mark.parametrize("explicit_override", [False, True])
def test_parse_ignores_user_data(monkeypatch, tmp_path, graph_classes,
                                explicit_override):
    from graph_databases import resource_path

    import sage.env
    from sage.graphs import isgci

    user_data = tmp_path / "db" / "graphs"
    legacy_data = tmp_path / "legacy"
    for directory in (user_data, legacy_data):
        directory.mkdir(parents=True)
        for name in ("graphs.db", isgci._XML_FILE, isgci._SMALLGRAPHS_FILE):
            (directory / name).write_text("invalid stale database")
    monkeypatch.setattr(sage.env, "DOT_SAGE", str(tmp_path))
    # Simulate both the old environment variable and its former Sage setting.
    override = str(legacy_data) if explicit_override else None
    if explicit_override:
        monkeypatch.setenv("GRAPHS_DATA_DIR", override)
    else:
        monkeypatch.delenv("GRAPHS_DATA_DIR", raising=False)
    monkeypatch.setattr(sage.env, "GRAPHS_DATA_DIR", override, raising=False)

    with resource_path("graphs.db") as packaged:
        assert Path(isgci.DatabaseGraphs().absolute_filename()).resolve() == packaged.resolve()
    graph_classes._parse_db()
    assert graph_classes.classes()["gc_151"]["name"] == "cograph"
    assert graph_classes.classes()["gc_151"]["problem"]["Clique"] == {"complexity": "Linear"}
    assert graph_classes.smallgraphs()["bull"].order() == 5
    assert graph_classes.Tree <= graph_classes.Chordal


def test_update_only_warns(monkeypatch, tmp_path, graph_classes):
    import urllib.request

    import sage.env
    from sage.graphs import isgci

    def unexpected_access(*args, **kwds):
        pytest.fail("the deprecated update method must not access a database or the network")

    monkeypatch.setattr(urllib.request, "urlopen", unexpected_access)
    monkeypatch.setattr(isgci, "DatabaseGraphs", unexpected_access)
    monkeypatch.setattr(sage.env, "DOT_SAGE", str(tmp_path))
    methods = (graph_classes.classes, graph_classes.inclusions,
               graph_classes.smallgraphs, graph_classes.inclusion_digraph)
    values = [object() for method in methods]
    for method, value in zip(methods, values):
        method.set_cache(value)

    with pytest.warns(DeprecationWarning, match="Update the graph-databases package"):
        graph_classes.update_db()

    assert all(method.is_in_cache() and method() is value
               for method, value in zip(methods, values))
    assert not list(tmp_path.iterdir())
