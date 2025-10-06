import pytest


def test_petersen_graph_is_projective_planar():
    r"""
    Ensure that the Petersen graph is projective planar without
    using the cache.
    """
    from sage.graphs.generators.smallgraphs import PetersenGraph
    P = PetersenGraph()
    P.is_projective_planar.clear_cache()
    assert not P.is_projective_planar.is_in_cache()
    assert P.is_projective_planar()
