import pytest


def test_shortened_000_111_extended_binary_Golay_code_graph():
    r"""
    Test that Sage produces a graph equal to the one that we get
    from this construction.

    The construction itself takes a long time.
    """
    from sage.coding import codes_catalog
    from sage.coding.linear_code import LinearCode
    from sage.graphs.generators.distance_regular import (
        shortened_000_111_extended_binary_Golay_code_graph
    )
    from sage.matrix.constructor import matrix
    from sage.rings.finite_rings.finite_field_constructor import FiniteField

    code = codes_catalog.GolayCode(FiniteField(2))
    C_basis = code.basis()

    # now special shortening
    v = C_basis[0] + C_basis[1] + C_basis[2]  # v has 111 at the start
    C_basis = C_basis[3:]
    C_basis.append(v)
    C_basis = [x[3:] for x in C_basis]

    code = LinearCode(matrix(FiniteField(2), C_basis))
    G = code.cosetGraph()
    G.name("Shortened 000 111 extended binary Golay code")
    assert G.is_distance_regular()

    H = shortened_000_111_extended_binary_Golay_code_graph()
    assert G == H
