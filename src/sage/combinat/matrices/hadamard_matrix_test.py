import pytest


@pytest.mark.longlong
@pytest.mark.parametrize("constructor", ["general", "cooper_wallis"])
def test_hadamard_matrix_2060(constructor):
    """
    Test the order-2060 construction through both entry points.

    Constructing the matrix and checking its orthogonality take a long time.
    """
    from sage.matrix.constructor import matrix

    from sage.combinat.matrices.hadamard_matrix import (
        hadamard_matrix_cooper_wallis_smallcases,
        is_hadamard_matrix,
    )

    if constructor == "general":
        H = matrix.hadamard(2060, check=False)
    else:
        H = hadamard_matrix_cooper_wallis_smallcases(2060, check=False)
    assert H.dimensions() == (2060, 2060)
    assert is_hadamard_matrix(H, normalized=True)
