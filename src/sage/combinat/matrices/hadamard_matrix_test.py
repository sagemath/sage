import pytest


@pytest.mark.parametrize("n", [40, 42, 44])
def test_hadamard_matrix_from_turyn_type_sequences(n):
    """Construct matrices using each new TT sequence from London--Kotsireas."""
    from sage.combinat.t_sequences import (
        T_sequences_construction_from_base_sequences,
        base_sequences_construction,
        turyn_type_sequences_smallcases,
    )
    from sage.combinat.matrices.hadamard_matrix import (
        hadamard_matrix,
        hadamard_matrix_cooper_wallis_construction,
        hadamard_matrix_cooper_wallis_smallcases,
        is_hadamard_matrix,
    )
    from sage.matrix.constructor import matrix

    seqs = turyn_type_sequences_smallcases(n)
    assert list(map(len, seqs)) == [n, n, n, n-1]
    # Check the TT identity directly, independently of the conversions below.
    for shift in range(1, n):
        assert sum(weight * seq[j] * seq[j + shift]
                   for weight, seq in zip((1, 1, 2, 2), seqs)
                   for j in range(len(seq) - shift)) == 0
    base = base_sequences_construction(seqs)
    tseqs = T_sequences_construction_from_base_sequences(base)
    unit = matrix([[1]])
    H = hadamard_matrix_cooper_wallis_construction(
        *tseqs, unit, unit, unit, unit, check=False
    )
    order = 4 * (3*n - 1)
    assert H.dimensions() == (order, order)
    assert is_hadamard_matrix(H)
    assert hadamard_matrix_cooper_wallis_smallcases(order, existence=True)
    assert hadamard_matrix(order, existence=True) is True


@pytest.mark.parametrize("constructor", ["general", "specialized"])
@pytest.mark.parametrize("n,construction", [
    (668, "spence"), (716, "spence"),
    (1676, "spence"), (1772, "spence"),
])
def test_hadamard_matrices_alpoge(n, construction, constructor):
    """Check the four public Spence seeds through both matrix entry points."""
    from sage.combinat.matrices.hadamard_matrix import (
        hadamard_matrix,
        hadamard_matrix_spence_construction,
        is_hadamard_matrix,
    )

    specialized = {
        "spence": hadamard_matrix_spence_construction,
    }
    build = hadamard_matrix if constructor == "general" else specialized[construction]
    assert build(n, existence=True) is True
    # Exercise the default check=True path in the specialized constructors.
    H = build(n, check=constructor == "specialized")
    assert H.dimensions() == (n, n)
    if constructor == "general":
        assert is_hadamard_matrix(H, normalized=True)


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
