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
    (892, "sds"), (1132, "sds"), (1244, "sds"),
    (1388, "bordered"), (1436, "bordered"),
    (1676, "spence"), (1772, "spence"),
    (1916, "bordered"), (1948, "sds"), (1964, "sds"),
])
def test_hadamard_matrices_alpoge(n, construction, constructor):
    """Check all twelve public seeds through both matrix entry points."""
    from sage.combinat.matrices.hadamard_matrix import (
        hadamard_matrix,
        hadamard_matrix_bordered_goethals_seidel,
        hadamard_matrix_from_sds,
        hadamard_matrix_spence_construction,
        is_hadamard_matrix,
    )

    specialized = {
        "spence": hadamard_matrix_spence_construction,
        "sds": hadamard_matrix_from_sds,
        "bordered": hadamard_matrix_bordered_goethals_seidel,
    }
    build = hadamard_matrix if constructor == "general" else specialized[construction]
    assert build(n, existence=True) is True
    # Exercise the default check=True path in the specialized constructors.
    H = build(n, check=constructor == "specialized")
    assert H.dimensions() == (n, n)
    if constructor == "general":
        assert is_hadamard_matrix(H, normalized=True)


def test_hadamard_alpoge_payload():
    """The stored records reproduce the original public payload exactly."""
    from hashlib import sha256
    from sage.combinat.matrices.hadamard_matrix_data import (
        ALPOGE_BORDERED_DATA,
        ALPOGE_GOETHALS_SEIDEL_DATA,
        ALPOGE_SPENCE_DATA,
    )

    records = {n: ''.join(rows) for n, rows in ALPOGE_SPENCE_DATA.items()}
    records.update({n: ''.join(rows) for n, rows in ALPOGE_GOETHALS_SEIDEL_DATA.items()})
    records.update({n: ''.join(data) for n, (_, _, _, data) in ALPOGE_BORDERED_DATA.items()})
    payload = ''.join(records[n] for n in sorted(records))
    assert len(payload) == 23828
    assert sha256(payload.encode('ascii')).hexdigest() == (
        '5b5fe8fa42f0d6a8b4e4c9926726d82a6aab8e1070c1ae4d1b430c1277e58db4'
    )


@pytest.mark.parametrize("n,digest", [
    (1388, 'a6b92584eb803b87026709d64fe892dec8f7182a120e13de9edd3065cf05bf0b'),
    (1436, 'e4d745a4d44f39a5671f9cd86f5c1d0aef93504dcfb2e253451cadc9e4086728'),
    (1916, 'be2073eeaa5399cfe104023829d2c6770b49dd2f07bf6347203f1cbd75577ae9'),
])
def test_hadamard_alpoge_decoder_conventions(n, digest):
    """Match the audited public decoder's larger-border outputs byte for byte."""
    from hashlib import sha256
    from sage.combinat.matrices.hadamard_matrix import hadamard_matrix_bordered_goethals_seidel

    H = hadamard_matrix_bordered_goethals_seidel(n, check=False)
    fingerprint = sha256()
    for row in H.rows():
        fingerprint.update((''.join('+' if c == 1 else '-' for c in row) + '\n').encode('ascii'))
    assert fingerprint.hexdigest() == digest


@pytest.mark.parametrize("n", [0, -4, 10, 668, 1389])
def test_hadamard_bordered_unavailable(n):
    """Availability queries and construction agree for unsupported orders."""
    from sage.combinat.matrices.hadamard_matrix import hadamard_matrix_bordered_goethals_seidel

    assert hadamard_matrix_bordered_goethals_seidel(n, existence=True) is False
    with pytest.raises(ValueError, match=f'No larger-border Goethals-Seidel data for order {n}'):
        hadamard_matrix_bordered_goethals_seidel(n)


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
