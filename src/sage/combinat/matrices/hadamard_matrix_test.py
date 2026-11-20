import pytest


@pytest.mark.parametrize("n", [40, 42, 44])
def test_hadamard_matrix_from_turyn_type_sequences(n):
    """Construct matrices using each new TT sequence from London--Kotsireas."""
    from sage.combinat.matrices.hadamard_matrix import (
        hadamard_matrix,
        hadamard_matrix_cooper_wallis_construction,
        hadamard_matrix_cooper_wallis_smallcases,
        is_hadamard_matrix,
    )
    from sage.combinat.t_sequences import (
        T_sequences_construction_from_base_sequences,
        base_sequences_construction,
        turyn_type_sequences_smallcases,
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


@pytest.mark.longlong
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


@pytest.mark.longlong
@pytest.mark.parametrize("n,digest", [
    (1388, 'a6b92584eb803b87026709d64fe892dec8f7182a120e13de9edd3065cf05bf0b'),
    (1436, 'e4d745a4d44f39a5671f9cd86f5c1d0aef93504dcfb2e253451cadc9e4086728'),
    (1916, 'be2073eeaa5399cfe104023829d2c6770b49dd2f07bf6347203f1cbd75577ae9'),
])
def test_hadamard_alpoge_decoder_conventions(n, digest):
    """Match the audited public decoder's larger-border outputs byte for byte."""
    from hashlib import sha256

    from sage.combinat.matrices.hadamard_matrix import (
        hadamard_matrix_bordered_goethals_seidel,
    )

    H = hadamard_matrix_bordered_goethals_seidel(n, check=False)
    fingerprint = sha256()
    for row in H.rows():
        fingerprint.update((''.join('+' if c == 1 else '-' for c in row) + '\n').encode('ascii'))
    assert fingerprint.hexdigest() == digest


@pytest.mark.parametrize("n", [0, -4, 10, 668, 1389])
def test_hadamard_bordered_unavailable(n):
    """Availability queries and construction agree for unsupported orders."""
    from sage.combinat.matrices.hadamard_matrix import (
        hadamard_matrix_bordered_goethals_seidel,
    )

    assert hadamard_matrix_bordered_goethals_seidel(n, existence=True) is False
    with pytest.raises(ValueError, match=f'No larger-border Goethals-Seidel data for order {n}'):
        hadamard_matrix_bordered_goethals_seidel(n)


@pytest.mark.longlong
def test_hadamard_matrix_2620_normalized():
    """Normalize Cooper--Wallis matrices enabled by the length-44 TT seed."""
    from sage.combinat.matrices.hadamard_matrix import (
        hadamard_matrix,
        is_hadamard_matrix,
    )

    assert hadamard_matrix(2620, existence=True, construction_name=True) == (
        'cooper-wallis small cases: 2620'
    )
    H = hadamard_matrix(2620, check=False)
    assert H.dimensions() == (2620, 2620)
    assert is_hadamard_matrix(H, normalized=True)


@pytest.mark.longlong
@pytest.mark.parametrize("constructor", ["general", "cooper_wallis"])
def test_hadamard_matrix_2060(constructor):
    """
    Test the order-2060 construction through both entry points.

    Constructing the matrix and checking its orthogonality take a long time.
    """
    from sage.combinat.matrices.hadamard_matrix import (
        hadamard_matrix_cooper_wallis_smallcases,
        is_hadamard_matrix,
    )
    from sage.matrix.constructor import matrix

    if constructor == "general":
        H = matrix.hadamard(2060, check=False)
    else:
        H = hadamard_matrix_cooper_wallis_smallcases(2060, check=False)
    assert H.dimensions() == (2060, 2060)
    assert is_hadamard_matrix(H, normalized=True)


def _assert_hadamard_by_hamming(H):
    """Check a positive square sign matrix using exact row-pair distances."""
    n = H.nrows()
    assert n > 0 and H.ncols() == n
    packed = []
    for row in H:
        entries = row.list()
        assert set(entries) <= {-1, 1}
        packed.append(int(''.join('1' if c == 1 else '0' for c in entries), 2))
    # For sign rows, the inner product is n minus twice the Hamming distance.
    # Each row has squared norm n, and every distinct pair must be orthogonal.
    for j, row in enumerate(packed):
        assert all((row ^ other).bit_count() == n // 2 for other in packed[:j])


def test_hadamard_hamming_verifier():
    """Check the exact verifier against valid and corrupted small matrices."""
    from sage.matrix.constructor import matrix

    _assert_hadamard_by_hamming(matrix([[1, 1], [1, -1]]))
    for rows in ([[1, 1], [1, 1]], [[1, 1], [1, 0]], [[1, 1, 1], [1, -1, 1]]):
        with pytest.raises(AssertionError):
            _assert_hadamard_by_hamming(matrix(rows))


@pytest.mark.parametrize(
    "n,digest",
    [
        (
            3732,
            "81f305ef4f38b4ec7d70d1316ba27628847cb278d8c31366319052aec593aaa9",
        ),
        (
            4492,
            "b694d19734a5eadd4d21fbbb62de95db6ccbfce42f30b2a1c16917dfc883ee6e",
        ),
        (
            5020,
            "abc51153c01dc99f2a0652546ad65e1cbf3dc55b796f5340848ac3ecf34c1d3b",
        ),
        (
            5932,
            "19c355158e4c4bf4a12c1305fa0442ab750fac24316ef1d579423618308e682c",
        ),
    ],
)
def test_hadamard_eklo_seeds(n, digest):
    """Check the published seeds exactly, independently of the matrix array."""
    from hashlib import sha256
    from math import gcd

    from sage.combinat.matrices.hadamard_matrix import (
        hadamard_matrix,
        hadamard_matrix_from_sds,
    )
    from sage.combinat.matrices.hadamard_matrix_data import (
        EKLO_GOETHALS_SEIDEL_DATA,
    )

    t = n // 4
    multipliers, representatives = EKLO_GOETHALS_SEIDEL_DATA[n]
    assert all(gcd(u, t) == 1 for u in multipliers)
    assert {u * v % t for u in multipliers for v in multipliers} == set(
        multipliers
    )
    supports = [
        {u * r % t for r in reps for u in multipliers}
        for reps in representatives
    ]
    rows = [
        [-1 if i in support else 1 for i in range(t)] for support in supports
    ]
    payload = ''.join('+' if a == 1 else '-' for row in rows for a in row)
    assert sha256(payload.encode('ascii')).hexdigest() == digest
    # Use integer sums, with no floating-point Fourier approximation.
    for shift in range(1, t):
        assert (
            sum(
                row[i] * row[(i + shift) % t] for row in rows for i in range(t)
            )
            == 0
        )
    assert hadamard_matrix_from_sds(n, existence=True) is True
    assert hadamard_matrix(n, existence=True) is True
    assert (
        hadamard_matrix(n, existence=True, construction_name=True)
        == f'SDS {n}'
    )


@pytest.mark.longlong
@pytest.mark.parametrize("n", [3732, 4492, 5020, 5932])
def test_hadamard_eklo_matrices(n):
    """Verify full matrices and exercise the general dispatcher."""
    from sage.combinat.matrices.hadamard_matrix import (
        hadamard_matrix,
        hadamard_matrix_from_sds,
        normalise_hadamard,
    )

    # The default check=True verifies the core identity and full Gram matrix.
    H = hadamard_matrix_from_sds(n)
    assert H.dimensions() == (n, n)
    assert hadamard_matrix(n, check=False) == normalise_hadamard(H)


def test_hadamard_bordered_5788_data():
    """Verify group autocorrelations and the Theorem 7 border identities."""
    from hashlib import sha256
    from itertools import product

    from sage.combinat.matrices.hadamard_matrix import (
        hadamard_matrix,
        hadamard_matrix_bordered_goethals_seidel,
    )
    from sage.combinat.matrices.hadamard_matrix_data import EKLO_BORDERED_DATA
    from sage.matrix.constructor import (
        block_matrix,
        identity_matrix,
        matrix,
        ones_matrix,
    )
    from sage.rings.integer_ring import ZZ

    factors, border, orbit, data = EKLO_BORDERED_DATA[5788]
    assert (factors, border, orbit) == ((2, 2, 19, 19), 12, 361)
    assert sha256(''.join(data).encode('ascii')).hexdigest() == (
        'f39d20d62fa95e936e948c2b508cebc3791675cd5307e27f355d86a2ca58a5ef'
    )
    assert set(''.join(data)) == {'+', '-'}
    size = 4 * orbit
    assert (
        list(map(len, data))
        == [4 * size, border * (border + 16)] + [4 * border] * 4
    )
    signs = [[1 if c == '+' else -1 for c in part] for part in data]
    rows = [signs[0][i * size : (i + 1) * size] for i in range(4)]
    coordinates = list(product(range(2), range(2), range(19), range(19)))
    for a, b, c, d in coordinates[1:]:
        shifted = [
            ((x ^ a) * 2 + (y ^ b)) * 361 + ((z + c) % 19) * 19 + (w + d) % 19
            for x, y, z, w in coordinates
        ]
        correlation = sum(
            row[i] * row[j] for row in rows for i, j in enumerate(shifted)
        )
        assert correlation == (-12 if a == b == 0 else 4)

    top = matrix(ZZ, border, border + 16, signs[1])
    K = top[:, :border]
    E = top[:, border:]
    F = matrix(ZZ, 16, border, sum(signs[2:], []))
    assert K * K.T == 12 * identity_matrix(ZZ, 12)
    assert E * E.T == 16 * identity_matrix(ZZ, 12)
    assert F * F.T == 16 * identity_matrix(ZZ, 16) - 4 * (
        identity_matrix(ZZ, 4).tensor_product(ones_matrix(ZZ, 4))
    )
    # The core acts on coset-constant vectors through this 16-by-16 quotient.
    sums = [
        [sum(row[i * orbit : (i + 1) * orbit]) for i in range(4)]
        for row in rows
    ]
    A, B, C, D = [matrix(ZZ, 4, 4, lambda i, j: row[i ^ j]) for row in sums]
    Q = block_matrix(
        [
            [A, B, C, D],
            [-B, A, D.T, -C.T],
            [-C, -D.T, A, B.T],
            [-D, C.T, -B.T, A],
        ]
    )
    assert Q * E.T + F * K.T == 0
    assert (
        hadamard_matrix_bordered_goethals_seidel(5788, existence=True) is True
    )
    assert hadamard_matrix(5788, existence=True) is True
    assert hadamard_matrix(5788, existence=True, construction_name=True) == (
        'bordered Goethals-Seidel 5788'
    )


@pytest.mark.longlong
def test_hadamard_bordered_5788_matrix():
    """Check full orthogonality, the source image fingerprint, and dispatch."""
    from hashlib import sha256

    from sage.combinat.matrices.hadamard_matrix import (
        hadamard_matrix,
        hadamard_matrix_bordered_goethals_seidel,
        normalise_hadamard,
    )

    # The default check=True computes the full exact product H H^T.
    H = hadamard_matrix_bordered_goethals_seidel(5788)
    assert H.dimensions() == (5788, 5788)
    fingerprint = sha256()
    for row in H.rows():
        fingerprint.update(
            (''.join('+' if c == 1 else '-' for c in row) + '\n').encode(
                'ascii'
            )
        )
    assert fingerprint.hexdigest() == (
        'dd5893a18e4fe0f8ab3a36df163e0e5643a64deb017726300f37ce4f9761b06d'
    )
    normalise_hadamard(H)
    assert set(H.row(0)) == set(H.column(0)) == {1}
    assert hadamard_matrix(5788, check=False) == H


def test_hadamard_bh86_seed_and_conversion():
    """Verify the recovered BH data and all Williamson-type identities."""
    from hashlib import sha256

    from sage.combinat.matrices.hadamard_matrix import (
        hadamard_matrix,
        hadamard_matrix_cooper_wallis_smallcases,
        symmetric_butson_hadamard_matrix_smallcases,
        williamson_type_from_symmetric_butson,
    )
    from sage.matrix.constructor import identity_matrix

    H = symmetric_butson_hadamard_matrix_smallcases(86)
    assert H.dimensions() == (86, 86)
    assert H == H.T
    assert H * H.conjugate_transpose() == 86 * identity_matrix(
        H.base_ring(), 86
    )
    i = H.base_ring().gen()
    symbols = {H.base_ring()(1): '+', H.base_ring()(-1): '-', i: 'i', -i: 'j'}
    payload = ''.join(symbols[c] for c in H.list())
    assert sha256(payload.encode('ascii')).hexdigest() == (
        '15da4d8dcaacbecd3f7128c44dcef8a2cc76d14b69e036c7d47271080989289e'
    )
    W = williamson_type_from_symmetric_butson(H)
    assert all(A.dimensions() == (173, 173) and A == A.T for A in W)
    assert all(set(A.list()) == {-1, 1} for A in W)
    assert all(A * B.T == B * A.T for j, A in enumerate(W) for B in W[j + 1 :])
    assert sum(A * A.T for A in W) == 692 * identity_matrix(173)
    for t in (3, 5, 7, 9, 11, 13):
        assert (
            hadamard_matrix_cooper_wallis_smallcases(692 * t, existence=True)
            is True
        )
        assert hadamard_matrix(
            692 * t, existence=True, construction_name=True
        ) == (f'cooper-wallis small cases: {692 * t}')


@pytest.mark.parametrize('n', [0, -2, 42, 84, 88])
def test_hadamard_bh_unavailable(n):
    from sage.combinat.matrices.hadamard_matrix import (
        symmetric_butson_hadamard_matrix_smallcases,
    )

    assert (
        symmetric_butson_hadamard_matrix_smallcases(n, existence=True) is False
    )
    with pytest.raises(ValueError, match='No symmetric BH'):
        symmetric_butson_hadamard_matrix_smallcases(n)


@pytest.mark.parametrize(
    'rows,message',
    [
        ([[1, 1, 1]], 'positive even order'),
        ([[1, 1], [-1, 1]], 'must be symmetric'),
        ([[2, 1], [1, 2]], 'fourth-root entries'),
        ([[1, 1], [1, 1]], 'must satisfy'),
    ],
)
def test_hadamard_bh_invalid_input(rows, message):
    from sage.combinat.matrices.hadamard_matrix import (
        williamson_type_from_symmetric_butson,
    )
    from sage.matrix.constructor import matrix

    with pytest.raises(ValueError, match=message):
        williamson_type_from_symmetric_butson(matrix(rows))


def test_hadamard_bh86_small_prime_power():
    """Exercise the general conversion over a non-prime finite field (q=9)."""
    from sage.combinat.matrices.hadamard_matrix import (
        hadamard_matrix,
        williamson_type_from_symmetric_butson,
    )
    from sage.matrix.constructor import identity_matrix

    W = williamson_type_from_symmetric_butson(hadamard_matrix(4))
    assert all(A.dimensions() == (9, 9) for A in W)
    assert sum(A * A.T for A in W) == 36 * identity_matrix(9)


@pytest.mark.longlong
@pytest.mark.parametrize('n', [2076, 3460, 4844, 6228, 7612, 8996])
def test_hadamard_bh86_full_matrices(n):
    """Verify all six full matrices by exact row-pair Hamming distances."""
    from sage.combinat.matrices.hadamard_matrix import (
        hadamard_matrix,
        hadamard_matrix_cooper_wallis_smallcases,
    )

    # Exercise default check=True (integer Gram product) at the first order.
    H = (
        hadamard_matrix_cooper_wallis_smallcases(n)
        if n == 2076
        else hadamard_matrix(n, check=False)
    )
    assert H.dimensions() == (n, n)
    assert set(H.row(0)) == set(H.column(0)) == {1}
    _assert_hadamard_by_hamming(H)
