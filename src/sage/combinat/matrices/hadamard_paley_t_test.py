import pytest


@pytest.mark.parametrize("n,q,t", [
    (3156, 1051, 3), (5260, 1051, 5), (6036, 2011, 3), (6356, 907, 7),
    (6756, 2251, 3), (7364, 1051, 7), (7660, 1531, 5), (8620, 1723, 5),
    (10060, 2011, 5), (10644, 3547, 3), (10932, 3643, 3), (11260, 2251, 5),
    (11572, 1051, 11), (11804, 907, 13),
])
def test_paley_t_periodic_autocorrelation(n, q, t):
    """Check the exact Goethals--Seidel condition without dense matrices."""
    from sage.combinat.matrices.hadamard_matrix import (
        _paley_t_sequences,
        hadamard_matrix,
    )

    assert hadamard_matrix(n, existence=True) is True
    rows = _paley_t_sequences(q, t)
    assert len(rows) == 4
    m = n // 4
    assert all(len(row) == m and set(row) <= {-1, 1} for row in rows)

    # The correlation at shift k is 4*m minus twice the total Hamming
    # distance from the four cyclically shifted sign vectors.
    mask = (1 << m) - 1
    words = [sum(1 << j for j, x in enumerate(row) if x == 1) for row in rows]
    for k in range(1, m):
        assert sum((w ^ (((w << k) & mask) | (w >> (m - k)))).bit_count()
                   for w in words) == 2 * m


@pytest.mark.longlong
def test_hadamard_matrix_paley_t_3156():
    """Build and verify a representative dense matrix."""
    from sage.combinat.matrices.hadamard_matrix import (
        _hadamard_matrix_paley_t,
        is_hadamard_matrix,
    )

    H = _hadamard_matrix_paley_t(1051, 3)
    assert H.dimensions() == (3156, 3156)
    assert is_hadamard_matrix(H, normalized=True)
