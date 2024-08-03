r"""
Generic Convolution

Asymptotically fast convolution of lists over any commutative ring
in which the multiply-by-two map is injective. (More precisely, if
`x \in R`, and `x = 2^k*y` for some `k \geq 0`, we require that
`R(x/2^k)` returns `y`.)

The main function to be exported is :func:`convolution`.

EXAMPLES::

    sage: convolution([1, 2, 3, 4, 5], [6, 7])
    [6, 19, 32, 45, 58, 35]

The convolution function is reasonably fast, even though it is written
in pure Python. For example, the following takes less than a second::

    sage: v = convolution(list(range(1000)), list(range(1000)))

ALGORITHM:

Converts the problem to multiplication in the ring
`S[x]/(x^M - 1)`, where `S = R[y]/(y^K + 1)` (where
`R` is the original base ring). Performs FFT with respect
to the roots of unity `1, y, y^2, \ldots, y^{2K-1}` in
`S`. The FFT/IFFT are accomplished with just additions and
subtractions and rotating python lists. (I think this algorithm is
essentially due to Schonhage, not completely sure.) The pointwise
multiplications are handled recursively, switching to a classical
algorithm at some point.

Complexity is O(n log(n) log(log(n))) additions/subtractions in R
and O(n log(n)) multiplications in R.

AUTHORS:

- David Harvey (2007-07): first implementation

- William Stein: editing the docstrings for inclusion in Sage.
"""
# ****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#                          David Harvey <dmharvey@math.harvard.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.structure.all import parent
from math import log, ceil

# -------------------------------------------------------------------
#      main exported routine


def convolution(L1, L2):
    """
    Return convolution of non-empty lists L1 and L2.

    L1 and L2 may have arbitrary lengths.

    EXAMPLES::

        sage: convolution([1, 2, 3], [4, 5, 6, 7])
        [4, 13, 28, 34, 32, 21]

    TESTS::

        sage: R = Integers(47)
        sage: L1 = [R.random_element() for _ in range(1000)]
        sage: L2 = [R.random_element() for _ in range(3756)]
        sage: L3 = convolution(L1, L2)
        sage: L3[2000] == sum([L1[i] * L2[2000-i] for i in range(1000)])
        True
        sage: len(L3) == 1000 + 3756 - 1
        True
    """
    if not L1 or not L2:
        raise ValueError("cannot compute convolution of empty lists")
    if len(L1) <= 100 and len(L2) <= 100:   # very arbitrary cutoff
        return _convolution_naive(L1, L2)
    return _convolution_fft(L1, L2)


# -------------------------------------------------------------------
#      naive convolution and negacyclic convolutions


def _convolution_naive(L1, L2):
    """
    Return convolution of non-empty lists L1 and L2, using naive algorithm.

    L1 and L2 may have arbitrary lengths.

    EXAMPLES::

        sage: from sage.rings.polynomial.convolution import _convolution_naive
        sage: _convolution_naive([2], [3])
        [6]
        sage: _convolution_naive([2, 5], [3])
        [6, 15]
        sage: _convolution_naive([2], [3, 6])
        [6, 12]
        sage: _convolution_naive([1, 2, 3], [4, 5, 6, 7])
        [4, 13, 28, 34, 32, 21]
        sage: _convolution_naive([4, 5, 6, 7], [1, 2, 3])
        [4, 13, 28, 34, 32, 21]
    """
    assert len(L1) and len(L2)

    m1 = len(L1)
    m2 = len(L2)

    return [sum([L1[i] * L2[k - i]
            for i in range(max(0, k - m2 + 1), min(k + 1, m1))])
            for k in range(m1 + m2 - 1)]


def _negaconvolution_naive(L1, L2):
    """
    Negacyclic convolution of L1 and L2, using naive algorithm.

    L1 and L2 must be the same length.

    EXAMPLES::

        sage: from sage.rings.polynomial.convolution import _negaconvolution_naive
        sage: from sage.rings.polynomial.convolution import _convolution_naive
        sage: _negaconvolution_naive([2], [3])
        [6]
        sage: _convolution_naive([1, 2, 3], [3, 4, 5])
        [3, 10, 22, 22, 15]
        sage: _negaconvolution_naive([1, 2, 3], [3, 4, 5])
        [-19, -5, 22]
    """
    assert len(L1)
    assert len(L1) == len(L2)

    N = len(L1)
    return [sum([L1[i] * L2[j - i] for i in range(j + 1)]) -
            sum([L1[i] * L2[N + j - i]
                 for i in range(j + 1, N)]) for j in range(N)]


# -------------------------------------------------------------------
#      FFT/IFFT routines


def _forward_butterfly(L1, L2, r):
    r"""
    L1 and L2 are both lists of length K, and
    `0 \leq r \leq K`. They represent polynomials in
    `S = R[y]/(y^K + 1)`. This function returns
    `(L_1 + y^r L_2, L_1 - y^r L_2)`, as a list.
    """
    assert len(L1) == len(L2)
    assert 0 <= r <= len(L1)

    K = len(L1)
    return [L1[i] - L2[i + K - r] for i in range(r)] + \
        [L1[i] + L2[i - r] for i in range(r, K)], \
        [L1[i] + L2[i + K - r] for i in range(r)] + \
        [L1[i] - L2[i - r] for i in range(r, K)]


def _inverse_butterfly(L1, L2, r):
    r"""
    L1 and L2 are both lists of length `K`, and
    `0 \leq r \leq K`. They represent polynomials in
    `S = R[y]/(y^K + 1)`. This function returns
    `(L_1 + L_2, y^{-r}*(L_1 - L_2))`, as a list.
    """
    assert len(L1) == len(L2)
    assert 0 <= r <= len(L1)

    K = len(L1)
    return [L1[i] + L2[i] for i in range(K)], \
        [L1[i] - L2[i] for i in range(r, K)] + \
        [L2[i] - L1[i] for i in range(r)]


def _fft(L, K, start, depth, root):
    r"""
    L is a list of length `M = 2^m`, each entry a list of
    length `K = 2^k`.

    This function only operates on the [start, start + D) portion of L,
    where `D = 2^\text{depth}`. This portion is interpreted as
    a polynomial in `S[x]/(x^D - y^(2*root))`, where
    `S = R[y]/(y^K + 1)`.

    This function performs an inplace FFT, i.e. evaluates the
    polynomial at x = each `D`-th root of unity in S (namely the powers
    of `y^{2K/D}`), with results in bit-reversed order.
    """
    half = 1 << (depth - 1)
    start2 = start + half

    # reduce mod (x^(D/2) - y^root) and mod (x^(D/2) + y^root)
    for i in range(half):
        L[start + i], L[start2 + i] = \
            _forward_butterfly(L[start + i], L[start2 + i], root)

    # recurse into each half
    if depth >= 2:
        _fft(L, K, start, depth - 1, root >> 1)
        _fft(L, K, start2, depth - 1, (root + K) >> 1)


def _ifft(L, K, start, depth, root):
    r"""
    Inverse operation of ``_fft_trunc()`` (except that
    result is a factor of ``2^depth`` too big)
    """
    half = 1 << (depth - 1)
    start2 = start + half

    # recurse into each half
    if depth >= 2:
        _ifft(L, K, start, depth - 1, root >> 1)
        _ifft(L, K, start2, depth - 1, (root + K) >> 1)

    # CRT together (x^(D/2) - y^root) and mod (x^(D/2) + y^root)
    for i in range(half):
        L[start + i], L[start2 + i] = \
            _inverse_butterfly(L[start + i], L[start2 + i], root)

# -------------------------------------------------------------------
#      splitting and recombining routines


def _split(L, m, k):
    """
    Assumes L is a list of length `2^{m+k-1}`. Splits it into
    `2^m` lists of length `2^{k-1}`, returned as a list
    of lists. Each list is zero padded up to length `2^k`.

    EXAMPLES::

        sage: from sage.rings.polynomial.convolution import _split
        sage: _split([1, 2, 3, 4, 5, 6, 7, 8], 2, 2)
        [[1, 2, 0, 0], [3, 4, 0, 0], [5, 6, 0, 0], [7, 8, 0, 0]]
        sage: _split([1, 2, 3, 4, 5, 6, 7, 8], 1, 3)
        [[1, 2, 3, 4, 0, 0, 0, 0], [5, 6, 7, 8, 0, 0, 0, 0]]
        sage: _split([1, 2, 3, 4, 5, 6, 7, 8], 3, 1)
        [[1, 0], [2, 0], [3, 0], [4, 0], [5, 0], [6, 0], [7, 0], [8, 0]]
    """
    K = 1 << (k - 1)
    zero = parent(L[0])(0)
    zeroes = [zero] * K
    return [[L[i + j] for j in range(K)] + zeroes for i in range(0, K << m, K)]


def _combine(L, m, k):
    r"""
    Assumes L is a list of length `2^m`, each entry a list of
    length `2^k`. Combines together into a single list,
    effectively inverting ``_split()``, but overlaying
    coefficients, i.e. list #i gets added in starting at position
    `2^{k-1} i`. Note that the second half of the last list is
    ignored.
    """
    M = 1 << m
    half_K = 1 << (k - 1)
    return [L[0][j] for j in range(half_K)] + \
        [L[i + 1][j] + L[i][j + half_K]
         for i in range(M - 1) for j in range(half_K)]


def _nega_combine(L, m, k):
    r"""
    Same as ``_combine()``, but doesn't ignore the second
    half of the last list; instead it makes that piece wrap around
    negacyclically.
    """
    M = 1 << m
    half_K = 1 << (k - 1)
    return [L[0][j] - L[M - 1][j + half_K] for j in range(half_K)] + \
        [L[i + 1][j] + L[i][j + half_K]
         for i in range(M - 1) for j in range(half_K)]


# -------------------------------------------------------------------
#      FFT-based convolution and negacyclic convolutions


def _negaconvolution(L1, L2, n):
    """
    Negacyclic convolution of L1 and L2.

    L1 and L2 must both have length `2^n`.
    """
    if n <= 3:    # arbitrary cutoff
        return _negaconvolution_naive(L1, L2)
    return _negaconvolution_fft(L1, L2, n)


def _negaconvolution_fft(L1, L2, n):
    r"""
    Return negacyclic convolution of lists L1 and L2, using FFT algorithm.

    L1 and L2 must both have length `2^n`, where `n \geq 3`.
    Assumes all entries of L1 and L2 belong to the same ring.

    EXAMPLES::

        sage: from sage.rings.polynomial.convolution import _negaconvolution_naive
        sage: from sage.rings.polynomial.convolution import _negaconvolution_fft
        sage: _negaconvolution_naive(list(range(8)), list(range(5, 13)))
        [-224, -234, -224, -192, -136, -54, 56, 196]
        sage: _negaconvolution_fft(list(range(8)), list(range(5, 13)), 3)
        [-224, -234, -224, -192, -136, -54, 56, 196]

    TESTS::

        sage: for n in range(3, 10):
        ....:    L1 = [ZZ.random_element(100) for _ in range(1 << n)]
        ....:    L2 = [ZZ.random_element(100) for _ in range(1 << n)]
        ....:    assert _negaconvolution_naive(L1, L2) == _negaconvolution_fft(L1, L2, n)
    """
    assert n >= 3

    R = parent(L1[0])

    # split into 2^m pieces of 2^(k-1) coefficients each, with k as small
    # as possible, subject to m <= k (so that the ring of Fourier coefficients
    # has enough roots of unity)
    m = (n + 1) >> 1
    k = n + 1 - m

    M = 1 << m
    K = 1 << k

    # split inputs into polynomials
    L1 = _split(L1, m, k)
    L2 = _split(L2, m, k)

    # fft each input
    _fft(L1, K, 0, m, K >> 1)
    _fft(L2, K, 0, m, K >> 1)

    # pointwise multiply
    L3 = [_negaconvolution(L1[i], L2[i], k) for i in range(M)]

    # inverse fft
    _ifft(L3, K, 0, m, K >> 1)

    # combine back into a single list
    L3 = _nega_combine(L3, m, k)

    # normalise
    return [R(x / M) for x in L3]


def _convolution_fft(L1, L2):
    r"""
    Return convolution of non-empty lists L1 and L2, using FFT algorithm.

    L1 and L2 may have arbitrary lengths `\geq 4`.
    Assumes all entries of L1 and L2 belong to the same ring.

    EXAMPLES::

        sage: from sage.rings.polynomial.convolution import _convolution_naive
        sage: from sage.rings.polynomial.convolution import _convolution_fft
        sage: _convolution_naive([1, 2, 3], [4, 5, 6])
        [4, 13, 28, 27, 18]
        sage: _convolution_fft([1, 2, 3], [4, 5, 6])
        [4, 13, 28, 27, 18]

    TESTS::

        sage: for len1 in range(4, 30):
        ....:    for len2 in range(4, 30):
        ....:       L1 = [ZZ.random_element(100) for _ in range(len1)]
        ....:       L2 = [ZZ.random_element(100) for _ in range(len2)]
        ....:       assert _convolution_naive(L1, L2) == _convolution_fft(L1, L2)
    """
    R = parent(L1[0])

    # choose n so that output convolution length is 2^n
    len1 = len(L1)
    len2 = len(L2)
    outlen = len1 + len2 - 1
    n = int(ceil(log(outlen) / log(2.0)))

    # split into 2^m pieces of 2^(k-1) coefficients each, with k as small
    # as possible, subject to m <= k + 1 (so that the ring of Fourier
    # coefficients has enough roots of unity)
    m = (n >> 1) + 1
    k = n + 1 - m

    N = 1 << n
    M = 1 << m
    K = 1 << k

    # zero pad inputs up to length N
    zero = R(0)
    L1 = L1 + [zero] * (N - len(L1))
    L2 = L2 + [zero] * (N - len(L2))

    # split inputs into polynomials
    L1 = _split(L1, m, k)
    L2 = _split(L2, m, k)

    # fft each input
    _fft(L1, K, 0, m, K)
    _fft(L2, K, 0, m, K)

    # pointwise multiply
    L3 = [_negaconvolution(L1[i], L2[i], k) for i in range(M)]

    # inverse fft
    _ifft(L3, K, 0, m, K)

    # combine back into a single list
    L3 = _combine(L3, m, k)

    # normalise, and truncate to correct length
    return [R(L3[i] / M) for i in range(outlen)]
