# sage.doctest: needs sage.graphs sage.modules
r"""
Hadamard matrices

A Hadamard matrix is an `n\times n` matrix `H` whose entries are either `+1` or `-1`
and whose rows are mutually orthogonal. For example, the matrix `H_2`
defined by

.. MATH::

    \left(\begin{array}{rr}
    1 & 1 \\
    1 & -1
    \end{array}\right)

is a Hadamard matrix. An `n\times n` matrix `H` whose entries are either `+1` or
`-1` is a Hadamard matrix if and only if:

(a) `|det(H)|=n^{n/2}` or

(b)  `H*H^t = n\cdot I_n`, where `I_n` is the identity matrix.

In general, the tensor product of an `m\times m` Hadamard matrix and an
`n\times n` Hadamard matrix is an `(mn)\times (mn)` matrix. In
particular, if there is an `n\times n` Hadamard matrix then there is a
`(2n)\times (2n)` Hadamard matrix (since one may tensor with `H_2`).
This particular case is sometimes called the Sylvester construction.

The Hadamard conjecture (possibly due to Paley) states that a Hadamard
matrix of order `n` exists if and only if `n= 1, 2` or `n` is a multiple
of `4`.

The module below implements constructions of Hadamard and skew Hadamard matrices
for all known orders `\le 1200`, plus some more greater than `1200`. It also
allows you to pull a Hadamard matrix from the database at [SloaHada]_.

The following code will test that a construction for all known orders `\le 4k`
is implemented. The assertion above can be verified by setting ``k=300``
(note that it will take a long time to run)::

    sage: from sage.combinat.matrices.hadamard_matrix import (hadamard_matrix,
    ....:                               skew_hadamard_matrix, is_hadamard_matrix,
    ....:                               is_skew_hadamard_matrix)
    sage: k = 20
    sage: unknown_hadamard = [668, 716, 892, 1132]
    sage: unknown_skew_hadamard = [356, 404, 428, 476, 596, 612, 668, 708, 712, 716,
    ....:                         764, 772, 804, 808, 820, 836, 856, 892, 900, 916,
    ....:                         932, 940, 952, 980, 996, 1004, 1012, 1028, 1036,
    ....:                         1044, 1060, 1076, 1100, 1108, 1132, 1140, 1148,
    ....:                         1156, 1180, 1192, 1196]
    sage: for n in range(1, k+1):
    ....:   if 4*n not in unknown_hadamard:
    ....:       H = hadamard_matrix(4*n, check=False)
    ....:       assert is_hadamard_matrix(H)
    ....:   if 4*n not in unknown_skew_hadamard:
    ....:       H = skew_hadamard_matrix(4*n, check=False)
    ....:       assert is_skew_hadamard_matrix(H)

AUTHORS:

- David Joyner (2009-05-17): initial version
- Matteo Cati (2023-03-18): implemented more constructions for Hadamard and skew
  Hadamard matrices, to cover all known orders up to 1200.

REFERENCES:

- [SloaHada]_

- [HadaWiki]_

- [Hora]_

- [CP2023]_
"""

# *****************************************************************************
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from math import sqrt
from urllib.request import urlopen

from sage.arith.misc import divisors, is_prime_power, is_square, is_prime
from sage.combinat.designs.difference_family import (get_fixed_relative_difference_set,
                                                     relative_difference_set_from_homomorphism,
                                                     skew_supplementary_difference_set,
                                                     complementary_difference_sets)
from sage.combinat.t_sequences import T_sequences_smallcases
from sage.cpython.string import bytes_to_str
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.matrix.constructor import (block_matrix,
                                     block_diagonal_matrix,
                                     diagonal_matrix,
                                     identity_matrix as I,
                                     ones_matrix as J,
                                     matrix,
                                     matrix_method,
                                     zero_matrix)
from sage.misc.unknown import Unknown
from sage.modules.free_module_element import vector


def normalise_hadamard(H, skew=False):
    r"""
    Return the normalised Hadamard matrix corresponding to ``H``.

    The normalised Hadamard matrix corresponding to a Hadamard matrix `H` is a
    matrix whose every entry in the first row and column is +1.

    If ``skew`` is True, the matrix returned will be skew-normal: a skew Hadamard
    matrix with first row of all `+1`.

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import normalise_hadamard, is_hadamard_matrix, skew_hadamard_matrix
        sage: H = normalise_hadamard(hadamard_matrix(4))
        sage: H == hadamard_matrix(4)
        True
        sage: H = normalise_hadamard(skew_hadamard_matrix(20, skew_normalize=False), skew=True)
        sage: is_hadamard_matrix(H, skew=True, normalized=True)
        True

    If ``skew`` is ``True`` but the Hadamard matrix is not skew, the matrix returned
    will not be normalized::

        sage: H = normalise_hadamard(hadamard_matrix(92), skew=True)
        sage: is_hadamard_matrix(H, normalized=True)
        False
    """

    if skew:
        dd = diagonal_matrix(H[0])
        return dd*H*dd
    else:
        for i in range(H.ncols()):
            if H[0, i] < 0:
                H.rescale_col(i, -1)
        for i in range(H.nrows()):
            if H[i, 0] < 0:
                H.rescale_row(i, -1)
        return H


def hadamard_matrix_paleyI(n, normalize=True):
    r"""
    Implement the Paley type I construction.

    The Paley type I case corresponds to the case `p=n-1 \cong 3 \mod{4}` for a
    prime power `p` (see [Hora]_).

    INPUT:

    - ``n`` -- the matrix size
    - ``normalize`` -- boolean (default: ``True``); whether to normalize the result

    EXAMPLES:

    We note that this method by default returns a normalised Hadamard matrix ::

        sage: from sage.combinat.matrices.hadamard_matrix import hadamard_matrix_paleyI
        sage: hadamard_matrix_paleyI(4)
        [ 1  1  1  1]
        [ 1 -1  1 -1]
        [ 1 -1 -1  1]
        [ 1  1 -1 -1]

    Otherwise, it returns a skew Hadamard matrix `H`, i.e. `H=S+I`, with
    `S=-S^\top`  ::

        sage: M = hadamard_matrix_paleyI(4, normalize=False); M
        [ 1  1  1  1]
        [-1  1  1 -1]
        [-1 -1  1  1]
        [-1  1 -1  1]
        sage: S = M - identity_matrix(4); -S == S.T
        True

    TESTS::

        sage: from sage.combinat.matrices.hadamard_matrix import is_hadamard_matrix
        sage: test_cases = [x+1 for x in range(100) if is_prime_power(x) and x%4==3]
        sage: all(is_hadamard_matrix(hadamard_matrix_paleyI(n),normalized=True,verbose=True)
        ....:     for n in test_cases)
        True
        sage: all(is_hadamard_matrix(hadamard_matrix_paleyI(n,normalize=False),verbose=True)
        ....:     for n in test_cases)
        True
    """
    p = n - 1
    if not (is_prime_power(p) and (p % 4 == 3)):
        raise ValueError("The order %s is not covered by the Paley type I construction." % n)

    from sage.rings.finite_rings.finite_field_constructor import FiniteField
    K = FiniteField(p, 'x')
    K_list = list(K)
    K_list.insert(0, K.zero())
    H = matrix(ZZ, [[(1 if (x-y).is_square() else -1)
                     for x in K_list]
                    for y in K_list])
    for i in range(n):
        H[i, 0] = -1
        H[0, i] = 1
    if normalize:
        for i in range(n):
            H[i, i] = -1
        H = normalise_hadamard(H)
    return H


def symmetric_conference_matrix_paley(n):
    r"""
    Construct a symmetric conference matrix of order n.

    A conference matrix is an `n\times n` matrix `C` with 0s on the main diagonal
    and 1s and -1s elsewhere, satisfying `CC^\top=(n-1)I`. This construction assumes
    that `q = n-1` is a prime power, with `q \cong 1 \mod 4`. See [Hora]_ or [Lon2013]_.

    These matrices are used in :func:`hadamard_matrix_paleyII`.

    INPUT:

    - ``n`` -- integer; the order of the symmetric conference matrix to construct

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import symmetric_conference_matrix_paley
        sage: symmetric_conference_matrix_paley(6)
        [ 0  1  1  1  1  1]
        [ 1  0  1 -1 -1  1]
        [ 1  1  0  1 -1 -1]
        [ 1 -1  1  0  1 -1]
        [ 1 -1 -1  1  0  1]
        [ 1  1 -1 -1  1  0]

    TESTS::

        sage: symmetric_conference_matrix_paley(5)
        Traceback (most recent call last):
        ...
        ValueError: The order 5 is not covered by Paley construction of symmetric conference matrices.
    """
    q = n - 1
    if not (is_prime_power(q) and (q % 4 == 1)):
        raise ValueError("The order %s is not covered by Paley construction of symmetric conference matrices." % n)

    from sage.rings.finite_rings.finite_field_constructor import FiniteField
    K = FiniteField(q, 'x')
    K_list = list(K)
    K_list.insert(0, K.zero())
    H = matrix(ZZ, [[(1 if (x-y).is_square() else -1)
                     for x in K_list]
                    for y in K_list])
    for i in range(n):
        H[0, i] = 1
        H[i, 0] = 1
        H[i, i] = 0
    return H


def hadamard_matrix_paleyII(n):
    r"""
    Implement the Paley type II construction.

    The Paley type II case corresponds to the case `p=n/2-1 \cong 1 \mod{4}` for a
    prime power `p` (see [Hora]_).

    EXAMPLES::

        sage: sage.combinat.matrices.hadamard_matrix.hadamard_matrix_paleyII(12).det()
        2985984
        sage: 12^6
        2985984

    We note that the method returns a normalised Hadamard matrix ::

        sage: sage.combinat.matrices.hadamard_matrix.hadamard_matrix_paleyII(12)
        [ 1  1| 1  1| 1  1| 1  1| 1  1| 1  1]
        [ 1 -1|-1  1|-1  1|-1  1|-1  1|-1  1]
        [-----+-----+-----+-----+-----+-----]
        [ 1 -1| 1 -1| 1  1|-1 -1|-1 -1| 1  1]
        [ 1  1|-1 -1| 1 -1|-1  1|-1  1| 1 -1]
        [-----+-----+-----+-----+-----+-----]
        [ 1 -1| 1  1| 1 -1| 1  1|-1 -1|-1 -1]
        [ 1  1| 1 -1|-1 -1| 1 -1|-1  1|-1  1]
        [-----+-----+-----+-----+-----+-----]
        [ 1 -1|-1 -1| 1  1| 1 -1| 1  1|-1 -1]
        [ 1  1|-1  1| 1 -1|-1 -1| 1 -1|-1  1]
        [-----+-----+-----+-----+-----+-----]
        [ 1 -1|-1 -1|-1 -1| 1  1| 1 -1| 1  1]
        [ 1  1|-1  1|-1  1| 1 -1|-1 -1| 1 -1]
        [-----+-----+-----+-----+-----+-----]
        [ 1 -1| 1  1|-1 -1|-1 -1| 1  1| 1 -1]
        [ 1  1| 1 -1|-1  1|-1  1| 1 -1|-1 -1]

    TESTS::

        sage: from sage.combinat.matrices.hadamard_matrix import (hadamard_matrix_paleyII, is_hadamard_matrix)
        sage: test_cases = [2*(x+1) for x in range(50) if is_prime_power(x) and x%4==1]
        sage: all(is_hadamard_matrix(hadamard_matrix_paleyII(n),normalized=True,verbose=True)
        ....:     for n in test_cases)
        True
    """
    q = n//2 - 1
    if not (n % 2 == 0 and is_prime_power(q) and (q % 4 == 1)):
        raise ValueError("The order %s is not covered by the Paley type II construction." % n)

    H = symmetric_conference_matrix_paley(q+1)

    tr = { 0: matrix(2, 2, [ 1, -1, -1, -1]),
           1: matrix(2, 2, [ 1,  1,  1, -1]),
          -1: matrix(2, 2, [-1, -1, -1,  1])}

    H = block_matrix(q+1, q+1, [tr[v] for r in H for v in r])

    return normalise_hadamard(H)


def hadamard_matrix_from_symmetric_conference_matrix(n, existence=False, check=True):
    r"""
    Construct a Hadamard matrix of order `n` from a symmetric conference matrix
    of order `n/2`.

    The construction is described in Theorem 4.3.24 of [IS2006]_.
    The symmetric conference matrices are obtained from
    :func:`sage.combinat.matrices.hadamard_matrix.symmetric_conference_matrix`.

    INPUT:

    - ``n`` -- integer; the order of the matrix to be constructed
    - ``existence`` -- boolean (default: ``False``); if ``True``, only check if
      the matrix exists
    - ``check`` -- boolean (default: ``True``); if ``True``, check that the matrix
      is a Hadamard before returning

    OUTPUT:

    If ``existence=False``, returns the Hadamard matrix of order `n`. It raises
    an error if no data is available to construct the matrix of the given order,
    or if `n` does not satisfies the constraints.
    If ``existence=True``, returns a boolean representing whether the matrix
    can be constructed or not.

    EXAMPLES:

    By default the function returns the Hadamard matrix ::

        sage: from sage.combinat.matrices.hadamard_matrix import hadamard_matrix_from_symmetric_conference_matrix
        sage: hadamard_matrix_from_symmetric_conference_matrix(20)
        20 x 20 dense matrix over Integer Ring...

    If ``existence`` is set to True, the function returns True if the matrix exists,
    False if the conference matrix does not exist, and Unknown if the conference
    matrix cannot be constructed yet ::

        sage: hadamard_matrix_from_symmetric_conference_matrix(12, existence=True)
        True
        sage: hadamard_matrix_from_symmetric_conference_matrix(4*787, existence=True)
        True

    TESTS::

        sage: from sage.combinat.matrices.hadamard_matrix import is_hadamard_matrix
        sage: is_hadamard_matrix(hadamard_matrix_from_symmetric_conference_matrix(60, check=False))
        True
        sage: hadamard_matrix_from_symmetric_conference_matrix(64, existence=True)
        False
        sage: hadamard_matrix_from_symmetric_conference_matrix(4*787, existence=True)
        True
        sage: hadamard_matrix_from_symmetric_conference_matrix(64)
        Traceback (most recent call last):
        ...
        ValueError: Cannot construct Hadamard matrix of order 64, a symmetric conference matrix of order 32 is not available in sage.
        sage: hadamard_matrix_from_symmetric_conference_matrix(14)
        Traceback (most recent call last):
        ...
        ValueError: No Hadamard matrix of order 14 exists.
    """
    if n < 0 or n % 4 != 0:
        raise ValueError(f'No Hadamard matrix of order {n} exists.')

    m = n//2
    exists = symmetric_conference_matrix(m, existence=True)

    if existence:
        return exists

    if not exists:
        raise ValueError(f'Cannot construct Hadamard matrix of order {n}, a symmetric conference matrix of order {m} is not available in sage.')

    C = symmetric_conference_matrix(m)

    H = block_matrix([[C + I(m), C - I(m)],
                      [C - I(m), -C - I(m)]])

    if check:
        assert is_hadamard_matrix(H)
    return H


def hadamard_matrix_miyamoto_construction(n, existence=False, check=True):
    r"""
    Construct Hadamard matrix using the Miyamoto construction.

    If `q = n/4` is a prime power, and there exists a Hadamard matrix of order
    `q-1`, then a Hadamard matrix of order `n` can be constructed (see [Miy1991]_).

    INPUT:

    - ``n`` -- integer; the order of the matrix to be constructed
    - ``check`` -- boolean (default: ``True``); if ``True``, check that the matrix
      is a Hadamard before returning
    - ``existence`` -- boolean (default: ``False``); if ``True``, only check if
      the matrix exists

    OUTPUT:

    If ``existence=False``, returns the Hadamard matrix of order `n`. It raises
    an error if no data is available to construct the matrix of the given order,
    or if `n` does not satisfies the constraints.
    If ``existence=True``, returns a boolean representing whether the matrix
    can be constructed or not.

    EXAMPLES:

    By default the function returns the Hadamard matrix ::

        sage: from sage.combinat.matrices.hadamard_matrix import hadamard_matrix_miyamoto_construction
        sage: hadamard_matrix_miyamoto_construction(20)
        20 x 20 dense matrix over Integer Ring...

    If ``existence`` is set to True, the function returns a boolean ::

        sage: hadamard_matrix_miyamoto_construction(36, existence=True)
        True

    TESTS::

        sage: from sage.combinat.matrices.hadamard_matrix import is_hadamard_matrix
        sage: is_hadamard_matrix(hadamard_matrix_miyamoto_construction(68, check=False))
        True
        sage: hadamard_matrix_miyamoto_construction(64, existence=True)
        False
        sage: hadamard_matrix_miyamoto_construction(4*65, existence=True)
        True
        sage: is_hadamard_matrix(hadamard_matrix_miyamoto_construction(4*65, check=False))
        True
        sage: hadamard_matrix_miyamoto_construction(64)
        Traceback (most recent call last):
        ...
        ValueError: The order 64 is not covered by Miyamoto construction.
        sage: hadamard_matrix_miyamoto_construction(14)
        Traceback (most recent call last):
        ...
        ValueError: No Hadamard matrix of order 14 exists.
    """
    if n < 0 or n % 4 != 0:
        raise ValueError(f'No Hadamard matrix of order {n} exists.')

    q = n // 4
    if existence:
        # return is_prime_power(q) and q % 4 == 1 and hadamard_matrix(q-1, existence=True) is True
        return symmetric_conference_matrix(q+1, existence=True) and hadamard_matrix(q-1, existence=True) is True

    # if not (is_prime_power(q) and q % 4 == 1 and hadamard_matrix(q-1, existence=True)):
    if not (symmetric_conference_matrix(q+1, existence=True) and hadamard_matrix(q-1, existence=True)):
        raise ValueError(f'The order {n} is not covered by Miyamoto construction.')

    m = (q-1) // 2

    C = symmetric_conference_matrix(q + 1)

    neg = [i for i in range(2, m+2) if C[1, i] == -1]
    pos = [i for i in range(m+2, 2*m+2) if C[1, i] == 1]

    for i, j in zip(neg, pos):
        C.swap_rows(i, j)
        C.swap_columns(i, j)

    C1 = -C.submatrix(row=2, col=2, nrows=m, ncols=m)
    C2 = C.submatrix(row=2, col=m+2, nrows=m, ncols=m)
    C4 = C.submatrix(row=m+2, col=m+2, nrows=m, ncols=m)

    K = hadamard_matrix(q - 1)
    K1 = K.submatrix(row=0, col=0, nrows=(q-1)//2, ncols=(q-1)//2)
    K2 = K.submatrix(row=0, col=(q-1)//2, nrows=(q-1)//2, ncols=(q-1)//2)
    K3 = -K.submatrix(row=(q-1)//2, col=0, nrows=(q-1)//2, ncols=(q-1)//2)
    K4 = K.submatrix(row=(q-1)//2, col=(q-1)//2, nrows=(q-1)//2, ncols=(q-1)//2)

    Zr = zero_matrix(m)
    Us = [[C1, C2, Zr, Zr], [C2.T, C4, Zr, Zr], [Zr, Zr, C1, C2], [Zr, Zr, C2.T, C4]]
    Vs = [[I(m), Zr, K1, K2], [Zr, I(m), K3, K4], [K1.T, K3.T, I(m), Zr], [K2.T, K4.T, Zr, I(m)]]

    def T(i, j):
        return block_matrix([[Us[i][j]+Vs[i][j], Us[i][j]-Vs[i][j]],
                             [Us[i][j]-Vs[i][j], Us[i][j]+Vs[i][j]]])

    e = matrix([[1] * (2*m)])
    one = matrix([1])
    H = block_matrix([[ one,       -e,  one,        e,  one,        e,  one,        e],
                      [-e.T,  T(0, 0),  e.T,  T(0, 1),  e.T,  T(0, 2),  e.T,  T(0, 3)],
                      [-one,       -e,  one,       -e,  one,        e, -one,       -e],
                      [-e.T, -T(1, 0), -e.T,  T(1, 1),  e.T,  T(1, 2), -e.T, -T(1, 3)],
                      [-one,       -e, -one,       -e,  one,       -e,  one,        e],
                      [-e.T, -T(2, 0), -e.T, -T(2, 1), -e.T,  T(2, 2),  e.T,  T(2, 3)],
                      [-one,       -e,  one,        e, -one,       -e,  one,       -e],
                      [-e.T, -T(3, 0),  e.T,  T(3, 1), -e.T, -T(3, 2), -e.T,  T(3, 3)]])

    if check:
        assert is_hadamard_matrix(H)
    return H


def hadamard_matrix_williamson_type(a, b, c, d, check=True):
    r"""
    Construction of Williamson type Hadamard matrix.

    Given `n\times n` circulant matrices `A`, `B`, `C`, `D` with 1,-1 entries,
    and satisfying `AA^\top + BB^\top + CC^\top + DD^\top = 4nI`,
    one can construct a  Hadamard matrix of order `4n`, cf. [Ha83]_.

    INPUT:

    - ``a`` -- (1,-1) list; the 1st row of `A`
    - ``b`` -- (1,-1) list; the 1st row of `B`
    - ``d`` -- (1,-1) list; the 1st row of `C`
    - ``c`` -- (1,-1) list; the 1st row of `D`
    - ``check`` -- boolean (default: ``True``); whether to check that the output
      is a Hadamard matrix before returning it

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import hadamard_matrix_williamson_type
        sage: a = [ 1,  1, 1]
        sage: b = [ 1, -1, -1]
        sage: c = [ 1, -1, -1]
        sage: d = [ 1, -1, -1]
        sage: M = hadamard_matrix_williamson_type(a,b,c,d,check=True)

    TESTS::

        sage: from sage.combinat.matrices.hadamard_matrix import hadamard_matrix_williamson_type, is_hadamard_matrix
        sage: a = [ 1,  1, 1]
        sage: b = [ 1, -1, -1]
        sage: c = [ 1, -1, -1]
        sage: d = [ 1, -1, -1]
        sage: is_hadamard_matrix(hadamard_matrix_williamson_type(a,b,c,d))
        True
        sage: e = [1, 1, 1]
        sage: hadamard_matrix_williamson_type(a,b,c,e, check=True)
        Traceback (most recent call last):
        ...
        AssertionError
        sage: f = [1, -1, 1, -1]
        sage: hadamard_matrix_williamson_type(a,b,c,f, check=True)
        Traceback (most recent call last):
        ...
        AssertionError
    """
    A, B, C, D = map(matrix.circulant, [a, b, c, d])

    n = len(a)
    assert len(a) == len(b) == len(c) == len(d)
    assert A*A.T+B*B.T+C*C.T+D*D.T == 4*n*I(n)

    M = block_matrix([[ A,  B,  C,  D],
                      [-B,  A, -D,  C],
                      [-C,  D,  A, -B],
                      [-D, -C,  B,  A]])
    if check:
        assert is_hadamard_matrix(M, normalized=False, skew=False)
    return M


def williamson_type_quadruples_smallcases(n, existence=False):
    r"""
    Quadruples of matrices that can be used to construct Williamson type Hadamard matrices.

    This function contains for some values of n, four `n\times n` matrices used in the
    Williamson construction of Hadamard matrices. Namely, the function returns the first row of
    4 `n\times n` circulant matrices with the properties described in
    :func:`sage.combinat.matrices.hadamard_matrix.hadamard_matrix_williamson_type`.
    The matrices for `n = 3, 5, ..., 29, 37, 43` are given in [Ha83]_. The matrices
    for `n = 31, 33, 39, 41, 45, 49, 51, 55, 57, 61, 63` are given in [Lon2013]_.

    INPUT:

    - ``n`` -- integer; the order of the matrices to be returned
    - ``existence`` -- boolean (dafault: ``False``); if ``True``, only check that
      we have the quadruple

    OUTPUT:

    If ``existence`` is false, returns a tuple containing four vectors, each being the first line
    of one of the four matrices. It raises an error if no such matrices are available.
    If ``existence`` is true, returns a boolean representing whether the matrices are available or not.

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import williamson_type_quadruples_smallcases
        sage: williamson_type_quadruples_smallcases(29)
        ((1, 1, 1, -1, -1, -1, 1, 1, -1, -1, 1, -1, 1, -1, -1, -1, -1, 1, -1, 1, -1, -1, 1, 1, -1, -1, -1, 1, 1),
         (1, -1, 1, -1, -1, -1, 1, 1, -1, -1, 1, -1, 1, 1, 1, 1, 1, 1, -1, 1, -1, -1, 1, 1, -1, -1, -1, 1, -1),
         (1, 1, 1, 1, -1, 1, 1, -1, 1, -1, -1, -1, 1, 1, 1, 1, 1, 1, -1, -1, -1, 1, -1, 1, 1, -1, 1, 1, 1),
         (1, 1, -1, -1, 1, -1, -1, 1, -1, 1, 1, 1, -1, 1, 1, 1, 1, -1, 1, 1, 1, -1, 1, -1, -1, 1, -1, -1, 1))
        sage: williamson_type_quadruples_smallcases(43, existence=True)
        True

    TESTS::

        sage: williamson_type_quadruples_smallcases(123, existence=True)
        False
        sage: williamson_type_quadruples_smallcases(123)
        Traceback (most recent call last):
        ...
        ValueError: The Williamson type quadruple of order 123 is not yet implemented.
    """
    db = {
        1: ('+', '+', '+', '+'),
        3: ('+++', '+--', '+--', '+--'),
        5: ('+-++-', '++--+', '+----', '+----'),
        7: ('+--++--', '+-+--+-', '++----+', '+------'),
        9: ('+---++---', '+--+--+--', '+-+----+-', '++------+'),
        11: ('++--------+', '++-+-++-+-+', '++-++--++-+', '+-++----++-'),
        13: ('++++-+--+-+++', '+---+-++-+---', '++---+--+---+', '++---+--+---+'),
        15: ('+-+---++++---+-', '++-++------++-+',
             '++-++++--++++-+', '++-++-+--+-++-+'),
        17: ('+---+++----+++---', '++-+---+--+---+-+',
             '+--+-++++++++-+--', '+-++-+++--+++-++-'),
        19: ('++--+++-+--+-+++--+', '++-++--+-++-+--++-+',
             '+-+---++++++++---+-', '++--+-++++++++-+--+'),
        21: ('+--++++---++---++++--', '++++-+---+--+---+-+++',
             '++--+-+-++--++-+-+--+', '++-+++++-+--+-+++++-+'),
        23: ('++---+---+-++-+---+---+', '+-++-++--++++++--++-++-',
             '+++---++-+-++-+-++---++', '+++-+++-+------+-+++-++'),
        25: ('++++-+-+-+--++--+-+-+-+++', '++--+--+-++++++++-+--+--+',
             '+++--+--++++--++++--+--++', '+-+--+++--++++++--+++--+-'),
        27: ('+--+--+-+++--++--+++-+--+--', '+++-++-+---++--++---+-++-++',
             '+---+++++-+-++++-+-+++++---', '+---+++++-+-++++-+-+++++---'),
        29: ('+++---++--+-+----+-+--++---++', '+-+---++--+-++++++-+--++---+-',
             '++++-++-+---++++++---+-++-+++', '++--+--+-+++-++++-+++-+--+--+'),
        31: ('++++++-+--+---++++---+--+-+++++', '+--++---+-+-++----++-+-+---++--',
             '+--++---+-+-++----++-+-+---++--', '+-----+-++-+++----+++-++-+-----'),
        33: ('++++++-+-+-+++------+++-+-+-+++++', '++-+-++-+----+++--+++----+-++-+-+',
             '++--++-+++-+--+-++-+--+-+++-++--+', '+--++--+++++-++----++-+++++--++--'),
        37: ('+--+-+-+-++---+--++++--+---++-+-+-+--', '+---++-++--+-+-++----++-+-+--++-++---',
             '+++++-+-----++----++----++-----+-++++', '+--+++-+-----+----++----+-----+-+++--'),
        39: ('+++--+-+-----+--++----++--+-----+-+--++', '+++--++-+---+-+--+----+--+-+---+-++--++',
             '++++---+--++----+-+--+-+----++--+---+++', '+---++-+-+-----+++-++-+++-----+-+-++---'),
        41: ('++++--+-++++-++--++----++--++-++++-+--+++', '++++--+-++++-++--++----++--++-++++-+--+++',
             '+++-++-+-+-+-----+++--+++-----+-+-+-++-++', '+--+--+-+-+-+++++---++---+++++-+-+-+--+--'),
        43: ('++---++++-+--+--++--------++--+--+-++++---+', '+++-+-++--+-+-++++-+----+-++++-+-+--++-+-++',
             '++-++++++----+-+--++-++-++--+-+----++++++-+', '+---++--++++-+-+++-++--++-+++-+-++++--++---'),
        45: ('+++++-++----+-++--++-++-++--++-+----++-++++', '+++---++--+-+-+-++--------++-+-+-+--++---++',
             '++-+-++++-+--+--+++--++--+++--+--+-++++-+-+', '+-++-----++++-+-+++-++++-+++-+-++++-----++-'),
        49: ('++++-++-+---++-+++---++-++-++---+++-++---+-++-+++', '++++-++-+---++-+++---++-++-++---+++-++---+-++-+++',
             '+----+-++++--+-+++-+-+++--+++-+-+++-+--++++-+----', '+++++-+----++-+---+-+---++---+-+---+-++----+-++++'),
        51: ('+---+++-++-+-+++--+++++--++--+++++--+++-+-++-+++---', '----+++-++-+-+++--+++++--++--+++++--+++-+-++-+++---',
             '-+--+----+-+++-+-+++++--+--+--+++++-+-+++-+----+--+', '-+--+----+-+++-+-+++++--+--+--+++++-+-+++-+----+--+'),
        55: ('+-+--+-+-++--+-+++++-+++--++++--+++-+++++-+--++-+-+--+-', '--+--+-+-++--+-+++++-+++--++++--+++-+++++-+--++-+-+--+-',
             '+++----++-++--++----+-+-++++++++-+-+----++--++-++----++', '+++----++-++--++----+-+-++++++++-+-+----++--++-++----++'),
        57: ('+---++-+--++++-+++-++---+-++++++-+---++-+++-++++--+-++---', '----++-+--++++-+++-++---+-++++++-+---++-+++-++++--+-++---',
             '--+-+-+++--+--+-++---+++++-++++-+++++---++-+--+--+++-+-+-', '--+-+-+++--+--+-++---+++++-++++-+++++---++-+--+--+++-+-+-'),
        61: ('++--+--++--+-+-++++--+-----+------+-----+--++++-+-+--++--+--+', '++--+--++--+-+-++++--+-----+------+-----+--++++-+-+--++--+--+',
             '+---+-+-++++---++--+-++-+---++++++---+-++-+--++---++++-+-+---', '++++-+-+----+++--++-+--+-+++------+++-+--+-++--+++----+-+-+++'),
        63: ('++-+++--++-++--+--+-++-+-+++--------+++-+-++-+--+--++-++--+++-+', '-+-+++--++-++--+--+-++-+-+++--------+++-+-++-+--+--++-++--+++-+',
             '++++-++-+-++++-+---+---+++---++++++---+++---+---+-++++-+-++-+++', '++++-++-+-++++-+---+---+++---++++++---+++---+---+-++++-+-++-+++'),
    }

    def pmtoZ(s):
        return [1 if x == '+' else -1 for x in s]

    if existence:
        return n in db

    if n not in db:
        raise ValueError("The Williamson type quadruple of order %s is not yet implemented." % n)

    a, b, c, d = (vector(pmtoZ(s)) for s in db[n])
    return a, b, c, d


def williamson_hadamard_matrix_smallcases(n, existence=False, check=True):
    r"""
    Construct Williamson type Hadamard matrices for some small values of n.

    This function uses the data contained in
    :func:`sage.combinat.matrices.hadamard_matrix.williamson_type_quadruples_smallcases`
    to create Hadamard matrices of the Williamson type, using the construction from
    :func:`sage.combinat.matrices.hadamard_matrix.hadamard_matrix_williamson_type`.

    INPUT:

    - ``n`` -- integer; the order of the matrix
    - ``existence`` -- boolean (dafault: ``False``); if ``True``, only check that
      we can do the construction
    - ``check`` -- boolean (default: ``True``); if ``True`` check the result

    TESTS::

        sage: from sage.combinat.matrices.hadamard_matrix import williamson_hadamard_matrix_smallcases
        sage: williamson_hadamard_matrix_smallcases(116)
        116 x 116 dense matrix over Integer Ring...
        sage: williamson_hadamard_matrix_smallcases(172)
        172 x 172 dense matrix over Integer Ring...
        sage: williamson_hadamard_matrix_smallcases(1000)
        Traceback (most recent call last):
        ...
        ValueError: The Williamson type Hadamard matrix of order 1000 is not yet implemented.
    """
    assert n % 4 == 0

    if not williamson_type_quadruples_smallcases(n // 4, existence=True):
        if existence:
            return False
        raise ValueError("The Williamson type Hadamard matrix of order %s is not yet implemented." % n)

    if existence:
        return True

    a, b, c, d = williamson_type_quadruples_smallcases(n//4)
    return hadamard_matrix_williamson_type(a, b, c, d, check=check)


def hadamard_matrix_156():
    r"""
    Construct a Hadamard matrix of order 156.

    The matrix is created using the construction detailed in [BH1965]_.
    This uses four circulant matrices of size `13\times 13`,
    which are composed into a `156\times 156` block matrix.

    TESTS::

        sage: from sage.combinat.matrices.hadamard_matrix import is_hadamard_matrix, hadamard_matrix_156
        sage: is_hadamard_matrix(hadamard_matrix_156())
        True
        sage: hadamard_matrix_156()
        156 x 156 dense matrix over Integer Ring...
    """
    a = [1, 1, -1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1]
    b = [1, -1, -1, -1, 1, 1, 1, 1, 1, 1, -1, -1, -1]
    c = [1, 1, 1, -1, 1, 1, -1, -1, 1, 1, -1, 1, 1]
    d = [1, 1, -1, 1, -1, 1, 1, 1, 1, -1, 1, -1, 1]

    A, B, C, D = map(matrix.circulant, [a, b, c, d])

    return block_matrix([[ A,  A,  A,  B, -B,  C, -C, -D,  B,  C, -D, -D],
                         [ A, -A,  B, -A, -B, -D,  D, -C, -B, -D, -C, -C],
                         [ A, -B, -A,  A, -D,  D, -B,  B, -C, -D,  C, -C],
                         [ B,  A, -A, -A,  D,  D,  D,  C,  C, -B, -B, -C],
                         [ B, -D,  D,  D,  A,  A,  A,  C, -C,  B, -C,  B],
                         [ B,  C, -D,  D,  A, -A,  C, -A, -D,  C,  B, -B],
                         [ D, -C,  B, -B,  A, -C, -A,  A,  B,  C,  D, -D],
                         [-C, -D, -C, -D,  C,  A, -A, -A, -D,  B, -B, -B],
                         [ D, -C, -B, -B, -B,  C,  C, -D,  A,  A,  A,  D],
                         [-D, -B,  C,  C,  C,  B,  B, -D,  A, -A,  D, -A],
                         [ C, -B, -C,  C,  D, -B, -D, -B,  A, -D, -A,  A],
                         [-C, -D, -D,  C, -C, -B,  B,  B,  D,  A, -A, -A]])


def construction_four_symbol_delta_code_I(X, Y, Z, W):
    r"""
    Construct 4-symbol `\delta` code of length `2n+1`.

    The 4-symbol `\delta` code is constructed from sequences `X, Y, Z, W` of
    length `n+1`, `n+1`, `n`, `n` satisfying for all `s > 0`:

    .. MATH::

        N_X(s) + N_Y(s) + N_Z(s) + N_W(s) = 0

    where `N_A(s)` is the nonperiodic correlation function:

    .. MATH::

        N_A(s) = \sum_{i=1}^{n-s}a_ia_{i+s}

    The construction (detailed in [Tur1974]_) is as follows:

    .. MATH::

        \begin{aligned}
        T_1 &= X;Z \\
        T_2 &= X;-Z \\
        T_3 &= Y;W \\
        T_4 &= Y;-W
        \end{aligned}

    INPUT:

    - ``X`` -- list; the first sequence (length `n+1`)
    - ``Y`` -- list; the second sequence (length `n+1`)
    - ``Z`` -- list; the third sequence (length `n`)
    - ``W`` -- list; the fourth sequence (length `n`)

    OUTPUT: tuple containing the 4-symbol `\delta` code of length `2n+1`

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import construction_four_symbol_delta_code_I
        sage: construction_four_symbol_delta_code_I([1, 1], [1, -1], [1], [1])
        ([1, 1, 1], [1, 1, -1], [1, -1, 1], [1, -1, -1])

    TESTS::

        sage: construction_four_symbol_delta_code_I([1, 1], [1, -1], [1, 1], [1])
        Traceback (most recent call last):
        ...
        AssertionError
        sage: construction_four_symbol_delta_code_I([1, 1], [1, 1], [-1], [1])
        Traceback (most recent call last):
        ...
        AssertionError
    """
    n = len(X)
    assert len(Y) == n and len(Z) == n-1 and len(W) == n-1

    def autocorrelation(seq, j):
        return sum([seq[i]*seq[i+j] for i in range(len(seq)-j)])
    for j in range(1, n):
        assert sum(autocorrelation(seq, j) for seq in [X, Y, Z, W]) == 0

    T1 = X + Z
    T2 = X + [-z for z in Z]
    T3 = Y + W
    T4 = Y + [-w for w in W]
    return T1, T2, T3, T4


def construction_four_symbol_delta_code_II(X, Y, Z, W):
    r"""
    Construct 4-symbol `\delta` code of length `4n+3`.

    The 4-symbol `\delta` code is constructed from sequences  `X, Y, Z, W` of
    length `n+1`, `n+1`, `n`, `n` satisfying for all `s > 0`:

    .. MATH::
        N_X(s) + N_Y(s) + N_Z(s) + N_W(s) = 0

    where `N_A(s)` is the nonperiodic correlation function:

    .. MATH::

        N_A(s) = \sum_{i=1}^{n-s}a_ia_{i+s}

    The construction (detailed in [Tur1974]_) is as follows (writing
    `A/B` to mean `A` alternated with `B`):

    .. MATH::

        \begin{aligned}
        T_1 &= X/Z;Y/W;1 \\
        T_2 &= X/Z;Y/-W;-1 \\
        T_3 &= X/Z;-Y/-W;1 \\
        T_4 &= X/Z;-Y/W;-1
        \end{aligned}

    INPUT:

    - ``X`` -- list; the first sequence (length `n+1`)
    - ``Y`` -- list; the second sequence (length `n+1`)
    - ``Z`` -- list; the third sequence (length `n`)
    - ``W`` -- list; the fourth sequence (length `n`)

    OUTPUT: tuple containing the four 4-symbol `\delta` code of length `4n+3`

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import construction_four_symbol_delta_code_II
        sage: construction_four_symbol_delta_code_II([1, 1], [1, -1], [1], [1])
        ([1, 1, 1, 1, 1, -1, 1],
         [1, 1, 1, 1, -1, -1, -1],
         [1, 1, 1, -1, -1, 1, 1],
         [1, 1, 1, -1, 1, 1, -1])

    TESTS::

        sage: construction_four_symbol_delta_code_II([1, 1], [1, -1], [1, 1], [1])
        Traceback (most recent call last):
        ...
        AssertionError
        sage: construction_four_symbol_delta_code_II([1, 1], [1, 1], [-1], [1, 1])
        Traceback (most recent call last):
        ...
        AssertionError
    """

    n = len(Z)
    assert len(X) == n+1 and len(Y) == n+1 and len(W) == n

    def autocorrelation(seq, j):
        return sum([seq[i]*seq[i+j] for i in range(len(seq)-j)])

    for j in range(1, n):
        assert sum(autocorrelation(seq, j) for seq in [X, Y, Z, W]) == 0

    def alternate(seq1, seq2):
        return [seq1[i//2] if i % 2 == 0 else seq2[(i-1)//2] for i in range(len(seq1)+len(seq2))]

    XaltZ = alternate(X, Z)
    Wneg = [-w for w in W]
    Yneg = [-y for y in Y]

    T1 = XaltZ + alternate(Y, W) + [1]
    T2 = XaltZ + alternate(Y, Wneg) + [-1]
    T3 = XaltZ + alternate(Yneg, Wneg) + [1]
    T4 = XaltZ + alternate(Yneg, W) + [-1]
    return T1, T2, T3, T4


def four_symbol_delta_code_smallcases(n, existence=False):
    r"""
    Return the 4-symobl `\delta` code of length ``n`` if available.

    The 4-symbol `\delta` codes are constructed using :func:`construction_four_symbol_delta_code_I`
    or :func:`construction_four_symbol_delta_code_II`.
    The base sequences used are taken from [Tur1974]_.

    INPUT:

    - ``n`` -- integer; the length of the desired 4-symbol `\delta` code
    - ``existence`` -- boolean (default: ``False``); if ``True``, only check if
      the sequences are available

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import four_symbol_delta_code_smallcases
        sage: four_symbol_delta_code_smallcases(3)
        ([1, -1, 1], [1, -1, -1], [1, 1, 1], [1, 1, -1])
        sage: four_symbol_delta_code_smallcases(3, existence=True)
        True

    TESTS::

        sage: four_symbol_delta_code_smallcases(17)
        Traceback (most recent call last):
        ...
        ValueError: The four-symbol delta code of length 17 have not yet been implemented
        sage: four_symbol_delta_code_smallcases(17, existence=True)
        False
    """
    db = {
        1: ([1, -1], [1, 1], [1], [1]),
        14: ([1, 1, -1, 1, 1, 1, -1, 1, -1, 1, 1, 1, -1, 1, 1],
             [1, 1, 1, -1, 1, 1, -1, -1, -1, 1, 1, -1, 1, 1, -1],
             [1, 1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1, -1, -1],
             [1, -1, -1, -1, -1, 1, -1, 1, -1, 1, 1, 1, 1, -1])
    }

    T1, T2, T3, T4 = None, None, None, None
    if n % 2 == 1 and (n-1)//2 in db:
        if existence:
            return True
        X, Y, Z, W = db[(n-1)//2]
        T1, T2, T3, T4 = construction_four_symbol_delta_code_I(X, Y, Z, W)
    elif n % 4 == 3 and (n-3) // 4 in db:
        if existence:
            return True
        X, Y, Z, W = db[(n-3)//4]
        T1, T2, T3, T4 = construction_four_symbol_delta_code_II(X, Y, Z, W)

    if existence:
        return False

    if T1 is None:
        raise ValueError("The four-symbol delta code of length %s have not yet been implemented" % n)

    return T1, T2, T3, T4


def _construction_goethals_seidel_matrix(A, B, C, D):
    r"""
    Construct the Goethals Seidel matrix.

    The matrix is described in [GS70s]_. Given matrices `A`, `B`, `C`, `D`
    the construction is:

    .. MATH::

        \left(\begin{array}{rrrr}
        A & BR & CR & DR \\
        -BR & A & -D\top R & C\top R \\
        -CR & D\top R & A & -B\top R \\
        -DR & -C\top R & B\top R & A
        \end{array}\right)

    Where `R` is the anti-diagonal matrix with all nonzero entries
    equal to one.

    INPUT:

    - ``A`` -- the first matrix used in the construction
    - ``B`` -- the second matrix used in the construction
    - ``C`` -- the third matrix used in the construction
    - ``D`` -- the fourth matrix used in the construction

    TESTS::

        sage: from sage.combinat.matrices.hadamard_matrix import _construction_goethals_seidel_matrix
        sage: A = matrix([[1, -1], [1, 1]])
        sage: B = matrix([[-1, -1], [-1, 1]])
        sage: C = matrix([[1, 1], [1, -1]])
        sage: D = matrix([[-1, 1], [-1, 1]])
        sage: _construction_goethals_seidel_matrix(A, B, C, D)
        [ 1 -1|-1 -1| 1  1| 1 -1]
        [ 1  1| 1 -1|-1  1| 1 -1]
        [-----+-----+-----+-----]
        [ 1  1| 1 -1| 1  1| 1  1]
        [-1  1| 1  1|-1 -1|-1  1]
        [-----+-----+-----+-----]
        [-1 -1|-1 -1| 1 -1| 1  1]
        [ 1 -1| 1  1| 1  1|-1  1]
        [-----+-----+-----+-----]
        [-1  1|-1 -1|-1 -1| 1 -1]
        [-1  1| 1 -1| 1 -1| 1  1]
    """
    n = len(A[0])
    R = matrix(ZZ, n, n, lambda i, j: 1 if i+j == n-1 else 0)
    return block_matrix([[   A,    B*R,    C*R,    D*R],
                         [-B*R,      A, -D.T*R,  C.T*R],
                         [-C*R,  D.T*R,      A, -B.T*R],
                         [-D*R, -C.T*R,  B.T*R,      A]])


def hadamard_matrix_from_sds(n, existence=False, check=True):
    r"""
    Construction of Hadamard matrices from supplementary difference sets.

    Given four SDS with parameters `4-\{n; n_1, n_2, n_3, n_4; \lambda\}` with
    `n_1 + n_2 + n_3 + n_4 = n+\lambda` we can construct four (-1,+1) sequences `a_i = (a_{i,0},...,a_{i,n-1})`
    where `a_{i,j} = -1` iff `j \in S_i`. These will be the fist rows of four circulant
    matrices `A_1, A_2, A_3, A_4` which, when plugged into the Goethals-Seidel array, create an
    Hadamard matrix of order `4n` (see [Djo1994b]_).

    The supplementary difference sets are taken from
    :func:`sage.combinat.designs.difference_family.supplementary_difference_set`.

    INPUT:

    - ``n`` -- integer; the order of the matrix to be constructed
    - ``check`` -- boolean (default: ``True``); if ``True``, check that the matrix
      is a Hadamard before returning
    - ``existence`` -- boolean (default: ``False``); if ``True``, only check if the
      matrix exists

    OUTPUT:

    If ``existence=False``, returns the Hadamard matrix of order `n`. It raises
    an error if no data is available to construct the matrix of the given order,
    or if `n` is not a multiple of `4`.
    If ``existence=True``, returns a boolean representing whether the matrix
    can be constructed or not.

    EXAMPLES:

    By default The function returns the Hadamard matrix ::

        sage: from sage.combinat.matrices.hadamard_matrix import hadamard_matrix_from_sds
        sage: hadamard_matrix_from_sds(148)
        148 x 148 dense matrix over Integer Ring...

    If ``existence`` is set to True, the function returns a boolean ::

        sage: hadamard_matrix_from_sds(764, existence=True)
        True

    TESTS::

        sage: from sage.combinat.matrices.hadamard_matrix import hadamard_matrix_from_sds, is_hadamard_matrix
        sage: is_hadamard_matrix(hadamard_matrix_from_sds(172))
        True
        sage: hadamard_matrix_from_sds(64, existence=True)
        False
        sage: hadamard_matrix_from_sds(64)
        Traceback (most recent call last):
        ...
        ValueError: SDS of order 16 not yet implemented.
        sage: hadamard_matrix_from_sds(14)
        Traceback (most recent call last):
        ...
        ValueError: n must be a positive multiple of four
    """
    from sage.combinat.designs.difference_family import supplementary_difference_set_hadamard

    if n <= 0 or n % 4 != 0:
        raise ValueError('n must be a positive multiple of four')
    t = n // 4

    if existence:
        return supplementary_difference_set_hadamard(t, existence=True)

    _, [S1, S2, S3, S4] = supplementary_difference_set_hadamard(t, check=False)
    a = [-1 if i in S1 else 1 for i in range(t)]
    b = [-1 if i in S2 else 1 for i in range(t)]
    c = [-1 if i in S3 else 1 for i in range(t)]
    d = [-1 if i in S4 else 1 for i in range(t)]

    if n == 956:
        a, b, c, d = [-el for el in d], a, b, c

    A, B, C, D = map(matrix.circulant, [a, b, c, d])
    if check:
        assert A*A.T+B*B.T+C*C.T+D*D.T == 4*t*I(t)

    H = _construction_goethals_seidel_matrix(A, B, C, D)
    if check:
        assert is_hadamard_matrix(H)
    return H


def hadamard_matrix_cooper_wallis_construction(x1, x2, x3, x4, A, B, C, D, check=True):
    r"""
    Create a Hadamard matrix using the construction detailed in [CW1972]_.

    Given four circulant matrices `X_1`, X_2, X_3, X_4` of order `n` with entries (0, 1, -1)
    such that the entrywise product of two distinct matrices is always equal to `0` and that
    `\sum_{i=1}^{4}X_iX_i^\top = nI_n` holds, and four matrices `A, B, C, D` of order `m` with
    elements (1, -1) such that `MN^\top = NM^\top` for all distinct `M`, `N` and
    `AA^\top + BB^\top + CC^\top + DD^\top =  4mI_n` holds, we construct a Hadamard matrix
    of order `4nm`.

    INPUT:

    - ``x1`` -- list or vector; the first row of the circulant matrix `X_1`
    - ``x2`` -- list or vector; the first row of the circulant matrix `X_2`
    - ``x3`` -- list or vector; the first row of the circulant matrix `X_3`
    - ``x4`` -- list or vector; the first row of the circulant matrix `X_4`
    - ``A`` -- the matrix described above
    - ``B`` -- the matrix described above
    - ``C`` -- the matrix described above
    - ``D`` -- the matrix described above
    - ``check`` -- boolean (default: ``True``); if ``True``, check that the resulting
      matrix is Hadamard before returing it

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import hadamard_matrix_cooper_wallis_construction
        sage: from sage.combinat.t_sequences import T_sequences_smallcases
        sage: seqs = T_sequences_smallcases(19)
        sage: hadamard_matrix_cooper_wallis_construction(seqs[0], seqs[1], seqs[2], seqs[3], matrix([1]), matrix([1]), matrix([1]), matrix([1]))
        76 x 76 dense matrix over Integer Ring...

    TESTS::

        sage: from sage.combinat.matrices.hadamard_matrix import hadamard_matrix_cooper_wallis_construction, is_hadamard_matrix
        sage: seqs = T_sequences_smallcases(13)
        sage: H = hadamard_matrix_cooper_wallis_construction(seqs[0], seqs[1], seqs[2], seqs[3], matrix([1]), matrix([1]), matrix([1]), matrix([1]))
        sage: is_hadamard_matrix(H)
        True
        sage: len(H[0]) == 13*4*1
        True
        sage: hadamard_matrix_cooper_wallis_construction(seqs[0], seqs[1], seqs[2], seqs[3], matrix([1]), matrix([1, -1]), matrix([1]), matrix([1]))
        Traceback (most recent call last):
        ...
        AssertionError
        sage: hadamard_matrix_cooper_wallis_construction([1,-1], [1, 1], [1,1], [1,1], matrix([1]), matrix([1]), matrix([1]), matrix([1]))
        Traceback (most recent call last):
        ...
        AssertionError
    """

    n = len(x1)
    assert n == len(x2) == len(x3) == len(x4)

    X1, X2, X3, X4 = map(matrix.circulant, [x1, x2, x3, x4])

    matrices = [X1, X2, X3, X4]
    for i in range(4):
        for j in range(i+1, 4):
            assert matrices[i].elementwise_product(matrices[j]) == zero_matrix(n)
    assert X1*X1.T + X2*X2.T + X3*X3.T + X4*X4.T == n*I(n)

    m = len(A[0])
    assert m == len(B[0]) == len(C[0]) == len(D[0])
    will_matrices = [A, B, C, D]
    for i in range(4):
        for j in range(i+1, 4):
            assert will_matrices[i]*will_matrices[j].T == will_matrices[j]*will_matrices[i].T
    assert A*A.T + B*B.T + C*C.T + D*D.T == 4*m*I(m)

    e1 = _construction_goethals_seidel_matrix(X1, X2, X3, X4)
    e2 = _construction_goethals_seidel_matrix(X2, -X1, X4, -X3)
    e3 = _construction_goethals_seidel_matrix(X3, -X4, -X1, X2)
    e4 = _construction_goethals_seidel_matrix(X4, X3, -X2, -X1)

    H = e1.tensor_product(A) + e2.tensor_product(B) + e3.tensor_product(C) + e4.tensor_product(D)
    if check:
        assert is_hadamard_matrix(H)
    return H


def hadamard_matrix_cooper_wallis_smallcases(n, check=True, existence=False):
    r"""
    Construct Hadamard matrices using the Cooper-Wallis construction for some small values of `n`.

    This function calls the function :func:`hadamard_matrix_cooper_wallis_construction`
    with the appropriate arguments.
    It constructs the matrices `X_1`, `X_2`, `X_3`, `X_4` using either
    T-matrices or the T-sequences from :func:`sage.combinat.t_sequences.T_sequences_smallcases`.
    The matrices `A`, `B`, `C`, `D` are taken from :func:`williamson_type_quadruples_smallcases`.

    Data for T-matrices of order 67 is taken from [Saw1985]_.

    INPUT:

    - ``n`` -- integer; the order of the matrix to be constructed
    - ``check`` -- boolean (default: ``True``); if ``True``, check that the matrix
      is a Hadamard matrix before returning
    - ``existence`` -- boolean (default: ``False``); if ``True``, only check if
      the matrix exists

    OUTPUT:

    If ``existence=False``, returns the Hadamard matrix of order `n`. It raises an error if no data
    is available to construct the matrix of the given order.
    If ``existence=True``, returns a boolean representing whether the matrix can be constructed or not.

    .. SEEALSO::

        :func:`hadamard_matrix_cooper_wallis_construction`

    EXAMPLES:

    By default The function returns the Hadamard matrix ::

        sage: from sage.combinat.matrices.hadamard_matrix import hadamard_matrix_cooper_wallis_smallcases
        sage: hadamard_matrix_cooper_wallis_smallcases(28)
        28 x 28 dense matrix over Integer Ring...

    If ``existence`` is set to True, the function returns a boolean ::

        sage: hadamard_matrix_cooper_wallis_smallcases(20, existence=True)
        True

    TESTS::

        sage: from sage.combinat.matrices.hadamard_matrix import hadamard_matrix_cooper_wallis_smallcases, is_hadamard_matrix
        sage: is_hadamard_matrix(hadamard_matrix_cooper_wallis_smallcases(188))
        True
        sage: hadamard_matrix_cooper_wallis_smallcases(64, existence=True)
        False
        sage: hadamard_matrix_cooper_wallis_smallcases(64)
        Traceback (most recent call last):
        ...
        ValueError: The Cooper-Wallis construction for Hadamard matrices of order 64 is not yet implemented.
        sage: hadamard_matrix_cooper_wallis_smallcases(14)
        Traceback (most recent call last):
        ...
        AssertionError
    """
    assert n % 4 == 0 and n > 0

    db = {
        67: (
            [1, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0, -1, -1, -1, 0, 0, -1, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, -1, 0, 1, -1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
            [0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 1, 0, 0, -1, -1, -1, 0, 0, 0, 0, 0, -1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0],
            [0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, 0, 0, 1, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 0, 1, 0, 0, 0, 0, 0, 0],
            [0, 0, -1, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 1, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, -1, 0, 0, -1, 0, 0, 0, 0, -1, 0, -1, 1, 1, 0, 0, 0]
        )
    }

    for T_seq_len in divisors(n//4):
        will_size = n // (4*T_seq_len)
        if (T_seq_len in db or T_sequences_smallcases(T_seq_len, existence=True)) and williamson_type_quadruples_smallcases(will_size, existence=True):
            if existence:
                return True

            e1, e2, e3, e4 = None, None, None, None
            if T_seq_len in db:
                e1, e2, e3, e4 = db[T_seq_len]
            else:
                e1, e2, e3, e4 = T_sequences_smallcases(T_seq_len, check=False)

            will_matrices = williamson_type_quadruples_smallcases(will_size)
            A, B, C, D = map(matrix.circulant, will_matrices)
            M = hadamard_matrix_cooper_wallis_construction(e1, e2, e3, e4, A, B, C, D, check=False)

            if check:
                assert is_hadamard_matrix(M)
            return M

    if existence:
        return False
    raise ValueError("The Cooper-Wallis construction for Hadamard matrices of order %s is not yet implemented." % n)


def _get_baumert_hall_units(n, existence=False):
    r"""
    Construct Baumert-Hall units of size `n` from available 4-symbol `\delta` codes.

    The construction is detailed in Theroem 2 from [Tur1974]_, and is based on the
    Goethals-Seidel construction of Hadamard matrices.
    We need a 4-symbol `\delta` code to detail the first row of circulant matrices M1, M2, M3, M4
    used in the construction.

    INPUT:

    - ``n`` -- integer; the size of the Baumert-Hall units
    - ``existence`` -- boolean (default: ``False``); if ``True``, only check whether
      the units can be contructed

    OUTPUT:

        If ``existence=True``, return a boolean representing whether the Baumert-Hall
        units can be constructed. Otherwise, return a tuple containing the four
        Baumert-Hall units.

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import _get_baumert_hall_units
        sage: _get_baumert_hall_units(28)
        (28 x 28 dense matrix over Integer Ring,
         28 x 28 dense matrix over Integer Ring,
         28 x 28 dense matrix over Integer Ring,
         28 x 28 dense matrix over Integer Ring)

    TESTS::

        sage: _get_baumert_hall_units(116, existence=True)
        True
        sage: _get_baumert_hall_units(200, existence=True)
        False
        sage: _get_baumert_hall_units(15)
        Traceback (most recent call last):
        ...
        AssertionError
        sage: _get_baumert_hall_units(200)
        Traceback (most recent call last):
        ...
        ValueError: The Baumert-Hall units of size 200 have not yet been implemented
    """
    assert n % 4 == 0 and n > 0

    delta_codes_len = n//4
    if not four_symbol_delta_code_smallcases(delta_codes_len, existence=True):
        if existence:
            return False
        raise ValueError("The Baumert-Hall units of size %s have not yet been implemented" % n)

    if existence:
        return True

    T1, T2, T3, T4 = four_symbol_delta_code_smallcases(delta_codes_len)
    M1 = matrix.circulant(T1)
    M2 = matrix.circulant(T2)
    M3 = matrix.circulant(T3)
    M4 = matrix.circulant(T4)

    M1hat = matrix(ZZ, 0.25*(M1+M2+M3+M4))
    M2hat = matrix(ZZ, 0.25*(M1-M2-M3+M4))
    M3hat = matrix(ZZ, 0.25*(M1+M2-M3-M4))
    M4hat = matrix(ZZ, 0.25*(M1-M2+M3-M4))

    e1 = _construction_goethals_seidel_matrix(M1hat, -M2hat, -M3hat, -M4hat)
    e2 = _construction_goethals_seidel_matrix(M2hat, M1hat, M4hat, -M3hat)
    e3 = _construction_goethals_seidel_matrix(M3hat, -M4hat, M1hat, M2hat)
    e4 = _construction_goethals_seidel_matrix(M4hat, M3hat, -M2hat, M1hat)
    return e1, e2, e3, e4


def hadamard_matrix_turyn_type(a, b, c, d, e1, e2, e3, e4, check=True):
    r"""
    Construction of Turyn type Hadamard matrix.

    Given `n\times n` circulant matrices `A`, `B`, `C`, `D` with 1,-1 entries,
    satisfying `AA^\top + BB^\top + CC^\top + DD^\top = 4nI`, and a set of
    Baumert-Hall units of order `4t`, one can construct a Hadamard matrix of order
    `4tn` as detailed by Turyn in [Tur1974]_.

    INPUT:

    - ``a`` -- 1,-1 list; the 1st row of `A`
    - ``b`` -- 1,-1 list; the 1st row of `B`
    - ``d`` -- 1,-1 list; the 1st row of `C`
    - ``c`` -- 1,-1 list; the 1st row of `D`
    - ``e1`` -- matrix; the first Baumert-Hall unit
    - ``e2`` -- matrix; the second Baumert-Hall unit
    - ``e3`` -- matrix; the third Baumert-Hall unit
    - ``e4`` -- matrix; the fourth Baumert-Hall unit
    - ``check`` -- boolean (default: ``True``); whether to check that the output
      is a Hadamard matrix before returning it

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import hadamard_matrix_turyn_type, _get_baumert_hall_units
        sage: A, B, C, D = _get_baumert_hall_units(28)
        sage: hadamard_matrix_turyn_type([1], [1], [1], [1], A, B, C, D)
        28 x 28 dense matrix over Integer Ring...

    TESTS::

        sage: from sage.combinat.matrices.hadamard_matrix import hadamard_matrix_turyn_type, _get_baumert_hall_units, is_hadamard_matrix
        sage: A, B, C, D = _get_baumert_hall_units(12)
        sage: is_hadamard_matrix(hadamard_matrix_turyn_type([1], [1], [1], [1], A, B, C, D))
        True
        sage: hadamard_matrix_turyn_type([1, -1], [1], [1], [1], A, B, C, D)
        Traceback (most recent call last):
        ...
        AssertionError
        sage: hadamard_matrix_turyn_type([1, -1], [1, 1], [1, 1], [1, 1], A, B, C, D)
        Traceback (most recent call last):
        ...
        AssertionError
    """
    A, B, C, D = map(matrix.circulant, [a, b, c, d])

    n = len(a)
    assert len(a) == len(b) == len(c) == len(d)
    assert A*A.T+B*B.T+C*C.T+D*D.T == 4*n*I(n)

    t4 = len(e1[0])
    assert t4 % 4 == 0
    t = t4 // 4

    # Check that e1, e2, e3, e4 are valid Baumert-Hall units
    for i in range(t4):
        for j in range(t4):
            assert abs(e1[i, j]) + abs(e2[i, j]) + abs(e3[i, j]) + abs(e4[i, j]) == 1

    assert e1*e1.T == t*I(t4) and e2*e2.T == t*I(t4) and e3*e3.T == t*I(t4) and e4*e4.T == t*I(t4)

    units = [e1, e2, e3, e4]
    for i in range(len(units)):
        for j in range(i+1, len(units)):
            assert units[i]*units[j].T + units[j]*units[i].T == 0*I(t4)

    H = e1.tensor_product(A) + e2.tensor_product(B) + e3.tensor_product(C) + e4.tensor_product(D)
    if check:
        assert is_hadamard_matrix(H)
    return H


def turyn_type_hadamard_matrix_smallcases(n, existence=False, check=True):
    r"""
    Construct a Hadamard matrix of order `n` from available 4-symbol `\delta` codes and Williamson quadruples.

    The function looks for Baumert-Hall units and Williamson type matrices from
    :func:`four_symbol_delta_code_smallcases` and :func:`williamson_type_quadruples_smallcases`
    and use them to construct a Hadamard matrix with the Turyn construction
    defined in :func:`hadamard_matrix_turyn_type`.

    INPUT:

    - ``n`` -- integer; the order of the matrix to be constructed
    - ``existence`` -- boolean (default: ``False``); if ``True``, only check if
      the matrix exists
    - ``check`` -- boolean (default: ``True``); if ``True``, check that the matrix
      is a Hadamard matrix before returning

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import turyn_type_hadamard_matrix_smallcases
        sage: turyn_type_hadamard_matrix_smallcases(28, existence=True)
        True
        sage: turyn_type_hadamard_matrix_smallcases(28)
        28 x 28 dense matrix over Integer Ring...

    TESTS::

        sage: from sage.combinat.matrices.hadamard_matrix import turyn_type_hadamard_matrix_smallcases, is_hadamard_matrix
        sage: is_hadamard_matrix(turyn_type_hadamard_matrix_smallcases(236)) # long time
        True
        sage: turyn_type_hadamard_matrix_smallcases(64, existence=True)
        False
        sage: turyn_type_hadamard_matrix_smallcases(64)
        Traceback (most recent call last):
        ...
        ValueError: The Turyn type construction for Hadamard matrices of order 64 is not yet implemented.
    """
    assert n % 4 == 0 and n > 0

    for delta_code_len in divisors(n//4):
        units_size = delta_code_len*4
        will_size = n//units_size
        if _get_baumert_hall_units(units_size, existence=True) and williamson_type_quadruples_smallcases(will_size, existence=True):
            if existence:
                return True

            e1, e2, e3, e4 = _get_baumert_hall_units(units_size)
            a, b, c, d = williamson_type_quadruples_smallcases(will_size)
            return hadamard_matrix_turyn_type(a, b, c, d, e1, e2, e3, e4, check=check)

    if existence:
        return False
    raise ValueError("The Turyn type construction for Hadamard matrices of order %s is not yet implemented." % n)


def hadamard_matrix_spence_construction(n, existence=False, check=True):
    r"""
    Create a Hadamard matrix of order `n` using the Spence construction.

    This construction (detailed in [Spe1975]_), uses supplementary difference sets implemented in
    :func:`sage.combinat.designs.difference_family.supplementary_difference_set_from_rel_diff_set` to create the
    desired matrix.

    INPUT:

    - ``n`` -- integer; the order of the matrix to be constructed
    - ``existence`` -- boolean (default: ``False``); if ``True``, only check if
      the matrix exists
    - ``check`` -- bolean (default: ``True``); if ``True``, check that the matrix
      is a Hadamard matrix before returning

    OUTPUT:

    If ``existence=True``, returns a boolean representing whether the Hadamard
    matrix can be constructed. Otherwise, returns the Hadamard matrix, or raises
    an error if it cannot be constructed.

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import hadamard_matrix_spence_construction
        sage: hadamard_matrix_spence_construction(36)
        36 x 36 dense matrix over Integer Ring...

    If ``existence`` is ``True``, the function returns a boolean ::

        sage: hadamard_matrix_spence_construction(52, existence=True)
        True

    TESTS::

        sage: from sage.combinat.matrices.hadamard_matrix import is_hadamard_matrix
        sage: is_hadamard_matrix(hadamard_matrix_spence_construction(100))
        True
        sage: hadamard_matrix_spence_construction(48, existence=True)
        False
        sage: hadamard_matrix_spence_construction(48)
        Traceback (most recent call last):
        ...
        ValueError: The order 48 is not covered by Spence construction.
        sage: hadamard_matrix_spence_construction(5)
        Traceback (most recent call last):
        ...
        AssertionError
        sage: hadamard_matrix_spence_construction(0)
        Traceback (most recent call last):
        ...
        AssertionError
    """
    from sage.combinat.designs.difference_family import supplementary_difference_set_from_rel_diff_set

    assert n % 4 == 0 and n > 0

    q = n//4

    if existence:
        return supplementary_difference_set_from_rel_diff_set(q, existence=True)

    if not supplementary_difference_set_from_rel_diff_set(q, existence=True):
        raise ValueError(f'The order {n} is not covered by Spence construction.')

    G, [S1, S2, S3, S4] = supplementary_difference_set_from_rel_diff_set(q, check=False)

    Glist = list(G)
    A1 = matrix.circulant([1 if j in S1 else -1 for j in Glist])
    A2 = matrix.circulant([1 if j in S4 else -1 for j in Glist])
    A3 = matrix.circulant([1 if j in S3 else -1 for j in Glist])
    A4 = matrix.circulant([1 if j in S2 else -1 for j in Glist])

    P = matrix(ZZ, [[1 if (i + j) % (q-1) == 0 else 0 for i in range(1, q)] for j in range(1, q)])

    e = matrix([1]*(q-1))
    m1 = matrix([-1])
    p1 = matrix([1])
    H = block_matrix([[  p1,   m1,   p1,   p1,     e,       e,       e,       e],
                      [  p1,   p1,   m1,   p1,    -e,       e,      -e,       e],
                      [  m1,   p1,   p1,   p1,    -e,       e,       e,      -e],
                      [  m1,   m1,   m1,   p1,    -e,      -e,       e,       e],
                      [-e.T,  e.T,  e.T, -e.T,    A1,    A2*P,    A3*P,    A4*P],
                      [-e.T, -e.T,  e.T,  e.T, -A2*P,      A1, -A4.T*P,  A3.T*P],
                      [-e.T, -e.T, -e.T, -e.T, -A3*P,  A4.T*P,      A1, -A2.T*P],
                      [ e.T, -e.T,  e.T, -e.T, -A4*P, -A3.T*P,  A2.T*P,      A1]])
    if check:
        assert is_hadamard_matrix(H, verbose=True)

    return H


def is_hadamard_matrix(M, normalized=False, skew=False, verbose=False):
    r"""
    Test if ``M`` is a Hadamard matrix.

    INPUT:

    - ``M`` -- a matrix
    - ``normalized`` -- boolean (default: ``False``); whether to test if ``M``
      is a normalized  Hadamard matrix, i.e. has its first row/column filled with +1
    - ``skew`` -- boolean (default: ``False``); whether to test if ``M`` is a skew
      Hadamard matrix, i.e. `M=S+I` for `-S=S^\top`, and `I` the identity matrix
    - ``verbose`` -- boolean (default: ``False``); whether to be verbose when
      the matrix is not Hadamard

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import is_hadamard_matrix
        sage: h = matrix.hadamard(12)
        sage: is_hadamard_matrix(h)
        True
        sage: from sage.combinat.matrices.hadamard_matrix import skew_hadamard_matrix
        sage: h = skew_hadamard_matrix(12)
        sage: is_hadamard_matrix(h, skew=True)
        True
        sage: h = matrix.hadamard(12)
        sage: h[0,0] = 2
        sage: is_hadamard_matrix(h,verbose=True)
        The matrix does not only contain +1 and -1 entries, e.g. 2
        False
        sage: h = matrix.hadamard(12)
        sage: for i in range(12):
        ....:     h[i,2] = -h[i,2]
        sage: is_hadamard_matrix(h,verbose=True,normalized=True)
        The matrix is not normalized
        False

    TESTS::

        sage: h = matrix.hadamard(12)
        sage: is_hadamard_matrix(h, skew=True)
        False
        sage: is_hadamard_matrix(h, skew=True, verbose=True)
        The matrix is not skew
        False
        sage: h = skew_hadamard_matrix(12)
        sage: is_hadamard_matrix(h, skew=True, verbose=True)
        True
        sage: is_hadamard_matrix(h, skew=False, verbose=True)
        True
        sage: h = -h
        sage: is_hadamard_matrix(h, skew=True, verbose=True)
        The matrix is not skew - diagonal entries must be all 1
        False
        sage: is_hadamard_matrix(h, skew=False, verbose=True)
        True
        sage: h = skew_hadamard_matrix(20, skew_normalize=False)
        sage: is_hadamard_matrix(h, skew=True, normalized=True, verbose=True)
        The matrix is not skew-normalized
        False
        sage: from sage.combinat.matrices.hadamard_matrix import normalise_hadamard
        sage: h = normalise_hadamard(h, skew=True)
        sage: is_hadamard_matrix(h, skew=True, normalized=True, verbose=True)
        True
    """
    n = M.ncols()
    if n != M.nrows():
        if verbose:
            print("The matrix is not square ({}x{})".format(M.nrows(), n))
        return False

    if n == 0:
        return True

    for r in M:
        for v in r:
            if v*v != 1:
                if verbose:
                    print("The matrix does not only contain +1 and -1 entries, e.g. " + str(v))
                return False

    prod = (M*M.transpose()).dict()
    if (len(prod) != n or set(prod.values()) != {n} or any((i, i) not in prod for i in range(n))):
        if verbose:
            print("The product M*M.transpose() is not equal to nI")
        return False

    if normalized:
        if skew and (set(M.row(0)) != {1}):
            if verbose:
                print("The matrix is not skew-normalized")
            return False
        elif not skew and (set(M.row(0)) != {1} or set(M.column(0)) != {1}):
            if verbose:
                print("The matrix is not normalized")
            return False

    if skew:
        for i in range(n-1):
            for j in range(i+1, n):
                if M[i, j] != -M[j, i]:
                    if verbose:
                        print("The matrix is not skew")
                    return False
        for i in range(n):
            if M[i, i] != 1:
                if verbose:
                    print("The matrix is not skew - diagonal entries must be all 1")
                return False
    return True


def is_skew_hadamard_matrix(M, normalized=False, verbose=False):
    r"""
    Test if ``M`` is a skew Hadamard matrix.

    this is a wrapper around the function :func:`is_hadamard_matrix`

    INPUT:

    - ``M`` -- a matrix
    - ``normalized`` -- boolean (default: ``False``); whether to test if ``M``
      is a skew-normalized Hadamard matrix, i.e. has its first row filled with +1
    - ``verbose`` -- boolean (default: ``False``); whether to be verbose when the
      matrix is not skew Hadamard

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import is_skew_hadamard_matrix, skew_hadamard_matrix
        sage: h = matrix.hadamard(12)
        sage: is_skew_hadamard_matrix(h, verbose=True)
        The matrix is not skew
        False
        sage: h = skew_hadamard_matrix(12)
        sage: is_skew_hadamard_matrix(h)
        True
        sage: from sage.combinat.matrices.hadamard_matrix import normalise_hadamard
        sage: h = normalise_hadamard(skew_hadamard_matrix(12), skew=True)
        sage: is_skew_hadamard_matrix(h, verbose=True, normalized=True)
        True
    """
    return is_hadamard_matrix(M, skew=True, normalized=normalized, verbose=verbose)


@matrix_method
def hadamard_matrix(n, existence=False, check=True, construction_name=False):
    r"""
    Try to construct a Hadamard matrix using the available methods.

    Currently all orders `\le 1200` for which a construction is
    known are implemented. For `n > 1200`, only some orders are available.

    INPUT:

    - ``n`` -- integer; dimension of the matrix
    - ``existence`` -- boolean (default: ``False``); whether to build the matrix
      or merely query if a construction is available in Sage. When set to ``True``,
      the function returns:

        - ``True`` -- meaning that Sage knows how to build the matrix
        - ``Unknown`` -- meaning that Sage does not know how to build the
          matrix, although the matrix may exist (see :mod:`sage.misc.unknown`).
        - ``False`` -- meaning that the matrix does not exist

    - ``check`` -- boolean (default: ``True``); whether to check that output is
      correct before returning it. As this is expected to be useless, you may
      want to disable it whenever you want speed.
    - ``construction_name`` -- boolean (default: ``False``); if it is ``True``,
      ``existence`` is ``True``, and a matrix exists, output the construction name.
      It has no effect if ``existence`` is set to ``False``.

    EXAMPLES::

        sage: hadamard_matrix(12).det()
        2985984
        sage: 12^6
        2985984
        sage: hadamard_matrix(1)
        [1]
        sage: hadamard_matrix(2)
        [ 1  1]
        [ 1 -1]
        sage: hadamard_matrix(8) # random
        [ 1  1  1  1  1  1  1  1]
        [ 1 -1  1 -1  1 -1  1 -1]
        [ 1  1 -1 -1  1  1 -1 -1]
        [ 1 -1 -1  1  1 -1 -1  1]
        [ 1  1  1  1 -1 -1 -1 -1]
        [ 1 -1  1 -1 -1  1 -1  1]
        [ 1  1 -1 -1 -1 -1  1  1]
        [ 1 -1 -1  1 -1  1  1 -1]
        sage: hadamard_matrix(8).det() == 8^4
        True

    We note that :func:`hadamard_matrix` returns a normalised Hadamard matrix
    (the entries in the first row and column are all +1) ::

        sage: hadamard_matrix(12) # random
        [ 1  1| 1  1| 1  1| 1  1| 1  1| 1  1]
        [ 1 -1|-1  1|-1  1|-1  1|-1  1|-1  1]
        [-----+-----+-----+-----+-----+-----]
        [ 1 -1| 1 -1| 1  1|-1 -1|-1 -1| 1  1]
        [ 1  1|-1 -1| 1 -1|-1  1|-1  1| 1 -1]
        [-----+-----+-----+-----+-----+-----]
        [ 1 -1| 1  1| 1 -1| 1  1|-1 -1|-1 -1]
        [ 1  1| 1 -1|-1 -1| 1 -1|-1  1|-1  1]
        [-----+-----+-----+-----+-----+-----]
        [ 1 -1|-1 -1| 1  1| 1 -1| 1  1|-1 -1]
        [ 1  1|-1  1| 1 -1|-1 -1| 1 -1|-1  1]
        [-----+-----+-----+-----+-----+-----]
        [ 1 -1|-1 -1|-1 -1| 1  1| 1 -1| 1  1]
        [ 1  1|-1  1|-1  1| 1 -1|-1 -1| 1 -1]
        [-----+-----+-----+-----+-----+-----]
        [ 1 -1| 1  1|-1 -1|-1 -1| 1  1| 1 -1]
        [ 1  1| 1 -1|-1  1|-1  1| 1 -1|-1 -1]

    To find how the matrix is obtained, use ``construction_name`` ::

        sage: matrix.hadamard(476, existence=True, construction_name=True)
        'cooper-wallis small cases: 476'

    TESTS::

        sage: matrix.hadamard(10,existence=True)
        False
        sage: matrix.hadamard(12,existence=True)
        True
        sage: matrix.hadamard(668,existence=True)
        Unknown
        sage: matrix.hadamard(10)
        Traceback (most recent call last):
        ...
        ValueError: The Hadamard matrix of order 10 does not exist
        sage: matrix.hadamard(312, existence=True)
        True
        sage: matrix.hadamard(1904, existence=True)
        True
        sage: matrix.hadamard(324, existence=True)
        True
    """
    name = str(n)
    if construction_name:
        def report_name(nam):
            return nam
    else:
        def report_name(nam):
            return True

    if not (n % 4 == 0) and (n > 2):
        if existence:
            return False
        raise ValueError("The Hadamard matrix of order %s does not exist" % n)
    if n == 2:
        if existence:
            return report_name(name)
        M = matrix([[1, 1], [1, -1]])
    elif n == 1:
        if existence:
            return report_name(name)
        M = matrix([1])
    elif is_prime_power(n//2 - 1) and (n//2 - 1) % 4 == 1:
        name = "paleyII " + name
        if existence:
            return report_name(name)
        M = hadamard_matrix_paleyII(n)
    elif n == 4 or n % 8 == 0 and hadamard_matrix(n//2, existence=True) is True:
        name = "doubling " + name
        if existence:
            return report_name(name)
        had = hadamard_matrix(n//2, check=False)
        chad1 = matrix([list(r) + list(r) for r in had.rows()])
        mhad = (-1) * had
        R = len(had.rows())
        chad2 = matrix([list(had.rows()[i]) + list(mhad.rows()[i])
                       for i in range(R)])
        M = chad1.stack(chad2)
    elif is_prime_power(n - 1) and (n - 1) % 4 == 3:
        name = "paleyI " + name
        if existence:
            return report_name(name)
        M = hadamard_matrix_paleyI(n)
    elif williamson_hadamard_matrix_smallcases(n, existence=True):
        name = "small cases: " + name
        if existence:
            return report_name(name)
        M = williamson_hadamard_matrix_smallcases(n, check=False)
    elif n == 156:
        if existence:
            return report_name(name)
        M = hadamard_matrix_156()
    elif hadamard_matrix_cooper_wallis_smallcases(n, existence=True):
        name = "cooper-wallis small cases: " + name
        if existence:
            return report_name(name)
        M = hadamard_matrix_cooper_wallis_smallcases(n, check=False)
    elif turyn_type_hadamard_matrix_smallcases(n, existence=True):
        name = "turyn type small cases: " + name
        if existence:
            return report_name(name)
        M = turyn_type_hadamard_matrix_smallcases(n, check=False)
    elif hadamard_matrix_miyamoto_construction(n, existence=True):
        name = "miyamoto " + name
        if existence:
            return report_name(name)
        M = hadamard_matrix_miyamoto_construction(n, check=False)
    elif hadamard_matrix_from_sds(n, existence=True):
        name = "SDS " + name
        if existence:
            return report_name(name)
        M = hadamard_matrix_from_sds(n, check=False)
    elif hadamard_matrix_spence_construction(n, existence=True):
        name = "spence " + name
        if existence:
            return report_name(name)
        M = hadamard_matrix_spence_construction(n, check=False)
    elif skew_hadamard_matrix(n, existence=True) is True:
        name = "skew " + name
        if existence:
            return report_name(name)
        M = skew_hadamard_matrix(n, check=False)
    elif regular_symmetric_hadamard_matrix_with_constant_diagonal(n, 1, existence=True) is True:
        name = "RSHCD " + name
        if existence:
            return report_name(name)
        M = regular_symmetric_hadamard_matrix_with_constant_diagonal(n, 1)
    elif hadamard_matrix_from_symmetric_conference_matrix(n, existence=True) is True:
        name = "Construction from symmetric conference matrix " + name
        if existence:
            return report_name(name)
        M = hadamard_matrix_from_symmetric_conference_matrix(n, check=False)
    else:
        if existence:
            return Unknown
        raise ValueError("The Hadamard matrix of order %s is not yet implemented." % n)

    if check:
        assert is_hadamard_matrix(M)

    return M


def hadamard_matrix_www(url_file, comments=False):
    r"""
    Pull file from Sloane's database and return the corresponding Hadamard
    matrix as a Sage matrix.

    You must input a filename of the form "had.n.xxx.txt" as described
    on the webpage http://neilsloane.com/hadamard/, where
    "xxx" could be empty or a number of some characters.

    If ``comments=True`` then the "Automorphism..." line of the had.n.xxx.txt
    file is printed if it exists. Otherwise nothing is done.

    EXAMPLES::

        sage: hadamard_matrix_www("had.4.txt")      # optional - internet
        [ 1  1  1  1]
        [ 1 -1  1 -1]
        [ 1  1 -1 -1]
        [ 1 -1 -1  1]
        sage: hadamard_matrix_www("had.16.2.txt",comments=True)   # optional - internet
        Automorphism group has order = 49152 = 2^14 * 3
        [ 1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1]
        [ 1 -1  1 -1  1 -1  1 -1  1 -1  1 -1  1 -1  1 -1]
        [ 1  1 -1 -1  1  1 -1 -1  1  1 -1 -1  1  1 -1 -1]
        [ 1 -1 -1  1  1 -1 -1  1  1 -1 -1  1  1 -1 -1  1]
        [ 1  1  1  1 -1 -1 -1 -1  1  1  1  1 -1 -1 -1 -1]
        [ 1 -1  1 -1 -1  1 -1  1  1 -1  1 -1 -1  1 -1  1]
        [ 1  1 -1 -1 -1 -1  1  1  1  1 -1 -1 -1 -1  1  1]
        [ 1 -1 -1  1 -1  1  1 -1  1 -1 -1  1 -1  1  1 -1]
        [ 1  1  1  1  1  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1]
        [ 1  1  1  1 -1 -1 -1 -1 -1 -1 -1 -1  1  1  1  1]
        [ 1  1 -1 -1  1 -1  1 -1 -1 -1  1  1 -1  1 -1  1]
        [ 1  1 -1 -1 -1  1 -1  1 -1 -1  1  1  1 -1  1 -1]
        [ 1 -1  1 -1  1 -1 -1  1 -1  1 -1  1 -1  1  1 -1]
        [ 1 -1  1 -1 -1  1  1 -1 -1  1 -1  1  1 -1 -1  1]
        [ 1 -1 -1  1  1  1 -1 -1 -1  1  1 -1 -1 -1  1  1]
        [ 1 -1 -1  1 -1 -1  1  1 -1  1  1 -1  1  1 -1 -1]
    """
    n = eval(url_file.split(".")[1])
    rws = []
    url = "http://neilsloane.com/hadamard/" + url_file
    with urlopen(url) as f:
        s = [bytes_to_str(line) for line in f.readlines()]
    for i in range(n):
        line = s[i]
        rws.append([1 if line[j] == "+" else -1 for j in range(n)])
    if comments:
        lastline = s[-1]
        if lastline[0] == "A":
            print(lastline)
    return matrix(rws)


_rshcd_cache = {}


def regular_symmetric_hadamard_matrix_with_constant_diagonal(n, e, existence=False):
    r"""
    Return a Regular Symmetric Hadamard Matrix with Constant Diagonal.

    A Hadamard matrix is said to be *regular* if its rows all sum to the same
    value.

    For `\epsilon\in\{-1,+1\}`, we say that `M` is a `(n,\epsilon)-RSHCD` if
    `M` is a regular symmetric Hadamard matrix with constant diagonal
    `\delta\in\{-1,+1\}` and row sums all equal to `\delta \epsilon
    \sqrt(n)`. For more information, see [HX2010]_ or 10.5.1 in
    [BH2012]_. For the case `n=324`, see :func:`RSHCD_324` and [CP2016]_.

    INPUT:

    - ``n`` -- integer; side of the matrix
    - ``e`` -- `-1` or `+1`; the value of `\epsilon`

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import regular_symmetric_hadamard_matrix_with_constant_diagonal
        sage: regular_symmetric_hadamard_matrix_with_constant_diagonal(4,1)
        [ 1  1  1 -1]
        [ 1  1 -1  1]
        [ 1 -1  1  1]
        [-1  1  1  1]
        sage: regular_symmetric_hadamard_matrix_with_constant_diagonal(4,-1)
        [ 1 -1 -1 -1]
        [-1  1 -1 -1]
        [-1 -1  1 -1]
        [-1 -1 -1  1]

    Other hardcoded values::

        sage: for n,e in [(36,1),(36,-1),(100,1),(100,-1),(196, 1)]:  # long time
        ....:     print(repr(regular_symmetric_hadamard_matrix_with_constant_diagonal(n,e)))
        36 x 36 dense matrix over Integer Ring
        36 x 36 dense matrix over Integer Ring
        100 x 100 dense matrix over Integer Ring
        100 x 100 dense matrix over Integer Ring
        196 x 196 dense matrix over Integer Ring

        sage: for n,e in [(324,1),(324,-1)]: # not tested - long time, tested in RSHCD_324
        ....:     print(repr(regular_symmetric_hadamard_matrix_with_constant_diagonal(n,e)))
        324 x 324 dense matrix over Integer Ring
        324 x 324 dense matrix over Integer Ring

    From two close prime powers::

        sage: regular_symmetric_hadamard_matrix_with_constant_diagonal(64,-1)
        64 x 64 dense matrix over Integer Ring (use the '.str()' method to see the entries)

    From a prime power and a conference matrix::

        sage: regular_symmetric_hadamard_matrix_with_constant_diagonal(676,1)  # long time
        676 x 676 dense matrix over Integer Ring (use the '.str()' method to see the entries)

    Recursive construction::

        sage: regular_symmetric_hadamard_matrix_with_constant_diagonal(144,-1)
        144 x 144 dense matrix over Integer Ring (use the '.str()' method to see the entries)

    REFERENCE:

    - [BH2012]_

    - [HX2010]_
    """
    if existence and (n, e) in _rshcd_cache:
        return _rshcd_cache[n, e]

    from sage.graphs.strongly_regular_db import strongly_regular_graph

    def true():
        _rshcd_cache[n, e] = True
        return True

    M = None
    if abs(e) != 1:
        raise ValueError
    sqn = None
    if is_square(n):
        sqn = int(sqrt(n))
    if n < 0:
        if existence:
            return False
        raise ValueError
    elif n == 4:
        if existence:
            return true()
        if e == 1:
            M = J(4)-2*matrix(4, [[int(i+j == 3) for i in range(4)] for j in range(4)])
        else:
            M = -J(4)+2*I(4)
    elif n == 36:
        if existence:
            return true()
        if e == 1:
            M = strongly_regular_graph(36, 15, 6, 6).adjacency_matrix()
            M = J(36) - 2*M
        else:
            M = strongly_regular_graph(36, 14, 4, 6).adjacency_matrix()
            M = -J(36) + 2*M + 2*I(36)
    elif n == 100:
        if existence:
            return true()
        if e == -1:
            M = strongly_regular_graph(100, 44, 18, 20).adjacency_matrix()
            M = 2*M - J(100) + 2*I(100)
        else:
            M = strongly_regular_graph(100, 45, 20, 20).adjacency_matrix()
            M = J(100) - 2*M
    elif n == 196 and e == 1:
        if existence:
            return true()
        M = strongly_regular_graph(196, 91, 42, 42).adjacency_matrix()
        M = J(196) - 2*M
    elif n == 324:
        if existence:
            return true()
        M = RSHCD_324(e)
    elif (e == 1 and
          n % 16 == 0 and
          sqn is not None and
          is_prime_power(sqn - 1) and
          is_prime_power(sqn + 1)):
        if existence:
            return true()
        M = -rshcd_from_close_prime_powers(sqn)

    elif (e == 1 and
          sqn is not None and
          sqn % 4 == 2 and
          strongly_regular_graph(sqn-1, (sqn-2)//2, (sqn-6)//4,
            existence=True) is True and
          is_prime_power(ZZ(sqn + 1))):
        if existence:
            return true()
        M = rshcd_from_prime_power_and_conference_matrix(sqn+1)

    # Recursive construction: the Kronecker product of two RSHCD is a RSHCD
    else:
        from itertools import product
        for n1, e1 in product(divisors(n)[1:-1], [-1, 1]):
            e2 = e1*e
            n2 = n//n1
            if (regular_symmetric_hadamard_matrix_with_constant_diagonal(n1, e1, existence=True) is True and
                regular_symmetric_hadamard_matrix_with_constant_diagonal(n2, e2, existence=True)) is True:
                if existence:
                    return true()
                M1 = regular_symmetric_hadamard_matrix_with_constant_diagonal(n1, e1)
                M2 = regular_symmetric_hadamard_matrix_with_constant_diagonal(n2, e2)
                M = M1.tensor_product(M2)
                break

    if M is None:
        from sage.misc.unknown import Unknown
        _rshcd_cache[n, e] = Unknown
        if existence:
            return Unknown
        raise ValueError("I do not know how to build a {}-RSHCD".format((n, e)))

    assert M*M.transpose() == n*I(n)
    assert set(map(sum, M)) == {ZZ(e*sqn)}

    return M


def RSHCD_324(e):
    r"""
    Return a size 324x324 Regular Symmetric Hadamard Matrix with Constant Diagonal.

    We build the matrix `M` for the case `n=324`, `\epsilon=1` directly from
    :meth:`JankoKharaghaniTonchevGraph
    <sage.graphs.graph_generators.GraphGenerators.JankoKharaghaniTonchevGraph>`
    and for the case `\epsilon=-1` from the "twist" `M'` of `M`, using Lemma 11
    in [HX2010]_. Namely, it turns out that the matrix

    .. MATH::

        M'=\begin{pmatrix} M_{12} & M_{11}\\ M_{11}^\top & M_{21} \end{pmatrix},
        \quad\text{where}\quad
        M=\begin{pmatrix} M_{11} & M_{12}\\ M_{21} & M_{22} \end{pmatrix},

    and the `M_{ij}` are 162x162-blocks, also RSHCD, its diagonal blocks having zero row
    sums, as needed by [loc.cit.]. Interestingly, the corresponding
    `(324,152,70,72)`-strongly regular graph
    has a vertex-transitive automorphism group of order 2592, twice the order of the
    (intransitive) automorphism group of the graph corresponding to `M`. Cf. [CP2016]_.

    INPUT:

    - ``e`` -- `-1` or `+1`; the value of `\epsilon`

    TESTS::

        sage: from sage.combinat.matrices.hadamard_matrix import RSHCD_324, is_hadamard_matrix
        sage: for e in [1,-1]:  # long time
        ....:     M = RSHCD_324(e)
        ....:     print("{} {} {}".format(M==M.T,is_hadamard_matrix(M),all(M[i,i]==1 for i in range(324))))
        ....:     print(list(set(sum(x) for x in M)))
        True True True
        [18]
        True True True
        [-18]

    REFERENCE:

    - [CP2016]_
    """
    from sage.graphs.generators.smallgraphs import JankoKharaghaniTonchevGraph as JKTG
    M = JKTG().adjacency_matrix()
    M = J(324) - 2*M
    if e == -1:
        M1 = M[:162].T
        M2 = M[162:].T
        M11 = M1[:162]
        M12 = M1[162:].T
        M21 = M2[:162].T
        M = block_matrix([[M12, -M11], [-M11.T, M21]])
    return M


def _helper_payley_matrix(n, zero_position=True):
    r"""
    Return the matrix constructed in Lemma 1.19 page 291 of [SWW1972]_.

    This function return a `n^2` matrix `M` whose rows/columns are indexed by
    the element of a finite field on `n` elements `x_1,...,x_n`. The value
    `M_{i,j}` is equal to `\chi(x_i-x_j)`.

    The elements `x_1,...,x_n` are ordered in such a way that the matrix
    (respectively, its submatrix obtained by removing first row and first column in the case
    ``zero_position=False``) is symmetric with respect to its second diagonal.
    The matrix is symmetric if `n=4k+1`, and skew-symmetric otherwise.

    INPUT:

    - ``n`` -- an odd prime power
    - ``zero_position`` -- boolean (default: ``True``); if it is true, place 0
      of ``F_n`` in the middle, otherwise place it first

    .. SEEALSO::

        :func:`rshcd_from_close_prime_powers`

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import _helper_payley_matrix
        sage: _helper_payley_matrix(5)
        [ 0  1  1 -1 -1]
        [ 1  0 -1  1 -1]
        [ 1 -1  0 -1  1]
        [-1  1 -1  0  1]
        [-1 -1  1  1  0]

    TESTS::

        sage: _helper_payley_matrix(11,zero_position=True)
        [ 0 -1  1 -1 -1  1 -1  1  1  1 -1]
        [ 1  0 -1  1 -1 -1 -1 -1  1  1  1]
        [-1  1  0 -1  1  1 -1 -1 -1  1  1]
        [ 1 -1  1  0 -1  1  1 -1 -1 -1  1]
        [ 1  1 -1  1  0  1 -1  1 -1 -1 -1]
        [-1  1 -1 -1 -1  0  1  1  1 -1  1]
        [ 1  1  1 -1  1 -1  0 -1  1 -1 -1]
        [-1  1  1  1 -1 -1  1  0 -1  1 -1]
        [-1 -1  1  1  1 -1 -1  1  0 -1  1]
        [-1 -1 -1  1  1  1  1 -1  1  0 -1]
        [ 1 -1 -1 -1  1 -1  1  1 -1  1  0]
        sage: _helper_payley_matrix(11,zero_position=False)
        [ 0 -1  1 -1 -1 -1  1  1  1 -1  1]
        [ 1  0 -1  1 -1 -1 -1  1  1  1 -1]
        [-1  1  0 -1  1 -1 -1 -1  1  1  1]
        [ 1 -1  1  0 -1  1 -1 -1 -1  1  1]
        [ 1  1 -1  1  0 -1  1 -1 -1 -1  1]
        [ 1  1  1 -1  1  0 -1  1 -1 -1 -1]
        [-1  1  1  1 -1  1  0 -1  1 -1 -1]
        [-1 -1  1  1  1 -1  1  0 -1  1 -1]
        [-1 -1 -1  1  1  1 -1  1  0 -1  1]
        [ 1 -1 -1 -1  1  1  1 -1  1  0 -1]
        [-1  1 -1 -1 -1  1  1  1 -1  1  0]
    """
    from sage.rings.finite_rings.finite_field_constructor import FiniteField
    K = FiniteField(n, prefix='x')

    # Order the elements of K in K_list
    # so that K_list[i] = -K_list[n - i - 1]
    K_pairs = set(tuple(sorted([x, -x])) for x in K if x)
    K_list = [K.zero()] * n
    shift = 0 if zero_position else 1

    for i, (x, y) in enumerate(sorted(K_pairs)):
        K_list[i + shift] = x
        K_list[-i - 1] = y

    M = matrix(ZZ, n, n, [(1 if (x - y).is_square() else -1)
                          for x in K_list for y in K_list])
    M -= I(n)
    assert (M * J(n)).is_zero()
    assert M * M.transpose() == n * I(n) - J(n)
    return M


def rshcd_from_close_prime_powers(n):
    r"""
    Return a `(n^2,1)`-RSHCD when `n-1` and `n+1` are odd prime powers and `n=0\pmod{4}`.

    The construction implemented here appears in Theorem 4.3 from [GS1970]_.

    Note that the authors of [SWW1972]_ claim in Corollary 5.12 (page 342) to have
    proved the same result without the `n=0\pmod{4}` restriction with a *very*
    similar construction. So far, however, I (Nathann Cohen) have not been able
    to make it work.

    INPUT:

    - ``n`` -- integer congruent to `0\pmod{4}`

    .. SEEALSO::

        :func:`regular_symmetric_hadamard_matrix_with_constant_diagonal`

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import rshcd_from_close_prime_powers
        sage: rshcd_from_close_prime_powers(4)
        [-1 -1  1 -1  1 -1 -1  1 -1  1 -1 -1  1 -1  1 -1]
        [-1 -1 -1  1  1 -1 -1  1 -1 -1  1 -1 -1  1 -1  1]
        [ 1 -1 -1  1  1  1  1 -1 -1 -1 -1 -1 -1 -1  1 -1]
        [-1  1  1 -1  1  1 -1 -1 -1 -1 -1  1 -1 -1 -1  1]
        [ 1  1  1  1 -1 -1 -1 -1 -1 -1  1 -1  1 -1 -1 -1]
        [-1 -1  1  1 -1 -1  1 -1 -1  1 -1  1 -1  1 -1 -1]
        [-1 -1  1 -1 -1  1 -1 -1  1 -1  1 -1 -1  1  1 -1]
        [ 1  1 -1 -1 -1 -1 -1 -1 -1  1 -1 -1 -1  1  1  1]
        [-1 -1 -1 -1 -1 -1  1 -1 -1 -1  1  1  1 -1  1  1]
        [ 1 -1 -1 -1 -1  1 -1  1 -1 -1 -1  1  1  1 -1 -1]
        [-1  1 -1 -1  1 -1  1 -1  1 -1 -1 -1  1  1 -1 -1]
        [-1 -1 -1  1 -1  1 -1 -1  1  1 -1 -1  1 -1 -1  1]
        [ 1 -1 -1 -1  1 -1 -1 -1  1  1  1  1 -1 -1 -1 -1]
        [-1  1 -1 -1 -1  1  1  1 -1  1  1 -1 -1 -1 -1 -1]
        [ 1 -1  1 -1 -1 -1  1  1  1 -1 -1 -1 -1 -1 -1  1]
        [-1  1 -1  1 -1 -1 -1  1  1 -1 -1  1 -1 -1  1 -1]

    REFERENCE:

    - [SWW1972]_
    """
    if n % 4:
        raise ValueError("n(={}) must be congruent to 0 mod 4")

    a, b = sorted([n-1, n+1], key=lambda x: -x % 4)
    Sa = _helper_payley_matrix(a)
    Sb = _helper_payley_matrix(b)
    U = matrix(a, [[int(i+j == a-1) for i in range(a)] for j in range(a)])

    K = (U*Sa).tensor_product(Sb) + U.tensor_product(J(b)-I(b)) - J(a).tensor_product(I(b))

    F = lambda x: diagonal_matrix([-(-1)**i for i in range(x)])
    G = block_diagonal_matrix([J(1), I(a).tensor_product(F(b))])
    e = matrix(a*b, [1]*(a*b))
    H = block_matrix(2, [-J(1), e.transpose(), e, K])

    HH = G*H*G
    assert len(set(map(sum, HH))) == 1
    assert HH**2 == n**2*I(n**2)
    return HH


def williamson_goethals_seidel_skew_hadamard_matrix(a, b, c, d, check=True):
    r"""
    Williamson-Goethals-Seidel construction of a skew Hadamard matrix.

    Given `n\times n` (anti)circulant matrices `A`, `B`, `C`, `D` with 1,-1 entries,
    and satisfying `A+A^\top = 2I`, `AA^\top + BB^\top + CC^\top + DD^\top = 4nI`,
    one can construct a skew Hadamard matrix of order `4n`, cf. [GS70s]_.

    INPUT:

    - ``a`` -- 1,-1 list; the 1st row of `A`
    - ``b`` -- 1,-1 list; the 1st row of `B`
    - ``d`` -- 1,-1 list; the 1st row of `C`
    - ``c`` -- 1,-1 list; the 1st row of `D`
    - ``check`` -- boolean (default: ``True``); if ``True``, check that the
      resulting matrix is skew Hadamard before returning it

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import williamson_goethals_seidel_skew_hadamard_matrix as WGS
        sage: a = [ 1,  1, 1, -1,  1, -1,  1, -1, -1]
        sage: b = [ 1, -1, 1,  1, -1, -1,  1,  1, -1]
        sage: c = [-1, -1]+[1]*6+[-1]
        sage: d = [ 1,  1, 1, -1,  1,  1, -1,  1,  1]
        sage: M = WGS(a,b,c,d,check=True)

    REFERENCES:

    - [GS70s]_

    - [Wall71]_

    - [KoSt08]_
    """
    n = len(a)
    A, B, C, D = map(matrix.circulant, [a, b, c, d])
    if check:
        assert A*A.T+B*B.T+C*C.T+D*D.T == 4*n*I(n)
        assert A+A.T == 2*I(n)

    M = _construction_goethals_seidel_matrix(A, B, C, D)
    if check:
        assert is_hadamard_matrix(M, normalized=False, skew=True)
    return M


def skew_hadamard_matrix_spence_construction(n, check=True):
    r"""
    Construct skew Hadamard matrix of order `n` using Spence construction.

    This function will construct skew Hadamard matrix of order `n=2(q+1)` where `q` is
    a prime power with `q = 5` (mod 8). The construction is taken from [Spe1977]_, and the
    relative difference sets are constructed using :func:`sage.combinat.designs.difference_family.relative_difference_set_from_homomorphism`.

    INPUT:

    - ``n`` -- positive integer
    - ``check`` -- boolean (default: ``True``); if ``True``, check that the
      resulting matrix is Hadamard before returning it

    OUTPUT:

    If ``n`` satisfies the requirements described above, the function returns a
    `n\times n` Hadamard matrix. Otherwise, an exception is raised.

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import skew_hadamard_matrix_spence_construction
        sage: skew_hadamard_matrix_spence_construction(28)
        28 x 28 dense matrix over Integer Ring...

    TESTS::

        sage: from sage.combinat.matrices.hadamard_matrix import is_hadamard_matrix
        sage: is_hadamard_matrix(skew_hadamard_matrix_spence_construction(12, check=False), skew=True)
        True
        sage: is_hadamard_matrix(skew_hadamard_matrix_spence_construction(60, check=False), skew=True)
        True
        sage: skew_hadamard_matrix_spence_construction(31)
        Traceback (most recent call last):
        ...
        ValueError: The order 31 is not covered by the Spence construction.
        sage: skew_hadamard_matrix_spence_construction(16)
        Traceback (most recent call last):
        ...
        ValueError: The order 16 is not covered by the Spence construction.
    """
    q = n//2 - 1
    m = (q+1)//2
    if n % 4 != 0 or not is_prime_power(q) or q % 8 != 5:
        raise ValueError(f'The order {n} is not covered by the Spence construction.')

    G, D = relative_difference_set_from_homomorphism(q, 2, (q-1)//4, check=False, return_group=True)
    D_fixed = get_fixed_relative_difference_set(G, D)
    D_union = D_fixed + [q+1+el for el in D_fixed]
    D_union = list(set([el % (4*(q+1)) for el in D_union]))

    def find_a(i):
        for a in range(8):
            if (a*(q+1)//2+i) % 8 == 0:
                return a

    ai = [find_a(0), find_a(1), find_a(2), find_a(3)]
    P = PolynomialRing(ZZ, 'x')

    Tm = 0
    for i in range(m):
        Tm += P.monomial(i)

    Ds = [[] for i in range(4)]

    for el in D_union:
        i = el % 8
        if i > 3:
            continue
        Ds[i].append(((el + ai[i] * m) // 8) % m)

    psis = [0 for i in range(4)]
    for i in range(4):
        for el in Ds[i]:
            psis[i] += P.monomial(el)

    diffs = [(2*psis[i] - Tm).mod(P.monomial(m)-1) for i in range(4)]
    a = [-el for el in diffs[1].coefficients()]
    b = diffs[0].coefficients()
    c = diffs[2].coefficients()
    d = diffs[3].coefficients()

    return williamson_goethals_seidel_skew_hadamard_matrix(a, b, c, d, check=check)


def skew_hadamard_matrix_spence_1975(n, existence=False, check=True):
    r"""
    Construct a skew Hadamard matrix of order `n = 4(1 + q + q^2)` using the
    Spence construction.

    If `n = 4(1 + q + q^2)` where `q` is a prime power such that either
    `1 + q + q^2` is a prime congruent to `3, 5, 7 \mod 8` or `3 + 2q + 2q^2` is
    a prime power, then a skew Hadamard matrix of order `n` can be constructed using
    the Goethals Seidel array. The four matrices `A, B, C, D` plugged into the
    GS-array are created using complementary difference sets of order `1 + q + q^2`
    (which exist if `q` satisfies the conditions above), and a cyclic planar
    difference set with parameters `(1 + q^2 + q^4, 1 + q^2, 1)`.
    These are constructed by the functions :func:`sage.combinat.designs.difference_family.complementary_difference_sets`
    and :func:`sage.combinat.designs.difference_family.difference_family`.

    For more details, see [Spe1975b]_.

    INPUT:

    - ``n`` -- positive integer; the order of the matrix to be constructed
    - ``existence`` -- boolean (default: ``False``); if ``True``, only return
      whether the Hadamard matrix can be constructed
    - ``check`` -- boolean (default: ``True``); check that the result
      is a skew Hadamard matrix before returning it

    OUTPUT:

    If ``existence=False``, returns the skew Hadamard matrix of order `n`. It
    raises an error if `n` does not satisfy the required conditions.
    If ``existence=True``, returns a boolean representing whether the matrix
    can be constructed or not.

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import skew_hadamard_matrix_spence_1975
        sage: skew_hadamard_matrix_spence_1975(52)
        52 x 52 dense matrix over Integer Ring...

    If ``existence`` is True, the function returns a boolean::

        sage: skew_hadamard_matrix_spence_1975(52, existence=True)
        True

    TESTS::

        sage: from sage.combinat.matrices.hadamard_matrix import is_skew_hadamard_matrix
        sage: is_skew_hadamard_matrix(skew_hadamard_matrix_spence_1975(52, check=False))
        True
        sage: skew_hadamard_matrix_spence_1975(100, existence=True)
        False
        sage: skew_hadamard_matrix_spence_1975(31)
        Traceback (most recent call last):
        ...
        ValueError: n is not in the form 4*(1+q+q^2)
        sage: skew_hadamard_matrix_spence_1975(16)
        Traceback (most recent call last):
        ...
        ValueError: n is not in the form 4*(1+q+q^2)
        sage: skew_hadamard_matrix_spence_1975(292)
        Traceback (most recent call last):
        ...
        ValueError: q=8 is not a valid parameter for this construction
    """
    from sage.combinat.designs.difference_family import is_fixed_relative_difference_set, difference_family, complementary_difference_sets

    q = None
    m = n // 4
    for i in range(m):
        if 1 + i + i**2 == m and is_prime_power(i):
            q = i
            break

    if n % 4 != 0 or q is None:
        if existence:
            return False
        raise ValueError('n is not in the form 4*(1+q+q^2)')

    is_valid = (is_prime(m) and m % 8 in [3, 5, 7]) or is_prime_power(3 + 2*q + 2*q**2)
    if existence:
        return is_valid

    if not is_valid:
        raise ValueError(f'q={q} is not a valid parameter for this construction')

    G, sets = difference_family(1 + q**2 + q**4, 1 + q**2, 1)

    def get_fixed_set(s, G, q):
        for el in G:
            fixed_set = [el + x for x in s]
            if is_fixed_relative_difference_set(fixed_set, q):
                return fixed_set
        assert False

    D = get_fixed_set(sets[0], G, q)

    D1 = [d for d in D if d % (1 - q + q**2) == 0]
    Dnot = [el for el in D if el not in D1]
    D2 = []

    indices = set()
    for i in range(len(Dnot)):
        if i in indices:
            continue
        indices.add(i)
        for j in range(i+1, len(Dnot)):
            if j not in indices and Dnot[i] % m == Dnot[j] % m:
                indices.add(j)
                D2.append(Dnot[i])

    D2 = [d % m for d in D2]
    D1 = [d % m for d in D1]

    G2, U, V = complementary_difference_sets(m, check=False)
    G2list = list(G2)
    p = [-1 if j in U else 1 for j in G2list]
    q = [-1 if j in V else 1 for j in G2list]
    r = [1 if j % m in D1 or j % m in D2 else -1 for j in range(m)]
    s = [1 if j % m in D2 else -1 for j in range(m)]

    return williamson_goethals_seidel_skew_hadamard_matrix(p, q, r, s, check=check)


def GS_skew_hadamard_smallcases(n, existence=False, check=True):
    r"""
    Data for Williamson-Goethals-Seidel construction of skew Hadamard matrices.

    Here we keep the data for this construction.
    Namely, it needs 4 circulant matrices with extra properties, as described in
    :func:`sage.combinat.matrices.hadamard_matrix.williamson_goethals_seidel_skew_hadamard_matrix`
    Matrices are taken from:

    * `n=36, 52`: [GS70s]_
    * `n=92`: [Wall71]_
    * `n=188`: [Djo2008a]_
    * `n=236`: [FKS2004]_
    * `n=276`: [Djo2023a]_

    Additional data is obtained from skew supplementary difference sets contained in
    :func:`sage.combinat.designs.difference_family.skew_supplementary_difference_set`, using the
    construction described in [Djo1992a]_.

    INPUT:

    - ``n`` -- integer; the order of the matrix
    - ``existence`` -- boolean (default: ``True``); if ``True``, only check that
      we can do the construction
    - ``check`` -- boolean (default: ``False``); if ``True``, check the result

    TESTS::

        sage: from sage.combinat.matrices.hadamard_matrix import GS_skew_hadamard_smallcases, is_skew_hadamard_matrix
        sage: GS_skew_hadamard_smallcases(36)
        36 x 36 dense matrix over Integer Ring...
        sage: GS_skew_hadamard_smallcases(52)
        52 x 52 dense matrix over Integer Ring...
        sage: GS_skew_hadamard_smallcases(92)
        92 x 92 dense matrix over Integer Ring...
        sage: GS_skew_hadamard_smallcases(100)
        sage: is_skew_hadamard_matrix(GS_skew_hadamard_smallcases(324, check=False))
        True
        sage: is_skew_hadamard_matrix(GS_skew_hadamard_smallcases(156, check=False))
        True
    """
    WGS = williamson_goethals_seidel_skew_hadamard_matrix

    def pmtoZ(s):
        return [1 if x == '+' else -1 for x in s]

    db = {
        36: ['+++-+-+--', '+-++--++-', '--++++++-', '+++-++-++'],
        52: ['++++-++--+---', '-+-++----++-+', '--+-+++++-+++', '--+-+++++-+++'],
        92: ['+-------++-+-+--+++++++', '++--+--++++--++++--+--+', '++---+-+-+-++-+-+-+---+', '+----+--+--++--+--+----'],
        188: ['+----+----++-+-+---++-++--+--+++-+-+--++++-++++',
              '++--+---+------++------++-+-++--+-+-+----+---++',
              '+-+-++---++-+---+++---++-++-++-++-+++++-+-+----',
              '+++-++-+-+---+-+++++--+-----++---+--+++++--++-+'],
        236: ['+-+---+-+-++-++---+----++-----+++++--++++-+++--+--+-+-+++-+',
              '+-+---+-+-++-++---+----++-----+++++--++++-+++--+--+-+-+++-+',
              '+++-++----+++-+-+++--+--++------+---+-----+--+-+--+---+----',
              '++++++--+++--+---++-+-+-+---+-+----++++-++-+--++-+--+------'],
        276: ['+--+++--+-+++--+---++-+++++-+++-++-+--+---+-----+--+++-++---+-++---++',
              '+-++--+-+----++-+---++++-+---+-++++++++-+---+-++++---+-++----+-+--++-',
              '--+--+-++---+--++--+-+-+++-+--++---++++-+-+-+--+-++-+++++++--+--+++++',
              '-+---+++-----++---+++-+++--+++++--+---+-+-++++-++++-++-++-+-+++++++++']
    }

    if existence:
        return n in db or skew_supplementary_difference_set(n//4, existence=True)

    if n in db:
        a, b, c, d = map(pmtoZ, db[n])
        return WGS(a, b, c, d, check=check)

    if skew_supplementary_difference_set(n//4, existence=True):
        t = n//4

        G, [S1, S2, S3, S4] = skew_supplementary_difference_set(t, check=False, return_group=True)
        Glist = list(G)

        A = matrix([[-1 if y-x in S1 else +1 for y in Glist] for x in Glist])
        B = matrix([[-1 if y-x in S2 else +1 for y in Glist] for x in Glist])
        C = matrix([[-1 if y-x in S3 else +1 for y in Glist] for x in Glist])
        D = matrix([[-1 if y-x in S4 else +1 for y in Glist] for x in Glist])

        H = _construction_goethals_seidel_matrix(A, B, C, D)
        if check:
            assert is_skew_hadamard_matrix(H)
        return H

    return None


def skew_hadamard_matrix_from_orthogonal_design(n, existence=False, check=True):
    r"""
    Construct skew Hadamard matrices of order `mn(n - 1)` if suitable orthogonal designs exist.

    In [Seb1978]_ is proved that if amicable Hadamard matrices of order `n` and an orthogonal
    design of type `(1, m, mn - m - 1)` in order `mn` exist, then a skew Hadamard matrix
    of order `mn(n - 1)` can be constructed. The paper uses amicable orthogonal designs
    instead of amicable Hadamard matrices, but the two are equivalent (see [Seb2017]_).

    Amicable Hadamard matrices are constructed using :func:`amicable_hadamard_matrices`,
    and the orthogonal designs are constructed using the Goethals-Seidel array,
    with data taken from [Seb2017]_.

    INPUT:

    - ``n`` -- positive integer; the order of the matrix to be constructed
    - ``existence`` -- boolean (default: ``False``); if ``True``, only return
      whether the skew Hadamard matrix can be constructed
    - ``check`` -- boolean (default: ``True``); if ``True``, check that the result
      is a skew Hadamard matrix before returning it

    OUTPUT:

    If ``existence=False``, returns the skew Hadamard matrix of order `n`. It
    raises an error if a construction for order `n` is not yet implemented, or if
    `n` does not satisfy the constraint.
    If ``existence=True``, returns a boolean representing whether the matrix
    can be constructed or not.

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import skew_hadamard_matrix_from_orthogonal_design
        sage: skew_hadamard_matrix_from_orthogonal_design(756)
        756 x 756 dense matrix over Integer Ring...

    If ``existence`` is True, the function returns a boolean::

        sage: skew_hadamard_matrix_from_orthogonal_design(200, existence=True)
        False

    TESTS::

        sage: from sage.combinat.matrices.hadamard_matrix import is_skew_hadamard_matrix
        sage: is_skew_hadamard_matrix(skew_hadamard_matrix_from_orthogonal_design(756, check=False))
        True
        sage: skew_hadamard_matrix_from_orthogonal_design(31)
        Traceback (most recent call last):
        ...
        ValueError: n must be a multiple of 4
        sage: skew_hadamard_matrix_from_orthogonal_design(16)
        Traceback (most recent call last):
        ...
        NotImplementedError: orthogonal designs for matrix of order 16 not yet implemented
    """
    # We use value i to represent entries where variable x_i should be, and -i for -x_i
    orthogonal_designs = {
        (1, 1, 26): [[1, 3, 3, -3, 3, -3, -3],  [2, 3, 3, -3, 3, -3, -3],
                     [3, 3, 3, -3, 3, 3, 3], [3, 3, -3, -3, -3, 3, -3]]
    }

    if n % 4 != 0:
        raise ValueError('n must be a multiple of 4')

    m1, m2 = None, None
    for d in divisors(n)[1:-1]:
        if (n//d) % (d-1) != 0:
            continue
        d1 = n // (d*(d - 1))
        if (1, d1, d1*d - d1 - 1) in orthogonal_designs and amicable_hadamard_matrices(d, existence=True):
            m1 = d1
            m2 = d

    if m2 is None or m1 is None:
        if existence:
            return False
        raise NotImplementedError(f'orthogonal designs for matrix of order {n} not yet implemented')

    if existence:
        return True

    M, N = amicable_hadamard_matrices(m2, check=False)
    M = normalise_hadamard(M, skew=True)
    N = normalise_hadamard(N)

    P = M[1:, 1:] - I(m2 - 1)
    D = N[1:, 1:]

    A1, A2, A3, A4 = map(matrix.circulant, orthogonal_designs[(1, m1, m1*m2 - m1 - 1)])
    OD = _construction_goethals_seidel_matrix(A1, A2, A3, A4)

    blocks = {1: P, -1: -P, 2: J(m2 - 1), -2: -J(m2 - 1), 3: D, -3: -D}
    H = block_matrix([[blocks[el] for el in row] for row in OD]) + I(n)

    if check:
        assert is_skew_hadamard_matrix(H)

    return H


def skew_hadamard_matrix_from_complementary_difference_sets(n, existence=False, check=True):
    r"""
    Construct a skew Hadamard matrix of order `n=4(m+1)` from complementary difference sets.

    If `A, B` are complementary difference sets over a group of order `2m+1`, then
    they can be used to construct a skew Hadamard matrix, as described in [BS1969]_.

    The complementary difference sets are constructed using the function
    :func:`sage.combinat.designs.difference_family.complementary_difference_sets`.

    INPUT:

    - ``n`` -- positive integer; the order of the matrix to be constructed
    - ``existence`` -- boolean (default: ``False``); if ``True``, only return
      whether the skew Hadamard matrix can be constructed
    - ``check`` -- boolean (default: ``True``);  if ``True``, check that the
      result is a skew Hadamard matrix before returning it

    OUTPUT:

    If ``existence=False``, returns the skew Hadamard matrix of order `n`. It
    raises an error if `n` does not satisfy the required conditions.
    If ``existence=True``, returns a boolean representing whether the matrix
    can be constructed or not.

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import skew_hadamard_matrix_from_complementary_difference_sets
        sage: skew_hadamard_matrix_from_complementary_difference_sets(20)
        20 x 20 dense matrix over Integer Ring...
        sage: skew_hadamard_matrix_from_complementary_difference_sets(52, existence=True)
        True

    TESTS::

        sage: from sage.combinat.matrices.hadamard_matrix import is_skew_hadamard_matrix
        sage: is_skew_hadamard_matrix(skew_hadamard_matrix_from_complementary_difference_sets(24, check=False))
        True
        sage: is_skew_hadamard_matrix(skew_hadamard_matrix_from_complementary_difference_sets(12, check=False))
        True
        sage: skew_hadamard_matrix_from_complementary_difference_sets(31)
        Traceback (most recent call last):
        ...
        ValueError: n must be 1, 2 or a multiple of four.
        sage: skew_hadamard_matrix_from_complementary_difference_sets(100)
        Traceback (most recent call last):
        ...
        NotImplementedError: hadamard matrix of order 100 from complementary difference sets is not implemented yet
        sage: skew_hadamard_matrix_from_complementary_difference_sets(100, existence=True)
        False
    """

    if n <= 0 or (n > 2 and n % 4 != 0):
        raise ValueError('n must be 1, 2 or a multiple of four.')

    m = n//4 - 1

    if existence:
        return complementary_difference_sets(2*m+1, existence=True)

    if not complementary_difference_sets(2*m+1, existence=True):
        raise NotImplementedError(f'hadamard matrix of order {n} from complementary difference sets is not implemented yet')

    G, A, B = complementary_difference_sets(2*m+1, check=False)

    m = n//4 - 1
    Glist = list(G)

    S = [[0 for i in range(n)] for j in range(n)]
    for i in range(2*m + 1):
        for j in range(2*m + 1):
            S[2*m + 1 + i][2*m + 1 + j] = -1 if Glist[j] - Glist[i] in A else 1
            S[i][j] = -S[2*m + 1 + i][2*m + 1 + j]
            S[2*m + 1 + j][i] = -1 if Glist[j] - Glist[i] in B else 1
            S[i][2*m + 1 + j] = -S[2*m + 1 + j][i]
        S[4*m + 2][i] = -1
        S[4*m + 2][2*m + 1 + i] = 1
        S[i][4*m + 2] = 1
        S[i + 2*m + 1][4*m + 2] = -1
    for i in range(4*m + 3):
        S[4*m + 3][i] = 1
        S[i][4*m + 3] = -1
    for i in range(4*m + 4):
        S[i][i] = 1

    H = matrix(S)

    if check:
        assert is_hadamard_matrix(H, skew=True)
    return H


def skew_hadamard_matrix_whiteman_construction(n, existence=False, check=True):
    r"""
    Construct a skew Hadamard matrix of order `n=2(q+1)` where `q=p^t` is a prime power with `p \cong 5 \mod 8` and `t \cong 2 \mod 4`.

    Assuming `n` satisfies the conditions above, it is possible to construct two supplementary difference sets
    `A, B` (see [Whi1971]_), and these can be used to construct a skew Hadamard matrix, as described in [BS1969]_.

    INPUT:

    - ``n`` -- positive integer; the order of the matrix to be constructed
    - ``existence`` -- boolean (default: ``False``); if ``True``, only return
      whether the Hadamard matrix can be constructed
    - ``check`` -- boolean (default: ``True``); if ``True``, check that the result
      is a skew Hadamard matrix before returning it

    OUTPUT:

    If ``existence=False``, returns the skew Hadamard matrix of order `n`. It
    raises an error if `n` does not satisfy the required conditions.
    If ``existence=True``, returns a boolean representing whether the matrix can
    be constructed or not.

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import skew_hadamard_matrix_whiteman_construction
        sage: skew_hadamard_matrix_whiteman_construction(52)
        52 x 52 dense matrix over Integer Ring...
        sage: skew_hadamard_matrix_whiteman_construction(52, existence=True)
        True

    TESTS::

        sage: from sage.combinat.matrices.hadamard_matrix import is_hadamard_matrix
        sage: is_hadamard_matrix(skew_hadamard_matrix_whiteman_construction(52, check=False), skew=True)
        True
        sage: is_hadamard_matrix(skew_hadamard_matrix_whiteman_construction(340, check=False), skew=True)
        True
        sage: skew_hadamard_matrix_whiteman_construction(31)
        Traceback (most recent call last):
        ...
        ValueError: The order 31 is not covered by the Whiteman construction.
        sage: skew_hadamard_matrix_whiteman_construction(100)
        Traceback (most recent call last):
        ...
        ValueError: The order 100 is not covered by the Whiteman construction.
        sage: skew_hadamard_matrix_whiteman_construction(100, existence=True)
        False

    .. NOTE::

        A more general version of this construction is :func:`skew_hadamard_matrix_from_complementary_difference_sets`.
    """

    q = n // 2 - 1
    p, t = is_prime_power(q, get_data=True)

    is_order_valid = n % 4 == 0 and t > 0 and p % 8 == 5 and t % 4 == 2
    if existence:
        return is_order_valid

    if not is_order_valid:
        raise ValueError(f'The order {n} is not covered by the Whiteman construction.')

    from sage.rings.finite_rings.finite_field_constructor import GF
    G = GF(q)
    f = (q-1) // 8
    Cs = {i: [G.gen()**(8*s+i) for s in range(f)] for i in [0, 1, 2, 3, 6, 7]}
    A = Cs[0] + Cs[1] + Cs[2] + Cs[3]
    B = Cs[0] + Cs[1] + Cs[6] + Cs[7]

    m = n//4 - 1
    Glist = list(G)

    S = [[0 for i in range(n)] for j in range(n)]
    for i in range(2*m + 1):
        for j in range(2*m + 1):
            S[2*m + 1 + i][2*m + 1 + j] = -1 if Glist[j] - Glist[i] in A else 1
            S[i][j] = -S[2*m + 1 + i][2*m + 1 + j]
            S[2*m + 1 + j][i] = -1 if Glist[j] - Glist[i] in B else 1
            S[i][2*m + 1 + j] = -S[2*m + 1 + j][i]
        S[4*m + 2][i] = -1
        S[4*m + 2][2*m + 1 + i] = 1
        S[i][4*m + 2] = 1
        S[i + 2*m + 1][4*m + 2] = -1
    for i in range(4*m + 3):
        S[4*m + 3][i] = 1
        S[i][4*m + 3] = -1
    for i in range(4*m + 4):
        S[i][i] = 1

    H = matrix(S)

    if check:
        assert is_hadamard_matrix(H, skew=True)
    return H


def skew_hadamard_matrix_from_good_matrices(a, b, c, d, check=True):
    r"""
    Construct skew Hadamard matrix from good matrices.

    Given good matrices `A`, `B`, `C`, `D` (`A` circulant, `B, C, D` back-circulant) they can be used to construct
    a skew Hadamard matrix using the following block matrix (as described in [Sze1988]_):

    .. MATH::

        \left(\begin{array}{rrrr}
        A & B & C & D \\
        -B & A & D & -C \\
        -C & -D & A & B \\
        -D & C & -B & A
        \end{array}\right)

    INPUT:

    - ``a`` -- (1,-1) list; the 1st row of `A`
    - ``b`` -- (1,-1) list; the 1st row of `B`
    - ``d`` -- (1,-1) list; the 1st row of `C`
    - ``c`` -- (1,-1) list; the 1st row of `D`
    - ``check`` -- boolean (default: ``True``); if ``True``, check that the matrix
      is a skew Hadamard matrix before returning it

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import skew_hadamard_matrix_from_good_matrices
        sage: a, b, c, d = ([1, 1, 1, -1, -1], [1, -1, 1, 1, -1], [1, -1, -1, -1, -1], [1, -1, -1, -1, -1])
        sage: skew_hadamard_matrix_from_good_matrices(a, b, c, d)
        20 x 20 dense matrix over Integer Ring...

    TESTS::

        sage: from sage.combinat.matrices.hadamard_matrix import is_hadamard_matrix
        sage: a, b, c, d = ([1, 1, 1, -1, -1], [1, -1, 1, 1, -1], [1, -1, -1, -1, -1], [1, -1, -1, -1, -1])
        sage: is_hadamard_matrix(skew_hadamard_matrix_from_good_matrices(a, b, c, d, check=False), skew=True)
        True
        sage: a, b, c, d = ([1, 1, 1, -1, -1], [1, -1, 1, 1, -1], [1, -1, -1, -1, -1], [1, -1, -1, -1, 1])
        sage: skew_hadamard_matrix_from_good_matrices(a, b, c, d)
        Traceback (most recent call last):
        ...
        AssertionError
        sage: a, b, c, d = ([1, 1, 1], [1, -1, 1, 1, -1], [1, -1, -1, -1, -1], [1, -1, -1, -1, -1])
        sage: skew_hadamard_matrix_from_good_matrices(a, b, c, d)
        Traceback (most recent call last):
        ...
        AssertionError
    """
    n = len(a)
    m = (n-1) // 2

    assert len(a) == len(b) == len(c) == len(d)
    assert a[0] == 1 and b[0] == 1 and c[0] == 1 and d[0] == 1
    for i in range(1, m+1):
        assert a[i] == -a[n-i] and b[i] == b[n-i] and c[i] == c[n-i] and d[i] == d[n-i]

    def back_circulant(row):
        length = len(row)
        return matrix([[row[(j+i) % length] for j in range(length)] for i in range(length)])

    A = matrix.circulant(a)
    B = back_circulant(b)
    C = back_circulant(c)
    D = back_circulant(d)

    H = block_matrix([[ A,  B,  C,  D],
                      [-B,  A,  D, -C],
                      [-C, -D,  A,  B],
                      [-D,  C, -B,  A]])

    if check:
        assert is_hadamard_matrix(H, skew=True)
    return H


def skew_hadamard_matrix_from_good_matrices_smallcases(n, existence=False, check=True):
    r"""
    Construct skew Hadamard matrices from good matrices for some small values of `n=4m`, with `m` odd.

    The function stores good matrices of odd orders `\le 31`, taken from [Sze1988]_.
    These are used to create skew Hadamard matrices of order `4m`, `1 \le m \le 31` (`m` odd), using the function
    :func:`skew_hadamard_matrix_from_good_matrices`.

    ALGORITHM:

    Given four sequences (stored in ``E_sequences``) of length `m`, they can be used to construct four `E-sequences`
    of length `n=2m+1`, as follows:

    .. MATH::

        \begin{aligned}
        a &= 1, a_0, a_1, ..., a_m, -a_m, -a_{m-1}, ..., -a_0 \\
        b &= 1, b_0, b_1, ..., b_m, b_m, b_{m-1}, ..., b_0 \\
        c &= 1, c_0, c_1, ..., c_m, c_m, c_{m-1}, ..., c_0 \\
        d &= 1, d_0, d_1, ..., d_m, d_m, d_{m-1}, ..., d_0 \\
        \end{aligned}

    These E-sequences will be the first rows of the four good matrices needed to construct a skew Hadamard matrix
    of order `4n`.

    INPUT:

    - ``n`` -- integer; the order of the skew Hadamard matrix to be constructed
    - ``existence`` -- boolean (default:  ``False``); if ``True``, only return
      whether the Hadamard matrix can be constructed
    - ``check`` -- boolean (default: ``True``); if ``True``, check that the matrix
      is a Hadamard matrix before returning it

    OUTPUT:

    If ``existence=False``, returns the skew Hadamard matrix of order `n`. It
    raises an error if no data is available to construct the matrix of the given
    order.
    If ``existence=True``, returns a boolean representing whether the matrix can
    be constructed or not.

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import skew_hadamard_matrix_from_good_matrices_smallcases
        sage: skew_hadamard_matrix_from_good_matrices_smallcases(20)
        20 x 20 dense matrix over Integer Ring...
        sage: skew_hadamard_matrix_from_good_matrices_smallcases(20, existence=True)
        True

    TESTS::

        sage: from sage.combinat.matrices.hadamard_matrix import is_hadamard_matrix
        sage: is_hadamard_matrix(skew_hadamard_matrix_from_good_matrices_smallcases(28, check=False), skew=True)
        True
        sage: skew_hadamard_matrix_from_good_matrices_smallcases(140)
        Traceback (most recent call last):
        ...
        ValueError: The Good matrices of order 35 are not yet implemented.
        sage: skew_hadamard_matrix_from_good_matrices_smallcases(14)
        Traceback (most recent call last):
        ...
        ValueError: The skew Hadamard matrix of order 14 from good matrices does not exist.
        sage: skew_hadamard_matrix_from_good_matrices_smallcases(140, existence=True)
        False
        sage: skew_hadamard_matrix_from_good_matrices_smallcases(14, existence=True)
        False
    """
    E_sequences = {
        0: ['', '', '', ''],
        1: ['+', '-', '-', '+'],
        2: ['++', '-+', '--', '--'],
        3: ['++-', '++-', '+-+', '-++'],
        4: ['+++-', '+-+-', '---+', '++-+'],
        5: ['+-+--', '+++--', '-+++-', '---+-'],
        6: ['+-+---', '---+++', '+-+--+', '----+-'],
        7: ['+++++--', '-++--++', '----+-+', '-+---+-'],
        8: ['+--++-+-', '+--+----', '++----+-', '+---+-+-'],
        9: ['-+-----++', '+-+++++--', '-+----++-', '--+-+-++-'],
        10: ['+--+++++++', '++--++++-+', '--++-+-+-+', '---+++-+-+'],
        11: ['++-+-------', '+----+--+--', '+-+--++---+', '--++-+-+-++'],
        12: ['+-----+-+---', '+-++++-+-++-', '---+--++++--', '--+-+++--+--'],
        13: ['+---+-+--++-+', '+++---++-++-+', '+++-+++-++---', '+---++++-+-+-'],
        14: ['+--+----+-+-++', '+---++++-++--+', '+-+----++-+--+', '++++++---+-+-+'],
        15: ['+--++----+---+-', '-++-+---+-+++--', '++---+--+--+++-', '-++++++++--+-+-']
    }

    def pm_to_good_matrix(s, sign=1):
        e1 = [1 if x == '+' else -1 for x in s]
        e2 = [sign * (1 if x == '+' else -1) for x in s]
        e2.reverse()
        return [1] + e1 + e2

    if not (n % 4 == 0 and (n//4) % 2 == 1):
        if existence:
            return False
        raise ValueError("The skew Hadamard matrix of order %s from good matrices does not exist." % n)

    m = n//4
    l = (m-1) // 2

    if existence:
        return l in E_sequences

    if l not in E_sequences:
        raise ValueError("The Good matrices of order %s are not yet implemented." % m)

    e1, e2, e3, e4 = E_sequences[l]
    a = pm_to_good_matrix(e1, sign=-1)
    b = pm_to_good_matrix(e2)
    c = pm_to_good_matrix(e3)
    d = pm_to_good_matrix(e4)
    return skew_hadamard_matrix_from_good_matrices(a, b, c, d, check=check)


_skew_had_cache = {}


def skew_hadamard_matrix(n, existence=False, skew_normalize=True, check=True,
                         construction_name=False):
    r"""
    Try to construct a skew Hadamard matrix.

    A Hadamard matrix `H` is called skew if `H=S-I`, for `I` the identity matrix
    and `-S=S^\top`. Currently all orders `\le 1200` for which a construction is
    known are implemented. For `n > 1200`, only some orders are available.

    INPUT:

    - ``n`` -- integer; dimension of the matrix
    - ``existence`` -- boolean (default: ``False``); whether to build the matrix
      or merely query if a construction is available in Sage. When set to ``True``,
      the function returns:

        - ``True`` -- meaning that Sage knows how to build the matrix
        - ``Unknown`` -- meaning that Sage does not know how to build the
          matrix, but that the design may exist (see :mod:`sage.misc.unknown`).
        - ``False`` -- meaning that the matrix does not exist

    - ``skew_normalize`` -- boolean (default: ``True``); whether to make the 1st
      row all-one, and adjust the 1st column accordingly
    - ``check`` -- boolean (default: ``True``); whether to check that output is
      correct before returning it. As this is expected to be useless, you may
      want to disable it whenever you want speed.
    - ``construction_name`` -- boolean (default: ``False``); if it is ``True``,
      ``existence`` is ``True``, and a matrix exists, output the construction name.
      It has no effect if ``existence`` is set to ``False``.

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import skew_hadamard_matrix
        sage: skew_hadamard_matrix(12).det()
        2985984
        sage: 12^6
        2985984
        sage: skew_hadamard_matrix(1)
        [1]
        sage: skew_hadamard_matrix(2)
        [ 1  1]
        [-1  1]
        sage: skew_hadamard_matrix(196, existence=True, construction_name=True)
        'Williamson-Goethals-Seidel: 196'

    TESTS::

        sage: skew_hadamard_matrix(10,existence=True)
        False
        sage: skew_hadamard_matrix(12,existence=True)
        True
        sage: skew_hadamard_matrix(784,existence=True)
        True
        sage: skew_hadamard_matrix(10)
        Traceback (most recent call last):
        ...
        ValueError: A skew Hadamard matrix of order 10 does not exist
        sage: skew_hadamard_matrix(36)
        36 x 36 dense matrix over Integer Ring...
        sage: skew_hadamard_matrix(36)==skew_hadamard_matrix(36,skew_normalize=False)
        False
        sage: skew_hadamard_matrix(52)
        52 x 52 dense matrix over Integer Ring...
        sage: skew_hadamard_matrix(92)
        92 x 92 dense matrix over Integer Ring...
        sage: skew_hadamard_matrix(816)     # long time
        816 x 816 dense matrix over Integer Ring...
        sage: skew_hadamard_matrix(356)
        Traceback (most recent call last):
        ...
        ValueError: A skew Hadamard matrix of order 356 is not yet implemented.
        sage: skew_hadamard_matrix(356,existence=True)
        Unknown

    Check that :issue:`28526` is fixed::

        sage: skew_hadamard_matrix(0)
        Traceback (most recent call last):
        ...
        ValueError: parameter n must be strictly positive

    REFERENCES:

    - [Ha83]_
    """
    if n < 1:
        raise ValueError("parameter n must be strictly positive")

    def true(nam):
        _skew_had_cache[n] = nam
        if construction_name:
            return nam+": "+str(n)
        return True
    M = None
    name = ''
    if existence and n in _skew_had_cache:
        return true(_skew_had_cache[n])
    if not (n % 4 == 0) and (n > 2):
        if existence:
            return False
        raise ValueError("A skew Hadamard matrix of order %s does not exist" % n)
    if n == 2:
        if existence:
            return true(name)
        M = matrix([[1, 1], [-1, 1]])
    elif n == 1:
        if existence:
            return true(name)
        M = matrix([1])
    elif skew_hadamard_matrix_from_good_matrices_smallcases(n, existence=True):
        name = "good matrices small cases"
        if existence:
            return true(name)
        M = skew_hadamard_matrix_from_good_matrices_smallcases(n, check=False)
    elif is_prime_power(n - 1) and ((n - 1) % 4 == 3):
        name = "paleyI"
        if existence:
            return true(name)
        M = hadamard_matrix_paleyI(n, normalize=False)
    elif is_prime_power(n//2 - 1) and (n//2 - 1) % 8 == 5:
        name = "spence"
        if existence:
            return true(name)
        M = skew_hadamard_matrix_spence_construction(n, check=False)
    elif skew_hadamard_matrix_from_complementary_difference_sets(n, existence=True):
        name = "complementary DS"
        if existence:
            return true(name)
        M = skew_hadamard_matrix_from_complementary_difference_sets(n, check=False)
    elif skew_hadamard_matrix_spence_1975(n, existence=True):
        name = "spence [Spe1975b]"
        if existence:
            return true(name)
        M = skew_hadamard_matrix_spence_1975(n, check=False)
    elif skew_hadamard_matrix_from_orthogonal_design(n, existence=True):
        name = "orthogonal design"
        if existence:
            return true(name)
        M = skew_hadamard_matrix_from_orthogonal_design(n, check=False)
    elif n % 8 == 0:
        if skew_hadamard_matrix(n//2, existence=True) is True:  # (Lemma 14.1.6 in [Ha83]_)
            name = "doubling"
            if existence:
                return true(name)
            H = skew_hadamard_matrix(n//2, check=False)
            M = block_matrix([[H, H], [-H.T, H.T]])

        else:  # try Williamson construction (Lemma 14.1.5 in [Ha83]_)
            for d in divisors(n)[2:-2]:  # skip 1, 2, n/2, and n
                n1 = n//d
                if is_prime_power(d - 1) and (d % 4 == 0) and (n1 % 4 == 0)\
                    and skew_hadamard_matrix(n1, existence=True) is True:
                    from sage.arith.misc import factor
                    name = "williamson - Lemma 14.1.5 [Ha83] ("+str(factor(d-1))+","+str(n1)+") "
                    if existence:
                        return true(name)
                    H = skew_hadamard_matrix(n1, check=False)-I(n1)
                    U = matrix(ZZ, d, lambda i, j: -1 if i == j == 0 else
                                        1 if i == j == 1 or (i > 1 and j-1 == d-i)
                                          else 0)
                    A = block_matrix([[matrix([0]), matrix(ZZ, 1, d-1, [1]*(d-1))],
                                      [matrix(ZZ, d-1, 1, [-1]*(d-1)),
                                       _helper_payley_matrix(d-1, zero_position=0)]])+I(d)
                    M = A.tensor_product(I(n1))+(U*A).tensor_product(H)
                    break
    if M is None:  # try Williamson-Goethals-Seidel construction
        if GS_skew_hadamard_smallcases(n, existence=True) is True:
            name = "Williamson-Goethals-Seidel"
            if existence:
                return true(name)
            M = GS_skew_hadamard_smallcases(n)

        else:
            if existence:
                return Unknown
            raise ValueError("A skew Hadamard matrix of order %s is not yet implemented." % n)
    if skew_normalize:
        M = normalise_hadamard(M, skew=True)
    if check:
        assert is_hadamard_matrix(M, normalized=skew_normalize, skew=True)
    _skew_had_cache[n] = name
    return M


def symmetric_conference_matrix(n, check=True, existence=False):
    r"""
    Try to construct a symmetric conference matrix.

    A conference matrix is an `n\times n` matrix `C` with 0s on the main diagonal
    and 1s and -1s elsewhere, satisfying `CC^\top=(n-1)I`.
    If `C=C^\top` then `n \cong 2 \mod 4` and `C` is Seidel adjacency matrix of
    a graph, whose descendent graphs are strongly regular graphs with parameters
    `(n-1,(n-2)/2,(n-6)/4,(n-2)/4)`, see Sec.10.4 of [BH2012]_. Thus we build `C`
    from the Seidel adjacency matrix of the latter by adding row and column of 1s.

    INPUT:

    - ``n`` -- integer;  dimension of the matrix
    - ``check`` -- boolean (default: ``True``); whether to check that output is
      correct before returning it. As this is expected to be useless, you may
      want to disable it whenever you want speed.
    - ``existence`` -- boolean (default: ``False``); if true, only check that such
      a matrix exists.

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import symmetric_conference_matrix
        sage: C = symmetric_conference_matrix(10); C
        [ 0  1  1  1  1  1  1  1  1  1]
        [ 1  0 -1 -1  1 -1  1  1  1 -1]
        [ 1 -1  0 -1  1  1 -1 -1  1  1]
        [ 1 -1 -1  0 -1  1  1  1 -1  1]
        [ 1  1  1 -1  0 -1 -1  1 -1  1]
        [ 1 -1  1  1 -1  0 -1  1  1 -1]
        [ 1  1 -1  1 -1 -1  0 -1  1  1]
        [ 1  1 -1  1  1  1 -1  0 -1 -1]
        [ 1  1  1 -1 -1  1  1 -1  0 -1]
        [ 1 -1  1  1  1 -1  1 -1 -1  0]
        sage: C^2 == 9*identity_matrix(10) and C == C.T
        True
    """
    from sage.graphs.strongly_regular_db import strongly_regular_graph as srg
    try:
        m = srg(n-1, (n-2)/2, (n-6)/4, (n-2)/4, existence=existence)
    except ValueError:
        raise
    if existence:
        return m
    C = matrix([0]+[1]*(n-1)).stack(matrix([1]*(n-1)).stack(m.seidel_adjacency_matrix()).T)
    if check:
        assert (C == C.T and C**2 == (n-1)*I(n))
    return C


def szekeres_difference_set_pair(m, check=True):
    r"""
    Construct Szekeres `(2m+1,m,1)`-cyclic difference family.

    Let `4m+3` be a prime power. Theorem 3 in [Sz1969]_ contains a construction of a pair
    of *complementary difference sets* `A`, `B` in the subgroup `G` of the quadratic
    residues in `F_{4m+3}^*`. Namely `|A|=|B|=m`, `a\in A` whenever `a-1\in G`, `b\in B`
    whenever `b+1 \in G`. See also Theorem 2.6 in [SWW1972]_ (there the formula for `B` is
    correct, as opposed to (4.2) in [Sz1969]_, where the sign before `1` is wrong.

    In modern terminology, for `m>1` the sets `A` and `B` form a
    :func:`difference family<sage.combinat.designs.difference_family>` with parameters `(2m+1,m,1)`.
    I.e. each non-identity `g \in G` can be expressed uniquely as `xy^{-1}` for `x,y \in A` or `x,y \in B`.
    Other, specific to this construction, properties of `A` and `B` are: for `a` in `A` one has
    `a^{-1}` not in `A`, whereas for `b` in `B` one has `b^{-1}` in `B`.

    INPUT:

    - ``m`` -- integer; dimension of the matrix
    - ``check`` -- boolean (default: ``True``); whether to check `A` and `B` for
      correctness

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import szekeres_difference_set_pair
        sage: G,A,B=szekeres_difference_set_pair(6)
        sage: G,A,B=szekeres_difference_set_pair(7)

    REFERENCE:

    - [Sz1969]_
    """
    from sage.rings.finite_rings.finite_field_constructor import GF
    F = GF(4*m+3)
    t = F.multiplicative_generator()**2
    G = F.cyclotomic_cosets(t, cosets=[F.one()])[0]
    sG = set(G)
    A = [a for a in G if a - F.one() in sG]
    B = [b for b in G if b + F.one() in sG]
    if check:
        from itertools import product, chain
        assert (len(A) == len(B) == m)
        if m > 1:
            assert (sG == set([xy[0] / xy[1]
                              for xy in chain(product(A, A), product(B, B))]))
        assert (all(F.one() / b + F.one() in sG for b in B))
        assert (not any(F.one() / a - F.one() in sG for a in A))
    return G, A, B


def typeI_matrix_difference_set(G, A):
    r"""
    (1,-1)-incidence type I matrix of a difference set `A` in `G`.

    Let `A` be a difference set in a group `G` of order `n`. Return `n\times n`
    matrix `M` with `M_{ij}=1` if `A_i A_j^{-1} \in A`, and `M_{ij}=-1` otherwise.

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import szekeres_difference_set_pair
        sage: from sage.combinat.matrices.hadamard_matrix import typeI_matrix_difference_set
        sage: G,A,B=szekeres_difference_set_pair(2)
        sage: typeI_matrix_difference_set(G,A)
        [-1  1 -1 -1  1]
        [-1 -1 -1  1  1]
        [ 1  1 -1 -1 -1]
        [ 1 -1  1 -1 -1]
        [-1 -1  1  1 -1]
    """
    n = len(G)
    return matrix(n, n, lambda i, j: 1 if G[i]/G[j] in A else -1)


def rshcd_from_prime_power_and_conference_matrix(n):
    r"""
    Return a `((n-1)^2,1)`-RSHCD if `n` is prime power, and symmetric `(n-1)`-conference matrix exists.

    The construction implemented here is Theorem 16 (and Corollary 17) from [WW1972]_.

    In [SWW1972]_ this construction (Theorem 5.15 and Corollary 5.16)
    is reproduced with a typo. Note that [WW1972]_ refers to [Sz1969]_ for the construction,
    provided by :func:`szekeres_difference_set_pair`,
    of complementary difference sets, and the latter has a typo.

    From a :func:`symmetric_conference_matrix`, we only need the Seidel
    adjacency matrix of the underlying strongly regular conference (i.e. Paley
    type) graph, which we construct directly.

    INPUT:

    - ``n`` -- integer

    .. SEEALSO::

        :func:`regular_symmetric_hadamard_matrix_with_constant_diagonal`

    EXAMPLES:

    A 36x36 example ::

        sage: from sage.combinat.matrices.hadamard_matrix import rshcd_from_prime_power_and_conference_matrix
        sage: from sage.combinat.matrices.hadamard_matrix import is_hadamard_matrix
        sage: H = rshcd_from_prime_power_and_conference_matrix(7); H
        36 x 36 dense matrix over Integer Ring (use the '.str()' method to see the entries)
        sage: H == H.T and is_hadamard_matrix(H) and H.diagonal() == [1]*36 and list(sum(H)) == [6]*36
        True

    Bigger examples, only provided by this construction ::

        sage: H = rshcd_from_prime_power_and_conference_matrix(27)  # long time
        sage: H == H.T and is_hadamard_matrix(H)                    # long time
        True
        sage: H.diagonal() == [1]*676 and list(sum(H)) == [26]*676  # long time
        True

    In this example the conference matrix is not Paley, as 45 is not a prime power ::

        sage: H = rshcd_from_prime_power_and_conference_matrix(47)  # not tested (long time)

    REFERENCE:

    - [WW1972]_
    """
    from sage.graphs.strongly_regular_db import strongly_regular_graph as srg
    if is_prime_power(n) and 2 == (n-1) % 4:
        try:
            M = srg(n-2, (n-3)//2, (n-7)//4)
        except ValueError:
            return
        m = (n-3)//4
        Q, X, Y = szekeres_difference_set_pair(m)
        B = typeI_matrix_difference_set(Q, X)
        A = -typeI_matrix_difference_set(Q, Y)  # must be symmetric
        W = M.seidel_adjacency_matrix()
        f = J(1, 4*m+1)
        e = J(1, 2*m+1)
        JJ = J(2*m+1, 2*m+1)
        II = I(n-2)
        Ib = I(2*m+1)
        J4m = J(4*m+1, 4*m+1)
        H34 = -(B+Ib).tensor_product(W)+Ib.tensor_product(J4m)+(Ib-JJ).tensor_product(II)
        A_t_W = A.tensor_product(W)
        e_t_f = e.tensor_product(f)
        H = block_matrix([
            [J(1, 1),                f,                      e_t_f,                  -e_t_f],
            [f.T,                  J4m,     e.tensor_product(W-II),  e.tensor_product(W+II)],
            [ e_t_f.T, (e.T).tensor_product(W-II), A_t_W+JJ.tensor_product(II),         H34],
            [-e_t_f.T, (e.T).tensor_product(W+II), H34.T,      -A_t_W+JJ.tensor_product(II)]])
        return H


def are_amicable_hadamard_matrices(M, N, verbose=False):
    r"""
    Check if ``M`` and ``N`` are amicable Hadamard matrices.

    Two matrices `M` and `N` of order `n` are called amicable if they
    satisfy the following conditions (see [Seb2017]_):

    * `M` is a skew Hadamard matrix
    * `N` is a symmetric Hadamard matrix
    * `MN^T = NM^T`

    INPUT:

    - ``M`` -- a square matrix
    - ``N`` -- a square matrix
    - ``verbose`` -- boolean (default: ``False``); whether to be verbose when the
      matrices are not amicable Hadamard matrices

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import are_amicable_hadamard_matrices
        sage: M = matrix([[1, 1], [-1, 1]])
        sage: N = matrix([[1, 1], [1, -1]])
        sage: are_amicable_hadamard_matrices(M, N)
        True

    If ``verbose`` is true, the function will be verbose when returning False::

        sage: N = matrix([[1, 1], [1, 1]])
        sage: are_amicable_hadamard_matrices(M, N, verbose=True)
        The second matrix is not Hadamard
        False

    TESTS::

        sage: N = matrix.hadamard(12)
        sage: are_amicable_hadamard_matrices(M, N)
        False
    """
    if not is_skew_hadamard_matrix(M):
        if verbose:
            print('The first matrix is not skew Hadamard')
        return False

    if not is_hadamard_matrix(N):
        if verbose:
            print('The second matrix is not Hadamard')
        return False

    if not N.is_symmetric():
        if verbose:
            print('The second matrix is not symmetric')
        return False

    if len(M[0]) != len(N[0]):
        if verbose:
            print('M*N.transpose() is not equal to N*M.transpose()')
        return False

    return True


def amicable_hadamard_matrices_wallis(n, check=True):
    r"""
    Construct amicable Hadamard matrices of order `n = q + 1` where `q` is a prime power.

    If `q` is a prime power `\equiv 3 \mod 4`, then amicable Hadamard matrices
    of order `q+1` can be constructed as described in [Wal1970b]_.

    INPUT:

    - ``n`` -- integer; the order of the matrices to be constructed
    - ``check`` -- boolean (default: ``True``); if ``True``, check that the
      resulting matrices are amicable Hadamard before returning them

    OUTPUT:

    The function returns two amicable Hadamard matrices, or raises an error if such
    matrices cannot be created using this construction.

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import amicable_hadamard_matrices_wallis
        sage: M, N = amicable_hadamard_matrices_wallis(28)

    TESTS::

        sage: from sage.combinat.matrices.hadamard_matrix import are_amicable_hadamard_matrices
        sage: M, N = amicable_hadamard_matrices_wallis(20, check=False)
        sage: are_amicable_hadamard_matrices(M, N)
        True
        sage: amicable_hadamard_matrices_wallis(18)
        Traceback (most recent call last):
        ...
        ValueError: n must be a positive multiple of 4
        sage: amicable_hadamard_matrices_wallis(16)
        Traceback (most recent call last):
        ...
        ValueError: q = n-1 must be a prime power
    """
    if n % 4 != 0 or n < 0:
        raise ValueError('n must be a positive multiple of 4')

    q = n - 1
    if not is_prime_power(q):
        raise ValueError('q = n-1 must be a prime power')

    from sage.rings.finite_rings.finite_field_constructor import GF

    G = GF(q)

    ls1, ls2 = [], []
    elements_added = set()
    for el in G:
        if el == 0 or el in elements_added:
            continue
        elements_added.add(el)
        ls1.append(el)
        elements_added.add(-el)
        ls2 = [-el] + ls2
    Glist = [0] + ls1 + ls2

    squares = []
    for el in Glist:
        squares.append(el*el)

    def chi(el):
        if el == 0:
            return 0
        if el in squares:
            return 1
        return -1

    S = matrix([[chi(Glist[i] - Glist[j]) for j in range(q)] for i in range(q)])
    R = matrix([[1 if (i, j) == (0, 0) else 1 if j == q-i else 0 for j in range(q)] for i in range(q)])

    P = S + I(q)
    D = R + R*S

    e = matrix([1 for _ in range(q)])
    one = matrix([1])

    M = block_matrix([[ one, e],
                      [-e.T, P]])
    N = block_matrix([[-one, -e],
                      [-e.T,  D]])

    if check:
        assert are_amicable_hadamard_matrices(M, N)
    return M, N


def amicable_hadamard_matrices(n, existence=False, check=True):
    r"""
    Construct amicable Hadamard matrices of order ``n`` using the available methods.

    INPUT:

    - ``n`` -- positive integer; the order of the amicable Hadamard matrices
    - ``existence`` -- boolean (default: ``False``); if ``True``, only return
      whether amicable Hadamard matrices of order `n` can be constructed
    - ``check`` -- boolean (default: ``True``); if ``True``, check that the
      matrices are amicable Hadamard matrices before returning them

    OUTPUT:

    If ``existence`` is true, the function returns a boolean representing whether
    amicable Hadamard matrices of order `n` can be constructed.
    If ``existence`` is false, returns two amicable Hadamard matrices, or raises
    an error if the matrices cannot be constructed.

    EXAMPLES::

        sage: from sage.combinat.matrices.hadamard_matrix import amicable_hadamard_matrices
        sage: amicable_hadamard_matrices(2)
        (
        [ 1  1]  [ 1  1]
        [-1  1], [ 1 -1]
        )

    If ``existence`` is true, the function returns a boolean::

        sage: amicable_hadamard_matrices(16, existence=True)
        False

    TESTS::

        sage: M, N = amicable_hadamard_matrices(20)
        sage: amicable_hadamard_matrices(18)
        Traceback (most recent call last):
        ...
        ValueError: Hadamard matrices of order 18 do not exist
        sage: amicable_hadamard_matrices(16)
        Traceback (most recent call last):
        ...
        NotImplementedError: construction for amicable Hadamard matrices of order 16 not yet implemented
    """
    if not ((n % 4 == 0 and n > 2) or n in [1, 2]):
        if existence:
            return False
        raise ValueError(f"Hadamard matrices of order {n} do not exist")

    M = None
    N = None
    if n == 1:
        if existence:
            return True
        M = matrix([1])
        N = matrix([1])
    elif n == 2:
        if existence:
            return True
        M = matrix([[1, 1], [-1, 1]])
        N = matrix([[1, 1], [1, -1]])
    elif is_prime_power(n-1):
        if existence:
            return True
        M, N = amicable_hadamard_matrices_wallis(n, check=False)

    if existence:
        return False

    if M is None:
        raise NotImplementedError(f'construction for amicable Hadamard matrices of order {n} not yet implemented')

    if check:
        assert are_amicable_hadamard_matrices(M, N)
    return M, N
