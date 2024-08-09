r"""
Calculate symplectic bases for matrices over fields and the integers.

This module finds a symplectic basis for an anti-symmetric,
alternating matrix M defined over a field or the integers.

Anti-symmetric means that `M = -M^t`, where `M^t` denotes the
transpose of `M`.  Alternating means that the diagonal of `M` is
identically zero.

A symplectic basis is a basis of the form `e_1,
\ldots, e_j, f_1, \ldots, f_j, z_1, \ldots, z_k` such that

* `z_i M v^t` = 0 for all vectors `v`;
* `e_i M {e_j}^t = 0` for all `i, j`;
* `f_i M {f_j}^t = 0` for all `i, j`;
* `e_i M {f_j}^t = 0` for all `i` not equal `j`;

and such that the nonzero terms

* `e_i M {f_i}^t` are "as nice as possible": 1 over fields, or
  integers satisfying divisibility properties otherwise.

REFERENCES:

Bourbaki gives a nice proof that can be made constructive but is
not efficient (see Section 5, Number 1, Theorem 1, page 79):

Bourbaki, N.  Elements of Mathematics, Algebra III, Springer
Verlag 2007.

Kuperberg gives a more efficient and constructive exposition (see
Theorem 18).

Kuperberg, Greg.  Kasteleyn Cokernels.  Electr. J. Comb. 9(1), 2002.

.. TODO::

    The routine over the integers applies over general principal ideal
    domains.

.. WARNING::

    This code is not a good candidate for conversion to Cython.  The
    majority of the execution time is spent adding multiples of
    columns and rows, which is already fast.  It would be better to
    devise a better algorithm, perhaps modular or based on a fast
    ``smith_form`` implementation.

AUTHOR:

- Nick Alexander: initial implementation
- David Loeffler (2008-12-08): changed conventions for consistency with smith_form
"""

######################################################################
#       Copyright (C) 2008 William Stein
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#                  https://www.gnu.org/licenses/
######################################################################

from sage.rings.integer_ring import ZZ
from sage.rings.infinity import Infinity


def _inplace_move_to_positive_pivot(G, row, col, B, pivot):
    r"""
    Modify G so that v = G[row, col] appears at G[pivot, pivot+1].
    Modify B so that B * G * B.transpose() has the value v appearing
    at G[pivot, pivot+1].

    WARNING: not intended for external use!

    EXAMPLES::

        sage: from sage.matrix.symplectic_basis import _inplace_move_to_positive_pivot
        sage: E = matrix(ZZ, 4, 4, [0, 16, 0, 2, -16, 0, 0, -4, 0, 0, 0, 0, -2, 4, 0, 0]); E
        [  0  16   0   2]
        [-16   0   0  -4]
        [  0   0   0   0]
        [ -2   4   0   0]
        sage: B = copy(E.parent().one()); G = copy(E)
        sage: _inplace_move_to_positive_pivot(G, 0, 3, B, 0)
        sage: E[0, 3] == G[0, 1]
        True
        sage: G
        [  0   2   0  16]
        [ -2   0   0   4]
        [  0   0   0   0]
        [-16  -4   0   0]
        sage: B * E * B.transpose()
        [  0   2   0  16]
        [ -2   0   0   4]
        [  0   0   0   0]
        [-16  -4   0   0]
    """
    v = G[row, col]

    if (row, col) == (pivot, pivot+1):
        pass
    elif (row, col) == (pivot+1, pivot):
        B.swap_rows(pivot, pivot+1)
        G.swap_rows(pivot, pivot+1)
        G.swap_columns(pivot, pivot+1)
    elif row != pivot and row != pivot+1 and col != pivot and col != pivot+1:
        B.swap_rows(pivot, row)
        B.swap_rows(pivot+1, col)

        G.swap_rows(pivot, row)
        G.swap_rows(pivot+1, col)
        G.swap_columns(pivot, row)
        G.swap_columns(pivot+1, col)
    elif row == pivot:
        B.swap_rows(pivot+1, col)
        G.swap_rows(pivot+1, col)
        G.swap_columns(pivot+1, col)
    elif row == pivot+1:
        B.swap_rows(pivot, col)
        G.swap_rows(pivot, col)
        G.swap_columns(pivot, col)
    elif col == pivot:
        B.swap_rows(pivot+1, row)
        G.swap_rows(pivot+1, row)
        G.swap_columns(pivot+1, row)
    elif col == pivot+1:
        B.swap_rows(pivot, row)
        G.swap_rows(pivot, row)
        G.swap_columns(pivot, row)

    # all that swapping can switch the sign of a row
    if G[pivot, pivot+1] != v:
        B.swap_rows(pivot, pivot+1)
        G.swap_rows(pivot, pivot+1)
        G.swap_columns(pivot, pivot+1)


def symplectic_basis_over_field(M):
    r"""
    Find a symplectic basis for an anti-symmetric, alternating matrix
    M defined over a field.

    Returns a pair (F, C) such that the rows of C form a symplectic
    basis for M and ``F = C * M * C.transpose()``.

    Anti-symmetric means that `M = -M^t`.  Alternating means that the
    diagonal of `M` is identically zero.

    A symplectic basis is a basis of the form `e_1,
    \ldots, e_j, f_1, \ldots f_j, z_1, \ldots, z_k` such that

    * `z_i M v^t` = 0 for all vectors `v`;
    * `e_i M {e_j}^t = 0` for all `i, j`;
    * `f_i M {f_j}^t = 0` for all `i, j`;
    * `e_i M {f_i}^t = 1` for all `i`;
    * `e_i M {f_j}^t = 0` for all `i` not equal `j`.

    See the examples for a pictorial description of such a basis.

    EXAMPLES::

        sage: from sage.matrix.symplectic_basis import symplectic_basis_over_field

    A full rank exact example::

        sage: E = matrix(QQ, 8, 8, [0, -1/2, -2, 1/2, 2, 0, -2, 1, 1/2, 0, -1, -3, 0, 2, 5/2, -3, 2, 1, 0, 3/2, -1, 0, -1, -2, -1/2, 3, -3/2, 0, 1, 3/2, -1/2, -1/2, -2, 0, 1, -1, 0, 0, 1, -1, 0, -2, 0, -3/2, 0, 0, 1/2, -2, 2, -5/2, 1, 1/2, -1, -1/2, 0, -1, -1, 3, 2, 1/2, 1, 2, 1, 0]); E
        [   0 -1/2   -2  1/2    2    0   -2    1]
        [ 1/2    0   -1   -3    0    2  5/2   -3]
        [   2    1    0  3/2   -1    0   -1   -2]
        [-1/2    3 -3/2    0    1  3/2 -1/2 -1/2]
        [  -2    0    1   -1    0    0    1   -1]
        [   0   -2    0 -3/2    0    0  1/2   -2]
        [   2 -5/2    1  1/2   -1 -1/2    0   -1]
        [  -1    3    2  1/2    1    2    1    0]
        sage: F, C = symplectic_basis_over_field(E); F
        [ 0  0  0  0  1  0  0  0]
        [ 0  0  0  0  0  1  0  0]
        [ 0  0  0  0  0  0  1  0]
        [ 0  0  0  0  0  0  0  1]
        [-1  0  0  0  0  0  0  0]
        [ 0 -1  0  0  0  0  0  0]
        [ 0  0 -1  0  0  0  0  0]
        [ 0  0  0 -1  0  0  0  0]
        sage: F == C * E * C.transpose()
        True

    An example over a finite field::

        sage: E = matrix(GF(7), 8, 8,
        ....:            [0, -1/2, -2, 1/2, 2, 0, -2, 1, 1/2, 0, -1, -3, 0, 2, 5/2,
        ....:             -3, 2, 1, 0, 3/2, -1, 0, -1, -2, -1/2, 3, -3/2, 0, 1, 3/2,
        ....:             -1/2, -1/2, -2, 0, 1, -1, 0, 0, 1, -1, 0, -2, 0, -3/2, 0,
        ....:             0, 1/2, -2, 2, -5/2, 1, 1/2, -1, -1/2, 0, -1, -1, 3, 2,
        ....:             1/2, 1, 2, 1, 0]); E
        [0 3 5 4 2 0 5 1]
        [4 0 6 4 0 2 6 4]
        [2 1 0 5 6 0 6 5]
        [3 3 2 0 1 5 3 3]
        [5 0 1 6 0 0 1 6]
        [0 5 0 2 0 0 4 5]
        [2 1 1 4 6 3 0 6]
        [6 3 2 4 1 2 1 0]
        sage: F, C = symplectic_basis_over_field(E); F
        [0 0 0 0 1 0 0 0]
        [0 0 0 0 0 1 0 0]
        [0 0 0 0 0 0 1 0]
        [0 0 0 0 0 0 0 1]
        [6 0 0 0 0 0 0 0]
        [0 6 0 0 0 0 0 0]
        [0 0 6 0 0 0 0 0]
        [0 0 0 6 0 0 0 0]
        sage: F == C * E * C.transpose()
        True

    The tricky case of characteristic 2::

        sage: E = matrix(GF(2), 8, 8,
        ....:            [0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
        ....:            0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0,
        ....:            1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1,
        ....:            0, 1, 0, 0]); E
        [0 0 1 1 0 1 0 1]
        [0 0 0 0 0 0 0 0]
        [1 0 0 0 0 0 1 1]
        [1 0 0 0 0 0 0 1]
        [0 0 0 0 0 1 1 0]
        [1 0 0 0 1 0 1 1]
        [0 0 1 0 1 1 0 0]
        [1 0 1 1 0 1 0 0]
        sage: F, C = symplectic_basis_over_field(E); F
        [0 0 0 1 0 0 0 0]
        [0 0 0 0 1 0 0 0]
        [0 0 0 0 0 1 0 0]
        [1 0 0 0 0 0 0 0]
        [0 1 0 0 0 0 0 0]
        [0 0 1 0 0 0 0 0]
        [0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0]
        sage: F == C * E * C.transpose()
        True

    An inexact example::

        sage: E = matrix(RR, 8, 8, [0.000000000000000, 0.420674846479344, -0.839702420666807, 0.658715385244413, 1.69467394825853, -1.14718543053828, 1.03076138152950, -0.227739521708484, -0.420674846479344, 0.000000000000000, 0.514381455379082, 0.282194064028260, -1.38977093018412, 0.278305070958958, -0.0781320488361574, -0.496003664217833, 0.839702420666807, -0.514381455379082, 0.000000000000000, -0.00618222322875384, -0.318386939149028, -0.0840205427053993, 1.28202592892333, -0.512563654267693, -0.658715385244413, -0.282194064028260, 0.00618222322875384, 0.000000000000000, 0.852525732369211, -0.356957405431611, -0.699960114607661, 0.0260496330859998, -1.69467394825853, 1.38977093018412, 0.318386939149028, -0.852525732369211, 0.000000000000000, -0.836072521423577, 0.450137632758469, -0.696145287292091, 1.14718543053828, -0.278305070958958, 0.0840205427053993, 0.356957405431611, 0.836072521423577, 0.000000000000000, 0.214878541347751, -1.20221688928379, -1.03076138152950, 0.0781320488361574, -1.28202592892333, 0.699960114607661, -0.450137632758469, -0.214878541347751, 0.000000000000000, 0.785074452163036, 0.227739521708484, 0.496003664217833, 0.512563654267693, -0.0260496330859998, 0.696145287292091, 1.20221688928379, -0.785074452163036, 0.000000000000000]); E
        [   0.000000000000000    0.420674846479344   -0.839702420666807    0.658715385244413     1.69467394825853    -1.14718543053828     1.03076138152950   -0.227739521708484]
        [  -0.420674846479344    0.000000000000000    0.514381455379082    0.282194064028260    -1.38977093018412    0.278305070958958  -0.0781320488361574   -0.496003664217833]
        [   0.839702420666807   -0.514381455379082    0.000000000000000 -0.00618222322875384   -0.318386939149028  -0.0840205427053993     1.28202592892333   -0.512563654267693]
        [  -0.658715385244413   -0.282194064028260  0.00618222322875384    0.000000000000000    0.852525732369211   -0.356957405431611   -0.699960114607661   0.0260496330859998]
        [   -1.69467394825853     1.38977093018412    0.318386939149028   -0.852525732369211    0.000000000000000   -0.836072521423577    0.450137632758469   -0.696145287292091]
        [    1.14718543053828   -0.278305070958958   0.0840205427053993    0.356957405431611    0.836072521423577    0.000000000000000    0.214878541347751    -1.20221688928379]
        [   -1.03076138152950   0.0781320488361574    -1.28202592892333    0.699960114607661   -0.450137632758469   -0.214878541347751    0.000000000000000    0.785074452163036]
        [   0.227739521708484    0.496003664217833    0.512563654267693  -0.0260496330859998    0.696145287292091     1.20221688928379   -0.785074452163036    0.000000000000000]
        sage: F, C = symplectic_basis_over_field(E); F  # random
        [    0.000000000000000     0.000000000000000  2.22044604925031e-16 -2.22044604925031e-16      1.00000000000000     0.000000000000000     0.000000000000000 -3.33066907387547e-16]
        [    0.000000000000000  8.14814392305203e-17 -1.66533453693773e-16 -1.11022302462516e-16     0.000000000000000      1.00000000000000 -1.11022302462516e-16     0.000000000000000]
        [-5.27829526256056e-16 -2.40004077757759e-16  1.28373418199470e-16 -1.11022302462516e-16     0.000000000000000 -3.15483812822081e-16      1.00000000000000 -4.44089209850063e-16]
        [ 1.31957381564014e-16  1.41622049084608e-16 -6.68515202578511e-17 -3.95597468756028e-17 -4.85722573273506e-17 -5.32388011580111e-17 -1.31328455615552e-16      1.00000000000000]
        [    -1.00000000000000     0.000000000000000     0.000000000000000  4.85722573273506e-17     0.000000000000000 -5.55111512312578e-17 -1.11022302462516e-16  2.22044604925031e-16]
        [    0.000000000000000     -1.00000000000000     0.000000000000000 -2.77555756156289e-17  5.55111512312578e-17 -8.69223574327834e-17     0.000000000000000 -4.44089209850063e-16]
        [    0.000000000000000 -1.05042437087238e-17     -1.00000000000000  3.33066907387547e-16  1.11022302462516e-16 -1.18333563634309e-16  4.40064433050777e-17  2.22044604925031e-16]
        [ 5.27829526256056e-16  1.99901485752317e-16  1.65710718121313e-17     -1.00000000000000 -2.22044604925031e-16  5.52150940090699e-16 -3.93560383111738e-16  1.01155762925061e-16]
        sage: F == C * E * C.transpose()
        True
        sage: abs(F[0, 4] - 1) < 1e-10
        True
        sage: abs(F[4, 0] + 1) < 1e-10
        True

        sage: F.parent()
        Full MatrixSpace of 8 by 8 dense matrices over Real Field with 53 bits of precision
        sage: C.parent()
        Full MatrixSpace of 8 by 8 dense matrices over Real Field with 53 bits of precision
    """
    if not M.base_ring().is_field():
        raise ValueError("Can only find symplectic bases for matrices over fields")
    if not M.is_square():
        raise ValueError("Can only find symplectic bases for square matrices")
    if not M.transpose() + M == 0:
        raise ValueError("Can only find symplectic bases for anti-symmetric matrices")

    E = M.__copy__()
    n = E.nrows()
    for i in range(n):
        if not E[i, i].is_zero():
            raise ValueError("Can only find symplectic bases for alternating matrices")
    B = E.parent().one().__copy__()

    zeroes = []
    es = []
    fs = []
    pivot = 0
    while pivot < n:
        # find nonzero element in row
        found_i = None
        for i in range(pivot, n):
            if E[pivot, i] != 0:
                found_i = i
                break
        if found_i is None:
            zeroes.append(pivot)
            pivot += 1
            continue

        # move non-trivial entry to pivot position
        _inplace_move_to_positive_pivot(E, pivot, found_i, B, pivot)

        # scale row and col
        v = ZZ(1)/E[pivot, pivot+1]
        E.rescale_row(pivot, v)
        E.rescale_col(pivot, v)
        B.rescale_row(pivot, v)

        # use nonzero element to clean row pivot
        for i in range(pivot+2, n):
            v = - E[i, pivot] / E[pivot+1, pivot]
            if v != 0:
                E.add_multiple_of_row(i, pivot+1, v)
                E.add_multiple_of_column(i, pivot+1, v)
                B.add_multiple_of_row(i, pivot+1, v)

        # use nonzero element to clean row pivot+1
        for i in range(pivot+2, n):
            v = - E[i, pivot+1] / E[pivot, pivot+1]
            if v != 0:
                E.add_multiple_of_row(i, pivot, v)
                E.add_multiple_of_column(i, pivot, v)
                B.add_multiple_of_row(i, pivot, v)

        # record for basis reconstruction
        es.append(pivot)
        fs.append(pivot+1)
        pivot += 2

    C = B.matrix_from_rows(es + fs + zeroes)
    F = C * M * C.transpose()
    return F, C


def _smallest_element_position_or_None(E, pivot):
    r"""
    Return a tuple (row, col) such that E[row, col] is the smallest
    positive element of E, or None if E has no positive elements, and
    ``row >= pivot`` and ``col >= pivot``.

    WARNING: not intended for external use!

    EXAMPLES::

        sage: from sage.matrix.symplectic_basis import _smallest_element_position_or_None

        sage: E = matrix(ZZ, 4, 4, [0, 16, 0, 2, -16, 0, 0, -4, 0, 0, 0, 0, -2, 4, 0, 0]); E
        [  0  16   0   2]
        [-16   0   0  -4]
        [  0   0   0   0]
        [ -2   4   0   0]
        sage: _smallest_element_position_or_None(E, 0)
        (0, 3)
        sage: _smallest_element_position_or_None(E, 1)
        (3, 1)
        sage: _smallest_element_position_or_None(E, 2) is None
        True
    """
    found = None
    min = Infinity
    n = E.nrows()
    for i in range(pivot, n):
        for j in range(pivot, n):
            v = E[j, i]
            if 0 < v and v < min:
                min = v
                found = (j, i)
    return found


def symplectic_basis_over_ZZ(M):
    r"""
    Find a symplectic basis for an anti-symmetric, alternating matrix
    M defined over the integers.

    Returns a pair (F, C) such that the rows of C form a symplectic
    basis for M and F = C * M * C.transpose().

    Anti-symmetric means that `M = -M^t`.  Alternating means that the
    diagonal of `M` is identically zero.

    A symplectic basis is a basis of the form `e_1,
    \ldots, e_j, f_1, \ldots, f_j, z_1, \ldots, z_k` such that

    * `z_i M v^t` = 0 for all vectors `v`;
    * `e_i M {e_j}^t = 0` for all `i, j`;
    * `f_i M {f_j}^t = 0` for all `i, j`;
    * `e_i M {f_i}^t = d_i` for all `i`, where d_i are positive integers such that `d_{i} | d_{i+1}` for all `i`;
    * `e_i M {f_j}^t = 0` for all `i` not equal `j`.

    The ordering for the factors `d_{i} | d_{i+1}` and for the
    placement of zeroes was chosen to agree with the output of
    ``smith_form``.

    See the examples for a pictorial description of such a basis.

    EXAMPLES::

        sage: from sage.matrix.symplectic_basis import symplectic_basis_over_ZZ

    An example which does not have full rank::

        sage: E = matrix(ZZ, 4, 4, [0, 16, 0, 2, -16, 0, 0, -4, 0, 0, 0, 0, -2, 4, 0, 0]); E
        [  0  16   0   2]
        [-16   0   0  -4]
        [  0   0   0   0]
        [ -2   4   0   0]
        sage: F, C = symplectic_basis_over_ZZ(E)
        sage: F
        [ 0  2  0  0]
        [-2  0  0  0]
        [ 0  0  0  0]
        [ 0  0  0  0]
        sage: C * E * C.transpose() == F
        True

    A larger example::

        sage: E = matrix(ZZ, 8, 8, [0, 25, 0, 0, -37, -3, 2, -5, -25, 0, 1, -5, -54, -3, 3, 3, 0, -1, 0, 7, 0, -4, -20, 0, 0, 5, -7, 0, 0, 14, 0, -3, 37, 54, 0, 0, 0, 2, 3, -12, 3, 3, 4, -14, -2, 0, -3, 2, -2, -3, 20, 0, -3, 3, 0, -2, 5, -3, 0, 3, 12, -2, 2, 0]); E
        [  0  25   0   0 -37  -3   2  -5]
        [-25   0   1  -5 -54  -3   3   3]
        [  0  -1   0   7   0  -4 -20   0]
        [  0   5  -7   0   0  14   0  -3]
        [ 37  54   0   0   0   2   3 -12]
        [  3   3   4 -14  -2   0  -3   2]
        [ -2  -3  20   0  -3   3   0  -2]
        [  5  -3   0   3  12  -2   2   0]
        sage: F, C = symplectic_basis_over_ZZ(E)
        sage: F
        [     0      0      0      0      1      0      0      0]
        [     0      0      0      0      0      1      0      0]
        [     0      0      0      0      0      0      1      0]
        [     0      0      0      0      0      0      0  20191]
        [    -1      0      0      0      0      0      0      0]
        [     0     -1      0      0      0      0      0      0]
        [     0      0     -1      0      0      0      0      0]
        [     0      0      0 -20191      0      0      0      0]
        sage: F == C * E * C.transpose()
        True
        sage: E.smith_form()[0]
        [    1     0     0     0     0     0     0     0]
        [    0     1     0     0     0     0     0     0]
        [    0     0     1     0     0     0     0     0]
        [    0     0     0     1     0     0     0     0]
        [    0     0     0     0     1     0     0     0]
        [    0     0     0     0     0     1     0     0]
        [    0     0     0     0     0     0 20191     0]
        [    0     0     0     0     0     0     0 20191]

    An odd dimensional example::

        sage: E = matrix(ZZ, 5, 5, [0, 14, 0, -8, -2, -14, 0, -3, -11, 4, 0, 3, 0, 0, 0, 8, 11, 0, 0, 8, 2, -4, 0, -8, 0]); E
        [  0  14   0  -8  -2]
        [-14   0  -3 -11   4]
        [  0   3   0   0   0]
        [  8  11   0   0   8]
        [  2  -4   0  -8   0]
        sage: F, C = symplectic_basis_over_ZZ(E)
        sage: F
        [ 0  0  1  0  0]
        [ 0  0  0  2  0]
        [-1  0  0  0  0]
        [ 0 -2  0  0  0]
        [ 0  0  0  0  0]
        sage: F == C * E * C.transpose()
        True
        sage: E.smith_form()[0]
        [1 0 0 0 0]
        [0 1 0 0 0]
        [0 0 2 0 0]
        [0 0 0 2 0]
        [0 0 0 0 0]

        sage: F.parent()
        Full MatrixSpace of 5 by 5 dense matrices over Integer Ring
        sage: C.parent()
        Full MatrixSpace of 5 by 5 dense matrices over Integer Ring
    """
    if not M.is_square():
        raise ValueError("Can only find symplectic bases for square matrices")
    if not M.transpose() + M == 0:
        raise ValueError("Can only find symplectic bases for anti-symmetric matrices")

    E = M.__copy__().change_ring(ZZ)
    n = E.nrows()
    for i in range(n):
        if not E[i, i].is_zero():
            raise ValueError("Can only find symplectic bases for alternating matrices")
    B = E.parent().one().__copy__()

    zeroes = []
    ps = []
    pivot = 0
    while pivot < n:
        # find smallest element in matrix
        found = _smallest_element_position_or_None(E, pivot)
        if found is None:
            zeroes.append(pivot)
            pivot += 1
            continue
        _inplace_move_to_positive_pivot(E, found[0], found[1], B, pivot)

        # use nonzero element to clean row pivot
        all_zero = True
        u = E[pivot+1, pivot]
        for i in range(pivot+2, n):
            v, r = (-E[i, pivot]).quo_rem(u)
            if v != 0:
                all_zero = False
                E.add_multiple_of_row(i, pivot+1, v)
                E.add_multiple_of_column(i, pivot+1, v)
                B.add_multiple_of_row(i, pivot+1, v)

        # use nonzero element to clean row pivot+1
        u = E[pivot, pivot+1]
        for i in range(pivot+2, n):
            v, r = (-E[i, pivot+1]).quo_rem(u)
            if v != 0:
                all_zero = False
                E.add_multiple_of_row(i, pivot, v)
                E.add_multiple_of_column(i, pivot, v)
                B.add_multiple_of_row(i, pivot, v)

        if all_zero:
            # record for basis reconstruction
            ps.append((E[pivot, pivot+1], pivot))
            pivot += 2

    ps.sort()
    es = [p[1] for p in ps]
    fs = [p[1] + 1 for p in ps]
    C = B.matrix_from_rows(es + fs + zeroes)
    F = C * M * C.transpose()
    return F, C
