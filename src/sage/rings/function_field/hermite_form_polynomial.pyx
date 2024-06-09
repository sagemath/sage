r"""
Hermite form computation for function fields

This module provides an optimized implementation of the algorithm computing
Hermite forms of matrices over polynomials. This is the workhorse of the
function field machinery of Sage.

EXAMPLES::

    sage: P.<x> = PolynomialRing(QQ)
    sage: A = matrix(P,3,[-(x-1)^((i-j+1) % 3) for i in range(3) for j in range(3)])
    sage: A
    [        -x + 1             -1 -x^2 + 2*x - 1]
    [-x^2 + 2*x - 1         -x + 1             -1]
    [            -1 -x^2 + 2*x - 1         -x + 1]
    sage: from sage.rings.function_field.hermite_form_polynomial import reversed_hermite_form
    sage: B = copy(A)
    sage: U = reversed_hermite_form(B, transformation=True)
    sage: U * A == B
    True
    sage: B
    [x^3 - 3*x^2 + 3*x - 2                     0                     0]
    [                    0 x^3 - 3*x^2 + 3*x - 2                     0]
    [        x^2 - 2*x + 1                 x - 1                     1]

The function :func:`reversed_hermite_form` computes the reversed hermite form,
which is reversed both row-wise and column-wise from the usual hermite form.
Let us check it::

    sage: A.reverse_rows_and_columns()
    sage: C = copy(A.hermite_form())
    sage: C.reverse_rows_and_columns()
    sage: C
    [x^3 - 3*x^2 + 3*x - 2                    0                     0]
    [                    0 x^3 - 3*x^2 + 3*x - 2                     0]
    [        x^2 - 2*x + 1                 x - 1                     1]
    sage: C == B
    True

AUTHORS:

- Kwankyu Lee (2021-05-21): initial version
"""

# ****************************************************************************
#       Copyright (C) 2021 Kwankyu Lee <ekwankyu@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.matrix.matrix cimport Matrix
from sage.rings.polynomial.polynomial_element cimport Polynomial
from sage.matrix.constructor import identity_matrix


def reversed_hermite_form(Matrix mat, bint transformation=False):
    """
    Transform the matrix in place to reversed hermite normal form and
    optionally return the transformation matrix.

    INPUT:

    - ``transformation`` -- boolean (default: ``False``); if ``True``,
      return the transformation matrix

    EXAMPLES::

        sage: from sage.rings.function_field.hermite_form_polynomial import reversed_hermite_form
        sage: P.<x> = PolynomialRing(QQ)
        sage: A = matrix(P,3,[-(x-1)^((i-2*j) % 4) for i in range(3) for j in range(3)])
        sage: A
        [                    -1         -x^2 + 2*x - 1                     -1]
        [                -x + 1 -x^3 + 3*x^2 - 3*x + 1                 -x + 1]
        [        -x^2 + 2*x - 1                     -1         -x^2 + 2*x - 1]
        sage: B = copy(A)
        sage: U = reversed_hermite_form(B, transformation=True)
        sage: U * A == B
        True
        sage: B
        [                        0                         0                         0]
        [                        0 x^4 - 4*x^3 + 6*x^2 - 4*x                         0]
        [                        1             x^2 - 2*x + 1                         1]
    """
    cdef Matrix A = mat
    cdef Matrix U

    cdef Py_ssize_t m = A.nrows()
    cdef Py_ssize_t n = A.ncols()
    cdef Py_ssize_t i = m - 1
    cdef Py_ssize_t j = n - 1
    cdef Py_ssize_t k, l, c, ip

    cdef Polynomial a, b, d, p, q, e, f, Aic, Alc, Uic, Ulc
    cdef Polynomial zero = A.base_ring().zero()

    cdef int di, dk

    cdef list pivot_cols = []

    if transformation:
        U = <Matrix> identity_matrix(A.base_ring(), m)

    while j >= 0:
        # find the first row with nonzero entry in jth column
        k = i
        while k >= 0 and not A.get_unsafe(k, j):
            k -= 1
        if k >= 0:
            # swap the kth row with the ith row
            if k < i:
                A.swap_rows_c(i, k)
                if transformation:
                    U.swap_rows_c(i, k)
            k -= 1
            # put the row with the smalllest degree to the ith row
            di = A.get_unsafe(i, j).degree()
            while k >= 0 and A.get_unsafe(k, j):
                dk = A.get_unsafe(k, j).degree()
                if dk < di:
                    A.swap_rows_c(i, k)
                    di = dk
                    if transformation:
                        U.swap_rows_c(i, k)
                k -= 1
            l = i - 1
            while True:
                # find a row with nonzero entry in the jth column
                while l >= 0 and not A.get_unsafe(l, j):
                    l -= 1
                if l < 0:
                    break

                a = <Polynomial> A.get_unsafe(i, j)
                b = <Polynomial> A.get_unsafe(l, j)
                d, p, q = a.xgcd(b) # p * a + q * b = d = gcd(a,b)
                e = a // d
                f = -(b // d)

                A.set_unsafe(i, j, d)
                A.set_unsafe(l, j, zero)

                for c in range(j):
                    Aic = <Polynomial> A.get_unsafe(i, c)
                    Alc = <Polynomial> A.get_unsafe(l, c)
                    A.set_unsafe(i, c, p * Aic + q * Alc)
                    A.set_unsafe(l, c, f * Aic + e * Alc)
                if transformation:
                    for c in range(m):
                        Uic = <Polynomial> U.get_unsafe(i, c)
                        Ulc = <Polynomial> U.get_unsafe(l, c)
                        U.set_unsafe(i, c, p * Uic + q * Ulc)
                        U.set_unsafe(l, c, f * Uic + e * Ulc)
            pivot_cols.append(j)
            i -= 1
        j -= 1

    # reduce entries below pivots
    for i in range(len(pivot_cols)):
        j = pivot_cols[i]
        ip = m - 1 - i
        pivot = <Polynomial> A.get_unsafe(ip, j)

        # normalize the leading coefficient to one
        coeff = ~pivot.lc()
        for c in range(j + 1):
            A.set_unsafe(ip, c, <Polynomial> A.get_unsafe(ip, c) * coeff)
        if transformation:
            for c in range(m):
                U.set_unsafe(ip, c, <Polynomial> U.get_unsafe(ip, c) * coeff)

        pivot = <Polynomial> A.get_unsafe(ip,j)
        for k in range(ip + 1, m):
            q = -(<Polynomial> A.get_unsafe(k, j) // pivot)
            if q:
                for c in range(j + 1):
                    A.set_unsafe(k, c, <Polynomial> A.get_unsafe(k, c)
                                 + q * <Polynomial> A.get_unsafe(ip, c))
                if transformation:
                    for c in range(m):
                        U.set_unsafe(k, c, <Polynomial> U.get_unsafe(k, c)
                                     + q * <Polynomial> U.get_unsafe(ip, c))

    if transformation:
        return U
