# sage_setup: distribution = sagemath-singular
# distutils: language = c++
# distutils: libraries = gmp m NTL_LIBRARIES
# distutils: extra_compile_args = NTL_CFLAGS
# distutils: include_dirs = NTL_INCDIR
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
"""
Optimized Cython code needed by quaternion algebras

This is a collection of miscellaneous routines that are in Cython for
speed purposes and are used by the quaternion algebra code.  For
example, there are functions for quickly constructing an n x 4 matrix
from a list of n rational quaternions.

AUTHORS:

- William Stein
"""

# ****************************************************************************
#       Copyright (C) 2009 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.integer cimport Integer
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.matrix_integer_dense cimport Matrix_integer_dense
from sage.matrix.matrix_rational_dense cimport Matrix_rational_dense

from sage.algebras.quatalg.quaternion_algebra_element cimport QuaternionAlgebraElement_rational_field

from sage.libs.gmp.mpz cimport mpz_t, mpz_lcm, mpz_init, mpz_set, mpz_clear, mpz_init_set, mpz_mul, mpz_fdiv_q, mpz_cmp_si

from sage.libs.flint.fmpz cimport fmpz_set_mpz
from sage.libs.flint.fmpq cimport fmpq_canonicalise
from sage.libs.flint.fmpq_mat cimport fmpq_mat_entry_num, fmpq_mat_entry_den, fmpq_mat_entry


def integral_matrix_and_denom_from_rational_quaternions(v, reverse=False):
    r"""
    Given a list of rational quaternions, return matrix `A` over `\ZZ`
    and denominator `d`, such that the rows of `(1/d)A` are the
    entries of the quaternions.

    INPUT:

    - ``v`` -- list of quaternions in a rational quaternion algebra
    - ``reverse`` -- whether order of the coordinates as well as the
      order of the list ``v`` should be reversed

    OUTPUT:

    - a matrix over `\ZZ`
    - an integer (the common denominator)

    EXAMPLES::

        sage: A.<i,j,k>=QuaternionAlgebra(-4,-5)
        sage: sage.algebras.quatalg.quaternion_algebra_cython.integral_matrix_and_denom_from_rational_quaternions([i/2,1/3+j+k])
        (
        [0 3 0 0]
        [2 0 6 6], 6
        )

        sage: sage.algebras.quatalg.quaternion_algebra_cython.integral_matrix_and_denom_from_rational_quaternions([i/2,1/3+j+k], reverse=True)
        (
        [6 6 0 2]
        [0 0 3 0], 6
        )
    """
    # This function is an optimized version of
    #   MatrixSpace(QQ,len(v),4)([x.coefficient_tuple() for x in v], coerce=False)._clear_denom

    cdef Py_ssize_t i, n=len(v)
    M = MatrixSpace(ZZ, n, 4)
    cdef Matrix_integer_dense A = M.zero_matrix().__copy__()
    if n == 0:
        return A

    # Find least common multiple of the denominators
    cdef QuaternionAlgebraElement_rational_field x
    cdef Integer d = Integer()
    # set denom to the denom of the first quaternion
    x = v[0]
    mpz_set(d.value, x.d)
    for x in v[1:]:
        mpz_lcm(d.value, d.value, x.d)

    # Now fill in each row x of A, multiplying it by q = d/denom(x)
    cdef mpz_t q
    cdef mpz_t tmp
    mpz_init(q)
    mpz_init(tmp)
    for i in range(n):
        x = v[i]
        mpz_fdiv_q(q, d.value, x.d)
        if reverse:
            mpz_mul(tmp, q, x.x)
            A.set_unsafe_mpz(n-i-1, 3, tmp)
            mpz_mul(tmp, q, x.y)
            A.set_unsafe_mpz(n-i-1, 2, tmp)
            mpz_mul(tmp, q, x.z)
            A.set_unsafe_mpz(n-i-1, 1, tmp)
            mpz_mul(tmp, q, x.w)
            A.set_unsafe_mpz(n-i-1, 0, tmp)
        else:
            mpz_mul(tmp, q, x.x)
            A.set_unsafe_mpz(i, 0, tmp)
            mpz_mul(tmp, q, x.y)
            A.set_unsafe_mpz(i, 1, tmp)
            mpz_mul(tmp, q, x.z)
            A.set_unsafe_mpz(i, 2, tmp)
            mpz_mul(tmp, q, x.w)
            A.set_unsafe_mpz(i, 3, tmp)
    mpz_clear(q)
    mpz_clear(tmp)
    return A, d


def rational_matrix_from_rational_quaternions(v, reverse=False):
    r"""
    Return matrix over the rationals whose rows have entries the
    coefficients of the rational quaternions in ``v``.

    INPUT:

    - ``v`` -- list of quaternions in a rational quaternion algebra
    - ``reverse`` -- whether order of the coordinates as well as the
      order of the list ``v`` should be reversed

    OUTPUT: a matrix over `\QQ`

    EXAMPLES::

        sage: A.<i,j,k>=QuaternionAlgebra(-4,-5)
        sage: sage.algebras.quatalg.quaternion_algebra_cython.rational_matrix_from_rational_quaternions([i/2,1/3+j+k])
        [  0 1/2   0   0]
        [1/3   0   1   1]

        sage: sage.algebras.quatalg.quaternion_algebra_cython.rational_matrix_from_rational_quaternions([i/2,1/3+j+k], reverse=True)
        [  1   1   0 1/3]
        [  0   0 1/2   0]
    """
    cdef Py_ssize_t i, j, n=len(v)
    M = MatrixSpace(QQ, n, 4)
    cdef Matrix_rational_dense A = M.zero_matrix().__copy__()
    if n == 0:
        return A

    cdef QuaternionAlgebraElement_rational_field x
    if reverse:
        for i in range(n):
            x = v[i]
            fmpz_set_mpz(fmpq_mat_entry_num(A._matrix, n-i-1, 3), x.x)
            fmpz_set_mpz(fmpq_mat_entry_num(A._matrix, n-i-1, 2), x.y)
            fmpz_set_mpz(fmpq_mat_entry_num(A._matrix, n-i-1, 1), x.z)
            fmpz_set_mpz(fmpq_mat_entry_num(A._matrix, n-i-1, 0), x.w)

            if mpz_cmp_si(x.d, 1):
                for j in range(4):
                    fmpz_set_mpz(fmpq_mat_entry_den(A._matrix, n-i-1, j), x.d)
                    fmpq_canonicalise(fmpq_mat_entry(A._matrix, n-i-1, j))
    else:
        for i in range(n):
            x = v[i]
            fmpz_set_mpz(fmpq_mat_entry_num(A._matrix, i, 0), x.x)
            fmpz_set_mpz(fmpq_mat_entry_num(A._matrix, i, 1), x.y)
            fmpz_set_mpz(fmpq_mat_entry_num(A._matrix, i, 2), x.z)
            fmpz_set_mpz(fmpq_mat_entry_num(A._matrix, i, 3), x.w)

            if mpz_cmp_si(x.d, 1):
                for j in range(4):
                    fmpz_set_mpz(fmpq_mat_entry_den(A._matrix, i, j), x.d)
                    fmpq_canonicalise(fmpq_mat_entry(A._matrix, i, j))

    return A


def rational_quaternions_from_integral_matrix_and_denom(A, Matrix_integer_dense H, Integer d, reverse=False):
    r"""
    Given an integral matrix and denominator, returns a list of
    rational quaternions.

    INPUT:

    - ``A`` -- rational quaternion algebra
    - ``H`` -- matrix over the integers
    - ``d`` -- integer
    - ``reverse`` -- whether order of the coordinates as well as the
      order of the list ``v`` should be reversed

    OUTPUT:

    - list of ``H.nrows()`` elements of ``A``

    EXAMPLES::

        sage: A.<i,j,k>=QuaternionAlgebra(-1,-2)
        sage: f = sage.algebras.quatalg.quaternion_algebra_cython.rational_quaternions_from_integral_matrix_and_denom
        sage: f(A, matrix([[1,2,3,4],[-1,2,-4,3]]), 3)
        [1/3 + 2/3*i + j + 4/3*k, -1/3 + 2/3*i - 4/3*j + k]

        sage: f(A, matrix([[3,-4,2,-1],[4,3,2,1]]), 3, reverse=True)
        [1/3 + 2/3*i + j + 4/3*k, -1/3 + 2/3*i - 4/3*j + k]
    """
    #
    # This is an optimized version of the following interpreted Python code.
    # H2 = H.change_ring(QQ)._rmul_(1/d)
    # return [A(v.list()) for v in H2.rows()]
    #
    cdef QuaternionAlgebraElement_rational_field x
    v = []
    cdef Integer a, b
    a = Integer(A.invariants()[0])
    b = Integer(A.invariants()[1])
    cdef Py_ssize_t i
    cdef mpz_t tmp
    mpz_init(tmp)

    if reverse:
        rng = range(H.nrows()-1, -1, -1)
    else:
        rng = range(H.nrows())

    for i in rng:
        x = <QuaternionAlgebraElement_rational_field> QuaternionAlgebraElement_rational_field.__new__(QuaternionAlgebraElement_rational_field)
        x._parent = A
        mpz_set(x.a, a.value)
        mpz_set(x.b, b.value)
        if reverse:
            H.get_unsafe_mpz(i, 3, tmp)
            mpz_init_set(x.x, tmp)
            H.get_unsafe_mpz(i, 2, tmp)
            mpz_init_set(x.y, tmp)
            H.get_unsafe_mpz(i, 1, tmp)
            mpz_init_set(x.z, tmp)
            H.get_unsafe_mpz(i, 0, tmp)
            mpz_init_set(x.w, tmp)
        else:
            H.get_unsafe_mpz(i, 0, tmp)
            mpz_init_set(x.x, tmp)
            H.get_unsafe_mpz(i, 1, tmp)
            mpz_init_set(x.y, tmp)
            H.get_unsafe_mpz(i, 2, tmp)
            mpz_init_set(x.z, tmp)
            H.get_unsafe_mpz(i, 3, tmp)
            mpz_init_set(x.w, tmp)
        mpz_init_set(x.d, d.value)
        # WARNING -- we do *not* canonicalize the entries in the quaternion.
        # This is I think _not_ needed for quaternion_element.pyx
        v.append(x)
    mpz_clear(tmp)
    return v
