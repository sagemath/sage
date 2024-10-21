# sage_setup: distribution = sagemath-flint
# cython: binding=True
# sage.doctest: needs sage.libs.flint sage.graphs sage.modules
r"""
Some fast computations for finite posets using FLINT matrices
"""
# ****************************************************************************
#       Copyright (C) 2020 Frédéric Chapoton <chapoton@unistra.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from cysignals.signals cimport sig_check
from cysignals.memory cimport sig_malloc, sig_free

from sage.libs.flint.fmpz cimport *
from sage.libs.flint.fmpz_mat cimport *
from sage.matrix.matrix_integer_dense cimport Matrix_integer_dense
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.integer_ring import ZZ
from sage.libs.flint.fmpz_poly_sage cimport Fmpz_poly


cpdef Fmpz_poly chain_poly(list positions):
    r"""
    Return the chain polynomial of a poset.

    INPUT:

    - ``positions`` -- a list of sets of integers describing the poset, as
      given by the lazy attribute ``_leq_storage`` of Hasse diagrams

    OUTPUT: a Flint polynomial in one variable over `\ZZ`.

    EXAMPLES::

        sage: from sage.combinat.posets.hasse_cython_flint import chain_poly
        sage: D = [{0, 1}, {1}]
        sage: chain_poly(D)
        3  1 2 1
        sage: P = posets.TamariLattice(5)
        sage: H = P._hasse_diagram
        sage: D = H._leq_storage
        sage: chain_poly(D)
        12  1 42 357 1385 3133 4635 4758 3468 1778 612 127 12
    """
    cdef Py_ssize_t n = len(positions)
    cdef Py_ssize_t i, j

    q = Fmpz_poly([0, 1])
    zero = Fmpz_poly(0)
    one = Fmpz_poly(1)

    cdef list chain_polys = [zero] * n

    # chain_polys[i] will be the generating function for the
    # chains with lowest vertex i (in the labelling of the
    # Hasse diagram).
    for i in range(n - 1, -1, -1):
        cpi = q
        for j in positions[i]:
            cpi += q * chain_polys[j]
        chain_polys[i] = cpi
    total = one
    for i in range(n):
        total += chain_polys[i]
    return total


cpdef Matrix_integer_dense moebius_matrix_fast(list positions):
    """
    Compute the Möbius matrix of a poset by a specific triangular inversion.

    INPUT:

    - ``positions`` -- a list of sets of integers describing the poset, as
      given by the lazy attribute ``_leq_storage`` of Hasse diagrams

    OUTPUT: a dense matrix

    EXAMPLES::

        sage: from sage.combinat.posets.hasse_cython_flint import moebius_matrix_fast
        sage: D = [{0,1},{1}]
        sage: moebius_matrix_fast(D)
        [ 1 -1]
        [ 0  1]
        sage: P = posets.TamariLattice(5)
        sage: H = P._hasse_diagram
        sage: D = H._leq_storage
        sage: moebius_matrix_fast(D)
        42 x 42 dense matrix over Integer Ring (...)
    """
    cdef Matrix_integer_dense A
    cdef Py_ssize_t n = len(positions)
    cdef Py_ssize_t i
    cdef int j, k
    A = Matrix_integer_dense.__new__(Matrix_integer_dense,
                                     MatrixSpace(ZZ, n, n), None, None, None)
    fmpz_mat_one(A._matrix)
    cdef Py_ssize_t *pos_lens = <Py_ssize_t*> sig_malloc(n*sizeof(Py_ssize_t))
    cdef int **pos_array = <int**> sig_malloc(n*sizeof(int*))
    cdef Py_ssize_t jind, kind
    for i in range(n):
        pos_lens[i] = len(positions[i])
        pos_array[i] = <int*> sig_malloc(pos_lens[i]*sizeof(int))
        for jind, j in enumerate(positions[i]):
            pos_array[i][jind] = j

    for i in range(n - 1, -1, -1):
        sig_check()
        for jind in range(pos_lens[i]):
            j = pos_array[i][jind]
            if j != i:
                for kind in range(pos_lens[j]):
                    k = pos_array[j][kind]
                    fmpz_sub(fmpz_mat_entry(A._matrix, i, k),
                             fmpz_mat_entry(A._matrix, i, k),
                             fmpz_mat_entry(A._matrix, j, k))
    for i in range(n):
        sig_free(pos_array[i])
    sig_free(pos_array)
    sig_free(pos_lens)
    return A


cpdef Matrix_integer_dense coxeter_matrix_fast(list positions):
    """
    Compute the Coxeter matrix of a poset by a specific algorithm.

    INPUT:

    - ``positions`` -- a list of sets of integers describing the poset, as
      given by the lazy attribute ``_leq_storage`` of Hasse diagrams

    OUTPUT: a dense matrix

    EXAMPLES::

        sage: from sage.combinat.posets.hasse_cython_flint import coxeter_matrix_fast
        sage: D = [{0,1},{1}]
        sage: coxeter_matrix_fast(D)
        [ 0 -1]
        [ 1 -1]
        sage: P = posets.TamariLattice(5)
        sage: H = P._hasse_diagram
        sage: D = H._leq_storage
        sage: coxeter_matrix_fast(D)
        42 x 42 dense matrix over Integer Ring (...)
    """
    cdef Matrix_integer_dense A
    cdef Py_ssize_t n = len(positions)
    cdef Py_ssize_t i
    cdef int j, k
    A = Matrix_integer_dense.__new__(Matrix_integer_dense,
                                     MatrixSpace(ZZ, n, n), None, None, None)
    fmpz_mat_one(A._matrix)
    cdef Py_ssize_t *pos_lens = <Py_ssize_t*> sig_malloc(n*sizeof(Py_ssize_t))
    cdef int **pos_array = <int**> sig_malloc(n*sizeof(int*))
    cdef Py_ssize_t jind, kind
    for i in range(n):
        pos_lens[i] = len(positions[i])
        pos_array[i] = <int*> sig_malloc(pos_lens[i]*sizeof(int))
        for jind, j in enumerate(positions[i]):
            pos_array[i][jind] = j

    for i in range(n - 1, -1, -1):
        sig_check()
        for jind in range(pos_lens[i]):
            j = pos_array[i][jind]
            if j != i:
                for kind in range(pos_lens[j]):
                    k = pos_array[j][kind]
                    fmpz_sub(fmpz_mat_entry(A._matrix, k, i),
                             fmpz_mat_entry(A._matrix, k, i),
                             fmpz_mat_entry(A._matrix, k, j))
    for i in range(n):
        sig_check()
        for jind in range(pos_lens[i]):
            j = pos_array[i][jind]
            if j != i:
                for k in range(n):
                    fmpz_add(fmpz_mat_entry(A._matrix, i, k),
                             fmpz_mat_entry(A._matrix, i, k),
                             fmpz_mat_entry(A._matrix, j, k))
    for i in range(n):
        sig_free(pos_array[i])
    sig_free(pos_array)
    sig_free(pos_lens)
    fmpz_mat_scalar_mul_si(A._matrix, A._matrix, -1)
    return A
