# -*- coding: utf-8 -*-
r"""
Some fast computations for finite posets
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

from sage.rings.integer_ring import ZZ
from sage.matrix.matrix_space import MatrixSpace
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet_forest


cpdef Matrix_integer_dense moebius_matrix_fast(list positions):
    """
    Compute the Möbius matrix of a poset by a specific triangular inversion.

    INPUT:

    a list of sets of integers describing the poset, as given by the
    lazy attribute ``_leq_storage`` of Hasse diagrams.

    OUTPUT:

    a dense matrix

    EXAMPLES::

        sage: from sage.combinat.posets.hasse_cython import moebius_matrix_fast
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
        for jind,j in enumerate(positions[i]):
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

    a list of sets of integers describing the poset, as given by the
    lazy attribute ``_leq_storage`` of Hasse diagrams.

    OUTPUT:

    a dense matrix

    EXAMPLES::

        sage: from sage.combinat.posets.hasse_cython import coxeter_matrix_fast
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
        for jind,j in enumerate(positions[i]):
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


class IncreasingChains(RecursivelyEnumeratedSet_forest):
    r"""
    The enumerated set of increasing chains.

    INPUT:

    - ``positions`` -- a list of sets of integers describing the poset,
      as given by the lazy attribute ``_leq_storage`` of Hasse diagrams.

    - ``element_constructor`` -- used to determine the type of chains,
      for example ``list`` or ``tuple``

    - ``exclude`` -- list of integers that should not belong to the chains

    - ``conversion`` -- optional list of elements of the poset

    If ``conversion`` is provided, it is used to convert chain elements
    to elements of this list.

    EXAMPLES::

        sage: from sage.combinat.posets.hasse_cython import IncreasingChains
        sage: D = IncreasingChains([{0,1},{1}], list, []); D
        An enumerated set with a forest structure
        sage: D.cardinality()
        4
        sage: list(D)
        [[], [0], [0, 1], [1]]
    """
    def __init__(self, list positions, element_constructor,
                 list exclude, conversion=None):
        """
        The enumerated set of increasing chains.

        TESTS::

            sage: from sage.combinat.posets.hasse_cython import IncreasingChains
            sage: D = IncreasingChains([{0,1},{1}], list, [])
            sage: TestSuite(D).run(skip='_test_pickling')
        """
        cdef int i
        cdef Py_ssize_t n = len(positions)
        self._n = n
        if exclude is not None:
            self._greater_than = [gti.difference(exclude) for gti in positions]
            self._vertices = [u for u in range(n) if u not in exclude]
        else:
            self._greater_than = positions
            self._vertices = list(range(n))

        self._element_constructor = element_constructor
        self._conversion = conversion
        if conversion is not None:
            self._from_poset = {elt: i for i, elt in enumerate(conversion)}

        self._roots = (tuple(),)

        RecursivelyEnumeratedSet_forest.__init__(self, algorithm='depth',
                    category=FiniteEnumeratedSets())

    def __contains__(self, tup):
        """
        Membership testing.

        If ``conversion`` was provided, it first converts elements of ``tup``
        to integers.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_cython import IncreasingChains
            sage: D = IncreasingChains([{0,1},{1}], list, [])
            sage: [x in D for x in D]
            [True, True, True, True]
            sage: [2] in D
            False
            sage: [1,1] in D
            False

            sage: P = Poset({'a':['b'],'b':[]})
            sage: ['a'] in P.chains()
            True
        """
        cdef int k
        cdef Py_ssize_t i, x, y
        if not tup:
            return True
        if self._conversion is not None:
            tup = [self._from_poset[elt] for elt in tup]
        for i in tup:
            if not(0 <= i < self._n):
                return False
        y = tup[0]
        for k in range(1, len(tup)):
            x = y
            y = tup[k]
            if x == y or y not in self._greater_than[x]:
                return False
        return True

    def post_process(self, chain):
        """
        Create a chain from the internal object.

        If ``conversion`` was provided, it first converts elements of the
        chain to elements of this list.

        Then the given ``element_constructor`` is applied to the chain.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_cython import IncreasingChains
            sage: D = IncreasingChains([{0,1},{1}], list, [])
            sage: D.post_process((0,1))
            [0, 1]

            sage: P = Poset({'a':['b'],'b':[]})
            sage: list(P.chains())
            [[], ['a'], ['a', 'b'], ['b']]
        """
        cdef Py_ssize_t i
        if self._conversion is not None:
            return self._element_constructor(self._conversion[i] for i in chain)
        return self._element_constructor(chain)

    def children(self, chain):
        """
        Return the children of a chain, by adding one largest element.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_cython import IncreasingChains
            sage: D = IncreasingChains([{0,1},{1}], list, [])
            sage: D.children((0,))
            [(0, 1)]

            sage: P = Poset({'a':['b'],'b':[]})
            sage: next(iter(P.chains()))
            []
        """
        cdef Py_ssize_t x, y
        if not chain:
            return [(x,) for x in self._vertices]
        else:
            x = chain[-1]
            return [chain + (y,) for y in self._greater_than[x] if x != y]
