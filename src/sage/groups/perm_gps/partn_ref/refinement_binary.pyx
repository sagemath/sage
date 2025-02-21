# sage.doctest: optional - sage.modules sage.rings.finite_rings
"""
Partition backtrack functions for binary codes

EXAMPLES::

    sage: import sage.groups.perm_gps.partn_ref.refinement_binary

REFERENCE:

- [1] McKay, Brendan D. Practical Graph Isomorphism. Congressus Numerantium,
  Vol. 30 (1981), pp. 45-87.

- [2] Leon, Jeffrey. Permutation Group Algorithms Based on Partitions, I:
  Theory and Algorithms. J. Symbolic Computation, Vol. 12 (1991), pp.
  533-583.
"""

#*****************************************************************************
#       Copyright (C) 2006 - 2011 Robert L. Miller <rlmillster@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from sage.data_structures.bitset_base cimport *
from sage.groups.perm_gps.partn_ref.data_structures cimport *
from sage.rings.integer cimport Integer
from sage.structure.element import Matrix
from sage.groups.perm_gps.partn_ref.double_coset cimport double_coset


cdef class LinearBinaryCodeStruct(BinaryCodeStruct):

    def __cinit__(self, matrix):
        cdef int i,j
        self.degree = matrix.ncols()
        self.dimension = matrix.nrows()
        if self.dimension >= <int>(sizeof(int) << 3):
            raise NotImplementedError
            # By the time the dimension gets this big, the computation is infeasible anyway...
        self.nwords = 1<<self.dimension

        self.basis = <bitset_s *> sig_malloc(self.dimension * sizeof(bitset_s))
        self.scratch_bitsets = <bitset_s *> sig_malloc((2*self.dimension+2) * sizeof(bitset_s))
        self.alpha_is_wd = <bitset_s *> sig_malloc(sizeof(bitset_s))
        self.word_ps = PS_new(self.nwords, 1)
        self.alpha = <int *> sig_malloc((self.nwords+self.degree) * sizeof(int))
        self.scratch = <int *> sig_malloc((3*self.nwords+3*self.degree+2) * sizeof(int))

        if self.basis       is NULL or self.scratch_bitsets is NULL \
        or self.alpha_is_wd is NULL or self.word_ps         is NULL \
        or self.alpha       is NULL or self.scratch         is NULL:
            sig_free(self.basis)
            sig_free(self.scratch_bitsets)
            sig_free(self.alpha_is_wd)
            PS_dealloc(self.word_ps)
            sig_free(self.alpha)
            sig_free(self.scratch)
            raise MemoryError

        cdef bint memerr = 0
        for i in range(self.dimension):
            try:
                bitset_init(&self.basis[i], self.degree)
            except MemoryError:
                for j from 0 <= j < i:
                    bitset_free(&self.basis[j])
                memerr = 1
                break
        if not memerr:
            for i in range(2*self.dimension+2):
                try:
                    bitset_init(&self.scratch_bitsets[i], self.degree)
                except MemoryError:
                    for j from 0 <= j < i:
                        bitset_free(&self.scratch_bitsets[j])
                    for j from 0 <= j < self.dimension:
                        bitset_free(&self.basis[j])
                    memerr = 1
                    break
        if not memerr:
            try:
                bitset_init(self.alpha_is_wd, self.nwords + self.degree)
            except MemoryError:
                for j from 0 <= j < 2*self.dimension+2:
                    bitset_free(&self.scratch_bitsets[j])
                for j from 0 <= j < self.dimension:
                    bitset_free(&self.basis[j])
                memerr = 1
        if memerr:
            sig_free(self.basis)
            sig_free(self.scratch_bitsets)
            sig_free(self.alpha_is_wd)
            PS_dealloc(self.word_ps)
            sig_free(self.alpha)
            sig_free(self.scratch)
            raise MemoryError
        else:
            bitset_zero(self.alpha_is_wd)
            for j from 0 <= j < self.dimension:
                bitset_zero(&self.basis[j])

        for i, j in matrix.nonzero_positions():
            bitset_set(&self.basis[i], j)

        self.output = NULL
        self.ith_word = &ith_word_linear

    def run(self, partition=None):
        """
        Perform the canonical labeling and automorphism group computation,
        storing results to ``self``.

        INPUT:

        - ``partition`` -- an optional list of lists partition of the columns;
          default is the unit partition

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref.refinement_binary import LinearBinaryCodeStruct

            sage: B = LinearBinaryCodeStruct(matrix(GF(2),[[1,0,1],[0,1,1]]))
            sage: B.run()
            sage: B.automorphism_group()
            ([[0, 2, 1], [1, 0, 2]], 6, [0, 1])
            sage: B.canonical_relabeling()
            [0, 1, 2]

            sage: B = LinearBinaryCodeStruct(matrix(GF(2),[[1,1,1,1]]))
            sage: B.automorphism_group()
            ([[0, 1, 3, 2], [0, 2, 1, 3], [1, 0, 2, 3]], 24, [0, 1, 2])
            sage: B.canonical_relabeling()
            [0, 1, 2, 3]

            sage: B = LinearBinaryCodeStruct(matrix(GF(2),[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]]))
            sage: B.automorphism_group()[1] == factorial(14)
            True

            sage: M = Matrix(GF(2),[\
            ....: [1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0],\
            ....: [0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0],\
            ....: [0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1],\
            ....: [0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1],\
            ....: [0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1]])
            sage: B = LinearBinaryCodeStruct(M)
            sage: B.automorphism_group()[1]
            322560

            sage: M = Matrix(GF(2),[\
            ....: [1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0],\
            ....: [0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0],\
            ....: [0,0,0,0,0,1,0,1,0,0,0,1,1,1,1,1,1],\
            ....: [0,0,0,1,1,0,0,0,0,1,1,0,1,1,0,1,1]])
            sage: B = LinearBinaryCodeStruct(M)
            sage: B.automorphism_group()[1]
            2304

            sage: M = Matrix(GF(2),[\
            ....: [1,0,0,1,1,1,1,0,0,1,0,0,0,0,0,0,0],\
            ....: [0,1,0,0,1,1,1,1,0,0,1,0,0,0,0,0,0],\
            ....: [0,0,1,0,0,1,1,1,1,0,0,1,0,0,0,0,0],\
            ....: [0,0,0,1,0,0,1,1,1,1,0,0,1,0,0,0,0],\
            ....: [0,0,0,0,1,0,0,1,1,1,1,0,0,1,0,0,0],\
            ....: [0,0,0,0,0,1,0,0,1,1,1,1,0,0,1,0,0],\
            ....: [0,0,0,0,0,0,1,0,0,1,1,1,1,0,0,1,0],\
            ....: [0,0,0,0,0,0,0,1,0,0,1,1,1,1,0,0,1]])
            sage: B = LinearBinaryCodeStruct(M)
            sage: B.automorphism_group()[1]
            136

            sage: M = Matrix(GF(2),[\
            ....: [1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0],
            ....: [0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0],
            ....: [0,0,0,0,1,1,0,0,0,0,0,0,1,1,1,1,1,1],
            ....: [0,0,1,1,0,0,0,0,0,0,1,1,1,1,0,0,1,1],
            ....: [0,0,0,1,0,0,0,1,0,1,0,1,0,1,1,1,0,1],
            ....: [0,1,0,0,0,1,0,0,0,1,1,1,0,1,0,1,1,0]])
            sage: B = LinearBinaryCodeStruct(M)
            sage: B.automorphism_group()[1]
            2160

            sage: M = Matrix(GF(2),[\
            ....: [0,1,0,1,1,1,0,0,0,1,0,0,0,1,0,0,0,1,1,1,0,1],\
            ....: [1,0,1,1,1,0,0,0,1,0,0,0,1,0,0,0,1,1,1,0,1,0],\
            ....: [0,1,1,1,0,0,0,1,0,0,1,1,0,0,0,1,1,1,0,1,0,0],\
            ....: [1,1,1,0,0,0,1,0,0,1,0,0,0,0,1,1,1,0,1,0,0,1],\
            ....: [1,1,0,0,0,1,0,0,1,0,1,0,0,1,1,1,0,1,0,0,1,0],\
            ....: [1,0,0,0,1,0,0,1,0,1,1,0,1,1,1,0,1,0,0,1,0,0],\
            ....: [0,0,0,1,0,0,1,0,1,1,1,1,1,1,0,1,0,0,1,0,0,0],\
            ....: [0,0,1,0,0,1,0,1,1,1,0,1,1,0,1,0,0,1,0,0,0,1],\
            ....: [0,1,0,0,1,0,1,1,1,0,0,1,0,1,0,0,1,0,0,0,1,1],\
            ....: [1,0,0,1,0,1,1,1,0,0,0,0,1,0,0,1,0,0,0,1,1,1],\
            ....: [0,0,1,0,1,1,1,0,0,0,1,1,0,0,1,0,0,0,1,1,1,0]])
            sage: B = LinearBinaryCodeStruct(M)
            sage: B.automorphism_group()[1]
            887040

            sage: M = Matrix(GF(2),[\
            ....: [1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
            ....: [0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
            ....: [0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0],
            ....: [0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0],
            ....: [0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0],
            ....: [0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0],
            ....: [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1],
            ....: [1,0,1,0,1,0,1,0,1,1,0,0,0,0,0,0,1,1,0,0],
            ....: [1,1,0,0,0,0,0,0,1,0,1,0,1,0,1,0,1,1,0,0],
            ....: [1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,1,0]])
            sage: B = LinearBinaryCodeStruct(M)
            sage: B.automorphism_group()[1]
            294912

            sage: M = Matrix(GF(2), [\
            ....: [1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
            ....: [0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
            ....: [0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0],
            ....: [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0],
            ....: [0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0],
            ....: [0,0,0,0,0,0,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0],
            ....: [0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,1,0,1,0,1,0,1,1],
            ....: [0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,1,0,0,1,1,1,0,1],
            ....: [0,0,0,0,0,1,0,1,0,0,0,1,0,0,0,1,1,1,1,0,0,0,1]])
            sage: B = LinearBinaryCodeStruct(M)
            sage: B.automorphism_group()[1]
            442368

            sage: M = Matrix(GF(2), [\
            ....: [1,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0],\
            ....: [1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,1]])
            sage: B = LinearBinaryCodeStruct(M)
            sage: B.automorphism_group()[1]
            17868913969917295853568000000
        """
        cdef int n = self.degree
        cdef PartitionStack *part
        if partition is None:
            part = PS_new(n, 1)
        else:
            part = PS_from_list(partition)
        if part is NULL:
            raise MemoryError

        self.first_time = 1

        self.output = get_aut_gp_and_can_lab(<void *>self, part, n, &all_children_are_equivalent, &refine_by_bip_degree, &compare_linear_codes, 1, NULL, NULL, NULL)

        PS_dealloc(part)

    def automorphism_group(self):
        """
        Return a list of generators of the automorphism group, along with its
        order and a base for which the list of generators is a strong generating
        set.

        EXAMPLES: (For more examples, see self.run())::

            sage: from sage.groups.perm_gps.partn_ref.refinement_binary import LinearBinaryCodeStruct

            sage: B = LinearBinaryCodeStruct(matrix(GF(2),[[1,1,1,1]]))
            sage: B.automorphism_group()
            ([[0, 1, 3, 2], [0, 2, 1, 3], [1, 0, 2, 3]], 24, [0, 1, 2])
        """
        cdef int i, j
        cdef list generators, base
        cdef Integer order
        if self.output is NULL:
            self.run()
        generators = []
        for i in range(self.output.num_gens):
            generators.append([self.output.generators[i*self.degree + j] for j from 0 <= j < self.degree])
        order = Integer()
        SC_order(self.output.group, 0, order.value)
        base = [self.output.group.base_orbits[i][0] for i in range(self.output.group.base_size)]
        return generators, order, base

    def canonical_relabeling(self):
        """
        Return a canonical relabeling (in list permutation format).

        EXAMPLES: (For more examples, see self.run())::

            sage: from sage.groups.perm_gps.partn_ref.refinement_binary import LinearBinaryCodeStruct

            sage: B = LinearBinaryCodeStruct(matrix(GF(2), [[1,1,0]]))
            sage: B.automorphism_group()
            ([[1, 0, 2]], 2, [0])
            sage: B.canonical_relabeling()
            [1, 2, 0]
            sage: B = LinearBinaryCodeStruct(matrix(GF(2), [[1,0,1]]))
            sage: B.automorphism_group()
            ([[2, 1, 0]], 2, [0])
            sage: B.canonical_relabeling()
            [1, 0, 2]
        """
        cdef int i
        if self.output is NULL:
            self.run()
        return [self.output.relabeling[i] for i in range(self.degree)]

    def is_isomorphic(self, LinearBinaryCodeStruct other):
        """
        Calculate whether ``self`` is isomorphic to ``other``.

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref.refinement_binary import LinearBinaryCodeStruct

            sage: B = LinearBinaryCodeStruct(Matrix(GF(2), [[1,1,1,1,0,0],[0,0,1,1,1,1]]))
            sage: C = LinearBinaryCodeStruct(Matrix(GF(2), [[1,1,1,0,0,1],[1,1,0,1,1,0]]))
            sage: B.is_isomorphic(C)
            [0, 1, 2, 5, 3, 4]
        """
        cdef int i, n = self.degree
        cdef int *output
        cdef int *ordering
        cdef PartitionStack *part
        part = PS_new(n, 1)
        ordering = <int *> sig_malloc(self.degree * sizeof(int))
        output = <int *> sig_malloc(self.degree * sizeof(int))
        if part is NULL or ordering is NULL or output is NULL:
            PS_dealloc(part)
            sig_free(ordering)
            sig_free(output)
            raise MemoryError
        for i in range(n):
            ordering[i] = i
        self.first_time = 1
        other.first_time = 1

        cdef bint isomorphic = double_coset(<void *> self, <void *> other, part, ordering, n, &all_children_are_equivalent, &refine_by_bip_degree, &compare_linear_codes, NULL, NULL, output)

        PS_dealloc(part)
        sig_free(ordering)
        if isomorphic:
            output_py = [output[i] for i in range(n)]
        else:
            output_py = False
        sig_free(output)
        return output_py

    def __dealloc__(self):
        cdef int j
        bitset_free(self.alpha_is_wd)
        for j from 0 <= j < 2*self.dimension+2:
            bitset_free(&self.scratch_bitsets[j])
        for j from 0 <= j < self.dimension:
            bitset_free(&self.basis[j])
        sig_free(self.basis)
        sig_free(self.scratch_bitsets)
        sig_free(self.alpha_is_wd)
        PS_dealloc(self.word_ps)
        sig_free(self.alpha)
        sig_free(self.scratch)
        if self.output is not NULL:
            deallocate_agcl_output(self.output)

cdef int ith_word_linear(BinaryCodeStruct self, int i, bitset_s *word) noexcept:
    cdef LinearBinaryCodeStruct LBCS = <LinearBinaryCodeStruct> self
    cdef int j
    bitset_zero(word)
    for j from 0 <= j < LBCS.dimension:
        if (1<<j)&i:
            bitset_xor(word, word, &LBCS.basis[j])
    return 0

cdef class NonlinearBinaryCodeStruct(BinaryCodeStruct):

    def __cinit__(self, arg):
        cdef int i,j
        if isinstance(arg, Matrix):
            self.degree = arg.ncols()
            self.nwords = arg.nrows()
        elif isinstance(arg, tuple):
            assert len(arg) == 2
            self.degree, self.nwords = arg
        else:
            raise NotImplementedError

        self.words = <bitset_s *> sig_malloc(self.nwords * sizeof(bitset_s))
        self.scratch_bitsets = <bitset_s *> sig_malloc((4*self.nwords+1) * sizeof(bitset_s))
        self.alpha_is_wd = <bitset_s *> sig_malloc(sizeof(bitset_s))
        self.word_ps = PS_new(self.nwords, 1)
        self.alpha = <int *> sig_malloc((self.nwords+self.degree) * sizeof(int))
        self.scratch = <int *> sig_malloc((3*self.nwords+3*self.degree+2) * sizeof(int))
        if self.words       is NULL or self.scratch_bitsets is NULL \
        or self.alpha_is_wd is NULL or self.word_ps         is NULL \
        or self.alpha       is NULL or self.scratch         is NULL:
            sig_free(self.words)
            sig_free(self.scratch_bitsets)
            sig_free(self.alpha_is_wd)
            PS_dealloc(self.word_ps)
            sig_free(self.alpha)
            sig_free(self.scratch)
            raise MemoryError

        cdef bint memerr = 0
        for i in range(self.nwords):
            try:
                bitset_init(&self.words[i], self.degree)
            except MemoryError:
                for j from 0 <= j < i:
                    bitset_free(&self.words[j])
                memerr = 1
                break
        if not memerr:
            for i in range(4*self.nwords):
                try:
                    bitset_init(&self.scratch_bitsets[i], self.degree)
                except MemoryError:
                    for j from 0 <= j < i:
                        bitset_free(&self.scratch_bitsets[j])
                    for j from 0 <= j < self.nwords:
                        bitset_free(&self.words[j])
                    memerr = 1
                    break
        if not memerr:
            try:
                bitset_init(&self.scratch_bitsets[4*self.nwords], self.nwords)
            except MemoryError:
                for j from 0 <= j < 4*self.nwords:
                    bitset_free(&self.scratch_bitsets[j])
                for j from 0 <= j < self.nwords:
                    bitset_free(&self.words[j])
                memerr = 1
        if not memerr:
            try:
                bitset_init(self.alpha_is_wd, self.nwords + self.degree)
            except MemoryError:
                for j from 0 <= j < 4*self.nwords + 1:
                    bitset_free(&self.scratch_bitsets[j])
                for j from 0 <= j < self.nwords:
                    bitset_free(&self.words[j])
                memerr = 1
        if memerr:
            sig_free(self.words)
            sig_free(self.scratch_bitsets)
            sig_free(self.alpha_is_wd)
            PS_dealloc(self.word_ps)
            sig_free(self.alpha)
            sig_free(self.scratch)
            raise MemoryError
        else:
            bitset_zero(self.alpha_is_wd)
            for j from 0 <= j < self.nwords:
                bitset_zero(&self.words[j])

        if isinstance(arg, Matrix):
            for i, j in arg.nonzero_positions():
                bitset_set(&self.words[i], j)

        self.output = NULL
        self.ith_word = &ith_word_nonlinear

    def __dealloc__(self):
        cdef int j
        bitset_free(self.alpha_is_wd)
        for j from 0 <= j < 4*self.nwords + 1:
            bitset_free(&self.scratch_bitsets[j])
        for j from 0 <= j < self.nwords:
            bitset_free(&self.words[j])
        sig_free(self.words)
        sig_free(self.scratch_bitsets)
        sig_free(self.alpha_is_wd)
        PS_dealloc(self.word_ps)
        sig_free(self.alpha)
        sig_free(self.scratch)
        if self.output is not NULL:
            deallocate_agcl_output(self.output)

    def run(self, partition=None):
        """
        Perform the canonical labeling and automorphism group computation,
        storing results to ``self``.

        INPUT:

        - ``partition`` -- an optional list of lists partition of the columns
          default is the unit partition

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref.refinement_binary import NonlinearBinaryCodeStruct

            sage: B = NonlinearBinaryCodeStruct(Matrix(GF(2), [[1,0,0,0],[0,0,1,0]]))
            sage: B.run()
            sage: B.automorphism_group()
            ([[2, 1, 0, 3], [0, 3, 2, 1]], 4, [1, 0])
            sage: B.canonical_relabeling()
            [2, 0, 3, 1]

            sage: B = NonlinearBinaryCodeStruct(Matrix(GF(2), [[1,1,1,0],[1,1,0,1],[1,0,1,1],[0,1,1,1]]))
            sage: B.run()
            sage: B.automorphism_group()
            ([[0, 1, 3, 2], [0, 2, 1, 3], [1, 0, 2, 3]], 24, [0, 1, 2])
            sage: B.canonical_relabeling()
            [0, 1, 2, 3]

            sage: B = NonlinearBinaryCodeStruct(Matrix(GF(2), [[1,1,1,0,0,0],[1,1,0,1,0,0],[1,0,1,1,0,0],[0,1,1,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]]))
            sage: B.run()
            sage: B.automorphism_group()
            ([[0, 1, 3, 2, 4, 5],
              [0, 2, 1, 3, 4, 5],
              [1, 0, 2, 3, 4, 5],
              [0, 1, 2, 3, 5, 4]],
             48,
             [4, 0, 1, 2])
            sage: B.canonical_relabeling()
            [2, 3, 4, 5, 0, 1]
        """
        cdef int n = self.degree
        cdef PartitionStack *part
        if partition is None:
            part = PS_new(n, 1)
        else:
            part = PS_from_list(partition)
        if part is NULL:
            raise MemoryError
        self.first_time = 1

        self.output = get_aut_gp_and_can_lab(<void *> self, part, self.degree, &all_children_are_equivalent, &refine_by_bip_degree, &compare_nonlinear_codes, 1, NULL, NULL, NULL)

        PS_dealloc(part)

    def automorphism_group(self):
        """
        Return a list of generators of the automorphism group, along with its
        order and a base for which the list of generators is a strong generating
        set.

        EXAMPLES: (For more examples, see self.run())::

            sage: from sage.groups.perm_gps.partn_ref.refinement_binary import NonlinearBinaryCodeStruct

            sage: B = NonlinearBinaryCodeStruct(Matrix(GF(2), [[1,1,1,0,0,0],[1,1,0,1,0,0],[1,0,1,1,0,0],[0,1,1,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]]))
            sage: B.run()
            sage: B.automorphism_group()
            ([[0, 1, 3, 2, 4, 5],
              [0, 2, 1, 3, 4, 5],
              [1, 0, 2, 3, 4, 5],
              [0, 1, 2, 3, 5, 4]],
             48,
             [4, 0, 1, 2])
        """
        cdef int i, j
        cdef list generators, base
        cdef Integer order
        if self.output is NULL:
            self.run()
        generators = []
        for i in range(self.output.num_gens):
            generators.append([self.output.generators[i*self.degree + j] for j from 0 <= j < self.degree])
        order = Integer()
        SC_order(self.output.group, 0, order.value)
        base = [self.output.group.base_orbits[i][0] for i in range(self.output.group.base_size)]
        return generators, order, base

    def canonical_relabeling(self):
        """
        Return a canonical relabeling (in list permutation format).

        EXAMPLES: (For more examples, see self.run())::

            sage: from sage.groups.perm_gps.partn_ref.refinement_binary import NonlinearBinaryCodeStruct

            sage: B = NonlinearBinaryCodeStruct(Matrix(GF(2), [[1,1,1,0,0,0],[1,1,0,1,0,0],[1,0,1,1,0,0],[0,1,1,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]]))
            sage: B.run()
            sage: B.canonical_relabeling()
            [2, 3, 4, 5, 0, 1]
        """
        cdef int i
        if self.output is NULL:
            self.run()
        return [self.output.relabeling[i] for i in range(self.degree)]

    def is_isomorphic(self, NonlinearBinaryCodeStruct other):
        """
        Calculate whether ``self`` is isomorphic to ``other``.

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref.refinement_binary import NonlinearBinaryCodeStruct

            sage: B = NonlinearBinaryCodeStruct(Matrix(GF(2), [[1,1,1,1,0,0],[0,0,1,1,1,1]]))
            sage: C = NonlinearBinaryCodeStruct(Matrix(GF(2), [[1,1,0,0,1,1],[1,1,1,1,0,0]]))
            sage: B.is_isomorphic(C)
            [2, 3, 0, 1, 4, 5]
        """
        cdef int i, n = self.degree
        cdef int *output
        cdef int *ordering
        cdef PartitionStack *part
        part = PS_new(n, 1)
        ordering = <int *> sig_malloc(n * sizeof(int))
        output = <int *> sig_malloc(n * sizeof(int))
        if part is NULL or ordering is NULL or output is NULL:
            PS_dealloc(part)
            sig_free(ordering)
            sig_free(output)
            raise MemoryError
        for i in range(n):
            ordering[i] = i
        self.first_time = 1
        other.first_time = 1

        cdef bint isomorphic = double_coset(<void *> self, <void *> other, part, ordering, n, &all_children_are_equivalent, &refine_by_bip_degree, &compare_nonlinear_codes, NULL, NULL, output)

        PS_dealloc(part)
        sig_free(ordering)
        if isomorphic:
            output_py = [output[i] for i in range(n)]
        else:
            output_py = False
        sig_free(output)
        return output_py

cdef int ith_word_nonlinear(BinaryCodeStruct self, int i, bitset_s *word) noexcept:
    cdef NonlinearBinaryCodeStruct NBCS = <NonlinearBinaryCodeStruct> self
    bitset_copy(word, &NBCS.words[i])
    return 0

cdef int refine_by_bip_degree(PartitionStack *col_ps, void *S, int *cells_to_refine_by, int ctrb_len) noexcept:
    r"""
    Refine the input partition by checking degrees of vertices to the given
    cells in the associated bipartite graph (vertices split into columns and
    words).

    INPUT:

    - ``col_ps`` -- a partition stack, whose finest partition is the partition to be refined
    - ``S`` -- a binary code struct object
    - ``cells_to_refine_by`` -- list of pointers to cells to check degrees against
      in refining the other cells (updated in place)
    - ``ctrb_len`` -- how many cells in cells_to_refine_by

    OUTPUT:

    An integer $I$ invariant under the orbits of $S_n$.  That is, if $\gamma$ is a
    permutation of the columns, then
    $$ I(G, PS, cells_to_refine_by) = I( \gamma(G), \gamma(PS), \gamma(cells_to_refine_by) ) .$$
    """
    cdef BinaryCodeStruct BCS = <BinaryCodeStruct> S
    cdef int current_cell_against = 0
    cdef int current_cell, i, r
    cdef int first_largest_subcell
    cdef int invariant = 0
    cdef PartitionStack *word_ps = BCS.word_ps
    cdef int *ctrb = BCS.alpha
    cdef bitset_s *ctrb_is_wd = BCS.alpha_is_wd

    word_ps.depth = col_ps.depth
    PS_clear(word_ps)
    bitset_zero(ctrb_is_wd)
    memcpy(ctrb, cells_to_refine_by, ctrb_len * sizeof(int))
    if BCS.first_time:
        BCS.first_time = 0
        ctrb[ctrb_len] = 0
        bitset_set(ctrb_is_wd, ctrb_len)
        ctrb_len += 1
    cdef int *col_degrees = BCS.scratch                                     # len degree
    cdef int *col_counts  = &BCS.scratch[BCS.degree]                        # len nwords+1
    cdef int *col_output  = &BCS.scratch[BCS.degree + BCS.nwords + 1]       # len degree
    cdef int *word_degrees  = &BCS.scratch[2*BCS.degree + BCS.nwords + 1]   # len nwords
    cdef int *word_counts   = &BCS.scratch[2*BCS.degree + 2*BCS.nwords + 1] # len degree+1
    cdef int *word_output   = &BCS.scratch[3*BCS.degree + 2*BCS.nwords + 2] # len nwords
    cdef bint necessary_to_split_cell
    cdef int against_index
    while not (PS_is_discrete(col_ps) and PS_is_discrete(word_ps)) and current_cell_against < ctrb_len:
        invariant += 1
        current_cell = 0
        if bitset_check(ctrb_is_wd, current_cell_against):
            while current_cell < col_ps.degree:
                invariant += 8
                i = current_cell
                necessary_to_split_cell = 0
                while True:
                    col_degrees[i-current_cell] = col_degree(col_ps, BCS, i, ctrb[current_cell_against], word_ps)
                    if col_degrees[i-current_cell] != col_degrees[0]:
                        necessary_to_split_cell = 1
                    i += 1
                    if col_ps.levels[i-1] <= col_ps.depth:
                        break
                # now, i points to the next cell (before refinement)
                if necessary_to_split_cell:
                    invariant += 8
                    first_largest_subcell = sort_by_function_codes(col_ps, current_cell, col_degrees, col_counts, col_output, BCS.nwords+1)
                    invariant += col_degree(col_ps, BCS, i-1, ctrb[current_cell_against], word_ps)
                    invariant += first_largest_subcell
                    against_index = current_cell_against
                    while against_index < ctrb_len:
                        if ctrb[against_index] == current_cell and not bitset_check(ctrb_is_wd, against_index):
                            ctrb[against_index] = first_largest_subcell
                            break
                        against_index += 1
                    r = current_cell
                    while True:
                        if r == current_cell or col_ps.levels[r-1] == col_ps.depth:
                            if r != first_largest_subcell:
                                ctrb[ctrb_len] = r
                                ctrb_len += 1
                        r += 1
                        if r >= i:
                            break
                    invariant += (i - current_cell)
                current_cell = i
        else:
            while current_cell < word_ps.degree:
                invariant += 64
                i = current_cell
                necessary_to_split_cell = 0
                while True:
                    word_degrees[i-current_cell] = word_degree(word_ps, BCS, i, ctrb[current_cell_against], col_ps)
                    if word_degrees[i-current_cell] != word_degrees[0]:
                        necessary_to_split_cell = 1
                    i += 1
                    if word_ps.levels[i-1] <= col_ps.depth:
                        break
                # now, i points to the next cell (before refinement)
                if necessary_to_split_cell:
                    invariant += 64
                    first_largest_subcell = sort_by_function_codes(word_ps, current_cell, word_degrees, word_counts, word_output, BCS.degree+1)
                    invariant += word_degree(word_ps, BCS, i-1, ctrb[current_cell_against], col_ps)
                    invariant += first_largest_subcell
                    against_index = current_cell_against
                    while against_index < ctrb_len:
                        if ctrb[against_index] == current_cell and bitset_check(ctrb_is_wd, against_index):
                            ctrb[against_index] = first_largest_subcell
                            break
                        against_index += 1
                    r = current_cell
                    while True:
                        if r == current_cell or word_ps.levels[r-1] == col_ps.depth:
                            if r != first_largest_subcell:
                                ctrb[ctrb_len] = r
                                bitset_set(ctrb_is_wd, ctrb_len)
                                ctrb_len += 1
                        r += 1
                        if r >= i:
                            break
                    invariant += (i - current_cell)
                current_cell = i
        current_cell_against += 1
    return invariant

cdef int compare_linear_codes(int *gamma_1, int *gamma_2, void *S1, void *S2, int degree) noexcept:
    r"""
    Compare gamma_1(S1) and gamma_2(S2).

    This returns:

    - -1 if gamma_1(S1) < gamma_2(S2),
    - 0 if gamma_1(S1) == gamma_2(S2),
    - 1 if gamma_1(S1) > gamma_2(S2).

    Abstractly, what this function does is relabel the basis of B by gamma_1 and
    gamma_2, run a row reduction on each, and verify that the matrices are the
    same, which holds if and only if the rowspan is the same. In practice, if
    the codes are not equal, the reductions (which are carried out in an
    interleaved way) will terminate as soon as this is discovered, and whichever
    code has a 1 in the entry in which they differ is reported as larger.

    INPUT:

    - gamma_1, gamma_2 -- list permutations (inverse)
    - S1, S2 -- binary code struct objects
    """
    cdef int i, piv_loc_1, piv_loc_2, cur_col, cur_row=0
    cdef bint is_pivot_1, is_pivot_2
    cdef LinearBinaryCodeStruct BCS1 = <LinearBinaryCodeStruct> S1
    cdef LinearBinaryCodeStruct BCS2 = <LinearBinaryCodeStruct> S2
    cdef bitset_s *basis_1 = BCS1.scratch_bitsets                   # len = dim
    cdef bitset_s *basis_2 = &BCS1.scratch_bitsets[BCS1.dimension]   # len = dim
    cdef bitset_s *pivots  = &BCS1.scratch_bitsets[2*BCS1.dimension] # len 1
    cdef bitset_s *temp  = &BCS1.scratch_bitsets[2*BCS1.dimension+1] # len 1
    for i from 0 <= i < BCS1.dimension:
        bitset_copy(&basis_1[i], &BCS1.basis[i])
        bitset_copy(&basis_2[i], &BCS2.basis[i])
    bitset_zero(pivots)
    for cur_col from 0 <= cur_col < BCS1.degree:
        is_pivot_1 = 0
        is_pivot_2 = 0
        for i from cur_row <= i < BCS1.dimension:
            if bitset_check(&basis_1[i], gamma_1[cur_col]):
                is_pivot_1 = 1
                piv_loc_1 = i
                break
        for i from cur_row <= i < BCS1.dimension:
            if bitset_check(&basis_2[i], gamma_2[cur_col]):
                is_pivot_2 = 1
                piv_loc_2 = i
                break
        if is_pivot_1 != is_pivot_2:
            return <int>is_pivot_2 - <int>is_pivot_1
        if is_pivot_1:
            bitset_set(pivots, cur_col)
            if piv_loc_1 != cur_row:
                bitset_copy(temp, &basis_1[piv_loc_1])
                bitset_copy(&basis_1[piv_loc_1], &basis_1[cur_row])
                bitset_copy(&basis_1[cur_row], temp)
            if piv_loc_2 != cur_row:
                bitset_copy(temp, &basis_2[piv_loc_2])
                bitset_copy(&basis_2[piv_loc_2], &basis_2[cur_row])
                bitset_copy(&basis_2[cur_row], temp)
            for i from 0 <= i < cur_row:
                if bitset_check(&basis_1[i], gamma_1[cur_col]):
                    bitset_xor(&basis_1[i], &basis_1[i], &basis_1[cur_row])
                if bitset_check(&basis_2[i], gamma_2[cur_col]):
                    bitset_xor(&basis_2[i], &basis_2[i], &basis_2[cur_row])
            for i from cur_row < i < BCS1.dimension:
                if bitset_check(&basis_1[i], gamma_1[cur_col]):
                    bitset_xor(&basis_1[i], &basis_1[i], &basis_1[cur_row])
                if bitset_check(&basis_2[i], gamma_2[cur_col]):
                    bitset_xor(&basis_2[i], &basis_2[i], &basis_2[cur_row])
            cur_row += 1
        else:
            for i from 0 <= i < cur_row:
                if bitset_check(&basis_1[i], gamma_1[cur_col]) != bitset_check(&basis_2[i], gamma_2[cur_col]):
                    return <int>bitset_check(&basis_2[i], gamma_2[cur_col]) - <int>bitset_check(&basis_1[i], gamma_1[cur_col])
    return 0

cdef int compare_nonlinear_codes(int *gamma_1, int *gamma_2, void *S1, void *S2, int degree) noexcept:
    r"""
    Compare gamma_1(S1) and gamma_2(S2).

    This returns:

    - -1 if gamma_1(S1) < gamma_2(S2),
    - 0 if gamma_1(S1) == gamma_2(S2),
    - 1 if gamma_1(S1) > gamma_2(S2).

    INPUT:

    - gamma_1, gamma_2 -- list permutations (inverse)
    - S1, S2 -- a binary code struct object
    """
    cdef int side=0, i, start, end, n_one_1, n_one_2, cur_col
    cdef int where_0, where_1
    cdef NonlinearBinaryCodeStruct BCS1 = <NonlinearBinaryCodeStruct> S1
    cdef NonlinearBinaryCodeStruct BCS2 = <NonlinearBinaryCodeStruct> S2
    cdef bitset_s *B_1_0 = BCS1.scratch_bitsets                   # nwords of len degree
    cdef bitset_s *B_1_1 = &BCS1.scratch_bitsets[BCS1.nwords]      # nwords of len degree
    cdef bitset_s *B_2_0 = &BCS1.scratch_bitsets[2*BCS1.nwords]    # nwords of len degree
    cdef bitset_s *B_2_1 = &BCS1.scratch_bitsets[3*BCS1.nwords]    # nwords of len degree
    cdef bitset_s *dividers = &BCS1.scratch_bitsets[4*BCS1.nwords] # 1 of len nwords
    cdef bitset_s *B_1_this
    cdef bitset_s *B_1_other
    cdef bitset_s *B_2_this
    cdef bitset_s *B_2_other
    for i from 0 <= i < BCS1.nwords:
        bitset_copy(&B_1_0[i], &BCS1.words[i])
        bitset_copy(&B_2_0[i], &BCS2.words[i])
    bitset_zero(dividers)
    bitset_set(dividers, BCS1.nwords-1)

    for cur_col from 0 <= cur_col < BCS1.degree:
        if side == 0:
            B_1_this  = B_1_0
            B_1_other = B_1_1
            B_2_this  = B_2_0
            B_2_other = B_2_1
        else:
            B_1_this  = B_1_1
            B_1_other = B_1_0
            B_2_this  = B_2_1
            B_2_other = B_2_0
        side ^= 1
        start = 0
        while start < BCS1.nwords:
            end = start
            while not bitset_check(dividers, end):
                end += 1
            end += 1
            n_one_1 = 0
            n_one_2 = 0
            for i from start <= i < end:
                n_one_1 += bitset_check(&B_1_this[i], gamma_1[cur_col])
                n_one_2 += bitset_check(&B_2_this[i], gamma_2[cur_col])
            if n_one_1 != n_one_2:
                if n_one_1 > n_one_2:
                    return 1
                else:
                    return -1
            where_0 = start
            where_1 = end - n_one_1
            if start < where_1 and where_1 < end:
                bitset_set(dividers, where_1 - 1)
            for i from start <= i < end:
                if bitset_check(&B_1_this[i], gamma_1[cur_col]):
                    bitset_copy(&B_1_other[where_1], &B_1_this[i])
                    where_1 += 1
                else:
                    bitset_copy(&B_1_other[where_0], &B_1_this[i])
                    where_0 += 1
            where_0 = start
            where_1 = end - n_one_2
            for i from start <= i < end:
                if bitset_check(&B_2_this[i], gamma_2[cur_col]):
                    bitset_copy(&B_2_other[where_1], &B_2_this[i])
                    where_1 += 1
                else:
                    bitset_copy(&B_2_other[where_0], &B_2_this[i])
                    where_0 += 1
            start = end

    return 0

cdef bint all_children_are_equivalent(PartitionStack *col_ps, void *S) noexcept:
    """
    Return ``True`` if any refinement of the current partition results in the
    same structure.

    INPUT:

    - ``col_ps`` -- the partition stack to be checked
    - ``S`` -- a binary code struct object

    .. WARNING::

        Converse does not hold in general! See Lemma 2.25 of [1] for details, noting
        that the binary code is interpreted as a bipartite graph (see module docs
        for details).
    """
    cdef BinaryCodeStruct BCS = <BinaryCodeStruct> S
    cdef PartitionStack *word_ps = BCS.word_ps
    cdef int i, n = col_ps.degree + BCS.nwords
    cdef bint in_cell = 0
    cdef int nontrivial_cells = 0
    cdef int total_cells = PS_num_cells(col_ps) + PS_num_cells(word_ps)
    if n <= total_cells + 4:
        return 1
    for i from 0 <= i < BCS.nwords:
        if word_ps.levels[i] <= col_ps.depth:
            if in_cell:
                nontrivial_cells += 1
            in_cell = 0
        else:
            in_cell = 1
    in_cell = 0
    for i from 0 <= i < BCS.degree:
        if col_ps.levels[i] <= col_ps.depth:
            if in_cell:
                nontrivial_cells += 1
            in_cell = 0
        else:
            in_cell = 1
    if n == total_cells + nontrivial_cells:
        return 1
    if n == total_cells + nontrivial_cells + 1:
        return 1
    return 0

cdef inline int word_degree(PartitionStack *word_ps, BinaryCodeStruct BCS, int entry, int cell_index, PartitionStack *col_ps) noexcept:
    """
    Return the number of edges from the vertex corresponding to entry to
    vertices in the cell corresponding to cell_index.

    INPUT:

    - ``word_ps`` -- the partition stack to be checked
    - ``col_ps`` -- corresponding partition stack on columns
    - ``BCS`` -- a binary code struct object
    - ``entry`` -- the position of the vertex in question in the entries of word_ps
    - ``cell_index`` -- the starting position of the cell in question in the entries of PS
    """
    cdef bitset_t cell, word
    cdef int h
    bitset_init(cell, BCS.degree)
    bitset_zero(cell)
    bitset_init(word, BCS.degree)
    entry = word_ps.entries[entry]
    bitset_set(cell, col_ps.entries[cell_index])
    while col_ps.levels[cell_index] > col_ps.depth:
        cell_index += 1
        bitset_set(cell, col_ps.entries[cell_index])
    BCS.ith_word(BCS, entry, word)
    bitset_and(cell, word, cell)
    h = bitset_hamming_weight(cell)
    bitset_free(cell)
    bitset_free(word)
    return h

cdef inline int col_degree(PartitionStack *col_ps, BinaryCodeStruct BCS, int entry, int cell_index, PartitionStack *word_ps) noexcept:
    """
    Return the number of edges from the vertex corresponding to entry to
    vertices in the cell corresponding to cell_index.

    INPUT:

    - ``col_ps`` -- the partition stack to be checked
    - ``word_ps`` -- corresponding partition stack on words
    - ``BCS`` -- a binary code struct object
    - ``entry`` -- the position of the vertex in question in the entries of word_ps
    - ``cell_index`` -- the starting position of the cell in question in the entries of PS
    """
    cdef bitset_t word
    bitset_init(word, BCS.degree)
    cdef int degree = 0
    entry = col_ps.entries[entry]
    while True:
        BCS.ith_word(BCS, word_ps.entries[cell_index], word)
        degree += bitset_check(word, entry)
        if not word_ps.levels[cell_index] > col_ps.depth:
            break
        cell_index += 1
    bitset_free(word)
    return degree

cdef inline int sort_by_function_codes(PartitionStack *PS, int start, int *degrees, int *counts, int *output, int count_max) noexcept:
    """
    A simple counting sort, given the degrees of vertices to a certain cell.

    INPUT:

    - ``PS`` -- the partition stack to be checked
    - ``start`` -- beginning index of the cell to be sorted
    - ``degrees`` -- the values to be sorted by
    - count, count_max, output -- scratch space
    """
    cdef int i, j, max, max_location
    for j in range(count_max):
        counts[j] = 0
    i = 0
    while PS.levels[i+start] > PS.depth:
        counts[degrees[i]] += 1
        i += 1
    counts[degrees[i]] += 1
    # i+start is the right endpoint of the cell now
    max = counts[0]
    max_location = 0
    for j from 0 < j < count_max:
        if counts[j] > max:
            max = counts[j]
            max_location = j
        counts[j] += counts[j - 1]
    for j from i >= j >= 0:
        counts[degrees[j]] -= 1
        output[counts[degrees[j]]] = PS.entries[start+j]
    max_location = counts[max_location]+start
    for j from 0 <= j <= i:
        PS.entries[start+j] = output[j]
    j = 1
    while j < count_max and counts[j] <= i:
        if counts[j] > 0:
            PS.levels[start + counts[j] - 1] = PS.depth
        PS_move_min_to_front(PS, start + counts[j-1], start + counts[j] - 1)
        j += 1
    return max_location


def random_tests(num=50, n_max=50, k_max=6, nwords_max=200, perms_per_code=10, density_range=(.1,.9)):
    """
    Test to make sure that ``C(gamma(B)) == C(B)`` for random permutations ``gamma``
    and random codes ``B``, and that :meth:`is_isomorphic` returns an isomorphism.

    INPUT:

    - ``num`` -- run tests for this many codes
    - ``n_max`` -- test codes with at most this many columns
    - ``k_max`` -- test codes with at most this for dimension
    - ``perms_per_code`` -- test each code with this many random permutations

    DISCUSSION:

    This code generates num random linear codes B on at most n_max columns with
    dimension at most k_max, and a random nonlinear code B2 on at most n_max
    columns with number of words at most nwords_max. The density of entries in
    the basis is chosen randomly between 0 and 1.

    For each code B (B2) generated, we uniformly generate perms_per_code random
    permutations and verify that the canonical labels of B and the image of B
    under the generated permutation are equal, and check that the double coset
    function returns an isomorphism.

    TESTS::

        sage: import sage.groups.perm_gps.partn_ref.refinement_binary
        sage: sage.groups.perm_gps.partn_ref.refinement_binary.random_tests()  # long time (up to 5s on sage.math, 2012)
        All passed: ... random tests on ... codes.
    """
    from sage.misc.prandom import random, randint
    from sage.combinat.permutation import Permutations
    from sage.matrix.constructor import random_matrix, matrix
    from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
    cdef int h, i, j, n, k, num_tests = 0, num_codes = 0
    cdef LinearBinaryCodeStruct B, C
    cdef NonlinearBinaryCodeStruct B_n, C_n
    for mmm in range(num):
        p = random()*(density_range[1]-density_range[0]) + density_range[0]
        n = randint(2, n_max)
        k = randint(1, min(n-1,k_max) )
        nwords = randint(1, min(n-1,nwords_max) )
        S = Permutations(n)

        M = random_matrix(GF(2), k, n, sparse=False, density=p).row_space().basis_matrix()
        M_n = random_matrix(GF(2), nwords, n, sparse=False, density=p)
        B = LinearBinaryCodeStruct( M )
        B_n = NonlinearBinaryCodeStruct( M_n )
        B.run()
        B_n.run()

        for i from 0 <= i < perms_per_code:
            perm = [a-1 for a in list(S.random_element())]
            C = LinearBinaryCodeStruct( matrix(GF(2), B.dimension, B.degree) )
            C_n = NonlinearBinaryCodeStruct( matrix(GF(2), B_n.nwords, B_n.degree) )
            for j from 0 <= j < B.dimension:
                for h from 0 <= h < B.degree:
                    bitset_set_to(&C.basis[j], <mp_bitcnt_t> perm[h], bitset_check(&B.basis[j], h))
            for j from 0 <= j < B_n.nwords:
                for h from 0 <= h < B_n.degree:
                    bitset_set_to(&C_n.words[j], <mp_bitcnt_t> perm[h], bitset_check(&B_n.words[j], h))
            # now C is a random permutation of B, and C_n of B_n
            C.run()
            C_n.run()
            B_relab = B.canonical_relabeling()
            C_relab = C.canonical_relabeling()
            B_n_relab = B_n.canonical_relabeling()
            C_n_relab = C_n.canonical_relabeling()
            B_M = matrix(GF(2), B.dimension, B.degree)
            C_M = matrix(GF(2), B.dimension, B.degree)
            B_n_M = matrix(GF(2), B_n.nwords, B_n.degree)
            C_n_M = matrix(GF(2), B_n.nwords, B_n.degree)
            for j from 0 <= j < B.dimension:
                for h from 0 <= h < B.degree:
                    B_M[j,B_relab[h]] = bitset_check(&B.basis[j], h)
                    C_M[j,C_relab[h]] = bitset_check(&C.basis[j], h)
            for j from 0 <= j < B_n.nwords:
                for h from 0 <= h < B_n.degree:
                    B_n_M[j,B_n_relab[h]] = bitset_check(&B_n.words[j], h)
                    C_n_M[j,C_n_relab[h]] = bitset_check(&C_n.words[j], h)
            if B_M.row_space() != C_M.row_space():
                print("can_lab error -- B:")
                for j from 0 <= j < B.dimension:
                    print(bitset_string(&B.basis[j]))
                print(perm)
                return
            if sorted(B_n_M.rows()) != sorted(C_n_M.rows()):
                print("can_lab error -- B_n:")
                for j from 0 <= j < B_n.nwords:
                    print(bitset_string(&B_n.words[j]))
                print(perm)
                return
            isom = B.is_isomorphic(C)
            if not isom:
                print("isom -- B:")
                for j from 0 <= j < B.dimension:
                    print(bitset_string(&B.basis[j]))
                print(perm)
                print(isom)
                return
            isom = B_n.is_isomorphic(C_n)
            if not isom:
                print("isom -- B_n:")
                for j from 0 <= j < B_n.nwords:
                    print(bitset_string(&B_n.words[j]))
                print(perm)
                print(isom)
                return

        num_tests += 4*perms_per_code
        num_codes += 2

    print("All passed: %d random tests on %d codes." % (num_tests, num_codes))
