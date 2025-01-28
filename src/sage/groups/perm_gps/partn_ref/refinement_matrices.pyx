"""
Partition backtrack functions for matrices

EXAMPLES::

    sage: import sage.groups.perm_gps.partn_ref.refinement_matrices

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

from libc.string cimport memcmp

from sage.groups.perm_gps.partn_ref.data_structures cimport *
from sage.data_structures.bitset_base cimport *
from sage.rings.integer cimport Integer
from sage.matrix.constructor import Matrix
from sage.groups.perm_gps.partn_ref.refinement_binary cimport NonlinearBinaryCodeStruct, refine_by_bip_degree
from sage.groups.perm_gps.partn_ref.double_coset cimport double_coset


cdef class MatrixStruct:

    def __cinit__(self, matrix):
        cdef int i, j
        cdef int *num_rows
        self.degree = matrix.ncols()
        self.nwords = matrix.nrows()
        cdef NonlinearBinaryCodeStruct S_temp
        self.matrix = matrix

        symbols = set(self.matrix.list())
        symbols.discard(0)
        self.symbols = sorted(symbols)
        self.nsymbols = len(symbols)

        self.symbol_structs = []
        num_rows = <int *> sig_malloc(self.nsymbols * sizeof(int))
        self.temp_col_ps = PS_new(self.degree, 1)
        if num_rows is NULL or self.temp_col_ps is NULL:
            sig_free(num_rows)
            PS_dealloc(self.temp_col_ps)
            raise MemoryError

        for i in range(self.nsymbols):
            num_rows[i] = 0
        for row in self.matrix.rows():
            row = set(row.list())
            if 0 in row:
                row.remove(0)
            for s in row:
                num_rows[self.symbols.index(s)] += 1
        for i in range(self.nsymbols):
            S_temp = NonlinearBinaryCodeStruct( (self.degree, num_rows[i]) )
            self.symbol_structs.append(S_temp)

        for i in range(self.nsymbols):
            num_rows[i] = 0
        for row in self.matrix.rows():
            row_list = row.list()
            row_list.reverse()
            for i in row.nonzero_positions():
                s = row[i]
                j = self.symbols.index(s)
                S_temp = <NonlinearBinaryCodeStruct>self.symbol_structs[j]
                bitset_set( &S_temp.words[num_rows[j]], i)
                if row_list.count(s) == 1 or row_list.index(s) == self.degree - i - 1:
                    num_rows[j] += 1
        sig_free(num_rows)
        self.output = NULL

    def __dealloc__(self):
        PS_dealloc(self.temp_col_ps)
        if self.output is not NULL:
            deallocate_agcl_output(self.output)

    def display(self):
        """
        Display the matrix, and associated data.

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref.refinement_matrices import MatrixStruct
            sage: M = MatrixStruct(Matrix(GF(5), [[0,1,1,4,4],[0,4,4,1,1]]))
            sage: M.display()
            [0 1 1 4 4]
            [0 4 4 1 1]
            <BLANKLINE>
            01100
            00011
            1
            <BLANKLINE>
            00011
            01100
            4
        """
        print(self.matrix)
        print("")
        cdef int i,j=0
        cdef NonlinearBinaryCodeStruct S_temp
        for S in self.symbol_structs:
            S_temp = <NonlinearBinaryCodeStruct>S
            for i in range(S_temp.nwords):
                print(bitset_string(&S_temp.words[i]))
            print(self.symbols[j])
            print("")
            j += 1

    def run(self, partition=None):
        """
        Perform the canonical labeling and automorphism group computation,
        storing results to ``self``.

        INPUT:

        - ``partition`` -- an optional list of lists partition of the columns;
          default is the unit partition

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref.refinement_matrices import MatrixStruct

            sage: M = MatrixStruct(matrix(GF(3),[[0,1,2],[0,2,1]]))
            sage: M.run()
            sage: M.automorphism_group()
            ([[0, 2, 1]], 2, [1])
            sage: M.canonical_relabeling()
            [0, 1, 2]

            sage: M = MatrixStruct(matrix(GF(3),[[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,0,1],[2,1,0]]))
            sage: M.automorphism_group()[1] == 6
            True

            sage: M = MatrixStruct(matrix(GF(3),[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,2]]))
            sage: M.automorphism_group()[1] == factorial(14)
            True
        """
        cdef int i, n = self.degree
        cdef PartitionStack *part
        cdef NonlinearBinaryCodeStruct S_temp
        for i in range(self.nsymbols):
            S_temp = <NonlinearBinaryCodeStruct> self.symbol_structs[i]
            S_temp.first_time = 1

        if partition is None:
            part = PS_new(n, 1)
        else:
            part = PS_from_list(partition)
        if part is NULL:
            raise MemoryError

        self.output = get_aut_gp_and_can_lab(<void *> self, part, self.degree, &all_matrix_children_are_equivalent, &refine_matrix, &compare_matrices, 1, NULL, NULL, NULL)

        PS_dealloc(part)

    def automorphism_group(self):
        """
        Return a list of generators of the automorphism group, along with its
        order and a base for which the list of generators is a strong generating
        set.

        For more examples, see self.run().

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref.refinement_matrices import MatrixStruct

            sage: M = MatrixStruct(matrix(GF(3),[[0,1,2],[0,2,1]]))
            sage: M.automorphism_group()
            ([[0, 2, 1]], 2, [1])
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

        For more examples, see self.run().

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref.refinement_matrices import MatrixStruct

            sage: M = MatrixStruct(matrix(GF(3),[[0,1,2],[0,2,1]]))
            sage: M.canonical_relabeling()
            [0, 1, 2]
        """
        cdef int i
        if self.output is NULL:
            self.run()
        return [self.output.relabeling[i] for i in range(self.degree)]

    def is_isomorphic(self, MatrixStruct other):
        """
        Calculate whether ``self`` is isomorphic to ``other``.

        EXAMPLES::

            sage: from sage.groups.perm_gps.partn_ref.refinement_matrices import MatrixStruct
            sage: M = MatrixStruct(Matrix(GF(11), [[1,2,3,0,0,0],[0,0,0,1,2,3]]))
            sage: N = MatrixStruct(Matrix(GF(11), [[0,1,0,2,0,3],[1,0,2,0,3,0]]))
            sage: M.is_isomorphic(N)
            [0, 2, 4, 1, 3, 5]
        """
        cdef int i, n = self.degree
        cdef int *output
        cdef int *ordering
        cdef PartitionStack *part
        cdef NonlinearBinaryCodeStruct S_temp
        for i in range(self.nsymbols):
            S_temp = self.symbol_structs[i]
            S_temp.first_time = 1
            S_temp = other.symbol_structs[i]
            S_temp.first_time = 1
        part = PS_new(n, 1)
        ordering = <int *> sig_malloc(self.degree * sizeof(int))
        output = <int *> sig_malloc(self.degree * sizeof(int))
        if part is NULL or ordering is NULL or output is NULL:
            PS_dealloc(part)
            sig_free(ordering)
            sig_free(output)
            raise MemoryError
        for i in range(self.degree):
            ordering[i] = i

        cdef bint isomorphic = double_coset(<void *> self, <void *> other, part, ordering, self.degree, &all_matrix_children_are_equivalent, &refine_matrix, &compare_matrices, NULL, NULL, output)

        PS_dealloc(part)
        sig_free(ordering)
        if isomorphic:
            output_py = [output[i] for i in range(self.degree)]
        else:
            output_py = False
        sig_free(output)
        return output_py

cdef int refine_matrix(PartitionStack *PS, void *S, int *cells_to_refine_by, int ctrb_len) noexcept:
    cdef MatrixStruct M = <MatrixStruct> S
    cdef int temp_inv, invariant = 1
    cdef bint changed = 1
    while changed:
        PS_copy_from_to(PS, M.temp_col_ps)
        for BCS in M.symbol_structs:
            temp_inv = refine_by_bip_degree(PS, <void *> BCS, cells_to_refine_by, ctrb_len)
            invariant *= temp_inv + 1
        if memcmp(PS.entries, M.temp_col_ps.entries, 2*M.degree * sizeof(int)) == 0:
            changed = 0
    return invariant

cdef int compare_matrices(int *gamma_1, int *gamma_2, void *S1, void *S2, int degree) noexcept:
    cdef MatrixStruct MS1 = <MatrixStruct> S1
    cdef MatrixStruct MS2 = <MatrixStruct> S2
    M1 = MS1.matrix
    M2 = MS2.matrix
    cdef int i
    MM1 = Matrix(M1.base_ring(), M1.nrows(), M1.ncols(), sparse=M1.is_sparse())
    MM2 = Matrix(M2.base_ring(), M2.nrows(), M2.ncols(), sparse=M2.is_sparse())
    for i in range(degree):
        MM1.set_column(i, M1.column(gamma_1[i]))
        MM2.set_column(i, M2.column(gamma_2[i]))
    rows1 = sorted(MM1.rows())
    rows2 = sorted(MM2.rows())
    if rows1 == rows2:
        return 0
    return -1 if rows1 < rows2 else 1

cdef bint all_matrix_children_are_equivalent(PartitionStack *PS, void *S) noexcept:
    return 0


def random_tests(n=10, nrows_max=50, ncols_max=50, nsymbols_max=10, perms_per_matrix=5, density_range=(.1,.9)):
    """
    Test to make sure that ``C(gamma(M)) == C(M)`` for random permutations ``gamma``
    and random matrices ``M``, and that ``M.is_isomorphic(gamma(M))`` returns an
    isomorphism.

    INPUT:

    - ``n`` -- run tests on this many matrices
    - ``nrows_max`` -- test matrices with at most this many rows
    - ``ncols_max`` -- test matrices with at most this many columns
    - ``perms_per_matrix`` -- test each matrix with this many random permutations
    - ``nsymbols_max`` -- maximum number of distinct symbols in the matrix

    This code generates n random matrices M on at most ncols_max columns and at
    most nrows_max rows. The density of entries in the basis is chosen randomly
    between 0 and 1.

    For each matrix M generated, we uniformly generate perms_per_matrix random
    permutations and verify that the canonical labels of M and the image of M
    under the generated permutation are equal, and that the isomorphism is
    discovered by the double coset function.

    TESTS::

        sage: import sage.groups.perm_gps.partn_ref.refinement_matrices
        sage: sage.groups.perm_gps.partn_ref.refinement_matrices.random_tests()  # long time (up to 30s on sage.math, 2011)
        All passed: ... random tests on ... matrices.
    """
    from sage.misc.prandom import random, randint
    from sage.combinat.permutation import Permutations
    from sage.matrix.constructor import random_matrix, matrix
    from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
    from sage.arith.misc import next_prime
    cdef int i, j, nrows, num_tests = 0, num_matrices = 0
    cdef MatrixStruct M, N
    for m in range(n):
        p = random()*(density_range[1]-density_range[0]) + density_range[0]
        nrows = randint(1, nrows_max)
        ncols = randint(1, ncols_max)
        nsymbols = next_prime(randint(1, nsymbols_max))
        S = Permutations(ncols)
        MM = random_matrix(GF(nsymbols), nrows, ncols, sparse=False, density=p)
        M = MatrixStruct( MM )
        M.run()

        for i in range(perms_per_matrix):
            perm = [a-1 for a in list(S.random_element())]
            NN = matrix(GF(nsymbols), nrows, ncols)
            for j from 0 <= j < ncols:
                NN.set_column(perm[j], MM.column(j))
            N = MatrixStruct(NN)
            # now N is a random permutation of M
            N.run()

            M_relab = M.canonical_relabeling()
            N_relab = N.canonical_relabeling()

            M_C = matrix(GF(nsymbols), nrows, ncols)
            N_C = matrix(GF(nsymbols), nrows, ncols)

            for j from 0 <= j < ncols:
                M_C.set_column(M_relab[j], MM.column(j))
                N_C.set_column(N_relab[j], NN.column(j))

            M_C = matrix(GF(nsymbols), sorted(M_C.rows()))
            N_C = matrix(GF(nsymbols), sorted(N_C.rows()))

            if M_C != N_C:
                print("M:")
                print(M.matrix.str())
                print("perm:")
                print(perm)
                return

            isom = M.is_isomorphic(N)
            if not isom:
                print("isom FAILURE: M:")
                print(M.matrix.str())
                print("isom FAILURE: N:")
                print(N.matrix.str())
                return

        num_tests += perms_per_matrix
        num_matrices += 2
    print("All passed: %d random tests on %d matrices." %
          (num_tests, num_matrices))
