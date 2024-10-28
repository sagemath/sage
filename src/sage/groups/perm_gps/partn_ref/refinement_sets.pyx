r"""
Partition backtrack functions for sets

EXAMPLES::

    sage: import sage.groups.perm_gps.partn_ref.refinement_sets

REFERENCE:

    [1] McKay, Brendan D. Practical Graph Isomorphism. Congressus Numerantium,
        Vol. 30 (1981), pp. 45-87.

    [2] Leon, Jeffrey. Permutation Group Algorithms Based on Partitions, I:
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
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.groups.perm_gps.partn_ref.data_structures cimport *
from sage.groups.perm_gps.partn_ref.double_coset cimport double_coset
from sage.data_structures.bitset_base cimport *


def set_stab_py(generators, sett, relab=False):
    r"""
    Compute the setwise stabilizer of a subset of [0..n-1] in a subgroup of S_n.

    EXAMPLES::

        sage: from sage.groups.perm_gps.partn_ref.refinement_sets import set_stab_py

    Degree 4 examples

    A four-cycle::

        sage: set_stab_py([[1,2,3,0]], [0])
        []
        sage: set_stab_py([[1,2,3,0]], [0,1])
        []
        sage: set_stab_py([[1,2,3,0]], [0,1,2])
        []
        sage: set_stab_py([[1,2,3,0]], [0,1,2,3])
        [[1, 2, 3, 0]]
        sage: set_stab_py([[1,2,3,0]], [0,2])
        [[2, 3, 0, 1]]

    Symmetric group::

        sage: set_stab_py([[1,0,2,3],[1,2,3,0]], [0])
        [[0, 1, 3, 2], [0, 2, 1, 3]]
        sage: set_stab_py([[1,0,2,3],[1,2,3,0]], [0,1])
        [[1, 0, 2, 3], [0, 1, 3, 2]]
        sage: set_stab_py([[1,0,2,3],[1,2,3,0]], [0,1,2,3])
        [[0, 1, 3, 2], [0, 2, 1, 3], [1, 0, 2, 3]]
        sage: set_stab_py([[1,0,2,3],[1,2,3,0]], [0,3])
        [[3, 1, 2, 0], [0, 2, 1, 3]]

    Klein 4-group::

        sage: set_stab_py([[1,0,2,3],[0,1,3,2]], [0])
        [[0, 1, 3, 2]]
        sage: set_stab_py([[1,0,2,3],[0,1,3,2]], [0,1])
        [[0, 1, 3, 2], [1, 0, 2, 3]]
        sage: set_stab_py([[1,0,2,3],[0,1,3,2]], [0,2])
        []

    Dihedral group::

        sage: set_stab_py([[1,2,3,0],[0,3,2,1]], [0])
        [[0, 3, 2, 1]]
        sage: set_stab_py([[1,2,3,0],[0,3,2,1]], [0,1])
        [[1, 0, 3, 2]]
        sage: set_stab_py([[1,2,3,0],[0,3,2,1]], [0,2])
        [[2, 1, 0, 3], [0, 3, 2, 1]]
        sage: set_stab_py([[1,2,3,0],[0,3,2,1]], [1])
        [[2, 1, 0, 3]]
        sage: set_stab_py([[1,2,3,0],[0,3,2,1]], [1,2,3])
        [[0, 3, 2, 1]]

    Alternating group::

        sage: set_stab_py([[1,2,0,3],[0,2,3,1]], [0])
        [[0, 2, 3, 1]]
        sage: set_stab_py([[1,2,0,3],[0,2,3,1]], [1])
        [[2, 1, 3, 0]]
        sage: set_stab_py([[1,2,0,3],[0,2,3,1]], [2])
        [[1, 3, 2, 0]]
        sage: set_stab_py([[1,2,0,3],[0,2,3,1]], [3])
        [[1, 2, 0, 3]]
        sage: set_stab_py([[1,2,0,3],[0,2,3,1]], [0,1])
        [[1, 0, 3, 2]]
        sage: set_stab_py([[1,2,0,3],[0,2,3,1]], [0,2])
        [[2, 3, 0, 1]]
        sage: set_stab_py([[1,2,0,3],[0,2,3,1]], [0,3])
        [[3, 2, 1, 0]]

    Larger degree examples

    Dihedral group of degree 5::

        sage: set_stab_py([[1,2,3,4,0],[0,4,3,2,1]], [])
        [[0, 4, 3, 2, 1], [1, 0, 4, 3, 2]]
        sage: set_stab_py([[1,2,3,4,0],[0,4,3,2,1]], [0])
        [[0, 4, 3, 2, 1]]
        sage: set_stab_py([[1,2,3,4,0],[0,4,3,2,1]], [0,2])
        [[2, 1, 0, 4, 3]]

    Dihedral group of degree 6::

        sage: set_stab_py([[1,2,3,4,5,0],[0,5,4,3,2,1]], [])
        [[0, 5, 4, 3, 2, 1], [1, 0, 5, 4, 3, 2]]
        sage: set_stab_py([[1,2,3,4,5,0],[0,5,4,3,2,1]], [0])
        [[0, 5, 4, 3, 2, 1]]
        sage: set_stab_py([[1,2,3,4,5,0],[0,5,4,3,2,1]], [0,1])
        [[1, 0, 5, 4, 3, 2]]
        sage: set_stab_py([[1,2,3,4,5,0],[0,5,4,3,2,1]], [0,2])
        [[2, 1, 0, 5, 4, 3]]
        sage: set_stab_py([[1,2,3,4,5,0],[0,5,4,3,2,1]], [0,3])
        [[0, 5, 4, 3, 2, 1], [3, 2, 1, 0, 5, 4]]
        sage: set_stab_py([[1,2,3,4,5,0],[0,5,4,3,2,1]], [0,2,4])
        [[2, 1, 0, 5, 4, 3], [4, 3, 2, 1, 0, 5]]

    Canonical labels::

        sage: set_stab_py([[0,2,1,4,3,5,8,7,6],[8,7,6,3,5,4,2,1,0]], [0,1,3,5,6], True)
        ([], [7, 8, 6, 3, 4, 5, 2, 0, 1])
        sage: set_stab_py([[0,2,1,4,3,5,8,7,6],[8,7,6,3,5,4,2,1,0]], [0,3,5,6,8], True)
        ([], [2, 1, 0, 5, 4, 3, 7, 6, 8])
    """
    if len(generators) == 0:
        return []
    cdef int i, j, n = len(generators[0]), n_gens = len(generators)
    cdef StabilizerChain *supergroup = SC_new(n)
    cdef aut_gp_and_can_lab *stabilizer
    cdef int *gens = <int *> sig_malloc(n*n_gens * sizeof(int))
    cdef subset *subset_sett = <subset *> sig_malloc(sizeof(subset))
    if gens is NULL or supergroup is NULL or subset_sett is NULL:
        SC_dealloc(supergroup)
        sig_free(gens)
        sig_free(subset_sett)
        raise MemoryError
    bitset_init(&subset_sett.bits, n)
    subset_sett.scratch = <int *> sig_malloc((3*n+1) * sizeof(int))
    for i from 0 <= i < len(generators):
        for j from 0 <= j < n:
            gens[n*i + j] = generators[i][j]
    if SC_insert(supergroup, 0, gens, n_gens):
        SC_dealloc(supergroup)
        sig_free(gens)
        sig_free(subset_sett)
        raise MemoryError
    sig_free(gens)
    bitset_clear(&subset_sett.bits)
    for i in sett:
        bitset_add(&subset_sett.bits, i)
    stabilizer = set_stab(supergroup, subset_sett, relab)
    SC_dealloc(supergroup)
    bitset_free(&subset_sett.bits)
    sig_free(subset_sett.scratch)
    sig_free(subset_sett)
    if stabilizer is NULL:
        raise MemoryError
    stab_gens = []
    for i from 0 <= i < stabilizer.num_gens:
        stab_gens.append([stabilizer.generators[i*n+j] for j from 0 <= j < n])
    if relab:
        relabeling = [stabilizer.relabeling[j] for j from 0 <= j < n]
    deallocate_agcl_output(stabilizer)
    if relab:
        return stab_gens, relabeling
    return stab_gens


cdef aut_gp_and_can_lab *set_stab(StabilizerChain *supergroup, subset *sett, bint relab) noexcept:
    r"""
    Compute the set stabilizer of ``sett`` within ``supergroup``. (Note that
    ``set`` is a reserved Python keyword.) If ``relab`` is specified then
    compute the canonical label of the set under the action of the group.
    """
    cdef aut_gp_and_can_lab *output
    cdef int n = supergroup.degree
    cdef PartitionStack *part = PS_new(n, 1)
    if part is NULL:
        return NULL
    output = get_aut_gp_and_can_lab(<void *> sett, part, n,
        &all_set_children_are_equivalent, &refine_set, &compare_sets, relab,
        supergroup, NULL, NULL)
    PS_dealloc(part)
    if output is NULL:
        return NULL
    return output


def sets_isom_py(generators, set1, set2):
    r"""
    Compute whether ``set1`` and ``set2`` are isomorphic under the action of
    the group generated by the generators given in list form.

    EXAMPLES::

        sage: from sage.groups.perm_gps.partn_ref.refinement_sets import sets_isom_py

    Degree 3 examples

    Trivial group::

        sage: sets_isom_py([], [0,1,2], [0,1])
        False
        sage: sets_isom_py([], [0,1,2], [0,2,1])
        [0, 1, 2]
        sage: sets_isom_py([[0,1,2]], [0,1,2], [0,2,1])
        [0, 1, 2]
        sage: sets_isom_py([[0,1,2]], [0], [0])
        [0, 1, 2]
        sage: sets_isom_py([[0,1,2]], [0], [1])
        False
        sage: sets_isom_py([[0,1,2]], [0], [2])
        False
        sage: sets_isom_py([[0,1,2]], [0,1], [1,0])
        [0, 1, 2]

    Three-cycle::

        sage: sets_isom_py([[1,2,0]], [0], [1])
        [1, 2, 0]
        sage: sets_isom_py([[1,2,0]], [0], [2])
        [2, 0, 1]
        sage: sets_isom_py([[1,2,0]], [0], [0])
        [0, 1, 2]
        sage: sets_isom_py([[1,2,0]], [0,1], [0])
        False
        sage: sets_isom_py([[1,2,0]], [0,1], [1])
        False
        sage: sets_isom_py([[1,2,0]], [0,1], [2])
        False
        sage: sets_isom_py([[1,2,0]], [0,1], [0,2])
        [2, 0, 1]
        sage: sets_isom_py([[1,2,0]], [0,1], [1,2])
        [1, 2, 0]
        sage: sets_isom_py([[1,2,0]], [0,1], [1,0])
        [0, 1, 2]
        sage: sets_isom_py([[1,2,0]], [0,2,1], [2,1,0])
        [0, 1, 2]

    Transposition::

        sage: sets_isom_py([[1,0,2]], [0], [])
        False
        sage: sets_isom_py([[1,0,2]], [0], [1,2])
        False
        sage: sets_isom_py([[1,0,2]], [0], [1])
        [1, 0, 2]
        sage: sets_isom_py([[1,0,2]], [0], [0])
        [0, 1, 2]
        sage: sets_isom_py([[1,0,2]], [0,1], [2])
        False
        sage: sets_isom_py([[1,0,2]], [0,2], [1,2])
        [1, 0, 2]
        sage: sets_isom_py([[1,0,2]], [0], [2])
        False
        sage: sets_isom_py([[1,0,2]], [0,1], [1,2])
        False

    Symmetric group S_3::

        sage: sets_isom_py([[1,0,2],[1,2,0]], [], [])
        [0, 1, 2]
        sage: sets_isom_py([[1,0,2],[1,2,0]], [0], [])
        False
        sage: sets_isom_py([[1,0,2],[1,2,0]], [0], [0])
        [0, 1, 2]
        sage: sets_isom_py([[1,0,2],[1,2,0]], [0], [1])
        [1, 0, 2]
        sage: sets_isom_py([[1,0,2],[1,2,0]], [0], [2])
        [2, 0, 1]
        sage: sets_isom_py([[1,0,2],[1,2,0]], [0,2], [2])
        False
        sage: sets_isom_py([[1,0,2],[1,2,0]], [0,2], [1,2])
        [1, 0, 2]
        sage: sets_isom_py([[1,0,2],[1,2,0]], [0,2], [0,1])
        [0, 2, 1]
        sage: sets_isom_py([[1,0,2],[1,2,0]], [0,2], [0,2])
        [0, 1, 2]
        sage: sets_isom_py([[1,0,2],[1,2,0]], [0,2], [0,1,2])
        False
        sage: sets_isom_py([[1,0,2],[1,2,0]], [0,2,1], [0,1,2])
        [0, 1, 2]

    Degree 4 examples

    Trivial group::

        sage: sets_isom_py([[0,1,2,3]], [], [])
        [0, 1, 2, 3]
        sage: sets_isom_py([[0,1,2,3]], [0], [])
        False
        sage: sets_isom_py([[0,1,2,3]], [0], [1])
        False
        sage: sets_isom_py([[0,1,2,3]], [0], [2])
        False
        sage: sets_isom_py([[0,1,2,3]], [0], [3])
        False
        sage: sets_isom_py([[0,1,2,3]], [0], [0])
        [0, 1, 2, 3]
        sage: sets_isom_py([[0,1,2,3]], [0,1], [1,2])
        False
        sage: sets_isom_py([[0,1,2,3]], [0,1], [0,1])
        [0, 1, 2, 3]
        sage: sets_isom_py([[0,1,2,3]], [0,1,2,3], [0,1])
        False
        sage: sets_isom_py([[0,1,2,3]], [0,1,2,3], [0,1,2,3])
        [0, 1, 2, 3]

    Four-cycle::

        sage: sets_isom_py([[1,2,3,0]], [0,1,2,3], [0,1,2,3])
        [0, 1, 2, 3]
        sage: sets_isom_py([[1,2,3,0]], [], [])
        [0, 1, 2, 3]
        sage: sets_isom_py([[1,2,3,0]], [0], [0])
        [0, 1, 2, 3]
        sage: sets_isom_py([[1,2,3,0]], [0], [1])
        [1, 2, 3, 0]
        sage: sets_isom_py([[1,2,3,0]], [0], [2])
        [2, 3, 0, 1]
        sage: sets_isom_py([[1,2,3,0]], [0], [3])
        [3, 0, 1, 2]
        sage: sets_isom_py([[1,2,3,0]], [0,1], [3])
        False
        sage: sets_isom_py([[1,2,3,0]], [0,1], [])
        False
        sage: sets_isom_py([[1,2,3,0]], [0,1], [1,2,3])
        False
        sage: sets_isom_py([[1,2,3,0]], [0,1], [1,2])
        [1, 2, 3, 0]
        sage: sets_isom_py([[1,2,3,0]], [0,1], [2,3])
        [2, 3, 0, 1]
        sage: sets_isom_py([[1,2,3,0]], [0,1], [0,3])
        [3, 0, 1, 2]
        sage: sets_isom_py([[1,2,3,0]], [0,2], [0,2])
        [0, 1, 2, 3]
        sage: sets_isom_py([[1,2,3,0]], [0,2], [1,3])
        [3, 0, 1, 2]
        sage: sets_isom_py([[1,2,3,0]], [0,1,2], [1,2,3])
        [1, 2, 3, 0]
        sage: sets_isom_py([[1,2,3,0]], [0,1,2], [0,2,3])
        [2, 3, 0, 1]
        sage: sets_isom_py([[1,2,3,0]], [0,1,2], [0,1,3])
        [3, 0, 1, 2]
        sage: sets_isom_py([[1,2,3,0]], [0,1,2], [0,1,2])
        [0, 1, 2, 3]

    Two transpositions::

        sage: from sage.groups.perm_gps.partn_ref.refinement_sets import sets_isom_py
        sage: sets_isom_py([[2,3,0,1]], [0], [2])
        [2, 3, 0, 1]
        sage: sets_isom_py([[2,3,0,1]], [1], [3])
        [2, 3, 0, 1]
        sage: sets_isom_py([[2,3,0,1]], [1], [1])
        [0, 1, 2, 3]
        sage: sets_isom_py([[2,3,0,1]], [1], [2])
        False
        sage: sets_isom_py([[2,3,0,1]], [0,3], [1,2])
        [2, 3, 0, 1]
        sage: sets_isom_py([[2,3,0,1]], [0,3], [3,0])
        [0, 1, 2, 3]
        sage: sets_isom_py([[2,3,0,1]], [0,1,3], [0,2,3])
        False
        sage: sets_isom_py([[2,3,0,1]], [0,1,3], [1,2,3])
        [2, 3, 0, 1]
    """
    set1 = set(set1)
    set2 = set(set2)
    if not generators:
        if set1 == set2:
            return list(range(max(set1) + 1))
        else:
            return False

    cdef int i, j, n = len(generators[0]), n_gens = len(generators)
    cdef StabilizerChain *supergroup = SC_new(n)
    cdef int *gens = <int *> sig_malloc(n*n_gens * sizeof(int))
    cdef int *isom = <int *> sig_malloc(n * sizeof(int))
    cdef subset *subset_sett1 = <subset *> sig_malloc(sizeof(subset))
    cdef subset *subset_sett2 = <subset *> sig_malloc(sizeof(subset))
    bitset_init(&subset_sett1.bits, n)
    bitset_init(&subset_sett2.bits, n)
    subset_sett1.scratch = <int *> sig_malloc((3*n+1) * sizeof(int))
    subset_sett2.scratch = <int *> sig_malloc((3*n+1) * sizeof(int))
    for i from 0 <= i < len(generators):
        for j from 0 <= j < n:
            gens[n*i + j] = generators[i][j]
    if SC_insert(supergroup, 0, gens, n_gens):
        raise MemoryError
    sig_free(gens)
    bitset_clear(&subset_sett1.bits)
    bitset_clear(&subset_sett2.bits)
    for i in set1:
        bitset_add(&subset_sett1.bits, i)
    for i in set2:
        bitset_add(&subset_sett2.bits, i)
    cdef bint isomorphic = sets_isom(supergroup, subset_sett1, subset_sett2, isom)
    SC_dealloc(supergroup)
    bitset_free(&subset_sett1.bits)
    bitset_free(&subset_sett2.bits)
    sig_free(subset_sett1.scratch)
    sig_free(subset_sett2.scratch)
    sig_free(subset_sett1)
    sig_free(subset_sett2)
    if isomorphic:
        output_py = [isom[i] for i from 0 <= i < n]
    else:
        output_py = False
    sig_free(isom)
    return output_py


cdef int sets_isom(StabilizerChain *supergroup, subset *set1, subset *set2, int *isom) except -1:
    r"""
    Underlying C function for testing two sets for isomorphism.
    """
    cdef int n = supergroup.degree
    cdef bint x
    cdef PartitionStack *part = PS_new(n, 1)
    if part is NULL:
        raise MemoryError
    x = double_coset(set1, set2, part, NULL, n,
        &all_set_children_are_equivalent, &refine_set, &compare_sets,
        supergroup, NULL, isom)
    PS_dealloc(part)
    return x

cdef bint all_set_children_are_equivalent(PartitionStack *PS, void *S) noexcept:
    return 0

cdef int refine_set(PartitionStack *PS, void *S, int *cells_to_refine_by, int ctrb_len) noexcept:
    """
    Given a set S, refine the partition stack PS so that each cell contains
    elements which are all either in the set or not in the set. If the depth is
    positive we do nothing since the cells will have already been split at an
    earlier level.
    """
    if PS.depth > 0:
        return 0
    cdef subset *subset1 = <subset *> S
    cdef int *scratch = subset1.scratch
    cdef int start, i, n = PS.degree
    start = 0
    while start < n:
        i = 0
        while True:
            scratch[i] = bitset_in(&subset1.bits, PS.entries[start+i])
            if PS.levels[start+i] <= PS.depth:
                break
            i += 1
        sort_by_function(PS, start, scratch)
        start += i + 1
    return 0

cdef inline int _bint_cmp(bint a, bint b) noexcept:
    return (<int> b) - (<int> a)

cdef int compare_sets(int *gamma_1, int *gamma_2, void *S1, void *S2, int degree) noexcept:
    r"""
    Compare two sets according to the lexicographic order.
    """
    cdef subset *subset1 = <subset *> S1
    cdef subset *subset2 = <subset *> S2
    cdef bitset_s set1 = subset1.bits
    cdef bitset_s set2 = subset2.bits
    cdef int i, j
    for i in range(degree):
        j = _bint_cmp(bitset_in(&set1, gamma_1[i]), bitset_in(&set2, gamma_2[i]))
        if j != 0:
            return j
    return 0

cdef void *allocate_subset(int n) noexcept:
    r"""
    Allocate a subset struct of degree n.
    """
    cdef subset *set1 = <subset *> sig_malloc(sizeof(subset))
    cdef int *scratch = <int *> sig_malloc((3*n+1) * sizeof(int))
    if set1 is NULL or scratch is NULL:
        sig_free(set1)
        sig_free(scratch)
        return NULL
    try:
        bitset_init(&set1.bits, n)
    except MemoryError:
        sig_free(set1)
        sig_free(scratch)
        return NULL
    set1.scratch = scratch
    return <void *> set1

cdef void free_subset(void *child) noexcept:
    r"""
    Deallocate a subset struct.
    """
    cdef subset *set1 = <subset *> child
    if set1 is not NULL:
        sig_free(set1.scratch)
        bitset_free(&set1.bits)
    sig_free(set1)

cdef void *allocate_sgd(int degree) noexcept:
    r"""
    Allocate the data part of an iterator which generates augmentations, i.e.,
    elements to add to the set.
    """
    cdef subset_generator_data *sgd = <subset_generator_data *> sig_malloc(sizeof(subset_generator_data))
    sgd.orbits = OP_new(degree)
    if sgd is NULL or sgd.orbits is NULL:
        deallocate_sgd(sgd)
        return NULL
    return <void *> sgd

cdef void deallocate_sgd(void *data) noexcept:
    r"""
    Deallocate the data part of the augmentation iterator.
    """
    cdef subset_generator_data *sgd = <subset_generator_data *> data
    if sgd is not NULL:
        OP_dealloc(sgd.orbits)
    sig_free(sgd)

cdef void *subset_generator_next(void *data, int *degree, bint *mem_err) noexcept:
    r"""
    Return the next element to consider adding to the set.
    """
    cdef subset_generator_data *sgd = <subset_generator_data *> data
    while True:
        sgd.cur_point += 1
        if sgd.cur_point == sgd.orbits.degree:
            break
        if OP_find(sgd.orbits, sgd.cur_point) == sgd.cur_point and \
          not bitset_in(&sgd.bits, sgd.cur_point):
            break
    if sgd.cur_point == sgd.orbits.degree or mem_err[0]:
        return NULL
    return <void *> &sgd.cur_point

cdef int generate_child_subsets(void *S, aut_gp_and_can_lab *group, iterator *child_iterator) noexcept:
    r"""
    Set up an iterator of augmentations, i.e., elements to add to the given set.
    """
    cdef subset *subset1 = <subset *> S
    cdef int i, j, n = group.group.degree
    cdef subset_generator_data *sgd = <subset_generator_data *> child_iterator.data
    OP_clear(sgd.orbits)
    for i in range(group.num_gens):
        for j in range(n):
            OP_join(sgd.orbits, j, group.generators[n*i + j])
    i = bitset_first(&subset1.bits)
    j = bitset_next(&subset1.bits, i+1)
    while j != -1:
        OP_join(sgd.orbits, i, j)
        j = bitset_next(&subset1.bits, j+1)
    sgd.cur_point = -1
    sgd.bits = subset1.bits
    return 0

cdef void *apply_subset_aug(void *parent, void *aug, void *child, int *degree, bint *mem_err) noexcept:
    r"""
    Add the element represented by ``aug`` to ``parent``, storing the result to
    ``child``.
    """
    cdef subset *set1 = <subset *> child
    cdef subset *par_set = <subset *> parent
    cdef bitset_s parbits = par_set.bits
    cdef int add_pt = (<int *> aug)[0], n = parbits.size
    bitset_copy(&set1.bits, &parbits)
    bitset_add(&set1.bits, add_pt)
    degree[0] = n
    return <void *> set1

cdef void free_subset_aug(void *aug) noexcept:
    return

cdef void *canonical_set_parent(void *child, void *parent, int *permutation, int *degree, bint *mem_err) noexcept:
    r"""
    Determine the canonical parent of the set ``child`` by applying
    ``permutation``, deleting the largest element in lexicographic order, and
    storing the result to ``parent``.
    """
    cdef subset *set1 = <subset *> child
    cdef int i, max_in_can_lab, max_loc, n = set1.bits.size
    cdef subset *par
    if parent is NULL:
        par = <subset *> allocate_subset(n)
        if par is NULL:
            mem_err[0] = 1
    else:
        par = <subset *> parent
    if mem_err[0]:
        return NULL
    i = bitset_first(&set1.bits)
    max_in_can_lab = permutation[i]
    max_loc = i
    while i != -1:
        if max_in_can_lab < permutation[i]:
            max_in_can_lab = permutation[i]
            max_loc = i
        i = bitset_next(&set1.bits, i+1)
    bitset_copy(&par.bits, &set1.bits)
    bitset_discard(&par.bits, max_loc)
    degree[0] = n
    return <void *> par

cdef iterator *allocate_subset_gen(int degree, int max_size) noexcept:
    r"""
    Allocate the generator of subsets.
    """
    cdef iterator *subset_gen = <iterator *> sig_malloc(sizeof(iterator))
    if subset_gen is not NULL:
        if allocate_subset_gen_2(degree, max_size, subset_gen):
            sig_free(subset_gen)
            subset_gen = NULL
    return subset_gen

cdef int allocate_subset_gen_2(int degree, int max_size, iterator *it) noexcept:
    r"""
    Given an already allocated iterator, allocate the generator of subsets.
    """
    cdef canonical_generator_data *cgd = allocate_cgd(max_size + 1, degree)
    if cgd is NULL:
        return 1
    cdef int i, j
    for i from 0 <= i < max_size + 1:
        cgd.object_stack[i] = allocate_subset(degree)
        cgd.parent_stack[i] = allocate_subset(degree)
        cgd.iterator_stack[i].data = allocate_sgd(degree)
        cgd.iterator_stack[i].next = &subset_generator_next
        if cgd.iterator_stack[i].data is NULL or \
           cgd.object_stack[i]        is NULL or \
           cgd.parent_stack[i]        is NULL:
            for j from 0 <= j <= i:
                deallocate_sgd(cgd.iterator_stack[i].data)
                free_subset(cgd.object_stack[i])
                free_subset(cgd.parent_stack[i])
            deallocate_cgd(cgd)
            return 1
    it.data = <void *> cgd
    it.next = canonical_generator_next
    return 0

cdef void free_subset_gen(iterator *subset_gen) noexcept:
    r"""
    Free the iterator of subsets.
    """
    if subset_gen is NULL:
        return
    cdef canonical_generator_data *cgd = <canonical_generator_data *> subset_gen.data
    deallocate_cgd(cgd)
    sig_free(subset_gen)

cdef iterator *setup_set_gen(iterator *subset_gen, int degree, int max_size) noexcept:
    r"""
    Initiate the iterator of subsets.
    """
    cdef subset *empty_set
    cdef iterator *subset_iterator = setup_canonical_generator(degree,
        &all_set_children_are_equivalent,
        &refine_set,
        &compare_sets,
        &generate_child_subsets,
        &apply_subset_aug,
        &free_subset,
        &deallocate_sgd,
        &free_subset_aug,
        &canonical_set_parent,
        max_size+1, 0, subset_gen)
    if subset_iterator is not NULL:
        empty_set = <subset *> (<canonical_generator_data *> subset_gen.data).object_stack[0]
        bitset_clear(&empty_set.bits)
    return subset_iterator


def sets_modulo_perm_group(list generators, int max_size,
                           bint indicate_mem_err=1):
    r"""
    Given generators of a permutation group, list subsets up to permutations in
    the group.

    INPUT:

    - ``generators`` -- list of generators in list form
    - ``max_size`` -- integer; maximum size of subsets to be generated
    - ``indicate_mem_err`` -- boolean; whether to raise an error.
      If we run out of memory, or simply append a :exc:`MemoryError`
      instance to the end of the output.

    EXAMPLES::

        sage: from sage.groups.perm_gps.partn_ref.refinement_sets import sets_modulo_perm_group
        sage: sets_modulo_perm_group([], 0)
        [[]]
        sage: sets_modulo_perm_group([], 1)
        [[0], []]
        sage: sets_modulo_perm_group([], 2)
        [[0, 1], [0], []]
        sage: sets_modulo_perm_group([], 3)
        [[0, 1, 2], [0, 1], [0], []]
        sage: sets_modulo_perm_group([], 4)
        [[0, 1, 2, 3], [0, 1, 2], [0, 1], [0], []]
        sage: len(sets_modulo_perm_group([], 99))
        100

    ::

        sage: sets_modulo_perm_group([[1,2,0]], 4)
        [[0, 1, 2], [0, 1], [0], []]
        sage: sets_modulo_perm_group([[1,2,0]], 3)
        [[0, 1, 2], [0, 1], [0], []]
        sage: sets_modulo_perm_group([[1,2,0]], 2)
        [[0, 1], [0], []]
        sage: sets_modulo_perm_group([[1,2,0]], 1)
        [[0], []]
        sage: sets_modulo_perm_group([[1,2,0]], 0)
        [[]]
        sage: sets_modulo_perm_group([[0,1,2]], 3)
        [[0, 1, 2], [0, 1], [0, 2], [0], [1, 2], [1], [2], []]
        sage: sets_modulo_perm_group([[1,0,2]], 3)
        [[0, 1, 2], [0, 1], [0, 2], [0], [2], []]
        sage: sets_modulo_perm_group([[1,0,2],[1,2,0]], 3)
        [[0, 1, 2], [0, 1], [1], []]

    ::

        sage: sets_modulo_perm_group([[1,2,3,0]], 4)
        [[0, 1, 2, 3], [0, 1, 2], [0, 1], [0, 2], [0], []]
        sage: sets_modulo_perm_group([[1,2,3,0]], 5)
        [[0, 1, 2, 3], [0, 1, 2], [0, 1], [0, 2], [0], []]
        sage: sets_modulo_perm_group([[1,2,3,0]], 3)
        [[0, 1, 2], [0, 1], [0, 2], [0], []]
        sage: sets_modulo_perm_group([[1,2,3,0]], 2)
        [[0, 1], [0, 2], [0], []]
        sage: sets_modulo_perm_group([[1,2,3,0]], 1)
        [[0], []]
        sage: sets_modulo_perm_group([[1,2,3,0]], 0)
        [[]]
        sage: sets_modulo_perm_group([[0,1,3,2],[1,0,2,3]], 4)
        [[0, 1, 2, 3], [0, 1, 2], [0, 1], [0, 2, 3], [0, 2], [0], [2, 3], [2], []]
        sage: sets_modulo_perm_group([[1,0,2,3],[1,2,0,3]], 4)
        [[0, 1, 2, 3], [0, 1, 2], [0, 1, 3], [0, 1], [1, 3], [1], [3], []]
        sage: sets_modulo_perm_group([[1,2,0,3],[0,2,3,1]], 4)
        [[0, 1, 2, 3], [0, 1, 2], [0, 1], [1], []]
        sage: sets_modulo_perm_group([[1,0,2,3],[1,2,3,0]], 4)
        [[0, 1, 2, 3], [0, 1, 2], [1, 2], [2], []]
        sage: L = list(powerset(range(4)))
        sage: L.sort()
        sage: L == sorted(sets_modulo_perm_group([[0,1,2,3]], 4))
        True

    ::

        sage: sets_modulo_perm_group([[1,2,3,4,0]], 5)
        [[0, 1, 2, 3, 4], [0, 1, 2, 3], [0, 1, 2], [0, 1], [0, 2, 3], [0, 2], [0], []]
        sage: sets_modulo_perm_group([[1,0,2,3,4],[0,1,3,4,2]], 5)
        [[0, 1, 2, 3, 4], [0, 1, 2, 3], [0, 1, 2], [0, 1], [0, 2, 3, 4], [0, 2, 3], [0, 2], [0], [2, 3, 4], [2, 3], [2], []]
        sage: L = list(powerset(range(5)))
        sage: L.sort()
        sage: L == sorted(sets_modulo_perm_group([[0,1,2,3,4]], 5))
        True
        sage: sets_modulo_perm_group([[1,0,2,3,4],[1,2,3,4,0]], 5)
        [[0, 1, 2, 3, 4], [0, 1, 2, 3], [1, 2, 3], [2, 3], [3], []]
        sage: sets_modulo_perm_group([[1,2,0,3,4],[0,2,3,1,4],[0,1,3,4,2]], 5)
        [[0, 1, 2, 3, 4], [0, 1, 2, 3], [1, 2, 3], [1, 2], [2], []]

    ::

        sage: X = sets_modulo_perm_group([[1,2,3,4,5,0]], 6)
        sage: [a for a in X if len(a) == 0]
        [[]]
        sage: [a for a in X if len(a) == 1]
        [[0]]
        sage: [a for a in X if len(a) == 2]
        [[0, 1], [0, 2], [0, 3]]
        sage: [a for a in X if len(a) == 3]
        [[0, 1, 2], [0, 1, 3], [0, 2, 3], [0, 2, 4]]
        sage: [a for a in X if len(a) == 4]
        [[0, 1, 2, 3], [0, 1, 3, 4], [0, 2, 3, 4]]
        sage: [a for a in X if len(a) == 5]
        [[0, 1, 2, 3, 4]]
        sage: [a for a in X if len(a) == 6]
        [[0, 1, 2, 3, 4, 5]]

    ::

        sage: X = sets_modulo_perm_group([[0,2,1,4,3,5,8,7,6],[8,7,6,3,5,4,2,1,0]], 9)
        sage: len(X)
        74
    """
    cdef list out_list = []
    cdef int i
    if max_size == 0:
        return [[]]
    if len(generators) == 0:
        ll = []
        for i in range(max_size,-1,-1):
            ll.append(list(range(i)))
        return ll
    cdef int n = len(generators[0]), n_gens = len(generators)
    cdef iterator *subset_iterator
    cdef subset *thing

    cdef StabilizerChain *group = SC_new(n)
    cdef int *gens = <int *> sig_malloc(n*n_gens * sizeof(int))
    if group is NULL or gens is NULL:
        SC_dealloc(group)
        sig_free(gens)
        raise MemoryError
    for i from 0 <= i < len(generators):
        for j from 0 <= j < n:
            gens[n*i + j] = generators[i][j]
    if SC_insert(group, 0, gens, n_gens):
        SC_dealloc(group)
        sig_free(gens)
        raise MemoryError
    sig_free(gens)

    cdef iterator *subset_gen = allocate_subset_gen(n, max_size)
    if subset_gen is NULL:
        SC_dealloc(group)
        raise MemoryError
    subset_iterator = setup_set_gen(subset_gen, n, max_size)
    cdef bint mem_err = 0
    if subset_iterator is NULL:
        SC_dealloc(group)
        free_subset_gen(subset_gen)
        mem_err = 1
    else:
        start_canonical_generator(group, NULL, n, subset_gen)
    while not mem_err:
        thing = <subset *> subset_iterator.next(subset_iterator.data, NULL, &mem_err)
        if thing is NULL:
            break
        out_list.append( bitset_list(&thing.bits) )
    free_subset_gen(subset_gen)
    SC_dealloc(group)
    if mem_err:
        if indicate_mem_err:
            raise MemoryError
        else:
            out_list.append(MemoryError())
    return out_list
