#*****************************************************************************
#       Copyright (C) 2006 - 2011 Robert L. Miller <rlmillster@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.data_structures.bitset_base cimport *

from libc.string cimport memcpy
from libc.stdlib cimport rand
from sage.libs.gmp.mpz cimport *


cdef enum:
    # The following is for the automorphism group computation, says what the
    # length of fixed point and minimal cell representative arrays should be,
    # which are used to prune the search tree under known symmetries.
    # Increasing this may make certain automorphism_group_and_canonical_label
    # computations faster for incredibly large groups.
    len_of_fp_and_mcr = 100

cdef struct OrbitPartition:
    # Disjoint set data structure for representing the orbits of the generators
    # so far found. Also keeps track of the minimum elements of the cells and
    # their sizes.
    int degree
    int num_cells
    int *parent
    int *rank
    int *mcr # minimum cell representatives - only valid at the root of a cell
    int *size # also only valid at the root of a cell

cdef struct PartitionStack:
    # Representation of a node of the search tree. A sequence of partitions of
    # length depth + 1, each of which is finer than the last. Partition k is
    # represented as PS.entries in order, broken up immediately after each
    # entry of levels which is at most k.
    int *entries
    int *levels
    int depth
    int degree

cdef struct StabilizerChain:
    # A representation of a permutation group acting on 0, 1, ..., degree-1.
    int degree
    int base_size

    int *orbit_sizes
    int *num_gens      # dimension of generator cube on each level
    int *array_size    # size of space to hold generators on each level (number of permutations)

    int **base_orbits #
    int **parents     # three n*n squares, orbits and tree structures
    int **labels      #

    int **generators   # generators for each level,
    int **gen_inverses # and their inverses

    bitset_s gen_used
    bitset_s gen_is_id
    int *perm_scratch
    OrbitPartition *OP_scratch


# OrbitPartition (OP)

cdef OrbitPartition *OP_new(int n) noexcept

cdef void OP_dealloc(OrbitPartition *OP) noexcept

cdef inline int OP_copy_from_to(OrbitPartition *OP, OrbitPartition *OP2) noexcept:
    """
    Copy all data from OP to OP2, we suppose that

    -   OP2.degree == OP.degree
    -   OP2.num_cells == OP.num_cells
    """
    memcpy(OP2.parent, OP.parent, 4*OP.degree * sizeof(int) )

cdef inline OrbitPartition *OP_copy(OrbitPartition *OP) noexcept:
    """
    Allocate and return a pointer to a copy of a OrbitPartition of degree n.

    Returns a null pointer in the case of an allocation failure.
    """
    cdef OrbitPartition *OP2 = OP_new(OP.degree)
    if OP is NULL:
        raise MemoryError("MemoryError allocating OrbitPartition in copy method")

    OP_copy_from_to(OP, OP2)
    return OP2

cdef OP_string(OrbitPartition *OP)

cdef inline void OP_clear(OrbitPartition *OP) noexcept:
    cdef int i, n = OP.degree
    for i from 0 <= i < n:
        OP.parent[i] = i
        OP.rank[i] = 0
        OP.mcr[i] = i
        OP.size[i] = 1

cdef inline int OP_find(OrbitPartition *OP, int n) noexcept:
    """
    Report the representative ("root") of the cell which contains n.
    """
    if OP.parent[n] == n:
        return n
    else:
        OP.parent[n] = OP_find(OP, OP.parent[n])
        return OP.parent[n]

cdef inline void OP_join(OrbitPartition *OP, int m, int n) noexcept:
    """
    Join the cells containing m and n, if they are different.
    """
    cdef int m_root = OP_find(OP, m)
    cdef int n_root = OP_find(OP, n)
    if OP.rank[m_root] > OP.rank[n_root]:
        OP.parent[n_root] = m_root
        OP.mcr[m_root] = min(OP.mcr[m_root], OP.mcr[n_root])
        OP.size[m_root] += OP.size[n_root]
    elif OP.rank[m_root] < OP.rank[n_root]:
        OP.parent[m_root] = n_root
        OP.mcr[n_root] = min(OP.mcr[m_root], OP.mcr[n_root])
        OP.size[n_root] += OP.size[m_root]
    elif m_root != n_root:
        OP.parent[n_root] = m_root
        OP.mcr[m_root] = min(OP.mcr[m_root], OP.mcr[n_root])
        OP.size[m_root] += OP.size[n_root]
        OP.rank[m_root] += 1
    if m_root != n_root:
        OP.num_cells -= 1


cdef inline void OP_make_set(OrbitPartition *OP) noexcept:
    cdef int i, n = OP.degree
    cdef int *new_parent, *new_rank, *new_mcr, *new_size

    cdef int *int_array = <int *> sig_malloc( 4*(n+1) * sizeof(int) )
    if not int_array is NULL:
        OP.degree = n + 1
        OP.num_cells = OP.num_cells + 1
        new_parent = int_array
        new_rank = int_array + (n + 1)
        new_mcr = int_array + (2*n + 2)
        new_size = int_array + (3 * n + 3)

        memcpy(new_parent, OP.parent, 4 * n * sizeof(int) )
        memcpy(new_rank, OP.rank, 4 * n * sizeof(int) )
        memcpy(new_mcr, OP.mcr, 4 * n * sizeof(int) )
        memcpy(new_size, OP.size, 4 * n * sizeof(int) )

        new_parent[n] = n
        new_rank[n] = 0
        new_mcr[n] = n
        new_size[n] = 1

        OP.parent = new_parent
        OP.rank   = new_rank
        OP.mcr    = new_mcr
        OP.size   = new_size

cdef inline int OP_merge_list_perm(OrbitPartition *OP, int *gamma) noexcept:
    """
    Joins the cells of OP which intersect the same orbit of gamma.

    INPUT:

    - ``gamma`` - an integer array representing ``i -> gamma[i]``.

    OUTPUT:

    - 1 -- something changed
    - 0 -- orbits of gamma all contained in cells of OP
    """
    cdef int i, i_root, gamma_i_root, changed = 0
    for i from 0 <= i < OP.degree:
        if gamma[i] == i:
            continue
        i_root = OP_find(OP, i)
        gamma_i_root = OP_find(OP, gamma[i])
        if i_root != gamma_i_root:
            changed = 1
            OP_join(OP, i_root, gamma_i_root)
    return changed


# PartitionStack (PS)

cdef inline int PS_copy_from_to(PartitionStack *PS, PartitionStack *PS2) noexcept:
    """
    Copy all data from PS to PS2.
    """
    PS2.depth = PS.depth
    PS2.degree = PS.degree
    memcpy(PS2.entries, PS.entries, 2*PS.degree * sizeof(int) )

cdef inline bint PS_is_discrete(PartitionStack *PS) noexcept:
    """
    Return whether the deepest partition consists only of singleton cells.
    """
    cdef int i
    for i from 0 <= i < PS.degree:
        if PS.levels[i] > PS.depth:
            return 0
    return 1

cdef inline int PS_num_cells(PartitionStack *PS) noexcept:
    """
    Return the number of cells.
    """
    cdef int i, ncells = 0
    for i from 0 <= i < PS.degree:
        if PS.levels[i] <= PS.depth:
            ncells += 1
    return ncells

cdef inline void PS_move_min_to_front(PartitionStack *PS, int start, int end) noexcept:
    """
    Makes sure that the first element of the segment of entries i with
    start <= i <= end is minimal.
    """
    cdef int i, min_loc = start, minimum = PS.entries[start]
    for i from start < i <= end:
        if PS.entries[i] < minimum:
            min_loc = i
            minimum = PS.entries[i]
    if min_loc != start:
        PS.entries[min_loc] = PS.entries[start]
        PS.entries[start] = minimum

cdef inline bint PS_is_mcr(PartitionStack *PS, int m) noexcept:
    """
    Return whether PS.elements[m] (not m!) is the smallest element of its cell.
    """
    return m == 0 or PS.levels[m-1] <= PS.depth

cdef inline bint PS_is_fixed(PartitionStack *PS, int m) noexcept:
    """
    Return whether PS.elements[m] (not m!) is in a singleton cell, assuming
    PS_is_mcr(PS, m) is already true.
    """
    return PS.levels[m] <= PS.depth

cdef inline int PS_clear(PartitionStack *PS) noexcept:
    """
    Set the current partition to the first shallower one, i.e. forgets about
    boundaries between cells that are new to the current level.
    """
    cdef int i, cur_start = 0
    for i from 0 <= i < PS.degree:
        if PS.levels[i] == PS.depth:
            PS.levels[i] += 1
        if PS.levels[i] < PS.depth:
            PS_move_min_to_front(PS, cur_start, i)
            cur_start = i+1

cdef inline int PS_move_all_mins_to_front(PartitionStack *PS) noexcept:
    """
    Move minimal cell elements to the front of each cell.
    """
    cdef int i, cur_start = 0
    for i from 0 <= i < PS.degree:
        if PS.levels[i] <= PS.depth:
            PS_move_min_to_front(PS, cur_start, i)
            cur_start = i+1

cdef inline int PS_get_perm_from(PartitionStack *PS1, PartitionStack *PS2, int *gamma) noexcept:
    """
    Store the permutation determined by PS2[i] -> PS1[i] for each i, where PS[i]
    denotes the entry of the `i`-th cell of the discrete partition PS.
    """
    cdef int i
    for i from 0 <= i < PS1.degree:
        gamma[PS2.entries[i]] = PS1.entries[i]

cdef PartitionStack *PS_new(int n, bint unit_partition) noexcept

cdef PartitionStack *PS_copy(PartitionStack *PS) noexcept

cdef void PS_dealloc(PartitionStack *PS) noexcept

cdef PS_print(PartitionStack *PS)

cdef void PS_unit_partition(PartitionStack *PS) noexcept

cdef int PS_first_smallest(PartitionStack *PS, bitset_t b, int *second_pos=?) noexcept

cdef PartitionStack *PS_from_list(list L) noexcept

cdef list PS_singletons(PartitionStack * part)

cdef int PS_all_new_cells(PartitionStack *PS, bitset_t** nonsingletons_ptr) noexcept

cdef inline bint stacks_are_equivalent(PartitionStack *PS1, PartitionStack *PS2) noexcept:
    cdef int i, j, depth = min(PS1.depth, PS2.depth)
    for i from 0 <= i < PS1.degree:
        if PS1.levels[i] == PS2.levels[i]:
            continue
        elif ((PS1.levels[i] <= depth < PS2.levels[i])
            or (PS2.levels[i] <= depth < PS1.levels[i])):
            return 0
    return 1

cdef int sort_by_function(PartitionStack *PS, int start, int *degrees) noexcept

cdef inline int PS_split_point(PartitionStack *PS, int v) noexcept:
    """
    Detaches the point v from the cell it is in, putting the singleton cell
    of just v in front. Returns the position where v is now located.
    """
    cdef int i = 0, index_of_v
    while PS.entries[i] != v:
        i += 1
    index_of_v = i
    while PS.levels[i] > PS.depth:
        i += 1
    if (index_of_v == 0 or PS.levels[index_of_v-1] <= PS.depth) \
       and PS.levels[index_of_v] > PS.depth:
        # if v is first (make sure v is not already alone)
        PS_move_min_to_front(PS, index_of_v+1, i)
        PS.levels[index_of_v] = PS.depth
        return index_of_v
    else:
        # If v is not at front, i.e. v is not minimal in its cell,
        # then move_min_to_front is not necessary since v will swap
        # with the first before making its own cell, leaving it at
        # the front of the other.
        i = index_of_v
        while i != 0 and PS.levels[i-1] > PS.depth:
            i -= 1
        PS.entries[index_of_v] = PS.entries[i+1] # move the second element to v
        PS.entries[i+1] = PS.entries[i] # move the first (min) to second
        PS.entries[i] = v # place v first
        PS.levels[i] = PS.depth
        return i

cdef inline int split_point_and_refine(PartitionStack *PS, int v, void *S,
    int (*refine_and_return_invariant)
         (PartitionStack *PS, void *S, int *cells_to_refine_by, int ctrb_len) noexcept,
    int *cells_to_refine_by) noexcept:
    """
    Make the partition stack one longer by copying the last partition in the
    stack, split off a given point, and refine. Return the invariant given by
    the refinement function.

    INPUT:
    - ``PS`` -- the partition stack to refine
    - ``v`` -- the point to split
    - ``S`` -- the structure
    - ``refine_and_return_invariant`` -- the refinement function provided
    - ``cells_to_refine_by`` -- an array, contents ignored
    - ``group`` -- the containing group, NULL for full S_n
    - ``perm_stack`` -- represents a partial traversal decomposition for group
    """
    PS.depth += 1
    PS_clear(PS)
    cells_to_refine_by[0] = PS_split_point(PS, v)
    return refine_and_return_invariant(PS, S, cells_to_refine_by, 1)


# StabilizerChain (SC)

cdef StabilizerChain *SC_new(int n, bint init_gens=?) noexcept

cdef int SC_realloc_gens(StabilizerChain *SC, int level, int size) noexcept

cdef void SC_dealloc(StabilizerChain *SC) noexcept

cdef int SC_copy_nomalloc(StabilizerChain *SC_dest, StabilizerChain *SC,
        int level) noexcept

cdef StabilizerChain *SC_alternating_group(int n) noexcept

cdef int SC_insert_and_sift(StabilizerChain *SC, int level, int *pi,
        int num_perms, bint sift) noexcept

cdef int SC_insert_base_point_nomalloc(StabilizerChain *SC_dest,
        StabilizerChain *SC, int level, int p) noexcept

cdef inline int SC_perm_is_identity(int *perm, int degree) noexcept:
    for i from 0 <= i < degree:
        if perm[i] != i:
            break
    else:
        return 1
    return 0

cdef inline void SC_mult_perms(int *out, int *first, int *second, int degree) noexcept:
    """
    DON'T DO THIS WITH out == second!
    """
    cdef int i
    for i from 0 <= i < degree:
        out[i] = second[first[i]]

cdef inline void SC_invert_perm(int *out, int *input, int degree) noexcept:
    """
    DON'T DO THIS WITH out == in!
    """
    cdef int i
    for i from 0 <= i < degree:
        out[input[i]] = i

cdef inline void SC_identify(int *perm, int degree) noexcept:
    cdef int i
    for i from 0 <= i < degree:
        perm[i] = i

cdef inline void SC_add_base_point(StabilizerChain *SC, int b) noexcept:
    """
    Add base point b to the end of SC. Assumes b is not already in the base.
    """
    cdef int i, n = SC.degree
    SC.orbit_sizes[SC.base_size] = 1
    SC.num_gens[SC.base_size] = 0
    SC.base_orbits[SC.base_size][0] = b
    for i from 0 <= i < n:
        SC.parents[SC.base_size][i] = -1
    SC.parents[SC.base_size][b] = b
    SC.labels[SC.base_size][b] = 0
    SC.base_size += 1

cdef inline int SC_cleanup(StabilizerChain *SC) noexcept:
    """
    Remove redundant base elements from SC.

    Returns 1 if nothing changed, and 2 in case of an allocation failure.
    """
    cdef int old, new = 0, i, n = SC.degree
    for old from 0 <= old < SC.base_size:
        if SC.orbit_sizes[old] != 1:
            if old != new:
                # copy row old to row new
                SC.orbit_sizes[new] = SC.orbit_sizes[old]
                SC.num_gens[new] = SC.num_gens[old]
                if SC.array_size[new] < SC.array_size[old]:
                    if SC_realloc_gens(SC, new, max(SC.array_size[old], 2*SC.array_size[new])):
                        return 2
                memcpy(SC.base_orbits[new],  SC.base_orbits[old],    n * sizeof(int))
                memcpy(SC.parents[new],      SC.parents[old],        n * sizeof(int))
                memcpy(SC.labels[new],       SC.labels[old],         n * sizeof(int))
                memcpy(SC.generators[new],   SC.generators[old],     n*SC.num_gens[old] * sizeof(int))
                memcpy(SC.gen_inverses[new], SC.gen_inverses[old],   n*SC.num_gens[old] * sizeof(int))
            new += 1
    old = SC.base_size
    SC.base_size = new
    return (old == new)

cdef inline void SC_compose_up_to_base(StabilizerChain *SC, int level, int x, int *perm) noexcept:
    """
    Repeatedly compose the given perm by labels on the Schreier tree, starting
    with x, until the base is reached. The composition is stored to perm.
    """
    cdef int b = SC.base_orbits[level][0], n = SC.degree
    cdef int *label
    cdef int label_no
    while x != b:
        label_no = SC.labels[level][x]
        if label_no < 0:
            label_no = -label_no - 1
            label = SC.gen_inverses[level] + n*label_no
        else:
            label_no = label_no - 1
            label = SC.generators[level] + n*label_no
        x = SC.parents[level][x]
        SC_mult_perms(perm, perm, label, n)

cdef inline void SC_scan(StabilizerChain *SC, int level, int x, int gen_index, int *gen, int sign) noexcept:
    """
    See whether the point x is moved to a point outside the
    tree by gen, and if so add it to the tree (arc label is gen_inv).

    gen_index - where in the generator array the generator is located
    gen - points to the generator
    gen_inv - points to the inverse
    sign - whether to take SC.generators or SC.gen_inverses to go *up* the tree
    """
    cdef int y = gen[x], n = SC.degree
    if SC.parents[level][y] == -1:
        SC.base_orbits[level][SC.orbit_sizes[level]] = y
        SC.orbit_sizes[level] += 1
        SC.parents[level][y] = x
        SC.labels[level][y] = sign*(gen_index+1)

cdef inline int SC_insert(StabilizerChain *SC, int level, int *pi, int num_perms) noexcept:
    """
    Add permutations in pi to the stabilizer chain. The array pi is a sequence
    of num_perms permutations, each in list representation, hence pi should be
    at least length SC.degree*num_perms. There must be at most SC.degree perms.
    (Simply call the function again if you want to add more.)

    The variable ``level`` is used for recursion. From the outside, should be
    set to zero. On the inside, used to bring the data structure up to date of
    level ``level``, given that it is up to date on ``level + 1``.

    Return values:
    0 - No errors.
    1 - Allocation failure.
    """
    return SC_insert_and_sift(SC, level, pi, num_perms, 1)

cdef inline int SC_update_tree(StabilizerChain *SC, int level, int *pi, int num_perms) noexcept:
    return SC_insert_and_sift(SC, level, pi, num_perms, 0)

cdef inline void SC_order(StabilizerChain *SC, int i, mpz_t order) noexcept:
    """
    Give the order of the stabilizer of base points up to but not including the
    `i`-th, storing it to ``order``, which must be already initialized.

    To get the order of the full group, let ``i = 0``.
    """
    cdef int k
    mpz_set_si(order, 1)
    for k from i <= k < SC.base_size:
        mpz_mul_si(order, order, SC.orbit_sizes[k])

cdef inline bint SC_contains(StabilizerChain *SC, int level, int *pi, bint modify) noexcept:
    """
    Test whether pi is in the level-th stabilizer.

    Assumes that pi stabilizes the first level base points.
    """
    cdef int b, i, j, x, n = SC.degree
    cdef int *perm
    if modify:
        perm = pi
    else:
        perm = SC.perm_scratch
        memcpy(perm, pi, n*sizeof(int))
    for i from level <= i < SC.base_size:
        b = SC.base_orbits[i][0]
        x = perm[b]
        if x == b:
            continue
        if SC.parents[i][x] == -1:
            return 0
        SC_compose_up_to_base(SC, i, x, perm)
    return SC_perm_is_identity(perm, n)

cdef inline void SC_random_element(StabilizerChain *SC, int level, int *perm) noexcept:
    """
    Give a random element of the level-th stabilizer. For a random element of the
    whole group, set level to 0. Must have level < SC.base_size.
    """
    cdef int i, x, n = SC.degree
    SC_identify(perm, n)
    for i from level <= i < SC.base_size:
        x = SC.base_orbits[i][rand()%SC.orbit_sizes[i]]
        SC_compose_up_to_base(SC, i, x, perm)

cdef int compute_relabeling(StabilizerChain *group,
        StabilizerChain *scratch_group,
        int *permutation, int *relabeling) noexcept

cdef inline void update_perm_stack(StabilizerChain *group, int level, int point,
    int *perm_stack) noexcept:
    """
    Ensure that perm_stack[level] is gamma_0^{-1}...gamma_{level-1}^{-1}, where
    each gamma_i represents the coset representative at the `i`-th level determined
    by our current position in the search tree.

    For internal use within the automorphism group, canonical label and double
    coset functions, to be called after refinement (level = depth after refinement).
    """
    cdef int n = group.degree
    memcpy(perm_stack + n*level, perm_stack + n*(level-1), n*sizeof(int))
    SC_compose_up_to_base(group, level-1, perm_stack[n*(level-1) + point], perm_stack + n*level)

cdef inline int split_point_and_refine_by_orbits(PartitionStack *PS, int v,
    void *S, int (*refine_and_return_invariant)
         (PartitionStack *PS, void *S, int *cells_to_refine_by, int ctrb_len) noexcept,
    int *cells_to_refine_by, StabilizerChain *SC, int *perm_stack) noexcept:
    """ """
    PS.depth += 1
    PS_clear(PS)
    cells_to_refine_by[0] = PS_split_point(PS, v)
    update_perm_stack(SC, PS.depth, v, perm_stack)
    return refine_also_by_orbits(PS, S, refine_and_return_invariant, cells_to_refine_by, 1, SC, perm_stack)

cdef int refine_by_orbits(PartitionStack *PS, StabilizerChain *SC,
        int *perm_stack, int *cells_to_refine_by, int *ctrb_len) noexcept

cdef inline int refine_also_by_orbits(PartitionStack *PS, void *S,
    int (*refine_and_return_invariant)
         (PartitionStack *PS, void *S, int *cells_to_refine_by, int ctrb_len) noexcept,
    int *cells_to_refine_by, int ctrb_len, StabilizerChain *SC, int *perm_stack) noexcept:
    """ """
    cdef int inv
    inv = refine_by_orbits(PS, SC, perm_stack, cells_to_refine_by, &ctrb_len)
    inv += refine_and_return_invariant(PS, S, cells_to_refine_by, ctrb_len)
    return inv
