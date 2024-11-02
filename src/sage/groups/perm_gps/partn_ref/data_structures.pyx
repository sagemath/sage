r"""
Data structures

This module implements basic data structures essential to the rest of the
partn_ref module.

REFERENCES:

    [1] McKay, Brendan D. Practical Graph Isomorphism. Congressus Numerantium,
        Vol. 30 (1981), pp. 45-87.
    [2] Fredman, M. and Saks, M. The cell probe complexity of dynamic data
        structures. Proceedings of the Twenty-First Annual ACM Symposium on
        Theory of Computing, pp. 345–354. May 1989.
    [3] Seress, Akos. Permutation Group Algorithms. Cambridge University Press,
        2003.
"""

#*****************************************************************************
#      Copyright (C) 2006 - 2011 Robert L. Miller <rlmillster@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from libc.math cimport log, ceil
from libc.string cimport memcpy, memset
from libc.stdlib cimport rand
from cysignals.memory cimport sig_malloc, sig_calloc, sig_realloc, sig_free

from sage.data_structures.bitset_base cimport *
from sage.rings.integer cimport Integer
# from sage.libs.flint.ulong_extras cimport n_is_prime
# -- avoid modularization obstruction -- function is only used for a doctest helper
from sage.arith.misc import is_prime as n_is_prime


# OrbitPartition (OP)

cdef inline OrbitPartition *OP_new(int n) noexcept:
    """
    Allocate and return a pointer to a new OrbitPartition of degree n. Return a
    null pointer in the case of an allocation failure.
    """
    cdef OrbitPartition *OP = <OrbitPartition *> \
                                sig_malloc(sizeof(OrbitPartition))
    cdef int *int_array = <int *> sig_malloc( 4*n * sizeof(int) )
    if OP is NULL or int_array is NULL:
        sig_free(OP)
        sig_free(int_array)
        return NULL
    OP.degree = n
    OP.num_cells = n
    OP.parent = int_array
    OP.rank   = int_array +   n
    OP.mcr    = int_array + 2*n
    OP.size   = int_array + 3*n
    OP_clear(OP)
    return OP

cdef inline void OP_dealloc(OrbitPartition *OP) noexcept:
    if OP is not NULL:
        sig_free(OP.parent)
        sig_free(OP)

cdef OP_string(OrbitPartition *OP):
    """
    Return a string representation of the OrbitPartition.
    """
    cdef i,j
    s = ""
    for i in range(OP.degree):
        s += " "
        j = OP_find(OP, i)
        s += "%d -> %d" % (i, j)
    return s


def OP_represent(int n, merges, perm):
    """
    Demonstration and testing.

    TESTS::

        sage: from sage.groups.perm_gps.partn_ref.data_structures import OP_represent
        sage: OP_represent(9, [(0,1),(2,3),(3,4)], [1,2,0,4,3,6,7,5,8])
        Allocating OrbitPartition...
        Allocation passed.
        Checking that each element reports itself as its root.
        Each element reports itself as its root.
        Merging:
        Merged 0 and 1.
        Merged 2 and 3.
        Merged 3 and 4.
        Done merging.
        Finding:
        0 -> 0, root: size=2, mcr=0, rank=1
        1 -> 0
        2 -> 2, root: size=3, mcr=2, rank=1
        3 -> 2
        4 -> 2
        5 -> 5, root: size=1, mcr=5, rank=0
        6 -> 6, root: size=1, mcr=6, rank=0
        7 -> 7, root: size=1, mcr=7, rank=0
        8 -> 8, root: size=1, mcr=8, rank=0
        Allocating array to test merge_perm.
        Allocation passed.
        Merging permutation: [1, 2, 0, 4, 3, 6, 7, 5, 8]
        Done merging.
        Finding:
        0 -> 0, root: size=5, mcr=0, rank=2
        1 -> 0
        2 -> 0
        3 -> 0
        4 -> 0
        5 -> 5, root: size=3, mcr=5, rank=1
        6 -> 5
        7 -> 5
        8 -> 8, root: size=1, mcr=8, rank=0
        Deallocating OrbitPartition.
        Done.
    """
    cdef int i
    print("Allocating OrbitPartition...")
    cdef OrbitPartition *OP = OP_new(n)
    if OP is NULL:
        print("Allocation failed!")
        return
    print("Allocation passed.")
    print("Checking that each element reports itself as its root.")
    good = True
    for i in range(n):
        if not OP_find(OP, i) == i:
            print("Failed at i = %d!" % i)
            good = False
    if good:
        print("Each element reports itself as its root.")
    print("Merging:")
    for i,j in merges:
        OP_join(OP, i, j)
        print("Merged %d and %d." % (i, j))
    print("Done merging.")
    print("Finding:")
    for i in range(n):
        j = OP_find(OP, i)
        s = "%d -> %d"%(i, j)
        if i == j:
            s += ", root: size=%d, mcr=%d, rank=%d"%\
                   (OP.size[i], OP.mcr[i], OP.rank[i])
        print(s)
    print("Allocating array to test merge_perm.")
    cdef int *gamma = <int *> sig_malloc( n * sizeof(int) )
    if gamma is NULL:
        print("Allocation failed!")
        OP_dealloc(OP)
        return
    print("Allocation passed.")
    for i in range(n):
        gamma[i] = perm[i]
    print("Merging permutation: %s" % perm)
    OP_merge_list_perm(OP, gamma)
    print("Done merging.")
    print("Finding:")
    for i in range(n):
        j = OP_find(OP, i)
        s = "%d -> %d"%(i, j)
        if i == j:
            s += ", root: size=%d, mcr=%d, rank=%d"%\
                   (OP.size[i], OP.mcr[i], OP.rank[i])
        print(s)
    print("Deallocating OrbitPartition.")
    sig_free(gamma)
    OP_dealloc(OP)
    print("Done.")


# PartitionStack (PS)

cdef inline PartitionStack *PS_new(int n, bint unit_partition) noexcept:
    """
    Allocate and return a pointer to a new PartitionStack of degree n. Return a
    null pointer in the case of an allocation failure.
    """
    cdef PartitionStack *PS = <PartitionStack *> \
                                sig_malloc(sizeof(PartitionStack))
    cdef int *int_array = <int *> sig_malloc( 2*n * sizeof(int) )
    if PS is NULL or int_array is NULL:
        sig_free(PS)
        sig_free(int_array)
        return NULL
    PS.entries = int_array
    PS.levels  = int_array + n
    PS.depth  = 0
    PS.degree = n
    if unit_partition:
        PS_unit_partition(PS)
    return PS

cdef void PS_unit_partition(PartitionStack *PS) noexcept:
    """
    Set partition stack to a single partition with a single cell.
    """
    cdef int i, n = PS.degree
    PS.depth = 0
    for i in range(n - 1):
        PS.entries[i] = i
        PS.levels[i] = n
    PS.entries[n-1] = n - 1
    PS.levels[n-1] = -1

cdef inline PartitionStack *PS_copy(PartitionStack *PS) noexcept:
    """
    Allocate and return a pointer to a copy of PartitionStack PS. Return a null
    pointer in the case of an allocation failure.
    """
    cdef int n = PS.degree

    cdef PartitionStack *PS2 = <PartitionStack *> \
                                sig_malloc(sizeof(PartitionStack))
    cdef int *int_array = <int *> sig_malloc( 2*n * sizeof(int) )
    if PS2 is NULL or int_array is NULL:
        sig_free(PS2)
        sig_free(int_array)
        return NULL
    PS2.entries = int_array
    PS2.levels  = int_array + n
    PS_copy_from_to(PS, PS2)
    return PS2

cdef inline void PS_dealloc(PartitionStack *PS) noexcept:
    if PS is not NULL:
        sig_free(PS.entries)
        sig_free(PS)

cdef PartitionStack *PS_from_list(list L) noexcept:
    """
    Allocate and return a pointer to a PartitionStack representing L. Return a
    null pointer in the case of an allocation failure.
    """
    cdef int cell, i, num_cells = len(L), cur_start = 0, cur_len, n = 0
    for cell in range(num_cells):
        n += len(L[cell])
    cdef PartitionStack *PS = PS_new(n, 0)
    if PS is NULL:
        return NULL
    for cell in range(num_cells):
        cur_len = len(L[cell])
        for i in range(cur_len):
            PS.entries[cur_start + i] = L[cell][i]
            PS.levels[cur_start + i] = n
        PS_move_min_to_front(PS, cur_start, cur_start+cur_len-1)
        cur_start += cur_len
        PS.levels[cur_start-1] = 0
    if n > 0:
        PS.levels[n-1] = -1
    PS.depth = 0
    PS.degree = n
    return PS

cdef PS_print(PartitionStack *PS):
    """
    Print a visual representation of PS.
    """
    cdef int i
    for i in range(PS.depth + 1):
        PS_print_partition(PS, i)

cdef PS_print_partition(PartitionStack *PS, int k):
    """
    Print the partition at depth k.
    """
    s = '('
    for i in range(PS.degree):
        s += str(PS.entries[i])
        if PS.levels[i] <= k:
            s += '|'
        else:
            s += ' '
    s = s[:-1] + ')'
    print(s)

cdef int PS_first_smallest(PartitionStack *PS, bitset_t b, int *second_pos=NULL) noexcept:
    """
    Find the first occurrence of the smallest cell of size greater than one,
    which is admissible (checked by the function ``test_allowance``).
    Its entries are stored to `b` and its minimum element is returned.
    """
    cdef int i = 0, j = 0, location = 0, n = PS.degree
    bitset_zero(b)
    while True:
        if PS.levels[i] <= PS.depth:
            if i != j and n > i - j + 1:
                n = i - j + 1
                location = j
            j = i + 1
        if PS.levels[i] == -1:
            break
        i += 1
    # location now points to the beginning of the first, smallest,
    # nontrivial cell
    i = location
    while True:
        bitset_flip(b, PS.entries[i])
        if PS.levels[i] <= PS.depth:
            break
        i += 1

    if second_pos != NULL:
        if n == 2:
            second_pos[0] = PS.entries[location + 1]
        else:
            second_pos[0] = -1

    return PS.entries[location]


cdef int PS_all_new_cells(PartitionStack *PS, bitset_t** nonsingletons_ptr) noexcept:
    """
    Suppose a cell ``C`` was split into ``a`` components at ``PS.level``.
    Set the rows of the matrix ``nonsingletons_ptr`` to the first
    ``a-1`` components of ``C`` for all those ``C``.
    Return the number of rows of ``nonsingletons_ptr``.
    """
    cdef int beg=0, end, n = PS.degree, count=0, i, n_1 = n-1
    cdef bitset_t scratch
    bitset_init(scratch, n)
    cdef bitset_t* nonsingletons = nonsingletons_ptr[0]

    while beg < n_1:
        end = beg
        while end!=n and PS.levels[end] > PS.depth:
            end+=1

        if end != n:
            if PS.levels[end] == PS.depth:
                bitset_zero(scratch)
                for i in range(beg, end + 1):
                    bitset_set(scratch, PS.entries[i])
                count +=1
                nonsingletons = <bitset_t*> sig_realloc(nonsingletons, count * sizeof(bitset_t))
                if nonsingletons is NULL:
                    raise MemoryError("Memory error in PS_all_new_cells")
                bitset_init(nonsingletons[count-1], n)
                bitset_copy(nonsingletons[count-1], scratch)
        else:
            if beg==0:
                nonsingletons = <bitset_t*> sig_realloc(nonsingletons, sizeof(bitset_t))
                if nonsingletons is NULL:
                    raise MemoryError("Memory error in PS_all_new_cells")
                bitset_init(nonsingletons[0], n)
                bitset_zero(scratch)
                bitset_complement(nonsingletons[0], scratch)
                count = 1
        beg = end+1
    nonsingletons_ptr[0] = nonsingletons
    return count

cdef int PS_find_element(PartitionStack *PS, bitset_t b, int x) except -1:
    """
    Find the cell containing x, store its entries to b and return the location
    of the beginning of the cell.
    """
    cdef int i, location, n = PS.degree
    bitset_zero(b)
    for i in range(n):
        if PS.entries[i] == x:
            location = i
            break
    else:
        raise ValueError("element not found")
    while location > 0 and PS.levels[location-1] > PS.depth:
        location -= 1
    i = 0
    while 1:
        bitset_set(b, PS.entries[location+i])
        if PS.levels[location+i] <= PS.depth:
            break
        i += 1
    return location

cdef list PS_singletons(PartitionStack * part):
    """
    Return the list of all singletons in the ``PartitionStack``.
    """
    cdef list l = []
    cdef int i

    if part.levels[0] <= part.depth:
        l.append(0)

    for i in range(1, part.degree):
        if part.levels[i] <= part.depth and part.levels[i - 1] <= part.depth:
            l.append(i)

    return l


def PS_represent(partition, splits):
    """
    Demonstration and testing.

    TESTS::

        sage: from sage.groups.perm_gps.partn_ref.data_structures import PS_represent
        sage: PS_represent([[6],[3,4,8,7],[1,9,5],[0,2]], [6,1,8,2])
        Allocating PartitionStack...
        Allocation passed:
        (0 1 2 3 4 5 6 7 8 9)
        Checking that entries are in order and correct level.
        Everything seems in order, deallocating.
        Deallocated.
        Creating PartitionStack from partition [[6], [3, 4, 8, 7], [1, 9, 5], [0, 2]].
        PartitionStack's data:
        entries -> [6, 3, 4, 8, 7, 1, 9, 5, 0, 2]
        levels -> [0, 10, 10, 10, 0, 10, 10, 0, 10, -1]
        depth = 0, degree = 10
        (6|3 4 8 7|1 9 5|0 2)
        Checking PS_is_discrete:
        False
        Checking PS_num_cells:
        4
        Checking PS_is_mcr, min cell reps are:
        [6, 3, 1, 0]
        Checking PS_is_fixed, fixed elements are:
        [6]
        Copying PartitionStack:
        (6|3 4 8 7|1 9 5|0 2)
        Checking for consistency.
        Everything is consistent.
        Clearing copy:
        (0 3 4 8 7 1 9 5 6 2)
        Splitting point 6 from original:
        0
        (6|3 4 8 7|1 9 5|0 2)
        Splitting point 1 from original:
        5
        (6|3 4 8 7|1|5 9|0 2)
        Splitting point 8 from original:
        1
        (6|8|3 4 7|1|5 9|0 2)
        Splitting point 2 from original:
        8
        (6|8|3 4 7|1|5 9|2|0)
        Getting permutation from PS2->PS:
        [6, 1, 0, 8, 3, 9, 2, 7, 4, 5]
        Finding first smallest:
        Minimal element is 5, bitset is:
        0000010001
        Finding element 1:
        Location is: 5
        Bitset is:
        0100000000
        Deallocating PartitionStacks.
        Done.
    """
    cdef int i, n = sum([len(cell) for cell in partition])
    cdef int *gamma
    cdef bitset_t b
    print("Allocating PartitionStack...")
    cdef PartitionStack *PS = PS_new(n, 1)
    cdef PartitionStack *PS2
    if PS is NULL:
        print("Allocation failed!")
        return
    print("Allocation passed:")
    PS_print(PS)
    print("Checking that entries are in order and correct level.")
    good = True
    for i in range(n - 1):
        if not (PS.entries[i] == i and PS.levels[i] == n):
            print("Failed at i = %d!" % i)
            print(PS.entries[i], PS.levels[i], i, n)
            good = False
    if not (PS.entries[n-1] == n-1 and PS.levels[n-1] == -1):
        print("Failed at i = %d!"%(n-1))
        good = False
    if not PS.degree == n or not PS.depth == 0:
        print("Incorrect degree or depth!")
        good = False
    if good:
        print("Everything seems in order, deallocating.")
    PS_dealloc(PS)
    print("Deallocated.")
    print("Creating PartitionStack from partition %s." % partition)
    PS = PS_from_list(partition)
    print("PartitionStack's data:")
    print("entries -> %s"%[PS.entries[i] for i in range(n)])
    print("levels -> %s"%[PS.levels[i] for i in range(n)])
    print("depth = %d, degree = %d" % (PS.depth,PS.degree))
    PS_print(PS)
    print("Checking PS_is_discrete:")
    print("True" if PS_is_discrete(PS) else "False")
    print("Checking PS_num_cells:")
    print(PS_num_cells(PS))
    print("Checking PS_is_mcr, min cell reps are:")
    L = [PS.entries[i] for i in range(n) if PS_is_mcr(PS, i)]
    print(L)
    print("Checking PS_is_fixed, fixed elements are:")
    print([PS.entries[l] for l in L if PS_is_fixed(PS, l)])
    print("Copying PartitionStack:")
    PS2 = PS_copy(PS)
    PS_print(PS2)
    print("Checking for consistency.")
    good = True
    for i in range(n):
        if PS.entries[i] != PS2.entries[i] or PS.levels[i] != PS2.levels[i]:
            print("Failed at i = %d!"%i)
            good = False
    if PS.degree != PS2.degree or PS.depth != PS2.depth:
        print("Failure with degree or depth!")
        good = False
    if good:
        print("Everything is consistent.")
    print("Clearing copy:")
    PS_clear(PS2)
    PS_print(PS2)
    for s in splits:
        print("Splitting point %d from original:"%s)
        print(PS_split_point(PS, s))
        PS_print(PS)
    print("Getting permutation from PS2->PS:")
    gamma = <int *> sig_malloc(n * sizeof(int))
    PS_get_perm_from(PS, PS2, gamma)
    print([gamma[i] for i in range(n)])
    sig_free(gamma)
    print("Finding first smallest:")
    bitset_init(b, n)
    i = PS_first_smallest(PS, b)
    print("Minimal element is %d, bitset is:"%i)
    print(bitset_string(b))
    bitset_free(b)
    print("Finding element 1:")
    bitset_init(b, n)
    print("Location is: {}".format(PS_find_element(PS, b, 1)))
    print("Bitset is:")
    print(bitset_string(b))
    bitset_free(b)
    print("Deallocating PartitionStacks.")
    PS_dealloc(PS)
    PS_dealloc(PS2)
    print("Done.")


# StabilizerChain (SC)

cdef enum:
    default_num_gens = 8
    default_num_bits = 64

cdef StabilizerChain *SC_new(int n, bint init_gens=True) noexcept:
    """
    Allocate and return a pointer to a new StabilizerChain of degree n. Return
    a null pointer in the case of an allocation failure.
    """
    cdef int i
    cdef StabilizerChain *SC = <StabilizerChain *> \
                                sig_calloc(1, sizeof(StabilizerChain))
    if SC is NULL:
        return NULL
    SC.degree = n
    SC.base_size = 0
    if n == 0:
        # All internal pointers have been initialized to NULL by sig_calloc
        return SC

    # first level allocations
    cdef int *int_array = <int *>  sig_malloc( (3*n*n + 6*n + 1) * sizeof(int) )
    cdef int **int_ptrs = <int **> sig_calloc( 5*n, sizeof(int *) )
    SC.OP_scratch = OP_new(n)
    # bitset_init without the MemoryError:
    cdef long limbs = (default_num_bits - 1)/(8*sizeof(unsigned long)) + 1
    SC.gen_used.size   = default_num_bits
    SC.gen_is_id.size  = default_num_bits
    SC.gen_used.limbs  = limbs
    SC.gen_is_id.limbs = limbs
    SC.gen_used.bits   = <mp_limb_t*>sig_malloc(limbs * sizeof(mp_limb_t))
    SC.gen_is_id.bits  = <mp_limb_t*>sig_malloc(limbs * sizeof(mp_limb_t))

    # check for allocation failures
    if int_array        is NULL or int_ptrs          is NULL or \
       SC.gen_used.bits is NULL or SC.gen_is_id.bits is NULL or \
       SC.OP_scratch    is NULL:
        sig_free(int_array)
        sig_free(int_ptrs)
        SC_dealloc(SC)
        return NULL

    SC.gen_used.bits[limbs-1] = 0
    SC.gen_is_id.bits[limbs-1] = 0

    SC.orbit_sizes  = int_array
    SC.num_gens     = int_array +   n
    SC.array_size   = int_array + 2*n
    SC.perm_scratch = int_array + 3*n # perm_scratch is length 3*n+1 for sorting
    int_array += 6*n + 1

    SC.generators   = int_ptrs
    SC.gen_inverses = int_ptrs +   n
    SC.base_orbits  = int_ptrs + 2*n
    SC.parents      = int_ptrs + 3*n
    SC.labels       = int_ptrs + 4*n
    for i in range(n):
        SC.base_orbits[i] = int_array
        SC.parents[i]     = int_array +   n
        SC.labels[i]      = int_array + 2*n
        int_array += 3*n

    # second level allocations
    if init_gens:
        for i in range(n):
            SC.array_size[i] = default_num_gens
            SC.generators[i]   = <int *> sig_malloc( default_num_gens*n * sizeof(int) )
            SC.gen_inverses[i] = <int *> sig_malloc( default_num_gens*n * sizeof(int) )
            if SC.generators[i] is NULL or SC.gen_inverses[i] is NULL:
                SC_dealloc(SC)
                return NULL

    return SC

cdef inline int SC_realloc_gens(StabilizerChain *SC, int level, int size) noexcept:
    """
    Reallocate generator array at level `level` to size `size`.

    Return 1 in case of an allocation failure.
    """
    cdef int *temp
    cdef int n = SC.degree

    temp = <int *> sig_realloc(SC.generators[level], n * size * sizeof(int))
    if temp is NULL:
        return 1
    SC.generators[level] = temp

    temp = <int *> sig_realloc(SC.gen_inverses[level], n * size * sizeof(int))
    if temp is NULL:
        return 1
    SC.gen_inverses[level] = temp

    SC.array_size[level] = size
    return 0

cdef inline void SC_dealloc(StabilizerChain *SC) noexcept:
    cdef int i, n
    if SC is not NULL:
        n =  SC.degree
        if SC.generators is not NULL:
            for i in range(n):
                sig_free(SC.generators[i])
                sig_free(SC.gen_inverses[i])
        sig_free(SC.generators) # frees int_ptrs
        sig_free(SC.orbit_sizes) # frees int_array
        sig_free(SC.gen_used.bits)
        sig_free(SC.gen_is_id.bits)
        OP_dealloc(SC.OP_scratch)
        sig_free(SC)

cdef StabilizerChain *SC_symmetric_group(int n) noexcept:
    """
    Return a stabilizer chain for the symmetric group on {0, 1, ..., n-1}.

    Return ``NULL`` in the case of an allocation failure.
    """
    cdef int i, j, b
    cdef StabilizerChain *SC = SC_new(n, False)
    if SC is NULL:
        return NULL
    SC.base_size = n-1
    for i in range(n - 1):
        SC.array_size[i] = n-i-1
    SC.array_size[n-1] = default_num_gens
    for i in range(n):
        SC.generators[i]   = <int *> sig_malloc( SC.array_size[i]*n * sizeof(int) )
        SC.gen_inverses[i] = <int *> sig_malloc( SC.array_size[i]*n * sizeof(int) )
        if SC.generators[i] is NULL or SC.gen_inverses[i] is NULL:
            SC_dealloc(SC)
            return NULL
    cdef int *id_perm = SC.perm_scratch
    for i in range(n):
        id_perm[i] = i
    for i in range(n - 1):
        b = i
        SC.orbit_sizes[i] = n-i
        SC.num_gens[i] = n-i-1
        for j in range(i):
            SC.parents[i][j] = -1
        for j in range(n - i):
            SC.base_orbits[i][j] = i+j
            SC.parents[i][i+j] = b
            SC.labels[i][i+j] = j
        for j in range(n - i - 1):
            #j-th generator sends i+j+1 to b
            memcpy(SC.generators[i] + n*j, id_perm, n * sizeof(int) )
            SC.generators[i][n*j + i+j+1] = b
            SC.generators[i][n*j + b] = i+j+1
            memcpy(SC.gen_inverses[i] + n*j, SC.generators[i] + n*j, n * sizeof(int) )
    return SC

cdef StabilizerChain *SC_alternating_group(int n) noexcept:
    """
    Return a stabilizer chain for the alternating group on {0, 1, ..., n-1}.

    Return ``NULL`` in the case of an allocation failure.
    """
    cdef int i, j, b
    cdef StabilizerChain *SC = SC_new(n, False)
    if SC is NULL:
        return NULL
    SC.base_size = n-2
    for i in range(n - 2):
        SC.array_size[i] = n-i-1
    SC.array_size[n-2] = default_num_gens
    SC.array_size[n-1] = default_num_gens
    for i in range(n):
        SC.generators[i]   = <int *> sig_malloc( SC.array_size[i]*n * sizeof(int) )
        SC.gen_inverses[i] = <int *> sig_malloc( SC.array_size[i]*n * sizeof(int) )
        if SC.generators[i] is NULL or SC.gen_inverses[i] is NULL:
            SC_dealloc(SC)
            return NULL
    cdef int *id_perm = SC.perm_scratch
    for i in range(n):
        id_perm[i] = i
    for i in range(n - 2):
        b = i
        SC.orbit_sizes[i] = n-i
        SC.num_gens[i] = n-i-2
        for j in range(i):
            SC.parents[i][j] = -1
        for j in range(n - i):
            SC.base_orbits[i][j] = i+j
            SC.parents[i][i+j] = b
            SC.labels[i][i+j] = j
        SC.labels[i][n-1] = -(n-i-2)
        for j in range(n - i - 2):
            #j-th generator sends i+j+1 to b, i+j+2 to i+j+1, and b to i+j+2
            memcpy(SC.generators[i] + n*j, id_perm, n * sizeof(int) )
            SC.generators[i][n*j + i+j+1] = b
            SC.generators[i][n*j + b    ] = i+j+2
            SC.generators[i][n*j + i+j+2] = i+j+1
            SC_invert_perm(SC.gen_inverses[i] + n*j, SC.generators[i] + n*j, n)
    return SC

cdef int SC_realloc_bitsets(StabilizerChain *SC, unsigned long size) noexcept:
    """
    If size is larger than current allocation, double the size of the bitsets
    until it is not.

    Return 1 in case of an allocation failure.
    """
    cdef unsigned long size_old = SC.gen_used.size
    if size <= size_old:
        return 0
    cdef unsigned long new_size = size_old
    while new_size < size:
        new_size *= 2
    cdef long limbs = (new_size - 1)/(8*sizeof(unsigned long)) + 1
    cdef mp_limb_t *tmp = <mp_limb_t*> sig_realloc(SC.gen_used.bits, limbs * sizeof(mp_limb_t))
    if tmp is not NULL:
        SC.gen_used.bits = tmp
    else:
        return 1
    tmp = <mp_limb_t*> sig_realloc(SC.gen_is_id.bits, limbs * sizeof(mp_limb_t))
    if tmp is not NULL:
        SC.gen_is_id.bits = tmp
    else:
        return 1
    SC.gen_used.limbs = limbs
    SC.gen_is_id.limbs = limbs
    SC.gen_used.size = new_size
    SC.gen_is_id.size = new_size
    SC.gen_used.bits[size_old >> index_shift] &= limb_lower_bits_down(size_old)
    memset(SC.gen_used.bits + (size_old >> index_shift) + 1, 0, (limbs - (size_old >> index_shift) - 1) * sizeof(unsigned long))
    SC.gen_is_id.bits[size_old >> index_shift] &= limb_lower_bits_down(size_old)
    memset(SC.gen_is_id.bits + (size_old >> index_shift) + 1, 0, (limbs - (size_old >> index_shift) - 1) * sizeof(unsigned long))
    return 0

cdef StabilizerChain *SC_copy(StabilizerChain *SC, int level) noexcept:
    """
    Create a copy of the first `level` levels of SC. Must have 0 < level.

    Return a null pointer in case of allocation failure.
    """
    cdef int i, n = SC.degree
    cdef StabilizerChain *SCC = SC_new(n, False)
    if SCC is NULL:
        return NULL
    level = min(level, SC.base_size)
    for i in range(level):
        SCC.generators[i]   = <int *> sig_malloc( SC.array_size[i]*n * sizeof(int) )
        SCC.gen_inverses[i] = <int *> sig_malloc( SC.array_size[i]*n * sizeof(int) )
        if SCC.generators[i] is NULL or SCC.gen_inverses[i] is NULL:
            SC_dealloc(SCC)
            return NULL
        SCC.array_size[i] = SC.array_size[i]
    for i in range(level, n):
        SCC.generators[i]   = <int *> sig_malloc( default_num_gens*n * sizeof(int) )
        SCC.gen_inverses[i] = <int *> sig_malloc( default_num_gens*n * sizeof(int) )
        if SCC.generators[i] is NULL or SCC.gen_inverses[i] is NULL:
            SC_dealloc(SCC)
            return NULL
        SCC.array_size[i] = default_num_gens
    SC_copy_nomalloc(SCC, SC, level) # no chance for memory error here...
    return SCC

cdef int SC_copy_nomalloc(StabilizerChain *SC_dest, StabilizerChain *SC, int level) noexcept:
    cdef int i, n = SC.degree
    level = min(level, SC.base_size)
    SC_dest.base_size = level
    memcpy(SC_dest.orbit_sizes, SC.orbit_sizes, 2*n * sizeof(int) ) # copies orbit_sizes, num_gens
    memcpy(SC_dest.base_orbits[0], SC.base_orbits[0], 3*n*n * sizeof(int) ) # copies base_orbits, parents, labels
    for i in range(level):
        if SC.num_gens[i] > SC_dest.array_size[i]:
            if SC_realloc_gens(SC_dest, i, max(SC.num_gens[i], 2*SC_dest.array_size[i])):
                return 1
        memcpy(SC_dest.generators[i], SC.generators[i], SC.num_gens[i]*n * sizeof(int) )
        memcpy(SC_dest.gen_inverses[i], SC.gen_inverses[i], SC.num_gens[i]*n * sizeof(int) )
    return 0

cdef SC_print_level(StabilizerChain *SC, int level):
    cdef int i, j, n = SC.degree
    if level < SC.base_size:
        print('/ level {}'.format(level))
        print('| orbit   {}'.format([SC.base_orbits[level][i]
                                     for i in range(SC.orbit_sizes[level])]))
        print('| parents {}'.format([SC.parents[level][i] for i in range(n)]))
        print('| labels  {}'.format([SC.labels[level][i] for i in range(n)]))
        print('|')
        print('| generators  {}'.format([[SC.generators[level][n*i + j]
                                          for j in range(n)]
                                         for i in range(SC.num_gens[level])]))
        print(r'\ inverses    {}'.format([[SC.gen_inverses[level][n*i + j]
                                           for j in range(n)]
                                          for i in range(SC.num_gens[level])]))
    else:
        print('/ level {}'.format(level))
        print('|')
        print(r'\ base_size {}'.format(SC.base_size))


cdef StabilizerChain *SC_new_base(StabilizerChain *SC, int *base, int base_len) noexcept:
    """
    Create a new stabilizer chain whose base starts with the given base, and
    which represents the same permutation group. Original StabilizerChain is
    unmodified.

    Use SC_cleanup to remove redundant base points.

    Return a null pointer in case of an allocation failure.
    """
    cdef StabilizerChain *NEW = SC_new(SC.degree)
    if NEW is NULL:
        return NULL
    if SC_new_base_nomalloc(NEW, SC, base, base_len):
        SC_dealloc(NEW)
        return NULL
    return NEW

cdef int SC_new_base_nomalloc(StabilizerChain *SC_dest, StabilizerChain *SC, int *base, int base_len) noexcept:
    cdef int i
    SC_dest.base_size = 0
    for i in range(base_len):
        SC_add_base_point(SC_dest, base[i])
    if SC_update(SC_dest, SC, 0):
        SC_dealloc(SC_dest)
        return 1
    return 0

cdef int SC_update(StabilizerChain *dest, StabilizerChain *source, int level) noexcept:
    cdef mpz_t src_order, dst_order
    cdef int *perm = dest.perm_scratch
    mpz_init(src_order)
    mpz_init(dst_order)
    SC_order(source, level, src_order)
    SC_order(dest, level, dst_order)
    cdef int i, first_moved, b
    while mpz_cmp(dst_order, src_order):
        SC_random_element(source, level, perm)
        i = level
        while i < dest.base_size:
            b = dest.base_orbits[i][0]
            if perm[b] != b:
                break
            i += 1
        else:
            for b in range(dest.degree):
                if perm[b] != b:
                    break
            else:
                continue
            SC_add_base_point(dest, b)
        first_moved = i
        for i in range(level, first_moved + 1):
            if SC_insert_and_sift(dest, i, perm, 1, 0): # don't sift!
                mpz_clear(dst_order)
                mpz_clear(src_order)
                return 1
        SC_order(dest, level, dst_order)
    mpz_clear(src_order)
    mpz_clear(dst_order)
    return 0

cdef StabilizerChain *SC_insert_base_point(StabilizerChain *SC, int level, int p) noexcept:
    """
    Insert the point ``p`` as a base point on level ``level``. Return a new
    StabilizerChain with this new base. Original StabilizerChain is unmodified.

    Use SC_cleanup to remove redundant base points.

    Return a null pointer in case of an allocation failure.
    """
    cdef int i, b, n = SC.degree
    cdef StabilizerChain *NEW
    if level > 0:
        NEW = SC_copy(SC, level)
    else:
        NEW = SC_new(n)
    if NEW is NULL:
        return NULL
    SC_add_base_point(NEW, p)
    for i in range(level, SC.base_size):
        b = SC.base_orbits[i][0]
        if b != p:
            SC_add_base_point(NEW, b)
    if SC_update(NEW, SC, level):
        SC_dealloc(NEW)
        return NULL
    return NEW

cdef int SC_insert_base_point_nomalloc(StabilizerChain *SC_dest, StabilizerChain *SC, int level, int p) noexcept:
    cdef int i, b
    SC_copy_nomalloc(SC_dest, SC, level)
    SC_add_base_point(SC_dest, p)
    for i in range(level, SC.base_size):
        b = SC.base_orbits[i][0]
        if b != p:
            SC_add_base_point(SC_dest, b)
    if SC_update(SC_dest, SC, level):
        return 1
    return 0

cdef int SC_re_tree(StabilizerChain *SC, int level, int *perm, int x) noexcept:
    """
    Return values:
    0 - No errors.
    1 - Allocation failure.
    """
    cdef int *gen
    cdef int *gen_inv
    cdef int i, b, gen_index, n = SC.degree

    # make sure we have room for the new generator:
    if SC.array_size[level] == SC.num_gens[level]:
        if SC_realloc_gens(SC, level, 2*SC.array_size[level]):
            return 1
    cdef int *new_gen = SC.generators[level] + n*SC.num_gens[level]
    cdef int *new_gen_inv = SC.gen_inverses[level] + n*SC.num_gens[level]

    # new generator is perm^(-1) * (path from x to base) (left to right composition)
    SC_invert_perm(new_gen, perm, n)
    SC_compose_up_to_base(SC, level, x, new_gen)
    SC_invert_perm(new_gen_inv, new_gen, n)
    SC.num_gens[level] += 1

    # now that we have our generators, regenerate the tree, breadth-first
    b = SC.base_orbits[level][0]
    for i in range(n):
        SC.parents[level][i] = -1
    SC.parents[level][b] = b
    i = 0
    SC.orbit_sizes[level] = 1
    while i < SC.orbit_sizes[level]:
        x = SC.base_orbits[level][i]
        for gen_index from SC.num_gens[level] > gen_index >= 0:
            gen_inv = SC.gen_inverses[level] + n*gen_index
            SC_scan(SC, level, x, gen_index, gen_inv, 1)
        for gen_index from 0 <= gen_index < SC.num_gens[level]:
            gen = SC.generators[level] + n*gen_index
            SC_scan(SC, level, x, gen_index, gen, -1)
        i += 1
    return 0

cdef int SC_sift(StabilizerChain *SC, int level, int x, int *gens, int num_gens, int *new_gens) noexcept:
    """
    Apply Schreier's subgroup lemma[1] as follows. Given a level, a point x, and
    a generator s, find the coset traversal element r coming from x.
    Try inserting rsR(rs)^(-1) (composing left to right) on level+1, where
    R(g) is the coset representative of g.

    level - which level we are currently working on
    x - the object representing r
    gens - points to the generators to sift by
    num_gens - how many of these there are
    new_gens - space of size at least num_gens*n for the sifted perms to go

    Return 1 in case of an allocation failure.
    """
    cdef int n = SC.degree
    if num_gens == 0:
        return 0

    # copy a representative taking base to the point x to each of these
    cdef int i
    cdef int *temp = SC.gen_inverses[level] + n*SC.num_gens[level] # one more scratch space
                                                                   # (available since num_gens > 0)
    cdef int *rep_inv = temp
    SC_identify(rep_inv, n)
    SC_compose_up_to_base(SC, level, x, rep_inv)
    SC_invert_perm(new_gens, rep_inv, n)
    for i from 0 < i < num_gens:
        memcpy(new_gens + n*i, new_gens, n*sizeof(int))

    # post-compose each one with a generator
    for i from 0 <= i < num_gens:
        SC_mult_perms(new_gens + n*i, new_gens + n*i, gens + n*i, n)

    # for each one, look up the representative it now gives, and post-compose with that inverse
    cdef int y, b = SC.base_orbits[level][0]
    cdef int *perm
    cdef int *perm_rep_inv = temp
    for i from 0 <= i < num_gens:
        perm = new_gens + n*i # this is now rs
        y = perm[b]
        SC_identify(perm_rep_inv, n)
        SC_compose_up_to_base(SC, level, y, perm_rep_inv)
        SC_mult_perms(perm, perm, perm_rep_inv, n)
    return SC_insert(SC, level+1, new_gens, num_gens)

cdef int SC_insert_and_sift(StabilizerChain *SC, int level, int *pi, int num_perms, bint sift) noexcept:
    cdef int i, j, b, n = SC.degree
    cdef int perm_gen_index
    if sift:
        if SC_realloc_bitsets(SC, num_perms):
            return 1
        bitset_clear(&SC.gen_used)
        bitset_clear(&SC.gen_is_id)
        b = -1
        for perm_gen_index from 0 <= perm_gen_index < num_perms:
            for i from 0 <= i < n:
                if pi[n*perm_gen_index + i] != i:
                    b = i
                    break
            else:
                bitset_set(&SC.gen_is_id, perm_gen_index)
            if b != -1:
                break
        if b == -1:
            return 0
    if sift and level == SC.base_size:
        SC_add_base_point(SC, b)
    else:
        b = SC.base_orbits[level][0]
    # Now b is the base element at level `level`

    # Record the old orbit elements and the old generators (see sifting phase)
    cdef int old_num_gens = SC.num_gens[level]

    # Add new points to the tree:
    cdef int x
    cdef int *perm
    cdef int start_over = 1
    cdef int re_treed = 0
    while start_over:
        start_over = 0
        for i from 0 <= i < SC.orbit_sizes[level]:
            x = SC.base_orbits[level][i]
            for perm_gen_index from 0 <= perm_gen_index < num_perms:
                if sift and bitset_check(&SC.gen_is_id, perm_gen_index):
                    continue
                perm = pi + n*perm_gen_index
                if SC.parents[level][perm[x]] == -1:
                    # now we have an x which maps to a new point under perm,
                    re_treed = 1
                    if sift:
                        bitset_set(&SC.gen_used, perm_gen_index)
                    if SC_re_tree(SC, level, perm, x):
                        return 1
                    start_over = 1 # we must look anew
                    break
            if start_over:
                break
            if not re_treed:
                continue
            for perm_gen_index from 0 <= perm_gen_index < old_num_gens:
                perm = SC.generators[level] + n*perm_gen_index
                if SC.parents[level][perm[x]] == -1:
                    # now we have an x which maps to a new point under perm,
                    if SC_re_tree(SC, level, perm, x):
                        return 1
                    start_over = 1 # we must look anew
                    break
            if start_over:
                break
            for j from level < j < SC.base_size:
                for perm_gen_index from 0 <= perm_gen_index < SC.num_gens[j]:
                    perm = SC.generators[j] + n*perm_gen_index
                    if SC.parents[level][perm[x]] == -1:
                        # now we have an x which maps to a new point under perm,
                        if SC_re_tree(SC, level, perm, x):
                            return 1
                        start_over = 1 # we must look anew
                        break
    if not sift:
        return 0

    # store the rest of the unused generators in the generator array, after the added ones
    cdef int new_size
    cdef int unused_gens = 0
    for perm_gen_index from 0 <= perm_gen_index < num_perms:
        if not bitset_check(&SC.gen_used, perm_gen_index) and not bitset_check(&SC.gen_is_id, perm_gen_index):
            unused_gens += 1
    if 2*(SC.num_gens[level] + unused_gens) > SC.array_size[level]:
        new_size = max(2*(SC.num_gens[level] + unused_gens), 2*SC.array_size[level])
        if SC_realloc_gens(SC, level, new_size):
            return 1
    i = 0
    for perm_gen_index from 0 <= perm_gen_index < num_perms:
        if not bitset_check(&SC.gen_used, perm_gen_index) and not bitset_check(&SC.gen_is_id, perm_gen_index):
            memcpy(SC.generators[level] + n*(SC.num_gens[level]+i), pi + n*perm_gen_index, n*sizeof(int))
            i += 1

    # Sift:
    cdef int *gens = SC.generators[level]
    cdef int total_gens = SC.num_gens[level] + unused_gens
    cdef int section, gens_in_section
    for i from 0 <= i < SC.orbit_sizes[level]:
        x = SC.base_orbits[level][i]
        section = 0
        while section*n < total_gens:
            gens_in_section = min(n, total_gens - section*n)
            if SC_sift(SC, level, x, gens + n*n*section, gens_in_section, gens + n*total_gens):
                return 1
            section += 1
    return 0

cdef bint SC_is_giant(int n, int num_perms, int *perms, float p, bitset_t support) noexcept:
    """
    Test whether the group generated by the input permutations is a giant, i.e.,
    the alternating or symmetric group.

    If the group is not a giant, this routine will return False. This could also
    indicate an allocation failure.

    If the group is a giant, this routine will return ``True`` with approximate
    probability ``p``. It will set `support' to the support of the group in this
    case. Use bitset_len to get the size of support.

    The bitset `support' must be initialized. Must have 0 <= p < 1.

    Running time is roughly O(-ln(1-p)*n*ln(m)) where m <= n is the size of the
    support of the group.
    """
    cdef int i, j, num_steps, support_root
    cdef size_t m = 1
    cdef unsigned long q
    cdef int *gen
    cdef int *perm = <int *> sig_malloc(n*sizeof(int))
    cdef OrbitPartition *OP = OP_new(n)
    if OP is NULL or perm is NULL:
        OP_dealloc(OP)
        sig_free(perm)
        return False

    # giants are transitive
    for j from 0 <= j < num_perms:
        gen = perms + n*j
        for i from 0 <= i < n:
            OP_join(OP, i, gen[i])
    for i from 0 <= i < n:
        if OP.parent[i] == i:
            if OP.size[i] != 1:
                if m != 1:
                    m = 1
                    break
                else:
                    m = OP.size[i]
                    support_root = i
    if m == 1:
        OP_dealloc(OP)
        sig_free(perm)
        return False
    bitset_zero(support)
    for i from 0 <= i < n:
        if OP_find(OP, i) == support_root:
            bitset_set(support, i)

    # get a bit lost in the group, so our random elements are more random:
    SC_identify(perm, n)
    for i from 0 <= i < 10:
        SC_mult_perms(perm, perm, perms + n*(rand()%num_perms), n)

    # look for elements with cycles of prime length q, m/2 < q < m-2
    num_steps = <int> ceil(-log(1-p)*log(m)/log(2))
    for j from 0 <= j < num_steps:
        OP_clear(OP)
        for i from 0 <= i < n:
            OP_join(OP, i, perm[i])
        for i from 0 <= i < n:
            if OP.parent[i] == i:
                q = OP.size[i]
                if m < 2*q and q < m-2:
                    if n_is_prime(q):
                        sig_free(perm)
                        OP_dealloc(OP)
                        return True
        SC_mult_perms(perm, perm, perms + n*(rand()%num_perms), n)
    OP_dealloc(OP)
    sig_free(perm)
    return False


def SC_test_list_perms(list L, int n, int limit, bint gap, bint limit_complain, bint test_contains):
    """
    Test that the permutation group generated by list perms in L of degree n
    is of the correct order, by comparing with GAP. Don't test if the group is
    of size greater than limit.

    TESTS::

        sage: # needs sage.groups
        sage: from sage.groups.perm_gps.partn_ref.data_structures import SC_test_list_perms
        sage: limit = 10^7
        sage: def test_Sn_on_m_points(n, m, gap, contains):
        ....:     perm1 = [1,0] + list(range(2, m))
        ....:     perm2 = [(i+1)%n for i in range( n )] + list(range(n,m))
        ....:     SC_test_list_perms([perm1, perm2], m, limit, gap, 0, contains)
        sage: for i in range(2,9):
        ....:     test_Sn_on_m_points(i,i,1,0)
        sage: for i in range(2,9):
        ....:     test_Sn_on_m_points(i,i,0,1)
        sage: for i in range(2,9):              # long time
        ....:     test_Sn_on_m_points(i,i,1,1)
        sage: test_Sn_on_m_points(8,8,1,1)
        sage: def test_stab_chain_fns_1(n, gap, contains):
        ....:     perm1 = sum([[2*i+1,2*i] for i in range(n)], [])
        ....:     perm2 = [(i+1)%(2*n) for i in range( 2*n)]
        ....:     SC_test_list_perms([perm1, perm2], 2*n, limit, gap, 0, contains)
        sage: for n in range(1,11):
        ....:     test_stab_chain_fns_1(n, 1, 0)
        sage: for n in range(1,11):
        ....:     test_stab_chain_fns_1(n, 0, 1)
        sage: for n in range(1,9):              # long time
        ....:     test_stab_chain_fns_1(n, 1, 1)
        sage: test_stab_chain_fns_1(11, 1, 1)
        sage: def test_stab_chain_fns_2(n, gap, contains):
        ....:     perms = []
        ....:     for p,e in factor(n):
        ....:         perm1 = [(p*(i//p)) + ((i+1)%p) for i in range(n)]
        ....:         perms.append(perm1)
        ....:     SC_test_list_perms(perms, n, limit, gap, 0, contains)
        sage: for n in range(2,11):
        ....:     test_stab_chain_fns_2(n, 1, 0)
        sage: for n in range(2,11):
        ....:     test_stab_chain_fns_2(n, 0, 1)
        sage: for n in range(2,11):             # long time
        ....:     test_stab_chain_fns_2(n, 1, 1)
        sage: test_stab_chain_fns_2(11, 1, 1)
        sage: def test_stab_chain_fns_3(n, gap, contains):
        ....:     perm1 = [(-i)%n for i in range( n )]
        ....:     perm2 = [(i+1)%n for i in range( n )]
        ....:     SC_test_list_perms([perm1, perm2], n, limit, gap, 0, contains)
        sage: for n in range(2,20):
        ....:     test_stab_chain_fns_3(n, 1, 0)
        sage: for n in range(2,20):
        ....:     test_stab_chain_fns_3(n, 0, 1)
        sage: for n in range(2,14):             # long time
        ....:     test_stab_chain_fns_3(n, 1, 1)
        sage: test_stab_chain_fns_3(20, 1, 1)
        sage: def test_stab_chain_fns_4(n, g, gap, contains):
        ....:     perms = []
        ....:     for _ in range(g):
        ....:         perm = list(range(n))
        ....:         shuffle(perm)
        ....:         perms.append(perm)
        ....:     SC_test_list_perms(perms, n, limit, gap, 0, contains)
        sage: for n in range(4,9):              # long time
        ....:     test_stab_chain_fns_4(n, 1, 1, 0)
        ....:     test_stab_chain_fns_4(n, 2, 1, 0)
        ....:     test_stab_chain_fns_4(n, 2, 1, 0)
        ....:     test_stab_chain_fns_4(n, 2, 1, 0)
        ....:     test_stab_chain_fns_4(n, 2, 1, 0)
        ....:     test_stab_chain_fns_4(n, 3, 1, 0)
        sage: for n in range(4,9):              # known bug (see :issue:`32187`), not tested
        ....:     test_stab_chain_fns_4(n, 1, 0, 1)
        ....:     for j in range(6):
        ....:         test_stab_chain_fns_4(n, 2, 0, 1)
        ....:     test_stab_chain_fns_4(n, 3, 0, 1)
        sage: for n in range(4,8):              # long time
        ....:     test_stab_chain_fns_4(n, 1, 1, 1)
        ....:     test_stab_chain_fns_4(n, 2, 1, 1)
        ....:     test_stab_chain_fns_4(n, 2, 1, 1)
        ....:     test_stab_chain_fns_4(n, 3, 1, 1)
        sage: test_stab_chain_fns_4(8, 2, 1, 1)
        sage: def test_stab_chain_fns_5(n, gap, contains):
        ....:     perms = []
        ....:     m = n//3
        ....:     perm1 = list(range(2*m))
        ....:     shuffle(perm1)
        ....:     perm1 += list(range(2*m,n))
        ....:     perm2 = list(range(m,n))
        ....:     shuffle(perm2)
        ....:     perm2 = list(range(m)) + perm2
        ....:     SC_test_list_perms([perm1, perm2], n, limit, gap, 0, contains)
        sage: for n in [4..9]:                  # long time
        ....:     for _ in range(2):
        ....:         test_stab_chain_fns_5(n, 1, 0)
        sage: for n in [4..8]:                  # long time
        ....:     test_stab_chain_fns_5(n, 0, 1)
        sage: for n in [4..9]:                  # long time
        ....:     test_stab_chain_fns_5(n, 1, 1)
        sage: def random_perm(x):
        ....:     shuffle(x)
        ....:     return x
        sage: def test_stab_chain_fns_6(m, n, k, gap, contains):
        ....:     perms = []
        ....:     for i in range(k):
        ....:         perm = sum([random_perm(list(range(i*(n//m),min(n,(i+1)*(n//m))))) for i in range(m)], [])
        ....:         perms.append(perm)
        ....:     SC_test_list_perms(perms, m*(n//m), limit, gap, 0, contains)
        sage: for m in range(2,9):              # long time
        ....:     for n in range(m,3*m):
        ....:         for k in range(1,3):
        ....:             test_stab_chain_fns_6(m,n,k, 1, 0)
        sage: for m in range(2,10):
        ....:     for n in range(m,4*m):
        ....:         for k in range(1,3):
        ....:             test_stab_chain_fns_6(m,n,k, 0, 1)
        sage: test_stab_chain_fns_6(10,20,2, 1, 1)
        sage: test_stab_chain_fns_6(8,16,2, 1, 1)
        sage: test_stab_chain_fns_6(6,36,2, 1, 1)
        sage: test_stab_chain_fns_6(4,40,3, 1, 1)
        sage: test_stab_chain_fns_6(4,40,2, 1, 1)
        sage: def test_stab_chain_fns_7(n, cop, gap, contains):
        ....:     perms = []
        ....:     for i in range(0,n//2,2):
        ....:         p = list(range(n))
        ....:         p[i] = i+1
        ....:         p[i+1] = i
        ....:     if cop:
        ....:         perms.append([c for c in p])
        ....:     else:
        ....:         perms.append(p)
        ....:     SC_test_list_perms(perms, n, limit, gap, 0, contains)
        sage: for n in [6..14]:
        ....:     test_stab_chain_fns_7(n, 1, 1, 0)
        ....:     test_stab_chain_fns_7(n, 0, 1, 0)
        sage: for n in [6..30]:
        ....:     test_stab_chain_fns_7(n, 1, 0, 1)
        ....:     test_stab_chain_fns_7(n, 0, 0, 1)
        sage: for n in [6..14]:                 # long time
        ....:     test_stab_chain_fns_7(n, 1, 1, 1)
        ....:     test_stab_chain_fns_7(n, 0, 1, 1)
        sage: test_stab_chain_fns_7(20, 1, 1, 1)
        sage: test_stab_chain_fns_7(20, 0, 1, 1)
    """
    if gap:
        from sage.groups.perm_gps.permgroup import PermutationGroup
        from sage.groups.perm_gps.constructor import PermutationGroupElement
        from sage.misc.prandom import shuffle

    cdef StabilizerChain *SC
    cdef StabilizerChain *SCC
    cdef StabilizerChain *SCCC
    cdef StabilizerChain *SC_nb
    cdef Integer order, order2
    cdef int i, j, m, SC_says
    cdef bitset_t giant_support
    if gap:
        G = PermutationGroup([[i+1 for i in p] for p in L])
        if G.order() > limit:
            if limit_complain:
                print('TOO BIG')
            return
    SC = SC_new(n)
    cdef int *perm = <int *>sig_malloc(n * (len(L)+3) * sizeof(int))
    try:
        bitset_init(giant_support, n)
    except MemoryError:
        sig_free(perm)
        SC_dealloc(SC)
        raise MemoryError
    if perm is NULL or SC is NULL:
        bitset_free(giant_support)
        sig_free(perm)
        SC_dealloc(SC)
        raise MemoryError
    cdef int *perm2 = perm +   n
    cdef int *perm3 = perm + 2*n
    for Lperm in L:
        for i from 0 <= i < n:
            perm[i] = Lperm[i]
        if SC_insert(SC, 0, perm, 1):
            bitset_free(giant_support)
            sig_free(perm)
            SC_dealloc(SC)
            raise MemoryError
    SCC = SC_copy(SC, n)
    SCCC = SC_insert_base_point(SC, 0, n-1)
    for i from 0 <= i < n:
        perm[i] = n-i-1
    SC_nb = SC_new_base(SC, perm, n)
    if SCC is NULL or SCCC is NULL or SC_nb is NULL:
        bitset_free(giant_support)
        sig_free(perm)
        SC_dealloc(SC)
        SC_dealloc(SCC)
        SC_dealloc(SCCC)
        SC_dealloc(SC_nb)
        raise MemoryError
    giant = False
    try:
        order = Integer(0)
        SC_order(SC,0,order.value)
        j = 0
        for Lperm in L:
            for i from 0 <= i < n:
                perm[n*j + i] = Lperm[i]
            j += 1
        if SC_is_giant(n, len(L), perm, 0.9, giant_support):
            giant = True
            m = bitset_len(giant_support)
            from sage.arith.misc import factorial
            if not (order == factorial(m) or order == factorial(m)/2):
                print("SC_is_giant failed: %s %s"%(str(L), order))
                raise AssertionError
            if order == factorial(n):
                SC_dealloc(SC)
                SC = SC_symmetric_group(n)
                SC_order(SC,0,order.value)
                if not order == factorial(n):
                    print("SC_symmetric_group failed: %s %s"%(str(L), order))
                    raise AssertionError
            elif order == factorial(n)/2:
                SC_dealloc(SC)
                SC = SC_alternating_group(n)
                SC_order(SC,0,order.value)
                if not order == factorial(n)/2:
                    print("SC_alternating_group failed: %s %s"%(str(L), order))
                    raise AssertionError
        order2 = Integer(0)
        SC_order(SCC,0,order2.value)
        if order != order2:
            print("FAIL {}".format(L))
            print('SC_copy(n) does not agree with order of original {} {}'.format(order, order2))
            raise AssertionError
        SC_order(SCCC,0,order2.value)
        if order != order2:
            print("FAIL {}".format(L))
            print('does not agree with order of inserted base point {} {}'.format(order, order2))
            raise AssertionError
        SC_order(SC_nb,0,order2.value)
        if order != order2:
            print("FAIL {}".format(L))
            print('does not agree with order of new base {} {}'.format(order, order2))
            raise AssertionError
        if test_contains:
            for i from 0 <= i < 3:
                SC_random_element(SC, 0, perm)
                if not SC_contains(SC, 0, perm, 0):
                    print("FAIL {}".format(L))
                    print('element {}'.format([perm[ii] for ii in range(n)]))
                    print('SC_random_element says it is an element, SC_contains(modify=0) does not')
                    raise AssertionError
                if not SC_contains(SC, 0, perm, 1):
                    print("FAIL {}".format(L))
                    print('element {}'.format([perm[ii] for ii in range(n)]))
                    print('SC_random_element says it is an element, SC_contains(modify=1) does not')
                    raise AssertionError
                if not SC_contains(SCC, 0, perm, 0):
                    print("FAIL {}".format(L))
                    print('element {}'.format([perm[ii] for ii in range(n)]))
                    print('SC_random_element says it is an element, SC_contains(modify=0) does not on copy')
                    raise AssertionError
                if not SC_contains(SCC, 0, perm, 1):
                    print("FAIL {}".format(L))
                    print('element {}'.format([perm[ii] for ii in range(n)]))
                    print('SC_random_element says it is an element, SC_contains(modify=1) does not on copy')
                    raise AssertionError
                if not SC_contains(SCCC, 0, perm, 0):
                    print("FAIL {}".format(L))
                    print('element {}'.format([perm[ii] for ii in range(n)]))
                    print('SC_random_element says it is an element, SC_contains(modify=0) does not on inserted base point')
                    raise AssertionError
                if not SC_contains(SCCC, 0, perm, 1):
                    print("FAIL {}".format(L))
                    print('element {}'.format([perm[ii] for ii in range(n)]))
                    print('SC_random_element says it is an element, SC_contains(modify=1) does not on inserted base point')
                    raise AssertionError
                if not SC_contains(SC_nb, 0, perm, 0):
                    print("FAIL {}".format(L))
                    print('element {}'.format([perm[ii] for ii in range(n)]))
                    print('SC_random_element says it is an element, SC_contains(modify=0) does not on new base')
                    raise AssertionError
                if not SC_contains(SC_nb, 0, perm, 1):
                    print("FAIL {}".format(L))
                    print('element {}'.format([perm[ii] for ii in range(n)]))
                    print('SC_random_element says it is an element, SC_contains(modify=1) does not on new base')
                    raise AssertionError
                SC_random_element(SCC, 0, perm2)
                if not SC_contains(SC, 0, perm2, 0):
                    print("FAIL {}".format(L))
                    print('element {}'.format([perm[ii] for ii in range(n)]))
                    print('SC_random_element says it is an element of copy, SC_contains(modify=0) does not')
                    raise AssertionError
                if not SC_contains(SC, 0, perm2, 1):
                    print("FAIL {}".format(L))
                    print('element {}'.format([perm[ii] for ii in range(n)]))
                    print('SC_random_element says it is an element of copy, SC_contains(modify=1) does not')
                    raise AssertionError
                SC_random_element(SCCC, 0, perm3)
                if not SC_contains(SC, 0, perm3, 0):
                    print("FAIL {}".format(L))
                    print('element {}'.format([perm[ii] for ii in range(n)]))
                    print('SC_random_element says it is an element of inserted base point, SC_contains(modify=0) does not')
                    raise AssertionError
                if not SC_contains(SC, 0, perm3, 1):
                    print("FAIL {}".format(L))
                    print('element {}'.format([perm[ii] for ii in range(n)]))
                    print('SC_random_element says it is an element of inserted base point, SC_contains(modify=1) does not')
                    raise AssertionError
                SC_random_element(SC_nb, 0, perm3)
                if not SC_contains(SC, 0, perm3, 0):
                    print("FAIL {}".format(L))
                    print('element {}'.format([perm[ii] for ii in range(n)]))
                    print('SC_random_element says it is an element of new base, SC_contains(modify=0) does not')
                    raise AssertionError
                if not SC_contains(SC, 0, perm3, 1):
                    print("FAIL {}".format(L))
                    print('element {}'.format([perm[ii] for ii in range(n)]))
                    print('SC_random_element says it is an element of new base, SC_contains(modify=1) does not')
                    raise AssertionError
        if gap:
            order = Integer(0)
            SC_order(SC,0,order.value)
            j = 0
            for Lperm in L:
                for i from 0 <= i < n:
                    perm[n*j + i] = Lperm[i]
                j += 1
            if SC_is_giant(n, len(L), perm, 0.9, giant_support):
                from sage.arith.misc import factorial
                m = bitset_len(giant_support)
                if order != factorial(m) and order != factorial(m)/2:
                    print("SC_is_giant failed: %s %s"%(str(L), order))
                    raise AssertionError
            if order != G.order():
                print("FAIL {}".format(L))
                print('order was computed to be {}'.format(order))
                print('GAP says it is {}'.format(G.order()))
                raise AssertionError
            if test_contains:
                for i from 0 <= i < 3:
                    permy = G.random_element()
                    for j from 0 <= j < n:
                        perm[j] = permy(j+1)-1
                    if not SC_contains(SC, 0, perm, 0):
                        print("FAIL {}".format(L))
                        print('element {}'.format(permy))
                        print('GAP says it is an element, SC_contains(modify=0) does not')
                        raise AssertionError
                    if not SC_contains(SC, 0, perm, 1):
                        print("FAIL {}".format(L))
                        print('element {}'.format(permy))
                        print('GAP says it is an element, SC_contains(modify=1) does not')
                        raise AssertionError
                    permy = list(range(1, n + 1))
                    shuffle(permy)
                    gap_says = (PermutationGroupElement(permy) in G)
                    for j from 0 <= j < n:
                        perm[j] = permy[j]-1
                    SC_says = SC_contains(SC, 0, perm, 0)
                    if bool(SC_says) != bool(gap_says):
                        print("FAIL {}".format(L))
                        print('element {}'.format(permy))
                        print('GAP says %d, SC_contains(modify=0) says %d'%(gap_says, SC_says))
                        raise AssertionError
                    SC_says = SC_contains(SC, 0, perm, 1)
                    if bool(SC_says) != bool(gap_says):
                        print("FAIL {}".format(L))
                        print('element {}'.format(permy))
                        print('GAP says %d, SC_contains(modify=0) says %d'%(gap_says, SC_says))
                        raise AssertionError
                    SC_random_element(SC, 0, perm)
                    for j from 0 <= j < n:
                        permy[j] = perm[j]+1
                    gap_says = (PermutationGroupElement(permy) in G)
                    if not SC_contains(SC, 0, perm, 0):
                        print("FAIL {}".format(L))
                        print('element {}'.format(permy))
                        print('random element not contained in group, modify=false')
                        raise AssertionError
                    if not SC_contains(SC, 0, perm, 1):
                        print("FAIL {}".format(L))
                        print('element {}'.format(permy))
                        print('random element not contained in group, modify=true')
                        raise AssertionError
                    if not gap_says:
                        print("FAIL {}".format(L))
                        print('element {}'.format(permy))
                        print('random element not contained in group, according to GAP')
                        raise AssertionError
    except Exception:
        if giant:
            print("detected giant!")
        raise
    finally:
        bitset_free(giant_support)
        sig_free(perm)
        SC_dealloc(SC)
        SC_dealloc(SCC)
        SC_dealloc(SCCC)
        SC_dealloc(SC_nb)


# Functions

cdef int sort_by_function(PartitionStack *PS, int start, int *degrees) noexcept:
    """
    A simple counting sort, given the degrees of vertices to a certain cell.

    INPUT:

    - ``PS`` -- the partition stack to be checked
    - ``start`` -- beginning index of the cell to be sorted
    - ``degrees`` -- the values to be sorted by, must have extra scratch space for a
      total of `3*n+1`
    """
    cdef int n = PS.degree
    cdef int i, j, max, max_location
    cdef int *counts = degrees + n
    cdef int *output = degrees + 2*n + 1
    for i from 0 <= i <= n:
        counts[i] = 0
    i = 0
    while PS.levels[i+start] > PS.depth:
        counts[degrees[i]] += 1
        i += 1
    counts[degrees[i]] += 1
    # i+start is the right endpoint of the cell now
    max = counts[0]
    max_location = 0
    for j from 0 < j <= n:
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
    while j <= n and counts[j] <= i:
        if counts[j] > 0:
            PS.levels[start + counts[j] - 1] = PS.depth
        PS_move_min_to_front(PS, start + counts[j-1], start + counts[j] - 1)
        j += 1
    return max_location

cdef int refine_by_orbits(PartitionStack *PS, StabilizerChain *SC, int *perm_stack, int *cells_to_refine_by, int *ctrb_len) noexcept:
    """
    Given a stabilizer chain SC, refine the partition stack PS so that each cell
    contains elements from at most one orbit, and sort the refined cells by
    their minimal cell representatives in the orbit of the group.
    """
    cdef int start, level, gen_index, i, j, k, x, n = SC.degree
    cdef int *gen
    cdef int *min_cell_reps = SC.perm_scratch
    cdef int *perm = perm_stack + n*PS.depth
    cdef OrbitPartition *OP = SC.OP_scratch
    cdef int invariant = 1
    OP_clear(OP)
    for level from PS.depth <= level < SC.base_size:
        for gen_index from 0 <= gen_index < SC.num_gens[level]:
            gen = SC.generators[level] + gen_index*n
            for i from 0 <= i < n:
                OP_join(OP, i, gen[i])
    start = 0
    while start < n:
        i = 0
        while 1:
            x = PS.entries[start+i]
            x = perm[x]
            min_cell_reps[i] = OP.mcr[OP_find(OP, x)]
            i += 1
            if PS.levels[start+i-1] <= PS.depth:
                break
        invariant += sort_by_function(PS, start, min_cell_reps)
        invariant += i
        # update cells_to_refine_by
        k = start
        for j from start <= j < start+i:
            if PS.levels[j] == PS.depth:
                k = j+1
                cells_to_refine_by[ctrb_len[0]] = k
                ctrb_len[0] += 1
        start += i
    return invariant

cdef int compute_relabeling(StabilizerChain *group, StabilizerChain *scratch_group,
    int *permutation, int *relabeling) noexcept:
    """
    Technically, compute the INVERSE of the relabeling
    """
    cdef int i, j, x, y, m, n = group.degree, orbit_element
    cdef int *scratch = group.perm_scratch
    if SC_new_base_nomalloc(scratch_group, group, permutation, n):
        return 1
    group = scratch_group
    SC_identify(relabeling, n)
    for i from 0 <= i < n:
        m = n
        for j from 0 <= j < group.orbit_sizes[i]:
            orbit_element = group.base_orbits[i][j]
            x = relabeling[orbit_element]
            if x < m:
                m = x
                y = orbit_element
        SC_invert_perm(scratch, relabeling, n)
        SC_compose_up_to_base(group, i, y, scratch)
        SC_invert_perm(relabeling, scratch, n)
    SC_invert_perm(scratch, relabeling, n)
    memcpy(relabeling, scratch, n*sizeof(int))
    return 0
