#*****************************************************************************
#       Copyright (C) 2006 - 2011 Robert L. Miller <rlmillster@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .data_structures cimport *
from sage.data_structures.bitset cimport bitset_t
from sage.rings.integer cimport Integer

cdef inline int int_cmp(int a, int b):
    if a < b:
        return -1
    elif a == b:
        return 0
    else:
        return 1

cdef struct dc_work_space:
    int degree
    # for nontrivial input groups
    int *perm_stack # n^2 ints
    StabilizerChain *group1 # degree n
    StabilizerChain *group2 # degree n
    # the rest
    PartitionStack *current_ps # degree n
    PartitionStack *first_ps # degree n
    int *int_array # 5*n ints
    bitset_t *bitset_array # (n + 2*len_of_fp_and_mcr + 1) bitset_t's, each of size n
    OrbitPartition *orbits_of_subgroup # degree n

cdef dc_work_space *allocate_dc_work_space(int)

cdef void deallocate_dc_work_space(dc_work_space *)

cdef int double_coset( void *, void *, PartitionStack *, int *, int,
    bint (*)(PartitionStack *, void *),
    int (*)(PartitionStack *, void *, int *, int),
    int (*)(int *, int *, void *, void *, int),
    StabilizerChain *, dc_work_space *, int *) except -1
