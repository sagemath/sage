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

cdef struct aut_gp_and_can_lab:
    int *generators
    int size_of_generator_array
    int num_gens
    StabilizerChain *group
    int *relabeling

cdef aut_gp_and_can_lab *allocate_agcl_output(int)

cdef void deallocate_agcl_output(aut_gp_and_can_lab *)

cdef struct agcl_work_space:
    int degree
    # for nontrivial input groups
    int *perm_stack # n^2 ints
    StabilizerChain *group1 # degree n
    StabilizerChain *group2 # degree n
    # for canonical labels
    int *label_indicators # n ints
    PartitionStack *label_ps # degree n
    # the rest
    int *int_array # 7*n ints
    bitset_t *bitset_array # (n + 2*len_of_fp_and_mcr + 1) bitset_t's, each of size n
    OrbitPartition *orbits_of_subgroup # degree n
    OrbitPartition *orbits_of_permutation # degree n
    PartitionStack *first_ps # degree n

cdef agcl_work_space *allocate_agcl_work_space(int)

cdef void deallocate_agcl_work_space(agcl_work_space *)

cdef aut_gp_and_can_lab *get_aut_gp_and_can_lab( void *,
    PartitionStack *, int,
    bint (*)(PartitionStack *, void *),
    int (*)(PartitionStack *, void *, int *, int),
    int (*)(int *, int *, void *, void *, int), bint, StabilizerChain *,
    agcl_work_space *, aut_gp_and_can_lab *) except NULL



