r"""
Declaration file for simple set datastructures

AUTHORS:

- Robert Miller (2011--2013): initial version

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

from .data_structures cimport *
from .automorphism_group_canonical_label cimport (
    get_aut_gp_and_can_lab, aut_gp_and_can_lab, agcl_work_space,
    allocate_agcl_output, deallocate_agcl_output,
    allocate_agcl_work_space, deallocate_agcl_work_space)
from .canonical_augmentation cimport (iterator,
    canonical_generator_data, allocate_cgd, deallocate_cgd,
    canonical_generator_next,
    setup_canonical_generator, start_canonical_generator)


cdef struct subset:
    bitset_s bits
    int *scratch # must be of size 3*n + 1

cdef int refine_set(PartitionStack *, void *, int *, int)
cdef int compare_sets(int *, int *, void *, void *, int)
cdef bint all_set_children_are_equivalent(PartitionStack *, void *)

cdef void *allocate_subset(int)

cdef struct subset_generator_data:
    OrbitPartition *orbits
    int cur_point
    bitset_s bits

cdef void *allocate_sgd(int)
cdef void deallocate_sgd(void *)

cdef void *subset_generator_next(void *, int *, bint *)

cdef int generate_child_subsets(void *S, aut_gp_and_can_lab *group, iterator *it)
cdef void *apply_subset_aug(void *parent, void *aug, void *child, int *degree, bint *mem_err)
cdef void free_subset(void *child)
cdef void free_subset_aug(void *)
cdef void *canonical_set_parent(void *child, void *parent, int *permutation, int *degree, bint *mem_err)

cdef iterator *allocate_subset_gen(int degree, int max_size)
cdef int allocate_subset_gen_2(int degree, int max_size, iterator *it)
cdef void free_subset_gen(iterator *subset_gen)
cdef iterator *setup_set_gen(iterator *subset_gen, int degree, int max_size)
