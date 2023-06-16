r"""
Declaration file for canonical augmentation

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
from .double_coset cimport (double_coset,
    dc_work_space, allocate_dc_work_space, deallocate_dc_work_space)


cdef struct iterator:
    void *data
    void *(*next)(void *data, int *degree, bint *mem_err)

cdef struct canonical_generator_data:
    StabilizerChain *group

    void **object_stack
    int *degree_stack
    iterator *iterator_stack
    aut_gp_and_can_lab **aut_gp_stack
    agcl_work_space **agcl_work_spaces
    dc_work_space **dc_work_spaces
    PartitionStack **ps_stack
    void **aug_stack
    void **parent_stack

    int level
    int max_level
    int allocd_levels
    bint reduce_children
    bint mem_err
    bint dealloc
    bint pr

    bint (*all_children_are_equivalent)(PartitionStack *, void *)
    int (*refine_and_return_invariant)(PartitionStack *, void *, int *, int)
    int (*compare_structures)(int *, int *, void *, void *, int)

    int (*generate_children)(void *, aut_gp_and_can_lab *, iterator *)
    void *(*apply_augmentation)(void *, void *, void *, int *, bint *)
    void (*free_object)(void *)
    void (* free_iter_data)(void *)
    void (*free_aug)(void *)
    void *(*canonical_parent)(void *child, void *parent, int *permutation, int *degree, bint *)

cdef canonical_generator_data *allocate_cgd(int, int)

cdef void deallocate_cgd(canonical_generator_data *)

cdef void *canonical_generator_next(void *, int *, bint *)

cdef iterator *setup_canonical_generator(int degree,
    bint (*all_children_are_equivalent)(PartitionStack *, void *),
    int (*refine_and_return_invariant)(PartitionStack *, void *, int *, int),
    int (*compare_structures)(int *, int *, void *, void *, int),
    int (*generate_children)(void *, aut_gp_and_can_lab *, iterator *),
    void *(*apply_augmentation)(void *, void *, void *, int *, bint *),
    void (*free_object)(void *),
    void (* free_iter_data)(void *),
    void (*free_aug)(void *),
    void *(*canonical_parent)(void *, void *, int *, int *, bint *),
    int max_depth, bint reduce_children,
    iterator *cangen_prealloc) except NULL

cdef iterator *start_canonical_generator(StabilizerChain *, void *, int, iterator *) except NULL

