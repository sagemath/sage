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
from sage.graphs.base.c_graph cimport CGraph
from sage.groups.perm_gps.partn_ref.automorphism_group_canonical_label cimport (
    get_aut_gp_and_can_lab, aut_gp_and_can_lab, agcl_work_space,
    allocate_agcl_output, deallocate_agcl_output,
    allocate_agcl_work_space, deallocate_agcl_work_space)
from sage.groups.perm_gps.partn_ref.canonical_augmentation cimport (iterator,
    canonical_generator_data, allocate_cgd, deallocate_cgd,
    canonical_generator_next,
    setup_canonical_generator, start_canonical_generator)
from sage.groups.perm_gps.partn_ref.refinement_sets cimport (subset, free_subset, all_set_children_are_equivalent,
    refine_set, compare_sets, generate_child_subsets, apply_subset_aug,
    canonical_set_parent, allocate_sgd, deallocate_sgd, allocate_subset_gen, free_subset_gen,
    setup_set_gen, subset_generator_next, subset_generator_data, allocate_subset_gen_2)


cdef class GraphStruct:
    cdef CGraph G
    cdef bint directed
    cdef bint loops
    cdef bint use_indicator
    cdef int *scratch # length 3n+1


cdef struct dg_edge_gen_data:
    iterator *edge_iterator
    void *graph
