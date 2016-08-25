#*****************************************************************************
#      Copyright (C) 2006 - 2011 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

include 'data_structures_pxd.pxi' # includes bitsets

from sage.graphs.base.c_graph cimport CGraph
from .automorphism_group_canonical_label cimport (
    get_aut_gp_and_can_lab, aut_gp_and_can_lab, agcl_work_space,
    allocate_agcl_output, deallocate_agcl_output,
    allocate_agcl_work_space, deallocate_agcl_work_space)
from .canonical_augmentation cimport (iterator,
    canonical_generator_data, allocate_cgd, deallocate_cgd,
    canonical_generator_next,
    setup_canonical_generator, start_canonical_generator)
from .refinement_sets cimport (subset, free_subset, all_set_children_are_equivalent,
    refine_set, compare_sets, generate_child_subsets, apply_subset_aug,
    canonical_set_parent, allocate_sgd, deallocate_sgd, allocate_subset_gen, free_subset_gen,
    setup_set_gen, subset_generator_next, subset_generator_data, allocate_subset_gen_2)


cdef class GraphStruct:
    cdef CGraph G
    cdef bint directed
    cdef bint loops
    cdef bint use_indicator
    cdef int *scratch # length 3n+1

cdef int refine_by_degree(PartitionStack *, void *, int *, int)
cdef int compare_graphs(int *, int *, void *, void *, int)
cdef bint all_children_are_equivalent(PartitionStack *, void *)
cdef inline int degree(PartitionStack *, CGraph, int, int, bint)

cdef struct dg_edge_gen_data:
    iterator *edge_iterator
    void *graph

cdef void *dg_edge_gen_next(void *, int *, bint *)

cdef int gen_children_dg_edge(void *S, aut_gp_and_can_lab *group, iterator *it)
cdef void *apply_dg_edge_aug(void *parent, void *aug, void *child, int *degree, bint *mem_err)
cdef void free_dg_edge(void *child)
cdef void *canonical_dg_edge_parent(void *child, void *parent, int *permutation, int *degree, bint *mem_err)

cdef int gen_children_dg_vert(void *S, aut_gp_and_can_lab *group, iterator *it)
cdef void *apply_dg_vert_aug(void *parent, void *aug, void *child, int *degree, bint *mem_err)
cdef void free_dg_vert(void *child)
cdef void *canonical_dg_vert_parent(void *child, void *parent, int *permutation, int *degree, bint *mem_err)



