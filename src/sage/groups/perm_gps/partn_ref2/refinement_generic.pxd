#*******************************************************************************
#       Copyright (C) 2012 Thomas Feulner <thomas.feulner@uni-bayreuth.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*******************************************************************************

from sage.data_structures.bitset_base cimport *
from sage.groups.perm_gps.partn_ref.data_structures cimport OrbitPartition, PartitionStack
from sage.libs.gap.element cimport GapElement, GapElement_Permutation
from sage.structure.parent cimport Parent

cdef extern from "refinement_generic.h":
    cdef long *global_refine_vals_array
    cdef int my_comp_func(void *a, void *b)
    cdef int BACKTRACK_WITHLATEX_DEBUG

cdef extern from *:
    void *qsort(void *base, size_t nmemb, size_t size,
                  int(*compar)(void *, void *))

cdef tuple PS_refinement(PartitionStack * part, long *refine_vals, long *best, 
                         int begin, int end,
                         bint * cand_initialized, bint *changed_partition)

cdef class _BestValStore:
    cdef int default_data_length
    cdef int storage_length
    cdef long *values
    cdef long * get_row(self, int i) noexcept

cdef class LabelledBranching:
    cdef int n
    cdef int count
    cdef int *father
    cdef int *act_perm
    cdef GapElement group, ClosureGroup
    cdef Parent sym_gp
    cdef bint has_empty_intersection(self, PartitionStack * part) noexcept
    cpdef add_gen(self, GapElement_Permutation gen)

cdef class PartitionRefinement_generic:
    cdef PartitionStack * _part
    cdef list _inner_min_order_best
    cdef list _fixed_minimized
    cdef list _fixed_not_minimized
    cdef int _n
    cdef long* _refine_vals_scratch
    cdef bint _is_candidate_initialized
    cdef list _allowance_best
    cdef int _nr_of_inner_min_unmin_calls

    cdef LabelledBranching _known_automorphisms
    cdef GapElement_Permutation _to_best, _to_best_inverse

    # the following allow us to debug the program via latex
    cdef object _latex_debug_string
    cdef void _init_latex(self) noexcept
    cdef void _finish_latex(self) noexcept
    cdef void _latex_new_lvl(self) noexcept
    cdef void _latex_finish_lvl(self) noexcept

     # methods which have to be implemented in derived classes
    cdef bint _inner_min_(self, int pos, bint* inner_group_changed) noexcept
    cdef bint _refine(self, bint *part_changed, bint inner_group_changed, bint first_step) noexcept
    cdef tuple _store_state_(self)
    cdef void _restore_state_(self, tuple act_state) noexcept
    cdef void _store_best_(self) noexcept
    cdef bint _minimization_allowed_on_col(self, int pos) noexcept
    cdef void _latex_act_node(self, str comment=*, int printlvl=*) noexcept # only if you want to allow
                                            # some debugging output

    # methods used in the main algorithm
    cdef void _init_partition_stack(self, list partition) noexcept
    cdef void _start_Sn_backtrack(self) noexcept
    cdef void _backtrack(self, bint first_step=?) noexcept
    cdef void _leaf_computations(self) noexcept
    cdef bint _cut_by_known_automs(self) noexcept
    cdef bint _inner_min_unminimized(self, bint *inner_group_changed) noexcept
    cdef bint _one_refinement(self, long *best, int begin, int end,
                              bint* inner_group_changed, bint* changed_partition,
                              str refine_name) noexcept
    cdef int len(self) noexcept


cdef int PS_first_smallest_PR(PartitionStack *PS, bitset_t b, int *second_pos=?,
        PartitionRefinement_generic partn_ref_alg=?) noexcept
