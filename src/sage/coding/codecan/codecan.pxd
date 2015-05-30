from sage.groups.perm_gps.permgroup_element cimport PermutationGroupElement
from sage.groups.semimonomial_transformations.semimonomial_transformation cimport SemimonomialTransformation
from sage.modules.free_module_element cimport FreeModuleElement
from sage.groups.perm_gps.partn_ref2.refinement_generic cimport *

cdef class InnerGroup:
    cdef int rank
    cdef OrbitPartition * row_partition
    cdef int frob_pow
    cdef bint permutational_only

    # for the transporter element computation we need
    cdef SemimonomialTransformation transporter
    cdef bint compute_transporter

    cdef inline int get_rep(self, int pos)
    cdef inline int join_rows(self, int rep1, int rep2)

    cdef InnerGroup _new_c(self)
    cdef void copy_from(self, InnerGroup other)
    cdef bint has_semilinear_action(self)
    cdef minimize_by_row_mult(self, FreeModuleElement v)
    cdef minimize_matrix_col(self, object m, int pos, list fixed_minimized_cols,
                             bint *group_changed)
    cdef void gaussian_elimination(self, object m, int pos, int pivot, list nz_pos)
    cdef void minimize_by_frobenius(self, object v, int *applied_frob, int *stab_pow)

    cdef SemimonomialTransformation get_transporter(self)

    cdef bint has_semilinear_action(self)
    cpdef int get_frob_pow(self)
    cpdef column_blocks(self, mat)

cdef class PartitionRefinementLinearCode(PartitionRefinement_generic):
    cdef int _k, _q
    cdef long *_hyp_refine_vals_scratch
    cdef object _inner_group_stabilizer_order
    cdef bitset_t *_hyp2points # hyperplanes to points
    cdef bitset_t *_points2hyp # points to hyperplanes, transpose of _hyp2points
    cdef PartitionStack *_hyp_part
    cdef object _matrix, _root_matrix
    cdef InnerGroup _inner_group
    cdef dict _stored_states

    # store the result of the refinements
    cdef _BestValStore _supp_refine_vals, _point_refine_vals, _hyp_refine_vals
    cdef int _nr_of_supp_refine_calls, _nr_of_point_refine_calls, _nr_of_hyp_refine_calls

    # the elements we want to compute
    cdef object _best_candidate
    cdef SemimonomialTransformation _transporter
    cdef list _autom_group_generators

    # specialized refine methods, called in refine
    cdef bint _inner_min_refine(self, bint *inner_stab_changed, bint *changed_partition)
    cdef bint _point_refine(self, bint *inner_stab_changed, bint *changed_partition)
    cdef bint _hyp_refine(self, bint *changed_partition)

    # some additional methods
    cdef _compute_group_element(self, SemimonomialTransformation trans, str algorithm_type)
    cdef _init_point_hyperplane_incidence(self)