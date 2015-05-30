
#*****************************************************************************
#      Copyright (C) 2006 - 2011 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

include 'sage/ext/cdefs.pxi'
include 'data_structures_pxd.pxi' # includes bitsets

from sage.rings.integer cimport Integer
from automorphism_group_canonical_label cimport \
    get_aut_gp_and_can_lab, aut_gp_and_can_lab, agcl_work_space, \
    allocate_agcl_output, deallocate_agcl_output, \
    allocate_agcl_work_space, deallocate_agcl_work_space
from refinement_binary cimport NonlinearBinaryCodeStruct, refine_by_bip_degree
from refinement_binary cimport all_children_are_equivalent as all_binary_children_are_equivalent
from double_coset cimport double_coset

cdef class MatrixStruct:
    cdef list symbol_structs
    cdef object matrix
    cdef int degree
    cdef int nwords
    cdef list symbols
    cdef int nsymbols
    cdef PartitionStack *temp_col_ps
    cdef aut_gp_and_can_lab *output

cdef int refine_matrix(PartitionStack *, void *, int *, int)
cdef int compare_matrices(int *, int *, void *, void *, int)
cdef bint all_matrix_children_are_equivalent(PartitionStack *, void *)

