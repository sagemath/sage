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

from sage.groups.perm_gps.partn_ref.automorphism_group_canonical_label cimport (
    get_aut_gp_and_can_lab, aut_gp_and_can_lab, agcl_work_space,
    allocate_agcl_output, deallocate_agcl_output,
    allocate_agcl_work_space, deallocate_agcl_work_space)

cdef class MatrixStruct:
    cdef list symbol_structs
    cdef object matrix
    cdef int degree
    cdef int nwords
    cdef list symbols
    cdef int nsymbols
    cdef PartitionStack *temp_col_ps
    cdef aut_gp_and_can_lab *output

cdef int refine_matrix(PartitionStack *, void *, int *, int) noexcept
cdef int compare_matrices(int *, int *, void *, void *, int) noexcept
cdef bint all_matrix_children_are_equivalent(PartitionStack *, void *) noexcept
