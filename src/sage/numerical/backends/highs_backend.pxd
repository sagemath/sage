#*****************************************************************************
#       Copyright (C) 2025 SageMath Developers
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.numerical.backends.generic_backend cimport GenericBackend

cdef class HiGHSBackend(GenericBackend):
    cdef void* highs
    cdef str prob_name
    cdef dict col_name_var
    cdef dict row_name_var
    cdef dict row_data_cache
    cdef int numcols
    cdef int numrows
    cdef void _get_col_bounds(self, int col, double* lb, double* ub) except *
    cdef void _get_row_bounds(self, int row, double* lb, double* ub) except *
    cpdef __copy__(self)
    cpdef get_row_prim(self, int i)
    cpdef double get_row_dual(self, int i) except? -1
    cpdef double get_col_dual(self, int j) except? -1
    cpdef int get_row_stat(self, int i) except? -1
    cpdef int get_col_stat(self, int j) except? -1
    cpdef set_row_stat(self, int i, int stat)
    cpdef set_col_stat(self, int j, int stat)
    cpdef int warm_up(self) noexcept
    cpdef int add_variable_with_type(self, int vtype, lower_bound=*, upper_bound=*, 
                                     obj=*, name=*) except -1