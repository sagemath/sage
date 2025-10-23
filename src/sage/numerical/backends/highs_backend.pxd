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
    cdef object highs_model
    cdef str prob_name
    cdef dict col_name_var
    cdef dict row_name_var
    cdef int numcols
    cdef int numrows
    cpdef __copy__(self)
