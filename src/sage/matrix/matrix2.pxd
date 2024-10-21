# sage_setup: distribution = sagemath-modules
"""
Generic matrices
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.matrix.matrix1 cimport Matrix as Matrix1

cdef class Matrix(Matrix1):
    cdef _det_by_minors(self, Py_ssize_t level)
    cdef _pf_bfl(self)
    cdef bint _is_positive_definite_or_semidefinite(self, bint semi) except -1
    cdef tuple _block_ldlt(self, bint classical)
    cpdef _echelon(self, str algorithm)
    cpdef _echelon_in_place(self, str algorithm)
    cpdef matrix_window(self, Py_ssize_t row=*, Py_ssize_t col=*, Py_ssize_t nrows=*, Py_ssize_t ncols=*, bint check=*)
