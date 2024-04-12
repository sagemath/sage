# sage_setup: distribution = sagemath-ntl
"""
Dense Matrix Template for C/C++ Library Interfaces
"""

from sage.ext.mod_int cimport *
from sage.matrix.matrix_dense cimport Matrix_dense

cdef class Matrix_modn_dense_template(Matrix_dense):
    cdef celement **_matrix
    cdef celement *_entries
    cdef mod_int p
    cdef xgcd_eliminate (self, celement * row1, celement* row2, Py_ssize_t start_col)
    cpdef _export_as_string(self)
    cdef int _copy_row_to_mod_int_array(self, mod_int *to, Py_ssize_t i) noexcept
