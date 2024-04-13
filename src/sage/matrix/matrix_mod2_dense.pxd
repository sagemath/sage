# sage_setup: distribution = sagemath-linbox
from .matrix_dense cimport Matrix_dense
from sage.libs.m4ri cimport *

cdef class Matrix_mod2_dense(Matrix_dense):
    cdef mzd_t *_entries
    cdef object _one
    cdef object _zero

    cpdef Matrix_mod2_dense _multiply_m4rm(Matrix_mod2_dense self, Matrix_mod2_dense right, int k)
    cpdef Matrix_mod2_dense _multiply_strassen(Matrix_mod2_dense self, Matrix_mod2_dense right, int cutoff)

    # For conversion to other systems
    cpdef _export_as_string(self)
