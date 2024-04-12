# sage_setup: distribution = sagemath-linbox
from sage.libs.m4rie cimport mzed_t
from sage.libs.m4ri cimport m4ri_word
from sage.matrix.matrix_dense cimport Matrix_dense


cdef class Matrix_gf2e_dense(Matrix_dense):
    cdef mzed_t *_entries
    cdef object _one
    cdef object _zero
    cdef m4ri_word _zero_word  # m4ri_word representation of _zero

    cpdef Matrix_gf2e_dense _multiply_newton_john(Matrix_gf2e_dense self, Matrix_gf2e_dense right)
    cpdef Matrix_gf2e_dense _multiply_karatsuba(Matrix_gf2e_dense self, Matrix_gf2e_dense right)
    cpdef Matrix_gf2e_dense _multiply_strassen(Matrix_gf2e_dense self, Matrix_gf2e_dense right, cutoff=*)
