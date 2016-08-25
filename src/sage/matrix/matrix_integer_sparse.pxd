include 'sage/modules/vector_integer_sparse_h.pxi'

from sage.libs.gmp.types cimport mpz_t
from sage.ext.mod_int cimport *
from .matrix_sparse cimport Matrix_sparse

cdef class Matrix_integer_sparse(Matrix_sparse):
    cdef mpz_vector* _matrix
    cdef int _initialized

    cdef _mod_int_c(self, mod_int p)
