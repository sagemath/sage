# sage_setup: distribution = sagemath-flint
from sage.libs.gmp.types cimport mpz_t
from sage.modules.vector_rational_sparse cimport mpq_vector
from sage.matrix.matrix_sparse cimport Matrix_sparse


cdef class Matrix_rational_sparse(Matrix_sparse):
    cdef mpq_vector* _matrix

    cdef int mpz_denom(self, mpz_t d) except -1
    cdef int mpz_height(self, mpz_t height) except -1
