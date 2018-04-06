from sage.libs.gmp.types cimport mpz_t
from sage.modules.vector_modn_sparse cimport *

from sage.matrix.matrix_integer_dense cimport mod_int

cdef class Linbox_modn_sparse:
    cdef c_vector_modint *rows
    cdef size_t nrows
    cdef size_t ncols
    cdef unsigned int modulus

    cdef set(self, int modulus, size_t nrows, size_t ncols, c_vector_modint *rows)
    cdef object rank(self, int gauss)
    cdef void solve(self, c_vector_modint **x, c_vector_modint *b, int method)

