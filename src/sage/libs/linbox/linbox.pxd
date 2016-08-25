# distutils: language = c++ 
# distutils: extra_compile_args = -std=c++98 

from sage.libs.gmp.types cimport mpz_t

include 'sage/modules/vector_modn_sparse_h.pxi'

from sage.matrix.matrix_integer_dense cimport mod_int

cdef class Linbox_modn_sparse:
    cdef c_vector_modint *rows
    cdef size_t nrows
    cdef size_t ncols
    cdef unsigned int modulus

    cdef set(self, int modulus, size_t nrows, size_t ncols, c_vector_modint *rows)
    cdef object rank(self, int gauss)
    cdef void solve(self, c_vector_modint **x, c_vector_modint *b, int method)


cdef class Linbox_integer_dense:
    cdef mpz_t** matrix
    cdef size_t nrows, ncols

    cdef set(self, mpz_t** matrix, size_t nrows, size_t ncols)
    cdef matrix_matrix_multiply(self,
                                mpz_t **ans,
                                mpz_t **B,
                                size_t B_nr, size_t B_nc)
    cdef unsigned long rank(self) except -1

