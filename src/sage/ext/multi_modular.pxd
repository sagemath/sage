# distutils: libraries = gmp
# distutils: extra_compile_args = -std=c99

include "sage/ext/cdefs.pxi"

from sage.ext.mod_int cimport *
from sage.rings.integer import Integer
from sage.rings.integer cimport Integer

cdef class MultiModularBasis_base:
    cdef int      n
    cdef mod_int* moduli
    cdef mpz_t*   partial_products
    cdef mod_int* C # precomputed values for CRT
    cdef mpz_t    product
    cdef mpz_t    half_product
    cdef unsigned long _l_bound
    cdef unsigned long _u_bound
    cdef unsigned long _num_primes

    cdef mod_int _new_random_prime(self, set known_primes) except 1
    cdef mod_int last_prime(self)
    cdef _realloc_to_new_count(self, new_count)
    cdef int _extend_moduli_to_height_c(self, mpz_t height) except -1
    cdef void _refresh_products(self, int start)
    cdef void _refresh_prod(self)
    cdef void _refresh_precomputations(self, int start) except *
    cdef int min_moduli_count(self, mpz_t height) except -1

    cdef int mpz_reduce_tail(self, mpz_t z, mod_int* b, int offset, int len) except -1
    cdef int mpz_reduce_vec_tail(self, mpz_t* z, mod_int** b, int vn, int offset, int len) except -1
    cdef int mpz_crt_tail(self, mpz_t z, mod_int* b, int offset, int len) except -1
    cdef int mpz_crt_vec_tail(self, mpz_t* z, mod_int** b, int vn, int offset, int len) except -1

cdef class MultiModularBasis(MultiModularBasis_base):
    cdef int mpz_reduce(self, mpz_t z, mod_int* b) except -1
    cdef int mpz_reduce_vec(self, mpz_t* z, mod_int** b, int vn) except -1
    cdef int mpz_crt(self, mpz_t z, mod_int* b) except -1
    cdef int mpz_crt_vec(self, mpz_t* z, mod_int** b, int vn) except -1

cdef class MutableMultiModularBasis(MultiModularBasis):
    cdef mod_int __last_prime
    cdef mod_int next_prime_c(self) except -1
    cdef mod_int replace_prime_c(self, int i) except -1

