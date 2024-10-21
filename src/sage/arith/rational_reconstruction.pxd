# sage_setup: distribution = sagemath-categories
from sage.libs.gmp.types cimport mpz_t, mpq_t

cdef int mpq_rational_reconstruction(mpq_t answer, mpz_t a, mpz_t m) except -1
