from sage.libs.gmp.types cimport mpz_srcptr, mpq_t

cdef int mpq_rational_reconstruction(mpq_t answer, mpz_srcptr a, mpz_srcptr m) except -1
