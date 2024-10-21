# sage_setup: distribution = sagemath-categories
from sage.libs.gmp.all cimport mpz_t

cdef mpz_stirling_s2(mpz_t s, unsigned long n, unsigned long k)

cdef list convert(Py_ssize_t* f, Py_ssize_t n)
