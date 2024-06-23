# sage_setup: distribution = sagemath-pari
from sage.rings.integer cimport Integer
from sage.rings.padics.pow_computer cimport PowComputer_class
from sage.libs.gmp.mpz cimport mpz_t, mpq_t

cdef long get_ordp(x, PowComputer_class prime_pow) except? -10000
cdef long get_preccap(x, PowComputer_class prime_pow) except? -10000
cdef long comb_prec(iprec, long prec) except? -10000
cdef int _process_args_and_kwds(long *aprec, long *rprec, args, kwds, bint absolute, PowComputer_class prime_pow) except -1

cdef long cconv_mpq_t_shared(mpz_t out, mpq_t x, long prec, bint absolute, PowComputer_class prime_pow) except? -10000
cdef int cconv_mpq_t_out_shared(mpq_t out, mpz_t x, long valshift, long prec, PowComputer_class prime_pow) except -1
cdef int cconv_shared(mpz_t out, x, long prec, long valshift, PowComputer_class prime_pow) except -2
cdef long cconv_mpz_t_shared(mpz_t out, mpz_t x, long prec, bint absolute, PowComputer_class prime_pow) except -2
cdef int cconv_mpz_t_out_shared(mpz_t out, mpz_t x, long valshift, long prec, PowComputer_class prime_pow) except -1
