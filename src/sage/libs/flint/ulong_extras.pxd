# distutils: libraries = flint
# distutils: depends = flint/ulong_extras.h

from sage.libs.flint.types cimport n_factor_t

# flint/ulong_extras.h
cdef extern from "flint_wrap.h":
    cdef int n_jacobi(long x, unsigned long y)

    cdef int n_is_prime(unsigned long n)

    cdef unsigned long n_gcd(long x, long y)

    cdef void n_factor(n_factor_t * factors, unsigned long n, int proved)
    cdef void n_factor_init(n_factor_t * factors)
