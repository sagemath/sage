# distutils: libraries = flint
# distutils: depends = flint/arith.h

from sage.libs.flint.types cimport fmpz_t, fmpq_t, ulong

# flint/arith.h
cdef extern from "flint_wrap.h":
    void arith_bell_number(fmpz_t b, ulong n)
    void arith_bernoulli_number(fmpq_t x, ulong n)
    void arith_euler_number ( fmpz_t res , ulong n )
    void arith_number_of_partitions(fmpz_t x, ulong n)
    void arith_dedekind_sum(fmpq_t, fmpz_t, fmpz_t)
    void arith_harmonic_number(fmpq_t, unsigned long n)
