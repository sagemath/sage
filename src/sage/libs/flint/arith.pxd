from sage.libs.flint.types cimport fmpz_t, fmpq_t, ulong

cdef extern from "flint/arith.h":
    void arith_bell_number(fmpz_t b, ulong n)
    void arith_number_of_partitions(fmpz_t x, ulong n)
    void arith_dedekind_sum(fmpq_t, fmpz_t, fmpz_t)
