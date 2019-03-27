# distutils: libraries = gmp flint arb

from ..flint.types cimport fmpq_t, ulong

cdef extern from "bernoulli.h":
    void bernoulli_fmpq_ui(fmpq_t b, ulong n)
