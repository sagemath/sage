# distutils: libraries = flint
# distutils: depends = flint/fmpz_factor.h

from sage.libs.flint.types cimport *

# flint/fmpz_factor.h
cdef extern from "flint_wrap.h":
    void fmpz_factor_clear(fmpz_factor_t)
    void fmpz_factor_init(fmpz_factor_t)
    void fmpz_factor(fmpz_factor_t, const fmpz_t)

cdef fmpz_factor_to_pairlist(const fmpz_factor_t)
