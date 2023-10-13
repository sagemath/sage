# Macros from fmpz_poly.h
# See https://github.com/flintlib/flint/issues/152

from .types cimport *

cdef extern from "flint_wrap.h":
    fmpz * fmpz_poly_get_coeff_ptr(fmpz_poly_t p, ulong n)
    # Macro giving a pointer to the n-th coefficient of p (or NULL)
