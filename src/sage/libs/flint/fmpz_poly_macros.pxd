# Macros from fmpz_poly.h
# See https://github.com/flintlib/flint/issues/1529

from .types cimport *

cdef extern from "flint_wrap.h":
    fmpz * fmpz_poly_get_coeff_ptr(fmpz_poly_t p, ulong n)
