# Macros from arb_poly.h
# See https://github.com/flintlib/flint/issues/1529

from .types cimport *

cdef extern from "flint_wrap.h":
    arb_ptr arb_poly_get_coeff_ptr(arb_poly_t p, ulong n)
