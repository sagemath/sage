# Macros from acb_poly.h
# See https://github.com/flintlib/flint/issues/1529

from .types cimport *

cdef extern from "flint_wrap.h":
    acb_ptr acb_poly_get_coeff_ptr(acb_poly_t p, ulong n)
