# Macros from fmpz_poly.h
# See https://github.com/flintlib/flint/issues/1529

from .types cimport *

cdef extern from "flint_wrap.h":
    bint COEFF_IS_MPZ(fmpz f)
