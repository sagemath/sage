# Macros from arb.h
# See https://github.com/flintlib/flint/issues/152

from .types cimport *

cdef extern from "flint_wrap.h":
    arf_ptr arb_midref(arb_t x)
    # Macro giving a pointer to the midpoint of x

    mag_ptr arb_radref(arb_t x)
    # Macro giving a pointer to the radius of x
