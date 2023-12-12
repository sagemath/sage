# Macros from arb.h
# See https://github.com/flintlib/flint/issues/1529

from .types cimport *

cdef extern from "flint_wrap.h":
    arf_ptr arb_midref(arb_t x)
    mag_ptr arb_radref(arb_t x)
