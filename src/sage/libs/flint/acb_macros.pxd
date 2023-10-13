# Macros from acb.h
# See https://github.com/flintlib/flint/issues/152

from .types cimport *

cdef extern from "flint_wrap.h":
    arb_ptr acb_realref(acb_t x)
    # Macro giving a pointer to the real part of x

    arb_ptr acb_imagref(acb_t x)
    # Macro giving a pointer to the imaginary part of x
