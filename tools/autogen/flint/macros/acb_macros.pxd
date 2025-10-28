# Macros from acb.h
# See https://github.com/flintlib/flint/issues/1529

from .types cimport *

cdef extern from "flint_wrap.h":
    arb_ptr acb_realref(acb_t x)
    arb_ptr acb_imagref(acb_t x)
