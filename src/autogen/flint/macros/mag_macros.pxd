# Macros from mag.h
# See https://github.com/flintlib/flint/issues/1529

from .types cimport *

cdef extern from "flint_wrap.h":
    long MAG_BITS
