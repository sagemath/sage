# sage_setup: distribution = sagemath-ntl
from sage.libs.ntl.types cimport GF2EContext_c
from sage.libs.ntl.ntl_GF2X cimport ntl_GF2X

cdef class ntl_GF2EContext_class():
    cdef GF2EContext_c x
    cdef ntl_GF2X m
    cdef void restore_c(self) noexcept
    cdef object __weakref__
