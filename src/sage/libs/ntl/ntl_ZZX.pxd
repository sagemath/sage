# sage_setup: distribution = sagemath-ntl
from sage.libs.ntl.types cimport ZZX_c

cdef class ntl_ZZX():
    cdef ZZX_c x
    cdef void setitem_from_int(ntl_ZZX self, long i, int value) noexcept
    cdef int getitem_as_int(ntl_ZZX self, long i) noexcept
