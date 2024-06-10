from sage.libs.ntl.types cimport ZZ_p_c
from sage.libs.ntl.ntl_ZZ_pContext cimport ntl_ZZ_pContext_class

cdef class ntl_ZZ_p():
    cdef ZZ_p_c x
    cdef ntl_ZZ_pContext_class c
    cdef int get_as_int(ntl_ZZ_p self) noexcept
    cdef void set_from_int(ntl_ZZ_p self, int value) noexcept
    cdef ntl_ZZ_p _new(self)
