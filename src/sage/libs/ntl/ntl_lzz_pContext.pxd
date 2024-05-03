from sage.libs.ntl.types cimport zz_pContext_c

cdef class ntl_zz_pContext_class():
    cdef zz_pContext_c x
    cdef void restore_c(self) noexcept
    cdef long p
