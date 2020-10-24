from .types cimport zz_pContext_c

cdef class ntl_zz_pContext_class(object):
    cdef zz_pContext_c x
    cdef void restore_c(self)
    cdef long p
