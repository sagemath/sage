from sage.libs.ntl.lzz_p cimport *
from sage.libs.ntl.ntl_lzz_pContext cimport ntl_zz_pContext_class

cdef class ntl_zz_p():
    cdef zz_p_c x
    cdef ntl_zz_pContext_class c
    cdef ntl_zz_p _new(ntl_zz_p self)
