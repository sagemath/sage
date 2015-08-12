include "decl.pxi"

from sage.libs.ntl.ntl_ZZ_pContext cimport ntl_ZZ_pContext_class

cdef class ntl_ZZ_p:
    cdef ZZ_p_c x
    cdef ntl_ZZ_pContext_class c
    cdef public int get_as_int(ntl_ZZ_p self)
    cdef public void set_from_int(ntl_ZZ_p self, int value)
    cdef ntl_ZZ_p _new(self)
