from sage.rings.padics.padic_generic_element cimport pAdicGenericElement
from sage.libs.ntl.types cimport ZZ_pX_c, ZZ_pE_c, ZZ_pEX_c, ZZ_p_c, ZZX_c
from sage.libs.ntl.ntl_ZZ_pContext cimport ntl_ZZ_pContext_class
from sage.libs.ntl.ntl_ZZ_pEContext cimport ntl_ZZ_pEContext_class

cdef class pAdicExtElement(pAdicGenericElement):
    cdef int _set_from_list(self, L) except -1
    cdef int _set_from_list_rel(self, L, long relprec) except -1
    cdef int _set_from_list_abs(self, L, long absprec) except -1
    cdef int _set_from_list_both(self, L, long absprec, long relprec) except -1

    cdef int _set_from_ZZX(self, ZZX_c poly) except -1
    cdef int _set_from_ZZX_rel(self, ZZX_c poly, long relprec) except -1
    cdef int _set_from_ZZX_abs(self, ZZX_c poly, long absprec) except -1
    cdef int _set_from_ZZX_both(self, ZZX_c poly, long absprec, long relprec) except -1

    cdef int _set_from_ZZ_pX(self, ZZ_pX_c* poly, ntl_ZZ_pContext_class ctx) except -1
    cdef int _set_from_ZZ_pX_rel(self, ZZ_pX_c* poly, ntl_ZZ_pContext_class ctx, long relprec) except -1
    cdef int _set_from_ZZ_pX_abs(self, ZZ_pX_c* poly, ntl_ZZ_pContext_class ctx, long absprec) except -1
    cdef int _set_from_ZZ_pX_both(self, ZZ_pX_c* poly, ntl_ZZ_pContext_class ctx, long absprec, long relprec) except -1

    cdef int _set_from_ZZ_pE(self, ZZ_pE_c* poly, ntl_ZZ_pEContext_class ctx) except -1
    cdef int _set_from_ZZ_pE_rel(self, ZZ_pE_c* poly, ntl_ZZ_pEContext_class ctx, long relprec) except -1
    cdef int _set_from_ZZ_pE_abs(self, ZZ_pE_c* poly, ntl_ZZ_pEContext_class ctx, long absprec) except -1
    cdef int _set_from_ZZ_pE_both(self, ZZ_pE_c* poly, ntl_ZZ_pEContext_class ctx, long absprec, long relprec) except -1

    cdef int _set_from_ZZ_pEX(self, ZZ_pEX_c* poly, ntl_ZZ_pEContext_class ctx) except -1
    cdef int _set_from_ZZ_pEX_rel(self, ZZ_pEX_c* poly, ntl_ZZ_pEContext_class ctx, long relprec) except -1
    cdef int _set_from_ZZ_pEX_abs(self, ZZ_pEX_c* poly, ntl_ZZ_pEContext_class ctx, long absprec) except -1
    cdef int _set_from_ZZ_pEX_both(self, ZZ_pEX_c* poly, ntl_ZZ_pEContext_class ctx, long absprec, long relprec) except -1

    cdef long _check_ZZ_pContext(self, ntl_ZZ_pContext_class ctx) except -1
    cdef long _check_ZZ_pEContext(self, ntl_ZZ_pEContext_class ctx) except -1

    cdef ext_p_list(self, bint pos)
    cdef ext_p_list_precs(self, bint pos, long prec)
    cdef ZZ_p_c _const_term(self)
