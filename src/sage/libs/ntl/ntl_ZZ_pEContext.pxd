include "decl.pxi"
include "sage/ext/cdefs.pxi"

from sage.libs.ntl.ntl_ZZ_pContext cimport ntl_ZZ_pContext_class

cdef struct ZZ_pEContext_ptrs:
    ZZ_pEContext_c *zzpec
    ZZ_pContext_c *zzpc

cdef class ntl_ZZ_pEContext_class:
    cdef ZZ_pEContext_ptrs ptrs
    cdef ZZ_pEContext_c x
    cdef ntl_ZZ_pContext_class pc
    cdef void restore_c(self)
    cdef f
