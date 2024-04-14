# sage_setup: distribution = sagemath-ntl
from sage.libs.ntl.types cimport ZZ_pEX_c
from sage.libs.ntl.ntl_ZZ_pEContext cimport ntl_ZZ_pEContext_class

cdef class ntl_ZZ_pEX():
    cdef ZZ_pEX_c x
    cdef ntl_ZZ_pEContext_class c
    #cdef void setitem_from_int(ntl_ZZ_pX self, long i, int value)
    #cdef int getitem_as_int(ntl_ZZ_pX self, long i)
    cdef ntl_ZZ_pEX _new(self)
