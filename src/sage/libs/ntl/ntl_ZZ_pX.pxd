# sage_setup: distribution = sagemath-ntl
from sage.libs.ntl.ZZ_pX cimport *
from sage.libs.ntl.ntl_ZZ_pContext cimport ntl_ZZ_pContext_class
from sage.rings.integer cimport Integer

cdef class ntl_ZZ_pX():
    cdef ZZ_pX_c x
    cdef ntl_ZZ_pContext_class c
    cdef void setitem_from_int(ntl_ZZ_pX self, long i, int value) noexcept
    cdef int getitem_as_int(ntl_ZZ_pX self, long i) noexcept
    cdef ntl_ZZ_pX _new(self)
    cdef ntl_ZZ_pX _pow(ntl_ZZ_pX self, long exp)
    cdef ntl_ZZ_pX _powmod(ntl_ZZ_pX self, Integer exp, ntl_ZZ_pX modulus)

cdef class ntl_ZZ_pX_Modulus():
    cdef ZZ_pX_Modulus_c x
    cdef ntl_ZZ_pX poly
