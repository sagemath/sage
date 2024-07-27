# sage_setup: distribution = sagemath-ntl
from sage.libs.ntl.types cimport ZZ_pE_c
from sage.libs.ntl.ntl_ZZ_pEContext cimport ntl_ZZ_pEContext_class
from sage.libs.ntl.ntl_ZZ_pX cimport ntl_ZZ_pX

cdef class ntl_ZZ_pE():
    cdef ZZ_pE_c x
    cdef ntl_ZZ_pEContext_class c
    cdef ntl_ZZ_pX get_as_ZZ_pX(ntl_ZZ_pE self)
    cdef void set_from_ZZ_pX(ntl_ZZ_pE self, ntl_ZZ_pX value) noexcept
    cdef ntl_ZZ_pE _new(self)
