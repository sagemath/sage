from sage.libs.ntl.types cimport ZZ_pContext_c, ZZ_pEContext_c
from sage.libs.ntl.ntl_ZZ_pContext cimport ntl_ZZ_pContext_class
from sage.libs.ntl.ntl_ZZ_pX cimport ntl_ZZ_pX
from sage.libs.ntl.types cimport ZZ_pX_Modulus_c


cdef struct ZZ_pEContext_ptrs:
    ZZ_pEContext_c *zzpec
    ZZ_pContext_c *zzpc


cdef class ntl_ZZ_pEContext_class():
    cdef ZZ_pEContext_ptrs ptrs
    cdef ZZ_pEContext_c x
    cdef ntl_ZZ_pContext_class pc
    cdef void restore_c(self) noexcept
    cdef ntl_ZZ_pX f
    cpdef void _assert_is_current_modulus(self) except *


cdef extern from "ntlwrap.h":
    cdef ZZ_pX_Modulus_c& ZZ_pE_current_modulus "ZZ_pE::modulus"()
