from sage.libs.ntl.types cimport ZZ_pContext_c
from sage.libs.ntl.ntl_ZZ cimport ntl_ZZ
from sage.libs.ntl.types cimport ZZ_c


cdef class ntl_ZZ_pContext_class():
    cdef ZZ_pContext_c x
    cdef void restore_c(self) noexcept
    cdef ntl_ZZ p
    cdef double p_bits
    cdef object __weakref__
    cpdef void _assert_is_current_modulus(self) except *


cdef class ntl_ZZ_pContext_factory():
    cdef object context_dict
    cdef ntl_ZZ_pContext_class make_c(self, ntl_ZZ v)


cdef extern from "ntlwrap.h":
    cdef const ZZ_c& ntl_ZZ_p_current_modulus "ZZ_p::modulus"()
