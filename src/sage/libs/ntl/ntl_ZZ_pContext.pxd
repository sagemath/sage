from .types cimport ZZ_pContext_c
from .ntl_ZZ cimport ntl_ZZ
from .types cimport ZZ_c


cdef class ntl_ZZ_pContext_class(object):
    cdef ZZ_pContext_c x
    cdef void restore_c(self)
    cdef ntl_ZZ p
    cdef double p_bits
    cdef object __weakref__
    cpdef void _assert_is_current_modulus(self) except *


cdef class ntl_ZZ_pContext_factory(object):
    cdef object context_dict
    cdef ntl_ZZ_pContext_class make_c(self, ntl_ZZ v)


cdef extern from "ntlwrap.h":
    cdef const ZZ_c& ntl_ZZ_p_current_modulus "ZZ_p::modulus"()
