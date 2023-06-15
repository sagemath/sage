# sage_setup: distribution = sagemath-ntl
from .types cimport GF2E_c
from .ntl_GF2EContext cimport ntl_GF2EContext_class

cdef class ntl_GF2E():
    cdef GF2E_c x
    cdef ntl_GF2EContext_class c
    cdef ntl_GF2E _new(self)

