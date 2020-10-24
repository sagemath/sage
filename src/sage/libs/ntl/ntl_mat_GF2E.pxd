from .types cimport mat_GF2E_c
from .ntl_GF2EContext cimport ntl_GF2EContext_class
from .ntl_GF2E cimport ntl_GF2E

cdef class ntl_mat_GF2E(object):
    cdef mat_GF2E_c x
    cdef ntl_GF2EContext_class c
    cdef ntl_GF2E _new_element(self)
    cdef ntl_mat_GF2E _new(self)
