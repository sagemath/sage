from .types cimport GF2EX_c
from .ntl_GF2EContext cimport ntl_GF2EContext_class
from .ntl_GF2E cimport ntl_GF2E

cdef class ntl_GF2EX(object):
    cdef GF2EX_c x
    cdef ntl_GF2EContext_class c
    cdef ntl_GF2E _new_element(self)
    cdef ntl_GF2EX _new(self)
