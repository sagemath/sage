from .types cimport mat_GF2_c
from .ntl_GF2 cimport ntl_GF2

cdef class ntl_mat_GF2(object):
    cdef mat_GF2_c x
    cdef ntl_GF2 _new_element(self)
    cdef ntl_mat_GF2 _new(self)
