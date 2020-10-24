from .types cimport GF2_c

cdef class ntl_GF2(object):
    cdef GF2_c x
