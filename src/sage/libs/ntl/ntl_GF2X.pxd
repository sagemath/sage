from .types cimport GF2X_c

cdef class ntl_GF2X(object):
    cdef GF2X_c x
