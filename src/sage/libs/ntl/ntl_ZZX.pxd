from .types cimport ZZX_c

cdef class ntl_ZZX(object):
    cdef ZZX_c x
    cdef void setitem_from_int(ntl_ZZX self, long i, int value)
    cdef int getitem_as_int(ntl_ZZX self, long i)
