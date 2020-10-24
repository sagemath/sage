from sage.libs.ntl.types cimport ZZ_c

cdef class ntl_ZZ(object):
    cdef ZZ_c x
    cdef int get_as_int(ntl_ZZ self)
    cdef void set_from_int(ntl_ZZ self, int value)
