include "decl.pxi"

cdef class ntl_ZZ:
    cdef ZZ_c x
    cdef int get_as_int(ntl_ZZ self)
    cdef void set_from_int(ntl_ZZ self, int value)

cdef void PyLong_to_ZZ(ZZ_c* z, value)
