from sage.libs.ntl.types cimport ZZ_c

cdef class ntl_ZZ():
    cdef ZZ_c x
    cdef int get_as_int(ntl_ZZ self) noexcept
    cdef void set_from_int(ntl_ZZ self, int value) noexcept
