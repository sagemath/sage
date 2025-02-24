from sage.libs.ntl.types cimport GF2_c

cdef class ntl_GF2():
    cdef GF2_c x
