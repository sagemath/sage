cdef class GSLDoubleArray:
    cdef size_t n
    cdef size_t stride
    cdef double * data
