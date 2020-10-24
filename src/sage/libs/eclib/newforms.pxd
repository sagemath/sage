from ..eclib cimport newforms

cdef class ECModularSymbol:
    cdef newforms* nfs
    cdef int n
    cdef object _E
    cdef int sign
