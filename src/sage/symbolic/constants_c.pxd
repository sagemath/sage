from ginac cimport GConstant

cdef class PynacConstant:
    cdef GConstant* pointer
    cdef GConstant object
    cdef object _name
