from .pynac cimport GConstant

cdef class PynacConstant:
    cdef GConstant* pointer
    cdef GConstant* _object
    cdef _name
