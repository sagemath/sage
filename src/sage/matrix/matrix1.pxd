from .matrix0 cimport Matrix as Matrix0

cdef class Matrix(Matrix0):
    cdef _stack_impl(self, bottom)
