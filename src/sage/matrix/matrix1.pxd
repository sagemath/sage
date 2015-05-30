cimport matrix0

cdef class Matrix(matrix0.Matrix):
    cdef _stack_impl(self, bottom)
