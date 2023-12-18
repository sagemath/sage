from sage.matrix.matrix0 cimport Matrix as Matrix0

cdef class Matrix(Matrix0):
    cdef _stack_impl(self, bottom) noexcept

    cpdef row_ambient_module(self, base_ring=*, sparse=*) noexcept
    cpdef column_ambient_module(self, base_ring=*, sparse=*) noexcept
