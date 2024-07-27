# sage_setup: distribution = sagemath-modules
from sage.matrix.matrix0 cimport Matrix as Matrix0

cdef class Matrix(Matrix0):
    cdef _stack_impl(self, bottom)

    cpdef row_ambient_module(self, base_ring=*, sparse=*)
    cpdef column_ambient_module(self, base_ring=*, sparse=*)
