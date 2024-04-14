# sage_setup: distribution = sagemath-modules
from .matrix cimport Matrix

cdef class Matrix_sparse(Matrix):
    pass
