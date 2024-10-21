# sage_setup: distribution = sagemath-modules

from sage.matrix.matrix cimport Matrix

cdef class Matrix_sparse(Matrix):
    pass
