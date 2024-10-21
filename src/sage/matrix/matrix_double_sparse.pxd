# sage_setup: distribution = sagemath-modules
from sage.matrix.matrix_generic_sparse cimport Matrix_generic_sparse

cdef class Matrix_double_sparse(Matrix_generic_sparse):
    pass
