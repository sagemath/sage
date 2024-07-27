# sage_setup: distribution = sagemath-symbolics

from sage.matrix.matrix_generic_sparse cimport Matrix_generic_sparse


cdef class Matrix_symbolic_sparse(Matrix_generic_sparse):
    pass
