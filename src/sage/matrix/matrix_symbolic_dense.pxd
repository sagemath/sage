# sage_setup: distribution = sagemath-symbolics
from sage.matrix.matrix_generic_dense cimport Matrix_generic_dense


cdef class Matrix_symbolic_dense(Matrix_generic_dense):
    pass
