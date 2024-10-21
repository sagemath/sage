# sage_setup: distribution = sagemath-modules
from sage.matrix.matrix_generic_dense cimport Matrix_generic_dense

cdef class Matrix_polynomial_dense(Matrix_generic_dense):
    pass
