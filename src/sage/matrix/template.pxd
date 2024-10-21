# sage_setup: distribution = sagemath-modules

from sage.matrix.matrix_dense cimport Matrix_dense


cdef class Matrix_generic_dense(Matrix_dense):
    pass
