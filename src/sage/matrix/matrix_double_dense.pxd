# sage_setup: distribution = sagemath-modules
from sage.matrix.matrix_numpy_dense cimport Matrix_numpy_dense


cdef class Matrix_double_dense(Matrix_numpy_dense):

    pass
