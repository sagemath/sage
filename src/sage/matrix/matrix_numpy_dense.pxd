from .matrix_dense cimport Matrix_dense
cimport numpy as cnumpy


cdef class Matrix_numpy_dense(Matrix_dense):
    cdef object _numpy_dtype
    cdef int _numpy_dtypeint
    cdef object _python_dtype
    cdef object _sage_dtype
    cdef object _sage_vector_dtype
    cdef Matrix_numpy_dense _new(self, int nrows=*, int ncols=*)
    cdef cnumpy.ndarray _matrix_numpy
