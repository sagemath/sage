from .matrix_dense cimport Matrix_dense
cimport numpy as cnumpy

cdef class Matrix_double_dense(Matrix_dense):
    cdef object _numpy_dtype
    # cdef cnumpy.NPY_TYPES _numpy_dtypeint
    cdef int _numpy_dtypeint
    cdef object _python_dtype
    cdef object _sage_dtype
    cdef object _sage_vector_dtype
    cdef Matrix_double_dense _new(self, int nrows=*, int ncols=*)
    cdef cnumpy.ndarray _matrix_numpy
