# sage_setup: distribution = sagemath-modules

from sage.matrix.matrix_dense cimport Matrix_dense


cdef class Matrix_generic_dense(Matrix_dense):
    cdef list _entries
    cdef Matrix_generic_dense _new(self, Py_ssize_t nrows, Py_ssize_t ncols)
