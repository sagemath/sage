from sage.matrix.matrix_double_dense cimport Matrix_double_dense


cdef class Matrix_real_double_dense(Matrix_double_dense):
    cdef set_unsafe_double(self, Py_ssize_t i, Py_ssize_t j, double value)
    cdef double get_unsafe_double(self, Py_ssize_t i, Py_ssize_t j) noexcept
