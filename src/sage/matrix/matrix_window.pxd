# sage_setup: distribution = sagemath-modules
from sage.matrix.matrix cimport Matrix


cdef class MatrixWindow:
    cdef Py_ssize_t _row, _col, _nrows, _ncols
    cdef Matrix _matrix
    cdef object _cached_zero

    # YOU *REALLY SHOULD* OVERRIDE THESE:
    cpdef add(MatrixWindow self, MatrixWindow A)
    cpdef subtract(MatrixWindow self, MatrixWindow A)
    cpdef set_to_sum(MatrixWindow self, MatrixWindow A, MatrixWindow B)
    cpdef set_to_diff(MatrixWindow self, MatrixWindow A, MatrixWindow B)
    cpdef set_to_prod(MatrixWindow self, MatrixWindow A, MatrixWindow B)
    cpdef add_prod(MatrixWindow self, MatrixWindow A, MatrixWindow B)
    cpdef subtract_prod(MatrixWindow self, MatrixWindow A, MatrixWindow B)

    cpdef bint element_is_zero(MatrixWindow self, Py_ssize_t i, Py_ssize_t j) noexcept
    cpdef set_to(MatrixWindow self, MatrixWindow A)
    cpdef set_to_zero(MatrixWindow self)

    # FOR BETTER SPEED, OVERRIDE ANY SUBSET OF THESE (OPTIONAL):
    cpdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, x)
    cpdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j)
    cpdef to_matrix(MatrixWindow self)
    cpdef new_empty_window(MatrixWindow self, Py_ssize_t nrows, Py_ssize_t ncols)

    # NO BENEFIT TO OVERRIDING THESE:
    cpdef MatrixWindow matrix_window(MatrixWindow self, Py_ssize_t row, Py_ssize_t col,
                                     Py_ssize_t n_rows, Py_ssize_t n_cols)
    cpdef MatrixWindow new_matrix_window(MatrixWindow self, Matrix matrix,
                                        Py_ssize_t row, Py_ssize_t col,
                                         Py_ssize_t n_rows, Py_ssize_t n_cols)
    cpdef matrix(MatrixWindow self)
    cpdef swap_rows(MatrixWindow self, Py_ssize_t a, Py_ssize_t b)
    cdef object _zero(self)
