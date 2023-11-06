from sage.matrix.matrix cimport Matrix


cdef class MatrixWindow:
    cdef Py_ssize_t _row, _col, _nrows, _ncols
    cdef Matrix _matrix
    cdef object _cached_zero

    # YOU *REALLY SHOULD* OVERRIDE THESE:
    cpdef add(MatrixWindow self, MatrixWindow A) noexcept
    cpdef subtract(MatrixWindow self, MatrixWindow A) noexcept
    cpdef set_to_sum(MatrixWindow self, MatrixWindow A, MatrixWindow B) noexcept
    cpdef set_to_diff(MatrixWindow self, MatrixWindow A, MatrixWindow B) noexcept
    cpdef set_to_prod(MatrixWindow self, MatrixWindow A, MatrixWindow B) noexcept
    cpdef add_prod(MatrixWindow self, MatrixWindow A, MatrixWindow B) noexcept
    cpdef subtract_prod(MatrixWindow self, MatrixWindow A, MatrixWindow B) noexcept

    cpdef bint element_is_zero(MatrixWindow self, Py_ssize_t i, Py_ssize_t j) noexcept
    cpdef set_to(MatrixWindow self, MatrixWindow A) noexcept
    cpdef set_to_zero(MatrixWindow self) noexcept

    # FOR BETTER SPEED, OVERRIDE ANY SUBSET OF THESE (OPTIONAL):
    cpdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, x) noexcept
    cpdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j) noexcept
    cpdef to_matrix(MatrixWindow self) noexcept
    cpdef new_empty_window(MatrixWindow self, Py_ssize_t nrows, Py_ssize_t ncols) noexcept

    # NO BENEFIT TO OVERRIDING THESE:
    cpdef MatrixWindow matrix_window(MatrixWindow self, Py_ssize_t row, Py_ssize_t col,
                                     Py_ssize_t n_rows, Py_ssize_t n_cols) noexcept
    cpdef MatrixWindow new_matrix_window(MatrixWindow self, Matrix matrix,
                                        Py_ssize_t row, Py_ssize_t col,
                                         Py_ssize_t n_rows, Py_ssize_t n_cols) noexcept
    cpdef matrix(MatrixWindow self) noexcept
    cpdef swap_rows(MatrixWindow self, Py_ssize_t a, Py_ssize_t b) noexcept
    cdef object _zero(self) noexcept
