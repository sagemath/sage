from sage.libs.flint.types cimport fmpz_t, fmpq_mat_t
from sage.matrix.matrix_dense cimport Matrix_dense

cdef class Matrix_rational_dense(Matrix_dense):
    cdef fmpq_mat_t _matrix

    cdef int fmpz_denom(self, fmpz_t d) except -1
    cdef int fmpz_height(self, fmpz_t height) except -1
#    cdef int _rescale(self, mpq_t a) except -1

    cdef _pickle_version0(self)
    cdef _unpickle_version0(self, data)
    cpdef _export_as_string(self, int base=?)

    cdef _add_ui_unsafe_assuming_int(self, Py_ssize_t i, Py_ssize_t j, unsigned long int n)
    cdef _sub_ui_unsafe_assuming_int(self, Py_ssize_t i, Py_ssize_t j, unsigned long int n)

    cdef inline Matrix_rational_dense _new_matrix(self, Py_ssize_t nrows, Py_ssize_t ncols)

cdef class MatrixWindow:
    cdef Matrix_rational_dense _matrix
    cdef int _row, _col, _nrows, _ncols
