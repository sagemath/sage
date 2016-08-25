from sage.libs.gmp.types cimport *
from .matrix_dense cimport Matrix_dense
from sage.rings.integer cimport Integer
from sage.ext.mod_int cimport *

cdef class Matrix_rational_dense(Matrix_dense):

    cdef mpq_t tmp
    cdef mpq_t *_entries
    cdef mpq_t ** _matrix
    cdef object __pivots

    cdef int mpz_denom(self, mpz_t d) except -1
    cdef int mpz_height(self, mpz_t height) except -1
    cdef int _rescale(self, mpq_t a) except -1

    cdef _pickle_version0(self)
    cdef _unpickle_version0(self, data)
    cpdef _export_as_string(self, int base=?)

    cdef _add_ui_unsafe_assuming_int(self, Py_ssize_t i, Py_ssize_t j, unsigned long int n)
    cdef _sub_ui_unsafe_assuming_int(self, Py_ssize_t i, Py_ssize_t j, unsigned long int n)

cdef class MatrixWindow:
    cdef Matrix_rational_dense _matrix
    cdef int _row, _col, _nrows, _ncols

################################################################
# fast conversion to pari on the stack
################################################################
ctypedef long* GEN
cdef inline GEN pari_GEN(Matrix_rational_dense B)
