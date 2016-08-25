from sage.libs.gmp.types cimport *
from sage.libs.flint.fmpz cimport *
from sage.libs.flint.fmpz_mat cimport *

from .matrix_dense cimport Matrix_dense
from sage.rings.integer cimport Integer
from sage.ext.mod_int cimport *

ctypedef long* GEN

cdef class Matrix_integer_dense(Matrix_dense):
    cdef fmpz_mat_t _matrix  # Always initialized in __cinit__
    cdef bint _initialized_mpz
    cdef mpz_t * _entries    # Only used if _initialized_mpz
    cdef mpz_t ** _rows      # Only used if _initialized_mpz
    cdef object _pivots
    cdef int mpz_height(self, mpz_t height) except -1
    cdef _mod_int_c(self, mod_int modulus)
    cdef _mod_two(self)
    cdef _pickle_version0(self)
    cdef _unpickle_version0(self, data)
    cpdef _export_as_string(self, int base=?)
    cdef inline int _init_mpz(self) except -1
    cdef int _init_mpz_impl(self) except -1
    cdef inline int _init_linbox(self) except -1
    cdef void _dealloc_mpz(self)
    cdef void set_unsafe_mpz(self, Py_ssize_t i, Py_ssize_t j, const mpz_t value)
    cdef void set_unsafe_si(self, Py_ssize_t i, Py_ssize_t j, long value)
    cdef void set_unsafe_double(self, Py_ssize_t i, Py_ssize_t j, double value)
    cdef inline void get_unsafe_mpz(self, Py_ssize_t i, Py_ssize_t j, mpz_t value)
    cdef inline double get_unsafe_double(self, Py_ssize_t i, Py_ssize_t j)

    # HNF Modn
    cdef int _hnf_modn(Matrix_integer_dense self, Matrix_integer_dense res,
            unsigned int det) except -1
    cdef int* _hnf_modn_impl(Matrix_integer_dense self, unsigned int det,
            Py_ssize_t nrows, Py_ssize_t ncols) except NULL

    cdef Matrix_integer_dense _new(self, Py_ssize_t nrows, Py_ssize_t ncols)

    cdef extract_hnf_from_pari_matrix(self, GEN H, int flag, bint include_zero_rows)

cpdef _lift_crt(Matrix_integer_dense M, residues, moduli=*)

################################################################
# fast conversion to PARI on the stack
################################################################
cdef inline GEN pari_GEN(Matrix_integer_dense B)
