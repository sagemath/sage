# sage_setup: distribution = sagemath-flint
from sage.libs.gmp.types cimport *
from sage.libs.flint.types cimport fmpz_mat_t

from sage.matrix.matrix_dense cimport Matrix_dense
from sage.rings.integer cimport Integer
from sage.ext.mod_int cimport *

ctypedef long* GEN

cdef class Matrix_integer_dense(Matrix_dense):
    cdef fmpz_mat_t _matrix
    cdef object _pivots
    cdef int mpz_height(self, mpz_t height) except -1
    cdef _mod_int_c(self, mod_int modulus)
    cdef _mod_two(self)
    cdef _pickle_version0(self)
    cdef _unpickle_version0(self, data)
    cpdef _export_as_string(self, int base=?)
    cdef void set_unsafe_mpz(self, Py_ssize_t i, Py_ssize_t j, const mpz_t value) noexcept
    cdef void set_unsafe_double(self, Py_ssize_t i, Py_ssize_t j, double value) noexcept
    cdef inline void get_unsafe_mpz(self, Py_ssize_t i, Py_ssize_t j, mpz_t value) noexcept
    cdef inline int get_unsafe_int(self, Py_ssize_t i, Py_ssize_t j) noexcept
    cdef inline double get_unsafe_double(self, Py_ssize_t i, Py_ssize_t j) noexcept

    # HNF Modn
    cdef int _hnf_modn(Matrix_integer_dense self, Matrix_integer_dense res,
            unsigned int det) except -1
    cdef int* _hnf_modn_impl(Matrix_integer_dense self, unsigned int det,
            Py_ssize_t nrows, Py_ssize_t ncols) except NULL

    cdef Matrix_integer_dense _new(self, Py_ssize_t nrows, Py_ssize_t ncols)


cpdef _lift_crt(Matrix_integer_dense M, residues, moduli=*)
