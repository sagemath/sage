from sage.libs.gmp.types cimport *
from sage.libs.flint.types cimport nmod_mat_t
from sage.ext.mod_int cimport mod_int

from .matrix_dense cimport Matrix_dense
from  sage.rings.finite_rings.integer_mod cimport NativeIntStruct

cdef class Matrix_modn_dense_flint(Matrix_dense):
    cdef nmod_mat_t _matrix
    cdef NativeIntStruct _modulus

    cdef Matrix_modn_dense_flint _new(self, Py_ssize_t nrows, Py_ssize_t ncols)
    cpdef _shift_mod(self, mp_limb_t modulus, mp_limb_t shift=*, bint mul=*, bint domod=*)
    cdef int _copy_row_to_mod_int_array(self, mod_int *to, Py_ssize_t i)
