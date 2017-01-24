from sage.libs.cypari2.types cimport GEN
from sage.libs.cypari2.gen cimport gen
from sage.libs.flint.types cimport fmpz_t, fmpz_mat_t

cdef GEN _new_GEN_from_fmpz_t(fmpz_t value)
cdef GEN _new_GEN_from_fmpz_mat_t(fmpz_mat_t B, Py_ssize_t nr, Py_ssize_t nc)
cdef GEN _new_GEN_from_fmpz_mat_t_rotate90(fmpz_mat_t B, Py_ssize_t nr, Py_ssize_t nc)
cdef gen integer_matrix(fmpz_mat_t B, Py_ssize_t nr, Py_ssize_t nc, bint permute_for_hnf)
