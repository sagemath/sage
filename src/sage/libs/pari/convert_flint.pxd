from cypari2.types cimport GEN
from cypari2.gen cimport Gen
from sage.libs.flint.types cimport fmpz_t, fmpz_mat_t, fmpq_t, fmpq_mat_t

cdef GEN _new_GEN_from_fmpz_t(fmpz_t value) noexcept
cdef GEN _new_GEN_from_fmpz_mat_t(fmpz_mat_t B) noexcept
cdef GEN _new_GEN_from_fmpz_mat_t_rotate90(fmpz_mat_t B) noexcept
cdef Gen integer_matrix(fmpz_mat_t B, bint rotate)

cdef GEN _new_GEN_from_fmpq_t(fmpq_t value) noexcept
cdef GEN _new_GEN_from_fmpq_mat_t(fmpq_mat_t B) noexcept
cdef GEN _new_GEN_from_fmpq_mat_t_rotate90(fmpq_mat_t B) noexcept
cdef Gen rational_matrix(fmpq_mat_t B, bint rotate)
