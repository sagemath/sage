from .types cimport *
from sage.libs.gmp.types cimport *
from sage.libs.flint.types cimport fmpz_t, fmpz_mat_t
from sage.structure.parent_base cimport ParentWithBase
cimport cython

from sage.libs.pari.gen cimport gen

cpdef long prec_bits_to_words(unsigned long prec_in_bits)
cpdef long prec_words_to_bits(long prec_in_words)
cpdef long default_bitprec()

cdef class PariInstance_auto(ParentWithBase):
    pass

@cython.final
cdef class PariInstance(PariInstance_auto):
    cdef long _real_precision
    cdef readonly gen PARI_ZERO, PARI_ONE, PARI_TWO
    cpdef gen zero(self)
    cpdef gen one(self)
    cdef inline gen new_gen(self, GEN x)
    cdef inline gen new_gen_noclear(self, GEN x)
    cdef gen new_gen_from_mpz_t(self, mpz_t value)
    cdef inline GEN _new_GEN_from_mpz_t(self, mpz_t value)
    cdef gen new_gen_from_mpq_t(self, mpq_t value)
    cdef inline GEN _new_GEN_from_mpq_t(self, mpq_t value)
    cdef inline GEN _new_GEN_from_fmpz_t(self, fmpz_t value)
    cdef gen new_gen_from_int(self, int value)
    cdef gen new_t_POL_from_int_star(self, int *vals, int length, long varnum)
    cdef gen new_gen_from_padic(self, long ordp, long relprec, mpz_t prime, mpz_t p_pow, mpz_t unit)
    cdef inline void clear_stack(self)
    cdef gen double_to_gen_c(self, double)
    cdef GEN double_to_GEN(self, double)
    cdef GEN deepcopy_to_python_heap(self, GEN x, pari_sp* address)
    cdef gen new_ref(self, GEN g, gen parent)
    cdef gen _empty_vector(self, long n)
    cdef long get_var(self, v) except -2
    cdef GEN _new_GEN_from_fmpz_mat_t(self, fmpz_mat_t B, Py_ssize_t nr, Py_ssize_t nc)
    cdef GEN _new_GEN_from_fmpz_mat_t_rotate90(self, fmpz_mat_t B, Py_ssize_t nr, Py_ssize_t nc)
    cdef gen integer_matrix(self, fmpz_mat_t B, Py_ssize_t nr, Py_ssize_t nc, bint permute_for_hnf)
    cdef GEN _new_GEN_from_mpq_t_matrix(self, mpq_t** B, Py_ssize_t nr, Py_ssize_t nc)
    cdef gen rational_matrix(self, mpq_t** B, Py_ssize_t nr, Py_ssize_t nc)

cdef PariInstance pari_instance

cdef void INT_to_mpz(mpz_ptr value, GEN g)
cdef void INTFRAC_to_mpq(mpq_ptr value, GEN g)
