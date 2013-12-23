include 'decl.pxi'

cimport sage.structure.parent_base
cimport cython

from sage.libs.pari.gen cimport gen

cpdef long prec_bits_to_words(unsigned long prec_in_bits)

@cython.final
cdef class PariInstance(sage.structure.parent_base.ParentWithBase):
    cdef gen PARI_ZERO, PARI_ONE, PARI_TWO
    cdef inline gen new_gen(self, GEN x)
    cdef inline gen new_gen_noclear(self, GEN x)
    cdef gen new_gen_from_mpz_t(self, mpz_t value)
    cdef inline GEN _new_GEN_from_mpz_t(self, mpz_t value)
    cdef gen new_gen_from_mpq_t(self, mpq_t value)
    cdef inline GEN _new_GEN_from_mpq_t(self, mpq_t value)
    cdef gen new_gen_from_int(self, int value)
    cdef gen new_t_POL_from_int_star(self, int *vals, int length, long varnum)
    cdef gen new_gen_from_padic(self, long ordp, long relprec, mpz_t prime, mpz_t p_pow, mpz_t unit)
    cdef inline void clear_stack(self)
    cdef gen double_to_gen_c(self, double)
    cdef GEN double_to_GEN(self, double)
    cdef GEN deepcopy_to_python_heap(self, GEN x, pari_sp* address)
    cdef gen new_ref(self, GEN g, gen parent)
    cdef gen _empty_vector(self, long n)
    cdef long get_var(self, v)
    cdef GEN _new_GEN_from_mpz_t_matrix(self, mpz_t** B, Py_ssize_t nr, Py_ssize_t nc)
    cdef GEN _new_GEN_from_mpz_t_matrix_rotate90(self, mpz_t** B, Py_ssize_t nr, Py_ssize_t nc)
    cdef gen integer_matrix(self, mpz_t** B, Py_ssize_t nr, Py_ssize_t nc, bint permute_for_hnf)
    cdef GEN _new_GEN_from_mpq_t_matrix(self, mpq_t** B, Py_ssize_t nr, Py_ssize_t nc)
    cdef gen rational_matrix(self, mpq_t** B, Py_ssize_t nr, Py_ssize_t nc)
