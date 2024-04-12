# sage_setup: distribution = sagemath-modules
#############################################################
#
#    Sparse Vector over mpz_t (the GMP integers)
#
#############################################################

from sage.libs.gmp.types cimport mpz_t

cdef struct mpz_vector:
    mpz_t *entries      # array of nonzero entries
    Py_ssize_t   *positions    # positions of those nonzero entries, starting at 0
    Py_ssize_t    degree       # the degree of this sparse vector
    Py_ssize_t    num_nonzero  # the number of nonzero entries of this vector.

cdef int allocate_mpz_vector(mpz_vector* v, Py_ssize_t num_nonzero) except -1
cdef int mpz_vector_init(mpz_vector* v, Py_ssize_t degree, Py_ssize_t num_nonzero) except -1
cdef void mpz_vector_clear(mpz_vector* v) noexcept
cdef Py_ssize_t mpz_binary_search0(mpz_t* v, Py_ssize_t n, mpz_t x) noexcept
cdef Py_ssize_t mpz_binary_search(mpz_t* v, Py_ssize_t n, mpz_t x, Py_ssize_t* ins) noexcept
cdef int mpz_vector_get_entry(mpz_t ans, mpz_vector* v, Py_ssize_t n) except -1
cdef bint mpz_vector_is_entry_zero_unsafe(mpz_vector* v, Py_ssize_t n) noexcept
cdef object mpz_vector_to_list(mpz_vector* v)
cdef int mpz_vector_set_entry(mpz_vector* v, Py_ssize_t n, mpz_t x) except -1
cdef int mpz_vector_set_entry_str(mpz_vector* v, Py_ssize_t n, char *x_str) except -1
cdef int add_mpz_vector_init(mpz_vector* sum, mpz_vector* v, mpz_vector* w, mpz_t multiple) except -1
cdef int mpz_vector_scale(mpz_vector* v, mpz_t scalar) except -1
cdef int mpz_vector_scalar_multiply(mpz_vector* v, mpz_vector* w, mpz_t scalar) except -1
cdef int mpz_vector_cmp(mpz_vector* v, mpz_vector* w) noexcept
