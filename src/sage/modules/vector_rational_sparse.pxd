# sage_setup: distribution = sagemath-modules
#############################################################
#
#    Sparse Vector over mpq_t (the GMP rationals)
#
#############################################################

from sage.libs.gmp.types cimport mpq_t

cdef struct mpq_vector:
    mpq_t *entries      # array of nonzero entries
    Py_ssize_t   *positions    # positions of those nonzero entries, starting at 0
    Py_ssize_t    degree       # the degree of this sparse vector
    Py_ssize_t    num_nonzero  # the number of nonzero entries of this vector.

cdef int reallocate_mpq_vector(mpq_vector* v, Py_ssize_t num_nonzero) except -1
cdef int allocate_mpq_vector(mpq_vector* v, Py_ssize_t num_nonzero) except -1
cdef int mpq_vector_init(mpq_vector* v, Py_ssize_t degree, Py_ssize_t num_nonzero) except -1
cdef void mpq_vector_clear(mpq_vector* v) noexcept
cdef Py_ssize_t mpq_binary_search0(mpq_t* v, Py_ssize_t n, mpq_t x) noexcept
cdef Py_ssize_t mpq_binary_search(mpq_t* v, Py_ssize_t n, mpq_t x, Py_ssize_t* ins) noexcept
cdef int mpq_vector_get_entry(mpq_t ans, mpq_vector* v, Py_ssize_t n) except -1
cdef bint mpq_vector_is_entry_zero_unsafe(mpq_vector* v, Py_ssize_t n) noexcept
cdef object mpq_vector_to_list(mpq_vector* v)
cdef int mpq_vector_set_entry(mpq_vector* v, Py_ssize_t n, mpq_t x) except -1
cdef int mpq_vector_set_entry_str(mpq_vector* v, Py_ssize_t n, char *x_str) except -1
cdef int add_mpq_vector_init(mpq_vector* sum, mpq_vector* v, mpq_vector* w, mpq_t multiple) except -1
cdef int mpq_vector_scale(mpq_vector* v, mpq_t scalar) except -1
cdef int mpq_vector_scalar_multiply(mpq_vector* v, mpq_vector* w, mpq_t scalar) except -1
cdef int mpq_vector_cmp(mpq_vector* v, mpq_vector* w) noexcept
