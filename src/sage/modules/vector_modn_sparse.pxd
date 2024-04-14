# sage_setup: distribution = sagemath-modules
from sage.rings.finite_rings.stdint cimport *

cdef struct c_vector_modint:
    int_fast64_t *entries
    int p
    Py_ssize_t *positions
    Py_ssize_t degree
    Py_ssize_t num_nonzero

cdef int allocate_c_vector_modint(c_vector_modint* v, Py_ssize_t num_nonzero) except -1
cdef int init_c_vector_modint(c_vector_modint* v, int p, Py_ssize_t degree, Py_ssize_t num_nonzero) except -1
cdef void clear_c_vector_modint(c_vector_modint* v) noexcept
cdef Py_ssize_t binary_search0_modn(Py_ssize_t* v, Py_ssize_t n, int_fast64_t x) noexcept
cdef Py_ssize_t binary_search_modn(Py_ssize_t* v, Py_ssize_t n, int_fast64_t x, Py_ssize_t* ins) noexcept
cdef int_fast64_t get_entry(c_vector_modint* v, Py_ssize_t n) except -1
cdef bint is_entry_zero_unsafe(c_vector_modint* v, Py_ssize_t n) noexcept
cdef object to_list(c_vector_modint* v)
cdef int set_entry(c_vector_modint* v, Py_ssize_t n, int_fast64_t x) except -1
cdef int add_c_vector_modint_init(c_vector_modint* sum, c_vector_modint* v, c_vector_modint* w, int multiple) except -1
cdef int scale_c_vector_modint(c_vector_modint* v, int_fast64_t scalar) except -1
