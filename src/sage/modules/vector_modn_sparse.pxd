cdef struct c_vector_modint:
    int *entries
    int p
    Py_ssize_t *positions
    Py_ssize_t degree
    Py_ssize_t num_nonzero

cdef int allocate_c_vector_modint(c_vector_modint* v, Py_ssize_t num_nonzero) except -1
cdef int init_c_vector_modint(c_vector_modint* v, int p, Py_ssize_t degree, Py_ssize_t num_nonzero) except -1
cdef void clear_c_vector_modint(c_vector_modint* v)
cdef Py_ssize_t binary_search0_modn(Py_ssize_t* v, Py_ssize_t n, int x)
cdef Py_ssize_t binary_search_modn(Py_ssize_t* v, Py_ssize_t n, int x, Py_ssize_t* ins)
cdef int get_entry(c_vector_modint* v, Py_ssize_t n) except -1
cdef object to_list(c_vector_modint* v)
cdef int set_entry(c_vector_modint* v, Py_ssize_t n, int x) except -1
cdef int add_c_vector_modint_init(c_vector_modint* sum, c_vector_modint* v, c_vector_modint* w, int multiple) except -1
cdef int scale_c_vector_modint(c_vector_modint* v, int scalar) except -1