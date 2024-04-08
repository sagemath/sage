# sage_setup: distribution = sagemath-objects
cdef class lazy_list_generic():
    cdef list cache                  # the cache
    cdef lazy_list_generic master   # a reference if self is a slice
    cdef Py_ssize_t start, stop, step

    cpdef get(self, Py_ssize_t i)
    cpdef int _fit(self, Py_ssize_t n) except -1
    cpdef int _update_cache_up_to(self, Py_ssize_t i) except -1
    cpdef list _get_cache_(self)

cdef class lazy_list_from_iterator(lazy_list_generic):
    cdef object iterator

cdef class lazy_list_from_function(lazy_list_generic):
    cdef object callable

cdef class lazy_list_from_update_function(lazy_list_generic):
    cdef object update_function
