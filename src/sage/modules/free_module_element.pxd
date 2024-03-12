from sage.structure.element cimport Vector

cdef class FreeModuleElement(Vector):
    cdef int set_unsafe(self, Py_ssize_t i, value) except -1
    cdef get_unsafe(self, Py_ssize_t i)
    cpdef int hamming_weight(self) noexcept

cdef class FreeModuleElement_generic_dense(FreeModuleElement):
    # data
    cdef list _entries

    # cdef'd methods
    cdef _new_c(self, object v)


cdef class FreeModuleElement_generic_sparse(FreeModuleElement):
    # data
    cdef dict _entries

    # cdef'd methods
    cdef _new_c(self, object v)

