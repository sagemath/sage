# sage_setup: distribution = sagemath-categories
cimport cython

cdef struct pair_s:
    size_t first
    size_t second


@cython.final
cdef class ListOfPairs:
    cdef pair_s** _lists
    cdef size_t length

    cdef inline int enlarge(self) except -1
    cdef inline int add(self, size_t first, size_t second) except -1
    cdef inline pair_s* get(self, size_t index) except NULL
