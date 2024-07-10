from libc.stdint cimport uint32_t, uint64_t
from sage.data_structures.binary_matrix cimport *


cdef class IndependentSets:
    cdef binary_matrix_t g
    cdef list vertices
    cdef dict vertex_to_int
    cdef int n
    cdef int i
    cdef int count_only
    cdef int maximal
