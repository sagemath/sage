from libc.stdint cimport uint8_t

cdef class FastDigraph:
    cdef uint8_t n
    cdef int * graph
    cdef list int_to_vertices
    cdef int * degree

cdef int compute_out_neighborhood_cardinality(FastDigraph, int) noexcept

cdef int popcount32(int) noexcept
cdef int slow_popcount32(int) noexcept
