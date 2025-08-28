from sage.data_structures.bitset cimport bitset_t
from sage.data_structures.binary_matrix cimport binary_matrix_t

cdef class ConvexityProperties:
    cdef int _n
    cdef list _list_integers_to_vertices
    cdef dict _dict_vertices_to_integers
    cdef binary_matrix_t _cache_hull_pairs

    cdef list _vertices_to_integers(self, vertices)
    cdef list _integers_to_vertices(self, list integers)
    cdef _bitset_convex_hull(self, bitset_t hull)
    cpdef hull(self, list vertices)
    cdef _greedy_increase(self, bitset_t bs)
    cpdef hull_number(self, value_only=*, verbose=*)
