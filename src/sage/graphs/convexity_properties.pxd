from sage.misc.bitset cimport bitset_t

cdef class ConvexityProperties:
    cdef int _n
    cdef list _list_integers_to_vertices
    cdef dict _dict_vertices_to_integers
    cdef bitset_t * _cache_hull_pairs

    cdef list _vertices_to_integers(self, vertices)
    cdef list _integers_to_vertices(self, integers)
    cdef _bitset_convex_hull(self, bitset_t hull)
    cpdef hull(self, list vertices)
    cdef _greedy_increase(self, bitset_t bs)
    cpdef hull_number(self, value_only = *, verbose = *)
