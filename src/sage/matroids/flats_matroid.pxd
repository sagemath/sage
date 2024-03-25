from sage.matroids.matroid cimport Matroid

cdef class FlatsMatroid(Matroid):
    cdef frozenset _groundset  # _E
    cdef int _matroid_rank  # _R
    cdef dict _F  # flats
    cpdef groundset(self)
    cpdef _rank(self, X)
    cpdef full_rank(self)
    cpdef _is_independent(self, F)

    # enumeration
    cpdef flats(self, k)
    cpdef whitney_numbers2(self)

    # isomorphism
    cpdef _is_isomorphic(self, other, certificate=*)

    # verification
    cpdef is_valid(self)
