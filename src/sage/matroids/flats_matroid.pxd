from sage.matroids.matroid cimport Matroid

cdef class FlatsMatroid(Matroid):
    cdef frozenset _groundset  # _E
    cdef int _matroid_rank  # _R
    cdef dict _F  # flats
    cpdef groundset(self) noexcept
    cpdef _rank(self, X) noexcept
    cpdef full_rank(self) noexcept
    cpdef _is_independent(self, F) noexcept

    # enumeration
    cpdef flats(self, k) noexcept
    cpdef whitney_numbers2(self) noexcept

    # isomorphism
    cpdef _is_isomorphic(self, other, certificate=*) noexcept

    # verification
    cpdef is_valid(self) noexcept
