from sage.matroids.matroid cimport Matroid

cdef class CircuitClosuresMatroid(Matroid):
    cdef frozenset _groundset  # _E
    cdef dict _circuit_closures  # _CC
    cdef int _matroid_rank  # _R
    cpdef frozenset groundset(self)
    cpdef int _rank(self, frozenset X) except? -1
    cpdef full_rank(self)
    cpdef bint _is_independent(self, frozenset F) noexcept
    cpdef frozenset _max_independent(self, frozenset F)
    cpdef frozenset _circuit(self, frozenset F)
    cpdef dict circuit_closures(self)
    cpdef _is_isomorphic(self, other, certificate=*)
    cpdef relabel(self, mapping)
