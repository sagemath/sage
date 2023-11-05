from sage.matroids.matroid cimport Matroid


cdef class CircuitClosuresMatroid(Matroid):
    cdef frozenset _groundset  # _E
    cdef dict _circuit_closures  # _CC
    cdef int _matroid_rank  # _R
    cpdef groundset(self) noexcept
    cpdef _rank(self, X) noexcept
    cpdef full_rank(self) noexcept
    cpdef _is_independent(self, F) noexcept
    cpdef _max_independent(self, F) noexcept
    cpdef _circuit(self, F) noexcept
    cpdef circuit_closures(self) noexcept
    cpdef _is_isomorphic(self, other, certificate=*) noexcept
