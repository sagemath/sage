from .matroid cimport Matroid

cdef class CircuitClosuresMatroid(Matroid):
    cdef frozenset _groundset  # _E
    cdef dict _circuit_closures  # _CC
    cdef int _matroid_rank  # _R
    cpdef groundset(self)
    cpdef _rank(self, X)
    cpdef full_rank(self)
    cpdef _is_independent(self, F)
    cpdef _max_independent(self, F)
    cpdef _circuit(self, F)
    cpdef circuit_closures(self)
    cpdef _is_isomorphic(self, other, certificate=*)
