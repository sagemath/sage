from sage.matroids.matroid cimport Matroid
from sage.matroids.set_system cimport SetSystem

cdef class CircuitsMatroid(Matroid):
    cdef frozenset _groundset  # _E
    cdef int _matroid_rank  # _R
    cdef SetSystem _C  # circuits
    cdef dict _k_C  # k-circuits (k=len)
    cdef bint _nsc_defined
    cpdef groundset(self) noexcept
    cpdef _rank(self, X) noexcept
    cpdef full_rank(self) noexcept
    cpdef _is_independent(self, F) noexcept
    cpdef _max_independent(self, F) noexcept
    cpdef _circuit(self, F) noexcept

    # enumeration
    cpdef bases(self) noexcept
    cpdef circuits(self, k=*) noexcept
    cpdef nonspanning_circuits(self) noexcept
    cpdef no_broken_circuits_sets(self, ordering=*) noexcept

    # properties
    cpdef girth(self) noexcept
    cpdef is_paving(self) noexcept

    # isomorphism
    cpdef _is_isomorphic(self, other, certificate=*) noexcept

    # verification
    cpdef is_valid(self) noexcept
