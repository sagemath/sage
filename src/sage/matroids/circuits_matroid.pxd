from sage.matroids.matroid cimport Matroid
from sage.matroids.set_system cimport SetSystem

cdef class CircuitsMatroid(Matroid):
    cdef frozenset _groundset  # _E
    cdef int _matroid_rank  # _R
    cdef SetSystem _C  # circuits
    cdef dict _k_C  # k-circuits (k=len)
    cdef bint _nsc_defined
    cpdef groundset(self)
    cpdef _rank(self, X)
    cpdef full_rank(self)
    cpdef _is_independent(self, F)
    cpdef _max_independent(self, F)
    cpdef _circuit(self, F)

    # enumeration
    cpdef bases(self)
    cpdef circuits(self, k=*)
    cpdef nonspanning_circuits(self)
    cpdef no_broken_circuits_sets(self, ordering=*)

    # properties
    cpdef girth(self)
    cpdef is_paving(self)

    # isomorphism and relabeling
    cpdef _is_isomorphic(self, other, certificate=*)
    cpdef relabel(self, mapping)

    # verification
    cpdef is_valid(self)
