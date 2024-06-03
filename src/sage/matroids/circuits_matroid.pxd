from sage.matroids.matroid cimport Matroid
from sage.matroids.set_system cimport SetSystem

cdef class CircuitsMatroid(Matroid):
    cdef frozenset _groundset
    cdef int _matroid_rank
    cdef set _C  # circuits
    cdef dict _k_C  # k-circuits (k=len)
    cdef list _sorted_C_lens
    cdef bint _nsc_defined
    cpdef frozenset groundset(self)
    cpdef int _rank(self, frozenset X)
    cpdef full_rank(self)
    cpdef bint _is_independent(self, frozenset X)
    cpdef frozenset _max_independent(self, frozenset X)
    cpdef frozenset _circuit(self, frozenset X)
    cpdef frozenset _closure(self, frozenset X)

    # enumeration
    cpdef SetSystem independent_sets(self, long k=*)
    cpdef SetSystem dependent_sets(self, long k)
    cpdef SetSystem circuits(self, k=*)
    cpdef SetSystem nonspanning_circuits(self)
    cpdef SetSystem no_broken_circuits_facets(self, ordering=*, reduced=*)
    cpdef SetSystem no_broken_circuits_sets(self, ordering=*, reduced=*)
    cpdef broken_circuit_complex(self, ordering=*, reduced=*)

    # properties
    cpdef girth(self)
    cpdef bint is_paving(self)

    # isomorphism and relabeling
    cpdef _is_isomorphic(self, other, certificate=*)
    cpdef relabel(self, mapping)

    # verification
    cpdef bint is_valid(self)
