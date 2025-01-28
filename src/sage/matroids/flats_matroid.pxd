from sage.matroids.matroid cimport Matroid
from sage.matroids.set_system cimport SetSystem

cdef class FlatsMatroid(Matroid):
    cdef frozenset _groundset
    cdef int _matroid_rank
    cdef dict _F  # flats
    cdef object _L  # lattice of flats
    cpdef frozenset groundset(self)

    cpdef int _rank(self, frozenset X) except? -1
    cpdef frozenset _closure(self, frozenset X)
    cpdef bint _is_closed(self, frozenset X) noexcept

    cpdef full_rank(self)

    # enumeration
    cpdef SetSystem flats(self, long k)
    cpdef list whitney_numbers(self)
    cpdef list whitney_numbers2(self)

    # isomorphism and relabeling
    cpdef _is_isomorphic(self, other, certificate=*)
    cpdef relabel(self, mapping)

    # verification
    cpdef is_valid(self, certificate=*)
