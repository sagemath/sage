from .matroid cimport Matroid

cdef class FlatsMatroid(Matroid):
    cdef frozenset _groundset  # _E
    cdef int _matroid_rank  # _R
    cdef dict _F  # flats
    cpdef groundset(self)
    cpdef _rank(self, X)
    cpdef full_rank(self)
    cpdef _closure(self, X)
    cpdef _is_closed(self, X)

    # enumeration
    cpdef flats(self, k)
    cpdef whitney_numbers2(self)

    # isomorphism and relabeling
    cpdef _is_isomorphic(self, other, certificate=*)
    cpdef relabel(self, mapping)

    # verification
    cpdef is_valid(self)

cdef class LatticeOfFlatsMatroid(FlatsMatroid):
    cdef object _L  # lattice_of_flats
    cpdef whitney_numbers(self)
    cpdef is_valid(self)
