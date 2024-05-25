from .matroid cimport Matroid

cdef class FlatsMatroid(Matroid):
    cdef frozenset _groundset
    cdef int _matroid_rank
    cdef dict _F  # flats
    cdef object _L  # lattice of flats
    cpdef groundset(self)
    cpdef _rank(self, X)
    cpdef full_rank(self)
    cpdef _closure(self, X)
    cpdef _is_closed(self, X)

    # enumeration
    cpdef flats(self, k)
    cpdef whitney_numbers(self)
    cpdef whitney_numbers2(self)

    # isomorphism and relabeling
    cpdef _is_isomorphic(self, other, certificate=*)
    cpdef relabel(self, mapping)

    # verification
    cpdef is_valid(self)
