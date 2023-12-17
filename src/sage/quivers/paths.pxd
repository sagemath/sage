from sage.structure.element cimport MonoidElement, Element
from sage.data_structures.bounded_integer_sequences cimport biseq_t

cdef class QuiverPath(MonoidElement):
    cdef biseq_t _path
    cdef int _start, _end
    cdef QuiverPath _new_(self, int start, int end) noexcept
    cpdef _mul_(self, other) noexcept
    cpdef _mod_(self, right) noexcept
    cpdef tuple complement(self, QuiverPath subpath) noexcept
    cpdef bint has_subpath(self, QuiverPath subpath) except -1
    cpdef bint has_prefix(self, QuiverPath subpath) except -1
