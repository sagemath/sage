from sage.data_structures.bitset cimport bitset_t
from sage.matroids.basis_matroid cimport BasisMatroid

cdef class CutNode:
    cdef LinearSubclasses _MC
    cdef bitset_t _p_free, _p_in, _l0, _l1
    cdef long _ml

    cdef CutNode copy(self)
    cdef bint insert_plane(self, long p0) noexcept
    cdef bint remove_plane(self, long p0) noexcept
    cdef select_plane(self)

    cdef list planes(self)

cdef class LinearSubclassesIter:
    cdef LinearSubclasses _MC
    cdef list _nodes

cdef class LinearSubclasses:
    cdef object _hyperlines, _hyperplanes
    cdef long _hyperlines_count, _hyperplanes_count
    cdef list _planes_on_line, _lines_on_plane

    cdef list _mandatory_planes
    cdef list _mandatory_lines
    cdef list _forbidden_planes
    cdef _line_length

cdef class MatroidExtensions(LinearSubclasses):
    cdef BasisMatroid _BX
    cdef list _BH
    cdef _orderly
