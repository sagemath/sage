from sage.matroids.matroid cimport Matroid


cdef class MatroidUnion(Matroid):
    cdef list matroids
    cdef frozenset _groundset
    cpdef groundset(self) noexcept
    cpdef _rank(self, X) noexcept

cdef class MatroidSum(Matroid):
    cdef list summands
    cdef frozenset _groundset
    cpdef groundset(self) noexcept
    cpdef _rank(self, X) noexcept

cdef class PartitionMatroid(Matroid):
    cdef dict p
    cdef frozenset _groundset
    cpdef groundset(self) noexcept
    cpdef _rank(self, X) noexcept
