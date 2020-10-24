from .matroid cimport Matroid

cdef class MatroidUnion(Matroid):
    cdef list matroids
    cdef frozenset _groundset
    cpdef groundset(self)
    cpdef _rank(self, X)

cdef class MatroidSum(Matroid):
    cdef list summands
    cdef frozenset _groundset
    cpdef groundset(self)
    cpdef _rank(self, X)

cdef class PartitionMatroid(Matroid):
    cdef dict p
    cdef frozenset _groundset
    cpdef groundset(self)
    cpdef _rank(self, X)
