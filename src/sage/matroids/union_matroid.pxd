from sage.matroids.matroid cimport Matroid

cdef class MatroidUnion(Matroid):
    cdef list matroids
    cdef frozenset _groundset
    cpdef frozenset groundset(self)
    cpdef int _rank(self, frozenset X) except? -1

cdef class MatroidSum(Matroid):
    cdef list summands
    cdef frozenset _groundset
    cpdef frozenset groundset(self)
    cpdef int _rank(self, frozenset X) except? -1

cdef class PartitionMatroid(Matroid):
    cdef dict p
    cdef frozenset _groundset
    cpdef frozenset groundset(self)
    cpdef int _rank(self, frozenset X) except? -1
