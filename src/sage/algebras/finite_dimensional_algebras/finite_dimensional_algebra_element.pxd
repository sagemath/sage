from sage.structure.element cimport AlgebraElement, Element, Vector
from sage.matrix.matrix cimport Matrix

cdef class FiniteDimensionalAlgebraElement(AlgebraElement):
    cdef public Matrix _vector
    cdef Matrix __matrix
    cdef FiniteDimensionalAlgebraElement __inverse
    cdef int _hash

cpdef FiniteDimensionalAlgebraElement unpickle_FiniteDimensionalAlgebraElement(A, vec, mat)
