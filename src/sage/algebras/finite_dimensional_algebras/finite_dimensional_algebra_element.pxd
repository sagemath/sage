from sage.structure.element cimport AlgebraElement, Element, Vector
from sage.matrix.matrix cimport Matrix

cdef class FiniteDimensionalAlgebraElement(AlgebraElement):
    cdef public Matrix _vector
    cdef Matrix __matrix
    cdef FiniteDimensionalAlgebraElement __inverse

    cpdef dict monomial_coefficients(self, bint copy=*)

cpdef FiniteDimensionalAlgebraElement unpickle_FiniteDimensionalAlgebraElement(A, vec, mat)
