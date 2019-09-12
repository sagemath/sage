from sage.structure.element cimport CommutativeAlgebraElement
from sage.structure.element cimport Element

cdef class AlgebraFMElement(CommutativeAlgebraElement):
    cdef Element _element

cdef class RingExtensionWithBasisElement(AlgebraFMElement):
    pass
