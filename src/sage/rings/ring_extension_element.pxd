from sage.structure.element cimport CommutativeAlgebraElement
from sage.structure.element cimport Element

cdef class RingExtensionElement(CommutativeAlgebraElement):
    cdef Element _element

cdef class RingExtensionWithBasisElement(RingExtensionElement):
    pass
