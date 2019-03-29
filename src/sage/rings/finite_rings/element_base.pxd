from sage.structure.element cimport CommutativeRingElement

cdef class FiniteRingElement(CommutativeRingElement):
    pass

cdef class FinitePolyExtElement(FiniteRingElement):
    pass
