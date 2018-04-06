from sage.structure.element cimport CommutativeRingElement

cdef class FiniteRingElement(CommutativeRingElement):
    cpdef _add_(self, other)
    cpdef _mul_(self, other)

cdef class FinitePolyExtElement(FiniteRingElement):
    pass
