from sage.rings.ring cimport CommutativeRing
from sage.structure.element cimport Element
from sage.structure.element cimport CommutativeAlgebraElement

cdef class RingExtensionElement(CommutativeAlgebraElement):
    cdef public Element _backend

cdef class RingExtensionFractionFieldElement(RingExtensionElement):
    pass

cdef class RingExtensionWithBasisElement(RingExtensionElement):
    cdef _vector(self, CommutativeRing base)
    cdef _matrix(self, CommutativeRing base)
    cdef _trace(self, CommutativeRing base)
    cdef _norm(self, CommutativeRing base)
    cpdef minpoly(self, var=*, base=*)


