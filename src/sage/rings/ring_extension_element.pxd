# sage_setup: distribution = sagemath-modules
from sage.rings.ring cimport CommutativeRing
from sage.structure.element cimport Element
from sage.structure.element cimport CommutativeAlgebraElement
from sage.rings.ring_extension cimport RingExtension_generic
from sage.rings.ring_extension cimport RingExtensionFractionField
from sage.rings.ring_extension cimport RingExtensionWithBasis


cdef class RingExtensionElement(CommutativeAlgebraElement):
    cdef Element _backend

cdef class RingExtensionFractionFieldElement(RingExtensionElement):
    pass

cdef class RingExtensionWithBasisElement(RingExtensionElement):
    cdef _vector(self, CommutativeRing base)
    cdef _matrix(self, CommutativeRing base)
    cdef _trace(self, CommutativeRing base)
    cdef _norm(self, CommutativeRing base)
    cpdef minpoly(self, base=*, var=*)


