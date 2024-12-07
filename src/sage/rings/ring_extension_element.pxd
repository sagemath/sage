from sage.structure.parent cimport Parent
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
    cdef _vector(self, Parent base)
    cdef _matrix(self, Parent base)
    cdef _trace(self, Parent base)
    cdef _norm(self, Parent base)
    cpdef minpoly(self, base=*, var=*)
