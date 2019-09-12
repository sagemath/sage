from sage.structure.element cimport Element
from sage.rings.morphism cimport RingHomomorphism

cdef class RingExtensionHomomorphism(RingHomomorphism):
    cdef _backend_morphism
    cpdef Element _call_(self, x)
