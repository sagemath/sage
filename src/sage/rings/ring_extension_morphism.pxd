from sage.structure.element cimport Element
from sage.rings.morphism cimport RingHomomorphism

cdef class RingExtensionHomomorphism(RingHomomorphism):
    cdef _backend_morphism
    cdef _gens
    cdef _im_gens
    cdef _base_map
    cpdef Element _call_(self, x)
