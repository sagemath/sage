from sage.rings.integer cimport Integer
from sage.structure.element cimport Element
from sage.structure.parent cimport Parent
from sage.categories.map cimport Map
from sage.categories.morphism cimport Morphism


cdef class RingMap(Morphism):
    pass

cdef class RingMap_lift(RingMap):
    cdef Parent S
    cdef Map to_S

cdef class RingHomomorphism(RingMap):
    cdef Morphism _lift

cdef class RingHomomorphism_im_gens(RingHomomorphism):
    cdef __im_gens

cdef class RingHomomorphism_from_base(RingHomomorphism):
    cdef __underlying

cdef class RingHomomorphism_cover(RingHomomorphism):
    pass

cdef class RingHomomorphism_from_quotient(RingHomomorphism):
    cdef phi

cdef class FrobeniusEndomorphism_generic(RingHomomorphism):
    cdef Integer _p
    cdef Integer _q
    cdef long _power
