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
    cdef public dict _cached_methods

cdef class RingHomomorphism_im_gens(RingHomomorphism):
    cdef _im_gens
    cdef _base_map

cdef class RingHomomorphism_from_base(RingHomomorphism):
    cdef _underlying

cdef class RingHomomorphism_from_fraction_field(RingHomomorphism):
    cdef _morphism

cdef class RingHomomorphism_cover(RingHomomorphism):
    pass

cdef class RingHomomorphism_from_quotient(RingHomomorphism):
    cdef phi

cdef class FrobeniusEndomorphism_generic(RingHomomorphism):
    cdef Integer _p
    cdef Integer _q
    cdef long _power
