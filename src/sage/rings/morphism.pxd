from sage.rings.integer cimport Integer
from sage.structure.element cimport Element
from sage.categories.morphism cimport Morphism
from sage.rings.ring cimport Ring
from sage.structure.category_object cimport CategoryObject

cdef class RingMap(Morphism):
    pass

cdef class RingMap_lift(RingMap):
    #cdef Ring S
    # S. King, trac ticket #11068: It is not used that it is a ring.
    # Thus, "CategoryObject" should be fine, and allows to
    # do computations in quotients of matrix algebras.
    cdef CategoryObject S

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
