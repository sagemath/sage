# sage_setup: distribution = sagemath-categories
from sage.rings.morphism cimport RingHomomorphism_im_gens, FrobeniusEndomorphism_generic
from sage.structure.element cimport Element
from sage.categories.map cimport Section


cdef class SectionFiniteFieldHomomorphism_generic(Section):
    pass


cdef class FiniteFieldHomomorphism_generic(RingHomomorphism_im_gens):
    cdef _gen
    cdef _section_class

    cpdef Element _call_(self, x)


cdef class FrobeniusEndomorphism_finite_field(FrobeniusEndomorphism_generic):
    cdef long _degree
    cdef long _degree_fixed
    cdef long _order

    cpdef Element _call_(self, x)
