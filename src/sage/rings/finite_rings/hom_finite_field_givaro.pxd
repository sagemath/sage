# sage_setup: distribution = sagemath-linbox

from sage.rings.finite_rings.hom_finite_field cimport (SectionFiniteFieldHomomorphism_generic,
    FiniteFieldHomomorphism_generic, FrobeniusEndomorphism_finite_field)

from sage.structure.element cimport Element
from sage.rings.finite_rings.element_givaro cimport Cache_givaro


cdef class SectionFiniteFieldHomomorphism_givaro(SectionFiniteFieldHomomorphism_generic):
    cdef long _order_codomain
    cdef long _gcd
    cdef long _power
    cdef Cache_givaro _codomain_cache

    cpdef Element _call_(self, x)


cdef class FiniteFieldHomomorphism_givaro(FiniteFieldHomomorphism_generic):
    cdef long _order_domain
    cdef long _order_codomain
    cdef long _power
    cdef Cache_givaro _codomain_cache

    cpdef Element _call_(self, x)


cdef class FrobeniusEndomorphism_givaro(FrobeniusEndomorphism_finite_field):
    pass
