from sage.rings.finite_rings.hom_finite_field cimport (SectionFiniteFieldHomomorphism_generic,
    FiniteFieldHomomorphism_generic, FrobeniusEndomorphism_finite_field)


cdef class SectionFiniteFieldHomomorphism_prime(SectionFiniteFieldHomomorphism_generic):
    pass


cdef class FiniteFieldHomomorphism_prime(FiniteFieldHomomorphism_generic):
    pass


cdef class FrobeniusEndomorphism_prime(FrobeniusEndomorphism_finite_field):
    pass
