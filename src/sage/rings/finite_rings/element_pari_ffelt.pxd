from cypari2.types cimport GEN
from sage.rings.finite_rings.element_base cimport FinitePolyExtElement


cdef class FiniteFieldElement_pari_ffelt(FinitePolyExtElement):
    # PARI t_FFELT describing the element.
    # This holds a reference to a PARI clone.
    cdef GEN val

    cdef FiniteFieldElement_pari_ffelt _new(self)
    cdef void construct(self, GEN g) noexcept
    cdef int construct_from(self, x) except -1
