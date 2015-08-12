from sage.libs.pari.types cimport GEN
from sage.rings.finite_rings.element_base cimport FinitePolyExtElement

cdef class FiniteFieldElement_pari_ffelt(FinitePolyExtElement):
    cdef GEN val        # PARI t_FFELT describing the element
    cdef void *block    # memory block containing the data
    cdef FiniteFieldElement_pari_ffelt _new(FiniteFieldElement_pari_ffelt self)
    cdef void construct(FiniteFieldElement_pari_ffelt self, GEN g)
    cdef void construct_from(FiniteFieldElement_pari_ffelt self, object x) except *
