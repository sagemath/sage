from sage.libs.flint.types cimport nmod_mpoly_t, nmod_mpoly_ctx_t

from sage.rings.polynomial.multi_polynomial cimport MPolynomial_flint as MPolynomial_flint_base
from sage.rings.polynomial.multi_polynomial_ring_base cimport MPolynomialRing_base


cdef class MPolynomialRing_zmod_flint(MPolynomialRing_base):
    cdef nmod_mpoly_ctx_t _ctx
    cdef unsigned long _modulus
    cdef MPolynomial_zmod_flint _new_element(self)


cdef class MPolynomial_zmod_flint(MPolynomial_flint_base):
    cdef nmod_mpoly_t _poly
    cdef MPolynomial_zmod_flint _new(self)
    cpdef _new_constant_poly(self, scalar, parent)
    cpdef long number_of_terms(self) noexcept
