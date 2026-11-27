from sage.libs.flint.types cimport fmpz_mpoly_t, fmpz_mpoly_ctx_t

from sage.rings.polynomial.multi_polynomial cimport MPolynomial_flint as MPolynomial_flint_base
from sage.rings.polynomial.multi_polynomial_ring_base cimport MPolynomialRing_base


cdef class MPolynomialRing_integer_flint(MPolynomialRing_base):
    cdef fmpz_mpoly_ctx_t _ctx
    cdef MPolynomial_integer_flint _new_element(self)


cdef class MPolynomial_integer_flint(MPolynomial_flint_base):
    cdef fmpz_mpoly_t _poly
    cdef MPolynomial_integer_flint _new(self)
    cpdef _new_constant_poly(self, scalar, parent)
    cpdef long number_of_terms(self) noexcept
