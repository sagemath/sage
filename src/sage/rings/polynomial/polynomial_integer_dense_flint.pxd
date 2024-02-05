from sage.libs.flint.types cimport fmpz_poly_t

from sage.rings.polynomial.polynomial_element cimport Polynomial
from sage.rings.integer cimport Integer
from sage.structure.parent cimport Parent

cdef class Polynomial_integer_dense_flint(Polynomial):
    cdef fmpz_poly_t _poly

    cdef Polynomial_integer_dense_flint _new(self) noexcept
    cpdef _unsafe_mutate(self, long n, value) noexcept
    cpdef Integer content(self) noexcept
