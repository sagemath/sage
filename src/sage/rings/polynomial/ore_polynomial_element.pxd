from sage.structure.element cimport AlgebraElement
from sage.structure.parent cimport Parent
from sage.rings.integer cimport Integer
from sage.rings.morphism cimport Morphism
from sage.structure.element cimport RingElement
from sage.rings.polynomial.polynomial_element cimport Polynomial_generic_dense

cdef class OrePolynomial(AlgebraElement):
    cdef _is_gen

    cdef long _hash_c(self) noexcept
    cdef OrePolynomial _new_c(self, list coeffs, Parent P, char check=*) noexcept
    cpdef OrePolynomial _new_constant_poly(self, RingElement a, Parent P, char check=*) noexcept
    cpdef _neg_(self) noexcept
    cpdef _floordiv_(self, right) noexcept
    cpdef _mod_(self, right) noexcept

    cpdef bint is_zero(self) noexcept
    cpdef bint is_one(self) noexcept
 
    cdef _left_quo_rem(self, OrePolynomial other) noexcept
    cdef _right_quo_rem(self, OrePolynomial other) noexcept
    cdef OrePolynomial _left_lcm_cofactor(self, OrePolynomial other) noexcept
    cdef OrePolynomial _right_lcm_cofactor(self, OrePolynomial other) noexcept

    # Abstract methods
    cpdef Integer degree(self) noexcept
    cpdef list coefficients(self, sparse=*) noexcept


cdef void lmul_gen(list A, Morphism m, d) noexcept

cdef class OrePolynomial_generic_dense(OrePolynomial):
    cdef list _coeffs

    cdef void _normalize(self) noexcept
    cpdef _add_(self, other) noexcept
    cdef list _mul_list(self, list A) noexcept
    cpdef _mul_(self, other) noexcept

    cpdef dict dict(self) noexcept
    cpdef list list(self, bint copy=*) noexcept


cdef class OrePolynomialBaseringInjection(Morphism):
    cdef RingElement _an_element
    cdef object _new_constant_poly_
