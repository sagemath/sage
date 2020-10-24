from sage.structure.element cimport AlgebraElement
from sage.structure.parent cimport Parent
from sage.rings.morphism cimport Morphism
from sage.structure.element cimport RingElement
from sage.rings.polynomial.polynomial_element cimport Polynomial_generic_dense

cdef class OrePolynomial(AlgebraElement):
    cdef _is_gen

    cdef long _hash_c(self)
    cdef OrePolynomial _new_c(self, list coeffs, Parent P, char check=*)
    cpdef OrePolynomial _new_constant_poly(self, RingElement a, Parent P, char check=*)
    cpdef _neg_(self)
    cpdef _floordiv_(self, right)
    cpdef _mod_(self, right)

    cpdef bint is_zero(self)
    cpdef bint is_one(self)
 
    cdef _left_quo_rem(self, OrePolynomial other)
    cdef _right_quo_rem(self, OrePolynomial other)
    cdef OrePolynomial _left_lcm_cofactor(self, OrePolynomial other)
    cdef OrePolynomial _right_lcm_cofactor(self, OrePolynomial other)

    # Abstract methods
    cpdef int degree(self)
    cpdef list coefficients(self, sparse=*)


cdef void lmul_gen(list A, Morphism m, d)

cdef class OrePolynomial_generic_dense(OrePolynomial):
    cdef list _coeffs

    cdef void __normalize(self)
    cpdef _add_(self, other)
    cdef list _mul_list(self, list A)
    cpdef _mul_(self, other)

    cpdef dict dict(self)
    cpdef list list(self, bint copy=*)


cdef class OrePolynomialBaseringInjection(Morphism):
    cdef RingElement _an_element
    cdef object _new_constant_poly_
