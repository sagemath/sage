# sage_setup: distribution = sagemath-singular
from sage.libs.singular.decl cimport poly, ring

from sage.rings.polynomial.multi_polynomial cimport MPolynomial_libsingular as MPolynomial_libsingular_base
from sage.rings.polynomial.multi_polynomial_ring_base cimport MPolynomialRing_base

cdef class MPolynomialRing_libsingular

cdef class MPolynomial_libsingular(MPolynomial_libsingular_base):
    cdef poly *_poly
    cdef ring *_parent_ring
    cpdef _add_(self, other)
    cpdef _mul_(self, other)
    cpdef _floordiv_(self, right)
    cpdef _repr_short_(self)
    cpdef is_constant(self)
    cpdef _homogenize(self, int var)
    cpdef MPolynomial_libsingular _new_constant_poly(self, x, MPolynomialRing_libsingular P)
    cpdef long number_of_terms(self) noexcept

cdef class MPolynomialRing_libsingular(MPolynomialRing_base):
    cdef object __singular
    cdef object __macaulay2
    cdef object __m2_set_ring_cache
    cdef object __minpoly
    cdef poly *_one_element_poly
    cdef ring *_ring

# new polynomials
cdef MPolynomial_libsingular new_MP(MPolynomialRing_libsingular parent, poly *p)
