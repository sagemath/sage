from sage.libs.singular.decl cimport poly, ring

from sage.rings.polynomial.multi_polynomial cimport MPolynomial
from sage.rings.polynomial.multi_polynomial_ring_generic cimport MPolynomialRing_generic

cdef class MPolynomialRing_libsingular

cdef class MPolynomial_libsingular(MPolynomial):
    cdef poly *_poly
    cdef ring *_parent_ring
    cpdef _repr_short_(self)
    cpdef is_constant(self)
    cpdef _homogenize(self, int var)
    cpdef MPolynomial_libsingular _new_constant_poly(self, x, MPolynomialRing_libsingular P)

cdef class MPolynomialRing_libsingular(MPolynomialRing_generic):
    cdef object __singular
    cdef object __macaulay2
    cdef object __m2_set_ring_cache
    cdef object __minpoly
    cdef poly *_one_element_poly
    cdef ring *_ring

# new polynomials
cdef MPolynomial_libsingular new_MP(MPolynomialRing_libsingular parent, poly *p)
