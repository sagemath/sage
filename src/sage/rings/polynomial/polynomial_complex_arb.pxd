from sage.libs.arb.acb_poly cimport *
from sage.rings.polynomial.polynomial_element cimport Polynomial

cdef class Polynomial_complex_arb(Polynomial):
    cdef acb_poly_t __poly
    cdef Polynomial_complex_arb _new(self)
