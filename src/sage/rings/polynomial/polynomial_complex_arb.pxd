from sage.libs.arb.acb_poly cimport *
from sage.rings.polynomial.polynomial_element cimport Polynomial

cdef class Polynomial_complex_arb(Polynomial):
    cdef acb_poly_struct[1] __poly # https://github.com/cython/cython/issues/1984
    cdef Polynomial_complex_arb _new(self)
