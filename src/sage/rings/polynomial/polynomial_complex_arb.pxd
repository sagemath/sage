# sage_setup: distribution = sagemath-flint

from sage.libs.flint.acb_poly cimport *
from sage.rings.polynomial.polynomial_element cimport Polynomial

cdef class Polynomial_complex_arb(Polynomial):
    cdef acb_poly_struct[1] _poly # https://github.com/cython/cython/issues/1984
    cdef Polynomial_complex_arb _new(self)
