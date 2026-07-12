from sage.libs.flint.arb_poly cimport *
from sage.rings.polynomial.polynomial_element cimport Polynomial

cdef class Polynomial_real_arb(Polynomial):
    cdef arb_poly_struct[1] _poly # https://github.com/cython/cython/issues/1984
    cdef Polynomial_real_arb _new(self)
