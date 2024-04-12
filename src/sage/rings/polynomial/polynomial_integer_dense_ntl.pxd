# sage_setup: distribution = sagemath-ntl
from sage.libs.ntl.types cimport ZZX_c
from sage.rings.polynomial.polynomial_element cimport Polynomial

cdef class Polynomial_integer_dense_ntl(Polynomial):
    cdef ZZX_c _poly

    cdef Polynomial_integer_dense_ntl _new(self)
