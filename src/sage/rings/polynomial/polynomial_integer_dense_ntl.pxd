include "sage/libs/ntl/decl.pxi"    # to get ZZX_c etc

from sage.rings.polynomial.polynomial_element cimport Polynomial

cdef class Polynomial_integer_dense_ntl(Polynomial):
    cdef ZZX_c __poly

    cdef Polynomial_integer_dense_ntl _new(self)
