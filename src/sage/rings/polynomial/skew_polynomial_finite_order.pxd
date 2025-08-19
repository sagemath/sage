from sage.rings.polynomial.skew_polynomial_element cimport SkewPolynomial_generic_dense

cdef class SkewPolynomial_finite_order_dense (SkewPolynomial_generic_dense):
    cdef _norm
    cdef _charpoly
    cdef _optbound

    cdef _matphir_c(self)
    cdef _matmul_c(self)
