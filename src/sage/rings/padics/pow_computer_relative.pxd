from sage.rings.padics.pow_computer cimport PowComputer_class
from sage.rings.polynomial.polynomial_element cimport Polynomial_generic_dense

cdef class PowComputer_relative(PowComputer_class):
    # p-adic elements are represented by polynomials in this ring
    cdef public object poly_ring
    # the p-adic ring is an extension of this base ring
    cdef public object base_ring
    # the modulus of the extension
    cdef public Polynomial_generic_dense modulus
    # storage for temporary variables used in the linkage files
    cdef Polynomial_generic_dense tmp_cconv_out
    cdef Polynomial_generic_dense tmp_ccoeffs
    cdef Polynomial_generic_dense tmp_ccoeffs_frac
    cdef Polynomial_generic_dense tmp_ccmp_a
    cdef Polynomial_generic_dense tmp_ccmp_b
    cdef Polynomial_generic_dense shift_rem
    cdef Polynomial_generic_dense aliasing
    # allow cached methods
    cdef public dict _cached_methods

    cdef unsigned long capdiv(self, unsigned long n) noexcept

cdef class PowComputer_relative_eis(PowComputer_relative):
    # (x^e - modulus)/p
    cdef public Polynomial_generic_dense _shift_seed
    cdef public Polynomial_generic_dense _inv_shift_seed
    cpdef Polynomial_generic_dense invert(self, Polynomial_generic_dense element, long prec)
