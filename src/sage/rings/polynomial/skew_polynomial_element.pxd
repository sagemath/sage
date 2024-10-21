# sage_setup: distribution = sagemath-modules
from sage.rings.polynomial.ore_polynomial_element cimport OrePolynomial_generic_dense

cdef class SkewPolynomial_generic_dense(OrePolynomial_generic_dense):
    cdef void _inplace_rmul(self, SkewPolynomial_generic_dense right) noexcept
    cdef void _inplace_pow(self, Py_ssize_t n) noexcept
    cpdef right_power_mod(self, exp, modulus)
    cpdef left_power_mod(self, exp, modulus)
    cpdef operator_eval(self, eval_pt)
