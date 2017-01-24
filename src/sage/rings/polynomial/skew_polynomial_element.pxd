from sage.structure.element cimport AlgebraElement
from sage.structure.parent cimport Parent
from sage.rings.morphism cimport Morphism
from sage.structure.element cimport RingElement
from sage.rings.polynomial.polynomial_element cimport Polynomial_generic_dense

cdef class SkewPolynomial(AlgebraElement):
    cdef _is_gen

    cdef long _hash_c(self)
    cdef SkewPolynomial _new_c(self,list coeffs,Parent P,char check=*)
    cpdef SkewPolynomial _new_constant_poly(self,RingElement a,Parent P,char check=*)
    cpdef _neg_(self)

    cpdef bint is_zero(self)
    cpdef bint is_one(self)

    cpdef operator_eval(self, eval_pt)

    # Abstract methods
    cdef void _inplace_rmul(self, SkewPolynomial_generic_dense right)
    cdef void _inplace_pow(self, Py_ssize_t n)
    cpdef int degree(self)
    cpdef list coefficients(self, sparse=*)

cdef class SkewPolynomial_generic_dense(SkewPolynomial):
    cdef list _coeffs

    cdef void __normalize(self)

    cpdef dict dict(self)
    cpdef list list(self)

    cpdef right_power_mod(self, exp, modulus)
    cpdef left_power_mod(self, exp, modulus)

cdef class SkewPolynomialBaseringInjection(Morphism):
    cdef RingElement _an_element
    cdef object _new_constant_poly_

