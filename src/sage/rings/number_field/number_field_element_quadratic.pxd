from sage.libs.gmp.types cimport mpz_t
from sage.libs.arb.types cimport arb_t
from sage.rings.integer cimport Integer
from .number_field_element cimport NumberFieldElement, NumberFieldElement_absolute


cdef class NumberFieldElement_quadratic(NumberFieldElement_absolute):
    # (a + b sqrt(D)) / denom
    cdef mpz_t a, b, denom
    cdef Integer D
    cdef bint standard_embedding
    cpdef NumberFieldElement galois_conjugate(self)
    cdef bint is_sqrt_disc(self)

    cpdef list _coefficients(self)

    cdef int _randomize(self, num_bound, den_bound, distribution) except -1
    cdef int arb_set_real(self, arb_t x, long prec) except -1
    cdef void arb_set_imag(self, arb_t x, long prec)

cdef class OrderElement_quadratic(NumberFieldElement_quadratic):
    pass
