from sage.libs.mpfr cimport *

cimport sage.structure.element
from .real_mpfr cimport RealNumber

cdef class ComplexNumber(sage.structure.element.FieldElement):
    cdef mpfr_t __re
    cdef mpfr_t __im
    cdef object _multiplicative_order
    cdef int _prec

    cdef RealNumber abs_c(ComplexNumber self)
    cdef RealNumber norm_c(ComplexNumber self)

    cdef ComplexNumber _new(self)
