from sage.libs.mpfr cimport *

cimport sage.structure.element
cimport real_mpfr

cdef class ComplexNumber(sage.structure.element.FieldElement):
    cdef mpfr_t __re
    cdef mpfr_t __im
    cdef object _multiplicative_order
    cdef int _prec

    cdef real_mpfr.RealNumber abs_c(ComplexNumber self)
    cdef real_mpfr.RealNumber norm_c(ComplexNumber self)

    cdef ComplexNumber _new(self)
