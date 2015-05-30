from sage.libs.mpfi cimport *

cimport sage.structure.element
cimport real_mpfi
cimport complex_number

cdef class ComplexIntervalFieldElement(sage.structure.element.FieldElement):
    cdef mpfi_t __re
    cdef mpfi_t __im
    cdef int _prec

    cdef real_mpfi.RealIntervalFieldElement abs_c(ComplexIntervalFieldElement self)
    cdef real_mpfi.RealIntervalFieldElement norm_c(ComplexIntervalFieldElement self)

    cdef ComplexIntervalFieldElement _new(self)
