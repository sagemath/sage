from sage.libs.mpfi cimport *

cimport sage.structure.element
from .real_mpfi cimport RealIntervalFieldElement

cdef class ComplexIntervalFieldElement(sage.structure.element.FieldElement):
    cdef mpfi_t __re
    cdef mpfi_t __im
    cdef int _prec

    cdef RealIntervalFieldElement abs_c(ComplexIntervalFieldElement self)
    cdef RealIntervalFieldElement norm_c(ComplexIntervalFieldElement self)

    cdef ComplexIntervalFieldElement _new(self)
