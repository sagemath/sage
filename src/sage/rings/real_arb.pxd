cimport sage.rings.number_field.number_field_element_quadratic as nfeq

from sage.libs.arb.arb cimport arb_t
from sage.libs.mpfi cimport mpfi_t
from sage.rings.real_mpfi cimport RealIntervalField_class, RealIntervalFieldElement
from sage.structure.parent cimport Parent
from sage.structure.element cimport RingElement

cdef void mpfi_to_arb(arb_t target, const mpfi_t source, const long precision)
cdef int arb_to_mpfi(mpfi_t target, arb_t source, const long precision) except -1
cdef int real_part_of_quadratic_element_to_arb(arb_t res, nfeq.NumberFieldElement_quadratic x, const long prec) except -1

cdef class RealBall(RingElement):
    cdef arb_t value
    cdef RealBall _new(self)
    cpdef _add_(self, other)
    cpdef _mul_(self, other)
    cpdef RealIntervalFieldElement _real_mpfi_(self, RealIntervalField_class parent)
    cpdef RealBall psi(self)
