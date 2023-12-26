from sage.libs.arb.acb cimport acb_t
from sage.rings.complex_interval cimport ComplexIntervalFieldElement
from sage.rings.real_arb cimport RealBall
from sage.structure.element cimport RingElement
from sage.rings.ring cimport Field

cdef void ComplexIntervalFieldElement_to_acb(
    acb_t target,
    ComplexIntervalFieldElement source) noexcept

cdef int acb_to_ComplexIntervalFieldElement(
    ComplexIntervalFieldElement target,
    const acb_t source) except -1

cdef class ComplexBall(RingElement):
    cdef acb_t value
    cdef ComplexBall _new(self) noexcept
    cpdef _add_(self, other) noexcept
    cpdef _mul_(self, other) noexcept
    cpdef ComplexIntervalFieldElement _complex_mpfi_(self, parent) noexcept
    cpdef RealBall real(self) noexcept
    cpdef RealBall imag(self) noexcept
    cpdef pow(self, expo, analytic=?) noexcept

    cdef inline ComplexBall _new(self) noexcept:
        cdef ComplexBall res = ComplexBall.__new__(ComplexBall)
        res._parent = self._parent
        return res
