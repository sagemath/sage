# sage_setup: distribution = sagemath-flint

from sage.libs.flint.acb cimport acb_t
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
    cdef ComplexBall _new(self)
    cpdef _add_(self, other)
    cpdef _mul_(self, other)
    cpdef ComplexIntervalFieldElement _complex_mpfi_(self, parent)
    cpdef RealBall real(self)
    cpdef RealBall imag(self)
    cpdef pow(self, expo, analytic=?)

    cdef inline ComplexBall _new(self):
        cdef ComplexBall res = ComplexBall.__new__(ComplexBall)
        res._parent = self._parent
        return res
