from sage.libs.arb.arb cimport arb_t
from sage.libs.mpfi.types cimport mpfi_t
from sage.rings.real_mpfi cimport RealIntervalField_class, RealIntervalFieldElement
from sage.structure.parent cimport Parent
from sage.structure.element cimport RingElement

cdef void mpfi_to_arb(arb_t target, const mpfi_t source, const long precision) noexcept
cdef int arb_to_mpfi(mpfi_t target, arb_t source, const long precision) except -1

cdef class RealBall(RingElement):
    cdef arb_t value
    cpdef _add_(self, other) noexcept
    cpdef _mul_(self, other) noexcept
    cpdef RealIntervalFieldElement _real_mpfi_(self, RealIntervalField_class parent) noexcept
    cpdef RealBall psi(self) noexcept

    cdef inline RealBall _new(self) noexcept:
        cdef RealBall res = RealBall.__new__(RealBall)
        res._parent = self._parent
        return res
