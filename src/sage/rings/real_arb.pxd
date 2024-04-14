# sage_setup: distribution = sagemath-flint

from sage.libs.flint.arb cimport arb_t
from sage.libs.mpfi.types cimport mpfi_t
from sage.rings.real_mpfi cimport RealIntervalField_class, RealIntervalFieldElement
from sage.structure.parent cimport Parent
from sage.structure.element cimport RingElement

cdef void mpfi_to_arb(arb_t target, const mpfi_t source, const long precision) noexcept
cdef int arb_to_mpfi(mpfi_t target, arb_t source, const long precision) except -1

cdef class RealBall(RingElement):
    cdef arb_t value
    cpdef _add_(self, other)
    cpdef _mul_(self, other)
    cpdef RealIntervalFieldElement _real_mpfi_(self, RealIntervalField_class parent)
    cpdef RealBall psi(self)

    cdef inline RealBall _new(self):
        cdef RealBall res = RealBall.__new__(RealBall)
        res._parent = self._parent
        return res
