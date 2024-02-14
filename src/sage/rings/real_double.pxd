from sage.structure.element cimport RingElement, ModuleElement, Element, FieldElement
from sage.rings.ring cimport Field
cimport sage.rings.abc

cdef class RealDoubleField_class(sage.rings.abc.RealDoubleField):
    cdef _new_c(self, double value) noexcept

cdef class RealDoubleElement(FieldElement):
    cdef double _value
    cdef _new_c(self, double value) noexcept
    cpdef _add_(self, other) noexcept
    cpdef _mul_(self, other) noexcept
    cpdef RealDoubleElement abs(RealDoubleElement self) noexcept

cdef double_repr(double x) noexcept
