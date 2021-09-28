from sage.structure.element cimport RingElement, ModuleElement, Element, FieldElement
from sage.rings.ring cimport Field
cimport sage.rings.abc

cdef class RealDoubleField_class(sage.rings.abc.RealDoubleField):
    cdef _new_c(self, double value)

cdef class RealDoubleElement(FieldElement):
    cdef double _value
    cdef _new_c(self, double value)
    cpdef _add_(self, other)
    cpdef _mul_(self, other)
    cpdef RealDoubleElement abs(RealDoubleElement self)
    cdef __pow_double(self, double exponent, double sign)
    cpdef _pow_(self, other)
    cdef _log_base(self, double log_of_base)

cdef double_repr(double x)
