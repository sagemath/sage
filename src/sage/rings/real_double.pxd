# sage_setup: distribution = sagemath-categories
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

cdef double_repr(double x)
