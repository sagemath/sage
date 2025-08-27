from sage.rings.ring cimport Field
from sage.structure.parent import Parent

cdef class RealField(Field):

    pass


cdef class RealIntervalField(Field):

    pass


cdef class RealDoubleField(Field):

    pass


cdef class ComplexField(Field):

    pass


cdef class ComplexDoubleField(Field):

    pass


cdef class SymbolicRing(Parent):

    pass
