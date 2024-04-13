from sage.rings.ring cimport CommutativeRing, Field

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


cdef class SymbolicRing(CommutativeRing):

    pass
