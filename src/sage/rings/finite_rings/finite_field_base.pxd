from sage.rings.ring cimport Field

cdef class FiniteField(Field):
    pass

cdef class FiniteFieldAbsolute(FiniteField):
    pass
