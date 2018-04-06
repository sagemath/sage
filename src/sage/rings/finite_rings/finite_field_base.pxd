from sage.rings.ring cimport Field

cdef class FiniteField(Field):
    cdef public object __factored_unit_order
    cdef public object __multiplicative_generator
    cdef public object __polynomial_ring
    cdef public object __vector_space
    cdef public object __interface

