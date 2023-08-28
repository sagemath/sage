# sage_setup: distribution = sagemath-categories
from sage.rings.ring cimport Field

cdef class FiniteField(Field):
    cdef public object __polynomial_ring
    cdef public object __vector_space
    cdef public object __interface
