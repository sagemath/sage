# sage_setup: distribution = sagemath-pari
from sage.rings.morphism cimport RingHomomorphism
from sage.structure.element cimport Element


cdef class FrobeniusEndomorphism_padics(RingHomomorphism):
    cdef long _degree
    cdef long _power
    cdef long _order

    cpdef Element _call_(self, x)
