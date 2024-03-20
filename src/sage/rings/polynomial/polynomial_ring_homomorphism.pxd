# sage_setup: distribution = sagemath-categories
from sage.rings.morphism cimport RingHomomorphism_from_base

cdef class PolynomialRingHomomorphism_from_base(RingHomomorphism_from_base):
    pass
