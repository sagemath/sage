# sage_setup: distribution = sagemath-modules
from sage.rings.real_double cimport RealDoubleElement


cdef class RealDoubleElement_gsl(RealDoubleElement):
    cdef __pow_double(self, double exponent, double sign)
    cpdef _pow_(self, other)
    cdef _log_base(self, double log_of_base)
