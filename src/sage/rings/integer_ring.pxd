from .ring cimport IntegralDomain
from .integer cimport Integer
from sage.libs.gmp.types cimport mpz_t

cdef class IntegerRing_class(IntegralDomain):
    cdef int _randomize_mpz(self, mpz_t value, x, y, distribution) except -1
    cdef object _zero
