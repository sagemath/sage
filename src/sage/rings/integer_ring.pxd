from sage.rings.ring cimport CommutativeRing
from sage.libs.gmp.types cimport mpz_t

cdef class IntegerRing_class(CommutativeRing):
    cdef int _randomize_mpz(self, mpz_t value, x, y, distribution) except -1
    cdef object _zero
