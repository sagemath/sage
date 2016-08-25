from .ring cimport PrincipalIdealDomain
from .integer cimport Integer
from sage.libs.gmp.types cimport mpz_t
from sage.libs.ntl.types cimport ZZ_c

cdef class IntegerRing_class(PrincipalIdealDomain):
    cdef Integer _coerce_ZZ(self, ZZ_c *z)
    cdef int _randomize_mpz(self, mpz_t value, x, y, distribution) except -1
    cdef object _zero
