from sage.libs.gmp.types cimport mpz_t, mpz_ptr
from sage.libs.ntl.types cimport ZZ_c

from sage.structure.element cimport EuclideanDomainElement, RingElement
from sage.categories.morphism cimport Morphism

cdef class Integer(EuclideanDomainElement):
    cdef mpz_t value

    cdef int _to_ZZ(self, ZZ_c *z) except -1
    cdef void set_from_mpz(self, mpz_t value)
    cdef hash_c(self)

    cpdef _pari_(self)

    cpdef _shift_helper(Integer self, y, int sign)
    cdef _and(Integer self, Integer other)
    cdef _or(Integer self, Integer other)
    cdef _xor(Integer self, Integer other)

    cpdef size_t _exact_log_log2_iter(self,Integer m)
    cpdef size_t _exact_log_mpfi_log(self,m)
    cpdef RingElement _valuation(Integer self, Integer p)
    cdef object _val_unit(Integer self, Integer p)
    cdef Integer _divide_knowing_divisible_by(Integer self, Integer right)
    cdef bint _is_power_of(Integer self, Integer n)

    cdef bint _pseudoprime_is_prime(self, proof) except -1
    cpdef list _pari_divisors_small(self)

    cdef _reduce_set(self, s) # do not use, since integers are immutable.

cdef int mpz_set_str_python(mpz_ptr z, char* s, int base) except -1

cdef Integer smallInteger(long value)

cdef class int_to_Z(Morphism):
    pass
