from sage.libs.gmp.types cimport *
from sage.rings.finite_rings.stdint cimport *
from sage.rings.finite_rings.element_base cimport FiniteRingElement
from sage.rings.integer cimport Integer


cdef class NativeIntStruct:
    cdef Integer sageInteger
    cdef int_fast32_t int32
    cdef int_fast64_t int64
    cdef readonly list table     # list of elements of IntegerModRing(n)
    cdef readonly list inverses  # list of inverses (None if not invertible)
    cdef inline type element_class(self):
        if self.int32 > 0:
            return IntegerMod_int
        elif self.int64 > 0:
            return IntegerMod_int64
        else:
            return IntegerMod_gmp


cdef class IntegerMod_abstract(FiniteRingElement):
    cdef NativeIntStruct _modulus
    cdef _new_c_from_long(self, long value)
    cdef IntegerMod_abstract _new_c_fast(self, unsigned long value)
    cdef void set_from_mpz(self, mpz_t value) noexcept
    cdef void set_from_long(self, long value) noexcept
    cdef void set_from_ulong_fast(self, unsigned long value) noexcept
    cdef bint is_square_c(self) except -2
    cpdef bint is_one(self) noexcept
    cpdef bint is_unit(self) noexcept
    cpdef _floordiv_(self, other)

cdef class IntegerMod_gmp(IntegerMod_abstract):
    cdef mpz_t value
    cdef IntegerMod_gmp _new_c(self)
    cdef shift(IntegerMod_gmp self, long k)

cdef class IntegerMod_int(IntegerMod_abstract):
    cdef int_fast32_t ivalue
    cdef void set_from_int(IntegerMod_int self, int_fast32_t value) noexcept
    cdef int_fast32_t get_int_value(IntegerMod_int self) noexcept
    cdef IntegerMod_int _new_c(self, int_fast32_t value)
    cdef shift(IntegerMod_int self, int k)

cdef class IntegerMod_int64(IntegerMod_abstract):
    cdef int_fast64_t ivalue
    cdef void set_from_int(IntegerMod_int64 self, int_fast64_t value) noexcept
    cdef int_fast64_t get_int_value(IntegerMod_int64 self) noexcept
    cdef IntegerMod_int64 _new_c(self, int_fast64_t value)
    cdef shift(IntegerMod_int64 self, int k)


cdef int_fast32_t mod_inverse_int(int_fast32_t x, int_fast32_t n) except 0
cdef bint use_32bit_type(int_fast64_t modulus) noexcept
