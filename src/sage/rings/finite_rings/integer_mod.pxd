include "sage/ext/cdefs.pxi"

from sage.rings.finite_rings.stdint cimport *
from sage.rings.finite_rings.element_base cimport FiniteRingElement
from sage.rings.integer cimport Integer

cdef class NativeIntStruct:
    cdef Integer sageInteger
    cdef int_fast32_t int32
    cdef int_fast64_t int64
    cdef object table # a list
    cdef object inverses # also a list
    cdef lookup(NativeIntStruct self, Py_ssize_t value)

cdef class IntegerMod_abstract(FiniteRingElement):
    cdef NativeIntStruct __modulus
    cdef _new_c_from_long(self, long value)
    cdef void set_from_mpz(self, mpz_t value)
    cdef void set_from_long(self, long value)
    cdef bint is_square_c(self) except -2
    cpdef bint is_one(self)
    cpdef bint is_unit(self)

cdef class IntegerMod_gmp(IntegerMod_abstract):
    cdef mpz_t value
    cdef IntegerMod_gmp _new_c(self)
    cdef shift(IntegerMod_gmp self, long k)

cdef class IntegerMod_int(IntegerMod_abstract):
    cdef int_fast32_t ivalue
    cdef void set_from_int(IntegerMod_int self, int_fast32_t value)
    cdef int_fast32_t get_int_value(IntegerMod_int self)
    cdef IntegerMod_int _new_c(self, int_fast32_t value)
    cdef shift(IntegerMod_int self, int k)
    #cdef Element _make_new_with_parent_c(self, ParentWithBase parent)

cdef class IntegerMod_int64(IntegerMod_abstract):
    cdef int_fast64_t ivalue
    cdef void set_from_int(IntegerMod_int64 self, int_fast64_t value)
    cdef int_fast64_t get_int_value(IntegerMod_int64 self)
    cdef IntegerMod_int64 _new_c(self, int_fast64_t value)
    cdef shift(IntegerMod_int64 self, int k)

cdef int jacobi_int(int_fast32_t a, int_fast32_t m) except -2
cdef int_fast32_t mod_inverse_int(int_fast32_t x, int_fast32_t n) except 0
cdef int_fast32_t mod_pow_int(int_fast32_t base, int_fast32_t exp, int_fast32_t n)


cdef bint use_32bit_type(int_fast64_t modulus)
