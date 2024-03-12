from sage.libs.gmp.types cimport mpz_t, mpz_srcptr
from sage.structure.sage_object cimport SageObject
from sage.rings.integer cimport Integer

cdef class PowComputer_class(SageObject):
    cdef Integer prime
    cdef Integer p2 # floor(p/2)
    cdef bint in_field
    cdef int _allocated
    cdef public object _prec_type

    cdef long ram_prec_cap # = prec_cap * e

    # the following constants should be set by the subclasses
    # the absolute degree of the p-adic ring
    cdef long deg
    # the absolute ramification index of the p-adic ring
    cdef long e
    # the absolute residual degree of the p-adic ring
    cdef long f

    cdef unsigned long cache_limit
    cdef unsigned long prec_cap

    cdef Integer pow_Integer(self, long n)
    cdef mpz_srcptr pow_mpz_t_top(self) noexcept
    cdef mpz_srcptr pow_mpz_t_tmp(self, long n) except NULL
    cdef mpz_t temp_m

cdef class PowComputer_base(PowComputer_class):
    cdef mpz_t* small_powers
    cdef mpz_t top_power
    cdef mpz_t powhelper_oneunit
    cdef mpz_t powhelper_teichdiff
    cdef mpz_t shift_rem
    cdef mpz_t aliasing
    cdef object __weakref__
