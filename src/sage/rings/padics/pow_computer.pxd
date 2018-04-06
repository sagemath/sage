from sage.libs.gmp.types cimport mpz_t, mpz_srcptr
from sage.structure.sage_object cimport SageObject
from sage.rings.integer cimport Integer

cdef class PowComputer_class(SageObject):
    cdef Integer prime
    cdef Integer p2 # floor(p/2)
    cdef bint in_field
    cdef int __allocated
    cdef public object _prec_type

    # the following three should be set by the subclasses
    cdef long ram_prec_cap # = prec_cap * e
    cdef long deg
    cdef long e
    cdef long f

    cdef unsigned long cache_limit
    cdef unsigned long prec_cap

    cdef Integer pow_Integer(self, long n)
    cdef mpz_srcptr pow_mpz_t_top(self)
    cdef mpz_srcptr pow_mpz_t_tmp(self, long n) except NULL
    cdef mpz_t temp_m

cdef class PowComputer_base(PowComputer_class):
    cdef mpz_t* small_powers
    cdef mpz_t top_power
    cdef object __weakref__
