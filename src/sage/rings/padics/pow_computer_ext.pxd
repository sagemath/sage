from sage.rings.padics.pow_computer cimport PowComputer_class
from sage.libs.ntl.ntl_ZZ_pContext cimport ntl_ZZ_pContext_class
from sage.libs.ntl.types cimport *
from sage.libs.gmp.types cimport mpz_t

cdef class PowComputer_ext(PowComputer_class):
    cdef ZZ_c* small_powers
    cdef ZZ_c top_power
    cdef ZZ_c temp_z
    cdef long _initialized
    cdef mpz_t temp_m2

    # the following are for unpickling
    cdef object _poly
    cdef object _shift_seed
    cdef object _ext_type

    cdef ZZ_c* pow_ZZ_tmp(self, long n) except NULL
    cdef ZZ_c* pow_ZZ_top(self) noexcept

    cdef void cleanup_ext(self) noexcept

cdef class PowComputer_ZZ_pX(PowComputer_ext):
    cdef ntl_ZZ_pContext_class get_context(self, long n)
    cdef ntl_ZZ_pContext_class get_context_capdiv(self, long n)
    cdef ntl_ZZ_pContext_class get_top_context(self)
    cdef restore_context(self, long n)
    cdef restore_context_capdiv(self, long n)
    cdef void restore_top_context(self) noexcept
    cdef ZZ_pX_Modulus_c* get_modulus(self, long n) noexcept
    cdef ZZ_pX_Modulus_c* get_modulus_capdiv(self, long n) noexcept
    cdef ZZ_pX_Modulus_c* get_top_modulus(self) noexcept
    cdef int eis_shift(self, ZZ_pX_c* x, ZZ_pX_c* a, long n, long finalprec) except -1
    cdef int eis_shift_capdiv(self, ZZ_pX_c* x, ZZ_pX_c* a, long n, long finalprec) except -1
    cdef long capdiv(self, long n) noexcept
    cdef int teichmuller_set_c (self, ZZ_pX_c* x, ZZ_pX_c* a, long absprec) except -1

cdef class PowComputer_ZZ_pX_FM(PowComputer_ZZ_pX):
    cdef ntl_ZZ_pContext_class c
    cdef ZZ_pX_Modulus_c mod

cdef class PowComputer_ZZ_pX_FM_Eis(PowComputer_ZZ_pX_FM):
    cdef int low_length
    cdef int high_length
    cdef ZZ_pX_Multiplier_c* low_shifter
    cdef ZZ_pX_Multiplier_c* high_shifter

    cdef void cleanup_ZZ_pX_FM_Eis(self) noexcept

cdef class PowComputer_ZZ_pX_small(PowComputer_ZZ_pX):
    cdef object c  # using a python list so that we can store ntl_ZZ_pContext_class objects
    cdef ZZ_pX_Modulus_c *mod

    cdef void cleanup_ZZ_pX_small(self) noexcept

cdef class PowComputer_ZZ_pX_small_Eis(PowComputer_ZZ_pX_small):
    cdef int low_length
    cdef int high_length
    cdef ZZ_pX_c* low_shifter
    cdef ZZ_pX_c* high_shifter

    cdef void cleanup_ZZ_pX_small_Eis(self) noexcept

cdef class PowComputer_ZZ_pX_big(PowComputer_ZZ_pX):
    cdef object context_list  # using a python list so that we can store ntl_ZZ_pContext_class objects
    cdef ZZ_pX_Modulus_c *modulus_list

    cdef ntl_ZZ_pContext_class top_context
    cdef ZZ_pX_Modulus_c top_mod

    cdef object context_dict  # currently using a dict, optimize for speed later
    cdef object modulus_dict  # currently using a dict, optimize for speed later

    cdef void cleanup_ZZ_pX_big(self) noexcept

cdef class PowComputer_ZZ_pX_big_Eis(PowComputer_ZZ_pX_big):
    cdef int low_length
    cdef int high_length
    cdef ZZ_pX_c* low_shifter
    cdef ZZ_pX_c* high_shifter

    cdef void cleanup_ZZ_pX_big_Eis(self) noexcept

# cdef class PowComputer_ZZ_pEX(PowComputer_ext):
#    cdef ntl_ZZ_pEContext get_context(self, long n)
#    cdef ntl_ZZ_pEContext get_top_context(self)
#    cdef void restore_context(self, long n)
#    cdef void restore_top_context(self)
#    cdef ZZ_pEX_Modulus_c get_modulus(self, long n)
#    cdef ZZ_pEX_Modulus_c get_top_modulus(self)
#
# cdef class PowComputer_ZZ_pEX_FM(PowComputer_ZZ_pEX):
#    cdef ZZ_pEX_c poly
#    cdef ntl_ZZ_pEContext c
#    cdef ZZ_pEX_Modulus_c mod
#
#
# cdef class PowComputer_ZZ_pEX_small(PowComputer_ZZ_pEX):
#    cdef ZZ_pEX_c *poly
#    cdef ntl_ZZ_pEContext *c
#    cdef ZZ_pEX_Modulus_c *mod
#
# cdef class PowComputer_ZZ_pEX_big(PowComputer_ZZ_pEX):
#    cdef ZZ_pEX_c poly
#    cdef context_dict  # currently using a dict, optimize for speed later
#    cdef modulus_dict  # currently using a dict, optimize for speed later
