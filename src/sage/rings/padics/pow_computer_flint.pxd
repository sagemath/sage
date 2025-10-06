from sage.libs.gmp.mpz cimport mpz_t
from sage.libs.flint.types cimport padic_ctx_t, fmpz_t, fmpz_poly_t
from sage.rings.padics.pow_computer cimport PowComputer_class

cdef class PowComputer_flint(PowComputer_class):
    cdef padic_ctx_t ctx
    cdef fmpz_t fprime
    cdef fmpz_t half_prime
    cdef fmpz_t _fpow_array
    cdef fmpz_t _fpow_variable
    cdef mpz_t top_power

    cdef fmpz_t* pow_fmpz_t_tmp(self, unsigned long n) except NULL
    cdef unsigned long capdiv(self, unsigned long n) noexcept

    cdef fmpz_t tfmpz

cdef class PowComputer_flint_1step(PowComputer_flint):
    cdef fmpz_t q
    cdef fmpz_poly_t modulus
    cdef fmpz_poly_t powhelper_oneunit
    cdef fmpz_poly_t powhelper_teichdiff
    cdef fmpz_poly_t* _moduli
    cdef fmpz_poly_t* get_modulus(self, unsigned long n) noexcept
    cdef fmpz_poly_t* get_modulus_capdiv(self, unsigned long n) noexcept
    cdef _new_fmpz_poly(self, fmpz_poly_t value, var=*)

cdef class PowComputer_flint_unram(PowComputer_flint_1step):
    # WARNING:
    # These variables are modified by the linkage and must no be used anywhere else
    # (other than in __(c)init__)
    cdef fmpz_t fmpz_ccmp, fmpz_cval, fmpz_cinv, fmpz_cinv2, fmpz_cexp, fmpz_ctm, fmpz_cconv
    cdef fmpz_poly_t poly_cconv, poly_ctm, poly_ccmp, poly_cinv, poly_cisunit, poly_cinv2, poly_flint_rep, poly_matmod, shift_rem, aliasing
    cdef mpz_t mpz_cpow, mpz_ctm, mpz_cconv, mpz_matmod

cdef class PowComputer_flint_eis(PowComputer_flint_1step):
    pass
