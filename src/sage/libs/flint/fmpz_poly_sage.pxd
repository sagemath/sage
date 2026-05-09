# distutils: libraries = flint
# distutils: depends = flint/fmpz_poly.h

from sage.libs.gmp.types cimport mpz_t
from .types cimport *

# functions removed from flint but still needed in sage
cdef void fmpz_poly_scalar_mul_mpz(fmpz_poly_t, const fmpz_poly_t, const mpz_t) noexcept
cdef void fmpz_poly_scalar_divexact_mpz(fmpz_poly_t, const fmpz_poly_t, const mpz_t) noexcept
cdef void fmpz_poly_scalar_fdiv_mpz(fmpz_poly_t, const fmpz_poly_t, const mpz_t) noexcept
cdef void fmpz_poly_set_coeff_mpz(fmpz_poly_t, slong, const mpz_t) noexcept
cdef void fmpz_poly_get_coeff_mpz(mpz_t, const fmpz_poly_t, slong) noexcept
cdef void fmpz_poly_set_mpz(fmpz_poly_t, const mpz_t) noexcept


# Wrapper Cython class
from sage.structure.sage_object cimport SageObject
cdef class Fmpz_poly(SageObject):
    cdef fmpz_poly_t poly
