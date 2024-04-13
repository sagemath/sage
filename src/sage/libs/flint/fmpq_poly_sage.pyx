# sage_setup: distribution = sagemath-flint
# Functions removed from flint but still needed in Sage. Code adapted from
# earlier versions of flint.

from sage.libs.gmp.mpq cimport *
from sage.libs.flint.fmpz cimport *
from sage.libs.flint.fmpq cimport *
from sage.libs.flint.fmpq_poly cimport *


cdef void fmpq_poly_scalar_mul_mpz(fmpq_poly_t rop, const fmpq_poly_t op, const mpz_t c) noexcept:
    cdef fmpz_t f
    fmpz_init_set_readonly(f, c)
    fmpq_poly_scalar_mul_fmpz(rop, op, f)
    fmpz_clear_readonly(f)

cdef void fmpq_poly_scalar_mul_mpq(fmpq_poly_t rop, const fmpq_poly_t op, const mpq_t c) noexcept:
    cdef fmpq_t f
    fmpq_init_set_readonly(f, c)
    fmpq_poly_scalar_mul_fmpq(rop, op, f)
    fmpq_clear_readonly(f)

cdef void fmpq_poly_set_coeff_mpq(fmpq_poly_t poly, slong n, const mpq_t x) noexcept:
    cdef fmpq_t t
    fmpq_init_set_readonly(t, x)
    fmpq_poly_set_coeff_fmpq(poly, n, t)
    fmpq_clear_readonly(t)

cdef void fmpq_poly_get_coeff_mpq(mpq_t x, const fmpq_poly_t poly, slong n) noexcept:
    cdef fmpq_t t
    fmpq_init(t)
    fmpq_poly_get_coeff_fmpq(t, poly, n)
    fmpq_get_mpq(x, t)
    fmpq_clear(t)

cdef void fmpq_poly_set_mpq(fmpq_poly_t poly, const mpq_t x) noexcept:
    fmpq_poly_fit_length(poly, 1)
    fmpz_set_mpz(fmpq_poly_numref(poly), mpq_numref(x))
    fmpz_set_mpz(fmpq_poly_denref(poly), mpq_denref(x))
    _fmpq_poly_set_length(poly, 1)
    _fmpq_poly_normalise(poly)

cdef void fmpq_poly_set_mpz(fmpq_poly_t poly, const mpz_t x) noexcept:
    fmpq_poly_fit_length(poly, 1)
    fmpz_set_mpz(fmpq_poly_numref(poly), x)
    fmpz_one(fmpq_poly_denref(poly))
    _fmpq_poly_set_length(poly, 1)
    _fmpq_poly_normalise(poly)
