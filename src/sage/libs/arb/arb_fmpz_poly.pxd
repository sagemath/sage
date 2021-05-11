# distutils: libraries = gmp flint ARB_LIBRARY
# distutils: depends = arb_fmpz_poly.h

from sage.libs.arb.types cimport *
from sage.libs.flint.types cimport fmpz, fmpz_poly_t

# acb.h
cdef extern from "arb_wrap.h":

    void _arb_fmpz_poly_evaluate_arb_horner(arb_t res, const fmpz * poly, long len, const arb_t x, long prec)
    void arb_fmpz_poly_evaluate_arb_horner(arb_t res, const fmpz_poly_t poly, const arb_t x, long prec)
    void _arb_fmpz_poly_evaluate_arb_rectangular(arb_t res, const fmpz * poly, long len, const arb_t x, long prec)
    void arb_fmpz_poly_evaluate_arb_rectangular(arb_t res, const fmpz_poly_t poly, const arb_t x, long prec)
    void _arb_fmpz_poly_evaluate_arb(arb_t res, const fmpz * poly, long len, const arb_t x, long prec)
    void arb_fmpz_poly_evaluate_arb(arb_t res, const fmpz_poly_t poly, const arb_t x, long prec)
    void _arb_fmpz_poly_evaluate_acb_horner(acb_t res, const fmpz * poly, long len, const acb_t x, long prec)
    void arb_fmpz_poly_evaluate_acb_horner(acb_t res, const fmpz_poly_t poly, const acb_t x, long prec)
    void _arb_fmpz_poly_evaluate_acb_rectangular(acb_t res, const fmpz * poly, long len, const acb_t x, long prec)
    void arb_fmpz_poly_evaluate_acb_rectangular(acb_t res, const fmpz_poly_t poly, const acb_t x, long prec)
    void _arb_fmpz_poly_evaluate_acb(acb_t res, const fmpz * poly, long len, const acb_t x, long prec)
    void arb_fmpz_poly_evaluate_acb(acb_t res, const fmpz_poly_t poly, const acb_t x, long prec)
    unsigned long arb_fmpz_poly_deflation(const fmpz_poly_t poly)
    void arb_fmpz_poly_deflate(fmpz_poly_t res, const fmpz_poly_t poly, unsigned long deflation)
    void arb_fmpz_poly_complex_roots(acb_ptr roots, const fmpz_poly_t poly, int flags, long prec)
    void arb_fmpz_poly_cos_minpoly(fmpz_poly_t res, unsigned long n)
    void arb_fmpz_poly_gauss_period_minpoly(fmpz_poly_t res, unsigned long q, unsigned long n)
