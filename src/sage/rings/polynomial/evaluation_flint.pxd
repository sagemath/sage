# sage_setup: distribution = sagemath-flint

from sage.libs.flint.types cimport fmpz_poly_t
from sage.libs.mpfr.types cimport mpfr_t
from sage.libs.mpfi.types cimport mpfi_t

cdef fmpz_poly_evaluation_mpfr(mpfr_t res, const fmpz_poly_t poly, const mpfr_t a)
cdef fmpz_poly_evaluation_mpfi(mpfi_t res, const fmpz_poly_t poly, const mpfi_t a)
