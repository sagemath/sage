from sage.libs.flint.types cimport fmpz_poly_t
from sage.libs.ntl.types cimport ZZX_c
from sage.libs.mpfr cimport mpfr_t
from sage.libs.mpfi cimport mpfi_t

cdef fmpz_poly_evaluation_mpfr(mpfr_t res, const fmpz_poly_t poly, const mpfr_t a)
cdef fmpz_poly_evaluation_mpfi(mpfi_t res, const fmpz_poly_t poly, const mpfi_t a)

cdef ZZX_evaluation_mpfr(mpfr_t res, ZZX_c poly, const mpfr_t a)
cdef ZZX_evaluation_mpfi(mpfi_t res, ZZX_c poly, const mpfi_t a)

