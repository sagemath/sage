# sage_setup: distribution = sagemath-ntl
from sage.libs.ntl.types cimport ZZX_c
from sage.libs.mpfr.types cimport mpfr_t
from sage.libs.mpfi.types cimport mpfi_t

cdef ZZX_evaluation_mpfr(mpfr_t res, ZZX_c poly, const mpfr_t a)
cdef ZZX_evaluation_mpfi(mpfi_t res, ZZX_c poly, const mpfi_t a)
