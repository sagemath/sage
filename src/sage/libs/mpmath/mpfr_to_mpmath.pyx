# sage_setup: distribution = sagemath-mpmath

from sage.ext.stdsage cimport PY_NEW
from sage.libs.mpfr cimport *
from sage.libs.gmp.all cimport *
from sage.rings.integer cimport Integer

cdef object finf, fnan, fninf, fzero


cdef inline late_import():
    global finf, fnan, fninf, fzero
    if finf is None:
        from sage.libs.mpmath._vendor.mpmath.libmp import finf, fnan, fninf, fzero


cdef mpfr_to_mpfval(mpfr_t value):
    """
    Given an MPFR value, return an mpmath mpf data tuple representing
    the same number.
    """
    late_import()
    if mpfr_nan_p(value):
        return fnan
    if mpfr_inf_p(value):
        if mpfr_sgn(value) > 0:
            return finf
        else:
            return fninf
    if mpfr_sgn(value) == 0:
        return fzero
    sign = 0
    cdef Integer man = PY_NEW(Integer)
    exp = mpfr_get_z_exp(man.value, value)
    if mpz_sgn(man.value) < 0:
        mpz_neg(man.value, man.value)
        sign = 1
    cdef unsigned long trailing
    trailing = mpz_scan1(man.value, 0)
    if trailing:
        mpz_tdiv_q_2exp(man.value, man.value, trailing)
        exp += trailing
    bc = mpz_sizeinbase(man.value, 2)
    return (sign, man, int(exp), bc)
