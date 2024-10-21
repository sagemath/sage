# sage_setup: distribution = sagemath-modules
from sage.libs.gmp.types cimport mp_limb_t

cdef extern from "mpfr.h":
    ctypedef int mpfr_sign_t
    ctypedef long mpfr_prec_t
    ctypedef long mpfr_exp_t
    ctypedef unsigned long mpfr_uexp_t

    ctypedef struct __mpfr_struct:
        mpfr_prec_t _mpfr_prec
        mpfr_sign_t _mpfr_sign
        mpfr_exp_t _mpfr_exp
        mp_limb_t* _mpfr_d

    ctypedef __mpfr_struct mpfr_t[1]
    ctypedef __mpfr_struct* mpfr_ptr
    ctypedef __mpfr_struct* mpfr_srcptr
    ctypedef enum mpfr_rnd_t:
        MPFR_RNDN
        MPFR_RNDZ
        MPFR_RNDU
        MPFR_RNDD
        MPFR_RNDA
        MPFR_RNDF
        MPFR_RNDNA
