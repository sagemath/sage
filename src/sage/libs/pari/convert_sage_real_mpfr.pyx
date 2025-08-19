# sage.doctest: needs sage.rings.real_mpfr

from cypari2.stack cimport new_gen
from cypari2.paridecl cimport *
from cysignals.signals cimport sig_on

from sage.libs.gmp.mpz cimport *
from sage.libs.mpfr cimport *
from sage.libs.mpfr.types cimport mpfr_prec_t
from sage.rings.real_mpfr cimport RealField_class, RealField


cpdef Gen new_gen_from_real_mpfr_element(RealNumber self):

    # This uses interfaces of MPFR and PARI which are documented
    # (and not marked subject-to-change).  It could be faster
    # by using internal interfaces of MPFR, which are documented
    # as subject-to-change.

    if mpfr_nan_p(self.value) or mpfr_inf_p(self.value):
        raise ValueError('cannot convert NaN or infinity to Pari float')

    # wordsize for PARI
    cdef unsigned long wordsize = sizeof(long)*8

    cdef mpfr_prec_t prec
    prec = (<RealField_class>self._parent)._prec

    # We round up the precision to the nearest multiple of wordsize.
    cdef int rounded_prec
    rounded_prec = nbits2prec(self.prec())

    # Yes, assigning to self works fine, even in Cython.
    if rounded_prec > prec:
        self = RealField(rounded_prec)(self)

    cdef mpz_t mantissa
    cdef mp_exp_t exponent
    cdef GEN pari_float

    sig_on()
    if mpfr_zero_p(self.value):
        pari_float = real_0_bit(-rounded_prec)
    else:
        # Now we can extract the mantissa, and it will be normalized
        # (the most significant bit of the most significant word will be 1).
        mpz_init(mantissa)
        exponent = mpfr_get_z_exp(mantissa, self.value)

        # Create a PARI REAL
        pari_float = cgetr(rounded_prec)
        pari_float[1] = evalexpo(exponent + rounded_prec - 1) + evalsigne(mpfr_sgn(self.value))
        mpz_export(&pari_float[2], NULL, 1, wordsize // 8, 0, 0, mantissa)
        mpz_clear(mantissa)

    return new_gen(pari_float)


cpdef bint set_real_mpfr_element_from_gen(RealNumber self, Gen x) noexcept:
    r"""
    EXAMPLES::

        sage: rt2 = sqrt(pari('2.0'))
        sage: rt2
        1.41421356237310
        sage: rt2.sage()
        1.41421356237309505
        sage: rt2.sage().prec()
        64
        sage: pari(rt2.sage()) == rt2
        True
        sage: for i in range(100, 200):
        ....:     assert(sqrt(pari(i)) == pari(sqrt(pari(i)).sage()))
        sage: (-3.1415).__pari__().sage()
        -3.14150000000000000
    """
    cdef GEN g = x.g

    if typ((<Gen>x).g) != t_REAL:
        return False

    cdef int sgn
    sgn = signe(g)

    if sgn == 0:
        mpfr_set_ui(self.value, 0, MPFR_RNDN)
        return True

    cdef int wordsize = 8 * sizeof(long)

    cdef mpz_t mantissa
    mpz_init(mantissa)
    mpz_import(mantissa, lg(g) - 2, 1, wordsize/8, 0, 0, &g[2])

    cdef mp_exp_t exponent = expo(g)

    # Round to nearest for best results when setting a low-precision
    # MPFR from a high-precision GEN
    mpfr_set_z(self.value, mantissa, MPFR_RNDN)
    mpfr_mul_2si(self.value, self.value, exponent - wordsize * (lg(g) - 2) + 1, MPFR_RNDN)

    if sgn < 0:
        mpfr_neg(self.value, self.value, MPFR_RNDN)

    mpz_clear(mantissa)

    return True
