# sage.doctest: needs sage.rings.complex_double

from cysignals.signals cimport sig_on, sig_off

from sage.libs.gsl.complex cimport *

from cypari2.paridecl cimport *
from cypari2.convert cimport new_gen_from_double, new_t_COMPLEX_from_double


cpdef ComplexDoubleElement pari_to_cdf(Gen g):
    """
    Create a CDF element from a PARI ``gen``.

    EXAMPLES::

        sage: CDF(pari("Pi"))  # indirect doctest
        3.141592653589793
        sage: CDF(pari("1 + I/2"))
        1.0 + 0.5*I

    TESTS:

    Check that we handle PARI errors gracefully, see :issue:`17329`::

        sage: CDF(-151.386325246 + 992.34771962*I).zeta()
        Traceback (most recent call last):
        ...
        PariError: overflow in t_REAL->double conversion
        sage: CDF(pari(x^2 + 5))                                                        # needs sage.symbolic
        Traceback (most recent call last):
        ...
        PariError: incorrect type in gtofp (t_POL)
    """
    cdef ComplexDoubleElement z = ComplexDoubleElement.__new__(ComplexDoubleElement)
    sig_on()
    if typ(g.g) == t_COMPLEX:
        z._complex = gsl_complex_rect(gtodouble(gel(g.g, 1)), gtodouble(gel(g.g, 2)))
    else:
        z._complex = gsl_complex_rect(gtodouble(g.g), 0.0)
    sig_off()
    return z


cpdef Gen new_gen_from_complex_double_element(ComplexDoubleElement self):
    """
    Return PARI version of ``self``, as ``t_COMPLEX`` or ``t_REAL``.

    EXAMPLES::

        sage: from sage.libs.pari.convert_sage_complex_double import new_gen_from_complex_double_element
        sage: new_gen_from_complex_double_element(CDF(1,0))
        1.00000000000000
        sage: new_gen_from_complex_double_element(CDF(0,1))
        1.00000000000000*I
        sage: new_gen_from_complex_double_element(CDF(1,1))
        1.00000000000000 + 1.00000000000000*I
    """
    if not GSL_IMAG(self._complex):
        return new_gen_from_double(GSL_REAL(self._complex))
    else:
        return new_t_COMPLEX_from_double(GSL_REAL(self._complex), GSL_IMAG(self._complex))


cpdef ComplexDoubleElement complex_double_element_eta(ComplexDoubleElement self, int flag):
    """
    TESTS::

        sage: from sage.libs.pari.convert_sage_complex_double import complex_double_element_eta
        sage: a = CDF(1,1)
        sage: complex_double_element_eta(a, 0)
        0.9981290699259585
        sage: complex_double_element_eta(a, 1)
        0.7420487758365647 + 0.1988313702299107*I
    """
    return pari_to_cdf(new_gen_from_complex_double_element(self).eta(flag))


cpdef ComplexDoubleElement complex_double_element_agm(ComplexDoubleElement self, right):
    """
    TESTS::

        sage: from sage.libs.pari.convert_sage_complex_double import complex_double_element_agm
        sage: complex_double_element_agm(CDF(1, 1), CDF(2, 2))
        1.4567910310469068 + 1.4567910310469068*I
    """
    return pari_to_cdf(new_gen_from_complex_double_element(self).agm(right))


cpdef ComplexDoubleElement complex_double_element_dilog(ComplexDoubleElement self):
    """
    TESTS::

        sage: from sage.libs.pari.convert_sage_complex_double import complex_double_element_dilog
        sage: complex_double_element_dilog(CDF(1, 1))
        0.6168502750680849 + 1.4603621167531196*I
    """
    return pari_to_cdf(new_gen_from_complex_double_element(self).dilog())


cpdef ComplexDoubleElement complex_double_element_gamma(ComplexDoubleElement self):
    """
    TESTS::

        sage: from sage.libs.pari.convert_sage_complex_double import complex_double_element_gamma
        sage: complex_double_element_gamma(CDF(1, 1))
        0.49801566811835607 - 0.15494982830181067*I
    """
    return pari_to_cdf(new_gen_from_complex_double_element(self).gamma())


cpdef ComplexDoubleElement complex_double_element_gamma_inc(ComplexDoubleElement self, t):
    """
    TESTS::

        sage: from sage.libs.pari.convert_sage_complex_double import complex_double_element_gamma_inc
        sage: complex_double_element_gamma_inc(CDF(1, 1), CDF(2, 2))
        0.054695987717541285 - 0.04676800059213122*I
    """
    return pari_to_cdf(new_gen_from_complex_double_element(self).incgam(t))


cpdef ComplexDoubleElement complex_double_element_zeta(ComplexDoubleElement self):
    """
    TESTS::

        sage: from sage.libs.pari.convert_sage_complex_double import complex_double_element_zeta
        sage: complex_double_element_zeta(CDF(1, 2))
        0.5981655697623818 - 0.35185474521784527*I
    """
    return pari_to_cdf(new_gen_from_complex_double_element(self).zeta())
