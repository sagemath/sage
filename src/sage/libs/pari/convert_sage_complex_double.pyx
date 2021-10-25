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

    Check that we handle PARI errors gracefully, see :trac:`17329`::

        sage: CDF(-151.386325246 + 992.34771962*I).zeta()
        Traceback (most recent call last):
        ...
        PariError: overflow in t_REAL->double conversion
        sage: CDF(pari(x^2 + 5))
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
    if not self._complex.imag:
        return new_gen_from_double(self._complex.real)
    else:
        return new_t_COMPLEX_from_double(self._complex.real, self._complex.imag)
