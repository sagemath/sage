# sage_setup: distribution = sagemath-modules
from sage.rings.complex_double cimport ComplexDoubleElement
from sage.rings.complex_mpfr cimport ComplexNumber
from sage.libs.mpfr cimport mpfr_get_d, MPFR_RNDN
from sage.libs.gsl.complex cimport GSL_SET_COMPLEX

cdef class CCtoCDF(Map):

    cpdef Element _call_(self, x):
        """
        EXAMPLES::
            sage: from sage.rings.complex_conversion import CCtoCDF
            sage: f = CCtoCDF(CC, CDF)  # indirect doctest
            sage: f(CC.0)
            1.0*I
            sage: f(exp(pi*CC.0/4))                                                     # needs sage.symbolic
            0.7071067811865476 + 0.7071067811865475*I
        """
        z = <ComplexDoubleElement>ComplexDoubleElement.__new__(ComplexDoubleElement)
        GSL_SET_COMPLEX(&z._complex,
                        mpfr_get_d((<ComplexNumber>x).__re, MPFR_RNDN),
                        mpfr_get_d((<ComplexNumber>x).__im, MPFR_RNDN))
        return z
