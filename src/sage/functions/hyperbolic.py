r"""
Hyperbolic functions

The full set of hyperbolic and inverse hyperbolic functions is
available:

 - hyperbolic sine: :class:`sinh() <Function_sinh>`
 - hyperbolic cosine: :class:`cosh() <Function_cosh>`
 - hyperbolic tangent: :class:`tanh() <Function_tanh>`
 - hyperbolic cotangent: :class:`coth() <Function_coth>`
 - hyperbolic secant: :class:`sech() <Function_sech>`
 - hyperbolic cosecant: :class:`csch() <Function_csch>`
 - inverse hyperbolic sine: :class:`asinh() <Function_arcsinh>`
 - inverse hyperbolic cosine: :class:`acosh() <Function_arccosh>`
 - inverse hyperbolic tangent: :class:`atanh() <Function_arctanh>`
 - inverse hyperbolic cotangent: :class:`acoth() <Function_arccoth>`
 - inverse hyperbolic secant: :class:`asech() <Function_arcsech>`
 - inverse hyperbolic cosecant: :class:`acsch() <Function_arccsch>`

REFERENCES:

- :wikipedia:`Hyperbolic function`
- :wikipedia:`Inverse hyperbolic functions`
- R. Roy, F. W. J. Olver, Elementary Functions, https://dlmf.nist.gov/4

EXAMPLES:

Inverse hyperbolic functions have logarithmic expressions,
so expressions of the form ``exp(c*f(x))`` simplify::

    sage: # needs sage.symbolic
    sage: exp(2*atanh(x))
    -(x + 1)/(x - 1)
    sage: exp(2*acoth(x))
    (x + 1)/(x - 1)
    sage: exp(2*asinh(x))
    (x + sqrt(x^2 + 1))^2
    sage: exp(2*acosh(x))
    (x + sqrt(x^2 - 1))^2
    sage: exp(2*asech(x))
    (sqrt(-x^2 + 1)/x + 1/x)^2
    sage: exp(2*acsch(x))
    (sqrt(1/x^2 + 1) + 1/x)^2
"""

from sage.symbolic.function import GinacFunction


class Function_sinh(GinacFunction):
    def __init__(self):
        r"""
        The hyperbolic sine function.

        EXAMPLES::

            sage: sinh(3.1415)                                                          # needs sage.rings.real_mpfr
            11.5476653707437

            sage: # needs sage.symbolic
            sage: sinh(pi)
            sinh(pi)
            sage: float(sinh(pi))
            11.54873935725774...
            sage: RR(sinh(pi))
            11.5487393572577
            sage: latex(sinh(x))
            \sinh\left(x\right)
            sage: sinh(x)._sympy_()                                                     # needs sympy
            sinh(x)

        To prevent automatic evaluation, use the ``hold`` parameter::

            sage: sinh(arccosh(x), hold=True)                                           # needs sage.symbolic
            sinh(arccosh(x))

        To then evaluate again, use the ``unhold`` method::

            sage: sinh(arccosh(x), hold=True).unhold()                                  # needs sage.symbolic
            sqrt(x + 1)*sqrt(x - 1)
        """
        GinacFunction.__init__(self, "sinh", latex_name=r"\sinh")


sinh = Function_sinh()


class Function_cosh(GinacFunction):
    def __init__(self):
        r"""
        The hyperbolic cosine function.

        EXAMPLES::

            sage: cosh(3.1415)                                                          # needs sage.rings.real_mpfr
            11.5908832931176

            sage: # needs sage.symbolic
            sage: cosh(pi)
            cosh(pi)
            sage: float(cosh(pi))
            11.591953275521519
            sage: RR(cosh(1/2))
            1.12762596520638
            sage: latex(cosh(x))
            \cosh\left(x\right)
            sage: cosh(x)._sympy_()                                                     # needs sympy
            cosh(x)

        To prevent automatic evaluation, use the ``hold`` parameter::

            sage: cosh(arcsinh(x), hold=True)                                           # needs sage.symbolic
            cosh(arcsinh(x))

        To then evaluate again, use the ``unhold`` method::

            sage: cosh(arcsinh(x), hold=True).unhold()                                  # needs sage.symbolic
            sqrt(x^2 + 1)
        """
        GinacFunction.__init__(self, "cosh", latex_name=r"\cosh")


cosh = Function_cosh()


class Function_tanh(GinacFunction):
    def __init__(self):
        r"""
        The hyperbolic tangent function.

        EXAMPLES::

            sage: tanh(3.1415)                                                          # needs sage.rings.real_mpfr
            0.996271386633702
            sage: tan(3.1415/4)                                                         # needs sage.rings.real_mpfr
            0.999953674278156

            sage: # needs sage.symbolic
            sage: tanh(pi)
            tanh(pi)
            sage: float(tanh(pi))
            0.99627207622075
            sage: tanh(pi/4)
            tanh(1/4*pi)
            sage: RR(tanh(1/2))
            0.462117157260010

        ::

            sage: CC(tanh(pi + I*e))                                                    # needs sage.rings.real_mpfr sage.symbolic
            0.997524731976164 - 0.00279068768100315*I
            sage: ComplexField(100)(tanh(pi + I*e))                                     # needs sage.rings.real_mpfr sage.symbolic
            0.99752473197616361034204366446 - 0.0027906876810031453884245163923*I
            sage: CDF(tanh(pi + I*e))  # rel tol 2e-15                                  # needs sage.rings.complex_double sage.symbolic
            0.9975247319761636 - 0.002790687681003147*I

        To prevent automatic evaluation, use the ``hold`` parameter::

            sage: tanh(arcsinh(x), hold=True)                                           # needs sage.symbolic
            tanh(arcsinh(x))

        To then evaluate again, use the ``unhold`` method::

            sage: tanh(arcsinh(x), hold=True).unhold()                                  # needs sage.symbolic
            x/sqrt(x^2 + 1)

        TESTS::

            sage: latex(tanh(x))                                                        # needs sage.symbolic
            \tanh\left(x\right)
            sage: tanh(x)._sympy_()                                                     # needs sympy sage.symbolic
            tanh(x)

        Check that real/imaginary parts are correct (:issue:`20098`)::

            sage: # needs sage.symbolic
            sage: tanh(1+2*I).n()
            1.16673625724092 - 0.243458201185725*I
            sage: tanh(1+2*I).real().n()
            1.16673625724092
            sage: tanh(1+2*I).imag().n()
            -0.243458201185725
            sage: tanh(x).real()
            sinh(2*real_part(x))/(cos(2*imag_part(x)) + cosh(2*real_part(x)))
            sage: tanh(x).imag()
            sin(2*imag_part(x))/(cos(2*imag_part(x)) + cosh(2*real_part(x)))
        """
        GinacFunction.__init__(self, "tanh", latex_name=r"\tanh")


tanh = Function_tanh()


class Function_coth(GinacFunction):
    def __init__(self):
        r"""
        The hyperbolic cotangent function.

        EXAMPLES::

            sage: coth(3.1415)                                                          # needs sage.rings.real_mpfr
            1.00374256795520
            sage: coth(complex(1, 2))  # abs tol 1e-15                                  # needs sage.rings.complex_double
            (0.8213297974938518+0.17138361290918508j)

            sage: # needs sage.symbolic
            sage: coth(pi)
            coth(pi)
            sage: coth(0)
            Infinity
            sage: coth(pi*I)
            Infinity
            sage: coth(pi*I/2)
            0
            sage: coth(7*pi*I/2)
            0
            sage: coth(8*pi*I/2)
            Infinity
            sage: coth(7.*pi*I/2)
            -I*cot(3.50000000000000*pi)
            sage: float(coth(pi))
            1.0037418731973213
            sage: RR(coth(pi))
            1.00374187319732

            sage: # needs sage.symbolic
            sage: bool(diff(coth(x), x) == diff(1/tanh(x), x))
            True
            sage: diff(coth(x), x)
            -1/sinh(x)^2
            sage: latex(coth(x))
            \coth\left(x\right)
            sage: coth(x)._sympy_()                                                     # needs sympy
            coth(x)
        """
        GinacFunction.__init__(self, "coth", latex_name=r"\coth")

    def _eval_numpy_(self, x):
        """
        EXAMPLES::

            sage: import numpy                                                          # needs numpy
            sage: a = numpy.arange(2, 5)                                                # needs numpy
            sage: coth(a)                                                               # needs numpy
            array([1.03731472, 1.00496982, 1.00067115])
        """
        return 1.0 / tanh(x)


coth = Function_coth()


class Function_sech(GinacFunction):
    def __init__(self):
        r"""
        The hyperbolic secant function.

        EXAMPLES::

            sage: sech(3.1415)                                                          # needs sage.rings.real_mpfr
            0.0862747018248192

            sage: # needs sage.symbolic
            sage: sech(pi)
            sech(pi)
            sage: float(sech(pi))
            0.0862667383340544...
            sage: RR(sech(pi))
            0.0862667383340544
            sage: sech(0)
            1
            sage: sech(pi*I)
            -1
            sage: sech(pi*I/2)
            Infinity
            sage: sech(7*pi*I/2)
            Infinity
            sage: sech(8*pi*I/2)
            1
            sage: sech(8.*pi*I/2)
            sec(4.00000000000000*pi)

            sage: # needs sage.symbolic
            sage: bool(diff(sech(x), x) == diff(1/cosh(x), x))
            True
            sage: diff(sech(x), x)
            -sech(x)*tanh(x)
            sage: latex(sech(x))
            \operatorname{sech}\left(x\right)
            sage: sech(x)._sympy_()                                                     # needs sympy
            sech(x)
        """
        GinacFunction.__init__(self, "sech", latex_name=r"\operatorname{sech}",)

    def _eval_numpy_(self, x):
        """
        EXAMPLES::

            sage: import numpy                                                          # needs numpy
            sage: a = numpy.arange(2, 5)                                                # needs numpy
            sage: sech(a)                                                               # needs numpy
            array([0.26580223, 0.09932793, 0.03661899])
        """
        return 1.0 / cosh(x)


sech = Function_sech()


class Function_csch(GinacFunction):
    def __init__(self):
        r"""
        The hyperbolic cosecant function.

        EXAMPLES::

            sage: csch(3.1415)                                                          # needs sage.rings.real_mpfr
            0.0865975907592133

            sage: # needs sage.symbolic
            sage: csch(pi)
            csch(pi)
            sage: float(csch(pi))
            0.0865895375300469...
            sage: RR(csch(pi))
            0.0865895375300470
            sage: csch(0)
            Infinity
            sage: csch(pi*I)
            Infinity
            sage: csch(pi*I/2)
            -I
            sage: csch(7*pi*I/2)
            I
            sage: csch(7.*pi*I/2)
            -I*csc(3.50000000000000*pi)

            sage: # needs sage.symbolic
            sage: bool(diff(csch(x), x) == diff(1/sinh(x), x))
            True
            sage: diff(csch(x), x)
            -coth(x)*csch(x)
            sage: latex(csch(x))
            \operatorname{csch}\left(x\right)
            sage: csch(x)._sympy_()                                                     # needs sympy
            csch(x)
        """
        GinacFunction.__init__(self, "csch", latex_name=r"\operatorname{csch}")

    def _eval_numpy_(self, x):
        """
        EXAMPLES::

            sage: import numpy                                                          # needs numpy
            sage: a = numpy.arange(2, 5)                                                # needs numpy
            sage: csch(a)                                                               # needs numpy
            array([0.27572056, 0.09982157, 0.03664357])
        """
        return 1.0 / sinh(x)


csch = Function_csch()


################################
# Inverse hyperbolic functions #
################################


class Function_arcsinh(GinacFunction):
    def __init__(self):
        r"""
        The inverse of the hyperbolic sine function.

        EXAMPLES::

            sage: asinh
            arcsinh
            sage: asinh(0.5)                                                            # needs sage.rings.real_mpfr
            0.481211825059603
            sage: asinh(1/2)                                                            # needs sage.symbolic
            arcsinh(1/2)
            sage: asinh(1 + I*1.0)                                                      # needs sage.symbolic
            1.06127506190504 + 0.666239432492515*I

        To prevent automatic evaluation use the ``hold`` argument::

            sage: asinh(-2, hold=True)                                                  # needs sage.symbolic
            arcsinh(-2)

        To then evaluate again, use the ``unhold`` method::

            sage: asinh(-2, hold=True).unhold()                                         # needs sage.symbolic
            -arcsinh(2)

        ``conjugate(asinh(x))==asinh(conjugate(x))`` unless on the branch
        cuts which run along the imaginary axis outside the interval [-I, +I].::

            sage: # needs sage.symbolic
            sage: conjugate(asinh(x))
            conjugate(arcsinh(x))
            sage: var('y', domain='positive')
            y
            sage: conjugate(asinh(y))
            arcsinh(y)
            sage: conjugate(asinh(y+I))
            conjugate(arcsinh(y + I))
            sage: conjugate(asinh(1/16))
            arcsinh(1/16)
            sage: conjugate(asinh(I/2))
            arcsinh(-1/2*I)
            sage: conjugate(asinh(2*I))
            conjugate(arcsinh(2*I))

        TESTS::

            sage: asinh(x).operator()                                                   # needs sage.symbolic
            arcsinh
            sage: latex(asinh(x))                                                       # needs sage.symbolic
            \operatorname{arsinh}\left(x\right)
            sage: asinh(x)._sympy_()                                                    # needs sympy sage.symbolic
            asinh(x)
        """
        GinacFunction.__init__(self, "arcsinh",
                latex_name=r"\operatorname{arsinh}",
                conversions=dict(maxima='asinh', sympy='asinh', fricas='asinh',
                                 giac='asinh', mathematica='ArcSinh'))


arcsinh = asinh = Function_arcsinh()


class Function_arccosh(GinacFunction):
    def __init__(self):
        r"""
        The inverse of the hyperbolic cosine function.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: acosh(1/2)
            arccosh(1/2)
            sage: acosh(1 + I*1.0)
            1.06127506190504 + 0.904556894302381*I
            sage: float(acosh(2))
            1.3169578969248168
            sage: cosh(float(acosh(2)))
            2.0

            sage: acosh(complex(1, 2))  # abs tol 1e-15                                 # needs sage.rings.complex_double
            (1.5285709194809982+1.1437177404024204j)

        .. warning::

            If the input is in the complex field or symbolic (which
            includes rational and integer input), the output will
            be complex.  However, if the input is a real decimal, the
            output will be real or `NaN`.  See the examples for details.

        ::

            sage: acosh(CC(0.5))                                                        # needs sage.rings.real_mpfr
            1.04719755119660*I

            sage: # needs sage.symbolic
            sage: acosh(0.5)
            NaN
            sage: acosh(1/2)
            arccosh(1/2)
            sage: acosh(1/2).n()
            NaN
            sage: acosh(0)
            1/2*I*pi
            sage: acosh(-1)
            I*pi

        To prevent automatic evaluation use the ``hold`` argument::

            sage: acosh(-1, hold=True)                                                  # needs sage.symbolic
            arccosh(-1)

        To then evaluate again, use the ``unhold`` method::

            sage: acosh(-1, hold=True).unhold()                                         # needs sage.symbolic
            I*pi

        ``conjugate(arccosh(x))==arccosh(conjugate(x))`` unless on the branch
        cut which runs along the real axis from +1 to -inf.::

            sage: # needs sage.symbolic
            sage: conjugate(acosh(x))
            conjugate(arccosh(x))
            sage: var('y', domain='positive')
            y
            sage: conjugate(acosh(y))
            conjugate(arccosh(y))
            sage: conjugate(acosh(y+I))
            conjugate(arccosh(y + I))
            sage: conjugate(acosh(1/16))
            conjugate(arccosh(1/16))
            sage: conjugate(acosh(2))
            arccosh(2)
            sage: conjugate(acosh(I/2))
            arccosh(-1/2*I)

        TESTS::

            sage: acosh(x).operator()                                                   # needs sage.symbolic
            arccosh
            sage: latex(acosh(x))                                                       # needs sage.symbolic
            \operatorname{arcosh}\left(x\right)
            sage: acosh(x)._sympy_()                                                    # needs sympy sage.symbolic
            acosh(x)
        """
        GinacFunction.__init__(self, "arccosh",
                latex_name=r"\operatorname{arcosh}",
                conversions=dict(maxima='acosh', sympy='acosh', fricas='acosh',
                                 giac='acosh', mathematica='ArcCosh'))


arccosh = acosh = Function_arccosh()


class Function_arctanh(GinacFunction):
    def __init__(self):
        r"""
        The inverse of the hyperbolic tangent function.

        EXAMPLES::

            sage: atanh(0.5)                                                            # needs sage.rings.real_mpfr
            0.549306144334055
            sage: atanh(1/2)                                                            # needs sage.symbolic
            1/2*log(3)
            sage: atanh(1 + I*1.0)                                                      # needs sage.symbolic
            0.402359478108525 + 1.01722196789785*I

        To prevent automatic evaluation use the ``hold`` argument::

            sage: atanh(-1/2, hold=True)                                                # needs sage.symbolic
            arctanh(-1/2)

        To then evaluate again, use the ``unhold`` method::

            sage: atanh(-1/2, hold=True).unhold()                                       # needs sage.symbolic
            -1/2*log(3)

        ``conjugate(arctanh(x)) == arctanh(conjugate(x))`` unless on the branch
        cuts which run along the real axis outside the interval [-1, +1]. ::

            sage: # needs sage.symbolic
            sage: conjugate(atanh(x))
            conjugate(arctanh(x))
            sage: var('y', domain='positive')
            y
            sage: conjugate(atanh(y))
            conjugate(arctanh(y))
            sage: conjugate(atanh(y + I))
            conjugate(arctanh(y + I))
            sage: conjugate(atanh(1/16))
            1/2*log(17/15)
            sage: conjugate(atanh(I/2))
            arctanh(-1/2*I)
            sage: conjugate(atanh(-2*I))
            arctanh(2*I)

        TESTS::

            sage: atanh(x).operator()                                                   # needs sage.symbolic
            arctanh
            sage: latex(atanh(x))                                                       # needs sage.symbolic
            \operatorname{artanh}\left(x\right)
            sage: atanh(x)._sympy_()                                                    # needs sympy sage.symbolic
            atanh(x)
        """
        GinacFunction.__init__(self, "arctanh",
                latex_name=r"\operatorname{artanh}",
                conversions=dict(maxima='atanh', sympy='atanh', fricas='atanh',
                                 giac='atanh', mathematica='ArcTanh'))


arctanh = atanh = Function_arctanh()


class Function_arccoth(GinacFunction):
    def __init__(self):
        r"""
        The inverse of the hyperbolic cotangent function.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: acoth(2.0)
            0.549306144334055
            sage: acoth(2)
            1/2*log(3)
            sage: acoth(1 + I*1.0)
            0.402359478108525 - 0.553574358897045*I
            sage: acoth(2).n(200)
            0.54930614433405484569762261846126285232374527891137472586735

            sage: bool(diff(acoth(x), x) == diff(atanh(x), x))                          # needs sage.symbolic
            True
            sage: diff(acoth(x), x)                                                     # needs sage.symbolic
            -1/(x^2 - 1)

            sage: float(acoth(2))                                                       # needs sage.symbolic
            0.5493061443340549
            sage: float(acoth(2).n(53))   # Correct result to 53 bits                   # needs sage.rings.real_mpfr sage.symbolic
            0.5493061443340549
            sage: float(acoth(2).n(100))  # Compute 100 bits and then round to 53       # needs sage.rings.real_mpfr sage.symbolic
            0.5493061443340549

        TESTS::

            sage: latex(acoth(x))                                                       # needs sage.symbolic
            \operatorname{arcoth}\left(x\right)
            sage: acoth(x)._sympy_()                                                    # needs sympy sage.symbolic
            acoth(x)

        Check that :issue:`23636` is fixed::

            sage: acoth(float(1.1))                                                     # needs sage.symbolic
            1.5222612188617113
        """
        GinacFunction.__init__(self, "arccoth",
                latex_name=r"\operatorname{arcoth}",
                conversions=dict(maxima='acoth', sympy='acoth',
                                 mathematica='ArcCoth',
                                 giac='acoth', fricas='acoth'))

    def _eval_numpy_(self, x):
        """
        EXAMPLES::

            sage: import numpy                                                          # needs numpy
            sage: a = numpy.arange(2,5)                                                 # needs numpy
            sage: acoth(a)                                                              # needs numpy
            array([0.54930614, 0.34657359, 0.25541281])
        """
        return arctanh(1.0 / x)


arccoth = acoth = Function_arccoth()


class Function_arcsech(GinacFunction):
    def __init__(self):
        r"""
        The inverse of the hyperbolic secant function.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: asech(0.5)
            1.31695789692482
            sage: asech(1/2)
            arcsech(1/2)
            sage: asech(1 + I*1.0)
            0.530637530952518 - 1.11851787964371*I
            sage: asech(1/2).n(200)
            1.3169578969248167086250463473079684440269819714675164797685
            sage: float(asech(1/2))
            1.3169578969248168

            sage: diff(asech(x), x)                                                     # needs sage.symbolic
            -1/(sqrt(-x^2 + 1)*x)
            sage: latex(asech(x))                                                       # needs sage.symbolic
            \operatorname{arsech}\left(x\right)
            sage: asech(x)._sympy_()                                                    # needs sympy sage.symbolic
            asech(x)
        """
        GinacFunction.__init__(self, "arcsech",
                latex_name=r"\operatorname{arsech}",
                conversions=dict(maxima='asech', sympy='asech',
                                 mathematica='ArcSech',
                                 fricas='asech'))

    def _eval_numpy_(self, x):
        """
        EXAMPLES::

            sage: import numpy                                                          # needs numpy
            sage: a = numpy.linspace(0,1,3)                                             # needs numpy
            sage: asech(a)                                                              # needs numpy
            doctest:...: RuntimeWarning: divide by zero encountered in ...divide
            array([       inf,  1.3169579,  0.       ])
        """
        return arccosh(1.0 / x)


arcsech = asech = Function_arcsech()


class Function_arccsch(GinacFunction):
    def __init__(self):
        r"""
        The inverse of the hyperbolic cosecant function.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: acsch(2.0)
            0.481211825059603
            sage: acsch(2)
            arccsch(2)
            sage: acsch(1 + I*1.0)
            0.530637530952518 - 0.452278447151191*I
            sage: acsch(1).n(200)
            0.88137358701954302523260932497979230902816032826163541075330
            sage: float(acsch(1))
            0.881373587019543

            sage: diff(acsch(x), x)                                                     # needs sage.symbolic
            -1/(sqrt(x^2 + 1)*x)
            sage: latex(acsch(x))                                                       # needs sage.symbolic
            \operatorname{arcsch}\left(x\right)

        TESTS:

        Check that :issue:`20818` is fixed::

            sage: acsch(float(0.1))                                                     # needs sage.symbolic
            2.99822295029797
            sage: acsch(x)._sympy_()                                                    # needs sympy sage.symbolic
            acsch(x)
        """
        GinacFunction.__init__(self, "arccsch",
                latex_name=r"\operatorname{arcsch}",
                conversions=dict(maxima='acsch',
                                 mathematica='ArcCsch',
                                 sympy='acsch', fricas='acsch'))

    def _eval_numpy_(self, x):
        """
        EXAMPLES::

            sage: import numpy                                                          # needs numpy
            sage: a = numpy.linspace(0,1,3)                                             # needs numpy
            sage: acsch(a)                                                              # needs numpy
            doctest:...: RuntimeWarning: divide by zero encountered in ...divide
            array([        inf,  1.44363548,  0.88137359])
        """
        return arcsinh(1.0 / x)


arccsch = acsch = Function_arccsch()
