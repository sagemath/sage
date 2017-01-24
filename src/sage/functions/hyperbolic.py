"""
Hyperbolic Functions
"""
from __future__ import division

from sage.symbolic.function import GinacFunction, BuiltinFunction
from sage.rings.infinity import Infinity
from sage.symbolic.expression import Expression
from sage.symbolic.constants import pi, I
from sage.rings.integer_ring import ZZ

import math


class HyperbolicFunction(BuiltinFunction):
    r"""
    Abstract base class for the functions defined in this file.

    EXAMPLES::

        sage: from sage.functions.hyperbolic import HyperbolicFunction
        sage: f = HyperbolicFunction('foo', latex_name='\\foo', conversions={'mathematica':'Foo'},evalf_float=lambda x: 2*x)
        sage: f(x)
        foo(x)
        sage: f(0.5r)
        1.0
        sage: latex(f(x))
        \foo\left(x\right)
        sage: f(x)._mathematica_init_()
        'Foo[x]'
    """
    def __init__(self, name, latex_name=None, conversions=None,
            evalf_float=None):
        """
        Note that subclasses of HyperbolicFunction should be instantiated only
        once, since they inherit from BuiltinFunction which only uses the
        name and class to check if a function was already registered.

        EXAMPLES::

            sage: from sage.functions.hyperbolic import HyperbolicFunction
            sage: class Barh(HyperbolicFunction):
            ....:     def __init__(self):
            ....:         HyperbolicFunction.__init__(self, 'barh')
            sage: barh = Barh()
            sage: barh(x)
            barh(x)
        """
        self._evalf_float = evalf_float
        BuiltinFunction.__init__(self, name, latex_name=latex_name,
                conversions=conversions)

    def _evalf_(self, x, parent=None, algorithm=None):
        """
        EXAMPLES::

            sage: from sage.functions.hyperbolic import HyperbolicFunction
            sage: class Fooh(HyperbolicFunction):
            ....:     def __init__(self):
            ....:         HyperbolicFunction.__init__(self, 'fooh',evalf_float=lambda x: 2*x)
            sage: fooh = Fooh()
            sage: fooh(float(5))
            10.0
            sage: fooh(0.5r)
            1.0
            sage: fooh(x).subs(x=.5r)
            1.0
            sage: fooh(x).n()
            Traceback (most recent call last):
            ...
            AttributeError: 'sage.symbolic.expression.Expression' object has no attribute 'fooh'
        """
        if parent is float:
            return self._evalf_float(x)
        return getattr(x, self.name())()


class Function_sinh(GinacFunction):
    def __init__(self):
        r"""
        The hyperbolic sine function.

        EXAMPLES::

            sage: sinh(pi)
            sinh(pi)
            sage: sinh(3.1415)
            11.5476653707437
            sage: float(sinh(pi))
            11.54873935725774...
            sage: RR(sinh(pi))
            11.5487393572577

            sage: latex(sinh(x))
            \sinh\left(x\right)
            sage: sinh(x)._sympy_()
            sinh(x)

        To prevent automatic evaluation, use the ``hold`` parameter::

            sage: sinh(arccosh(x),hold=True)
            sinh(arccosh(x))

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: sinh(arccosh(x),hold=True).simplify()
            sqrt(x + 1)*sqrt(x - 1)

        """
        GinacFunction.__init__(self, "sinh", latex_name=r"\sinh")

sinh = Function_sinh()


class Function_cosh(GinacFunction):
    def __init__(self):
        r"""
        The hyperbolic cosine function.

        EXAMPLES::

            sage: cosh(pi)
            cosh(pi)
            sage: cosh(3.1415)
            11.5908832931176
            sage: float(cosh(pi))
            11.591953275521519
            sage: RR(cosh(1/2))
            1.12762596520638

            sage: latex(cosh(x))
            \cosh\left(x\right)
            sage: cosh(x)._sympy_()
            cosh(x)

        To prevent automatic evaluation, use the ``hold`` parameter::

            sage: cosh(arcsinh(x),hold=True)
            cosh(arcsinh(x))

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: cosh(arcsinh(x),hold=True).simplify()
            sqrt(x^2 + 1)

        """
        GinacFunction.__init__(self, "cosh", latex_name=r"\cosh")

cosh = Function_cosh()


class Function_tanh(GinacFunction):
    def __init__(self):
        r"""
        The hyperbolic tangent function.

        EXAMPLES::

            sage: tanh(pi)
            tanh(pi)
            sage: tanh(3.1415)
            0.996271386633702
            sage: float(tanh(pi))
            0.99627207622075
            sage: tan(3.1415/4)
            0.999953674278156
            sage: tanh(pi/4)
            tanh(1/4*pi)
            sage: RR(tanh(1/2))
            0.462117157260010

        ::

            sage: CC(tanh(pi + I*e))
            0.997524731976164 - 0.00279068768100315*I
            sage: ComplexField(100)(tanh(pi + I*e))
            0.99752473197616361034204366446 - 0.0027906876810031453884245163923*I
            sage: CDF(tanh(pi + I*e))  # rel tol 2e-15
            0.9975247319761636 - 0.002790687681003147*I

        To prevent automatic evaluation, use the ``hold`` parameter::

            sage: tanh(arcsinh(x),hold=True)
            tanh(arcsinh(x))

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: tanh(arcsinh(x),hold=True).simplify()
            x/sqrt(x^2 + 1)

        TESTS::

            sage: latex(tanh(x))
            \tanh\left(x\right)
            sage: tanh(x)._sympy_()
            tanh(x)

        Check that real/imaginary parts are correct (:trac:`20098`)::

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
            sage: coth(3.1415)
            1.00374256795520
            sage: float(coth(pi))
            1.0037418731973213
            sage: RR(coth(pi))
            1.00374187319732
            sage: coth(complex(1, 2))  # abs tol 1e-15
            (0.8213297974938518+0.17138361290918508j)

            sage: bool(diff(coth(x), x) == diff(1/tanh(x), x))
            True
            sage: diff(coth(x), x)
            -1/sinh(x)^2
            sage: latex(coth(x))
            \operatorname{coth}\left(x\right)
            sage: coth(x)._sympy_()
            coth(x)
        """
        GinacFunction.__init__(self, "coth", latex_name=r"\operatorname{coth}")

    def _eval_numpy_(self, x):
        """
        EXAMPLES::

            sage: import numpy
            sage: a = numpy.arange(2, 5)
            sage: coth(a)
            array([ 1.03731472,  1.00496982,  1.00067115])
        """
        return 1.0 / tanh(x)

coth = Function_coth()


class Function_sech(GinacFunction):
    def __init__(self):
        r"""
        The hyperbolic secant function.

        EXAMPLES::

            sage: sech(pi)
            sech(pi)
            sage: sech(3.1415)
            0.0862747018248192
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

            sage: bool(diff(sech(x), x) == diff(1/cosh(x), x))
            True
            sage: diff(sech(x), x)
            -sech(x)*tanh(x)
            sage: latex(sech(x))
            \operatorname{sech}\left(x\right)
            sage: sech(x)._sympy_()
            sech(x)
        """
        GinacFunction.__init__(self, "sech", latex_name=r"\operatorname{sech}",)

    def _eval_numpy_(self, x):
        """
        EXAMPLES::

            sage: import numpy
            sage: a = numpy.arange(2, 5)
            sage: sech(a)
            array([ 0.26580223,  0.09932793,  0.03661899])
        """
        return 1.0 / cosh(x)

sech = Function_sech()


class Function_csch(GinacFunction):
    def __init__(self):
        r"""
        The hyperbolic cosecant function.

        EXAMPLES::

            sage: csch(pi)
            csch(pi)
            sage: csch(3.1415)
            0.0865975907592133
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

            sage: bool(diff(csch(x), x) == diff(1/sinh(x), x))
            True
            sage: diff(csch(x), x)
            -coth(x)*csch(x)
            sage: latex(csch(x))
            {\rm csch}\left(x\right)
            sage: csch(x)._sympy_()
            csch(x)
        """
        GinacFunction.__init__(self, "csch", latex_name=r"{\rm csch}")

    def _eval_numpy_(self, x):
        """
        EXAMPLES::

            sage: import numpy
            sage: a = numpy.arange(2, 5)
            sage: csch(a)
            array([ 0.27572056,  0.09982157,  0.03664357])
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

            sage: arcsinh
            arcsinh
            sage: arcsinh(0.5)
            0.481211825059603
            sage: arcsinh(1/2)
            arcsinh(1/2)
            sage: arcsinh(1 + I*1.0)
            1.06127506190504 + 0.666239432492515*I

        To prevent automatic evaluation use the ``hold`` argument::

            sage: arcsinh(-2,hold=True)
            arcsinh(-2)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: arcsinh(-2,hold=True).simplify()
            -arcsinh(2)

        ``conjugate(arcsinh(x))==arcsinh(conjugate(x))`` unless on the branch
        cuts which run along the imaginary axis outside the interval [-I, +I].::

            sage: conjugate(arcsinh(x))
            conjugate(arcsinh(x))
            sage: var('y', domain='positive')
            y
            sage: conjugate(arcsinh(y))
            arcsinh(y)
            sage: conjugate(arcsinh(y+I))
            conjugate(arcsinh(y + I))
            sage: conjugate(arcsinh(1/16))
            arcsinh(1/16)
            sage: conjugate(arcsinh(I/2))
            arcsinh(-1/2*I)
            sage: conjugate(arcsinh(2*I))
            conjugate(arcsinh(2*I))

        TESTS::

            sage: arcsinh(x).operator()
            arcsinh
            sage: latex(arcsinh(x))
            {\rm arcsinh}\left(x\right)
            sage: asinh(x)._sympy_()
            asinh(x)
        """
        GinacFunction.__init__(self, "arcsinh", latex_name=r"{\rm arcsinh}",
                conversions=dict(maxima='asinh', sympy='asinh', fricas='asinh'))

arcsinh = asinh = Function_arcsinh()


class Function_arccosh(GinacFunction):
    def __init__(self):
        r"""
        The inverse of the hyperbolic cosine function.

        EXAMPLES::

            sage: arccosh(1/2)
            arccosh(1/2)
            sage: arccosh(1 + I*1.0)
            1.06127506190504 + 0.904556894302381*I
            sage: float(arccosh(2))
            1.3169578969248168
            sage: cosh(float(arccosh(2)))
            2.0
            sage: arccosh(complex(1, 2))  # abs tol 1e-15
            (1.5285709194809982+1.1437177404024204j)

        .. warning::

            If the input is in the complex field or symbolic (which
            includes rational and integer input), the output will
            be complex.  However, if the input is a real decimal, the
            output will be real or `NaN`.  See the examples for details.

        ::

            sage: arccosh(0.5)
            NaN
            sage: arccosh(1/2)
            arccosh(1/2)
            sage: arccosh(1/2).n()
            NaN
            sage: arccosh(CC(0.5))
            1.04719755119660*I
            sage: arccosh(0)
            1/2*I*pi
            sage: arccosh(-1)
            I*pi

        To prevent automatic evaluation use the ``hold`` argument::

            sage: arccosh(-1,hold=True)
            arccosh(-1)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: arccosh(-1,hold=True).simplify()
            I*pi

        ``conjugate(arccosh(x))==arccosh(conjugate(x))`` unless on the branch
        cut which runs along the real axis from +1 to -inf.::

            sage: conjugate(arccosh(x))
            conjugate(arccosh(x))
            sage: var('y', domain='positive')
            y
            sage: conjugate(arccosh(y))
            conjugate(arccosh(y))
            sage: conjugate(arccosh(y+I))
            conjugate(arccosh(y + I))
            sage: conjugate(arccosh(1/16))
            conjugate(arccosh(1/16))
            sage: conjugate(arccosh(2))
            arccosh(2)
            sage: conjugate(arccosh(I/2))
            arccosh(-1/2*I)

        TESTS::

            sage: arccosh(x).operator()
            arccosh
            sage: latex(arccosh(x))
            {\rm arccosh}\left(x\right)
            sage: acosh(x)._sympy_()
            acosh(x)
        """
        GinacFunction.__init__(self, "arccosh", latex_name=r"{\rm arccosh}",
                conversions=dict(maxima='acosh', sympy='acosh', fricas='acosh'))

arccosh = acosh = Function_arccosh()


class Function_arctanh(GinacFunction):
    def __init__(self):
        r"""
        The inverse of the hyperbolic tangent function.

        EXAMPLES::

            sage: arctanh(0.5)
            0.549306144334055
            sage: arctanh(1/2)
            arctanh(1/2)
            sage: arctanh(1 + I*1.0)
            0.402359478108525 + 1.01722196789785*I

        To prevent automatic evaluation use the ``hold`` argument::

            sage: arctanh(-1/2,hold=True)
            arctanh(-1/2)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: arctanh(-1/2,hold=True).simplify()
            -arctanh(1/2)

        ``conjugate(arctanh(x))==arctanh(conjugate(x))`` unless on the branch
        cuts which run along the real axis outside the interval [-1, +1].::

            sage: conjugate(arctanh(x))
            conjugate(arctanh(x))
            sage: var('y', domain='positive')
            y
            sage: conjugate(arctanh(y))
            conjugate(arctanh(y))
            sage: conjugate(arctanh(y+I))
            conjugate(arctanh(y + I))
            sage: conjugate(arctanh(1/16))
            arctanh(1/16)
            sage: conjugate(arctanh(I/2))
            arctanh(-1/2*I)
            sage: conjugate(arctanh(-2*I))
            arctanh(2*I)

        TESTS::

            sage: arctanh(x).operator()
            arctanh
            sage: latex(arctanh(x))
            {\rm arctanh}\left(x\right)
            sage: atanh(x)._sympy_()
            atanh(x)
        """
        GinacFunction.__init__(self, "arctanh", latex_name=r"{\rm arctanh}",
                conversions=dict(maxima='atanh', sympy='atanh', fricas='atanh'))

arctanh = atanh = Function_arctanh()


class Function_arccoth(GinacFunction):
    def __init__(self):
        r"""
        The inverse of the hyperbolic cotangent function.

        EXAMPLES::

            sage: arccoth(2.0)
            0.549306144334055
            sage: arccoth(2)
            arccoth(2)
            sage: arccoth(1 + I*1.0)
            0.402359478108525 - 0.553574358897045*I
            sage: arccoth(2).n(200)
            0.54930614433405484569762261846126285232374527891137472586735

            sage: bool(diff(acoth(x), x) == diff(atanh(x), x))
            True
            sage: diff(acoth(x), x)
            -1/(x^2 - 1)

        Using first the `.n(53)` method is slightly more precise than
        converting directly to a ``float``::

            sage: float(arccoth(2))  # abs tol 1e-16
            0.5493061443340548
            sage: float(arccoth(2).n(53))   # Correct result to 53 bits
            0.5493061443340549
            sage: float(arccoth(2).n(100))  # Compute 100 bits and then round to 53
            0.5493061443340549

        TESTS::

            sage: latex(arccoth(x))
            \operatorname{arccoth}\left(x\right)
            sage: acoth(x)._sympy_()
            acoth(x)
        """
        GinacFunction.__init__(self, "arccoth",
                latex_name=r"\operatorname{arccoth}",
                conversions=dict(maxima='acoth', sympy='acoth', fricas='acoth'))

    def _eval_numpy_(self, x):
        """
        EXAMPLES::

            sage: import numpy
            sage: a = numpy.arange(2,5)
            sage: arccoth(a)
            array([ 0.54930614,  0.34657359,  0.25541281])
        """
        return arctanh(1.0 / x)

arccoth = acoth = Function_arccoth()


class Function_arcsech(GinacFunction):
    def __init__(self):
        r"""
        The inverse of the hyperbolic secant function.

        EXAMPLES::

            sage: arcsech(0.5)
            1.31695789692482
            sage: arcsech(1/2)
            arcsech(1/2)
            sage: arcsech(1 + I*1.0)
            0.530637530952518 - 1.11851787964371*I
            sage: arcsech(1/2).n(200)
            1.3169578969248167086250463473079684440269819714675164797685
            sage: float(arcsech(1/2))
            1.3169578969248168

            sage: diff(asech(x), x)
            -1/(sqrt(-x^2 + 1)*x)
            sage: latex(arcsech(x))
            \operatorname{arcsech}\left(x\right)
            sage: asech(x)._sympy_()
            asech(x)
        """
        GinacFunction.__init__(self, "arcsech",
                latex_name=r"\operatorname{arcsech}",
                conversions=dict(maxima='asech', sympy='asech', fricas='asech'))

    def _eval_numpy_(self, x):
        """
        EXAMPLES::

            sage: import numpy
            sage: a = numpy.linspace(0,1,3)
            sage: arcsech(a)
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

            sage: arccsch(2.0)
            0.481211825059603
            sage: arccsch(2)
            arccsch(2)
            sage: arccsch(1 + I*1.0)
            0.530637530952518 - 0.452278447151191*I
            sage: arccsch(1).n(200)
            0.88137358701954302523260932497979230902816032826163541075330
            sage: float(arccsch(1))
            0.881373587019543

            sage: diff(acsch(x), x)
            -1/(sqrt(x^2 + 1)*x)
            sage: latex(arccsch(x))
            \operatorname{arccsch}\left(x\right)

        TESTS:

        Check if :trac:`20818` is fixed::

            sage: arccsch(float(0.1))
            2.99822295029797
        """
        GinacFunction.__init__(self, "arccsch",
                latex_name=r"\operatorname{arccsch}",
                conversions=dict(maxima='acsch', sympy='asech', fricas='acsch'))

    def _eval_numpy_(self, x):
        """
        EXAMPLES::

            sage: import numpy
            sage: a = numpy.linspace(0,1,3)
            sage: arccsch(a)
            doctest:...: RuntimeWarning: divide by zero encountered in ...divide
            array([        inf,  1.44363548,  0.88137359])
        """
        return arcsinh(1.0 / x)

arccsch = acsch = Function_arccsch()
