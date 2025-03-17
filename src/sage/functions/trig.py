r"""
Trigonometric functions
"""
import math

from sage.symbolic.function import GinacFunction


class Function_sin(GinacFunction):
    def __init__(self):
        r"""
        The sine function.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: sin(0)
            0
            sage: sin(x).subs(x==0)
            0
            sage: sin(2).n(100)
            0.90929742682568169539601986591
            sage: sin(x)._sympy_()                                                      # needs sympy
            sin(x)

        We can prevent evaluation using the ``hold`` parameter::

            sage: sin(0, hold=True)                                                     # needs sage.symbolic
            sin(0)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: a = sin(0, hold=True); a.simplify()                                   # needs sage.symbolic
            0

        If possible, the argument is also reduced modulo the
        period length `2\pi`, and well-known identities are
        directly evaluated::

            sage: k = var('k', domain='integer')                                        # needs sage.symbolic
            sage: sin(1 + 2*k*pi)                                                       # needs sage.symbolic
            sin(1)
            sage: sin(k*pi)                                                             # needs sage.symbolic
            0

        TESTS::

            sage: loads(dumps(sin))
            sin

            sage: conjugate(sin(x))                                                     # needs sage.symbolic
            sin(conjugate(x))
            sage: sin(complex(1,1))     # rel tol 1e-15
            (1.2984575814159773+0.6349639147847361j)

            sage: # needs sage.symbolic
            sage: sin(pi/5)
            1/4*sqrt(-2*sqrt(5) + 10)
            sage: sin(pi/8)
            1/2*sqrt(-sqrt(2) + 2)
            sage: sin(pi/24)
            1/4*sqrt(-2*sqrt(6) - 2*sqrt(2) + 8)
            sage: sin(pi/30)
            -1/8*sqrt(5) + 1/4*sqrt(-3/2*sqrt(5) + 15/2) - 1/8
            sage: sin(104*pi/105)
            sin(1/105*pi)
            sage: cos(pi/8)
            1/2*sqrt(sqrt(2) + 2)
            sage: cos(pi/10)
            1/4*sqrt(2*sqrt(5) + 10)
            sage: cos(pi/12)
            1/4*sqrt(6) + 1/4*sqrt(2)
            sage: cos(pi/15)
            1/8*sqrt(5) + 1/4*sqrt(3/2*sqrt(5) + 15/2) - 1/8
            sage: cos(pi/24)
            1/4*sqrt(2*sqrt(6) + 2*sqrt(2) + 8)
            sage: cos(104*pi/105)
            -cos(1/105*pi)
            sage: tan(pi/5)
            sqrt(-2*sqrt(5) + 5)
            sage: tan(pi/8)
            sqrt(2) - 1
            sage: tan(pi/10)
            1/5*sqrt(-10*sqrt(5) + 25)
            sage: tan(pi/16)
            -sqrt(2) + sqrt(2*sqrt(2) + 4) - 1
            sage: tan(pi/20)
            sqrt(5) - sqrt(2*sqrt(5) + 5) + 1
            sage: tan(pi/24)
            sqrt(6) - sqrt(3) + sqrt(2) - 2
            sage: tan(104*pi/105)
            -tan(1/105*pi)
            sage: cot(104*pi/105)
            -cot(1/105*pi)
            sage: sec(104*pi/105)
            -sec(1/105*pi)
            sage: csc(104*pi/105)
            csc(1/105*pi)

            sage: all(sin(rat*pi).n(200) - sin(rat*pi, hold=True).n(200) < 1e-30        # needs sage.symbolic
            ....:     for rat in [1/5, 2/5, 1/30, 7/30, 11/30, 13/30,
            ....:                 1/8, 3/8, 1/24, 5/24, 7/24, 11/24])
            True
            sage: all(cos(rat*pi).n(200)-cos(rat*pi, hold=True).n(200) < 1e-30          # needs sage.symbolic
            ....:     for rat in [1/10, 3/10, 1/12, 5/12, 1/15, 2/15, 4/15, 7/15,
            ....:                 1/8, 3/8, 1/24, 5/24, 11/24])
            True
            sage: all(tan(rat*pi).n(200)-tan(rat*pi, hold=True).n(200) < 1e-30          # needs sage.symbolic
            ....:     for rat in [1/5, 2/5, 1/10, 3/10, 1/20, 3/20, 7/20, 9/20,
            ....:                 1/8, 3/8, 1/16, 3/16, 5/16, 7/16, 1/24, 5/24, 7/24, 11/24])
            True

        Check that :issue:`20456` is fixed::

            sage: assume(x > 0)                                                         # needs sage.symbolic
            sage: sin(pi*x)                                                             # needs sage.symbolic
            sin(pi*x)
            sage: forget()                                                              # needs sage.symbolic

        Check that :issue:`20752` is fixed::

            sage: sin(3*pi + 41/42*pi)                                                  # needs sage.symbolic
            -sin(1/42*pi)
            sage: sin(-5*pi + 1/42*pi)                                                  # needs sage.symbolic
            -sin(1/42*pi)
            sage: sin(pi - 1/42*pi)                                                     # needs sage.symbolic
            sin(1/42*pi)
        """
        GinacFunction.__init__(self, 'sin', latex_name=r"\sin",
                conversions=dict(maxima='sin', mathematica='Sin',
                                 giac='sin', fricas='sin', sympy='sin'))


sin = Function_sin()


class Function_cos(GinacFunction):
    def __init__(self):
        r"""
        The cosine function.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: cos(pi)
            -1
            sage: cos(x).subs(x==pi)
            -1
            sage: cos(2).n(100)
            -0.41614683654714238699756822950
            sage: cos(x)._sympy_()                                                      # needs sympy
            cos(x)

        We can prevent evaluation using the ``hold`` parameter::

            sage: cos(0, hold=True)                                                     # needs sage.symbolic
            cos(0)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: a = cos(0, hold=True); a.simplify()                                   # needs sage.symbolic
            1

        If possible, the argument is also reduced modulo the
        period length `2\pi`, and well-known identities are
        directly evaluated::

            sage: # needs sage.symbolic
            sage: k = var('k', domain='integer')
            sage: cos(1 + 2*k*pi)
            cos(1)
            sage: cos(k*pi)
            cos(pi*k)
            sage: cos(pi/3 + 2*k*pi)
            1/2

        TESTS::

            sage: loads(dumps(cos))
            cos

            sage: conjugate(cos(x))                                                     # needs sage.symbolic
            cos(conjugate(x))
            sage: cos(complex(1,1))     # rel tol 1e-15
            (0.8337300251311491-0.9888977057628651j)

        Check that :issue:`20752` is fixed::

            sage: cos(3*pi + 41/42*pi)                                                  # needs sage.symbolic
            cos(1/42*pi)
            sage: cos(-5*pi + 1/42*pi)                                                  # needs sage.symbolic
            -cos(1/42*pi)
            sage: cos(pi - 1/42*pi)                                                     # needs sage.symbolic
            -cos(1/42*pi)
        """
        GinacFunction.__init__(self, 'cos', latex_name=r"\cos",
                conversions=dict(maxima='cos', mathematica='Cos',
                                 giac='cos', fricas='cos', sympy='cos'))


cos = Function_cos()


class Function_tan(GinacFunction):
    def __init__(self):
        r"""
        The tangent function.

        EXAMPLES::

            sage: # needs sage.rings.real_mpfr
            sage: tan(3.1415)
            -0.0000926535900581913
            sage: tan(3.1415/4)
            0.999953674278156

            sage: # needs sage.symbolic
            sage: tan(pi)
            0
            sage: tan(pi/4)
            1
            sage: tan(1/2)
            tan(1/2)
            sage: RR(tan(1/2))
            0.546302489843790

        We can prevent evaluation using the ``hold`` parameter::

            sage: tan(pi/4, hold=True)                                                  # needs sage.symbolic
            tan(1/4*pi)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: a = tan(pi/4, hold=True); a.simplify()                                # needs sage.symbolic
            1

        If possible, the argument is also reduced modulo the
        period length `\pi`, and well-known identities are
        directly evaluated::

            sage: k = var('k', domain='integer')                                        # needs sage.symbolic
            sage: tan(1 + 2*k*pi)                                                       # needs sage.symbolic
            tan(1)
            sage: tan(k*pi)                                                             # needs sage.symbolic
            0

        TESTS::

            sage: tan(x)._sympy_()                                                      # needs sympy sage.symbolic
            tan(x)
            sage: conjugate(tan(x))                                                     # needs sage.symbolic
            tan(conjugate(x))
            sage: tan(complex(1,1))     # rel tol 1e-15
            (0.2717525853195118+1.0839233273386946j)

        Check that :issue:`19791` is fixed::

            sage: tan(2+I).imag().n()                                                   # needs sage.symbolic
            1.16673625724092
        """
        GinacFunction.__init__(self, 'tan', latex_name=r"\tan",
                               conversions=dict(maxima='tan', mathematica='Tan',
                                 giac='tan', fricas='tan', sympy='tan'))


tan = Function_tan()


class Function_cot(GinacFunction):
    def __init__(self):
        r"""
        The cotangent function.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: cot(pi/4)
            1
            sage: RR(cot(pi/4))
            1.00000000000000
            sage: cot(1/2)
            cot(1/2)
            sage: cot(0.5)
            1.83048772171245

            sage: latex(cot(x))                                                         # needs sage.symbolic
            \cot\left(x\right)
            sage: cot(x)._sympy_()                                                      # needs sympy sage.symbolic
            cot(x)

        We can prevent evaluation using the ``hold`` parameter::

            sage: cot(pi/4, hold=True)                                                  # needs sage.symbolic
            cot(1/4*pi)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: a = cot(pi/4, hold=True); a.simplify()                                # needs sage.symbolic
            1

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: cot(pi/4)
            1
            sage: cot(x).subs(x==pi/4)
            1
            sage: cot(pi/7)
            cot(1/7*pi)
            sage: cot(x)
            cot(x)

            sage: # needs sage.symbolic
            sage: n(cot(pi/4), 100)
            1.0000000000000000000000000000
            sage: float(cot(1))
            0.64209261593433...
            sage: bool(diff(cot(x), x) == diff(1/tan(x), x))
            True
            sage: diff(cot(x), x)
            -cot(x)^2 - 1

        TESTS::

            sage: # needs sage.symbolic
            sage: cot(float(0))
            Infinity
            sage: cot(SR(0))
            Infinity
            sage: cot(float(0.1))
            9.966644423259238
            sage: type(_)
            <... 'float'>
            sage: cot(float(0))
            Infinity
            sage: cot(SR(0))
            Infinity
            sage: cot(float(0.1))
            9.966644423259238
            sage: type(_)
            <... 'float'>

        Test complex input::

            sage: cot(complex(1,1))     # rel tol 1e-15                                 # needs sage.rings.complex_double
            (0.21762156185440273-0.8680141428959249j)
            sage: cot(1.+I)                                                             # needs sage.symbolic
            0.217621561854403 - 0.868014142895925*I
        """
        GinacFunction.__init__(self, 'cot', latex_name=r"\cot",
                               conversions=dict(maxima='cot', mathematica='Cot',
                                 giac='cot', fricas='cot', sympy='cot'))

    def _eval_numpy_(self, x):
        """
        EXAMPLES::

             sage: import numpy                                                         # needs numpy
             sage: a = numpy.arange(2, 5)                                               # needs numpy
             sage: cot(a)                                                               # needs numpy
             array([-0.45765755, -7.01525255,  0.86369115])
        """
        return 1.0 / tan(x)


cot = Function_cot()


class Function_sec(GinacFunction):
    def __init__(self):
        r"""
        The secant function.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: sec(pi/4)
            sqrt(2)
            sage: sec(x).subs(x==pi/4)
            sqrt(2)
            sage: sec(pi/7)
            sec(1/7*pi)
            sage: sec(x)
            sec(x)
            sage: RR(sec(pi/4))
            1.41421356237310
            sage: n(sec(pi/4),100)
            1.4142135623730950488016887242
            sage: float(sec(pi/4))
            1.4142135623730951
            sage: sec(1/2)
            sec(1/2)
            sage: sec(0.5)
            1.13949392732455

            sage: # needs sage.symbolic
            sage: bool(diff(sec(x), x) == diff(1/cos(x), x))
            True
            sage: diff(sec(x), x)
            sec(x)*tan(x)
            sage: latex(sec(x))
            \sec\left(x\right)
            sage: sec(x)._sympy_()                                                      # needs sympy
            sec(x)

        We can prevent evaluation using the ``hold`` parameter::

            sage: sec(pi/4, hold=True)                                                  # needs sage.symbolic
            sec(1/4*pi)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: a = sec(pi/4, hold=True); a.simplify()                                # needs sage.symbolic
            sqrt(2)

        TESTS:

        Test complex input::

            sage: sec(complex(1,1))     # rel tol 1e-15                                 # needs sage.rings.complex_double
            (0.49833703055518686+0.5910838417210451j)
        """
        GinacFunction.__init__(self, 'sec', latex_name=r"\sec",
                               conversions=dict(maxima='sec', mathematica='Sec',
                                 giac='sec', fricas='sec', sympy='sec'))

    def _eval_numpy_(self, x):
        """
        EXAMPLES::

            sage: import numpy                                                          # needs numpy
            sage: a = numpy.arange(2, 5)                                                # needs numpy
            sage: sec(a)                                                                # needs numpy
            array([-2.40299796, -1.01010867, -1.52988566])
        """
        return 1 / cos(x)


sec = Function_sec()


class Function_csc(GinacFunction):
    def __init__(self):
        r"""
        The cosecant function.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: csc(pi/4)
            sqrt(2)
            sage: csc(x).subs(x==pi/4)
            sqrt(2)
            sage: csc(pi/7)
            csc(1/7*pi)
            sage: csc(x)
            csc(x)
            sage: RR(csc(pi/4))
            1.41421356237310
            sage: n(csc(pi/4), 100)
            1.4142135623730950488016887242
            sage: float(csc(pi/4))
            1.4142135623730951
            sage: csc(1/2)
            csc(1/2)
            sage: csc(0.5)
            2.08582964293349

            sage: # needs sage.symbolic
            sage: bool(diff(csc(x), x) == diff(1/sin(x), x))
            True
            sage: diff(csc(x), x)
            -cot(x)*csc(x)
            sage: latex(csc(x))
            \csc\left(x\right)
            sage: csc(x)._sympy_()                                                      # needs sympy
            csc(x)

        We can prevent evaluation using the ``hold`` parameter::

            sage: csc(pi/4, hold=True)                                                  # needs sage.symbolic
            csc(1/4*pi)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: a = csc(pi/4,hold=True); a.simplify()                                 # needs sage.symbolic
            sqrt(2)

        TESTS:

        Test complex input::

            sage: csc(complex(1,1))     # rel tol 1e-15                                 # needs sage.rings.complex_double
            (0.6215180171704284-0.30393100162842646j)
        """
        GinacFunction.__init__(self, 'csc', latex_name=r"\csc",
                               conversions=dict(maxima='csc', mathematica='Csc',
                                 giac='csc', fricas='csc', sympy='csc'))

    def _eval_numpy_(self, x):
        """
        EXAMPLES::

            sage: import numpy                                                          # needs numpy
            sage: a = numpy.arange(2, 5)                                                # needs numpy
            sage: csc(a)                                                                # needs numpy
            array([ 1.09975017,  7.0861674 , -1.32134871])
        """
        return 1 / sin(x)


csc = Function_csc()


###################################
# Inverse Trigonometric Functions #
###################################

class Function_arcsin(GinacFunction):
    def __init__(self):
        """
        The arcsine function.

        EXAMPLES::

            sage: arcsin(0.5)                                                           # needs sage.rings.real_mpfr
            0.523598775598299
            sage: arcsin(1/2)                                                           # needs sage.symbolic
            1/6*pi
            sage: arcsin(1 + 1.0*I)                                                     # needs sage.symbolic
            0.666239432492515 + 1.06127506190504*I

        We can delay evaluation using the ``hold`` parameter::

            sage: arcsin(0, hold=True)                                                  # needs sage.symbolic
            arcsin(0)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: a = arcsin(0, hold=True); a.simplify()                                # needs sage.symbolic
            0

        ``conjugate(arcsin(x))==arcsin(conjugate(x))``, unless on the branch
        cuts which run along the real axis outside the interval [-1, +1].::

            sage: # needs sage.symbolic
            sage: conjugate(arcsin(x))
            conjugate(arcsin(x))
            sage: var('y', domain='positive')
            y
            sage: conjugate(arcsin(y))
            conjugate(arcsin(y))
            sage: conjugate(arcsin(y+I))
            conjugate(arcsin(y + I))
            sage: conjugate(arcsin(1/16))
            arcsin(1/16)
            sage: conjugate(arcsin(2))
            conjugate(arcsin(2))
            sage: conjugate(arcsin(-2))
            -conjugate(arcsin(2))

        TESTS::

            sage: arcsin(x)._sympy_()                                                   # needs sympy sage.symbolic
            asin(x)
            sage: arcsin(x).operator()                                                  # needs sage.symbolic
            arcsin
            sage: asin(complex(1,1))                                                    # needs sage.rings.complex_double
            (0.6662394324925152+1.0612750619050357j)
            sage: asin(SR(2.1))                                                         # needs sage.symbolic
            1.57079632679490 - 1.37285914424258*I
        """
        GinacFunction.__init__(self, 'arcsin', latex_name=r"\arcsin",
                conversions=dict(maxima='asin', sympy='asin',
                                 mathematica='ArcSin',
                                 fricas='asin', giac='asin'))


arcsin = asin = Function_arcsin()


class Function_arccos(GinacFunction):
    def __init__(self):
        """
        The arccosine function.

        EXAMPLES::

            sage: arccos(0.5)                                                           # needs sage.rings.real_mpfr
            1.04719755119660
            sage: arccos(1/2)                                                           # needs sage.symbolic
            1/3*pi
            sage: arccos(1 + 1.0*I)                                                     # needs sage.symbolic
            0.904556894302381 - 1.06127506190504*I
            sage: arccos(3/4).n(100)                                                    # needs sage.symbolic
            0.72273424781341561117837735264

        We can delay evaluation using the ``hold`` parameter::

            sage: arccos(0, hold=True)                                                  # needs sage.symbolic
            arccos(0)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: a = arccos(0, hold=True); a.simplify()                                # needs sage.symbolic
            1/2*pi

        ``conjugate(arccos(x))==arccos(conjugate(x))``, unless on the branch
        cuts, which run along the real axis outside the interval [-1, +1].::

            sage: # needs sage.symbolic
            sage: conjugate(arccos(x))
            conjugate(arccos(x))
            sage: var('y', domain='positive')
            y
            sage: conjugate(arccos(y))
            conjugate(arccos(y))
            sage: conjugate(arccos(y+I))
            conjugate(arccos(y + I))
            sage: conjugate(arccos(1/16))
            arccos(1/16)
            sage: conjugate(arccos(2))
            conjugate(arccos(2))
            sage: conjugate(arccos(-2))
            pi - conjugate(arccos(2))

        TESTS::

            sage: arccos(x)._sympy_()                                                   # needs sympy sage.symbolic
            acos(x)
            sage: arccos(x).operator()                                                  # needs sage.symbolic
            arccos
            sage: acos(complex(1,1))                                                    # needs sage.rings.complex_double
            (0.9045568943023814-1.0612750619050357j)
            sage: acos(SR(2.1))                                                         # needs sage.symbolic
            1.37285914424258*I

            sage: arcsin(sqrt(2)/2)
            1/4*pi
        """
        GinacFunction.__init__(self, 'arccos', latex_name=r"\arccos",
                conversions=dict(maxima='acos', sympy='acos',
                                 mathematica='ArcCos',
                                 fricas='acos', giac='acos'))


arccos = acos = Function_arccos()


class Function_arctan(GinacFunction):
    def __init__(self):
        """
        The arctangent function.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: arctan(1/2)
            arctan(1/2)
            sage: RDF(arctan(1/2))  # rel tol 1e-15
            0.46364760900080615
            sage: arctan(1 + I)
            arctan(I + 1)
            sage: arctan(1/2).n(100)
            0.46364760900080611621425623146

        We can delay evaluation using the ``hold`` parameter::

            sage: arctan(0, hold=True)                                                  # needs sage.symbolic
            arctan(0)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: a = arctan(0, hold=True); a.simplify()                                # needs sage.symbolic
            0

        ``conjugate(arctan(x))==arctan(conjugate(x))``, unless on the branch
        cuts which run along the imaginary axis outside the interval [-I, +I].::

            sage: # needs sage.symbolic
            sage: conjugate(arctan(x))
            conjugate(arctan(x))
            sage: var('y', domain='positive')
            y
            sage: conjugate(arctan(y))
            arctan(y)
            sage: conjugate(arctan(y+I))
            conjugate(arctan(y + I))
            sage: conjugate(arctan(1/16))
            arctan(1/16)
            sage: conjugate(arctan(-2*I))
            conjugate(arctan(-2*I))
            sage: conjugate(arctan(2*I))
            conjugate(arctan(2*I))
            sage: conjugate(arctan(I/2))
            arctan(-1/2*I)

        TESTS::

            sage: arctan(x)._sympy_()                                                   # needs sympy sage.symbolic
            atan(x)
            sage: arctan(x).operator()                                                  # needs sage.symbolic
            arctan
            sage: atan(complex(1,1))                                                    # needs sage.rings.complex_double
            (1.0172219678978514+0.4023594781085251j)

        Check that :issue:`19918` is fixed::

            sage: arctan(-x).subs(x=oo)                                                 # needs sage.symbolic
            -1/2*pi
            sage: arctan(-x).subs(x=-oo)                                                # needs sage.symbolic
            1/2*pi
        """
        GinacFunction.__init__(self, 'arctan', latex_name=r"\arctan",
                conversions=dict(maxima='atan', sympy='atan',
                                 mathematica='ArcTan',
                                 fricas='atan', giac='atan'))


arctan = atan = Function_arctan()


class Function_arccot(GinacFunction):
    def __init__(self):
        """
        The arccotangent function.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: arccot(1/2)
            arccot(1/2)
            sage: RDF(arccot(1/2))  # abs tol 2e-16
            1.1071487177940906
            sage: arccot(1 + I)
            arccot(I + 1)
            sage: arccot(1/2).n(100)
            1.1071487177940905030170654602
            sage: float(arccot(1/2))  # abs tol 2e-16
            1.1071487177940906
            sage: bool(diff(acot(x), x) == -diff(atan(x), x))
            True
            sage: diff(acot(x), x)
            -1/(x^2 + 1)

        We can delay evaluation using the ``hold`` parameter::

            sage: arccot(1, hold=True)                                                  # needs sage.symbolic
            arccot(1)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: a = arccot(1, hold=True); a.simplify()                                # needs sage.symbolic
            1/4*pi

        TESTS:

        Test complex input::

            sage: arccot(x)._sympy_()                                                   # needs sympy sage.symbolic
            acot(x)
            sage: arccot(complex(1,1))  # rel tol 1e-15                                 # needs sage.rings.complex_double
            (0.5535743588970452-0.4023594781085251j)
            sage: arccot(1.+I)                                                          # needs sage.symbolic
            0.553574358897045 - 0.402359478108525*I
        """
        GinacFunction.__init__(self, 'arccot', latex_name=r"\operatorname{arccot}",
                conversions=dict(maxima='acot', sympy='acot',
                                 mathematica='ArcCot',
                                 fricas='acot', giac='acot'))

    def _eval_numpy_(self, x):
        """
        EXAMPLES::

            sage: import numpy                                                          # needs numpy
            sage: a = numpy.arange(2, 5)                                                # needs numpy
            sage: arccot(a)                                                             # needs numpy
            array([0.46364761, 0.32175055, 0.24497866])
        """
        return math.pi / 2 - arctan(x)


arccot = acot = Function_arccot()


class Function_arccsc(GinacFunction):
    def __init__(self):
        """
        The arccosecant function.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: arccsc(2)
            arccsc(2)
            sage: RDF(arccsc(2))  # rel tol 1e-15
            0.5235987755982988
            sage: arccsc(2).n(100)
            0.52359877559829887307710723055
            sage: float(arccsc(2))
            0.52359877559829...
            sage: arccsc(1 + I)
            arccsc(I + 1)
            sage: diff(acsc(x), x)
            -1/(sqrt(x^2 - 1)*x)
            sage: arccsc(x)._sympy_()                                                   # needs sympy
            acsc(x)

        We can delay evaluation using the ``hold`` parameter::

            sage: arccsc(1, hold=True)                                                  # needs sage.symbolic
            arccsc(1)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: a = arccsc(1, hold=True); a.simplify()                                # needs sage.symbolic
            1/2*pi

        TESTS:

        Test complex input::

            sage: arccsc(complex(1,1))  # rel tol 1e-15                                 # needs sage.rings.complex_double
            (0.45227844715119064-0.5306375309525178j)
        """
        GinacFunction.__init__(self, 'arccsc', latex_name=r"\operatorname{arccsc}",
                               conversions=dict(maxima='acsc', sympy='acsc',
                                                mathematica='ArcCsc',
                                                fricas='acsc', giac='acsc'))

    def _eval_numpy_(self, x):
        """
        EXAMPLES::

            sage: import numpy                                                          # needs numpy
            sage: a = numpy.arange(2, 5)                                                # needs numpy
            sage: arccsc(a)                                                             # needs numpy
            array([0.52359878, 0.33983691, 0.25268026])
        """
        return arcsin(1.0 / x)


arccsc = acsc = Function_arccsc()


class Function_arcsec(GinacFunction):
    def __init__(self):
        """
        The arcsecant function.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: arcsec(2)
            arcsec(2)
            sage: arcsec(2.0)
            1.04719755119660
            sage: arcsec(2).n(100)
            1.0471975511965977461542144611
            sage: arcsec(1/2).n(100)
            1.3169578969248167086250463473*I
            sage: RDF(arcsec(2))  # abs tol 1e-15
            1.0471975511965976
            sage: arcsec(1 + I)
            arcsec(I + 1)
            sage: diff(asec(x), x)
            1/(sqrt(x^2 - 1)*x)
            sage: arcsec(x)._sympy_()                                                   # needs sympy
            asec(x)

        We can delay evaluation using the ``hold`` parameter::

            sage: arcsec(1, hold=True)                                                  # needs sage.symbolic
            arcsec(1)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: a = arcsec(1, hold=True); a.simplify()                                # needs sage.symbolic
            0

        TESTS:

        Test complex input::

            sage: arcsec(complex(1,1))  # rel tol 1e-15                                 # needs sage.rings.complex_double
            (1.118517879643706+0.5306375309525178j)
        """
        GinacFunction.__init__(self, 'arcsec', latex_name=r"\operatorname{arcsec}",
                               conversions=dict(maxima='asec', sympy='asec',
                                                mathematica='ArcSec',
                                                fricas='asec', giac='asec'))

    def _eval_numpy_(self, x):
        """
        EXAMPLES::

            sage: import numpy                                                          # needs numpy
            sage: a = numpy.arange(2, 5)                                                # needs numpy
            sage: arcsec(a)                                                             # needs numpy
            array([1.04719755, 1.23095942, 1.31811607])
        """
        return arccos(1.0 / x)


arcsec = asec = Function_arcsec()


class Function_arctan2(GinacFunction):
    def __init__(self):
        r"""
        The modified arctangent function.

        Returns the arc tangent (measured in radians) of `y/x`, where
        unlike ``arctan(y/x)``, the signs of both ``x`` and ``y`` are
        considered.  In particular, this function measures the angle
        of a ray through the origin and `(x,y)`, with the positive
        `x`-axis the zero mark, and with output angle `\theta`
        being between `-\pi<\theta<=\pi`.

        Hence, ``arctan2(y,x) = arctan(y/x)`` only for `x>0`.  One
        may consider the usual arctan to measure angles of lines
        through the origin, while the modified function measures
        rays through the origin.

        Note that the `y`-coordinate is by convention the first input.

        EXAMPLES:

        Note the difference between the two functions::

            sage: arctan2(1, -1)                                                        # needs sage.symbolic
            3/4*pi
            sage: arctan(1/-1)                                                          # needs sage.symbolic
            -1/4*pi

        This is consistent with Python and Maxima::

            sage: maxima.atan2(1, -1)                                                   # needs sage.symbolic
            (3*%pi)/4
            sage: math.atan2(1, -1)
            2.356194490192345

        More examples::

            sage: arctan2(1, 0)                                                         # needs sage.symbolic
            1/2*pi
            sage: arctan2(2, 3)                                                         # needs sage.symbolic
            arctan(2/3)
            sage: arctan2(-1, -1)                                                       # needs sage.symbolic
            -3/4*pi

        Of course we can approximate as well::

            sage: arctan2(-1/2, 1).n(100)                                               # needs sage.symbolic
            -0.46364760900080611621425623146
            sage: arctan2(2, 3).n(100)                                                  # needs sage.symbolic
            0.58800260354756755124561108063

        We can delay evaluation using the ``hold`` parameter::

            sage: arctan2(-1/2, 1, hold=True)                                           # needs sage.symbolic
            arctan2(-1/2, 1)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: arctan2(-1/2, 1, hold=True).simplify()                                # needs sage.symbolic
            -arctan(1/2)

        The function also works with numpy arrays as input::

            sage: # needs numpy
            sage: import numpy
            sage: a = numpy.linspace(1, 3, 3)
            sage: b = numpy.linspace(3, 6, 3)
            sage: atan2(a, b)
            array([0.32175055, 0.41822433, 0.46364761])

            sage: atan2(1,a)                                                            # needs numpy
            array([0.78539816, 0.46364761, 0.32175055])

            sage: atan2(a, 1)                                                           # needs numpy
            array([0.78539816, 1.10714872, 1.24904577])

        TESTS::

            sage: x,y = var('x,y')                                                      # needs sage.symbolic
            sage: arctan2(y, x).operator()                                              # needs sage.symbolic
            arctan2

        Check if :issue:`8565` is fixed::

            sage: atan2(-pi, 0)                                                         # needs sage.symbolic
            -1/2*pi

        Check if :issue:`8564` is fixed::

            sage: arctan2(x,x)._sympy_()                                                # needs sympy sage.symbolic
            atan2(x, x)

        Check if numerical evaluation works :issue:`9913`::

            sage: arctan2(0, -log(2)).n()                                               # needs sage.symbolic
            3.14159265358979

        Check that atan2(0,0) returns NaN :issue:`21614`::

            sage: # needs sage.symbolic
            sage: atan2(0, 0)
            NaN
            sage: atan2(0, 0).n()
            NaN
            sage: atan2(0, 0, hold=True)
            arctan2(0, 0)
            sage: atan2(0, 0, hold=True).n()
            Traceback (most recent call last):
            ...
            RuntimeError: atan2(): division by zero

        Check if :issue:`10062` is fixed, this was caused by
        ``(I*I).is_positive()`` returning ``True``::

            sage: arctan2(0, I*I)                                                       # needs sage.symbolic
            pi
        """
        GinacFunction.__init__(self, 'arctan2', nargs=2, latex_name=r"\arctan",
                conversions=dict(maxima='atan2', sympy='atan2', giac='atan2'))


arctan2 = atan2 = Function_arctan2()
