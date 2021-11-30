r"""
Exponential Integrals

AUTHORS:

- Benjamin Jones (2011-06-12)

This module provides easy access to many exponential integral
special functions. It utilizes Maxima's `special functions package`_ and
the `mpmath library`_.

REFERENCES:

- [AS1964]_ Abramowitz and Stegun: *Handbook of Mathematical Functions*
- :wikipedia:`Exponential_integral`
- Online Encyclopedia of Special Function: http://algo.inria.fr/esf/index.html
- NIST Digital Library of Mathematical Functions: https://dlmf.nist.gov/
- Maxima `special functions package`_
- `mpmath library`_

.. _`special functions package`: http://maxima.sourceforge.net/docs/manual/en/maxima_15.html
.. _`mpmath library`: https://github.com/fredrik-johansson/mpmath/

AUTHORS:

- Benjamin Jones

    Implementations of the classes ``Function_exp_integral_*``.

- David Joyner and William Stein

    Authors of the code which was moved from special.py and trans.py.
    Implementation of :meth:`exp_int` (from sage/functions/special.py).
    Implementation of :meth:`exponential_integral_1` (from
    sage/functions/transcendental.py).

"""

#*****************************************************************************
#       Copyright (C) 2011 Benjamin Jones <benjaminfjones@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.symbolic.function import BuiltinFunction
from sage.symbolic.expression import Expression
from sage.structure.all import parent
from sage.misc.latex import latex
from sage.libs.mpmath import utils as mpmath_utils
mpmath_utils_call = mpmath_utils.call # eliminate some overhead in _evalf_

from sage.rings.real_mpfr import RealField
from sage.rings.integer_ring import ZZ
from sage.functions.log import exp, log
from sage.functions.trig import sin, cos
from sage.functions.hyperbolic import sinh, cosh
import math


class Function_exp_integral_e(BuiltinFunction):
    r"""
    The generalized complex exponential integral `E_n(z)` defined by

    .. MATH::

        E_n(z) = \int_1^{\infty} \frac{e^{-z t}}{t^n} \; dt

    for complex numbers `n` and `z`, see [AS1964]_ 5.1.4.

    The special case where `n = 1` is denoted in Sage by
    ``exp_integral_e1``.

    EXAMPLES:

    Numerical evaluation is handled using mpmath::

        sage: N(exp_integral_e(1,1))
        0.219383934395520
        sage: exp_integral_e(1, RealField(100)(1))
        0.21938393439552027367716377546

    We can compare this to PARI's evaluation of
    :meth:`exponential_integral_1`::

        sage: N(exponential_integral_1(1))
        0.219383934395520

    We can verify one case of [AS1964]_ 5.1.45, i.e.
    `E_n(z) = z^{n-1}\Gamma(1-n,z)`::

        sage: N(exp_integral_e(2, 3+I))
        0.00354575823814662 - 0.00973200528288687*I
        sage: N((3+I)*gamma(-1, 3+I))
        0.00354575823814662 - 0.00973200528288687*I

    Maxima returns the following improper integral as a multiple of
    ``exp_integral_e(1,1)``::

        sage: uu = integral(e^(-x)*log(x+1),x,0,oo)
        sage: uu
        e*exp_integral_e(1, 1)
        sage: uu.n(digits=30)
        0.596347362323194074341078499369

    Symbolic derivatives and integrals are handled by Sage and Maxima::

        sage: x = var('x')
        sage: f = exp_integral_e(2,x)
        sage: f.diff(x)
        -exp_integral_e(1, x)

        sage: f.integrate(x)
        -exp_integral_e(3, x)

        sage: f = exp_integral_e(-1,x)
        sage: f.integrate(x)
        Ei(-x) - gamma(-1, x)

    Some special values of ``exp_integral_e`` can be simplified.
    [AS1964]_ 5.1.23::

        sage: exp_integral_e(0,x)
        e^(-x)/x

    [AS1964]_ 5.1.24::

        sage: exp_integral_e(6,0)
        1/5
        sage: nn = var('nn')
        sage: assume(nn > 1)
        sage: f = exp_integral_e(nn,0)
        sage: f.simplify()
        1/(nn - 1)


    ALGORITHM:

    Numerical evaluation is handled using mpmath, but symbolics are handled
    by Sage and Maxima.

    """
    def __init__(self):
        """
        See the docstring for :meth:`Function_exp_integral_e`.

        EXAMPLES::

            sage: exp_integral_e(1, 0)
            exp_integral_e(1, 0)
            sage: exp_integral_e(1, x)._sympy_()
            expint(1, x)

        """
        BuiltinFunction.__init__(self, "exp_integral_e", nargs=2,
                                 conversions=dict(maxima='expintegral_e',
                                                  sympy='expint'))

    def _eval_(self, n, z):
        """
        EXAMPLES::

            sage: exp_integral_e(1.0, x)
            exp_integral_e(1.00000000000000, x)
            sage: exp_integral_e(x, 1.0)
            exp_integral_e(x, 1.00000000000000)
            sage: exp_integral_e(3, 0)
            1/2

        TESTS:

        Check that Python ints work (:trac:`14766`)::

            sage: exp_integral_e(int(3), 0)
            1/2
        """
        z_zero = False
        # special case: z == 0 and n > 1
        if isinstance(z, Expression):
            if z.is_trivial_zero():
                z_zero = True # for later
                if n > 1:
                    return 1/(n-1)
        else:
            if not z:
                z_zero = True
                if n > 1:
                    return 1/(n-1)

        # special case: n == 0
        if isinstance(n, Expression):
            if n.is_trivial_zero():
                if z_zero:
                    return None
                else:
                    return exp(-z)/z
        else:
            if not n:
                if z_zero:
                    return None
                else:
                    return exp(-z)/z

        return None # leaves the expression unevaluated

    def _evalf_(self, n, z, parent=None, algorithm=None):
        """
        EXAMPLES::

            sage: exp_integral_e(1.0, 1.0)
            0.219383934395520
            sage: N(exp_integral_e(1, 1+I))
            0.000281624451981418 - 0.179324535039359*I
            sage: exp_integral_e(1, RealField(100)(1))
            0.21938393439552027367716377546

        """
        import mpmath
        return mpmath_utils.call(mpmath.expint, n, z, parent=parent)

    def _print_latex_(self, n, z):
        r"""
        Custom ``_print_latex_`` method.

        EXAMPLES::

            sage: latex(exp_integral_e(1, -x - 1))
            E_{1}\left(-x - 1\right)
        """
        return r"E_{{{}}}\left({}\right)".format(latex(n), latex(z))

    def _derivative_(self, n, z, diff_param=None):
        """
        If `n` is an integer strictly larger than 0, then the derivative of
        `E_n(z)` with respect to `z` is
        `-E_{n-1}(z)`. See [AS1964]_ 5.1.26.

        EXAMPLES::

            sage: x = var('x')
            sage: f = exp_integral_e(2,x)
            sage: f.diff(x)
            -exp_integral_e(1, x)

            sage: f = exp_integral_e(2,sqrt(x))
            sage: f.diff(x)
            -1/2*exp_integral_e(1, sqrt(x))/sqrt(x)

        """
        if n in ZZ and n > 0:
            return -1*exp_integral_e(n-1,z)
        else:
            raise NotImplementedError("The derivative of this function is only implemented for n = 1, 2, 3, ...")

exp_integral_e = Function_exp_integral_e()


class Function_exp_integral_e1(BuiltinFunction):
    r"""
    The generalized complex exponential integral `E_1(z)` defined by

    .. MATH::

        E_1(z) = \int_z^\infty \frac{e^{-t}}{t} \; dt

    see [AS1964]_ 5.1.4.

    EXAMPLES::

        sage: exp_integral_e1(x)
        exp_integral_e1(x)
        sage: exp_integral_e1(1.0)
        0.219383934395520

    Numerical evaluation is handled using mpmath::

        sage: N(exp_integral_e1(1))
        0.219383934395520
        sage: exp_integral_e1(RealField(100)(1))
        0.21938393439552027367716377546

    We can compare this to PARI's evaluation of
    :meth:`exponential_integral_1`::

        sage: N(exp_integral_e1(2.0))
        0.0489005107080611
        sage: N(exponential_integral_1(2.0))
        0.0489005107080611

    Symbolic derivatives and integrals are handled by Sage and Maxima::

        sage: x = var('x')
        sage: f = exp_integral_e1(x)
        sage: f.diff(x)
        -e^(-x)/x

        sage: f.integrate(x)
        -exp_integral_e(2, x)

    ALGORITHM:

    Numerical evaluation is handled using mpmath, but symbolics are handled
    by Sage and Maxima.

    """
    def __init__(self):
        """
        See the docstring for :class:`Function_exp_integral_e1`.

        EXAMPLES::

            sage: exp_integral_e1(1)
            exp_integral_e1(1)
            sage: exp_integral_e1(x)._sympy_()
            expint(1, x)

        """
        BuiltinFunction.__init__(self, "exp_integral_e1", nargs=1,
                                 conversions=dict(maxima='expintegral_e1',
                                                  sympy='E1'))

    def _evalf_(self, z, parent=None, algorithm=None):
        """
        EXAMPLES::

            sage: N(exp_integral_e1(1+I))
            0.000281624451981418 - 0.179324535039359*I
            sage: exp_integral_e1(RealField(200)(0.5))
            0.55977359477616081174679593931508523522684689031635351524829

        """
        import mpmath
        return mpmath_utils_call(mpmath.e1, z, parent=parent)


    def _print_latex_(self, z):
        r"""
        Custom ``_print_latex_`` method.

        EXAMPLES::

            sage: latex(exp_integral_e1(2))
            E_{1}\left(2\right)
        """
        return r"E_{{1}}\left({}\right)".format(latex(z))

    def _derivative_(self, z, diff_param=None):
        """
        The derivative of `E_1(z)` is `-e^{-z}/z`.

        See [AS1964]_ 5.1.26.

        EXAMPLES::

            sage: x = var('x')
            sage: f = exp_integral_e1(x)
            sage: f.diff(x)
            -e^(-x)/x

            sage: f = exp_integral_e1(x^2)
            sage: f.diff(x)
            -2*e^(-x^2)/x

        """
        return -exp(-z)/z

exp_integral_e1 = Function_exp_integral_e1()


class Function_log_integral(BuiltinFunction):
    r"""
    The logarithmic integral `\operatorname{li}(z)` defined by

    .. MATH::

        \operatorname{li}(x) = \int_0^z \frac{dt}{\ln(t)} = \operatorname{Ei}(\ln(x))

    for x > 1 and by analytic continuation for complex arguments z (see [AS1964]_ 5.1.3).

    EXAMPLES:

    Numerical evaluation for real and complex arguments is handled using mpmath::

        sage: N(log_integral(3))
        2.16358859466719
        sage: N(log_integral(3), digits=30)
        2.16358859466719197287692236735
        sage: log_integral(ComplexField(100)(3+I))
        2.2879892769816826157078450911 + 0.87232935488528370139883806779*I
        sage: log_integral(0)
        0

    Symbolic derivatives and integrals are handled by Sage and Maxima::

        sage: x = var('x')
        sage: f = log_integral(x)
        sage: f.diff(x)
        1/log(x)

        sage: f.integrate(x)
        x*log_integral(x) - Ei(2*log(x))

    Here is a test from the mpmath documentation. There are
    1,925,320,391,606,803,968,923 many prime numbers less than 1e23. The
    value of ``log_integral(1e23)`` is very close to this::

        sage: log_integral(1e23)
        1.92532039161405e21

    ALGORITHM:

    Numerical evaluation is handled using mpmath, but symbolics are handled
    by Sage and Maxima.

    REFERENCES:

    - :wikipedia:`Logarithmic_integral_function`
    - mpmath documentation: `logarithmic-integral`_

    .. _`logarithmic-integral`: http://mpmath.org/doc/current/functions/expintegrals.html#logarithmic-integral


    """
    def __init__(self):
        r"""
        See the docstring for ``Function_log_integral``.

        EXAMPLES::

            sage: log_integral(3)
            log_integral(3)
            sage: log_integral(x)._sympy_()
            li(x)
            sage: log_integral(x)._fricas_init_()
            'li(x)'

        TESTS:

        Verify that :trac:`28917` is fixed::

            sage: latex(log_integral(x))
            \operatorname{log\_integral}\left(x\right)
        """
        BuiltinFunction.__init__(self, "log_integral", nargs=1,
                                 latex_name=r'\operatorname{log\_integral}',
                                 conversions=dict(maxima='expintegral_li',
                                                  sympy='li',
                                                  fricas='li'))

    def _eval_(self, z):
        """
        EXAMPLES::

            sage: z = var('z')
            sage: log_integral(z)
            log_integral(z)
            sage: log_integral(3.0)
            2.16358859466719
            sage: log_integral(0)
            0

        """
        # Special case z = 0
        if isinstance(z, Expression):
            if z.is_trivial_zero():
                return z
        elif not z:
            return z

    def _evalf_(self, z, parent=None, algorithm=None):
        """
        EXAMPLES::

            sage: N(log_integral(1e6))
            78627.5491594622
            sage: log_integral(RealField(200)(1e6))
            78627.549159462181919862910747947261161321874382421767074759

        """
        import mpmath
        return mpmath_utils_call(mpmath.li, z, parent=parent)

    def _derivative_(self, z, diff_param=None):
        r"""
        The derivative of `\operatorname{li}(z) is `1/log(z)`.

        EXAMPLES::

            sage: x = var('x')
            sage: f = log_integral(x)
            sage: f.diff(x)
            1/log(x)

            sage: f = log_integral(x^2)
            sage: f.diff(x)
            2*x/log(x^2)

        """
        return 1/log(z)

li = log_integral = Function_log_integral()

class Function_log_integral_offset(BuiltinFunction):
    r"""
    The offset logarithmic integral, or Eulerian logarithmic integral,
    `\operatorname{Li}(x)` is defined by

    .. MATH::

        \operatorname{Li}(x) = \int_2^x \frac{dt}{\ln(t)} =
        \operatorname{li}(x)-\operatorname{li}(2)

    for `x \ge 2`.

    The offset logarithmic integral should also not be confused with the
    polylogarithm (also denoted by `\operatorname{Li}(x)` ), which is
    implemented as :class:`sage.functions.log.Function_polylog`.

    `\operatorname{Li}(x)` is identical to `\operatorname{li}(x)` except that
    the lower limit of integration is `2` rather than `0` to avoid the
    singularity at `x = 1` of

    .. MATH::

        \frac{1}{\ln(t)}

    See :class:`Function_log_integral` for details of `\operatorname{li}(x)`.
    Thus `\operatorname{Li}(x)` can also be represented by

    .. MATH::

        \operatorname{Li}(x) = \operatorname{li}(x)-\operatorname{li}(2)

    So we have::

        sage: li(4.5)-li(2.0)-Li(4.5)
        0.000000000000000

    `\operatorname{Li}(x)` is extended to complex arguments `z`
    by analytic continuation (see [AS1964]_ 5.1.3)::

        sage: Li(6.6+5.4*I)
        3.97032201503632 + 2.62311237593572*I

    The function `\operatorname{Li}` is an approximation for the number of
    primes up to `x`. In fact, the famous Riemann Hypothesis is

    .. MATH::

        |\pi(x) - \operatorname{Li}(x)| \leq \sqrt{x} \log(x).

    For "small" `x`, `\operatorname{Li}(x)` is always slightly bigger
    than `\pi(x)`. However it is a theorem that there are very
    large values of `x` (e.g., around `10^{316}`), such that
    `\exists x: \pi(x) > \operatorname{Li}(x)`.  See "A new bound for the
    smallest x with `\pi(x) > \operatorname{li}(x)`",
    Bays and Hudson, Mathematics of Computation, 69 (2000) 1285-1296.

    .. NOTE::

        Definite integration returns a part symbolic and part
        numerical result.  This is because when Li(x) is evaluated it is
        passed as li(x)-li(2).

    EXAMPLES:

    Numerical evaluation for real and complex arguments is handled using mpmath::

        sage: N(log_integral_offset(3))
        1.11842481454970
        sage: N(log_integral_offset(3), digits=30)
        1.11842481454969918803233347815
        sage: log_integral_offset(ComplexField(100)(3+I))
        1.2428254968641898308632562019 + 0.87232935488528370139883806779*I
        sage: log_integral_offset(2)
        0
        sage: for n in range(1,7):
        ....:  print('%-10s%-10s%-20s'%(10^n, prime_pi(10^n), N(Li(10^n))))
        10        4         5.12043572466980
        100       25        29.0809778039621
        1000      168       176.564494210035
        10000     1229      1245.09205211927
        100000    9592      9628.76383727068
        1000000   78498     78626.5039956821

    Here is a test from the mpmath documentation.
    There are 1,925,320,391,606,803,968,923 prime numbers less than 1e23.
    The value of ``log_integral_offset(1e23)`` is very close to this::

        sage: log_integral_offset(1e23)
        1.92532039161405e21

    Symbolic derivatives are handled by Sage and integration by Maxima::

        sage: x = var('x')
        sage: f = log_integral_offset(x)
        sage: f.diff(x)
        1/log(x)
        sage: f.integrate(x)
        -x*log_integral(2) + x*log_integral(x) - Ei(2*log(x))
        sage: Li(x).integrate(x,2.0,4.5)
        -2.5*log_integral(2) + 5.799321147411334
        sage: N(f.integrate(x,2.0,3.0))   # abs tol 1e-15
        0.601621785860587

    ALGORITHM:

    Numerical evaluation is handled using mpmath, but symbolics are handled
    by Sage and Maxima.

    REFERENCES:

    - :wikipedia:`Logarithmic_integral_function`
    - mpmath documentation: `logarithmic-integral`_

    .. _`logarithmic-integral`: http://mpmath.org/doc/current/functions/expintegrals.html#logarithmic-integral
    """

    def __init__(self):
        r"""
        See the docstring for ``Function_log_integral_offset``.

        EXAMPLES::

            sage: log_integral_offset(3)
            log_integral(3) - log_integral(2)
            sage: log_integral_offset(x, hold=True)._sympy_()
            Li(x)

        TESTS:

        Verify that the problem described in :trac:`28917` no longer appears here::

            sage: latex(log_integral_offset)
            \operatorname{log\_integral\_offset}

        """
        BuiltinFunction.__init__(self, "log_integral_offset", nargs=1,
                                 latex_name=r'\operatorname{log\_integral\_offset}',
                                 conversions=dict(sympy='Li'))

    def _eval_(self,z):
        """
        EXAMPLES::

            sage: z = var('z')
            sage: log_integral_offset(z)
            -log_integral(2) + log_integral(z)
            sage: log_integral_offset(3.0)
            1.11842481454970
            sage: log_integral_offset(2)
            0

        """
        if z == 2:
            import sage.symbolic.ring
            return sage.symbolic.ring.SR(0)
        return li(z)-li(2)
        # If we return:(li(z)-li(2)) we get correct symbolic integration.
        # But on definite integration it returns x.xxxx-li(2).

    def _evalf_(self, z, parent=None, algorithm=None):
        """
        EXAMPLES::

            sage: N(log_integral_offset(1e6))
            78626.5039956821
            sage: log_integral_offset(RealField(200)(1e6))
            78626.503995682064427078066159058066548185351766843615873183
            sage: li(4.5)-li(2.0)-Li(4.5)
            0.000000000000000

        """
        import mpmath
        return mpmath_utils_call(mpmath.li, z, offset=True, parent=parent)

    def _derivative_(self, z, diff_param=None):
        r"""
        The derivative of `\operatorname{Li}(z) is `1/log(z)`.

        EXAMPLES::

            sage: x = var('x')
            sage: f = log_integral_offset(x)
            sage: f.diff(x)
            1/log(x)

            sage: f = log_integral_offset(x^2)
            sage: f.diff(x)
            2*x/log(x^2)
        """
        return 1/log(z)

Li = log_integral_offset = Function_log_integral_offset()

class Function_sin_integral(BuiltinFunction):
    r"""
    The trigonometric integral `\operatorname{Si}(z)` defined by

    .. MATH::

        \operatorname{Si}(z) = \int_0^z \frac{\sin(t)}{t} \; dt,

    see [AS1964]_ 5.2.1.

    EXAMPLES:

    Numerical evaluation for real and complex arguments is handled using mpmath::

        sage: sin_integral(0)
        0
        sage: sin_integral(0.0)
        0.000000000000000
        sage: sin_integral(3.0)
        1.84865252799947
        sage: N(sin_integral(3), digits=30)
        1.84865252799946825639773025111
        sage: sin_integral(ComplexField(100)(3+I))
        2.0277151656451253616038525998 + 0.015210926166954211913653130271*I

    The alias ``Si`` can be used instead of ``sin_integral``::

        sage: Si(3.0)
        1.84865252799947

    The limit of `\operatorname{Si}(z)` as `z \to \infty` is `\pi/2`::

        sage: N(sin_integral(1e23))
        1.57079632679490
        sage: N(pi/2)
        1.57079632679490

    At 200 bits of precision `\operatorname{Si}(10^{23})` agrees with `\pi/2` up to
    `10^{-24}`::

        sage: sin_integral(RealField(200)(1e23))
        1.5707963267948966192313288218697837425815368604836679189519
        sage: N(pi/2, prec=200)
        1.5707963267948966192313216916397514420985846996875529104875

    The exponential sine integral is analytic everywhere::

        sage: sin_integral(-1.0)
        -0.946083070367183
        sage: sin_integral(-2.0)
        -1.60541297680269
        sage: sin_integral(-1e23)
        -1.57079632679490

    Symbolic derivatives and integrals are handled by Sage and Maxima::

        sage: x = var('x')
        sage: f = sin_integral(x)
        sage: f.diff(x)
        sin(x)/x

        sage: f.integrate(x)
        x*sin_integral(x) + cos(x)

        sage: integrate(sin(x)/x, x)
        -1/2*I*Ei(I*x) + 1/2*I*Ei(-I*x)


    Compare values of the functions `\operatorname{Si}(x)` and
    `f(x) = (1/2)i \cdot \operatorname{Ei}(-ix) - (1/2)i \cdot
    \operatorname{Ei}(ix) - \pi/2`, which are both anti-derivatives of
    `\sin(x)/x`, at some random positive real numbers::

        sage: f(x) = 1/2*I*Ei(-I*x) - 1/2*I*Ei(I*x) - pi/2
        sage: g(x) = sin_integral(x)
        sage: R = [ abs(RDF.random_element()) for i in range(100) ]
        sage: all(abs(f(x) - g(x)) < 1e-10 for x in R)
        True

    The Nielsen spiral is the parametric plot of (Si(t), Ci(t))::

        sage: x=var('x')
        sage: f(x) = sin_integral(x)
        sage: g(x) = cos_integral(x)
        sage: P = parametric_plot([f, g], (x, 0.5 ,20))
        sage: show(P, frame=True, axes=False)

    ALGORITHM:

    Numerical evaluation is handled using mpmath, but symbolics are handled
    by Sage and Maxima.

    REFERENCES:

    - :wikipedia:`Trigonometric_integral`
    - mpmath documentation: `si`_

    .. _`si`: http://mpmath.org/doc/current/functions/expintegrals.html#si

    """
    def __init__(self):
        """
        See the docstring for ``Function_sin_integral``.

        EXAMPLES::

            sage: sin_integral(1)
            sin_integral(1)
            sage: sin_integral(x)._sympy_()
            Si(x)
            sage: sin_integral(x)._fricas_init_()
            'Si(x)'
            sage: sin_integral(x)._giac_()
            Si(sageVARx)
        """
        BuiltinFunction.__init__(self, "sin_integral", nargs=1,
                                 latex_name=r'\operatorname{Si}',
                                 conversions=dict(maxima='expintegral_si',
                                                  sympy='Si',
                                                  fricas='Si', giac='Si'))

    def _eval_(self, z):
        """
        EXAMPLES::

            sage: z = var('z')
            sage: sin_integral(z)
            sin_integral(z)
            sage: sin_integral(3.0)
            1.84865252799947
            sage: sin_integral(0)
            0

        """
        if isinstance(z, Expression):
            if z.is_trivial_zero():
                return z
        elif not z:
            return z

    def _evalf_(self, z, parent=None, algorithm=None):
        r"""
        EXAMPLES:

        The limit `\operatorname{Si}(z)` as `z \to \infty`  is `\pi/2`::

            sage: N(sin_integral(1e23) - pi/2)
            0.000000000000000

        At 200 bits of precision `\operatorname{Si}(10^{23})` agrees with `\pi/2` up to
        `10^{-24}`::

            sage: sin_integral(RealField(200)(1e23))
            1.5707963267948966192313288218697837425815368604836679189519
            sage: N(pi/2, prec=200)
            1.5707963267948966192313216916397514420985846996875529104875

        The exponential sine integral is analytic everywhere, even on the
        negative real axis::

            sage: sin_integral(-1.0)
            -0.946083070367183
            sage: sin_integral(-2.0)
            -1.60541297680269
            sage: sin_integral(-1e23)
            -1.57079632679490
        """
        import mpmath
        return mpmath_utils_call(mpmath.si, z, parent=parent)

    def _derivative_(self, z, diff_param=None):
        r"""
        The derivative of `\operatorname{Si}(z)` is `\sin(z)/z` if `z`
        is not zero. The derivative at `z = 0` is `1` (but this
        exception is not currently implemented).

        EXAMPLES::

            sage: x = var('x')
            sage: f = sin_integral(x)
            sage: f.diff(x)
            sin(x)/x

            sage: f = sin_integral(x^2)
            sage: f.diff(x)
            2*sin(x^2)/x

        """
        return sin(z)/z

Si = sin_integral = Function_sin_integral()


class Function_cos_integral(BuiltinFunction):
    r"""
    The trigonometric integral `\operatorname{Ci}(z)` defined by

    .. MATH::

        \operatorname{Ci}(z) = \gamma + \log(z) + \int_0^z \frac{\cos(t)-1}{t} \; dt,

    where `\gamma` is the Euler gamma constant (``euler_gamma`` in Sage),
    see [AS1964]_ 5.2.1.

    EXAMPLES::

        sage: z = var('z')
        sage: cos_integral(z)
        cos_integral(z)
        sage: cos_integral(3.0)
        0.119629786008000
        sage: cos_integral(0)
        cos_integral(0)
        sage: N(cos_integral(0))
        -infinity

    Numerical evaluation for real and complex arguments is handled using mpmath::

        sage: cos_integral(3.0)
        0.119629786008000

    The alias ``Ci`` can be used instead of ``cos_integral``::

        sage: Ci(3.0)
        0.119629786008000

    Compare ``cos_integral(3.0)`` to the definition of the value using
    numerical integration::

        sage: a = numerical_integral((cos(x)-1)/x, 0, 3)[0]
        sage: abs(N(euler_gamma + log(3)) + a - N(cos_integral(3.0))) < 1e-14
        True

    Arbitrary precision and complex arguments are handled::

        sage: N(cos_integral(3), digits=30)
        0.119629786008000327626472281177
        sage: cos_integral(ComplexField(100)(3+I))
        0.078134230477495714401983633057 - 0.37814733904787920181190368789*I

    The limit `\operatorname{Ci}(z)` as `z \to \infty` is zero::

        sage: N(cos_integral(1e23))
        -3.24053937643003e-24

    Symbolic derivatives and integrals are handled by Sage and Maxima::

        sage: x = var('x')
        sage: f = cos_integral(x)
        sage: f.diff(x)
        cos(x)/x

        sage: f.integrate(x)
        x*cos_integral(x) - sin(x)

    The Nielsen spiral is the parametric plot of (Si(t), Ci(t))::

        sage: t=var('t')
        sage: f(t) = sin_integral(t)
        sage: g(t) = cos_integral(t)
        sage: P = parametric_plot([f, g], (t, 0.5 ,20))
        sage: show(P, frame=True, axes=False)

    ALGORITHM:

    Numerical evaluation is handled using mpmath, but symbolics are handled
    by Sage and Maxima.

    REFERENCES:

    - :wikipedia:`Trigonometric_integral`
    - mpmath documentation: `ci`_

    .. _`ci`: http://mpmath.org/doc/current/functions/expintegrals.html#ci

    """
    def __init__(self):
        """
        See the docstring for :class:`Function_cos_integral`.

        EXAMPLES::

            sage: cos_integral(1)
            cos_integral(1)
            sage: cos_integral(x)._sympy_()
            Ci(x)
            sage: cos_integral(x)._fricas_init_()
            'Ci(x)'
            sage: cos_integral(x)._giac_()
            Ci(sageVARx)
        """
        BuiltinFunction.__init__(self, "cos_integral", nargs=1,
                                 latex_name=r'\operatorname{Ci}',
                                 conversions=dict(maxima='expintegral_ci',
                                                  sympy='Ci',
                                                  fricas='Ci', giac='Ci'))

    def _evalf_(self, z, parent=None, algorithm=None):
        """
        EXAMPLES::

            sage: N(cos_integral(1e23)) < 1e-20
            True
            sage: N(cos_integral(10^-10), digits=30)
            -22.4486352650389239795759024568
            sage: cos_integral(ComplexField(100)(I))
            0.83786694098020824089467857943 + 1.5707963267948966192313216916*I

        """
        import mpmath
        return mpmath_utils_call(mpmath.ci, z, parent=parent)

    def _derivative_(self, z, diff_param=None):
        r"""
        The derivative of `\operatorname{Ci}(z)` is `\cos(z)/z` if `z` is not zero.

        EXAMPLES::

            sage: x = var('x')
            sage: f = cos_integral(x)
            sage: f.diff(x)
            cos(x)/x

            sage: f = cos_integral(x^2)
            sage: f.diff(x)
            2*cos(x^2)/x

        """
        return cos(z)/z

Ci = cos_integral = Function_cos_integral()


class Function_sinh_integral(BuiltinFunction):
    r"""
    The trigonometric integral `\operatorname{Shi}(z)` defined by

    .. MATH::

        \operatorname{Shi}(z) = \int_0^z \frac{\sinh(t)}{t} \; dt,

    see [AS1964]_ 5.2.3.

    EXAMPLES:

    Numerical evaluation for real and complex arguments is handled using mpmath::

        sage: sinh_integral(3.0)
        4.97344047585981
        sage: sinh_integral(1.0)
        1.05725087537573
        sage: sinh_integral(-1.0)
        -1.05725087537573

    The alias ``Shi`` can be used instead of ``sinh_integral``::

        sage: Shi(3.0)
        4.97344047585981

    Compare ``sinh_integral(3.0)`` to the definition of the value using
    numerical integration::

        sage: a = numerical_integral(sinh(x)/x, 0, 3)[0]
        sage: abs(a - N(sinh_integral(3))) < 1e-14
        True

    Arbitrary precision and complex arguments are handled::

        sage: N(sinh_integral(3), digits=30)
        4.97344047585980679771041838252
        sage: sinh_integral(ComplexField(100)(3+I))
        3.9134623660329374406788354078 + 3.0427678212908839256360163759*I

    The limit `\operatorname{Shi}(z)` as `z \to \infty` is `\infty`::

        sage: N(sinh_integral(Infinity))
        +infinity

    Symbolic derivatives and integrals are handled by Sage and Maxima::

        sage: x = var('x')
        sage: f = sinh_integral(x)
        sage: f.diff(x)
        sinh(x)/x

        sage: f.integrate(x)
        x*sinh_integral(x) - cosh(x)

    Note that due to some problems with the way Maxima handles these
    expressions, definite integrals can sometimes give unexpected
    results (typically when using inexact endpoints) due to
    inconsistent branching::

        sage: integrate(sinh_integral(x), x, 0, 1/2)
        -cosh(1/2) + 1/2*sinh_integral(1/2) + 1
        sage: integrate(sinh_integral(x), x, 0, 1/2).n() # correct
        0.125872409703453
        sage: integrate(sinh_integral(x), x, 0, 0.5).n() # fixed in maxima 5.29.1
        0.125872409703453

    ALGORITHM:

    Numerical evaluation is handled using mpmath, but symbolics are handled
    by Sage and Maxima.

    REFERENCES:

    - :wikipedia:`Trigonometric_integral`
    - mpmath documentation: `shi`_

    .. _`shi`: http://mpmath.org/doc/current/functions/expintegrals.html#shi

    """
    def __init__(self):
        """
        See the docstring for ``Function_sinh_integral``.

        EXAMPLES::

            sage: sinh_integral(1)
            sinh_integral(1)
            sage: sinh_integral(x)._sympy_()
            Shi(x)

        """
        BuiltinFunction.__init__(self, "sinh_integral", nargs=1,
                                 latex_name=r'\operatorname{Shi}',
                                 conversions=dict(maxima='expintegral_shi',
                                                  sympy='Shi',
                                                  fricas='Shi'))

    def _eval_(self, z):
        """
        EXAMPLES::

            sage: z = var('z')
            sage: sinh_integral(z)
            sinh_integral(z)
            sage: sinh_integral(3.0)
            4.97344047585981
            sage: sinh_integral(0)
            0

        """
        # special case: z = 0
        if isinstance(z, Expression):
            if z.is_trivial_zero():
                return z
        elif not z:
            return z

    def _evalf_(self, z, parent=None, algorithm=None):
        """
        EXAMPLES::

            sage: N(sinh_integral(10^-10), digits=30)
            1.00000000000000000000055555556e-10
            sage: sinh_integral(ComplexField(100)(I))
            0.94608307036718301494135331382*I

        """
        import mpmath
        return mpmath_utils_call(mpmath.shi, z, parent=parent)

    def _derivative_(self, z, diff_param=None):
        r"""
        The derivative of `\operatorname{Shi}(z)` is `\sinh(z)/z`.

        EXAMPLES::

            sage: x = var('x')
            sage: f = sinh_integral(x)
            sage: f.diff(x)
            sinh(x)/x

            sage: f = sinh_integral(ln(x))
            sage: f.diff(x)
            1/2*(x^2 - 1)/(x^2*log(x))

        """
        return sinh(z)/z

Shi = sinh_integral = Function_sinh_integral()


class Function_cosh_integral(BuiltinFunction):
    r"""
    The trigonometric integral `\operatorname{Chi}(z)` defined by

    .. MATH::

        \operatorname{Chi}(z) = \gamma + \log(z) + \int_0^z \frac{\cosh(t)-1}{t} \; dt,

    see [AS1964]_ 5.2.4.

    EXAMPLES::

        sage: z = var('z')
        sage: cosh_integral(z)
        cosh_integral(z)
        sage: cosh_integral(3.0)
        4.96039209476561

    Numerical evaluation for real and complex arguments is handled using mpmath::

        sage: cosh_integral(1.0)
        0.837866940980208

    The alias ``Chi`` can be used instead of ``cosh_integral``::

        sage: Chi(1.0)
        0.837866940980208

    Here is an example from the mpmath documentation::

        sage: f(x) = cosh_integral(x)
        sage: find_root(f, 0.1, 1.0)
        0.523822571389...

    Compare ``cosh_integral(3.0)`` to the definition of the value using
    numerical integration::

        sage: a = numerical_integral((cosh(x)-1)/x, 0, 3)[0]
        sage: abs(N(euler_gamma + log(3)) + a - N(cosh_integral(3.0))) < 1e-14
        True

    Arbitrary precision and complex arguments are handled::

        sage: N(cosh_integral(3), digits=30)
        4.96039209476560976029791763669
        sage: cosh_integral(ComplexField(100)(3+I))
        3.9096723099686417127843516794 + 3.0547519627014217273323873274*I

    The limit of `\operatorname{Chi}(z)` as `z \to \infty` is `\infty`::

        sage: N(cosh_integral(Infinity))
        +infinity

    Symbolic derivatives and integrals are handled by Sage and Maxima::

        sage: x = var('x')
        sage: f = cosh_integral(x)
        sage: f.diff(x)
        cosh(x)/x

        sage: f.integrate(x)
        x*cosh_integral(x) - sinh(x)

    ALGORITHM:

    Numerical evaluation is handled using mpmath, but symbolics are handled
    by Sage and Maxima.

    REFERENCES:

    - :wikipedia:`Trigonometric_integral`
    - mpmath documentation: `chi`_

    .. _`chi`: http://mpmath.org/doc/current/functions/expintegrals.html#chi

    """
    def __init__(self):
        """
        See the docstring for ``Function_cosh_integral``.

        EXAMPLES::

            sage: cosh_integral(1)
            cosh_integral(1)
            sage: cosh_integral(x)._sympy_()
            Chi(x)

        """
        BuiltinFunction.__init__(self, "cosh_integral", nargs=1,
                                 latex_name=r'\operatorname{Chi}',
                                 conversions=dict(maxima='expintegral_chi',
                                                  sympy='Chi',
                                                  fricas='Chi'))

    def _evalf_(self, z, parent=None, algorithm=None):
        """
        EXAMPLES::

            sage: N(cosh_integral(10^-10), digits=30)
            -22.4486352650389239795709024568
            sage: cosh_integral(ComplexField(100)(I))
            0.33740392290096813466264620389 + 1.5707963267948966192313216916*I

        """
        import mpmath
        return mpmath_utils_call(mpmath.chi, z, parent=parent)

    def _derivative_(self, z, diff_param=None):
        r"""
        The derivative of `\operatorname{Chi}(z)` is `\cosh(z)/z`.

        EXAMPLES::

            sage: x = var('x')
            sage: f = cosh_integral(x)
            sage: f.diff(x)
            cosh(x)/x

            sage: f = cosh_integral(ln(x))
            sage: f.diff(x)
            1/2*(x^2 + 1)/(x^2*log(x))

        """
        return cosh(z)/z

Chi = cosh_integral = Function_cosh_integral()


###################################################################
# Code below here was moved from sage/functions/transcendental.py
# This occurred as part of Trac #11143.
###################################################################

# This class has a name which is not specific enough
# see Function_exp_integral_e above, for example, which
# is the "generalized" exponential integral function. We
# are leaving the name the same for backwards compatibility
# purposes.
class Function_exp_integral(BuiltinFunction):
    r"""
    The generalized complex exponential integral Ei(z) defined by

    .. MATH::

        \operatorname{Ei}(x) = \int_{-\infty}^x \frac{e^t}{t} \; dt

    for x > 0 and for complex arguments by analytic continuation,
    see [AS1964]_ 5.1.2.

    EXAMPLES::

        sage: Ei(10)
        Ei(10)
        sage: Ei(I)
        Ei(I)
        sage: Ei(3+I)
        Ei(I + 3)
        sage: Ei(1.3)
        2.72139888023202
        sage: Ei(10r)
        Ei(10)
        sage: Ei(1.3r)
        2.7213988802320235

    The branch cut for this function is along the negative real axis::

        sage: Ei(-3 + 0.1*I)
        -0.0129379427181693 + 3.13993830250942*I
        sage: Ei(-3 - 0.1*I)
        -0.0129379427181693 - 3.13993830250942*I

    The precision for the result is deduced from the precision of the
    input. Convert the input to a higher precision explicitly if a
    result with higher precision is desired::

        sage: Ei(RealField(300)(1.1))
        2.16737827956340282358378734233807621497112737591639704719499002090327541763352339357795426

    ALGORITHM: Uses mpmath.

    TESTS:

    Show that the evaluation and limit issue in :trac:`13271` is fixed::

        sage: var('Z')
        Z
        sage: (Ei(-Z)).limit(Z=oo)
        0
        sage: (Ei(-Z)).limit(Z=1000)
        Ei(-1000)
        sage: (Ei(-Z)).limit(Z=1000).n()
        -5.07089306023517e-438
    """
    def __init__(self):
        """
        TESTS::

            sage: Ei(10)
            Ei(10)
            sage: Ei(x)._sympy_()
            Ei(x)
        """
        BuiltinFunction.__init__(self, "Ei",
                                 conversions=dict(maxima='expintegral_ei',
                                                  sympy='Ei',
                                                  fricas='Ei'))

    def _evalf_(self, x, parent=None, algorithm=None):
        """
        EXAMPLES::

            sage: Ei(10).n()
            2492.22897624188
            sage: Ei(20).n()
            2.56156526640566e7
            sage: Ei(I).n()
            0.337403922900968 + 2.51687939716208*I
            sage: Ei(3+I).n()
            7.82313467600158 + 6.09751978399231*I
        """
        import mpmath
        return mpmath_utils_call(mpmath.ei, x, parent=parent)

    def _derivative_(self, x, diff_param=None):
        """
        EXAMPLES::

            sage: Ei(x).diff(x)
            e^x/x
            sage: Ei(x).diff(x).subs(x=1)
            e
            sage: Ei(x^2).diff(x)
            2*e^(x^2)/x
            sage: f = function('f')
            sage: Ei(f(x)).diff(x)
            e^f(x)*diff(f(x), x)/f(x)
        """
        return exp(x)/x

Ei = exp_integral_ei = Function_exp_integral()


# moved here from sage/functions/transcendental.py
def exponential_integral_1(x, n=0):
    r"""
    Returns the exponential integral `E_1(x)`. If the optional
    argument `n` is given, computes list of the first
    `n` values of the exponential integral
    `E_1(x m)`.

    The exponential integral `E_1(x)` is

    .. MATH::

                      E_1(x) = \int_{x}^{\infty} \frac{e^{-t}}{t} \; dt

    INPUT:

    - ``x`` -- a positive real number

    - ``n`` -- (default: 0) a nonnegative integer; if
      nonzero, then return a list of values ``E_1(x*m)`` for m =
      1,2,3,...,n. This is useful, e.g., when computing derivatives of
      L-functions.


    OUTPUT:

    A real number if n is 0 (the default) or a list of reals if n > 0.
    The precision is the same as the input, with a default of 53 bits
    in case the input is exact.

    EXAMPLES::

        sage: exponential_integral_1(2)
        0.0489005107080611
        sage: exponential_integral_1(2, 4)  # abs tol 1e-18
        [0.0489005107080611, 0.00377935240984891, 0.000360082452162659, 0.0000376656228439245]
        sage: exponential_integral_1(40, 5)
        [0.000000000000000, 2.22854325868847e-37, 6.33732515501151e-55, 2.02336191509997e-72, 6.88522610630764e-90]
        sage: exponential_integral_1(0)
        +Infinity
        sage: r = exponential_integral_1(RealField(150)(1))
        sage: r
        0.21938393439552027367716377546012164903104729
        sage: parent(r)
        Real Field with 150 bits of precision
        sage: exponential_integral_1(RealField(150)(100))
        3.6835977616820321802351926205081189876552201e-46

    TESTS:

    The relative error for a single value should be less than 1 ulp::

        sage: for prec in [20..1000]:  # long time (22s on sage.math, 2013)
        ....:     R = RealField(prec)
        ....:     S = RealField(prec+64)
        ....:     for t in range(8):  # Try 8 values for each precision
        ....:         a = R.random_element(-15,10).exp()
        ....:         x = exponential_integral_1(a)
        ....:         y = exponential_integral_1(S(a))
        ....:         e = float(abs(S(x) - y)/x.ulp())
        ....:         if e >= 1.0:
        ....:             print("exponential_integral_1(%s) with precision %s has error of %s ulp"%(a, prec, e))

    The absolute error for a vector should be less than `2^{-p} c`, where
    `p` is the precision in bits of `x` and `c = 2` ``max(1, exponential_integral_1(x))``::

        sage: for prec in [20..128]:  # long time (15s on sage.math, 2013)
        ....:     R = RealField(prec)
        ....:     S = RealField(prec+64)
        ....:     a = R.random_element(-15,10).exp()
        ....:     n = 2^ZZ.random_element(14)
        ....:     x = exponential_integral_1(a, n)
        ....:     y = exponential_integral_1(S(a), n)
        ....:     c = RDF(4 * max(1.0, y[0]))
        ....:     for i in range(n):
        ....:         e = float(abs(S(x[i]) - y[i]) << prec)
        ....:         if e >= c:
        ....:             print("exponential_integral_1(%s, %s)[%s] with precision %s has error of %s >= %s"%(a, n, i, prec, e, c))

    ALGORITHM: use the PARI C-library function ``eint1``.

    REFERENCE:

    - See Proposition 5.6.12 of Cohen's book "A Course in
      Computational Algebraic Number Theory".
    """
    if isinstance(x, Expression):
        if x.is_trivial_zero():
            from sage.rings.infinity import Infinity
            return Infinity
        else:
            raise NotImplementedError("Use the symbolic exponential integral " +
                                      "function: exp_integral_e1.")

    # x == 0  =>  return Infinity
    if not x:
        from sage.rings.infinity import Infinity
        return Infinity

    # Figure out output precision
    try:
        prec = parent(x).precision()
    except AttributeError:
        prec = 53

    R = RealField(prec)
    if n <= 0:
        # Add extra bits to the input.
        # (experimentally verified -- Jeroen Demeyer)
        inprec = prec + 5 + math.ceil(math.log(prec))
        x = RealField(inprec)(x).__pari__()
        return R(x.eint1())
    else:
        # PARI's algorithm is less precise as n grows larger:
        # add extra bits.
        # (experimentally verified -- Jeroen Demeyer)
        inprec = prec + 1 + math.ceil(1.4427 * math.log(n))
        x = RealField(inprec)(x).__pari__()
        return [R(z) for z in x.eint1(n)]
