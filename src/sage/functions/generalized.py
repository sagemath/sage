r"""
Generalized functions

Sage implements several generalized functions (also known as
distributions) such as Dirac delta, Heaviside step functions. These
generalized functions can be manipulated within Sage like any other
symbolic functions.


AUTHORS:

- Golam Mortuza Hossain (2009-06-26): initial version

EXAMPLES:

Dirac delta function::

    sage: dirac_delta(x)                                                                # needs sage.symbolic
    dirac_delta(x)

Heaviside step function::

    sage: heaviside(x)                                                                  # needs sage.symbolic
    heaviside(x)

Unit step function::

    sage: unit_step(x)                                                                  # needs sage.symbolic
    unit_step(x)

Signum (sgn) function::

    sage: sgn(x)                                                                        # needs sage.symbolic
    sgn(x)

Kronecker delta function::

    sage: m, n = var('m,n')                                                             # needs sage.symbolic
    sage: kronecker_delta(m, n)                                                         # needs sage.symbolic
    kronecker_delta(m, n)
"""

##############################################################################
#
#       Copyright (C) 2009 Golam Mortuza Hossain <gmhossain@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL v2+)
#                  http://www.gnu.org/licenses/
#
##############################################################################

from sage.misc.lazy_import import lazy_import
from sage.rings.integer_ring import ZZ
from sage.symbolic.function import BuiltinFunction, GinacFunction

lazy_import('sage.misc.latex', 'latex')
lazy_import('sage.rings.complex_interval_field', 'ComplexIntervalField')


class FunctionDiracDelta(BuiltinFunction):
    r"""
    The Dirac delta (generalized) function, `\delta(x)` (``dirac_delta(x)``).

    INPUT:

    - ``x`` -- a real number or a symbolic expression

    DEFINITION:

    Dirac delta function `\delta(x)`, is defined in Sage as:

        `\delta(x) = 0` for real `x \ne 0` and
        `\int_{-\infty}^{\infty} \delta(x) dx = 1`

    Its alternate definition with respect to an arbitrary test
    function `f(x)` is

        `\int_{-\infty}^{\infty} f(x) \delta(x-a) dx = f(a)`

    EXAMPLES::

        sage: # needs sage.symbolic
        sage: dirac_delta(1)
        0
        sage: dirac_delta(0)
        dirac_delta(0)
        sage: dirac_delta(x)
        dirac_delta(x)
        sage: integrate(dirac_delta(x), x, -1, 1, algorithm='sympy')                    # needs sympy
        1

    REFERENCES:

    - :wikipedia:`Dirac_delta_function`
    """
    def __init__(self):
        r"""
        The Dirac delta (generalized) function, ``dirac_delta(x)``.

        INPUT:

        - ``x`` -- a real number or a symbolic expression

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: dirac_delta(1)
            0
            sage: dirac_delta(0)
            dirac_delta(0)
            sage: dirac_delta(x)
            dirac_delta(x)
            sage: latex(dirac_delta(x))
            \delta\left(x\right)

            sage: loads(dumps(dirac_delta(x)))                                          # needs sage.symbolic
            dirac_delta(x)
            sage: dirac_delta(x)._sympy_()                                              # needs sympy sage.symbolic
            DiracDelta(x)
        """
        BuiltinFunction.__init__(self, "dirac_delta", latex_name=r"\delta",
                                 conversions=dict(maxima='delta',
                                                  mathematica='DiracDelta',
                                                  sympy='DiracDelta',
                                                  giac='Dirac'))

    def _eval_(self, x):
        """
        INPUT:

        - ``x`` -- a real number or a symbolic expression

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: dirac_delta(1)
            0
            sage: dirac_delta(0)
            dirac_delta(0)
            sage: dirac_delta(x)
            dirac_delta(x)
            sage: dirac_delta(exp(-10000000000000000000))
            0

        Evaluation test::

            sage: dirac_delta(x).subs(x=1)                                              # needs sage.symbolic
            0
        """
        try:
            return self._evalf_(x)
        except (TypeError, ValueError):      # x is symbolic
            pass
        return None

    def _evalf_(self, x, **kwds):
        """
        TESTS::

            sage: h(x) = dirac_delta(x)                                                 # needs sage.symbolic
            sage: h(pi).numerical_approx()                                              # needs sage.symbolic
            0.000000000000000
        """
        approx_x = ComplexIntervalField()(x)
        if bool(approx_x.imag() == 0):      # x is real
            if bool(approx_x.real() == 0):  # x is zero
                return None
            else:
                return 0
        raise ValueError("Numeric evaluation of symbolic expression")


dirac_delta = FunctionDiracDelta()


class FunctionHeaviside(GinacFunction):
    r"""
    The Heaviside step function, `H(x)` (``heaviside(x)``).

    INPUT:

    - ``x`` -- a real number or a symbolic expression

    DEFINITION:

    The Heaviside step function, `H(x)` is defined in Sage as:

        `H(x) = 0` for `x < 0` and `H(x) = 1` for `x > 0`

    .. SEEALSO:: :func:`unit_step()<sage.functions.generalized.FunctionUnitStep>`

    EXAMPLES::

        sage: # needs sage.symbolic
        sage: heaviside(-1)
        0
        sage: heaviside(1)
        1
        sage: heaviside(0)
        heaviside(0)
        sage: heaviside(x)
        heaviside(x)

        sage: heaviside(-1/2)                                                           # needs sage.symbolic
        0
        sage: heaviside(exp(-1000000000000000000000))                                   # needs sage.symbolic
        1

    TESTS::

        sage: heaviside(x)._sympy_()                                                    # needs sympy sage.symbolic
        Heaviside(x)
        sage: heaviside(x).subs(x=1)                                                    # needs sage.symbolic
        1
        sage: heaviside(x).subs(x=-1)                                                   # needs sage.symbolic
        0

    ::

        sage: # needs sage.symbolic
        sage: ex = heaviside(x) + 1
        sage: t = loads(dumps(ex)); t
        heaviside(x) + 1
        sage: bool(t == ex)
        True
        sage: t.subs(x=1)
        2

    REFERENCES:

    -  :wikipedia:`Heaviside_function`
    """
    def __init__(self):
        r"""
        The Heaviside step function, ``heaviside(x)``.

        INPUT:

        - ``x`` -- a real number or a symbolic expression

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: heaviside(-1)
            0
            sage: heaviside(1)
            1
            sage: heaviside(0)
            heaviside(0)
            sage: heaviside(x)
            heaviside(x)
            sage: latex(heaviside(x))
            H\left(x\right)
            sage: heaviside(x)._sympy_()                                                # needs sympy
            Heaviside(x)
            sage: heaviside(x)._giac_()                                                 # needs giac
            Heaviside(sageVARx)
            sage: h(x) = heaviside(x)
            sage: h(pi).numerical_approx()
            1.00000000000000
        """
        GinacFunction.__init__(self, "heaviside", latex_name='H',
                               conversions=dict(maxima='hstep',
                                                mathematica='HeavisideTheta',
                                                sympy='Heaviside',
                                                giac='Heaviside'))

    def _derivative_(self, x, diff_param=None):
        """
        Derivative of Heaviside step function.

        EXAMPLES::

            sage: heaviside(x).diff(x)                                                  # needs sage.symbolic
            dirac_delta(x)
        """
        return dirac_delta(x)


heaviside = FunctionHeaviside()


class FunctionUnitStep(GinacFunction):
    r"""
    The unit step function, `\mathrm{u}(x)` (``unit_step(x)``).

    INPUT:

    - ``x`` -- a real number or a symbolic expression

    DEFINITION:

    The unit step function, `\mathrm{u}(x)` is defined in Sage as:

        `\mathrm{u}(x) = 0` for `x < 0` and `\mathrm{u}(x) = 1` for `x \geq 0`

    .. SEEALSO:: :func:`heaviside()<sage.functions.generalized.FunctionHeaviside>`

    EXAMPLES::

        sage: # needs sage.symbolic
        sage: unit_step(-1)
        0
        sage: unit_step(1)
        1
        sage: unit_step(0)
        1
        sage: unit_step(x)
        unit_step(x)
        sage: unit_step(-exp(-10000000000000000000))
        0

    TESTS::

        sage: # needs sage.symbolic
        sage: unit_step(x).subs(x=1)
        1
        sage: unit_step(x).subs(x=0)
        1
        sage: h(x) = unit_step(x)
        sage: h(pi).numerical_approx()
        1.00000000000000
    """
    def __init__(self):
        r"""
        The unit step function, ``unit_step(x)``.

        INPUT:

        - ``x`` -- a real number or a symbolic expression

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: unit_step(-1)
            0
            sage: unit_step(1)
            1
            sage: unit_step(0)
            1
            sage: unit_step(x)
            unit_step(x)
            sage: latex(unit_step(x))
            \mathrm{u}\left(x\right)

        TESTS::

            sage: t = loads(dumps(unit_step(x) + 1)); t                                 # needs sage.symbolic
            unit_step(x) + 1
            sage: t.subs(x=0)                                                           # needs sage.symbolic
            2
        """
        GinacFunction.__init__(self, "unit_step", latex_name=r"\mathrm{u}",
                               conversions=dict(mathematica='UnitStep'))

    def _derivative_(self, x, diff_param=None):
        """
        Derivative of unit step function.

        EXAMPLES::

            sage: unit_step(x).diff(x)                                                  # needs sage.symbolic
            dirac_delta(x)
        """
        return dirac_delta(x)


unit_step = FunctionUnitStep()


class FunctionSignum(BuiltinFunction):
    r"""
    The signum or sgn function `\mathrm{sgn}(x)` (``sgn(x)``).

    INPUT:

    - ``x`` -- a real number or a symbolic expression

    DEFINITION:

    The sgn function, `\mathrm{sgn}(x)` is defined as:

        `\mathrm{sgn}(x) =  1` for `x > 0`,
        `\mathrm{sgn}(x) =  0` for `x = 0` and
        `\mathrm{sgn}(x) = -1` for `x < 0`

    EXAMPLES::

        sage: sgn(-1)
        -1
        sage: sgn(1)
        1
        sage: sgn(0)
        0
        sage: sgn(x)                                                                    # needs sage.symbolic
        sgn(x)

    We can also use ``sign``::

        sage: sign(1)
        1
        sage: sign(0)
        0
        sage: a = AA(-5).nth_root(7)                                                    # needs sage.rings.number_field
        sage: sign(a)                                                                   # needs sage.rings.number_field
        -1

    TESTS:

    Check if conversions to sympy and others work (:issue:`11921`)::

        sage: sgn(x)._sympy_()                                                          # needs sympy sage.symbolic
        sign(x)
        sage: sgn(x)._fricas_init_()                                                    # needs sage.symbolic
        '(x+->abs(x)/x)(x)'
        sage: sgn(x)._giac_()                                                           # needs giac sage.symbolic
        sign(sageVARx)

    Test for :issue:`31085`::

        sage: fricas(sign(x)).eval(x=-3)        # optional - fricas                     # needs sage.symbolic
        - 1

    REFERENCES:

    - :wikipedia:`Sign_function`
    """
    def __init__(self):
        r"""
        The sgn function, ``sgn(x)``.

        EXAMPLES::

            sage: sgn(-1)
            -1
            sage: sgn(1)
            1
            sage: sgn(0)
            0
            sage: sgn(x)                                                                # needs sage.symbolic
            sgn(x)
            sage: sgn(x)._sympy_()                                                      # needs sympy sage.symbolic
            sign(x)
        """
        BuiltinFunction.__init__(self, "sgn", latex_name=r"\mathrm{sgn}",
                conversions=dict(maxima='signum', mathematica='Sign',
                                 sympy='sign', giac='sign',
                                 fricas='(x+->abs(x)/x)'),
                alt_name='sign')

    def _eval_(self, x):
        """

        EXAMPLES::

            sage: sgn(-1)
            -1
            sage: sgn(1)
            1
            sage: sgn(0)
            0
            sage: sgn(x)                                                                # needs sage.symbolic
            sgn(x)
            sage: sgn(-exp(-10000000000000000000))                                      # needs sage.symbolic
            -1

        Evaluation test::

            sage: sgn(x).subs(x=1)                                                      # needs sage.symbolic
            1
            sage: sgn(x).subs(x=0)                                                      # needs sage.symbolic
            0
            sage: sgn(x).subs(x=-1)                                                     # needs sage.symbolic
            -1

        More tests::

            sage: sign(RR(2))
            1
            sage: sign(RDF(2))
            1
            sage: sign(AA(-2))                                                          # needs sage.rings.number_field
            -1
            sage: sign(AA(0))                                                           # needs sage.rings.number_field
            0
        """
        try:
            return self._evalf_(x)
        except (TypeError, ValueError):      # x is symbolic
            pass
        return None

    def _evalf_(self, x, **kwds):
        """
        TESTS:

        Check that :issue:`16587` is fixed::

            sage: # needs sage.symbolic
            sage: M = sgn(3/2, hold=True); M
            sgn(3/2)
            sage: M.n()
            1
            sage: h(x) = sgn(x)
            sage: h(pi).numerical_approx()
            1.00000000000000
        """
        if hasattr(x, 'sign'):  # First check if x has a sign method
            return x.sign()
        if hasattr(x, 'sgn'):  # or a sgn method
            return x.sgn()
        approx_x = ComplexIntervalField()(x)
        if bool(approx_x.imag() == 0):      # x is real
            if bool(approx_x.real() == 0):  # x is zero
                return ZZ(0)
            # Now we have a nonzero real
            if bool((approx_x**(0.5)).imag() == 0):  # Check: x > 0
                return ZZ(1)
            else:
                return ZZ(-1)
        raise ValueError("Numeric evaluation of symbolic expression")

    def _derivative_(self, x, diff_param=None):
        """
        Derivative of sgn function.

        EXAMPLES::

            sage: sgn(x).diff(x)                                                        # needs sage.symbolic
            2*dirac_delta(x)
        """
        assert diff_param == 0
        return 2 * dirac_delta(x)


sgn = FunctionSignum()
sign = sgn


class FunctionKroneckerDelta(BuiltinFunction):
    r"""
    The Kronecker delta function `\delta_{m,n}` (``kronecker_delta(m, n)``).

    INPUT:

    - ``m`` -- a number or a symbolic expression
    - ``n`` -- a number or a symbolic expression

    DEFINITION:

    Kronecker delta function `\delta_{m,n}` is defined as:

        `\delta_{m,n} = 0` for `m \ne n` and
        `\delta_{m,n} = 1` for `m = n`

    EXAMPLES::

        sage: kronecker_delta(1, 2)                                                     # needs sage.rings.complex_interval_field
        0
        sage: kronecker_delta(1, 1)                                                     # needs sage.rings.complex_interval_field
        1
        sage: m, n = var('m,n')                                                         # needs sage.symbolic
        sage: kronecker_delta(m, n)                                                     # needs sage.symbolic
        kronecker_delta(m, n)

    REFERENCES:

    - :wikipedia:`Kronecker_delta`
    """
    def __init__(self):
        r"""
        The Kronecker delta function.

        EXAMPLES::

            sage: kronecker_delta(1, 2)                                                 # needs sage.rings.complex_interval_field
            0
            sage: kronecker_delta(1, 1)                                                 # needs sage.rings.complex_interval_field
            1
            sage: y = var('y')                                                          # needs sage.symbolic
            sage: kronecker_delta(x, y)._sympy_()                                       # needs sympy sage.symbolic
            KroneckerDelta(x, y)
        """
        BuiltinFunction.__init__(self, "kronecker_delta", nargs=2,
                                 conversions=dict(maxima='kron_delta',
                                                  mathematica='KroneckerDelta',
                                                  sympy='KroneckerDelta'))

    def _eval_(self, m, n):
        """
        The Kronecker delta function.

        EXAMPLES::

            sage: kronecker_delta(1, 2)                                                 # needs sage.rings.complex_interval_field
            0
            sage: kronecker_delta(1, 1)                                                 # needs sage.rings.complex_interval_field
            1

        Kronecker delta is a symmetric function. We keep arguments sorted to
        ensure that k_d(m, n) - k_d(n, m) cancels automatically::

            sage: # needs sage.symbolic
            sage: x, y = var('x,y')
            sage: kronecker_delta(x, y)
            kronecker_delta(x, y)
            sage: kronecker_delta(y, x)
            kronecker_delta(x, y)
            sage: kronecker_delta(x, 2*x)
            kronecker_delta(2*x, x)

        Evaluation test::

            sage: kronecker_delta(1, x).subs(x=1)                                       # needs sage.symbolic
            1
        """
        try:
            return self._evalf_(m, n)
        except (TypeError, ValueError):      # x is symbolic
            pass
        return None

    def _evalf_(self, m, n, **kwds):
        """
        TESTS::

            sage: h(x) = kronecker_delta(3, x)                                          # needs sage.symbolic
            sage: h(pi).numerical_approx()                                              # needs sage.symbolic
            0.000000000000000
        """
        if bool(repr(m) > repr(n)):
            return kronecker_delta(n, m)
        x = m - n
        approx_x = ComplexIntervalField()(x)
        if approx_x.imag() == 0:      # x is real
            if approx_x.real() == 0:  # x is zero
                return 1
            else:
                return 0
        return 0            # x is complex

    def _derivative_(self, *args, **kwds):
        """
        Derivative of Kronecker delta.

        EXAMPLES::

            sage: kronecker_delta(x, 1).diff(x)                                         # needs sage.symbolic
            0
        """
        # Kronecker delta is nonzero (but finite) only in the set of
        # zero-measure unlike Dirac delta. Consequently, it is null
        # for the purpose of integration/differentiation. For *discrete sum*
        # Kronecker delta is however non-trivial.
        return 0

    def _print_latex_(self, m, n, **kwds):
        r"""
        Return latex expression.

        EXAMPLES::

            sage: from sage.misc.latex import latex
            sage: m, n = var('m,n')                                                     # needs sage.symbolic
            sage: latex(kronecker_delta(m, n))                                          # needs sage.symbolic
            \delta_{m,n}
        """
        return r"\delta_{%s,%s}" % (latex(m), latex(n))


kronecker_delta = FunctionKroneckerDelta()
