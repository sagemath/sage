"""
Gamma and related functions
"""
from sage.misc.lazy_import import lazy_import
from sage.rings.infinity import Infinity
from sage.rings.rational import Rational
from sage.structure.element import parent as s_parent
from sage.symbolic.function import GinacFunction, BuiltinFunction
from sage.symbolic.symbols import register_symbol, symbol_table

lazy_import('sage.rings.complex_mpfr', 'ComplexField')

lazy_import('sage.functions.error', 'erf')
lazy_import('sage.functions.exp_integral', 'Ei')
lazy_import('sage.functions.log', 'exp')
lazy_import('sage.functions.other', 'sqrt')

lazy_import('sage.symbolic.constants', 'pi')

lazy_import('sage.libs.mpmath.utils', 'call', as_='_mpmath_utils_call')
lazy_import('mpmath', 'gammainc', as_='_mpmath_gammainc')


class Function_gamma(GinacFunction):
    def __init__(self):
        r"""
        The Gamma function.  This is defined by

        .. MATH::

            \Gamma(z) = \int_0^\infty t^{z-1}e^{-t} dt

        for complex input `z` with real part greater than zero, and by
        analytic continuation on the rest of the complex plane (except
        for negative integers, which are poles).

        It is computed by various libraries within Sage, depending on
        the input type.

        EXAMPLES::

            sage: from sage.functions.gamma import gamma1
            sage: gamma1(CDF(0.5, 14))                                                  # needs sage.libs.pari sage.rings.complex_double
            -4.0537030780372815e-10 - 5.773299834553605e-10*I
            sage: gamma1(CDF(I))                                                        # needs sage.libs.pari sage.rings.complex_double sage.symbolic
            -0.15494982830181067 - 0.49801566811835607*I

        Recall that `\Gamma(n)` is `n-1` factorial::

            sage: gamma1(11) == factorial(10)
            True
            sage: gamma1(6)
            120
            sage: gamma1(1/2)                                                           # needs sage.symbolic
            sqrt(pi)
            sage: gamma1(-1)
            Infinity
            sage: gamma1(I)                                                             # needs sage.symbolic
            gamma(I)
            sage: gamma1(x/2)(x=5)                                                      # needs sage.symbolic
            3/4*sqrt(pi)

            sage: gamma1(float(6))  # For ARM: rel tol 3e-16
            120.0
            sage: gamma(6.)                                                             # needs sage.symbolic
            120.000000000000
            sage: gamma1(x)                                                             # needs sage.symbolic
            gamma(x)

        ::

            sage: gamma1(pi)                                                            # needs sage.symbolic
            gamma(pi)
            sage: gamma1(i)                                                             # needs sage.symbolic
            gamma(I)
            sage: gamma1(i).n()                                                         # needs sage.symbolic
            -0.154949828301811 - 0.498015668118356*I
            sage: gamma1(int(5))
            24

        ::

            sage: conjugate(gamma(x))                                                   # needs sage.symbolic
            gamma(conjugate(x))

        ::

            sage: plot(gamma1(x), (x,1,5))                                              # needs sage.plot sage.symbolic
            Graphics object consisting of 1 graphics primitive

        We are also able to compute the Laurent expansion of the
        Gamma function (as well as of functions containing
        the Gamma function)::

            sage: gamma(x).series(x==0, 2)                                              # needs sage.symbolic
            1*x^(-1) + (-euler_gamma)
            + (1/2*euler_gamma^2 + 1/12*pi^2)*x + Order(x^2)
            sage: (gamma(x)^2).series(x==0, 1)                                          # needs sage.symbolic
            1*x^(-2) + (-2*euler_gamma)*x^(-1)
            + (2*euler_gamma^2 + 1/6*pi^2) + Order(x)

        To prevent automatic evaluation, use the ``hold`` argument::

            sage: gamma1(1/2, hold=True)                                                # needs sage.symbolic
            gamma(1/2)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: gamma1(1/2, hold=True).simplify()                                     # needs sage.symbolic
            sqrt(pi)

        TESTS:

            sage: gamma(x)._sympy_()                                                    # needs sympy sage.symbolic
            gamma(x)

        We verify that we can convert this function to Maxima and
        convert back to Sage::

            sage: z = var('z')                                                          # needs sage.symbolic
            sage: maxima(gamma1(z)).sage()                                              # needs sage.symbolic
            gamma(z)
            sage: latex(gamma1(z))                                                      # needs sage.symbolic
            \Gamma\left(z\right)

        Test that :issue:`5556` is fixed::

            sage: gamma1(3/4)                                                           # needs sage.symbolic
            gamma(3/4)

            sage: gamma1(3/4).n(100)                                                    # needs sage.symbolic
            1.2254167024651776451290983034

        Check that negative integer input works::

            sage: (-1).gamma()
            Infinity
            sage: (-1.).gamma()                                                         # needs sage.rings.real_mpfr
            NaN
            sage: CC(-1).gamma()                                                        # needs sage.libs.pari sage.rings.real_mpfr
            Infinity
            sage: RDF(-1).gamma()                                                       # needs sage.rings.real_mpfr
            NaN
            sage: CDF(-1).gamma()                                                       # needs sage.libs.pari sage.rings.complex_double
            Infinity

        Check if :issue:`8297` is fixed::

            sage: latex(gamma(1/4))                                                     # needs sage.symbolic
            \Gamma\left(\frac{1}{4}\right)

        Test pickling::

            sage: loads(dumps(gamma(x)))                                                # needs sage.symbolic
            gamma(x)

        Check that the implementations roughly agrees (note there might be
        difference of several ulp on more complicated entries)::

            sage: import mpmath                                                         # needs mpmath
            sage: float(gamma(10.)) == gamma(10.r) == float(gamma(mpmath.mpf(10)))      # needs mpmath
            True
            sage: float(gamma(8.5)) == gamma(8.5r) == float(gamma(mpmath.mpf(8.5)))     # needs mpmath
            True

        Check that ``QQbar`` half integers work with the ``pi`` formula::

            sage: gamma(QQbar(1/2))                                                     # needs sage.rings.number_field sage.symbolic
            sqrt(pi)
            sage: gamma(QQbar(-9/2))                                                    # needs sage.rings.number_field sage.symbolic
            -32/945*sqrt(pi)

        .. SEEALSO::

            :meth:`gamma`
        """
        GinacFunction.__init__(self, 'gamma', latex_name=r"\Gamma",
                               ginac_name='gamma',
                               conversions={'mathematica':'Gamma',
                                            'maple':'GAMMA',
                                            'sympy':'gamma',
                                            'fricas':'Gamma',
                                            'giac':'Gamma'})


gamma1 = Function_gamma()


class Function_log_gamma(GinacFunction):
    def __init__(self):
        r"""
        The principal branch of the log gamma function. Note that for
        `x < 0`, ``log(gamma(x))`` is not, in general, equal to
        ``log_gamma(x)``.

        It is computed by the ``log_gamma`` function for the number type,
        or by ``lgamma`` in Ginac, failing that.

        Gamma is defined for complex input `z` with real part greater
        than zero, and by analytic continuation on the rest of the
        complex plane (except for negative integers, which are poles).

        EXAMPLES:

        Numerical evaluation happens when appropriate, to the
        appropriate accuracy (see :issue:`10072`)::

            sage: # needs sage.symbolic
            sage: log_gamma(6)
            log(120)
            sage: log_gamma(6.)
            4.78749174278205
            sage: log_gamma(6).n()
            4.78749174278205
            sage: log_gamma(2.4 + I)
            -0.0308566579348816 + 0.693427705955790*I
            sage: log_gamma(-3.1)
            0.400311696703985 - 12.5663706143592*I
            sage: log_gamma(-1.1) == log(gamma(-1.1))
            False

            sage: log_gamma(RealField(100)(6))                                          # needs sage.rings.real_mpfr
            4.7874917427820459942477009345

        Symbolic input works (see :issue:`10075`)::

            sage: log_gamma(3*x)                                                        # needs sage.symbolic
            log_gamma(3*x)
            sage: log_gamma(3 + I)                                                      # needs sage.symbolic
            log_gamma(I + 3)
            sage: log_gamma(3 + I + x)                                                  # needs sage.symbolic
            log_gamma(x + I + 3)

        Check that :issue:`12521` is fixed::

            sage: # needs sage.symbolic
            sage: log_gamma(-2.1)
            1.53171380819509 - 9.42477796076938*I
            sage: log_gamma(CC(-2.1))
            1.53171380819509 - 9.42477796076938*I
            sage: log_gamma(-21/10).n()
            1.53171380819509 - 9.42477796076938*I
            sage: exp(log_gamma(-1.3) + log_gamma(-0.4) -
            ....:     log_gamma(-1.3 - 0.4)).real_part()  # beta(-1.3, -0.4)
            -4.92909641669610

        In order to prevent evaluation, use the ``hold`` argument;
        to evaluate a held expression, use the ``n()`` numerical
        evaluation method::

            sage: log_gamma(SR(5), hold=True)                                           # needs sage.symbolic
            log_gamma(5)
            sage: log_gamma(SR(5), hold=True).n()                                       # needs sage.symbolic
            3.17805383034795

        TESTS::

            sage: log_gamma(pari(6))                                                    # needs sage.libs.pari
            4.78749174278205

            sage: # needs sage.symbolic
            sage: log_gamma(-2.1 + I)
            -1.90373724496982 - 7.18482377077183*I
            sage: log_gamma(x)._sympy_()                                                # needs sympy
            loggamma(x)
            sage: log_gamma(CC(6))                                                      # needs sage.rings.real_mpfr
            4.78749174278205
            sage: log_gamma(CC(-2.5))                                                   # needs sage.rings.real_mpfr
            -0.0562437164976741 - 9.42477796076938*I
            sage: log_gamma(RDF(-2.5))
            -0.056243716497674054 - 9.42477796076938*I
            sage: log_gamma(CDF(-2.5))                                                  # needs sage.rings.complex_double
            -0.056243716497674054 - 9.42477796076938*I
            sage: log_gamma(float(-2.5))
            (-0.056243716497674054-9.42477796076938j)
            sage: log_gamma(complex(-2.5))                                              # needs sage.rings.complex_double
            (-0.056243716497674054-9.42477796076938j)

        ``conjugate(log_gamma(x)) == log_gamma(conjugate(x))`` unless on the
        branch cut, which runs along the negative real axis.::

            sage: # needs sage.symbolic
            sage: conjugate(log_gamma(x))
            conjugate(log_gamma(x))
            sage: var('y', domain='positive')
            y
            sage: conjugate(log_gamma(y))
            log_gamma(y)
            sage: conjugate(log_gamma(y + I))
            conjugate(log_gamma(y + I))
            sage: log_gamma(-2)
            +Infinity
            sage: conjugate(log_gamma(-2))
            +Infinity
        """
        GinacFunction.__init__(self, "log_gamma", latex_name=r'\log\Gamma',
                               conversions=dict(mathematica='LogGamma',
                                                maxima='log_gamma',
                                                sympy='loggamma',
                                                fricas='logGamma'))


log_gamma = Function_log_gamma()


class Function_gamma_inc(BuiltinFunction):
    def __init__(self):
        r"""
        The upper incomplete gamma function.

        It is defined by the integral

        .. MATH::

            \Gamma(a,z)=\int_z^\infty t^{a-1}e^{-t}\,\mathrm{d}t

        EXAMPLES::

            sage: gamma_inc(CDF(0,1), 3)                                                # needs sage.libs.pari sage.rings.complex_double
            0.0032085749933691158 + 0.012406185811871568*I
            sage: gamma_inc(RDF(1), 3)                                                  # needs sage.rings.complex_double
            0.049787068367863944

            sage: # needs sage.symbolic
            sage: gamma_inc(3, 2)
            gamma(3, 2)
            sage: gamma_inc(x, 0)
            gamma(x)
            sage: latex(gamma_inc(3, 2))
            \Gamma\left(3, 2\right)
            sage: loads(dumps((gamma_inc(3, 2))))
            gamma(3, 2)

            sage: i = ComplexField(30).0; gamma_inc(2, 1 + i)                           # needs sage.rings.real_mpfr
            0.70709210 - 0.42035364*I
            sage: gamma_inc(2., 5)                                                      # needs sage.rings.complex_double
            0.0404276819945128

            sage: x, y = var('x,y')                                                     # needs sage.symbolic
            sage: gamma_inc(x,y).diff(x)                                                # needs sage.symbolic
            diff(gamma(x, y), x)
            sage: (gamma_inc(x, x + 1).diff(x)).simplify()                              # needs sage.symbolic
            -(x + 1)^(x - 1)*e^(-x - 1) + D[0](gamma)(x, x + 1)

        TESTS:

        Check that :issue:`21407` is fixed::

            sage: gamma(-1, 5)._sympy_()                                                # needs sympy sage.symbolic
            expint(2, 5)/5
            sage: gamma(-3/2, 5)._sympy_()                                              # needs sympy sage.symbolic
            -6*sqrt(5)*exp(-5)/25 + 4*sqrt(pi)*erfc(sqrt(5))/3

        Check that :issue:`25597` is fixed::

            sage: gamma(-1, 5)._fricas_()                                       # optional - fricas, needs sage.symbolic
            Gamma(- 1,5)

            sage: var('t')                                                              # needs sage.symbolic
            t
            sage: integrate(-exp(-x)*x^(t-1), x, algorithm='fricas')            # optional - fricas, needs sage.symbolic
            gamma(t, x)

        .. SEEALSO::

            :meth:`gamma`
        """
        BuiltinFunction.__init__(self, "gamma", nargs=2, latex_name=r"\Gamma",
                conversions={'maxima':'gamma_incomplete', 'mathematica':'Gamma',
                             'maple':'GAMMA', 'sympy':'uppergamma', 'fricas':'Gamma',
                             'giac':'ugamma'})

    def _method_arguments(self, x, y):
        r"""
        TESTS::

            sage: # needs sage.libs.flint
            sage: b = RBF(1, 1e-10)
            sage: gamma(b) # abs tol 1e-9
            [1.0000000000 +/- 5.78e-11]
            sage: gamma(CBF(b)) # abs tol 1e-9
            [1.0000000000 +/- 5.78e-11]
            sage: gamma(CBF(b), 4) # abs tol 2e-9
            [0.018315639 +/- 9.00e-10]
            sage: gamma(CBF(1), b)
            [0.3678794412 +/- 6.54e-11]
        """
        return [x, y]

    def _eval_(self, x, y):
        """
        EXAMPLES::

            sage: gamma_inc(2., 0)                                                      # needs sage.rings.complex_double
            1.00000000000000
            sage: gamma_inc(2, 0)                                                       # needs sage.rings.real_mpfr
            1

            sage: # needs sage.symbolic
            sage: gamma_inc(1/2, 2)
            -sqrt(pi)*(erf(sqrt(2)) - 1)
            sage: gamma_inc(1/2, 1)
            -sqrt(pi)*(erf(1) - 1)
            sage: gamma_inc(1/2, 0)
            sqrt(pi)
            sage: gamma_inc(x, 0)
            gamma(x)
            sage: gamma_inc(1, 2)
            e^(-2)
            sage: gamma_inc(0, 2)
            -Ei(-2)
        """
        if y == 0:
            return gamma(x)
        if x == 1:
            return exp(-y)
        if x == 0:
            return -Ei(-y)
        if x == Rational((1, 2)):  # only for x>0
            return sqrt(pi) * (1 - erf(sqrt(y)))
        return None

    def _evalf_(self, x, y, parent=None, algorithm='pari'):
        """
        EXAMPLES::

            sage: gamma_inc(0, 2)                                                       # needs sage.symbolic
            -Ei(-2)
            sage: gamma_inc(0, 2.)                                                      # needs sage.rings.real_mpfr
            0.0489005107080611
            sage: gamma_inc(0, 2).n(algorithm='pari')                                   # needs sage.libs.pari sage.symbolic
            0.0489005107080611
            sage: gamma_inc(0, 2).n(200)                                                # needs sage.symbolic
            0.048900510708061119567239835228...
            sage: gamma_inc(3, 2).n()                                                   # needs sage.symbolic
            1.35335283236613

        TESTS:

        Check that :issue:`7099` is fixed::

            sage: R = RealField(1024)                                                   # needs sage.rings.real_mpfr
            sage: gamma(R(9), R(10^-3))  # rel tol 1e-308                               # needs sage.rings.real_mpfr
            40319.99999999999999999999999999988898884344822911869926361916294165058203634104838326009191542490601781777105678829520585311300510347676330951251563007679436243294653538925717144381702105700908686088851362675381239820118402497959018315224423868693918493033078310647199219674433536605771315869983788442389633
            sage: numerical_approx(gamma(9, 10^(-3)) - gamma(9), digits=40)  # abs tol 1e-36        # needs sage.symbolic
            -1.110111598370794007949063502542063148294e-28

        Check that :issue:`17328` is fixed::

            sage: gamma_inc(float(-1), float(-1))                                       # needs sage.rings.real_mpfr
            (-0.8231640121031085+3.141592653589793j)
            sage: gamma_inc(RR(-1), RR(-1))                                             # needs sage.rings.complex_double
            -0.823164012103109 + 3.14159265358979*I
            sage: gamma_inc(-1, float(-log(3))) - gamma_inc(-1, float(-log(2)))  # abs tol 1e-15    # needs sage.symbolic
            (1.2730972164471142+0j)

        Check that :issue:`17130` is fixed::

            sage: r = gamma_inc(float(0), float(1)); r                                  # needs sage.rings.real_mpfr
            0.21938393439552029
            sage: type(r)                                                               # needs sage.rings.real_mpfr
            <... 'float'>
        """
        R = parent or s_parent(x)
        # C is the complex version of R
        # prec is the precision of R
        if R is float:
            prec = 53
            C = complex
        else:
            try:
                prec = R.precision()
            except AttributeError:
                prec = 53
            try:
                C = R.complex_field()
            except AttributeError:
                C = R

        if algorithm == 'pari':
            v = ComplexField(prec)(x).gamma_inc(y)
        else:
            v = ComplexField(prec)(_mpmath_utils_call(_mpmath_gammainc, x, y, parent=R))
        if v.is_real():
            return R(v)
        else:
            return C(v)


# shorter alias
gamma_inc = Function_gamma_inc()


class Function_gamma_inc_lower(BuiltinFunction):
    def __init__(self):
        r"""
        The lower incomplete gamma function.

        It is defined by the integral

        .. MATH::

            \Gamma(a,z)=\int_0^z t^{a-1}e^{-t}\,\mathrm{d}t

        EXAMPLES::

            sage: gamma_inc_lower(CDF(0,1), 3)                                          # needs sage.rings.complex_double
            -0.1581584032951798 - 0.5104218539302277*I
            sage: gamma_inc_lower(RDF(1), 3)                                            # needs sage.rings.complex_double
            0.950212931632136

            sage: # needs sage.symbolic
            sage: gamma_inc_lower(3, 2, hold=True)
            gamma_inc_lower(3, 2)
            sage: gamma_inc_lower(3, 2)
            -10*e^(-2) + 2
            sage: gamma_inc_lower(x, 0)
            0
            sage: latex(gamma_inc_lower(x, x))
            \gamma\left(x, x\right)
            sage: loads(dumps((gamma_inc_lower(x, x))))
            gamma_inc_lower(x, x)

            sage: i = ComplexField(30).0; gamma_inc_lower(2, 1 + i)                     # needs sage.rings.real_mpfr
            0.29290790 + 0.42035364*I
            sage: gamma_inc_lower(2., 5)                                                # needs sage.rings.complex_double
            0.959572318005487

        Interfaces to other software::

            sage: gamma_inc_lower(x, x)._sympy_()                                       # needs sympy sage.symbolic
            lowergamma(x, x)
            sage: maxima(gamma_inc_lower(x, x))                                         # needs sage.symbolic
            gamma_incomplete_lower(_SAGE_VAR_x,_SAGE_VAR_x)

        .. SEEALSO::

        :class:`Function_gamma_inc`
        """
        BuiltinFunction.__init__(self, "gamma_inc_lower", nargs=2, latex_name=r"\gamma",
                conversions={'maxima':'gamma_incomplete_lower',
                    'maple':'GAMMA', 'sympy':'lowergamma', 'giac':'igamma'})

    def _eval_(self, x, y):
        """
        EXAMPLES::

            sage: gamma_inc_lower(2., 0)                                                # needs sage.rings.complex_double
            0.000000000000000

            sage: # needs sage.symbolic
            sage: gamma_inc_lower(2, 0)
            0
            sage: gamma_inc_lower(1/2, 2)
            sqrt(pi)*erf(sqrt(2))
            sage: gamma_inc_lower(1/2, 1)
            sqrt(pi)*erf(1)
            sage: gamma_inc_lower(1/2, 0)
            0
            sage: gamma_inc_lower(x, 0)
            0
            sage: gamma_inc_lower(1, 2)
            -e^(-2) + 1
            sage: gamma_inc_lower(0, 2)
            +Infinity
            sage: gamma_inc_lower(2, 377/79)
            -456/79*e^(-377/79) + 1
            sage: gamma_inc_lower(3, x)
            -(x^2 + 2*x + 2)*e^(-x) + 2
            sage: gamma_inc_lower(9/2, 37/7)
            -1/38416*sqrt(pi)*(1672946*sqrt(259)*e^(-37/7)/sqrt(pi) - 252105*erf(1/7*sqrt(259)))
        """
        if y == 0:
            return 0
        if x == 0:
            return Infinity
        elif x == 1:
            return 1 - exp(-y)
        elif (2 * x).is_integer():
            return self(x, y, hold=True)._sympy_()
        else:
            return None

    def _evalf_(self, x, y, parent=None, algorithm='mpmath'):
        """
        EXAMPLES::

            sage: gamma_inc_lower(3, 2.)                                                # needs sage.rings.real_mpfr
            0.646647167633873
            sage: gamma_inc_lower(3, 2).n(200)                                          # needs sage.symbolic
            0.646647167633873081060005050275155...
            sage: gamma_inc_lower(0, 2.)                                                # needs sage.rings.real_mpfr
            +infinity
        """
        R = parent or s_parent(x)
        # C is the complex version of R
        # prec is the precision of R
        if R is float:
            prec = 53
            C = complex
        else:
            try:
                prec = R.precision()
            except AttributeError:
                prec = 53
            try:
                C = R.complex_field()
            except AttributeError:
                C = R
        if algorithm == 'pari':
            Cx = ComplexField(prec)(x)
            v = Cx.gamma() - Cx.gamma_inc(y)
        else:
            v = ComplexField(prec)(_mpmath_utils_call(_mpmath_gammainc, x, 0, y, parent=R))
        return R(v) if v.is_real() else C(v)

    def _derivative_(self, x, y, diff_param=None):
        """
        EXAMPLES::

            sage: x, y = var('x,y')                                                     # needs sage.symbolic
            sage: gamma_inc_lower(x, y).diff(y)                                         # needs sage.symbolic
            y^(x - 1)*e^(-y)
            sage: gamma_inc_lower(x, y).diff(x)                                         # needs sage.symbolic
            Traceback (most recent call last):
            ...
            NotImplementedError: cannot differentiate gamma_inc_lower in the first parameter
        """
        if diff_param == 0:
            raise NotImplementedError("cannot differentiate gamma_inc_lower in the"
                                      " first parameter")
        else:
            return exp(-y) * y**(x - 1)

    def _mathematica_init_evaled_(self, *args):
        r"""
        EXAMPLES::

            sage: gamma_inc_lower(4/3, 1)._mathematica_()  # indirect doctest   # optional - mathematica, needs sage.symbolic
            Gamma[4/3, 0, 1]
        """
        args_mathematica = []
        for a in args:
            if isinstance(a, str):
                args_mathematica.append(a)
            elif hasattr(a, '_mathematica_init_'):
                args_mathematica.append(a._mathematica_init_())
            else:
                args_mathematica.append(str(a))
        x, z = args_mathematica
        return "Gamma[%s,0,%s]" % (x, z)


# shorter alias
gamma_inc_lower = Function_gamma_inc_lower()


def gamma(a, *args, **kwds):
    r"""
    Gamma and upper incomplete gamma functions in one symbol.

    Recall that `\Gamma(n)` is `n-1` factorial::

        sage: gamma(11) == factorial(10)
        True
        sage: gamma(6)
        120
        sage: gamma(1/2)                                                                # needs sage.symbolic
        sqrt(pi)
        sage: gamma(-4/3)                                                               # needs sage.symbolic
        gamma(-4/3)
        sage: gamma(-1)
        Infinity
        sage: gamma(0)
        Infinity

    ::

        sage: gamma_inc(3, 2)                                                           # needs sage.symbolic
        gamma(3, 2)
        sage: gamma_inc(x,0)                                                            # needs sage.symbolic
        gamma(x)

    ::

        sage: gamma(5, hold=True)                                                       # needs sage.symbolic
        gamma(5)
        sage: gamma(x, 0, hold=True)                                                    # needs sage.symbolic
        gamma(x, 0)

    ::

        sage: gamma(CDF(I))                                                             # needs sage.libs.pari sage.rings.complex_double sage.symbolic
        -0.15494982830181067 - 0.49801566811835607*I
        sage: gamma(CDF(0.5, 14))                                                       # needs sage.libs.pari sage.rings.complex_double
        -4.0537030780372815e-10 - 5.773299834553605e-10*I

    Use ``numerical_approx`` to get higher precision from
    symbolic expressions::

        sage: gamma(pi).n(100)                                                          # needs sage.symbolic
        2.2880377953400324179595889091
        sage: gamma(3/4).n(100)                                                         # needs sage.symbolic
        1.2254167024651776451290983034

    The precision for the result is also deduced from the precision of the
    input. Convert the input to a higher precision explicitly if a result
    with higher precision is desired.::

        sage: t = gamma(RealField(100)(2.5)); t                                         # needs sage.rings.real_mpfr
        1.3293403881791370204736256125
        sage: t.prec()                                                                  # needs sage.rings.real_mpfr
        100

    The gamma function only works with input that can be coerced to the
    Symbolic Ring::

        sage: x = polygen(ZZ, 'x')
        sage: Q.<i> = NumberField(x^2 + 1)                                              # needs sage.rings.number_field
        sage: gamma(i)                                                                  # needs sage.rings.number_field sage.symbolic
        Traceback (most recent call last):
        ...
        TypeError: cannot coerce arguments: no canonical coercion
        from Number Field in i with defining polynomial x^2 + 1 to Symbolic Ring

    .. SEEALSO::

        :class:`Function_gamma`
    """
    if not args:
        return gamma1(a, **kwds)
    if len(args) > 1:
        raise TypeError("Symbolic function gamma takes at most 2 arguments (%s given)" % (len(args) + 1))
    return gamma_inc(a, args[0], **kwds)


# We have to add the wrapper function manually to the symbol_table when we have
# two functions with different number of arguments and the same name
symbol_table['functions']['gamma'] = gamma


def _mathematica_gamma3(*args):
    r"""
    EXAMPLES::

        sage: gamma(4/3)._mathematica_().sage()       # indirect doctest, optional - mathematica
        gamma(4/3)
        sage: gamma(4/3, 1)._mathematica_().sage()    # indirect doctest, optional - mathematica
        gamma(4/3, 1)
        sage: mathematica('Gamma[4/3, 0, 1]').sage()  # indirect doctest, optional - mathematica
        gamma(4/3) - gamma(4/3, 1)
    """
    assert len(args) == 3
    return gamma_inc(args[0], args[1]) - gamma_inc(args[0], args[2])


register_symbol(_mathematica_gamma3, dict(mathematica='Gamma'), 3)


class Function_psi1(GinacFunction):
    def __init__(self):
        r"""
        The digamma function, `\psi(x)`, is the logarithmic derivative of the
        gamma function.

        .. MATH::

            \psi(x) = \frac{d}{dx} \log(\Gamma(x)) = \frac{\Gamma'(x)}{\Gamma(x)}

        EXAMPLES::

            sage: from sage.functions.gamma import psi1
            sage: psi1(x)                                                               # needs sage.symbolic
            psi(x)
            sage: psi1(x).derivative(x)                                                 # needs sage.symbolic
            psi(1, x)

        ::

            sage: psi1(3)                                                               # needs sage.symbolic
            -euler_gamma + 3/2

        ::

            sage: psi(.5)                                                               # needs sage.symbolic
            -1.96351002602142
            sage: psi(RealField(100)(.5))                                               # needs sage.rings.real_mpfr sage.symbolic
            -1.9635100260214234794409763330

        TESTS::

            sage: latex(psi1(x))                                                        # needs sage.symbolic
            \psi\left(x\right)
            sage: loads(dumps(psi1(x) + 1))                                             # needs sage.symbolic
            psi(x) + 1

            sage: # needs sage.symbolic
            sage: t = psi1(x); t
            psi(x)
            sage: t.subs(x=.2)
            -5.28903989659219
            sage: psi(x)._sympy_()                                                      # needs sympy
            polygamma(0, x)
            sage: psi(x)._fricas_()             # optional - fricas
            digamma(x)
        """
        GinacFunction.__init__(self, "psi", nargs=1, latex_name=r'\psi',
                               conversions=dict(mathematica='PolyGamma',
                                                maxima='psi[0]',
                                                maple='Psi',
                                                sympy='digamma',
                                                fricas='digamma'))


class Function_psi2(GinacFunction):
    def __init__(self):
        r"""
        Derivatives of the digamma function `\psi(x)`.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: from sage.functions.gamma import psi2
            sage: psi2(2, x)
            psi(2, x)
            sage: psi2(2, x).derivative(x)
            psi(3, x)
            sage: n = var('n')
            sage: psi2(n, x).derivative(x)
            psi(n + 1, x)

        ::

            sage: psi2(0, x)                                                            # needs sage.symbolic
            psi(x)
            sage: psi2(-1, x)                                                           # needs sage.symbolic
            log(gamma(x))
            sage: psi2(3, 1)                                                            # needs sage.symbolic
            1/15*pi^4

        ::

            sage: psi2(2, .5).n()                                                       # needs sage.symbolic
            -16.8287966442343
            sage: psi2(2, .5).n(100)                                                    # needs sage.symbolic
            -16.828796644234319995596334261

        TESTS::

            sage: psi2(n, x).derivative(n)                                              # needs sage.symbolic
            Traceback (most recent call last):
            ...
            RuntimeError: cannot diff psi(n,x) with respect to n

            sage: # needs sage.symbolic
            sage: latex(psi2(2,x))
            \psi\left(2, x\right)
            sage: loads(dumps(psi2(2,x) + 1))
            psi(2, x) + 1
            sage: psi(2, x)._sympy_()                                                   # needs sympy
            polygamma(2, x)
            sage: psi(2, x)._fricas_()          # optional - fricas
            polygamma(2,x)

        Fixed conversion::

            sage: psi(2, x)._maple_init_()                                              # needs sage.symbolic
            'Psi(2,x)'
        """
        GinacFunction.__init__(self, "psi", nargs=2, latex_name=r'\psi',
                               conversions=dict(mathematica='PolyGamma',
                                                sympy='polygamma',
                                                maple='Psi',
                                                giac='Psi',
                                                fricas='polygamma'))

    def _maxima_init_evaled_(self, *args):
        """
        EXAMPLES:

        These are indirect doctests for this function.::

            sage: from sage.functions.gamma import psi2
            sage: psi2(2, x)._maxima_()                                                 # needs sage.symbolic
            psi[2](_SAGE_VAR_x)
            sage: psi2(4, x)._maxima_()                                                 # needs sage.symbolic
            psi[4](_SAGE_VAR_x)
        """
        args_maxima = []
        for a in args:
            if isinstance(a, str):
                args_maxima.append(a)
            elif hasattr(a, '_maxima_init_'):
                args_maxima.append(a._maxima_init_())
            else:
                args_maxima.append(str(a))
        n, x = args_maxima
        return "psi[%s](%s)" % (n, x)


psi1 = Function_psi1()
psi2 = Function_psi2()


def psi(x, *args, **kwds):
    r"""
    The digamma function, `\psi(x)`, is the logarithmic derivative of the
    gamma function.

    .. MATH::

        \psi(x) = \frac{d}{dx} \log(\Gamma(x)) = \frac{\Gamma'(x)}{\Gamma(x)}

    We represent the `n`-th derivative of the digamma function with
    `\psi(n, x)` or `psi(n, x)`.

    EXAMPLES::

        sage: # needs sage.symbolic
        sage: psi(x)
        psi(x)
        sage: psi(.5)
        -1.96351002602142
        sage: psi(3)
        -euler_gamma + 3/2
        sage: psi(1, 5)
        1/6*pi^2 - 205/144
        sage: psi(1, x)
        psi(1, x)
        sage: psi(1, x).derivative(x)
        psi(2, x)

    ::

        sage: psi(3, hold=True)                                                         # needs sage.symbolic
        psi(3)
        sage: psi(1, 5, hold=True)                                                      # needs sage.symbolic
        psi(1, 5)

    TESTS::

        sage: psi(2, x, 3)                                                              # needs sage.symbolic
        Traceback (most recent call last):
        ...
        TypeError: Symbolic function psi takes at most 2 arguments (3 given)
    """
    if not args:
        return psi1(x, **kwds)
    if len(args) > 1:
        raise TypeError("Symbolic function psi takes at most 2 arguments (%s given)" % (len(args) + 1))
    return psi2(x, args[0], **kwds)


# We have to add the wrapper function manually to the symbol_table when we have
# two functions with different number of arguments and the same name
symbol_table['functions']['psi'] = psi


def _swap_psi(a, b):
    return psi(b, a)


register_symbol(_swap_psi, {'giac': 'Psi'}, 2)


class Function_beta(GinacFunction):
    def __init__(self):
        r"""
        Return the beta function.  This is defined by

        .. MATH::

            \operatorname{B}(p,q) = \int_0^1 t^{p-1}(1-t)^{q-1} dt

        for complex or symbolic input `p` and `q`.
        Note that the order of inputs does not matter:
        `\operatorname{B}(p,q)=\operatorname{B}(q,p)`.

        GiNaC is used to compute `\operatorname{B}(p,q)`.  However, complex inputs
        are not yet handled in general.  When GiNaC raises an error on
        such inputs, we raise a :exc:`NotImplementedError`.

        If either input is 1, GiNaC returns the reciprocal of the
        other.  In other cases, GiNaC uses one of the following
        formulas:

        .. MATH::

            \operatorname{B}(p,q) = \frac{\Gamma(p)\Gamma(q)}{\Gamma(p+q)}

        or

        .. MATH::

            \operatorname{B}(p,q) = (-1)^q \operatorname{B}(1-p-q, q).


        For numerical inputs, GiNaC uses the formula

        .. MATH::

            \operatorname{B}(p,q) =  \exp[\log\Gamma(p)+\log\Gamma(q)-\log\Gamma(p+q)]


        INPUT:

        - ``p`` -- number or symbolic expression

        - ``q`` -- number or symbolic expression

        OUTPUT: number or symbolic expression (if input is symbolic)

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: beta(3, 2)
            1/12
            sage: beta(3, 1)
            1/3
            sage: beta(1/2, 1/2)
            beta(1/2, 1/2)
            sage: beta(-1, 1)
            -1
            sage: beta(-1/2, -1/2)
            0
            sage: ex = beta(x/2, 3)
            sage: set(ex.operands()) == set([1/2*x, 3])
            True
            sage: beta(.5, .5)
            3.14159265358979
            sage: beta(1, 2.0+I)
            0.400000000000000 - 0.200000000000000*I
            sage: ex = beta(3, x+I)
            sage: set(ex.operands()) == set([x+I, 3])
            True

        The result is symbolic if exact input is given::

            sage: # needs sage.symbolic
            sage: ex = beta(2, 1+5*I); ex
            beta(...
            sage: set(ex.operands()) == set([1+5*I, 2])
            True
            sage: beta(2, 2.)
            0.166666666666667
            sage: beta(I, 2.)
            -0.500000000000000 - 0.500000000000000*I
            sage: beta(2., 2)
            0.166666666666667
            sage: beta(2., I)
            -0.500000000000000 - 0.500000000000000*I

            sage: beta(x, x)._sympy_()                                                  # needs sympy sage.symbolic
            beta(x, x)

        Test pickling::

            sage: loads(dumps(beta))
            beta

        Check that :issue:`15196` is fixed::

            sage: beta(-1.3, -0.4)                                                      # needs sage.symbolic
            -4.92909641669610
        """
        GinacFunction.__init__(self, 'beta', nargs=2,
                               latex_name=r"\operatorname{B}",
                               conversions=dict(maxima='beta',
                                                mathematica='Beta',
                                                maple='Beta',
                                                sympy='beta',
                                                fricas='Beta',
                                                giac='Beta'))

    def _method_arguments(self, x, y):
        r"""
        TESTS::

            sage: RBF(beta(sin(3), sqrt(RBF(2).add_error(1e-8)/3)))  # abs tol 6e-7     # needs sage.libs.flint sage.symbolic
            [7.407662 +/- 6.17e-7]
        """
        return [x, y]


beta = Function_beta()
