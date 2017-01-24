"""
Logarithmic functions

AUTHORS:

- Yoora Yi Tenen (2012-11-16): Add documentation for :meth:`log()` (:trac:`12113`)

- Tomas Kalvoda (2015-04-01): Add :meth:`exp_polar()` (:trac:`18085`)

"""
from six.moves import range

from sage.symbolic.function import GinacFunction, BuiltinFunction
from sage.symbolic.constants import e as const_e
from sage.symbolic.constants import pi as const_pi

from sage.libs.mpmath import utils as mpmath_utils
from sage.structure.all import parent as s_parent
from sage.symbolic.expression import Expression
from sage.rings.real_double import RDF
from sage.rings.complex_double import CDF
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.rational import Rational

class Function_exp(GinacFunction):
    def __init__(self):
        r"""
        The exponential function, `\exp(x) = e^x`.

        EXAMPLES::

            sage: exp(-1)
            e^(-1)
            sage: exp(2)
            e^2
            sage: exp(2).n(100)
            7.3890560989306502272304274606
            sage: exp(x^2 + log(x))
            e^(x^2 + log(x))
            sage: exp(x^2 + log(x)).simplify()
            x*e^(x^2)
            sage: exp(2.5)
            12.1824939607035
            sage: exp(float(2.5))
            12.182493960703473
            sage: exp(RDF('2.5'))
            12.182493960703473
            sage: exp(I*pi/12)
            (1/4*I + 1/4)*sqrt(6) - (1/4*I - 1/4)*sqrt(2)

        To prevent automatic evaluation, use the ``hold`` parameter::

            sage: exp(I*pi,hold=True)
            e^(I*pi)
            sage: exp(0,hold=True)
            e^0

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: exp(0,hold=True).simplify()
            1

        ::

            sage: exp(pi*I/2)
            I
            sage: exp(pi*I)
            -1
            sage: exp(8*pi*I)
            1
            sage: exp(7*pi*I/2)
            -I

        The precision for the result is deduced from the precision of
        the input. Convert the input to a higher precision explicitly
        if a result with higher precision is desired::

            sage: t = exp(RealField(100)(2)); t
            7.3890560989306502272304274606
            sage: t.prec()
            100
            sage: exp(2).n(100)
            7.3890560989306502272304274606

        TEST::

            sage: latex(exp(x))
            e^{x}
            sage: latex(exp(sqrt(x)))
            e^{\sqrt{x}}
            sage: latex(exp)
            \exp
            sage: latex(exp(sqrt(x))^x)
            \left(e^{\sqrt{x}}\right)^{x}
            sage: latex(exp(sqrt(x)^x))
            e^{\left(\sqrt{x}^{x}\right)}
            sage: exp(x)._sympy_()
            exp(x)

        Test conjugates::

            sage: conjugate(exp(x))
            e^conjugate(x)

        Test simplifications when taking powers of exp, #7264::

            sage: var('a,b,c,II')
            (a, b, c, II)
            sage: model_exp = exp(II)**a*(b)
            sage: sol1_l={b: 5.0, a: 1.1}
            sage: model_exp.subs(sol1_l)
            5.00000000000000*(e^II)^1.10000000000000

        ::

            sage: exp(3)^II*exp(x)
            (e^3)^II*e^x
            sage: exp(x)*exp(x)
            e^(2*x)
            sage: exp(x)*exp(a)
            e^(a + x)
            sage: exp(x)*exp(a)^2
            e^(2*a + x)

        Another instance of the same problem, #7394::

            sage: 2*sqrt(e)
            2*sqrt(e)

        Check that :trac:`19918` is fixed::

            sage: exp(-x^2).subs(x=oo)
            0
            sage: exp(-x).subs(x=-oo)
            +Infinity
        """
        GinacFunction.__init__(self, "exp", latex_name=r"\exp",
                                   conversions=dict(maxima='exp', fricas='exp'))

exp = Function_exp()

class Function_log(GinacFunction):
    def __init__(self):
        r"""
        The natural logarithm of x.  See `log?` for
        more information about its behavior.

        EXAMPLES::

            sage: ln(e^2)
            2
            sage: ln(2)
            log(2)
            sage: ln(10)
            log(10)

        ::

            sage: ln(RDF(10))
            2.302585092994046
            sage: ln(2.718)
            0.999896315728952
            sage: ln(2.0)
            0.693147180559945
            sage: ln(float(-1))
            3.141592653589793j
            sage: ln(complex(-1))
            3.141592653589793j

        The ``hold`` parameter can be used to prevent automatic evaluation::

            sage: log(-1,hold=True)
            log(-1)
            sage: log(-1)
            I*pi
            sage: I.log(hold=True)
            log(I)
            sage: I.log(hold=True).simplify()
            1/2*I*pi

        TESTS::

            sage: latex(x.log())
            \log\left(x\right)
            sage: latex(log(1/4))
            \log\left(\frac{1}{4}\right)
            sage: log(x)._sympy_()
            log(x)
            sage: loads(dumps(ln(x)+1))
            log(x) + 1

        ``conjugate(log(x))==log(conjugate(x))`` unless on the branch cut which
        runs along the negative real axis.::

            sage: conjugate(log(x))
            conjugate(log(x))
            sage: var('y', domain='positive')
            y
            sage: conjugate(log(y))
            log(y)
            sage: conjugate(log(y+I))
            conjugate(log(y + I))
            sage: conjugate(log(-1))
            -I*pi
            sage: log(conjugate(-1))
            I*pi

        Check if float arguments are handled properly.::

            sage: from sage.functions.log import function_log as log
            sage: log(float(5))
            1.6094379124341003
            sage: log(float(0))
            -inf
            sage: log(float(-1))
            3.141592653589793j
            sage: log(x).subs(x=float(-1))
            3.141592653589793j
        """
        GinacFunction.__init__(self, 'log', latex_name=r'\log',
                               conversions=dict(maxima='log', fricas='log'))

    def __call__(self, *args, **kwds):
        """
        Return the logarithm of x to the given base.

        Calls the ``log`` method of the object x when computing
        the logarithm, thus allowing use of logarithm on any object
        containing a ``log`` method. In other words, log works
        on more than just real numbers.

        EXAMPLES::

            sage: log(e^2)
            2

        To change the base of the logarithm, add a second parameter::

            sage: log(1000,10)
            3

        You can use
        :class:`RDF<sage.rings.real_double.RealDoubleField_class>`,
        :class:`~sage.rings.real_mpfr.RealField` or ``n`` to get a
        numerical real approximation::

            sage: log(1024, 2)
            10
            sage: RDF(log(1024, 2))
            10.0
            sage: log(10, 4)
            log(10)/log(4)
            sage: RDF(log(10, 4))
            1.6609640474436813
            sage: log(10, 2)
            log(10)/log(2)
            sage: n(log(10, 2))
            3.32192809488736
            sage: log(10, e)
            log(10)
            sage: n(log(10, e))
            2.30258509299405

        The log function works for negative numbers, complex
        numbers, and symbolic numbers too, picking the branch
        with angle between `-pi` and `pi`::

            sage: log(-1+0*I)
            I*pi
            sage: log(CC(-1))
            3.14159265358979*I
            sage: log(-1.0)
            3.14159265358979*I

        For input zero, the following behavior occurs::

            sage: log(0)
            -Infinity
            sage: log(CC(0))
            -infinity
            sage: log(0.0)
            -infinity

        The log function also works in finite fields as long as the
        argument lies in the multiplicative group generated by the base::

            sage: F = GF(13); g = F.multiplicative_generator(); g
            2
            sage: a = F(8)
            sage: log(a,g); g^log(a,g)
            3
            8
            sage: log(a,3)
            Traceback (most recent call last):
            ...
            ValueError: No discrete log of 8 found to base 3
            sage: log(F(9), 3)
            2

        The log function also works for p-adics (see documentation for
        p-adics for more information)::

            sage: R = Zp(5); R
            5-adic Ring with capped relative precision 20
            sage: a = R(16); a
            1 + 3*5 + O(5^20)
            sage: log(a)
            3*5 + 3*5^2 + 3*5^4 + 3*5^5 + 3*5^6 + 4*5^7 + 2*5^8 + 5^9 +
            5^11 + 2*5^12 + 5^13 + 3*5^15 + 2*5^16 + 4*5^17 + 3*5^18 +
            3*5^19 + O(5^20)


        TESTS:

        Check if :trac:`10136` is fixed::

            sage: log(x).operator() is log
            True
            sage: log(x).operator() is ln
            True

            sage: log(1000, 10, base=5)
            Traceback (most recent call last):
            ...
            TypeError: Symbolic function log must be called as log(x),
            log(x, base=b) or log(x, b)
        """
        base = kwds.pop('base', None)
        if base is None:
            if len(args) == 1:
                return GinacFunction.__call__(self, *args, **kwds)
            # second argument is base
            base = args[1]
            args = args[:1]

        if len(args) != 1:
            raise TypeError("Symbolic function log must be called as "
                    "log(x), log(x, base=b) or log(x, b)")

        try:
            return args[0].log(base)
        except (AttributeError, TypeError):
            return GinacFunction.__call__(self, *args, **kwds) / \
                GinacFunction.__call__(self, base, **kwds)

ln = log = function_log = Function_log()


class Function_polylog(GinacFunction):
    def __init__(self):
        r"""
        The polylog function
        `\text{Li}_s(z) = \sum_{k=1}^{\infty} z^k / k^s`.

        This definition is valid for arbitrary complex order `s` and for
        all complex arguments `z` with `|z| < 1`; it can be extended to
        `|z| \ge 1` by the process of analytic continuation. So the
        function may have a discontinuity at `z=1` which can cause a
        `NaN` value returned for floating point arguments.

        EXAMPLES::

            sage: polylog(2.7, 0)
            0
            sage: polylog(2, 1)
            1/6*pi^2
            sage: polylog(2, -1)
            -1/12*pi^2
            sage: polylog(3, -1)
            -3/4*zeta(3)
            sage: polylog(2, I)
            I*catalan - 1/48*pi^2
            sage: polylog(4, 1/2)
            polylog(4, 1/2)
            sage: polylog(4, 0.5)
            0.517479061673899

            sage: polylog(1, x)
            -log(-x + 1)
            sage: polylog(2,x^2+1)
            dilog(x^2 + 1)

            sage: f = polylog(4, 1); f
            1/90*pi^4
            sage: f.n()
            1.08232323371114

            sage: polylog(4, 2).n()
            2.42786280675470 - 0.174371300025453*I
            sage: complex(polylog(4,2))
            (2.4278628067547032-0.17437130002545306j)
            sage: float(polylog(4,0.5))
            0.5174790616738993

            sage: z = var('z')
            sage: polylog(2,z).series(z==0, 5)
            1*z + 1/4*z^2 + 1/9*z^3 + 1/16*z^4 + Order(z^5)

            sage: loads(dumps(polylog))
            polylog

            sage: latex(polylog(5, x))
            {\rm Li}_{5}(x)
            sage: polylog(x, x)._sympy_()
            polylog(x, x)

        TESTS:

        Check if :trac:`8459` is fixed::

            sage: t = maxima(polylog(5,x)).sage(); t
            polylog(5, x)
            sage: t.operator() == polylog
            True
            sage: t.subs(x=.5).n()
            0.50840057924226...

        Check if :trac:`18386` is fixed::

            sage: polylog(2.0, 1)
            1.64493406684823
            sage: polylog(2, 1.0)
            NaN - NaN*I
            sage: polylog(2.0, 1.0)
            NaN - NaN*I
        """
        GinacFunction.__init__(self, "polylog", nargs=2)

    def _maxima_init_evaled_(self, *args):
        """
        EXAMPLES:

        These are indirect doctests for this function.::

            sage: polylog(2, x)._maxima_()
            li[2](_SAGE_VAR_x)
            sage: polylog(4, x)._maxima_()
            polylog(4,_SAGE_VAR_x)
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
        if int(n) in [1,2,3]:
            return 'li[%s](%s)'%(n, x)
        else:
            return 'polylog(%s, %s)'%(n, x)


polylog = Function_polylog()

class Function_dilog(GinacFunction):
    def __init__(self):
        r"""
        The dilogarithm function
        `\text{Li}_2(z) = \sum_{k=1}^{\infty} z^k / k^2`.

        This is simply an alias for polylog(2, z).

        EXAMPLES::

            sage: dilog(1)
            1/6*pi^2
            sage: dilog(1/2)
            1/12*pi^2 - 1/2*log(2)^2
            sage: dilog(x^2+1)
            dilog(x^2 + 1)
            sage: dilog(-1)
            -1/12*pi^2
            sage: dilog(-1.1)
            -0.890838090262283
            sage: float(dilog(1))
            1.6449340668482262
            sage: var('z')
            z
            sage: dilog(z).diff(z, 2)
            log(-z + 1)/z^2 - 1/((z - 1)*z)
            sage: dilog(z).series(z==1/2, 3)
            (1/12*pi^2 - 1/2*log(2)^2) + (-2*log(1/2))*(z - 1/2) + (2*log(1/2) + 2)*(z - 1/2)^2 + Order(1/8*(2*z - 1)^3)

            sage: latex(dilog(z))
            {\rm Li}_2\left(z\right)

        TESTS:

        ``conjugate(dilog(x))==dilog(conjugate(x))`` unless on the branch cuts
        which run along the positive real axis beginning at 1.::

            sage: conjugate(dilog(x))
            conjugate(dilog(x))
            sage: var('y',domain='positive')
            y
            sage: conjugate(dilog(y))
            conjugate(dilog(y))
            sage: conjugate(dilog(1/19))
            dilog(1/19)
            sage: conjugate(dilog(1/2*I))
            dilog(-1/2*I)
            sage: dilog(conjugate(1/2*I))
            dilog(-1/2*I)
            sage: conjugate(dilog(2))
            conjugate(dilog(2))
        """
        GinacFunction.__init__(self, 'dilog',
                conversions=dict(maxima='li[2]'))

dilog = Function_dilog()


class Function_lambert_w(BuiltinFunction):
    r"""
    The integral branches of the Lambert W function `W_n(z)`.

    This function satisfies the equation

    .. MATH::

        z = W_n(z) e^{W_n(z)}

    INPUT:

    - ``n`` - an integer. `n=0` corresponds to the principal branch.

    - ``z`` - a complex number

    If called with a single argument, that argument is ``z`` and the branch ``n`` is
    assumed to be 0 (the principal branch).

    ALGORITHM:

    Numerical evaluation is handled using the mpmath and SciPy libraries.

    REFERENCES:

    - :wikipedia:`Lambert_W_function`

    EXAMPLES:

    Evaluation of the principal branch::

        sage: lambert_w(1.0)
        0.567143290409784
        sage: lambert_w(-1).n()
        -0.318131505204764 + 1.33723570143069*I
        sage: lambert_w(-1.5 + 5*I)
        1.17418016254171 + 1.10651494102011*I

    Evaluation of other branches::

        sage: lambert_w(2, 1.0)
        -2.40158510486800 + 10.7762995161151*I

    Solutions to certain exponential equations are returned in terms of lambert_w::

        sage: S = solve(e^(5*x)+x==0, x, to_poly_solve=True)
        sage: z = S[0].rhs(); z
        -1/5*lambert_w(5)
        sage: N(z)
        -0.265344933048440

    Check the defining equation numerically at `z=5`::

        sage: N(lambert_w(5)*exp(lambert_w(5)) - 5)
        0.000000000000000

    There are several special values of the principal branch which
    are automatically simplified::

        sage: lambert_w(0)
        0
        sage: lambert_w(e)
        1
        sage: lambert_w(-1/e)
        -1

    Integration (of the principal branch) is evaluated using Maxima::

        sage: integrate(lambert_w(x), x)
        (lambert_w(x)^2 - lambert_w(x) + 1)*x/lambert_w(x)
        sage: integrate(lambert_w(x), x, 0, 1)
        (lambert_w(1)^2 - lambert_w(1) + 1)/lambert_w(1) - 1
        sage: integrate(lambert_w(x), x, 0, 1.0)
        0.3303661247616807

    Warning: The integral of a non-principal branch is not implemented,
    neither is numerical integration using GSL. The :meth:`numerical_integral`
    function does work if you pass a lambda function::

        sage: numerical_integral(lambda x: lambert_w(x), 0, 1)
        (0.33036612476168054, 3.667800782666048e-15)
    """

    def __init__(self):
        r"""
        See the docstring for :meth:`Function_lambert_w`.

        EXAMPLES::

            sage: lambert_w(0, 1.0)
            0.567143290409784
            sage: lambert_w(x, x)._sympy_()
            LambertW(x, x)
        """
        BuiltinFunction.__init__(self, "lambert_w", nargs=2,
                                 conversions={'mathematica': 'ProductLog',
                                              'maple': 'LambertW',
                                              'matlab': 'lambertw',
                                              'maxima': 'generalized_lambert_w',
                                              'sympy': 'LambertW'})

    def __call__(self, *args, **kwds):
        r"""
        Custom call method allows the user to pass one argument or two. If
        one argument is passed, we assume it is ``z`` and that ``n=0``.

        EXAMPLES::

            sage: lambert_w(1)
            lambert_w(1)
            sage: lambert_w(1, 2)
            lambert_w(1, 2)
        """
        if len(args) == 2:
            return BuiltinFunction.__call__(self, *args, **kwds)
        elif len(args) == 1:
            return BuiltinFunction.__call__(self, 0, args[0], **kwds)
        else:
            raise TypeError("lambert_w takes either one or two arguments.")

    def _eval_(self, n, z):
        """
        EXAMPLES::

            sage: lambert_w(6.0)
            1.43240477589830
            sage: lambert_w(1)
            lambert_w(1)
            sage: lambert_w(x+1)
            lambert_w(x + 1)

        There are three special values which are automatically simplified::

            sage: lambert_w(0)
            0
            sage: lambert_w(e)
            1
            sage: lambert_w(-1/e)
            -1
            sage: lambert_w(SR(0))
            0

        The special values only hold on the principal branch::

            sage: lambert_w(1,e)
            lambert_w(1, e)
            sage: lambert_w(1, e.n())
            -0.532092121986380 + 4.59715801330257*I

        TESTS:

        When automatic simplication occurs, the parent of the output value should be
        either the same as the parent of the input, or a Sage type::

            sage: parent(lambert_w(int(0)))
            <type 'int'>
            sage: parent(lambert_w(Integer(0)))
            Integer Ring
            sage: parent(lambert_w(e))
            Symbolic Ring
        """
        if not isinstance(z, Expression):
            if n == 0 and z == 0:
                return s_parent(z)(0)
        elif n == 0:
            if z.is_trivial_zero():
                return s_parent(z)(Integer(0))
            elif (z-const_e).is_trivial_zero():
                return s_parent(z)(Integer(1))
            elif (z+1/const_e).is_trivial_zero():
                return s_parent(z)(Integer(-1))

    def _evalf_(self, n, z, parent=None, algorithm=None):
        """
        EXAMPLES::

            sage: N(lambert_w(1))
            0.567143290409784
            sage: lambert_w(RealField(100)(1))
            0.56714329040978387299996866221

        SciPy is used to evaluate for float, RDF, and CDF inputs::

            sage: lambert_w(RDF(1))
            0.5671432904097838
            sage: lambert_w(float(1))
            0.5671432904097838
            sage: lambert_w(CDF(1))
            0.5671432904097838
            sage: lambert_w(complex(1))
            (0.5671432904097838+0j)
            sage: lambert_w(RDF(-1))  # abs tol 2e-16
            -0.31813150520476413 + 1.3372357014306895*I
            sage: lambert_w(float(-1))  # abs tol 2e-16
            (-0.31813150520476413+1.3372357014306895j)
        """
        R = parent or s_parent(z)
        if R is float or R is RDF:
            from scipy.special import lambertw
            res = lambertw(z, n)
            # SciPy always returns a complex value, make it real if possible
            if not res.imag:
                return R(res.real)
            elif R is float:
                return complex(res)
            else:
                return CDF(res)
        elif R is complex or R is CDF:
            from scipy.special import lambertw
            return R(lambertw(z, n))
        else:
            import mpmath
            return mpmath_utils.call(mpmath.lambertw, z, n, parent=R)

    def _derivative_(self, n, z, diff_param=None):
        """
        The derivative of `W_n(x)` is `W_n(x)/(x \cdot W_n(x) + x)`.

        EXAMPLES::

            sage: x = var('x')
            sage: derivative(lambert_w(x), x)
            lambert_w(x)/(x*lambert_w(x) + x)

            sage: derivative(lambert_w(2, exp(x)), x)
            e^x*lambert_w(2, e^x)/(e^x*lambert_w(2, e^x) + e^x)

        TESTS:

        Differentiation in the first parameter raises an error :trac:`14788`::

            sage: n = var('n')
            sage: lambert_w(n, x).diff(n)
            Traceback (most recent call last):
            ...
            ValueError: cannot differentiate lambert_w in the first parameter
        """
        if diff_param == 0:
            raise ValueError("cannot differentiate lambert_w in the first parameter")

        return lambert_w(n, z)/(z*lambert_w(n, z)+z)

    def _maxima_init_evaled_(self, n, z):
        """
        EXAMPLES:

        These are indirect doctests for this function.::

            sage: lambert_w(0, x)._maxima_()
            lambert_w(_SAGE_VAR_x)
            sage: lambert_w(1, x)._maxima_()
            generalized_lambert_w(1,_SAGE_VAR_x)

        TESTS::

            sage: lambert_w(x)._maxima_()._sage_()
            lambert_w(x)
            sage: lambert_w(2, x)._maxima_()._sage_()
            lambert_w(2, x)
        """
        if isinstance(z, str):
            maxima_z = z
        elif hasattr(z, '_maxima_init_'):
            maxima_z = z._maxima_init_()
        else:
            maxima_z = str(z)
        if n == 0:
            return "lambert_w(%s)" % maxima_z
        else:
            return "generalized_lambert_w(%s,%s)" % (n, maxima_z)

    def _print_(self, n, z):
        """
        Custom _print_ method to avoid printing the branch number if
        it is zero.

        EXAMPLES::

            sage: lambert_w(1)
            lambert_w(1)
            sage: lambert_w(0,x)
            lambert_w(x)
        """
        if n == 0:
            return "lambert_w(%s)" % z
        else:
            return "lambert_w(%s, %s)" % (n, z)

    def _print_latex_(self, n, z):
        """
        Custom _print_latex_ method to avoid printing the branch
        number if it is zero.

        EXAMPLES::

            sage: latex(lambert_w(1))
            \operatorname{W}({1})
            sage: latex(lambert_w(0,x))
            \operatorname{W}({x})
            sage: latex(lambert_w(1,x))
            \operatorname{W_{1}}({x})
            sage: latex(lambert_w(1,x+exp(x)))
            \operatorname{W_{1}}({x + e^{x}})
        """
        if n == 0:
            return r"\operatorname{W}({%s})" % z._latex_()
        else:
            return r"\operatorname{W_{%s}}({%s})" % (n, z._latex_())

lambert_w = Function_lambert_w()


class Function_exp_polar(BuiltinFunction):
    def __init__(self):
        r"""
        Representation of a complex number in a polar form.

        INPUT:

        - ``z`` - a complex number `z = a + ib`.

        OUTPUT:

        A complex number with modulus `\exp(a)` and argument `b`.

        If `-\pi < b \leq \pi` then `\operatorname{exp\_polar}(z)=\exp(z)`.
        For other values of `b` the function is left unevaluated.

        EXAMPLES:

        The following expressions are evaluated using the exponential
        function::

            sage: exp_polar(pi*I/2)
            I
            sage: x = var('x', domain='real')
            sage: exp_polar(-1/2*I*pi + x)
            e^(-1/2*I*pi + x)

        The function is left unevaluated when the imaginary part of the
        input `z` does not satisfy `-\pi < \Im(z) \leq \pi`::

            sage: exp_polar(2*pi*I)
            exp_polar(2*I*pi)
            sage: exp_polar(-4*pi*I)
            exp_polar(-4*I*pi)

        This fixes :trac:`18085`::

            sage: integrate(1/sqrt(1+x^3),x,algorithm='sympy')
            1/3*x*hypergeometric((1/3, 1/2), (4/3,), -x^3)*gamma(1/3)/gamma(4/3)

        SEEALSO:

            `Examples in Sympy documentation <http://docs.sympy.org/latest/modules/functions/special.html?highlight=exp_polar>`_,
            `Sympy source code of exp_polar <http://docs.sympy.org/0.7.4/_modules/sympy/functions/elementary/exponential.html>`_

        REFERENCES:

            :wikipedia:`Complex_number#Polar_form`
        """
        BuiltinFunction.__init__(self, "exp_polar",
                                latex_name=r"\operatorname{exp\_polar}",
                                conversions=dict(sympy='exp_polar'))

    def _evalf_(self, z, parent=None, algorithm=None):
        r"""
        EXAMPLES:

        If the imaginary part of `z` obeys `-\pi < z \leq \pi`, then
        `\operatorname{exp\_polar}(z)` is evaluated as `\exp(z)`::

            sage: exp_polar(1.0 + 2.0*I)
            -1.13120438375681 + 2.47172667200482*I

        If the imaginary part of `z` is outside of that interval the
        expression is left unevaluated::

            sage: exp_polar(-5.0 + 8.0*I)
            exp_polar(-5.00000000000000 + 8.00000000000000*I)

        An attempt to numerically evaluate such an expression raises an error::

            sage: exp_polar(-5.0 + 8.0*I).n()
            Traceback (most recent call last):
            ...
            ValueError: invalid attempt to numerically evaluate exp_polar()

        """
        from sage.functions.other import imag

        if (not isinstance(z, Expression)
            and bool(-const_pi < imag(z) <= const_pi)):
            return exp(z)
        else:
            raise ValueError("invalid attempt to numerically evaluate exp_polar()")

    def _eval_(self, z):
        """
        EXAMPLES::

            sage: exp_polar(3*I*pi)
            exp_polar(3*I*pi)
            sage: x = var('x', domain='real')
            sage: exp_polar(4*I*pi + x)
            exp_polar(4*I*pi + x)

        """
        if (isinstance(z, Expression)
            and bool(-const_pi < z.imag_part() <= const_pi)):
            return exp(z)
        else:
            return None

exp_polar = Function_exp_polar()


class Function_harmonic_number_generalized(BuiltinFunction):
    r"""
    Harmonic and generalized harmonic number functions,
    defined by:

    .. MATH::

        H_{n}=H_{n,1}=\sum_{k=1}^n\frac1k
        
        H_{n,m}=\sum_{k=1}^n\frac1{k^m}

    They are also well-defined for complex argument, through:

    .. MATH::

        H_{s}=\int_0^1\frac{1-x^s}{1-x}

        H_{s,m}=\zeta(m)-\zeta(m,s-1)

    If called with a single argument, that argument is ``s`` and ``m`` is
    assumed to be 1 (the normal harmonic numbers ``H_s``).

    ALGORITHM:

    Numerical evaluation is handled using the mpmath and FLINT libraries.

    REFERENCES:

    - :wikipedia:`/Harmonic_number`

    EXAMPLES:

    Evaluation of integer, rational, or complex argument::

        sage: harmonic_number(5)
        137/60
        sage: harmonic_number(3,3)
        251/216
        sage: harmonic_number(5/2)
        -2*log(2) + 46/15
        sage: harmonic_number(3.,3)
        zeta(3) - 0.0400198661225573
        sage: harmonic_number(3.,3.)
        1.16203703703704
        sage: harmonic_number(3,3).n(200)
        1.16203703703703703703703...
        sage: harmonic_number(1+I,5)
        harmonic_number(I + 1, 5)
        sage: harmonic_number(5,1.+I)
        1.57436810798989 - 1.06194728851357*I

    Solutions to certain sums are returned in terms of harmonic numbers::

        sage: k=var('k')
        sage: sum(1/k^7,k,1,x)
        harmonic_number(x, 7)

    Check the defining integral at a random integer::

        sage: n=randint(10,100)
        sage: bool(SR(integrate((1-x^n)/(1-x),x,0,1)) == harmonic_number(n))
        True

    There are several special values which are automatically simplified::

        sage: harmonic_number(0)
        0
        sage: harmonic_number(1)
        1
        sage: harmonic_number(x,1)
        harmonic_number(x)

    Arguments are swapped with respect to the same functions in
    Maxima::

        sage: maxima(harmonic_number(x,2)) # maxima expect interface
        gen_harmonic_number(2,_SAGE_VAR_x)
        sage: from sage.calculus.calculus import symbolic_expression_from_maxima_string as sefms
        sage: sefms('gen_harmonic_number(3,x)')
        harmonic_number(x, 3)
        sage: from sage.interfaces.maxima_lib import maxima_lib, max_to_sr
        sage: c=maxima_lib(harmonic_number(x,2)); c
        gen_harmonic_number(2,_SAGE_VAR_x)
        sage: max_to_sr(c.ecl())
        harmonic_number(x, 2)
    """

    def __init__(self):
        r"""
        EXAMPLES::

            sage: loads(dumps(harmonic_number(x,5)))
            harmonic_number(x, 5)
            sage: harmonic_number(x, x)._sympy_()
            harmonic(x, x)
        """
        BuiltinFunction.__init__(self, "harmonic_number", nargs=2,
                conversions={'sympy':'harmonic'})

    def __call__(self, z, m=1, **kwds):
        r"""
        Custom call method allows the user to pass one argument or two. If
        one argument is passed, we assume it is ``z`` and that ``m=1``.

        EXAMPLES::

            sage: harmonic_number(x)
            harmonic_number(x)
            sage: harmonic_number(x,1)
            harmonic_number(x)
            sage: harmonic_number(x,2)
            harmonic_number(x, 2)
        """
        return BuiltinFunction.__call__(self, z, m, **kwds)

    def _eval_(self, z, m):
        """
        EXAMPLES::

            sage: harmonic_number(x,0)
            x
            sage: harmonic_number(x,1)
            harmonic_number(x)
            sage: harmonic_number(5)
            137/60
            sage: harmonic_number(3,3)
            251/216
            sage: harmonic_number(3,3).n() # this goes from rational to float
            1.16203703703704
            sage: harmonic_number(3,3.) # the following uses zeta functions
            1.16203703703704
            sage: harmonic_number(3.,3)
            zeta(3) - 0.0400198661225573
            sage: harmonic_number(0.1,5)
            zeta(5) - 0.650300133161038
            sage: harmonic_number(0.1,5).n()
            0.386627621982332
            sage: harmonic_number(3,5/2)
            1/27*sqrt(3) + 1/8*sqrt(2) + 1
        """
        if m == 0:
            return z
        elif m == 1:
            return harmonic_m1._eval_(z)

        if z in ZZ and z >= 0:
            return sum(1 / (k ** m) for k in range(1, z + 1))

    def _evalf_(self, z, m, parent=None, algorithm=None):
        """
        EXAMPLES::

            sage: harmonic_number(3.,3)
            zeta(3) - 0.0400198661225573
            sage: harmonic_number(3.,3.)
            1.16203703703704
            sage: harmonic_number(3,3).n(200)
            1.16203703703703703703703...
            sage: harmonic_number(5,I).n()
            2.36889632899995 - 3.51181956521611*I
        """
        if m == 0:
            if parent is None:
                return z
            return parent(z)
        elif m == 1:
            return harmonic_m1._evalf_(z, parent, algorithm)

        from sage.functions.transcendental import zeta, hurwitz_zeta
        return zeta(m) - hurwitz_zeta(m,z+1)

    def _maxima_init_evaled_(self, n, z):
        """
        EXAMPLES:

            sage: maxima_calculus(harmonic_number(x,2))
            gen_harmonic_number(2,_SAGE_VAR_x)
            sage: maxima_calculus(harmonic_number(3,harmonic_number(x,3),hold=True))
            1/3^gen_harmonic_number(3,_SAGE_VAR_x)+1/2^gen_harmonic_number(3,_SAGE_VAR_x)+1
        """
        if isinstance(n,str):
            maxima_n=n
        elif hasattr(n,'_maxima_init_'):
            maxima_n=n._maxima_init_()
        else:
            maxima_n=str(n)
        if isinstance(z,str):
            maxima_z=z
        elif hasattr(z,'_maxima_init_'):
            maxima_z=z._maxima_init_()
        else:
            maxima_z=str(z)
        return "gen_harmonic_number(%s,%s)" % (maxima_z, maxima_n) # swap arguments

    def _derivative_(self, n, m, diff_param=None):
        """
        The derivative of `H_{n,m}`.

        EXAMPLES::

            sage: k,m,n = var('k,m,n')
            sage: sum(1/k, k, 1, x).diff(x)
            1/6*pi^2 - harmonic_number(x, 2)
            sage: harmonic_number(x, 1).diff(x)
            1/6*pi^2 - harmonic_number(x, 2)
            sage: harmonic_number(n, m).diff(n)
            -m*(harmonic_number(n, m + 1) - zeta(m + 1))
            sage: harmonic_number(n, m).diff(m)
            Traceback (most recent call last):
            ...
            ValueError: cannot differentiate harmonic_number in the second parameter
        """
        from sage.functions.transcendental import zeta
        if diff_param == 1:
            raise ValueError("cannot differentiate harmonic_number in the second parameter")
        if m==1:
            return harmonic_m1(n).diff()
        else:
            return m*(zeta(m+1) - harmonic_number(n, m+1))

    def _print_(self, z, m):
        """
        EXAMPLES::

            sage: harmonic_number(x)
            harmonic_number(x)
            sage: harmonic_number(x,2)
            harmonic_number(x, 2)
        """
        if m == 1:
            return "harmonic_number(%s)" % z
        else:
            return "harmonic_number(%s, %s)" % (z, m)

    def _print_latex_(self, z, m):
        """
        EXAMPLES::

            sage: latex(harmonic_number(x))
            H_{x}
            sage: latex(harmonic_number(x,2))
            H_{{x},{2}}
        """
        if m == 1:
            return r"H_{%s}" % z
        else:
            return r"H_{{%s},{%s}}" % (z, m)

harmonic_number = Function_harmonic_number_generalized()

def _swap_harmonic(a,b): return harmonic_number(b,a)

from sage.libs.pynac.pynac import register_symbol

register_symbol(_swap_harmonic,{'maxima':'gen_harmonic_number'})
register_symbol(_swap_harmonic,{'maple':'harmonic'})

class Function_harmonic_number(BuiltinFunction):
    r"""
    Harmonic number function, defined by:

    .. MATH::

        H_{n}=H_{n,1}=\sum_{k=1}^n\frac1k

        H_{s}=\int_0^1\frac{1-x^s}{1-x}

    See the docstring for :meth:`Function_harmonic_number_generalized`.

    This class exists as callback for ``harmonic_number`` returned by Maxima.
    """

    def __init__(self):
        r"""
        EXAMPLES::

            sage: k=var('k')
            sage: loads(dumps(sum(1/k,k,1,x)))
            harmonic_number(x)
            sage: harmonic_number(x)._sympy_()
            harmonic(x)
        """
        BuiltinFunction.__init__(self, "harmonic_number", nargs=1,
                                 conversions={'mathematica':'HarmonicNumber',
                                              'maple':'harmonic',
                                              'maxima':'harmonic_number',
                                              'sympy':'harmonic'})

    def _eval_(self, z, **kwds):
        """
        EXAMPLES::

            sage: harmonic_number(0)
            0
            sage: harmonic_number(1)
            1
            sage: harmonic_number(20)
            55835135/15519504
            sage: harmonic_number(5/2)
            -2*log(2) + 46/15
            sage: harmonic_number(2*x)
            harmonic_number(2*x)
        """
        if z in ZZ:
            if z == 0:
                return Integer(0)
            elif z == 1:
                return Integer(1)
            elif z > 1:
                import sage.libs.flint.arith as flint_arith
                return flint_arith.harmonic_number(z)
        elif z in QQ:
            from sage.functions.other import psi1
            return psi1(z+1) - psi1(1)

    def _evalf_(self, z, parent=None, algorithm='mpmath'):
        """
        EXAMPLES::

            sage: harmonic_number(20).n() # this goes from rational to float
            3.59773965714368
            sage: harmonic_number(20).n(200)
            3.59773965714368191148376906...
            sage: harmonic_number(20.) # this computes the integral with mpmath
            3.59773965714368
            sage: harmonic_number(1.0*I)
            0.671865985524010 + 1.07667404746858*I
        """
        from sage.libs.mpmath import utils as mpmath_utils
        import mpmath
        return mpmath_utils.call(mpmath.harmonic, z, parent=parent)

    def _derivative_(self, z, diff_param=None):
        """
        The derivative of `H_x`.

        EXAMPLES::

            sage: k=var('k')
            sage: sum(1/k,k,1,x).diff(x)
            1/6*pi^2 - harmonic_number(x, 2)
        """
        from sage.functions.transcendental import zeta
        return zeta(2)-harmonic_number(z,2)

    def _print_latex_(self, z):
        """
        EXAMPLES::

            sage: k=var('k')
            sage: latex(sum(1/k,k,1,x))
            H_{x}
        """
        return r"H_{%s}" % z

harmonic_m1 = Function_harmonic_number()
