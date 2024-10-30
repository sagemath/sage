# sage.doctest: optional - sympy
"""
SymPy --> Sage conversion

The file consists of ``_sage_()`` methods that are added lazily to
the respective SymPy objects. Any call of the ``_sympy_()`` method
of a symbolic expression will trigger the addition. See
:class:`sage.symbolic.expression_conversion.SymPyConverter` for the
conversion to SymPy.

Only ``Function`` objects where the names differ need their own ``_sage()_``
method. There are several functions with differing name that have an alias
in Sage that is the same as the name in SymPy, so no explicit translation
is needed for them::

    sage: from sympy import Symbol, Si, Ci, Shi, Chi, sign
    sage: sx = Symbol('x')
    sage: assert sin_integral(x)._sympy_() == Si(sx)
    sage: assert sin_integral(x) == Si(sx)._sage_()
    sage: assert sinh_integral(x)._sympy_() == Shi(sx)
    sage: assert sinh_integral(x) == Shi(sx)._sage_()
    sage: assert cos_integral(x)._sympy_() == Ci(sx)
    sage: assert cos_integral(x) == Ci(sx)._sage_()
    sage: assert cosh_integral(x)._sympy_() == Chi(sx)
    sage: assert cosh_integral(x) == Chi(sx)._sage_()
    sage: assert sgn(x)._sympy_() == sign(sx)
    sage: assert sgn(x) == sign(sx)._sage_()

TESTS:

Check that :issue:`24212` is fixed::

    sage: integrate(sin(x^2), x, algorithm='sympy')
    3/8*sqrt(2)*sqrt(pi)*fresnel_sin(sqrt(2)*x/sqrt(pi))*gamma(3/4)/gamma(7/4)

Test that conversion of symbolic functions with latex names works (:issue:`31047`)::

    sage: var('phi')
    phi
    sage: function('Cp', latex_name='C_+')
    Cp
    sage: test = Cp(phi)._sympy_()._sage_()
    sage: test.operator() == Cp
    True
    sage: test.operator()._latex_() == 'C_+'
    True

AUTHORS:

- Ralf Stephan (2017-10)
"""
################################################################
#   Distributed under GNU GPL3, see www.gnu.org
################################################################
from sage.misc.misc import run_once


# ################         numbers and constants      ##############

def _sympysage_float(self):
    """
    EXAMPLES::

        sage: from sympy.core.numbers import RealNumber as RN
        sage: assert SR(-1.34)._sympy_() == RN('-1.34')
        sage: assert SR(-1.34) == RN('-1.34')._sage_()
    """
    from sage.rings.real_mpfr import create_RealNumber
    return create_RealNumber(str(self))


def _sympysage_integer_ring(self):
    r"""
    EXAMPLES::

        sage: import sympy
        sage: sympy.ZZ._sage_()
        Integer Ring
    """
    from sage.rings.integer_ring import ZZ
    return ZZ


def _sympysage_integer(self):
    """
    EXAMPLES::

        sage: from sympy.core.numbers import Integer as SympyInt
        sage: assert SR(2)._sympy_() == SympyInt(int(2))
        sage: assert SR(2) == SympyInt(int(2))._sage_()
        sage: type(SympyInt(int(2))._sage_())
        <class 'sage.rings.integer.Integer'>
    """
    from sage.rings.integer import Integer
    return Integer(self.p)


def _sympysage_rational(self):
    """
    EXAMPLES::

        sage: from sympy.core.numbers import Rational
        sage: assert SR(-5/7)._sympy_() == Rational(int(-5),int(7))
        sage: assert SR(-5/7) == Rational(int(-5),int(7))._sage_()
    """
    from sage.rings.integer import Integer
    from sage.rings.rational import Rational
    return Rational((Integer(self.p), Integer(self.q)))


def _sympysage_rational_field(self):
    r"""
    EXAMPLES::

        sage: import sympy
        sage: sympy.QQ._sage_()
        Rational Field
    """
    from sage.rings.rational_field import QQ
    return QQ


def _sympysage_real_interval(self):
    r"""
    EXAMPLES::

        sage: from sage.interfaces.sympy import sympy_init
        sage: sympy_init()

        sage: from sympy import CRootOf
        sage: from sympy.abc import x
        sage: root = CRootOf(x**3 - x^2 - x - 1, 0)
        sage: interval = root._get_interval()
        sage: interval._sage_()
        2.?
        sage: interval._sage_().parent()
        Real Interval Field with 1024 bits of precision
    """
    # NOTE: this is a very approximate conversion as we do not consider any
    # potential issue with precision
    # Just to be (a little bit) safe, we set it to 1024
    from sage.rings.real_mpfi import RealIntervalField
    RIF = RealIntervalField(1024)
    # NOTE: we call fraction_field since sympy stores mpq even
    # for integral entries
    domain = self.dom._sage_().fraction_field()
    return RIF(domain(self.a)).union(RIF(domain(self.b)))


def _sympysage_complex_interval(self):
    r"""
    EXAMPLES::

        sage: from sage.interfaces.sympy import sympy_init
        sage: sympy_init()

        sage: from sympy import CRootOf
        sage: from sympy.abc import x
        sage: root = CRootOf(x**10 - 2*x + 3, 9)
        sage: interval = root._get_interval()
        sage: interval._sage_()
        0.1? + 1.2?*I
        sage: interval._sage_().parent()
        Complex Interval Field with 1024 bits of precision
    """
    # NOTE: this is a very approximate conversion as we do not consider any
    # potential issue with precision
    # Just to be (a little bit) safe, we set it to 1024
    from sage.rings.complex_interval_field import ComplexIntervalField
    CIF = ComplexIntervalField(1024)
    # NOTE: we call fraction_field since sympy stores mpq even
    # for integral entries
    domain = self.dom._sage_().fraction_field()
    return CIF(domain(self.ax), domain(self.ay)).union(CIF(domain(self.bx), domain(self.by)))


def _sympysage_polynomial_ring(self):
    r"""
    EXAMPLES::

        sage: from sage.interfaces.sympy import sympy_init
        sage: sympy_init()

        sage: import sympy
        sage: ZZx = sympy.PolynomialRing(sympy.ZZ, 'x')
        sage: ZZx._sage_()
        Univariate Polynomial Ring in x over Integer Ring

        sage: ZZxy = sympy.PolynomialRing(ZZx, 'y')
        sage: ZZxy._sage_()
        Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Integer Ring
    """
    base_ring = self.domain._sage_()
    variables = ','.join(map(str, self.gens))
    return base_ring[variables]


def _sympysage_polynomial(self):
    r"""
    EXAMPLES::

        sage: from sage.interfaces.sympy import sympy_init
        sage: sympy_init()

        sage: import sympy
        sage: from sympy.abc import x, y
        sage: p = sympy.Poly(x*(x**2 + x - 1)**2)
        sage: p._sage_()
        x^5 + 2*x^4 - x^3 - 2*x^2 + x
        sage: p._sage_().parent()
        Univariate Polynomial Ring in x over Integer Ring

        sage: p = sympy.Poly(y*x**2 + x*y + 1)
        sage: p._sage_()
        x^2*y + x*y + 1
        sage: p._sage_().parent()
        Multivariate Polynomial Ring in x, y over Integer Ring

        sage: p = sympy.Poly(y*x**2 + x*y + 1, x)
        sage: p._sage_()
        y*x^2 + y*x + 1
        sage: p._sage_().parent()
        Univariate Polynomial Ring in x over Univariate Polynomial Ring in y over Integer Ring
    """
    base_ring = self.domain._sage_()
    variables = ','.join(map(str, self.gens))
    R = base_ring[variables]
    return R.sum(base_ring(coeff) * R.monomial(*exp) for exp, coeff in self.rep.terms(order=None))


def _sympysage_pinfty(self):
    """
    EXAMPLES::

        sage: from sympy.core.numbers import oo as sinf
        sage: assert SR(oo)._sympy_() == sinf
        sage: assert SR(oo) == sinf._sage_()
    """
    from sage.rings.infinity import PlusInfinity
    return PlusInfinity()


def _sympysage_ninfty(self):
    """
    EXAMPLES::

        sage: from sympy.core.numbers import oo as sinf
        sage: assert SR(-oo)._sympy_() == -sinf
        sage: assert SR(-oo) == (-sinf)._sage_()
    """
    from sage.rings.infinity import MinusInfinity
    return MinusInfinity()


def _sympysage_uinfty(self):
    """
    EXAMPLES::

        sage: from sympy.core.numbers import zoo
        sage: assert unsigned_infinity._sympy_() == zoo
        sage: assert unsigned_infinity == zoo._sage_()
    """
    from sage.rings.infinity import unsigned_infinity
    return unsigned_infinity


def _sympysage_nan(self):
    """
    EXAMPLES::

        sage: from sympy.core.numbers import nan as snan
        sage: assert NaN._sympy_() == snan
        sage: assert NaN == snan._sage_()
    """
    from sage.symbolic.constants import NaN
    return NaN


def _sympysage_e(self):
    """
    EXAMPLES::

        sage: from sympy.core.numbers import E
        sage: assert e._sympy_() == E
        sage: assert e == E._sage_()
    """
    from sage.symbolic.constants import e
    return e


def _sympysage_pi(self):
    """
    EXAMPLES::

        sage: from sympy.core.numbers import pi as spi
        sage: assert pi._sympy_() == spi
        sage: assert pi == spi._sage_()
    """
    from sage.symbolic.constants import pi
    return pi


def _sympysage_golden_ratio(self):
    """
    EXAMPLES::

        sage: from sympy.core.singleton import S
        sage: assert golden_ratio._sympy_() == S.GoldenRatio
        sage: assert golden_ratio == S.GoldenRatio._sage_()
    """
    from sage.symbolic.constants import golden_ratio
    return golden_ratio


def _sympysage_eulerg(self):
    """
    EXAMPLES::

        sage: from sympy.core.singleton import S
        sage: assert euler_gamma._sympy_() == S.EulerGamma
        sage: assert euler_gamma == S.EulerGamma._sage_()
    """
    from sage.symbolic.constants import euler_gamma
    return euler_gamma


def _sympysage_catalan(self):
    """
    EXAMPLES::

        sage: from sympy.core.singleton import S
        sage: assert catalan._sympy_() == S.Catalan
        sage: assert catalan == S.Catalan._sage_()
    """
    from sage.symbolic.constants import catalan
    return catalan


def _sympysage_i(self):
    """
    EXAMPLES::

        sage: from sympy.core.singleton import S
        sage: assert I._sympy_() == S.ImaginaryUnit
        sage: assert I == S.ImaginaryUnit._sage_()
    """
    from sage.symbolic.constants import I
    return I


# #################       basic operators         ##############

def _sympysage_add(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol
        sage: from sympy.core.singleton import S
        sage: assert (x-pi+1)._sympy_() == Symbol('x')-S.Pi+1
        sage: assert x-pi+1 == (Symbol('x')-S.Pi+1)._sage_()
    """
    s = 0
    for x in self.args:
        s += x._sage_()
    return s


def _sympysage_mul(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol
        sage: from sympy.core.singleton import S
        sage: assert (-x*pi*5)._sympy_() == -Symbol('x')*S.Pi*5
        sage: assert -x*pi*5 == (-Symbol('x')*S.Pi*5)._sage_()
    """
    s = 1
    for x in self.args:
        s *= x._sage_()
    return s


def _sympysage_pow(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol
        sage: from sympy.core.singleton import S
        sage: assert (x^pi^5)._sympy_() == Symbol('x')**S.Pi**5
        sage: assert x^pi^5 == (Symbol('x')**S.Pi**5)._sage_()
    """
    return self.args[0]._sage_()**self.args[1]._sage_()


def _sympysage_symbol(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol
        sage: assert x._sympy_() == Symbol('x')
        sage: assert x == Symbol('x')._sage_()
    """
    from sage.symbolic.ring import SR
    try:
        return SR.var(self.name)
    except ValueError:
        # sympy sometimes returns dummy variables
        # with name = 'None', str rep = '_None'
        # in particular in inverse Laplace and inverse Mellin transforms
        return SR.var(str(self))


def _sympysage_Subs(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol
        sage: from sympy.core.singleton import S
    """
    args = self.args
    substi = {args[1][i]._sage_(): args[2][i]._sage_()
              for i in range(len(args[1]))}
    return args[0]._sage_().subs(substi)


# #############       functions       ###############

def _sympysage_function_by_name(fname):
    """
    Given a sympy function with name ``fname`` find the corresponding
    sage function or create a new one with the given name.

    EXAMPLES::

        sage: from sympy import Function
        sage: f = function('f')
        sage: F = Function('f')
        sage: assert f._sympy_() == F
        sage: assert f == F._sage_()
    """
    from sage.functions import all as sagefuncs
    func = getattr(sagefuncs, fname, None)
    # In the case the function is not known in sage:
    if func is None:
        import sympy
        if getattr(sympy, fname, None) is None:
            # symbolic function
            from sage.symbolic.expression import symbol_table
            func = symbol_table['functions'].get(fname)
            if func is None:
                from sage.calculus.var import function
                return function(fname)

        else:
            # the function defined in sympy is not known in sage
            raise AttributeError
    return func


# the convoluted class structure with metaclasses and stuff sympy uses
# to implement undefined functions makes things a bit harder for us
# here
class UndefSageHelper:
    """
    Helper class to convert sympy function objects to sage functions

    EXAMPLES::

        sage: from sympy import Function
        sage: f = function('f')
        sage: F = Function('f')
        sage: assert f._sympy_() == F
        sage: assert f == F._sage_()
    """
    def __get__(self, ins, typ):
        if ins is None:
            return lambda: _sympysage_function_by_name(typ.__name__)
        else:
            args = [arg._sage_() for arg in ins.args]
            return lambda: _sympysage_function_by_name(ins.__class__.__name__)(*args)


def _sympysage_function(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, Function, sin as Sin
        sage: assert sin(x)._sympy_() == Sin(Symbol('x'))
        sage: assert sin(x) == Sin(Symbol('x'))._sage_()

        sage: f = function('f')
        sage: F = Function('f')
        sage: assert f(x)._sympy_() == F(x)
        sage: assert f(x) == F(x)._sage_()
        sage: assert f(x+3)._sympy_() == F(x+3)
        sage: assert f(x+3) == F(x+3)._sage_()
        sage: assert (3*f(x))._sympy_() == 3*F(x)
        sage: assert 3*f(x) == (3*F(x))._sage_()

    Test that functions unknown to Sage raise an exception::

        sage: from sympy.functions.combinatorial.numbers import lucas
        sage: lucas(Symbol('x'))._sage_()
        Traceback (most recent call last):
        ...
        AttributeError...
    """
    fname = self.func.__name__
    func = _sympysage_function_by_name(fname)
    args = [arg._sage_() for arg in self.args]

    return func(*args)


def _sympysage_integral(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, Integral
        sage: sx = Symbol('x')
        sage: assert integral(x, x, hold=True)._sympy_() == Integral(sx, sx)
        sage: assert integral(x, x, hold=True) == Integral(sx, sx)._sage_()
        sage: assert integral(x, x, 0, 1, hold=True)._sympy_() == Integral(sx, (sx,0,1))
        sage: assert integral(x, x, 0, 1, hold=True) == Integral(sx, (sx,0,1))._sage_()
    """
    from sage.misc.functional import integral
    f, limits = self.function._sage_(), list(self.limits)
    for limit in limits:
        if len(limit) == 1:
            x = limit[0]
            f = integral(f, x._sage_(), hold=True)
        elif len(limit) == 2:
            x, b = limit
            f = integral(f, x._sage_(), b._sage_(), hold=True)
        else:
            x, a, b = limit
            f = integral(f, (x._sage_(), a._sage_(), b._sage_()), hold=True)
    return f


def _sympysage_derivative(self):
    """
    EXAMPLES::

        sage: from sympy import Derivative
        sage: f = function('f')
        sage: sympy_diff = Derivative(f(x)._sympy_(), x._sympy_())
        sage: assert diff(f(x),x)._sympy_() == sympy_diff
        sage: assert diff(f(x),x) == sympy_diff._sage_()

    TESTS:

    Check that :issue:`28964` is fixed::

        sage: f = function('f')
        sage: _ = var('x,t')
        sage: assert diff(f(x, t), t)._sympy_()._sage_() == diff(f(x, t), t)
        sage: assert diff(f(x, t), x, 2, t)._sympy_()._sage_() == diff(f(x, t), x, 2, t)

        sage: diff(f(x, t), x).integrate(x)
        f(x, t)
        sage: diff(f(x, t), x).integrate(t, algorithm='maxima')
        integrate(diff(f(x, t), x), t)
        sage: diff(f(x, t), x).integrate(t, algorithm='sympy')
        integrate(diff(f(x, t), x), t)
        sage: result = integrate(f(x, t), x).diff(t)
        ...
        sage: result
        integrate(diff(f(x, t), t), x)
    """
    from sage.calculus.functional import derivative
    from sympy.core.containers import Tuple
    f = self.args[0]._sage_()
    args = [a._sage_() for arg in self.args[1:]
            for a in (arg if isinstance(arg, (tuple, Tuple)) else [arg])]
    return derivative(f, *args)


def _sympysage_order(self):
    """
    EXAMPLES::

        sage: from sage.functions.other import Order
        sage: from sympy.series import Order as SOrder
        sage: assert Order(1)._sympy_() == SOrder(1)
        sage: assert Order(1) == SOrder(1)._sage_()
    """
    from sage.functions.other import Order
    return Order(self.args[0])._sage_()


def _sympysage_lambertw(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, LambertW
        sage: assert lambert_w(x)._sympy_() == LambertW(0, Symbol('x'))
        sage: assert lambert_w(x) == LambertW(Symbol('x'))._sage_()
    """
    from sage.functions.log import lambert_w
    return lambert_w(self.args[0]._sage_())


def _sympysage_rf(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, rf
        sage: _ = var('x, y')
        sage: rfxy = rf(Symbol('x'), Symbol('y'))
        sage: assert rising_factorial(x,y)._sympy_() == rfxy.rewrite('gamma', piecewise=False)
        sage: assert rising_factorial(x,y) == rfxy._sage_()
    """
    from sage.arith.misc import rising_factorial
    return rising_factorial(self.args[0]._sage_(), self.args[1]._sage_())


def _sympysage_ff(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, ff
        sage: _ = var('x, y')
        sage: ffxy = ff(Symbol('x'), Symbol('y'))
        sage: assert falling_factorial(x,y)._sympy_() == ffxy.rewrite('gamma') # known bug
        sage: assert falling_factorial(x,y) == ffxy._sage_()
    """
    from sage.arith.misc import falling_factorial
    return falling_factorial(self.args[0]._sage_(), self.args[1]._sage_())


def _sympysage_lgamma(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, loggamma
        sage: assert log_gamma(x)._sympy_() == loggamma(Symbol('x'))
        sage: assert log_gamma(x) == loggamma(Symbol('x'))._sage_()
    """
    from sage.functions.gamma import log_gamma
    return log_gamma(self.args[0]._sage_())


def _sympysage_polygamma(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, polygamma as pg
        sage: _ = var('x, y')
        sage: pgxy = pg(Symbol('x'), Symbol('y'))
        sage: assert psi(x)._sympy_() == pg(0, Symbol('x'))
        sage: assert psi(x) == pg(0, Symbol('x'))._sage_()
        sage: assert psi(x,y)._sympy_() == pgxy
        sage: assert psi(x,y) == pgxy._sage_()
        sage: integrate(psi(x), x, algorithm='sympy')
        integrate(psi(x), x)
    """
    from sage.functions.gamma import psi
    return psi(self.args[0]._sage_(), self.args[1]._sage_())


def _sympysage_dirac_delta(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, DiracDelta
        sage: assert dirac_delta(x)._sympy_() == DiracDelta(Symbol('x'))
        sage: assert dirac_delta(x) == DiracDelta(Symbol('x'))._sage_()
    """
    from sage.functions.generalized import dirac_delta
    return dirac_delta(self.args[0]._sage_())


def _sympysage_heaviside(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, Heaviside
        sage: assert heaviside(x)._sympy_() == Heaviside(Symbol('x'))
        sage: assert heaviside(x) == Heaviside(Symbol('x'))._sage_()
    """
    from sage.functions.generalized import heaviside
    return heaviside(self.args[0]._sage_())


def _sympysage_expint(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, expint
        sage: _ = var('x, y')
        sage: sy = expint(Symbol('x'), Symbol('y'))
        sage: assert exp_integral_e(x,y)._sympy_() == sy
        sage: assert exp_integral_e(x,y) == sy._sage_()
    """
    from sage.functions.exp_integral import exp_integral_e
    return exp_integral_e(self.args[0]._sage_(), self.args[1]._sage_())


def _sympysage_hyp(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, hyper
        sage: _ = var('a,b,p,q,x')
        sage: sy = hyper((Symbol('a'), Symbol('b')), (Symbol('p'), Symbol('q')), Symbol('x'))
        sage: assert hypergeometric((a,b),(p,q),x)._sympy_() == sy
        sage: assert hypergeometric((a,b),(p,q),x) == sy._sage_()
    """
    from sage.functions.hypergeometric import hypergeometric
    ap = [arg._sage_() for arg in self.args[0]]
    bq = [arg._sage_() for arg in self.args[1]]
    return hypergeometric(ap, bq, self.argument._sage_())


def _sympysage_elliptic_k(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, elliptic_k
        sage: assert elliptic_kc(x)._sympy_() == elliptic_k(Symbol('x'))
        sage: assert elliptic_kc(x) == elliptic_k(Symbol('x'))._sage_()
    """
    from sage.functions.special import elliptic_kc
    return elliptic_kc(self.args[0]._sage_())


def _sympysage_kronecker_delta(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, KroneckerDelta
        sage: _ = var('x, y')
        sage: sy = KroneckerDelta(Symbol('x'), Symbol('y'))
        sage: assert kronecker_delta(x,y)._sympy_() == sy
        sage: assert kronecker_delta(x,y) == sy._sage_()
    """
    from sage.functions.generalized import kronecker_delta
    return kronecker_delta(self.args[0]._sage_(), self.args[1]._sage_())


def _sympysage_ceiling(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, ceiling
        sage: assert ceil(x)._sympy_() == ceiling(Symbol('x'))
        sage: assert ceil(x) == ceiling(Symbol('x'))._sage_()
        sage: integrate(ceil(x), x, 0, infinity, algorithm='sympy')
        integrate(ceil(x), x, 0, +Infinity)
    """
    from sage.functions.other import ceil
    return ceil(self.args[0]._sage_())


def _sympysage_piecewise(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, pi as spi, Eq, Lt, Piecewise
        sage: sx = Symbol('x')
        sage: sp = Piecewise((spi, Lt(sx,0)), (1, Eq(sx,1)), (0, True))
        sage: ex = cases(((x<0, pi), (x==1, 1), (True, 0)))
        sage: assert ex._sympy_() == sp
        sage: assert ex == sp._sage_()

        sage: _ = var('y, z')
        sage: (x^y - z).integrate(y, algorithm='sympy')
        -y*z + cases(((log(x) != 0, x^y/log(x)), (1, y)))
    """
    from sage.functions.other import cases
    return cases([(p.cond._sage_(), p.expr._sage_()) for p in self.args])


def _sympysage_fresnels(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, pi as spi, fresnels
        sage: sx = Symbol('x')
        sage: sp = fresnels(sx)
        sage: ex =  fresnel_sin(x)
        sage: assert ex._sympy_() == sp
        sage: assert ex == sp._sage_()
    """
    from sage.functions.error import fresnel_sin
    return fresnel_sin(self.args[0]._sage_())


def _sympysage_fresnelc(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, pi as spi, fresnelc
        sage: sx = Symbol('x')
        sage: sp = fresnelc(sx)
        sage: ex =  fresnel_cos(x)
        sage: assert ex._sympy_() == sp
        sage: assert ex == sp._sage_()
    """
    from sage.functions.error import fresnel_cos
    return fresnel_cos(self.args[0]._sage_())


def _sympysage_besselj(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, besselj
        sage: _ = var('x, y')
        sage: sy = besselj(Symbol('x'), Symbol('y'))
        sage: assert bessel_J(x,y)._sympy_() == sy
        sage: assert bessel_J(x,y) == sy._sage_()
    """
    from sage.functions.bessel import bessel_J
    return bessel_J(self.args[0]._sage_(), self.args[1]._sage_())


def _sympysage_bessely(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, bessely
        sage: _ = var('x, y')
        sage: sy = bessely(Symbol('x'), Symbol('y'))
        sage: assert bessel_Y(x,y)._sympy_() == sy
        sage: assert bessel_Y(x,y) == sy._sage_()
    """
    from sage.functions.bessel import bessel_Y
    return bessel_Y(self.args[0]._sage_(), self.args[1]._sage_())


def _sympysage_besseli(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, besseli
        sage: _ = var('x, y')
        sage: sy = besseli(Symbol('x'), Symbol('y'))
        sage: assert bessel_I(x,y)._sympy_() == sy
        sage: assert bessel_I(x,y) == sy._sage_()
    """
    from sage.functions.bessel import bessel_I
    return bessel_I(self.args[0]._sage_(), self.args[1]._sage_())


def _sympysage_besselk(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, besselk
        sage: _ = var('x, y')
        sage: sy = besselk(Symbol('x'), Symbol('y'))
        sage: assert bessel_K(x,y)._sympy_() == sy
        sage: assert bessel_K(x,y) == sy._sage_()
    """
    from sage.functions.bessel import bessel_K
    return bessel_K(self.args[0]._sage_(), self.args[1]._sage_())


def _sympysage_ynm(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, Ynm
        sage: _ = var('n,m,t,p')
        sage: sy = Ynm(Symbol('n'), Symbol('m'), Symbol('t'), Symbol('p'))
        sage: assert spherical_harmonic(n,m,t,p)._sympy_() == sy
        sage: assert spherical_harmonic(n,m,t,p) == sy._sage_()
    """
    from sage.functions.special import spherical_harmonic
    return spherical_harmonic(self.args[0]._sage_(),
                              self.args[1]._sage_(),
                              self.args[2]._sage_(),
                              self.args[3]._sage_())


def _sympysage_re(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, re
        sage: assert real_part(x)._sympy_() == re(Symbol('x'))
        sage: assert real_part(x) == re(Symbol('x'))._sage_()
    """
    from sage.functions.other import real_part
    return real_part(self.args[0]._sage_())


def _sympysage_im(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, im
        sage: assert imag_part(x)._sympy_() == im(Symbol('x'))
        sage: assert imag_part(x) == im(Symbol('x'))._sage_()
    """
    from sage.functions.other import imag_part
    return imag_part(self.args[0]._sage_())


def _sympysage_abs(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, Abs
        sage: assert abs(x)._sympy_() == Abs(Symbol('x'))
        sage: assert abs(x) == Abs(Symbol('x'))._sage_()
    """
    from sage.functions.other import abs_symbolic
    return abs_symbolic(self.args[0]._sage_())


def _sympysage_crootof(self):
    """
    EXAMPLES::

        sage: from sympy import Symbol, CRootOf
        sage: sobj = CRootOf(Symbol('x')**2 - 2, 1)
        sage: assert complex_root_of(x^2-2, 1)._sympy_() == sobj
        sage: assert complex_root_of(x^2-2, 1) == sobj._sage_()

        sage: from sympy import solve as ssolve
        sage: sols = ssolve(x^6+x+1, x)
        sage: (sols[0]+1)._sage_().n()
        0.209332811185582 - 0.300506920309552*I
    """
    from sage.functions.other import complex_root_of
    from sage.symbolic.ring import SR
    return complex_root_of(self.args[0]._sage_(), SR(self.args[1]))


def _sympysage_matrix(self):
    """
    Convert SymPy matrix ``self`` to Sage.

    EXAMPLES::

        sage: from sympy.matrices import Matrix, SparseMatrix, ImmutableMatrix
        sage: from sage.interfaces.sympy import sympy_init
        sage: from sympy.abc import x
        sage: sympy_init()
        sage: sM = Matrix([[1, x + 1], [x - 1, 1]]); sM
        Matrix([
        [    1, x + 1],
        [x - 1,     1]])
        sage: M = sM._sage_(); M
        [    1 x + 1]
        [x - 1     1]
        sage: M.parent()
        Full MatrixSpace of 2 by 2 dense matrices over Symbolic Ring

        sage: sN = SparseMatrix.eye(3); sN
        Matrix([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1]])
        sage: N = sN._sage_(); N
        [1 0 0]
        [0 1 0]
        [0 0 1]
        sage: N.parent()
        Full MatrixSpace of 3 by 3 sparse matrices over Integer Ring

        sage: sO = SparseMatrix.zeros(3); sO
        Matrix([
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0]])
        sage: O = sO._sage_(); O
        [0 0 0]
        [0 0 0]
        [0 0 0]
        sage: O.parent()
        Full MatrixSpace of 3 by 3 sparse matrices over Integer Ring

    If ``self`` is immutable, the result is cached::

        sage: sImmM = ImmutableMatrix([[1, x + 1], [x - 1, 1]]); sImmM
        Matrix([
        [    1, x + 1],
        [x - 1,     1]])
        sage: ImmM = sImmM._sage_(); ImmM
        [    1 x + 1]
        [x - 1     1]
        sage: ImmM is sImmM._sage_()
        True

    If ``self`` is mutable, the conversion is redone every time::

        sage: sM[0, 0] = 1000
        sage: MutatedM = sM._sage_(); MutatedM
        [ 1000 x + 1]
        [x - 1     1]
        sage: M == MutatedM
        False
    """
    try:
        return self._sage_object
    except AttributeError:
        from sympy.matrices import SparseMatrix, ImmutableMatrix
        from sage.matrix.constructor import matrix

        rows, cols = self.shape
        d = {row_col: value._sage_()
             for row_col, value in self.todok().items()}
        if not d:
            from sage.rings.integer_ring import ZZ
            base_ring = ZZ
        else:
            from sage.structure.element import get_coercion_model
            from sage.symbolic.ring import SR
            coercion_model = get_coercion_model()
            try:
                base_ring = coercion_model.common_parent(*d.values())
            except TypeError:  # no common canonical parent
                base_ring = SR
        result = matrix(base_ring, rows, cols, d,
                        sparse=isinstance(self, SparseMatrix),
                        immutable=True)
        if isinstance(self, ImmutableMatrix):
            self._sage_object = result
        return result


def _sympysage_relational(self):
    """
    EXAMPLES::

        sage: from sympy import Eq, Ne, Gt, Ge, Lt, Le, Symbol
        sage: sx = Symbol('x')
        sage: assert (x == 0)._sympy_() == Eq(sx, 0)
        sage: assert (x == 0) == Eq(x, 0)._sage_()
        sage: assert (x != 0)._sympy_() == Ne(sx, 0)
        sage: assert (x != 0) == Ne(x, 0)._sage_()
        sage: assert (x > 0)._sympy_() == Gt(sx, 0)
        sage: assert (x > 0) == Gt(x, 0)._sage_()
        sage: assert (x >= 0)._sympy_() == Ge(sx, 0)
        sage: assert (x >= 0) == Ge(x, 0)._sage_()
        sage: assert (x < 0)._sympy_() == Lt(sx, 0)
        sage: assert (x < 0) == Lt(x, 0)._sage_()
        sage: assert (x <= 0)._sympy_() == Le(sx, 0)
        sage: assert (x <= 0) == Le(x, 0)._sage_()
    """
    from operator import eq, ne, gt, lt, ge, le
    from sympy import Eq, Ne, Gt, Ge, Lt, Le
    ops = {Eq: eq, Ne: ne, Gt: gt, Lt: lt, Ge: ge, Le: le}
    return ops.get(self.func)(self.lhs._sage_(), self.rhs._sage_())


def _sympysage_false(self):
    """
    EXAMPLES::

        sage: from sympy.logic.boolalg import BooleanFalse
        sage: assert SR(False)._sympy_() == BooleanFalse()      # known bug
        sage: assert SR(False) == BooleanFalse()._sage_()
    """
    from sage.symbolic.ring import SR
    return SR(False)


def _sympysage_true(self):
    """
    EXAMPLES::

        sage: from sympy.logic.boolalg import BooleanTrue
        sage: assert SR(True)._sympy_() == BooleanTrue()      # known bug
        sage: assert SR(True) == BooleanTrue()._sage_()
    """
    from sage.symbolic.ring import SR
    return SR(True)


# ------------------------------------------------------------------


@run_once
def sympy_init():
    """
    Add ``_sage_()`` methods to SymPy objects where needed.

    This gets called with every call to ``Expression._sympy_()``
    so there is only need to call it if you bypass ``_sympy_()`` to
    create SymPy objects. Note that SymPy objects have ``_sage_()``
    methods hard installed but having them inside Sage as
    one file makes them easier to maintain for Sage developers.

    EXAMPLES::

        sage: from sage.interfaces.sympy import sympy_init
        sage: from sympy import Symbol, Abs
        sage: sympy_init()
        sage: assert abs(x) == Abs(Symbol('x'))._sage_()
    """
    from sympy import Add
    if Add._sage_ == _sympysage_add:
        return

    from sympy import Mul, Pow, Symbol, Subs
    from sympy.core.function import (Function, AppliedUndef, Derivative)
    from sympy.core.numbers import (Float, Integer, Rational, Infinity,
                                    NegativeInfinity, ComplexInfinity,
                                    Exp1, Pi, GoldenRatio,
                                    EulerGamma, Catalan, ImaginaryUnit)
    from sympy.core.numbers import NaN as sympy_nan
    from sympy.core.relational import Relational
    from sympy.functions.combinatorial.factorials import (RisingFactorial,
                                                          FallingFactorial)
    from sympy.functions.elementary.complexes import (re, im, Abs)
    from sympy.functions.elementary.exponential import LambertW
    from sympy.functions.elementary.integers import ceiling
    from sympy.functions.elementary.piecewise import Piecewise
    from sympy.functions.special.error_functions import fresnels, fresnelc
    from sympy.functions.special.bessel import (besselj, bessely, besseli, besselk)
    from sympy.functions.special.delta_functions import (DiracDelta, Heaviside)
    from sympy.functions.special.error_functions import expint
    from sympy.functions.special.elliptic_integrals import elliptic_k
    from sympy.functions.special.gamma_functions import loggamma, polygamma
    from sympy.functions.special.hyper import hyper
    from sympy.functions.special.spherical_harmonics import Ynm
    from sympy.functions.special.tensor_functions import KroneckerDelta
    from sympy.logic.boolalg import BooleanTrue, BooleanFalse
    from sympy.integrals.integrals import Integral
    from sympy.polys import Poly
    from sympy.polys.domains.integerring import IntegerRing
    from sympy.polys.domains.rationalfield import RationalField
    from sympy.polys.domains.polynomialring import PolynomialRing
    from sympy.polys.rootoftools import CRootOf
    from sympy.polys.rootisolation import RealInterval, ComplexInterval
    from sympy.series.order import Order
    from sympy.matrices import ImmutableMatrix, ImmutableSparseMatrix, Matrix, SparseMatrix

    Float._sage_ = _sympysage_float
    Integer._sage_ = _sympysage_integer
    Rational._sage_ = _sympysage_rational
    RealInterval._sage_ = _sympysage_real_interval
    ComplexInterval._sage_ = _sympysage_complex_interval
    IntegerRing._sage_ = _sympysage_integer_ring
    RationalField._sage_ = _sympysage_rational_field
    PolynomialRing._sage_ = _sympysage_polynomial_ring
    Poly._sage_ = _sympysage_polynomial
    Infinity._sage_ = _sympysage_pinfty
    NegativeInfinity._sage_ = _sympysage_ninfty
    ComplexInfinity._sage_ = _sympysage_uinfty
    sympy_nan._sage_ = _sympysage_nan
    ImmutableMatrix._sage_ = _sympysage_matrix
    ImmutableSparseMatrix._sage_ = _sympysage_matrix
    Matrix._sage_ = _sympysage_matrix
    SparseMatrix._sage_ = _sympysage_matrix
    Relational._sage_ = _sympysage_relational
    Exp1._sage_ = _sympysage_e
    Pi._sage_ = _sympysage_pi
    GoldenRatio._sage_ = _sympysage_golden_ratio
    EulerGamma._sage_ = _sympysage_eulerg
    Catalan._sage_ = _sympysage_catalan
    ImaginaryUnit._sage_ = _sympysage_i
    Add._sage_ = _sympysage_add
    Mul._sage_ = _sympysage_mul
    Pow._sage_ = _sympysage_pow
    Symbol._sage_ = _sympysage_symbol
    Subs._sage_ = _sympysage_Subs
    Function._sage_ = _sympysage_function
    AppliedUndef._sage_ = _sympysage_function
    import sympy.core.function
    sympy.core.function._undef_sage_helper = UndefSageHelper()
    Integral._sage_ = _sympysage_integral
    Derivative._sage_ = _sympysage_derivative
    Order._sage_ = _sympysage_order
    LambertW._sage_ = _sympysage_lambertw
    RisingFactorial._sage_ = _sympysage_rf
    FallingFactorial._sage_ = _sympysage_ff
    loggamma._sage_ = _sympysage_lgamma
    polygamma._sage_ = _sympysage_polygamma
    DiracDelta._sage_ = _sympysage_dirac_delta
    Heaviside._sage_ = _sympysage_heaviside
    expint._sage_ = _sympysage_expint
    hyper._sage_ = _sympysage_hyp
    elliptic_k._sage_ = _sympysage_elliptic_k
    KroneckerDelta._sage_ = _sympysage_kronecker_delta
    Piecewise._sage_ = _sympysage_piecewise
    fresnels._sage_ = _sympysage_fresnels
    fresnelc._sage_ = _sympysage_fresnelc
    besselj._sage_ = _sympysage_besselj
    bessely._sage_ = _sympysage_bessely
    besseli._sage_ = _sympysage_besseli
    besselk._sage_ = _sympysage_besselk
    Ynm._sage_ = _sympysage_ynm
    re._sage_ = _sympysage_re
    im._sage_ = _sympysage_im
    Abs._sage_ = _sympysage_abs
    CRootOf._sage_ = _sympysage_crootof
    BooleanFalse._sage_ = _sympysage_false
    BooleanTrue._sage_ = _sympysage_true
    ceiling._sage_ = _sympysage_ceiling


def check_expression(expr, var_symbols, only_from_sympy=False):
    """
    Do ``eval(expr)`` both in Sage and SymPy and other checks.

    EXAMPLES::

        sage: from sage.interfaces.sympy import check_expression
        sage: check_expression("1.123*x", "x")
    """
    from sage.symbolic.ring import SR
    from sympy import (__dict__ as sympydict, Basic, S, var as svar)
    # evaluate the expression in the context of Sage:
    if var_symbols:
        SR.var(var_symbols)
    is_different = False
    try:
        e_sage = SR(expr)
        assert not isinstance(e_sage, Basic)
    except (NameError, TypeError):
        is_different = True

    # evaluate the expression in the context of SymPy:
    if var_symbols:
        svar(var_symbols)
    b = globals().copy()
    b.update(sympydict)
    assert "sin" in b
    b.update(sympydict)
    e_sympy = eval(expr, b)
    assert isinstance(e_sympy, Basic)

    # Sympy func may have specific _sage_ method
    if is_different:
        _sage_method = getattr(e_sympy.func, "_sage_")
        e_sage = _sage_method(S(e_sympy))

    # Do the actual checks:
    if not only_from_sympy:
        assert S(e_sage) == e_sympy
    assert e_sage == SR(e_sympy)


def test_all():
    """
    Call some tests that were originally in SymPy.

    EXAMPLES::

        sage: from sage.interfaces.sympy import test_all
        sage: test_all()
    """
    def test_basics():
        check_expression("x", "x")
        check_expression("x**2", "x")
        check_expression("x**2+y**3", "x y")
        check_expression("1/(x+y)**2-x**3/4", "x y")

    def test_complex():
        check_expression("I", "")
        check_expression("23+I*4", "x")

    def test_complex_fail():
        # Sage doesn't properly implement _sympy_ on I
        check_expression("I*y", "y")
        check_expression("x+I*y", "x y")

    def test_integer():
        check_expression("4*x", "x")
        check_expression("-4*x", "x")

    def test_real():
        check_expression("1.123*x", "x")
        check_expression("-18.22*x", "x")

    def test_functions():
        # Test at least one Function without own _sage_ method
        from sympy import factorial
        assert "_sage_" not in factorial.__dict__
        check_expression("factorial(x)", "x")
        check_expression("sin(x)", "x")
        check_expression("cos(x)", "x")
        check_expression("tan(x)", "x")
        check_expression("cot(x)", "x")
        check_expression("asin(x)", "x")
        check_expression("acos(x)", "x")
        check_expression("atan(x)", "x")
        check_expression("atan2(y, x)", "x, y")
        check_expression("acot(x)", "x")
        check_expression("sinh(x)", "x")
        check_expression("cosh(x)", "x")
        check_expression("tanh(x)", "x")
        check_expression("coth(x)", "x")
        check_expression("asinh(x)", "x")
        check_expression("acosh(x)", "x")
        check_expression("atanh(x)", "x")
        check_expression("acoth(x)", "x")
        check_expression("exp(x)", "x")
        check_expression("log(x)", "x")
        check_expression("abs(x)", "x")
        check_expression("arg(x)", "x")
        check_expression("conjugate(x)", "x")

    def test_issue_4023():
        from sage.symbolic.ring import SR
        from sage.functions.all import log
        from sympy import integrate, simplify
        a, x = SR.var("a x")
        i = integrate(log(x) / a, (x, a, a + 1))
        i2 = simplify(i)
        s = SR(i2)
        assert s == (a * log(1 + a) - a * log(a) + log(1 + a) - 1) / a

    def test_integral():
        # test Sympy-->Sage
        check_expression("Integral(x, (x,))", "x", only_from_sympy=True)
        check_expression("Integral(x, (x, 0, 1))", "x", only_from_sympy=True)
        check_expression("Integral(x*y, (x,), (y, ))", "x,y", only_from_sympy=True)
        check_expression("Integral(x*y, (x,), (y, 0, 1))", "x,y", only_from_sympy=True)
        check_expression("Integral(x*y, (x, 0, 1), (y,))", "x,y", only_from_sympy=True)
        check_expression("Integral(x*y, (x, 0, 1), (y, 0, 1))", "x,y", only_from_sympy=True)
        check_expression("Integral(x*y*z, (x, 0, 1), (y, 0, 1), (z, 0, 1))", "x,y,z", only_from_sympy=True)

    def test_integral_failing():
        # Note: sage may attempt to turn this into Integral(x, (x, x, 0))
        check_expression("Integral(x, (x, 0))", "x", only_from_sympy=True)
        check_expression("Integral(x*y, (x,), (y, 0))", "x,y", only_from_sympy=True)
        check_expression("Integral(x*y, (x, 0, 1), (y, 0))", "x,y", only_from_sympy=True)

    def test_undefined_function():
        from sage.symbolic.ring import SR
        from sage.calculus.var import function
        from sympy import Symbol, Function
        f = function('f')
        sf = Function('f')
        x, y = SR.var('x y')
        sx = Symbol('x')
        sy = Symbol('y')
        assert f(x)._sympy_() == sf(sx)
        assert f(x) == sf(sx)._sage_()
        assert f(x, y)._sympy_() == sf(sx, sy)
        assert f(x, y) == sf(sx, sy)._sage_()
        assert f._sympy_() == sf
        assert f == sf._sage_()

    test_basics()
    test_complex()
    test_complex_fail()
    test_integer()
    test_real()
    test_functions()
    test_issue_4023()
    test_integral()
    # test_integral_failing()
    test_undefined_function()


def sympy_set_to_list(set, vars):
    """
    Convert all set objects that can be returned by SymPy's solvers.
    """
    from sage.rings.infinity import UnsignedInfinity
    from sympy import (FiniteSet, And, Or, Union, Interval, oo, S)
    from sympy.core.relational import Relational
    if set == S.Reals:
        return [x._sage_() < oo for x in vars]
    elif set == S.Complexes:
        return [x._sage_() != UnsignedInfinity for x in vars]
    elif set is None or set == S.EmptySet:
        return []
    if isinstance(set, (And, Or, Relational)):
        if isinstance(set, And):
            return [[item for rel in set._args[0]
                    for item in sympy_set_to_list(rel, vars)]]
        elif isinstance(set, Or):
            return [sympy_set_to_list(iv, vars) for iv in set._args[0]]
        elif isinstance(set, Relational):
            return [set._sage_()]
    elif isinstance(set, FiniteSet):
        x = vars[0]
        return [x._sage_() == arg._sage_() for arg in set.args]
    elif isinstance(set, (Union, Interval)):
        x = vars[0]
        if isinstance(set, Interval):
            left, right, lclosed, rclosed = set._args
            if lclosed:
                rel1 = [x._sage_() > left._sage_()]
            else:
                rel1 = [x._sage_() >= left._sage_()]
            if rclosed:
                rel2 = [x._sage_() < right._sage_()]
            else:
                return [x._sage_() <= right._sage_()]
            if right == oo:
                return rel1
            if left == -oo:
                return rel2
            return [rel1, rel2]
        if isinstance(set, Union):
            return [sympy_set_to_list(iv, vars) for iv in set._args]
    return set
