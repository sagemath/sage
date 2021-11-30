"""
Other functions

TESTS:

Check that gamma function imports are deprecated (:trac:`24411`)::

    sage: from sage.functions.other import beta
    sage: beta(x, x)
    doctest:...: DeprecationWarning:
    Importing beta from here is deprecated. If you need to use it, please import it directly from sage.functions.gamma
    See http://trac.sagemath.org/24411 for details.
    beta(x, x)
"""

from sage.misc.lazy_import import lazy_import
lazy_import('sage.functions.gamma',
            ('gamma', 'log_gamma', 'gamma_inc',
             'gamma_inc_lower', 'psi', 'beta'), deprecation=24411)

from sage.symbolic.function import GinacFunction, BuiltinFunction
from sage.symbolic.expression import Expression, register_symbol, symbol_table
from sage.symbolic.ring import SR, SymbolicRing
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.rational import Rational
from sage.rings.complex_mpfr import ComplexField
from sage.rings.real_mpfr import RealField
from sage.misc.latex import latex
from sage.structure.element import Element
import math

from sage.structure.element import coercion_model

# avoid name conflicts with `parent` as a function parameter
from sage.structure.all import parent as s_parent

from sage.functions.trig import arctan2

from sage.arith.all import binomial as arith_binomial

from sage.misc.functional import sqrt


class Function_abs(GinacFunction):
    def __init__(self):
        r"""
        The absolute value function.

        EXAMPLES::

            sage: var('x y')
            (x, y)
            sage: abs(x)
            abs(x)
            sage: abs(x^2 + y^2)
            abs(x^2 + y^2)
            sage: abs(-2)
            2
            sage: sqrt(x^2)
            sqrt(x^2)
            sage: abs(sqrt(x))
            sqrt(abs(x))
            sage: complex(abs(3*I))
            (3+0j)

            sage: f = sage.functions.other.Function_abs()
            sage: latex(f)
            \mathrm{abs}
            sage: latex(abs(x))
            {\left| x \right|}
            sage: abs(x)._sympy_()
            Abs(x)

        Test pickling::

            sage: loads(dumps(abs(x)))
            abs(x)

        TESTS:

        Check that :trac:`12588` is fixed::

            sage: abs(pi*I)
            pi
            sage: abs(pi*I*catalan)
            catalan*pi
            sage: abs(pi*catalan*x)
            catalan*pi*abs(x)
            sage: abs(pi*I*catalan*x)
            catalan*pi*abs(x)
            sage: abs(1.0j*pi)
            1.00000000000000*pi
            sage: abs(I*x)
            abs(x)
            sage: abs(I*pi)
            pi
            sage: abs(I*log(2))
            log(2)
            sage: abs(I*e^5)
            e^5
            sage: abs(log(1/2))
            -log(1/2)
            sage: abs(log(3/2))
            log(3/2)
            sage: abs(log(1/2)*log(1/3))
            log(1/2)*log(1/3)
            sage: abs(log(1/2)*log(1/3)*log(1/4))
            -log(1/2)*log(1/3)*log(1/4)
            sage: abs(log(1/2)*log(1/3)*log(1/4)*i)
            -log(1/2)*log(1/3)*log(1/4)
            sage: abs(log(x))
            abs(log(x))
            sage: abs(zeta(I))
            abs(zeta(I))
            sage: abs(e^2*x)
            abs(x)*e^2
            sage: abs((pi+e)*x)
            (pi + e)*abs(x)

            sage: fricas(abs(x)).sage().derivative()  # optional - fricas
            1/2*(x + conjugate(x))/abs(x)
        """
        GinacFunction.__init__(self, "abs", latex_name=r"\mathrm{abs}",
                               conversions=dict(sympy='Abs',
                                                mathematica='Abs',
                                                giac='abs',
                                                fricas='abs'))

abs = abs_symbolic = Function_abs()


def _eval_floor_ceil(self, x, method, bits=0, **kwds):
    """
    Helper function to compute ``floor(x)`` or ``ceil(x)``.

    INPUT:

    - ``x`` -- a number

    - ``method`` -- should be either ``"floor"`` or ``"ceil"``

    - ``bits`` -- how many bits to use before giving up

    See :class:`Function_floor` and :class:`Function_ceil` for examples
    and tests.

    TESTS::

        sage: numbers = [SR(10^100 + exp(-100)), SR(10^100 - exp(-100)), SR(10^100)]
        sage: numbers += [-n for n in numbers]
        sage: for n in numbers:
        ....:     f = floor(n)
        ....:     c = ceil(n)
        ....:     if f == c:
        ....:         assert n in ZZ
        ....:     else:
        ....:         assert f + 1 == c

    A test from :trac:`12121`::

        sage: e1 = pi - continued_fraction(pi).convergent(2785)
        sage: e2 = e - continued_fraction(e).convergent(1500)
        sage: f = e1/e2
        sage: f = 1 / (f - continued_fraction(f).convergent(1000))
        sage: f = f - continued_fraction(f).convergent(1)
        sage: floor(f, bits=10000)
        -1
        sage: ceil(f, bits=10000)
        0

    These do not work but fail gracefully::

        sage: ceil(Infinity)
        Traceback (most recent call last):
        ...
        ValueError: Calling ceil() on infinity or NaN
        sage: ceil(NaN)
        Traceback (most recent call last):
        ...
        ValueError: Calling ceil() on infinity or NaN

    Test that elements of symbolic subrings work in the same way as
    elements of ``SR``, :trac:`32724`::

        sage: SCR = SR.subring(no_variables=True)
        sage: floor(log(2^(3/2)) / log(2) + 1/2)
        2
        sage: floor(SCR(log(2^(-3/2)) / log(2) + 1/2))
        -1
    """
    # First, some obvious things...
    try:
        m = getattr(x, method)
    except AttributeError:
        pass
    else:
        return m()

    if isinstance(x, int):
        return Integer(x)
    if isinstance(x, (float, complex)):
        m = getattr(math, method)
        return Integer(m(x))
    if type(x).__module__ == 'numpy':
        import numpy
        m = getattr(numpy, method)
        return m(x)

    # The strategy is to convert the number to an interval field and
    # hope that this interval will have a unique floor/ceiling.
    #
    # There are 2 reasons why this could fail:
    # (A) The expression is very complicated and we simply require
    #     more bits.
    # (B) The expression is a non-obvious exact integer. In this
    #     case, adding bits will not help since an interval around
    #     an integer will not have a unique floor/ceiling, no matter
    #     how many bits are used.
    #
    # The strategy is to first reduce the absolute diameter of the
    # interval until its size is at most 10^(-6). Then we check for
    # (B) by simplifying the expression.
    from sage.rings.all import RealIntervalField

    # Might it be needed to simplify x? This only applies for
    # elements of SR (or its subrings)
    need_to_simplify = isinstance(s_parent(x), SymbolicRing)

    # An integer which is close to x. We use this to increase precision
    # by subtracting this guess before converting to an interval field.
    # This mostly helps with the case that x is close to, but not equal
    # to, an exact integer.
    guess = Integer(0)

    # We do not use the target number of bits immediately, we just use
    # it as indication of when to stop.
    target_bits = bits
    bits = 32
    attempts = 5
    while attempts:
        attempts -= 1
        if not attempts and bits < target_bits:
            # Add one more attempt as long as the precision is less
            # than requested
            attempts = 1

        RIF = RealIntervalField(bits)
        if guess:
            y = x - guess
        else:
            y = x
        try:
            y_interval = RIF(y)
        except TypeError:
            # If we cannot compute a numerical enclosure, leave the
            # expression unevaluated.
            return BuiltinFunction.__call__(self, SR(x))
        diam = y_interval.absolute_diameter()
        if diam.is_infinity():
            # We have a very bad approximation => increase the number
            # of bits a lot
            bits *= 4
            continue
        fdiam = float(diam)
        if fdiam >= 1.0:
            # Increase number of bits to get to a diameter less than
            # 2^(-32), assuming that the diameter scales as 2^(-bits)
            bits += 32 + int(diam.log2())
            continue

        # Compute ceil/floor of both ends of the interval:
        # if these match, we are done!
        a = getattr(y_interval.lower(), method)()
        b = getattr(y_interval.upper(), method)()
        if a == b:
            return a + guess

        # Compute a better guess for the next attempt. Since diam < 1,
        # there is a unique integer in our interval. This integer equals
        # the ceil of the lower bound and the floor of the upper bound.
        if self is floor:
            guess += b
        else:
            assert self is ceil
            guess += a

        if need_to_simplify and fdiam <= 1e-6:
            x = x.full_simplify().canonicalize_radical()
            need_to_simplify = False
            continue

        bits *= 2

    raise ValueError("cannot compute {}({!r}) using {} bits of precision".format(method, x, RIF.precision()))


class Function_ceil(BuiltinFunction):
    def __init__(self):
        r"""
        The ceiling function.

        The ceiling of `x` is computed in the following manner.


        #. The ``x.ceil()`` method is called and returned if it
           is there. If it is not, then Sage checks if `x` is one of
           Python's native numeric data types. If so, then it calls and
           returns ``Integer(math.ceil(x))``.

        #. Sage tries to convert `x` into a
           ``RealIntervalField`` with 53 bits of precision. Next,
           the ceilings of the endpoints are computed. If they are the same,
           then that value is returned. Otherwise, the precision of the
           ``RealIntervalField`` is increased until they do match
           up or it reaches ``bits`` of precision.

        #. If none of the above work, Sage returns a
           ``Expression`` object.


        EXAMPLES::

            sage: a = ceil(2/5 + x)
            sage: a
            ceil(x + 2/5)
            sage: a(x=4)
            5
            sage: a(x=4.0)
            5
            sage: ZZ(a(x=3))
            4
            sage: a = ceil(x^3 + x + 5/2); a
            ceil(x^3 + x + 5/2)
            sage: a.simplify()
            ceil(x^3 + x + 1/2) + 2
            sage: a(x=2)
            13

        ::

            sage: ceil(sin(8)/sin(2))
            2

        ::

            sage: ceil(5.4)
            6
            sage: type(ceil(5.4))
            <class 'sage.rings.integer.Integer'>

        ::

            sage: ceil(factorial(50)/exp(1))
            11188719610782480504630258070757734324011354208865721592720336801
            sage: ceil(SR(10^50 + 10^(-50)))
            100000000000000000000000000000000000000000000000001
            sage: ceil(SR(10^50 - 10^(-50)))
            100000000000000000000000000000000000000000000000000

        Small numbers which are extremely close to an integer are hard to
        deal with::

            sage: ceil((33^100 + 1)^(1/100))
            Traceback (most recent call last):
            ...
            ValueError: cannot compute ceil(...) using 256 bits of precision

        This can be fixed by giving a sufficiently large ``bits`` argument::

            sage: ceil((33^100 + 1)^(1/100), bits=500)
            Traceback (most recent call last):
            ...
            ValueError: cannot compute ceil(...) using 512 bits of precision
            sage: ceil((33^100 + 1)^(1/100), bits=1000)
            34

        ::

            sage: ceil(sec(e))
            -1

            sage: latex(ceil(x))
            \left \lceil x \right \rceil
            sage: ceil(x)._sympy_()
            ceiling(x)

        ::

            sage: import numpy
            sage: a = numpy.linspace(0,2,6)
            sage: ceil(a)
            array([0., 1., 1., 2., 2., 2.])

        Test pickling::

            sage: loads(dumps(ceil))
            ceil
        """
        BuiltinFunction.__init__(self, "ceil",
                                   conversions=dict(maxima='ceiling',
                                                    sympy='ceiling',
                                                    giac='ceil'))

    def _print_latex_(self, x):
        r"""
        EXAMPLES::

            sage: latex(ceil(x)) # indirect doctest
            \left \lceil x \right \rceil
        """
        return r"\left \lceil %s \right \rceil"%latex(x)

    #FIXME: this should be moved to _eval_
    def __call__(self, x, **kwds):
        """
        Allows an object of this class to behave like a function. If
        ``ceil`` is an instance of this class, we can do ``ceil(n)`` to get
        the ceiling of ``n``.

        TESTS::

            sage: ceil(SR(10^50 + 10^(-50)))
            100000000000000000000000000000000000000000000000001
            sage: ceil(SR(10^50 - 10^(-50)))
            100000000000000000000000000000000000000000000000000
            sage: ceil(int(10^50))
            100000000000000000000000000000000000000000000000000
            sage: ceil((1725033*pi - 5419351)/(25510582*pi - 80143857))
            -2
        """
        return _eval_floor_ceil(self, x, "ceil", **kwds)

    def _eval_(self, x):
        """
        EXAMPLES::

            sage: ceil(x).subs(x==7.5)
            8
            sage: ceil(x)
            ceil(x)

            sage: var('x',domain='integer')
            x
            sage: ceil(x)
            x
            sage: ceil(factorial(x) + binomial(x^2, x))
            binomial(x^2, x) + factorial(x)
            sage: ceil(gamma(abs(2*x)+1) * real(x))
            x*gamma(2*abs(x) + 1)
            sage: forget()
        """
        try:
            if SR(x).variables() and x.is_integer():
                return x
        except TypeError:
            pass
        try:
            return x.ceil()
        except AttributeError:
            if isinstance(x, int):
                return Integer(x)
            elif isinstance(x, (float, complex)):
                return Integer(math.ceil(x))
        return None

ceil = Function_ceil()


class Function_floor(BuiltinFunction):
    def __init__(self):
        r"""
        The floor function.

        The floor of `x` is computed in the following manner.


        #. The ``x.floor()`` method is called and returned if
           it is there. If it is not, then Sage checks if `x` is one
           of Python's native numeric data types. If so, then it calls and
           returns ``Integer(math.floor(x))``.

        #. Sage tries to convert `x` into a
           ``RealIntervalField`` with 53 bits of precision. Next,
           the floors of the endpoints are computed. If they are the same,
           then that value is returned. Otherwise, the precision of the
           ``RealIntervalField`` is increased until they do match
           up or it reaches ``bits`` of precision.

        #. If none of the above work, Sage returns a
           symbolic ``Expression`` object.


        EXAMPLES::

            sage: floor(5.4)
            5
            sage: type(floor(5.4))
            <class 'sage.rings.integer.Integer'>
            sage: var('x')
            x
            sage: a = floor(5.4 + x); a
            floor(x + 5.40000000000000)
            sage: a.simplify()
            floor(x + 0.4000000000000004) + 5
            sage: a(x=2)
            7

        ::

            sage: floor(cos(8) / cos(2))
            0
            sage: floor(log(4) / log(2))
            2
            sage: a = floor(5.4 + x); a
            floor(x + 5.40000000000000)
            sage: a.subs(x==2)
            7
            sage: floor(log(2^(3/2)) / log(2) + 1/2)
            2
            sage: floor(log(2^(-3/2)) / log(2) + 1/2)
            -1

        ::

            sage: floor(factorial(50)/exp(1))
            11188719610782480504630258070757734324011354208865721592720336800
            sage: floor(SR(10^50 + 10^(-50)))
            100000000000000000000000000000000000000000000000000
            sage: floor(SR(10^50 - 10^(-50)))
            99999999999999999999999999999999999999999999999999
            sage: floor(int(10^50))
            100000000000000000000000000000000000000000000000000

        Small numbers which are extremely close to an integer are hard to
        deal with::

            sage: floor((33^100 + 1)^(1/100))
            Traceback (most recent call last):
            ...
            ValueError: cannot compute floor(...) using 256 bits of precision

        This can be fixed by giving a sufficiently large ``bits`` argument::

            sage: floor((33^100 + 1)^(1/100), bits=500)
            Traceback (most recent call last):
            ...
            ValueError: cannot compute floor(...) using 512 bits of precision
            sage: floor((33^100 + 1)^(1/100), bits=1000)
            33

        ::

            sage: import numpy
            sage: a = numpy.linspace(0,2,6)
            sage: floor(a)
            array([0., 0., 0., 1., 1., 2.])
            sage: floor(x)._sympy_()
            floor(x)

        Test pickling::

            sage: loads(dumps(floor))
            floor
        """
        BuiltinFunction.__init__(self, "floor",
                                 conversions=dict(sympy='floor', giac='floor'))

    def _print_latex_(self, x):
        r"""
        EXAMPLES::

            sage: latex(floor(x))
            \left \lfloor x \right \rfloor
        """
        return r"\left \lfloor %s \right \rfloor"%latex(x)

    #FIXME: this should be moved to _eval_
    def __call__(self, x, **kwds):
        """
        Allows an object of this class to behave like a function. If
        ``floor`` is an instance of this class, we can do ``floor(n)`` to
        obtain the floor of ``n``.

        TESTS::

            sage: floor(SR(10^50 + 10^(-50)))
            100000000000000000000000000000000000000000000000000
            sage: floor(SR(10^50 - 10^(-50)))
            99999999999999999999999999999999999999999999999999
            sage: floor(int(10^50))
            100000000000000000000000000000000000000000000000000
            sage: floor((1725033*pi - 5419351)/(25510582*pi - 80143857))
            -3
        """
        return _eval_floor_ceil(self, x, "floor", **kwds)

    def _eval_(self, x):
        """
        EXAMPLES::

            sage: floor(x).subs(x==7.5)
            7
            sage: floor(x)
            floor(x)

            sage: var('x',domain='integer')
            x
            sage: floor(x)
            x
            sage: floor(factorial(x) + binomial(x^2, x))
            binomial(x^2, x) + factorial(x)
            sage: floor(gamma(abs(2*x)+1) * real(x))
            x*gamma(2*abs(x) + 1)
            sage: forget()
        """
        try:
            if SR(x).variables() and x.is_integer():
                return x
        except TypeError:
            pass
        try:
            return x.floor()
        except AttributeError:
            if isinstance(x, int):
                return Integer(x)
            elif isinstance(x, (float, complex)):
                return Integer(math.floor(x))
        return None

floor = Function_floor()


class Function_Order(GinacFunction):
    def __init__(self):
        r"""
        The order function.

        This function gives the order of magnitude of some expression,
        similar to `O`-terms.

        .. SEEALSO::

            :meth:`~sage.symbolic.expression.Expression.Order`,
            :mod:`~sage.rings.big_oh`

        EXAMPLES::

            sage: x = SR('x')
            sage: x.Order()
            Order(x)
            sage: (x^2 + x).Order()
            Order(x^2 + x)

        TESTS:

        Check that :trac:`19425` is resolved::

            sage: x.Order().operator()
            Order
        """
        GinacFunction.__init__(self, "Order",
                conversions=dict(),
                latex_name=r"\mathcal{O}")

    def _sympy_(self, arg):
        """
        EXAMPLES::

            sage: x.Order()._sympy_()
            O(x)
            sage: SR(1).Order()._sympy_()
            O(1)
            sage: ((x-1)^3).Order()._sympy_()
            O((x - 1)**3, (x, 1))
            sage: exp(x).series(x==1, 3)._sympy_()
            E + E*(x - 1) + E*(x - 1)**2/2 + O((x - 1)**3, (x, 1))

            sage: (-(pi-x)^3).Order()._sympy_()
            O((x - pi)**3, (x, pi))
            sage: cos(x).series(x==pi, 3)._sympy_()
            -1 + (pi - x)**2/2 + O((x - pi)**3, (x, pi))
        """
        roots = arg.solve(arg.default_variable(), algorithm='sympy',
                          multiplicities=False, explicit_solutions=True)
        if len(roots) == 1:
            arg = (arg, (roots[0].lhs(), roots[0].rhs()))
        elif len(roots) > 1:
            raise ValueError("order term %s has multiple roots" % arg)
        # else there are no roots, e.g. O(1), so we leave arg unchanged
        import sympy
        return sympy.O(*sympy.sympify(arg, evaluate=False))

Order = Function_Order()


class Function_frac(BuiltinFunction):
    def __init__(self):
        r"""
        The fractional part function `\{x\}`.

        ``frac(x)`` is defined as `\{x\} = x - \lfloor x\rfloor`.

        EXAMPLES::

            sage: frac(5.4)
            0.400000000000000
            sage: type(frac(5.4))
            <class 'sage.rings.real_mpfr.RealNumber'>
            sage: frac(456/123)
            29/41
            sage: var('x')
            x
            sage: a = frac(5.4 + x); a
            frac(x + 5.40000000000000)
            sage: frac(cos(8)/cos(2))
            cos(8)/cos(2)
            sage: latex(frac(x))
            \operatorname{frac}\left(x\right)
            sage: frac(x)._sympy_()
            frac(x)

        Test pickling::

            sage: loads(dumps(floor))
            floor
        """
        BuiltinFunction.__init__(self, "frac",
                                 conversions=dict(sympy='frac'),
                                 latex_name=r"\operatorname{frac}")

    def _evalf_(self, x, **kwds):
        """
        EXAMPLES::

            sage: frac(pi).n()
            0.141592653589793
            sage: frac(pi).n(200)
            0.14159265358979323846264338327950288419716939937510582097494
        """
        return x - floor(x)

    def _eval_(self, x):
        """
        EXAMPLES::

            sage: frac(x).subs(x==7.5)
            0.500000000000000
            sage: frac(x)
            frac(x)
        """
        try:
            return x - x.floor()
        except AttributeError:
            if isinstance(x, int):
                return Integer(0)
            elif isinstance(x, (float, complex)):
                return x - Integer(math.floor(x))
            elif isinstance(x, Expression):
                ret = floor(x)
                if not hasattr(ret, "operator") or not ret.operator() == floor:
                    return x - ret
        return None

frac = Function_frac()


# register sqrt in pynac symbol_table for conversion back from other systems
register_symbol(sqrt, dict(mathematica='Sqrt'))
symbol_table['functions']['sqrt'] = sqrt

Function_sqrt = type('deprecated_sqrt', (),
        {'__call__': staticmethod(sqrt),
            '__setstate__': lambda x, y: None})


class Function_real_nth_root(BuiltinFunction):
    r"""
    Real `n`-th root function `x^\frac{1}{n}`.

    The function assumes positive integer `n` and real number `x`.

    EXAMPLES::

        sage: real_nth_root(2, 3)
        2^(1/3)
        sage: real_nth_root(-2, 3)
        -2^(1/3)
        sage: real_nth_root(8, 3)
        2
        sage: real_nth_root(-8, 3)
        -2

        sage: real_nth_root(-2, 4)
        Traceback (most recent call last):
        ...
        ValueError: no real nth root of negative real number with even n

    For numeric input, it gives a numerical approximation. ::

        sage: real_nth_root(2., 3)
        1.25992104989487
        sage: real_nth_root(-2., 3)
        -1.25992104989487

    Some symbolic calculus::

        sage: f = real_nth_root(x, 5)^3
        sage: f
        real_nth_root(x^3, 5)
        sage: f.diff()
        3/5*x^2*real_nth_root(x^(-12), 5)
        sage: result = f.integrate(x)
        ...
        sage: result
        integrate((abs(x)^3)^(1/5)*sgn(x^3), x)
        sage: _.diff()
        (abs(x)^3)^(1/5)*sgn(x^3)
    """
    def __init__(self):
        r"""
        Initialize.

        TESTS::

            sage: cube_root = real_nth_root(x, 3)
            sage: loads(dumps(cube_root))
            real_nth_root(x, 3)

        ::

            sage: f = real_nth_root(x, 3)
            sage: f._sympy_()
            Piecewise((Abs(x)**(1/3)*sign(x), Eq(im(x), 0)), (x**(1/3), True))

        """
        BuiltinFunction.__init__(self, "real_nth_root", nargs=2,
                                 conversions=dict(sympy='real_root',
                                                  mathematica='Surd',
                                                  maple='surd'))

    def _print_latex_(self, base, exp):
        r"""
        TESTS::

            sage: latex(real_nth_root(x, 3))
            x^{\frac{1}{3}}
            sage: latex(real_nth_root(x^2 + x, 3))
            {\left(x^{2} + x\right)}^{\frac{1}{3}}
        """
        return latex(base**(1/exp))

    def _evalf_(self, base, exp, parent=None):
        """
        TESTS::

            sage: real_nth_root(RDF(-2), 3)
            -1.25992104989487...
            sage: real_nth_root(Reals(100)(2), 2)
            1.4142135623730950488016887242
        """
        negative = base < 0

        if negative:
            if exp % 2 == 0:
                raise ValueError('no real nth root of negative real number with even n')
            base = -base

        r = base**(1/exp)

        if negative:
            return -r
        else:
            return r

    def _eval_(self, base, exp):
        """
        TESTS::

            sage: real_nth_root(x, 1)
            x
            sage: real_nth_root(x, 3)
            real_nth_root(x, 3)

            sage: real_nth_root(RIF(2), 3)
            1.259921049894873?
            sage: real_nth_root(RBF(2), 3)
            [1.259921049894873 +/- 3.92e-16]
        """
        if not isinstance(base, Expression) and not isinstance(exp, Expression):
            if isinstance(base, Integer):
                try:
                    return base.nth_root(exp)
                except ValueError:
                    pass
            return self._evalf_(base, exp, parent=s_parent(base))

        if isinstance(exp, Integer) and exp.is_one():
            return base

    def _power_(self, base, exp, power_param=None):
        """
        TESTS::

            sage: f = real_nth_root(x, 3)
            sage: f^5
            real_nth_root(x^5, 3)
        """
        return self(base**power_param, exp)

    def _derivative_(self, base, exp, diff_param=None):
        """
        TESTS::

            sage: f = real_nth_root(x, 3)
            sage: f.diff()
            1/3*real_nth_root(x^(-2), 3)
            sage: f = real_nth_root(-x, 3)
            sage: f.diff()
            -1/3*real_nth_root(x^(-2), 3)
            sage: f = real_nth_root(x, 4)
            sage: f.diff()
            1/4*real_nth_root(x^(-3), 4)
            sage: f = real_nth_root(-x, 4)
            sage: f.diff()
            -1/4*real_nth_root(-1/x^3, 4)
        """
        return 1/exp * self(base, exp)**(1-exp)

real_nth_root = Function_real_nth_root()


class Function_arg(BuiltinFunction):
    def __init__(self):
        r"""
        The argument function for complex numbers.

        EXAMPLES::

            sage: arg(3+i)
            arctan(1/3)
            sage: arg(-1+i)
            3/4*pi
            sage: arg(2+2*i)
            1/4*pi
            sage: arg(2+x)
            arg(x + 2)
            sage: arg(2.0+i+x)
            arg(x + 2.00000000000000 + 1.00000000000000*I)
            sage: arg(-3)
            pi
            sage: arg(3)
            0
            sage: arg(0)
            0

            sage: latex(arg(x))
            {\rm arg}\left(x\right)
            sage: maxima(arg(x))
            atan2(0,_SAGE_VAR_x)
            sage: maxima(arg(2+i))
            atan(1/2)
            sage: maxima(arg(sqrt(2)+i))
            atan(1/sqrt(2))
            sage: arg(x)._sympy_()
            arg(x)

            sage: arg(2+i)
            arctan(1/2)
            sage: arg(sqrt(2)+i)
            arg(sqrt(2) + I)
            sage: arg(sqrt(2)+i).simplify()
            arctan(1/2*sqrt(2))

        TESTS::

            sage: arg(0.0)
            0.000000000000000
            sage: arg(3.0)
            0.000000000000000
            sage: arg(-2.5)
            3.14159265358979
            sage: arg(2.0+3*i)
            0.982793723247329
        """
        BuiltinFunction.__init__(self, "arg",
                conversions=dict(maxima='carg',
                                 mathematica='Arg',
                                 sympy='arg',
                                 giac='arg'))

    def _eval_(self, x):
        """
        EXAMPLES::

            sage: arg(3+i)
            arctan(1/3)
            sage: arg(-1+i)
            3/4*pi
            sage: arg(2+2*i)
            1/4*pi
            sage: arg(2+x)
            arg(x + 2)
            sage: arg(2.0+i+x)
            arg(x + 2.00000000000000 + 1.00000000000000*I)
            sage: arg(-3)
            pi
            sage: arg(3)
            0
            sage: arg(0)
            0
            sage: arg(sqrt(2)+i)
            arg(sqrt(2) + I)

        """
        if isinstance(x,Expression):
            if x.is_trivial_zero():
                return x
        else:
            if not x:
                return x
            else:
                return arctan2(imag_part(x),real_part(x))

    def _evalf_(self, x, parent=None, algorithm=None):
        """
        EXAMPLES::

            sage: arg(0.0)
            0.000000000000000
            sage: arg(3.0)
            0.000000000000000
            sage: arg(3.00000000000000000000000000)
            0.00000000000000000000000000
            sage: arg(3.00000000000000000000000000).prec()
            90
            sage: arg(ComplexIntervalField(90)(3)).prec()
            90
            sage: arg(ComplexIntervalField(90)(3)).parent()
            Real Interval Field with 90 bits of precision
            sage: arg(3.0r)
            0.0
            sage: arg(RDF(3))
            0.0
            sage: arg(RDF(3)).parent()
            Real Double Field
            sage: arg(-2.5)
            3.14159265358979
            sage: arg(2.0+3*i)
            0.982793723247329

        TESTS:

        Make sure that the ``_evalf_`` method works when it receives a
        keyword argument ``parent`` :trac:`12289`::

            sage: arg(5+I, hold=True).n()
            0.197395559849881
        """
        try:
            return x.arg()
        except AttributeError:
            pass
        # try to find a parent that support .arg()
        if parent is None:
            parent = s_parent(x)
        try:
            parent = parent.complex_field()
        except AttributeError:
            try:
                parent = ComplexField(x.prec())
            except AttributeError:
                parent = ComplexField()

        return parent(x).arg()

arg=Function_arg()


############################
# Real and Imaginary Parts #
############################
class Function_real_part(GinacFunction):
    def __init__(self):
        r"""
        Returns the real part of the (possibly complex) input.

        It is possible to prevent automatic evaluation using the
        ``hold`` parameter::

            sage: real_part(I,hold=True)
            real_part(I)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: real_part(I,hold=True).simplify()
            0

        EXAMPLES::

            sage: z = 1+2*I
            sage: real(z)
            1
            sage: real(5/3)
            5/3
            sage: a = 2.5
            sage: real(a)
            2.50000000000000
            sage: type(real(a))
            <class 'sage.rings.real_mpfr.RealLiteral'>
            sage: real(1.0r)
            1.0
            sage: real(complex(3, 4))
            3.0

        Sage can recognize some expressions as real and accordingly
        return the identical argument::

            sage: SR.var('x', domain='integer').real_part()
            x
            sage: SR.var('x', domain='integer').imag_part()
            0
            sage: real_part(sin(x)+x)
            x + sin(x)
            sage: real_part(x*exp(x))
            x*e^x
            sage: imag_part(sin(x)+x)
            0
            sage: real_part(real_part(x))
            x
            sage: forget()

        TESTS::

            sage: loads(dumps(real_part))
            real_part
            sage: real_part(x)._sympy_()
            re(x)

        Check if :trac:`6401` is fixed::

            sage: latex(x.real())
            \Re \left( x \right)

            sage: f(x) = function('f')(x)
            sage: latex( f(x).real())
            \Re \left( f\left(x\right) \right)

        Check that some real part expansions evaluate correctly
        (:trac:`21614`)::

            sage: real(sqrt(sin(x))).subs(x==0)
            0
        """
        GinacFunction.__init__(self, "real_part",
                               conversions=dict(maxima='realpart',
                                                sympy='re',
                                                giac='re'),
                               alt_name="real")

    def __call__(self, x, **kwargs):
        r"""
        TESTS::

            sage: type(real(complex(3, 4)))
            <... 'float'>
        """
        if isinstance(x, complex):
            return x.real
        else:
            return GinacFunction.__call__(self, x, **kwargs)

real = real_part = Function_real_part()


class Function_imag_part(GinacFunction):
    def __init__(self):
        r"""
        Returns the imaginary part of the (possibly complex) input.

        It is possible to prevent automatic evaluation using the
        ``hold`` parameter::

            sage: imag_part(I,hold=True)
            imag_part(I)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: imag_part(I,hold=True).simplify()
            1

        TESTS::

            sage: z = 1+2*I
            sage: imaginary(z)
            2
            sage: imag(z)
            2
            sage: imag(complex(3, 4))
            4.0
            sage: loads(dumps(imag_part))
            imag_part
            sage: imag_part(x)._sympy_()
            im(x)

        Check if :trac:`6401` is fixed::

            sage: latex(x.imag())
            \Im \left( x \right)

            sage: f(x) = function('f')(x)
            sage: latex( f(x).imag())
            \Im \left( f\left(x\right) \right)
        """
        GinacFunction.__init__(self, "imag_part",
                               conversions=dict(maxima='imagpart',
                                                sympy='im',
                                                giac='im'),
                               alt_name="imag")

    def __call__(self, x, **kwargs):
        r"""
        TESTS::

            sage: type(imag(complex(3, 4)))
            <... 'float'>
        """
        if isinstance(x, complex):
            return x.imag
        else:
            return GinacFunction.__call__(self, x, **kwargs)

imag = imag_part = imaginary = Function_imag_part()


############################
# Complex Conjugate        #
############################
class Function_conjugate(GinacFunction):
    def __init__(self):
        r"""
        Returns the complex conjugate of the input.

        It is possible to prevent automatic evaluation using the
        ``hold`` parameter::

            sage: conjugate(I,hold=True)
            conjugate(I)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: conjugate(I,hold=True).simplify()
            -I

        TESTS::

            sage: x,y = var('x,y')
            sage: x.conjugate()
            conjugate(x)
            sage: _._sympy_()
            conjugate(x)
            sage: latex(conjugate(x))
            \overline{x}
            sage: f = function('f')
            sage: latex(f(x).conjugate())
            \overline{f\left(x\right)}
            sage: f = function('psi')(x,y)
            sage: latex(f.conjugate())
            \overline{\psi\left(x, y\right)}
            sage: x.conjugate().conjugate()
            x
            sage: x.conjugate().operator()
            conjugate
            sage: x.conjugate().operator() == conjugate
            True

        Check if :trac:`8755` is fixed::

            sage: conjugate(sqrt(-3))
            conjugate(sqrt(-3))
            sage: conjugate(sqrt(3))
            sqrt(3)
            sage: conjugate(sqrt(x))
            conjugate(sqrt(x))
            sage: conjugate(x^2)
            conjugate(x)^2
            sage: var('y',domain='positive')
            y
            sage: conjugate(sqrt(y))
            sqrt(y)

        Check if :trac:`10964` is fixed::

            sage: z= I*sqrt(-3); z
            I*sqrt(-3)
            sage: conjugate(z)
            -I*conjugate(sqrt(-3))
            sage: var('a')
            a
            sage: conjugate(a*sqrt(-2)*sqrt(-3))
            conjugate(sqrt(-2))*conjugate(sqrt(-3))*conjugate(a)

        Check that sums are handled correctly::

            sage: y = var('y', domain='real')
            sage: conjugate(y + I)
            y - I

        Test pickling::

            sage: loads(dumps(conjugate))
            conjugate
        """
        GinacFunction.__init__(self, "conjugate",
                               conversions=dict(sympy='conjugate',
                                                giac='conj'))

conjugate = Function_conjugate()


class Function_factorial(GinacFunction):
    def __init__(self):
        r"""
        Returns the factorial of `n`.

        INPUT:

        -  ``n`` - a non-negative integer, a complex number (except negative
           integers) or any symbolic expression


        OUTPUT: an integer or symbolic expression

        EXAMPLES::

            sage: factorial(0)
            1
            sage: factorial(4)
            24
            sage: factorial(10)
            3628800
            sage: factorial(6) == 6*5*4*3*2
            True

            sage: x = SR.var('x')
            sage: f = factorial(x + factorial(x)); f
            factorial(x + factorial(x))
            sage: f(x=3)
            362880
            sage: factorial(x)^2
            factorial(x)^2

        To prevent automatic evaluation use the ``hold`` argument::

            sage: factorial(5, hold=True)
            factorial(5)

        To then evaluate again, we currently must use Maxima via
        :meth:`sage.symbolic.expression.Expression.simplify`::

            sage: factorial(5, hold=True).simplify()
            120

        We can also give input other than nonnegative integers.  For
        other nonnegative numbers, the :func:`sage.functions.gamma.gamma`
        function is used::

            sage: factorial(1/2)
            1/2*sqrt(pi)
            sage: factorial(3/4)
            gamma(7/4)
            sage: factorial(2.3)
            2.68343738195577

        But negative input always fails::

            sage: factorial(-32)
            Traceback (most recent call last):
            ...
            ValueError: factorial only defined for non-negative integers

        And very large integers remain unevaluated::

            sage: factorial(2**64)
            factorial(18446744073709551616)
            sage: SR(2**64).factorial()
            factorial(18446744073709551616)

        TESTS:

        We verify that we can convert this function to Maxima and
        bring it back into Sage.::

            sage: z = var('z')
            sage: factorial._maxima_init_()
            'factorial'
            sage: maxima(factorial(z))
            factorial(_SAGE_VAR_z)
            sage: _.sage()
            factorial(z)
            sage: _._sympy_()
            factorial(z)
            sage: k = var('k')
            sage: factorial(k)
            factorial(k)

            sage: factorial(3.14)
            7.173269190187...

        Test latex typesetting::

            sage: latex(factorial(x))
            x!
            sage: latex(factorial(2*x))
            \left(2 \, x\right)!
            sage: latex(factorial(sin(x)))
            \sin\left(x\right)!
            sage: latex(factorial(sqrt(x+1)))
            \left(\sqrt{x + 1}\right)!
            sage: latex(factorial(sqrt(x)))
            \sqrt{x}!
            sage: latex(factorial(x^(2/3)))
            \left(x^{\frac{2}{3}}\right)!

            sage: latex(factorial)
            {\rm factorial}

        Check that :trac:`11539` is fixed::

            sage: (factorial(x) == 0).simplify()
            factorial(x) == 0
            sage: maxima(factorial(x) == 0).sage()
            factorial(x) == 0
            sage: y = var('y')
            sage: (factorial(x) == y).solve(x)
            [factorial(x) == y]

        Check that :trac:`16166` is fixed::

            sage: RBF = RealBallField(53)
            sage: factorial(RBF(4.2)) # abs tol 1e-13
            [32.5780960503314 +/- 6.06e-14]

        Test pickling::

            sage: loads(dumps(factorial))
            factorial
        """
        GinacFunction.__init__(self, "factorial", latex_name='{\\rm factorial}',
                conversions=dict(maxima='factorial',
                                 mathematica='Factorial',
                                 sympy='factorial',
                                 fricas='factorial',
                                 giac='factorial'))

    def _eval_(self, x):
        """
        Evaluate the factorial function.

        Note that this method overrides the eval method defined in GiNaC
        which calls numeric evaluation on all numeric input. We preserve
        exact results if the input is a rational number.

        EXAMPLES::

            sage: k = var('k')
            sage: k.factorial()
            factorial(k)
            sage: SR(1/2).factorial()
            1/2*sqrt(pi)
            sage: SR(3/4).factorial()
            gamma(7/4)
            sage: SR(5).factorial()
            120
            sage: SR(3245908723049857203948572398475r).factorial()
            factorial(3245908723049857203948572398475)
            sage: SR(3245908723049857203948572398475).factorial()
            factorial(3245908723049857203948572398475)

        TESTS:

        Check that :trac:`25421` is fixed::

            sage: factorial(RBF(2)**64)
            [+/- 2.30e+347382171326740403407]

        Check that :trac:`26749` is fixed::

            sage: factorial(float(3.2))        # abs tol 1e-14
            7.7566895357931776
            sage: type(factorial(float(3.2)))
            <class 'float'>
        """
        if isinstance(x, Integer):
            try:
                return x.factorial()
            except OverflowError:
                return
        elif isinstance(x, Rational):
            from sage.functions.gamma import gamma
            return gamma(x + 1)
        elif isinstance(x, Element) and hasattr(x.parent(), 'precision'):
            return (x + 1).gamma()
        elif self._is_numerical(x):
            from sage.functions.gamma import gamma
            return gamma(x + 1)

factorial = Function_factorial()


class Function_binomial(GinacFunction):
    def __init__(self):
        r"""
        Return the binomial coefficient

        .. MATH::

            \binom{x}{m} = x (x-1) \cdots (x-m+1) / m!


        which is defined for `m \in \ZZ` and any
        `x`. We extend this definition to include cases when
        `x-m` is an integer but `m` is not by

        .. MATH::

            \binom{x}{m}= \binom{x}{x-m}

        If `m < 0`, return `0`.

        INPUT:

        -  ``x``, ``m`` - numbers or symbolic expressions. Either ``m``
           or ``x-m`` must be an integer, else the output is symbolic.

        OUTPUT: number or symbolic expression (if input is symbolic)

        EXAMPLES::

            sage: binomial(5,2)
            10
            sage: binomial(2,0)
            1
            sage: binomial(1/2, 0)
            1
            sage: binomial(3,-1)
            0
            sage: binomial(20,10)
            184756
            sage: binomial(-2, 5)
            -6
            sage: binomial(RealField()('2.5'), 2)
            1.87500000000000
            sage: n=var('n'); binomial(n,2)
            1/2*(n - 1)*n
            sage: n=var('n'); binomial(n,n)
            1
            sage: n=var('n'); binomial(n,n-1)
            n
            sage: binomial(2^100, 2^100)
            1

        ::

            sage: k, i = var('k,i')
            sage: binomial(k,i)
            binomial(k, i)

        We can use a ``hold`` parameter to prevent automatic evaluation::

            sage: SR(5).binomial(3, hold=True)
            binomial(5, 3)
            sage: SR(5).binomial(3, hold=True).simplify()
            10

        TESTS:

        We verify that we can convert this function to Maxima and
        bring it back into Sage.

        ::

            sage: n,k = var('n,k')
            sage: maxima(binomial(n,k))
            binomial(_SAGE_VAR_n,_SAGE_VAR_k)
            sage: _.sage()
            binomial(n, k)
            sage: _._sympy_()
            binomial(n, k)
            sage: binomial._maxima_init_()
            'binomial'

        For polynomials::

            sage: y = polygen(QQ, 'y')
            sage: binomial(y, 2).parent()
            Univariate Polynomial Ring in y over Rational Field

        :trac:`16726`::

            sage: binomial(CIF(1), 2)
            0
            sage: binomial(CIF(3), 2)
            3

        Test pickling::

            sage: loads(dumps(binomial(n,k)))
            binomial(n, k)
        """
        GinacFunction.__init__(self, "binomial", nargs=2, preserved_arg=1,
                conversions=dict(maxima='binomial',
                                 mathematica='Binomial',
                                 sympy='binomial',
                                 fricas='binomial',
                                 giac='comb'))

    def _binomial_sym(self, n, k):
        """
        Expand the binomial formula symbolically when the second argument
        is an integer.

        EXAMPLES::

            sage: binomial._binomial_sym(x, 3)
            1/6*(x - 1)*(x - 2)*x
            sage: binomial._binomial_sym(x, x)
            Traceback (most recent call last):
            ...
            ValueError: second argument must be an integer
            sage: binomial._binomial_sym(x, SR(3))
            1/6*(x - 1)*(x - 2)*x

            sage: binomial._binomial_sym(x, 0r)
            1
            sage: binomial._binomial_sym(x, -1)
            0

            sage: y = polygen(QQ, 'y')
            sage: binomial._binomial_sym(y, 2).parent()
            Univariate Polynomial Ring in y over Rational Field
        """
        if isinstance(k, Expression):
            if k.is_integer():
                k = k.pyobject()
            else:
                raise ValueError("second argument must be an integer")

        if k < 0:
            return s_parent(k)(0)
        if k == 0:
            return s_parent(k)(1)
        if k == 1:
            return n

        from sage.misc.misc_c import prod
        return prod(n - i for i in range(k)) / factorial(k)

    def _eval_(self, n, k):
        """
        EXAMPLES::

            sage: binomial._eval_(5, 3)
            10
            sage: type(binomial._eval_(5, 3))
            <class 'sage.rings.integer.Integer'>
            sage: type(binomial._eval_(5., 3))
            <class 'sage.rings.real_mpfr.RealNumber'>
            sage: binomial._eval_(x, 3)
            1/6*(x - 1)*(x - 2)*x
            sage: binomial._eval_(x, x-2)
            1/2*(x - 1)*x
            sage: n = var('n')
            sage: binomial._eval_(x, n) is None
            True

        TESTS::

            sage: y = polygen(QQ, 'y')
            sage: binomial._eval_(y, 2).parent()
            Univariate Polynomial Ring in y over Rational Field
        """
        if not isinstance(k, Expression):
            if not isinstance(n, Expression):
                n, k = coercion_model.canonical_coercion(n, k)
                return self._evalf_(n, k)
        if k in ZZ:
            return self._binomial_sym(n, k)
        if (n - k) in ZZ:
            return self._binomial_sym(n, n - k)

        return None

    def _evalf_(self, n, k, parent=None, algorithm=None):
        """
        EXAMPLES::

            sage: binomial._evalf_(5.r, 3)
            10.0
            sage: type(binomial._evalf_(5.r, 3))
            <... 'float'>
            sage: binomial._evalf_(1/2,1/1)
            1/2
            sage: binomial._evalf_(10^20+1/1,10^20)
            100000000000000000001
            sage: binomial._evalf_(SR(10**7),10**7)
            1
            sage: binomial._evalf_(3/2,SR(1/1))
            3/2
        """
        return arith_binomial(n, k)

binomial = Function_binomial()


class Function_sum(BuiltinFunction):
    """
    Placeholder symbolic sum function that is only accessible internally.

    EXAMPLES::

        sage: from sage.functions.other import symbolic_sum as ssum
        sage: r = ssum(x, x, 1, 10); r
        sum(x, x, 1, 10)
        sage: r.unhold()
        55
    """
    def __init__(self):
        """
        EXAMPLES::

            sage: from sage.functions.other import symbolic_sum as ssum
            sage: maxima(ssum(x, x, 1, 10))
            55
        """
        BuiltinFunction.__init__(self, "sum", nargs=4,
                               conversions=dict(maxima='sum'))

    def _print_latex_(self, x, var, a, b):
        r"""
        EXAMPLES::

            sage: from sage.functions.other import symbolic_sum as ssum
            sage: latex(ssum(x^2, x, 1, 10))
            {\sum_{x=1}^{10} x^{2}}
        """
        return r"{{\sum_{{{}={}}}^{{{}}} {}}}".format(latex(var), latex(a),
                                                      latex(b), latex(x))

    def _sympy_(self, term, k, a, n):
        """
        Convert to sympy Sum.

        EXAMPLES::

            sage: var('k, n')
            (k, n)
            sage: s = sum(k, k, 1, n, hold=True)
            sage: s
            sum(k, k, 1, n)
            sage: s._sympy_() # indirect test
            Sum(k, (k, 1, n))
            sage: s._sympy_().doit()
            n**2/2 + n/2

        """
        import sympy
        return sympy.Sum(term, (k, a, n))

symbolic_sum = Function_sum()


class Function_prod(BuiltinFunction):
    """
    Placeholder symbolic product function that is only accessible internally.

    EXAMPLES::

        sage: from sage.functions.other import symbolic_product as sprod
        sage: r = sprod(x, x, 1, 10); r
        product(x, x, 1, 10)
        sage: r.unhold()
        3628800
    """
    def __init__(self):
        """
        EXAMPLES::

            sage: from sage.functions.other import symbolic_product as sprod
            sage: _ = var('m n', domain='integer')
            sage: r = maxima(sprod(sin(m), m, 1, n)).sage(); r
            product(sin(m), m, 1, n)
            sage: isinstance(r.operator(), sage.functions.other.Function_prod)
            True
            sage: r = sympy(sprod(sin(m), m, 1, n)).sage(); r # known bug
            product(sin(m), m, 1, n)
            sage: isinstance(r.operator(),
            ....:     sage.functions.other.Function_prod) # known bug
            True
            sage: giac(sprod(m, m, 1, n)).sage()
            factorial(n)
        """
        BuiltinFunction.__init__(self, "product", nargs=4,
                               conversions=dict(maxima='product',
                                   sympy='Product', giac='product'))

    def _print_latex_(self, x, var, a, b):
        r"""
        EXAMPLES::

            sage: from sage.functions.other import symbolic_product as sprod
            sage: latex(sprod(x^2, x, 1, 10))
            {\prod_{x=1}^{10} x^{2}}
        """
        return r"{{\prod_{{{}={}}}^{{{}}} {}}}".format(latex(var), latex(a),
                                                       latex(b), latex(x))

    def _sympy_(self, term, k, a, n):
        """
        Convert to sympy Product.

        EXAMPLES::

            sage: var('k, n')
            (k, n)
            sage: p = product(k^2+k+1,k,1,n, hold=True)
            sage: p._sympy_() # indirect test
            Product(k**2 + k + 1, (k, 1, n))
        """
        import sympy
        return sympy.Product(term, (k, a, n))

symbolic_product = Function_prod()


class Function_limit(BuiltinFunction):
    """
    Placeholder symbolic limit function that is only accessible internally.

    This function is called to create formal wrappers of limits that
    Maxima can't compute::

        sage: a = lim(exp(x^2)*(1-erf(x)), x=infinity); a
        -limit((erf(x) - 1)*e^(x^2), x, +Infinity)

    EXAMPLES::

        sage: from sage.functions.other import symbolic_limit as slimit
        sage: slimit(1/x, x, +oo)
        limit(1/x, x, +Infinity)
        sage: var('minus,plus')
        (minus, plus)
        sage: slimit(1/x, x, +oo)
        limit(1/x, x, +Infinity)
        sage: slimit(1/x, x, 0, plus)
        limit(1/x, x, 0, plus)
        sage: slimit(1/x, x, 0, minus)
        limit(1/x, x, 0, minus)
    """
    def __init__(self):
        """
        EXAMPLES::

            sage: from sage.functions.other import symbolic_limit as slimit
            sage: maxima(slimit(1/x, x, +oo))
            0
        """
        BuiltinFunction.__init__(self, "limit", nargs=0,
                               conversions=dict(maxima='limit'))

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: from sage.functions.other import symbolic_limit as slimit
            sage: latex(slimit)
            \lim
        """
        return r'\lim'

    def _print_latex_(self, ex, var, to, direction=''):
        r"""
        EXAMPLES::

            sage: from sage.functions.other import symbolic_limit as slimit
            sage: var('x,a')
            (x, a)
            sage: f = function('f')
            sage: latex(slimit(f(x), x, a))
            \lim_{x \to a}\, f\left(x\right)
            sage: latex(limit(f(x), x=oo))
            \lim_{x \to +\infty}\, f\left(x\right)

        TESTS:

        When one-sided limits are converted back from maxima, the direction
        argument becomes a symbolic variable. We check if typesetting these works::

            sage: from sage.functions.other import symbolic_limit as slimit
            sage: var('minus,plus')
            (minus, plus)
            sage: latex(slimit(f(x), x, a, minus))
            \lim_{x \to a^-}\, f\left(x\right)
            sage: latex(slimit(f(x), x, a, plus))
            \lim_{x \to a^+}\, f\left(x\right)
            sage: latex(limit(f(x),x=a,dir='+'))
            \lim_{x \to a^+}\, f\left(x\right)
            sage: latex(limit(f(x),x=a,dir='right'))
            \lim_{x \to a^+}\, f\left(x\right)
            sage: latex(limit(f(x),x=a,dir='-'))
            \lim_{x \to a^-}\, f\left(x\right)
            sage: latex(limit(f(x),x=a,dir='left'))
            \lim_{x \to a^-}\, f\left(x\right)

        Check if :trac:`13181` is fixed::

            sage: t = var('t')
            sage: latex(limit(exp_integral_e(1/2, I*t - I*x)*sqrt(-t + x),t=x,dir='-'))
            \lim_{t \to x^-}\, \sqrt{-t + x} E_{\frac{1}{2}}\left(i \, t - i \, x\right)
            sage: latex(limit(exp_integral_e(1/2, I*t - I*x)*sqrt(-t + x),t=x,dir='+'))
            \lim_{t \to x^+}\, \sqrt{-t + x} E_{\frac{1}{2}}\left(i \, t - i \, x\right)
            sage: latex(limit(exp_integral_e(1/2, I*t - I*x)*sqrt(-t + x),t=x))
            \lim_{t \to x}\, \sqrt{-t + x} E_{\frac{1}{2}}\left(i \, t - i \, x\right)
        """
        if repr(direction) == 'minus':
            dir_str = '^-'
        elif repr(direction) == 'plus':
            dir_str = '^+'
        else:
            dir_str = ''
        return r"\lim_{{{} \to {}{}}}\, {}".format(latex(var),
                latex(to), dir_str, latex(ex))

symbolic_limit = Function_limit()


class Function_cases(GinacFunction):
    """
    Formal function holding ``(condition, expression)`` pairs.

    Numbers are considered conditions with zero being ``False``.
    A true condition marks a default value. The function is not
    evaluated as long as it contains a relation that cannot be
    decided by Pynac.

    EXAMPLES::

        sage: ex = cases([(x==0, pi), (True, 0)]); ex
        cases(((x == 0, pi), (1, 0)))
        sage: ex.subs(x==0)
        pi
        sage: ex.subs(x==2)
        0
        sage: ex + 1
        cases(((x == 0, pi), (1, 0))) + 1
        sage: _.subs(x==0)
        pi + 1

    The first encountered default is used, as well as the first relation
    that can be trivially decided::

        sage: cases(((True, pi), (True, 0)))
        pi

        sage: _ = var('y')
        sage: ex = cases(((x==0, pi), (y==1, 0))); ex
        cases(((x == 0, pi), (y == 1, 0)))
        sage: ex.subs(x==0)
        pi
        sage: ex.subs(x==0, y==1)
        pi
    """
    def __init__(self):
        """
        EXAMPLES::

            sage: loads(dumps(cases))
            cases
        """
        GinacFunction.__init__(self, "cases")

    def __call__(self, l, **kwargs):
        """
        EXAMPLES::

            sage: ex = cases([(x==0, pi), (True, 0)]); ex
            cases(((x == 0, pi), (1, 0)))

        TESTS::

            sage: cases()
            Traceback (most recent call last):
            ...
            TypeError: ...__call__() missing 1 required positional argument: 'l'

            sage: cases(x)
            Traceback (most recent call last):
            ...
            RuntimeError: cases argument not a sequence
        """
        return GinacFunction.__call__(self,
                SR._force_pyobject(l), **kwargs)

    def _print_latex_(self, l, **kwargs):
        r"""
        EXAMPLES::

            sage: ex = cases([(x==0, pi), (True, 0)]); ex
            cases(((x == 0, pi), (1, 0)))
            sage: latex(ex)
            \begin{cases}{\pi} & {x = 0}\\{0} & {1}\end{cases}

        TESTS:

        Verify that :trac:`25624` is fixed::

            sage: L = latex(cases([(x == 0, 0), (1, 1)]))
            sage: L
            \begin{cases}{0} & {x = 0}\\{1} & {1}\end{cases}
        """
        if not isinstance(l, (list, tuple)):
            raise ValueError("cases() argument must be a list")
        str = r"\begin{cases}"
        for pair in l:
            left = None
            if (isinstance(pair, tuple)):
                right,left = pair
            else:
                right = pair
            str += r"{%s} & {%s}\\" % (latex(left), latex(right))
        return str[:-2] + r"\end{cases}"

    def _sympy_(self, l):
        """
        Convert this cases expression to its SymPy equivalent.

        EXAMPLES::

            sage: ex = cases(((x<0, pi), (x==1, 1), (True, 0)))
            sage: assert ex == ex._sympy_()._sage_()
        """
        from sage.symbolic.ring import SR
        from sympy import Piecewise as pw
        args = []
        for tup in l.operands():
            cond,expr = tup.operands()
            if SR(cond).is_numeric():
                args.append((SR(expr)._sympy_(), bool(SR(cond)._sympy_())))
            else:
                args.append((SR(expr)._sympy_(), SR(cond)._sympy_()))
        return pw(*args)

cases = Function_cases()


class Function_crootof(BuiltinFunction):
    """
    Formal function holding ``(polynomial, index)`` pairs.

    The expression evaluates to a floating point value that is an
    approximation to a specific complex root of the polynomial. The
    ordering is fixed so you always get the same root.

    The functionality is imported from SymPy, see
    http://docs.sympy.org/latest/_modules/sympy/polys/rootoftools.html

    EXAMPLES::

        sage: c = complex_root_of(x^6 + x + 1, 1); c
        complex_root_of(x^6 + x + 1, 1)
        sage: c.n()
        -0.790667188814418 + 0.300506920309552*I
        sage: c.n(100)
        -0.79066718881441764449859281847 + 0.30050692030955162512001002521*I
        sage: (c^6 + c + 1).n(100) < 1e-25
        True
    """
    def __init__(self):
        """
        EXAMPLES::

            sage: loads(dumps(complex_root_of))
            complex_root_of
        """
        BuiltinFunction.__init__(self, "complex_root_of", nargs=2,
                                   conversions=dict(sympy='CRootOf'),
                                   evalf_params_first=False)

    def _eval_(self, poly, index):
        """
        TESTS::

            sage: _ = var('y')
            sage: complex_root_of(1, 1)
            Traceback (most recent call last):
            ...
            ValueError: polynomial in one variable required
            sage: complex_root_of(x+y, 1)
            Traceback (most recent call last):
            ...
            ValueError: polynomial in one variable required
            sage: complex_root_of(sin(x), 1)
            Traceback (most recent call last):
            ...
            ValueError: polynomial in one variable required
        """
        try:
            vars = poly.variables()
        except AttributeError:
            raise ValueError('polynomial in one variable required')
        if len(vars) != 1 or not poly.is_polynomial(vars[0]):
            raise ValueError('polynomial in one variable required')

    def _evalf_(self, poly, index, parent=None, algorithm=None):
        """
        EXAMPLES::

            sage: complex_root_of(x^2-2, 1).n()
            1.41421356237309
            sage: complex_root_of(x^2-2, 3).n()
            Traceback (most recent call last):
            ...
            IndexError: root index out of [-2, 1] range, got 3

        TESTS:

        Check that low precision is handled (:trac:`24378`)::

            sage: complex_root_of(x^8-1, 7).n(2)
            0.75 + 0.75*I
            sage: complex_root_of(x^8-1, 7).n(20)
            0.70711 + 0.70711*I
        """
        from sympy.core.evalf import prec_to_dps
        from sympy.polys import CRootOf, Poly
        try:
            prec = parent.precision()
        except AttributeError:
            prec = 53
        sobj = CRootOf(Poly(poly._sympy_()), int(index))
        return parent(sobj.n(1 + prec_to_dps(prec))._sage_())

complex_root_of = Function_crootof()


class Function_elementof(BuiltinFunction):
    """
    Formal set membership function that is only accessible internally.

    This function is called to express a set membership statement,
    usually as part of a solution set returned by ``solve()``.
    See :class:`sage.sets.set.Set` and :class:`sage.sets.real_set.RealSet`
    for possible set arguments.

    EXAMPLES::

        sage: from sage.functions.other import element_of
        sage: element_of(x, SR(ZZ))
        element_of(x, Integer Ring)
        sage: element_of(sin(x), SR(QQ))
        element_of(sin(x), Rational Field)
        sage: element_of(x, SR(RealSet.open_closed(0,1)))
        element_of(x, (0, 1])
        sage: element_of(x, SR(Set([4,6,8])))
        element_of(x, {8, 4, 6})
    """
    def __init__(self):
        """
        EXAMPLES::

            sage: from sage.functions.other import element_of
            sage: loads(dumps(element_of))
            element_of
        """
        BuiltinFunction.__init__(self, "element_of", nargs=2,
                                 conversions=dict(sympy='Contains'))

    def _eval_(self, x, s):
        """
        EXAMPLES::

            sage: from sage.functions.other import element_of
            sage: element_of(x, SR(RealSet(-oo, oo)))
            element_of(x, (-oo, +oo))
            sage: element_of(x, 0)
            Traceback (most recent call last):
            ...
            ValueError: not a set: 0
        """
        from sage.categories.sets_cat import Sets
        if s not in Sets():
            raise ValueError("not a set: {}".format(s))

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: from sage.functions.other import element_of
            sage: latex(element_of)
            \in
        """
        return r'\in'

    def _print_latex_(self, ex, s):
        r"""
        EXAMPLES::

            sage: from sage.functions.other import element_of
            sage: latex(element_of(x, SR(ZZ)))
            x \in \Bold{Z}
            sage: latex(element_of(x, SR(Set([4,6,8]))))
            x \in \left\{8, 4, 6\right\}
        """
        return r"{} \in {}".format(latex(ex), latex(s))

element_of = Function_elementof()
