"""
Symbolic Integration
"""
# ****************************************************************************
#       Copyright (C) 2009 Golam Mortuza Hossain <gmhossain@gmail.com>
#       Copyright (C) 2010 Burcin Erocal <burcin@erocal.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************`
from sage.structure.element import Expression
from sage.symbolic.ring import SR
from sage.symbolic.function import BuiltinFunction

##################################################################
#  Table of available integration routines
##################################################################

import sage.symbolic.integration.external as external

# Add new integration routines to the dictionary below. This will make them
# accessible with the 'algorithm' keyword parameter of top level integrate().
available_integrators = {}
available_integrators['maxima'] = external.maxima_integrator
available_integrators['sympy'] = external.sympy_integrator
available_integrators['mathematica_free'] = external.mma_free_integrator
available_integrators['fricas'] = external.fricas_integrator
available_integrators['giac'] = external.libgiac_integrator
available_integrators['libgiac'] = external.libgiac_integrator

######################################################
#
# Class implementing symbolic integration
#
######################################################


class IndefiniteIntegral(BuiltinFunction):
    def __init__(self):
        """
        Class to represent an indefinite integral.

        EXAMPLES::

            sage: from sage.symbolic.integration.integral import indefinite_integral
            sage: indefinite_integral(log(x), x) #indirect doctest
            x*log(x) - x
            sage: indefinite_integral(x^2, x)
            1/3*x^3
            sage: indefinite_integral(4*x*log(x), x)
            2*x^2*log(x) - x^2
            sage: indefinite_integral(exp(x), 2*x)
            2*e^x

        TESTS:

        Check for :issue:`28913`::

            sage: # needs sage.libs.giac
            sage: Ex = (1-2*x^(1/3))^(3/4)/x
            sage: integrate(Ex, x, algorithm='giac')  # long time
            4*(-2*x^(1/3) + 1)^(3/4) + 6*arctan((-2*x^(1/3) + 1)^(1/4)) - 3*log((-2*x^(1/3) + 1)^(1/4) + 1) + 3*log(abs((-2*x^(1/3) + 1)^(1/4) - 1))

        Check for :issue:`29833`::

            sage: (x,a,b)=var('x a b')
            sage: assume(b > 0)
            sage: f = (exp((x-a)/b) + 1)**(-1)
            sage: (f*f).integrate(x, algorithm='mathematica_free') # optional -- internet
            -b*log(e^(a/b) + e^(x/b)) + x + b/(e^(-(a - x)/b) + 1)

        After :issue:`25119` we can integrate the following function,
        although giac and sympy give different-looking answers::

            sage: result = integrate(sqrt(x^2)/x,x)
            ...
            sage: result in [x*sgn(x), sqrt(x^2)]
            True
        """
        # The automatic evaluation routine will try these integrators
        # in the given order. This is an attribute of the class instead of
        # a global variable in this module to enable customization by
        # creating a subclasses which define a different set of integrators
        #
        # The libgiac integrator may immediately return a symbolic
        # (unevaluated) answer if libgiac is unavailable. This essentially
        # causes it to be skipped.
        self.integrators = [external.maxima_integrator,
                            external.libgiac_integrator,
                            external.sympy_integrator]

        BuiltinFunction.__init__(self, "integrate", nargs=2, conversions={'sympy': 'Integral',
                                                                          'giac': 'integrate'})

    def _eval_(self, f, x):
        """
        EXAMPLES::

            sage: from sage.symbolic.integration.integral import indefinite_integral
            sage: indefinite_integral(exp(x), x) # indirect doctest
            e^x
            sage: indefinite_integral(exp(x), x^2)
            2*(x - 1)*e^x

        TESTS:

        Check that :issue:`28842` is fixed::

            sage: integrate(1/(x^4 + x^3 + 1), x)
            integrate(1/(x^4 + x^3 + 1), x)

        Check that :issue:`32002` is fixed. This needs giac since only
        giac can integrate it in any case::

            sage: # needs sage.libs.giac
            sage: result = integral(2*min_symbolic(x,2*x),x)
            ...
            sage: result
            -1/2*x^2*sgn(x) + 3/2*x^2
        """
        # Check for x
        if not (isinstance(x, Expression) and x.is_symbol()):
            if len(x.variables()) == 1:
                nx = x.variables()[0]
                f = f * x.diff(nx)
                x = nx
            else:
                return None

        # we try all listed integration algorithms
        A = None
        for integrator in self.integrators:
            try:
                A = integrator(f, x)
            except (NotImplementedError, TypeError,
                    AttributeError, RuntimeError):
                pass
            except ValueError:
                # maxima is telling us something
                raise
            else:
                if not hasattr(A, 'operator'):
                    return A
                else:
                    uneval = integral(SR.wild(0), x, hold=True)
                    if not A.has(uneval):
                        return A
        return A

    def _tderivative_(self, f, x, diff_param=None):
        """
        EXAMPLES::

            sage: from sage.symbolic.integration.integral import indefinite_integral
            sage: f = function('f'); a,b=var('a,b')
            sage: h = indefinite_integral(f(x), x)
            sage: h.diff(x) # indirect doctest
            f(x)
            sage: h.diff(a)
            0
        """
        if x.has(diff_param):
            return f * x.derivative(diff_param)
        else:
            return f.derivative(diff_param).integral(x)

    def _print_latex_(self, f, x):
        r"""
        EXAMPLES::

            sage: from sage.symbolic.integration.integral import indefinite_integral
            sage: print_latex = indefinite_integral._print_latex_
            sage: var('x,a,b')
            (x, a, b)
            sage: f = function('f')
            sage: print_latex(f(x),x)
            '\\int f\\left(x\\right)\\,{d x}'
            sage: latex(integrate(tan(x)/x, x))
            \int \frac{\tan\left(x\right)}{x}\,{d x}
        """
        from sage.misc.latex import latex
        if not (isinstance(x, Expression) and x.is_symbol()):
            dx_str = "{d \\left(%s\\right)}" % latex(x)
        else:
            dx_str = "{d %s}" % latex(x)

        return "\\int %s\\,%s" % (latex(f), dx_str)


indefinite_integral = IndefiniteIntegral()


class DefiniteIntegral(BuiltinFunction):
    def __init__(self):
        """
        The symbolic function representing a definite integral.

        EXAMPLES::

            sage: from sage.symbolic.integration.integral import definite_integral
            sage: definite_integral(sin(x),x,0,pi)
            2

        TESTS:

        Check for :issue:`32354`::

            sage: # needs sage.libs.giac
            sage: ex = 1/max_symbolic(x, 1)**2
            sage: integral(ex, x, 0, 2, algorithm='giac')
            3/2
            sage: integral(1/max_symbolic(x, 1)**2, x, 0, oo, algorithm='giac')
            2
        """
        # The automatic evaluation routine will try these integrators
        # in the given order. This is an attribute of the class instead of
        # a global variable in this module to enable customization by
        # creating a subclasses which define a different set of integrators
        self.integrators = [external.maxima_integrator,
                            external.libgiac_integrator,
                            external.sympy_integrator]

        BuiltinFunction.__init__(self, "integrate", nargs=4, conversions={'sympy': 'Integral',
                                                                          'giac': 'integrate'})

    def _eval_(self, f, x, a, b):
        """
        Return the results of symbolic evaluation of the integral.

        EXAMPLES::

            sage: from sage.symbolic.integration.integral import definite_integral
            sage: definite_integral(exp(x),x,0,1) # indirect doctest
            e - 1

        TESTS:

        Check that :issue:`32002` is fixed::

            sage: result = integral(2*min_symbolic(x,2*x),x,-1,1)
            ...
            sage: result
            -1
        """
        # Check for x
        if not (isinstance(x, Expression) and x.is_symbol()):
            if len(x.variables()) == 1:
                nx = x.variables()[0]
                f = f * x.diff(nx)
                x = nx
            else:
                return None

        args = (f, x, a, b)

        # we try all listed integration algorithms
        A = None
        for integrator in self.integrators:
            try:
                A = integrator(*args)
            except (NotImplementedError, TypeError,
                    AttributeError, RuntimeError):
                pass
            except ValueError:
                # maxima is telling us something
                raise
            else:
                if not hasattr(A, 'operator'):
                    return A
                else:
                    uneval = integral(SR.wild(0), x, a, b, hold=True)
                    if not A.has(uneval):
                        return A
        return A

    def _evalf_(self, f, x, a, b, parent=None, algorithm=None):
        """
        Return a numerical approximation of the integral.

        EXAMPLES::

            sage: from sage.symbolic.integration.integral import definite_integral
            sage: f = sin(x)*log(x)/x^2
            sage: h = definite_integral(f, x, 1, 2, hold=True); h
            integrate(log(x)*sin(x)/x^2, x, 1, 2)
            sage: h.n() # indirect doctest
            0.14839875208053...

        TESTS:

        Check if :issue:`3863` is fixed::

            sage: integrate(x^2.7 * e^(-2.4*x), x, 0, 3).n()
            0.154572952320790
        """
        from sage.calculus.integration import numerical_integral
        # The gsl routine returns a tuple, which also contains the error.
        # We only return the result.
        return numerical_integral(f, a, b)[0]

    def _tderivative_(self, f, x, a, b, diff_param=None):
        """
        Return the derivative of symbolic integration.

        EXAMPLES::

            sage: from sage.symbolic.integration.integral import definite_integral
            sage: f = function('f'); a,b=var('a,b')
            sage: h = definite_integral(f(x), x,a,b)
            ...
            sage: h.diff(x) # indirect doctest
            0
            sage: h.diff(a)
            -f(a)
            sage: h.diff(b)
            f(b)

        TESTS:

        Check for :issue:`28656`::

            sage: t = var("t")
            sage: f = function("f")
            sage: F(x) = integrate(f(t),t,0,x)
            sage: F(x).diff(x)
            f(x)
        """
        if not x.has(diff_param):
            # integration variable != differentiation variable
            ans = definite_integral(f.diff(diff_param), x, a, b)
        else:
            ans = SR.zero()
        if hasattr(b, 'diff'):
            ans += f.subs(x == b) * b.diff(diff_param)
        if hasattr(a, 'diff'):
            ans -= f.subs(x == a) * a.diff(diff_param)
        return ans

    def _print_latex_(self, f, x, a, b):
        r"""
        Convert this integral to LaTeX notation.

        EXAMPLES::

            sage: from sage.symbolic.integration.integral import definite_integral
            sage: print_latex = definite_integral._print_latex_
            sage: var('x,a,b')
            (x, a, b)
            sage: f = function('f')
            sage: print_latex(f(x),x,0,1)
            '\\int_{0}^{1} f\\left(x\\right)\\,{d x}'
            sage: latex(integrate(tan(x)/x, x, 0, 1))
            \int_{0}^{1} \frac{\tan\left(x\right)}{x}\,{d x}
        """
        from sage.misc.latex import latex
        if not (isinstance(x, Expression) and x.is_symbol()):
            dx_str = "{d \\left(%s\\right)}" % latex(x)
        else:
            dx_str = "{d %s}" % latex(x)
        return "\\int_{%s}^{%s} %s\\,%s" % (latex(a), latex(b),
                                            latex(f), dx_str)

    def _sympy_(self, f, x, a, b):
        """
        Convert this integral to the equivalent SymPy object.

        The resulting SymPy integral can be evaluated using ``doit()``.

        EXAMPLES::

            sage: integral(x, x, 0, 1, hold=True)._sympy_()                             # needs sympy
            Integral(x, (x, 0, 1))
            sage: _.doit()                                                              # needs sympy
            1/2
        """
        from sympy.integrals import Integral
        return Integral(f, (x, a, b))


definite_integral = DefiniteIntegral()


def _normalize_integral_input(f, v, a=None, b=None):
    r"""
    Validate and return variable and endpoints for an integral.

    INPUT:

    - ``f`` -- an expression to integrate

    - ``v`` -- a variable of integration or a triple

    - ``a`` -- (optional) the left endpoint of integration

    - ``b`` -- (optional) the right endpoint of integration

    It is also possible to pass the last three parameters in ``v`` as a triple.

    If the input contains endpoints, both endpoints must be given.

    OUTPUT: a tuple of ``f``, ``v``, ``a``, and ``b``

    EXAMPLES::

        sage: from sage.symbolic.integration.integral import \
        ....:     _normalize_integral_input
        sage: _normalize_integral_input(x^2, x, 0, 3)
        (x^2, x, 0, 3)
        sage: _normalize_integral_input(x^2, [x, 0, 3], None, None)
        (x^2, x, 0, 3)
        sage: _normalize_integral_input(x^2, [x], None, None)
        (x^2, x, None, None)

    TESTS::

        sage: _normalize_integral_input(x^2, [0, 3], None, None)
        Traceback (most recent call last):
        ...
        TypeError: invalid input [0, 3] - please use variable,
        with or without two endpoints
        sage: _normalize_integral_input(x^2, x, 0, None)
        Traceback (most recent call last):
        ...
        TypeError: only one endpoint was given!
    """
    if isinstance(v, (list, tuple)) and a is None and b is None:
        if len(v) == 1:  # bare variable in a tuple
            v = v[0]
        elif len(v) == 3:  # variable and two endpoints
            v, a, b = v
        else:
            raise TypeError("invalid input %s - please use variable, "
                            "with or without two endpoints" % repr(v))

    if (a is None) ^ (b is None):
        raise TypeError('only one endpoint was given!')

    return f, v, a, b


def integrate(expression, v=None, a=None, b=None, algorithm=None, hold=False):
    r"""
    Return the indefinite integral with respect to the variable
    `v`, ignoring the constant of integration. Or, if endpoints
    `a` and `b` are specified, returns the definite
    integral over the interval `[a, b]`.

    If ``self`` has only one variable, then it returns the
    integral with respect to that variable.

    If definite integration fails, it could be still possible to
    evaluate the definite integral using indefinite integration with
    the Newton - Leibniz theorem (however, the user has to ensure that the
    indefinite integral is continuous on the compact interval `[a,b]` and
    this theorem can be applied).

    INPUT:

    - ``v`` -- a variable or variable name; this can also be a tuple of
      the variable (optional) and endpoints (i.e., ``(x,0,1)`` or ``(0,1)``)

    - ``a`` -- (optional) lower endpoint of definite integral

    - ``b`` -- (optional) upper endpoint of definite integral

    - ``algorithm`` -- (default: ``'maxima'``, ``'libgiac'`` and ``'sympy'``) one of

      - ``'maxima'`` -- use maxima

      - ``'sympy'`` -- use sympy (also in Sage)

      - ``'mathematica_free'`` -- use http://integrals.wolfram.com/

      - ``'fricas'`` -- use FriCAS (the optional fricas spkg has to be installed)

      - ``'giac'`` -- use libgiac

      - ``'libgiac'`` -- use libgiac (alias for ``'giac'``)

    To prevent automatic evaluation, use the ``hold`` argument.

    .. SEEALSO::

        To integrate a polynomial over a polytope, use the optional
        ``latte_int`` package
        :meth:`sage.geometry.polyhedron.base.Polyhedron_base.integrate`.

    EXAMPLES::

        sage: x = var('x')
        sage: h = sin(x)/(cos(x))^2
        sage: h.integral(x)
        1/cos(x)

    ::

        sage: f = x^2/(x+1)^3
        sage: f.integral(x)
        1/2*(4*x + 3)/(x^2 + 2*x + 1) + log(x + 1)

    ::

        sage: f = x*cos(x^2)
        sage: f.integral(x, 0, sqrt(pi))
        0
        sage: f.integral(x, a=-pi, b=pi)
        0

    ::

        sage: f(x) = sin(x)
        sage: f.integral(x, 0, pi/2)
        1

    The variable is required, but the endpoints are optional::

        sage: y = var('y')
        sage: integral(sin(x), x)
        -cos(x)
        sage: integral(sin(x), y)
        y*sin(x)
        sage: integral(sin(x), x, pi, 2*pi)
        -2
        sage: integral(sin(x), y, pi, 2*pi)
        pi*sin(x)
        sage: integral(sin(x), (x, pi, 2*pi))
        -2
        sage: integral(sin(x), (y, pi, 2*pi))
        pi*sin(x)

    Using the ``hold`` parameter it is possible to prevent automatic
    evaluation, which can then be evaluated via :meth:`simplify`::

        sage: integral(x^2, x, 0, 3)
        9
        sage: a = integral(x^2, x, 0, 3, hold=True) ; a
        integrate(x^2, x, 0, 3)
        sage: a.simplify()
        9

    Constraints are sometimes needed::

        sage: var('x, n')
        (x, n)
        sage: integral(x^n,x)
        Traceback (most recent call last):
        ...
        ValueError: Computation failed since Maxima requested additional
        constraints; using the 'assume' command before evaluation
        *may* help (example of legal syntax is 'assume(n>0)', see `assume?`
        for more details)
        Is n equal to -1?
        sage: assume(n > 0)
        sage: integral(x^n,x)
        x^(n + 1)/(n + 1)
        sage: forget()

    Usually the constraints are of sign, but others are possible::

        sage: assume(n==-1)
        sage: integral(x^n,x)
        log(x)

    Note that an exception is raised when a definite integral is
    divergent::

        sage: forget() # always remember to forget assumptions you no longer need
        sage: integrate(1/x^3,(x,0,1))
        Traceback (most recent call last):
        ...
        ValueError: Integral is divergent.
        sage: integrate(1/x^3,x,-1,3)
        Traceback (most recent call last):
        ...
        ValueError: Integral is divergent.

    But Sage can calculate the convergent improper integral of
    this function::

        sage: integrate(1/x^3,x,1,infinity)
        1/2

    The examples in the Maxima documentation::

        sage: var('x, y, z, b')
        (x, y, z, b)
        sage: integral(sin(x)^3, x)
        1/3*cos(x)^3 - cos(x)
        sage: integral(x/sqrt(b^2-x^2), b)
        x*log(2*b + 2*sqrt(b^2 - x^2))
        sage: integral(x/sqrt(b^2-x^2), x)
        -sqrt(b^2 - x^2)
        sage: integral(cos(x)^2 * exp(x), x, 0, pi)
        3/5*e^pi - 3/5
        sage: integral(x^2 * exp(-x^2), x, -oo, oo)
        1/2*sqrt(pi)

    We integrate the same function in both Mathematica and Sage (via
    Maxima)::

        sage: _ = var('x, y, z')
        sage: f = sin(x^2) + y^z
        sage: g = mathematica(f)                            # optional - mathematica
        sage: print(g)                                      # optional - mathematica
                  z        2
                 y  + Sin[x ]
        sage: print(g.Integrate(x))                         # optional - mathematica
                    z        Pi                2
                 x y  + Sqrt[--] FresnelS[Sqrt[--] x]
                             2                 Pi
        sage: print(f.integral(x))
        x*y^z + 1/16*sqrt(pi)*((I + 1)*sqrt(2)*erf((1/2*I + 1/2)*sqrt(2)*x)
                                + (I - 1)*sqrt(2)*erf((1/2*I - 1/2)*sqrt(2)*x)
                                - (I - 1)*sqrt(2)*erf(sqrt(-I)*x)
                                + (I + 1)*sqrt(2)*erf((-1)^(1/4)*x))

    Alternatively, just use algorithm='mathematica_free' to integrate via Mathematica
    over the internet (does NOT require a Mathematica license!)::

        sage: _ = var('x, y, z')  # optional - internet
        sage: f = sin(x^2) + y^z   # optional - internet
        sage: f.integrate(x, algorithm='mathematica_free')   # optional - internet
        x*y^z + sqrt(1/2)*sqrt(pi)*fresnel_sin(sqrt(2)*x/sqrt(pi))

    We can also use Sympy::

        sage: integrate(x*sin(log(x)), x)
        -1/5*x^2*(cos(log(x)) - 2*sin(log(x)))
        sage: integrate(x*sin(log(x)), x, algorithm='sympy')                            # needs sympy
        -1/5*x^2*cos(log(x)) + 2/5*x^2*sin(log(x))
        sage: _ = var('y, z')
        sage: (x^y - z).integrate(y)
        -y*z + x^y/log(x)
        sage: (x^y - z).integrate(y, algorithm='sympy')                                 # needs sympy
        -y*z + cases(((log(x) != 0, x^y/log(x)), (1, y)))

    We integrate the above function in Maple now::

        sage: g = maple(f); g.sort()         # optional - maple
        y^z+sin(x^2)
        sage: g.integrate(x).sort()          # optional - maple
        x*y^z+1/2*2^(1/2)*Pi^(1/2)*FresnelS(2^(1/2)/Pi^(1/2)*x)

    We next integrate a function with no closed form integral. Notice
    that the answer comes back as an expression that contains an
    integral itself. ::

        sage: A = integral(1/ ((x-4) * (x^4+x+1)), x); A
        integrate(1/((x^4 + x + 1)*(x - 4)), x)

    Sometimes, in this situation, using the algorithm "maxima"
    gives instead a partially integrated answer::

        sage: integral(1/(x**7-1),x,algorithm='maxima')
        -1/7*integrate((x^5 + 2*x^4 + 3*x^3 + 4*x^2 + 5*x + 6)/(x^6 + x^5 + x^4 + x^3 + x^2 + x + 1), x) + 1/7*log(x - 1)

    We now show that floats are not converted to rationals
    automatically since we by default have keepfloat: true in maxima.

    ::

        sage: integral(e^(-x^2),(x, 0, 0.1))
        0.05623145800914245*sqrt(pi)

    An example of an integral that fricas can integrate::

        sage: f(x) = sqrt(x+sqrt(1+x^2))/x
        sage: integrate(f(x), x, algorithm='fricas')      # optional - fricas
        2*sqrt(x + sqrt(x^2 + 1)) - 2*arctan(sqrt(x + sqrt(x^2 + 1)))
         - log(sqrt(x + sqrt(x^2 + 1)) + 1) + log(sqrt(x + sqrt(x^2 + 1)) - 1)

    where the default integrator obtains another answer::

        sage: integrate(f(x), x)  # long time
        1/8*sqrt(x)*gamma(1/4)*gamma(-1/4)^2*hypergeometric((-1/4, -1/4, 1/4),
                                                            (1/2, 3/4),
                                                            -1/x^2)/(pi*gamma(3/4))

    The following definite integral is not found by maxima::

        sage: f(x) = (x^4 - 3*x^2 + 6) / (x^6 - 5*x^4 + 5*x^2 + 4)
        sage: integrate(f(x), x, 1, 2, algorithm='maxima')
        integrate((x^4 - 3*x^2 + 6)/(x^6 - 5*x^4 + 5*x^2 + 4), x, 1, 2)

    but is nevertheless computed::

        sage: integrate(f(x), x, 1, 2)
        -1/2*pi + arctan(8) + arctan(5) + arctan(2) + arctan(1/2)

    Both fricas and sympy give the correct result::

        sage: integrate(f(x), x, 1, 2, algorithm='fricas')  # optional - fricas
        -1/2*pi + arctan(8) + arctan(5) + arctan(2) + arctan(1/2)
        sage: integrate(f(x), x, 1, 2, algorithm='sympy')                               # needs sympy
        -1/2*pi + arctan(8) + arctan(5) + arctan(2) + arctan(1/2)

    Using Giac to integrate the absolute value of a trigonometric
    expression. If Giac is installed, this will be attempted
    automatically in the event that Maxima is unable to integrate the
    expression::

        sage: # needs sage.libs.giac
        sage: result = integrate(abs(cos(x)), x, 0, 2*pi, algorithm='giac')
        ...
        sage: result
        4
        sage: result = integrate(abs(cos(x)), x, 0, 2*pi)
        ...
        sage: result
        4

    ALIASES: integral() and integrate() are the same.

    EXAMPLES:

    Here is an example where we have to use assume::

        sage: a,b = var('a,b')
        sage: integrate(1/(x^3 *(a+b*x)^(1/3)), x)
        Traceback (most recent call last):
        ...
        ValueError: Computation failed since Maxima requested additional
        constraints; using the 'assume' command before evaluation
        *may* help (example of legal syntax is 'assume(a>0)', see `assume?`
        for more details)
        Is a positive or negative?

    So we just assume that `a>0` and the integral works::

        sage: assume(a>0)
        sage: integrate(1/(x^3 *(a+b*x)^(1/3)), x)
        2/9*sqrt(3)*b^2*arctan(1/3*sqrt(3)*(2*(b*x + a)^(1/3) + a^(1/3))/a^(1/3))/a^(7/3)
         - 1/9*b^2*log((b*x + a)^(2/3) + (b*x + a)^(1/3)*a^(1/3) + a^(2/3))/a^(7/3)
         + 2/9*b^2*log((b*x + a)^(1/3) - a^(1/3))/a^(7/3) + 1/6*(4*(b*x + a)^(5/3)*b^2
         - 7*(b*x + a)^(2/3)*a*b^2)/((b*x + a)^2*a^2 - 2*(b*x + a)*a^3 + a^4)

    TESTS:

    The following integral was broken prior to Maxima 5.15.0 -
    see :issue:`3013`::

        sage: integrate(sin(x)*cos(10*x)*log(x), x)
        -1/198*(9*cos(11*x) - 11*cos(9*x))*log(x) + 1/44*Ei(11*I*x) - 1/36*Ei(9*I*x) - 1/36*Ei(-9*I*x) + 1/44*Ei(-11*I*x)

    It is no longer possible to use certain functions without an
    explicit variable.  Instead, evaluate the function at a variable,
    and then take the integral::

        sage: integrate(sin)
        Traceback (most recent call last):
        ...
        TypeError: unable to convert sin to a symbolic expression

        sage: integrate(sin(x), x)
        -cos(x)
        sage: integrate(sin(x), x, 0, 1)
        -cos(1) + 1

    Check if :issue:`780` is fixed::

        sage: _ = var('x,y')
        sage: f = log(x^2+y^2)
        sage: res = integral(f,x,1414/10^7, 1); res
        Traceback (most recent call last):
        ...
        ValueError: Computation failed since Maxima requested additional constraints; using the 'assume' command before evaluation *may* help ...
        Is ... positive, negative or zero?
        sage: assume(y>1)
        sage: res = integral(f,x,1414/10^7, 1); res
        -2*y*arctan(707/5000000/y) + 2*y*arctan(1/y) + log(y^2 + 1) - 707/5000000*log(y^2 + 499849/25000000000000) - 4999293/2500000
        sage: nres = numerical_integral(f.subs(y=2), 0.0001414, 1.); nres
        (1.4638323264144..., 1.6251803529759...e-14)
        sage: res.subs(y=2).n()
        1.46383232641443
        sage: nres = numerical_integral(f.subs(y=.5), 0.0001414, 1.); nres
        (-0.669511708872807, 7.768678110854711e-15)
        sage: res.subs(y=.5).n()
        -0.669511708872807

    Check if :issue:`6816` is fixed::

        sage: var('t,theta')
        (t, theta)
        sage: integrate(t*cos(-theta*t),t,0,pi)
        (pi*theta*sin(pi*theta) + cos(pi*theta))/theta^2 - 1/theta^2
        sage: integrate(t*cos(-theta*t),(t,0,pi))
        (pi*theta*sin(pi*theta) + cos(pi*theta))/theta^2 - 1/theta^2
        sage: integrate(t*cos(-theta*t),t)
        (t*theta*sin(t*theta) + cos(t*theta))/theta^2
        sage: integrate(x^2,(x)) # this worked before
        1/3*x^3
        sage: integrate(x^2,(x,)) # this didn't
        1/3*x^3
        sage: integrate(x^2,(x,1,2))
        7/3
        sage: integrate(x^2,(x,1,2,3))
        Traceback (most recent call last):
        ...
        TypeError: invalid input (x, 1, 2, 3) - please use variable, with or without two endpoints

    Note that this used to be the test, but it is actually divergent
    (Maxima currently asks for assumptions on theta)::

        sage: integrate(t*cos(-theta*t),(t,-oo,oo))
        Traceback (most recent call last):
        ...
        ValueError: Computation failed since Maxima requested additional constraints;...

    Check if :issue:`6189` is fixed::

        sage: n = N; n
        <function numerical_approx at ...>
        sage: F(x) = 1/sqrt(2*pi*1^2)*exp(-1/(2*1^2)*(x-0)^2)
        sage: G(x) = 1/sqrt(2*pi*n(1)^2)*exp(-1/(2*n(1)^2)*(x-n(0))^2)
        sage: integrate( (F(x)-F(x))^2, x, -infinity, infinity).n()
        0.000000000000000
        sage: integrate( ((F(x)-G(x))^2).expand(), x, -infinity, infinity).n()
        -6.26376265908397e-17
        sage: integrate( (F(x)-G(x))^2, x, -infinity, infinity).n()# abstol 1e-6
        0

    This was broken before Maxima 5.20::

        sage: exp(-x*i).integral(x,0,1)
        I*e^(-I) - I

    Test that :issue:`8729` is fixed::

        sage: t = var('t')
        sage: a = sqrt((sin(t))^2 + (cos(t))^2)
        sage: integrate(a, t, 0, 2*pi)
        2*pi
        sage: a.simplify_full().simplify_trig()
        1

    Maxima uses Cauchy Principal Value calculations to
    integrate certain convergent integrals.  Here we test
    that this does not raise an error message (see :issue:`11987`)::

        sage: integrate(sin(x)*sin(x/3)/x^2, x, 0, oo)
        1/6*pi

    Maxima returned a negative value for this integral prior to
    maxima-5.24 (:issue:`10923`). Ideally we would get an answer in terms
    of the gamma function; however, we get something equivalent::

        sage: actual_result = integral(e^(-1/x^2), x, 0, 1)
        sage: actual_result.canonicalize_radical()
        (sqrt(pi)*(erf(1)*e - e) + 1)*e^(-1)
        sage: ideal_result = 1/2*gamma(-1/2, 1)
        sage: error = actual_result - ideal_result
        sage: error.numerical_approx() # abs tol 1e-10
        0

    We will get a correct answer here, which is better than
    the previous (wrong) answer of zero. See :issue:`10914`::

        sage: f = abs(sin(x))
        sage: result = integrate(f, x, 0, 2*pi)
        ...
        sage: result
        4

    Another incorrect integral fixed upstream in Maxima, from
    :issue:`11233`::

        sage: a,t = var('a,t')
        sage: assume(a>0)
        sage: assume(x>0)
        sage: f = log(1 + a/(x * t)^2)
        sage: F = integrate(f, t, 1, Infinity)
        sage: F(x=1, a=7).numerical_approx() # abs tol 1e-10
        4.32025625668262
        sage: forget()

    Verify that MinusInfinity works with sympy (:issue:`12345`)::

        sage: integral(1/x^2, x, -infinity, -1, algorithm='sympy')                      # needs sympy
        1

    Check that :issue:`11737` is fixed::

        sage: N(integrate(sin(x^2)/(x^2), x, 1, infinity), prec=54)
        0.285736646322853
        sage: N(integrate(sin(x^2)/(x^2), x, 1, infinity))  # known bug (nonzero imag part)
        0.285736646322853

    Check that :issue:`14209` is fixed::

        sage: integral(e^(-abs(x))/cosh(x),x,-infinity,infinity)
        2*log(2)
        sage: integral(e^(-abs(x))/cosh(x),x,-infinity,infinity)
        2*log(2)

    Check that :issue:`12628` is fixed::

        sage: var('z,n')
        (z, n)
        sage: f(z, n) = sin(n*z) / (n*z)
        sage: integrate(f(z,1)*f(z,3)*f(z,5)*f(z,7),z,0,oo)
        22/315*pi
        sage: for k in srange(1, 16, 2):
        ....:     print(integrate(prod(f(z, ell)
        ....:                     for ell in srange(1, k+1, 2)), z, 0, oo))
        1/2*pi
        1/6*pi
        1/10*pi
        22/315*pi
        3677/72576*pi
        48481/1247400*pi
        193359161/6227020800*pi
        5799919/227026800*pi

    Check that :issue:`12628` is fixed::

        sage: integrate(1/(sqrt(x)*((1+sqrt(x))^2)),x,1,9)
        1/2

    Check that :issue:`8728` is fixed::

        sage: forget()
        sage: c,w,T = var('c,w,T')
        sage: assume(1-c^2 > 0)
        sage: assume(abs(c) - sqrt(1-c^2) - 1 > 0)
        sage: assume(abs(sqrt(1-c^2)-1) - abs(c) > 0)
        sage: integrate(cos(w+T) / (1+c*cos(T))^2, T, 0, 2*pi)
        2*pi*sqrt(-c^2 + 1)*c*cos(w)/(c^4 - 2*c^2 + 1)

    Check that :issue:`13733` is fixed (but the bug reappeared, see :issue:`30063`)::

        sage: a = integral(log(cot(x) - 1), x, 0, pi/4); a  # long time (about 6 s) # known bug
        -1/4*pi*log(2) - 1/2*I*dilog(I + 1) + 1/2*I*dilog(-I + 1) + 1/2*I*dilog(1/2*I + 1/2) - 1/2*I*dilog(-1/2*I + 1/2)
        sage: abs(N(a - pi*log(2)/8)) < 1e-15  # long time # known bug
        True

    Check that :issue:`17968` is fixed::

        sage: a = N(integrate(exp(x^3), (x, 1, 2)), prec=54)
        sage: a.real_part()    # abs tol 1e-13
        275.510983763312
        sage: a.imag_part()    # abs tol 1e-13
        0.0

    This used to be solved by the ``abs_integrate`` Maxima package
    but can be solved now without it::

        sage: integrate(abs(x), x)
        1/2*x*abs(x)
        sage: integral(abs(cos(x))*sin(x),(x,pi/2,pi))
        1/2
        sage: f = (x^2)*exp(x) / (1+exp(x))^2
        sage: integrate(f, (x, -infinity, infinity))
        1/3*pi^2

    Some integrals are now working (:issue:`27958`, using giac or sympy)::

        sage: result = integrate(abs(x^2 - 1), x, -2, 2)
        ...
        sage: result
        4

        sage: result = integrate(sgn(x) - sgn(1-x), x)
        ...
        sage: result
        abs(x - 1) + abs(x)

        sage: result = integrate(1 / (1 + abs(x-5)), x, -5, 6)
        ...
        sage: result
        log(11) + log(2)

        sage: result = integrate(abs(x^2 - 1), x, -2, 2)
        ...
        sage: result
        4

    Examples that only giac can correctly integrate (for now)::

        sage: # needs sage.libs.giac
        sage: g = abs(sin(x)*cos(x))
        sage: result = g.integrate(x, 0, 2*pi)
        ...
        sage: result
        2

    ::

        sage: # needs sage.libs.giac
        sage: f = sqrt(x + 1/x^2)
        sage: actual = integrate(f, x)
        ...
        sage: expected = (1/3*(2*sqrt(x^3 + 1) - log(sqrt(x^3 + 1) + 1)
        ....:             + log(abs(sqrt(x^3 + 1) - 1)))*sgn(x))
        sage: bool(actual == expected)
        True

    ::

        sage: # needs sage.libs.giac
        sage: result = integrate(cos(x + abs(x)), x)
        ...
        sage: result
        sin(x*sgn(x) + x)/(sgn(x) + 1)

    ::

        sage: # needs sage.libs.giac
        sage: result = integrate(1/(1 + abs(x)), x)
        ...
        sage: result
        log(abs(x*sgn(x) + 1))/sgn(x)

    ::

        sage: # needs sage.libs.giac
        sage: result = integrate(1/sqrt(abs(x)), x)
        ...
        sage: result
        2*sqrt(x*sgn(x))/sgn(x)

    ::

        sage: # needs sage.libs.giac
        sage: result = integrate(1/(1 + abs(x)), x)
        ...
        sage: result
        log(abs(x*sgn(x) + 1))/sgn(x)

    ::

        sage: # needs sage.libs.giac
        sage: result = integrate(cos(x + abs(x)), x)
        ...
        sage: result
        sin(x*sgn(x) + x)/(sgn(x) + 1)

    Some tests for :issue:`17468`::

        sage: integral(log(abs(2*sin(x))), x, 0, pi/3)
        1/36*I*pi^2 + I*dilog(1/2*I*sqrt(3) + 1/2) + I*dilog(-1/2*I*sqrt(3) - 1/2)
        sage: integral(log(abs(sin(x))), x, 0, pi/2)
        -1/2*pi*log(2)

    Check that :issue:`25823` is fixed::

        sage: f = log(sin(x))*sin(x)^2
        sage: g = integrate(f, x) ; g
        1/4*I*x^2
        - 1/2*I*x*arctan2(sin(x), cos(x) + 1)
        + 1/2*I*x*arctan2(sin(x), -cos(x) + 1)
        - 1/4*x*log(cos(x)^2 + sin(x)^2 + 2*cos(x) + 1)
        - 1/4*x*log(cos(x)^2 + sin(x)^2 - 2*cos(x) + 1)
        + 1/4*(2*x - sin(2*x))*log(sin(x))
        + 1/4*x
        + 1/2*I*dilog(-e^(I*x))
        + 1/2*I*dilog(e^(I*x)) + 1/8*sin(2*x)

    Indeed::

        sage: (g.derivative() - f).full_simplify().full_simplify()
        0

    Test for :issue:`24117`::

        sage: integrate(sqrt(1-4*sin(x)^2),x, algorithm='maxima')
        integrate(sqrt(-4*sin(x)^2 + 1), x)

    Check that :issue:`30353` is fixed::

        sage: a = SR.var('a')
        sage: assume(a > 0)
        sage: assume(a < 1)
        sage: integrate(x*log(1/(a*x+(1-x)^2)), x, 0, 1, algorithm='maxima')
        1/4*a^2*log(a) + 1/2*sqrt(-a^2 + 4*a)*a*arctan(sqrt(-a^2 + 4*a)*(a - 2)/(a^2 - 4*a)) - 1/2*sqrt(-a^2 + 4*a)*a*arctan(sqrt(-a^2 + 4*a)/(a - 4)) - a*log(a) - sqrt(-a^2 + 4*a)*arctan(sqrt(-a^2 + 4*a)*(a - 2)/(a^2 - 4*a)) + sqrt(-a^2 + 4*a)*arctan(sqrt(-a^2 + 4*a)/(a - 4)) - 1/2*a + 3/2

    Check that :issue:`25905` is fixed::

        sage: var('a d x c')
        (a, d, x, c)
        sage: f = (I*a*tan(d*x + c) + a)^3*tan(d*x + c)
        sage: integrate(f, x, algorithm='fricas')  # optional - fricas
        -2/3*(24*a^3*e^(4*I*d*x + 4*I*c) + 33*a^3*e^(2*I*d*x + 2*I*c) + 13*a^3 + 6*(a^3*e^(6*I*d*x + 6*I*c) + 3*a^3*e^(4*I*d*x + 4*I*c) + 3*a^3*e^(2*I*d*x + 2*I*c) + a^3)*log(e^(2*I*d*x + 2*I*c) + 1))/(d*e^(6*I*d*x + 6*I*c) + 3*d*e^(4*I*d*x + 4*I*c) + 3*d*e^(2*I*d*x + 2*I*c) + d)

    The fundamental theorem of calculus holds for elliptic integrals
    of the second kind, :issue:`26563`::

        sage: x,m = SR.var('x,m', domain='real')    # long time
        sage: integrate(elliptic_e(x,m).diff(x), x) # long time
        elliptic_e(x, m)

    Check that :issue:`20467` is fixed::

        sage: k = var('k')
        sage: integral(sin(k*x)/x*erf(x^2), x, 0, oo, algorithm='maxima')
        integrate(erf(x^2)*sin(k*x)/x, x, 0, +Infinity)
    """
    expression, v, a, b = _normalize_integral_input(expression, v, a, b)
    if algorithm is not None:
        integrator = available_integrators.get(algorithm)
        if not integrator:
            raise ValueError("Unknown algorithm: %s" % algorithm)
        return integrator(expression, v, a, b)
    if a is None:
        return indefinite_integral(expression, v, hold=hold)
    else:
        return definite_integral(expression, v, a, b, hold=hold)


integral = integrate
