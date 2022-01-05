
from .calculus import maxima as maxima_calculus
from .calculus import (laplace, inverse_laplace,
                      limit, lim)

from .integration import numerical_integral, monte_carlo_integral
integral_numerical = numerical_integral

from .interpolation import spline, Spline

from .functional import (diff, derivative,
                        expand,
                        taylor, simplify)

from .functions import (wronskian,jacobian)

from .ode import ode_solver, ode_system

from .desolvers import (desolve, desolve_laplace, desolve_system,
                       eulers_method, eulers_method_2x2,
                       eulers_method_2x2_plot, desolve_rk4, desolve_system_rk4,
                       desolve_odeint, desolve_mintides, desolve_tides_mpfr)

from .var import (var, function, clear_vars)

from .transforms.all import *

# We lazy_import the following modules since they import numpy which slows down sage startup
from sage.misc.lazy_import import lazy_import
lazy_import("sage.calculus.riemann",["Riemann_Map"])
lazy_import("sage.calculus.interpolators",["polygon_spline","complex_cubic_spline"])

from sage.modules.all import vector

def symbolic_expression(x):
    """
    Create a symbolic expression or vector of symbolic expressions from x.

    INPUT:

    - ``x`` - an object

    OUTPUT:

    - a symbolic expression.

    EXAMPLES::

        sage: a = symbolic_expression(3/2); a
        3/2
        sage: type(a)
        <type 'sage.symbolic.expression.Expression'>
        sage: R.<x> = QQ[]; type(x)
        <type 'sage.rings.polynomial.polynomial_rational_flint.Polynomial_rational_flint'>
        sage: a = symbolic_expression(2*x^2 + 3); a
        2*x^2 + 3
        sage: type(a)
        <type 'sage.symbolic.expression.Expression'>
        sage: from sage.symbolic.expression import is_Expression
        sage: is_Expression(a)
        True
        sage: a in SR
        True
        sage: a.parent()
        Symbolic Ring

    Note that equations exist in the symbolic ring::

        sage: E = EllipticCurve('15a'); E
        Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 10*x - 10 over Rational Field
        sage: symbolic_expression(E)
        x*y + y^2 + y == x^3 + x^2 - 10*x - 10
        sage: symbolic_expression(E) in SR
        True

    If ``x`` is a list or tuple, create a vector of symbolic expressions::

        sage: v=symbolic_expression([x,1]); v
        (x, 1)
        sage: v.base_ring()
        Symbolic Ring
        sage: v=symbolic_expression((x,1)); v
        (x, 1)
        sage: v.base_ring()
        Symbolic Ring
        sage: v=symbolic_expression((3,1)); v
        (3, 1)
        sage: v.base_ring()
        Symbolic Ring
        sage: E = EllipticCurve('15a'); E
        Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 10*x - 10 over Rational Field
        sage: v=symbolic_expression([E,E]); v
        (x*y + y^2 + y == x^3 + x^2 - 10*x - 10, x*y + y^2 + y == x^3 + x^2 - 10*x - 10)
        sage: v.base_ring()
        Symbolic Ring

    If ``x`` is a function, for example defined by a ``lambda`` expression, create a
    symbolic function::

        sage: f = symbolic_expression(lambda z: z^2 + 1); f
        z |--> z^2 + 1
        sage: f.parent()
        Callable function ring with argument z
        sage: f(7)
        50

    If ``x`` is a list or tuple of functions, or if ``x`` is a function that returns a list
    or tuple, create a callable symbolic vector::

        sage: symbolic_expression([lambda mu, nu: mu^2 + nu^2, lambda mu, nu: mu^2 - nu^2])
        (mu, nu) |--> (mu^2 + nu^2, mu^2 - nu^2)
        sage: f = symbolic_expression(lambda uwu: [1, uwu, uwu^2]); f
        uwu |--> (1, uwu, uwu^2)
        sage: f.parent()
        Vector space of dimension 3 over Callable function ring with argument uwu
        sage: f(5)
        (1, 5, 25)
        sage: f(5).parent()
        Vector space of dimension 3 over Symbolic Ring

    TESTS:

    Also functions defined using ``def`` can be used, but we do not advertise it as a use case::

        sage: def sos(x, y):
        ....:     return x^2 + y^2
        sage: symbolic_expression(sos)
        (x, y) |--> x^2 + y^2

    Functions that take a varying number of arguments or keyword-only arguments are not accepted::

        sage: def variadic(x, *y):
        ....:     return x
        sage: symbolic_expression(variadic)
        Traceback (most recent call last):
        ...
        TypeError: unable to convert <function variadic at 0x...> to a symbolic expression

        sage: def function_with_keyword_only_arg(x, *, sign=1):
        ....:     return sign * x
        sage: symbolic_expression(function_with_keyword_only_arg)
        Traceback (most recent call last):
        ...
        TypeError: unable to convert <function function_with_keyword_only_arg at 0x...>
        to a symbolic expression
    """
    from sage.symbolic.expression import Expression
    from sage.symbolic.ring import SR
    if isinstance(x, Expression):
        return x
    elif hasattr(x, '_symbolic_'):
        return x._symbolic_(SR)
    elif isinstance(x, (tuple, list)):
        return vector([symbolic_expression(item) for item in x])
    elif callable(x):
        from inspect import signature, Parameter
        try:
            s = signature(x)
        except ValueError:
            pass
        else:
            if all(param.kind in (Parameter.POSITIONAL_ONLY, Parameter.POSITIONAL_OR_KEYWORD)
                   for param in s.parameters.values()):
                vars = [SR.var(name) for name in s.parameters.keys()]
                result = x(*vars)
                if isinstance(result, (tuple, list)):
                    return vector(SR, result).function(*vars)
                else:
                    return SR(result).function(*vars)
    return SR(x)

from . import desolvers
