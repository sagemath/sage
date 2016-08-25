from __future__ import absolute_import

from .calculus import maxima as maxima_calculus
from .calculus import (laplace, inverse_laplace,
                      limit, lim)

from .functional import (diff, derivative,
                        expand,
                        taylor, simplify)

from .functions import (wronskian,jacobian)

from .desolvers import (desolve, desolve_laplace, desolve_system,
                       eulers_method, eulers_method_2x2,
                       eulers_method_2x2_plot, desolve_rk4, desolve_system_rk4,
                       desolve_odeint, desolve_mintides, desolve_tides_mpfr)

from .var import (var, function, clear_vars)

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

    If x is a list or tuple, create a vector of symbolic expressions::

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
    """
    from sage.symbolic.expression import Expression
    from sage.symbolic.ring import SR
    if isinstance(x, Expression):
        return x
    elif hasattr(x, '_symbolic_'):
        return x._symbolic_(SR)
    elif isinstance(x, (tuple,list)):
        return vector(SR,x)
    else:
        return SR(x)

from . import desolvers
