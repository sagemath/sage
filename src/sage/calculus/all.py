
from .calculus import maxima as maxima_calculus
from .calculus import (laplace, inverse_laplace,
                       limit, lim)

from .integration import numerical_integral, monte_carlo_integral
integral_numerical = numerical_integral

from .interpolation import spline, Spline

from .functional import (diff, derivative,
                         expand,
                         taylor, simplify)

from .functions import (wronskian, jacobian)

from .ode import ode_solver, ode_system

from .desolvers import (desolve, desolve_laplace, desolve_system,
                        eulers_method, eulers_method_2x2,
                        eulers_method_2x2_plot, desolve_rk4, desolve_system_rk4,
                        desolve_odeint, desolve_mintides, desolve_tides_mpfr)

from sage.calculus.expr import symbolic_expression
from sage.calculus.var import (var, function, clear_vars)

from .transforms.all import *

# We lazy_import the following modules since they import numpy which slows down sage startup
from sage.misc.lazy_import import lazy_import
lazy_import("sage.calculus.riemann", ["Riemann_Map"])
lazy_import("sage.calculus.interpolators", ["polygon_spline", "complex_cubic_spline"])

from . import desolvers
