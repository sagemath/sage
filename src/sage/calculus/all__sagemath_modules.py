# sage_setup: distribution = sagemath-modules
from sage.calculus.all__sagemath_categories import *

from sage.calculus.integration import numerical_integral, monte_carlo_integral
integral_numerical = numerical_integral

from sage.calculus.interpolation import spline, Spline

from sage.calculus.functions import wronskian, jacobian

from sage.calculus.ode import ode_solver, ode_system

# We lazy_import the following modules since they import numpy which slows down sage startup
from sage.misc.lazy_import import lazy_import
lazy_import("sage.calculus.riemann", ["Riemann_Map"])
lazy_import("sage.calculus.interpolators", ["polygon_spline", "complex_cubic_spline"])

from sage.calculus.transforms.all import *
del lazy_import
