from .all__sagemath_categories import *

from .integration import numerical_integral, monte_carlo_integral
integral_numerical = numerical_integral

from .interpolation import spline, Spline

from .functions import wronskian, jacobian

from .ode import ode_solver, ode_system

# We lazy_import the following modules since they import numpy which slows down sage startup
from sage.misc.lazy_import import lazy_import
lazy_import("sage.calculus.riemann", ["Riemann_Map"])
lazy_import("sage.calculus.interpolators", ["polygon_spline","complex_cubic_spline"])

from .transforms.all import *
