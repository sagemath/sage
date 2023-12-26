# sage_setup: distribution = sagemath-plot
from .graphics import Graphics
from .plot import (plot, graphics_array, multi_graphics, list_plot,
                   parametric_plot, polar_plot, plot_loglog, plot_semilogx,
                   plot_semilogy, list_plot_loglog, list_plot_semilogx,
                   list_plot_semilogy)
from .line import line, line2d
from .arrow import arrow, arrow2d
from .bar_chart import bar_chart
from .histogram import histogram
from .bezier_path import bezier_path
from .scatter_plot import scatter_plot
from .disk import disk
from .point import point, points, point2d
from .matrix_plot import matrix_plot
from .plot_field import plot_vector_field, plot_slope_field
from .text import text
from .polygon import polygon, polygon2d
from .circle import circle
from .ellipse import ellipse
from .contour_plot import contour_plot, implicit_plot, region_plot
from .density_plot import density_plot
from .streamline_plot import streamline_plot

from sage.misc.lazy_import import lazy_import
lazy_import("sage.plot.complex_plot", ["complex_plot"])

from sage.plot.arc import arc

from sage.plot.animate import animate

from sage.plot.plot3d.tachyon import Tachyon

from sage.plot.colors import Color, hue, rainbow, colors, colormaps

from sage.plot.step import plot_step_function

lazy_import("sage.plot.hyperbolic_arc", "hyperbolic_arc")
lazy_import("sage.plot.hyperbolic_polygon", [
            "hyperbolic_triangle", "hyperbolic_polygon"])
lazy_import("sage.plot.hyperbolic_regular_polygon", "hyperbolic_regular_polygon")
del lazy_import
