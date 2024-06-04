from sage.plot.graphics import Graphics
from sage.plot.plot import (plot, graphics_array, multi_graphics, list_plot,
                            parametric_plot, polar_plot, plot_loglog, plot_semilogx,
                            plot_semilogy, list_plot_loglog, list_plot_semilogx,
                            list_plot_semilogy)
from sage.plot.line import line, line2d
from sage.plot.arrow import arrow, arrow2d
from sage.plot.bar_chart import bar_chart
from sage.plot.histogram import histogram
from sage.plot.bezier_path import bezier_path
from sage.plot.scatter_plot import scatter_plot
from sage.plot.disk import disk
from sage.plot.point import point, points, point2d
from sage.plot.matrix_plot import matrix_plot
from sage.plot.plot_field import plot_vector_field, plot_slope_field
from sage.plot.text import text
from sage.plot.polygon import polygon, polygon2d
from sage.plot.circle import circle
from sage.plot.ellipse import ellipse
from sage.plot.contour_plot import contour_plot, implicit_plot, region_plot
from sage.plot.density_plot import density_plot
from sage.plot.streamline_plot import streamline_plot

from sage.misc.lazy_import import lazy_import
lazy_import("sage.plot.complex_plot", ["complex_plot"])

from sage.plot.arc import arc

from sage.plot.animate import animate

from sage.plot.plot3d.tachyon import Tachyon

from sage.plot.colors import Color, hue, rainbow, colors, colormaps

from sage.plot.step import plot_step_function

from sage.plot.hyperbolic_arc import hyperbolic_arc
from sage.plot.hyperbolic_polygon import hyperbolic_triangle, hyperbolic_polygon
lazy_import("sage.plot.hyperbolic_regular_polygon", "hyperbolic_regular_polygon")
