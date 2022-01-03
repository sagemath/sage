"""
Complex Plots
"""
# ****************************************************************************
#       Copyright (C) 2009 Robert Bradshaw <robertwb@math.washington.edu>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

# TODO: use NumPy buffers and complex fast_callable (when supported)
from cysignals.signals cimport sig_on, sig_off, sig_check

cimport numpy as cnumpy

from sage.plot.primitive import GraphicPrimitive
from sage.misc.decorators import options
from sage.rings.complex_double cimport ComplexDoubleElement
from sage.arith.srange import srange

from libc.math cimport hypot, atan2, atan, log, sqrt
from sage.arith.constants cimport M_PI as PI

from sage.libs.gsl.complex cimport *

cdef inline ComplexDoubleElement new_CDF_element(double x, double y):
    z = <ComplexDoubleElement>ComplexDoubleElement.__new__(ComplexDoubleElement)
    GSL_SET_COMPLEX(&z._complex, x, y)
    return z


cdef inline double mag_to_lightness(double r):
    """
    Tweak this to adjust how the magnitude affects the color.
    For instance, changing ``sqrt(r)`` to ``r`` will cause
    anything near a zero to be much darker and poles to be
    much lighter, while ``r**(.25)`` would cause the reverse
    effect.

    INPUT:

    - ``r`` - a non-negative real number

    OUTPUT:

    A value between `-1` (black) and `+1` (white), inclusive.

    EXAMPLES:

    This tests it implicitly::

        sage: from sage.plot.complex_plot import complex_to_rgb
        sage: complex_to_rgb([[0, 1, 10]])
        array([[[0.        , 0.        , 0.        ],
                [0.77172568, 0.        , 0.        ],
                [1.        , 0.22134776, 0.22134776]]])
    """
    return atan(log(sqrt(r)+1)) * (4/PI) - 1


def complex_to_rgb(z_values):
    """
    INPUT:

    - ``z_values`` -- A grid of complex numbers, as a list of lists

    OUTPUT:

    An `N \\times M \\times 3` floating point Numpy array ``X``, where
    ``X[i,j]`` is an (r,g,b) tuple.

    EXAMPLES::

        sage: from sage.plot.complex_plot import complex_to_rgb
        sage: complex_to_rgb([[0, 1, 1000]])
        array([[[0.        , 0.        , 0.        ],
                [0.77172568, 0.        , 0.        ],
                [1.        , 0.64421177, 0.64421177]]])
        sage: complex_to_rgb([[0, 1j, 1000j]])
        array([[[0.        , 0.        , 0.        ],
                [0.38586284, 0.77172568, 0.        ],
                [0.82210588, 1.        , 0.64421177]]])
    """
    import numpy
    cdef unsigned int i, j, imax, jmax
    cdef double x, y, mag, arg
    cdef double lightness, hue, top, bot
    cdef double r, g, b
    cdef int ihue
    cdef ComplexDoubleElement z
    from sage.rings.complex_double import CDF

    imax = len(z_values)
    jmax = len(z_values[0])
    cdef cnumpy.ndarray[cnumpy.float_t, ndim=3, mode='c'] rgb = numpy.empty(dtype=float, shape=(imax, jmax, 3))

    sig_on()
    for i from 0 <= i < imax:

        row = z_values[i]
        for j from 0 <= j < jmax:

            zz = row[j]
            if type(zz) is ComplexDoubleElement:
                z = <ComplexDoubleElement>zz
            else:
                z = CDF(zz)
            x = GSL_REAL(z._complex)
            y = GSL_IMAG(z._complex)
            mag = hypot(x, y)
            arg = atan2(y, x) # math module arctan has range from -pi to pi, so cut along negative x-axis

            lightness = mag_to_lightness(mag)
            if lightness < 0: # in hsv, variable value, full saturation (s=1, v=1+lightness)
                bot = 0
                top = (1+lightness)
            else: # in hsv, variable saturation, full value (v=1, s=1-lightness)
                bot = lightness
                top = 1

            hue = 3*arg/PI # Note that does same thing as colorsys module hsv_to_rgb for this setup, but in Cython
            if hue < 0: hue += 6 # usual hsv hue is thus h=arg/(2*pi) for positive, h=arg/(2*PI)+1 for negative
            ihue = <int>hue
            if ihue == 0:
                r = top
                g = bot + hue * (top-bot)
                b = bot
            elif ihue == 1:
                r = bot + (2-hue) * (top-bot)
                g = top
                b = bot
            elif ihue == 2:
                r = bot
                g = top
                b = bot + (hue-2) * (top-bot)
            elif ihue == 3:
                r = bot
                g = bot + (4-hue) * (top-bot)
                b = top
            elif ihue == 4:
                r = bot + (hue-4) * (top-bot)
                g = bot
                b = top
            else:
                r = top
                g = bot
                b = bot + (6-hue) * (top-bot)

            rgb[i, j, 0] = r
            rgb[i, j, 1] = g
            rgb[i, j, 2] = b

    sig_off()
    return rgb


class ComplexPlot(GraphicPrimitive):
    """
    The GraphicsPrimitive to display complex functions in using the domain
    coloring method

    INPUT:

    - ``rgb_data`` -- An array of colored points to be plotted.

    - ``x_range`` -- A minimum and maximum x value for the plot.

    - ``y_range`` -- A minimum and maximum y value for the plot.

    TESTS::

        sage: p = complex_plot(lambda z: z^2-1, (-2, 2), (-2, 2))
    """
    def __init__(self, rgb_data, x_range, y_range, options):
        """
        TESTS::

            sage: p = complex_plot(lambda z: z^2-1, (-2, 2), (-2, 2))
        """
        self.x_range = x_range
        self.y_range = y_range
        self.x_count = len(rgb_data)
        self.y_count = len(rgb_data[0])
        self.rgb_data = rgb_data
        GraphicPrimitive.__init__(self, options)

    def get_minmax_data(self):
        """
        Return a dictionary with the bounding box data.

        EXAMPLES::

            sage: p = complex_plot(lambda z: z, (-1, 2), (-3, 4))
            sage: sorted(p.get_minmax_data().items())
            [('xmax', 2.0), ('xmin', -1.0), ('ymax', 4.0), ('ymin', -3.0)]
            sage: p = complex_plot(lambda z: z, (1, 2), (3, 4))
            sage: sorted(p.get_minmax_data().items())
            [('xmax', 2.0), ('xmin', 1.0), ('ymax', 4.0), ('ymin', 3.0)]
        """
        from sage.plot.plot import minmax_data
        return minmax_data(self.x_range, self.y_range, dict=True)

    def _allowed_options(self):
        """
        TESTS::

            sage: isinstance(complex_plot(lambda z: z, (-1,1), (-1,1))[0]._allowed_options(), dict)
            True
        """
        return {'plot_points':'How many points to use for plotting precision',
                'interpolation':'What interpolation method to use'}

    def _repr_(self):
        """
        TESTS::

            sage: isinstance(complex_plot(lambda z: z, (-1,1), (-1,1))[0]._repr_(), str)
            True
        """
        return "ComplexPlot defined by a %s x %s data grid"%(self.x_count, self.y_count)

    def _render_on_subplot(self, subplot):
        """
        TESTS::

            sage: complex_plot(lambda x: x^2, (-5, 5), (-5, 5))
            Graphics object consisting of 1 graphics primitive
        """
        options = self.options()
        x0, x1 = float(self.x_range[0]), float(self.x_range[1])
        y0, y1 = float(self.y_range[0]), float(self.y_range[1])
        subplot.imshow(self.rgb_data, origin='lower', extent=(x0, x1, y0, y1),
                       interpolation=options['interpolation'])


@options(plot_points=100, interpolation='catrom')
def complex_plot(f, x_range, y_range, **options):
    r"""
    ``complex_plot`` takes a complex function of one variable,
    `f(z)` and plots output of the function over the specified
    ``x_range`` and ``y_range`` as demonstrated below. The magnitude of the
    output is indicated by the brightness (with zero being black and
    infinity being white) while the argument is represented by the
    hue (with red being positive real, and increasing through orange,
    yellow, ... as the argument increases).

    ``complex_plot(f, (xmin, xmax), (ymin, ymax), ...)``

    INPUT:

    - ``f`` -- a function of a single complex value `x + iy`

    - ``(xmin, xmax)`` -- 2-tuple, the range of ``x`` values

    - ``(ymin, ymax)`` -- 2-tuple, the range of ``y`` values

    The following inputs must all be passed in as named parameters:

    - ``plot_points`` -- integer (default: 100); number of points to plot
      in each direction of the grid

    - ``interpolation`` -- string (default: ``'catrom'``), the interpolation
      method to use: ``'bilinear'``, ``'bicubic'``, ``'spline16'``,
      ``'spline36'``, ``'quadric'``, ``'gaussian'``, ``'sinc'``,
      ``'bessel'``, ``'mitchell'``, ``'lanczos'``, ``'catrom'``,
      ``'hermite'``, ``'hanning'``, ``'hamming'``, ``'kaiser'``


    EXAMPLES:

    Here we plot a couple of simple functions::

        sage: complex_plot(sqrt(x), (-5, 5), (-5, 5))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(complex_plot(sqrt(x), (-5, 5), (-5, 5)))

    ::

        sage: complex_plot(sin(x), (-5, 5), (-5, 5))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(complex_plot(sin(x), (-5, 5), (-5, 5)))

    ::

        sage: complex_plot(log(x), (-10, 10), (-10, 10))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(complex_plot(log(x), (-10, 10), (-10, 10)))

    ::

        sage: complex_plot(exp(x), (-10, 10), (-10, 10))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(complex_plot(exp(x), (-10, 10), (-10, 10)))

    A function with some nice zeros and a pole::

        sage: f(z) = z^5 + z - 1 + 1/z
        sage: complex_plot(f, (-3, 3), (-3, 3))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        def f(z): return z**5 + z - 1 + 1/z
        sphinx_plot(complex_plot(f, (-3, 3), (-3, 3)))

    Here is the identity, useful for seeing what values map to what colors::

        sage: complex_plot(lambda z: z, (-3, 3), (-3, 3))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(complex_plot(lambda z: z, (-3, 3), (-3, 3)))

    The Riemann Zeta function::

        sage: complex_plot(zeta, (-30,30), (-30,30))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(complex_plot(zeta, (-30,30), (-30,30)))

    Extra options will get passed on to show(), as long as they are valid::

        sage: complex_plot(lambda z: z, (-3, 3), (-3, 3), figsize=[1,1])
        Graphics object consisting of 1 graphics primitive

    ::

        sage: complex_plot(lambda z: z, (-3, 3), (-3, 3)).show(figsize=[1,1]) # These are equivalent

    TESTS:

    Test to make sure that using fast_callable functions works::

        sage: f(x) = x^2
        sage: g = fast_callable(f, domain=CC, vars='x')
        sage: h = fast_callable(f, domain=CDF, vars='x')
        sage: P = complex_plot(f, (-10, 10), (-10, 10))
        sage: Q = complex_plot(g, (-10, 10), (-10, 10))
        sage: R = complex_plot(h, (-10, 10), (-10, 10))
        sage: S = complex_plot(exp(x)-sin(x), (-10, 10), (-10, 10))
        sage: P; Q; R; S
        Graphics object consisting of 1 graphics primitive
        Graphics object consisting of 1 graphics primitive
        Graphics object consisting of 1 graphics primitive
        Graphics object consisting of 1 graphics primitive

    Test to make sure symbolic functions still work without declaring
    a variable.  (We don't do this in practice because it doesn't use
    fast_callable, so it is much slower.)

    ::

        sage: complex_plot(sqrt, (-5, 5), (-5, 5))
        Graphics object consisting of 1 graphics primitive
    """
    from sage.plot.all import Graphics
    from sage.plot.misc import setup_for_eval_on_grid
    from sage.ext.fast_callable import fast_callable
    from sage.rings.complex_double import CDF

    try:
        f = fast_callable(f, domain=CDF, expect_one_var=True)
    except (AttributeError, TypeError, ValueError):
        pass

    cdef double x, y
    _, ranges = setup_for_eval_on_grid([], [x_range, y_range],
                                       options['plot_points'])
    x_range = ranges[0]
    y_range = ranges[1]
    cdef list z_values = []
    cdef list row
    for y in srange(*y_range, include_endpoint=True):
        row = []
        for x in srange(*x_range, include_endpoint=True):
            sig_check()
            row.append(f(new_CDF_element(x, y)))
        z_values.append(row)

    g = Graphics()
    g._set_extra_kwds(Graphics._extract_kwds_for_show(options, ignore=['xmin', 'xmax']))
    g.add_primitive(ComplexPlot(complex_to_rgb(z_values),
                                x_range[:2], y_range[:2], options))
    return g
