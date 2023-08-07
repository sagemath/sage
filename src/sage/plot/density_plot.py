# sage.doctest: needs sage.symbolic
"""
Density plots
"""

#*****************************************************************************
#       Copyright (C) 2006 Alex Clemesha <clemesha@gmail.com>,
#                          William Stein <wstein@gmail.com>,
#                     2008 Mike Hansen <mhansen@gmail.com>,
#                          Arnaud Bergeron <abergeron@gmail.com>
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
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.plot.primitive import GraphicPrimitive
from sage.misc.decorators import options
from sage.plot.colors import get_cmap
from sage.arith.srange import xsrange


class DensityPlot(GraphicPrimitive):
    """
    Primitive class for the density plot graphics type.  See
    ``density_plot?`` for help actually doing density plots.

    INPUT:

    - ``xy_data_array`` - list of lists giving evaluated values of the
      function on the grid

    - ``xrange`` - tuple of 2 floats indicating range for horizontal direction

    - ``yrange`` - tuple of 2 floats indicating range for vertical direction

    - ``options`` - dict of valid plot options to pass to constructor

    EXAMPLES:

    Note this should normally be used indirectly via ``density_plot``::

        sage: from sage.plot.density_plot import DensityPlot
        sage: D = DensityPlot([[1,3],[2,4]], (1,2), (2,3),options={})
        sage: D
        DensityPlot defined by a 2 x 2 data grid
        sage: D.yrange
        (2, 3)
        sage: D.options()
        {}

    TESTS:

    We test creating a density plot::

        sage: x,y = var('x,y')
        sage: density_plot(x^2 - y^3 + 10*sin(x*y), (x,-4,4), (y,-4,4), plot_points=121, cmap='hsv')
        Graphics object consisting of 1 graphics primitive
    """
    def __init__(self, xy_data_array, xrange, yrange, options):
        """
        Initializes base class DensityPlot.

        EXAMPLES::

            sage: x,y = var('x,y')
            sage: D = density_plot(x^2 - y^3 + 10*sin(x*y), (x,-4,4), (y,-4,4), plot_points=121, cmap='hsv')
            sage: D[0].xrange
            (-4.0, 4.0)
            sage: D[0].options()['plot_points']
            121
        """
        self.xrange = xrange
        self.yrange = yrange
        self.xy_data_array = xy_data_array
        self.xy_array_row = len(xy_data_array)
        self.xy_array_col = len(xy_data_array[0])
        GraphicPrimitive.__init__(self, options)

    def get_minmax_data(self):
        """
        Returns a dictionary with the bounding box data.

        EXAMPLES::

            sage: x,y = var('x,y')
            sage: f(x, y) = x^2 + y^2
            sage: d = density_plot(f, (3,6), (3,6))[0].get_minmax_data()
            sage: d['xmin']
            3.0
            sage: d['ymin']
            3.0
        """
        from sage.plot.plot import minmax_data
        return minmax_data(self.xrange, self.yrange, dict=True)

    def _allowed_options(self):
        """
        Return the allowed options for the DensityPlot class.

        TESTS::

            sage: isinstance(density_plot(x, (-2,3), (1,10))[0]._allowed_options(), dict)
            True
        """
        return {'plot_points': 'How many points to use for plotting precision',
                'cmap': """the name of a predefined colormap,
                       a list of colors or an instance of a
                       matplotlib Colormap. Type: import matplotlib.cm; matplotlib.cm.datad.keys()
                       for available colormap names.""",
                'interpolation': 'What interpolation method to use'}

    def _repr_(self):
        """
        String representation of DensityrPlot primitive.

        EXAMPLES::

            sage: x,y = var('x,y')
            sage: D = density_plot(x^2 - y^2, (x,-2,2), (y,-2,2))
            sage: d = D[0]; d
            DensityPlot defined by a 25 x 25 data grid
        """
        return "DensityPlot defined by a %s x %s data grid"%(self.xy_array_row, self.xy_array_col)

    def _render_on_subplot(self, subplot):
        """
        TESTS:

        A somewhat random plot, but fun to look at::

            sage: x,y = var('x,y')
            sage: density_plot(x^2 - y^3 + 10*sin(x*y), (x,-4,4), (y,-4,4), plot_points=121, cmap='hsv')
            Graphics object consisting of 1 graphics primitive
        """
        options = self.options()
        cmap = get_cmap(options['cmap'])

        x0, x1 = float(self.xrange[0]), float(self.xrange[1])
        y0, y1 = float(self.yrange[0]), float(self.yrange[1])

        subplot.imshow(self.xy_data_array, origin='lower',
                       cmap=cmap, extent=(x0,x1,y0,y1),
                       interpolation=options['interpolation'])


@options(plot_points=25, cmap='gray', interpolation='catrom')
def density_plot(f, xrange, yrange, **options):
    r"""
    ``density_plot`` takes a function of two variables, `f(x,y)`
    and plots the height of the function over the specified
    ``xrange`` and ``yrange`` as demonstrated below.

    ``density_plot(f, (xmin,xmax), (ymin,ymax), ...)``

    INPUT:

    - ``f`` -- a function of two variables

    - ``(xmin,xmax)`` -- 2-tuple, the range of ``x`` values OR 3-tuple
      ``(x,xmin,xmax)``

    - ``(ymin,ymax)`` -- 2-tuple, the range of ``y`` values OR 3-tuple
      ``(y,ymin,ymax)``

    The following inputs must all be passed in as named parameters:

    - ``plot_points`` -- integer (default: 25); number of points to plot
      in each direction of the grid

    - ``cmap`` -- a colormap (default: ``'gray'``), the name of
      a predefined colormap, a list of colors or an instance of a matplotlib
      Colormap. Type: ``import matplotlib.cm; matplotlib.cm.datad.keys()``
      for available colormap names.

    - ``interpolation`` -- string (default: ``'catrom'``), the interpolation
      method to use: ``'bilinear'``, ``'bicubic'``, ``'spline16'``,
      ``'spline36'``, ``'quadric'``, ``'gaussian'``, ``'sinc'``,
      ``'bessel'``, ``'mitchell'``, ``'lanczos'``, ``'catrom'``,
      ``'hermite'``, ``'hanning'``, ``'hamming'``, ``'kaiser'``


    EXAMPLES:

    Here we plot a simple function of two variables.  Note that
    since the input function is an expression, we need to explicitly
    declare the variables in 3-tuples for the range::

        sage: x,y = var('x,y')
        sage: density_plot(sin(x) * sin(y), (x,-2,2), (y,-2,2))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        g = density_plot(sin(x) * sin(y), (x,-2,2), (y,-2,2))
        sphinx_plot(g)

    Here we change the ranges and add some options; note that here
    ``f`` is callable (has variables declared), so we can use 2-tuple ranges::

        sage: x,y = var('x,y')
        sage: f(x,y) = x^2 * cos(x*y)
        sage: density_plot(f, (x,-10,5), (y,-5,5), interpolation='sinc', plot_points=100)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        def f(x,y): return x**2 * cos(x*y)
        g = density_plot(f, (x,-10,5), (y,-5,5), interpolation='sinc', plot_points=100)
        sphinx_plot(g)

    An even more complicated plot::

        sage: x,y = var('x,y')
        sage: density_plot(sin(x^2+y^2) * cos(x) * sin(y), (x,-4,4), (y,-4,4), cmap='jet', plot_points=100)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        g = density_plot(sin(x**2 + y**2)*cos(x)*sin(y), (x,-4,4), (y,-4,4), cmap='jet', plot_points=100)
        sphinx_plot(g)

    This should show a "spotlight" right on the origin::

        sage: x,y = var('x,y')
        sage: density_plot(1/(x^10 + y^10), (x,-10,10), (y,-10,10))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        g = density_plot(1/(x**10 + y**10), (x,-10,10), (y,-10,10))
        sphinx_plot(g)

    Some elliptic curves, but with symbolic endpoints.  In the first
    example, the plot is rotated 90 degrees because we switch the
    variables `x`, `y`::

        sage: density_plot(y^2 + 1 - x^3 - x, (y,-pi,pi), (x,-pi,pi))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        g = density_plot(y**2 + 1 - x**3 - x, (y,-pi,pi), (x,-pi,pi))
        sphinx_plot(g)

    ::

        sage: density_plot(y^2 + 1 - x^3 - x, (x,-pi,pi), (y,-pi,pi))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        g = density_plot(y**2 + 1 - x**3 - x, (x,-pi,pi), (y,-pi,pi))
        sphinx_plot(g)

    Extra options will get passed on to show(), as long as they are valid::

        sage: density_plot(log(x) + log(y), (x,1,10), (y,1,10), aspect_ratio=1)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        g = density_plot(log(x) + log(y), (x,1,10), (y,1,10), aspect_ratio=1)
        sphinx_plot(g)

    ::

        sage: density_plot(log(x) + log(y), (x,1,10), (y,1,10)).show(aspect_ratio=1) # These are equivalent

    TESTS:

    Check that :trac:`15315` is fixed, i.e., density_plot respects the
    ``aspect_ratio`` parameter. Without the fix, it looks like a thin line
    of width a few mm. With the fix it should look like a nice fat layered
    image::

        sage: density_plot((x*y)^(1/2), (x,0,3), (y,0,500), aspect_ratio=.01)
        Graphics object consisting of 1 graphics primitive

    Default ``aspect_ratio`` is ``"automatic"``, and that should work too::

        sage: density_plot((x*y)^(1/2), (x,0,3), (y,0,500))
        Graphics object consisting of 1 graphics primitive

    Check that :trac:`17684` is fixed, i.e., symbolic values can be plotted::

        sage: def f(x,y):
        ....:     return SR(x)
        sage: density_plot(f, (0,1), (0,1))
        Graphics object consisting of 1 graphics primitive
    """
    from sage.plot.all import Graphics
    from sage.plot.misc import setup_for_eval_on_grid
    from sage.rings.real_double import RDF
    g, ranges = setup_for_eval_on_grid([f], [xrange, yrange], options['plot_points'])
    g = g[0]
    xrange, yrange = [r[:2] for r in ranges]

    xy_data_array = [[RDF(g(x,y)) for x in xsrange(*ranges[0], include_endpoint=True)]
                            for y in xsrange(*ranges[1], include_endpoint=True)]

    g = Graphics()
    g._set_extra_kwds(Graphics._extract_kwds_for_show(options, ignore=['xmin','xmax']))
    g.add_primitive(DensityPlot(xy_data_array, xrange, yrange, options))
    return g
