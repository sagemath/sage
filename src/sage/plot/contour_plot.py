# sage.doctest: needs sage.symbolic
"""
Contour plots
"""
# ****************************************************************************
#       Copyright (C) 2006 Alex Clemesha <clemesha@gmail.com>,
#                          William Stein <wstein@gmail.com>,
#                     2008 Mike Hansen <mhansen@gmail.com>,
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
import operator
from sage.plot.primitive import GraphicPrimitive
from sage.misc.decorators import options, suboptions
from sage.plot.colors import rgbcolor, get_cmap
from sage.arith.srange import xsrange


class ContourPlot(GraphicPrimitive):
    """
    Primitive class for the contour plot graphics type.

    See ``contour_plot?`` for help actually doing contour plots.

    INPUT:

    - ``xy_data_array`` -- list of lists giving evaluated values of the function
      on the grid

    - ``xrange`` -- tuple of 2 floats indicating range for horizontal direction

    - ``yrange`` -- tuple of 2 floats indicating range for vertical direction

    - ``options`` -- dictionary of valid plot options to pass to constructor

    EXAMPLES:

    Note this should normally be used indirectly via ``contour_plot``::

        sage: from sage.plot.contour_plot import ContourPlot
        sage: C = ContourPlot([[1,3],[2,4]], (1,2), (2,3), options={})
        sage: C
        ContourPlot defined by a 2 x 2 data grid
        sage: C.xrange
        (1, 2)

    TESTS:

    We test creating a contour plot::

        sage: x,y = var('x,y')
        sage: contour_plot(x^2-y^3+10*sin(x*y), (x,-4,4), (y,-4,4),
        ....:              plot_points=121, cmap='hsv')
        Graphics object consisting of 1 graphics primitive
    """
    def __init__(self, xy_data_array, xrange, yrange, options):
        """
        Initialize base class ``ContourPlot``.

        EXAMPLES::

            sage: x,y = var('x,y')
            sage: C = contour_plot(x^2-y^3+10*sin(x*y), (x,-4,4), (y,-4,4),
            ....:                  plot_points=121, cmap='hsv')
            sage: C[0].xrange
            (-4.0, 4.0)
            sage: C[0].options()['plot_points']
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
        Return a dictionary with the bounding box data.

        EXAMPLES::

            sage: x,y = var('x,y')
            sage: f(x,y) = x^2 + y^2
            sage: d = contour_plot(f, (3,6), (3,6))[0].get_minmax_data()
            sage: d['xmin']
            3.0
            sage: d['ymin']
            3.0
        """
        from sage.plot.plot import minmax_data
        return minmax_data(self.xrange, self.yrange, dict=True)

    def _allowed_options(self):
        """
        Return the allowed options for the ContourPlot class.

        EXAMPLES::

            sage: x,y = var('x,y')
            sage: C = contour_plot(x^2 - y^2, (x,-2,2), (y,-2,2))
            sage: isinstance(C[0]._allowed_options(), dict)
            True
        """
        return {'plot_points': 'How many points to use for plotting precision',
                'cmap': """the name of a predefined colormap,
                        a list of colors, or an instance of a
                        matplotlib Colormap. Type: import matplotlib.cm;
                        matplotlib.cm.datad.keys()
                        for available colormap names.""",
                'colorbar': "Include a colorbar indicating the levels",
                'colorbar_options': "a dictionary of options for colorbars",
                'fill': 'Fill contours or not',
                'legend_label': 'The label for this item in the legend.',
                'contours': """Either an integer specifying the number of
                        contour levels, or a sequence of numbers giving
                        the actual contours to use.""",
                'linewidths': 'the width of the lines to be plotted',
                'linestyles': 'the style of the lines to be plotted',
                'labels': 'show line labels or not',
                'label_options': 'a dictionary of options for the labels',
                'zorder': 'The layer level in which to draw'}

    def _repr_(self):
        """
        String representation of ``ContourPlot`` primitive.

        EXAMPLES::

            sage: x,y = var('x,y')
            sage: C = contour_plot(x^2 - y^2, (x,-2,2), (y,-2,2))
            sage: c = C[0]; c
            ContourPlot defined by a 100 x 100 data grid
        """
        msg = "ContourPlot defined by a %s x %s data grid"
        return msg % (self.xy_array_row, self.xy_array_col)

    def _render_on_subplot(self, subplot):
        """
        TESTS:

        A somewhat random plot, but fun to look at::

            sage: x,y = var('x,y')
            sage: contour_plot(x^2 - y^3 + 10*sin(x*y), (x,-4,4), (y,-4,4),
            ....:              plot_points=121, cmap='hsv')
            Graphics object consisting of 1 graphics primitive
        """
        from sage.rings.integer import Integer
        options = self.options()
        fill = options['fill']
        contours = options['contours']
        if 'cmap' in options:
            cmap = get_cmap(options['cmap'])
        elif fill or contours is None:
            cmap = get_cmap('gray')
        else:
            if isinstance(contours, (int, Integer)):
                cmap = get_cmap([(i, i, i)
                                 for i in xsrange(0, 1, 1 / contours)])
            else:
                step = 1 / Integer(len(contours))
                cmap = get_cmap([(i, i, i) for i in xsrange(0, 1, step)])

        x0, x1 = float(self.xrange[0]), float(self.xrange[1])
        y0, y1 = float(self.yrange[0]), float(self.yrange[1])

        if isinstance(contours, (int, Integer)):
            contours = int(contours)

        CSF = None
        if fill:
            if contours is None:
                CSF = subplot.contourf(self.xy_data_array, cmap=cmap,
                                       extent=(x0, x1, y0, y1))
            else:
                CSF = subplot.contourf(self.xy_data_array, contours, cmap=cmap,
                                       extent=(x0, x1, y0, y1), extend='both')

        linewidths = options.get('linewidths', None)
        if isinstance(linewidths, (int, Integer)):
            linewidths = int(linewidths)
        elif isinstance(linewidths, (list, tuple)):
            linewidths = tuple(int(x) for x in linewidths)

        from sage.plot.misc import get_matplotlib_linestyle
        linestyles = options.get('linestyles', None)
        if isinstance(linestyles, (list, tuple)):
            linestyles = [get_matplotlib_linestyle(i, 'long')
                          for i in linestyles]
        else:
            linestyles = get_matplotlib_linestyle(linestyles, 'long')
        if contours is None:
            CS = subplot.contour(self.xy_data_array, cmap=cmap,
                                 extent=(x0, x1, y0, y1),
                                 linewidths=linewidths, linestyles=linestyles)
        else:
            CS = subplot.contour(self.xy_data_array, contours, cmap=cmap,
                                 extent=(x0, x1, y0, y1),
                                 linewidths=linewidths, linestyles=linestyles)
        if options.get('labels', False):
            label_options = options['label_options']
            label_options['fontsize'] = int(label_options['fontsize'])
            if fill and label_options is None:
                label_options['inline'] = False
            subplot.clabel(CS, **label_options)
        if options.get('colorbar', False):
            colorbar_options = options['colorbar_options']
            from matplotlib import colorbar
            cax, kwds = colorbar.make_axes_gridspec(subplot, **colorbar_options)
            if CSF is None:
                cb = colorbar.Colorbar(cax, CS, **kwds)
            else:
                cb = colorbar.Colorbar(cax, CSF, **kwds)
                cb.add_lines(CS)


@suboptions('colorbar', orientation='vertical', format=None, spacing='uniform')
@suboptions('label', fontsize=9, colors='blue', inline=None, inline_spacing=3,
            fmt="%1.2f")
@options(plot_points=100, fill=True, contours=None, linewidths=None,
         linestyles=None, labels=False, frame=True, axes=False, colorbar=False,
         legend_label=None, aspect_ratio=1, region=None)
def contour_plot(f, xrange, yrange, **options):
    r"""
    ``contour_plot`` takes a function of two variables, `f(x,y)`
    and plots contour lines of the function over the specified
    ``xrange`` and ``yrange`` as demonstrated below.

    ``contour_plot(f, (xmin,xmax), (ymin,ymax), ...)``

    INPUT:

    - ``f`` -- a function of two variables

    - ``(xmin, xmax)`` -- 2-tuple, the range of ``x`` values OR 3-tuple
      ``(x,xmin,xmax)``

    - ``(ymin, ymax)`` -- 2-tuple, the range of ``y`` values OR 3-tuple
      ``(y,ymin,ymax)``

    The following inputs must all be passed in as named parameters:

    - ``plot_points`` -- integer (default: 100); number of points to plot
      in each direction of the grid.  For old computers, 25 is fine, but
      should not be used to verify specific intersection points.

    - ``fill`` -- boolean (default: ``True``); whether to color in the area
      between contour lines

    - ``cmap`` -- a colormap (default: ``'gray'``), the name of
      a predefined colormap, a list of colors or an instance of a matplotlib
      Colormap. Type: ``import matplotlib.cm; matplotlib.cm.datad.keys()``
      for available colormap names.

    - ``contours`` -- integer or list of numbers (default: ``None``):
      If a list of numbers is given, then this specifies the contour levels
      to use.  If an integer is given, then this many contour lines are
      used, but the exact levels are determined automatically. If ``None``
      is passed (or the option is not given), then the number of contour
      lines is determined automatically, and is usually about 5.

    - ``linewidths`` -- integer or list of integer (default: ``None``), if
      a single integer all levels will be of the width given,
      otherwise the levels will be plotted with the width in the order
      given.  If the list is shorter than the number of contours, then
      the widths will be repeated cyclically.

    - ``linestyles`` -- string or list of strings (default: ``None``), the
      style of the lines to be plotted, one of: ``'solid'``, ``'dashed'``,
      ``'dashdot'``, ``'dotted'``, respectively ``'-'``, ``'--'``,
      ``'-.'``, ``':'``.  If the list is shorter than the number of
      contours, then the styles will be repeated cyclically.

    - ``labels`` -- boolean (default: ``False``); show level labels or not

      The following options are to adjust the style and placement of
      labels, they have no effect if no labels are shown.

      - ``label_fontsize`` -- integer (default: 9); the font size of
        the labels

      - ``label_colors`` -- string or sequence of colors (default:
        None) If a string, gives the name of a single color with which
        to draw all labels.  If a sequence, gives the colors of the
        labels.  A color is a string giving the name of one or a
        3-tuple of floats.

      - ``label_inline`` -- boolean (default: ``False`` if fill is True,
        otherwise True), controls whether the underlying contour is
        removed or not.

      - ``label_inline_spacing`` -- integer (default: 3); when inline,
        this is the amount of contour that is removed from each side,
        in pixels.

      - ``label_fmt`` -- a format string (default: ``"%1.2f"``), this is
        used to get the label text from the level.  This can also be a
        dictionary with the contour levels as keys and corresponding
        text string labels as values.  It can also be any callable which
        returns a string when called with a numeric contour level.

    - ``colorbar`` -- boolean (default: ``False``); show a colorbar or not

      The following options are to adjust the style and placement of
      colorbars.  They have no effect if a colorbar is not shown.

      - ``colorbar_orientation`` -- string (default: ``'vertical'``),
        controls placement of the colorbar, can be either 'vertical'
        or 'horizontal'

      - ``colorbar_format`` -- a format string, this is used to format
        the colorbar labels

      - ``colorbar_spacing`` -- string (default: ``'proportional'``); if
        'proportional', make the contour divisions proportional to
        values.  If 'uniform', space the colorbar divisions uniformly,
        without regard for numeric values.

    - ``legend_label`` -- the label for this item in the legend

    - ``region`` -- (default: ``None``) if region is given, it must be a function
        of two variables. Only segments of the surface where region(x,y)
        returns a number >0 will be included in the plot.

    .. WARNING::

        Due to an implementation detail in matplotlib, single-contour
        plots whose data all lie on one side of the sole contour may
        not be plotted correctly. We attempt to detect this situation
        and to produce something better than an empty plot when it
        happens; a ``UserWarning`` is emitted in that case.

    EXAMPLES:

    Here we plot a simple function of two variables.  Note that
    since the input function is an expression, we need to explicitly
    declare the variables in 3-tuples for the range::

        sage: x,y = var('x,y')
        sage: contour_plot(cos(x^2 + y^2), (x,-4,4), (y,-4,4))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        g = contour_plot(cos(x**2 + y**2), (x,-4,4), (y,-4,4))
        sphinx_plot(g)

    Here we change the ranges and add some options::

        sage: x,y = var('x,y')
        sage: contour_plot((x^2) * cos(x*y), (x,-10,5), (y,-5,5), fill=False, plot_points=150)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        g = contour_plot((x**2) * cos(x*y), (x,-10,5), (y,-5,5), fill=False, plot_points=150)
        sphinx_plot(g)

    An even more complicated plot::

        sage: x,y = var('x,y')
        sage: contour_plot(sin(x^2+y^2) * cos(x) * sin(y), (x,-4,4), (y,-4,4), plot_points=150)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        g = contour_plot(sin(x**2+y**2) * cos(x) * sin(y), (x,-4,4), (y,-4,4),plot_points=150)
        sphinx_plot(g)

    Some elliptic curves, but with symbolic endpoints.  In the first
    example, the plot is rotated 90 degrees because we switch the
    variables `x`, `y`::

        sage: x,y = var('x,y')
        sage: contour_plot(y^2 + 1 - x^3 - x, (y,-pi,pi), (x,-pi,pi))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        g = contour_plot(y**2 + 1 - x**3 - x, (y,-pi,pi), (x,-pi,pi))
        sphinx_plot(g)

    ::

        sage: contour_plot(y^2 + 1 - x^3 - x, (x,-pi,pi), (y,-pi,pi))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        g = contour_plot(y**2 + 1 - x**3 - x, (x,-pi,pi), (y,-pi,pi))
        sphinx_plot(g)

    We can play with the contour levels::

        sage: x,y = var('x,y')
        sage: f(x,y) = x^2 + y^2
        sage: contour_plot(f, (-2,2), (-2,2))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        def f(x,y): return x**2 + y**2
        g = contour_plot(f, (-2,2), (-2,2))
        sphinx_plot(g)

    ::

        sage: contour_plot(f, (-2,2), (-2,2), contours=2, cmap=[(1,0,0), (0,1,0), (0,0,1)])
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        def f(x,y): return x**2 + y**2
        g = contour_plot(f, (-2, 2), (-2, 2), contours=2, cmap=[(1,0,0), (0,1,0), (0,0,1)])
        sphinx_plot(g)

    ::

        sage: contour_plot(f, (-2,2), (-2,2),
        ....:              contours=(0.1,1.0,1.2,1.4), cmap='hsv')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        def f(x,y): return x**2 + y**2
        g = contour_plot(f, (-2,2), (-2,2), contours=(0.1,1.0,1.2,1.4), cmap='hsv')
        sphinx_plot(g)

    ::

        sage: contour_plot(f, (-2,2), (-2,2), contours=(1.0,), fill=False)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        def f(x,y): return x**2 + y**2
        g = contour_plot(f, (-2,2), (-2,2), contours=(1.0,), fill=False)
        sphinx_plot(g)

    ::

        sage: contour_plot(x - y^2, (x,-5,5), (y,-3,3), contours=[-4,0,1])
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        g = contour_plot(x - y**2, (x,-5,5), (y,-3,3), contours=[-4,0,1])
        sphinx_plot(g)

    We can change the style of the lines::

        sage: contour_plot(f, (-2,2), (-2,2), fill=False, linewidths=10)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        def f(x,y): return x**2 + y**2
        g = contour_plot(f, (-2,2), (-2,2), fill=False, linewidths=10)
        sphinx_plot(g)


    ::

        sage: contour_plot(f, (-2,2), (-2,2), fill=False, linestyles='dashdot')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        def f(x,y): return x**2 + y**2
        g = contour_plot(f, (-2,2), (-2,2), fill=False, linestyles='dashdot')
        sphinx_plot(g)

    ::

        sage: P = contour_plot(x^2 - y^2, (x,-3,3), (y,-3,3),
        ....:                  contours=[0,1,2,3,4], linewidths=[1,5],
        ....:                  linestyles=['solid','dashed'], fill=False)
        sage: P
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        P = contour_plot(x**2 - y**2, (x,-3,3), (y,-3,3),
                         contours=[0,1,2,3,4], linewidths=[1,5],
                         linestyles=['solid','dashed'], fill=False)
        sphinx_plot(P)

    ::

        sage: P = contour_plot(x^2 - y^2, (x,-3,3), (y,-3,3),
        ....:                  contours=[0,1,2,3,4], linewidths=[1,5],
        ....:                  linestyles=['solid','dashed'])
        sage: P
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        P = contour_plot(x**2 - y**2, (x,-3,3), (y,-3,3),
                         contours=[0,1,2,3,4], linewidths=[1,5],
                         linestyles=['solid','dashed'])
        sphinx_plot(P)

    ::

        sage: P = contour_plot(x^2 - y^2, (x,-3,3), (y,-3,3),
        ....:                  contours=[0,1,2,3,4], linewidths=[1,5],
        ....:                  linestyles=['-',':'])
        sage: P
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        P = contour_plot(x**2 - y**2, (x,-3,3), (y,-3,3),
                         contours=[0,1,2,3,4], linewidths=[1,5],
                         linestyles=['-',':'])
        sphinx_plot(P)

    We can add labels and play with them::

        sage: contour_plot(y^2 + 1 - x^3 - x, (x,-pi,pi), (y,-pi,pi),
        ....:              fill=False, cmap='hsv', labels=True)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        P = contour_plot(y**2 + 1 - x**3 - x, (x,-pi,pi), (y,-pi,pi),
                         fill=False, cmap='hsv', labels=True)
        sphinx_plot(P)

    ::

        sage: P=contour_plot(y^2 + 1 - x^3 - x, (x,-pi,pi), (y,-pi,pi),
        ....:                fill=False, cmap='hsv',
        ....:                labels=True, label_fmt="%1.0f",
        ....:                label_colors='black')
        sage: P
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        P=contour_plot(y**2 + 1 - x**3 - x, (x,-pi,pi), (y,-pi,pi),
                       fill=False, cmap='hsv',
                       labels=True, label_fmt="%1.0f",
                        label_colors='black')
        sphinx_plot(P)

    ::

        sage: P = contour_plot(y^2 + 1 - x^3 - x, (x,-pi,pi), (y,-pi,pi),
        ....:                  fill=False, cmap='hsv', labels=True,
        ....:                  contours=[-4,0,4],
        ....:                  label_fmt={-4:"low", 0:"medium", 4: "hi"},
        ....:                  label_colors='black')
        sage: P
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        P = contour_plot(y**2 + 1 - x**3 - x, (x,-pi,pi), (y,-pi,pi),
                         fill=False, cmap='hsv', labels=True,
                         contours=[-4,0,4],
                         label_fmt={-4:"low", 0:"medium", 4: "hi"},
                         label_colors='black')
        sphinx_plot(P)

    ::

        sage: P = contour_plot(y^2 + 1 - x^3 - x, (x,-pi,pi), (y,-pi,pi),
        ....:                fill=False, cmap='hsv', labels=True,
        ....:                contours=[-4,0,4], label_fmt=lambda x: "$z=%s$"%x,
        ....:                label_colors='black', label_inline=True,
        ....:                label_fontsize=12)
        sage: P
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        P = contour_plot(y**2 + 1 - x**3 - x, (x,-pi,pi), (y,-pi,pi),
                         fill=False, cmap='hsv', labels=True,
                         contours=[-4,0,4], label_fmt=lambda x: "$z=%s$"%x,
                         label_colors='black', label_inline=True,
                         label_fontsize=12)
        sphinx_plot(P)

    ::

        sage: P = contour_plot(y^2 + 1 - x^3 - x, (x,-pi,pi), (y,-pi,pi),
        ....:                  fill=False, cmap='hsv', labels=True,
        ....:                  label_fontsize=18)
        sage: P
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        P = contour_plot(y**2 + 1 - x**3 - x, (x,-pi,pi), (y,-pi,pi),
                         fill=False, cmap='hsv', labels=True,
                         label_fontsize=18)
        sphinx_plot(P)

    ::

        sage: P = contour_plot(y^2 + 1 - x^3 - x, (x,-pi,pi), (y,-pi,pi),
        ....:                fill=False, cmap='hsv', labels=True,
        ....:                label_inline_spacing=1)
        sage: P
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        P = contour_plot(y**2 + 1 - x**3 - x, (x,-pi,pi), (y,-pi,pi),
                         fill=False, cmap='hsv', labels=True,
                         label_inline_spacing=1)
        sphinx_plot(P)

    ::

        sage: P = contour_plot(y^2 + 1 - x^3 - x, (x,-pi,pi), (y,-pi,pi),
        ....:                  fill=False, cmap='hsv', labels=True,
        ....:                  label_inline=False)
        sage: P
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        P = contour_plot(y**2 + 1 - x**3 - x, (x,-pi,pi), (y,-pi,pi),
                         fill=False, cmap='hsv', labels=True,
                         label_inline=False)
        sphinx_plot(P)

    We can change the color of the labels if so desired::

        sage: contour_plot(f, (-2,2), (-2,2), labels=True, label_colors='red')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        def f(x,y): return x**2 + y**2
        g = contour_plot(f, (-2,2), (-2,2), labels=True, label_colors='red')
        sphinx_plot(g)

    We can add a colorbar as well::

        sage: f(x, y) = x^2 + y^2
        sage: contour_plot(f, (x,-3,3), (y,-3,3), colorbar=True)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        def f(x,y): return x**2 + y**2
        g = contour_plot(f, (x,-3,3), (y,-3,3), colorbar=True)
        sphinx_plot(g)

    ::

        sage: contour_plot(f, (x,-3,3), (y,-3,3), colorbar=True, colorbar_orientation='horizontal')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        def f(x,y): return x**2 + y**2
        g = contour_plot(f, (x,-3,3), (y,-3,3), colorbar=True, colorbar_orientation='horizontal')
        sphinx_plot(g)

    ::

        sage: contour_plot(f, (x,-3,3), (y,-3,3), contours=[-2,-1,4], colorbar=True)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        def f(x,y): return x**2 + y**2
        g = contour_plot(f, (x,-3,3), (y,-3,3), contours=[-2,-1,4],
                        colorbar=True)
        sphinx_plot(g)

    ::

        sage: contour_plot(f, (x,-3,3), (y,-3,3), contours=[-2,-1,4],
        ....:              colorbar=True, colorbar_spacing='uniform')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        def f(x,y): return x**2 + y**2
        g = contour_plot(f, (x,-3,3), (y,-3,3), contours=[-2,-1,4],
                         colorbar=True, colorbar_spacing='uniform')
        sphinx_plot(g)

    ::

        sage: contour_plot(f, (x,-3,3), (y,-3,3), contours=[0,2,3,6],
        ....:              colorbar=True, colorbar_format='%.3f')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        def f(x,y): return x**2 + y**2
        g = contour_plot(f, (x,-3,3), (y,-3,3), contours=[0,2,3,6],
                         colorbar=True, colorbar_format='%.3f')
        sphinx_plot(g)

    ::

        sage: contour_plot(f, (x,-3,3), (y,-3,3), labels=True,
        ....:              label_colors='red', contours=[0,2,3,6],
        ....:              colorbar=True)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        def f(x,y): return x**2 + y**2
        g = contour_plot(f, (x,-3,3), (y,-3,3), labels=True,
                         label_colors='red', contours=[0,2,3,6],
                         colorbar=True)
        sphinx_plot(g)

    ::

        sage: contour_plot(f, (x,-3,3), (y,-3,3), cmap='winter',
        ....:              contours=20, fill=False, colorbar=True)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        def f(x,y): return x**2 + y**2
        g = contour_plot(f, (x,-3,3), (y,-3,3), cmap='winter',
                         contours=20, fill=False, colorbar=True)
        sphinx_plot(g)

    This should plot concentric circles centered at the origin::

        sage: x,y = var('x,y')
        sage: contour_plot(x^2 + y^2-2,(x,-1,1), (y,-1,1))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        g = contour_plot(x**2 + y**2-2,(x,-1,1), (y,-1,1))
        sphinx_plot(g)


    Extra options will get passed on to show(), as long as they are valid::

        sage: f(x,y) = cos(x) + sin(y)
        sage: contour_plot(f, (0,pi), (0,pi), axes=True)
        Graphics object consisting of 1 graphics primitive

    ::

        sage: contour_plot(f, (0,pi), (0,pi)).show(axes=True) # These are equivalent

    .. PLOT::

        x,y = var('x,y')
        def f(x,y): return cos(x) + sin(y)
        g = contour_plot(f, (0,pi), (0,pi), axes=True)
        sphinx_plot(g)

    One can also plot over a reduced region::

        sage: contour_plot(x**2 - y**2, (x,-2,2), (y,-2,2), region=x - y, plot_points=300)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        g = contour_plot(x**2 - y**2, (x,-2,2), (y,-2,2), region=x - y,
                         plot_points=300)
        sphinx_plot(g)

    Note that with ``fill=False`` and grayscale contours, there is the
    possibility of confusion between the contours and the axes, so use
    ``fill=False`` together with ``axes=True`` with caution::

        sage: contour_plot(f, (-pi,pi), (-pi,pi), fill=False, axes=True)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        def f(x,y): return cos(x) + sin(y)
        g = contour_plot(f, (-pi,pi), (-pi,pi), fill=False, axes=True)
        sphinx_plot(g)

    If you are plotting a sole contour and if all of your data lie on
    one side of it, then (as part of :issue:`21042`) a heuristic may be
    used to improve the result; in that case, a warning is emitted::

        sage: contour_plot(lambda x,y: abs(x^2-y^2), (-1,1), (-1,1),
        ....:              contours=[0], fill=False, cmap=['blue'])
        ...
        UserWarning: pathological contour plot of a function whose values
        all lie on one side of the sole contour; we are adding more plot
        points and perturbing your function values.
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            g = contour_plot(lambda x,y: abs(x**2-y**2), (-1,1), (-1,1),
                             contours=[0], fill=False, cmap=['blue'])
        sphinx_plot(g)

    Constant functions (with a single contour) can be plotted as well;
    this was not possible before :issue:`21042`::

        sage: contour_plot(lambda x,y: 0, (-1,1), (-1,1),
        ....:              contours=[0], fill=False, cmap=['blue'])
        ...Graphics object consisting of 1 graphics primitive

    .. PLOT::

        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            g = contour_plot(lambda x,y: 0, (-1,1), (-1,1),
                             contours=[0], fill=False, cmap=['blue'])
        sphinx_plot(g)

    TESTS:

    To check that :issue:`5221` is fixed, note that this has three curves, not
    two::

        sage: x,y = var('x,y')
        sage: contour_plot(x - y^2, (x,-5,5), (y,-3,3),
        ....:              contours=[-4,-2,0], fill=False)
        Graphics object consisting of 1 graphics primitive

    Check that :issue:`18074` is fixed::

        sage: contour_plot(0, (0,1), (0,1))
        ...Graphics object consisting of 1 graphics primitive

    Domain points in :issue:`11648` with complex output are now skipped::

        sage: x,y = SR.var('x,y', domain='real')
        sage: contour_plot(log(x) + log(y), (-1, 5), (-1, 5))
        Graphics object consisting of 1 graphics primitive
    """
    from sage.plot.all import Graphics
    from sage.plot.misc import setup_for_eval_on_grid

    region = options.pop('region')
    ev = [f] if region is None else [f, region]

    F, ranges = setup_for_eval_on_grid(ev, [xrange, yrange],
                                       options['plot_points'])
    h = F[0]
    xrange, yrange = (r[:2] for r in ranges)

    xy_data_array = [[h(x, y) for x in xsrange(*ranges[0],
                                               include_endpoint=True)]
                     for y in xsrange(*ranges[1], include_endpoint=True)]

    g = Graphics()

    # Reset aspect_ratio to 'automatic' in case scale is 'semilog[xy]'.
    # Otherwise matplotlib complains.
    scale = options.get('scale', None)
    if isinstance(scale, (list, tuple)):
        scale = scale[0]
    if scale in ('semilogy', 'semilogx'):
        options['aspect_ratio'] = 'automatic'

    g._set_extra_kwds(Graphics._extract_kwds_for_show(options,
                                                      ignore=['xmin', 'xmax']))

    # Was a single contour level explicitly given? If "contours" is
    # the integer 1, then there will be a single level, but we can't
    # know what it is because it's determined within matplotlib's
    # "contour" or "contourf" function. So we punt in that case. If
    # there's a single contour and fill=True, we fall through to let
    # matplotlib complain that "Filled contours require at least 2
    # levels."
    if (isinstance(options["contours"], (list, tuple))
        and len(options["contours"]) == 1
        and options.get("fill") is False):
        # When there's only one level (say, zero), matplotlib doesn't
        # handle it well. If all of the data lie on one side of that
        # level -- for example, if f(x,y) >= 0 for all x,y -- then it
        # will fail to plot the points where f(x,y) == 0. This is
        # especially catastrophic for implicit_plot(), which tries to
        # do just that. Here we handle that special case: if there's
        # only one level, and if all of the data lie on one side of
        # it, we perturb the data a bit so that they don't. The resulting
        # plots don't look great, but they're not empty, which is an
        # improvement.
        import numpy as np
        dx = ranges[0][2]
        dy = ranges[1][2]
        z0 = options["contours"][0]

        # This works OK for the examples in the doctests, but basing
        # it off the plot scale rather than how fast the function
        # changes can never be truly satisfactory.
        tol = max(dx, dy) / 4.0
        xy_data_array = np.ma.asarray(xy_data_array, dtype=float)

        # Special case for constant functions. This is needed because
        # otherwise the forthcoming perturbation trick will take values
        # like 0,0,0... and perturb them to -tol, -tol, -tol... which
        # doesn't work for the same reason 0,0,0... doesn't work.
        if np.all(np.abs(xy_data_array - z0) <= tol):
            # Up to our tolerance, this is the const_z0 function.
            # ...make it actually the const_z0 function.
            xy_data_array.fill(z0)

            # We're going to set fill=True in a moment, so we need to
            # prepend an entry to the cmap so that the user's original
            # cmap winds up in the right place.
            if "cmap" in options:
                if isinstance(options["cmap"], (list, tuple)):
                    oldcmap = options["cmap"][0]
                else:
                    oldcmap = options["cmap"]
            else:
                # The docs promise this as the default.
                oldcmap = "gray"

            # Trick matplotlib into plotting all of the points (minus
            # those masked) by using a single, filled contour that
            # covers the entire plotting surface.
            options["cmap"] = ["white", oldcmap]
            options["contours"] = (z0 - 1, z0)
            options["fill"] = True
        else:
            # The "c" constant is set to plus/minus one to handle both
            # of the "all values greater than z0" and "all values less
            # than z0" cases at once.
            c = 1
            if np.all(xy_data_array <= z0):
                xy_data_array *= -1
                c = -1
            # Now we check if (a) all of the data lie on one side of
            # z0, and (b) if perturbing the data will actually help by
            # moving anything across z0.
            if (np.all(xy_data_array >= z0) and
                np.any(xy_data_array - z0 < tol)):

                from warnings import warn
                warn("pathological contour plot of a function whose "
                     "values all lie on one side of the sole contour; "
                     "we are adding more plot points and perturbing "
                     "your function values.")

                # The choice of "4" here is not based on much of anything.
                # It works well enough for the examples in the doctests.
                if not isinstance(options["plot_points"], (list, tuple)):
                    options["plot_points"] = (options["plot_points"],
                                              options["plot_points"])
                    options["plot_points"] = (options["plot_points"][0] * 4,
                                              options["plot_points"][1] * 4)

                # Re-plot with more points...
                F, ranges = setup_for_eval_on_grid(ev, [xrange, yrange],
                                                   options['plot_points'])
                h = F[0]
                xrange, yrange = (r[:2] for r in ranges)

                # ...and a function whose values are shifted towards
                # z0 by "tol".
                xy_data_array = [[h(x, y) - c * tol
                                  for x in xsrange(*ranges[0],
                                                   include_endpoint=True)]
                                 for y in xsrange(*ranges[1],
                                                  include_endpoint=True)]

    if region is not None:
        import numpy

        xy_data_array = numpy.ma.asarray(xy_data_array, dtype=float)

        m = F[1]

        mask = numpy.asarray([[m(x, y) <= 0
                               for x in xsrange(*ranges[0],
                                                include_endpoint=True)]
                              for y in xsrange(*ranges[1],
                                               include_endpoint=True)],
                             dtype=bool)

        xy_data_array[mask] = numpy.ma.masked

    g.add_primitive(ContourPlot(xy_data_array, xrange, yrange, options))

    return g


@options(plot_points=150, contours=(0,), fill=False, cmap=["blue"])
def implicit_plot(f, xrange, yrange, **options):
    r"""
    ``implicit_plot`` takes a function of two variables, `f(x, y)`
    and plots the curve `f(x,y) = 0` over the specified
    ``xrange`` and ``yrange`` as demonstrated below.

    ``implicit_plot(f, (xmin,xmax), (ymin,ymax), ...)``

    ``implicit_plot(f, (x,xmin,xmax), (y,ymin,ymax), ...)``

    INPUT:

    - ``f`` -- a function of two variables or equation in two variables

    - ``(xmin, xmax)`` -- 2-tuple, the range of ``x``
      values or ``(x,xmin,xmax)``

    - ``(ymin, ymax)`` -- 2-tuple, the range of ``y``
      values or ``(y,ymin,ymax)``

    The following inputs must all be passed in as named parameters:

    - ``plot_points`` -- integer (default: 150); number of points to plot
      in each direction of the grid

    - ``fill`` -- boolean (default: ``False``); if ``True``, fill the region
      `f(x, y) < 0`

    - ``fillcolor`` -- string (default: ``'blue'``); the color of the region
      where `f(x,y) < 0` if ``fill = True``. Colors are defined in
      :mod:`sage.plot.colors`; try ``colors?`` to see them all.

    - ``linewidth`` -- integer (default: ``None``); if a single integer all levels
      will be of the width given, otherwise the levels will be plotted with the
      widths in the order given.

    - ``linestyle`` -- string (default: ``None``); the style of the line to be
      plotted, one of: ``'solid'``, ``'dashed'``, ``'dashdot'`` or
      ``'dotted'``, respectively ``'-'``, ``'--'``, ``'-.'``, or ``':'``.

    - ``color`` -- string (default: ``'blue'``); the color of the plot. Colors
      are defined in :mod:`sage.plot.colors`; try ``colors?`` to see them all.
      If ``fill = True``, then this sets only the color of the border of the
      plot. See ``fillcolor`` for setting the color of the fill region.

    - ``legend_label`` -- the label for this item in the legend

    - ``base`` -- (default: 10) the base of the logarithm if
      a logarithmic scale is set. This must be greater than 1. The base
      can be also given as a list or tuple ``(basex, basey)``.
      ``basex`` sets the base of the logarithm along the horizontal
      axis and ``basey`` sets the base along the vertical axis.

    - ``scale`` -- (default: ``'linear'``) string. The scale of the axes.
      Possible values are ``'linear'``, ``'loglog'``, ``'semilogx'``,
      ``'semilogy'``.

      The scale can be also be given as single argument that is a list
      or tuple ``(scale, base)`` or ``(scale, basex, basey)``.

      The ``'loglog'`` scale sets both the horizontal and vertical axes to
      logarithmic scale. The ``'semilogx'`` scale sets the horizontal axis
      to logarithmic scale. The ``'semilogy'`` scale sets the vertical axis
      to logarithmic scale. The ``'linear'`` scale is the default value
      when :class:`~sage.plot.graphics.Graphics` is initialized.

    .. WARNING::

        Due to an implementation detail in matplotlib, implicit plots
        whose data are all nonpositive or nonnegative may not be
        plotted correctly. We attempt to detect this situation and to
        produce something better than an empty plot when it happens; a
        ``UserWarning`` is emitted in that case.

    EXAMPLES:

    A simple circle with a radius of 2. Note that
    since the input function is an expression, we need to explicitly
    declare the variables in 3-tuples for the range::

        sage: var("x y")
        (x, y)
        sage: implicit_plot(x^2 + y^2 - 2, (x,-3,3), (y,-3,3))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x, y =var("x y")
        g = implicit_plot(x**2 + y**2 - 2, (x,-3,3), (y,-3,3))
        sphinx_plot(g)

    We can do the same thing, but using a callable function so we do not
    need to explicitly define the variables in the ranges. We also fill
    the inside::

        sage: f(x,y) = x^2 + y^2 - 2
        sage: implicit_plot(f, (-3,3), (-3,3), fill=True, plot_points=500) # long time
        Graphics object consisting of 2 graphics primitives

    .. PLOT::

        def f(x,y): return x**2 + y**2 - 2
        g = implicit_plot(f, (-3,3), (-3,3), fill=True, plot_points=500)
        sphinx_plot(g)

    The same circle but with a different line width::

        sage: implicit_plot(f, (-3,3), (-3,3), linewidth=6)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        def f(x,y): return x**2 + y**2 - 2
        g = implicit_plot(f, (-3,3), (-3,3), linewidth=6)
        sphinx_plot(g)

    Again the same circle but this time with a dashdot border::

        sage: implicit_plot(f, (-3,3), (-3,3), linestyle='dashdot')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x, y =var("x y")
        def f(x,y): return x**2 + y**2 - 2
        g = implicit_plot(f, (-3,3), (-3,3), linestyle='dashdot')
        sphinx_plot(g)

    The same circle with different line and fill colors::

        sage: implicit_plot(f, (-3,3), (-3,3), color='red',  # long time
        ....:               fill=True, fillcolor='green',
        ....:               plot_points=500)
        Graphics object consisting of 2 graphics primitives

    .. PLOT::

        def f(x,y): return x**2 + y**2 - 2
        g = implicit_plot(f, (-3,3), (-3,3), color='red', fill=True, fillcolor='green',
                                             plot_points=500)
        sphinx_plot(g)

    You can also plot an equation::

        sage: var("x y")
        (x, y)
        sage: implicit_plot(x^2 + y^2 == 2, (x,-3,3), (y,-3,3))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x, y =var("x y")
        g = implicit_plot(x**2 + y**2 == 2, (x,-3,3), (y,-3,3))
        sphinx_plot(g)

    You can even change the color of the plot::

        sage: implicit_plot(x^2 + y^2 == 2, (x,-3,3), (y,-3,3), color='red')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x, y =var("x y")
        g = implicit_plot(x**2 + y**2 == 2, (x,-3,3), (y,-3,3), color='red')
        sphinx_plot(g)

    The color of the fill region can be changed::

        sage: implicit_plot(x**2 + y**2 == 2, (x,-3,3), (y,-3,3), fill=True, fillcolor='red')
        Graphics object consisting of 2 graphics primitives

    .. PLOT::

        x, y =var("x y")
        g = implicit_plot(x**2 + y**2 == 2, (x,-3,3), (y,-3,3), fill=True, fillcolor='red')
        sphinx_plot(g)

    Here is a beautiful (and long) example which also tests that all
    colors work with this::

        sage: G = Graphics()
        sage: counter = 0
        sage: for col in colors.keys():  # long time
        ....:     G += implicit_plot(x^2 + y^2 == 1 + counter*.1, (x,-4,4),(y,-4,4), color=col)
        ....:     counter += 1
        sage: G  # long time
        Graphics object consisting of 148 graphics primitives

    .. PLOT::

        x, y = var("x y")
        G = Graphics()
        counter = 0
        for col in colors.keys():
            G += implicit_plot(x**2 + y**2 == 1 + counter*.1, (x,-4,4), (y,-4,4), color=col)
            counter += 1
        sphinx_plot(G)

    We can define a level-`n` approximation of the boundary of the
    Mandelbrot set::

        sage: def mandel(n):
        ....:     c = polygen(CDF, 'c')
        ....:     z = 0
        ....:     for i in range(n):
        ....:         z = z*z + c
        ....:     def f(x,y):
        ....:         val = z(CDF(x, y))
        ....:         return val.norm() - 4
        ....:     return f

    The first-level approximation is just a circle::

        sage: implicit_plot(mandel(1), (-3,3), (-3,3))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        def mandel(n):
             c = polygen(CDF, 'c')
             z = 0
             for i in range(n):
                 z = z*z + c
             def f(x,y):
                 val = z(CDF(x, y))
                 return val.norm() - 4
             return f
        g = implicit_plot(mandel(1), (-3,3), (-3,3))
        sphinx_plot(g)

    A third-level approximation starts to get interesting::

        sage: implicit_plot(mandel(3), (-2,1), (-1.5,1.5))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        def mandel(n):
             c = polygen(CDF, 'c')
             z = 0
             for i in range(n):
                 z = z*z + c
             def f(x,y):
                 val = z(CDF(x, y))
                 return val.norm() - 4
             return f
        g = implicit_plot(mandel(3), (-2,1), (-1.5,1.5))
        sphinx_plot(g)

    The seventh-level approximation is a degree 64 polynomial, and
    ``implicit_plot`` does a pretty good job on this part of the curve.
    (``plot_points=200`` looks even better, but it takes over a second.)

    ::

        sage: implicit_plot(mandel(7), (-0.3, 0.05), (-1.15, -0.9), plot_points=50)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        def mandel(n):
             c = polygen(CDF, 'c')
             z = 0
             for i in range(n):
                 z = z*z + c
             def f(x,y):
                 val = z(CDF(x, y))
                 return val.norm() - 4
             return f
        g = implicit_plot(mandel(7), (-0.3,0.05), (-1.15,-0.9),
                          plot_points=50)
        sphinx_plot(g)

    When making a filled implicit plot using a python function rather than a
    symbolic expression the user should increase the number of plot points to
    avoid artifacts::

        sage: implicit_plot(lambda x, y: x^2 + y^2 - 2, (x,-3,3),  # long time
        ....:               (y,-3,3), fill=True, plot_points=500)
        Graphics object consisting of 2 graphics primitives

    .. PLOT::

        x, y = var("x y")
        g = implicit_plot(lambda x, y: x**2 + y**2 - 2, (x,-3,3), (y,-3,3),
                          fill=True, plot_points=500)
        sphinx_plot(g)

    An example of an implicit plot on 'loglog' scale::

        sage: implicit_plot(x^2 + y^2 == 200, (x,1,200), (y,1,200), scale='loglog')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x, y = var("x y")
        g = implicit_plot(x**2 + y**2 == 200, (x,1,200), (y,1,200), scale='loglog')
        sphinx_plot(g)

    TESTS::

        sage: f(x,y) = x^2 + y^2 - 2
        sage: implicit_plot(f, (-3,3), (-3,3), fill=5)
        Traceback (most recent call last):
        ...
        ValueError: fill=5 is not supported

    To check that :issue:`9654` is fixed::

        sage: f(x,y) = x^2 + y^2 - 2
        sage: implicit_plot(f, (-3,3), (-3,3), rgbcolor=(1,0,0))
        Graphics object consisting of 1 graphics primitive
        sage: implicit_plot(f, (-3,3), (-3,3), color='green')
        Graphics object consisting of 1 graphics primitive
        sage: implicit_plot(f, (-3,3), (-3,3), rgbcolor=(1,0,0), color='green')
        Traceback (most recent call last):
        ...
        ValueError: only one of color or rgbcolor should be specified
    """
    from sage.structure.element import Expression
    if isinstance(f, Expression) and f.is_relational():
        if f.operator() != operator.eq:
            raise ValueError("input to implicit plot must be function "
                             "or equation")
        f = f.lhs() - f.rhs()
    linewidths = options.pop('linewidth', None)
    linestyles = options.pop('linestyle', None)

    if 'color' in options and 'rgbcolor' in options:
        raise ValueError('only one of color or rgbcolor should be specified')

    if 'color' in options:
        options['cmap'] = [options.pop('color', None)]
    elif 'rgbcolor' in options:
        options['cmap'] = [rgbcolor(options.pop('rgbcolor', None))]

    if options['fill'] is True:
        options.pop('fill')
        options.pop('contours', None)
        incol = options.pop('fillcolor', 'blue')
        bordercol = options.pop('cmap', [None])[0]
        from sage.structure.element import Expression
        if not isinstance(f, Expression):
            return region_plot(lambda x, y: f(x, y) < 0, xrange, yrange,
                               borderwidth=linewidths, borderstyle=linestyles,
                               incol=incol, bordercol=bordercol,
                               **options)
        return region_plot(f < 0, xrange, yrange, borderwidth=linewidths,
                           borderstyle=linestyles,
                           incol=incol, bordercol=bordercol,
                           **options)
    elif options['fill'] is False:
        options.pop('fillcolor', None)
        return contour_plot(f, xrange, yrange, linewidths=linewidths,
                            linestyles=linestyles, **options)
    else:
        raise ValueError("fill=%s is not supported" % options['fill'])


@options(plot_points=100, incol='blue', outcol=None, bordercol=None,
         borderstyle=None, borderwidth=None, frame=False, axes=True,
         legend_label=None, aspect_ratio=1, alpha=1)
def region_plot(f, xrange, yrange, **options):
    r"""
    ``region_plot`` takes a boolean function of two variables, `f(x, y)`
    and plots the region where f is ``True`` over the specified
    ``xrange`` and ``yrange`` as demonstrated below.

    ``region_plot(f, (xmin,xmax), (ymin,ymax), ...)``

    INPUT:

    - ``f`` -- boolean function or a list of boolean functions of
      two variables

    - ``(xmin, xmax)`` -- 2-tuple, the range of ``x`` values OR 3-tuple
      ``(x,xmin,xmax)``

    - ``(ymin, ymax)`` -- 2-tuple, the range of ``y`` values OR 3-tuple
      ``(y,ymin,ymax)``

    - ``plot_points`` -- integer (default: 100); number of points to plot
      in each direction of the grid

    - ``incol`` -- a color (default: ``'blue'``), the color inside the region

    - ``outcol`` -- a color (default: ``None``), the color of the outside
      of the region

    If any of these options are specified, the border will be shown as
    indicated, otherwise it is only implicit (with color ``incol``) as the
    border of the inside of the region.

    - ``bordercol`` -- a color (default: ``None``), the color of the border
       (``'black'`` if ``borderwidth`` or ``borderstyle`` is specified but
       not ``bordercol``)

    - ``borderstyle`` -- string (default: ``'solid'``); one of ``'solid'``,
      ``'dashed'``, ``'dotted'``, ``'dashdot'``, respectively ``'-'``,
      ``'--'``, ``':'``, ``'-.'``

    - ``borderwidth`` -- integer (default: ``None``); the width of the
      border in pixels

    - ``alpha`` -- (default: 1) how transparent the fill is; a number
      between 0 and 1

    - ``legend_label`` -- the label for this item in the legend

    - ``base`` -- (default: 10) the base of the logarithm if
      a logarithmic scale is set. This must be greater than 1. The base
      can be also given as a list or tuple ``(basex, basey)``.
      ``basex`` sets the base of the logarithm along the horizontal
      axis and ``basey`` sets the base along the vertical axis.

    - ``scale`` -- string (default: ``'linear'``); the scale of the axes.
      Possible values are ``'linear'``, ``'loglog'``, ``'semilogx'``,
      ``'semilogy'``.

      The scale can be also be given as single argument that is a list
      or tuple ``(scale, base)`` or ``(scale, basex, basey)``.

      The ``'loglog'`` scale sets both the horizontal and vertical axes to
      logarithmic scale. The ``'semilogx'`` scale sets the horizontal axis
      to logarithmic scale. The ``'semilogy'`` scale sets the vertical axis
      to logarithmic scale. The ``'linear'`` scale is the default value
      when :class:`~sage.plot.graphics.Graphics` is initialized.

    EXAMPLES:

    Here we plot a simple function of two variables::

        sage: x,y = var('x,y')
        sage: region_plot(cos(x^2 + y^2) <= 0, (x,-3,3), (y,-3,3))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        g = region_plot(cos(x**2 + y**2) <= 0, (x,-3,3), (y,-3,3))
        sphinx_plot(g)

    Here we play with the colors::

        sage: region_plot(x^2 + y^3 < 2, (x,-2,2), (y,-2,2), incol='lightblue', bordercol='gray')
        Graphics object consisting of 2 graphics primitives

    .. PLOT::

        x,y = var('x,y')
        g = region_plot(x**2 + y**3 < 2, (x,-2,2), (y,-2,2), incol='lightblue', bordercol='gray')
        sphinx_plot(g)

    An even more complicated plot, with dashed borders::

        sage: region_plot(sin(x) * sin(y) >= 1/4, (x,-10,10), (y,-10,10),
        ....:             incol='yellow', bordercol='black',
        ....:             borderstyle='dashed', plot_points=250)
        Graphics object consisting of 2 graphics primitives

    .. PLOT::

        x,y = var('x,y')
        g = region_plot(sin(x) * sin(y) >= 1/4, (x,-10,10), (y,-10,10),
                        incol='yellow', bordercol='black',
                        borderstyle='dashed', plot_points=250)
        sphinx_plot(g)

    A disk centered at the origin::

        sage: region_plot(x^2 + y^2 < 1, (x,-1,1), (y,-1,1))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        g = region_plot(x**2 + y**2 < 1, (x,-1,1), (y,-1,1))
        sphinx_plot(g)

    A plot with more than one condition (all conditions must be true for the
    statement to be true)::

        sage: region_plot([x^2 + y^2 < 1, x < y], (x,-2,2), (y,-2,2))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        g = region_plot([x**2 + y**2 < 1, x < y], (x,-2,2), (y,-2,2))
        sphinx_plot(g)

    Since it does not look very good, let us increase ``plot_points``::

        sage: region_plot([x^2 + y^2 < 1, x< y], (x,-2,2), (y,-2,2), plot_points=400)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        g = region_plot([x**2 + y**2 < 1, x < y], (x,-2,2), (y,-2,2), plot_points=400)
        sphinx_plot(g)

    To get plots where only one condition needs to be true, use a function.
    Using lambda functions, we definitely need the extra ``plot_points``::

        sage: region_plot(lambda x, y: x^2 + y^2 < 1 or x < y, (x,-2,2), (y,-2,2), plot_points=400)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x,y = var('x,y')
        g = region_plot(lambda x, y: x**2 + y**2 < 1 or x < y, (x,-2,2), (y,-2,2), plot_points=400)
        sphinx_plot(g)

    The first quadrant of the unit circle::

        sage: region_plot([y > 0, x > 0, x^2 + y^2 < 1], (x,-1.1,1.1), (y,-1.1,1.1), plot_points=400)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x, y = var("x y")
        g = region_plot([y > 0, x > 0, x**2 + y**2 < 1], (x,-1.1,1.1), (y,-1.1,1.1), plot_points=400)
        sphinx_plot(g)

    Here is another plot, with a huge border::

        sage: region_plot(x*(x-1)*(x+1) + y^2 < 0, (x,-3,2), (y,-3,3),
        ....:             incol='lightblue', bordercol='gray', borderwidth=10,
        ....:             plot_points=50)
        Graphics object consisting of 2 graphics primitives

    .. PLOT::

        x, y = var("x y")
        g = region_plot(x*(x-1)*(x+1) + y**2 < 0, (x,-3,2), (y,-3,3),
                          incol='lightblue', bordercol='gray', borderwidth=10,
                          plot_points=50)
        sphinx_plot(g)

    If we want to keep only the region where x is positive::

        sage: region_plot([x*(x-1)*(x+1) + y^2 < 0, x > -1], (x,-3,2), (y,-3,3),
        ....:             incol='lightblue', plot_points=50)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x, y =var("x y")
        g = region_plot([x*(x-1)*(x+1) + y**2 < 0, x > -1], (x,-3,2), (y,-3,3),
                        incol='lightblue', plot_points=50)
        sphinx_plot(g)

    Here we have a cut circle::

        sage: region_plot([x^2 + y^2 < 4, x > -1], (x,-2,2), (y,-2,2),
        ....:             incol='lightblue', bordercol='gray', plot_points=200)
        Graphics object consisting of 2 graphics primitives

    .. PLOT::

        x, y =var("x y")
        g = region_plot([x**2 + y**2 < 4, x > -1], (x,-2,2), (y,-2,2),
                        incol='lightblue', bordercol='gray', plot_points=200)
        sphinx_plot(g)

    The first variable range corresponds to the horizontal axis and
    the second variable range corresponds to the vertical axis::

        sage: s, t = var('s, t')
        sage: region_plot(s > 0, (t,-2,2), (s,-2,2))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        s, t = var('s, t')
        g = region_plot(s > 0, (t,-2,2), (s,-2,2))
        sphinx_plot(g)

    ::

        sage: region_plot(s>0,(s,-2,2),(t,-2,2))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        s, t = var('s, t')
        g = region_plot(s > 0, (s,-2,2), (t,-2,2))
        sphinx_plot(g)

    An example of a region plot in 'loglog' scale::

        sage: region_plot(x^2 + y^2 < 100, (x,1,10), (y,1,10), scale='loglog')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x, y = var("x y")
        g = region_plot(x**2 + y**2 < 100, (x,1,10), (y,1,10), scale='loglog')
        sphinx_plot(g)

    TESTS:

    To check that :issue:`16907` is fixed::

        sage: x, y = var('x, y')
        sage: disc1 = region_plot(x^2 + y^2 < 1, (x,-1,1), (y,-1,1), alpha=0.5)
        sage: disc2 = region_plot((x-0.7)^2 + (y-0.7)^2 < 0.5, (x,-2,2), (y,-2,2), incol='red', alpha=0.5)
        sage: disc1 + disc2
        Graphics object consisting of 2 graphics primitives

    To check that :issue:`18286` is fixed::

        sage: x, y = var('x, y')
        sage: region_plot([x == 0], (x,-1,1), (y,-1,1))
        Graphics object consisting of 1 graphics primitive
        sage: region_plot([x^2 + y^2 == 1, x < y], (x,-1,1), (y,-1,1))
        Graphics object consisting of 1 graphics primitive
    """
    from sage.plot.all import Graphics
    from sage.plot.misc import setup_for_eval_on_grid
    from sage.structure.element import Expression
    from warnings import warn
    import numpy

    plot_points = options['plot_points']
    incol = options.pop('incol')
    outcol = options.pop('outcol')
    bordercol = options.pop('bordercol')
    borderstyle = options.pop('borderstyle')
    borderwidth = options.pop('borderwidth')
    alpha = options.pop('alpha')

    if not isinstance(f, (list, tuple)):
        f = [f]

    feqs = [equify(g) for g in f
            if isinstance(g, Expression) and g.operator() is operator.eq
            and not equify(g).is_zero()]
    f = [equify(g) for g in f
         if not (isinstance(g, Expression) and g.operator() is operator.eq)]
    neqs = len(feqs)
    if neqs > 1:
        warn("There are at least 2 equations; "
             "If the region is degenerated to points, "
             "plotting might show nothing.")
        feqs = [sum([fn**2 for fn in feqs])]
        neqs = 1
    if neqs and not bordercol:
        bordercol = incol
    if not f:
        return implicit_plot(feqs[0], xrange, yrange, fill=False,
                             linewidth=borderwidth, linestyle=borderstyle,
                             color=bordercol, **options)
    f_all, ranges = setup_for_eval_on_grid(feqs + f,
                                           [xrange, yrange],
                                           plot_points)
    xrange, yrange = (r[:2] for r in ranges)

    xy_data_arrays = numpy.asarray([[[func(x, y)
                                      for x in xsrange(*ranges[0],
                                                       include_endpoint=True)]
                                     for y in xsrange(*ranges[1],
                                                      include_endpoint=True)]
                                    for func in f_all[neqs::]], dtype=float)
    xy_data_array = numpy.abs(xy_data_arrays.prod(axis=0))
    # Now we need to set entries to negative iff all
    # functions were negative at that point.
    neg_indices = (xy_data_arrays < 0).all(axis=0)
    xy_data_array[neg_indices] = -xy_data_array[neg_indices]

    from matplotlib.colors import ListedColormap
    incol = rgbcolor(incol)
    if outcol:
        outcol = rgbcolor(outcol)
        cmap = ListedColormap([incol, outcol])
        cmap.set_over(outcol, alpha=alpha)
    else:
        outcol = rgbcolor('white')
        cmap = ListedColormap([incol, outcol])
        cmap.set_over(outcol, alpha=0)
    cmap.set_under(incol, alpha=alpha)

    g = Graphics()

    # Reset aspect_ratio to 'automatic' in case scale is 'semilog[xy]'.
    # Otherwise matplotlib complains.
    scale = options.get('scale', None)
    if isinstance(scale, (list, tuple)):
        scale = scale[0]
    if scale in ('semilogy', 'semilogx'):
        options['aspect_ratio'] = 'automatic'

    g._set_extra_kwds(Graphics._extract_kwds_for_show(options,
                                                      ignore=['xmin', 'xmax']))

    if neqs == 0:
        g.add_primitive(ContourPlot(xy_data_array, xrange, yrange,
                                    dict(contours=[-1e-20, 0, 1e-20],
                                         cmap=cmap,
                                         fill=True, **options)))
    else:
        mask = numpy.asarray([[elt > 0 for elt in rows]
                              for rows in xy_data_array],
                             dtype=bool)
        xy_data_array = numpy.asarray([[f_all[0](x, y)
                                        for x in xsrange(*ranges[0],
                                                         include_endpoint=True)]
                                       for y in xsrange(*ranges[1],
                                                        include_endpoint=True)],
                                      dtype=float)
        xy_data_array[mask] = None
    if bordercol or borderstyle or borderwidth:
        cmap = [rgbcolor(bordercol)] if bordercol else ['black']
        linestyles = [borderstyle] if borderstyle else None
        linewidths = [borderwidth] if borderwidth else None
        g.add_primitive(ContourPlot(xy_data_array, xrange, yrange,
                                    dict(linestyles=linestyles,
                                         linewidths=linewidths,
                                         contours=[0], cmap=[bordercol],
                                         fill=False, **options)))

    return g


def equify(f):
    """
    Return the equation rewritten as a symbolic function to give
    negative values when ``True``, positive when ``False``.

    EXAMPLES::

        sage: from sage.plot.contour_plot import equify
        sage: var('x, y')
        (x, y)
        sage: equify(x^2 < 2)
        x^2 - 2
        sage: equify(x^2 > 2)
        -x^2 + 2
        sage: equify(x*y > 1)
        -x*y + 1
        sage: equify(y > 0)
        -y
        sage: f = equify(lambda x, y: x > y)
        sage: f(1, 2)
        1
        sage: f(2, 1)
        -1
    """
    from sage.calculus.all import symbolic_expression
    from sage.structure.element import Expression
    if not isinstance(f, Expression):
        return lambda x, y: -1 if f(x, y) else 1

    op = f.operator()
    if op is operator.gt or op is operator.ge:
        return symbolic_expression(f.rhs() - f.lhs())
    return symbolic_expression(f.lhs() - f.rhs())
