# sage_setup: distribution = sagemath-plot
# sage.doctest: needs sage.symbolic
r"""
2D plotting

Sage provides extensive 2D plotting functionality. The underlying
rendering is done using the matplotlib Python library.

The following graphics primitives are supported:

-  :func:`~sage.plot.arrow.arrow` -- an arrow from a min point to a max point

-  :func:`~sage.plot.circle.circle` -- a circle with given radius

-  :func:`~sage.plot.ellipse.ellipse` -- an ellipse with given radii
   and angle

-  :func:`~sage.plot.arc.arc` -- an arc of a circle or an ellipse

-  :func:`~sage.plot.disk.disk` -- a filled disk (i.e. a sector or wedge of a circle)

-  :func:`~sage.plot.line.line` -- a line determined by a sequence of points (this need not
   be straight!)

-  :func:`~sage.plot.point.point` -- a point

-  :func:`~sage.plot.text.text` -- some text

-  :func:`~sage.plot.polygon.polygon` -- a filled polygon

The following plotting functions are supported:

-  :func:`plot` -- plot of a function or other Sage object (e.g., elliptic
   curve)

-  :func:`parametric_plot`

-  :func:`~sage.plot.contour_plot.implicit_plot`

-  :func:`polar_plot`

-  :func:`~sage.plot.contour_plot.region_plot`

-  :func:`list_plot`

-  :func:`~sage.plot.scatter_plot.scatter_plot`

-  :func:`~sage.plot.bar_chart.bar_chart`

-  :func:`~sage.plot.contour_plot.contour_plot`

-  :func:`~sage.plot.density_plot.density_plot`

-  :func:`~sage.plot.plot_field.plot_vector_field`

-  :func:`~sage.plot.plot_field.plot_slope_field`

-  :func:`~sage.plot.matrix_plot.matrix_plot`

-  :func:`~sage.plot.complex_plot.complex_plot`

-  :func:`graphics_array`

-  :func:`multi_graphics`

-  The following log plotting functions:

   - :func:`plot_loglog`

   - :func:`plot_semilogx` and :func:`plot_semilogy`

   - :func:`list_plot_loglog`

   - :func:`list_plot_semilogx` and :func:`list_plot_semilogy`

The following miscellaneous Graphics functions are included:

-  :func:`~sage.plot.graphics.Graphics`

-  :func:`~sage.plot.graphics.is_Graphics`

-  :func:`~sage.plot.colors.hue`

Type ``?`` after each primitive in Sage for help and examples.

EXAMPLES:

We draw a curve::

    sage: plot(x^2, (x,0,5))
    Graphics object consisting of 1 graphics primitive

.. PLOT::

    g = plot(x**2, (x,0,5))
    sphinx_plot(g)

We draw a circle and a curve::

    sage: circle((1,1), 1) + plot(x^2, (x,0,5))
    Graphics object consisting of 2 graphics primitives

.. PLOT::

    g = circle((1,1), 1) + plot(x**2, (x,0,5))
    sphinx_plot(g)

Notice that the aspect ratio of the above plot makes the plot very tall
because the plot adopts the default aspect ratio of the circle (to make
the circle appear like a circle).  We can change the aspect ratio to be
what we normally expect for a plot by explicitly asking for an
'automatic' aspect ratio::

    sage: show(circle((1,1), 1) + plot(x^2, (x,0,5)), aspect_ratio='automatic')

The aspect ratio describes the apparently height/width ratio of a unit
square.  If you want the vertical units to be twice as big as the
horizontal units, specify an aspect ratio of 2::

    sage: show(circle((1,1), 1) + plot(x^2, (x,0,5)), aspect_ratio=2)

The ``figsize`` option adjusts the figure size.  The default figsize is
4.  To make a figure that is roughly twice as big, use ``figsize=8``::

    sage: show(circle((1,1), 1) + plot(x^2, (x,0,5)), figsize=8)

You can also give separate horizontal and vertical dimensions.  Both
will be measured in inches::

    sage: show(circle((1,1), 1) + plot(x^2, (x,0,5)), figsize=[4,8])

However, do not make the figsize too big (e.g. one dimension greater
than 327 or both in the mid-200s) as this will lead to errors or crashes.
See :meth:`~sage.plot.graphics.Graphics.show` for full details.

Note that the axes will not cross if the data is not on both sides of
both axes, even if it is quite close::

    sage: plot(x^3, (x,1,10))
    Graphics object consisting of 1 graphics primitive

.. PLOT::

    g = plot(x**3, (x,1,10))
    sphinx_plot(g)

When the labels have quite different orders of magnitude or are very
large, scientific notation (the `e` notation for powers of ten) is used::

    sage: plot(x^2, (x,480,500))  # no scientific notation
    Graphics object consisting of 1 graphics primitive

.. PLOT::

    g = plot(x**2, (x,480,500))
    sphinx_plot(g)

::

    sage: plot(x^2, (x,300,500))  # scientific notation on y-axis
    Graphics object consisting of 1 graphics primitive

.. PLOT::

    g = plot(x**2, (x,300,500))
    sphinx_plot(g)

But you can fix your own tick labels, if you know what to expect and
have a preference::

    sage: plot(x^2, (x,300,500), ticks=[100,50000])
    Graphics object consisting of 1 graphics primitive

.. PLOT::

    g = plot(x**2, (x,300,500), ticks=[100,50000])
    sphinx_plot(g)

To change the ticks on one axis only, use the following notation::

    sage: plot(x^2, (x,300,500), ticks=[None,50000])
    Graphics object consisting of 1 graphics primitive

.. PLOT::

    g = plot(x**2, (x,300,500), ticks=[None,50000])
    sphinx_plot(g)

You can even have custom tick labels along with custom positioning. ::

    sage: plot(x^2, (x,0,3), ticks=[[1,2.5],pi/2], tick_formatter=[["$x_1$","$x_2$"],pi])  # long time
    Graphics object consisting of 1 graphics primitive

.. PLOT::

    g = plot(x**2, (x,0,3), ticks=[[1,2.5],pi/2], tick_formatter=[["$x_1$","$x_2$"],pi])
    sphinx_plot(g)

We construct a plot involving several graphics objects::

    sage: G = plot(cos(x), (x, -5, 5), thickness=5, color='green', title='A plot')
    sage: P = polygon([[1,2], [5,6], [5,0]], color='red')
    sage: G + P
    Graphics object consisting of 2 graphics primitives

.. PLOT::

    G = plot(cos(x), (x, -5, 5), thickness=5, color='green', title='A plot')
    P = polygon([[1,2], [5,6], [5,0]], color='red')
    sphinx_plot(G + P)

Next we construct the reflection of the above polygon about the
`y`-axis by iterating over the list of first-coordinates of
the first graphic element of ``P`` (which is the actual
Polygon; note that ``P`` is a Graphics object, which consists
of a single polygon)::

    sage: Q = polygon([(-x,y) for x,y in P[0]], color='blue')
    sage: Q   # show it
    Graphics object consisting of 1 graphics primitive

.. PLOT::

    P = polygon([[1,2], [5,6], [5,0]], color='red')
    Q = polygon([(-x,y) for x,y in P[0]], color='blue')
    sphinx_plot(Q)

We combine together different graphics objects using "+"::

    sage: H = G + P + Q
    sage: print(H)
    Graphics object consisting of 3 graphics primitives
    sage: type(H)
    <class 'sage.plot.graphics.Graphics'>
    sage: H[1]
    Polygon defined by 3 points
    sage: list(H[1])
    [(1.0, 2.0), (5.0, 6.0), (5.0, 0.0)]
    sage: H       # show it
    Graphics object consisting of 3 graphics primitives

.. PLOT::

    G = plot(cos(x), (x, -5, 5), thickness=5, color='green', title='A plot')
    P = polygon([[1,2], [5,6], [5,0]], color='red')
    Q = polygon([(-x,y) for x,y in P[0]], color='blue')
    H = G + P + Q
    sphinx_plot(H)

We can put text in a graph::

    sage: L = [[cos(pi*i/100)^3,sin(pi*i/100)] for i in range(200)]
    sage: p = line(L, rgbcolor=(1/4,1/8,3/4))
    sage: tt = text('A Bulb', (1.5, 0.25))
    sage: tx = text('x axis', (1.5,-0.2))
    sage: ty = text('y axis', (0.4,0.9))
    sage: g = p + tt + tx + ty
    sage: g.show(xmin=-1.5, xmax=2, ymin=-1, ymax=1)

.. PLOT::

    L = [[cos(pi*i/100)**3,sin(pi*i/100)] for i in range(200)]
    p = line(L, rgbcolor=(1.0/4.0,1.0/8.0,3.0/4.0))
    t = text('A Bulb', (1.5, 0.25))
    x = text('x axis', (1.5,-0.2))
    y = text('y axis', (0.4,0.9))
    g = p+t+x+y
    g.xmin(-1.5)
    g.xmax(2)
    g.ymin(-1)
    g.ymax(1)
    sphinx_plot(g)

We can add a graphics object to another one as an inset::

    sage: g1 = plot(x^2*sin(1/x), (x, -2, 2), axes_labels=['$x$', '$y$'])
    sage: g2 = plot(x^2*sin(1/x), (x, -0.3, 0.3), axes_labels=['$x$', '$y$'],
    ....:           frame=True)
    sage: g1.inset(g2, pos=(0.15, 0.7, 0.25, 0.25))
    Multigraphics with 2 elements

.. PLOT::

    g1 = plot(x**2*sin(1/x), (x, -2, 2), axes_labels=['$x$', '$y$'])
    g2 = plot(x**2*sin(1/x), (x, -0.3, 0.3), axes_labels=['$x$', '$y$'], \
              frame=True)
    sphinx_plot(g1.inset(g2, pos=(0.15, 0.7, 0.25, 0.25)))

We can add a title to a graph::

    sage: plot(x^2, (x,-2,2), title='A plot of $x^2$')
    Graphics object consisting of 1 graphics primitive

.. PLOT::

    g=plot(x**2, (x,-2,2), title='A plot of $x^2$')
    sphinx_plot(g)

We can set the position of the title::

    sage: plot(x^2, (-2,2), title='Plot of $x^2$', title_pos=(0.5,-0.05))
    Graphics object consisting of 1 graphics primitive

.. PLOT::

    g=plot(x**2, (-2,2), title='Plot of $x^2$', title_pos=(0.5,-0.05))
    sphinx_plot(g)

We plot the Riemann zeta function along the critical line and see
the first few zeros::

    sage: i = CDF.0      # define i this way for maximum speed.
    sage: p1 = plot(lambda t: arg(zeta(0.5+t*i)), 1, 27, rgbcolor=(0.8,0,0))
    sage: p2 = plot(lambda t: abs(zeta(0.5+t*i)), 1, 27, color=hue(0.7))
    sage: print(p1 + p2)
    Graphics object consisting of 2 graphics primitives
    sage: p1 + p2    # display it
    Graphics object consisting of 2 graphics primitives

.. PLOT::

    from sage.rings.complex_double import ComplexDoubleElement
    i = ComplexDoubleElement(0,1)      # define i this way for maximum speed.
    p1 = plot(lambda t: arg(zeta(0.5+t*i)), 1, 27, rgbcolor=(0.8,0,0))
    p2 = plot(lambda t: abs(zeta(0.5+t*i)), 1, 27, color=hue(0.7))
    g = p1 + p2
    sphinx_plot(g)

.. NOTE::

    Not all functions in Sage are symbolic. When plotting non-symbolic functions
    they should be wrapped in ``lambda``::

        sage: plot(lambda x:fibonacci(round(x)), (x,1,10))
        Graphics object consisting of 1 graphics primitive

.. PLOT::

    g=plot(lambda x:fibonacci(round(x)), (x,1,10))
    sphinx_plot(g)

Many concentric circles shrinking toward the origin::

    sage: show(sum(circle((i,0), i, hue=sin(i/10)) for i in [10,9.9,..,0]))  # long time

.. PLOT::

    g = sum(circle((i,0), i, hue=sin(i/10)) for i in srange(0,10,0.1))
    sphinx_plot(g)

Here is a pretty graph::

    sage: g = Graphics()
    sage: for i in range(60):
    ....:    p = polygon([(i*cos(i),i*sin(i)), (0,i), (i,0)],\
    ....:                color=hue(i/40+0.4), alpha=0.2)
    ....:    g = g + p
    sage: g.show(dpi=200, axes=False)

.. PLOT::

    g=Graphics()
    for i in range(60):
        # i/40 doesn't convert to real number
        p = polygon([(i*cos(i),i*sin(i)), (0,i), (i,0)],\
                    color=hue(0.025*i+0.4), alpha=0.2)
        g = g + p
    g.axes(False)
    sphinx_plot(g)

Another graph::

    sage: x = var('x')
    sage: P = plot(sin(x)/x, -4, 4, color='blue') + \
    ....:     plot(x*cos(x), -4, 4, color='red') + \
    ....:     plot(tan(x), -4, 4, color='green')
    sage: P.show(ymin=-pi, ymax=pi)

.. PLOT::

    g = plot(sin(x)/x, -4, 4, color='blue') + \
        plot(x*cos(x), -4, 4, color='red') + \
        plot(tan(x), -4, 4, color='green')
    g.ymin(-pi)
    g.ymax(pi)
    sphinx_plot(g)

PYX EXAMPLES: These are some examples of plots similar to some of
the plots in the PyX (http://pyx.sourceforge.net) documentation:

Symbolline::

    sage: y(x) = x*sin(x^2)
    sage: v = [(x, y(x)) for x in [-3,-2.95,..,3]]
    sage: show(points(v, rgbcolor=(0.2,0.6, 0.1), pointsize=30) + plot(spline(v), -3.1, 3))

.. PLOT::

    #y(x)=x*sin(x**2) gave SyntaxError: can't assign to function call
    def y(x): return x*sin(x**2)
    v=list()
    for x in srange(-3,3,0.05):
        v.append((x, y(x)))
    g = points(v, rgbcolor=(0.2,0.6, 0.1), pointsize=30) + plot(spline(v), -3.1, 3)
    sphinx_plot(g)

Cycliclink::

    sage: g1 = plot(cos(20*x)*exp(-2*x), 0, 1)
    sage: g2 = plot(2*exp(-30*x) - exp(-3*x), 0, 1)
    sage: show(graphics_array([g1, g2], 2, 1))

.. PLOT::

    g1 = plot(cos(20*x)*exp(-2*x), 0, 1)
    g2 = plot(2*exp(-30*x) - exp(-3*x), 0, 1)
    g = graphics_array([g1, g2], 2, 1)
    sphinx_plot(g)

Pi Axis::

    sage: g1 = plot(sin(x), 0, 2*pi)
    sage: g2 = plot(cos(x), 0, 2*pi, linestyle='--')
    sage: (g1 + g2).show(ticks=pi/6,        # show their sum, nicely formatted  # long time
    ....:                tick_formatter=pi)

.. PLOT::

    g1 = plot(sin(x), 0, 2*pi, ticks=pi/6, tick_formatter=pi)
    g2 = plot(cos(x), 0, 2*pi, linestyle='--', ticks=pi/6, tick_formatter=pi)
    sphinx_plot(g1+g2)

An illustration of integration::

    sage: f(x) = (x-3)*(x-5)*(x-7)+40
    sage: P = line([(2,0),(2,f(2))], color='black')
    sage: P += line([(8,0),(8,f(8))], color='black')
    sage: P += polygon([(2,0),(2,f(2))] + [(x, f(x)) for x in [2,2.1,..,8]] + [(8,0),(2,0)],
    ....:              rgbcolor=(0.8,0.8,0.8), aspect_ratio='automatic')
    sage: P += text("$\\int_{a}^b f(x) dx$", (5, 20), fontsize=16, color='black')
    sage: P += plot(f, (1, 8.5), thickness=3)
    sage: P    # show the result
    Graphics object consisting of 5 graphics primitives

.. PLOT::

    #inline f substitution to avoid SyntaxError: can't assign to function call in sphinx_plot
    def f(x): return (x-3)*(x-5)*(x-7)+40
    P = line([(2,0),(2,f(2))], color='black')
    P = P + line([(8,0),(8,f(8))], color='black')
    L = list(((2,0), (2,f(2))))
    for i in srange(2,8.1,0.1):
        L.append((i,f(i)))
    L.append((8,0))
    L.append((2,0))
    P = P + polygon(L, rgbcolor=(0.8,0.8,0.8), aspect_ratio='automatic')
    P = P + text("$\\int_{a}^b f(x) dx$", (5, 20), fontsize=16, color='black')
    P = P + plot(f, (1, 8.5), thickness=3)
    sphinx_plot(P)

NUMERICAL PLOTTING:

Sage includes Matplotlib, which provides 2D plotting with an interface
that is a likely very familiar to people doing numerical
computation.
You can use ``plt.clf()`` to clear the current image frame
and ``plt.close()`` to close it.
For example,

::

    sage: import pylab as plt
    sage: t = plt.arange(0.0, 2.0, 0.01)
    sage: s = sin(2*pi*t)
    sage: P = plt.plot(t, s, linewidth=1.0)
    sage: xl = plt.xlabel('time (s)')
    sage: yl = plt.ylabel('voltage (mV)')
    sage: t = plt.title('About as simple as it gets, folks')
    sage: plt.grid(True)
    sage: import tempfile
    sage: with tempfile.NamedTemporaryFile(suffix='.png') as f1:
    ....:     plt.savefig(f1.name)
    sage: plt.clf()
    sage: with tempfile.NamedTemporaryFile(suffix='.png') as f2:
    ....:     plt.savefig(f2.name)
    sage: plt.close()
    sage: plt.imshow([[1,2],[0,1]])
    <matplotlib.image.AxesImage object at ...>

We test that ``imshow`` works as well, verifying that
:issue:`2900` is fixed (in Matplotlib).

::

    sage: plt.imshow([[(0.0,0.0,0.0)]])
    <matplotlib.image.AxesImage object at ...>
    sage: import tempfile
    sage: with tempfile.NamedTemporaryFile(suffix='.png') as f:
    ....:     plt.savefig(f.name)

Since the above overwrites many Sage plotting functions, we reset
the state of Sage, so that the examples below work!

::

    sage: reset()

See http://matplotlib.sourceforge.net for complete documentation
about how to use Matplotlib.

TESTS:

We test dumping and loading a plot.

::

    sage: p = plot(sin(x), (x, 0,2*pi))
    sage: Q = loads(dumps(p))

Verify that a clean sage startup does *not* import matplotlib::

    sage: os.system("sage -c \"if 'matplotlib' in sys.modules: sys.exit(1)\"")  # long time
    0

Verify that :issue:`10980` is fixed::

    sage: plot(x,0,2,gridlines=([sqrt(2)],[]))
    Graphics object consisting of 1 graphics primitive

AUTHORS:

- Alex Clemesha and William Stein (2006-04-10): initial version

- David Joyner: examples

- Alex Clemesha (2006-05-04) major update

- William Stein (2006-05-29): fine tuning, bug fixes, better server
  integration

- William Stein (2006-07-01): misc polish

- Alex Clemesha (2006-09-29): added contour_plot, frame axes, misc
  polishing

- Robert Miller (2006-10-30): tuning, NetworkX primitive

- Alex Clemesha (2006-11-25): added plot_vector_field, matrix_plot,
  arrow, bar_chart, Axes class usage (see axes.py)

- Bobby Moretti and William Stein (2008-01): Change plot to specify
  ranges using the (varname, min, max) notation.

- William Stein (2008-01-19): raised the documentation coverage from a
  miserable 12 percent to a 'wopping' 35 percent, and fixed and
  clarified numerous small issues.

- Jason Grout (2009-09-05): shifted axes and grid functionality over
  to matplotlib; fixed a number of smaller issues.

- Jason Grout (2010-10): rewrote aspect ratio portions of the code

- Jeroen Demeyer (2012-04-19): move parts of this file to graphics.py (:issue:`12857`)

- Aaron Lauve (2016-07-13): reworked handling of 'color' when passed
  a list of functions; now more in-line with other CAS's. Added list functionality
  to linestyle and legend_label options as well. (:issue:`12962`)

- Eric Gourgoulhon (2019-04-24): add :func:`multi_graphics` and insets
"""
# ****************************************************************************
#       Copyright (C) 2006 Alex Clemesha <clemesha@gmail.com>
#       Copyright (C) 2006-2008 William Stein <wstein@gmail.com>
#       Copyright (C) 2010 Jason Grout
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from functools import reduce

# IMPORTANT: Do *not* import matplotlib at module scope.  It takes a
# surprisingly long time to initialize itself.  It's better if it is
# imported in functions, so it only gets started if it is actually
# going to be used.

# DEFAULT_FIGSIZE=(6, 3.70820393249937)
import sage.misc.verbose
from sage.arith.srange import srange

from sage.misc.randstate import current_randstate  # for plot adaptive refinement
from math import sin, cos, pi, log, exp  # for polar_plot and log scaling

from sage.ext.fast_eval import fast_float, is_fast_float
from sage.structure.element import Expression
from sage.misc.decorators import options

from .graphics import Graphics
from .multigraphics import GraphicsArray, MultiGraphics
from sage.plot.polygon import polygon

# import of line2d below is only for redirection of imports
from sage.plot.line import line
from sage.misc.lazy_import import lazy_import
lazy_import('sage.plot.line', 'line2d', deprecation=28717)

# Currently not used - see comment immediately above about
# figure.canvas.mpl_connect('draw_event', pad_for_tick_labels)
# TODO - figure out how to use this, add documentation
# def pad_for_tick_labels(event):
#    import matplotlib.transforms as mtransforms
#    figure=event.canvas.figure
#    bboxes = []
#    for ax in figure.axes:
#        bbox = ax.xaxis.get_label().get_window_extent()
#        # the figure transform goes from relative coords->pixels and we
#        # want the inverse of that
#        bboxi = bbox.inverse_transformed(figure.transFigure)
#        bboxes.append(bboxi)
#
#        bbox = ax.yaxis.get_label().get_window_extent()
#        bboxi = bbox.inverse_transformed(figure.transFigure)
#        bboxes.append(bboxi)
#        for label in (ax.get_xticklabels()+ax.get_yticklabels() \
#                          + ax.get_xticklabels(minor=True) \
#                          +ax.get_yticklabels(minor=True)):
#            bbox = label.get_window_extent()
#            bboxi = bbox.inverse_transformed(figure.transFigure)
#            bboxes.append(bboxi)
#
#    # this is the bbox that bounds all the bboxes, again in relative
#    # figure coords
#    bbox = mtransforms.Bbox.union(bboxes)
#    adjusted=adjust_figure_to_contain_bbox(figure,bbox)
#
#    if adjusted:
#        figure.canvas.draw()
#    return False
#
# Currently not used - see comment above about
# figure.canvas.mpl_connect('draw_event', pad_for_tick_labels)
# TODO - figure out how to use this, add documentation
# def adjust_figure_to_contain_bbox(fig, bbox,pad=1.1):
#    """
#    For each amount we are over (in axes coordinates), we adjust by over*pad
#    to give ourselves a bit of padding.
#    """
#    left=fig.subplotpars.left
#    bottom=fig.subplotpars.bottom
#    right=fig.subplotpars.right
#    top=fig.subplotpars.top
#
#    adjusted=False
#    if bbox.xmin<0:
#        left-=bbox.xmin*pad
#        adjusted=True
#    if bbox.ymin<0:
#        bottom-=bbox.ymin*pad
#        adjusted=True
#    if bbox.xmax>1:
#        right-=(bbox.xmax-1)*pad
#        adjusted=True
#    if bbox.ymax>1:
#        top-=(bbox.ymax-1)*pad
#        adjusted=True
#
#    if left<right and bottom<top:
#        fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top)
#        return adjusted
#    else:
#        return False

_SelectiveFormatterClass = None


def SelectiveFormatter(formatter, skip_values):
    """
    This matplotlib formatter selectively omits some tick values and
    passes the rest on to a specified formatter.

    EXAMPLES:

    This example is almost straight from a matplotlib example.

    ::

        sage: # needs numpy
        sage: from sage.plot.plot import SelectiveFormatter
        sage: import matplotlib.pyplot as plt
        sage: import numpy
        sage: fig = plt.figure()
        sage: ax = fig.add_subplot(111)
        sage: t = numpy.arange(0.0, 2.0, 0.01)
        sage: s = numpy.sin(2*numpy.pi*t)
        sage: p = ax.plot(t, s)
        sage: formatter = SelectiveFormatter(ax.xaxis.get_major_formatter(),
        ....:                                skip_values=[0,1])
        sage: ax.xaxis.set_major_formatter(formatter)
        sage: import tempfile
        sage: with tempfile.NamedTemporaryFile(suffix='.png') as f:
        ....:     fig.savefig(f.name)
    """
    global _SelectiveFormatterClass
    if _SelectiveFormatterClass is None:

        from matplotlib.ticker import Formatter

        class _SelectiveFormatterClass(Formatter):
            def __init__(self, formatter, skip_values):
                """
                Initialize a SelectiveFormatter object.

                INPUT:

                  - ``formatter`` -- the formatter object to which we should pass labels

                  - ``skip_values`` -- list of values that we should skip when
                    formatting the tick labels

                EXAMPLES::

                    sage: # needs numpy
                    sage: from sage.plot.plot import SelectiveFormatter
                    sage: import matplotlib.pyplot as plt
                    sage: import numpy
                    sage: fig = plt.figure()
                    sage: ax = fig.add_subplot(111)
                    sage: t = numpy.arange(0.0, 2.0, 0.01)
                    sage: s = numpy.sin(2*numpy.pi*t)
                    sage: line = ax.plot(t, s)
                    sage: formatter = SelectiveFormatter(ax.xaxis.get_major_formatter(),
                    ....:                                skip_values=[0,1])
                    sage: ax.xaxis.set_major_formatter(formatter)
                    sage: from tempfile import NamedTemporaryFile
                    sage: with NamedTemporaryFile(suffix='.png') as f:
                    ....:     fig.savefig(f.name)
                """
                self.formatter = formatter
                self.skip_values = skip_values

            def set_locs(self, locs):
                """
                Set the locations for the ticks that are not skipped.

                EXAMPLES::

                    sage: from sage.plot.plot import SelectiveFormatter
                    sage: import matplotlib.ticker
                    sage: formatter = SelectiveFormatter(matplotlib.ticker.Formatter(),
                    ....:                                skip_values=[0,200])
                    sage: formatter.set_locs([i*100 for i in range(10)])
                """
                self.formatter.set_locs([l for l in locs if l not in self.skip_values])

            def __call__(self, x, *args, **kwds):
                """
                Return the format for tick val *x* at position *pos*.

                EXAMPLES::

                    sage: from sage.plot.plot import SelectiveFormatter
                    sage: import matplotlib.ticker
                    sage: formatter = SelectiveFormatter(matplotlib.ticker.FixedFormatter(['a','b']),
                    ....:                                skip_values=[0,2])
                    sage: [formatter(i,1) for i in range(10)]
                    ['', 'b', '', 'b', 'b', 'b', 'b', 'b', 'b', 'b']
                """
                if x in self.skip_values:
                    return ''
                else:
                    return self.formatter(x, *args, **kwds)

    return _SelectiveFormatterClass(formatter, skip_values)


def xydata_from_point_list(points):
    r"""
    Return two lists (xdata, ydata), each coerced to a list of floats,
    which correspond to the x-coordinates and the y-coordinates of the
    points.

    The points parameter can be a list of 2-tuples or some object that
    yields a list of one or two numbers.

    This function can potentially be very slow for large point sets.

    TESTS::

        sage: from sage.plot.plot import xydata_from_point_list
        sage: xydata_from_point_list([CC(0), CC(1)])   # issue 8082
        ([0.0, 1.0], [0.0, 0.0])

    This function should work for anything than can be turned into a
    list, such as iterators and such (see :issue:`10478`)::

        sage: xydata_from_point_list(iter([(0,0), (sqrt(3), 2)]))
        ([0.0, 1.7320508075688772], [0.0, 2.0])
        sage: xydata_from_point_list((x, x^2) for x in range(5))
        ([0.0, 1.0, 2.0, 3.0, 4.0], [0.0, 1.0, 4.0, 9.0, 16.0])
        sage: xydata_from_point_list(enumerate(prime_range(1, 15)))
        ([0.0, 1.0, 2.0, 3.0, 4.0, 5.0], [2.0, 3.0, 5.0, 7.0, 11.0, 13.0])
        sage: from builtins import zip
        sage: xydata_from_point_list(list(zip([2,3,5,7], [11, 13, 17, 19])))
        ([2.0, 3.0, 5.0, 7.0], [11.0, 13.0, 17.0, 19.0])

    The code now accepts mixed lists of complex and real numbers::

        sage: xydata_from_point_list(map(N,[0,1,1+I,I,I-1,-1,-1-I,-I,1-I]))
        ([0.0, 1.0, 1.0, 0.0, -1.0, -1.0, -1.0, 0.0, 1.0],
        [0.0, 0.0, 1.0, 1.0, 1.0, 0.0, -1.0, -1.0, -1.0])
        sage: point2d([0, 1., CC(0,1)])
        Graphics object consisting of 1 graphics primitive
        sage: point2d((x^5-1).roots(multiplicities=False))
        Graphics object consisting of 1 graphics primitive
    """
    import numbers
    zero = float(0)

    xdata = []
    ydata = []
    for xy in points:
        if isinstance(xy, Expression):
            xy = xy.n()
        if isinstance(xy, numbers.Real):
            xdata.append(float(xy))
            ydata.append(zero)
        elif isinstance(xy, numbers.Complex):
            xdata.append(float(xy.real()))
            ydata.append(float(xy.imag()))
        else:
            try:
                x, y = xy
            except TypeError:
                raise TypeError('invalid input for 2D point')
            else:
                xdata.append(float(x))
                ydata.append(float(y))
    return xdata, ydata


@options(alpha=1, thickness=1, fill=False, fillcolor='automatic',
         fillalpha=0.5, plot_points=200, adaptive_tolerance=0.01,
         adaptive_recursion=5, detect_poles=False, exclude=None,
         legend_label=None, __original_opts=True,
         aspect_ratio='automatic', imaginary_tolerance=1e-8)
def plot(funcs, *args, **kwds):
    r"""
    Use plot by writing.

    ``plot(X, ...)``

    where `X` is a Sage object (or list of Sage objects) that
    either is callable and returns numbers that can be coerced to
    floats, or has a plot method that returns a
    ``GraphicPrimitive`` object.

    There are many other specialized 2D plot commands available
    in Sage, such as ``plot_slope_field``, as well as various
    graphics primitives like :class:`~sage.plot.arrow.Arrow`;
    type ``sage.plot.plot?`` for a current list.

    Type ``plot.options`` for a dictionary of the default
    options for plots. You can change this to change the defaults for
    all future plots. Use ``plot.reset()`` to reset to the
    default options.

    PLOT OPTIONS:

    - ``plot_points`` -- (default: 200) the minimal number of plot points

    - ``adaptive_recursion`` -- (default: 5) how many levels of recursion to go
      before giving up when doing adaptive refinement.  Setting this to 0
      disables adaptive refinement.

    - ``adaptive_tolerance`` -- (default: 0.01) how large a difference should be
      before the adaptive refinement code considers it significant.  See the
      documentation further below for more information, starting at "the
      algorithm used to insert".

    - ``imaginary_tolerance`` -- (default: ``1e-8``) if an imaginary
      number arises (due, for example, to numerical issues), this
      tolerance specifies how large it has to be in magnitude before
      we raise an error.  In other words, imaginary parts smaller than
      this are ignored in your plot points.

    - ``base`` -- (default: `10`) the base of the logarithm if
      a logarithmic scale is set. This must be greater than 1. The base
      can be also given as a list or tuple ``(basex, basey)``.
      ``basex`` sets the base of the logarithm along the horizontal
      axis and ``basey`` sets the base along the vertical axis.

    - ``scale`` -- string (default: ``'linear'``); scale of the axes.
      Possible values are ``'linear'``, ``'loglog'``, ``'semilogx'``,
      ``'semilogy'``.

      The scale can be also be given as single argument that is a list
      or tuple ``(scale, base)`` or ``(scale, basex, basey)``.

      The ``'loglog'`` scale sets both the horizontal and vertical axes to
      logarithmic scale. The ``'semilogx'`` scale sets the horizontal axis
      to logarithmic scale. The ``'semilogy'`` scale sets the vertical axis
      to logarithmic scale. The ``'linear'`` scale is the default value
      when :class:`~sage.plot.graphics.Graphics` is initialized.

    - ``xmin`` -- starting x value in the rendered figure. This parameter is
      passed directly to the ``show`` procedure and it could be overwritten.

    - ``xmax`` -- ending x value in the rendered figure. This parameter is passed
      directly to the ``show`` procedure and it could be overwritten.

    - ``ymin`` -- starting y value in the rendered figure. This parameter is
      passed directly to the ``show`` procedure and it could be overwritten.

    - ``ymax`` -- ending y value in the rendered figure. This parameter is passed
      directly to the ``show`` procedure and it could be overwritten.

    - ``detect_poles`` -- boolean (default: ``False``); if set to ``True`` poles are detected.
      If set to "show" vertical asymptotes are drawn.

    - ``legend_label`` -- a (TeX) string serving as the label for `X` in the legend.
      If `X` is a list, then this option can be a single string, or a list or dictionary
      with strings as entries/values. If a dictionary, then keys are taken from ``range(len(X))``.

    .. NOTE::

        - If the ``scale`` is ``'linear'``, then irrespective of what
          ``base`` is set to, it will default to 10 and will remain unused.

        - If you want to limit the plot along the horizontal axis in the
          final rendered figure, then pass the ``xmin`` and ``xmax``
          keywords to the :meth:`~sage.plot.graphics.Graphics.show` method.
          To limit the plot along the vertical axis, ``ymin`` and ``ymax``
          keywords can be provided to either this ``plot`` command or to
          the ``show`` command.

        - This function does NOT simply sample equally spaced points
          between xmin and xmax. Instead it computes equally spaced points
          and adds small perturbations to them. This reduces the possibility
          of, e.g., sampling `\sin` only at multiples of `2\pi`, which would
          yield a very misleading graph.

        - If there is a range of consecutive points where the function has
          no value, then those points will be excluded from the plot. See
          the example below on automatic exclusion of points.

        - For the other keyword options that the ``plot`` function can
          take, refer to the method :meth:`~sage.plot.graphics.Graphics.show`
          and the further options below.

    COLOR OPTIONS:

    - ``color`` -- (default: ``'blue'``) one of:

      - an RGB tuple (r,g,b) with each of r,g,b between 0 and 1.

      - a color name as a string (e.g., ``'purple'``).

      - an HTML color such as '#aaff0b'.

      - a list or dictionary of colors (valid only if `X` is a list):
        if a dictionary, keys are taken from ``range(len(X))``;
        the entries/values of the list/dictionary may be any of the options above.

      - ``'automatic'`` -- maps to default ('blue') if `X` is a single Sage object; and
        maps to a fixed sequence of regularly spaced colors if `X` is a list

    - ``legend_color`` -- the color of the text for `X` (or each item in `X`) in the legend.
      Default color is 'black'. Options are as in ``color`` above, except that the choice 'automatic' maps to 'black' if `X` is a single Sage object

    - ``fillcolor`` -- the color of the fill for the plot of `X` (or each item in `X`).
      Default color is 'gray' if `X` is a single Sage object or if ``color`` is a single color. Otherwise, options are as in ``color`` above

    APPEARANCE OPTIONS:

    The following options affect the appearance of
    the line through the points on the graph of `X` (these are
    the same as for the line function):

    INPUT:

    - ``alpha`` -- how transparent the line is

    - ``thickness`` -- how thick the line is

    - ``rgbcolor`` -- the color as an RGB tuple

    - ``hue`` -- the color given as a hue

    LINE OPTIONS:

    Any MATPLOTLIB line option may also be passed in.  E.g.,

    - ``linestyle`` -- (default: ``'-'``) the style of the line, which is one of

      - ``'-'`` or ``'solid'``
      - ``'--'`` or ``'dashed'``
      - ``'-.'`` or ``'dash dot'``
      - ``':'`` or ``'dotted'``
      - ``"None"`` or ``" "`` or ``""`` (nothing)
      - a list or dictionary (see below)

      The linestyle can also be prefixed with a drawing style (e.g., ``'steps--'``)

      - ``'default'`` (connect the points with straight lines)
      - ``'steps'`` or ``'steps-pre'`` (step function; horizontal
        line is to the left of point)
      - ``'steps-mid'`` (step function; points are in the middle of
        horizontal lines)
      - ``'steps-post'`` (step function; horizontal line is to the
        right of point)

      If `X` is a list, then ``linestyle`` may be a list (with entries
      taken from the strings above) or a dictionary (with keys in ``range(len(X))``
      and values taken from the strings above).

    - ``marker`` -- the style of the markers, which is one of

      - ``"None"`` or ``" "`` or ``""`` (nothing) -- default
      - ``","`` (pixel), ``"."`` (point)
      - ``"_"`` (horizontal line), ``"|"`` (vertical line)
      - ``"o"`` (circle), ``"p"`` (pentagon), ``"s"`` (square), ``"x"`` (x), ``"+"`` (plus), ``"*"`` (star)
      - ``"D"`` (diamond), ``"d"`` (thin diamond)
      - ``"H"`` (hexagon), ``"h"`` (alternative hexagon)
      - ``"<"`` (triangle left), ``">"`` (triangle right), ``"^"`` (triangle up), ``"v"`` (triangle down)
      - ``"1"`` (tri down), ``"2"`` (tri up), ``"3"`` (tri left), ``"4"`` (tri right)
      - ``0`` (tick left), ``1`` (tick right), ``2`` (tick up), ``3`` (tick down)
      - ``4`` (caret left), ``5`` (caret right), ``6`` (caret up), ``7`` (caret down), ``8`` (octagon)
      - ``"$...$"`` (math TeX string)
      - ``(numsides, style, angle)`` to create a custom, regular symbol

        - ``numsides`` -- the number of sides
        - ``style`` -- ``0`` (regular polygon), ``1`` (star shape), ``2`` (asterisk), ``3`` (circle)
        - ``angle`` -- the angular rotation in degrees

    - ``markersize`` -- the size of the marker in points

    - ``markeredgecolor`` -- the color of the marker edge

    - ``markerfacecolor`` -- the color of the marker face

    - ``markeredgewidth`` -- the size of the marker edge in points

    - ``exclude`` -- (default: ``None``) values which are excluded from the plot range.
      Either a list of real numbers, or an equation in one variable.

    FILLING OPTIONS:

    - ``fill`` -- boolean (default: ``False``); one of:

      - "axis" or ``True``: Fill the area between the function and the x-axis.

      - "min": Fill the area between the function and its minimal value.

      - "max": Fill the area between the function and its maximal value.

      - a number c: Fill the area between the function and the horizontal line y = c.

      - a function g: Fill the area between the function that is plotted and g.

      - a dictionary ``d`` (only if a list of functions are plotted):
        The keys of the dictionary should be integers.
        The value of ``d[i]`` specifies the fill options for the i-th function
        in the list.  If ``d[i] == [j]``: Fill the area between the i-th and
        the j-th function in the list.  (But if ``d[i] == j``: Fill the area
        between the i-th function in the list and the horizontal line y = j.)

    - ``fillalpha`` -- (default: 0.5) how transparent the fill is;
      a number between 0 and 1

    MATPLOTLIB STYLE SHEET OPTION:

    - ``stylesheet`` -- (default: classic) support for loading a full matplotlib style sheet.
      Any style sheet listed in ``matplotlib.pyplot.style.available`` is acceptable. If a
      non-existing style is provided the default classic is applied.

    EXAMPLES:

    We plot the `\sin` function::

        sage: P = plot(sin, (0,10)); print(P)
        Graphics object consisting of 1 graphics primitive
        sage: len(P)     # number of graphics primitives
        1
        sage: len(P[0])  # how many points were computed (random)
        225
        sage: P          # render
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        P = plot(sin, (0,10))
        sphinx_plot(P)

    ::

        sage: P = plot(sin, (0,10), plot_points=10); print(P)
        Graphics object consisting of 1 graphics primitive
        sage: len(P[0])  # random output
        32
        sage: P          # render
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        P = plot(sin, (0,10), plot_points=10)
        sphinx_plot(P)

    We plot with ``randomize=False``, which makes the initial sample points
    evenly spaced (hence always the same). Adaptive plotting might
    insert other points, however, unless ``adaptive_recursion=0``.

    ::

        sage: p = plot(1, (x,0,3), plot_points=4, randomize=False, adaptive_recursion=0)
        sage: list(p[0])
        [(0.0, 1.0), (1.0, 1.0), (2.0, 1.0), (3.0, 1.0)]

    Some colored functions::

        sage: plot(sin, 0, 10, color='purple')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        P=plot(sin,0,10,color='purple')
        sphinx_plot(P)

    ::

        sage: plot(sin, 0, 10, color='#ff00ff')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        P=plot(sin, 0, 10, color='#ff00ff')
        sphinx_plot(P)

    We plot several functions together by passing a list of functions
    as input::

        sage: plot([x*exp(-n*x^2)/.4 for n in [1..5]], (0, 2), aspect_ratio=.8)
        Graphics object consisting of 5 graphics primitives

    .. PLOT::

        g = plot([x*exp(-n*x**2)/.4 for n in range(1,6)], (0, 2), aspect_ratio=.8)
        sphinx_plot(g)

    By default, color will change from one primitive to the next.
    This may be controlled by modifying ``color`` option::

        sage: g1 = plot([x*exp(-n*x^2)/.4 for n in [1..3]], (0, 2),
        ....:           color='blue', aspect_ratio=.8); g1
        Graphics object consisting of 3 graphics primitives
        sage: g2 = plot([x*exp(-n*x^2)/.4 for n in [1..3]], (0, 2),
        ....:           color=['red','red','green'], linestyle=['-','--','-.'],
        ....:           aspect_ratio=.8); g2
        Graphics object consisting of 3 graphics primitives

    .. PLOT::

        g1 = plot([x*exp(-n*x**2)/.4 for n in range(1,4)], (0, 2), color='blue', aspect_ratio=.8)
        g2 = plot([x*exp(-n*x**2)/.4 for n in range(1,4)], (0, 2), color=['red','red','green'], linestyle=['-','--','-.'], aspect_ratio=.8)
        sphinx_plot(graphics_array([[g1], [g2]]))

    While plotting real functions, imaginary numbers that are "almost
    real" will inevitably arise due to numerical issues. By tweaking
    the ``imaginary_tolerance``, you can decide how large of an
    imaginary part you're willing to sweep under the rug in order to
    plot the corresponding point. If a particular value's imaginary
    part has magnitude larger than ``imaginary_tolerance``, that point
    will not be plotted. The default tolerance is ``1e-8``, so the
    imaginary part in the first example below is ignored, but the
    second example "fails," emits a warning, and produces an empty
    graph::

        sage: f = x + I*1e-12
        sage: plot(f, x, -1, 1)
        Graphics object consisting of 1 graphics primitive
        sage: plot(f, x, -1, 1, imaginary_tolerance=0)
        ...WARNING: ...Unable to compute ...
        Graphics object consisting of 0 graphics primitives

    We can also build a plot step by step from an empty plot::

        sage: a = plot([]); a       # passing an empty list returns an empty plot (Graphics() object)
        Graphics object consisting of 0 graphics primitives
        sage: a += plot(x**2); a    # append another plot
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        a = plot([])
        a = a + plot(x**2)
        sphinx_plot(a)

    ::

        sage: a += plot(x**3); a    # append yet another plot
        Graphics object consisting of 2 graphics primitives

    .. PLOT::

        a = plot([])
        a = a + plot(x**2)
        a = a + plot(x**3)
        sphinx_plot(a)

    The function `\sin(1/x)` wiggles wildly near `0`.
    Sage adapts to this and plots extra points near the origin.

    ::

        sage: plot(sin(1/x), (x, -1, 1))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = plot(sin(1/x), (x, -1, 1))
        sphinx_plot(g)

    Via the matplotlib library, Sage makes it easy to tell whether
    a graph is on both sides of both axes, as the axes only cross
    if the origin is actually part of the viewing area::

        sage: plot(x^3, (x,0,2))  # this one has the origin
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = plot(x**3,(x,0,2))
        sphinx_plot(g)

    ::

        sage: plot(x^3, (x,1,2))  # this one does not
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = plot(x**3,(x,1,2))
        sphinx_plot(g)

    Another thing to be aware of with axis labeling is that when
    the labels have quite different orders of magnitude or are very
    large, scientific notation (the `e` notation for powers of ten) is used::

        sage: plot(x^2, (x,480,500))  # this one has no scientific notation
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = plot(x**2,(x,480,500))
        sphinx_plot(g)

    ::

        sage: plot(x^2, (x,300,500))  # this one has scientific notation on y-axis
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = plot(x**2,(x,300,500))
        sphinx_plot(g)

    You can put a legend with ``legend_label`` (the legend is only put
    once in the case of multiple functions)::

        sage: plot(exp(x), 0, 2, legend_label='$e^x$')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = plot(exp(x), 0, 2, legend_label='$e^x$')
        sphinx_plot(g)

    Sage understands TeX, so these all are slightly different, and you can choose
    one based on your needs::

        sage: plot(sin, legend_label='sin')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = plot(sin,legend_label='sin')
        sphinx_plot(g)

    ::

        sage: plot(sin, legend_label='$sin$')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = plot(sin,legend_label='$sin$')
        sphinx_plot(g)

    ::

        sage: plot(sin, legend_label=r'$\sin$')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = plot(sin,legend_label=r'$\sin$')
        sphinx_plot(g)

    It is possible to use a different color for the text of each label::

        sage: p1 = plot(sin, legend_label='sin', legend_color='red')
        sage: p2 = plot(cos, legend_label='cos', legend_color='green')
        sage: p1 + p2
        Graphics object consisting of 2 graphics primitives

    .. PLOT::

        p1 = plot(sin, legend_label='sin', legend_color='red')
        p2 = plot(cos, legend_label='cos', legend_color='green')
        g = p1 + p2
        sphinx_plot(g)

    Prior to :issue:`19485`, legends by default had a shadowless gray
    background. This behavior can be recovered by setting the legend
    options on your plot object::

        sage: p = plot(sin(x), legend_label=r'$\sin(x)$')
        sage: p.set_legend_options(back_color=(0.9,0.9,0.9), shadow=False)

    .. PLOT::

        g = plot(sin(x), legend_label=r'$\sin(x)$')
        g.set_legend_options(back_color=(0.9,0.9,0.9), shadow=False)
        sphinx_plot(g)

    If `X` is a list of Sage objects and ``legend_label`` is 'automatic', then Sage will
    create labels for each function according to their internal representation::

        sage: plot([sin(x), tan(x), 1 - x^2], legend_label='automatic')
        Graphics object consisting of 3 graphics primitives

    .. PLOT::

        g = plot([sin(x), tan(x), 1-x**2], legend_label='automatic')
        sphinx_plot(g)

    If ``legend_label`` is any single string other than 'automatic',
    then it is repeated for all members of `X`::

        sage: plot([sin(x), tan(x)], color='blue', legend_label='trig')
        Graphics object consisting of 2 graphics primitives

    .. PLOT::

        g = plot([sin(x), tan(x)], color='blue', legend_label='trig')
        sphinx_plot(g)

    Note that the independent variable may be omitted if there is no
    ambiguity::

        sage: plot(sin(1.0/x), (-1, 1))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = plot(sin(1.0/x), (-1, 1))
        sphinx_plot(g)

    Plotting in logarithmic scale is possible for 2D plots.  There
    are two different syntaxes supported::

        sage: plot(exp, (1, 10), scale='semilogy') # log axis on vertical
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = plot(exp, (1, 10), scale='semilogy')
        sphinx_plot(g)

    ::

        sage: plot_semilogy(exp, (1, 10))  # same thing
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = plot_semilogy(exp, (1, 10))
        sphinx_plot(g)

    ::

        sage: plot_loglog(exp, (1, 10))   # both axes are log
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = plot_loglog(exp, (1, 10))
        sphinx_plot(g)

    ::

        sage: plot(exp, (1, 10), scale='loglog', base=2)  # base of log is 2    # long time
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = plot(exp, (1, 10), scale='loglog', base=2)
        sphinx_plot(g)

    We can also change the scale of the axes in the graphics just before
    displaying::

        sage: G = plot(exp, 1, 10)                                              # long time
        sage: G.show(scale=('semilogy', 2))                                     # long time

    The algorithm used to insert extra points is actually pretty
    simple. On the picture drawn by the lines below::

        sage: p = plot(x^2, (-0.5, 1.4)) + line([(0,0), (1,1)], color='green')
        sage: p += line([(0.5, 0.5), (0.5, 0.5^2)], color='purple')
        sage: p += point(((0, 0), (0.5, 0.5), (0.5, 0.5^2), (1, 1)),
        ....:            color='red', pointsize=20)
        sage: p += text('A', (-0.05, 0.1), color='red')
        sage: p += text('B', (1.01, 1.1), color='red')
        sage: p += text('C', (0.48, 0.57), color='red')
        sage: p += text('D', (0.53, 0.18), color='red')
        sage: p.show(axes=False, xmin=-0.5, xmax=1.4, ymin=0, ymax=2)

    .. PLOT::

        g = plot(x**2, (-0.5, 1.4)) + line([(0,0), (1,1)], color='green')
        g = g + line([(0.5, 0.5), (0.5, 0.5**2)], color='purple')
        g = g + point(((0, 0), (0.5, 0.5), (0.5, 0.5**2), (1, 1)), color='red', pointsize=20)
        g = g + text('A', (-0.05, 0.1), color='red')
        g = g + text('B', (1.01, 1.1), color='red')
        g = g + text('C', (0.48, 0.57), color='red')
        g = g + text('D', (0.53, 0.18), color='red')
        g.axes(False)
        g.xmin(-0.5)
        g.xmax(1.4)
        g.ymin(0)
        g.ymax(2)
        sphinx_plot(g)

    You have the function (in blue) and its approximation (in green)
    passing through the points A and B. The algorithm finds the
    midpoint C of AB and computes the distance between C and D. If that
    distance exceeds the ``adaptive_tolerance`` threshold (*relative* to
    the size of the initial plot subintervals), the point D is
    added to the curve.  If D is added to the curve, then the
    algorithm is applied recursively to the points A and D, and D and
    B. It is repeated ``adaptive_recursion`` times (5, by default).

    The actual sample points are slightly randomized, so the above
    plots may look slightly different each time you draw them.

    We draw the graph of an elliptic curve as the union of graphs of 2
    functions.

    ::

        sage: def h1(x): return abs(sqrt(x^3 - 1))
        sage: def h2(x): return -abs(sqrt(x^3 - 1))
        sage: P = plot([h1, h2], 1,4)
        sage: P          # show the result
        Graphics object consisting of 2 graphics primitives

    .. PLOT::

        def h1(x): return abs(sqrt(x**3  - 1))
        def h2(x): return -abs(sqrt(x**3  - 1))
        P = plot([h1, h2], 1,4)
        sphinx_plot(P)

    It is important to mention that when we draw several graphs at the same time,
    parameters ``xmin``, ``xmax``, ``ymin`` and ``ymax`` are just passed directly
    to the ``show`` procedure. In fact, these parameters would be overwritten::

        sage: p=plot(x^3, x, xmin=-1, xmax=1,ymin=-1, ymax=1)
        sage: q=plot(exp(x), x, xmin=-2, xmax=2, ymin=0, ymax=4)
        sage: (p+q).show()

    As a workaround, we can perform the trick::

        sage: p1 = line([(a,b) for a, b in zip(p[0].xdata, p[0].ydata)
        ....:            if b>=-1 and b<=1])
        sage: q1 = line([(a,b) for a, b in zip(q[0].xdata, q[0].ydata)
        ....:            if b>=0 and b<=4])
        sage: (p1 + q1).show()

    We can also directly plot the elliptic curve::

        sage: E = EllipticCurve([0,-1])                                                 # needs sage.schemes
        sage: plot(E, (1, 4), color=hue(0.6))                                           # needs sage.schemes
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        E = EllipticCurve([0,-1])
        g = plot(E, (1, 4), color=hue(0.6))
        sphinx_plot(g)

    We can change the line style as well::

        sage: plot(sin(x), (x, 0, 10), linestyle='-.')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = plot(sin(x), (x, 0, 10), linestyle='-.')
        sphinx_plot(g)

    If we have an empty linestyle and specify a marker, we can see the
    points that are actually being plotted::

        sage: plot(sin(x), (x,0,10), plot_points=20, linestyle='', marker='.')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = plot(sin(x), (x, 0, 10), plot_points=20, linestyle='', marker='.')
        sphinx_plot(g)

    The marker can be a TeX symbol as well::

        sage: plot(sin(x), (x, 0, 10), plot_points=20, linestyle='', marker=r'$\checkmark$')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = plot(sin(x), (x, 0, 10), plot_points=20, linestyle='', marker=r'$\checkmark$')
        sphinx_plot(g)

    Sage currently ignores points that cannot be evaluated

    ::

        sage: from sage.misc.verbose import set_verbose
        sage: set_verbose(-1)
        sage: plot(-x*log(x), (x, 0, 1))  # this works fine since the failed endpoint is just skipped.
        Graphics object consisting of 1 graphics primitive
        sage: set_verbose(0)

    This prints out a warning and plots where it can (we turn off the
    warning by setting the verbose mode temporarily to -1.)

    ::

        sage: set_verbose(-1)
        sage: plot(x^(1/3), (x, -1, 1))
        Graphics object consisting of 1 graphics primitive
        sage: set_verbose(0)

    .. PLOT::

        set_verbose(-1)
        g = plot(x**(1.0/3.0), (x, -1, 1))
        sphinx_plot(g)
        set_verbose(0)

    Plotting the real cube root function for negative input requires avoiding
    the complex numbers one would usually get.  The easiest way is to use
    :class:`real_nth_root(x, n)<sage.functions.other.Function_real_nth_root>` ::

        sage: plot(real_nth_root(x, 3), (x, -1, 1))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

       g = plot(real_nth_root(x, 3), (x, -1, 1))
       sphinx_plot(g)

    We can also get the same plot in the following way::

        sage: plot(sign(x)*abs(x)^(1/3), (x, -1, 1))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

       g = plot(sign(x)*abs(x)**(1./3.), (x, -1, 1))
       sphinx_plot(g)

    A way to plot other functions without symbolic variants is to use lambda
    functions::

        sage: plot(lambda x : RR(x).nth_root(3), (x,-1, 1))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(plot(lambda x : RR(x).nth_root(3), (x,-1, 1)))

    We can detect the poles of a function::

        sage: plot(gamma, (-3, 4), detect_poles=True).show(ymin=-5, ymax=5)

    .. PLOT::

        g = plot(gamma, (-3, 4), detect_poles=True)
        g.ymin(-5)
        g.ymax(5)
        sphinx_plot(g)

    We draw the Gamma-Function with its poles highlighted::

        sage: plot(gamma, (-3, 4), detect_poles='show').show(ymin=-5, ymax=5)

    .. PLOT::

        g = plot(gamma, (-3, 4), detect_poles='show')
        g.ymin(-5)
        g.ymax(5)
        sphinx_plot(g)

    The basic options for filling a plot::

        sage: p1 = plot(sin(x), -pi, pi, fill='axis')
        sage: p2 = plot(sin(x), -pi, pi, fill='min', fillalpha=1)
        sage: p3 = plot(sin(x), -pi, pi, fill='max')
        sage: p4 = plot(sin(x), -pi, pi, fill=(1-x)/3,
        ....:           fillcolor='blue', fillalpha=.2)
        sage: graphics_array([[p1, p2],                                             # long time
        ....:                 [p3, p4]]).show(frame=True, axes=False)

    .. PLOT::

        p1 = plot(sin(x), -pi, pi, fill='axis'); print(p1)
        p2 = plot(sin(x), -pi, pi, fill='min', fillalpha=1)
        p3 = plot(sin(x), -pi, pi, fill='max')
        p4 = plot(sin(x), -pi, pi, fill=(1-x)/3, fillcolor='blue', fillalpha=.2)
        g = graphics_array([[p1, p2], [p3, p4]])
        sphinx_plot(g, frame=True, axes=False)

    The basic options for filling a list of plots::

        sage: (f1, f2) = x*exp(-1*x^2)/.35, x*exp(-2*x^2)/.35
        sage: p1 = plot([f1, f2], -pi, pi, fill={1: [0]},
        ....:           fillcolor='blue', fillalpha=.25, color='blue')
        sage: p2 = plot([f1, f2], -pi, pi, fill={0: x/3, 1:[0]}, color=['blue'])
        sage: p3 = plot([f1, f2], -pi, pi, fill=[0, [0]],
        ....:           fillcolor=['orange','red'], fillalpha=1, color={1: 'blue'})
        sage: p4 = plot([f1, f2], (x,-pi, pi), fill=[x/3, 0],
        ....:           fillcolor=['grey'], color=['red', 'blue'])
        sage: graphics_array([[p1, p2],                                             # long time
        ....:                 [p3, p4]]).show(frame=True, axes=False)

    .. PLOT::

        (f1, f2) = x*exp(-1*x**2)/.35, x*exp(-2*x**2)/.35
        p1 = plot([f1, f2], -pi, pi, fill={1: [0]}, fillcolor='blue', fillalpha=.25, color='blue')
        p2 = plot([f1, f2], -pi, pi, fill={0: x/3, 1:[0]}, color=['blue'])
        p3 = plot([f1, f2], -pi, pi, fill=[0, [0]], fillcolor=['orange','red'], fillalpha=1, color={1: 'blue'})
        p4 = plot([f1, f2], (x,-pi, pi), fill=[x/3, 0], fillcolor=['grey'], color=['red', 'blue'])
        g = graphics_array([[p1, p2], [p3, p4]])
        sphinx_plot(g, frame=True, axes=False)

    A example about the growth of prime numbers::

        sage: plot(1.13*log(x), 1, 100,
        ....:      fill=lambda x: nth_prime(x)/floor(x), fillcolor='red')
        Graphics object consisting of 2 graphics primitives

    .. PLOT::

        sphinx_plot(plot(1.13*log(x), 1, 100, fill=lambda x: nth_prime(x)/floor(x), fillcolor='red'))

    Fill the area between a function and its asymptote::

        sage: f = (2*x^3+2*x-1)/((x-2)*(x+1))
        sage: plot([f, 2*x+2], -7, 7, fill={0: [1]}, fillcolor='#ccc').show(ymin=-20, ymax=20)

    .. PLOT::

        f = (2*x**3+2*x-1)/((x-2)*(x+1))
        g = plot([f, 2*x+2], -7,7, fill={0: [1]}, fillcolor='#ccc')
        g.ymin(-20)
        g.ymax(20)
        sphinx_plot(g)

    Fill the area between a list of functions and the x-axis::

        sage: def b(n): return lambda x: bessel_J(n, x)
        sage: plot([b(n) for n in [1..5]], 0, 20, fill='axis')
        Graphics object consisting of 10 graphics primitives

    .. PLOT::

        def b(n): return lambda x: bessel_J(n, x)
        g = plot([b(n) for n in range(1,6)], 0, 20, fill='axis')
        sphinx_plot(g)

    Note that to fill between the ith and jth functions, you must use
    the dictionary key-value syntax ``i:[j]``; using key-value pairs
    like ``i:j`` will fill between the ith function and the line y=j::

        sage: def b(n): return lambda x: bessel_J(n, x) + 0.5*(n-1)
        sage: plot([b(c) for c in [1..5]], 0, 20, fill={i: [i-1] for i in [1..4]},
        ....:      color={i: 'blue' for i in [1..5]}, aspect_ratio=3, ymax=3)
        Graphics object consisting of 9 graphics primitives
        sage: plot([b(c) for c in [1..5]], 0, 20, fill={i: i-1 for i in [1..4]},        # long time
        ....:      color='blue', aspect_ratio=3)
        Graphics object consisting of 9 graphics primitives

    .. PLOT::

        def b(n): return lambda x: bessel_J(n, x) + 0.5*(n-1)
        g1 = plot([b(n) for n in range(1,6)], 0, 20, fill={i:[i-1] for i in range(1,5)}, color={i:'blue' for i in range(1,6)}, aspect_ratio=3)
        g2 = plot([b(n) for n in range(1,6)], 0, 20, fill={i:i-1 for i in range(1,5)}, color='blue', aspect_ratio=3)
        g1.ymax(3)
        g = graphics_array([[g1], [g2]])
        sphinx_plot(g)

    Extra options will get passed on to :meth:`~sage.plot.graphics.Graphics.show`,
    as long as they are valid::

        sage: plot(sin(x^2), (x, -3, 3),  # These labels will be nicely typeset
        ....:      title=r'Plot of $\sin(x^2)$', axes_labels=['$x$','$y$'])
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = plot(sin(x**2), (x, -3, 3), title=r'Plot of $\sin(x^2)$', axes_labels=['$x$','$y$'])
        sphinx_plot(g)

    ::

        sage: plot(sin(x^2), (x, -3, 3),  # These will not
        ....:      title='Plot of sin(x^2)', axes_labels=['x','y'])
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = plot(sin(x**2), (x, -3, 3), title='Plot of sin(x^2)', axes_labels=['x','y'])
        sphinx_plot(g)

    ::

        sage: plot(sin(x^2), (x, -3, 3),  # Large axes labels (w.r.t. the tick marks)
        ....:      axes_labels=['x','y'], axes_labels_size=2.5)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = plot(sin(x**2), (x, -3, 3), axes_labels=['x','y'], axes_labels_size=2.5)
        sphinx_plot(g)

    ::

        sage: plot(sin(x^2), (x, -3, 3), figsize=[8,2])
        Graphics object consisting of 1 graphics primitive
        sage: plot(sin(x^2), (x, -3, 3)).show(figsize=[8,2])  # These are equivalent

    .. PLOT::

        g = plot(sin(x**2), (x, -3, 3), figsize=[8,2])
        sphinx_plot(g)

    This includes options for custom ticks and formatting.  See documentation
    for :meth:`show` for more details.

    ::

        sage: plot(sin(pi*x), (x, -8, 8), ticks=[[-7,-3,0,3,7], [-1/2,0,1/2]])
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = plot(sin(pi*x), (x, -8, 8), ticks=[[-7,-3,0,3,7], [-1/2,0,1/2]])
        sphinx_plot(g)

    ::

        sage: plot(2*x + 1, (x, 0, 5),
        ....:      ticks=[[0, 1, e, pi, sqrt(20)],
        ....:             [1, 3, 2*e + 1, 2*pi + 1, 2*sqrt(20) + 1]],
        ....:      tick_formatter='latex')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = plot(2*x + 1, (x, 0, 5), ticks=[[0, 1, e, pi, sqrt(20)], [1, 3, 2*e + 1, 2*pi + 1, 2*sqrt(20) + 1]], tick_formatter='latex')
        sphinx_plot(g)

    This is particularly useful when setting custom ticks in multiples of `\pi`.

    ::

        sage: plot(sin(x), (x,0,2*pi), ticks=pi/3, tick_formatter=pi)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = plot(sin(x),(x,0,2*pi),ticks=pi/3,tick_formatter=pi)
        sphinx_plot(g)

    You can even have custom tick labels along with custom positioning. ::

        sage: plot(x**2, (x,0,3), ticks=[[1,2.5], [0.5,1,2]],
        ....:      tick_formatter=[["$x_1$","$x_2$"],["$y_1$","$y_2$","$y_3$"]])
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = plot(x**2, (x,0,3), ticks=[[1,2.5],[0.5,1,2]], tick_formatter=[["$x_1$","$x_2$"], ["$y_1$","$y_2$","$y_3$"]])
        sphinx_plot(g)

    You can force Type 1 fonts in your figures by providing the relevant
    option as shown below. This also requires that LaTeX, dvipng and
    Ghostscript be installed::

        sage: plot(x, typeset='type1')                          # optional - latex
        Graphics object consisting of 1 graphics primitive

    A example with excluded values::

        sage: plot(floor(x), (x, 1, 10), exclude=[1..10])
        Graphics object consisting of 11 graphics primitives

    .. PLOT::

        g = plot(floor(x), (x, 1, 10), exclude=list(range(1,11)))
        sphinx_plot(g)

    We exclude all points where :class:`~sage.functions.prime_pi.PrimePi`
    makes a jump::

        sage: jumps = [n for n in [1..100] if prime_pi(n) != prime_pi(n-1)]
        sage: plot(lambda x: prime_pi(x), (x, 1, 100), exclude=jumps)
        Graphics object consisting of 26 graphics primitives

    .. PLOT::

        #jumps = [n for n in [1..100] if prime_pi(n) != prime_pi(n-1)]
        #syntaxError: invalid syntax, so we need more code
        jumps=list()
        for n in range(1,101):
            if prime_pi(n) != prime_pi(n-1):
                jumps.append(n)
        g = plot(lambda x: prime_pi(x), (x, 1, 100), exclude=jumps)
        sphinx_plot(g)

    Excluded points can also be given by an equation::

        sage: g(x) = x^2 - 2*x - 2
        sage: plot(1/g(x), (x, -3, 4), exclude=g(x)==0, ymin=-5, ymax=5)    # long time
        Graphics object consisting of 3 graphics primitives

    .. PLOT::

        def g(x): return x**2-2*x-2
        G = plot(1/g(x), (x, -3, 4), exclude=g(x)==0, ymin=-5, ymax=5)
        sphinx_plot(G)

    ``exclude`` and ``detect_poles`` can be used together::

        sage: f(x) = (floor(x)+0.5) / (1-(x-0.5)^2)
        sage: plot(f, (x, -3.5, 3.5), detect_poles='show', exclude=[-3..3],
        ....:      ymin=-5, ymax=5)
        Graphics object consisting of 12 graphics primitives

    .. PLOT::

        def f(x): return (floor(x)+0.5) / (1-(x-0.5)**2)
        g = plot(f, (x, -3.5, 3.5), detect_poles='show', exclude=list(range(-3,4)), ymin=-5, ymax=5)
        sphinx_plot(g)

    Regions in which the plot has no values are automatically excluded. The
    regions thus excluded are in addition to the exclusion points present
    in the ``exclude`` keyword argument.::

        sage: set_verbose(-1)
        sage: plot(arcsec, (x, -2, 2))  # [-1, 1] is excluded automatically
        Graphics object consisting of 2 graphics primitives

    .. PLOT::

        set_verbose(-1)
        g = plot(arcsec, (x, -2, 2))
        sphinx_plot(g)

    ::

        sage: plot(arcsec, (x, -2, 2), exclude=[1.5])  # x=1.5 is also excluded
        Graphics object consisting of 3 graphics primitives

    .. PLOT::

        set_verbose(-1)
        g = plot(arcsec, (x, -2, 2), exclude=[1.5])
        sphinx_plot(g)

    ::

        sage: plot(arcsec(x/2), -2, 2)  # plot should be empty; no valid points
        Graphics object consisting of 0 graphics primitives
        sage: plot(sqrt(x^2 - 1), -2, 2)  # [-1, 1] is excluded automatically
        Graphics object consisting of 2 graphics primitives

    .. PLOT::

        set_verbose(-1)
        g = plot(sqrt(x**2-1), -2, 2)
        sphinx_plot(g)

    ::

        sage: plot(arccsc, -2, 2)       # [-1, 1] is excluded automatically
        Graphics object consisting of 2 graphics primitives
        sage: set_verbose(0)

    .. PLOT::

        set_verbose(-1)
        g = plot(arccsc, -2, 2)
        sphinx_plot(g)

    TESTS:

    We do not randomize the endpoints::

        sage: p = plot(x, (x,-1,1))
        sage: p[0].xdata[0] == -1
        True
        sage: p[0].xdata[-1] == 1
        True

    We check to make sure that the x/y min/max data get set correctly
    when there are multiple functions.

    ::

        sage: d = plot([sin(x), cos(x)], 100, 120).get_minmax_data()
        sage: d['xmin']
        100.0
        sage: d['xmax']
        120.0

    We check various combinations of tuples and functions, ending with
    tests that lambda functions work properly with explicit variable
    declaration, without a tuple.

    ::

        sage: p = plot(lambda x: x,(x,-1,1))
        sage: p = plot(lambda x: x,-1,1)
        sage: p = plot(x,x,-1,1)
        sage: p = plot(x,-1,1)
        sage: p = plot(x^2,x,-1,1)
        sage: p = plot(x^2,xmin=-1,xmax=2)
        sage: p = plot(lambda x: x,x,-1,1)
        sage: p = plot(lambda x: x^2,x,-1,1)
        sage: p = plot(lambda x: 1/x,x,-1,1)
        sage: f(x) = sin(x+3)-.1*x^3
        sage: p = plot(lambda x: f(x),x,-1,1)

    We check to handle cases where the function gets evaluated at a
    point which causes an 'inf' or '-inf' result to be produced.

    ::

        sage: p = plot(1/x, 0, 1)
        sage: p = plot(-1/x, 0, 1)

    Bad options now give better errors::

        sage: P = plot(sin(1/x), (x,-1,3), foo=10)
        Traceback (most recent call last):
        ...
        RuntimeError: error in line(): option 'foo' not valid
        sage: P = plot(x, (x,1,1)) # github issue #11753
        Traceback (most recent call last):
        ...
        ValueError: plot start point and end point must be different

    We test that we can plot `f(x)=x` (see :issue:`10246`)::

        sage: f(x)=x; f
        x |--> x
        sage: plot(f,(x,-1,1))
        Graphics object consisting of 1 graphics primitive

    Check that :issue:`15030` is fixed::

        sage: plot(abs(log(x)), x)
        Graphics object consisting of 1 graphics primitive

    Check that if excluded points are less than xmin then the exclusion
    still works for polar and parametric plots. The following should
    show two excluded points::

        sage: set_verbose(-1)
        sage: polar_plot(sin(sqrt(x^2-1)), (x,0,2*pi), exclude=[1/2,2,3])
        Graphics object consisting of 3 graphics primitives

        sage: parametric_plot((sqrt(x^2-1),sqrt(x^2-1/2)), (x,0,5), exclude=[1,2,3])
        Graphics object consisting of 3 graphics primitives

        sage: set_verbose(0)

    Legends can contain variables with long names, :issue:`13543`::

        sage: hello = var('hello')
        sage: label = '$' + latex(hello) + '$'
        sage: plot(x, x, 0, 1, legend_label=label)
        Graphics object consisting of 1 graphics primitive

    Extra keywords should be saved if object has a plot method, :issue:`20924`::

        sage: G = graphs.PetersenGraph()
        sage: p = G.plot()
        sage: p.aspect_ratio()
        1.0
        sage: pp = plot(G)
        sage: pp.aspect_ratio()
        1.0
    """
    if 'color' in kwds and 'rgbcolor' in kwds:
        raise ValueError('only one of color or rgbcolor should be specified')
    elif 'color' in kwds:
        kwds['rgbcolor'] = kwds.pop('color', (0,0,1)) # take blue as default ``rgbcolor``
    G_kwds = Graphics._extract_kwds_for_show(kwds, ignore=['xmin', 'xmax'])
    if 'scale' in G_kwds:
        kwds['scale'] = G_kwds['scale'] # pass scaling information to _plot too

    original_opts = kwds.pop('__original_opts', {})
    do_show = kwds.pop('show',False)

    from sage.structure.element import Vector
    if kwds.get('parametric',False) and isinstance(funcs, Vector):
        funcs = tuple(funcs)

    if hasattr(funcs, 'plot'):
        G = funcs.plot(*args, **original_opts)

        # If we have extra keywords already set, then update them
        for ext in G._extra_kwds:
            if ext in G_kwds:
                G_kwds[ext] = G._extra_kwds[ext]

    # if we are using the generic plotting method
    else:
        n = len(args)
        # if there are no extra args, try to get xmin,xmax from
        # keyword arguments or pick some silly default
        if n == 0:
            xmin = kwds.pop('xmin', -1)
            xmax = kwds.pop('xmax', 1)
            G = _plot(funcs, (xmin, xmax), **kwds)

        # if there is one extra arg, then it had better be a tuple
        elif n == 1:
            G = _plot(funcs, *args, **kwds)
        elif n == 2:
        # if there are two extra args, then pull them out and pass them as a tuple
            xmin = args[0]
            xmax = args[1]
            args = args[2:]
            G = _plot(funcs, (xmin, xmax), *args, **kwds)
        elif n == 3:
        # if there are three extra args, then pull them out and pass them as a tuple
            var = args[0]
            xmin = args[1]
            xmax = args[2]
            args = args[3:]
            G = _plot(funcs, (var, xmin, xmax), *args, **kwds)
        elif ('xmin' in kwds) or ('xmax' in kwds):
            xmin = kwds.pop('xmin', -1)
            xmax = kwds.pop('xmax', 1)
            G = _plot(funcs, (xmin, xmax), *args, **kwds)
        else:
            sage.misc.verbose.verbose(f"there were {n} extra arguments (besides {funcs})", level=0)

    G._set_extra_kwds(G_kwds)
    if do_show:
        G.show()
    return G


def _plot(funcs, xrange, parametric=False,
          polar=False, fill=False, label='', randomize=True, **options):
    """
    Internal function which does the actual plotting.

    INPUT:

    - ``funcs`` -- function or list of functions to be plotted
    - ``xrange`` -- two or three tuple of [input variable], min and max
    - ``parametric`` -- boolean (default: ``False``); whether
      this is a parametric plot
    - ``polar`` -- boolean (default: ``False``); whether
      this is a polar plot
    - ``fill`` -- boolean (default: ``False``); whether
      this plot is filled
    - ``randomize`` -- boolean (default: ``True``); whether
      to use random plot points

    The following option is deprecated in favor of ``legend_label``:

    - ``label`` -- (default: ``''``) a string for the label

    All other usual plot options are also accepted, and a number
    are required (see the example below) which are normally passed
    through the options decorator to :func:`plot`.

    OUTPUT: a ``Graphics`` object

    EXAMPLES:

    See :func:`plot` for many, many implicit examples.
    Here is an explicit one::

        sage: from sage.plot.plot import _plot
        sage: P = _plot(e^(-x^2),(-3,3), fill=True,
        ....:                            color='red',
        ....:                            plot_points=50,
        ....:                            adaptive_tolerance=2,
        ....:                            adaptive_recursion=True,
        ....:                            exclude=None,
        ....:                            imaginary_tolerance=1e-8)
        sage: P.show(aspect_ratio='automatic')

    TESTS:

    Make sure that we get the right number of legend entries as the number of
    functions varies (:issue:`10514`)::

        sage: p1 = plot(1*x, legend_label='1x')
        sage: p2 = plot(2*x, legend_label='2x', color='green')
        sage: p1+p2
        Graphics object consisting of 2 graphics primitives

    ::

        sage: len(p1.matplotlib().axes[0].legend().texts)
        1
        sage: len((p1+p2).matplotlib().axes[0].legend().texts)
        2
        sage: q1 = plot([sin(x), tan(x)], color='blue', legend_label='trig')
        sage: len(q1.matplotlib().axes[0].legend().texts)
        2
        sage: q1
        Graphics object consisting of 2 graphics primitives
        sage: q2 = plot([sin(x), tan(x)], legend_label={1:'tan'})
        sage: len(q2.matplotlib().axes[0].legend().texts)
        2
        sage: q2
        Graphics object consisting of 2 graphics primitives

    Make sure that we don't get multiple legend labels for plot segments
    (:issue:`11998`)::

        sage: p1 = plot(1/(x^2-1),(x,-2,2),legend_label='foo',detect_poles=True)
        sage: len(p1.matplotlib().axes[0].legend().texts)
        1
        sage: p1.show(ymin=-10,ymax=10) # should be one legend

    Parametric plots that get evaluated at invalid points should still
    plot properly (:issue:`13246`)::

        sage: parametric_plot((x, arcsec(x)), (x, -2, 2))
        Graphics object consisting of 1 graphics primitive

    Verify that :issue:`31089` is fixed::

        sage: plot(x, -1, 1, detect_poles=True)
        Graphics object consisting of 1 graphics primitive
    """
    from sage.plot.colors import Color
    from sage.plot.misc import setup_for_eval_on_grid
    if funcs == []:
        return Graphics()
    orig_funcs = funcs # keep the original functions (for use in legend labels)
    excluded_points = []
    imag_tol = options["imaginary_tolerance"]
    funcs, ranges = setup_for_eval_on_grid(funcs,
                                           [xrange],
                                           options['plot_points'],
                                           imaginary_tolerance=imag_tol)
    xmin, xmax, delta = ranges[0]
    xrange = ranges[0][:2]
    # parametric_plot will be a list or tuple of two functions (f,g)
    # and will plotted as (f(x), g(x)) for all x in the given range
    if parametric:
        f, g = funcs
    # or we have only a single function to be plotted:
    else:
        f = funcs

    # Keep ``rgbcolor`` option 'automatic' only for lists of functions.
    # Otherwise, let the plot color be 'blue'.
    if parametric or not isinstance(funcs, (list, tuple)):
        if 'rgbcolor' in options.keys() and options['rgbcolor'] == 'automatic':
            options['rgbcolor'] = (0, 0, 1) # default color for a single curve.

    # Check to see if funcs is a list of functions that will be all plotted together.
    if isinstance(funcs, (list, tuple)) and not parametric:
        def golden_rainbow(i,lightness=0.4):
            # note: sage's "blue" has hue-saturation-lightness values (2/3, 1, 1/2).
            g = 0.61803398875
            return Color((0.66666666666666 + i*g) % 1, 1, lightness, space='hsl')

        default_line_styles = ("-", "--", "-.", ":")*len(funcs)

        G = Graphics()
        for i, h in enumerate(funcs):
            options_temp = options.copy()
            color_temp = options_temp.pop('rgbcolor', 'automatic')
            fill_temp = options_temp.pop('fill', fill)
            fillcolor_temp = options_temp.pop('fillcolor', 'automatic') # perhaps the 2nd argument should be ``options_temp['color']``
            linestyle_temp = options_temp.pop('linestyle', None)
            legend_label_temp = options_temp.pop('legend_label', None)
            legend_color_temp = options_temp.pop('legend_color', None)

            # passed more than one color directive?
            one_plot_color = False
            if isinstance(color_temp, dict):
                if i in color_temp:
                    color_entry = color_temp[i]
                else:
                    color_entry = golden_rainbow(i)
            elif isinstance(color_temp, (list, tuple)) and isinstance(color_temp[0], (str, list, tuple)):
                if i < len(color_temp):
                    color_entry = color_temp[i]
                else:
                    color_entry = golden_rainbow(i)
            elif color_temp == 'automatic':
                color_entry = golden_rainbow(i)
            else:
                # assume a single color was assigned for all plots
                one_plot_color = True
                color_entry = color_temp

            # passed more than one fill directive?
            if isinstance(fill_temp, dict):
                if i in fill_temp:
                    fill_entry_listQ = fill_temp[i]
                    if isinstance(fill_entry_listQ, list):
                        if fill_entry_listQ[0] < len(funcs):
                            fill_entry = funcs[fill_entry_listQ[0]]
                        else:
                            fill_entry = False
                    else:
                        fill_entry = fill_entry_listQ
                else:
                    fill_entry = False
            elif isinstance(fill_temp, (list, tuple)):
                if i < len(fill_temp):
                    fill_entry_listQ = fill_temp[i]
                    if isinstance(fill_entry_listQ, list):
                        if fill_entry_listQ[0] < len(funcs):
                            fill_entry = funcs[fill_entry_listQ[0]]
                        else:
                            fill_entry = False
                    else:
                        fill_entry = fill_entry_listQ
                else:
                    fill_entry = False
            else:
                fill_entry = fill_temp

            # passed more than one fillcolor directive?
            fillcolor_entry = 'gray' # the default choice
            if isinstance(fillcolor_temp, dict):
                if i in fillcolor_temp:
                    fillcolor_entry = fillcolor_temp[i]
                else:
                    fillcolor_entry = golden_rainbow(i,0.85)
            elif isinstance(fillcolor_temp, (list, tuple)):
                if i < len(fillcolor_temp):
                    fillcolor_entry = fillcolor_temp[i]
                else:
                    fillcolor_entry = golden_rainbow(i,0.85)
            elif fillcolor_temp == 'automatic':
                # check that we haven't overwritten automatic multi-colors in color_temp
                if len(funcs) > 1 and not one_plot_color:
                    fillcolor_entry = golden_rainbow(i,0.85)
            elif fillcolor_temp is not None:
                fillcolor_entry = fillcolor_temp

            # passed more than one linestyle directive?
            if isinstance(linestyle_temp, dict):
                if i in linestyle_temp:
                    linestyle_entry = linestyle_temp[i]
                else:
                    linestyle_entry = default_line_styles[i]
            elif isinstance(linestyle_temp, (list, tuple)):
                if i < len(linestyle_temp):
                    linestyle_entry = linestyle_temp[i]
                else:
                    linestyle_entry = default_line_styles[i]
            elif linestyle_temp == 'automatic':
                linestyle_entry = default_line_styles[i]
            elif linestyle_temp is None:
                linestyle_entry = 'solid'
            else:
                linestyle_entry = linestyle_temp

            # passed more than one legend_label directive?
            legend_label_entry = orig_funcs[i].__repr__() # the 'automatic' choice
            if isinstance(legend_label_temp, dict):
                if i in legend_label_temp:
                    legend_label_entry = legend_label_temp[i]
            elif isinstance(legend_label_temp, (list, tuple)):
                if i < len(legend_label_temp):
                    legend_label_entry = legend_label_temp[i]
            elif legend_label_temp != 'automatic':
                # assume it is None or str.
                legend_label_entry = legend_label_temp

            # passed more than one legend_color directive?
            legend_color_entry = 'black'  # the default choice
            if isinstance(legend_color_temp, dict):
                if i in legend_color_temp:
                    legend_color_entry = legend_color_temp[i]
            elif isinstance(legend_color_temp, (list, tuple)) and isinstance(legend_color_temp[0], (str, list, tuple)):
                if i < len(legend_color_temp):
                    legend_color_entry = legend_color_temp[i]
            elif legend_color_temp == 'automatic':
                if len(funcs) > 1:
                    legend_color_entry = golden_rainbow(i)
            elif legend_color_temp is not None:
                legend_color_entry = legend_color_temp

            G += plot(h, xrange, polar=polar, fill=fill_entry, fillcolor=fillcolor_entry,
                      rgbcolor=color_entry, linestyle=linestyle_entry,
                      legend_label=legend_label_entry, legend_color=legend_color_entry, **options_temp)
        return G

    adaptive_tolerance = options.pop('adaptive_tolerance')
    adaptive_recursion = options.pop('adaptive_recursion')
    imag_tol = options.pop('imaginary_tolerance')
    plot_points = int(options.pop('plot_points'))

    exclude = options.pop('exclude')
    if exclude is not None:
        if isinstance(exclude, Expression) and exclude.is_relational():
            if len(exclude.variables()) > 1:
                raise ValueError('exclude has to be an equation of only one variable')
            v = exclude.variables()[0]
            points = [e.right() for e in exclude.solve(v) if e.left() == v and (v not in e.right().variables())]
            # We are only interested in real solutions
            for x in points:
                try:
                    excluded_points.append(float(x))
                except TypeError:
                    pass
            excluded_points.sort()

        # We should either have a list in excluded points or exclude
        # itself must be a list
        elif isinstance(exclude, (list, tuple)):
            excluded_points = sorted(exclude)
        else:
            raise ValueError('exclude needs to be a list of numbers or an equation')

        # We make sure that points plot points close to the excluded points are computed
        epsilon = 0.001*(xmax - xmin)
        initial_points = reduce(lambda a,b: a+b,
                                [[x - epsilon, x + epsilon]
                                 for x in excluded_points], [])
    else:
        initial_points = None

    # If we are a log scale plot on the x axis, do a change of variables
    # so we sample the range in log scale
    is_log_scale = ('scale' in options.keys() and
                    not parametric and
                    options['scale'] in ['loglog', 'semilogx'])
    if is_log_scale:
        f_exp = lambda x: f(exp(x))
        log_xrange = (log(xrange[0]), log(xrange[1]))
        if initial_points is None:
            log_initial_points = None
        else:
            log_initial_points = [log(x) for x in initial_points]
        data, extra_excluded = generate_plot_points(
            f_exp, log_xrange, plot_points,
            adaptive_tolerance, adaptive_recursion,
            randomize, log_initial_points,
            excluded=True, imaginary_tolerance=imag_tol)
    else:
        data, extra_excluded = generate_plot_points(
            f, xrange, plot_points,
            adaptive_tolerance, adaptive_recursion,
            randomize, initial_points,
            excluded=True, imaginary_tolerance=imag_tol)

    excluded_points += extra_excluded

    # If we did a change in variables, undo it now
    if is_log_scale:
        for i,(a,fa) in enumerate(data):
            data[i] = (exp(a), fa)
        for i,p in enumerate(excluded_points):
            excluded_points[i] = exp(p)

    if parametric:
        # We need the original x-values to be able to exclude points in parametric plots
        exclude_data = data
        newdata = []
        for x,fdata in data:
            try:
                newdata.append((fdata, g(x)))
            except (ValueError, TypeError):
                newdata.append((fdata, 0)) # append a dummy value 0
                excluded_points.append(x)
        data = newdata

    excluded_points.sort(reverse=True)
    G = Graphics()

    fillcolor = options.pop('fillcolor', 'automatic')
    fillalpha = options.pop('fillalpha', 0.5)

    # TODO: Use matplotlib's fill and fill_between commands.
    if fill is not False and fill is not None:
        if parametric:
            filldata = data
        else:
            if fill == 'axis' or fill is True:
                base_level = 0
            elif fill == 'min':
                base_level = min(t[1] for t in data)
            elif fill == 'max':
                base_level = max(t[1] for t in data)
            elif callable(fill):
                if fill == max or fill == min:
                    if fill == max:
                        fstr = 'max'
                    else:
                        fstr = 'min'
                    msg = "WARNING: You use the built-in function {} for filling. You probably wanted the string '{}'.".format(fstr, fstr)
                    sage.misc.verbose.verbose(msg, level=0)
                if not is_fast_float(fill):
                    fill_f = fast_float(fill, expect_one_var=True)
                else:
                    fill_f = fill

                filldata = generate_plot_points(fill_f,
                                                xrange,
                                                plot_points,
                                                adaptive_tolerance,
                                                adaptive_recursion,
                                                randomize,
                                                imaginary_tolerance=imag_tol)
                filldata.reverse()
                filldata += data
            else:
                try:
                    base_level = float(fill)
                except TypeError:
                    base_level = 0

            if not callable(fill) and polar:
                filldata = generate_plot_points(lambda x: base_level,
                                                xrange,
                                                plot_points,
                                                adaptive_tolerance,
                                                adaptive_recursion,
                                                randomize,
                                                imaginary_tolerance=imag_tol)
                filldata.reverse()
                filldata += data
            if not callable(fill) and not polar:
                filldata = [(data[0][0], base_level)] + data + [(data[-1][0], base_level)]

        if fillcolor == 'automatic':
            fillcolor = (0.5, 0.5, 0.5)
        fill_options = {}
        fill_options['rgbcolor'] = fillcolor
        fill_options['alpha'] = fillalpha
        fill_options['thickness'] = 0
        if polar:
            filldata = [(y*cos(x), y*sin(x)) for x, y in filldata]
        G += polygon(filldata, **fill_options)

    # We need the original data to be able to exclude points in polar plots
    if not parametric:
        exclude_data = data
    if polar:
        data = [(y*cos(x), y*sin(x)) for x, y in data]

    detect_poles = options.pop('detect_poles', False)
    legend_label = options.pop('legend_label', None)
    if excluded_points or detect_poles:
        start_index = 0
        # setup for pole detection
        from sage.rings.real_double import RDF
        epsilon = 0.0001
        pole_options = {}
        pole_options['linestyle'] = '--'
        pole_options['thickness'] = 1
        pole_options['rgbcolor'] = '#ccc'

        # setup for exclusion points
        if excluded_points:
            exclusion_point = excluded_points.pop()
        else:
            # default value of exclusion point must be outside the plot interval
            # (see :issue:`31089`)
            exclusion_point = xmax + 1

        flag = True
        for i in range(len(data)-1):
            x0, y0 = exclude_data[i]
            x1, y1 = exclude_data[i+1]

            # detect poles
            if (not (polar or parametric)) and detect_poles \
               and ((y1 > 0 and y0 < 0) or (y1 < 0 and y0 > 0)):
                # calculate the slope of the line segment
                dy = abs(y1-y0)
                dx = x1 - x0
                alpha = (RDF(dy)/RDF(dx)).arctan()
                if alpha >= RDF(pi/2) - epsilon:
                    G += line(data[start_index:i], **options)
                    if detect_poles == 'show':
                        # draw a vertical asymptote
                        G += line([(x0, y0), (x1, y1)], **pole_options)
                    start_index = i+2

            # exclude points
            if x0 > exclusion_point:
                while exclusion_point <= x1:
                    try:
                        exclusion_point = excluded_points.pop()
                    except IndexError:
                        # all excluded points were considered
                        flag = False
                        break

            elif flag and (x0 <= exclusion_point <= x1):
                G += line(data[start_index:i], **options)
                start_index = i + 2
                while exclusion_point <= x1:
                    try:
                        exclusion_point = excluded_points.pop()
                    except IndexError:
                        # all excluded points were considered
                        flag = False
                        break

        G += line(data[start_index:], legend_label=legend_label, **options)
    else:
        G += line(data, legend_label=legend_label, **options)

    return G


########## misc functions ###################

@options(aspect_ratio=1.0)
def parametric_plot(funcs, *args, **kwargs):
    r"""
    Plot a parametric curve or surface in 2d or 3d.

    :func:`parametric_plot` takes two or three functions as a
    list or a tuple and makes a plot with the first function giving the
    `x` coordinates, the second function giving the `y`
    coordinates, and the third function (if present) giving the
    `z` coordinates.

    In the 2d case, :func:`parametric_plot` is equivalent to the :func:`plot` command
    with the option ``parametric=True``.  In the 3d case, :func:`parametric_plot`
    is equivalent to :func:`~sage.plot.plot3d.parametric_plot3d.parametric_plot3d`.
    See each of these functions for more help and examples.

    INPUT:

    - ``funcs`` -- 2 or 3-tuple of functions, or a vector of dimension 2 or 3

    - ``other options`` -- passed to :func:`plot` or :func:`~sage.plot.plot3d.parametric_plot3d.parametric_plot3d`

    EXAMPLES: We draw some 2d parametric plots.  Note that the default aspect ratio
    is 1, so that circles look like circles. ::

        sage: t = var('t')
        sage: parametric_plot((cos(t), sin(t)), (t, 0, 2*pi))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        t = var('t')
        g = parametric_plot( (cos(t), sin(t)), (t, 0, 2*pi))
        sphinx_plot(g)

    ::

        sage: parametric_plot((sin(t), sin(2*t)), (t, 0, 2*pi), color=hue(0.6))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        t = var('t')
        g = parametric_plot( (sin(t), sin(2*t)), (t, 0, 2*pi), color=hue(0.6))
        sphinx_plot(g)

    ::

        sage: parametric_plot((1, t), (t, 0, 4))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        t =var('t')
        g = parametric_plot((1, t), (t, 0, 4))
        sphinx_plot(g)

    Note that in parametric_plot, there is only fill or no fill.

    ::

        sage: parametric_plot((t, t^2), (t, -4, 4), fill=True)
        Graphics object consisting of 2 graphics primitives

    .. PLOT::

        t = var('t')
        g = parametric_plot((t, t**2), (t, -4, 4), fill=True)
        sphinx_plot(g)

    A filled Hypotrochoid::

        sage: parametric_plot([cos(x) + 2 * cos(x/4), sin(x) - 2 * sin(x/4)],
        ....:                 (x, 0, 8*pi), fill=True)
        Graphics object consisting of 2 graphics primitives

    .. PLOT::

        g = parametric_plot([cos(x) + 2 * cos(x/4), sin(x) - 2 * sin(x/4)], (x, 0, 8*pi), fill=True)
        sphinx_plot(g)

    ::

        sage: parametric_plot((5*cos(x), 5*sin(x), x), (x, -12, 12),  # long time
        ....:                 plot_points=150, color='red')
        Graphics3d Object

    .. PLOT::

        #AttributeError: 'Line' object has no attribute 'plot'
        #g = parametric_plot( (5*cos(x), 5*sin(x), x), (x,-12, 12), plot_points=150, color='red') # long time
        #sphinx_plot(g)

    ::

        sage: y = var('y')
        sage: parametric_plot((5*cos(x), x*y, cos(x*y)), (x, -4, 4), (y, -4, 4))  # long time
        Graphics3d Object

    .. PLOT::

        #AttributeError: 'sage.plot.plot3d.parametric_surface.ParametricSurf' object has no attribute 'plot'
        #y = var('y')
        #g = parametric_plot( (5*cos(x), x*y, cos(x*y)), (x, -4,4), (y,-4,4))  # long time
        #sphinx_plot(g)

    ::

        sage: t = var('t')
        sage: parametric_plot(vector((sin(t), sin(2*t))), (t, 0, 2*pi), color='green')  # long time
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        t = var('t')
        g = parametric_plot( vector((sin(t), sin(2*t))), (t, 0, 2*pi), color='green')  # long time
        sphinx_plot(g)

    ::

        sage: t = var('t')
        sage: parametric_plot( vector([t, t+1, t^2]), (t, 0, 1)) # long time
        Graphics3d Object

    .. PLOT::

        #t = var('t')
        #g = parametric_plot( vector([t, t+1, t**2]), (t, 0, 1)) # long time
        #sphinx_plot(g)

    Plotting in logarithmic scale is possible with 2D plots. The keyword
    ``aspect_ratio`` will be ignored if the scale is not ``'loglog'`` or
    ``'linear'``.::

        sage: parametric_plot((x, x**2), (x, 1, 10), scale='loglog')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = parametric_plot((x, x**2), (x, 1, 10), scale='loglog')
        sphinx_plot(g)

    We can also change the scale of the axes in the graphics just before
    displaying. In this case, the ``aspect_ratio`` must be specified as
    ``'automatic'`` if the ``scale`` is set to ``'semilogx'`` or ``'semilogy'``. For
    other values of the ``scale`` parameter, any ``aspect_ratio`` can be
    used, or the keyword need not be provided.::

        sage: p = parametric_plot((x, x**2), (x, 1, 10))
        sage: p.show(scale='semilogy', aspect_ratio='automatic')

    TESTS::

        sage: parametric_plot((x, t^2), (x, -4, 4))
        Traceback (most recent call last):
        ...
        ValueError: there are more variables than variable ranges

        sage: parametric_plot((1, x+t), (x, -4, 4))
        Traceback (most recent call last):
        ...
        ValueError: there are more variables than variable ranges

        sage: parametric_plot((-t, x+t), (x, -4, 4))
        Traceback (most recent call last):
        ...
        ValueError: there are more variables than variable ranges

        sage: parametric_plot((1, x+t, y), (x, -4, 4), (t, -4, 4))
        Traceback (most recent call last):
        ...
        ValueError: there are more variables than variable ranges

        sage: parametric_plot((1, x, y), 0, 4)
        Traceback (most recent call last):
        ...
        ValueError: there are more variables than variable ranges

    One test for :issue:`7165`::

        sage: m = SR.var('m')
        sage: parametric_plot([real(exp(i*m)), imaginary(exp(i*m))], (m, 0, 7))
        Graphics object consisting of 1 graphics primitive
    """
    num_ranges = 0
    for i in args:
        if isinstance(i, (list, tuple)):
            num_ranges += 1
        else:
            break

    num_funcs = len(funcs)

    num_vars = len(sage.plot.misc.unify_arguments(funcs)[0])
    if num_vars > num_ranges:
        raise ValueError("there are more variables than variable ranges")

    # Reset aspect_ratio to 'automatic' in case scale is 'semilog[xy]'.
    # Otherwise matplotlib complains.
    scale = kwargs.get('scale', None)
    if isinstance(scale, (list, tuple)):
        scale = scale[0]
    if scale == 'semilogy' or scale == 'semilogx':
        kwargs['aspect_ratio'] = 'automatic'

    if num_funcs == 2 and num_ranges == 1:
        kwargs['parametric'] = True
        return plot(funcs, *args, **kwargs)
    elif (num_funcs == 3 and num_ranges <= 2):
        return sage.plot.plot3d.parametric_plot3d.parametric_plot3d(funcs, *args, **kwargs)
    else:
        raise ValueError("the number of functions and the number of variable ranges is not a supported combination for a 2d or 3d parametric plots")


@options(aspect_ratio=1.0)
def polar_plot(funcs, *args, **kwds):
    r"""
    ``polar_plot`` takes a single function or a list or
    tuple of functions and plots them with polar coordinates in the given
    domain.

    This function is equivalent to the :func:`plot` command with the options
    ``polar=True`` and ``aspect_ratio=1``. For more help on options,
    see the documentation for :func:`plot`.

    INPUT:

    - ``funcs`` -- a function

    - other options are passed to plot

    EXAMPLES:

    Here is a blue 8-leaved petal::

        sage: polar_plot(sin(5*x)^2, (x, 0, 2*pi), color='blue')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = polar_plot(sin(5*x)**2, (x, 0, 2*pi), color='blue')
        sphinx_plot(g)

    A red figure-8::

        sage: polar_plot(abs(sqrt(1 - sin(x)^2)), (x, 0, 2*pi), color='red')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = polar_plot(abs(sqrt(1 - sin(x)**2)), (x, 0, 2*pi), color='red')
        sphinx_plot(g)

    A green limacon of Pascal::

        sage: polar_plot(2 + 2*cos(x), (x, 0, 2*pi), color=hue(0.3))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = polar_plot(2 + 2*cos(x), (x, 0, 2*pi), color=hue(0.3))
        sphinx_plot(g)

    Several polar plots::

        sage: polar_plot([2*sin(x), 2*cos(x)], (x, 0, 2*pi))
        Graphics object consisting of 2 graphics primitives

    .. PLOT::

        g = polar_plot([2*sin(x), 2*cos(x)], (x, 0, 2*pi))
        sphinx_plot(g)

    A filled spiral::

        sage: polar_plot(sqrt, 0, 2 * pi, fill=True)
        Graphics object consisting of 2 graphics primitives

    .. PLOT::

        g = polar_plot(sqrt, 0, 2 * pi, fill=True)
        sphinx_plot(g)

    Fill the area between two functions::

        sage: polar_plot(cos(4*x) + 1.5, 0, 2*pi, fill=0.5 * cos(4*x) + 2.5,
        ....:            fillcolor='orange')
        Graphics object consisting of 2 graphics primitives

    .. PLOT::

        g = polar_plot(cos(4*x) + 1.5, 0, 2*pi, fill=0.5 * cos(4*x) + 2.5, fillcolor='orange')
        sphinx_plot(g)

    Fill the area between several spirals::

        sage: polar_plot([(1.2+k*0.2)*log(x) for k in range(6)], 1, 3 * pi,
        ....:            fill={0: [1], 2: [3], 4: [5]})
        Graphics object consisting of 9 graphics primitives

    .. PLOT::

        g = polar_plot([(1.2+k*0.2)*log(x) for k in range(6)], 1, 3 * pi, fill={0: [1], 2: [3], 4: [5]})
        sphinx_plot(g)

    Exclude points at discontinuities::

        sage: polar_plot(log(floor(x)), (x, 1, 4*pi), exclude=[1..12])
        Graphics object consisting of 12 graphics primitives

    .. PLOT::

        g = polar_plot(log(floor(x)), (x, 1, 4*pi), exclude=list(range(1,13)))
        sphinx_plot(g)
    """
    kwds['polar'] = True
    return plot(funcs, *args, **kwds)


@options(aspect_ratio='automatic')
def list_plot(data, plotjoined=False, **kwargs):
    r"""
    ``list_plot`` takes either a list of numbers, a list of tuples, a numpy
    array, or a dictionary and plots the corresponding points.

    If given a list of numbers (that is, not a list of tuples or lists),
    ``list_plot`` forms a list of tuples ``(i, x_i)`` where ``i`` goes from
    0 to ``len(data)-1`` and ``x_i`` is the ``i``-th data value, and puts
    points at those tuple values.

    ``list_plot`` will plot a list of complex numbers in the obvious
    way; any numbers for which
    :func:`CC()<sage.rings.complex_mpfr.ComplexField>` makes sense will
    work.

    ``list_plot`` also takes a list of tuples ``(x_i, y_i)`` where ``x_i``
    and ``y_i`` are the ``i``-th values representing the ``x``- and
    ``y``-values, respectively.

    If given a dictionary, ``list_plot`` interprets the keys as
    `x`-values and the values as `y`-values.

    The ``plotjoined=True`` option tells ``list_plot`` to plot a line
    joining all the data.

    For other keyword options that the ``list_plot`` function can
    take, refer to :func:`~sage.plot.plot.plot`.

    It is possible to pass empty dictionaries, lists, or tuples to
    ``list_plot``. Doing so will plot nothing (returning an empty plot).

    EXAMPLES::

        sage: list_plot([i^2 for i in range(5)])  # long time
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = list_plot([i**2 for i in range(5)]) # long time
        sphinx_plot(g)

    Here are a bunch of random red points::

        sage: r = [(random(),random()) for _ in range(20)]
        sage: list_plot(r, color='red')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        r = [(random(),random()) for _ in range(20)]
        g = list_plot(r, color='red')
        sphinx_plot(g)

    This gives all the random points joined in a purple line::

        sage: list_plot(r, plotjoined=True, color='purple')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        r = [(random(),random()) for _ in range(20)]
        g = list_plot(r, plotjoined=True, color='purple')
        sphinx_plot(g)

    You can provide a numpy array.::

        sage: import numpy                                                              # needs numpy
        sage: list_plot(numpy.arange(10))                                               # needs numpy
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        import numpy
        g = list_plot(numpy.arange(10))
        sphinx_plot(g)

    ::

        sage: list_plot(numpy.array([[1,2], [2,3], [3,4]]))                             # needs numpy
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        import numpy
        g = list_plot(numpy.array([[1,2], [2,3], [3,4]]))
        sphinx_plot(g)

    Plot a list of complex numbers::

        sage: list_plot([1, I, pi + I/2, CC(.25, .25)])
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = list_plot([1, I, pi + I/2, CC(.25, .25)])
        sphinx_plot(g)

    ::

        sage: list_plot([exp(I*theta) for theta in [0, .2..pi]])
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = list_plot([exp(I*theta) for theta in srange(0,pi,0.2)])
        sphinx_plot(g)

    Note that if your list of complex numbers are all actually real,
    they get plotted as real values, so this

    ::

        sage: list_plot([CDF(1), CDF(1/2), CDF(1/3)])
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = list_plot([CDF(1), CDF(1/2), CDF(1/3)])
        sphinx_plot(g)

    is the same as ``list_plot([1, 1/2, 1/3])`` -- it produces a plot of
    the points `(0,1)`, `(1,1/2)`, and `(2,1/3)`.

    If you have separate lists of `x` values and `y` values which you
    want to plot against each other, use the ``zip`` command to make a
    single list whose entries are pairs of `(x,y)` values, and feed
    the result into ``list_plot``::

        sage: x_coords = [cos(t)^3 for t in srange(0, 2*pi, 0.02)]
        sage: y_coords = [sin(t)^3 for t in srange(0, 2*pi, 0.02)]
        sage: list_plot(list(zip(x_coords, y_coords)))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        x_coords = [cos(t)**3 for t in srange(0, 2*pi, 0.02)]
        y_coords = [sin(t)**3 for t in srange(0, 2*pi, 0.02)]
        g = list_plot(list(zip(x_coords, y_coords)))
        sphinx_plot(g)

    If instead you try to pass the two lists as separate arguments,
    you will get an error message::

        sage: list_plot(x_coords, y_coords)
        Traceback (most recent call last):
        ...
        TypeError: The second argument 'plotjoined' should be boolean (True or False).
        If you meant to plot two lists 'x' and 'y' against each other,
        use 'list_plot(list(zip(x,y)))'.

    Dictionaries with numeric keys and values can be plotted::

        sage: list_plot({22: 3365, 27: 3295, 37: 3135, 42: 3020, 47: 2880, 52: 2735, 57: 2550})
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = list_plot({22: 3365, 27: 3295, 37: 3135, 42: 3020, 47: 2880, 52: 2735, 57: 2550})
        sphinx_plot(g)

    Plotting in logarithmic scale is possible for 2D list plots.
    There are two different syntaxes available::

        sage: yl = [2**k for k in range(20)]
        sage: list_plot(yl, scale='semilogy')  # long time  # log axis on vertical
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        yl = [2**k for k in range(20)]
        g = list_plot(yl, scale='semilogy')  # long time  # log axis on vertical
        sphinx_plot(g)

    ::

        sage: list_plot_semilogy(yl)       # same
        Graphics object consisting of 1 graphics primitive

    .. warning::

        If ``plotjoined`` is ``False`` then the axis that is in log scale
        must have all points strictly positive. For instance, the following
        plot will show no points in the figure since the points in the
        horizontal axis starts from `(0,1)`. Further, matplotlib will display
        a user warning.

        ::

            sage: list_plot(yl, scale='loglog')         # both axes are log
            doctest:warning
            ...
            Graphics object consisting of 1 graphics primitive

        Instead this will work. We drop the point `(0,1)`.::

            sage: list_plot(list(zip(range(1,len(yl)), yl[1:])), scale='loglog')  # long time
            Graphics object consisting of 1 graphics primitive

    We use :func:`list_plot_loglog` and plot in a different base.::

        sage: list_plot_loglog(list(zip(range(1,len(yl)), yl[1:])), base=2)  # long time
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        yl = [2**k for k in range(20)]
        g = list_plot_loglog(list(zip(range(1,len(yl)), yl[1:])), base=2) # long time
        sphinx_plot(g)

    We can also change the scale of the axes in the graphics just before
    displaying::

        sage: G = list_plot(yl) # long time
        sage: G.show(scale=('semilogy', 2)) # long time

    TESTS:

    We check to see whether elements of the Symbolic Ring are properly
    handled; see :issue:`16378` ::

        sage: list_plot([1+I, 2+I])
        Graphics object consisting of 1 graphics primitive

        sage: list_plot([1+I, 2, CC(3+I)])
        Graphics object consisting of 1 graphics primitive

        sage: list_plot([2, SR(1), CC(1+i)])
        Graphics object consisting of 1 graphics primitive

    We check to see that the x/y min/max data are set correctly::

        sage: d = list_plot([(100,100), (120, 120)]).get_minmax_data()
        sage: d['xmin']
        100.0
        sage: d['ymin']
        100.0
    """
    from sage.plot.all import point
    try:
        if not data:
            return Graphics()
    except ValueError:  # numpy raises ValueError if it is not empty
        pass
    if not isinstance(plotjoined, bool):
        raise TypeError("The second argument 'plotjoined' should be boolean "
                    "(True or False).  If you meant to plot two lists 'x' "
                    "and 'y' against each other, use 'list_plot(list(zip(x,y)))'.")
    if isinstance(data, dict):
        if plotjoined:
            list_data = sorted(data.items())
        else:
            list_data = list(data.items())
        return list_plot(list_data, plotjoined=plotjoined, **kwargs)
    try:
        from sage.rings.real_double import RDF
        RDF(data[0])
        data = list(enumerate(data))
    except TypeError: # we can get this TypeError if the element is a list
                      # or tuple or numpy array, or an element of CC, CDF
        # We also want to avoid doing CC(data[0]) here since it will go
        # through if data[0] is really a tuple and every element of the
        # data will be converted to a complex and later converted back to
        # a tuple.
        # So, the only other check we need to do is whether data[0] is an
        # element of the Symbolic Ring.
        if isinstance(data[0], Expression):
            data = list(enumerate(data))

    try:
        if plotjoined:
            return line(data, **kwargs)
        else:
            return point(data, **kwargs)
    except (TypeError, IndexError):
        # Assume we have complex-valued input and plot real and imaginary parts.
        # Need to catch IndexError because if data is, say, [(0, 1), (1, I)],
        # point3d() throws an IndexError on the (0,1) before it ever
        # gets to (1, I).
        from sage.rings.cc import CC
        # if we get here, we already did "list(enumerate(data))",
        # so look at z[1] in inner list
        data = [(z.real(), z.imag()) for z in [CC(z[1]) for z in data]]
        if plotjoined:
            return line(data, **kwargs)
        else:
            return point(data, **kwargs)

#------------------------ Graphs on log scale ---------------------------#


@options(base=10)
def plot_loglog(funcs, *args, **kwds):
    """
    Plot graphics in 'loglog' scale, that is, both the horizontal and the
    vertical axes will be in logarithmic scale.

    INPUT:

    - ``base`` -- (default: `10`) the base of the logarithm; this must be
      greater than 1. The base can be also given as a list or tuple
      ``(basex, basey)``.  ``basex`` sets the base of the logarithm along the
      horizontal axis and ``basey`` sets the base along the vertical axis.

    - ``funcs`` -- any Sage object which is acceptable to the :func:`plot`

    For all other inputs, look at the documentation of :func:`plot`.

    EXAMPLES::

        sage: plot_loglog(exp, (1,10))  # plot in loglog scale with base 10
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = plot_loglog(exp, (1,10)) # plot in loglog scale with base 10
        sphinx_plot(g)

    ::

        sage: plot_loglog(exp, (1,10), base=2.1)  # with base 2.1 on both axes  # long time
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = plot_loglog(exp, (1,10), base=2.1)
        sphinx_plot(g)

    ::

        sage: plot_loglog(exp, (1,10), base=(2,3))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = plot_loglog(exp, (1,10), base=(2,3))
        sphinx_plot(g)
    """
    return plot(funcs, *args, scale='loglog', **kwds)


@options(base=10)
def plot_semilogx(funcs, *args, **kwds):
    """
    Plot graphics in 'semilogx' scale, that is, the horizontal axis will be
    in logarithmic scale.

    INPUT:

    - ``base`` -- (default: `10`) the base of the logarithm; this must be
      greater than 1

    - ``funcs`` -- any Sage object which is acceptable to the :func:`plot`

    For all other inputs, look at the documentation of :func:`plot`.

    EXAMPLES::

        sage: plot_semilogx(exp, (1,10))  # plot in semilogx scale, base 10     # long time
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = plot_semilogx(exp, (1,10))
        sphinx_plot(g)

    ::

        sage: plot_semilogx(exp, (1,10), base=2)  # with base 2
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = plot_semilogx(exp, (1,10), base=2)
        sphinx_plot(g)

    ::

        sage: s = var('s')  # Samples points logarithmically so graph is smooth
        sage: f = 4000000/(4000000 + 4000*s*i - s*s)
        sage: plot_semilogx(20*log(abs(f), 10), (s, 1, 1e6))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        s = var('s') # Samples points logarithmically so graph is smooth
        f = 4000000/(4000000 + 4000*s*i - s*s)
        g = plot_semilogx(20*log(abs(f), 10), (s, 1, 1e6))
        sphinx_plot(g)
    """
    return plot(funcs, *args, scale='semilogx', **kwds)


@options(base=10)
def plot_semilogy(funcs, *args, **kwds):
    """
    Plot graphics in 'semilogy' scale, that is, the vertical axis will be
    in logarithmic scale.

    INPUT:

    - ``base`` -- (default: `10`) the base of the logarithm; this must be
      greater than 1

    - ``funcs`` -- any Sage object which is acceptable to the :func:`plot`

    For all other inputs, look at the documentation of :func:`plot`.

    EXAMPLES::

        sage: plot_semilogy(exp, (1, 10))  # long time # plot in semilogy scale, base 10
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = plot_semilogy(exp, (1,10))  # long time # plot in semilogy scale, base 10
        sphinx_plot(g)

    ::

        sage: plot_semilogy(exp, (1, 10), base=2)  # long time # with base 2
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = plot_semilogy(exp, (1,10), base=2)  # long time # with base 2
        sphinx_plot(g)
    """
    return plot(funcs, *args, scale='semilogy', **kwds)


@options(base=10)
def list_plot_loglog(data, plotjoined=False, **kwds):
    """
    Plot the ``data`` in 'loglog' scale, that is, both the horizontal and the
    vertical axes will be in logarithmic scale.

    INPUT:

    - ``base`` -- (default: `10`) the base of the logarithm; this must be
      greater than 1. The base can be also given as a list or tuple
      ``(basex, basey)``.  ``basex`` sets the base of the logarithm along the
      horizontal axis and ``basey`` sets the base along the vertical axis.

    For all other inputs, look at the documentation of :func:`list_plot`.

    EXAMPLES::

        sage: yl = [5**k for k in range(10)]; xl = [2**k for k in range(10)]
        sage: list_plot_loglog(list(zip(xl, yl)))  # use loglog scale with base 10      # long time
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        yl = [5**k for k in range(10)]
        xl = [2**k for k in range(10)]
        g = list_plot_loglog(list(zip(xl, yl)))
        sphinx_plot(g)

    ::

        sage: list_plot_loglog(list(zip(xl, yl)),  # with base 2.1 on both axes         # long time
        ....:                  base=2.1)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        yl = [5**k for k in range(10)]
        xl = [2**k for k in range(10)]
        g = list_plot_loglog(list(zip(xl, yl)), base=2.1)
        sphinx_plot(g)

    ::

        sage: list_plot_loglog(list(zip(xl, yl)), base=(2,5))                           # long time
        Graphics object consisting of 1 graphics primitive

    .. warning::

        If ``plotjoined`` is ``False`` then the axis that is in log scale
        must have all points strictly positive. For instance, the following
        plot will show no points in the figure since the points in the
        horizontal axis starts from `(0,1)`.

        ::

            sage: yl = [2**k for k in range(20)]
            sage: list_plot_loglog(yl)
            Graphics object consisting of 1 graphics primitive

        Instead this will work. We drop the point `(0,1)`.::

            sage: list_plot_loglog(list(zip(range(1,len(yl)), yl[1:])))
            Graphics object consisting of 1 graphics primitive
    """
    return list_plot(data, plotjoined=plotjoined, scale='loglog', **kwds)


@options(base=10)
def list_plot_semilogx(data, plotjoined=False, **kwds):
    """
    Plot ``data`` in 'semilogx' scale, that is, the horizontal axis will be
    in logarithmic scale.

    INPUT:

    - ``base`` -- (default: `10`) the base of the logarithm; this must be
      greater than 1

    For all other inputs, look at the documentation of :func:`list_plot`.

    EXAMPLES::

        sage: yl = [2**k for k in range(12)]
        sage: list_plot_semilogx(list(zip(yl,yl)))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        yl = [2**k for k in range(12)]
        g = list_plot_semilogx(list(zip(yl,yl)))
        sphinx_plot(g)

    .. warning::

        If ``plotjoined`` is ``False`` then the horizontal axis must have all
        points strictly positive. Otherwise the plot will come up empty.
        For instance the following plot contains a point at `(0,1)`.

        ::

            sage: yl = [2**k for k in range(12)]
            sage: list_plot_semilogx(yl)  # plot empty due to (0,1)
            Graphics object consisting of 1 graphics primitive

        We remove `(0,1)` to fix this.::

            sage: list_plot_semilogx(list(zip(range(1, len(yl)), yl[1:])))
            Graphics object consisting of 1 graphics primitive

    ::

        sage: list_plot_semilogx([(1,2),(3,4),(3,-1),(25,3)], base=2)  # with base 2
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = list_plot_semilogx([(1,2),(3,4),(3,-1),(25,3)], base=2)
        sphinx_plot(g)
    """
    return list_plot(data, plotjoined=plotjoined, scale='semilogx', **kwds)


@options(base=10)
def list_plot_semilogy(data, plotjoined=False, **kwds):
    """
    Plot ``data`` in 'semilogy' scale, that is, the vertical axis will be
    in logarithmic scale.

    INPUT:

    - ``base`` -- (default: `10`) the base of the logarithm; this must be
      greater than 1

    For all other inputs, look at the documentation of :func:`list_plot`.

    EXAMPLES::

        sage: yl = [2**k for k in range(12)]
        sage: list_plot_semilogy(yl)  # plot in semilogy scale, base 10
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        yl = [2**k for k in range(12)]
        g = list_plot_semilogy(yl)  # plot in semilogy scale, base 10
        sphinx_plot(g)

    .. warning::

        If ``plotjoined`` is ``False`` then the vertical axis must have all
        points strictly positive. Otherwise the plot will come up empty.
        For instance the following plot contains a point at `(1,0)`. Further,
        matplotlib will display a user warning.

        ::

            sage: xl = [2**k for k in range(12)]; yl = range(len(xl))
            sage: list_plot_semilogy(list(zip(xl, yl)))  # plot empty due to (1,0)
            doctest:warning
            ...
            Graphics object consisting of 1 graphics primitive

        We remove `(1,0)` to fix this.::

            sage: list_plot_semilogy(list(zip(xl[1:],yl[1:])))
            Graphics object consisting of 1 graphics primitive

    ::

        sage: list_plot_semilogy([2, 4, 6, 8, 16, 31], base=2)  # with base 2
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        g = list_plot_semilogy([2, 4, 6, 8, 16, 31], base=2)
        sphinx_plot(g)
    """
    return list_plot(data, plotjoined=plotjoined, scale='semilogy', **kwds)


def to_float_list(v):
    """
    Given a list or tuple or iterable v, coerce each element of v to a
    float and make a list out of the result.

    EXAMPLES::

        sage: from sage.plot.plot import to_float_list
        sage: to_float_list([1,1/2,3])
        [1.0, 0.5, 3.0]
    """
    return [float(x) for x in v]


def reshape(v, n, m):
    """
    Helper function for creating graphics arrays.

    The input array is flattened and turned into an `n\times m`
    array, with blank graphics object padded at the end, if
    necessary.

    INPUT:

    - ``v`` -- list of lists or tuples

    - ``n``, ``m`` -- integers

    OUTPUT: list of lists of graphics objects

    EXAMPLES::

        sage: L = [plot(sin(k*x), (x,-pi,pi)) for k in range(10)]
        sage: graphics_array(L,3,4)  # long time (up to 4s on sage.math, 2012)
        Graphics Array of size 3 x 4

    ::

        sage: M = [[plot(sin(k*x), (x,-pi,pi)) for k in range(3)],
        ....:      [plot(cos(j*x), (x,-pi,pi)) for j in [3..5]]]
        sage: graphics_array(M,6,1)  # long time (up to 4s on sage.math, 2012)
        Graphics Array of size 6 x 1

    TESTS::

        sage: L = [plot(sin(k*x), (x,-pi,pi)) for k in [1..3]]
        sage: graphics_array(L,0,-1) # indirect doctest
        Traceback (most recent call last):
        ...
        ValueError: array sizes must be positive
    """
    if not (n > 0 and m > 0):
        raise ValueError('array sizes must be positive')
    G = Graphics()
    G.axes(False)
    if len(v) == 0:
        return [[G]*m]*n

    if not isinstance(v[0], Graphics):
        # a list of lists -- flatten it
        v = sum([list(x) for x in v], [])

    # Now v should be a single list.
    # First, make it have the right length.
    v = list(v)   # do not mutate the argument
    for i in range(n * m - len(v)):
        v.append(G)

    # Next, create a list of lists out of it.
    L = []
    k = 0
    for i in range(n):
        w = []
        for j in range(m):
            w.append(v[k])
            k += 1
        L.append(w)

    return L


def graphics_array(array, nrows=None, ncols=None):
    r"""
    Plot a list of lists (or tuples) of graphics objects on one canvas,
    arranged as an array.

    INPUT:

    - ``array`` -- either a list of lists of
      :class:`~sage.plot.graphics.Graphics` elements or a
      single list of :class:`~sage.plot.graphics.Graphics` elements

    - ``nrows``, ``ncols`` -- (optional) integers. If both are given then
      the input array is flattened and turned into an ``nrows`` x
      ``ncols`` array, with blank graphics objects padded at the end,
      if necessary. If only one is specified, the other is chosen
      automatically.

    OUTPUT: an instance of :class:`~sage.plot.multigraphics.GraphicsArray`

    EXAMPLES:

    Make some plots of `\sin` functions::

        sage: # long time
        sage: f(x) = sin(x)
        sage: g(x) = sin(2*x)
        sage: h(x) = sin(4*x)
        sage: p1 = plot(f, (-2*pi,2*pi), color=hue(0.5))
        sage: p2 = plot(g, (-2*pi,2*pi), color=hue(0.9))
        sage: p3 = parametric_plot((f,g), (0,2*pi), color=hue(0.6))
        sage: p4 = parametric_plot((f,h), (0,2*pi), color=hue(1.0))

    Now make a graphics array out of the plots::

        sage: graphics_array(((p1,p2), (p3,p4)))        # long time
        Graphics Array of size 2 x 2

    .. PLOT::

        def f(x): return sin(x)
        def g(x): return sin(2*x)
        def h(x): return sin(4*x)
        p1 = plot(f, (-2*pi,2*pi), color=hue(0.5))
        p2 = plot(g, (-2*pi,2*pi), color=hue(0.9))
        p3 = parametric_plot((f,g), (0,2*pi), color=hue(0.6))
        p4 = parametric_plot((f,h), (0,2*pi), color=hue(1.0))
        g = graphics_array(((p1, p2), (p3, p4)))
        sphinx_plot(g)

    One can also name the array, and then use
    :meth:`~sage.plot.multigraphics.MultiGraphics.show`
    or :meth:`~sage.plot.multigraphics.MultiGraphics.save`::

        sage: ga = graphics_array(((p1,p2), (p3,p4)))   # long time
        sage: ga.show() # long time; same output as above

    Here we give only one row::

        sage: p1 = plot(sin, (-4,4))
        sage: p2 = plot(cos, (-4,4))
        sage: ga = graphics_array([p1, p2]); ga
        Graphics Array of size 1 x 2
        sage: ga.show()

    .. PLOT::

        p1 = plot(sin,(-4,4))
        p2 = plot(cos,(-4,4))
        ga = graphics_array([p1, p2])
        sphinx_plot(ga)

    It is possible to use ``figsize`` to change the size of the plot
    as a whole::

        sage: x = var('x')
        sage: L = [plot(sin(k*x), (x,-pi,pi)) for k in [1..3]]
        sage: ga = graphics_array(L)
        sage: ga.show(figsize=[5,3])  # smallish and compact

    .. PLOT::

        ga = graphics_array([plot(sin(k*x), (x,-pi,pi)) for k in range(1, 4)])
        sphinx_plot(ga, figsize=[5,3])

    ::

        sage: ga.show(figsize=[5,7])  # tall and thin; long time

    .. PLOT::

        ga = graphics_array([plot(sin(k*x), (x,-pi,pi)) for k in range(1, 4)])
        sphinx_plot(ga, figsize=[5,7])

    ::

        sage: ga.show(figsize=4)  # width=4 inches, height fixed from default aspect ratio

    .. PLOT::

        ga = graphics_array([plot(sin(k*x), (x,-pi,pi)) for k in range(1, 4)])
        sphinx_plot(ga, figsize=4)

    Specifying only the number of rows or the number of columns
    computes the other dimension automatically::

        sage: ga = graphics_array([plot(sin)] * 10, nrows=3)
        sage: ga.nrows(), ga.ncols()
        (3, 4)
        sage: ga = graphics_array([plot(sin)] * 10, ncols=3)
        sage: ga.nrows(), ga.ncols()
        (4, 3)
        sage: ga = graphics_array([plot(sin)] * 4, nrows=2)
        sage: ga.nrows(), ga.ncols()
        (2, 2)
        sage: ga = graphics_array([plot(sin)] * 6, ncols=2)
        sage: ga.nrows(), ga.ncols()
        (3, 2)

    The options like ``fontsize``, ``scale`` or ``frame`` passed to individual
    plots are preserved::

        sage: p1 = plot(sin(x^2), (x, 0, 6),
        ....:           axes_labels=[r'$\theta$', r'$\sin(\theta^2)$'], fontsize=16)
        sage: p2 = plot(x^3, (x, 1, 100), axes_labels=[r'$x$', r'$y$'],
        ....:           scale='semilogy', frame=True, gridlines='minor')
        sage: ga = graphics_array([p1, p2])
        sage: ga.show()

    .. PLOT::

        p1 = plot(sin(x**2), (x, 0, 6), \
                  axes_labels=[r'$\theta$', r'$\sin(\theta^2)$'], fontsize=16)
        p2 = plot(x**3, (x, 1, 100), axes_labels=[r'$x$', r'$y$'], \
                  scale='semilogy', frame=True, gridlines='minor')
        ga = graphics_array([p1, p2])
        sphinx_plot(ga)

    .. SEEALSO::

        :class:`~sage.plot.multigraphics.GraphicsArray` for more examples
    """
    # TODO: refactor the whole array flattening and reshaping into a class
    if nrows is None and ncols is None:
        pass
    elif nrows is not None and ncols is not None:
        nrows = int(nrows)
        ncols = int(ncols)
        array = reshape(array, nrows, ncols)
    else:
        # nrows is None xor ncols is None
        if len(array) > 0 and isinstance(array[0], Graphics):
            length = len(array)
        else:
            length = sum(map(len, array))
        if nrows is None:
            ncols = int(ncols)
            nrows = length // ncols
            if nrows*ncols < length or nrows == 0:
                nrows += 1
        elif ncols is None:
            nrows = int(nrows)
            ncols = length // nrows
            if nrows*ncols < length or ncols == 0:
                ncols += 1
        else:
            assert False
        array = reshape(array, nrows, ncols)
    return GraphicsArray(array)


def multi_graphics(graphics_list):
    r"""
    Plot a list of graphics at specified positions on a single canvas.

    If the graphics positions define a regular array, use
    :func:`graphics_array` instead.

    INPUT:

    - ``graphics_list`` -- list of graphics along with their
      positions on the canvas; each element of ``graphics_list`` is either

      - a pair ``(graphics, position)``, where ``graphics`` is a
        :class:`~sage.plot.graphics.Graphics` object and ``position`` is
        the 4-tuple ``(left, bottom, width, height)`` specifying the location
        and size of the graphics on the canvas, all quantities being in
        fractions of the canvas width and height

      - or a single :class:`~sage.plot.graphics.Graphics` object; its position
        is then assumed to occupy the whole canvas, except for some padding;
        this corresponds to the default position
        ``(left, bottom, width, height) = (0.125, 0.11, 0.775, 0.77)``

    OUTPUT: an instance of :class:`~sage.plot.multigraphics.MultiGraphics`

    EXAMPLES:

    ``multi_graphics`` is to be used for plot arrangements that cannot be
    achieved with :func:`graphics_array`, for instance::

        sage: g1 = plot(sin(x), (x, -10, 10), frame=True)
        sage: g2 = EllipticCurve([0,0,1,-1,0]).plot(color='red', thickness=2,
        ....:                    axes_labels=['$x$', '$y$']) \
        ....:      + text(r"$y^2 + y = x^3 - x$", (1.2, 2), color='red')
        sage: g3 = matrix_plot(matrix([[1,3,5,1], [2,4,5,6], [1,3,5,7]]))
        sage: G = multi_graphics([(g1, (0.125, 0.65, 0.775, 0.3)),
        ....:                     (g2, (0.125, 0.11, 0.4, 0.4)),
        ....:                     (g3, (0.55, 0.18, 0.4, 0.3))])
        sage: G
        Multigraphics with 3 elements

    .. PLOT::

        g1 = plot(sin(x), (x, -10, 10), frame=True)
        g2 = EllipticCurve([0,0,1,-1,0]).plot(color='red', thickness=2, \
                           axes_labels=['$x$', '$y$']) \
             + text(r"$y^2 + y = x^3 - x$", (1.2, 2), color='red')
        g3 = matrix_plot(matrix([[1,3,5,1], [2,4,5,6], [1,3,5,7]]))
        G = multi_graphics([(g1, (0.125, 0.65, 0.775, 0.3)), \
                            (g2, (0.125, 0.11, 0.4, 0.4)), \
                            (g3, (0.55, 0.18, 0.4, 0.3))])
        sphinx_plot(G)

    An example with a list containing a graphics object without any specified
    position (the graphics, here ``g3``, occupies then the whole canvas)::

        sage: G = multi_graphics([g3, (g1, (0.4, 0.4, 0.2, 0.2))])
        sage: G
        Multigraphics with 2 elements

    .. PLOT::

        g1 = plot(sin(x), (x, -10, 10), frame=True)
        g3 = matrix_plot(matrix([[1,3,5,1], [2,4,5,6], [1,3,5,7]]))
        G = multi_graphics([g3, (g1, (0.4, 0.4, 0.2, 0.2))])
        sphinx_plot(G)

    .. SEEALSO::

        :class:`~sage.plot.multigraphics.MultiGraphics` for more examples
    """
    return MultiGraphics(graphics_list)


def minmax_data(xdata, ydata, dict=False):
    """
    Return the minimums and maximums of ``xdata`` and ``ydata``.

    If dict is False, then minmax_data returns the tuple (xmin, xmax,
    ymin, ymax); otherwise, it returns a dictionary whose keys are
    'xmin', 'xmax', 'ymin', and 'ymax' and whose values are the
    corresponding values.

    EXAMPLES::

        sage: from sage.plot.plot import minmax_data
        sage: minmax_data([], [])
        (-1, 1, -1, 1)
        sage: minmax_data([-1, 2], [4, -3])
        (-1, 2, -3, 4)
        sage: minmax_data([1, 2], [4, -3])
        (1, 2, -3, 4)
        sage: d = minmax_data([-1, 2], [4, -3], dict=True)
        sage: list(sorted(d.items()))
        [('xmax', 2), ('xmin', -1), ('ymax', 4), ('ymin', -3)]
        sage: d = minmax_data([1, 2], [3, 4], dict=True)
        sage: list(sorted(d.items()))
        [('xmax', 2), ('xmin', 1), ('ymax', 4), ('ymin', 3)]
    """
    xmin = min(xdata) if len(xdata) else -1
    xmax = max(xdata) if len(xdata) else 1
    ymin = min(ydata) if len(ydata) else -1
    ymax = max(ydata) if len(ydata) else 1
    if dict:
        return {'xmin': xmin, 'xmax': xmax,
                'ymin': ymin, 'ymax': ymax}
    else:
        return xmin, xmax, ymin, ymax


def adaptive_refinement(f, p1, p2, adaptive_tolerance=0.01,
                        adaptive_recursion=5, level=0, *, excluded=False):
    r"""
    The adaptive refinement algorithm for plotting a function ``f``. See
    the docstring for plot for a description of the algorithm.

    INPUT:

    - ``f`` -- a function of one variable

    - ``p1``, ``p2`` -- two points to refine between

    - ``adaptive_recursion`` -- (default: `5`) how many
      levels of recursion to go before giving up when doing adaptive
      refinement. Setting this to 0 disables adaptive refinement.

    - ``adaptive_tolerance`` -- (default: `0.01`) how large
      a relative difference should be before the adaptive refinement
      code considers it significant; see documentation for generate_plot_points
      for more information.  See the documentation for :func:`plot` for more
      information on how the adaptive refinement algorithm works.

    - ``excluded`` -- (default: ``False``) also return locations where it has been
      discovered that the function is not defined
      (y-value will be ``'NaN'`` in this case)

    OUTPUT:

    A list of points to insert between ``p1`` and
    ``p2`` to get a better linear approximation between them.
    If ``excluded``, also x-values for which the calculation failed are given
    with ``'NaN'`` as y-value.

    TESTS::

        sage: from sage.plot.plot import adaptive_refinement
        sage: adaptive_refinement(sin, (0,0), (pi,0), adaptive_tolerance=0.01,
        ....:                     adaptive_recursion=0)
        []
        sage: adaptive_refinement(sin, (0,0), (pi,0), adaptive_tolerance=0.01)
        [(0.125*pi, 0.3826834323650898), (0.1875*pi, 0.5555702330196022),
         (0.25*pi, 0.7071067811865475),  (0.3125*pi, 0.8314696123025452),
         (0.375*pi, 0.9238795325112867), (0.4375*pi, 0.9807852804032304),
         (0.5*pi, 1.0), (0.5625*pi, 0.9807852804032304), (0.625*pi, 0.9238795325112867),
         (0.6875*pi, 0.8314696123025455), (0.75*pi, 0.7071067811865476),
         (0.8125*pi, 0.5555702330196022), (0.875*pi, 0.3826834323650899)]

    This shows that lowering ``adaptive_tolerance`` and raising
    ``adaptive_recursion`` both increase the number of subdivision
    points, though which one creates more points is heavily
    dependent upon the function being plotted.

    ::

        sage: x = var('x')
        sage: f(x) = sin(1/x)
        sage: n1 = len(adaptive_refinement(f, (0,0), (pi,0), adaptive_tolerance=0.01)); n1
        15
        sage: n2 = len(adaptive_refinement(f, (0,0), (pi,0), adaptive_recursion=10,
        ....:                              adaptive_tolerance=0.01)); n2
        79
        sage: n3 = len(adaptive_refinement(f, (0,0), (pi,0), adaptive_tolerance=0.001)); n3
        26

    Exclusion points will be added if ``excluded`` is set::

        sage: f(x) = 1/x
        sage: adaptive_refinement(f, (-1, -1), (3, 1/3), adaptive_recursion=2, excluded=False)
        [(1.0, 1.0), (2.0, 0.5)]
        sage: adaptive_refinement(f, (-1, -1), (3, 1/3), adaptive_recursion=2, excluded=True)
        [(0.0, 'NaN'), (1.0, 1.0), (2.0, 0.5)]
    """
    if level >= adaptive_recursion:
        return []

    x = (p1[0] + p2[0])/2.0
    msg = ''

    try:
        y = float(f(x))
        if str(y) in ['nan', 'NaN', 'inf', '-inf']:
            sage.misc.verbose.verbose(f"{msg}\nUnable to compute f({x})",1)
            # give up for this branch
            if excluded:
                return [(x, 'NaN')]
            return []

    except (ZeroDivisionError, TypeError, ValueError, OverflowError) as msg:
        sage.misc.verbose.verbose(f"{msg}\nUnable to compute f({x})", 1)
        # give up for this branch
        if excluded:
            return [(x, 'NaN')]
        return []

    # this distance calculation is not perfect.
    if abs((p1[1] + p2[1])/2.0 - y) > adaptive_tolerance:
        return adaptive_refinement(f, p1, (x, y),
                    adaptive_tolerance=adaptive_tolerance,
                    adaptive_recursion=adaptive_recursion,
                    level=level+1,
                    excluded=excluded) \
                    + [(x, y)] + \
            adaptive_refinement(f, (x, y), p2,
                    adaptive_tolerance=adaptive_tolerance,
                    adaptive_recursion=adaptive_recursion,
                    level=level+1,
                    excluded=excluded)
    else:
        return []


def generate_plot_points(f, xrange, plot_points=5, adaptive_tolerance=0.01,
                         adaptive_recursion=5, randomize=True,
                         initial_points=None, *, excluded=False,
                         imaginary_tolerance=1e-8):
    r"""
    Calculate plot points for a function f in the interval xrange.  The
    adaptive refinement algorithm is also automatically invoked with a
    *relative* adaptive tolerance of adaptive_tolerance; see below.

    INPUT:

    - ``f`` -- a function of one variable

    - ``p1``, ``p2`` -- two points to refine between

    - ``plot_points`` -- (default: 5) the minimal number of plot points. (Note
      however that in any actual plot a number is passed to this, with default
      value 200.)

    - ``adaptive_recursion`` -- (default: 5) how many levels of recursion to go
      before giving up when doing adaptive refinement.  Setting this to 0
      disables adaptive refinement.

    - ``adaptive_tolerance`` -- (default: 0.01) how large the relative difference
      should be before the adaptive refinement code considers it significant.  If
      the actual difference is greater than adaptive_tolerance*delta, where delta
      is the initial subinterval size for the given xrange and plot_points, then
      the algorithm will consider it significant.

    - ``initial_points`` -- (default: ``None``) a list of x-values that should be evaluated

    - ``excluded`` -- (default: ``False``) add a list of discovered x-values, for
      which ``f`` is not defined

    - ``imaginary_tolerance`` -- (default: ``1e-8``) if an imaginary
      number arises (due, for example, to numerical issues), this
      tolerance specifies how large it has to be in magnitude before
      we raise an error.  In other words, imaginary parts smaller than
      this are ignored in your plot points.

    OUTPUT:

    - a list of points (x, f(x)) in the interval xrange, which approximate
      the function f.

    - if ``excluded`` a tuple consisting of the above and a list of x-values
      at which ``f`` is not defined

    TESTS::

        sage: from sage.plot.plot import generate_plot_points
        sage: generate_plot_points(sin, (0, pi), plot_points=2, adaptive_recursion=0)
        [(0.0, 0.0), (3.141592653589793, 1.2246...e-16)]

        sage: from sage.plot.plot import generate_plot_points
        sage: generate_plot_points(lambda x: x^2, (0, 6), plot_points=2, adaptive_recursion=0, initial_points=[1,2,3])
        [(0.0, 0.0), (1.0, 1.0), (2.0, 4.0), (3.0, 9.0), (6.0, 36.0)]

    The delta remains consistent with ``randomize=False`` and no
    adaptive recursion::

        sage: pps = generate_plot_points(sin(x).function(x),
        ....:                            (-pi, pi),
        ....:                            randomize=False,
        ....:                            adaptive_recursion=0)
        sage: [pps[k][0] - pps[k-1][0] for k in range(1,len(pps))] # abs tol 1e-10
        [1.5707963267948966,
         1.5707963267948966,
         1.5707963267948966,
         1.5707963267948966]

    This shows that lowering adaptive_tolerance and raising
    adaptive_recursion both increase the number of subdivision points.
    (Note that which creates more points is heavily dependent on the
    particular function plotted.)

    ::

        sage: x = var('x')
        sage: f(x) = sin(1/x)
        sage: [len(generate_plot_points(f, (-pi, pi), plot_points=16, adaptive_tolerance=i, randomize=False)) for i in [0.01, 0.001, 0.0001]]
        [97, 161, 275]

        sage: [len(generate_plot_points(f, (-pi, pi), plot_points=16, adaptive_recursion=i, randomize=False)) for i in [5, 10, 15]]
        [97, 499, 2681]

    Excluded x-values will be added, if ``exclusion`` is set::

        sage: generate_plot_points(log, (0, 1), plot_points=2, adaptive_recursion=0)
        [(1.0, 0.0)]
        sage: generate_plot_points(log, (0, 1), plot_points=2, adaptive_recursion=0, excluded=True)
        ([(1.0, 0.0)], [0.0])
    """
    from sage.plot.misc import setup_for_eval_on_grid
    f, ranges = setup_for_eval_on_grid(f,
                                       [xrange],
                                       plot_points,
                                       imaginary_tolerance=imaginary_tolerance)
    xmin, xmax, delta = ranges[0]
    x_values = srange(*ranges[0], include_endpoint=True)

    random = current_randstate().python_random().random

    for i in range(len(x_values)):
        xi = x_values[i]
        # Slightly randomize the interior sample points if
        # randomize is true
        if randomize and i > 0 and i < plot_points-1:
            xi += delta*(random() - 0.5)
            x_values[i] = xi

    # add initial points
    if isinstance(initial_points, list):
        x_values = sorted(x_values + initial_points)

    data = [None]*len(x_values)

    exceptions = 0
    exception_indices = []
    for i in range(len(x_values)):
        xi = x_values[i]

        try:
            data[i] = (float(xi), float(f(xi)))
            if str(data[i][1]) in ['nan', 'NaN', 'inf', '-inf']:
                msg = "Unable to compute f(%s)" % xi
                sage.misc.verbose.verbose(msg, 1)
                exceptions += 1
                exception_indices.append(i)

        except (ArithmeticError, TypeError, ValueError) as m:
            sage.misc.verbose.verbose(f"{m}\nUnable to compute f({xi})", 1)

            if i == 0: # Given an error for left endpoint, try to move it in slightly
                for j in range(1, 99):
                    xj = xi + delta*j/100.0
                    try:
                        data[i] = (float(xj), float(f(xj)))
                        # nan != nan
                        if data[i][1] != data[i][1]:
                            continue
                        break
                    except (ArithmeticError, TypeError, ValueError):
                        pass
                else:
                    msg = m
                    exceptions += 1
                    exception_indices.append(i)

            elif i == plot_points-1: # Given an error for right endpoint, try to move it in slightly
                for j in range(1, 99):
                    xj = xi - delta*j/100.0
                    try:
                        data[i] = (float(xj), float(f(xj)))
                        # nan != nan
                        if data[i][1] != data[i][1]:
                            continue
                        break
                    except (ArithmeticError, TypeError, ValueError):
                        pass
                else:
                    msg = m
                    exceptions += 1
                    exception_indices.append(i)
            else:
                msg = m
                exceptions += 1
                exception_indices.append(i)

    data = [data[i] for i in range(len(data)) if i not in exception_indices]
    excluded_points = [x_values[i] for i in exception_indices]

    # calls adaptive refinement
    i, j = 0, 0
    adaptive_tolerance = delta * float(adaptive_tolerance)
    adaptive_recursion = int(adaptive_recursion)

    while i < len(data) - 1:
        for p in adaptive_refinement(f, data[i], data[i+1],
                                     adaptive_tolerance=adaptive_tolerance,
                                     adaptive_recursion=adaptive_recursion,
                                     excluded=True):
            if p[1] == "NaN":
                excluded_points.append(p[0])
            else:
                data.insert(i+1, p)
                i += 1
        i += 1

    if (len(data) == 0 and exceptions > 0) or exceptions > 10:
        sage.misc.verbose.verbose("WARNING: When plotting, failed to evaluate function at %s points." % exceptions, level=0)
        sage.misc.verbose.verbose("Last error message: '%s'" % msg, level=0)

    if excluded:
        return data, excluded_points
    return data
