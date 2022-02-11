"""
Line Plots
"""

#*****************************************************************************
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
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.plot.primitive import GraphicPrimitive_xydata
from sage.misc.decorators import options, rename_keyword
from sage.plot.colors import to_mpl_color

class Line(GraphicPrimitive_xydata):
    """
    Primitive class that initializes the line graphics type.

    EXAMPLES::

        sage: from sage.plot.line import Line
        sage: Line([1,2,7], [1,5,-1], {})
        Line defined by 3 points
    """
    def __init__(self, xdata, ydata, options):
        """
        Initialize a line graphics primitive.

        EXAMPLES::

            sage: from sage.plot.line import Line
            sage: Line([-1,2], [17,4], {'thickness':2})
            Line defined by 2 points
        """
        valid_options = self._allowed_options()
        for opt in options:
            if opt not in valid_options:
                raise RuntimeError("Error in line(): option '%s' not valid." % opt)
        self.xdata = xdata
        self.ydata = ydata
        GraphicPrimitive_xydata.__init__(self, options)

    def _allowed_options(self):
        """
        Displayed the list of allowed line options.

        EXAMPLES::

            sage: from sage.plot.line import Line
            sage: list(sorted(Line([-1,2], [17,4], {})._allowed_options().items()))
            [('alpha', 'How transparent the line is.'),
             ('hue', 'The color given as a hue.'),
             ('legend_color', 'The color of the legend text.'),
             ('legend_label', 'The label for this item in the legend.'),
             ('linestyle',
              "The style of the line, which is one of '--' (dashed), '-.' (dash dot), '-' (solid), 'steps', ':' (dotted)."),
             ('marker', 'the marker symbol (see documentation for line2d for details)'),
             ('markeredgecolor', 'the color of the marker edge'),
             ('markeredgewidth', 'the size of the marker edge in points'),
             ('markerfacecolor', 'the color of the marker face'),
             ('markersize', 'the size of the marker in points'),
             ('rgbcolor', 'The color as an RGB tuple.'),
             ('thickness', 'How thick the line is.'),
             ('zorder', 'The layer level in which to draw')]
        """
        return {'alpha':'How transparent the line is.',
                'legend_color':'The color of the legend text.',
                'legend_label':'The label for this item in the legend.',
                'thickness':'How thick the line is.',
                'rgbcolor':'The color as an RGB tuple.',
                'hue':'The color given as a hue.',
                'linestyle':"The style of the line, which is one of '--' (dashed), '-.' (dash dot), '-' (solid), 'steps', ':' (dotted).",
                'marker':"the marker symbol (see documentation for line2d for details)",
                'markersize':'the size of the marker in points',
                'markeredgecolor':'the color of the marker edge',
                'markeredgewidth':'the size of the marker edge in points',
                'markerfacecolor':'the color of the marker face',
                'zorder':'The layer level in which to draw'
                }

    def _plot3d_options(self, options=None):
        """
        Translate 2D plot options into 3D plot options.

        EXAMPLES::

            sage: L = line([(1,1), (1,2), (2,2), (2,1)], alpha=.5, thickness=10, zorder=2)
            sage: l=L[0]; l
            Line defined by 4 points
            sage: m=l.plot3d(z=2)
            sage: m.texture.opacity
            0.5
            sage: m.thickness
            10
            sage: L = line([(1,1), (1,2), (2,2), (2,1)], linestyle=":")
            sage: L.plot3d()
            Traceback (most recent call last):
            ...
            NotImplementedError: Invalid 3d line style: ':'
        """
        if options is None:
            options = dict(self.options())
        options_3d = {}
        if 'thickness' in options:
            options_3d['thickness'] = options['thickness']
            del options['thickness']
        if 'zorder' in options:
            del options['zorder']
        if 'linestyle' in options:
            if options['linestyle'] not in ('-', 'solid'):
                raise NotImplementedError("Invalid 3d line style: '%s'"%
                                          (options['linestyle']))
            del options['linestyle']
        options_3d.update(GraphicPrimitive_xydata._plot3d_options(self, options))
        return options_3d

    def plot3d(self, z=0, **kwds):
        """
        Plots a 2D line in 3D, with default height zero.

        EXAMPLES::

            sage: E = EllipticCurve('37a').plot(thickness=5).plot3d()
            sage: F = EllipticCurve('37a').plot(thickness=5).plot3d(z=2)
            sage: E + F  # long time (5s on sage.math, 2012)
            Graphics3d Object
        """
        from sage.plot.plot3d.shapes2 import line3d
        options = self._plot3d_options()
        options.update(kwds)
        return line3d([(x, y, z) for x, y in zip(self.xdata, self.ydata)], **options)

    def _repr_(self):
        """
        String representation of a line primitive.

        EXAMPLES::

            sage: from sage.plot.line import Line
            sage: Line([-1,2,3,3], [17,4,0,2], {})._repr_()
            'Line defined by 4 points'
        """
        return "Line defined by %s points"%len(self)

    def __getitem__(self, i):
        """
        Extract the i-th element of the line (which is stored as a list of points).

        INPUT:

        - ``i`` -- an integer between 0 and the number of points minus 1

        OUTPUT:

        A 2-tuple of floats.

        EXAMPLES::

            sage: L = line([(1,2), (3,-4), (2, 5), (1,2)])
            sage: line_primitive = L[0]; line_primitive
            Line defined by 4 points
            sage: line_primitive[0]
            (1.0, 2.0)
            sage: line_primitive[2]
            (2.0, 5.0)
            sage: list(line_primitive)
            [(1.0, 2.0), (3.0, -4.0), (2.0, 5.0), (1.0, 2.0)]
        """
        return self.xdata[i], self.ydata[i]

    def __setitem__(self, i, point):
        """
        Set the i-th element of this line (really a sequence of lines
        through given points).

        INPUT:

        - ``i`` -- an integer between 0 and the number of points on the
          line minus 1

        - ``point`` -- a 2-tuple of floats

        EXAMPLES:

        We create a line graphics object $L$ and get ahold of the
        corresponding line graphics primitive::

            sage: L = line([(1,2), (3,-4), (2, 5), (1,2)])
            sage: line_primitive = L[0]; line_primitive
            Line defined by 4 points

        We then set the 0th point to `(0,0)` instead of `(1,2)`::

            sage: line_primitive[0] = (0,0)
            sage: line_primitive[0]
            (0.0, 0.0)

        Plotting we visibly see the change -- now the line starts at `(0,0)`::

            sage: L
            Graphics object consisting of 1 graphics primitive
        """
        self.xdata[i] = float(point[0])
        self.ydata[i] = float(point[1])

    def __len__(self):
        r"""
        Return the number of points on this line (where a line is really a sequence
        of line segments through a given list of points).

        EXAMPLES:

        We create a line, then grab the line primitive as ``L[0]`` and compute
        its length::

            sage: L = line([(1,2), (3,-4), (2, 5), (1,2)])
            sage: len(L[0])
            4
        """
        return len(self.xdata)

    def _render_on_subplot(self, subplot):
        """
        Render this line on a matplotlib subplot.

        INPUT:

        - ``subplot`` -- a matplotlib subplot

        EXAMPLES:

        This implicitly calls this function::

            sage: line([(1,2), (3,-4), (2, 5), (1,2)])
            Graphics object consisting of 1 graphics primitive
        """
        import matplotlib.lines as lines
        options = dict(self.options())
        for o in ('alpha', 'legend_color', 'legend_label', 'linestyle',
                  'rgbcolor', 'thickness'):
            if o in options:
                del options[o]
        p = lines.Line2D(self.xdata, self.ydata, **options)
        options = self.options()
        a = float(options['alpha'])
        p.set_alpha(a)
        p.set_linewidth(float(options['thickness']))
        p.set_color(to_mpl_color(options['rgbcolor']))
        p.set_label(options['legend_label'])
        # we don't pass linestyle in directly since the drawstyles aren't
        # pulled off automatically.  This (I think) is a bug in matplotlib 1.0.1
        if 'linestyle' in options:
            from sage.plot.misc import get_matplotlib_linestyle
            p.set_linestyle(get_matplotlib_linestyle(options['linestyle'],
                                                     return_type='short'))
        subplot.add_line(p)

def line(points, **kwds):
    """
    Returns either a 2-dimensional or 3-dimensional line depending
    on value of points.

    INPUT:

    -  ``points`` - either a single point (as a tuple), a list of
       points, a single complex number, or a list of complex numbers.

    For information regarding additional arguments, see either line2d?
    or line3d?.

    EXAMPLES::

        sage: line([(0,0), (1,1)])
        Graphics object consisting of 1 graphics primitive

    ::

        sage: line([(0,0,1), (1,1,1)])
        Graphics3d Object
    """
    try:
        return line2d(points, **kwds)
    except ValueError:
        from sage.plot.plot3d.shapes2 import line3d
        return line3d(points, **kwds)


@rename_keyword(color='rgbcolor')
@options(alpha=1, rgbcolor=(0,0,1), thickness=1, legend_label=None,
         legend_color=None, aspect_ratio ='automatic')
def line2d(points, **options):
    r"""
    Create the line through the given list of points.

    INPUT:

    -  ``points`` - either a single point (as a tuple), a list of
       points, a single complex number, or a list of complex numbers.

    Type ``line2d.options`` for a dictionary of the default options for
    lines.  You can change this to change the defaults for all future
    lines.  Use ``line2d.reset()`` to reset to the default options.

    INPUT:

    - ``alpha`` -- How transparent the line is

    - ``thickness`` -- How thick the line is

    - ``rgbcolor`` -- The color as an RGB tuple

    - ``hue`` -- The color given as a hue

    - ``legend_color`` -- The color of the text in the legend

    - ``legend_label`` -- the label for this item in the legend


    Any MATPLOTLIB line option may also be passed in.  E.g.,

    - ``linestyle`` - (default: "-") The style of the line, which is one of
       - ``"-"`` or ``"solid"``
       - ``"--"`` or ``"dashed"``
       - ``"-."`` or ``"dash dot"``
       - ``":"`` or ``"dotted"``
       - ``"None"`` or ``" "`` or ``""`` (nothing)

       The linestyle can also be prefixed with a drawing style (e.g., ``"steps--"``)

       - ``"default"`` (connect the points with straight lines)
       - ``"steps"`` or ``"steps-pre"`` (step function; horizontal
         line is to the left of point)
       - ``"steps-mid"`` (step function; points are in the middle of
         horizontal lines)
       - ``"steps-post"`` (step function; horizontal line is to the
         right of point)

    - ``marker``  - The style of the markers, which is one of
       - ``"None"`` or ``" "`` or ``""`` (nothing) -- default
       - ``","`` (pixel), ``"."`` (point)
       - ``"_"`` (horizontal line), ``"|"`` (vertical line)
       - ``"o"`` (circle), ``"p"`` (pentagon), ``"s"`` (square), ``"x"`` (x), ``"+"`` (plus), ``"*"`` (star)
       - ``"D"`` (diamond), ``"d"`` (thin diamond)
       - ``"H"`` (hexagon), ``"h"`` (alternative hexagon)
       - ``"<"`` (triangle left), ``">"`` (triangle right), ``"^"`` (triangle up), ``"v"`` (triangle down)
       - ``"1"`` (tri down), ``"2"`` (tri up), ``"3"`` (tri left), ``"4"`` (tri right)
       - ``0`` (tick left), ``1`` (tick right), ``2`` (tick up), ``3`` (tick down)
       - ``4`` (caret left), ``5`` (caret right), ``6`` (caret up), ``7`` (caret down)
       - ``"$...$"`` (math TeX string)

    - ``markersize`` -- the size of the marker in points

    - ``markeredgecolor`` -- the color of the marker edge

    - ``markerfacecolor`` -- the color of the marker face

    - ``markeredgewidth`` -- the size of the marker edge in points

    EXAMPLES:

    A line with no points or one point::

        sage: line([])      #returns an empty plot
        Graphics object consisting of 0 graphics primitives
        sage: import numpy; line(numpy.array([]))
        Graphics object consisting of 0 graphics primitives
        sage: line([(1,1)])
        Graphics object consisting of 1 graphics primitive

    A line with numpy arrays::

        sage: line(numpy.array([[1,2], [3,4]]))
        Graphics object consisting of 1 graphics primitive

    A line with a legend::

        sage: line([(0,0),(1,1)], legend_label='line')
        Graphics object consisting of 1 graphics primitive

    Lines with different colors in the legend text::

        sage: p1 = line([(0,0),(1,1)], legend_label='line')
        sage: p2 = line([(1,1),(2,4)], legend_label='squared', legend_color='red')
        sage: p1 + p2
        Graphics object consisting of 2 graphics primitives

    Extra options will get passed on to show(), as long as they are valid::

        sage: line([(0,1), (3,4)], figsize=[10, 2])
        Graphics object consisting of 1 graphics primitive
        sage: line([(0,1), (3,4)]).show(figsize=[10, 2]) # These are equivalent

    We can also use a logarithmic scale if the data will support it::

        sage: line([(1,2),(2,4),(3,4),(4,8),(4.5,32)],scale='loglog',base=2)
        Graphics object consisting of 1 graphics primitive

    Many more examples below!

    A blue conchoid of Nicomedes::

        sage: L = [[1+5*cos(pi/2+pi*i/100), tan(pi/2+pi*i/100)*(1+5*cos(pi/2+pi*i/100))] for i in range(1,100)]
        sage: line(L, rgbcolor=(1/4,1/8,3/4))
        Graphics object consisting of 1 graphics primitive

    A line with 2 complex points::

        sage: i = CC.0
        sage: line([1+i, 2+3*i])
        Graphics object consisting of 1 graphics primitive

    A blue hypotrochoid (3 leaves)::

        sage: n = 4; h = 3; b = 2
        sage: L = [[n*cos(pi*i/100)+h*cos((n/b)*pi*i/100),n*sin(pi*i/100)-h*sin((n/b)*pi*i/100)] for i in range(200)]
        sage: line(L, rgbcolor=(1/4,1/4,3/4))
        Graphics object consisting of 1 graphics primitive

    A blue hypotrochoid (4 leaves)::

        sage: n = 6; h = 5; b = 2
        sage: L = [[n*cos(pi*i/100)+h*cos((n/b)*pi*i/100),n*sin(pi*i/100)-h*sin((n/b)*pi*i/100)] for i in range(200)]
        sage: line(L, rgbcolor=(1/4,1/4,3/4))
        Graphics object consisting of 1 graphics primitive

    A red limacon of Pascal::

        sage: L = [[sin(pi*i/100)+sin(pi*i/50),-(1+cos(pi*i/100)+cos(pi*i/50))] for i in range(-100,101)]
        sage: line(L, rgbcolor=(1,1/4,1/2))
        Graphics object consisting of 1 graphics primitive

    A light green trisectrix of Maclaurin::

        sage: L = [[2*(1-4*cos(-pi/2+pi*i/100)^2),10*tan(-pi/2+pi*i/100)*(1-4*cos(-pi/2+pi*i/100)^2)] for i in range(1,100)]
        sage: line(L, rgbcolor=(1/4,1,1/8))
        Graphics object consisting of 1 graphics primitive

    A green lemniscate of Bernoulli::

        sage: cosines = [cos(-pi/2+pi*i/100) for i in range(201)]
        sage: v = [(1/c, tan(-pi/2+pi*i/100)) for i,c in enumerate(cosines) if c != 0]
        sage: L = [(a/(a^2+b^2), b/(a^2+b^2)) for a,b in v]
        sage: line(L, rgbcolor=(1/4,3/4,1/8))
        Graphics object consisting of 1 graphics primitive

    A red plot of the Jacobi elliptic function `\text{sn}(x,2)`, `-3 < x < 3`::

        sage: L = [(i/100.0, real_part(jacobi('sn', i/100.0, 2.0))) for i in
        ....:      range(-300, 300, 30)]
        sage: line(L, rgbcolor=(3/4, 1/4, 1/8))
        Graphics object consisting of 1 graphics primitive

    A red plot of `J`-Bessel function `J_2(x)`, `0 < x < 10`::

        sage: L = [(i/10.0, bessel_J(2,i/10.0)) for i in range(100)]
        sage: line(L, rgbcolor=(3/4,1/4,5/8))
        Graphics object consisting of 1 graphics primitive


    A purple plot of the Riemann zeta function `\zeta(1/2 + it)`, `0 < t < 30`::

        sage: i = CDF.gen()
        sage: v = [zeta(0.5 + n/10 * i) for n in range(300)]
        sage: L = [(z.real(), z.imag()) for z in v]
        sage: line(L, rgbcolor=(3/4,1/2,5/8))
        Graphics object consisting of 1 graphics primitive

    A purple plot of the Hasse-Weil `L`-function `L(E, 1 + it)`, `-1 < t < 10`::

        sage: E = EllipticCurve('37a')
        sage: vals = E.lseries().values_along_line(1-I, 1+10*I, 100) # critical line
        sage: L = [(z[1].real(), z[1].imag()) for z in vals]
        sage: line(L, rgbcolor=(3/4,1/2,5/8))
        Graphics object consisting of 1 graphics primitive

    A red, blue, and green "cool cat"::

        sage: G = plot(-cos(x), -2, 2, thickness=5, rgbcolor=(0.5,1,0.5))
        sage: P = polygon([[1,2], [5,6], [5,0]], rgbcolor=(1,0,0))
        sage: Q = polygon([(-x,y) for x,y in P[0]], rgbcolor=(0,0,1))
        sage: G + P + Q   # show the plot
        Graphics object consisting of 3 graphics primitives

    TESTS:

    Check that :trac:`13690` is fixed. The legend label should have circles
    as markers.::

        sage: line(enumerate(range(2)), marker='o', legend_label='circle')
        Graphics object consisting of 1 graphics primitive
    """
    from sage.plot.all import Graphics
    from sage.plot.plot import xydata_from_point_list
    points = list(points) # make sure points is a python list
    if not points:
        return Graphics()
    xdata, ydata = xydata_from_point_list(points)
    g = Graphics()
    g._set_extra_kwds(Graphics._extract_kwds_for_show(options))
    g.add_primitive(Line(xdata, ydata, options))
    if options['legend_label']:
        g.legend(True)
        g._legend_colors = [options['legend_color']]
    return g
