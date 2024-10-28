"""
Circles
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
from .primitive import GraphicPrimitive
from sage.misc.decorators import options, rename_keyword
from sage.plot.colors import to_mpl_color
from math import sin, cos, pi


class Circle(GraphicPrimitive):
    """
    Primitive class for the :class:`Circle` graphics type.  See ``circle?`` for information
    about actually plotting circles.

    INPUT:

    - ``x`` -- `x`-coordinate of center of Circle

    - ``y`` -- `y`-coordinate of center of Circle

    - ``r`` -- radius of Circle object

    - ``options`` -- dictionary of valid plot options to pass to constructor

    EXAMPLES:

    Note this should normally be used indirectly via ``circle``::

        sage: from sage.plot.circle import Circle
        sage: C = Circle(2,3,5,{'zorder':2})
        sage: C
        Circle defined by (2.0,3.0) with r=5.0
        sage: C.options()['zorder']
        2
        sage: C.r
        5.0

    TESTS:

    We test creating a circle::

        sage: C = circle((2,3), 5)
    """
    def __init__(self, x, y, r, options):
        """
        Initialize base class Circle.

        EXAMPLES::

            sage: C = circle((2,3), 5, edgecolor='red', alpha=.5, fill=True)
            sage: C[0].x
            2.0
            sage: C[0].r
            5.0
            sage: C[0].options()['edgecolor']
            'red'
            sage: C[0].options()['alpha']
            0.500000000000000
        """
        self.x = float(x)
        self.y = float(y)
        self.r = float(r)
        GraphicPrimitive.__init__(self, options)

    def get_minmax_data(self):
        """
        Return a dictionary with the bounding box data.

        EXAMPLES::

            sage: p = circle((3, 3), 1)
            sage: d = p.get_minmax_data()
            sage: d['xmin']
            2.0
            sage: d['ymin']
            2.0
        """
        from sage.plot.plot import minmax_data
        return minmax_data([self.x - self.r, self.x + self.r],
                           [self.y - self.r, self.y + self.r],
                           dict=True)

    def _allowed_options(self):
        """
        Return the allowed options for the Circle class.

        EXAMPLES::

            sage: p = circle((3, 3), 1)
            sage: p[0]._allowed_options()['alpha']
            'How transparent the figure is.'
            sage: p[0]._allowed_options()['facecolor']
            '2D only: The color of the face as an RGB tuple.'
        """
        return {'alpha': 'How transparent the figure is.',
                'fill': 'Whether or not to fill the circle.',
                'legend_label': 'The label for this item in the legend.',
                'legend_color': 'The color of the legend text.',
                'thickness': 'How thick the border of the circle is.',
                'edgecolor': '2D only: The color of the edge as an RGB tuple.',
                'facecolor': '2D only: The color of the face as an RGB tuple.',
                'rgbcolor': 'The color (edge and face) as an RGB tuple.',
                'hue': 'The color given as a hue.',
                'zorder': '2D only: The layer level in which to draw',
                'linestyle': "2D only: The style of the line, which is one of "
                "'dashed', 'dotted', 'solid', 'dashdot', or '--', ':', '-', '-.', "
                "respectively.",
                'clip': 'Whether or not to clip the circle.'}

    def _repr_(self):
        """
        String representation of Circle primitive.

        EXAMPLES::

            sage: C = circle((2,3), 5)
            sage: c = C[0]; c
            Circle defined by (2.0,3.0) with r=5.0
        """
        return f"Circle defined by ({self.x},{self.y}) with r={self.r}"

    def _render_on_subplot(self, subplot):
        """
        TESTS::

            sage: from math import pi
            sage: C = circle((2,pi), 2, edgecolor='black', facecolor='green', fill=True)
        """
        import matplotlib.patches as patches
        from sage.plot.misc import get_matplotlib_linestyle

        options = self.options()
        p = patches.Circle((float(self.x), float(self.y)), float(self.r), clip_on=options['clip'])
        if not options['clip']:
            self._bbox_extra_artists = [p]
        p.set_linewidth(float(options['thickness']))
        p.set_fill(options['fill'])
        a = float(options['alpha'])
        p.set_alpha(a)
        ec = to_mpl_color(options['edgecolor'])
        fc = to_mpl_color(options['facecolor'])
        if 'rgbcolor' in options:
            ec = fc = to_mpl_color(options['rgbcolor'])
        p.set_edgecolor(ec)
        p.set_facecolor(fc)
        p.set_linestyle(get_matplotlib_linestyle(options['linestyle'],return_type='long'))
        p.set_label(options['legend_label'])
        z = int(options.pop('zorder', 0))
        p.set_zorder(z)
        subplot.add_patch(p)

    def plot3d(self, z=0, **kwds):
        """
        Plots a 2D circle (actually a 50-gon) in 3D,
        with default height zero.

        INPUT:

        - ``z`` -- (optional) 3D height above `xy`-plane

        EXAMPLES::

            sage: circle((0,0), 1).plot3d()
            Graphics3d Object

        This example uses this method implicitly, but does not pass
        the optional parameter z to this method::

            sage: sum(circle((random(),random()), random()).plot3d(z=random())
            ....:     for _ in range(20))
            Graphics3d Object

        .. PLOT::

            P = sum([circle((random(),random()), random()).plot3d(z=random()) for _ in range(20)])
            sphinx_plot(P)

        These examples are explicit, and pass z to this method::

            sage: from math import pi
            sage: C = circle((2,pi), 2, hue=.8, alpha=.3, fill=True)
            sage: c = C[0]
            sage: d = c.plot3d(z=2)
            sage: d.texture.opacity
            0.3

        ::

            sage: C = circle((2,pi), 2, hue=.8, alpha=.3, linestyle='dotted')
            sage: c = C[0]
            sage: d = c.plot3d(z=2)
            sage: d.jmol_repr(d.testing_render_params())[0][-1]
            'color $line_1 translucent 0.7 [204,0,255]'
        """
        options = dict(self.options())
        fill = options['fill']
        for s in ['clip', 'edgecolor', 'facecolor', 'fill', 'linestyle',
                'zorder']:
            if s in options:
                del options[s]

        n = 50
        dt = float(2*pi/n)
        x, y, r = self.x, self.y, self.r
        xdata = [x+r*cos(t*dt) for t in range(n+1)]
        ydata = [y+r*sin(t*dt) for t in range(n+1)]
        if fill:
            from .polygon import Polygon
            return Polygon(xdata, ydata, options).plot3d(z)
        else:
            from .line import Line
            return Line(xdata, ydata, options).plot3d().translate((0,0,z))


@rename_keyword(color='rgbcolor')
@options(alpha=1, fill=False, thickness=1, edgecolor='blue', facecolor='blue', linestyle='solid',
         zorder=5, legend_label=None, legend_color=None, clip=True, aspect_ratio=1.0)
def circle(center, radius, **options):
    """
    Return a circle at a point center = `(x,y)` (or `(x,y,z)` and
    parallel to the `xy`-plane) with radius = `r`.  Type
    ``circle.options`` to see all options.

    OPTIONS:

    - ``alpha`` -- (default: 1)

    - ``fill`` -- (default: ``False``)

    - ``thickness`` -- (default: 1)

    - ``linestyle`` -- (default: ``'solid'``) (2D plotting only) the style of the
      line, which is one of ``'dashed'``, ``'dotted'``, ``'solid'``, ``'dashdot'``,
      or ``'--'``, ``':'``, ``'-'``, ``'-.'``, respectively

    - ``edgecolor`` -- (default: ``'blue'``) 2D plotting only

    - ``facecolor`` -- (default: ``'blue'``) 2D plotting only, useful only
      if ``fill=True``

    - ``rgbcolor`` -- 2D or 3D plotting.  This option overrides
      ``edgecolor`` and ``facecolor`` for 2D plotting

    - ``legend_label`` -- the label for this item in the legend

    - ``legend_color`` -- the color for the legend label

    EXAMPLES:

    The default color is blue, the default linestyle is solid, but this is easy to change::

        sage: c = circle((1,1), 1)
        sage: c
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(circle((1,1), 1))

    ::

        sage: c = circle((1,1), 1, rgbcolor=(1,0,0), linestyle='-.')
        sage: c
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        c = circle((1,1), 1, rgbcolor=(1,0,0), linestyle='-.')
        sphinx_plot(c)

    We can also use this command to plot three-dimensional circles parallel
    to the `xy`-plane::

        sage: c = circle((1,1,3), 1, rgbcolor=(1,0,0))
        sage: c
        Graphics3d Object
        sage: type(c)
        <class 'sage.plot.plot3d.base.TransformGroup'>

    .. PLOT::

        c = circle((1,1,3), 1, rgbcolor=(1,0,0))
        sphinx_plot(c)

    To correct the aspect ratio of certain graphics, it is necessary
    to show with a ``figsize`` of square dimensions::

        sage: c.show(figsize=[5,5],xmin=-1,xmax=3,ymin=-1,ymax=3)

    Here we make a more complicated plot, with many circles of different colors::

        sage: g = Graphics()
        sage: step = 6; ocur = 1/5; paths = 16
        sage: PI = math.pi    # numerical for speed -- fine for graphics
        sage: for r in range(1,paths+1):
        ....:     for x,y in [((r+ocur)*math.cos(n), (r+ocur)*math.sin(n))
        ....:                 for n in srange(0, 2*PI+PI/step, PI/step)]:
        ....:         g += circle((x,y), ocur, rgbcolor=hue(r/paths))
        ....:     rnext = (r+1)^2
        ....:     ocur = (rnext-r)-ocur
        sage: g.show(xmin=-(paths+1)^2, xmax=(paths+1)^2,
        ....:        ymin=-(paths+1)^2, ymax=(paths+1)^2, figsize=[6,6])

    .. PLOT::

        g = Graphics()
        step=6; ocur=1/5; paths=16;
        PI = math.pi    # numerical for speed -- fine for graphics
        for r in range(1,paths+1):
             for x,y in [((r+ocur)*math.cos(n), (r+ocur)*math.sin(n)) for n in srange(0, 2*PI+PI/step, PI/step)]:
                 g += circle((x,y), ocur, rgbcolor=hue(r*1.0/paths))
             rnext = (r+1)**2
             ocur = (rnext-r)-ocur
        g.set_axes_range(-(paths+1)**2,(paths+1)**2,-(paths+1)**2,(paths+1)**2)
        sphinx_plot(g)

    Note that the ``rgbcolor`` option overrides the other coloring options.
    This produces red fill in a blue circle::

        sage: circle((2,3), 1, fill=True, edgecolor='blue', facecolor='red')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(circle((2,3), 1, fill=True, edgecolor='blue', facecolor='red'))

    This produces an all-green filled circle::

        sage: circle((2,3), 1, fill=True, edgecolor='blue', rgbcolor='green')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(circle((2,3), 1, fill=True, edgecolor='blue', rgbcolor='green'))

    The option ``hue`` overrides *all* other options, so be careful with its use.
    This produces a purplish filled circle::

        sage: circle((2,3), 1, fill=True, edgecolor='blue', rgbcolor='green', hue=.8)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        C = circle((2,3), 1, fill=True, edgecolor='blue', rgbcolor='green', hue=.8)
        sphinx_plot(C)

    And circles with legends::

        sage: circle((4,5), 1, rgbcolor='yellow', fill=True,
        ....:        legend_label='the sun').show(xmin=0, ymin=0)

    .. PLOT::

        C = circle((4,5), 1, rgbcolor='yellow', fill=True, legend_label='the sun')
        C.set_axes_range(xmin=0, ymin=0)
        sphinx_plot(C)

    ::

        sage: circle((4,5), 1,
        ....:        legend_label='the sun', legend_color='yellow').show(xmin=0, ymin=0)

    .. PLOT::

        C = circle((4,5), 1, legend_label='the sun', legend_color='yellow')
        C.set_axes_range(xmin=0, ymin=0)
        sphinx_plot(C)

    Extra options will get passed on to show(), as long as they are valid::

        sage: circle((0, 0), 2, figsize=[10,10])  # That circle is huge!
        Graphics object consisting of 1 graphics primitive

    ::

        sage: circle((0, 0), 2).show(figsize=[10,10])  # These are equivalent

    TESTS:

    We cannot currently plot circles in more than three dimensions::

        sage: circle((1,1,1,1), 1, rgbcolor=(1,0,0))
        Traceback (most recent call last):
        ...
        ValueError: the center of a plotted circle should have two or three coordinates

    The default aspect ratio for a circle is 1.0::

        sage: P = circle((1,1), 1)
        sage: P.aspect_ratio()
        1.0

    Verify that :issue:`36153` does not arise::

        sage: C = circle((1,1), 1, legend_label='test')
    """
    from sage.plot.all import Graphics

    # Reset aspect_ratio to 'automatic' in case scale is 'semilog[xy]'.
    # Otherwise matplotlib complains.
    scale = options.get('scale', None)
    if isinstance(scale, (list, tuple)):
        scale = scale[0]
    if scale == 'semilogy' or scale == 'semilogx':
        options['aspect_ratio'] = 'automatic'

    g = Graphics()
    g._set_extra_kwds(Graphics._extract_kwds_for_show(options))
    g.add_primitive(Circle(center[0], center[1], radius, options))
    if options['legend_label']:
        g.legend(True)
        g._legend_colors = [options['legend_color']]
    if len(center) == 2:
        return g
    elif len(center) == 3:
        return g[0].plot3d(z=center[2])
    raise ValueError('the center of a plotted circle should have '
                     'two or three coordinates')
