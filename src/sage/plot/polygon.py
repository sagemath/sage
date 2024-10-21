# sage_setup: distribution = sagemath-plot
"""
Polygons
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
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.plot.primitive import GraphicPrimitive_xydata
from sage.misc.decorators import options, rename_keyword
from sage.plot.colors import to_mpl_color


class Polygon(GraphicPrimitive_xydata):
    """
    Primitive class for the Polygon graphics type.  For information
    on actual plotting, please see :func:`polygon`, :func:`polygon2d`,
    or :func:`~sage.plot.plot3d.shapes2.polygon3d`.

    INPUT:

    - ``xdata`` -- list of `x`-coordinates of points defining Polygon

    - ``ydata`` -- list of `y`-coordinates of points defining Polygon

    - ``options`` -- dictionary of valid plot options to pass to constructor

    EXAMPLES:

    Note this should normally be used indirectly via :func:`polygon`::

        sage: from sage.plot.polygon import Polygon
        sage: P = Polygon([1,2,3],[2,3,2],{'alpha':.5})
        sage: P
        Polygon defined by 3 points
        sage: P.options()['alpha']
        0.500000000000000
        sage: P.ydata
        [2, 3, 2]

    TESTS:

    We test creating polygons::

        sage: polygon([(0,0), (1,1), (0,1)])
        Graphics object consisting of 1 graphics primitive

    ::

        sage: polygon([(0,0,1), (1,1,1), (2,0,1)])
        Graphics3d Object

    ::

        sage: polygon2d([(1, 1), (0, 1), (1, 0)], fill=False, linestyle='dashed')
        Graphics object consisting of 1 graphics primitive
    """
    def __init__(self, xdata, ydata, options):
        """
        Initialize base class Polygon.

        EXAMPLES::

            sage: P = polygon([(0,0), (1,1), (-1,3)], thickness=2)
            sage: P[0].xdata
            [0.0, 1.0, -1.0]
            sage: P[0].options()['thickness']
            2
        """
        self.xdata = xdata
        self.ydata = ydata
        GraphicPrimitive_xydata.__init__(self, options)

    def _repr_(self):
        """
        String representation of Polygon primitive.

        EXAMPLES::

            sage: P = polygon([(0,0), (1,1), (-1,3)])
            sage: p=P[0]; p
            Polygon defined by 3 points
        """
        return "Polygon defined by %s points" % len(self)

    def __getitem__(self, i):
        """
        Return `i`-th vertex of Polygon primitive.

        It is starting count from 0th vertex.

        EXAMPLES::

            sage: P = polygon([(0,0), (1,1), (-1,3)])
            sage: p=P[0]
            sage: p[0]
            (0.0, 0.0)
        """
        return self.xdata[i], self.ydata[i]

    def __setitem__(self, i, point):
        """
        Change `i`-th vertex of Polygon primitive.

        It is starting count from 0th vertex.

        Note that this only changes a vertex, but does not create new vertices.

        EXAMPLES::

            sage: P = polygon([(0,0), (1,2), (0,1), (-1,2)])
            sage: p=P[0]
            sage: [p[i] for i in range(4)]
            [(0.0, 0.0), (1.0, 2.0), (0.0, 1.0), (-1.0, 2.0)]
            sage: p[2]=(0,.5)
            sage: p[2]
            (0.0, 0.5)
        """
        i = int(i)
        self.xdata[i] = float(point[0])
        self.ydata[i] = float(point[1])

    def __len__(self):
        """
        Return number of vertices of Polygon primitive.

        EXAMPLES::

            sage: P = polygon([(0,0), (1,2), (0,1), (-1,2)])
            sage: p=P[0]
            sage: len(p)
            4
        """
        return len(self.xdata)

    def _allowed_options(self):
        """
        Return the allowed options for the Polygon class.

        EXAMPLES::

            sage: P = polygon([(1,1), (1,2), (2,2), (2,1)], alpha=.5)
            sage: P[0]._allowed_options()['alpha']
            'How transparent the figure is.'
        """
        return {'alpha': 'How transparent the figure is.',
                'thickness': 'How thick the border line is.',
                'edgecolor': 'The color for the border of filled polygons.',
                'fill': 'Whether or not to fill the polygon.',
                'legend_label': 'The label for this item in the legend.',
                'legend_color': 'The color of the legend text.',
                'linestyle': 'The style of the enclosing line.',
                'rgbcolor': 'The color as an RGB tuple.',
                'hue': 'The color given as a hue.',
                'zorder': 'The layer level in which to draw'}

    def _plot3d_options(self, options=None):
        """
        Translate 2d plot options into 3d plot options.

        EXAMPLES::

            sage: P = polygon([(1,1), (1,2), (2,2), (2,1)], alpha=.5)
            sage: p=P[0]; p
            Polygon defined by 4 points
            sage: q=p.plot3d()
            sage: q.texture.opacity
            0.5
        """
        if options is None:
            options = dict(self.options())
        for o in ['thickness', 'zorder', 'legend_label', 'fill', 'edgecolor']:
            options.pop(o, None)
        return GraphicPrimitive_xydata._plot3d_options(self, options)

    def plot3d(self, z=0, **kwds):
        """
        Plots a 2D polygon in 3D, with default height zero.

        INPUT:

        - ``z`` -- (optional) 3D height above `xy`-plane, or a list of
          heights corresponding to the list of 2D polygon points

        EXAMPLES:

        A pentagon::

            sage: polygon([(cos(t), sin(t))                                             # needs sage.symbolic
            ....:          for t in srange(0, 2*pi, 2*pi/5)]).plot3d()
            Graphics3d Object

        .. PLOT::

            L = polygon([(cos(t), sin(t)) for t in srange(0, 2*pi, 2*pi/5)]).plot3d()
            sphinx_plot(L)

        Showing behavior of the optional parameter z::

            sage: P = polygon([(0,0), (1,2), (0,1), (-1,2)])
            sage: p = P[0]; p
            Polygon defined by 4 points
            sage: q = p.plot3d()
            sage: q.obj_repr(q.testing_render_params())[2]
            ['v 0 0 0', 'v 1 2 0', 'v 0 1 0', 'v -1 2 0']
            sage: r = p.plot3d(z=3)
            sage: r.obj_repr(r.testing_render_params())[2]
            ['v 0 0 3', 'v 1 2 3', 'v 0 1 3', 'v -1 2 3']
            sage: s = p.plot3d(z=[0,1,2,3])
            sage: s.obj_repr(s.testing_render_params())[2]
            ['v 0 0 0', 'v 1 2 1', 'v 0 1 2', 'v -1 2 3']

        TESTS:

        Heights passed as a list should have same length as
        number of points::

            sage: P = polygon([(0,0), (1,2), (0,1), (-1,2)])
            sage: p = P[0]
            sage: q = p.plot3d(z=[2,-2])
            Traceback (most recent call last):
            ...
            ValueError: Incorrect number of heights given
        """
        from sage.plot.plot3d.index_face_set import IndexFaceSet
        options = self._plot3d_options()
        options.update(kwds)
        zdata = []
        if isinstance(z, list):
            zdata = z
        else:
            zdata = [z] * len(self.xdata)
        if len(zdata) == len(self.xdata):
            return IndexFaceSet([list(zip(self.xdata, self.ydata, zdata))],
                                **options)
        else:
            raise ValueError('Incorrect number of heights given')

    def _render_on_subplot(self, subplot):
        """
        TESTS::

            sage: P = polygon([(0,0), (1,2), (0,1), (-1,2)])
        """
        import matplotlib.patches as patches
        options = self.options()
        p = patches.Polygon([(self.xdata[i], self.ydata[i])
                             for i in range(len(self.xdata))])
        p.set_linewidth(float(options['thickness']))
        if 'linestyle' in options:
            p.set_linestyle(options['linestyle'])
        a = float(options['alpha'])
        z = int(options.pop('zorder', 1))
        p.set_alpha(a)
        f = options.pop('fill')
        p.set_fill(f)
        c = to_mpl_color(options['rgbcolor'])
        if f:
            ec = options['edgecolor']
            if ec is None:
                p.set_color(c)
            else:
                p.set_facecolor(c)
                p.set_edgecolor(to_mpl_color(ec))
        else:
            p.set_color(c)
        p.set_label(options['legend_label'])
        p.set_zorder(z)
        subplot.add_patch(p)


def polygon(points, **options):
    """
    Return either a 2-dimensional or 3-dimensional polygon depending
    on value of points.

    For information regarding additional arguments, see either
    :func:`polygon2d` or :func:`~sage.plot.plot3d.shapes2.polygon3d`.
    Options may be found and set using the dictionaries ``polygon2d.options``
    and ``polygon3d.options``.

    EXAMPLES::

        sage: polygon([(0,0), (1,1), (0,1)])
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(polygon([(0,0), (1,1), (0,1)]))

    ::

        sage: polygon([(0,0,1), (1,1,1), (2,0,1)])
        Graphics3d Object

    Extra options will get passed on to show(), as long as they are valid::

        sage: polygon([(0,0), (1,1), (0,1)], axes=False)
        Graphics object consisting of 1 graphics primitive
        sage: polygon([(0,0), (1,1), (0,1)]).show(axes=False) # These are equivalent
    """
    try:
        return polygon2d(points, **options)
    except ValueError:
        from sage.plot.plot3d.shapes2 import polygon3d
        return polygon3d(points, **options)


@rename_keyword(color='rgbcolor')
@options(alpha=1, rgbcolor=(0, 0, 1), edgecolor=None, thickness=None,
         legend_label=None, legend_color=None,
         aspect_ratio=1.0, fill=True)
def polygon2d(points, **options):
    r"""
    Return a 2-dimensional polygon defined by ``points``.

    Type ``polygon2d.options`` for a dictionary of the default
    options for polygons.  You can change this to change the
    defaults for all future polygons.  Use ``polygon2d.reset()``
    to reset to the default options.

    EXAMPLES:

    We create a purple-ish polygon::

        sage: polygon2d([[1,2], [5,6], [5,0]], rgbcolor=(1,0,1))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(polygon2d([[1,2], [5,6], [5,0]], rgbcolor=(1,0,1)))

    By default, polygons are filled in, but we can make them
    without a fill as well::

        sage: polygon2d([[1,2], [5,6], [5,0]], fill=False)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(polygon2d([[1,2], [5,6], [5,0]], fill=False))

    In either case, the thickness of the border can be controlled::

        sage: polygon2d([[1,2], [5,6], [5,0]], fill=False, thickness=4, color='orange')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        P = polygon2d([[1,2], [5,6], [5,0]], fill=False, thickness=4, color='orange')
        sphinx_plot(P)

    For filled polygons, one can use different colors for the border
    and the interior as follows::

        sage: L = [[0,0]]+[[i/100, 1.1+cos(i/20)] for i in range(100)]+[[1,0]]          # needs sage.symbolic
        sage: polygon2d(L, color='limegreen', edgecolor='black', axes=False)            # needs sage.symbolic
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        L = [[0,0]]+[[i*0.01, 1.1+cos(i*0.05)] for i in range(100)]+[[1,0]]
        P = polygon2d(L, color='limegreen', edgecolor='black', axes=False)
        sphinx_plot(P)

    Some modern art -- a random polygon, with legend::

        sage: v = [(randrange(-5,5), randrange(-5,5)) for _ in range(10)]
        sage: polygon2d(v, legend_label='some form')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        v = [(randrange(-5,5), randrange(-5,5)) for _ in range(10)]
        P = polygon2d(v, legend_label='some form')
        sphinx_plot(P)

    An aperiodic monotile, [Smi2023]_::

        sage: s = sqrt(3)                                                               # needs sage.symbolic
        sage: polygon2d([[0, 0], [0, s], [1, s], [3/2, 3/2*s], [3, s], [3, 0], [4, 0],  # needs sage.symbolic
        ....:            [9/2, -1/2*s], [3, -s], [3/2, -1/2*s], [1, -s], [-1, -s],
        ....:            [-3/2, -1/2*s]], axes=False)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        s = sqrt(3)
        P = polygon2d([[0, 0], [0, s], [1, s], [3/2, 3/2*s], [3, s], [3, 0], [4, 0], [9/2, -1/2*s], [3, -s], \
                       [3/2, -1/2*s], [1, -s], [-1, -s], [-3/2, -1/2*s]], axes=False)
        sphinx_plot(P)

    A purple hexagon::

        sage: L = [[cos(pi*i/3),sin(pi*i/3)] for i in range(6)]                         # needs sage.symbolic
        sage: polygon2d(L, rgbcolor=(1,0,1))                                            # needs sage.symbolic
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        L = [[cos(pi*i/3.0),sin(pi*i/3.0)] for i in range(6)]
        P = polygon2d(L, rgbcolor=(1,0,1))
        sphinx_plot(P)

    A green deltoid::

        sage: L = [[-1+cos(pi*i/100)*(1+cos(pi*i/100)),                                 # needs sage.symbolic
        ....:       2*sin(pi*i/100)*(1-cos(pi*i/100))] for i in range(200)]
        sage: polygon2d(L, rgbcolor=(1/8,3/4,1/2))                                      # needs sage.symbolic
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        L = [[-1+cos(pi*i*0.01)*(1+cos(pi*i*0.01)),2*sin(pi*i*0.01)*(1-cos(pi*i*0.01))] for i in range(200)]
        P = polygon2d(L, rgbcolor=(0.125,0.75,0.5))
        sphinx_plot(P)

    A blue hypotrochoid::

        sage: L = [[6*cos(pi*i/100)+5*cos((6/2)*pi*i/100),                              # needs sage.symbolic
        ....:       6*sin(pi*i/100)-5*sin((6/2)*pi*i/100)] for i in range(200)]
        sage: polygon2d(L, rgbcolor=(1/8,1/4,1/2))                                      # needs sage.symbolic
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        L = [[6*cos(pi*i*0.01)+5*cos(3*pi*i*0.01),6*sin(pi*i*0.01)-5*sin(3*pi*i*0.01)] for i in range(200)]
        P = polygon2d(L, rgbcolor=(0.125,0.25,0.5))
        sphinx_plot(P)

    Another one::

        sage: n = 4; h = 5; b = 2
        sage: L = [[n*cos(pi*i/100)+h*cos((n/b)*pi*i/100),                              # needs sage.symbolic
        ....:       n*sin(pi*i/100)-h*sin((n/b)*pi*i/100)] for i in range(200)]
        sage: polygon2d(L, rgbcolor=(1/8,1/4,3/4))                                      # needs sage.symbolic
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        n = 4.0; h = 5.0; b = 2.0
        L = [[n*cos(pi*i*0.01)+h*cos((n/b)*pi*i*0.01),n*sin(pi*i*0.01)-h*sin((n/b)*pi*i*0.01)] for i in range(200)]
        P = polygon2d(L, rgbcolor=(0.125,0.25,0.75))
        sphinx_plot(P)

    A purple epicycloid::

        sage: m = 9; b = 1
        sage: L = [[m*cos(pi*i/100)+b*cos((m/b)*pi*i/100),                              # needs sage.symbolic
        ....:       m*sin(pi*i/100)-b*sin((m/b)*pi*i/100)] for i in range(200)]
        sage: polygon2d(L, rgbcolor=(7/8,1/4,3/4))                                      # needs sage.symbolic
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        m = 9.0; b = 1
        L = [[m*cos(pi*i*0.01)+b*cos((m/b)*pi*i*0.01),m*sin(pi*i*0.01)-b*sin((m/b)*pi*i*0.01)] for i in range(200)]
        P = polygon2d(L, rgbcolor=(0.875,0.25,0.75))
        sphinx_plot(P)

    A brown astroid::

        sage: L = [[cos(pi*i/100)^3, sin(pi*i/100)^3] for i in range(200)]              # needs sage.symbolic
        sage: polygon2d(L, rgbcolor=(3/4,1/4,1/4))                                      # needs sage.symbolic
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        L = [[cos(pi*i*0.01)**3,sin(pi*i*0.01)**3] for i in range(200)]
        P = polygon2d(L, rgbcolor=(0.75,0.25,0.25))
        sphinx_plot(P)

    And, my favorite, a greenish blob::

        sage: L = [[cos(pi*i/100)*(1+cos(pi*i/50)),                                     # needs sage.symbolic
        ....:       sin(pi*i/100)*(1+sin(pi*i/50))] for i in range(200)]
        sage: polygon2d(L, rgbcolor=(1/8,3/4,1/2))                                      # needs sage.symbolic
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        L = [[cos(pi*i*0.01)*(1+cos(pi*i*0.02)), sin(pi*i*0.01)*(1+sin(pi*i*0.02))] for i in range(200)]
        P = polygon2d(L, rgbcolor=(0.125,0.75,0.5))
        sphinx_plot(P)

    This one is for my wife::

        sage: L = [[sin(pi*i/100)+sin(pi*i/50),                                         # needs sage.symbolic
        ....:       -(1+cos(pi*i/100)+cos(pi*i/50))] for i in range(-100,100)]
        sage: polygon2d(L, rgbcolor=(1,1/4,1/2))                                        # needs sage.symbolic
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        L = [[sin(pi*i*0.01)+sin(pi*i*0.02),-(1+cos(pi*i*0.01)+cos(pi*i*0.02))] for i in range(-100,100)]
        P = polygon2d(L, rgbcolor=(1,0.25,0.5))
        sphinx_plot(P)

    One can do the same one with a colored legend label::

        sage: polygon2d(L, color='red', legend_label='For you!', legend_color='red')    # needs sage.symbolic
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        L = [[sin(pi*i*0.01)+sin(pi*i*0.02),-(1+cos(pi*i*0.01)+cos(pi*i*0.02))] for i in range(-100,100)]
        P = polygon2d(L, color='red', legend_label='For you!', legend_color='red')
        sphinx_plot(P)

    Polygons have a default aspect ratio of 1.0::

        sage: polygon2d([[1,2], [5,6], [5,0]]).aspect_ratio()
        1.0

    TESTS:

    Verify that :issue:`36153` does not arise::

        sage: P = polygon2d([[1,2], [5,6], [5,0]], legend_label='test')

    AUTHORS:

    - David Joyner (2006-04-14): the long list of examples above.
    """
    from sage.plot.plot import xydata_from_point_list
    from sage.plot.all import Graphics
    if options["thickness"] is None:   # If the user did not specify thickness
        if options["fill"] and options["edgecolor"] is None:
            # If the user chose fill
            options["thickness"] = 0
        else:
            options["thickness"] = 1
    xdata, ydata = xydata_from_point_list(points)
    g = Graphics()

    # Reset aspect_ratio to 'automatic' in case scale is 'semilog[xy]'.
    # Otherwise matplotlib complains.
    scale = options.get('scale', None)
    if isinstance(scale, (list, tuple)):
        scale = scale[0]
    if scale == 'semilogy' or scale == 'semilogx':
        options['aspect_ratio'] = 'automatic'

    g._set_extra_kwds(Graphics._extract_kwds_for_show(options))
    g.add_primitive(Polygon(xdata, ydata, options))
    if options['legend_label']:
        g.legend(True)
        g._legend_colors = [options['legend_color']]
    return g
