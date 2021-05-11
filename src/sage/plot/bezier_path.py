r"""
Bezier Paths
"""
#*****************************************************************************
#       Copyright (C) 2006 Alex Clemesha <clemesha@gmail.com>,
#                          William Stein <wstein@gmail.com>,
#                     2008 Mike Hansen <mhansen@gmail.com>,
#                     2009 Emily Kirkman
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


class BezierPath(GraphicPrimitive_xydata):
    """
    Path of Bezier Curves graphics primitive.

    The input to this constructor is a list of curves, each a list of points,
    along which to create the curves, along with a dict of any options passed.

    EXAMPLES::

        sage: from sage.plot.bezier_path import BezierPath
        sage: BezierPath([[(0,0), (.5,.5),(1,0)],[(.5,1),(0,0)]], {'linestyle':'dashed'})
        Bezier path from (0.0, 0.0) to (0.0, 0.0)

    We use :func:`bezier_path` to actually plot Bezier curves::

        sage: bezier_path([[(0,0),(.5,.5),(1,0)],[(.5,1),(0,0)]], linestyle="dashed")
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

         P = bezier_path([[(0,0),(.5,.5),(1,0)],[(.5,1),(0,0)]], linestyle="dashed")
         sphinx_plot(P)

    """
    def __init__(self, path, options):
        """
        Returns a graphics primitive of a path of Bezier curves.

        EXAMPLES::

            sage: from sage.plot.bezier_path import BezierPath
            sage: BezierPath([[(0,0),(.5,.5),(1,0)],[(.5,1),(0,0)]], {'linestyle':'dashed'})
            Bezier path from (0.0, 0.0) to (0.0, 0.0)

            sage: BezierPath([[(0,0), (1,2), (3,6), (2,-1), (3,3)]], {})
            Traceback (most recent call last):
            ...
            ValueError: invalid input for BezierPath

        TESTS:

        Check :trac:`31646`::

            sage: from sage.plot.bezier_path import BezierPath
            sage: p2d = [[(3,0),(4,1),(2,1),(3,0)], [(2,2),(3,1),(2,1)]]
            sage: P = BezierPath(p2d, {})
            sage: P.path
            [array([[3., 0.], [4., 1.], [2., 1.], [3., 0.]]),
             array([[2., 2.], [3., 1.], [2., 1.]])]
        """
        import numpy as np

        self.path = [np.array(l, float) for l in path]

        # In oder to feed later to matplotlib.path.Path we convert in the following form
        # - vertices: an Nx2 float array of vertices
        # - codes: an N-length uint8 array of vertex types, or None
        #   where each code could be MOVETO (=1), LINETO (=2), CURVE3 (=3), CURVE4 (=4)
        self.vertices = np.concatenate(self.path)
        N, _ = self.vertices.shape
        codes = np.zeros((N,), np.uint8)
        k = 0
        for i, curve in enumerate(self.path):
            code = len(curve) + (i > 0)
            if code < 2 or code > 4:
                raise ValueError('invalid input for BezierPath')
            codes[k:k+len(curve)] = code
            k += len(curve)
        codes[0] = 1 # MOVETO
        self.codes = codes
        GraphicPrimitive_xydata.__init__(self, options)

    def _allowed_options(self):
        """
        Returns a dict of allowed options for ``bezier_path``.

        EXAMPLES::

            sage: from sage.plot.bezier_path import BezierPath
            sage: list(sorted(BezierPath([[[-1,2], [14,2.3], [17,4]]], {})._allowed_options().items()))
            [('alpha', 'How transparent the line is.'),
            ('fill', 'Whether or not to fill the polygon.'),
            ('linestyle',
            "The style of the line, which is one of 'dashed', 'dotted', 'solid',
            'dashdot', or '--', ':', '-', '-.', respectively."),
            ('rgbcolor', 'The color as an RGB tuple.'),
            ('thickness', 'How thick the border of the polygon is.'),
            ('zorder', 'The layer level in which to draw')]

        """
        return {'alpha': 'How transparent the line is.',
                'fill': 'Whether or not to fill the polygon.',
                'thickness': 'How thick the border of the polygon is.',
                'rgbcolor': 'The color as an RGB tuple.',
                'zorder': 'The layer level in which to draw',
                'linestyle': "The style of the line, which is one of 'dashed',"
                " 'dotted', 'solid', 'dashdot', or '--', ':', '-', '-.',"
                " respectively."}

    def _plot3d_options(self, options=None):
        """
        Updates ``BezierPath`` options to those allowed by 3D implementation.

        EXAMPLES::

            sage: from sage.plot.bezier_path import BezierPath
            sage: B = BezierPath([[(0,0),(.5,.5),(1,0)],[(.5,1),(0,0)]], {'linestyle':'dashed'})
            sage: B._plot3d_options()
            Traceback (most recent call last):
            ...
            NotImplementedError: Invalid 3d line style: 'dashed'
            sage: B = BezierPath([[(0,0),(.5,.5),(1,0)],[(.5,1),(0,0)]], {'fill':False, 'thickness':2})
            sage: B._plot3d_options()
            {'thickness': 2}
        """
        if options is None:
            options = dict(self.options())
        options_3d = {}
        if 'thickness' in options:
            options_3d['thickness'] = options['thickness']
            del options['thickness']
        if 'fill' in options:
            if options['fill']:
                raise NotImplementedError("Invalid 3d fill style. Must set fill to False.")
            del options['fill']
        if 'linestyle' in options:
            if options['linestyle'] not in ('solid', '-'):
                raise NotImplementedError("Invalid 3d line style: '%s'" %
                                          (options['linestyle']))
            del options['linestyle']
        options_3d.update(GraphicPrimitive_xydata._plot3d_options(self, options))
        return options_3d

    def plot3d(self, z=0, **kwds):
        """
        Returns a 3D plot (Jmol) of the Bezier path.  Since a ``BezierPath``
        primitive contains only `x,y` coordinates, the path will be drawn in
        some plane (default is `z=0`).  To create a Bezier path with nonzero
        (and nonidentical) `z` coordinates in the path and control points, use
        the function :func:`~sage.plot.plot3d.shapes2.bezier3d` instead of
        :func:`bezier_path`.

        EXAMPLES::

            sage: b = bezier_path([[(0,0),(0,1),(1,0)]])
            sage: A = b.plot3d()
            sage: B = b.plot3d(z=2)
            sage: A + B
            Graphics3d Object

        .. PLOT::

            b = bezier_path([[(0,0),(0,1),(1,0)]])
            A = b.plot3d()
            B = b.plot3d(z=2)
            sphinx_plot(A + B)

        ::

            sage: bezier3d([[(0,0,0),(1,0,0),(0,1,0),(0,1,1)]])
            Graphics3d Object

        .. PLOT::

            sphinx_plot(bezier3d([[(0,0,0),(1,0,0),(0,1,0),(0,1,1)]]))

        """
        from sage.plot.plot3d.shapes2 import bezier3d
        options = self._plot3d_options()
        options.update(kwds)
        return bezier3d([[(x,y,0) for x,y in self.path[i]] for i in range(len(self.path))], **options)

    def _repr_(self):
        """
        Return text representation of this Bezier path graphics primitive.

        EXAMPLES::

            sage: from sage.plot.bezier_path import BezierPath
            sage: B = BezierPath([[(0,0),(.5,.5),(1,0)],[(.5,1),(0,0)]], {'linestyle':'dashed'})
            sage: B._repr_()
            'Bezier path from (0.0, 0.0) to (0.0, 0.0)'
        """
        x0, y0 = self.vertices[0]
        x1, y1 = self.vertices[-1]
        return "Bezier path from (%s, %s) to (%s, %s)" % (x0, y0, x1, y1)

    def _render_on_subplot(self, subplot):
        """
        Render this Bezier path in a subplot.  This is the key function that
        defines how this Bezier path graphics primitive is rendered in matplotlib's
        library.

        TESTS::

            sage: bezier_path([[(0,1),(.5,0),(1,1)]])
            Graphics object consisting of 1 graphics primitive

        ::

            sage: bezier_path([[(0,1),(.5,0),(1,1),(-3,5)]])
            Graphics object consisting of 1 graphics primitive
        """
        from matplotlib.patches import PathPatch
        from matplotlib.path import Path
        from sage.plot.misc import get_matplotlib_linestyle

        options = dict(self.options())

        del options['alpha']
        del options['thickness']
        del options['rgbcolor']
        del options['zorder']
        del options['fill']
        del options['linestyle']

        bpath = Path(self.vertices, self.codes)
        bpatch = PathPatch(bpath, **options)
        options = self.options()
        bpatch.set_linewidth(float(options['thickness']))
        bpatch.set_fill(options['fill'])
        bpatch.set_zorder(options['zorder'])
        a = float(options['alpha'])
        bpatch.set_alpha(a)
        c = to_mpl_color(options['rgbcolor'])
        bpatch.set_edgecolor(c)
        bpatch.set_facecolor(c)
        bpatch.set_linestyle(get_matplotlib_linestyle(options['linestyle'], return_type='long'))
        subplot.add_patch(bpatch)

    def get_minmax_data(self):
        """
        Returns a dictionary with the bounding box data.

        EXAMPLES::

            sage: b = bezier_path([[(0,0),(.5,.5),(1,0)],[(.5,1),(0,0)]])
            sage: d = b.get_minmax_data()
            sage: d['xmin']
            0.0
            sage: d['xmax']
            1.0
        """
        return {'xmin': self.vertices[:,0].min(),
                'xmax': self.vertices[:,0].max(),
                'ymin': self.vertices[:,1].min(),
                'ymax': self.vertices[:,1].max()}


@rename_keyword(color='rgbcolor')
@options(alpha=1, fill=False, thickness=1, rgbcolor=(0,0,0), zorder=2, linestyle='solid')
def bezier_path(path, **options):
    """
    Returns a Graphics object of a Bezier path corresponding to the
    path parameter.  The path is a list of curves, and each curve is
    a list of points.  Each point is a tuple ``(x,y)``.

    The first curve contains the endpoints as the first and last point
    in the list.  All other curves assume a starting point given by the
    last entry in the preceding list, and take the last point in the list
    as their opposite endpoint.  A curve can have 0, 1 or 2 control points
    listed between the endpoints.  In the input example for path below,
    the first and second curves have 2 control points, the third has one,
    and the fourth has no control points:

    path = [[p1, c1, c2, p2], [c3, c4, p3], [c5, p4], [p5], ...]

    In the case of no control points, a straight line will be drawn
    between the two endpoints.  If one control point is supplied, then
    the curve at each of the endpoints will be tangent to the line from
    that endpoint to the control point.  Similarly, in the case of two
    control points, at each endpoint the curve will be tangent to the line
    connecting that endpoint with the control point immediately after or
    immediately preceding it in the list.

    .. PLOT::

        p1 = (0,0)
        c1 = (1,1)
        c2 = (1.5,0.5)
        p2 = (4,-1)
        c3 = (3.5,0)
        c4 = (2,1)
        p3 = (0,2)
        c5 = (0.5,3)
        p4 = (1.5,2)
        p5 = (0,4)
        path = [[p1, c1, c2, p2], [c3, c4, p3], [c5, p4], [p5]]
        P = bezier_path(path)
        P += line([p1,c1], color="red", linestyle="dashed")
        P += line([p2,c2], color="red", linestyle="dashed")
        P += line([p2,c3], color="red", linestyle="dashed")
        P += line([p3,c4], color="red", linestyle="dashed")
        P += line([p3,c5], color="red", linestyle="dashed")
        P += text("c1", c1, horizontal_alignment='left')
        P += text("c2", c2, horizontal_alignment='right')
        P += text("c3", c3, horizontal_alignment='left', vertical_alignment='bottom')
        P += text("c4", c4, horizontal_alignment='left')
        P += text("c5", c5, horizontal_alignment='left')
        P += text("p1", p1, horizontal_alignment='left', vertical_alignment='top')
        P += text("p2", p2, horizontal_alignment='left')
        P += text("p3", p3, horizontal_alignment='right', vertical_alignment='top')
        P += text("p4", p4, horizontal_alignment='left')
        P += text("p5", p5, horizontal_alignment='left', vertical_alignment='bottom')
        P += point([c1, c2, c3, c4, c5])
        sphinx_plot(P)

    So in our example above, the curve between p1 and p2 is tangent to the
    line through p1 and c1 at p1, and tangent to the line through p2 and c2
    at p2.  Similarly, the curve between p2 and p3 is tangent to line(p2,c3)
    at p2 and tangent to line(p3,c4) at p3.  Curve(p3,p4) is tangent to
    line(p3,c5) at p3 and tangent to line(p4,c5) at p4.  Curve(p4,p5) is a
    straight line.

    INPUT:

    - ``path`` -- a list of lists of tuples (see above)
    - ``alpha`` -- default: 1
    - ``fill`` -- default: False
    - ``thickness`` -- default: 1
    - ``linestyle`` -- default: ``'solid'``, The style of the line, which is one
       of ``'dashed'``, ``'dotted'``, ``'solid'``, ``'dashdot'``, or ``'--'``,
       ``':'``, ``'-'``, ``'-.'``, respectively.
    - ``rgbcolor`` -- default: (0,0,0)
    - ``zorder`` -- the layer in which to draw

    EXAMPLES::

        sage: path = [[(0,0),(.5,.1),(.75,3),(1,0)],[(.5,1),(.5,0)],[(.2,.5)]]
        sage: b = bezier_path(path, linestyle='dashed', rgbcolor='green')
        sage: b
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        path = [[(0,0),(.5,.1),(.75,3),(1,0)],[(.5,1),(.5,0)],[(.2,.5)]]
        b = bezier_path(path, linestyle='dashed', rgbcolor='green')
        sphinx_plot(b)

    To construct a simple curve, create a list containing a single list::

        sage: path = [[(0,0),(.5,1),(1,0)]]
        sage: curve = bezier_path(path, linestyle='dashed', rgbcolor='green')
        sage: curve
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        path = [[(0,0),(.5,1),(1,0)]]
        curve = bezier_path(path, linestyle='dashed', rgbcolor='green')
        sphinx_plot(curve)

    Extra options will get passed on to :meth:`~Graphics.show`, as long as they are valid::

        sage: bezier_path([[(0,1),(.5,0),(1,1)]], fontsize=50)
        Graphics object consisting of 1 graphics primitive
        sage: bezier_path([[(0,1),(.5,0),(1,1)]]).show(fontsize=50) # These are equivalent

    .. PLOT::

        sphinx_plot(bezier_path([[(0,1),(.5,0),(1,1)]], fontsize=50))

    TESTS:

    We shouldn't modify our argument, :trac:`13822`::

        sage: bp = [[(1,1),(2,3),(3,3)], [(4,4),(5,5)]]
        sage: foo = bezier_path(bp)
        sage: bp
        [[(1, 1), (2, 3), (3, 3)], [(4, 4), (5, 5)]]

    """
    from sage.plot.all import Graphics
    g = Graphics()
    g._set_extra_kwds(g._extract_kwds_for_show(options))
    g.add_primitive(BezierPath(path, options))
    return g
