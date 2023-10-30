"""
Disks
"""
# *****************************************************************************
#       Copyright (C) 2006 Alex Clemesha <clemesha@gmail.com>,
#                          William Stein <wstein@gmail.com>,
#                     2008 Mike Hansen <mhansen@gmail.com>,
#                     2021 Javier Honrubia <jhonrubia6@alumno.uned.es>
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
# *****************************************************************************

from sage.plot.primitive import GraphicPrimitive
from sage.misc.decorators import options, rename_keyword
from sage.plot.colors import to_mpl_color
from math import sin, cos, pi


class Disk(GraphicPrimitive):
    """
    Primitive class for the ``Disk`` graphics type.  See ``disk?`` for
    information about actually plotting a disk (the Sage term for a sector
    or wedge of a circle).

    INPUT:

    - ``point`` - coordinates of center of disk

    - ``r`` - radius of disk

    - ``angle`` - beginning and ending angles of disk (i.e.
      angle extent of sector/wedge)

    - ``options`` - dict of valid plot options to pass to constructor

    EXAMPLES:

    Note this should normally be used indirectly via ``disk``::

        sage: from math import pi
        sage: from sage.plot.disk import Disk
        sage: D = Disk((1,2), 2, (pi/2,pi), {'zorder':3})
        sage: D
        Disk defined by (1.0,2.0) with r=2.0
         spanning (1.5707963267..., 3.1415926535...) radians
        sage: D.options()['zorder']
        3
        sage: D.x
        1.0

    TESTS:

    We test creating a disk::

        sage: disk((2,3), 2, (0,pi/2))
        Graphics object consisting of 1 graphics primitive
    """
    def __init__(self, point, r, angle, options):
        """
        Initializes base class ``Disk``.

        EXAMPLES::

            sage: from math import pi
            sage: D = disk((2,3), 1, (pi/2, pi), fill=False, color='red', thickness=1, alpha=.5)
            sage: D[0].x
            2.0
            sage: D[0].r
            1.0
            sage: D[0].rad1
            1.5707963267948966
            sage: D[0].options()['rgbcolor']
            'red'
            sage: D[0].options()['alpha']
            0.500000000000000
            sage: print(loads(dumps(D)))
            Graphics object consisting of 1 graphics primitive
        """
        self.x = float(point[0])
        self.y = float(point[1])
        self.r = float(r)
        self.rad1 = float(angle[0])
        self.rad2 = float(angle[1])
        GraphicPrimitive.__init__(self, options)

    def get_minmax_data(self):
        """
        Returns a dictionary with the bounding box data.

        EXAMPLES::

            sage: from math import pi
            sage: D = disk((5,4), 1, (pi/2, pi))
            sage: d = D.get_minmax_data()
            sage: d['xmin']
            4.0
            sage: d['ymin']
            3.0
            sage: d['xmax']
            6.0
            sage: d['ymax']
            5.0

        """
        from sage.plot.plot import minmax_data
        return minmax_data([self.x - self.r, self.x + self.r],
                           [self.y - self.r, self.y + self.r],
                           dict=True)

    def _allowed_options(self):
        """
        Return the allowed options for the ``Disk`` class.

        EXAMPLES::

            sage: from math import pi
            sage: p = disk((3, 3), 1, (0, pi/2))
            sage: p[0]._allowed_options()['alpha']
            'How transparent the figure is.'
            sage: p[0]._allowed_options()['zorder']
            'The layer level in which to draw'
        """
        return {'alpha': 'How transparent the figure is.',
                'fill': 'Whether or not to fill the disk.',
                'legend_label': 'The label for this item in the legend.',
                'legend_color': 'The color of the legend text.',
                'thickness': 'How thick the border of the disk is.',
                'rgbcolor': 'The color as an RGB tuple.',
                'hue': 'The color given as a hue.',
                'zorder': 'The layer level in which to draw'}

    def _repr_(self):
        """
        String representation of ``Disk`` primitive.

        EXAMPLES::

            sage: from math import pi
            sage: P = disk((3, 3), 1, (0, pi/2))
            sage: p = P[0]; p
            Disk defined by (3.0,3.0) with r=1.0 spanning (0.0, 1.5707963267...) radians
        """
        return "Disk defined by (%s,%s) with r=%s spanning (%s, %s) radians" % (self.x, self.y, self.r, self.rad1, self.rad2)

    def _render_on_subplot(self, subplot):
        """
        TESTS::

            sage: from math import pi
            sage: D = disk((2,-1), 2, (0, pi), color='black', thickness=3, fill=False); D
            Graphics object consisting of 1 graphics primitive

        Save alpha information in pdf (see :trac:`13732`)::

            sage: f = tmp_filename(ext='.pdf')
            sage: p = disk((0,0), 5, (0, pi/4), alpha=0.5)
            sage: p.save(f)

        """
        import matplotlib.patches as patches
        options = self.options()
        deg1 = self.rad1*(180./pi)  # convert radians to degrees
        deg2 = self.rad2*(180./pi)
        z = int(options.pop('zorder', 0))
        p = patches.Wedge((float(self.x), float(self.y)), float(self.r), float(deg1),
                          float(deg2), zorder=z)
        a = float(options['alpha'])
        p.set_alpha(a)
        p.set_linewidth(float(options['thickness']))
        p.set_fill(options['fill'])
        c = to_mpl_color(options['rgbcolor'])
        p.set_edgecolor(c)
        p.set_facecolor(c)
        p.set_label(options['legend_label'])
        subplot.add_patch(p)

    def plot3d(self, z=0, **kwds):
        """
        Plots a 2D disk (actually a 52-gon) in 3D,
        with default height zero.

        INPUT:


        -  ``z`` - optional 3D height above `xy`-plane.

        AUTHORS:

        - Karl-Dieter Crisman (05-09)

        EXAMPLES::

            sage: from math import pi
            sage: disk((0,0), 1, (0, pi/2)).plot3d()
            Graphics3d Object
            sage: disk((0,0), 1, (0, pi/2)).plot3d(z=2)
            Graphics3d Object
            sage: disk((0,0), 1, (pi/2, 0), fill=False).plot3d(3)
            Graphics3d Object


        These examples show that the appropriate options are passed::

            sage: D = disk((2,3), 1, (pi/4,pi/3), hue=.8, alpha=.3, fill=True)
            sage: d = D[0]
            sage: d.plot3d(z=2).texture.opacity
            0.3

        ::

            sage: D = disk((2,3), 1, (pi/4,pi/3), hue=.8, alpha=.3, fill=False)
            sage: d = D[0]
            sage: dd = d.plot3d(z=2)
            sage: dd.jmol_repr(dd.testing_render_params())[0][-1]
            'color $line_4 translucent 0.7 [204,0,255]'
        """
        options = dict(self.options())
        fill = options['fill']
        del options['fill']
        if 'zorder' in options:
            del options['zorder']
        n = 50
        x, y, r, rad1, rad2 = self.x, self.y, self.r, self.rad1, self.rad2
        dt = float((rad2-rad1)/n)
        xdata = [x]
        ydata = [y]
        xdata.extend([x+r*cos(t*dt+rad1) for t in range(n+1)])
        ydata.extend([y+r*sin(t*dt+rad1) for t in range(n+1)])
        xdata.append(x)
        ydata.append(y)
        if fill:
            from .polygon import Polygon
            return Polygon(xdata, ydata, options).plot3d(z)
        else:
            from .line import Line
            return Line(xdata, ydata, options).plot3d().translate((0, 0, z))


@rename_keyword(color='rgbcolor')
@options(alpha=1, fill=True, rgbcolor=(0, 0, 1), thickness=0, legend_label=None,
         aspect_ratio=1.0)
def disk(point, radius, angle, **options):
    r"""
    A disk (that is, a sector or wedge of a circle) with center
    at a point = `(x,y)` (or `(x,y,z)` and parallel to the
    `xy`-plane) with radius = `r` spanning (in radians)
    angle=`(rad1, rad2)`.

    Type ``disk.options`` to see all options.

    EXAMPLES:

    Make some dangerous disks::

        sage: from math import pi
        sage: bl = disk((0.0,0.0), 1, (pi, 3*pi/2), color='yellow')
        sage: tr = disk((0.0,0.0), 1, (0, pi/2), color='yellow')
        sage: tl = disk((0.0,0.0), 1, (pi/2, pi), color='black')
        sage: br = disk((0.0,0.0), 1, (3*pi/2, 2*pi), color='black')
        sage: P  = tl + tr + bl + br
        sage: P.show(xmin=-2, xmax=2, ymin=-2, ymax=2)

    .. PLOT::

        from sage.plot.disk import Disk
        bl = disk((0.0,0.0), 1, (pi, 3*pi/2), color='yellow')
        tr = disk((0.0,0.0), 1, (0, pi/2), color='yellow')
        tl = disk((0.0,0.0), 1, (pi/2, pi), color='black')
        br = disk((0.0,0.0), 1, (3*pi/2, 2*pi), color='black')
        P  = tl+tr+bl+br
        sphinx_plot(P)

    The default aspect ratio is 1.0::

        sage: disk((0.0,0.0), 1, (pi, 3*pi/2)).aspect_ratio()
        1.0

    Another example of a disk::

        sage: bl = disk((0.0,0.0), 1, (pi, 3*pi/2), rgbcolor=(1,1,0))
        sage: bl.show(figsize=[5,5])

    .. PLOT::

        from sage.plot.disk import Disk
        bl = disk((0.0,0.0), 1, (pi, 3*pi/2), rgbcolor=(1,1,0))
        sphinx_plot(bl)

    Note that since ``thickness`` defaults to zero, it is best to change
    that option when using ``fill=False``::

        sage: disk((2,3), 1, (pi/4,pi/3), hue=.8, alpha=.3, fill=False, thickness=2)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        from sage.plot.disk import Disk
        D = disk((2,3), 1, (pi/4,pi/3), hue=.8, alpha=.3, fill=False, thickness=2)
        sphinx_plot(D)

    The previous two examples also illustrate using ``hue`` and ``rgbcolor``
    as ways of specifying the color of the graphic.

    We can also use this command to plot three-dimensional disks parallel
    to the `xy`-plane::

        sage: d = disk((1,1,3), 1, (pi,3*pi/2), rgbcolor=(1,0,0))
        sage: d
        Graphics3d Object
        sage: type(d)
        <... 'sage.plot.plot3d.index_face_set.IndexFaceSet'>

    .. PLOT::

        from sage.plot.disk import Disk
        d = disk((1,1,3), 1, (pi,3*pi/2), rgbcolor=(1,0,0))
        sphinx_plot(d)

    Extra options will get passed on to ``show()``, as long as they are valid::

        sage: disk((0, 0), 5, (0, pi/2), rgbcolor=(1, 0, 1),
        ....:      xmin=0, xmax=5, ymin=0, ymax=5, figsize=(2,2))
        Graphics object consisting of 1 graphics primitive
        sage: disk((0, 0), 5, (0, pi/2), rgbcolor=(1, 0, 1)).show(  # These are equivalent
        ....:     xmin=0, xmax=5, ymin=0, ymax=5, figsize=(2,2))

    TESTS:

    Testing that legend labels work right::

        sage: disk((2,4), 3, (pi/8, pi/4), hue=1, legend_label='disk', legend_color='blue')
        Graphics object consisting of 1 graphics primitive

    We cannot currently plot disks in more than three dimensions::

        sage: d = disk((1,1,1,1), 1, (0,pi))
        Traceback (most recent call last):
        ...
        ValueError: the center point of a plotted disk should have two or three coordinates
    """
    from sage.plot.all import Graphics
    g = Graphics()

    # Reset aspect_ratio to 'automatic' in case scale is 'semilog[xy]'.
    # Otherwise matplotlib complains.
    scale = options.get('scale', None)
    if isinstance(scale, (list, tuple)):
        scale = scale[0]
    if scale == 'semilogy' or scale == 'semilogx':
        options['aspect_ratio'] = 'automatic'

    g._set_extra_kwds(Graphics._extract_kwds_for_show(options))
    g.add_primitive(Disk(point, radius, angle, options))
    if options['legend_label']:
        g.legend(True)
        g._legend_colors = [options['legend_color']]
    if len(point) == 2:
        return g
    elif len(point) == 3:
        return g[0].plot3d(z=point[2])
    raise ValueError('the center point of a plotted disk should have '
                     'two or three coordinates')
