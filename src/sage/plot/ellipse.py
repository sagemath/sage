"""
Ellipses
"""
# ****************************************************************************
#       Copyright (C) 2010 Vincent Delecroix <20100.delecroix@gmail.com>
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

from .primitive import GraphicPrimitive
from sage.misc.decorators import options, rename_keyword
from sage.plot.colors import to_mpl_color
from math import sin, cos, sqrt, pi, fmod


class Ellipse(GraphicPrimitive):
    """
    Primitive class for the ``Ellipse`` graphics type.  See ``ellipse?`` for
    information about actually plotting ellipses.

    INPUT:

    - ``x``, ``y`` -- coordinates of the center of the ellipse

    - ``r1``, ``r2`` -- radii of the ellipse

    - ``angle`` -- angle

    - ``options`` -- dictionary of options

    EXAMPLES:

    Note that this construction should be done using ``ellipse``::

        sage: from math import pi
        sage: from sage.plot.ellipse import Ellipse
        sage: Ellipse(0, 0, 2, 1, pi/4, {})
        Ellipse centered at (0.0, 0.0) with radii (2.0, 1.0) and angle 0.78539816339...
    """
    def __init__(self, x, y, r1, r2, angle, options):
        """
        Initialize base class ``Ellipse``.

        TESTS::

            sage: from sage.plot.ellipse import Ellipse
            sage: e = Ellipse(0, 0, 1, 1, 0, {})
            sage: print(loads(dumps(e)))
            Ellipse centered at (0.0, 0.0) with radii (1.0, 1.0) and angle 0.0
            sage: ellipse((0,0),0,1)
            Traceback (most recent call last):
            ...
            ValueError: both radii must be positive
        """
        self.x = float(x)
        self.y = float(y)
        self.r1 = float(r1)
        self.r2 = float(r2)
        if self.r1 <= 0 or self.r2 <= 0:
            raise ValueError("both radii must be positive")
        self.angle = fmod(angle, 2 * pi)
        if self.angle < 0:
            self.angle += 2 * pi
        GraphicPrimitive.__init__(self, options)

    def get_minmax_data(self):
        r"""
        Return a dictionary with the bounding box data.

        The bounding box is computed to be as minimal as possible.

        EXAMPLES:

        An example without an angle::

            sage: p = ellipse((-2, 3), 1, 2)
            sage: d = p.get_minmax_data()
            sage: d['xmin']
            -3.0
            sage: d['xmax']
            -1.0
            sage: d['ymin']
            1.0
            sage: d['ymax']
            5.0

        The same example with a rotation of angle `\pi/2`::

            sage: from math import pi
            sage: p = ellipse((-2, 3), 1, 2, pi/2)
            sage: d = p.get_minmax_data()
            sage: d['xmin']
            -4.0
            sage: d['xmax']
            0.0
            sage: d['ymin']
            2.0
            sage: d['ymax']
            4.0
        """
        from sage.plot.plot import minmax_data

        epsilon = 0.000001
        cos_angle = cos(self.angle)

        if abs(cos_angle) > 1-epsilon:
            xmax = self.r1
            ymax = self.r2
        elif abs(cos_angle) < epsilon:
            xmax = self.r2
            ymax = self.r1
        else:
            sin_angle = sin(self.angle)
            tan_angle = sin_angle / cos_angle
            sxmax = ((self.r2*tan_angle)/self.r1)**2
            symax = (self.r2/(self.r1*tan_angle))**2
            xmax = (
                abs(self.r1 * cos_angle / sqrt(sxmax+1.)) +
                abs(self.r2 * sin_angle / sqrt(1./sxmax+1.)))
            ymax = (
                abs(self.r1 * sin_angle / sqrt(symax+1.)) +
                abs(self.r2 * cos_angle / sqrt(1./symax+1.)))

        return minmax_data([self.x - xmax, self.x + xmax],
                           [self.y - ymax, self.y + ymax],
                           dict=True)

    def _allowed_options(self):
        """
        Return the allowed options for the ``Ellipse`` class.

        EXAMPLES::

            sage: p = ellipse((3, 3), 2, 1)
            sage: p[0]._allowed_options()['alpha']
            'How transparent the figure is.'
            sage: p[0]._allowed_options()['facecolor']
            '2D only: The color of the face as an RGB tuple.'
        """
        return {'alpha':'How transparent the figure is.',
                'fill': 'Whether or not to fill the ellipse.',
                'legend_label':'The label for this item in the legend.',
                'legend_color':'The color of the legend text.',
                'thickness':'How thick the border of the ellipse is.',
                'edgecolor':'2D only: The color of the edge as an RGB tuple.',
                'facecolor':'2D only: The color of the face as an RGB tuple.',
                'rgbcolor':'The color (edge and face) as an RGB tuple.',
                'hue':'The color given as a hue.',
                'zorder':'2D only: The layer level in which to draw',
                'linestyle':"2D only: The style of the line, which is one of "
                "'dashed', 'dotted', 'solid', 'dashdot', or '--', ':', '-', '-.', "
                "respectively."}

    def _repr_(self):
        """
        String representation of ``Ellipse`` primitive.

        TESTS::

            sage: from sage.plot.ellipse import Ellipse
            sage: Ellipse(0,0,2,1,0,{})._repr_()
            'Ellipse centered at (0.0, 0.0) with radii (2.0, 1.0) and angle 0.0'
        """
        return "Ellipse centered at ({}, {}) with radii ({}, {}) and angle {}".format(self.x, self.y, self.r1, self.r2, self.angle)

    def _render_on_subplot(self, subplot):
        """
        Render this ellipse in a subplot.  This is the key function that
        defines how this ellipse graphics primitive is rendered in matplotlib's
        library.

        TESTS::

            sage: from math import pi
            sage: ellipse((0,0),3,1,pi/6,fill=True,alpha=0.3)
            Graphics object consisting of 1 graphics primitive

        ::

            sage: ellipse((3,2),1,2)
            Graphics object consisting of 1 graphics primitive
        """
        import matplotlib.patches as patches
        from sage.plot.misc import get_matplotlib_linestyle

        options = self.options()
        p = patches.Ellipse(
                (self.x,self.y),
                self.r1*2.,self.r2*2.,
                angle=self.angle/pi*180.)
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

    def plot3d(self):
        r"""
        Plotting in 3D is not implemented.

        TESTS::

            sage: from sage.plot.ellipse import Ellipse
            sage: Ellipse(0,0,2,1,pi/4,{}).plot3d()                                     # needs sage.symbolic
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


@rename_keyword(color='rgbcolor')
@options(alpha=1, fill=False, thickness=1, edgecolor='blue', facecolor='blue', linestyle='solid', zorder=5,
         aspect_ratio=1.0, legend_label=None, legend_color=None)
def ellipse(center, r1, r2, angle=0, **options):
    """
    Return an ellipse centered at a point center = ``(x,y)`` with radii =
    ``r1,r2`` and angle ``angle``.  Type ``ellipse.options`` to see all
    options.

    INPUT:

    - ``center`` -- 2-tuple of real numbers; coordinates of the center

    - ``r1``, ``r2`` -- positive real numbers; the radii of the ellipse

    - ``angle`` -- real number (default: 0) -- the angle between the first axis
      and the horizontal

    OPTIONS:

    - ``alpha`` -- (default: 1) transparency

    - ``fill`` -- (default: ``False``) whether to fill the ellipse or not

    - ``thickness`` -- (default: 1) thickness of the line

    - ``linestyle`` -- (default: ``'solid'``) the style of the line, which is one
      of ``'dashed'``, ``'dotted'``, ``'solid'``, ``'dashdot'``, or ``'--'``,
      ``':'``, ``'-'``, ``'-.'``,  respectively

    - ``edgecolor`` -- (default: ``'black'``) color of the contour

    - ``facecolor`` -- (default: ``'red'``) color of the filling

    - ``rgbcolor`` -- 2D or 3D plotting.  This option overrides
      ``edgecolor`` and ``facecolor`` for 2D plotting

    - ``legend_label`` -- the label for this item in the legend

    - ``legend_color`` -- the color for the legend label

    EXAMPLES:

    An ellipse centered at (0,0) with major and minor axes of lengths 2 and 1.
    Note that the default color is blue::

        sage: ellipse((0,0),2,1)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        E=ellipse((0,0),2,1)
        sphinx_plot(E)

    More complicated examples with tilted axes and drawing options::

        sage: from math import pi
        sage: ellipse((0,0), 3, 1, pi/6, fill=True, alpha=0.3, linestyle='dashed')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        E = ellipse((0,0),3,1,pi/6,fill=True,alpha=0.3,linestyle='dashed')
        sphinx_plot(E)

    other way to indicate dashed linestyle::

        sage: ellipse((0,0),3,1,pi/6,fill=True,alpha=0.3,linestyle='--')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        E =ellipse((0,0),3,1,pi/6,fill=True,alpha=0.3,linestyle='--')
        sphinx_plot(E)

    with colors ::

        sage: ellipse((0,0),3,1,pi/6,fill=True,edgecolor='black',facecolor='red')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        E=ellipse((0,0),3,1,pi/6,fill=True,edgecolor='black',facecolor='red')
        sphinx_plot(E)

    We see that ``rgbcolor`` overrides these other options, as this plot
    is green::

        sage: ellipse((0,0),3,1,pi/6,fill=True,edgecolor='black',facecolor='red',rgbcolor='green')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        E=ellipse((0,0),3,1,pi/6,fill=True,edgecolor='black',facecolor='red',rgbcolor='green')
        sphinx_plot(E)

    The default aspect ratio for ellipses is 1.0::

        sage: ellipse((0,0),2,1).aspect_ratio()
        1.0

    One cannot yet plot ellipses in 3D::

        sage: ellipse((0,0,0),2,1)
        Traceback (most recent call last):
        ...
        NotImplementedError: plotting ellipse in 3D is not implemented

    We can also give ellipses a legend::

        sage: ellipse((0,0),2,1,legend_label="My ellipse", legend_color='green')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        E=ellipse((0,0),2,1,legend_label="My ellipse", legend_color='green')
        sphinx_plot(E)

    TESTS:

    Verify that :issue:`36153` does not arise::

        sage: E = ellipse((0,0), 2, 1, legend_label='test')
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
    g.add_primitive(Ellipse(center[0],center[1],r1,r2,angle,options))
    if options['legend_label']:
        g.legend(True)
        g._legend_colors = [options['legend_color']]
    if len(center) == 2:
        return g
    elif len(center) == 3:
        raise NotImplementedError("plotting ellipse in 3D is not implemented")
