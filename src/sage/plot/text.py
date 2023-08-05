"""
Text in plots
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
from sage.plot.primitive import GraphicPrimitive
from sage.misc.decorators import options, rename_keyword
from sage.plot.colors import to_mpl_color


class Text(GraphicPrimitive):
    """
    Base class for Text graphics primitive.

    TESTS:

    We test creating some text::

        sage: text("I like Fibonacci",(3,5))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(text("I like Fibonacci",(3,5)))

    """
    def __init__(self, string, point, options):
        """
        Initialize base class Text.

        EXAMPLES::

            sage: T = text("I like Fibonacci", (3,5))
            sage: t = T[0]
            sage: t.string
            'I like Fibonacci'
            sage: t.x
            3.0
            sage: t.options()['fontsize']
            10
        """
        self.string = string
        self.x = float(point[0])
        self.y = float(point[1])
        GraphicPrimitive.__init__(self, options)

    def get_minmax_data(self):
        """
        Return a dictionary with the bounding box data. Notice
        that, for text, the box is just the location itself.

        EXAMPLES::

            sage: T = text("Where am I?",(1,1))
            sage: t=T[0]
            sage: t.get_minmax_data()['ymin']
            1.0
            sage: t.get_minmax_data()['ymax']
            1.0
        """
        from sage.plot.plot import minmax_data
        return minmax_data([self.x], [self.y], dict=True)

    def _repr_(self):
        """
        String representation of Text primitive.

        EXAMPLES::

            sage: T = text("I like cool constants", (pi,e))                             # needs sage.symbolic
            sage: t = T[0];t                                                            # needs sage.symbolic
            Text 'I like cool constants' at the point (3.1415926535...,2.7182818284...)
        """
        return "Text '%s' at the point (%s,%s)" % (self.string, self.x, self.y)

    def _allowed_options(self):
        """
        Return the allowed options for the Text class.

        EXAMPLES::

            sage: T = text("ABC",(1,1),zorder=3)
            sage: T[0]._allowed_options()['fontsize']
            "How big the text is. Either the size in points or a relative size, e.g. 'smaller', 'x-large', etc"
            sage: T[0]._allowed_options()['zorder']
            'The layer level in which to draw'
            sage: T[0]._allowed_options()['rotation']
            'How to rotate the text: angle in degrees, vertical, horizontal'

        """
        return {'fontsize': 'How big the text is. Either the size in points or a relative size, e.g. \'smaller\', \'x-large\', etc',
                'fontstyle': 'A string either \'normal\', \'italic\' or \'oblique\'',
                'fontweight': 'A numeric value in the range 0-1000 or a string'
                              '\'ultralight\', \'light\', \'normal\', \'regular\', \'book\','
                              '\'medium\', \'roman\', \'semibold\', \'demibold\', \'demi\','
                              '\'bold,\', \'heavy\', \'extra bold\', \'black\'',
                'rgbcolor': 'The color as an RGB tuple',
                'background_color': 'The background color',
                'bounding_box': 'A dictionary specifying a bounding box',
                'hue': 'The color given as a hue',
                'alpha': 'A float (0.0 transparent through 1.0 opaque)',
                'axis_coords': 'If True use axis coordinates: (0,0) lower left and (1,1) upper right',
                'rotation': 'How to rotate the text: angle in degrees, vertical, horizontal',
                'vertical_alignment': 'How to align vertically: top, center, bottom',
                'horizontal_alignment': 'How to align horizontally: left, center, right',
                'zorder': 'The layer level in which to draw',
                'clip': 'Whether to clip or not'}

    def _plot3d_options(self, options=None):
        """
        Translate 2D plot options into 3D plot options.

        EXAMPLES::

            sage: T = text("ABC",(1,1))
            sage: t = T[0]
            sage: t.options()['rgbcolor']
            (0.0, 0.0, 1.0)
            sage: s=t.plot3d()
            sage: s.jmol_repr(s.testing_render_params())[0][1]
            'color atom  [0,0,255]'

        """
        if options is None:
            options = dict(self.options())
        options_3d = {}
        for s in ['fontfamily', 'fontsize', 'fontstyle', 'fontweight']:
            if s in options:
                options_3d[s] = options.pop(s)
        # TODO: figure out how to implement rather than ignore
        for s in ['axis_coords', 'clip', 'horizontal_alignment',
                  'rotation', 'vertical_alignment']:
            if s in options:
                del options[s]
        options_3d.update(GraphicPrimitive._plot3d_options(self, options))
        return options_3d

    def plot3d(self, **kwds):
        """
        Plot 2D text in 3D.

        EXAMPLES::

            sage: T = text("ABC", (1, 1))
            sage: t = T[0]
            sage: s = t.plot3d()
            sage: s.jmol_repr(s.testing_render_params())[0][2]
            'label "ABC"'
            sage: s._trans
            (1.0, 1.0, 0)

        """
        from sage.plot.plot3d.shapes2 import text3d
        options = self._plot3d_options()
        options.update(kwds)
        return text3d(self.string, (self.x, self.y, 0), **options)

    def _render_on_subplot(self, subplot):
        """
        TESTS::

            sage: t1 = text("Hello",(1,1), vertical_alignment="top", fontsize=30, rgbcolor='black')
            sage: t2 = text("World", (1,1), horizontal_alignment="left", fontsize=20, zorder=-1)
            sage: t1 + t2   # render the sum
            Graphics object consisting of 2 graphics primitives

        """
        options = self.options()
        opts = {}
        opts['color'] = options['rgbcolor']
        opts['verticalalignment'] = options['vertical_alignment']
        opts['horizontalalignment'] = options['horizontal_alignment']
        if 'background_color' in options:
            opts['backgroundcolor'] = options['background_color']
        if 'fontweight' in options:
            opts['fontweight'] = options['fontweight']
        if 'alpha' in options:
            opts['alpha'] = options['alpha']
        if 'fontstyle' in options:
            opts['fontstyle'] = options['fontstyle']
        if 'bounding_box' in options:
            opts['bbox'] = options['bounding_box']
        if 'zorder' in options:
            opts['zorder'] = options['zorder']
        if options['axis_coords']:
            opts['transform'] = subplot.transAxes
        if 'fontsize' in options:
            val = options['fontsize']
            if isinstance(val, str):
                opts['fontsize'] = val
            else:
                opts['fontsize'] = int(val)
        if 'rotation' in options:
            val = options['rotation']
            if isinstance(val, str):
                opts['rotation'] = options['rotation']
            else:
                opts['rotation'] = float(options['rotation'])

        p = subplot.text(self.x, self.y, self.string, clip_on=options['clip'], **opts)
        if not options['clip']:
            self._bbox_extra_artists = [p]


@rename_keyword(color='rgbcolor')
@options(fontsize=10, rgbcolor=(0,0,1), horizontal_alignment='center',
         vertical_alignment='center', axis_coords=False, clip=False)
def text(string, xy, **options):
    r"""
    Return a 2D text graphics object at the point `(x, y)`.

    Type ``text.options`` for a dictionary of options for 2D text.

    2D OPTIONS:

    - ``fontsize`` - How big the text is. Either an integer that
      specifies the size in points or a string which specifies a size (one of
      'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large')

    - ``fontstyle`` - A string either 'normal', 'italic' or 'oblique'

    - ``fontweight`` - A numeric value in the range 0-1000 or a string (one of
      'ultralight', 'light', 'normal', 'regular', 'book',' 'medium', 'roman',
      'semibold', 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black')

    - ``rgbcolor`` - The color as an RGB tuple

    - ``hue`` - The color given as a hue

    - ``alpha`` - A float (0.0 transparent through 1.0 opaque)

    - ``background_color`` - The background color

    - ``rotation`` - How to rotate the text: angle in degrees, vertical, horizontal

    - ``vertical_alignment`` - How to align vertically: top, center, bottom

    - ``horizontal_alignment`` - How to align horizontally: left, center, right

    - ``zorder`` - The layer level in which to draw

    - ``clip`` - (default: False) Whether to clip or not

    - ``axis_coords`` - (default: False) If True, use axis coordinates, so that
      (0,0) is the lower left and (1,1) upper right, regardless of the x and y
      range of plotted values.

    - ``bounding_box`` - A dictionary specifying a bounding box. Currently the text location.

    EXAMPLES::

        sage: text("Sage graphics are really neat because they use matplotlib!", (2,12))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        t = "Sage graphics are really neat because they use matplotlib!"
        sphinx_plot(text(t,(2,12)))

    Larger font, bold, colored red and transparent text::

        sage: text("I had a dream!", (2,12), alpha=0.3,
        ....:      fontsize='large', fontweight='bold', color='red')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(text("I had a dream!", (2,12), alpha=0.3, fontsize='large', fontweight='bold', color='red'))

    By setting ``horizontal_alignment`` to 'left' the text is guaranteed to be
    in the lower left no matter what::

        sage: text("I got a horse and he lives in a tree", (0,0),
        ....:      axis_coords=True, horizontal_alignment='left')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        t = "I got a horse and he lives in a tree"
        sphinx_plot(text(t, (0,0), axis_coords=True, horizontal_alignment='left'))

    Various rotations::

        sage: text("noitator", (0,0), rotation=45.0,
        ....:      horizontal_alignment='left', vertical_alignment='bottom')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(text("noitator", (0,0), rotation=45.0, horizontal_alignment='left', vertical_alignment='bottom'))

    ::

        sage: text("Sage is really neat!!",(0,0), rotation="vertical")
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(text("Sage is really neat!!",(0,0), rotation="vertical"))

    You can also align text differently::

        sage: t1 = text("Hello", (1,1), vertical_alignment="top")
        sage: t2 = text("World", (1,0.5), horizontal_alignment="left")
        sage: t1 + t2   # render the sum
        Graphics object consisting of 2 graphics primitives

    .. PLOT::

        t1 = text("Hello",(1,1), vertical_alignment="top")
        t2 = text("World", (1,0.5), horizontal_alignment="left")
        sphinx_plot(t1 + t2)

    You can save text as part of PDF output::

        sage: import tempfile
        sage: with tempfile.NamedTemporaryFile(suffix=".pdf") as f:
        ....:     text("sage", (0,0), rgbcolor=(0,0,0)).save(f.name)

    Some examples of bounding box::

        sage: bbox = {'boxstyle': "rarrow,pad=0.3", 'fc': "cyan", 'ec': "b", 'lw': 2}
        sage: text("I feel good", (1,2), bounding_box=bbox)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

         bbox = {'boxstyle':"rarrow,pad=0.3", 'fc':"cyan", 'ec':"b", 'lw':2}
         sphinx_plot(text("I feel good", (1,2), bounding_box=bbox))

    ::

        sage: text("So good", (0,0), bounding_box={'boxstyle': 'round', 'fc': 'w'})
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        bbox = {'boxstyle':'round', 'fc':'w'}
        sphinx_plot(text("So good", (0,0), bounding_box=bbox))

    The possible options of the bounding box are 'boxstyle' (one of 'larrow',
    'rarrow', 'round', 'round4', 'roundtooth', 'sawtooth', 'square'), 'fc' or
    'facecolor', 'ec' or 'edgecolor', 'ha' or 'horizontalalignment', 'va' or
    'verticalalignment', 'lw' or 'linewidth'.

    A text with a background color::

        sage: text("So good", (-2,2), background_color='red')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(text("So good", (-2,2), background_color='red'))

    Use dollar signs for LaTeX and raw strings to avoid having to
    escape backslash characters::

        sage: A = arc((0, 0), 1, sector=(0.0, RDF.pi()))
        sage: a = sqrt(1./2.)
        sage: PQ = point2d([(-a, a), (a, a)])
        sage: botleft = dict(horizontal_alignment='left', vertical_alignment='bottom')
        sage: botright = dict(horizontal_alignment='right', vertical_alignment='bottom')
        sage: tp = text(r'$z_P = e^{3i\pi/4}$',
        ....:           (-a, a), **botright)
        sage: tq = text(r'$Q = (\frac{\sqrt{2}}{2}, \frac{\sqrt{2}}{2})$',
        ....:           (a, a), **botleft)
        sage: A + PQ + tp + tq
        Graphics object consisting of 4 graphics primitives

    .. PLOT::

        A = arc((0, 0), 1, sector=(0.0, RDF.pi()))
        a = sqrt(1./2.)
        PQ = point2d([(-a, a), (a, a)])
        botleft = dict(horizontal_alignment='left', vertical_alignment='bottom')
        botright = dict(horizontal_alignment='right', vertical_alignment='bottom')
        tp = text(r'$z_P = e^{3i\pi/4}$', (-a, a), **botright)
        tq = text(r'$Q = (\frac{\sqrt{2}}{2}, \frac{\sqrt{2}}{2})$', (a, a), **botleft)
        sphinx_plot(A + PQ + tp + tq)

    Text coordinates must be 2D, an error is raised if 3D coordinates are passed::

        sage: t = text("hi", (1, 2, 3))
        Traceback (most recent call last):
        ...
        ValueError: use text3d instead for text in 3d

    Use the :func:`text3d <sage.plot.plot3d.shapes2.text3d>` function for 3D text::

        sage: t = text3d("hi", (1, 2, 3))

    Or produce 2D text with coordinates `(x, y)` and plot it in 3D (at `z = 0`)::

        sage: t = text("hi", (1, 2))
        sage: t.plot3d()  # text at position (1, 2, 0)
        Graphics3d Object

    Extra options will get passed on to ``show()``, as long as they are valid. Hence this ::

        sage: text("MATH IS AWESOME", (0, 0), fontsize=40, axes=False)
        Graphics object consisting of 1 graphics primitive

    is equivalent to ::

        sage: text("MATH IS AWESOME", (0, 0), fontsize=40).show(axes=False)
    """
    try:
        x, y = xy
    except ValueError:
        if isinstance(xy, (list, tuple)) and len(xy) == 3:
            raise ValueError("use text3d instead for text in 3d")
        raise
    from sage.plot.all import Graphics
    options['rgbcolor'] = to_mpl_color(options['rgbcolor'])
    point = (float(x), float(y))
    g = Graphics()
    g._set_extra_kwds(Graphics._extract_kwds_for_show(options, ignore='fontsize'))
    g.add_primitive(Text(string, point, options))
    return g
