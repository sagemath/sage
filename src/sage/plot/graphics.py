r"""
Graphics objects

This file contains the definition of the class :class:`Graphics`.
Usually, you don't call the constructor of this class directly
(although you can do it), you would use :func:`plot` instead.

AUTHORS:

- Jeroen Demeyer (2012-04-19): split off this file from plot.py (:issue:`12857`)

- Punarbasu Purkayastha (2012-05-20): Add logarithmic scale (:issue:`4529`)

- Emily Chen (2013-01-05): Add documentation for
  :meth:`~sage.plot.graphics.Graphics.show` figsize parameter (:issue:`5956`)

- Eric Gourgoulhon (2015-03-19): Add parameter axes_labels_size (:issue:`18004`)

- Eric Gourgoulhon (2019-05-24): :class:`~sage.plot.multigraphics.GraphicsArray`
  moved to new module :mod:`~sage.plot.multigraphics`; various improvements and
  fixes in :meth:`Graphics.matplotlib` and ``Graphics._set_scale``; new method
  :meth:`Graphics.inset`
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

import os
from numbers import Integral
from collections.abc import Iterable
from math import isnan
import sage.misc.verbose
from sage.misc.temporary_file import tmp_filename
from sage.misc.fast_methods import WithEqualityById
from sage.structure.sage_object import SageObject
from sage.misc.decorators import suboptions
from .colors import rgbcolor

ALLOWED_EXTENSIONS = ['.eps', '.pdf', '.pgf', '.png', '.ps', '.sobj', '.svg']
DEFAULT_DPI = 100


# If do_verify is True, options are checked when drawing a
# GraphicsPrimitive.  See primitive.py
do_verify = True


def is_Graphics(x):
    """
    Return ``True`` if `x` is a Graphics object.

    EXAMPLES::

        sage: from sage.plot.graphics import is_Graphics
        sage: is_Graphics(1)
        doctest:warning...
        DeprecationWarning: The function is_Graphics is deprecated;
        use 'isinstance(..., Graphics)' instead.
        See https://github.com/sagemath/sage/issues/38184 for details.
        False
        sage: is_Graphics(disk((0.0, 0.0), 1, (0, pi/2)))                               # needs sage.symbolic
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(38184,
                "The function is_Graphics is deprecated; "
                "use 'isinstance(..., Graphics)' instead.")
    return isinstance(x, Graphics)


def _parse_figsize(figsize):
    r"""
    Helper function to get a figure size in matplotlib format.

    INPUT:

    - ``figsize`` -- width or [width, height] in inches; if only the width is
      provided, the height is computed from matplotlib's default aspect ratio

    OUTPUT:

    - a pair of ``float``'s representing ``(width, height)``

    EXAMPLES::

        sage: from sage.plot.graphics import _parse_figsize
        sage: _parse_figsize([5, 4])
        (5.0, 4.0)

    The default aspect ratio is 4/3::

        sage: _parse_figsize(5)  # tol 1.0e-13
        (5.0, 3.75)
    """
    from matplotlib import rcParams
    if isinstance(figsize, (list, tuple)):
        # figsize should be a pair of positive numbers
        if len(figsize) != 2:
            raise ValueError("figsize should be a positive number or a list "
                             f"of two positive numbers, not {figsize}")
        figsize = (float(figsize[0]), float(figsize[1]))  # floats for mpl
        if not (figsize[0] > 0 and figsize[1] > 0):
            raise ValueError("figsize should be positive numbers, "
                             f"not {figsize[0]} and {figsize[1]}")
    else:
        # in this case, figsize is a single number representing the width and
        # should be positive
        try:
            figsize = float(figsize)  # to pass to mpl
        except TypeError:
            raise TypeError(f"figsize should be a positive number, not {figsize}")
        if figsize > 0:
            default_width, default_height = rcParams['figure.figsize']
            figsize = (figsize, default_height * figsize / default_width)
        else:
            raise ValueError(f"figsize should be positive, not {figsize}")
    return figsize


class Graphics(WithEqualityById, SageObject):
    """
    The Graphics object is an empty list of graphics objects. It is
    useful to use this object when initializing a for loop where
    different graphics object will be added to the empty object.

    EXAMPLES::

        sage: G = Graphics(); print(G)
        Graphics object consisting of 0 graphics primitives
        sage: c = circle((1,1), 1)
        sage: G += c; print(G)
        Graphics object consisting of 1 graphics primitive

    Here we make a graphic of embedded isosceles triangles, coloring
    each one with a different color as we go::

        sage: h = 10; c = 0.4; p = 0.5
        sage: G = Graphics()
        sage: for x in srange(1, h+1):                                                  # needs sage.symbolic
        ....:     l = [[0,x*sqrt(3)],[-x/2,-x*sqrt(3)/2],[x/2,-x*sqrt(3)/2],[0,x*sqrt(3)]]
        ....:     G += line(l, color=hue(c + p*(x/h)))
        sage: G.show(figsize=[5,5])                                                     # needs sage.symbolic

    We can change the scale of the axes in the graphics before displaying.::

        sage: G = plot(exp, 1, 10)              # long time                             # needs sage.symbolic
        sage: G.show(scale='semilogy')          # long time                             # needs sage.symbolic

    TESTS:

    From :issue:`4604`, ensure Graphics can handle 3d objects::

        sage: g = Graphics()
        sage: g += sphere((1, 1, 1), 2)
        sage: g.show()

    We check that graphics can be pickled (we can't use equality on
    graphics so we just check that the load/dump cycle gives a
    :class:`Graphics` instance)::

        sage: g = Graphics()
        sage: g2 = loads(dumps(g))
        sage: g2.show()

    ::

        sage: isinstance(g2, Graphics)
        True

        sage: hash(Graphics()) # random
        42

    .. automethod:: _rich_repr_
    """

    def __init__(self):
        """
        Create a new empty Graphics objects with all the defaults.

        EXAMPLES::

            sage: G = Graphics()
        """
        self._axes_color = (0, 0, 0)
        self._axes_label_color = (0, 0, 0)
        self._axes_width = 0.8
        self._bbox_extra_artists = []
        self._extra_kwds = {}
        self._fontsize = 10
        self._axes_labels_size = 1.6
        self._legend_colors = []
        self._legend_opts = {}
        self._objects = []
        self._show_axes = True
        self._show_legend = False
        self._tick_label_color = (0, 0, 0)

    def set_aspect_ratio(self, ratio):
        """
        Set the aspect ratio, which is the ratio of height and width
        of a unit square (i.e., height/width of a unit square), or
        'automatic' (expand to fill the figure).

        INPUT:

        - ``ratio`` -- a positive real number or 'automatic'

        EXAMPLES: We create a plot of the upper half of a circle, but it
        doesn't look round because the aspect ratio is off::

            sage: P = plot(sqrt(1-x^2),(x,-1,1)); P                                     # needs sage.symbolic
            Graphics object consisting of 1 graphics primitive

        So we set the aspect ratio and now it is round::

            sage: P.set_aspect_ratio(1)                                                 # needs sage.symbolic
            sage: P.aspect_ratio()                                                      # needs sage.symbolic
            1.0
            sage: P                                                                     # needs sage.symbolic
            Graphics object consisting of 1 graphics primitive

        Note that the aspect ratio is inherited upon addition (which takes
        the max of aspect ratios of objects whose aspect ratio has been
        set)::

            sage: P + plot(sqrt(4-x^2),(x,-2,2))                                        # needs sage.symbolic
            Graphics object consisting of 2 graphics primitives

        In the following example, both plots produce a circle that looks
        twice as tall as wide::

            sage: Q = circle((0,0), 0.5); Q.set_aspect_ratio(2)
            sage: (P + Q).aspect_ratio(); P + Q                                         # needs sage.symbolic
            2.0
            Graphics object consisting of 2 graphics primitives
            sage: (Q + P).aspect_ratio(); Q + P                                         # needs sage.symbolic
            2.0
            Graphics object consisting of 2 graphics primitives
        """
        if ratio != 'auto' and ratio != 'automatic':
            ratio = float(ratio)
            if ratio <= 0:
                raise ValueError("the aspect ratio must be positive or 'automatic'")
        else:
            ratio = 'automatic'
        self._extra_kwds['aspect_ratio'] = ratio

    def aspect_ratio(self):
        """
        Get the current aspect ratio, which is the ratio of height to
        width of a unit square, or ``'automatic'``.

        OUTPUT: a positive float (height/width of a unit square), or ``'automatic'``
        (expand to fill the figure).

        EXAMPLES:

        The default aspect ratio for a new blank :class:`Graphics` object is ``'automatic'``::

            sage: P = Graphics()
            sage: P.aspect_ratio()
            'automatic'

        The aspect ratio can be explicitly set different from the object's default::

            sage: P = circle((1,1), 1)
            sage: P.aspect_ratio()
            1.0
            sage: P.set_aspect_ratio(2)
            sage: P.aspect_ratio()
            2.0
            sage: P.set_aspect_ratio('automatic')
            sage: P.aspect_ratio()
            'automatic'
        """
        return self._extra_kwds.get('aspect_ratio', 'automatic')

    def legend(self, show=None):
        r"""
        Set whether or not the legend is shown by default.

        INPUT:

        - ``show`` -- (default: ``None``) a boolean

        If called with no input, return the current legend setting.

        EXAMPLES:

        By default no legend is displayed::

            sage: P = plot(sin)                                                         # needs sage.symbolic
            sage: P.legend()                                                            # needs sage.symbolic
            False

        But if we put a label then the legend is shown::

            sage: P = plot(sin, legend_label='sin')                                     # needs sage.symbolic
            sage: P.legend()                                                            # needs sage.symbolic
            True

        We can turn it on or off::

            sage: # needs sage.symbolic
            sage: P.legend(False)
            sage: P.legend()
            False
            sage: P.legend(True)
            sage: P  # show with the legend
            Graphics object consisting of 1 graphics primitive
        """
        if show is None:
            return self._show_legend
        else:
            self._show_legend = bool(show)

    def set_legend_options(self, **kwds):
        r"""
        Set various legend options.

        INPUT:

        - ``title`` -- (default: ``None``) string, the legend title

        - ``ncol`` -- (default: 1) positive integer, the number of columns

        - ``columnspacing`` -- (default: ``None``) the spacing between columns

        - ``borderaxespad`` -- (default: ``None``) float, length between the axes and the legend

        - ``back_color`` -- (default: ``'white'``) this parameter can be a string
          denoting a color or an RGB tuple. The string can be a color name
          as in ('red', 'green', 'yellow', ...) or a floating point number
          like '0.8' which gets expanded to (0.8, 0.8, 0.8). The
          tuple form is just a floating point RGB tuple with all values ranging
          from 0 to 1.

        - ``handlelength`` -- (default: 0.05) float, the length of the legend handles

        - ``handletextpad`` -- (default: 0.5) float, the pad between the legend handle and text

        - ``labelspacing`` -- (default: 0.02) float, vertical space between legend entries

        - ``loc`` -- (default: ``'best'``) may be a string, an integer or a tuple. String or
              integer inputs must be one of the following:

          - 0, 'best'

          - 1, 'upper right'

          - 2, 'upper left'

          - 3, 'lower left'

          - 4, 'lower right'

          - 5, 'right'

          - 6, 'center left'

          - 7, 'center right'

          - 8, 'lower center'

          - 9, 'upper center'

          - 10, 'center'

          - Tuple arguments represent an absolute (x, y) position on the plot
            in axes coordinates (meaning from 0 to 1 in each direction).

        - ``markerscale`` -- (default: 0.6) float, how much to scale the markers in the legend

        - ``numpoints`` -- (default: 2) integer, the number of points in the legend for line

        - ``borderpad`` -- (default: 0.6) float, the fractional whitespace inside the legend border
          (between 0 and 1)

        - ``font_family`` -- (default: ``'sans-serif'``) string, one of
          ``'serif'``, ``'sans-serif'``, ``'cursive'``, ``'fantasy'``,
          ``'monospace'``

        - ``font_style`` -- (default: ``'normal'``) string, one of
          ``'normal'``, ``'italic'``, ``'oblique'``

        - ``font_variant`` -- (default: ``'normal'``) string, one of
          ``'normal'``, ``'small-caps'``

        - ``font_weight`` -- (default: ``'medium'``) string, one of
          ``'black'``, ``'extra bold'``, ``'bold'``, ``'semibold'``,
          ``'medium'``, ``'normal'``, ``'light'``

        - ``font_size`` -- (default: ``'medium'``) string, one of
          ``'xx-small'``, ``'x-small'``, ``'small'``, ``'medium'``,
          ``'large'``, ``'x-large'``, ``'xx-large'``, or an absolute font size
          (e.g. 12)

        - ``shadow`` -- boolean (default: ``True``); draw a shadow behind the legend

        - ``fancybox`` -- boolean (default: ``False``); if
          ``True``, draws a frame with a round fancybox

        These are all keyword arguments.

        OUTPUT: a dictionary of all current legend options

        EXAMPLES:

        By default, no options are set::

            sage: p = plot(tan, legend_label='tan')                                     # needs sage.symbolic
            sage: p.set_legend_options()                                                # needs sage.symbolic
            {}

        We build a legend without a shadow::

            sage: p.set_legend_options(shadow=False)                                    # needs sage.symbolic
            sage: p.set_legend_options()['shadow']                                      # needs sage.symbolic
            False

        To set the legend position to the center of the plot, all these
        methods are roughly equivalent::

            sage: p.set_legend_options(loc='center'); p                                 # needs sage.symbolic
            Graphics object consisting of 1 graphics primitive

        ::

            sage: p.set_legend_options(loc=10); p                                       # needs sage.symbolic
            Graphics object consisting of 1 graphics primitive

        ::

            sage: p.set_legend_options(loc=(0.5,0.5)); p  # aligns the bottom of the box to the center                # needs sage.symbolic
            Graphics object consisting of 1 graphics primitive

        The parameters ``loc`` and ``borderaxespad`` can be altered
        in order to place the legend below the x-axis label or to
        the left of the y-axis label::

            sage: p = line([(0, 0), (1, 1)], legend_label='test')
            sage: p.axes_labels(['X-Label', 'Y-Label'])  # adding labels for axes
            sage: p.set_legend_options(loc=8, borderaxespad=-7.5-0.01*p.fontsize())
            sage: p
            Graphics object consisting of 1 graphics primitive
        """
        if len(kwds) == 0:
            return self._legend_opts
        else:
            self._legend_opts.update(kwds)

    def get_axes_range(self):
        """
        Return a dictionary of the range of the axes for this graphics
        object.  This is fall back to the ranges in ``get_minmax_data()`` for
        any value which the user has not explicitly set.

        .. warning::

           Changing the dictionary returned by this function does not
           change the axes range for this object.  To do that, use the
           :meth:`set_axes_range` method.

        EXAMPLES::

            sage: L = line([(1,2), (3,-4), (2, 5), (1,2)])
            sage: list(sorted(L.get_axes_range().items()))
            [('xmax', 3.0), ('xmin', 1.0), ('ymax', 5.0), ('ymin', -4.0)]
            sage: L.set_axes_range(xmin=-1)
            sage: list(sorted(L.get_axes_range().items()))
            [('xmax', 3.0), ('xmin', -1.0), ('ymax', 5.0), ('ymin', -4.0)]
        """
        axes_range = self.get_minmax_data()
        axes_range.update(self._get_axes_range_dict())
        return axes_range

    def set_axes_range(self, xmin=None, xmax=None, ymin=None, ymax=None):
        """
        Set the ranges of the `x` and `y` axes.

        INPUT:

        - ``xmin``, ``xmax``, ``ymin``, ``ymax`` -- floats

        EXAMPLES::

            sage: L = line([(1,2), (3,-4), (2, 5), (1,2)])
            sage: L.set_axes_range(-1, 20, 0, 2)
            sage: d = L.get_axes_range()
            sage: d['xmin'], d['xmax'], d['ymin'], d['ymax']
            (-1.0, 20.0, 0.0, 2.0)
        """
        l = locals()
        axes_range = self._get_axes_range_dict()
        for name in ['xmin', 'xmax', 'ymin', 'ymax']:
            if l[name] is not None:
                axes_range[name] = float(l[name])

    axes_range = set_axes_range

    def _get_axes_range_dict(self):
        """
        Return the underlying dictionary used to store the user's
        custom ranges for the axes on this object.

        EXAMPLES::

            sage: L = line([(1,2), (3,-4), (2, 5), (1,2)])
            sage: L._get_axes_range_dict()
            {}
            sage: L.set_axes_range(xmin=-1)
            sage: L._get_axes_range_dict()
            {'xmin': -1.0}
        """
        try:
            return self._axes_range
        except AttributeError:
            self._axes_range = {}
            return self._axes_range

    def set_flip(self, flip_x=None, flip_y=None):
        """
        Set the flip options for this graphics object.

        INPUT:

        - ``flip_x`` -- boolean (default: ``None``); if not ``None``, set the
          ``flip_x`` option to this value
        - ``flip_y`` -- boolean (default: ``None``); if not ``None``, set the
          ``flip_y`` option to this value

        EXAMPLES::

            sage: L = line([(1, 0), (2, 3)])
            sage: L.set_flip(flip_y=True)
            sage: L.flip()
            (False, True)
            sage: L.set_flip(True, False)
            sage: L.flip()
            (True, False)
        """
        if flip_x is not None:
            self._extra_kwds['flip_x'] = flip_x
        if flip_y is not None:
            self._extra_kwds['flip_y'] = flip_y

    def flip(self, flip_x=False, flip_y=False):
        """
        Get the flip options and optionally mirror this graphics object.

        INPUT:

        - ``flip_x`` -- boolean (default: ``False``); if ``True``, replace the
          current ``flip_x`` option by its opposite
        - ``flip_y`` -- boolean (default: ``False``); if ``True``, replace the
          current ``flip_y`` option by its opposite

        OUTPUT: a tuple containing the new flip options

        EXAMPLES:

        When called without arguments, this just returns the current flip
        options::

            sage: L = line([(1, 0), (2, 3)])
            sage: L.flip()
            (False, False)

        Otherwise, the specified options are changed and the new options are
        returned::

            sage: L.flip(flip_y=True)
            (False, True)
            sage: L.flip(True, True)
            (True, False)
        """
        a = self._extra_kwds.get('flip_x', self.SHOW_OPTIONS['flip_x'])
        b = self._extra_kwds.get('flip_y', self.SHOW_OPTIONS['flip_y'])
        if flip_x:
            a = not a
            self._extra_kwds['flip_x'] = a
        if flip_y:
            b = not b
            self._extra_kwds['flip_y'] = b
        return (a, b)

    def fontsize(self, s=None):
        """
        Set the font size of axes labels and tick marks.

        Note that the relative size of the axes labels font w.r.t. the tick
        marks font can be adjusted via :meth:`axes_labels_size`.

        INPUT:

        - ``s`` -- integer, a font size in points


        If called with no input, return the current fontsize.

        EXAMPLES::

            sage: L = line([(1,2), (3,-4), (2, 5), (1,2)])
            sage: L.fontsize()
            10
            sage: L.fontsize(20)
            sage: L.fontsize()
            20

        All the numbers on the axes will be very large in this plot::

            sage: L
            Graphics object consisting of 1 graphics primitive
        """
        if s is None:
            try:
                return self._fontsize
            except AttributeError:
                self._fontsize = 10
                return self._fontsize
        self._fontsize = int(s)

    def axes_labels_size(self, s=None):
        """
        Set the relative size of axes labels w.r.t. the axes tick marks.

        INPUT:

        - ``s`` -- float, relative size of axes labels w.r.t. to the tick marks,
          the size of the tick marks being set by :meth:`fontsize`

        If called with no input, return the current relative size.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: p = plot(sin(x^2), (x, -3, 3), axes_labels=['$x$','$y$'])
            sage: p.axes_labels_size()  # default value
            1.6
            sage: p.axes_labels_size(2.5)
            sage: p.axes_labels_size()
            2.5

        Now the axes labels are large w.r.t. the tick marks::

            sage: p                                                                     # needs sage.symbolic
            Graphics object consisting of 1 graphics primitive
        """
        if s is None:
            try:
                return self._axes_labels_size
            except AttributeError:
                self._axes_labels_size = 1.6
                return self._axes_labels_size
        self._axes_labels_size = float(s)

    def axes(self, show=None):
        """
        Set whether or not the `x` and `y` axes are shown
        by default.

        INPUT:

        - ``show`` -- boolean

        If called with no input, return the current axes setting.

        EXAMPLES::

            sage: L = line([(1,2), (3,-4), (2, 5), (1,2)])

        By default the axes are displayed.

        ::

            sage: L.axes()
            True

        But we turn them off, and verify that they are off

        ::

            sage: L.axes(False)
            sage: L.axes()
            False

        Displaying L now shows a triangle but no axes.

        ::

            sage: L
            Graphics object consisting of 1 graphics primitive
        """
        if show is None:
            try:
                return self._show_axes
            except AttributeError:
                self._show_axes = True
                return self._show_axes
        self._show_axes = bool(show)

    def axes_color(self, c=None):
        """
        Set the axes color.

        If called with no input, return the current axes_color setting.

        INPUT:

        - ``c`` -- an RGB color 3-tuple, where each tuple entry
          is a float between 0 and 1

        EXAMPLES: We create a line, which has like everything a default
        axes color of black.

        ::

            sage: L = line([(1,2), (3,-4), (2, 5), (1,2)])
            sage: L.axes_color()
            (0, 0, 0)

        We change the axes color to red and verify the change.

        ::

            sage: L.axes_color((1,0,0))
            sage: L.axes_color()
            (1.0, 0.0, 0.0)

        When we display the plot, we'll see a blue triangle and bright red
        axes.

        ::

            sage: L
            Graphics object consisting of 1 graphics primitive
        """
        if c is None:
            try:
                return self._axes_color

            except AttributeError:
                self._axes_color = (0.0, 0.0, 0.0)
                return self._axes_color
        self._axes_color = rgbcolor(c)

    def axes_labels(self, l=None):
        """
        Set the axes labels.

        INPUT:

        - ``l`` -- (default: ``None``) a list of two strings or
          ``None``

        OUTPUT: a 2-tuple of strings

        If l is ``None``, returns the current ``axes_labels``,
        which is itself by default ``None``. The default labels are both
        empty.

        EXAMPLES: We create a plot and put x and y axes labels on it.

        ::

            sage: p = plot(sin(x), (x, 0, 10))                                          # needs sage.symbolic
            sage: p.axes_labels(['$x$','$y$'])                                          # needs sage.symbolic
            sage: p.axes_labels()                                                       # needs sage.symbolic
            ('$x$', '$y$')

        Now when you plot p, you see x and y axes labels::

            sage: p                                                                     # needs sage.symbolic
            Graphics object consisting of 1 graphics primitive

        Notice that some may prefer axes labels which are not
        typeset::

            sage: plot(sin(x), (x, 0, 10), axes_labels=['x','y'])                       # needs sage.symbolic
            Graphics object consisting of 1 graphics primitive

        TESTS:

        Unicode strings are acceptable; see :issue:`13161`. Note that
        this does not guarantee that matplotlib will handle the strings
        properly, although it should.

        ::

            sage: c = circle((0,0), 1)
            sage: c.axes_labels(['axe des abscisses', 'axe des ordonnées'])
            sage: c._axes_labels
            ('axe des abscisses', 'axe des ordonnées')
        """
        if l is None:
            try:
                return self._axes_labels
            except AttributeError:
                self._axes_labels = None
                return self._axes_labels
        if not isinstance(l, (list, tuple)):
            raise TypeError("l must be a list or tuple")
        if len(l) != 2:
            raise ValueError("l must have length 2")
        self._axes_labels = tuple(l)

    def axes_label_color(self, c=None):
        r"""
        Set the color of the axes labels.

        The axes labels are placed at the edge of the x and y axes, and are
        not on by default (use the :meth:`axes_labels` command to
        set them; see the example below). This function just changes their
        color.

        INPUT:

        - ``c`` -- an RGB 3-tuple of numbers between 0 and 1


        If called with no input, return the current axes_label_color
        setting.

        EXAMPLES: We create a plot, which by default has axes label color
        black.

        ::

            sage: p = plot(sin, (-1,1))                                                 # needs sage.symbolic
            sage: p.axes_label_color()                                                  # needs sage.symbolic
            (0, 0, 0)

        We change the labels to be red, and confirm this::

            sage: p.axes_label_color((1,0,0))                                           # needs sage.symbolic
            sage: p.axes_label_color()                                                  # needs sage.symbolic
            (1.0, 0.0, 0.0)

        We set labels, since otherwise we won't see anything.

        ::

            sage: p.axes_labels(['$x$ axis', '$y$ axis'])                               # needs sage.symbolic

        In the plot below, notice that the labels are red::

            sage: p                                                                     # needs sage.symbolic
            Graphics object consisting of 1 graphics primitive
        """
        if c is None:
            try:
                return self._axes_label_color
            except AttributeError:
                self._axes_label_color = (0, 0, 0)
                return self._axes_label_color
        self._axes_label_color = rgbcolor(c)

    def axes_width(self, w=None):
        r"""
        Set the axes width. Use this to draw a plot with really fat or
        really thin axes.

        INPUT:

        - ``w`` -- a float


        If called with no input, return the current
        ``axes_width`` setting.

        EXAMPLES: We create a plot, see the default axes width (with funny
        Python float rounding), then reset the width to 10 (very fat).

        ::

            sage: # needs sage.symbolic
            sage: p = plot(cos, (-3,3))
            sage: p.axes_width()
            0.8
            sage: p.axes_width(10)
            sage: p.axes_width()
            10.0

        Finally we plot the result, which is a graph with very fat axes.

        ::

            sage: p                                                                     # needs sage.symbolic
            Graphics object consisting of 1 graphics primitive
        """
        if w is None:
            try:
                return self._axes_width
            except AttributeError:
                self._axes_width = True
                return self._axes_width
        self._axes_width = float(w)

    def tick_label_color(self, c=None):
        """
        Set the color of the axes tick labels.

        INPUT:

        - ``c`` -- an RGB 3-tuple of numbers between 0 and 1


        If called with no input, return the current ``tick_label_color``
        setting.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: p = plot(cos, (-3,3))
            sage: p.tick_label_color()
            (0, 0, 0)
            sage: p.tick_label_color((1,0,0))
            sage: p.tick_label_color()
            (1.0, 0.0, 0.0)
            sage: p
            Graphics object consisting of 1 graphics primitive
        """
        if c is None:
            try:
                return self._tick_label_color
            except AttributeError:
                self._tick_label_color = (0, 0, 0)
                return self._tick_label_color
        self._tick_label_color = rgbcolor(c)

    def _repr_(self):
        r"""
        Return a string representation of the graphics objects.

        OUTPUT: string

        EXAMPLES:

        We create a plot and call :meth:`show` on it, which causes it
        to be displayed as a plot::

            sage: P = plot(cos, (-1,1))                                                 # needs sage.symbolic
            sage: P.show()                                                              # needs sage.symbolic

        Just doing this also displays the plot::

            sage: P                                                                     # needs sage.symbolic
            Graphics object consisting of 1 graphics primitive

        Using the Python `repr` or `str` commands do not display the
        plot::

            sage: repr(P)                                                               # needs sage.symbolic
            'Graphics object consisting of 1 graphics primitive'
            sage: str(P)                                                                # needs sage.symbolic
            'Graphics object consisting of 1 graphics primitive'
            sage: print(P)                                                              # needs sage.symbolic
            Graphics object consisting of 1 graphics primitive

        TESTS::

            sage: P._repr_()                                                            # needs sage.symbolic
            'Graphics object consisting of 1 graphics primitive'
        """
        return str(self)

    def _rich_repr_(self, display_manager, **kwds):
        """
        Rich Output Magic Method.

        See :mod:`sage.repl.rich_output` for details.

        EXAMPLES::

            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: g = Graphics()
            sage: g._rich_repr_(dm)
            OutputImagePng container
        """
        types = display_manager.types
        prefer_raster = (
            ('.png', types.OutputImagePng),
            ('.jpg', types.OutputImageJpg),
            ('.gif', types.OutputImageGif),
        )
        prefer_vector = (
            ('.svg', types.OutputImageSvg),
            ('.pdf', types.OutputImagePdf),
        )
        graphics = display_manager.preferences.graphics
        if graphics == 'disable':
            return
        elif graphics == 'raster' or graphics is None:
            preferred = prefer_raster + prefer_vector
        elif graphics == 'vector':
            preferred = prefer_vector + prefer_raster
        else:
            raise ValueError('unknown graphics output preference')
        for file_ext, output_container in preferred:
            if output_container in display_manager.supported_output():
                return display_manager.graphics_from_save(
                    self.save, kwds, file_ext, output_container)

    def __str__(self):
        r"""
        Return string representation of this plot.

        OUTPUT: string

        EXAMPLES::

            sage: S = circle((0,0), 2); S.__str__()
            'Graphics object consisting of 1 graphics primitive'
            sage: str(S)
            'Graphics object consisting of 1 graphics primitive'
            sage: print(S)
            Graphics object consisting of 1 graphics primitive
        """
        s = "Graphics object consisting of %s graphics primitives" % (len(self))
        if len(self) == 1:
            s = s[:-1]
        return s

    def __getitem__(self, i):
        """
        Return the i-th graphics primitive object:

        EXAMPLES::

            sage: G = circle((1,1),2) + circle((2,2),5); print(G)
            Graphics object consisting of 2 graphics primitives
            sage: G[1]
            Circle defined by (2.0,2.0) with r=5.0
        """
        return self._objects[i]

    def __len__(self):
        """
        If G is of type Graphics, then len(G) gives the number of distinct
        graphics primitives making up that object.

        EXAMPLES::

            sage: G = circle((1,1),1) + circle((1,2),1) + circle((1,2),5); print(G)
            Graphics object consisting of 3 graphics primitives
            sage: len(G)
            3
        """
        return len(self._objects)

    def __delitem__(self, i):
        """
        If G is of type Graphics, then del(G[i]) removes the i-th distinct
        graphic primitive making up that object.

        EXAMPLES::

            sage: G = circle((1,1),1) + circle((1,2),1) + circle((1,2),5); print(G)
            Graphics object consisting of 3 graphics primitives
            sage: len(G)
            3
            sage: del(G[2])
            sage: print(G)
            Graphics object consisting of 2 graphics primitives
            sage: len(G)
            2
        """
        del self._objects[int(i)]

    def __setitem__(self, i, x):
        """
        You can replace a :class:`GraphicPrimitive` (point, line, circle, etc...) in
        a :class:`Graphics` object G with any other :class:`GraphicPrimitive`

        EXAMPLES::

            sage: G = circle((1,1),1) + circle((1,2),1) + circle((1,2),5); print(G)
            Graphics object consisting of 3 graphics primitives

        ::

            sage: p = polygon([[1,3],[2,-2],[1,1],[1,3]]); print(p)
            Graphics object consisting of 1 graphics primitive

        ::

            sage: G[1] = p[0]
            sage: G    # show the plot
            Graphics object consisting of 3 graphics primitives
        """
        from sage.plot.primitive import GraphicPrimitive
        if not isinstance(x, GraphicPrimitive):
            raise TypeError("x must be a GraphicPrimitive")
        self._objects[int(i)] = x

    def __radd__(self, other):
        """
        Compute and return other + this graphics object.

        This only works when other is a Python int equal to 0. In all other
        cases a :exc:`TypeError` is raised. The main reason for this
        function is to make summing a list of graphics objects easier.

        EXAMPLES::

            sage: S = circle((0,0), 2)
            sage: print(int(0) + S)
            Graphics object consisting of 1 graphics primitive
            sage: print(S + int(0))
            Graphics object consisting of 1 graphics primitive

        The following would fail were it not for this function::

            sage: v = [circle((0,0), 2), circle((2,3), 1)]
            sage: print(sum(v))
            Graphics object consisting of 2 graphics primitives
        """
        if isinstance(other, int) and other == 0:
            return self
        raise TypeError

    def __add__(self, other):
        """
        If you have any Graphics object G1, you can always add any other
        amount of Graphics objects G2,G3,... to form a new Graphics object:
        ``G4 = G1 + G2 + G3``.

        The xmin, xmax, ymin, and ymax properties of the graphics objects
        are expanded to include all objects in both scenes. If the aspect
        ratio property of either or both objects are set, then the larger
        aspect ratio is chosen, with 'automatic' being overridden by a
        numeric aspect ratio.

        If one of the graphics object is set to show a legend, then
        the resulting object will also be set to show a legend. Legend
        options are propagated if set. If the same legend option is
        present in both arguments, the latter value is used.

        EXAMPLES::

            sage: g1 = plot(abs(sqrt(x^3-1)), (x,1,5), frame=True)                      # needs sage.symbolic
            sage: g2 = plot(-abs(sqrt(x^3-1)), (x,1,5), color='red')                    # needs sage.symbolic
            sage: g1 + g2  # displays the plot                                          # needs sage.symbolic
            Graphics object consisting of 2 graphics primitives

        TESTS:

        Extra keywords to show are propagated::

            sage: # needs sage.symbolic
            sage: (g1 + g2)._extra_kwds=={'aspect_ratio': 'automatic', 'frame': True}
            True
            sage: g1.set_aspect_ratio(2)
            sage: (g1+g2).aspect_ratio()
            2.0
            sage: g2.set_aspect_ratio(3)
            sage: (g1+g2).aspect_ratio()
            3.0

        As are legend options, :issue:`12936`::

            sage: # needs sage.symbolic
            sage: p1 = plot(x, x, 0, 1)
            sage: p2 = p1
            sage: p1.set_legend_options(back_color='black')
            sage: p2.set_legend_options(shadow=False)
            sage: p3 = p1 + p2
            sage: p3._legend_opts
            {'back_color': 'black', 'shadow': False}

        If the same legend option is specified more than once, the
        latter takes precedence::

            sage: # needs sage.symbolic
            sage: p1 = plot(x, x, 0, 1)
            sage: p2 = p1
            sage: p1.set_legend_options(shadow=True)
            sage: p2.set_legend_options(shadow=False)
            sage: p3 = p1 + p2
            sage: p3._legend_opts
            {'shadow': False}

        Flipped axes take precedence over non-flipped axes::

            sage: p1 = plot(x, x, 0, 1, flip_x=True, flip_y=True)                       # needs sage.symbolic
            sage: p2 = plot(x^2, x, 0, 1)                                               # needs sage.symbolic
            sage: [p._extra_kwds[k] for p in [p1 + p2, p2 + p1]                         # needs sage.symbolic
            ....:  for k in ['flip_x', 'flip_y']]
            [True, True, True, True]
        """
        if isinstance(other, int) and other == 0:
            return self
        if not isinstance(other, Graphics):
            from sage.plot.plot3d.base import Graphics3d
            if isinstance(other, Graphics3d):
                return self.plot3d() + other
            raise TypeError("other (=%s) must be a Graphics objects" % other)
        g = Graphics()
        g._objects = self._objects + other._objects
        g._show_legend = self._show_legend or other._show_legend
        g._extra_kwds.update(self._extra_kwds)
        g._extra_kwds.update(other._extra_kwds)
        g._legend_colors = self._legend_colors + other._legend_colors
        g._legend_opts.update(self._legend_opts)
        g._legend_opts.update(other._legend_opts)
        if 'flip_x' in self._extra_kwds and 'flip_x' in other._extra_kwds:
            g._extra_kwds['flip_x'] = (self._extra_kwds['flip_x']
                                       or other._extra_kwds['flip_x'])
        if 'flip_y' in self._extra_kwds and 'flip_y' in other._extra_kwds:
            g._extra_kwds['flip_y'] = (self._extra_kwds['flip_y']
                                       or other._extra_kwds['flip_y'])
        if self.aspect_ratio() == 'automatic':
            g.set_aspect_ratio(other.aspect_ratio())
        elif other.aspect_ratio() == 'automatic':
            g.set_aspect_ratio(self.aspect_ratio())
        else:
            g.set_aspect_ratio(max(self.aspect_ratio(), other.aspect_ratio()))
        return g

    def add_primitive(self, primitive):
        """
        Add a primitive to this graphics object.

        EXAMPLES:

        We give a very explicit example::

            sage: G = Graphics()
            sage: from math import e
            sage: from sage.plot.line import Line
            sage: from sage.plot.arrow import Arrow
            sage: L = Line([3,4,2,7,-2], [1,2,e,4,5.],
            ....:          {'alpha': 1, 'thickness': 2, 'rgbcolor': (0,1,1),
            ....:           'legend_label': ''})
            sage: A = Arrow(2, -5, .1, .2,
            ....:           {'width': 3, 'head': 0, 'rgbcolor': (1,0,0),
            ....:            'linestyle': 'dashed', 'zorder': 8, 'legend_label': ''})
            sage: G.add_primitive(L)
            sage: G.add_primitive(A)
            sage: G
            Graphics object consisting of 2 graphics primitives
        """
        self._objects.append(primitive)

    def plot(self):
        """
        Draw a 2D plot of this graphics object, which just returns this
        object since this is already a 2D graphics object.

        EXAMPLES::

            sage: S = circle((0,0), 2)
            sage: S.plot() is S
            True

        It does not accept any argument (:issue:`19539`)::

            sage: S.plot(1)
            Traceback (most recent call last):
            ...
            TypeError: ...plot() takes 1 positional argument but 2 were given

            sage: S.plot(hey='hou')
            Traceback (most recent call last):
            ...
            TypeError: ...plot() got an unexpected keyword argument 'hey'
        """
        return self

    def plot3d(self, z=0, **kwds):
        """
        Return an embedding of this 2D plot into the xy-plane of 3D space,
        as a 3D plot object. An optional parameter ``z`` can be given to
        specify the z-coordinate.

        EXAMPLES::

            sage: sum(plot(z*sin(x), 0, 10).plot3d(z)   # long time                     # needs sage.symbolic
            ....:     for z in range(6))
            Graphics3d Object
        """
        from sage.plot.plot3d.base import Graphics3dGroup
        g = Graphics3dGroup([g.plot3d(**kwds) for g in self._objects])
        if z:
            g = g.translate(0, 0, z)
        return g

    @classmethod
    def _extract_kwds_for_show(cls, kwds, ignore=[]):
        """
        Extract keywords relevant to show() from the provided dictionary.

        EXAMPLES::

            sage: kwds = {'f': lambda x: x, 'xmin': 0, 'figsize': [1,1], 'plot_points': (40, 40)}
            sage: G_kwds = Graphics._extract_kwds_for_show(kwds, ignore='xmin')
            sage: kwds  # Note how this action modifies the passed dictionary
            {'f': <function <lambda> at 0x...>,
             'plot_points': (40, 40),
             'xmin': 0}
            sage: G_kwds
            {'figsize': [1, 1]}

        This method is intended to be used with _set_extra_kwds(). Here is an
        idiom to ensure the correct keywords will get passed on to show()::

            sage: options = {}  # Usually this will come from an argument
            sage: g = Graphics()
            sage: g._set_extra_kwds(Graphics._extract_kwds_for_show(options))
        """
        result = {}
        for option in cls.SHOW_OPTIONS:
            if option not in ignore:
                try:
                    result[option] = kwds.pop(option)
                except KeyError:
                    pass
        return result

    def _set_extra_kwds(self, kwds):
        """
        Set a dictionary of keywords that will get passed on to show().

        TESTS::

            sage: g = Graphics()
            sage: g._extra_kwds
            {}
            sage: g._set_extra_kwds({'figsize': [10,10]})
            sage: g._extra_kwds
            {'figsize': [10, 10]}
            sage: g.show()  # Now the (blank) plot will be extra large
        """
        self._extra_kwds = kwds

    def _set_scale(self, subplot, scale=None, base=None):
        """
        Set the scale of the axes in the current subplot. This function is
        only for internal use.

        INPUT:

        - ``subplot`` -- matplotlib Axes instance
        - ``scale`` -- the scale of the figure. Values it can take are
          ``'linear'``, ``'loglog'``, ``'semilogx'``, ``'semilogy'``. See
          :meth:`show` for other options it can take.
        - ``base`` -- the base of the logarithm if a logarithmic scale is
          set. See :meth:`show` for the options it can take

        OUTPUT:

        The scale in the form of a tuple: (xscale, yscale, basex, basey)

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: p = plot(x, 1, 10)
            sage: fig = p.matplotlib()
            sage: ax = fig.get_axes()[0]
            sage: p._set_scale(ax, scale='linear', base=2)
            ('linear', 'linear', 10, 10)
            sage: p._set_scale(ax, scale='semilogy', base=2)
            ('linear', 'log', 10, 2)
            sage: p._set_scale(ax, scale=('loglog', 2, 3))
            ('log', 'log', 2, 3)
            sage: p._set_scale(ax, scale=['semilogx', 2])
            ('log', 'linear', 2, 10)

        TESTS::

            sage: # needs sage.symbolic
            sage: p._set_scale(ax, 'log')
            Traceback (most recent call last):
            ...
            ValueError: The scale must be one of 'linear', 'loglog', 'semilogx' or 'semilogy' -- got 'log'
            sage: p._set_scale(ax, ('loglog', 1))
            Traceback (most recent call last):
            ...
            ValueError: The base of the logarithm must be greater than 1
        """
        if scale is None:
            return ('linear', 'linear', 10, 10)
        if isinstance(scale, (list, tuple)):
            if len(scale) != 2 and len(scale) != 3:
                raise ValueError("If the input is a tuple, it must be of "
                                 "the form (scale, base) or (scale, basex, basey)")
            if len(scale) == 2:
                base = scale[1]
            else:
                base = scale[1:]
            scale = scale[0]

        if scale not in ('linear', 'loglog', 'semilogx', 'semilogy'):
            raise ValueError("The scale must be one of 'linear', 'loglog',"
                             f" 'semilogx' or 'semilogy' -- got '{scale}'")

        if isinstance(base, (list, tuple)):
            basex, basey = base
        elif base is None:
            basex = basey = 10
        else:
            basex = basey = base

        if basex <= 1 or basey <= 1:
            raise ValueError("The base of the logarithm must be greater "
                             "than 1")

        xscale = yscale = 'linear'
        if scale == 'linear':
            basex = basey = 10
        elif scale == 'loglog':
            subplot.set_xscale('log', base=basex)
            subplot.set_yscale('log', base=basey)
            xscale = yscale = 'log'
        elif scale == 'semilogx':
            subplot.set_xscale('log', base=basex)
            basey = 10
            xscale = 'log'
        elif scale == 'semilogy':
            subplot.set_yscale('log', base=basey)
            basex = 10
            yscale = 'log'

        return (xscale, yscale, basex, basey)

    # This dictionary has the default values for the keywords to show(). When
    # show is invoked with keyword arguments, those arguments are merged with
    # this dictionary to create a set of keywords with the defaults filled in.
    # Then, those keywords are passed on to save().

    # NOTE: If you intend to use a new parameter in show(), you should update
    # this dictionary to contain the default value for that parameter.

    SHOW_OPTIONS = {  # axes options
        'axes': None, 'axes_labels': None, 'axes_labels_size': None,
        'axes_pad': None, 'base': None, 'scale': None,
        'xmin': None, 'xmax': None, 'ymin': None, 'ymax': None,
        'flip_x': False, 'flip_y': False,
        # Figure options
        'aspect_ratio': None, 'dpi': DEFAULT_DPI, 'fig_tight': True,
        'figsize': None, 'fontsize': None, 'frame': False,
        'title': None, 'title_pos': None, 'transparent': False,
        # Grid options
        'gridlines': None, 'gridlinesstyle': None,
        'hgridlinesstyle': None, 'vgridlinesstyle': None,
        # Legend options
        'legend_options': {}, 'show_legend': None,
        # Ticks options
        'ticks': None, 'tick_formatter': None, 'ticks_integer': False,
        # Text options
        'typeset': 'default'}

    # Default options for the legends:

    LEGEND_OPTIONS = {'back_color': 'white', 'borderpad': 0.6,
                      'borderaxespad': None,
                      'columnspacing': None,
                      'fancybox': False, 'font_family': 'sans-serif',
                      'font_size': 'medium', 'font_style': 'normal',
                      'font_variant': 'normal', 'font_weight': 'medium',
                      'handlelength': 0.05, 'handletextpad': 0.5,
                      'labelspacing': 0.02, 'loc': 'best',
                      'markerscale': 0.6, 'ncol': 1, 'numpoints': 2,
                      'shadow': True, 'title': None}

    @suboptions('legend', **LEGEND_OPTIONS)
    def show(self, **kwds):
        r"""
        Show this graphics image immediately.

        This method attempts to display the graphics immediately,
        without waiting for the currently running code (if any) to
        return to the command line. Be careful, calling it from within
        a loop will potentially launch a large number of external
        viewer programs.

        OPTIONAL INPUT:

        - ``dpi`` -- (default: 100) dots per inch

        - ``figsize`` -- (default: [6.4, 4.8]) [width, height] inches. The
          maximum value of each of the width and the height can be 327
          inches, at the default ``dpi`` of 100 dpi, which is just shy of
          the maximum allowed value of 32768 dots (pixels).

        - ``fig_tight`` -- boolean (default: ``True``); whether to clip the drawing
          tightly around drawn objects.  If ``True``, then the resulting
          image will usually not have dimensions corresponding to
          ``figsize``.  If ``False``, the resulting image will have
          dimensions corresponding to ``figsize``.

        - ``aspect_ratio`` -- the perceived height divided by the
          perceived width. For example, if the aspect ratio is set to ``1``, circles
          will look round and a unit square will appear to have sides
          of equal length, and if the aspect ratio is set ``2``, vertical units will be
          twice as long as horizontal units, so a unit square will be twice as
          high as it is wide.  If set to ``'automatic'``, the aspect ratio
          is determined by ``figsize`` and the picture fills the figure.

        - ``axes`` -- (default: ``True``)

        - ``axes_labels`` -- (default: ``None``) list (or tuple) of two
          strings; the first is used as the label for the horizontal
          axis, and the second for the vertical axis.

        - ``axes_labels_size`` -- (default: current setting -- 1.6) scale factor
          relating the size of the axes labels with respect to the size of the
          tick marks.

        - ``fontsize`` -- (default: current setting -- 10) positive
          integer; used for axes labels; if you make this very large,
          you may have to increase figsize to see all labels.

        - ``frame`` -- boolean (default: ``False``); draw a frame around the image

        - ``gridlines`` -- (default: ``None``) can be any of the following:

          - None, False: do not add grid lines.

          - True, "automatic", "major": add grid lines at major ticks of the axes.

          - "minor": add grid at major and minor ticks.

          - [xlist,ylist]: a tuple or list containing
            two elements, where xlist (or ylist) can be
            any of the following.


            - None, False: don't add horizontal (or vertical) lines.

            - True, "automatic", "major": add horizontal (or vertical) grid lines at
              the major ticks of the axes.

            - "minor": add horizontal (or vertical) grid lines at major and minor ticks of
              axes.

            - an iterable yielding numbers n or pairs (n,opts), where n
              is the coordinate of the line and opt is a dictionary of
              MATPLOTLIB options for rendering the line.


        - ``gridlinesstyle, hgridlinesstyle, vgridlinesstyle`` -
          (default: ``None``); a dictionary of MATPLOTLIB options for the
          rendering of the grid lines, the horizontal grid lines or the
          vertical grid lines, respectively.

        - ``transparent`` -- boolean (default: ``False``); if True, make the
          background transparent

        - ``axes_pad`` -- (default: 0.02 on ``'linear'`` scale, 1 on
          ``'log'`` scale)

          - In the ``'linear'`` scale, it determines the percentage of the
            axis range that is added to each end of each axis. This helps
            avoid problems like clipping lines because of line-width, etc.
            To get axes that are exactly the specified limits, set
            ``axes_pad`` to zero.

          - On the ``'log'`` scale, it determines the exponent of the
            fraction of the minimum (resp. maximum) that is subtracted from
            the minimum (resp. added to the maximum) value of the axis. For
            instance if the minimum is `m` and the base of the axis is `b`
            then the new minimum after padding the axis will be
            `m - m/b^{\mathrm{axes\_pad}}`.

        - ``ticks_integer`` -- boolean (default: ``False``); guarantee that the ticks
          are integers (the ``ticks`` option, if specified, will
          override this)

        - ``ticks`` -- a matplotlib locator for the major ticks, or
          a number. There are several options.  For more information about
          locators, type ``from matplotlib import ticker`` and then
          ``ticker?``.

          - If this is a locator object, then it is the locator for
            the horizontal axis.  A value of None means use the default
            locator.

          - If it is a list of two locators, then the first is for the
            horizontal axis and one for the vertical axis.  A value of
            None means use the default locator (so a value of
            [None, my_locator] uses my_locator for the vertical axis and
            the default for the horizontal axis).

          - If in either case above one of the entries is a number `m`
            (something which can be coerced to a float), it will be
            replaced by a MultipleLocator which places major ticks at
            integer multiples of `m`.  See examples.

          - If in either case above one of the entries is a list of
            numbers, it will be replaced by a FixedLocator which places
            ticks at the locations specified.  This includes the case of
            of the empty list, which will give no ticks.  See examples.

        - ``tick_formatter`` -- a matplotlib formatter for the major
          ticks. There are several options.  For more information about
          formatters, type ``from matplotlib import ticker`` and then
          ``ticker?``.

          If the value of this keyword is a single item, then this will
          give the formatting for the horizontal axis *only* (except for
          the ``'latex'`` option).  If it is a list or tuple, the first
          is for the horizontal axis, the second for the vertical axis.
          The options are below:

          - If one of the entries is a formatter object, then it used.
            A value of None means to use the default locator (so using
            ``tick_formatter=[None, my_formatter]`` uses my_formatter
            for the vertical axis and the default for the horizontal axis).

          - If one of the entries is a symbolic constant such as `\pi`,
            `e`, or `sqrt(2)`, ticks will be formatted nicely at rational
            multiples of this constant.

          .. warning::

             This should only be used with the ``ticks`` option using nice
             rational multiples of that constant!

          - If one of the entries is the string ``'latex'``, then the
            formatting will be nice typesetting of the ticks.  This is
            intended to be used when the tick locator for at least one of
            the axes is a list including some symbolic elements. This uses
            matplotlib's internal LaTeX rendering engine. If you want to
            use an external LaTeX compiler, then set the keyword option
            ``typeset``.  See examples.

        - ``title`` -- (default: ``None``) the title for the plot

        - ``title_pos`` -- (default: ``None``) the position of the title for the
            plot. It must be a tuple or a list of two real numbers
            ``(x_pos, y_pos)`` which indicate the relative position of the
            title within the plot. The plot itself can be considered to
            occupy, in relative terms, the region within a unit square
            `[0, 1] \times [0, 1]`.  The title text is centered around the
            horizontal factor ``x_pos`` of the plot. The baseline of the
            title text is present at the vertical factor ``y_pos`` of the
            plot. Hence, ``title_pos=(0.5, 0.5)`` will center the title in
            the plot, whereas ``title_pos=(0.5, 1.1)`` will center the
            title along the horizontal direction, but will place the title
            a fraction `0.1` times above the plot.

          - If the first entry is a list of strings (or numbers), then the
            formatting for the horizontal axis will be typeset with the strings
            present in the list. Each entry of the list of strings must be
            provided with a corresponding number in the first entry of
            ``ticks`` to indicate its position on the axis. To typeset the
            strings with ``'latex'`` enclose them within ``'$'`` symbols. To
            have similar custom formatting of the labels along the vertical
            axis, the second entry must be a list of strings and the second
            entry of ``ticks`` must also be a list of numbers which give the
            positions of the labels. See the examples below.

        - ``show_legend`` -- (default: ``None``) if ``True``, show the legend

        - ``legend_*`` -- all the options valid for :meth:`set_legend_options`
            prefixed with ``legend_``

        - ``base`` -- (default: 10) the base of the logarithm if
          a logarithmic scale is set. This must be greater than 1. The base
          can be also given as a list or tuple ``(basex, basey)``.
          ``basex`` sets the base of the logarithm along the horizontal
          axis and ``basey`` sets the base along the vertical axis.

        - ``scale`` -- (default: ``'linear'``) string. The scale of the axes.
          Possible values are

          - ``'linear'`` -- linear scaling of both the axes
          - ``'loglog'`` -- sets both the horizontal and vertical axes to
            logarithmic scale
          - ``'semilogx'`` -- sets only the horizontal axis to logarithmic
            scale
          - ``'semilogy'`` -- sets only the vertical axis to logarithmic
            scale

          The scale can be also be given as single argument that is a list
          or tuple ``(scale, base)`` or ``(scale, basex, basey)``.

          .. NOTE::

            - If the ``scale`` is ``'linear'``, then irrespective of what
              ``base`` is set to, it will default to 10 and will remain
              unused.

        - ``xmin`` -- starting x value in the rendered figure

        - ``xmax`` -- ending x value in the rendered figure

        - ``ymin`` -- starting y value in the rendered figure

        - ``ymax`` -- ending y value in the rendered figure

        - ``flip_x`` -- boolean (default: ``False``); if ``True``, flip the horizontal
          axis

        - ``flip_y`` -- boolean (default: ``False``); if ``True``, flip the vertical
          axis

        - ``typeset`` -- (default: ``'default'``) string. The type of
          font rendering that should be used for the text. The possible
          values are

          - ``'default'`` -- uses matplotlib's internal text rendering
            engine called Mathtext ( see
            https://matplotlib.org/users/mathtext.html ). If you have
            modified the default matplotlib settings, for instance via
            a matplotlibrc file, then this option will not change any of
            those settings.
          - ``'latex'`` -- laTeX is used for rendering the fonts. This
            requires LaTeX, dvipng and Ghostscript to be installed
          - ``'type1'`` -- type 1 fonts are used by matplotlib in the text
            in the figure.  This requires LaTeX, dvipng and Ghostscript to
            be installed.

        OUTPUT:

        This method does not return anything. Use :meth:`save` if you
        want to save the figure as an image.

        EXAMPLES::

            sage: c = circle((1,1), 1, color='red')
            sage: c.show(xmin=-1, xmax=3, ymin=-1, ymax=3)

        You can make the picture larger by changing ``figsize`` with width,
        height each having a maximum value of 327 inches at default dpi::

            sage: p = ellipse((0,0),4,1)
            sage: p.show(figsize=[327,10], dpi=100)
            sage: p.show(figsize=[328,10], dpi=80)

        You can turn off the drawing of the axes::

            sage: show(plot(sin,-4,4), axes=False)                                      # needs sage.symbolic

        You can also label the axes.  Putting something in dollar
        signs formats it as a mathematical expression::

            sage: show(plot(sin,-4,4), axes_labels=('$x$','$y$'))                       # needs sage.symbolic

        You can add a title to a plot::

            sage: show(plot(sin,-4,4), title=r'A plot of $\sin(x)$')                    # needs sage.symbolic

        You can also provide the position for the title to the plot. In the
        plot below the title is placed on the bottom left of the figure.::

            sage: plot(sin, -4, 4, title='Plot sin(x)', title_pos=(0.05,-0.05))         # needs sage.symbolic
            Graphics object consisting of 1 graphics primitive

        If you want all the text to be rendered by using an external LaTeX
        installation then set the ``typeset`` to ``'latex'``. This
        requires that LaTeX, dvipng and Ghostscript be installed::

            sage: plot(x, typeset='latex')                                      # optional - latex, needs sage.symbolic
            Graphics object consisting of 1 graphics primitive

        If you want all the text in your plot to use Type 1 fonts, then
        set the ``typeset`` option to ``'type1'``. This requires that
        LaTeX, dvipng and Ghostscript be installed::

            sage: plot(x, typeset='type1')                                      # optional - latex, needs sage.symbolic
            Graphics object consisting of 1 graphics primitive

        You can turn on the drawing of a frame around the plots::

            sage: show(plot(sin,-4,4), frame=True)                                      # needs sage.symbolic

        You can make the background transparent::

            sage: plot(sin(x), (x, -4, 4), transparent=True)                            # needs sage.symbolic
            Graphics object consisting of 1 graphics primitive

        Prior to :issue:`19485`, legends by default had a shadowless gray
        background. This behavior can be recovered by passing in certain
        ``legend_options``::

            sage: p = plot(sin(x), legend_label=r'$\sin(x)$')                           # needs sage.symbolic
            sage: p.show(legend_options={'back_color': (0.9,0.9,0.9),                   # needs sage.symbolic
            ....:                        'shadow': False})

        We can change the scale of the axes in the graphics before
        displaying::

            sage: G = plot(exp, 1, 10)                                                  # needs sage.symbolic
            sage: G.show(scale='semilogy')                                              # needs sage.symbolic

        We can change the base of the logarithm too. The following changes
        the vertical axis to be on log scale, and with base 2. Note that
        the ``base`` argument will ignore any changes to the axis which is
        in linear scale.::

            sage: G.show(scale='semilogy', base=2)  # y axis as powers of 2         # long time, needs sage.symbolic

        ::

            sage: G.show(scale='semilogy', base=(3,2))  # base ignored for x-axis       # needs sage.symbolic

        The scale can be also given as a 2-tuple or a 3-tuple.::

            sage: G.show(scale=('loglog', 2.1))   # both x and y axes in base 2.1   # long time, needs sage.symbolic

        ::

            sage: G.show(scale=('loglog', 2, 3))  # x in base 2, y in base 3        # long time, needs sage.symbolic

        The base need not be an integer, though it does have to be made
        a float.::

            sage: G.show(scale='semilogx', base=float(e))  # base is e                  # needs sage.symbolic

        Logarithmic scale can be used for various kinds of plots. Here are
        some examples.::

            sage: G = list_plot([10**i for i in range(10)])                         # long time, needs sage.symbolic
            sage: G.show(scale='semilogy')                                          # long time, needs sage.symbolic

        ::

            sage: G = parametric_plot((x, x**2), (x, 1, 10))                            # needs sage.symbolic
            sage: G.show(scale='loglog')                                                # needs sage.symbolic

        ::

            sage: disk((5,5), 4, (0, 3*pi/2)).show(scale='loglog',base=2)               # needs sage.symbolic

        ::

            sage: x, y = var('x, y')                                                    # needs sage.symbolic
            sage: G = plot_vector_field((2^x,y^2), (x,1,10), (y,1,100))                 # needs sage.symbolic
            sage: G.show(scale='semilogx',base=2)                                       # needs sage.symbolic

        Flip the horizontal or vertical axis.

        ::

            sage: G = plot(x^3, -2, 3)                                                  # needs sage.symbolic
            sage: G.show(flip_x=True)                                                   # needs sage.symbolic
            sage: G.show(flip_y=True)                                                   # needs sage.symbolic

        Add grid lines at the major ticks of the axes.

        ::

            sage: c = circle((0,0), 1)
            sage: c.show(gridlines=True)
            sage: c.show(gridlines='automatic')
            sage: c.show(gridlines='major')

        Add grid lines at the major and minor ticks of the axes.

        ::

            sage: # needs sage.symbolic
            sage: u,v = var('u v')
            sage: f = exp(-(u^2+v^2))
            sage: p = plot_vector_field(f.gradient(), (u,-2,2), (v,-2,2))
            sage: p.show(gridlines='minor')

        Add only horizontal or vertical grid lines.

        ::

            sage: p = plot(sin, -10, 20)                                                # needs sage.symbolic
            sage: p.show(gridlines=[None, "automatic"])                                 # needs sage.symbolic
            sage: p.show(gridlines=["minor", False])                                    # needs sage.symbolic

        Add grid lines at specific positions (using lists/tuples).

        ::

            sage: x, y = var('x, y')                                                    # needs sage.symbolic
            sage: p = implicit_plot((y^2-x^2)*(x-1)*(2*x-3) - 4*(x^2+y^2-2*x)^2,        # needs sage.symbolic
            ....:                   (x,-2,2), (y,-2,2), plot_points=1000)
            sage: p.show(gridlines=[[1,0],[-1,0,1]])                                    # needs sage.symbolic

        Add grid lines at specific positions (using iterators).

        ::

            sage: def maple_leaf(t):
            ....:     return (100/(100+(t-pi/2)^8))*(2-sin(7*t)-cos(30*t)/2)
            sage: p = polar_plot(maple_leaf, -pi/4, 3*pi/2,                         # long time, needs sage.symbolic
            ....:                color='red',plot_points=1000)
            sage: p.show(gridlines=([-3,-2.75,..,3], range(-1,5,2)))                # long time, needs sage.symbolic

        Add grid lines at specific positions (using functions).

        ::

            sage: # needs sage.symbolic
            sage: y = x^5 + 4*x^4 - 10*x^3 - 40*x^2 + 9*x + 36
            sage: p = plot(y, -4.1, 1.1)
            sage: xlines = lambda a, b: [z for z, m in y.roots()]
            sage: p.show(gridlines=[xlines, [0]], frame=True, axes=False)

        Change the style of all the grid lines.

        ::

            sage: b = bar_chart([-3,5,-6,11], color='red')
            sage: b.show(gridlines=([-1,-0.5,..,4], True),
            ....:        gridlinesstyle=dict(color='blue', linestyle=':'))

        Change the style of the horizontal or vertical grid lines
        separately.

        ::

            sage: p = polar_plot(2 + 2*cos(x), 0, 2*pi, color=hue(0.3))                 # needs sage.symbolic
            sage: p.show(gridlines=True,                                                # needs sage.symbolic
            ....:        hgridlinesstyle=dict(color='orange', linewidth=1.0),
            ....:        vgridlinesstyle=dict(color='blue', linestyle=':'))

        Change the style of each grid line individually.

        ::

            sage: x, y = var('x, y')                                                    # needs sage.symbolic
            sage: p = implicit_plot((y^2-x^2)*(x-1)*(2*x-3) - 4*(x^2+y^2-2*x)^2,        # needs sage.symbolic
            ....:                   (x,-2,2), (y,-2,2), plot_points=1000)
            sage: p.show(gridlines=(                                                    # needs sage.symbolic
            ....:    [
            ....:     (1,{"color":"red","linestyle":":"}),
            ....:     (0,{"color":"blue","linestyle":"--"})
            ....:    ],
            ....:    [
            ....:     (-1,{"color":"red","linestyle":":"}),
            ....:     (0,{"color":"blue","linestyle":"--"}),
            ....:     (1,{"color":"red","linestyle":":"}),
            ....:    ]
            ....:    ),
            ....:    gridlinesstyle=dict(marker='x',color='black'))

        Grid lines can be added to contour plots.

        ::

            sage: f = sin(x^2 + y^2)*cos(x)*sin(y)                                      # needs sage.symbolic
            sage: c = contour_plot(f, (x, -4, 4), (y, -4, 4), plot_points=100)          # needs sage.symbolic
            sage: c.show(gridlines=True,                                                # needs sage.symbolic
            ....:        gridlinesstyle={'linestyle': ':', 'linewidth': 1, 'color': 'red'})

        Grid lines can be added to matrix plots.

        ::

            sage: M = MatrixSpace(QQ,10).random_element()
            sage: matrix_plot(M).show(gridlines=True)

        By default, Sage increases the horizontal and vertical axes
        limits by a certain percentage in all directions.  This is
        controlled by the ``axes_pad`` parameter.  Increasing the range
        of the axes helps avoid problems with lines and dots being
        clipped because the linewidth extends beyond the axes.  To get
        axes limits that are exactly what is specified, set
        ``axes_pad`` to zero.  Compare the following two examples

        ::

            sage: (plot(sin(x), (x, -pi, pi), thickness=2)                              # needs sage.symbolic
            ....:   + point((pi, -1), pointsize=15))
            Graphics object consisting of 2 graphics primitives
            sage: (plot(sin(x), (x, -pi, pi), thickness=2, axes_pad=0)                  # needs sage.symbolic
            ....:   + point((pi, -1), pointsize=15))
            Graphics object consisting of 2 graphics primitives

        The behavior of the ``axes_pad`` parameter is different if the axis
        is in the ``'log'`` scale. If `b` is the base of the axis, the
        minimum value of the axis, is decreased by the factor
        `1/b^{\mathrm{axes\_pad}}` of the minimum and the maximum value of the axis
        is increased by the same factor of the maximum value.  Compare the
        axes in the following two plots to see the difference.

        ::

            sage: plot_loglog(x, (1.1*10**-2, 9990))                                    # needs sage.symbolic
            Graphics object consisting of 1 graphics primitive

            sage: plot_loglog(x, (1.1*10**-2, 9990), axes_pad=0)                        # needs sage.symbolic
            Graphics object consisting of 1 graphics primitive

        Via matplotlib, Sage allows setting of custom ticks.  See above
        for more details.

        Here the labels are not so useful::

            sage: plot(sin(pi*x), (x, -8, 8))                                           # needs sage.symbolic
            Graphics object consisting of 1 graphics primitive

        Now put ticks at multiples of 2::

            sage: plot(sin(pi*x), (x, -8, 8), ticks=2)                                  # needs sage.symbolic
            Graphics object consisting of 1 graphics primitive

        Or just choose where you want the ticks::

            sage: plot(sin(pi*x), (x, -8, 8), ticks=[[-7,-3,0,3,7], [-1/2,0,1/2]])      # needs sage.symbolic
            Graphics object consisting of 1 graphics primitive

        Or no ticks at all::

            sage: plot(sin(pi*x), (x, -8, 8), ticks=[[], []])                           # needs sage.symbolic
            Graphics object consisting of 1 graphics primitive

        This can be very helpful in showing certain features of plots. ::

            sage: plot(1.5/(1+e^(-x)), (x, -10, 10))  # doesn't quite show value of inflection point                    # needs sage.symbolic
            Graphics object consisting of 1 graphics primitive

        ::

            sage: plot(1.5/(1+e^(-x)), (x, -10, 10),  # It's right at f(x)=0.75!        # needs sage.symbolic
            ....:      ticks=[None, 1.5/4])
            Graphics object consisting of 1 graphics primitive

        But be careful to leave enough room for at least two major ticks, so that
        the user can tell what the scale is::

            sage: plot(x^2, (x,1,8), ticks=6).show()                                    # needs sage.symbolic
            Traceback (most recent call last):
            ...
            ValueError: Expand the range of the independent variable to
            allow two multiples of your tick locator (option `ticks`).

        We can also do custom formatting if you need it.  See above for full
        details::

            sage: plot(2*x + 1, (x,0,5),        # not tested (broken with matplotlib 3.6), needs sage.symbolic
            ....:      ticks=[[0,1,e,pi,sqrt(20)], 2],
            ....:      tick_formatter='latex')
            Graphics object consisting of 1 graphics primitive

        This is particularly useful when setting custom ticks in multiples
        of `\pi`.

        ::

            sage: plot(sin(x), (x,0,2*pi), ticks=pi/3, tick_formatter=pi)               # needs sage.symbolic
            Graphics object consisting of 1 graphics primitive

        But keep in mind that you will get exactly the formatting you asked
        for if you specify both formatters.  The first syntax is recommended
        for best style in that case. ::

            sage: plot(arcsin(x), (x,-1,1), ticks=[None, pi/6],  # Nice-looking!        # needs sage.symbolic
            ....:      tick_formatter=["latex", pi])
            Graphics object consisting of 1 graphics primitive

        ::

            sage: plot(arcsin(x), (x,-1,1), ticks=[None, pi/6],  # Not so nice-looking  # needs sage.symbolic
            ....:      tick_formatter=[None, pi])
            Graphics object consisting of 1 graphics primitive

        Custom tick labels can be provided by providing the keyword
        ``tick_formatter`` with the list of labels, and simultaneously
        providing the keyword ``ticks`` with the positions of the labels. ::

            sage: plot(x, (x,0,3), ticks=[[1,2.5], [0.5,1,2]],                          # needs sage.symbolic
            ....:      tick_formatter=[["$x_1$","$x_2$"], ["$y_1$","$y_2$","$y_3$"]])
            Graphics object consisting of 1 graphics primitive

        The following sets the custom tick labels only along the horizontal
        axis. ::

            sage: plot(x**2, (x,0,2), ticks=[[1,2], None],                              # needs sage.symbolic
            ....:      tick_formatter=[["$x_1$","$x_2$"], None])
            Graphics object consisting of 1 graphics primitive

        If the number of tick labels do not match the number of positions of
        tick labels, then it results in an error.::

            sage: plot(x**2, (x,0,2), ticks=[[2], None],                                # needs sage.symbolic
            ....:      tick_formatter=[["$x_1$","$x_2$"], None]).show()
            Traceback (most recent call last):
            ...
            ValueError: If the first component of the list `tick_formatter` is a list
            then the first component of `ticks` must also be a list of equal length.

        When using logarithmic scale along the axis, make sure to have
        enough room for two ticks so that the user can tell what the scale
        is. This can be effected by increasing the range of the independent
        variable, or by changing the ``base``, or by providing enough tick
        locations by using the ``ticks`` parameter.

        By default, Sage will expand the variable range so that at least two
        ticks are included along the logarithmic axis. However, if you
        specify ``ticks`` manually, this safety measure can be defeated::

            sage: list_plot_loglog([(1,2),(2,3)], plotjoined=True, ticks=[[1],[1]])
            doctest:...: UserWarning: The x-axis contains fewer than 2 ticks;
            the logarithmic scale of the plot may not be apparent to the reader.
            doctest:...: UserWarning: The y-axis contains fewer than 2 ticks;
            the logarithmic scale of the plot may not be apparent to the reader.
            Graphics object consisting of 1 graphics primitive

        This one works, since the horizontal axis is automatically expanded
        to contain two ticks and the vertical axis is provided with two ticks::

            sage: list_plot_loglog([(1,2),(2,3)], plotjoined=True, ticks=[None,[1,10]])
            Graphics object consisting of 1 graphics primitive

        Another example in the log scale where both the axes are automatically
        expanded to show two major ticks::

            sage: list_plot_loglog([(2,0.5), (3, 4)], plotjoined=True)
            Graphics object consisting of 1 graphics primitive

        When using ``title_pos``, it must be ensured that a list or a tuple
        of length two is used. Otherwise, a warning is raised::

            sage: plot(x, -4, 4, title='Plot x', title_pos=0.05)                        # needs sage.symbolic
            doctest:...: ...RichReprWarning: Exception in _rich_repr_ while displaying
            object: 'title_pos' must be a list or tuple of two real numbers.
            Graphics object consisting of 1 graphics primitive

        TESTS:

        The following tests result in a segmentation fault and should not
        be run or doctested::

            sage: p = ellipse((0,0),4,1)
            sage: p.show(figsize=[232,232],dpi=100)  # not tested
            ------------------------------------------------------------------------
            Unhandled SIGSEGV: A segmentation fault occurred.
            This probably occurred because a *compiled* module has a bug
            in it and is not properly wrapped with sig_on(), sig_off().
            Python will now terminate.
            ------------------------------------------------------------------------
            sage: p.show(figsize=[327,181],dpi=100)  # not tested
            ------------------------------------------------------------------------
            Unhandled SIGSEGV: A segmentation fault occurred.
            This probably occurred because a *compiled* module has a bug
            in it and is not properly wrapped with sig_on(), sig_off().
            Python will now terminate.
            ------------------------------------------------------------------------

        The following tests ensure we give a good error message for
        negative figsizes::

            sage: # needs sage.symbolic
            sage: P = plot(x^2,(x,0,1))
            sage: P.show(figsize=[-1,1])
            Traceback (most recent call last):
            ...
            ValueError: figsize should be positive numbers, not -1.0 and 1.0
            sage: P.show(figsize=-1)
            Traceback (most recent call last):
            ...
            ValueError: figsize should be positive, not -1.0
            sage: P.show(figsize=x^2)
            Traceback (most recent call last):
            ...
            TypeError: figsize should be a positive number, not x^2
            sage: P.show(figsize=[2,3,4])
            Traceback (most recent call last):
            ...
            ValueError: figsize should be a positive number or
            a list of two positive numbers, not [2, 3, 4]
            sage: P.show(figsize=[sqrt(2),sqrt(3)])
        """
        from sage.repl.rich_output import get_display_manager
        dm = get_display_manager()
        dm.display_immediately(self, **kwds)

    def xmin(self, xmin=None):
        """
        EXAMPLES::

            sage: g = line([(-1,1), (3,2)])
            sage: g.xmin()
            -1.0
            sage: g.xmin(-3)
            sage: g.xmin()
            -3.0
        """
        if xmin is None:
            return self.get_axes_range()['xmin']
        else:
            self.set_axes_range(xmin=xmin)

    def xmax(self, xmax=None):
        """
        EXAMPLES::

            sage: g = line([(-1,1), (3,2)])
            sage: g.xmax()
            3.0
            sage: g.xmax(10)
            sage: g.xmax()
            10.0
        """
        if xmax is None:
            return self.get_axes_range()['xmax']
        else:
            self.set_axes_range(xmax=xmax)

    def ymin(self, ymin=None):
        """
        EXAMPLES::

            sage: g = line([(-1,1), (3,2)])
            sage: g.ymin()
            1.0
            sage: g.ymin(-3)
            sage: g.ymin()
            -3.0
        """
        if ymin is None:
            return self.get_axes_range()['ymin']
        else:
            self.set_axes_range(ymin=ymin)

    def ymax(self, ymax=None):
        """
        EXAMPLES::

            sage: g = line([(-1,1), (3,2)])
            sage: g.ymax()
            2.0
            sage: g.ymax(10)
            sage: g.ymax()
            10.0
        """
        if ymax is None:
            return self.get_axes_range()['ymax']
        else:
            self.set_axes_range(ymax=ymax)

    def get_minmax_data(self):
        r"""
        Return the x and y coordinate minimum and maximum.

        .. warning::

           The returned dictionary is mutable, but changing it does
           not change the xmin/xmax/ymin/ymax data.  The minmax data is a function
           of the primitives which make up this Graphics object.  To change the
           range of the axes, call methods :meth:`xmin`, :meth:`xmax`,
           :meth:`ymin`, :meth:`ymax`, or :meth:`set_axes_range`.

        OUTPUT:

        A dictionary whose keys give the xmin, xmax, ymin, and ymax
        data for this graphic.

        EXAMPLES::

            sage: g = line([(-1,1), (3,2)])
            sage: list(sorted(g.get_minmax_data().items()))
            [('xmax', 3.0), ('xmin', -1.0), ('ymax', 2.0), ('ymin', 1.0)]

        Note that changing ymax doesn't change the output of get_minmax_data::

            sage: g.ymax(10)
            sage: list(sorted(g.get_minmax_data().items()))
            [('xmax', 3.0), ('xmin', -1.0), ('ymax', 2.0), ('ymin', 1.0)]

        The width/height ratio (in output units, after factoring in the
        chosen aspect ratio) of the plot is limited to `10^{-15}\dots
        10^{15}`, otherwise floating point errors cause problems in
        matplotlib::

            sage: l = line([(1e-19,-1), (-1e-19,+1)], aspect_ratio=1.0)
            sage: l.get_minmax_data()
            {'xmax': 1.00010000000000e-15,
             'xmin': -9.99900000000000e-16,
             'ymax': 1.0,
             'ymin': -1.0}
            sage: l = line([(0,0), (1,1)], aspect_ratio=1e19)
            sage: l.get_minmax_data()
            {'xmax': 5000.50000000000, 'xmin': -4999.50000000000,
             'ymax': 1.0, 'ymin': 0.0}
        """
        objects = self._objects
        if objects:
            minmax_data = [o.get_minmax_data() for o in objects]
            xmin = min(d['xmin'] for d in minmax_data)
            xmax = max(d['xmax'] for d in minmax_data)
            ymin = min(d['ymin'] for d in minmax_data)
            ymax = max(d['ymax'] for d in minmax_data)
            if isnan(xmin):
                xmin = 0
                sage.misc.verbose.verbose("xmin was NaN (setting to 0)", level=0)
            if isnan(xmax):
                xmax = 0
                sage.misc.verbose.verbose("xmax was NaN (setting to 0)", level=0)
            if isnan(ymin):
                ymin = 0
                sage.misc.verbose.verbose("ymin was NaN (setting to 0)", level=0)
            if isnan(ymax):
                ymax = 0
                sage.misc.verbose.verbose("ymax was NaN (setting to 0)", level=0)
        else:
            xmin = xmax = ymin = ymax = 0

        if xmin == xmax:
            xmin -= 1
            xmax += 1
        if ymin == ymax:
            ymin -= 1
            ymax += 1
        return self._limit_output_aspect_ratio(xmin, xmax, ymin, ymax)

    def _limit_output_aspect_ratio(self, xmin, xmax, ymin, ymax):
        r"""
        Private helper function for :meth:`get_minmax_data`.

        INPUT:

        - ``xmin``, ``xmax``, ``ymin``, ``ymax`` -- bounding box for
          the graphics

        OUTPUT:

        A dictionary whose keys give the xmin, xmax, ymin, and ymax
        data for this graphic. Possibly enlarged in order to keep the
        width/height ratio (in output units, after factoring in the
        chosen aspect ratio) of the plot is limited to `10^{-15}\dots
        10^{15}` to avoid floating point issues in matplotlib.

        EXAMPLES::

            sage: l = line([(0,0), (1,1)], aspect_ratio=1.0)
            sage: l._limit_output_aspect_ratio(1, 2, 1e19, 3)
            {'xmax': -4999.50000000000,
             'xmin': 5000.50000000000,
             'ymax': 3,
             'ymin': 1.00000000000000e19}
            sage: l._limit_output_aspect_ratio(1, 2, 3, 1e19)
            {'xmax': 5000.50000000000,
             'xmin': -4999.50000000000,
             'ymax': 1.00000000000000e19,
             'ymin': 3}
            sage: l = line([(0,0), (1,1)], aspect_ratio=1e16)
            sage: l._limit_output_aspect_ratio(0, 1, 2, 3)
            {'xmax': 5.50000000000000, 'xmin': -4.50000000000000, 'ymax': 3, 'ymin': 2}
        """
        aspect_ratio = self.aspect_ratio()
        if aspect_ratio != 'automatic':
            width = xmax - xmin
            height = ymax - ymin
            output_aspect = abs(width / height / aspect_ratio)
            if output_aspect > 1e15:
                height = 1e15 * width / aspect_ratio
                ycenter = (ymax - ymin) / 2
                ymin = ycenter - height / 2
                ymax = ycenter + height / 2
            if output_aspect < 1e-15:
                width = 1e-15 * height * aspect_ratio
                xcenter = (xmax - xmin) / 2
                xmin = xcenter - width / 2
                xmax = xcenter + width / 2
        return {'xmin': xmin, 'xmax': xmax, 'ymin': ymin, 'ymax': ymax}

    def _matplotlib_tick_formatter(self, subplot, base=(10, 10),
                                   locator_options={}, scale=('linear', 'linear'),
                                   tick_formatter=(None, None), ticks=(None, None),
                                   xmax=None, xmin=None, ymax=None, ymin=None):
        r"""
        Take a matplotlib subplot instance representing the graphic and set
        the ticks formatting. This function is only for internal use.

        INPUT:

        - ``subplot`` -- the subplot instance

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: from matplotlib.figure import Figure
            sage: p = plot(x); d = p.get_minmax_data()
            sage: subplot = Figure().add_subplot(111)
            sage: p._objects[0]._render_on_subplot(subplot)
            sage: p._matplotlib_tick_formatter(subplot, **d)
            (<Axes...>,
            <matplotlib.ticker.MaxNLocator object at ...>,
            <matplotlib.ticker.MaxNLocator object at ...>,
            <matplotlib.ticker.ScalarFormatter object at ...>,
            <matplotlib.ticker.ScalarFormatter object at ...>)
        """
        # This function is created to refactor some code that is repeated
        # in the matplotlib function
        from matplotlib.ticker import (FixedLocator, Locator,
                                       LogFormatterMathtext,
                                       LogLocator, MaxNLocator,
                                       MultipleLocator,
                                       NullLocator, ScalarFormatter)

        x_locator, y_locator = ticks
        # ---------------------- Location of x-ticks ---------------------

        if x_locator is None:
            if scale[0] == 'log':
                x_locator = LogLocator(base=base[0])
            else:
                x_locator = MaxNLocator(**locator_options)
        elif isinstance(x_locator, Locator):
            pass
        elif x_locator == []:
            x_locator = NullLocator()
        elif isinstance(x_locator, list):
            x_locator = FixedLocator([float(x) for x in x_locator])
        else:  # x_locator is a number which can be made a float
            from sage.functions.other import ceil, floor
            if floor(xmax / x_locator) - ceil(xmin / x_locator) > 1:
                x_locator = MultipleLocator(float(x_locator))
            else:  # not enough room for two major ticks
                raise ValueError('Expand the range of the independent '
                                 'variable to allow two multiples of your tick locator '
                                 '(option `ticks`).')

        # ---------------------- Location of y-ticks ---------------------
        if y_locator is None:
            if scale[1] == 'log':
                y_locator = LogLocator(base=base[1])
            else:
                y_locator = MaxNLocator(**locator_options)
        elif isinstance(y_locator, Locator):
            pass
        elif y_locator == []:
            y_locator = NullLocator()
        elif isinstance(y_locator, list):
            y_locator = FixedLocator([float(y) for y in y_locator])
        else:  # y_locator is a number which can be made a float
            from sage.functions.other import ceil, floor
            if floor(ymax / y_locator) - ceil(ymin / y_locator) > 1:
                y_locator = MultipleLocator(float(y_locator))
            else:  # not enough room for two major ticks
                raise ValueError('Expand the range of the dependent '
                                 'variable to allow two multiples of your tick locator '
                                 '(option `ticks`).')

        x_formatter, y_formatter = tick_formatter
        from matplotlib.ticker import FuncFormatter, FixedFormatter
        from sage.misc.latex import latex
        from sage.structure.element import Expression
        from .misc import _multiple_of_constant
        # ---------------------- Formatting x-ticks ----------------------
        if x_formatter is None:
            if scale[0] == 'log':
                x_formatter = LogFormatterMathtext(base=base[0])
            else:
                x_formatter = ScalarFormatter()
        elif isinstance(x_formatter, Expression):
            x_const = x_formatter
            x_formatter = FuncFormatter(lambda n, pos:
                                        _multiple_of_constant(n, pos, x_const))
        elif x_formatter == "latex":
            if scale[0] == 'log':
                # We need to strip out '\\mathdefault' from the string
                x_formatter = FuncFormatter(lambda n, pos:
                                            LogFormatterMathtext(base=base[0])(n, pos).replace(
                                                "\\mathdefault", ""))
            else:
                # circumvent the problem of symbolic tick values (trac #34693)
                if isinstance(x_locator, FixedLocator):
                    x_formatter = FixedFormatter(['$%s$' % latex(n) for n in ticks[0]])
                else:
                    x_formatter = FuncFormatter(lambda n, pos: '$%s$' % latex(n))
        elif isinstance(x_formatter, (list, tuple)):
            if (not isinstance(ticks[0], (list, tuple)) or
                    len(ticks[0]) != len(x_formatter)):
                raise ValueError("If the first component of the list "
                                 "`tick_formatter` is a list then the first component "
                                 "of `ticks` must also be a list of equal length.")
            x_formatter = FixedFormatter(x_formatter)
        # ---------------------- Formatting y-ticks ----------------------
        if y_formatter is None:
            if scale[1] == 'log':
                y_formatter = LogFormatterMathtext(base=base[1])
            else:
                y_formatter = ScalarFormatter()
        elif isinstance(y_formatter, Expression):
            y_const = y_formatter
            y_formatter = FuncFormatter(lambda n, pos:
                                        _multiple_of_constant(n, pos, y_const))
        elif y_formatter == "latex":
            if scale[1] == 'log':
                # We need to strip out '\\mathdefault' from the string
                y_formatter = FuncFormatter(lambda n, pos:
                                            LogFormatterMathtext(base=base[1])(n, pos).replace(
                                                "\\mathdefault", ""))
            else:
                # circumvent the problem of symbolic tick values (trac #34693)
                if isinstance(y_locator, FixedLocator):
                    y_formatter = FixedFormatter(['$%s$' % latex(n) for n in ticks[1]])
                else:
                    y_formatter = FuncFormatter(lambda n, pos: '$%s$' % latex(n))
        elif isinstance(y_formatter, (list, tuple)):
            if (not isinstance(ticks[1], (list, tuple)) or
                    len(ticks[1]) != len(y_formatter)):
                raise ValueError("If the second component of the list "
                                 "`tick_formatter` is a list then the second component "
                                 "of `ticks` must also be a list of equal length.")
            y_formatter = FixedFormatter(y_formatter)

        subplot.xaxis.set_major_locator(x_locator)
        subplot.yaxis.set_major_locator(y_locator)
        subplot.xaxis.set_major_formatter(x_formatter)
        subplot.yaxis.set_major_formatter(y_formatter)

        # Check for whether there will be too few ticks in the log scale case.
        # If there are not enough ticks (2 or more) to determine that the scale
        # is non-linear, we throw a warning.
        from warnings import warn
        tickwarnmsg = 'The %s-axis contains fewer than 2 ticks; '
        tickwarnmsg += 'the logarithmic scale of the plot may not be apparent '
        tickwarnmsg += 'to the reader.'

        if (scale[0] == 'log' and not isinstance(x_locator, NullLocator) and
                len(subplot.xaxis.get_ticklocs()) < 2):
            warn(tickwarnmsg % 'x')

        if (scale[1] == 'log' and not isinstance(y_locator, NullLocator) and
                len(subplot.yaxis.get_ticklocs()) < 2):
            warn(tickwarnmsg % 'y')

        return (subplot, x_locator, y_locator, x_formatter, y_formatter)

    def _get_vmin_vmax(self, vmin, vmax, basev, axes_pad):
        r"""
        Determine the min/max value for a variable plotted on a logarithmic
        scale. The motivation is that we desire at least two ticks for a log
        plot; otherwise the reader may assume that the scale is linear. For
        internal use only.

        We check if this case occurs (for e.g. assuming xmin < xmax)::

           floor(logxmin)              ceil(logxmax)
           ----|---------+----------+----------|----------------------|--
                      logxmin     logxmax

        Or if this case occurs (assuming xmin < xmax)::

           floor(logxmin)             floor(logxmax)         ceil(logxmax)
           ----|---------+---------------------|-----+----------------|--
                      logxmin                     logxmax


        INPUT:

        - ``vmin`` -- the current min for this variable (e.g. xmin or ymin)

        - ``vmax`` -- the current max for this variable (e.g. xmax or ymax)

        - ``basev`` -- the base of the logarithmic scale for this variable

        - ``axes_pad`` -- the padding for the axis. It determines the
          exponent of the fraction of the minimum (resp. maximum) that is
          subtracted from the minimum (resp. added to the maximum) value of
          the axis. For instance if the minimum is `m` and the base of the
          axis is `b` then the new minimum after padding the axis will be
          `m - m/b^{\mathrm{axes\_pad}}`.

        OUTPUT:

        A new (min,max) pair for this variable, suitable for its logarithmic
        scale.

        EXAMPLES:

        On a base-10 logarithmic scale, we should have ``vmin``/``vmax``
        at least 10 units apart::

            sage: p = Graphics()
            sage: p._get_vmin_vmax(1, 2, 10, None) == (9/10, 10)
            True
            sage: p._get_vmin_vmax(1, 5, 10, None) == (9/10, 10)
            True
            sage: p._get_vmin_vmax(1, 10, 10, None)
            (9/10, 11)
            sage: p._get_vmin_vmax(1, 11, 10, None)
            (9/10, 121/10)
            sage: p._get_vmin_vmax(1, 50, 10, None)
            (9/10, 55)

        We can set the ``axes_pad`` separately::

            sage: p._get_vmin_vmax(1, 50, 2, 2)
            (0.75, 62.5)

        Nonpositive values of ``vmin`` are not accepted due to the domain
        of the logarithm function::

            sage: p = Graphics()
            sage: p._get_vmin_vmax(-1,2,10, None)
            Traceback (most recent call last):
            ...
            ValueError: vmin must be positive

        And ``vmax`` must be greater than ``vmin``::

            sage: p._get_vmin_vmax(1,-2,10, None)
            Traceback (most recent call last):
            ...
            ValueError: vmin must be less than vmax
        """
        if vmin <= 0:
            raise ValueError('vmin must be positive')

        if vmin >= vmax:
            raise ValueError('vmin must be less than vmax')

        import math
        if axes_pad is None:
            axes_pad = 1
        else:
            axes_pad = float(abs(axes_pad))

        logvmin = math.log(vmin) / math.log(basev)
        logvmax = math.log(vmax) / math.log(basev)

        if math.floor(logvmax) - math.ceil(logvmin) < 0:
            vmax = basev**math.ceil(logvmax)
            vmin = basev**math.floor(logvmin)
        elif math.floor(logvmax) - math.ceil(logvmin) < 1:
            if logvmax - math.floor(logvmax) > math.ceil(logvmin) - logvmin:
                vmax = basev**math.ceil(logvmax)
                if axes_pad > 0:
                    vmin -= vmin * basev**(-axes_pad)
            else:
                vmin = basev**math.floor(logvmin)
                if axes_pad > 0:
                    vmax += vmax * basev**(-axes_pad)
        elif axes_pad > 0:
            # pad the axes if we haven't expanded the axes earlier.
            vmin -= vmin * basev**(-axes_pad)
            vmax += vmax * basev**(-axes_pad)

        return vmin, vmax

    def matplotlib(self, filename=None,
                   xmin=None, xmax=None, ymin=None, ymax=None,
                   figsize=None, figure=None, sub=None,
                   axes=None, axes_labels=None, axes_labels_size=None,
                   flip_x=False, flip_y=False,
                   fontsize=None, frame=False, verify=True,
                   aspect_ratio=None,
                   gridlines=None, gridlinesstyle=None,
                   vgridlinesstyle=None, hgridlinesstyle=None,
                   show_legend=None, legend_options=None,
                   axes_pad=None, ticks_integer=None,
                   tick_formatter=None, ticks=None, title=None,
                   title_pos=None, base=None, scale=None,
                   stylesheet=None,
                   typeset='default'):
        r"""
        Construct or modify a Matplotlib figure by drawing ``self`` on it.

        INPUT (partial description, involving only Matplotlib objects; see
        :meth:`show` for the other arguments):

        - ``figure`` -- (default: ``None``) Matplotlib figure (class
          ``matplotlib.figure.Figure``) on which ``self`` is to be displayed;
          if ``None``, the figure will be created from the parameter
          ``figsize``

        - ``figsize`` -- (default: ``None``) width or [width, height] in inches
          of the Matplotlib figure in case ``figure`` is ``None``; if
          ``figsize`` is ``None``, Matplotlib's default (6.4 x 4.8 inches) is
          used

        - ``sub`` -- (default: ``None``) subpart of the figure, as an
          instance of Matplotlib "axes" (class ``matplotlib.axes.Axes``) on
          which ``self`` is to be drawn; if ``None``, the subpart will be
          created so as to cover the whole figure

        OUTPUT:

        - a ``matplotlib.figure.Figure`` object; if the argument ``figure`` is
          provided, this is the same object as ``figure``.

        EXAMPLES::

            sage: c = circle((1,1),1)
            sage: print(c.matplotlib())
            Figure(640x480)

        To obtain the first Matplotlib ``Axes`` object inside of the
        figure, you can do something like the following.

        ::

            sage: p = plot(sin(x), (x, -2*pi, 2*pi))                                    # needs sage.symbolic
            sage: figure = p.matplotlib()                                               # needs sage.symbolic
            sage: axes = figure.axes[0]                                                 # needs sage.symbolic

        Scientific notation is automatically used for large enough function values, cf. :issue:`34233` ::

            sage: plot(90000*x/(100-x), (x,0,100), ymin=0, ymax=1000000)                # needs sage.symbolic
            ...

        .. PLOT::

            from sage.symbolic.ring import var
            p = var('p')
            g = plot(90000*p/(100-p), (p,0,100), ymin=0, ymax=1000000)
            sphinx_plot(g)

        TESTS:

        We verify that :issue:`10291` is fixed::

            sage: # needs sage.symbolic
            sage: p = plot(sin(x), (x, -2*pi, 2*pi))
            sage: figure = p.matplotlib()
            sage: axes_range = p.get_axes_range()
            sage: figure = p.matplotlib()
            sage: axes_range2 = p.get_axes_range()
            sage: axes_range == axes_range2
            True

        We verify that legend options are properly handled (:issue:`12960`).
        First, we test with no options, and next with an incomplete set of
        options.::

            sage: # needs sage.symbolic
            sage: p = plot(x, legend_label='aha')
            sage: p.legend(True)
            sage: pm = p.matplotlib()
            sage: pm = p.matplotlib(legend_options={'font_size': 'small'})

        The title should not overlap with the axes labels nor the frame in
        the following plot (see :issue:`10512`)::

            sage: plot(sin(x^2), (x, -3, 3), title='Plot of sin(x^2)',                  # needs sage.symbolic
            ....:      axes_labels=['x','y'], frame=True)
            Graphics object consisting of 1 graphics primitive

        ``typeset`` must not be set to an arbitrary string::

            sage: plot(x, typeset='garbage')                                            # needs sage.symbolic
            doctest:...: ...RichReprWarning: Exception in _rich_repr_ while
            displaying object: typeset must be set to one of 'default',
            'latex', or 'type1'; got 'garbage'.
            Graphics object consisting of 1 graphics primitive

        We verify that numerical options are changed to float before saving (:issue:`14741`).
        By default, Sage 5.10 changes float objects to the `RealLiteral` type.
        The patch changes them to float before creating `matplotlib` objects.::

            sage: # long time, needs sage.symbolic
            sage: f = lambda x, y: abs(cos((x + I * y) ** 4)) - 1
            sage: g = implicit_plot(f, (-4, 4), (-3, 3), linewidth=0.6)
            sage: gm = g.matplotlib()

        If the axes are flipped, the limits of the axes get swapped::

            sage: # needs sage.symbolic
            sage: p = plot(2*x, 1, 2)
            sage: sub, = p.matplotlib(flip_y=True, flip_x=True).axes
            sage: xmin, xmax = sub.get_xlim()
            sage: ymin, ymax = sub.get_ylim()
            sage: xmin > xmax, ymin > ymax
            (...True..., ...True...)
        """
        if not isinstance(ticks, (list, tuple)):
            ticks = (ticks, None)
        if legend_options is None:
            legend_options = {}
        # as discussed in trac #25799 and #23696, Sage prefers the computer
        # modern fonts of TeX for math texts such as axes labels, but otherwise
        # adopts the default style of matplotlib
        from matplotlib import rcParams
        rcParams['mathtext.fontset'] = 'cm'
        rcParams['mathtext.rm'] = 'serif'
        rcParams['axes.formatter.use_mathtext'] = True

        import matplotlib.pyplot as plt
        if stylesheet in plt.style.available:
            plt.style.use(stylesheet)

        from sage.structure.element import Expression
        # make sure both formatters typeset or both don't
        if not isinstance(tick_formatter, (list, tuple)):
            if tick_formatter == "latex" or isinstance(tick_formatter, Expression):
                tick_formatter = (tick_formatter, "latex")
            else:
                tick_formatter = (tick_formatter, None)

        global do_verify
        do_verify = verify

        if axes is None:
            axes = self._show_axes

        from matplotlib.figure import Figure
        if typeset == 'type1':  # Requires LaTeX, dvipng, gs to be installed.
            rcParams['ps.useafm'] = True
            rcParams['pdf.use14corefonts'] = True
            rcParams['text.usetex'] = True
        elif typeset == 'latex':  # Requires LaTeX, dvipng, gs to be installed.
            rcParams['ps.useafm'] = False
            rcParams['pdf.use14corefonts'] = False
            rcParams['text.usetex'] = True
        elif typeset != 'default':  # We won't change (maybe user-set) defaults
            raise ValueError("typeset must be set to one of 'default', 'latex',"
                             f" or 'type1'; got '{typeset}'.")

        self.fontsize(fontsize)
        self.axes_labels(l=axes_labels)
        self.axes_labels_size(s=axes_labels_size)

        # If no matplotlib figure is provided, it is created here:
        if figure is None:
            if figsize is not None:
                figsize = _parse_figsize(figsize)
            figure = Figure(figsize=figsize)

        # The incoming subplot instance
        subplot = sub
        if not subplot:
            subplot = figure.add_subplot(111)
        # Add all the primitives to the subplot
        old_opts = {}
        for g in self._objects:
            opts, old_opts[g] = g.options(), g.options()
            for k, v in opts.items():
                try:
                    if v.parent() in sage.categories.fields.Fields():
                        opts[k] = float(v)
                except (AttributeError, TypeError):
                    pass
            g.set_options(opts)
            g._render_on_subplot(subplot)
            if hasattr(g, '_bbox_extra_artists'):
                self._bbox_extra_artists.extend(g._bbox_extra_artists)
        # Set the aspect ratio
        if aspect_ratio is None:
            aspect_ratio = self.aspect_ratio()
        if aspect_ratio == 'automatic':
            subplot.set_aspect('auto', adjustable='box')
        else:
            subplot.set_aspect(aspect_ratio, adjustable='box')

        # ---------------- Set the axes limits and scale ------------------
        self.set_axes_range(xmin, xmax, ymin, ymax)
        d = self.get_axes_range()
        xmin = d['xmax' if flip_x else 'xmin']
        xmax = d['xmin' if flip_x else 'xmax']
        ymin = d['ymax' if flip_y else 'ymin']
        ymax = d['ymin' if flip_y else 'ymax']

        xscale, yscale, basex, basey = self._set_scale(subplot, scale=scale,
                                                       base=base)

        # If any of the x-data are negative, we leave the min/max alone.
        if xscale == 'log' and min(xmin, xmax) > 0:
            if xmin < xmax:
                xmin, xmax = self._get_vmin_vmax(xmin, xmax, basex, axes_pad)
            else:
                xmax, xmin = self._get_vmin_vmax(xmax, xmin, basex, axes_pad)
        else:
            xpad = 0.02 if axes_pad is None else axes_pad
            xpad = (xmax - xmin) * float(xpad)
            xmax += xpad
            xmin -= xpad

        # Likewise for the y-data.
        if yscale == 'log' and min(ymin, ymax) > 0:
            if ymin < ymax:
                ymin, ymax = self._get_vmin_vmax(ymin, ymax, basey, axes_pad)
            else:
                ymax, ymin = self._get_vmin_vmax(ymax, ymin, basey, axes_pad)
        else:
            ypad = 0.02 if axes_pad is None else axes_pad
            ypad = (ymax - ymin) * float(ypad)
            ymax += ypad
            ymin -= ypad

        # -------------------------- Set the legend -----------------------
        if show_legend is None:
            show_legend = self._show_legend

        if show_legend:
            from matplotlib.font_manager import FontProperties
            lopts = {}
            lopts.update(legend_options)
            lopts.update(self._legend_opts)
            prop = FontProperties(
                family=lopts.pop('font_family', 'sans-serif'),
                size=lopts.pop('font_size', 'medium'),
                style=lopts.pop('font_style', 'normal'),
                weight=lopts.pop('font_weight', 'medium'),
                variant=lopts.pop('font_variant', 'normal'))
            color = lopts.pop('back_color', 'white')
            if 'loc' in lopts:
                loc = lopts['loc']
                if isinstance(loc, Integral):
                    # matplotlib 3.8 doesn't support sage integers
                    lopts['loc'] = int(loc)
            leg = subplot.legend(prop=prop, **lopts)
            if leg is None:
                from warnings import warn
                warn("legend requested but no items are labeled")
            else:
                # color
                lframe = leg.get_frame()
                lframe.set_facecolor(color)
                from sage.plot.colors import to_mpl_color
                for txt, color in zip(leg.get_texts(), self._legend_colors):
                    if color is not None:
                        txt.set_color(to_mpl_color(color))

        subplot.set_xlim([xmin, xmax])
        subplot.set_ylim([ymin, ymax])

        locator_options = {'nbins': 9, 'steps': [1, 2, 5, 10],
                           'integer': ticks_integer}

        if axes is None:
            axes = self._show_axes

        for spine in subplot.spines.values():
            spine.set_color(self._axes_color)
            spine.set_linewidth(self._axes_width)

        if frame:
            # For now, set the formatter to the old one, since that is
            # sort of what we are used to.  We should eventually look at
            # the default one to see if we like it better.

            (subplot, x_locator, y_locator,
                x_formatter, y_formatter) = self._matplotlib_tick_formatter(
                    subplot, base=(basex, basey),
                    locator_options=locator_options,
                    scale=(xscale, yscale),
                    tick_formatter=tick_formatter, ticks=ticks,
                    xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin)

            subplot.set_frame_on(True)
            if axes and xscale == 'linear' and yscale == 'linear':
                if (ymin <= 0 and ymax >= 0) or (ymax <= 0 and ymin >= 0):
                    subplot.axhline(color=self._axes_color,
                                    linewidth=self._axes_width)
                if (xmin <= 0 and xmax >= 0) or (xmax <= 0 and xmin >= 0):
                    subplot.axvline(color=self._axes_color,
                                    linewidth=self._axes_width)

        elif axes:
            ymiddle = False
            xmiddle = False
            # Note that the user may specify a custom xmin and xmax which
            # flips the axis horizontally. Hence we need to check for both
            # the possibilities in the if statements below. Similar
            # comments hold for ymin and ymax.
            if xscale == 'log':
                if xmax > xmin:
                    subplot.spines['right'].set_visible(False)
                    subplot.spines['left'].set_position(('outward', 10))
                    subplot.yaxis.set_ticks_position('left')
                    subplot.yaxis.set_label_position('left')
                    yaxis = 'left'
                elif xmax < xmin:
                    subplot.spines['left'].set_visible(False)
                    subplot.spines['right'].set_position(('outward', 10))
                    subplot.yaxis.set_ticks_position('right')
                    subplot.yaxis.set_label_position('right')
                    yaxis = 'right'
            elif (xmin > 0 and xmax > xmin) or (xmax > 0 and xmin > xmax):
                subplot.spines['right'].set_visible(False)
                subplot.spines['left'].set_position(('outward', 10))
                subplot.yaxis.set_ticks_position('left')
                subplot.yaxis.set_label_position('left')
                yaxis = 'left'
            elif (xmax < 0 and xmax > xmin) or (xmin < 0 and xmin > xmax):
                subplot.spines['left'].set_visible(False)
                subplot.spines['right'].set_position(('outward', 10))
                subplot.yaxis.set_ticks_position('right')
                subplot.yaxis.set_label_position('right')
                yaxis = 'right'
            else:
                subplot.spines['left'].set_position('zero')
                subplot.yaxis.set_ticks_position('left')
                subplot.yaxis.set_label_position('left')
                subplot.spines['right'].set_visible(False)
                ymiddle = True
                yaxis = 'left'

            if yscale == 'log':
                if ymax > ymin:
                    subplot.spines['top'].set_visible(False)
                    subplot.spines['bottom'].set_position(('outward', 10))
                    subplot.xaxis.set_ticks_position('bottom')
                    subplot.xaxis.set_label_position('bottom')
                    xaxis = 'bottom'
                elif ymax < ymin:
                    subplot.spines['bottom'].set_visible(False)
                    subplot.spines['top'].set_position(('outward', 10))
                    subplot.xaxis.set_ticks_position('top')
                    subplot.xaxis.set_label_position('top')
                    xaxis = 'top'
            elif (ymin > 0 and ymax > ymin) or (ymax > 0 and ymin > ymax):
                subplot.spines['top'].set_visible(False)
                subplot.spines['bottom'].set_position(('outward', 10))
                subplot.xaxis.set_ticks_position('bottom')
                subplot.xaxis.set_label_position('bottom')
                xaxis = 'bottom'
            elif (ymax < 0 and ymax > ymin) or (ymin < 0 and ymin > ymax):
                subplot.spines['bottom'].set_visible(False)
                subplot.spines['top'].set_position(('outward', 10))
                subplot.xaxis.set_ticks_position('top')
                subplot.xaxis.set_label_position('top')
                xaxis = 'top'
            else:
                subplot.spines['bottom'].set_position('zero')
                subplot.xaxis.set_ticks_position('bottom')
                subplot.xaxis.set_label_position('bottom')
                subplot.spines['top'].set_visible(False)
                xmiddle = True
                xaxis = 'bottom'

            # For now, set the formatter to the old one, since that is
            # sort of what we are used to.  We should eventually look at
            # the default one to see if we like it better.

            (subplot, x_locator, y_locator,
                x_formatter, y_formatter) = self._matplotlib_tick_formatter(
                    subplot, base=(basex, basey),
                    locator_options=locator_options,
                    scale=(xscale, yscale),
                    tick_formatter=tick_formatter, ticks=ticks,
                    xmax=xmax, xmin=xmin, ymax=ymax, ymin=ymin)

            # Make ticklines go on both sides of the axes
            #             if xmiddle:
            #                 for t in subplot.xaxis.get_majorticklines():
            #                     t.set_marker("|")
            #                     t.set_markersize(8)
            #                 for t in subplot.xaxis.get_minorticklines():
            #                     t.set_marker("|")
            #                     t.set_markersize(4)

            #             if ymiddle:
            #                 for t in subplot.yaxis.get_majorticklines():
            #                     t.set_marker("|")
            #                     t.set_markersize(8)
            #                 for t in subplot.yaxis.get_minorticklines():
            #                     t.set_marker("|")
            #                     t.set_markersize(4)

            # Make the zero tick labels disappear if the axes cross
            # inside the picture, but only if log scale is not used
            if (xmiddle and ymiddle and xscale == 'linear' == yscale):
                ## get the ticks, which are the numeric location of the ticks, and eliminate 0 from them
                if type(ticks[0]) is not list:
                    subplot.set_xticks([i for i in subplot.get_xticks() if i != 0][1:-1])
                if type(ticks[1]) is not list:
                    subplot.set_yticks([i for i in subplot.get_yticks() if i != 0][1:-1])
        else:
            for spine in subplot.spines.values():
                spine.set_visible(False)
            from matplotlib.ticker import NullFormatter, NullLocator
            subplot.xaxis.set_major_formatter(NullFormatter())
            subplot.yaxis.set_major_formatter(NullFormatter())
            subplot.xaxis.set_major_locator(NullLocator())
            subplot.yaxis.set_major_locator(NullLocator())

        if frame or axes:
            # Make minor tickmarks, unless we specify fixed ticks or no ticks
            # We do this change only on linear scale, otherwise matplotlib
            # errors out with a memory error.
            from matplotlib.ticker import (AutoMinorLocator, FixedLocator,
                                           LogLocator, NullLocator)
            if isinstance(x_locator, (NullLocator, FixedLocator)):
                subplot.xaxis.set_minor_locator(NullLocator())
            elif xscale == 'linear':
                subplot.xaxis.set_minor_locator(AutoMinorLocator())
            else:  # log scale
                from sage.arith.srange import srange
                base_inv = 1.0 / basex
                subs = [float(_) for _ in srange(2 * base_inv, 1, base_inv)]
                subplot.xaxis.set_minor_locator(LogLocator(base=basex,
                                                           subs=subs))
            if isinstance(y_locator, (NullLocator, FixedLocator)):
                subplot.yaxis.set_minor_locator(NullLocator())
            elif yscale == 'linear':
                subplot.yaxis.set_minor_locator(AutoMinorLocator())
            else:  # log scale
                from sage.arith.srange import srange
                base_inv = 1.0 / basey
                subs = [float(_) for _ in srange(2 * base_inv, 1, base_inv)]
                subplot.yaxis.set_minor_locator(LogLocator(base=basey,
                                                           subs=subs))
            # Set the color and fontsize of ticks
            subplot.tick_params(color=self._axes_color,
                                labelcolor=self._tick_label_color,
                                labelsize=self._fontsize, which='both')

        if gridlines is not None:
            if isinstance(gridlines, (list, tuple)):
                vgridlines, hgridlines = gridlines
            else:
                hgridlines = gridlines
                vgridlines = gridlines

            if gridlinesstyle is None:
                # Set up the default grid style
                gridlinesstyle = {'color': 'black', 'linestyle': ':',
                                  'linewidth': 0.5}

            vgridstyle = gridlinesstyle.copy()
            if vgridlinesstyle is not None:
                vgridstyle.update(vgridlinesstyle)

            hgridstyle = gridlinesstyle.copy()
            if hgridlinesstyle is not None:
                hgridstyle.update(hgridlinesstyle)

            if hgridlines == 'minor':
                hgridstyle['which'] = 'both'
            if vgridlines == 'minor':
                vgridstyle['which'] = 'both'

            if not isinstance(hgridlines, str) and isinstance(hgridlines, Iterable):
                hlines = iter(hgridlines)
                hgridstyle.pop("minor", None)
                for hline in hlines:
                    if isinstance(hline, (list, tuple)):
                        hl, style = hline
                        st = hgridstyle.copy()
                        st.update(style)
                    else:
                        hl = hline
                        st = hgridstyle
                    subplot.axhline(hl, **st)
            else:
                if hgridlines not in (None, False):
                    subplot.yaxis.grid(True, **hgridstyle)

            if not isinstance(vgridlines, str) and isinstance(vgridlines, Iterable):
                vlines = iter(vgridlines)
                vgridstyle.pop("minor", None)
                for vline in vlines:
                    if isinstance(vline, (list, tuple)):
                        vl, style = vline
                        st = vgridstyle.copy()
                        st.update(style)
                    else:
                        vl = vline
                        st = vgridstyle
                    subplot.axvline(vl, **st)
            else:
                if vgridlines not in (None, False):
                    subplot.xaxis.grid(True, **vgridstyle)

        if self._axes_labels is not None:
            label_options = {}
            label_options['color'] = self._axes_label_color
            label_options['size'] = int(self._axes_labels_size * self._fontsize)
            subplot.set_xlabel(self._axes_labels[0], **label_options)
            subplot.set_ylabel(self._axes_labels[1], **label_options)

            if axes is True and frame is False:
                # We set the label positions according to where we are
                # drawing the axes.
                if xaxis == 'bottom':
                    yaxis_labely = subplot.get_ylim()[1]
                    yaxis_labeloffset = 8
                    yaxis_vert = 'bottom'
                    xaxis_labely = 0
                    xaxis_vert = 'baseline'
                else:
                    yaxis_labely = subplot.get_ylim()[0]
                    yaxis_labeloffset = -8
                    yaxis_vert = 'top'
                    xaxis_labely = 1
                    xaxis_vert = 'top'

                if yaxis == 'left':
                    xaxis_labelx = subplot.get_xlim()[1]
                    xaxis_labeloffset = 8
                    xaxis_horiz = 'left'
                    yaxis_labelx = 0
                else:
                    xaxis_labelx = subplot.get_xlim()[0]
                    xaxis_labeloffset = -8
                    xaxis_horiz = 'right'
                    yaxis_labelx = 1

                from matplotlib.transforms import offset_copy
                xlabel = subplot.xaxis.get_label()
                xlabel.set_horizontalalignment(xaxis_horiz)
                xlabel.set_verticalalignment(xaxis_vert)
                trans = subplot.spines[xaxis].get_transform()
                labeltrans = offset_copy(trans, figure, x=xaxis_labeloffset,
                                         y=0, units='points')
                subplot.xaxis.set_label_coords(x=xaxis_labelx,
                                               y=xaxis_labely, transform=labeltrans)

                ylabel = subplot.yaxis.get_label()
                ylabel.set_horizontalalignment('center')
                ylabel.set_verticalalignment(yaxis_vert)
                ylabel.set_rotation('horizontal')
                trans = subplot.spines[yaxis].get_transform()
                labeltrans = offset_copy(trans, figure, x=0,
                                         y=yaxis_labeloffset, units='points')
                subplot.yaxis.set_label_coords(x=yaxis_labelx,
                                               y=yaxis_labely, transform=labeltrans)

        # This option makes the xlim and ylim limits not take effect
        # todo: figure out which limits were specified, and let the
        # free limits autoscale
        # subplot.autoscale_view(tight=True)
        if title is not None:
            if title_pos is not None:
                if (not isinstance(title_pos, (list, tuple)) or
                        len(title_pos) != 2):
                    raise ValueError("'title_pos' must be a list or tuple "
                                     "of two real numbers.")
                title_pos = (float(title_pos[0]), float(title_pos[1]))

            if (frame) or (axes_labels is None):
                if title_pos is not None:
                    subplot.set_title(title, fontsize=fontsize,
                                      position=title_pos)
                else:
                    subplot.set_title(title, fontsize=fontsize)
            else:
                # frame is false axes is not None, and neither is axes_labels
                # Then, the title is moved up to avoid overlap with axes labels
                if title_pos is None:
                    title_pos = (0.5, 1.05)
                subplot.set_title(title, fontsize=fontsize, position=title_pos)

        for g in self._objects:
            g.set_options(old_opts[g])

        return figure

    def save_image(self, filename=None, *args, **kwds):
        r"""
        Save an image representation of ``self``.

        The image type is determined by the extension of the filename.
        For example, this could be ``.png``, ``.jpg``, ``.gif``,
        ``.pdf``, ``.svg``.  Currently this is implemented by calling
        the :meth:`save` method of self, passing along all arguments
        and keywords.

        .. NOTE::

            Not all image types are necessarily implemented for all
            graphics types.  See :meth:`save` for more details.

        EXAMPLES::

            sage: import tempfile
            sage: c = circle((1,1), 1, color='red')
            sage: with tempfile.NamedTemporaryFile(suffix='.png') as f:
            ....:     c.save_image(f.name, xmin=-1, xmax=3,
            ....:                  ymin=-1, ymax=3)
        """
        self.save(filename, *args, **kwds)

    # filename argument is written explicitly so that it can be used as a
    # positional one, which is a very likely usage for this function.

    @suboptions('legend', **LEGEND_OPTIONS)
    def save(self, filename, **kwds):
        r"""
        Save the graphics to an image file.

        INPUT:

        - ``filename`` -- string. The filename and the image format
          given by the extension, which can be one of the following:

            * ``.eps``,

            * ``.pdf``,

            * ``.pgf``,

            * ``.png``,

            * ``.ps``,

            * ``.sobj`` (for a Sage object you can load later),

            * ``.svg``,

            * empty extension will be treated as ``.sobj``.

        All other keyword arguments will be passed to the plotter.

        OUTPUT: none

        EXAMPLES::

            sage: c = circle((1,1), 1, color='red')
            sage: from tempfile import NamedTemporaryFile
            sage: with NamedTemporaryFile(suffix='.png') as f:
            ....:     c.save(f.name, xmin=-1, xmax=3, ymin=-1, ymax=3)

        To make a figure bigger or smaller, use ``figsize``::

            sage: c.save(f.name, figsize=5, xmin=-1, xmax=3, ymin=-1, ymax=3)

        By default, the figure grows to include all of the graphics and text,
        so the final image may not be exactly the figure size you specified.
        If you want a figure to be exactly a certain size, specify the keyword
        ``fig_tight=False``::

            sage: c.save(f.name, figsize=[8,4], fig_tight=False,
            ....:        xmin=-1, xmax=3, ymin=-1, ymax=3)

        You can also pass extra options to the plot command instead of this
        method, e.g. ::

            sage: plot(x^2 - 5, (x, 0, 5), ymin=0).save(tmp_filename(ext='.png'))       # needs sage.symbolic

        will save the same plot as the one shown by this command::

            sage: plot(x^2 - 5, (x, 0, 5), ymin=0)                                      # needs sage.symbolic
            Graphics object consisting of 1 graphics primitive

        (This test verifies that :issue:`8632` is fixed.)

        TESTS:

        Legend labels should save correctly::

            sage: # needs sage.symbolic
            sage: P = plot(x,(x,0,1),legend_label='$xyz$')
            sage: P.set_legend_options(back_color=(1,0,0))
            sage: P.set_legend_options(loc=7)
            sage: import tempfile
            sage: with tempfile.NamedTemporaryFile(suffix='.png') as f:
            ....:     P.save(f.name)

        This plot should save with the frame shown, showing :issue:`7524`
        is fixed (same issue as :issue:`7981` and :issue:`8632`)::

            sage: var('x,y')                                                            # needs sage.symbolic
            (x, y)
            sage: a = plot_vector_field((x,-y),(x,-1,1),(y,-1,1))                       # needs sage.symbolic
            sage: import tempfile
            sage: with tempfile.NamedTemporaryFile(suffix='.png') as f:                 # needs sage.symbolic
            ....:     a.save(f.name)

        The following plot should show the axes; fixes :issue:`14782` ::

            sage: plot(x^2, (x, 1, 2), ticks=[[], []])                                  # needs sage.symbolic
            Graphics object consisting of 1 graphics primitive
        """
        options = {}
        options.update(self.SHOW_OPTIONS)
        options.update(self._extra_kwds)
        options.update(kwds)
        dpi = options.pop('dpi')
        transparent = options.pop('transparent')
        fig_tight = options.pop('fig_tight')

        ext = os.path.splitext(filename)[1].lower()

        if ext in ['', '.sobj']:
            SageObject.save(self, filename)
        elif ext not in ALLOWED_EXTENSIONS:
            raise ValueError("allowed file extensions for images are '" +
                             "', '".join(ALLOWED_EXTENSIONS) + "'!")
        else:
            from matplotlib import rcParams
            rc_backup = (rcParams['ps.useafm'], rcParams['pdf.use14corefonts'],
                         rcParams['text.usetex'])  # save the rcParams
            figure = self.matplotlib(**options)
            # You can output in PNG, PS, EPS, PDF, PGF, or SVG format, depending
            # on the file extension.
            # PGF is handled by a different backend
            if ext == '.pgf':
                from sage.features.latex import xelatex, pdflatex, lualatex
                latex_implementations = []
                if xelatex().is_present():
                    latex_implementations.append('xelatex')
                if pdflatex().is_present():
                    latex_implementations.append('pdflatex')
                if lualatex().is_present():
                    latex_implementations.append('lualatex')
                if not latex_implementations:
                    raise ValueError("Matplotlib requires either xelatex, "
                                     "lualatex, or pdflatex.")
                if latex_implementations[0] == "pdflatex":
                    # use pdflatex and set font encoding as per
                    # matplotlib documentation:
                    # https://matplotlib.org/users/pgf.html#pgf-tutorial
                    pgf_options = {"pgf.texsystem": "pdflatex",
                                   "pgf.preamble": [
                                       r"\usepackage[utf8x]{inputenc}",
                                       r"\usepackage[T1]{fontenc}"]}
                else:
                    pgf_options = {
                        "pgf.texsystem": latex_implementations[0],
                    }
                from matplotlib import rcParams
                rcParams.update(pgf_options)
                from matplotlib.backends.backend_pgf import FigureCanvasPgf
                figure.set_canvas(FigureCanvasPgf(figure))

            # matplotlib looks at the file extension to see what the renderer should be.
            # The default is FigureCanvasAgg for PNG's because this is by far the most
            # common type of files rendered, like in the notebook, for example.
            # if the file extension is not '.png', then matplotlib will handle it.
            else:
                from matplotlib.backends.backend_agg import FigureCanvasAgg
                figure.set_canvas(FigureCanvasAgg(figure))
            # this messes up the aspect ratio!
            # figure.canvas.mpl_connect('draw_event', pad_for_tick_labels)

            # tight_layout adjusts the *subplot* parameters so ticks aren't cut off, etc.
            figure.tight_layout()

            opts = {'dpi': dpi, 'transparent': transparent}
            if fig_tight is True:
                opts['bbox_inches'] = 'tight'
            if self._bbox_extra_artists:
                opts['bbox_extra_artists'] = self._bbox_extra_artists

            figure.savefig(filename, **opts)

            # Restore the rcParams to the original, possibly user-set values
            (rcParams['ps.useafm'], rcParams['pdf.use14corefonts'],
             rcParams['text.usetex']) = rc_backup

    def _latex_(self, **kwds):
        """
        Return a string plotting ``self`` with PGF.

        INPUT:

        - ``**kwds`` -- all keyword arguments will be passed to the plotter

        OUTPUT: string of PGF commands to plot ``self``

        EXAMPLES::

            sage: L = line([(0,0), (1,1)], axes=False)
            sage: L._latex_()     # not tested
            '%% Creator: Matplotlib, PGF backend...
        """
        tmpfilename = tmp_filename(ext='.pgf')
        self.save(filename=tmpfilename, **kwds)
        with open(tmpfilename) as tmpfile:
            latex_list = tmpfile.readlines()
        from sage.misc.latex import latex
        latex.add_package_to_preamble_if_available('pgf')
        return ''.join(latex_list)

    def description(self):
        r"""
        Print a textual description to stdout.

        This method is mostly used for doctests.

        EXAMPLES::

            sage: print(polytopes.hypercube(2).plot().description())                    # needs sage.geometry.polyhedron
            Polygon defined by 4 points: [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)]
            Line defined by 2 points: [(-1.0, 1.0), (-1.0, -1.0)]
            Line defined by 2 points: [(1.0, -1.0), (-1.0, -1.0)]
            Line defined by 2 points: [(1.0, -1.0), (1.0, 1.0)]
            Line defined by 2 points: [(1.0, 1.0), (-1.0, 1.0)]
            Point set defined by 4 point(s): [(1.0, -1.0), (1.0, 1.0), (-1.0, 1.0), (-1.0, -1.0)]
        """
        data = []
        for g in self:
            g_zorder = g.options().get('zorder', 0)
            if hasattr(g, 'xdata'):
                g_str = f'{g}:\t{list(zip(g.xdata, g.ydata))}'
            else:
                g_str = repr(g)
            data.append([g_zorder, g_str, g])
        data.sort()
        return '\n'.join(g[1] for g in data)

    def inset(self, graphics, pos=None, fontsize=None):
        r"""
        Add a graphics object as an inset.

        INPUT:

        - ``graphics`` -- the graphics object (instance of :class:`Graphics`)
          to be added as an inset to the current graphics

        - ``pos`` -- (default: ``None``) 4-tuple
          ``(left, bottom, width, height)``
          specifying the location and size of the inset on the final figure,
          all quantities being in fractions of the figure width and height; if
          ``None``, the value ``(0.7, 0.7, 0.2, 0.2)`` is used

        - ``fontsize`` -- (default: ``None``)  integer, font size (in points)
          for the inset; if ``None``, the value of 6 points is used, unless
          ``fontsize`` has been explicitly set in the construction of
          ``graphics`` (in this case, it is not overwritten here)

        OUTPUT: instance of :class:`~sage.plot.multigraphics.MultiGraphics`

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: f(x) = x^2*sin(1/x)
            sage: g1 = plot(f(x), (x, -2, 2), axes_labels=['$x$', '$y$'])
            sage: g2 = plot(f(x), (x, -0.3, 0.3), axes_labels=['$x$', '$y$'],
            ....:           frame=True)
            sage: g1.inset(g2)
            Multigraphics with 2 elements

        .. PLOT::

            f = (x**2*sin(1/x)).function(x)
            g1 = plot(f(x), (x, -2, 2), axes_labels=['$x$', '$y$'])
            g2 = plot(f(x), (x, -0.3, 0.3), axes_labels=['$x$', '$y$'], \
                      frame=True)
            sphinx_plot(g1.inset(g2))

        Using non-default values for the position/size and the font size::

            sage: g1.inset(g2, pos=(0.15, 0.7, 0.25, 0.25), fontsize=8)                 # needs sage.symbolic
            Multigraphics with 2 elements

        .. PLOT::

            f = (x**2*sin(1/x)).function(x)
            g1 = plot(f(x), (x, -2, 2), axes_labels=['$x$', '$y$'])
            g2 = plot(f(x), (x, -0.3, 0.3), axes_labels=['$x$', '$y$'], \
                      frame=True)
            sphinx_plot(g1.inset(g2, pos=(0.15, 0.7, 0.25, 0.25), fontsize=8))

        We can add another inset by invoking ``inset`` on the last output::

            sage: g1g2 = _                                                              # needs sage.symbolic
            sage: g3 = plot(f(x), (x, -0.05, 0.05), axes_labels=['$x$', '$y$'],         # needs sage.symbolic
            ....:           frame=True)
            sage: g1g2.inset(g3, pos=(0.65, 0.12, 0.25, 0.25))                          # needs sage.symbolic
            Multigraphics with 3 elements

        .. PLOT::

            f = (x**2*sin(1/x)).function(x)
            g1 = plot(f(x), (x, -2, 2), axes_labels=['$x$', '$y$'])
            g2 = plot(f(x), (x, -0.3, 0.3), axes_labels=['$x$', '$y$'], \
                      frame=True)
            g1g2 = g1.inset(g2, pos=(0.15, 0.7, 0.25, 0.25), fontsize=8)
            g3 = plot(f(x), (x, -0.05, 0.05), axes_labels=['$x$', '$y$'], \
                      frame=True)
            sphinx_plot(g1g2.inset(g3, pos=(0.65, 0.12, 0.25, 0.25)))
        """
        from .multigraphics import MultiGraphics
        if pos is None:
            pos = (0.7, 0.7, 0.2, 0.2)
        pos0 = (0.05, 0.05, 0.9, 0.9)
        if fontsize is not None:
            graphics._extra_kwds['fontsize'] = fontsize
        elif 'fontsize' not in graphics._extra_kwds:
            graphics._extra_kwds['fontsize'] = 6
        return MultiGraphics([(self, pos0), (graphics, pos)])


# Deprecation notice for GraphicsArray import
def GraphicsArray(*args, **kwargs):
    r"""
    This is deprecated (see :issue:`28675`).
    Use :class:`sage.plot.multigraphics.GraphicsArray` instead.

    TESTS::

        sage: from sage.plot.graphics import GraphicsArray
        sage: c = circle((0,0), 1)
        sage: G = GraphicsArray([c, c])
        doctest:...: DeprecationWarning: GraphicsArray must be imported from
        sage.plot.multigraphics and no longer from sage.plot.graphics.
        See https://github.com/sagemath/sage/issues/28675 for details.
        sage: G
        Graphics Array of size 1 x 2
    """
    from .multigraphics import GraphicsArray as NewGraphicsArray
    from sage.misc.superseded import deprecation
    deprecation(28675, "GraphicsArray must be imported from "
                "sage.plot.multigraphics and no longer from "
                "sage.plot.graphics.")
    return NewGraphicsArray(*args, **kwargs)
