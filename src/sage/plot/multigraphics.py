# sage_setup: distribution = sagemath-plot
# sage.doctest: needs sage.symbolic
r"""
Graphics arrays and insets

This module defines the classes :class:`MultiGraphics` and
:class:`GraphicsArray`. The class :class:`MultiGraphics` is the base class
for 2-dimensional graphical objects that are composed of various
:class:`~sage.plot.graphics.Graphics` objects, arranged in a given canvas.
The subclass :class:`GraphicsArray` is for
:class:`~sage.plot.graphics.Graphics` objects arranged in a regular array.

AUTHORS:

- Eric Gourgoulhon (2019-05-24): initial version, refactoring the class
  ``GraphicsArray`` that was defined in the module :mod:`~sage.plot.graphics`.
"""
import os
from sage.misc.fast_methods import WithEqualityById
from sage.structure.sage_object import SageObject
from sage.misc.temporary_file import tmp_filename
from .graphics import Graphics, ALLOWED_EXTENSIONS, _parse_figsize


class MultiGraphics(WithEqualityById, SageObject):
    r"""
    Base class for objects composed of :class:`~sage.plot.graphics.Graphics`
    objects.

    Both the display and the output to a file of ``MultiGraphics`` objects
    are governed by the method :meth:`save`, which is called by the rich output
    display manager, via
    :meth:`~sage.repl.rich_output.display_manager.DisplayManager.graphics_from_save`.

    The user interface is through the functions
    :func:`~sage.plot.plot.multi_graphics` (generic multi-graphics) and
    :func:`~sage.plot.plot.graphics_array` (subclass :class:`GraphicsArray`).

    INPUT:

    - ``graphics_list`` -- list of graphics along with their positions on the
      common canvas; each element of ``graphics_list`` is either

      - a pair ``(graphics, position)``, where ``graphics`` is a
        :class:`~sage.plot.graphics.Graphics` object and ``position`` is the
        4-tuple ``(left, bottom, width, height)`` specifying the location and
        size of the graphics on the canvas, all quantities being in fractions
        of the canvas width and height

      - or a single :class:`~sage.plot.graphics.Graphics` object; its position
        is then assumed to occupy the whole canvas, except for some padding;
        this corresponds to the default position
        ``(left, bottom, width, height) = (0.125, 0.11, 0.775, 0.77)``

    EXAMPLES:

    A multi-graphics made from two graphics objects::

        sage: g1 = plot(sin(x^3), (x, -pi, pi))
        sage: g2 = circle((0,0), 1, color='red')
        sage: G = multi_graphics([g1, (g2, (0.2, 0.55, 0.3, 0.3))])
        sage: G
        Multigraphics with 2 elements

    .. PLOT::

        g1 = plot(sin(x**3), (x, -pi, pi))
        g2 = circle((0,0), 1, color='red')
        G = multi_graphics([g1, (g2, (0.2, 0.55, 0.3, 0.3))])
        sphinx_plot(G)

    Since no position was given for ``g1``, it occupies the whole canvas.
    Moreover, we note that ``g2`` has been drawn over ``g1`` with a white
    background. To have a transparent background instead, one has to construct
    ``g2`` with the keyword ``transparent`` set to ``True``::

        sage: g2 = circle((0,0), 1, color='red', transparent=True)
        sage: G = multi_graphics([g1, (g2, (0.2, 0.55, 0.3, 0.3))])
        sage: G
        Multigraphics with 2 elements

    .. PLOT::

        g1 = plot(sin(x**3), (x, -pi, pi))
        g2 = circle((0,0), 1, color='red', transparent=True)
        G = multi_graphics([g1, (g2, (0.2, 0.55, 0.3, 0.3))])
        sphinx_plot(G)

    We can add a new graphics object to G via the method :meth:`append`::

        sage: g3 = complex_plot(zeta, (-20, 10), (-20, 20),
        ....:                   axes_labels=['$x$', '$y$'], frame=True)
        sage: G.append(g3, pos=(0.63, 0.12, 0.3, 0.3))
        sage: G
        Multigraphics with 3 elements

    .. PLOT::

        g1 = plot(sin(x**3), (x, -pi, pi))
        g2 = circle((0,0), 1, color='red', transparent=True)
        G = multi_graphics([g1, (g2, (0.2, 0.55, 0.3, 0.3))])
        g3 = complex_plot(zeta, (-20, 10), (-20, 20), \
                          axes_labels=['$x$', '$y$'], frame=True)
        G.append(g3, pos=(0.63, 0.12, 0.3, 0.3))
        sphinx_plot(G)

    We can access the individual elements composing ``G`` with the
    square-bracket operator::

        sage: print(G[0])
        Graphics object consisting of 1 graphics primitive
        sage: G[0] is g1
        True
        sage: G[1] is g2
        True
        sage: G[2] is g3
        True

    ``G[:]`` returns the full list of graphics objects composing ``G``::

        sage: G[:]
        [Graphics object consisting of 1 graphics primitive,
         Graphics object consisting of 1 graphics primitive,
         Graphics object consisting of 1 graphics primitive]
        sage: len(G)
        3
    """
    def __init__(self, graphics_list):
        r"""
        Initialize the attributes common to all MultiGraphics objects.

        TESTS::

            sage: from sage.plot.multigraphics import MultiGraphics
            sage: G = MultiGraphics([])
            sage: print(G)
            Multigraphics with 0 element
            sage: c = circle((0,0), 1)
            sage: G = MultiGraphics([c, (c, (0.7, 0.6, 0.2, 0.2))])
            sage: print(G)
            Multigraphics with 2 elements
        """
        self._glist = []
        self._positions = []

        for ins in graphics_list:
            if isinstance(ins, Graphics):
                self.append(ins)  # default position
            else:
                if not isinstance(ins, (list, tuple)) or len(ins) != 2:
                    raise TypeError("a pair (Graphics, position) is "
                                    f"expected, not {ins}")
                self.append(ins[0], pos=ins[1])

    def _repr_(self):
        r"""
        Representation of ``self``.

        EXAMPLES::

            sage: c = circle((0,0), 1)
            sage: G = graphics_array([c, c, c])
            sage: G._repr_()
            'Graphics Array of size 1 x 3'
            sage: G
            Graphics Array of size 1 x 3
        """
        return str(self)

    def _rich_repr_(self, display_manager, **kwds):
        r"""
        Rich Output Magic Method.

        See :mod:`sage.repl.rich_output` for details.

        .. TODO::

           This method is identical to Graphics._rich_repr_ so it could be
           inherited from a common base class

        EXAMPLES::

            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: G = graphics_array([Graphics(), Graphics()], 1, 2)
            sage: G._rich_repr_(dm)
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

    def __getitem__(self, i):
        r"""
        Return the ``i``-th element of the list of graphics composing ``self``.

        EXAMPLES:

        We can access and view individual plots::

            sage: L = [[plot(x^2)], [plot(x^3)]]
            sage: G = graphics_array(L)
            sage: G[1]
            Graphics object consisting of 1 graphics primitive

        Another example::

            sage: L = [plot(sin(k*x), (x,-pi,pi)) + circle((k,k), 1,
            ....:           color='red') for k in range(10)]
            sage: G = graphics_array(L, 5, 2)
            sage: G[3]
            Graphics object consisting of 2 graphics primitives
        """
        return self._glist[i]

    def __setitem__(self, i, g):
        r"""
        Set the ``i``-th element of the list of graphics composing ``self``.

        EXAMPLES::

            sage: L = [[plot(x^2)], [plot(x^3)]]
            sage: G = graphics_array(L)
            sage: G[1] # the plot of x^3
            Graphics object consisting of 1 graphics primitive

        Now we change it::

            sage: G[1] = circle((1,1), 2) + points([(1,2), (3,2), (5,5)],
            ....:                                  color='purple')
            sage: G[1] # a circle and some purple points
            Graphics object consisting of 2 graphics primitives
        """
        self._glist[i] = g

    def __len__(self):
        r"""
        Total number of Graphics objects composing ``self``.

        EXAMPLES::

            sage: L = [circle((0,0), n) for n in range(6)]
            sage: G = graphics_array(L, 2, 3)
            sage: len(G)
            6
        """
        return len(self._glist)

    def matplotlib(self, figure=None, figsize=None, **kwds):
        r"""
        Construct or modify a Matplotlib figure by drawing ``self`` on it.

        INPUT:

        - ``figure`` -- (default: ``None``) Matplotlib figure (class
          ``matplotlib.figure.Figure``) on which ``self`` is to be displayed;
          if ``None``, the figure will be created from the parameter
          ``figsize``

        - ``figsize`` -- (default: ``None``) width or [width, height] in inches
          of the Matplotlib figure in case ``figure`` is ``None``; if
          ``figsize`` is ``None``, Matplotlib's default (6.4 x 4.8 inches) is
          used

        - ``kwds`` -- options passed to the
          :meth:`~sage.plot.graphics.Graphics.matplotlib` method of
          each graphics object constituting ``self``

        OUTPUT:

        - a ``matplotlib.figure.Figure`` object; if the argument ``figure`` is
          provided, this is the same object as ``figure``.

        EXAMPLES:

        Let us consider a :class:`GraphicsArray` object with 3 elements::

            sage: G = graphics_array([plot(sin(x^k), (x, 0, 3))
            ....:                     for k in range(1, 4)])

        If ``matplotlib()`` is invoked without any argument, a Matplotlib
        figure is created and contains the 3 graphics element of the array
        as 3 Matplotlib ``Axes``::

            sage: fig = G.matplotlib()
            sage: fig
            <Figure size 640x480 with 3 Axes>
            sage: type(fig)
            <class 'matplotlib.figure.Figure'>

        Specifying the figure size (in inches)::

            sage: G.matplotlib(figsize=(8., 5.))
            <Figure size 800x500 with 3 Axes>

        If a single number is provided for ``figsize``, it is considered to be
        the width; the height is then computed according to Matplotlib's
        default aspect ratio (4/3)::

            sage: G.matplotlib(figsize=8.)
            <Figure size 800x600 with 3 Axes>

        An example of use with a preexisting created figure, created by
        ``pyplot``::

            sage: import matplotlib.pyplot as plt
            sage: fig1 = plt.figure(1)
            sage: fig1
            <Figure size 640x480 with 0 Axes>
            sage: fig_out = G.matplotlib(figure=fig1)
            sage: fig_out
            <Figure size 640x480 with 3 Axes>

        Note that the output figure is the same object as the input one::

            sage: fig_out is fig1
            True

        It has however been modified by ``G.matplotlib(figure=fig1)``, which
        has added 3 new ``Axes`` to it.

        Another example, with a figure created from scratch, via Matplolib's
        ``Figure``::

            sage: from matplotlib.figure import Figure
            sage: fig2 = Figure()
            sage: fig2
            <Figure size 640x480 with 0 Axes>
            sage: G.matplotlib(figure=fig2)
            <Figure size 640x480 with 3 Axes>
            sage: fig2
            <Figure size 640x480 with 3 Axes>
        """
        from matplotlib.figure import Figure
        glist = self._glist
        if len(glist) == 0:       # for an empty MultiGraphics, we create
            glist = [Graphics()]  # a 1-element list with an empty graphics
        # If no Matplotlib figure is provided, it is created here:
        if figure is None:
            if figsize is not None:
                figsize = _parse_figsize(figsize)
            figure = Figure(figsize=figsize)
        global do_verify
        do_verify = True
        for i, g in enumerate(glist):
            # Options for g.matplotlib():
            options = {}
            options.update(Graphics.SHOW_OPTIONS)  # default options for show()
            options['legend_options'] = Graphics.LEGEND_OPTIONS  # default leg.
            options.update(g._extra_kwds)  # options set in g
            options.update(kwds)
            # We get rid of options that are not relevant for g.matplotlib():
            options.pop('dpi', None)
            options.pop('fig_tight', None)
            transparent = options.pop('transparent', None)
            # Creating the Matplotlib Axes object "subplot" on the figure:
            subplot = self._add_subplot(figure, i)
            # and drawing g on it:
            g.matplotlib(figure=figure, sub=subplot, verify=do_verify,
                         **options)
            if transparent:
                subplot.set_facecolor('none')
        return figure

    def save(self, filename, figsize=None, **kwds):
        r"""
        Save ``self`` to a file, in various formats.

        INPUT:

        - ``filename`` -- string; the file name. The image format is given by
          the extension, which can be one of the following:

            * ``.eps``,

            * ``.pdf``,

            * ``.png``,

            * ``.ps``,

            * ``.sobj`` (for a Sage object you can load later),

            * ``.svg``,

            * empty extension will be treated as ``.sobj``.

        - ``figsize`` -- (default: ``None``) width or [width, height] in inches
          of the Matplotlib figure; if none is provided, Matplotlib's default
          (6.4 x 4.8 inches) is used

        - ``kwds`` -- keyword arguments, like ``dpi=...``, passed to the
          plotter, see :meth:`show`

        EXAMPLES::

            sage: F = tmp_filename(ext='.png')
            sage: L = [plot(sin(k*x), (x,-pi,pi)) for k in [1..3]]
            sage: G = graphics_array(L)
            sage: G.save(F, dpi=500, axes=False)

        TESTS::

            sage: graphics_array([]).save(F)
            sage: graphics_array([[]]).save(F)
        """
        from matplotlib import rcParams
        ext = os.path.splitext(filename)[1].lower()
        if ext in ['', '.sobj']:
            SageObject.save(self, filename)
        elif ext not in ALLOWED_EXTENSIONS:
            raise ValueError("allowed file extensions for images are '" +
                             "', '".join(ALLOWED_EXTENSIONS) + "'!")
        else:
            rc_backup = (rcParams['ps.useafm'], rcParams['pdf.use14corefonts'],
                         rcParams['text.usetex'])  # save the rcParams
            figure = self.matplotlib(figsize=figsize, **kwds)
            transparent = kwds.get('transparent',
                                   Graphics.SHOW_OPTIONS['transparent'])
            fig_tight = kwds.get('fig_tight',
                                 Graphics.SHOW_OPTIONS['fig_tight'])
            dpi = kwds.get('dpi', Graphics.SHOW_OPTIONS['dpi'])
            # One can output in PNG, PS, EPS, PDF, PGF, or SVG format,
            # depending on the file extension.
            # PGF is handled by a different backend
            if ext == '.pgf':
                from sage.features.latex import xelatex,pdflatex,lualatex
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
                    # Matplotlib documentation:
                    # https://matplotlib.org/users/pgf.html#pgf-tutorial
                    pgf_options = {"pgf.texsystem": "pdflatex",
                                   "pgf.preamble": [
                                      r"\usepackage[utf8x]{inputenc}",
                                      r"\usepackage[T1]{fontenc}"
                                      ]}
                else:
                    pgf_options = {"pgf.texsystem": latex_implementations[0]}
                rcParams.update(pgf_options)
                from matplotlib.backends.backend_pgf import FigureCanvasPgf
                figure.set_canvas(FigureCanvasPgf(figure))
            # Matplotlib looks at the file extension to see what the renderer
            # should be. The default is FigureCanvasAgg for PNG's because this
            # is by far the most common type of files rendered, like in the
            # notebook, for example. If the file extension is not '.png', then
            # Matplotlib will handle it.
            else:
                from matplotlib.backends.backend_agg import FigureCanvasAgg
                figure.set_canvas(FigureCanvasAgg(figure))
            if isinstance(self, GraphicsArray):
                # tight_layout adjusts the *subplot* parameters so ticks aren't
                # cut off, etc.
                figure.tight_layout()
            opts = {"dpi": dpi, "transparent": transparent}
            if fig_tight is True:
                opts['bbox_inches'] = 'tight'
            figure.savefig(filename, **opts)
            # Restore the rcParams to the original, possibly user-set values
            (rcParams['ps.useafm'], rcParams['pdf.use14corefonts'],
             rcParams['text.usetex']) = rc_backup

    def save_image(self, filename=None, *args, **kwds):
        r"""
        Save an image representation of ``self``.  The image type is
        determined by the extension of the filename.  For example,
        this could be ``.png``, ``.jpg``, ``.gif``, ``.pdf``,
        ``.svg``.  Currently this is implemented by calling the
        :meth:`save` method of self, passing along all arguments and
        keywords.

        .. NOTE::

            Not all image types are necessarily implemented for all
            graphics types.  See :meth:`save` for more details.

        EXAMPLES::

            sage: plots = [[plot(m*cos(x + n*pi/4), (x, 0, 2*pi))
            ....:           for n in range(3)] for m in range(1,3)]
            sage: G = graphics_array(plots)
            sage: G.save_image(tmp_filename(ext='.png'))
        """
        self.save(filename, *args, **kwds)

    def _latex_(self, **kwds):
        r"""
        Return a string plotting ``self`` with PGF.

        INPUT:

        - ``**kwds`` -- all keyword arguments will be passed to the plotter

        OUTPUT: string of PGF commands to plot ``self``

        EXAMPLES::

            sage: A = graphics_array([plot(sin), plot(cos)])
            sage: A._latex_()[:40]  # not tested (see comment below)
            '%% Creator: Matplotlib, PGF backend\n%%\n%'

        The above doctest fails on macOS due to the following Matplotlib issue: https://github.com/matplotlib/matplotlib/issues/10307
        """
        tmpfilename = tmp_filename(ext='.pgf')
        self.save(filename=tmpfilename, **kwds)
        with open(tmpfilename) as tmpfile:
            latex_list = tmpfile.readlines()
        return ''.join(latex_list)

    def show(self, **kwds):
        r"""
        Show ``self`` immediately.

        This method attempts to display the graphics immediately,
        without waiting for the currently running code (if any) to
        return to the command line. Be careful, calling it from within
        a loop will potentially launch a large number of external
        viewer programs.

        OPTIONAL INPUT:

        - ``dpi`` -- dots per inch

        - ``figsize`` -- width or [width, height] of the figure, in inches; the
          default is 6.4 x 4.8 inches

        - ``axes`` -- boolean; if ``True``, all individual graphics are
          endowed with axes; if ``False``, all axes are removed (this overrides
          the ``axes`` option set in each graphics)

        - ``frame`` -- boolean; if ``True``, all individual graphics are
          drawn with a frame around them; if ``False``, all frames are removed
          (this overrides the ``frame`` option set in each graphics)

        - ``fontsize`` -- positive integer, the size of fonts for the axes
          labels (this overrides the ``fontsize`` option set in each graphics)

        OUTPUT:

        This method does not return anything. Use :meth:`save` if you
        want to save the figure as an image.

        EXAMPLES:

        This draws a graphics array with four trig plots and no axes in any of
        the plots and a figure width of 4 inches::

            sage: G = graphics_array([[plot(sin), plot(cos)],
            ....:                     [plot(tan), plot(sec)]])
            sage: G.show(axes=False, figsize=4)

        .. PLOT::

            G = graphics_array([[plot(sin), plot(cos)], \
                                [plot(tan), plot(sec)]])
            sphinx_plot(G, axes=False, figsize=4)

        Same thing with a frame around each individual graphics::

            sage: G.show(axes=False, frame=True, figsize=4)

        .. PLOT::

            G = graphics_array([[plot(sin), plot(cos)], \
                                [plot(tan), plot(sec)]])
            sphinx_plot(G, axes=False, frame=True, figsize=4)

        Actually, many options are possible; for instance, we may set
        ``fontsize`` and ``gridlines``::

            sage: G.show(axes=False, frame=True, figsize=4, fontsize=8,
            ....:        gridlines='major')

        .. PLOT::

            G = graphics_array([[plot(sin), plot(cos)], \
                                [plot(tan), plot(sec)]])
            sphinx_plot(G, axes=False, frame=True, figsize=4, fontsize=8, \
                        gridlines='major')
        """
        from sage.repl.rich_output import get_display_manager
        dm = get_display_manager()
        dm.display_immediately(self, **kwds)

    def plot(self):
        r"""
        Return ``self`` since ``self`` is already a graphics object.

        EXAMPLES::

            sage: g1 = plot(cos, 0, 1)
            sage: g2 = circle((0,0), 1)
            sage: G = multi_graphics([g1, g2])
            sage: G.plot() is G
            True
        """
        return self

    def inset(self, graphics, pos=None, fontsize=None):
        r"""
        Add a graphics object as an inset.

        INPUT:

        - ``graphics`` -- the graphics object (instance of :class:`Graphics`)
          to be added as an inset

        - ``pos`` -- (default: ``None``) 4-tuple
          ``(left, bottom, width, height)`` specifying the location and
          relative size of the inset on the canvas, all quantities being
          expressed in fractions of the canvas width and height; if ``None``,
          the value ``(0.7, 0.7, 0.2, 0.2)`` is used

        - ``fontsize`` -- (default: ``None``)  integer, font size (in points)
          for the inset; if ``None``, the value of 6 points is used, unless
          ``fontsize`` has been explicitly set in the construction of
          ``graphics`` (in this case, it is not overwritten here)

        OUTPUT: instance of :class:`~sage.plot.multigraphics.MultiGraphics`

        EXAMPLES:

        Let us consider a graphics array of 2 elements::

            sage: G = graphics_array([plot(sin, (0, 2*pi)),
            ....:                     plot(cos, (0, 2*pi))])
            sage: G
            Graphics Array of size 1 x 2

        .. PLOT::

            G = graphics_array([plot(sin, (0, 2*pi)), plot(cos, (0, 2*pi))])
            sphinx_plot(G)

        and add some inset at the default position::

            sage: c = circle((0,0), 1, color='red', thickness=2, frame=True)
            sage: G.inset(c)
            Multigraphics with 3 elements

        .. PLOT::

            G = graphics_array([plot(sin, (0, 2*pi)), plot(cos, (0, 2*pi))])
            c = circle((0,0), 1, color='red', thickness=2, frame=True)
            sphinx_plot(G.inset(c))

        We may customize the position and font size of the inset::

            sage: G.inset(c, pos=(0.3, 0.7, 0.2, 0.2), fontsize=8)
            Multigraphics with 3 elements

        .. PLOT::

            G = graphics_array([plot(sin, (0, 2*pi)), plot(cos, (0, 2*pi))])
            c = circle((0,0), 1, color='red', thickness=2, frame=True)
            sphinx_plot(G.inset(c, pos=(0.3, 0.7, 0.2, 0.2), fontsize=8))
        """
        if pos is None:
            pos = (0.7, 0.7, 0.2, 0.2)
        if fontsize is not None:
            graphics._extra_kwds['fontsize'] = fontsize
        elif 'fontsize' not in graphics._extra_kwds:
            graphics._extra_kwds['fontsize'] = 6
        current = []  # list of current pairs (graphics, position)
        for i, g in enumerate(self._glist):
            current.append((g, self.position(i)))
        resu = MultiGraphics(current)
        resu.append(graphics, pos=pos)
        return resu

    #
    # Methods to reimplement in derived classes:
    #
    def __str__(self):
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: from sage.plot.multigraphics import MultiGraphics
            sage: G = MultiGraphics([])
            sage: G.__str__()
            'Multigraphics with 0 element'
            sage: str(G)
            'Multigraphics with 0 element'
            sage: c = circle((0,0), 1)
            sage: G = MultiGraphics([c])
            sage: str(G)
            'Multigraphics with 1 element'
            sage: G = MultiGraphics([c, c, c])
            sage: str(G)
            'Multigraphics with 3 elements'
        """
        n = len(self._glist)
        if n <= 1:
            return f"Multigraphics with {n} element"
        return f"Multigraphics with {n} elements"

    def _add_subplot(self, figure, index, **options):
        r"""
        Add a subplot to a given Matplotlib ``Figure``, the position of
        which is governed by a given element of ``self``.

        This method encapsulates the Matplotlib method ``Figure.add_axes``
        and is intended to be called by :meth:`MultiGraphics.save`.

        INPUT:

        - ``figure`` -- a Matplotlib ``Figure`` object
        - ``index`` -- integer specifying the element of ``self``
        - ``options`` -- extra options to be passed to ``Figure.add_axes``

        OUTPUT: a Matplotlib ``Axes`` object

        EXAMPLES::

            sage: g0 = circle((0,0), 1)
            sage: g1 = plot(sin)
            sage: G = multi_graphics([g0, (g1, (0.2, 0.3, 0.4, 0.1))])
            sage: from matplotlib.figure import Figure
            sage: fig = Figure()
            sage: fig
            <Figure size 640x480 with 0 Axes>
            sage: ax0 = G._add_subplot(fig, 0)
            sage: type(ax0)
            <class 'matplotlib.axes._axes.Axes'>
            sage: fig
            <Figure size 640x480 with 1 Axes>
            sage: ax1 = G._add_subplot(fig, 1)
            sage: fig
            <Figure size 640x480 with 2 Axes>

        TESTS::

            sage: [ax0, ax1] == fig.get_axes()
            True
            sage: G.position(1)
            (0.2, 0.3, 0.4, 0.1)
            sage: import numpy  # to ensure numpy 2.0 compatibility
            sage: if int(numpy.version.short_version[0]) > 1:
            ....:     numpy.set_printoptions(legacy="1.25")
            sage: ax1.get_position().bounds  # tol 1.0e-13
            (0.2, 0.3, 0.4000000000000001, 0.10000000000000003)
        """
        # Note: using label=str(index) ensures that a new Axes is generated
        # for each element of ``self``, even if some elements share the same
        # positions
        return figure.add_axes(self._positions[index], label=str(index),
                               **options)

    def position(self, index):
        r"""
        Return the position and relative size of an element of ``self`` on the
        canvas.

        INPUT:

        - ``index`` -- integer specifying which element of ``self``

        OUTPUT:

        - a 4-tuple ``(left, bottom, width, height)`` giving the location and
          relative size of the element on the canvas, all quantities being
          expressed in fractions of the canvas width and height

        EXAMPLES::

            sage: g1 = plot(sin(x^2), (x, 0, 4))
            sage: g2 = circle((0,0), 1, rgbcolor='red', fill=True, axes=False)
            sage: G = multi_graphics([g1, (g2, (0.15, 0.2, 0.1, 0.15))])
            sage: G.position(0)  # tol 1.0e-13
            (0.125, 0.11, 0.775, 0.77)
            sage: G.position(1)  # tol 1.0e-13
            (0.15, 0.2, 0.1, 0.15)
        """
        return self._positions[index]

    def append(self, graphics, pos=None):
        r"""
        Append a graphics object to ``self``.

        INPUT:

        - ``graphics`` -- the graphics object (instance of :class:`Graphics`)
          to be added to ``self``

        - ``pos`` -- (default: ``None``) 4-tuple
          ``(left, bottom, width, height)`` specifying the location and size
          of ``graphics`` on the canvas, all quantities being in fractions of
          the canvas width and height; if ``None``, ``graphics`` is assumed to
          occupy the whole canvas, except for some padding; this corresponds to
          the default position
          ``(left, bottom, width, height) = (0.125, 0.11, 0.775, 0.77)``

        EXAMPLES:

        Let us consider a multigraphics with 2 elements::

            sage: g1 = plot(chebyshev_T(4, x), (x, -1, 1), title='n=4')
            sage: g2 = plot(chebyshev_T(8, x), (x, -1, 1), title='n=8',
            ....:           color='red')
            sage: G = multi_graphics([(g1, (0.125, 0.2, 0.4, 0.4)),
            ....:                     (g2, (0.55, 0.4, 0.4, 0.4))])
            sage: G
            Multigraphics with 2 elements

        .. PLOT::

            g1 = plot(chebyshev_T(4, x), (x, -1, 1), title='n=4')
            g2 = plot(chebyshev_T(8, x), (x, -1, 1), title='n=8', color='red')
            G = multi_graphics([(g1, (0.125, 0.2, 0.4, 0.4)), \
                                (g2, (0.55, 0.4, 0.4, 0.4))])
            sphinx_plot(G)

        We append a third plot to it::

            sage: g3 = plot(chebyshev_T(16, x), (x, -1, 1), title='n=16',
            ....:           color='brown')
            sage: G.append(g3, pos=(0.55, 0.11, 0.4, 0.15))
            sage: G
            Multigraphics with 3 elements

        .. PLOT::

            g1 = plot(chebyshev_T(4, x), (x, -1, 1), title='n=4')
            g2 = plot(chebyshev_T(8, x), (x, -1, 1), title='n=8', color='red')
            G = multi_graphics([(g1, (0.125, 0.2, 0.4, 0.4)), \
                                (g2, (0.55, 0.4, 0.4, 0.4))])
            g3 = plot(chebyshev_T(16, x), (x, -1, 1), title='n=16', \
                      color='brown')
            G.append(g3, pos=(0.55, 0.11, 0.4, 0.15))
            sphinx_plot(G)

        We may use ``append`` to add a title::

            sage: title = text("Chebyshev polynomials", (0, 0), fontsize=16,
            ....:              axes=False)
            sage: G.append(title, pos=(0.18, 0.8, 0.7, 0.1))
            sage: G
            Multigraphics with 4 elements

        .. PLOT::

            g1 = plot(chebyshev_T(4, x), (x, -1, 1), title='n=4')
            g2 = plot(chebyshev_T(8, x), (x, -1, 1), title='n=8', color='red')
            G = multi_graphics([(g1, (0.125, 0.2, 0.4, 0.4)), \
                                (g2, (0.55, 0.4, 0.4, 0.4))])
            g3 = plot(chebyshev_T(16, x), (x, -1, 1), title='n=16', \
                      color='brown')
            G.append(g3, pos=(0.55, 0.11, 0.4, 0.15))
            title = text("Chebyshev polynomials", (0, 0), fontsize=16, \
                         axes=False)
            G.append(title, pos=(0.18, 0.8, 0.7, 0.1))
            sphinx_plot(G)

        .. SEEALSO::

            :meth:`inset`
        """
        from matplotlib import rcParams
        if not isinstance(graphics, Graphics):
            raise TypeError("a Graphics object is expected, "
                            f"not {graphics}")
        if pos is None:
            # Default position:
            left = rcParams['figure.subplot.left']
            bottom = rcParams['figure.subplot.bottom']
            width = rcParams['figure.subplot.right'] - left
            height = rcParams['figure.subplot.top'] - bottom
            pos = (left, bottom, width, height)
        elif not isinstance(pos, (list, tuple)) or len(pos) != 4:
            raise TypeError(f"pos must be a 4-tuple, not {pos}")
        pos = tuple(float(p) for p in pos)
        self._glist.append(graphics)
        self._positions.append(pos)


# ****************************************************************************


class GraphicsArray(MultiGraphics):
    r"""
    This class implements 2-dimensional graphical objects that constitute
    an array of :class:`~sage.plot.graphics.Graphics` drawn on a single
    canvas.

    The user interface is through the function
    :func:`~sage.plot.plot.graphics_array`.

    INPUT:

    - ``array`` -- either a list of lists of
      :class:`~sage.plot.graphics.Graphics` elements (generic case) or a
      single list of :class:`~sage.plot.graphics.Graphics` elements (case of a
      single-row array)

    EXAMPLES:

    An array made of four graphics objects::

        sage: g1 = plot(sin(x^2), (x, 0, 6), axes_labels=['$x$', '$y$'],
        ....:           axes=False, frame=True, gridlines='minor')
        sage: y = var('y')
        sage: g2 = streamline_plot((sin(x), cos(y)), (x,-3,3), (y,-3,3),
        ....:                      aspect_ratio=1)
        sage: g3 = graphs.DodecahedralGraph().plot()
        sage: g4 = polar_plot(sin(5*x)^2, (x, 0, 2*pi), color='green',
        ....:                 fontsize=8) \
        ....:      + circle((0,0), 0.5, rgbcolor='red', fill=True, alpha=0.1,
        ....:               legend_label='pink')
        sage: g4.set_legend_options(loc='upper right')
        sage: G = graphics_array([[g1, g2], [g3, g4]])
        sage: G
        Graphics Array of size 2 x 2

    .. PLOT::

        g1 = plot(sin(x**2), (x, 0, 6), axes_labels=['$x$', '$y$'], \
                  axes=False, frame=True, gridlines='minor')
        y = var('y')
        g2 = streamline_plot((sin(x), cos(y)), (x,-3,3), (y,-3,3), \
                             aspect_ratio=1)
        g3 = graphs.DodecahedralGraph().plot()
        g4 = polar_plot(sin(5*x)**2, (x, 0, 2*pi), color='green', fontsize=8) \
             + circle((0,0), 0.5, rgbcolor='red', fill=True, alpha=0.1, \
                      legend_label='pink')
        g4.set_legend_options(loc='upper right')
        G = graphics_array([[g1, g2], [g3, g4]])
        sphinx_plot(G)

    If one constructs the graphics array from a single list of graphics
    objects, one obtains a single-row array::

        sage: G = graphics_array([g1, g2, g3, g4])
        sage: G
        Graphics Array of size 1 x 4

    .. PLOT::

        g1 = plot(sin(x**2), (x, 0, 6), axes_labels=['$x$', '$y$'], \
                  axes=False, frame=True, gridlines='minor')
        y = var('y')
        g2 = streamline_plot((sin(x), cos(y)), (x,-3,3), (y,-3,3), \
                             aspect_ratio=1)
        g3 = graphs.DodecahedralGraph().plot()
        g4 = polar_plot(sin(5*x)**2, (x, 0, 2*pi), color='green', fontsize=8) \
             + circle((0,0), 0.5, rgbcolor='red', fill=True, alpha=0.1, \
                      legend_label='pink')
        g4.set_legend_options(loc='upper right')
        G = graphics_array([g1, g2, g3, g4])
        sphinx_plot(G)

    We note that the overall aspect ratio of the figure is 4/3 (the default),
    which makes ``g1`` elongated, while the aspect ratio of ``g2``, which has
    been specified with the parameter ``aspect_ratio=1`` is preserved. To get
    a better aspect ratio for the whole figure, one can use the option
    ``figsize`` in the method :meth:`~MultiGraphics.show`::

        sage: G.show(figsize=[8, 3])

    .. PLOT::

        g1 = plot(sin(x**2), (x, 0, 6), axes_labels=['$x$', '$y$'], \
                  axes=False, frame=True, gridlines='minor')
        y = var('y')
        g2 = streamline_plot((sin(x), cos(y)), (x,-3,3), (y,-3,3), \
                             aspect_ratio=1)
        g3 = graphs.DodecahedralGraph().plot()
        g4 = polar_plot(sin(5*x)**2, (x, 0, 2*pi), color='green', fontsize=8) \
             + circle((0,0), 0.5, rgbcolor='red', fill=True, alpha=0.1, \
                      legend_label='pink')
        g4.set_legend_options(loc='upper right')
        G = graphics_array([g1, g2, g3, g4])
        sphinx_plot(G, figsize=[8, 3])

    We can access individual elements of the graphics array with the
    square-bracket operator::

        sage: G = graphics_array([[g1, g2], [g3, g4]])  # back to the 2x2 array
        sage: print(G)
        Graphics Array of size 2 x 2
        sage: G[0] is g1
        True
        sage: G[1] is g2
        True
        sage: G[2] is g3
        True
        sage: G[3] is g4
        True

    Note that with respect to the square-bracket operator, ``G`` is considered
    as a flattened list of graphics objects, not as an array. For instance,
    ``G[0, 1]`` throws an error::

        sage: G[0, 1]
        Traceback (most recent call last):
        ...
        TypeError: list indices must be integers or slices, not tuple

    ``G[:]`` returns the full (flattened) list of graphics objects composing
    ``G``::

        sage: G[:]
        [Graphics object consisting of 1 graphics primitive,
        Graphics object consisting of 1 graphics primitive,
        Graphics object consisting of 51 graphics primitives,
        Graphics object consisting of 2 graphics primitives]

    The total number of Graphics objects composing the array is returned
    by the function ``len``::

        sage: len(G)
        4

    The square-bracket operator can be used to replace elements in the array::

        sage: G[0] = g4
        sage: G
        Graphics Array of size 2 x 2

    .. PLOT::

        g1 = plot(sin(x**2), (x, 0, 6), axes_labels=['$x$', '$y$'], \
                  axes=False, frame=True, gridlines='minor')
        y = var('y')
        g2 = streamline_plot((sin(x), cos(y)), (x,-3,3), (y,-3,3), \
                             aspect_ratio=1)
        g3 = graphs.DodecahedralGraph().plot()
        g4 = polar_plot(sin(5*x)**2, (x, 0, 2*pi), color='green', fontsize=8) \
             + circle((0,0), 0.5, rgbcolor='red', fill=True, alpha=0.1, \
                      legend_label='pink')
        g4.set_legend_options(loc='upper right')
        G = graphics_array([[g1, g2], [g3, g4]])
        G[0] = g4
        sphinx_plot(G)
    """
    def __init__(self, array):
        r"""
        Construct a ``GraphicsArray``.

        TESTS::

            sage: from sage.plot.multigraphics import GraphicsArray
            sage: g = circle((0,0), 1)  # a Graphics object
            sage: G = GraphicsArray([[g, g], [g, g]])
            sage: print(G)
            Graphics Array of size 2 x 2

        Construction from a single list ==> 1-row array::

            sage: G = GraphicsArray([g, g, g])
            sage: print(G)
            Graphics Array of size 1 x 3
            sage: G = GraphicsArray([g])
            sage: print(G)
            Graphics Array of size 1 x 1

        Empty array::

            sage: G = GraphicsArray([])
            sage: print(G)
            Graphics Array of size 0 x 0
            sage: len(G)
            0
            sage: G = GraphicsArray([[]])
            sage: print(G)
            Graphics Array of size 1 x 0
            sage: len(G)
            0

        Check treatment of wrong inputs::

            sage: G = GraphicsArray([[g, g], [g]])
            Traceback (most recent call last):
            ...
            TypeError: array must be a list of equal-size lists of Graphics
             objects, not [[Graphics object consisting of 1 graphics primitive,
             Graphics object consisting of 1 graphics primitive],
             [Graphics object consisting of 1 graphics primitive]]
            sage: G = GraphicsArray(g)
            Traceback (most recent call last):
            ...
            TypeError: array must be a list of lists of Graphics objects, not
             Graphics object consisting of 1 graphics primitive
            sage: G = GraphicsArray([g, x])
            Traceback (most recent call last):
            ...
            TypeError: every element of array must be a Graphics object
        """
        MultiGraphics.__init__(self, [])
        if not isinstance(array, (list, tuple)):
            raise TypeError("array must be a list of lists of Graphics "
                            f"objects, not {array}")
        array = list(array)
        self._rows = len(array)
        if self._rows > 0:
            if not isinstance(array[0], (list, tuple)):
                array = [array]
                self._rows = 1
            self._cols = len(array[0])
        else:
            self._cols = 0
        for row in array:  # basically flatten the list
            if not isinstance(row, (list, tuple)) or len(row) != self._cols:
                raise TypeError("array must be a list of equal-size lists of "
                                f"Graphics objects, not {array}")
            for g in row:
                if not isinstance(g, Graphics):
                    raise TypeError("every element of array must be a "
                                    "Graphics object")
                self._glist.append(g)
        # self._positions is not initialized since most of the time, it is not
        # not used. It is required only by the method inset(); it is then
        # initialized by the method position().

    def __str__(self):
        r"""
        String representation of the graphics array.

        EXAMPLES::

            sage: c = circle((0,0), 1)
            sage: G = graphics_array([c]*6, 2, 3)
            sage: G.__str__()
            'Graphics Array of size 2 x 3'
            sage: str(G)
            'Graphics Array of size 2 x 3'
        """
        return f"Graphics Array of size {self._rows} x {self._cols}"

    def _add_subplot(self, figure, index, **options):
        r"""
        Add a subplot to a given Matplotlib ``Figure``, the position of
        which is governed by a given element of ``self``.

        This method encapsulates the Matplotlib method ``Figure.add_subplot``
        and is intended to be called by :meth:`MultiGraphics.save`.

        INPUT:

        - ``figure`` -- a Matplotlib ``Figure`` object
        - ``index`` -- integer specifying the element of ``self``
        - ``options`` -- extra options to be passed to ``Figure.add_subplot``

        OUTPUT: a Matplotlib ``Axes`` object

        EXAMPLES::

            sage: c = circle((0,0), 1)
            sage: G = graphics_array([c, c])
            sage: from matplotlib.figure import Figure
            sage: fig = Figure()
            sage: ax1 = G._add_subplot(fig, 0)
            sage: type(ax1)
            <class 'matplotlib.axes...'>
            sage: ax2 = G._add_subplot(fig, 1)
            sage: fig.get_axes() == [ax1, ax2]
            True
        """
        if self._rows == 0 or self._cols == 0:
            rows = 1
            cols = 1
        else:
            rows = self._rows
            cols = self._cols
        # index --> index + 1 for Figure.add_subplot:
        return figure.add_subplot(rows, cols, index + 1, **options)

    def nrows(self):
        r"""
        Number of rows of the graphics array.

        EXAMPLES::

            sage: R = rainbow(6)
            sage: L = [plot(x^n, (x,0,1), color=R[n]) for n in range(6)]
            sage: G = graphics_array(L, 2, 3)
            sage: G.nrows()
            2
            sage: graphics_array(L).nrows()
            1
        """
        return self._rows

    def ncols(self):
        r"""
        Number of columns of the graphics array.

        EXAMPLES::

            sage: R = rainbow(6)
            sage: L = [plot(x^n, (x,0,1), color=R[n]) for n in range(6)]
            sage: G = graphics_array(L, 2, 3)
            sage: G.ncols()
            3
            sage: graphics_array(L).ncols()
            6
        """
        return self._cols

    def append(self, g):
        r"""
        Append a graphics to the array.

        Currently not implemented.

        TESTS::

            sage: from sage.plot.multigraphics import GraphicsArray
            sage: c = circle((0,0), 1)
            sage: G = GraphicsArray([c, c])
            sage: G.append(c)
            Traceback (most recent call last):
            ...
            NotImplementedError: Appending to a graphics array is not yet
             implemented
        """
        # Not clear if there is a way to do this
        raise NotImplementedError('Appending to a graphics array is not '
                                  'yet implemented')

    def position(self, index):
        r"""
        Return the position and relative size of an element of ``self`` on the
        canvas.

        INPUT:

        - ``index`` -- integer specifying which element of ``self``

        OUTPUT:

        - a 4-tuple ``(left, bottom, width, height)`` giving the location and
          relative size of the element on the canvas, all quantities being
          expressed in fractions of the canvas width and height

        EXAMPLES::

            sage: g1 = plot(sin(x), (x, -pi, pi))
            sage: g2 = circle((0,1), 1.)
            sage: G = graphics_array([g1, g2])
            sage: import numpy  # to ensure numpy 2.0 compatibility
            sage: if int(numpy.version.short_version[0]) > 1:
            ....:     numpy.set_printoptions(legacy="1.25")
            sage: G.position(0)  # tol 5.0e-3
            (0.025045451349937315,
             0.03415488992713045,
             0.4489880779745068,
             0.9345951100728696)
            sage: G.position(1)  # tol 5.0e-3
            (0.5170637412999687,
             0.20212705964722733,
             0.4489880779745068,
             0.5986507706326758)
        """
        if not self._positions:
            # self._positions must be generated, by invoking get_position() on
            # each of the Axes of the Matplotlib figure corresponding to self:
            from matplotlib.backends.backend_agg import FigureCanvasAgg
            figure = self.matplotlib()
            figure.set_canvas(FigureCanvasAgg(figure))
            figure.tight_layout()
            axes = figure.get_axes()
            self._positions = [ax.get_position().bounds for ax in axes]
        return self._positions[index]
