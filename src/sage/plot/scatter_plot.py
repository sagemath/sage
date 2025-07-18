"""
Scatter plots
"""

# *****************************************************************************
#       Copyright (C) 2006 Alex Clemesha <clemesha@gmail.com>,
#                          William Stein <wstein@gmail.com>,
#                     2008 Mike Hansen <mhansen@gmail.com>,
#                     2009 Emily Kirkman
#                     2022 Javier Honrubia <jhonrubia6@alumno.uned.es>
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
from sage.misc.decorators import options


class ScatterPlot(GraphicPrimitive):
    """
    Scatter plot graphics primitive.

    Input consists of two lists/arrays of the same length, whose
    values give the horizontal and vertical coordinates of each
    point in the scatter plot.  Options may be passed in
    dictionary format.

    EXAMPLES::

        sage: from sage.plot.scatter_plot import ScatterPlot
        sage: ScatterPlot([0,1,2], [3.5,2,5.1], {'facecolor':'white', 'marker':'s'})
        Scatter plot graphics primitive on 3 data points
    """
    def __init__(self, xdata, ydata, options):
        """
        Scatter plot graphics primitive.

        EXAMPLES::

            sage: import numpy
            sage: from sage.plot.scatter_plot import ScatterPlot
            sage: ScatterPlot(numpy.array([0,1,2]), numpy.array([3.5,2,5.1]), {'facecolor':'white', 'marker':'s'})
            Scatter plot graphics primitive on 3 data points
        """
        self.xdata = xdata
        self.ydata = ydata
        GraphicPrimitive.__init__(self, options)

    def get_minmax_data(self):
        """
        Return a dictionary with the bounding box data.

        EXAMPLES::

            sage: s = scatter_plot([[0,1],[2,4],[3.2,6]])
            sage: d = s.get_minmax_data()
            sage: d['xmin']
            ...0.0...
            sage: d['ymin']
            ...1.0...
        """
        return {'xmin': self.xdata.min(),
                'xmax': self.xdata.max(),
                'ymin': self.ydata.min(),
                'ymax': self.ydata.max()}

    def _allowed_options(self):
        """
        Return the dictionary of allowed options for the scatter plot
        graphics primitive.

        EXAMPLES::

            sage: from sage.plot.scatter_plot import ScatterPlot
            sage: list(sorted(ScatterPlot([-1,2], [17,4], {})._allowed_options().items()))
            [('alpha', 'How transparent the marker border is.'),
            ('clip', 'Whether or not to clip.'),
            ('edgecolor', 'The color of the marker border.'),
            ('facecolor', 'The color of the marker face.'),
            ('hue', 'The color given as a hue.'),
            ('marker', 'What shape to plot the points. See the documentation of plot() for the full list of markers.'),
            ('markersize', 'the size of the markers.'),
            ('rgbcolor', 'The color as an RGB tuple.'),
            ('zorder', 'The layer level in which to draw.')]
        """
        return {'markersize': 'the size of the markers.',
                'marker': 'What shape to plot the points. See the documentation of plot() for the full list of markers.',
                'alpha': 'How transparent the marker border is.',
                'rgbcolor': 'The color as an RGB tuple.',
                'hue': 'The color given as a hue.',
                'facecolor': 'The color of the marker face.',
                'edgecolor': 'The color of the marker border.',
                'zorder': 'The layer level in which to draw.',
                'clip': 'Whether or not to clip.'}

    def _repr_(self):
        """
        Text representation of a scatter plot graphics primitive.

        EXAMPLES::

            sage: import numpy
            sage: from sage.plot.scatter_plot import ScatterPlot
            sage: ScatterPlot(numpy.array([0,1,2]), numpy.array([3.5,2,5.1]), {})
            Scatter plot graphics primitive on 3 data points
        """
        return 'Scatter plot graphics primitive on %s data points' % len(self.xdata)

    def _render_on_subplot(self, subplot):
        """
        Render this scatter plot in a subplot.  This is the key function that
        defines how this scatter plot graphics primitive is rendered in
        matplotlib's library.

        EXAMPLES::

            sage: scatter_plot([[0,1],[2,2],[4.3,1.1]], marker='s')
            Graphics object consisting of 1 graphics primitive

        ::

            sage: scatter_plot([[n,n] for n in range(5)])
            Graphics object consisting of 1 graphics primitive
        """
        options = self.options()
        p = subplot.scatter(self.xdata, self.ydata, alpha=options['alpha'],
                            zorder=options['zorder'], marker=options['marker'],
                            s=options['markersize'], facecolors=options['facecolor'],
                            edgecolors=options['edgecolor'], clip_on=options['clip'])
        if not options['clip']:
            self._bbox_extra_artists = [p]


@options(alpha=1, markersize=50, marker='o', zorder=5, facecolor='#fec7b8', edgecolor='black', clip=True, aspect_ratio='automatic')
def scatter_plot(datalist, **options):
    """
    Return a Graphics object of a scatter plot containing all points in
    the datalist.  Type ``scatter_plot.options`` to see all available
    plotting options.

    INPUT:

    - ``datalist`` -- list of tuples ``(x,y)``

    - ``alpha`` -- (default: 1)

    - ``markersize`` -- (default: 50)

    - ``marker`` -- the style of the markers (default: ``'o'``); see the
      documentation of :func:`plot` for the full list of markers

    - ``facecolor`` -- (default: ``'#fec7b8'``)

    - ``edgecolor`` -- (default: ``'black'``)

    - ``zorder`` -- (default: 5)

    EXAMPLES::

        sage: scatter_plot([[0,1],[2,2],[4.3,1.1]], marker='s')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        from sage.plot.scatter_plot import ScatterPlot
        S = scatter_plot([[0,1],[2,2],[4.3,1.1]], marker='s')
        sphinx_plot(S)

    Extra options will get passed on to :meth:`~Graphics.show`, as long as they are valid::

        sage: scatter_plot([(0, 0), (1, 1)], markersize=100, facecolor='green', ymax=100)
        Graphics object consisting of 1 graphics primitive
        sage: scatter_plot([(0, 0), (1, 1)], markersize=100, facecolor='green').show(ymax=100) # These are equivalent

    .. PLOT::

        from sage.plot.scatter_plot import ScatterPlot
        S = scatter_plot([(0, 0), (1, 1)], markersize=100, facecolor='green', ymax=100)
        sphinx_plot(S)
    """
    import numpy
    from sage.plot.all import Graphics
    g = Graphics()
    g._set_extra_kwds(Graphics._extract_kwds_for_show(options))
    data = numpy.array(datalist, dtype='float')
    if len(data) != 0:
        xdata = data[:, 0]
        ydata = data[:, 1]
        g.add_primitive(ScatterPlot(xdata, ydata, options=options))
    return g
