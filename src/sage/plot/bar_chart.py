# sage_setup: distribution = sagemath-plot
"""
Bar charts
"""

# *****************************************************************************
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
# *****************************************************************************
from sage.plot.primitive import GraphicPrimitive
from sage.plot.plot import minmax_data
from sage.plot.graphics import Graphics
from sage.misc.decorators import options, rename_keyword

# TODO: make bar_chart more general


class BarChart(GraphicPrimitive):
    """
    Graphics primitive that represents a bar chart.

    EXAMPLES::

        sage: from sage.plot.bar_chart import BarChart
        sage: g = BarChart(list(range(4)), [1,3,2,0], {}); g
        BarChart defined by a 4 datalist
        sage: type(g)
        <class 'sage.plot.bar_chart.BarChart'>
    """
    def __init__(self, ind, datalist, options):
        """
        Initialize a ``BarChart`` primitive.

        EXAMPLES::

            sage: from sage.plot.bar_chart import BarChart
            sage: BarChart(list(range(3)), [10,3,5], {'width':0.7})
            BarChart defined by a 3 datalist
        """
        self.datalist = datalist
        self.ind = ind
        GraphicPrimitive.__init__(self, options)

    def get_minmax_data(self):
        """
        Return a dictionary with the bounding box data.

        EXAMPLES::

            sage: b = bar_chart([-2.3,5,-6,12])
            sage: d = b.get_minmax_data()
            sage: d['xmin']
            0
            sage: d['xmax']
            4
        """
        return minmax_data([0, len(self.datalist)], self.datalist, dict=True)

    def _allowed_options(self):
        """
        Return the allowed options with descriptions for this graphics
        primitive. This is used in displaying an error message when the
        user gives an option that doesn't make sense.

        EXAMPLES::

            sage: from sage.plot.bar_chart import BarChart
            sage: g = BarChart(list(range(4)), [1,3,2,0], {})
            sage: list(sorted(g._allowed_options().items()))
            [('hue', 'The color given as a hue.'), ('legend_label', 'The label for this item in the legend.'), ('rgbcolor', 'The color as an RGB tuple.'), ('width', 'The width of the bars'), ('zorder', 'The layer level in which to draw')]
        """
        return {'rgbcolor': 'The color as an RGB tuple.',
                'hue': 'The color given as a hue.',
                'legend_label': 'The label for this item in the legend.',
                'width': 'The width of the bars',
                'zorder': 'The layer level in which to draw'}

    def _repr_(self):
        """
        Return text representation of this bar chart graphics primitive.

        EXAMPLES::

            sage: from sage.plot.bar_chart import BarChart
            sage: g = BarChart(list(range(4)), [1,3,2,0], {})
            sage: g._repr_()
            'BarChart defined by a 4 datalist'
        """
        return "BarChart defined by a %s datalist" % (len(self.datalist))

    def _render_on_subplot(self, subplot):
        """
        Render this bar chart graphics primitive on a matplotlib subplot
        object.

        EXAMPLES:

        This rendering happens implicitly when the following command
        is executed::

            sage: bar_chart([1,2,10])
            Graphics object consisting of 1 graphics primitive
        """
        options = self.options()
        color = options['rgbcolor']
        width = float(options['width'])
        # it is critical to make NumPy arrays of type float below,
        # or bar will go boom:
        import numpy
        ind = numpy.array(self.ind, dtype=float)
        datalist = numpy.array(self.datalist, dtype=float)
        subplot.bar(ind, datalist, color=color, width=width, label=options['legend_label'])


@rename_keyword(color='rgbcolor')
@options(width=0.5, rgbcolor=(0, 0, 1), legend_label=None, aspect_ratio='automatic')
def bar_chart(datalist, **options):
    """
    A bar chart of (currently) one list of numerical data.
    Support for more data lists in progress.

    EXAMPLES:

    A bar_chart with blue bars::

        sage: bar_chart([1,2,3,4])
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(bar_chart([1,2,3,4]))

    A bar_chart with thinner bars::

        sage: bar_chart([x^2 for x in range(1,20)], width=0.2)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(bar_chart([x**2 for x in range(1,20)], width=0.2))

    A bar_chart with negative values and red bars::

        sage: bar_chart([-3,5,-6,11], rgbcolor=(1,0,0))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(bar_chart([-3,5,-6,11], rgbcolor=(1,0,0)))

    A bar chart with a legend (it's possible, not necessarily useful)::

        sage: bar_chart([-1,1,-1,1], legend_label='wave')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(bar_chart([-1,1,-1,1], legend_label='wave'))

    Extra options will get passed on to show(), as long as they are valid::

        sage: bar_chart([-2,8,-7,3], rgbcolor=(1,0,0), axes=False)
        Graphics object consisting of 1 graphics primitive
        sage: bar_chart([-2,8,-7,3], rgbcolor=(1,0,0)).show(axes=False) # These are equivalent

    .. PLOT::

        sphinx_plot(bar_chart([-2,8,-7,3], rgbcolor=(1,0,0), axes=False))
    """
    dl = len(datalist)
    if dl == 3:
        datalist = datalist + [0]
    # bardata = []
    # cnt = 1
    # for pnts in datalist:
        # ind = [i+cnt/dl for i in range(len(pnts))]
        # bardata.append([ind, pnts, xrange, yrange])
        # cnt += 1

    g = Graphics()
    g._set_extra_kwds(Graphics._extract_kwds_for_show(options))
    # TODO: improve below for multiple data sets!
    # cnt = 1
    # for ind, pnts, xrange, yrange in bardata:
    #    options={'rgbcolor':hue(cnt/dl),'width':0.5/dl}
    #    g._bar_chart(ind, pnts, xrange, yrange, options=options)
    #    cnt += 1
    # else:
    ind = list(range(len(datalist)))
    g.add_primitive(BarChart(ind, datalist, options=options))
    if options['legend_label']:
        g.legend(True)
    return g
