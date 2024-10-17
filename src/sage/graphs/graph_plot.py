# sage.doctest: optional - sage.plot
r"""
Graph plotting

*(For LaTeX drawings of graphs, see the* :mod:`~sage.graphs.graph_latex` *module.)*

All graphs have an associated Sage graphics object, which you can display::

    sage: G = graphs.WheelGraph(15)
    sage: P = G.plot()
    sage: P.show()  # long time

.. PLOT::

    sphinx_plot(graphs.WheelGraph(15))

When plotting a graph created using Sage's ``Graph`` command,
node positions are determined using the spring-layout algorithm.
Special graphs available from ``graphs.*`` have preset positions.
For example, compare the two plots of the Petersen graph,
as obtained using ``Graph`` or as obtained from that database::

    sage: petersen_spring = Graph(':I`ES@obGkqegW~')
    sage: petersen_spring.show()  # long time

.. PLOT::

    petersen_spring = Graph(':I`ES@obGkqegW~')
    sphinx_plot(petersen_spring)

::

    sage: petersen_database = graphs.PetersenGraph()
    sage: petersen_database.show()  # long time

.. PLOT::

    petersen_database = graphs.PetersenGraph()
    sphinx_plot(petersen_database)

All constructors in this database (except some random graphs) prefill
the position dictionary, bypassing the spring-layout positioning algorithm.

**Plot options**

Here is the list of options accepted by
:meth:`~sage.graphs.generic_graph.GenericGraph.plot` and the constructor of
:class:`GraphPlot`. Those two functions also accept all options of
:meth:`sage.plot.graphics.Graphics.show`.

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

{PLOT_OPTIONS_TABLE}

**Default options**

This module defines two dictionaries containing default options for the
:meth:`~sage.graphs.generic_graph.GenericGraph.plot` and
:meth:`~sage.graphs.generic_graph.GenericGraph.show` methods. These two
dictionaries are ``sage.graphs.graph_plot.DEFAULT_PLOT_OPTIONS`` and
``sage.graphs.graph_plot.DEFAULT_SHOW_OPTIONS``, respectively.

Obviously, these values are overruled when arguments are given explicitly.

Here is how to define the default size of a graph drawing to be ``(6, 6)``.
The first two calls to :meth:`~sage.graphs.generic_graph.GenericGraph.show`
use this option, while the third does not (a value for ``figsize``
is explicitly given)::

    sage: import sage.graphs.graph_plot
    sage: sage.graphs.graph_plot.DEFAULT_SHOW_OPTIONS['figsize'] = (6, 6)
    sage: graphs.PetersenGraph().show()  # long time
    sage: graphs.ChvatalGraph().show()  # long time
    sage: graphs.PetersenGraph().show(figsize=(4, 4))  # long time

We can now reset the default to its initial value, and now display graphs as
previously::

    sage: sage.graphs.graph_plot.DEFAULT_SHOW_OPTIONS['figsize'] = (4, 4)
    sage: graphs.PetersenGraph().show()  # long time
    sage: graphs.ChvatalGraph().show()  # long time

.. NOTE::

    * While ``DEFAULT_PLOT_OPTIONS`` affects both ``G.show()`` and ``G.plot()``,
      settings from ``DEFAULT_SHOW_OPTIONS`` only affects ``G.show()``.

    * In order to define a default value permanently, you can add a couple of
      lines to `Sage's startup scripts <../../../repl/startup.html>`_. Example::

       sage: import sage.graphs.graph_plot
       sage: sage.graphs.graph_plot.DEFAULT_SHOW_OPTIONS['figsize'] = (4, 4)

**Index of methods and functions**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`GraphPlot.set_pos` | Set the position plotting parameters for this GraphPlot.
    :meth:`GraphPlot.set_vertices` | Set the vertex plotting parameters for this GraphPlot.
    :meth:`GraphPlot.set_edges` | Set the edge (or arrow) plotting parameters for the GraphPlot object.
    :meth:`GraphPlot.show` | Show the (Di)Graph associated with this GraphPlot object.
    :meth:`GraphPlot.plot` | Return a graphics object representing the (di)graph.
    :meth:`GraphPlot.layout_tree` | Compute a nice layout of a tree.
"""

# ****************************************************************************
#      Copyright (C) 2009   Emily Kirkman
#                    2009   Robert L. Miller <rlmillster@gmail.com>
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

from collections import defaultdict
from math import sqrt, cos, sin, atan, pi
from sage.structure.sage_object import SageObject
from sage.misc.lazy_import import lazy_import
lazy_import("sage.plot.all", [
    "Graphics", "scatter_plot", "bezier_path", "line", "arrow", "text", "circle"])


layout_options = {
    'layout':
        'A layout algorithm -- one of : "acyclic", "circular" (plots the '
        'graph with vertices evenly distributed on a circle), "ranked", '
        '"graphviz", "planar", "spring" (traditional spring layout, using '
        'the graph\'s current positions as initial positions), or "tree" '
        '(the tree will be plotted in levels, depending on minimum distance '
        'for the root).',
    'iterations':
        'The number of times to execute the spring layout algorithm.',
    'heights':
        'A dictionary mapping heights to the list of vertices at this height.',
    'spring':
        'Use spring layout to finalize the current layout.',
    'tree_root':
        'A vertex designation for drawing trees. A vertex of the tree to '
        'be used as the root for the ``layout=\'tree\'`` option. If no root '
        'is specified, then one is chosen close to the center of the tree. '
        'Ignored unless ``layout=\'tree\'``.',
    'forest_roots':
        'An iterable specifying which vertices to use as roots for the '
        '``layout=\'forest\'`` option. If no root is specified for a tree, '
        'then one is chosen close to the center of the tree. '
        'Ignored unless ``layout=\'forest\'``.',
    'tree_orientation':
        'The direction of tree branches -- \'up\', \'down\', '
        '\'left\' or \'right\'.',
    'save_pos':
        'Whether or not to save the computed position for the graph.',
    'dim':
        'The dimension of the layout -- 2 or 3.',
    'prog':
        'Which graphviz layout program to use -- one of '
        '"circo", "dot", "fdp", "neato", or "twopi".',
    'by_component':
        'Whether to do the spring layout by connected component -- boolean.'}

graphplot_options = layout_options.copy()

graphplot_options.update({
    'pos':
        'The position dictionary of vertices.',
    'vertex_labels':
        'Vertex labels to draw. This can be ``True``/``False`` to indicate '
        'whether to print the vertex string representation of not, '
        'a dictionary keyed by vertices and associating to each vertex '
        'a label string, or a function taking as input a vertex and returning '
        'a label string.',
    'vertex_label_shift':
        'If layout is circular and we have vertex labels, will shift vertices '
        'away from center of circle in coordinate fashion `(x, y)`.',
    'vertex_color':
        'Default color for vertices not listed '
        'in vertex_colors dictionary.',
    'vertex_colors':
        'A dictionary specifying vertex colors: '
        'each key is a color recognizable by matplotlib, '
        'and each corresponding value is a list of vertices.',
    'vertex_size':
        'The size to draw the vertices.',
    'vertex_shape':
        'The shape to draw the vertices. '
        'Currently unavailable for Multi-edged DiGraphs.',
    'edge_labels':
        'Whether or not to draw edge labels.',
    'edge_style':
        'The linestyle of the edges. It should be '
        'one of "solid", "dashed", "dotted", "dashdot", '
        'or "-", "--", ":", "-.", respectively. ',
    'edge_styles':
        'A dictionary specifying edge styles: '
        'each key is an edge or a label (all same) and value is the linestyle '
        'of the edge. It should be one of "solid", "dashed", "dotted", '
        '"dashdot", or "-", "--", ":", "-.", respectively.',
    'edge_thickness':
        'The thickness of the edges.',
    'edge_thicknesses':
        'A dictionary specifying edge thicknesses: '
        'each key is an edge or a label (all same) and thickness of the '
        'corresponding edge.',
    'edge_color':
        'The default color for edges not listed in edge_colors.',
    'edge_colors':
        'A dictionary specifying edge colors: '
        'each key is a color recognized by matplotlib, '
        'and each corresponding value is a list of edges.',
    'color_by_label':
        'Whether to color the edges according to their labels. This also '
        'accepts a function or dictionary mapping labels to colors.',
    'partition':
        'A partition of the vertex set. If specified, plot will show each '
        'cell in a different color; vertex_colors takes precedence.',
    'loop_size':
        'The radius of the smallest loop.',
    'arrowsize':
        'Size of arrows.',
    'dist':
        'The distance between multiedges.',
    'max_dist':
        'The max distance range to allow multiedges.',
    'talk':
        'Whether to display the vertices in talk mode (larger and white).',
    'label_fontsize':
        'font size of all labels',
    'graph_border':
        'Whether or not to draw a frame around the graph.',
    'edge_labels_background':
        'The color of the background of the edge labels.'})

_PLOT_OPTIONS_TABLE = ""

for key, value in graphplot_options.items():
    _PLOT_OPTIONS_TABLE += f"    ``{key}`` | {value}\n"

__doc__ = __doc__.format(PLOT_OPTIONS_TABLE=_PLOT_OPTIONS_TABLE)

DEFAULT_SHOW_OPTIONS = {'figsize': (4, 4)}

DEFAULT_PLOT_OPTIONS = {
    'vertex_size'               : 200,
    'vertex_labels'             : True,
    'vertex_label_shift'        : None,
    'layout'                    : None,
    'edge_style'                : 'solid',
    'edge_styles'               : None,
    'edge_thickness'            : 1,
    'edge_thicknesses'          : None,
    'edge_color'                : 'black',
    'edge_colors'               : None,
    'edge_labels'               : False,
    'iterations'                : 50,
    'tree_orientation'          : 'down',
    'heights'                   : None,
    'graph_border'              : False,
    'talk'                      : False,
    'color_by_label'            : False,
    'partition'                 : None,
    'dist'                      : .075,
    'max_dist'                  : 1.5,
    'label_fontsize'            : 12,
    'loop_size'                 : .075,
    'edge_labels_background'    : 'white'}


class GraphPlot(SageObject):
    def __init__(self, graph, options):
        """
        Return a ``GraphPlot`` object, which stores all the parameters needed
        for plotting (Di)Graphs.

        A ``GraphPlot`` has a plot and show function, as well as some functions
        to set parameters for vertices and edges.  This constructor assumes
        default options are set.  Defaults are shown in the example below.

        EXAMPLES::

            sage: from sage.graphs.graph_plot import GraphPlot
            sage: options = {
            ....:     'vertex_size': 200,
            ....:     'vertex_labels': True,
            ....:     'layout': None,
            ....:     'edge_style': 'solid',
            ....:     'edge_color': 'black',
            ....:     'edge_colors': None,
            ....:     'edge_labels': False,
            ....:     'iterations': 50,
            ....:     'tree_orientation': 'down',
            ....:     'heights': None,
            ....:     'graph_border': False,
            ....:     'talk': False,
            ....:     'color_by_label': False,
            ....:     'partition': None,
            ....:     'dist': .075,
            ....:     'max_dist': 1.5,
            ....:     'loop_size': .075,
            ....:     'edge_labels_background': 'transparent'}
            sage: g = Graph({0: [1, 2], 2: [3], 4: [0, 1]})
            sage: GP = GraphPlot(g, options)
        """
        # Setting the default values if needed
        for k, value in DEFAULT_PLOT_OPTIONS.items():
            if k not in options:
                options[k] = value
        self._plot_components = {}
        self._nodelist = list(graph)
        self._graph = graph
        self._options = options  # contains both plot and show options
        self._arcs = self._graph.has_multiple_edges(to_undirected=True)
        self._loops = self._graph.has_loops()
        self._arcdigraph = self._graph.is_directed() and self._arcs

        self.set_pos()
        self.set_vertices()
        self.set_edges()

    def _repr_(self):
        """
        Return a string representation of a ``GraphPlot`` object.

        EXAMPLES:

        This function is called implicitly by the code below::

            sage: g = Graph({0: [1, 2], 2: [3], 4: [0, 1]})
            sage: g.graphplot()  # indirect doctest
            GraphPlot object for Graph on 5 vertices
        """
        return f"GraphPlot object for {self._graph}"

    def set_pos(self):
        """
        Set the position plotting parameters for this GraphPlot.

        EXAMPLES:

        This function is called implicitly by the code below::

            sage: g = Graph({0: [1, 2], 2: [3], 4: [0, 1]})
            sage: g.graphplot(save_pos=True, layout='circular')  # indirect doctest
            GraphPlot object for Graph on 5 vertices

        The following illustrates the format of a position dictionary, but due
        to numerical noise we do not check the values themselves::

            sage: g.get_pos()
            {0: (0.0, 1.0),
             1: (-0.951..., 0.309...),
             2: (-0.587..., -0.809...),
             3: (0.587..., -0.809...),
             4: (0.951..., 0.309...)}

        ::

            sage: T = list(graphs.trees(7))
            sage: t = T[3]
            sage: t.plot(heights={0: [0], 1: [4, 5, 1], 2: [2], 3: [3, 6]})
            Graphics object consisting of 14 graphics primitives

        .. PLOT::

            g = Graph({0: [1, 2], 2: [3], 4: [0, 1]})
            g.graphplot(save_pos=True, layout='circular')  # indirect doctest
            T = list(graphs.trees(7))
            t = T[3]
            P = t.plot(heights={0: [0], 1: [4, 5, 1], 2: [2], 3: [3, 6]})
            sphinx_plot(P)

        TESTS:

        Make sure that vertex locations are floats.  Not being floats isn't
        a bug in itself but made it too easy to accidentally introduce a bug
        elsewhere, such as in :meth:`set_edges` (:issue:`10124`), via silent
        truncating division of Python 2 integers::

            sage: g = graphs.FruchtGraph()
            sage: gp = g.graphplot()
            sage: set(map(type, flatten(gp._pos.values())))
            {<... 'float'>}
            sage: g = graphs.BullGraph()
            sage: gp = g.graphplot(save_pos=True)
            sage: set(map(type, flatten(gp._pos.values())))
            {<... 'float'>}

        Non-ascii labels are also possible using unicode (:issue:`21008`)::

            sage: Graph({u'où': [u'là', u'ici']}).plot()
            Graphics object consisting of 6 graphics primitives
        """
        self._pos = self._graph.layout(**self._options)
        # Make sure the positions are floats (trac #10124)
        self._pos = {k: (float(v[0]), float(v[1]))
                     for k, v in self._pos.items()}

    def set_vertices(self, **vertex_options):
        """
        Set the vertex plotting parameters for this ``GraphPlot``.

        This function is called by the constructor but can also be
        called to make updates to the vertex options of an existing
        ``GraphPlot`` object. Note that the changes are cumulative.

        EXAMPLES::

            sage: g = Graph({}, loops=True, multiedges=True, sparse=True)
            sage: g.add_edges([(0, 0, 'a'), (0, 0, 'b'), (0, 1, 'c'),
            ....:              (0, 1, 'd'), (0, 1, 'e'), (0, 1, 'f'),
            ....:              (0, 1, 'f'), (2, 1, 'g'), (2, 2, 'h')])
            sage: GP = g.graphplot(vertex_size=100, edge_labels=True,
            ....:                  color_by_label=True, edge_style='dashed')
            sage: GP.set_vertices(talk=True)
            sage: GP.plot()
            Graphics object consisting of 22 graphics primitives
            sage: GP.set_vertices(vertex_color='green', vertex_shape='^')
            sage: GP.plot()
            Graphics object consisting of 22 graphics primitives

        .. PLOT::

            g = Graph({}, loops=True, multiedges=True, sparse=True)
            g.add_edges([(0, 0, 'a'), (0, 0, 'b'), (0, 1, 'c'),
                         (0, 1, 'd'), (0, 1, 'e'), (0, 1, 'f'),
                         (0, 1, 'f'), (2, 1, 'g'), (2, 2, 'h')])
            GP = g.graphplot(vertex_size=100, edge_labels=True,
                             color_by_label=True, edge_style='dashed')
            GP.set_vertices(talk=True)
            sphinx_plot(GP)

        .. PLOT::

            g = Graph({}, loops=True, multiedges=True, sparse=True)
            g.add_edges([(0, 0, 'a'), (0, 0, 'b'), (0, 1, 'c'),
                         (0, 1, 'd'), (0, 1, 'e'), (0, 1, 'f'),
                         (0, 1, 'f'), (2, 1, 'g'), (2, 2, 'h')])
            GP = g.graphplot(vertex_size=100, edge_labels=True,
                             color_by_label=True, edge_style='dashed')
            GP.set_vertices(talk=True)
            GP.set_vertices(vertex_color='green', vertex_shape='^')
            sphinx_plot(GP)

        Vertex labels are flexible::

            sage: g = graphs.PathGraph(4)
            sage: g.plot(vertex_labels=False)
            Graphics object consisting of 4 graphics primitives

        .. PLOT::

            g = graphs.PathGraph(4)
            P = g.graphplot(vertex_labels=False)
            sphinx_plot(P)

        ::

            sage: g = graphs.PathGraph(4)
            sage: g.plot(vertex_labels=True)
            Graphics object consisting of 8 graphics primitives

        .. PLOT::

            g = graphs.PathGraph(4)
            P = g.graphplot(vertex_labels=True)
            sphinx_plot(P)

        ::

            sage: g = graphs.PathGraph(4)
            sage: g.plot(vertex_labels=dict(zip(g, ['+', '-', '/', '*'])))
            Graphics object consisting of 8 graphics primitives

        .. PLOT::

            g = graphs.PathGraph(4)
            P = g.graphplot(vertex_labels=dict(zip(g, ['+', '-', '/', '*'])))
            sphinx_plot(P)

        ::

            sage: g = graphs.PathGraph(4)
            sage: g.plot(vertex_labels=lambda x: str(x % 2))
            Graphics object consisting of 8 graphics primitives

        .. PLOT::

            g = graphs.PathGraph(4)
            P = g.graphplot(vertex_labels=lambda x: str(x % 2))
            sphinx_plot(P)
        """
        # Handle base vertex options
        voptions = {}

        for arg in vertex_options:
            self._options[arg] = vertex_options[arg]

        # First set defaults for styles
        vertex_colors = None
        if self._options['talk']:
            voptions['markersize'] = 500
            if self._options['partition'] is None:
                vertex_colors = '#ffffff'
        else:
            voptions['markersize'] = self._options['vertex_size']

        if ('vertex_color' not in self._options
                or self._options['vertex_color'] is None):
            vertex_color = '#fec7b8'
        else:
            vertex_color = self._options['vertex_color']

        if ('vertex_colors' not in self._options
                or self._options['vertex_colors'] is None):
            if self._options['partition'] is not None:
                from sage.plot.colors import rainbow
                partition = self._options['partition']
                length = len(partition)
                R = rainbow(length)
                vertex_colors = {R[i]: partition[i] for i in range(length)}
            elif not vertex_colors:
                vertex_colors = vertex_color
        else:
            vertex_colors = self._options['vertex_colors']

        if 'vertex_shape' in self._options:
            voptions['marker'] = self._options['vertex_shape']

        if self._graph.is_directed():
            self._vertex_radius = sqrt(voptions['markersize'] / pi)
            self._arrowshorten = 2 * self._vertex_radius
            if self._arcdigraph:
                self._vertex_radius = sqrt(voptions['markersize'] / (20500 * pi))

        voptions['zorder'] = 7

        if not isinstance(vertex_colors, dict):
            voptions['facecolor'] = vertex_colors
            pos = list(self._pos.values())
            if self._arcdigraph:
                self._plot_components['vertices'] = [
                    circle(p, self._vertex_radius, fill=True, clip=False,
                           edgecolor='black', facecolor=vertex_colors)
                    for p in pos]
            else:
                self._plot_components['vertices'] = (
                    scatter_plot(pos, clip=False, **voptions))
        else:
            # Color list must be ordered:
            pos = []
            colors = []
            for i in vertex_colors:
                pos.extend([self._pos[j] for j in vertex_colors[i]])
                colors.extend([i] * len(vertex_colors[i]))

            # If all the vertices have not been assigned a color
            if len(self._pos) != len(pos):
                leftovers = [j for j in self._pos.values() if j not in pos]
                pos.extend(leftovers)
                colors.extend([vertex_color] * len(leftovers))

            if self._arcdigraph:
                self._plot_components['vertices'] = [
                    circle(p, self._vertex_radius, fill=True, clip=False,
                           facecolor=colors[i], edgecolor='black')
                    for i, p in enumerate(pos)]
            else:
                self._plot_components['vertices'] = scatter_plot(
                    pos, facecolor=colors, clip=False, **voptions)

        vlabels = self._options['vertex_labels']
        if vlabels:
            if vlabels is True:
                vfun = str
            elif isinstance(vlabels, dict):
                def vfun(x):
                    return vlabels.get(x, "")
            else:
                vfun = vlabels

            # TODO: allow text options
            if self._options['layout'] == 'circular' and self._options['vertex_label_shift'] is not None:
                def pos_shift(v, shift):
                    return (v[0] + (v[0] * shift[0])/100, v[1] + (v[1] * shift[1])/100)
                self._plot_components['vertex_labels'] = [
                    text(
                        vfun(v),
                        pos_shift(self._pos[v], self._options['vertex_label_shift']),
                        fontsize=self._options['label_fontsize'],
                        color='black',
                        zorder=8
                    )
                    for v in self._nodelist
                    ]
            else:
                self._plot_components['vertex_labels'] = [
                    text(vfun(v), self._pos[v], color='black', zorder=8, fontsize=self._options['label_fontsize'])
                    for v in self._nodelist
                    ]

    def set_edges(self, **edge_options):
        """
        Set edge plotting parameters for the ``GraphPlot`` object.

        This function is called by the constructor but can also be called to
        update the edge options of an existing ``GraphPlot`` object.
        Note that the changes are cumulative.

        EXAMPLES::

            sage: g = Graph(loops=True, multiedges=True, sparse=True)
            sage: g.add_edges([(0, 0, 'a'), (0, 0, 'b'), (0, 1, 'c'),
            ....:              (0, 1, 'd'), (0, 1, 'e'), (0, 1, 'f'),
            ....:              (0, 1, 'f'), (2, 1, 'g'), (2, 2, 'h')])
            sage: GP = g.graphplot(vertex_size=100, edge_labels=True,
            ....:                  color_by_label=True, edge_style='dashed')
            sage: GP.set_edges(edge_style='solid')
            sage: GP.plot()
            Graphics object consisting of 22 graphics primitives

        .. PLOT::

            g = Graph(loops=True, multiedges=True, sparse=True)
            g.add_edges([(0, 0, 'a'), (0, 0, 'b'), (0, 1, 'c'),
                         (0, 1, 'd'), (0, 1, 'e'), (0, 1, 'f'),
                         (0, 1, 'f'), (2, 1, 'g'), (2, 2, 'h')])
            GP = g.graphplot(vertex_size=100, edge_labels=True,
                             color_by_label=True, edge_style='dashed')
            GP.set_edges(edge_style='solid')
            sphinx_plot(GP)

        ::

            sage: GP.set_edges(edge_color='black')
            sage: GP.plot()
            Graphics object consisting of 22 graphics primitives

        .. PLOT::

            g = Graph(loops=True, multiedges=True, sparse=True)
            g.add_edges([(0, 0, 'a'), (0, 0, 'b'), (0, 1, 'c'),
                         (0, 1, 'd'), (0, 1, 'e'), (0, 1, 'f'),
                         (0, 1, 'f'), (2, 1, 'g'), (2, 2, 'h')])
            GP = g.graphplot(vertex_size=100, edge_labels=True,
                             color_by_label=True, edge_style='dashed')
            GP.set_edges(edge_style='solid')
            GP.set_edges(edge_color='black')
            sphinx_plot(GP)

        ::

            sage: d = DiGraph(loops=True, multiedges=True, sparse=True)
            sage: d.add_edges([(0, 0, 'a'), (0, 0, 'b'), (0, 1, 'c'),
            ....:              (0, 1, 'd'), (0, 1, 'e'), (0, 1, 'f'),
            ....:              (0, 1, 'f'), (2, 1, 'g'), (2, 2, 'h')])
            sage: GP = d.graphplot(vertex_size=100, edge_labels=True,
            ....:                  color_by_label=True, edge_style='dashed')
            sage: GP.set_edges(edge_style='solid')
            sage: GP.plot()
            Graphics object consisting of 24 graphics primitives

        .. PLOT::

            d = DiGraph(loops=True, multiedges=True, sparse=True)
            d.add_edges([(0, 0, 'a'), (0, 0, 'b'), (0, 1, 'c'),
                         (0, 1, 'd'), (0, 1, 'e'), (0, 1, 'f'),
                         (0, 1, 'f'), (2, 1, 'g'), (2, 2, 'h')])
            GP = d.graphplot(vertex_size=100, edge_labels=True,
                             color_by_label=True, edge_style='dashed')
            GP.set_edges(edge_style='solid')
            sphinx_plot(GP)

        ::

            sage: GP.set_edges(edge_color='black')
            sage: GP.plot()
            Graphics object consisting of 24 graphics primitives

        .. PLOT::

            d = DiGraph(loops=True, multiedges=True, sparse=True)
            d.add_edges([(0, 0, 'a'), (0, 0, 'b'), (0, 1, 'c'),
                         (0, 1, 'd'), (0, 1, 'e'), (0, 1, 'f'),
                         (0, 1, 'f'), (2, 1, 'g'), (2, 2, 'h')])
            GP = d.graphplot(vertex_size=100, edge_labels=True,
                             color_by_label=True, edge_style='dashed')
            GP.set_edges(edge_style='solid')
            GP.set_edges(edge_color='black')
            sphinx_plot(GP)

        TESTS::

            sage: G = Graph("Fooba")
            sage: G.show(edge_colors={'red':[(3, 6), (2, 5)]})

        Check default edge labels are pretty close to halfway between
        the vertices in some cases where they weren't due to Python 2
        truncating division (:issue:`10124`)::

            sage: test_graphs = graphs.FruchtGraph(), graphs.BullGraph()
            sage: tol = 0.001
            sage: for G in test_graphs:
            ....:     E = G.edges(sort=True)
            ....:     for e0, e1, elab in E:
            ....:         G.set_edge_label(e0, e1, '%d %d' % (e0, e1))
            ....:     gp = G.graphplot(save_pos=True, edge_labels=True)
            ....:     vx = gp._plot_components['vertices'][0].xdata
            ....:     vy = gp._plot_components['vertices'][0].ydata
            ....:     for elab in gp._plot_components['edge_labels']:
            ....:         textobj = elab[0]
            ....:         x, y, s = textobj.x, textobj.y, textobj.string
            ....:         v0, v1 = map(int, s.split())
            ....:         m = sum(vector((vx[v], vy[v])) for v in (v0, v1))/2
            ....:         assert (vector((x, y)) - m).norm() < tol

        Issue :issue:`24051` is fixed::

            sage: G = Graph([(0, 1), (0, 1)], multiedges=True)
            sage: G.plot(edge_colors={"red": [(1, 0)]})
            Graphics object consisting of 5 graphics primitives

        Issue :issue:`31542` is fixed::

            sage: s = 'ABCCCCDABCDABCDA'
            sage: g = DiGraph({}, loops=True, multiedges=True)
            sage: for a, b in [(s[i], u) for i, u in enumerate(s[1:])]:
            ....:     g.add_edge(a, b, b)
            sage: g.plot(color_by_label=True, layout='circular')
            Graphics object consisting of 23 graphics primitives
        """
        for arg in edge_options:
            self._options[arg] = edge_options[arg]
        if 'edge_colors' in edge_options:
            self._options['color_by_label'] = False
        if self._options['edge_labels_background'] == "transparent":
            self._options['edge_labels_background'] = "None"

        # Whether a key is an edge or not:
        #   None => edge_x is not set
        #   True => keys are edges
        #   False => keys are labels
        style_key_edges = None
        thickness_key_edges = None
        if isinstance(self._options['edge_styles'], dict):
            for k in self._options['edge_styles']:
                style_key_edges = k in self._graph.edges()
                break
        if isinstance(self._options['edge_thicknesses'], dict):
            for k in self._options['edge_thicknesses']:
                thickness_key_edges = k in self._graph.edges()
                break

        eoptions = {}
        if 'arrowsize' in self._options:
            eoptions['arrowsize'] = self._options['arrowsize']

        # Set labels param to add labels on the fly
        labels = False
        if self._options['edge_labels']:
            labels = True
            self._plot_components['edge_labels'] = []

        # Make dict collection of all edges (keep label and edge color)
        edges_to_draw = defaultdict(list)

        v_to_int = {v: i for i, v in enumerate(self._graph)}

        if (self._options['color_by_label']
                or isinstance(self._options['edge_colors'], dict)):
            if self._options['color_by_label']:
                edge_colors = self._graph._color_by_label(
                    format=self._options['color_by_label'])
            else:
                edge_colors = self._options['edge_colors']
            edges_drawn = []
            for color in edge_colors:
                for edge in edge_colors[color]:
                    a, b = edge[0], edge[1]
                    if v_to_int[a] < v_to_int[b]:
                        key = (a, b)
                        head = 1
                    else:
                        key = (b, a)
                        head = 0
                    if len(edge) < 3:
                        label = self._graph.edge_label(a, b)
                        if isinstance(label, list):
                            edges_to_draw[key].append((label[-1], color, head))
                            edges_drawn.append((a, b, label[-1]))
                            for lab in label[:-1]:
                                edges_to_draw[key].append((lab, color, head))
                                edges_drawn.append((a, b, lab))
                        else:
                            edges_to_draw[key].append((label, color, head))
                            edges_drawn.append((a, b, label))
                    else:
                        label = edge[2]
                        edges_to_draw[key].append((label, color, head))
                        edges_drawn.append((a, b, label))

            # Add unspecified edges (default color black set in DEFAULT_PLOT_OPTIONS)
            for a, b, c in self._graph.edge_iterator():
                if ((a, b, c) not in edges_drawn
                        and (self._graph.is_directed()
                             or (b, a, c) not in edges_drawn)):
                    if v_to_int[a] < v_to_int[b]:
                        key = (a, b)
                        head = 1
                    else:
                        key = (b, a)
                        head = 0
                    edges_to_draw[key].append((c, self._options['edge_color'], head))

        else:
            for a, b, c in self._graph.edge_iterator():
                if v_to_int[a] < v_to_int[b]:
                    key = (a, b)
                    head = 1
                else:
                    key = (b, a)
                    head = 0
                edges_to_draw[key].append((c, self._options['edge_color'], head))

        if edges_to_draw:
            self._plot_components['edges'] = []
        else:
            return

        # Check for multi-edges or loops
        if self._arcs or self._loops:
            tmp = edges_to_draw.copy()
            dist = self._options['dist'] * 2
            min_loop_size = self._options['loop_size']
            max_dist = self._options['max_dist']
            from sage.misc.functional import sqrt
            for a, b in tmp:
                if a == b:
                    # Multiple loops need varying loop radius starting at
                    # minimum loop size and respecting other distances
                    loop_size = min_loop_size  # current loop radius
                    distance = dist
                    local_labels = edges_to_draw.pop((a, b))
                    len_local_labels = len(local_labels)
                    if len_local_labels * dist > max_dist:
                        distance = float(max_dist) / len_local_labels
                    loop_size_increment = distance / 4
                    # Now add all the loops at this vertex, varying their size
                    for lab, col, _ in local_labels:
                        x, y = self._pos[a][0], self._pos[a][1] - loop_size

                        estyle = self._options['edge_style']
                        ethickness = self._options['edge_thickness']
                        if style_key_edges is not None and ((style_key_edges and (x, y) in self._options['edge_styles']) or (not style_key_edges and lab in self._options['edge_styles'])):
                            estyle = style_key_edges and self._options['edge_styles'][(x, y)] or self._options['edge_styles'][lab]
                        if thickness_key_edges is not None and ((thickness_key_edges and (x, y) in self._options['edge_thicknesses']) or (not thickness_key_edges and lab in self._options['edge_thicknesses'])):
                            ethickness = thickness_key_edges and self._options['edge_thicknesses'][(x, y)] or self._options['edge_thicknesses'][lab]
                            
                        c = circle((x, y), loop_size, rgbcolor=col, linestyle=estyle, thickness=ethickness)
                        self._plot_components['edges'].append(c)
                        if labels:
                            bg = self._options['edge_labels_background']
                            y -= loop_size  # place label at bottom of loop
                            t = text(lab, (x, y), background_color=bg, fontsize=self._options['label_fontsize'])
                            self._plot_components['edge_labels'].append(t)
                        loop_size += loop_size_increment
                elif len(edges_to_draw[a, b]) > 1:
                    # Multi-edge
                    local_labels = edges_to_draw.pop((a, b))

                    # Compute perpendicular bisector
                    p1 = self._pos[a]
                    p2 = self._pos[b]
                    m = ((p1[0] + p2[0]) / 2., (p1[1] + p2[1]) / 2.)  # midpoint
                    if not p1[1] == p2[1]:
                        s = (p1[0] - p2[0]) / (p2[1] - p1[1])  # perp slope

                        def y(x):
                            return s * (x - m[0]) + m[1]  # perp bisector line

                        # f, g are functions to determine x-values of point
                        # on line y at distance d from point m (on each side)
                        def f(d):
                            return sqrt(d**2 / (1. + s**2)) + m[0]

                        def g(d):
                            return -sqrt(d**2 / (1. + s**2)) + m[0]

                        odd_x = f
                        even_x = g
                        if p1[0] == p2[0]:
                            def odd_y(d):
                                return m[1]

                            even_y = odd_y
                        else:
                            def odd_y(x):
                                return y(f(x))

                            def even_y(x):
                                return y(g(x))
                    else:
                        def odd_x(d):
                            return m[0]

                        even_x = odd_x

                        def odd_y(d):
                            return m[1] + d

                        def even_y(d):
                            return m[1] - d

                    def odd_xy(d):
                        return (odd_x(d), odd_y(d))

                    def even_xy(d):
                        return (even_x(d), even_y(d))

                    # We now have the control points for each Bezier curve
                    # in terms of distance parameter d.
                    # Also note each edge label should be drawn at d/2.
                    # (This is because we're using the perp bisectors).
                    distance = dist
                    len_local_labels = len(local_labels)
                    if len_local_labels * dist > max_dist:
                        distance = float(max_dist) / len_local_labels
                    for i in range(len_local_labels // 2):
                        k = (i + 1.0) * distance
                        estyle = self._options['edge_style']
                        ethickness = self._options['edge_thickness']

                        if self._arcdigraph:
                            vr = self._vertex_radius
                            ph = self._polar_hack_for_multidigraph
                            odd_start = ph(p1, odd_xy(k), vr)[0]
                            odd_end = ph(odd_xy(k), p2, vr)[1]
                            even_start = ph(p1, even_xy(k), vr)[0]
                            even_end = ph(even_xy(k), p2, vr)[1]

                            self._plot_components['edges'].append(
                                arrow(path=[[odd_start, odd_xy(k), odd_end]],
                                      head=local_labels[2 * i][2], zorder=1,
                                      rgbcolor=local_labels[2 * i][1],
                                      linestyle=estyle,
                                      width=ethickness,
                                      **eoptions
                                      ))
                            self._plot_components['edges'].append(
                                arrow(path=[[even_start, even_xy(k), even_end]],
                                      head=local_labels[2 * i + 1][2], zorder=1,
                                      rgbcolor=local_labels[2 * i + 1][1],
                                      linestyle=estyle,
                                      width=ethickness,
                                      **eoptions
                                      ))
                        else:
                            self._plot_components['edges'].append(
                                bezier_path([[p1, odd_xy(k), p2]], zorder=1,
                                            rgbcolor=local_labels[2 * i][1],
                                            linestyle=estyle,
                                            thickness=ethickness
                                            ))
                            self._plot_components['edges'].append(
                                bezier_path([[p1, even_xy(k), p2]], zorder=1,
                                            rgbcolor=local_labels[2 * i + 1][1],
                                            linestyle=estyle,
                                            thickness=ethickness
                                            ))
                        if labels:
                            j = k / 2.0
                            bg = self._options['edge_labels_background']
                            self._plot_components['edge_labels'].append(
                                text(local_labels[2 * i][0], odd_xy(j),
                                     background_color=bg, fontsize=self._options['label_fontsize']))
                            self._plot_components['edge_labels'].append(
                                text(local_labels[2 * i + 1][0], even_xy(j),
                                     background_color=bg, fontsize=self._options['label_fontsize']))
                    if len_local_labels % 2:
                        # draw line for last odd
                        edges_to_draw[a, b] = [local_labels[-1]]

        is_directed = self._graph.is_directed()
        for a, b in edges_to_draw:
            elabel = edges_to_draw[a, b][0][0]
            ecolor = edges_to_draw[a, b][0][1]
            ehead = edges_to_draw[a, b][0][2]
            e = (a, b, elabel)

            estyle = self._options['edge_style']
            ethickness = self._options['edge_thickness']
            if style_key_edges is not None and ((style_key_edges and e in self._options['edge_styles']) or (not style_key_edges and elabel in self._options['edge_styles'])):
                estyle = style_key_edges and self._options['edge_styles'][e] or self._options['edge_styles'][elabel]
            if thickness_key_edges is not None and ((thickness_key_edges and e in self._options['edge_thicknesses']) or (not thickness_key_edges and elabel in self._options['edge_thicknesses'])):
                ethickness = thickness_key_edges and self._options['edge_thicknesses'][e] or self._options['edge_thicknesses'][elabel]

            if self._arcdigraph:
                ph = self._polar_hack_for_multidigraph
                C, D = ph(self._pos[a], self._pos[b], self._vertex_radius)
                self._plot_components['edges'].append(
                    arrow(C, D,
                          rgbcolor=ecolor,
                          head=ehead,
                          linestyle=estyle,
                          width=ethickness,
                          **eoptions
                          ))
                if labels:
                    bg = self._options['edge_labels_background']
                    self._plot_components['edge_labels'].append(
                        text(str(elabel),
                             [(C[0] + D[0]) / 2., (C[1] + D[1]) / 2.],
                             background_color=bg,
                             fontsize=self._options['label_fontsize']))
            elif is_directed:
                self._plot_components['edges'].append(
                    arrow(self._pos[a], self._pos[b],
                          rgbcolor=ecolor,
                          arrowshorten=self._arrowshorten,
                          head=ehead,
                          linestyle=estyle,
                          width=ethickness,
                          **eoptions
                          ))
            else:
                self._plot_components['edges'].append(
                    line([self._pos[a], self._pos[b]],
                         rgbcolor=ecolor,
                         linestyle=estyle,
                         thickness=ethickness))
            if labels and not self._arcdigraph:
                bg = self._options['edge_labels_background']
                self._plot_components['edge_labels'].append(
                    text(str(edges_to_draw[a, b][0][0]),
                         [(self._pos[a][0] + self._pos[b][0]) / 2.,
                         (self._pos[a][1] + self._pos[b][1]) / 2.],
                         background_color=bg,
                         fontsize=self._options['label_fontsize']))

    def _polar_hack_for_multidigraph(self, A, B, VR):
        """
        Helper function to quickly compute the two points of intersection
        of a line segment from ``A`` to ``B`` (provided as xy pairs) and
        circles centered at ``A`` and ``B``, both with radius ``VR``.
        Returns a pair of xy pairs representing the two points.

        EXAMPLES::

            sage: d = DiGraph(loops=True, multiedges=True, sparse=True)
            sage: d.add_edges([(0, 0, 'a'), (0, 0, 'b'), (0, 1, 'c'),
            ....:              (0, 1, 'd'), (0, 1, 'e'), (0, 1, 'f'),
            ....:              (0, 1, 'f'), (2, 1, 'g'), (2, 2, 'h')])
            sage: GP = d.graphplot(vertex_size=100, edge_labels=True,
            ....:                  color_by_label=True, edge_style='dashed')
            sage: GP._polar_hack_for_multidigraph((0, 1), (1, 1), .1)
            ([0.10..., 1.00...], [0.90..., 1.00...])

        TESTS:

        Make sure that Python ints are acceptable arguments (:issue:`10124`)::

            sage: GP = DiGraph().graphplot()
            sage: GP._polar_hack_for_multidigraph((0, 1), (2, 2), .1)
            ([0.08..., 1.04...], [1.91..., 1.95...])
            sage: GP._polar_hack_for_multidigraph((int(0), int(1)),
            ....:                                 (int(2), int(2)), .1)
            ([0.08..., 1.04...], [1.91..., 1.95...])
        """
        D = [float(B[i] - A[i]) for i in range(2)]
        R = sqrt(D[0]**2 + D[1]**2)
        theta = 3 * pi / 2
        if D[0] > 0:
            theta = atan(D[1] / D[0])
            if D[1] < 0:
                theta += 2 * pi
        elif D[0] < 0:
            theta = atan(D[1] / D[0]) + pi
        elif D[1] > 0:
            theta = pi / 2
        cos_theta = cos(theta)
        sin_theta = sin(theta)
        return ([VR * cos_theta + A[0], VR * sin_theta + A[1]],
                [(R - VR) * cos_theta + A[0], (R - VR) * sin_theta + A[1]])

    def show(self, **kwds):
        """
        Show the (di)graph associated with this ``GraphPlot`` object.

        INPUT:

        This method accepts all parameters of
        :meth:`sage.plot.graphics.Graphics.show`.

        .. NOTE::

            - See :mod:`the module's documentation <sage.graphs.graph_plot>`
              for information on default values of this method.

            - Any options not used by plot will be passed on to the
              :meth:`~sage.plot.graphics.Graphics.show` method.

        EXAMPLES::

            sage: C = graphs.CubeGraph(8)
            sage: P = C.graphplot(vertex_labels=False, vertex_size=0,
            ....:                 graph_border=True)
            sage: P.show()

        .. PLOT::

            C = graphs.CubeGraph(8)
            P = C.graphplot(vertex_labels=False, vertex_size=0,
                            graph_border=True)
            sphinx_plot(P)
        """
        # Setting the default values if needed
        for k, value in DEFAULT_SHOW_OPTIONS.items():
            if k not in kwds:
                kwds[k] = value

        self.plot().show(**kwds)

    def plot(self, **kwds):
        """
        Return a graphics object representing the (di)graph.

        INPUT:

        The options accepted by this method are to be found in the
        documentation of the :mod:`sage.graphs.graph_plot` module,
        and the :meth:`~sage.plot.graphics.Graphics.show` method.

        .. NOTE::

            See :mod:`the module's documentation <sage.graphs.graph_plot>` for
            information on default values of this method.

        We can specify some pretty precise plotting of familiar graphs::

            sage: from math import sin, cos, pi
            sage: P = graphs.PetersenGraph()
            sage: d = {'#FF0000': [0, 5], '#FF9900': [1, 6], '#FFFF00': [2, 7],
            ....:      '#00FF00': [3, 8], '#0000FF': [4,9]}
            sage: pos_dict = {}
            sage: for i in range(5):
            ....:  x = float(cos(pi/2 + ((2*pi)/5)*i))
            ....:  y = float(sin(pi/2 + ((2*pi)/5)*i))
            ....:  pos_dict[i] = [x,y]
            ...
            sage: for i in range(5, 10):
            ....:  x = float(0.5*cos(pi/2 + ((2*pi)/5)*i))
            ....:  y = float(0.5*sin(pi/2 + ((2*pi)/5)*i))
            ....:  pos_dict[i] = [x,y]
            ...
            sage: pl = P.graphplot(pos=pos_dict, vertex_colors=d)
            sage: pl.show()

        .. PLOT::

            from math import sin, cos, pi
            P = graphs.PetersenGraph()
            d = {'#FF0000': [0, 5], '#FF9900': [1, 6], '#FFFF00': [2, 7],
                 '#00FF00': [3, 8], '#0000FF': [4,9]}
            pos_dict = {}
            for i in range(5):
                x = float(cos(pi/2 + ((2*pi)/5)*i))
                y = float(sin(pi/2 + ((2*pi)/5)*i))
                pos_dict[i] = [x,y]

            for i in range(5, 10):
                x = float(0.5*cos(pi/2 + ((2*pi)/5)*i))
                y = float(0.5*sin(pi/2 + ((2*pi)/5)*i))
                pos_dict[i] = [x,y]

            pl = P.graphplot(pos=pos_dict, vertex_colors=d)
            sphinx_plot(pl)

        Here are some more common graphs with typical options::

            sage: C = graphs.CubeGraph(8)
            sage: P = C.graphplot(vertex_labels=False, vertex_size=0,
            ....:                 graph_border=True)
            sage: P.show()

        .. PLOT::

            C = graphs.CubeGraph(8)
            P = C.graphplot(vertex_labels=False, vertex_size=0,
                            graph_border=True)
            sphinx_plot(P)

        ::

            sage: G = graphs.HeawoodGraph().copy(sparse=True)
            sage: for u, v, l in G.edges(sort=True):
            ....:     G.set_edge_label(u, v, f'({u},{v})')
            sage: G.graphplot(edge_labels=True).show()

        .. PLOT::

            G = graphs.HeawoodGraph().copy(sparse=True)
            for u, v, l in G.edges(sort=True):
                G.set_edge_label(u, v, f'({u},{v})')
            sphinx_plot(G.graphplot(edge_labels=True))

        The options for plotting also work with directed graphs::

            sage: D = DiGraph({
            ....:     0: [1, 10, 19], 1: [8, 2], 2: [3, 6], 3: [19, 4],
            ....:     4: [17, 5], 5: [6, 15], 6: [7], 7: [8, 14], 8: [9],
            ....:     9: [10, 13], 10: [11], 11: [12, 18], 12: [16, 13],
            ....:     13: [14], 14: [15], 15: [16], 16: [17], 17: [18],
            ....:     18: [19], 19: []})
            sage: for u, v, l in D.edges(sort=True):
            ....:     D.set_edge_label(u, v, f'({u},{v})')
            sage: D.graphplot(edge_labels=True, layout='circular').show()

        .. PLOT::

            D = DiGraph({
                0: [1, 10, 19], 1: [8, 2], 2: [3, 6], 3: [19, 4],
                4: [17, 5], 5: [6, 15], 6: [7], 7: [8, 14], 8: [9],
                9: [10, 13], 10: [11], 11: [12, 18], 12: [16, 13],
                13: [14], 14: [15], 15: [16], 16: [17], 17: [18],
                18: [19], 19: []})
            for u, v, l in D.edges(sort=True):
                D.set_edge_label(u, v, f'({u},{v})')
            sphinx_plot(D.graphplot(edge_labels=True, layout='circular'))
        
        For graphs with ``circular`` layouts, one may shift the vertex labels by
        specifying coordinates to shift by::

            sage: D = DiGraph({
            ....:     0: [1, 10, 19], 1: [8, 2], 2: [3, 6], 3: [19, 4],
            ....:     4: [17, 5], 5: [6, 15], 6: [7], 7: [8, 14], 8: [9],
            ....:     9: [10, 13], 10: [11], 11: [12, 18], 12: [16, 13],
            ....:     13: [14], 14: [15], 15: [16], 16: [17], 17: [18],
            ....:     18: [19], 19: []})
            sage: for u, v, l in D.edges(sort=True):
            ....:     D.set_edge_label(u, v, f'({u},{v})')
            sage: D.graphplot(edge_labels=True, layout='circular', vertex_label_shift=(15,10)).show()

        .. PLOT::

            D = DiGraph({
                0: [1, 10, 19], 1: [8, 2], 2: [3, 6], 3: [19, 4],
                4: [17, 5], 5: [6, 15], 6: [7], 7: [8, 14], 8: [9],
                9: [10, 13], 10: [11], 11: [12, 18], 12: [16, 13],
                13: [14], 14: [15], 15: [16], 16: [17], 17: [18],
                18: [19], 19: []})
            for u, v, l in D.edges(sort=True):
                D.set_edge_label(u, v, f'({u},{v})')
            sphinx_plot(D.graphplot(edge_labels=True, layout='circular', vertex_label_shift=(15,10)))

        This example shows off the coloring of edges::

            sage: from sage.plot.colors import rainbow
            sage: C = graphs.CubeGraph(5)
            sage: R = rainbow(5)
            sage: edge_colors = {}
            sage: for i in range(5):
            ....:     edge_colors[R[i]] = []
            sage: for u, v, l in C.edges(sort=True):
            ....:     for i in range(5):
            ....:         if u[i] != v[i]:
            ....:             edge_colors[R[i]].append((u, v, l))
            sage: C.graphplot(vertex_labels=False, vertex_size=0,
            ....:             edge_colors=edge_colors).show()

        .. PLOT::

            from sage.plot.colors import rainbow
            C = graphs.CubeGraph(5)
            R = rainbow(5)
            edge_colors = {}
            for i in range(5):
                edge_colors[R[i]] = []
            for u, v, l in C.edges(sort=True):
                for i in range(5):
                    if u[i] != v[i]:
                        edge_colors[R[i]].append((u, v, l))
            sphinx_plot(C.graphplot(vertex_labels=False, vertex_size=0,
                                    edge_colors=edge_colors))

        With the ``partition`` option, we can separate out same-color groups
        of vertices::

            sage: D = graphs.DodecahedralGraph()
            sage: Pi = [[6, 5, 15, 14, 7], [16, 13, 8, 2, 4],
            ....:       [12, 17, 9, 3, 1], [0, 19, 18, 10, 11]]
            sage: D.show(partition=Pi)

        .. PLOT::

            D = graphs.DodecahedralGraph()
            Pi = [[6, 5, 15, 14, 7], [16, 13, 8, 2, 4],
                  [12, 17, 9, 3, 1], [0, 19, 18, 10, 11]]
            sphinx_plot(D.plot(partition=Pi))

        Loops are also plotted correctly::

            sage: G = graphs.PetersenGraph()
            sage: G.allow_loops(True)
            sage: G.add_edge(0,0)
            sage: G.show()

        .. PLOT::

            G = graphs.PetersenGraph()
            G.allow_loops(True)
            G.add_edge(0,0)
            sphinx_plot(G)

        ::

            sage: D = DiGraph({0:[0,1], 1:[2], 2:[3]}, loops=True)
            sage: D.show()
            sage: D.show(edge_colors={(0, 1, 0): [(0, 1, None), (1, 2, None)],
            ....:                     (0, 0, 0): [(2, 3, None)]})

        .. PLOT::

            D = DiGraph({0:[0,1], 1:[2], 2:[3]}, loops=True)
            P = D.plot(edge_colors={(0, 1, 0): [(0, 1, None), (1, 2, None)],
                                    (0, 0, 0): [(2, 3, None)]})
            sphinx_plot(P)

        More options::

            sage: pos = {0: [0.0, 1.5], 1: [-0.8, 0.3], 2: [-0.6, -0.8],
            ....:        3:[0.6, -0.8], 4:[0.8, 0.3]}
            sage: g = Graph({0: [1], 1: [2], 2: [3], 3: [4], 4: [0]})
            sage: g.graphplot(pos=pos, layout='spring', iterations=0).plot()
            Graphics object consisting of 11 graphics primitives

        .. PLOT::

            pos = {0: [0.0, 1.5], 1: [-0.8, 0.3], 2: [-0.6, -0.8],
                   3: [0.6, -0.8], 4:[0.8, 0.3]}
            g = Graph({0: [1], 1: [2], 2: [3], 3: [4], 4: [0]})
            P = g.graphplot(pos=pos, layout='spring', iterations=0).plot()
            sphinx_plot(P)

        ::

            sage: D = graphs.CubeGraph(3)
            sage: D.graphplot(layout='planar').plot()
            Launched png viewer for Graphics object consisting of 21 graphics primitives
        
        .. PLOT::

            D = graphs.CubeGraph(3)
            sphinx_plot(D.graphplot(layout='planar'))

        ::

            sage: G = Graph()
            sage: P = G.graphplot().plot()
            sage: P.axes()
            False
            sage: G = DiGraph()
            sage: P = G.graphplot().plot()
            sage: P.axes()
            False

        We can plot multiple graphs::

            sage: T = list(graphs.trees(7))
            sage: t = T[3]
            sage: t.graphplot(heights={0: [0], 1: [4, 5, 1],
            ....:                      2: [2], 3: [3, 6]}
            ....:            ).plot()
            Graphics object consisting of 14 graphics primitives

        .. PLOT::

            T = list(graphs.trees(7))
            t = T[3]
            sphinx_plot(t.graphplot(heights={0: [0], 1: [4, 5, 1],
                                             2: [2], 3: [3, 6]}))

        ::

            sage: T = list(graphs.trees(7))
            sage: t = T[3]
            sage: t.graphplot(heights={0: [0], 1: [4, 5, 1],
            ....:                      2: [2], 3: [3, 6]}
            ....:            ).plot()
            Graphics object consisting of 14 graphics primitives

        .. PLOT::

            T = list(graphs.trees(7))
            t = T[3]
            sphinx_plot(t.graphplot(heights={0: [0], 1: [4, 5, 1],
                                             2: [2], 3: [3, 6]}))

        ::

            sage: t.set_edge_label(0, 1, -7)
            sage: t.set_edge_label(0, 5, 3)
            sage: t.set_edge_label(0, 5, 99)
            sage: t.set_edge_label(1, 2, 1000)
            sage: t.set_edge_label(3, 2, 'spam')
            sage: t.set_edge_label(2, 6, 3/2)
            sage: t.set_edge_label(0, 4, 66)
            sage: t.graphplot(heights={0: [0], 1: [4, 5, 1],
            ....:                      2: [2], 3: [3, 6]},
            ....:             edge_labels=True
            ....:            ).plot()
            Graphics object consisting of 20 graphics primitives

        .. PLOT::

            T = list(graphs.trees(7))
            t = T[3]
            t.set_edge_label(0, 1, -7)
            t.set_edge_label(0, 5, 3)
            t.set_edge_label(0, 5, 99)
            t.set_edge_label(1, 2, 1000)
            t.set_edge_label(3, 2, 'spam')
            t.set_edge_label(2, 6, 3/2)
            t.set_edge_label(0, 4, 66)
            sphinx_plot(t.graphplot(heights={0: [0], 1: [4, 5, 1],
                                             2: [2], 3: [3, 6]},
                                    edge_labels=True))

        ::

            sage: T = list(graphs.trees(7))
            sage: t = T[3]
            sage: t.graphplot(layout='tree').show()

        .. PLOT::

            T = list(graphs.trees(7))
            t = T[3]
            sphinx_plot(t.graphplot(layout='tree'))

        The tree layout is also useful::

            sage: t = DiGraph('JCC???@A??GO??CO??GO??')
            sage: t.graphplot(layout='tree', tree_root=0,
            ....:             tree_orientation="up"
            ....:            ).show()

        .. PLOT::

            t = DiGraph('JCC???@A??GO??CO??GO??')
            sphinx_plot(t.graphplot(layout='tree', tree_root=0,
                                    tree_orientation='up'))

        More examples::

            sage: D = DiGraph({0:[1,2,3], 2:[1,4], 3:[0]})
            sage: D.graphplot().show()

        .. PLOT::

            D = DiGraph({0:[1,2,3], 2:[1,4], 3:[0]})
            sphinx_plot(D.graphplot())

        ::

            sage: D = DiGraph({0:[1,2,3], 2:[1,4], 3:[0]})
            sage: D.graphplot(label_fontsize=20).show()
            Graphics object consisting of 8 graphics primitives

        .. PLOT::

            D = DiGraph({0:[1,2,3], 2:[1,4], 3:[0]})
            sphinx_plot(D.graphplot(label_fontsize=20))


        ::

            sage: D = DiGraph(multiedges=True, sparse=True)
            sage: for i in range(5):
            ....:   D.add_edge((i, i + 1, 'a'))
            ....:   D.add_edge((i, i - 1, 'b'))
            sage: D.graphplot(edge_labels=True,
            ....:             edge_colors=D._color_by_label()
            ....:            ).plot()
            Graphics object consisting of 34 graphics primitives

        .. PLOT::

            D = DiGraph(multiedges=True, sparse=True)
            for i in range(5):
                 D.add_edge((i, i + 1, 'a'))
                 D.add_edge((i, i - 1, 'b'))
            sphinx_plot(D.graphplot(edge_labels=True,
                                    edge_colors=D._color_by_label()))

        ::

            sage: g = Graph({}, loops=True, multiedges=True, sparse=True)
            sage: g.add_edges([(0, 0, 'a'), (0, 0, 'b'), (0, 1, 'c'),
            ....:              (0, 1, 'd'), (0, 1, 'e'), (0, 1, 'f'),
            ....:              (0, 1, 'f'), (2, 1, 'g'), (2, 2, 'h')])
            sage: g.graphplot(edge_labels=True,
            ....:             color_by_label=True,
            ....:             edge_style='dashed'
            ....:            ).plot()
            Graphics object consisting of 22 graphics primitives

        .. PLOT::

            g = Graph({}, loops=True, multiedges=True, sparse=True)
            g.add_edges([(0, 0, 'a'), (0, 0, 'b'), (0, 1, 'c'),
                         (0, 1, 'd'), (0, 1, 'e'), (0, 1, 'f'),
                         (0, 1, 'f'), (2, 1, 'g'), (2, 2, 'h')])
            sphinx_plot(g.graphplot(edge_labels=True,
                                    color_by_label=True,
                                    edge_style='dashed'))

        The ``edge_style`` option may be provided in the short format too::

            sage: g.graphplot(edge_labels=True,
            ....:             color_by_label=True,
            ....:             edge_style='--'
            ....:            ).plot()
            Graphics object consisting of 22 graphics primitives

        The ``edge_styles`` option may be provided if you need only certain edges
        to have certain styles::

            sage: GP.set_edges(edge_styles={'a':'dashed', 'g':'dotted'})
            sage: GP.plot()
            Graphics object consisting of 22 graphics primitives

        .. PLOT::

            g = Graph(loops=True, multiedges=True, sparse=True)
            g.add_edges([(0, 0, 'a'), (0, 0, 'b'), (0, 1, 'c'),
                         (0, 1, 'd'), (0, 1, 'e'), (0, 1, 'f'),
                         (0, 1, 'f'), (2, 1, 'g'), (2, 2, 'h')])
            GP = g.graphplot(vertex_size=100, edge_labels=True,
                             color_by_label=True, edge_style='dashed')
            GP.set_edges(edge_style='solid')
            GP.set_edges(edge_color='black')
            GP.set_edges(edge_styles={'a':'dashed', 'g':'dotted'})
            sphinx_plot(GP)

        TESTS:

        Make sure that show options work with plot also::

            sage: g = Graph({})
            sage: g.plot(title='empty graph', axes=True)
            Graphics object consisting of 0 graphics primitives

        Check for invalid inputs::

            sage: p = graphs.PetersenGraph().plot(egabrag='garbage')
            Traceback (most recent call last):
            ...
            ValueError: invalid input 'egabrag=garbage'

        Make sure that no graphics primitive is clipped::

            sage: tadpole = Graph({0: [0, 1]}).plot()
            sage: bbox = tadpole.get_minmax_data()
            sage: for part in tadpole:
            ....:      part_bbox = part.get_minmax_data()
            ....:      assert (bbox['xmin'] <= part_bbox['xmin']
            ....:              <= part_bbox['xmax'] <= bbox['xmax'])
            ....:      assert (bbox['ymin'] <= part_bbox['ymin']
            ....:              <= part_bbox['ymax'] <= bbox['ymax'])

        Check that one can plot immutable graphs (:issue:`17340`)::

            sage: Graph({0: [0]}, immutable=True).plot()
            Graphics object consisting of 3 graphics primitives
        """
        G = Graphics()
        options = self._options.copy()
        options.update(kwds)
        G._set_extra_kwds(Graphics._extract_kwds_for_show(options))

        # Check the arguments
        for o in options:
            if o not in graphplot_options and o not in G._extra_kwds:
                raise ValueError("invalid input '{}={}'".format(o, options[o]))

        for comp in self._plot_components.values():
            if not isinstance(comp, list):
                G += comp
            else:
                for item in comp:
                    G += item

        if self._options['graph_border']:
            xmin = G.xmin()
            xmax = G.xmax()
            ymin = G.ymin()
            ymax = G.ymax()
            dx = (xmax - xmin) / 10.0
            dy = (ymax - ymin) / 10.0
            border = (line([(xmin - dx, ymin - dy), (xmin - dx, ymax + dy),
                            (xmax + dx, ymax + dy), (xmax + dx, ymin - dy),
                            (xmin - dx, ymin - dy)], thickness=1.3))
            border.axes_range(xmin=(xmin - dx), xmax=(xmax + dx),
                              ymin=(ymin - dy), ymax=(ymax + dy))
            G += border
        G.set_aspect_ratio(1)
        G.axes(False)
        return G

    def layout_tree(self, root, orientation):
        """
        Compute a nice layout of a tree.

        INPUT:

        - ``root`` -- the root vertex

        - ``orientation`` -- whether to place the root at the top or at the
          bottom:

          * ``orientation="down"`` -- children are placed below their parent
          * ``orientation="top"`` -- children are placed above their parent

        EXAMPLES::

            sage: from sage.graphs.graph_plot import GraphPlot
            sage: G = graphs.HoffmanSingletonGraph()
            sage: T = Graph()
            sage: T.add_edges(G.min_spanning_tree(starting_vertex=0))
            sage: T.show(layout='tree', tree_root=0)  # indirect doctest
        """
        T = self._graph

        if not self._graph.is_tree():
            raise RuntimeError("cannot use tree layout on this graph: "
                               "self.is_tree() returns False")

        children = {root: T.neighbors(root)}

        # Always make a copy of the children because they get eaten
        stack = [list(children[root])]
        stick = [root]
        parent = {u: root for u in children[root]}
        pos = {}
        obstruction = [0.0] * T.num_verts()
        if orientation == 'down':
            o = -1
        else:
            o = 1

        def slide(v, dx):
            """
            Shift the vertex v and its descendants to the right by dx.

            Precondition: v and its descendents have already had their
            positions computed.
            """
            level = [v]
            while level:
                nextlevel = []
                for u in level:
                    x, y = pos[u]
                    x += dx
                    obstruction[y] = max(x + 1, obstruction[y])
                    pos[u] = x, y
                    nextlevel += children[u]
                level = nextlevel

        while stack:
            C = stack[-1]
            if not C:
                p = stick.pop()
                stack.pop()
                cp = children[p]
                y = o * len(stack)
                if not cp:
                    x = obstruction[y]
                    pos[p] = x, y
                else:
                    x = sum([pos[c][0] for c in cp]) / float(len(cp))
                    pos[p] = x, y
                    ox = obstruction[y]
                    if x < ox:
                        slide(p, ox - x)
                        x = ox
                obstruction[y] = x + 1
                continue

            t = C.pop()
            pt = parent[t]

            ct = [u for u in T.neighbors(t) if u != pt]
            for c in ct:
                parent[c] = t
            children[t] = ct

            stack.append(ct)
            stick.append(t)

        return pos
