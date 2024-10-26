# sage.doctest: needs sage.plot
r"""
Graph plotting in Javascript with d3.js

This module implements everything that can be used to draw graphs with `d3.js
<https://d3js.org/>`_ in Sage.

On Python's side, this is mainly done by wrapping a graph's edges and vertices
in a structure that can then be used in the javascript code. This javascript
code is then inserted into a .html file to be opened by a browser.

In the browser, the displayed page contains at the bottom right a menu
that allows to save the picture under the svg file format.

What Sage feeds javascript with is a "graph" object with the following content:

- ``vertices`` -- each vertex is a dictionary defining :

    - ``name`` -- the vertex's label

    - ``group`` -- the vertex' color (integer)

The ID of a vertex is its index in the vertex list.

- ``edges`` -- each edge is a dictionary defining :

    - ``source`` -- the ID (int) of the edge's source

    - ``target`` -- the ID (int) of the edge's destination

    - ``color`` -- the edge's color (integer)

    - ``value`` -- thickness of the edge

    - ``strength`` -- the edge's strength in the automatic layout

    - ``color`` -- color (hexadecimal code)

    - ``curve`` -- distance from the barycenter of the two endpoints and the
      center of the edge. It defines the curve of the edge, which can be useful
      for multigraphs.

- ``pos`` -- list whose `i` th element is a dictionary defining the position
  of the `i` th vertex

It also contains the definition of some numerical/boolean variables whose
definition can be found in the documentation of
:meth:`~sage.graphs.generic_graph.GenericGraph.show` : ``directed``, ``charge``,
``link_distance``, ``link_strength``, ``gravity``, ``vertex_size``,
``edge_thickness``.

.. WARNING::

    Since the d3js package is not standard yet, the javascript is fetched from
    d3js.org website by the browser. If you want to avoid that (e.g.  to protect
    your privacy or by lack of internet connection), you can install the d3js
    package for offline use by running ``sage -i d3js`` from the command line.

.. TODO::

    - Add tooltip like in `<http://bl.ocks.org/bentwonk/2514276>`_.

    - Add a zoom through scrolling (`<http://bl.ocks.org/mbostock/3681006>`_).

Authors:

- Nathann Cohen, Brice Onfroy -- July 2013 --
  Initial version of the Sage code,
  Javascript code, using examples from `d3.js <https://d3js.org/>`_.

- Thierry Monteil (June 2014): allow offline use of d3.js provided by d3js spkg.

Functions
---------
"""
from pathlib import Path
from sage.misc.temporary_file import tmp_filename
from sage.misc.lazy_import import lazy_import
lazy_import("sage.plot.colors", "rainbow")

# ****************************************************************************
#       Copyright (C) 2013 Nathann Cohen <nathann.cohen@gmail.com>
#                          Brice Onfroy  <onfroy.brice@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


def gen_html_code(G,
                  vertex_labels=True,
                  edge_labels=False,
                  vertex_partition=[],
                  vertex_colors=None,
                  edge_partition=[],
                  force_spring_layout=False,
                  charge=-120,
                  link_distance=30,
                  link_strength=2,
                  gravity=.04,
                  vertex_size=7,
                  edge_thickness=4):
    r"""
    Create a .html file showing the graph using `d3.js <https://d3js.org/>`_.

    This function returns the name of the .html file. If you want to visualize
    the actual graph use :meth:`~sage.graphs.generic_graph.GenericGraph.show`.

    INPUT:

    - ``G`` -- the graph

    - ``vertex_labels`` -- boolean (default: ``False``); whether to display
      vertex labels

    - ``edge_labels`` -- boolean (default: ``False``); whether to display edge
      labels

    - ``vertex_partition`` -- list (default: ``[]``); a list of lists
      representing a partition of the vertex set. Vertices are then colored in
      the graph according to the partition.

    - ``vertex_colors`` -- dictionary (default: ``None``); a dictionary
      representing a partition of the vertex set. Keys are colors (ignored) and
      values are lists of vertices. Vertices are then colored in the graph
      according to the partition.

    - ``edge_partition`` -- list (default: ``[]``); same as
      ``vertex_partition``, with edges instead

    - ``force_spring_layout`` -- boolean (default: ``False``); whether to take
      previously computed position of nodes into account if there is one, or to
      compute a spring layout

    - ``vertex_size`` -- integer (default: 7); the size of a vertex' circle

    - ``edge_thickness`` -- integer (default: 4); thickness of an edge

    - ``charge`` -- integer (default: -120); the vertices' charge. Defines
      how they repulse each other. See
      `<https://github.com/mbostock/d3/wiki/Force-Layout>`_ for more
      information

    - ``link_distance`` -- integer (default: 30); see
      `<https://github.com/mbostock/d3/wiki/Force-Layout>`_ for more
      information

    - ``link_strength`` -- integer (default: 2); see
      `<https://github.com/mbostock/d3/wiki/Force-Layout>`_ for more
      information

    - ``gravity`` -- float (default: 0.04); see
      `<https://github.com/mbostock/d3/wiki/Force-Layout>`_ for more
      information

    .. WARNING::

        Since the d3js package is not standard yet, the javascript is fetched
        from d3js.org website by the browser. If you want to avoid that (e.g.
        to protect your privacy or by lack of internet connection), you can
        install the d3js package for offline use by running ``sage -i d3js``
        from the command line.

    EXAMPLES::

        sage: graphs.RandomTree(50).show(method='js')                           # optional - internet

        sage: g = graphs.PetersenGraph()
        sage: g.show(method='js', vertex_partition=g.coloring())                # optional - internet

        sage: graphs.DodecahedralGraph().show(method='js',                      # optional - internet
        ....:                                 force_spring_layout=True)

        sage: graphs.DodecahedralGraph().show(method='js')                      # optional - internet

        sage: # needs sage.combinat
        sage: g = digraphs.DeBruijn(2, 2)
        sage: g.allow_multiple_edges(True)
        sage: g.add_edge("10", "10", "a")
        sage: g.add_edge("10", "10", "b")
        sage: g.add_edge("10", "10", "c")
        sage: g.add_edge("10", "10", "d")
        sage: g.add_edge("01", "11", "1")
        sage: g.show(method='js', vertex_labels=True, edge_labels=True,         # optional - internet
        ....:        link_distance=200, gravity=.05, charge=-500,
        ....:        edge_partition=[[("11", "12", "2"), ("21", "21", "a")]],
        ....:        edge_thickness=4)

    TESTS::

        sage: from sage.graphs.graph_plot_js import gen_html_code
        sage: filename = gen_html_code(graphs.PetersenGraph())

    :issue:`17370`::

        sage: filename = gen_html_code(graphs.CompleteBipartiteGraph(4, 5))

    In the generated html code, the source (resp. target) of a link is the index
    of the node in the list defining the names of the nodes. We check that the
    order is correct (:issue:`27460`)::

        sage: filename = gen_html_code(DiGraph({1: [10]}))
        sage: with open(filename, 'r') as f:
        ....:     data = f.read()
        sage: nodes = data.partition('"nodes":')[2]; nodes
        ...[{..."name": "10"...}, {..."name": "1"...}]...
        sage: links = data.partition('"links":')[2]
        sage: '"source": 1' in links and '"target": 0' in links
        True
    """
    directed = G.is_directed()
    multiple_edges = G.has_multiple_edges()

    # Associated an integer to each vertex
    v_to_id = {v: i for i, v in enumerate(G)}

    # Vertex colors
    if vertex_colors is not None:
        vertex_partition = list(vertex_colors.values())
    len_vertex_partition = len(vertex_partition)
    color = {i: len_vertex_partition for i in range(G.order())}
    for i, l in enumerate(vertex_partition):
        for v in l:
            color[v_to_id[v]] = i

    # Vertex list
    # Data for vertex v must be at position v_to_id[v] in list nodes
    nodes = [{"name": str(v), "group": str(color[v_to_id[v]])} for v in G]

    # Edge colors.
    edge_color_default = "#aaa"
    color_list = rainbow(len(edge_partition))
    edge_color = {}
    for i, l in enumerate(edge_partition):
        for e in l:
            u, v, label = e if len(e) == 3 else e+(None,)
            edge_color[u, v, label] = color_list[i]
            if not directed:
                edge_color[v, u, label] = color_list[i]

    # Edge list
    edges = []
    seen = {}  # How many times has this edge been seen ?

    for u, v, l in G.edge_iterator():

        # Edge color
        color = edge_color.get((u, v, l), edge_color_default)

        # Computes the curve of the edge
        curve = 0

        # Loop ?
        if u == v:
            seen[u, v] = seen.get((u, v), 0) + 1
            curve = seen[u, v] * 10 + 10

        # For directed graphs, one also has to take into accounts
        # edges in the opposite direction
        elif directed:
            if G.has_edge(v, u):
                seen[u, v] = seen.get((u, v), 0) + 1
                curve = seen[u, v] * 15
            else:
                if multiple_edges and len(G.edge_label(u, v)) != 1:
                    # Multiple edges. The first one has curve 15, then
                    # -15, then 30, then -30, ...
                    seen[u, v] = seen.get((u, v), 0) + 1
                    curve = (1 if seen[u, v] % 2 else -1) * (seen[u, v] // 2) * 15

        elif not directed and multiple_edges:
            # Same formula as above for multiple edges
            if len(G.edge_label(u, v)) != 1:
                seen[u, v] = seen.get((u, v), 0) + 1
                curve = (1 if seen[u, v] % 2 else -1) * (seen[u, v] // 2) * 15

        # Adding the edge to the list
        # The source (resp. target) is the index of u (resp. v) in list nodes
        edges.append({"source": v_to_id[u],
                      "target": v_to_id[v],
                      "strength": 0,
                      "color": color,
                      "curve": curve,
                      "name": str(l) if edge_labels else ""})

    loops = [e for e in edges if e["source"] == e["target"]]
    edges = [e for e in edges if e["source"] != e["target"]]

    # Defines the vertices' layout if possible
    Gpos = G.get_pos()
    pos = []
    if Gpos is not None and force_spring_layout is False:
        charge = 0
        link_strength = 0
        gravity = 0

        for v in G:
            x, y = Gpos[v]
            pos.append([float(x), float(-y)])

    # Encodes the data as a JSON string
    from json import JSONEncoder
    string = JSONEncoder().encode({"nodes": nodes,
                                   "links": edges,
                                   "loops": loops,
                                   "pos": pos,
                                   "directed": G.is_directed(),
                                   "charge": int(charge),
                                   "link_distance": int(link_distance),
                                   "link_strength": int(link_strength),
                                   "gravity": float(gravity),
                                   "vertex_labels": bool(vertex_labels),
                                   "edge_labels": bool(edge_labels),
                                   "vertex_size": int(vertex_size),
                                   "edge_thickness": int(edge_thickness)})

    from sage.env import SAGE_EXTCODE, SAGE_SHARE
    with open(Path(SAGE_EXTCODE) / "graphs" / "graph_plot_js.html") as f:
        js_code = f.read().replace("// GRAPH_DATA_HEREEEEEEEEEEE", string)

    # Add d3.js script depending on whether d3js package is installed.
    d3js_filepath = Path(SAGE_SHARE) / 'd3js' / 'd3.min.js'
    if d3js_filepath.exists():
        with open(d3js_filepath) as d3js_code_file:
            d3js_script = '<script>' + d3js_code_file.read() + '</script>'
    else:
        d3js_script = '<script src="http://d3js.org/d3.v3.min.js"></script>'
    js_code = js_code.replace('// D3JS_SCRIPT_HEREEEEEEEEEEE', d3js_script)

    # Writes the temporary .html file
    filename = tmp_filename(ext='.html')
    with open(filename, 'w') as f:
        f.write(js_code)

    return filename
