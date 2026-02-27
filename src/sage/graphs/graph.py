r"""
Undirected graphs

This module implements functions and operations involving undirected graphs.

{INDEX_OF_METHODS}

AUTHORS:

- Robert L. Miller (2006-10-22): initial version

- William Stein (2006-12-05): Editing

- Robert L. Miller (2007-01-13): refactoring, adjusting for NetworkX-0.33, fixed
   plotting bugs (2007-01-23): basic tutorial, edge labels, loops, multiple
   edges and arcs (2007-02-07): graph6 and sparse6 formats, matrix input

- Emily Kirkmann (2007-02-11): added graph_border option to plot and show

- Robert L. Miller (2007-02-12): vertex color-maps, graph boundaries, graph6
   helper functions in Cython

- Robert L. Miller Sage Days 3 (2007-02-17-21): 3d plotting in Tachyon

- Robert L. Miller (2007-02-25): display a partition

- Robert L. Miller (2007-02-28): associate arbitrary objects to vertices, edge
   and arc label display (in 2d), edge coloring

- Robert L. Miller (2007-03-21): Automorphism group, isomorphism check,
   canonical label

- Robert L. Miller (2007-06-07-09): NetworkX function wrapping

- Michael W. Hansen (2007-06-09): Topological sort generation

- Emily Kirkman, Robert L. Miller Sage Days 4: Finished wrapping NetworkX

- Emily Kirkman (2007-07-21): Genus (including circular planar, all embeddings
   and all planar embeddings), all paths, interior paths

- Bobby Moretti (2007-08-12): fixed up plotting of graphs with edge colors
   differentiated by label

- Jason Grout (2007-09-25): Added functions, bug fixes, and general enhancements

- Robert L. Miller (Sage Days 7): Edge labeled graph isomorphism

- Tom Boothby (Sage Days 7): Miscellaneous awesomeness

- Tom Boothby (2008-01-09): Added graphviz output

- David Joyner (2009-2): Fixed docstring bug related to GAP.

- Stephen Hartke (2009-07-26): Fixed bug in blocks_and_cut_vertices() that
   caused an incorrect result when the vertex 0 was a cut vertex.

- Stephen Hartke (2009-08-22): Fixed bug in blocks_and_cut_vertices() where the
   list of cut_vertices is not treated as a set.

- Anders Jonsson (2009-10-10): Counting of spanning trees and out-trees added.

- Nathann Cohen (2009-09) : Cliquer, Connectivity, Flows and everything that
                             uses Linear Programming and class numerical.MIP

- Nicolas M. Thiery (2010-02): graph layout code refactoring, dot2tex/graphviz
  interface

- David Coudert (2012-04) : Reduction rules in vertex_cover.

- Birk Eisermann (2012-06): added recognition of weakly chordal graphs and
                            long-hole-free / long-antihole-free graphs

- Alexandre P. Zuge (2013-07): added join operation.

- Amritanshu Prasad (2014-08): added clique polynomial

- Julian Rüth (2018-06-21): upgrade to NetworkX 2

- David Coudert (2018-10-07): cleaning

- Amanda Francis, Caitlin Lienkaemper, Kate Collins, Rajat Mittal (2019-03-10):
  methods for computing effective resistance

- Amanda Francis, Caitlin Lienkaemper, Kate Collins, Rajat Mittal (2019-03-19):
  most_common_neighbors and common_neighbors_matrix added.

- Jean-Florent Raymond (2019-04): is_redundant, is_dominating,
   private_neighbors

- Cyril Bouvier (2024-11): is_module

Graph Format
------------

Supported formats
~~~~~~~~~~~~~~~~~

Sage Graphs can be created from a wide range of inputs. A few examples are
covered here.

- NetworkX dictionary format:

   ::

       sage: d = {0: [1,4,5], 1: [2,6], 2: [3,7], 3: [4,8], 4: [9], \
       ....: 5: [7, 8], 6: [8,9], 7: [9]}
       sage: G = Graph(d); G
       Graph on 10 vertices
       sage: G.plot().show()    # or G.show()                                           # needs sage.plot

- A NetworkX graph:

   ::

       sage: # needs networkx
       sage: import networkx
       sage: K = networkx.complete_bipartite_graph(12,7)
       sage: G = Graph(K)
       sage: G.degree()
       [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 12, 12, 12, 12, 12, 12, 12]

- graph6 or sparse6 format:

   ::

       sage: s = ':I`AKGsaOs`cI]Gb~'
       sage: G = Graph(s, sparse=True); G
       Looped multi-graph on 10 vertices
       sage: G.plot().show()    # or G.show()                                           # needs sage.plot

   Note that the ``\`` character is an escape character in Python, and also a
   character used by graph6 strings:

   ::

       sage: G = Graph('Ihe\n@GUA')
       Traceback (most recent call last):
       ...
       RuntimeError: the string (Ihe) seems corrupt: for n = 10, the string is too short

   In Python, the escaped character ``\`` is represented by ``\\``:

   ::

       sage: G = Graph('Ihe\\n@GUA')
       sage: G.plot().show()    # or G.show()                                           # needs sage.plot

- adjacency matrix: In an adjacency matrix, each column and each row represent a
   vertex. If a 1 shows up in row `i`, column `j`, there is an edge `(i,j)`.

   ::

       sage: # needs sage.modules
       sage: M = Matrix([(0,1,0,0,1,1,0,0,0,0), (1,0,1,0,0,0,1,0,0,0),
       ....:             (0,1,0,1,0,0,0,1,0,0), (0,0,1,0,1,0,0,0,1,0),
       ....:             (1,0,0,1,0,0,0,0,0,1), (1,0,0,0,0,0,0,1,1,0), (0,1,0,0,0,0,0,0,1,1),
       ....:             (0,0,1,0,0,1,0,0,0,1), (0,0,0,1,0,1,1,0,0,0), (0,0,0,0,1,0,1,1,0,0)])
       sage: M
       [0 1 0 0 1 1 0 0 0 0]
       [1 0 1 0 0 0 1 0 0 0]
       [0 1 0 1 0 0 0 1 0 0]
       [0 0 1 0 1 0 0 0 1 0]
       [1 0 0 1 0 0 0 0 0 1]
       [1 0 0 0 0 0 0 1 1 0]
       [0 1 0 0 0 0 0 0 1 1]
       [0 0 1 0 0 1 0 0 0 1]
       [0 0 0 1 0 1 1 0 0 0]
       [0 0 0 0 1 0 1 1 0 0]
       sage: G = Graph(M); G
       Graph on 10 vertices
       sage: G.plot().show()    # or G.show()                                           # needs sage.plot

- incidence matrix: In an incidence matrix, each row represents a vertex and
   each column represents an edge.

   ::

       sage: # needs sage.modules
       sage: M = Matrix([(-1, 0, 0, 0, 1, 0, 0, 0, 0, 0,-1, 0, 0, 0, 0),
       ....:             ( 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 0),
       ....:             ( 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0),
       ....:             ( 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0),
       ....:             ( 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1),
       ....:             ( 0, 0, 0, 0, 0,-1, 0, 0, 0, 1, 1, 0, 0, 0, 0),
       ....:             ( 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 1, 0, 0, 0),
       ....:             ( 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 1, 0, 0),
       ....:             ( 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 1, 0),
       ....:             ( 0, 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0, 1)])
       sage: M
       [-1  0  0  0  1  0  0  0  0  0 -1  0  0  0  0]
       [ 1 -1  0  0  0  0  0  0  0  0  0 -1  0  0  0]
       [ 0  1 -1  0  0  0  0  0  0  0  0  0 -1  0  0]
       [ 0  0  1 -1  0  0  0  0  0  0  0  0  0 -1  0]
       [ 0  0  0  1 -1  0  0  0  0  0  0  0  0  0 -1]
       [ 0  0  0  0  0 -1  0  0  0  1  1  0  0  0  0]
       [ 0  0  0  0  0  0  0  1 -1  0  0  1  0  0  0]
       [ 0  0  0  0  0  1 -1  0  0  0  0  0  1  0  0]
       [ 0  0  0  0  0  0  0  0  1 -1  0  0  0  1  0]
       [ 0  0  0  0  0  0  1 -1  0  0  0  0  0  0  1]
       sage: G = Graph(M); G
       Graph on 10 vertices
       sage: G.plot().show()    # or G.show()                                           # needs sage.plot
       sage: DiGraph(matrix(2, [0,0,-1,1]), format='incidence_matrix')
       Traceback (most recent call last):
       ...
       ValueError: there must be two nonzero entries (-1 & 1) per column

- a list of edges::

       sage: g = Graph([(1, 3), (3, 8), (5, 2)]); g
       Graph on 5 vertices

- an igraph Graph::

       sage: import igraph                                 # optional - python_igraph
       sage: g = Graph(igraph.Graph([(1,3),(3,2),(0,2)]))  # optional - python_igraph
       sage: g                                             # optional - python_igraph
       Graph on 4 vertices

Generators
----------

Use ``graphs(n)`` to iterate through all non-isomorphic graphs of given size::

    sage: for g in graphs(4):
    ....:     print(g.degree_sequence())
    [0, 0, 0, 0]
    [1, 1, 0, 0]
    [2, 1, 1, 0]
    [3, 1, 1, 1]
    [1, 1, 1, 1]
    [2, 2, 1, 1]
    [2, 2, 2, 0]
    [3, 2, 2, 1]
    [2, 2, 2, 2]
    [3, 3, 2, 2]
    [3, 3, 3, 3]

Similarly ``graphs()`` will iterate through all graphs. The complete graph of 4
vertices is of course the smallest graph with chromatic number bigger than
three::

    sage: for g in graphs():
    ....:     if g.chromatic_number() > 3:
    ....:         break
    sage: g.is_isomorphic(graphs.CompleteGraph(4))
    True

For some commonly used graphs to play with, type::

    sage: graphs.[tab]          # not tested

and hit {tab}. Most of these graphs come with their own custom plot, so you can
see how people usually visualize these graphs.

::

    sage: G = graphs.PetersenGraph()
    sage: G.plot().show()    # or G.show()                                              # needs sage.plot
    sage: G.degree_histogram()
    [0, 0, 0, 10]
    sage: G.adjacency_matrix()                                                          # needs sage.modules
    [0 1 0 0 1 1 0 0 0 0]
    [1 0 1 0 0 0 1 0 0 0]
    [0 1 0 1 0 0 0 1 0 0]
    [0 0 1 0 1 0 0 0 1 0]
    [1 0 0 1 0 0 0 0 0 1]
    [1 0 0 0 0 0 0 1 1 0]
    [0 1 0 0 0 0 0 0 1 1]
    [0 0 1 0 0 1 0 0 0 1]
    [0 0 0 1 0 1 1 0 0 0]
    [0 0 0 0 1 0 1 1 0 0]

::

    sage: S = G.subgraph([0,1,2,3])
    sage: S.plot().show()    # or S.show()                                              # needs sage.plot
    sage: S.density()
    1/2

::

    sage: G = GraphQuery(display_cols=['graph6'], num_vertices=7, diameter=5)
    sage: L = G.get_graphs_list()
    sage: graphs_list.show_graphs(L)                                                    # needs sage.plot

.. _Graph:labels:

Labels
------

Each vertex can have any hashable object as a label. These are things like
strings, numbers, and tuples. Each edge is given a default label of ``None``,
but if specified, edges can have any label at all. Edges between vertices `u`
and `v` are represented typically as ``(u, v, l)``, where ``l`` is the label for
the edge.

Note that vertex labels themselves cannot be mutable items::

    sage: M = Matrix([[0,0], [0,0]])                                                    # needs sage.modules
    sage: G = Graph({ 0 : { M : None } })                                               # needs sage.modules
    Traceback (most recent call last):
    ...
    TypeError: ...mutable matrices are unhashable...

However, if one wants to define a dictionary, with the same keys and arbitrary
objects for entries, one can make that association::

    sage: d = {0 : graphs.DodecahedralGraph(), 1 : graphs.FlowerSnark(), \
    ....: 2 : graphs.MoebiusKantorGraph(), 3 : graphs.PetersenGraph() }
    sage: d[2]
    Moebius-Kantor Graph: Graph on 16 vertices
    sage: T = graphs.TetrahedralGraph()
    sage: T.vertices(sort=True)
    [0, 1, 2, 3]
    sage: T.set_vertices(d)
    sage: T.get_vertex(1)
    Flower Snark: Graph on 20 vertices

Database
--------

There is a database available for searching for graphs that satisfy a certain
set of parameters, including number of vertices and edges, density, maximum and
minimum degree, diameter, radius, and connectivity. To see a list of all search
parameter keywords broken down by their designated table names, type ::

    sage: graph_db_info()
    {...}

For more details on data types or keyword input, enter ::

    sage: GraphQuery?    # not tested

The results of a query can be viewed with the show method, or can be viewed
individually by iterating through the results ::

    sage: Q = GraphQuery(display_cols=['graph6'],num_vertices=7, diameter=5)
    sage: Q.show()
    Graph6
    --------------------
    F?`po
    F?gqg
    F@?]O
    F@OKg
    F@R@o
    FA_pW
    FEOhW
    FGC{o
    FIAHo

Show each graph as you iterate through the results::

    sage: for g in Q:                                                                   # needs sage.plot
    ....:     show(g)

Visualization
-------------

To see a graph `G` you are working with, there are three main options. You can
view the graph in two dimensions via matplotlib with ``show()``. ::

    sage: G = graphs.RandomGNP(15,.3)
    sage: G.show()                                                                      # needs sage.plot

And you can view it in three dimensions with ``show3d()``. ::

    sage: G.show3d()                                                                    # needs sage.plot

Or it can be rendered with `\LaTeX`.  This requires the right additions to a
standard `\mbox{\rm\TeX}` installation.  Then standard Sage commands, such as
``view(G)`` will display the graph, or ``latex(G)`` will produce a string
suitable for inclusion in a `\LaTeX` document.  More details on this are at the
:mod:`sage.graphs.graph_latex` module. ::

    sage: from sage.graphs.graph_latex import check_tkz_graph
    sage: check_tkz_graph()  # random - depends on TeX installation
    sage: latex(G)
    \begin{tikzpicture}
    ...
    \end{tikzpicture}

Mutability
----------

Graphs are mutable, and thus unusable as dictionary keys, unless
``data_structure="static_sparse"`` is used::

    sage: G = graphs.PetersenGraph()
    sage: {G:1}[G]
    Traceback (most recent call last):
    ...
    TypeError: ...This graph is mutable, and thus not hashable...
    sage: G_immutable = Graph(G, immutable=True)
    sage: G_immutable == G
    True
    sage: {G_immutable:1}[G_immutable]
    1

Methods
-------
"""


# ****************************************************************************
#       Copyright (C) 2006-2007 Robert L. Miller <rlmillster@gmail.com>
#                          2018 Julian Rüth <julian.rueth@fsfe.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import itertools

from copy import copy
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
import sage.graphs.generic_graph_pyx as generic_graph_pyx
from sage.graphs.generic_graph import GenericGraph
from sage.graphs.independent_sets import IndependentSets
from sage.misc.rest_index_of_methods import doc_index, gen_thematic_rest_table_index
from sage.graphs.views import EdgesView
from sage.parallel.decorate import parallel
from sage.misc.lazy_import import lazy_import, LazyImport
from sage.features.mcqd import Mcqd
from sage.misc.cachefunc import cached_method

lazy_import('sage.graphs.mcqd', ['mcqd'],
            feature=Mcqd())


class Graph(GenericGraph):
    r"""
    Undirected graph.

    A graph is a set of vertices connected by edges. See the
    :wikipedia:`Graph_(mathematics)` for more information. For a collection of
    pre-defined graphs, see the :mod:`~sage.graphs.graph_generators` module.

    A :class:`Graph` object has many methods whose list can be obtained by
    typing ``g.<tab>`` (i.e. hit the :kbd:`Tab` key) or by reading the documentation
    of :mod:`~sage.graphs.graph`, :mod:`~sage.graphs.generic_graph`, and
    :mod:`~sage.graphs.digraph`.

    INPUT:

    By default, a :class:`Graph` object is simple (i.e. no *loops* nor *multiple
    edges*) and unweighted. This can be easily tuned with the appropriate flags
    (see below).

    - ``data`` -- can be any of the following (see the ``format`` argument):

      #. ``Graph()`` -- build a graph on 0 vertices.

      #. ``Graph(5)`` -- return an edgeless graph on the 5 vertices 0,...,4.

      #. ``Graph([list_of_vertices, list_of_edges])`` -- returns a graph with
         given vertices/edges.

         To bypass auto-detection, prefer the more explicit
         ``Graph([V, E], format='vertices_and_edges')``.

      #. ``Graph(list_of_edges)`` -- return a graph with a given list of edges
         (see documentation of
         :meth:`~sage.graphs.generic_graph.GenericGraph.add_edges`).

         To bypass auto-detection, prefer the more explicit
         ``Graph(L, format='list_of_edges')``.

      #. ``Graph({1: [2, 3, 4], 3: [4]})`` -- return a graph by associating to
         each vertex the list of its neighbors.

         To bypass auto-detection, prefer the more explicit
         ``Graph(D, format='dict_of_lists')``.

      #. ``Graph({1: {2: 'a', 3:'b'} ,3:{2:'c'}})`` -- return a graph by
         associating a list of neighbors to each vertex and providing its edge
         label.

         To bypass auto-detection, prefer the more explicit
         ``Graph(D, format='dict_of_dicts')``.

         For graphs with multiple edges, you can provide a list of labels
         instead, e.g.: ``Graph({1: {2: ['a1', 'a2'], 3:['b']} ,3:{2:['c']}})``.

      #. ``Graph(a_symmetric_matrix)`` -- return a graph with given (weighted)
         adjacency matrix (see documentation of
         :meth:`~sage.graphs.generic_graph.GenericGraph.adjacency_matrix`).

         To bypass auto-detection, prefer the more explicit ``Graph(M,
         format='adjacency_matrix')``. To take weights into account, use
         ``format='weighted_adjacency_matrix'`` instead.

      #. ``Graph(a_nonsymmetric_matrix)`` -- return a graph with given incidence
         matrix (see documentation of
         :meth:`~sage.graphs.generic_graph.GenericGraph.incidence_matrix`).

         To bypass auto-detection, prefer the more explicit
         ``Graph(M, format='incidence_matrix')``.

      #. ``Graph([V, f])`` -- return a graph from a vertex set ``V`` and a
         *symmetric* function ``f``. The graph contains an edge `u,v` whenever
         ``f(u,v)`` is ``True``.. Example: ``Graph([ [1..10], lambda x,y:
         abs(x-y).is_square()])``

      #. ``Graph(':I`ES@obGkqegW~')`` -- return a graph from a graph6 or sparse6
         string (see documentation of :meth:`graph6_string` or
         :meth:`sparse6_string`).

      #. ``Graph(a_seidel_matrix, format='seidel_adjacency_matrix')`` -- return
         a graph with a given Seidel adjacency matrix (see documentation of
         :meth:`seidel_adjacency_matrix`).

      #. ``Graph(another_graph)`` -- return a graph from a Sage (di)graph,
         `pygraphviz <https://pygraphviz.github.io/>`__ graph, `NetworkX
         <https://networkx.github.io/>`__ graph, or `igraph
         <http://igraph.org/python/>`__ graph.

    - ``pos`` -- a positioning dictionary (cf. documentation of
      :meth:`~sage.graphs.generic_graph.GenericGraph.layout`). For example, to
      draw 4 vertices on a square::

         {0: [-1,-1],
          1: [ 1,-1],
          2: [ 1, 1],
          3: [-1, 1]}

    - ``name`` -- (must be an explicitly named parameter, i.e.,
       ``name='complete')`` gives the graph a name

    - ``loops`` -- boolean (default: ``None``); whether to allow loops (ignored
      if data is an instance of the ``Graph`` class)

    - ``multiedges`` -- boolean (default: ``None``); whether to allow multiple
      edges (ignored if data is an instance of the ``Graph`` class)

    - ``weighted`` -- boolean (default: ``None``); whether graph thinks of
      itself as weighted or not. See
      :meth:`~sage.graphs.generic_graph.GenericGraph.weighted`.

    - ``format`` -- if set to ``None`` (default), :class:`Graph` tries to guess
      input's format. To avoid this possibly time-consuming step, one of the
      following values can be specified (see description above): ``'int'``,
      ``'graph6'``, ``'sparse6'``, ``'rule'``, ``'list_of_edges'``,
      ``'dict_of_lists'``, ``'dict_of_dicts'``, ``'adjacency_matrix'``,
      ``'weighted_adjacency_matrix'``, ``'seidel_adjacency_matrix'``,
      ``'incidence_matrix'``, ``"NX"``, ``'igraph'``.

    - ``sparse`` -- boolean (default: ``True``); ``sparse=True`` is an alias for
      ``data_structure="sparse"``, and ``sparse=False`` is an alias for
      ``data_structure="dense"``.

    - ``data_structure`` -- one of the following (for more information, see
      :mod:`~sage.graphs.base.overview`)

       * ``'dense'`` -- selects the :mod:`~sage.graphs.base.dense_graph`
         backend.

       * ``'sparse'`` -- selects the :mod:`~sage.graphs.base.sparse_graph`
         backend.

       * ``'static_sparse'`` -- selects the
         :mod:`~sage.graphs.base.static_sparse_backend` (this backend is faster
         than the sparse backend and smaller in memory, and it is immutable, so
         that the resulting graphs can be used as dictionary keys).

    - ``immutable`` -- boolean (default: ``False``); whether to create a
      immutable graph. Note that ``immutable=True`` is actually a shortcut for
      ``data_structure='static_sparse'``. Set to ``False`` by default.

    - ``hash_labels`` -- boolean (default: ``None``); whether to include edge
      labels during hashing. This parameter defaults to ``True`` if the graph is
      weighted. This parameter is ignored if the graph is mutable.
      Beware that trying to hash unhashable labels will raise an error.

    - ``vertex_labels`` -- boolean (default: ``True``); whether to allow any
      object as a vertex (slower), or only the integers `0,...,n-1`, where `n`
      is the number of vertices.

    - ``convert_empty_dict_labels_to_None`` -- this arguments sets the default
      edge labels used by NetworkX (empty dictionaries) to be replaced by
      ``None``, the default Sage edge label. It is set to ``True`` iff a
      NetworkX graph is on the input.

    EXAMPLES:

    We illustrate the first seven input formats (the other two involve packages
    that are currently not standard in Sage):

    #. An integer giving the number of vertices::

        sage: g = Graph(5); g
        Graph on 5 vertices
        sage: g.vertices(sort=True)
        [0, 1, 2, 3, 4]
        sage: g.edges(sort=False)
        []

    #. A dictionary of dictionaries::

        sage: g = Graph({0:{1:'x',2:'z',3:'a'}, 2:{5:'out'}}); g
        Graph on 5 vertices

       The labels ('x', 'z', 'a', 'out') are labels for edges. For example,
       'out' is the label for the edge on 2 and 5. Labels can be used as
       weights, if all the labels share some common parent.::

        sage: a, b, c, d, e, f = sorted(SymmetricGroup(3))                              # needs sage.groups
        sage: Graph({b: {d: 'c', e: 'p'}, c: {d: 'p', e: 'c'}})                         # needs sage.groups
        Graph on 4 vertices

    #. A dictionary of lists::

        sage: g = Graph({0:[1,2,3], 2:[4]}); g
        Graph on 5 vertices

    #. A list of vertices and a function describing adjacencies. Note that the
       list of vertices and the function must be enclosed in a list (i.e., [list
       of vertices, function]).

       Construct the Paley graph over GF(13).::

          sage: g = Graph([GF(13), lambda i,j: i!=j and (i-j).is_square()])             # needs sage.rings.finite_rings
          sage: g.vertices(sort=True)                                                   # needs sage.rings.finite_rings
          [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
          sage: g.adjacency_matrix()                                                    # needs sage.modules sage.rings.finite_rings
          [0 1 0 1 1 0 0 0 0 1 1 0 1]
          [1 0 1 0 1 1 0 0 0 0 1 1 0]
          [0 1 0 1 0 1 1 0 0 0 0 1 1]
          [1 0 1 0 1 0 1 1 0 0 0 0 1]
          [1 1 0 1 0 1 0 1 1 0 0 0 0]
          [0 1 1 0 1 0 1 0 1 1 0 0 0]
          [0 0 1 1 0 1 0 1 0 1 1 0 0]
          [0 0 0 1 1 0 1 0 1 0 1 1 0]
          [0 0 0 0 1 1 0 1 0 1 0 1 1]
          [1 0 0 0 0 1 1 0 1 0 1 0 1]
          [1 1 0 0 0 0 1 1 0 1 0 1 0]
          [0 1 1 0 0 0 0 1 1 0 1 0 1]
          [1 0 1 1 0 0 0 0 1 1 0 1 0]

       Construct the line graph of a complete graph.::

          sage: g = graphs.CompleteGraph(4)
          sage: line_graph = Graph([g.edges(sort=True, labels=false),
          ....:                     lambda i,j: len(set(i).intersection(set(j)))>0],
          ....:                    loops=False)
          sage: line_graph.vertices(sort=True)
          [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
          sage: line_graph.adjacency_matrix()                                           # needs sage.modules
          [0 1 1 1 1 0]
          [1 0 1 1 0 1]
          [1 1 0 0 1 1]
          [1 1 0 0 1 1]
          [1 0 1 1 0 1]
          [0 1 1 1 1 0]

    #. A graph6 or sparse6 string: Sage automatically recognizes whether a
       string is in graph6 or sparse6 format::

           sage: s = ':I`AKGsaOs`cI]Gb~'
           sage: Graph(s, sparse=True)
           Looped multi-graph on 10 vertices

       ::

           sage: G = Graph('G?????')
           sage: G = Graph("G'?G?C")
           Traceback (most recent call last):
           ...
           RuntimeError: the string seems corrupt: valid characters are
           ?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
           sage: G = Graph('G??????')
           Traceback (most recent call last):
           ...
           RuntimeError: the string (G??????) seems corrupt: for n = 8, the string is too long

       ::

          sage: G = Graph(":I'AKGsaOs`cI]Gb~")
          Traceback (most recent call last):
          ...
          RuntimeError: the string seems corrupt: valid characters are
          ?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~

       There are also list functions to take care of lists of graphs::

           sage: s = ':IgMoqoCUOqeb\n:I`AKGsaOs`cI]Gb~\n:I`EDOAEQ?PccSsge\\N\n'
           sage: graphs_list.from_sparse6(s)
           [Looped multi-graph on 10 vertices,
            Looped multi-graph on 10 vertices,
            Looped multi-graph on 10 vertices]

    #. A Sage matrix:
       Note: If format is not specified, then Sage assumes a symmetric square
       matrix is an adjacency matrix, otherwise an incidence matrix.

       - an adjacency matrix::

            sage: M = graphs.PetersenGraph().am(); M                                    # needs sage.modules
            [0 1 0 0 1 1 0 0 0 0]
            [1 0 1 0 0 0 1 0 0 0]
            [0 1 0 1 0 0 0 1 0 0]
            [0 0 1 0 1 0 0 0 1 0]
            [1 0 0 1 0 0 0 0 0 1]
            [1 0 0 0 0 0 0 1 1 0]
            [0 1 0 0 0 0 0 0 1 1]
            [0 0 1 0 0 1 0 0 0 1]
            [0 0 0 1 0 1 1 0 0 0]
            [0 0 0 0 1 0 1 1 0 0]
            sage: Graph(M)                                                              # needs sage.modules
            Graph on 10 vertices

         ::

            sage: Graph(matrix([[1,2], [2,4]]), loops=True, sparse=True)                # needs sage.modules
            Looped multi-graph on 2 vertices

            sage: M = Matrix([[0,1,-1], [1,0,-1/2], [-1,-1/2,0]]); M                    # needs sage.modules
            [   0    1   -1]
            [   1    0 -1/2]
            [  -1 -1/2    0]
            sage: G = Graph(M, sparse=True); G                                          # needs sage.modules
            Graph on 3 vertices
            sage: G.weighted()                                                          # needs sage.modules
            True

       - an incidence matrix::

            sage: M = Matrix(6, [-1,0,0,0,1, 1,-1,0,0,0, 0,1,-1,0,0,                    # needs sage.modules
            ....:                0,0,1,-1,0, 0,0,0,1,-1, 0,0,0,0,0]); M
            [-1  0  0  0  1]
            [ 1 -1  0  0  0]
            [ 0  1 -1  0  0]
            [ 0  0  1 -1  0]
            [ 0  0  0  1 -1]
            [ 0  0  0  0  0]
            sage: Graph(M)                                                              # needs sage.modules
            Graph on 6 vertices

            sage: Graph(Matrix([[1],[1],[1]]))                                          # needs sage.modules
            Traceback (most recent call last):
            ...
            ValueError: there must be one or two nonzero entries per column
            in an incidence matrix, got entries [1, 1, 1] in column 0
            sage: Graph(Matrix([[1],[1],[0]]))                                          # needs sage.modules
            Graph on 3 vertices

            sage: M = Matrix([[0,1,-1], [1,0,-1], [-1,-1,0]]); M                        # needs sage.modules
            [ 0  1 -1]
            [ 1  0 -1]
            [-1 -1  0]
            sage: Graph(M, sparse=True)                                                 # needs sage.modules
            Graph on 3 vertices

            sage: M = Matrix([[0,1,1], [1,0,1], [-1,-1,0]]); M                          # needs sage.modules
            [ 0  1  1]
            [ 1  0  1]
            [-1 -1  0]
            sage: Graph(M)                                                              # needs sage.modules
            Traceback (most recent call last):
            ...
            ValueError: there must be one or two nonzero entries per column
            in an incidence matrix, got entries [1, 1] in column 2

        Check that :issue:`9714` is fixed::

            sage: # needs sage.modules
            sage: MA = Matrix([[1,2,0], [0,2,0], [0,0,1]])
            sage: GA = Graph(MA, format='adjacency_matrix')
            sage: MI = GA.incidence_matrix(oriented=False); MI
            [2 1 1 0 0 0]
            [0 1 1 2 2 0]
            [0 0 0 0 0 2]
            sage: Graph(MI).edges(sort=True, labels=None)
            [(0, 0), (0, 1), (0, 1), (1, 1), (1, 1), (2, 2)]

            sage: M = Matrix([[1], [-1]]); M                                            # needs sage.modules
            [ 1]
            [-1]
            sage: Graph(M).edges(sort=True)                                             # needs sage.modules
            [(0, 1, None)]

    #. A Seidel adjacency matrix::

          sage: from sage.combinat.matrices.hadamard_matrix import (                    # needs sage.combinat sage.modules
          ....:  regular_symmetric_hadamard_matrix_with_constant_diagonal as rshcd)
          sage: m = rshcd(16,1) - matrix.identity(16)                                   # needs sage.combinat sage.modules
          sage: Graph(m,                                                                # needs sage.combinat sage.modules
          ....:       format='seidel_adjacency_matrix').is_strongly_regular(parameters=True)
          (16, 6, 2, 2)

    #. List of edges, or labelled edges::

          sage: g = Graph([(1, 3), (3, 8), (5, 2)]); g
          Graph on 5 vertices

          sage: g = Graph([(1, 2, "Peace"), (7, -9, "and"), (77, 2, "Love")]); g
          Graph on 5 vertices
          sage: g = Graph([(0, 2, '0'), (0, 2, '1'), (3, 3, '2')],
          ....:           loops=True, multiedges=True)
          sage: g.loops()
          [(3, 3, '2')]

    #. A NetworkX MultiGraph::

          sage: import networkx                                                         # needs networkx
          sage: g = networkx.MultiGraph({0:[1,2,3], 2:[4]})                             # needs networkx
          sage: Graph(g)                                                                # needs networkx
          Multi-graph on 5 vertices

    #. A NetworkX graph::

           sage: import networkx                                                        # needs networkx
           sage: g = networkx.Graph({0:[1,2,3], 2:[4]})                                 # needs networkx
           sage: DiGraph(g)                                                             # needs networkx
           Digraph on 5 vertices

    #. An igraph Graph (see also
       :meth:`~sage.graphs.generic_graph.GenericGraph.igraph_graph`)::

           sage: import igraph                       # optional - python_igraph
           sage: g = igraph.Graph([(0, 1), (0, 2)])  # optional - python_igraph
           sage: Graph(g)                            # optional - python_igraph
           Graph on 3 vertices

       If ``vertex_labels`` is ``True``, the names of the vertices are given by
       the vertex attribute ``'name'``, if available::

           sage: # optional - python_igraph
           sage: g = igraph.Graph([(0,1),(0,2)], vertex_attrs={'name':['a','b','c']})
           sage: Graph(g).vertices(sort=True)
           ['a', 'b', 'c']
           sage: g = igraph.Graph([(0,1),(0,2)], vertex_attrs={'label':['a','b','c']})
           sage: Graph(g).vertices(sort=True)
           [0, 1, 2]

       If the igraph Graph has edge attributes, they are used as edge labels::

           sage: g = igraph.Graph([(0, 1), (0, 2)],                             # optional - python_igraph
           ....:                  edge_attrs={'name': ['a', 'b'], 'weight': [1, 3]})
           sage: Graph(g).edges(sort=True)                                      # optional - python_igraph
           [(0, 1, {'name': 'a', 'weight': 1}), (0, 2, {'name': 'b', 'weight': 3})]


    When defining an undirected graph from a function ``f``, it is *very*
    important that ``f`` be symmetric. If it is not, anything can happen::

        sage: f_sym = lambda x,y: abs(x-y) == 1
        sage: f_nonsym = lambda x,y: (x-y) == 1
        sage: G_sym = Graph([[4,6,1,5,3,7,2,0], f_sym])
        sage: G_sym.is_isomorphic(graphs.PathGraph(8))
        True
        sage: G_nonsym = Graph([[4,6,1,5,3,7,2,0], f_nonsym])
        sage: G_nonsym.size()
        4
        sage: G_nonsym.is_isomorphic(G_sym)
        False

    By default, graphs are mutable and can thus not be used as a dictionary
    key::

          sage: G = graphs.PetersenGraph()
          sage: {G:1}[G]
          Traceback (most recent call last):
          ...
          TypeError: ...This graph is mutable, and thus not hashable...

    When providing the optional arguments ``data_structure="static_sparse"`` or
    ``immutable=True`` (both mean the same), then an immutable graph results::

          sage: G_imm = Graph(G, immutable=True)
          sage: H_imm = Graph(G, data_structure='static_sparse')
          sage: G_imm == H_imm == G
          True
          sage: {G_imm:1}[H_imm]
          1

    TESTS::

        sage: Graph(4, format='HeyHeyHey')
        Traceback (most recent call last):
        ...
        ValueError: Unknown input format 'HeyHeyHey'

        sage: Graph(igraph.Graph(directed=True))  # optional - python_igraph
        Traceback (most recent call last):
        ...
        ValueError: An *undirected* igraph graph was expected.
        To build a directed graph, call the DiGraph constructor.

        sage: # needs sage.modules
        sage: m = matrix([[0, -1], [-1, 0]])
        sage: Graph(m, format='seidel_adjacency_matrix')
        Graph on 2 vertices
        sage: m[0,1] = 1
        sage: Graph(m, format='seidel_adjacency_matrix')
        Traceback (most recent call last):
        ...
        ValueError: the adjacency matrix of a Seidel graph must be symmetric

        sage: m[0,1] = -1; m[1,1] = 1                                                   # needs sage.modules
        sage: Graph(m, format='seidel_adjacency_matrix')                                # needs sage.modules
        Traceback (most recent call last):
        ...
        ValueError: the adjacency matrix of a Seidel graph must have 0s on the main diagonal

    From a list of vertices and a list of edges::

        sage: G = Graph([[1,2,3], [(1,2)]]); G
        Graph on 3 vertices
        sage: G.edges(sort=True)
        [(1, 2, None)]

    Check that :issue:`27505` is fixed::

        sage: Graph(Graph().networkx_graph(), weighted=None, format='NX')               # needs networkx
        Graph on 0 vertices
    """
    _directed = False

    def __init__(self, data=None, pos=None, loops=None, format=None,
                 weighted=None, data_structure='sparse',
                 vertex_labels=True, name=None,
                 multiedges=None, convert_empty_dict_labels_to_None=None,
                 sparse=True, immutable=False, hash_labels=None):
        """
        TESTS::

            sage: G = Graph()
            sage: loads(dumps(G)) == G
            True
            sage: a = matrix(2,2,[1,0,0,1])                                             # needs sage.modules
            sage: Graph(a).adjacency_matrix() == a                                      # needs sage.modules
            True

            sage: a = matrix(2,2,[2,0,0,1])                                             # needs sage.modules
            sage: Graph(a,sparse=True).adjacency_matrix() == a                          # needs sage.modules
            True

        The positions are copied when the graph is built from another graph ::

            sage: g = graphs.PetersenGraph()
            sage: h = Graph(g)
            sage: g.get_pos() == h.get_pos()
            True

        The position dictionary is not the input one (:issue:`22424`)::

            sage: my_pos = {0:(0,0), 1:(1,1)}
            sage: G = Graph([[0,1], [(0,1)]], pos=my_pos)
            sage: my_pos == G._pos
            True
            sage: my_pos is G._pos
            False

        Or from a DiGraph ::

            sage: d = DiGraph(g)
            sage: h = Graph(d)
            sage: g.get_pos() == h.get_pos()
            True

        Loops are not counted as multiedges (see :issue:`11693`) and edges are
        not counted twice ::

            sage: Graph({1:[1]}).n_edges()
            1
            sage: Graph({1:[2,2]}).n_edges()
            2

        An empty list or dictionary defines a simple graph
        (:issue:`10441` and :issue:`12910`)::

            sage: Graph([])
            Graph on 0 vertices
            sage: Graph({})
            Graph on 0 vertices
            sage: # not "Multi-graph on 0 vertices"

        Verify that the int format works as expected (:issue:`12557`)::

            sage: Graph(2).adjacency_matrix()                                           # needs sage.modules
            [0 0]
            [0 0]
            sage: Graph(3) == Graph(3, format='int')
            True

        Problem with weighted adjacency matrix (:issue:`13919`)::

            sage: B = {0:{1:2,2:5,3:4},1:{2:2,4:7},2:{3:1,4:4,5:3},3:{5:4},4:{5:1,6:5},5:{6:7}}
            sage: grafo3 = Graph(B, weighted=True)
            sage: matad = grafo3.weighted_adjacency_matrix()                            # needs sage.modules
            sage: grafo4 = Graph(matad, format='adjacency_matrix', weighted=True)       # needs sage.modules
            sage: grafo4.shortest_path(0, 6, by_weight=True)                            # needs sage.modules
            [0, 1, 2, 5, 4, 6]

        Graphs returned when setting ``immutable=False`` are mutable::

            sage: g = graphs.PetersenGraph()
            sage: g = Graph(g.edges(sort=True), immutable=False)
            sage: g.add_edge("Hey", "Heyyyyyyy")

        And their name is set::

            sage: g = graphs.PetersenGraph()
            sage: Graph(g, immutable=True)
            Petersen graph: Graph on 10 vertices

        Check error messages for graphs built from incidence matrices (see
        :issue:`18440`)::

            sage: Graph(matrix([[-1, 1, 0],[1, 0, 0]]))                                 # needs sage.modules
            Traceback (most recent call last):
            ...
            ValueError: column 1 of the (oriented) incidence matrix
            contains only one nonzero value
            sage: Graph(matrix([[1,1],[1,1],[1,0]]))                                    # needs sage.modules
            Traceback (most recent call last):
            ...
            ValueError: there must be one or two nonzero entries per column
            in an incidence matrix, got entries [1, 1, 1] in column 0
            sage: Graph(matrix([[3,1,1],[0,1,1]]))                                      # needs sage.modules
            Traceback (most recent call last):
            ...
            ValueError: each column of a non-oriented incidence matrix
            must sum to 2, but column 0 does not

        Vertex labels are retained in the graph (:issue:`14708`)::

            sage: g = Graph()
            sage: g.add_vertex(0)
            sage: g.set_vertex(0, 'foo')
            sage: g.get_vertices()
            {0: 'foo'}
            sage: Graph(g).get_vertices()
            {0: 'foo'}
        """
        GenericGraph.__init__(self)

        from sage.structure.element import Matrix

        if sparse is False:
            if data_structure != "sparse":
                raise ValueError("The 'sparse' argument is an alias for "
                                 "'data_structure'. Please do not define both.")
            data_structure = "dense"

        if multiedges or weighted:
            if data_structure == "dense":
                raise RuntimeError("Multiedge and weighted c_graphs must be sparse.")
        if immutable:
            data_structure = 'static_sparse'

        # If the data structure is static_sparse, we first build a graph
        # using the sparse data structure, then re-encode the resulting graph
        # as a static sparse graph.
        from sage.graphs.base.sparse_graph import SparseGraphBackend
        from sage.graphs.base.dense_graph import DenseGraphBackend
        if data_structure in ["sparse", "static_sparse"]:
            CGB = SparseGraphBackend
        elif data_structure == "dense":
            CGB = DenseGraphBackend
        else:
            raise ValueError("data_structure must be equal to 'sparse', "
                             "'static_sparse' or 'dense'")
        self._backend = CGB(0, directed=False)

        if format is None and isinstance(data, str):
            if data.startswith(">>graph6<<"):
                data = data[10:]
                format = 'graph6'
            elif data.startswith(">>sparse6<<"):
                data = data[11:]
                format = 'sparse6'
            elif data[0] == ':':
                format = 'sparse6'
            else:
                format = 'graph6'
        if format is None and isinstance(data, Matrix):
            if data.is_symmetric():
                format = 'adjacency_matrix'
            else:
                format = 'incidence_matrix'
        if format is None and isinstance(data, Graph):
            format = 'Graph'
        from sage.graphs.digraph import DiGraph
        if format is None and isinstance(data, DiGraph):
            data = data.to_undirected()
            format = 'Graph'
        if (format is None and
                isinstance(data, list) and
                len(data) >= 2 and
                callable(data[1])):
            format = 'rule'

        if (format is None and
                isinstance(data, list) and
                len(data) == 2 and
                isinstance(data[0], list) and    # a list of two lists, the second of
                ((isinstance(data[1], list) and  # which contains iterables (the edges)
                 (not data[1] or callable(getattr(data[1][0], "__iter__", None)))) or
                 isinstance(data[1], EdgesView))):
            format = "vertices_and_edges"

        if format is None and isinstance(data, dict):
            if not data:
                format = 'dict_of_dicts'
            else:
                val = next(iter(data.values()))
                if isinstance(val, (list, EdgesView)):
                    format = 'dict_of_lists'
                elif isinstance(val, dict):
                    format = 'dict_of_dicts'
        if format is None and hasattr(data, 'adj'):
            # the input is a networkx (Multi)(Di)Graph
            format = 'NX'

        if (format is None and
                hasattr(data, 'vcount') and
                hasattr(data, 'get_edgelist')):
            try:
                import igraph
            except ImportError:
                raise ImportError("The data seems to be a igraph object, but "
                                  "igraph is not installed in Sage. To install "
                                  "it, run 'sage -i python_igraph'")
            if format is None and isinstance(data, igraph.Graph):
                format = 'igraph'
        if format is None and isinstance(data, (int, Integer)):
            format = 'int'
        if format is None and data is None:
            format = 'int'
            data = 0

        # Input is a list of edges or an EdgesView
        if format is None and isinstance(data, (list, EdgesView)):
            format = "list_of_edges"
            if weighted is None:
                weighted = False

        if format is None:
            raise ValueError("This input cannot be turned into a graph")

        if format == 'weighted_adjacency_matrix':
            if weighted is False:
                raise ValueError("Format was weighted_adjacency_matrix but weighted was False.")
            if weighted is None:
                weighted = True
            if multiedges is None:
                multiedges = False
            format = 'adjacency_matrix'

        # At this point, 'format' has been set. We build the graph

        if format == 'graph6':
            if weighted is None:
                weighted = False
            self.allow_loops(loops if loops else False, check=False)
            self.allow_multiple_edges(multiedges if multiedges else False, check=False)
            from .graph_input import from_graph6
            from_graph6(self, data)

        elif format == 'sparse6':
            if weighted is None:
                weighted = False
            self.allow_loops(False if loops is False else True, check=False)
            self.allow_multiple_edges(False if multiedges is False else True, check=False)
            from .graph_input import from_sparse6
            from_sparse6(self, data)

        elif format == 'adjacency_matrix':
            from .graph_input import from_adjacency_matrix
            from_adjacency_matrix(self, data, loops=loops, multiedges=multiedges, weighted=weighted)

        elif format == 'incidence_matrix':
            from .graph_input import from_incidence_matrix
            from_incidence_matrix(self, data, loops=loops, multiedges=multiedges, weighted=weighted)

        elif format == 'seidel_adjacency_matrix':
            weighted = False
            self.allow_loops(False)
            self.allow_multiple_edges(False)
            from .graph_input import from_seidel_adjacency_matrix
            from_seidel_adjacency_matrix(self, data)
        elif format == 'Graph':
            if loops is None:
                loops = data.allows_loops()
            if multiedges is None:
                multiedges = data.allows_multiple_edges()
            if weighted is None:
                weighted = data.weighted()
            self.allow_loops(loops, check=False)
            self.allow_multiple_edges(multiedges, check=False)
            if data.get_pos() is not None:
                pos = data.get_pos()
            self.name(data.name())
            self.set_vertices(data.get_vertices())
            data._backend.subgraph_given_vertices(self._backend, data)

        elif format == 'NX':
            from sage.graphs.graph_input import from_networkx_graph
            from_networkx_graph(self, data,
                                weighted=weighted, multiedges=multiedges, loops=loops,
                                convert_empty_dict_labels_to_None=convert_empty_dict_labels_to_None)
            if weighted is None:
                weighted = self.allows_multiple_edges()

        elif format == 'igraph':
            if data.is_directed():
                raise ValueError("An *undirected* igraph graph was expected. "
                                 "To build a directed graph, call the DiGraph "
                                 "constructor.")

            self.add_vertices(range(data.vcount()))
            self.add_edges((e.source, e.target, e.attributes()) for e in data.es())

            if vertex_labels and 'name' in data.vertex_attributes():
                vs = data.vs()
                self.relabel({v: vs[v]['name'] for v in self})

        elif format == 'rule':
            f = data[1]
            verts = data[0]
            if loops is None:
                loops = any(f(v, v) for v in verts)
            if weighted is None:
                weighted = False
            self.allow_loops(loops, check=False)
            self.allow_multiple_edges(bool(multiedges), check=False)
            self.add_vertices(verts)
            self.add_edges(e for e in itertools.combinations(verts, 2) if f(*e))
            if loops:
                self.add_edges((v, v) for v in verts if f(v, v))

        elif format == "vertices_and_edges":
            self.allow_multiple_edges(bool(multiedges), check=False)
            self.allow_loops(bool(loops), check=False)
            self.add_vertices(data[0])
            self.add_edges(data[1])

        elif format == 'dict_of_dicts':
            from .graph_input import from_dict_of_dicts
            from_dict_of_dicts(self, data, loops=loops, multiedges=multiedges, weighted=weighted,
                               convert_empty_dict_labels_to_None=(False if convert_empty_dict_labels_to_None is None
                                                                  else convert_empty_dict_labels_to_None))

        elif format == 'dict_of_lists':
            from .graph_input import from_dict_of_lists
            from_dict_of_lists(self, data, loops=loops, multiedges=multiedges, weighted=weighted)

        elif format == 'int':
            self.allow_loops(loops if loops else False, check=False)
            self.allow_multiple_edges(multiedges if multiedges else False, check=False)
            if data < 0:
                raise ValueError("The number of vertices cannot be strictly negative!")
            if data:
                self.add_vertices(range(data))

        elif format == 'list_of_edges':
            self.allow_multiple_edges(bool(multiedges),
                                      check=False)
            self.allow_loops(bool(loops), check=False)
            self.add_edges(data)
        else:
            raise ValueError("Unknown input format '{}'".format(format))

        if weighted is None:
            weighted = False
        self._weighted = getattr(self, '_weighted', weighted)

        if hash_labels is None and hasattr(data, '_hash_labels'):
            hash_labels = data._hash_labels
        self._hash_labels = hash_labels

        self._pos = copy(pos)

        if format != 'Graph' or name is not None:
            self.name(name)

        if data_structure == "static_sparse":
            from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            ib = StaticSparseBackend(self,
                                     loops=self.allows_loops(),
                                     multiedges=self.allows_multiple_edges())
            self._backend = ib
            self._immutable = True

    # Formats

    @doc_index("Basic methods")
    def graph6_string(self):
        r"""
        Return the graph6 representation of the graph as an ASCII string.

        This is only valid for simple (no loops, no multiple edges) graphs
        on at most `2^{18}-1=262143` vertices.

        .. NOTE::

            As the graph6 format only handles graphs with vertex set
            `\{0,...,n-1\}`, a :meth:`relabelled copy
            <sage.graphs.generic_graph.GenericGraph.relabel>` will
            be encoded, if necessary.

        .. SEEALSO::

            * :meth:`~sage.graphs.digraph.DiGraph.dig6_string` --
              a similar string format for directed graphs

        EXAMPLES::

            sage: G = graphs.KrackhardtKiteGraph()
            sage: G.graph6_string()
            'IvUqwK@?G'

        TESTS::

            sage: Graph().graph6_string()
            '?'
        """
        n = self.order()
        if n > 262143:
            raise ValueError('graph6 format supports graphs on 0 to 262143 vertices only.')
        elif self.has_loops() or self.has_multiple_edges():
            raise ValueError('graph6 format supports only simple graphs (no loops, no multiple edges)')
        return generic_graph_pyx.small_integer_to_graph6(n) + generic_graph_pyx.binary_string_to_graph6(self._bit_vector())

    @doc_index("Basic methods")
    def sparse6_string(self):
        r"""
        Return the sparse6 representation of the graph as an ASCII string.

        Only valid for undirected graphs on 0 to 262143 vertices, but loops
        and multiple edges are permitted.

        .. NOTE::

            As the sparse6 format only handles graphs whose vertex set is
            `\{0,...,n-1\}`, a :meth:`relabelled copy
            <sage.graphs.generic_graph.GenericGraph.relabel>` of your graph will
            be encoded if necessary.

        EXAMPLES::

            sage: G = graphs.BullGraph()
            sage: G.sparse6_string()
            ':Da@en'

        ::

            sage: G = Graph(loops=True, multiedges=True, data_structure='sparse')
            sage: Graph(':?', data_structure='sparse') == G
            True

        TESTS::

            sage: G = Graph()
            sage: G.sparse6_string()
            ':?'

        Check that :issue:`18445` is fixed::

            sage: Graph(graphs.KneserGraph(5,2).sparse6_string()).size()
            15

        Graphs with 1 vertex are correctly handled (:issue:`24923`)::

            sage: Graph([(0, 0)], loops=True).sparse6_string()
            ':@^'
            sage: G = Graph(_)
            sage: G.order(), G.size()
            (1, 1)
            sage: Graph([(0, 0), (0, 0)], loops=True, multiedges=True).sparse6_string()
            ':@N'
            sage: H = Graph(_)
            sage: H.order(), H.size()
            (1, 2)

        Sparse6 encoding of canonical graph is unique (:issue:`31026`)::

            sage: G = Graph([(0,1),(1,2),(2,3),(3,0),(0,2)])
            sage: H = Graph([(0,1),(1,2),(2,3),(3,0),(1,3)])
            sage: G == H
            False
            sage: G.is_isomorphic(H)
            True
            sage: G.sparse6_string() == H.sparse6_string()
            False
            sage: G_ = G.canonical_label()
            sage: H_ = H.canonical_label()
            sage: G_ == H_
            True
            sage: G_.sparse6_string() == H_.sparse6_string()
            True

        The method can handle vertices with different types (:issue:`31026`)::

            sage: G = Graph([(1, 'a')])
            sage: H = Graph(G.sparse6_string())
            sage: G.is_isomorphic(H)
            True
            sage: set(G) == set(H)
            False
        """
        n = self.order()
        if not n:
            return ':?'
        if n > 262143:
            raise ValueError('sparse6 format supports graphs on 0 to 262143 vertices only.')
        if n == 1:
            s = '0' * self.size()
        else:
            try:
                V = sorted(self)
            except TypeError:
                V = self
            v_to_int = {v: i for i, v in enumerate(V)}
            edges = [sorted((v_to_int[u], v_to_int[v])) for u, v in self.edge_iterator(labels=False)]
            edges.sort(key=lambda e: (e[1], e[0]))  # reverse lexicographic order

            # encode bit vector
            k = int((ZZ(n) - 1).nbits())
            v = 0
            i = 0
            m = 0
            s = ''
            while m < len(edges):
                if edges[m][1] > v + 1:
                    sp = generic_graph_pyx.int_to_binary_string(edges[m][1])
                    sp = '0'*(k-len(sp)) + sp
                    s += '1' + sp
                    v = edges[m][1]
                elif edges[m][1] == v + 1:
                    sp = generic_graph_pyx.int_to_binary_string(edges[m][0])
                    sp = '0'*(k-len(sp)) + sp
                    s += '1' + sp
                    v += 1
                    m += 1
                else:
                    sp = generic_graph_pyx.int_to_binary_string(edges[m][0])
                    sp = '0'*(k-len(sp)) + sp
                    s += '0' + sp
                    m += 1

        # encode s as a 6-string, as in R(x), but padding with 1's
        # pad on the right to make a multiple of 6
        s = s + ('1' * ((6 - len(s)) % 6))

        # split into groups of 6, and convert numbers to decimal, adding 63
        six_bits = ''
        for i in range(0, len(s), 6):
            six_bits += chr(int(s[i:i+6], 2) + 63)
        return ':' + generic_graph_pyx.small_integer_to_graph6(n) + six_bits

    # Attributes

    @doc_index("Basic methods")
    def is_directed(self):
        """
        Since graph is undirected, returns False.

        EXAMPLES::

            sage: Graph().is_directed()
            False
        """
        return False

    # Properties

    @doc_index("Graph properties")
    def is_tree(self, certificate=False, output='vertex'):
        r"""
        Test if the graph is a tree.

        The empty graph is defined to be not a tree.

        INPUT:

        - ``certificate`` -- boolean (default: ``False``); whether to return a
          certificate. The method only returns boolean answers when
          ``certificate = False`` (default). When it is set to ``True``, it
          either answers ``(True, None)`` when the graph is a tree or ``(False,
          cycle)`` when it contains a cycle. It returns ``(False, None)`` when
          the graph is empty or not connected.

        - ``output`` -- either ``'vertex'`` (default) or ``'edge'``; whether the
          certificate is given as a list of vertices (``output = 'vertex'``) or
          a list of edges (``output = 'edge'``).

        When the certificate cycle is given as a list of edges, the edges are
        given as `(v_i, v_{i+1}, l)` where `v_1, v_2, \dots, v_n` are the
        vertices of the cycles (in their cyclic order).

        EXAMPLES::

            sage: all(T.is_tree() for T in graphs.trees(15))
            True

        With certificates::

            sage: g = graphs.RandomTree(30)
            sage: g.is_tree(certificate=True)
            (True, None)
            sage: g.add_edge(10,-1)
            sage: g.add_edge(11,-1)
            sage: isit, cycle = g.is_tree(certificate=True)
            sage: isit
            False
            sage: -1 in cycle
            True

        One can also ask for the certificate as a list of edges::

            sage: g = graphs.CycleGraph(4)
            sage: g.is_tree(certificate=True, output='edge')
            (False, [(3, 2, None), (2, 1, None), (1, 0, None), (0, 3, None)])

        This is useful for graphs with multiple edges::

            sage: G = Graph([(1, 2, 'a'), (1, 2, 'b')], multiedges=True)
            sage: G.is_tree(certificate=True)
            (False, [1, 2])
            sage: G.is_tree(certificate=True, output='edge')
            (False, [(1, 2, 'b'), (2, 1, 'a')])

        TESTS:

        :issue:`14434` is fixed::

            sage: g = Graph({0:[1,4,5],3:[4,8,9],4:[9],5:[7,8],7:[9]})
            sage: _,cycle = g.is_tree(certificate=True)
            sage: g.size()
            10
            sage: g.add_cycle(cycle)
            sage: g.size()
            10

        The empty graph::

            sage: graphs.EmptyGraph().is_tree()
            False
            sage: graphs.EmptyGraph().is_tree(certificate=True)
            (False, None)

        :issue:`22912` is fixed::

            sage: G = Graph([(0,0), (0,1)], loops=True)
            sage: G.is_tree(certificate=True)
            (False, [0])
            sage: G.is_tree(certificate=True, output='edge')
            (False, [(0, 0, None)])

        Case of edges with incomparable types (see :issue:`35903`)::

            sage: G = Graph(multiedges=True)
            sage: G.add_cycle(['A', 1, 2, 3])
            sage: G.add_cycle(['A', 1, 2, 3])
            sage: G.is_tree(certificate=True, output='vertex')
            (False, ['A', 1])
            sage: G.is_tree(certificate=True, output='edge')
            (False, [('A', 1, None), (1, 'A', None)])
        """
        if output not in ['vertex', 'edge']:
            raise ValueError('output must be either vertex or edge')

        if not self.order() or not self.is_connected():
            return (False, None) if certificate else False

        if certificate:
            if self.order() == self.size() + 1:
                return (True, None)

            if self.allows_loops():
                L = self.loop_edges() if output == 'edge' else self.loop_vertices()
                if L:
                    return False, L[:1]

            if self.has_multiple_edges():
                multiple_edges = self.multiple_edges(sort=False)
                if output == 'vertex':
                    return (False, list(multiple_edges[0][:2]))
                # Search for 2 edges between u and v.
                # We do this way to handle the case of edges with incomparable
                # types
                u1, v1, w1 = multiple_edges[0]
                for u2, v2, w2 in multiple_edges[1:]:
                    if u1 == u2 and v1 == v2:
                        return (False, [(u1, v1, w1), (v2, u2, w2)])
                    elif u1 == v2 and v1 == u2:
                        return (False, [(u1, v1, w1), (u2, v2, w2)])

            if output == 'edge':
                if self.allows_multiple_edges():
                    def vertices_to_edges(x):
                        return [(u[0], u[1], self.edge_label(u[0], u[1])[0])
                                for u in zip(x, x[1:] + [x[0]])]
                else:
                    def vertices_to_edges(x):
                        return [(u[0], u[1], self.edge_label(u[0], u[1]))
                                for u in zip(x, x[1:] + [x[0]])]

            # This code is a depth-first search that looks for a cycle in the
            # graph. We *know* it exists as there are too many edges around.
            seen = {}
            u = next(self.vertex_iterator())
            seen[u] = u
            stack = [(u, v) for v in self.neighbor_iterator(u)]
            while stack:
                u, v = stack.pop()
                if v in seen:
                    continue
                for w in self.neighbor_iterator(v):
                    if u == w:
                        continue
                    elif w in seen:
                        cycle = [w, v]
                        while u != w:
                            cycle.append(u)
                            u = seen[u]
                        cycle.reverse()
                        if output == 'vertex':
                            return (False, cycle)
                        return (False, vertices_to_edges(cycle))
                    else:
                        stack.append((v, w))
                seen[v] = u

        return self.order() == self.size() + 1

    @doc_index("Graph properties")
    def is_forest(self, certificate=False, output='vertex'):
        """
        Test if the graph is a forest, i.e. a disjoint union of trees.

        INPUT:

        - ``certificate`` -- boolean (default: ``False``); whether to return a
          certificate. The method only returns boolean answers when
          ``certificate = False`` (default). When it is set to ``True``, it
          either answers ``(True, None)`` when the graph is a forest or
          ``(False, cycle)`` when it contains a cycle.

        - ``output`` -- either ``'vertex'`` (default) or ``'edge'``; whether the
          certificate is given as a list of vertices (``output = 'vertex'``) or
          a list of edges (``output = 'edge'``).

        EXAMPLES::

            sage: seven_acre_wood = sum(graphs.trees(7), Graph())
            sage: seven_acre_wood.is_forest()
            True

        With certificates::

            sage: g = graphs.RandomTree(30)
            sage: g.is_forest(certificate=True)
            (True, None)
            sage: (2*g + graphs.PetersenGraph() + g).is_forest(certificate=True)
            (False, [64, 69, 67, 65, 60])
        """
        connected_components = self.connected_components(sort=False)
        number_of_connected_components = len(connected_components)
        isit = (self.order() ==
                self.size() + number_of_connected_components)

        if not certificate:
            return isit
        if isit:
            return (True, None)

        # The graph contains a cycle, and the user wants to see it.
        if number_of_connected_components == 1:
            return self.is_tree(certificate=True, output=output)
        # We try to find a cycle in each connected component
        for cc in connected_components:
            isit, cycle = self.subgraph(cc).is_tree(certificate=True, output=output)
            if not isit:
                return (False, cycle)

    @doc_index("Graph properties")
    def is_cactus(self):
        """
        Check whether the graph is cactus graph.

        A graph is called *cactus graph* if it is connected and every pair of
        simple cycles have at most one common vertex.

        There are other definitions, see the :wikipedia:`Cactus_graph`.

        EXAMPLES::

            sage: g = Graph({1: [2], 2: [3, 4], 3: [4, 5, 6, 7], 8: [3, 5], 9: [6, 7]})
            sage: g.is_cactus()
            True

            sage: c6 = graphs.CycleGraph(6)
            sage: naphthalene = c6 + c6
            sage: naphthalene.is_cactus()  # Not connected
            False
            sage: naphthalene.merge_vertices([0, 6])
            sage: naphthalene.is_cactus()
            True
            sage: naphthalene.merge_vertices([1, 7])
            sage: naphthalene.is_cactus()
            False

        TESTS::

            sage: all(graphs.PathGraph(i).is_cactus() for i in range(5))
            True

            sage: Graph('Fli@?').is_cactus()
            False

            sage: Graph('BG').is_cactus()
            False
            sage: (Graph(0).is_cactus(), Graph(1).is_cactus())
            (True, True)
            sage: (Graph(2).is_cactus(), Graph(3).is_cactus())
            (False, False)

        Test a graph that is not outerplanar, see :issue:`24480`::

            sage: graphs.Balaban10Cage().is_cactus()
            False
        """
        self._scream_if_not_simple()

        if not self.is_connected():
            return False

        # Special cases
        if self.order() < 4:
            return True

        # trees are cacti
        if self.order() == self.size() + 1:
            return True

        if self.size() > 3 * (self.order() - 1) / 2:
            return False

        # Every cactus graph is outerplanar
        if not self.is_circular_planar():
            return False

        # the number of faces is 1 plus the number of blocks of order > 2
        B = self.blocks_and_cut_vertices()[0]
        return len(self.faces()) == sum(1 for b in B if len(b) > 2) + 1

    @doc_index("Graph properties")
    def is_block_graph(self):
        r"""
        Return whether this graph is a block graph.

        A block graph is a connected graph in which every biconnected component
        (block) is a clique.

        .. SEEALSO::

            - :wikipedia:`Block_graph` for more details on these graphs
            - :meth:`~sage.graphs.graph_generators.GraphGenerators.RandomBlockGraph`
              -- generator of random block graphs
            - :meth:`~sage.graphs.generic_graph.GenericGraph.blocks_and_cut_vertices`
            - :meth:`~sage.graphs.generic_graph.GenericGraph.blocks_and_cuts_tree`

        EXAMPLES::

            sage: G = graphs.RandomBlockGraph(6, 2, kmax=4)
            sage: G.is_block_graph()
            True
            sage: from sage.graphs.isgci import graph_classes
            sage: G in graph_classes.Block
            True
            sage: graphs.CompleteGraph(4).is_block_graph()
            True
            sage: graphs.RandomTree(6).is_block_graph()
            True
            sage: graphs.PetersenGraph().is_block_graph()
            False
            sage: Graph(4).is_block_graph()
            False
        """
        if not self.is_connected():
            return False
        if self.is_clique():
            return True

        B, C = self.blocks_and_cut_vertices()
        return all(self.is_clique(vertices=block) for block in B)

    @doc_index("Graph properties")
    def is_cograph(self):
        """
        Check whether the graph is cograph.

        A cograph is defined recursively: the single-vertex graph is
        cograph, complement of cograph is cograph, and disjoint union
        of two cographs is cograph. There are many other
        characterizations, see the :wikipedia:`Cograph`.

        EXAMPLES::

            sage: graphs.HouseXGraph().is_cograph()
            True
            sage: graphs.HouseGraph().is_cograph()                                      # needs sage.modules
            False

        .. TODO::

            Implement faster recognition algorithm, as for instance
            the linear time recognition algorithm using LexBFS proposed
            in [Bre2008]_.

        TESTS::

            sage: [graphs.PathGraph(i).is_cograph() for i in range(6)]                  # needs sage.modules
            [True, True, True, True, False, False]
            sage: graphs.CycleGraph(5).is_cograph()  # Self-complemented                # needs sage.modules
            False
        """
        # A cograph has no 4-vertex path as an induced subgraph.
        # We will first try to "decompose" graph by complements and
        # split to connected components, and use fairly slow
        # subgraph search if that fails.
        self._scream_if_not_simple()
        if self.order() < 4:
            return True
        if self.density()*2 > 1:
            return self.complement().is_cograph()
        if not self.is_connected():
            return all(part.is_cograph() for part in self.connected_components_subgraphs())
        P4 = Graph({0: [1], 1: [2], 2: [3]})
        return self.subgraph_search(P4, induced=True) is None

    @doc_index("Graph properties")
    def is_apex(self):
        r"""
        Test if the graph is apex.

        A graph is apex if it can be made planar by the removal of a single
        vertex. The deleted vertex is called ``an apex`` of the graph, and a
        graph may have more than one apex. For instance, in the minimal
        nonplanar graphs `K_5` or `K_{3,3}`, every vertex is an apex. The apex
        graphs include graphs that are themselves planar, in which case again
        every vertex is an apex. The null graph is also counted as an apex graph
        even though it has no vertex to remove.  If the graph is not connected,
        we say that it is apex if it has at most one non planar connected
        component and that this component is apex.  See the :wikipedia:`Apex_graph`
        for more information.

        .. SEEALSO::

          - :meth:`~Graph.apex_vertices`
          - :meth:`~sage.graphs.generic_graph.GenericGraph.is_planar`

        EXAMPLES:

        `K_5` and `K_{3,3}` are apex graphs, and each of their vertices is an
        apex::

            sage: G = graphs.CompleteGraph(5)
            sage: G.is_apex()
            True
            sage: G = graphs.CompleteBipartiteGraph(3,3)
            sage: G.is_apex()
            True

        The Petersen graph is not apex::

            sage: G = graphs.PetersenGraph()
            sage: G.is_apex()
            False

        A graph is apex if all its connected components are apex, but at most
        one is not planar::

            sage: M = graphs.Grid2dGraph(3,3)
            sage: K5 = graphs.CompleteGraph(5)
            sage: (M+K5).is_apex()
            True
            sage: (M+K5+K5).is_apex()
            False

        TESTS:

        The null graph is apex::

            sage: G = Graph()
            sage: G.is_apex()
            True

        The graph might be mutable or immutable::

            sage: G = Graph(M+K5, immutable=True)
            sage: G.is_apex()
            True
        """
        # Easy cases: null graph, subgraphs of K_5 and K_3,3
        if self.order() <= 5 or (self.order() <= 6 and self.is_bipartite()):
            return True

        return len(self.apex_vertices(k=1)) > 0

    @doc_index("Graph properties")
    def apex_vertices(self, k=None):
        r"""
        Return the list of apex vertices.

        A graph is apex if it can be made planar by the removal of a single
        vertex. The deleted vertex is called ``an apex`` of the graph, and a
        graph may have more than one apex. For instance, in the minimal
        nonplanar graphs `K_5` or `K_{3,3}`, every vertex is an apex. The apex
        graphs include graphs that are themselves planar, in which case again
        every vertex is an apex. The null graph is also counted as an apex graph
        even though it has no vertex to remove.  If the graph is not connected,
        we say that it is apex if it has at most one non planar connected
        component and that this component is apex.  See the
        :wikipedia:`Apex_graph` for more information.

        .. SEEALSO::

          - :meth:`~Graph.is_apex`
          - :meth:`~sage.graphs.generic_graph.GenericGraph.is_planar`

        INPUT:

        - ``k`` -- integer (default: ``None``); when set to ``None``, the method
          returns the list of all apex of the graph, possibly empty if the graph
          is not apex. When set to a positive integer, the method ends as soon
          as `k` apex vertices are found.

        OUTPUT:

        By default, the method returns the list of all apex of the graph. When
        parameter ``k`` is set to a positive integer, the returned list is
        bounded to `k` apex vertices.

        EXAMPLES:

        `K_5` and `K_{3,3}` are apex graphs, and each of their vertices is an
        apex::

            sage: G = graphs.CompleteGraph(5)
            sage: G.apex_vertices()
            [0, 1, 2, 3, 4]
            sage: G = graphs.CompleteBipartiteGraph(3,3)
            sage: G.is_apex()
            True
            sage: G.apex_vertices()
            [0, 1, 2, 3, 4, 5]
            sage: G.apex_vertices(k=3)
            [0, 1, 2]

        A `4\\times 4`-grid is apex and each of its vertices is an apex. When
        adding a universal vertex, the resulting graph is apex and the universal
        vertex is the unique apex vertex ::

            sage: G = graphs.Grid2dGraph(4,4)
            sage: set(G.apex_vertices()) == set(G.vertices(sort=False))
            True
            sage: G.add_edges([('universal',v) for v in G])
            sage: G.apex_vertices()
            ['universal']

        The Petersen graph is not apex::

            sage: G = graphs.PetersenGraph()
            sage: G.apex_vertices()
            []

        A graph is apex if all its connected components are apex, but at most
        one is not planar::

            sage: M = graphs.Grid2dGraph(3,3)
            sage: K5 = graphs.CompleteGraph(5)
            sage: (M+K5).apex_vertices()
            [9, 10, 11, 12, 13]
            sage: (M+K5+K5).apex_vertices()
            []

        Neighbors of an apex of degree 2 are apex::

            sage: G = graphs.Grid2dGraph(5,5)
            sage: v = (666, 666)
            sage: G.add_path([(1, 1), v, (3, 3)])
            sage: G.is_planar()
            False
            sage: G.degree(v)
            2
            sage: sorted(G.apex_vertices())
            [(1, 1), (2, 2), (3, 3), (666, 666)]


        TESTS:

        The null graph is apex although it has no apex vertex::

            sage: G = Graph()
            sage: G.apex_vertices()
            []

        Parameter ``k`` cannot be a negative integer::

            sage: G.apex_vertices(k=-1)
            Traceback (most recent call last):
            ...
            ValueError: parameter k must be a nonnegative integer

        The graph might be mutable or immutable::

            sage: G = Graph(M+K5, immutable=True)
            sage: G.apex_vertices()
            [9, 10, 11, 12, 13]
        """
        if k is None:
            k = self.order()
        elif k < 0:
            raise ValueError("parameter k must be a nonnegative integer")

        # Easy cases: null graph, subgraphs of K_5 and K_3,3
        if self.order() <= 5 or (self.order() <= 6 and self.is_bipartite()):
            it = self.vertex_iterator()
            return [next(it) for _ in range(k)]

        if not self.is_connected():
            # We search for its non planar connected components. If it has more
            # than one such component, the graph is not apex. It is apex if
            # either it has no such component, in which case the graph is
            # planar, or if its unique non planar component is apex.

            P = [H for H in self.connected_components_subgraphs() if not H.is_planar()]
            if not P:  # The graph is planar
                it = self.vertex_iterator()
                return [next(it) for _ in range(k)]
            if len(P) > 1:
                return []

            # We proceed with the non planar component
            if P[0].is_immutable():
                H = Graph(P[0].edges(labels=0, sort=False), immutable=False, loops=False, multiedges=False)
            else:
                H = P[0]

        elif self.is_planar():
            # A planar graph is apex.
            it = self.vertex_iterator()
            return [next(it) for _ in range(k)]

        else:
            # We make a basic copy of the graph since we will modify it
            H = Graph(self.edges(labels=0, sort=False), immutable=False, loops=False, multiedges=False)

        # General case: basic implementation
        #
        # Test for each vertex if its removal makes the graph planar.
        # Obviously, we don't test vertices of degree one. Furthermore, if a
        # vertex of degree 2 is an apex, its neighbors also are. So we start
        # with vertices of degree 2.
        V = {}
        for u in H:
            d = H.degree(u)
            if d > 1:
                if d in V:
                    V[d].append(u)
                else:
                    V[d] = [u]
        apex = set()
        for deg in sorted(V):
            for u in V[deg]:
                if u in apex:  # True if neighbor of an apex of degree 2
                    if deg == 2:
                        # We ensure that its neighbors are known apex
                        apex.update(H.neighbor_iterator(u))
                        if len(apex) >= k:
                            return list(apex)[:k]
                    continue

                E = H.edges_incident(u, labels=0)
                H.delete_vertex(u)
                if H.is_planar():
                    apex.add(u)
                    if deg == 2:
                        # The neighbors of an apex of degree 2 also are
                        apex.update(self.neighbor_iterator(u))

                    if len(apex) >= k:
                        return list(apex)[:k]

                H.add_edges(E)

        return list(apex)

    @doc_index("Graph properties")
    def is_overfull(self):
        r"""
        Test whether the current graph is overfull.

        A graph `G` on `n` vertices and `m` edges is said to be overfull if:

        - `n` is odd

        - It satisfies `2m > (n-1)\Delta(G)`, where `\Delta(G)` denotes the
          maximum degree among all vertices in `G`.

        An overfull graph must have a chromatic index of `\Delta(G)+1`.

        EXAMPLES:

        A complete graph of order `n > 1` is overfull if and only if `n` is
        odd::

            sage: graphs.CompleteGraph(6).is_overfull()
            False
            sage: graphs.CompleteGraph(7).is_overfull()
            True
            sage: graphs.CompleteGraph(1).is_overfull()
            False

        The claw graph is not overfull::

            sage: from sage.graphs.graph_coloring import edge_coloring
            sage: g = graphs.ClawGraph()
            sage: g
            Claw graph: Graph on 4 vertices
            sage: edge_coloring(g, value_only=True)                                     # needs sage.numerical_mip
            3
            sage: g.is_overfull()
            False

        The Holt graph is an example of a overfull graph::

            sage: G = graphs.HoltGraph()
            sage: G.is_overfull()
            True

        Checking that all complete graphs `K_n` for even `0 \leq n \leq 100`
        are not overfull::

            sage: def check_overfull_Kn_even(n):
            ....:     i = 0
            ....:     while i <= n:
            ....:         if graphs.CompleteGraph(i).is_overfull():
            ....:             print("A complete graph of even order cannot be overfull.")
            ....:             return
            ....:         i += 2
            ....:     print("Complete graphs of even order up to %s are not overfull." % n)
            ...
            sage: check_overfull_Kn_even(100)  # long time
            Complete graphs of even order up to 100 are not overfull.

        The null graph, i.e. the graph with no vertices, is not overfull::

            sage: Graph().is_overfull()
            False
            sage: graphs.CompleteGraph(0).is_overfull()
            False

        Checking that all complete graphs `K_n` for odd `1 < n \leq 100`
        are overfull::

            sage: def check_overfull_Kn_odd(n):
            ....:     i = 3
            ....:     while i <= n:
            ....:         if not graphs.CompleteGraph(i).is_overfull():
            ....:             print("A complete graph of odd order > 1 must be overfull.")
            ....:             return
            ....:         i += 2
            ....:     print("Complete graphs of odd order > 1 up to %s are overfull." % n)
            ...
            sage: check_overfull_Kn_odd(100)  # long time
            Complete graphs of odd order > 1 up to 100 are overfull.

        The Petersen Graph, though, is not overfull while
        its chromatic index is `\Delta+1`::

            sage: g = graphs.PetersenGraph()
            sage: g.is_overfull()
            False
            sage: from sage.graphs.graph_coloring import edge_coloring
            sage: max(g.degree()) + 1 ==  edge_coloring(g, value_only=True)             # needs sage.numerical_mip
            True
        """
        # # A possible optimized version. But the gain in speed is very little.
        # return bool(self._backend.n_vertices() & 1) and (  # odd order n
        #     2 * self._backend.n_edges(self._directed) > #2m > \Delta(G)*(n-1)
        #     max(self.degree()) * (self._backend.n_vertices() - 1))
        # unoptimized version
        return (self.order() % 2 == 1) and (
            2 * self.size() > max(self.degree()) * (self.order() - 1))

    @doc_index("Graph properties")
    def is_even_hole_free(self, certificate=False):
        r"""
        Test whether ``self`` contains an induced even hole.

        A Hole is a cycle of length at least 4 (included). It is said to be even
        (resp. odd) if its length is even (resp. odd).

        Even-hole-free graphs always contain a bisimplicial vertex, which
        ensures that their chromatic number is at most twice their clique number
        [ACHRS2008]_.

        INPUT:

        - ``certificate`` -- boolean (default: ``False``); when ``certificate =
          False``, this method only returns ``True`` or ``False``. If
          ``certificate = True``, the subgraph found is returned instead of
          ``False``.

        EXAMPLES:

        Is the Petersen Graph even-hole-free ::

            sage: g = graphs.PetersenGraph()
            sage: g.is_even_hole_free()                                                 # needs sage.modules
            False

        As any chordal graph is hole-free, interval graphs behave the same way::

            sage: g = graphs.RandomIntervalGraph(20)
            sage: g.is_even_hole_free()                                                 # needs sage.modules
            True

        It is clear, though, that a random Bipartite Graph which is not a forest
        has an even hole::

            sage: g = graphs.RandomBipartite(10, 10, .5)                                # needs numpy
            sage: g.is_even_hole_free() and not g.is_forest()                           # needs numpy sage.modules
            False

        We can check the certificate returned is indeed an even cycle::

            sage: if not g.is_forest():                                                 # needs numpy sage.modules
            ....:    cycle = g.is_even_hole_free(certificate=True)
            ....:    if cycle.order() % 2 == 1:
            ....:        print("Error !")
            ....:    if not cycle.is_isomorphic(
            ....:           graphs.CycleGraph(cycle.order())):
            ....:        print("Error !")
            ...
            sage: print("Everything is Fine !")
            Everything is Fine !

        TESTS:

        Bug reported in :issue:`9925`, and fixed by :issue:`9420`::

            sage: g = Graph(':SiBFGaCEF_@CE`DEGH`CEFGaCDGaCDEHaDEF`CEH`ABCDEF',
            ....:           loops=False, multiedges=False)
            sage: g.is_even_hole_free()                                                 # needs sage.modules
            False
            sage: g.is_even_hole_free(certificate=True)                                 # needs sage.modules
            Subgraph of (): Graph on 4 vertices

        Making sure there are no other counter-examples around ::

            sage: t = lambda x: (Graph(x).is_forest() or
            ....:       isinstance(Graph(x).is_even_hole_free(certificate=True), Graph))
            sage: all(t(graphs.RandomBipartite(10, 10, .5)) for i in range(100))        # needs numpy sage.modules
            True
        """
        girth = self.girth()

        if girth > self.order():
            start = 4

        elif not girth % 2:
            if not certificate:
                return False
            start = girth

        else:
            start = girth + 1

        from sage.graphs.generators.basic import CycleGraph

        while start <= self.order():

            subgraph = self.subgraph_search(CycleGraph(start), induced=True)

            if subgraph is not None:
                if certificate:
                    return subgraph
                return False

            start += 2

        return True

    @doc_index("Graph properties")
    def is_odd_hole_free(self, certificate=False):
        r"""
        Test whether ``self`` contains an induced odd hole.

        A Hole is a cycle of length at least 4 (included). It is said to be even
        (resp. odd) if its length is even (resp. odd).

        It is interesting to notice that while it is polynomial to check whether
        a graph has an odd hole or an odd antihole [CCLSV2005]_, it is not known
        whether testing for one of these two cases independently is polynomial
        too.

        INPUT:

        - ``certificate`` -- boolean (default: ``False``); when ``certificate =
          False``, this method only returns ``True`` or ``False``. If
          ``certificate = True``, the subgraph found is returned instead of
          ``False``.

        EXAMPLES:

        Is the Petersen Graph odd-hole-free ::

            sage: g = graphs.PetersenGraph()
            sage: g.is_odd_hole_free()                                                  # needs sage.modules
            False

        Which was to be expected, as its girth is 5 ::

            sage: g.girth()
            5

        We can check the certificate returned is indeed a 5-cycle::

            sage: cycle = g.is_odd_hole_free(certificate=True)                          # needs sage.modules
            sage: cycle.is_isomorphic(graphs.CycleGraph(5))                             # needs sage.modules
            True

        As any chordal graph is hole-free, no interval graph has an odd hole::

            sage: g = graphs.RandomIntervalGraph(20)
            sage: g.is_odd_hole_free()                                                  # needs sage.modules
            True
        """
        girth = self.odd_girth()

        if girth > self.order():
            return True
        if girth == 3:
            start = 5
        else:
            if not certificate:
                return False
            start = girth

        from sage.graphs.generators.basic import CycleGraph

        while start <= self.order():

            subgraph = self.subgraph_search(CycleGraph(start), induced=True)

            if subgraph is not None:
                if certificate:
                    return subgraph
                return False

            start += 2

        return True

    @doc_index("Graph properties")
    def is_triangle_free(self, algorithm='dense_graph', certificate=False):
        r"""
        Check whether ``self`` is triangle-free.

        INPUT:

        - ``algorithm`` -- (default: ``'dense_graph'``) specifies the algorithm
          to use among:

          - ``'matrix'`` -- tests if the trace of the adjacency matrix is
            positive

          - ``'bitset'`` -- encodes adjacencies into bitsets and uses fast
            bitset operations to test if the input graph contains a
            triangle. This method is generally faster than standard matrix
            multiplication.

          - ``'dense_graph'`` -- use the implementation of
            :mod:`sage.graphs.base.static_dense_graph`

        - ``certificate`` -- boolean (default: ``False``); whether to return a
          triangle if one is found. This parameter is ignored when ``algorithm``
          is ``'matrix'``.

        EXAMPLES:

        The Petersen Graph is triangle-free::

            sage: g = graphs.PetersenGraph()
            sage: g.is_triangle_free()
            True

        or a complete Bipartite Graph::

            sage: G = graphs.CompleteBipartiteGraph(5,6)
            sage: G.is_triangle_free(algorithm='matrix')                                # needs sage.modules
            True
            sage: G.is_triangle_free(algorithm='bitset')
            True
            sage: G.is_triangle_free(algorithm='dense_graph')
            True

        a tripartite graph, though, contains many triangles::

            sage: G = (3 * graphs.CompleteGraph(5)).complement()
            sage: G.is_triangle_free(algorithm='matrix')                                # needs sage.modules
            False
            sage: G.is_triangle_free(algorithm='bitset')
            False
            sage: G.is_triangle_free(algorithm='dense_graph')
            False

        Asking for a certificate::

            sage: K4 = graphs.CompleteGraph(4)
            sage: K4.is_triangle_free(algorithm='dense_graph', certificate=True)
            (False, [0, 1, 2])
            sage: K4.is_triangle_free(algorithm='bitset', certificate=True)
            (False, [0, 1, 2])

        TESTS:

        Comparison of algorithms::

            sage: for i in range(10):           # long time                             # needs networkx
            ....:     G = graphs.RandomBarabasiAlbert(50,2)
            ....:     bm = G.is_triangle_free(algorithm='matrix')
            ....:     bb = G.is_triangle_free(algorithm='bitset')
            ....:     bd = G.is_triangle_free(algorithm='dense_graph')
            ....:     if bm != bb or bm != bd:
            ....:        print("That's not good!")

        Asking for an unknown algorithm::

            sage: g.is_triangle_free(algorithm='tip top')
            Traceback (most recent call last):
            ...
            ValueError: Algorithm 'tip top' not yet implemented. Please contribute.

        Check the empty graph::

            sage: graphs.EmptyGraph().is_triangle_free()
            True
        """
        if algorithm == 'dense_graph':
            from sage.graphs.base.static_dense_graph import is_triangle_free
            return is_triangle_free(self, certificate=certificate)

        if algorithm == 'bitset':
            if self.order() < 3:
                return (True, []) if certificate else True
            from sage.data_structures.bitset import Bitset
            N = self.order()
            vertex_to_int = {}
            B = {}
            for i, u in enumerate(self):
                vertex_to_int[u] = i
                B[u] = Bitset(capacity=N)
            # map adjacency to bitsets
            for u, v in self.edge_iterator(labels=None):
                if u != v:
                    B[u].add(vertex_to_int[v])
                    B[v].add(vertex_to_int[u])
            # Search for a triangle
            for u, v in self.edge_iterator(labels=None):
                BB = B[u] & B[v]
                if BB:
                    if certificate:
                        for w in self.neighbor_iterator(u):
                            if vertex_to_int[w] in BB:
                                return False, [u, v, w]
                    return False
            return (True, []) if certificate else True

        elif algorithm == 'matrix':
            if self.order() < 3:
                return True
            return (self.adjacency_matrix()**3).trace() == 0

        raise ValueError("Algorithm '%s' not yet implemented. Please contribute." % (algorithm))

    @doc_index("Graph properties")
    def is_split(self):
        r"""
        Return ``True`` if the graph is a Split graph, ``False`` otherwise.

        A Graph `G` is said to be a split graph if its vertices `V(G)` can be
        partitioned into two sets `K` and `I` such that the vertices of `K`
        induce a complete graph, and those of `I` are an independent set.

        There is a simple test to check whether a graph is a split graph (see,
        for instance, the book "Graph Classes, a survey" [BLS1999]_ page
        203) :

        Given the degree sequence `d_1 \geq ... \geq d_n` of `G`, a graph is a
        split graph if and only if :

        .. MATH::

            \sum_{i=1}^\omega d_i = \omega (\omega - 1) + \sum_{i=\omega + 1}^nd_i

        where `\omega = max \{i:d_i\geq i-1\}`.

        EXAMPLES:

        Split graphs are, in particular, chordal graphs. Hence, The Petersen
        graph can not be split::

            sage: graphs.PetersenGraph().is_split()
            False

        We can easily build some "random" split graph by creating a complete
        graph, and adding vertices only connected to some random vertices of the
        clique::

            sage: g = graphs.CompleteGraph(10)
            sage: sets = Subsets(Set(range(10)))
            sage: for i in range(10, 25):
            ....:    g.add_edges([(i,k) for k in sets.random_element()])
            sage: g.is_split()
            True

        Another characterisation of split graph states that a graph is a split
        graph if and only if does not contain the 4-cycle, 5-cycle or `2K_2` as
        an induced subgraph. Hence for the above graph we have::

            sage: forbidden_subgraphs = [graphs.CycleGraph(4),
            ....:                        graphs.CycleGraph(5),
            ....:                        2 * graphs.CompleteGraph(2)]
            sage: sum(g.subgraph_search_count(H, induced=True)                          # needs sage.modules
            ....:     for H in forbidden_subgraphs)
            0
        """
        self._scream_if_not_simple()
        # our degree sequence is numbered from 0 to n-1, so to avoid
        # any mistake, let's fix it :-)
        degree_sequence = [0] + sorted(self.degree(), reverse=True)

        for i, d in enumerate(degree_sequence):
            if d >= i - 1:
                omega = i
            else:
                break

        left = sum(degree_sequence[:omega + 1])
        right = omega * (omega - 1) + sum(degree_sequence[omega + 1:])

        return left == right

    @doc_index("Algorithmically hard stuff")
    def is_perfect(self, certificate=False):
        r"""
