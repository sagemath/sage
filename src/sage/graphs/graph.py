# -*- coding: utf-8 -*-
r"""
Undirected graphs

This module implements functions and operations involving undirected
graphs.

{INDEX_OF_METHODS}

AUTHORS:

-  Robert L. Miller (2006-10-22): initial version

-  William Stein (2006-12-05): Editing

-  Robert L. Miller (2007-01-13): refactoring, adjusting for
   NetworkX-0.33, fixed plotting bugs (2007-01-23): basic tutorial,
   edge labels, loops, multiple edges and arcs (2007-02-07): graph6
   and sparse6 formats, matrix input

-  Emily Kirkmann (2007-02-11): added graph_border option to plot
   and show

-  Robert L. Miller (2007-02-12): vertex color-maps, graph
   boundaries, graph6 helper functions in Cython

-  Robert L. Miller Sage Days 3 (2007-02-17-21): 3d plotting in
   Tachyon

-  Robert L. Miller (2007-02-25): display a partition

-  Robert L. Miller (2007-02-28): associate arbitrary objects to
   vertices, edge and arc label display (in 2d), edge coloring

-  Robert L. Miller (2007-03-21): Automorphism group, isomorphism
   check, canonical label

-  Robert L. Miller (2007-06-07-09): NetworkX function wrapping

-  Michael W. Hansen (2007-06-09): Topological sort generation

-  Emily Kirkman, Robert L. Miller Sage Days 4: Finished wrapping
   NetworkX

-  Emily Kirkman (2007-07-21): Genus (including circular planar,
   all embeddings and all planar embeddings), all paths, interior
   paths

-  Bobby Moretti (2007-08-12): fixed up plotting of graphs with
   edge colors differentiated by label

-  Jason Grout (2007-09-25): Added functions, bug fixes, and
   general enhancements

-  Robert L. Miller (Sage Days 7): Edge labeled graph isomorphism

-  Tom Boothby (Sage Days 7): Miscellaneous awesomeness

-  Tom Boothby (2008-01-09): Added graphviz output

-  David Joyner (2009-2): Fixed docstring bug related to GAP.

-  Stephen Hartke (2009-07-26): Fixed bug in blocks_and_cut_vertices()
   that caused an incorrect result when the vertex 0 was a cut vertex.

-  Stephen Hartke (2009-08-22): Fixed bug in blocks_and_cut_vertices()
   where the list of cut_vertices is not treated as a set.

-  Anders Jonsson (2009-10-10): Counting of spanning trees and out-trees added.

-  Nathann Cohen (2009-09) : Cliquer, Connectivity, Flows
                             and everything that uses Linear Programming
                             and class numerical.MIP

- Nicolas M. Thiery (2010-02): graph layout code refactoring, dot2tex/graphviz interface

- David Coudert (2012-04) : Reduction rules in vertex_cover.

- Birk Eisermann (2012-06): added recognition of weakly chordal graphs and
                            long-hole-free / long-antihole-free graphs

- Alexandre P. Zuge (2013-07): added join operation.

- Amritanshu Prasad (2014-08): added clique polynomial

Graph Format
------------

Supported formats
~~~~~~~~~~~~~~~~~

Sage Graphs can be created from a wide range of inputs. A few
examples are covered here.


-  NetworkX dictionary format:

   ::

       sage: d = {0: [1,4,5], 1: [2,6], 2: [3,7], 3: [4,8], 4: [9], \
             5: [7, 8], 6: [8,9], 7: [9]}
       sage: G = Graph(d); G
       Graph on 10 vertices
       sage: G.plot().show()    # or G.show()

-  A NetworkX graph:

   ::

       sage: import networkx
       sage: K = networkx.complete_bipartite_graph(12,7)
       sage: G = Graph(K)
       sage: G.degree()
       [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 12, 12, 12, 12, 12, 12, 12]

-  graph6 or sparse6 format:

   ::

       sage: s = ':I`AKGsaOs`cI]Gb~'
       sage: G = Graph(s, sparse=True); G
       Looped multi-graph on 10 vertices
       sage: G.plot().show()    # or G.show()

   Note that the ``\`` character is an escape character in Python, and
   also a character used by graph6 strings:

   ::

       sage: G = Graph('Ihe\n@GUA')
       Traceback (most recent call last):
       ...
       RuntimeError: The string (Ihe) seems corrupt: for n = 10, the string is too short.

   In Python, the escaped character ``\`` is represented by ``\\``:

   ::

       sage: G = Graph('Ihe\\n@GUA')
       sage: G.plot().show()    # or G.show()

-  adjacency matrix: In an adjacency matrix, each column and each
   row represent a vertex. If a 1 shows up in row `i`, column
   `j`, there is an edge `(i,j)`.

   ::

       sage: M = Matrix([(0,1,0,0,1,1,0,0,0,0),(1,0,1,0,0,0,1,0,0,0), \
       (0,1,0,1,0,0,0,1,0,0), (0,0,1,0,1,0,0,0,1,0),(1,0,0,1,0,0,0,0,0,1), \
       (1,0,0,0,0,0,0,1,1,0), (0,1,0,0,0,0,0,0,1,1),(0,0,1,0,0,1,0,0,0,1), \
       (0,0,0,1,0,1,1,0,0,0), (0,0,0,0,1,0,1,1,0,0)])
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
       sage: G.plot().show()    # or G.show()

-  incidence matrix: In an incidence matrix, each row represents a
   vertex and each column represents an edge.

   ::

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
       sage: G.plot().show()    # or G.show()
       sage: DiGraph(matrix(2,[0,0,-1,1]), format="incidence_matrix")
       Traceback (most recent call last):
       ...
       ValueError: There must be two nonzero entries (-1 & 1) per column.

-  a list of edges::

       sage: g = Graph([(1,3),(3,8),(5,2)])
       sage: g
       Graph on 5 vertices

-  an igraph Graph::

       sage: import igraph                                # optional - python_igraph
       sage: g = Graph(igraph.Graph([(1,3),(3,2),(0,2)])) # optional - python_igraph
       sage: g                                            # optional - python_igraph
       Graph on 4 vertices

Generators
----------

Use ``graphs(n)`` to iterate through all non-isomorphic graphs of given size::

    sage: for g in graphs(4):
    ....:     print(g.spectrum())
    [0, 0, 0, 0]
    [1, 0, 0, -1]
    [1.4142135623..., 0, 0, -1.4142135623...]
    [2, 0, -1, -1]
    [1.7320508075..., 0, 0, -1.7320508075...]
    [1, 1, -1, -1]
    [1.6180339887..., 0.6180339887..., -0.6180339887..., -1.6180339887...]
    [2.1700864866..., 0.3111078174..., -1, -1.4811943040...]
    [2, 0, 0, -2]
    [2.5615528128..., 0, -1, -1.5615528128...]
    [3, -1, -1, -1]

Similarly ``graphs()`` will iterate through all graphs. The complete
graph of 4 vertices is of course the smallest graph with chromatic number
bigger than three::

    sage: for g in graphs():
    ....:     if g.chromatic_number() > 3:
    ....:         break
    sage: g.is_isomorphic(graphs.CompleteGraph(4))
    True

For some commonly used graphs to play with, type

::

    sage: graphs.[tab]          # not tested

and hit {tab}. Most of these graphs come with their own custom
plot, so you can see how people usually visualize these graphs.

::

    sage: G = graphs.PetersenGraph()
    sage: G.plot().show()    # or G.show()
    sage: G.degree_histogram()
    [0, 0, 0, 10]
    sage: G.adjacency_matrix()
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
    sage: S.plot().show()    # or S.show()
    sage: S.density()
    1/2

::

    sage: G = GraphQuery(display_cols=['graph6'], num_vertices=7, diameter=5)
    sage: L = G.get_graphs_list()
    sage: graphs_list.show_graphs(L)

.. _Graph:labels:

Labels
------

Each vertex can have any hashable object as a label. These are
things like strings, numbers, and tuples. Each edge is given a
default label of ``None``, but if specified, edges can
have any label at all. Edges between vertices `u` and
`v` are represented typically as ``(u, v, l)``, where
``l`` is the label for the edge.

Note that vertex labels themselves cannot be mutable items::

    sage: M = Matrix( [[0,0],[0,0]] )
    sage: G = Graph({ 0 : { M : None } })
    Traceback (most recent call last):
    ...
    TypeError: mutable matrices are unhashable

However, if one wants to define a dictionary, with the same keys
and arbitrary objects for entries, one can make that association::

    sage: d = {0 : graphs.DodecahedralGraph(), 1 : graphs.FlowerSnark(), \
          2 : graphs.MoebiusKantorGraph(), 3 : graphs.PetersenGraph() }
    sage: d[2]
    Moebius-Kantor Graph: Graph on 16 vertices
    sage: T = graphs.TetrahedralGraph()
    sage: T.vertices()
    [0, 1, 2, 3]
    sage: T.set_vertices(d)
    sage: T.get_vertex(1)
    Flower Snark: Graph on 20 vertices

Database
--------

There is a database available for searching for graphs that satisfy
a certain set of parameters, including number of vertices and
edges, density, maximum and minimum degree, diameter, radius, and
connectivity. To see a list of all search parameter keywords broken
down by their designated table names, type

::

    sage: graph_db_info()
    {...}

For more details on data types or keyword input, enter

::

    sage: GraphQuery?    # not tested

The results of a query can be viewed with the show method, or can be
viewed individually by iterating through the results:

::

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

Show each graph as you iterate through the results:

::

    sage: for g in Q:
    ....:     show(g)

Visualization
-------------

To see a graph `G` you are working with, there
are three main options. You can view the graph in two dimensions via
matplotlib with ``show()``. ::

    sage: G = graphs.RandomGNP(15,.3)
    sage: G.show()

And you can view it in three dimensions via jmol with ``show3d()``. ::

    sage: G.show3d()

Or it can be rendered with `\LaTeX`.  This requires the right
additions to a standard `\mbox{\rm\TeX}` installation.  Then standard
Sage commands, such as ``view(G)`` will display the graph, or
``latex(G)`` will produce a string suitable for inclusion in a
`\LaTeX` document.  More details on this are at
the :mod:`sage.graphs.graph_latex` module. ::

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
    TypeError: This graph is mutable, and thus not hashable. Create an immutable copy by `g.copy(immutable=True)`
    sage: G_immutable = Graph(G, immutable=True)
    sage: G_immutable == G
    True
    sage: {G_immutable:1}[G_immutable]
    1

Methods
-------
"""

#*****************************************************************************
#      Copyright (C) 2006 - 2007 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function
from __future__ import absolute_import

from copy import copy
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.misc.superseded import deprecation
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
import sage.graphs.generic_graph_pyx as generic_graph_pyx
from sage.graphs.generic_graph import GenericGraph
from sage.graphs.digraph import DiGraph
from sage.graphs.independent_sets import IndependentSets
from sage.combinat.combinatorial_map import combinatorial_map
from sage.misc.rest_index_of_methods import doc_index, gen_thematic_rest_table_index
from sage.misc.decorators import rename_keyword

class Graph(GenericGraph):
    r"""
    Undirected graph.

    A graph is a set of vertices connected by edges. See also the
    :wikipedia:`Wikipedia article on graphs <Graph_(mathematics)>`. For a
    collection of pre-defined graphs, see the
    :mod:`~sage.graphs.graph_generators` module.

    A :class:`Graph` object has many methods whose list can be obtained by
    typing ``g.<tab>`` (i.e. hit the 'tab' key) or by reading the documentation
    of :mod:`~sage.graphs.graph`, :mod:`~sage.graphs.generic_graph`, and
    :mod:`~sage.graphs.digraph`.

    INPUT:

    By default, a :class:`Graph` object is simple (i.e. no *loops* nor *multiple
    edges*) and unweighted. This can be easily tuned with the appropriate flags
    (see below).

    -  ``data`` -- can be any of the following (see the ``format`` argument):

      #. ``Graph()`` -- build a graph on 0 vertices.

      #. ``Graph(5)`` -- return an edgeless graph on the 5 vertices 0,...,4.

      #. ``Graph([list_of_vertices,list_of_edges])`` -- returns a graph with
         given vertices/edges.

         To bypass auto-detection, prefer the more explicit
         ``Graph([V,E],format='vertices_and_edges')``.

      #. ``Graph(list_of_edges)`` -- return a graph with a given list of edges
         (see documentation of
         :meth:`~sage.graphs.generic_graph.GenericGraph.add_edges`).

         To bypass auto-detection, prefer the more explicit ``Graph(L,
         format='list_of_edges')``.

      #. ``Graph({1:[2,3,4],3:[4]})`` -- return a graph by associating to each
         vertex the list of its neighbors.

         To bypass auto-detection, prefer the more explicit ``Graph(D,
         format='dict_of_lists')``.

      #. ``Graph({1: {2: 'a', 3:'b'} ,3:{2:'c'}})`` -- return a graph by
         associating a list of neighbors to each vertex and providing its edge
         label.

         To bypass auto-detection, prefer the more explicit ``Graph(D,
         format='dict_of_dicts')``.

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

         To bypass auto-detection, prefer the more explicit ``Graph(M,
         format='incidence_matrix')``.

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

    - ``pos`` - a positioning dictionary (cf. documentation of
      :meth:`~sage.graphs.generic_graph.GenericGraph.layout`). For example, to
      draw 4 vertices on a square::

         {0: [-1,-1],
          1: [ 1,-1],
          2: [ 1, 1],
          3: [-1, 1]}

    -  ``name`` - (must be an explicitly named parameter,
       i.e., ``name="complete")`` gives the graph a name

    -  ``loops`` - boolean, whether to allow loops (ignored
       if data is an instance of the ``Graph`` class)

    -  ``multiedges`` - boolean, whether to allow multiple
       edges (ignored if data is an instance of the ``Graph`` class).

    - ``weighted`` - whether graph thinks of itself as weighted or not. See
      :meth:`~sage.graphs.generic_graph.GenericGraph.weighted`.

    - ``format`` - if set to ``None`` (default), :class:`Graph` tries to guess
      input's format. To avoid this possibly time-consuming step, one of the
      following values can be specified (see description above): ``"int"``,
      ``"graph6"``, ``"sparse6"``, ``"rule"``, ``"list_of_edges"``,
      ``"dict_of_lists"``, ``"dict_of_dicts"``, ``"adjacency_matrix"``,
      ``"weighted_adjacency_matrix"``, ``"seidel_adjacency_matrix"``,
      ``"incidence_matrix"``, ``"NX"``, ``"igraph"``.

    - ``sparse`` (boolean) -- ``sparse=True`` is an alias for
      ``data_structure="sparse"``, and ``sparse=False`` is an alias for
      ``data_structure="dense"``.

    - ``data_structure`` -- one of the following (for more information, see
      :mod:`~sage.graphs.base.overview`)

       * ``"dense"`` -- selects the :mod:`~sage.graphs.base.dense_graph`
         backend.

       * ``"sparse"`` -- selects the :mod:`~sage.graphs.base.sparse_graph`
         backend.

       * ``"static_sparse"`` -- selects the
         :mod:`~sage.graphs.base.static_sparse_backend` (this backend is faster
         than the sparse backend and smaller in memory, and it is immutable, so
         that the resulting graphs can be used as dictionary keys).

    - ``immutable`` (boolean) -- whether to create a immutable graph. Note that
      ``immutable=True`` is actually a shortcut for
      ``data_structure='static_sparse'``. Set to ``False`` by default.

    - ``vertex_labels`` - Whether to allow any object as a vertex (slower), or
      only the integers `0,...,n-1`, where `n` is the number of vertices.

    -  ``convert_empty_dict_labels_to_None`` - this arguments sets
       the default edge labels used by NetworkX (empty dictionaries)
       to be replaced by None, the default Sage edge label. It is
       set to ``True`` iff a NetworkX graph is on the input.

    EXAMPLES:

    We illustrate the first seven input formats (the other two
    involve packages that are currently not standard in Sage):

    #. An integer giving the number of vertices::

        sage: g = Graph(5); g
        Graph on 5 vertices
        sage: g.vertices()
        [0, 1, 2, 3, 4]
        sage: g.edges()
        []

    #. A dictionary of dictionaries::

        sage: g = Graph({0:{1:'x',2:'z',3:'a'}, 2:{5:'out'}}); g
        Graph on 5 vertices

       The labels ('x', 'z', 'a', 'out') are labels for edges. For
       example, 'out' is the label for the edge on 2 and 5. Labels can be
       used as weights, if all the labels share some common parent.

       ::

        sage: a,b,c,d,e,f = sorted(SymmetricGroup(3))
        sage: Graph({b:{d:'c',e:'p'}, c:{d:'p',e:'c'}})
        Graph on 4 vertices

    #. A dictionary of lists::

        sage: g = Graph({0:[1,2,3], 2:[4]}); g
        Graph on 5 vertices

    #. A list of vertices and a function describing adjacencies. Note
       that the list of vertices and the function must be enclosed in a
       list (i.e., [list of vertices, function]).

       Construct the Paley graph over GF(13).

       ::

          sage: g=Graph([GF(13), lambda i,j: i!=j and (i-j).is_square()])
          sage: g.vertices()
          [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
          sage: g.adjacency_matrix()
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

       Construct the line graph of a complete graph.

       ::

          sage: g=graphs.CompleteGraph(4)
          sage: line_graph=Graph([g.edges(labels=false), \
                 lambda i,j: len(set(i).intersection(set(j)))>0], \
                 loops=False)
          sage: line_graph.vertices()
          [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
          sage: line_graph.adjacency_matrix()
          [0 1 1 1 1 0]
          [1 0 1 1 0 1]
          [1 1 0 0 1 1]
          [1 1 0 0 1 1]
          [1 0 1 1 0 1]
          [0 1 1 1 1 0]

    #. A graph6 or sparse6 string: Sage automatically recognizes
       whether a string is in graph6 or sparse6 format::

           sage: s = ':I`AKGsaOs`cI]Gb~'
           sage: Graph(s,sparse=True)
           Looped multi-graph on 10 vertices

       ::

           sage: G = Graph('G?????')
           sage: G = Graph("G'?G?C")
           Traceback (most recent call last):
           ...
           RuntimeError: The string seems corrupt: valid characters are
           ?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
           sage: G = Graph('G??????')
           Traceback (most recent call last):
           ...
           RuntimeError: The string (G??????) seems corrupt: for n = 8, the string is too long.

       ::

          sage: G = Graph(":I'AKGsaOs`cI]Gb~")
          Traceback (most recent call last):
          ...
          RuntimeError: The string seems corrupt: valid characters are
          ?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~

       There are also list functions to take care of lists of graphs::

           sage: s = ':IgMoqoCUOqeb\n:I`AKGsaOs`cI]Gb~\n:I`EDOAEQ?PccSsge\N\n'
           sage: graphs_list.from_sparse6(s)
           [Looped multi-graph on 10 vertices, Looped multi-graph on 10 vertices, Looped multi-graph on 10 vertices]

    #. A Sage matrix:
       Note: If format is not specified, then Sage assumes a symmetric square
       matrix is an adjacency matrix, otherwise an incidence matrix.

       - an adjacency matrix::

            sage: M = graphs.PetersenGraph().am(); M
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
            sage: Graph(M)
            Graph on 10 vertices

         ::

            sage: Graph(matrix([[1,2],[2,4]]),loops=True,sparse=True)
            Looped multi-graph on 2 vertices

            sage: M = Matrix([[0,1,-1],[1,0,-1/2],[-1,-1/2,0]]); M
            [   0    1   -1]
            [   1    0 -1/2]
            [  -1 -1/2    0]
            sage: G = Graph(M,sparse=True); G
            Graph on 3 vertices
            sage: G.weighted()
            True

       - an incidence matrix::

            sage: M = Matrix(6, [-1,0,0,0,1, 1,-1,0,0,0, 0,1,-1,0,0, 0,0,1,-1,0, 0,0,0,1,-1, 0,0,0,0,0]); M
            [-1  0  0  0  1]
            [ 1 -1  0  0  0]
            [ 0  1 -1  0  0]
            [ 0  0  1 -1  0]
            [ 0  0  0  1 -1]
            [ 0  0  0  0  0]
            sage: Graph(M)
            Graph on 6 vertices

            sage: Graph(Matrix([[1],[1],[1]]))
            Traceback (most recent call last):
            ...
            ValueError: There must be one or two nonzero entries per column in an incidence matrix. Got entries [1, 1, 1] in column 0
            sage: Graph(Matrix([[1],[1],[0]]))
            Graph on 3 vertices

            sage: M = Matrix([[0,1,-1],[1,0,-1],[-1,-1,0]]); M
            [ 0  1 -1]
            [ 1  0 -1]
            [-1 -1  0]
            sage: Graph(M,sparse=True)
            Graph on 3 vertices

            sage: M = Matrix([[0,1,1],[1,0,1],[-1,-1,0]]); M
            [ 0  1  1]
            [ 1  0  1]
            [-1 -1  0]
            sage: Graph(M)
            Traceback (most recent call last):
            ...
            ValueError: There must be one or two nonzero entries per column in an incidence matrix. Got entries [1, 1] in column 2

        Check that :trac:`9714` is fixed::

            sage: MA = Matrix([[1,2,0], [0,2,0], [0,0,1]])
            sage: GA = Graph(MA, format='adjacency_matrix')
            sage: MI = GA.incidence_matrix(oriented=False)
            sage: MI
            [2 1 1 0 0 0]
            [0 1 1 2 2 0]
            [0 0 0 0 0 2]
            sage: Graph(MI).edges(labels=None)
            [(0, 0), (0, 1), (0, 1), (1, 1), (1, 1), (2, 2)]

            sage: M = Matrix([[1], [-1]]); M
            [ 1]
            [-1]
            sage: Graph(M).edges()
            [(0, 1, None)]

    #.  A Seidel adjacency matrix::

          sage: from sage.combinat.matrices.hadamard_matrix import \
          ....:  regular_symmetric_hadamard_matrix_with_constant_diagonal as rshcd
          sage: m=rshcd(16,1)- matrix.identity(16)
          sage: Graph(m,format="seidel_adjacency_matrix").is_strongly_regular(parameters=True)
          (16, 6, 2, 2)

    #. a list of edges, or labelled edges::

          sage: g = Graph([(1,3),(3,8),(5,2)])
          sage: g
          Graph on 5 vertices

          sage: g = Graph([(1,2,"Peace"),(7,-9,"and"),(77,2, "Love")])
          sage: g
          Graph on 5 vertices
          sage: g = Graph([(0, 2, '0'), (0, 2, '1'), (3, 3, '2')], loops=True, multiedges=True)
          sage: g.loops()
          [(3, 3, '2')]

    #. A NetworkX MultiGraph::

          sage: import networkx
          sage: g = networkx.MultiGraph({0:[1,2,3], 2:[4]})
          sage: Graph(g)
          Graph on 5 vertices

    #. A NetworkX graph::

           sage: import networkx
           sage: g = networkx.Graph({0:[1,2,3], 2:[4]})
           sage: DiGraph(g)
           Digraph on 5 vertices

    #. An igraph Graph (see also
       :meth:`~sage.graphs.generic_graph.GenericGraph.igraph_graph`)::

           sage: import igraph                   # optional - python_igraph
           sage: g = igraph.Graph([(0,1),(0,2)]) # optional - python_igraph
           sage: Graph(g)                        # optional - python_igraph
           Graph on 3 vertices

       If ``vertex_labels`` is ``True``, the names of the vertices are given by
       the vertex attribute ``'name'``, if available::

           sage: g = igraph.Graph([(0,1),(0,2)], vertex_attrs={'name':['a','b','c']})  # optional - python_igraph
           sage: Graph(g).vertices()                                                   # optional - python_igraph
           ['a', 'b', 'c']
           sage: g = igraph.Graph([(0,1),(0,2)], vertex_attrs={'label':['a','b','c']}) # optional - python_igraph
           sage: Graph(g).vertices()                                                   # optional - python_igraph
           [0, 1, 2]

       If the igraph Graph has edge attributes, they are used as edge labels::

           sage: g = igraph.Graph([(0,1),(0,2)], edge_attrs={'name':['a','b'], 'weight':[1,3]}) # optional - python_igraph
           sage: Graph(g).edges()                                                               # optional - python_igraph
           [(0, 1, {'name': 'a', 'weight': 1}), (0, 2, {'name': 'b', 'weight': 3})]


    When defining an undirected graph from a function ``f``, it is *very*
    important that ``f`` be symmetric. If it is not, anything can happen::

        sage: f_sym = lambda x,y : abs(x-y) == 1
        sage: f_nonsym = lambda x,y : (x-y) == 1
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
          TypeError: This graph is mutable, and thus not hashable. Create an immutable copy by `g.copy(immutable=True)`

    When providing the optional arguments ``data_structure="static_sparse"``
    or ``immutable=True`` (both mean the same), then an immutable graph
    results. ::

          sage: G_imm = Graph(G, immutable=True)
          sage: H_imm = Graph(G, data_structure='static_sparse')
          sage: G_imm == H_imm == G
          True
          sage: {G_imm:1}[H_imm]
          1

    TESTS::

        sage: Graph(4,format="HeyHeyHey")
        Traceback (most recent call last):
        ...
        ValueError: Unknown input format 'HeyHeyHey'

        sage: Graph(igraph.Graph(directed=True)) # optional - python_igraph
        Traceback (most recent call last):
        ...
        ValueError: An *undirected* igraph graph was expected. To build an directed graph, call the DiGraph constructor.

        sage: m = matrix([[0,-1],[-1,0]])
        sage: Graph(m,format="seidel_adjacency_matrix")
        Graph on 2 vertices
        sage: m[0,1]=1
        sage: Graph(m,format="seidel_adjacency_matrix")
        Traceback (most recent call last):
        ...
        ValueError: Graph's Seidel adjacency matrix must be symmetric

        sage: m[0,1]=-1; m[1,1]=1
        sage: Graph(m,format="seidel_adjacency_matrix")
        Traceback (most recent call last):
        ...
        ValueError: Graph's Seidel adjacency matrix must have 0s on the main diagonal

    From a a list of vertices and a list of edges::

        sage: G = Graph([[1,2,3],[(1,2)]]); G
        Graph on 3 vertices
        sage: G.edges()
        [(1, 2, None)]
    """
    _directed = False

    def __init__(self, data=None, pos=None, loops=None, format=None,
                 weighted=None, implementation='c_graph',
                 data_structure="sparse", vertex_labels=True, name=None,
                 multiedges=None, convert_empty_dict_labels_to_None=None,
                 sparse=True, immutable=False):
        """
        TESTS::

            sage: G = Graph()
            sage: loads(dumps(G)) == G
            True
            sage: a = matrix(2,2,[1,0,0,1])
            sage: Graph(a).adjacency_matrix() == a
            True

            sage: a = matrix(2,2,[2,0,0,1])
            sage: Graph(a,sparse=True).adjacency_matrix() == a
            True

        The positions are copied when the graph is built from
        another graph ::

            sage: g = graphs.PetersenGraph()
            sage: h = Graph(g)
            sage: g.get_pos() == h.get_pos()
            True

        Or from a DiGraph ::

            sage: d = DiGraph(g)
            sage: h = Graph(d)
            sage: g.get_pos() == h.get_pos()
            True

        Loops are not counted as multiedges (see :trac:`11693`) and edges are
        not counted twice ::

            sage: Graph({1:[1]}).num_edges()
            1
            sage: Graph({1:[2,2]}).num_edges()
            2

        An empty list or dictionary defines a simple graph
        (:trac:`10441` and :trac:`12910`)::

            sage: Graph([])
            Graph on 0 vertices
            sage: Graph({})
            Graph on 0 vertices
            sage: # not "Multi-graph on 0 vertices"

        Verify that the int format works as expected (:trac:`12557`)::

            sage: Graph(2).adjacency_matrix()
            [0 0]
            [0 0]
            sage: Graph(3) == Graph(3,format='int')
            True

        Problem with weighted adjacency matrix (:trac:`13919`)::

            sage: B = {0:{1:2,2:5,3:4},1:{2:2,4:7},2:{3:1,4:4,5:3},3:{5:4},4:{5:1,6:5},5:{6:7}}
            sage: grafo3 = Graph(B,weighted=True)
            sage: matad = grafo3.weighted_adjacency_matrix()
            sage: grafo4 = Graph(matad,format = "adjacency_matrix", weighted=True)
            sage: grafo4.shortest_path(0,6,by_weight=True)
            [0, 1, 2, 5, 4, 6]

        Graphs returned when setting ``immutable=False`` are mutable::

            sage: g = graphs.PetersenGraph()
            sage: g = Graph(g.edges(),immutable=False)
            sage: g.add_edge("Hey", "Heyyyyyyy")

        And their name is set::

            sage: g = graphs.PetersenGraph()
            sage: Graph(g, immutable=True)
            Petersen graph: Graph on 10 vertices

        Check error messages for graphs built from incidence matrices (see
        :trac:`18440`)::

            sage: Graph(matrix([[-1, 1, 0],[1, 0, 0]]))
            Traceback (most recent call last):
            ...
            ValueError: Column 1 of the (oriented) incidence matrix contains
            only one nonzero value
            sage: Graph(matrix([[1,1],[1,1],[1,0]]))
            Traceback (most recent call last):
            ...
            ValueError: There must be one or two nonzero entries per column in an incidence matrix. Got entries [1, 1, 1] in column 0
            sage: Graph(matrix([[3,1,1],[0,1,1]]))
            Traceback (most recent call last):
            ...
            ValueError: Each column of a non-oriented incidence matrix must sum
            to 2, but column 0 does not
        """
        GenericGraph.__init__(self)

        from sage.structure.element import is_Matrix

        if sparse is False:
            if data_structure != "sparse":
                raise ValueError("The 'sparse' argument is an alias for "
                                 "'data_structure'. Please do not define both.")
            data_structure = "dense"

        # Choice of the backend

        if implementation != 'c_graph':
            deprecation(18375,"The 'implementation' keyword is deprecated, "
                        "and the graphs has been stored as a 'c_graph'")

        if multiedges or weighted:
            if data_structure == "dense":
                raise RuntimeError("Multiedge and weighted c_graphs must be sparse.")
        if immutable:
            data_structure = 'static_sparse'

        # If the data structure is static_sparse, we first build a graph
        # using the sparse data structure, then reencode the resulting graph
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
        if format is None and is_Matrix(data):
            if data.is_symmetric():
                format = 'adjacency_matrix'
            else:
                format = 'incidence_matrix'
        if format is None and isinstance(data, Graph):
            format = 'Graph'
        from sage.graphs.all import DiGraph
        if format is None and isinstance(data, DiGraph):
            data = data.to_undirected()
            format = 'Graph'
        if (format is None        and
            isinstance(data,list) and
            len(data)>=2          and
            callable(data[1])):
            format = 'rule'

        if (format is None           and
            isinstance(data,list)    and
            len(data) == 2           and
            isinstance(data[0],list) and # a list of two lists, the second of
            isinstance(data[1],list) and # which contains iterables (the edges)
            (not data[1] or callable(getattr(data[1][0],"__iter__",None)))):
            format = "vertices_and_edges"

        if format is None and isinstance(data,dict):
            keys = data.keys()
            if len(keys) == 0: format = 'dict_of_dicts'
            else:
                if isinstance(data[keys[0]], list):
                    format = 'dict_of_lists'
                elif isinstance(data[keys[0]], dict):
                    format = 'dict_of_dicts'
        if format is None and hasattr(data, 'adj'):
            import networkx
            if isinstance(data, (networkx.DiGraph, networkx.MultiDiGraph)):
                data = data.to_undirected()
            elif isinstance(data, (networkx.Graph, networkx.MultiGraph)):
                format = 'NX'

        if (format is None          and
            hasattr(data, 'vcount') and
            hasattr(data, 'get_edgelist')):
            try:
                import igraph
            except ImportError:
                raise ImportError("The data seems to be a igraph object, but "+
                                  "igraph is not installed in Sage. To install "+
                                  "it, run 'sage -i python_igraph'")
            if format is None and isinstance(data, igraph.Graph):
                format = 'igraph'
        if format is None and isinstance(data, (int, Integer)):
            format = 'int'
        if format is None and data is None:
            format = 'int'
            data = 0

        # Input is a list of edges
        if format is None and isinstance(data,list):
            format = "list_of_edges"
            if weighted is None: weighted = False
            num_verts=0

        if format is None:
            raise ValueError("This input cannot be turned into a graph")

        if format == 'weighted_adjacency_matrix':
            if weighted is False:
                raise ValueError("Format was weighted_adjacency_matrix but weighted was False.")
            if weighted   is None: weighted   = True
            if multiedges is None: multiedges = False
            format = 'adjacency_matrix'

        # At this point, 'format' has been set. We build the graph

        if format == 'graph6':
            if weighted   is None: weighted   = False
            self.allow_loops(loops if loops else False, check=False)
            self.allow_multiple_edges(multiedges if multiedges else False, check=False)
            from .graph_input import from_graph6
            from_graph6(self, data)

        elif format == 'sparse6':
            if weighted   is None: weighted   = False
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
            multiedges = False
            weighted = False
            loops = False
            self.allow_loops(False)
            self.allow_multiple_edges(False)
            from .graph_input import from_seidel_adjacency_matrix
            from_seidel_adjacency_matrix(self, data)
        elif format == 'Graph':
            if loops is None:      loops      = data.allows_loops()
            if multiedges is None: multiedges = data.allows_multiple_edges()
            if weighted is None:   weighted   = data.weighted()
            self.allow_loops(loops, check=False)
            self.allow_multiple_edges(multiedges, check=False)
            if data.get_pos() is not None:
                pos = data.get_pos().copy()
            self.name(data.name())
            self.add_vertices(data.vertex_iterator())
            self.add_edges(data.edge_iterator())
        elif format == 'NX':
            if convert_empty_dict_labels_to_None is not False:
                r = lambda x:None if x=={} else x
            else:
                r = lambda x:x
            if weighted is None:
                if isinstance(data, networkx.Graph):
                    weighted = False
                    if multiedges is None:
                        multiedges = False
                    if loops is None:
                        loops = False
                else:
                    weighted = True
                    if multiedges is None:
                        multiedges = True
                    if loops is None:
                        loops = True
            self.allow_loops(loops, check=False)
            self.allow_multiple_edges(multiedges, check=False)
            self.add_vertices(data.nodes())
            self.add_edges((u,v,r(l)) for u,v,l in data.edges_iter(data=True))
        elif format == 'igraph':
            if data.is_directed():
                raise ValueError("An *undirected* igraph graph was expected. "+
                                 "To build an directed graph, call the DiGraph "+
                                 "constructor.")

            self.add_vertices(range(data.vcount()))
            self.add_edges([(e.source, e.target, e.attributes()) for e in data.es()])

            if vertex_labels and 'name' in data.vertex_attributes():
                vs = data.vs()
                self.relabel({v:vs[v]['name'] for v in self})

        elif format == 'rule':
            f = data[1]
            verts = data[0]
            if loops is None: loops = any(f(v,v) for v in verts)
            if weighted is None: weighted = False
            self.allow_loops(loops, check=False)
            self.allow_multiple_edges(True if multiedges else False, check=False)
            from itertools import combinations
            self.add_vertices(verts)
            self.add_edges(e for e in combinations(verts,2) if f(*e))
            self.add_edges((v,v) for v in verts if f(v,v))

        elif format == "vertices_and_edges":
            self.allow_multiple_edges(bool(multiedges), check=False)
            self.allow_loops(bool(loops), check=False)
            self.add_vertices(data[0])
            self.add_edges(data[1])

        elif format == 'dict_of_dicts':
            from .graph_input import from_dict_of_dicts
            from_dict_of_dicts(self, data, loops=loops, multiedges=multiedges, weighted=weighted,
                               convert_empty_dict_labels_to_None = False if convert_empty_dict_labels_to_None is None else convert_empty_dict_labels_to_None)

        elif format == 'dict_of_lists':
            from .graph_input import from_dict_of_lists
            from_dict_of_lists(self, data, loops=loops, multiedges=multiedges, weighted=weighted)

        elif format == 'int':
            self.allow_loops(loops if loops else False, check=False)
            self.allow_multiple_edges(multiedges if multiedges else False, check=False)
            if data<0:
                raise ValueError("The number of vertices cannot be strictly negative!")
            if data:
                self.add_vertices(range(data))

        elif format == 'list_of_edges':
            self.allow_multiple_edges(False if multiedges is False else True, check=False)
            self.allow_loops(False if loops is False else True, check=False)
            self.add_edges(data)
            if multiedges is not True and self.has_multiple_edges():
                deprecation(15706, "You created a graph with multiple edges "
                            "from a list. Please set 'multiedges' to 'True' "
                            "when you do so, as in the future the default "
                            "behaviour will be to ignore those edges")
            elif multiedges is None:
                self.allow_multiple_edges(False, check=False)

            if loops is not True and self.has_loops():
                deprecation(15706, "You created a graph with loops from a list. "+
                            "Please set 'loops' to 'True' when you do so, as in "+
                            "the future the default behaviour will be to ignore "+
                            "those edges")
            elif loops is None:
                self.allow_loops(False, check=False)
        else:
            raise ValueError("Unknown input format '{}'".format(format))

        if weighted   is None: weighted   = False
        self._weighted = getattr(self,'_weighted',weighted)

        self._pos = pos

        if format != 'Graph' or name is not None:
            self.name(name)

        if data_structure == "static_sparse":
            from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            ib = StaticSparseBackend(self,
                                     loops = self.allows_loops(),
                                     multiedges = self.allows_multiple_edges())
            self._backend = ib
            self._immutable = True

    ### Formats

    @doc_index("Basic methods")
    def graph6_string(self):
        """
        Returns the graph6 representation of the graph as an ASCII string.
        Only valid for simple (no loops, multiple edges) graphs on 0 to
        262143 vertices.

        .. NOTE::

            As the graph6 format only handles graphs whose vertex set is
            `\{0,...,n-1\}`, a :meth:`relabelled copy
            <sage.graphs.generic_graph.GenericGraph.relabel>` of your graph will
            be encoded if necessary.

        EXAMPLES::

            sage: G = graphs.KrackhardtKiteGraph()
            sage: G.graph6_string()
            'IvUqwK@?G'
        """
        n = self.order()
        if n > 262143:
            raise ValueError('graph6 format supports graphs on 0 to 262143 vertices only.')
        elif self.has_loops() or self.has_multiple_edges():
            raise ValueError('graph6 format supports only simple graphs (no loops, no multiple edges)')
        else:
            return generic_graph_pyx.small_integer_to_graph6(n) + generic_graph_pyx.binary_string_to_graph6(self._bit_vector())

    @doc_index("Basic methods")
    def sparse6_string(self):
        r"""
        Returns the sparse6 representation of the graph as an ASCII string.
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

            sage: G = Graph()
            sage: G.sparse6_string()
            ':?'

        ::

            sage: G = Graph(loops=True, multiedges=True,data_structure="sparse")
            sage: Graph(':?',data_structure="sparse") == G
            True

        TEST:

        Check that :trac:`18445` is fixed::

            sage: Graph(graphs.KneserGraph(5,2).sparse6_string()).size()
            15

        """
        n = self.order()
        if n == 0:
            return ':?'
        if n > 262143:
            raise ValueError('sparse6 format supports graphs on 0 to 262143 vertices only.')
        else:
            v_to_int = {v:i for i,v in enumerate(self.vertices())}
            edges = [sorted((v_to_int[u],v_to_int[v])) for u,v in self.edge_iterator(labels=False)]
            edges.sort(key=lambda e: (e[1],e[0])) # reverse lexicographic order

            # encode bit vector
            from math import ceil
            from sage.misc.functional import log
            k = int(ceil(log(n,2)))
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
            s = s + ( '1' * ((6 - len(s))%6) )

            # split into groups of 6, and convert numbers to decimal, adding 63
            six_bits = ''
            for i in range(len(s)//6):
                six_bits += chr( int( s[6*i:6*(i+1)], 2) + 63 )
            return ':' + generic_graph_pyx.small_integer_to_graph6(n) + six_bits

    ### Attributes

    @combinatorial_map(name="partition of connected components")
    @doc_index("Deprecated")
    def to_partition(self):
        """
        Return the partition of connected components of ``self``.

        EXAMPLES::

            sage: for x in graphs(3):    print(x.to_partition())
            doctest:...: DeprecationWarning: Please use G.connected_components_sizes() instead
            See http://trac.sagemath.org/17449 for details.
            [1, 1, 1]
            [2, 1]
            [3]
            [3]
        """
        from sage.misc.superseded import deprecation
        deprecation(17449, "Please use G.connected_components_sizes() instead")

        from sage.combinat.partition import Partition
        return Partition(sorted([len(y) for y in self.connected_components()], reverse=True))

    @doc_index("Basic methods")
    def is_directed(self):
        """
        Since graph is undirected, returns False.

        EXAMPLES::

            sage: Graph().is_directed()
            False
        """
        return False

    @doc_index("Connectivity, orientations, trees")
    def bridges(self):
        r"""
        Returns a list of the bridges (or cut edges).

        A bridge is an edge so that deleting it disconnects the graph.

        .. NOTE::

            This method assumes the graph is connected.

        EXAMPLES::

             sage: g = 2*graphs.PetersenGraph()
             sage: g.add_edge(1,10)
             sage: g.is_connected()
             True
             sage: g.bridges()
             [(1, 10, None)]
        """
        gs = self.strong_orientation()
        bridges = []
        for scc in gs.strongly_connected_components():
            bridges.extend(gs.edge_boundary(scc))
        return bridges

    @doc_index("Connectivity, orientations, trees")
    def spanning_trees(self):
        """
        Returns a list of all spanning trees.

        If the graph is disconnected, returns the empty list.

        Uses the Read-Tarjan backtracking algorithm [RT75]_.

        EXAMPLES::

            sage: G = Graph([(1,2),(1,2),(1,3),(1,3),(2,3),(1,4)],multiedges=True)
            sage: len(G.spanning_trees())
            8
            sage: G.spanning_trees_count()
            8
            sage: G = Graph([(1,2),(2,3),(3,1),(3,4),(4,5),(4,5),(4,6)],multiedges=True)
            sage: len(G.spanning_trees())
            6
            sage: G.spanning_trees_count()
            6

        .. SEEALSO::

            - :meth:`~sage.graphs.generic_graph.GenericGraph.spanning_trees_count`
              -- counts the number of spanning trees.

            - :meth:`~sage.graphs.graph.Graph.random_spanning_tree`
              -- returns a random spanning tree.

        TESTS:

        Works with looped graphs::

            sage: g = Graph({i:[i,(i+1)%6] for i in range(6)})
            sage: g.spanning_trees()
            [Graph on 6 vertices,
             Graph on 6 vertices,
             Graph on 6 vertices,
             Graph on 6 vertices,
             Graph on 6 vertices,
             Graph on 6 vertices]

        REFERENCES:

        .. [RT75] Read, R. C. and Tarjan, R. E.
          Bounds on Backtrack Algorithms for Listing Cycles, Paths, and Spanning Trees
          Networks, Volume 5 (1975), numer 3, pages 237-252.
        """

        def _recursive_spanning_trees(G,forest):
            """
            Returns all the spanning trees of G containing forest
            """
            if not G.is_connected():
                return []

            if G.size() == forest.size():
                return [forest.copy()]
            else:
                # Pick an edge e from G-forest
                for e in G.edge_iterator(labels=False):
                    if not forest.has_edge(e):
                        break

                # 1) Recursive call with e removed from G
                G.delete_edge(e)
                trees = _recursive_spanning_trees(G,forest)
                G.add_edge(e)

                # 2) Recursive call with e include in forest
                #
                # e=xy links the CC (connected component) of forest containing x
                # with the CC containing y. Any other edge which does that
                # cannot be added to forest anymore, and B is the list of them
                c1 = forest.connected_component_containing_vertex(e[0])
                c2 = forest.connected_component_containing_vertex(e[1])
                G.delete_edge(e)
                B = G.edge_boundary(c1,c2,sort=False)
                G.add_edge(e)

                # Actual call
                forest.add_edge(e)
                G.delete_edges(B)
                trees.extend(_recursive_spanning_trees(G,forest))
                G.add_edges(B)
                forest.delete_edge(e)

                return trees

        if self.is_connected() and len(self):
            forest = Graph([])
            forest.add_vertices(self.vertices())
            forest.add_edges(self.bridges())
            return _recursive_spanning_trees(Graph(self,immutable=False,loops=False), forest)
        else:
            return []

    ### Properties
    @doc_index("Graph properties")
    def is_tree(self, certificate=False, output='vertex'):
        """
        Tests if the graph is a tree

        INPUT:

        - ``certificate`` (boolean) -- whether to return a certificate. The
          method only returns boolean answers when ``certificate = False``
          (default). When it is set to ``True``, it either answers ``(True,
          None)`` when the graph is a tree and ``(False, cycle)`` when it
          contains a cycle. It returns ``(False, None)`` when the graph is not
          connected.

        - ``output`` (``'vertex'`` (default) or ``'edge'``) -- whether the
          certificate is given as a list of vertices or a list of
          edges.

        When the certificate cycle is given as a list of edges, the
        edges are given as `(v_i, v_{i+1}, l)` where `v_1, v_2, \dots,
        v_n` are the vertices of the cycles (in their cyclic order).

        EXAMPLES::

            sage: all(T.is_tree() for T in graphs.trees(15))
            True

        The empty graph is not considered to be a tree::

            sage: graphs.EmptyGraph().is_tree()
            False

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
            (False, [(1, 2, 'a'), (2, 1, 'b')])

        TESTS:

        :trac:`14434` is fixed::

            sage: g = Graph({0:[1,4,5],3:[4,8,9],4:[9],5:[7,8],7:[9]})
            sage: _,cycle = g.is_tree(certificate=True)
            sage: g.size()
            10
            sage: g.add_cycle(cycle)
            sage: g.size()
            10
        """
        if not output in ['vertex', 'edge']:
            raise ValueError('output must be either vertex or edge')

        if self.order() == 0:
            return False

        if not self.is_connected():
            return (False, None) if certificate else False

        if certificate:
            if self.num_verts() == self.num_edges() + 1:
                return (True, None)

            if self.has_multiple_edges():
                if output == 'vertex':
                    return (False, list(self.multiple_edges()[0][:2]))
                edge1, edge2 = self.multiple_edges()[:2]
                if edge1[0] != edge2[0]:
                    return (False, [edge1, edge2])
                return (False, [edge1, (edge2[1], edge2[0], edge2[2])])

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
            n = self.order()
            seen = {}
            u = next(self.vertex_iterator())
            seen[u] = u
            stack = [(u, v) for v in self.neighbor_iterator(u)]
            while stack:
                u, v = stack.pop(-1)
                if v in seen:
                    continue
                for w in self.neighbors(v):
                    if u == w:
                        continue
                    elif w in seen:
                        cycle = [v, w]
                        while u != w:
                            cycle.insert(0, u)
                            u = seen[u]
                        if output == 'vertex':
                            return (False, cycle)
                        return (False, vertices_to_edges(cycle))
                    else:
                        stack.append((v, w))
                seen[v] = u

        else:
            return self.num_verts() == self.num_edges() + 1

    @doc_index("Graph properties")
    def is_forest(self, certificate=False, output='vertex'):
        """
        Tests if the graph is a forest, i.e. a disjoint union of trees.

        INPUT:

        - ``certificate`` (boolean) -- whether to return a certificate. The
          method only returns boolean answers when ``certificate = False``
          (default). When it is set to ``True``, it either answers ``(True,
          None)`` when the graph is a forest and ``(False, cycle)`` when it
          contains a cycle.

        - ``output`` (``'vertex'`` (default) or ``'edge'``) -- whether the
          certificate is given as a list of vertices or a list of
          edges.

        EXAMPLES::

            sage: seven_acre_wood = sum(graphs.trees(7), Graph())
            sage: seven_acre_wood.is_forest()
            True

        With certificates::

            sage: g = graphs.RandomTree(30)
            sage: g.is_forest(certificate=True)
            (True, None)
            sage: (2*g + graphs.PetersenGraph() + g).is_forest(certificate=True)
            (False, [63, 62, 61, 60, 64])
        """
        number_of_connected_components = len(self.connected_components())
        isit = (self.num_verts() ==
                self.num_edges() + number_of_connected_components)

        if not certificate:
            return isit
        else:
            if isit:
                return (True, None)
            # The graph contains a cycle, and the user wants to see it.

            # No need to copy the graph
            if number_of_connected_components == 1:
                return self.is_tree(certificate=True, output=output)

            # We try to find a cycle in each connected component
            for gg in self.connected_components_subgraphs():
                isit, cycle = gg.is_tree(certificate=True, output=output)
                if not isit:
                    return (False, cycle)

    @doc_index("Graph properties")
    def is_overfull(self):
        r"""
        Tests whether the current graph is overfull.

        A graph `G` on `n` vertices and `m` edges is said to
        be overfull if:

        - `n` is odd

        - It satisfies `2m > (n-1)\Delta(G)`, where
          `\Delta(G)` denotes the maximum degree
          among all vertices in `G`.

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
            sage: edge_coloring(g, value_only=True)
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
            ...       i = 0
            ...       while i <= n:
            ...           if graphs.CompleteGraph(i).is_overfull():
            ...               print("A complete graph of even order cannot be overfull.")
            ...               return
            ...           i += 2
            ...       print("Complete graphs of even order up to %s are not overfull." % n)
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
            ...       i = 3
            ...       while i <= n:
            ...           if not graphs.CompleteGraph(i).is_overfull():
            ...               print("A complete graph of odd order > 1 must be overfull.")
            ...               return
            ...           i += 2
            ...       print("Complete graphs of odd order > 1 up to %s are overfull." % n)
            ...
            sage: check_overfull_Kn_odd(100)  # long time
            Complete graphs of odd order > 1 up to 100 are overfull.

        The Petersen Graph, though, is not overfull while
        its chromatic index is `\Delta+1`::

            sage: g = graphs.PetersenGraph()
            sage: g.is_overfull()
            False
            sage: from sage.graphs.graph_coloring import edge_coloring
            sage: max(g.degree()) + 1 ==  edge_coloring(g, value_only=True)
            True
        """
        # # A possible optimized version. But the gain in speed is very little.
        # return bool(self._backend.num_verts() & 1) and (  # odd order n
        #     2 * self._backend.num_edges(self._directed) > #2m > \Delta(G)*(n-1)
        #     max(self.degree()) * (self._backend.num_verts() - 1))
        # unoptimized version
        return (self.order() % 2 == 1) and (
            2 * self.size() > max(self.degree()) * (self.order() - 1))

    @doc_index("Graph properties")
    def is_even_hole_free(self, certificate = False):
        r"""
        Tests whether ``self`` contains an induced even hole.

        A Hole is a cycle of length at least 4 (included). It is said
        to be even (resp. odd) if its length is even (resp. odd).

        Even-hole-free graphs always contain a bisimplicial vertex,
        which ensures that their chromatic number is at most twice
        their clique number [ABCHRS08]_.

        INPUT:

        - ``certificate`` (boolean) -- When ``certificate = False``,
          this method only returns ``True`` or ``False``. If
          ``certificate = True``, the subgraph found is returned
          instead of ``False``.

        EXAMPLE:

        Is the Petersen Graph even-hole-free ::

            sage: g = graphs.PetersenGraph()
            sage: g.is_even_hole_free()
            False

        As any chordal graph is hole-free, interval graphs behave the
        same way::

            sage: g = graphs.RandomIntervalGraph(20)
            sage: g.is_even_hole_free()
            True

        It is clear, though, that a random Bipartite Graph which is
        not a forest has an even hole::

            sage: g = graphs.RandomBipartite(10, 10, .5)
            sage: g.is_even_hole_free() and not g.is_forest()
            False

        We can check the certificate returned is indeed an even
        cycle::

            sage: if not g.is_forest():
            ...      cycle = g.is_even_hole_free(certificate = True)
            ...      if cycle.order() % 2 == 1:
            ...          print("Error !")
            ...      if not cycle.is_isomorphic(
            ...             graphs.CycleGraph(cycle.order())):
            ...          print("Error !")
            ...
            sage: print("Everything is Fine !")
            Everything is Fine !

        TESTS:

        Bug reported in :trac:`9925`, and fixed by :trac:`9420`::

            sage: g = Graph(':SiBFGaCEF_@CE`DEGH`CEFGaCDGaCDEHaDEF`CEH`ABCDEF', loops=False, multiedges=False)
            sage: g.is_even_hole_free()
            False
            sage: g.is_even_hole_free(certificate = True)
            Subgraph of (): Graph on 4 vertices

        Making sure there are no other counter-examples around ::

            sage: t = lambda x : (Graph(x).is_forest() or
            ...         isinstance(Graph(x).is_even_hole_free(certificate = True),Graph))
            sage: all( t(graphs.RandomBipartite(10,10,.5)) for i in range(100) )
            True

        REFERENCE:

        .. [ABCHRS08] \L. Addario-Berry, M. Chudnovsky, F. Havet, B. Reed, P. Seymour
          Bisimplicial vertices in even-hole-free graphs
          Journal of Combinatorial Theory, Series B
          vol 98, n.6 pp 1119-1164, 2008
        """
        from sage.graphs.graph_generators import GraphGenerators

        girth = self.girth()

        if girth > self.order():
            start = 4

        elif girth % 2 == 0:
            if not certificate:
                return False
            start = girth

        else:
            start = girth + 1

        while start <= self.order():


            subgraph = self.subgraph_search(GraphGenerators().CycleGraph(start), induced = True)

            if not subgraph is None:
                if certificate:
                    return subgraph
                else:
                    return False

            start = start + 2

        return True

    @doc_index("Graph properties")
    def is_odd_hole_free(self, certificate = False):
        r"""
        Tests whether ``self`` contains an induced odd hole.

        A Hole is a cycle of length at least 4 (included). It is said
        to be even (resp. odd) if its length is even (resp. odd).

        It is interesting to notice that while it is polynomial to
        check whether a graph has an odd hole or an odd antihole [CRST06]_, it is
        not known whether testing for one of these two cases
        independently is polynomial too.

        INPUT:

        - ``certificate`` (boolean) -- When ``certificate = False``,
          this method only returns ``True`` or ``False``. If
          ``certificate = True``, the subgraph found is returned
          instead of ``False``.

        EXAMPLE:

        Is the Petersen Graph odd-hole-free ::

            sage: g = graphs.PetersenGraph()
            sage: g.is_odd_hole_free()
            False

        Which was to be expected, as its girth is 5 ::

            sage: g.girth()
            5

        We can check the certificate returned is indeed a 5-cycle::

            sage: cycle = g.is_odd_hole_free(certificate = True)
            sage: cycle.is_isomorphic(graphs.CycleGraph(5))
            True

        As any chordal graph is hole-free, no interval graph has an odd hole::

            sage: g = graphs.RandomIntervalGraph(20)
            sage: g.is_odd_hole_free()
            True

        REFERENCES:

        .. [CRST06] \M. Chudnovsky, G. Cornuejols, X. Liu, P. Seymour, K. Vuskovic
          Recognizing berge graphs
          Combinatorica vol 25, n 2, pages 143--186
          2005
        """
        from sage.graphs.graph_generators import GraphGenerators

        girth = self.odd_girth()

        if girth > self.order():
            return True
        if girth == 3:
            start = 5

        else:
            if not certificate:
                return False
            start = girth

        while start <= self.order():

            subgraph = self.subgraph_search(GraphGenerators().CycleGraph(start), induced = True)

            if not subgraph is None:
                if certificate:
                    return subgraph
                else:
                    return False

            start += 2

        return True

    @doc_index("Graph properties")
    def is_bipartite(self, certificate = False):
        """
        Returns ``True`` if graph `G` is bipartite, ``False`` if not.

        Traverse the graph G with breadth-first-search and color nodes.

        INPUT:

        - ``certificate`` -- whether to return a certificate (``False`` by
          default). If set to ``True``, the certificate returned in a proper
          2-coloring when `G` is bipartite, and an odd cycle otherwise.

        EXAMPLES::

            sage: graphs.CycleGraph(4).is_bipartite()
            True
            sage: graphs.CycleGraph(5).is_bipartite()
            False
            sage: graphs.RandomBipartite(100,100,0.7).is_bipartite()
            True

        A random graph is very rarely bipartite::

            sage: g = graphs.PetersenGraph()
            sage: g.is_bipartite()
            False
            sage: false, oddcycle = g.is_bipartite(certificate = True)
            sage: len(oddcycle) % 2
            1
        """
        color = {}

        # For any uncolored vertex in the graph (to ensure we do the right job
        # when the graph is not connected !)
        for u in self:
            if u in color:
                continue

            # Let us run a BFS starting from u
            queue = [u]
            color[u] = 1
            while queue:
                v = queue.pop(0)
                c = 1-color[v]
                for w in self.neighbor_iterator(v):

                    # If the vertex has already been colored
                    if w in color:

                        # The graph is not bipartite !
                        if color[w] == color[v]:

                            # Should we return an odd cycle ?
                            if certificate:

                                # We build the first half of the cycle, i.e. a
                                # u-w path
                                cycle = self.shortest_path(u,w)

                                # The second half is a v-u path, but there may
                                # be common vertices in the two paths. But we
                                # can avoid that !

                                for v in self.shortest_path(v,u):
                                    if v in cycle:
                                        return False, cycle[cycle.index(v):]
                                    else:
                                        cycle.append(v)
                            else:
                                return False

                    # We color a new vertex
                    else:
                        color[w] = c
                        queue.append(w)
        if certificate:
            return True, color
        else:
            return True

    @doc_index("Graph properties")
    def is_triangle_free(self, algorithm='bitset'):
        r"""
        Returns whether ``self`` is triangle-free

        INPUT:

        - ``algorithm`` -- (default: ``'bitset'``) specifies the algorithm to
          use among:

          - ``'matrix'`` -- tests if the trace of the adjacency matrix is
            positive.

          - ``'bitset'`` -- encodes adjacencies into bitsets and uses fast
            bitset operations to test if the input graph contains a
            triangle. This method is generaly faster than stantard matrix
            multiplication.

        EXAMPLE:

        The Petersen Graph is triangle-free::

            sage: g = graphs.PetersenGraph()
            sage: g.is_triangle_free()
            True

        or a complete Bipartite Graph::

            sage: G = graphs.CompleteBipartiteGraph(5,6)
            sage: G.is_triangle_free(algorithm='matrix')
            True
            sage: G.is_triangle_free(algorithm='bitset')
            True

        a tripartite graph, though, contains many triangles::

            sage: G = (3 * graphs.CompleteGraph(5)).complement()
            sage: G.is_triangle_free(algorithm='matrix')
            False
            sage: G.is_triangle_free(algorithm='bitset')
            False

        TESTS:

        Comparison of algorithms::

            sage: for i in xrange(10): # long test
            ...       G = graphs.RandomBarabasiAlbert(50,2)
            ...       bm = G.is_triangle_free(algorithm='matrix')
            ...       bb = G.is_triangle_free(algorithm='bitset')
            ...       if bm != bb:
            ...          print("That's not good!")

        Asking for an unknown algorithm::

            sage: g.is_triangle_free(algorithm='tip top')
            Traceback (most recent call last):
            ...
            ValueError: Algorithm 'tip top' not yet implemented. Please contribute.
        """
        if algorithm=='bitset':
            from sage.data_structures.bitset import Bitset
            N = self.num_verts()
            map = {}
            i = 0
            B = {}
            for u in self.vertex_iterator():
                map[u] = i
                i += 1
                B[u] = Bitset(capacity=N)
            # map adjacency to bitsets
            for u,v in self.edge_iterator(labels=None):
                B[u].add(map[v])
                B[v].add(map[u])
            # map lengths 2 paths to bitsets
            BB = Bitset(capacity=N)
            for u in self.vertex_iterator():
                BB.clear()
                for v in self.vertex_iterator():
                    if B[u]&B[v]:
                        BB.add(map[v])
                # search for triangles
                if B[u]&BB:
                    return False
            return True

        elif algorithm=='matrix':
            return (self.adjacency_matrix()**3).trace() == 0

        else:
            raise ValueError("Algorithm '%s' not yet implemented. Please contribute." %(algorithm))

    @doc_index("Graph properties")
    def is_split(self):
        r"""
        Returns ``True`` if the graph is a Split graph, ``False`` otherwise.

        A Graph `G` is said to be a split graph if its vertices `V(G)`
        can be partitioned into two sets `K` and `I` such that the
        vertices of `K` induce a complete graph, and those of `I` are
        an independent set.

        There is a simple test to check whether a graph is a split
        graph (see, for instance, the book "Graph Classes, a survey"
        [GraphClasses]_ page 203) :

        Given the degree sequence `d_1 \geq ... \geq d_n` of `G`, a graph
        is a split graph if and only if :

        .. MATH::

            \sum_{i=1}^\omega d_i = \omega (\omega - 1) + \sum_{i=\omega + 1}^nd_i

        where `\omega = max \{i:d_i\geq i-1\}`.


        EXAMPLES:

        Split graphs are, in particular, chordal graphs. Hence, The Petersen graph
        can not be split::

            sage: graphs.PetersenGraph().is_split()
            False

        We can easily build some "random" split graph by creating a
        complete graph, and adding vertices only connected
        to some random vertices of the clique::

            sage: g = graphs.CompleteGraph(10)
            sage: sets = Subsets(Set(range(10)))
            sage: for i in range(10, 25):
            ...      g.add_edges([(i,k) for k in sets.random_element()])
            sage: g.is_split()
            True

        Another caracterisation of split graph states that a graph is a split graph
        if and only if does not contain the 4-cycle, 5-cycle or 2K_2 as an induced
        subgraph. Hence for the above graph we have::

            sage: sum([g.subgraph_search_count(H,induced=True) for H in [graphs.CycleGraph(4),graphs.CycleGraph(5), 2*graphs.CompleteGraph(2)]])
            0


        REFERENCES:

        .. [GraphClasses] \A. Brandstadt, VB Le and JP Spinrad
          Graph classes: a survey
          SIAM Monographs on Discrete Mathematics and Applications},
          1999
        """
        self._scream_if_not_simple()
        # our degree sequence is numbered from 0 to n-1, so to avoid
        # any mistake, let's fix it :-)
        degree_sequence = [0] + sorted(self.degree(), reverse = True)

        for (i, d) in enumerate(degree_sequence):
            if d >= i - 1:
                omega = i
            else:
                break

        left = sum(degree_sequence[:omega + 1])
        right = omega * (omega - 1) + sum(degree_sequence[omega + 1:])

        return left == right

    @doc_index("Algorithmically hard stuff")
    def treewidth(self,k=None,certificate=False,algorithm=None):
        r"""
        Computes the tree-width of `G` (and provides a decomposition)

        INPUT:

        - ``k`` (integer) -- the width to be considered. When ``k`` is an
          integer, the method checks that the graph has treewidth `\leq k`. If
          ``k`` is ``None`` (default), the method computes the optimal
          tree-width.

        - ``certificate`` -- whether to return the tree-decomposition itself.

        - ``algorithm`` -- whether to use ``"sage"`` or ``"tdlib"`` (requires
          the installation of the 'tdlib' package). The default behaviour is to
          use 'tdlib' if it is available, and Sage's own algorithm when it is
          not.

        OUTPUT:

            ``g.treewidth()`` returns the treewidth of ``g``. When ``k`` is
             specified, it returns ``False`` when no tree-decomposition of width
             `\leq k` exists or ``True`` otherwise. When ``certificate=True``,
             the tree-decomposition is also returned.

        ALGORITHM:

            This function virtually explores the graph of all pairs
            ``(vertex_cut,cc)``, where ``vertex_cut`` is a vertex cut of the
            graph of cardinality `\leq k+1`, and ``connected_component`` is a
            connected component of the graph induced by ``G-vertex_cut``.

            We deduce that the pair ``(vertex_cut,cc)`` is feasible with
            tree-width `k` if ``cc`` is empty, or if a vertex ``v`` from
            ``vertex_cut`` can be replaced with a vertex from ``cc``, such that
            the pair ``(vertex_cut+v,cc-v)`` is feasible.

        .. NOTE::

            The implementation would be much faster if ``cc``, the argument of the
            recursive function, was a bitset. It would also be very nice to not copy
            the graph in order to compute connected components, for this is really a
            waste of time.

        .. SEEALSO::

            :meth:`~sage.graphs.graph_decompositions.vertex_separation.path_decomposition`
            computes the pathwidth of a graph. See also the
            :mod:`~sage.graphs.graph_decompositions.vertex_separation` module.

        EXAMPLES:

        The PetersenGraph has treewidth 4::

            sage: graphs.PetersenGraph().treewidth()
            4
            sage: graphs.PetersenGraph().treewidth(certificate=True)
            Tree decomposition: Graph on 6 vertices

        The treewidth of a 2d grid is its smallest side::

            sage: graphs.Grid2dGraph(2,5).treewidth()
            2
            sage: graphs.Grid2dGraph(3,5).treewidth()
            3

        TESTS::

            sage: g = graphs.PathGraph(3)
            sage: g.treewidth()
            1
            sage: g = 2*graphs.PathGraph(3)
            sage: g.treewidth()
            1
            sage: g.treewidth(certificate=True)
            Tree decomposition: Graph on 4 vertices
            sage: g.treewidth(2)
            True
            sage: g.treewidth(1)
            True
            sage: Graph(1).treewidth()
            0
            sage: Graph(0).treewidth()
            -1
            sage: graphs.PetersenGraph().treewidth(k=2)
            False
            sage: graphs.PetersenGraph().treewidth(k=6)
            True
            sage: graphs.PetersenGraph().treewidth(certificate=True).is_tree()
            True
            sage: graphs.PetersenGraph().treewidth(k=3,certificate=True)
            False
            sage: graphs.PetersenGraph().treewidth(k=4,certificate=True)
            Tree decomposition: Graph on 6 vertices

        All edges do appear (:trac:`17893`)::

            sage: from itertools import combinations
            sage: g = graphs.PathGraph(10)
            sage: td = g.treewidth(certificate=True)
            sage: for bag in td:
            ....:    g.delete_edges(list(combinations(bag,2)))
            sage: g.size()
            0

        :trac:`19358`::

            sage: g = Graph()
            sage: for i in range(3):
            ....:     for j in range(2):
            ....:         g.add_path([i,(i,j),(i+1)%3])
            sage: g.treewidth()
            2

        Trivially true::

            sage: graphs.PetersenGraph().treewidth(k=35)
            True
            sage: graphs.PetersenGraph().treewidth(k=35,certificate=True)
            Tree decomposition: Graph on 1 vertex

        Bad input:

            sage: graphs.PetersenGraph().treewidth(k=-3)
            Traceback (most recent call last):
            ...
            ValueError: k(=-3) must be a nonnegative integer
        """
        g = self

        # Check Input
        if algorithm is None:
            try:
                import sage.graphs.graph_decompositions.tdlib as tdlib
                algorithm = "tdlib"
            except ImportError:
                algorithm = "sage"
        elif (algorithm != "sage"   and
              algorithm != "tdlib"):
            raise ValueError("'algorithm' must be equal to 'tdlib', 'sage', or None")

        if k is not None and k<0:
            raise ValueError("k(={}) must be a nonnegative integer".format(k))

        # Stupid cases
        if g.order() == 0:
            if certificate: return Graph()
            elif k is None: return -1
            else:           return True

        if k is not None and k >= g.order()-1:
            if certificate:
                return Graph({sage.sets.set.Set(g.vertices()):[]},
                             name="Tree decomposition")
            return True

        # TDLIB
        if algorithm == 'tdlib':
            try:
                import sage.graphs.graph_decompositions.tdlib as tdlib
            except ImportError:
                from sage.misc.package import PackageNotFoundError
                raise PackageNotFoundError("tdlib")

            T = tdlib.treedecomposition_exact(g, -1 if k is None else k)
            width = tdlib.get_width(T)

            if certificate:
                return T if (k is None or width <= k) else False
            elif k is None:
                return width
            else:
                return (width <= k)

        # Disconnected cases
        if not g.is_connected():
            if certificate is False:
                if k is None:
                    return max(cc.treewidth() for cc in g.connected_components_subgraphs())
                else:
                    return all(cc.treewidth(k) for cc in g.connected_components_subgraphs())
            else:
                return Graph(sum([cc.treewidth(certificate=True).edges(labels=False)
                                  for cc in g.connected_components_subgraphs()],[]),
                             name="Tree decomposition")

        # Forcing k to be defined
        if k is None:
            for i in range(max(0,g.clique_number()-1,min(g.degree())),
                           g.order()+1):
                ans = g.treewidth(k=i, certificate=certificate)
                if ans:
                    return ans if certificate else i

        # This is the recursion described in the method's documentation. All
        # computations are cached, and depends on the pair ``cut,
        # connected_component`` only.
        #
        # It returns either a boolean or the corresponding tree-decomposition, as a
        # list of edges between vertex cuts (as it is done for the complete
        # tree-decomposition at the end of the main function.
        from sage.misc.cachefunc import cached_function
        @cached_function
        def rec(cut,cc):
            # Easy cases
            if len(cut) > k:
                return False
            if len(cc)+len(cut) <= k+1:
                return [(cut,cut.union(cc))] if certificate else True

            # We explore all possible extensions of the cut
            for v in cc:

                # New cuts and connected components, with v respectively added and
                # removed
                cutv = cut.union([v])
                ccv = cc.difference([v])

                # The values returned by the recursive calls.
                sons = []

                # Removing v may have disconnected cc. We iterate on its connected
                # components
                for cci in g.subgraph(ccv).connected_components():

                    # The recursive subcalls. We remove on-the-fly the vertices from
                    # the cut which play no role in separating the connected
                    # component from the rest of the graph.
                    reduced_cut = frozenset([x for x in cutv if any(xx in cci for xx in g.neighbors(x))])
                    son = rec(reduced_cut,frozenset(cci))
                    if son is False:
                        break

                    if certificate:
                        sons.extend(son)
                        sons.append((cut,cutv))
                        sons.append((cutv,reduced_cut))

                # Weird Python syntax which is useful once in a lifetime : if break
                # was never called in the loop above, we return "sons".
                else:
                    return sons if certificate else True

            return False

        # Main call to rec function, i.e. rec({v},V-{v})
        V = g.vertices()
        v = frozenset([V.pop(0)])
        TD = rec(v,frozenset(V))

        if TD is False:
            return False

        if not certificate:
            return True

        # Building the Tree-Decomposition graph. Its vertices are cuts of the
        # decomposition, and there is an edge from a cut C1 to a cut C2 if C2 is an
        # immediate subcall of C1
        from sage.sets.set import Set
        G = Graph(name="Tree decomposition")
        G.add_edges([(Set(x),Set(y)) for x,y in TD])

        # The Tree-Decomposition contains a lot of useless nodes.
        #
        # We merge all edges between two sets S,S' where S is a subset of S'
        changed = True
        while changed:
            changed=False
            for v in G.vertices():
                for u in G.neighbors(v):
                    if u.issuperset(v):
                        G.merge_vertices([u,v]) # the new vertex is named 'u'
                        changed = True
                        break

        return G

    @doc_index("Algorithmically hard stuff")
    def is_perfect(self, certificate = False):
        r"""
        Tests whether the graph is perfect.

        A graph `G` is said to be perfect if `\chi(H)=\omega(H)` hold
        for any induced subgraph `H\subseteq_i G` (and so for `G`
        itself, too), where `\chi(H)` represents the chromatic number
        of `H`, and `\omega(H)` its clique number. The Strong Perfect
        Graph Theorem [SPGT]_ gives another characterization of
        perfect graphs:

        A graph is perfect if and only if it contains no odd hole
        (cycle on an odd number `k` of vertices, `k>3`) nor any odd
        antihole (complement of a hole) as an induced subgraph.

        INPUT:

        - ``certificate`` (boolean) -- whether to return
          a certificate (default : ``False``)

        OUTPUT:

        When ``certificate = False``, this function returns
        a boolean value. When ``certificate = True``, it returns
        a subgraph of ``self`` isomorphic to an odd hole or an odd
        antihole if any, and ``None`` otherwise.

        EXAMPLE:

        A Bipartite Graph is always perfect ::

            sage: g = graphs.RandomBipartite(8,4,.5)
            sage: g.is_perfect()
            True

        So is the line graph of a bipartite graph::

            sage: g = graphs.RandomBipartite(4,3,0.7)
            sage: g.line_graph().is_perfect() # long time
            True

        As well as the Cartesian product of two complete graphs::

            sage: g = graphs.CompleteGraph(3).cartesian_product(graphs.CompleteGraph(3))
            sage: g.is_perfect()
            True

        Interval Graphs, which are chordal graphs, too ::

            sage: g =  graphs.RandomIntervalGraph(7)
            sage: g.is_perfect()
            True

        The PetersenGraph, which is triangle-free and
        has chromatic number 3 is obviously not perfect::

            sage: g = graphs.PetersenGraph()
            sage: g.is_perfect()
            False

        We can obtain an induced 5-cycle as a certificate::

            sage: g.is_perfect(certificate = True)
            Subgraph of (Petersen graph): Graph on 5 vertices

        TEST:

        Check that :trac:`13546` has been fixed::

            sage: Graph(':FgGE@I@GxGs', loops=False, multiedges=False).is_perfect()
            False
            sage: g = Graph({0: [2, 3, 4, 5],
            ...              1: [3, 4, 5, 6],
            ...              2: [0, 4, 5, 6],
            ...              3: [0, 1, 5, 6],
            ...              4: [0, 1, 2, 6],
            ...              5: [0, 1, 2, 3],
            ...              6: [1, 2, 3, 4]})
            sage: g.is_perfect()
            False

        REFERENCES:

        .. [SPGT] \M. Chudnovsky, N. Robertson, P. Seymour, R. Thomas.
          The strong perfect graph theorem
          Annals of Mathematics
          vol 164, number 1, pages 51--230
          2006

        TESTS::

            sage: Graph(':Ab').is_perfect()
            Traceback (most recent call last):
            ...
            ValueError: This method is only defined for simple graphs, and yours is not one of them !
            sage: g = Graph()
            sage: g.allow_loops(True)
            sage: g.add_edge(0,0)
            sage: g.edges()
            [(0, 0, None)]
            sage: g.is_perfect()
            Traceback (most recent call last):
            ...
            ValueError: This method is only defined for simple graphs, and yours is not one of them !

        """

        if self.has_multiple_edges() or self.has_loops():
            raise ValueError("This method is only defined for simple graphs,"
                             " and yours is not one of them !")
        if self.is_bipartite():

            return True if not certificate else None

        self_complement = self.complement()

        self_complement.remove_loops()
        self_complement.remove_multiple_edges()

        if self_complement.is_bipartite():
            return True if not certificate else None

        answer = self.is_odd_hole_free(certificate = certificate)
        if not (answer is True):
            return answer

        return self_complement.is_odd_hole_free(certificate = certificate)

    @doc_index("Graph properties")
    def odd_girth(self):
        r"""
        Returns the odd girth of self.

        The odd girth of a graph is defined as the smallest cycle of odd length.

        OUTPUT:

        The odd girth of ``self``.

        EXAMPLES:

        The McGee graph has girth 7 and therefore its odd girth is 7 as well. ::

            sage: G = graphs.McGeeGraph()
            sage: G.odd_girth()
            7

        Any complete graph on more than 2 vertices contains a triangle and has
        thus odd girth 3. ::

            sage: G = graphs.CompleteGraph(10)
            sage: G.odd_girth()
            3

        Every bipartite graph has no odd cycles and consequently odd girth of
        infinity. ::

            sage: G = graphs.CompleteBipartiteGraph(100,100)
            sage: G.odd_girth()
            +Infinity

        .. SEEALSO::

            * :meth:`~sage.graphs.generic_graph.GenericGraph.girth` -- computes
              the girth of a graph.

        REFERENCES:

        The property relating the odd girth to the coefficients of the
        characteristic polynomial is an old result from algebraic graph theory
        see

        .. [Har62] Harary, F (1962). The determinant of the adjacency matrix of
          a graph, SIAM Review 4, 202-210

        .. [Biggs93] Biggs, N. L. Algebraic Graph Theory, 2nd ed. Cambridge,
          England: Cambridge University Press, pp. 45, 1993.

        TESTS::

            sage: graphs.CycleGraph(5).odd_girth()
            5
            sage: graphs.CycleGraph(11).odd_girth()
            11
        """
        ch = ((self.am()).charpoly()).coefficients(sparse=False)
        n = self.order()

        for i in xrange(n-1,-1,-2):
            if ch[i] != 0:
                return n-i

        from sage.rings.infinity import Infinity

        return Infinity

    @doc_index("Graph properties")
    def is_edge_transitive(self):
        """
        Returns true if self is an edge transitive graph.

        A graph is edge-transitive if its automorphism group acts transitively
        on its edge set.

        Equivalently, if there exists for any pair of edges `uv,u'v'\in E(G)` an
        automorphism `\phi` of `G` such that `\phi(uv)=u'v'` (note this does not
        necessarily mean that `\phi(u)=u'` and `\phi(v)=v'`).

        See :wikipedia:`the wikipedia article on edge-transitive graphs
        <Edge-transitive_graph>` for more information.

        .. SEEALSO::

          - :meth:`~Graph.is_arc_transitive`
          - :meth:`~Graph.is_half_transitive`
          - :meth:`~Graph.is_semi_symmetric`

        EXAMPLES::

            sage: P = graphs.PetersenGraph()
            sage: P.is_edge_transitive()
            True
            sage: C = graphs.CubeGraph(3)
            sage: C.is_edge_transitive()
            True
            sage: G = graphs.GrayGraph()
            sage: G.is_edge_transitive()
            True
            sage: P = graphs.PathGraph(4)
            sage: P.is_edge_transitive()
            False
        """
        from sage.interfaces.gap import gap

        if self.size() == 0:
            return True

        A = self.automorphism_group()
        e = next(self.edge_iterator(labels=False))
        e = [A._domain_to_gap[e[0]], A._domain_to_gap[e[1]]]

        return gap("OrbitLength("+str(A._gap_())+",Set(" + str(e) + "),OnSets);") == self.size()

    @doc_index("Graph properties")
    def is_arc_transitive(self):
        """
        Returns true if self is an arc-transitive graph

        A graph is arc-transitive if its automorphism group acts transitively on
        its pairs of adjacent vertices.

        Equivalently, if there exists for any pair of edges `uv,u'v'\in E(G)` an
        automorphism `\phi_1` of `G` such that `\phi_1(u)=u'` and
        `\phi_1(v)=v'`, as well as another automorphism `\phi_2` of `G` such
        that `\phi_2(u)=v'` and `\phi_2(v)=u'`

        See :wikipedia:`the wikipedia article on arc-transitive graphs
        <arc-transitive_graph>` for more information.

        .. SEEALSO::

          - :meth:`~Graph.is_edge_transitive`
          - :meth:`~Graph.is_half_transitive`
          - :meth:`~Graph.is_semi_symmetric`

        EXAMPLES::

            sage: P = graphs.PetersenGraph()
            sage: P.is_arc_transitive()
            True
            sage: G = graphs.GrayGraph()
            sage: G.is_arc_transitive()
            False
        """

        from sage.interfaces.gap import gap

        if self.size() == 0:
            return True

        A = self.automorphism_group()
        e = next(self.edge_iterator(labels=False))
        e = [A._domain_to_gap[e[0]], A._domain_to_gap[e[1]]]

        return gap("OrbitLength("+str(A._gap_())+",Set(" + str(e) + "),OnTuples);") == 2*self.size()

    @doc_index("Graph properties")
    def is_half_transitive(self):
        """
        Returns true if self is a half-transitive graph.

        A graph is half-transitive if it is both vertex and edge transitive
        but not arc-transitive.

        See :wikipedia:`the wikipedia article on half-transitive graphs
        <half-transitive_graph>` for more information.

        .. SEEALSO::

          - :meth:`~Graph.is_edge_transitive`
          - :meth:`~Graph.is_arc_transitive`
          - :meth:`~Graph.is_semi_symmetric`

        EXAMPLES:

        The Petersen Graph is not half-transitive::

            sage: P = graphs.PetersenGraph()
            sage: P.is_half_transitive()
            False

        The smallest half-transitive graph is the Holt Graph::

            sage: H = graphs.HoltGraph()
            sage: H.is_half_transitive()
            True
        """

        # A half-transitive graph always has only vertices of even degree
        if not all(d%2 == 0 for d in self.degree_iterator()):
            return False

        return (self.is_edge_transitive() and
                self.is_vertex_transitive() and
                not self.is_arc_transitive())

    @doc_index("Graph properties")
    def is_semi_symmetric(self):
        """
        Returns true if self is semi-symmetric.

        A graph is semi-symmetric if it is regular, edge-transitve but not
        vertex-transitive.

        See :wikipedia:`the wikipedia article on semi-symmetric graphs
        <Semi-symmetric_graph>` for more information.

        .. SEEALSO::

          - :meth:`~Graph.is_edge_transitive`
          - :meth:`~Graph.is_arc_transitive`
          - :meth:`~Graph.is_half_transitive`

        EXAMPLES:

        The Petersen graph is not semi-symmetric::

            sage: P = graphs.PetersenGraph()
            sage: P.is_semi_symmetric()
            False

        The Gray graph is the smallest possible cubic semi-symmetric graph::

            sage: G = graphs.GrayGraph()
            sage: G.is_semi_symmetric()
            True

        Another well known semi-symmetric graph is the Ljubljana graph::

            sage: L = graphs.LjubljanaGraph()
            sage: L.is_semi_symmetric()
            True
        """
        # A semi-symmetric graph is always bipartite
        if not self.is_bipartite():
            return False

        return (self.is_regular() and
                self.is_edge_transitive() and not
                self.is_vertex_transitive())

    @doc_index("Connectivity, orientations, trees")
    def degree_constrained_subgraph(self, bounds=None, solver=None, verbose=0):
        r"""
        Returns a degree-constrained subgraph.

        Given a graph `G` and two functions `f, g:V(G)\rightarrow \mathbb Z`
        such that `f \leq g`, a degree-constrained subgraph in `G` is
        a subgraph `G' \subseteq G` such that for any vertex `v \in G`,
        `f(v) \leq d_{G'}(v) \leq g(v)`.

        INPUT:

        - ``bounds`` -- (default: ``None``) Two possibilities:

          - A dictionary whose keys are the vertices, and values a pair of
            real values ``(min,max)`` corresponding to the values
            `(f(v),g(v))`.

          - A function associating to each vertex a pair of
            real values ``(min,max)`` corresponding to the values
            `(f(v),g(v))`.


        - ``solver`` -- (default: ``None``) Specify a Linear Program (LP)
          solver to be used. If set to ``None``, the default one is used. For
          more information on LP solvers and which default solver is used, see
          the method
          :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>`
          of the class
          :class:`MixedIntegerLinearProgram <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``). Sets the level of
          verbosity. Set to 0 by default, which means quiet.

        OUTPUT:

        - When a solution exists, this method outputs the degree-constained
          subgraph as a Graph object.

        - When no solution exists, returns ``False``.

        .. NOTE::

            - This algorithm computes the degree-constrained subgraph of minimum weight.
            - If the graph's edges are weighted, these are taken into account.
            - This problem can be solved in polynomial time.

        EXAMPLES:

        Is there a perfect matching in an even cycle? ::

            sage: g = graphs.CycleGraph(6)
            sage: bounds = lambda x: [1,1]
            sage: m = g.degree_constrained_subgraph(bounds=bounds)
            sage: m.size()
            3
        """
        self._scream_if_not_simple()
        from sage.numerical.mip import MixedIntegerLinearProgram, MIPSolverException

        p = MixedIntegerLinearProgram(maximization=False, solver=solver)
        b = p.new_variable(binary=True)

        reorder = lambda x,y: (x,y) if x<y else (y,x)

        if bounds is None:
            raise ValueError("The `bounds` keyword can not be equal to None")
        elif isinstance(bounds,dict):
            f_bounds = lambda x: bounds[x]
        else:
            f_bounds = bounds


        if self.weighted():
            from sage.rings.real_mpfr import RR
            weight = lambda x: x if x in RR else 1
        else:
            weight = lambda x: 1

        for v in self:
            minimum,maximum = f_bounds(v)
            p.add_constraint(p.sum( b[reorder(x,y)]*weight(l) for x,y,l in self.edges_incident(v)), min=minimum, max=maximum)

        p.set_objective(p.sum( b[reorder(x,y)]*weight(l) for x,y,l in self.edge_iterator()))

        try:
            p.solve(log=verbose)
            g = copy(self)
            b = p.get_values(b)
            g.delete_edges([(x,y) for x,y,_ in g.edge_iterator() if b[reorder(x,y)] < 0.5])
            return g


        except MIPSolverException:
            return False


    ### Orientations

    @doc_index("Connectivity, orientations, trees")
    def strong_orientation(self):
        r"""
        Returns a strongly connected orientation of the current graph.

        An orientation of an undirected graph is a digraph obtained by giving an
        unique direction to each of its edges. An orientation is said to be
        strong if there is a directed path between each pair of vertices.  See
        also the :wikipedia:`Strongly_connected_component`.

        If the graph is 2-edge-connected, a strongly connected orientation
        can be found in linear time. If the given graph is not 2-connected,
        the orientation returned will ensure that each 2-connected component
        has a strongly connected orientation.

        OUTPUT:

        A digraph representing an orientation of the current graph.

        .. NOTE::

            - This method assumes the graph is connected.
            - This algorithm works in O(m).

        EXAMPLE:

        For a 2-regular graph, a strong orientation gives to each vertex
        an out-degree equal to 1::

            sage: g = graphs.CycleGraph(5)
            sage: g.strong_orientation().out_degree()
            [1, 1, 1, 1, 1]

        The Petersen Graph is 2-edge connected. It then has a strongly
        connected orientation::

            sage: g = graphs.PetersenGraph()
            sage: o = g.strong_orientation()
            sage: len(o.strongly_connected_components())
            1

        The same goes for the CubeGraph in any dimension ::

            sage: all(len(graphs.CubeGraph(i).strong_orientation().strongly_connected_components()) == 1 for i in xrange(2,6))
            True

        A multigraph also has a strong orientation ::

            sage: g = Graph([(1,2),(1,2)],multiedges=True)
            sage: g.strong_orientation()
            Multi-digraph on 2 vertices

        """
        from sage.graphs.all import DiGraph
        d = DiGraph(multiedges=self.allows_multiple_edges())

        id = {}
        i = 0

        # The algorithm works through a depth-first search. Any edge
        # used in the depth-first search is oriented in the direction
        # in which it has been used. All the other edges are oriented
        # backward

        v = next(self.vertex_iterator())
        seen = {}
        i=1

        # Time at which the vertices have been discovered
        seen[v] = i

        # indicates the stack of edges to explore
        next_ = self.edges_incident(v)

        while next_:
            e = next_.pop(-1)
            # We assume e[0] to be a `seen` vertex
            e = e if seen.get(e[0],False) is not False else (e[1],e[0],e[2])

            # If we discovered a new vertex
            if seen.get(e[1],False) is False:
                d.add_edge(e)
                next_.extend([ee for ee in self.edges_incident(e[1]) if (((e[0],e[1]) != (ee[0],ee[1])) and ((e[0],e[1]) != (ee[1],ee[0])))])
                i+=1
                seen[e[1]]=i

            # Else, we orient the edges backward
            else:
                if seen[e[0]] < seen[e[1]]:
                    d.add_edge((e[1],e[0],e[2]))
                else:
                    d.add_edge(e)

        # Case of multiple edges. If another edge has already been inserted, we add the new one
        # in the opposite direction.
        tmp = None
        for e in self.multiple_edges():
            if tmp == (e[0],e[1]):
                if d.has_edge(e[0],e[1]):
                    d.add_edge(e[1],e[0],e[2])
                else:
                    d.add_edge(e)
            tmp = (e[0],e[1])

        return d

    @doc_index("Connectivity, orientations, trees")
    def minimum_outdegree_orientation(self, use_edge_labels=False, solver=None, verbose=0):
        r"""
        Returns an orientation of ``self`` with the smallest possible maximum
        outdegree.

        Given a Graph `G`, it is polynomial to compute an orientation
        `D` of the edges of `G` such that the maximum out-degree in
        `D` is minimized. This problem, though, is NP-complete in the
        weighted case [AMOZ06]_.

        INPUT:

        - ``use_edge_labels`` -- boolean (default: ``False``)

          - When set to ``True``, uses edge labels as weights to
            compute the orientation and assumes a weight of `1`
            when there is no value available for a given edge.

          - When set to ``False`` (default), gives a weight of 1
            to all the edges.

        - ``solver`` -- (default: ``None``) Specify a Linear Program (LP)
          solver to be used. If set to ``None``, the default one is used. For
          more information on LP solvers and which default solver is used, see
          the method
          :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>`
          of the class
          :class:`MixedIntegerLinearProgram <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``). Sets the level of
          verbosity. Set to 0 by default, which means quiet.

        EXAMPLE:

        Given a complete bipartite graph `K_{n,m}`, the maximum out-degree
        of an optimal orientation is `\left\lceil \frac {nm} {n+m}\right\rceil`::

            sage: g = graphs.CompleteBipartiteGraph(3,4)
            sage: o = g.minimum_outdegree_orientation()
            sage: max(o.out_degree()) == ceil((4*3)/(3+4))
            True

        REFERENCES:

        .. [AMOZ06] Asahiro, Y. and Miyano, E. and Ono, H. and Zenmyo, K.
          Graph orientation algorithms to minimize the maximum outdegree
          Proceedings of the 12th Computing: The Australasian Theory Symposium
          Volume 51, page 20
          Australian Computer Society, Inc. 2006
        """
        self._scream_if_not_simple()
        if self.is_directed():
            raise ValueError("Cannot compute an orientation of a DiGraph. "+\
                                 "Please convert it to a Graph if you really mean it.")

        if use_edge_labels:
            from sage.rings.real_mpfr import RR
            weight = lambda u,v : self.edge_label(u,v) if self.edge_label(u,v) in RR else 1
        else:
            weight = lambda u,v : 1

        from sage.numerical.mip import MixedIntegerLinearProgram

        p = MixedIntegerLinearProgram(maximization=False, solver=solver)

        # The orientation of an edge is boolean
        # and indicates whether the edge uv
        # with u<v goes from u to v ( equal to 0 )
        # or from v to u ( equal to 1)
        orientation = p.new_variable(binary=True)

        degree = p.new_variable(nonnegative=True)

        # Whether an edge adjacent to a vertex u counts
        # positively or negatively
        outgoing = lambda u,v,variable : (1-variable) if u>v else variable

        for u in self:
            p.add_constraint(p.sum(weight(u,v)*outgoing(u,v,orientation[min(u,v),max(u,v)]) for v in self.neighbors(u))-degree['max'], max=0)

        p.set_objective(degree['max'])

        p.solve(log=verbose)

        orientation = p.get_values(orientation)

        # All the edges from self are doubled in O
        # ( one in each direction )
        from sage.graphs.digraph import DiGraph
        O = DiGraph(self)

        # Builds the list of edges that should be removed
        edges=[]

        for u,v in self.edge_iterator(labels=None):
            # assumes u<v
            if u>v:
                u,v=v,u

            if orientation[min(u,v),max(u,v)] == 1:
                edges.append((max(u,v),min(u,v)))
            else:
                edges.append((min(u,v),max(u,v)))

        O.delete_edges(edges)

        return O

    @doc_index("Connectivity, orientations, trees")
    def bounded_outdegree_orientation(self, bound):
        r"""
        Computes an orientation of ``self`` such that every vertex `v`
        has out-degree less than `b(v)`

        INPUT:

        - ``bound`` -- Maximum bound on the out-degree. Can be of
          three different types :

         * An integer `k`. In this case, computes an orientation
           whose maximum out-degree is less than `k`.

         * A dictionary associating to each vertex its associated
           maximum out-degree.

         * A function associating to each vertex its associated
           maximum out-degree.

        OUTPUT:

        A DiGraph representing the orientation if it exists. A
        ``ValueError`` exception is raised otherwise.

        ALGORITHM:

        The problem is solved through a maximum flow :

        Given a graph `G`, we create a ``DiGraph`` `D` defined on
        `E(G)\cup V(G)\cup \{s,t\}`. We then link `s` to all of `V(G)`
        (these edges having a capacity equal to the bound associated
        to each element of `V(G)`), and all the elements of `E(G)` to
        `t` . We then link each `v \in V(G)` to each of its incident
        edges in `G`. A maximum integer flow of value `|E(G)|`
        corresponds to an admissible orientation of `G`. Otherwise,
        none exists.

        EXAMPLES:

        There is always an orientation of a graph `G` such that a
        vertex `v` has out-degree at most `\lceil \frac {d(v)} 2
        \rceil`::

            sage: g = graphs.RandomGNP(40, .4)
            sage: b = lambda v : ceil(g.degree(v)/2)
            sage: D = g.bounded_outdegree_orientation(b)
            sage: all( D.out_degree(v) <= b(v) for v in g )
            True


        Chvatal's graph, being 4-regular, can be oriented in such a
        way that its maximum out-degree is 2::

            sage: g = graphs.ChvatalGraph()
            sage: D = g.bounded_outdegree_orientation(2)
            sage: max(D.out_degree())
            2

        For any graph `G`, it is possible to compute an orientation
        such that the maximum out-degree is at most the maximum
        average degree of `G` divided by 2. Anything less, though, is
        impossible.

            sage: g = graphs.RandomGNP(40, .4)
            sage: mad = g.maximum_average_degree()

        Hence this is possible ::

            sage: d = g.bounded_outdegree_orientation(ceil(mad/2))

        While this is not::

            sage: try:
            ...      g.bounded_outdegree_orientation(ceil(mad/2-1))
            ...      print("Error")
            ... except ValueError:
            ...       pass

        TESTS:

        As previously for random graphs, but more intensively::

            sage: for i in xrange(30):      # long time (up to 6s on sage.math, 2012)
            ...       g = graphs.RandomGNP(40, .4)
            ...       b = lambda v : ceil(g.degree(v)/2)
            ...       D = g.bounded_outdegree_orientation(b)
            ...       if not (
            ...            all( D.out_degree(v) <= b(v) for v in g ) or
            ...            D.size() != g.size()):
            ...           print("Something wrong happened")

        """
        self._scream_if_not_simple()
        from sage.graphs.all import DiGraph
        n = self.order()

        if n == 0:
            return DiGraph()

        vertices = self.vertices()
        vertices_id = dict((y, x) for x,y in enumerate(vertices))

        b = {}


        # Checking the input type. We make a dictionary out of it
        if isinstance(bound, dict):
            b = bound
        else:
            try:
                b = dict(zip(vertices,map(bound, vertices)))

            except TypeError:
                b = dict(zip(vertices, [bound]*n))

        d = DiGraph()

        # Adding the edges (s,v) and ((u,v),t)
        d.add_edges( ('s', vertices_id[v], b[v]) for v in vertices)

        d.add_edges( ((vertices_id[u], vertices_id[v]), 't', 1)
                     for (u,v) in self.edges(labels=None) )

        # each v is linked to its incident edges

        for u,v in self.edges(labels = None):
            u,v = vertices_id[u], vertices_id[v]
            d.add_edge(u, (u,v), 1)
            d.add_edge(v, (u,v), 1)

        # Solving the maximum flow
        value, flow = d.flow('s','t', value_only = False, integer = True, use_edge_labels = True)

        if value != self.size():
            raise ValueError("No orientation exists for the given bound")

        D = DiGraph()
        D.add_vertices(vertices)

        # The flow graph may not contain all the vertices, if they are
        # not part of the flow...

        for u in [x for x in range(n) if x in flow]:

            for (uu,vv) in flow.neighbors_out(u):
                v = vv if vv != u else uu
                D.add_edge(vertices[u], vertices[v])

        # I do not like when a method destroys the embedding ;-)

        D.set_pos(self.get_pos())

        return D


    ### Coloring

    @doc_index("Basic methods")
    def bipartite_color(self):
        """
        Returns a dictionary with vertices as the keys and the color class
        as the values. Fails with an error if the graph is not bipartite.

        EXAMPLES::

            sage: graphs.CycleGraph(4).bipartite_color()
            {0: 1, 1: 0, 2: 1, 3: 0}
            sage: graphs.CycleGraph(5).bipartite_color()
            Traceback (most recent call last):
            ...
            RuntimeError: Graph is not bipartite.
        """
        isit, certificate = self.is_bipartite(certificate = True)

        if isit:
            return certificate
        else:
            raise RuntimeError("Graph is not bipartite.")

    @doc_index("Basic methods")
    def bipartite_sets(self):
        """
        Returns `(X,Y)` where `X` and `Y` are the nodes in each bipartite set of
        graph `G`. Fails with an error if graph is not bipartite.

        EXAMPLES::

            sage: graphs.CycleGraph(4).bipartite_sets()
            ({0, 2}, {1, 3})
            sage: graphs.CycleGraph(5).bipartite_sets()
            Traceback (most recent call last):
            ...
            RuntimeError: Graph is not bipartite.
        """
        color = self.bipartite_color()
        left = set([])
        right = set([])

        for u,s in color.iteritems():
            if s:
                left.add(u)
            else:
                right.add(u)

        return left, right

    @doc_index("Algorithmically hard stuff")
    def chromatic_number(self, algorithm="DLX", verbose = 0):
        r"""
        Returns the minimal number of colors needed to color the vertices
        of the graph `G`.

        INPUT:

        - ``algorithm`` -- Select an algorithm from the following supported
          algorithms:

          - If ``algorithm="DLX"`` (default), the chromatic number is
            computed using the dancing link algorithm. It is
            inefficient speedwise to compute the chromatic number through
            the dancing link algorithm because this algorithm computes
            *all* the possible colorings to check that one exists.

          - If ``algorithm="CP"``, the chromatic number is computed
            using the coefficients of the chromatic polynomial. Again, this
            method is inefficient in terms of speed and it only useful for
            small graphs.

          - If ``algorithm="MILP"``, the chromatic number is computed using a
            mixed integer linear program. The performance of this implementation
            is affected by whether optional MILP solvers have been installed
            (see the :mod:`MILP module <sage.numerical.mip>`, or Sage's tutorial
            on Linear Programming).

        - ``verbose`` -- integer (default: ``0``). Sets the level of verbosity
          for the MILP algorithm. Its default value is 0, which means *quiet*.

        .. SEEALSO::

            For more functions related to graph coloring, see the
            module :mod:`sage.graphs.graph_coloring`.

        EXAMPLES::

            sage: G = Graph({0: [1, 2, 3], 1: [2]})
            sage: G.chromatic_number(algorithm="DLX")
            3
            sage: G.chromatic_number(algorithm="MILP")
            3
            sage: G.chromatic_number(algorithm="CP")
            3

        A bipartite graph has (by definition) chromatic number 2::

            sage: graphs.RandomBipartite(50,50,0.7).chromatic_number()
            2

        A complete multipartite graph with k parts has chromatic number k::

            sage: all(graphs.CompleteMultipartiteGraph([5]*i).chromatic_number() == i for i in xrange(2,5))
            True

        The complete graph has the largest chromatic number from all the graphs
        of order n. Namely its chromatic number is n::

            sage: all(graphs.CompleteGraph(i).chromatic_number() == i for i in xrange(10))
            True

        The Kneser graph with parameters (n,2) for n > 3 has chromatic number n-2::

            sage: all(graphs.KneserGraph(i,2).chromatic_number() == i-2 for i in xrange(4,6))
            True

        A snark has chromatic index 4 hence its line graph has chromatic number 4::

            sage: graphs.FlowerSnark().line_graph().chromatic_number()
            4

        TESTS::

            sage: G = Graph({0: [1, 2, 3], 1: [2]})
            sage: G.chromatic_number(algorithm="foo")
            Traceback (most recent call last):
            ...
            ValueError: The 'algorithm' keyword must be set to either 'DLX', 'MILP' or 'CP'.
        """
        self._scream_if_not_simple(allow_multiple_edges=True)
        # default built-in algorithm; bad performance
        if algorithm == "DLX":
            from sage.graphs.graph_coloring import chromatic_number
            return chromatic_number(self)
        # Algorithm with good performance, but requires an optional
        # package: choose any of GLPK or CBC.
        elif algorithm == "MILP":
            from sage.graphs.graph_coloring import vertex_coloring
            return vertex_coloring(self, value_only=True, verbose = verbose)
        # another algorithm with bad performance; only good for small graphs
        elif algorithm == "CP":
            f = self.chromatic_polynomial()
            i = 0
            while f(i) == 0:
                i += 1
            return i
        else:
            raise ValueError("The 'algorithm' keyword must be set to either 'DLX', 'MILP' or 'CP'.")

    @doc_index("Algorithmically hard stuff")
    def coloring(self, algorithm="DLX", hex_colors=False, verbose = 0):
        r"""
        Returns the first (optimal) proper vertex-coloring found.

        INPUT:

        - ``algorithm`` -- Select an algorithm from the following supported
          algorithms:

          - If ``algorithm="DLX"`` (default), the coloring is computed using the
            dancing link algorithm.

          - If ``algorithm="MILP"``, the coloring is computed using a mixed
            integer linear program. The performance of this implementation is
            affected by whether optional MILP solvers have been installed (see
            the :mod:`MILP module <sage.numerical.mip>`).

        - ``hex_colors`` -- (default: ``False``) if ``True``, return a
          dictionary which can easily be used for plotting.

        - ``verbose`` -- integer (default: ``0``). Sets the level of verbosity
          for the MILP algorithm. Its default value is 0, which means *quiet*.

        .. SEEALSO::

            For more functions related to graph coloring, see the
            module :mod:`sage.graphs.graph_coloring`.

        EXAMPLES::

            sage: G = Graph("Fooba")
            sage: P = G.coloring(algorithm="MILP"); P
            [[2, 1, 3], [0, 6, 5], [4]]
            sage: P = G.coloring(algorithm="DLX"); P
            [[1, 2, 3], [0, 5, 6], [4]]
            sage: G.plot(partition=P)
            Graphics object consisting of 16 graphics primitives
            sage: H = G.coloring(hex_colors=True, algorithm="MILP")
            sage: for c in sorted(H.keys()):
            ....:     print("{} {}".format(c, H[c]))
            #0000ff [4]
            #00ff00 [0, 6, 5]
            #ff0000 [2, 1, 3]
            sage: H = G.coloring(hex_colors=True, algorithm="DLX")
            sage: for c in sorted(H.keys()):
            ....:     print("{} {}".format(c, H[c]))
            #0000ff [4]
            #00ff00 [1, 2, 3]
            #ff0000 [0, 5, 6]
            sage: G.plot(vertex_colors=H)
            Graphics object consisting of 16 graphics primitives

        .. PLOT::

            g = Graph("Fooba")
            sphinx_plot(g.plot(partition=g.coloring()))

        TESTS::

            sage: G.coloring(algorithm="foo")
            Traceback (most recent call last):
            ...
            ValueError: The 'algorithm' keyword must be set to either 'DLX' or 'MILP'.
        """
        self._scream_if_not_simple(allow_multiple_edges=True)
        if algorithm == "MILP":
            from sage.graphs.graph_coloring import vertex_coloring
            return vertex_coloring(self, hex_colors=hex_colors, verbose = verbose)
        elif algorithm == "DLX":
            from sage.graphs.graph_coloring import first_coloring
            return first_coloring(self, hex_colors=hex_colors)
        else:
            raise ValueError("The 'algorithm' keyword must be set to either 'DLX' or 'MILP'.")

    @doc_index("Algorithmically hard stuff")
    def chromatic_symmetric_function(self, R=None):
        r"""
        Return the chromatic symmetric function of ``self``.

        Let `G` be a graph. The chromatic symmetric function `X_G` was
        described in [Stanley95]_, specifically Theorem 2.5 states that

        .. MATH::

            X_G = \sum_{F \subseteq E(G)} (-1)^{|F|} p_{\lambda(F)},

        where `\lambda(F)` is the partition of the sizes of the connected
        components of the subgraph induced by the edges `F` and `p_{\mu}`
        is the powersum symmetric function.

        INPUT:

        - ``R`` -- (optional) the base ring for the symmetric functions;
          this uses `\ZZ` by default

        EXAMPLES::

            sage: s = SymmetricFunctions(ZZ).s()
            sage: G = graphs.CycleGraph(5)
            sage: XG = G.chromatic_symmetric_function(); XG
            p[1, 1, 1, 1, 1] - 5*p[2, 1, 1, 1] + 5*p[2, 2, 1]
             + 5*p[3, 1, 1] - 5*p[3, 2] - 5*p[4, 1] + 4*p[5]
            sage: s(XG)
            30*s[1, 1, 1, 1, 1] + 10*s[2, 1, 1, 1] + 10*s[2, 2, 1]

        Not all graphs have a positive Schur expansion::

            sage: G = graphs.ClawGraph()
            sage: XG = G.chromatic_symmetric_function(); XG
            p[1, 1, 1, 1] - 3*p[2, 1, 1] + 3*p[3, 1] - p[4]
            sage: s(XG)
            8*s[1, 1, 1, 1] + 5*s[2, 1, 1] - s[2, 2] + s[3, 1]

        We show that given a triangle `\{e_1, e_2, e_3\}`, we have
        `X_G = X_{G - e_1} + X_{G - e_2} - X_{G - e_1 - e_2}`::

            sage: G = Graph([[1,2],[1,3],[2,3]])
            sage: XG = G.chromatic_symmetric_function()
            sage: G1 = copy(G)
            sage: G1.delete_edge([1,2])
            sage: XG1 = G1.chromatic_symmetric_function()
            sage: G2 = copy(G)
            sage: G2.delete_edge([1,3])
            sage: XG2 = G2.chromatic_symmetric_function()
            sage: G3 = copy(G1)
            sage: G3.delete_edge([1,3])
            sage: XG3 = G3.chromatic_symmetric_function()
            sage: XG == XG1 + XG2 - XG3
            True

        REFERENCES:

        .. [Stanley95] \R. P. Stanley, *A symmetric function generalization
           of the chromatic polynomial of a graph*, Adv. Math., ***111***
           no.1 (1995), 166-194.
        """
        from sage.combinat.sf.sf import SymmetricFunctions
        from sage.combinat.partition import _Partitions
        from sage.misc.misc import powerset
        if R is None:
            R = ZZ
        p = SymmetricFunctions(R).p()
        ret = p.zero()
        for F in powerset(self.edges()):
            la = _Partitions(self.subgraph(edges=F).connected_components_sizes())
            ret += (-1)**len(F) * p[la]
        return ret

    @doc_index("Algorithmically hard stuff")
    def chromatic_quasisymmetric_function(self, t=None, R=None):
        r"""
        Return the chromatic quasisymmetric function of ``self``.

        Let `G` be a graph whose vertex set is totally ordered. The
        chromatic quasisymmetric function `X_G(t)` was first
        described in [SW12]_. We use the equivalent definition
        given in [BC15]_:

        .. MATH::

            X_G(t) = \sum_{\sigma=(\sigma_1,\ldots,\sigma_n)}
            t^{\operatorname{asc}(\sigma)}
            M_{|\sigma_1|,\ldots,|\sigma_n|},

        where we sum over all ordered set partitions of the vertex
        set of `G` such that each block `\sigma_i` is an independent
        (i.e., stable) set of `G`, and where
        `\operatorname{asc}(\sigma)` denotes the number of edges
        `\{u, v\}` of `G` such that `u < v` and `v` appears in a
        later part of `\sigma` than `u`.

        INPUT:

        - ``t`` -- (optional) the parameter `t`; uses the variable `t`
          in `\ZZ[t]` by default
        - ``R`` -- (optional) the base ring for the quasisymmetric
          functions; uses the parent of `t` by default

        EXAMPLES::

            sage: G = Graph([[1,2,3], [[1,3], [2,3]]])
            sage: G.chromatic_quasisymmetric_function()
            (2*t^2+2*t+2)*M[1, 1, 1] + M[1, 2] + t^2*M[2, 1]
            sage: G = graphs.PathGraph(4)
            sage: XG = G.chromatic_quasisymmetric_function(); XG
            (t^3+11*t^2+11*t+1)*M[1, 1, 1, 1] + (3*t^2+3*t)*M[1, 1, 2]
             + (3*t^2+3*t)*M[1, 2, 1] + (3*t^2+3*t)*M[2, 1, 1]
             + (t^2+t)*M[2, 2]
            sage: XG.to_symmetric_function()
            (t^3+11*t^2+11*t+1)*m[1, 1, 1, 1] + (3*t^2+3*t)*m[2, 1, 1]
             + (t^2+t)*m[2, 2]
            sage: G = graphs.CompleteGraph(4)
            sage: G.chromatic_quasisymmetric_function()
            (t^6+3*t^5+5*t^4+6*t^3+5*t^2+3*t+1)*M[1, 1, 1, 1]

        Not all chromatic quasisymmetric functions are symmetric::

            sage: G = Graph([[1,2], [1,5], [3,4], [3,5]])
            sage: G.chromatic_quasisymmetric_function().is_symmetric()
            False

        We check that at `t = 1`, we recover the usual chromatic
        symmetric function::

            sage: p = SymmetricFunctions(QQ).p()
            sage: G = graphs.CycleGraph(5)
            sage: XG = G.chromatic_quasisymmetric_function(t=1); XG
            120*M[1, 1, 1, 1, 1] + 30*M[1, 1, 1, 2] + 30*M[1, 1, 2, 1]
             + 30*M[1, 2, 1, 1] + 10*M[1, 2, 2] + 30*M[2, 1, 1, 1]
             + 10*M[2, 1, 2] + 10*M[2, 2, 1]
            sage: p(XG.to_symmetric_function())
            p[1, 1, 1, 1, 1] - 5*p[2, 1, 1, 1] + 5*p[2, 2, 1]
             + 5*p[3, 1, 1] - 5*p[3, 2] - 5*p[4, 1] + 4*p[5]

            sage: G = graphs.ClawGraph()
            sage: XG = G.chromatic_quasisymmetric_function(t=1); XG
            24*M[1, 1, 1, 1] + 6*M[1, 1, 2] + 6*M[1, 2, 1] + M[1, 3]
             + 6*M[2, 1, 1] + M[3, 1]
            sage: p(XG.to_symmetric_function())
            p[1, 1, 1, 1] - 3*p[2, 1, 1] + 3*p[3, 1] - p[4]

        REFERENCES:

        .. [SW12] John Shareshian and Michelle Wachs.
           *Chromatic quasisymmetric functions and Hessenberg varieties*.
           Configuration Spaces. CRM Series. Scuola Normale Superiore.
           (2012) pp. 433-460.
           http://www.math.miami.edu/~wachs/papers/chrom.pdf

        .. [BC15] Patrick Brosnan and Timothy Y. Chow.
           *Unit interval orders and the dot action on the cohomology
           of regular semisimple Hessenberg varieties*.
           (2015) :arxiv:`1511.00773v1`.
        """
        from sage.combinat.ncsf_qsym.qsym import QuasiSymmetricFunctions
        from sage.combinat.composition import Compositions
        from sage.combinat.set_partition_ordered import OrderedSetPartitions
        if t is None:
            t = ZZ['t'].gen()
        if R is None:
            R = t.parent()
        M = QuasiSymmetricFunctions(R).M()
        ret = M.zero()
        V = self.vertices()
        def asc(sigma):
            stat = 0
            for i, s in enumerate(sigma):
                for u in s:
                    stat += sum(1 for p in sigma[i+1:] for v in p
                                if v > u and self.has_edge(u, v))
            return stat
        for sigma in OrderedSetPartitions(V):
            if any(not self.is_independent_set(s) for s in sigma):
                continue
            ret += M.term(sigma.to_composition(), t**asc(sigma))
        return ret

    @doc_index("Leftovers")
    def matching(self, value_only=False, algorithm="Edmonds", use_edge_labels=True, solver=None, verbose=0):
        r"""
        Returns a maximum weighted matching of the graph
        represented by the list of its edges. For more information, see the
        `Wikipedia article on matchings
        <http://en.wikipedia.org/wiki/Matching_%28graph_theory%29>`_.

        Given a graph `G` such that each edge `e` has a weight `w_e`,
        a maximum matching is a subset `S` of the edges of `G` of
        maximum weight such that no two edges of `S` are incident
        with each other.

        As an optimization problem, it can be expressed as:

        .. math::

            \mbox{Maximize : }&\sum_{e\in G.edges()} w_e b_e\\
            \mbox{Such that : }&\forall v \in G, \sum_{(u,v)\in G.edges()} b_{(u,v)}\leq 1\\
            &\forall x\in G, b_x\mbox{ is a binary variable}

        INPUT:

        - ``value_only`` -- boolean (default: ``False``). When set to
          ``True``, only the cardinal (or the weight) of the matching is
          returned.

        - ``algorithm`` -- string (default: ``"Edmonds"``)

          - ``"Edmonds"`` selects Edmonds' algorithm as implemented in NetworkX

          - ``"LP"`` uses a Linear Program formulation of the matching problem

        - ``use_edge_labels`` -- boolean (default: ``False``)

          - When set to ``True``, computes a weighted matching where each edge
            is weighted by its label. (If an edge has no label, `1` is assumed.)

          - When set to ``False``, each edge has weight `1`.

        - ``solver`` -- (default: ``None``) Specify a Linear Program (LP)
          solver to be used. If set to ``None``, the default one is used. For
          more information on LP solvers and which default solver is used, see
          the method
          :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>`
          of the class
          :class:`MixedIntegerLinearProgram <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``). Sets the level of
          verbosity. Set to 0 by default, which means quiet.
          Only useful when ``algorithm == "LP"``.

        ALGORITHM:

        The problem is solved using Edmond's algorithm implemented in
        NetworkX, or using Linear Programming depending on the value of
        ``algorithm``.

        EXAMPLES:

        Maximum matching in a Pappus Graph::

           sage: g = graphs.PappusGraph()
           sage: g.matching(value_only=True)
           9.0

        Same test with the Linear Program formulation::

           sage: g = graphs.PappusGraph()
           sage: g.matching(algorithm="LP", value_only=True)
           9.0

        .. PLOT::

            g = graphs.PappusGraph()
            sphinx_plot(g.plot(edge_colors={"red":g.matching()}))

        TESTS:

        If ``algorithm`` is set to anything different from ``"Edmonds"`` or
        ``"LP"``, an exception is raised::

           sage: g = graphs.PappusGraph()
           sage: g.matching(algorithm="somethingdifferent")
           Traceback (most recent call last):
           ...
           ValueError: algorithm must be set to either "Edmonds" or "LP"
        """
        self._scream_if_not_simple(allow_loops=True)
        from sage.rings.real_mpfr import RR
        weight = lambda x: x if x in RR else 1

        if algorithm == "Edmonds":
            import networkx
            if use_edge_labels:
                g = networkx.Graph()
                for u, v, l in self.edges():
                    g.add_edge(u, v, attr_dict={"weight": weight(l)})
            else:
                g = self.networkx_graph(copy=False)
            d = networkx.max_weight_matching(g)
            if value_only:
                if use_edge_labels:
                    return sum(weight(self.edge_label(u, v))
                                for u, v in d.iteritems()) * 0.5
                else:
                    return Integer(len(d) // 2)
            else:
                return [(u, v, self.edge_label(u, v))
                        for u, v in d.iteritems() if u < v]

        elif algorithm == "LP":
            from sage.numerical.mip import MixedIntegerLinearProgram
            g = self
            # returns the weight of an edge considering it may not be
            # weighted ...
            p = MixedIntegerLinearProgram(maximization=True, solver=solver)
            b = p.new_variable(binary = True)
            p.set_objective(
                p.sum(weight(w) * b[min(u, v),max(u, v)]
                     for u, v, w in g.edges()))
            # for any vertex v, there is at most one edge incident to v in
            # the maximum matching
            for v in g.vertex_iterator():
                p.add_constraint(
                    p.sum(b[min(u, v),max(u, v)]
                         for u in g.neighbors(v)), max=1)
            if value_only:
                if use_edge_labels:
                    return p.solve(objective_only=True, log=verbose)
                else:
                    return Integer(round(p.solve(objective_only=True, log=verbose)))
            else:
                p.solve(log=verbose)
                b = p.get_values(b)
                return [(u, v, w) for u, v, w in g.edges()
                        if b[min(u, v),max(u, v)] == 1]

        else:
            raise ValueError('algorithm must be set to either "Edmonds" or "LP"')

    @doc_index("Algorithmically hard stuff")
    def has_homomorphism_to(self, H, core = False, solver = None, verbose = 0):
        r"""
        Checks whether there is a homomorphism between two graphs.

        A homomorphism from a graph `G` to a graph `H` is a function
        `\phi:V(G)\mapsto V(H)` such that for any edge `uv \in E(G)` the pair
        `\phi(u)\phi(v)` is an edge of `H`.

        Saying that a graph can be `k`-colored is equivalent to saying that it
        has a homomorphism to `K_k`, the complete graph on `k` elements.

        For more information, see the `Wikipedia article on graph homomorphisms
        <Graph_homomorphism>`_.

        INPUT:

        - ``H`` -- the graph to which ``self`` should be sent.

        - ``core`` (boolean) -- whether to minimize the size of the mapping's
          image (see note below). This is set to ``False`` by default.

        - ``solver`` -- (default: ``None``) Specify a Linear Program (LP)
          solver to be used. If set to ``None``, the default one is used. For
          more information on LP solvers and which default solver is used, see
          the method
          :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>`
          of the class
          :class:`MixedIntegerLinearProgram <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``). Sets the level of
          verbosity. Set to 0 by default, which means quiet.

        .. NOTE::

           One can compute the core of a graph (with respect to homomorphism)
           with this method ::

               sage: g = graphs.CycleGraph(10)
               sage: mapping = g.has_homomorphism_to(g, core = True)
               sage: print("The size of the core is {}".format(len(set(mapping.values()))))
               The size of the core is 2

        OUTPUT:

        This method returns ``False`` when the homomorphism does not exist, and
        returns the homomorphism otherwise as a dictionary associating a vertex
        of `H` to a vertex of `G`.

        EXAMPLE:

        Is Petersen's graph 3-colorable::

            sage: P = graphs.PetersenGraph()
            sage: P.has_homomorphism_to(graphs.CompleteGraph(3)) is not False
            True

        An odd cycle admits a homomorphism to a smaller odd cycle, but not to an
        even cycle::

            sage: g = graphs.CycleGraph(9)
            sage: g.has_homomorphism_to(graphs.CycleGraph(5)) is not False
            True
            sage: g.has_homomorphism_to(graphs.CycleGraph(7)) is not False
            True
            sage: g.has_homomorphism_to(graphs.CycleGraph(4)) is not False
            False
        """
        self._scream_if_not_simple()
        from sage.numerical.mip import MixedIntegerLinearProgram, MIPSolverException
        p = MixedIntegerLinearProgram(solver=solver, maximization = False)
        b = p.new_variable(binary = True)

        # Each vertex has an image
        for ug in self:
            p.add_constraint(p.sum(b[ug,uh] for uh in H) == 1)

        nonedges = H.complement().edges(labels = False)
        for ug,vg in self.edges(labels = False):
            # Two adjacent vertices cannot be mapped to the same element
            for uh in H:
                p.add_constraint(b[ug,uh] + b[vg,uh] <= 1)

            # Two adjacent vertices cannot be mapped to no adjacent vertices
            for uh,vh in nonedges:
                p.add_constraint(b[ug,uh] + b[vg,vh] <= 1)
                p.add_constraint(b[ug,vh] + b[vg,uh] <= 1)

        # Minimize the mapping's size
        if core:

            # the value of m is one if the corresponding vertex of h is used.
            m = p.new_variable(nonnegative=True)
            for uh in H:
                for ug in self:
                    p.add_constraint(b[ug,uh] <= m[uh])

            p.set_objective(p.sum(m[vh] for vh in H))

        try:
            p.solve(log = verbose)
            b = p.get_values(b)
            mapping = dict(x[0] for x in b.items() if x[1])
            return mapping

        except MIPSolverException:
            return False

    @doc_index("Leftovers")
    def fractional_chromatic_index(self, solver = None, verbose_constraints = 0, verbose = 0):
        r"""
        Computes the fractional chromatic index of ``self``

        The fractional chromatic index is a relaxed version of edge-coloring. An
        edge coloring of a graph being actually a covering of its edges into the
        smallest possible number of matchings, the fractional chromatic index of
        a graph `G` is the smallest real value `\chi_f(G)` such that there
        exists a list of matchings `M_1, ..., M_k` of `G` and coefficients
        `\alpha_1, ..., \alpha_k` with the property that each edge is covered by
        the matchings in the following relaxed way

        .. MATH::

            \forall e \in E(G), \sum_{e \in M_i} \alpha_i \geq 1

        For more information, see the `Wikipedia article on fractional coloring
        <http://en.wikipedia.org/wiki/Fractional_coloring>`_.

        ALGORITHM:

        The fractional chromatic index is computed through Linear Programming
        through its dual. The LP solved by sage is actually:

        .. MATH::

            \mbox{Maximize : }&\sum_{e\in E(G)} r_{e}\\
            \mbox{Such that : }&\\
            &\forall M\text{ matching }\subseteq G, \sum_{e\in M}r_{v}\leq 1\\

        INPUT:

        - ``solver`` -- (default: ``None``) Specify a Linear Program (LP)
          solver to be used. If set to ``None``, the default one is used. For
          more information on LP solvers and which default solver is used, see
          the method
          :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>`
          of the class
          :class:`MixedIntegerLinearProgram <sage.numerical.mip.MixedIntegerLinearProgram>`.

          .. NOTE::

              If you want exact results, i.e. a rational number, use
              ``solver="PPL"``. This may be slower, though.

        - ``verbose_constraints`` -- whether to display which constraints are
          being generated.

        - ``verbose`` -- level of verbosity required from the LP solver

        .. NOTE::

            This implementation can be improved by computing matchings through a
            LP formulation, and not using the Python implementation of Edmonds'
            algorithm (which requires to copy the graph, etc). It may be more
            efficient to write the matching problem as a LP, as we would then
            just have to update the weights on the edges between each call to
            ``solve`` (and so avoiding the generation of all the constraints).

        EXAMPLE:

        The fractional chromatic index of a `C_5` is `5/2`::

            sage: g = graphs.CycleGraph(5)
            sage: g.fractional_chromatic_index()
            2.5

        With PPL::

            sage: g.fractional_chromatic_index(solver="PPL")
            5/2
        """
        self._scream_if_not_simple()
        from sage.numerical.mip import MixedIntegerLinearProgram

        g = copy(self)
        p = MixedIntegerLinearProgram(solver=solver, constraint_generation = True)

        # One variable per edge
        r = p.new_variable(nonnegative=True)
        R = lambda x,y : r[x,y] if x<y else r[y,x]

        # We want to maximize the sum of weights on the edges
        p.set_objective( p.sum( R(u,v) for u,v in g.edges(labels = False)))

        # Each edge being by itself a matching, its weight can not be more than
        # 1

        for u,v in g.edges(labels = False):
            p.add_constraint( R(u,v), max = 1)

        obj = p.solve(log = verbose)

        while True:

            # Updating the value on the edges of g
            for u,v in g.edges(labels = False):
                g.set_edge_label(u,v,p.get_values(R(u,v)))

            # Computing a matching of maximum weight...

            matching = g.matching()

            # If the maximum matching has weight at most 1, we are done !
            if sum((x[2] for x in matching)) <= 1:
                break

            # Otherwise, we add a new constraint

            if verbose_constraints:
                print("Adding a constraint on matching : {}".format(matching))

            p.add_constraint( p.sum( R(u,v) for u,v,_ in matching), max = 1)

            # And solve again
            obj = p.solve(log = verbose)

        # Accomplished !
        return obj

    @doc_index("Leftovers")
    def maximum_average_degree(self, value_only=True, solver = None, verbose = 0):
        r"""
        Returns the Maximum Average Degree (MAD) of the current graph.

        The Maximum Average Degree (MAD) of a graph is defined as
        the average degree of its densest subgraph. More formally,
        ``Mad(G) = \max_{H\subseteq G} Ad(H)``, where `Ad(G)` denotes
        the average degree of `G`.

        This can be computed in polynomial time.

        INPUT:

        - ``value_only`` (boolean) -- ``True`` by default

          - If ``value_only=True``, only the numerical
            value of the `MAD` is returned.

          - Else, the subgraph of `G` realizing the `MAD`
            is returned.

        - ``solver`` -- (default: ``None``) Specify a Linear Program (LP)
          solver to be used. If set to ``None``, the default one is used. For
          more information on LP solvers and which default solver is used, see
          the method
          :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>`
          of the class
          :class:`MixedIntegerLinearProgram <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``). Sets the level of
          verbosity. Set to 0 by default, which means quiet.

        EXAMPLES:

        In any graph, the `Mad` is always larger than the average
        degree::

            sage: g = graphs.RandomGNP(20,.3)
            sage: mad_g = g.maximum_average_degree()
            sage: g.average_degree() <= mad_g
            True

        Unlike the average degree, the `Mad` of the disjoint
        union of two graphs is the maximum of the `Mad` of each
        graphs::

            sage: h = graphs.RandomGNP(20,.3)
            sage: mad_h = h.maximum_average_degree()
            sage: (g+h).maximum_average_degree() == max(mad_g, mad_h)
            True

        The subgraph of a regular graph realizing the maximum
        average degree is always the whole graph ::

            sage: g = graphs.CompleteGraph(5)
            sage: mad_g = g.maximum_average_degree(value_only=False)
            sage: g.is_isomorphic(mad_g)
            True

        This also works for complete bipartite graphs ::

            sage: g = graphs.CompleteBipartiteGraph(3,4)
            sage: mad_g = g.maximum_average_degree(value_only=False)
            sage: g.is_isomorphic(mad_g)
            True
        """
        self._scream_if_not_simple()
        g = self
        from sage.numerical.mip import MixedIntegerLinearProgram

        p = MixedIntegerLinearProgram(maximization=True, solver = solver)

        d = p.new_variable(nonnegative=True)
        one = p.new_variable(nonnegative=True)

        # Reorders u and v so that uv and vu are not considered
        # to be different edges
        reorder = lambda u,v : (min(u,v),max(u,v))

        for u,v in g.edge_iterator(labels=False):
            p.add_constraint( one[ reorder(u,v) ] - 2*d[u] , max = 0 )
            p.add_constraint( one[ reorder(u,v) ] - 2*d[v] , max = 0 )

        p.add_constraint( p.sum(d[v] for v in g), max = 1)

        p.set_objective( p.sum( one[reorder(u,v)] for u,v in g.edge_iterator(labels=False)) )

        obj = p.solve(log = verbose)

        # Paying attention to numerical error :
        # The zero values could be something like 0.000000000001
        # so I can not write l > 0
        # And the non-zero, though they should be equal to
        # 1/(order of the optimal subgraph) may be a bit lower

        # setting the minimum to 1/(10 * size of the whole graph )
        # should be safe :-)
        m = 1/(10 *Integer(g.order()))
        g_mad = g.subgraph([v for v,l in p.get_values(d).iteritems() if l>m ])

        if value_only:
            return g_mad.average_degree()
        else:
            return g_mad

    @doc_index("Algorithmically hard stuff")
    def independent_set_of_representatives(self, family, solver=None, verbose=0):
        r"""
        Returns an independent set of representatives.

        Given a graph `G` and and a family `F=\{F_i:i\in [1,...,k]\}` of
        subsets of ``g.vertices()``, an Independent Set of Representatives
        (ISR) is an assignation of a vertex `v_i\in F_i` to each set `F_i`
        such that `v_i != v_j` if `i<j` (they are representatives) and the
        set `\cup_{i}v_i` is an independent set in `G`.

        It generalizes, for example, graph coloring and graph list coloring.

        (See [AhaBerZiv07]_ for more information.)

        INPUT:

        - ``family`` -- A list of lists defining the family `F`
          (actually, a Family of subsets of ``G.vertices()``).

        - ``solver`` -- (default: ``None``) Specify a Linear Program (LP)
          solver to be used. If set to ``None``, the default one is used. For
          more information on LP solvers and which default solver is used, see
          the method
          :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>`
          of the class
          :class:`MixedIntegerLinearProgram <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``). Sets the level of
          verbosity. Set to 0 by default, which means quiet.

        OUTPUT:

        - A list whose `i^{\mbox{th}}` element is the representative of the
          `i^{\mbox{th}}` element of the ``family`` list. If there is no ISR,
          ``None`` is returned.

        EXAMPLES:

        For a bipartite graph missing one edge, the solution is as expected::

           sage: g = graphs.CompleteBipartiteGraph(3,3)
           sage: g.delete_edge(1,4)
           sage: g.independent_set_of_representatives([[0,1,2],[3,4,5]])
           [1, 4]

        The Petersen Graph is 3-colorable, which can be expressed as an
        independent set of representatives problem : take 3 disjoint copies
        of the Petersen Graph, each one representing one color. Then take
        as a partition of the set of vertices the family defined by the three
        copies of each vertex. The ISR of such a family
        defines a 3-coloring::

            sage: g = 3 * graphs.PetersenGraph()
            sage: n = g.order()/3
            sage: f = [[i,i+n,i+2*n] for i in xrange(n)]
            sage: isr = g.independent_set_of_representatives(f)
            sage: c = [floor(i/n) for i in isr]
            sage: color_classes = [[],[],[]]
            sage: for v,i in enumerate(c):
            ...     color_classes[i].append(v)
            sage: for classs in color_classes:
            ...     g.subgraph(classs).size() == 0
            True
            True
            True

        REFERENCE:

        .. [AhaBerZiv07] \R. Aharoni and E. Berger and R. Ziv
          Independent systems of representatives in weighted graphs
          Combinatorica vol 27, num 3, p253--267
          2007

        """

        from sage.numerical.mip import MixedIntegerLinearProgram
        p=MixedIntegerLinearProgram(solver=solver)

        # Boolean variable indicating whether the vertex
        # is the representative of some set
        vertex_taken=p.new_variable(binary=True)

        # Boolean variable in two dimension whose first
        # element is a vertex and whose second element
        # is one of the sets given as arguments.
        # When true, indicated that the vertex is the representant
        # of the corresponding set

        classss=p.new_variable(binary = True)

        # Associates to the vertices the classes
        # to which they belong

        lists=dict([(v,[]) for v in self.vertex_iterator()])
        for i,f in enumerate(family):
            [lists[v].append(i) for v in f]

            # a classss has exactly one representant
            p.add_constraint(p.sum(classss[v,i] for v in f), max=1, min=1)

        # A vertex represents at most one classss (vertex_taken is binary), and
        # vertex_taken[v]==1 if v is the representative of some classss

        [p.add_constraint(p.sum(classss[v,i] for i in lists[v]) - vertex_taken[v], max=0) for v in self.vertex_iterator()]

        # Two adjacent vertices can not both be representants of a set

        for (u,v) in self.edges(labels=None):
            p.add_constraint(vertex_taken[u]+vertex_taken[v],max=1)

        p.set_objective(None)

        try:
            p.solve(log=verbose)
        except Exception:
            return None

        classss=p.get_values(classss)

        repr=[]
        for i,f in enumerate(family):
            for v in f:
                if classss[v,i]==1:
                    repr.append(v)
                    break

        return repr

    @doc_index("Algorithmically hard stuff")
    def minor(self, H, solver=None, verbose=0):
        r"""
        Returns the vertices of a minor isomorphic to `H` in the current graph.

        We say that a graph `G` has a `H`-minor (or that it has
        a graph isomorphic to `H` as a minor), if for all `h\in H`,
        there exist disjoint sets `S_h \subseteq V(G)` such that
        once the vertices of each `S_h` have been merged to create
        a new graph `G'`, this new graph contains `H` as a subgraph.

        For more information, see the
        `Wikipedia article on graph minor <http://en.wikipedia.org/wiki/Minor_%28graph_theory%29>`_.

        INPUT:

        - ``H`` -- The minor to find for in the current graph.

        - ``solver`` -- (default: ``None``) Specify a Linear Program (LP)
          solver to be used. If set to ``None``, the default one is used. For
          more information on LP solvers and which default solver is used, see
          the method
          :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>`
          of the class
          :class:`MixedIntegerLinearProgram <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``). Sets the level of
          verbosity. Set to 0 by default, which means quiet.

        OUTPUT:

        A dictionary associating to each vertex of `H` the set of vertices
        in the current graph representing it.

        ALGORITHM:

        Mixed Integer Linear Programming

        COMPLEXITY:

        Theoretically, when `H` is fixed, testing for the existence of
        a `H`-minor is polynomial. The known algorithms are highly
        exponential in `H`, though.

        .. NOTE::

            This function can be expected to be *very* slow, especially
            where the minor does not exist.

        EXAMPLES:

        Trying to find a minor isomorphic to `K_4` in
        the `4\times 4` grid::

            sage: g = graphs.GridGraph([4,4])
            sage: h = graphs.CompleteGraph(4)
            sage: L = g.minor(h)
            sage: gg = g.subgraph(flatten(L.values(), max_level = 1))
            sage: _ = [gg.merge_vertices(l) for l in L.values() if len(l)>1]
            sage: gg.is_isomorphic(h)
            True

        We can also try to prove this way that the Petersen graph
        is not planar, as it has a `K_5` minor::

            sage: g = graphs.PetersenGraph()
            sage: K5_minor = g.minor(graphs.CompleteGraph(5))                    # long time

        And even a `K_{3,3}` minor::

            sage: K33_minor = g.minor(graphs.CompleteBipartiteGraph(3,3))        # long time

        (It is much faster to use the linear-time test of
        planarity in this situation, though.)

        As there is no cycle in a tree, looking for a `K_3` minor is useless.
        This function will raise an exception in this case::

            sage: g = graphs.RandomGNP(20,.5)
            sage: g = g.subgraph(edges = g.min_spanning_tree())
            sage: g.is_tree()
            True
            sage: L = g.minor(graphs.CompleteGraph(3))
            Traceback (most recent call last):
            ...
            ValueError: This graph has no minor isomorphic to H !
        """
        self._scream_if_not_simple()
        H._scream_if_not_simple()
        from sage.numerical.mip import MixedIntegerLinearProgram, MIPSolverException
        p = MixedIntegerLinearProgram(solver=solver)

        # sorts an edge
        S = lambda x_y: x_y if x_y[0] < x_y[1] else (x_y[1], x_y[0])

        # rs = Representative set of a vertex
        # for h in H, v in G is such that rs[h,v] == 1 if and only if v
        # is a representant of h in self
        rs = p.new_variable(binary = True)

        for v in self:
            p.add_constraint(p.sum(rs[h,v] for h in H), max = 1)

        # We ensure that the set of representatives of a
        # vertex h contains a tree, and thus is connected

        # edges represents the edges of the tree
        edges = p.new_variable(binary = True)

        # there can be a edge for h between two vertices
        # only if those vertices represent h
        for u,v in self.edges(labels=None):
            for h in H:
                p.add_constraint(edges[h,S((u,v))] - rs[h,u], max = 0 )
                p.add_constraint(edges[h,S((u,v))] - rs[h,v], max = 0 )

        # The number of edges of the tree in h is exactly the cardinal
        # of its representative set minus 1

        for h in H:
            p.add_constraint(p.sum(edges[h,S(e)] for e in self.edges(labels=None))-p.sum(rs[h,v] for v in self), min=-1, max=-1)

        # a tree  has no cycle
        epsilon = 1/(5*Integer(self.order()))
        r_edges = p.new_variable(nonnegative=True)

        for h in H:
            for u,v in self.edges(labels=None):
                p.add_constraint(r_edges[h,(u,v)] + r_edges[h,(v,u)] - edges[h,S((u,v))], min = 0)

            for v in self:
                p.add_constraint(p.sum(r_edges[h,(u,v)] for u in self.neighbors(v)), max = 1-epsilon)

        # Once the representative sets are described, we must ensure
        # there are arcs corresponding to those of H between them
        h_edges = p.new_variable(nonnegative=True)

        for h1, h2 in H.edges(labels=None):

            for v1, v2 in self.edges(labels=None):

                p.add_constraint(h_edges[(h1,h2),S((v1,v2))] - rs[h2,v2], max = 0)
                p.add_constraint(h_edges[(h1,h2),S((v1,v2))] - rs[h1,v1], max = 0)

                p.add_constraint(h_edges[(h2,h1),S((v1,v2))] - rs[h1,v2], max = 0)
                p.add_constraint(h_edges[(h2,h1),S((v1,v2))] - rs[h2,v1], max = 0)

            p.add_constraint(p.sum(h_edges[(h1,h2),S(e)] + h_edges[(h2,h1),S(e)] for e in self.edges(labels=None) ), min = 1)

        p.set_objective(None)

        try:
            p.solve(log=verbose)
        except MIPSolverException:
            raise ValueError("This graph has no minor isomorphic to H !")

        rs = p.get_values(rs)

        rs_dict = {}
        for h in H:
            rs_dict[h] = [v for v in self if rs[h,v]==1]

        return rs_dict

    ### Convexity

    @doc_index("Algorithmically hard stuff")
    def convexity_properties(self):
        r"""
        Returns a ``ConvexityProperties`` object corresponding to ``self``.

        This object contains the methods related to convexity in graphs (convex
        hull, hull number) and caches useful information so that it becomes
        comparatively cheaper to compute the convex hull of many different sets
        of the same graph.

        .. SEEALSO::

            In order to know what can be done through this object, please refer
            to module :mod:`sage.graphs.convexity_properties`.

        .. NOTE::

            If you want to compute many convex hulls, keep this object in memory
            ! When it is created, it builds a table of useful information to
            compute convex hulls. As a result ::

                sage: g = graphs.PetersenGraph()
                sage: g.convexity_properties().hull([1, 3])
                [1, 2, 3]
                sage: g.convexity_properties().hull([3, 7])
                [2, 3, 7]

            Is a terrible waste of computations, while ::

                sage: g = graphs.PetersenGraph()
                sage: CP = g.convexity_properties()
                sage: CP.hull([1, 3])
                [1, 2, 3]
                sage: CP.hull([3, 7])
                [2, 3, 7]

            Makes perfect sense.
        """
        from sage.graphs.convexity_properties import ConvexityProperties
        return ConvexityProperties(self)

    # Centrality
    @doc_index("Distances")
    def centrality_degree(self, v=None):
        r"""
        Returns the degree centrality of a vertex.

        The degree centrality of a vertex `v` is its degree, divided by
        `|V(G)|-1`. For more information, see the :wikipedia:`Centrality`.

        INPUT:

        - ``v`` - a vertex. Set to ``None`` (default) to get a dictionary
          associating each vertex with its centrality degree.

        .. SEEALSO::

            - :meth:`~sage.graphs.generic_graph.GenericGraph.centrality_closeness`
            - :meth:`~sage.graphs.generic_graph.GenericGraph.centrality_betweenness`

        EXAMPLES::

            sage: (graphs.ChvatalGraph()).centrality_degree()
            {0: 4/11, 1: 4/11, 2: 4/11, 3: 4/11,  4: 4/11,  5: 4/11,
             6: 4/11, 7: 4/11, 8: 4/11, 9: 4/11, 10: 4/11, 11: 4/11}
            sage: D = graphs.DiamondGraph()
            sage: D.centrality_degree()
            {0: 2/3, 1: 1, 2: 1, 3: 2/3}
            sage: D.centrality_degree(v=1)
            1

        TESTS::

            sage: Graph(1).centrality_degree()
            Traceback (most recent call last):
            ...
            ValueError: The centrality degree is not defined on graphs with only one vertex
        """
        from sage.rings.integer import Integer
        n_minus_one = Integer(self.order()-1)
        if n_minus_one == 0:
            raise ValueError("The centrality degree is not defined "
                             "on graphs with only one vertex")
        if v is None:
            return {v:self.degree(v)/n_minus_one for v in self}
        else:
            return self.degree(v)/n_minus_one

    ### Constructors

    @doc_index("Basic methods")
    def to_directed(self, implementation='c_graph', data_structure=None,
                    sparse=None):
        """
        Returns a directed version of the graph. A single edge becomes two
        edges, one in each direction.

        INPUT:

         - ``data_structure`` -- one of ``"sparse"``, ``"static_sparse"``, or
           ``"dense"``. See the documentation of :class:`Graph` or
           :class:`DiGraph`.

         - ``sparse`` (boolean) -- ``sparse=True`` is an alias for
           ``data_structure="sparse"``, and ``sparse=False`` is an alias for
           ``data_structure="dense"``.

        EXAMPLES::

            sage: graphs.PetersenGraph().to_directed()
            Petersen graph: Digraph on 10 vertices

        TESTS:

        Immutable graphs yield immutable graphs::

            sage: Graph([[1, 2]], immutable=True).to_directed()._backend
            <type 'sage.graphs.base.static_sparse_backend.StaticSparseBackend'>

        :trac:`17005`::

            sage: Graph([[1,2]], immutable=True).to_directed()
            Digraph on 2 vertices
        """
        if sparse is not None:
            if data_structure is not None:
                raise ValueError("The 'sparse' argument is an alias for "
                                 "'data_structure'. Please do not define both.")
            data_structure = "sparse" if sparse else "dense"

        if data_structure is None:
            from sage.graphs.base.dense_graph import DenseGraphBackend
            from sage.graphs.base.sparse_graph import SparseGraphBackend
            if isinstance(self._backend, DenseGraphBackend):
                data_structure = "dense"
            elif isinstance(self._backend, SparseGraphBackend):
                data_structure = "sparse"
            else:
                data_structure = "static_sparse"
        from sage.graphs.all import DiGraph
        D = DiGraph(name           = self.name(),
                    pos            = self._pos,
                    multiedges     = self.allows_multiple_edges(),
                    loops          = self.allows_loops(),
                    implementation = implementation,
                    data_structure = (data_structure if data_structure!="static_sparse"
                                      else "sparse")) # we need a mutable copy

        D.add_vertices(self.vertex_iterator())
        for u,v,l in self.edge_iterator():
            D.add_edge(u,v,l)
            D.add_edge(v,u,l)
        if hasattr(self, '_embedding'):
            D._embedding = copy(self._embedding)
        D._weighted = self._weighted

        if data_structure == "static_sparse":
            D = D.copy(data_structure=data_structure)

        return D

    @doc_index("Basic methods")
    def to_undirected(self):
        """
        Since the graph is already undirected, simply returns a copy of
        itself.

        EXAMPLES::

            sage: graphs.PetersenGraph().to_undirected()
            Petersen graph: Graph on 10 vertices
        """
        return self.copy()

    @doc_index("Basic methods")
    def join(self, other, verbose_relabel=None, labels="pairs", immutable=None):
        """
        Returns the join of ``self`` and ``other``.

        INPUT:

        - ``verbose_relabel`` - deprecated.

        - ``labels`` - (defaults to 'pairs') If set to 'pairs', each
          element ``v`` in the first graph will be named ``(0,v)`` and
          each element ``u`` in ``other`` will be named ``(1,u)`` in
          the result. If set to 'integers', the elements of the result
          will be relabeled with consecutive integers.

        - ``immutable`` (boolean) -- whether to create a mutable/immutable
          join. ``immutable=None`` (default) means that the graphs and their
          join will behave the same way.

        .. SEEALSO::

            * :meth:`~sage.graphs.generic_graph.GenericGraph.union`

            * :meth:`~sage.graphs.generic_graph.GenericGraph.disjoint_union`

        EXAMPLES::

            sage: G = graphs.CycleGraph(3)
            sage: H = Graph(2)
            sage: J = G.join(H); J
            Cycle graph join : Graph on 5 vertices
            sage: J.vertices()
            [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1)]
            sage: J = G.join(H, labels='integers'); J
            Cycle graph join : Graph on 5 vertices
            sage: J.vertices()
            [0, 1, 2, 3, 4]
            sage: J.edges()
            [(0, 1, None), (0, 2, None), (0, 3, None), (0, 4, None), (1, 2, None), (1, 3, None), (1, 4, None), (2, 3, None), (2, 4, None)]

        ::

            sage: G = Graph(3)
            sage: G.name("Graph on 3 vertices")
            sage: H = Graph(2)
            sage: H.name("Graph on 2 vertices")
            sage: J = G.join(H); J
            Graph on 3 vertices join Graph on 2 vertices: Graph on 5 vertices
            sage: J.vertices()
            [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1)]
            sage: J = G.join(H, labels='integers'); J
            Graph on 3 vertices join Graph on 2 vertices: Graph on 5 vertices
            sage: J.edges()
            [(0, 3, None), (0, 4, None), (1, 3, None), (1, 4, None), (2, 3, None), (2, 4, None)]
        """
        if verbose_relabel is not None:
            deprecation(17053, "Instead of verbose_relabel=True/False use labels='pairs'/'integers'.")
            if verbose_relabel is True:
                labels="pairs"
            if verbose_relabel is False:
                labels="integers"

        G = self.disjoint_union(other, labels=labels, immutable=False)
        if labels=="integers":
            G.add_edges((u,v) for u in range(self.order())
                        for v in range(self.order(), self.order()+other.order()))
        else:
            G.add_edges(((0,u), (1,v)) for u in self.vertices()
                        for v in other.vertices())

        G.name('%s join %s'%(self.name(), other.name()))

        if immutable is None:
            immutable = self.is_immutable() and other.is_immutable()
        if immutable:
            G = G.copy(immutable=True)

        return G

    @doc_index("Leftovers")
    def seidel_adjacency_matrix(self, vertices=None):
        r"""
        Returns the Seidel adjacency matrix of ``self``.

        Returns `J-I-2A`, for `A` the (ordinary) :meth:`adjacency
        matrix
        <sage.graphs.generic_graph.GenericGraph.adjacency_matrix>` of
        ``self``, `I` the identity matrix, and `J` the all-1 matrix.
        It is closely related to :meth:`twograph`.

        The matrix returned is over the integers. If a different ring is
        desired, use either :meth:`sage.matrix.matrix0.Matrix.change_ring`
        method or :class:`matrix <sage.matrix.constructor.MatrixFactory>` function.

        INPUT:

        - ``vertices`` (list) -- the ordering of the vertices defining
          how they should appear in the matrix. By default, the
          ordering given by
          :meth:`~sage.graphs.generic_graph.GenericGraph.vertices` is
          used.

        EXAMPLES::

            sage: G = graphs.CycleGraph(5)
            sage: G = G.disjoint_union(graphs.CompleteGraph(1))
            sage: G.seidel_adjacency_matrix().minpoly()
            x^2 - 5
        """

        return -self.adjacency_matrix(sparse=False, vertices=vertices)+ \
                  self.complement().adjacency_matrix(sparse=False, \
                                            vertices=vertices)

    @doc_index("Leftovers")
    def seidel_switching(self, s, inplace=True):
        r"""
        Returns the Seidel switching of ``self`` w.r.t. subset of vertices ``s``.

        Returns the graph obtained by Seidel switching of ``self``
        with respect to the subset of vertices ``s``. This is the graph
        given by Seidel adjacency matrix `DSD`, for `S` the Seidel
        adjacency matrix of ``self``, and `D` the diagonal matrix with -1s
        at positions corresponding to ``s``, and 1s elsewhere.

        INPUT:

         - ``s`` -- a list of vertices of ``self``

        - ``inplace`` (boolean) -- whether to do the modification inplace, or to
          return a copy of the graph after switching.

        EXAMPLES::

            sage: G = graphs.CycleGraph(5)
            sage: G = G.disjoint_union(graphs.CompleteGraph(1))
            sage: G.seidel_switching([(0,1),(1,0),(0,0)])
            sage: G.seidel_adjacency_matrix().minpoly()
            x^2 - 5
            sage: G.is_connected()
            True

        TESTS::

            sage: H = G.seidel_switching([1,4,5],inplace=False)
            sage: G.seidel_switching([1,4,5])
            sage: G == H
            True
        """
        from itertools import product
        G = self if inplace else copy(self)
        boundary = self.edge_boundary(s)
        G.add_edges(product(s, set(self).difference(s)))
        G.delete_edges(boundary)
        if not inplace:
            return G

    @doc_index("Leftovers")
    def twograph(self):
        r"""
        Returns the two-graph of ``self``

        Returns the :class:`two-graph <sage.combinat.designs.twographs.TwoGraph>`
        with the triples
        `T=\{t \in \binom {V}{3} : \left| \binom {t}{2} \cap E \right| \text{odd} \}`
        where `V` and `E` are vertices and edges of ``self``, respectively.

        EXAMPLES::

            sage: p=graphs.PetersenGraph()
            sage: p.twograph()
            Incidence structure with 10 points and 60 blocks
            sage: p=graphs.chang_graphs()
            sage: T8 = graphs.CompleteGraph(8).line_graph()
            sage: C = T8.seidel_switching([(0,1,None),(2,3,None),(4,5,None),(6,7,None)],inplace=False)
            sage: T8.twograph()==C.twograph()
            True
            sage: T8.is_isomorphic(C)
            False

        TESTS::

            sage: from sage.combinat.designs.twographs import TwoGraph
            sage: p=graphs.PetersenGraph().twograph()
            sage: TwoGraph(p, check=True)
            Incidence structure with 10 points and 60 blocks

        .. SEEALSO::

            - :meth:`~sage.combinat.designs.twographs.TwoGraph.descendant`
              -- computes the descendant graph of the two-graph of self at a vertex

            - :func:`~sage.combinat.designs.twographs.twograph_descendant`
              -- ditto, but much faster.
        """
        from sage.combinat.designs.twographs import TwoGraph
        G = self.relabel(inplace=False)
        T = []

        # Triangles
        for x,y,z in G.subgraph_search_iterator(Graph({1:[2,3],2:[3]})):
            if x < y and y < z:
                T.append([x,y,z])

        # Triples with just one edge
        for x,y,z in G.subgraph_search_iterator(Graph({1:[2],3:[]}),induced=True):
            if x < y:
                T.append([x,y,z])

        T = TwoGraph(T)
        T.relabel({i:v for i,v in enumerate(self.vertices())})

        return T

    ### Visualization

    @doc_index("Basic methods")
    def write_to_eps(self, filename, **options):
        r"""
        Writes a plot of the graph to ``filename`` in ``eps`` format.

        INPUT:

         - ``filename`` -- a string
         - ``**options`` -- same layout options as :meth:`.layout`

        EXAMPLES::

            sage: P = graphs.PetersenGraph()
            sage: P.write_to_eps(tmp_filename(ext='.eps'))

        It is relatively simple to include this file in a LaTeX
        document.  ``\usepackage{graphics}`` must appear in the
        preamble, and ``\includegraphics{filename}`` will include
        the file. To compile the document to ``pdf`` with ``pdflatex`` or ``xelatex``
        the file needs first to be converted to ``pdf``, for example
        with ``ps2pdf filename.eps filename.pdf``.
        """
        from sage.graphs.print_graphs import print_graph_eps
        pos = self.layout(**options)
        [xmin, xmax, ymin, ymax] = self._layout_bounding_box(pos)
        for v in pos:
            pos[v] = (1.8*(pos[v][0] - xmin)/(xmax - xmin) - 0.9, 1.8*(pos[v][1] - ymin)/(ymax - ymin) - 0.9)
        if filename[-4:] != '.eps':
            filename += '.eps'
        f = open(filename, 'w')
        f.write( print_graph_eps(self.vertices(), self.edge_iterator(), pos) )
        f.close()

    @doc_index("Algorithmically hard stuff")
    def topological_minor(self, H, vertices = False, paths = False, solver=None, verbose=0):
        r"""
        Returns a topological `H`-minor from ``self`` if one exists.

        We say that a graph `G` has a topological `H`-minor (or that
        it has a graph isomorphic to `H` as a topological minor), if
        `G` contains a subdivision of a graph isomorphic to `H` (i.e.
        obtained from `H` through arbitrary subdivision of its edges)
        as a subgraph.

        For more information, see the :wikipedia:`Minor_(graph_theory)`.

        INPUT:

        - ``H`` -- The topological minor to find in the current graph.

        - ``solver`` -- (default: ``None``) Specify a Linear Program (LP)
          solver to be used. If set to ``None``, the default one is used. For
          more information on LP solvers and which default solver is used, see
          the method
          :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>`
          of the class
          :class:`MixedIntegerLinearProgram <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``). Sets the level of
          verbosity. Set to 0 by default, which means quiet.

        OUTPUT:

        The topological `H`-minor found is returned as a subgraph `M`
        of ``self``, such that the vertex `v` of `M` that represents a
        vertex `h\in H` has ``h`` as a label (see
        :meth:`get_vertex <sage.graphs.generic_graph.GenericGraph.get_vertex>`
        and
        :meth:`set_vertex <sage.graphs.generic_graph.GenericGraph.set_vertex>`),
        and such that every edge of `M` has as a label the edge of `H`
        it (partially) represents.

        If no topological minor is found, this method returns
        ``False``.

        ALGORITHM:

        Mixed Integer Linear Programming.

        COMPLEXITY:

        Theoretically, when `H` is fixed, testing for the existence of
        a topological `H`-minor is polynomial. The known algorithms
        are highly exponential in `H`, though.

        .. NOTE::

            This function can be expected to be *very* slow, especially where
            the topological minor does not exist.

            (CPLEX seems to be *much* more efficient than GLPK on this kind of
            problem)

        EXAMPLES:

        Petersen's graph has a topological `K_4`-minor::

            sage: g = graphs.PetersenGraph()
            sage: g.topological_minor(graphs.CompleteGraph(4))
            Subgraph of (Petersen graph): Graph on ...

        And a topological `K_{3,3}`-minor::

            sage: g.topological_minor(graphs.CompleteBipartiteGraph(3,3))
            Subgraph of (Petersen graph): Graph on ...

        And of course, a tree has no topological `C_3`-minor::

            sage: g = graphs.RandomGNP(15,.3)
            sage: g = g.subgraph(edges = g.min_spanning_tree())
            sage: g.topological_minor(graphs.CycleGraph(3))
            False
        """
        self._scream_if_not_simple()
        H._scream_if_not_simple()
        # Useful alias ...
        G = self

        from sage.numerical.mip import MixedIntegerLinearProgram, MIPSolverException
        p = MixedIntegerLinearProgram(solver=solver)

        # This is an existence problem
        p.set_objective(None)

        #######################
        # Vertex representant #
        #######################
        #
        # v_repr[h,g] = 1 if vertex h from H is represented by vertex
        # g from G, 0 otherwise

        v_repr = p.new_variable(binary = True)

        # Exactly one representant per vertex of H
        for h in H:
            p.add_constraint( p.sum( v_repr[h,g] for g in G), min = 1, max = 1)

        # A vertex of G can only represent one vertex of H
        for g in G:
            p.add_constraint( p.sum( v_repr[h,g] for h in H), max = 1)

        ###################
        # Is representent #
        ###################
        #
        # is_repr[v] = 1 if v represents some vertex of H

        is_repr = p.new_variable(binary = True)

        for g in G:
            for h in H:
                p.add_constraint( v_repr[h,g] - is_repr[g], max = 0)

        ###################################
        # paths between the representents #
        ###################################
        #
        # For any edge (h1,h2) in H, we have a corresponding path in G
        # between the representants of h1 and h2. Which means there is
        # a flow of intensity 1 from one to the other.
        # We are then writing a flow problem for each edge of H.
        #
        # The variable flow[(h1,h2),(g1,g2)] indicates the amount of
        # flow on the edge (g1,g2) representing the edge (h1,h2).

        flow = p.new_variable(binary = True)

        # This lambda function returns the balance of flow
        # corresponding to commodity C at vertex v v

        flow_in = lambda C, v : p.sum( flow[C,(v,u)] for u in G.neighbors(v) )
        flow_out = lambda C, v : p.sum( flow[C,(u,v)] for u in G.neighbors(v) )

        flow_balance = lambda C, v : flow_in(C,v) - flow_out(C,v)

        for h1,h2 in H.edges(labels = False):

            for v in G:

                # The flow balance depends on whether the vertex v is
                # a representant of h1 or h2 in G, or a reprensentant
                # of none

                p.add_constraint( flow_balance((h1,h2),v) == v_repr[h1,v] - v_repr[h2,v] )

        #############################
        # Internal vertex of a path #
        #############################
        #
        # is_internal[C][g] = 1 if a vertex v from G is located on the
        # path representing the edge (=commodity) C

        is_internal = p.new_variable(binary = True)

        # When is a vertex internal for a commodity ?
        for C in H.edges(labels = False):
            for g in G:
                p.add_constraint( flow_in(C,g) + flow_out(C,g) - is_internal[C,g], max = 1)

        ############################
        # Two paths do not cross ! #
        ############################

        # A vertex can only be internal for one commodity, and zero if
        # the vertex is a representent

        for g in G:
            p.add_constraint( p.sum( is_internal[C,g] for C in H.edges(labels = False))
                              + is_repr[g], max = 1 )

        # (The following inequalities are not necessary, but they seem
        # to be of help (the solvers find the answer quicker when they
        # are added)

        # The flow on one edge can go in only one direction. Besides,
        # it can belong to at most one commodity and has a maximum
        # intensity of 1.

        for g1,g2 in G.edges(labels = None):

            p.add_constraint(   p.sum( flow[C,(g1,g2)] for C in H.edges(labels = False) )
                              + p.sum( flow[C,(g2,g1)] for C in H.edges(labels = False) ),
                                max = 1)


        # Now we can solve the problem itself !

        try:
            p.solve(log = verbose)

        except MIPSolverException:
            return False


        minor = G.subgraph(immutable=False)

        is_repr = p.get_values(is_repr)
        v_repr = p.get_values(v_repr)
        flow = p.get_values(flow)

        for u,v in minor.edges(labels = False):
            used = False
            for C in H.edges(labels = False):

                if flow[C,(u,v)] + flow[C,(v,u)] > .5:
                    used = True
                    minor.set_edge_label(u,v,C)
                    break
            if not used:
                minor.delete_edge(u,v)

        minor.delete_vertices( [v for v in minor
                                if minor.degree(v) == 0 ] )

        for g in minor:
            if is_repr[g] > .5:
                for h in H:
                    if v_repr[h,v] > .5:
                        minor.set_vertex(g,h)
                        break

        return minor

    ### Cliques

    @doc_index("Clique-related methods")
    def cliques_maximal(self, algorithm = "native"):
        """
        Returns the list of all maximal cliques, with each clique represented
        by a list of vertices. A clique is an induced complete subgraph, and a
        maximal clique is one not contained in a larger one.

        INPUT:

        - ``algorithm`` -- can be set to ``"native"`` (default) to use Sage's
          own implementation, or to ``"NetworkX"`` to use NetworkX'
          implementation of the Bron and Kerbosch Algorithm [BroKer1973]_.


        .. NOTE::

            This method sorts its output before returning it. If you prefer to
            save the extra time, you can call
            :class:`sage.graphs.independent_sets.IndependentSets` directly.

        .. NOTE::

            Sage's implementation of the enumeration of *maximal* independent
            sets is not much faster than NetworkX' (expect a 2x speedup), which
            is surprising as it is written in Cython. This being said, the
            algorithm from NetworkX appears to be sligthly different from this
            one, and that would be a good thing to explore if one wants to
            improve the implementation.

        ALGORITHM:

        This function is based on NetworkX's implementation of the Bron and
        Kerbosch Algorithm [BroKer1973]_.

        REFERENCE:

        .. [BroKer1973] Coen Bron and Joep Kerbosch. (1973). Algorithm 457:
          Finding All Cliques of an Undirected Graph. Commun. ACM. v
          16. n 9.  pages 575-577. ACM Press. [Online] Available:
          http://www.ram.org/computing/rambin/rambin.html

        EXAMPLES::

            sage: graphs.ChvatalGraph().cliques_maximal()
            [[0, 1], [0, 4], [0, 6], [0, 9], [1, 2], [1, 5], [1, 7], [2, 3],
             [2, 6], [2, 8], [3, 4], [3, 7], [3, 9], [4, 5], [4, 8], [5, 10],
             [5, 11], [6, 10], [6, 11], [7, 8], [7, 11], [8, 10], [9, 10], [9, 11]]
            sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage: G.show(figsize=[2,2])
            sage: G.cliques_maximal()
            [[0, 1, 2], [0, 1, 3]]
            sage: C=graphs.PetersenGraph()
            sage: C.cliques_maximal()
            [[0, 1], [0, 4], [0, 5], [1, 2], [1, 6], [2, 3], [2, 7], [3, 4],
             [3, 8], [4, 9], [5, 7], [5, 8], [6, 8], [6, 9], [7, 9]]
            sage: C = Graph('DJ{')
            sage: C.cliques_maximal()
            [[0, 4], [1, 2, 3, 4]]

        Comparing the two implementations::

            sage: g = graphs.RandomGNP(20,.7)
            sage: s1 = Set(map(Set, g.cliques_maximal(algorithm="NetworkX")))
            sage: s2 = Set(map(Set, g.cliques_maximal(algorithm="native")))
            sage: s1 == s2
            True
        """
        if algorithm == "native":
            from sage.graphs.independent_sets import IndependentSets
            return sorted(IndependentSets(self, maximal = True, complement = True))
        elif algorithm == "NetworkX":
            import networkx
            return sorted(networkx.find_cliques(self.networkx_graph(copy=False)))
        else:
            raise ValueError("Algorithm must be equal to 'native' or to 'NetworkX'.")

    @doc_index("Clique-related methods")
    def clique_maximum(self,  algorithm="Cliquer"):
        """
        Returns the vertex set of a maximal order complete subgraph.

        INPUT:

        - ``algorithm`` -- the algorithm to be used :

          - If ``algorithm = "Cliquer"`` (default) - This wraps the C program
            Cliquer [NisOst2003]_.

          - If ``algorithm = "MILP"``, the problem is solved through a Mixed
            Integer Linear Program.

            (see :class:`~sage.numerical.mip.MixedIntegerLinearProgram`)

          - If ``algorithm = "mcqd"`` - Uses the MCQD solver
            (`<http://www.sicmm.org/~konc/maxclique/>`_). Note that the MCQD
            package must be installed.

        .. NOTE::

            Currently only implemented for undirected graphs. Use to_undirected
            to convert a digraph to an undirected graph.

        ALGORITHM:

        This function is based on Cliquer [NisOst2003]_.

        EXAMPLES:

        Using Cliquer (default)::

            sage: C=graphs.PetersenGraph()
            sage: C.clique_maximum()
            [7, 9]
            sage: C = Graph('DJ{')
            sage: C.clique_maximum()
            [1, 2, 3, 4]

        Through a Linear Program::

            sage: len(C.clique_maximum(algorithm = "MILP"))
            4

        TESTS:

        Wrong algorithm::

            sage: C.clique_maximum(algorithm = "BFS")
            Traceback (most recent call last):
            ...
            NotImplementedError: Only 'MILP', 'Cliquer' and 'mcqd' are supported.

        """
        self._scream_if_not_simple(allow_multiple_edges=True)
        if algorithm=="Cliquer":
            from sage.graphs.cliquer import max_clique
            return max_clique(self)
        elif algorithm == "MILP":
            return self.complement().independent_set(algorithm = algorithm)
        elif algorithm == "mcqd":
            try:
                from sage.graphs.mcqd import mcqd
            except ImportError:
                raise ImportError("Please install the mcqd package")
            return mcqd(self)
        else:
            raise NotImplementedError("Only 'MILP', 'Cliquer' and 'mcqd' are supported.")

    @doc_index("Clique-related methods")
    def clique_number(self, algorithm="Cliquer", cliques=None):
        r"""
        Returns the order of the largest clique of the graph (the clique
        number).

        .. NOTE::

            Currently only implemented for undirected graphs. Use ``to_undirected``
            to convert a digraph to an undirected graph.

        INPUT:

        - ``algorithm`` -- the algorithm to be used :

          - If ``algorithm = "Cliquer"`` - This wraps the C program Cliquer
            [NisOst2003]_.

          - If ``algorithm = "networkx"`` - This function is based on
            NetworkX's implementation of the Bron and Kerbosch Algorithm
            [BroKer1973]_.

          - If ``algorithm = "MILP"``, the problem is solved through a Mixed
            Integer Linear Program.

            (see :class:`~sage.numerical.mip.MixedIntegerLinearProgram`)

          - If ``algorithm = "mcqd"`` - Uses the MCQD solver
            (`<http://www.sicmm.org/~konc/maxclique/>`_). Note that the MCQD
            package must be installed.

        - ``cliques`` - an optional list of cliques that can be input if
          already computed. Ignored unless ``algorithm=="networkx"``.

        ALGORITHM:

        This function is based on Cliquer [NisOst2003]_ and [BroKer1973]_.

        EXAMPLES::

            sage: C = Graph('DJ{')
            sage: C.clique_number()
            4
            sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage: G.show(figsize=[2,2])
            sage: G.clique_number()
            3

        By definition the clique number of a complete graph is its order::

            sage: all(graphs.CompleteGraph(i).clique_number() == i for i in xrange(1,15))
            True

        A non-empty graph without edges has a clique number of 1::

            sage: all((i*graphs.CompleteGraph(1)).clique_number() == 1 for i in xrange(1,15))
            True

        A complete multipartite graph with k parts has clique number k::

            sage: all((i*graphs.CompleteMultipartiteGraph(i*[5])).clique_number() == i for i in xrange(1,6))
            True

        TESTS::

            sage: g = graphs.PetersenGraph()
            sage: g.clique_number(algorithm="MILP")
            2
            sage: for i in range(10):                                            # optional - mcqd
            ...       g = graphs.RandomGNP(15,.5)                                # optional - mcqd
            ...       if g.clique_number() != g.clique_number(algorithm="mcqd"): # optional - mcqd
            ...           print("This is dead wrong !")                          # optional - mcqd
        """
        self._scream_if_not_simple(allow_loops=False)
        if algorithm=="Cliquer":
            from sage.graphs.cliquer import clique_number
            return clique_number(self)
        elif algorithm=="networkx":
            import networkx
            return networkx.graph_clique_number(self.networkx_graph(copy=False),cliques)
        elif algorithm == "MILP":
            return len(self.complement().independent_set(algorithm = algorithm))
        elif algorithm == "mcqd":
            try:
                from sage.graphs.mcqd import mcqd
            except ImportError:
                raise ImportError("Please install the mcqd package")
            return len(mcqd(self))
        else:
            raise NotImplementedError("Only 'networkx' 'MILP' 'Cliquer' and 'mcqd' are supported.")

    @doc_index("Clique-related methods")
    def cliques_number_of(self, vertices=None, cliques=None):
        """
        Returns a dictionary of the number of maximal cliques containing each
        vertex, keyed by vertex. (Returns a single value if
        only one input vertex).

        .. NOTE::

            Currently only implemented for undirected graphs. Use to_undirected
            to convert a digraph to an undirected graph.

        INPUT:

        -  ``vertices`` - the vertices to inspect (default is
           entire graph)

        -  ``cliques`` - list of cliques (if already
           computed)


        EXAMPLES::

            sage: C = Graph('DJ{')
            sage: C.cliques_number_of()
            {0: 1, 1: 1, 2: 1, 3: 1, 4: 2}
            sage: E = C.cliques_maximal()
            sage: E
            [[0, 4], [1, 2, 3, 4]]
            sage: C.cliques_number_of(cliques=E)
            {0: 1, 1: 1, 2: 1, 3: 1, 4: 2}
            sage: F = graphs.Grid2dGraph(2,3)
            sage: X = F.cliques_number_of()
            sage: for v in sorted(X.iterkeys()):
            ....:     print("{} {}".format(v, X[v]))
            (0, 0) 2
            (0, 1) 3
            (0, 2) 2
            (1, 0) 2
            (1, 1) 3
            (1, 2) 2
            sage: F.cliques_number_of(vertices=[(0, 1), (1, 2)])
            {(0, 1): 3, (1, 2): 2}
            sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage: G.show(figsize=[2,2])
            sage: G.cliques_number_of()
            {0: 2, 1: 2, 2: 1, 3: 1}
        """
        import networkx
        return networkx.number_of_cliques(self.networkx_graph(copy=False), vertices, cliques)

    @doc_index("Clique-related methods")
    def cliques_get_max_clique_graph(self, name=''):
        """
        Returns a graph constructed with maximal cliques as vertices, and
        edges between maximal cliques with common members in the original
        graph.

        For more information, see the :wikipedia:`Clique_graph`.

        .. NOTE::

            Currently only implemented for undirected graphs. Use to_undirected
            to convert a digraph to an undirected graph.

        INPUT:

        -  ``name`` - The name of the new graph.

        EXAMPLES::

            sage: (graphs.ChvatalGraph()).cliques_get_max_clique_graph()
            Graph on 24 vertices
            sage: ((graphs.ChvatalGraph()).cliques_get_max_clique_graph()).show(figsize=[2,2], vertex_size=20, vertex_labels=False)
            sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage: G.show(figsize=[2,2])
            sage: G.cliques_get_max_clique_graph()
            Graph on 2 vertices
            sage: (G.cliques_get_max_clique_graph()).show(figsize=[2,2])
        """
        import networkx
        return Graph(networkx.make_max_clique_graph(self.networkx_graph(copy=False), name=name, create_using=networkx.MultiGraph()))

    @doc_index("Clique-related methods")
    def cliques_get_clique_bipartite(self, **kwds):
        """
        Returns a bipartite graph constructed such that maximal cliques are the
        right vertices and the left vertices are retained from the given
        graph. Right and left vertices are connected if the bottom vertex
        belongs to the clique represented by a top vertex.

        .. NOTE::

            Currently only implemented for undirected graphs. Use to_undirected
            to convert a digraph to an undirected graph.

        EXAMPLES::

            sage: (graphs.ChvatalGraph()).cliques_get_clique_bipartite()
            Bipartite graph on 36 vertices
            sage: ((graphs.ChvatalGraph()).cliques_get_clique_bipartite()).show(figsize=[2,2], vertex_size=20, vertex_labels=False)
            sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage: G.show(figsize=[2,2])
            sage: G.cliques_get_clique_bipartite()
            Bipartite graph on 6 vertices
            sage: (G.cliques_get_clique_bipartite()).show(figsize=[2,2])
        """
        from .bipartite_graph import BipartiteGraph
        import networkx
        return BipartiteGraph(networkx.make_clique_bipartite(self.networkx_graph(copy=False), **kwds))

    @doc_index("Algorithmically hard stuff")
    def independent_set(self, algorithm = "Cliquer", value_only = False, reduction_rules = True, solver = None, verbosity = 0):
        r"""
        Returns a maximum independent set.

        An independent set of a graph is a set of pairwise non-adjacent
        vertices. A maximum independent set is an independent set of maximum
        cardinality.  It induces an empty subgraph.

        Equivalently, an independent set is defined as the complement of a
        vertex cover.

        For more information, see the
        :wikipedia:`Independent_set_(graph_theory)` and the
        :wikipedia:`Vertex_cover`.

        INPUT:

        - ``algorithm`` -- the algorithm to be used

          * If ``algorithm = "Cliquer"`` (default), the problem is solved
            using Cliquer [NisOst2003]_.

            (see the :mod:`Cliquer modules <sage.graphs.cliquer>`)

          * If ``algorithm = "MILP"``, the problem is solved through a Mixed
            Integer Linear Program.

            (see :class:`~sage.numerical.mip.MixedIntegerLinearProgram`)

         * If ``algorithm = "mcqd"`` - Uses the MCQD solver
           (`<http://www.sicmm.org/~konc/maxclique/>`_). Note that the MCQD
           package must be installed.

        - ``value_only`` -- boolean (default: ``False``). If set to ``True``,
          only the size of a maximum independent set is returned. Otherwise,
          a maximum independent set is returned as a list of vertices.

        - ``reduction_rules`` -- (default: ``True``) Specify if the reductions
          rules from kernelization must be applied as pre-processing or not.
          See [ACFLSS04]_ for more details. Note that depending on the
          instance, it might be faster to disable reduction rules.

        - ``solver`` -- (default: ``None``) Specify a Linear Program (LP)
          solver to be used. If set to ``None``, the default one is used. For
          more information on LP solvers and which default solver is used, see
          the method
          :meth:`~sage.numerical.mip.MixedIntegerLinearProgram.solve`
          of the class
          :class:`~sage.numerical.mip.MixedIntegerLinearProgram`.

        - ``verbosity`` -- non-negative integer (default: ``0``). Set the level
          of verbosity you want from the linear program solver. Since the
          problem of computing an independent set is `NP`-complete, its solving
          may take some time depending on the graph. A value of 0 means that
          there will be no message printed by the solver. This option is only
          useful if ``algorithm="MILP"``.

        .. NOTE::

            While Cliquer/MCAD are usually (and by far) the most efficient
            implementations, the MILP formulation sometimes proves faster on
            very "symmetrical" graphs.

        EXAMPLES:

        Using Cliquer::

            sage: C = graphs.PetersenGraph()
            sage: C.independent_set()
            [0, 3, 6, 7]

        As a linear program::

            sage: C = graphs.PetersenGraph()
            sage: len(C.independent_set(algorithm = "MILP"))
            4

        .. PLOT::

            g = graphs.PetersenGraph()
            sphinx_plot(g.plot(partition=[g.independent_set()]))
        """
        my_cover = self.vertex_cover(algorithm=algorithm, value_only=value_only, reduction_rules=reduction_rules, solver=solver, verbosity=verbosity)
        if value_only:
            return self.order() - my_cover
        else:
            return [u for u in self.vertices() if not u in my_cover]


    @doc_index("Algorithmically hard stuff")
    def vertex_cover(self, algorithm = "Cliquer", value_only = False,
                     reduction_rules = True, solver = None, verbosity = 0):
        r"""
        Returns a minimum vertex cover of self represented by a set of vertices.

        A minimum vertex cover of a graph is a set `S` of vertices such that
        each edge is incident to at least one element of `S`, and such that `S`
        is of minimum cardinality. For more information, see the
        :wikipedia:`Wikipedia article on vertex cover <Vertex_cover>`.

        Equivalently, a vertex cover is defined as the complement of an
        independent set.

        As an optimization problem, it can be expressed as follows:

        .. MATH::

            \mbox{Minimize : }&\sum_{v\in G} b_v\\
            \mbox{Such that : }&\forall (u,v) \in G.edges(), b_u+b_v\geq 1\\
            &\forall x\in G, b_x\mbox{ is a binary variable}

        INPUT:

        - ``algorithm`` -- string (default: ``"Cliquer"``). Indicating
          which algorithm to use. It can be one of those two values.

          - ``"Cliquer"`` will compute a minimum vertex cover
            using the Cliquer package.

          - ``"MILP"`` will compute a minimum vertex cover through a mixed
            integer linear program.

          - If ``algorithm = "mcqd"`` - Uses the MCQD solver
            (`<http://www.sicmm.org/~konc/maxclique/>`_). Note that the MCQD
            package must be installed.

        - ``value_only`` -- boolean (default: ``False``). If set to ``True``,
          only the size of a minimum vertex cover is returned. Otherwise,
          a minimum vertex cover is returned as a list of vertices.

        - ``reduction_rules`` -- (default: ``True``) Specify if the reductions
          rules from kernelization must be applied as pre-processing or not.
          See [ACFLSS04]_ for more details. Note that depending on the
          instance, it might be faster to disable reduction rules.

        - ``solver`` -- (default: ``None``) Specify a Linear Program (LP)
          solver to be used. If set to ``None``, the default one is used. For
          more information on LP solvers and which default solver is used, see
          the method
          :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>`
          of the class
          :class:`MixedIntegerLinearProgram <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbosity`` -- non-negative integer (default: ``0``). Set the level
          of verbosity you want from the linear program solver. Since the
          problem of computing a vertex cover is `NP`-complete, its solving may
          take some time depending on the graph. A value of 0 means that there
          will be no message printed by the solver. This option is only useful
          if ``algorithm="MILP"``.

        EXAMPLES:

        On the Pappus graph::

           sage: g = graphs.PappusGraph()
           sage: g.vertex_cover(value_only=True)
           9

        .. PLOT::

            g = graphs.PappusGraph()
            sphinx_plot(g.plot(partition=[g.vertex_cover()]))

        TESTS:

        The two algorithms should return the same result::

           sage: g = graphs.RandomGNP(10,.5)
           sage: vc1 = g.vertex_cover(algorithm="MILP")
           sage: vc2 = g.vertex_cover(algorithm="Cliquer")
           sage: len(vc1) == len(vc2)
           True

        The cardinality of the vertex cover is unchanged when reduction rules are used. First for trees::

           sage: for i in range(20):
           ...       g = graphs.RandomTree(20)
           ...       vc1_set = g.vertex_cover()
           ...       vc1 = len(vc1_set)
           ...       vc2 = g.vertex_cover(value_only = True, reduction_rules = False)
           ...       if vc1 != vc2:
           ...           print("Error :", vc1, vc2)
           ...           print("With reduction rules :", vc1)
           ...           print("Without reduction rules :", vc2)
           ...           break
           ...       g.delete_vertices(vc1_set)
           ...       if g.size() != 0:
           ...           print("This thing is not a vertex cover !")

        Then for random GNP graphs::

           sage: for i in range(20):
           ...       g = graphs.RandomGNP(50,4/50)
           ...       vc1_set = g.vertex_cover()
           ...       vc1 = len(vc1_set)
           ...       vc2 = g.vertex_cover(value_only = True, reduction_rules = False)
           ...       if vc1 != vc2:
           ...           print("Error :", vc1, vc2)
           ...           print("With reduction rules :", vc1)
           ...           print("Without reduction rules :", vc2)
           ...           break
           ...       g.delete_vertices(vc1_set)
           ...       if g.size() != 0:
           ...           print("This thing is not a vertex cover !")

        Testing mcqd::

            sage: graphs.PetersenGraph().vertex_cover(algorithm="mcqd",value_only=True) # optional - mcqd
            6

        Given a wrong algorithm::

            sage: graphs.PetersenGraph().vertex_cover(algorithm = "guess")
            Traceback (most recent call last):
            ...
            ValueError: The algorithm must be "Cliquer" "MILP" or "mcqd".

        REFERENCE:

        .. [ACFLSS04] \F. N. Abu-Khzam, R. L. Collins, M. R. Fellows, M. A.
          Langston, W. H. Suters, and C. T. Symons: Kernelization Algorithm for
          the Vertex Cover Problem: Theory and Experiments. *SIAM ALENEX/ANALCO*
          2004: 62-69.
        """
        self._scream_if_not_simple(allow_multiple_edges=True)
        g = self

        ppset = []
        folded_vertices = []

        ###################
        # Reduction rules #
        ###################

        if reduction_rules:
            # We apply simple reduction rules allowing to identify vertices that
            # belongs to an optimal vertex cover

            # We first create manually a copy of the graph to prevent creating
            # multi-edges when merging vertices, if edges have labels (e.g., weights).
            g = copy(self)

            degree_at_most_two = set([u for u,du in g.degree(labels = True).items() if du <= 2])

            while degree_at_most_two:

                u = degree_at_most_two.pop()
                du = g.degree(u)

                if du == 0:
                    # RULE 1: isolated vertices are not part of the cover. We
                    # simply remove them from the graph. The degree of such
                    # vertices may have been reduced to 0 while applying other
                    # reduction rules
                    g.delete_vertex(u)

                elif du == 1:
                    # RULE 2: If a vertex u has degree 1, we select its neighbor
                    # v and remove both u and v from g.
                    v = g.neighbors(u)[0]
                    ppset.append(v)
                    g.delete_vertex(u)

                    for w in g.neighbors(v):
                        if g.degree(w) <= 3:
                            # The degree of w will be at most two after the
                            # deletion of v
                            degree_at_most_two.add(w)

                    g.delete_vertex(v)
                    degree_at_most_two.discard(v)

                elif du == 2:
                    v,w  = g.neighbors(u)

                    if g.has_edge(v,w):
                        # RULE 3: If the neighbors v and w of a degree 2 vertex
                        # u are incident, then we select both v and w and remove
                        # u, v, and w from g.
                        ppset.append(v)
                        ppset.append(w)
                        g.delete_vertex(u)
                        neigh = set(g.neighbors(v) + g.neighbors(w)).difference(set([v,w]))
                        g.delete_vertex(v)
                        g.delete_vertex(w)

                        for z in neigh:
                            if g.degree(z) <= 2:
                                degree_at_most_two.add(z)

                    else:
                        # RULE 4, folded vertices: If the neighbors v and w of a
                        # degree 2 vertex u are not incident, then we contract
                        # edges (u, v), (u,w). Then, if the solution contains u,
                        # we replace it with v and w. Otherwise, we let u in the
                        # solution.
                        neigh = set(g.neighbors(v) + g.neighbors(w)).difference(set([u,v,w]))
                        g.delete_vertex(v)
                        g.delete_vertex(w)
                        for z in neigh:
                            g.add_edge(u,z)

                        folded_vertices += [(u,v,w)]

                        if g.degree(u) <= 2:
                            degree_at_most_two.add(u)

                    degree_at_most_two.discard(v)
                    degree_at_most_two.discard(w)


                # RULE 5:
                # TODO: add extra reduction rules


        ##################
        # Main Algorithm #
        ##################

        if g.order() == 0:
            # Reduction rules were sufficients to get the solution
            size_cover_g = 0
            cover_g = []

        elif algorithm == "Cliquer" or algorithm == "mcqd":
            independent = g.complement().clique_maximum(algorithm=algorithm)
            if value_only:
                size_cover_g = g.order() - len(independent)
            else:
                cover_g = [u for u in g.vertices() if not u in independent]

        elif algorithm == "MILP":

            from sage.numerical.mip import MixedIntegerLinearProgram
            p = MixedIntegerLinearProgram(maximization=False, solver=solver)
            b = p.new_variable(binary=True)

            # minimizes the number of vertices in the set
            p.set_objective(p.sum(b[v] for v in g.vertices()))

            # an edge contains at least one vertex of the minimum vertex cover
            for (u,v) in g.edges(labels=None):
                p.add_constraint(b[u] + b[v], min=1)

            if value_only:
                size_cover_g = p.solve(objective_only=True, log=verbosity)
            else:
                p.solve(log=verbosity)
                b = p.get_values(b)
                cover_g = [v for v in g.vertices() if b[v] == 1]
        else:
            raise ValueError("The algorithm must be \"Cliquer\" \"MILP\" or \"mcqd\".")

        #########################
        # Returning the results #
        #########################

        # We finally reconstruct the solution according the reduction rules
        if value_only:
            return len(ppset) + len(folded_vertices) + size_cover_g
        else:
            # RULES 2 and 3:
            cover_g.extend(ppset)
            # RULE 4:
            folded_vertices.reverse()
            for u,v,w in folded_vertices:
                if u in cover_g:
                    cover_g.remove(u)
                    cover_g += [v,w]
                else:
                    cover_g += [u]
            cover_g.sort()
            return cover_g

    @doc_index("Clique-related methods")
    def cliques_vertex_clique_number(self, algorithm="cliquer", vertices=None,
                                     cliques=None):
        """
        Returns a dictionary of sizes of the largest maximal cliques containing
        each vertex, keyed by vertex. (Returns a single value if only one
        input vertex).

        .. NOTE::

            Currently only implemented for undirected graphs. Use to_undirected
            to convert a digraph to an undirected graph.

        INPUT:

         - ``algorithm`` - either ``cliquer`` or ``networkx``

           - ``cliquer`` - This wraps the C program Cliquer [NisOst2003]_.

           - ``networkx`` - This function is based on NetworkX's implementation
             of the Bron and Kerbosch Algorithm [BroKer1973]_.

        -  ``vertices`` - the vertices to inspect (default is entire graph).
           Ignored unless ``algorithm=='networkx'``.

        -  ``cliques`` - list of cliques (if already computed).  Ignored unless
           ``algorithm=='networkx'``.

        EXAMPLES::

            sage: C = Graph('DJ{')
            sage: C.cliques_vertex_clique_number()
            {0: 2, 1: 4, 2: 4, 3: 4, 4: 4}
            sage: E = C.cliques_maximal()
            sage: E
            [[0, 4], [1, 2, 3, 4]]
            sage: C.cliques_vertex_clique_number(cliques=E,algorithm="networkx")
            {0: 2, 1: 4, 2: 4, 3: 4, 4: 4}
            sage: F = graphs.Grid2dGraph(2,3)
            sage: X = F.cliques_vertex_clique_number(algorithm="networkx")
            sage: for v in sorted(X.iterkeys()):
            ....:     print("{} {}".format(v, X[v]))
            (0, 0) 2
            (0, 1) 2
            (0, 2) 2
            (1, 0) 2
            (1, 1) 2
            (1, 2) 2
            sage: F.cliques_vertex_clique_number(vertices=[(0, 1), (1, 2)])
            {(0, 1): 2, (1, 2): 2}
            sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage: G.show(figsize=[2,2])
            sage: G.cliques_vertex_clique_number()
            {0: 3, 1: 3, 2: 3, 3: 3}

        """

        if algorithm=="cliquer":
            from sage.graphs.cliquer import clique_number
            if vertices is None:
                vertices=self
            value={}
            for v in vertices:
                value[v] = 1+clique_number(self.subgraph(self.neighbors(v)))
                self.subgraph(self.neighbors(v)).plot()
            return value
        elif algorithm=="networkx":
            import networkx
            return networkx.node_clique_number(self.networkx_graph(copy=False),vertices, cliques)
        else:
            raise NotImplementedError("Only 'networkx' and 'cliquer' are supported.")

    @doc_index("Clique-related methods")
    def cliques_containing_vertex(self, vertices=None, cliques=None):
        """
        Returns the cliques containing each vertex, represented as a dictionary
        of lists of lists, keyed by vertex. (Returns a single list if only one
        input vertex).

        .. NOTE::

            Currently only implemented for undirected graphs. Use to_undirected
            to convert a digraph to an undirected graph.

        INPUT:

        -  ``vertices`` - the vertices to inspect (default is
           entire graph)

        -  ``cliques`` - list of cliques (if already
           computed)

        EXAMPLES::

            sage: C = Graph('DJ{')
            sage: C.cliques_containing_vertex()
            {0: [[4, 0]], 1: [[4, 1, 2, 3]], 2: [[4, 1, 2, 3]], 3: [[4, 1, 2, 3]], 4: [[4, 0], [4, 1, 2, 3]]}
            sage: E = C.cliques_maximal()
            sage: E
            [[0, 4], [1, 2, 3, 4]]
            sage: C.cliques_containing_vertex(cliques=E)
            {0: [[0, 4]], 1: [[1, 2, 3, 4]], 2: [[1, 2, 3, 4]], 3: [[1, 2, 3, 4]], 4: [[0, 4], [1, 2, 3, 4]]}
            sage: F = graphs.Grid2dGraph(2,3)
            sage: X = F.cliques_containing_vertex()
            sage: for v in sorted(X.iterkeys()):
            ....:     print("{} {}".format(v, X[v]))
            (0, 0) [[(0, 1), (0, 0)], [(1, 0), (0, 0)]]
            (0, 1) [[(0, 1), (0, 0)], [(0, 1), (0, 2)], [(0, 1), (1, 1)]]
            (0, 2) [[(0, 1), (0, 2)], [(1, 2), (0, 2)]]
            (1, 0) [[(1, 0), (0, 0)], [(1, 0), (1, 1)]]
            (1, 1) [[(0, 1), (1, 1)], [(1, 2), (1, 1)], [(1, 0), (1, 1)]]
            (1, 2) [[(1, 2), (0, 2)], [(1, 2), (1, 1)]]
            sage: F.cliques_containing_vertex(vertices=[(0, 1), (1, 2)])
            {(0, 1): [[(0, 1), (0, 0)], [(0, 1), (0, 2)], [(0, 1), (1, 1)]], (1, 2): [[(1, 2), (0, 2)], [(1, 2), (1, 1)]]}
            sage: G = Graph({0:[1,2,3], 1:[2], 3:[0,1]})
            sage: G.show(figsize=[2,2])
            sage: G.cliques_containing_vertex()
            {0: [[0, 1, 2], [0, 1, 3]], 1: [[0, 1, 2], [0, 1, 3]], 2: [[0, 1, 2]], 3: [[0, 1, 3]]}

        """
        import networkx
        return networkx.cliques_containing_node(self.networkx_graph(copy=False),vertices, cliques)

    @doc_index("Clique-related methods")
    def clique_complex(self):
        """
        Returns the clique complex of self. This is the largest simplicial complex on
        the vertices of self whose 1-skeleton is self.

        This is only makes sense for undirected simple graphs.

        EXAMPLES::

            sage: g = Graph({0:[1,2],1:[2],4:[]})
            sage: g.clique_complex()
            Simplicial complex with vertex set (0, 1, 2, 4) and facets {(4,), (0, 1, 2)}

            sage: h = Graph({0:[1,2,3,4],1:[2,3,4],2:[3]})
            sage: x = h.clique_complex()
            sage: x
            Simplicial complex with vertex set (0, 1, 2, 3, 4) and facets {(0, 1, 4), (0, 1, 2, 3)}
            sage: i = x.graph()
            sage: i==h
            True
            sage: x==i.clique_complex()
            True

        """
        if self.is_directed() or self.has_loops() or self.has_multiple_edges():
            raise ValueError("Self must be an undirected simple graph to have a clique complex.")
        import sage.homology.simplicial_complex
        C = sage.homology.simplicial_complex.SimplicialComplex(self.cliques_maximal(), maximality_check=True)
        C._graph = self
        return C

    @doc_index("Clique-related methods")
    def clique_polynomial(self, t = None):
        """
        Returns the clique polynomial of self.

        This is the polynomial where the coefficient of `t^n` is the number of
        cliques in the graph with `n` vertices. The constant term of the
        clique polynomial is always taken to be one.

        EXAMPLES::

            sage: g = Graph()
            sage: g.clique_polynomial()
            1
            sage: g = Graph({0:[1]})
            sage: g.clique_polynomial()
            t^2 + 2*t + 1
            sage: g = graphs.CycleGraph(4)
            sage: g.clique_polynomial()
            4*t^2 + 4*t + 1

        """
        if t is None:
            R = PolynomialRing(ZZ, 't')
            t = R.gen()
        number_of = [0]*(self.order() + 1)
        for x in IndependentSets(self, complement = True):
            number_of[len(x)] += 1
        return sum(coeff*t**i for i,coeff in enumerate(number_of) if coeff)

    ### Miscellaneous

    @doc_index("Leftovers")
    def cores(self, k = None, with_labels=False):
        """
        Returns the core number for each vertex in an ordered list.

        (for homomorphisms cores, see the :meth:`Graph.has_homomorphism_to`
        method)

        **DEFINITIONS**

        * *K-cores* in graph theory were introduced by Seidman in 1983 and by
          Bollobas in 1984 as a method of (destructively) simplifying graph
          topology to aid in analysis and visualization. They have been more
          recently defined as the following by Batagelj et al:

          *Given a graph `G` with vertices set `V` and edges set `E`, the
          `k`-core of `G` is the graph obtained from `G` by recursively removing
          the vertices with degree less than `k`, for as long as there are any.*

          This operation can be useful to filter or to study some properties of
          the graphs. For instance, when you compute the 2-core of graph G, you
          are cutting all the vertices which are in a tree part of graph.  (A
          tree is a graph with no loops). [WPkcore]_

          [PSW1996]_ defines a `k`-core of `G` as the largest subgraph (it is
          unique) of `G` with minimum degree at least `k`.

        * Core number of a vertex

          The core number of a vertex `v` is the largest integer `k` such that
          `v` belongs to the `k`-core of `G`.

        * Degeneracy

          The *degeneracy* of a graph `G`, usually denoted `\delta^*(G)`, is the
          smallest integer `k` such that the graph `G` can be reduced to the
          empty graph by iteratively removing vertices of degree `\leq
          k`. Equivalently, `\delta^*(G)=k` if `k` is the smallest integer such
          that the `k`-core of `G` is empty.

        **IMPLEMENTATION**

        This implementation is based on the NetworkX implementation of
        the algorithm described in [BZ]_.

        **INPUT**

        - ``k`` (integer)

            * If ``k = None`` (default), returns the core number for each vertex.

            * If ``k`` is an integer, returns a pair ``(ordering, core)``, where
              ``core`` is the list of vertices in the `k`-core of ``self``, and
              ``ordering`` is an elimination order for the other vertices such
              that each vertex is of degree strictly less than `k` when it is to
              be eliminated from the graph.

        - ``with_labels`` (boolean)

           * When set to ``False``, and ``k = None``, the method returns a list
             whose `i` th element is the core number of the `i` th vertex. When
             set to ``True``, the method returns a dictionary whose keys are
             vertices, and whose values are the corresponding core numbers.

             By default, ``with_labels = False``.

        .. SEEALSO::

           * Graph cores is also a notion related to graph homomorphisms. For
             this second meaning, see :meth:`Graph.has_homomorphism_to`.

        REFERENCE:

        .. [WPkcore] K-core. Wikipedia. (2007). [Online] Available:
          :wikipedia:`K-core`

        .. [PSW1996] Boris Pittel, Joel Spencer and Nicholas Wormald. Sudden
          Emergence of a Giant k-Core in a Random
          Graph. (1996). J. Combinatorial Theory. Ser B 67. pages
          111-151. [Online] Available:
          http://cs.nyu.edu/cs/faculty/spencer/papers/k-core.pdf

        .. [BZ] Vladimir Batagelj and Matjaz Zaversnik. An `O(m)`
          Algorithm for Cores Decomposition of
          Networks. :arxiv:`cs/0310049v1`.

        EXAMPLES::

            sage: (graphs.FruchtGraph()).cores()
            [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
            sage: (graphs.FruchtGraph()).cores(with_labels=True)
            {0: 3, 1: 3, 2: 3, 3: 3, 4: 3, 5: 3, 6: 3, 7: 3, 8: 3, 9: 3, 10: 3, 11: 3}
            sage: a=random_matrix(ZZ,20,x=2,sparse=True, density=.1)
            sage: b=Graph(20)
            sage: b.add_edges(a.nonzero_positions())
            sage: cores=b.cores(with_labels=True); cores
            {0: 3, 1: 3, 2: 3, 3: 3, 4: 2, 5: 2, 6: 3, 7: 1, 8: 3, 9: 3, 10: 3, 11: 3, 12: 3, 13: 3, 14: 2, 15: 3, 16: 3, 17: 3, 18: 3, 19: 3}
            sage: [v for v,c in cores.items() if c>=2] # the vertices in the 2-core
            [0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]

        Checking the 2-core of a random lobster is indeed the empty set::

            sage: g = graphs.RandomLobster(20,.5,.5)
            sage: ordering, core = g.cores(2)
            sage: len(core) == 0
            True
        """
        self._scream_if_not_simple()
        # compute the degrees of each vertex
        degrees=self.degree(labels=True)

        # sort vertices by degree.  Store in a list and keep track of
        # where a specific degree starts (effectively, the list is
        # sorted by bins).
        verts= sorted( degrees.keys(), key=lambda x: degrees[x])
        bin_boundaries=[0]
        curr_degree=0
        for i,v in enumerate(verts):
            if degrees[v]>curr_degree:
                bin_boundaries.extend([i]*(degrees[v]-curr_degree))
                curr_degree=degrees[v]
        vert_pos = dict((v,pos) for pos,v in enumerate(verts))
        # Set up initial guesses for core and lists of neighbors.
        core= degrees
        nbrs=dict((v,set(self.neighbors(v))) for v in self)
        # form vertex core building up from smallest
        for v in verts:

            # If all the vertices have a degree larger than k, we can
            # return our answer if k is not None
            if k is not None and core[v] >= k:
                return verts[:vert_pos[v]], verts[vert_pos[v]:]

            for u in nbrs[v]:
                if core[u] > core[v]:
                    nbrs[u].remove(v)

                    # cleverly move u to the end of the next smallest
                    # bin (i.e., subtract one from the degree of u).
                    # We do this by swapping u with the first vertex
                    # in the bin that contains u, then incrementing
                    # the bin boundary for the bin that contains u.
                    pos=vert_pos[u]
                    bin_start=bin_boundaries[core[u]]
                    vert_pos[u]=bin_start
                    vert_pos[verts[bin_start]]=pos
                    verts[bin_start],verts[pos]=verts[pos],verts[bin_start]
                    bin_boundaries[core[u]]+=1
                    core[u] -= 1

        if k is not None:
            return verts, []

        if with_labels:
            return core
        else:
            return core.values()

    @doc_index("Leftovers")
    def modular_decomposition(self):
        r"""
        Returns the modular decomposition of the current graph.

        .. NOTE::

            In order to use this method you must install the
            ``modular_decomposition`` optional package. See
            :mod:`sage.misc.package`.

        Crash course on modular decomposition:

        A module `M` of a graph `G` is a proper subset of its vertices
        such that for all `u \in V(G)-M, v,w\in M` the relation `u
        \sim v \Leftrightarrow u \sim w` holds, where `\sim` denotes
        the adjacency relation in `G`. Equivalently, `M \subset V(G)`
        is a module if all its vertices have the same adjacency
        relations with each vertex outside of the module (vertex by
        vertex).

        Hence, for a set like a module, it is very easy to encode the
        information of the adjacencies between the vertices inside and
        outside the module -- we can actually add a new vertex `v_M`
        to our graph representing our module `M`, and let `v_M` be
        adjacent to `u\in V(G)-M` if and only if some `v\in M` (and
        hence all the vertices contained in the module) is adjacent to
        `u`. We can now independently (and recursively) study the
        structure of our module `M` and the new graph `G-M+\{v_M\}`,
        without any loss of information.

        Here are two very simple modules :

            * A connected component `C` (or the union of some --but
              not all-- of them) of a disconnected graph `G`, for
              instance, is a module, as no vertex of `C` has a
              neighbor outside of it.

            * An anticomponent `C` (or the union of some --but not
              all-- of them) of an non-anticonnected graph `G`, for
              the same reason (it is just the complement of the
              previous graph !).

        These modules being of special interest, the disjoint union of
        graphs is called a Parallel composition, and the complement of
        a disjoint union is called a Parallel composition. A graph
        whose only modules are singletons is called Prime.

        For more information on modular decomposition, in particular
        for an explanation of the terms "Parallel," "Prime" and
        "Serie," see the `Wikipedia article on modular decomposition
        <http://en.wikipedia.org/wiki/Modular_decomposition>`_.

        You may also be interested in the survey from Michel Habib and
        Christophe Paul entitled "A survey on Algorithmic aspects of
        modular decomposition" [HabPau10]_.

        OUTPUT:

        A pair of two values (recursively encoding the decomposition) :

        * The type of the current module :

          * ``"Parallel"``
          * ``"Prime"``
          * ``"Serie"``

        * The list of submodules (as list of pairs ``(type, list)``,
          recursively...) or the vertex's name if the module is a
          singleton.

        EXAMPLES:

        The Bull Graph is prime::

            sage: graphs.BullGraph().modular_decomposition() # optional -- modular_decomposition
            ('Prime', [3, 4, 0, 1, 2])

        The Petersen Graph too::

            sage: graphs.PetersenGraph().modular_decomposition() # optional -- modular_decomposition
            ('Prime', [2, 6, 3, 9, 7, 8, 0, 1, 5, 4])

        This a clique on 5 vertices with 2 pendant edges, though, has a more
        interesting decomposition ::

            sage: g = graphs.CompleteGraph(5)
            sage: g.add_edge(0,5)
            sage: g.add_edge(0,6)
            sage: g.modular_decomposition() # optional -- modular_decomposition
            ('Serie', [0, ('Parallel', [5, ('Serie', [1, 4, 3, 2]), 6])])

        ALGORITHM:

        This function uses a C implementation of a 2-step algorithm
        implemented by Fabien de Montgolfier [FMDec]_ :

            * Computation of a factorizing permutation [HabibViennot1999]_.

            * Computation of the tree itself [CapHabMont02]_.

        .. SEEALSO::

            - :meth:`is_prime` -- Tests whether a graph is prime.

        REFERENCE:

        .. [FMDec] Fabien de Montgolfier
          http://www.liafa.jussieu.fr/~fm/algos/index.html

        .. [HabibViennot1999] Michel Habib, Christiphe Paul, Laurent Viennot
          Partition refinement techniques: An interesting algorithmic tool kit
          International Journal of Foundations of Computer Science
          vol. 10 n2 pp.147--170, 1999

        .. [CapHabMont02] \C. Capelle, M. Habib et F. de Montgolfier
          Graph decomposition and Factorising Permutations
          Discrete Mathematics and Theoretical Computer Sciences, vol 5 no. 1 , 2002.

        .. [HabPau10] Michel Habib and Christophe Paul
          A survey of the algorithmic aspects of modular decomposition
          Computer Science Review
          vol 4, number 1, pages 41--59, 2010
          http://www.lirmm.fr/~paul/md-survey.pdf
        """
        try:
            from sage.graphs.modular_decomposition import modular_decomposition
        except ImportError:
            raise RuntimeError("In order to use this method you must "
                               "install the modular_decomposition package")

        self._scream_if_not_simple()
        from sage.misc.stopgap import stopgap
        stopgap("Graph.modular_decomposition is known to return wrong results",13744)

        D = modular_decomposition(self)

        id_label = dict(enumerate(self.vertices()))

        relabel = lambda x : (x[0], [relabel(_) for _ in x[1]]) if isinstance(x,tuple) else id_label[x]

        return relabel(D)

    @doc_index("Graph properties")
    def is_prime(self):
        r"""
        Tests whether the current graph is prime.

        A graph is prime if all its modules are trivial (i.e. empty, all of the
        graph or singletons) -- see :meth:`modular_decomposition`.

        .. NOTE::

            In order to use this method you must install the
            ``modular_decomposition`` optional package. See
            :mod:`sage.misc.package`.

        EXAMPLE:

        The Petersen Graph and the Bull Graph are both prime::

            sage: graphs.PetersenGraph().is_prime() # optional - modular_decomposition
            True
            sage: graphs.BullGraph().is_prime()     # optional - modular_decomposition
            True

        Though quite obviously, the disjoint union of them is not::

            sage: (graphs.PetersenGraph() + graphs.BullGraph()).is_prime() # optional - modular_decomposition
            False
        """

        D = self.modular_decomposition()

        return D[0] == "Prime" and len(D[1]) == self.order()

    @rename_keyword(deprecation=19550, method='algorithm')
    def _gomory_hu_tree(self, vertices, algorithm=None):
        r"""
        Return a Gomory-Hu tree associated to self.

        This function is the private counterpart of ``gomory_hu_tree()``,
        with the difference that it has an optional argument
        needed for recursive computations, which the user is not
        interested in defining himself.

        See the documentation of ``gomory_hu_tree()`` for more information.

        INPUT:

        - ``vertices`` - a set of "real" vertices, as opposed to the
          fakes one introduced during the computations. This variable is
          useful for the algorithm and for recursion purposes.

        - ``algorithm`` -- select the algorithm used by the :meth:`edge_cut`
          method. Refer to its documentation for allowed values and default
          behaviour.

        EXAMPLE:

        This function is actually tested in ``gomory_hu_tree()``, this
        example is only present to have a doctest coverage of 100%.

            sage: g = graphs.PetersenGraph()
            sage: t = g._gomory_hu_tree(frozenset(g.vertices()))
        """
        self._scream_if_not_simple()

        # Small case, not really a problem ;-)
        if len(vertices) == 1:
            g = Graph()
            g.add_vertices(vertices)
            return g

        # Take any two vertices (u,v)
        it = iter(vertices)
        u,v = next(it),next(it)

        # Compute a uv min-edge-cut.
        #
        # The graph is split into U,V with u \in U and v\in V.
        flow,edges,[U,V] = self.edge_cut(u, v, use_edge_labels=True, vertices=True, algorithm=algorithm)

        # One graph for each part of the previous one
        gU,gV = self.subgraph(U, immutable=False), self.subgraph(V, immutable=False)

        # A fake vertex fU (resp. fV) to represent U (resp. V)
        fU = frozenset(U)
        fV = frozenset(V)

        # Each edge (uu,vv) with uu \in U and vv\in V yields:
        # - an edge (uu,fV) in gU
        # - an edge (vv,fU) in gV
        #
        # If the same edge is added several times their capacities add up.

        from sage.rings.real_mpfr import RR
        for uu,vv,capacity in edges:
            capacity = capacity if capacity in RR else 1

            # Assume uu is in gU
            if uu in V:
                uu,vv = vv,uu

            # Create the new edges if necessary
            if not gU.has_edge(uu, fV):
                gU.add_edge(uu, fV, 0)
            if not gV.has_edge(vv, fU):
                gV.add_edge(vv, fU, 0)

            # update the capacities
            gU.set_edge_label(uu, fV, gU.edge_label(uu, fV) + capacity)
            gV.set_edge_label(vv, fU, gV.edge_label(vv, fU) + capacity)

        # Recursion on each side
        gU_tree = gU._gomory_hu_tree(vertices & frozenset(gU), algorithm=algorithm)
        gV_tree = gV._gomory_hu_tree(vertices & frozenset(gV), algorithm=algorithm)

        # Union of the two partial trees
        g = gU_tree.union(gV_tree)

        # An edge to connect them, with the appropriate label
        g.add_edge(u, v, flow)

        return g

    @doc_index("Connectivity, orientations, trees")
    @rename_keyword(deprecation=19550, method='algorithm')
    def gomory_hu_tree(self, algorithm=None):
        r"""
        Returns a Gomory-Hu tree of self.

        Given a tree `T` with labeled edges representing capacities, it is very
        easy to determine the maximum flow between any pair of vertices :
        it is the minimal label on the edges of the unique path between them.

        Given a graph `G`, a Gomory-Hu tree `T` of `G` is a tree
        with the same set of vertices, and such that the maximum flow
        between any two vertices is the same in `G` as in `T`. See the
        `Wikipedia article on Gomory-Hu tree <http://en.wikipedia.org/wiki/Gomory%E2%80%93Hu_tree>`_.
        Note that, in general, a graph admits more than one Gomory-Hu tree.

        See also 15.4 (Gomory-Hu trees) from [SchrijverCombOpt]_.

        INPUT:

        - ``algorithm`` -- select the algorithm used by the :meth:`edge_cut`
          method. Refer to its documentation for allowed values and default
          behaviour.

        OUTPUT:

        A graph with labeled edges

        EXAMPLE:

        Taking the Petersen graph::

            sage: g = graphs.PetersenGraph()
            sage: t = g.gomory_hu_tree()

        Obviously, this graph is a tree::

            sage: t.is_tree()
            True

        Note that if the original graph is not connected, then the
        Gomory-Hu tree is in fact a forest::

            sage: (2*g).gomory_hu_tree().is_forest()
            True
            sage: (2*g).gomory_hu_tree().is_connected()
            False

        On the other hand, such a tree has lost nothing of the initial
        graph connectedness::

            sage: all([ t.flow(u,v) == g.flow(u,v) for u,v in Subsets( g.vertices(), 2 ) ])
            True

        Just to make sure, we can check that the same is true for two vertices
        in a random graph::

            sage: g = graphs.RandomGNP(20,.3)
            sage: t = g.gomory_hu_tree()
            sage: g.flow(0,1) == t.flow(0,1)
            True

        And also the min cut::

            sage: g.edge_connectivity() == min(t.edge_labels())
            True

        TESTS:

        :trac:`16475`::

            sage: G = graphs.PetersenGraph()
            sage: for u,v in G.edge_iterator(labels=False):
            ....:     G.set_edge_label(u, v, 1)
            sage: for u, v in [(0, 1), (0, 4), (0, 5), (1, 2), (1, 6), (3, 4), (5, 7), (5, 8)]:
            ....:     G.set_edge_label(u, v, 2)
            sage: T = G.gomory_hu_tree()
            sage: from itertools import combinations
            sage: for u,v in combinations(G,2):
            ....:     assert T.flow(u,v,use_edge_labels=True) == G.flow(u,v,use_edge_labels=True)
        """
        if not self.is_connected():
            g = Graph()
            for cc in self.connected_components_subgraphs():
                g = g.union(cc._gomory_hu_tree(frozenset(cc.vertices()), algorithm=algorithm))
        else:
            g = self._gomory_hu_tree(frozenset(self.vertices()), algorithm=algorithm)

        if self.get_pos() is not None:
            g.set_pos(dict(self.get_pos()))
        return g

    @doc_index("Leftovers")
    def two_factor_petersen(self):
        r"""
        Returns a decomposition of the graph into 2-factors.

        Petersen's 2-factor decomposition theorem asserts that any
        `2r`-regular graph `G` can be decomposed into 2-factors.
        Equivalently, it means that the edges of any `2r`-regular
        graphs can be partitionned in `r` sets `C_1,\dots,C_r` such
        that for all `i`, the set `C_i` is a disjoint union of cycles
        ( a 2-regular graph ).

        As any graph of maximal degree `\Delta` can be completed into
        a regular graph of degree `2\lceil\frac\Delta 2\rceil`, this
        result also means that the edges of any graph of degree `\Delta`
        can be partitionned in `r=2\lceil\frac\Delta 2\rceil` sets
        `C_1,\dots,C_r` such that for all `i`, the set `C_i` is a
        graph of maximal degree `2` ( a disjoint union of paths
        and cycles ).

        EXAMPLE:

        The Complete Graph on `7` vertices is a `6`-regular graph, so it can
        be edge-partitionned into `2`-regular graphs::

            sage: g = graphs.CompleteGraph(7)
            sage: classes = g.two_factor_petersen()
            sage: for c in classes:
            ....:     gg = Graph()
            ....:     gg.add_edges(c)
            ....:     print(max(gg.degree())<=2)
            True
            True
            True
            sage: Set(set(classes[0]) | set(classes[1]) | set(classes[2])).cardinality() == g.size()
            True

        ::

            sage: g = graphs.CirculantGraph(24, [7, 11])
            sage: cl = g.two_factor_petersen()
            sage: g.plot(edge_colors={'black':cl[0], 'red':cl[1]})
            Graphics object consisting of 73 graphics primitives

        """
        self._scream_if_not_simple()
        d = self.eulerian_orientation()

        # This new graph is bipartite, and built the following way :
        #
        # To each vertex v of the digraph are associated two vertices,
        # a sink (-1,v) and a source (1,v)
        # Any edge (u,v) in the digraph is then added as ((-1,u),(1,v))

        from sage.graphs.graph import Graph
        g = Graph()
        g.add_edges([((-1,u),(1,v)) for (u,v) in d.edge_iterator(labels=None)])

        # This new bipartite graph is now edge_colored
        from sage.graphs.graph_coloring import edge_coloring
        classes = edge_coloring(g)

        # The edges in the classes are of the form ((-1,u),(1,v))
        # and have to be translated back to (u,v)
        classes_b = []
        for c in classes:
            classes_b.append([(u,v) for ((uu,u),(vv,v)) in c])

        return classes_b

    @doc_index("Leftovers")
    def kirchhoff_symanzik_polynomial(self, name='t'):
        """
        Return the Kirchhoff-Symanzik polynomial of a graph.

        This is a polynomial in variables `t_e` (each of them representing an
        edge of the graph `G`) defined as a sum over all spanning trees:

        .. MATH::

            \Psi_G(t) = \sum_{T\subseteq V\\atop{\\text{a spanning tree}}} \prod_{e \\not\in E(T)} t_e

        This is also called the first Symanzik polynomial or the Kirchhoff
        polynomial.

        INPUT:

        - ``name``: name of the variables (default: ``'t'``)

        OUTPUT:

        - a polynomial with integer coefficients

        ALGORITHM:

            This is computed here using a determinant, as explained in Section
            3.1 of [Marcolli2009]_.

            As an intermediate step, one computes a cycle basis `\mathcal C` of
            `G` and a rectangular `|\mathcal C| \\times |E(G)|` matrix with
            entries in `\{-1,0,1\}`, which describes which edge belong to which
            cycle of `\mathcal C` and their respective orientations.

            More precisely, after fixing an arbitrary orientation for each edge
            `e\in E(G)` and each cycle `C\in\mathcal C`, one gets a sign for
            every incident pair (edge, cycle) which is `1` if the orientation
            coincide and `-1` otherwise.

        EXAMPLES:

        For the cycle of length 5::

            sage: G = graphs.CycleGraph(5)
            sage: G.kirchhoff_symanzik_polynomial()
            t0 + t1 + t2 + t3 + t4

        One can use another letter for variables::

            sage: G.kirchhoff_symanzik_polynomial(name='u')
            u0 + u1 + u2 + u3 + u4

        For the 'coffee bean' graph::

            sage: G = Graph([(0,1,'a'),(0,1,'b'),(0,1,'c')],multiedges=True)
            sage: G.kirchhoff_symanzik_polynomial()
            t0*t1 + t0*t2 + t1*t2

        For the 'parachute' graph::

            sage: G = Graph([(0,2,'a'),(0,2,'b'),(0,1,'c'),(1,2,'d')], multiedges=True)
            sage: G.kirchhoff_symanzik_polynomial()
            t0*t1 + t0*t2 + t1*t2 + t1*t3 + t2*t3

        For the complete graph with 4 vertices::

            sage: G = graphs.CompleteGraph(4)
            sage: G.kirchhoff_symanzik_polynomial()
            t0*t1*t3 + t0*t2*t3 + t1*t2*t3 + t0*t1*t4 + t0*t2*t4 + t1*t2*t4
            + t1*t3*t4 + t2*t3*t4 + t0*t1*t5 + t0*t2*t5 + t1*t2*t5 + t0*t3*t5
            + t2*t3*t5 + t0*t4*t5 + t1*t4*t5 + t3*t4*t5

        REFERENCES:

        .. [Marcolli2009] Matilde Marcolli, Feynman Motives, Chapter 3,
           Feynman integrals and algebraic varieties,
           http://www.its.caltech.edu/~matilde/LectureN3.pdf

        .. [Brown2011] Francis Brown, Multiple zeta values and periods: From
           moduli spaces to Feynman integrals, in Contemporary Mathematics vol
           539
        """
        from sage.matrix.constructor import matrix
        from sage.rings.integer_ring import ZZ
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

        edges = self.edges()
        cycles = self.cycle_basis(output='edge')

        edge2int = {e: j for j, e in enumerate(edges)}
        circuit_mtrx = matrix(ZZ, self.size(), len(cycles))
        for i, cycle in enumerate(cycles):
            for edge in cycle:
                if edge in edges:
                    circuit_mtrx[edge2int[edge], i] = +1
                else:
                    circuit_mtrx[edge2int[(edge[1], edge[0], edge[2])], i] = -1

        D = matrix.diagonal(PolynomialRing(ZZ, name, self.size()).gens())
        return (circuit_mtrx.transpose() * D * circuit_mtrx).determinant()

    @doc_index("Leftovers")
    def ihara_zeta_function_inverse(self):
        """
        Compute the inverse of the Ihara zeta function of the graph.

        This is a polynomial in one variable with integer coefficients. The
        Ihara zeta function itself is the inverse of this polynomial.

        See :wikipedia:`Ihara zeta function`.

        ALGORITHM:

        This is computed here as the (reversed) characteristic
        polynomial of a square matrix of size twice the number of edges,
        related to the adjacency matrix of the line graph, see for example
        Proposition 9 in [ScottStorm]_ and Def. 4.1 in [Terras]_.

        The graph is first replaced by its 2-core, as this does not change
        the Ihara zeta function.

        EXAMPLES::

            sage: G = graphs.CompleteGraph(4)
            sage: factor(G.ihara_zeta_function_inverse())
            (2*t - 1) * (t + 1)^2 * (t - 1)^3 * (2*t^2 + t + 1)^3

            sage: G = graphs.CompleteGraph(5)
            sage: factor(G.ihara_zeta_function_inverse())
            (-1) * (3*t - 1) * (t + 1)^5 * (t - 1)^6 * (3*t^2 + t + 1)^4

            sage: G = graphs.PetersenGraph()
            sage: factor(G.ihara_zeta_function_inverse())
            (-1) * (2*t - 1) * (t + 1)^5 * (t - 1)^6 * (2*t^2 + 2*t + 1)^4
            * (2*t^2 - t + 1)^5

            sage: G = graphs.RandomTree(10)
            sage: G.ihara_zeta_function_inverse()
            1

        REFERENCES:

        .. [HST] Matthew D. Horton, H. M. Stark, and Audrey A. Terras,
           What are zeta functions of graphs and what are they good for?
           in Quantum graphs and their applications, 173-189,
           Contemp. Math., Vol. 415

        .. [Terras] Audrey Terras, Zeta functions of graphs: a stroll through
           the garden, Cambridge Studies in Advanced Mathematics, Vol. 128

        .. [ScottStorm] Geoffrey Scott and Christopher Storm, The coefficients
           of the Ihara zeta function, Involve (http://msp.org/involve/2008/1-2/involve-v1-n2-p08-p.pdf)
        """
        from sage.matrix.constructor import matrix

        H = self.subgraph(vertices=self.cores(k=2)[1])
        E = H.edges()
        m = len(E)
        # compute (Hashimoto) edge matrix T
        T = matrix(ZZ, 2 * m, 2 * m, 0)
        for i in range(m):
            for j in range(m):
                if i != j:
                    if E[i][1] == E[j][0]:  # same orientation
                        T[2 * i, 2 * j] = 1
                        T[2 * j + 1, 2 * i + 1] = 1
                    elif E[i][1] == E[j][1]:  # opposite orientation (towards)
                        T[2 * i, 2 * j + 1] = 1
                        T[2 * j, 2 * i + 1] = 1
                    elif E[i][0] == E[j][0]:  # opposite orientation (away)
                        T[2 * i + 1, 2 * j] = 1
                        T[2 * j + 1, 2 * i] = 1
        return T.charpoly('t').reverse()

    @doc_index("Leftovers")
    def perfect_matchings(self, labels=False):
        """
        Return an interator over all perfect matchings of the graph.

        ALGORITHM:

        Choose a vertex `v`, then recurse through all edges incident to `v`,
        removing one edge at a time whenever an edge is added to a matching.

        INPUT:

        - ``labels`` -- boolean (default: ``False``); when ``True``, the
          edges in each perfect matching are triples (containing the label
          as the third element), otherwise the edges are pairs

        .. SEEALSO::

            :meth:`matching`

        EXAMPLES::

            sage: G=graphs.GridGraph([2,3])
            sage: list(G.perfect_matchings())
            [[((0, 0), (0, 1)), ((0, 2), (1, 2)), ((1, 0), (1, 1))],
             [((0, 1), (0, 2)), ((1, 1), (1, 2)), ((0, 0), (1, 0))],
             [((0, 1), (1, 1)), ((0, 2), (1, 2)), ((0, 0), (1, 0))]]

            sage: G = graphs.CompleteGraph(4)
            sage: list(G.perfect_matchings(labels=True))
            [[(0, 1, None), (2, 3, None)],
             [(0, 2, None), (1, 3, None)],
             [(0, 3, None), (1, 2, None)]]

            sage: G = Graph([[1,-1,'a'], [2,-2, 'b'], [1,-2,'x'], [2,-1,'y']])
            sage: list(G.perfect_matchings(labels=True))
            [[(-2, 1, 'x'), (-1, 2, 'y')],
             [(-1, 1, 'a'), (-2, 2, 'b')]]

            sage: G = graphs.CompleteGraph(8)
            sage: mpc = G.matching_polynomial().coefficients(sparse=False)[0]
            sage: len(list(G.perfect_matchings())) == mpc
            True

            sage: G = graphs.PetersenGraph().copy(immutable=True)
            sage: list(G.perfect_matchings())
            [[(0, 1), (2, 3), (4, 9), (6, 8), (5, 7)],
             [(0, 1), (2, 7), (3, 4), (5, 8), (6, 9)],
             [(0, 4), (1, 2), (3, 8), (6, 9), (5, 7)],
             [(0, 4), (1, 6), (2, 3), (5, 8), (7, 9)],
             [(0, 5), (1, 2), (3, 4), (6, 8), (7, 9)],
             [(0, 5), (1, 6), (2, 7), (3, 8), (4, 9)]]

            sage: list(Graph().perfect_matchings())
            [[]]

            sage: G = graphs.CompleteGraph(5)
            sage: list(G.perfect_matchings())
            []
        """
        if not self:
            yield []
            return
        # if every connected component has an even number of vertices
        if all(len(cc) % 2 == 0 for cc in self.connected_components()):
            v = next(self.vertex_iterator())
            for e in self.edges_incident(v, labels=labels):
                Gp = self.copy(immutable=False)
                Gp.delete_vertices([e[0], e[1]])
                for mat in Gp.perfect_matchings(labels):
                    yield [e] + mat


# Aliases to functions defined in Cython modules
import types

import sage.graphs.weakly_chordal
Graph.is_long_hole_free         = types.MethodType(sage.graphs.weakly_chordal.is_long_hole_free, None, Graph)
Graph.is_long_antihole_free     = types.MethodType(sage.graphs.weakly_chordal.is_long_antihole_free, None, Graph)
Graph.is_weakly_chordal         = types.MethodType(sage.graphs.weakly_chordal.is_weakly_chordal, None, Graph)

import sage.graphs.asteroidal_triples
Graph.is_asteroidal_triple_free = types.MethodType(sage.graphs.asteroidal_triples.is_asteroidal_triple_free, None, Graph)

import sage.graphs.chrompoly
Graph.chromatic_polynomial      = types.MethodType(sage.graphs.chrompoly.chromatic_polynomial, None, Graph)

import sage.graphs.graph_decompositions.rankwidth
Graph.rank_decomposition        = types.MethodType(sage.graphs.graph_decompositions.rankwidth.rank_decomposition, None, Graph)

import sage.graphs.matchpoly
Graph.matching_polynomial       = types.MethodType(sage.graphs.matchpoly.matching_polynomial, None, Graph)

import sage.graphs.cliquer
Graph.cliques_maximum           = types.MethodType(sage.graphs.cliquer.all_max_clique, None, Graph)

import sage.graphs.spanning_tree
Graph.random_spanning_tree      = types.MethodType(sage.graphs.spanning_tree.random_spanning_tree, None, Graph)

import sage.graphs.graph_decompositions.graph_products
Graph.is_cartesian_product      = types.MethodType(sage.graphs.graph_decompositions.graph_products.is_cartesian_product, None, Graph)

import sage.graphs.distances_all_pairs
Graph.is_distance_regular       = types.MethodType(sage.graphs.distances_all_pairs.is_distance_regular, None, Graph)

import sage.graphs.base.static_dense_graph
Graph.is_strongly_regular       = types.MethodType(sage.graphs.base.static_dense_graph.is_strongly_regular, None, Graph)

# From Python modules
import sage.graphs.line_graph
Graph.is_line_graph             = sage.graphs.line_graph.is_line_graph

from sage.graphs.tutte_polynomial import tutte_polynomial
Graph.tutte_polynomial          = tutte_polynomial

from sage.graphs.lovasz_theta import lovasz_theta
Graph.lovasz_theta              = lovasz_theta

from sage.graphs.partial_cube import is_partial_cube
Graph.is_partial_cube           = is_partial_cube

_additional_categories = {
    Graph.is_long_hole_free         : "Graph properties",
    Graph.is_long_antihole_free     : "Graph properties",
    Graph.is_weakly_chordal         : "Graph properties",
    Graph.is_asteroidal_triple_free : "Graph properties",
    Graph.chromatic_polynomial      : "Algorithmically hard stuff",
    Graph.rank_decomposition        : "Algorithmically hard stuff",
    Graph.matching_polynomial       : "Algorithmically hard stuff",
    Graph.cliques_maximum           : "Clique-related methods",
    Graph.random_spanning_tree      : "Connectivity, orientations, trees",
    Graph.is_cartesian_product      : "Graph properties",
    Graph.is_distance_regular       : "Graph properties",
    Graph.is_strongly_regular       : "Graph properties",
    Graph.is_line_graph             : "Graph properties",
    Graph.is_partial_cube           : "Graph properties",
    Graph.tutte_polynomial          : "Algorithmically hard stuff",
    Graph.lovasz_theta              : "Leftovers",
    }

__doc__ = __doc__.replace("{INDEX_OF_METHODS}",gen_thematic_rest_table_index(Graph,_additional_categories))
