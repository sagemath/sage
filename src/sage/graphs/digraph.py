# sage_setup: distribution = sagemath-graphs
r"""
Directed graphs

This module implements functions and operations involving directed
graphs. Here is what they can do

**Graph basic operations:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~DiGraph.layout_acyclic_dummy` | Compute a (dummy) ranked layout so that all edges point upward.
    :meth:`~DiGraph.layout_acyclic` | Compute a ranked layout so that all edges point upward.
    :meth:`~DiGraph.reverse` | Return a copy of digraph with edges reversed in direction.
    :meth:`~DiGraph.reverse_edge` | Reverse an edge.
    :meth:`~DiGraph.reverse_edges` | Reverse a list of edges.
    :meth:`~DiGraph.out_degree_sequence` | Return the outdegree sequence.
    :meth:`~DiGraph.out_degree_iterator` | Same as degree_iterator, but for out degree.
    :meth:`~DiGraph.out_degree` | Same as degree, but for out degree.
    :meth:`~DiGraph.in_degree_sequence` | Return the indegree sequence of this digraph.
    :meth:`~DiGraph.in_degree_iterator` | Same as degree_iterator, but for in degree.
    :meth:`~DiGraph.in_degree` | Same as degree, but for in-degree.
    :meth:`~DiGraph.neighbors_out` | Return the list of the out-neighbors of a given vertex.
    :meth:`~DiGraph.neighbor_out_iterator` | Return an iterator over the out-neighbors of a given vertex.
    :meth:`~DiGraph.neighbors_in` | Return the list of the in-neighbors of a given vertex.
    :meth:`~DiGraph.neighbor_in_iterator` | Return an iterator over the in-neighbors of vertex.
    :meth:`~DiGraph.outgoing_edges` | Return a list of edges departing from vertices.
    :meth:`~DiGraph.outgoing_edge_iterator` | Return an iterator over all departing edges from vertices
    :meth:`~DiGraph.incoming_edges` | Return a list of edges arriving at vertices.
    :meth:`~DiGraph.incoming_edge_iterator` | Return an iterator over all arriving edges from vertices
    :meth:`~DiGraph.sources` | Return the list of all sources (vertices without incoming edges) of this digraph.
    :meth:`~DiGraph.sinks` | Return the list of all sinks (vertices without outgoing edges) of this digraph.
    :meth:`~DiGraph.to_undirected` | Return an undirected version of the graph.
    :meth:`~DiGraph.to_directed` | Since the graph is already directed, simply returns a copy of itself.
    :meth:`~DiGraph.is_directed` | Since digraph is directed, returns True.
    :meth:`~DiGraph.dig6_string` | Return the ``dig6`` representation of the digraph as an ASCII string.

**Distances:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~DiGraph.eccentricity` | Return the eccentricity of vertex (or vertices) ``v``.
    :meth:`~DiGraph.radius` | Return the radius of the DiGraph.
    :meth:`~DiGraph.diameter` | Return the diameter of the DiGraph.
    :meth:`~DiGraph.center` | Return the set of vertices in the center of the DiGraph.
    :meth:`~DiGraph.periphery` | Return the set of vertices in the periphery of the DiGraph.

**Paths and cycles:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~DiGraph.all_cycles_iterator` | Return an iterator over all the cycles of ``self`` starting with one of the given vertices.
    :meth:`~DiGraph.all_simple_cycles` | Return a list of all simple cycles of ``self``.

**Representation theory:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~DiGraph.path_semigroup` | Return the (partial) semigroup formed by the paths of the digraph.
    :meth:`~DiGraph.auslander_reiten_quiver` | Return the Auslander-Reiten quiver of ``self``.

**Connectivity:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~DiGraph.is_strongly_connected` | Check whether the current ``DiGraph`` is strongly connected.
    :meth:`~DiGraph.strongly_connected_components_digraph` | Return the digraph of the strongly connected components
    :meth:`~DiGraph.strongly_connected_components_subgraphs` | Return the strongly connected components as a list of subgraphs.
    :meth:`~DiGraph.strongly_connected_component_containing_vertex` | Return the strongly connected component containing a given vertex
    :meth:`~DiGraph.strongly_connected_components` | Return the list of strongly connected components.
    :meth:`~DiGraph.strong_articulation_points` | Return the strong articulation points of this digraph.


**Acyclicity:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~DiGraph.is_directed_acyclic` | Check whether the digraph is acyclic or not.
    :meth:`~DiGraph.is_transitive` | Check whether the digraph is transitive or not.
    :meth:`~DiGraph.is_aperiodic` | Check whether the digraph is aperiodic or not.
    :meth:`~DiGraph.is_tournament` | Check whether the digraph is a tournament.
    :meth:`~DiGraph.period` | Return the period of the digraph.
    :meth:`~DiGraph.level_sets` | Return the level set decomposition of the digraph.
    :meth:`~DiGraph.topological_sort_generator` | Return a list of all topological sorts of the digraph if it is acyclic
    :meth:`~DiGraph.topological_sort` | Return a topological sort of the digraph if it is acyclic

**Hard stuff:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~DiGraph.feedback_edge_set` | Compute the minimum feedback edge (arc) set of a digraph

**Miscellaneous:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~DiGraph.flow_polytope` | Compute the flow polytope of a digraph
    :meth:`~DiGraph.degree_polynomial` | Return the generating polynomial of degrees of vertices in ``self``.
    :meth:`~DiGraph.out_branchings` | Return an iterator over the out branchings rooted at given vertex in ``self``.
    :meth:`~DiGraph.in_branchings` | Return an iterator over the in branchings rooted at given vertex in ``self``.

Methods
-------
"""

# ****************************************************************************
#       Copyright (C) 2010      Alexandre Blondin Masse <alexandre.blondin.masse at gmail.com>
#                               Carl Witty <cwitty@newtonlabs.com>
#                               Gregory McWhirter <gmcwhirt@uci.edu>
#                               Minh Van Nguyen <nguyenminh2@gmail.com>
#                     2010-2011 Robert L. Miller <rlm@rlmiller.org>
#                     2010-2015 Nathann Cohen <nathann.cohen@gmail.com>
#                               Nicolas M. Thiery <nthiery@users.sf.net>
#                     2011      Johannes Klaus Fichte <fichte@kr.tuwien.ac.at>
#                     2012      Javier López Peña <vengoroso@gmail.com>
#                     2012      Jim Stark <jstarx@gmail.com>
#                     2012      Karl-Dieter Crisman <kcrisman@gmail.com>
#                     2012      Keshav Kini <keshav.kini@gmail.com>
#                     2012      Lukas Lansky <lansky@kam.mff.cuni.cz>
#                     2012-2015 Volker Braun <vbraun.name@gmail.com>
#                     2012-2017 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#                     2012-2018 David Coudert <david.coudert@inria.fr>
#                     2013      Emily Gunawan <egunawan@umn.edu>
#                     2013      Gregg Musiker <musiker@math.mit.edu>
#                     2013      Mathieu Guay-Paquet <mathieu.guaypaquet@gmail.com>
#                     2013-2014 Simon King <simon.king@uni-jena.de>
#                     2014      Clemens Heuberger <clemens.heuberger@aau.at>
#                               Erik Massop <e.massop@hccnet.nl>
#                               R. Andrew Ohana <andrew.ohana@gmail.com>
#                               Wilfried Luebbe <wluebbe@gmail.com>
#                     2014-2015 André Apitzsch <andre.apitzsch@etit.tu-chemnitz.de>
#                               Darij Grinberg <darijgrinberg@gmail.com>
#                               Travis Scrimshaw <tscrim at ucdavis.edu>
#                               Vincent Delecroix <20100.delecroix@gmail.com>
#                     2014-2017 Frédéric Chapoton <chapoton@math.univ-lyon1.fr>
#                     2015      Michele Borassi <michele.borassi@imtlucca.it>
#                     2015-2017 John H. Palmieri <palmieri@math.washington.edu>
#                               Jori Mäntysalo <jori.mantysalo@uta.fi>
#                     2016      Dima Pasechnik <dimpase@gmail.com>
#                     2018      Meghana M Reddy <mreddymeghana@gmail.com>
#                               Julian Rüth <julian.rueth@fsfe.org>
#                     2019      Rajat Mittal <rajat.mttl@gmail.com>
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from copy import copy
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from itertools import product
import sage.graphs.generic_graph_pyx as generic_graph_pyx
from sage.graphs.generic_graph import GenericGraph
from sage.graphs.dot2tex_utils import have_dot2tex
from sage.graphs.views import EdgesView


class DiGraph(GenericGraph):
    r"""
    Directed graph.

    A digraph or directed graph is a set of vertices connected by oriented
    edges. See also the :wikipedia:`Directed_graph`. For a collection of
    pre-defined digraphs, see the :mod:`~sage.graphs.digraph_generators` module.

    A :class:`DiGraph` object has many methods whose list can be obtained by
    typing ``g.<tab>`` (i.e. hit the :kbd:`Tab` key) or by reading the documentation
    of :mod:`~sage.graphs.digraph`, :mod:`~sage.graphs.generic_graph`, and
    :mod:`~sage.graphs.graph`.

    INPUT:

    By default, a :class:`DiGraph` object is simple (i.e. no *loops* nor
    *multiple edges*) and unweighted. This can be easily tuned with the
    appropriate flags (see below).

    - ``data`` -- can be any of the following (see the ``format`` argument):

      #. ``DiGraph()`` -- build a digraph on 0 vertices

      #. ``DiGraph(5)`` -- return an edgeless digraph on the 5 vertices 0,...,4

      #. ``DiGraph([list_of_vertices, list_of_edges])`` -- return a digraph with
         given vertices/edges

         To bypass auto-detection, prefer the more explicit
         ``DiGraph([V, E], format='vertices_and_edges')``.

      #. ``DiGraph(list_of_edges)`` -- return a digraph with a given list of
         edges (see documentation of
         :meth:`~sage.graphs.generic_graph.GenericGraph.add_edges`).

         To bypass auto-detection, prefer the more explicit
         ``DiGraph(L, format='list_of_edges')``.

      #. ``DiGraph({1: [2,3,4], 3: [4]})`` -- return a digraph by associating to
         each vertex the list of its out-neighbors.

         To bypass auto-detection, prefer the more explicit
         ``DiGraph(D, format='dict_of_lists')``.

      #. ``DiGraph({1: {2: 'a', 3: 'b'}, 3: {2: 'c'}})`` -- return a digraph by
         associating a list of out-neighbors to each vertex and providing its
         edge label.

         To bypass auto-detection, prefer the more explicit
         ``DiGraph(D, format='dict_of_dicts')``.

         For digraphs with multiple edges, you can provide a list of labels
         instead, e.g.: ``DiGraph({1: {2: ['a1', 'a2'], 3:['b']},
         3:{2:['c']}})``.

      #. ``DiGraph(a_matrix)`` -- return a digraph with given (weighted)
         adjacency matrix (see documentation of
         :meth:`~sage.graphs.generic_graph.GenericGraph.adjacency_matrix`).

         To bypass auto-detection, prefer the more explicit ``DiGraph(M,
         format='adjacency_matrix')``. To take weights into account, use
         ``format='weighted_adjacency_matrix'`` instead.

      #. ``DiGraph(a_nonsquare_matrix)`` -- return a digraph with given
         incidence matrix (see documentation of
         :meth:`~sage.graphs.generic_graph.GenericGraph.incidence_matrix`).

         To bypass auto-detection, prefer the more explicit ``DiGraph(M,
         format='incidence_matrix')``.

      #. ``DiGraph([V, f])`` -- return a digraph with a vertex set ``V`` and an
         edge `u,v` whenever `f(u, v)` is ``True``. Example: ``DiGraph([
         [1..10], lambda x,y: abs(x - y).is_square()])``

      #. ``DiGraph('FOC@?OC@_?')`` -- return a digraph from a ``dig6`` string
         (see documentation of :meth:`~dig6_string`).

      #. ``DiGraph(another_digraph)`` -- return a digraph from a Sage (di)graph,
         `pygraphviz <https://pygraphviz.github.io/>`__ digraph, `NetworkX
         <https://networkx.github.io/>`__ digraph, or `igraph
         <http://igraph.org/python/>`__ digraph.

    - ``pos`` -- dictionary (default: ``None``); a positioning dictionary. For
      example, the spring layout from NetworkX for the 5-cycle is::

         {0: [-0.91679746, 0.88169588],
          1: [ 0.47294849, 1.125     ],
          2: [ 1.125     ,-0.12867615],
          3: [ 0.12743933,-1.125     ],
          4: [-1.125     ,-0.50118505]}

    - ``name`` -- string (default: ``None``); gives the graph a name (e.g.,
      name='complete')

    - ``loops`` -- boolean (default: ``None``); whether to allow loops (ignored
      if data is an instance of the DiGraph class)

    - ``multiedges`` -- boolean (default: ``None``); whether to allow multiple
      edges (ignored if data is an instance of the DiGraph class)

    - ``weighted`` -- boolean (default: ``None``); whether digraph thinks of
      itself as weighted or not. See ``self.weighted()``

    - ``format`` -- string (default: ``None``); if set to ``None``,
      :class:`DiGraph` tries to guess input's format. To avoid this possibly
      time-consuming step, one of the following values can be specified (see
      description above): ``'int'``, ``'dig6'``, ``'rule'``,
      ``'list_of_edges'``, ``'dict_of_lists'``, ``'dict_of_dicts'``,
      ``'adjacency_matrix'``, ``'weighted_adjacency_matrix'``,
      ``'incidence_matrix'``, ``"NX"``, ``'igraph'``.

    - ``sparse`` -- boolean (default: ``True``); ``sparse=True`` is an alias for
      ``data_structure="sparse"``, and ``sparse=False`` is an alias for
      ``data_structure="dense"``

    - ``data_structure`` -- string (default: ``'sparse'``); one of the following
      (for more information, see :mod:`~sage.graphs.base.overview`):

      * ``'dense'`` -- selects the :mod:`~sage.graphs.base.dense_graph` backend

      * ``'sparse'`` -- selects the :mod:`~sage.graphs.base.sparse_graph`
        backend

      * ``'static_sparse'`` -- selects the
        :mod:`~sage.graphs.base.static_sparse_backend` (this backend is faster
        than the sparse backend and smaller in memory, and it is immutable, so
        that the resulting graphs can be used as dictionary keys).

    - ``immutable`` -- boolean (default: ``False``); whether to create a
      immutable digraph. Note that ``immutable=True`` is actually a shortcut for
      ``data_structure='static_sparse'``.

    - ``hash_labels`` -- boolean (default: ``None``); whether to include edge
      labels during hashing. This parameter defaults to ``True`` if the digraph
      is weighted. This parameter is ignored if the digraph is mutable.
      Beware that trying to hash unhashable labels will raise an error.

    - ``vertex_labels`` -- boolean (default: ``True``); whether to allow any
      object as a vertex (slower), or only the integers `0,...,n-1`, where `n`
      is the number of vertices.

    - ``convert_empty_dict_labels_to_None`` -- boolean (default: ``None``); this
      arguments sets the default edge labels used by NetworkX (empty
      dictionaries) to be replaced by ``None``, the default Sage edge label. It
      is set to ``True`` iff a NetworkX graph is on the input.

    EXAMPLES:

    #. A dictionary of dictionaries::

            sage: g = DiGraph({0: {1: 'x', 2: 'z', 3: 'a'}, 2: {5: 'out'}}); g
            Digraph on 5 vertices

       The labels ('x', 'z', 'a', 'out') are labels for edges. For example,
       'out' is the label for the edge from 2 to 5. Labels can be used as
       weights, if all the labels share some common parent.

    #. A dictionary of lists (or iterables)::

            sage: g = DiGraph({0: [1, 2, 3], 2: [4]}); g
            Digraph on 5 vertices
            sage: g = DiGraph({0: (1, 2, 3), 2: (4,)}); g
            Digraph on 5 vertices

    #. A list of vertices and a function describing adjacencies. Note that the
       list of vertices and the function must be enclosed in a list (i.e.,
       ``[list of vertices, function]``).

       We construct a graph on the integers 1 through 12 such that there is a
       directed edge from `i` to `j` if and only if `i` divides `j`::

            sage: g = DiGraph([[1..12], lambda i,j: i != j and i.divides(j)])
            sage: g.vertices(sort=True)
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
            sage: g.adjacency_matrix()                                                  # needs sage.modules
            [0 1 1 1 1 1 1 1 1 1 1 1]
            [0 0 0 1 0 1 0 1 0 1 0 1]
            [0 0 0 0 0 1 0 0 1 0 0 1]
            [0 0 0 0 0 0 0 1 0 0 0 1]
            [0 0 0 0 0 0 0 0 0 1 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 1]
            [0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0]

    #. A Sage matrix: Note: If format is not specified, then Sage assumes a
       square matrix is an adjacency matrix, and a nonsquare matrix is an
       incidence matrix.

       - an adjacency matrix::

            sage: M = Matrix([[0, 1, 1, 1, 0], [0, 0, 0, 0, 0],                         # needs sage.modules
            ....:             [0, 0, 0, 0, 1], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]]); M
            [0 1 1 1 0]
            [0 0 0 0 0]
            [0 0 0 0 1]
            [0 0 0 0 0]
            [0 0 0 0 0]
            sage: DiGraph(M)                                                            # needs sage.modules
            Digraph on 5 vertices

            sage: M = Matrix([[0,1,-1], [-1,0,-1/2], [1,1/2,0]]); M                     # needs sage.modules
            [   0    1   -1]
            [  -1    0 -1/2]
            [   1  1/2    0]
            sage: G = DiGraph(M, sparse=True, weighted=True); G                         # needs sage.modules
            Digraph on 3 vertices
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
            sage: DiGraph(M)                                                            # needs sage.modules
            Digraph on 6 vertices

    #. A ``dig6`` string: Sage automatically recognizes whether a string is in
       ``dig6`` format, which is a directed version of ``graph6``::

            sage: D = DiGraph('IRAaDCIIOWEOKcPWAo')
            sage: D
            Digraph on 10 vertices

            sage: D = DiGraph('IRAaDCIIOEOKcPWAo')
            Traceback (most recent call last):
            ...
            RuntimeError: the string (IRAaDCIIOEOKcPWAo) seems corrupt: for n = 10, the string is too short

            sage: D = DiGraph("IRAaDCI'OWEOKcPWAo")
            Traceback (most recent call last):
            ...
            RuntimeError: the string seems corrupt: valid characters are
            ?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~

    #. A NetworkX MultiDiGraph::

            sage: import networkx                                                       # needs networkx
            sage: g = networkx.MultiDiGraph({0: [1, 2, 3], 2: [4]})                     # needs networkx
            sage: DiGraph(g)                                                            # needs networkx
            Multi-digraph on 5 vertices


    #. A NetworkX digraph::

            sage: import networkx                                                       # needs networkx
            sage: g = networkx.DiGraph({0: [1, 2, 3], 2: [4]})                          # needs networkx
            sage: DiGraph(g)                                                            # needs networkx
            Digraph on 5 vertices

    #. An igraph directed Graph (see also
       :meth:`~sage.graphs.generic_graph.GenericGraph.igraph_graph`)::

           sage: import igraph                                   # optional - python_igraph
           sage: g = igraph.Graph([(0,1),(0,2)], directed=True)  # optional - python_igraph
           sage: DiGraph(g)                                      # optional - python_igraph
           Digraph on 3 vertices

       If ``vertex_labels`` is ``True``, the names of the vertices are given by
       the vertex attribute ``'name'``, if available::

           sage: # optional - python_igraph
           sage: g = igraph.Graph([(0,1),(0,2)], directed=True, vertex_attrs={'name':['a','b','c']})
           sage: DiGraph(g).vertices(sort=True)
           ['a', 'b', 'c']
           sage: g = igraph.Graph([(0,1),(0,2)], directed=True, vertex_attrs={'label':['a','b','c']})
           sage: DiGraph(g).vertices(sort=True)
           [0, 1, 2]

       If the igraph Graph has edge attributes, they are used as edge labels::

           sage: g = igraph.Graph([(0, 1), (0, 2)], directed=True,                  # optional - python_igraph
           ....:                  edge_attrs={'name':['a', 'b'], 'weight':[1, 3]})
           sage: DiGraph(g).edges(sort=True)                                        # optional - python_igraph
           [(0, 1, {'name': 'a', 'weight': 1}), (0, 2, {'name': 'b', 'weight': 3})]


    TESTS::

        sage: DiGraph({0:[1,2,3], 2:[4]}).edges(sort=True)
        [(0, 1, None), (0, 2, None), (0, 3, None), (2, 4, None)]
        sage: DiGraph({0:(1,2,3), 2:(4,)}).edges(sort=True)
        [(0, 1, None), (0, 2, None), (0, 3, None), (2, 4, None)]
        sage: DiGraph({0:Set([1,2,3]), 2:Set([4])}).edges(sort=True)
        [(0, 1, None), (0, 2, None), (0, 3, None), (2, 4, None)]

    Demonstrate that digraphs using the static backend are equal to mutable
    graphs but can be used as dictionary keys::

        sage: # needs networkx
        sage: import networkx
        sage: g = networkx.DiGraph({0:[1,2,3], 2:[4]})
        sage: G = DiGraph(g)
        sage: G_imm = DiGraph(G, data_structure='static_sparse')
        sage: H_imm = DiGraph(G, data_structure='static_sparse')
        sage: H_imm is G_imm
        False
        sage: H_imm == G_imm == G
        True
        sage: {G_imm:1}[H_imm]
        1
        sage: {G_imm:1}[G]
        Traceback (most recent call last):
        ...
        TypeError: This graph is mutable, and thus not hashable. Create an
        immutable copy by `g.copy(immutable=True)`

    The error message states that one can also create immutable graphs by
    specifying the ``immutable`` optional argument (not only by
    ``data_structure='static_sparse'`` as above)::

        sage: J_imm = DiGraph(G, immutable=True)                                        # needs networkx
        sage: J_imm == G_imm                                                            # needs networkx
        True
        sage: type(J_imm._backend) == type(G_imm._backend)                              # needs networkx
        True

    From a list of vertices and a list of edges::

        sage: G = DiGraph([[1,2,3],[(1,2)]]); G
        Digraph on 3 vertices
        sage: G.edges(sort=True)
        [(1, 2, None)]

    Check that :issue:`27505` is fixed::

        sage: DiGraph(DiGraph().networkx_graph(), weighted=None, format='NX')           # needs networkx
        Digraph on 0 vertices
    """
    _directed = True

    def __init__(self, data=None, pos=None, loops=None, format=None,
                 weighted=None, data_structure='sparse',
                 vertex_labels=True, name=None,
                 multiedges=None, convert_empty_dict_labels_to_None=None,
                 sparse=True, immutable=False, hash_labels=None):
        """
        TESTS::

            sage: D = DiGraph()
            sage: loads(dumps(D)) == D
            True

            sage: a = matrix(2,2,[1,2,0,1])                                             # needs sage.modules
            sage: DiGraph(a, sparse=True).adjacency_matrix() == a                       # needs sage.modules
            True

            sage: a = matrix(2,2,[3,2,0,1])                                             # needs sage.modules
            sage: DiGraph(a, sparse=True).adjacency_matrix() == a                       # needs sage.modules
            True

        The positions are copied when the DiGraph is built from another DiGraph
        or from a Graph ::

            sage: g = DiGraph(graphs.PetersenGraph())
            sage: h = DiGraph(g)
            sage: g.get_pos() == h.get_pos()
            True
            sage: g.get_pos() == graphs.PetersenGraph().get_pos()
            True

        The position dictionary is not the input one (:issue:`22424`)::

            sage: my_pos = {0:(0,0), 1:(1,1)}
            sage: D = DiGraph([[0,1], [(0,1)]], pos=my_pos)
            sage: my_pos == D._pos
            True
            sage: my_pos is D._pos
            False

        Detection of multiple edges::

            sage: DiGraph({1:{2:[0,1]}})
            Multi-digraph on 2 vertices
            sage: DiGraph({1:{2:0}})
            Digraph on 2 vertices

        An empty list or dictionary defines a simple graph (:issue:`10441` and
        :issue:`12910`)::

            sage: DiGraph([])
            Digraph on 0 vertices
            sage: DiGraph({})
            Digraph on 0 vertices
            sage: # not "Multi-digraph on 0 vertices"

        Problem with weighted adjacency matrix (:issue:`13919`)::

            sage: B = {0:{1:2,2:5,3:4},1:{2:2,4:7},2:{3:1,4:4,5:3},
            ....:      3:{5:4},4:{5:1,6:5},5:{4:1,6:7,5:1}}
            sage: grafo3 = DiGraph(B, weighted=True)
            sage: matad = grafo3.weighted_adjacency_matrix()                            # needs sage.modules
            sage: grafo4 = DiGraph(matad, format='adjacency_matrix', weighted=True)     # needs sage.modules
            sage: grafo4.shortest_path(0, 6, by_weight=True)                            # needs sage.modules
            [0, 1, 2, 5, 4, 6]

        Building a DiGraph with ``immutable=False`` returns a mutable graph::

            sage: g = graphs.PetersenGraph()
            sage: g = DiGraph(g.edges(sort=True), immutable=False)
            sage: g.add_edge("Hey", "Heyyyyyyy")
            sage: {g:1}[g]
            Traceback (most recent call last):
            ...
            TypeError: This graph is mutable, and thus not hashable. Create an immutable copy by `g.copy(immutable=True)`
            sage: copy(g) is g
            False
            sage: {g.copy(immutable=True):1}[g.copy(immutable=True)]
            1

        But building it with ``immutable=True`` returns an immutable graph::

            sage: g = DiGraph(graphs.PetersenGraph(), immutable=True)
            sage: g.add_edge("Hey", "Heyyyyyyy")
            Traceback (most recent call last):
            ...
            ValueError: graph is immutable; please change a copy instead (use function copy())
            sage: {g:1}[g]
            1
            sage: copy(g) is g    # copy is mutable again
            False

        Unknown input format::

            sage: DiGraph(4, format='HeyHeyHey')
            Traceback (most recent call last):
            ...
            ValueError: unknown input format 'HeyHeyHey'

        Sage DiGraph from igraph undirected graph::

            sage: import igraph            # optional - python_igraph
            sage: DiGraph(igraph.Graph())  # optional - python_igraph
            Traceback (most recent call last):
            ...
            ValueError: a *directed* igraph graph was expected. To build an undirected graph, call the Graph constructor

        Vertex labels are retained in the graph (:issue:`14708`)::

            sage: g = DiGraph()
            sage: g.add_vertex(0)
            sage: g.set_vertex(0, 'foo')
            sage: g.get_vertices()
            {0: 'foo'}
            sage: DiGraph(g).get_vertices()
            {0: 'foo'}
        """
        msg = ''
        GenericGraph.__init__(self)
        from sage.structure.element import Matrix

        if sparse is False:
            if data_structure != "sparse":
                raise ValueError("the 'sparse' argument is an alias for "
                                 "'data_structure', please do not define both")
            data_structure = "dense"

        if multiedges or weighted:
            if data_structure == "dense":
                raise RuntimeError("multiedge and weighted c_graphs must be sparse")

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
        self._backend = CGB(0, directed=True)

        if format is None and isinstance(data, str):
            format = 'dig6'
            if data[:8] == ">>dig6<<":
                data = data[8:]
        if format is None and isinstance(data, Matrix):
            if data.is_square():
                format = 'adjacency_matrix'
            else:
                format = 'incidence_matrix'
                msg += "Non-symmetric or non-square matrix assumed to be an incidence matrix: "
        if format is None and isinstance(data, DiGraph):
            format = 'DiGraph'
        from sage.graphs.graph import Graph
        if format is None and isinstance(data, Graph):
            data = data.to_directed()
            format = 'DiGraph'
        if format is None and isinstance(data, list) and \
           len(data) >= 2 and callable(data[1]):
            format = 'rule'

        if (format is None and
            isinstance(data, list) and
            len(data) == 2 and
            isinstance(data[0], list) and    # a list of two lists, the second of
            ((isinstance(data[1], list) and  # which contains iterables (the edges)
              (not data[1] or callable(getattr(data[1][0], "__iter__", None)))) or
             (isinstance(data[1], EdgesView)))):
            format = "vertices_and_edges"

        if format is None and isinstance(data, dict):
            if not data:
                format = 'dict_of_dicts'
            else:
                val = next(iter(data.values()))
                if isinstance(val, dict):
                    format = 'dict_of_dicts'
                else:
                    format = 'dict_of_lists'
        if format is None and hasattr(data, 'adj'):
            # the input is a networkx (Multi)(Di)Graph
            format = 'NX'

        if (format is None and hasattr(data, 'vcount') and
                hasattr(data, 'get_edgelist')):
            try:
                import igraph
            except ImportError:
                raise ImportError("the data seems to be a igraph object, but "
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

        if format == 'weighted_adjacency_matrix':
            if weighted is False:
                raise ValueError("format was weighted_adjacency_matrix but weighted was False")
            if weighted is None:
                weighted = True
            if multiedges is None:
                multiedges = False
            format = 'adjacency_matrix'

        if format is None:
            raise ValueError("This input cannot be turned into a graph")

        # At this point, format has been set. We build the graph

        if format == 'dig6':
            if weighted is None:
                self._weighted = False
            self.allow_loops(bool(loops), check=False)
            self.allow_multiple_edges(bool(multiedges), check=False)
            from .graph_input import from_dig6
            from_dig6(self, data)

        elif format == 'adjacency_matrix':
            from .graph_input import from_adjacency_matrix
            from_adjacency_matrix(self, data, loops=loops, multiedges=multiedges, weighted=weighted)

        elif format == 'incidence_matrix':
            from .graph_input import from_oriented_incidence_matrix
            from_oriented_incidence_matrix(self, data, loops=loops, multiedges=multiedges, weighted=weighted)

        elif format == 'DiGraph':
            if loops is None:
                loops = data.allows_loops()
            elif not loops and data.has_loops():
                raise ValueError("the digraph was built with loops=False but input data has a loop")
            if multiedges is None:
                multiedges = data.allows_multiple_edges()
            elif not multiedges:
                e = data.edges(labels=False, sort=False)
                if len(e) != len(set(e)):
                    raise ValueError("no multiple edges but input digraph"
                                     " has multiple edges")
            self.allow_multiple_edges(multiedges, check=False)
            self.allow_loops(loops, check=False)
            if weighted is None:
                weighted = data.weighted()
            if data.get_pos() is not None:
                pos = data.get_pos()
            self.set_vertices(data.get_vertices())
            data._backend.subgraph_given_vertices(self._backend, data)
            self.name(data.name())
        elif format == 'rule':
            f = data[1]
            if loops is None:
                loops = any(f(v, v) for v in data[0])
            if weighted is None:
                weighted = False
            self.allow_multiple_edges(bool(multiedges), check=False)
            self.allow_loops(loops, check=False)
            self.add_vertices(data[0])
            self.add_edges((u, v) for u in data[0] for v in data[0] if f(u, v))

        elif format == "vertices_and_edges":
            self.allow_multiple_edges(bool(multiedges), check=False)
            self.allow_loops(bool(loops), check=False)
            self.add_vertices(data[0])
            self.add_edges(data[1])

        elif format == 'dict_of_dicts':
            from .graph_input import from_dict_of_dicts
            from_dict_of_dicts(self, data, loops=loops, multiedges=multiedges, weighted=weighted,
                               convert_empty_dict_labels_to_None=False if convert_empty_dict_labels_to_None is None else convert_empty_dict_labels_to_None)

        elif format == 'dict_of_lists':
            from .graph_input import from_dict_of_lists
            from_dict_of_lists(self, data, loops=loops, multiedges=multiedges, weighted=weighted)

        elif format == 'NX':
            from sage.graphs.graph_input import from_networkx_graph
            from_networkx_graph(self, data,
                                weighted=weighted, multiedges=multiedges, loops=loops,
                                convert_empty_dict_labels_to_None=convert_empty_dict_labels_to_None)
            if weighted is None:
                weighted = self.allows_multiple_edges()

        elif format == 'igraph':
            if not data.is_directed():
                raise ValueError("a *directed* igraph graph was expected. To "
                                 "build an undirected graph, call the Graph "
                                 "constructor")

            self.add_vertices(range(data.vcount()))
            self.add_edges((e.source, e.target, e.attributes()) for e in data.es())

            if vertex_labels and 'name' in data.vertex_attributes():
                vs = data.vs()
                self.relabel({v: vs[v]['name'] for v in self})

        elif format == 'int':
            if weighted is None:
                weighted = False
            self.allow_loops(bool(loops), check=False)
            self.allow_multiple_edges(bool(multiedges),
                                      check=False)
            if data < 0:
                raise ValueError("the number of vertices cannot be strictly negative")
            elif data:
                self.add_vertices(range(data))
        elif format == 'list_of_edges':
            self.allow_multiple_edges(bool(multiedges),
                                      check=False)
            self.allow_loops(bool(loops), check=False)
            self.add_edges(data)
        else:
            raise ValueError("unknown input format '{}'".format(format))

        # weighted, multiedges, loops, verts and num_verts should now be set
        self._weighted = weighted

        if hash_labels is None and hasattr(data, '_hash_labels'):
            hash_labels = data._hash_labels
        self._hash_labels = hash_labels

        self._pos = copy(pos)

        if format != 'DiGraph' or name is not None:
            self.name(name)

        if data_structure == "static_sparse":
            from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            ib = StaticSparseBackend(self,
                                     loops=self.allows_loops(),
                                     multiedges=self.allows_multiple_edges())
            self._backend = ib
            self._immutable = True

    # Formats

    def dig6_string(self):
        r"""
        Return the ``dig6`` representation of the digraph as an ASCII string.

        This is only valid for simple (no multiple edges) digraphs on at most
        `2^{18} - 1 = 262143` vertices.

        .. NOTE::

            As the ``dig6`` format only handles graphs with vertex set `\{0,
            \ldots, n-1\}`, a :meth:`relabelled copy
            <sage.graphs.generic_graph.GenericGraph.relabel>` will be encoded,
            if necessary.

        .. SEEALSO::

            * :meth:`~sage.graphs.graph.Graph.graph6_string` -- a similar string
              format for undirected graphs

        EXAMPLES::

            sage: D = DiGraph({0: [1, 2], 1: [2], 2: [3], 3: [0]})
            sage: D.dig6_string()
            'CW`_'
            sage: L = DiGraph({0: [1, 2], 1: [2], 2: [3], 3: [3]})
            sage: L.dig6_string()
            'CW`C'

        TESTS::

            sage: DiGraph().dig6_string()
            '?'
        """
        n = self.order()
        if n > 262143:
            raise ValueError('dig6 format supports graphs on 0 to 262143 vertices only')
        elif self.has_multiple_edges():
            raise ValueError('dig6 format does not support multiple edges')
        return generic_graph_pyx.small_integer_to_graph6(n) + generic_graph_pyx.binary_string_to_graph6(self._bit_vector())

    # Attributes

    def is_directed(self):
        """
        Since digraph is directed, return ``True``.

        EXAMPLES::

            sage: DiGraph().is_directed()
            True
        """
        return True

    # Properties

    def is_directed_acyclic(self, certificate=False):
        r"""
        Check whether the digraph is acyclic or not.

        A directed graph is acyclic if for any vertex `v`, there is no directed
        path that starts and ends at `v`. Every directed acyclic graph (DAG)
        corresponds to a partial ordering of its vertices, however multiple dags
        may lead to the same partial ordering.

        INPUT:

        - ``certificate`` -- boolean (default: ``False``); whether to return a
          certificate

        OUTPUT:

        * When ``certificate=False``, returns a boolean value

        * When ``certificate=True``:

          * If the graph is acyclic, returns a pair ``(True, ordering)`` where
            ``ordering`` is a list of the vertices such that `u` appears
            before `v` in ``ordering`` if `uv` is an edge.

          * Else, returns a pair ``(False, cycle)`` where ``cycle`` is a list of
            vertices representing a circuit in the graph.

        EXAMPLES:

        At first, the following graph is acyclic::

            sage: D = DiGraph({0:[1, 2, 3], 4:[2, 5], 1:[8], 2:[7], 3:[7], 5:[6,7], 7:[8], 6:[9], 8:[10], 9:[10]})
            sage: D.plot(layout='circular').show()                                      # needs sage.plot
            sage: D.is_directed_acyclic()
            True

        Adding an edge from `9` to `7` does not change it::

            sage: D.add_edge(9, 7)
            sage: D.is_directed_acyclic()
            True

        We can obtain as a proof an ordering of the vertices such that `u`
        appears before `v` if `uv` is an edge of the graph::

            sage: D.is_directed_acyclic(certificate=True)
            (True, [4, 5, 6, 9, 0, 1, 2, 3, 7, 8, 10])

        Adding an edge from 7 to 4, though, makes a difference::

            sage: D.add_edge(7, 4)
            sage: D.is_directed_acyclic()
            False

        Indeed, it creates a circuit `7, 4, 5`::

            sage: D.is_directed_acyclic(certificate=True)
            (False, [7, 4, 5])

        Checking acyclic graphs are indeed acyclic ::

            sage: def random_acyclic(n, p):
            ....:  g = graphs.RandomGNP(n, p)
            ....:  h = DiGraph()
            ....:  h.add_edges(((u, v) if u < v else (v, u)) for u, v in g.edge_iterator(labels=False))
            ....:  return h
            ...
            sage: all(random_acyclic(100, .2).is_directed_acyclic()    # long time
            ....:      for i in range(50))
            True

        TESTS:

        What about loops? ::

            sage: g = digraphs.ButterflyGraph(3)
            sage: g.allow_loops(True)
            sage: g.is_directed_acyclic()
            True
            sage: g.add_edge(0, 0)
            sage: g.is_directed_acyclic()
            False
        """
        return self._backend.is_directed_acyclic(certificate=certificate)

    def to_directed(self):
        """
        Since the graph is already directed, simply returns a copy of itself.

        EXAMPLES::

            sage: DiGraph({0: [1, 2, 3], 4: [5, 1]}).to_directed()
            Digraph on 6 vertices
        """
        return self.copy()

    def to_undirected(self, data_structure=None, sparse=None):
        """
        Return an undirected version of the graph.

        Every directed edge becomes an edge.

        INPUT:

        - ``data_structure`` -- string (default: ``None``); one of
          ``'sparse'``, ``'static_sparse'``, or ``'dense'``. See the
          documentation of :class:`Graph` or :class:`DiGraph`.

        - ``sparse`` -- boolean (default: ``None``); ``sparse=True`` is an
          alias for ``data_structure="sparse"``, and ``sparse=False`` is an
          alias for ``data_structure="dense"``.

        EXAMPLES::

            sage: D = DiGraph({0: [1, 2], 1: [0]})
            sage: G = D.to_undirected()
            sage: D.edges(sort=True, labels=False)
            [(0, 1), (0, 2), (1, 0)]
            sage: G.edges(sort=True, labels=False)
            [(0, 1), (0, 2)]

        TESTS:

        Immutable graphs yield immutable graphs (:issue:`17005`)::

            sage: DiGraph([[1, 2]], immutable=True).to_undirected()._backend
            <sage.graphs.base.static_sparse_backend.StaticSparseBackend object at ...>

        Vertex labels will be retained (:issue:`14708`)::

            sage: D.set_vertex(0, 'foo')
            sage: G = D.to_undirected()
            sage: D.get_vertices()
            {0: 'foo', 1: None, 2: None}
            sage: G.get_vertices()
            {0: 'foo', 1: None, 2: None}
        """
        if sparse is not None:
            if data_structure is not None:
                raise ValueError("the 'sparse' argument is an alias for "
                                 "'data_structure'. Please do not define both")
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
        from sage.graphs.graph import Graph
        G = Graph(name=self.name(),
                  pos=self._pos,
                  multiedges=self.allows_multiple_edges(),
                  loops=self.allows_loops(),
                  data_structure=(data_structure if data_structure != "static_sparse"
                                  else "sparse"))  # we need a mutable copy first

        G.add_vertices(self.vertex_iterator())
        G.set_vertices(self.get_vertices())
        G.add_edges(self.edge_iterator())
        if hasattr(self, '_embedding'):
            G._embedding = copy(self._embedding)
        G._weighted = self._weighted

        if data_structure == "static_sparse":
            G = G.copy(data_structure=data_structure)

        return G

    # Edge Handlers

    def incoming_edge_iterator(self, vertices, labels=True):
        """
        Return an iterator over all arriving edges from vertices.

        INPUT:

        - ``vertices`` -- a vertex or a list of vertices

        - ``labels`` -- boolean (default: ``True``); whether to return edges as
          pairs of vertices, or as triples containing the labels

        EXAMPLES::

            sage: D = DiGraph({0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1]})
            sage: for a in D.incoming_edge_iterator([0]):
            ....:     print(a)
            (1, 0, None)
            (4, 0, None)
        """
        if vertices is None:
            vertices = self
        elif vertices in self:
            vertices = [vertices]
        else:
            vertices = [v for v in vertices if v in self]
        return self._backend.iterator_in_edges(vertices, labels)

    def incoming_edges(self, vertices, labels=True):
        """
        Return a list of edges arriving at vertices.

        INPUT:

        - ``vertices`` -- a vertex or a list of vertices

        - ``labels`` -- boolean (default: ``True``); whether to return edges as
          pairs of vertices, or as triples containing the labels

        EXAMPLES::

            sage: D = DiGraph({0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1]})
            sage: D.incoming_edges([0])
            [(1, 0, None), (4, 0, None)]
        """
        return list(self.incoming_edge_iterator(vertices, labels=labels))

    def outgoing_edge_iterator(self, vertices, labels=True):
        """
        Return an iterator over all departing edges from vertices.

        INPUT:

        - ``vertices`` -- a vertex or a list of vertices

        - ``labels`` -- boolean (default: ``True``); whether to return edges as
          pairs of vertices, or as triples containing the labels

        EXAMPLES::

            sage: D = DiGraph({0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1]})
            sage: for a in D.outgoing_edge_iterator([0]):
            ....:     print(a)
            (0, 1, None)
            (0, 2, None)
            (0, 3, None)
        """
        if vertices is None:
            vertices = self
        elif vertices in self:
            vertices = [vertices]
        else:
            vertices = [v for v in vertices if v in self]
        return self._backend.iterator_out_edges(vertices, labels)

    def outgoing_edges(self, vertices, labels=True):
        """
        Return a list of edges departing from vertices.

        INPUT:

        - ``vertices`` -- a vertex or a list of vertices

        - ``labels`` -- boolean (default: ``True``); whether to return edges as
          pairs of vertices, or as triples containing the labels

        EXAMPLES::

            sage: D = DiGraph({0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1]})
            sage: D.outgoing_edges([0])
            [(0, 1, None), (0, 2, None), (0, 3, None)]
        """
        return list(self.outgoing_edge_iterator(vertices, labels=labels))

    def neighbor_in_iterator(self, vertex):
        """
        Return an iterator over the in-neighbors of ``vertex``.

        A vertex `u` is an in-neighbor of a vertex `v` if `uv` in an edge.

        EXAMPLES::

            sage: D = DiGraph({0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1]})
            sage: for a in D.neighbor_in_iterator(0):
            ....:     print(a)
            1
            4

        TESTS:

        With multiple edges, check that the neighbors are listed only once::

            sage: D = DiGraph([[0, 1, 2], [(0, 1), (0, 1), (1, 2) ]],multiedges=True)
            sage: list(D.neighbor_in_iterator(0))
            []
            sage: list(D.neighbor_in_iterator(1))
            [0]
            sage: list(D.neighbor_in_iterator(2))
            [1]

        Check that the iterator lists ``vertex`` in the presence of loop(s):

            sage: D = DiGraph([[0, 1, 2], [(0, 0), (0, 0), (1, 1)]],multiedges=True, loops=True)
            sage: list(D.neighbor_in_iterator(0))
            [0]
            sage: list(D.neighbor_in_iterator(1))
            [1]
            sage: list(D.neighbor_in_iterator(2))
            []

        """
        yield from self._backend.iterator_in_nbrs(vertex)

    def neighbors_in(self, vertex):
        """
        Return the list of the in-neighbors of a given vertex.

        A vertex `u` is an in-neighbor of a vertex `v` if `uv` in an edge.

        EXAMPLES::

            sage: D = DiGraph({0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1]})
            sage: D.neighbors_in(0)
            [1, 4]
        """
        return list(self.neighbor_in_iterator(vertex))

    def neighbor_out_iterator(self, vertex):
        """
        Return an iterator over the out-neighbors of a given vertex.

        A vertex `u` is an out-neighbor of a vertex `v` if `vu` in an edge.

        EXAMPLES::

            sage: D = DiGraph({0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1]})
            sage: for a in D.neighbor_out_iterator(0):
            ....:     print(a)
            1
            2
            3

        TESTS:

        With multiple edges, check that the neighbors are listed only once::

            sage: D = DiGraph([[0, 1, 2], [(0, 1), (0, 1), (1, 2) ]],multiedges=True)
            sage: list(D.neighbor_out_iterator(0))
            [1]
            sage: list(D.neighbor_out_iterator(1))
            [2]
            sage: list(D.neighbor_out_iterator(2))
            []

        Check that the iterator lists ``vertex`` in the presence of loop(s):

            sage: D = DiGraph([[0, 1, 2], [(0, 0), (0, 0), (1, 1)]],multiedges=True, loops=True)
            sage: list(D.neighbor_out_iterator(0))
            [0]
            sage: list(D.neighbor_out_iterator(1))
            [1]
            sage: list(D.neighbor_out_iterator(2))
            []

        """
        yield from self._backend.iterator_out_nbrs(vertex)

    def neighbors_out(self, vertex):
        """
        Return the list of the out-neighbors of a given vertex.

        A vertex `u` is an out-neighbor of a vertex `v` if `vu` in an edge.

        EXAMPLES::

            sage: D = DiGraph({0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1]})
            sage: D.neighbors_out(0)
            [1, 2, 3]
        """
        return list(self.neighbor_out_iterator(vertex))

    # Degree functions

    def in_degree(self, vertices=None, labels=False):
        """
        Same as degree, but for in degree.

        EXAMPLES::

            sage: D = DiGraph({0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1]})
            sage: D.in_degree(vertices=[0, 1, 2], labels=True)
            {0: 2, 1: 2, 2: 2}
            sage: D.in_degree()
            [2, 2, 2, 2, 1, 1]
            sage: G = graphs.PetersenGraph().to_directed()
            sage: G.in_degree(0)
            3
        """
        if vertices in self:
            return self._backend.in_degree(vertices)
        elif labels:
            return dict(self.in_degree_iterator(vertices, labels=labels))
        return list(self.in_degree_iterator(vertices, labels=labels))

    def in_degree_iterator(self, vertices=None, labels=False):
        """
        Same as degree_iterator, but for in degree.

        EXAMPLES::

            sage: D = graphs.Grid2dGraph(2,4).to_directed()
            sage: sorted(D.in_degree_iterator())
            [2, 2, 2, 2, 3, 3, 3, 3]
            sage: sorted(D.in_degree_iterator(labels=True))
            [((0, 0), 2),
             ((0, 1), 3),
             ((0, 2), 3),
             ((0, 3), 2),
             ((1, 0), 2),
             ((1, 1), 3),
             ((1, 2), 3),
             ((1, 3), 2)]
        """
        if vertices is None:
            vertices = self.vertex_iterator()
        if labels:
            for v in vertices:
                yield (v, self.in_degree(v))
        else:
            for v in vertices:
                yield self.in_degree(v)

    def in_degree_sequence(self):
        r"""
        Return the in-degree sequence.

        EXAMPLES:

        The in-degree sequences of two digraphs::

            sage: g = DiGraph({1: [2, 5, 6], 2: [3, 6], 3: [4, 6], 4: [6], 5: [4, 6]})
            sage: g.in_degree_sequence()
            [5, 2, 1, 1, 1, 0]

        ::

            sage: V = [2, 3, 5, 7, 8, 9, 10, 11]
            sage: E = [[], [8, 10], [11], [8, 11], [9], [], [], [2, 9, 10]]
            sage: g = DiGraph(dict(zip(V, E)))
            sage: g.in_degree_sequence()
            [2, 2, 2, 2, 1, 0, 0, 0]
        """
        return sorted(self.in_degree_iterator(), reverse=True)

    def out_degree(self, vertices=None, labels=False):
        """
        Same as degree, but for out degree.

        EXAMPLES::

            sage: D = DiGraph({0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1]})
            sage: D.out_degree(vertices=[0, 1 ,2], labels=True)
            {0: 3, 1: 2, 2: 1}
            sage: D.out_degree()
            [3, 2, 1, 1, 2, 1]
            sage: D.out_degree(2)
            1
        """
        if vertices in self:
            return self._backend.out_degree(vertices)
        elif labels:
            return dict(self.out_degree_iterator(vertices, labels=labels))
        return list(self.out_degree_iterator(vertices, labels=labels))

    def out_degree_iterator(self, vertices=None, labels=False):
        """
        Same as degree_iterator, but for out degree.

        EXAMPLES::

            sage: D = graphs.Grid2dGraph(2,4).to_directed()
            sage: sorted(D.out_degree_iterator())
            [2, 2, 2, 2, 3, 3, 3, 3]
            sage: sorted(D.out_degree_iterator(labels=True))
            [((0, 0), 2),
             ((0, 1), 3),
             ((0, 2), 3),
             ((0, 3), 2),
             ((1, 0), 2),
             ((1, 1), 3),
             ((1, 2), 3),
             ((1, 3), 2)]
        """
        if vertices is None:
            vertices = self.vertex_iterator()
        if labels:
            for v in vertices:
                yield (v, self.out_degree(v))
        else:
            for v in vertices:
                yield self.out_degree(v)

    def out_degree_sequence(self):
        r"""
        Return the outdegree sequence of this digraph.

        EXAMPLES:

        The outdegree sequences of two digraphs::

            sage: g = DiGraph({1: [2, 5, 6], 2: [3, 6], 3: [4, 6], 4: [6], 5: [4, 6]})
            sage: g.out_degree_sequence()
            [3, 2, 2, 2, 1, 0]

        ::

            sage: V = [2, 3, 5, 7, 8, 9, 10, 11]
            sage: E = [[], [8, 10], [11], [8, 11], [9], [], [], [2, 9, 10]]
            sage: g = DiGraph(dict(zip(V, E)))
            sage: g.out_degree_sequence()
            [3, 2, 2, 1, 1, 0, 0, 0]
        """
        return sorted(self.out_degree_iterator(), reverse=True)

    def sources(self):
        r"""
        Return a list of sources of the digraph.

        OUTPUT: list of the vertices of the digraph that have no edges going into them

        EXAMPLES::

            sage: G = DiGraph({1: {3: ['a']}, 2: {3: ['b']}})
            sage: G.sources()
            [1, 2]
            sage: T = DiGraph({1: {}})
            sage: T.sources()
            [1]
        """
        return [x for x in self if not self.in_degree(x)]

    def sinks(self):
        """
        Return a list of sinks of the digraph.

        OUTPUT: list of the vertices of the digraph that have no edges beginning at them

        EXAMPLES::

            sage: G = DiGraph({1: {3: ['a']}, 2: {3: ['b']}})
            sage: G.sinks()
            [3]
            sage: T = DiGraph({1: {}})
            sage: T.sinks()
            [1]
        """
        return [x for x in self if not self.out_degree(x)]

    def degree_polynomial(self):
        r"""
        Return the generating polynomial of degrees of vertices in ``self``.

        This is the sum

        .. MATH::

            \sum_{v \in G} x^{\operatorname{in}(v)} y^{\operatorname{out}(v)},

        where ``in(v)`` and ``out(v)`` are the number of incoming and outgoing
        edges at vertex `v` in the digraph `G`.

        Because this polynomial is multiplicative for Cartesian product of
        digraphs, it is useful to help see if the digraph can be isomorphic to a
        Cartesian product.

        .. SEEALSO::

            :meth:`num_verts` for the value at `(x, y) = (1, 1)`

        EXAMPLES::

            sage: G = posets.PentagonPoset().hasse_diagram()                            # needs sage.modules
            sage: G.degree_polynomial()                                                 # needs sage.modules
            x^2 + 3*x*y + y^2

            sage: G = posets.BooleanLattice(4).hasse_diagram()
            sage: G.degree_polynomial().factor()                                        # needs sage.libs.pari
            (x + y)^4
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        R = PolynomialRing(ZZ, 'x,y')
        x, y = R.gens()
        return R.sum(x ** self.in_degree(v) * y ** self.out_degree(v) for v in self)

    def feedback_edge_set(self, constraint_generation=True, value_only=False,
                          solver=None, verbose=0, *, integrality_tolerance=1e-3):
        r"""
        Compute the minimum feedback edge set of a digraph (also called
        feedback arc set).

        The minimum feedback edge set of a digraph is a set of edges that
        intersect all the circuits of the digraph.  Equivalently, a minimum
        feedback arc set of a DiGraph is a set `S` of arcs such that the digraph
        `G - S` is acyclic. For more information, see the
        :wikipedia:`Feedback_arc_set`.

        INPUT:

        - ``value_only`` -- boolean (default: ``False``)

          - When set to ``True``, only the minimum cardinal of a minimum edge
            set is returned.

          - When set to ``False``, the ``Set`` of edges of a minimal edge set is
            returned.

        - ``constraint_generation`` -- boolean (default: ``True``); whether to
          use constraint generation when solving the Mixed Integer Linear
          Program.

        - ``solver`` -- string (default: ``None``); specifies a Mixed Integer
          Linear Programming (MILP) solver to be used. If set to ``None``, the
          default one is used. For more information on MILP solvers and which
          default solver is used, see the method :meth:`solve
          <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
          :class:`MixedIntegerLinearProgram
          <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: 0); sets the level of
          verbosity. Set to 0 by default, which means quiet.

        - ``integrality_tolerance`` -- float; parameter for use with MILP
          solvers over an inexact base ring; see
          :meth:`MixedIntegerLinearProgram.get_values`.

        ALGORITHM:

        This problem is solved using Linear Programming, in two different
        ways. The first one is to solve the following formulation:

        .. MATH::

            \mbox{Minimize : }&\sum_{(u,v)\in G} b_{(u,v)}\\
            \mbox{Such that : }&\\
            &\forall (u,v)\in G, d_u-d_v+ n \cdot b_{(u,v)}\geq 0\\
            &\forall u\in G, 0\leq d_u\leq |G|\\

        An explanation:

        An acyclic digraph can be seen as a poset, and every poset has a linear
        extension. This means that in any acyclic digraph the vertices can be
        ordered with a total order `<` in such a way that if `(u,v) \in G`, then
        `u < v`.

        Thus, this linear program is built in order to assign to each vertex `v`
        a number `d_v \in [0,\dots,n-1]` such that if there exists an edge
        `(u, v) \in G` such that `d_v < d_u`, then the edge `(u,v)` is removed.

        The number of edges removed is then minimized, which is the objective.

        (Constraint Generation)

        If the parameter ``constraint_generation`` is enabled, a more efficient
        formulation is used :

        .. MATH::

            \mbox{Minimize : }&\sum_{(u,v)\in G} b_{(u,v)}\\
            \mbox{Such that : }&\\
            &\forall C\text{ circuits }\subseteq G, \sum_{uv\in C}b_{(u,v)}\geq 1\\

        As the number of circuits contained in a graph is exponential, this LP
        is solved through constraint generation. This means that the solver is
        sequentially asked to solved the problem, knowing only a portion of the
        circuits contained in `G`, each time adding to the list of its
        constraints the circuit which its last answer had left intact.

        EXAMPLES:

        If the digraph is created from a graph, and hence is symmetric (if `uv`
        is an edge, then `vu` is an edge too), then obviously the cardinality of
        its feedback arc set is the number of edges in the first graph::

            sage: cycle = graphs.CycleGraph(5)
            sage: dcycle = DiGraph(cycle)
            sage: cycle.size()
            5
            sage: dcycle.feedback_edge_set(value_only=True)                             # needs sage.numerical.mip
            5

        And in this situation, for any edge `uv` of the first graph, `uv` of
        `vu` is in the returned feedback arc set::

           sage: g = graphs.RandomGNP(5,.3)
           sage: while not g.num_edges():
           ....:     g = graphs.RandomGNP(5,.3)
           sage: dg = DiGraph(g)
           sage: feedback = dg.feedback_edge_set()                                      # needs sage.numerical.mip
           sage: u,v,l = next(g.edge_iterator())
           sage: (u,v) in feedback or (v,u) in feedback                                 # needs sage.numerical.mip
           True

        TESTS:

        Comparing with/without constraint generation. Also double-checks issue
        :issue:`12833`::

            sage: for i in range(20):                                                   # needs sage.numerical.mip
            ....:     g = digraphs.RandomDirectedGNP(10, .3)
            ....:     x = g.feedback_edge_set(value_only=True)
            ....:     y = g.feedback_edge_set(value_only=True,
            ....:            constraint_generation=False)
            ....:     if x != y:
            ....:         print("Oh my, oh my !")
            ....:         break

        Loops are part of the feedback edge set (:issue:`23989`)::

            sage: # needs sage.combinat
            sage: D = digraphs.DeBruijn(2, 2)
            sage: sorted(D.loops(labels=None))
            [('00', '00'), ('11', '11')]
            sage: FAS = D.feedback_edge_set(value_only=False)                           # needs sage.numerical.mip
            sage: all(l in FAS for l in D.loops(labels=None))                           # needs sage.numerical.mip
            True
            sage: FAS2 = D.feedback_edge_set(value_only=False,                          # needs sage.numerical.mip
            ....:                            constraint_generation=False)
            sage: len(FAS) == len(FAS2)                                                 # needs sage.numerical.mip
            True

        Check that multi-edges are properly taken into account::

            sage: cycle = graphs.CycleGraph(5)
            sage: dcycle = DiGraph(cycle)
            sage: dcycle.feedback_edge_set(value_only=True)                             # needs sage.numerical.mip
            5
            sage: dcycle.allow_multiple_edges(True)
            sage: dcycle.add_edges(dcycle.edges(sort=True))
            sage: dcycle.feedback_edge_set(value_only=True)                             # needs sage.numerical.mip
            10
            sage: dcycle.feedback_edge_set(value_only=True,                             # needs sage.numerical.mip
            ....:                          constraint_generation=False)
            10

        Strongly connected components are well handled (:issue:`23989`)::

            sage: g = digraphs.Circuit(3) * 2
            sage: g.add_edge(0, 3)
            sage: g.feedback_edge_set(value_only=True)                                  # needs sage.numerical.mip
            2
        """
        # It would be a pity to start a LP if the digraph is already acyclic
        if self.is_directed_acyclic():
            return 0 if value_only else []

        if self.has_loops():
            # We solve the problem on a copy without loops of the digraph
            D = DiGraph(self.edges(sort=False), multiedges=self.allows_multiple_edges(), loops=True)
            loops = D.loops(labels=None)
            D.delete_edges(loops)
            D.allow_loops(False, check=False)
            FAS = D.feedback_edge_set(constraint_generation=constraint_generation,
                                      value_only=value_only, solver=solver, verbose=verbose,
                                      integrality_tolerance=integrality_tolerance)
            if value_only:
                return FAS + len(loops)
            return FAS + loops

        if not self.is_strongly_connected():
            # If the digraph is not strongly connected, we solve the problem on
            # each of its strongly connected components

            FAS = 0 if value_only else []

            for h in self.strongly_connected_components_subgraphs():
                if not h.size():
                    continue
                if value_only:
                    FAS += h.feedback_edge_set(constraint_generation=constraint_generation,
                                               value_only=True, solver=solver, verbose=verbose,
                                               integrality_tolerance=integrality_tolerance)
                else:
                    FAS.extend(h.feedback_edge_set(constraint_generation=constraint_generation,
                                                   value_only=False, solver=solver, verbose=verbose,
                                                   integrality_tolerance=integrality_tolerance))
            return FAS

        from sage.numerical.mip import MixedIntegerLinearProgram

        ########################################
        # Constraint Generation Implementation #
        ########################################
        if constraint_generation:

            p = MixedIntegerLinearProgram(constraint_generation=True,
                                          maximization=False, solver=solver)

            # A variable for each edge
            b = p.new_variable(binary=True)

            # Variables are binary, and their coefficient in the objective is
            # the number of occurrences of the corresponding edge, so 1 if the
            # graph is simple
            p.set_objective(p.sum(b[e] for e in self.edge_iterator(labels=False)))

            # For as long as we do not break because the digraph is acyclic....
            while True:

                # Building the graph without the edges removed by the MILP
                p.solve(log=verbose)
                val = p.get_values(b, convert=bool, tolerance=integrality_tolerance)
                h = DiGraph([e for e in self.edge_iterator(labels=False) if not val[e]],
                            format='list_of_edges')

                # Is the digraph acyclic ?
                isok, certificate = h.is_directed_acyclic(certificate=True)

                # If so, we are done !
                if isok:
                    if value_only:
                        return sum(1 for e in self.edge_iterator(labels=False) if val[e])
                    # listing the edges contained in the MFAS
                    return [e for e in self.edge_iterator(labels=False) if val[e]]

                # There is a circuit left. Let's add the corresponding
                # constraint !
                while not isok:

                    if verbose:
                        print("Adding a constraint on circuit : {}".format(certificate))

                    edges = zip(certificate, certificate[1:] + [certificate[0]])
                    p.add_constraint(p.sum(b[e] for e in edges), min=1)

                    # Is there another edge disjoint circuit ?
                    # for python3, we need to recreate the zip iterator
                    edges = zip(certificate, certificate[1:] + [certificate[0]])
                    h.delete_edges(edges)
                    isok, certificate = h.is_directed_acyclic(certificate=True)

        ######################################
        # Ordering-based MILP Implementation #
        ######################################
        else:
            p = MixedIntegerLinearProgram(maximization=False, solver=solver)

            b = p.new_variable(binary=True)
            d = p.new_variable(integer=True, nonnegative=True)

            n = self.order()

            for u, v in self.edge_iterator(labels=None):
                p.add_constraint(d[u] - d[v] + n * b[u, v], min=1)

            for v in self:
                p.add_constraint(d[v] <= n)

            p.set_objective(p.sum(b[e] for e in self.edge_iterator(labels=False)))

            p.solve(log=verbose)

            b_sol = p.get_values(b, convert=bool, tolerance=integrality_tolerance)

            if value_only:
                return sum(1 for e in self.edge_iterator(labels=False) if b_sol[e])
            return [e for e in self.edge_iterator(labels=False) if b_sol[e]]

    # Construction

    def reverse(self, immutable=None):
        """
        Return a copy of digraph with edges reversed in direction.

        INPUT:

        - ``immutable`` -- boolean (default: ``None``); whether to return an
          immutable digraph or not. By default (``None``), the returned digraph
          has the same setting than ``self``. That is, if ``self`` is immutable,
          the returned digraph also is.

        EXAMPLES::

            sage: adj = {0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1]}
            sage: D = DiGraph(adj)
            sage: R = D.reverse(); R
            Reverse of (): Digraph on 6 vertices
            sage: H = R.reverse()
            sage: adj == H.to_dictionary()
            True

        TESTS::

            sage: adj = {0: [1, 1], 1: [1]}
            sage: D = DiGraph(adj, immutable=True, multiedges=True, loops=True)
            sage: R = D.reverse()
            sage: R.is_immutable() and R.allows_loops() and R.allows_multiple_edges()
            True
            sage: adj == R.reverse().to_dictionary(multiple_edges=True)
            True

        Check the behavior of parameter ``immutable``::

            sage: D = DiGraph([(0, 1)], immutable=False)
            sage: R = D.reverse()
            sage: R.is_immutable()
            False
            sage: R = D.reverse(immutable=True)
            sage: R.is_immutable()
            True
            sage: H = R.reverse()
            sage: H.is_immutable()
            True
            sage: H = R.reverse(immutable=False)
            sage: H.is_immutable()
            False
        """
        from sage.graphs.base.dense_graph import DenseGraphBackend
        if isinstance(self._backend, DenseGraphBackend):
            data_structure = "dense"
        else:
            data_structure = "sparse"

        H = DiGraph(data_structure=data_structure,
                    multiedges=self.allows_multiple_edges(), loops=self.allows_loops(),
                    pos=copy(self._pos), weighted=self.weighted(),
                    hash_labels=self._hash_labels)
        H.add_vertices(self)
        H.add_edges((v, u, d) for u, v, d in self.edge_iterator())
        name = self.name()
        if name is None:
            name = ''
        H.name("Reverse of (%s)" % name)

        attributes_to_copy = ('_assoc', '_embedding')
        for attr in attributes_to_copy:
            if hasattr(self, attr):
                copy_attr = {}
                old_attr = getattr(self, attr)
                if isinstance(old_attr, dict):
                    for v, value in old_attr.items():
                        try:
                            copy_attr[v] = value.copy()
                        except AttributeError:
                            copy_attr[v] = copy(value)
                    setattr(H, attr, copy_attr)
                else:
                    setattr(H, attr, copy(old_attr))

        if immutable or (immutable is None and self.is_immutable()):
            return H.copy(immutable=True)

        return H

    def reverse_edge(self, u, v=None, label=None, inplace=True, multiedges=None):
        """
        Reverse the edge from `u` to `v`.

        INPUT:

        - ``inplace`` -- boolean (default: ``True``); if ``False``, a new
          digraph is created and returned as output, otherwise ``self`` is
          modified.

        - ``multiedges`` -- boolean (default: ``None``); how to decide what
          should be done in case of doubt (for instance when edge `(1,2)` is to
          be reversed in a graph while `(2,1)` already exists):

          - If set to ``True``, input graph will be forced to allow parallel
            edges if necessary and edge `(1,2)` will appear twice in the graph.

          - If set to ``False``, only one edge `(1,2)` will remain in the graph
            after `(2,1)` is reversed. Besides, the label of edge `(1,2)` will
            be overwritten with the label of edge `(2,1)`.

          The default behaviour (``multiedges = None``) will raise an exception
          each time a subjective decision (setting ``multiedges`` to ``True`` or
          ``False``) is necessary to perform the operation.

        The following forms are all accepted:

        - D.reverse_edge( 1, 2 )
        - D.reverse_edge( (1, 2) )
        - D.reverse_edge( [1, 2] )
        - D.reverse_edge( 1, 2, 'label' )
        - D.reverse_edge( ( 1, 2, 'label') )
        - D.reverse_edge( [1, 2, 'label'] )
        - D.reverse_edge( ( 1, 2), label='label' )

        EXAMPLES:

        If ``inplace`` is ``True`` (default), ``self`` is modified::

            sage: D = DiGraph([(0, 1 ,2)])
            sage: D.reverse_edge(0, 1)
            sage: D.edges(sort=True)
            [(1, 0, 2)]

        If ``inplace`` is ``False``, ``self`` is not modified and a new digraph
        is returned::

            sage: D = DiGraph([(0, 1, 2)])
            sage: re = D.reverse_edge(0, 1, inplace=False)
            sage: re.edges(sort=True)
            [(1, 0, 2)]
            sage: D.edges(sort=True)
            [(0, 1, 2)]

        If ``multiedges`` is ``True``, ``self`` will be forced to allow parallel
        edges when and only when it is necessary::

            sage: D = DiGraph([(1, 2, 'A'), (2, 1, 'A'), (2, 3, None)])
            sage: D.reverse_edge(1, 2, multiedges=True)
            sage: D.edges(sort=True)
            [(2, 1, 'A'), (2, 1, 'A'), (2, 3, None)]
            sage: D.allows_multiple_edges()
            True

        Even if ``multiedges`` is ``True``, ``self`` will not be forced to allow
        parallel edges when it is not necessary::

            sage: D = DiGraph( [(1, 2, 'A'), (2, 1, 'A'), (2, 3, None)] )
            sage: D.reverse_edge(2, 3, multiedges=True)
            sage: D.edges(sort=True)
            [(1, 2, 'A'), (2, 1, 'A'), (3, 2, None)]
            sage: D.allows_multiple_edges()
            False

        If user specifies ``multiedges = False``, ``self`` will not be forced to
        allow parallel edges and a parallel edge will get deleted::

            sage: D = DiGraph( [(1, 2, 'A'), (2, 1, 'A'), (2, 3, None)] )
            sage: D.edges(sort=True)
            [(1, 2, 'A'), (2, 1, 'A'), (2, 3, None)]
            sage: D.reverse_edge(1, 2, multiedges=False)
            sage: D.edges(sort=True)
            [(2, 1, 'A'), (2, 3, None)]

        Note that in the following graph, specifying ``multiedges = False`` will
        result in overwriting the label of `(1, 2)` with the label of `(2, 1)`::

            sage: D = DiGraph( [(1, 2, 'B'), (2, 1, 'A'), (2, 3, None)] )
            sage: D.edges(sort=True)
            [(1, 2, 'B'), (2, 1, 'A'), (2, 3, None)]
            sage: D.reverse_edge(2, 1, multiedges=False)
            sage: D.edges(sort=True)
            [(1, 2, 'A'), (2, 3, None)]

        If input edge in digraph has weight/label, then the weight/label should
        be preserved in the output digraph.  User does not need to specify the
        weight/label when calling function::

            sage: D = DiGraph([[0, 1, 2], [1, 2, 1]], weighted=True)
            sage: D.reverse_edge(0, 1)
            sage: D.edges(sort=True)
            [(1, 0, 2), (1, 2, 1)]
            sage: re = D.reverse_edge([1, 2], inplace=False)
            sage: re.edges(sort=True)
            [(1, 0, 2), (2, 1, 1)]

        If ``self`` has multiple copies (parallel edges) of the input edge, only
        1 of the parallel edges is reversed::

            sage: D = DiGraph([(0, 1, '01'), (0, 1, '01'), (0, 1, 'cat'), (1, 2, '12')], weighted=True, multiedges=True)
            sage: re = D.reverse_edge([0, 1, '01'], inplace=False)
            sage: re.edges(sort=True)
            [(0, 1, '01'), (0, 1, 'cat'), (1, 0, '01'), (1, 2, '12')]

        If ``self`` has multiple copies (parallel edges) of the input edge but
        with distinct labels and no input label is specified, only 1 of the
        parallel edges is reversed (the edge that is labeled by the first label
        on the list returned by :meth:`.edge_label`)::

            sage: D = DiGraph([(0, 1, 'A'), (0, 1, 'B'), (0, 1, 'mouse'), (0, 1, 'cat')], multiedges=true)
            sage: D.edge_label(0, 1)
            ['cat', 'mouse', 'B', 'A']
            sage: D.reverse_edge(0, 1)
            sage: D.edges(sort=True)
            [(0, 1, 'A'), (0, 1, 'B'), (0, 1, 'mouse'), (1, 0, 'cat')]

        Finally, an exception is raised when Sage does not know how to choose
        between allowing multiple edges and losing some data::

            sage: D = DiGraph([(0, 1, 'A'), (1, 0, 'B')])
            sage: D.reverse_edge(0, 1)
            Traceback (most recent call last):
            ...
            ValueError: reversing the given edge is about to create two parallel
            edges but input digraph doesn't allow them - User needs to specify
            multiedges is True or False.

        The following syntax is supported, but note that you must use the
        ``label`` keyword::

            sage: D = DiGraph()
            sage: D.add_edge((1, 2), label='label')
            sage: D.edges(sort=True)
            [(1, 2, 'label')]
            sage: D.reverse_edge((1, 2), label='label')
            sage: D.edges(sort=True)
            [(2, 1, 'label')]
            sage: D.add_edge((1, 2), 'label')
            sage: D.edges(sort=False)
            [((1, 2), 'label', None), (2, 1, 'label')]
            sage: D.reverse_edge((1, 2), 'label')
            sage: D.edges(sort=False)
            [('label', (1, 2), None), (2, 1, 'label')]

        TESTS::

            sage: D = DiGraph([(0, 1, None)])
            sage: D.reverse_edge(0, 1, 'mylabel')
            Traceback (most recent call last):
            ...
            ValueError: input edge must exist in the digraph
        """
        # Assigns the expected values to u,v, and label depending on the input.
        if label is None:
            if v is None:
                try:
                    u, v, label = u
                except Exception:
                    try:
                        u, v = u
                    except Exception:
                        pass
        else:
            if v is None:
                try:
                    u, v = u
                except Exception:
                    pass

        if not self.has_edge(u, v, label):
            raise ValueError("input edge must exist in the digraph")

        tempG = self if inplace else copy(self)

        if label is None:
            if not tempG.allows_multiple_edges():
                label = tempG.edge_label(u, v)
            else:
                # If digraph has parallel edges for input edge, pick the first
                # from the labels on the list
                label = tempG.edge_label(u, v)[0]

        if ((not tempG.allows_multiple_edges()) and (tempG.has_edge(v, u))):
            # If user wants to force digraph to allow parallel edges
            if multiedges:
                tempG.allow_multiple_edges(True)
                tempG.delete_edge(u, v, label)
                tempG.add_edge(v, u, label)

            # If user does not want to force digraph to allow parallel edges,
            # we delete edge u to v and overwrite v,u with the label of u,v
            elif multiedges is False:
                tempG.delete_edge(u, v, label)
                tempG.set_edge_label(v, u, label)

            # User is supposed to specify multiedges True or False
            else:
                raise ValueError("reversing the given edge is about to "
                                 "create two parallel edges but input digraph "
                                 "doesn't allow them - User needs to specify "
                                 "multiedges is True or False.")
        else:
            tempG.delete_edge(u, v, label)
            tempG.add_edge(v, u, label)

        if not inplace:
            return tempG

    def reverse_edges(self, edges, inplace=True, multiedges=None):
        """
        Reverse a list of edges.

        INPUT:

        - ``edges`` -- list of edges in the DiGraph

        - ``inplace`` -- boolean (default: ``True``); if ``False``, a new
          digraph is created and returned as output, otherwise ``self`` is
          modified

        - ``multiedges`` -- boolean (default: ``None``); if ``True``, input
          graph will be forced to allow parallel edges when necessary (for more
          information see the documentation of :meth:`~DiGraph.reverse_edge`)

        .. SEEALSO::

            :meth:`~DiGraph.reverse_edge` -- reverses a single edge

        EXAMPLES:

        If ``inplace`` is ``True`` (default), ``self`` is modified::

            sage: D = DiGraph({ 0: [1, 1, 3], 2: [3, 3], 4: [1, 5]}, multiedges=true)
            sage: D.reverse_edges([[0, 1], [0, 3]])
            sage: D.reverse_edges([(2, 3), (4, 5)])
            sage: D.edges(sort=True)
            [(0, 1, None), (1, 0, None), (2, 3, None), (3, 0, None),
             (3, 2, None), (4, 1, None), (5, 4, None)]

        If ``inplace`` is ``False``, ``self`` is not modified and a new digraph
        is returned::

            sage: D = DiGraph([(0, 1, 'A'), (1, 0, 'B'), (1, 2, 'C')])
            sage: re = D.reverse_edges([(0, 1), (1, 2)],
            ....:                       inplace=False,
            ....:                       multiedges=True)
            sage: re.edges(sort=True)
            [(1, 0, 'A'), (1, 0, 'B'), (2, 1, 'C')]
            sage: D.edges(sort=True)
            [(0, 1, 'A'), (1, 0, 'B'), (1, 2, 'C')]
            sage: D.allows_multiple_edges()
            False
            sage: re.allows_multiple_edges()
            True

        If ``multiedges`` is ``True``, ``self`` will be forced to allow parallel
        edges when and only when it is necessary::

            sage: D = DiGraph([(1, 2, 'A'), (2, 1, 'A'), (2, 3, None)])
            sage: D.reverse_edges([(1, 2), (2, 3)], multiedges=True)
            sage: D.edges(sort=True)
            [(2, 1, 'A'), (2, 1, 'A'), (3, 2, None)]
            sage: D.allows_multiple_edges()
            True

        Even if ``multiedges`` is ``True``, ``self`` will not be forced to allow
        parallel edges when it is not necessary::

            sage: D = DiGraph([(1, 2, 'A'), (2, 1, 'A'), (2, 3, None)])
            sage: D.reverse_edges([(2, 3)], multiedges=True)
            sage: D.edges(sort=True)
            [(1, 2, 'A'), (2, 1, 'A'), (3, 2, None)]
            sage: D.allows_multiple_edges()
            False

        If ``multiedges`` is ``False``, ``self`` will not be forced to allow
        parallel edges and an edge will get deleted::

            sage: D = DiGraph([(1, 2), (2, 1)])
            sage: D.edges(sort=True)
            [(1, 2, None), (2, 1, None)]
            sage: D.reverse_edges([(1, 2)], multiedges=False)
            sage: D.edges(sort=True)
            [(2, 1, None)]

        If input edge in digraph has weight/label, then the weight/label should
        be preserved in the output digraph.  User does not need to specify the
        weight/label when calling function::

            sage: D = DiGraph([(0, 1, '01'), (1, 2, 1), (2, 3, '23')], weighted=True)
            sage: D.reverse_edges([(0, 1, '01'), (1, 2), (2, 3)])
            sage: D.edges(sort=True)
            [(1, 0, '01'), (2, 1, 1), (3, 2, '23')]

        TESTS::

            sage: D = digraphs.Circuit(6)
            sage: D.reverse_edges(D.edges(sort=True), inplace=False).edges(sort=True)
            [(0, 5, None), (1, 0, None), (2, 1, None),
             (3, 2, None), (4, 3, None), (5, 4, None)]

            sage: D = digraphs.Kautz(2, 3)                                              # needs sage.combinat
            sage: Dr = D.reverse_edges(D.edges(sort=True), inplace=False,               # needs sage.combinat
            ....:                      multiedges=True)
            sage: Dr.edges(sort=True) == D.reverse().edges(sort=True)                   # needs sage.combinat
            True
        """
        tempG = self if inplace else copy(self)
        for e in edges:
            tempG.reverse_edge(e, inplace=True, multiedges=multiedges)
        if not inplace:
            return tempG

    # Distances

    def eccentricity(self, v=None, by_weight=False, algorithm=None,
                     weight_function=None, check_weight=True, dist_dict=None,
                     with_labels=False):
        """
        Return the eccentricity of vertex (or vertices) ``v``.

        The eccentricity of a vertex is the maximum distance to any other
        vertex.

        For more information and examples on how to use input variables, see
        :meth:`~GenericGraph.shortest_path_all_pairs`,
        :meth:`~GenericGraph.shortest_path_lengths` and
        :meth:`~GenericGraph.shortest_paths`

        INPUT:

        - ``v`` -- either a single vertex or a list of vertices. If it is not
          specified, then it is taken to be all vertices

        - ``by_weight`` -- boolean (default: ``False``); if ``True``, edge
          weights are taken into account; if ``False``, all edges have weight 1

        - ``algorithm`` -- string (default: ``None``); one of the following
          algorithms:

          - ``'BFS'`` -- the computation is done through a BFS centered on each
            vertex successively. Works only if ``by_weight==False``

          - ``'Floyd-Warshall-Cython'`` -- a Cython implementation of the
            Floyd-Warshall algorithm. Works only if ``by_weight==False`` and
            ``v is None`` or ``v`` should contain all vertices of ``self``.

          - ``'Floyd-Warshall-Python'`` -- a Python implementation of the
            Floyd-Warshall algorithm. Works also with weighted graphs, even with
            negative weights (but no negative cycle is allowed). However, ``v``
            must be ``None`` or ``v`` should contain all vertices of ``self``.

          - ``'Dijkstra_NetworkX'`` -- the Dijkstra algorithm, implemented in
            NetworkX. It works with weighted graphs, but no negative weight is
            allowed.

          - ``'Dijkstra_Boost'`` -- the Dijkstra algorithm, implemented in Boost
            (works only with positive weights)

          - ``'Johnson_Boost'`` -- the Johnson algorithm, implemented in
            Boost (works also with negative weights, if there is no negative
            cycle). Works only if ``v is None`` or ``v`` should contain all
            vertices of ``self``.

          - ``'From_Dictionary'`` -- uses the (already computed) distances, that
            are provided by input variable ``dist_dict``

          - ``None`` (default): Sage chooses the best algorithm:
            ``'From_Dictionary'`` if ``dist_dict`` is not None, ``'BFS'`` for
            unweighted graphs, ``'Dijkstra_Boost'`` if all weights are
            positive, ``'Johnson_Boost'`` otherwise.

        - ``weight_function`` -- function (default: ``None``); a function that
          takes as input an edge ``(u, v, l)`` and outputs its weight. If not
          ``None``, ``by_weight`` is automatically set to ``True``. If ``None``
          and ``by_weight`` is ``True``, we use the edge label ``l``, if ``l``
          is not ``None``, else ``1`` as a weight.

        - ``check_weight`` -- boolean (default: ``True``); if ``True``, we check
          that the ``weight_function`` outputs a number for each edge

        - ``dist_dict`` -- dictionary (default: ``None``); a dict of dicts of
          distances (used only if ``algorithm=='From_Dictionary'``)

        - ``with_labels`` -- boolean (default: ``False``); whether to return a
          list or a dictionary keyed by vertices

        EXAMPLES::

            sage: G = graphs.KrackhardtKiteGraph().to_directed()
            sage: G.eccentricity()
            [4, 4, 4, 4, 4, 3, 3, 2, 3, 4]
            sage: G.vertices(sort=True)
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            sage: G.eccentricity(7)
            2
            sage: G.eccentricity([7,8,9])
            [2, 3, 4]
            sage: G.eccentricity([7,8,9], with_labels=True) == {8: 3, 9: 4, 7: 2}
            True
            sage: G = DiGraph(3)
            sage: G.eccentricity(with_labels=True)
            {0: +Infinity, 1: +Infinity, 2: +Infinity}
            sage: G = DiGraph({0:[]})
            sage: G.eccentricity(with_labels=True)
            {0: 0}
            sage: G = DiGraph([(0,1,2), (1,2,3), (2,0,2)])
            sage: G.eccentricity(algorithm='BFS')
            [2, 2, 2]
            sage: G.eccentricity(algorithm='Floyd-Warshall-Cython')
            [2, 2, 2]
            sage: G.eccentricity(by_weight=True, algorithm='Dijkstra_NetworkX')         # needs networkx
            [5, 5, 4]
            sage: G.eccentricity(by_weight=True, algorithm='Dijkstra_Boost')
            [5, 5, 4]
            sage: G.eccentricity(by_weight=True, algorithm='Johnson_Boost')
            [5, 5, 4]
            sage: G.eccentricity(by_weight=True, algorithm='Floyd-Warshall-Python')
            [5, 5, 4]
            sage: G.eccentricity(dist_dict=G.shortest_path_all_pairs(by_weight=True)[0])
            [5, 5, 4]

        TESTS:

        A non-implemented algorithm::

            sage: G.eccentricity(algorithm='boh')
            Traceback (most recent call last):
            ...
            ValueError: unknown algorithm "boh"

        An algorithm that does not work with edge weights::

            sage: G.eccentricity(by_weight=True, algorithm='BFS')
            Traceback (most recent call last):
            ...
            ValueError: algorithm 'BFS' does not work with weights
            sage: G.eccentricity(by_weight=True, algorithm='Floyd-Warshall-Cython')
            Traceback (most recent call last):
            ...
            ValueError: algorithm 'Floyd-Warshall-Cython' does not work with weights

        An algorithm that computes the all-pair-shortest-paths when not all
        vertices are needed::

            sage: G.eccentricity(0, algorithm='Floyd-Warshall-Cython')
            Traceback (most recent call last):
            ...
            ValueError: algorithm 'Floyd-Warshall-Cython' works only if all eccentricities are needed
            sage: G.eccentricity(0, algorithm='Floyd-Warshall-Python')
            Traceback (most recent call last):
            ...
            ValueError: algorithm 'Floyd-Warshall-Python' works only if all eccentricities are needed
            sage: G.eccentricity(0, algorithm='Johnson_Boost')
            Traceback (most recent call last):
            ...
            ValueError: algorithm 'Johnson_Boost' works only if all eccentricities are needed
        """
        by_weight, weight_function = self._get_weight_function(by_weight=by_weight,
                                                               weight_function=weight_function,
                                                               check_weight=check_weight)

        if not by_weight:
            # We don't want the default weight function
            weight_function = None
        elif algorithm in ['BFS', 'Floyd-Warshall-Cython']:
            raise ValueError("algorithm '{}' does not work with weights".format(algorithm))
        if algorithm is None:
            if dist_dict is not None:
                algorithm = 'From_Dictionary'
            elif not by_weight:
                algorithm = 'BFS'
            elif any(float(weight_function(e)) < 0 for e in self.edge_iterator()):
                algorithm = 'Johnson_Boost'
            if algorithm is None:
                algorithm = 'Dijkstra_Boost'

        if v is not None:
            if not isinstance(v, list):
                v = [v]
            v_set = set(v)

        if v is None or all(u in v_set for u in self):
            if v is None:
                v = list(self)

            # If we want to use BFS, we use the Cython routine
            if algorithm == 'BFS':
                from sage.graphs.distances_all_pairs import eccentricity
                algo = 'standard'
                if with_labels:
                    return dict(zip(v, eccentricity(self, algorithm=algo, vertex_list=v)))
                else:
                    return eccentricity(self, algorithm=algo, vertex_list=v)

            if algorithm in ['Floyd-Warshall-Python', 'Floyd-Warshall-Cython', 'Johnson_Boost']:
                dist_dict = self.shortest_path_all_pairs(by_weight=by_weight, algorithm=algorithm,
                                                         weight_function=weight_function,
                                                         check_weight=False)[0]
                algorithm = 'From_Dictionary'

        elif algorithm in ['Floyd-Warshall-Python', 'Floyd-Warshall-Cython', 'Johnson_Boost']:
            raise ValueError("algorithm '" + algorithm + "' works only if all" +
                             " eccentricities are needed")

        ecc = {}

        from sage.rings.infinity import Infinity

        for u in v:
            if algorithm == 'From_Dictionary':
                length = dist_dict[u]
            else:
                # If algorithm is wrong, the error is raised by the
                # shortest_path_lengths function
                length = self.shortest_path_lengths(u, by_weight=by_weight,
                                                    algorithm=algorithm,
                                                    weight_function=weight_function,
                                                    check_weight=False)

            if len(length) != self.num_verts():
                ecc[u] = Infinity
            else:
                ecc[u] = max(length.values())

        if with_labels:
            return ecc
        if len(ecc) == 1:
            # return single value
            v, = ecc.values()
            return v
        return [ecc[u] for u in v]

    def radius(self, by_weight=False, algorithm=None, weight_function=None,
               check_weight=True):
        r"""
        Return the radius of the DiGraph.

        The radius is defined to be the minimum eccentricity of any vertex,
        where the eccentricity is the maximum distance to any other
        vertex. For more information and examples on how to use input variables,
        see :meth:`~GenericGraph.shortest_paths` and
        :meth:`~DiGraph.eccentricity`

        INPUT:

        - ``by_weight`` -- boolean (default: ``False``); if ``True``, edge
          weights are taken into account; if ``False``, all edges have weight 1

        - ``algorithm`` -- string (default: ``None``); see method
          :meth:`eccentricity` for the list of available algorithms

        - ``weight_function`` -- function (default: ``None``); a function that
          takes as input an edge ``(u, v, l)`` and outputs its weight. If not
          ``None``, ``by_weight`` is automatically set to ``True``. If ``None``
          and ``by_weight`` is ``True``, we use the edge label ``l``, if ``l``
          is not ``None``, else ``1`` as a weight.

        - ``check_weight`` -- boolean (default: ``True``); if ``True``, we check
          that the ``weight_function`` outputs a number for each edge

        EXAMPLES:

        The more symmetric a DiGraph is, the smaller (diameter - radius) is::

            sage: G = graphs.BarbellGraph(9, 3).to_directed()
            sage: G.radius()
            3
            sage: G.diameter()
            6

        ::

            sage: G = digraphs.Circuit(9)
            sage: G.radius()
            8
            sage: G.diameter()
            8

        TESTS::

            sage: G = DiGraph()
            sage: G.radius()
            Traceback (most recent call last):
            ...
            ValueError: radius is not defined for the empty DiGraph

        Check that :issue:`35300` is fixed::

            sage: H = DiGraph([[42, 'John'], [(42, 'John')]])
            sage: H.radius()
            1
        """
        if not self.order():
            raise ValueError("radius is not defined for the empty DiGraph")

        return min(self.eccentricity(v=list(self), by_weight=by_weight,
                                     weight_function=weight_function,
                                     check_weight=check_weight,
                                     algorithm=algorithm))

    def diameter(self, by_weight=False, algorithm=None, weight_function=None,
                 check_weight=True):
        r"""
        Return the diameter of the DiGraph.

        The diameter is defined to be the maximum distance between two vertices.
        It is infinite if the DiGraph is not strongly connected.

        For more information and examples on how to use input variables, see
        :meth:`~GenericGraph.shortest_paths` and
        :meth:`~DiGraph.eccentricity`

        INPUT:

        - ``by_weight`` -- boolean (default: ``False``); if ``True``, edge
          weights are taken into account; if ``False``, all edges have weight 1

        - ``algorithm`` -- string (default: ``None``); one of the following
          algorithms:

          - ``'BFS'``: the computation is done through a BFS centered on each
            vertex successively. Works only if ``by_weight==False``. It computes
            all the eccentricities and return the maximum value.

          - ``'Floyd-Warshall-Cython'``: a Cython implementation of the
            Floyd-Warshall algorithm. Works only if ``by_weight==False``. It
            computes all the eccentricities and return the maximum value.

          - ``'Floyd-Warshall-Python'``: a Python implementation of the
            Floyd-Warshall algorithm. Works also with weighted graphs, even with
            negative weights (but no negative cycle is allowed). It computes all
            the eccentricities and return the maximum value.

          - ``'Dijkstra_NetworkX'``: the Dijkstra algorithm, implemented in
            NetworkX. It works with weighted graphs, but no negative weight is
            allowed. It computes all the eccentricities and return the maximum
            value.

          - ``'DiFUB'``, ``'2Dsweep'``: these algorithms are
            implemented in :func:`sage.graphs.distances_all_pairs.diameter` and
            :func:`sage.graphs.base.boost_graph.diameter`. ``'2Dsweep'`` returns
            lower bound on the diameter, while ``'DiFUB'`` returns the exact
            computed diameter. They also work with negative weight, if there is
            no negative cycle. See the functions documentation for more
            information.

          - ``'standard'`` : the standard algorithm is implemented in
            :func:`sage.graphs.distances_all_pairs.diameter`. It works only
            if ``by_weight==False``. See the function documentation for more
            information. It computes all the eccentricities and return the
            maximum value.

          - ``'Dijkstra_Boost'``: the Dijkstra algorithm, implemented in Boost
            (works only with positive weights). It computes all the
            eccentricities and return the maximum value.

          - ``'Johnson_Boost'``: the Johnson algorithm, implemented in
            Boost (works also with negative weights, if there is no negative
            cycle). It computes all the eccentricities and return the maximum
            value.

          - ``None`` (default): Sage chooses the best algorithm: ``'DiFUB'``.

        - ``weight_function`` -- function (default: ``None``); a function that
          takes as input an edge ``(u, v, l)`` and outputs its weight. If not
          ``None``, ``by_weight`` is automatically set to ``True``. If ``None``
          and ``by_weight`` is ``True``, we use the edge label ``l``, if ``l``
          is not ``None``, else ``1`` as weight.

        - ``check_weight`` -- boolean (default: ``True``); if ``True``, we check
          that the ``weight_function`` outputs a number for each edge

        EXAMPLES::

            sage: # needs sage.combinat
            sage: G = digraphs.DeBruijn(5,4)
            sage: G.diameter()
            4
            sage: G = digraphs.GeneralizedDeBruijn(9, 3)
            sage: G.diameter()
            2

        TESTS::

            sage: G = graphs.RandomGNP(40, 0.4).to_directed()
            sage: d1 = G.diameter(algorithm='DiFUB', by_weight=True)
            sage: d2 = max(G.eccentricity(algorithm='Dijkstra_Boost', by_weight=True))
            sage: d1 == d2
            True
            sage: G.diameter(algorithm='BFS', by_weight=True)
            Traceback (most recent call last):
            ...
            ValueError: algorithm 'BFS' does not work with weights
            sage: G.diameter(algorithm='Floyd-Warshall-Cython', by_weight=True)
            Traceback (most recent call last):
            ...
            ValueError: algorithm 'Floyd-Warshall-Cython' does not work with weights
            sage: G = digraphs.Path(5)
            sage: G.diameter(algorithm = 'DiFUB')
            +Infinity
            sage: G = DiGraph([(1,2,4), (2,1,7)])
            sage: G.diameter(algorithm='2Dsweep', by_weight=True)
            7.0
            sage: G.delete_edge(2,1,7)
            sage: G.add_edge(2,1,-5)
            sage: G.diameter(algorithm='2Dsweep', by_weight=True)
            Traceback (most recent call last):
            ...
            ValueError: the graph contains a negative cycle
            sage: G = DiGraph()
            sage: G.diameter()
            Traceback (most recent call last):
            ...
            ValueError: diameter is not defined for the empty DiGraph

        :issue:`32095` is fixed::

            sage: g6 = 'guQOUOQCW[IaDBCVP_IE\\RfxV@WMSaeHgheEIA@tfOJkB~@EpGLCrs'
            sage: g6 += 'aPIpwgQI_`Abs_x?VWxNJAo@w\\hffCDAW]bYGMIZGC_PYOrIw[Gp['
            sage: g6 += '@FTgc_O}E?fXAnGCB{gSaUcD'
            sage: G = Graph(g6).to_directed()
            sage: G.diameter(algorithm='DiFUB', by_weight=False)
            3
            sage: G.diameter(algorithm='DiFUB', by_weight=True)
            3.0

        Check that :issue:`35300` is fixed::

            sage: H = DiGraph([[42, 'John'], [(42, 'John')]])
            sage: H.diameter()
            +Infinity
            sage: H.add_edge('John', 42)
            sage: H.diameter()
            1
        """
        if not self.order():
            raise ValueError("diameter is not defined for the empty DiGraph")

        by_weight, weight_function = self._get_weight_function(by_weight=by_weight,
                                                               weight_function=weight_function,
                                                               check_weight=check_weight)

        if not by_weight:
            # We don't want the default weight function
            weight_function = None
        elif algorithm in ['BFS', 'Floyd-Warshall-Cython']:
            raise ValueError("algorithm '{}' does not work with weights".format(algorithm))

        if algorithm is None:
            algorithm = 'DiFUB'

        if algorithm in ['2Dsweep', 'DiFUB']:
            if not by_weight:
                from sage.graphs.distances_all_pairs import diameter
                return diameter(self, algorithm=algorithm)
            else:
                from sage.graphs.base.boost_graph import diameter
                return diameter(self, algorithm=algorithm,
                                weight_function=weight_function,
                                check_weight=False)

        if algorithm == 'BFS':
            from sage.graphs.distances_all_pairs import diameter
            return diameter(self, algorithm='standard')

        return max(self.eccentricity(v=list(self), by_weight=by_weight,
                                     weight_function=weight_function,
                                     check_weight=False,
                                     algorithm=algorithm))

    def center(self, by_weight=False, algorithm=None, weight_function=None,
               check_weight=True):
        r"""
        Return the set of vertices in the center of the DiGraph.

        The center is the set of vertices whose eccentricity is equal to the
        radius of the DiGraph, i.e., achieving the minimum eccentricity.

        For more information and examples on how to use input variables,
        see :meth:`~GenericGraph.shortest_paths` and
        :meth:`~DiGraph.eccentricity`

        INPUT:

        - ``by_weight`` -- boolean (default: ``False``); if ``True``, edge
          weights are taken into account; if ``False``, all edges have weight 1

        - ``algorithm`` -- string (default: ``None``); see method
          :meth:`eccentricity` for the list of available algorithms

        - ``weight_function`` -- function (default: ``None``); a function that
          takes as input an edge ``(u, v, l)`` and outputs its weight. If not
          ``None``, ``by_weight`` is automatically set to ``True``. If ``None``
          and ``by_weight`` is ``True``, we use the edge label ``l`` as a
          weight, if ``l`` is not ``None``, else ``1`` as a weight.

        - ``check_weight`` -- boolean (default: ``True``); if ``True``, we check
          that the ``weight_function`` outputs a number for each edge

        EXAMPLES:

        Every vertex is a center in a Circuit-DiGraph::

            sage: G = digraphs.Circuit(9)
            sage: G.center()
            [0, 1, 2, 3, 4, 5, 6, 7, 8]

        Center can be the whole graph::

            sage: G.subgraph(G.center()) == G
            True

        Some other graphs::

            sage: G = digraphs.Path(5)
            sage: G.center()
            [0]
            sage: G = DiGraph([(0,1,2), (1,2,3), (2,0,2)])
            sage: G.center(by_weight=True)
            [2]

        TESTS::

            sage: G = DiGraph()
            sage: G.center()
            []
            sage: G = DiGraph(3)
            sage: G.center()
            [0, 1, 2]
        """
        ecc = self.eccentricity(v=list(self), by_weight=by_weight,
                                weight_function=weight_function,
                                algorithm=algorithm,
                                check_weight=check_weight,
                                with_labels=True)
        try:
            r = min(ecc.values())
        except Exception:
            return []
        return [v for v in self if ecc[v] == r]

    def periphery(self, by_weight=False, algorithm=None, weight_function=None,
                  check_weight=True):
        r"""
        Return the set of vertices in the periphery of the DiGraph.

        The periphery is the set of vertices whose eccentricity is equal to the
        diameter of the DiGraph, i.e., achieving the maximum eccentricity.

        For more information and examples on how to use input variables,
        see :meth:`~GenericGraph.shortest_paths` and
        :meth:`~DiGraph.eccentricity`

        INPUT:

        - ``by_weight`` -- boolean (default: ``False``); if ``True``, edge
          weights are taken into account; if ``False``, all edges have weight 1

        - ``algorithm`` -- string (default: ``None``); see method
          :meth:`eccentricity` for the list of available algorithms

        - ``weight_function`` -- function (default: ``None``); a function that
          takes as input an edge ``(u, v, l)`` and outputs its weight. If not
          ``None``, ``by_weight`` is automatically set to ``True``. If ``None``
          and ``by_weight`` is ``True``, we use the edge label ``l`` as a
          weight, if ``l`` is not ``None``, else ``1`` as a weight.

        - ``check_weight`` -- boolean (default: ``True``); if ``True``, we check
          that the ``weight_function`` outputs a number for each edge

        EXAMPLES::

            sage: G = graphs.DiamondGraph().to_directed()
            sage: G.periphery()
            [0, 3]
            sage: P = digraphs.Path(5)
            sage: P.periphery()
            [1, 2, 3, 4]
            sage: G = digraphs.Complete(5)
            sage: G.subgraph(G.periphery()) == G
            True

        TESTS::

            sage: G = DiGraph()
            sage: G.periphery()
            []
            sage: G.add_vertex()
            0
            sage: G.periphery()
            [0]
        """
        ecc = self.eccentricity(v=list(self), by_weight=by_weight,
                                weight_function=weight_function,
                                algorithm=algorithm,
                                check_weight=check_weight,
                                with_labels=True)
        try:
            d = max(ecc.values())
        except Exception:
            return []
        return [v for v in self if ecc[v] == d]

    # Paths and cycles iterators

    def _all_cycles_iterator_vertex(self, vertex, starting_vertices=None, simple=False,
                                    rooted=False, max_length=None, trivial=False,
                                    remove_acyclic_edges=True):
        r"""
        Return an iterator over the cycles of ``self`` starting with the given
        vertex.

        INPUT:

        - ``vertex`` -- the starting vertex of the cycle

        - ``starting_vertices`` -- iterable (default: ``None``); vertices from
          which the cycles must start. If ``None``, then all vertices of the
          graph can be starting points. This argument is necessary if ``rooted``
          is set to ``True``.

        - ``simple`` -- boolean (default: ``False``); if set to ``True``, then
          only simple cycles are considered. A cycle is simple if the only
          vertex occurring twice in it is the starting and ending one.

        - ``rooted`` -- boolean (default: ``False``); if set to False, then
          cycles differing only by their starting vertex are considered the same
          (e.g. ``['a', 'b', 'c', 'a']`` and ``['b', 'c', 'a',
          'b']``). Otherwise, all cycles are enumerated.

        - ``max_length`` -- nonnegative integer (default: ``None``); the
          maximum length of the enumerated paths. If set to ``None``, then all
          lengths are allowed.

        - ``trivial`` -- boolean (default: ``False``); if set to ``True``, then
          the empty paths are also enumerated

        - ``remove_acyclic_edges`` -- boolean (default: ``True``); whether
          acyclic edges must be removed from the graph.  Used to avoid
          recomputing it for each vertex

        OUTPUT: iterator

        EXAMPLES::

            sage: g = DiGraph({'a': ['a', 'b'], 'b': ['c'], 'c': ['d'], 'd': ['c']}, loops=True)
            sage: it = g._all_cycles_iterator_vertex('a', simple=False, max_length=None)
            sage: for i in range(5): print(next(it))
            ['a', 'a']
            ['a', 'a', 'a']
            ['a', 'a', 'a', 'a']
            ['a', 'a', 'a', 'a', 'a']
            ['a', 'a', 'a', 'a', 'a', 'a']
            sage: it = g._all_cycles_iterator_vertex('c', simple=False, max_length=None)
            sage: for i in range(5): print(next(it))
            ['c', 'd', 'c']
            ['c', 'd', 'c', 'd', 'c']
            ['c', 'd', 'c', 'd', 'c', 'd', 'c']
            ['c', 'd', 'c', 'd', 'c', 'd', 'c', 'd', 'c']
            ['c', 'd', 'c', 'd', 'c', 'd', 'c', 'd', 'c', 'd', 'c']

            sage: it = g._all_cycles_iterator_vertex('d', simple=False, max_length=None)
            sage: for i in range(5): print(next(it))
            ['d', 'c', 'd']
            ['d', 'c', 'd', 'c', 'd']
            ['d', 'c', 'd', 'c', 'd', 'c', 'd']
            ['d', 'c', 'd', 'c', 'd', 'c', 'd', 'c', 'd']
            ['d', 'c', 'd', 'c', 'd', 'c', 'd', 'c', 'd', 'c', 'd']

        It is possible to set a maximum length so that the number of cycles is
        finite::

            sage: it = g._all_cycles_iterator_vertex('d', simple=False, max_length=6)
            sage: list(it)
            [['d', 'c', 'd'], ['d', 'c', 'd', 'c', 'd'], ['d', 'c', 'd', 'c', 'd', 'c', 'd']]

        When ``simple`` is set to True, the number of cycles is finite since no vertex
        but the first one can occur more than once::

            sage: it = g._all_cycles_iterator_vertex('d', simple=True, max_length=None)
            sage: list(it)
            [['d', 'c', 'd']]

        By default, the empty cycle is not enumerated::

            sage: it = g._all_cycles_iterator_vertex('d', simple=True, trivial=True)
            sage: list(it)
            [['d'], ['d', 'c', 'd']]
        """
        if starting_vertices is None:
            starting_vertices = [vertex]
        # First enumerate the empty cycle
        if trivial:
            yield [vertex]
        # First we remove vertices and edges that are not part of any cycle
        if remove_acyclic_edges:
            sccs = self.strongly_connected_components()
            d = {}
            for id, component in enumerate(sccs):
                for v in component:
                    d[v] = id
            h = copy(self)
            h.delete_edges((u, v) for u, v in h.edge_iterator(labels=False) if d[u] != d[v])
        else:
            h = self
        queue = [[vertex]]
        if max_length is None:
            from sage.rings.infinity import Infinity
            max_length = Infinity
        while queue:
            path = queue.pop(0)
            # Checks if a cycle has been found
            if len(path) > 1 and path[0] == path[-1]:
                yield path
            # Makes sure that the current cycle is not too long
            # Also if a cycle has been encountered and only simple cycles are
            # allowed, Then it discards the current path
            if len(path) <= max_length and (not simple or path.count(path[-1]) == 1):
                for neighbor in h.neighbor_out_iterator(path[-1]):
                    # If cycles are not rooted, makes sure to keep only the
                    # minimum cycle according to the lexicographic order
                    if rooted or neighbor not in starting_vertices or path[0] <= neighbor:
                        queue.append(path + [neighbor])

    def all_cycles_iterator(self, starting_vertices=None, simple=False,
                            rooted=False, max_length=None, trivial=False):
        r"""
        Return an iterator over all the cycles of ``self`` starting with one of
        the given vertices.

        The cycles are enumerated in increasing length order.

        INPUT:

        - ``starting_vertices`` -- iterable (default: ``None``); vertices from
          which the cycles must start. If ``None``, then all vertices of the
          graph can be starting points. This argument is necessary if ``rooted``
          is set to ``True``.

        - ``simple`` -- boolean (default: ``False``); if set to ``True``, then
          only simple cycles are considered. A cycle is simple if the only
          vertex occurring twice in it is the starting and ending one.

        - ``rooted`` -- boolean (default: ``False``); if set to False, then
          cycles differing only by their starting vertex are considered the same
          (e.g. ``['a', 'b', 'c', 'a']`` and ``['b', 'c', 'a',
          'b']``). Otherwise, all cycles are enumerated.

        - ``max_length`` -- nonnegative integer (default: ``None``); the
          maximum length of the enumerated paths. If set to ``None``, then all
          lengths are allowed.

        - ``trivial`` -- boolean (default: ``False``); if set to ``True``, then
          the empty paths are also enumerated

        OUTPUT: iterator

        .. SEEALSO::

            - :meth:`all_simple_cycles`

        AUTHOR:

            Alexandre Blondin Masse

        EXAMPLES::

            sage: g = DiGraph({'a': ['a', 'b'], 'b': ['c'], 'c': ['d'], 'd': ['c']}, loops=True)
            sage: it = g.all_cycles_iterator()
            sage: for _ in range(7): print(next(it))
            ['a', 'a']
            ['a', 'a', 'a']
            ['c', 'd', 'c']
            ['a', 'a', 'a', 'a']
            ['a', 'a', 'a', 'a', 'a']
            ['c', 'd', 'c', 'd', 'c']
            ['a', 'a', 'a', 'a', 'a', 'a']

        There are no cycles in the empty graph and in acyclic graphs::

            sage: g = DiGraph()
            sage: it = g.all_cycles_iterator()
            sage: list(it)
            []
            sage: g = DiGraph({0:[1]})
            sage: it = g.all_cycles_iterator()
            sage: list(it)
            []

        It is possible to restrict the starting vertices of the cycles::

            sage: g = DiGraph({'a': ['a', 'b'], 'b': ['c'], 'c': ['d'], 'd': ['c']}, loops=True)
            sage: it = g.all_cycles_iterator(starting_vertices=['b', 'c'])
            sage: for _ in range(3): print(next(it))
            ['c', 'd', 'c']
            ['c', 'd', 'c', 'd', 'c']
            ['c', 'd', 'c', 'd', 'c', 'd', 'c']

        Also, one can bound the length of the cycles::

            sage: it = g.all_cycles_iterator(max_length=3)
            sage: list(it)
            [['a', 'a'], ['a', 'a', 'a'], ['c', 'd', 'c'],
             ['a', 'a', 'a', 'a']]

        By default, cycles differing only by their starting point are not all
        enumerated, but this may be parametrized::

            sage: it = g.all_cycles_iterator(max_length=3, rooted=False)
            sage: list(it)
            [['a', 'a'], ['a', 'a', 'a'], ['c', 'd', 'c'],
             ['a', 'a', 'a', 'a']]
            sage: it = g.all_cycles_iterator(max_length=3, rooted=True)
            sage: list(it)
            [['a', 'a'], ['a', 'a', 'a'], ['c', 'd', 'c'], ['d', 'c', 'd'],
             ['a', 'a', 'a', 'a']]

        One may prefer to enumerate simple cycles, i.e. cycles such that the only
        vertex occurring twice in it is the starting and ending one (see also
        :meth:`all_simple_cycles`)::

            sage: it = g.all_cycles_iterator(simple=True)
            sage: list(it)
            [['a', 'a'], ['c', 'd', 'c']]
            sage: g = digraphs.Circuit(4)
            sage: list(g.all_cycles_iterator(simple=True))
            [[0, 1, 2, 3, 0]]
        """
        if starting_vertices is None:
            starting_vertices = self
        # Since a cycle is always included in a given strongly connected
        # component, we may remove edges from the graph
        sccs = self.strongly_connected_components()
        d = {}
        for id, component in enumerate(sccs):
            for v in component:
                d[v] = id
        h = copy(self)
        h.delete_edges((u, v) for u, v in h.edge_iterator(labels=False) if d[u] != d[v])

        # We create one cycles iterator per vertex. This is necessary if we
        # want to iterate over cycles with increasing length.
        def cycle_iter(v):
            return h._all_cycles_iterator_vertex(v,
                                                 starting_vertices=starting_vertices,
                                                 simple=simple,
                                                 rooted=rooted,
                                                 max_length=max_length,
                                                 trivial=trivial,
                                                 remove_acyclic_edges=False)

        vertex_iterators = {v: cycle_iter(v) for v in starting_vertices}

        cycles = []
        for vi in vertex_iterators.values():
            try:
                cycle = next(vi)
                cycles.append((len(cycle), cycle))
            except StopIteration:
                pass
        # Since we always extract a shortest path, using a heap
        # can speed up the algorithm
        from heapq import heapify, heappop, heappush
        heapify(cycles)
        while cycles:
            # We choose the shortest available cycle
            _, shortest_cycle = heappop(cycles)
            yield shortest_cycle
            # We update the cycle iterator to its next available cycle if it
            # exists
            try:
                cycle = next(vertex_iterators[shortest_cycle[0]])
                heappush(cycles, (len(cycle), cycle))
            except StopIteration:
                pass

    def all_simple_cycles(self, starting_vertices=None, rooted=False,
                          max_length=None, trivial=False):
        r"""
        Return a list of all simple cycles of ``self``.

        INPUT:

        - ``starting_vertices`` -- iterable (default: ``None``); vertices from
          which the cycles must start. If ``None``, then all vertices of the
          graph can be starting points. This argument is necessary if ``rooted``
          is set to ``True``.

        - ``rooted`` -- boolean (default: ``False``); if set to False, then
          cycles differing only by their starting vertex are considered the same
          (e.g. ``['a', 'b', 'c', 'a']`` and ``['b', 'c', 'a',
          'b']``). Otherwise, all cycles are enumerated.

        - ``max_length`` -- nonnegative integer (default: ``None``); the
          maximum length of the enumerated paths. If set to ``None``, then all
          lengths are allowed.

        - ``trivial`` -- boolean (default: ``False``); if set to ``True``, then
          the empty paths are also enumerated

        OUTPUT: list

        .. NOTE::

            Although the number of simple cycles of a finite graph is always
            finite, computing all its cycles may take a very long time.

        EXAMPLES::

            sage: g = DiGraph({'a': ['a', 'b'], 'b': ['c'], 'c': ['d'], 'd': ['c']}, loops=True)
            sage: g.all_simple_cycles()
            [['a', 'a'], ['c', 'd', 'c']]

        The directed version of the Petersen graph::

            sage: g = graphs.PetersenGraph().to_directed()
            sage: g.all_simple_cycles(max_length=4)
            [[0, 1, 0], [0, 4, 0], [0, 5, 0], [1, 2, 1], [1, 6, 1], [2, 3, 2],
             [2, 7, 2], [3, 4, 3], [3, 8, 3], [4, 9, 4], [5, 7, 5], [5, 8, 5],
             [6, 8, 6], [6, 9, 6], [7, 9, 7]]
            sage: g.all_simple_cycles(max_length=6)
            [[0, 1, 0], [0, 4, 0], [0, 5, 0], [1, 2, 1], [1, 6, 1], [2, 3, 2],
             [2, 7, 2], [3, 4, 3], [3, 8, 3], [4, 9, 4], [5, 7, 5], [5, 8, 5],
             [6, 8, 6], [6, 9, 6], [7, 9, 7], [0, 1, 2, 3, 4, 0],
             [0, 1, 2, 7, 5, 0], [0, 1, 6, 8, 5, 0], [0, 1, 6, 9, 4, 0],
             [0, 4, 3, 2, 1, 0], [0, 4, 3, 8, 5, 0], [0, 4, 9, 6, 1, 0],
             [0, 4, 9, 7, 5, 0], [0, 5, 7, 2, 1, 0], [0, 5, 7, 9, 4, 0],
             [0, 5, 8, 3, 4, 0], [0, 5, 8, 6, 1, 0], [1, 2, 3, 8, 6, 1],
             [1, 2, 7, 9, 6, 1], [1, 6, 8, 3, 2, 1], [1, 6, 9, 7, 2, 1],
             [2, 3, 4, 9, 7, 2], [2, 3, 8, 5, 7, 2], [2, 7, 5, 8, 3, 2],
             [2, 7, 9, 4, 3, 2], [3, 4, 9, 6, 8, 3], [3, 8, 6, 9, 4, 3],
             [5, 7, 9, 6, 8, 5], [5, 8, 6, 9, 7, 5], [0, 1, 2, 3, 8, 5, 0],
             [0, 1, 2, 7, 9, 4, 0], [0, 1, 6, 8, 3, 4, 0],
             [0, 1, 6, 9, 7, 5, 0], [0, 4, 3, 2, 7, 5, 0],
             [0, 4, 3, 8, 6, 1, 0], [0, 4, 9, 6, 8, 5, 0],
             [0, 4, 9, 7, 2, 1, 0], [0, 5, 7, 2, 3, 4, 0],
             [0, 5, 7, 9, 6, 1, 0], [0, 5, 8, 3, 2, 1, 0],
             [0, 5, 8, 6, 9, 4, 0], [1, 2, 3, 4, 9, 6, 1],
             [1, 2, 7, 5, 8, 6, 1], [1, 6, 8, 5, 7, 2, 1],
             [1, 6, 9, 4, 3, 2, 1], [2, 3, 8, 6, 9, 7, 2],
             [2, 7, 9, 6, 8, 3, 2], [3, 4, 9, 7, 5, 8, 3],
             [3, 8, 5, 7, 9, 4, 3]]

        The complete graph (without loops) on `4` vertices::

            sage: g = graphs.CompleteGraph(4).to_directed()
            sage: g.all_simple_cycles()
            [[0, 1, 0], [0, 2, 0], [0, 3, 0], [1, 2, 1], [1, 3, 1], [2, 3, 2],
             [0, 1, 2, 0], [0, 1, 3, 0], [0, 2, 1, 0], [0, 2, 3, 0],
             [0, 3, 1, 0], [0, 3, 2, 0], [1, 2, 3, 1], [1, 3, 2, 1],
             [0, 1, 2, 3, 0], [0, 1, 3, 2, 0], [0, 2, 1, 3, 0],
             [0, 2, 3, 1, 0], [0, 3, 1, 2, 0], [0, 3, 2, 1, 0]]

        If the graph contains a large number of cycles, one can bound the length
        of the cycles, or simply restrict the possible starting vertices of the
        cycles::

            sage: g = graphs.CompleteGraph(20).to_directed()
            sage: g.all_simple_cycles(max_length=2)
            [[0, 16, 0], [0, 1, 0], [0, 17, 0], [0, 2, 0], [0, 18, 0],
             [0, 3, 0], [0, 19, 0], [0, 4, 0], [0, 5, 0], [0, 6, 0], [0, 7, 0],
             [0, 8, 0], [0, 9, 0], [0, 10, 0], [0, 11, 0], [0, 12, 0],
             [0, 13, 0], [0, 14, 0], [0, 15, 0], [1, 16, 1], [1, 17, 1],
             [1, 2, 1], [1, 18, 1], [1, 3, 1], [1, 19, 1], [1, 4, 1], [1, 5, 1],
             [1, 6, 1], [1, 7, 1], [1, 8, 1], [1, 9, 1], [1, 10, 1], [1, 11, 1],
             [1, 12, 1], [1, 13, 1], [1, 14, 1], [1, 15, 1], [2, 16, 2],
             [2, 17, 2], [2, 18, 2], [2, 3, 2], [2, 19, 2], [2, 4, 2],
             [2, 5, 2], [2, 6, 2], [2, 7, 2], [2, 8, 2], [2, 9, 2], [2, 10, 2],
             [2, 11, 2], [2, 12, 2], [2, 13, 2], [2, 14, 2], [2, 15, 2],
             [3, 16, 3], [3, 17, 3], [3, 18, 3], [3, 19, 3], [3, 4, 3],
             [3, 5, 3], [3, 6, 3], [3, 7, 3], [3, 8, 3], [3, 9, 3], [3, 10, 3],
             [3, 11, 3], [3, 12, 3], [3, 13, 3], [3, 14, 3], [3, 15, 3],
             [4, 16, 4], [4, 17, 4], [4, 18, 4], [4, 19, 4], [4, 5, 4],
             [4, 6, 4], [4, 7, 4], [4, 8, 4], [4, 9, 4], [4, 10, 4], [4, 11, 4],
             [4, 12, 4], [4, 13, 4], [4, 14, 4], [4, 15, 4], [5, 16, 5],
             [5, 17, 5], [5, 18, 5], [5, 19, 5], [5, 6, 5], [5, 7, 5],
             [5, 8, 5], [5, 9, 5], [5, 10, 5], [5, 11, 5], [5, 12, 5],
             [5, 13, 5], [5, 14, 5], [5, 15, 5], [6, 16, 6], [6, 17, 6],
             [6, 18, 6], [6, 19, 6], [6, 7, 6], [6, 8, 6], [6, 9, 6],
             [6, 10, 6], [6, 11, 6], [6, 12, 6], [6, 13, 6], [6, 14, 6],
             [6, 15, 6], [7, 16, 7], [7, 17, 7], [7, 18, 7], [7, 19, 7],
             [7, 8, 7], [7, 9, 7], [7, 10, 7], [7, 11, 7], [7, 12, 7],
             [7, 13, 7], [7, 14, 7], [7, 15, 7], [8, 16, 8], [8, 17, 8],
             [8, 18, 8], [8, 19, 8], [8, 9, 8], [8, 10, 8], [8, 11, 8],
             [8, 12, 8], [8, 13, 8], [8, 14, 8], [8, 15, 8], [9, 16, 9],
             [9, 17, 9], [9, 18, 9], [9, 19, 9], [9, 10, 9], [9, 11, 9],
             [9, 12, 9], [9, 13, 9], [9, 14, 9], [9, 15, 9], [10, 16, 10],
             [10, 17, 10], [10, 18, 10], [10, 19, 10], [10, 11, 10],
             [10, 12, 10], [10, 13, 10], [10, 14, 10], [10, 15, 10],
             [11, 16, 11], [11, 17, 11], [11, 18, 11], [11, 19, 11],
             [11, 12, 11], [11, 13, 11], [11, 14, 11], [11, 15, 11],
             [12, 16, 12], [12, 17, 12], [12, 18, 12], [12, 19, 12],
             [12, 13, 12], [12, 14, 12], [12, 15, 12], [13, 16, 13],
             [13, 17, 13], [13, 18, 13], [13, 19, 13], [13, 14, 13],
             [13, 15, 13], [14, 16, 14], [14, 17, 14], [14, 18, 14],
             [14, 19, 14], [14, 15, 14], [15, 16, 15], [15, 17, 15],
             [15, 18, 15], [15, 19, 15], [16, 17, 16], [16, 18, 16],
             [16, 19, 16], [17, 18, 17], [17, 19, 17], [18, 19, 18]]

            sage: g = graphs.CompleteGraph(20).to_directed()
            sage: g.all_simple_cycles(max_length=2, starting_vertices=[0])
            [[0, 16, 0], [0, 1, 0], [0, 17, 0], [0, 2, 0], [0, 18, 0],
             [0, 3, 0], [0, 19, 0], [0, 4, 0], [0, 5, 0], [0, 6, 0],
             [0, 7, 0], [0, 8, 0], [0, 9, 0], [0, 10, 0], [0, 11, 0],
             [0, 12, 0], [0, 13, 0], [0, 14, 0], [0, 15, 0]]

        One may prefer to distinguish equivalent cycles having distinct starting
        vertices (compare the following examples)::

            sage: g = graphs.CompleteGraph(4).to_directed()
            sage: g.all_simple_cycles(max_length=2, rooted=False)
            [[0, 1, 0], [0, 2, 0], [0, 3, 0], [1, 2, 1], [1, 3, 1], [2, 3, 2]]
            sage: g.all_simple_cycles(max_length=2, rooted=True)
            [[0, 1, 0], [0, 2, 0], [0, 3, 0], [1, 0, 1], [1, 2, 1], [1, 3, 1],
             [2, 0, 2], [2, 1, 2], [2, 3, 2], [3, 0, 3], [3, 1, 3], [3, 2, 3]]
        """
        return list(self.all_cycles_iterator(starting_vertices=starting_vertices,
                                             simple=True, rooted=rooted,
                                             max_length=max_length, trivial=trivial))

    def path_semigroup(self):
        """
        The partial semigroup formed by the paths of this quiver.

        EXAMPLES::

            sage: Q = DiGraph({1: {2: ['a', 'c']}, 2: {3: ['b']}})
            sage: F = Q.path_semigroup(); F                                             # needs sage.libs.flint
            Partial semigroup formed by the directed paths of Multi-digraph on 3 vertices
            sage: list(F)                                                               # needs sage.libs.flint
            [e_1, e_2, e_3, a, c, b, a*b, c*b]
        """
        from sage.quivers.path_semigroup import PathSemigroup
        return PathSemigroup(self)

    def auslander_reiten_quiver(self):
        r"""
        Return the Auslander-Reiten quiver of ``self``.

        .. SEEALSO::

            :class:`~sage.quivers.ar_quiver.AuslanderReitenQuiver`

        EXAMPLES::

            sage: D = DiGraph([[1,2,'a'], [1,2,'b']], multiedges=True)
            sage: D.auslander_reiten_quiver()
            Auslander-Reiten quiver of Multi-digraph on 2 vertices
        """
        from sage.quivers.ar_quiver import AuslanderReitenQuiver
        return AuslanderReitenQuiver(self)

    # Directed Acyclic Graphs (DAGs)

    def topological_sort(self, implementation='default'):
        """
        Return a topological sort of the digraph if it is acyclic.

        If the digraph contains a directed cycle, a :exc:`TypeError`
        is raised. As topological sorts are not necessarily unique,
        different implementations may yield different results.

        A topological sort is an ordering of the vertices of the digraph such
        that each vertex comes before all of its successors. That is, if `u`
        comes before `v` in the sort, then there may be a directed path from `u`
        to `v`, but there will be no directed path from `v` to `u`.

        INPUT:

        - ``implementation`` -- string (default: ``'default'``); either use the
          default Cython implementation, or the default NetworkX library
          (``implementation = "NetworkX"``)

        .. SEEALSO::

            - :meth:`is_directed_acyclic` -- tests whether a directed graph is
              acyclic (can also join a certificate; a topological sort or a
              circuit in the graph)

        EXAMPLES::

            sage: D = DiGraph({0: [1, 2, 3], 4: [2, 5], 1: [8], 2: [7], 3: [7],
            ....:   5: [6, 7], 7: [8], 6: [9], 8: [10], 9: [10]})
            sage: D.plot(layout='circular').show()                                      # needs sage.plot
            sage: D.topological_sort()
            [4, 5, 6, 9, 0, 1, 2, 3, 7, 8, 10]

        ::

            sage: D.add_edge(9, 7)
            sage: D.topological_sort()
            [4, 5, 6, 9, 0, 1, 2, 3, 7, 8, 10]

        Using the NetworkX implementation ::

            sage: s = list(D.topological_sort(implementation="NetworkX")); s  # random  # needs networkx
            [0, 4, 1, 3, 2, 5, 6, 9, 7, 8, 10]
            sage: all(s.index(u) < s.index(v)                                           # needs networkx
            ....:     for u, v in D.edges(sort=False, labels=False))
            True

        ::

            sage: D.add_edge(7, 4)
            sage: D.topological_sort()
            Traceback (most recent call last):
            ...
            TypeError: digraph is not acyclic; there is no topological sort

        TESTS:

        A wrong value for the ``implementation`` keyword::

            sage: D.topological_sort(implementation = "cloud-reading")
            Traceback (most recent call last):
            ...
            ValueError: implementation must be set to one of "default" or "NetworkX"
        """
        if implementation == "default":
            b, ordering = self._backend.is_directed_acyclic(certificate=True)
            if b:
                return ordering
            else:
                raise TypeError('digraph is not acyclic; there is no topological sort')

        elif implementation == "NetworkX":
            import networkx
            S = networkx.topological_sort(self.networkx_graph())
            if S is None:
                raise TypeError('digraph is not acyclic; there is no topological sort')
            else:
                return S

        raise ValueError("implementation must be set to one of \"default\" or \"NetworkX\"")

    def topological_sort_generator(self):
        """
        Return an iterator over all topological sorts of the digraph if
        it is acyclic.

        If the digraph contains a directed cycle, a :exc:`TypeError`
        is raised.

        A topological sort is an ordering of the vertices of the digraph such
        that each vertex comes before all of its successors. That is, if u comes
        before v in the sort, then there may be a directed path from u to v, but
        there will be no directed path from v to u. See also
        :meth:`topological_sort`.

        AUTHORS:

        - Mike Hansen - original implementation

        - Robert L. Miller: wrapping, documentation

        REFERENCE:

        - [1] Pruesse, Gara and Ruskey, Frank. Generating Linear
          Extensions Fast. SIAM J. Comput., Vol. 23 (1994), no. 2, pp.
          373-386.

        EXAMPLES::

            sage: D = DiGraph({0: [1, 2], 1: [3], 2: [3, 4]})
            sage: D.plot(layout='circular').show()                                      # needs sage.plot
            sage: list(D.topological_sort_generator())                                  # needs sage.modules sage.rings.finite_rings
            [[0, 1, 2, 3, 4], [0, 2, 1, 3, 4], [0, 2, 1, 4, 3],
             [0, 2, 4, 1, 3], [0, 1, 2, 4, 3]]

        ::

            sage: for sort in D.topological_sort_generator():                           # needs sage.modules sage.rings.finite_rings
            ....:     for u, v in D.edge_iterator(labels=False):
            ....:         if sort.index(u) > sort.index(v):
            ....:             print("this should never happen")
        """
        from sage.combinat.posets.posets import Poset
        return Poset(self).linear_extensions()

    # Visualization

    def layout_acyclic(self, rankdir='up', **options):
        """
        Return a ranked layout so that all edges point upward.

        To this end, the heights of the vertices are set according to the level
        set decomposition of the graph (see :meth:`.level_sets`).

        This is achieved by calling ``graphviz`` and ``dot2tex`` if available
        (see :meth:`.layout_graphviz`), and using a spring layout with fixed
        vertical placement of the vertices otherwise (see
        :meth:`.layout_acyclic_dummy` and
        :meth:`~sage.graphs.generic_graph.GenericGraph.layout_ranked`).

        Non acyclic graphs are partially supported by ``graphviz``, which then
        chooses some edges to point down.

        INPUT:

        - ``rankdir`` -- string (default: ``'up'``); indicates which direction
          the edges should point toward among ``'up'``, ``'down'``, ``'left'``,
          or ``'right'``

        - ``**options`` -- passed down to
          :meth:`~sage.graphs.generic_graph.GenericGraph.layout_ranked` or
          :meth:`~sage.graphs.generic_graph.GenericGraph.layout_graphviz`

        EXAMPLES::

            sage: H = DiGraph({0: [1, 2], 1: [3], 2: [3], 3: [], 5: [1, 6], 6: [2, 3]})

        The actual layout computed depends on whether dot2tex and graphviz are
        installed, so we don't test its relative values::

            sage: H.layout_acyclic()
            {0: [..., ...], 1: [..., ...], 2: [..., ...], 3: [..., ...], 5: [..., ...], 6: [..., ...]}

            sage: H = DiGraph({0: [1]})
            sage: pos = H.layout_acyclic(rankdir='up')
            sage: pos[1][1] > pos[0][1] + .5
            True
            sage: pos = H.layout_acyclic(rankdir='down')
            sage: pos[1][1] < pos[0][1] - .5
            True
            sage: pos = H.layout_acyclic(rankdir='right')
            sage: pos[1][0] > pos[0][0] + .5
            True
            sage: pos = H.layout_acyclic(rankdir='left')
            sage: pos[1][0] < pos[0][0] - .5
            True
        """
        if have_dot2tex():
            return self.layout_graphviz(rankdir=rankdir, **options)
        return self.layout_acyclic_dummy(rankdir=rankdir, **options)

    def layout_acyclic_dummy(self, heights=None, rankdir='up', **options):
        """
        Return a ranked layout so that all edges point upward.

        To this end, the heights of the vertices are set according to the level
        set decomposition of the graph (see :meth:`level_sets`). This is
        achieved by a spring layout with fixed vertical placement of the
        vertices otherwise (see :meth:`layout_acyclic_dummy` and
        :meth:`~sage.graphs.generic_graph.GenericGraph.layout_ranked`).

        INPUT:

        - ``rankdir`` -- string (default: ``'up'``); indicates which direction
          the edges should point toward among ``'up'``, ``'down'``, ``'left'``,
          or ``'right'``

        - ``**options`` -- passed down to
          :meth:`~sage.graphs.generic_graph.GenericGraph.layout_ranked`

        EXAMPLES::

            sage: H = DiGraph({0: [1, 2], 1: [3], 2: [3], 3: [], 5: [1, 6], 6: [2, 3]})
            sage: H.layout_acyclic_dummy()
            {0: [1.0..., 0], 1: [1.0..., 1], 2: [1.5..., 2], 3: [1.5..., 3], 5: [2.0..., 0], 6: [2.0..., 1]}

            sage: H = DiGraph({0: [1]})
            sage: H.layout_acyclic_dummy(rankdir='up')
            {0: [0.5..., 0], 1: [0.5..., 1]}
            sage: H.layout_acyclic_dummy(rankdir='down')
            {0: [0.5..., 1], 1: [0.5..., 0]}
            sage: H.layout_acyclic_dummy(rankdir='left')
            {0: [1, 0.5...], 1: [0, 0.5...]}
            sage: H.layout_acyclic_dummy(rankdir='right')
            {0: [0, 0.5...], 1: [1, 0.5...]}
            sage: H = DiGraph({0: [1, 2], 1: [3], 2: [3], 3: [1], 5: [1, 6], 6: [2, 3]})
            sage: H.layout_acyclic_dummy()
            Traceback (most recent call last):
            ...
            ValueError: `self` should be an acyclic graph

        TESTS:

        :issue:`31681` is fixed::

            sage: H = DiGraph({0: [1], 'X': [1]}, format='dict_of_lists')
            sage: pos = H.layout_acyclic_dummy(rankdir='up')
            sage: pos['X'][1] == 0 and pos[0][1] == 0
            True
            sage: pos[1][1] == 1
            True
        """
        if heights is None:
            if not self.is_directed_acyclic():
                raise ValueError("`self` should be an acyclic graph")
            levels = self.level_sets()
            # Sort vertices in each level in best effort mode
            for i in range(len(levels)):
                try:
                    label = sorted(levels[i])
                    levels[i] = label
                except TypeError:
                    continue
            if rankdir == 'down' or rankdir == 'left':
                levels.reverse()
            heights = {i: levels[i] for i in range(len(levels))}
        positions = self.layout_ranked(heights=heights, **options)
        if rankdir == 'left' or rankdir == 'right':
            for coordinates in positions.values():
                coordinates.reverse()
        return positions

    def level_sets(self):
        r"""
        Return the level set decomposition of the digraph.

        OUTPUT: list of non empty lists of vertices of this graph

        The level set decomposition of the digraph is a list `l` such that the
        level `l[i]` contains all the vertices having all their predecessors in
        the levels `l[j]` for `j < i`, and at least one in level `l[i-1]`
        (unless `i = 0`).

        The level decomposition contains exactly the vertices not occurring in
        any cycle of the graph. In particular, the graph is acyclic if and only
        if the decomposition forms a set partition of its vertices, and we
        recover the usual level set decomposition of the corresponding poset.

        EXAMPLES::

            sage: H = DiGraph({0: [1, 2], 1: [3], 2: [3], 3: [], 5: [1, 6], 6: [2, 3]})
            sage: H.level_sets()
            [[0, 5], [1, 6], [2], [3]]

            sage: H = DiGraph({0: [1, 2], 1: [3], 2: [3], 3: [1], 5: [1, 6], 6: [2, 3]})
            sage: H.level_sets()
            [[0, 5], [6], [2]]

        This routine is mostly used for Hasse diagrams of posets::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0: [1, 2], 1: [3], 2: [3], 3: []})
            sage: [len(x) for x in H.level_sets()]
            [1, 2, 1]

        ::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0: [1, 2], 1: [3], 2: [4], 3: [4]})
            sage: [len(x) for x in H.level_sets()]
            [1, 2, 1, 1]

        Complexity: `O(n+m)` in time and `O(n)` in memory (besides the storage
        of the graph itself), where `n` and `m` are respectively the number of
        vertices and edges (assuming that appending to a list is constant time,
        which it is not quite).
        """
        in_degrees = self.in_degree(labels=True)
        level = [x for x in in_degrees if not in_degrees[x]]
        Levels = []
        while level:
            Levels.append(level)
            new_level = []
            for x in level:
                for y in self.neighbor_out_iterator(x):
                    in_degrees[y] -= 1
                    if not in_degrees[y]:
                        new_level.append(y)
            level = new_level
        return Levels

    def is_aperiodic(self):
        r"""
        Return whether the current ``DiGraph`` is aperiodic.

        A directed graph is aperiodic if there is no integer `k > 1` that
        divides the length of every cycle in the graph. See the
        :wikipedia:`Aperiodic_graph` for more information.

        EXAMPLES:

        The following graph has period ``2``, so it is not aperiodic::

            sage: g = DiGraph({0: [1], 1: [0]})
            sage: g.is_aperiodic()
            False

        The following graph has a cycle of length 2 and a cycle of length 3,
        so it is aperiodic::

            sage: g = DiGraph({0: [1, 4], 1: [2], 2: [0], 4: [0]})
            sage: g.is_aperiodic()
            True

        .. SEEALSO::

            :meth:`period`
        """
        return self.period() == 1

    def period(self):
        r"""
        Return the period of the current ``DiGraph``.

        The period of a directed graph is the largest integer that divides the
        length of every cycle in the graph. See the :wikipedia:`Aperiodic_graph`
        for more information.

        EXAMPLES:

        The following graph has period ``2``::

            sage: g = DiGraph({0: [1], 1: [0]})
            sage: g.period()
            2

        The following graph has a cycle of length 2 and a cycle of length 3,
        so it has period ``1``::

            sage: g = DiGraph({0: [1, 4], 1: [2], 2: [0], 4: [0]})
            sage: g.period()
            1

        Here is an example of computing the period of a digraph which is not
        strongly connected. By definition, it is the :func:`gcd` of the periods
        of its strongly connected components::

            sage: g = DiGraph({-1: [-2], -2: [-3], -3: [-1],
            ....:     1: [2], 2: [1]})
            sage: g.period()
            1
            sage: sorted([s.period() for s
            ....:         in g.strongly_connected_components_subgraphs()])
            [2, 3]

        ALGORITHM:

        See the networkX implementation of ``is_aperiodic``, that is based on
        breadth first search.

        .. SEEALSO::

            :meth:`is_aperiodic`
        """
        from sage.arith.misc import GCD as gcd

        g = 0

        for component in self.strongly_connected_components():
            levels = {s: None for s in component}
            vertices_in_scc = levels  # considers level as a set
            s = component[0]
            levels[s] = 0
            this_level = [s]
            idx = 1
            while this_level:
                next_level = []
                for u in this_level:
                    # we have levels[u] == idx - 1
                    for v in self.neighbor_out_iterator(u):
                        # ignore edges leaving the component
                        if v not in vertices_in_scc:
                            continue
                        level_v = levels[v]
                        if level_v is not None:  # Non-Tree Edge
                            g = gcd(g, idx - level_v)
                            if g == 1:
                                return 1
                        else:  # Tree Edge
                            next_level.append(v)
                            levels[v] = idx
                this_level = next_level
                idx += 1

        return g

    def flow_polytope(self, edges=None, ends=None, backend=None):
        r"""
        Return the flow polytope of a digraph.

        The flow polytope of a directed graph is the polytope consisting of all
        nonnegative flows on the graph with a given set `S` of sources and a
        given set `T` of sinks.

        A *flow* on a directed graph `G` with a given set `S` of sources and a
        given set `T` of sinks means an assignment of a nonnegative real to each
        edge of `G` such that the flow is conserved in each vertex outside of
        `S` and `T`, and there is a unit of flow entering each vertex in `S` and
        a unit of flow leaving each vertex in `T`. These flows clearly form a
        polytope in the space of all assignments of reals to the edges of `G`.

        The polytope is empty unless the sets `S` and `T` are equinumerous.

        By default, `S` is taken to be the set of all sources (i.e., vertices of
        indegree `0`) of `G`, and `T` is taken to be the set of all sinks (i.e.,
        vertices of outdegree `0`) of `G`. If a different choice of `S` and `T`
        is desired, it can be specified using the optional ``ends`` parameter.

        The polytope is returned as a polytope in `\RR^m`, where `m` is the
        number of edges of the digraph ``self``. The `k`-th coordinate of a
        point in the polytope is the real assigned to the `k`-th edge of
        ``self``. The order of the edges is the one returned by
        ``self.edges(sort=True)``. If a different order is desired, it can be specified
        using the optional ``edges`` parameter.

        The faces and volume of these polytopes are of interest. Examples of
        these polytopes are the Chan-Robbins-Yuen polytope and the
        Pitman-Stanley polytope [PS2002]_.

        INPUT:

        - ``edges`` -- list (default: ``None``); a list of edges of ``self``. If
          not specified, the list of all edges of ``self`` is used with the
          default ordering of ``self.edges(sort=True)``. This determines which coordinate
          of a point in the polytope will correspond to which edge of
          ``self``. It is also possible to specify a list which contains not all
          edges of ``self``; this results in a polytope corresponding to the
          flows which are `0` on all remaining edges. Notice that the edges
          entered here must be in the precisely same format as outputted by
          ``self.edges()``; so, if ``self.edges()`` outputs an edge in the form
          ``(1, 3, None)``, then ``(1, 3)`` will not do!

        - ``ends`` -- (default: ``(self.sources(), self.sinks())``) a
          pair `(S, T)` of an iterable `S` and an iterable `T`

        - ``backend`` -- string or ``None`` (default); the backend to use;
          see :meth:`sage.geometry.polyhedron.constructor.Polyhedron`

        .. NOTE::

            Flow polytopes can also be built through the ``polytopes.<tab>``
            object::

                sage: polytopes.flow_polytope(digraphs.Path(5))                         # needs sage.geometry.polyhedron
                A 0-dimensional polyhedron in QQ^4 defined as the convex hull of 1 vertex

        EXAMPLES:

        A commutative square::

            sage: G = DiGraph({1: [2, 3], 2: [4], 3: [4]})
            sage: fl = G.flow_polytope(); fl                                            # needs sage.geometry.polyhedron
            A 1-dimensional polyhedron in QQ^4 defined as the convex hull
            of 2 vertices
            sage: fl.vertices()                                                         # needs sage.geometry.polyhedron
            (A vertex at (0, 1, 0, 1), A vertex at (1, 0, 1, 0))

        Using a different order for the edges of the graph::

            sage: ordered_edges = G.edges(sort=True, key=lambda x: x[0] - x[1])
            sage: fl = G.flow_polytope(edges=ordered_edges); fl                         # needs sage.geometry.polyhedron
            A 1-dimensional polyhedron in QQ^4 defined as the convex hull of 2 vertices
            sage: fl.vertices()                                                         # needs sage.geometry.polyhedron
            (A vertex at (0, 1, 1, 0), A vertex at (1, 0, 0, 1))

        A tournament on 4 vertices::

            sage: H = digraphs.TransitiveTournament(4)
            sage: fl = H.flow_polytope(); fl                                            # needs sage.geometry.polyhedron
            A 3-dimensional polyhedron in QQ^6 defined as the convex hull
            of 4 vertices
            sage: fl.vertices()                                                         # needs sage.geometry.polyhedron
            (A vertex at (0, 0, 1, 0, 0, 0),
             A vertex at (0, 1, 0, 0, 0, 1),
             A vertex at (1, 0, 0, 0, 1, 0),
             A vertex at (1, 0, 0, 1, 0, 1))

        Restricting to a subset of the edges::

            sage: fl = H.flow_polytope(edges=[(0, 1, None), (1, 2, None),               # needs sage.geometry.polyhedron
            ....:                             (2, 3, None), (0, 3, None)]); fl
            A 1-dimensional polyhedron in QQ^4 defined as the convex hull
            of 2 vertices
            sage: fl.vertices()                                                         # needs sage.geometry.polyhedron
            (A vertex at (0, 0, 0, 1), A vertex at (1, 1, 1, 0))

        Using a different choice of sources and sinks::

            sage: # needs sage.geometry.polyhedron
            sage: fl = H.flow_polytope(ends=([1], [3])); fl
            A 1-dimensional polyhedron in QQ^6 defined as the convex hull
            of 2 vertices
            sage: fl.vertices()
            (A vertex at (0, 0, 0, 1, 0, 1), A vertex at (0, 0, 0, 0, 1, 0))
            sage: fl = H.flow_polytope(ends=([0, 1], [3])); fl
            The empty polyhedron in QQ^6
            sage: fl = H.flow_polytope(ends=([3], [0])); fl
            The empty polyhedron in QQ^6
            sage: fl = H.flow_polytope(ends=([0, 1], [2, 3])); fl
            A 3-dimensional polyhedron in QQ^6 defined as the convex hull
            of 5 vertices
            sage: fl.vertices()
            (A vertex at (0, 0, 1, 1, 0, 0),
             A vertex at (0, 1, 0, 0, 1, 0),
             A vertex at (1, 0, 0, 2, 0, 1),
             A vertex at (1, 0, 0, 1, 1, 0),
             A vertex at (0, 1, 0, 1, 0, 1))
            sage: fl = H.flow_polytope(edges=[(0, 1, None), (1, 2, None),
            ....:                             (2, 3, None), (0, 2, None),
            ....:                             (1, 3, None)],
            ....:                      ends=([0, 1], [2, 3])); fl
            A 2-dimensional polyhedron in QQ^5 defined as the convex hull
            of 4 vertices
            sage: fl.vertices()
            (A vertex at (0, 0, 0, 1, 1),
             A vertex at (1, 2, 1, 0, 0),
             A vertex at (1, 1, 0, 0, 1),
             A vertex at (0, 1, 1, 1, 0))

        A digraph with one source and two sinks::

            sage: Y = DiGraph({1: [2], 2: [3, 4]})
            sage: Y.flow_polytope()                                                     # needs sage.geometry.polyhedron
            The empty polyhedron in QQ^3

        A digraph with one vertex and no edge::

            sage: Z = DiGraph({1: []})
            sage: Z.flow_polytope()                                                     # needs sage.geometry.polyhedron
            A 0-dimensional polyhedron in QQ^0 defined as the convex hull
            of 1 vertex

        A digraph with multiple edges (:issue:`28837`)::

            sage: G = DiGraph([(0, 1), (0,1)], multiedges=True); G
            Multi-digraph on 2 vertices
            sage: P = G.flow_polytope(); P                                              # needs sage.geometry.polyhedron
            A 1-dimensional polyhedron in QQ^2 defined as the convex hull of 2 vertices
            sage: P.vertices()                                                          # needs sage.geometry.polyhedron
            (A vertex at (1, 0), A vertex at (0, 1))
            sage: P.lines()                                                             # needs sage.geometry.polyhedron
            ()
        """
        from sage.geometry.polyhedron.constructor import Polyhedron
        if edges is None:
            edges = self.edges(sort=False)
        m = len(edges)
        ineqs = [[0] * (i + 1) + [1] + [0] * (m - i - 1) for i in range(m)]

        eqs = []
        for u in self:
            ins = set(self.incoming_edge_iterator(u))
            outs = set(self.outgoing_edge_iterator(u))
            eq = [Integer(j in ins) - Integer(j in outs) for j in edges]

            const = 0
            if ends is None:
                if not ins:  # sources (indegree 0)
                    const += 1
                if not outs:  # sinks (outdegree 0)
                    const -= 1
            else:
                if u in ends[0]:  # chosen sources
                    const += 1
                if u in ends[1]:  # chosen sinks
                    const -= 1

            eq = [const] + eq
            eqs.append(eq)

        return Polyhedron(ieqs=ineqs, eqns=eqs, backend=backend)

    def is_tournament(self):
        r"""
        Check whether the digraph is a tournament.

        A tournament is a digraph in which each pair of distinct vertices is
        connected by a single arc.

        EXAMPLES::

            sage: g = digraphs.RandomTournament(6)
            sage: g.is_tournament()
            True
            sage: u,v = next(g.edge_iterator(labels=False))
            sage: g.add_edge(v, u)
            sage: g.is_tournament()
            False
            sage: g.add_edges([(u, v), (v, u)])
            sage: g.is_tournament()
            False

        .. SEEALSO::

            - :wikipedia:`Tournament_(graph_theory)`
            - :meth:`~sage.graphs.digraph_generators.DiGraphGenerators.RandomTournament`
            - :meth:`~sage.graphs.digraph_generators.DiGraphGenerators.TransitiveTournament`
        """
        self._scream_if_not_simple()

        if self.size() != self.order() * (self.order() - 1) // 2:
            return False

        import itertools
        return not any(self.has_edge(u, v) == self.has_edge(v, u)
                       for u, v in itertools.combinations(self, 2))

    def _girth_bfs(self, odd=False, certificate=False):
        r"""
        Return the girth of the digraph using breadth-first search.

        Loops are ignored, so the returned value is at least 2.

        INPUT:

        - ``odd`` -- boolean (default: ``False``); whether to compute the odd
          girth

        - ``certificate`` -- boolean (default: ``False``); whether to return
          ``(g, c)``, where ``g`` is the (odd) girth and ``c`` is a list
          of vertices of a directed cycle of length ``g`` in the graph,
          thus providing a certificate that the (odd) girth is at most ``g``,
          or ``None`` if ``g`` is infinite

        EXAMPLES:

        A digraph with girth 4 and odd girth 5::

            sage: G = DiGraph([(0, 1), (1, 2), (1, 3), (2, 3), (3, 4), (4, 0)])
            sage: G._girth_bfs(certificate=True)  # random
            (4, [1, 3, 4, 0])
            sage: G._girth_bfs(odd=True)
            5

        .. SEEALSO::

            * :meth:`~sage.graphs.GenericGraph.girth` -- return the girth of the
              graph
            * :meth:`~sage.graphs.GenericGraph.odd_girth` -- return the odd
              girth of the graph
        """
        n = self.num_verts()
        best = n + 1
        seen = set()
        for w in self:
            seen.add(w)
            inSpan, outSpan = {w: None}, {w: None}
            depth = 1
            outList, inList = set([w]), set([w])
            while 2 * depth <= best:
                nextOutList, nextInList = set(), set()
                for v in outList:
                    for u in self.neighbor_out_iterator(v):
                        if u in seen:
                            continue
                        if u not in outSpan:
                            outSpan[u] = v
                            nextOutList.add(u)
                        if u in inList:
                            best = depth * 2 - 1
                            ends = (v, u)
                            bestSpans = (outSpan, inSpan)
                            break
                    if best == 2 * depth - 1:
                        break
                if best == 2 * depth - 1:
                    break
                for v in inList:
                    for u in self.neighbor_in_iterator(v):
                        if u in seen:
                            continue
                        if u not in inSpan:
                            inSpan[u] = v
                            nextInList.add(u)
                        if not odd and u in nextOutList:
                            best = depth * 2
                            ends = (u, v)
                            bestSpans = (outSpan, inSpan)
                            break
                    if best == 2 * depth:
                        break
                if best <= 2:
                    break
                outList = nextOutList
                inList = nextInList
                depth += 1
        if best == n + 1:
            from sage.rings.infinity import Infinity
            return (Infinity, None) if certificate else Infinity
        if certificate:
            cycles = {}
            for x, span in zip(ends, bestSpans):
                cycles[x] = []
                y = x
                while span[y] is not None:
                    cycles[x].append(y)
                    y = span[y]
            cycles[x].append(y)
            u, v = ends
            return (best, list(reversed(cycles[u])) + cycles[v])
        return best

    def out_branchings(self, source, spanning=True):
        r"""
        Return an iterator over the out branchings rooted at given vertex in
        ``self``.

        An out-branching is a directed tree rooted at ``source`` whose arcs are
        directed from source to leaves. An out-branching is spanning if it
        contains all vertices of the digraph.

        If no spanning out branching rooted at ``source`` exist, raises
        :exc:`ValueError` or return non spanning out branching rooted at
        ``source``, depending on the value of ``spanning``.

        INPUT:

        - ``source`` -- vertex used as the source for all out branchings

        - ``spanning`` -- boolean (default: ``True``); if ``False`` return
          maximum out branching from ``source``. Otherwise, return spanning out
          branching if exists.

        OUTPUT: an iterator over the out branchings rooted in the given source

        .. SEEALSO::

            - :meth:`~sage.graphs.digraph.DiGraph.in_branchings`
              -- iterator over in-branchings rooted at given vertex.
            - :meth:`~sage.graphs.graph.Graph.spanning_trees`
              -- returns all spanning trees.
            - :meth:`~sage.graphs.generic_graph.GenericGraph.spanning_trees_count`
              -- counts the number of spanning trees.

        ALGORITHM:

        Recursively computes all out branchings.

        At each step:

            0. clean the graph (see below)
            1. pick an edge e out of source
            2. find all out branchings that do not contain e by first
               removing it
            3. find all out branchings that do contain e by first
               merging the end vertices of e

        Cleaning the graph implies to remove loops and replace multiedges by a
        single one with an appropriate label since these lead to similar steps
        of computation.

        EXAMPLES:

        A bidirectional 4-cycle::

            sage: G = DiGraph({1:[2,3], 2:[1,4], 3:[1,4], 4:[2,3]}, format='dict_of_lists')
            sage: list(G.out_branchings(1))
            [Digraph on 4 vertices,
             Digraph on 4 vertices,
             Digraph on 4 vertices,
             Digraph on 4 vertices]

        With the Petersen graph turned into a symmetric directed graph::

            sage: G = graphs.PetersenGraph().to_directed()
            sage: len(list(G.out_branchings(0)))
            2000

        With a non connected ``DiGraph`` and ``spanning = True``::

            sage: G = graphs.PetersenGraph().to_directed() + graphs.PetersenGraph().to_directed()
            sage: G.out_branchings(0, spanning=True)
            Traceback (most recent call last):
            ...
            ValueError: no spanning out branching from vertex (0) exist

        With a non connected ``DiGraph`` and ``spanning = False``::

            sage: g=DiGraph([(0,1), (0,1), (1,2), (3,4)],multiedges=True)
            sage: list(g.out_branchings(0, spanning=False))
            [Digraph on 3 vertices, Digraph on 3 vertices]

        With multiedges::

            sage: G = DiGraph({0:[1,1,1], 1:[2,2]}, format='dict_of_lists', multiedges=True)
            sage: len(list(G.out_branchings(0)))
            6

        With a DiGraph already being a spanning out branching::

            sage: G = DiGraph({0:[1,2], 1:[3,4], 2:[5], 3:[], 4:[], 5:[]}, format='dict_of_lists')
            sage: next(G.out_branchings(0)) == G
            True

        TESTS:

        The empty ``DiGraph``::

            sage: G = DiGraph()
            sage: G.out_branchings(0)
            Traceback (most recent call last):
            ...
            ValueError: vertex (0) is not a vertex of the digraph

            sage: edges = [(0,0,'x'), (0,0,'y')]
            sage: G = DiGraph(edges, multiedges=True, loops=True, weighted=True)
            sage: list(G.out_branchings(0))
            [Digraph on 1 vertex]

            sage: edges = [(0,1,'x'), (0,1,'y'), (1,2,'z'), (2,0,'w')]
            sage: G = DiGraph(edges, multiedges=True, loops=True, weighted=True)
            sage: len(list(G.out_branchings(0)))
            2
        """
        def _rec_out_branchings(depth):
            r"""
            The recursive function used to enumerate out branchings.

            This function makes use of the following to keep track of partial
            out branchings:

            - ``list_edges`` -- list of edges in self
            - ``list_merged_edges`` -- list of edges that are currently merged
            - ``graph`` -- a copy of self where edges have an appropriate label
            """
            if not depth:
                # We have enough merged edges to form a out_branching
                # We iterate over the lists of labels in list_merged_edges and
                # yield the corresponding out_branchings
                for indexes in product(*list_merged_edges):
                    yield DiGraph([list_edges[index] for index in indexes],
                                  format='list_of_edges', pos=self.get_pos())

            # 1) Clean the graph
            # delete loops on source if any
            D.delete_edges(D.incoming_edge_iterator(source))

            # merge multi-edges if any by concatenating their labels
            if D.has_multiple_edges():
                merged_multiple_edges = {}
                for u, v, label in D.multiple_edges():
                    D.delete_edge(u, v, label)
                    if (u, v) not in merged_multiple_edges:
                        merged_multiple_edges[(u, v)] = label
                    else:
                        merged_multiple_edges[(u, v)] += label
                D.add_edges([(u, v, label) for (u, v), label in merged_multiple_edges.items()])

            # 2) Pick an edge e outgoing from the source
            try:
                s, x, label = next(D.outgoing_edge_iterator(source))
            except StopIteration:
                return
            # 3) Find all out_branchings that do not contain e
            # by first removing it
            D.delete_edge(s, x, label)
            if len(list(D.depth_first_search(source))) == depth + 1:
                for out_branch in _rec_out_branchings(depth):
                    yield out_branch
            D.add_edge(s, x, label)

            # 4) Find all out_branchings that do contain e by merging
            # the end vertices of e
            # store different edges to unmerged the end vertices of e
            saved_edges = D.outgoing_edges(source)
            saved_edges.remove((s, x, label))
            saved_edges += D.outgoing_edges(x)
            saved_edges += D.incoming_edges(x)

            D.merge_vertices((source, x))

            list_merged_edges.add(label)

            for out_branch in _rec_out_branchings(depth - 1):
                yield out_branch

            list_merged_edges.remove(label)

            # unmerge the end vertices of e
            D.delete_vertex(source)
            D.add_edges(saved_edges)

        def _singleton_out_branching():
            r"""
            Return a DiGraph containing only ``source`` and no edges.
            """
            D = DiGraph()
            D.add_vertex(source)
            yield D

        if not self.has_vertex(source):
            raise ValueError("vertex ({0}) is not a vertex of the digraph".format(source))

        # check if self.order == 1
        if self.order() == 1:
            return _singleton_out_branching()

        # check if the source can access to every other vertex
        if spanning:
            depth = self.order() - 1
            if len(list(self.depth_first_search(source))) < self.order():
                raise ValueError("no spanning out branching from vertex ({0}) exist".format(source))
        else:
            depth = len(list(self.depth_first_search(source))) - 1
            # if vertex is isolated
            if not depth:
                return _singleton_out_branching()

        # We build a copy of self in which each edge has a distinct label.
        # On the way, we remove loops and edges incoming to source.
        D = DiGraph(multiedges=True, loops=True)
        list_edges = list(self.edges(sort=False))
        for i, (u, v, _) in enumerate(list_edges):
            if u != v and v != source:
                D.add_edge(u, v, (i,))
        list_merged_edges = set()
        return _rec_out_branchings(depth)

    def in_branchings(self, source, spanning=True):
        r"""
        Return an iterator over the in branchings rooted at given vertex in
        ``self``.

        An in-branching is a directed tree rooted at ``source`` whose arcs are
        directed to source from leaves. An in-branching is spanning if it
        contains all vertices of the digraph.

        If no spanning in branching rooted at ``source`` exist, raises
        :exc:`ValueError` or return non spanning in branching rooted at
        ``source``, depending on the value of ``spanning``.

        INPUT:

        - ``source`` -- vertex used as the source for all in branchings

        - ``spanning`` -- boolean (default: ``True``); if ``False`` return
          maximum in branching to ``source``. Otherwise, return spanning in
          branching if exists.

        OUTPUT: an iterator over the in branchings rooted in the given source

        .. SEEALSO::

            - :meth:`~sage.graphs.digraph.DiGraph.out_branchings`
              -- iterator over out-branchings rooted at given vertex.
            - :meth:`~sage.graphs.graph.Graph.spanning_trees`
              -- returns all spanning trees.
            - :meth:`~sage.graphs.generic_graph.GenericGraph.spanning_trees_count`
              -- counts the number of spanning trees.

        ALGORITHM:

        Recursively computes all in branchings.

        At each step:

            0. clean the graph (see below)
            1. pick an edge e incoming to source
            2. find all in branchings that do not contain e by first
               removing it
            3. find all in branchings that do contain e by first
               merging the end vertices of e

        Cleaning the graph implies to remove loops and replace multiedges by a
        single one with an appropriate label since these lead to similar steps
        of computation.

        EXAMPLES:

        A bidirectional 4-cycle::

            sage: G = DiGraph({1:[2,3], 2:[1,4], 3:[1,4], 4:[2,3]}, format='dict_of_lists')
            sage: list(G.in_branchings(1))
            [Digraph on 4 vertices,
             Digraph on 4 vertices,
             Digraph on 4 vertices,
             Digraph on 4 vertices]

        With the Petersen graph turned into a symmetric directed graph::

            sage: G = graphs.PetersenGraph().to_directed()
            sage: len(list(G.in_branchings(0)))
            2000

        With a non connected ``DiGraph`` and ``spanning = True``::

            sage: G = graphs.PetersenGraph().to_directed() + graphs.PetersenGraph().to_directed()
            sage: G.in_branchings(0)
            Traceback (most recent call last):
            ...
            ValueError: no spanning in branching to vertex (0) exist

        With a non connected ``DiGraph`` and ``spanning = False``::

            sage: g=DiGraph([(1,0), (1,0), (2,1), (3,4)],multiedges=True)
            sage: list(g.in_branchings(0,spanning=False))
            [Digraph on 3 vertices, Digraph on 3 vertices]

        With multiedges::

            sage: G = DiGraph({0:[1,1,1], 1:[2,2]}, format='dict_of_lists', multiedges=True)
            sage: len(list(G.in_branchings(2)))
            6

        With a DiGraph already being a spanning in branching::

            sage: G = DiGraph({0:[], 1:[0], 2:[0], 3:[1], 4:[1], 5:[2]}, format='dict_of_lists')
            sage: next(G.in_branchings(0)) == G
            True

        TESTS:

        The empty ``DiGraph``::

            sage: G = DiGraph()
            sage: G.in_branchings(0)
            Traceback (most recent call last):
            ...
            ValueError: vertex (0) is not a vertex of the digraph

            sage: edges = [(0,0,'x'), (0,0,'y')]
            sage: G = DiGraph(edges, multiedges=True, loops=True, weighted=True)
            sage: list(G.in_branchings(0))
            [Digraph on 1 vertex]

            sage: edges = [(0,1,'x'), (0,1,'y'), (1,2,'z'), (2,0,'w')]
            sage: G = DiGraph(edges, multiedges=True, loops=True, weighted=True)
            sage: len(list(G.in_branchings(0)))
            1
        """
        def _rec_in_branchings(depth):
            r"""
            The recursive function used to enumerate in branchings.

            This function makes use of the following to keep track of partial in
            branchings:

            - ``list_edges`` -- list of edges in self
            - ``list_merged_edges`` -- list of edges that are currently merged
            - ``graph`` -- a copy of self where edges have an appropriate label
            """
            if not depth:
                # We have enough merged edges to form a in_branching
                # We iterate over the lists of labels in list_merged_edges and
                # yield the corresponding in_branchings
                for indexes in product(*list_merged_edges):
                    yield DiGraph([list_edges[index] for index in indexes],
                                  format='list_of_edges', pos=self.get_pos())

            # 1) Clean the graph
            # delete loops on source if any
            D.delete_edges(D.outgoing_edge_iterator(source))

            # merge multi-edges if any by concatenating their labels
            if D.has_multiple_edges():
                merged_multiple_edges = {}
                for u, v, label in D.multiple_edges():
                    D.delete_edge(u, v, label)
                    if (u, v) not in merged_multiple_edges:
                        merged_multiple_edges[(u, v)] = label
                    else:
                        merged_multiple_edges[(u, v)] += label
                D.add_edges([(u, v, label) for (u, v), label in merged_multiple_edges.items()])

            # 2) Pick an edge e incoming to the source
            try:
                x, s, label = next(D.incoming_edge_iterator(source))
            except StopIteration:
                return
            # 3) Find all in_branchings that do not contain e
            # by first removing it
            D.delete_edge(x, s, label)
            if len(list(D.depth_first_search(source, neighbors=D.neighbor_in_iterator))) == depth + 1:
                for in_branch in _rec_in_branchings(depth):
                    yield in_branch
            D.add_edge(x, s, label)

            # 4) Find all in_branchings that do contain e by merging
            # the end vertices of e
            # store different edges to unmerged the end vertices of e
            saved_edges = D.incoming_edges(source)
            saved_edges.remove((x, s, label))
            saved_edges += D.outgoing_edges(x)
            saved_edges += D.incoming_edges(x)

            D.merge_vertices((source, x))

            list_merged_edges.add(label)

            for in_branch in _rec_in_branchings(depth - 1):
                yield in_branch

            list_merged_edges.remove(label)

            # unmerge the end vertices of e
            D.delete_vertex(source)
            D.add_edges(saved_edges)

        def _singleton_in_branching():
            r"""
            Return a DiGraph containing only ``source`` and no edges.
            """
            D = DiGraph()
            D.add_vertex(source)
            yield D

        if not self.has_vertex(source):
            raise ValueError("vertex ({0}) is not a vertex of the digraph".format(source))

        # check if self.order == 1
        if self.order() == 1:
            return _singleton_in_branching()

        # check if the source can access to every other vertex
        if spanning:
            depth = self.order() - 1
            if len(list(self.depth_first_search(source, neighbors=self.neighbor_in_iterator))) < self.order():
                raise ValueError("no spanning in branching to vertex ({0}) exist".format(source))
        else:
            depth = len(list(self.depth_first_search(source, neighbors=self.neighbor_in_iterator))) - 1
            # if vertex is isolated
            if not depth:
                return _singleton_in_branching()

        # We build a copy of self in which each edge has a distinct label.
        # On the way, we remove loops and edges incoming to source.
        D = DiGraph(multiedges=True, loops=True)
        list_edges = list(self.edges(sort=False))
        for i, (u, v, _) in enumerate(list_edges):
            if u != v and u != source:
                D.add_edge(u, v, (i,))
        list_merged_edges = set()
        return _rec_in_branchings(depth)

    # Aliases to functions defined in other modules
    from sage.graphs.comparability import is_transitive
    from sage.graphs.base.static_sparse_graph import tarjan_strongly_connected_components as strongly_connected_components
    from sage.graphs.connectivity import is_strongly_connected
    from sage.graphs.connectivity import strongly_connected_components_digraph
    from sage.graphs.connectivity import strongly_connected_components_subgraphs
    from sage.graphs.connectivity import strongly_connected_component_containing_vertex
    from sage.graphs.connectivity import strong_articulation_points
