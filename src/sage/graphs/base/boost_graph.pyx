r"""
Interface to run Boost algorithms

Wrapper for a Boost graph. The Boost graphs are Cython C++ variables, and they
cannot be converted to Python objects: as a consequence, only functions defined
with cdef are able to create, read, modify, and delete these graphs.

A very important feature of Boost graph library is that all object are generic:
for instance, adjacency lists can be stored using different data structures,
and (most of) the functions work with all implementations provided. This feature
is implemented in our interface using fused types: however, Cython's support for
fused types is still experimental, and some features are missing. For instance,
there cannot be nested generic function calls, and no variable can have a
generic type, apart from the arguments of a generic function.

All the input functions use pointers, because otherwise we might have problems
with ``delete()``.

**Basic Boost Graph operations:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :func:`clustering_coeff`  | Return the clustering coefficient of all vertices in the graph.
    :func:`edge_connectivity` | Return the edge connectivity of the graph.
    :func:`dominator_tree` | Return a dominator tree of the graph.
    :func:`bandwidth_heuristics` | Use heuristics to approximate the bandwidth of the graph.
    :func:`min_spanning_tree` | Compute a minimum spanning tree of a (weighted) graph.
    :func:`shortest_paths` | Use Dijkstra or Bellman-Ford algorithm to compute the single-source shortest paths.
    :func:`johnson_shortest_paths` | Use Johnson algorithm to compute the all-pairs shortest paths.
    :func:`floyd_warshall_shortest_paths` | Use Floyd-Warshall algorithm to compute the all-pairs shortest paths.
    :func:`johnson_closeness_centrality` | Use Johnson algorithm to compute the closeness centrality of all vertices.
    :func:`blocks_and_cut_vertices` | Use Tarjan's algorithm to compute the blocks and cut vertices of the graph.
    :func:`min_cycle_basis` | Return a minimum weight cycle basis of the input graph.

Functions
---------
"""

#*****************************************************************************
#       Copyright (C) 2015 Michele Borassi michele.borassi@imtlucca.it
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import

cimport cython
from cysignals.signals cimport sig_check, sig_on, sig_off
from libcpp.set cimport set as cset
from libcpp.pair cimport pair

cdef boost_graph_from_sage_graph(BoostGenGraph *g, g_sage, vertex_to_int, reverse=False):
    r"""
    Initialize the Boost graph ``g`` to be equal to ``g_sage``.

    INPUT:

    - ``g`` -- a Boost graph; it must represent an empty graph (an exception is
      raised otherwise)

    - ``g_sage`` -- a Sage graph

    - ``vertex_to_int`` -- a dictionary; it is a mapping from the vertex set of
      ``g_sage`` to `(0, \ldots, n-1)`

    - ``reverse`` -- boolean (default: ``False``); when set to ``True``, the
      Boost graph is initialized with reversed edges
    """
    from sage.graphs.generic_graph import GenericGraph

    if not isinstance(g_sage, GenericGraph):
        raise TypeError("the input must be a Sage graph")

    if g.num_verts():
        raise AssertionError("the given Boost graph must be empty")

    cdef int N = g_sage.num_verts()
    cdef int i

    for i in range(N):
        g.add_vertex()

    if reverse:
        for u,v in g_sage.edge_iterator(labels=None):
            g.add_edge(vertex_to_int[v], vertex_to_int[u])
    else:
        for u,v in g_sage.edge_iterator(labels=None):
            g.add_edge(vertex_to_int[u], vertex_to_int[v])


cdef boost_weighted_graph_from_sage_graph(BoostWeightedGraph *g,
                                          g_sage,
                                          vertex_to_int,
                                          weight_function=None,
                                          reverse=False):
    r"""
    Initialize the Boost weighted graph ``g`` to be equal to ``g_sage``.

    INPUT:

    - ``g`` -- a Boost weighted graph; it must represent an empty weighted graph
      (an exception is raised otherwise)

    - ``g_sage`` -- a Sage graph

    - ``vertex_to_int`` -- a dictionary; it is a mapping from the vertex set of
      ``g_sage`` to `(0, \ldots, n-1)`

    - ``weight_function`` -- function (default: ``None``); a function which
      inputs an edge ``e`` and outputs a number. The edge weights are chosen as
      follows, and they must be convertible to floats, otherwise an error is
      raised.

      - If ``weight_function`` is not ``None``, this function is used

      - If ``weight_function`` is ``None`` and ``g`` is weighted, the edge
        labels of ``g`` are used; in other words, the weight of an edge
        ``e = (u, v, l)`` is ``l``

      - Otherwise, all weights are set to 1

    - ``reverse`` -- boolean (default: ``False``); when set to ``True``, the
      Boost graph is initialized with reversed edges
    """
    from sage.graphs.generic_graph import GenericGraph

    if not isinstance(g_sage, GenericGraph):
        raise TypeError("the input must be a Sage graph")

    if g.num_verts():
        raise AssertionError("the given Boost graph must be empty")

    cdef int N = g_sage.num_verts()
    cdef int i

    for i in range(N):
        g.add_vertex()

    if weight_function is not None:
        if reverse:
            for e in g_sage.edge_iterator():
                g.add_edge(vertex_to_int[e[1]],
                        vertex_to_int[e[0]],
                        float(weight_function(e)))
        else:
            for e in g_sage.edge_iterator():
                g.add_edge(vertex_to_int[e[0]],
                        vertex_to_int[e[1]],
                        float(weight_function(e)))
    elif g_sage.weighted():
        if reverse:
            for u,v,w in g_sage.edge_iterator():
                g.add_edge(vertex_to_int[v], vertex_to_int[u], float(w))
        else:
            for u,v,w in g_sage.edge_iterator():
                g.add_edge(vertex_to_int[u], vertex_to_int[v], float(w))
    else:
        if reverse:
            for u,v in g_sage.edge_iterator(labels=False):
                g.add_edge(vertex_to_int[v], vertex_to_int[u], 1)
        else:
            for u,v in g_sage.edge_iterator(labels=False):
                g.add_edge(vertex_to_int[u], vertex_to_int[v], 1)


cdef boost_edge_connectivity(BoostVecGenGraph *g):
    r"""
    Compute the edge connectivity of the input Boost graph.

    The output is a pair ``[ec,edges]``, where ``ec`` is the edge connectivity,
    ``edges`` is the list of edges in a minimum cut.
    """
    cdef result_ec result

    sig_on()
    result = g[0].edge_connectivity()
    sig_off()

    cdef size_t i
    edges = [(result.edges[i], result.edges[i+1])
             for i in range(0, result.edges.size(), 2)]

    return (result.ec, edges)


cpdef edge_connectivity(g):
    r"""
    Compute the edge connectivity of the input graph, using Boost.

    OUTPUT: a pair ``(ec, edges)``, where ``ec`` is the edge
    connectivity, ``edges`` is the list of edges in a minimum cut.

    .. SEEALSO::

        :meth:`sage.graphs.generic_graph.GenericGraph.edge_connectivity`

    EXAMPLES:

    Computing the edge connectivity of a clique::

        sage: from sage.graphs.base.boost_graph import edge_connectivity
        sage: g = graphs.CompleteGraph(5)
        sage: edge_connectivity(g)
        (4, [(0, 1), (0, 2), (0, 3), (0, 4)])

    Vertex-labeled graphs::

        sage: from sage.graphs.base.boost_graph import edge_connectivity
        sage: g = graphs.GridGraph([2,2])
        sage: edge_connectivity(g)  # py2
        (2, [((0, 1), (1, 1)), ((0, 1), (0, 0))])
        sage: edge_connectivity(g)  # py3
        (2, [((0, 0), (0, 1)), ((0, 0), (1, 0))])
    """
    from sage.graphs.graph import Graph
    from sage.graphs.digraph import DiGraph

    # These variables are automatically deleted when the function terminates.
    cdef BoostVecGraph g_boost_und
    cdef BoostVecDiGraph g_boost_dir
    cdef v_index i
    cdef list int_to_vertex = list(g)
    cdef dict vertex_to_int = {v: i for i, v in enumerate(int_to_vertex)}

    if isinstance(g, Graph):
        boost_graph_from_sage_graph(&g_boost_und, g, vertex_to_int)
        ec, edges = boost_edge_connectivity(&g_boost_und)

    elif isinstance(g, DiGraph):
        from sage.misc.stopgap import stopgap
        stopgap("The edge connectivity of directed graphs is not implemented " +
                "in Boost. The result may be mathematically unreliable.",18753)

        boost_graph_from_sage_graph(&g_boost_dir, g, vertex_to_int)
        ec, edges = boost_edge_connectivity(&g_boost_dir)

    else:
        raise TypeError("the input must be a Sage graph")

    return (ec, [(int_to_vertex[u], int_to_vertex[v]) for u, v in edges])


cdef boost_clustering_coeff(BoostGenGraph *g, vertices):
    r"""
    Compute the clustering coefficient of all vertices in the list provided.

    The output is a pair ``[average_clustering_coefficient, clust_of_v]``, where
    ``average_clustering_coefficient`` is the average clustering of the vertices
    in variable ``vertices``, ``clust_of_v`` is a dictionary that associates to
    each vertex (stored as an integer) its clustering coefficient.
    """
    cdef result_cc result
    cdef double result_d
    cdef v_index vi
    cdef dict clust_of_v

    if len(vertices) == g.num_verts():
        sig_on()
        result = g[0].clustering_coeff_all()
        sig_off()
        clust_of_v = {v: result.clust_of_v[v] for v in range(g.num_verts())}
        return (result.average_clustering_coefficient, clust_of_v)

    else:
        clust_of_v = {}
        for v in vertices:
            vi = v
            sig_on()
            result_d = g[0].clustering_coeff(vi)
            sig_off()
            clust_of_v[v] = result_d
        return ((sum(clust_of_v.itervalues()) / len(clust_of_v)), clust_of_v)


cpdef clustering_coeff(g, vertices=None):
    r"""
    Compute the clustering coefficient of the input graph, using Boost.

    .. SEEALSO::

        :meth:`sage.graphs.generic_graph.GenericGraph.clustering_coeff`

    INPUT:

    - ``g`` -- the input Sage Graph

    - ``vertices`` -- list (default: ``None``); the list of vertices to analyze
      (if ``None``, compute the clustering coefficient of all vertices)

    OUTPUT: a pair ``(average_clustering_coefficient, clust_of_v)``, where
    ``average_clustering_coefficient`` is the average clustering of the vertices
    in variable ``vertices``, ``clust_of_v`` is a dictionary that associates to
    each vertex its clustering coefficient. If ``vertices`` is ``None``, all
    vertices are considered.

    EXAMPLES:

    Computing the clustering coefficient of a clique::

        sage: from sage.graphs.base.boost_graph import clustering_coeff
        sage: g = graphs.CompleteGraph(5)
        sage: clustering_coeff(g)
        (1.0, {0: 1.0, 1: 1.0, 2: 1.0, 3: 1.0, 4: 1.0})
        sage: clustering_coeff(g, vertices = [0,1,2])
        (1.0, {0: 1.0, 1: 1.0, 2: 1.0})

    Of a non-clique graph with triangles::

        sage: g = graphs.IcosahedralGraph()
        sage: clustering_coeff(g, vertices=[1,2,3])
        (0.5, {1: 0.5, 2: 0.5, 3: 0.5})

    With labels::

        sage: g.relabel(list("abcdefghiklm"))
        sage: clustering_coeff(g, vertices="abde")
        (0.5, {'a': 0.5, 'b': 0.5, 'd': 0.5, 'e': 0.5})
    """
    from sage.graphs.graph import Graph

    # These variables are automatically deleted when the function terminates.
    cdef BoostVecGraph g_boost
    cdef v_index i
    cdef list g_vertices = list(g)
    cdef dict vertex_to_int = {v: i for i, v in enumerate(g_vertices)}

    if not isinstance(g, Graph):
        raise TypeError("the input must be a Sage Graph")

    boost_graph_from_sage_graph(&g_boost, g, vertex_to_int)

    if vertices is None:
        vertices = g_vertices

    cdef list vertices_boost = [vertex_to_int[v] for v in vertices]
    average_clustering, clust_v_int = boost_clustering_coeff(&g_boost, vertices_boost)
    cdef dict clust_v_sage = {g_vertices[v]: clust_v_int[v] for v in vertices_boost}
    return (average_clustering, clust_v_sage)


@cython.binding(True)
cpdef dominator_tree(g, root, return_dict=False, reverse=False):
    r"""
    Use Boost to compute the dominator tree of ``g``, rooted at ``root``.

    A node `d` dominates a node `n` if every path from the entry node
    ``root`` to `n` must go through `d`. The immediate dominator of a node
    `n` is the unique node that strictly dominates `n` but does not dominate
    any other node that dominates `n`. A dominator tree is a tree where each
    node's children are those nodes it immediately dominates. For more
    information, see the :wikipedia:`Dominator_(graph_theory)`.

    If the graph is connected and undirected, the parent of a vertex `v` is:

     - the root if `v` is in the same biconnected component as the root;

     - the first cut vertex in a path from `v` to the root, otherwise.

    If the graph is not connected, the dominator tree of the whole graph is
    equal to the dominator tree of the connected component of the root.

    If the graph is directed, computing a dominator tree is more complicated,
    and it needs time `O(m\log m)`, where `m` is the number of edges. The
    implementation provided by Boost is the most general one, so it needs time
    `O(m\log m)` even for undirected graphs.

    INPUT:

    - ``g`` -- the input Sage (Di)Graph

    - ``root`` -- the root of the dominator tree

    - ``return_dict`` -- boolean (default: ``False``); if ``True``, the function
      returns a dictionary associating to each vertex its parent in the
      dominator tree. If ``False`` (default), it returns the whole tree, as a
      ``Graph`` or a ``DiGraph``.

    - ``reverse`` -- boolean (default: ``False``); when set to ``True``,
      computes the dominator tree in the reverse graph

    OUTPUT:

    The dominator tree, as a graph or as a dictionary, depending on the
    value of ``return_dict``. If the output is a dictionary, it will contain
    ``None`` in correspondence of ``root`` and of vertices that are not
    reachable from ``root``. If the output is a graph, it will not contain
    vertices that are not reachable from ``root``.

    EXAMPLES:

    An undirected grid is biconnected, and its dominator tree is a star
    (everyone's parent is the root)::

        sage: g = graphs.GridGraph([2,2]).dominator_tree((0,0))
        sage: g.to_dictionary()
        {(0, 0): [(0, 1), (1, 0), (1, 1)], (0, 1): [(0, 0)], (1, 0): [(0, 0)], (1, 1): [(0, 0)]}

    If the graph is made by two 3-cycles `C_1,C_2` connected by an edge `(v,w)`,
    with `v \in C_1`, `w \in C_2`, the cut vertices are `v` and `w`, the
    biconnected components are `C_1`, `C_2`, and the edge `(v,w)`. If the root
    is in `C_1`, the parent of each vertex in `C_1` is the root, the parent of
    `w` is `v`, and the parent of each vertex in `C_2` is `w`::

         sage: G = 2 * graphs.CycleGraph(3)
         sage: v = 0
         sage: w = 3
         sage: G.add_edge(v,w)
         sage: G.dominator_tree(1, return_dict=True)
         {0: 1, 1: None, 2: 1, 3: 0, 4: 3, 5: 3}

    An example with a directed graph::

        sage: g = digraphs.Circuit(10).dominator_tree(5)
        sage: g.to_dictionary()
        {0: [1], 1: [2], 2: [3], 3: [4], 4: [], 5: [6], 6: [7], 7: [8], 8: [9], 9: [0]}
        sage: g = digraphs.Circuit(10).dominator_tree(5, reverse=True)
        sage: g.to_dictionary()
        {0: [9], 1: [0], 2: [1], 3: [2], 4: [3], 5: [4], 6: [], 7: [6], 8: [7], 9: [8]}

    If the output is a dictionary::

        sage: graphs.GridGraph([2,2]).dominator_tree((0,0), return_dict=True)
        {(0, 0): None, (0, 1): (0, 0), (1, 0): (0, 0), (1, 1): (0, 0)}

    TESTS:

    If ``g`` is not a graph, an error is raised::

        sage: from sage.graphs.base.boost_graph import dominator_tree
        sage: dominator_tree('I am not a graph', 0)
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage Graph or DiGraph

    If ``root`` is not a vertex, an error is raised::

        sage: digraphs.TransitiveTournament(10).dominator_tree('Not a vertex!')
        Traceback (most recent call last):
        ...
        ValueError: the input root must be a vertex of the given graph
        sage: graphs.GridGraph([2,2]).dominator_tree(0)
        Traceback (most recent call last):
        ...
        ValueError: the input root must be a vertex of the given graph
    """
    from sage.graphs.graph import Graph
    from sage.graphs.digraph import DiGraph

    if not isinstance(g, (Graph, DiGraph)):
        raise TypeError("the input must be a Sage Graph or DiGraph")
    if root not in g:
        raise ValueError("the input root must be a vertex of the given graph")

    # These variables are automatically deleted when the function terminates.
    cdef BoostVecGraph g_boost_und
    cdef BoostVecDiGraph g_boost_dir
    cdef vector[v_index] result
    cdef v_index i
    cdef list int_to_vertex = list(g)
    cdef dict vertex_to_int = {v: i for i, v in enumerate(int_to_vertex)}

    if isinstance(g, Graph):
        boost_graph_from_sage_graph(&g_boost_und, g, vertex_to_int, reverse)
        i = vertex_to_int[root]
        sig_on()
        result = g_boost_und.dominator_tree(i)
        sig_off()

    elif isinstance(g, DiGraph):
        boost_graph_from_sage_graph(&g_boost_dir, g, vertex_to_int, reverse)
        i = vertex_to_int[root]
        sig_on()
        result = g_boost_dir.dominator_tree(i)
        sig_off()

    cdef v_index no_parent = -1

    if return_dict:
        return {v: (None if result[i] == no_parent else int_to_vertex[<int> result[i]]) for i, v in enumerate(int_to_vertex)}

    cdef list edges = [[int_to_vertex[<int> result[i]], v] for i, v in enumerate(int_to_vertex) if result[i] != no_parent]

    if g.is_directed():
        if not edges:
            g = DiGraph()
            g.add_vertex(root)
            return g
        else:
            return DiGraph(edges)
    else:
        if not edges:
            g = Graph()
            g.add_vertex(root)
            return g
        else:
            return Graph(edges)


cpdef bandwidth_heuristics(g, algorithm='cuthill_mckee'):
    r"""
    Use Boost heuristics to approximate the bandwidth of the input graph.

    The bandwidth `bw(M)` of a matrix `M` is the smallest integer `k` such that
    all non-zero entries of `M` are at distance `k` from the diagonal. The
    bandwidth `bw(g)` of an undirected graph `g` is the minimum bandwidth of
    the adjacency matrix of `g`, over all possible relabellings of its vertices
    (for more information, see the
    :mod:`~sage.graphs.graph_decompositions.bandwidth`
    module).

    Unfortunately, exactly computing the bandwidth is NP-hard (and an
    exponential algorithm is implemented in Sagemath in routine
    :func:`~sage.graphs.graph_decompositions.bandwidth.bandwidth`). Here, we
    implement two heuristics to find good orderings: Cuthill-McKee, and King.

    This function works only in undirected graphs, and its running time is
    `O(md_{max}\log d_{max})` for the Cuthill-McKee ordering, and
    `O(md_{max}^2\log d_{max})` for the King ordering, where `m` is the number
    of edges, and `d_{max}` is the maximum degree in the graph.

    INPUT:

    - ``g`` -- the input Sage graph

    - ``algorithm`` -- string (default: ``'cuthill_mckee'``); the heuristic used
      to compute the ordering among ``'cuthill_mckee'`` and ``'king'``

    OUTPUT:

    A pair ``[bandwidth, ordering]``, where ``ordering`` is the ordering of
    vertices, ``bandwidth`` is the bandwidth of that specific ordering (which
    is not necessarily the bandwidth of the graph, because this is a heuristic).

    EXAMPLES::

        sage: from sage.graphs.base.boost_graph import bandwidth_heuristics
        sage: bandwidth_heuristics(graphs.PathGraph(10))
        (1, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
        sage: bandwidth_heuristics(graphs.GridGraph([3,3]))  # py2
        (3, [(2, 2), (2, 1), (1, 2), (2, 0), (1, 1), (0, 2), (1, 0), (0, 1), (0, 0)])
        sage: bandwidth_heuristics(graphs.GridGraph([3,3]), algorithm='king')  # py2
        (3, [(2, 2), (2, 1), (1, 2), (2, 0), (1, 1), (0, 2), (1, 0), (0, 1), (0, 0)])
        sage: bandwidth_heuristics(graphs.GridGraph([3,3]))  # py3
        (3, [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2), (2, 1), (1, 2), (2, 2)])
        sage: bandwidth_heuristics(graphs.GridGraph([3,3]), algorithm='king')  # py3
        (3, [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2), (2, 1), (1, 2), (2, 2)])

    TESTS:

    Given an input which is not a graph::

        sage: from sage.graphs.base.boost_graph import bandwidth_heuristics
        sage: bandwidth_heuristics(digraphs.Path(10))
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage Graph
        sage: bandwidth_heuristics("I am not a graph!")
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage Graph

    Given a wrong algorithm::

        from sage.graphs.base.boost_graph import bandwidth_heuristics
        sage: bandwidth_heuristics(graphs.PathGraph(3), algorithm='tip top')
        Traceback (most recent call last):
        ...
        ValueError: unknown algorithm 'tip top'

    Given a graph with no edges::

        from sage.graphs.base.boost_graph import bandwidth_heuristics
        sage: bandwidth_heuristics(Graph())
        (0, [])
        sage: bandwidth_heuristics(graphs.RandomGNM(10,0))
        (0, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])

    """
    from sage.graphs.graph import Graph

    # Tests for errors and trivial cases
    if not isinstance(g, Graph):
        raise TypeError("the input must be a Sage Graph")
    if not algorithm in ['cuthill_mckee', 'king']:
        raise ValueError(f"unknown algorithm {algorithm!r}")
    if not g.num_edges():
        return (0, list(g))

    # These variables are automatically deleted when the function terminates.
    cdef BoostVecGraph g_boost
    cdef vector[v_index] result
    cdef v_index i
    cdef list int_to_vertex = list(g)
    cdef dict vertex_to_int = {v: i for i, v in enumerate(int_to_vertex)}

    boost_graph_from_sage_graph(&g_boost, g, vertex_to_int)
    cdef bint use_cuthill_mckee = (algorithm == 'cuthill_mckee')
    sig_on()
    result = g_boost.bandwidth_ordering(use_cuthill_mckee)
    sig_off()

    cdef int n = g.num_verts()
    cdef dict pos = {int_to_vertex[<int> result[i]]: i for i in range(n)}
    cdef int bandwidth = max([abs(pos[u] - pos[v]) for u, v in g.edge_iterator(labels=False)])

    return (bandwidth, [int_to_vertex[<int> result[i]] for i in range(n)])


cpdef min_spanning_tree(g,
                        weight_function=None,
                        algorithm='Kruskal'):
    r"""
    Use Boost to compute the minimum spanning tree of the input graph.

    INPUT:

    - ``g`` -- the input Sage graph

    - ``weight_function`` -- function (default: ``None``); a function that
      inputs an edge ``e`` and outputs its weight. An edge has the form
      ``(u,v,l)``, where ``u`` and ``v`` are vertices, ``l`` is a label (that
      can be of any kind). The ``weight_function`` can be used to transform the
      label into a weight (see the example below). In particular:

      - if ``weight_function`` is not ``None``, the weight of an edge ``e`` is
        ``weight_function(e)``;

      - if ``weight_function`` is ``None`` (default) and ``g`` is weighted (that
        is, ``g.weighted()==True``), for each edge ``e=(u,v,l)``, we set weight
        ``l``;

      - if ``weight_function`` is ``None`` and ``g`` is not weighted, we set all
        weights to 1 (hence, the output can be any spanning tree).

      Note that, if the weight is not convertible to a number with function
      ``float()``, an error is raised (see tests below).

    - ``algorithm`` -- string (default: ``'Kruskal'``); the algorithm to use
      among ``'Kruskal'`` and ``'Prim'``

    OUTPUT:

    The edges of a minimum spanning tree of ``g``, if one exists, otherwise
    the empty list.

    .. SEEALSO::

        - :meth:`sage.graphs.generic_graph.GenericGraph.min_spanning_tree`

    EXAMPLES::

        sage: from sage.graphs.base.boost_graph import min_spanning_tree
        sage: min_spanning_tree(graphs.PathGraph(4))
        [(0, 1, None), (1, 2, None), (2, 3, None)]

        sage: G = Graph([(0,1,{'name':'a','weight':1}), (0,2,{'name':'b','weight':3}), (1,2,{'name':'b','weight':1})])
        sage: min_spanning_tree(G, weight_function=lambda e: e[2]['weight'])
        [(0, 1, {'name': 'a', 'weight': 1}), (1, 2, {'name': 'b', 'weight': 1})]

    TESTS:

    Given an input which is not a graph::

        sage: min_spanning_tree("I am not a graph!")
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage Graph

    Given a wrong algorithm::

        sage: min_spanning_tree(graphs.PathGraph(3), algorithm='tip top')
        Traceback (most recent call last):
        ...
        ValueError: algorithm 'tip top' not yet implemented, please contribute

    If the weight is not a number::

        sage: g = Graph([(0,1,1), (1,2,'a')], weighted=True)
        sage: min_spanning_tree(g)
        Traceback (most recent call last):
        ...
        ValueError: could not convert string to float:...

        sage: g = Graph([(0,1,1), (1,2,[1,2,3])], weighted=True)
        sage: min_spanning_tree(g)
        Traceback (most recent call last):
        ...
        TypeError: float() argument must be a string or a number...
    """
    from sage.graphs.graph import Graph

    if not isinstance(g, Graph):
        raise TypeError("the input must be a Sage Graph")
    if not algorithm in ['Kruskal', 'Prim']:
        raise ValueError("algorithm '%s' not yet implemented, please contribute" %(algorithm))

    if g.allows_loops() or g.allows_multiple_edges():
        g = g.to_simple()
    # Now g has no self loops and no multiple edges.
    # These variables are automatically deleted when the function terminates.
    cdef BoostVecWeightedGraph g_boost
    cdef vector[v_index] result
    cdef v_index i
    cdef list int_to_vertex = list(g)
    cdef dict vertex_to_int = {v: i for i, v in enumerate(int_to_vertex)}

    boost_weighted_graph_from_sage_graph(&g_boost, g, vertex_to_int, weight_function)

    if algorithm == 'Kruskal':
        sig_on()
        result = g_boost.kruskal_min_spanning_tree()
        sig_off()
    elif algorithm == 'Prim':
        sig_on()
        result = g_boost.prim_min_spanning_tree()
        sig_off()

    cdef size_t n = g.num_verts()

    if result.size() != 2 * (n - 1):
        return []
    else:
        edges = [(int_to_vertex[<int> result[2*i]], int_to_vertex[<int> result[2*i+1]]) for i in range(n-1)]
        return [(min(e[0],e[1]), max(e[0],e[1]), g.edge_label(e[0], e[1])) for e in edges]


cpdef blocks_and_cut_vertices(g):
    r"""
    Compute the blocks and cut vertices of the graph.

    This method uses the implementation of Tarjan's algorithm available in the
    Boost library .

    INPUT:

    - ``g`` -- the input Sage graph

    OUTPUT:

    A 2-dimensional vector with m+1 rows (m is the number of biconnected
    components), where each of the first m rows correspond to vertices in a
    block, and the last row is the list of cut vertices.

    .. SEEALSO::

        - :meth:`sage.graphs.generic_graph.GenericGraph.blocks_and_cut_vertices`

    EXAMPLES::

        sage: from sage.graphs.base.boost_graph import blocks_and_cut_vertices
        sage: g = graphs.KrackhardtKiteGraph()
        sage: blocks_and_cut_vertices(g)
        ([[8, 9], [7, 8], [0, 1, 2, 3, 5, 4, 6, 7]], [8, 7])

        sage: G = Graph([(0,1,{'name':'a','weight':1}), (0,2,{'name':'b','weight':3}), (1,2,{'name':'b','weight':1})])
        sage: blocks_and_cut_vertices(G)
        ([[0, 1, 2]], [])

    TESTS:

    Given an input which is not a graph::

        sage: blocks_and_cut_vertices("I am not a graph!")
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage graph
    """
    from sage.graphs.generic_graph import GenericGraph

    if not isinstance(g, GenericGraph):
        raise TypeError("the input must be a Sage graph")

    if g.allows_loops() or g.allows_multiple_edges() or g.is_directed():
        g = g.to_simple()

    cdef BoostVecGraph g_boost
    cdef vector[vector[v_index]] result
    cdef v_index vi
    cdef list int_to_vertex = list(g)
    cdef dict vertex_to_int = {vv: vi for vi, vv in enumerate(int_to_vertex)}
    cdef list vertex_status = [-1] * g.order()

    boost_graph_from_sage_graph(&g_boost, g, vertex_to_int)
    sig_on()
    result = g_boost.blocks_and_cut_vertices()
    sig_off()

    cdef list result_blocks = []
    cdef set result_cut = set()
    cdef list result_temp = []
    cdef v_index v

    # We iterate over the vertices in the blocks and find articulation points
    for i in range(len(result)):
        for v in result[i]:
            # The vertex is seen for the first time
            if vertex_status[v] == -1:
                result_temp.append(int_to_vertex[<int> v])
                vertex_status[v] = i
            # Vertex belongs to a previous block also, must be a cut vertex
            elif vertex_status[v] < i:
                result_cut.add(int_to_vertex[<int> v])
                result_temp.append(int_to_vertex[<int> v])
                # Change the block number to avoid adding the vertex twice 
                # as a cut vertex if it is repeated in block i
                vertex_status[v] = i
            # elif vertex_status[v] == i:
            # Nothing to do since we have already added the vertex to block i

        result_blocks.append(result_temp)
        result_temp = []

    # If a vertex does not belong to any block, it must be an isolated vertex.
    # Hence, it is considered a block.
    for i in range(g.order()):
        if vertex_status[i] == -1:
            result_blocks.append([int_to_vertex[<int> i]])

    return (result_blocks, list(result_cut))


cpdef shortest_paths(g, start, weight_function=None, algorithm=None):
    r"""
    Compute the shortest paths from ``start`` to all other vertices.

    This routine outputs all shortest paths from node ``start`` to any other
    node in the graph. The input graph can be weighted: if the algorithm is
    Dijkstra, no negative weights are allowed, while if the algorithm is
    Bellman-Ford, negative weights are allowed, but there must be no negative
    cycle (otherwise, the shortest paths might not exist).

    However, Dijkstra algorithm is more efficient: for this reason, we suggest
    to use Bellman-Ford only if necessary (which is also the default option).
    Note that, if the graph is undirected, a negative edge automatically creates
    a negative cycle: for this reason, in this case, Dijkstra algorithm is
    always better.

    The running-time is `O(n \log n+m)` for Dijkstra algorithm and `O(mn)` for
    Bellman-Ford algorithm, where `n` is the number of nodes and `m` is the
    number of edges.

    INPUT:

    - ``g`` -- the input Sage graph

    - ``start`` -- the starting vertex to compute shortest paths

    - ``weight_function`` -- function (default: ``None``); a function that
      associates a weight to each edge. If ``None`` (default), the weights of
      ``g`` are used, if available, otherwise all edges have weight 1.

    - ``algorithm`` -- string (default: ``None``); one of the following
      algorithms:

      - ``'Dijkstra'``, ``'Dijkstra_Boost'``: the Dijkstra algorithm implemented
        in Boost (works only with positive weights)

      - ``'Bellman-Ford'``, ``'Bellman-Ford_Boost'``: the Bellman-Ford algorithm
        implemented in Boost (works also with negative weights, if there is no
        negative cycle)

    OUTPUT:

    A pair of dictionaries ``(distances, predecessors)`` such that, for each
    vertex ``v``, ``distances[v]`` is the distance from ``start`` to ``v``,
    ``predecessors[v]`` is the last vertex in a shortest path from ``start`` to
    ``v``.

    EXAMPLES:

    Undirected graphs::

        sage: from sage.graphs.base.boost_graph import shortest_paths
        sage: g = Graph([(0,1,1),(1,2,2),(1,3,4),(2,3,1)], weighted=True)
        sage: shortest_paths(g, 1)
        ({0: 1, 1: 0, 2: 2, 3: 3}, {0: 1, 1: None, 2: 1, 3: 2})
        sage: g = graphs.GridGraph([2,2])
        sage: shortest_paths(g,(0,0),weight_function=lambda e:2)
        ({(0, 0): 0, (0, 1): 2, (1, 0): 2, (1, 1): 4},
         {(0, 0): None, (0, 1): (0, 0), (1, 0): (0, 0), (1, 1): (0, 1)})

    Directed graphs::

        sage: g = DiGraph([(0,1,1),(1,2,2),(1,3,4),(2,3,1)], weighted=True)
        sage: shortest_paths(g, 1)
        ({1: 0, 2: 2, 3: 3}, {1: None, 2: 1, 3: 2})

    TESTS:

    Given an input which is not a graph::

        sage: shortest_paths("I am not a graph!", 1)
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage graph

    If there is a negative cycle::

        sage: from sage.graphs.base.boost_graph import shortest_paths
        sage: g = DiGraph([(0,1,1),(1,2,-2),(2,0,0.5),(2,3,1)], weighted=True)
        sage: shortest_paths(g, 1)
        Traceback (most recent call last):
        ...
        ValueError: the graph contains a negative cycle

    If Dijkstra is used with negative weights::

        sage: from sage.graphs.base.boost_graph import shortest_paths
        sage: g = DiGraph([(0,1,1),(1,2,-2),(1,3,4),(2,3,1)], weighted=True)
        sage: shortest_paths(g, 1, algorithm='Dijkstra')
        Traceback (most recent call last):
        ...
        RuntimeError: Dijkstra algorithm does not work with negative weights, use Bellman-Ford instead

    Wrong starting vertex::

        sage: shortest_paths(g, 55)
        Traceback (most recent call last):
        ...
        ValueError: the starting vertex 55 is not in the graph
    """
    from sage.graphs.generic_graph import GenericGraph

    if not isinstance(g, GenericGraph):
        raise TypeError("the input must be a Sage graph")

    if start not in g:
        raise ValueError("the starting vertex " + str(start) + " is not in " +
                         "the graph")

    if not g.num_edges():
        return ({start:0}, {start:None})

    # These variables are automatically deleted when the function terminates.
    cdef v_index vi
    cdef dict int_to_v = dict(enumerate(g))
    cdef dict v_to_int = {vv: vi for vi, vv in enumerate(g)}
    cdef BoostVecWeightedDiGraphU g_boost_dir
    cdef BoostVecWeightedGraph g_boost_und
    cdef result_distances result

    if algorithm is None:
        # Check if there are edges with negative weights
        if weight_function is not None:
            for e in g.edge_iterator():
                if float(weight_function(e)) < 0:
                    algorithm = 'Bellman-Ford'
                    break
        else:
            for _,_,w in g.edge_iterator():
                if float(w) < 0:
                    algorithm = 'Bellman-Ford'
                    break

        if algorithm is None:
            algorithm = 'Dijkstra'

    if algorithm in ['Bellman-Ford', 'Bellman-Ford_Boost']:
        if g.is_directed():
            boost_weighted_graph_from_sage_graph(&g_boost_dir, g, v_to_int, weight_function)
            vi = v_to_int[start]
            sig_on()
            result = g_boost_dir.bellman_ford_shortest_paths(vi)
            sig_off()
        else:
            boost_weighted_graph_from_sage_graph(&g_boost_und, g, v_to_int, weight_function)
            vi = v_to_int[start]
            sig_on()
            result = g_boost_und.bellman_ford_shortest_paths(vi)
            sig_off()
        if not result.distances.size():
            raise ValueError("the graph contains a negative cycle")

    elif algorithm in ['Dijkstra', 'Dijkstra_Boost']:
        if g.is_directed():
            boost_weighted_graph_from_sage_graph(&g_boost_dir, g, v_to_int, weight_function)
            vi = v_to_int[start]
            sig_on()
            result = g_boost_dir.dijkstra_shortest_paths(vi)
            sig_off()
        else:
            boost_weighted_graph_from_sage_graph(&g_boost_und, g, v_to_int, weight_function)
            vi = v_to_int[start]
            sig_on()
            result = g_boost_und.dijkstra_shortest_paths(vi)
            sig_off()
        if not result.distances.size():
            raise RuntimeError("Dijkstra algorithm does not work with negative weights, use Bellman-Ford instead")

    else:
        raise ValueError(f"unknown algorithm {algorithm!r}")

    dist = {}
    pred = {}

    if weight_function is not None:
        correct_type = type(weight_function(next(g.edge_iterator())))
    elif g.weighted():
        correct_type = type(next(g.edge_iterator())[2])
    else:
        correct_type = int
    # Needed for rational curves.
    from sage.rings.real_mpfr import RealNumber, RR
    if correct_type == RealNumber:
        correct_type = RR

    import sys
    for v in range(g.num_verts()):
        if result.distances[v] != sys.float_info.max:
            w = int_to_v[v]
            dist[w] = correct_type(result.distances[v])
            pred[w] = int_to_v[result.predecessors[v]] if result.predecessors[v] != v else None
    return (dist, pred)


cdef get_predecessors(BoostWeightedGraph g, result, int_to_v, directed, weight_type):
    r"""
    Return the predecessor matrix from the distance matrix of the graph.
    
    INPUT:
    
    - ``g`` -- the input boost graph
    
    - ``result`` -- the matrix of shortest distances
    
    - ``int_to_v`` -- a list; it is a mapping from `(0, \ldots, n-1)`
      to the vertex set of the original sage graph.
    
    - ``directed`` -- boolean; whether the input graph is directed
    
    - ``weight_type`` -- correct data type for edge weights
    
    OUTPUT:
    
    A dictionary of dictionaries ``pred`` such that ``pred[u][v]`` indicates 
    the predecessor  of `v` in the shortest path from `u` to `v`.

    TESTS::

        sage: from sage.graphs.base.boost_graph import johnson_shortest_paths
        sage: g = Graph([(0,1,1),(1,2,2),(1,3,4),(2,3,1)], weighted=True)
        sage: expected = {0: {0: None, 1: 0, 2: 1, 3: 2},
        ....:             1: {0: 1, 1: None, 2: 1, 3: 2},
        ....:             2: {0: 1, 1: 2, 2: None, 3: 2},
        ....:             3: {0: 1, 1: 2, 2: 3, 3: None}}
        sage: johnson_shortest_paths(g, distances=False, predecessors=True) == expected
        True
    """
    cdef vector[pair[int, pair[int, double]]] edges
    sig_on()
    edges = g.edge_list()
    sig_off()
    cdef int N = g.num_verts()
    cdef dict pred = {v: {v: None} for v in int_to_v}
    import sys
    for p in edges:
        dst = weight_type(p.second.second)
        # dst is the weight of the edge (u, v)
        u = p.first
        v = p.second.first
        for k in range(N):
            if result[k][u] == sys.float_info.max or result[k][v] == sys.float_info.max:
                continue
            if weight_type(result[k][u]) + dst == weight_type(result[k][v]):
                pred[int_to_v[k]][int_to_v[v]] = int_to_v[u]
            if directed:
                continue
            if weight_type(result[k][u]) == weight_type(result[k][v]) + dst:
                pred[int_to_v[k]][int_to_v[u]] = int_to_v[v]
    return pred

cpdef johnson_shortest_paths(g, weight_function=None, distances=True, predecessors=False):
    r"""
    Use Johnson algorithm to solve the all-pairs-shortest-paths.

    This routine outputs the distance between each pair of vertices and the 
    predecessors matrix (depending on the values of boolean ``distances`` and
    ``predecessors``) using a dictionary of dictionaries. It works on all kinds
    of graphs, but it is designed specifically for graphs with negative weights
    (otherwise there are more efficient algorithms, like Dijkstra).

    The time-complexity is `O(mn\log n)`, where `n` is the number of nodes and
    `m` is the number of edges.

    INPUT:

    - ``g`` -- the input Sage graph

    - ``weight_function`` -- function (default: ``None``); a function that
      associates a weight to each edge. If ``None`` (default), the weights of
      ``g`` are used, if available, otherwise all edges have weight 1.
      
    - ``distances`` -- boolean (default: ``True``); whether to return the
      dictionary of shortest distances
      
    - ``predecessors`` -- boolean (default: ``False``); whether to return the
      predecessors matrix 

    OUTPUT:

    Depending on the input, this function return the dictionary of
    predecessors, the dictionary of distances, or a pair of dictionaries
    ``(distances, predecessors)`` where ``distance[u][v]`` denotes the distance
    of a shortest path from `u` to `v` and ``predecessors[u][v]`` indicates the
    predecessor of `w` on a shortest path from `u` to `v`.

    EXAMPLES:

    Undirected graphs::

        sage: from sage.graphs.base.boost_graph import johnson_shortest_paths
        sage: g = Graph([(0,1,1),(1,2,2),(1,3,4),(2,3,1)], weighted=True)
        sage: johnson_shortest_paths(g)
        {0: {0: 0, 1: 1, 2: 3, 3: 4},
         1: {0: 1, 1: 0, 2: 2, 3: 3},
         2: {0: 3, 1: 2, 2: 0, 3: 1},
         3: {0: 4, 1: 3, 2: 1, 3: 0}}
        sage: expected = {0: {0: None, 1: 0, 2: 1, 3: 2},
        ....:             1: {0: 1, 1: None, 2: 1, 3: 2},
        ....:             2: {0: 1, 1: 2, 2: None, 3: 2},
        ....:             3: {0: 1, 1: 2, 2: 3, 3: None}}
        sage: johnson_shortest_paths(g, distances=False, predecessors=True) == expected
        True

    Directed graphs::

        sage: g = DiGraph([(0,1,1),(1,2,-2),(1,3,4),(2,3,1)], weighted=True)
        sage: johnson_shortest_paths(g)
        {0: {0: 0, 1: 1, 2: -1, 3: 0},
         1: {1: 0, 2: -2, 3: -1},
         2: {2: 0, 3: 1},
         3: {3: 0}}
        sage: g = DiGraph([(1,2,3),(2,3,2),(1,4,1),(4,2,1)], weighted=True)
        sage: johnson_shortest_paths(g, distances=False, predecessors=True)
        {1: {1: None, 2: 4, 3: 2, 4: 1},
         2: {2: None, 3: 2},
         3: {3: None},
         4: {2: 4, 3: 2, 4: None}}

    TESTS:

    Given an input which is not a graph::

        sage: johnson_shortest_paths("I am not a graph!")
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage graph

    If there is a negative cycle::

        sage: g = DiGraph([(0,1,1),(1,2,-2),(2,0,0.5),(2,3,1)], weighted=True)
        sage: johnson_shortest_paths(g)
        Traceback (most recent call last):
        ...
        ValueError: the graph contains a negative cycle
    """
    from sage.graphs.generic_graph import GenericGraph
    cdef dict dist = {}
    cdef dict pred = {}

    if not isinstance(g, GenericGraph):
        raise TypeError("the input must be a Sage graph")
    elif not g.num_edges():
        dist = {v: {v: 0} for v in g}
        pred = {v: {v: None} for v in g}
        if distances and predecessors:
            return (dist, pred)
        if distances:
            return dist
        if predecessors:
            return pred
    # These variables are automatically deleted when the function terminates.
    cdef v_index i
    cdef list int_to_v = list(g)
    cdef dict v_to_int = {v: i for i, v in enumerate(int_to_v)}
    cdef BoostVecWeightedDiGraphU g_boost_dir
    cdef BoostVecWeightedGraph g_boost_und
    cdef int N = g.num_verts()
    cdef vector[vector[double]] result
    cdef int u_int, v_int

    if g.is_directed():
        boost_weighted_graph_from_sage_graph(&g_boost_dir, g, v_to_int, weight_function)
        sig_on()
        result = g_boost_dir.johnson_shortest_paths()
        sig_off()
    else:
        boost_weighted_graph_from_sage_graph(&g_boost_und, g, v_to_int, weight_function)
        sig_on()
        result = g_boost_und.johnson_shortest_paths()
        sig_off()

    if not result.size():
        raise ValueError("the graph contains a negative cycle")

    if weight_function is not None:
        correct_type = type(weight_function(next(g.edge_iterator())))
    elif g.weighted():
        correct_type = type(next(g.edge_iterator())[2])
    else:
        correct_type = int
    # Needed for rational curves.
    from sage.rings.real_mpfr import RealNumber, RR
    if correct_type == RealNumber:
        correct_type = RR

    import sys
    if distances:
        dist = {int_to_v[v]: {int_to_v[w]: correct_type(result[v][w])
                              for w in range(N) if result[v][w] != sys.float_info.max}
                for v in range(N)}

    if predecessors:
        if g.is_directed():
            pred = get_predecessors(g_boost_dir, result, int_to_v, True, correct_type)
        else:
            pred = get_predecessors(g_boost_und, result, int_to_v, False, correct_type)

    if distances and predecessors:
        return (dist, pred)
    if distances:
        return dist
    if predecessors:
        return pred

cpdef floyd_warshall_shortest_paths(g, weight_function=None, distances=True, predecessors=False):
    r"""
    Use Floyd-Warshall algorithm to solve the all-pairs-shortest-paths.

    This routine outputs the distance between each pair of vertices and the 
    predecessors matrix (depending on the values of boolean ``distances`` and
    ``predecessors``) using a dictionary of dictionaries. This method should be
    preferred only if the graph is dense. If the graph is sparse the much
    faster johnson_shortest_paths should be used.
    
    The time-complexity is `O(n^3 + nm)`, where `n` is the number of nodes and
    `m` the number of edges. The factor `nm` in the complexity is added only
    when ``predecessors`` is set to ``True``.
    
    INPUT:

    - ``g`` -- the input Sage graph

    - ``weight_function`` -- function (default: ``None``); a function that
      associates a weight to each edge. If ``None`` (default), the weights of
      ``g`` are used, if available, otherwise all edges have weight 1.
      
    - ``distances`` -- boolean (default: ``True``); whether to return
      the dictionary of shortest distances
      
    - ``predecessors`` -- boolean (default: ``False``); whether to return the
      predecessors matrix 

    OUTPUT:

    Depending on the input, this function return the dictionary of
    predecessors, the dictionary of distances, or a pair of dictionaries
    ``(distances, predecessors)`` where ``distance[u][v]`` denotes the distance
    of a shortest path from `u` to `v` and ``predecessors[u][v]`` indicates the
    predecessor of `w` on a shortest path from `u` to `v`.

    EXAMPLES:

    Undirected graphs::

        sage: from sage.graphs.base.boost_graph import floyd_warshall_shortest_paths
        sage: g = Graph([(0,1,1),(1,2,2),(1,3,4),(2,3,1)], weighted=True)
        sage: floyd_warshall_shortest_paths(g)
        {0: {0: 0, 1: 1, 2: 3, 3: 4},
         1: {0: 1, 1: 0, 2: 2, 3: 3},
         2: {0: 3, 1: 2, 2: 0, 3: 1},
         3: {0: 4, 1: 3, 2: 1, 3: 0}}
        sage: expected = {0: {0: None, 1: 0, 2: 1, 3: 2},
        ....:             1: {0: 1, 1: None, 2: 1, 3: 2},
        ....:             2: {0: 1, 1: 2, 2: None, 3: 2},
        ....:             3: {0: 1, 1: 2, 2: 3, 3: None}}
        sage: floyd_warshall_shortest_paths(g, distances=False, predecessors=True) == expected
        True

    Directed graphs::

        sage: g = DiGraph([(0,1,1),(1,2,-2),(1,3,4),(2,3,1)], weighted=True)
        sage: floyd_warshall_shortest_paths(g)
        {0: {0: 0, 1: 1, 2: -1, 3: 0},
         1: {1: 0, 2: -2, 3: -1},
         2: {2: 0, 3: 1},
         3: {3: 0}}
        sage: g = DiGraph([(1,2,3),(2,3,2),(1,4,1),(4,2,1)], weighted=True)
        sage: floyd_warshall_shortest_paths(g, distances=False, predecessors=True)
        {1: {1: None, 2: 4, 3: 2, 4: 1},
         2: {2: None, 3: 2},
         3: {3: None},
         4: {2: 4, 3: 2, 4: None}}

    TESTS:

    Given an input which is not a graph::

        sage: floyd_warshall_shortest_paths("I am not a graph!")
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage graph

    If there is a negative cycle:

        sage: g = DiGraph([(0,1,1),(1,2,-2),(2,0,0.5),(2,3,1)], weighted=True)
        sage: floyd_warshall_shortest_paths(g)
        Traceback (most recent call last):
        ...
        ValueError: the graph contains a negative cycle
    """
    from sage.graphs.generic_graph import GenericGraph
    cdef dict dist = {}
    cdef dict pred = {}

    if not isinstance(g, GenericGraph):
        raise TypeError("the input must be a Sage graph")
    elif not g.num_edges():
        dist = {v: {v: 0} for v in g}
        pred = {v: {v: None} for v in g}
        if distances and predecessors:
            return (dist, pred)
        if distances:
            return dist
        if predecessors:
            return pred
    # These variables are automatically deleted when the function terminates.
    cdef v_index i
    cdef list int_to_v = list(g)
    cdef dict v_to_int = {v: i for i, v in enumerate(int_to_v)}
    cdef BoostVecWeightedDiGraphU g_boost_dir
    cdef BoostVecWeightedGraph g_boost_und
    cdef int N = g.num_verts()
    cdef vector[vector[double]] result
    cdef int u_int, v_int

    if g.is_directed():
        boost_weighted_graph_from_sage_graph(&g_boost_dir, g, v_to_int, weight_function)
        sig_on()
        result = g_boost_dir.floyd_warshall_shortest_paths()
        sig_off()
    else:
        boost_weighted_graph_from_sage_graph(&g_boost_und, g, v_to_int, weight_function)
        sig_on()
        result = g_boost_und.floyd_warshall_shortest_paths()
        sig_off()

    if not result.size():
        raise ValueError("the graph contains a negative cycle")

    if weight_function is not None:
        correct_type = type(weight_function(next(g.edge_iterator())))
    elif g.weighted():
        correct_type = type(next(g.edge_iterator())[2])
    else:
        correct_type = int
    # Needed for rational curves.
    from sage.rings.real_mpfr import RealNumber, RR
    if correct_type == RealNumber:
        correct_type = RR

    import sys
    if distances:
        dist = {int_to_v[v]: {int_to_v[w]: correct_type(result[v][w])
                              for w in range(N) if result[v][w] != sys.float_info.max}
                for v in range(N)}

    if predecessors:
        if g.is_directed():
            pred = get_predecessors(g_boost_dir, result, int_to_v, True, correct_type)
        else:
            pred = get_predecessors(g_boost_und, result, int_to_v, False, correct_type)

    if distances and predecessors:
        return (dist, pred)
    if distances:
        return dist
    if predecessors:
        return pred


cpdef johnson_closeness_centrality(g, weight_function=None):
    r"""
    Use Johnson algorithm to compute the closeness centrality of all vertices.

    This routine is preferrable to :func:`~johnson_shortest_paths` because it
    does not create a doubly indexed dictionary of distances, saving memory.

    The time-complexity is `O(mn\log n)`, where `n` is the number of nodes and
    `m` is the number of edges.

    INPUT:

    - ``g`` -- the input Sage graph

    - ``weight_function`` -- function (default: ``None``); a function that
      associates a weight to each edge. If ``None`` (default), the weights of
      ``g`` are used, if available, otherwise all edges have weight 1.

    OUTPUT:

    A dictionary associating each vertex ``v`` to its closeness centrality.

    EXAMPLES:

    Undirected graphs::

        sage: from sage.graphs.base.boost_graph import johnson_closeness_centrality
        sage: g = Graph([(0,1,1),(1,2,2),(1,3,4),(2,3,1)], weighted=True)
        sage: johnson_closeness_centrality(g)
        {0: 0.375, 1: 0.5, 2: 0.5, 3: 0.375}

    Directed graphs::

        sage: from sage.graphs.base.boost_graph import johnson_closeness_centrality
        sage: g = DiGraph([(0,1,1),(1,2,-2),(1,3,4),(2,3,1)], weighted=True)
        sage: johnson_closeness_centrality(g)
        {0: inf, 1: -0.4444444444444444, 2: 0.3333333333333333}

    TESTS:

    Given an input which is not a graph::

        sage: from sage.graphs.base.boost_graph import johnson_closeness_centrality
        sage: johnson_closeness_centrality("I am not a graph!")
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage graph

    If there is a negative cycle::

        sage: from sage.graphs.base.boost_graph import johnson_closeness_centrality
        sage: g = DiGraph([(0,1,1),(1,2,-2),(2,0,0.5),(2,3,1)], weighted=True)
        sage: johnson_closeness_centrality(g)
        Traceback (most recent call last):
        ...
        ValueError: the graph contains a negative cycle
    """
    from sage.graphs.generic_graph import GenericGraph

    if not isinstance(g, GenericGraph):
        raise TypeError("the input must be a Sage graph")
    elif not g.num_edges():
        return {}
    # These variables are automatically deleted when the function terminates.
    cdef BoostVecWeightedDiGraphU g_boost_dir
    cdef BoostVecWeightedGraph g_boost_und
    cdef int N = g.num_verts()
    cdef vector[vector[double]] result
    cdef vector[double] closeness
    cdef double farness
    cdef int i, j, reach
    cdef list int_to_v = list(g)
    cdef dict v_to_int = {v: i for i, v in enumerate(int_to_v)}

    if g.is_directed():
        boost_weighted_graph_from_sage_graph(&g_boost_dir, g, v_to_int, weight_function)
        sig_on()
        result = g_boost_dir.johnson_shortest_paths()
        sig_off()
    else:
        boost_weighted_graph_from_sage_graph(&g_boost_und, g, v_to_int, weight_function)
        sig_on()
        result = g_boost_und.johnson_shortest_paths()
        sig_off()

    if not result.size():
        raise ValueError("the graph contains a negative cycle")

    import sys
    for i in range(N):
        farness = 0
        reach = 0
        for j in range(N):
            if result[i][j] != sys.float_info.max:
                farness += result[i][j]
                reach += 1
        if reach > 1:
            closeness.push_back((<double>reach-1) * (reach-1) / ((N-1) * farness))
        else:
            closeness.push_back(sys.float_info.max)
        sig_check()
    return {v: closeness[i] for i,v in enumerate(int_to_v) if closeness[i] != sys.float_info.max}

cpdef min_cycle_basis(g_sage, weight_function=None, by_weight=False):
    r"""
    Return a minimum weight cycle basis of the input graph ``g_sage``.

    A cycle basis is a list of cycles (list of vertices forming a cycle) of
    ``g_sage``. Note that the vertices are not necessarily returned in the order
    in which they appear in the cycle.

    A minimum weight cycle basis is a cycle basis that minimizes the sum of the
    weights (length for unweighted graphs) of its cycles.

    Not implemented for directed graphs and multigraphs.

    INPUT:

    - ``g_sage`` -- a Sage Graph

    - ``weight_function`` -- function (default: ``None``); a function that takes
      as input an edge ``(u, v, l)`` and outputs its weight. If not ``None``,
      ``by_weight`` is automatically set to ``True``. If ``None`` and
      ``by_weight`` is ``True``, we use the edge label ``l`` as a weight.

    - ``by_weight`` -- boolean (default: ``False``); if ``True``, the edges in
      the graph are weighted, otherwise all edges have weight 1

    EXAMPLES::

        sage: g = Graph([(1, 2, 3), (2, 3, 5), (3, 4, 8), (4, 1, 13), (1, 3, 250), (5, 6, 9), (6, 7, 17), (7, 5, 20)])
        sage: sorted(g.minimum_cycle_basis(by_weight=True))
        [[1, 2, 3], [1, 2, 3, 4], [5, 6, 7]]
        sage: sorted(g.minimum_cycle_basis())
        [[1, 2, 3], [1, 3, 4], [5, 6, 7]]

    .. SEEALSO::

        * :wikipedia:`Cycle_basis`
    """
    cdef Py_ssize_t u_int, v_int, i, j
    cdef object u, v
    cdef Py_ssize_t n = g_sage.num_verts()
    cdef list int_to_vertex = list(g_sage)
    cdef dict vertex_to_int = {u: u_int for u_int, u in enumerate(int_to_vertex)}
    cdef list edgelist
    if by_weight and weight_function is not None:
        edgelist = [(vertex_to_int[e[0]], vertex_to_int[e[1]], float(weight_function(e)))
                        for e in g_sage.edge_iterator()]
    if by_weight:
        edgelist = [(vertex_to_int[e[0]], vertex_to_int[e[1]], float(e[2]))
                        for e in g_sage.edge_iterator()]
    else:
        edgelist = [(vertex_to_int[u], vertex_to_int[v], 1) for u, v in g_sage.edge_iterator(labels=False)]

    # We just need the edges of any spanning tree here not necessarily a
    # minimum spanning tree.
    
    cdef list sp_edges = min_spanning_tree(g_sage)
    cdef cset[pair[int, int]] edges_s
    for a, b, c in sp_edges:
        edges_s.insert((a, b))
    # Edges of self that are not in the spanning tree
    cdef list edges_c = [e for e in g_sage.edge_iterator(labels=False) if not edges_s.count(e)]
    cdef list edges_complement = [frozenset((vertex_to_int[u], vertex_to_int[v])) for u, v in edges_c]
    cdef Py_ssize_t l = len(edges_complement)
    cdef list orth_set = [set([e]) for e in edges_complement]
    cdef list cycle_basis = []
    cdef set base
    cdef BoostVecWeightedGraph g_boost_und
    cdef vector[vector[double]] all_pair_shortest_pathlens
    cdef result_distances min_path
    cdef list min_path_nodes
    cdef dict cross_paths_lens
    cdef float edge_w

    for i in range(l):
        base = orth_set[i]
        # For each edge in g, add 2 edges to g_boost_und: "cross" edges if edge
        # is in base, otherwise "in-plane" edges
        g_boost_und = BoostVecWeightedGraph()
        for u_int in range(2 * n):
            g_boost_und.add_vertex()
        for u_int, v_int, edge_w in edgelist:
            # mapping the nodes in g from 0 to n-1
            if frozenset((u_int, v_int)) in base:
                g_boost_und.add_edge(u_int, n + v_int, edge_w)
                g_boost_und.add_edge(n + u_int, v_int, edge_w)
            else:
                g_boost_und.add_edge(u_int, v_int, edge_w)
                g_boost_und.add_edge(n + u_int, n + v_int, edge_w)

        sig_on()
        all_pair_shortest_pathlens = g_boost_und.johnson_shortest_paths()
        sig_off()
        cross_paths_lens = {j: all_pair_shortest_pathlens[j][n + j] for j in range(n)}
        u_int = min(cross_paths_lens, key=cross_paths_lens.get)
        v_int = n + u_int
        sig_on()
        min_path = g_boost_und.dijkstra_shortest_paths(u_int)
        sig_off()
        # Mapping the nodes in G to nodes in g
        min_path_nodes = [v_int if v_int < n else v_int - n]
        while v_int != u_int:
            v_int = min_path.predecessors[v_int]
            min_path_nodes.append(v_int if v_int < n else v_int - n)

        # removal of edges occuring even number of times
        edges = set()
        for edge in zip(min_path_nodes[:-1], min_path_nodes[1:]):
            edges ^= {edge}
        new_cycle = {frozenset(e) for e in edges}
        cycle_basis.append([int_to_vertex[u_int] for u_int in set().union(*new_cycle)])
        # updating orth_set so that i+1, i+2, ...th elements are orthogonal
        # to the newly found cycle
        for j in range(i + 1, l):
            if len(orth_set[j] & new_cycle) % 2:
                orth_set[j] = orth_set[j] ^ base
    return cycle_basis
