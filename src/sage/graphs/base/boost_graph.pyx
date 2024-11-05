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

# ****************************************************************************
#       Copyright (C) 2015 Michele Borassi michele.borassi@imtlucca.it
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

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

    - ``vertex_to_int`` -- dictionary; it is a mapping from the vertex set of
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
        for u, v in g_sage.edge_iterator(labels=None):
            g.add_edge(vertex_to_int[v], vertex_to_int[u])
    else:
        for u, v in g_sage.edge_iterator(labels=None):
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

    - ``vertex_to_int`` -- dictionary; it is a mapping from the vertex set of
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
            for u, v, w in g_sage.edge_iterator():
                g.add_edge(vertex_to_int[v], vertex_to_int[u], float(w))
        else:
            for u, v, w in g_sage.edge_iterator():
                g.add_edge(vertex_to_int[u], vertex_to_int[v], float(w))
    else:
        if reverse:
            for u, v in g_sage.edge_iterator(labels=False):
                g.add_edge(vertex_to_int[v], vertex_to_int[u], 1)
        else:
            for u, v in g_sage.edge_iterator(labels=False):
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
        sage: edge_connectivity(g)
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
        stopgap("The edge connectivity of directed graphs is not implemented "
                "in Boost. The result may be mathematically unreliable.", 18753)

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
        sage: clustering_coeff(g, vertices='abde')
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
    all nonzero entries of `M` are at distance `k` from the diagonal. The
    bandwidth `bw(g)` of an undirected graph `g` is the minimum bandwidth of
    the adjacency matrix of `g`, over all possible relabellings of its vertices
    (for more information, see the
    :mod:`~sage.graphs.graph_decompositions.bandwidth`
    module).

    Unfortunately, exactly computing the bandwidth is NP-hard (and an
    exponential algorithm is implemented in Sagemath in routine
    :func:`~sage.graphs.graph_decompositions.bandwidth.bandwidth`). Here, we
    implement two heuristics to find good orderings: Cuthill-McKee, reverse
    Cuthill-McKee (also known as ``RCM``) and King.

    This function works only in undirected graphs, and its running time is
    `O(md_{max}\log d_{max})` for the Cuthill-McKee ordering, and
    `O(md_{max}^2\log d_{max})` for the King ordering, where `m` is the number
    of edges, and `d_{max}` is the maximum degree in the graph.

    INPUT:

    - ``g`` -- the input Sage graph

    - ``algorithm`` -- string (default: ``'cuthill_mckee'``); the heuristic used
      to compute the ordering among ``'cuthill_mckee'``,
      ``'reverse_cuthill_mckee'`` and ``'king'``

    OUTPUT:

    A pair ``[bandwidth, ordering]``, where ``ordering`` is the ordering of
    vertices, ``bandwidth`` is the bandwidth of that specific ordering (which
    is not necessarily the bandwidth of the graph, because this is a heuristic).

    EXAMPLES::

        sage: from sage.graphs.base.boost_graph import bandwidth_heuristics
        sage: bandwidth_heuristics(graphs.PathGraph(10))
        (1, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
        sage: bandwidth_heuristics(graphs.GridGraph([3,3]))
        (3, [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2), (2, 1), (1, 2), (2, 2)])
        sage: bandwidth_heuristics(graphs.GridGraph([3,3]), algorithm='reverse_cuthill_mckee')
        (3, [(2, 2), (1, 2), (2, 1), (0, 2), (1, 1), (2, 0), (0, 1), (1, 0), (0, 0)])
        sage: bandwidth_heuristics(graphs.GridGraph([3,3]), algorithm='king')
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

        sage: from sage.graphs.base.boost_graph import bandwidth_heuristics
        sage: bandwidth_heuristics(graphs.PathGraph(3), algorithm='tip top')
        Traceback (most recent call last):
        ...
        ValueError: unknown algorithm 'tip top'

    Given a graph with no edges::

        sage: from sage.graphs.base.boost_graph import bandwidth_heuristics
        sage: bandwidth_heuristics(Graph())
        (0, [])
        sage: bandwidth_heuristics(graphs.RandomGNM(10,0))                              # needs networkx
        (0, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    """
    from sage.graphs.graph import Graph

    # Tests for errors and trivial cases
    if not isinstance(g, Graph):
        raise TypeError("the input must be a Sage Graph")
    if algorithm not in ['cuthill_mckee', 'reverse_cuthill_mckee', 'king']:
        raise ValueError(f"unknown algorithm {algorithm!r}")
    if not g.num_edges():
        return (0, list(g))

    cdef bint reverse = False
    if algorithm == 'reverse_cuthill_mckee':
        reverse = True
        algorithm = 'cuthill_mckee'

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
    cdef int bandwidth = max([abs(pos[u] - pos[v])
                              for u, v in g.edge_iterator(labels=False)])

    if reverse:
        return (bandwidth, [int_to_vertex[<int> result[i]]
                            for i in range(n - 1, -1, -1)])
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
        TypeError: float() argument must be a string or a... number...

    Check that the method is robust to incomparable vertices::

        sage: G = Graph([(1, 2, 10), (1, 'a', 1), ('a', 'b', 1), ('b', 2, 1)], weighted=True)
        sage: E = min_spanning_tree(G, algorithm='Kruskal')
        sage: sum(w for _, _, w in E)
        3
        sage: F = min_spanning_tree(G, algorithm='Prim')
        sage: sum(w for _, _, w in F)
        3
    """
    from sage.graphs.graph import Graph

    if not isinstance(g, Graph):
        raise TypeError("the input must be a Sage Graph")
    if algorithm not in ['Kruskal', 'Prim']:
        raise ValueError("algorithm '%s' not yet implemented, please contribute" % (algorithm))

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

    cdef v_index n = g.num_verts()

    if <v_index> result.size() != 2 * (n - 1):
        return []
    edges = [(int_to_vertex[<int> result[2*i]], int_to_vertex[<int> result[2*i + 1]]) for i in range(n - 1)]
    return [(u, v, g.edge_label(u, v)) for u, v in edges]


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
      ``g`` are used, if ``g.weighted()==True``, otherwise all edges have
      weight 1.

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
        return ({start: 0}, {start: None})

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
        elif g.weighted():
            for _, _, w in g.edge_iterator():
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
        try:
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
        except RuntimeError:
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

    - ``int_to_v`` -- list; it is a mapping from `(0, \ldots, n-1)`
      to the vertex set of the original sage graph

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
      ``g`` are used, if ``g.weighted()==True``, otherwise all edges have
      weight 1.

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
      ``g`` are used, if ``g.weighted()==True``, otherwise all edges have
      weight 1.

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

    This routine is preferable to :func:`~johnson_shortest_paths` because it
    does not create a doubly indexed dictionary of distances, saving memory.

    The time-complexity is `O(mn\log n)`, where `n` is the number of nodes and
    `m` is the number of edges.

    INPUT:

    - ``g`` -- the input Sage graph

    - ``weight_function`` -- function (default: ``None``); a function that
      associates a weight to each edge. If ``None`` (default), the weights of
      ``g`` are used, if ``g.weighted()==True``, otherwise all edges have
      weight 1.

    OUTPUT: a dictionary associating each vertex ``v`` to its closeness centrality

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
            closeness.push_back((<double>reach - 1) * (reach - 1) / ((N - 1) * farness))
        else:
            closeness.push_back(sys.float_info.max)
        sig_check()
    return {v: closeness[i] for i, v in enumerate(int_to_v) if closeness[i] != sys.float_info.max}


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
      ``by_weight`` is ``True``, the weights of ``g_sage`` are used, if
      ``g_sage.weighted()==True``, otherwise all edges have weight 1.

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
    cdef Py_ssize_t lec = len(edges_complement)
    cdef list orth_set = [set([e]) for e in edges_complement]
    cdef list cycle_basis = []
    cdef set base
    cdef BoostVecWeightedGraph g_boost_und
    cdef vector[vector[double]] all_pair_shortest_pathlens
    cdef result_distances min_path
    cdef list min_path_nodes
    cdef dict cross_paths_lens
    cdef float edge_w

    for i in range(lec):
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

        # removal of edges occurring even number of times
        edges = set()
        for edge in zip(min_path_nodes[:-1], min_path_nodes[1:]):
            edges ^= {edge}
        new_cycle = {frozenset(e) for e in edges}
        cycle_basis.append([int_to_vertex[u_int] for u_int in set().union(*new_cycle)])
        # updating orth_set so that i+1, i+2, ...th elements are orthogonal
        # to the newly found cycle
        for j in range(i + 1, lec):
            if len(orth_set[j] & new_cycle) % 2:
                orth_set[j] = orth_set[j] ^ base
    return cycle_basis


cpdef eccentricity_DHV(g, vertex_list=None, weight_function=None, check_weight=True):
    r"""
    Return the vector of eccentricities using the algorithm of [Dragan2018]_.

    The array returned is of length `n`, and by default its `i`-th component is
    the eccentricity of the `i`-th vertex in ``g.vertices(sort=True)``,
    if ``vertex_list is None``, otherwise ``ecc[i]`` is the eccentricity of
    vertex ``vertex_list[i]``.

    The algorithm proposed in [Dragan2018]_ is based on the observation that for
    all nodes `v,w\in V`, we have `\max(ecc[v]-d(v,w), d(v,w))\leq ecc[w] \leq
    ecc[v] + d(v,w)`. Also the algorithm iteratively improves upper and lower
    bounds on the eccentricity of each vertex until no further improvements can
    be done.

    INPUT:

    - ``g`` -- the input Sage graph

    - ``vertex_list`` -- list (default: ``None``); a list of `n` vertices
      specifying a mapping from `(0, \ldots, n-1)` to vertex labels in `g`. When
      set, ``ecc[i]`` is the eccentricity of vertex ``vertex_list[i]``. When
      ``vertex_list`` is ``None``, ``ecc[i]`` is the eccentricity of vertex
      ``g.vertices(sort=True)[i]``.

    - ``weight_function`` -- function (default: ``None``); a function that
      associates a weight to each edge. If ``None`` (default), the weights of
      ``g`` are used, if ``g.weighted()==True``, otherwise all edges have
      weight 1.

    - ``check_weight`` -- boolean (default: ``True``); if ``True``, we check
      that the ``weight_function`` outputs a number for each edge

    EXAMPLES::

        sage: from sage.graphs.base.boost_graph import eccentricity_DHV
        sage: G = graphs.BullGraph()
        sage: eccentricity_DHV(G)
        [2.0, 2.0, 2.0, 3.0, 3.0]

    TESTS:

        sage: G = Graph(2)
        sage: eccentricity_DHV(G)
        [+Infinity, +Infinity]
        sage: G = graphs.RandomGNP(20, 0.7)
        sage: eccentricity_DHV(G) == G.eccentricity()
        True
        sage: G = Graph([(0,1,-1)], weighted=True)
        sage: eccentricity_DHV(G)
        Traceback (most recent call last):
        ...
        ValueError: graph contains negative edge weights, use Johnson_Boost instead
    """
    if g.is_directed():
        raise TypeError("the 'DHV' algorithm only works on undirected graphs")

    cdef v_index n = g.order()
    if not n:
        return []
    if n == 1:
        return [0]

    if weight_function and check_weight:
        g._check_weight_function(weight_function)

    if weight_function is not None:
        for e in g.edge_iterator():
            if float(weight_function(e)) < 0:
                raise ValueError("graph contains negative edge weights, use Johnson_Boost instead")
    elif g.weighted():
        for _, _, w in g.edge_iterator():
            if w and float(w) < 0:
                raise ValueError("graph contains negative edge weights, use Johnson_Boost instead")

    if vertex_list is None:
        vertex_list = g.vertices(sort=True)
    elif not len(vertex_list) == n or not set(vertex_list) == set(g):
        raise ValueError("parameter vertex_list is incorrect for this graph")

    # These variables are automatically deleted when the function terminates.
    cdef dict v_to_int = {vv: vi for vi, vv in enumerate(vertex_list)}
    cdef BoostVecWeightedGraph g_boost
    boost_weighted_graph_from_sage_graph(&g_boost, g, v_to_int, weight_function)

    import sys
    cdef v_index u, antipode, v
    cdef double ecc_u, ecc_antipode, tmp
    cdef v_index i, idx

    cdef list active = list(range(n))
    cdef vector[double] ecc_lower_bound
    cdef vector[double] ecc_upper_bound
    cdef vector[double] distances

    ecc_lower_bound.assign(n, 0)
    ecc_upper_bound.assign(n, sys.float_info.max)

    # Algorithm
    while active:
        # Select vertex with minimum eccentricity in active and update
        # eccentricity upper bounds.
        # For this, we select u with minimum eccentricity lower bound in active
        # if ecc_u == ecc_lb[u], we are done. Otherwise, we update eccentricity
        # lower bounds and repeat

        tmp = sys.float_info.max
        for i, v in enumerate(active):
            if ecc_lower_bound[v] < tmp:
                tmp = ecc_lower_bound[v]
                idx = i
        active[idx], active[-1] = active[-1], active[idx]
        u = active.pop()

        # compute distances from u
        sig_on()
        distances = g_boost.dijkstra_shortest_paths(u).distances
        sig_off()

        # Compute eccentricity of u
        ecc_u = 0
        for v in range(n):
            if ecc_u < distances[v]:
                ecc_u = distances[v]
                antipode = v
        ecc_upper_bound[u] = ecc_u

        if ecc_u == sys.float_info.max:  # Disconnected graph
            break

        if ecc_u == ecc_lower_bound[u]:
            # We found the good vertex.
            # Update eccentricity upper bounds and remove from active those
            # vertices for which gap is closed
            i = 0
            while i < len(active):
                v = active[i]
                ecc_upper_bound[v] = min(ecc_upper_bound[v], distances[v] + ecc_u)
                if ecc_upper_bound[v] == ecc_lower_bound[v]:
                    active[i] = active[-1]
                    active.pop()
                else:
                    i += 1

        else:
            # u was not a good choice.
            # We use its antipode to update eccentricity lower bounds.
            # Observe that this antipode might have already been seen.
            for i, v in enumerate(active):
                if v == antipode:
                    active[i] = active[-1]
                    active.pop()
                    break

            # Compute distances from antipode
            sig_on()
            distances = g_boost.dijkstra_shortest_paths(antipode).distances
            sig_off()

            # Compute eccentricity of antipode
            ecc_antipode = 0
            for v in range(n):
                ecc_antipode = max(ecc_antipode, distances[v])
            ecc_upper_bound[antipode] = ecc_antipode

            # Update eccentricity lower bounds and remove from active those
            # vertices for which the gap is closed
            i = 0
            while i < len(active):
                v = active[i]
                ecc_lower_bound[v] = max(ecc_lower_bound[v], distances[v])
                if ecc_upper_bound[v] == ecc_lower_bound[v]:
                    active[i] = active[-1]
                    active.pop()
                else:
                    i += 1

    from sage.rings.infinity import Infinity
    cdef list eccentricity = []
    for i in range(n):
        if ecc_upper_bound[i] != sys.float_info.max:
            eccentricity.append(ecc_upper_bound[i])
        else:
            eccentricity.append(+Infinity)

    return eccentricity


cpdef radius_DHV(g, weight_function=None, check_weight=True):
    r"""
    Return the radius of weighted graph `g`.

    This method computes the radius of undirected graph using the algorithm
    given in [Dragan2018]_.

    This method returns Infinity if graph is not connected.

    INPUT:

    - ``g`` -- the input Sage graph

    - ``weight_function`` -- function (default: ``None``); a function that
      associates a weight to each edge. If ``None`` (default), the weights of
      ``g`` are used, if ``g.weighted()==True``, otherwise all edges have
      weight 1.

    - ``check_weight`` -- boolean (default: ``True``); if ``True``, we check
      that the ``weight_function`` outputs a number for each edge

    EXAMPLES::

        sage: from sage.graphs.base.boost_graph import radius_DHV
        sage: G = Graph([(0,1,1), (1,2,1), (0,2,3)])
        sage: radius_DHV(G)
        1.0
        sage: G = graphs.PathGraph(7)
        sage: radius_DHV(G) == G.radius(algorithm='Dijkstra_Boost')
        True

    TESTS::

        sage: G = Graph()
        sage: radius_DHV(G)
        0
        sage: G = Graph(1)
        sage: radius_DHV(G)
        0
        sage: G = Graph(2)
        sage: radius_DHV(G)
        +Infinity
        sage: G = Graph([(0, 1, 2)],weighted=True)
        sage: radius_DHV(G)
        2.0
        sage: G = DiGraph(1)
        sage: radius_DHV(G)
        Traceback (most recent call last):
        ...
        TypeError: this method works for undirected graphs only
    """
    if g.is_directed():
        raise TypeError("this method works for undirected graphs only")

    cdef int n = g.order()
    if n <= 1:
        return 0

    if weight_function and check_weight:
        g._check_weight_function(weight_function)

    if weight_function is not None:
        for e in g.edge_iterator():
            if float(weight_function(e)) < 0:
                raise ValueError("graphs contains negative weights, use Johnson_Boost instead")
    elif g.weighted():
        for _, _, w in g.edge_iterator():
            if w and float(w) < 0:
                raise ValueError("graphs contains negative weights, use Johnson_Boost instead")

    # These variables are automatically deleted when the function terminates.
    cdef dict v_to_int = {vv: vi for vi, vv in enumerate(g)}
    cdef BoostVecWeightedGraph g_boost
    boost_weighted_graph_from_sage_graph(&g_boost, g, v_to_int, weight_function)

    import sys
    cdef v_index source = 0
    cdef v_index antipode
    cdef v_index v
    cdef double ecc_source
    cdef double UB = sys.float_info.max
    cdef double LB = 0
    # For storing distances of all nodes from source
    cdef vector[double] distances
    # For storing lower bound on eccentricity of nodes
    cdef vector[double] ecc_lower_bound

    # Initializing
    for i in range(n):
        ecc_lower_bound.push_back(0)

    # Algorithm
    while LB < UB:
        # 1) pick vertex with minimum eccentricity lower bound
        # and compute its eccentricity
        sig_on()
        distances = g_boost.dijkstra_shortest_paths(source).distances
        sig_off()

        # Determine the eccentricity of source and its antipode, that is a
        # vertex at largest distance from source
        ecc_source = 0
        for v in range(n):
            if ecc_source < distances[v]:
                ecc_source = distances[v]
                antipode = v

        if ecc_source == sys.float_info.max:  # Disconnected graph
            break

        UB = min(UB, ecc_source)  # minimum among exact computed eccentricities
        if ecc_source == ecc_lower_bound[source]:
            # we have found minimum eccentricity vertex and hence the radius
            break

        # 2) Compute distances from antipode
        sig_on()
        distances = g_boost.dijkstra_shortest_paths(antipode).distances
        sig_off()

        # 3) Use distances from antipode to improve eccentricity lower bounds.
        # We also determine the next source
        LB = sys.float_info.max
        for v in range(n):
            ecc_lower_bound[v] = max(ecc_lower_bound[v], distances[v])
            if LB > ecc_lower_bound[v]:
                LB = ecc_lower_bound[v]
                source = v  # vertex with minimum eccentricity lower bound

    if UB == sys.float_info.max:
        from sage.rings.infinity import Infinity
        return +Infinity

    return UB


cpdef diameter_DHV(g, weight_function=None, check_weight=True):
    r"""
    Return the diameter of weighted graph `g`.

    This method computes the diameter of undirected graph using the
    algorithm proposed in [Dragan2018]_.

    This method returns Infinity if graph is not connected.

    INPUT:

    - ``g`` -- the input Sage graph

    - ``weight_function`` -- function (default: ``None``); a function that
      associates a weight to each edge. If ``None`` (default), the weights of
      ``g`` are used, if ``g.weighted()==True``, otherwise all edges have
      weight 1.

    - ``check_weight`` -- boolean (default: ``True``); if ``True``, we check
      that the ``weight_function`` outputs a number for each edge

    EXAMPLES::

        sage: from sage.graphs.base.boost_graph import diameter_DHV
        sage: G = graphs.ButterflyGraph()
        sage: diameter_DHV(G)
        2.0

    TESTS::

        sage: G = graphs.RandomBarabasiAlbert(17,6)                                     # needs networkx
        sage: diameter_DHV(G) == G.diameter(algorithm = 'Dijkstra_Boost')               # needs networkx
        True
        sage: G = Graph([(0,1,-1)], weighted=True)
        sage: diameter_DHV(G)
        Traceback (most recent call last):
        ...
        ValueError: graph contains negative edge weights, use Johnson_Boost instead
    """
    if g.is_directed():
        raise TypeError("this method works for undirected graphs only")

    cdef int n = g.order()
    if n <= 1:
        return 0

    if weight_function and check_weight:
        g._check_weight_function(weight_function)

    if weight_function is not None:
        for e in g.edges(sort=False):
            if float(weight_function(e)) < 0:
                raise ValueError("graph contains negative edge weights, use Johnson_Boost instead")
    elif g.weighted():
        for _, _, w in g.edges(sort=False):
            if w and float(w) < 0:
                raise ValueError("graph contains negative edge weights, use Johnson_Boost instead")

    # These variables are automatically deleted when the function terminates.
    cdef dict v_to_int = {vv: vi for vi, vv in enumerate(g)}
    cdef BoostVecWeightedGraph g_boost
    boost_weighted_graph_from_sage_graph(&g_boost, g, v_to_int, weight_function)

    import sys
    cdef v_index u, x, antipode
    cdef double ecc_u, ecc_x, ecc_antipode
    cdef double LB = 0
    cdef double UB = sys.float_info.max
    cdef v_index v
    cdef double tmp
    cdef size_t i, idx

    cdef list active = list(range(n))
    cdef vector[double] ecc_lower_bound, ecc_upper_bound, distances

    for i in range(n):
        ecc_lower_bound.push_back(0)
        ecc_upper_bound.push_back(sys.float_info.max)

    # Algorithm
    while LB < UB and active:
        # 1. Select vertex u with maximum eccentricity upper bound
        tmp = 0
        for i, v in enumerate(active):
            if ecc_upper_bound[v] > tmp:
                tmp = ecc_upper_bound[v]
                idx = i
        active[idx], active[-1] = active[-1], active[idx]
        u = active.pop()

        # Compute the distances from u
        sig_on()
        distances = g_boost.dijkstra_shortest_paths(u).distances
        sig_off()

        # compute the eccentricity of u and update eccentricity lower bounds
        ecc_u = 0
        for v in range(n):
            ecc_lower_bound[v] = max(ecc_lower_bound[v], distances[v])
            ecc_u = max(ecc_u, distances[v])

        LB = max(LB, ecc_u)

        if LB == sys.float_info.max:  # Disconnected graph
            break

        # 2. Select x such that dist(u, x) + ecc[x] == ecc[u].
        # Since we don't know ecc[x], we select x with minimum eccentricity
        # lower bound.  If ecc[x] == ecc_lb[x], we are done. Otherwise, we
        # update eccentricity lower bounds and repeat
        while active:
            # Select v with minimum eccentricity lower bound
            tmp = sys.float_info.max
            for i, v in enumerate(active):
                if ecc_lower_bound[v] < tmp:
                    tmp = ecc_lower_bound[v]
                    idx = i
            active[idx], active[-1] = active[-1], active[idx]
            x = active.pop()

            # compute the distances from x
            sig_on()
            distances = g_boost.dijkstra_shortest_paths(x).distances
            sig_off()

            # compute the eccentricity of x and its antipode
            ecc_x = 0
            for v in range(n):
                if distances[v] > ecc_x:
                    ecc_x = distances[v]
                    antipode = v
            LB = max(LB, ecc_x)

            if ecc_x == ecc_lower_bound[x]:
                # We found the good vertex x
                # We update eccentricity upper bounds and break
                UB = ecc_x
                for v in active:
                    ecc_upper_bound[v] = min(ecc_upper_bound[v], distances[v] + ecc_x)
                    UB = max(UB, ecc_upper_bound[v])
                break
            else:
                # x was not a good choice
                # We use its antipode to update eccentricity lower bounds.
                # Observe that this antipode might have already been seen.
                for i, v in enumerate(active):
                    if v == antipode:
                        active[i] = active[-1]
                        active.pop()
                        break

                # compute the distances from antipode
                sig_on()
                distances = g_boost.dijkstra_shortest_paths(antipode).distances
                sig_off()

                # compute the eccentricity of antipode and update
                # eccentricity lower bounds
                ecc_antipode = 0
                for v in range(n):
                    ecc_antipode = max(ecc_antipode, distances[v])
                    ecc_lower_bound[v] = max(ecc_lower_bound[v], distances[v])
                LB = max(LB, ecc_antipode)

    if LB == sys.float_info.max:
        from sage.rings.infinity import Infinity
        return +Infinity

    return LB

cdef tuple diameter_lower_bound_2Dsweep(BoostVecWeightedDiGraphU g_boost,
                                        BoostVecWeightedDiGraphU rev_g_boost,
                                        v_index source,
                                        str algorithm):
    r"""
    Return a lower bound on the diameter of `G`.

    This method implements the weighted version of the algorithm proposed
    in [Broder2000]_ to compute a lower bound on the diameter of the
    weighted digraph `G`.

    If the digraph is not strongly connected, the returned value is infinity.

    Firstly, this method computes forward distances from `source` and selects a
    vertex `vf` at maximum forward distance from `source` (i.e. an
    antipode). Then, it computes backward eccentricity of `vf`. Observe that the
    backward eccentricity of `vf` is at least the forward eccentricity of
    `source`.

    Secondly, this method computes backward distances from `source` and selects
    a vertex `vb` at maximum backward distance from `source`. Then, it computes
    the forward eccentricity of `vb`, which is at least the backward
    eccentricity of `source`.

    The lower bound on the diameter is the maximum among the backward
    eccentricity of `vf` and forward eccentricity of `vb`.

    The method returns `(LB, s, m, d)`, where `LB` is best found lower bound on
    diameter, `s` is a vertex whose forward / backward eccentricity is `LB`, `d`
    is a vertex at a distance `LB` from / to `s` , and `m` is a vertex at
    distance `LB/2` from / to both `s` and `d`.

    INPUT:

    - ``g_boost`` -- a boost weighted digraph

    - ``rev_g_boost`` -- a copy of ``g_boost`` with edges reversed

    - ``source`` -- starting node for forward and backward distance computation

    - ``algorithm`` -- string; algorithm for computing single source shortest
      distances. If ``g_boost`` contains negative edge weights then it will be
      ``Bellman-Ford``, otherwise it will be ``Dijkstra_Boost``.

    TESTS::

        sage: from sage.graphs.base.boost_graph import diameter
        sage: G = DiGraph()
        sage: diameter(DiGraph(), algorithm='2Dsweep')
        0
        sage: diameter(DiGraph(1), algorithm='2Dsweep')
        0
        sage: diameter(DiGraph(2), algorithm='2Dsweep')
        +Infinity
    """
    import sys

    cdef int n = g_boost.num_verts()

    if n <= 1:
        return (0, 0, 0, 0)

    cdef v_index source_1, source_2, m, s, d, antipode_1, antipode_2, v
    cdef double LB_1, LB_2, LB, LB_m
    cdef result_distances result_1, result_2
    source_1 = source_2 = source

    # Algorithm

    # 1) Compute forward distances from source_1.
    # Get forward eccentricity and antipode (vertex at maximum forward distance)
    if algorithm == 'Bellman-Ford':
        sig_on()
        result_1 = g_boost.bellman_ford_shortest_paths(source_1)
        sig_off()
    else:
        sig_on()
        result_1 = g_boost.dijkstra_shortest_paths(source_1)
        sig_off()
    if not result_1.distances.size():
        raise ValueError("the graph contains a negative cycle")

    LB_1 = -sys.float_info.max
    for v in range(n):
        if result_1.distances[v] > LB_1:
            LB_1 = result_1.distances[v]
            antipode_1 = v

    if LB_1 == sys.float_info.max:
        return (LB_1, 0, 0, 0)

    # 2) Compute backward distances from antipode_1.
    source_1, antipode_1 = antipode_1, source_1
    if algorithm == 'Bellman-Ford':
        sig_on()
        result_1 = rev_g_boost.bellman_ford_shortest_paths(source_1)
        sig_off()
    else:
        sig_on()
        result_1 = rev_g_boost.dijkstra_shortest_paths(source_1)
        sig_off()
    if not result_1.distances.size():
        raise ValueError("the graph contains a negative cycle")

    for v in range(n):
        if result_1.distances[v] > LB_1:
            LB_1 = result_1.distances[v]
            antipode_1 = v

    if LB_1 == sys.float_info.max:
        return (LB_1, 0, 0, 0)

    # 3) Compute backward distances from source_2.
    # Get backward eccentricity and antipode.
    if algorithm == 'Bellman-Ford':
        sig_on()
        result_2 = rev_g_boost.bellman_ford_shortest_paths(source_2)
        sig_off()
    else:
        sig_on()
        result_2 = rev_g_boost.dijkstra_shortest_paths(source_2)
        sig_off()
    if not result_2.distances.size():
        raise ValueError("the graph contains a negative cycle")

    LB_2 = -sys.float_info.max
    for v in range(n):
        if result_2.distances[v] > LB_2:
            LB_2 = result_2.distances[v]
            antipode_2 = v

    if LB_2 == sys.float_info.max:
        return (LB_2, 0, 0, 0)

    # 4) Compute forward distances from antipode_2.
    source_2, antipode_2 = antipode_2, source_2
    if algorithm == 'Bellman-Ford':
        sig_on()
        result_2 = g_boost.bellman_ford_shortest_paths(source_2)
        sig_off()
    else:
        sig_on()
        result_2 = g_boost.dijkstra_shortest_paths(source_2)
        sig_off()
    if not result_2.distances.size():
        raise ValueError("the graph contains a negative cycle")

    for v in range(n):
        if result_2.distances[v] > LB_2:
            LB_2 = result_2.distances[v]
            antipode_2 = v

    if LB_2 == sys.float_info.max:
        return (LB_2, 0, 0, 0)

    # 5) Select the best found lower bound as LB with corresponding source s and
    # antipode d. Then find a vertex m at a distance LB/2 from/to both s and d.
    if LB_1 < LB_2:
        LB = LB_2
        s = source_2
        d = antipode_2
        LB_m = LB_2 / 2
        m = d
        while result_2.distances[m] > LB_m:
            m = result_2.predecessors[m]
    else:
        LB = LB_1
        s = source_1
        d = antipode_1
        LB_m = LB_1 / 2
        m = d
        while result_1.distances[m] > LB_m:
            m = result_1.predecessors[m]

    return (LB, s, m, d)


cdef double diameter_DiFUB(BoostVecWeightedDiGraphU g_boost,
                           BoostVecWeightedDiGraphU rev_g_boost,
                           v_index source,
                           str algorithm) except? -1:
    r"""
    Return the diameter of a weighted directed graph.

    The ``DiFUB`` (Directed iterative Fringe Upper Bound) algorithm calculates
    the exact value of the diameter of a weighted directed graph [CGLM2012]_.

    This algorithm starts from a vertex found through a 2Dsweep call (a directed
    version of the 2sweep method). The worst case time complexity of the DiFUB
    algorithm is `O(nm)`, but it can be very fast in practice. See the code's
    documentation and [CGLM2012]_ for more details.

    If the digraph is not strongly connected, the returned value is infinity.

    INPUT:

    - ``g_boost`` -- a boost weighted digraph

    - ``rev_g_boost`` -- a copy of ``g_boost`` with edges reversed

    - ``source`` -- starting node for forward and backward distance computation

    - ``algorithm`` -- string; algorithm for computing single source shortest
      distances. If ``g_boost`` contains negative edge weights then it will be
      ``Bellman-Ford``, otherwise it will be ``Dijkstra_Boost``.

    TESTS::

        sage: from sage.graphs.base.boost_graph import diameter
        sage: G = DiGraph()
        sage: diameter(DiGraph(), algorithm='DiFUB')
        0
        sage: diameter(DiGraph(1), algorithm='DiFUB')
        0
        sage: diameter(DiGraph(2), algorithm='DiFUB')
        +Infinity
    """
    cdef v_index n = g_boost.num_verts()
    if n <= 1:
        return 0

    import sys
    # These variables are automatically deleted when the function terminates.
    cdef double LB, LB_1, LB_2, UB
    cdef v_index m, v, tmp
    cdef v_index i
    cdef vector[double] distances
    cdef vector[pair[double, v_index]] order_1, order_2

    # We select a vertex with low eccentricity using 2Dsweep
    LB, _, m, _ = diameter_lower_bound_2Dsweep(g_boost, rev_g_boost,
                                               source, algorithm)

    # If the lower bound is a very large number, it means that the digraph is
    # not strongly connected and so the diameter is infinite.
    if LB == sys.float_info.max:
        return LB

    # Compute Forward distances from `m`.
    if algorithm == 'Bellman-Ford':
        sig_on()
        distances = g_boost.bellman_ford_shortest_paths(m).distances
        sig_off()
    else:
        sig_on()
        distances = g_boost.dijkstra_shortest_paths(m).distances
        sig_off()
    if not distances.size():
        raise ValueError("the graph contains a negative cycle")

    # Obtain Forward eccentricity of `m` and store pair of
    # forward distances, vertex in order_1
    LB_1 = sys.float_info.min
    for v in range(n):
        LB_1 = max(LB_1, distances[v])
        order_1.push_back(pair[double, v_index](distances[v], v))
    # Compute Backward distances from `m`.
    if algorithm == 'Bellman-Ford':
        sig_on()
        distances = rev_g_boost.bellman_ford_shortest_paths(m).distances
        sig_off()
    else:
        sig_on()
        distances = rev_g_boost.dijkstra_shortest_paths(m).distances
        sig_off()
    if not distances.size():
        raise ValueError("the graph contains a negative cycle")

    # Obtain Backward eccentricity of `m` and store pair of
    # backward distances, vertex in order_2.
    LB_2 = sys.float_info.min
    for v in range(n):
        LB_2 = max(LB_2, distances[v])
        order_2.push_back(pair[double, v_index](distances[v], v))

    # Now sort order_1 / order_2 in decreasing order of forward / backward
    # distances respectively.
    # Now order_1 and order_2 will contain order of vertices in which
    # further distance computations will be done.
    order_1 = sorted(order_1, reverse=True)
    order_2 = sorted(order_2, reverse=True)

    LB = max(LB, LB_1, LB_2)
    if LB == sys.float_info.max:
        return LB

    # The algorithm:
    #
    # The diameter of the digraph is equal to the maximum forward or backward
    # eccentricity of a vertex. Let `\[db_1, db_2,..., db_i\]` represents the
    # different backward distances from `m` containing at least one vertex at
    # that distance. Similarly, let `\[df_1, df_2,..., df_i\]` represents the
    # different forward distances from `m` containing at least one vertex at
    # that distance.
    #
    # The algorithm is based on the following two observations:
    #
    # 1). All the nodes `x` at a backward distance greater than `\[db_i\]` from
    # `m` having forward eccentricity greater than `\[2db_{i-1}\]` have a
    # corresponding node `y` whose backward eccentricity is greater than or
    # equal to the forward eccentricity of `x`, at a forward distance greater
    # than `\[db_i\]` from `m`.
    #
    # 2). All the nodes `x` at a forward distance greater than `\[df_i\]` from
    # `m` having backward eccentricity greater than `\[2df_{i-1}\]` have a
    # corresponding node `y` whose forward eccentricity is greater than or equal
    # to the backward eccentricity of `x`, at a backward distance greater than
    # `\[df_i\]` from `m`.
    #
    # Therefore, we calculate backward / forward eccentricity of all nodes at
    # forward / backward distance `\[df_i / db_i\]` from `m` respectively. And
    # their maximum is `LB`. If `LB` is greater than `2(next maximum forward /
    # backward distance)` then we are done, else we proceed further.

    i = 0
    UB = max(2 * order_1[i].first, 2 * order_2[i].first)

    while LB < UB:
        v = order_1[i].second
        if algorithm == 'Bellman-Ford':
            sig_on()
            distances = rev_g_boost.bellman_ford_shortest_paths(v).distances
            sig_off()
        else:
            sig_on()
            distances = rev_g_boost.dijkstra_shortest_paths(v).distances
            sig_off()
        if not distances.size():
            raise ValueError("the graph contains a negative cycle")

        LB_1 = sys.float_info.min
        for tmp in range(n):
            LB_1 = max(LB_1, distances[tmp])

        v = order_2[i].second
        if algorithm == 'Bellman-Ford':
            sig_on()
            distances = g_boost.bellman_ford_shortest_paths(v).distances
            sig_off()
        else:
            sig_on()
            distances = g_boost.dijkstra_shortest_paths(v).distances
            sig_off()
        if not distances.size():
            raise ValueError("the graph contains a negative cycle")

        LB_2 = sys.float_info.min
        for tmp in range(n):
            LB_2 = max(LB_2, distances[tmp])

        # Update the lower bound
        LB = max(LB, LB_1, LB_2)
        i += 1

        if LB == sys.float_info.max or i == n:
            break

        # next maximum forward / backward distance
        UB = max(2 * order_1[i].first, 2 * order_2[i].first)

    # Finally return the computed diameter
    return LB

cpdef diameter(G, algorithm=None, source=None,
               weight_function=None, check_weight=True):
    r"""
    Return the diameter of `G`.

    This method returns Infinity if the digraph is not strongly connected. It
    can also quickly return a lower bound on the diameter using the ``2Dsweep``
    scheme.

    INPUT:

    - ``G`` -- the input sage digraph

    - ``algorithm`` -- string (default: ``None``); specifies the algorithm to
      use among:

      - ``'2Dsweep'`` -- computes lower bound on the diameter of a weighted
        directed graph using the weighted version of the algorithm proposed in
        [Broder2000]_. See the code's documentation for more details.

      - ``'DiFUB'`` -- computes the diameter of a weighted directed graph
        using the weighted version of the algorithm proposed in [CGLM2012]_.
        See the code's documentation for more details.

    - ``source`` -- (default: ``None``) vertex from which to start the
      computation. If ``source==None``, an arbitrary vertex of the graph is
      chosen. Raise an error if the initial vertex is not in `G`.

    - ``weight_function`` -- function (default: ``None``); a function that
      associates a weight to each edge. If ``None`` (default), the weights of
      ``G`` are used, if ``G.weighted()==True``, otherwise all edges have
      weight 1.

    - ``check_weight`` -- boolean (default: ``True``); if ``True``, we check
      that the ``weight_function`` outputs a number for each edge

    EXAMPLES::

        sage: from sage.graphs.base.boost_graph import diameter
        sage: G = DiGraph([(0, 1, 2), (1, 0, -1)])
        sage: diameter(G, algorithm='DiFUB')
        1.0
        sage: diameter(G, algorithm='DiFUB', weight_function=lambda e:e[2])
        2.0
        sage: G = DiGraph([(0, 1, -1), (1, 0, 2)])
        sage: diameter(G, algorithm='DiFUB', weight_function=lambda e:e[2])
        2.0

    TESTS:

    Diameter of weakly connected digraph is Infinity::

        sage: G = DiGraph(2)
        sage: diameter(G, algorithm='DiFUB')
        +Infinity
        sage: diameter(G, algorithm='2Dsweep')
        +Infinity

    DiGraph containing negative cycle::

        sage: G = DiGraph([(0,1,-2), (1,0,1)])
        sage: diameter(G, algorithm='2Dsweep', weight_function=lambda e:e[2])
        Traceback (most recent call last):
        ...
        ValueError: the graph contains a negative cycle
        sage: diameter(G, algorithm='DiFUB', weight_function=lambda e:e[2])
        Traceback (most recent call last):
        ...
        ValueError: the graph contains a negative cycle
    """
    import sys

    if not G.is_directed():
        raise TypeError("this method works only for digraphs")

    cdef int n = G.order()

    if n <= 1:
        return 0

    if weight_function and check_weight:
        G._check_weight_function(weight_function)

    # Algorithm for single source shortest distance computations.
    cdef str algo = 'Dijkstra_Boost'

    # If digraph contains negative edge weight then
    # algo is set to `Bellman-Ford`
    if weight_function is not None:
        for e in G.edges(sort=False):
            if float(weight_function(e)) < 0:
                algo = 'Bellman-Ford'
                break
    elif G.weighted():
        for _, _, w in G.edges(sort=False):
            if w and float(w) < 0:
                algo = 'Bellman-Ford'
                break

    if algorithm is None:  # default algorithm for diameter computation
        algorithm = 'DiFUB'

    if algorithm not in ['2Dsweep', 'DiFUB']:
        raise ValueError("unknown algorithm for computing the diameter of directed graph")

    if source is None:
        source = next(G.vertex_iterator())
    elif not G.has_vertex(source):
        raise ValueError("the specified source is not a vertex of the input Graph")

    # These variables are automatically deleted when the function terminates.
    cdef dict v_to_int = {vv: vi for vi, vv in enumerate(G)}

    # boost copy of G
    cdef BoostVecWeightedDiGraphU g_boost
    # boost copy of G with edges reversed
    cdef BoostVecWeightedDiGraphU rev_g_boost

    # Initializing
    boost_weighted_graph_from_sage_graph(&g_boost, G, v_to_int, weight_function)
    boost_weighted_graph_from_sage_graph(&rev_g_boost, G, v_to_int, weight_function, reverse=True)

    cdef v_index isource = 0 if source is None else v_to_int[source]
    cdef double LB

    if algorithm == '2Dsweep':
        LB = diameter_lower_bound_2Dsweep(g_boost, rev_g_boost, isource, algo)[0]
    else:
        LB = diameter_DiFUB(g_boost, rev_g_boost, isource, algo)

    if LB == sys.float_info.max:
        from sage.rings.infinity import Infinity
        return +Infinity
    else:
        return LB

cpdef shortest_paths_from_vertices(g, vertex_list=None, order=None,
                                   weight_function=None, algorithm=None):
    r"""
    Compute the shortest paths to all vertices from each vertex in
    ``vertex_list``.

    The input graph can be weighted: if the algorithm is Dijkstra, no negative
    weights are allowed, while if the algorithm is Bellman-Ford, negative
    weights are allowed, but there must be no negative cycle (otherwise, the
    shortest paths might not exist).

    However, Dijkstra algorithm is more efficient: for this reason, we suggest
    to use Bellman-Ford only if necessary (which is also the default option).

    The running-time for each vertex is `O(n \log n+m)` for Dijkstra algorithm
    and `O(mn)` for Bellman-Ford algorithm, where `n` is the number of nodes and
    `m` is the number of edges.

    INPUT:

    - ``g`` -- the input Sage graph

    - ``vertex_list`` -- list (default: ``None``); list of vertices to compute
      shortest paths from. By default (``None``), compute shortest paths from
      all vertices.

    - ``order`` -- list (default: ``None``); order of vertices of `g`

    - ``weight_function`` -- function (default: ``None``); a function that
      associates a weight to each edge. If ``None`` (default), the weights of
      ``g`` are used, if ``g.weighted()==True``, otherwise all edges have
      weight 1.

    - ``algorithm`` -- string (default: ``None``); one of the following
      algorithms:

      - ``'Dijkstra'``, ``'Dijkstra_Boost'`` -- the Dijkstra algorithm
        implemented in Boost (works only with positive weights)

      - ``'Bellman-Ford'``, ``'Bellman-Ford_Boost'`` -- the Bellman-Ford
        algorithm implemented in Boost (works also with negative weights,
        if there is no negative cycle)

    OUTPUT:

    The type of output depends on the input. More precisely -

    - A pair of dictionaries of list ``(distances, predecessors)``, when
      ``order is not None``, such that for each vertex ``v`` in ``vertex_list``,
      ``distances[v][i]`` store the shortest distance between ``v`` and
      ``order[i]`` and ``predecessors[v][i]`` store the last vertex in the
      shortest path from ``v`` to ``order[i]``.

    - A pair of dictionaries of dictionaries ``(distances, predecessors)`` such
      that for each vertex ``v`` in ``vertex_list``, ``distances[v]`` store the
      shortest distances of all the other vertices from ``v``,
      ``predecessors[v]`` store the last vertices in the shortest path from
      ``v`` to all the other vertices.

    EXAMPLES:

    Undirected graphs::

        sage: from sage.graphs.base.boost_graph import shortest_paths_from_vertices
        sage: g = Graph([(0,1,1),(1,2,2),(1,3,4),(2,3,1)], weighted=True)
        sage: shortest_paths_from_vertices(g,[1,2])
        ({1: {0: 1.0, 1: 0.0, 2: 2.0, 3: 3.0}, 2: {0: 3.0, 1: 2.0, 2: 0.0, 3: 1.0}},
         {1: {0: 1, 1: None, 2: 1, 3: 2}, 2: {0: 1, 1: 2, 2: None, 3: 2}})

    Directed graphs::

        sage: g = DiGraph([(0,1,1),(1,2,-1),(2,0,2),(2,3,1)], weighted=True)
        sage: shortest_paths_from_vertices(g,1)
        ({1: {0: 1.0, 1: 0.0, 2: -1.0, 3: 0.0}}, {1: {0: 2, 1: None, 2: 1, 3: 2}})
        sage: shortest_paths_from_vertices(g, 1, [0,1,2,3])
        ({1: [1.0, 0.0, -1.0, 0.0]}, {1: [2, None, 1, 2]})

    TESTS:

    Given an input which is not a graph::

        sage: shortest_paths_from_vertices("X-AE A-12", 1)
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage graph

    If there is a negative cycle::

        sage: g = DiGraph([(0,1,1),(1,2,-2),(2,0,0.5),(2,3,1)], weighted=True)
        sage: shortest_paths_from_vertices(g, 1)
        Traceback (most recent call last):
        ...
        ValueError: the graph contains a negative cycle

    If the given ordering is not valid::

        sage: g = DiGraph([(0,1,1),(1,2,2),(2,0,0.5),(2,3,1)], weighted=True)
        sage: shortest_paths_from_vertices(g,1,[0,1])
        Traceback (most recent call last):
        ...
        ValueError: Given ordering is not valid

    If Dijkstra is used with negative weights::

        sage: g = Graph([(0,1,1),(1,2,-2),(1,3,4)], weighted=True)
        sage: shortest_paths_from_vertices(g, 1, algorithm='Dijkstra')
        Traceback (most recent call last):
        ...
        RuntimeError: Dijkstra algorithm does not work with negative weights, use Bellman-Ford instead

    Wrong starting vertex::

        sage: shortest_paths_from_vertices(g, 55)
        Traceback (most recent call last):
        ...
        ValueError: the starting vertex 55 is not in the graph
    """
    import sys
    from sage.rings.infinity import Infinity
    from sage.graphs.generic_graph import GenericGraph

    if not isinstance(g, GenericGraph):
        raise TypeError("the input must be a Sage graph")

    if vertex_list is None:
        vertex_list = g

    else:
        if not isinstance(vertex_list, list):
            vertex_list = [vertex_list]

        for vertex in vertex_list:
            if vertex not in g:
                raise ValueError(f"the starting vertex {vertex} is not in the graph")

    if order is not None:
        if len(g) == len(order):
            for vertex in order:
                if vertex not in g:
                    raise ValueError("Given ordering is not valid")
        else:
            raise ValueError("Given ordering is not valid")

    cdef bint use_Bellman_Ford = algorithm in ['Bellman-Ford', 'Bellman-Ford_Boost']
    if not use_Bellman_Ford:
        # Check if there are edges with negative weights
        if weight_function is not None:
            for e in g.edges(sort=False):
                if float(weight_function(e)) < 0:
                    use_Bellman_Ford = True
                    break
        elif g.weighted():
            for _, _, wt in g.edges(sort=False):
                if float(wt) < 0:
                    use_Bellman_Ford = True
                    break

        if algorithm in ['Dijkstra', 'Dijkstra_Boost']:
            if use_Bellman_Ford:
                raise RuntimeError("Dijkstra algorithm does not work with "
                                   "negative weights, use Bellman-Ford instead")
        elif algorithm is not None:
            raise ValueError(f"unknown algorithm {algorithm!r}")

    # These variables are automatically deleted when the function terminates.
    cdef v_index vi, v, vert, pred, w
    cdef list int_to_v = list(g)
    cdef dict v_to_int = {vv: vi for vi, vv in enumerate(g)}
    cdef result_distances result
    cdef BoostVecWeightedDiGraphU g_boost_dir
    cdef BoostVecWeightedGraph g_boost_und
    cdef dict dist_v_dict, pred_v_dict
    cdef list dist_v_list, pred_v_list

    if g.is_directed():
        boost_weighted_graph_from_sage_graph(&g_boost_dir, g, v_to_int, weight_function)
    else:
        boost_weighted_graph_from_sage_graph(&g_boost_und, g, v_to_int, weight_function)

    distances = {}
    predecessors = {}

    for v in vertex_list:
        vi = v_to_int[v]
        if use_Bellman_Ford:
            if g.is_directed():
                sig_on()
                result = g_boost_dir.bellman_ford_shortest_paths(vi)
                sig_off()
            else:
                sig_on()
                result = g_boost_und.bellman_ford_shortest_paths(vi)
                sig_off()
            if not result.distances.size():
                raise ValueError("the graph contains a negative cycle")
        else:
            if g.is_directed():
                sig_on()
                result = g_boost_dir.dijkstra_shortest_paths(vi)
                sig_off()
            else:
                sig_on()
                result = g_boost_und.dijkstra_shortest_paths(vi)
                sig_off()
            if not result.distances.size():
                # This situation should never happen
                raise RuntimeError("something goes wrong. Please report the "
                                   "bug on sage-devel@googlegroups.com")

        if order is None:
            dist_v_dict = {}
            pred_v_dict = {}

            for vert in range(g.num_verts()):
                if result.distances[vert] != sys.float_info.max:
                    w = int_to_v[vert]
                    dist_v_dict[w] = result.distances[vert]
                    pred = result.predecessors[vert]
                    if pred == vert:
                        pred_v_dict[w] = None
                    else:
                        pred_v_dict[w] = int_to_v[pred]

            distances[v] = dist_v_dict
            predecessors[v] = pred_v_dict
        else:
            dist_v_list = []
            pred_v_list = []

            for w in order:
                vert = v_to_int[w]
                if result.distances[vert] != sys.float_info.max:
                    dist_v_list.append(result.distances[vert])
                    pred = result.predecessors[vert]
                    if pred == vert:
                        pred_v_list.append(None)
                    else:
                        pred_v_list.append(int_to_v[pred])

            distances[v] = dist_v_list
            predecessors[v] = pred_v_list

    return distances, predecessors

cpdef wiener_index(g, algorithm=None, weight_function=None, check_weight=True):
    r"""
    Return the Wiener index of the graph.

    The Wiener index of an undirected graph `G` is defined as
    `W(G) = \frac{1}{2} \sum_{u,v\in G} d(u,v)` where `d(u,v)` denotes the
    distance between vertices `u` and `v` (see [KRG1996]_).

    The Wiener index of a directed graph `G` is defined as the sum of the
    distances between each pairs of vertices, `W(G) = \sum_{u,v\in G} d(u,v)`.

    INPUT:

    - ``g`` -- the input Sage graph

    - ``algorithm`` -- string (default: ``None``); one of the following
      algorithms:

      - ``'Dijkstra'``, ``'Dijkstra_Boost'``: the Dijkstra algorithm implemented
        in Boost (works only with positive weights)

      - ``'Bellman-Ford'``, ``'Bellman-Ford_Boost'``: the Bellman-Ford algorithm
        implemented in Boost (works also with negative weights, if there is no
        negative cycle)

    - ``weight_function`` -- function (default: ``None``); a function that
      associates a weight to each edge. If ``None`` (default), the weights of
      ``g`` are used, if ``g.weighted()==True``, otherwise all edges have
      weight 1.

    - ``check_weight`` -- boolean (default: ``True``); if ``True``, we check
      that the ``weight_function`` outputs a number for each edge

    EXAMPLES:

        sage: from sage.graphs.base.boost_graph import wiener_index
        sage: g = Graph([(0,1,9), (1,2,7), (2,3,4), (3,0,3)])
        sage: wiener_index(g)
        8.0
        sage: g.weighted(True)
        sage: wiener_index(g)
        41.0

    Wiener index of circuit digraphs::

        sage: n = 10
        sage: g = digraphs.Circuit(n)
        sage: w = lambda x: (x*x*(x-1))/2
        sage: wiener_index(g) == w(n)
        True

    Wiener index of a graph of order 1::

        sage: wiener_index(Graph(1))
        0

    The Wiener index is not defined on the empty graph::

        sage: wiener_index(Graph())
        Traceback (most recent call last):
        ...
        ValueError: Wiener index is not defined for the empty graph

    TESTS:

    Using ``'Dijkstra'`` on a graph with negative weights::

        sage: g = Graph([(0, 1, -1), (1, 2, 1)])
        sage: def weight_of(e):
        ....:     return e[2]
        sage: wiener_index(g, algorithm='Dijkstra', weight_function=weight_of)
        Traceback (most recent call last):
        ...
        RuntimeError: Dijkstra algorithm does not work with negative weights, use Bellman-Ford instead

    Directed graph with a negative weight cycle::

        sage: g = DiGraph([(0, 1, -1), (1, 2, -1), (2, 0, -1)])
        sage: def weight_of(e):
        ....:     return e[2]
        sage: wiener_index(g, algorithm='Bellman-Ford', weight_function=weight_of)
        Traceback (most recent call last):
        ...
        ValueError: the graph contains a negative cycle
    """
    if not g:
        raise ValueError("Wiener index is not defined for the empty graph")

    cdef unsigned int n = g.order()
    if n == 1:
        return 0

    import sys

    if weight_function and check_weight:
        g._check_weight_function(weight_function)

    cdef bint use_Bellman_Ford = algorithm in ['Bellman-Ford', 'Bellman-Ford_Boost']
    if not use_Bellman_Ford:
        # Check if there are edges with negative weights
        if weight_function is not None:
            for e in g.edges(sort=False):
                if float(weight_function(e)) < 0:
                    use_Bellman_Ford = True
                    break
        elif g.weighted():
            for _, _, w in g.edges(sort=False):
                if float(w) < 0:
                    use_Bellman_Ford = True
                    break

        if algorithm in ['Dijkstra', 'Dijkstra_Boost']:
            if use_Bellman_Ford:
                raise RuntimeError("Dijkstra algorithm does not work with "
                                   "negative weights, use Bellman-Ford instead")
        elif algorithm is not None:
            raise ValueError(f"unknown algorithm {algorithm!r}")

    # These variables are automatically deleted when the function terminates.
    cdef v_index vi, u, v
    cdef dict v_to_int = {vv: vi for vi, vv in enumerate(g)}
    cdef BoostVecWeightedDiGraphU g_boost_dir
    cdef BoostVecWeightedGraph g_boost_und
    cdef vector[double] distances
    cdef double s = 0

    if g.is_directed():
        boost_weighted_graph_from_sage_graph(&g_boost_dir, g, v_to_int, weight_function)
    else:
        boost_weighted_graph_from_sage_graph(&g_boost_und, g, v_to_int, weight_function)

    for u in range(n):
        if use_Bellman_Ford:
            if g.is_directed():
                sig_on()
                distances = g_boost_dir.bellman_ford_shortest_paths(u).distances
                sig_off()
            else:
                sig_on()
                distances = g_boost_und.bellman_ford_shortest_paths(u).distances
                sig_off()
            if not distances.size():
                raise ValueError("the graph contains a negative cycle")
        else:
            if g.is_directed():
                sig_on()
                distances = g_boost_dir.dijkstra_shortest_paths(u).distances
                sig_off()
            else:
                sig_on()
                distances = g_boost_und.dijkstra_shortest_paths(u).distances
                sig_off()
            if not distances.size():
                # This situation should never happen
                raise RuntimeError("something goes wrong. Please report the "
                                   "bug on sage-devel@googlegroups.com")

        for v in range(0 if g.is_directed() else (u + 1), n):
            if distances[v] == sys.float_info.max:
                from sage.rings.infinity import Infinity
                return +Infinity
            else:
                s += distances[v]

    return s
