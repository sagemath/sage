#*****************************************************************************
#       Copyright (C) 2015 Michele Borassi michele.borassi@imtlucca.it
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
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

    :func:`clustering_coeff`  | Returns the clustering coefficient of all vertices in the graph.
    :func:`edge_connectivity` | Returns the edge connectivity of the graph.
    :func:`dominator_tree` | Returns a dominator tree of the graph.
    :func:`bandwidth_heuristics` | Uses heuristics to approximate the bandwidth of the graph.
    :func:`min_spanning_tree` | Computes a minimum spanning tree of a (weighted) graph.
    :func:`shortest_paths` | Uses Dijkstra or Bellman-Ford algorithm to compute the single-source shortest paths.
    :func:`johnson_shortest_paths` | Uses Johnson algorithm to compute the all-pairs shortest paths.
    :func:`johnson_closeness_centrality` | Uses Johnson algorithm to compute the closeness centrality of all vertices.

Functions
---------
"""

include "cysignals/signals.pxi"

cdef boost_graph_from_sage_graph(BoostGenGraph *g, g_sage):
    r"""
    Initializes the Boost graph ``g`` to be equal to ``g_sage``.

    The Boost graph ``*g`` must represent an empty graph (an exception is raised
    otherwise).
    """

    from sage.graphs.generic_graph import GenericGraph

    if not isinstance(g_sage, GenericGraph):
        raise ValueError("The input parameter must be a Sage graph.")

    if g.num_verts() > 0:
        raise ValueError("The Boost graph in input must be empty")

    N = g_sage.num_verts()
    cdef dict vertex_to_int = {v:i for i,v in enumerate(g_sage.vertices())}

    for i in range(N):
        g.add_vertex()

    for u,v in g_sage.edge_iterator(labels=None):
        g.add_edge(vertex_to_int[u], vertex_to_int[v])

cdef boost_weighted_graph_from_sage_graph(BoostWeightedGraph *g,
                                          g_sage,
                                          weight_function = None):
    r"""
    Initializes the Boost weighted graph ``g`` to be equal to ``g_sage``.

    The Boost graph ``*g`` must represent an empty weighted graph. The edge
    weights are chosen as follows, and they must be convertible to floats,
    otherwise an error is raised.

    - If ``weight_function`` is not ``None``, this function is used.

    - If ``weight_function`` is ``None`` and ``g`` is weighted, the edge labels
      of ``g`` are used; in other words, the weight of an edge ``e=(u,v,l)`` is
      ``l``.

    - Otherwise, all weights are set to 1.

    In particular, the ``weight_function`` must be a function which inputs an
    edge ``e`` and outputs a number.
    """

    from sage.graphs.generic_graph import GenericGraph

    if not isinstance(g_sage, GenericGraph):
        raise ValueError("The input parameter must be a Sage graph.")

    if g.num_verts() > 0:
        raise ValueError("The Boost graph in input must be empty")

    N = g_sage.num_verts()
    cdef dict vertex_to_int = {v:i for i,v in enumerate(g_sage.vertices())}

    for i in range(N):
        g.add_vertex()

    if weight_function is not None:
        for e in g_sage.edge_iterator():
            try:
                g.add_edge(vertex_to_int[e[0]],
                           vertex_to_int[e[1]],
                           float(weight_function(e)))

            except (ValueError, TypeError):
                raise ValueError("The weight function cannot find the" +
                                         " weight of " + str(e) + ".")

    elif g_sage.weighted():
        for u,v,w in g_sage.edge_iterator():
            try:
                g.add_edge(vertex_to_int[u], vertex_to_int[v], float(w))
            except (ValueError, TypeError):
                raise ValueError("The weight function cannot find the" +
                                         " weight of " + str((u,v,w)) + ".")
    else:
        for u,v in g_sage.edge_iterator(labels=False):
            g.add_edge(vertex_to_int[u], vertex_to_int[v], 1)



cdef boost_edge_connectivity(BoostVecGenGraph *g):
    r"""
    Computes the edge connectivity of the input Boost graph.

    The output is a pair ``[ec,edges]``, where ``ec`` is the edge connectivity,
    ``edges`` is the list of edges in a minimum cut.
    """
    result = g[0].edge_connectivity()

    cdef int i
    edges = [(result.edges[i], result.edges[i+1])
             for i in range(0, result.edges.size(), 2)]

    return [result.ec, edges]

cpdef edge_connectivity(g):
    r"""
    Computes the edge connectivity of the input graph, using Boost.

    The output is a pair ``[ec,edges]``, where ``ec`` is the edge connectivity,
    ``edges`` is the list of edges in a minimum cut.

    .. SEEALSO::

        :meth:`sage.graphs.generic_graph.GenericGraph.edge_connectivity`

    EXAMPLES:

    Computing the edge connectivity of a clique::

        sage: from sage.graphs.base.boost_graph import edge_connectivity
        sage: g = graphs.CompleteGraph(5)
        sage: edge_connectivity(g)
        [4, [(0, 1), (0, 2), (0, 3), (0, 4)]]

    Vertex-labeled graphs::

        sage: from sage.graphs.base.boost_graph import edge_connectivity
        sage: g = graphs.GridGraph([2,2])
        sage: edge_connectivity(g)
        [2, [((0, 0), (0, 1)), ((0, 0), (1, 0))]]

    """
    from sage.graphs.graph import Graph
    from sage.graphs.digraph import DiGraph

    # These variables are automatically deleted when the function terminates.
    cdef BoostVecGraph g_boost_und
    cdef BoostVecDiGraph g_boost_dir
    cdef list int_to_vertex = g.vertices()

    if isinstance(g, Graph):
        boost_graph_from_sage_graph(&g_boost_und, g)
        sig_check()
        ec, edges = boost_edge_connectivity(&g_boost_und)

    elif isinstance(g, DiGraph):
        from sage.misc.stopgap import stopgap
        stopgap("The edge connectivity of directed graphs is not implemented " +
                "in Boost. The result may be mathematically unreliable.",18753)

        boost_graph_from_sage_graph(&g_boost_dir, g)
        sig_check()
        ec, edges = boost_edge_connectivity(&g_boost_dir)

    else:
        raise ValueError("The input must be a Sage graph.")

    return [ec, [(int_to_vertex[u], int_to_vertex[v]) for (u,v) in edges]]

cdef boost_clustering_coeff(BoostGenGraph *g, vertices):
    r"""
    Computes the clustering coefficient of all vertices in the list provided.

    The output is a pair ``[average_clustering_coefficient, clust_of_v]``, where
    ``average_clustering_coefficient`` is the average clustering of the vertices
    in variable ``vertices``, ``clust_of_v`` is a dictionary that associates to
    each vertex (stored as an integer) its clustering coefficient.
    """
    cdef result_cc result
    cdef dict clust_of_v

    if len(vertices) == g.num_verts():
        result = g[0].clustering_coeff_all()
        clust_of_v = {v:result.clust_of_v[v] for v in range(g.num_verts())}
        return [result.average_clustering_coefficient, clust_of_v]

    else:
        clust_of_v = {v:g[0].clustering_coeff(v) for v in vertices}
        return [(sum(clust_of_v.itervalues()) / len(clust_of_v)), clust_of_v]


cpdef clustering_coeff(g, vertices = None):
    r"""
    Computes the clustering coefficient of the input graph, using Boost.

    The output is a pair ``[average_clustering_coefficient, clust_of_v]``, where
    ``average_clustering_coefficient`` is the average clustering of the vertices
    in variable ``vertices``, ``clust_of_v`` is a dictionary that associates to
    each vertex its clustering coefficient. If ``vertices`` is ``None``, all
    vertices are considered.

    .. SEEALSO::

        :meth:`sage.graphs.generic_graph.GenericGraph.clustering_coeff`

    INPUT:

    - ``g`` (Graph) - the input graph.

    - ``vertices`` (list) - the list of vertices we need to analyze (if
      ``None``, we will compute the clustering coefficient of all vertices).

    EXAMPLES:

    Computing the clustering coefficient of a clique::

        sage: from sage.graphs.base.boost_graph import clustering_coeff
        sage: g = graphs.CompleteGraph(5)
        sage: clustering_coeff(g)
        [1.0, {0: 1.0, 1: 1.0, 2: 1.0, 3: 1.0, 4: 1.0}]
        sage: clustering_coeff(g, vertices = [0,1,2])
        [1.0, {0: 1.0, 1: 1.0, 2: 1.0}]

    Of a non-clique graph with triangles::

        sage: g = graphs.IcosahedralGraph()
        sage: clustering_coeff(g, vertices=[1,2,3])
        [0.5, {1: 0.5, 2: 0.5, 3: 0.5}]

    With labels::

        sage: g.relabel(list("abcdefghiklm"))
        sage: clustering_coeff(g, vertices="abde")
        [0.5, {'a': 0.5, 'b': 0.5, 'd': 0.5, 'e': 0.5}]
    """
    from sage.graphs.graph import Graph

    sig_on()
    # These variables are automatically deleted when the function terminates.
    cdef BoostVecGraph g_boost
    cdef list g_vertices = g.vertices()
    cdef dict vertex_to_int = {v:i for i,v in enumerate(g_vertices)}

    if not isinstance(g, Graph):
        sig_off()
        raise ValueError("The input must be a Sage graph.")

    boost_graph_from_sage_graph(&g_boost, g)

    if vertices is None:
        vertices = g_vertices

    vertices_boost = [vertex_to_int[v] for v in vertices]
    [average_clustering, clust_v_int] = boost_clustering_coeff(&g_boost, vertices_boost)
    clust_v_sage = {g_vertices[v]: clust_v_int[v] for v in vertices_boost}
    sig_off()
    return [average_clustering, clust_v_sage]


cpdef dominator_tree(g, root, return_dict = False):
    r"""
    Uses Boost to compute the dominator tree of ``g``, rooted at ``root``.

    A node `d` dominates a node `n` if every path from the entry node
    ``root`` to `n` must go through `d`. The immediate dominator of a node
    `n` is the unique node that strictly dominates `n` but does not dominate
    any other node that dominates `n`. A dominator tree is a tree where each
    node's children are those nodes it immediately dominates. For more
    information, see :wikipedia:`Dominator_(graph_theory)`.

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

    - ``g`` (generic_graph) - the input graph.

    - ``root`` (vertex) - the root of the dominator tree.

    - ``return_dict`` (boolean) - if ``True``, the function returns a
      dictionary associating to each vertex its parent in the dominator
      tree. If ``False`` (default), it returns the whole tree, as a ``Graph``
      or a ``DiGraph``.

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

    If the output is a dictionary::

        sage: graphs.GridGraph([2,2]).dominator_tree((0,0), return_dict = True)
        {(0, 0): None, (0, 1): (0, 0), (1, 0): (0, 0), (1, 1): (0, 0)}

    TESTS:

    If ``g`` is not a graph, an error is raised::

        sage: from sage.graphs.base.boost_graph import dominator_tree
        sage: dominator_tree('I am not a graph', 0)
        Traceback (most recent call last):
        ...
        ValueError: The input g must be a Sage graph.

    If ``root`` is not a vertex, an error is raised::

        sage: digraphs.TransitiveTournament(10).dominator_tree('Not a vertex!')
        Traceback (most recent call last):
        ...
        ValueError: The input root must be a vertex of g.
        sage: graphs.GridGraph([2,2]).dominator_tree(0)
        Traceback (most recent call last):
        ...
        ValueError: The input root must be a vertex of g.

    """
    from sage.graphs.graph import Graph
    from sage.graphs.digraph import DiGraph

    if not isinstance(g, Graph) and not isinstance(g, DiGraph):
        raise ValueError("The input g must be a Sage graph.")
    if not root in g.vertices():
        raise ValueError("The input root must be a vertex of g.")

    sig_on()
    # These variables are automatically deleted when the function terminates.
    cdef BoostVecGraph g_boost_und
    cdef BoostVecDiGraph g_boost_dir
    cdef vector[v_index] result
    cdef dict vertex_to_int = {v:i for i,v in enumerate(g.vertices())}
    cdef list int_to_vertex = g.vertices()

    if isinstance(g, Graph):
        boost_graph_from_sage_graph(&g_boost_und, g)
        result = <vector[v_index]> g_boost_und.dominator_tree(vertex_to_int[root])

    elif isinstance(g, DiGraph):
        boost_graph_from_sage_graph(&g_boost_dir, g)
        result = <vector[v_index]> g_boost_dir.dominator_tree(vertex_to_int[root])

    sig_off()

    cdef v_index no_parent = -1

    if return_dict:
        return {v:(None if result[vertex_to_int[v]] == no_parent else int_to_vertex[<int> result[vertex_to_int[v]]]) for v in g.vertices()};

    edges = [[int_to_vertex[<int> result[vertex_to_int[v]]], v] for v in g.vertices() if result[vertex_to_int[v]] != no_parent]

    if g.is_directed():
        if len(edges) == 0:
            g = DiGraph()
            g.add_vertex(root)
            return g
        else:
            return DiGraph(edges)
    else:
        if len(edges) == 0:
            g = Graph()
            g.add_vertex(root)
            return g
        else:
            return Graph(edges)


cpdef bandwidth_heuristics(g, algorithm = 'cuthill_mckee'):
    r"""
    Uses Boost heuristics to approximate the bandwidth of the input graph.

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

    - ``g`` (``Graph``) - the input graph.

    - ``algorithm`` (``'cuthill_mckee'`` or ``'king'``) - the heuristic used to
      compute the ordering: Cuthill-McKee, or King.

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
        sage: bandwidth_heuristics(graphs.GridGraph([3,3]), algorithm='king')
        (3, [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2), (2, 1), (1, 2), (2, 2)])

    TESTS:

    Given an input which is not a graph::

        sage: from sage.graphs.base.boost_graph import bandwidth_heuristics
        sage: bandwidth_heuristics(digraphs.Path(10))
        Traceback (most recent call last):
        ...
        ValueError: The input g must be a Graph.
        sage: bandwidth_heuristics("I am not a graph!")
        Traceback (most recent call last):
        ...
        ValueError: The input g must be a Graph.

    Given a wrong algorithm::

        from sage.graphs.base.boost_graph import bandwidth_heuristics
        sage: bandwidth_heuristics(graphs.PathGraph(3), algorithm='tip top')
        Traceback (most recent call last):
        ...
        ValueError: Algorithm 'tip top' not yet implemented. Please contribute.

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
        raise ValueError("The input g must be a Graph.")
    if not algorithm in ['cuthill_mckee', 'king']:
        raise ValueError("Algorithm '%s' not yet implemented. Please contribute." %(algorithm))
    if g.num_edges()==0:
        return (0, g.vertices());

    sig_on()
    # These variables are automatically deleted when the function terminates.
    cdef BoostVecGraph g_boost
    cdef vector[v_index] result
    cdef dict vertex_to_int = {v:i for i,v in enumerate(g.vertices())}
    cdef list int_to_vertex = g.vertices()

    boost_graph_from_sage_graph(&g_boost, g)
    result = <vector[v_index]> g_boost.bandwidth_ordering(algorithm=='cuthill_mckee')

    cdef int n = g.num_verts()
    cdef dict pos = {int_to_vertex[<int> result[i]]:i for i in range(n)}
    cdef int bandwidth = max([abs(pos[u]-pos[v]) for u,v in g.edges(labels=False)])

    sig_off()
    return (bandwidth, [int_to_vertex[<int> result[i]] for i in range(n)])

cpdef min_spanning_tree(g,
                        weight_function=None,
                        algorithm='Kruskal'):
    r"""
    Uses Boost to compute the minimum spanning tree of the input graph.

    INPUT:

    - ``g`` (``Graph``) - the input graph.

    - ``weight_function`` (function) - a function that inputs an edge ``e`` and
      outputs its weight. An edge has the form ``(u,v,l)``, where ``u`` and
      ``v`` are vertices, ``l`` is a label (that can be of any kind). The
      ``weight_function`` can be used to transform the label into a weight (see
      the example below). In particular:

      - if ``weight_function`` is not ``None``, the weight of an edge ``e`` is
        ``weight_function(e)``;

      - if ``weight_function`` is ``None`` (default) and ``g`` is weighted (that
        is, ``g.weighted()==True``), for each edge ``e=(u,v,l)``, we set weight
        ``l``;

      - if ``weight_function`` is ``None`` and ``g`` is not weighted, we set all
        weights to 1 (hence, the output can be any spanning tree).

      Note that, if the weight is not convertible to a number with function
      ``float()``, an error is raised (see tests below).

    - ``algorithm`` (``'Kruskal'`` or ``'Prim'``) - the algorithm used.

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
        ValueError: The input g must be a Sage Graph.

    Given a wrong algorithm::

        sage: min_spanning_tree(graphs.PathGraph(3), algorithm='tip top')
        Traceback (most recent call last):
        ...
        ValueError: Algorithm 'tip top' not yet implemented. Please contribute.

    If the weight is not a number::

        sage: g = Graph([(0,1,1), (1,2,'a')], weighted=True)
        sage: min_spanning_tree(g)
        Traceback (most recent call last):
        ...
        ValueError: The weight function cannot find the weight of (1, 2, 'a').

        sage: g = Graph([(0,1,1), (1,2,[1,2,3])], weighted=True)
        sage: min_spanning_tree(g)
        Traceback (most recent call last):
        ...
        ValueError: The weight function cannot find the weight of (1, 2, [1, 2, 3]).

    """
    from sage.graphs.graph import Graph

    if not isinstance(g, Graph):
        raise ValueError("The input g must be a Sage Graph.")
    if not algorithm in ['Kruskal', 'Prim']:
        raise ValueError("Algorithm '%s' not yet implemented. Please contribute." %(algorithm))

    if g.allows_loops() or g.allows_multiple_edges():
        g = g.to_simple()
    # Now g has no self loops and no multiple edges.
    sig_on()
    # These variables are automatically deleted when the function terminates.
    cdef BoostVecWeightedGraph g_boost
    cdef vector[v_index] result
    cdef dict vertex_to_int = {v:i for i,v in enumerate(g.vertices())}
    cdef list int_to_vertex = g.vertices()

    try:
        boost_weighted_graph_from_sage_graph(&g_boost, g, weight_function)
    except Exception as e:
        sig_off()
        raise e

    if algorithm=='Kruskal':
        result = <vector[v_index]> g_boost.kruskal_min_spanning_tree()
    elif algorithm=='Prim':
        result = <vector[v_index]> g_boost.prim_min_spanning_tree()

    cdef int n = g.num_verts()
    sig_off()

    if result.size() != 2 * (n - 1):
        return []
    else:
        edges = [(int_to_vertex[<int> result[2*i]], int_to_vertex[<int> result[2*i+1]]) for i in range(n-1)]
        return sorted([(min(e[0],e[1]), max(e[0],e[1]), g.edge_label(e[0], e[1])) for e in edges])


cpdef shortest_paths(g, start, weight_function=None, algorithm=None):
    r"""
    Computes the shortest paths from ``start`` to all other vertices.

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

    - ``g`` (generic_graph) - the input graph.

    - ``start`` (vertex) - the starting vertex to compute shortest paths.

    - ``weight_function`` (function) - a function that associates a weight to
      each edge. If ``None`` (default), the weights of ``g`` are used, if
      available, otherwise all edges have weight 1.

    - ``algorithm`` (string) - one of the following algorithms:

      - ``'Dijkstra','Dijkstra_Boost'``: the Dijkstra algorithm implemented in
        Boost (works only with positive weights).

      - ``'Bellman-Ford','Bellman-Ford_Boost'``: the Bellman-Ford algorithm
        implemented in Boost (works also with negative weights, if there is no
        negative cycle).

    OUTPUT:

    A pair of dictionaries ``[distances, predecessors]`` such that, for each
    vertex ``v``, ``distances[v]`` is the distance from ``start`` to ``v``,
    ``predecessors[v]`` is the last vertex in a shortest path from ``start`` to
    ``v``.

    EXAMPLES:

    Undirected graphs::

        sage: from sage.graphs.base.boost_graph import shortest_paths
        sage: g = Graph([(0,1,1),(1,2,2),(1,3,4),(2,3,1)], weighted=True)
        sage: shortest_paths(g, 1)
        [{0: 1, 1: 0, 2: 2, 3: 3}, {0: 1, 1: None, 2: 1, 3: 2}]
        sage: g = graphs.GridGraph([2,2])
        sage: shortest_paths(g,(0,0),weight_function=lambda e:2)
        [{(0, 0): 0, (0, 1): 2, (1, 0): 2, (1, 1): 4},
         {(0, 0): None, (0, 1): (0, 0), (1, 0): (0, 0), (1, 1): (0, 1)}]

    Directed graphs::

        sage: g = DiGraph([(0,1,1),(1,2,2),(1,3,4),(2,3,1)], weighted=True)
        sage: shortest_paths(g, 1)
        [{1: 0, 2: 2, 3: 3}, {1: None, 2: 1, 3: 2}]

    TESTS:

    Given an input which is not a graph::

        sage: shortest_paths("I am not a graph!", 1)
        Traceback (most recent call last):
        ...
        ValueError: The input g must be a Sage Graph or DiGraph.

    If there is a negative cycle::

        sage: from sage.graphs.base.boost_graph import shortest_paths
        sage: g = DiGraph([(0,1,1),(1,2,-2),(2,0,0.5),(2,3,1)], weighted=True)
        sage: shortest_paths(g, 1)
        Traceback (most recent call last):
        ...
        ValueError: The graph contains a negative cycle.

    If Dijkstra is used with negative weights::

        sage: from sage.graphs.base.boost_graph import shortest_paths
        sage: g = DiGraph([(0,1,1),(1,2,-2),(1,3,4),(2,3,1)], weighted=True)
        sage: shortest_paths(g, 1, algorithm='Dijkstra')
        Traceback (most recent call last):
        ...
        ValueError: Dijkstra algorithm does not work with negative weights. Please, use Bellman-Ford.

    Wrong starting vartex::

        sage: shortest_paths(g, [])
        Traceback (most recent call last):
        ...
        ValueError: The starting vertex [] is not in the graph.
    """
    from sage.graphs.generic_graph import GenericGraph

    if not isinstance(g, GenericGraph):
        raise ValueError("The input g must be a Sage Graph or DiGraph.")
    elif g.num_edges() == 0:
        from sage.rings.infinity import Infinity
        return [{start:0}, {start:None}]


    # These variables are automatically deleted when the function terminates.
    cdef dict v_to_int = {v:i for i,v in enumerate(g.vertices())}
    cdef dict int_to_v = {i:v for i,v in enumerate(g.vertices())}
    cdef BoostVecWeightedDiGraphU g_boost_dir
    cdef BoostVecWeightedGraph g_boost_und
    cdef result_distances result

    if start not in v_to_int.keys():
        raise ValueError("The starting vertex " + str(start) + " is not in " +
                         "the graph.")

    if algorithm is None:
        # Check if there are edges with negative weights
        if weight_function is not None:
            for e in g.edge_iterator():
                try:
                    if float(weight_function(e)) < 0:
                        algorithm = 'Bellman-Ford'
                        break
                except (ValueError, TypeError):
                    raise ValueError("I cannot find the weight of edge " +
                                     str(e) + ".")

        else:
            for _,_,w in g.edge_iterator():
                try:
                    if float(w) < 0:
                        algorithm = 'Bellman-Ford'
                        break
                except (ValueError, TypeError):
                    raise ValueError("The label '", str(w), "' is not convertible " +
                                     "to a float.")

        if algorithm is None:
            algorithm = 'Dijkstra'

    if algorithm in ['Bellman-Ford', 'Bellman-Ford_Boost']:
        if g.is_directed():
            boost_weighted_graph_from_sage_graph(&g_boost_dir, g, weight_function)
            sig_on()
            result = g_boost_dir.bellman_ford_shortest_paths(v_to_int[start])
            sig_off()
        else:
            boost_weighted_graph_from_sage_graph(&g_boost_und, g, weight_function)
            sig_on()
            result = g_boost_und.bellman_ford_shortest_paths(v_to_int[start])
            sig_off()
        if result.distances.size() == 0:
            raise ValueError("The graph contains a negative cycle.");

    elif algorithm in ['Dijkstra', 'Dijkstra_Boost']:
        if g.is_directed():
            boost_weighted_graph_from_sage_graph(&g_boost_dir, g, weight_function)
            sig_on()
            result = g_boost_dir.dijkstra_shortest_paths(v_to_int[start])
            sig_off()
        else:
            boost_weighted_graph_from_sage_graph(&g_boost_und, g, weight_function)
            sig_on()
            result = g_boost_und.dijkstra_shortest_paths(v_to_int[start])
            sig_off()
        if result.distances.size() == 0:
            raise ValueError("Dijkstra algorithm does not work with negative weights. Please, use Bellman-Ford.");

    else:
        raise ValueError("Algorithm '%s' not yet implemented. Please contribute." %(algorithm))


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
    return [dist, pred]

cpdef johnson_shortest_paths(g, weight_function = None):
    r"""
    Uses Johnson algorithm to solve the all-pairs-shortest-paths.

    This routine outputs the distance between each pair of vertices, using a
    dictionary of dictionaries. It works on all kinds of graphs, but it is
    designed specifically for graphs with negative weights (otherwise there are
    more efficient algorithms, like Dijkstra).

    The time-complexity is `O(mn\log n)`, where `n` is the number of nodes and
    `m` is the number of edges.

    INPUT:

    - ``g`` (generic_graph) - the input graph.

    - ``weight_function`` (function) - a function that inputs an edge
      ``(u, v, l)`` and outputs its weight. If not ``None``, ``by_weight``
      is automatically set to ``True``. If ``None`` and ``by_weight`` is
      ``True``, we use the edge label ``l`` as a weight.

    OUTPUT:

    A dictionary of dictionary ``distances`` such that ``distances[v][w]`` is
    the distance between vertex ``v`` and vertex ``w``.

    EXAMPLES:

    Undirected graphs::

        sage: from sage.graphs.base.boost_graph import johnson_shortest_paths
        sage: g = Graph([(0,1,1),(1,2,2),(1,3,4),(2,3,1)], weighted=True)
        sage: johnson_shortest_paths(g)
        {0: {0: 0, 1: 1, 2: 3, 3: 4},
         1: {0: 1, 1: 0, 2: 2, 3: 3},
         2: {0: 3, 1: 2, 2: 0, 3: 1},
         3: {0: 4, 1: 3, 2: 1, 3: 0}}

    Directed graphs::

        sage: g = DiGraph([(0,1,1),(1,2,-2),(1,3,4),(2,3,1)], weighted=True)
        sage: johnson_shortest_paths(g)
        {0: {0: 0, 1: 1, 2: -1, 3: 0},
         1: {1: 0, 2: -2, 3: -1},
         2: {2: 0, 3: 1},
         3: {3: 0}}

    TESTS:

    Given an input which is not a graph::

        sage: johnson_shortest_paths("I am not a graph!")
        Traceback (most recent call last):
        ...
        ValueError: The input g must be a Sage Graph or DiGraph.

    If there is a negative cycle::

        sage: g = DiGraph([(0,1,1),(1,2,-2),(2,0,0.5),(2,3,1)], weighted=True)
        sage: johnson_shortest_paths(g)
        Traceback (most recent call last):
        ...
        ValueError: The graph contains a negative cycle.

    """
    from sage.graphs.generic_graph import GenericGraph

    if not isinstance(g, GenericGraph):
        raise ValueError("The input g must be a Sage Graph or DiGraph.")
    elif g.num_edges() == 0:
        from sage.rings.infinity import Infinity
        return {v:{v:0} for v in g.vertices()}
    # These variables are automatically deleted when the function terminates.
    cdef dict v_to_int = {v:i for i,v in enumerate(g.vertices())}
    cdef dict int_to_v = {i:v for i,v in enumerate(g.vertices())}
    cdef BoostVecWeightedDiGraphU g_boost_dir
    cdef BoostVecWeightedGraph g_boost_und
    cdef int N = g.num_verts()
    cdef vector[vector[double]] result

    if g.is_directed():
        boost_weighted_graph_from_sage_graph(&g_boost_dir, g, weight_function)
        sig_on()
        result = g_boost_dir.johnson_shortest_paths()
        sig_off()
    else:
        boost_weighted_graph_from_sage_graph(&g_boost_und, g, weight_function)
        sig_on()
        result = g_boost_und.johnson_shortest_paths()
        sig_off()

    if result.size() == 0:
        raise ValueError("The graph contains a negative cycle.")

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
    return {int_to_v[v]:{int_to_v[w]:correct_type(result[v][w])
                    for w in range(N) if result[v][w] != sys.float_info.max}
            for v in range(N)}

cpdef johnson_closeness_centrality(g, weight_function = None):
    r"""
    Uses Johnson algorithm to compute the closeness centrality of all vertices.

    This routine is preferrable to :func:`~johnson_shortest_paths` because it
    does not create a doubly indexed dictionary of distances, saving memory.

    The time-complexity is `O(mn\log n)`, where `n` is the number of nodes and
    `m` is the number of edges.

    INPUT:

    - ``g`` (generic_graph) - the input graph.

    - ``weight_function`` (function) - a function that inputs an edge
      ``(u, v, l)`` and outputs its weight. If not ``None``, ``by_weight``
      is automatically set to ``True``. If ``None`` and ``by_weight`` is
      ``True``, we use the edge label ``l`` as a weight.

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
        ValueError: The input g must be a Sage Graph or DiGraph.

    If there is a negative cycle::

        sage: from sage.graphs.base.boost_graph import johnson_closeness_centrality
        sage: g = DiGraph([(0,1,1),(1,2,-2),(2,0,0.5),(2,3,1)], weighted=True)
        sage: johnson_closeness_centrality(g)
        Traceback (most recent call last):
        ...
        ValueError: The graph contains a negative cycle.

    """
    from sage.graphs.generic_graph import GenericGraph

    if not isinstance(g, GenericGraph):
        raise ValueError("The input g must be a Sage Graph or DiGraph.")
    elif g.num_edges() == 0:
        from sage.rings.infinity import Infinity
        return {}
    sig_on()
    # These variables are automatically deleted when the function terminates.
    cdef BoostVecWeightedDiGraphU g_boost_dir
    cdef BoostVecWeightedGraph g_boost_und
    cdef int N = g.num_verts()
    cdef vector[vector[double]] result
    cdef vector[double] closeness
    cdef double farness
    cdef int i, j, reach

    if g.is_directed():
        boost_weighted_graph_from_sage_graph(&g_boost_dir, g, weight_function)
        result = g_boost_dir.johnson_shortest_paths()
    else:
        boost_weighted_graph_from_sage_graph(&g_boost_und, g, weight_function)
        result = g_boost_und.johnson_shortest_paths()

    if result.size() == 0:
        sig_off()
        raise ValueError("The graph contains a negative cycle.")

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
    sig_off()
    return {v: closeness[i] for i,v in enumerate(g.vertices()) if closeness[i] != sys.float_info.max}
