# cython: binding=True
r"""
Spanning trees

This module is a collection of algorithms on spanning trees. Also included in
the collection are algorithms for minimum spanning trees. See the book
[JNC2010]_ for descriptions of spanning tree algorithms,
including minimum spanning trees.

.. SEEALSO::

    - :meth:`GenericGraph.min_spanning_tree
      <sage.graphs.generic_graph.GenericGraph.min_spanning_tree>`.

.. TODO::

    - Parallel version of Boruvka's algorithm.


Methods
-------
"""

# ****************************************************************************
#       Copyright (c) 2007 Jason Grout <jason-sage@creativetrax.com>
#       Copyright (c) 2009 Mike Hansen <mhansen@gmail.com>
#       Copyright (c) 2010 Gregory McWhirter <gmcwhirt@uci.edu>
#       Copyright (c) 2010 Minh Van Nguyen <nguyenminh2@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from memory_allocator cimport MemoryAllocator
from sage.sets.disjoint_set cimport DisjointSet_of_hashables


def kruskal(G, by_weight=True, weight_function=None, check_weight=False, check=False):
    r"""
    Minimum spanning tree using Kruskal's algorithm.

    This function assumes that we can only compute minimum spanning trees for
    undirected graphs. Such graphs can be weighted or unweighted, and they can
    have multiple edges (since we are computing the minimum spanning tree, only
    the minimum weight among all `(u,v)`-edges is considered, for each pair
    of vertices `u`, `v`).

    INPUT:

    - ``G`` -- an undirected graph

    - ``by_weight`` -- boolean (default: ``True``); if ``True``, the edges in
      the graph are weighted. If ``False``, all edges have weight 1.

    - ``weight_function`` -- function (default: ``None``); a function that takes
      as input an edge ``(u, v, l)`` and outputs its weight. If not ``None``,
      ``by_weight`` is automatically set to ``True``. If ``None`` and
      ``by_weight`` is ``True``, we use the edge label ``l``, if ``l`` is not
      ``None``, else ``1`` as a weight.

    - ``check_weight`` -- boolean (default: ``False``); whether to check that
      the ``weight_function`` outputs a number for each edge

    - ``check`` -- boolean (default: ``False``); whether to first perform sanity
      checks on the input graph ``G``. Default: ``check=False``. If we toggle
      ``check=True``, the following sanity checks are first performed on ``G``
      prior to running Kruskal's algorithm on that input graph:

      - Is ``G`` the null graph?
      - Is ``G`` disconnected?
      - Is ``G`` a tree?
      - Does ``G`` have self-loops?
      - Does ``G`` have multiple edges?

      By default, we turn off the sanity checks for performance reasons. This
      means that by default the function assumes that its input graph is
      connected, and has at least one vertex. Otherwise, you should set
      ``check=True`` to perform some sanity checks and preprocessing on the
      input graph. If ``G`` has multiple edges or self-loops, the algorithm
      still works, but the running-time can be improved if these edges are
      removed. To further improve the runtime of this function, you should call
      it directly instead of using it indirectly via
      :meth:`sage.graphs.generic_graph.GenericGraph.min_spanning_tree`.

    OUTPUT:

    The edges of a minimum spanning tree of ``G``, if one exists, otherwise
    returns the empty list.

    .. SEEALSO::

        - :meth:`sage.graphs.generic_graph.GenericGraph.min_spanning_tree`
        - :func:`kruskal_iterator`
        - :func:`filter_kruskal` and :func:`filter_kruskal_iterator`

    EXAMPLES:

    An example from pages 727--728 in [Sah2000]_. ::

        sage: from sage.graphs.spanning_tree import kruskal
        sage: G = Graph({1:{2:28, 6:10}, 2:{3:16, 7:14}, 3:{4:12}, 4:{5:22, 7:18}, 5:{6:25, 7:24}})
        sage: G.weighted(True)
        sage: E = kruskal(G, check=True); E
        [(1, 6, 10), (3, 4, 12), (2, 7, 14), (2, 3, 16), (4, 5, 22), (5, 6, 25)]

    Variants of the previous example. ::

        sage: H = Graph(G.edges(sort=True, labels=False))
        sage: kruskal(H, check=True)
        [(1, 2, None), (1, 6, None), (2, 3, None), (2, 7, None), (3, 4, None), (4, 5, None)]
        sage: G.allow_loops(True)
        sage: G.allow_multiple_edges(True)
        sage: G
        Looped multi-graph on 7 vertices
        sage: for i in range(20):
        ....:     u = randint(1, 7)
        ....:     v = randint(1, 7)
        ....:     w = randint(0, 20)
        ....:     G.add_edge(u, v, w)
        sage: H = copy(G)
        sage: H
        Looped multi-graph on 7 vertices
        sage: def sanitize(G):
        ....:     G.allow_loops(False)
        ....:     G.allow_multiple_edges(False, keep_label='min')
        sage: sanitize(H)
        sage: H
        Graph on 7 vertices
        sage: sum(e[2] for e in kruskal(G, check=True)) == sum(e[2] for e in kruskal(H, check=True))
        True

    An example from pages 599--601 in [GT2001]_. ::

        sage: G = Graph({"SFO":{"BOS":2704, "ORD":1846, "DFW":1464, "LAX":337},
        ....: "BOS":{"ORD":867, "JFK":187, "MIA":1258},
        ....: "ORD":{"PVD":849, "JFK":740, "BWI":621, "DFW":802},
        ....: "DFW":{"JFK":1391, "MIA":1121, "LAX":1235},
        ....: "LAX":{"MIA":2342},
        ....: "PVD":{"JFK":144},
        ....: "JFK":{"MIA":1090, "BWI":184},
        ....: "BWI":{"MIA":946}})
        sage: G.weighted(True)
        sage: kruskal(G, check=True)
        [('JFK', 'PVD', 144),
         ('BWI', 'JFK', 184),
         ('BOS', 'JFK', 187),
         ('LAX', 'SFO', 337),
         ('BWI', 'ORD', 621),
         ('DFW', 'ORD', 802),
         ('BWI', 'MIA', 946),
         ('DFW', 'LAX', 1235)]

    An example from pages 568--569 in [CLRS2001]_. ::

        sage: G = Graph({"a":{"b":4, "h":8}, "b":{"c":8, "h":11},
        ....: "c":{"d":7, "f":4, "i":2}, "d":{"e":9, "f":14},
        ....: "e":{"f":10}, "f":{"g":2}, "g":{"h":1, "i":6}, "h":{"i":7}})
        sage: G.weighted(True)
        sage: T = Graph(kruskal(G, check=True), format='list_of_edges')
        sage: sum(T.edge_labels())
        37
        sage: T.is_tree()
        True

    An example with custom edge labels::

        sage: G = Graph([[0,1,1],[1,2,1],[2,0,10]], weighted=True)
        sage: weight = lambda e:3-e[0]-e[1]
        sage: sorted(kruskal(G, check=True))
        [(0, 1, 1), (1, 2, 1)]
        sage: sorted(kruskal(G, weight_function=weight, check=True))
        [(0, 2, 10), (1, 2, 1)]
        sage: sorted(kruskal(G, weight_function=weight, check=False))
        [(0, 2, 10), (1, 2, 1)]

    TESTS:

    The input graph must not be empty. ::

        sage: from sage.graphs.spanning_tree import kruskal
        sage: kruskal(graphs.EmptyGraph(), check=True)
        []
        sage: kruskal(Graph(), check=True)
        []
        sage: kruskal(Graph(multiedges=True), check=True)
        []
        sage: kruskal(Graph(loops=True), check=True)
        []
        sage: kruskal(Graph(multiedges=True, loops=True), check=True)
        []

    The input graph must be connected. ::

        sage: # long time
        sage: def my_disconnected_graph(n, ntries, directed=False, multiedges=False, loops=False):
        ....:     G = Graph()
        ....:     k = randint(2, n)
        ....:     G.add_vertices(range(k))
        ....:     if directed:
        ....:         G = G.to_directed()
        ....:     if multiedges:
        ....:         G.allow_multiple_edges(True)
        ....:     if loops:
        ....:         G.allow_loops(True)
        ....:     for i in range(ntries):
        ....:         u = randint(0, k-1)
        ....:         v = randint(0, k-1)
        ....:         if u != v or loops:
        ....:             G.add_edge(u, v)
        ....:     while G.is_connected():
        ....:         u = randint(0, k-1)
        ....:         v = randint(0, k-1)
        ....:         G.delete_edge(u, v)
        ....:     return G
        sage: G = my_disconnected_graph(100, 50, directed=False, multiedges=False, loops=False)
        sage: kruskal(G, check=True)
        []
        sage: G = my_disconnected_graph(100, 50, directed=False, multiedges=True, loops=False)
        sage: kruskal(G, check=True)
        []
        sage: G = my_disconnected_graph(100, 50, directed=False, multiedges=True, loops=True)
        sage: kruskal(G, check=True)
        []

    If the input graph is a tree, then return its edges::

        sage: T = graphs.RandomTree(randint(1, 50))  # long time
        sage: sorted(T.edge_iterator()) == sorted(kruskal(T, check=True))  # long time
        True

    If the input is not a Graph::

        sage: kruskal("I am not a graph")
        Traceback (most recent call last):
        ...
        ValueError: the input graph must be undirected
        sage: kruskal(digraphs.Path(10))
        Traceback (most recent call last):
        ...
        ValueError: the input graph must be undirected

    Check that the method is robust to incomparable vertices::

        sage: G = Graph([(1, 2, 10), (1, 'a', 1), ('a', 'b', 1), ('b', 2, 1)])
        sage: E = kruskal(G, by_weight=True)
        sage: sum(w for _, _, w in E)
        3
    """
    return list(kruskal_iterator(G, by_weight=by_weight, weight_function=weight_function,
                                 check_weight=check_weight, check=check))


def kruskal_iterator(G, by_weight=True, weight_function=None, check_weight=False, bint check=False):
    """
    Return an iterator implementation of Kruskal algorithm.

    INPUT:

    - ``G`` -- an undirected graph

    - ``by_weight`` -- boolean (default: ``True``); if ``True``, the edges in
      the graph are weighted. If ``False``, all edges have weight 1.

    - ``weight_function`` -- function (default: ``None``); a function that takes
      as input an edge ``(u, v, l)`` and outputs its weight. If not ``None``,
      ``by_weight`` is automatically set to ``True``. If ``None`` and
      ``by_weight`` is ``True``, we use the edge label ``l``, if ``l`` is not
      ``None``, else ``1`` as a weight.

    - ``check_weight`` -- boolean (default: ``False``); whether to check that
      the ``weight_function`` outputs a number for each edge

    - ``check`` -- boolean (default: ``False``); whether to first perform sanity
      checks on the input graph ``G``. Default: ``check=False``. If we toggle
      ``check=True``, the following sanity checks are first performed on ``G``
      prior to running Kruskal's algorithm on that input graph:

      - Is ``G`` the null graph?
      - Is ``G`` disconnected?
      - Is ``G`` a tree?
      - Does ``G`` have self-loops?
      - Does ``G`` have multiple edges?

      By default, we turn off the sanity checks for performance reasons. This
      means that by default the function assumes that its input graph is
      connected, and has at least one vertex. Otherwise, you should set
      ``check=True`` to perform some sanity checks and preprocessing on the
      input graph. If ``G`` has multiple edges or self-loops, the algorithm
      still works, but the running-time can be improved if these edges are
      removed. To further improve the runtime of this function, you should call
      it directly instead of using it indirectly via
      :meth:`sage.graphs.generic_graph.GenericGraph.min_spanning_tree`.

    OUTPUT: the edges of a minimum spanning tree of ``G``, one by one

    .. SEEALSO:: :func:`kruskal`

    EXAMPLES::

        sage: from sage.graphs.spanning_tree import kruskal_iterator
        sage: G = Graph({1:{2:28, 6:10}, 2:{3:16, 7:14}, 3:{4:12}, 4:{5:22, 7:18}, 5:{6:25, 7:24}})
        sage: G.weighted(True)
        sage: next(kruskal_iterator(G, check=True))
        (1, 6, 10)

    TESTS:

    If the input is not a Graph::

        sage: list(kruskal_iterator("I am not a graph"))
        Traceback (most recent call last):
        ...
        ValueError: the input graph must be undirected
        sage: list(kruskal_iterator(digraphs.Path(2)))
        Traceback (most recent call last):
        ...
        ValueError: the input graph must be undirected

    Check that the method is robust to incomparable vertices::

        sage: G = Graph([(1, 2, 10), (1, 'a', 1), ('a', 'b', 1), ('b', 2, 1)])
        sage: E = list(kruskal_iterator(G, by_weight=True))
        sage: sum(w for _, _, w in E)
        3
    """
    from sage.graphs.graph import Graph
    if not isinstance(G, Graph):
        raise ValueError("the input graph must be undirected")

    # sanity checks
    if check:
        if not G.order():
            return
        if not G.is_connected():
            return
        # G is now assumed to be a nonempty connected graph
        if G.num_verts() == G.num_edges() + 1:
            # G is a tree
            yield from G.edge_iterator()
            return

    cdef DisjointSet_of_hashables union_find = DisjointSet_of_hashables(G)
    by_weight, weight_function = G._get_weight_function(by_weight=by_weight,
                                                        weight_function=weight_function,
                                                        check_weight=check_weight)
    yield from kruskal_iterator_from_edges(G.edge_iterator(), union_find,
                                           by_weight=by_weight,
                                           weight_function=weight_function,
                                           check_weight=False)


def kruskal_iterator_from_edges(edges, union_find, by_weight=True,
                                weight_function=None, check_weight=False):
    """
    Return an iterator implementation of Kruskal algorithm on list of edges.

    INPUT:

    - ``edges`` -- list of edges

    - ``union_find`` -- a
      :class:`~sage.sets.disjoint_set.DisjointSet_of_hashables` encoding a
      forest

    - ``by_weight`` -- boolean (default: ``True``); if ``True``, the edges in
      the graph are weighted. If ``False``, all edges have weight 1.

    - ``weight_function`` -- function (default: ``None``); a function that takes
      as input an edge ``(u, v, l)`` and outputs its weight. If not ``None``,
      ``by_weight`` is automatically set to ``True``. If ``None`` and
      ``by_weight`` is ``True``, we use the edge label ``l``, if ``l`` is not
      ``None``, else ``1`` as a weight.

    - ``check_weight`` -- boolean (default: ``False``); whether to check that
      the ``weight_function`` outputs a number for each edge

    OUTPUT: the edges of a minimum spanning tree of ``G``, one by one

    .. SEEALSO::

        - :func:`kruskal`
        - :func:`filter_kruskal`

    EXAMPLES::

        sage: from sage.graphs.spanning_tree import kruskal_iterator_from_edges
        sage: G = Graph({1:{2:28, 6:10}, 2:{3:16, 7:14}, 3:{4:12}, 4:{5:22, 7:18}, 5:{6:25, 7:24}})
        sage: G.weighted(True)
        sage: union_set = DisjointSet(G)
        sage: next(kruskal_iterator_from_edges(G.edges(sort=False), union_set, by_weight=G.weighted()))
        (1, 6, 10)

    Check that the method is robust to incomparable vertices::

        sage: G = Graph([(1, 2, 10), (1, 'a', 1), ('a', 'b', 1), ('b', 2, 1)])
        sage: union_set = DisjointSet(G)
        sage: E = list(kruskal_iterator_from_edges(G.edges(sort=False), union_set, by_weight=True))
        sage: sum(w for _, _, w in E)
        3
    """
    # We sort edges, as specified.
    if weight_function is not None:
        edges = sorted(edges, key=weight_function)
    elif by_weight:
        from operator import itemgetter
        edges = sorted(edges, key=itemgetter(2))

    # Kruskal's algorithm
    for e in edges:
        # acyclic test via union-find
        u = union_find.find(e[0])
        v = union_find.find(e[1])
        if u != v:
            yield e
            # merge the trees
            union_find.union(u, v)
            if union_find.number_of_subsets() == 1:
                return


def filter_kruskal(G, threshold=10000, by_weight=True, weight_function=None,
                   check_weight=True, bint check=False):
    """
    Minimum spanning tree using Filter Kruskal algorithm.

    This function implements the variant of Kruskal's algorithm proposed in
    [OSS2009]_. Instead of directly sorting the whole set of edges, it
    partitions it in a similar way to quicksort and filter out edges that
    connect vertices of the same tree to reduce the cost of sorting.

    This function assumes that we can only compute minimum spanning trees for
    undirected graphs. Such graphs can be weighted or unweighted, and they can
    have multiple edges (since we are computing the minimum spanning tree, only
    the minimum weight among all `(u,v)`-edges is considered, for each pair of
    vertices `u`, `v`).

    INPUT:

    - ``G`` -- an undirected graph

    - ``threshold`` -- integer (default: 10000); maximum number of edges on
      which to run kruskal algorithm. Above that value, edges are partitioned
      into sets of size at most ``threshold``

    - ``by_weight`` -- boolean (default: ``True``); if ``True``, the edges in
      the graph are weighted. If ``False``, all edges have weight 1.

    - ``weight_function`` -- function (default: ``None``); a function that takes
      as input an edge ``(u, v, l)`` and outputs its weight. If not ``None``,
      ``by_weight`` is automatically set to ``True``. If ``None`` and
      ``by_weight`` is ``True``, we use the edge label ``l``, if ``l`` is not
      ``None``, else ``1`` as a weight.

    - ``check_weight`` -- boolean (default: ``False``); whether to check that
      the ``weight_function`` outputs a number for each edge

    - ``check`` -- boolean (default: ``False``); whether to first perform sanity
      checks on the input graph ``G``. Default: ``check=False``. If we toggle
      ``check=True``, the following sanity checks are first performed on ``G``
      prior to running Kruskal's algorithm on that input graph:

      - Is ``G`` the null graph?
      - Is ``G`` disconnected?
      - Is ``G`` a tree?
      - Does ``G`` have self-loops?
      - Does ``G`` have multiple edges?

    OUTPUT:

    The edges of a minimum spanning tree of ``G``, if one exists, otherwise
    returns the empty list.

    .. SEEALSO::

        - :meth:`sage.graphs.generic_graph.GenericGraph.min_spanning_tree`
        - :wikipedia:`Kruskal%27s_algorithm`
        - :func:`kruskal`
        - :func:`filter_kruskal_iterator`

    EXAMPLES::

        sage: from sage.graphs.spanning_tree import filter_kruskal
        sage: G = Graph({1:{2:28, 6:10}, 2:{3:16, 7:14}, 3:{4:12}, 4:{5:22, 7:18}, 5:{6:25, 7:24}})
        sage: G.weighted(True)
        sage: filter_kruskal(G, check=True)
        [(1, 6, 10), (3, 4, 12), (2, 7, 14), (2, 3, 16), (4, 5, 22), (5, 6, 25)]

        sage: filter_kruskal(Graph(2), check=True)
        []

    TESTS:

    Check that the method is robust to incomparable vertices::

        sage: G = Graph([(1, 2, 10), (1, 'a', 1), ('a', 'b', 1), ('b', 2, 1)])
        sage: E = filter_kruskal(G, by_weight=True)
        sage: sum(w for _, _, w in E)
        3
    """
    return list(filter_kruskal_iterator(G, threshold=threshold,
                                        by_weight=by_weight, weight_function=weight_function,
                                        check_weight=check_weight, check=check))


def filter_kruskal_iterator(G, threshold=10000, by_weight=True, weight_function=None,
                            check_weight=True, bint check=False):
    r"""
    Return an iterator implementation of Filter Kruskal's algorithm.

    INPUT:

    - ``G`` -- an undirected graph

    - ``threshold`` -- integer (default: 10000); maximum number of edges on
      which to run kruskal algorithm. Above that value, edges are partitioned
      into sets of size at most ``threshold``

    - ``by_weight`` -- boolean (default: ``True``); if ``True``, the edges in
      the graph are weighted. If ``False``, all edges have weight 1.

    - ``weight_function`` -- function (default: ``None``); a function that takes
      as input an edge ``(u, v, l)`` and outputs its weight. If not ``None``,
      ``by_weight`` is automatically set to ``True``. If ``None`` and
      ``by_weight`` is ``True``, we use the edge label ``l``, if ``l`` is not
      ``None``, else ``1`` as a weight.

    - ``check_weight`` -- boolean (default: ``False``); whether to check that
      the ``weight_function`` outputs a number for each edge

    - ``check`` -- boolean (default: ``False``); whether to first perform sanity
      checks on the input graph ``G``. Default: ``check=False``. If we toggle
      ``check=True``, the following sanity checks are first performed on ``G``
      prior to running Kruskal's algorithm on that input graph:

      - Is ``G`` the null graph?
      - Is ``G`` disconnected?
      - Is ``G`` a tree?
      - Does ``G`` have self-loops?
      - Does ``G`` have multiple edges?

    OUTPUT: the edges of a minimum spanning tree of ``G``, one by one

    .. SEEALSO::

        - :meth:`sage.graphs.generic_graph.GenericGraph.min_spanning_tree`
        - :wikipedia:`Kruskal%27s_algorithm`
        - :func:`kruskal`
        - :func:`filter_kruskal`

    EXAMPLES:

    The edges of a minimum spanning tree of ``G``, if one exists, otherwise
    returns the empty list. ::

        sage: from sage.graphs.spanning_tree import filter_kruskal_iterator
        sage: G = Graph({1:{2:28, 6:10}, 2:{3:16, 7:14}, 3:{4:12}, 4:{5:22, 7:18}, 5:{6:25, 7:24}})
        sage: G.weighted(True)
        sage: list(filter_kruskal_iterator(G, threshold=3, check=True))
        [(1, 6, 10), (3, 4, 12), (2, 7, 14), (2, 3, 16), (4, 5, 22), (5, 6, 25)]

    The weights of the spanning trees returned by :func:`kruskal_iterator` and
    :func:`filter_kruskal_iterator` are the same::

        sage: # needs networkx
        sage: from sage.graphs.spanning_tree import kruskal_iterator
        sage: G = graphs.RandomBarabasiAlbert(50, 2)
        sage: for u, v in G.edge_iterator(labels=False):
        ....:     G.set_edge_label(u, v, randint(1, 10))
        sage: G.weighted(True)
        sage: sum(e[2] for e in kruskal_iterator(G)) == sum(e[2]
        ....:     for e in filter_kruskal_iterator(G, threshold=20))
        True

    TESTS:

    The threshold must be at least 1::

        sage: from sage.graphs.spanning_tree import filter_kruskal_iterator
        sage: next(filter_kruskal_iterator(Graph(), threshold=0))
        Traceback (most recent call last):
        ...
        ValueError: the threshold mut be at least 1

    Check that a threshold of 1 is accepted::

        sage: len(list(filter_kruskal_iterator(graphs.HouseGraph(), threshold=1)))
        4

    Check that the method is robust to incomparable vertices::

        sage: G = Graph([(1, 2, 10), (1, 'a', 1), ('a', 'b', 1), ('b', 2, 1)])
        sage: E = list(filter_kruskal_iterator(G, by_weight=True))
        sage: sum(w for _, _, w in E)
        3
    """
    from sage.graphs.graph import Graph
    if not isinstance(G, Graph):
        raise ValueError("the input graph must be undirected")
    if threshold < 1:
        raise ValueError("the threshold mut be at least 1")
    if check:
        if not G.order() or not G.is_connected():
            return
        # G is now assumed to be a nonempty connected graph
        if G.order() == G.size() + 1:
            # G is a tree
            yield from G.edge_iterator()
            return

        g = G.to_simple(to_undirected=False, keep_label='min')
    else:
        g = G

    cdef int m = g.size()
    if m <= threshold:
        yield from kruskal_iterator_from_edges(g.edge_iterator(),
                                               DisjointSet_of_hashables(g),
                                               by_weight=by_weight,
                                               weight_function=weight_function,
                                               check_weight=check_weight)
        return

    #
    # Initialize some data structure
    #
    cdef list edges = list(g.edge_iterator())
    # Precompute edge weights to avoid frequent calls to weight_function
    cdef list weight
    _, weight_function = G._get_weight_function(by_weight=by_weight,
                                                weight_function=weight_function,
                                                check_weight=check_weight)
    if weight_function is None:
        weight = [1 for _ in range(m)]
    else:
        weight = [weight_function(e) for e in edges]

    cdef MemoryAllocator mem = MemoryAllocator()
    # Array storing a permutation of the edges.
    # e_index[i] is the position of edge i in list edges
    cdef int* e_index = <int*> mem.allocarray(m, sizeof(int))
    cdef int i, j
    for i in range(m):
        e_index[i] = i
    # Stack of range of edge partitions
    cdef list stack = [(0, m - 1)]
    cdef int begin, end
    # Parameter to  equally divide edges with weight equal the to pivot
    cdef bint ch = True
    # Data structure to record the vertices in each tree of the forest
    cdef DisjointSet_of_hashables union_find = DisjointSet_of_hashables(g)

    #
    # Iteratively partition the list of edges
    #
    while stack:
        begin, end = stack.pop()

        if end - begin < threshold:
            # Filter edges connecting vertices of a same tree
            L = [edges[e_index[i]] for i in range(begin, end + 1)
                 if union_find.find(edges[e_index[i]][0]) != union_find.find(edges[e_index[i]][1])]
            yield from kruskal_iterator_from_edges(L, union_find,
                                                   by_weight=by_weight,
                                                   weight_function=weight_function,
                                                   check_weight=False)
            if union_find.number_of_subsets() == 1:
                return
            continue

        # Choose a pivot
        pivot = weight[e_index[(begin + end) // 2]]

        # Partition edges with respect to pivot, as in quicksort
        i, j = begin, end
        while i < j:
            while weight[e_index[i]] < pivot and i < j:
                i += 1
            if ch and weight[e_index[i]] == pivot and i < j:
                i += 1
                ch = False
                continue
            while weight[e_index[j]] > pivot and i < j:
                j -= 1
            if not ch and weight[e_index[j]] == pivot and i < j:
                j -= 1
                ch = True
                continue
            if i < j:
                e_index[i], e_index[j] = e_index[j], e_index[i]

        # Record range of edge partitions
        if weight[e_index[i]] <= pivot:
            stack.append((i + 1, end))
            stack.append((begin, i))
        else:
            stack.append((i, end))
            stack.append((begin, i - 1))


def boruvka(G, by_weight=True, weight_function=None, check_weight=True, check=False):
    r"""
    Minimum spanning tree using Boruvka's algorithm.

    This function assumes that we can only compute minimum spanning trees for
    undirected graphs. Such graphs can be weighted or unweighted, and they can
    have multiple edges (since we are computing the minimum spanning tree, only
    the minimum weight among all `(u,v)`-edges is considered, for each pair of
    vertices `u`, `v`).

    INPUT:

    - ``G`` -- an undirected graph

    - ``by_weight`` -- boolean (default: ``True``); if ``True``, the edges in
      the graph are weighted. If ``False``, all edges have weight 1

    - ``weight_function`` -- function (default: ``None``); a function that takes
      as input an edge ``(u, v, l)`` and outputs its weight. If not ``None``,
      ``by_weight`` is automatically set to ``True``. If ``None`` and
      ``by_weight`` is ``True``, we use the edge label ``l``, if ``l`` is not
      ``None``, else ``1`` as a weight.

    - ``check_weight`` -- boolean (default: ``False``); whether to check that
      the ``weight_function`` outputs a number for each edge

    - ``check`` -- boolean (default: ``False``); whether to first perform sanity
      checks on the input graph ``G``. Default: ``check=False``. If we toggle
      ``check=True``, the following sanity checks are first performed on ``G``
      prior to running Boruvka's algorithm on that input graph:

      - Is ``G`` the null graph or graph on one vertex?
      - Is ``G`` disconnected?
      - Is ``G`` a tree?

      By default, we turn off the sanity checks for performance reasons. This
      means that by default the function assumes that its input graph is
      connected, and has at least one vertex. Otherwise, you should set
      ``check=True`` to perform some sanity checks and preprocessing on the
      input graph.

    OUTPUT:

    The edges of a minimum spanning tree of ``G``, if one exists, otherwise
    returns the empty list.

    .. SEEALSO::

        - :meth:`~sage.graphs.generic_graph.GenericGraph.min_spanning_tree`

    EXAMPLES:

    An example from pages 727--728 in [Sah2000]_::

        sage: from sage.graphs.spanning_tree import boruvka
        sage: G = Graph({1:{2:28, 6:10}, 2:{3:16, 7:14}, 3:{4:12}, 4:{5:22, 7:18}, 5:{6:25, 7:24}})
        sage: G.weighted(True)
        sage: E = boruvka(G, check=True); E
        [(1, 6, 10), (2, 7, 14), (3, 4, 12), (4, 5, 22), (5, 6, 25), (2, 3, 16)]
        sage: boruvka(G, by_weight=True)
        [(1, 6, 10), (2, 7, 14), (3, 4, 12), (4, 5, 22), (5, 6, 25), (2, 3, 16)]
        sage: sorted(boruvka(G, by_weight=False))
        [(1, 2, 28), (1, 6, 10), (2, 3, 16), (2, 7, 14), (3, 4, 12), (4, 5, 22)]

    An example with custom edge labels::

        sage: G = Graph([[0,1,1],[1,2,1],[2,0,10]], weighted=True)
        sage: weight = lambda e:3-e[0]-e[1]
        sage: boruvka(G, weight_function=lambda e:3-e[0]-e[1], by_weight=True)
        [(0, 2, 10), (1, 2, 1)]
        sage: boruvka(G, weight_function=lambda e:float(1/e[2]), by_weight=True)
        [(0, 2, 10), (0, 1, 1)]

    An example of disconnected graph with ``check`` disabled::

        sage: from sage.graphs.spanning_tree import boruvka
        sage: G = Graph({1:{2:28}, 3:{4:16}}, weighted=True)
        sage: boruvka(G, check=False)
        []

    TESTS:

    If the input graph is a tree, then return its edges::

        sage: T = graphs.RandomTree(randint(1, 10))
        sage: list(T.edges(sort=True)) == sorted(boruvka(T, check=True))
        True

    Check if the weight of MST returned by Prim's and Boruvka's is the same::

        sage: G = Graph([(u,v,randint(1,5)) for u,v in graphs.CompleteGraph(4).edges(sort=True, labels=0)], weighted=True)
        sage: G.weighted()
        True
        sage: E1 = G.min_spanning_tree(algorithm='Boruvka')
        sage: E2 = G.min_spanning_tree(algorithm='Prim_Boost')
        sage: sum(e[2] for e in E1) == sum(e[2] for e in E2)
        True

    If the input is not a Graph::

        sage: boruvka("I am not a graph")
        Traceback (most recent call last):
        ...
        ValueError: the input graph must be undirected
        sage: boruvka(digraphs.Path(10))
        Traceback (most recent call last):
        ...
        ValueError: the input graph must be undirected

    Check that the method is robust to incomparable vertices::

        sage: G = Graph([(1, 2, 10), (1, 'a', 1), ('a', 'b', 1), ('b', 2, 1)])
        sage: E = boruvka(G, by_weight=True)
        sage: sum(w for _, _, w in E)
        3
    """
    from sage.graphs.graph import Graph
    if not isinstance(G, Graph):
        raise ValueError("the input graph must be undirected")

    if G.order() <= 1:
        return []

    # sanity checks
    if check:
        if not G.is_connected():
            return []
        # G is now assumed to be a nonempty connected graph
        if G.num_verts() == G.num_edges() + 1:
            # G is a tree
            return G.edges(sort=False)

    by_weight, weight_function = G._get_weight_function(by_weight=by_weight,
                                                        weight_function=weight_function,
                                                        check_weight=check_weight)

    # Boruvka's algorithm

    # Store the list of active edges as (e, e_weight) in a list
    if weight_function is not None:
        edge_list = [(e, weight_function(e)) for e in G.edge_iterator()]
    else:
        edge_list = [(e, 1) for e in G.edge_iterator()]

    # initially, each vertex is a connected component
    cdef DisjointSet_of_hashables partitions = DisjointSet_of_hashables(G)
    # a dictionary to store the least weight outgoing edge for each component
    cdef dict cheapest = {}
    cdef list T = []  # stores the edges in minimum spanning tree
    cdef int numConComp = G.order()
    cdef int numConCompPrevIter = numConComp + 1

    # Dictionary to maintain active cheapest edges between pairs of components
    cdef dict components_dict = {}

    while numConComp > 1:
        # Check if number of connected components decreased.
        # Otherwise, the graph is not connected.
        if numConCompPrevIter == numConComp:
            return []
        else:
            numConCompPrevIter = numConComp

        # Iterate over all active edges to identify the cheapest edge between
        # each pair of components (trees of the forest), as well as cheapest
        # active edge incident to a component.
        for e, e_weight in edge_list:
            component1 = partitions.find(e[0])
            component2 = partitions.find(e[1])

            if component1 != component2:
                if component1 in cheapest:
                    if cheapest[component1][1] > e_weight:
                        cheapest[component1] = (e, e_weight)
                else:
                    cheapest[component1] = (e, e_weight)

                if component2 in cheapest:
                    if cheapest[component2][1] > e_weight:
                        cheapest[component2] = (e, e_weight)
                else:
                    cheapest[component2] = (e, e_weight)
                # store the cheapest edge between the two components
                pair = frozenset((component1, component2))
                if pair in components_dict:
                    if components_dict[pair][1] > e_weight:
                        components_dict[pair] = (e, e_weight)
                else:
                    components_dict[pair] = (e, e_weight)

        # Update the list of active edges
        edge_list = components_dict.values()

        # Go through all the current connected components and merge wherever
        # possible
        for v in cheapest:
            e, e_weight = cheapest[v]
            component1 = partitions.find(e[0])
            component2 = partitions.find(e[1])

            if component1 != component2:
                partitions.union(component1, component2)
                T.append(e)
                numConComp = numConComp - 1

        # reset the dictionaries for next iteration
        cheapest = {}
        components_dict = {}

    return T


def random_spanning_tree(G, output_as_graph=False, by_weight=False, weight_function=None, check_weight=True):
    r"""
    Return a random spanning tree of the graph.

    This uses the Aldous-Broder algorithm ([Bro1989]_, [Ald1990]_) to generate
    a random spanning tree with the uniform distribution, as follows.

    Start from any vertex. Perform a random walk by choosing at every step one
    neighbor uniformly at random. Every time a new vertex `j` is met, add the
    edge `(i, j)` to the spanning tree, where `i` is the previous vertex in the
    random walk.

    When ``by_weight`` is ``True`` or a weight function is given, the selection
    of the neighbor is done proportionaly to the edge weights.

    INPUT:

    - ``G`` -- an undirected graph

    - ``output_as_graph`` -- boolean (default: ``False``); whether to return a
      list of edges or a graph

    - ``by_weight`` -- boolean (default: ``False``); if ``True``, the edges in
      the graph are weighted, otherwise all edges have weight 1

    - ``weight_function`` -- function (default: ``None``); a function that takes
      as input an edge ``(u, v, l)`` and outputs its weight. If not ``None``,
      ``by_weight`` is automatically set to ``True``. If ``None`` and
      ``by_weight`` is ``True``, we use the edge label ``l`` , if ``l`` is not
      ``None``, else ``1`` as a weight. The ``weight_function`` can be used to
      transform the label into a weight (note that, if the weight returned is
      not convertible to a float, an error is raised)

    - ``check_weight`` -- boolean (default: ``True``); whether to check that
      the ``weight_function`` outputs a number for each edge

    .. SEEALSO::

        :meth:`~sage.graphs.generic_graph.GenericGraph.spanning_trees_count`
        and :meth:`~sage.graphs.graph.Graph.spanning_trees`

    EXAMPLES::

        sage: G = graphs.TietzeGraph()
        sage: G.random_spanning_tree(output_as_graph=True)
        Graph on 12 vertices
        sage: rg = G.random_spanning_tree(); rg # random
        [(0, 9),
        (9, 11),
        (0, 8),
        (8, 7),
        (7, 6),
        (7, 2),
        (2, 1),
        (1, 5),
        (9, 10),
        (5, 4),
        (2, 3)]
        sage: Graph(rg).is_tree()
        True

    A visual example for the grid graph::

        sage: G = graphs.Grid2dGraph(6, 6)
        sage: pos = G.get_pos()
        sage: T = G.random_spanning_tree(True)
        sage: T.set_pos(pos)
        sage: T.show(vertex_labels=False)                                               # needs sage.plot

    We can also use edge weights to change the probability of returning a
    spanning tree::

        sage: def foo(G, k):
        ....:     S = set()
        ....:     for _ in range(k):
        ....:         E = G.random_spanning_tree(by_weight=True)
        ....:         S.add(Graph(E).graph6_string())
        ....:     return S
        sage: K3 = graphs.CompleteGraph(3)
        sage: for u, v in K3.edges(sort=True, labels=False):
        ....:     K3.set_edge_label(u, v, randint(1, 2))
        sage: foo(K3, 100) == {'BW', 'Bg', 'Bo'}  # random
        True
        sage: K4 = graphs.CompleteGraph(4)
        sage: for u, v in K4.edges(sort=True, labels=False):
        ....:     K4.set_edge_label(u, v, randint(1, 2))
        sage: print(len(foo(K4, 100)))  # random
        16

    Check that the spanning tree returned when using weights is a tree::

        sage: # needs networkx
        sage: G = graphs.RandomBarabasiAlbert(50, 2)
        sage: for u, v in G.edge_iterator(labels=False):
        ....:     G.set_edge_label(u, v, randint(1, 10))
        sage: T = G.random_spanning_tree(by_weight=True, output_as_graph=True)
        sage: T.is_tree()
        True

    TESTS::

        sage: G = Graph()
        sage: G.random_spanning_tree()
        Traceback (most recent call last):
        ...
        ValueError: works only for non-empty connected graphs

        sage: G = graphs.CompleteGraph(3).complement()
        sage: G.random_spanning_tree()
        Traceback (most recent call last):
        ...
        ValueError: works only for non-empty connected graphs

    Check that the method is robust to incomparable vertices::

        sage: G = Graph([(1, 2, 10), (1, 'a', 1), ('a', 'b', 1), ('b', 2, 1)])
        sage: T = G.random_spanning_tree(by_weight=True, output_as_graph=True)
        sage: T.is_tree()
        True
    """
    from sage.misc.prandom import randint
    from sage.misc.prandom import random
    from sage.graphs.graph import Graph

    cdef int N = G.order()

    if not N or not G.is_connected():
        raise ValueError('works only for non-empty connected graphs')

    if G.order() == G.size() + 1:
        # G is a tree
        if output_as_graph:
            return G.copy()
        return list(G.edge_iterator(label=False))

    by_weight, weight_function = G._get_weight_function(by_weight=by_weight,
                                                        weight_function=weight_function,
                                                        check_weight=check_weight)

    if by_weight:
        def next_neighbor(s):
            p = random() * sum(weight_function(e)
                               for e in G.edge_iterator(s, sort_vertices=False))
            for e in G.edge_iterator(s, sort_vertices=False):
                p -= weight_function(e)
                if p <= 0:
                    break
            return e[1] if e[0] == s else e[0]
    else:
        def next_neighbor(s):
            return G.neighbors(s)[randint(0, G.degree(s) - 1)]

    s = next(G.vertex_iterator())
    cdef set found = set([s])
    cdef int found_nr = 1
    cdef list tree_edges = []
    while found_nr < N:
        new_s = next_neighbor(s)
        if new_s not in found:
            found.add(new_s)
            found_nr += 1
            tree_edges.append((s, new_s))
        s = new_s

    if not output_as_graph:
        return tree_edges
    return Graph(tree_edges)


def spanning_trees(g, labels=False):
    r"""
    Return an iterator over all spanning trees of the graph `g`.

    A disconnected graph has no spanning tree.

    Uses the Read-Tarjan backtracking algorithm [RT1975a]_.

    INPUT:

    - ``labels`` -- boolean (default: ``False``); whether to return edges labels
      in the spanning trees or not

    EXAMPLES::

        sage: G = Graph([(1,2),(1,2),(1,3),(1,3),(2,3),(1,4)], multiedges=True)
        sage: len(list(G.spanning_trees()))
        8
        sage: G.spanning_trees_count()                                                  # needs sage.modules
        8
        sage: G = Graph([(1,2),(2,3),(3,1),(3,4),(4,5),(4,5),(4,6)], multiedges=True)
        sage: len(list(G.spanning_trees()))
        6
        sage: G.spanning_trees_count()                                                  # needs sage.modules
        6

    .. SEEALSO::

        - :meth:`~sage.graphs.generic_graph.GenericGraph.spanning_trees_count`
          -- counts the number of spanning trees

        - :meth:`~sage.graphs.graph.Graph.random_spanning_tree`
          -- returns a random spanning tree

    TESTS:

    Works with looped graphs::

        sage: g = Graph({i:[i,(i+1)%6] for i in range(6)})
        sage: list(g.spanning_trees())
        [Graph on 6 vertices,
         Graph on 6 vertices,
         Graph on 6 vertices,
         Graph on 6 vertices,
         Graph on 6 vertices,
         Graph on 6 vertices]

    Edges of the spanning trees can be labeled or unlabeled (:issue:`27557`)::

        sage: g = Graph([(1,2,2),(1,2,1),(1,2,4),(1,4,5)],multiedges=True)
        sage: l = list(g.spanning_trees(labels=True))
        sage: l[0].edges(sort=True)
        [(1, 2, 4), (1, 4, 5)]
        sage: l[1].edges(sort=True)
        [(1, 2, 1), (1, 4, 5)]
        sage: l[2].edges(sort=True)
        [(1, 2, 2), (1, 4, 5)]

    Small cases::

        sage: list(Graph().spanning_trees())
        []
        sage: list(Graph(1).spanning_trees())
        [Graph on 1 vertex]
        sage: list(Graph(2).spanning_trees())
        []

    Giving anything else than a graph::

        sage: from sage.graphs.spanning_tree import spanning_trees
        sage: list(spanning_trees(DiGraph()))
        Traceback (most recent call last):
        ...
        ValueError: this method is for undirected graphs only
        sage: list(spanning_trees("bike"))
        Traceback (most recent call last):
        ...
        ValueError: this method is for undirected graphs only

    Check that the method is robust to incomparable vertices::

        sage: G = Graph([(1, 2, 10), (1, 'a', 1), ('a', 'b', 1), ('b', 2, 1)])
        sage: len(list(G.spanning_trees(labels=False)))
        4
    """
    from sage.graphs.graph import Graph
    if not isinstance(g, Graph):
        raise ValueError("this method is for undirected graphs only")

    def _recursive_spanning_trees(G, forest, labels):
        """
        Return an iterator over all the spanning trees of G containing forest
        """
        if not G.is_connected():
            return

        if G.size() == forest.size():
            yield forest.copy()
        else:
            # Pick an edge e from G-forest
            for e in G.edge_iterator(labels=labels):
                if not forest.has_edge(e):
                    break

            # 1) Recursive call with e removed from G
            G.delete_edge(e)
            yield from _recursive_spanning_trees(G, forest, labels)
            G.add_edge(e)

            # 2) Recursive call with e include in forest
            #
            # e=xy links the CC (connected component) of forest containing x
            # with the CC containing y. Any other edge which does that cannot be
            # added to forest anymore, and B is the list of them
            c1 = forest.connected_component_containing_vertex(e[0], sort=False)
            c2 = forest.connected_component_containing_vertex(e[1], sort=False)
            G.delete_edge(e)
            B = G.edge_boundary(c1, c2, sort=False)
            G.add_edge(e)

            # Actual call
            forest.add_edge(e)
            G.delete_edges(B)
            yield from _recursive_spanning_trees(G, forest, labels)
            G.add_edges(B)
            forest.delete_edge(e)

    if g.order() and g.is_connected():
        forest = Graph([g, g.bridges()], format='vertices_and_edges')
        yield from _recursive_spanning_trees(Graph(g, immutable=False, loops=False), forest, labels)


def edge_disjoint_spanning_trees(G, k, by_weight=False, weight_function=None, check_weight=True):
    r"""
    Return `k` edge-disjoint spanning trees of minimum cost.

    This method implements the Roskind-Tarjan algorithm for finding `k`
    minimum-cost edge-disjoint spanning trees in simple undirected graphs
    [RT1985]_. When edge weights are taken into account, the algorithm ensures
    that the sum of the weights of the returned spanning trees is minimized. The
    time complexity of the algorithm is in `O(k^2n^2)` for the unweighted case
    and otherwise in `O(m\log{m} + k^2n^2)`.

    This method raises an error if the graph does not contain the requested
    number of spanning trees.

    INPUT:

    - ``G`` -- a simple undirected graph

    - ``k`` -- the requested number of edge-disjoint spanning trees

    - ``by_weight`` -- boolean (default: ``False``); if ``True``, the edges in
      the graph are weighted, otherwise all edges have weight 1

    - ``weight_function`` -- function (default: ``None``); a function that takes
      as input an edge ``(u, v, l)`` and outputs its weight. If not ``None``,
      ``by_weight`` is automatically set to ``True``. If ``None`` and
      ``by_weight`` is ``True``, we use the edge label ``l``, if ``l`` is not
      ``None``, else ``1`` as a weight.

    - ``check_weight`` -- boolean (default: ``True``); if ``True``, we check
      that the ``weight_function`` outputs a number for each edge

    EXAMPLES:

    Example from [RT1985]_::

        sage: from sage.graphs.spanning_tree import edge_disjoint_spanning_trees
        sage: G = Graph({'a': ['b', 'c', 'd', 'e'], 'b': ['c', 'e'], 'c': ['d'], 'd': ['e']})
        sage: F = edge_disjoint_spanning_trees(G, 2)
        sage: F
        [Graph on 5 vertices, Graph on 5 vertices]
        sage: [f.is_tree() for f in F]
        [True, True]

    This method raises an error if the graph does not contain the required
    number of trees::

        sage: edge_disjoint_spanning_trees(G, 3)
        Traceback (most recent call last):
        ...
        EmptySetError: this graph does not contain the required number of trees/arborescences

    A clique of order `n` has `\lfloor n/2 \rfloor` edge disjoint spanning
    trees::

        sage: for n in range(1, 10):
        ....:     g = graphs.CompleteGraph(n)
        ....:     F = edge_disjoint_spanning_trees(g, n//2)

    The sum of the weights of the returned spanning trees is minimum::

        sage: g = graphs.CompleteGraph(5)
        sage: for u, v in g.edges(sort=True, labels=False):
        ....:     g.set_edge_label(u, v, 1)
        sage: g.set_edge_label(0, 1, 33)
        sage: g.set_edge_label(1, 3, 33)
        sage: F = edge_disjoint_spanning_trees(g, 2, by_weight=True)
        sage: sum(F[0].edge_labels()) + sum(F[1].edge_labels())
        8

    TESTS:

    A graph with a single vertex has a spanning tree::

        sage: from sage.graphs.spanning_tree import edge_disjoint_spanning_trees
        sage: edge_disjoint_spanning_trees(Graph(1), 1)
        [Graph on 1 vertex]

    Check parameter `k`::

        sage: G = graphs.CompleteGraph(4)
        sage: edge_disjoint_spanning_trees(G, -1)
        Traceback (most recent call last):
        ...
        ValueError: parameter k must be a nonnegative integer
        sage: edge_disjoint_spanning_trees(G, 0)
        []
        sage: edge_disjoint_spanning_trees(G, 1)
        [Graph on 4 vertices]

    This method is for undirected graphs only::

        sage: edge_disjoint_spanning_trees(DiGraph(), 1)
        Traceback (most recent call last):
        ...
        ValueError: this method is for undirected graphs only

    Check that the method is robust to incomparable vertices::

        sage: G = Graph()
        sage: G.add_clique([0, 1, 2, 'a', 'b'])
        sage: F = G.edge_disjoint_spanning_trees(k=2)
        sage: len(F)
        2
    """
    if G.is_directed():
        raise ValueError("this method is for undirected graphs only")
    G._scream_if_not_simple()

    from sage.categories.sets_cat import EmptySetError
    from sage.graphs.graph import Graph
    msg_no_solution = "this graph does not contain the required number of trees/arborescences"
    if k < 0:
        raise ValueError("parameter k must be a nonnegative integer")
    elif not k:
        return []
    elif k == 1:
        E = G.min_spanning_tree()
        if not E and G.order() != 1:
            raise EmptySetError(msg_no_solution)
        return [Graph([G, E], format='vertices_and_edges')]
    elif k > 1 + min(G.degree()) // 2:
        raise EmptySetError(msg_no_solution)

    # Initialization of data structures

    # - partition[0] is used to maintain known clumps.
    # - partition[i], 1 <= i <= k, is used to check if a given edge has both its
    #   endpoints in the same tree of forest Fi.
    partition = [DisjointSet_of_hashables(G) for _ in range(k + 1)]

    # Mapping from edge to forests:
    # - edge_index[e] == i if edge e is in Fi, and 0 if not in any Fi
    # This mapping is sufficient to extract the spanning trees.
    edge_index = {frozenset(e): 0 for e in G.edge_iterator(labels=False)}

    # Data structure to maintain the edge sets of each forest.
    # This is not a requirement of the algorithm as we can use the mapping
    # edge_index. However, it is convenient to maintain the forest as graphs to
    # simplify some operations.
    H = Graph([G, []], format='vertices_and_edges')
    F = [H.copy() for _ in range(k + 1)]

    # We consider the edges by increasing weight
    by_weight, weight_function = G._get_weight_function(by_weight=by_weight,
                                                        weight_function=weight_function,
                                                        check_weight=check_weight)
    if not by_weight:
        weight_function = None

    for x, y, _ in G.edges(sort=by_weight, key=weight_function):
        # {x, y} is edge e0 in the algorithm

        if partition[0].find(x) == partition[0].find(y):
            # x and y are in a same clump. That is x and y are in a same tree
            # in every forest Fi. We proceed with the next edge.
            continue

        # else, we apply the labeling algorithm

        # Label assigned to each edge by the labeling algorithm
        edge_label = {}

        # We use a queue of edges
        queue = [(x, y)]
        queue_begin = 0
        queue_end = 1

        # We find the tree Ti in Fi containing x, root Ti at x and
        # compute the parent pi(v) of every vertex in Ti
        p = [{x: x} for _ in range(k + 1)]
        for i in range(1, k + 1):
            # BFS will consider only vertices of the tree Ti of Fi containing x
            for u, v in F[i].breadth_first_search(x, edges=True):
                p[i][v] = u

        # and we search for an augmenting sequence
        augmenting_sequence_found = False
        while queue_begin < queue_end:
            e = queue[queue_begin]
            queue_begin += 1
            fe = frozenset(e)
            i = (edge_index[fe] % k) + 1
            v, w = e
            if partition[i].find(v) != partition[i].find(w):
                # v and w are in different subtrees of Fi. We have detected an
                # augmenting sequence since we can join the two subtrees.
                augmenting_sequence_found = True
                break
            else:
                # One of v and w is in the subtree of labeled edges in Fi
                if v == x or (v in p[i] and frozenset((v, p[i][v])) in edge_label):
                    u = w
                else:
                    u = v

                # Let F(e) be the unique path joining v and w.
                # We find the unlabeled edges of Fi(e) by ascending through the
                # tree one vertex at a time from z toward x, until reaching
                # either x or a previously labeled edge.

                # Stack of edges to be labeled
                edges_to_label = []
                while u != x and (u in p[i] and frozenset((u, p[i][u])) not in edge_label):
                    edges_to_label.append((u, p[i][u]))
                    u = p[i][u]

                # We now label edges
                while edges_to_label:
                    ep = edges_to_label.pop()
                    edge_label[frozenset(ep)] = fe
                    queue.append(ep)
                    queue_end += 1

        if augmenting_sequence_found:
            # We perform the corresponding augmentation
            partition[i].union(v, w)

            while fe in edge_label:
                F[edge_index[fe]].delete_edge(fe)
                F[i].add_edge(fe)
                e, edge_index[fe], i = edge_label[fe], i, edge_index[fe]
                fe = frozenset(e)

            # Finally, add edge e = e0 = (x, y) to Fi
            F[i].add_edge(e)
            edge_index[fe] = i

        else:
            # x and y are in a same tree in every Fi, so in a same clump
            partition[0].union(x, y)

    res = [F[i] for i in range(1, k + 1) if F[i].size() == G.order() - 1]
    if len(res) != k:
        raise EmptySetError(msg_no_solution)

    for f in res:
        for u, v in f.edges(sort=False, labels=False):
            f.set_edge_label(u, v, G.edge_label(u, v))
    return res
