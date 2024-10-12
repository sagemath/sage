r"""
Orientations

This module implements several methods to compute orientations of undirected
graphs subject to specific constraints (e.g., acyclic, strongly connected,
etc.). It also implements some iterators over all these orientations.

**This module contains the following methods**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`orient` | Return an oriented version of `G` according the input function `f`.
    :meth:`acyclic_orientations` | Return an iterator over all acyclic orientations of an undirected graph `G`.
    :meth:`strong_orientation` | Return a strongly connected orientation of the graph `G`.
    :meth:`strong_orientations_iterator` | Return an iterator over all strong orientations of a graph `G`
    :meth:`random_orientation` | Return a random orientation of a graph `G`


Authors
-------

- Kolja Knauer, Petru Valicov (2017-01-10) -- initial version


Methods
-------
"""
# ****************************************************************************
#       Copyright (C)      2017 Kolja Knauer <kolja.knauer@gmail.com>
#                          2017 Petru Valicov <petru.valicov@lirmm.fr>
#                     2017-2024 David Coudert <david.coudert@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from copy import copy
from sage.graphs.digraph import DiGraph


def orient(G, f, weighted=None, data_structure=None, sparse=None,
           immutable=None, hash_labels=None):
    r"""
    Return an oriented version of `G` according the input function `f`.

    INPUT:

    - ``G`` -- an undirected graph

    - ``f`` -- a function that inputs an edge and outputs an orientation of this
      edge

    - ``weighted`` -- boolean (default: ``None``); weightedness for the oriented
      digraph. By default (``None``), the graph and its orientation will behave
      the same.

    - ``sparse`` -- boolean (default: ``None``); ``sparse=True`` is an alias for
      ``data_structure="sparse"``, and ``sparse=False`` is an alias for
      ``data_structure="dense"``. Only used when ``data_structure=None``.

    - ``data_structure`` -- string (default: ``None``); one of ``'sparse'``,
      ``'static_sparse'``, or ``'dense'``. See the documentation of
      :class:`DiGraph`.

    - ``immutable`` -- boolean (default: ``None``); whether to create a
      mutable/immutable digraph. Only used when ``data_structure=None``.

      * ``immutable=None`` (default) means that the graph and its orientation
        will behave the same way.

      * ``immutable=True`` is a shortcut for ``data_structure='static_sparse'``

      * ``immutable=False`` means that the created digraph is mutable. When used
        to orient an immutable graph, the data structure used is ``'sparse'``
        unless anything else is specified.

    - ``hash_labels`` -- boolean (default: ``None``); whether to include edge
      labels during hashing of the oriented digraph. This parameter defaults to
      ``True`` if the graph is weighted. This parameter is ignored when
      parameter ``immutable`` is not ``True``. Beware that trying to hash
      unhashable labels will raise an error.

    OUTPUT: a :class:`DiGraph` object

    .. NOTE::

        This method behaves similarly to method
        :meth:`~sage.graphs.generic_graph.GenericGraph.copy`. That is, the
        returned digraph uses the same data structure by default, unless the
        user asks to use another data structure, and the attributes of the input
        graph are copied.

    EXAMPLES::

        sage: G = graphs.CycleGraph(4); G
        Cycle graph: Graph on 4 vertices
        sage: D = G.orient(lambda e:e if e[0] < e[1] else (e[1], e[0], e[2])); D
        Orientation of Cycle graph: Digraph on 4 vertices
        sage: sorted(D.edges(labels=False))
        [(0, 1), (0, 3), (1, 2), (2, 3)]

    TESTS:

    We make sure that one can get an immutable orientation by providing the
    ``data_structure`` optional argument::

        sage: def foo(e):
        ....:     return e if e[0] < e[1] else (e[1], e[0], e[2])
        sage: G = graphs.CycleGraph(4)
        sage: D = G.orient(foo, data_structure='static_sparse')
        sage: D.is_immutable()
        True
        sage: D = G.orient(foo, immutable=True)
        sage: D.is_immutable()
        True

    Bad input::

        sage: G.orient(foo, data_structure='sparse', sparse=False)
        Traceback (most recent call last):
        ...
        ValueError: you cannot define 'immutable' or 'sparse' when 'data_structure' has a value
        sage: G.orient(foo, data_structure='sparse', immutable=True)
        Traceback (most recent call last):
        ...
        ValueError: you cannot define 'immutable' or 'sparse' when 'data_structure' has a value
        sage: G.orient(foo, immutable=True, sparse=False)
        Traceback (most recent call last):
        ...
        ValueError: there is no dense immutable backend at the moment

    Which backend? ::

        sage: G.orient(foo, data_structure='sparse')._backend
        <sage.graphs.base.sparse_graph.SparseGraphBackend object at ...>
        sage: G.orient(foo, data_structure='dense')._backend
        <sage.graphs.base.dense_graph.DenseGraphBackend object at ...>
        sage: G.orient(foo, data_structure='static_sparse')._backend
        <sage.graphs.base.static_sparse_backend.StaticSparseBackend object at ...>
        sage: G.orient(foo, immutable=True)._backend
        <sage.graphs.base.static_sparse_backend.StaticSparseBackend object at ...>
        sage: G.orient(foo, immutable=True, sparse=True)._backend
        <sage.graphs.base.static_sparse_backend.StaticSparseBackend object at ...>
        sage: G.orient(foo, immutable=False, sparse=True)._backend
        <sage.graphs.base.sparse_graph.SparseGraphBackend object at ...>
        sage: G.orient(foo, immutable=False, sparse=False)._backend
        <sage.graphs.base.sparse_graph.SparseGraphBackend object at ...>
        sage: G.orient(foo, data_structure=None, immutable=None, sparse=True)._backend
        <sage.graphs.base.sparse_graph.SparseGraphBackend object at ...>
        sage: G.orient(foo, data_structure=None, immutable=None, sparse=False)._backend
        <sage.graphs.base.dense_graph.DenseGraphBackend object at ...>
        sage: G.orient(foo, data_structure=None, immutable=None, sparse=None)._backend
        <sage.graphs.base.sparse_graph.SparseGraphBackend object at ...>
        sage: H = Graph(data_structure='dense')
        sage: H.orient(foo, data_structure=None, immutable=None, sparse=None)._backend
        <sage.graphs.base.dense_graph.DenseGraphBackend object at ...>
    """
    # Which data structure should be used ?
    if data_structure is not None:
        # data_structure is already defined so there is nothing left to do
        # here. Did the user try to define too much ?
        if immutable is not None or sparse is not None:
            raise ValueError("you cannot define 'immutable' or 'sparse' "
                             "when 'data_structure' has a value")
    # At this point, data_structure is None.
    elif immutable is True:
        data_structure = 'static_sparse'
        if sparse is False:
            raise ValueError("there is no dense immutable backend at the moment")
    elif immutable is False:
        # If the user requests a mutable digraph and input is immutable, we
        # choose the 'sparse' cgraph backend. Unless the user explicitly
        # asked for something different.
        if G.is_immutable():
            data_structure = 'dense' if sparse is False else 'sparse'
    elif sparse is True:
        data_structure = "sparse"
    elif sparse is False:
        data_structure = "dense"

    if data_structure is None:
        from sage.graphs.base.dense_graph import DenseGraphBackend
        if isinstance(G._backend, DenseGraphBackend):
            data_structure = "dense"
        else:
            data_structure = "sparse"

    if weighted is None:
        weighted = G.weighted()

    edges = (f(e) for e in G.edge_iterator())
    D = DiGraph([G, edges], format='vertices_and_edges',
                data_structure=data_structure,
                loops=G.allows_loops(),
                multiedges=G.allows_multiple_edges(),
                name=f"Orientation of {G.name()}",
                pos=copy(G._pos), weighted=weighted,
                hash_labels=hash_labels)

    attributes_to_copy = ('_assoc', '_embedding')
    for attr in attributes_to_copy:
        if hasattr(G, attr):
            copy_attr = {}
            old_attr = getattr(G, attr)
            if isinstance(old_attr, dict):
                for v, value in old_attr.items():
                    try:
                        copy_attr[v] = value.copy()
                    except AttributeError:
                        copy_attr[v] = copy(value)
                setattr(D, attr, copy_attr)
            else:
                setattr(D, attr, copy(old_attr))

    return D


def acyclic_orientations(G):
    r"""
    Return an iterator over all acyclic orientations of an undirected graph `G`.

    ALGORITHM:

    The algorithm is based on [Sq1998]_.
    It presents an efficient algorithm for listing the acyclic orientations of a
    graph. The algorithm is shown to require O(n) time per acyclic orientation
    generated, making it the most efficient known algorithm for generating acyclic
    orientations.

    The function uses a recursive approach to generate acyclic orientations of the
    graph. It reorders the vertices and edges of the graph, creating a new graph
    with updated labels. Then, it iteratively generates acyclic orientations by
    considering subsets of edges and checking whether they form upsets in a
    corresponding poset.

    INPUT:

    - ``G`` -- an undirected graph

    OUTPUT: an iterator over all acyclic orientations of the input graph

    .. NOTE::

        The function assumes that the input graph is undirected and the edges
        are unlabelled.

    EXAMPLES:

    To count the number of acyclic orientations for a graph::

        sage: g = Graph([(0, 3), (0, 4), (3, 4), (1, 3), (1, 2), (2, 3), (2, 4)])
        sage: it = g.acyclic_orientations()
        sage: len(list(it))
        54

    Test for arbitrary vertex labels::

        sage: g_str = Graph([('abc', 'def'), ('ghi', 'def'), ('xyz', 'abc'),
        ....:                ('xyz', 'uvw'), ('uvw', 'abc'), ('uvw', 'ghi')])
        sage: it = g_str.acyclic_orientations()
        sage: len(list(it))
        42

    Check that the method returns properly relabeled acyclic digraphs::

        sage: g = Graph([(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)])
        sage: orientations = set([frozenset(d.edges(labels=false)) for d in g.acyclic_orientations()])
        sage: len(orientations)
        18
        sage: all(d.is_directed_acyclic() for d in g.acyclic_orientations())
        True

    TESTS:

    To count the number of acyclic orientations for a graph with 0 vertices::

        sage: list(Graph().acyclic_orientations())
        []

    To count the number of acyclic orientations for a graph with 1 vertex::

        sage: list(Graph(1).acyclic_orientations())
        []

    To count the number of acyclic orientations for a graph with 2 vertices::

        sage: list(Graph(2).acyclic_orientations())
        []

    Acyclic orientations of a complete graph::

        sage: g = graphs.CompleteGraph(5)
        sage: it = g.acyclic_orientations()
        sage: len(list(it))
        120

    Graph with one edge::

        sage: list(Graph([(0, 1)]).acyclic_orientations())
        [Digraph on 2 vertices, Digraph on 2 vertices]

    Graph with two edges::

        sage: len(list(Graph([(0, 1), (1, 2)]).acyclic_orientations()))
        4

    Cycle graph::

        sage: len(list(Graph([(0, 1), (1, 2), (2, 0)]).acyclic_orientations()))
        6
    """
    if not G.size():
        # A graph without edge cannot be oriented
        return

    from sage.rings.infinity import Infinity
    from sage.combinat.subset import Subsets

    def reorder_vertices(G):
        n = G.order()
        ko = n
        k = n
        G_copy = G.copy()
        vertex_labels = {v: None for v in G_copy.vertices()}

        while G_copy.size() > 0:
            min_val = float('inf')
            uv = None
            for u, v, _ in G_copy.edges():
                du = G_copy.degree(u)
                dv = G_copy.degree(v)
                val = (du + dv) / (du * dv)
                if val < min_val:
                    min_val = val
                    uv = (u, v)

            if uv:
                u, v = uv
                vertex_labels[u] = ko
                vertex_labels[v] = ko - 1
                G_copy.delete_vertex(u)
                G_copy.delete_vertex(v)
                ko -= 2

            if G_copy.size() == 0:
                break

        for vertex, label in vertex_labels.items():
            if label is None:
                vertex_labels[vertex] = ko
                ko -= 1

        return vertex_labels

    def order_edges(G, vertex_labels):
        n = len(vertex_labels)
        m = 1
        edge_labels = {}

        for j in range(2, n + 1):
            for i in range(1, j):
                if G.has_edge(i, j):
                    edge_labels[(i, j)] = m
                    m += 1

        return edge_labels

    def is_upset_of_poset(Poset, subset, keys):
        for (u, v) in subset:
            for (w, x) in keys:
                if (Poset[(u, v), (w, x)] == 1 and (w, x) not in subset):
                    return False
        return True

    def generate_orientations(globO, starting_of_Ek, m, k, keys):
        # Creating a poset
        Poset = {}
        for i in range(starting_of_Ek, m - 1):
            for j in range(starting_of_Ek, m - 1):
                u, v = keys[i]
                w, x = keys[j]
                Poset[(u, v), (w, x)] = 0

        # Create a new graph to determine reachable vertices
        new_G = DiGraph()

        # Process vertices up to starting_of_Ek
        new_G.add_edges([(v, u) if globO[(u, v)] == 1 else (u, v) for u, v in keys[:starting_of_Ek]])

        # Process vertices starting from starting_of_Ek
        new_G.add_vertices([u for u, _ in keys[starting_of_Ek:]] + [v for _, v in keys[starting_of_Ek:]])

        if (globO[(k-1, k)] == 1):
            new_G.add_edge(k, k - 1)
        else:
            new_G.add_edge(k-1, k)

        for i in range(starting_of_Ek, m - 1):
            for j in range(starting_of_Ek, m - 1):
                u, v = keys[i]
                w, x = keys[j]
                # w should be reachable from u and v should be reachable from x
                if w in new_G.depth_first_search(u) and v in new_G.depth_first_search(x):
                    Poset[(u, v), (w, x)] = 1

        # For each subset of the base set of E_k, check if it is an upset or not
        upsets = []
        for subset in Subsets(keys[starting_of_Ek:m-1]):
            if (is_upset_of_poset(Poset, subset, keys[starting_of_Ek:m-1])):
                upsets.append(list(subset))

        for upset in upsets:
            for i in range(starting_of_Ek, m - 1):
                u, v = keys[i]
                if (u, v) in upset:
                    globO[(u, v)] = 1
                else:
                    globO[(u, v)] = 0

            yield globO.copy()

    def helper(G, globO, m, k):
        keys = list(globO.keys())
        keys = keys[0:m]

        if m <= 0:
            yield {}
            return

        starting_of_Ek = 0
        for (u, v) in keys:
            if u >= k - 1 or v >= k - 1:
                break
            else:
                starting_of_Ek += 1

        # s is the size of E_k
        s = m - 1 - starting_of_Ek

        # Recursively generate acyclic orientations
        orientations_G_small = helper(G, globO, starting_of_Ek, k - 2)

        # For each orientation of G_k-2, yield acyclic orientations
        for alpha in orientations_G_small:
            for (u, v) in alpha:
                globO[(u, v)] = alpha[(u, v)]

            # Orienting H_k as 1
            globO[(k-1, k)] = 1
            yield from generate_orientations(globO, starting_of_Ek, m, k, keys)

            # Orienting H_k as 0
            globO[(k-1, k)] = 0
            yield from generate_orientations(globO, starting_of_Ek, m, k, keys)

    # Reorder vertices based on the logic in reorder_vertices function
    vertex_labels = reorder_vertices(G)

    # Create a new graph with updated vertex labels using SageMath, Assuming the graph edges are unlabelled
    new_G = G.relabel(perm=vertex_labels, inplace=False)

    G = new_G

    # Order the edges based on the logic in order_edges function
    edge_labels = order_edges(G, vertex_labels)

    # Create globO array
    globO = {uv: 0 for uv in edge_labels}

    m = len(edge_labels)
    k = len(vertex_labels)
    orientations = helper(G, globO, m, k)

    # Create a mapping between original and new vertex labels
    reverse_vertex_labels = {label: vertex for vertex, label in vertex_labels.items()}

    # Iterate over acyclic orientations and create relabeled graphs
    for orientation in orientations:
        D = DiGraph([(u, v) if label else (v, u) for (u, v), label in orientation.items()])
        D.relabel(perm=reverse_vertex_labels, inplace=True)
        yield D


def strong_orientation(G):
    r"""
    Return a strongly connected orientation of the graph `G`.

    An orientation of an undirected graph is a digraph obtained by giving an
    unique direction to each of its edges. An orientation is said to be strong
    if there is a directed path between each pair of vertices. See also the
    :wikipedia:`Strongly_connected_component`.

    If the graph is 2-edge-connected, a strongly connected orientation can be
    found in linear time. If the given graph is not 2-connected, the orientation
    returned will ensure that each 2-connected component has a strongly
    connected orientation.

    OUTPUT: a digraph representing an orientation of the current graph

    .. NOTE::

        - This method assumes that the input the graph is connected.
        - The time complexity is `O(n+m)` for ``SparseGraph`` and `O(n^2)` for
          ``DenseGraph`` .

    .. SEEALSO::

        - :meth:`~sage.graphs.graph.Graph.orientations`
        - :meth:`~sage.graphs.orientations.strong_orientations_iterator`
        - :meth:`~sage.graphs.digraph_generators.DiGraphGenerators.nauty_directg`
        - :meth:`~sage.graphs.orientations.random_orientation`

    EXAMPLES:

    For a 2-regular graph, a strong orientation gives to each vertex an
    out-degree equal to 1::

        sage: g = graphs.CycleGraph(5)
        sage: g.strong_orientation().out_degree()
        [1, 1, 1, 1, 1]

    The Petersen Graph is 2-edge connected. It then has a strongly connected
    orientation::

        sage: g = graphs.PetersenGraph()
        sage: o = g.strong_orientation()
        sage: len(o.strongly_connected_components())
        1

    The same goes for the CubeGraph in any dimension::

        sage: all(len(graphs.CubeGraph(i).strong_orientation().strongly_connected_components()) == 1
        ....:     for i in range(2,6))
        True

    A multigraph also has a strong orientation::

        sage: g = Graph([(1,2),(1,2)], multiedges=True)
        sage: g.strong_orientation()
        Multi-digraph on 2 vertices
    """
    d = DiGraph(multiedges=G.allows_multiple_edges())
    i = 0

    # The algorithm works through a depth-first search. Any edge used in the
    # depth-first search is oriented in the direction in which it has been
    # used. All the other edges are oriented backward

    v = next(G.vertex_iterator())
    seen = {}
    i = 1

    # Time at which the vertices have been discovered
    seen[v] = i

    # indicates the stack of edges to explore
    next_ = G.edges_incident(v)

    while next_:
        e = next_.pop()

        # Ignore loops
        if e[0] == e[1]:
            continue

        # We assume e[0] to be a `seen` vertex
        e = e if seen.get(e[0], False) is not False else (e[1], e[0], e[2])

        # If we discovered a new vertex
        if seen.get(e[1], False) is False:
            d.add_edge(e)
            next_.extend(ee for ee in G.edges_incident(e[1])
                         if ((e[0], e[1]) != (ee[0], ee[1])) and ((e[0], e[1]) != (ee[1], ee[0])))
            i += 1
            seen[e[1]] = i

        # Else, we orient the edges backward
        else:
            if seen[e[0]] < seen[e[1]]:
                d.add_edge(e[1], e[0], e[2])
            else:
                d.add_edge(e)

    # Case of multiple edges. If another edge has already been inserted, we add
    # the new one in the opposite direction.
    tmp = None
    for e in G.multiple_edges():
        if tmp == (e[0], e[1]):
            if d.has_edge(e[0], e[1]):
                d.add_edge(e[1], e[0], e[2])
            else:
                d.add_edge(e)
        tmp = (e[0], e[1])

    return d


def strong_orientations_iterator(G):
    r"""
    Return an iterator over all strong orientations of a graph `G`.

    A strong orientation of a graph is an orientation of its edges such that the
    obtained digraph is strongly connected (i.e. there exist a directed path
    between each pair of vertices). According to Robbins' theorem (see the
    :wikipedia:`Robbins_theorem`), the graphs that have strong orientations are
    exactly the 2-edge-connected graphs (i.e., the bridgeless graphs).

    ALGORITHM:

    It is an adaptation of the algorithm published in [CGMRV16]_.
    It runs in `O(mn)` amortized time, where `m` is the number of edges and
    `n` is the number of vertices. The amortized time can be improved to `O(m)`
    with a more involved method.
    In this function, first the graph is preprocessed and a spanning tree is
    generated. Then every orientation of the non-tree edges of the graph can be
    extended to at least one new strong orientation by orienting properly
    the edges of the spanning tree (this property is proved in [CGMRV16]_).
    Therefore, this function generates all partial orientations of the non-tree
    edges and then launches a helper function corresponding to the generation
    algorithm described in [CGMRV16]_.
    In order to avoid trivial symmetries, the orientation of an arbitrary edge
    is fixed before the start of the enumeration process.

    INPUT:

    - ``G`` -- an undirected graph

    OUTPUT: an iterator which will produce all strong orientations of this graph

    .. NOTE::

        Works only for simple graphs (no multiple edges).
        To avoid symmetries an orientation of an arbitrary edge is fixed.

    .. SEEALSO::

        - :meth:`~Graph.orientations`
        - :meth:`~Graph.strong_orientation`
        - :meth:`~sage.graphs.digraph_generators.DiGraphGenerators.nauty_directg`
        - :meth:`~sage.graphs.orientations.random_orientation`

    EXAMPLES:

    A cycle has one possible (non-symmetric) strong orientation::

        sage: g = graphs.CycleGraph(4)
        sage: it = g.strong_orientations_iterator()
        sage: len(list(it))
        1

    A tree cannot be strongly oriented::

        sage: g = graphs.RandomTree(10)
        sage: len(list(g.strong_orientations_iterator()))
        0

    Neither can be a disconnected graph::

        sage: g = graphs.CompleteGraph(6)
        sage: g.add_vertex(7)
        sage: len(list(g.strong_orientations_iterator()))
        0

    TESTS:

        sage: g = graphs.CompleteGraph(2)
        sage: len(list(g.strong_orientations_iterator()))
        0

        sage: g = graphs.CubeGraph(3)
        sage: b = True
        sage: for orientedGraph in g.strong_orientations_iterator():
        ....:     if not orientedGraph.is_strongly_connected():
        ....:         b = False
        sage: b
        True

    The total number of strong orientations of a graph can be counted using
    the Tutte polynomial evaluated at points (0,2)::

        sage: g = graphs.PetersenGraph()
        sage: nr1 = len(list(g.strong_orientations_iterator()))
        sage: nr2 = g.tutte_polynomial()(0,2)
        sage: nr1 == nr2/2  # The Tutte polynomial counts also the symmetrical orientations
        True
    """
    # if the graph has a bridge or is disconnected,
    # then it cannot be strongly oriented
    if G.order() < 3 or not G.is_connected() or any(G.bridges(labels=False)):
        return

    V = list(G)

    # compute an arbitrary spanning tree of the undirected graph
    T = G.subgraph(vertices=G, edges=G.min_spanning_tree(), inplace=False)
    treeEdges = list(T.edges(labels=False, sort=False))
    A = [edge for edge in G.edge_iterator(labels=False) if not T.has_edge(edge)]

    # Initialize a digraph with the edges of the spanning tree doubly oriented
    Dg = T.to_directed(sparse=True)
    Dg.add_edges(A)

    # initialization of the first binary word 00...0
    # corresponding to the current orientation of the non-tree edges
    existingAedges = [0] * len(A)

    # Generate all orientations for non-tree edges (using Gray code)
    # Each of these orientations can be extended to a strong orientation
    # of G by orienting properly the tree-edges
    previousWord = 0

    # the orientation of one edge is fixed so we consider one edge less
    nr = 2**(len(A) - 1)
    for i in range(nr):
        word = (i >> 1) ^ i
        bitChanged = word ^ previousWord

        bit = 0
        while bitChanged > 1:
            bitChanged >>= 1
            bit += 1

        previousWord = word
        if not existingAedges[bit]:
            Dg.reverse_edge(A[bit])
            existingAedges[bit] = 1
        else:
            Dg.reverse_edge(A[bit][1], A[bit][0])
            existingAedges[bit] = 0
        # launch the algorithm for enumeration of the solutions
        yield from _strong_orientations_of_a_mixed_graph(Dg, V, treeEdges)


def _strong_orientations_of_a_mixed_graph(Dg, V, E):
    r"""
    Helper function for the generation of all strong orientations.

    Generates all strong orientations of a given partially directed graph
    (also called mixed graph). The algorithm finds bound edges i.e undirected
    edges whose orientation is forced and tries all possible orientations for
    the other edges. See [CGMRV16]_ for more details.

    INPUT:

    - ``Dg`` -- the mixed graph. The undirected edges are doubly oriented

    - ``V`` -- the set of vertices

    - ``E`` -- the set of undirected edges (they are oriented in both ways);
      no labels are allowed

    OUTPUT: an iterator which will produce all strong orientations of the input
    partially directed graph

    EXAMPLES::

        sage: from sage.graphs.orientations import _strong_orientations_of_a_mixed_graph
        sage: g = graphs.CycleGraph(5)
        sage: Dg = DiGraph(g)  # all edges of g will be doubly oriented
        sage: it = _strong_orientations_of_a_mixed_graph(Dg, list(g), list(g.edges(labels=False, sort=False)))
        sage: len(list(it))  # there are two orientations of this multigraph
        2
    """
    length = len(E)
    i = 0
    boundEdges = []
    while i < length:
        u, v = E[i]
        Dg.delete_edge(u, v)
        if v not in Dg.depth_first_search(u):
            # del E[i] in constant time
            E[i] = E[-1]
            E.pop()
            length -= 1
            Dg.add_edge(u, v)
            Dg.delete_edge(v, u)
            boundEdges.append((v, u))
        else:
            Dg.add_edge(u, v)
            Dg.delete_edge(v, u)
            if u not in Dg.depth_first_search(v):
                # del E[i] in constant time
                E[i] = E[-1]
                E.pop()
                length -= 1
                boundEdges.append((u, v))
                Dg.delete_edge(u, v)
            else:
                i += 1
            Dg.add_edge(v, u)

    # if true the obtained orientation is strong
    if not E:
        yield Dg.copy()
    else:
        u, v = E.pop()
        Dg.delete_edge(v, u)
        for orientation in _strong_orientations_of_a_mixed_graph(Dg, V, E):
            yield orientation
        Dg.add_edge(v, u)
        Dg.delete_edge(u, v)
        for orientation in _strong_orientations_of_a_mixed_graph(Dg, V, E):
            yield orientation
        Dg.add_edge(u, v)
        E.append((u, v))
    Dg.add_edges(boundEdges)
    E.extend(boundEdges)


def random_orientation(G):
    r"""
    Return a random orientation of a graph `G`.

    An *orientation* of an undirected graph is a directed graph such that every
    edge is assigned a direction. Hence there are `2^m` oriented digraphs for a
    simple graph with `m` edges.

    INPUT:

    - ``G`` -- a Graph

    EXAMPLES::

        sage: from sage.graphs.orientations import random_orientation
        sage: G = graphs.PetersenGraph()
        sage: D = random_orientation(G)
        sage: D.order() == G.order(), D.size() == G.size()
        (True, True)

    TESTS:

    Giving anything else than a Graph::

        sage: random_orientation([])
        Traceback (most recent call last):
        ...
        ValueError: the input parameter must be a Graph

    .. SEEALSO::

        - :meth:`~Graph.orientations`
        - :meth:`~Graph.strong_orientation`
        - :meth:`~sage.graphs.orientations.strong_orientations_iterator`
        - :meth:`~sage.graphs.digraph_generators.DiGraphGenerators.nauty_directg`
    """
    from sage.graphs.graph import Graph
    if not isinstance(G, Graph):
        raise ValueError("the input parameter must be a Graph")

    D = DiGraph(data=[G.vertices(sort=False), []],
                format='vertices_and_edges',
                multiedges=G.allows_multiple_edges(),
                loops=G.allows_loops(),
                weighted=G.weighted(),
                pos=G.get_pos(),
                name="Random orientation of {}".format(G.name()))
    if hasattr(G, '_embedding'):
        D._embedding = copy(G._embedding)

    from sage.misc.prandom import getrandbits
    rbits = getrandbits(G.size())
    for u, v, l in G.edge_iterator():
        if rbits % 2:
            D.add_edge(u, v, l)
        else:
            D.add_edge(v, u, l)
        rbits >>= 1
    return D
