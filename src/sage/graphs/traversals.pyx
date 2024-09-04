# cython: binding=True
# distutils: language = c++
r"""
Graph traversals

**This module implements the following graph traversals**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~lex_BFS` | Perform a lexicographic breadth first search (LexBFS) on the graph.
    :meth:`~lex_DFS` | Perform a lexicographic depth first search (LexDFS) on the graph.
    :meth:`~lex_UP` | Perform a lexicographic UP search (LexUP) on the graph.
    :meth:`~lex_DOWN` | Perform a lexicographic DOWN search (LexDOWN) on the graph.
    :meth:`~lex_M` |     Return an ordering of the vertices according the LexM graph traversal.
    :meth:`~lex_M_slow` | Return an ordering of the vertices according the LexM graph traversal.
    :meth:`~lex_M_fast` | Return an ordering of the vertices according the LexM graph traversal.
    :meth:`~maximum_cardinality_search` | Return an ordering of the vertices according a maximum cardinality search.
    :meth:`~maximum_cardinality_search_M` | Return the ordering and the edges of the triangulation produced by MCS-M.


ALGORITHM:

For :meth:`~lex_BFS` with ``algorithm="slow"``, :meth:`~lex_DFS`,
:meth:`~lex_UP` and :meth:`~lex_DOWN` the same generic implementation is used.
It corresponds to an implementation the generic algorithm described in
"Algorithm 1" of [Mil2017]_.

This algorithm maintains for each vertex left in the graph a lexicographic label
corresponding to the vertices already removed. The vertex of maximal
lexicographic label is then removed, and the lexicographic labels of its
neighbors are updated. Depending on how the update is done, it corresponds to
LexBFS, LexUP, LexDFS or LexDOWN: during the `i`-th iteration of the algorithm
`n-i` (for LexBFS and LexDOWN) or `i` (for LexDFS and LexUP) is appended (for
LexBFS and LexUP) or prepended (for LexDFS and LexDOWN) to the lexicographic
labels of all neighbors of the selected vertex that are left in the graph.

The time complexity of the algorithm is `O(mn)` for ``SparseGraph`` and
`O(\max\{mn, n^2\})` for ``DenseGraph``, where `n` is the number of vertices
and `m` is the number of edges.

See [CK2008]_ and [Mil2017]_ for more details on the algorithm and graphs
searching.

Methods
-------
"""
# ****************************************************************************
# Copyright (C) 2019 Georgios Giapitzakis Tzintanos <giorgosgiapis@mail.com>
#                    David Coudert <david.coudert@inria.fr>
#               2024 Cyril Bouvier <cyril.bouvier@lirmm.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from collections import deque

from libc.string cimport memset
from libc.stdint cimport uint32_t
from libcpp.queue cimport priority_queue
from libcpp.pair cimport pair
from libcpp.vector cimport vector
from cysignals.signals cimport sig_on, sig_off
from memory_allocator cimport MemoryAllocator

from sage.graphs.base.static_sparse_graph cimport init_short_digraph
from sage.graphs.base.static_sparse_graph cimport free_short_digraph
from sage.graphs.base.static_sparse_graph cimport out_degree
from sage.graphs.base.c_graph cimport CGraph, CGraphBackend
from sage.graphs.graph_decompositions.slice_decomposition cimport \
        extended_lex_BFS


def _lex_order_common(G, algo, reverse, tree, initial_vertex):
    r"""
    Perform a lexicographic search (LexBFS, LexUP, LexDFS or LexDOWN) on the
    graph.

    INPUT:

    - ``G`` -- a sage graph

    - ``algo`` -- string; the name of the actual algorithm among:

        - ``"lex_BFS"``

        - ``"lex_UP"``

        - ``"lex_DFS"``

        - ``"lex_DOWN"``

    - ``reverse`` -- whether to return the vertices in discovery order, or the
      reverse

    - ``tree`` -- whether to return the discovery directed tree (each vertex
      being linked to the one that saw it last)

    - ``initial_vertex`` -- the first vertex to consider, can be None

    .. NOTE::

        Loops and multiple edges are ignored and directed graphs are considered
        as undirected graphs.

    ALGORITHM:

    See the documentation of the :mod:`~sage.graphs.traversals` module.

    TESTS:

    Lex ordering of a graph on one vertex::

        sage: Graph(1).lex_BFS(tree=True, algorithm="slow")
        ([0], Digraph on 1 vertex)
        sage: Graph(1).lex_UP(tree=True)
        ([0], Digraph on 1 vertex)
        sage: Graph(1).lex_DFS(tree=True)
        ([0], Digraph on 1 vertex)
        sage: Graph(1).lex_DOWN(tree=True)
        ([0], Digraph on 1 vertex)

    Lex ordering of an empty (di)graph is an empty sequence::

        sage: g = Graph()
        sage: g.lex_BFS(algorithm="slow")
        []
        sage: g.lex_BFS(algorithm="slow", tree=True)
        ([], Digraph on 0 vertices)
        sage: g.lex_UP()
        []
        sage: g.lex_UP(tree=True)
        ([], Digraph on 0 vertices)
        sage: g.lex_DFS()
        []
        sage: g.lex_DFS(tree=True)
        ([], Digraph on 0 vertices)
        sage: g.lex_DOWN()
        []
        sage: g.lex_DFS(tree=True)
        ([], Digraph on 0 vertices)

    Lex UP ordering of a symmetric digraph should be the same as the Lex UP
    ordering of the corresponding undirected graph::

        sage: G = Graph([(1, 2), (1, 3), (2, 3), (2, 4), (2, 5), (3, 5), (3, 6), (4, 5), (5, 6)])
        sage: H = DiGraph(G)
        sage: G.lex_BFS(algorithm="slow") == H.lex_BFS(algorithm="slow")
        True
        sage: G.lex_UP() == H.lex_UP()
        True
        sage: G.lex_DFS() == H.lex_DFS()
        True
        sage: G.lex_DOWN() == H.lex_DOWN()
        True

    ``initial_vertex`` should be a valid graph vertex::

        sage: G = graphs.CompleteGraph(6)
        sage: G.lex_BFS(initial_vertex='foo', algorithm="slow")
        Traceback (most recent call last):
        ...
        ValueError: 'foo' is not a graph vertex
        sage: G.lex_UP(initial_vertex='foo')
        Traceback (most recent call last):
        ...
        ValueError: 'foo' is not a graph vertex
        sage: G.lex_DFS(initial_vertex='foo')
        Traceback (most recent call last):
        ...
        ValueError: 'foo' is not a graph vertex
        sage: G.lex_DOWN(initial_vertex='foo')
        Traceback (most recent call last):
        ...
        ValueError: 'foo' is not a graph vertex
    """
    if initial_vertex is not None and initial_vertex not in G:
        raise ValueError(f"'{initial_vertex}' is not a graph vertex")

    if algo not in ("lex_BFS", "lex_UP", "lex_DFS", "lex_DOWN"):
        raise ValueError(f"unknown algorithm '{algo}'")

    cdef size_t n = G.order()
    cdef list sigma = []
    cdef dict predecessor = {}

    cdef bint right = algo in ("lex_BFS", "lex_UP")
    cdef bint decr = algo in ("lex_BFS", "lex_DOWN")

    cdef size_t cur_label = n if decr else -1
    cdef int label_incr = -1 if decr else 1

    # Perform the search
    lexicographic_label = {u: deque() for u in G}
    if initial_vertex is not None:
        # append or appendleft does not matter here, as the deque is empty
        lexicographic_label[initial_vertex].append(cur_label)
    while lexicographic_label:
        u = max(lexicographic_label, key=lexicographic_label.get)
        lexicographic_label.pop(u)
        sigma.append(u)
        cur_label += label_incr
        for v in G.neighbor_iterator(u):  # graphs are considered undirected
            if v in lexicographic_label:
                if right:
                    lexicographic_label[v].append(cur_label)
                else:
                    lexicographic_label[v].appendleft(cur_label)
                predecessor[v] = u

    if reverse:
        sigma.reverse()

    if tree:
        from sage.graphs.digraph import DiGraph
        edges = predecessor.items()
        g = DiGraph([G, edges], format='vertices_and_edges', sparse=True)
        return sigma, g
    return sigma


def _is_valid_lex_BFS_order(G, L):
    r"""
    Check whether ``L`` is a valid LexBFS ordering of the vertices of ``G``.

    Given two vertices `a` and `b` of `G = (V, E)`, we write `a < b` if `a` has
    a smaller label than `b`, and so if `a` is after `b` in the ordering `L`.
    It is proved in [DNB1996]_ that any LexBFS ordering satisfies that,
    if `a < b < c` and `ac \in E` and `bc \not\in E`, then there exists `d\in V`
    such that `c < d`, `db \in E` and `da \not\in E`.

    INPUT:

    - ``G`` -- a sage Graph

    - ``L`` -- list; an ordering of the vertices of `G`

    OUTPUT:

    - ``True`` if ``L`` is a LexBFS ordering of ``G``; ``False`` otherwise

    .. NOTE::

        Loops and multiple edges are ignored for LexBFS ordering and directed
        graphs are considered as undirected graphs.

    .. SEEALSO::

        * :wikipedia:`Lexicographic_breadth-first_search`
        * :meth:`~sage.graphs.generic_graph.GenericGraph.lex_BFS` -- perform a
          lexicographic depth first search (LexBFS) on the graph

    TESTS::

        sage: from sage.graphs.traversals import _is_valid_lex_BFS_order
        sage: G = graphs.PathGraph(3)
        sage: _is_valid_lex_BFS_order(G, [0, 1, 2])
        True
        sage: _is_valid_lex_BFS_order(G, [0, 2, 1])
        False

        sage: G = DiGraph("I?O@??A?CCA?A??C??")
        sage: _is_valid_lex_BFS_order(G, [0, 7, 1, 2, 3, 4, 5, 8, 6, 9])
        False
        sage: _is_valid_lex_BFS_order(G, G.lex_BFS())
        True
        sage: H = G.to_undirected()
        sage: _is_valid_lex_BFS_order(H, G.lex_BFS())
        True
        sage: _is_valid_lex_BFS_order(G, H.lex_BFS())
        True
    """
    # Convert G to a simple undirected graph
    if G.has_loops() or G.has_multiple_edges() or G.is_directed():
        G = G.to_simple(immutable=False, to_undirected=True)

    cdef int n = G.order()

    if set(L) != set(G):
        return False

    cdef dict L_inv = {u: i for i, u in enumerate(L)}
    cdef int pos_a, pos_b, pos_c

    for pos_a in range(n - 1, -1, -1):
        a = L[pos_a]
        for c in G.neighbor_iterator(a):
            pos_c = L_inv[c]
            if pos_c > pos_a:
                continue
            for pos_b in range(pos_c + 1, pos_a):
                b = L[pos_b]
                if G.has_edge(c, b):
                    continue
                if any(L_inv[d] < pos_c and not G.has_edge(d, a)
                       for d in G.neighbor_iterator(b)):
                    # The condition is satisfied for a < b < c
                    continue
                return False
    return True


def lex_BFS(G, reverse=False, tree=False, initial_vertex=None, algorithm="fast"):
    r"""
    Perform a lexicographic breadth first search (LexBFS) on the graph.

    INPUT:

    - ``G`` -- a sage graph

    - ``reverse`` -- boolean (default: ``False``); whether to return the
      vertices in discovery order, or the reverse

    - ``tree`` -- boolean (default: ``False``); whether to return the
      discovery directed tree (each vertex being linked to the one that saw
      it last)

    - ``initial_vertex`` -- (default: ``None``) the first vertex to
      consider

    - ``algorithm`` -- string (default: ``'fast'``); algorithm to use among:

      - ``'slow'`` -- it use the generic algorithm for all the lexicographic
        searchs. See the documentation of the :mod:`~sage.graphs.traversals`
        module for more details.

      - ``'fast'`` -- this algorithm uses the notion of *partition refinement*
        to determine the position of the vertices in the ordering. The time
        complexity of this algorithm is in `O(n + m)`, and our implementation
        follows that complexity for ``SparseGraph``. For ``DenseGraph``,
        the complexity is `O(n^2)`. See [HMPV2000]_ and [TCHP2008]_ for more
        details. This algorithm is also used to compute slice decompositions of
        undirected graphs, a more thorough description can be found in the
        documentation of the
        :mod:`~sage.graphs.graph_decompositions.slice_decomposition` module.

    .. NOTE::

        Loops and multiple edges are ignored during the computation of
        ``lex_BFS`` and directed graphs are converted to undirected graphs.

    .. SEEALSO::

        * :wikipedia:`Lexicographic_breadth-first_search`
        * :mod:`~sage.graphs.graph_decompositions.slice_decomposition` module
          and :meth:`~sage.graphs.graph.Graph.slice_decomposition` -- compute a
          slice decomposition of the graph using an extended lex BFS algorithm
        * :meth:`~sage.graphs.generic_graph.GenericGraph.lex_DFS` -- perform a
          lexicographic depth first search (LexDFS) on the graph
        * :meth:`~sage.graphs.generic_graph.GenericGraph.lex_UP` -- perform a
          lexicographic UP search (LexUP) on the graph
        * :meth:`~sage.graphs.generic_graph.GenericGraph.lex_DOWN` -- perform a
          lexicographic DOWN search (LexDOWN) on the graph

    EXAMPLES:

    A Lex BFS is obviously an ordering of the vertices::

        sage: g = graphs.CompleteGraph(6)
        sage: set(g.lex_BFS()) == set(g)
        True

    Lex BFS ordering of the 3-sun graph::

        sage: g = Graph([(1, 2), (1, 3), (2, 3), (2, 4), (2, 5), (3, 5), (3, 6), (4, 5), (5, 6)])
        sage: g.lex_BFS()
        [1, 2, 3, 5, 4, 6]

    The method also works for directed graphs::

        sage: G = DiGraph([(1, 2), (2, 3), (1, 3)])
        sage: correct_anwsers = [[2, 1, 3], [2, 3, 1]]
        sage: G.lex_BFS(initial_vertex=2, algorithm='slow') in correct_anwsers
        True
        sage: G.lex_BFS(initial_vertex=2, algorithm='fast') in correct_anwsers
        True

    For a Chordal Graph, a reversed Lex BFS is a Perfect Elimination Order::

        sage: g = graphs.PathGraph(3).lexicographic_product(graphs.CompleteGraph(2))
        sage: g.lex_BFS(reverse=True)
        [(2, 1), (2, 0), (1, 1), (1, 0), (0, 1), (0, 0)]

    And the vertices at the end of the tree of discovery are, for chordal
    graphs, simplicial vertices (their neighborhood is a complete graph)::

        sage: g = graphs.ClawGraph().lexicographic_product(graphs.CompleteGraph(2))
        sage: v = g.lex_BFS()[-1]
        sage: peo, tree = g.lex_BFS(initial_vertex = v,  tree=True)
        sage: leaves = [v for v in tree if tree.in_degree(v) ==0]
        sage: all(g.subgraph(g.neighbors(v)).is_clique() for v in leaves)
        True

    Different orderings for different traversals::

        sage: # needs sage.combinat
        sage: G = digraphs.DeBruijn(2,3)
        sage: G.lex_BFS(initial_vertex='000', algorithm='fast')
        ['000', '001', '100', '010', '011', '110', '101', '111']
        sage: G.lex_BFS(initial_vertex='000', algorithm='slow')
        ['000', '001', '100', '010', '011', '110', '101', '111']
        sage: G.lex_DFS(initial_vertex='000')
        ['000', '001', '100', '010', '101', '110', '011', '111']
        sage: G.lex_UP(initial_vertex='000')
        ['000', '001', '010', '101', '110', '111', '011', '100']
        sage: G.lex_DOWN(initial_vertex='000')
        ['000', '001', '100', '011', '010', '110', '111', '101']

    TESTS:

    Computed orderings are valid::

        sage: from sage.graphs.traversals import _is_valid_lex_BFS_order
        sage: G = graphs.RandomChordalGraph(15)
        sage: v0 = ZZ.random_element(G.order())
        sage: L = G.lex_BFS(initial_vertex=v0, algorithm='fast')
        sage: _is_valid_lex_BFS_order(G, L)
        True
        sage: L = G.lex_BFS(initial_vertex=v0, algorithm='slow')
        sage: _is_valid_lex_BFS_order(G, L)
        True
        sage: G = digraphs.RandomDirectedGNP(15, .3)
        sage: v0 = ZZ.random_element(G.order())
        sage: L = G.lex_BFS(initial_vertex=v0, algorithm='fast')
        sage: _is_valid_lex_BFS_order(G, L)
        True
        sage: L = G.lex_BFS(initial_vertex=v0, algorithm='slow')
        sage: _is_valid_lex_BFS_order(G, L)
        True

    Lex BFS ordering of a graph on one vertex::

        sage: Graph(1).lex_BFS(algorithm="fast", tree=True)
        ([0], Digraph on 1 vertex)

    Lex BFS ordering of an empty (di)graph is an empty sequence::

        sage: g = Graph()
        sage: g.lex_BFS(algorithm="fast")
        []
        sage: g.lex_BFS(algorithm="fast", tree=True)
        ([], Digraph on 0 vertices)

    Lex BFS ordering of a symmetric digraph should be the same as the Lex BFS
    ordering of the corresponding undirected graph::

        sage: G = Graph([(1, 2), (1, 3), (2, 3), (2, 4), (2, 5), (3, 5), (3, 6), (4, 5), (5, 6)])
        sage: H = DiGraph(G)
        sage: G.lex_BFS(algorithm="fast") == H.lex_BFS(algorithm="fast")
        True

    ``initial_vertex`` should be a valid graph vertex::

        sage: G = graphs.CompleteGraph(6)
        sage: G.lex_BFS(algorithm="fast", initial_vertex='foo')
        Traceback (most recent call last):
        ...
        ValueError: 'foo' is not a graph vertex
    """
    if initial_vertex is not None and initial_vertex not in G:
        raise ValueError(f"'{initial_vertex}' is not a graph vertex")

    if algorithm == "slow":
        return _lex_order_common(G, "lex_BFS", reverse, tree, initial_vertex)

    if algorithm != "fast":
        raise ValueError(f"unknown algorithm '{algorithm}'")

    cdef size_t n = G.order()

    # For algorithm "fast" we need to convert G to an undirected graph
    if G.is_directed():
        G = G.to_undirected()

    # Initialize variables needed by the fast and slow algorithms
    cdef CGraphBackend Gbackend = <CGraphBackend> G._backend
    cdef CGraph cg = Gbackend.cg()
    cdef list sigma = []
    cdef dict predecessor = {}
    # Initialize variables needed by the fast algorithm
    cdef vector[int] sigma_int
    cdef vector[int] pred
    # Initialize variables needed by the slow algorithm
    cdef dict lexicographic_label
    # Temporary variables
    cdef int vi, i, initial_v_int

    # Perform Lex BFS
    if initial_vertex is not None:
        # we already checked that initial_vertex is in G
        initial_v_int = Gbackend.get_vertex(initial_vertex)
    else:
        initial_v_int = -1
    extended_lex_BFS(cg, sigma_int, NULL, initial_v_int, &pred, NULL, NULL)
    sigma = [ Gbackend.vertex_label(vi) for vi in sigma_int ]
    predecessor = { u: sigma[i] for u, i in zip(sigma, pred) if i != -1 }

    if reverse:
        sigma.reverse()

    if tree:
        from sage.graphs.digraph import DiGraph
        edges = predecessor.items()
        g = DiGraph([G, edges], format='vertices_and_edges', sparse=True)
        return sigma, g
    return sigma


def lex_UP(G, reverse=False, tree=False, initial_vertex=None):
    r"""
    Perform a lexicographic UP search (LexUP) on the graph.

    INPUT:

    - ``G`` -- a sage graph

    - ``reverse`` -- boolean (default: ``False``); whether to return the
      vertices in discovery order, or the reverse

    - ``tree`` -- boolean (default: ``False``); whether to return the
      discovery directed tree (each vertex being linked to the one that saw
      it last)

    - ``initial_vertex`` -- (default: ``None``) the first vertex to
      consider

    .. NOTE::

        Loops and multiple edges are ignored during the computation of
        ``lex_UP`` and directed graphs are converted to undirected graphs.

    .. SEEALSO::

        * :meth:`~sage.graphs.generic_graph.GenericGraph.lex_BFS` -- perform a
          lexicographic breadth first search (LexBFS) on the graph
        * :meth:`~sage.graphs.generic_graph.GenericGraph.lex_DFS` -- perform a
          lexicographic depth first search (LexDFS) on the graph
        * :meth:`~sage.graphs.generic_graph.GenericGraph.lex_DOWN` -- perform a
          lexicographic DOWN search (LexDOWN) on the graph

    ALGORITHM:

    See the documentation of the :mod:`~sage.graphs.traversals` module.

    EXAMPLES:

    A Lex UP is obviously an ordering of the vertices::

        sage: g = graphs.CompleteGraph(6)
        sage: set(g.lex_UP()) == set(g)
        True

    Lex UP ordering of the 3-sun graph::

        sage: g = Graph([(1, 2), (1, 3), (2, 3), (2, 4), (2, 5), (3, 5), (3, 6), (4, 5), (5, 6)])
        sage: g.lex_UP()
        [1, 2, 4, 5, 6, 3]

    The method also works for directed graphs::

        sage: G = DiGraph([(1, 2), (2, 3), (1, 3)])
        sage: correct_anwsers = [[2, 1, 3], [2, 3, 1]]
        sage: G.lex_UP(initial_vertex=2) in correct_anwsers
        True

    Different orderings for different traversals::

        sage: # needs sage.combinat
        sage: G = digraphs.DeBruijn(2,3)
        sage: G.lex_BFS(initial_vertex='000')
        ['000', '001', '100', '010', '011', '110', '101', '111']
        sage: G.lex_DFS(initial_vertex='000')
        ['000', '001', '100', '010', '101', '110', '011', '111']
        sage: G.lex_UP(initial_vertex='000')
        ['000', '001', '010', '101', '110', '111', '011', '100']
        sage: G.lex_DOWN(initial_vertex='000')
        ['000', '001', '100', '011', '010', '110', '111', '101']
    """
    return _lex_order_common(G, "lex_UP", reverse, tree, initial_vertex)


def lex_DFS(G, reverse=False, tree=False, initial_vertex=None):
    r"""
    Perform a lexicographic depth first search (LexDFS) on the graph.

    INPUT:

    - ``G`` -- a sage graph

    - ``reverse`` -- boolean (default: ``False``); whether to return the
      vertices in discovery order, or the reverse

    - ``tree`` -- boolean (default: ``False``); whether to return the
      discovery directed tree (each vertex being linked to the one that saw
      it last)

    - ``initial_vertex`` -- (default: ``None``) the first vertex to
      consider

    .. NOTE::

        Loops and multiple edges are ignored during the computation of
        ``lex_DFS`` and directed graphs are converted to undirected graphs.

    .. SEEALSO::

        * :meth:`~sage.graphs.generic_graph.GenericGraph.lex_BFS` -- perform a
          lexicographic breadth first search (LexBFS) on the graph
        * :meth:`~sage.graphs.generic_graph.GenericGraph.lex_UP` -- perform a
          lexicographic UP search (LexUP) on the graph
        * :meth:`~sage.graphs.generic_graph.GenericGraph.lex_DOWN` -- perform a
          lexicographic DOWN search (LexDOWN) on the graph

    ALGORITHM:

    See the documentation of the :mod:`~sage.graphs.traversals` module.

    EXAMPLES:

    A Lex DFS is obviously an ordering of the vertices::

        sage: g = graphs.CompleteGraph(6)
        sage: set(g.lex_DFS()) == set(g)
        True

    Lex DFS ordering of the 3-sun graph::

        sage: g = Graph([(1, 2), (1, 3), (2, 3), (2, 4), (2, 5), (3, 5), (3, 6), (4, 5), (5, 6)])
        sage: g.lex_DFS()
        [1, 2, 3, 5, 6, 4]

    The method also works for directed graphs::

        sage: G = DiGraph([(1, 2), (2, 3), (1, 3)])
        sage: correct_anwsers = [[2, 1, 3], [2, 3, 1]]
        sage: G.lex_DFS(initial_vertex=2) in correct_anwsers
        True

    Different orderings for different traversals::

        sage: # needs sage.combinat
        sage: G = digraphs.DeBruijn(2,3)
        sage: G.lex_BFS(initial_vertex='000')
        ['000', '001', '100', '010', '011', '110', '101', '111']
        sage: G.lex_DFS(initial_vertex='000')
        ['000', '001', '100', '010', '101', '110', '011', '111']
        sage: G.lex_UP(initial_vertex='000')
        ['000', '001', '010', '101', '110', '111', '011', '100']
        sage: G.lex_DOWN(initial_vertex='000')
        ['000', '001', '100', '011', '010', '110', '111', '101']
    """
    return _lex_order_common(G, "lex_DFS", reverse, tree, initial_vertex)


def lex_DOWN(G, reverse=False, tree=False, initial_vertex=None):
    r"""
    Perform a lexicographic DOWN search (LexDOWN) on the graph.

    INPUT:

    - ``G`` -- a sage graph

    - ``reverse`` -- boolean (default: ``False``); whether to return the
      vertices in discovery order, or the reverse

    - ``tree`` -- boolean (default: ``False``); whether to return the
      discovery directed tree (each vertex being linked to the one that saw
      it)

    - ``initial_vertex`` -- (default: ``None``) the first vertex to
      consider

    .. NOTE::

        Loops and multiple edges are ignored during the computation of
        ``lex_DOWN`` and directed graphs are converted to undirected graphs.

    .. SEEALSO::

        * :meth:`~sage.graphs.generic_graph.GenericGraph.lex_BFS` -- perform a
          lexicographic breadth first search (LexBFS) on the graph
        * :meth:`~sage.graphs.generic_graph.GenericGraph.lex_DFS` -- perform a
          lexicographic depth first search (LexDFS) on the graph
        * :meth:`~sage.graphs.generic_graph.GenericGraph.lex_UP` -- perform a
          lexicographic UP search (LexUP) on the graph

    ALGORITHM:

    See the documentation of the :mod:`~sage.graphs.traversals` module.

    EXAMPLES:

    A Lex DOWN is obviously an ordering of the vertices::

        sage: g = graphs.CompleteGraph(6)
        sage: set(g.lex_DOWN()) == set(g)
        True

    Lex DOWN ordering of the 3-sun graph::

        sage: g = Graph([(1, 2), (1, 3), (2, 3), (2, 4), (2, 5), (3, 5), (3, 6), (4, 5), (5, 6)])
        sage: g.lex_DOWN()
        [1, 2, 3, 4, 6, 5]

    The method also works for directed graphs::

        sage: G = DiGraph([(1, 2), (2, 3), (1, 3)])
        sage: correct_anwsers = [[2, 1, 3], [2, 3, 1]]
        sage: G.lex_DOWN(initial_vertex=2) in correct_anwsers
        True

    Different orderings for different traversals::

        sage: # needs sage.combinat
        sage: G = digraphs.DeBruijn(2,3)
        sage: G.lex_BFS(initial_vertex='000')
        ['000', '001', '100', '010', '011', '110', '101', '111']
        sage: G.lex_DFS(initial_vertex='000')
        ['000', '001', '100', '010', '101', '110', '011', '111']
        sage: G.lex_UP(initial_vertex='000')
        ['000', '001', '010', '101', '110', '111', '011', '100']
        sage: G.lex_DOWN(initial_vertex='000')
        ['000', '001', '100', '011', '010', '110', '111', '101']
    """
    return _lex_order_common(G, "lex_DOWN", reverse, tree, initial_vertex)


def lex_M(self, triangulation=False, labels=False, initial_vertex=None, algorithm=None):
    r"""
    Return an ordering of the vertices according the LexM graph traversal.

    LexM is a lexicographic ordering scheme that is a special type of
    breadth-first-search. LexM can also produce a triangulation of the
    given graph. This functionality is implemented in this method. For
    more details on the algorithms used see Sections 4 (``'lex_M_slow'``)
    and 5.3 (``'lex_M_fast'``) of [RTL76]_.

    .. NOTE::

        This method works only for undirected graphs.

    INPUT:

    - ``triangulation`` -- boolean (default: ``False``); whether to return a
      list of edges that need to be added in order to triangulate the graph

    - ``labels`` -- boolean (default: ``False``); whether to return the labels
      assigned to each vertex

    - ``initial_vertex`` -- (default: ``None``) the first vertex to
      consider

    - ``algorithm`` -- string (default: ``None``); one of the following
      algorithms:

      - ``'lex_M_slow'``: slower implementation of LexM traversal

      - ``'lex_M_fast'``: faster implementation of LexM traversal (works only
        when ``labels`` is set to ``False``)

      - ``None``: Sage chooses the best algorithm: ``'lex_M_slow'`` if
        ``labels`` is set to ``True``, ``'lex_M_fast'`` otherwise.

    OUTPUT:

    Depending on the values of the parameters ``triangulation`` and ``labels``
    the method will return one or more of the following (in that order):

    - an ordering of vertices of the graph according to LexM ordering scheme

    - the labels assigned to each vertex

    - a list of edges that when added to the graph will triangulate it

    EXAMPLES:

    LexM produces an ordering of the vertices::

        sage: g = graphs.CompleteGraph(6)
        sage: ord = g.lex_M(algorithm='lex_M_fast')
        sage: len(ord) == g.order()
        True
        sage: set(ord) == set(g.vertices(sort=False))
        True
        sage: ord = g.lex_M(algorithm='lex_M_slow')
        sage: len(ord) == g.order()
        True
        sage: set(ord) == set(g.vertices(sort=False))
        True

    Both algorithms produce a valid LexM ordering `\alpha` (i.e the neighbors of
    `\alpha(i)` in `G[\{\alpha(i), ..., \alpha(n)\}]` induce a clique)::

        sage: from sage.graphs.traversals import is_valid_lex_M_order
        sage: G = graphs.PetersenGraph()
        sage: ord, F = G.lex_M(triangulation=True, algorithm='lex_M_slow')
        sage: is_valid_lex_M_order(G, ord, F)
        True
        sage: ord, F = G.lex_M(triangulation=True, algorithm='lex_M_fast')
        sage: is_valid_lex_M_order(G, ord, F)
        True

    LexM produces a triangulation of given graph::

        sage: G = graphs.PetersenGraph()
        sage: _, F = G.lex_M(triangulation=True)
        sage: H = Graph(F, format='list_of_edges')
        sage: H.is_chordal()
        True

    LexM ordering of the 3-sun graph::

        sage: g = Graph([(1, 2), (1, 3), (2, 3), (2, 4), (2, 5), (3, 5), (3, 6), (4, 5), (5, 6)])
        sage: g.lex_M()
        [6, 4, 5, 3, 2, 1]

    TESTS:

    ``'lex_M_fast'`` cannot return labels::

        sage: Graph().lex_M(labels=True, algorithm='lex_M_fast')
        Traceback (most recent call last):
        ...
        ValueError: 'lex_M_fast' cannot return labels assigned to vertices

    The method works only for undirected graphs::

        sage: from sage.graphs.traversals import lex_M
        sage: lex_M(DiGraph())
        Traceback (most recent call last):
        ...
        ValueError: input graph must be undirected

    LexM ordering of empty graph::

        sage: G = Graph()
        sage: G.lex_M()
        []

    Parameter ``algorithm`` must be either ``'lex_M_slow'``,
    ``'lex_M_fast'`` or ``None``::

        sage: G = graphs.CompleteGraph(6)
        sage: G.lex_M(algorithm='Bob')
        Traceback (most recent call last):
        ...
        ValueError: unknown algorithm 'Bob'

    ``initial_vertex`` should be a valid graph vertex::

        sage: Graph().lex_M(initial_vertex='foo')
        Traceback (most recent call last):
        ...
        ValueError: 'foo' is not a graph vertex
    """
    if initial_vertex is not None and initial_vertex not in self:
        raise ValueError("'{}' is not a graph vertex".format(initial_vertex))

    if self.is_directed():
        raise ValueError("input graph must be undirected")

    if not algorithm:
        if labels:
            algorithm = "lex_M_slow"
        else:
            algorithm = "lex_M_fast"

    elif algorithm not in ["lex_M_slow", "lex_M_fast"]:
        raise ValueError("unknown algorithm '{}'".format(algorithm))

    if algorithm == "lex_M_slow":
        return lex_M_slow(self, triangulation=triangulation, labels=labels, initial_vertex=initial_vertex)
    if labels:
        raise ValueError("'{}' cannot return labels assigned to vertices".format(algorithm))
    return lex_M_fast(self, triangulation=triangulation, initial_vertex=initial_vertex)


def lex_M_slow(G, triangulation=False, labels=False, initial_vertex=None):
    r"""
    Return an ordering of the vertices according the LexM graph traversal.

    LexM is a lexicographic ordering scheme that is a special type of
    breadth-first-search. This function implements the algorithm described in
    Section 4 of [RTL76]_.

    During the search, the vertices are numbered from `n` to `1`. Let
    `\alpha(i)` denote the vertex numbered `i` and let `\alpha^{-1}(u)` denote
    the number assigned to `u`. Each vertex `u` has also a label, denoted by
    `label(u)`, consisting of a list of numbers selected from `[1,n]` and
    ordered in decreasing order. Given two labels `L_1=[p_1, p_2,\ldots, p_k]`
    and `L_1=[q_1, q_2,\ldots, q_l]`, we define `L_1<L_2` if, for some `j`,
    `p_i==q_i` for `i=1,\ldots,j-1` and `p_j<q_j`, or if `p_i==q_i` for
    `i=1,\ldots,k` and `k<l`. Observe that this is exactly how Python compares
    two lists.

    .. NOTE::

        This method works only for undirected graphs.

    INPUT:

    - ``G`` -- a sage graph

    - ``triangulation`` -- boolean (default: ``False``); whether to return the
      triangulation of the graph produced by the method

    - ``labels`` -- boolean (default: ``False``); whether to return the labels
      assigned to each vertex

    - ``initial_vertex`` -- (default: ``None``) the first vertex to
      consider. If not specified, an arbitrary vertex is chosen.

    OUTPUT:

    Depending on the values of the parameters ``triangulation`` and ``labels``
    the method will return one or more of the following (in that order):

    - the ordering of vertices of `G`

    - the labels assigned to each vertex

    - a list of edges that when added to `G` will produce a triangulation of `G`

    EXAMPLES:

    A LexM ordering is obviously an ordering of the vertices::

        sage: from sage.graphs.traversals import lex_M_slow
        sage: g = graphs.CompleteGraph(6)
        sage: len(lex_M_slow(g)) == g.order()
        True

    LexM ordering and label assignments on the vertices of the 3-sun graph::

        sage: from sage.graphs.traversals import lex_M_slow
        sage: g = Graph([(1, 2), (1, 3), (2, 3), (2, 4), (2, 5), (3, 5), (3, 6), (4, 5), (5, 6)])
        sage: lex_M_slow(g, labels=True)
        ([6, 4, 5, 3, 2, 1],
         {1: [], 2: [5], 3: [5, 4], 4: [4, 2], 5: [4, 3], 6: [3, 2]})

    LexM produces a triangulation of given graph::

        sage: from sage.graphs.traversals import lex_M_slow
        sage: G = graphs.PetersenGraph()
        sage: _, F = lex_M_slow(G, triangulation=True)
        sage: H = G.copy()
        sage: H.add_edges(F)
        sage: H.is_chordal()
        True

    TESTS:

    LexM ordering of empty graph::

        sage: from sage.graphs.traversals import lex_M_slow
        sage: G = Graph()
        sage: lex_M_slow(G)
        []

    The method works only for undirected graphs::

        sage: from sage.graphs.traversals import lex_M_slow
        sage: G = digraphs.Circuit(15)
        sage: lex_M_slow(G)
        Traceback (most recent call last):
        ...
        ValueError: input graph must be undirected

    ``initial_vertex`` should be a valid graph vertex::

        sage: G = graphs.CompleteGraph(6)
        sage: from sage.graphs.traversals import lex_M_slow
        sage: lex_M_slow(G, initial_vertex='foo')
        Traceback (most recent call last):
        ...
        ValueError: 'foo' is not a graph vertex
    """
    if initial_vertex is not None and initial_vertex not in G:
        raise ValueError("'{}' is not a graph vertex".format(initial_vertex))

    if G.is_directed():
        raise ValueError("input graph must be undirected")

    # ==>Initialization
    # Assign empty label to all vertices of G and empty list to F
    cdef list unnumbered_vertices = list(G)
    cdef int n = G.order()
    cdef list alpha = [0] * n
    cdef dict label = {v: [] for v in unnumbered_vertices}
    cdef list F = []
    cdef int i
    cdef set active, reach

    if initial_vertex is not None:
        i = unnumbered_vertices.index(initial_vertex)
        unnumbered_vertices[0], unnumbered_vertices[i] = unnumbered_vertices[i], unnumbered_vertices[0]

    for i in range(n-1, -1, -1):
        # Select: pick an unnumbered vertex u with largest label
        u = unnumbered_vertices[0]
        for v in unnumbered_vertices[1:]:
            if label[u] < label[v]:
                u = v

        unnumbered_vertices.remove(u)
        alpha[i] = u

        # Update: for each vertex v in unnumbered_vertices such that there is a
        # chain u = w_1, w_2, ..., w_{p+1} = v with w_j unnumbered and
        # label(w_j) < label(v) for all j in {2,...,p}. If so, we add i to the
        # label of v and add edge {u,v} to F.
        for v in unnumbered_vertices:

            # We check if there is a chain u = w_1, w_2, ..., w_{p+1} = v with
            # w_j unnumbered and label(w_j) < label(v) for all j in {2, ..., p}
            active = set([w for w in unnumbered_vertices if label[w] < label[v]])
            active.add(v)
            reach = set([u])
            while active and reach and v not in reach:
                w = reach.pop()
                for x in G.neighbor_iterator(w):
                    if x in active:
                        reach.add(x)
                        active.discard(x)

            if v in reach:
                label[v].append(i)
                if triangulation:
                    F.append((u, v))

    if triangulation and labels:
        return alpha, label, F
    elif triangulation:
        return alpha, F
    elif labels:
        return alpha, label
    return alpha


def lex_M_fast(G, triangulation=False, initial_vertex=None):
    r"""
    Return an ordering of the vertices according the LexM graph traversal.

    LexM is a lexicographic ordering scheme that is a special type of
    breadth-first-search. This function implements the algorithm described in
    Section 5.3 of [RTL76]_.

    Note that instead of using labels `1, 2, \ldots, k` and adding `1/2`, we
    use labels `2, 4, \ldots, k` and add `1`, thus avoiding to use floats or
    rationals.

    .. NOTE::

        This method works only for undirected graphs.

    INPUT:

    - ``G`` -- a sage graph

    - ``triangulation`` -- boolean (default: ``False``); whether to return the
      triangulation of given graph produced by the method

    - ``initial_vertex`` -- (default: ``None``) the first vertex to consider

    OUTPUT:

    This method will return an ordering of the vertices of ``G`` according to
    the LexM ordering scheme. Furthermore, if ``triangulation`` is set to
    ``True`` the method also returns a list of edges ``F`` such that when added
    to ``G`` the resulting graph is a triangulation of ``G``.

    EXAMPLES:

    A LexM ordering is obviously an ordering of the vertices::

        sage: from sage.graphs.traversals import lex_M_fast
        sage: g = graphs.CompleteGraph(6)
        sage: len(lex_M_fast(g)) == g.order()
        True

    LexM ordering of the 3-sun graph::

        sage: from sage.graphs.traversals import lex_M_fast
        sage: g = Graph([(1, 2), (1, 3), (2, 3), (2, 4), (2, 5), (3, 5), (3, 6), (4, 5), (5, 6)])
        sage: lex_M_fast(g)
        [6, 4, 5, 3, 2, 1]

    LexM produces a triangulation of given graph::

        sage: from sage.graphs.traversals import lex_M_fast
        sage: G = graphs.PetersenGraph()
        sage: _, F = lex_M_fast(G, triangulation=True)
        sage: H = G.copy()
        sage: H.add_edges(F)
        sage: H.is_chordal()
        True

    TESTS:

    LexM ordering of empty graph::

        sage: from sage.graphs.traversals import lex_M_fast
        sage: G = Graph()
        sage: lex_M_fast(G)
        []

    The method works only for undirected graphs::

        sage: from sage.graphs.traversals import lex_M_fast
        sage: G = digraphs.Circuit(15)
        sage: lex_M_fast(G)
        Traceback (most recent call last):
        ...
        ValueError: input graph must be undirected

    ``initial_vertex`` should be a valid graph vertex::

        sage: G = graphs.CompleteGraph(6)
        sage: from sage.graphs.traversals import lex_M_fast
        sage: lex_M_fast(G, initial_vertex='foo')
        Traceback (most recent call last):
        ...
        ValueError: 'foo' is not a graph vertex
    """
    if initial_vertex is not None and initial_vertex not in G:
        raise ValueError("'{}' is not a graph vertex".format(initial_vertex))

    if G.is_directed():
        raise ValueError("input graph must be undirected")

    # ==> Initialization

    cdef list int_to_v = list(G)
    cdef int i, j, k, v, w, z

    if initial_vertex is not None:
        # We put the initial vertex at first place in the ordering
        i = int_to_v.index(initial_vertex)
        int_to_v[0], int_to_v[i] = int_to_v[i], int_to_v[0]

    cdef short_digraph sd
    init_short_digraph(sd, G, edge_labelled=False, vertex_list=int_to_v,
                       sort_neighbors=False)
    cdef uint32_t* p_tmp
    cdef uint32_t* p_end

    cdef int n = G.order()

    cdef list unnumbered_vertices = list(range(n))

    cdef MemoryAllocator mem = MemoryAllocator()
    cdef int* label = <int*>mem.allocarray(n, sizeof(int))
    cdef int* alpha = <int*>mem.allocarray(n, sizeof(int))
    cdef int* alphainv = <int*>mem.allocarray(n, sizeof(int))
    cdef bint* reached = <bint*>mem.allocarray(n, sizeof(bint))

    for i in range(n):
        label[i] = 2
        alpha[i] = 0
        alphainv[i] = 0
        reached[i] = False

    cdef list F = list()
    cdef dict reach

    k = 2
    for i in range(n - 1, -1, -1):

        # Select: pick an unnumbered vertex v with label(v)==k and assign it
        # number i
        for v in unnumbered_vertices:
            if label[v] == k:
                alpha[i] = v
                alphainv[v] = i
                reached[v] = True
                unnumbered_vertices.remove(v)
                break
        else:
            raise ValueError('unable to find an unnumbered vertex v with label[v] == k')

        # Mark all unnumbered vertices unreached
        for w in unnumbered_vertices:
            reached[w] = False

        reach = dict()
        for j in range(2, k + 1, 2):
            reach[j] = set()

        p_tmp = sd.neighbors[v]
        p_end = sd.neighbors[v + 1]
        while p_tmp < p_end:
            w = p_tmp[0]
            p_tmp += 1
            if alphainv[w]:
                continue
            reach[label[w]].add(w)
            reached[w] = True
            label[w] += 1
            if triangulation:
                F.append((int_to_v[v], int_to_v[w]))

        # Search
        for j in range(2, k + 1, 2):
            while reach[j]:
                w = reach[j].pop()
                p_tmp = sd.neighbors[w]
                p_end = sd.neighbors[w + 1]
                while p_tmp < p_end:
                    z = p_tmp[0]
                    p_tmp += 1
                    if reached[z]:
                        continue
                    reached[z] = True
                    if label[z] > j:
                        reach[label[z]].add(z)
                        label[z] += 1
                        if triangulation:
                            F.append((int_to_v[v], int_to_v[z]))
                    else:
                        reach[j].add(z)

        if unnumbered_vertices:
            # Sort: sort unnumbered vertices by label(w) value
            order = sorted((label[w], w) for w in unnumbered_vertices)

            # Reassign labels as integers from 2 to k, redefining k appropriately
            k = 2
            l, _ = order[0]
            for ll, w in order:
                if l != ll:
                    l = ll
                    k += 2
                label[w] = k

    free_short_digraph(sd)

    cdef list ordering = [int_to_v[alpha[i]] for i in range(n)]

    if triangulation:
        return ordering, F
    return ordering


def is_valid_lex_M_order(G, alpha, F):
    r"""
    Check whether the ordering alpha and the triangulation F are valid for G.

    Given the graph `G = (V, E)` with vertex set `V` and edge set `E`, and the
    set `F` of edges of a triangulation of `G`, let `H = (V, E\cup F)`.
    By induction one can see that for every `i \in \{1, ..., n - 1\}` the
    neighbors of `\alpha(i)` in `H[\{\alpha(i), ..., \alpha(n)\}]` induce a
    clique. The ordering `\alpha` is a perfect elimination ordering of `H`, so
    `H` is chordal. See [RTL76]_ for more details.

    INPUT:

    - ``G`` -- a Graph

    - ``alpha`` -- list; an ordering of the vertices of `G`

    - ``F`` -- an iterable of edges given either as ``(u, v)`` or ``(u, v,
      label)``, the edges of the triangulation of `G`


    TESTS::

        sage: from sage.graphs.traversals import lex_M_slow, is_valid_lex_M_order
        sage: G = graphs.PetersenGraph()
        sage: alpha, F = lex_M_slow(G, triangulation=True)
        sage: is_valid_lex_M_order(G, alpha, F)
        True
        sage: H = Graph(G.edges(sort=False))
        sage: H.add_edges(F)
        sage: H.is_chordal()
        True
        sage: from sage.graphs.traversals import lex_M_fast
        sage: alpha, F = lex_M_fast(G, triangulation=True)
        sage: is_valid_lex_M_order(G, alpha, F)
        True
        sage: H = Graph(G.edges(sort=False))
        sage: H.add_edges(F)
        sage: H.is_chordal()
        True
    """
    H = G.copy()
    H.add_edges(F)
    s_alpha = set(alpha)
    for u in alpha:
        K = H.subgraph(H.neighbors(u))
        s_alpha.discard(u)
        K.delete_vertices([v for v in K if v not in s_alpha])
        if not K.is_clique():
            return False
    return True


def maximum_cardinality_search(G, reverse=False, tree=False, initial_vertex=None):
    r"""
    Return an ordering of the vertices according a maximum cardinality search.

    Maximum cardinality search (MCS) is a graph traversal introduced in
    [TY1984]_. It starts by assigning an arbitrary vertex (or the specified
    ``initial_vertex``) of `G` the last position in the ordering `\alpha`. Every
    vertex keeps a weight equal to the number of its already processed neighbors
    (i.e., already added to `\alpha`), and a vertex of largest such number is
    chosen at each step `i` to be placed in position `n - i` in `\alpha`. This
    ordering can be computed in time `O(n + m)`.

    Time complexity is `O(n+m)` for ``SparseGraph`` and `O(n^2)` for
    ``DenseGraph`` where `n` is the number of vertices and `m` is the number of
    edges.

    When the graph is chordal, the ordering returned by MCS is a *perfect
    elimination ordering*, like :meth:`~sage.graphs.traversals.lex_BFS`. So
    this ordering can be used to recognize chordal graphs. See [He2006]_ for
    more details.

    .. NOTE::

        The current implementation is for connected graphs only.

    INPUT:

    - ``G`` -- a Sage graph

    - ``reverse`` -- boolean (default: ``False``); whether to return the
      vertices in discovery order, or the reverse

    - ``tree`` -- boolean (default: ``False``); whether to also return the
      discovery directed tree (each vertex being linked to the one that saw
      it for the first time)

    - ``initial_vertex`` -- (default: ``None``) the first vertex to consider

    OUTPUT:

    By default, return the ordering `\alpha` as a list. When ``tree`` is
    ``True``, the method returns a tuple `(\alpha, T)`, where `T` is a directed
    tree with the same set of vertices as `G` and a directed edge from `u` to `v`
    if `u` was the first vertex to see `v`.

    EXAMPLES:

    When specified, the ``initial_vertex`` is placed at the end of the ordering,
    unless parameter ``reverse`` is ``True``, in which case it is placed at the
    beginning::

        sage: G = graphs.PathGraph(4)
        sage: G.maximum_cardinality_search(initial_vertex=0)
        [3, 2, 1, 0]
        sage: G.maximum_cardinality_search(initial_vertex=1)
        [0, 3, 2, 1]
        sage: G.maximum_cardinality_search(initial_vertex=2)
        [0, 1, 3, 2]
        sage: G.maximum_cardinality_search(initial_vertex=3)
        [0, 1, 2, 3]
        sage: G.maximum_cardinality_search(initial_vertex=3, reverse=True)
        [3, 2, 1, 0]

    Returning the discovery tree::

        sage: G = graphs.PathGraph(4)
        sage: _, T = G.maximum_cardinality_search(tree=True, initial_vertex=0)
        sage: T.order(), T.size()
        (4, 3)
        sage: T.edges(labels=False, sort=True)
        [(1, 0), (2, 1), (3, 2)]
        sage: _, T = G.maximum_cardinality_search(tree=True, initial_vertex=3)
        sage: T.edges(labels=False, sort=True)
        [(0, 1), (1, 2), (2, 3)]

    TESTS::

        sage: Graph().maximum_cardinality_search()
        []
        sage: Graph(1).maximum_cardinality_search()
        [0]
        sage: Graph(2).maximum_cardinality_search()
        Traceback (most recent call last):
        ...
        ValueError: the input graph is not connected
        sage: graphs.PathGraph(2).maximum_cardinality_search(initial_vertex=17)
        Traceback (most recent call last):
        ...
        ValueError: vertex (17) is not a vertex of the graph
    """
    if tree:
        from sage.graphs.digraph import DiGraph

    cdef int N = G.order()
    if not N:
        return ([], DiGraph()) if tree else []
    if N == 1:
        return (list(G), DiGraph(G)) if tree else list(G)

    cdef list int_to_vertex = list(G)

    if initial_vertex is None:
        initial_vertex = 0
    elif initial_vertex in G:
        initial_vertex = int_to_vertex.index(initial_vertex)
    else:
        raise ValueError("vertex ({0}) is not a vertex of the graph".format(initial_vertex))

    cdef short_digraph sd
    init_short_digraph(sd, G, edge_labelled=False, vertex_list=int_to_vertex,
                       sort_neighbors=False)
    cdef uint32_t** p_vertices = sd.neighbors
    cdef uint32_t* p_tmp
    cdef uint32_t* p_end

    cdef MemoryAllocator mem = MemoryAllocator()
    cdef int* weight = <int*>mem.calloc(N, sizeof(int))
    cdef bint* seen = <bint*>mem.calloc(N, sizeof(bint))
    cdef int* pred = <int *>mem.allocarray(N, sizeof(int))

    cdef int i, u, v
    for i in range(N):
        weight[i] = 0
        seen[i] = False
        pred[i] = i

    # We emulate a heap with decrease key operation using a priority queue.
    # A vertex can be inserted multiple times (up to its degree), but only the
    # first extraction (with maximum weight) matters. The size of the queue will
    # never exceed O(m).
    cdef priority_queue[pair[int, int]] pq
    pq.push((0, initial_vertex))

    # The ordering alpha is feed in reversed order and revert afterword
    cdef list alpha = []

    while not pq.empty():
        _, u = pq.top()
        pq.pop()
        if seen[u]:
            # We use a lazy decrease key mode, so u can be several times in pq
            continue

        alpha.append(int_to_vertex[u])
        seen[u] = True

        p_tmp = p_vertices[u]
        p_end = p_vertices[u + 1]
        while p_tmp < p_end:
            v = p_tmp[0]
            if not seen[v]:
                weight[v] += 1
                pq.push((weight[v], v))
                if pred[v] == v:
                    pred[v] = u
            p_tmp += 1

    free_short_digraph(sd)

    if len(alpha) < N:
        raise ValueError("the input graph is not connected")

    if not reverse:
        alpha.reverse()

    if tree:
        D = DiGraph([int_to_vertex, [(int_to_vertex[i], int_to_vertex[pred[i]])
                                     for i in range(N) if pred[i] != i]],
                    format='vertices_and_edges')
        return alpha, D

    return alpha


cdef inline int swap(int* alpha, int* alpha_inv, int u, int new_pos_u) noexcept:
    """
    Swap positions of u and v in alpha, where v is be the vertex occupying cell
    new_pos_u in alpha.
    """
    cdef int v = alpha[new_pos_u]
    alpha[new_pos_u], alpha[alpha_inv[u]] = u, v
    alpha_inv[u], alpha_inv[v] = alpha_inv[v], alpha_inv[u]
    return v


cdef maximum_cardinality_search_M_short_digraph(short_digraph sd, int initial_vertex,
                                                int* alpha, int* alpha_inv, list F, bint* X):
    r"""
    Compute the ordering and the edges of the triangulation produced by MCS-M.

    Maximum cardinality search M (MCS-M) is an extension of MCS
    (:meth:`~sage.graphs.traversals.maximum_cardinality_search`) in the same way
    that Lex-M (:meth:`~sage.graphs.traversals.lex_M`) is an extension of
    Lex-BFS (:meth:`~sage.graphs.traversalslex_BFS`). That is, in MCS-M when `u`
    receives number `i` at step `n - i + 1`, it increments the weight of all
    unnumbered vertices `v` for which there exists a path between `u` and `v`
    consisting only of unnumbered vertices with weight strictly less than
    `w^-(u)` and `w^-(v)`, where `w^-` is the number of times a vertex has been
    reached during previous iterations. See [BBHP2004]_ for the details of this
    `O(nm)` time algorithm.

    If `G` is not connected, the orderings of each of its connected components
    are added consecutively.

    This method is the core of
    :meth:`~sage.graphs.traversals.maximum_cardinality_search_M`.

    INPUT:

    - ``sd`` -- a ``short_digraph`` as documented in
      :mod:`~sage.graphs.base.static_sparse_graph`

    - ``initial_vertex`` -- integer; initial vertex for the search

    - ``alpha`` -- int array of size `N`; the computed ordering of MCS-M

    - ``alpha_inv`` -- int array of size `N`; the position of vertex ``u`` in
      ``alpha``, that is the inverse function of alpha. So we have
      ``alpha[alpha_inv[u]] == u`` for all `0 \leq u < N - 1`.

    - ``F`` -- list; to be filled with the edges of the triangulation

    - ``X`` -- boolean array of size `N`; ``X[u]`` is set to ``True`` if the
      neighborhood of `u` is a separator of the graph

    TESTS::

        sage: Graph().maximum_cardinality_search_M()
        ([], [], [])
        sage: Graph(1).maximum_cardinality_search_M()
        ([0], [], [])
        sage: graphs.PathGraph(2).maximum_cardinality_search_M(initial_vertex=17)
        Traceback (most recent call last):
        ...
        ValueError: vertex (17) is not a vertex of the graph

    .. TODO::

        Use a fast heap data structure with decrease-key operation.
    """
    # Initialization of data structures
    cdef int N = sd.n
    cdef MemoryAllocator mem = MemoryAllocator()
    # number of times a vertex is reached, initially 0
    cdef int* weight = <int*>mem.calloc(N, sizeof(int))
    # has a vertex been reached, initially False
    cdef bint* reached = <bint*>mem.calloc(N, sizeof(bint))

    cdef int i, u, v, xi
    for i in range(N):
        weight[i] = 0
        alpha[i] = i
        alpha_inv[i] = i
        X[i] = False

    # If an initial vertex is specified, we put it at position 0 in alpha.
    # This way, it will be the first vertex to be considered.
    if initial_vertex:
        swap(alpha, alpha_inv, initial_vertex, 0)

    # variables for the manipulation of the short digraph
    cdef uint32_t** p_vertices = sd.neighbors
    cdef uint32_t* p_tmp
    cdef uint32_t* p_end

    cdef vector[vector[int]] reach
    cdef int s = -1
    cdef int current_pos = N

    while current_pos:

        # Choose an unnumbered vertex of maximum weight.
        # This could be done faster if we had a heap data structure with
        # decrease key operation.
        u = alpha[0]
        for i in range(current_pos):
            v = alpha[i]
            if weight[u] < weight[v]:
                u = v

        # Swap u and the vertex v occupying position current_pos in alpha
        current_pos -= 1
        v = swap(alpha, alpha_inv, u, current_pos)
        reached[u] = True

        # If the weight decreases, the neighborhood of u is a separator
        if weight[u] <= s:
            X[u] = True
        s = weight[u]

        # Search for new edges of the triangulation.
        # We add an edge to the triangulation between u and any unnumbered
        # vertex v such that there is a path (v, x1, x2,... , xk, u) through
        # unnumbered vertices such that count-[xi] < count-[v], 1 <= i <= k. If
        # such an edge is found, we increase the count of v for next round.

        # Mark all unnumbered vertices unreached. These vertices occupy
        # positions 0,..,current_pos-1 in alpha
        reach.clear()
        reach.resize(N)
        for i in range(current_pos):
            v = alpha[i]
            reached[v] = False

        # Initialize reach with unnumbered neighbors of u
        p_tmp = p_vertices[u]
        p_end = p_vertices[u + 1]
        while p_tmp < p_end:
            v = p_tmp[0]
            p_tmp += 1
            if not reached[v]:
                reach[weight[v]].push_back(v)
                reached[v] = True
                weight[v] += 1

        # Search
        for i in range(N):
            while not reach[i].empty():
                xi = reach[i].back()
                reach[i].pop_back()
                p_tmp = p_vertices[xi]
                p_end = p_vertices[xi + 1]
                while p_tmp < p_end:
                    v = p_tmp[0]
                    p_tmp += 1
                    if reached[v]:
                        continue
                    reached[v] = True
                    if i < weight[v]:
                        reach[weight[v]].push_back(v)
                        weight[v] += 1
                        F.append((u, v))
                    else:
                        reach[i].push_back(v)

    reach.clear()


def maximum_cardinality_search_M(G, initial_vertex=None):
    r"""
    Return the ordering and the edges of the triangulation produced by MCS-M.

    Maximum cardinality search M (MCS-M) is an extension of MCS
    (:meth:`~sage.graphs.traversals.maximum_cardinality_search`) in the same way
    that Lex-M (:meth:`~sage.graphs.traversals.lex_M`) is an extension of
    Lex-BFS (:meth:`~sage.graphs.traversals.lex_BFS`). That is, in MCS-M when
    `u` receives number `i` at step `n - i + 1`, it increments the weight of all
    unnumbered vertices `v` for which there exists a path between `u` and `v`
    consisting only of unnumbered vertices with weight strictly less than
    `w^-(u)` and `w^-(v)`, where `w^-` is the number of times a vertex has been
    reached during previous iterations. See [BBHP2004]_ for the details of this
    `O(nm)` time algorithm.

    If `G` is not connected, the orderings of each of its connected components
    are added consecutively. Furthermore, if `G` has `k` connected components
    `C_i` for `0 \leq i < k`, `X` contains at least one vertex of `C_i` for each
    `i \geq 1`. Hence, `|X| \geq k - 1`. In particular, some isolated vertices
    (i.e., of degree 0) can appear in `X` as for such a vertex `x`, we have that
    `G \setminus N(x) = G` is not connected.

    INPUT:

    - ``G`` -- a Sage graph

    - ``initial_vertex`` -- (default: ``None``) the first vertex to consider

    OUTPUT: a tuple `(\alpha, F, X)`, where

    - `\alpha` is the resulting ordering of the vertices. If an initial vertex
      is specified, it gets the last position in the ordering `\alpha`.

    - `F` is the list of edges of a minimal triangulation of `G` according
      `\alpha`

    - `X` is a list of vertices such that for each `x \in X`, the
      neighborhood of `x` in `G` is a separator (i.e., `G \setminus N(x)` is not
      connected). Note that we may have `N(x) = \emptyset` if `G` is not
      connected and `x` has degree 0.

    EXAMPLES:

    Chordal graphs have a perfect elimination ordering, and so the set `F` of
    edges of the triangulation is empty::

        sage: G = graphs.RandomChordalGraph(20)
        sage: alpha, F, X = G.maximum_cardinality_search_M(); F
        []

    The cycle of order 4 is not chordal and so the triangulation has one edge::

        sage: G = graphs.CycleGraph(4)
        sage: alpha, F, X = G.maximum_cardinality_search_M(); len(F)
        1

    The number of edges needed to triangulate of a cycle graph or order `n` is
    `n - 3`, independently of the initial vertex::

        sage: n = randint(3, 20)
        sage: C = graphs.CycleGraph(n)
        sage: _, F, X = C.maximum_cardinality_search_M()
        sage: len(F) == n - 3
        True
        sage: _, F, X = C.maximum_cardinality_search_M(initial_vertex=C.random_vertex())
        sage: len(F) == n - 3
        True

    When an initial vertex is specified, it gets the last position in the
    ordering::

        sage: G = graphs.PathGraph(4)
        sage: G.maximum_cardinality_search_M(initial_vertex=0)
        ([3, 2, 1, 0], [], [2, 3])
        sage: G.maximum_cardinality_search_M(initial_vertex=1)
        ([3, 2, 0, 1], [], [2, 3])
        sage: G.maximum_cardinality_search_M(initial_vertex=2)
        ([0, 1, 3, 2], [], [0, 1])
        sage: G.maximum_cardinality_search_M(initial_vertex=3)
        ([0, 1, 2, 3], [], [0, 1])


    When `G` is not connected, the orderings of each of its connected components
    are added consecutively, the vertices of the component containing the
    initial vertex occupying the last positions::

        sage: G = graphs.CycleGraph(4) * 2
        sage: G.maximum_cardinality_search_M()[0]
        [5, 4, 6, 7, 2, 3, 1, 0]
        sage: G.maximum_cardinality_search_M(initial_vertex=7)[0]
        [2, 1, 3, 0, 5, 6, 4, 7]

    Furthermore, if `G` has `k` connected components, `X` contains at least one
    vertex per connected component, except for the first one, and so at least `k
    - 1` vertices::

        sage: for k in range(1, 5):
        ....:     _, _, X = Graph(k).maximum_cardinality_search_M()
        ....:     if len(X) < k - 1:
        ....:         raise ValueError("something goes wrong")
        sage: G = graphs.RandomGNP(10, .2)
        sage: cc = G.connected_components(sort=False)
        sage: _, _, X = G.maximum_cardinality_search_M()
        sage: len(X) >= len(cc) - 1
        True

    In the example of [BPS2010]_, the triangulation has 3 edges::

        sage: G = Graph({'a': ['b', 'k'], 'b': ['c'], 'c': ['d', 'j', 'k'],
        ....:            'd': ['e', 'f', 'j', 'k'], 'e': ['g'],
        ....:            'f': ['g', 'j', 'k'], 'g': ['j', 'k'], 'h': ['i', 'j'],
        ....:            'i': ['k'], 'j': ['k']})
        sage: _, F, _ = G.maximum_cardinality_search_M(initial_vertex='a')
        sage: len(F)
        3

    TESTS::

        sage: Graph().maximum_cardinality_search_M()
        ([], [], [])
        sage: Graph(1).maximum_cardinality_search_M()
        ([0], [], [])
        sage: graphs.PathGraph(2).maximum_cardinality_search_M(initial_vertex=17)
        Traceback (most recent call last):
        ...
        ValueError: vertex (17) is not a vertex of the graph
    """
    cdef list int_to_vertex = list(G)

    if initial_vertex is None:
        initial_vertex = 0
    elif initial_vertex in G:
        initial_vertex = int_to_vertex.index(initial_vertex)
    else:
        raise ValueError("vertex ({0}) is not a vertex of the graph".format(initial_vertex))

    cdef int N = G.order()
    if not N:
        return ([], [], [])
    if N == 1:
        return (list(G), [], [])

    # Copying the whole graph to obtain the list of neighbors quicker than by
    # calling out_neighbors. This data structure is well documented in the
    # module sage.graphs.base.static_sparse_graph
    cdef short_digraph sd
    init_short_digraph(sd, G, edge_labelled=False, vertex_list=int_to_vertex,
                       sort_neighbors=False)

    cdef MemoryAllocator mem = MemoryAllocator()
    cdef int* alpha = <int*>mem.calloc(N, sizeof(int))
    cdef int* alpha_inv = <int*>mem.calloc(N, sizeof(int))
    cdef bint* X = <bint*>mem.calloc(N, sizeof(bint))
    cdef list F = []

    sig_on()
    maximum_cardinality_search_M_short_digraph(sd, initial_vertex, alpha, alpha_inv, F, X)
    sig_off()

    free_short_digraph(sd)

    cdef int u, v
    return ([int_to_vertex[alpha[u]] for u in range(N)],
            [(int_to_vertex[u], int_to_vertex[v]) for u, v in F],
            [int_to_vertex[u] for u in range(N) if X[u]])
