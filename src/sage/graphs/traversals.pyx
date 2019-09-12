# -*- coding: utf-8 -*-
# cython: binding=True
# distutils: language = c++
r"""
Graph traversals.

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

Methods
-------
"""
# ****************************************************************************
# Copyright (C) 2019 Georgios Giapitzakis Tzintanos <giorgosgiapis@mail.com>
#                    David Coudert <david.coudert@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import collections

from libc.string cimport memset
from sage.ext.memory_allocator cimport MemoryAllocator
from sage.graphs.base.static_sparse_graph cimport short_digraph
from sage.graphs.base.static_sparse_graph cimport init_short_digraph
from sage.graphs.base.static_sparse_graph cimport free_short_digraph
from sage.graphs.base.static_sparse_graph cimport out_degree, has_edge
from libc.stdint cimport uint32_t

def lex_BFS(G, reverse=False, tree=False, initial_vertex=None):
    r"""
    Perform a lexicographic breadth first search (LexBFS) on the graph.

    INPUT:

    - ``G`` -- a sage graph

    - ``reverse`` -- boolean (default: ``False``); whether to return the
      vertices in discovery order, or the reverse

    - ``tree`` -- boolean (default: ``False``); whether to return the
      discovery directed tree (each vertex being linked to the one that saw
      it for the first time)

    - ``initial_vertex`` -- (default: ``None``); the first vertex to
      consider

    ALGORITHM:

    This algorithm maintains for each vertex left in the graph a code
    corresponding to the vertices already removed. The vertex of maximal
    code (according to the lexicographic order) is then removed, and the
    codes are updated.

    Time complexity is `O(n+m)` where `n` is the number of vertices and `m` is
    the number of edges.

    See [CK2008]_ for more details on the algorithm.

    .. SEEALSO::

        * :wikipedia:`Lexicographic_breadth-first_search`
        * :meth:`~sage.graphs.generic_graph.GenericGraph.lex_DFS` -- perform a
          lexicographic depth first search (LexDFS) on the graph
        * :meth:`~sage.graphs.generic_graph.GenericGraph.lex_UP` -- perform a
          lexicographic UP search (LexUP) on the graph
        * :meth:`~sage.graphs.generic_graph.GenericGraph.lex_DOWN` -- perform a
          lexicographic DOWN search (LexDOWN) on the graph

    EXAMPLES:

    A Lex BFS is obviously an ordering of the vertices::

        sage: g = graphs.CompleteGraph(6)
        sage: len(g.lex_BFS()) == g.order()
        True

    Lex BFS ordering of the 3-sun graph::

        sage: g = Graph([(1, 2), (1, 3), (2, 3), (2, 4), (2, 5), (3, 5), (3, 6), (4, 5), (5, 6)])
        sage: g.lex_BFS()
        [1, 2, 3, 5, 4, 6]

    The method also works for directed graphs::

        sage: G = DiGraph([(1, 2), (2, 3), (1, 3)])
        sage: G.lex_BFS(initial_vertex=2)
        [2, 3, 1]

    For a Chordal Graph, a reversed Lex BFS is a Perfect Elimination Order::

        sage: g = graphs.PathGraph(3).lexicographic_product(graphs.CompleteGraph(2))
        sage: g.lex_BFS(reverse=True)  # py2
        [(2, 0), (2, 1), (1, 1), (1, 0), (0, 0), (0, 1)]
        sage: g.lex_BFS(reverse=True)  # py3
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

        sage: G = digraphs.DeBruijn(2,3)
        sage: G.lex_BFS(initial_vertex='000')
        ['000', '001', '100', '010', '011', '110', '101', '111']
        sage: G.lex_DFS(initial_vertex='000')
        ['000', '001', '100', '010', '101', '110', '011', '111']
        sage: G.lex_UP(initial_vertex='000')
        ['000', '001', '010', '101', '110', '111', '011', '100']
        sage: G.lex_DOWN(initial_vertex='000')
        ['000', '001', '100', '011', '010', '110', '111', '101']

    TESTS:

    Lex BFS ordering of a graph on one vertex::

        sage: Graph(1).lex_BFS(tree=True)
        ([0], Digraph on 1 vertex)

    Lex BFS ordering of an empty (di)graph is an empty sequence::

        sage: g = Graph()
        sage: g.lex_BFS()
        []

    Lex BFS ordering of a symmetric digraph should be the same as the Lex BFS
    ordering of the corresponding undirected graph::

        sage: G = Graph([(1, 2), (1, 3), (2, 3), (2, 4), (2, 5), (3, 5), (3, 6), (4, 5), (5, 6)])
        sage: H = DiGraph(G)
        sage: G.lex_BFS() == H.lex_BFS()
        True

    ``initial_vertex`` should be a valid graph vertex::

        sage: G = graphs.CompleteGraph(6)
        sage: G.lex_BFS(initial_vertex='foo')
        Traceback (most recent call last):
        ...
        ValueError: 'foo' is not a graph vertex

    """
    if initial_vertex is not None and initial_vertex not in G:
        raise ValueError("'{}' is not a graph vertex".format(initial_vertex))

    # Loops and multiple edges are not needed in Lex BFS
    if G.allows_loops() or G.allows_multiple_edges():
        G = G.to_simple(immutable=False)

    cdef int nV = G.order()

    if not nV:
        if tree:
            from sage.graphs.digraph import DiGraph
            g = DiGraph(sparse=True)
            return [], g
        else:
            return []

    # Build adjacency list of G
    cdef list int_to_v = list(G)

    cdef short_digraph sd
    init_short_digraph(sd, G, edge_labelled=False, vertex_list=int_to_v)

    # Perform Lex BFS

    cdef list code = [[] for i in range(nV)]

    def l_func(x):
        return code[x]

    cdef list value = []

    # Initialize the predecessors array
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef int *pred = <int *>mem.allocarray(nV, sizeof(int))
    memset(pred, -1, nV * sizeof(int))

    cdef set vertices = set(range(nV))

    cdef int source = 0 if initial_vertex is None else int_to_v.index(initial_vertex)
    code[source].append(nV + 1)

    cdef int now = 1, v, int_neighbor
    while vertices:
        v = max(vertices, key=l_func)
        vertices.remove(v)
        for i in range(0, out_degree(sd, v)):
            int_neighbor = sd.neighbors[v][i]
            if int_neighbor in vertices:
                code[int_neighbor].append(nV - now)
                pred[int_neighbor] = v
        value.append(int_to_v[v])
        now += 1

    free_short_digraph(sd)

    if reverse:
        value.reverse()

    if tree:
        from sage.graphs.digraph import DiGraph
        g = DiGraph(sparse=True)
        g.add_vertices(G)
        edges = [(int_to_v[i], int_to_v[pred[i]]) for i in range(nV) if pred[i] != -1]
        g.add_edges(edges)
        return value, g

    else:
        return value

def lex_UP(G, reverse=False, tree=False, initial_vertex=None):
    r"""
    Perform a lexicographic UP search (LexUP) on the graph.

    INPUT:

    - ``G`` -- a sage graph

    - ``reverse`` -- boolean (default: ``False``); whether to return the
      vertices in discovery order, or the reverse

    - ``tree`` -- boolean (default: ``False``); whether to return the
      discovery directed tree (each vertex being linked to the one that saw
      it for the first time)

    - ``initial_vertex`` -- (default: ``None``); the first vertex to
      consider

    ALGORITHM:

    This algorithm maintains for each vertex left in the graph a code
    corresponding to the vertices already removed. The vertex of maximal
    code (according to the lexicographic order) is then removed, and the
    codes are updated. During the `i`-th iteration of the algorithm `i` is
    appended to the codes of all neighbors of the selected vertex that are left
    in the graph.

    Time complexity is `O(n+m)` where `n` is the number of vertices and `m` is
    the number of edges.

    See [Mil2017]_ for more details on the algorithm.

    .. SEEALSO::

        * :meth:`~sage.graphs.generic_graph.GenericGraph.lex_BFS` -- perform a
          lexicographic breadth first search (LexBFS) on the graph
        * :meth:`~sage.graphs.generic_graph.GenericGraph.lex_DFS` -- perform a
          lexicographic depth first search (LexDFS) on the graph
        * :meth:`~sage.graphs.generic_graph.GenericGraph.lex_DOWN` -- perform a
          lexicographic DOWN search (LexDOWN) on the graph

    EXAMPLES:

    A Lex UP is obviously an ordering of the vertices::

        sage: g = graphs.CompleteGraph(6)
        sage: len(g.lex_UP()) == g.order()
        True

    Lex UP ordering of the 3-sun graph::

        sage: g = Graph([(1, 2), (1, 3), (2, 3), (2, 4), (2, 5), (3, 5), (3, 6), (4, 5), (5, 6)])
        sage: g.lex_UP()
        [1, 2, 4, 5, 6, 3]

    The method also works for directed graphs::

        sage: G = DiGraph([(1, 2), (2, 3), (1, 3)])
        sage: G.lex_UP(initial_vertex=2)
        [2, 3, 1]

    Different orderings for different traversals::

        sage: G = digraphs.DeBruijn(2,3)
        sage: G.lex_BFS(initial_vertex='000')
        ['000', '001', '100', '010', '011', '110', '101', '111']
        sage: G.lex_DFS(initial_vertex='000')
        ['000', '001', '100', '010', '101', '110', '011', '111']
        sage: G.lex_UP(initial_vertex='000')
        ['000', '001', '010', '101', '110', '111', '011', '100']
        sage: G.lex_DOWN(initial_vertex='000')
        ['000', '001', '100', '011', '010', '110', '111', '101']

    TESTS:

    Lex UP ordering of a graph on one vertex::

        sage: Graph(1).lex_UP(tree=True)
        ([0], Digraph on 1 vertex)

    Lex UP ordering of an empty (di)graph is an empty sequence::

        sage: g = Graph()
        sage: g.lex_UP()
        []

    Lex UP ordering of a symmetric digraph should be the same as the Lex UP
    ordering of the corresponding undirected graph::

        sage: G = Graph([(1, 2), (1, 3), (2, 3), (2, 4), (2, 5), (3, 5), (3, 6), (4, 5), (5, 6)])
        sage: H = DiGraph(G)
        sage: G.lex_UP() == H.lex_UP()
        True

    ``initial_vertex`` should be a valid graph vertex::

        sage: G = graphs.CompleteGraph(6)
        sage: G.lex_UP(initial_vertex='foo')
        Traceback (most recent call last):
        ...
        ValueError: 'foo' is not a graph vertex

    """
    if initial_vertex is not None and initial_vertex not in G:
        raise ValueError("'{}' is not a graph vertex".format(initial_vertex))

    # Loops and multiple edges are not needed in Lex UP
    if G.allows_loops() or G.allows_multiple_edges():
        G = G.to_simple(immutable=False)

    cdef int nV = G.order()

    if not nV:
        if tree:
            from sage.graphs.digraph import DiGraph
            g = DiGraph(sparse=True)
            return [], g
        else:
            return []

    # Build adjacency list of G
    cdef list int_to_v = list(G)

    cdef short_digraph sd
    init_short_digraph(sd, G, edge_labelled=False, vertex_list=int_to_v)

    # Perform Lex UP

    cdef list code = [[] for i in range(nV)]

    def l_func(x):
        return code[x]

    cdef list value = []

    # Initialize the predecessors array
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef int *pred = <int *>mem.allocarray(nV, sizeof(int))
    memset(pred, -1, nV * sizeof(int))

    cdef set vertices = set(range(nV))

    cdef int source = 0 if initial_vertex is None else int_to_v.index(initial_vertex)
    code[source].append(nV + 1)

    cdef int now = 1, v, int_neighbor
    while vertices:
        v = max(vertices, key=l_func)
        vertices.remove(v)
        for i in range(0, out_degree(sd, v)):
            int_neighbor = sd.neighbors[v][i]
            if int_neighbor in vertices:
                code[int_neighbor].append(now)
                pred[int_neighbor] = v
        value.append(int_to_v[v])
        now += 1

    free_short_digraph(sd)

    if reverse:
        value.reverse()

    if tree:
        from sage.graphs.digraph import DiGraph
        g = DiGraph(sparse=True)
        g.add_vertices(G)
        edges = [(int_to_v[i], int_to_v[pred[i]]) for i in range(nV) if pred[i] != -1]
        g.add_edges(edges)
        return value, g

    else:
        return value

def lex_DFS(G, reverse=False, tree=False, initial_vertex=None):
    r"""
    Perform a lexicographic depth first search (LexDFS) on the graph.

    INPUT:

    - ``G`` -- a sage graph

    - ``reverse`` -- boolean (default: ``False``); whether to return the
      vertices in discovery order, or the reverse

    - ``tree`` -- boolean (default: ``False``); whether to return the
      discovery directed tree (each vertex being linked to the one that saw
      it for the first time)

    - ``initial_vertex`` -- (default: ``None``); the first vertex to
      consider

    ALGORITHM:

    This algorithm maintains for each vertex left in the graph a code
    corresponding to the vertices already removed. The vertex of maximal
    code (according to the lexicographic order) is then removed, and the
    codes are updated. Lex DFS differs from Lex BFS only in the way codes are
    updated after each iteration.

    Time complexity is `O(n+m)` where `n` is the number of vertices and `m` is
    the number of edges.

    See [CK2008]_ for more details on the algorithm.

    .. SEEALSO::

        * :meth:`~sage.graphs.generic_graph.GenericGraph.lex_BFS` -- perform a
          lexicographic breadth first search (LexBFS) on the graph
        * :meth:`~sage.graphs.generic_graph.GenericGraph.lex_UP` -- perform a
          lexicographic UP search (LexUP) on the graph
        * :meth:`~sage.graphs.generic_graph.GenericGraph.lex_DOWN` -- perform a
          lexicographic DOWN search (LexDOWN) on the graph

    EXAMPLES:

    A Lex DFS is obviously an ordering of the vertices::

        sage: g = graphs.CompleteGraph(6)
        sage: len(g.lex_DFS()) == g.order()
        True

    Lex DFS ordering of the 3-sun graph::

        sage: g = Graph([(1, 2), (1, 3), (2, 3), (2, 4), (2, 5), (3, 5), (3, 6), (4, 5), (5, 6)])
        sage: g.lex_DFS()
        [1, 2, 3, 5, 6, 4]

    The method also works for directed graphs::

        sage: G = DiGraph([(1, 2), (2, 3), (1, 3)])
        sage: G.lex_DFS(initial_vertex=2)
        [2, 3, 1]

    Different orderings for different traversals::

        sage: G = digraphs.DeBruijn(2,3)
        sage: G.lex_BFS(initial_vertex='000')
        ['000', '001', '100', '010', '011', '110', '101', '111']
        sage: G.lex_DFS(initial_vertex='000')
        ['000', '001', '100', '010', '101', '110', '011', '111']
        sage: G.lex_UP(initial_vertex='000')
        ['000', '001', '010', '101', '110', '111', '011', '100']
        sage: G.lex_DOWN(initial_vertex='000')
        ['000', '001', '100', '011', '010', '110', '111', '101']

    TESTS:

    Lex DFS ordering of a graph on one vertex::

        sage: Graph(1).lex_DFS(tree=True)
        ([0], Digraph on 1 vertex)

    Lex DFS ordering of an empty (di)graph is an empty sequence::

        sage: g = Graph()
        sage: g.lex_DFS()
        []

    Lex DFS ordering of a symmetric digraph should be the same as the Lex DFS
    ordering of the corresponding undirected graph::

        sage: G = Graph([(1, 2), (1, 3), (2, 3), (2, 4), (2, 5), (3, 5), (3, 6), (4, 5), (5, 6)])
        sage: H = DiGraph(G)
        sage: G.lex_DFS() == H.lex_DFS()
        True

    ``initial_vertex`` should be a valid graph vertex::

        sage: G = graphs.CompleteGraph(6)
        sage: G.lex_DFS(initial_vertex='foo')
        Traceback (most recent call last):
        ...
        ValueError: 'foo' is not a graph vertex

    """
    if initial_vertex is not None and initial_vertex not in G:
        raise ValueError("'{}' is not a graph vertex".format(initial_vertex))

    # Loops and multiple edges are not needed in Lex DFS
    if G.allows_loops() or G.allows_multiple_edges():
        G = G.to_simple(immutable=False)

    cdef int nV = G.order()

    if not nV:
        if tree:
            from sage.graphs.digraph import DiGraph
            g = DiGraph(sparse=True)
            return [], g
        else:
            return []

    # Build adjacency list of G
    cdef list int_to_v = list(G)

    cdef short_digraph sd
    init_short_digraph(sd, G, edge_labelled=False, vertex_list=int_to_v)

    # Perform Lex DFS

    # We are using deque in order to prepend items in list efficiently
    cdef list code = [collections.deque([]) for i in range(nV)]

    def l_func(x):
        return code[x]

    cdef list value = []

    # initialize the predecessors array
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef int *pred = <int *>mem.allocarray(nV, sizeof(int))
    memset(pred, -1, nV * sizeof(int))

    cdef set vertices = set(range(nV))

    cdef int source = 0 if initial_vertex is None else int_to_v.index(initial_vertex)
    code[source].appendleft(0)

    cdef int now = 1, v, int_neighbor
    while vertices:
        v = max(vertices, key=l_func)
        vertices.remove(v)
        for i in range(0, out_degree(sd, v)):
            int_neighbor = sd.neighbors[v][i]
            if int_neighbor in vertices:
                code[int_neighbor].appendleft(now)
                pred[int_neighbor] = v
        value.append(int_to_v[v])
        now += 1

    free_short_digraph(sd)

    if reverse:
        value.reverse()

    if tree:
        from sage.graphs.digraph import DiGraph
        g = DiGraph(sparse=True)
        g.add_vertices(G)
        edges = [(int_to_v[i], int_to_v[pred[i]]) for i in range(nV) if pred[i] != -1]
        g.add_edges(edges)
        return value, g

    else:
        return value

def lex_DOWN(G, reverse=False, tree=False, initial_vertex=None):
    r"""
    Perform a lexicographic DOWN search (LexDOWN) on the graph.

    INPUT:

    - ``G`` -- a sage graph

    - ``reverse`` -- boolean (default: ``False``); whether to return the
      vertices in discovery order, or the reverse

    - ``tree`` -- boolean (default: ``False``); whether to return the
      discovery directed tree (each vertex being linked to the one that saw
      it for the first time)

    - ``initial_vertex`` -- (default: ``None``); the first vertex to
      consider

    ALGORITHM:

    This algorithm maintains for each vertex left in the graph a code
    corresponding to the vertices already removed. The vertex of maximal
    code (according to the lexicographic order) is then removed, and the
    codes are updated. During the `i`-th iteration of the algorithm `n-i` is
    prepended to the codes of all neighbors of the selected vertex that are left
    in the graph.

    Time complexity is `O(n+m)` where `n` is the number of vertices and `m` is
    the number of edges.

    See [Mil2017]_ for more details on the algorithm.

    .. SEEALSO::

        * :meth:`~sage.graphs.generic_graph.GenericGraph.lex_BFS` -- perform a
          lexicographic breadth first search (LexBFS) on the graph
        * :meth:`~sage.graphs.generic_graph.GenericGraph.lex_DFS` -- perform a
          lexicographic depth first search (LexDFS) on the graph
        * :meth:`~sage.graphs.generic_graph.GenericGraph.lex_UP` -- perform a
          lexicographic UP search (LexUP) on the graph

    EXAMPLES:

    A Lex DOWN is obviously an ordering of the vertices::

        sage: g = graphs.CompleteGraph(6)
        sage: len(g.lex_DOWN()) == g.order()
        True

    Lex DOWN ordering of the 3-sun graph::

        sage: g = Graph([(1, 2), (1, 3), (2, 3), (2, 4), (2, 5), (3, 5), (3, 6), (4, 5), (5, 6)])
        sage: g.lex_DOWN()
        [1, 2, 3, 4, 6, 5]

    The method also works for directed graphs::

        sage: G = DiGraph([(1, 2), (2, 3), (1, 3)])
        sage: G.lex_DOWN(initial_vertex=2)
        [2, 3, 1]

    Different orderings for different traversals::

        sage: G = digraphs.DeBruijn(2,3)
        sage: G.lex_BFS(initial_vertex='000')
        ['000', '001', '100', '010', '011', '110', '101', '111']
        sage: G.lex_DFS(initial_vertex='000')
        ['000', '001', '100', '010', '101', '110', '011', '111']
        sage: G.lex_UP(initial_vertex='000')
        ['000', '001', '010', '101', '110', '111', '011', '100']
        sage: G.lex_DOWN(initial_vertex='000')
        ['000', '001', '100', '011', '010', '110', '111', '101']

    TESTS:

    Lex DOWN ordering of a graph on one vertex::

        sage: Graph(1).lex_DOWN(tree=True)
        ([0], Digraph on 1 vertex)

    Lex DOWN ordering of an empty (di)graph is an empty sequence::

        sage: g = Graph()
        sage: g.lex_DOWN()
        []

    Lex DOWN ordering of a symmetric digraph should be the same as the Lex DOWN
    ordering of the corresponding undirected graph::

        sage: G = Graph([(1, 2), (1, 3), (2, 3), (2, 4), (2, 5), (3, 5), (3, 6), (4, 5), (5, 6)])
        sage: H = DiGraph(G)
        sage: G.lex_DOWN() == H.lex_DOWN()
        True

    ``initial_vertex`` should be a valid graph vertex::

        sage: G = graphs.CompleteGraph(6)
        sage: G.lex_DOWN(initial_vertex='foo')
        Traceback (most recent call last):
        ...
        ValueError: 'foo' is not a graph vertex

    """
    if initial_vertex is not None and initial_vertex not in G:
        raise ValueError("'{}' is not a graph vertex".format(initial_vertex))

    # Loops and multiple edges are not needed in Lex DOWN
    if G.allows_loops() or G.allows_multiple_edges():
        G = G.to_simple(immutable=False)

    cdef int nV = G.order()

    if not nV:
        if tree:
            from sage.graphs.digraph import DiGraph
            g = DiGraph(sparse=True)
            return [], g
        else:
            return []

    # Build adjacency list of G
    cdef list int_to_v = list(G)

    cdef short_digraph sd
    init_short_digraph(sd, G, edge_labelled=False, vertex_list=int_to_v)

    # Perform Lex DOWN

    # We are using deque in order to prepend items in list efficiently
    cdef list code = [collections.deque([]) for i in range(nV)]

    def l_func(x):
        return code[x]

    cdef list value = []

    # initialize the predecessors array
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef int *pred = <int *>mem.allocarray(nV, sizeof(int))
    memset(pred, -1, nV * sizeof(int))

    cdef set vertices = set(range(nV))

    cdef int source = 0 if initial_vertex is None else int_to_v.index(initial_vertex)
    code[source].appendleft(0)

    cdef int now = 1, v, int_neighbor
    while vertices:
        v = max(vertices, key=l_func)
        vertices.remove(v)
        for i in range(0, out_degree(sd, v)):
            int_neighbor = sd.neighbors[v][i]
            if int_neighbor in vertices:
                code[int_neighbor].appendleft(nV - now)
                pred[int_neighbor] = v
        value.append(int_to_v[v])
        now += 1

    free_short_digraph(sd)

    if reverse:
        value.reverse()

    if tree:
        from sage.graphs.digraph import DiGraph
        g = DiGraph(sparse=True)
        g.add_vertices(G)
        edges = [(int_to_v[i], int_to_v[pred[i]]) for i in range(nV) if pred[i] != -1]
        g.add_edges(edges)
        return value, g

    else:
        return value

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

    - ``initial_vertex`` -- (default: ``None``); the first vertex to
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
        sage: set(ord) == set(g.vertices())
        True
        sage: ord = g.lex_M(algorithm='lex_M_slow')
        sage: len(ord) == g.order()
        True
        sage: set(ord) == set(g.vertices())
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

        sage: G = graphs.CompleteGraph(6)
        sage: G.lex_M(labels=True, algorithm='lex_M_fast')
        Traceback (most recent call last):
        ...
        ValueError: 'lex_M_fast' cannot return labels assigned to vertices

    The method works only for undirected graphs::

        sage: G = digraphs.Circuit(15)
        sage: G.lex_M()
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

        sage: G = graphs.CompleteGraph(6)
        sage: G.lex_M(initial_vertex='foo')
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
    else:
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

    - ``initial_vertex`` -- (default: ``None``); the first vertex to
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
    else:
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

    - ``initial_vertex`` -- (default: ``None``); the first vertex to consider

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
    init_short_digraph(sd, G, edge_labelled=False, vertex_list=int_to_v)
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
    for i in range(n-1, -1, -1):

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
            order = sorted( (label[w], w) for w in unnumbered_vertices )

            # Reassign labels as integers from 2 to k, redefining k appropriately
            k = 2
            l,_ = order[0]
            for ll,w in order:
                if l != ll:
                    l = ll
                    k += 2
                label[w] = k

    free_short_digraph(sd)

    cdef list ordering = [int_to_v[alpha[i]] for i in range(n)]

    if triangulation:
        return ordering, F
    else:
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

    INPUTS:

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
