# cython: binding=True
r"""
Weakly chordal graphs

This module deals with everything related to weakly chordal graphs. It currently
contains the following functions:

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~sage.graphs.weakly_chordal.is_long_hole_free` | Test whether ``g`` contains an induced cycle of length at least 5.
    :meth:`~sage.graphs.weakly_chordal.is_long_antihole_free` | Test whether ``g`` contains an induced anticycle of length at least 5.
    :meth:`~sage.graphs.weakly_chordal.is_weakly_chordal` | Test whether ``g`` is weakly chordal.

Author:

- Birk Eisermann (initial implementation)
- Nathann Cohen (some doc and optimization)
- David Coudert (remove recursion)


Methods
-------
"""

# ****************************************************************************
#       Copyright (c) 2012 Birk Eisermann <eisermbi@fastmail.fm>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from memory_allocator cimport MemoryAllocator
from sage.data_structures.bitset_base cimport *
from sage.graphs.base.static_sparse_backend cimport StaticSparseCGraph
from sage.graphs.base.static_sparse_backend cimport StaticSparseBackend
from sage.graphs.base.static_sparse_graph cimport short_digraph
from sage.graphs.base.static_sparse_graph cimport init_short_digraph
from sage.graphs.base.static_sparse_graph cimport free_short_digraph
from sage.graphs.base.static_sparse_graph cimport out_degree


cdef inline int has_edge(bitset_t bs, int u, int v, int n) noexcept:
    return bitset_in(bs, u * n + v)


cdef inline is_long_hole_free_process(g, short_digraph sd, bitset_t dense_graph,
                                      list id_label, int* path, int* InPath,
                                      int* neighbor_index, set VisitedP3,
                                      bint certificate,
                                      int a, int b, int c, int n):
    """
    This method is part of method ``is_long_hole_free``.

    EXAMPLES::

        sage: g = graphs.PetersenGraph()
        sage: g.is_long_hole_free()
        False
    """
    cdef int d, u, v, i
    cdef Py_ssize_t path_top = 2
    cdef list C
    cdef dict C_index

    path[2] = c
    InPath[c] = path_top  # c is the (i+1)-th vertex at position i
    neighbor_index[c] = 0
    VisitedP3.add((a, b, c))
    VisitedP3.add((c, b, a))

    while path_top >= 2:
        a = path[path_top - 2]
        b = path[path_top - 1]
        c = path[path_top]

        if neighbor_index[c] < out_degree(sd, c):
            d = sd.neighbors[c][neighbor_index[c]]
            neighbor_index[c] += 1
            if not has_edge(dense_graph, d, a, n) and not has_edge(dense_graph, d, b, n):
                # a-b-c-d form an induced path P_4

                if InPath[d] != -1:
                    # d is already contained in InPath
                    # HOLE FOUND !!!
                    if certificate:
                        # We extract the hole and relabel it on-the-fly with the
                        # vertices' real name
                        C = [id_label[path[i]] for i in range(InPath[d], path_top + 1)]
                        C_index = {label: i for i, label in enumerate(C)}

                        # At this step C[0]C[1]..... is a cycle such that any 4
                        # consecutive vertices induce a P4. C may not be an
                        # induced cycle, so we extract one from it.

                        # To do so, we look for the *shortest* edge C[i]C[j]
                        # between two nonconsecutive vertices of C, where the
                        # length is the difference |i-j|.
                        #
                        # C[i]...C[j] is necessarily an induced cycle.

                        gg = g.subgraph(C, immutable=False)
                        gg.delete_edges(zip(C[:-1], C[1:]))

                        def dist(X):
                            return abs(C_index[X[0]] - C_index[X[1]])

                        label_u, label_v = min(gg.edge_iterator(labels=False), key=dist)
                        u, v = C_index[label_u], C_index[label_v]

                        # Return the answer
                        return False, g.subgraph(C[min(u, v): max(u, v) + 1])

                    else:
                        return False, None

                elif (b, c, d) not in VisitedP3:
                    # search for another P_4
                    path_top += 1
                    path[path_top] = d
                    InPath[d] = path_top
                    neighbor_index[d] = 0
                    VisitedP3.add((b, c, d))
                    VisitedP3.add((d, c, b))

        else:
            # We are done with c. We trackback
            path_top -= 1
            InPath[c] = -1

    return True, []


def is_long_hole_free(g, certificate=False):
    r"""
    Test whether ``g`` contains an induced cycle of length at least 5.

    INPUT:

    - ``certificate`` -- boolean (default: ``False``)

      Whether to return a certificate. When ``certificate = True``, then
      the function returns

      * ``(True, [])`` if ``g`` does not contain such a cycle.
        For this case, it is not known how to provide a certificate.
      * ``(False, Hole)`` if ``g`` contains an induced cycle of length at
        least 5. ``Hole`` returns this cycle.

      If ``certificate = False``, the function returns just ``True`` or
      ``False`` accordingly.

    ALGORITHM:

    This algorithm tries to find a cycle in the graph of all induced `P_4` of
    `g`, where two copies `P` and `P'` of `P_4` are adjacent if there exists a
    (not necessarily induced) copy of `P_5=u_1u_2u_3u_4u_5` such that
    `P=u_1u_2u_3u_4` and `P'=u_2u_3u_4u_5`.

    This is done through a depth-first-search. For efficiency, the auxiliary
    graph is constructed on-the-fly and never stored in memory.

    The run time of this algorithm is `O(n+m^2)` for ``SparseGraph`` and
    `O(n^2 + m^2)` for ``DenseGraph`` [NP2007]_ (where `n` is the number of
    vertices and `m` is the number of edges of the graph).

    EXAMPLES:

    The Petersen Graph contains a hole::

        sage: g = graphs.PetersenGraph()
        sage: g.is_long_hole_free()
        False

    The following graph contains a hole, which we want to display::

        sage: g = graphs.FlowerSnark()
        sage: r,h = g.is_long_hole_free(certificate=True)
        sage: r
        False
        sage: Graph(h).is_isomorphic(graphs.CycleGraph(h.order()))
        True

    TESTS:

    Another graph with vertices 2, ..., 8, 10::

        sage: g = Graph({2:[3,8],3:[2,4],4:[3,8,10],5:[6,10],6:[5,7],7:[6,8],8:[2,4,7,10],10:[4,5,8]})
        sage: r,hole = g.is_long_hole_free(certificate=True)
        sage: r
        False
        sage: hole
        Subgraph of (): Graph on 5 vertices
        sage: hole.is_isomorphic(graphs.CycleGraph(hole.order()))
        True

        sage: graphs.EmptyGraph().is_long_hole_free()
        True

    Immutable graphs::

        sage: G = graphs.RandomGNP(10, .7)
        sage: G._backend
        <sage.graphs.base.sparse_graph.SparseGraphBackend ...>
        sage: H = Graph(G, immutable=True)
        sage: H._backend
        <sage.graphs.base.static_sparse_backend.StaticSparseBackend ...>
        sage: G.is_long_hole_free() == H.is_long_hole_free()
        True
    """
    g._scream_if_not_simple()

    if g.order() < 5:
        return (True, []) if certificate else True

    cdef int u, v, w, vv, ww

    # Make a copy of the graph as a short_digraph. This data structure is well
    # documented in the module sage.graphs.base.static_sparse_graph.
    # Vertices are relabeled in 0..n-1
    cdef int n = g.order()
    cdef list id_label
    cdef StaticSparseCGraph cg
    cdef short_digraph sd
    if isinstance(g, StaticSparseBackend):
        cg = <StaticSparseCGraph> g._cg
        sd = <short_digraph> cg.g
        id_label = cg._vertex_to_labels
    else:
        id_label = list(g)
        init_short_digraph(sd, g, edge_labelled=False, vertex_list=id_label)

    # Make a dense copy of the graph for quick adjacency tests
    cdef bitset_t dense_graph
    bitset_init(dense_graph, n * n)
    bitset_set_first_n(dense_graph, 0)
    for u in range(n):
        for vv in range(out_degree(sd, u)):
            v = sd.neighbors[u][vv]
            bitset_add(dense_graph, u * n + v)
            bitset_add(dense_graph, v * n + u)

    # Allocate some data structures
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef int* path = <int*> mem.allocarray(n, sizeof(int))
    cdef int* InPath = <int*> mem.allocarray(n, sizeof(int))
    for u in range(n):
        InPath[u] = -1

    cdef int* neighbor_index = <int*> mem.allocarray(n, sizeof(int))

    cdef set VisitedP3 = set()  # stores triples (u,v,w) which represent visited paths of length 3

    # main algorithm
    # For all triples u,v,w of vertices such that uvw is a P_3
    for u in range(n):
        # u is the first vertex of the path, at position 0
        path[0] = u
        InPath[u] = 0
        for vv in range(out_degree(sd, u)):
            v = sd.neighbors[u][vv]
            # v is the second vertex of the path, at position 1
            path[1] = v
            InPath[v] = 1
            for ww in range(out_degree(sd, v)):
                w = sd.neighbors[v][ww]
                if u != w and not has_edge(dense_graph, u, w, n) and (u, v, w) not in VisitedP3:

                    res, hole = is_long_hole_free_process(g, sd, dense_graph, id_label,
                                                          path, InPath, neighbor_index, VisitedP3,
                                                          certificate, u, v, w, n)

                    if not res:
                        # We release memory before returning the result
                        if isinstance(g, StaticSparseBackend):
                            free_short_digraph(sd)
                        bitset_free(dense_graph)

                        if certificate:
                            return False, hole
                        return False

            InPath[v] = -1
        InPath[u] = -1

    # Release memory
    if not isinstance(g, StaticSparseBackend):
        free_short_digraph(sd)
    bitset_free(dense_graph)

    if certificate:
        return True, []
    return True


cdef inline is_long_antihole_free_process(g, short_digraph sd, bitset_t dense_graph,
                                          list id_label, int* path, int* InPath,
                                          int* neighbor_index, set VisitedP3,
                                          bint certificate,
                                          int a, int b, int c, int n):
    """
    This method is part of method ``is_long_antihole_free``.

    EXAMPLES::

        sage: g = graphs.PetersenGraph()
        sage: g.is_long_antihole_free()
        False
    """
    cdef int d, u, v, i
    cdef Py_ssize_t path_top = 2
    cdef list C
    cdef dict C_index

    path[2] = c
    InPath[c] = 2  # c is the (i+1)-th vertex at position i
    neighbor_index[b] = 0
    VisitedP3.add((a, b, c))
    VisitedP3.add((c, b, a))

    while path_top >= 2:
        # We consider the antichain a,b,c
        a = path[path_top - 2]
        b = path[path_top - 1]
        c = path[path_top]

        if neighbor_index[b] < out_degree(sd, b):
            d = sd.neighbors[b][neighbor_index[b]]
            neighbor_index[b] += 1
            if has_edge(dense_graph, d, a, n) and not has_edge(dense_graph, d, c, n):
                # We found a neighbor of a and b that is not adjacent to c
                if InPath[d] != -1:
                    if certificate:
                        # Calculation of induced cycle in complement
                        # Relabel it on-the-fly with the vertices' real name
                        C = [id_label[path[i]] for i in range(InPath[d], path_top + 1)]
                        C_index = {label: i for i, label in enumerate(C)}

                        # At this step C[0]C[1]..... is an anticycle such that
                        # any 4 consecutive vertices induce the complement of a
                        # P4. C may not be an induced anticycle, so we extract
                        # one from it.

                        # To do so, we look for the *shortest* nonedge C[i]C[j]
                        # between two nonconsecutive vertices of C, where the
                        # length is the difference |i-j|.
                        #
                        # C[i]...C[j] is necessarily an induced anticycle.

                        gg = g.subgraph(C, immutable=False).complement()
                        gg.delete_edges(zip(C[:-1], C[1:]))

                        def dist(X):
                            return abs(C_index[X[0]] - C_index[X[1]])

                        label_u, label_v = min(gg.edge_iterator(labels=False), key=dist)
                        u, v = C_index[label_u], C_index[label_v]

                        # Return the answer
                        return False, g.subgraph(C[min(u, v): max(u, v) + 1])

                    else:
                        return False, []

                elif (b, c, d) not in VisitedP3:
                    path_top += 1
                    path[path_top] = d
                    InPath[d] = path_top
                    neighbor_index[c] = 0
                    VisitedP3.add((b, c, d))
                    VisitedP3.add((d, c, b))

        else:
            # We trackback
            path_top -= 1
            InPath[c] = -1

    return True, []


def is_long_antihole_free(g, certificate=False):
    r"""
    Test whether the given graph contains an induced subgraph that is
    isomorphic to the complement of a cycle of length at least 5.

    INPUT:

    - ``certificate`` -- boolean (default: ``False``)

      Whether to return a certificate. When ``certificate = True``, then
      the function returns

      * ``(False, Antihole)`` if ``g`` contains an induced complement
        of a cycle of length at least 5 returned as ``Antihole``.
      * ``(True, [])`` if ``g`` does not contain an induced complement of
        a cycle of length at least 5.
        For this case it is not known how to provide a certificate.

      When ``certificate = False``, the function returns just ``True`` or
      ``False`` accordingly.

    ALGORITHM:

    This algorithm tries to find a cycle in the graph of all induced
    `\overline{P_4}` of `g`, where two copies `\overline{P}` and `\overline{P'}`
    of `\overline{P_4}` are adjacent if there exists a (not necessarily induced)
    copy of `\overline{P_5}=u_1u_2u_3u_4u_5` such that
    `\overline{P}=u_1u_2u_3u_4` and `\overline{P'}=u_2u_3u_4u_5`.

    This is done through a depth-first-search. For efficiency, the auxiliary
    graph is constructed on-the-fly and never stored in memory.

    The run time of this algorithm is `O(n+m^2)` for ``SparseGraph`` and
    `O(n^2\log{m} + m^2)` for ``DenseGraph`` [NP2007]_ (where `n` is the number
    of vertices and `m` is the number of edges of the graph).

    EXAMPLES:

    The Petersen Graph contains an antihole::

        sage: g = graphs.PetersenGraph()
        sage: g.is_long_antihole_free()
        False

    The complement of a cycle is an antihole::

        sage: g = graphs.CycleGraph(6).complement()
        sage: r,a = g.is_long_antihole_free(certificate=True)
        sage: r
        False
        sage: a.complement().is_isomorphic(graphs.CycleGraph(6))
        True

    TESTS:

    Further tests::

        sage: g = Graph({0:[6,7],1:[7,8],2:[8,9],3:[9,10],4:[10,11],5:[11,6],6:[0,5,7],7:[0,1,6],8:[1,2,9],9:[2,3,8],10:[3,4,11],11:[4,5,10]}).complement()
        sage: r,a = g.is_long_antihole_free(certificate=True)
        sage: r
        False
        sage: a.complement().is_isomorphic(graphs.CycleGraph(9))
        True

        sage: graphs.EmptyGraph().is_long_hole_free()
        True

    Immutable graphs::

        sage: G = graphs.RandomGNP(10, .7)
        sage: G._backend
        <sage.graphs.base.sparse_graph.SparseGraphBackend ...>
        sage: H = Graph(G, immutable=True)
        sage: H._backend
        <sage.graphs.base.static_sparse_backend.StaticSparseBackend ...>
        sage: G.is_long_antihole_free() == H.is_long_antihole_free()
        True
    """
    g._scream_if_not_simple()

    if g.order() < 5:
        return (True, []) if certificate else True

    cdef int u, v, w, vv, ww

    # Make a copy of the graph as a short_digraph. This data structure is well
    # documented in the module sage.graphs.base.static_sparse_graph.
    # Vertices are relabeled in 0..n-1
    cdef int n = g.order()
    cdef list id_label
    cdef StaticSparseCGraph cg
    cdef short_digraph sd
    if isinstance(g, StaticSparseBackend):
        cg = <StaticSparseCGraph> g._cg
        sd = <short_digraph> cg.g
        id_label = cg._vertex_to_labels
    else:
        id_label = list(g)
        init_short_digraph(sd, g, edge_labelled=False, vertex_list=id_label)

    # Make a dense copy of the graph for quick adjacency tests
    cdef bitset_t dense_graph
    bitset_init(dense_graph, n * n)
    bitset_set_first_n(dense_graph, 0)
    for u in range(n):
        for vv in range(out_degree(sd, u)):
            v = sd.neighbors[u][vv]
            bitset_add(dense_graph, u * n + v)
            bitset_add(dense_graph, v * n + u)

    # Allocate some data structures
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef int* path = <int*> mem.allocarray(n, sizeof(int))
    cdef int* InPath = <int*> mem.allocarray(n, sizeof(int))
    for u in range(n):
        InPath[u] = -1

    cdef int* neighbor_index = <int*> mem.allocarray(n, sizeof(int))

    cdef set VisitedP3 = set()  # stores triples (u,v,w) which represent visited paths of length 3

    # main algorithm
    # For all triples u,v,w of vertices such that uvw is a complement of P_3
    for u in range(n):
        # u is the first vertex of the path, at position 0
        path[1] = u
        InPath[u] = 1
        for v in range(n):
            if v == u or has_edge(dense_graph, u, v, n):
                continue
            path[0] = v
            InPath[v] = 0
            for ww in range(out_degree(sd, v)):
                w = sd.neighbors[v][ww]
                if v < w and not has_edge(dense_graph, u, w, n) and (v, u, w) not in VisitedP3:

                    res, antihole = is_long_antihole_free_process(g, sd, dense_graph, id_label,
                                                                  path, InPath, neighbor_index,
                                                                  VisitedP3, certificate,
                                                                  v, u, w, n)

                    if not res:
                        # We release memory before returning the result
                        if isinstance(g, StaticSparseBackend):
                            free_short_digraph(sd)
                        bitset_free(dense_graph)

                        if certificate:
                            return False, antihole
                        return False

            InPath[v] = -1
        InPath[u] = -1

    # Release memory
    if not isinstance(g, StaticSparseBackend):
        free_short_digraph(sd)
    bitset_free(dense_graph)

    if certificate:
        return True, []
    return True


def is_weakly_chordal(g, certificate=False):
    r"""
    Test whether the given graph is weakly chordal, i.e., the graph and its
    complement have no induced cycle of length at least 5.

    INPUT:

    - ``certificate`` -- boolean (default: ``False``); whether to
      return a certificate. If ``certificate = False``, return ``True`` or
      ``False`` according to the graph. If ``certificate = True``, return

      * ``(False, forbidden_subgraph)`` when the graph contains a
        forbidden subgraph H, this graph is returned.
      * ``(True, [])`` when the graph is weakly chordal.
          For this case, it is not known how to provide a certificate.

    ALGORITHM:

    This algorithm checks whether the graph ``g`` or its complement
    contain an induced cycle of length at least 5.

    Using is_long_hole_free() and is_long_antihole_free() yields a run time
    of `O(n+m^2)` for ``SparseGraph`` and `O(n^2 + m^2)` for ``DenseGraph``
    (where `n` is the number of vertices and `m` is the number of edges of the
    graph).

    EXAMPLES:

    The Petersen Graph is not weakly chordal and contains a hole::

        sage: g = graphs.PetersenGraph()
        sage: r,s = g.is_weakly_chordal(certificate=True)
        sage: r
        False
        sage: l = s.order()
        sage: s.is_isomorphic(graphs.CycleGraph(l))
        True

    TESTS::

        sage: graphs.EmptyGraph().is_weakly_chordal()
        True
    """
    if g.order() < 5:
        return (True, []) if certificate else True

    if certificate:
        r, forbid_subgr = g.is_long_hole_free(certificate=True)
        if not r:
            return False, forbid_subgr

        return g.is_long_antihole_free(certificate=True)

    return g.is_long_hole_free() and g.is_long_antihole_free()
