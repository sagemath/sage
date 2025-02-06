# cython: binding=True
# distutils: language = c++
r"""
Path enumeration

This module is meant for all functions related to path enumeration in graphs.

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :func:`all_paths` | Return the list of all paths between a pair of vertices.
    :func:`yen_k_shortest_simple_paths` | Return an iterator over the simple paths between a pair of vertices in increasing order of weights.
    :func:`feng_k_shortest_simple_paths` | Return an iterator over the simple paths between a pair of vertices in increasing order of weights.
    :func:`all_paths_iterator` | Return an iterator over the paths of ``self``.
    :func:`all_simple_paths` | Return a list of all the simple paths of ``self`` starting with one of the given vertices.
    :func:`shortest_simple_paths` | Return an iterator over the simple paths between a pair of vertices.

Functions
---------
"""
# ****************************************************************************
# Copyright (C) 2019 Rajat Mittal <rajat.mttl@gmail.com>
#                    David Coudert <david.coudert@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from itertools import product

from sage.misc.misc_c import prod
from libcpp.queue cimport priority_queue
from libcpp.pair cimport pair
from sage.rings.integer_ring import ZZ
import copy


def all_paths(G, start, end, use_multiedges=False, report_edges=False, labels=False):
    """
    Return the list of all paths between a pair of vertices.

    If ``start`` is the same vertex as ``end``, then ``[[start]]`` is returned
    -- a list containing the 1-vertex, 0-edge path "``start``".

    If ``G`` has multiple edges, a path will be returned as many times as the
    product of the multiplicity of the edges along that path depending on the
    value of the flag ``use_multiedges``.

    INPUT:

    - ``start`` -- a vertex of a graph, where to start

    - ``end`` -- a vertex of a graph, where to end

    - ``use_multiedges`` -- boolean (default: ``False``); this parameter is
      used only if the graph has multiple edges

        - If ``False``, the graph is considered as simple and an edge label
          is arbitrarily selected for each edge as in
          :meth:`sage.graphs.generic_graph.GenericGraph.to_simple` if
          ``report_edges`` is ``True``

        - If ``True``, a path will be reported as many times as the edges
          multiplicities along that path (when ``report_edges = False`` or
          ``labels = False``), or with all possible combinations of edge
          labels (when ``report_edges = True`` and ``labels = True``)

    - ``report_edges`` -- boolean (default: ``False``); whether to report
      paths as list of vertices (default) or list of edges, if ``False``
      then ``labels`` parameter is ignored

    - ``labels`` -- boolean (default: ``False``); if ``False``, each edge
      is simply a pair ``(u, v)`` of vertices. Otherwise a list of edges
      along with its edge labels are used to represent the path.

    EXAMPLES::

        sage: eg1 = Graph({0:[1, 2], 1:[4], 2:[3, 4], 4:[5], 5:[6]})
        sage: eg1.all_paths(0, 6)
        [[0, 1, 4, 5, 6], [0, 2, 4, 5, 6]]
        sage: eg2 = graphs.PetersenGraph()
        sage: sorted(eg2.all_paths(1, 4))
        [[1, 0, 4],
         [1, 0, 5, 7, 2, 3, 4],
         [1, 0, 5, 7, 2, 3, 8, 6, 9, 4],
         [1, 0, 5, 7, 9, 4],
         [1, 0, 5, 7, 9, 6, 8, 3, 4],
         [1, 0, 5, 8, 3, 2, 7, 9, 4],
         [1, 0, 5, 8, 3, 4],
         [1, 0, 5, 8, 6, 9, 4],
         [1, 0, 5, 8, 6, 9, 7, 2, 3, 4],
         [1, 2, 3, 4],
         [1, 2, 3, 8, 5, 0, 4],
         [1, 2, 3, 8, 5, 7, 9, 4],
         [1, 2, 3, 8, 6, 9, 4],
         [1, 2, 3, 8, 6, 9, 7, 5, 0, 4],
         [1, 2, 7, 5, 0, 4],
         [1, 2, 7, 5, 8, 3, 4],
         [1, 2, 7, 5, 8, 6, 9, 4],
         [1, 2, 7, 9, 4],
         [1, 2, 7, 9, 6, 8, 3, 4],
         [1, 2, 7, 9, 6, 8, 5, 0, 4],
         [1, 6, 8, 3, 2, 7, 5, 0, 4],
         [1, 6, 8, 3, 2, 7, 9, 4],
         [1, 6, 8, 3, 4],
         [1, 6, 8, 5, 0, 4],
         [1, 6, 8, 5, 7, 2, 3, 4],
         [1, 6, 8, 5, 7, 9, 4],
         [1, 6, 9, 4],
         [1, 6, 9, 7, 2, 3, 4],
         [1, 6, 9, 7, 2, 3, 8, 5, 0, 4],
         [1, 6, 9, 7, 5, 0, 4],
         [1, 6, 9, 7, 5, 8, 3, 4]]
        sage: dg = DiGraph({0:[1, 3], 1:[3], 2:[0, 3]})
        sage: sorted(dg.all_paths(0, 3))
        [[0, 1, 3], [0, 3]]
        sage: ug = dg.to_undirected()
        sage: sorted(ug.all_paths(0, 3))
        [[0, 1, 3], [0, 2, 3], [0, 3]]

        sage: g = Graph([(0, 1), (0, 1), (1, 2), (1, 2)], multiedges=True)
        sage: g.all_paths(0, 2, use_multiedges=True)
        [[0, 1, 2], [0, 1, 2], [0, 1, 2], [0, 1, 2]]

        sage: dg = DiGraph({0:[1, 2, 1], 3:[0, 0]}, multiedges=True)
        sage: dg.all_paths(3, 1, use_multiedges=True)
        [[3, 0, 1], [3, 0, 1], [3, 0, 1], [3, 0, 1]]

        sage: g = Graph([(0, 1, 'a'), (0, 1, 'b'), (1, 2, 'c'), (1, 2, 'd')], multiedges=True)
        sage: g.all_paths(0, 2, use_multiedges=False)
        [[0, 1, 2]]
        sage: g.all_paths(0, 2, use_multiedges=True)
        [[0, 1, 2], [0, 1, 2], [0, 1, 2], [0, 1, 2]]
        sage: g.all_paths(0, 2, use_multiedges=True, report_edges=True)
        [[(0, 1), (1, 2)], [(0, 1), (1, 2)], [(0, 1), (1, 2)], [(0, 1), (1, 2)]]
        sage: g.all_paths(0, 2, use_multiedges=True, report_edges=True, labels=True)
        [((0, 1, 'b'), (1, 2, 'd')),
         ((0, 1, 'b'), (1, 2, 'c')),
         ((0, 1, 'a'), (1, 2, 'd')),
         ((0, 1, 'a'), (1, 2, 'c'))]
        sage: g.all_paths(0, 2, use_multiedges=False, report_edges=True, labels=True)
        [((0, 1, 'b'), (1, 2, 'd'))]
        sage: g.all_paths(0, 2, use_multiedges=False, report_edges=False, labels=True)
        [[0, 1, 2]]
        sage: g.all_paths(0, 2, use_multiedges=True, report_edges=False, labels=True)
        [[0, 1, 2], [0, 1, 2], [0, 1, 2], [0, 1, 2]]

    TESTS:

    Starting and ending at the same vertex (see :issue:`13006`)::

        sage: graphs.CompleteGraph(4).all_paths(2, 2)
        [[2]]

    Non-existing vertex as end vertex (see :issue:`24495`)::

        sage: g = graphs.PathGraph(5)
        sage: g.all_paths(1, 'junk')
        Traceback (most recent call last):
        ...
        LookupError: end vertex (junk) is not a vertex of the graph

    Distinguishing between multiedged paths (see :issue:`27501`)::

        sage: g = Graph(multiedges=True)
        sage: g.add_edge(0, 3, 1)
        sage: g.add_edge(0, 2, 3)
        sage: g.add_edge(0, 1, 3)
        sage: g.add_edge(2, 3, 5)
        sage: g.add_edge(2, 3, 15)
        sage: g.add_edge(2, 4, 12)
        sage: g.add_edge(3, 5, 7)
        sage: g.all_paths(0, 5, use_multiedges=True)
        [[0, 2, 3, 5], [0, 2, 3, 5], [0, 3, 5]]

        sage: g = Graph(multiedges=True)
        sage: g.add_edge(0, 1, 1)
        sage: g.add_edge(0, 2, 3)
        sage: g.add_edge(1, 4, 3)
        sage: g.add_edge(2, 3, 5)
        sage: g.add_edge(2, 4, 15)
        sage: g.add_edge(2, 4, 12)
        sage: g.add_edge(4, 5, 7)
        sage: g.add_edge(4, 5, 8)
        sage: g.add_edge(5, 6, 2)
        sage: g.all_paths(0, 6, use_multiedges=True)
        [[0, 1, 4, 5, 6],
         [0, 1, 4, 5, 6],
         [0, 2, 4, 5, 6],
         [0, 2, 4, 5, 6],
         [0, 2, 4, 5, 6],
         [0, 2, 4, 5, 6]]

    Added reporting of edges (see :issue:`27501`)::

        sage: G = DiGraph(multiedges=True)
        sage: G.add_edges([(0, 2), (0, 3), (0, 4), (1, 2), (1, 2), (1, 5), (3, 5), (3, 5)])
        sage: G.all_paths(0, 5, report_edges=True)
        [[(0, 3), (3, 5)]]
        sage: G.all_paths(0, 5, report_edges=True, use_multiedges=True)
        [[(0, 3), (3, 5)], [(0, 3), (3, 5)]]
    """
    if start not in G:
        raise LookupError("start vertex ({0}) is not a vertex of the graph".format(start))
    if end not in G:
        raise LookupError("end vertex ({0}) is not a vertex of the graph".format(end))

    if G.is_directed():
        iterator = G.neighbor_out_iterator
    else:
        iterator = G.neighbor_iterator

    if report_edges and labels:
        edge_labels = {}
        if use_multiedges:
            for e in G.edge_iterator():
                if (e[0], e[1]) in edge_labels:
                    edge_labels[(e[0], e[1])].append(e)
                else:
                    edge_labels[(e[0], e[1])] = [e]
        else:
            for e in G.edge_iterator():
                if (e[0], e[1]) not in edge_labels:
                    edge_labels[(e[0], e[1])] = [e]
        if not G.is_directed():
            for u, v in list(edge_labels):
                edge_labels[v, u] = edge_labels[u, v]
    elif use_multiedges and G.has_multiple_edges():
        from collections import Counter
        edge_multiplicity = Counter(G.edge_iterator(labels=False))

    if start == end:
        return [[start]]

    all_paths = []      # list of
    act_path = []       # the current path
    act_path_iter = []  # the neighbor/successor-iterators of the current path
    done = False
    s = start
    while not done:
        if s == end:    # if path completes, add to list
            all_paths.append(act_path + [s])
        else:
            if s not in act_path:   # we want vertices just once in a path
                act_path.append(s)  # extend current path
                act_path_iter.append(iterator(s))  # save the state of the neighbor/successor-iterator of the current vertex
        s = None
        while (s is None) and not done:
            try:
                s = next(act_path_iter[-1])  # try to get the next neighbor/successor, ...
            except (StopIteration):          # ... if there is none ...
                act_path.pop()               # ... go one step back
                act_path_iter.pop()
            if not act_path:                 # there is no other vertex ...
                done = True                  # ... so we are done

    if report_edges and labels:
        path_with_labels = []
        for p in all_paths:
            path_with_labels.extend(product(*[edge_labels[e] for e in zip(p[:-1], p[1:])]))
        return path_with_labels
    elif use_multiedges and G.has_multiple_edges():
        multiple_all_paths = []
        for p in all_paths:
            m = prod(edge_multiplicity[e] for e in zip(p[:-1], p[1:]))
            if report_edges:
                ep = list(zip(p[:-1], p[1:]))
            for _ in range(m):
                if report_edges:
                    multiple_all_paths.append(ep)
                else:
                    multiple_all_paths.append(p)
        return multiple_all_paths
    elif report_edges:
        return [list(zip(p[:-1], p[1:])) for p in all_paths]
    return all_paths


def shortest_simple_paths(self, source, target, weight_function=None,
                          by_weight=False, check_weight=True,
                          algorithm=None, report_edges=False,
                          labels=False, report_weight=False):
    r"""
    Return an iterator over the simple paths between a pair of vertices.

    This method returns an iterator over the simple paths (i.e., without
    repetition) from ``source`` to ``target``. By default (``by_weight`` is
    ``False``), the paths are reported by increasing number of edges. When
    ``by_weight`` is ``True``, the paths are reported by increasing weights.

    In case of weighted graphs negative weights are not allowed.

    If ``source`` is the same vertex as ``target``, then ``[[source]]`` is
    returned -- a list containing the 1-vertex, 0-edge path ``source``.

    By default ``Yen's`` algorithm [Yen1970]_ is used for undirected graphs and
    ``Feng's`` algorithm is used for directed graphs [Feng2014]_.

    The loops and the multiedges if present in the given graph are ignored and
    only minimum of the edge labels is kept in case of multiedges.

    INPUT:

    - ``source`` -- a vertex of the graph, where to start

    - ``target`` -- a vertex of the graph, where to end

    - ``weight_function`` -- function (default: ``None``); a function that
      takes as input an edge ``(u, v, l)`` and outputs its weight. If not
      ``None``, ``by_weight`` is automatically set to ``True``. If ``None``
      and ``by_weight`` is ``True``, we use the edge label ``l`` as a
      weight.

    - ``by_weight`` -- boolean (default: ``False``); if ``True``, the edges
      in the graph are weighted, otherwise all edges have weight 1

    - ``check_weight`` -- boolean (default: ``True``); whether to check that the
      ``weight_function`` outputs a number for each edge

    - ``algorithm`` -- string (default: ``None``); the algorithm to use in
      computing ``k`` shortest paths of ``self``. The following algorithms are
      supported:

      - ``'Yen'`` -- Yen's algorithm [Yen1970]_

      - ``'Feng'`` -- an improved version of Yen's algorithm but that works only
        for directed graphs [Feng2014]_

    - ``report_edges`` -- boolean (default: ``False``); whether to report paths
      as list of vertices (default) or list of edges. When set to ``False``, the
      ``labels`` parameter is ignored.

    - ``labels`` -- boolean (default: ``False``); if ``False``, each edge is
      simply a pair ``(u, v)`` of vertices. Otherwise a list of edges along
      with its edge labels are used to represent the path.

    - ``report_weight`` -- boolean (default: ``False``); if ``False``, just the
      path between ``source`` and ``target`` is returned. Otherwise a tuple of
      path length and path is returned.

    EXAMPLES::

        sage: g = DiGraph([(1, 2, 20), (1, 3, 10), (1, 4, 30),
        ....:              (2, 5, 20), (3, 5, 10), (4, 5, 30)])
        sage: list(g.shortest_simple_paths(1, 5, by_weight=True, algorithm='Yen'))
        [[1, 3, 5], [1, 2, 5], [1, 4, 5]]
        sage: list(g.shortest_simple_paths(1, 5, algorithm='Yen'))
        [[1, 2, 5], [1, 3, 5], [1, 4, 5]]
        sage: list(g.shortest_simple_paths(1, 1))
        [[1]]
        sage: list(g.shortest_simple_paths(1, 5, by_weight=True,
        ....:                              report_edges=True, report_weight=True, labels=True))
        [(20, [(1, 3, 10), (3, 5, 10)]),
         (40, [(1, 2, 20), (2, 5, 20)]),
         (60, [(1, 4, 30), (4, 5, 30)])]
        sage: list(g.shortest_simple_paths(1, 5, by_weight=True, algorithm='Feng',
        ....:                              report_edges=True, report_weight=True))
        [(20, [(1, 3), (3, 5)]), (40, [(1, 2), (2, 5)]), (60, [(1, 4), (4, 5)])]
        sage: list(g.shortest_simple_paths(1, 5, report_edges=True, report_weight=True))
        [(2, [(1, 2), (2, 5)]), (2, [(1, 3), (3, 5)]), (2, [(1, 4), (4, 5)])]
        sage: list(g.shortest_simple_paths(1, 5, by_weight=True, report_edges=True))
        [[(1, 3), (3, 5)], [(1, 2), (2, 5)], [(1, 4), (4, 5)]]
        sage: list(g.shortest_simple_paths(1, 5, by_weight=True, algorithm='Feng',
        ....:                              report_edges=True, labels=True))
        [[(1, 3, 10), (3, 5, 10)], [(1, 2, 20), (2, 5, 20)], [(1, 4, 30), (4, 5, 30)]]
        sage: g = Graph([(1, 2, 20), (1, 3, 10), (1, 4, 30), (2, 5, 20),
        ....:            (3, 5, 10), (4, 5, 30), (1, 6, 100), (5, 6, 5)])
        sage: list(g.shortest_simple_paths(1, 6, by_weight = True))
        [[1, 3, 5, 6], [1, 2, 5, 6], [1, 4, 5, 6], [1, 6]]
        sage: list(g.shortest_simple_paths(1, 6, algorithm='Yen'))
        [[1, 6], [1, 2, 5, 6], [1, 3, 5, 6], [1, 4, 5, 6]]
        sage: list(g.shortest_simple_paths(1, 6,
        ....:                              report_edges=True, report_weight=True, labels=True))
        [(1, [(1, 6, 100)]),
         (3, [(1, 2, 20), (2, 5, 20), (5, 6, 5)]),
         (3, [(1, 3, 10), (3, 5, 10), (5, 6, 5)]),
         (3, [(1, 4, 30), (4, 5, 30), (5, 6, 5)])]
        sage: list(g.shortest_simple_paths(1, 6, by_weight=True,
        ....:                              report_edges=True, report_weight=True, labels=True))
        [(25, [(1, 3, 10), (3, 5, 10), (5, 6, 5)]),
         (45, [(1, 2, 20), (2, 5, 20), (5, 6, 5)]),
         (65, [(1, 4, 30), (4, 5, 30), (5, 6, 5)]),
         (100, [(1, 6, 100)])]
        sage: list(g.shortest_simple_paths(1, 6, by_weight=True,
        ....:                              report_edges=True, labels=True))
        [[(1, 3, 10), (3, 5, 10), (5, 6, 5)],
         [(1, 2, 20), (2, 5, 20), (5, 6, 5)],
         [(1, 4, 30), (4, 5, 30), (5, 6, 5)],
         [(1, 6, 100)]]
        sage: list(g.shortest_simple_paths(1, 6, report_edges=True, labels=True))
        [[(1, 6, 100)],
         [(1, 2, 20), (2, 5, 20), (5, 6, 5)],
         [(1, 3, 10), (3, 5, 10), (5, 6, 5)],
         [(1, 4, 30), (4, 5, 30), (5, 6, 5)]]

    TESTS::

        sage: g = Graph([(1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 2),
        ....:            (5, 6, 100), (4, 7, 3), (7, 6, 4), (3, 8, 5),
        ....:            (8, 9, 2), (9, 6, 2), (9, 10, 7), (9, 11, 10),
        ....:            (11, 6, 8), (10, 6, 2)])
        sage: list(g.shortest_simple_paths(1, 6, algorithm='Yen', by_weight=True))
        [[1, 2, 3, 4, 7, 6],
         [1, 2, 3, 8, 9, 6],
         [1, 2, 3, 8, 9, 10, 6],
         [1, 2, 3, 8, 9, 11, 6],
         [1, 2, 3, 4, 5, 6]]
        sage: list(g.shortest_simple_paths(1, 6, report_edges=True, labels=True, by_weight=True))
        [[(1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 7, 3), (6, 7, 4)],
         [(1, 2, 1), (2, 3, 1), (3, 8, 5), (8, 9, 2), (6, 9, 2)],
         [(1, 2, 1), (2, 3, 1), (3, 8, 5), (8, 9, 2), (9, 10, 7), (6, 10, 2)],
         [(1, 2, 1), (2, 3, 1), (3, 8, 5), (8, 9, 2), (9, 11, 10), (6, 11, 8)],
         [(1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 2), (5, 6, 100)]]
        sage: list(g.shortest_simple_paths(1, 6, report_edges=True, labels=True, by_weight=True, report_weight=True))
        [(10, [(1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 7, 3), (6, 7, 4)]),
         (11, [(1, 2, 1), (2, 3, 1), (3, 8, 5), (8, 9, 2), (6, 9, 2)]),
         (18, [(1, 2, 1), (2, 3, 1), (3, 8, 5), (8, 9, 2), (9, 10, 7), (6, 10, 2)]),
         (27, [(1, 2, 1), (2, 3, 1), (3, 8, 5), (8, 9, 2), (9, 11, 10), (6, 11, 8)]),
         (105, [(1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 2), (5, 6, 100)])]
        sage: list(g.shortest_simple_paths(1, 6, algorithm='Yen'))
        [[1, 2, 3, 4, 5, 6],
         [1, 2, 3, 4, 7, 6],
         [1, 2, 3, 8, 9, 6],
         [1, 2, 3, 8, 9, 10, 6],
         [1, 2, 3, 8, 9, 11, 6]]
        sage: g = DiGraph([(1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1),
        ....:              (1, 7, 1), (7, 8, 1), (8, 5, 1), (1, 6, 1),
        ....:              (6, 9, 1), (9, 5, 1), (4, 2, 1), (9, 3, 1),
        ....:              (9, 10, 1), (10, 5, 1), (9, 11, 1), (11, 10, 1)])
        sage: list(g.shortest_simple_paths(1, 5, algorithm='Feng'))
        [[1, 6, 9, 5],
         [1, 7, 8, 5],
         [1, 2, 3, 4, 5],
         [1, 6, 9, 10, 5],
         [1, 6, 9, 11, 10, 5],
         [1, 6, 9, 3, 4, 5]]

        sage: # needs sage.combinat
        sage: G = digraphs.DeBruijn(2, 3)
        sage: for u,v in G.edges(sort=True, labels=False):
        ....:     G.set_edge_label(u, v, 1)
        sage: G.allow_multiple_edges(True)
        sage: for u,v in G.edges(sort=True, labels=False):
        ....:     G.add_edge(u, v, 2)
        sage: list(G.shortest_simple_paths('000', '111'))
        [['000', '001', '011', '111'], ['000', '001', '010', '101', '011', '111']]
        sage: list(G.shortest_simple_paths('000', '111', by_weight=True))
        [['000', '001', '011', '111'], ['000', '001', '010', '101', '011', '111']]
        sage: list(G.shortest_simple_paths('000', '111', by_weight=True, report_weight=True))
        [(3, ['000', '001', '011', '111']),
         (5, ['000', '001', '010', '101', '011', '111'])]
        sage: list(G.shortest_simple_paths('000', '111', by_weight=True, report_weight=True, report_edges=True, labels=True))
        [(3, [('000', '001', 1), ('001', '011', 1), ('011', '111', 1)]),
         (5,
          [('000', '001', 1),
           ('001', '010', 1),
           ('010', '101', 1),
           ('101', '011', 1),
           ('011', '111', 1)])]

    Feng's algorithm cannot be used on undirected graphs::

        sage: list(graphs.PathGraph(2).shortest_simple_paths(0, 1, algorithm='Feng'))
        Traceback (most recent call last):
        ...
        ValueError: Feng's algorithm works only for directed graphs

    If the algorithm is not implemented::

        sage: list(g.shortest_simple_paths(1, 5, algorithm='tip top'))
        Traceback (most recent call last):
        ...
        ValueError: unknown algorithm "tip top"

    Check for consistency of results of Yen's and Feng's::

        sage: # needs sage.combinat
        sage: G = digraphs.DeBruijn(2, 4)
        sage: s = set()
        sage: for p in G.shortest_simple_paths('0000', '1111', by_weight=False, algorithm='Yen'):
        ....:     s.add(tuple(p))
        sage: k = set()
        sage: for p in G.shortest_simple_paths('0000', '1111', by_weight=False, algorithm='Feng'):
        ....:     k.add(tuple(p))
        sage: k == s
        True

        sage: G = DiGraph(graphs.Grid2dGraph(3, 3))
        sage: s = set()
        sage: for i, p in enumerate(G.shortest_simple_paths((0, 0), (0, 1), by_weight=False, algorithm='Feng')):
        ....:     s.add(tuple(p))
        sage: k = set()
        sage: for i, p in enumerate(G.shortest_simple_paths((0, 0), (0, 1), by_weight=False, algorithm='Yen')):
        ....:     k.add(tuple(p))
        sage: s == k
        True

        sage: G = DiGraph('SL{Sa??B[??iSOBIgA_K?a?@H??aGCsc??_oGCC__AA?H????c@_GA?C@?A_?_C???a?')
        sage: s = set()
        sage: for i, p in enumerate(G.shortest_simple_paths(0, 1, by_weight=False, algorithm='Yen')):
        ....:     s.add(tuple(p))
        sage: t = set()
        sage: for i, p in enumerate(G.shortest_simple_paths(0, 1, by_weight=False, algorithm='Feng')):
        ....:     t.add(tuple(p))
        sage: s == t
        True

        sage: G = digraphs.Circulant(10, [2, 3])
        sage: s = set()
        sage: for i, p in enumerate(G.shortest_simple_paths(1, 7, by_weight=False, algorithm='Yen')):
        ....:     s.add(tuple(p))
        sage: t = set()
        sage: for i, p in enumerate(G.shortest_simple_paths(1, 7, by_weight=False, algorithm='Feng')):
        ....:     t.add(tuple(p))
        sage: s == t
        True

    Check that "Yen" and "Feng" provide same results on random digraphs::

        sage: G = digraphs.RandomDirectedGNP(30, .05)
        sage: while not G.is_strongly_connected():
        ....:     G = digraphs.RandomDirectedGNP(30, .1)
        sage: for u, v in list(G.edges(labels=False, sort=False)):
        ....:     G.set_edge_label(u, v, randint(1, 10))
        sage: V = G.vertices(sort=False)
        sage: shuffle(V)
        sage: u, v = V[:2]
        sage: it_Y = G.shortest_simple_paths(u, v, by_weight=True, report_weight=True, algorithm='Yen')
        sage: it_F = G.shortest_simple_paths(u, v, by_weight=True, report_weight=True, algorithm='Feng')
        sage: for i, (y, f) in enumerate(zip(it_Y, it_F)):
        ....:     if y[0] != f[0]:
        ....:         raise ValueError("something goes wrong !")
        ....:     if i == 100:
        ....:         break
    """
    if source not in self:
        raise ValueError("vertex '{}' is not in the graph".format(source))

    if target not in self:
        raise ValueError("vertex '{}' is not in the graph".format(target))

    if source == target:
        if report_edges:
            yield []
        elif report_weight:
            yield (0, [source])
        else:
            yield [source]
        return

    if self.has_loops() or self.allows_multiple_edges():
        self = self.to_simple(to_undirected=False, keep_label='min', immutable=False)

    if algorithm is None:
        algorithm = "Feng" if self.is_directed() else "Yen"

    if algorithm == "Feng":
        if not self.is_directed():
            raise ValueError("Feng's algorithm works only for directed graphs")

        yield from feng_k_shortest_simple_paths(self, source=source, target=target,
                                                weight_function=weight_function,
                                                by_weight=by_weight, check_weight=check_weight,
                                                report_edges=report_edges,
                                                labels=labels, report_weight=report_weight)

    elif algorithm == "Yen":
        yield from yen_k_shortest_simple_paths(self, source=source, target=target,
                                               weight_function=weight_function,
                                               by_weight=by_weight, check_weight=check_weight,
                                               report_edges=report_edges,
                                               labels=labels, report_weight=report_weight)

    else:
        raise ValueError('unknown algorithm "{}"'.format(algorithm))


def yen_k_shortest_simple_paths(self, source, target, weight_function=None,
                                by_weight=False, check_weight=True,
                                report_edges=False,
                                labels=False, report_weight=False):
    r"""
    Return an iterator over the simple paths between a pair of vertices in
    increasing order of weights.

    For unweighted graphs paths are returned in order of increasing number
    of edges.

    In case of weighted graphs negative weights are not allowed.

    If ``source`` is the same vertex as ``target``, then ``[[source]]`` is
    returned -- a list containing the 1-vertex, 0-edge path ``source``.

    The loops and the multiedges if present in the given graph are ignored and
    only minimum of the edge labels is kept in case of multiedges.

    INPUT:

    - ``source`` -- a vertex of the graph, where to start

    - ``target`` -- a vertex of the graph, where to end

    - ``weight_function`` -- function (default: ``None``); a function that
      takes as input an edge ``(u, v, l)`` and outputs its weight. If not
      ``None``, ``by_weight`` is automatically set to ``True``. If ``None``
      and ``by_weight`` is ``True``, we use the edge label ``l`` as a
      weight.

    - ``by_weight`` -- boolean (default: ``False``); if ``True``, the edges
      in the graph are weighted, otherwise all edges have weight 1

    - ``check_weight`` -- boolean (default: ``True``); whether to check that
      the ``weight_function`` outputs a number for each edge

    - ``report_edges`` -- boolean (default: ``False``); whether to report
      paths as list of vertices (default) or list of edges, if ``False``
      then ``labels`` parameter is ignored

    - ``labels`` -- boolean (default: ``False``); if ``False``, each edge
      is simply a pair ``(u, v)`` of vertices. Otherwise a list of edges
      along with its edge labels are used to represent the path.

    - ``report_weight`` -- boolean (default: ``False``); if ``False``, just
      the path between ``source`` and ``target`` is returned. Otherwise a
      tuple of path length and path is returned.

    ALGORITHM:

    This algorithm can be divided into two parts. Firstly, it determines a
    shortest path from ``source`` to ``target``. Then, it determines all the
    other `k`-shortest paths.  This algorithm finds the deviations of previous
    shortest paths to determine the next shortest paths.

    Time complexity is `O(kn(m+n\log{n}))` where `n` is the number of vertices
    and `m` is the number of edges and `k` is the number of shortest paths
    needed to find.

    See [Yen1970]_ and the :wikipedia:`Yen%27s_algorithm` for more details on the
    algorithm.

    EXAMPLES::

        sage: from sage.graphs.path_enumeration import yen_k_shortest_simple_paths
        sage: g = DiGraph([(1, 2, 20), (1, 3, 10), (1, 4, 30), (2, 5, 20), (3, 5, 10), (4, 5, 30)])
        sage: list(yen_k_shortest_simple_paths(g, 1, 5, by_weight=True))
        [[1, 3, 5], [1, 2, 5], [1, 4, 5]]
        sage: list(yen_k_shortest_simple_paths(g, 1, 5))
        [[1, 2, 5], [1, 3, 5], [1, 4, 5]]
        sage: list(yen_k_shortest_simple_paths(g, 1, 1))
        [[1]]
        sage: list(yen_k_shortest_simple_paths(g, 1, 5, by_weight=True, report_edges=True, report_weight=True, labels=True))
        [(20, [(1, 3, 10), (3, 5, 10)]),
         (40, [(1, 2, 20), (2, 5, 20)]),
         (60, [(1, 4, 30), (4, 5, 30)])]
        sage: list(yen_k_shortest_simple_paths(g, 1, 5, by_weight=True, report_edges=True, report_weight=True))
        [(20, [(1, 3), (3, 5)]), (40, [(1, 2), (2, 5)]), (60, [(1, 4), (4, 5)])]
        sage: list(yen_k_shortest_simple_paths(g, 1, 5, report_edges=True, report_weight=True))
        [(2, [(1, 2), (2, 5)]), (2, [(1, 3), (3, 5)]), (2, [(1, 4), (4, 5)])]
        sage: list(yen_k_shortest_simple_paths(g, 1, 5, by_weight=True, report_edges=True))
        [[(1, 3), (3, 5)], [(1, 2), (2, 5)], [(1, 4), (4, 5)]]
        sage: list(yen_k_shortest_simple_paths(g, 1, 5, by_weight=True, report_edges=True, labels=True))
        [[(1, 3, 10), (3, 5, 10)], [(1, 2, 20), (2, 5, 20)], [(1, 4, 30), (4, 5, 30)]]
        sage: from sage.graphs.path_enumeration import yen_k_shortest_simple_paths
        sage: g = Graph([(1, 2, 20), (1, 3, 10), (1, 4, 30), (2, 5, 20), (3, 5, 10), (4, 5, 30), (1, 6, 100), (5, 6, 5)])
        sage: list(yen_k_shortest_simple_paths(g, 1, 6, by_weight = True))
        [[1, 3, 5, 6], [1, 2, 5, 6], [1, 4, 5, 6], [1, 6]]
        sage: list(yen_k_shortest_simple_paths(g, 1, 6))
        [[1, 6], [1, 2, 5, 6], [1, 3, 5, 6], [1, 4, 5, 6]]
        sage: list(yen_k_shortest_simple_paths(g, 1, 6, report_edges=True, report_weight=True, labels=True))
        [(1, [(1, 6, 100)]),
         (3, [(1, 2, 20), (2, 5, 20), (5, 6, 5)]),
         (3, [(1, 3, 10), (3, 5, 10), (5, 6, 5)]),
         (3, [(1, 4, 30), (4, 5, 30), (5, 6, 5)])]
        sage: list(yen_k_shortest_simple_paths(g, 1, 6, report_edges=True, report_weight=True, labels=True, by_weight=True))
        [(25, [(1, 3, 10), (3, 5, 10), (5, 6, 5)]),
         (45, [(1, 2, 20), (2, 5, 20), (5, 6, 5)]),
         (65, [(1, 4, 30), (4, 5, 30), (5, 6, 5)]),
         (100, [(1, 6, 100)])]
        sage: list(yen_k_shortest_simple_paths(g, 1, 6, report_edges=True, labels=True, by_weight=True))
        [[(1, 3, 10), (3, 5, 10), (5, 6, 5)],
         [(1, 2, 20), (2, 5, 20), (5, 6, 5)],
         [(1, 4, 30), (4, 5, 30), (5, 6, 5)],
         [(1, 6, 100)]]
        sage: list(yen_k_shortest_simple_paths(g, 1, 6, report_edges=True, labels=True))
        [[(1, 6, 100)],
         [(1, 2, 20), (2, 5, 20), (5, 6, 5)],
         [(1, 3, 10), (3, 5, 10), (5, 6, 5)],
         [(1, 4, 30), (4, 5, 30), (5, 6, 5)]]

    TESTS::

        sage: from sage.graphs.path_enumeration import yen_k_shortest_simple_paths
        sage: g = Graph([(1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 2),
        ....:            (5, 6, 100), (4, 7, 3), (7, 6, 4), (3, 8, 5),
        ....:            (8, 9, 2), (9, 6, 2), (9, 10, 7), (9, 11, 10),
        ....:            (11, 6, 8), (10, 6, 2)])
        sage: list(yen_k_shortest_simple_paths(g, 1, 6, by_weight=True))
        [[1, 2, 3, 4, 7, 6],
         [1, 2, 3, 8, 9, 6],
         [1, 2, 3, 8, 9, 10, 6],
         [1, 2, 3, 8, 9, 11, 6],
         [1, 2, 3, 4, 5, 6]]
        sage: list(yen_k_shortest_simple_paths(g, 1, 6, report_edges=True, labels=True, by_weight=True))
        [[(1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 7, 3), (6, 7, 4)],
         [(1, 2, 1), (2, 3, 1), (3, 8, 5), (8, 9, 2), (6, 9, 2)],
         [(1, 2, 1), (2, 3, 1), (3, 8, 5), (8, 9, 2), (9, 10, 7), (6, 10, 2)],
         [(1, 2, 1), (2, 3, 1), (3, 8, 5), (8, 9, 2), (9, 11, 10), (6, 11, 8)],
         [(1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 2), (5, 6, 100)]]
        sage: list(yen_k_shortest_simple_paths(g, 1, 6, report_edges=True, labels=True, by_weight=True, report_weight=True))
        [(10, [(1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 7, 3), (6, 7, 4)]),
         (11, [(1, 2, 1), (2, 3, 1), (3, 8, 5), (8, 9, 2), (6, 9, 2)]),
         (18, [(1, 2, 1), (2, 3, 1), (3, 8, 5), (8, 9, 2), (9, 10, 7), (6, 10, 2)]),
         (27, [(1, 2, 1), (2, 3, 1), (3, 8, 5), (8, 9, 2), (9, 11, 10), (6, 11, 8)]),
         (105, [(1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 2), (5, 6, 100)])]
        sage: list(yen_k_shortest_simple_paths(g, 1, 6))
        [[1, 2, 3, 4, 5, 6],
         [1, 2, 3, 4, 7, 6],
         [1, 2, 3, 8, 9, 6],
         [1, 2, 3, 8, 9, 10, 6],
         [1, 2, 3, 8, 9, 11, 6]]
        sage: from sage.graphs.path_enumeration import yen_k_shortest_simple_paths
        sage: g = DiGraph([(1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1),
        ....:              (1, 7, 1), (7, 8, 1), (8, 5, 1), (1, 6, 1),
        ....:              (6, 9, 1), (9, 5, 1), (4, 2, 1), (9, 3, 1),
        ....:              (9, 10, 1), (10, 5, 1), (9, 11, 1), (11, 10, 1)])
        sage: list(yen_k_shortest_simple_paths(g, 1, 5))
        [[1, 6, 9, 5],
         [1, 7, 8, 5],
         [1, 2, 3, 4, 5],
         [1, 6, 9, 10, 5],
         [1, 6, 9, 3, 4, 5],
         [1, 6, 9, 11, 10, 5]]
    """
    if source not in self:
        raise ValueError("vertex '{}' is not in the graph".format(source))
    if target not in self:
        raise ValueError("vertex '{}' is not in the graph".format(target))

    if source == target:
        if report_edges:
            yield []
        elif report_weight:
            yield (0, [source])
        else:
            yield [source]
        return

    if self.has_loops() or self.allows_multiple_edges():
        G = self.to_simple(to_undirected=False, keep_label='min', immutable=False)
    else:
        G = self

    by_weight, weight_function = self._get_weight_function(by_weight=by_weight,
                                                           weight_function=weight_function,
                                                           check_weight=check_weight)

    cdef dict edge_wt
    if by_weight:
        # dictionary to get weight of the edges
        edge_wt = {(e[0], e[1]): weight_function(e) for e in G.edge_iterator()}
        if not G.is_directed():
            for u, v in G.edge_iterator(labels=False):
                edge_wt[v, u] = edge_wt[u, v]

        def length_func(path):
            return sum(edge_wt[e] for e in zip(path[:-1], path[1:]))
        # shortest path function for weighted graph
        shortest_path_func = G._backend.bidirectional_dijkstra_special
    else:
        def length_func(path):
            return len(path) - 1
        # shortest path function for unweighted graph
        shortest_path_func = G._backend.shortest_path_special

    # compute the shortest path between the source and the target
    cdef list path
    if by_weight:
        path = shortest_path_func(source, target, weight_function=weight_function)
    else:
        path = shortest_path_func(source, target)
    # corner case
    if not path:
        if report_weight:
            yield (0, [])
        else:
            yield []
        return

    cdef dict edge_labels
    if report_edges and labels:
        edge_labels = {(e[0], e[1]): e for e in G.edge_iterator()}
        if not G.is_directed():
            for u, v in G.edge_iterator(labels=False):
                edge_labels[v, u] = edge_labels[u, v]

    # heap data structure containing the candidate paths
    cdef priority_queue[pair[double, pair[int, int]]] heap_sorted_paths
    cdef int idx = 0
    heap_sorted_paths.push((-length_func(path), (idx, 0)))
    cdef dict idx_to_path = {idx: path}
    idx = idx + 1
    # list of all paths already yielded
    cdef list listA = list()

    cdef set exclude_vertices
    cdef set exclude_edges
    cdef list prev_path, new_path, root
    cdef int path_idx, dev_idx

    while idx_to_path:
        # extracting the next best path from the heap
        cost, (path_idx, dev_idx) = heap_sorted_paths.top()
        heap_sorted_paths.pop()
        prev_path = idx_to_path[path_idx]
        del idx_to_path[path_idx]
        if report_weight:
            cost = -cost
            if cost in ZZ:
                cost = int(cost)
            if report_edges and labels:
                yield (cost, [edge_labels[e] for e in zip(prev_path[:-1], prev_path[1:])])
            elif report_edges:
                yield (cost, list(zip(prev_path[:-1], prev_path[1:])))
            else:
                yield (cost, prev_path)
        else:
            if report_edges and labels:
                yield [edge_labels[e] for e in zip(prev_path[:-1], prev_path[1:])]
            elif report_edges:
                yield list(zip(prev_path[:-1], prev_path[1:]))
            else:
                yield prev_path

        listA.append(prev_path)
        exclude_vertices = set(prev_path[:dev_idx])
        exclude_edges = set()
        root = prev_path[:dev_idx]

        # deviating from the previous path to find the candidate paths
        for i in range(dev_idx + 1, len(prev_path)):
            # root part of the previous path
            root.append(prev_path[i - 1])
            for path in listA:
                if path[:i] == root:
                    exclude_edges.add((path[i - 1], path[i]))
                    if not G.is_directed():
                        exclude_edges.add((path[i], path[i - 1]))
            try:
                # finding the spur part of the path after excluding certain
                # vertices and edges
                if by_weight:
                    spur = shortest_path_func(root[-1], target,
                                              exclude_vertices=exclude_vertices,
                                              exclude_edges=exclude_edges,
                                              weight_function=weight_function)
                else:
                    spur = shortest_path_func(root[-1], target,
                                              exclude_vertices=exclude_vertices,
                                              exclude_edges=exclude_edges)
                if not spur:
                    continue
                # concatenating the root and the spur paths
                new_path = root[:-1] + spur
                # push operation
                idx_to_path[idx] = new_path
                heap_sorted_paths.push((-length_func(new_path), (idx, i - 1)))
                idx = idx + 1
            except Exception:
                pass
            exclude_vertices.add(root[-1])


def feng_k_shortest_simple_paths(self, source, target, weight_function=None,
                                 by_weight=False, check_weight=True,
                                 report_edges=False,
                                 labels=False, report_weight=False):
    r"""
    Return an iterator over the simple paths between a pair of vertices in
    increasing order of weights.

    Works only for directed graphs.

    For unweighted graphs, paths are returned in order of increasing number
    of edges.

    In case of weighted graphs, negative weights are not allowed.

    If ``source`` is the same vertex as ``target``, then ``[[source]]`` is
    returned -- a list containing the 1-vertex, 0-edge path ``source``.

    The loops and the multiedges if present in the given graph are ignored and
    only minimum of the edge labels is kept in case of multiedges.

    INPUT:

    - ``source`` -- a vertex of the graph, where to start

    - ``target`` -- a vertex of the graph, where to end

    - ``weight_function`` -- function (default: ``None``); a function that
      takes as input an edge ``(u, v, l)`` and outputs its weight. If not
      ``None``, ``by_weight`` is automatically set to ``True``. If ``None``
      and ``by_weight`` is ``True``, we use the edge label ``l`` as a
      weight.

    - ``by_weight`` -- boolean (default: ``False``); if ``True``, the edges
      in the graph are weighted, otherwise all edges have weight 1

    - ``check_weight`` -- boolean (default: ``True``); whether to check that
      the ``weight_function`` outputs a number for each edge

    - ``report_edges`` -- boolean (default: ``False``); whether to report
      paths as list of vertices (default) or list of edges, if ``False``
      then ``labels`` parameter is ignored

    - ``labels`` -- boolean (default: ``False``); if ``False``, each edge
      is simply a pair ``(u, v)`` of vertices. Otherwise a list of edges
      along with its edge labels are used to represent the path.

    - ``report_weight`` -- boolean (default: ``False``); if ``False``, just
      the path between ``source`` and ``target`` is returned. Otherwise a
      tuple of path length and path is returned.

    ALGORITHM:

    This algorithm can be divided into two parts. Firstly, it determines the
    shortest path from ``source`` to ``target``. Then, it determines all the
    other `k`-shortest paths. This algorithm finds the deviations of previous
    shortest paths to determine the next shortest paths. This algorithm finds
    the candidate paths more efficiently using a node classification
    technique. At first the candidate path is separated by its deviation node
    as prefix and suffix. Then the algorithm classify the nodes as red, yellow
    and green. A node on the prefix is assigned a red color, a node that can
    reach t (the destination node) through a shortest path without visiting a
    red node is assigned a green color, and all other nodes are assigned a
    yellow color. When searching for the suffix of a candidate path, all green
    nodes are bypassed, and ``Dijkstra’s algorithm`` is applied to find an
    all-yellow-node subpath.  Since on average the number of yellow nodes is
    much smaller than n, this algorithm has a much lower average-case running
    time.

    Time complexity is `O(kn(m+n\log{n}))` where `n` is the number of vertices
    and `m` is the number of edges and `k` is the number of shortest paths
    needed to find. Its average running time is much smaller as compared to
    `Yen's` algorithm.

    See [Feng2014]_ for more details on this algorithm.

    EXAMPLES::

        sage: from sage.graphs.path_enumeration import feng_k_shortest_simple_paths
        sage: g = DiGraph([(1, 2, 20), (1, 3, 10), (1, 4, 30), (2, 5, 20), (3, 5, 10), (4, 5, 30)])
        sage: list(feng_k_shortest_simple_paths(g, 1, 5, by_weight=True))
        [[1, 3, 5], [1, 2, 5], [1, 4, 5]]
        sage: list(feng_k_shortest_simple_paths(g, 1, 5))
        [[1, 2, 5], [1, 3, 5], [1, 4, 5]]
        sage: list(feng_k_shortest_simple_paths(g, 1, 1))
        [[1]]
        sage: list(feng_k_shortest_simple_paths(g, 1, 5, report_edges=True, labels=True))
        [[(1, 2, 20), (2, 5, 20)], [(1, 3, 10), (3, 5, 10)], [(1, 4, 30), (4, 5, 30)]]
        sage: list(feng_k_shortest_simple_paths(g, 1, 5, report_edges=True, labels=True, by_weight=True))
        [[(1, 3, 10), (3, 5, 10)], [(1, 2, 20), (2, 5, 20)], [(1, 4, 30), (4, 5, 30)]]
        sage: list(feng_k_shortest_simple_paths(g, 1, 5, report_edges=True, labels=True, by_weight=True, report_weight=True))
        [(20, [(1, 3, 10), (3, 5, 10)]),
         (40, [(1, 2, 20), (2, 5, 20)]),
         (60, [(1, 4, 30), (4, 5, 30)])]

        sage: from sage.graphs.path_enumeration import feng_k_shortest_simple_paths
        sage: g = DiGraph([(1, 2, 20), (1, 3, 10), (1, 4, 30), (2, 5, 20), (3, 5, 10), (4, 5, 30), (1, 6, 100), (5, 6, 5)])
        sage: list(feng_k_shortest_simple_paths(g, 1, 6, by_weight = True))
        [[1, 3, 5, 6], [1, 2, 5, 6], [1, 4, 5, 6], [1, 6]]
        sage: list(feng_k_shortest_simple_paths(g, 1, 6))
        [[1, 6], [1, 2, 5, 6], [1, 3, 5, 6], [1, 4, 5, 6]]
        sage: list(feng_k_shortest_simple_paths(g, 1, 6, report_edges=True, labels=True, by_weight=True, report_weight=True))
        [(25, [(1, 3, 10), (3, 5, 10), (5, 6, 5)]),
         (45, [(1, 2, 20), (2, 5, 20), (5, 6, 5)]),
         (65, [(1, 4, 30), (4, 5, 30), (5, 6, 5)]),
         (100, [(1, 6, 100)])]
        sage: list(feng_k_shortest_simple_paths(g, 1, 6, report_edges=True, labels=True, report_weight=True))
        [(1, [(1, 6, 100)]),
         (3, [(1, 2, 20), (2, 5, 20), (5, 6, 5)]),
         (3, [(1, 3, 10), (3, 5, 10), (5, 6, 5)]),
         (3, [(1, 4, 30), (4, 5, 30), (5, 6, 5)])]
        sage: from sage.graphs.path_enumeration import feng_k_shortest_simple_paths
        sage: g = DiGraph([(1, 2, 5), (2, 3, 0), (1, 4, 2), (4, 5, 1), (5, 3, 0)])
        sage: list(feng_k_shortest_simple_paths(g, 1, 3, by_weight=True))
        [[1, 4, 5, 3], [1, 2, 3]]
        sage: list(feng_k_shortest_simple_paths(g, 1, 3))
        [[1, 2, 3], [1, 4, 5, 3]]
        sage: list(feng_k_shortest_simple_paths(g, 1, 3, report_weight=True))
        [(2, [1, 2, 3]), (3, [1, 4, 5, 3])]
        sage: list(feng_k_shortest_simple_paths(g, 1, 3, report_weight=True, report_edges=True))
        [(2, [(1, 2), (2, 3)]), (3, [(1, 4), (4, 5), (5, 3)])]
        sage: list(feng_k_shortest_simple_paths(g, 1, 3, report_weight=True, report_edges=True, by_weight=True))
        [(3, [(1, 4), (4, 5), (5, 3)]), (5, [(1, 2), (2, 3)])]
        sage: list(feng_k_shortest_simple_paths(g, 1, 3, report_weight=True, report_edges=True, by_weight=True, labels=True))
        [(3, [(1, 4, 2), (4, 5, 1), (5, 3, 0)]), (5, [(1, 2, 5), (2, 3, 0)])]

    TESTS::

        sage: from sage.graphs.path_enumeration import feng_k_shortest_simple_paths
        sage: g = DiGraph([(1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 2), (5, 6, 100),
        ....:              (4, 7, 3), (7, 6, 4), (3, 8, 5), (8, 9, 2), (9, 6, 2),
        ....:              (9, 10, 7), (9, 11, 10), (11, 6, 8), (10, 6, 2)])
        sage: list(feng_k_shortest_simple_paths(g, 1, 6, by_weight=True))
        [[1, 2, 3, 4, 7, 6],
         [1, 2, 3, 8, 9, 6],
         [1, 2, 3, 8, 9, 10, 6],
         [1, 2, 3, 8, 9, 11, 6],
         [1, 2, 3, 4, 5, 6]]
        sage: list(feng_k_shortest_simple_paths(g, 1, 6, by_weight=True, report_edges=True))
        [[(1, 2), (2, 3), (3, 4), (4, 7), (7, 6)],
         [(1, 2), (2, 3), (3, 8), (8, 9), (9, 6)],
         [(1, 2), (2, 3), (3, 8), (8, 9), (9, 10), (10, 6)],
         [(1, 2), (2, 3), (3, 8), (8, 9), (9, 11), (11, 6)],
         [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)]]
        sage: list(feng_k_shortest_simple_paths(g, 1, 6, by_weight=True, report_edges=True, report_weight=True))
        [(10, [(1, 2), (2, 3), (3, 4), (4, 7), (7, 6)]),
         (11, [(1, 2), (2, 3), (3, 8), (8, 9), (9, 6)]),
         (18, [(1, 2), (2, 3), (3, 8), (8, 9), (9, 10), (10, 6)]),
         (27, [(1, 2), (2, 3), (3, 8), (8, 9), (9, 11), (11, 6)]),
         (105, [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)])]
        sage: list(feng_k_shortest_simple_paths(g, 1, 6, by_weight=True, report_edges=True, report_weight=True, labels=True))
        [(10, [(1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 7, 3), (7, 6, 4)]),
         (11, [(1, 2, 1), (2, 3, 1), (3, 8, 5), (8, 9, 2), (9, 6, 2)]),
         (18, [(1, 2, 1), (2, 3, 1), (3, 8, 5), (8, 9, 2), (9, 10, 7), (10, 6, 2)]),
         (27, [(1, 2, 1), (2, 3, 1), (3, 8, 5), (8, 9, 2), (9, 11, 10), (11, 6, 8)]),
         (105, [(1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 2), (5, 6, 100)])]
        sage: list(feng_k_shortest_simple_paths(g, 1, 6))
        [[1, 2, 3, 4, 5, 6],
         [1, 2, 3, 4, 7, 6],
         [1, 2, 3, 8, 9, 6],
         [1, 2, 3, 8, 9, 11, 6],
         [1, 2, 3, 8, 9, 10, 6]]
        sage: from sage.graphs.path_enumeration import feng_k_shortest_simple_paths
        sage: g = DiGraph([(1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1),
        ....:              (1, 7, 1), (7, 8, 1), (8, 5, 1), (1, 6, 1),
        ....:              (6, 9, 1), (9, 5, 1), (4, 2, 1), (9, 3, 1),
        ....:              (9, 10, 1), (10, 5, 1), (9, 11, 1), (11, 10, 1)])
        sage: list(feng_k_shortest_simple_paths(g, 1, 5))
        [[1, 6, 9, 5],
         [1, 7, 8, 5],
         [1, 2, 3, 4, 5],
         [1, 6, 9, 10, 5],
         [1, 6, 9, 11, 10, 5],
         [1, 6, 9, 3, 4, 5]]
        sage: list(feng_k_shortest_simple_paths(g, 1, 5, by_weight=True))
        [[1, 6, 9, 5],
         [1, 7, 8, 5],
         [1, 2, 3, 4, 5],
         [1, 6, 9, 10, 5],
         [1, 6, 9, 11, 10, 5],
         [1, 6, 9, 3, 4, 5]]
        sage: from sage.graphs.path_enumeration import feng_k_shortest_simple_paths
        sage: g = DiGraph([(1, 2, 5), (6, 3, 0), (2, 6, 6), (1, 4, 15),
        ....:              (4, 5, 1), (4, 3, 0), (7, 1, 2), (8, 7, 1)])
        sage: list(feng_k_shortest_simple_paths(g, 1, 3))
        [[1, 4, 3], [1, 2, 6, 3]]
        sage: list(feng_k_shortest_simple_paths(g, 1, 3, by_weight=True, report_edges=True, report_weight=True, labels=True))
        [(11, [(1, 2, 5), (2, 6, 6), (6, 3, 0)]), (15, [(1, 4, 15), (4, 3, 0)])]
        sage: list(feng_k_shortest_simple_paths(g, 1, 3, by_weight=True))
        [[1, 2, 6, 3], [1, 4, 3]]
        sage: from sage.graphs.path_enumeration import feng_k_shortest_simple_paths
        sage: G = DiGraph([(0, 1, 9), (0, 3, 1), (0, 4, 2), (1, 6, 4),
        ....:              (1, 7, 1), (2, 0, 5), (2, 1, 4), (2, 7, 1),
        ....:              (3, 1, 7), (3, 2, 4), (3, 4, 2), (4, 0, 8),
        ....:              (4, 1, 10), (4, 3, 3), (4, 7, 10), (5, 2, 5),
        ....:              (5, 4, 9), (6, 2, 9)], weighted=True)
        sage: list(feng_k_shortest_simple_paths(G, 2, 1, by_weight=True, report_weight=True, report_edges=True, labels=True))
        [(4, [(2, 1, 4)]),
         (13, [(2, 0, 5), (0, 3, 1), (3, 1, 7)]),
         (14, [(2, 0, 5), (0, 1, 9)]),
         (17, [(2, 0, 5), (0, 4, 2), (4, 1, 10)]),
         (17, [(2, 0, 5), (0, 4, 2), (4, 3, 3), (3, 1, 7)]),
         (18, [(2, 0, 5), (0, 3, 1), (3, 4, 2), (4, 1, 10)])]
    """
    if not self.is_directed():
        raise ValueError("this algorithm works only for directed graphs")

    if source not in self:
        raise ValueError("vertex '{}' is not in the graph".format(source))

    if target not in self:
        raise ValueError("vertex '{}' is not in the graph".format(target))

    if source == target:
        if report_edges:
            yield []
        elif report_weight:
            yield (0, [source])
        else:
            yield [source]
        return

    if self.has_loops() or self.allows_multiple_edges():
        G = self.to_simple(to_undirected=False, keep_label='min', immutable=False)
    else:
        G = self.copy()

    # removing the incoming edges to source and outgoing edges from target as
    # they do not contribute towards the k shortest simple paths
    G.delete_edges(G.incoming_edges(source, labels=False))
    G.delete_edges(G.outgoing_edges(target, labels=False))

    if weight_function is not None:
        by_weight = True

    if weight_function is None and by_weight:
        def weight_function(e):
            return e[2]

    by_weight, weight_function = self._get_weight_function(by_weight=by_weight,
                                                           weight_function=weight_function,
                                                           check_weight=check_weight)
    if by_weight:
        def reverse_weight_function(e):
            return weight_function((e[1], e[0], e[2]))
    else:
        def reverse_weight_function(e):
            return 1

    cdef dict edge_labels
    if report_edges and labels:
        edge_labels = {(e[0], e[1]): e for e in G.edge_iterator()}
        if not G.is_directed():
            for u, v in G.edge_iterator(labels=False):
                edge_labels[v, u] = edge_labels[u, v]

    from sage.graphs.base.boost_graph import shortest_paths
    # dictionary of parent node in the shortest path tree of the target vertex
    cdef dict parent = {}
    # assign color to each vertex as green, red or yellow
    cdef dict color = {}
    # express edges are the edges with head node as green and tail node as
    # yellow or tail node is a deviation node
    cdef dict expressEdges = {}
    # a dictionary of the new edges added to the graph used for restoring the
    # graph after the iteration
    cdef dict dic = {}
    # used to keep track of temporary edges added to the graph
    cdef dict temp_dict = {}
    # father of the path
    cdef dict father = {}

    def getUpStreamNodes(v):
        """
        If there exist a path in shortest path subtree of target node from u to
        v then u is said to be an upstream node of v
        """
        cdef list ver = list()
        S = [v]
        while S:
            u = S.pop(0)
            if u in parent:
                for u_node in parent[u]:
                    # if node color is green
                    if color[u_node] == 0:
                        S.append(u_node)
                        ver.append(u_node)
        return ver

    def findExpressEdges(Y):
        """
        Find the express edges whose tail nodes belong to a set Y and update the
        head node of each express edge
        """
        for v in Y:
            for w in G.neighbors_out(v):
                # if node color is green
                if color[w] == 0:
                    if w not in expressEdges:
                        expressEdges[w] = []
                    expressEdges[w].append((v, w, G.edge_label(v, w)))
                    if w != target and not G.has_edge(w, target):
                        G.add_edge(w, target, 0)
                        dic[w, target] = 1
                        reduced_cost[w, target] = 0
                    elif w != target and reduced_cost[w, target]:
                        temp_dict[w, target] = reduced_cost[w, target]
                        reduced_cost[w, target] = 0
                    include_vertices.add(w)

    reverse_graph = G.reverse()
    cdef dict dist
    cdef dict successor
    dist, successor = shortest_paths(reverse_graph, target, weight_function=reverse_weight_function,
                                     algorithm='Dijkstra_Boost')

    # successor is a child node in the shortest path subtree
    cdef dict reduced_cost = {(e[0], e[1]): weight_function(e) + dist[e[1]] - dist[e[0]]
                              for e in G.edge_iterator()
                              if e[0] in dist and e[1] in dist}

    cdef set exclude_vert_set = set(G) - set(dist)
    # finding the parent information from successor
    for key in successor:
        if successor[key] and successor[key] not in parent:
            parent[successor[key]] = [key]
        elif successor[key]:
            parent[successor[key]].append(key)

    def length_func(path):
        return sum(reduced_cost[e] for e in zip(path[:-1], path[1:]))

    # shortest path function for weighted/unweighted graph using reduced weights
    shortest_path_func = G._backend.bidirectional_dijkstra_special

    try:
        # compute the shortest path between the source and the target
        path = shortest_path_func(source, target, exclude_vertices=exclude_vert_set,
                                  weight_function=weight_function, reduced_weight=reduced_cost)
        # corner case
        if not path:
            if report_weight:
                yield (0, [])
            else:
                yield []
            return
    except Exception:
        if report_weight:
            yield (0, [])
        else:
            yield []
        return

    hash_path = tuple(path)
    father[hash_path] = None

    # heap data structure containing the candidate paths
    cdef priority_queue[pair[double, pair[int, int]]] heap_sorted_paths
    cdef int idx = 0
    cdef dict idx_to_path = {idx: path}
    heap_sorted_paths.push((-length_func(path), (idx, 0)))
    idx = idx + 1
    shortest_path_len = dist[source]

    cdef set exclude_edges
    cdef set exclude_vertices
    cdef set include_vertices
    cdef list allY
    cdef list prev_path, new_path, root
    cdef int path_idx, dev_idx

    while idx_to_path:
        # extracting the next best path from the heap
        cost, (path_idx, dev_idx) = heap_sorted_paths.top()
        heap_sorted_paths.pop()
        prev_path = idx_to_path[path_idx]
        # removing the path from dictionary
        del idx_to_path[path_idx]
        if len(set(prev_path)) == len(prev_path):
            if report_weight:
                cost = -cost
                if cost in ZZ:
                    cost = int(cost)
                if report_edges and labels:
                    yield (cost + shortest_path_len, [edge_labels[e] for e in zip(prev_path[:-1], prev_path[1:])])
                elif report_edges:
                    yield (cost + shortest_path_len, list(zip(prev_path[:-1], prev_path[1:])))
                else:
                    yield (cost + shortest_path_len, prev_path)
            else:
                if report_edges and labels:
                    yield [edge_labels[e] for e in zip(prev_path[:-1], prev_path[1:])]
                elif report_edges:
                    yield list(zip(prev_path[:-1], prev_path[1:]))
                else:
                    yield prev_path

        # deep copy of the exclude vertices set
        exclude_vertices = copy.deepcopy(exclude_vert_set)
        exclude_edges = set()
        include_vertices = set()
        expressEdges = {}
        dic = {}
        temp_dict = {}
        root = prev_path[:dev_idx]
        # coloring all the nodes as green initially
        color = {v: 0 for v in G}
        # list of yellow nodes
        allY = list()
        for j in range(dev_idx):
            exclude_vertices.add(prev_path[j])
        for j in range(dev_idx + 1):
            # coloring red
            color[prev_path[j]] = 1
            Yv = getUpStreamNodes(prev_path[j])
            allY += Yv
            for y in Yv:
                # color yellow for upstream nodes
                color[y] = 2
                include_vertices.add(y)
        # adding the deviation node to find the express edges
        allY.append(prev_path[dev_idx])
        findExpressEdges(allY)
        color[target] = 2
        include_vertices.add(target)
        # deviating from the previous path to find the candidate paths
        for i in range(dev_idx + 1, len(prev_path)):
            # root part of the previous path
            root.append(prev_path[i - 1])
            # if it is the deviation node
            if i == dev_idx + 1:
                p = father[tuple(prev_path)]
                # comparing the deviation nodes
                while p and len(p) > dev_idx + 1 and p[dev_idx] == root[dev_idx]:
                    # using fatherly approach to filter the edges to be removed
                    exclude_edges.add((p[i - 1], p[i]))
                    p = father[tuple(p)]
            else:
                # coloring it red
                color[root[-1]] = 1
                Yu = getUpStreamNodes(root[-1])
                for y in Yu:
                    # coloring upstream nodes as yellow
                    color[y] = 2
                    include_vertices.add(y)
                Yu.append(root[-1])
                for n in Yu:
                    if n in expressEdges and expressEdges[n]:
                        # recovering the express edges incident to a node n
                        for e in expressEdges[n]:
                            if (e[1], target) in dic:
                                # restoration of edges in the original graph
                                G.delete_edge(e[1], target)
                                del dic[(e[1], target)]
                                del reduced_cost[(e[1], target)]
                            if (e[1], target) in temp_dict:
                                # restoration of cost function
                                reduced_cost[e[1], target] = temp_dict[e[1], target]
                                del temp_dict[e[1], target]
                        # resetting the expressEdges for node n
                        expressEdges[n] = []
                findExpressEdges(Yu)
            # removing the edge in the previous shortest path to find a new
            # candidate path
            exclude_edges.add((prev_path[i - 1], prev_path[i]))
            try:
                # finding the spur part of the path after excluding certain
                # vertices and edges, this spur path is an all yellow subpath
                # so the shortest path algorithm is applied only on the all
                # yellow node subtree
                spur = shortest_path_func(root[-1], target,
                                          exclude_vertices=exclude_vertices,
                                          exclude_edges=exclude_edges,
                                          include_vertices=include_vertices,
                                          weight_function=weight_function,
                                          reduced_weight=reduced_cost)
                # finding the spur path in the original graph
                if spur and ((spur[-2], target) in dic or (spur[-2], target) in temp_dict):
                    spur.pop()
                    st = spur[-1]
                    while st != target:
                        st = successor[st]
                        spur.append(st)
                if not spur:
                    exclude_vertices.add(root[-1])
                    continue
                # concatenating the root and the spur path
                new_path = root[:-1] + spur
                # push operation
                hash_path = tuple(new_path)
                father[hash_path] = prev_path
                idx_to_path[idx] = new_path
                heap_sorted_paths.push((-length_func(new_path), (idx, i - 1)))
                idx = idx + 1
            except Exception:
                pass
            exclude_vertices.add(root[-1])
        # restoring the original graph here
        for e in dic:
            G.delete_edge(e)
            del reduced_cost[(e[0], e[1])]
        for e in temp_dict:
            reduced_cost[e[0], e[1]] = temp_dict[e[0], e[1]]


def _all_paths_iterator(self, vertex, ending_vertices=None,
                        simple=False, max_length=None, trivial=False,
                        use_multiedges=False, report_edges=False,
                        labels=False, data=None):
    r"""
    Return an iterator over the paths of ``self`` starting with the
    given vertex.

    INPUT:

    - ``vertex`` -- the starting vertex of the paths

    - ``ending_vertices`` -- iterable (default: ``None``); allowed ending
      vertices of the paths. If ``None``, then all vertices are allowed.

    - ``simple`` -- boolean (default: ``False``); if set to ``True``, then
      only simple paths are considered. Simple paths are paths in which no
      two arcs share a head or share a tail, i.e. every vertex in the path
      is entered at most once and exited at most once.

    - ``max_length`` -- nonnegative integer (default: ``None``); the
      maximum length of the enumerated paths. If set to ``None``, then all
      lengths are allowed.

    - ``trivial`` -- boolean (default: ``False``); if set to ``True``, then
      the empty paths are also enumerated

    - ``use_multiedges`` -- boolean (default: ``False``); this parameter is
      used only if the graph has multiple edges

        - If ``False``, the graph is considered as simple and an edge label
          is arbitrarily selected for each edge as in
          :meth:`sage.graphs.generic_graph.GenericGraph.to_simple` if
          ``report_edges`` is ``True``

        - If ``True``, a path will be reported as many times as the edges
          multiplicities along that path (when ``report_edges = False`` or
          ``labels = False``), or with all possible combinations of edge
          labels (when ``report_edges = True`` and ``labels = True``)

    - ``report_edges`` -- boolean (default: ``False``); whether to report
      paths as list of vertices (default) or list of edges, if ``False``
      then ``labels`` parameter is ignored

    - ``labels`` -- boolean (default: ``False``); if ``False``, each edge
      is simply a pair ``(u, v)`` of vertices. Otherwise a list of edges
      along with its edge labels are used to represent the path.

    - ``data`` -- dictionary (default: ``None``); optional parameter to
      pass information about edge multiplicities of the graph, if ``None``
      edge multiplicity values are computed inside the method.

    OUTPUT: iterator

    EXAMPLES::

        sage: G = graphs.CompleteGraph(4)
        sage: list(G._all_paths_iterator(1, ending_vertices=[3], simple=True))
        [[1, 3], [1, 0, 3], [1, 2, 3], [1, 0, 2, 3], [1, 2, 0, 3]]
        sage: list(G.shortest_simple_paths(1, 3))
        [[1, 3], [1, 0, 3], [1, 2, 3], [1, 2, 0, 3], [1, 0, 2, 3]]
        sage: pi = G._all_paths_iterator(1, ending_vertices=[3])
        sage: for _ in range(6):
        ....:     print(next(pi))
        [1, 3]
        [1, 0, 3]
        [1, 2, 3]
        [1, 0, 1, 3]
        [1, 0, 2, 3]
        [1, 2, 0, 3]

        sage: g = DiGraph({'a': ['a', 'b'], 'b': ['c'], 'c': ['d'], 'd': ['c']}, loops=True)
        sage: pi = g._all_paths_iterator('a', ending_vertices=['d'], report_edges=True, simple=True)
        sage: list(pi)
        [[('a', 'b'), ('b', 'c'), ('c', 'd')]]

        sage: g = DiGraph([(0, 1, 'a'), (0, 1, 'b'), (1, 2,'c'), (1, 2,'d')], multiedges=True)
        sage: pi =  g._all_paths_iterator(0, use_multiedges=True)
        sage: for _ in range(6):
        ....:     print(next(pi))
        [0, 1]
        [0, 1]
        [0, 1, 2]
        [0, 1, 2]
        [0, 1, 2]
        [0, 1, 2]
        sage: pi =  g._all_paths_iterator(0, use_multiedges=True, report_edges=True, labels=True)
        sage: for _ in range(6):
        ....:     print(next(pi))
        [(0, 1, 'b')]
        [(0, 1, 'a')]
        [(0, 1, 'b'), (1, 2, 'd')]
        [(0, 1, 'b'), (1, 2, 'c')]
        [(0, 1, 'a'), (1, 2, 'd')]
        [(0, 1, 'a'), (1, 2, 'c')]
        sage: list(g._all_paths_iterator(1, ending_vertices=[2], use_multiedges=False, report_edges=True, labels=True, simple=True))
        [[(1, 2, 'd')]]
        sage: list(g._all_paths_iterator(0, ending_vertices=[2], use_multiedges=False, report_edges=False, labels=True))
        [[0, 1, 2]]
        sage: list(g._all_paths_iterator(0, use_multiedges=True, report_edges=False, labels=True, max_length=1))
        [[0, 1], [0, 1]]
        sage: list(g._all_paths_iterator(0, use_multiedges=True, report_edges=True, labels=True, max_length=1))
        [[(0, 1, 'b')], [(0, 1, 'a')]]

        sage: g = DiGraph({'a': ['a', 'b'], 'b': ['c'], 'c': ['d'], 'd': ['c']}, loops=True)
        sage: pi = g._all_paths_iterator('a')
        sage: [len(next(pi)) - 1 for _ in range(5)]
        [1, 1, 2, 2, 2]

    ::

        sage: pi = g._all_paths_iterator('b')
        sage: for _ in range(5):
        ....:     print(next(pi))
        ['b', 'c']
        ['b', 'c', 'd']
        ['b', 'c', 'd', 'c']
        ['b', 'c', 'd', 'c', 'd']
        ['b', 'c', 'd', 'c', 'd', 'c']

    One may wish to enumerate simple paths, which are paths in which no two
    arcs share a head or share a tail, i.e. every vertex in the path is
    entered at most once and exited at most once. The result is always
    finite but may take a long time to compute::

        sage: pi = g._all_paths_iterator('a', simple=True)
        sage: sorted(pi)
        [['a', 'a'], ['a', 'b'], ['a', 'b', 'c'], ['a', 'b', 'c', 'd']]
        sage: pi = g._all_paths_iterator('d', simple=True)
        sage: sorted(pi)
        [['d', 'c'], ['d', 'c', 'd']]

    It is possible to specify the allowed ending vertices::

        sage: pi = g._all_paths_iterator('a', ending_vertices=['c'])
        sage: [len(next(pi)) - 1 for _ in range(5)]
        [2, 3, 4, 4, 5]
        sage: pi = g._all_paths_iterator('a', ending_vertices=['a', 'b'])
        sage: [len(next(pi)) - 1 for _ in range(5)]
        [1, 1, 2, 2, 3]

    One can bound the length of the paths::

        sage: pi = g._all_paths_iterator('d', max_length=3)
        sage: sorted(pi)
        [['d', 'c'], ['d', 'c', 'd'], ['d', 'c', 'd', 'c']]

    Or include the trivial empty path::

        sage: pi = g._all_paths_iterator('a', max_length=3, trivial=True)
        sage: sorted(list(pi), key=lambda x:(len(x), x))
        [['a'], ['a', 'a'], ['a', 'b'], ['a', 'a', 'a'], ['a', 'a', 'b'],
            ['a', 'b', 'c'], ['a', 'a', 'a', 'a'], ['a', 'a', 'a', 'b'],
            ['a', 'a', 'b', 'c'], ['a', 'b', 'c', 'd']]
        sage: pi = g._all_paths_iterator('a', max_length=3, trivial=True)
        sage: [len(p) - 1 for p in pi]
        [0, 1, 1, 2, 2, 2, 3, 3, 3, 3]
    """
    if ending_vertices is None:
        ending_vertices = self
    else:
        ending_vertices = frozenset(ending_vertices)
    if max_length is None:
        from sage.rings.infinity import Infinity
        max_length = Infinity
    if max_length < 1:
        return

    if self.is_directed():
        neighbor_iterator = self.neighbor_out_iterator
    else:
        neighbor_iterator = self.neighbor_iterator

    cdef dict my_dict = {}
    cdef dict edge_multiplicity
    if not data:
        if report_edges and labels:
            if use_multiedges:
                for e in self.edge_iterator():
                    if (e[0], e[1]) in my_dict:
                        my_dict[(e[0], e[1])].append(e)
                    else:
                        my_dict[(e[0], e[1])] = [e]
            else:
                for e in self.edge_iterator():
                    if (e[0], e[1]) not in my_dict:
                        my_dict[(e[0], e[1])] = [e]
        elif use_multiedges and self.has_multiple_edges():
            from collections import Counter
            edge_multiplicity = dict(Counter(self.edge_iterator(labels=False)))
    else:
        if report_edges and labels:
            my_dict = data
        elif use_multiedges and self.has_multiple_edges():
            edge_multiplicity = data
    # Start with the empty path; we will try all extensions of it
    cdef list queue = []
    cdef list path = [vertex]
    cdef list newpath
    cdef int m
    if trivial and not report_edges and vertex in ending_vertices:
        yield path
    while True:
        # Build next generation of paths, one arc longer; max_length refers
        # to edges and not vertices, hence <= and not <
        if len(path) <= max_length:
            # We try all possible extensions
            if simple:
                # We only keep simple extensions. An extension is simple iff
                # the new vertex being entered has not previously occurred
                # in the path, or has occurred but only been exited (i.e. is
                # the first vertex in the path). In this latter case we must
                # not exit the new vertex again, so we do not consider it
                # for further extension, but just yield it immediately. See
                # Issue #12385.
                frozen_path = frozenset(path)
                for neighbor in neighbor_iterator(path[-1]):
                    if neighbor not in frozen_path:
                        queue.append(path + [neighbor])
                    elif (neighbor == path[0] and
                          neighbor in ending_vertices):
                        newpath = path + [neighbor]
                        if report_edges and labels:
                            for p in product(*[my_dict[e] for e in zip(newpath[:-1], newpath[1:])]):
                                yield list(p)
                        elif use_multiedges and self.has_multiple_edges():
                            m = prod(edge_multiplicity[e] for e in zip(newpath[:-1], newpath[1:]))
                            if report_edges:
                                newpath = list(zip(newpath[:-1], newpath[1:]))
                            for _ in range(m):
                                yield newpath
                        elif report_edges:
                            yield list(zip(newpath[:-1], newpath[1:]))
                        else:
                            yield newpath
            else:
                # Non-simple paths requested: we add all of them
                for neighbor in neighbor_iterator(path[-1]):
                    queue.append(path + [neighbor])

        if not queue:
            break
        path = queue.pop(0)     # get the next path

        if path[-1] in ending_vertices:
            # yield good path
            if report_edges and labels:
                for p in product(*[my_dict[e] for e in zip(path[:-1], path[1:])]):
                    yield list(p)
            elif use_multiedges and self.has_multiple_edges():
                m = prod(edge_multiplicity[e] for e in zip(path[:-1], path[1:]))
                if report_edges:
                    newpath = list(zip(path[:-1], path[1:]))
                else:
                    newpath = path
                for _ in range(m):
                    yield newpath
            elif report_edges:
                yield list(zip(path[:-1], path[1:]))
            else:
                yield path


def all_paths_iterator(self, starting_vertices=None, ending_vertices=None,
                       simple=False, max_length=None, trivial=False,
                       use_multiedges=False, report_edges=False, labels=False):
    r"""
    Return an iterator over the paths of ``self``.

    The paths are enumerated in increasing length order.

    INPUT:

    - ``starting_vertices`` -- iterable (default: ``None``); vertices from
      which the paths must start. If ``None``, then all vertices of the
      graph can be starting points.

    - ``ending_vertices`` -- iterable (default: ``None``); allowed ending
      vertices of the paths. If ``None``, then all vertices are allowed.

    - ``simple`` -- boolean (default: ``False``); if set to ``True``, then
      only simple paths are considered. Simple paths are paths in which no
      two arcs share a head or share a tail, i.e. every vertex in the path
      is entered at most once and exited at most once.

    - ``max_length`` -- nonnegative integer (default: ``None``); the
      maximum length of the enumerated paths. If set to ``None``, then all
      lengths are allowed.

    - ``trivial`` -- boolean (default: ``False``); if set to ``True``, then
      the empty paths are also enumerated

    - ``use_multiedges`` -- boolean (default: ``False``); this parameter is
      used only if the graph has multiple edges

        - If ``False``, the graph is considered as simple and an edge label
          is arbitrarily selected for each edge as in
          :meth:`sage.graphs.generic_graph.GenericGraph.to_simple` if
          ``report_edges`` is ``True``

        - If ``True``, a path will be reported as many times as the edges
          multiplicities along that path (when ``report_edges = False`` or
          ``labels = False``), or with all possible combinations of edge
          labels (when ``report_edges = True`` and ``labels = True``)

    - ``report_edges`` -- boolean (default: ``False``); whether to report
      paths as list of vertices (default) or list of edges, if ``False``
      then ``labels`` parameter is ignored

    - ``labels`` -- boolean (default: ``False``); if ``False``, each edge
      is simply a pair ``(u, v)`` of vertices. Otherwise a list of edges
      along with its edge labels are used to represent the path.

    OUTPUT: iterator

    AUTHOR:

        Alexandre Blondin Masse

    EXAMPLES::

        sage: G = graphs.CompleteGraph(4)
        sage: list(G.all_paths_iterator(starting_vertices=[1], ending_vertices=[3], simple=True))
        [[1, 3], [1, 0, 3], [1, 2, 3], [1, 0, 2, 3], [1, 2, 0, 3]]
        sage: list(G.shortest_simple_paths(1, 3))
        [[1, 3], [1, 0, 3], [1, 2, 3], [1, 2, 0, 3], [1, 0, 2, 3]]
        sage: pi = G.all_paths_iterator(starting_vertices=[1], ending_vertices=[3])
        sage: for _ in range(6):
        ....:     print(next(pi))
        [1, 3]
        [1, 0, 3]
        [1, 2, 3]
        [1, 0, 1, 3]
        [1, 0, 2, 3]
        [1, 2, 0, 3]

        sage: g = DiGraph({'a': ['a', 'b'], 'b': ['c'], 'c': ['d'], 'd': ['c']}, loops=True)
        sage: pi = g.all_paths_iterator(starting_vertices=['a'], ending_vertices=['d'], report_edges=True, simple=True)
        sage: list(pi)
        [[('a', 'b'), ('b', 'c'), ('c', 'd')]]

        sage: g = DiGraph([(0, 1, 'a'), (0, 1, 'b'), (1, 2,'c'), (1, 2,'d')], multiedges=True)
        sage: pi =  g.all_paths_iterator(starting_vertices=[0], use_multiedges=True)
        sage: for _ in range(6):
        ....:     print(next(pi))
        [0, 1]
        [0, 1]
        [0, 1, 2]
        [0, 1, 2]
        [0, 1, 2]
        [0, 1, 2]
        sage: pi =  g.all_paths_iterator(starting_vertices=[0], use_multiedges=True, report_edges=True, labels=True)
        sage: for _ in range(6):
        ....:     print(next(pi))
        [(0, 1, 'b')]
        [(0, 1, 'a')]
        [(0, 1, 'b'), (1, 2, 'd')]
        [(0, 1, 'b'), (1, 2, 'c')]
        [(0, 1, 'a'), (1, 2, 'd')]
        [(0, 1, 'a'), (1, 2, 'c')]
        sage: list(g.all_paths_iterator(starting_vertices=[0, 1], ending_vertices=[2], use_multiedges=False, report_edges=True, labels=True, simple=True))
        [[(1, 2, 'd')], [(0, 1, 'b'), (1, 2, 'd')]]
        sage: list(g.all_paths_iterator(starting_vertices=[0, 1], ending_vertices=[2], use_multiedges=False, report_edges=False, labels=True))
        [[1, 2], [0, 1, 2]]
        sage: list(g.all_paths_iterator(use_multiedges=True, report_edges=False, labels=True, max_length=1))
        [[1, 2], [1, 2], [0, 1], [0, 1]]
        sage: list(g.all_paths_iterator(use_multiedges=True, report_edges=True, labels=True, max_length=1))
        [[(1, 2, 'd')], [(1, 2, 'c')], [(0, 1, 'b')], [(0, 1, 'a')]]

        sage: g = DiGraph({'a': ['a', 'b'], 'b': ['c'], 'c': ['d'], 'd': ['c']}, loops=True)
        sage: pi = g.all_paths_iterator()
        sage: [len(next(pi)) - 1 for _ in range(7)]
        [1, 1, 1, 1, 1, 2, 2]

    It is possible to precise the allowed starting and/or ending vertices::

        sage: pi = g.all_paths_iterator(starting_vertices=['a'])
        sage: [len(next(pi)) - 1 for _ in range(5)]
        [1, 1, 2, 2, 2]
        sage: pi = g.all_paths_iterator(starting_vertices=['a'], ending_vertices=['b'])
        sage: for _ in range(5):
        ....:     print(next(pi))
        ['a', 'b']
        ['a', 'a', 'b']
        ['a', 'a', 'a', 'b']
        ['a', 'a', 'a', 'a', 'b']
        ['a', 'a', 'a', 'a', 'a', 'b']

    One may prefer to enumerate only simple paths (see
    :meth:`all_simple_paths`)::

        sage: pi = g.all_paths_iterator(simple=True)
        sage: sorted(list(pi), key=lambda x:(len(x), x))
        [['a', 'a'], ['a', 'b'], ['b', 'c'], ['c', 'd'], ['d', 'c'],
            ['a', 'b', 'c'], ['b', 'c', 'd'], ['c', 'd', 'c'], ['d', 'c', 'd'],
            ['a', 'b', 'c', 'd']]
        sage: pi = g.all_paths_iterator(simple=True)
        sage: [len(p) - 1 for p in pi]
        [1, 1, 1, 1, 1, 2, 2, 2, 2, 3]

    Or simply bound the length of the enumerated paths::

        sage: pi = g.all_paths_iterator(starting_vertices=['a'], ending_vertices=['b', 'c'], max_length=6)
        sage: sorted(list(pi), key=lambda x:(len(x), x))
        [['a', 'b'], ['a', 'a', 'b'], ['a', 'b', 'c'], ['a', 'a', 'a', 'b'],
            ['a', 'a', 'b', 'c'], ['a', 'a', 'a', 'a', 'b'],
            ['a', 'a', 'a', 'b', 'c'], ['a', 'b', 'c', 'd', 'c'],
            ['a', 'a', 'a', 'a', 'a', 'b'], ['a', 'a', 'a', 'a', 'b', 'c'],
            ['a', 'a', 'b', 'c', 'd', 'c'],
            ['a', 'a', 'a', 'a', 'a', 'a', 'b'],
            ['a', 'a', 'a', 'a', 'a', 'b', 'c'],
            ['a', 'a', 'a', 'b', 'c', 'd', 'c'],
            ['a', 'b', 'c', 'd', 'c', 'd', 'c']]
        sage: pi = g.all_paths_iterator(starting_vertices=['a'], ending_vertices=['b', 'c'], max_length=6)
        sage: [len(p) - 1 for p in pi]
        [1, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6]

    By default, empty paths are not enumerated, but it may be parametrized::

        sage: pi = g.all_paths_iterator(simple=True, trivial=True)
        sage: sorted(list(pi), key=lambda x:(len(x), x))
        [['a'], ['b'], ['c'], ['d'], ['a', 'a'], ['a', 'b'], ['b', 'c'],
            ['c', 'd'], ['d', 'c'], ['a', 'b', 'c'], ['b', 'c', 'd'],
            ['c', 'd', 'c'], ['d', 'c', 'd'], ['a', 'b', 'c', 'd']]
        sage: pi = g.all_paths_iterator(simple=True, trivial=True)
        sage: [len(p) - 1 for p in pi]
        [0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3]
        sage: pi = g.all_paths_iterator(simple=True, trivial=False)
        sage: sorted(list(pi), key=lambda x:(len(x), x))
        [['a', 'a'], ['a', 'b'], ['b', 'c'], ['c', 'd'], ['d', 'c'],
            ['a', 'b', 'c'], ['b', 'c', 'd'], ['c', 'd', 'c'], ['d', 'c', 'd'],
            ['a', 'b', 'c', 'd']]
        sage: pi = g.all_paths_iterator(simple=True, trivial=False)
        sage: [len(p) - 1 for p in pi]
        [1, 1, 1, 1, 1, 2, 2, 2, 2, 3]
    """
    if starting_vertices is None:
        starting_vertices = self
    cdef dict data = {}
    cdef dict vertex_iterators
    cdef list path

    if report_edges and labels:
        if use_multiedges:
            for e in self.edge_iterator():
                if (e[0], e[1]) in data:
                    data[(e[0], e[1])].append(e)
                else:
                    data[(e[0], e[1])] = [e]
        else:
            for e in self.edge_iterator():
                if (e[0], e[1]) not in data:
                    data[(e[0], e[1])] = [e]
    elif use_multiedges and self.has_multiple_edges():
        from collections import Counter
        edge_multiplicity = Counter(self.edge_iterator(labels=False))
        data = dict(edge_multiplicity)

    # We create one paths iterator per vertex
    # This is necessary if we want to iterate over paths
    # with increasing length
    vertex_iterators = {v: self._all_paths_iterator(v, ending_vertices=ending_vertices,
                                                    simple=simple, max_length=max_length,
                                                    trivial=trivial, use_multiedges=use_multiedges,
                                                    report_edges=report_edges, labels=labels, data=data)
                        for v in starting_vertices}

    cdef priority_queue[pair[int, int]] pq
    cdef int idx = 0
    cdef dict idx_to_path = {}
    for vi in vertex_iterators.values():
        try:
            path = next(vi)
            idx_to_path[idx] = path
            pq.push((-len(path), idx))
            idx = idx + 1
        except(StopIteration):
            pass
    # Since we always extract a shortest path, using a heap
    # can speed up the algorithm
    while not pq.empty():
        # We choose the shortest available path
        _, shortest_path_idx = pq.top()
        pq.pop()
        prev_path = idx_to_path[shortest_path_idx]
        yield prev_path
        del idx_to_path[shortest_path_idx]
        # We update the path iterator to its next available path if it exists
        try:
            if report_edges:
                path = next(vertex_iterators[prev_path[0][0]])
            else:
                path = next(vertex_iterators[prev_path[0]])
            idx_to_path[idx] = path
            pq.push((-len(path), idx))
            idx = idx + 1
        except(StopIteration):
            pass


def all_simple_paths(self, starting_vertices=None, ending_vertices=None,
                     max_length=None, trivial=False, use_multiedges=False,
                     report_edges=False, labels=False):
    r"""
    Return a list of all the simple paths of ``self`` starting with one of
    the given vertices.

    Simple paths are paths in which no two arcs share a head or share a
    tail, i.e. every vertex in the path is entered at most once and exited
    at most once.

    INPUT:

    - ``starting_vertices`` -- list (default: ``None``); vertices from which
      the paths must start. If ``None``, then all vertices of the graph can
      be starting points.

    - ``ending_vertices`` -- iterable (default: ``None``); allowed ending
      vertices of the paths. If ``None``, then all vertices are allowed.

    - ``max_length`` -- nonnegative integer (default: ``None``); the
      maximum length of the enumerated paths. If set to ``None``, then all
      lengths are allowed.

    - ``trivial`` -- boolean (default: ``False``); if set to ``True``, then
      the empty paths are also enumerated

    - ``use_multiedges`` -- boolean (default: ``False``); this parameter is
      used only if the graph has multiple edges

        - If ``False``, the graph is considered as simple and an edge label
          is arbitrarily selected for each edge as in
          :meth:`sage.graphs.generic_graph.GenericGraph.to_simple` if
          ``report_edges`` is ``True``

        - If ``True``, a path will be reported as many times as the edges
          multiplicities along that path (when ``report_edges = False`` or
          ``labels = False``), or with all possible combinations of edge
          labels (when ``report_edges = True`` and ``labels = True``)

    - ``report_edges`` -- boolean (default: ``False``); whether to report
      paths as list of vertices (default) or list of edges, if ``False``
      then ``labels`` parameter is ignored

    - ``labels`` -- boolean (default: ``False``); if ``False``, each edge
      is simply a pair ``(u, v)`` of vertices. Otherwise a list of edges
      along with its edge labels are used to represent the path.

    OUTPUT: list

    .. NOTE::

        Although the number of simple paths of a finite graph is always
        finite, computing all its paths may take a very long time.

    EXAMPLES::

        sage: G = graphs.CompleteGraph(4)
        sage: G.all_simple_paths([1], [3])
        [[1, 3], [1, 0, 3], [1, 2, 3], [1, 0, 2, 3], [1, 2, 0, 3]]
        sage: list(G.shortest_simple_paths(1, 3))
        [[1, 3], [1, 0, 3], [1, 2, 3], [1, 2, 0, 3], [1, 0, 2, 3]]
        sage: G.all_simple_paths([0, 1], [2, 3])
        [[1, 2],
         [1, 3],
         [0, 2],
         [0, 3],
         [0, 1, 2],
         [0, 1, 3],
         [0, 2, 3],
         [0, 3, 2],
         [1, 0, 2],
         [1, 0, 3],
         [1, 2, 3],
         [1, 3, 2],
         [1, 0, 2, 3],
         [1, 0, 3, 2],
         [1, 2, 0, 3],
         [1, 3, 0, 2],
         [0, 1, 2, 3],
         [0, 1, 3, 2],
         [0, 2, 1, 3],
         [0, 3, 1, 2]]

        sage: g = DiGraph({0: [0, 1], 1: [2], 2: [3], 3: [2]}, loops=True)
        sage: g.all_simple_paths()
        [[3, 2],
         [2, 3],
         [1, 2],
         [0, 0],
         [0, 1],
         [0, 1, 2],
         [1, 2, 3],
         [2, 3, 2],
         [3, 2, 3],
         [0, 1, 2, 3]]

        sage: g = DiGraph([(0, 1, 'a'), (0, 1, 'b'), (1, 2,'c'), (1, 2,'d')], multiedges=True)
        sage: g.all_simple_paths(starting_vertices=[0], ending_vertices=[2], use_multiedges=False)
        [[0, 1, 2]]
        sage: g.all_simple_paths(starting_vertices=[0], ending_vertices=[2], use_multiedges=True)
        [[0, 1, 2], [0, 1, 2], [0, 1, 2], [0, 1, 2]]
        sage: g.all_simple_paths(starting_vertices=[0], ending_vertices=[2], use_multiedges=True, report_edges=True)
        [[(0, 1), (1, 2)], [(0, 1), (1, 2)], [(0, 1), (1, 2)], [(0, 1), (1, 2)]]
        sage: g.all_simple_paths(starting_vertices=[0], ending_vertices=[2], use_multiedges=True, report_edges=True, labels=True)
        [[(0, 1, 'b'), (1, 2, 'd')],
         [(0, 1, 'b'), (1, 2, 'c')],
         [(0, 1, 'a'), (1, 2, 'd')],
         [(0, 1, 'a'), (1, 2, 'c')]]
        sage: g.all_simple_paths(starting_vertices=[0, 1], ending_vertices=[2], use_multiedges=False, report_edges=True, labels=True)
        [[(1, 2, 'd')], [(0, 1, 'b'), (1, 2, 'd')]]
        sage: g.all_simple_paths(starting_vertices=[0, 1], ending_vertices=[2], use_multiedges=False, report_edges=False, labels=True)
        [[1, 2], [0, 1, 2]]
        sage: g.all_simple_paths(use_multiedges=True, report_edges=False, labels=True)
        [[1, 2], [1, 2], [0, 1], [0, 1], [0, 1, 2], [0, 1, 2], [0, 1, 2], [0, 1, 2]]
        sage: g.all_simple_paths(starting_vertices=[0, 1], ending_vertices=[2], use_multiedges=False, report_edges=True, labels=True, trivial=True)
        [[(1, 2, 'd')], [(0, 1, 'b'), (1, 2, 'd')]]

    One may compute all paths having specific starting and/or ending
    vertices::

        sage: g = DiGraph({'a': ['a', 'b'], 'b': ['c'], 'c': ['d'], 'd': ['c']}, loops=True)
        sage: g.all_simple_paths(starting_vertices=['a'])
        [['a', 'a'], ['a', 'b'], ['a', 'b', 'c'], ['a', 'b', 'c', 'd']]
        sage: g.all_simple_paths(starting_vertices=['a'], ending_vertices=['c'])
        [['a', 'b', 'c']]
        sage: g.all_simple_paths(starting_vertices=['a'], ending_vertices=['b', 'c'])
        [['a', 'b'], ['a', 'b', 'c']]

    It is also possible to bound the length of the paths::

        sage: g = DiGraph({0: [0, 1], 1: [2], 2: [3], 3: [2]}, loops=True)
        sage: g.all_simple_paths(max_length=2)
        [[3, 2],
         [2, 3],
         [1, 2],
         [0, 0],
         [0, 1],
         [0, 1, 2],
         [1, 2, 3],
         [2, 3, 2],
         [3, 2, 3]]

    By default, empty paths are not enumerated, but this can be
    parametrized::

        sage: g = DiGraph({'a': ['a', 'b'], 'b': ['c'], 'c': ['d'], 'd': ['c']}, loops=True)
        sage: g.all_simple_paths(starting_vertices=['a'], trivial=True)
        [['a'], ['a', 'a'], ['a', 'b'], ['a', 'b', 'c'],
         ['a', 'b', 'c', 'd']]
        sage: g.all_simple_paths(starting_vertices=['a'], trivial=False)
        [['a', 'a'], ['a', 'b'], ['a', 'b', 'c'], ['a', 'b', 'c', 'd']]
    """
    return list(self.all_paths_iterator(starting_vertices=starting_vertices,
                                        ending_vertices=ending_vertices,
                                        simple=True, max_length=max_length,
                                        trivial=trivial, use_multiedges=use_multiedges,
                                        report_edges=report_edges, labels=labels))
