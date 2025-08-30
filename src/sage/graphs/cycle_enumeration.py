r"""
Cycle enumeration

This module is meant for all functions related to cycle enumeration in graphs.

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :func:`all_cycles_iterator` | Return an iterator over the cycles of ``self``.
    :func:`all_simple_cycles` | Return an iterator over the simple cycles of ``self``.

Functions
---------
"""
# ****************************************************************************
# Copyright (C) 2025 Yuta Inoue <yutainoue888@gmail.com>
#                    David Coudert <david.coudert@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from copy import copy


def _all_cycles_iterator_vertex(self, vertex, starting_vertices=None, simple=False,
                                rooted=False, max_length=None, trivial=False,
                                remove_acyclic_edges=True,
                                weight_function=None, by_weight=False,
                                check_weight=True, report_weight=False):
    r"""
    Return an iterator over the cycles of ``self`` starting with the given
    vertex in increasing length order. Each edge must have a positive weight.

    INPUT:

    - ``vertex`` -- the starting vertex of the cycle

    - ``starting_vertices`` -- iterable (default: ``None``); vertices from
      which the cycles must start. If ``None``, then all vertices of the
      graph can be starting points. This argument is necessary if ``rooted``
      is set to ``True``.

    - ``simple`` -- boolean (default: ``False``); if set to ``True``, then
      only simple cycles are considered. A cycle is simple if the only
      vertex occurring twice in it is the starting and ending one and no edge
      occurs twice.

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
      acyclic edges must be removed from the graph if ``self`` is directed.
      Used to avoid recomputing it for each vertex

    - ``weight_function`` -- function (default: ``None``); a function that
      takes as input an edge ``(u, v, l)`` and outputs its weight. If not
      ``None``, ``by_weight`` is automatically set to ``True``. If ``None``
      and ``by_weight`` is ``True``, we use the edge label ``l`` as a
      weight.

    - ``by_weight`` -- boolean (default: ``False``); if ``True``, the edges
      in the graph are weighted, otherwise all edges have weight 1

    - ``check_weight`` -- boolean (default: ``True``); whether to check that
      the ``weight_function`` outputs a number for each edge

    - ``report_weight`` -- boolean (default: ``False``); if ``False``, just
      a cycle is returned. Otherwise a tuple of cycle length and cycle is
      returned.

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

    A cycle is enumerated in increasing length order for a weighted graph::

        sage: g = DiGraph()
        sage: g.add_edges([('a', 'b', 2), ('a', 'e', 2), ('b', 'a', 1), ('b', 'c', 1),
        ....:              ('c', 'd', 2), ('d', 'b', 1), ('e', 'a', 2)])
        sage: it = g._all_cycles_iterator_vertex('a', simple=False, max_length=None,
        ....:                                    by_weight=True, report_weight=True)
        sage: for i in range(7): print(next(it))
        (3, ['a', 'b', 'a'])
        (4, ['a', 'e', 'a'])
        (6, ['a', 'b', 'a', 'b', 'a'])
        (7, ['a', 'b', 'a', 'e', 'a'])
        (7, ['a', 'b', 'c', 'd', 'b', 'a'])
        (7, ['a', 'e', 'a', 'b', 'a'])
        (8, ['a', 'e', 'a', 'e', 'a'])

    Each edge must have a positive weight::

        sage: g = DiGraph()
        sage: g.add_edges([('a', 'b', -2), ('b', 'a', 1)])
        sage: next(g._all_cycles_iterator_vertex('a', simple=False, max_length=None,
        ....:                                     by_weight=True, report_weight=True))
        Traceback (most recent call last):
        ...
        ValueError: negative weight is not allowed

    The function works for an undirected graph. Specifically, each cycle is
    enumerated exactly once, meaning a cycle and its reverse are not listed separately::

        sage: g = Graph({0: [1, 2], 1: [0, 2], 2: [0, 1]})
        sage: it = g._all_cycles_iterator_vertex(0, simple=False)
        sage: for i in range(7): print(next(it))
        [0, 1, 0]
        [0, 2, 0]
        [0, 1, 2, 0]
        [0, 1, 0, 1, 0]
        [0, 1, 0, 2, 0]
        [0, 1, 2, 1, 0]
        [0, 2, 0, 2, 0]
        sage: for cycle in g._all_cycles_iterator_vertex(0, simple=True):
        ....:     print(cycle)
        [0, 1, 2, 0]
    """
    if starting_vertices is None:
        starting_vertices = [vertex]
    # First enumerate the empty cycle
    if trivial:
        if report_weight:
            yield (0, [vertex])
        else:
            yield [vertex]
    # First we remove vertices and edges that are not part of any cycle
    if self.is_directed() and remove_acyclic_edges:
        sccs = self.strongly_connected_components()
        if len(sccs) == 1:
            h = self
        else:
            d = {}
            for id, component in enumerate(sccs):
                for v in component:
                    d[v] = id
            h = copy(self)
            h.delete_edges((u, v) for u, v in h.edge_iterator(labels=False) if d[u] != d[v])
    else:
        h = self
        int_to_vertex = list(h)
        vertex_to_int = {v: i for i, v in enumerate(int_to_vertex)}

    by_weight, weight_function = self._get_weight_function(by_weight=by_weight,
                                                           weight_function=weight_function,
                                                           check_weight=check_weight)
    if by_weight:
        for e in h.edge_iterator():
            if weight_function(e) < 0:
                raise ValueError("negative weight is not allowed")

    from heapq import heapify, heappop, heappush
    heap_queue = [(0, [vertex])]
    if max_length is None:
        from sage.rings.infinity import Infinity
        max_length = Infinity
    while heap_queue:
        length, path = heappop(heap_queue)
        # Checks if a cycle has been found
        if len(path) > 1 and path[0] == path[-1]:
            report = True
            if not self.is_directed():
                if simple:
                    report = len(path) > 3 and vertex_to_int[path[1]] < vertex_to_int[path[-2]]
                else:
                    L = len(path)
                    for i in range(1, L // 2):
                        if vertex_to_int[path[i]] > vertex_to_int[path[L - i - 1]]:
                            report = False
                            break
                        if vertex_to_int[path[i]] < vertex_to_int[path[L - i - 1]]:
                            break
            if report:
                if report_weight:
                    yield (length, path)
                else:
                    yield path
        # If simple is set to True, only simple cycles are
        # allowed, Then it discards the current path
        if (not simple or path.count(path[-1]) == 1):
            for e in h.edge_iterator(vertices=[path[-1]]):
                neighbor = e[1] if e[0] == path[-1] else e[0]
                # Makes sure that the current cycle is not too long.
                # If cycles are not rooted, makes sure to keep only the
                # minimum cycle according to the lexicographic order
                if length + weight_function(e) <= max_length and \
                   (rooted or neighbor not in starting_vertices or path[0] <= neighbor):
                    heappush(heap_queue, (length + weight_function(e), path + [neighbor]))


def _all_simple_cycles_iterator_edge(self, edge, max_length=None,
                                     remove_unnecessary_edges=True,
                                     weight_function=None, by_weight=False,
                                     check_weight=True, report_weight=False):
    r"""
    Return an iterator over the **simple** cycles of ``self`` starting with the
    given edge in increasing length order. Each edge must have a positive weight.

    INPUT:

    - ``edge`` -- the starting edge of the cycle.

    - ``max_length`` -- nonnegative integer (default: ``None``); the
      maximum length of the enumerated paths. If set to ``None``, then all
      lengths are allowed.

    - ``remove_unnecessary_edges`` -- boolean (default: ``True``); whether
      unnecessary edges for enumerating simple cycles must be removed from the
      graph. If ``self`` is directed, edges not in the strongly connected
      component that contains ``edge`` are unnecessary. Otherwise, edges not in
      the 2-connected component that contains ``edge`` are unnecessary.

    - ``weight_function`` -- function (default: ``None``); a function that
      takes as input an edge ``(u, v, l)`` and outputs its weight. If not
      ``None``, ``by_weight`` is automatically set to ``True``. If ``None``
      and ``by_weight`` is ``True``, we use the edge label ``l`` as a
      weight.

    - ``by_weight`` -- boolean (default: ``False``); if ``True``, the edges
      in the graph are weighted, otherwise all edges have weight 1

    - ``check_weight`` -- boolean (default: ``True``); whether to check that
      the ``weight_function`` outputs a number for each edge

    - ``report_weight`` -- boolean (default: ``False``); if ``False``, just
      a cycle is returned. Otherwise a tuple of cycle length and cycle is
      returned.

    OUTPUT: iterator

    ALGORITHM:

    Given an edge `uv`, this algorithm extracts k-shortest `vu`-path in `G-uv`
    in increasing length order by using k-shortest path algorithm. Thus, it
    extracts only simple cycles. See
    :math:`~sage.graphs.path_enumeration.shortest_simple_paths` for more
    information.

    EXAMPLES:

        sage: g = graphs.Grid2dGraph(2, 5).to_directed()
        sage: it = g._all_simple_cycles_iterator_edge(((0, 0), (0, 1), None), report_weight=True)
        sage: for i in range(5): print(next(it))
        (2.0, [(0, 0), (0, 1), (0, 0)])
        (4.0, [(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)])
        (6.0, [(0, 0), (0, 1), (0, 2), (1, 2), (1, 1), (1, 0), (0, 0)])
        (8.0, [(0, 0), (0, 1), (0, 2), (0, 3), (1, 3), (1, 2), (1, 1), (1, 0), (0, 0)])
        (10.0, [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (1, 4), (1, 3), (1, 2), (1, 1), (1, 0), (0, 0)])

    The function works for undirected graphs as well::

        sage: g = graphs.Grid2dGraph(2, 5)
        sage: it = g._all_simple_cycles_iterator_edge(((0, 0), (0, 1), None), report_weight=True)
        sage: for i in range(4): print(next(it))
        (4, [(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)])
        (6, [(0, 0), (0, 1), (0, 2), (1, 2), (1, 1), (1, 0), (0, 0)])
        (8, [(0, 0), (0, 1), (0, 2), (0, 3), (1, 3), (1, 2), (1, 1), (1, 0), (0, 0)])
        (10, [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (1, 4), (1, 3), (1, 2), (1, 1), (1, 0), (0, 0)])

    Each edge must have a positive weight::

        sage: g = DiGraph()
        sage: g.add_edges([('a', 'b', -2), ('b', 'a', 1)])
        sage: next(g._all_simple_cycles_iterator_edge(('a', 'b', -2), max_length=None, by_weight=True))
        Traceback (most recent call last):
        ...
        ValueError: negative weight is not allowed

    TESTS:

    A graph with a loop::

        sage: g = DiGraph(loops=True)
        sage: g.add_edges([('a', 'b'), ('b', 'c'), ('c', 'a'), ('c', 'c')])
        sage: it = g._all_simple_cycles_iterator_edge(('a', 'b'), remove_unnecessary_edges=True)
        sage: next(it)
        ['a', 'b', 'c', 'a']
        sage: g = DiGraph(loops=True)
        sage: g.add_edges([('a', 'b'), ('b', 'c'), ('c', 'c')])
        sage: it = g._all_simple_cycles_iterator_edge(('c', 'c'), remove_unnecessary_edges=True)
        sage: next(it)
        ['c', 'c']

    An undirected graph with 2 biconnected components::

        sage: g = Graph()
        sage: g.add_edges([('a', 'b'), ('b', 'c'), ('c', 'a'), ('a', 'd'), ('d', 'e'), ('e', 'a')])
        sage: it = g._all_simple_cycles_iterator_edge(('a', 'b'), remove_unnecessary_edges=True)
        sage: next(it)
        ['a', 'b', 'c', 'a']
    """
    # First we remove vertices and edges that are not part of any cycle
    if remove_unnecessary_edges:
        if self.is_directed():
            components = self.strongly_connected_components_subgraphs()
        else:
            components = self.biconnected_components_subgraphs()
        h = None
        for component in components:
            if component.has_edge(edge):
                h = component
                break
        else:
            # edge connects two strongly connected components, so
            # no simple cycle starting with edge exists.
            return
    else:
        h = copy(self)
    # delete edge
    h.delete_edge(edge)

    by_weight, weight_function = self._get_weight_function(by_weight=by_weight,
                                                           weight_function=weight_function,
                                                           check_weight=check_weight)

    if by_weight:
        for e in self.edge_iterator():
            if weight_function(e) < 0:
                raise ValueError("negative weight is not allowed")

    it = h.shortest_simple_paths(source=edge[1], target=edge[0],
                                 weight_function=weight_function,
                                 by_weight=by_weight,
                                 check_weight=check_weight,
                                 algorithm=('Feng' if h.is_directed() else 'Yen'),
                                 report_edges=False,
                                 report_weight=True)

    edge_weight = weight_function(edge)

    if max_length is None:
        from sage.rings.infinity import Infinity
        max_length = Infinity
    for length, path in it:
        if length + edge_weight > max_length:
            break

        if report_weight:
            yield length + edge_weight, [edge[0]] + path
        else:
            yield [edge[0]] + path


def all_cycles_iterator(self, starting_vertices=None, simple=False,
                        rooted=False, max_length=None, trivial=False,
                        weight_function=None, by_weight=False,
                        check_weight=True, report_weight=False,
                        algorithm='A'):
    r"""
    Return an iterator over all the cycles of ``self`` starting with one of
    the given vertices. Each edge must have a positive weight.

    The cycles are enumerated in increasing length order. Here, a cycle
    means a closed walk `v_0e_0v_1e_1 \ldots v_{k-1}e_{k-1}v_0` where
    `v_i` are vertices and `e_i` are edges. `v_i` and `e_i`
    (`0 \leq i < k`) may not be distinct. If ``simple=True``, then cycles
    with distinct `v_i` and `e_i` (`0 \leq i < k`) are enumerated.

    INPUT:

    - ``starting_vertices`` -- iterable (default: ``None``); vertices from
      which the cycles must start. If ``None``, then all vertices of the
      graph can be starting points. This argument is necessary if ``rooted``
      is set to ``True``.

    - ``simple`` -- boolean (default: ``False``); if set to ``True``, then
      only simple cycles are considered. A cycle is simple if the only
      vertex occurring twice in it is the starting and ending one and no edge
      occurs twice.

    - ``rooted`` -- boolean (default: ``False``); if set to False, then
      cycles differing only by their starting vertex are considered the same
      (e.g. ``['a', 'b', 'c', 'a']`` and ``['b', 'c', 'a',
      'b']``). Otherwise, all cycles are enumerated.

    - ``max_length`` -- nonnegative integer (default: ``None``); the
      maximum length of the enumerated paths. If set to ``None``, then all
      lengths are allowed.

    - ``trivial`` -- boolean (default: ``False``); if set to ``True``, then
      the empty paths are also enumerated

    - ``weight_function`` -- function (default: ``None``); a function that
      takes as input an edge ``(u, v, l)`` and outputs its weight. If not
      ``None``, ``by_weight`` is automatically set to ``True``. If ``None``
      and ``by_weight`` is ``True``, we use the edge label ``l`` as a
      weight.

    - ``by_weight`` -- boolean (default: ``False``); if ``True``, the edges
      in the graph are weighted, otherwise all edges have weight 1

    - ``check_weight`` -- boolean (default: ``True``); whether to check that
      the ``weight_function`` outputs a number for each edge

    - ``report_weight`` -- boolean (default: ``False``); if ``False``, just
      a cycle is returned. Otherwise a tuple of cycle length and cycle is
      returned.

    - ``algorithm`` -- string (default: ``'A'``); the algorithm used to
      enumerate the cycles.

        - The algorithm ``'A'`` holds cycle iterators starting with each vertex,
          and output them in increasing length order.

        - The algorithm ``'B'`` holds cycle iterators starting with each edge,
          and output them in increasing length order. It depends on the k-shortest
          simple paths algorithm. Thus, it is not available if ``simple=False``.

    OUTPUT: iterator

    .. SEEALSO::

        - :meth:`all_simple_cycles`
        - :meth:`~sage.graphs.path_enumeration.shortest_simple_paths`

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

    A cycle is enumerated in increasing length order for a weighted graph::

        sage: g = DiGraph()
        sage: g.add_edges([('a', 'b', 2), ('a', 'e', 2), ('b', 'a', 1), ('b', 'c', 1),
        ....:              ('c', 'd', 2), ('d', 'b', 1), ('e', 'a', 2)])
        sage: it = g.all_cycles_iterator(simple=False, max_length=None,
        ....:                            by_weight=True, report_weight=True)
        sage: for i in range(9): print(next(it))
        (3, ['a', 'b', 'a'])
        (4, ['a', 'e', 'a'])
        (4, ['b', 'c', 'd', 'b'])
        (6, ['a', 'b', 'a', 'b', 'a'])
        (7, ['a', 'b', 'a', 'e', 'a'])
        (7, ['a', 'b', 'c', 'd', 'b', 'a'])
        (7, ['a', 'e', 'a', 'b', 'a'])
        (8, ['a', 'e', 'a', 'e', 'a'])
        (8, ['b', 'c', 'd', 'b', 'c', 'd', 'b'])

    Each edge must have a positive weight::

        sage: g = DiGraph()
        sage: g.add_edges([('a', 'b', -2), ('b', 'a', 1)])
        sage: next(g.all_cycles_iterator(simple=False, max_length=None,
        ....:                            by_weight=True, report_weight=True))
        Traceback (most recent call last):
        ...
        ValueError: negative weight is not allowed

    The algorithm ``'B'`` works for undirected graphs as well::

        sage: g = Graph({0: [1, 2], 1: [0, 2], 2: [0, 1, 3, 5], 3: [2, 4], 4: [3, 5], 5: [4, 2]})
        sage: for cycle in g.all_cycles_iterator(algorithm='B', simple=True, by_weight=True):
        ....:     print(cycle)
        [0, 1, 2, 0]
        [2, 3, 4, 5, 2]

    The algorithm ``'B'`` is available only when ``simple=True``::

        sage: g = DiGraph()
        sage: g.add_edges([('a', 'b', 1), ('b', 'a', 1)])
        sage: next(g.all_cycles_iterator(algorithm='B', simple=False))
        ....:
        Traceback (most recent call last):
        ...
        ValueError: The algorithm 'B' is available only when simple=True.

    The algorithm ``'A'`` works for undirected graphs as well. Specifically, each cycle is
    enumerated exactly once, meaning a cycle and its reverse are not listed separately::

        sage: g = Graph({0: [1, 2], 1: [0, 2], 2: [0, 1]})
        sage: for cycle in g.all_cycles_iterator(algorithm='A', simple=True):
        ....:     print(cycle)
        [0, 1, 2, 0]
    """
    if starting_vertices is None:
        starting_vertices = self

    if algorithm == 'B' and not simple:
        raise ValueError("The algorithm 'B' is available only when simple=True.")

    by_weight, weight_function = self._get_weight_function(by_weight=by_weight,
                                                           weight_function=weight_function,
                                                           check_weight=check_weight)

    if by_weight:
        for e in self.edge_iterator():
            if weight_function(e) < 0:
                raise ValueError("negative weight is not allowed")

    if algorithm == 'A':
        if self.is_directed():
            # Since a cycle is always included in a given strongly connected
            # component, we may remove edges from the graph
            sccs = self.strongly_connected_components()
            if len(sccs) == 1:
                h = self
            else:
                d = {}
                for id, component in enumerate(sccs):
                    for v in component:
                        d[v] = id
                h = copy(self)
                h.delete_edges((u, v) for u, v in h.edge_iterator(labels=False) if d[u] != d[v])
        else:
            h = self

        # We create one cycles iterator per vertex. This is necessary if we
        # want to iterate over cycles with increasing length.
        def cycle_iter(v):
            return h._all_cycles_iterator_vertex(v,
                                                 starting_vertices=starting_vertices,
                                                 simple=simple,
                                                 rooted=rooted,
                                                 max_length=max_length,
                                                 trivial=trivial,
                                                 remove_acyclic_edges=False,
                                                 weight_function=weight_function,
                                                 by_weight=by_weight,
                                                 check_weight=check_weight,
                                                 report_weight=True)

        iterators = {v: cycle_iter(v) for v in starting_vertices}
    elif algorithm == 'B':
        def simple_cycle_iter(hh, e):
            return hh._all_simple_cycles_iterator_edge(e,
                                                       max_length=max_length,
                                                       remove_unnecessary_edges=False,
                                                       weight_function=weight_function,
                                                       by_weight=by_weight,
                                                       check_weight=check_weight,
                                                       report_weight=True)
        if self.is_directed():
            def decompose(hh):
                return hh.strongly_connected_components_subgraphs()
        else:
            def decompose(hh):
                return hh.biconnected_components_subgraphs()

        components = decompose(self)
        iterators = dict()
        while components:
            hh = components.pop()
            if not hh.size():
                continue
            e = next(hh.edge_iterator())
            iterators[e] = simple_cycle_iter(hh, e)
            hh.delete_edge(e)
            components.extend(decompose(hh))
    else:
        raise ValueError(f"The algorithm {algorithm} is not valid. \
                            Use the algorithm 'A' or 'B'.")

    cycles = []
    for key, it in iterators.items():
        try:
            length, cycle = next(it)
            cycles.append((length, cycle, key))
        except StopIteration:
            pass
    # Since we always extract a shortest path, using a heap
    # can speed up the algorithm
    from heapq import heapify, heappop, heappush
    heapify(cycles)
    while cycles:
        # We choose the shortest available cycle
        length, shortest_cycle, key = heappop(cycles)
        if report_weight:
            yield (length, shortest_cycle)
        else:
            yield shortest_cycle
        # We update the cycle iterator to its next available cycle if it
        # exists
        try:
            length, cycle = next(iterators[key])
            heappush(cycles, (length, cycle, key))
        except StopIteration:
            pass


def all_simple_cycles(self, starting_vertices=None, rooted=False,
                      max_length=None, trivial=False,
                      weight_function=None, by_weight=False,
                      check_weight=True, report_weight=False,
                      algorithm='A'):
    r"""
    Return a list of all simple cycles of ``self``. The cycles are
    enumerated in increasing length order. Each edge must have a
    positive weight.

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

    - ``weight_function`` -- function (default: ``None``); a function that
      takes as input an edge ``(u, v, l)`` and outputs its weight. If not
      ``None``, ``by_weight`` is automatically set to ``True``. If ``None``
      and ``by_weight`` is ``True``, we use the edge label ``l`` as a
      weight.

    - ``by_weight`` -- boolean (default: ``False``); if ``True``, the edges
      in the graph are weighted, otherwise all edges have weight 1

    - ``check_weight`` -- boolean (default: ``True``); whether to check that
      the ``weight_function`` outputs a number for each edge

    - ``report_weight`` -- boolean (default: ``False``); if ``False``, just
      a cycle is returned. Otherwise a tuple of cycle length and cycle is
      returned.

    - ``algorithm`` -- string (default: ``'A'``); the algorithm used to
      enumerate the cycles.

        - The algorithm ``'A'`` holds cycle iterators starting with each vertex,
          and output them in increasing length order.

        - The algorithm ``'B'`` holds cycle iterators starting with each edge,
          and output them in increasing length order.

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
        [[0, 1, 0], [0, 2, 0], [0, 3, 0], [0, 4, 0], [0, 5, 0], [0, 6, 0], [0, 7, 0],
         [0, 8, 0], [0, 9, 0], [0, 10, 0], [0, 11, 0], [0, 12, 0], [0, 13, 0],
         [0, 14, 0], [0, 15, 0], [0, 16, 0], [0, 17, 0], [0, 18, 0], [0, 19, 0],
         [1, 2, 1], [1, 3, 1], [1, 4, 1], [1, 5, 1], [1, 6, 1], [1, 7, 1], [1, 8, 1],
         [1, 9, 1], [1, 10, 1], [1, 11, 1], [1, 12, 1], [1, 13, 1], [1, 14, 1],
         [1, 15, 1], [1, 16, 1], [1, 17, 1], [1, 18, 1], [1, 19, 1], [2, 3, 2],
         [2, 4, 2], [2, 5, 2], [2, 6, 2], [2, 7, 2], [2, 8, 2], [2, 9, 2], [2, 10, 2],
         [2, 11, 2], [2, 12, 2], [2, 13, 2], [2, 14, 2], [2, 15, 2], [2, 16, 2],
         [2, 17, 2], [2, 18, 2], [2, 19, 2], [3, 4, 3], [3, 5, 3], [3, 6, 3],
         [3, 7, 3], [3, 8, 3], [3, 9, 3], [3, 10, 3], [3, 11, 3], [3, 12, 3],
         [3, 13, 3], [3, 14, 3], [3, 15, 3], [3, 16, 3], [3, 17, 3], [3, 18, 3],
         [3, 19, 3], [4, 5, 4], [4, 6, 4], [4, 7, 4], [4, 8, 4], [4, 9, 4], [4, 10, 4],
         [4, 11, 4], [4, 12, 4], [4, 13, 4], [4, 14, 4], [4, 15, 4], [4, 16, 4],
         [4, 17, 4], [4, 18, 4], [4, 19, 4], [5, 6, 5], [5, 7, 5], [5, 8, 5],
         [5, 9, 5], [5, 10, 5], [5, 11, 5], [5, 12, 5], [5, 13, 5], [5, 14, 5],
         [5, 15, 5], [5, 16, 5], [5, 17, 5], [5, 18, 5], [5, 19, 5], [6, 7, 6],
         [6, 8, 6], [6, 9, 6], [6, 10, 6], [6, 11, 6], [6, 12, 6], [6, 13, 6],
         [6, 14, 6], [6, 15, 6], [6, 16, 6], [6, 17, 6], [6, 18, 6], [6, 19, 6],
         [7, 8, 7], [7, 9, 7], [7, 10, 7], [7, 11, 7], [7, 12, 7], [7, 13, 7],
         [7, 14, 7], [7, 15, 7], [7, 16, 7], [7, 17, 7], [7, 18, 7], [7, 19, 7],
         [8, 9, 8], [8, 10, 8], [8, 11, 8], [8, 12, 8], [8, 13, 8], [8, 14, 8],
         [8, 15, 8], [8, 16, 8], [8, 17, 8], [8, 18, 8], [8, 19, 8], [9, 10, 9],
         [9, 11, 9], [9, 12, 9], [9, 13, 9], [9, 14, 9], [9, 15, 9], [9, 16, 9],
         [9, 17, 9], [9, 18, 9], [9, 19, 9], [10, 11, 10], [10, 12, 10], [10, 13, 10],
         [10, 14, 10], [10, 15, 10], [10, 16, 10], [10, 17, 10], [10, 18, 10],
         [10, 19, 10], [11, 12, 11], [11, 13, 11], [11, 14, 11], [11, 15, 11],
         [11, 16, 11], [11, 17, 11], [11, 18, 11], [11, 19, 11], [12, 13, 12],
         [12, 14, 12], [12, 15, 12], [12, 16, 12], [12, 17, 12], [12, 18, 12],
         [12, 19, 12], [13, 14, 13], [13, 15, 13], [13, 16, 13], [13, 17, 13],
         [13, 18, 13], [13, 19, 13], [14, 15, 14], [14, 16, 14], [14, 17, 14],
         [14, 18, 14], [14, 19, 14], [15, 16, 15], [15, 17, 15], [15, 18, 15],
         [15, 19, 15], [16, 17, 16], [16, 18, 16], [16, 19, 16], [17, 18, 17],
         [17, 19, 17], [18, 19, 18]]

        sage: g = graphs.CompleteGraph(20).to_directed()
        sage: g.all_simple_cycles(max_length=2, starting_vertices=[0])
        [[0, 1, 0], [0, 2, 0], [0, 3, 0], [0, 4, 0], [0, 5, 0],
         [0, 6, 0], [0, 7, 0], [0, 8, 0], [0, 9, 0], [0, 10, 0],
         [0, 11, 0], [0, 12, 0], [0, 13, 0], [0, 14, 0], [0, 15, 0],
         [0, 16, 0], [0, 17, 0], [0, 18, 0], [0, 19, 0]]

    One may prefer to distinguish equivalent cycles having distinct starting
    vertices (compare the following examples)::

        sage: g = graphs.CompleteGraph(4).to_directed()
        sage: g.all_simple_cycles(max_length=2, rooted=False)
        [[0, 1, 0], [0, 2, 0], [0, 3, 0], [1, 2, 1], [1, 3, 1], [2, 3, 2]]
        sage: g.all_simple_cycles(max_length=2, rooted=True)
        [[0, 1, 0], [0, 2, 0], [0, 3, 0], [1, 0, 1], [1, 2, 1], [1, 3, 1],
         [2, 0, 2], [2, 1, 2], [2, 3, 2], [3, 0, 3], [3, 1, 3], [3, 2, 3]]

    A cycle is enumerated in increasing length order for a weighted graph::

        sage: cycles = g.all_simple_cycles(weight_function=lambda e:e[0]+e[1],
        ....:                              by_weight=True, report_weight=True)
        sage: cycles
        [(2, [0, 1, 0]), (4, [0, 2, 0]), (6, [0, 1, 2, 0]), (6, [0, 2, 1, 0]),
         (6, [0, 3, 0]), (6, [1, 2, 1]), (8, [0, 1, 3, 0]), (8, [0, 3, 1, 0]),
         (8, [1, 3, 1]), (10, [0, 2, 3, 0]), (10, [0, 3, 2, 0]), (10, [2, 3, 2]),
         (12, [0, 1, 2, 3, 0]), (12, [0, 1, 3, 2, 0]), (12, [0, 2, 1, 3, 0]),
         (12, [0, 2, 3, 1, 0]), (12, [0, 3, 1, 2, 0]), (12, [0, 3, 2, 1, 0]),
         (12, [1, 2, 3, 1]), (12, [1, 3, 2, 1])]

    The algorithm ``'B'`` can be used::

        sage: cycles_B = g.all_simple_cycles(weight_function=lambda e:e[0]+e[1], by_weight=True,
        ....:                                report_weight=True, algorithm='B')
        sage: cycles_B
        [(2.0, [0, 1, 0]), (4.0, [0, 2, 0]), (6.0, [0, 1, 2, 0]), (6.0, [0, 2, 1, 0]),
         (6.0, [0, 3, 0]), (6.0, [1, 2, 1]), (8.0, [0, 1, 3, 0]), (8.0, [0, 3, 1, 0]),
         (8.0, [1, 3, 1]), (10.0, [0, 2, 3, 0]), (10.0, [0, 3, 2, 0]), (10.0, [2, 3, 2]),
         (12.0, [0, 1, 3, 2, 0]), (12.0, [0, 1, 2, 3, 0]), (12.0, [0, 2, 3, 1, 0]),
         (12.0, [0, 2, 1, 3, 0]), (12.0, [0, 3, 2, 1, 0]), (12.0, [0, 3, 1, 2, 0]),
         (12.0, [1, 2, 3, 1]), (12.0, [1, 3, 2, 1])]
        sage: cycles.sort() == cycles_B.sort()
        True

    The algorithm ``'A'`` is available for undirected graphs. Specifically, each cycle is
    enumerated exactly once, meaning a cycle and its reverse are not listed separately::

        sage: g = Graph({0: [1, 2], 1: [0, 2], 2: [0, 1]})
        sage: g.all_simple_cycles(algorithm='A')
        [[0, 1, 2, 0]]
    """
    return list(self.all_cycles_iterator(starting_vertices=starting_vertices,
                                         simple=True, rooted=rooted,
                                         max_length=max_length, trivial=trivial,
                                         weight_function=weight_function,
                                         by_weight=by_weight,
                                         check_weight=check_weight,
                                         report_weight=report_weight,
                                         algorithm=algorithm))
