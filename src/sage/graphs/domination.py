# -*- coding: utf-8 -*-
r"""
Domination

This module implements methods related to the notion of domination in graphs,
and more precisely:

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~dominating_set` | Return a minimum distance-`k` dominating set of the graph.
    :meth:`~minimal_dominating_sets` | Return an iterator over the minimal dominating sets of a graph.
    :meth:`~is_dominating` | Check whether a set of vertices dominates a graph.
    :meth:`~is_redundant` | Check whether a set of vertices has redundant vertices (with respect to domination).
    :meth:`~private_neighbors` | Return the private neighbors of a vertex with respect to other vertices.
    :meth:`~greedy_dominating_set` | Return a greedy distance-`k` dominating set of of the graph.


EXAMPLES:

We compute the size of a minimum dominating set of the Petersen graph::

    sage: g = graphs.PetersenGraph()
    sage: g.dominating_set(value_only=True)
    3

We enumerate the minimal dominating sets of the 5-star graph::

    sage: g = graphs.StarGraph(5)
    sage: list(g.minimal_dominating_sets())
    [{0}, {1, 2, 3, 4, 5}]

Now only those that dominate the middle vertex::

    sage: list(g.minimal_dominating_sets([0]))
    [{0}, {1}, {2}, {3}, {4}, {5}]

Now the minimal dominating sets of the 5-path graph::

    sage: g = graphs.PathGraph(5)
    sage: list(g.minimal_dominating_sets())
    [{0, 2, 4}, {1, 4}, {0, 3}, {1, 3}]

We count the minimal dominating sets of the Petersen graph::

    sage: sum(1 for _ in graphs.PetersenGraph().minimal_dominating_sets())
    27

Methods
-------
"""

# ****************************************************************************
#       Copyright (C) 2009 Nathann Cohen <nathann.cohen@gmail.com>
#                     2019 Jean-Florent Raymond  <j-florent.raymond@uca.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from copy import copy
from sage.rings.integer import Integer

def is_dominating(G, dom, focus=None):
    r"""
    Check whether ``dom`` is a dominating set of ``G``.

    We say that a set `D` of vertices of a graph `G` dominates a set `S` if
    every vertex of `S` either belongs to `D` or is adjacent to a vertex of `D`.
    Also, `D` is a dominating set of `G` if it dominates `V(G)`.

    INPUT:

    - ``dom`` -- iterable of vertices of ``G``; the vertices of the supposed
      dominating set.

    - ``focus`` -- iterable of vertices of ``G`` (default: ``None``); if
      specified, this method checks instead if ``dom`` dominates the vertices in
      ``focus``.

    EXAMPLES::

        sage: g = graphs.CycleGraph(5)
        sage: g.is_dominating([0,1], [4, 2])
        True

        sage: g.is_dominating([0,1])
        False
    """
    to_dom = set(G) if focus is None else set(focus)

    for v in dom:
        if not to_dom:
            return True
        to_dom.difference_update(G.neighbor_iterator(v, closed=True))

    return not to_dom

def is_redundant(G, dom, focus=None):
    r"""
    Check whether ``dom`` has redundant vertices.

    For a graph `G` and sets `D` and `S` of vertices, we say that a vertex `v
    \in D` is *redundant* in `S` if `v` has no private neighbor with respect to
    `D` in `S`.  In other words, there is no vertex in `S` that is dominated by
    `v` but not by `D \setminus \{v\}`.

    INPUT:

    - ``dom`` -- iterable of vertices of ``G``; where we look for redundant
      vertices.

    - ``focus`` -- iterable of vertices of ``G`` (default: ``None``); if
      specified, this method checks instead whether ``dom`` has a redundant
      vertex in ``focus``.

    .. WARNING::

        The assumption is made that ``focus`` (if provided) does not contain
        repeated vertices.

    EXAMPLES::

        sage: G = graphs.CubeGraph(3)
        sage: G.is_redundant(['000', '101'], ['011'])
        True
        sage: G.is_redundant(['000', '101'])
        False
    """
    dom = list(dom)
    focus = list(G) if focus is None else list(focus)

    # dominator[v] (for v in focus) will be equal to:
    #  - (0, None) if v has no neighbor in dom
    #  - (1, u) if v has a unique neighbor in dom, u
    #  - (2, None) if v has >= 2 neighbors in dom

    # Initialization
    dominator = {v: (0, None) for v in focus}

    for x in dom:
        for v in G.neighbor_iterator(x, closed=True):
            if v in focus:
                # remember about x only if we never encountered
                # neighbors of v so far
                if dominator[v][0] == 0:
                    dominator[v] = (1, x)
                elif dominator[v][0] == 1:
                    dominator[v] = (2, None)

    # We now compute the subset of vertices of dom that have a private neighbor
    # in focus. A vertex v in dom has a private neighbor in focus if it is
    # dominated by a unique vertex, that is if dominator[v][0] == 1
    with_private = set()
    for v in focus:
        if dominator[v][0] == 1:
            with_private.add(dominator[v][1])

    # By construction with_private is a subset of dom and we assume the elements
    # of dom to be unique, so the following is equivalent to checking
    # with_private != set(dom)
    return len(with_private) != len(dom)

def private_neighbors(G, vertex, dom):
    r"""
    Return the private neighbors of a vertex with respect to other vertices.

    A private neighbor of a vertex `v` with respect to a vertex subset `D`
    is a closed neighbor of `v` that is not dominated by a vertex of `D
    \setminus \{v\}`.

    INPUT:

    - ``vertex`` -- a vertex of ``G``.

    - ``dom`` -- iterable of vertices of ``G``; the vertices possibly stealing
      private neighbors from ``vertex``.

    OUTPUT:

    Return the closed neighbors of ``vertex`` that are not closed neighbors
    of any other vertex of ``dom``.

    EXAMPLES::

        sage: g = graphs.PathGraph(5)
        sage: list(g.private_neighbors(1, [1, 3, 4]))
        [1, 0]

        sage: list(g.private_neighbors(1, [3, 4]))
        [1, 0]

        sage: list(g.private_neighbors(1, [3, 4, 0]))
        []
    """
    # The set of all vertices that are dominated by vertex_subset - vertex:
    closed_neighborhood_vs = set()
    for u in dom:
        if u != vertex:
            closed_neighborhood_vs.update(
                G.neighbor_iterator(u, closed=True))

    return (neighbor
            for neighbor in G.neighbor_iterator(vertex, closed=True)
            if neighbor not in closed_neighborhood_vs)


# ==============================================================================
# Computation of minimum dominating
# ==============================================================================

def dominating_set(g, k=1, independent=False, total=False, value_only=False,
                   solver=None, verbose=0, *, integrality_tolerance=1e-3):
    r"""
    Return a minimum distance-`k` dominating set of the graph.

    A minimum dominating set `S` of a graph `G` is a set of its vertices of
    minimal cardinality such that any vertex of `G` is in `S` or has one of its
    neighbors in `S`. See the :wikipedia:`Dominating_set`.

    A minimum distance-`k` dominating set is a set `S` of vertices of `G` of
    minimal cardinality such that any vertex of `G` is in `S` or at distance at
    most `k` from a vertex in `S`. A distance-`0` dominating set is the set of
    vertices itself, and when `k` is the radius of the graph, any vertex
    dominates all the other vertices.

    As an optimization problem, it can be expressed as follows, where `N^k(u)`
    denotes the set of vertices at distance at most `k` from `u` (the set of
    neighbors when `k=1`):

    .. MATH::

        \mbox{Minimize : }&\sum_{v\in G} b_v\\
        \mbox{Such that : }&\forall v \in G, b_v+\sum_{u \in N^k(v)} b_u\geq 1\\
        &\forall x\in G, b_x\mbox{ is a binary variable}

    INPUT:

    - ``k`` -- a non-negative integer (default: ``1``); the domination distance

    - ``independent`` -- boolean (default: ``False``); when ``True``, computes a
      minimum independent dominating set, that is a minimum dominating set that
      is also an independent set (see also
      :meth:`~sage.graphs.graph.independent_set`)

    - ``total`` -- boolean (default: ``False``); when ``True``, computes a total
      dominating set (see the See the :wikipedia:`Dominating_set`)

    - ``value_only`` -- boolean (default: ``False``); whether to only return the
      cardinality of the computed dominating set, or to return its list of
      vertices (default)

    - ``solver`` -- string (default: ``None``); specify a Mixed Integer Linear
      Programming (MILP) solver to be used. If set to ``None``, the default one
      is used. For more information on MILP solvers and which default solver is
      used, see the method :meth:`solve
      <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
      :class:`MixedIntegerLinearProgram
      <sage.numerical.mip.MixedIntegerLinearProgram>`.

    - ``verbose`` -- integer (default: ``0``); sets the level of verbosity. Set
      to 0 by default, which means quiet.

    - ``integrality_tolerance`` -- float; parameter for use with MILP solvers
      over an inexact base ring; see
      :meth:`MixedIntegerLinearProgram.get_values`.

    EXAMPLES:

    A basic illustration on a ``PappusGraph``::

        sage: g = graphs.PappusGraph()
        sage: g.dominating_set(value_only=True)
        5

    If we build a graph from two disjoint stars, then link their centers we will
    find a difference between the cardinality of an independent set and a stable
    independent set::

        sage: g = 2 * graphs.StarGraph(5)
        sage: g.add_edge(0, 6)
        sage: len(g.dominating_set())
        2
        sage: len(g.dominating_set(independent=True))
        6

    The total dominating set of the Petersen graph has cardinality 4::

        sage: G = graphs.PetersenGraph()
        sage: G.dominating_set(total=True, value_only=True)
        4

    The dominating set is calculated for both the directed and undirected graphs
    (modification introduced in :trac:`17905`)::

        sage: g = digraphs.Path(3)
        sage: g.dominating_set(value_only=True)
        2
        sage: g = graphs.PathGraph(3)
        sage: g.dominating_set(value_only=True)
        1

    Cardinality of distance-`k` dominating sets::

        sage: G = graphs.PetersenGraph()
        sage: [G.dominating_set(k=k, value_only=True) for k in range(G.radius() + 1)]
        [10, 3, 1]
        sage: G = graphs.PathGraph(5)
        sage: [G.dominating_set(k=k, value_only=True) for k in range(G.radius() + 1)]
        [5, 2, 1]
    """
    g._scream_if_not_simple(allow_multiple_edges=True, allow_loops=not total)

    if not k:
        return g.order() if value_only else list(g)
    elif k < 0:
        raise ValueError("the domination distance must be a non-negative integer")

    from sage.numerical.mip import MixedIntegerLinearProgram
    p = MixedIntegerLinearProgram(maximization=False, solver=solver)
    b = p.new_variable(binary=True)

    if k == 1:
        # For any vertex v, one of its neighbors or v itself is in the minimum
        # dominating set. If g is directed, we use the in neighbors of v
        # instead.
        neighbors_iter = g.neighbor_in_iterator if g.is_directed() else g.neighbor_iterator
    else:
        # When k > 1, we use BFS to determine the vertices that can reach v
        # through a path of length at most k
        gg = g.reverse() if g.is_directed() else g

        def neighbors_iter(x):
            it = gg.breadth_first_search(x, distance=k)
            _ = next(it)
            yield from it

    if total:
        # We want a total dominating set
        for v in g:
            p.add_constraint(p.sum(b[u] for u in neighbors_iter(v)), min=1)
    else:
        for v in g:
            p.add_constraint(b[v] + p.sum(b[u] for u in neighbors_iter(v)), min=1)

    if independent:
        # no two adjacent vertices are in the set
        for u, v in g.edge_iterator(labels=None):
            p.add_constraint(b[u] + b[v], max=1)

    # Minimizes the number of vertices used
    p.set_objective(p.sum(b[v] for v in g))

    p.solve(log=verbose)
    b = p.get_values(b, convert=bool, tolerance=integrality_tolerance)
    dom = [v for v in g if b[v]]
    return Integer(len(dom)) if value_only else dom

# ==============================================================================
# Enumeration of minimal dominating set as described in [BDHPR2019]_
# ==============================================================================

def _parent(G, dom, V_prev):
    r"""
    Return a subset of dom that is irredundant in ``V_prev``.

    For internal use.

    INPUT:

    - ``G`` -- a graph

    - ``dom`` -- an iterable of vertices of ``G``

    - ``V_prev`` -- an iterable of vertices of ``G``

    OUTPUT:

    Return the list obtained from ``dom`` by iteratively removing those
    vertices of minimum index that have no private neighbor in ``V_prev``.

    TESTS::

        sage: from sage.graphs.domination import _parent
        sage: G = graphs.PathGraph(4)
        sage: G.add_vertices([4, 5])
        sage: G.add_edges([(4, 1), (5, 2)])
        sage: _parent(G, [0, 2, 4, 5], [1, 2])
        [4, 5]
        sage: _parent(G, [0, 2, 4, 5], [1, 3])
        [2]
    """
    # The list where we search vertices
    D_start = sorted(dom, reverse=True)

    # The list to be output at the end, that we construct:
    D_end = []

    while D_start:
        v = D_start.pop()  # element of min index
        priv = set(G.neighbor_iterator(v, closed=True))
        # We remove the vertices already dominated
        # by other vertices of (D_end union D_start)
        priv.difference_update(*(G.neighbor_iterator(u, closed=True)
                                 for u in D_start if u != v))
        priv.difference_update(*(G.neighbor_iterator(u, closed=True)
                                 for u in D_end if u != v))
        # Now priv is the private neighborhood of v in G wrt D_start + D_end
        if priv.intersection(V_prev) != set():
            # if v has a private in V_prev, we keep it
            D_end.append(v)

    return D_end


def _peel(G, A):
    r"""
    Return a peeling of a vertex iterable of a graph.

    For internal use.
    Given a graph `G` and a subset `A` of its vertices, a peeling of `(G,A)` is
    a list `[(u_0, V_0), \dots, (u_p, V_p)]` such that `u_0` is ``None``,
    `V_0` is the empty set, `V_p = A` and for every `i \in \{1, \dots, p\}`,
    `V_{i-1} = V_i \setminus N[v_i]`, for some vertex `u_i` of `V_i`.

    INPUT:

    - ``G`` -- a graph

    - ``A`` -- a set of vertices of `G`

    OUTPUT:

    A peeling of `(G, A)`.

    TESTS::

        sage: from sage.graphs.domination import _peel
        sage: G = Graph(10); _peel(G, {0, 1, 2, 3, 4})
        [(None, set()),
        (4, {4}),
        (3, {3, 4}),
        (2, {2, 3, 4}),
        (1, {1, 2, 3, 4}),
        (0, {0, 1, 2, 3, 4})]


        sage: from sage.graphs.domination import _peel
        sage: G = graphs.PathGraph(10); _peel(G, set((i for i in range(10) if i%2==0)))
        [(None, set()),
        (8, {8}),
        (6, {6, 8}),
        (4, {4, 6, 8}),
        (2, {2, 4, 6, 8}),
        (0, {0, 2, 4, 6, 8})]

    """
    Acomp = set(G)
    Acomp.difference_update(A)  # Acomp  = V - A

    peeling = []
    H = copy(G)
    H.delete_vertices(list(Acomp))
    del Acomp

    while H:
        ui = next(H.vertex_iterator())  # pick some vertex of H
        Vi = set(H)
        peeling.append((ui, Vi))
        H.delete_vertices(H.neighbor_iterator(ui, closed=True))
    peeling.append((None, set()))
    peeling.reverse()
    return peeling


def _cand_ext_enum(G, to_dom, u_next):
    r"""
    Return the minimal dominating sets of ``to_dom``.

    For internal use.
    Assumption: ``u_next`` dominates ``to_dom``.

    INPUT:

    - ``G`` -- a graph

    - ``to_dom`` -- a ``set()`` of vertices of ``G``

    - ``u_next`` -- a vertex of ``G`` that dominates ``to_dom``

    OUTPUT:

    An iterator over the minimal dominating sets of ``to_dom``.

    TESTS::

        sage: from sage.graphs.domination import _cand_ext_enum
        sage: g = graphs.DiamondGraph()
        sage: l = list(_cand_ext_enum(g, {0, 2, 3}, 1,))
        sage: len(l) == 3 and {1} in l and {2} in l and {0, 3} in l
        True
    """

    def _aux_with_rep(H, to_dom, u_next):
        """
        Return the minimal dominating sets of ``to_dom``, with the
        assumption that ``u_next`` dominates ``to_som``.

        .. WARNING::

            The same output may be output several times (up to |H| times).

        In order to later remove duplicates, we here output pairs ``(ext, i)``
        where ``ext`` is the output candidate extension and ``i`` counts how
        many elements have already been output.
        """
        if u_next not in to_dom:
            # In this case, enumerating the minimal DSs of the subset
            # to_dom is a smaller instance as it excludes u_next:

            cand_ext_index = 0

            for ext in H.minimal_dominating_sets(to_dom):
                yield (ext, cand_ext_index)
                cand_ext_index += 1

        elif to_dom == {u_next}:
            # In this case, only u_next has to be dominated
            cand_ext_index = 0
            for w in H.neighbor_iterator(u_next, closed=True):
                # Notice that the case w = u_next is included
                yield ({w}, cand_ext_index)
                cand_ext_index += 1

        else:
            # In this case, both u_next and to_dom-u_next have to be dominated

            # We first output the trivial output
            # (as to_dom is subset of N(u_next)):
            yield ({u_next}, 0)
            # Start from 1 because we already output the 0-th elt:
            cand_ext_index = 1

            # When u_next is not in the DS, one of its neighbors w should be:
            for w in H.neighbor_iterator(u_next):

                remains_to_dom = set(to_dom)
                remains_to_dom.difference_update(
                    H.neighbor_iterator(w, closed=True))
                # Here again we recurse on a smaller instance at it
                # excludes u_next (and w)
                for Q in H.minimal_dominating_sets(remains_to_dom):
                    ext = set(Q)
                    ext.add(w)
                    # By construction w dominates u_next and Q dominates
                    # to_dom - N[w], so ext dominates to_dom: it is a
                    # valid output iff it is not redundant
                    if not H.is_redundant(ext):
                        yield (ext, cand_ext_index)
                        cand_ext_index += 1
    #
    # End of aux_with_rep routine

    # Here we use aux_with_rep twice to enumerate the minimal
    # dominating sets while avoiding repeated outputs
    for (X, i) in _aux_with_rep(G, to_dom, u_next):
        for (Y, j) in _aux_with_rep(G, to_dom, u_next):
            if j >= i:
                # This is the first time we meet X: we output it
                yield X
                break
            elif Y == X:  # These are sets
                # X has already been output in the past: we ignore it
                break


def minimal_dominating_sets(G, to_dominate=None, work_on_copy=True):
    r"""
    Return an iterator over the minimal dominating sets of a graph.

    INPUT:

    - ``G`` -- a graph.

    - ``to_dominate`` -- vertex iterable or ``None`` (default: ``None``);
      the set of vertices to be dominated.

    - ``work_on_copy`` -- boolean (default: ``True``); whether or not to work on
      a copy of the input graph; if set to ``False``, the input graph will be
      modified (relabeled).
    
    OUTPUT:

    An iterator over the inclusion-minimal sets of vertices of ``G``.
    If ``to_dominate`` is provided, return an iterator over the
    inclusion-minimal sets of vertices that dominate the vertices of
    ``to_dominate``.

    ALGORITHM: The algorithm described in [BDHPR2019]_.

    AUTHOR: Jean-Florent Raymond (2019-03-04) -- initial version.

    EXAMPLES::

        sage: G = graphs.ButterflyGraph()
        sage: ll = list(G.minimal_dominating_sets())
        sage: pp = [{0, 1}, {1, 3}, {0, 2}, {2, 3}, {4}]
        sage: len(ll) == len(pp) and all(x in pp for x in ll) and all(x in ll for x in pp)
        True

        sage: ll = list(G.minimal_dominating_sets([0,3]))
        sage: pp = [{0}, {3}, {4}]
        sage: len(ll) == len(pp) and all(x in pp for x in ll) and all(x in ll for x in pp)
        True

        sage: ll = list(G.minimal_dominating_sets([4]))
        sage: pp = [{4}, {0}, {1}, {2}, {3}]
        sage: len(ll) == len(pp) and all(x in pp for x in ll) and all(x in ll for x in pp)
        True

    ::

        sage: ll = list(graphs.PetersenGraph().minimal_dominating_sets())
        sage: pp = [{0, 2, 6},
        ....: {0, 9, 3},
        ....: {0, 8, 7},
        ....: {1, 3, 7},
        ....: {1, 4, 5},
        ....: {8, 1, 9},
        ....: {8, 2, 4},
        ....: {9, 2, 5},
        ....: {3, 5, 6},
        ....: {4, 6, 7},
        ....: {0, 8, 2, 9},
        ....: {0, 3, 6, 7},
        ....: {1, 3, 5, 9},
        ....: {8, 1, 4, 7},
        ....: {2, 4, 5, 6},
        ....: {0, 1, 2, 3, 4},
        ....: {0, 1, 2, 5, 7},
        ....: {0, 1, 4, 6, 9},
        ....: {0, 1, 5, 6, 8},
        ....: {0, 8, 3, 4, 5},
        ....: {0, 9, 4, 5, 7},
        ....: {8, 1, 2, 3, 6},
        ....: {1, 2, 9, 6, 7},
        ....: {9, 2, 3, 4, 7},
        ....: {8, 2, 3, 5, 7},
        ....: {8, 9, 3, 4, 6},
        ....: {8, 9, 5, 6, 7}]
        sage: len(ll) == len(pp) and all(x in pp for x in ll) and all(x in ll for x in pp)
        True

    TESTS:

    The empty graph is handled correctly::

        sage: list(Graph().minimal_dominating_sets())
        [set()]

    Test on all graphs on 6 vertices::

        sage: from sage.combinat.subset import Subsets
        sage: def minimal_dominating_sets_naive(G):
        ....:     return (S for S in Subsets(G.vertices())
        ....:             if not(G.is_redundant(S)) and G.is_dominating(S))
        sage: def big_check(n):
        ....:     for G in graphs(n):
        ....:         ll = list(G.minimal_dominating_sets())
        ....:         pp = list(minimal_dominating_sets_naive(G))
        ....:         if (len(pp) != len(pp)
        ....:             or any(x not in pp for x in ll)
        ....:             or any(x not in ll for x in pp)):
        ....:             return False
        ....:     return True
        sage: big_check(6)  # long time
        True

    Outputs are unique::

        sage: def check_uniqueness(g):
        ....:     counter_1 = 0
        ....:     for dom_1 in g.minimal_dominating_sets():
        ....:         counter_1 += 1
        ....:         counter_2 = 0
        ....:         for dom_2 in g.minimal_dominating_sets():
        ....:             counter_2 += 1
        ....:             if counter_2 >= counter_1:
        ....:                 break
        ....:             if dom_1 == dom_2:
        ....:                 return False
        ....:     return True
        sage: check_uniqueness(graphs.RandomGNP(9, 0.5))
        True
    """

    def tree_search(H, plng, dom, i):
        r"""
        Enumerate minimal dominating sets recursively.

        INPUT:

        - ``H`` -- a graph

        - ``plng`` -- a peeling of H (result of :func:`_peel`)

        - ``dom`` -- a minimal dominating set of ``plng[i][1]``

        - ``i`` -- an integer, the current position in ``plng``

        OUTPUT:

        An iterator over those minimal dominating sets (in ``H``) of
        ``plng[-1][1]`` that are children of ``dom`` (with respect to the
        :func:`parent` function).

        ALGORITHM:

        We iterate over those minimal dominating sets of ``plng[i + 1][1]`` that
        are children of dom and call recursively on each. The fact that we
        iterate over children (with respect to the `parent` function) ensures
        that we do not have repeated outputs.
        """
        if i == len(plng) - 1:
            # we reached a leaf, i.e. dom is a minimal dominating set
            # of plng[i][1] = plng[-1][1]
            yield dom
            return

        u_next, V_next = plng[i + 1]

        if H.is_dominating(dom, V_next):
            # if dom dominates V_next
            # then dom is its unique extension: we recurse on it
            for Di in tree_search(H, plng, dom, i + 1):
                yield Di
            return

        # Otherwise, V_next - <what dom dominates> is what we have to dominate
        to_dom = V_next - set().union(
            *(G.neighbor_iterator(vert, closed=True)
              for vert in dom))

        for can_ext in _cand_ext_enum(H, to_dom, u_next):

            # We complete dom with can_ext -> canD
            canD = set().union(can_ext, dom)

            if (not H.is_redundant(canD, V_next)) and set(dom) == set(_parent(H, canD, plng[i][1])):
                # By construction, can_ext is a dominating set of
                # `V_next - N[dom]`, so canD dominates V_next.
                # If canD is a legitimate child of dom and is not redundant, we
                # recurse on it:
                for Di in tree_search(H, plng, canD, i + 1):
                    yield Di
    ##
    # end of tree-search routine

    int_to_vertex = list(G)
    vertex_to_int = {u: i for i, u in enumerate(int_to_vertex)}

    if work_on_copy:
        G.relabel(perm=vertex_to_int)
    else:
        G = G.relabel(perm=vertex_to_int, inplace=False)

    if to_dominate is None:
        vertices_to_dominate = set(G)
    else:
        vertices_to_dominate = set(to_dominate)

    if not vertices_to_dominate:
        # base case: vertices_to_dominate is empty
        # the empty set/list is the only minimal DS of the empty set
        yield set()
        return

    peeling = _peel(G, vertices_to_dominate)

    for dom in tree_search(G, peeling, set(), 0):
        yield {int_to_vertex[v] for v in dom}

# ==============================================================================
# Greedy heuristic for dominating set
# ==============================================================================

def greedy_dominating_set(G, k=1, vertices=None, ordering=None, return_sets=False, closest=False):
    r"""
    Return a greedy distance-`k` dominating set of the graph.

    A distance-`k` dominating set `S` of a graph `G` is a set of its vertices of
    minimal cardinality such that any vertex of `G` is in `S` or is at distance
    at most `k` from a vertex in `S`. See the :wikipedia:`Dominating_set`.

    When `G` is directed, vertex `u` can be a dominator of vertex `v` if there
    is a directed path of length at most `k` from `u` to `v`.

    This method implements a greedy heuristic to find a minimal dominatic set.

    INPUT:

    - ``G`` -- a Graph

    - ``k`` -- integer (default: ``1``); the domination distance to consider

    - ``vertices`` -- iterable container of vertices (default: ``None``); when
      specified, return a dominating set of the specified vertices only

    - ``ordering`` -- string (default: ``None``); specify the order in which to
      consider the vertices

      - ``None`` -- if ``vertices`` is ``None``, then consider the vertices in
        the order given by ``list(G)``. Otherwise, consider the vertices in the
        order of iteration of ``vertices``.

      - ``"degree_min"`` -- consider the vertices by increasing degree

      - ``"degree_max"`` -- consider the vertices by decreasing degree

    - ``return_sets`` -- boolean (default: ``False``); whether to return the
      vertices of the dominating set only (default), or a dictionary mapping
      each vertex of the dominating set to the set of vertices it dominates.

    - ``closest`` -- boolean (default: ``False``); whether to attach a vertex to
      its closest dominator or not. This parameter is use only when
      ``return_sets`` is ``True``.

    EXAMPLES:

    Dominating sets of a path::

        sage: from sage.graphs.domination import greedy_dominating_set
        sage: G = graphs.PathGraph(5)
        sage: sorted(greedy_dominating_set(G, ordering=None))
        [0, 2, 4]
        sage: sorted(greedy_dominating_set(G, ordering="degree_min"))
        [0, 2, 4]
        sage: sorted(greedy_dominating_set(G, ordering="degree_max"))
        [1, 3]
        sage: sorted(greedy_dominating_set(G, k=2, ordering=None))
        [0, 3]
        sage: sorted(greedy_dominating_set(G, k=2, ordering="degree_min"))
        [0, 4]
        sage: sorted(greedy_dominating_set(G, k=2, ordering="degree_max"))
        [1, 4]
        sage: greedy_dominating_set(G, k=3, ordering="degree_min", return_sets=True, closest=False)
        {0: {0, 1, 2, 3}, 4: {4}}
        sage: greedy_dominating_set(G, k=3, ordering="degree_min", return_sets=True, closest=True)
        {0: {0, 2, 3}, 4: {1, 4}}

    Asking for a dominating set of a subset of vertices::

        sage: from sage.graphs.domination import greedy_dominating_set
        sage: from sage.graphs.domination import is_dominating
        sage: G = graphs.PetersenGraph()
        sage: vertices = {0, 1, 2, 3, 4, 5}
        sage: dom = greedy_dominating_set(G, vertices=vertices, return_sets=True)
        sage: sorted(dom)
        [0, 2]
        sage: is_dominating(G, dom, focus=vertices)
        True
        sage: is_dominating(G, dom)
        False
        sage: dominated = [u for v in dom for u in dom[v]]
        sage: sorted(dominated) == sorted(vertices)
        True

    Influence of the ordering of the vertices on the result::

        sage: from sage.graphs.domination import greedy_dominating_set
        sage: G = graphs.StarGraph(4)
        sage: greedy_dominating_set(G, vertices=[0, 1, 2, 3, 4])
        [0]
        sage: sorted(greedy_dominating_set(G, vertices=[1, 2, 3, 4, 0]))
        [1, 2, 3, 4]

    Dominating set of a directed graph::

        sage: from sage.graphs.domination import greedy_dominating_set
        sage: D = digraphs.Path(3)
        sage: sorted(greedy_dominating_set(D, vertices=[0, 1, 2]))
        [0, 2]

    TESTS:

    Random tests::

        sage: from sage.graphs.domination import greedy_dominating_set
        sage: from sage.graphs.domination import is_dominating
        sage: G = graphs.RandomGNP(15, .2)
        sage: for o in [None, "degree_min", "degree_max"]:
        ....:     for c in [True, False]:
        ....:         dom = greedy_dominating_set(G, ordering=o, closest=c)
        ....:         if not is_dominating(G, dom):
        ....:             print("something goes wrong")

    Corner cases::

        sage: greedy_dominating_set(Graph())
        []
        sage: greedy_dominating_set(Graph(1))
        [0]
        sage: greedy_dominating_set(Graph(2))
        [0, 1]
        sage: G = graphs.PathGraph(5)
        sage: dom = greedy_dominating_set(G, vertices=[0, 1, 3, 4])

    Check parameters::

        sage: greedy_dominating_set(G, ordering="foo")
        Traceback (most recent call last):
        ...
        ValueError: ordering must be None, "degree_min" or "degree_max"
    """
    if vertices is None:
        vertices = list(G)
    else:
        vertices = [u for u in vertices if u in G]

    if ordering in ["degree_min", "degree_max"]:
        vertices = sorted(vertices, key=G.degree, reverse=ordering.endswith("max"))
    elif ordering is not None:
        raise ValueError('ordering must be None, "degree_min" or "degree_max"')

    if not G:
        return dict() if return_sets else []
    if not k:
        return vertices

    n = G.order()
    dom = dict()
    seen = set()
    # We want to dominate only the set vertices
    to_avoid = set() if len(vertices) == n else set(G).difference(vertices)

    if closest:
        # Attach each dominated vertex to its closest dominator
        from sage.rings.infinity import Infinity
        dominator = {u: (u, +Infinity) for u in vertices}
        for u in vertices:
            if u in seen:
                continue
            dom[u] = set()
            for v, d in G.breadth_first_search(u, distance=k, report_distance=True):
                if v in to_avoid:
                    continue
                if v not in seen:
                    dom[u].add(v)
                    seen.add(v)
                    dominator[v] = (u, d)
                else:
                    x, dx = dominator[v]
                    if dx < d:
                        dom[x].discard(v)
                        dom[u].add(v)
                        dominator[v] = (u, d)

    else:
        for u in vertices:
            if u in seen:
                continue
            dom[u] = set()
            for v in G.breadth_first_search(u, distance=k):
                if v not in to_avoid and v not in seen:
                    dom[u].add(v)
                    seen.add(v)

    if return_sets:
        return dom
    else:
        return list(dom)
