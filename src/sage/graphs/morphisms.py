r"""
Morphisms

This module gathers methods related to homeomorphims, homomorphisms,
isomorphisms, etc. in (di)graphs.

**This module contains the following methods**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~is_homeomorphic` | Check whether ``G`` and ``H`` are homeomorphic.
    :meth:`~reduced_homeomorphic_graph` | Return the smallest graph homeomorphic to ``G``.
    :meth:`~has_homomorphism_to` | Check whether there is a homomorphism between two graphs.
    :meth:`~is_homomorphism` | Check whether a mapping is a valid graph homomorphism.
    :meth:`~homomorphisms` | Return an iterator that lazily yields homomorphisms from ``G`` to ``H``.
    :meth:`~homomorphism_count` | Return the total number of homomorphisms from ``G`` to ``H``.

.. TODO::

    - Move methods related to graph automorphisms to this module
    - Move methods related to graph isomorphisms to this module

Methods
-------
"""
# ****************************************************************************
#       Copyright (C) 2025 David Coudert <david.coudert@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


def reduced_homeomorphic_graph(G, allow_multiple_edges=False, allow_loops=False,
                               return_steps=False, immutable=None):
    r"""
    Return the smallest graph homeomorphic to `G`.

    Two graphs `G` and `H` are homeomorphic if there is an isomorphism from some
    subdivision of `G` to some subdivision of `H`. For more details, see the
    :wikipedia:`Homeomorphism_(graph_theory)`.

    By default (i.e., when ``allow_multiple_edges == False`` and ``allow_loops
    == False``), given a graph `G`, a vertex `u` of degree two and its neighbors
    `x` and `y`, with `x \neq y`, this methods replaces the path `(x, u, y)`
    with the edge `(x, y)` unless the graph already has edge `(x, y)`. This
    process is repeated for each vertex of degree two. The resulting graph `H`
    is the smallest graph that is homeomorphic to `G`.

    When ``allow_multiple_edges == True`` and ``allow_loops == False``, this
    method always replaces the path `(x, u, y)` with a new edge `(x, y)`. Hence,
    the resulting graph may have several edges between `x` and `y`. This
    operation is performed only if `x \neq y`.

    When ``allow_loops == True``, this method also assumes that
    ``allow_multiple_edges == True``. If a vertex `u` of degree two is connected
    by two edges to a vertex `x`, this method replaces the two edges by a loop
    edge on `x`.

    For digraphs, the method considers the vertices with in and out degree one.

    INPUT:

    - ``G`` -- a graph or a digraph

    - ``allow_multiple_edges`` -- boolean (default: ``False``); whether to allow
      the creation of new multiple edges.  This parameter is considered ``True``
      when ``allow_loops`` is ``True``.

    - ``allow_loops`` -- boolean (default: ``False``); whether to allow the
      creation of new loops

    - ``return_steps`` -- boolean (default: ``False``); whether to return the
      steps of the reduction as a list of triples `(x, u, y)` indicating that
      path `(x, u, y)` has been replaced by edge `(x, y)`. The original graph
      can be reconstructed by using this list in reverse order.

    - ``immutable`` -- boolean (default: ``None``); whether to create a
      mutable/immutable (di)graph. ``immutable=None`` (default) means that the
      (di)graph and its reduced (di)graph will behave the same way.

    OUTPUT: When ``return_steps`` is ``False``, this method returns the reduced
    graph. When ``return_steps`` is ``True``, this method returns both the
    reduced graph and the ordered list of reduction operations. Each reduction
    operation is a triple `(x, u, y)` indicating that the path `(x, u, y)`, with
    `u` of degree two (or with in and out degree one for digraphs), has been
    replaced by edge `(x, y)`.

    EXAMPLES:

    Reduction of a Cycle Graph::

        sage: G = graphs.CycleGraph(4)
        sage: G.reduced_homeomorphic_graph()
        Graph on 3 vertices
        sage: G.reduced_homeomorphic_graph(allow_multiple_edges=True)
        Multi-graph on 2 vertices
        sage: G.reduced_homeomorphic_graph(allow_loops=True)
        Looped multi-graph on 1 vertex

    Reduction of a Circuit::

        sage: G = digraphs.Circuit(4)
        sage: G.reduced_homeomorphic_graph()
        Digraph on 2 vertices
        sage: G.reduced_homeomorphic_graph(allow_multiple_edges=True)
        Multi-digraph on 2 vertices
        sage: G.reduced_homeomorphic_graph(allow_loops=True)
        Looped multi-digraph on 1 vertex

    Check that the construction is reversible::

        sage: def revert_steps(g, steps):
        ....:     h = g.copy(immutable=False)
        ....:     for P in reversed(steps):
        ....:         h.add_path(P)
        ....:         h.delete_edge(P[0], P[2])
        ....:     return h

        sage: G = graphs.WindmillGraph(3, 5)
        sage: G.order(), G.size()
        (11, 15)
        sage: H, steps = G.reduced_homeomorphic_graph(return_steps=True)
        sage: H.order(), H.size()
        (11, 15)
        sage: G.is_isomorphic(revert_steps(H, steps))
        True
        sage: H, steps = G.reduced_homeomorphic_graph(allow_multiple_edges=True, return_steps=True)
        sage: H.order(), H.size()
        (6, 10)
        sage: G.is_isomorphic(revert_steps(H, steps))
        True
        sage: H, steps = G.reduced_homeomorphic_graph(allow_loops=True, return_steps=True)
        sage: H.order(), H.size()
        (1, 5)
        sage: len(H.loop_edges())
        5
        sage: G.is_isomorphic(revert_steps(H, steps))
        True

    Random digraph::

        sage: G = digraphs.RandomDirectedGNP(20, 0.05, loops=True)
        sage: H, steps = G.reduced_homeomorphic_graph(return_steps=True)
        sage: G.is_isomorphic(revert_steps(H, steps))
        True
        sage: H, steps = G.reduced_homeomorphic_graph(allow_multiple_edges=True, return_steps=True)
        sage: G.is_isomorphic(revert_steps(H, steps))
        True
        sage: H, steps = G.reduced_homeomorphic_graph(allow_loops=True, return_steps=True)
        sage: G.is_isomorphic(revert_steps(H, steps))
        True

    TESTS:

    Check the behavior of parameter ``immutable``::

        sage: G = graphs.CycleGraph(3)
        sage: G.reduced_homeomorphic_graph().is_immutable()
        False
        sage: G.reduced_homeomorphic_graph(immutable=True).is_immutable()
        True
        sage: G = G.copy(immutable=True)
        sage: G.reduced_homeomorphic_graph().is_immutable()
        True
        sage: G.reduced_homeomorphic_graph(immutable=False).is_immutable()
        False
    """
    if allow_loops:
        allow_multiple_edges = True

    if G.is_directed():
        from sage.graphs.digraph import DiGraph as MyGraph

        # candidates is the list of vertices with in and out degree 1
        out_degree_one = (u for u, d in G.out_degree_iterator(labels=True) if d == 1)
        candidates = (u for u, d in G.in_degree_iterator(vertices=out_degree_one, labels=True) if d == 1)

        def get_neighbors(g, u):
            return (next(g.neighbor_in_iterator(u)),
                    next(g.neighbor_out_iterator(u)))

    else:
        from sage.graphs.graph import Graph as MyGraph

        # candidates is the list of vertices with degree 2
        candidates = (u for u, d in G.degree_iterator(labels=True) if d == 2)

        if allow_multiple_edges:

            def get_neighbors(g, u):
                N = g.neighbors(u)
                if len(N) == 1:
                    return N * 2
                return N

        else:

            def get_neighbors(g, u):
                return g.neighbors(u)

    # Copy of the (di)graph with required settings for loops and multiple edges
    H = MyGraph([G, G.edge_iterator(labels=False)], format='vertices_and_edges',
                multiedges=G.allows_multiple_edges() or allow_multiple_edges,
                loops=G.allows_loops() or allow_loops, immutable=False)

    steps = []
    for u in candidates:
        x, y = get_neighbors(H, u)
        if (not allow_loops and x == y) or x == u:
            # The case x = u = y may occur when contracting a cycle
            continue
        if not allow_multiple_edges and H.has_edge(x, y):
            continue
        # Replace path (x, u, y) with edge (x, y)
        H.delete_vertex(u)
        H.add_edge(x, y)
        steps.append((x, u, y))

    if immutable is None:
        immutable = G.is_immutable()
    if immutable:
        H = H.copy(immutable=True)

    if return_steps:
        return H, steps
    return H


def is_homeomorphic(G, H):
    r"""
    Check whether ``G`` and ``H`` are homeomorphic.

    Two graphs `G` and `H` are homeomorphic if there is an isomorphism from some
    subdivision of `G` to some subdivision of `H`. To check whether `G` and `H`
    are homeomorphic, it suffices to check whether their reduced homeomorphic
    (di)graphs are isomorphic. For more details, see the
    :wikipedia:`Homeomorphism_(graph_theory)`.

    INPUT:

    - ``G``, ``H`` -- two (di)graphs

    EXAMPLES::

        sage: G = graphs.RandomGNP(10, .2)
        sage: H = G.copy()
        sage: for e in list(G.edges()):
        ....:     G.subdivide_edge(e, randint(0, 5))
        ....:     H.subdivide_edge(e, randint(0, 5))
        sage: G.is_homeomorphic(H)
        True
        sage: G = graphs.RandomGNP(10, .2)
        sage: G.allow_multiple_edges(True)
        sage: G.add_edges(G.edges())
        sage: H = G.copy()
        sage: for e in list(G.edges()):
        ....:     G.subdivide_edge(e, randint(0, 5))
        ....:     H.subdivide_edge(e, randint(0, 5))
        sage: G.is_homeomorphic(H)
        True

        sage: G = digraphs.RandomDirectedGNP(10, .2)
        sage: H = G.copy()
        sage: for e in list(G.edges()):
        ....:     G.subdivide_edge(e, randint(0, 5))
        ....:     H.subdivide_edge(e, randint(0, 5))
        sage: G.is_homeomorphic(H)
        True
        sage: G = digraphs.RandomDirectedGNP(10, .2)
        sage: G.allow_multiple_edges(True)
        sage: G.add_edges(G.edges())
        sage: H = G.copy()
        sage: for e in list(G.edges()):
        ....:     G.subdivide_edge(e, randint(0, 5))
        ....:     H.subdivide_edge(e, randint(0, 5))
        sage: G.is_homeomorphic(H)
        True

        sage: G = digraphs.RandomDirectedGNP(10, .2)
        sage: G.allow_loops(True)
        sage: G.add_edges((u, u) for u in G if randint(0, 1))
        sage: G.allow_multiple_edges(True)
        sage: G.add_edges(G.edges())
        sage: H = G.copy()
        sage: for e in list(G.edges()):
        ....:     G.subdivide_edge(e, randint(0, 5))
        ....:     H.subdivide_edge(e, randint(0, 5))
        sage: G.is_homeomorphic(H)
        True

    TESTS::

        sage: Graph(1).is_homeomorphic(DiGraph(1))
        False
    """
    if G.is_directed() is not H.is_directed():
        return False
    X = G.reduced_homeomorphic_graph(allow_loops=True, immutable=False)
    Y = H.reduced_homeomorphic_graph(allow_loops=True, immutable=False)
    return X.is_isomorphic(Y)


def has_homomorphism_to(G, H, core=False, solver=None, verbose=0,
                        *, integrality_tolerance=1e-3):
    r"""
    Check whether there is a homomorphism between two graphs.

    A homomorphism from a graph `G` to a graph `H` is a function
    `\phi:V(G)\mapsto V(H)` such that for any edge `uv \in E(G)` the pair
    `\phi(u)\phi(v)` is an edge of `H`.

    Saying that a graph can be `k`-colored is equivalent to saying that it has a
    homomorphism to `K_k`, the complete graph of order `k`.

    For more information, see the :wikipedia:`Graph_homomorphism`.

    INPUT:

    - ``G`` -- the graph to map

    - ``H`` -- the graph to which ``G`` should be sent

    - ``core`` -- boolean (default: ``False``); whether to minimize the size of
      the mapping's image (see examples below). This is set to ``False`` by
      default.

    - ``solver`` -- string (default: ``None``); specifies a Mixed Integer Linear
      Programming (MILP) solver to be used. If set to ``None``, the default one
      is used. For more information on MILP solvers and which default solver is
      used, see the method :meth:`solve
      <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
      :class:`MixedIntegerLinearProgram
      <sage.numerical.mip.MixedIntegerLinearProgram>`.

    - ``verbose`` -- integer (default: 0); sets the level of verbosity. Set to 0
      by default, which means quiet.

    - ``integrality_tolerance`` -- float; parameter for use with MILP solvers
      over an inexact base ring; see
      :meth:`MixedIntegerLinearProgram.get_values`.

    OUTPUT:

    This method returns ``False`` when the homomorphism does not exist, and
    returns the homomorphism otherwise as a dictionary associating a vertex of
    `H` to a vertex of `G`.

    EXAMPLES:

    Is Petersen's graph 3-colorable::

        sage: P = graphs.PetersenGraph()
        sage: P.has_homomorphism_to(graphs.CompleteGraph(3)) is not False               # needs sage.numerical.mip
        True

    An odd cycle admits a homomorphism to a smaller odd cycle, but not to an
    even cycle::

        sage: g = graphs.CycleGraph(9)
        sage: g.has_homomorphism_to(graphs.CycleGraph(5)) is not False                  # needs sage.numerical.mip
        True
        sage: g.has_homomorphism_to(graphs.CycleGraph(7)) is not False                  # needs sage.numerical.mip
        True
        sage: g.has_homomorphism_to(graphs.CycleGraph(4)) is not False                  # needs sage.numerical.mip
        False

    One can compute the core of a graph (with respect to homomorphism)
    with this method::

        sage: g = graphs.CycleGraph(8)
        sage: mapping = g.has_homomorphism_to(g, core=True)
        sage: print(f"The size of the core is {len(set(mapping.values()))}")
        The size of the core is 2
        sage: g = graphs.CycleGraph(9)
        sage: mapping = g.has_homomorphism_to(g, core=True)
        sage: print(f"The size of the core is {len(set(mapping.values()))}")
        The size of the core is 9

    The chromatic number of a graph is the order of the smallest clique to which
    it has an homomorphism::

        sage: g = graphs.CycleGraph(9)
        sage: g.chromatic_number()
        3
        sage: g.has_homomorphism_to(graphs.CompleteGraph(3)) is not False
        True
        sage: g.has_homomorphism_to(graphs.CompleteGraph(2)) is not False
        False
        sage: K6 = graphs.CompleteGraph(6)
        sage: g.has_homomorphism_to(K6) is not False
        True
        sage: mapping = g.has_homomorphism_to(K6, core=True)
        sage: print(f"The size of the core is {len(set(mapping.values()))}")
        The size of the core is 3

    A circuit of order `n` admits a homomorphism to smaller circuit of order `p
    \leq n` if `p` is a divisor of `n`::

        sage: g = digraphs.Circuit(12)
        sage: [i for i in range(2, g.order() + 1)                                       # needs sage.numerical.mip
        ....:  if g.has_homomorphism_to(digraphs.Circuit(i)) is not False]
        [2, 3, 4, 6, 12]

    TESTS::

        sage: Graph(1).has_homomorphism_to(DiGraph(1))
        False
    """
    G._scream_if_not_simple()
    H._scream_if_not_simple()
    if G.is_directed() is not H.is_directed():
        return False
    undirected = not G.is_directed()

    from sage.numerical.mip import MIPSolverException, MixedIntegerLinearProgram
    p = MixedIntegerLinearProgram(solver=solver, maximization=False)
    b = p.new_variable(binary=True)

    # Each vertex has an image
    for ug in G:
        p.add_constraint(p.sum(b[ug, uh] for uh in H) == 1)

    nonedges = H.complement().edges(sort=False, labels=False)
    for ug, vg in G.edges(sort=False, labels=False):
        # Two adjacent vertices cannot be mapped to the same element
        for uh in H:
            p.add_constraint(b[ug, uh] + b[vg, uh] <= 1)

        # Two adjacent vertices cannot be mapped to no adjacent vertices
        for uh, vh in nonedges:
            p.add_constraint(b[ug, uh] + b[vg, vh] <= 1)

        if undirected:
            # Both directions of edges must be considered for undirected graphs
            for uh, vh in nonedges:
                p.add_constraint(b[ug, vh] + b[vg, uh] <= 1)

    # Minimize the mapping's size
    if core:

        # The value of m is one if the corresponding vertex of H is used
        m = p.new_variable(nonnegative=True)
        for uh in H:
            for ug in G:
                p.add_constraint(b[ug, uh] <= m[uh])

        # Minimize the number of used vertices of H
        p.set_objective(p.sum(m[vh] for vh in H))

    try:
        p.solve(log=verbose)
    except MIPSolverException:
        return False

    b = p.get_values(b, convert=bool, tolerance=integrality_tolerance)
    return dict(x[0] for x in b.items() if x[1])


def is_homomorphism(G, H, mapping):
    r"""
    Check whether ``mapping`` is a valid graph homomorphism from ``G`` to ``H``.

    Let `G` and `H` be simple undirected graphs. A homomorphism from `G` to `H`
    is a vertex mapping `\phi: V(G) \to V(H)` such that for every edge
    `\{u, v\} \in E(G)`, the image `\{\phi(u), \phi(v)\}` is an edge in `H`.

    Every vertex of `G` must have an assigned image in ``mapping``, and each
    such image must be a vertex of `H`. Extra keys in ``mapping`` that do not
    belong to `V(G)` are ignored. Mappings are not required to be injective:
    non-injective mappings are valid homomorphisms as long as edge preservation
    holds.

    INPUT:

    - ``G`` -- a simple undirected graph (domain)

    - ``H`` -- a simple undirected graph (codomain)

    - ``mapping`` -- a mapping (e.g. :class:`dict`) from vertices of ``G`` to
      vertices of ``H``

    OUTPUT: boolean

    EXAMPLES:

    A valid 3-coloring of a 5-cycle `C_5` is a homomorphism to `K_3`::

        sage: from sage.graphs.morphisms import is_homomorphism
        sage: C5 = graphs.CycleGraph(5)
        sage: K3 = graphs.CompleteGraph(3)
        sage: phi = {0: 0, 1: 1, 2: 0, 3: 1, 4: 2}
        sage: is_homomorphism(C5, K3, phi)
        True

    Adjacent vertices in `G` cannot be mapped to non-adjacent vertices in `H`::

        sage: P3 = graphs.PathGraph(3)
        sage: H = Graph([(0, 1), (1, 2)])
        sage: is_homomorphism(P3, H, {0: 0, 1: 1, 2: 0})
        True
        sage: is_homomorphism(P3, H, {0: 0, 1: 2, 2: 1})
        False

    Adjacent vertices cannot be mapped to the same vertex in a simple graph::

        sage: is_homomorphism(C5, K3, {0: 0, 1: 0, 2: 1, 3: 2, 4: 1})
        False

    Non-injective homomorphisms are valid::

        sage: C4 = graphs.CycleGraph(4)
        sage: K2 = graphs.CompleteGraph(2)
        sage: is_homomorphism(C4, K2, {0: 0, 1: 1, 2: 0, 3: 1})
        True

    Extra mapping keys not in `V(G)` are accepted and ignored::

        sage: is_homomorphism(C4, K2, {0: 0, 1: 1, 2: 0, 3: 1, 'extra': 99})
        True

    Missing vertices of `G` or target values not in `V(H)` return ``False``::

        sage: is_homomorphism(C4, K2, {0: 0, 1: 1, 2: 0})
        False
        sage: is_homomorphism(C4, K2, {0: 0, 1: 1, 2: 0, 3: 99})
        False

    Empty source graph has exactly one homomorphism (the empty mapping)::

        sage: is_homomorphism(Graph(), K3, {})
        True
        sage: is_homomorphism(Graph(), Graph(), {})
        True
        sage: is_homomorphism(Graph(), Graph(), {'irrelevant': 42, 'extra': 'ignored'})
        True

    A non-empty source graph to an empty target has no homomorphisms::

        sage: is_homomorphism(K3, Graph(), {0: 0, 1: 1, 2: 2})
        False

    Edgeless source graphs with isolated vertices::

        sage: E3 = Graph(3)
        sage: is_homomorphism(E3, K3, {0: 0, 1: 0, 2: 0})
        True

    Non-mapping inputs raise :class:`TypeError`::

        sage: is_homomorphism(C5, K3, [0, 1, 0, 1, 2])
        Traceback (most recent call last):
        ...
        TypeError: mapping must be a Mapping, got list
        sage: is_homomorphism(C5, K3, 42)
        Traceback (most recent call last):
        ...
        TypeError: mapping must be a Mapping, got Integer
        sage: is_homomorphism(C5, K3, int(42))
        Traceback (most recent call last):
        ...
        TypeError: mapping must be a Mapping, got int

    Directed graphs and non-simple graphs are rejected::

        sage: is_homomorphism(DiGraph(1), Graph(1), {0: 0})
        Traceback (most recent call last):
        ...
        ValueError: directed graphs are not supported
        sage: is_homomorphism(Graph(1), DiGraph(1), {0: 0})
        Traceback (most recent call last):
        ...
        ValueError: directed graphs are not supported
        sage: is_homomorphism(Graph([(0, 0)], loops=True), Graph(1), {0: 0})
        Traceback (most recent call last):
        ...
        ValueError: This method is not known to work on graphs with loops...
        sage: is_homomorphism(Graph([(0, 1), (0, 1)], multiedges=True), Graph(2), {0: 0, 1: 1})
        Traceback (most recent call last):
        ...
        ValueError: This method is not known to work on graphs with multiedges...
    """
    from collections.abc import Mapping

    if not isinstance(mapping, Mapping):
        raise TypeError(f"mapping must be a Mapping, got {type(mapping).__name__}")

    if G.is_directed() or H.is_directed():
        raise ValueError("directed graphs are not supported")

    G._scream_if_not_simple()
    H._scream_if_not_simple()

    # Every vertex of G must be in mapping and map to a vertex of H
    for u in G:
        if u not in mapping:
            return False
        if mapping[u] not in H:
            return False

    # Every edge of G must map to an edge in H
    for u, v in G.edges(sort=False, labels=False):
        if not H.has_edge(mapping[u], mapping[v]):
            return False

    return True


def _vertex_ordering(C):
    r"""
    Return a vertex ordering of connected graph ``C`` for backtracking search.

    The ordering ensures that each vertex after the first has at least one
    neighbor among earlier vertices in the ordering.

    TESTS:

    Verify ordering property on a connected cycle::

        sage: from sage.graphs.morphisms import _vertex_ordering
        sage: C = graphs.CycleGraph(5)
        sage: order = _vertex_ordering(C)
        sage: len(order) == 5 and set(order) == set(C)
        True
        sage: all(any(C.has_edge(order[i], order[j]) for j in range(i))
        ....:     for i in range(1, len(order)))
        True
        sage: _vertex_ordering(Graph())
        Traceback (most recent call last):
        ...
        ValueError: graph must be non-empty and connected
        sage: _vertex_ordering(graphs.PathGraph(2) + graphs.PathGraph(2))
        Traceback (most recent call last):
        ...
        ValueError: graph must be non-empty and connected
    """
    if not C or not C.is_connected():
        raise ValueError("graph must be non-empty and connected")

    n = C.order()
    if n == 1:
        return list(C)

    u0 = max(C, key=lambda v: (C.degree(v), str(v)))
    order = [u0]
    visited = {u0}
    frontier = set(C.neighbors(u0))

    while len(order) < n:
        v = max(frontier, key=lambda u: (sum(1 for w in C.neighbors(u) if w in visited), C.degree(u), str(u)))
        order.append(v)
        visited.add(v)
        frontier.remove(v)
        for w in C.neighbors(v):
            if w not in visited:
                frontier.add(w)

    return order


def _homomorphisms_connected(C, H):
    r"""
    Yield homomorphisms from connected graph C to H.
    """
    if not C:
        yield {}
        return
    if not H:
        return
    if C.order() == 1:
        u = next(iter(C))
        for v in H:
            yield {u: v}
        return

    order = _vertex_ordering(C)
    n = len(order)
    pos = {u: i for i, u in enumerate(order)}
    pred = [[p for p in C.neighbors(order[i]) if pos[p] < i] for i in range(n)]

    adj_H = {v: set(H.neighbors(v)) for v in H}

    u0_deg = C.degree(order[0])
    if u0_deg > 0:
        init_cand = [v for v in H if H.degree(v) > 0]
    else:
        init_cand = list(H)

    current_phi = [None] * n

    def backtrack(idx):
        if idx == n:
            yield {order[i]: current_phi[i] for i in range(n)}
            return

        if idx == 0:
            candidates = init_cand
        else:
            p_list = pred[idx]
            candidates = adj_H[current_phi[pos[p_list[0]]]]
            if len(p_list) > 1:
                candidates = candidates.intersection(*(adj_H[current_phi[pos[p]]] for p in p_list[1:]))

        for v in candidates:
            current_phi[idx] = v
            yield from backtrack(idx + 1)

    yield from backtrack(0)


def _homomorphisms_generator(G, H):
    r"""
    Internal generator yielding homomorphisms from G to H.
    """
    if not G:
        yield {}
        return

    if not H:
        return

    components = G.connected_components(sort=False)
    if len(components) == 1:
        yield from _homomorphisms_connected(G, H)
        return

    subgraphs = [G.subgraph(c) for c in components]

    def combine(comp_idx, current_mapping):
        if comp_idx == len(subgraphs):
            yield dict(current_mapping)
            return

        for phi in _homomorphisms_connected(subgraphs[comp_idx], H):
            current_mapping.update(phi)
            yield from combine(comp_idx + 1, current_mapping)
            for k in phi:
                del current_mapping[k]

    current_mapping = {}
    yield from combine(0, current_mapping)


def homomorphisms(G, H, max_count=None):
    r"""
    Return an iterator that lazily yields homomorphisms from ``G`` to ``H``.

    Let `G` and `H` be simple undirected graphs. A homomorphism from `G` to `H`
    is a vertex mapping `\phi: V(G) \to V(H)` such that for every edge
    `\{u, v\} \in E(G)`, `\{\phi(u), \phi(v)\} \in E(H)`.

    This function returns an iterator that searches for homomorphisms
    using depth-first backtracking with degree-ordered connected vertex
    traversal and connected component factorization. Mappings are yielded
    lazily as they are found. The enumeration order is unspecified and
    depends on internal vertex orderings.

    ALGORITHMIC COMPLEXITY:

    For a fixed simple undirected target graph `H`, deciding whether there exists
    a homomorphism from `G` to `H` is solvable in polynomial time if `H` is
    bipartite, and is NP-complete otherwise. Counting the number of homomorphisms
    is #P-complete even for `H = K_3`.

    In the worst case, the total number of candidate vertex mappings is
    `|V(H)|^{|V(G)|}`. Backtracking can explore exponentially many partial
    mappings before finding a solution or certifying non-existence. Actual
    running time depends heavily on the structure of `G` and `H`, pruning
    efficiency from degree-ordered neighbor constraints, and the total size of
    the output. While counting can exploit component factorization without
    materializing product combinations, enumeration runtime is necessarily
    lower-bounded by the number of homomorphisms yielded.

    INPUT:

    - ``G`` -- a simple undirected graph (domain)

    - ``H`` -- a simple undirected graph (codomain)

    - ``max_count`` -- integer or ``None`` (default: ``None``); maximum number
      of homomorphisms to yield. If ``None``, all homomorphisms are yielded.

    OUTPUT: an iterator yielding dictionaries representing valid homomorphisms

    EXAMPLES:

    Homomorphisms from `C_5` to `K_3` correspond to proper 3-colorings::

        sage: from sage.graphs.morphisms import homomorphisms, is_homomorphism
        sage: C5 = graphs.CycleGraph(5)
        sage: K3 = graphs.CompleteGraph(3)
        sage: h_list = list(homomorphisms(C5, K3))
        sage: len(h_list)
        30

    There are no homomorphisms from an odd cycle to a bipartite graph::

        sage: K2 = graphs.CompleteGraph(2)
        sage: list(homomorphisms(K3, K2))
        []

    Homomorphisms from `K_2` to `K_3` (directed edges / ordered pairs)::

        sage: len(list(homomorphisms(K2, K3)))
        6

    Homomorphisms from bipartite `C_4` to `K_2` (2 proper 2-colorings)::

        sage: len(list(homomorphisms(graphs.CycleGraph(4), K2)))
        2

    Using ``max_count`` to limit the number of yielded homomorphisms::

        sage: len(list(homomorphisms(C5, K3, max_count=1)))
        1
        sage: res1 = list(homomorphisms(C5, K3, max_count=1))[0]
        sage: is_homomorphism(C5, K3, res1)
        True
        sage: len(list(homomorphisms(C5, K3, max_count=5)))
        5
        sage: list(homomorphisms(C5, K3, max_count=0))
        []
        sage: len(list(homomorphisms(C5, K3, max_count=100)))
        30
        sage: len(list(homomorphisms(C5, K3, max_count=Integer(5))))
        5

    Validation of ``max_count`` parameter::

        sage: list(homomorphisms(C5, K3, max_count=-1))
        Traceback (most recent call last):
        ...
        ValueError: max_count must be a non-negative integer or None
        sage: list(homomorphisms(C5, K3, max_count=True))
        Traceback (most recent call last):
        ...
        TypeError: max_count must be a non-negative integer or None
        sage: list(homomorphisms(C5, K3, max_count=2.5))
        Traceback (most recent call last):
        ...
        TypeError: max_count must be a non-negative integer or None
        sage: list(homomorphisms(C5, K3, max_count="5"))
        Traceback (most recent call last):
        ...
        TypeError: max_count must be a non-negative integer or None

    Degenerate boundary cases::

        sage: list(homomorphisms(Graph(), K3))
        [{}]
        sage: list(homomorphisms(Graph(), Graph()))
        [{}]
        sage: list(homomorphisms(Graph(), K3, max_count=0))
        []
        sage: list(homomorphisms(K3, Graph()))
        []
        sage: list(homomorphisms(Graph(1), Graph()))
        []

    Edgeless source graphs with isolated vertices::

        sage: E3 = Graph(3)
        sage: len(list(homomorphisms(E3, K3)))
        27
        sage: len(list(homomorphisms(E3, K3, max_count=4)))
        4

    Disconnected graphs factor across connected components and obey ``max_count``::

        sage: G = graphs.PathGraph(2) + graphs.PathGraph(2)
        sage: len(list(homomorphisms(G, K3)))
        36
        sage: len(list(homomorphisms(G, K3, max_count=5)))
        5
        sage: all(is_homomorphism(G, K3, phi) for phi in homomorphisms(G, K3, max_count=5))
        True
        sage: list(homomorphisms(G, K3, max_count=0))
        []

    Laziness verification: halting on demand without computing all solutions::

        sage: it = homomorphisms(C5, K3)
        sage: is_homomorphism(C5, K3, next(it))
        True

    Directed and non-simple graphs are rejected::

        sage: list(homomorphisms(DiGraph(1), Graph(1)))
        Traceback (most recent call last):
        ...
        ValueError: directed graphs are not supported
        sage: list(homomorphisms(Graph(1), DiGraph(1)))
        Traceback (most recent call last):
        ...
        ValueError: directed graphs are not supported
        sage: list(homomorphisms(Graph([(0, 0)], loops=True), Graph(1)))
        Traceback (most recent call last):
        ...
        ValueError: This method is not known to work on graphs with loops...
        sage: list(homomorphisms(Graph([(0, 1), (0, 1)], multiedges=True), Graph(2)))
        Traceback (most recent call last):
        ...
        ValueError: This method is not known to work on graphs with multiedges...
    """
    from sage.rings.integer import Integer

    if G.is_directed() or H.is_directed():
        raise ValueError("directed graphs are not supported")

    G._scream_if_not_simple()
    H._scream_if_not_simple()

    if max_count is not None:
        if isinstance(max_count, bool) or not isinstance(max_count, (int, Integer)):
            raise TypeError("max_count must be a non-negative integer or None")
        if max_count < 0:
            raise ValueError("max_count must be a non-negative integer or None")

    gen = _homomorphisms_generator(G, H)
    if max_count is None:
        return gen
    from itertools import islice
    return islice(gen, int(max_count))


def _homomorphism_count_connected(C, H):
    r"""
    Return the number of homomorphisms from connected graph C to H.
    """
    if not C:
        return 1
    if not H:
        return 0
    if C.order() == 1:
        return H.order()

    order = _vertex_ordering(C)
    n = len(order)
    pos = {u: i for i, u in enumerate(order)}
    pred = [[p for p in C.neighbors(order[i]) if pos[p] < i] for i in range(n)]

    adj_H = {v: set(H.neighbors(v)) for v in H}

    u0_deg = C.degree(order[0])
    if u0_deg > 0:
        init_cand = [v for v in H if H.degree(v) > 0]
    else:
        init_cand = list(H)

    current_phi = [None] * n
    total = 0

    def backtrack(idx):
        nonlocal total
        if idx == n:
            total += 1
            return

        if idx == 0:
            candidates = init_cand
        else:
            p_list = pred[idx]
            candidates = adj_H[current_phi[pos[p_list[0]]]]
            if len(p_list) > 1:
                candidates = candidates.intersection(*(adj_H[current_phi[pos[p]]] for p in p_list[1:]))

        for v in candidates:
            current_phi[idx] = v
            backtrack(idx + 1)

    backtrack(0)
    return total


def homomorphism_count(G, H):
    r"""
    Return the total number of homomorphisms from ``G`` to ``H``.

    Let `G` and `H` be simple undirected graphs. A homomorphism from `G` to `H`
    is a vertex mapping `\phi: V(G) \to V(H)` such that for every edge
    `\{u, v\} \in E(G)`, `\{\phi(u), \phi(v)\} \in E(H)`.

    Homomorphism counting is multiplicative across connected components of
    undirected graphs:

    .. MATH::

        |\mathrm{Hom}(G_1 \sqcup G_2, H)| = |\mathrm{Hom}(G_1, H)| \times |\mathrm{Hom}(G_2, H)|

    For an edgeless graph `E_n` on `n` isolated vertices into `H` with `m`
    vertices, edge constraints are vacuously satisfied, yielding exactly `m^n`
    homomorphisms.

    ALGORITHMIC COMPLEXITY:

    For a fixed simple undirected target graph `H`, deciding whether there exists
    a homomorphism from `G` to `H` is solvable in polynomial time if `H` is
    bipartite, and is NP-complete otherwise. Counting homomorphisms is
    #P-complete even for `H = K_3`.

    In the worst case, the total number of candidate vertex mappings is
    `|V(H)|^{|V(G)|}`, and backtracking can explore an exponential search
    space. Connected component factorization allows counting disconnected
    graphs by multiplying component counts directly without materializing the
    Cartesian product of assignments across components. However, on large or
    dense connected components, worst-case exponential complexity remains
    unavoidable.

    INPUT:

    - ``G`` -- a simple undirected graph (domain)

    - ``H`` -- a simple undirected graph (codomain)

    OUTPUT: an :class:`Integer` representing the exact number of homomorphisms

    EXAMPLES:

    Counting homomorphisms from the Petersen graph to `K_3` (3-colorings)::

        sage: from sage.graphs.morphisms import homomorphism_count
        sage: P = graphs.PetersenGraph()
        sage: K3 = graphs.CompleteGraph(3)
        sage: homomorphism_count(P, K3)
        120

    Multiplicative component factorization on disconnected graphs::

        sage: homomorphism_count(P + P, K3)
        14400
        sage: 120 * 120
        14400

    If any component admits zero homomorphisms, the total count is zero::

        sage: K2 = graphs.CompleteGraph(2)
        sage: homomorphism_count(P + graphs.CompleteGraph(3), K2)
        0

    Edgeless source graphs with isolated vertices follow the `m^n` rule::

        sage: E3 = Graph(3)
        sage: K4 = graphs.CompleteGraph(4)
        sage: homomorphism_count(E3, K4)
        64
        sage: 4^3
        64

    Boundary cases::

        sage: homomorphism_count(Graph(), K3)
        1
        sage: homomorphism_count(K3, Graph())
        0
        sage: homomorphism_count(Graph(), Graph())
        1

    Endomorphisms of complete graphs (`K_q \to K_q` is `q!`)::

        sage: factorial(4)
        24
        sage: homomorphism_count(K4, K4)
        24

    Directed and non-simple graphs are rejected::

        sage: homomorphism_count(DiGraph(1), Graph(1))
        Traceback (most recent call last):
        ...
        ValueError: directed graphs are not supported
        sage: homomorphism_count(Graph(1), DiGraph(1))
        Traceback (most recent call last):
        ...
        ValueError: directed graphs are not supported
        sage: homomorphism_count(Graph([(0, 0)], loops=True), Graph(1))
        Traceback (most recent call last):
        ...
        ValueError: This method is not known to work on graphs with loops...
        sage: homomorphism_count(Graph([(0, 1), (0, 1)], multiedges=True), Graph(2))
        Traceback (most recent call last):
        ...
        ValueError: This method is not known to work on graphs with multiedges...

    TESTS:

    Cross-verification across a finite graph corpus::

        sage: from sage.graphs.morphisms import (has_homomorphism_to,
        ....:                                    homomorphisms, homomorphism_count,
        ....:                                    is_homomorphism)
        sage: corpus = [
        ....:     graphs.CompleteGraph(1),
        ....:     graphs.CompleteGraph(2),
        ....:     graphs.CompleteGraph(3),
        ....:     graphs.CycleGraph(4),
        ....:     graphs.CycleGraph(5),
        ....:     graphs.PathGraph(3),
        ....:     graphs.PathGraph(4),
        ....:     graphs.CompleteBipartiteGraph(2, 2),
        ....: ]
        sage: for g in corpus:
        ....:     for h in corpus[:4]:
        ....:         c = homomorphism_count(g, h)
        ....:         h_list = list(homomorphisms(g, h))
        ....:         assert c == len(h_list)
        ....:         assert all(is_homomorphism(g, h, phi) for phi in h_list)

    Cross-verification with chromatic polynomials::

        sage: P = graphs.PetersenGraph()
        sage: P_poly = P.chromatic_polynomial()
        sage: all(homomorphism_count(P, graphs.CompleteGraph(q)) == P_poly(q)
        ....:     for q in [3, 4, 5])
        True

        sage: C5 = graphs.CycleGraph(5)
        sage: C5_poly = C5.chromatic_polynomial()
        sage: all(homomorphism_count(C5, graphs.CompleteGraph(q)) == C5_poly(q)
        ....:     for q in [2, 3, 4])
        True

        sage: T = graphs.PathGraph(6)
        sage: q = 3
        sage: homomorphism_count(T, graphs.CompleteGraph(q)) == q * (q - 1)^5
        True

    Agreement with :func:`has_homomorphism_to`::

        sage: has_homomorphism_to(graphs.CompleteGraph(3), graphs.CompleteGraph(2)) is False  # needs sage.numerical.mip
        True
        sage: homomorphism_count(graphs.CompleteGraph(3), graphs.CompleteGraph(2)) == 0
        True
        sage: has_homomorphism_to(graphs.CycleGraph(4), graphs.CompleteGraph(2)) is not False  # needs sage.numerical.mip
        True
        sage: homomorphism_count(graphs.CycleGraph(4), graphs.CompleteGraph(2)) > 0
        True

    Arbitrary non-integer vertex labels::

        sage: G_str = Graph([('a', 'b'), ('b', 'c')])
        sage: H_str = Graph([('x', 'y'), ('y', 'z')])
        sage: homomorphism_count(G_str, H_str) == len(list(homomorphisms(G_str, H_str)))
        True
        sage: all(is_homomorphism(G_str, H_str, phi) for phi in homomorphisms(G_str, H_str))
        True
    """
    from sage.rings.integer import Integer

    if G.is_directed() or H.is_directed():
        raise ValueError("directed graphs are not supported")

    G._scream_if_not_simple()
    H._scream_if_not_simple()

    if not G:
        return Integer(1)
    if not H:
        return Integer(0)

    m = H.order()
    components = G.connected_components(sort=False)

    num_isolated = 0
    total_count = 1

    for c in components:
        if len(c) == 1:
            num_isolated += 1
        else:
            subG = G.subgraph(c)
            c_count = _homomorphism_count_connected(subG, H)
            if c_count == 0:
                return Integer(0)
            total_count *= c_count

    if num_isolated > 0:
        total_count *= m ** num_isolated

    return Integer(total_count)
