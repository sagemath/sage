r"""
`*`morphisms

This module gathers methods related to homeomorphims, homomorphisms,
isomorphisms, etc. in (di)graphs.

{INDEX_OF_METHODS}

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

from sage.misc.rest_index_of_methods import doc_index, gen_thematic_rest_table_index


@doc_index("Homeomorphism")
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

        sage: G = graphs.WindmillGraph(3, 10)
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


@doc_index("Homeomorphism")
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


_additional_categories = {}
__doc__ = __doc__.replace("{INDEX_OF_METHODS}", gen_thematic_rest_table_index(Graph, _additional_categories))
