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
    :meth:`~is_2isomorphic` | Check whether two graphs have isomorphic cycle matroids.
    :meth:`~verify_2isomorphism_certificate` | Verify a Whitney-operation certificate.
    :meth:`~reduced_homeomorphic_graph` | Return the smallest graph homeomorphic to ``G``.
    :meth:`~has_homomorphism_to` | Check whether there is a homomorphism between two graphs.

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


_TWO_ISOMORPHISM_COLORS = (
    'loop edge',
    'coloop edge',
    'real edge',
    'virtual port',
    'Q component',
    'S component',
    'P component',
    'R component',
    'R vertex',
)

_TWO_ISOMORPHISM_EDGE_MAPPING_ONLY = object()


def _two_isomorphism_work_graph(G):
    """
    Return a uniquely edge-labeled loopless copy of ``G``.

    The labels created here identify edge *occurrences*.  In particular, they
    remain distinct when parallel edges have identical endpoints and labels.
    """
    from sage.graphs.graph import Graph

    edges = list(G.edge_iterator())
    vertices = list(G)
    vertex_to_int = {vertex: i for i, vertex in enumerate(vertices)}
    marker = object()
    labels = [(marker, i) for i in range(len(edges))]
    label_to_edge = {label: i for i, label in enumerate(labels)}
    loop_edges = []
    nonloop_edges = []
    active_vertices = set()
    for i, (u, v, _) in enumerate(edges):
        if u == v:
            loop_edges.append(i)
        else:
            u_int = vertex_to_int[u]
            v_int = vertex_to_int[v]
            active_vertices.add(u_int)
            active_vertices.add(v_int)
            nonloop_edges.append((u_int, v_int, labels[i]))

    work = Graph([active_vertices, nonloop_edges],
                 format='vertices_and_edges',
                 loops=False, multiedges=True)
    return edges, work, loop_edges, label_to_edge


def _two_isomorphism_encoding(work, loop_edges, label_to_edge, number_of_edges,
                              decomposition=False):
    """
    Encode the cycle matroid of ``work`` as a vertex-colored simple graph.

    The encoding retains incidences only in rigid SPQR components.  Edges and
    virtual ports in series and parallel components are deliberately unordered.
    """
    from sage.graphs.graph import Graph

    auxiliary = Graph(multiedges=False, loops=False)
    colors = {color: [] for color in _TWO_ISOMORPHISM_COLORS}
    seen_real_edges = set()
    decompositions = []

    def add_vertex(vertex, color):
        auxiliary.add_vertex(vertex)
        colors[color].append(vertex)

    def add_unordered_component(kind, block_index, vertices, edges):
        """
        Encode a block that is itself a series or parallel component.
        """
        anchor = ('component', block_index, 0)
        add_vertex(anchor, '{} component'.format(kind))
        skeleton_edges = []
        real_edges = []
        for u, v, label in edges:
            edge = label_to_edge[label]
            edge_node = ('edge', edge)
            add_vertex(edge_node, 'real edge')
            seen_real_edges.add(edge)
            auxiliary.add_edge(anchor, edge_node)
            if decomposition:
                skeleton_edges.append((u, v, edge_node))
                real_edges.append(edge)
        if decomposition:
            skeleton = Graph([vertices, skeleton_edges],
                             format='vertices_and_edges', loops=False,
                             multiedges=True)
            decompositions.append({
                'block_index': block_index,
                'edge_ids': frozenset(real_edges),
                'components': [{
                    'kind': kind,
                    'anchor': anchor,
                    'skeleton': skeleton,
                }],
                'tree_edges': [],
            })

    for edge in loop_edges:
        colors['loop edge'].append(('edge', edge))
        seen_real_edges.add(edge)

    blocks = work.blocks_and_cut_vertices()[0] if work else []
    vertex_blocks = {}
    for block_index, vertices in enumerate(blocks):
        for vertex in vertices:
            vertex_blocks.setdefault(vertex, []).append(block_index)
    vertex_block_sets = {vertex: set(indices)
                         for vertex, indices in vertex_blocks.items()}
    block_edges = [[] for _ in blocks]
    endpoint_block = {}
    for edge in work.edge_iterator():
        u, v, _ = edge
        endpoints = (u, v) if u < v else (v, u)
        try:
            block_index = endpoint_block[endpoints]
        except KeyError:
            blocks_u = vertex_blocks[u]
            blocks_v = vertex_blocks[v]
            if len(blocks_u) == 1:
                block_index = blocks_u[0]
            elif len(blocks_v) == 1:
                block_index = blocks_v[0]
            else:
                common = vertex_block_sets[u].intersection(
                    vertex_block_sets[v])
                if len(common) != 1:
                    raise RuntimeError(
                        "an edge must belong to exactly one block")
                block_index = common.pop()
            endpoint_block[endpoints] = block_index
        block_edges[block_index].append(edge)

    for block_index, (vertices, edges) in enumerate(zip(blocks, block_edges)):
        if not edges:
            continue

        # A one-edge block is a coloop.  It must have the same representation
        # whether it is an isolated K2 or a bridge attached to another block.
        if len(edges) == 1:
            label = edges[0][2]
            edge = label_to_edge[label]
            colors['coloop edge'].append(('edge', edge))
            seen_real_edges.add(edge)
            if decomposition:
                decompositions.append({
                    'block_index': block_index,
                    'edge_ids': frozenset([edge]),
                    'components': [],
                    'tree_edges': [],
                })
            continue

        # Avoid rebuilding large parallel classes and cycles as mutable
        # multigraphs: their matroid encoding is just one unordered component.
        if len(vertices) == 2:
            add_unordered_component('P', block_index, vertices, edges)
            continue
        if len(vertices) == len(edges):
            degrees = dict.fromkeys(vertices, 0)
            for u, v, _ in edges:
                degrees[u] += 1
                degrees[v] += 1
            if all(degree == 2 for degree in degrees.values()):
                add_unordered_component('S', block_index, vertices, edges)
                continue

        block = Graph([vertices, edges], format='vertices_and_edges',
                      loops=False, multiedges=True)
        tree = block.spqr_tree()
        components = list(tree)
        tree_edges = list(tree.edge_iterator(labels=False))

        component_index = {component: i
                           for i, component in enumerate(components)}
        component_labels = []
        for _, skeleton in components:
            labels = {label for _, _, label in skeleton.edge_iterator()}
            component_labels.append(labels)

        # Recover the two occurrences of every virtual edge from the public
        # SPQR tree.  Real labels are private and globally unique, so the sole
        # label shared by adjacent skeletons is their virtual port.
        virtual_port = {}
        for port, (left, right) in enumerate(tree_edges):
            i = component_index[left]
            j = component_index[right]
            shared = component_labels[i].intersection(component_labels[j])
            if len(shared) != 1:
                raise RuntimeError("invalid SPQR virtual-edge representation")
            label = shared.pop()
            if label in label_to_edge:
                raise RuntimeError("a real edge occurs in two SPQR components")
            virtual_port[i, label] = port
            virtual_port[j, label] = port

        port_nodes = {port: [] for port in range(len(tree_edges))}
        component_records = []
        for i, (kind, skeleton) in enumerate(components):
            if kind not in 'QSPR':
                raise RuntimeError("unknown SPQR component type")
            anchor = ('component', block_index, i)
            add_vertex(anchor, '{} component'.format(kind))

            rigid_vertices = {}
            skeleton_edges = []
            if kind == 'R':
                for j, vertex in enumerate(skeleton.vertex_iterator()):
                    rigid_vertex = ('R vertex', block_index, i, j)
                    rigid_vertices[vertex] = rigid_vertex
                    add_vertex(rigid_vertex, 'R vertex')

            for j, (u, v, label) in enumerate(skeleton.edge_iterator()):
                if label in label_to_edge:
                    edge = label_to_edge[label]
                    edge_node = ('edge', edge)
                    add_vertex(edge_node, 'real edge')
                    seen_real_edges.add(edge)
                else:
                    try:
                        port = virtual_port[i, label]
                    except KeyError:
                        raise RuntimeError("unmatched SPQR virtual edge") from None
                    edge_node = ('port', block_index, i, j)
                    add_vertex(edge_node, 'virtual port')
                    port_nodes[port].append(edge_node)

                auxiliary.add_edge(anchor, edge_node)
                if decomposition:
                    skeleton_edges.append((u, v, edge_node))
                if kind == 'R':
                    auxiliary.add_edge(edge_node, rigid_vertices[u])
                    auxiliary.add_edge(edge_node, rigid_vertices[v])

            if decomposition:
                component_records.append({
                    'kind': kind,
                    'anchor': anchor,
                    'skeleton': Graph(
                        [list(skeleton.vertex_iterator()), skeleton_edges],
                        format='vertices_and_edges', loops=False,
                        multiedges=True),
                })

        decomposition_tree_edges = []
        for nodes in port_nodes.values():
            if len(nodes) != 2:
                raise RuntimeError("an SPQR virtual edge must have two ports")
            auxiliary.add_edge(nodes[0], nodes[1])
            if decomposition:
                decomposition_tree_edges.append(
                    (nodes[0][2], nodes[1][2], nodes[0], nodes[1]))

        if decomposition:
            decompositions.append({
                'block_index': block_index,
                'edge_ids': frozenset(
                    label_to_edge[label] for _, _, label in edges),
                'components': component_records,
                'tree_edges': decomposition_tree_edges,
            })

    if seen_real_edges != set(range(number_of_edges)):
        raise RuntimeError("some graph edges are missing from the block decomposition")
    return auxiliary, colors, decompositions


def _two_isomorphism_component_forms(auxiliary, colors, algorithm,
                                     certificate):
    """
    Canonicalize every nontrivial block gadget separately.

    Blocks are direct summands of a cycle matroid.  Treating their gadgets as
    a multiset avoids introducing a large symmetric disconnected graph when a
    graph has many equal blocks.
    """
    color_of = {vertex: i
                for i, color in enumerate(_TWO_ISOMORPHISM_COLORS)
                for vertex in colors[color]}
    forms = {}
    for vertices in auxiliary.connected_components(sort=False):
        # Singletons do not need graph canonization.  Loops and coloops are
        # kept out of ``auxiliary`` altogether and matched directly by
        # ``is_2isomorphic``.
        if len(vertices) == 1:
            continue
        cells = [[] for _ in _TWO_ISOMORPHISM_COLORS]
        for vertex in vertices:
            cells[color_of[vertex]].append(vertex)
        sizes = tuple(map(len, cells))
        real_edges = [vertex for vertex in vertices
                      if vertex[0] == 'edge']

        # A block consisting of one Q, S, or P component represents a coloop,
        # circuit, or parallel class.  Its colored isomorphism type depends
        # only on the component type and number of real elements, so avoid a
        # graph-canonization call for this common case.
        if (not sizes[3] and not sizes[7] and not sizes[8] and
                sum(sizes[4:7]) == 1):
            if certificate:
                mapping = {
                    vertex: (color, i)
                    for color, cell in enumerate(cells)
                    for i, vertex in enumerate(cell)
                }
                entry = mapping, list(vertices), real_edges
            else:
                entry = None
            forms.setdefault((sizes, None), []).append(entry)
            continue

        partition = [cell for cell in cells if cell]
        component = auxiliary.subgraph(vertices)
        if certificate:
            canonical, mapping = component.canonical_label(
                partition=partition, algorithm=algorithm,
                certificate=True, immutable=True)
            entry = mapping, list(vertices), real_edges
        else:
            canonical = component.canonical_label(
                partition=partition, algorithm=algorithm, immutable=True)
            entry = None
        forms.setdefault((sizes, canonical), []).append(entry)
    return forms


def _two_isomorphism_normalized_graph(G, target_to_source=None):
    """
    Return a mutable copy with integer vertices and occurrence edge labels.
    """
    from sage.graphs.graph import Graph

    vertices = list(G)
    vertex_to_int = {vertex: i for i, vertex in enumerate(vertices)}
    edges = []
    for i, (u, v, _) in enumerate(G.edge_iterator()):
        edge = i if target_to_source is None else target_to_source[i]
        edges.append((vertex_to_int[u], vertex_to_int[v], edge))
    return Graph([range(len(vertices)), edges], format='vertices_and_edges',
                 loops=True, multiedges=True)


def _two_isomorphism_normalized_edge_table(G, target_to_source=None):
    """
    Return normalized edge endpoints and isolated vertex identifiers.
    """
    vertices = list(G)
    vertex_to_int = {vertex: i for i, vertex in enumerate(vertices)}
    degrees = [0] * len(vertices)
    table = {}
    for i, (source_u, source_v, _) in enumerate(G.edge_iterator()):
        edge = i if target_to_source is None else target_to_source[i]
        u = vertex_to_int[source_u]
        v = vertex_to_int[source_v]
        table[edge] = (u, v)
        degrees[u] += 1
        degrees[v] += 1
    isolates = [vertex for vertex, degree in enumerate(degrees) if not degree]
    return table, isolates, len(vertices)


def _two_isomorphism_edge_table(G):
    """
    Return the endpoints indexed by the private edge occurrence labels.
    """
    table = {}
    for u, v, edge in G.edge_iterator():
        if edge in table:
            raise RuntimeError("edge occurrence labels must be unique")
        table[edge] = (u, v)
    return table


def _two_isomorphism_graph_state(G):
    """
    Return mutable endpoint and incidence indices for certificate replay.
    """
    table = _two_isomorphism_edge_table(G)
    return _two_isomorphism_endpoint_state(table, set(G))


def _two_isomorphism_endpoint_state(table, vertices):
    """
    Return mutable endpoint and incidence indices from normalized data.
    """
    table = dict(table)
    vertices = set(vertices)
    incidence = {vertex: set() for vertex in vertices}
    for edge, endpoints in table.items():
        for vertex in set(endpoints):
            incidence[vertex].add(edge)
    return {'table': table, 'incidence': incidence, 'vertices': vertices}


def _two_isomorphism_update_edge_state(state, edge, old, new):
    """
    Update endpoint and incidence indices after rewiring one edge.
    """
    if state is None:
        return
    edge_atoms = state.get('edge_atoms')
    atom = edge_atoms[edge] if edge_atoms is not None else None
    for vertex in set(old):
        state['incidence'][vertex].remove(edge)
        if edge_atoms is not None:
            atom_edges = state['atom_incidence'][vertex][atom]
            atom_edges.remove(edge)
            if not atom_edges:
                del state['atom_incidence'][vertex][atom]
    for vertex in set(new):
        state['incidence'].setdefault(vertex, set()).add(edge)
        if edge_atoms is not None:
            state['atom_incidence'].setdefault(vertex, {}).setdefault(
                atom, set()).add(edge)
    state['table'][edge] = new


def _two_isomorphism_edge_subgraph(G, edge_ids, table=None):
    """
    Return the edge-induced subgraph on ``edge_ids`` including its endpoints.
    """
    from sage.graphs.graph import Graph

    if table is None:
        table = _two_isomorphism_edge_table(G)
    edges = [(table[edge][0], table[edge][1], edge) for edge in edge_ids]
    vertices = {vertex for u, v, _ in edges for vertex in (u, v)}
    return Graph([vertices, edges], format='vertices_and_edges',
                 loops=True, multiedges=True)


def _two_isomorphism_fixed_edge_vertex_mapping(G, H, table_G=None,
                                                table_H=None):
    """
    Return the vertex map induced by equal unique edge labels, if it exists.

    Every nontrivial connected component has only two possible orientations.
    Propagating either orientation along its labeled edges is linear in the
    size of the component and avoids general graph canonization.
    """
    if table_G is None:
        table_G = _two_isomorphism_edge_table(G)
    if table_H is None:
        table_H = _two_isomorphism_edge_table(H)
    if set(table_G) != set(table_H):
        return None
    adjacency = {}
    for edge, endpoints in table_G.items():
        for vertex in set(endpoints):
            adjacency.setdefault(vertex, []).append(edge)
    target_vertices = {vertex for endpoints in table_H.values()
                       for vertex in endpoints}
    mapping = {}
    inverse = {}

    for start, incident in adjacency.items():
        if start in mapping:
            continue
        loop = next((edge for edge in incident
                     if table_G[edge][0] == table_G[edge][1]), None)
        if loop is not None:
            target_u, target_v = table_H[loop]
            if target_u != target_v:
                return None
            candidates = [target_u]
        else:
            target_u, target_v = table_H[incident[0]]
            if target_u == target_v:
                return None
            candidates = [target_u, target_v]

        found = None
        for image in candidates:
            local = {}
            local_inverse = {}
            stack = [(start, image)]
            valid = True
            while stack and valid:
                vertex, target_vertex = stack.pop()
                if vertex in local:
                    valid = local[vertex] == target_vertex
                    continue
                if (target_vertex in local_inverse and
                        local_inverse[target_vertex] != vertex):
                    valid = False
                    break
                if (target_vertex in inverse and
                        inverse[target_vertex] != vertex):
                    valid = False
                    break
                local[vertex] = target_vertex
                local_inverse[target_vertex] = vertex
                for edge in adjacency[vertex]:
                    u, v = table_G[edge]
                    target_u, target_v = table_H[edge]
                    if u == v:
                        if target_u != target_v or target_vertex != target_u:
                            valid = False
                            break
                        continue
                    if target_u == target_v:
                        valid = False
                        break
                    if target_vertex == target_u:
                        other_target = target_v
                    elif target_vertex == target_v:
                        other_target = target_u
                    else:
                        valid = False
                        break
                    other = v if u == vertex else u
                    stack.append((other, other_target))
            if valid:
                found = local
                break
        if found is None:
            return None
        mapping.update(found)
        inverse.update((target, source) for source, target in found.items())

    if set(mapping.values()) != target_vertices:
        return None
    return mapping


def _two_isomorphism_atoms(G):
    """
    Return the cycle-matroid components of an occurrence-labeled graph.

    Loops and coloops are individual atoms.  Every other atom is a nontrivial
    block.  The atoms are sorted by their smallest edge occurrence identifier.
    """
    return _two_isomorphism_atoms_from_table(_two_isomorphism_edge_table(G))


def _two_isomorphism_atoms_from_table(table):
    """
    Return cycle-matroid components from an occurrence endpoint table.

    The block decomposition only needs the underlying simple graph.  Keeping
    parallel occurrences in a side table avoids constructing a large mutable
    multigraph merely to discover that they belong to the same block.
    """
    from sage.graphs.graph import Graph

    table = dict(table)
    atoms = [{edge} for edge, (u, v) in table.items() if u == v]
    endpoint_edges = {}
    for edge, (u, v) in table.items():
        if u == v:
            continue
        endpoints = (u, v) if u < v else (v, u)
        endpoint_edges.setdefault(endpoints, set()).add(edge)
    active_vertices = {vertex for endpoints in endpoint_edges
                       for vertex in endpoints}
    work = Graph([active_vertices, endpoint_edges],
                 format='vertices_and_edges', loops=False, multiedges=False)
    blocks = work.blocks_and_cut_vertices()[0] if work else []
    vertex_blocks = {}
    for block_index, vertices in enumerate(blocks):
        for vertex in vertices:
            vertex_blocks.setdefault(vertex, []).append(block_index)
    vertex_block_sets = {vertex: set(indices)
                         for vertex, indices in vertex_blocks.items()}
    block_edges = [set() for _ in blocks]
    for (u, v), edges in endpoint_edges.items():
        blocks_u = vertex_blocks[u]
        blocks_v = vertex_blocks[v]
        if len(blocks_u) == 1:
            block_index = blocks_u[0]
        elif len(blocks_v) == 1:
            block_index = blocks_v[0]
        else:
            common = vertex_block_sets[u].intersection(vertex_block_sets[v])
            if len(common) != 1:
                raise RuntimeError("an edge must belong to exactly one block")
            block_index = common.pop()
        block_edges[block_index].update(edges)
    atoms.extend(block_edges)
    atoms = [frozenset(atom) for atom in atoms if atom]
    covered_edges = set().union(*atoms) if atoms else set()
    if covered_edges != set(table):
        raise RuntimeError("some edges are missing from the atom decomposition")
    return sorted(atoms, key=lambda atom: min(atom))


def _two_isomorphism_twist(G, edge_ids, operations=None, state=None,
                           mutate=True):
    """
    Apply one Whitney twist and return its two boundary vertices.
    """
    side = frozenset(edge_ids)
    table = state['table'] if state is not None else _two_isomorphism_edge_table(G)
    if not side:
        raise ValueError("a Whitney twist needs two nonempty edge sides")
    if any(edge not in table for edge in side):
        raise ValueError("a Whitney twist references an unknown edge")
    if len(side) == len(table):
        raise ValueError("a Whitney twist needs two nonempty edge sides")
    side_incidence = {}
    for edge in side:
        for vertex in set(table[edge]):
            side_incidence[vertex] = side_incidence.get(vertex, 0) + 1
    if state is None:
        total_incidence = {}
        for endpoints in table.values():
            for vertex in set(endpoints):
                total_incidence[vertex] = total_incidence.get(vertex, 0) + 1
        boundary = {vertex for vertex, count in side_incidence.items()
                    if total_incidence[vertex] > count}
    else:
        boundary = {vertex for vertex, count in side_incidence.items()
                    if len(state['incidence'][vertex]) > count}
    if len(boundary) != 2:
        raise ValueError("the two sides of a Whitney twist must share two vertices")
    a, b = sorted(boundary)
    for edge in sorted(side):
        u, v = table[edge]
        new_u = b if u == a else a if u == b else u
        new_v = b if v == a else a if v == b else v
        if mutate:
            G.delete_edge(u, v, edge)
            G.add_edge(new_u, new_v, edge)
        elif state is None:
            raise RuntimeError("state-only replay needs endpoint indices")
        _two_isomorphism_update_edge_state(
            state, edge, (u, v), (new_u, new_v))
    if operations is not None:
        operations.append({
            'operation': 'whitney_twist',
            'vertices': (a, b),
            'edges': tuple(sorted(side)),
        })
    return a, b


def _two_isomorphism_split_is_valid(G, vertex, edge_ids, edge_atoms=None,
                                    state=None):
    """
    Test whether moving ``edge_ids`` is a valid vertex cleaving.
    """
    move = set(edge_ids)
    table = state['table'] if state is not None else _two_isomorphism_edge_table(G)
    if state is None:
        incident = {edge for edge, endpoints in table.items()
                    if vertex in endpoints}
    else:
        incident = state['incidence'][vertex]
    if not move or not move < incident:
        return False

    if edge_atoms is not None:
        if state is not None and 'atom_incidence' in state:
            atom_incidence = state['atom_incidence'][vertex]
            for edge in move:
                atom = edge_atoms[edge]
                if not atom_incidence.get(atom, set()).issubset(move):
                    return False
            return True
        groups = {}
        for edge in incident:
            groups.setdefault(edge_atoms[edge], set()).add(edge)
        return all(group.issubset(move) or group.isdisjoint(move)
                   for group in groups.values())

    remainder = G.copy(immutable=False)
    remainder.delete_vertex(vertex)
    component = {}
    for i, vertices in enumerate(remainder.connected_components(sort=False)):
        component.update((v, i) for v in vertices)
    groups = {}
    for edge in incident:
        u, v = table[edge]
        if u == v == vertex:
            key = ('loop', edge)
        else:
            other = v if u == vertex else u
            key = ('component', component[other])
        groups.setdefault(key, set()).add(edge)
    return all(group.issubset(move) or group.isdisjoint(move)
               for group in groups.values())


def _two_isomorphism_split(G, vertex, new_vertex, edge_ids, operations=None,
                           check=True, edge_atoms=None, state=None,
                           mutate=True):
    """
    Cleave ``vertex`` by moving complete edge incidences to ``new_vertex``.
    """
    move = tuple(sorted(edge_ids))
    vertices = G if state is None else state['vertices']
    if new_vertex in vertices:
        raise ValueError("the new vertex of a cleaving must be fresh")
    if check and not _two_isomorphism_split_is_valid(
            G, vertex, move, edge_atoms=edge_atoms, state=state):
        raise ValueError("the edge partition is not a valid vertex cleaving")
    table = state['table'] if state is not None else _two_isomorphism_edge_table(G)
    if mutate:
        G.add_vertex(new_vertex)
    elif state is None:
        raise RuntimeError("state-only replay needs endpoint indices")
    if state is not None:
        state['vertices'].add(new_vertex)
        state['incidence'][new_vertex] = set()
        if 'atom_incidence' in state:
            state['atom_incidence'][new_vertex] = {}
    for edge in move:
        try:
            u, v = table[edge]
        except KeyError:
            raise ValueError("a vertex cleaving references an unknown edge") from None
        if vertex not in (u, v):
            raise ValueError("a moved edge is not incident with the cleaved vertex")
        new_u = new_vertex if u == vertex else u
        new_v = new_vertex if v == vertex else v
        if mutate:
            G.delete_edge(u, v, edge)
            G.add_edge(new_u, new_v, edge)
        _two_isomorphism_update_edge_state(
            state, edge, (u, v), (new_u, new_v))
    if operations is not None:
        operations.append({
            'operation': 'vertex_cleaving',
            'vertex': vertex,
            'new_vertex': new_vertex,
            'edges': move,
        })


def _two_isomorphism_identify(G, keep, drop, operations=None, check=True,
                              state=None, mutate=True):
    """
    Identify vertices in different connected components.
    """
    vertices = G if state is None else state['vertices']
    if keep == drop or keep not in vertices or drop not in vertices:
        raise ValueError("vertex identification needs two distinct vertices")
    if check:
        if not mutate:
            raise RuntimeError("state-only identification needs a component check")
        if drop in G.connected_component_containing_vertex(keep, sort=False):
            raise ValueError(
                "identified vertices must lie in different components")
    table = state['table'] if state is not None else _two_isomorphism_edge_table(G)
    if state is None:
        incident = [edge for edge, endpoints in table.items()
                    if drop in endpoints]
    else:
        incident = list(state['incidence'][drop])
    for edge in incident:
        u, v = table[edge]
        new_u = keep if u == drop else u
        new_v = keep if v == drop else v
        if mutate:
            G.delete_edge(u, v, edge)
            G.add_edge(new_u, new_v, edge)
        _two_isomorphism_update_edge_state(
            state, edge, (u, v), (new_u, new_v))
    if mutate:
        G.delete_vertex(drop)
    elif state is None:
        raise RuntimeError("state-only replay needs endpoint indices")
    if state is not None:
        state['vertices'].remove(drop)
        state['incidence'].pop(drop)
        if 'atom_incidence' in state:
            state['atom_incidence'].pop(drop)
    if operations is not None:
        operations.append({
            'operation': 'vertex_identification',
            'keep': keep,
            'drop': drop,
        })


def _two_isomorphism_cycle_order(skeleton):
    """
    Return the edge labels in cyclic order in an S skeleton.
    """
    table = _two_isomorphism_edge_table(skeleton)
    if len(table) < 3:
        raise RuntimeError("an S component must contain a cycle")
    adjacency = {vertex: [] for vertex in skeleton}
    for edge, (u, v) in table.items():
        adjacency[u].append(edge)
        adjacency[v].append(edge)
    if any(len(edges) != 2 for edges in adjacency.values()):
        raise RuntimeError("an S component must be a cycle")
    first = min(table, key=repr)
    start, current = table[first]
    order = [first]
    used = {first}
    while len(order) < len(table):
        candidates = [edge for edge in adjacency[current] if edge not in used]
        if len(candidates) != 1:
            raise RuntimeError("cannot traverse an S component")
        edge = candidates[0]
        order.append(edge)
        used.add(edge)
        u, v = table[edge]
        current = v if u == current else u
    if current != start:
        raise RuntimeError("an S component is not cyclic")
    return order


def _two_isomorphism_reversal_plan(source, target):
    """
    Return interval reversals changing one anchored cyclic word into another.
    """
    current = list(source)
    moves = []
    for i in range(1, len(current) - 1):
        j = current.index(target[i], i)
        if i != j:
            moves.append((i, j))
            current[i:j + 1] = reversed(current[i:j + 1])
    if current != list(target):
        raise RuntimeError("failed to align an S component")
    return moves


def _two_isomorphism_mapped_skeleton(component, inverse_auxiliary_mapping):
    """
    Relabel target skeleton edges by their matched source auxiliary vertices.
    """
    from sage.graphs.graph import Graph

    skeleton = component['skeleton']
    edges = [(u, v, inverse_auxiliary_mapping[label])
             for u, v, label in skeleton.edge_iterator()]
    return Graph([list(skeleton), edges], format='vertices_and_edges',
                 loops=False, multiedges=True)


def _two_isomorphism_align_block(G, source, target, auxiliary_mapping,
                                 source_vertex_map, operations, state=None):
    """
    Transform one nontrivial source block into its matched target block.
    """
    source_components = [{
        'kind': component['kind'],
        'anchor': component['anchor'],
        'skeleton': component['skeleton'].copy(immutable=False),
    } for component in source['components']]
    target_components = target['components']
    inverse_auxiliary_mapping = {
        target_vertex: source_vertex
        for source_vertex, target_vertex in auxiliary_mapping.items()
    }
    component_match = {}
    for i, component in enumerate(source_components):
        target_anchor = auxiliary_mapping[component['anchor']]
        component_match[i] = target_anchor[2]
        if target_components[target_anchor[2]]['anchor'] != target_anchor:
            raise RuntimeError("invalid SPQR component certificate")
        if component['kind'] != target_components[target_anchor[2]]['kind']:
            raise RuntimeError("matched SPQR components have different types")

    tree_adjacency = {i: [] for i in range(len(source_components))}
    port_neighbor = {}
    for left, right, left_item, right_item in source['tree_edges']:
        tree_adjacency[left].append((right, left_item, right_item))
        tree_adjacency[right].append((left, right_item, left_item))
        port_neighbor[left, left_item] = right
        port_neighbor[right, right_item] = left

    direct_real_edges = []
    for component in source_components:
        direct_real_edges.append({label[1]
                                  for _, _, label in component['skeleton'].edge_iterator()
                                  if label[0] == 'edge'})

    def component_side(component, item):
        neighbor = port_neighbor[component, item]
        side = set()
        stack = [(neighbor, component)]
        while stack:
            vertex, parent = stack.pop()
            side.add(vertex)
            stack.extend((other, vertex)
                         for other, _, _ in tree_adjacency[vertex]
                         if other != parent)
        return side

    expansion = {}

    def item_expansion(component, item):
        """
        Return the real edges and components represented by one skeleton item.

        Most virtual ports never participate in an S-node reversal.  Compute
        their SPQR-tree sides lazily instead of materializing both sides of
        every tree edge up front.
        """
        key = component, item
        try:
            return expansion[key]
        except KeyError:
            if item[0] == 'edge':
                value = {item[1]}, set()
            else:
                side = component_side(component, item)
                edges = set().union(
                    *(direct_real_edges[j] for j in side))
                value = edges, side
            expansion[key] = value
            return value

    embeddings = [{vertex: source_vertex_map[vertex]
                   for vertex in component['skeleton']}
                  for component in source_components]

    # In an S component an arbitrary permutation of pieces need not be a
    # dihedral graph automorphism.  Anchored interval reversals are Whitney
    # twists and generate every permutation.
    for i, component in enumerate(source_components):
        if component['kind'] != 'S':
            continue
        source_order = _two_isomorphism_cycle_order(component['skeleton'])
        target_component = target_components[component_match[i]]
        target_order = [inverse_auxiliary_mapping[item]
                        for item in _two_isomorphism_cycle_order(
                            target_component['skeleton'])]
        anchor = source_order[0]
        anchor_index = target_order.index(anchor)
        target_order = (target_order[anchor_index:] +
                        target_order[:anchor_index])
        reflected = [target_order[0]] + list(reversed(target_order[1:]))
        plans = [(_two_isomorphism_reversal_plan(source_order, target_order),
                  target_order),
                 (_two_isomorphism_reversal_plan(source_order, reflected),
                  reflected)]
        plan, desired = min(plans, key=lambda item: len(item[0]))
        current_order = list(source_order)
        for first, last in plan:
            interval = current_order[first:last + 1]
            edge_side = set().union(
                *(item_expansion(i, item)[0] for item in interval))
            affected_components = set().union(
                *(item_expansion(i, item)[1] for item in interval))
            abstract_boundary = _two_isomorphism_twist(
                component['skeleton'], interval)
            expected_boundary = {embeddings[i][vertex]
                                 for vertex in abstract_boundary}
            actual_boundary = _two_isomorphism_twist(
                G, edge_side, operations, state=state)
            if set(actual_boundary) != expected_boundary:
                raise RuntimeError("an S-component twist has wrong boundary")
            a, b = actual_boundary
            for j in affected_components:
                embeddings[j] = {
                    vertex: b if image == a else a if image == b else image
                    for vertex, image in embeddings[j].items()
                }
            current_order[first:last + 1] = reversed(
                current_order[first:last + 1])
        if current_order != desired:
            raise RuntimeError("failed to execute S-component reversals")

    local_maps = {}
    for i, component in enumerate(source_components):
        target_skeleton = _two_isomorphism_mapped_skeleton(
            target_components[component_match[i]],
            inverse_auxiliary_mapping)
        mapping = _two_isomorphism_fixed_edge_vertex_mapping(
            component['skeleton'], target_skeleton)
        if mapping is None:
            raise RuntimeError("matched SPQR skeletons are not isomorphic")
        local_maps[i] = mapping

    # Root the SPQR tree.  A crossed pair of virtual ports is repaired by one
    # twist of the complete descendant side.
    if source_components:
        parent = {0: None}
        parent_port = {}
        order = [0]
        for vertex in order:
            for other, item, other_item in tree_adjacency[vertex]:
                if other in parent:
                    continue
                parent[other] = vertex
                parent_port[other] = (item, other_item)
                order.append(other)

        for child in order[1:]:
            parent_vertex = parent[child]
            parent_item, child_item = parent_port[child]
            parent_endpoints = _two_isomorphism_edge_table(
                source_components[parent_vertex]['skeleton'])[parent_item]
            child_endpoints = _two_isomorphism_edge_table(
                source_components[child]['skeleton'])[child_item]
            parent_targets = {
                embeddings[parent_vertex][vertex]:
                local_maps[parent_vertex][vertex]
                for vertex in parent_endpoints
            }
            child_targets = {
                embeddings[child][vertex]: local_maps[child][vertex]
                for vertex in child_endpoints
            }
            if set(parent_targets) != set(child_targets):
                raise RuntimeError("SPQR ports are not glued at the same vertices")
            if set(parent_targets.values()) != set(child_targets.values()):
                raise RuntimeError("matched target SPQR ports have different endpoints")
            if parent_targets != child_targets:
                side = set()
                stack = [(child, parent_vertex)]
                while stack:
                    vertex, previous = stack.pop()
                    side.add(vertex)
                    stack.extend((other, vertex)
                                 for other, _, _ in tree_adjacency[vertex]
                                 if other != previous)
                edge_side = set().union(*(direct_real_edges[j] for j in side))
                boundary = _two_isomorphism_twist(
                    G, edge_side, operations, state=state)
                if set(boundary) != set(parent_targets):
                    raise RuntimeError("an SPQR port twist has wrong boundary")
                a, b = boundary
                for j in side:
                    embeddings[j] = {
                        vertex: b if image == a else a if image == b else image
                        for vertex, image in embeddings[j].items()
                    }
                child_targets = {
                    embeddings[child][vertex]: local_maps[child][vertex]
                    for vertex in child_endpoints
                }
                if parent_targets != child_targets:
                    raise RuntimeError("failed to align an SPQR virtual port")


def _two_isomorphism_replay_certificate(G, H, witness):
    """
    Replay and verify a Whitney-operation certificate.
    """
    if not isinstance(witness, dict) or witness.get('version') != 1:
        return False
    edge_mapping = witness.get('edge_mapping')
    if not isinstance(edge_mapping, dict):
        return False
    operations = witness.get('operations')
    if not isinstance(operations, (list, tuple)):
        return False
    size = G.size()
    if (set(edge_mapping) != set(range(size)) or
            set(edge_mapping.values()) != set(range(H.size()))):
        return False
    target_to_source = {target: source
                        for source, target in edge_mapping.items()}
    current_table, _, current_order = _two_isomorphism_normalized_edge_table(G)
    target_table, _, target_order = _two_isomorphism_normalized_edge_table(
        H, target_to_source)
    state = _two_isomorphism_endpoint_state(
        current_table, range(current_order))
    atoms = _two_isomorphism_atoms_from_table(current_table)
    edge_atoms = {edge: atom_index
                  for atom_index, atom in enumerate(atoms)
                  for edge in atom}
    state['edge_atoms'] = edge_atoms
    state['atom_incidence'] = {vertex: {} for vertex in state['vertices']}
    for vertex, edges in state['incidence'].items():
        for edge in edges:
            state['atom_incidence'][vertex].setdefault(
                edge_atoms[edge], set()).add(edge)
    phases = {
        'delete_isolated_vertex': 0,
        'vertex_cleaving': 1,
        'whitney_twist': 2,
        'vertex_identification': 3,
        'add_isolated_vertex': 4,
    }
    phase = -1
    component_parent = None
    component_of = None

    def component_root(component):
        while component_parent[component] != component:
            component_parent[component] = component_parent[
                component_parent[component]]
            component = component_parent[component]
        return component

    try:
        for step in operations:
            if not isinstance(step, dict):
                return False
            operation = step['operation']
            step_phase = phases[operation]
            if step_phase < phase:
                return False
            phase = step_phase
            if operation == 'vertex_cleaving':
                if not isinstance(step.get('edges'), (list, tuple)):
                    return False
                _two_isomorphism_split(
                    None, step['vertex'], step['new_vertex'],
                    step['edges'], check=True, edge_atoms=edge_atoms,
                    state=state, mutate=False)
            elif operation == 'vertex_identification':
                if component_parent is None:
                    vertex_parent = {vertex: vertex
                                     for vertex in state['vertices']}

                    def vertex_root(vertex):
                        while vertex_parent[vertex] != vertex:
                            vertex_parent[vertex] = vertex_parent[
                                vertex_parent[vertex]]
                            vertex = vertex_parent[vertex]
                        return vertex

                    for u, v in state['table'].values():
                        root_u = vertex_root(u)
                        root_v = vertex_root(v)
                        if root_u != root_v:
                            vertex_parent[root_v] = root_u
                    roots = {}
                    component_of = {}
                    for vertex in state['vertices']:
                        root = vertex_root(vertex)
                        component_of[vertex] = roots.setdefault(
                            root, len(roots))
                    component_parent = list(range(len(roots)))
                keep = step['keep']
                drop = step['drop']
                keep_component = component_root(component_of[keep])
                drop_component = component_root(component_of[drop])
                if keep_component == drop_component:
                    return False
                component_parent[drop_component] = keep_component
                _two_isomorphism_identify(
                    None, keep, drop, check=False, state=state,
                    mutate=False)
                component_of.pop(drop)
            elif operation == 'whitney_twist':
                if (not isinstance(step.get('edges'), (list, tuple)) or
                        not isinstance(step.get('vertices'), (list, tuple)) or
                        len(step['vertices']) != 2):
                    return False
                boundary = _two_isomorphism_twist(
                    None, step['edges'], state=state, mutate=False)
                if set(boundary) != set(step['vertices']):
                    return False
            elif operation == 'delete_isolated_vertex':
                vertex = step['vertex']
                if (vertex not in state['vertices'] or
                        state['incidence'][vertex]):
                    return False
                state['vertices'].remove(vertex)
                state['incidence'].pop(vertex)
                state['atom_incidence'].pop(vertex)
            elif operation == 'add_isolated_vertex':
                vertex = step['vertex']
                if vertex in state['vertices']:
                    return False
                state['vertices'].add(vertex)
                state['incidence'][vertex] = set()
                state['atom_incidence'][vertex] = {}
            else:
                return False
    except (IndexError, KeyError, TypeError, ValueError, RuntimeError):
        return False

    vertex_mapping = witness.get('vertex_mapping')
    if (not isinstance(vertex_mapping, dict) or
            set(vertex_mapping) != state['vertices'] or
            set(vertex_mapping.values()) != set(range(target_order)) or
            len(vertex_mapping) != target_order):
        return False
    current_table = state['table']
    if set(current_table) != set(target_table):
        return False
    for edge, (u, v) in current_table.items():
        target_u, target_v = target_table[edge]
        if {vertex_mapping[u], vertex_mapping[v]} != {target_u, target_v}:
            return False
        if (u == v) != (target_u == target_v):
            return False
    return True


def verify_2isomorphism_certificate(G, H, witness):
    """
    Verify and replay a certificate returned by
    :func:`~sage.graphs.morphisms.is_2isomorphic`.

    INPUT:

    - ``G``, ``H`` -- undirected graphs

    - ``witness`` -- a Whitney-operation certificate

    OUTPUT: boolean

    The input graphs are not modified.  Every vertex cleaving, Whitney twist,
    and vertex identification is checked before it is replayed.  Finally the
    edge and vertex mappings are checked against ``H``.

    EXAMPLES::

        sage: G = graphs.PathGraph(5)
        sage: H = graphs.StarGraph(4)
        sage: ok, witness = G.is_2isomorphic(H, certificate=True)
        sage: ok and G.verify_2isomorphism_certificate(H, witness)
        True

    A modified operation is rejected::

        sage: witness['operations'][0]['edges'] = (99,)
        sage: G.verify_2isomorphism_certificate(H, witness)
        False

    Malformed mapping data is rejected as well::

        sage: witness['edge_mapping'][0] = []
        sage: G.verify_2isomorphism_certificate(H, witness)
        False

    The final vertex map must be a bijection; it cannot stand in for a missing
    vertex identification::

        sage: G = Graph([(0, 1), (2, 3)])
        sage: H = graphs.PathGraph(3)
        sage: fake = {'version': 1,
        ....:         'edge_mapping': {0: 0, 1: 1},
        ....:         'operations': [],
        ....:         'vertex_mapping': {0: 0, 1: 1, 2: 1, 3: 2}}
        sage: G.verify_2isomorphism_certificate(H, fake)
        False
    """
    from sage.graphs.graph import Graph

    if not isinstance(G, Graph) or not isinstance(H, Graph):
        raise TypeError("can only verify 2-isomorphism between undirected graphs")
    try:
        return _two_isomorphism_replay_certificate(G, H, witness)
    except (IndexError, KeyError, TypeError, ValueError, RuntimeError):
        # A verifier must reject malformed external data rather than expose
        # implementation errors from, for example, unhashable mapping values.
        return False


def _two_isomorphism_operation_certificate(G, H, edge_mapping,
                                           decompositions_G=(),
                                           decompositions_H=(),
                                           auxiliary_mapping=None):
    """
    Compile an edge bijection into explicit Whitney operations.
    """
    auxiliary_mapping = auxiliary_mapping or {}
    target_to_source = {target: source
                        for source, target in edge_mapping.items()}
    source_table, source_isolates, source_order = (
        _two_isomorphism_normalized_edge_table(G))
    target_table, target_isolates, _ = (
        _two_isomorphism_normalized_edge_table(H, target_to_source))
    direct_vertex_map = _two_isomorphism_fixed_edge_vertex_mapping(
        None, None, source_table, target_table)
    if direct_vertex_map is not None:
        operations = [{
            'operation': 'delete_isolated_vertex',
            'vertex': vertex,
        } for vertex in source_isolates]
        next_vertex = source_order
        for target_vertex in target_isolates:
            operations.append({
                'operation': 'add_isolated_vertex',
                'vertex': next_vertex,
            })
            direct_vertex_map[next_vertex] = target_vertex
            next_vertex += 1
        return {
            'version': 1,
            'edge_mapping': dict(edge_mapping),
            'operations': operations,
            'vertex_mapping': direct_vertex_map,
        }

    current = _two_isomorphism_normalized_graph(G)
    target = _two_isomorphism_normalized_graph(H, target_to_source)
    current_state = _two_isomorphism_graph_state(current)
    operations = []
    next_vertex = max(current, default=-1) + 1

    for vertex in sorted([v for v in current if not current.degree(v)]):
        current.delete_vertex(vertex)
        current_state['vertices'].remove(vertex)
        current_state['incidence'].pop(vertex)
        operations.append({
            'operation': 'delete_isolated_vertex',
            'vertex': vertex,
        })

    source_atoms = _two_isomorphism_atoms(current)
    target_atoms = _two_isomorphism_atoms(target)
    if set(source_atoms) != set(target_atoms):
        raise RuntimeError("the edge map does not preserve matroid components")

    atom_vertex_maps = []
    edge_table = current_state['table']
    for atom in source_atoms:
        vertices = {vertex for edge in atom for vertex in edge_table[edge]}
        atom_vertex_maps.append({vertex: vertex for vertex in vertices})

    vertex_atoms = {}
    for atom_index, atom in enumerate(source_atoms):
        for edge in atom:
            for vertex in edge_table[edge]:
                vertex_atoms.setdefault(vertex, set()).add(atom_index)
    for vertex in sorted(vertex_atoms):
        atom_indices = sorted(vertex_atoms[vertex])
        for atom_index in atom_indices[1:]:
            move = [edge for edge in source_atoms[atom_index]
                    if vertex in current_state['table'][edge]]
            _two_isomorphism_split(
                current, vertex, next_vertex, move, operations, check=False,
                state=current_state)
            atom_vertex_maps[atom_index][vertex] = next_vertex
            next_vertex += 1

    source_decomposition = {
        decomposition['edge_ids']: decomposition
        for decomposition in decompositions_G
    }
    target_decomposition = {
        frozenset(target_to_source[edge]
                  for edge in decomposition['edge_ids']): decomposition
        for decomposition in decompositions_H
    }

    atom_maps = []
    for atom_index, atom in enumerate(source_atoms):
        if len(atom) == 1:
            edge = next(iter(atom))
            source_u, source_v = current_state['table'][edge]
            target_u, target_v = target_table[edge]
            if (source_u == source_v) != (target_u == target_v):
                raise RuntimeError("an edge map confuses a loop and a coloop")
            if source_u == source_v:
                atom_maps.append({source_u: target_u})
            else:
                atom_maps.append({source_u: target_u, source_v: target_v})
            continue
        source_graph = _two_isomorphism_edge_subgraph(
            current, atom, current_state['table'])
        target_graph = _two_isomorphism_edge_subgraph(
            target, atom, target_table)
        vertex_map = _two_isomorphism_fixed_edge_vertex_mapping(
            source_graph, target_graph)
        if vertex_map is None:
            try:
                source_record = source_decomposition[atom]
                target_record = target_decomposition[atom]
            except KeyError:
                raise RuntimeError(
                    "missing SPQR data for a non-isomorphic block pair") from None
            _two_isomorphism_align_block(
                current, source_record, target_record, auxiliary_mapping,
                atom_vertex_maps[atom_index], operations,
                state=current_state)
            source_graph = _two_isomorphism_edge_subgraph(
                current, atom, current_state['table'])
            vertex_map = _two_isomorphism_fixed_edge_vertex_mapping(
                source_graph, target_graph)
            if vertex_map is None:
                raise RuntimeError("Whitney twists did not align a block")
        atom_maps.append(vertex_map)

    current_to_target = {}
    for vertex_map in atom_maps:
        for source_vertex, target_vertex in vertex_map.items():
            if (source_vertex in current_to_target and
                    current_to_target[source_vertex] != target_vertex):
                raise RuntimeError("inconsistent local block isomorphisms")
            current_to_target[source_vertex] = target_vertex

    target_to_current = {}
    for source_vertex, target_vertex in current_to_target.items():
        target_to_current.setdefault(target_vertex, []).append(source_vertex)
    for target_vertex in sorted(target_to_current):
        vertices = target_to_current[target_vertex]
        keep = vertices[0]
        for drop in vertices[1:]:
            _two_isomorphism_identify(
                current, keep, drop, operations, check=False,
                state=current_state)
            current_to_target.pop(drop)

    target_isolates = sorted(vertex for vertex in target
                             if not target.degree(vertex))
    for target_vertex in target_isolates:
        current.add_vertex(next_vertex)
        current_state['vertices'].add(next_vertex)
        current_state['incidence'][next_vertex] = set()
        operations.append({
            'operation': 'add_isolated_vertex',
            'vertex': next_vertex,
        })
        current_to_target[next_vertex] = target_vertex
        next_vertex += 1

    witness = {
        'version': 1,
        'edge_mapping': dict(edge_mapping),
        'operations': operations,
        'vertex_mapping': current_to_target,
    }
    return witness


def is_2isomorphic(G, H, certificate=False):
    r"""
    Test whether ``G`` and ``H`` are Whitney 2-isomorphic.

    Two graphs are *2-isomorphic* if their cycle matroids are isomorphic.  In
    other words, there is a bijection between their edge occurrences that
    preserves cycles.  Edge labels and isolated vertices are ignored.

    INPUT:

    - ``G``, ``H`` -- undirected graphs, possibly with loops and multiple edges

    - ``certificate`` -- boolean (default: ``False``); whether to return a
      constructive Whitney-operation certificate

    OUTPUT:

    If ``certificate`` is ``False``, return a boolean.  Otherwise return a pair
    ``(result, witness)``.  On failure, ``witness`` is ``None``.  On success it
    is a dictionary containing:

    - ``'version'`` -- the certificate format version, currently ``1``;

    - ``'edge_mapping'`` -- a bijection between positions in
      ``list(G.edge_iterator())`` and ``list(H.edge_iterator())``;

    - ``'operations'`` -- a replayable sequence of vertex cleavings, Whitney
      twists, vertex identifications, and isolated-vertex adjustments; and

    - ``'vertex_mapping'`` -- the resulting graph isomorphism.

    The operations act on a normalized mutable copy of ``G``: its initial
    vertices and edges are numbered by their positions in ``list(G)`` and
    ``list(G.edge_iterator())``, respectively.  Fresh vertices introduced by
    cleaving have larger integer identifiers.  The vertex mapping takes the
    graph obtained after all operations to the similarly normalized ``H``.
    Cleaving steps contain ``vertex``, ``new_vertex``, and ``edges``; twist
    steps contain the two ``vertices`` and one ``edges`` side; identification
    steps contain ``keep`` and ``drop``.  Edge identifiers remain unchanged
    throughout the sequence.

    ALGORITHM:

    Loops and coloops are separated first.  The remaining cycle matroid is a
    direct sum over the nontrivial blocks.  Each block is decomposed into its
    SPQR tree and converted into a vertex-colored simple graph.  Series and
    parallel components are represented as unordered sets of ports, while the
    incidence graph of every rigid component is retained.  Canonical labeling
    of this auxiliary graph decides 2-isomorphism and supplies the edge map.
    This is a port-aware implementation of the SPQR-tree approach described
    in [RS2008]_.  The auxiliary graphs have linear size.  Runtime also
    includes block and SPQR decomposition and colored graph canonization of
    the rigid components.

    For a certificate, the matched block and SPQR decompositions are compiled
    into vertex cleavings, interval-reversal and virtual-port Whitney twists,
    and vertex identifications.  The resulting deterministic sequence is not
    intended to contain the minimum possible number of operations.  Since
    every twist explicitly lists one edge side, a deeply nested decomposition
    can produce a witness of quadratic total serialized size.

    EXAMPLES:

    A Whitney twist can change the graph isomorphism type without changing its
    cycle matroid::

        sage: G = Graph([('u','a'), ('a','v'), ('u','b'), ('b','a'),
        ....:            ('u','c'), ('c','v'), ('u','d'), ('d','c')])
        sage: H = Graph([('u','a'), ('a','v'), ('u','b'), ('b','a'),
        ....:            ('v','c'), ('c','u'), ('v','d'), ('d','c')])
        sage: G.is_isomorphic(H)
        False
        sage: G.is_2isomorphic(H)
        True
        sage: ok, witness = G.is_2isomorphic(H, certificate=True)
        sage: any(step['operation'] == 'whitney_twist'
        ....:     for step in witness['operations'])
        True
        sage: ok and G.verify_2isomorphism_certificate(H, witness)
        True

    Vertex identifications and cleavings between blocks are also ignored::

        sage: G = graphs.CycleGraph(3).disjoint_union(graphs.CycleGraph(4))
        sage: H = Graph([(0,1), (1,2), (2,0),
        ....:            (0,3), (3,4), (4,5), (5,0)])
        sage: G.is_2isomorphic(H)
        True
        sage: G = graphs.PathGraph(5)
        sage: H = graphs.StarGraph(4)
        sage: ok, witness = G.is_2isomorphic(H, certificate=True)
        sage: {'vertex_cleaving', 'vertex_identification'} <= {
        ....:     step['operation'] for step in witness['operations']}
        True
        sage: ok and G.verify_2isomorphism_certificate(H, witness)
        True

    Rigid components and loops are still distinguished::

        sage: graphs.CompleteBipartiteGraph(3, 3).is_2isomorphic(
        ....:     graphs.CircularLadderGraph(3))
        False
        sage: G = Graph([(0,0), (0,1)], loops=True)
        sage: H = Graph([(0,1), (0,1)], multiedges=True)
        sage: G.is_2isomorphic(H)
        False

    Certificates contain explicit Whitney operations and refer to edge
    occurrences, so repeated edge triples are safe::

        sage: G = Graph(multiedges=True)
        sage: G.add_edges([(0, 1, None), (0, 1, None), (1, 2, 'newVEdge0')])
        sage: H = Graph(multiedges=True)
        sage: H.add_edges([(3, 4, None), (3, 4, None), (4, 5, None)])
        sage: ok, witness = G.is_2isomorphic(H, certificate=True)
        sage: ok
        True
        sage: mapping = witness['edge_mapping']
        sage: set(mapping) == set(range(G.size()))
        True
        sage: set(mapping.values()) == set(range(H.size()))
        True
        sage: G.verify_2isomorphism_certificate(H, witness)
        True

    Isolated vertices do not affect the answer::

        sage: Graph(100).is_2isomorphic(Graph())
        True

    TESTS:

    Check empty and incompatible inputs::

        sage: ok, witness = Graph().is_2isomorphic(Graph(), certificate=True)
        sage: (ok, witness['edge_mapping'], witness['operations'])
        (True, {}, [])
        sage: Graph([(0, 1)]).is_2isomorphic(Graph(), certificate=True)
        (False, None)
        sage: Graph().is_2isomorphic(DiGraph())
        Traceback (most recent call last):
        ...
        TypeError: can only test 2-isomorphism between undirected graphs
        sage: from sage.graphs.morphisms import is_2isomorphic
        sage: is_2isomorphic(DiGraph(), Graph())
        Traceback (most recent call last):
        ...
        TypeError: can only test 2-isomorphism between undirected graphs

    Series ports may be permuted beyond the dihedral automorphisms of their
    skeleton cycle::

        sage: def decorated_cycle(lengths):
        ....:     graph = Graph(len(lengths))
        ....:     new_vertex = len(lengths)
        ....:     for i, length in enumerate(lengths):
        ....:         u, v = i, (i + 1) % len(lengths)
        ....:         graph.add_edge(u, v)
        ....:         path = [u]
        ....:         for _ in range(length - 1):
        ....:             path.append(new_vertex)
        ....:             new_vertex += 1
        ....:         graph.add_path(path + [v])
        ....:     return graph
        sage: G = decorated_cycle([2, 3, 4, 5])
        sage: H = decorated_cycle([2, 4, 3, 5])
        sage: G.is_isomorphic(H), G.is_2isomorphic(H)
        (False, True)

    In contrast, ports in different edge orbits of a rigid skeleton cannot be
    exchanged::

        sage: def decorated_wheel(rim_length, spoke_length):
        ....:     graph = Graph()
        ....:     graph.add_cycle([0, 1, 2, 3])
        ....:     graph.add_edges((4, i) for i in range(4))
        ....:     new_vertex = 5
        ....:     for (u, v), length in [((0, 1), rim_length),
        ....:                            ((4, 0), spoke_length)]:
        ....:         internal = range(new_vertex, new_vertex + length - 1)
        ....:         graph.add_path([u] + list(internal) + [v])
        ....:         new_vertex += length - 1
        ....:     return graph
        sage: G = decorated_wheel(2, 3)
        sage: H = decorated_wheel(3, 2)
        sage: G.is_2isomorphic(H)
        False

    User labels are not confused with the implementation's virtual-edge
    labels::

        sage: G.set_edge_label(0, 1, 'newVEdge0')
        sage: H = G.relabel({v: ('v', v) for v in G}, inplace=False)
        sage: G.is_2isomorphic(H)
        True

    Labels need not be hashable because the test works with private edge
    occurrence identifiers::

        sage: G = Graph()
        sage: G.add_edge(0, 1, ['unhashable'])
        sage: H = Graph()
        sage: H.add_edge('a', 'b', {'also': 'unhashable'})
        sage: G.is_2isomorphic(H)
        True
    """
    from sage.graphs.graph import Graph

    if not isinstance(G, Graph) or not isinstance(H, Graph):
        raise TypeError("can only test 2-isomorphism between undirected graphs")
    edge_mapping_only = certificate is _TWO_ISOMORPHISM_EDGE_MAPPING_ONLY
    if edge_mapping_only:
        certificate = True
    operation_certificate = certificate and not edge_mapping_only

    def certified(mapping, decompositions_G=(), decompositions_H=(),
                  auxiliary_mapping=None):
        if edge_mapping_only:
            return True, dict(mapping)
        witness = _two_isomorphism_operation_certificate(
            G, H, mapping, decompositions_G, decompositions_H,
            auxiliary_mapping)
        return True, witness

    if G.size() != H.size():
        return (False, None) if certificate else False
    if not G.size():
        return certified({}) if certificate else True
    if G is H:
        mapping = {i: i for i in range(G.size())}
        if not certificate:
            return True
        if edge_mapping_only:
            return True, mapping
        witness = {
            'version': 1,
            'edge_mapping': mapping,
            'operations': [],
            'vertex_mapping': {i: i for i in range(G.order())},
        }
        return True, witness

    edges_G, work_G, loops_G, labels_G = _two_isomorphism_work_graph(G)
    edges_H, work_H, loops_H, labels_H = _two_isomorphism_work_graph(H)
    if len(loops_G) != len(loops_H):
        return (False, None) if certificate else False
    rank_G = work_G.order() - work_G.connected_components_number()
    rank_H = work_H.order() - work_H.connected_components_number()
    if rank_G != rank_H:
        return (False, None) if certificate else False

    # Forests have only coloops once the loops have been removed, so any
    # occurrence bijection that preserves these two classes is a certificate.
    if work_G.size() == rank_G:
        if not certificate:
            return True
        loop_set_G = set(loops_G)
        loop_set_H = set(loops_H)
        nonloops_G = [i for i in range(len(edges_G)) if i not in loop_set_G]
        nonloops_H = [i for i in range(len(edges_H)) if i not in loop_set_H]
        mapping = dict(zip(loops_G, loops_H))
        mapping.update(zip(nonloops_G, nonloops_H))
        return certified(mapping)

    # Whitney's theorem reduces 2-isomorphism to ordinary graph isomorphism
    # for 3-connected simple graphs.  Keep this important fast path instead of
    # replacing a graph by its larger incidence encoding.
    if (not work_G.has_multiple_edges() and
            not work_H.has_multiple_edges() and
            work_G.is_triconnected() and work_H.is_triconnected()):
        if not certificate:
            return work_G.is_isomorphic(work_H)
        isomorphic, vertex_mapping = work_G.is_isomorphic(work_H,
                                                          certificate=True)
        if not isomorphic:
            return False, None
        mapping = dict(zip(loops_G, loops_H))
        target_edges = {(min(u, v), max(u, v)): label
                        for u, v, label in work_H.edge_iterator()}
        for u, v, label in work_G.edge_iterator():
            target_u = vertex_mapping[u]
            target_v = vertex_mapping[v]
            target_label = target_edges[(min(target_u, target_v),
                                         max(target_u, target_v))]
            mapping[labels_G[label]] = labels_H[target_label]
        return certified(mapping)

    auxiliary_G, colors_G, decompositions_G = _two_isomorphism_encoding(
        work_G, loops_G, labels_G, len(edges_G),
        decomposition=operation_certificate)
    auxiliary_H, colors_H, decompositions_H = _two_isomorphism_encoding(
        work_H, loops_H, labels_H, len(edges_H),
        decomposition=operation_certificate)
    color_sizes_G = tuple(len(colors_G[color])
                          for color in _TWO_ISOMORPHISM_COLORS)
    color_sizes_H = tuple(len(colors_H[color])
                          for color in _TWO_ISOMORPHISM_COLORS)
    if color_sizes_G != color_sizes_H:
        return (False, None) if certificate else False
    try:
        from sage.graphs.bliss import canonical_form  # noqa: F401
        algorithm = 'bliss'
    except ImportError:
        algorithm = 'sage'

    forms_G = _two_isomorphism_component_forms(
        auxiliary_G, colors_G, algorithm, certificate)
    forms_H = _two_isomorphism_component_forms(
        auxiliary_H, colors_H, algorithm, certificate)
    if (forms_G.keys() != forms_H.keys() or
            any(len(forms_G[key]) != len(forms_H[key]) for key in forms_G)):
        return (False, None) if certificate else False
    if not certificate:
        return True

    mapping = {}
    auxiliary_mapping = {}
    for color in ('loop edge', 'coloop edge'):
        source = [vertex[1] for vertex in colors_G[color]]
        target = [vertex[1] for vertex in colors_H[color]]
        mapping.update(zip(source, target))
    for key in forms_G:
        for (map_G, vertices_G, real_G), (map_H, vertices_H, _) in zip(
                forms_G[key], forms_H[key]):
            inverse_H = {image: vertex for vertex, image in map_H.items()}
            component_mapping = {
                vertex: inverse_H[map_G[vertex]] for vertex in vertices_G
            }
            if set(component_mapping.values()) != set(vertices_H):
                raise RuntimeError("canonical block mappings are inconsistent")
            auxiliary_mapping.update(component_mapping)
            mapping.update((vertex[1], component_mapping[vertex][1])
                           for vertex in real_G)
    return certified(mapping, decompositions_G, decompositions_H,
                     auxiliary_mapping)


def _two_isomorphism_edge_mapping(G, H):
    """
    Return only the edge-occurrence mapping used by ``GraphicMatroid``.

    The matroid API translates this mapping to its own groundset certificate
    and does not expose graph operations.  This private path avoids compiling
    a potentially large Whitney sequence that would immediately be discarded.
    """
    return is_2isomorphic(
        G, H, certificate=_TWO_ISOMORPHISM_EDGE_MAPPING_ONLY)


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

    from sage.numerical.mip import MixedIntegerLinearProgram, MIPSolverException
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
