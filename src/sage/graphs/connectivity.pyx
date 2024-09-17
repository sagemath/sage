# cython: binding=True
r"""
Connectivity related functions

This module implements the connectivity based functions for graphs and digraphs.
The methods in this module are also available as part of GenericGraph, DiGraph
or Graph classes as aliases, and these methods can be accessed through this
module or as class methods.
Here is what the module can do:

**For both directed and undirected graphs:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`is_connected` | Check whether the (di)graph is connected.
    :meth:`connected_components` | Return the list of connected components
    :meth:`connected_components_number` | Return the number of connected components.
    :meth:`connected_components_subgraphs` | Return a list of connected components as graph objects.
    :meth:`connected_component_containing_vertex` | Return a list of the vertices connected to vertex.
    :meth:`connected_components_sizes` | Return the sizes of the connected components as a list.
    :meth:`blocks_and_cut_vertices` | Return the blocks and cut vertices of the graph.
    :meth:`blocks_and_cuts_tree` | Return the blocks-and-cuts tree of the graph.
    :meth:`is_cut_edge` | Check whether the input edge is a cut-edge or a bridge.
    :meth:`is_edge_cut` | Check whether the input edges form an edge cut.
    :meth:`is_cut_vertex` | Check whether the input vertex is a cut-vertex.
    :meth:`is_vertex_cut` | Check whether the input vertices form a vertex cut.
    :meth:`edge_connectivity` | Return the edge connectivity of the graph.
    :meth:`vertex_connectivity` | Return the vertex connectivity of the graph.

**For DiGraph:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`is_strongly_connected` | Check whether the current ``DiGraph`` is strongly connected.
    :meth:`strongly_connected_components_digraph` | Return the digraph of the strongly connected components
    :meth:`strongly_connected_components_subgraphs` | Return the strongly connected components as a list of subgraphs.
    :meth:`strongly_connected_component_containing_vertex` | Return the strongly connected component containing a given vertex.
    :meth:`strong_articulation_points` | Return the strong articulation points of this digraph.

**For undirected graphs:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`bridges` | Return an iterator over the bridges (or cut edges) of given undirected graph.
    :meth:`cleave` | Return the connected subgraphs separated by the input vertex cut.
    :meth:`is_triconnected` | Check whether the graph is triconnected.
    :meth:`spqr_tree` | Return a SPQR-tree representing the triconnected components of the graph.
    :meth:`spqr_tree_to_graph` | Return the graph represented by the SPQR-tree `T`.

Methods
-------
"""

# ****************************************************************************
#
#       Copyright (C) 2023      David Coudert <david.coudert@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.superseded import deprecation
from sage.sets.disjoint_set cimport DisjointSet


def is_connected(G):
    """
    Check whether the (di)graph is connected.

    Note that in a graph, path connected is equivalent to connected.

    INPUT:

    - ``G`` -- the input graph

    .. SEEALSO::

        - :meth:`~Graph.is_biconnected`

    EXAMPLES::

        sage: from sage.graphs.connectivity import is_connected
        sage: G = Graph({0: [1, 2], 1: [2], 3: [4, 5], 4: [5]})
        sage: is_connected(G)
        False
        sage: G.is_connected()
        False
        sage: G.add_edge(0,3)
        sage: is_connected(G)
        True
        sage: D = DiGraph({0: [1, 2], 1: [2], 3: [4, 5], 4: [5]})
        sage: is_connected(D)
        False
        sage: D.add_edge(0, 3)
        sage: is_connected(D)
        True
        sage: D = DiGraph({1: [0], 2: [0]})
        sage: is_connected(D)
        True

    TESTS:

    If ``G`` is not a Sage graph, an error is raised::

        sage: from sage.graphs.connectivity import is_connected
        sage: is_connected('I am not a graph')
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage graph
    """
    from sage.graphs.generic_graph import GenericGraph
    if not isinstance(G, GenericGraph):
        raise TypeError("the input must be a Sage graph")

    if not G.order():
        return True

    try:
        return G._backend.is_connected()
    except AttributeError:
        v = next(G.vertex_iterator())
        conn_verts = list(G.depth_first_search(v, ignore_direction=True))
        return len(conn_verts) == G.num_verts()


def connected_components(G, sort=None, key=None):
    """
    Return the list of connected components.

    This returns a list of lists of vertices, each list representing a connected
    component. The list is ordered from largest to smallest component.

    INPUT:

    - ``G`` -- the input graph

    - ``sort`` -- boolean (default: ``None``); if ``True``, vertices inside each
      component are sorted according to the default ordering

      As of :issue:`35889`, this argument must be explicitly specified (unless a
      ``key`` is given); otherwise a warning is printed and ``sort=True`` is
      used. The default will eventually be changed to ``False``.

    - ``key`` -- a function (default: ``None``); a function that takes a
      vertex as its one argument and returns a value that can be used for
      comparisons in the sorting algorithm (we must have ``sort=True``)

    EXAMPLES::

        sage: from sage.graphs.connectivity import connected_components
        sage: G = Graph({0: [1, 3], 1: [2], 2: [3], 4: [5, 6], 5: [6]})
        sage: connected_components(G, sort=True)
        [[0, 1, 2, 3], [4, 5, 6]]
        sage: G.connected_components(sort=True)
        [[0, 1, 2, 3], [4, 5, 6]]
        sage: D = DiGraph({0: [1, 3], 1: [2], 2: [3], 4: [5, 6], 5: [6]})
        sage: connected_components(D, sort=True)
        [[0, 1, 2, 3], [4, 5, 6]]
        sage: connected_components(D, sort=True, key=lambda x: -x)
        [[3, 2, 1, 0], [6, 5, 4]]

    TESTS:

    If ``G`` is not a Sage graph, an error is raised::

        sage: from sage.graphs.connectivity import connected_components
        sage: connected_components('I am not a graph')
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage graph

    When parameter ``key`` is set, parameter ``sort`` must be ``True``::

        sage: G = Graph(2)
        sage: G.connected_components(sort=False, key=lambda x: x)
        Traceback (most recent call last):
        ...
        ValueError: sort keyword is False, yet a key function is given

    Deprecation warning for ``sort=None`` (:issue:`35889`)::

        sage: G = graphs.HouseGraph()
        sage: G.connected_components()
        doctest:...: DeprecationWarning: parameter 'sort' will be set to False by default in the future
        See https://github.com/sagemath/sage/issues/35889 for details.
        [[0, 1, 2, 3, 4]]
    """
    from sage.graphs.generic_graph import GenericGraph
    if not isinstance(G, GenericGraph):
        raise TypeError("the input must be a Sage graph")

    if sort is None:
        if key is None:
            deprecation(35889, "parameter 'sort' will be set to False by default in the future")
        sort = True

    if (not sort) and key:
        raise ValueError('sort keyword is False, yet a key function is given')

    cdef set seen = set()
    cdef list components = []
    for v in G:
        if v not in seen:
            c = connected_component_containing_vertex(G, v, sort=sort, key=key)
            seen.update(c)
            components.append(c)
    components.sort(key=lambda comp: -len(comp))
    return components


def connected_components_number(G):
    """
    Return the number of connected components.

    INPUT:

    - ``G`` -- the input graph

    EXAMPLES::

        sage: from sage.graphs.connectivity import connected_components_number
        sage: G = Graph({0: [1, 3], 1: [2], 2: [3], 4: [5, 6], 5: [6]})
        sage: connected_components_number(G)
        2
        sage: G.connected_components_number()
        2
        sage: D = DiGraph({0: [1, 3], 1: [2], 2: [3], 4: [5, 6], 5: [6]})
        sage: connected_components_number(D)
        2

    TESTS:

    If ``G`` is not a Sage graph, an error is raised::

        sage: from sage.graphs.connectivity import connected_components_number
        sage: connected_components_number('I am not a graph')
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage graph
    """
    return len(connected_components(G, sort=False))


def connected_components_subgraphs(G):
    """
    Return a list of connected components as graph objects.

    EXAMPLES::

        sage: from sage.graphs.connectivity import connected_components_subgraphs
        sage: G = Graph({0: [1, 3], 1: [2], 2: [3], 4: [5, 6], 5: [6]})
        sage: L = connected_components_subgraphs(G)
        sage: graphs_list.show_graphs(L)                                                # needs sage.plot
        sage: D = DiGraph({0: [1, 3], 1: [2], 2: [3], 4: [5, 6], 5: [6]})
        sage: L = connected_components_subgraphs(D)
        sage: graphs_list.show_graphs(L)                                                # needs sage.plot
        sage: L = D.connected_components_subgraphs()
        sage: graphs_list.show_graphs(L)                                                # needs sage.plot

    TESTS:

    If ``G`` is not a Sage graph, an error is raised::

        sage: from sage.graphs.connectivity import connected_components_subgraphs
        sage: connected_components_subgraphs('I am not a graph')
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage graph
    """
    from sage.graphs.generic_graph import GenericGraph
    if not isinstance(G, GenericGraph):
        raise TypeError("the input must be a Sage graph")

    return [G.subgraph(c, inplace=False) for c in connected_components(G, sort=False)]


def connected_component_containing_vertex(G, vertex, sort=None, key=None):
    """
    Return a list of the vertices connected to vertex.

    INPUT:

    - ``G`` -- the input graph

    - ``v`` -- the vertex to search for

    - ``sort`` -- boolean (default: ``None``); if ``True``, vertices inside the
      component are sorted according to the default ordering

      As of :issue:`35889`, this argument must be explicitly specified (unless a
      ``key`` is given); otherwise a warning is printed and ``sort=True`` is
      used. The default will eventually be changed to ``False``.

    - ``key`` -- a function (default: ``None``); a function that takes a
      vertex as its one argument and returns a value that can be used for
      comparisons in the sorting algorithm (we must have ``sort=True``)

    EXAMPLES::

        sage: from sage.graphs.connectivity import connected_component_containing_vertex
        sage: G = Graph({0: [1, 3], 1: [2], 2: [3], 4: [5, 6], 5: [6]})
        sage: connected_component_containing_vertex(G, 0, sort=True)
        [0, 1, 2, 3]
        sage: G.connected_component_containing_vertex(0, sort=True)
        [0, 1, 2, 3]
        sage: D = DiGraph({0: [1, 3], 1: [2], 2: [3], 4: [5, 6], 5: [6]})
        sage: connected_component_containing_vertex(D, 0, sort=True)
        [0, 1, 2, 3]
        sage: connected_component_containing_vertex(D, 0, sort=True, key=lambda x: -x)
        [3, 2, 1, 0]

    TESTS:

    If ``G`` is not a Sage graph, an error is raised::

        sage: from sage.graphs.connectivity import connected_component_containing_vertex
        sage: connected_component_containing_vertex('I am not a graph', 0)
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage graph

    :issue:`35889` is fixed::

        sage: G = Graph([('A', 1)])
        sage: G.connected_component_containing_vertex(1, sort=False)
        [1, 'A']
        sage: G.connected_component_containing_vertex(1, sort=True)
        Traceback (most recent call last):
        ...
        TypeError: '<' not supported between instances of 'str' and 'int'

    When parameter ``key`` is set, parameter ``sort`` must be ``True``::

        sage: G = Graph(2)
        sage: G.connected_component_containing_vertex(1, sort=False, key=lambda x: x)
        Traceback (most recent call last):
        ...
        ValueError: sort keyword is False, yet a key function is given

    Deprecation warning for ``sort=None`` (:issue:`35889`)::

        sage: G = graphs.HouseGraph()
        sage: G.connected_component_containing_vertex(1)
        doctest:...: DeprecationWarning: parameter 'sort' will be set to False by default in the future
        See https://github.com/sagemath/sage/issues/35889 for details.
        [0, 1, 2, 3, 4]
    """
    from sage.graphs.generic_graph import GenericGraph
    if not isinstance(G, GenericGraph):
        raise TypeError("the input must be a Sage graph")

    if sort is None:
        if key is None:
            deprecation(35889, "parameter 'sort' will be set to False by default in the future")
        sort = True

    if (not sort) and key:
        raise ValueError('sort keyword is False, yet a key function is given')

    try:
        c = list(G._backend.depth_first_search(vertex, ignore_direction=True))
    except AttributeError:
        c = list(G.depth_first_search(vertex, ignore_direction=True))

    if sort:
        return sorted(c, key=key)
    return c


def connected_components_sizes(G):
    """
    Return the sizes of the connected components as a list.

    The list is sorted from largest to lower values.

    EXAMPLES::

        sage: from sage.graphs.connectivity import connected_components_sizes
        sage: for x in graphs(3):
        ....:     print(connected_components_sizes(x))
        [1, 1, 1]
        [2, 1]
        [3]
        [3]
        sage: for x in graphs(3):
        ....:     print(x.connected_components_sizes())
        [1, 1, 1]
        [2, 1]
        [3]
        [3]

    TESTS:

    If ``G`` is not a Sage graph, an error is raised::

        sage: from sage.graphs.connectivity import connected_components_sizes
        sage: connected_components_sizes('I am not a graph')
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage graph
    """
    from sage.graphs.generic_graph import GenericGraph
    if not isinstance(G, GenericGraph):
        raise TypeError("the input must be a Sage graph")

    # connected components are sorted from largest to smallest
    return [len(cc) for cc in connected_components(G, sort=False)]


def blocks_and_cut_vertices(G, algorithm='Tarjan_Boost', sort=False, key=None):
    """
    Return the blocks and cut vertices of the graph.

    In the case of a digraph, this computation is done on the underlying
    graph.

    A cut vertex is one whose deletion increases the number of connected
    components. A block is a maximal induced subgraph which itself has no
    cut vertices. Two distinct blocks cannot overlap in more than a single
    cut vertex.

    INPUT:

    - ``algorithm`` -- string (default: ``'Tarjan_Boost'``); the algorithm to
      use among:

      - ``'Tarjan_Boost'`` -- default; Tarjan's algorithm (Boost
        implementation)

      - ``'Tarjan_Sage'`` -- Tarjan's algorithm (Sage implementation)

    - ``sort`` -- boolean (default: ``False``); whether to sort vertices inside
      the components and the list of cut vertices
      **currently only available for ``'Tarjan_Sage'``**

    - ``key`` -- a function (default: ``None``); a function that takes a
      vertex as its one argument and returns a value that can be used for
      comparisons in the sorting algorithm (we must have ``sort=True``)

    OUTPUT: ``(B, C)``, where ``B`` is a list of blocks - each is a list of
    vertices and the blocks are the corresponding induced subgraphs - and
    ``C`` is a list of cut vertices.

    ALGORITHM:

      We implement the algorithm proposed by Tarjan in [Tarjan72]_. The
      original version is recursive. We emulate the recursion using a stack.

    .. SEEALSO::

        - :meth:`blocks_and_cuts_tree`
        - :func:`sage.graphs.base.boost_graph.blocks_and_cut_vertices`
        - :meth:`~Graph.is_biconnected`
        - :meth:`~Graph.bridges`

    EXAMPLES:

    We construct a trivial example of a graph with one cut vertex::

        sage: from sage.graphs.connectivity import blocks_and_cut_vertices
        sage: rings = graphs.CycleGraph(10)
        sage: rings.merge_vertices([0, 5])
        sage: blocks_and_cut_vertices(rings)
        ([[0, 1, 4, 2, 3], [0, 6, 9, 7, 8]], [0])
        sage: rings.blocks_and_cut_vertices()
        ([[0, 1, 4, 2, 3], [0, 6, 9, 7, 8]], [0])
        sage: B, C = blocks_and_cut_vertices(rings, algorithm='Tarjan_Sage', sort=True)
        sage: B, C
        ([[0, 1, 2, 3, 4], [0, 6, 7, 8, 9]], [0])
        sage: B2, C2 = blocks_and_cut_vertices(rings, algorithm='Tarjan_Sage', sort=False)
        sage: Set(map(Set, B)) == Set(map(Set, B2)) and set(C) == set(C2)
        True

    The Petersen graph is biconnected, hence has no cut vertices::

        sage: blocks_and_cut_vertices(graphs.PetersenGraph())
        ([[0, 1, 4, 5, 2, 6, 3, 7, 8, 9]], [])

    Decomposing paths to pairs::

        sage: g = graphs.PathGraph(4) + graphs.PathGraph(5)
        sage: blocks_and_cut_vertices(g)
        ([[2, 3], [1, 2], [0, 1], [7, 8], [6, 7], [5, 6], [4, 5]], [1, 2, 5, 6, 7])

    A disconnected graph::

        sage: g = Graph({1: {2: 28, 3: 10}, 2: {1: 10, 3: 16}, 4: {}, 5: {6: 3, 7: 10, 8: 4}})
        sage: blocks_and_cut_vertices(g)
        ([[1, 2, 3], [5, 6], [5, 7], [5, 8], [4]], [5])

    A directed graph with Boost's algorithm (:issue:`25994`)::

        sage: rings = graphs.CycleGraph(10)
        sage: rings.merge_vertices([0, 5])
        sage: rings = rings.to_directed()
        sage: blocks_and_cut_vertices(rings, algorithm='Tarjan_Boost')
        ([[0, 1, 4, 2, 3], [0, 6, 9, 7, 8]], [0])

    TESTS::

        sage: blocks_and_cut_vertices(Graph(0))
        ([], [])
        sage: blocks_and_cut_vertices(Graph(1))
        ([[0]], [])
        sage: blocks_and_cut_vertices(Graph(2))
        ([[0], [1]], [])

    If ``G`` is not a Sage graph, an error is raised::

        sage: from sage.graphs.connectivity import connected_components_sizes
        sage: connected_components_sizes('I am not a graph')
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage graph
    """
    from sage.graphs.generic_graph import GenericGraph
    if not isinstance(G, GenericGraph):
        raise TypeError("the input must be a Sage graph")

    if algorithm == "Tarjan_Boost":
        from sage.graphs.base.boost_graph import blocks_and_cut_vertices
        return blocks_and_cut_vertices(G)

    if algorithm != "Tarjan_Sage":
        raise NotImplementedError("blocks and cut vertices algorithm '%s' is not implemented" % algorithm)

    # If algorithm is "Tarjan_Sage"
    if (not sort) and key:
        raise ValueError('sort keyword is False, yet a key function is given')

    blocks = []
    cut_vertices = set()

    # We iterate over all vertices to ensure that we visit each connected
    # component of the graph
    seen = set()
    for start in G.vertex_iterator():
        if start in seen:
            continue

        # Special case of an isolated vertex
        if not G.degree(start):
            blocks.append([start])
            seen.add(start)
            continue

        # Each vertex is numbered with an integer from 1...|V(G)|,
        # corresponding to the order in which it is discovered during the
        # DFS.
        number = {}
        num = 1

        # Associates to each vertex v the smallest number of a vertex that
        # can be reached from v in the orientation of the graph that the
        # algorithm creates.
        low_point = {}

        # Associates to each vertex an iterator over its neighbors
        neighbors = {}

        stack = [start]
        edge_stack = []
        start_already_seen = False

        while stack:
            v = stack[-1]
            seen.add(v)

            # The first time we meet v
            if v not in number:
                # We number the vertices in the order they are reached
                # during DFS
                number[v] = num
                neighbors[v] = G.neighbor_iterator(v)
                low_point[v] = num
                num += 1

            try:
                # We consider the next of its neighbors
                w = next(neighbors[v])

                # If we never met w before, we remember the direction of
                # edge vw, and add w to the stack.
                if w not in number:
                    edge_stack.append((v, w))
                    stack.append(w)

                # If w is an ancestor of v in the DFS tree, we remember the
                # direction of edge vw
                elif number[w] < number[v]:
                    edge_stack.append((v, w))
                    low_point[v] = min(low_point[v], number[w])

            # We went through all of v's neighbors
            except StopIteration:
                # We trackback, so w takes the value of v and we pop the
                # stack
                w = stack.pop()

                # Test termination of the algorithm
                if not stack:
                    break

                v = stack[-1]

                # Propagating the information : low_point[v] indicates the
                # smallest vertex (the vertex x with smallest number[x])
                # that can be reached from v
                low_point[v] = min(low_point[v], low_point[w])

                # The situation in which there is no path from w to an
                # ancestor of v : we have identified a new biconnected
                # component
                if low_point[w] >= number[v]:
                    new_block = set()
                    nw = number[w]
                    u1, u2 = edge_stack.pop()
                    while number[u1] >= nw:
                        new_block.add(u1)
                        u1, u2 = edge_stack.pop()
                    new_block.add(u1)
                    if sort:
                        this_block = sorted(new_block, key=key)
                    else:
                        this_block = list(new_block)
                    blocks.append(this_block)

                    # We update the set of cut vertices.
                    #
                    # If v is start, then we add it only if it belongs to
                    # several blocks.
                    if (v is not start) or start_already_seen:
                        cut_vertices.add(v)
                    else:
                        start_already_seen = True

    if sort:
        return blocks, sorted(cut_vertices, key=key)
    return blocks, list(cut_vertices)


def blocks_and_cuts_tree(G):
    """
    Return the blocks-and-cuts tree of ``self``.

    This new graph has two different kinds of vertices, some representing the
    blocks (type B) and some other the cut vertices of the graph (type C).

    There is an edge between a vertex `u` of type B and a vertex `v` of type C
    if the cut-vertex corresponding to `v` is in the block corresponding to `u`.

    The resulting graph is a tree, with the additional characteristic property
    that the distance between two leaves is even. When ``self`` is not
    connected, the resulting graph is a forest.

    When ``self`` is biconnected, the tree is reduced to a single node of
    type `B`.

    We referred to [HarPri]_ and [Gallai]_ for blocks and cuts tree.

    .. SEEALSO::

        - :meth:`~sage.graphs.generic_graph.GenericGraph.blocks_and_cut_vertices`
        - :meth:`~Graph.is_biconnected`

    EXAMPLES::

        sage: from sage.graphs.connectivity import blocks_and_cuts_tree
        sage: T = blocks_and_cuts_tree(graphs.KrackhardtKiteGraph()); T
        Graph on 5 vertices
        sage: T.is_isomorphic(graphs.PathGraph(5))
        True
        sage: from sage.graphs.connectivity import blocks_and_cuts_tree
        sage: T = graphs.KrackhardtKiteGraph().blocks_and_cuts_tree(); T
        Graph on 5 vertices

    The distance between two leaves is even::

        sage: T = blocks_and_cuts_tree(graphs.RandomTree(40))
        sage: T.is_tree()
        True
        sage: leaves = [v for v in T if T.degree(v) == 1]
        sage: all(T.distance(u,v) % 2 == 0 for u in leaves for v in leaves)
        True

    The tree of a biconnected graph has a single vertex, of type `B`::

        sage: T = blocks_and_cuts_tree(graphs.PetersenGraph())
        sage: T.vertices(sort=True)
        [('B', (0, 1, 4, 5, 2, 6, 3, 7, 8, 9))]

    TESTS:

    When ``self`` is not connected, the resulting graph is a forest
    (:issue:`24163`)::

        sage: from sage.graphs.connectivity import blocks_and_cuts_tree
        sage: T = blocks_and_cuts_tree(Graph(2))
        sage: T.is_forest()
        True

    If ``G`` is not a Sage graph, an error is raised::

        sage: blocks_and_cuts_tree('I am not a graph')
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage graph
    """
    from sage.graphs.generic_graph import GenericGraph
    if not isinstance(G, GenericGraph):
        raise TypeError("the input must be a Sage graph")

    from sage.graphs.graph import Graph
    B, C = G.blocks_and_cut_vertices()
    B = map(tuple, B)
    set_C = set(C)
    g = Graph()
    for bloc in B:
        g.add_vertex(('B', bloc))
        for c in bloc:
            if c in set_C:
                g.add_edge(('B', bloc), ('C', c))
    return g


def is_edge_cut(G, edges):
    """
    Check whether ``edges`` form an edge cut.

    A set of edges is an edge cut of a graph if its removal increases the number
    of connected components. In a digraph, we consider the number of (weakly)
    connected components.

    This method is not working for (di)graphs with multiple edges. Furthermore,
    edge labels are ignored.

    INPUT:

    - ``G`` -- a (di)graph

    - ``edges`` -- a set of edges

    EXAMPLES:

    A cycle graph of order 4::

        sage: from sage.graphs.connectivity import is_edge_cut
        sage: G = graphs.CycleGraph(4)
        sage: is_edge_cut(G, [(1, 2)])
        False
        sage: is_edge_cut(G, [(1, 2), (2, 3)])
        True
        sage: is_edge_cut(G, [(1, 2), (3, 0)])
        True

    A pending edge is a cut-edge::

        sage: G.add_edge((0, 5, 'silly'))
        sage: is_edge_cut(G, [(0, 5, 'silly')])
        True

    Edge labels are ignored, even if specified::

        sage: G.add_edge((2, 5, 'xyz'))
        sage: is_edge_cut(G, [(0, 5), (2, 5)])
        True
        sage: is_edge_cut(G, [(0, 5), (2, 5, 'xyz')])
        True
        sage: is_edge_cut(G, [(0, 5, 'silly'), (2, 5)])
        True
        sage: is_edge_cut(G, [(0, 5, 'aa'), (2, 5, 'bb')])
        True

    The graph can have loops::

        sage: G.allow_loops(True)
        sage: G.add_edge(0, 0)
        sage: is_edge_cut(G, [(0, 5), (2, 5)])
        True
        sage: is_edge_cut(G, [(0, 0), (0, 5), (2, 5)])
        True

    Multiple edges are not allowed::

        sage: G.allow_multiple_edges(True)
        sage: is_edge_cut(G, [(0, 5), (2, 5)])
        Traceback (most recent call last):
        ...
        ValueError: This method is not known to work on graphs with
         multiedges. Perhaps this method can be updated to handle them, but in
         the meantime if you want to use it please disallow multiedges using
         allow_multiple_edges().

    An error is raised if an element of ``edges`` is not an edge of `G`::

        sage: G = graphs.CycleGraph(4)
        sage: is_edge_cut(G, [(0, 2)])
        Traceback (most recent call last):
        ...
        ValueError: edge (0, 2) is not an edge of the graph

    For digraphs, this method considers the number of (weakly) connected
    components::

        sage: G = digraphs.Circuit(4)
        sage: is_edge_cut(G, [(0, 1)])
        False
        sage: G = digraphs.Circuit(4)
        sage: is_edge_cut(G, [(0, 1), (1, 2)])
        True

    For disconnected (di)graphs, the method checks if the number of (weakly)
    connected components increases::

        sage: G = graphs.CycleGraph(4) * 2
        sage: is_edge_cut(G, [(1, 2), (2, 3)])
        True
        sage: G = digraphs.Circuit(4) * 2
        sage: is_edge_cut(G, [(0, 1), (1, 2)])
        True
    """
    G._scream_if_not_simple(allow_loops=True)

    cdef set C = set()  # set of edges of the potential cut
    cdef set S = set()  # set of incident vertices
    for e in edges:
        u, v = e[0], e[1]
        if not G.has_edge(u, v):
            raise ValueError("edge {0} is not an edge of the graph".format(repr(e)))
        if u == v:
            # We ignore loops
            continue
        if G.degree(u) == 1 or G.degree(v) == 1:
            # e is a pending edge and so a cut-edge
            return True
        S.add(u)
        S.add(v)
        C.add((u, v))
        if not G.is_directed():
            C.add((v, u))

    cdef list queue
    cdef set seen
    DS = DisjointSet(G)

    for comp in G.connected_components():
        if not S.intersection(comp):
            # This component is not involved in the cut
            continue

        # We run a DFS in comp from any vertex and avoid edges in C
        start = comp[0]
        queue = [start]
        seen = set(queue)
        while queue:
            v = queue.pop()
            for e in G.edge_iterator(vertices=[v], labels=False, ignore_direction=True, sort_vertices=False):
                if e in C:
                    continue
                w = e[1] if e[0] == v else e[0]
                if w not in seen:
                    seen.add(w)
                    DS.union(v, w)
                    queue.append(w)

        # We now check if some vertices of comp have not been reached
        if len(set(DS.find(v) for v in comp)) > 1:
            return True

    return False


def is_cut_edge(G, u, v=None, label=None):
    """
    Check whether the edge ``(u, v)`` is a cut-edge or a bridge of graph ``G``.

    A cut edge (or bridge) is an edge that when removed increases
    the number of connected components. This function works with
    simple graphs as well as graphs with loops and multiedges. In
    a digraph, a cut edge is an edge that when removed increases
    the number of (weakly) connected components.

    INPUT: The following forms are accepted

    - is_cut_edge(G, 1, 2 )

    - is_cut_edge(G, (1, 2) )

    - is_cut_edge(G, 1, 2, 'label' )

    - is_cut_edge(G, (1, 2, 'label') )

    OUTPUT:

    - Returns ``True`` if (u,v) is a cut edge, False otherwise

    EXAMPLES::

        sage: from sage.graphs.connectivity import is_cut_edge
        sage: G = graphs.CompleteGraph(4)
        sage: is_cut_edge(G,0,2)
        False
        sage: G.is_cut_edge(0,2)
        False

        sage: G = graphs.CompleteGraph(4)
        sage: G.add_edge((0,5,'silly'))
        sage: is_cut_edge(G,(0,5,'silly'))
        True

        sage: G = Graph([[0,1],[0,2],[3,4],[4,5],[3,5]])
        sage: is_cut_edge(G,(0,1))
        True

        sage: G = Graph([[0,1],[0,2],[1,1]], loops = True)
        sage: is_cut_edge(G,(1,1))
        False

        sage: G = digraphs.Circuit(5)
        sage: is_cut_edge(G,(0,1))
        False

        sage: G = graphs.CompleteGraph(6)
        sage: is_cut_edge(G,(0,7))
        Traceback (most recent call last):
        ...
        ValueError: edge not in graph
    """
    if label is None:
        if v is None:
            try:
                u, v, label = u
            except ValueError:
                u, v = u
                label = None

    if not G.has_edge(u, v):
        raise ValueError('edge not in graph')

    # If edge (u,v) is a pending edge, it is also a cut-edge
    if G.degree(u) == 1 or G.degree(v) == 1:
        return True
    elif G.allows_multiple_edges():
        # If we have two or more edges between u and v, it is not a cut-edge
        if len([(uu, vv) for uu, vv, ll in G.edges_incident(u) if uu == v or vv == v]) > 1:
            return False

    g = G.copy(immutable=False) if G.is_immutable() else G
    g.delete_edge(u, v, label)
    if g.is_directed():
        # (u,v) is a cut-edge if u is not in the connected
        # component containing v of self-(u,v)
        sol = u not in connected_component_containing_vertex(g, v)
    else:
        # (u,v) is a cut-edge if there is no path from u to v in
        # self-(u,v)
        sol = not g.distance(u, v) < g.order()

    g.add_edge(u, v, label)
    return sol


def is_vertex_cut(G, cut, weak=False):
    r"""
    Check whether the input vertices form a vertex cut.

    A set of vertices is a vertex cut if its removal from the (di)graph
    increases the number of (strongly) connected components. This function works
    with simple graphs as well as graphs with loops and multiple edges.

    INPUT:

    - ``G`` -- a Sage (Di)Graph

    - ``cut`` -- a set of vertices

    - ``weak`` -- boolean (default: ``False``); whether the connectivity of
      directed graphs is to be taken in the weak sense, that is ignoring edges
      orientations

    EXAMPLES:

    Giving a cycle graph of order 4::

        sage: from sage.graphs.connectivity import is_vertex_cut
        sage: G = graphs.CycleGraph(4)
        sage: is_vertex_cut(G, [0, 1])
        False
        sage: is_vertex_cut(G, [0, 2])
        True

    Giving a disconnected graph::

        sage: from sage.graphs.connectivity import is_vertex_cut
        sage: G = graphs.CycleGraph(4) * 2
        sage: G.connected_components()
        [[0, 1, 2, 3], [4, 5, 6, 7]]
        sage: is_vertex_cut(G, [0, 2])
        True
        sage: is_vertex_cut(G, [4, 6])
        True
        sage: is_vertex_cut(G, [0, 6])
        False
        sage: is_vertex_cut(G, [0, 4, 6])
        True

    Comparing the weak and strong connectivity of a digraph::

        sage: D = digraphs.Circuit(6)
        sage: D.is_strongly_connected()
        True
        sage: is_vertex_cut(D, [2])
        True
        sage: is_vertex_cut(D, [2], weak=True)
        False

    Giving a vertex that is not in the graph::

        sage: G = graphs.CompleteGraph(4)
        sage: is_vertex_cut(G, [7])
        Traceback (most recent call last):
        ...
        ValueError: vertex (7) is not a vertex of the graph

    TESTS:

    If ``G`` is not a Sage graph, an error is raised::

        sage: is_vertex_cut('I am not a graph', [0])
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage graph
    """
    from sage.graphs.generic_graph import GenericGraph
    if not isinstance(G, GenericGraph):
        raise TypeError("the input must be a Sage graph")

    cdef set cutset = set(cut)
    for u in cutset:
        if u not in G:
            raise ValueError("vertex ({0}) is not a vertex of the graph".format(repr(u)))

    if len(cutset) >= G.order() - 1:
        # A vertex cut must be of size at most n - 2
        return False

    # We deal with graphs with multiple (strongly) connected components
    cdef list CC
    if G.is_directed() and not weak:
        CC = G.strongly_connected_components()
    else:
        CC = G.connected_components(sort=False)
    if len(CC) > 1:
        for comp in CC:
            subcut = cutset.intersection(comp)
            if subcut and is_vertex_cut(G.subgraph(comp), subcut, weak=weak):
                return True
        return False

    cdef list boundary = G.vertex_boundary(cutset)
    if not boundary:
        # We need at least 1 vertex in the boundary of the cut
        return False

    cdef list cases = [(G.neighbor_iterator, boundary)]
    if not weak and G.is_directed():
        # Strong connectivity for digraphs.
        # We perform two DFS starting from an out neighbor of cut and avoiding
        # cut. The first DFS follows the edges directions, and the second is
        # in the reverse order. If both allow to reach all neighbors of cut,
        # then it is not a vertex cut.
        # We set data for the reverse order
        in_boundary = set()
        for u in cutset:
            in_boundary.update(G.neighbor_in_iterator(u))
        in_boundary.difference_update(cutset)
        if not in_boundary:
            return False
        cases.append((G.neighbor_in_iterator, list(in_boundary)))

    cdef list queue
    cdef set seen
    cdef set targets
    start = boundary[0]

    for neighbors, this_boundary in cases:

        # We perform a DFS starting from start and avoiding cut
        queue = [start]
        seen = set(cutset)
        seen.add(start)
        targets = set(this_boundary)
        targets.discard(start)
        while queue:
            v = queue.pop()
            for w in neighbors(v):
                if w not in seen:
                    seen.add(w)
                    queue.append(w)
                    targets.discard(w)

        # If some neighbors cannot be reached, we have a vertex cut
        if targets:
            return True

    return False


def is_cut_vertex(G, u, weak=False):
    r"""
    Check whether the input vertex is a cut-vertex.

    A vertex is a cut-vertex if its removal from the (di)graph increases the
    number of (strongly) connected components. Isolated vertices or leaves are
    not cut-vertices. This function works with simple graphs as well as graphs
    with loops and multiple edges.

    INPUT:

    - ``G`` -- a Sage (Di)Graph

    - ``u`` -- a vertex

    - ``weak`` -- boolean (default: ``False``); whether the connectivity of
      directed graphs is to be taken in the weak sense, that is ignoring edges
      orientations

    OUTPUT:

    Return ``True`` if ``u`` is a cut-vertex, and ``False`` otherwise.

    EXAMPLES:

    Giving a LollipopGraph(4,2), that is a complete graph with 4 vertices with a
    pending edge::

        sage: from sage.graphs.connectivity import is_cut_vertex
        sage: G = graphs.LollipopGraph(4, 2)
        sage: is_cut_vertex(G, 0)
        False
        sage: is_cut_vertex(G, 3)
        True
        sage: G.is_cut_vertex(3)
        True

    Comparing the weak and strong connectivity of a digraph::

        sage: D = digraphs.Circuit(6)
        sage: D.is_strongly_connected()
        True
        sage: is_cut_vertex(D, 2)
        True
        sage: is_cut_vertex(D, 2, weak=True)
        False

    Giving a vertex that is not in the graph::

        sage: G = graphs.CompleteGraph(4)
        sage: is_cut_vertex(G, 7)
        Traceback (most recent call last):
        ...
        ValueError: vertex (7) is not a vertex of the graph

    TESTS:

    If ``G`` is not a Sage graph, an error is raised::

        sage: is_cut_vertex('I am not a graph', 0)
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage graph
    """
    return is_vertex_cut(G, [u], weak=weak)


def edge_connectivity(G,
                      value_only=True,
                      implementation=None,
                      use_edge_labels=False,
                      vertices=False,
                      solver=None,
                      verbose=0,
                      *, integrality_tolerance=1e-3):
    r"""
    Return the edge connectivity of the graph.

    For more information, see the :wikipedia:`Connectivity_(graph_theory)`.

    .. NOTE::

        When the graph is a directed graph, this method actually computes the
        *strong* connectivity, (i.e. a directed graph is strongly `k`-connected
        if there are `k` disjoint paths between any two vertices `u, v`). If you
        do not want to consider strong connectivity, the best is probably to
        convert your ``DiGraph`` object to a ``Graph`` object, and compute the
        connectivity of this other graph.

    INPUT:

    - ``G`` -- the input Sage (Di)Graph

    - ``value_only`` -- boolean (default: ``True``)

      - When set to ``True`` (default), only the value is returned.

      - When set to ``False``, both the value and a minimum vertex cut are
        returned.

    - ``implementation`` -- string (default: ``None``); selects an
      implementation:

      - ``None`` -- default; selects the best implementation available

      - ``'boost'`` -- use the Boost graph library (which is much more
        efficient). It is not available when ``edge_labels=True``, and it is
        unreliable for directed graphs (see :issue:`18753`).

      - ``'Sage'`` -- use Sage's implementation based on integer linear
        programming

    - ``use_edge_labels`` -- boolean (default: ``False``)

      - When set to ``True``, computes a weighted minimum cut where each edge
        has a weight defined by its label. (If an edge has no label, `1` is
        assumed.). Implies ``boost`` = ``False``.

      - When set to ``False``, each edge has weight `1`.

    - ``vertices`` -- boolean (default: ``False``)

      - When set to ``True``, also returns the two sets of vertices that are
        disconnected by the cut. Implies ``value_only=False``.

    - ``solver`` -- string (default: ``None``); specifies a Mixed Integer Linear
      Programming (MILP) solver to be used. If set to ``None``, the default one
      is used. For more information on MILP solvers and which default solver is
      used, see the method :meth:`solve
      <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
      :class:`MixedIntegerLinearProgram
      <sage.numerical.mip.MixedIntegerLinearProgram>`.

    - ``verbose`` -- integer (default: 0); sets the level of verbosity. Set
      to 0 by default, which means quiet.

    - ``integrality_tolerance`` -- float; parameter for use with MILP solvers
      over an inexact base ring; see
      :meth:`MixedIntegerLinearProgram.get_values`.

    EXAMPLES:

    A basic application on the PappusGraph::

        sage: from sage.graphs.connectivity import edge_connectivity
        sage: g = graphs.PappusGraph()
        sage: edge_connectivity(g)
        3
        sage: g.edge_connectivity()
        3

    The edge connectivity of a complete graph is its minimum degree, and one of
    the two parts of the bipartition is reduced to only one vertex. The graph of
    the cut edges is isomorphic to a Star graph::

        sage: g = graphs.CompleteGraph(5)
        sage: [ value, edges, [ setA, setB ]] = edge_connectivity(g,vertices=True)
        sage: value
        4
        sage: len(setA) == 1 or len(setB) == 1
        True
        sage: cut = Graph()
        sage: cut.add_edges(edges)
        sage: cut.is_isomorphic(graphs.StarGraph(4))
        True

    Even if obviously in any graph we know that the edge connectivity is less
    than the minimum degree of the graph::

        sage: g = graphs.RandomGNP(10,.3)
        sage: min(g.degree()) >= edge_connectivity(g)
        True

    If we build a tree then assign to its edges a random value, the minimum cut
    will be the edge with minimum value::

        sage: tree = graphs.RandomTree(10)
        sage: for u,v in tree.edge_iterator(labels=None):
        ....:      tree.set_edge_label(u, v, random())
        sage: minimum = min(tree.edge_labels())
        sage: [_, [(_, _, l)]] = edge_connectivity(tree, value_only=False,              # needs sage.numerical.mip
        ....:                                      use_edge_labels=True)
        sage: l == minimum                                                              # needs sage.numerical.mip
        True

    When ``value_only=True`` and ``implementation="sage"``, this function is
    optimized for small connectivity values and does not need to build a linear
    program.

    It is the case for graphs which are not connected ::

        sage: g = 2 * graphs.PetersenGraph()
        sage: edge_connectivity(g, implementation='sage')
        0.0

    For directed graphs, the strong connectivity is tested through the dedicated
    function::

        sage: g = digraphs.ButterflyGraph(3)
        sage: edge_connectivity(g, implementation='sage')
        0.0

    We check that the result with Boost is the same as the result without Boost::

        sage: g = graphs.RandomGNP(15, .3)
        sage: (edge_connectivity(g, implementation='boost')                             # needs sage.numerical.mip
        ....:    == edge_connectivity(g, implementation='sage'))
        True

    Boost interface also works with directed graphs::

        sage: edge_connectivity(digraphs.Circuit(10), implementation='boost',
        ....:                   vertices=True)
        [1, [(0, 1)], [{0}, {1, 2, 3, 4, 5, 6, 7, 8, 9}]]

    However, the Boost algorithm is not reliable if the input is directed
    (see :issue:`18753`)::

        sage: g = digraphs.Path(3)
        sage: edge_connectivity(g)
        0.0
        sage: edge_connectivity(g, implementation='boost')
        1
        sage: g.add_edge(1, 0)
        sage: edge_connectivity(g)
        0.0
        sage: edge_connectivity(g, implementation='boost')
        0

    TESTS:

    Checking that the two implementations agree::

        sage: for i in range(10):                                                       # needs sage.numerical.mip
        ....:     g = graphs.RandomGNP(30, 0.3)
        ....:     e1 = edge_connectivity(g, implementation='boost')
        ....:     e2 = edge_connectivity(g, implementation='sage')
        ....:     assert (e1 == e2)

    Disconnected graphs and ``vertices=True``::

        sage: g = graphs.PetersenGraph()
        sage: edge_connectivity((2 * g), vertices=True)
        [0, [], [{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, {10, 11, 12, 13, 14, 15, 16, 17, 18, 19}]]
        sage: edge_connectivity(Graph(), vertices=True)
        [0, [], [{}, {}]]

    If ``G`` is not a Sage graph, an error is raised::

        sage: edge_connectivity('I am not a graph')
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage graph
    """
    from sage.graphs.generic_graph import GenericGraph
    if not isinstance(G, GenericGraph):
        raise TypeError("the input must be a Sage graph")

    G._scream_if_not_simple(allow_loops=True)
    g = G

    if vertices:
        value_only = False

    if implementation is None:
        if use_edge_labels or g.is_directed():
            implementation = "sage"
        else:
            implementation = "boost"

    implementation = implementation.lower()
    if implementation not in ["boost", "sage"]:
        raise ValueError("'implementation' must be set to 'boost', 'sage' or None.")
    elif implementation == "boost" and use_edge_labels:
        raise ValueError("the Boost implementation is currently not able to handle edge labels")

    # Otherwise, an error is created
    if not g.num_edges() or not g.num_verts():
        if value_only:
            return 0
        elif vertices:
            return [0, [], [{}, {}]]
        return [0, []]

    if implementation == "boost":
        from sage.graphs.base.boost_graph import edge_connectivity

        [obj, edges] = edge_connectivity(g)

        if value_only:
            return obj

        val = [obj, edges]

        if vertices and obj:
            H = G.copy()
            H.delete_edges(edges)

            if H.is_directed():
                a = set(H.breadth_first_search([x for x, y in edges]))
                b = set(H).difference(a)
                val.append([a, b])
            else:
                val.append([set(c) for c in connected_components(H, sort=False)])
        elif vertices:
            val.append([set(c) for c in connected_components(G, sort=False)])

        return val

    if use_edge_labels:
        from sage.rings.real_mpfr import RR

        def weight(x):
            return x if x in RR else 1
    else:
        def weight(x):
            return 1

    # Better methods for small connectivity tests, when one is not interested in
    # cuts...
    if value_only and not use_edge_labels:

        if G.is_directed():
            if not is_strongly_connected(G):
                return 0.0

        else:
            if not is_connected(G):
                return 0.0

            h = G.strong_orientation()
            if not is_strongly_connected(h):
                return 1.0

    from sage.numerical.mip import MixedIntegerLinearProgram

    p = MixedIntegerLinearProgram(maximization=False, solver=solver)

    in_set = p.new_variable(binary=True)
    in_cut = p.new_variable(binary=True)

    # A vertex has to be in some set
    for v in g:
        p.add_constraint(in_set[0, v] + in_set[1, v], max=1, min=1)

    # There is no empty set
    p.add_constraint(p.sum(in_set[1, v] for v in g), min=1)
    p.add_constraint(p.sum(in_set[0, v] for v in g), min=1)

    if g.is_directed():
        # There is no edge from set 0 to set 1 which is not in the cut
        for u, v in g.edge_iterator(labels=None):
            p.add_constraint(in_set[0, u] + in_set[1, v] - in_cut[u, v], max=1)

        p.set_objective(p.sum(weight(l) * in_cut[u, v] for u, v, l in g.edge_iterator()))

    else:

        # Two adjacent vertices are in different sets if and only if
        # the edge between them is in the cut
        for u, v in g.edge_iterator(labels=None):
            p.add_constraint(in_set[0, u] + in_set[1, v] - in_cut[frozenset((u, v))], max=1)
            p.add_constraint(in_set[1, u] + in_set[0, v] - in_cut[frozenset((u, v))], max=1)

        p.set_objective(p.sum(weight(l) * in_cut[frozenset((u, v))] for u, v, l in g.edge_iterator()))

    obj = p.solve(log=verbose)

    in_cut = p.get_values(in_cut, convert=bool, tolerance=integrality_tolerance)

    if use_edge_labels is False:
        if g.is_directed():
            obj = sum(1 for u, v in g.edge_iterator(labels=False) if in_cut[u, v])
        else:
            obj = sum(1 for u, v in g.edge_iterator(labels=False) if in_cut[frozenset((u, v))])

    if value_only:
        return obj

    val = [obj]
    in_set = p.get_values(in_set, convert=bool, tolerance=integrality_tolerance)

    if g.is_directed():
        edges = [(u, v, l) for u, v, l in g.edge_iterator() if in_cut[u, v]]
    else:
        edges = [(u, v, l) for u, v, l in g.edge_iterator() if in_cut[frozenset((u, v))]]

    val.append(edges)

    if vertices:
        a = {}
        b = {}
        for v in g:
            if in_set[0, v]:
                a.add(v)
            else:
                b.add(v)
        val.append([a, b])

    return val


def vertex_connectivity(G, value_only=True, sets=False, k=None, solver=None, verbose=0,
                        *, integrality_tolerance=1e-3):
    r"""
    Return the vertex connectivity of the graph.

    For more information, see the :wikipedia:`Connectivity_(graph_theory)` and
    the :wikipedia:`K-vertex-connected_graph`.

    .. NOTE::

        * When the graph is directed, this method actually computes the *strong*
          connectivity, (i.e. a directed graph is strongly `k`-connected if
          there are `k` vertex disjoint paths between any two vertices `u,
          v`). If you do not want to consider strong connectivity, the best is
          probably to convert your ``DiGraph`` object to a ``Graph`` object, and
          compute the connectivity of this other graph.

        * By convention, a complete graph on `n` vertices is `n-1` connected. In
          this case, no certificate can be given as there is no pair of vertices
          split by a cut of order `k-1`. For this reason, the certificates
          returned in this situation are empty.

    INPUT:

    - ``G`` -- the input Sage (Di)Graph

    - ``value_only`` -- boolean (default: ``True``)

      - When set to ``True`` (default), only the value is returned.

      - When set to ``False``, both the value and a minimum vertex cut are
        returned.

    - ``sets`` -- boolean (default: ``False``); whether to also return the two
        sets of vertices that are disconnected by the cut (implies
        ``value_only=False``)

    - ``k`` -- integer (default: ``None``); when specified, check if the vertex
      connectivity of the (di)graph is larger or equal to `k`. The method thus
      outputs a boolean only.

    - ``solver`` -- string (default: ``None``); specifies a Mixed Integer Linear
      Programming (MILP) solver to be used. If set to ``None``, the default one
      is used. For more information on MILP solvers and which default solver is
      used, see the method :meth:`solve
      <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
      :class:`MixedIntegerLinearProgram
      <sage.numerical.mip.MixedIntegerLinearProgram>`.

    - ``verbose`` -- integer (default: 0); sets the level of verbosity. Set
      to 0 by default, which means quiet.

    - ``integrality_tolerance`` -- float; parameter for use with MILP solvers
      over an inexact base ring; see
      :meth:`MixedIntegerLinearProgram.get_values`.

    EXAMPLES:

    A basic application on a ``PappusGraph``::

       sage: from sage.graphs.connectivity import vertex_connectivity
       sage: g = graphs.PappusGraph()
       sage: vertex_connectivity(g)                                                     # needs sage.numerical.mip
       3
       sage: g.vertex_connectivity()                                                    # needs sage.numerical.mip
       3

    In a grid, the vertex connectivity is equal to the minimum degree, in which
    case one of the two sets is of cardinality `1`::

       sage: g = graphs.GridGraph([ 3,3 ])
       sage: [value, cut, [ setA, setB ]] = vertex_connectivity(g, sets=True)           # needs sage.numerical.mip
       sage: len(setA) == 1 or len(setB) == 1                                           # needs sage.numerical.mip
       True

    A vertex cut in a tree is any internal vertex::

       sage: tree = graphs.RandomTree(15)
       sage: val, [cut_vertex] = vertex_connectivity(tree, value_only=False)            # needs sage.numerical.mip
       sage: tree.degree(cut_vertex) > 1                                                # needs sage.numerical.mip
       True

    When ``value_only = True``, this function is optimized for small
    connectivity values and does not need to build a linear program.

    It is the case for connected graphs which are not connected::

       sage: g = 2 * graphs.PetersenGraph()
       sage: vertex_connectivity(g)                                                     # needs sage.numerical.mip
       0

    Or if they are just 1-connected::

       sage: g = graphs.PathGraph(10)
       sage: vertex_connectivity(g)                                                     # needs sage.numerical.mip
       1

    For directed graphs, the strong connectivity is tested through the dedicated
    function::

       sage: g = digraphs.ButterflyGraph(3)
       sage: vertex_connectivity(g)                                                     # needs sage.numerical.mip
       0

    A complete graph on `10` vertices is `9`-connected::

       sage: g = graphs.CompleteGraph(10)
       sage: vertex_connectivity(g)                                                     # needs sage.numerical.mip
       9

    A complete digraph on `10` vertices is `9`-connected::

       sage: g = DiGraph(graphs.CompleteGraph(10))
       sage: vertex_connectivity(g)                                                     # needs sage.numerical.mip
       9

    When parameter ``k`` is set, we only check for the existence of a vertex cut
    of order at least ``k``::

       sage: g = graphs.PappusGraph()
       sage: vertex_connectivity(g, k=3)                                                # needs sage.numerical.mip
       True
       sage: vertex_connectivity(g, k=4)                                                # needs sage.numerical.mip
       False

    TESTS:

    Giving negative value to parameter ``k``::

       sage: g = graphs.PappusGraph()
       sage: vertex_connectivity(g, k=-1)
       Traceback (most recent call last):
       ...
       ValueError: parameter k must be strictly positive

    The empty graph has vertex connectivity 0, is considered connected but not
    biconnected. The empty digraph is considered strongly connected::

       sage: from sage.graphs.connectivity import is_strongly_connected
       sage: from sage.graphs.connectivity import is_connected
       sage: empty = Graph()
       sage: vertex_connectivity(empty)                                                 # needs sage.numerical.mip
       0
       sage: vertex_connectivity(empty, k=1) == is_connected(empty)                     # needs sage.numerical.mip
       True
       sage: vertex_connectivity(Graph(), k=2) == empty.is_biconnected()                # needs sage.numerical.mip
       True
       sage: vertex_connectivity(DiGraph(), k=1) == is_strongly_connected(DiGraph())    # needs sage.numerical.mip
       True

    If ``G`` is not a Sage (Di)Graph, an error is raised::

        sage: vertex_connectivity('I am not a graph')
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage graph

    Complete Graph with loops or multiple edges (:issue:`25589`)::

        sage: G = Graph([(0, 1), (0, 1)], multiedges=True)
        sage: G.vertex_connectivity()                                                   # needs sage.numerical.mip
        1
        sage: G = graphs.CompleteGraph(4)
        sage: G.allow_loops(True)
        sage: G.add_edge(0, 0)
        sage: G.vertex_connectivity(value_only=False, verbose=1)                        # needs sage.numerical.mip
        (3, [])
        sage: G.allow_multiple_edges(True)
        sage: G.add_edge(0, 1)
        sage: G.vertex_connectivity(value_only=False, verbose=1)                        # needs sage.numerical.mip
        (3, [])
    """
    from sage.graphs.generic_graph import GenericGraph
    if not isinstance(G, GenericGraph):
        raise TypeError("the input must be a Sage graph")

    g = G

    if k is not None:
        if k < 1:
            raise ValueError("parameter k must be strictly positive")
        if not g.order():
            # We follow the convention of is_connected, is_biconnected and
            # is_strongly_connected
            return k == 1
        if (g.is_directed() and k > min(min(g.in_degree()), min(g.out_degree()))) \
           or (not g.is_directed() and (k > min(g.degree()))):
            return False
        value_only = True
        sets = False

    elif sets:
        value_only = False

    # When the graph is complete, the MILP below is infeasible.
    if (g.is_clique(directed_clique=g.is_directed())
            or (not g.is_directed() and g.to_simple().is_clique())):
        if k is not None:
            return g.order() > k
        if value_only:
            return max(g.order() - 1, 0)
        elif not sets:
            return max(g.order() - 1, 0), []
        return max(g.order() - 1, 0), [], [[], []]

    if value_only:
        if G.is_directed():
            if not is_strongly_connected(G):
                return 0 if k is None else False

        else:
            if not is_connected(G):
                return 0 if k is None else False

            if G.blocks_and_cut_vertices()[1]:
                return 1 if k is None else (k == 1)

            if not G.is_triconnected():
                return 2 if k is None else (k == 2)
            elif k == 3:
                return True

        if k == 1:
            # We know that the (di)graph is (strongly) connected
            return True

    from sage.numerical.mip import MixedIntegerLinearProgram, MIPSolverException

    p = MixedIntegerLinearProgram(maximization=False, solver=solver)

    # Sets 0 and 2 are "real" sets while set 1 represents the cut
    in_set = p.new_variable(binary=True)

    # A vertex has to be in some set
    for v in g:
        p.add_constraint(in_set[0, v] + in_set[1, v] + in_set[2, v], max=1, min=1)

    # There is no empty set
    p.add_constraint(p.sum(in_set[0, v] for v in g), min=1)
    p.add_constraint(p.sum(in_set[2, v] for v in g), min=1)

    if g.is_directed():
        # There is no edge from set 0 to set 1 which is not in the cut
        for u, v in g.edge_iterator(labels=None):
            p.add_constraint(in_set[0, u] + in_set[2, v], max=1)
    else:
        # Two adjacent vertices are in different sets if and only if
        # the edge between them is in the cut
        for u, v in g.edge_iterator(labels=None):
            p.add_constraint(in_set[0, u] + in_set[2, v], max=1)
            p.add_constraint(in_set[2, u] + in_set[0, v], max=1)

    if k is not None:
        # To check if the vertex connectivity is at least k, we check if
        # there exists a cut of order at most k-1. If the ILP is infeasible,
        # the vertex connectivity is >= k.
        p.add_constraint(p.sum(in_set[1, v] for v in g) <= k-1)
        try:
            p.solve(log=verbose)
            return False
        except MIPSolverException:
            return True

    p.set_objective(p.sum(in_set[1, v] for v in g))

    val = p.solve(log=verbose)

    in_set = p.get_values(in_set, convert=bool, tolerance=integrality_tolerance)

    if value_only:
        return sum(1 for v in g if in_set[1, v])

    cut = []
    a = []
    b = []

    for v in g:
        if in_set[0, v]:
            a.append(v)
        elif in_set[1, v]:
            cut.append(v)
        else:
            b.append(v)

    if sets:
        return val, cut, [a, b]

    return val, cut


def is_strongly_connected(G):
    r"""
    Check whether the current ``DiGraph`` is strongly connected.

    EXAMPLES:

    The circuit is obviously strongly connected::

        sage: from sage.graphs.connectivity import is_strongly_connected
        sage: g = digraphs.Circuit(5)
        sage: is_strongly_connected(g)
        True
        sage: g.is_strongly_connected()
        True

    But a transitive triangle is not::

        sage: g = DiGraph({0: [1, 2], 1: [2]})
        sage: is_strongly_connected(g)
        False

    TESTS:

    If ``G`` is not a Sage DiGraph, an error is raised::

        sage: is_strongly_connected('I am not a graph')
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage DiGraph
    """
    from sage.graphs.digraph import DiGraph
    if not isinstance(G, DiGraph):
        raise TypeError("the input must be a Sage DiGraph")

    if G.order() <= 1:
        return True

    try:
        return G._backend.is_strongly_connected()

    except AttributeError:
        return len(G.strongly_connected_components()) == 1


def strongly_connected_components_digraph(G, keep_labels=False):
    r"""
    Return the digraph of the strongly connected components.

    The digraph of the strongly connected components of a graph `G` has a vertex
    per strongly connected component included in `G`. There is an edge from a
    component `C_1` to a component `C_2` if there is an edge in `G` from a
    vertex `u_1 \in C_1` to a vertex `u_2 \in C_2`.

    INPUT:

    - ``G`` -- the input DiGraph

    - ``keep_labels`` -- boolean (default: ``False``); when
      ``keep_labels=True``, the resulting digraph has an edge from a component
      `C_i` to a component `C_j` for each edge in `G` from a vertex `u_i \in
      C_i` to a vertex `u_j \in C_j`. Hence the resulting digraph may have loops
      and multiple edges. However, edges in the result with same source, target,
      and label are not duplicated (see examples below). When
      ``keep_labels=False``, the return digraph is simple, so without loops nor
      multiple edges, and edges are unlabelled.

    EXAMPLES:

    Such a digraph is always acyclic::

        sage: from sage.graphs.connectivity import strongly_connected_components_digraph
        sage: g = digraphs.RandomDirectedGNP(15, .1)
        sage: scc_digraph = strongly_connected_components_digraph(g)
        sage: scc_digraph.is_directed_acyclic()
        True
        sage: scc_digraph = g.strongly_connected_components_digraph()
        sage: scc_digraph.is_directed_acyclic()
        True

    The vertices of the digraph of strongly connected components are exactly the
    strongly connected components::

        sage: g = digraphs.ButterflyGraph(2)
        sage: scc_digraph = strongly_connected_components_digraph(g)
        sage: g.is_directed_acyclic()
        True
        sage: V_scc = list(scc_digraph)
        sage: all(Set(scc) in V_scc for scc in g.strongly_connected_components())
        True

    The following digraph has three strongly connected components, and the
    digraph of those is a
    :meth:`~sage.graphs.digraph_generators.TransitiveTournament`::

        sage: g = DiGraph({0: {1: "01", 2: "02", 3: "03"}, 1: {2: "12"}, 2:{1: "21", 3: "23"}})
        sage: scc_digraph = strongly_connected_components_digraph(g)
        sage: scc_digraph.is_isomorphic(digraphs.TransitiveTournament(3))
        True

    By default, the labels are discarded, and the result has no loops nor
    multiple edges. If ``keep_labels`` is ``True``, then the labels are kept,
    and the result is a multi digraph, possibly with multiple edges and
    loops. However, edges in the result with same source, target, and label are
    not duplicated (see the edges from 0 to the strongly connected component
    `\{1,2\}` below)::

        sage: g = DiGraph({0: {1: "0-12", 2: "0-12", 3: "0-3"}, 1: {2: "1-2", 3: "1-3"}, 2: {1: "2-1", 3: "2-3"}})
        sage: g.order(), g.size()
        (4, 7)
        sage: scc_digraph = strongly_connected_components_digraph(g, keep_labels=True)
        sage: (scc_digraph.order(), scc_digraph.size())
        (3, 6)
        sage: set(g.edge_labels()) == set(scc_digraph.edge_labels())
        True

    TESTS:

    If ``G`` is not a Sage DiGraph, an error is raised::

        sage: strongly_connected_components_digraph('I am not a graph')
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage DiGraph
    """
    from sage.graphs.digraph import DiGraph
    if not isinstance(G, DiGraph):
        raise TypeError("the input must be a Sage DiGraph")

    from sage.sets.set import Set

    cdef list scc = G.strongly_connected_components()
    cdef list scc_set = [Set(_) for _ in scc]
    cdef dict d = {v: i for i, c in enumerate(scc) for v in c}

    if keep_labels:
        g = DiGraph(len(scc), multiedges=True, loops=True)
        g.add_edges(set((d[u], d[v], label) for u, v, label in G.edge_iterator()))

    else:
        g = DiGraph(len(scc), multiedges=False, loops=False)
        g.add_edges(((d[u], d[v]) for u, v in G.edge_iterator(labels=False)), loops=False)

    g.relabel(scc_set, inplace=True)
    return g


def strongly_connected_components_subgraphs(G):
    r"""
    Return the strongly connected components as a list of subgraphs.

    EXAMPLES:

    In the symmetric digraph of a graph, the strongly connected components are
    the connected components::

        sage: from sage.graphs.connectivity import strongly_connected_components_subgraphs
        sage: g = graphs.PetersenGraph()
        sage: d = DiGraph(g)
        sage: strongly_connected_components_subgraphs(d)
        [Subgraph of (Petersen graph): Digraph on 10 vertices]
        sage: d.strongly_connected_components_subgraphs()
        [Subgraph of (Petersen graph): Digraph on 10 vertices]

    ::

        sage: g = DiGraph([(0, 1), (1, 0), (1, 2), (2, 3), (3, 2)])
        sage: strongly_connected_components_subgraphs(g)
        [Subgraph of (): Digraph on 2 vertices, Subgraph of (): Digraph on 2 vertices]

    TESTS:

    If ``G`` is not a Sage DiGraph, an error is raised::

        sage: strongly_connected_components_subgraphs('I am not a graph')
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage DiGraph
    """
    from sage.graphs.digraph import DiGraph
    if not isinstance(G, DiGraph):
        raise TypeError("the input must be a Sage DiGraph")

    return [G.subgraph(_) for _ in G.strongly_connected_components()]


def strongly_connected_component_containing_vertex(G, v):
    """
    Return the strongly connected component containing a given vertex.

    INPUT:

    - ``G`` -- the input DiGraph

    - ``v`` -- a vertex

    EXAMPLES:

    In the symmetric digraph of a graph, the strongly connected components are
    the connected components::

        sage: from sage.graphs.connectivity import strongly_connected_component_containing_vertex
        sage: g = graphs.PetersenGraph()
        sage: d = DiGraph(g)
        sage: strongly_connected_component_containing_vertex(d, 0)
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        sage: d.strongly_connected_component_containing_vertex(0)
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    ::

        sage: g = DiGraph([(0, 1), (1, 0), (1, 2), (2, 3), (3, 2)])
        sage: strongly_connected_component_containing_vertex(g, 0)
        [0, 1]

    TESTS:

    If ``G`` is not a Sage DiGraph, an error is raised::

        sage: strongly_connected_component_containing_vertex('I am not a graph', 0)
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage DiGraph

    If the vertex is not in the DiGraph::

        sage: strongly_connected_component_containing_vertex(DiGraph(1), 'z')
        Traceback (most recent call last):
        ...
        ValueError: vertex ('z') is not a vertex of the DiGraph
    """
    from sage.graphs.digraph import DiGraph
    if not isinstance(G, DiGraph):
        raise TypeError("the input must be a Sage DiGraph")

    if v not in G:
        raise ValueError("vertex ({0}) is not a vertex of the DiGraph".format(repr(v)))

    if G.order() == 1:
        return [v]

    try:
        return G._backend.strongly_connected_component_containing_vertex(v)

    except AttributeError:
        raise AttributeError("this function is only defined for C graphs")


def strong_articulation_points(G):
    r"""
    Return the strong articulation points of this digraph.

    A vertex is a strong articulation point if its deletion increases the
    number of strongly connected components. This method implements the
    algorithm described in [ILS2012]_. The time complexity is dominated by
    the time complexity of the immediate dominators finding algorithm.

    OUTPUT: the list of strong articulation points

    EXAMPLES:

    Two cliques sharing a vertex::

        sage: from sage.graphs.connectivity import strong_articulation_points
        sage: D = digraphs.Complete(4)
        sage: D.add_clique([3, 4, 5, 6])
        sage: strong_articulation_points(D)
        [3]
        sage: D.strong_articulation_points()
        [3]

    Two cliques connected by some arcs::

        sage: D = digraphs.Complete(4) * 2
        sage: D.add_edges([(0, 4), (7, 3)])
        sage: sorted(strong_articulation_points(D))
        [0, 3, 4, 7]
        sage: D.add_edge(1, 5)
        sage: sorted(strong_articulation_points(D))
        [3, 7]
        sage: D.add_edge(6, 2)
        sage: strong_articulation_points(D)
        []

    .. SEEALSO::

        - :meth:`~sage.graphs.digraph.DiGraph.strongly_connected_components`
        - :meth:`~sage.graphs.base.boost_graph.dominator_tree`

    TESTS:

    All strong articulation points are found::

        sage: from sage.graphs.connectivity import strong_articulation_points
        sage: def sap_naive(G):
        ....:     nscc = len(G.strongly_connected_components())
        ....:     S = []
        ....:     for u in G:
        ....:         H = copy(G)
        ....:         H.delete_vertex(u)
        ....:         if len(H.strongly_connected_components()) > nscc:
        ....:             S.append(u)
        ....:     return S
        sage: D = digraphs.RandomDirectedGNP(20, 0.1)
        sage: X = sap_naive(D)
        sage: SAP = strong_articulation_points(D)
        sage: set(X) == set(SAP)
        True

    Trivial cases::

        sage: strong_articulation_points(DiGraph())
        []
        sage: strong_articulation_points(DiGraph(1))
        []
        sage: strong_articulation_points(DiGraph(2))
        []

    If ``G`` is not a Sage DiGraph, an error is raised::

        sage: strong_articulation_points('I am not a graph')
        Traceback (most recent call last):
        ...
        TypeError: the input must be a Sage DiGraph

    Issue :issue:`29958` is fixed::

        sage: D = DiGraph('SA?GA??_??a???@?@OH_?@?I??b??G?AgGGCO??AC????a?????A@????AOCOQ?d??I?')
        sage: SAP = strong_articulation_points(D)
        sage: set(SAP) == {1, 2, 4, 17, 18}
        True
    """
    from sage.graphs.digraph import DiGraph
    if not isinstance(G, DiGraph):
        raise TypeError("the input must be a Sage DiGraph")

    # The method is applied on each strongly connected component
    if is_strongly_connected(G):
        # Make a mutable copy of self
        L = [DiGraph([(u, v) for u, v in G.edge_iterator(labels=0) if u != v],
                     data_structure='sparse', immutable=False)]
    else:
        # Get the list of strongly connected components of self as mutable
        # subgraphs
        L = [G.subgraph(scc, immutable=False) for scc in G.strongly_connected_components()]

    SAP = []
    for g in L:
        n = g.order()
        if n <= 2:
            continue

        # 1. Choose arbitrarily a vertex r, and test whether r is a strong
        # articulation point.
        r = next(g.vertex_iterator())
        E = g.incoming_edges(r) + g.outgoing_edges(r)
        g.delete_vertex(r)
        if not is_strongly_connected(g):
            SAP.append(r)
        g.add_edges(E)

        # 2. Compute the set of non-trivial immediate dominators in g
        Dr = set(g.dominator_tree(r, return_dict=True).values())

        # 3. Compute the set of non-trivial immediate dominators in the
        # reverse digraph
        DRr = set(g.dominator_tree(r, return_dict=True, reverse=True).values())

        # 4. Store D(r) + DR(r) - r
        SAP.extend(Dr.union(DRr).difference([r, None]))

    return SAP


def bridges(G, labels=True):
    r"""
    Return an iterator over the bridges (or cut edges).

    A bridge is an edge whose deletion disconnects the undirected graph.
    A disconnected graph has no bridge.

    INPUT:

    - ``labels`` -- boolean (default: ``True``); if ``False``, each bridge is a
      tuple `(u, v)` of vertices

    EXAMPLES::

        sage: from sage.graphs.connectivity import bridges
        sage: from sage.graphs.connectivity import is_connected
        sage: g = 2 * graphs.PetersenGraph()
        sage: g.add_edge(1, 10)
        sage: is_connected(g)
        True
        sage: list(bridges(g))
        [(1, 10, None)]
        sage: list(g.bridges())
        [(1, 10, None)]

    Every edge of a tree is a bridge::

        sage: g = graphs.RandomTree(100)
        sage: sum(1 for _ in g.bridges()) == 99
        True

    TESTS:

    Disconnected graphs have no bridges::

        sage: g = 2*graphs.PetersenGraph()
        sage: next(g.bridges())
        Traceback (most recent call last):
        ...
        StopIteration

    Graph with multiple edges and edge labels::

        sage: g = 2 * graphs.CycleGraph(3)
        sage: g.allow_multiple_edges(True)
        sage: g.add_edges(g.edges(sort=False))
        sage: g.add_edge(2, 3, "label")
        sage: list(bridges(g, labels=True))
        [(2, 3, 'label')]

    Issue :issue:`23817` is solved::

        sage: G = Graph()
        sage: G.add_edge(0, 1)
        sage: list(bridges(G))
        [(0, 1, None)]
        sage: G.allow_loops(True)
        sage: G.add_edge(0, 0)
        sage: G.add_edge(1, 1)
        sage: list(bridges(G))
        [(0, 1, None)]

    If ``G`` is not a Sage Graph, an error is raised::

        sage: next(bridges('I am not a graph'))
        Traceback (most recent call last):
        ...
        TypeError: the input must be an undirected Sage graph
    """
    from sage.graphs.graph import Graph
    if not isinstance(G, Graph):
        raise TypeError("the input must be an undirected Sage graph")

    # Small graphs and disconnected graphs have no bridge
    if G.order() < 2 or not is_connected(G):
        return

    B, _ = G.blocks_and_cut_vertices()

    # A block of size 2 is a bridge, unless the vertices are connected with
    # multiple edges.
    cdef bint multiple_edges = G.allows_multiple_edges()
    cdef set ME = set(G.multiple_edges(labels=False)) if multiple_edges else set()
    for b in B:
        if len(b) == 2 and not tuple(b) in ME:
            if labels:
                if multiple_edges:
                    [label] = G.edge_label(b[0], b[1])
                else:
                    label = G.edge_label(b[0], b[1])
                yield (b[0], b[1], label)
            else:
                yield tuple(b)


# ==============================================================================
# Methods for finding 3-vertex-connected components and building SPQR-tree
# ==============================================================================

def cleave(G, cut_vertices=None, virtual_edges=True, solver=None, verbose=0,
           *, integrality_tolerance=1e-3):
    r"""
    Return the connected subgraphs separated by the input vertex cut.

    Given a connected (multi)graph `G` and a vertex cut `X`, this method
    computes the list of subgraphs of `G` induced by each connected component
    `c` of `G\setminus X` plus `X`, i.e., `G[c\cup X]`.

    INPUT:

    - ``G`` -- a Graph

    - ``cut_vertices`` -- iterable container of vertices (default: ``None``); a
      set of vertices representing a vertex cut of ``G``. If no vertex cut is
      given, the method will compute one via a call to
      :meth:`~sage.graphs.connectivity.vertex_connectivity`.

    - ``virtual_edges`` -- boolean (default: ``True``); whether to add virtual
      edges to the sides of the cut or not. A virtual edge is an edge between a
      pair of vertices of the cut that are not connected by an edge in ``G``.

    - ``solver`` -- string (default: ``None``); specifies a Mixed Integer Linear
      Programming (MILP) solver to be used. If set to ``None``, the default one
      is used. For more information on MILP solvers and which default solver is
      used, see the method :meth:`solve
      <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
      :class:`MixedIntegerLinearProgram
      <sage.numerical.mip.MixedIntegerLinearProgram>`.

    - ``verbose`` -- integer (default: 0); sets the level of verbosity. Set
      to 0 by default, which means quiet.

    - ``integrality_tolerance`` -- float; parameter for use with MILP solvers
      over an inexact base ring; see
      :meth:`MixedIntegerLinearProgram.get_values`.

    OUTPUT: a triple `(S, C, f)`, where

    - `S` is a list of the graphs that are sides of the vertex cut.

    - `C` is the graph of the cocycles. For each pair of vertices of the cut,
      if there exists an edge between them, `C` has one copy of each edge
      connecting them in ``G`` per sides of the cut plus one extra copy.
      Furthermore, when ``virtual_edges == True``, if a pair of vertices of the
      cut is not connected by an edge in ``G``, then it has one virtual edge
      between them per sides of the cut.

    - `f` is the complement of the subgraph of ``G`` induced by the vertex
      cut. Hence, its vertex set is the vertex cut, and its edge set is the set
      of virtual edges (i.e., edges between pairs of vertices of the cut that
      are not connected by an edge in ``G``). When ``virtual_edges == False``,
      the edge set is empty.

    EXAMPLES:

    If there is an edge between cut vertices::

        sage: from sage.graphs.connectivity import cleave
        sage: G = Graph(2)
        sage: for _ in range(3):
        ....:     G.add_clique([0, 1, G.add_vertex(), G.add_vertex()])
        sage: S1,C1,f1 = cleave(G, cut_vertices=[0, 1])
        sage: [g.order() for g in S1]
        [4, 4, 4]
        sage: C1.order(), C1.size()
        (2, 4)
        sage: f1.vertices(sort=True), f1.edges(sort=True)
        ([0, 1], [])

    If ``virtual_edges == False`` and there is an edge between cut vertices::

        sage: G.subgraph([0, 1]).complement() == Graph([[0, 1], []])
        True
        sage: S2,C2,f2 = cleave(G, cut_vertices=[0, 1], virtual_edges=False)
        sage: (S1 == S2, C1 == C2, f1 == f2)
        (True, True, True)

    If cut vertices doesn't have edge between them::

        sage: G.delete_edge(0, 1)
        sage: S1,C1,f1 = cleave(G, cut_vertices=[0, 1])
        sage: [g.order() for g in S1]
        [4, 4, 4]
        sage: C1.order(), C1.size()
        (2, 3)
        sage: f1.vertices(sort=True), f1.edges(sort=True)
        ([0, 1], [(0, 1, None)])

    If ``virtual_edges == False`` and the cut vertices are not connected by an
    edge::

        sage: G.subgraph([0, 1]).complement() == Graph([[0, 1], []])
        False
        sage: S2,C2,f2 = cleave(G, cut_vertices=[0, 1], virtual_edges=False)
        sage: [g.order() for g in S2]
        [4, 4, 4]
        sage: C2.order(), C2.size()
        (2, 0)
        sage: f2.vertices(sort=True), f2.edges(sort=True)
        ([0, 1], [])
        sage: (S1 == S2, C1 == C2, f1 == f2)
        (False, False, False)

    If `G` is a biconnected multigraph::

        sage: G = graphs.CompleteBipartiteGraph(2, 3)
        sage: G.add_edge(2, 3)
        sage: G.allow_multiple_edges(True)
        sage: G.add_edges(G.edge_iterator())
        sage: G.add_edges([(0, 1), (0, 1), (0, 1)])
        sage: S,C,f = cleave(G, cut_vertices=[0, 1])
        sage: for g in S:
        ....:     print(g.edges(sort=True, labels=0))
        [(0, 1), (0, 1), (0, 1), (0, 2), (0, 2), (0, 3), (0, 3), (1, 2), (1, 2), (1, 3), (1, 3), (2, 3), (2, 3)]
        [(0, 1), (0, 1), (0, 1), (0, 4), (0, 4), (1, 4), (1, 4)]

    TESTS::

        sage: cleave(Graph(2))
        Traceback (most recent call last):
        ...
        ValueError: this method is designed for connected graphs only
        sage: cleave(graphs.PathGraph(3), cut_vertices=[5])
        Traceback (most recent call last):
        ...
        ValueError: vertex 5 is not a vertex of the input graph
        sage: cleave(Graph())
        Traceback (most recent call last):
        ...
        ValueError: the input graph has no vertex cut
        sage: cleave(graphs.CompleteGraph(5))
        Traceback (most recent call last):
        ...
        ValueError: the input graph has no vertex cut
        sage: cleave(graphs.CompleteGraph(5), cut_vertices=[3, 4])
        Traceback (most recent call last):
        ...
        ValueError: the set cut_vertices is not a vertex cut of the graph
    """
    if not G.is_connected():
        raise ValueError("this method is designed for connected graphs only")

    # If a vertex cut is given, we check that it is valid. Otherwise, we compute
    # a small vertex cut
    if cut_vertices is None:
        _, cut_vertices = G.vertex_connectivity(value_only=False,
                                                solver=solver, verbose=verbose,
                                                integrality_tolerance=integrality_tolerance)
        if not cut_vertices:
            # Typical example is a clique
            raise ValueError("the input graph has no vertex cut")
    else:
        cut_vertices = list(cut_vertices)
        for u in cut_vertices:
            if u not in G:
                raise ValueError(f"vertex {u} is not a vertex of the input graph")

    H = G.copy(immutable=False)
    H.delete_vertices(cut_vertices)
    CC = H.connected_components(sort=False)
    if len(CC) == 1:
        raise ValueError("the set cut_vertices is not a vertex cut of the graph")

    # We identify the virtual edges, i.e., pairs of vertices of the vertex cut
    # that are not connected by an edge in G
    from sage.graphs.graph import Graph
    K = G.subgraph(cut_vertices)
    if virtual_edges:
        if K.allows_multiple_edges():
            virtual_cut_graph = K.to_simple().complement()
        else:
            virtual_cut_graph = K.complement()
    else:
        virtual_cut_graph = Graph([cut_vertices, []])

    # We now build the graphs in each side of the cut, including the vertices
    # from the vertex cut
    cut_sides = []
    for comp in CC:
        h = G.subgraph(comp + cut_vertices)
        if virtual_edges:
            h.add_edges(virtual_cut_graph.edge_iterator())
        cut_sides.append(h)

    # We build the cocycles for re-assembly. For each edge between a pair of
    # vertices of the cut in the original graph G, a bond with one edge more
    # than the number of cut sides is needed. For pairs of vertices of the cut
    # that are not connected by an edge in G, a bond with one edge per cut side
    # is needed.
    cocycles = Graph([cut_vertices, []], multiedges=True)
    if K.size():
        cocycles.add_edges(K.edges(sort=False) * (len(cut_sides) + 1))
    if virtual_edges and virtual_cut_graph:
        cocycles.add_edges(virtual_cut_graph.edges(sort=False) * len(cut_sides))

    return cut_sides, cocycles, virtual_cut_graph


def spqr_tree(G, algorithm='Hopcroft_Tarjan', solver=None, verbose=0,
              *, integrality_tolerance=1e-3):
    r"""
    Return an SPQR-tree representing the triconnected components of the graph.

    An SPQR-tree is a tree data structure used to represent the triconnected
    components of a biconnected (multi)graph and the 2-vertex cuts separating
    them. A node of a SPQR-tree, and the graph associated with it, can be one of
    the following four types:

    - ``'S'`` -- the associated graph is a cycle with at least three vertices
      ``'S'`` stands for ``series``

    - ``'P'`` -- the associated graph is a dipole graph, a multigraph with two
      vertices and three or more edges. ``'P'`` stands for ``parallel``

    - ``'Q'`` -- the associated graph has a single real edge. This trivial case
      is necessary to handle the graph that has only one edge

    - ``'R'`` -- the associated graph is a 3-connected graph that is not a cycle
      or dipole. ``'R'`` stands for ``rigid``

    This method decomposes a biconnected graph into cycles, cocycles, and
    3-connected blocks summed over cocycles, and arranges them as a SPQR-tree.
    More precisely, it splits the graph at each of its 2-vertex cuts, giving a
    unique decomposition into 3-connected blocks, cycles and cocycles. The
    cocycles are dipole graphs with one edge per real edge between the included
    vertices and one additional (virtual) edge per connected component resulting
    from deletion of the vertices in the cut. See the :wikipedia:`SPQR_tree`.

    INPUT:

    - ``G`` -- the input graph

    - ``algorithm`` -- string (default: ``'Hopcroft_Tarjan'``); the algorithm to
      use among:

      - ``'Hopcroft_Tarjan'`` -- default; use the algorithm proposed by
        Hopcroft and Tarjan in [Hopcroft1973]_ and later corrected by Gutwenger
        and Mutzel in [Gut2001]_. See
        :class:`~sage.graphs.connectivity.TriconnectivitySPQR`.

      - ``'cleave'`` -- using method :meth:`~sage.graphs.connectivity.cleave`

    - ``solver`` -- string (default: ``None``); specifies a Mixed Integer Linear
      Programming (MILP) solver to be used. If set to ``None``, the default one
      is used. For more information on MILP solvers and which default solver is
      used, see the method :meth:`solve
      <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
      :class:`MixedIntegerLinearProgram
      <sage.numerical.mip.MixedIntegerLinearProgram>`.

    - ``verbose`` -- integer (default: 0); sets the level of verbosity. Set
      to 0 by default, which means quiet.

    - ``integrality_tolerance`` -- float; parameter for use with MILP solvers
      over an inexact base ring; see
      :meth:`MixedIntegerLinearProgram.get_values`.

    OUTPUT: ``SPQR-tree`` a tree whose vertices are labeled with the block's type
    and the subgraph of three-blocks in the decomposition.

    EXAMPLES::

        sage: from sage.graphs.connectivity import spqr_tree
        sage: G = Graph(2)
        sage: for i in range(3):
        ....:     G.add_clique([0, 1, G.add_vertex(), G.add_vertex()])
        sage: Tree = spqr_tree(G)
        sage: Tree.order()
        4
        sage: K4 = graphs.CompleteGraph(4)
        sage: all(u[1].is_isomorphic(K4) for u in Tree if u[0] == 'R')
        True
        sage: from sage.graphs.connectivity import spqr_tree_to_graph
        sage: G.is_isomorphic(spqr_tree_to_graph(Tree))
        True

        sage: G = Graph(2)
        sage: for i in range(3):
        ....:     G.add_path([0, G.add_vertex(), G.add_vertex(), 1])
        sage: Tree = spqr_tree(G)
        sage: Tree.order()
        4
        sage: C4 = graphs.CycleGraph(4)
        sage: all(u[1].is_isomorphic(C4) for u in Tree if u[0] == 'S')
        True
        sage: G.is_isomorphic(spqr_tree_to_graph(Tree))
        True

        sage: G.allow_multiple_edges(True)
        sage: G.add_edges(G.edge_iterator())
        sage: Tree = spqr_tree(G)
        sage: Tree.order()
        13
        sage: all(u[1].is_isomorphic(C4) for u in Tree if u[0] == 'S')
        True
        sage: G.is_isomorphic(spqr_tree_to_graph(Tree))
        True

        sage: G = graphs.CycleGraph(6)
        sage: Tree = spqr_tree(G)
        sage: Tree.order()
        1
        sage: G.is_isomorphic(spqr_tree_to_graph(Tree))
        True
        sage: G.add_edge(0, 3)
        sage: Tree = spqr_tree(G)
        sage: Tree.order()
        3
        sage: G.is_isomorphic(spqr_tree_to_graph(Tree))
        True

        sage: G = Graph('LlCG{O@?GBoMw?')
        sage: T = spqr_tree(G, algorithm='Hopcroft_Tarjan')
        sage: G.is_isomorphic(spqr_tree_to_graph(T))
        True
        sage: T2 = spqr_tree(G, algorithm='cleave')                                     # needs sage.numerical.mip
        sage: G.is_isomorphic(spqr_tree_to_graph(T2))                                   # needs sage.numerical.mip
        True

        sage: G = Graph([(0, 1)], multiedges=True)
        sage: T = spqr_tree(G, algorithm='cleave')                                      # needs sage.numerical.mip
        sage: T.vertices(sort=True)                                                     # needs sage.numerical.mip
        [('Q', Multi-graph on 2 vertices)]
        sage: G.is_isomorphic(spqr_tree_to_graph(T))                                    # needs sage.numerical.mip
        True
        sage: T = spqr_tree(G, algorithm='Hopcroft_Tarjan')
        sage: T.vertices(sort=True)
        [('Q', Multi-graph on 2 vertices)]
        sage: G.add_edge(0, 1)
        sage: spqr_tree(G, algorithm='cleave').vertices(sort=True)                      # needs sage.numerical.mip
        [('P', Multi-graph on 2 vertices)]

        sage: from collections import Counter
        sage: G = graphs.PetersenGraph()
        sage: T = G.spqr_tree(algorithm='Hopcroft_Tarjan')
        sage: Counter(u[0] for u in T)
        Counter({'R': 1})
        sage: T = G.spqr_tree(algorithm='cleave')                                       # needs sage.numerical.mip
        sage: Counter(u[0] for u in T)                                                  # needs sage.numerical.mip
        Counter({'R': 1})
        sage: for u,v in list(G.edges(labels=False, sort=False)):
        ....:     G.add_path([u, G.add_vertex(), G.add_vertex(), v])
        sage: T = G.spqr_tree(algorithm='Hopcroft_Tarjan')
        sage: sorted(Counter(u[0] for u in T).items())
        [('P', 15), ('R', 1), ('S', 15)]
        sage: T = G.spqr_tree(algorithm='cleave')                                       # needs sage.numerical.mip
        sage: sorted(Counter(u[0] for u in T).items())                                  # needs sage.numerical.mip
        [('P', 15), ('R', 1), ('S', 15)]
        sage: for u,v in list(G.edges(labels=False, sort=False)):
        ....:     G.add_path([u, G.add_vertex(), G.add_vertex(), v])
        sage: T = G.spqr_tree(algorithm='Hopcroft_Tarjan')
        sage: sorted(Counter(u[0] for u in T).items())
        [('P', 60), ('R', 1), ('S', 75)]
        sage: T = G.spqr_tree(algorithm='cleave')       # long time                     # needs sage.numerical.mip
        sage: sorted(Counter(u[0] for u in T).items())  # long time                     # needs sage.numerical.mip
        [('P', 60), ('R', 1), ('S', 75)]

    TESTS::

        sage: G = graphs.PathGraph(4)
        sage: spqr_tree(G)
        Traceback (most recent call last):
        ...
        ValueError: graph is not biconnected

        sage: G = Graph([(0, 0)], loops=True)
        sage: spqr_tree(G)
        Traceback (most recent call last):
        ...
        ValueError: graph is not biconnected

        sage: spqr_tree(Graph(), algorithm='easy')
        Traceback (most recent call last):
        ...
        NotImplementedError: SPQR tree algorithm 'easy' is not implemented
    """
    from sage.graphs.generic_graph import GenericGraph
    if not isinstance(G, GenericGraph):
        raise TypeError("the input must be a Sage graph")

    if algorithm == "Hopcroft_Tarjan":
        tric = TriconnectivitySPQR(G)
        return tric.get_spqr_tree()

    if algorithm != "cleave":
        raise NotImplementedError("SPQR tree algorithm '{}' is not implemented".format(algorithm))

    from sage.graphs.graph import Graph
    from collections import Counter

    if G.has_loops():
        raise ValueError("generation of SPQR-trees is only implemented for graphs without loops")

    if G.order() == 2 and G.size():
        return Graph({('Q' if G.size() == 1 else 'P', Graph(G, immutable=True, multiedges=True)): []},
                     name='SPQR-tree of {}'.format(G.name()))

    cut_size, cut_vertices = G.vertex_connectivity(value_only=False, solver=solver, verbose=verbose,
                                                   integrality_tolerance=integrality_tolerance)

    if cut_size < 2:
        raise ValueError("generation of SPQR-trees is only implemented for 2-connected graphs")
    elif cut_size > 2:
        return Graph({('R', Graph(G, immutable=True)): []}, name='SPQR-tree of {}'.format(G.name()))

    # Split_multiple_edge Algorithm. If the input graph has multiple edges, we
    # make SG a simple graph while recording virtual edges that will be needed
    # to 2-sum with cocycles during reconstruction. After this step, next steps
    # will be finding S and R blocks.
    if G.has_multiple_edges():
        SG = G.to_simple()
        counter_multiedges = Counter([frozenset(e) for e in G.multiple_edges(labels=False)])
    else:
        SG = G
        counter_multiedges = Counter()

    # If SG simplifies to a cycle, we return the corresponding SPQR tree
    if SG.is_cycle():
        T = Graph(name='SPQR-tree of {}'.format(G.name()))
        S_node = ('S', Graph(SG, immutable=True))
        T.add_vertex(S_node)
        for e, num in counter_multiedges.items():
            P_node = ('P', Graph([e] * (num + 1), multiedges=True, immutable=True))
            T.add_edge(S_node, P_node)
        return T

    # We now search for 2-vertex cuts. If a cut is found, we split the graph and
    # check each side for S or R blocks or another 2-vertex cut
    R_blocks = []
    two_blocks = [(SG, cut_vertices)]
    cocycles_count = Counter()
    cycles_list = []
    virtual_edge_to_cycles = dict()

    while two_blocks:
        B, B_cut = two_blocks.pop()
        # B will be always simple graph.
        S, C, f = cleave(B, cut_vertices=B_cut)
        # Store the number of edges of the cocycle (P block)
        fe = frozenset(B_cut)
        cocycles_count[fe] += C.size()
        if f.size():
            virtual_edge_to_cycles[fe] = []
        # Check the sides of the cut
        for K in S:
            if K.is_cycle():
                # Add this cycle to the list of cycles
                cycles_list.append(K)
            else:
                K_cut_size, K_cut_vertices = K.vertex_connectivity(value_only=False, solver=solver, verbose=verbose,
                                                                   integrality_tolerance=integrality_tolerance)
                if K_cut_size == 2:
                    # The graph has a 2-vertex cut. We add it to the stack
                    two_blocks.append((K, K_cut_vertices))
                else:
                    # The graph is 3-vertex connected
                    R_blocks.append(('R', Graph(K, immutable=True)))

    # Cycles of order > 3 may have been triangulated; We undo this to reduce the
    # number of S-blocks. Two cycles can be merged if they share a virtual edge
    # that is not shared by any other block, i.e., cocycles_count[e] == 2. We
    # first associate cycles to virtual edges. Then, we use a DisjointSet to
    # form the groups of cycles to be merged.
    for K_index, K in enumerate(cycles_list):
        for e in K.edge_iterator(labels=False):
            fe = frozenset(e)
            if fe in virtual_edge_to_cycles:
                virtual_edge_to_cycles[fe].append(K_index)
    from sage.sets.disjoint_set import DisjointSet
    DS = DisjointSet(range(len(cycles_list)))
    for fe in virtual_edge_to_cycles:
        if cocycles_count[fe] == 2 and len(virtual_edge_to_cycles[fe]) == 2:
            # This virtual edge is only between 2 cycles
            C1, C2 = virtual_edge_to_cycles[fe]
            DS.union(C1, C2)
            cycles_list[C1].delete_edge(fe)
            cycles_list[C2].delete_edge(fe)
            cocycles_count[fe] -= 2

    # We finalize the creation of S_blocks.
    S_blocks = []
    for root, indexes in DS.root_to_elements_dict().items():
        E = []
        for i in indexes:
            E.extend(cycles_list[i].edge_iterator(labels=False))
        S_blocks.append(('S', Graph(E, immutable=True)))

    # We now build the SPQR tree
    Tree = Graph(name='SPQR tree of {}'.format(G.name()))
    SR_blocks = S_blocks + R_blocks
    Tree.add_vertices(SR_blocks)
    P2 = []
    for e, num in cocycles_count.items():
        if num:
            P_block = ('P', Graph([e] * (num + max(0, counter_multiedges[e] - 1)), multiedges=True, immutable=True))
            for block in SR_blocks:
                # Note: here we use a try...except statement since the immutable
                # graph backend raises an error if an end vertex of the edge is
                # not in the graph.
                try:
                    if block[1].has_edge(e):
                        Tree.add_edge(block, P_block)
                except LookupError:
                    continue
            if num == 2:
                # When 2 S or R blocks are separated by a 2-cut without edge, we
                # have added a P block with only 2 edges. We must remove them
                # and connect neighbors by an edge. So we record these blocks
                P2.append(P_block)

    # We now remove the P blocks with only 2 edges.
    for P_block in P2:
        u, v = Tree.neighbors(P_block)
        Tree.add_edge(u, v)
        Tree.delete_vertex(P_block)

    # We finally add P blocks to account for multiple edges of the input graph
    # that are not involved in any separator of the graph
    for e, num in counter_multiedges.items():
        if not cocycles_count[e]:
            P_block = ('P', Graph([e] * (num + 1), multiedges=True, immutable=True))
            for block in SR_blocks:
                try:
                    if block[1].has_edge(e):
                        Tree.add_edge(block, P_block)
                        break
                except LookupError:
                    continue

    return Tree


def spqr_tree_to_graph(T):
    r"""
    Return the graph represented by the SPQR-tree `T`.

    The main purpose of this method is to test :meth:`spqr_tree`.

    INPUT:

    - ``T`` -- a SPQR tree as returned by :meth:`spqr_tree`

    OUTPUT: a (multi) graph

    EXAMPLES:

    :wikipedia:`SPQR_tree` reference paper example::

        sage: from sage.graphs.connectivity import spqr_tree
        sage: from sage.graphs.connectivity import spqr_tree_to_graph
        sage: G = Graph([(1, 2), (1, 4), (1, 8), (1, 12), (3, 4), (2, 3),
        ....: (2, 13), (3, 13), (4, 5), (4, 7), (5, 6), (5, 8), (5, 7), (6, 7),
        ....: (8, 11), (8, 9), (8, 12), (9, 10), (9, 11), (9, 12), (10, 12)])
        sage: T = spqr_tree(G)
        sage: H = spqr_tree_to_graph(T)
        sage: H.is_isomorphic(G)
        True

    A small multigraph ::

        sage: G = Graph([(0, 2), (0, 2), (1, 3), (2, 3)], multiedges=True)
        sage: for i in range(3):
        ....:     G.add_clique([0, 1, G.add_vertex(), G.add_vertex()])
        sage: for i in range(3):
        ....:     G.add_clique([2, 3, G.add_vertex(), G.add_vertex()])
        sage: T = spqr_tree(G)
        sage: H = spqr_tree_to_graph(T)
        sage: H.is_isomorphic(G)
        True

    TESTS:

    Check that the method is working for the empty SPQR tree::

        sage: H = spqr_tree_to_graph(Graph())
        sage: H.is_isomorphic(Graph())
        True

    Check that :issue:`38527` is fixed::

        sage: from sage.graphs.connectivity import spqr_tree, spqr_tree_to_graph
        sage: G = Graph('LlCG{O@?GBoMw?')
        sage: T1 = spqr_tree(G, algorithm="Hopcroft_Tarjan")
        sage: T2 = spqr_tree(G, algorithm="cleave")
        sage: T1.is_isomorphic(T2)
        True
        sage: G1 = spqr_tree_to_graph(T1)
        sage: G2 = spqr_tree_to_graph(T2)
        sage: G.is_isomorphic(G1)
        True
        sage: G.is_isomorphic(G2)
        True
    """
    from sage.graphs.graph import Graph
    from collections import Counter

    count_G = Counter()
    count_P = Counter()
    vertex_to_int = dict()
    for t, g in T:
        for u in g:
            if u not in vertex_to_int:
                vertex_to_int[u] = len(vertex_to_int)
        if t in ['P', 'Q']:
            for u, v, label in g.edge_iterator():
                if vertex_to_int[u] < vertex_to_int[v]:
                    count_P[u, v, label] += 1
                else:
                    count_P[v, u, label] += 1
        else:
            for u, v, label in g.edge_iterator():
                if vertex_to_int[u] < vertex_to_int[v]:
                    count_G[u, v, label] += 1
                else:
                    count_G[v, u, label] += 1

    G = Graph(multiedges=True)
    for e, num in count_G.items():
        if e in count_P:
            num = abs(count_P[e] - count_G[e])
        elif num == 2:
            # Case of 2 S or R blocks separated by a 2-cut without edge.
            # No P-block was created as a P-block has at least 3 edges.
            continue
        for _ in range(num):
            G.add_edge(e)

    # Some edges might only be in P_blocks. Such edges are true edges of the
    # graph. This happen when virtual edges have distinct labels.
    for e, num in count_P.items():
        if e not in count_G:
            for _ in range(num):
                G.add_edge(e)

    return G


# Helper methods for ``TriconnectivitySPQR``.
# Define a doubly linked list

cdef inline _LinkedListNode_initialize(_LinkedListNode * node, Py_ssize_t data):
    """
    Initialize the ``_LinkedListNode`` with value data.
    """
    node.prev = NULL
    node.next = NULL
    node.data = data


cdef inline _LinkedList_initialize(_LinkedList * ll):
    """
    Initialize the ``_LinkedList``.
    """
    ll.head = NULL
    ll.tail = NULL
    ll.length = 0

cdef _LinkedList_set_head(_LinkedList * ll, _LinkedListNode * h):
    """
    Set the node ``h`` as the head and tail of the linked list ``ll``.
    """
    ll.head = h
    ll.tail = h
    ll.length = 1

cdef inline _LinkedListNode * _LinkedList_get_head(_LinkedList * ll) noexcept:
    """
    Return the head of the linked list ``ll``.
    """
    return ll.head

cdef inline Py_ssize_t _LinkedList_get_length(_LinkedList * ll) noexcept:
    """
    Return the length of the linked list ``ll``.
    """
    return ll.length

cdef _LinkedList_append(_LinkedList * ll, _LinkedListNode * node):
    """
    Append the node ``node`` to the linked list ``ll``.
    """
    if not ll.head:
        _LinkedList_set_head(ll, node)
    else:
        ll.tail.next = node
        node.prev = ll.tail
        ll.tail = node
        ll.length += 1

cdef _LinkedList_remove(_LinkedList * ll, _LinkedListNode * node):
    """
    Remove the node ``node`` from the linked list ``ll``.
    """
    if not node.prev and not node.next:
        ll.head = NULL
        ll.tail = NULL
    elif not node.prev:  # node is head
        ll.head = node.next
        node.next.prev = NULL
    elif not node.next:  # node is tail
        node.prev.next = NULL
        ll.tail = node.prev
    else:
        node.prev.next = node.next
        node.next.prev = node.prev
    ll.length -= 1

cdef _LinkedList_push_front(_LinkedList * ll, _LinkedListNode * node):
    """
    Add node ``node`` to the beginning of the linked list ``ll``.
    """
    if not ll.head:
        _LinkedList_set_head(ll, node)
    else:
        ll.head.prev = node
        node.next = ll.head
        ll.head = node
        ll.length += 1

cdef _LinkedList_concatenate(_LinkedList * lst1, _LinkedList * lst2):
    """
    Concatenate lst2 to lst1.

    Makes lst2 empty.
    """
    lst1.tail.next = lst2.head
    lst2.head.prev = lst1.tail
    lst1.tail = lst2.tail
    lst1.length += lst2.length
    lst2.head = NULL
    lst2.length = 0

cdef str _LinkedList_to_string(_LinkedList * ll):
    """
    Return a string representation of ``self``.
    """
    cdef _LinkedListNode * temp = ll.head
    cdef list s = []
    while temp:
        s.append(str(temp.data))
        temp = temp.next
    return " ".join(s)

cdef class _Component:
    """
    Connected component class.

    This is a helper class for ``TriconnectivitySPQR``.

    This class is used to store a connected component. It contains:

    - ``edge_list`` -- list of edges belonging to the component,
      stored as a :class:`_LinkedList`

    - ``component_type`` -- the type of the component

      - 0 if bond.
      - 1 if polygon.
      - 2 is triconnected component.
    """
    def __init__(self, list edge_list, int type_c):
        """
        Initialize this component.

        INPUT:

        - ``edge_list`` -- list of edges to be added to the component

        - ``type_c`` -- type of the component (0, 1, or 2)

        TESTS::

            sage: cython_code = [
            ....: 'from sage.graphs.connectivity cimport _Component',
            ....: 'cdef _Component comp = _Component([], 0)',
            ....: 'comp.add_edge(2)',
            ....: 'comp.add_edge(3)',
            ....: 'comp.finish_tric_or_poly(4)',
            ....: 'print(comp)']
            sage: cython(os.linesep.join(cython_code))                                  # needs sage.misc.cython
            Polygon: 2 3 4
        """
        self.mem = MemoryAllocator()
        self.edge_list = <_LinkedList *> self.mem.malloc(sizeof(_LinkedList))
        _LinkedList_initialize(self.edge_list)

        cdef Py_ssize_t e_index
        for e_index in edge_list:
            self.add_edge(e_index)
        self.component_type = type_c

    cdef add_edge(self, Py_ssize_t e_index):
        """
        Add edge index ``e_index`` to the component.
        """
        cdef _LinkedListNode * node = <_LinkedListNode *> self.mem.malloc(sizeof(_LinkedListNode))
        _LinkedListNode_initialize(node, e_index)
        _LinkedList_append(self.edge_list, node)

    cdef finish_tric_or_poly(self, Py_ssize_t e_index):
        r"""
        Finalize the component by adding edge ``e``.

        Edge ``e`` is the last edge to be added to the component.
        Classify the component as a polygon or triconnected component
        depending on the number of edges belonging to it.
        """
        self.add_edge(e_index)
        if _LinkedList_get_length(self.edge_list) > 3:
            self.component_type = 2
        else:
            self.component_type = 1

    def __str__(self):
        """
        Return a string representation of the component.

        TESTS::

            sage: cython_code = [
            ....: 'from sage.graphs.connectivity cimport _Component',
            ....: 'cdef _Component comp = _Component([], 0)',
            ....: 'comp.add_edge(2)',
            ....: 'comp.add_edge(3)',
            ....: 'comp.finish_tric_or_poly(4)',
            ....: 'print(comp)']
            sage: cython(os.linesep.join(cython_code))                                  # needs sage.misc.cython
            Polygon: 2 3 4
        """
        if self.component_type == 0:
            type_str = "Bond: "
        elif self.component_type == 1:
            type_str = "Polygon: "
        else:
            type_str = "Triconnected: "
        return type_str + _LinkedList_to_string(self.edge_list)

    cdef list get_edge_list(self):
        """
        Return the list of edges belonging to the component.
        """
        cdef list e_list = []
        cdef _LinkedListNode * e_node = _LinkedList_get_head(self.edge_list)
        while e_node:
            e_list.append(e_node.data)
            e_node = e_node.next
        return e_list


cdef class TriconnectivitySPQR:
    r"""
    Decompose a graph into triconnected components and build SPQR-tree.

    This class implements the algorithm proposed by Hopcroft and Tarjan in
    [Hopcroft1973]_, and later corrected by Gutwenger and Mutzel in [Gut2001]_,
    for finding the triconnected components of a biconnected graph. It then
    organizes these components into a SPQR-tree. See the:wikipedia:`SPQR_tree`.

    A SPQR-tree is a tree data structure used to represent the triconnected
    components of a biconnected (multi)graph and the 2-vertex cuts separating
    them. A node of a SPQR-tree, and the graph associated with it, can be one of
    the following four types:

    - ``'S'`` -- the associated graph is a cycle with at least three vertices
      ``'S'`` stands for ``series`` and is also called a ``polygon``

    - ``'P'`` -- the associated graph is a dipole graph, a multigraph with two
      vertices and three or more edges. ``'P'`` stands for ``parallel`` and the
      node is called a ``bond``.

    - ``'Q'`` -- the associated graph has a single real edge. This trivial case
      is necessary to handle the graph that has only one edge.

    - ``'R'`` -- the associated graph is a 3-vertex-connected graph that is not
      a cycle or dipole. ``'R'`` stands for ``rigid``.

    The edges of the tree indicate the 2-vertex cuts of the graph.

    INPUT:

    - ``G`` -- graph; if ``G`` is a :class:`DiGraph`, the computation is done on
      the underlying :class:`Graph` (i.e., ignoring edge orientation)

    - ``check`` -- boolean (default: ``True``); indicates whether ``G`` needs to
      be tested for biconnectivity

    .. SEEALSO::

        - :meth:`sage.graphs.connectivity.spqr_tree`
        - :meth:`~Graph.is_biconnected`
        - :wikipedia:`SPQR_tree`

    EXAMPLES:

    Example from the :wikipedia:`SPQR_tree`::

        sage: from sage.graphs.connectivity import TriconnectivitySPQR
        sage: from sage.graphs.connectivity import spqr_tree_to_graph
        sage: G = Graph([(1, 2), (1, 4), (1, 8), (1, 12), (3, 4), (2, 3),
        ....: (2, 13), (3, 13), (4, 5), (4, 7), (5, 6), (5, 8), (5, 7), (6, 7),
        ....: (8, 11), (8, 9), (8, 12), (9, 10), (9, 11), (9, 12), (10, 12)])
        sage: tric = TriconnectivitySPQR(G)
        sage: T = tric.get_spqr_tree()
        sage: G.is_isomorphic(spqr_tree_to_graph(T))
        True

    An example from [Hopcroft1973]_::

        sage: G = Graph([(1, 2), (1, 4), (1, 8), (1, 12), (1, 13), (2, 3),
        ....: (2, 13), (3, 4), (3, 13), (4, 5), (4, 7), (5, 6), (5, 7), (5, 8),
        ....: (6, 7), (8, 9), (8, 11), (8, 12), (9, 10), (9, 11), (9, 12),
        ....: (10, 11), (10, 12)])
        sage: tric = TriconnectivitySPQR(G)
        sage: T = tric.get_spqr_tree()
        sage: G.is_isomorphic(spqr_tree_to_graph(T))
        True
        sage: tric.print_triconnected_components()
        Triconnected: [(8, 9, None), (9, 12, None), (9, 11, None), (8, 11, None), (10, 11, None), (9, 10, None), (10, 12, None), (8, 12, 'newVEdge0')]
        Bond: [(8, 12, None), (8, 12, 'newVEdge0'), (8, 12, 'newVEdge1')]
        Polygon: [(6, 7, None), (5, 6, None), (7, 5, 'newVEdge2')]
        Bond: [(7, 5, 'newVEdge2'), (5, 7, 'newVEdge3'), (5, 7, None)]
        Polygon: [(5, 7, 'newVEdge3'), (4, 7, None), (5, 4, 'newVEdge4')]
        Bond: [(5, 4, 'newVEdge4'), (4, 5, 'newVEdge5'), (4, 5, None)]
        Polygon: [(4, 5, 'newVEdge5'), (5, 8, None), (1, 4, 'newVEdge9'), (1, 8, 'newVEdge10')]
        Triconnected: [(1, 2, None), (2, 13, None), (1, 13, None), (3, 13, None), (2, 3, None), (1, 3, 'newVEdge7')]
        Polygon: [(1, 3, 'newVEdge7'), (3, 4, None), (1, 4, 'newVEdge8')]
        Bond: [(1, 4, None), (1, 4, 'newVEdge8'), (1, 4, 'newVEdge9')]
        Bond: [(1, 8, None), (1, 8, 'newVEdge10'), (1, 8, 'newVEdge11')]
        Polygon: [(8, 12, 'newVEdge1'), (1, 8, 'newVEdge11'), (1, 12, None)]

    An example from [Gut2001]_::

        sage: G = Graph([(1, 2), (1, 4), (2, 3), (2, 5), (3, 4), (3, 5), (4, 5),
        ....: (4, 6), (5, 7), (5, 8), (5, 14), (6, 8), (7, 14), (8, 9), (8, 10),
        ....: (8, 11), (8, 12), (9, 10), (10, 13), (10, 14), (10, 15), (10, 16),
        ....: (11, 12), (11, 13), (12, 13), (14, 15), (14, 16), (15, 16)])
        sage: T = TriconnectivitySPQR(G).get_spqr_tree()
        sage: G.is_isomorphic(spqr_tree_to_graph(T))
        True

    An example with multi-edges and accessing the triconnected components::

        sage: G = Graph([(1, 2), (1, 5), (1, 5), (2, 3), (2, 3), (3, 4), (4, 5)], multiedges=True)
        sage: tric = TriconnectivitySPQR(G)
        sage: T = tric.get_spqr_tree()
        sage: G.is_isomorphic(spqr_tree_to_graph(T))
        True
        sage: tric.print_triconnected_components()
        Bond:  [(1, 5, None), (1, 5, None), (1, 5, 'newVEdge0')]
        Bond:  [(2, 3, None), (2, 3, None), (2, 3, 'newVEdge1')]
        Polygon:  [(4, 5, None), (1, 5, 'newVEdge0'), (3, 4, None), (2, 3, 'newVEdge1'), (1, 2, None)]

    An example of a triconnected graph::

        sage: G = Graph([('a', 'b'), ('a', 'c'), ('a', 'd'), ('b', 'c'), ('b', 'd'), ('c', 'd')])
        sage: T = TriconnectivitySPQR(G).get_spqr_tree()
        sage: print(T.vertices(sort=True))
        [('R', Multi-graph on 4 vertices)]
        sage: G.is_isomorphic(spqr_tree_to_graph(T))
        True

    An example of a directed graph with multi-edges::

        sage: G = DiGraph([(1, 2), (2, 3), (3, 4), (4, 5), (1, 5), (5, 1)])
        sage: tric = TriconnectivitySPQR(G)
        sage: tric.print_triconnected_components()
        Bond:  [(1, 5, None), (5, 1, None), (1, 5, 'newVEdge0')]
        Polygon:  [(4, 5, None), (1, 5, 'newVEdge0'), (3, 4, None), (2, 3, None), (1, 2, None)]

    Edge labels are preserved by the construction::

        sage: G = Graph([(0, 1, '01'), (0, 4, '04'), (1, 2, '12'), (1, 5, '15'),
        ....: (2, 3, '23'), (2, 6, '26'), (3, 7, '37'), (4, 5, '45'),
        ....: (5, 6, '56'), (6, 7, 67)])
        sage: T = TriconnectivitySPQR(G).get_spqr_tree()
        sage: H = spqr_tree_to_graph(T)
        sage: all(G.has_edge(e) for e in H.edge_iterator())
        True
        sage: all(H.has_edge(e) for e in G.edge_iterator())
        True

    TESTS:

    A disconnected graph::

        sage: from sage.graphs.connectivity import TriconnectivitySPQR
        sage: G = Graph([(1,2),(3,5)])
        sage: tric = TriconnectivitySPQR(G)
        Traceback (most recent call last):
        ...
        ValueError: graph is not connected

    A graph with a cut vertex::

        sage: from sage.graphs.connectivity import TriconnectivitySPQR
        sage: G = Graph([(1,2),(1,3),(2,3),(3,4),(3,5),(4,5)])
        sage: tric = TriconnectivitySPQR(G)
        Traceback (most recent call last):
        ...
        ValueError: graph has a cut vertex
    """
    def __init__(self, G, check=True):
        """
        Initialize this object, decompose the graph and build SPQR-tree.

        INPUT:

        - ``G`` -- graph; if ``G`` is a :class:`DiGraph`, the computation is
          done on the underlying :class:`Graph` (i.e., ignoring edge
          orientation)

        - ``check`` -- boolean (default: ``True``); indicates whether ``G``
          needs to be tested for biconnectivity

        EXAMPLES:

        Example from the :wikipedia:`SPQR_tree`::

            sage: from sage.graphs.connectivity import TriconnectivitySPQR
            sage: from sage.graphs.connectivity import spqr_tree_to_graph
            sage: G = Graph([(1, 2), (1, 4), (1, 8), (1, 12), (3, 4), (2, 3),
            ....: (2, 13), (3, 13), (4, 5), (4, 7), (5, 6), (5, 8), (5, 7), (6, 7),
            ....: (8, 11), (8, 9), (8, 12), (9, 10), (9, 11), (9, 12), (10, 12)])
            sage: tric = TriconnectivitySPQR(G)
            sage: T = tric.get_spqr_tree()
            sage: G.is_isomorphic(spqr_tree_to_graph(T))
            True
        """
        self.n = G.order()
        self.m = G.size()
        self.graph_name = G.name()
        self.mem = MemoryAllocator()

        # We set the largest possible index of an edge to 2 * m + 1
        # The algorithm creates at most n virtual edges, so this is large enough
        self.max_number_of_edges = 2 * self.m + 1

        # Trivial cases
        if self.n < 2:
            raise ValueError("graph is not biconnected")
        elif self.n == 2 and self.m:
            # a P block with at least 1 edge
            self.comp_final_edge_list = [G.edges(sort=False)]
            self.comp_type = [0]
            self.__build_spqr_tree()
            return
        elif self.m < self.n - 1:
            # less edges than a tree
            raise ValueError("graph is not connected")
        elif self.m < self.n:
            # less edges than a cycle
            raise ValueError("graph is not biconnected")

        cdef Py_ssize_t i, j
        cdef Py_ssize_t e_index

        # We relabel the graph and store it in different arrays:
        # - Vertices are relabeled as integers in [0..n-1]
        # - Edges are relabeled with distinct labels in [0..m-1] to distinguish
        #   between multi-edges
        # - Virtual edges created by the algorithm have labels >= m
        # - We use these edge labels as unique edge identifiers. Each of these
        #   edge labels is also the index of the edge extremities and original
        #   edge label in appropriate arrays
        # - The status of an edge is: unseen=0, tree=1, frond=2, inactive=-1
        self.int_to_vertex = list(G)
        self.vertex_to_int = {u: i for i, u in enumerate(self.int_to_vertex)}
        self.edge_extremity_first = <int * > self.mem.allocarray(self.max_number_of_edges, sizeof(int))
        self.edge_extremity_second = <int * > self.mem.allocarray(self.max_number_of_edges, sizeof(int))
        self.int_to_original_edge_label = []  # to associate original edge label
        self.edge_status = <int *> self.mem.allocarray(self.max_number_of_edges, sizeof(int))
        for e_index, (u, v, l) in enumerate(G.edge_iterator()):
            self.int_to_original_edge_label.append(l)
            self.edge_extremity_first[e_index] = self.vertex_to_int[u]
            self.edge_extremity_second[e_index] = self.vertex_to_int[v]
            self.edge_status[e_index] = 0

        # Label used for virtual edges, incremented at every new virtual edge
        self.virtual_edge_num = 0

        #
        # Initialize data structures needed for the algorithm
        #

        # Edges of the graph which are in the reverse direction in palm tree
        self.reverse_edges = <bint *> self.mem.allocarray(self.max_number_of_edges, sizeof(bint))
        for i in range(self.max_number_of_edges):
            self.reverse_edges[i] = False

        # DFS number of vertex i
        self.dfs_number = <int *> self.mem.allocarray(self.n, sizeof(int))
        for i in range(self.n):
            self.dfs_number[i] = 0

        # Linked list of fronds entering vertex i in the order they are visited
        self.highpt = <_LinkedList **> self.mem.allocarray(self.n, sizeof(_LinkedList *))
        for i in range(self.n):
            self.highpt[i] = <_LinkedList *> self.mem.malloc(sizeof(_LinkedList))
            _LinkedList_initialize(self.highpt[i])

        # A dictionary whose key is an edge e, value is a pointer to element in
        # self.highpt containing the edge e. Used in the `path_search` function.
        self.in_high = <_LinkedListNode **> self.mem.allocarray(self.max_number_of_edges, sizeof(_LinkedListNode *))
        for i in range(self.max_number_of_edges):
            self.in_high[i] = NULL

        # Translates DFS number of a vertex to its new number
        self.old_to_new = <int *> self.mem.allocarray(self.n + 1, sizeof(int))
        self.newnum = <int *> self.mem.allocarray(self.n + 1, sizeof(int))
        self.node_at = <int *> self.mem.allocarray(self.n + 1, sizeof(int))
        self.lowpt1 = <int *> self.mem.allocarray(self.n + 1, sizeof(int))
        self.lowpt2 = <int *> self.mem.allocarray(self.n + 1, sizeof(int))
        for i in range(self.n + 1):
            self.old_to_new[i] = 0
            self.newnum[i] = 0
            self.node_at[i] = 0
            self.lowpt1[i] = -1
            self.lowpt2[i] = -1

        # i^th value contains a LinkedList of incident edges of vertex i
        self.adj = <_LinkedList **> self.mem.allocarray(self.n, sizeof(_LinkedList *))
        for i in range(self.n):
            self.adj[i] = <_LinkedList *> self.mem.malloc(sizeof(_LinkedList))
            _LinkedList_initialize(self.adj[i])

        # A dictionary whose key is an edge, value is a pointer to element in
        # self.adj containing the edge. Used in the `path_search` function.
        self.in_adj = <_LinkedListNode **> self.mem.allocarray(self.max_number_of_edges, sizeof(_LinkedListNode *))
        for i in range(self.max_number_of_edges):
            self.in_adj[i] = NULL

        self.nd = <int *> self.mem.allocarray(self.n, sizeof(int))
        for i in range(self.n):
            self.nd[i] = 0

        # Parent vertex of vertex i in the palm tree
        self.parent = <int *> self.mem.allocarray(self.n, sizeof(int))
        self.degree = <int *> self.mem.allocarray(self.n, sizeof(int))
        self.tree_arc = <int *> self.mem.allocarray(self.n, sizeof(int))
        self.vertex_at = <int *> self.mem.allocarray(self.n, sizeof(int))
        for i in range(self.n):
            self.parent[i] = -1
            self.degree[i] = 0
            self.tree_arc[i] = -1
            self.vertex_at[i] = 1

        self.dfs_counter = 0
        self.components_list = []  # list of components of `graph_copy`
        self.graph_copy_adjacency = [[] for i in range(self.n)]  # Stores adjacency list

        # Dictionary of (e, True/False) to denote if edge e starts a path
        self.starts_path = <bint *> self.mem.allocarray(self.max_number_of_edges, sizeof(bint))

        # Stacks used in `path_search` function
        self.e_stack = []
        self.t_stack_h = <int *> self.mem.allocarray(self.max_number_of_edges, sizeof(int))
        self.t_stack_a = <int *> self.mem.allocarray(self.max_number_of_edges, sizeof(int))
        self.t_stack_b = <int *> self.mem.allocarray(self.max_number_of_edges, sizeof(int))
        self.t_stack_top = 0
        self.t_stack_a[self.t_stack_top] = -1

        # The final triconnected components are stored
        self.comp_final_edge_list = []  # i^th entry is list of edges in i^th component
        self.comp_type = []  # i^th entry is type of i^th component
        # associate final edge e to its internal index
        self.final_edge_to_edge_index = {}
        # The final SPQR tree is stored
        self.spqr_tree = None  # Graph

        # Arrays used in different methods. We allocate them only once
        self.tmp_array_n_int_1 = <int *> self.mem.allocarray(self.n, sizeof(int))
        self.tmp_array_n_int_2 = <int *> self.mem.allocarray(self.n, sizeof(int))
        self.tmp_array_n_int_3 = <int *> self.mem.allocarray(self.n, sizeof(int))
        self.tmp_array_n_bint_1 = <bint *> self.mem.allocarray(self.n, sizeof(bint))

        #
        # Triconnectivity algorithm
        #

        # Deal with multiple edges
        self.__split_multiple_edges()

        # Build adjacency list
        for e_index in range(self.m + self.virtual_edge_num):
            if self.edge_status[e_index] == -1:
                continue
            i = self.edge_extremity_first[e_index]
            j = self.edge_extremity_second[e_index]
            self.graph_copy_adjacency[i].append(e_index)
            self.graph_copy_adjacency[j].append(e_index)
            self.degree[i] += 1
            self.degree[j] += 1

        self.dfs_counter = 0  # Initialisation for dfs1()
        self.start_vertex = 0  # Initialisation for dfs1()
        cdef int cut_vertex = self.__dfs1(self.start_vertex, check=check)

        if check:
            # If graph is disconnected
            if self.dfs_counter < self.n:
                raise ValueError("graph is not connected")

            # If graph has a cut vertex
            if cut_vertex != -1:
                raise ValueError("graph has a cut vertex")

        # Identify reversed edges to reflect the palm tree arcs and fronds
        cdef bint up
        for e_index in range(self.m + self.virtual_edge_num):
            if self.edge_status[e_index] == -1:
                continue
            i = self.edge_extremity_first[e_index]
            j = self.edge_extremity_second[e_index]
            up = (self.dfs_number[j] - self.dfs_number[i]) > 0
            if (up and self.edge_status[e_index] == 2) or (not up and self.edge_status[e_index] == 1):
                # Add edge to the set reverse_edges
                self.reverse_edges[e_index] = True

        self.__build_acceptable_adj_struct()
        self.__dfs2()

        self.__path_search(self.start_vertex)

        # last split component
        cdef _Component c
        if self.e_stack:
            e_index = self.__estack_pop()
            c = _Component(self.e_stack, 0)
            c.finish_tric_or_poly(e_index)
            self.components_list.append(c)

        self.__assemble_triconnected_components()

        self.__build_spqr_tree()

    cdef int __new_virtual_edge(self, int u, int v) noexcept:
        """
        Return a new virtual edge between ``u`` and ``v``.
        """
        cdef Py_ssize_t e_index = self.m + self.virtual_edge_num
        self.int_to_original_edge_label.append("newVEdge"+str(self.virtual_edge_num))
        self.virtual_edge_num += 1
        self.edge_extremity_first[e_index] = u
        self.edge_extremity_second[e_index] = v
        self.edge_status[e_index] = 0
        return e_index

    cdef _LinkedListNode * __new_LinkedListNode(self, Py_ssize_t e_index) noexcept:
        """
        Create a new ``_LinkedListNode`` initialized with value ``e_index``.
        """
        cdef _LinkedListNode * node = <_LinkedListNode *> self.mem.malloc(sizeof(_LinkedListNode))
        _LinkedListNode_initialize(node, e_index)
        return node

    cdef Py_ssize_t __high(self, Py_ssize_t v) noexcept:
        """
        Return the ``high(v)`` value, which is the first value in
        ``highpt`` list of ``v``.
        """
        cdef _LinkedListNode * head = _LinkedList_get_head(self.highpt[v])
        if head:
            return head.data
        return 0

    cdef __del_high(self, int e_index):
        """
        Delete edge ``e`` from the ``highpt`` list of the endpoint ``v``
        it belongs to.
        """
        cdef int v
        cdef _LinkedListNode * it = self.in_high[e_index]
        if it:
            if self.reverse_edges[e_index]:
                v = self.edge_extremity_first[e_index]
            else:
                v = self.edge_extremity_second[e_index]
            _LinkedList_remove(self.highpt[v], it)

    cdef __split_multiple_edges(self):
        """
        Make the graph simple and build bonds recording multiple edges.

        If there are `k` multiple edges between `u` and `v`, then a new
        component (a bond) with `k+1` edges (one of them is a virtual edge) will
        be created, all the `k` edges are deleted from the graph and the virtual
        edge between `u` and `v` is added to the graph.
        """
        cdef dict sub_bucket
        cdef list b, sb
        cdef int u, v, e_index, virtual_e_index
        cdef list bucket = [[] for u in range(self.n)]

        # We form buckets of edges with same min(e[0], e[1])
        for e_index in range(self.m):
            u = min(self.edge_extremity_first[e_index], self.edge_extremity_second[e_index])
            bucket[u].append(e_index)

        # We split each bucket into sub-buckets with same max(e[0], e[1]) thus
        # identifying groups of multiple edges
        for u, b in enumerate(bucket):
            if not b or len(b) == 1:
                # Nothing to do
                continue
            sub_bucket = {}
            for e_index in b:
                v = self.__edge_other_extremity(e_index, u)
                if v in sub_bucket:
                    sub_bucket[v].append(e_index)
                else:
                    sub_bucket[v] = [e_index]

            for v, sb in sub_bucket.items():
                if len(sb) == 1:
                    continue

                # We have multiple edges. We remove them from graph_copy, add a
                # virtual edge to graph_copy, and create a component containing
                # all removed multiple edges and the virtual edge.
                for e_index in sb:
                    self.edge_status[e_index] = -1

                virtual_e_index = self.__new_virtual_edge(u, v)
                self.edge_status[virtual_e_index] = 0

                sb.append(virtual_e_index)
                self.__new_component(sb, 0)

    cdef int __dfs1(self, int start, bint check=True) noexcept:
        """
        Build the palm-tree of the graph using a dfs traversal.

        Also populates the lists ``lowpt1``, ``lowpt2``, ``nd``, ``parent``,
        and ``dfs_number``.  It updates the dict ``edge_status`` to reflect
        palm tree arcs and fronds.

        INPUT:

        - ``start`` -- the start vertex for DFS

        - ``check`` -- if ``True``, the graph is tested for biconnectivity; if
          the graph has a cut vertex, the cut vertex is returned; otherwise
          the graph is assumed to be biconnected, function returns ``None``

        OUTPUT:

        - If ``check`` is set to ``True`` and a cut vertex is found, the cut
          vertex is returned. If no cut vertex is found, return ``-1``.
        - If ``check`` is set to ``False``, ``-1`` is returned.
        """
        cdef Py_ssize_t v, w
        cdef Py_ssize_t e_index
        cdef int cut_vertex = -1  # Storing the cut vertex, if any
        cdef int* adjacency = self.tmp_array_n_int_3
        cdef list cur_adj
        cdef Py_ssize_t len_cur_adj
        for v in range(self.n):
            adjacency[v] = 0

        # Defining a stack. stack_top == -1 means empty stack
        cdef int* stack = self.tmp_array_n_int_1
        cdef Py_ssize_t stack_top = 0
        stack[stack_top] = start

        # Used for testing biconnectivity
        cdef int* first_son = self.tmp_array_n_int_2
        for v in range(self.n):
            first_son[v] = -1

        while stack_top != -1:
            v = stack[stack_top]

            if not self.dfs_number[v]:
                self.dfs_counter += 1
                self.dfs_number[v] = self.dfs_counter
                self.lowpt1[v] = self.lowpt2[v] = self.dfs_number[v]
                self.nd[v] = 1

            cur_adj = self.graph_copy_adjacency[v]
            len_cur_adj = len(cur_adj)
            # Find the next e_index such that self.edge_status[e_index] is False
            if adjacency[v] == len_cur_adj:
                adjacency[v] = -1
            elif adjacency[v] != -1:
                e_index = cur_adj[adjacency[v]]
                adjacency[v] += 1
                while self.edge_status[e_index] > 0:
                    if adjacency[v] == len_cur_adj:
                        adjacency[v] = -1
                        break
                    e_index = cur_adj[adjacency[v]]
                    adjacency[v] += 1

            if adjacency[v] != -1:
                # Opposite vertex of edge e
                w = self.__edge_other_extremity(e_index, v)
                if not self.dfs_number[w]:
                    self.edge_status[e_index] = 1  # tree edge
                    if first_son[v] == -1:
                        first_son[v] = w
                    self.tree_arc[w] = e_index

                    stack_top += 1
                    stack[stack_top] = w
                    self.parent[w] = v

                else:
                    self.edge_status[e_index] = 2  # frond
                    if self.dfs_number[w] < self.lowpt1[v]:
                        self.lowpt2[v] = self.lowpt1[v]
                        self.lowpt1[v] = self.dfs_number[w]
                    elif self.dfs_number[w] > self.lowpt1[v]:
                        self.lowpt2[v] = min(self.lowpt2[v], self.dfs_number[w])

            else:
                # We trackback, so w takes the value of v and we pop the stack
                w = stack[stack_top]
                stack_top -= 1

                # Test termination
                if stack_top == -1:
                    break

                v = stack[stack_top]

                if check:
                    # Check for cut vertex.
                    # The situation in which there is no path from w to an
                    # ancestor of v : we have identified a cut vertex
                    if ((self.lowpt1[w] >= self.dfs_number[v])
                            and (w != first_son[v] or self.parent[v] != -1)):
                        cut_vertex = v

                # Calculate the `lowpt1` and `lowpt2` values.
                # `lowpt1` is the smallest vertex (the vertex x with smallest
                # dfs_number[x]) that can be reached from v.
                # `lowpt2` is the next smallest vertex that can be reached from v.
                if self.lowpt1[w] < self.lowpt1[v]:
                    self.lowpt2[v] = min(self.lowpt1[v], self.lowpt2[w])
                    self.lowpt1[v] = self.lowpt1[w]

                elif self.lowpt1[w] == self.lowpt1[v]:
                    self.lowpt2[v] = min(self.lowpt2[v], self.lowpt2[w])

                else:
                    self.lowpt2[v] = min(self.lowpt2[v], self.lowpt1[w])

                self.nd[v] += self.nd[w]

        return cut_vertex  # cut_vertex is -1 if graph does not have a cut vertex

    cdef __build_acceptable_adj_struct(self):
        """
        Build the adjacency lists for each vertex with certain properties of
        the ordering, using the ``lowpt1`` and ``lowpt2`` values.

        The list ``adj`` and the dictionary ``in_adj`` are populated.

        ``phi`` values of each edge are calculated using the ``lowpt`` values of
        incident vertices. The edges are then sorted by the ``phi`` values and
        added to adjacency list.
        """
        cdef Py_ssize_t max_size = 3 * self.n + 2
        cdef Py_ssize_t i, u, v
        cdef int e_index, edge_type, phi
        cdef list bucket = [[] for i in range(max_size + 1)]
        cdef _LinkedListNode * node

        for e_index in range(self.m + self.virtual_edge_num):
            edge_type = self.edge_status[e_index]
            if edge_type == -1:
                continue
            u = self.edge_extremity_first[e_index]
            v = self.edge_extremity_second[e_index]

            # Compute phi value
            # bucket sort adjacency list by phi values
            if self.reverse_edges[e_index]:
                if edge_type == 1:  # tree arc
                    if self.lowpt2[u] < self.dfs_number[v]:
                        phi = 3 * self.lowpt1[u]
                    else:
                        phi = 3 * self.lowpt1[u] + 2
                else:  # tree frond
                    phi = 3 * self.dfs_number[u] + 1
            else:
                if edge_type == 1:  # tree arc
                    if self.lowpt2[v] < self.dfs_number[u]:
                        phi = 3 * self.lowpt1[v]
                    else:
                        phi = 3 * self.lowpt1[v] + 2
                else:  # tree frond
                    phi = 3 * self.dfs_number[v] + 1

            bucket[phi].append(e_index)

        # Populate `adj` and `in_adj` with the sorted edges
        for i in range(1, max_size + 1):
            for e_index in bucket[i]:
                node = self.__new_LinkedListNode(e_index)
                if self.reverse_edges[e_index]:
                    _LinkedList_append(self.adj[self.edge_extremity_second[e_index]], node)
                else:
                    _LinkedList_append(self.adj[self.edge_extremity_first[e_index]], node)
                self.in_adj[e_index] = node

    cdef __path_finder(self, int start):
        """
        This function is a helper function for :meth:`__dfs2` function.

        Calculate ``newnum[v]`` and identify the edges which start a new path.

        INPUT:

        - ``start`` -- the start vertex
        """
        cdef bint new_path = True
        cdef Py_ssize_t v, w
        cdef Py_ssize_t e_index
        cdef _LinkedListNode * e_node
        cdef _LinkedListNode * highpt_node

        # Defining a stack. stack_top == -1 means empty stack
        cdef int* stack = self.tmp_array_n_int_1
        cdef Py_ssize_t stack_top = 0
        stack[stack_top] = start

        cdef bint * seen = self.tmp_array_n_bint_1
        for v in range(self.n):
            seen[v] = False

        cdef _LinkedListNode ** pointer_e_node = <_LinkedListNode ** > self.mem.allocarray(self.n, sizeof(_LinkedListNode *))
        for v in range(self.n):
            pointer_e_node[v] = _LinkedList_get_head(self.adj[v])

        while stack_top != -1:
            v = stack[stack_top]
            if not seen[v]:
                self.newnum[v] = self.dfs_counter - self.nd[v] + 1
                seen[v] = True
            e_node = pointer_e_node[v]

            if e_node:
                e_index = e_node.data
                pointer_e_node[v] = e_node.next
                # opposite vertex of e
                w = self.__edge_other_extremity(e_index, v)
                if new_path:
                    new_path = False
                    self.starts_path[e_index] = True
                if self.edge_status[e_index] == 1:  # tree arc
                    stack_top += 1
                    stack[stack_top] = w
                else:
                    # Identified a new frond that enters `w`. Add to `highpt[w]`.
                    highpt_node = self.__new_LinkedListNode(self.newnum[v])
                    _LinkedList_append(self.highpt[w], highpt_node)
                    self.in_high[e_index] = highpt_node
                    new_path = True

            else:
                # We trackback
                self.dfs_counter -= 1
                stack_top -= 1

    cdef __dfs2(self):
        """
        Update the values of ``lowpt1`` and ``lowpt2`` lists with the
        help of new numbering obtained from :meth:`__path_finder`.
        Populate ``highpt`` values.
        """
        cdef Py_ssize_t v
        cdef Py_ssize_t e_index

        self.dfs_counter = self.n
        for e_index in range(self.m + self.virtual_edge_num):
            self.in_high[e_index] = NULL
            self.starts_path[e_index] = False

        # We call the pathFinder function with the start vertex
        self.__path_finder(self.start_vertex)

        # Update `old_to_new` values with the calculated `newnum` values
        for v in range(self.n):
            self.old_to_new[self.dfs_number[v]] = self.newnum[v]

        # Update lowpt values according to `newnum` values.
        for v in range(self.n):
            self.node_at[self.newnum[v]] = v
            self.lowpt1[v] = self.old_to_new[self.lowpt1[v]]
            self.lowpt2[v] = self.old_to_new[self.lowpt2[v]]

    cdef int __path_search(self, int start) except -1:
        """
        Find the separation pairs and construct the split components.

        Check for type-1 and type-2 separation pairs, and construct the split
        components while also creating new virtual edges wherever required.

        INPUT:

        - ``start`` -- the start vertex
        """
        cdef int e_index, e_virt_index
        cdef int x, y, h, xx
        cdef int v, vnum, outv
        cdef int w, wnum
        cdef int temp_index, temp_target
        cdef int a, b, e_ab_index, e_ab_source
        cdef int e1_index, e2_index, e2_source
        cdef int xy_index, xy_target
        cdef int eh_index, eh_source
        cdef _LinkedListNode * it
        cdef _LinkedListNode * e_node
        cdef _LinkedListNode * temp_node
        cdef _LinkedListNode * vnum_node
        cdef _LinkedListNode * e_virt_node
        cdef _Component comp

        # Defining a stack. stack_v_top == -1 means empty stack
        cdef int* stack_v = self.tmp_array_n_int_1
        cdef Py_ssize_t stack_v_top = 0
        stack_v[stack_v_top] = start

        cdef int* y_dict = self.tmp_array_n_int_2
        y_dict[start] = 0

        cdef int* outv_dict = self.tmp_array_n_int_3
        outv_dict[start] = _LinkedList_get_length(self.adj[start])

        cdef _LinkedListNode ** e_node_dict = <_LinkedListNode **> self.mem.allocarray(self.n, sizeof(_LinkedListNode *))
        e_node_dict[start] = _LinkedList_get_head(self.adj[start])

        while stack_v_top != -1:
            v = stack_v[stack_v_top]
            e_node = e_node_dict[v]

            if e_node:
                # Restore values of variables
                y = y_dict[v]
                vnum = self.newnum[v]
                outv = outv_dict[v]
                e_index = e_node.data
                it = e_node
                if self.reverse_edges[e_index]:
                    w = self.edge_extremity_first[e_index]  # target
                else:
                    w = self.edge_extremity_second[e_index]
                wnum = self.newnum[w]

                if self.edge_status[e_index] == 1:  # e is a tree arc
                    if self.starts_path[e_index]:  # if a new path starts at edge e
                        # Pop all (h,a,b) from tstack where a > lowpt1[w]
                        if self.t_stack_a[self.t_stack_top] > self.lowpt1[w]:
                            while self.t_stack_a[self.t_stack_top] > self.lowpt1[w]:
                                y = max(y, self.t_stack_h[self.t_stack_top])
                                b = self.t_stack_b[self.t_stack_top]
                                self.t_stack_top -= 1
                            self.__tstack_push(y, self.lowpt1[w], b)

                        else:
                            self.__tstack_push(wnum + self.nd[w] - 1, self.lowpt1[w], vnum)
                        self.__tstack_push_eos()

                    # We emulate the recursive call on w using a stack
                    stack_v_top += 1
                    stack_v[stack_v_top] = w
                    y_dict[w] = 0
                    outv_dict[w] = _LinkedList_get_length(self.adj[w])
                    e_node_dict[w] = _LinkedList_get_head(self.adj[w])
                    y_dict[v] = y
                    continue

                else:  # e is a frond
                    if self.starts_path[e_index]:
                        # pop all (h,a,b) from tstack where a > w
                        if self.t_stack_a[self.t_stack_top] > wnum:
                            while self.t_stack_a[self.t_stack_top] > wnum:
                                y = max(y, self.t_stack_h[self.t_stack_top])
                                b = self.t_stack_b[self.t_stack_top]
                                self.t_stack_top -= 1
                            self.__tstack_push(y, wnum, b)

                        else:
                            self.__tstack_push(vnum, wnum, vnum)
                    self.e_stack.append(e_index)  # add edge (v,w) to ESTACK

            else:
                # We are done with v, so we trackback
                stack_v_top -= 1

                # Test termination
                if stack_v_top == -1:
                    continue

                # Restore state of variables
                v = stack_v[stack_v_top]
                e_node = e_node_dict[v]
                y = y_dict[v]
                vnum = self.newnum[v]
                outv = outv_dict[v]
                e_index = e_node.data
                it = e_node
                if self.reverse_edges[e_index]:
                    w = self.edge_extremity_first[e_index]  # target
                else:
                    w = self.edge_extremity_second[e_index]
                wnum = self.newnum[w]

                # Continue operations with tree arc e

                self.e_stack.append(self.tree_arc[w])
                temp_node = _LinkedList_get_head(self.adj[w])
                temp_index = temp_node.data
                if self.reverse_edges[temp_index]:
                    temp_target = self.edge_extremity_first[temp_index]
                else:
                    temp_target = self.edge_extremity_second[temp_index]

                # Type-2 separation pair check
                # while v is not the start_vertex
                while vnum != 1 and ((self.t_stack_a[self.t_stack_top] == vnum)
                                     or (self.degree[w] == 2 and self.newnum[temp_target] > wnum)):
                    a = self.t_stack_a[self.t_stack_top]
                    b = self.t_stack_b[self.t_stack_top]
                    if a == vnum and self.parent[self.node_at[b]] == self.node_at[a]:
                        self.t_stack_top -= 1

                    else:
                        e_ab_index = -1
                        if self.degree[w] == 2 and self.newnum[temp_target] > wnum:
                            # found type-2 separation pair - (v, temp_target)
                            e1_index = self.__estack_pop()
                            e2_index = self.__estack_pop()
                            _LinkedList_remove(self.adj[w], self.in_adj[e2_index])

                            if self.reverse_edges[e2_index]:
                                x = self.edge_extremity_first[e2_index]  # target
                            else:
                                x = self.edge_extremity_second[e2_index]  # target

                            e_virt_index = self.__new_virtual_edge(v, x)
                            self.degree[v] -= 1
                            self.degree[x] -= 1

                            if self.reverse_edges[e2_index]:
                                e2_source = self.edge_extremity_second[e2_index]  # target
                            else:
                                e2_source = self.edge_extremity_first[e2_index]
                            if e2_source != w:
                                raise ValueError("graph is not biconnected")

                            self.__new_component([e1_index, e2_index, e_virt_index], 1)

                            if self.e_stack:
                                e1_index = self.e_stack[-1]
                                if self.reverse_edges[e1_index]:
                                    if (self.edge_extremity_first[e1_index] == v
                                            and self.edge_extremity_second[e1_index] == x):
                                        e_ab_index = self.__estack_pop()
                                        _LinkedList_remove(self.adj[x], self.in_adj[e_ab_index])
                                        self.__del_high(e_ab_index)
                                else:
                                    if (self.edge_extremity_first[e1_index] == x
                                            and self.edge_extremity_second[e1_index] == v):
                                        e_ab_index = self.__estack_pop()
                                        _LinkedList_remove(self.adj[x], self.in_adj[e_ab_index])
                                        self.__del_high(e_ab_index)

                        else:  # found type-2 separation pair - (self.node_at[a], self.node_at[b])
                            h = self.t_stack_h[self.t_stack_top]
                            self.t_stack_top -= 1

                            comp = _Component([], 0)
                            while True:
                                xy_index = self.e_stack[-1]
                                if self.reverse_edges[xy_index]:
                                    x = self.edge_extremity_second[xy_index]
                                    xy_target = self.edge_extremity_first[xy_index]
                                else:
                                    x = self.edge_extremity_first[xy_index]
                                    xy_target = self.edge_extremity_second[xy_index]
                                if not (a <= self.newnum[x] and self.newnum[x] <= h
                                        and a <= self.newnum[xy_target] and self.newnum[xy_target] <= h):
                                    break
                                if ((self.newnum[x] == a and self.newnum[xy_target] == b)
                                        or (self.newnum[xy_target] == a and self.newnum[x] == b)):
                                    e_ab_index = self.__estack_pop()
                                    if self.reverse_edges[e_ab_index]:
                                        e_ab_source = self.edge_extremity_second[e_ab_index]  # source
                                    else:
                                        e_ab_source = self.edge_extremity_first[e_ab_index]  # source
                                    _LinkedList_remove(self.adj[e_ab_source], self.in_adj[e_ab_index])
                                    self.__del_high(e_ab_index)

                                else:
                                    eh_index = self.__estack_pop()
                                    if self.reverse_edges[eh_index]:
                                        eh_source = self.edge_extremity_second[eh_index]
                                    else:
                                        eh_source = self.edge_extremity_first[eh_index]
                                    if it != self.in_adj[eh_index]:
                                        _LinkedList_remove(self.adj[eh_source], self.in_adj[eh_index])
                                        self.__del_high(eh_index)

                                    comp.add_edge(eh_index)
                                    self.degree[x] -= 1
                                    self.degree[xy_target] -= 1

                            e_virt_index = self.__new_virtual_edge(self.node_at[a], self.node_at[b])
                            comp.finish_tric_or_poly(e_virt_index)
                            self.components_list.append(comp)
                            comp = None
                            x = self.node_at[b]

                        if e_ab_index != -1:
                            comp = _Component([e_ab_index, e_virt_index], 0)
                            e_virt_index = self.__new_virtual_edge(v, x)
                            comp.add_edge(e_virt_index)
                            self.degree[x] -= 1
                            self.degree[v] -= 1
                            self.components_list.append(comp)
                            comp = None

                        self.e_stack.append(e_virt_index)
                        # Replace the edge in `it` with `e_virt`
                        it.data = e_virt_index

                        self.in_adj[e_virt_index] = it
                        self.degree[x] += 1
                        self.degree[v] += 1
                        self.parent[x] = v
                        self.tree_arc[x] = e_virt_index
                        self.edge_status[e_virt_index] = 1
                        w = x
                        wnum = self.newnum[w]

                    # update the values used in the while loop check
                    temp_node = _LinkedList_get_head(self.adj[w])
                    temp_index = temp_node.data
                    if self.reverse_edges[temp_index]:
                        temp_target = self.edge_extremity_first[temp_index]
                    else:
                        temp_target = self.edge_extremity_second[temp_index]

                # start type-1 check
                if (self.lowpt2[w] >= vnum and self.lowpt1[w] < vnum
                        and (self.parent[v] != self.start_vertex or outv >= 2)):
                    # type-1 separation pair - (self.node_at[self.lowpt1[w]], v)
                    # Create a new component and add edges to it
                    comp = _Component([], 0)
                    if not self.e_stack:
                        raise ValueError("stack is empty")
                    while self.e_stack:
                        xy_index = self.e_stack[-1]
                        if self.reverse_edges[xy_index]:
                            xx = self.newnum[self.edge_extremity_second[xy_index]]  # source
                            y = self.newnum[self.edge_extremity_first[xy_index]]  # target
                        else:
                            xx = self.newnum[self.edge_extremity_first[xy_index]]  # source
                            y = self.newnum[self.edge_extremity_second[xy_index]]  # target

                        if not ((wnum <= xx and xx < wnum + self.nd[w])
                                or (wnum <= y and y < wnum + self.nd[w])):
                            break

                        comp.add_edge(self.__estack_pop())
                        self.__del_high(xy_index)
                        self.degree[self.node_at[xx]] -= 1
                        self.degree[self.node_at[y]] -= 1

                    e_virt_index = self.__new_virtual_edge(v, self.node_at[self.lowpt1[w]])
                    comp.finish_tric_or_poly(e_virt_index)  # Add virtual edge to component
                    self.components_list.append(comp)
                    comp = None

                    if ((xx == vnum and y == self.lowpt1[w])
                            or (y == vnum and xx == self.lowpt1[w])):
                        comp_bond = _Component([], 0)  # new triple bond
                        eh_index = self.__estack_pop()
                        if self.in_adj[eh_index] != it:
                            if self.reverse_edges[eh_index]:
                                _LinkedList_remove(self.adj[self.edge_extremity_second[eh_index]], self.in_adj[eh_index])
                            else:
                                _LinkedList_remove(self.adj[self.edge_extremity_first[eh_index]], self.in_adj[eh_index])

                        comp_bond.add_edge(eh_index)
                        comp_bond.add_edge(e_virt_index)
                        e_virt_index = self.__new_virtual_edge(v, self.node_at[self.lowpt1[w]])
                        comp_bond.add_edge(e_virt_index)
                        if self.in_high[eh_index]:
                            self.in_high[e_virt_index] = self.in_high[eh_index]
                        self.degree[v] -= 1
                        self.degree[self.node_at[self.lowpt1[w]]] -= 1

                        self.components_list.append(comp_bond)
                        comp_bond = None

                    if self.node_at[self.lowpt1[w]] != self.parent[v]:
                        self.e_stack.append(e_virt_index)

                        # replace edge in `it` with `e_virt`
                        it.data = e_virt_index

                        self.in_adj[e_virt_index] = it
                        if not self.in_high[e_virt_index] and self.__high(self.node_at[self.lowpt1[w]]) < vnum:
                            vnum_node = self.__new_LinkedListNode(vnum)
                            _LinkedList_push_front(self.highpt[self.node_at[self.lowpt1[w]]], vnum_node)
                            self.in_high[e_virt_index] = vnum_node

                        self.degree[v] += 1
                        self.degree[self.node_at[self.lowpt1[w]]] += 1

                    else:
                        _LinkedList_remove(self.adj[v], it)
                        comp_bond = _Component([e_virt_index], 0)
                        e_virt_index = self.__new_virtual_edge(self.node_at[self.lowpt1[w]], v)
                        comp_bond.add_edge(e_virt_index)

                        eh_index = self.tree_arc[v]
                        comp_bond.add_edge(eh_index)

                        self.components_list.append(comp_bond)
                        comp_bond = None

                        self.tree_arc[v] = e_virt_index
                        self.edge_status[e_virt_index] = 1
                        if self.in_adj[eh_index]:
                            self.in_adj[e_virt_index] = self.in_adj[eh_index]
                        e_virt_node = self.__new_LinkedListNode(e_virt_index)
                        self.in_adj[eh_index] = e_virt_node
                        # end type-1 search

                # if a path starts at edge e, empty the tstack.
                if self.starts_path[e_index]:
                    while self.__tstack_not_eos():
                        self.t_stack_top -= 1
                    self.t_stack_top -= 1

                while (self.__tstack_not_eos() and self.t_stack_b[self.t_stack_top] != vnum
                       and self.__high(v) > self.t_stack_h[self.t_stack_top]):
                    self.t_stack_top -= 1

                outv_dict[v] -= 1

            # Go to next edge in adjacency list
            e_node_dict[v] = e_node.next

    cdef __assemble_triconnected_components(self):
        """
        Iterate through all the split components built by :meth:`__path_finder`
        and merges two bonds or two polygons that share an edge for constructing
        the final triconnected components.

        Subsequently, convert the edges in triconnected components into original
        vertices and edges. The triconnected components are stored in
        ``self.comp_final_edge_list`` and ``self.comp_type``.
        """
        cdef Py_ssize_t i, j
        cdef Py_ssize_t e_index
        cdef _Component c1, c2
        cdef int c1_type
        cdef _LinkedListNode * e_node
        cdef _LinkedListNode * e_node_next
        cdef _LinkedList * l1
        cdef _LinkedList * l2

        cdef Py_ssize_t num_components = len(self.components_list)
        cdef bint* visited = <bint*> self.mem.allocarray(num_components, sizeof(bint))
        for i in range(num_components):
            visited[i] = False

        # The index of first (second) component that an edge belongs to
        cdef int* comp1 = <int*> self.mem.allocarray(self.m + self.virtual_edge_num, sizeof(int))
        cdef int* comp2 = <int*> self.mem.allocarray(self.m + self.virtual_edge_num, sizeof(int))

        # Pointer to the edge node in first (second) component
        cdef _LinkedListNode ** item1 = <_LinkedListNode **> self.mem.allocarray(self.m + self.virtual_edge_num, sizeof(_LinkedListNode *))
        cdef _LinkedListNode ** item2 = <_LinkedListNode **> self.mem.allocarray(self.m + self.virtual_edge_num, sizeof(_LinkedListNode *))
        for i in range(self.m + self.virtual_edge_num):
            item1[i] = NULL
            item2[i] = NULL

        # For each edge, we populate the comp1, comp2, item1 and item2 values
        for i in range(num_components):  # for each component
            e_node = _LinkedList_get_head((<_Component> self.components_list[i]).edge_list)
            while e_node:  # for each edge
                e_index = e_node.data
                if not item1[e_index]:
                    comp1[e_index] = i
                    item1[e_index] = e_node
                else:
                    comp2[e_index] = i
                    item2[e_index] = e_node

                e_node = e_node.next

        # For each edge in a component, if the edge is a virtual edge, merge
        # the two components the edge belongs to
        for i in range(num_components):
            c1 = <_Component> self.components_list[i]
            c1_type = c1.component_type
            l1 = c1.edge_list
            visited[i] = True

            if not _LinkedList_get_length(l1):
                continue

            if c1_type == 0 or c1_type == 1:
                e_node = _LinkedList_get_head((<_Component> self.components_list[i]).edge_list)
                # Iterate through each edge in the component
                while e_node:
                    e_index = e_node.data
                    e_node_next = e_node.next
                    if not self.__is_virtual_edge(e_index):
                        e_node = e_node_next
                        continue

                    j = comp1[e_index]
                    if visited[j]:
                        j = comp2[e_index]
                        if visited[j]:
                            e_node = e_node_next
                            continue
                        e_node2 = item2[e_index]
                    else:
                        e_node2 = item1[e_index]

                    c2 = <_Component> self.components_list[j]

                    # If the two components are not the same type, do not merge
                    if c1_type != c2.component_type:
                        e_node = e_node_next  # Go to next edge
                        continue

                    visited[j] = True
                    l2 = c2.edge_list

                    # Remove the corresponding virtual edges in both the components
                    # and merge the components
                    _LinkedList_remove(l2, e_node2)
                    _LinkedList_concatenate(l1, l2)

                    # if `e_node_next` was empty, after merging two components,
                    # more edges are added to the component.
                    if not e_node_next:
                        e_node_next = e_node.next  # Go to next edge

                    _LinkedList_remove(l1, e_node)

                    e_node = e_node_next

        # Convert connected components into original graph vertices and edges
        cdef list e_list_new
        cdef list e_index_list
        cdef tuple e_new
        self.comp_final_edge_list = []
        self.comp_type = []
        self.final_edge_to_edge_index = {}
        for comp in self.components_list:
            if _LinkedList_get_length((<_Component> comp).edge_list):
                e_index_list = (<_Component> comp).get_edge_list()
                e_list_new = []
                # For each edge, get the original source, target and label
                for e_index in e_index_list:
                    source = self.int_to_vertex[self.edge_extremity_first[e_index]]
                    target = self.int_to_vertex[self.edge_extremity_second[e_index]]
                    label = self.int_to_original_edge_label[e_index]
                    e_new = (source, target, label)
                    e_list_new.append(e_new)
                    self.final_edge_to_edge_index[e_new] = e_index
                # Add the component data to `comp_final_edge_list` and `comp_type`
                self.comp_type.append((<_Component> comp).component_type)
                self.comp_final_edge_list.append(e_list_new)

    cdef __build_spqr_tree(self):
        """
        Build the SPQR-tree of the graph and store it in variable
        ``self.spqr_tree``. See
        :meth:`~sage.graphs.connectivity.TriconnectivitySPQR.get_spqr_tree`.
        """
        # Types of components 0: "P", 1: "S", 2: "R"
        cdef list component_type = ["P", "S", "R"]

        from sage.graphs.graph import Graph
        self.spqr_tree = Graph(multiedges=False, name='SPQR-tree of {}'.format(self.graph_name))

        if len(self.comp_final_edge_list) == 1 and self.comp_type[0] == 0:
            self.spqr_tree.add_vertex(('Q' if len(self.comp_final_edge_list[0]) == 1 else 'P',
                                       Graph(self.comp_final_edge_list[0], immutable=True, multiedges=True)))
            return

        cdef list int_to_vertex = []
        cdef dict partner_nodes = {}
        cdef Py_ssize_t i, j
        cdef Py_ssize_t e_index

        for i in range(len(self.comp_final_edge_list)):
            # Create a new tree vertex
            u = (component_type[self.comp_type[i]],
                 Graph(self.comp_final_edge_list[i], immutable=True, multiedges=True))
            self.spqr_tree.add_vertex(u)
            int_to_vertex.append(u)

            # Add an edge to each node containing the same virtual edge
            for e in self.comp_final_edge_list[i]:
                e_index = self.final_edge_to_edge_index[e]
                if self.__is_virtual_edge(e_index):
                    if e_index in partner_nodes:
                        for j in partner_nodes[e_index]:
                            self.spqr_tree.add_edge(int_to_vertex[i], int_to_vertex[j])
                        partner_nodes[e_index].append(i)
                    else:
                        partner_nodes[e_index] = [i]

    def print_triconnected_components(self):
        """
        Print the type and list of edges of each component.

        EXAMPLES:

        An example from [Hopcroft1973]_::

            sage: from sage.graphs.connectivity import TriconnectivitySPQR
            sage: from sage.graphs.connectivity import spqr_tree_to_graph
            sage: G = Graph([(1, 2), (1, 4), (1, 8), (1, 12), (1, 13), (2, 3),
            ....: (2, 13), (3, 4), (3, 13), (4, 5), (4, 7), (5, 6), (5, 7), (5, 8),
            ....: (6, 7), (8, 9), (8, 11), (8, 12), (9, 10), (9, 11), (9, 12),
            ....: (10, 11), (10, 12)])
            sage: tric = TriconnectivitySPQR(G)
            sage: T = tric.get_spqr_tree()
            sage: G.is_isomorphic(spqr_tree_to_graph(T))
            True
            sage: tric.print_triconnected_components()
            Triconnected: [(8, 9, None), (9, 12, None), (9, 11, None), (8, 11, None), (10, 11, None), (9, 10, None), (10, 12, None), (8, 12, 'newVEdge0')]
            Bond: [(8, 12, None), (8, 12, 'newVEdge0'), (8, 12, 'newVEdge1')]
            Polygon: [(6, 7, None), (5, 6, None), (7, 5, 'newVEdge2')]
            Bond: [(7, 5, 'newVEdge2'), (5, 7, 'newVEdge3'), (5, 7, None)]
            Polygon: [(5, 7, 'newVEdge3'), (4, 7, None), (5, 4, 'newVEdge4')]
            Bond: [(5, 4, 'newVEdge4'), (4, 5, 'newVEdge5'), (4, 5, None)]
            Polygon: [(4, 5, 'newVEdge5'), (5, 8, None), (1, 4, 'newVEdge9'), (1, 8, 'newVEdge10')]
            Triconnected: [(1, 2, None), (2, 13, None), (1, 13, None), (3, 13, None), (2, 3, None), (1, 3, 'newVEdge7')]
            Polygon: [(1, 3, 'newVEdge7'), (3, 4, None), (1, 4, 'newVEdge8')]
            Bond: [(1, 4, None), (1, 4, 'newVEdge8'), (1, 4, 'newVEdge9')]
            Bond: [(1, 8, None), (1, 8, 'newVEdge10'), (1, 8, 'newVEdge11')]
            Polygon: [(8, 12, 'newVEdge1'), (1, 8, 'newVEdge11'), (1, 12, None)]
        """
        # The types are {0: "Bond", 1: "Polygon", 2: "Triconnected"}
        cdef list prefix = ["Bond", "Polygon", "Triconnected"]
        cdef Py_ssize_t i
        for i in range(len(self.comp_final_edge_list)):
            print("{}: {}".format(prefix[self.comp_type[i]], self.comp_final_edge_list[i]))

    def get_triconnected_components(self):
        r"""
        Return the triconnected components as a list of tuples.

        Each component is represented as a tuple of the type of the component
        and the list of edges of the component.

        EXAMPLES::

            sage: from sage.graphs.connectivity import TriconnectivitySPQR
            sage: G = Graph(2)
            sage: for i in range(3):
            ....:     G.add_path([0, G.add_vertex(), G.add_vertex(), 1])
            sage: tric = TriconnectivitySPQR(G)
            sage: tric.get_triconnected_components()
            [('Polygon', [(4, 5, None), (0, 4, None), (1, 5, None), (1, 0, 'newVEdge1')]),
            ('Polygon', [(6, 7, None), (0, 6, None), (1, 7, None), (1, 0, 'newVEdge3')]),
            ('Bond', [(1, 0, 'newVEdge1'), (1, 0, 'newVEdge3'), (1, 0, 'newVEdge4')]),
            ('Polygon', [(1, 3, None), (1, 0, 'newVEdge4'), (2, 3, None), (0, 2, None)])]
        """
        cdef list comps = []
        cdef Py_ssize_t i
        # The types are {0: "Bond", 1: "Polygon", 2: "Triconnected"}
        cdef list prefix = ["Bond", "Polygon", "Triconnected"]
        for i in range(len(self.comp_final_edge_list)):
            comps.append((prefix[self.comp_type[i]], self.comp_final_edge_list[i]))
        return comps

    def get_spqr_tree(self):
        r"""
        Return an SPQR-tree representing the triconnected components of the
        graph.

        An SPQR-tree is a tree data structure used to represent the triconnected
        components of a biconnected (multi)graph and the 2-vertex cuts
        separating them. A node of a SPQR-tree, and the graph associated with
        it, can be one of the following four types:

        - ``'S'`` -- the associated graph is a cycle with at least three vertices.
          ``'S'`` stands for ``series``.

        - ``'P'`` -- the associated graph is a dipole graph, a multigraph with
          two vertices and three or more edges. ``'P'`` stands for ``parallel``.

        - ``'Q'`` -- the associated graph has a single real edge. This trivial
          case is necessary to handle the graph that has only one edge.

        - ``'R'`` -- the associated graph is a 3-connected graph that is not a
          cycle or dipole. ``'R'`` stands for ``rigid``.

        The edges of the tree indicate the 2-vertex cuts of the graph.

        OUTPUT:

        ``SPQR-tree`` a tree whose vertices are labeled with the block's
        type and the subgraph of three-blocks in the decomposition.

        EXAMPLES::

            sage: from sage.graphs.connectivity import TriconnectivitySPQR
            sage: G = Graph(2)
            sage: for i in range(3):
            ....:     G.add_clique([0, 1, G.add_vertex(), G.add_vertex()])
            sage: tric = TriconnectivitySPQR(G)
            sage: Tree = tric.get_spqr_tree()
            sage: K4 = graphs.CompleteGraph(4)
            sage: all(u[1].is_isomorphic(K4) for u in Tree if u[0] == 'R')
            True
            sage: from sage.graphs.connectivity import spqr_tree_to_graph
            sage: G.is_isomorphic(spqr_tree_to_graph(Tree))
            True

            sage: G = Graph(2)
            sage: for i in range(3):
            ....:     G.add_path([0, G.add_vertex(), G.add_vertex(), 1])
            sage: tric = TriconnectivitySPQR(G)
            sage: Tree = tric.get_spqr_tree()
            sage: C4 = graphs.CycleGraph(4)
            sage: all(u[1].is_isomorphic(C4) for u in Tree if u[0] == 'S')
            True
            sage: G.is_isomorphic(spqr_tree_to_graph(Tree))
            True

            sage: G.allow_multiple_edges(True)
            sage: G.add_edges(G.edge_iterator())
            sage: tric = TriconnectivitySPQR(G)
            sage: Tree = tric.get_spqr_tree()
            sage: all(u[1].is_isomorphic(C4) for u in Tree if u[0] == 'S')
            True
            sage: G.is_isomorphic(spqr_tree_to_graph(Tree))
            True

            sage: G = graphs.CycleGraph(6)
            sage: tric = TriconnectivitySPQR(G)
            sage: Tree = tric.get_spqr_tree()
            sage: Tree.order()
            1
            sage: G.is_isomorphic(spqr_tree_to_graph(Tree))
            True
            sage: G.add_edge(0, 3)
            sage: tric = TriconnectivitySPQR(G)
            sage: Tree = tric.get_spqr_tree()
            sage: Tree.order()
            3
            sage: G.is_isomorphic(spqr_tree_to_graph(Tree))
            True

            sage: G = Graph([(0, 1)], multiedges=True)
            sage: tric = TriconnectivitySPQR(G)
            sage: Tree = tric.get_spqr_tree()
            sage: Tree.vertices(sort=True)
            [('Q', Multi-graph on 2 vertices)]
            sage: G.add_edge(0, 1)
            sage: Tree = TriconnectivitySPQR(G).get_spqr_tree()
            sage: Tree.vertices(sort=True)
            [('P', Multi-graph on 2 vertices)]
        """
        return self.spqr_tree


def is_triconnected(G):
    r"""
    Check whether the graph is triconnected.

    A triconnected graph is a connected graph on 3 or more vertices that is not
    broken into disconnected pieces by deleting any pair of vertices.

    EXAMPLES:

    The Petersen graph is triconnected::

        sage: G = graphs.PetersenGraph()
        sage: G.is_triconnected()
        True

    But a 2D grid is not::

        sage: G = graphs.Grid2dGraph(3, 3)
        sage: G.is_triconnected()
        False

    By convention, a cycle of order 3 is triconnected::

        sage: G = graphs.CycleGraph(3)
        sage: G.is_triconnected()
        True

    But cycles of order 4 and more are not::

        sage: [graphs.CycleGraph(i).is_triconnected() for i in range(4, 8)]
        [False, False, False, False]

    Comparing different methods on random graphs that are not always
    triconnected::

        sage: G = graphs.RandomBarabasiAlbert(50, 3)                                    # needs networkx
        sage: G.is_triconnected() == G.vertex_connectivity(k=3)                         # needs networkx
        True

    .. SEEALSO::

        - :meth:`~sage.graphs.generic_graph.GenericGraph.is_connected`
        - :meth:`~Graph.is_biconnected`
        - :meth:`~sage.graphs.connectivity.spqr_tree`
        - :wikipedia:`SPQR_tree`

    TESTS::

        sage: [Graph(i).is_triconnected() for i in range(4)]
        [False, False, False, False]
        sage: [graphs.CompleteGraph(i).is_triconnected() for i in range(3, 6)]
        [True, True, True]
    """
    if G.order() < 3:
        return False

    try:
        T = G.spqr_tree()
    except ValueError:
        # The graph is not biconnected
        return False

    from collections import Counter
    C = Counter(v[0] for v in T)
    if 'S' in C:
        return G.order() == 3
    # Since the graph has order >= 3, is biconnected and has no 'S' block, it
    # has at least one 'R' block. A triconnected graph has only one such block.
    return C['R'] == 1
