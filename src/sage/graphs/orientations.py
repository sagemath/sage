r"""
Orientations

This module implements several methods to compute orientations of undirected
graphs subject to specific constraints (e.g., acyclic, strongly connected,
etc.). It also implements some iterators over all these orientations.

**This module contains the following methods**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`strong_orientations_iterator` | Return an iterator over all strong orientations of a graph `G`
    :meth:`random_orientation` | Return a random orientation of a graph `G`


Authors
-------

- Kolja Knauer, Petru Valicov (2017-01-10) -- initial version


Methods
-------
"""
from copy import copy
from sage.graphs.digraph import DiGraph


def strong_orientations_iterator(G):
    r"""
    Returns an iterator over all strong orientations of a graph `G`.

    A strong orientation of a graph is an orientation of its edges such that
    the obtained digraph is strongly connected (i.e. there exist a directed path
    between each pair of vertices).

    ALGORITHM:

    It is an adaptation of the algorithm published in [CGMRV16]_.
    It runs in `O(mn)` amortized time, where `m` is the number of edges and
    `n` is the number of vertices. The amortized time can be improved to `O(m)`
    with a more involved method.
    In this function, first the graph is preprocessed and a spanning tree is
    generated. Then every orientation of the non-tree edges of the graph can be
    extended to at least one new strong orientation by orienting properly
    the edges of the spanning tree (this property is proved in [CGMRV16]_).
    Therefore, this function generates all partial orientations of the non-tree
    edges and then launches a helper function corresponding to the generation
    algorithm described in [CGMRV16]_.
    In order to avoid trivial symmetries, the orientation of an arbitrary edge
    is fixed before the start of the enumeration process.

    INPUT:

    - ``G`` -- an undirected graph.

    OUTPUT:

    - an iterator which will produce all strong orientations of this graph.

    .. NOTE::

        Works only for simple graphs (no multiple edges).
        To avoid symmetries an orientation of an arbitrary edge is fixed.


    EXAMPLES:

    A cycle has one possible (non-symmetric) strong orientation::

        sage: g = graphs.CycleGraph(4)
        sage: it = g.strong_orientations_iterator()
        sage: len(list(it))
        1

    A tree cannot be strongly oriented::

        sage: g = graphs.RandomTree(100)
        sage: len(list(g.strong_orientations_iterator()))
        0

    Neither can be a disconnected graph::

        sage: g = graphs.CompleteGraph(6)
        sage: g.add_vertex(7)
        sage: len(list(g.strong_orientations_iterator()))
        0

    TESTS:

        sage: g = graphs.CompleteGraph(2)
        sage: len(list(g.strong_orientations_iterator()))
        0

        sage: g = graphs.CubeGraph(3)
        sage: b = True
        sage: for orientedGraph in g.strong_orientations_iterator():
        ....:     if not orientedGraph.is_strongly_connected():
        ....:         b = False
        sage: b
        True

    The total number of strong orientations of a graph can be counted using
    the Tutte polynomial evaluated at points (0,2)::

        sage: g = graphs.PetersenGraph()
        sage: nr1 = len(list(g.strong_orientations_iterator()))
        sage: nr2 = g.tutte_polynomial()(0,2)
        sage: nr1 == nr2/2 # The Tutte polynomial counts also the symmetrical orientations
        True

    """
    # if the graph has a bridge or is disconnected,
    # then it cannot be strongly oriented
    if G.order() < 3 or not G.is_biconnected():
        return

    V = G.vertices()
    Dg = DiGraph([G.vertices(), G.edges()], pos=G.get_pos())

    # compute an arbitrary spanning tree of the undirected graph
    te = G.min_spanning_tree()
    treeEdges = [(u,v) for u,v,_ in te]
    tree_edges_set = set(treeEdges)
    A = [edge for edge in G.edge_iterator(labels=False) if edge not in tree_edges_set]

    # initialization of the first binary word 00...0
    # corresponding to the current orientation of the non-tree edges
    existingAedges = [0] * len(A)

    # Make the edges of the spanning tree doubly oriented
    for e in treeEdges:
        if Dg.has_edge(e):
            Dg.add_edge(e[1], e[0])
        else:
            Dg.add_edge(e)

    # Generate all orientations for non-tree edges (using Gray code)
    # Each of these orientations can be extended to a strong orientation
    # of G by orienting properly the tree-edges
    previousWord = 0
    i = 0

    # the orientation of one edge is fixed so we consider one edge less
    nr = 2**(len(A)-1)
    while i < nr:
        word = (i >> 1) ^ i
        bitChanged = word ^ previousWord
        
        bit = 0
        while bitChanged > 1:
            bitChanged >>= 1
            bit += 1

        previousWord = word
        if not existingAedges[bit]:
            Dg.reverse_edge(A[bit])
            existingAedges[bit] = 1
        else:
            Dg.reverse_edge(A[bit][1], A[bit][0])
            existingAedges[bit] = 0
        # launch the algorithm for enumeration of the solutions
        for sol in _strong_orientations_of_a_mixed_graph(Dg, V, treeEdges):
            yield sol
        i = i + 1


def _strong_orientations_of_a_mixed_graph(Dg, V, E):
    r"""
    Helper function for the generation of all strong orientations.

    Generates all strong orientations of a given partially directed graph
    (also called mixed graph). The algorithm finds bound edges i.e undirected
    edges whose orientation is forced and tries all possible orientations for
    the other edges. See [CGMRV16]_ for more details.

    INPUT:

    - ``Dg`` -- the mixed graph. The undirected edges are doubly oriented.

    - ``V`` -- the set of vertices

    - ``E`` -- the set of undirected edges (they are oriented in both ways);
      No labels are allowed.

    OUTPUT:

    - an iterator which will produce all strong orientations of the input
      partially directed graph.

    EXAMPLES::

        sage: from sage.graphs.orientations import _strong_orientations_of_a_mixed_graph
        sage: g = graphs.CycleGraph(5)
        sage: Dg = DiGraph(g) # all edges of g will be doubly oriented
        sage: it = _strong_orientations_of_a_mixed_graph(Dg, list(g), list(g.edges(labels=False, sort=False)))
        sage: len(list(it)) # there are two orientations of this multigraph
        2
    """
    length = len(E)
    i = 0
    boundEdges = []
    while i < length:
        u, v = E[i]
        Dg.delete_edge(u, v)
        if not (v in Dg.depth_first_search(u)):
            del E[i]
            length -= 1
            Dg.add_edge(u, v)
            Dg.delete_edge(v, u)
            boundEdges.append((v, u))
        else:
            Dg.add_edge(u, v)
            Dg.delete_edge(v, u)
            if not (u in Dg.depth_first_search(v)):
                del E[i]
                length -= 1
                boundEdges.append((u, v))
                Dg.delete_edge(u, v)
            else:
                i += 1
            Dg.add_edge(v, u)

    # if true the obtained orientation is strong
    if not E:
        yield Dg.copy()
    else:
        u, v = E.pop()
        Dg.delete_edge(v, u)
        for orientation in _strong_orientations_of_a_mixed_graph(Dg, V, E):
            yield orientation
        Dg.add_edge(v, u)
        Dg.delete_edge(u, v)
        for orientation in _strong_orientations_of_a_mixed_graph(Dg, V, E):
            yield orientation
        Dg.add_edge(u, v)
        E.append((u, v))
    Dg.add_edges(boundEdges)
    E.extend(boundEdges)


def random_orientation(G):
    r"""
    Return a random orientation of a graph `G`.

    An *orientation* of an undirected graph is a directed graph such that every
    edge is assigned a direction. Hence there are `2^m` oriented digraphs for a
    simple graph with `m` edges.

    INPUT:

    - ``G`` -- a Graph.

    EXAMPLES::

        sage: from sage.graphs.orientations import random_orientation
        sage: G = graphs.PetersenGraph()
        sage: D = random_orientation(G)
        sage: D.order() == G.order(), D.size() == G.size()
        (True, True)

    TESTS:

    Giving anything else than a Graph::

        sage: random_orientation([])
        Traceback (most recent call last):
        ...
        ValueError: the input parameter must be a Graph

    .. SEEALSO::

        - :meth:`~Graph.orientations`
    """
    from sage.graphs.graph import Graph
    if not isinstance(G, Graph):
        raise ValueError("the input parameter must be a Graph")

    D = DiGraph(data=[G.vertices(), []],
                format='vertices_and_edges',
                multiedges=G.allows_multiple_edges(),
                loops=G.allows_loops(),
                weighted=G.weighted(),
                pos=G.get_pos(),
                name="Random orientation of {}".format(G.name()))
    if hasattr(G, '_embedding'):
        D._embedding = copy(G._embedding)

    from sage.misc.prandom import getrandbits
    rbits = getrandbits(G.size())
    for u,v,l in G.edge_iterator():
        if rbits % 2:
            D.add_edge(u, v, l)
        else:
            D.add_edge(v, u, l)
        rbits >>= 1
    return D
