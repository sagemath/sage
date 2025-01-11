r"""
Basic graphs

The methods defined here appear in :mod:`sage.graphs.graph_generators`.
"""
# ****************************************************************************
#           Copyright (C) 2006 Robert L. Miller <rlmillster@gmail.com>
#                              and Emily A. Kirkman
#                         2009 Michael C. Yurko <myurko@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

# import from Sage library
from sage.graphs.graph import Graph
from math import sin, cos, pi


def BullGraph(immutable=False):
    r"""
    Return a bull graph with 5 nodes.

    A bull graph is named for its shape. It's a triangle with horns.
    See the :wikipedia:`Bull_graph` for more information.

    INPUT:

    - ``immutable`` -- boolean (default: ``False``); whether to return an
      immutable or a mutable graph

    PLOTTING:

    Upon construction, the position dictionary is filled to override the
    spring-layout algorithm. By convention, the bull graph is drawn as a
    triangle with the first node (0) on the bottom. The second and third nodes
    (1 and 2) complete the triangle. Node 3 is the horn connected to 1 and node
    4 is the horn connected to node 2.

    EXAMPLES:

    Construct and show a bull graph::

        sage: g = graphs.BullGraph(); g
        Bull graph: Graph on 5 vertices
        sage: g.show()                          # long time                             # needs sage.plot

    The bull graph has 5 vertices and 5 edges. Its radius is 2, its
    diameter 3, and its girth 3. The bull graph is planar with chromatic
    number 3 and chromatic index also 3::

        sage: g.order(); g.size()
        5
        5
        sage: g.radius(); g.diameter(); g.girth()
        2
        3
        3
        sage: g.chromatic_number()
        3

    The bull graph has chromatic polynomial `x(x - 2)(x - 1)^3` and
    Tutte polynomial `x^4 + x^3 + x^2 y`. Its characteristic polynomial
    is `x(x^2 - x - 3)(x^2 + x - 1)`, which follows from the definition of
    characteristic polynomials for graphs, i.e. `\det(xI - A)`, where
    `x` is a variable, `A` the adjacency matrix of the graph, and `I`
    the identity matrix of the same dimensions as `A`::

        sage: # needs sage.libs.flint
        sage: chrompoly = g.chromatic_polynomial()
        sage: x = chrompoly.parent()('x')
        sage: x * (x - 2) * (x - 1)^3 == chrompoly
        True

        sage: # needs sage.libs.flint sage.modules
        sage: charpoly = g.characteristic_polynomial()
        sage: M = g.adjacency_matrix(); M
        [0 1 1 0 0]
        [1 0 1 1 0]
        [1 1 0 0 1]
        [0 1 0 0 0]
        [0 0 1 0 0]
        sage: Id = identity_matrix(ZZ, M.nrows())
        sage: D = x*Id - M
        sage: D.determinant() == charpoly                                               # needs sage.symbolic
        True
        sage: x * (x^2 - x - 3) * (x^2 + x - 1) == charpoly
        True

    TESTS:

    Check the behavior of parameter ``immutable``::

        sage: graphs.BullGraph(immutable=True).is_immutable()
        True
    """
    edge_list = [(0, 1), (0, 2), (1, 2), (1, 3), (2, 4)]
    pos_dict = {0: (0, 0), 1: (-1, 1), 2: (1, 1), 3: (-2, 2), 4: (2, 2)}
    return Graph([range(5), edge_list], format='vertices_and_edges',
                 immutable=immutable, pos=pos_dict, name="Bull graph")


def ButterflyGraph(immutable=False):
    r"""
    Return the butterfly graph.

    Let `C_3` be the cycle graph on 3 vertices. The butterfly or bowtie
    graph is obtained by joining two copies of `C_3` at a common vertex,
    resulting in a graph that is isomorphic to the friendship graph `F_2`.
    See the :wikipedia:`Butterfly_graph` for more information.

    .. SEEALSO::

        - :meth:`~sage.graphs.graph_generators.GraphGenerators.FriendshipGraph`

    INPUT:

    - ``immutable`` -- boolean (default: ``False``); whether to return an
      immutable or a mutable graph

    EXAMPLES:

    The butterfly graph is a planar graph on 5 vertices and having 6 edges::

        sage: G = graphs.ButterflyGraph(); G
        Butterfly graph: Graph on 5 vertices
        sage: G.show()                          # long time                             # needs sage.plot
        sage: G.is_planar()
        True
        sage: G.order()
        5
        sage: G.size()
        6

    It has diameter 2, girth 3, and radius 1::

        sage: G.diameter()
        2
        sage: G.girth()
        3
        sage: G.radius()
        1

    The butterfly graph is Eulerian, with chromatic number 3::

        sage: G.is_eulerian()
        True
        sage: G.chromatic_number()
        3

    TESTS:

    Check the behavior of parameter ``immutable``::

        sage: graphs.ButterflyGraph(immutable=True).is_immutable()
        True
    """
    edge_dict = {
        0: [3, 4],
        1: [2, 4],
        2: [4],
        3: [4]}
    pos_dict = {
        0: [-1, 1],
        1: [1, 1],
        2: [1, -1],
        3: [-1, -1],
        4: [0, 0]}
    return Graph(edge_dict, format='dict_of_lists',
                 immutable=immutable, pos=pos_dict, name="Butterfly graph")


def CircularLadderGraph(n, immutable=False):
    r"""
    Return a circular ladder graph with `2 * n` nodes.

    A Circular ladder graph is a ladder graph that is connected at the ends,
    i.e.: a ladder bent around so that top meets bottom. Thus it can be
    described as two parallel cycle graphs connected at each corresponding node
    pair.

    PLOTTING: Upon construction, the position dictionary is filled to override
    the spring-layout algorithm. By convention, the circular ladder graph is
    displayed as an inner and outer cycle pair, with the first `n` nodes drawn
    on the inner circle. The first (0) node is drawn at the top of the
    inner-circle, moving clockwise after that. The outer circle is drawn with
    the `(n+1)`-th node at the top, then counterclockwise as well.
    When `n == 2`, we rotate the outer circle by an angle of `\pi/8` to ensure
    that all edges are visible (otherwise the 4 vertices of the graph would be
    placed on a single line).

    INPUT:

    - ``n`` -- nonnegative integer

    - ``immutable`` -- boolean (default: ``False``); whether to return an
      immutable or a mutable graph

    EXAMPLES:

    Construct and show a circular ladder graph with 26 nodes::

        sage: g = graphs.CircularLadderGraph(13)
        sage: g.show()                          # long time                             # needs sage.plot

    Create several circular ladder graphs in a Sage graphics array::

        sage: g = []
        sage: j = []
        sage: for i in range(9):
        ....:    k = graphs.CircularLadderGraph(i+3)
        ....:    g.append(k)
        sage: for i in range(3):                                                        # needs sage.plot
        ....:    n = []
        ....:    for m in range(3):
        ....:        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ....:    j.append(n)
        sage: G = graphics_array(j)                                                     # needs sage.plot
        sage: G.show()                          # long time                             # needs sage.plot

    TESTS::

        sage: G = graphs.CircularLadderGraph(4, immutable=True)
        sage: G.is_immutable()
        True
    """
    from itertools import chain
    edges_1 = zip(range(n), chain(range(1, n), (0,)))
    edges_2 = zip(range(n, 2 * n), chain(range(n + 1, 2 * n), (n,)))
    edges_3 = ((i, i + n) for i in range(n))
    G = Graph([range(2 * n), chain(edges_1, edges_2, edges_3)],
              format='vertices_and_edges', immutable=immutable,
              name="Circular Ladder graph")
    G._circle_embedding(list(range(n)), radius=1, angle=pi/2)
    if n == 2:
        G._circle_embedding(list(range(4)), radius=1, angle=pi/2 + pi/8)
    else:
        G._circle_embedding(list(range(n, 2*n)), radius=2, angle=pi/2)
    return G


def ClawGraph(immutable=False):
    """
    Return a claw graph.

    A claw graph is named for its shape. It is actually a complete
    bipartite graph with ``(n1, n2) = (1, 3)``.

    PLOTTING: See :meth:`CompleteBipartiteGraph`.

    INPUT:

    - ``immutable`` -- boolean (default: ``False``); whether to return an
      immutable or a mutable graph

    EXAMPLES:

    Show a Claw graph::

        sage: (graphs.ClawGraph()).show()       # long time                             # needs sage.plot

    Inspect a Claw graph::

        sage: G = graphs.ClawGraph()
        sage: G
        Claw graph: Graph on 4 vertices

    TESTS:

    Check the behavior of parameter ``immutable``::

        sage: graphs.ClawGraph(immutable=True).is_immutable()
        True
    """
    edge_list = [(0, 1), (0, 2), (0, 3)]
    pos_dict = {0: (0, 1), 1: (-1, 0), 2: (0, 0), 3: (1, 0)}
    return Graph([range(4), edge_list], format='vertices_and_edges',
                 immutable=immutable, pos=pos_dict, name="Claw graph")


def CycleGraph(n, immutable=False):
    r"""
    Return a cycle graph with `n` nodes.

    A cycle graph is a basic structure which is also typically called an
    `n`-gon.

    INPUT:

    - ``n`` -- nonnegative integer; the number of vertices of the cycle graph

    - ``immutable`` -- boolean (default: ``False``); whether to return an
      immutable or a mutable graph

    PLOTTING: Upon construction, the position dictionary is filled to override
    the spring-layout algorithm. By convention, each cycle graph will be
    displayed with the first (0) node at the top, with the rest following in a
    counterclockwise manner.

    The cycle graph is a good opportunity to compare efficiency of filling a
    position dictionary vs. using the spring-layout algorithm for
    plotting. Because the cycle graph is very symmetric, the resulting plots
    should be similar (in cases of small `n`).

    Filling the position dictionary in advance adds `O(n)` to the constructor.

    EXAMPLES:

    Compare plotting using the predefined layout and networkx::

        sage: # needs networkx sage.plot
        sage: import networkx
        sage: n = networkx.cycle_graph(23)
        sage: spring23 = Graph(n)
        sage: posdict23 = graphs.CycleGraph(23)
        sage: spring23.show()                   # long time
        sage: posdict23.show()                  # long time

    We next view many cycle graphs as a Sage graphics array. First we use the
    ``CycleGraph`` constructor, which fills in the position dictionary::

        sage: # needs networkx sage.plot
        sage: g = []
        sage: j = []
        sage: for i in range(9):
        ....:     k = graphs.CycleGraph(i+3)
        ....:     g.append(k)
        sage: for i in range(3):
        ....:     n = []
        ....:     for m in range(3):
        ....:         n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ....:     j.append(n)
        sage: G = graphics_array(j)
        sage: G.show()                          # long time

    Compare to plotting with the spring-layout algorithm::

        sage: # needs networkx sage.plot
        sage: g = []
        sage: j = []
        sage: for i in range(9):
        ....:     spr = networkx.cycle_graph(i+3)
        ....:     k = Graph(spr)
        ....:     g.append(k)
        sage: for i in range(3):
        ....:     n = []
        ....:     for m in range(3):
        ....:         n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ....:     j.append(n)
        sage: G = graphics_array(j)
        sage: G.show()                          # long time

    TESTS:

    The input parameter must be a positive integer::

        sage: G = graphs.CycleGraph(-1)
        Traceback (most recent call last):
        ...
        ValueError: parameter n must be a positive integer

    Check the behavior of parameter ``immutable``::

        sage: graphs.CycleGraph(4, immutable=True).is_immutable()
        True
    """
    if n < 0:
        raise ValueError("parameter n must be a positive integer")

    from itertools import chain
    edges = zip(range(n), chain(range(1, n), (0,))) if n > 1 else []
    G = Graph([range(n), edges], format='vertices_and_edges',
              immutable=immutable, name="Cycle graph")
    if n == 1:
        G.set_pos({0: (0, 0)})
    else:
        G._circle_embedding(list(range(n)), angle=pi/2)
    return G


def CompleteGraph(n, immutable=False):
    r"""
    Return a complete graph on `n` nodes.

    A Complete Graph is a graph in which all nodes are connected to all
    other nodes.

    INPUT:

    - ``n`` -- nonnegative integer; the number of vertices of the complete graph

    - ``immutable`` -- boolean (default: ``False``); whether to return an
      immutable or a mutable graph

    PLOTTING: Upon construction, the position dictionary is filled to
    override the spring-layout algorithm. By convention, each complete
    graph will be displayed with the first (0) node at the top, with
    the rest following in a counterclockwise manner.

    In the complete graph, there is a big difference visually in using
    the spring-layout algorithm vs. the position dictionary used in
    this constructor. The position dictionary flattens the graph,
    making it clear which nodes an edge is connected to. But the
    complete graph offers a good example of how the spring-layout
    works. The edges push outward (everything is connected), causing
    the graph to appear as a 3-dimensional pointy ball. (See examples
    below).

    EXAMPLES:

    We view many Complete graphs with a Sage Graphics Array, first with this
    constructor (i.e., the position dictionary filled)::

        sage: # needs sage.plot
        sage: g = []
        sage: j = []
        sage: for i in range(9):
        ....:     k = graphs.CompleteGraph(i+3)
        ....:     g.append(k)
        sage: for i in range(3):
        ....:     n = []
        ....:     for m in range(3):
        ....:         n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ....:     j.append(n)
        sage: G = graphics_array(j)
        sage: G.show()                          # long time

    We compare to plotting with the spring-layout algorithm::

        sage: # needs networkx sage.plot
        sage: import networkx
        sage: g = []
        sage: j = []
        sage: for i in range(9):
        ....:     spr = networkx.complete_graph(i+3)
        ....:     k = Graph(spr)
        ....:     g.append(k)
        sage: for i in range(3):
        ....:     n = []
        ....:     for m in range(3):
        ....:         n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ....:     j.append(n)
        sage: G = graphics_array(j)
        sage: G.show()                          # long time

    Compare the constructors (results will vary)::

        sage: # needs networkx
        sage: import networkx
        sage: t = cputime()
        sage: n = networkx.complete_graph(389); spring389 = Graph(n)
        sage: cputime(t)  # random
        0.59203700000000126
        sage: t = cputime()
        sage: posdict389 = graphs.CompleteGraph(389)
        sage: cputime(t)  # random
        0.6680419999999998

    We compare plotting::

        sage: # needs networkx
        sage: import networkx
        sage: n = networkx.complete_graph(23)
        sage: spring23 = Graph(n)
        sage: posdict23 = graphs.CompleteGraph(23)
        sage: spring23.show()                   # long time                             # needs sage.plot
        sage: posdict23.show()                  # long time                             # needs sage.plot

    TESTS:

    Check the behavior of parameter ``immutable``::

        sage: graphs.CompleteGraph(4, immutable=True).is_immutable()
        True
    """
    from itertools import combinations
    G = Graph([range(n), combinations(range(n), 2)],
              format='vertices_and_edges', immutable=immutable,
              name="Complete graph")
    if n == 1:
        G.set_pos({0: (0, 0)})
    else:
        G._circle_embedding(list(range(n)), angle=pi/2)
    return G


def CorrelationGraph(seqs, alpha, include_anticorrelation, immutable=False):
    r"""
    Return a correlation graph with a node per sequence in ``seqs``.

    Edges are added between nodes where the corresponding sequences have a
    correlation coefficient greater than alpha.

    If ``include_anticorrelation`` is ``True``, then edges are also added
    between nodes with correlation coefficient less than ``-alpha``.

    INPUT:

    - ``seqs`` -- list of sequences, that is a list of lists

    - ``alpha`` -- float; threshold on the correlation coefficient between two
      sequences for adding an edge

    - ``include_anticorrelation`` -- boolean; whether to add edges between nodes
      with correlation coefficient less than ``-alpha`` or not

    - ``immutable`` -- boolean (default: ``False``); whether to return an
      immutable or a mutable graph

    EXAMPLES::

        sage: # needs numpy
        sage: from sage.graphs.generators.basic import CorrelationGraph
        sage: data = [[1,2,3], [4,5,6], [7,8,9999]]
        sage: CG1 = CorrelationGraph(data, 0.9, False)
        sage: CG2 = CorrelationGraph(data, 0.9, True)
        sage: CG3 = CorrelationGraph(data, 0.1, True)
        sage: CG1.edges(sort=False)
        [(0, 0, None), (0, 1, None), (1, 1, None), (2, 2, None)]
        sage: CG2.edges(sort=False)
        [(0, 0, None), (0, 1, None), (1, 1, None), (2, 2, None)]
        sage: CG3.edges(sort=False)
        [(0, 0, None), (0, 1, None), (0, 2, None), (1, 1, None), (1, 2, None), (2, 2, None)]

    TESTS:

    Check the behavior of parameter ``immutable``::

        sage: # needs numpy
        sage: from sage.graphs.generators.basic import CorrelationGraph
        sage: CorrelationGraph(data, 0.9, False, immutable=True).is_immutable()
        True
        sage: CorrelationGraph(data, 0.9, True, immutable=True).is_immutable()
        True
    """
    from numpy import corrcoef
    from sage.matrix.constructor import Matrix

    # compute pairwise correlation coefficients
    corrs = corrcoef(seqs)

    # compare against alpha to get adjacency matrix
    if include_anticorrelation:
        boolean_adjacency_matrix = abs(corrs) >= alpha
    else:
        boolean_adjacency_matrix = corrs >= alpha

    adjacency_matrix = Matrix(boolean_adjacency_matrix.astype(int))

    # call graph constructor
    return Graph(adjacency_matrix, format='adjacency_matrix',
                 immutable=immutable, name="Correlation Graph")


def CompleteBipartiteGraph(p, q, set_position=True, immutable=False):
    r"""
    Return a Complete Bipartite Graph on `p + q` vertices.

    A Complete Bipartite Graph is a graph with its vertices partitioned into two
    groups, `V_1 = \{0,...,p-1\}` and `V_2 = \{p,...,p+q-1\}`. Each `u \in
    V_1` is connected to every `v \in V_2`.

    INPUT:

    - ``p``, ``q`` -- number of vertices in each side

    - ``set_position`` -- boolean (default: ``True``); if set to ``True``, we
      assign positions to the vertices so that the set of cardinality `p` is
      on the line `y=1` and the set of cardinality `q` is on the line `y=0`.

    - ``immutable`` -- boolean (default: ``False``); whether to return an
      immutable or a mutable graph

    PLOTTING: Upon construction, the position dictionary is filled to override
    the spring-layout algorithm. By convention, each complete bipartite graph
    will be displayed with the first `p` nodes on the top row (at `y=1`) from
    left to right. The remaining `q` nodes appear at `y=0`, also from left to
    right. The shorter row (partition with fewer nodes) is stretched to the same
    length as the longer row, unless the shorter row has 1 node; in which case
    it is centered. The `x` values in the plot are in domain `[0, \max(p, q)]`.

    In the Complete Bipartite graph, there is a visual difference in using the
    spring-layout algorithm vs. the position dictionary used in this
    constructor. The position dictionary flattens the graph and separates the
    partitioned nodes, making it clear which nodes an edge is connected to. The
    Complete Bipartite graph plotted with the spring-layout algorithm tends to
    center the nodes in `p` (see ``spring_med`` in examples below), thus
    overlapping its nodes and edges, making it typically hard to decipher.

    Filling the position dictionary in advance adds `O(n)` to the constructor.
    Feel free to race the constructors below in the examples section. The much
    larger difference is the time added by the spring-layout algorithm when
    plotting. (Also shown in the example below). The spring model is typically
    described as `O(n^3)`, as appears to be the case in the NetworkX source
    code.

    EXAMPLES:

    Two ways of constructing the complete bipartite graph, using different
    layout algorithms::

        sage: # needs networkx
        sage: import networkx
        sage: n = networkx.complete_bipartite_graph(389, 157)   # long time
        sage: spring_big = Graph(n)             # long time
        sage: posdict_big = graphs.CompleteBipartiteGraph(389, 157)        # long time

    Compare the plotting::

        sage: n = networkx.complete_bipartite_graph(11, 17)                             # needs networkx
        sage: spring_med = Graph(n)                                                     # needs networkx
        sage: posdict_med = graphs.CompleteBipartiteGraph(11, 17)

    Notice here how the spring-layout tends to center the nodes of `n1`::

        sage: spring_med.show()                 # long time                             # needs networkx
        sage: posdict_med.show()                # long time                             # needs sage.plot

    View many complete bipartite graphs with a Sage Graphics Array, with this
    constructor (i.e., the position dictionary filled)::

        sage: g = []
        sage: j = []
        sage: for i in range(9):
        ....:     k = graphs.CompleteBipartiteGraph(i+1,4)
        ....:     g.append(k)
        sage: for i in range(3):                                                        # needs sage.plot
        ....:     n = []
        ....:     for m in range(3):
        ....:         n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ....:     j.append(n)
        sage: G = graphics_array(j)                                                     # needs sage.plot
        sage: G.show()                          # long time                             # needs sage.plot

    We compare to plotting with the spring-layout algorithm::

        sage: # needs networkx sage.plot
        sage: g = []
        sage: j = []
        sage: for i in range(9):
        ....:     spr = networkx.complete_bipartite_graph(i+1,4)
        ....:     k = Graph(spr)
        ....:     g.append(k)
        sage: for i in range(3):
        ....:     n = []
        ....:     for m in range(3):
        ....:         n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ....:     j.append(n)
        sage: G = graphics_array(j)
        sage: G.show()                          # long time

    :issue:`12155`::

        sage: graphs.CompleteBipartiteGraph(5,6).complement()
        complement(Complete bipartite graph of order 5+6): Graph on 11 vertices

    TESTS:

    Prevent negative dimensions (:issue:`18530`)::

        sage: graphs.CompleteBipartiteGraph(-1,1)
        Traceback (most recent call last):
        ...
        ValueError: the arguments p(=-1) and q(=1) must be positive integers
        sage: graphs.CompleteBipartiteGraph(1,-1)
        Traceback (most recent call last):
        ...
        ValueError: the arguments p(=1) and q(=-1) must be positive integers

    Check the behavior of parameter ``immutable``::

        sage: graphs.CompleteBipartiteGraph(1, 2, immutable=True).is_immutable()
        True
    """
    if p < 0 or q < 0:
        raise ValueError('the arguments p(={}) and q(={}) must be positive integers'.format(p, q))

    edges = ((i, j) for i in range(p) for j in range(p, p + q))
    G = Graph([range(p + q), edges], format='vertices_and_edges',
              immutable=immutable,
              name=f"Complete bipartite graph of order {p}+{q}")

    # We now assign positions to vertices:
    # - vertices 0,..,p-1 are placed on the line (0, 1) to (max(p, q), 1)
    # - vertices p,..,p+q-1 are placed on the line (0, 0) to (max(p, q), 0)
    # If p (or q) is 1, the vertex is centered in the line.
    if set_position:
        nmax = max(p, q)
        G._line_embedding(list(range(p)), first=(0, 1), last=(nmax, 1))
        G._line_embedding(list(range(p, p + q)), first=(0, 0), last=(nmax, 0))

    return G


def CompleteMultipartiteGraph(L, immutable=False):
    r"""
    Return a complete multipartite graph.

    INPUT:

    - ``L`` -- list of integers; the respective sizes of the components

    - ``immutable`` -- boolean (default: ``False``); whether to return an
      immutable or a mutable graph

    PLOTTING: Produce a layout of the vertices so that vertices in the same
    vertex set are adjacent and clearly separated from vertices in other vertex
    sets.

    This is done by calculating the vertices of an `r`-gon then calculating the
    slope between adjacent vertices. We then 'walk' around the `r`-gon placing
    graph vertices in regular intervals between adjacent vertices of the
    `r`-gon.

    Makes a nicely organized graph like in this picture:
    https://commons.wikimedia.org/wiki/File:Turan_13-4.svg

    EXAMPLES:

    A complete tripartite graph with sets of sizes `5, 6, 8`::

        sage: g = graphs.CompleteMultipartiteGraph([5, 6, 8]); g
        Multipartite Graph with set sizes [5, 6, 8]: Graph on 19 vertices

    It clearly has a chromatic number of 3::

        sage: g.chromatic_number()
        3

    TESTS:

    Prevent negative dimensions::

        sage: graphs.CompleteMultipartiteGraph([1, -1, 2])
        Traceback (most recent call last):
        ...
        ValueError: the sizes of the components must be positive integers

    Check the bahavior of parameter ``immutable``::

        sage: graphs.CompleteMultipartiteGraph([1], immutable=True).is_immutable()
        True
        sage: graphs.CompleteMultipartiteGraph([1, 2], immutable=True).is_immutable()
        True
        sage: graphs.CompleteMultipartiteGraph([1, 2, 3], immutable=True).is_immutable()
        True
    """
    if any(p < 0 for p in L):
        raise ValueError("the sizes of the components must be positive integers")

    r = len(L)  # getting the number of partitions
    name = "Multipartite Graph with set sizes {}".format(L)

    if not r:
        g = Graph(name=name, immutable=immutable)
    elif r == 1:
        g = Graph(L[0], immutable=immutable, name=name)
        g._line_embedding(range(L[0]), first=(0, 0), last=(L[0], 0))
    elif r == 2:
        p, q = L
        edges = ((i, j) for i in range(p) for j in range(p, p + q))
        g = Graph([range(p + q), edges], format='vertices_and_edges',
                  immutable=immutable,
                  name=name)
        nmax = max(p, q)
        g._line_embedding(list(range(p)), first=(0, 1), last=(nmax, 1))
        g._line_embedding(list(range(p, p + q)), first=(0, 0), last=(nmax, 0))
    else:
        # This position code gives bad results on bipartite or isolated graphs
        points = [(cos(2 * pi * i / r), sin(2 * pi * i / r)) for i in range(r)]
        slopes = [(points[(i + 1) % r][0] - points[i % r][0],
                   points[(i + 1) % r][1] - points[i % r][1]) for i in range(r)]

        counter = 0
        positions = {}
        for i in range(r):
            vertex_set_size = L[i] + 1
            for j in range(1, vertex_set_size):
                x = points[i][0] + slopes[i][0] * j / vertex_set_size
                y = points[i][1] + slopes[i][1] * j / vertex_set_size
                positions[counter] = (x, y)
                counter += 1

        g = Graph(sum(L))
        s = 0
        for i in L:
            g.add_clique(range(s, s + i))
            s += i

        g = Graph([range(g.order()), g.complement().edges(sort_vertices=False)],
                  format='vertices_and_edges', immutable=immutable, name=name,
                  pos=positions)

    return g


def DiamondGraph(immutable=False):
    """
    Return a diamond graph with 4 nodes.

    A diamond graph is a square with one pair of diagonal nodes connected.

    INPUT:

    - ``immutable`` -- boolean (default: ``False``); whether to return an
      immutable or a mutable graph

    PLOTTING: Upon construction, the position dictionary is filled to override
    the spring-layout algorithm. By convention, the diamond graph is drawn as a
    diamond, with the first node on top, second on the left, third on the right,
    and fourth on the bottom; with the second and third node connected.

    EXAMPLES:

    Construct and show a diamond graph::

        sage: g = graphs.DiamondGraph()
        sage: g.show()                          # long time                             # needs sage.plot

    TESTS:

    Check the behavior of parameter ``immutable``::

        sage: graphs.DiamondGraph(immutable=True).is_immutable()
        True
    """
    pos_dict = {0: (0, 1), 1: (-1, 0), 2: (1, 0), 3: (0, -1)}
    edges = [(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)]
    return Graph([range(4), edges], format='vertices_and_edges',
                 immutable=immutable, pos=pos_dict, name="Diamond Graph")


def GemGraph(immutable=False):
    """
    Return a gem graph with 5 nodes.

    A gem graph is a fan graph (4,1).

    INPUT:

    - ``immutable`` -- boolean (default: ``False``); whether to return an
      immutable or a mutable graph

    PLOTTING: Upon construction, the position dictionary is filled to override
    the spring-layout algorithm. By convention, the gem graph is drawn as a gem,
    with the sharp part on the bottom.

    EXAMPLES:

    Construct and show a gem graph::

        sage: g = graphs.GemGraph()
        sage: g.show()                          # long time                             # needs sage.plot

    TESTS:

    Check the behavior of parameter ``immutable``::

        sage: graphs.GemGraph(immutable=True).is_immutable()
        True
    """
    pos_dict = {0: (0.5, 0), 1: (0, 0.75), 2: (0.25, 1), 3: (0.75, 1), 4: (1, 0.75)}
    edges = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (2, 3), (3, 4)]
    return Graph([range(5), edges], format='vertices_and_edges',
                 immutable=immutable, pos=pos_dict, name="Gem Graph")


def ForkGraph(immutable=False):
    """
    Return a fork graph with 5 nodes.

    A fork graph, sometimes also called chair graph, is 5 vertex tree.

    INPUT:

    - ``immutable`` -- boolean (default: ``False``); whether to return an
      immutable or a mutable graph

    PLOTTING: Upon construction, the position dictionary is filled to override
    the spring-layout algorithm. By convention, the fork graph is drawn as a
    fork, with the sharp part on the bottom.

    EXAMPLES:

    Construct and show a fork graph::

        sage: g = graphs.ForkGraph()
        sage: g.show()                          # long time                             # needs sage.plot

    TESTS:

    Check the behavior of parameter ``immutable``::

        sage: graphs.ForkGraph(immutable=True).is_immutable()
        True
    """
    pos_dict = {0: (0, 0), 1: (1, 0), 2: (0, 1), 3: (1, 1), 4: (0, 2)}
    edges = [(0, 2), (2, 3), (3, 1), (2, 4)]
    return Graph([range(5), edges], format='vertices_and_edges',
                 immutable=immutable, pos=pos_dict, name="Fork Graph")


def DartGraph(immutable=False):
    """
    Return a dart graph with 5 nodes.

    INPUT:

    - ``immutable`` -- boolean (default: ``False``); whether to return an
      immutable or a mutable graph

    PLOTTING: Upon construction, the position dictionary is filled to override
    the spring-layout algorithm. By convention, the dart graph is drawn as a
    dart, with the sharp part on the bottom.

    EXAMPLES:

    Construct and show a dart graph::

        sage: g = graphs.DartGraph()
        sage: g.show()                          # long time                             # needs sage.plot

    TESTS:

    Check the behavior of parameter ``immutable``::

        sage: graphs.DartGraph(immutable=True).is_immutable()
        True
    """
    pos_dict = {0: (0, 1), 1: (-1, 0), 2: (1, 0), 3: (0, -1), 4: (0, 0)}
    edges = [(0, 1), (0, 2), (1, 4), (2, 4), (0, 4), (3, 4)]
    return Graph([range(5), edges], format='vertices_and_edges',
                 immutable=immutable, pos=pos_dict, name="Dart Graph")


def EmptyGraph(immutable=False):
    """
    Return an empty graph (0 nodes and 0 edges).

    This is useful for constructing graphs by adding edges and vertices
    individually or in a loop.

    INPUT:

    - ``immutable`` -- boolean (default: ``False``); whether to return an
      immutable or a mutable graph

    PLOTTING: When plotting, this graph will use the default
    spring-layout algorithm, unless a position dictionary is
    specified.

    EXAMPLES:

    Add one vertex to an empty graph and then show::

        sage: empty1 = graphs.EmptyGraph()
        sage: empty1.add_vertex()
        0
        sage: empty1.show()                     # long time                             # needs sage.plot

    Use for loops to build a graph from an empty graph::

        sage: empty2 = graphs.EmptyGraph()
        sage: for i in range(5):
        ....:     empty2.add_vertex()  # add 5 nodes, labeled 0-4
        0
        1
        2
        3
        4
        sage: for i in range(3):
        ....:     empty2.add_edge(i,i+1)  # add edges {[0:1],[1:2],[2:3]}
        sage: for i in range(1, 4):
        ....:     empty2.add_edge(4,i)  # add edges {[1:4],[2:4],[3:4]}
        sage: empty2.show()                     # long time                             # needs sage.plot

    TESTS:

    Check the behavior of parameter ``immutable``::

        sage: graphs.EmptyGraph(immutable=True).is_immutable()
        True
    """
    return Graph(sparse=True, immutable=immutable)


def ToroidalGrid2dGraph(p, q, immutable=False):
    r"""
    Return a toroidal 2-dimensional grid graph with `p \times q` nodes (`p` rows
    and `q` columns).

    The toroidal 2-dimensional grid with parameters `p,q` is the 2-dimensional
    grid graph with identical parameters to which are added the edges
    `((i, 0), (i, q - 1))` and `((0, i), (p - 1, i))`.

    INPUT:

    - ``p, q`` -- nonnegative integers; the sides of the toroidal grid

    - ``immutable`` -- boolean (default: ``False``); whether to return an
      immutable or a mutable graph

    EXAMPLES:

    The toroidal 2-dimensional grid is a regular graph, while the usual
    2-dimensional grid is not ::

        sage: tgrid = graphs.ToroidalGrid2dGraph(8,9)
        sage: print(tgrid)
        Toroidal 2D Grid Graph with parameters 8,9
        sage: grid = graphs.Grid2dGraph(8,9)
        sage: grid.is_regular()
        False
        sage: tgrid.is_regular()
        True

    TESTS:

    Check the behavior of parameter ``immutable``::

        sage: graphs.ToroidalGrid2dGraph(2, 3, immutable=True).is_immutable()
        True
    """
    name = f"Toroidal 2D Grid Graph with parameters {p},{q}"
    g = Grid2dGraph(p, q, set_positions=True, immutable=False)
    g.name(name)

    g.add_edges([((i, 0), (i, q - 1)) for i in range(p)])
    g.add_edges([((0, i), (p - 1, i)) for i in range(q)])

    pos = g._pos
    p += 0.
    q += 0.
    uf = (p / 2) * (p / 2)
    vf = (q / 2) * (q / 2)
    for u, v in g:
        x, y = pos[u, v]
        x += 0.25 * (1.0 + u * (u - p + 1) / uf)
        y += 0.25 * (1.0 + v * (v - q + 1) / vf)
        pos[u, v] = (x, y)

    if immutable:
        return Graph(g, immutable=True, pos=g.get_pos(), name=name)
    return g


def Toroidal6RegularGrid2dGraph(p, q, immutable=False):
    r"""
    Return a toroidal 6-regular grid.

    The toroidal 6-regular grid is a 6-regular graph on `p \times q` vertices
    and its elements have coordinates `(i, j)` for `i \in \{0...p-1\}` and
    `j \in \{0...q-1\}`.

    Its edges are those of the :meth:`ToroidalGrid2dGraph`, to which are added
    the edges between `(i, j)` and `((i + 1) \% p, (j + 1) \% q)`.

    INPUT:

    - ``p``, ``q`` -- integers

    - ``immutable`` -- boolean (default: ``False``); whether to return an
      immutable or a mutable graph

    EXAMPLES:

    The toroidal 6-regular grid on `25` elements::

        sage: g = graphs.Toroidal6RegularGrid2dGraph(5,5)
        sage: g.is_regular(k=6)
        True
        sage: g.is_vertex_transitive()                                                  # needs sage.groups
        True
        sage: g.line_graph().is_vertex_transitive()                                     # needs sage.groups
        True
        sage: g.automorphism_group().cardinality()                                      # needs sage.groups
        300
        sage: g.is_hamiltonian()                                                        # needs sage.numerical.mip
        True

    TESTS:

    Senseless input::

        sage: graphs.Toroidal6RegularGrid2dGraph(5,2)
        Traceback (most recent call last):
        ...
        ValueError: parameters p and q must be integers larger than 3
        sage: graphs.Toroidal6RegularGrid2dGraph(2,0)
        Traceback (most recent call last):
        ...
        ValueError: parameters p and q must be integers larger than 3

    Check the behavior of parameter ``immutable``::

        sage: graphs.Toroidal6RegularGrid2dGraph(4, 4, immutable=True).is_immutable()
        True
    """
    if p <= 3 or q <= 3:
        raise ValueError("parameters p and q must be integers larger than 3")

    g = ToroidalGrid2dGraph(p, q, immutable=False)
    for u, v in g:
        g.add_edge((u, v), ((u + 1) % p, (v + 1) % q))

    name = f"Toroidal Hexagonal Grid graph on {p}x{q} elements"
    if immutable:
        return Graph(g, immutable=True, pos=g.get_pos(), name=name)
    g.name(name)
    return g


def Grid2dGraph(p, q, set_positions=True, immutable=False, name=None):
    r"""
    Return a `2`-dimensional grid graph with `p \times q` nodes (`p` rows and
    `q` columns).

    A 2d grid graph resembles a `2` dimensional grid. All inner nodes are
    connected to their `4` neighbors. Outer (non-corner) nodes are connected to
    their `3` neighbors. Corner nodes are connected to their 2 neighbors.

    INPUT:

    - ``p``, ``q`` -- two positive integers

    - ``set_positions`` -- boolean (default: ``True``); whether to set the
      position of the nodes

    - ``immutable`` -- boolean (default: ``False``); whether to return an
      immutable or a mutable graph

    - ``name`` -- string (default: ``None``); used as the name of the returned
      graph when set

    PLOTTING: Upon construction, the position dictionary is filled to override
    the spring-layout algorithm. By convention, nodes are labelled in (row,
    column) pairs with `(0, 0)` in the top left corner. Edges will always be
    horizontal and vertical - another advantage of filling the position
    dictionary.

    EXAMPLES:

    Construct and show a grid 2d graph Rows = `5`, Columns = `7`::

        sage: g = graphs.Grid2dGraph(5,7)
        sage: g.show()                          # long time                             # needs sage.plot

    TESTS:

    Senseless input::

        sage: graphs.Grid2dGraph(5,0)
        Traceback (most recent call last):
        ...
        ValueError: parameters p and q must be positive integers
        sage: graphs.Grid2dGraph(-1,0)
        Traceback (most recent call last):
        ...
        ValueError: parameters p and q must be positive integers

    The graph name contains the dimension::

        sage: g = graphs.Grid2dGraph(5,7)
        sage: g.name()
        '2D Grid Graph for [5, 7]'

    Check the behavior of parameter ``ìmmutable``::

        sage: graphs.Grid2dGraph(2, 3, immutable=True).is_immutable()
        True

    Check the behavior of parameter ``name``::

        sage: graphs.Grid2dGraph(2, 3, name='foo')
        foo: Graph on 6 vertices
    """
    if p <= 0 or q <= 0:
        raise ValueError("parameters p and q must be positive integers")

    vertices = ((i, j) for i in range(p) for j in range(q))
    from itertools import chain
    edges = chain((((i, j), (i + 1, j)) for i in range(p - 1) for j in range(q)),
                  (((i, j), (i, j + 1)) for i in range(p) for j in range(q - 1)))
    pos_dict = None
    if set_positions:
        pos_dict = {(i, j): (j, -i) for i in range(p) for j in range(q)}
    return Graph([vertices, edges], format='vertices_and_edges',
                 immutable=immutable, pos=pos_dict,
                 name=f"2D Grid Graph for [{p}, {q}]" if name is None else name)


def GridGraph(dim_list, immutable=False):
    r"""
    Return an `n`-dimensional grid graph.

    INPUT:

    - ``dim_list`` -- list of integers representing the number of nodes to
      extend in each dimension

    - ``immutable`` -- boolean (default: ``False``); whether to return an
      immutable or a mutable graph

    PLOTTING: When plotting, this graph will use the default spring-layout
    algorithm, unless a position dictionary is specified.

    EXAMPLES::

        sage: G = graphs.GridGraph([2,3,4])
        sage: G.show()                          # long time                             # needs sage.plot

    ::

        sage: C = graphs.CubeGraph(4)
        sage: G = graphs.GridGraph([2,2,2,2])
        sage: C.show()                          # long time                             # needs sage.plot
        sage: G.show()                          # long time                             # needs sage.plot

    TESTS:

    The graph name contains the dimension::

        sage: g = graphs.GridGraph([5, 7])
        sage: g.name()
        'Grid Graph for [5, 7]'
        sage: g = graphs.GridGraph([2, 3, 4])
        sage: g.name()
        'Grid Graph for [2, 3, 4]'
        sage: g = graphs.GridGraph([2, 4, 3])
        sage: g.name()
        'Grid Graph for [2, 4, 3]'

    One dimensional grids (i.e., path) have simple vertex labels::

        sage: g = graphs.GridGraph([5])
        sage: g.vertices(sort=True)
        [0, 1, 2, 3, 4]

    The graph is correct::

        sage: dim = [randint(1,4) for i in range(4)]
        sage: g = graphs.GridGraph(dim)
        sage: import networkx                                                           # needs networkx
        sage: h = Graph(networkx.grid_graph(list(dim)))                                 # needs networkx
        sage: g.is_isomorphic(h)                                                        # needs networkx
        True

    Trivial cases::

        sage: g = graphs.GridGraph([]); g; g.vertices(sort=False)
        Grid Graph for []: Graph on 0 vertices
        []
        sage: g = graphs.GridGraph([1]); g; g.vertices(sort=False)
        Grid Graph for [1]: Graph on 1 vertex
        [0]
        sage: g = graphs.GridGraph([2]); g; g.vertices(sort=True)
        Grid Graph for [2]: Graph on 2 vertices
        [0, 1]
        sage: g = graphs.GridGraph([1,1]); g; g.vertices(sort=False)
        Grid Graph for [1, 1]: Graph on 1 vertex
        [(0, 0)]
        sage: g = graphs.GridGraph([1, 1, 1]); g; g.vertices(sort=False)
        Grid Graph for [1, 1, 1]: Graph on 1 vertex
        [(0, 0, 0)]
        sage: g = graphs.GridGraph([1,1,2]); g; g.vertices(sort=True)
        Grid Graph for [1, 1, 2]: Graph on 2 vertices
        [(0, 0, 0), (0, 0, 1)]

    All dimensions must be positive integers::

        sage: g = graphs.GridGraph([2,-1,3])
        Traceback (most recent call last):
        ...
        ValueError: all dimensions must be positive integers

    Check the behavior of parameter ``ìmmutable``::

        sage: graphs.GridGraph([], immutable=True).is_immutable()
        True
        sage: graphs.GridGraph([2], immutable=True).is_immutable()
        True
        sage: graphs.GridGraph([2, 2], immutable=True).is_immutable()
        True
        sage: graphs.GridGraph([2, 2, 2], immutable=True).is_immutable()
        True
    """
    dim = [int(a) for a in dim_list]
    if any(a <= 0 for a in dim):
        raise ValueError("all dimensions must be positive integers")

    name = "Grid Graph for {}".format(dim)
    n_dim = len(dim)
    if not n_dim:
        return Graph(name=name, immutable=immutable)
    if n_dim == 1:
        # Vertices are labeled from 0 to dim[0]-1
        return PathGraph(dim[0], immutable=immutable, name=name)
    if n_dim == 2:
        # We use the Grid2dGraph generator to also get the positions
        return Grid2dGraph(*dim, immutable=immutable, name=name)

    # Now, n_dim > 2 and we don't set positions
    # Vertices are tuples of dimension n_dim, and the graph contains at
    # least vertex (0, 0, ..., 0)
    V = [tuple([0] * n_dim)]
    def edges():
        from itertools import product
        for u in product(*[range(d) for d in dim]):
            for i in range(n_dim):
                if u[i] + 1 < dim[i]:
                    v = list(u)
                    v[i] = u[i] + 1
                    yield (u, tuple(v))

    return Graph([V, edges()], format='vertices_and_edges',
                 immutable=immutable, name=name)


def HouseGraph(immutable=False):
    """
    Return a house graph with 5 nodes.

    A house graph is named for its shape. It is a triangle (roof) over a
    square (walls).

    PLOTTING: Upon construction, the position dictionary is filled to override
    the spring-layout algorithm. By convention, the house graph is drawn with
    the first node in the lower-left corner of the house, the second in the
    lower-right corner of the house. The third node is in the upper-left corner
    connecting the roof to the wall, and the fourth is in the upper-right corner
    connecting the roof to the wall. The fifth node is the top of the roof,
    connected only to the third and fourth.

    INPUT:

    - ``immutable`` -- boolean (default: ``False``); whether to return an
      immutable or a mutable graph

    EXAMPLES:

    Construct and show a house graph::

        sage: g = graphs.HouseGraph()
        sage: g.show()                          # long time                             # needs sage.plot

    TESTS:

    Check the behavior of parameter ``immutable``::

        sage: graphs.HouseGraph(immutable=True).is_immutable()
        True
    """
    pos_dict = {0: (-1, 0), 1: (1, 0), 2: (-1, 1), 3: (1, 1), 4: (0, 2)}
    edges = [(0, 1), (0, 2), (1, 3), (2, 3), (2, 4), (3, 4)]
    return Graph([range(5), edges], format='vertices_and_edges',
                 immutable=immutable, pos=pos_dict, name="House Graph")


def HouseXGraph(immutable=False):
    """
    Return a house X graph with 5 nodes.

    A house X graph is a house graph with two additional edges. The upper-right
    corner is connected to the lower-left. And the upper-left corner is
    connected to the lower-right.

    PLOTTING: Upon construction, the position dictionary is filled to override
    the spring-layout algorithm. By convention, the house X graph is drawn with
    the first node in the lower-left corner of the house, the second in the
    lower-right corner of the house. The third node is in the upper-left corner
    connecting the roof to the wall, and the fourth is in the upper-right corner
    connecting the roof to the wall. The fifth node is the top of the roof,
    connected only to the third and fourth.

    INPUT:

    - ``immutable`` -- boolean (default: ``False``); whether to return an
      immutable or a mutable graph

    EXAMPLES:

    Construct and show a house X graph::

        sage: g = graphs.HouseXGraph()
        sage: g.show()                          # long time                             # needs sage.plot

    TESTS:

    Check the behavior of parameter ``immutable``::

        sage: graphs.HouseXGraph(immutable=True).is_immutable()
        True
    """
    pos_dict = {0: (-1, 0), 1: (1, 0), 2: (-1, 1), 3: (1, 1), 4: (0, 2)}
    edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3), (2, 4), (3, 4)]
    return Graph([range(5), edges], format='vertices_and_edges',
                 immutable=immutable, pos=pos_dict, name="House Graph")


def LadderGraph(n, immutable=False):
    r"""
    Return a ladder graph with `2 * n` nodes.

    A ladder graph is a basic structure that is typically displayed as a ladder,
    i.e.: two parallel path graphs connected at each corresponding node pair.

    PLOTTING: Upon construction, the position dictionary is filled to override
    the spring-layout algorithm. By convention, each ladder graph will be
    displayed horizontally, with the first n nodes displayed left to right on
    the top horizontal line.

    INPUT:

    - ``n`` -- a nonnegative integer; number of nodes is `2n`

    - ``immutable`` -- boolean (default: ``False``); whether to return an
      immutable or a mutable graph

    EXAMPLES:

    Construct and show a ladder graph with 14 nodes::

        sage: g = graphs.LadderGraph(7)
        sage: g.show()                          # long time                             # needs sage.plot

    Create several ladder graphs in a Sage graphics array::

        sage: # needs sage.plot
        sage: g = []
        sage: j = []
        sage: for i in range(9):
        ....:     k = graphs.LadderGraph(i+2)
        ....:     g.append(k)
        sage: for i in range(3):
        ....:     n = []
        ....:     for m in range(3):
        ....:         n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ....:     j.append(n)
        sage: G = graphics_array(j)
        sage: G.show()                          # long time

    TESTS:

    Check the behavior of parameter ``immutable``::

        sage: graphs.LadderGraph(4, immutable=True).is_immutable()
        True
    """
    pos_dict = {i: (i, 1) for i in range(n)}
    for i in range(n, 2 * n):
        x = i - n
        pos_dict[i] = (x, 0)
    from itertools import chain
    edges_1 = zip(range(n), range(1, n))
    edges_2 = zip(range(n, 2 * n), range(n + 1, 2 * n))
    edges_3 = ((i, i + n) for i in range(n))
    return Graph([range(2 * n), chain(edges_1, edges_2, edges_3)],
                 format='vertices_and_edges', immutable=immutable,
                 pos=pos_dict, name="Ladder graph")


def MoebiusLadderGraph(n, immutable=False):
    r"""
    Return a Möbius ladder graph with `2n` nodes

    A Möbius ladder graph of order `2n` is a ladder graph of the same order
    that is connected at the ends with a single twist, i.e., a ladder graph
    bent around so that top meets bottom with a single twist. Alternatively,
    it can be described as a single cycle graph (of order `2n`) with the
    addition of edges (called `rungs`) joining the antipodal pairs of nodes.
    Also, note that the Möbius ladder graph ``graphs.MoebiusLadderGraph(n)`` is
    precisely the same graph as the circulant graph
    ``graphs.CirculantGraph(2 * n, [1, n])``.

    PLOTTING:

    Upon construction, the position dictionary is filled to override the
    spring-layout algorithm. By convention, each Möbius ladder graph will be
    displayed with the first (0) node at the top, with the rest following in a
    counterclockwise manner.

    INPUT:

    - ``n`` -- a nonnegative integer; number of nodes is `2n`

    - ``immutable`` -- boolean (default: ``False``); whether to return an
      immutable or a mutable graph

    OUTPUT:

    - ``G`` -- a Möbius ladder graph of order `2n`; note that a
      :exc:`ValueError` is returned if `n < 0`

    EXAMPLES:

    Construct and show a Möbius ladder graph with 26 nodes::

        sage: g = graphs.MoebiusLadderGraph(13)
        sage: g.show()                          # long time                             # needs sage.plot

    Create several Möbius ladder graphs in a Sage graphics array::

        sage: # needs sage.plots
        sage: g = []
        sage: j = []
        sage: for i in range(9):
        ....:    k = graphs.MoebiusLadderGraph(i+3)
        ....:    g.append(k)
        sage: for i in range(3):
        ....:    n = []
        ....:    for m in range(3):
        ....:        n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ....:    j.append(n)
        sage: G = graphics_array(j)
        sage: G.show()                          # long time

    TESTS:

    The input parameter must be a nonnegative integer::

        sage: G = graphs.MoebiusLadderGraph(-1)
        Traceback (most recent call last):
        ...
        ValueError: parameter n must be a nonnegative integer

    Check the behavior of parameter ``immutable``::

        sage: graphs.MoebiusLadderGraph(4, immutable=True).is_immutable()
        True

    REFERENCES:

    - :wikipedia:`Möbius_ladder`

    .. SEEALSO::
        :meth:`~sage.graphs.graph_generators.GraphGenerators.LadderGraph`,
        :meth:`~sage.graphs.graph_generators.GraphGenerators.CircularLadderGraph`,
        :meth:`~sage.graphs.graph_generators.GraphGenerators.CirculantGraph`

    AUTHORS:

    - Janmenjaya Panda (2024-05-26)
    """
    if n < 0:
        raise ValueError("parameter n must be a nonnegative integer")

    from itertools import chain
    edges_1 = zip(range(2 * n), chain(range(1, 2 * n), (0,)))
    edges_2 = ((i, i + n) for i in range(n))
    G = Graph([range(2 * n), chain(edges_1, edges_2)],
              format='vertices_and_edges', immutable=immutable,
              name="Moebius ladder graph")
    G._circle_embedding(list(range(2 * n)), angle=pi/2)
    return G


def PathGraph(n, pos=None, immutable=False, name=None):
    r"""
    Return a path graph with `n` nodes.

    A path graph is a graph where all inner nodes are connected to their two
    neighbors and the two end-nodes are connected to their one inner neighbors
    (i.e.: a cycle graph without the first and last node connected).

    INPUT:

    - ``n`` -- nonnegative integer; number of nodes of the path graph

    - ``pos`` -- string (default: ``None``); indicates the embedding to use
      between 'circle', 'line' or the default algorithm. See the plotting
      section below for more detail.

    - ``immutable`` -- boolean (default: ``False``); whether to return an
      immutable or a mutable graph

    - ``name`` -- string (default: ``None``); used as the name of the returned
      graph when set

    PLOTTING: Upon construction, the position dictionary is filled to override
    the spring-layout algorithm. By convention, the graph may be drawn in one of
    two ways: The 'line' argument will draw the graph in a horizontal line (left
    to right) if there are less than 11 nodes. Otherwise the 'line' argument
    will append horizontal lines of length 10 nodes below, alternating left to
    right and right to left. The 'circle' argument will cause the graph to be
    drawn in a cycle-shape, with the first node at the top and then about the
    circle in a clockwise manner. By default (without an appropriate string
    argument) the graph will be drawn as a 'circle' if `10 < n < 41` and as a
    'line' for all other `n`.

    EXAMPLES:

    Show default drawing by size: 'line': `n \leq 10`::

        sage: p = graphs.PathGraph(10)
        sage: p.show()                          # long time                             # needs sage.plot

    'circle': `10 < n < 41`::

        sage: q = graphs.PathGraph(25)
        sage: q.show()                          # long time                             # needs sage.plot

    'line': `n \geq 41`::

        sage: r = graphs.PathGraph(55)
        sage: r.show()                          # long time                             # needs sage.plot

    Override the default drawing::

        sage: s = graphs.PathGraph(5,'circle')
        sage: s.show()                          # long time                             # needs sage.plot

    TESTS:

    Check the behavior of parameter ``immutable``::

        sage: graphs.PathGraph(4, immutable=True).is_immutable()
        True

    Check the behavior of parameter ``name``::

        sage: graphs.PathGraph(4, name='foo')
        foo: Graph on 4 vertices
    """
    edges = ((i, i + 1) for i in range(n - 1))
    G = Graph([range(n), edges], format='vertices_and_edges',
              immutable=immutable, name="Path graph" if name is None else name)

    pos_dict = {}

    # Choose appropriate drawing pattern
    circle = False
    if pos == "circle":
        circle = True
    elif pos == "line":
        circle = False
    # Otherwise use default by size of n
    elif 10 < n < 41:
        circle = True

    # Draw 'circle'
    if circle:
        if n == 1:
            G.set_pos({0: (0, 0)})
        else:
            G._circle_embedding(list(range(n)), angle=pi/2)
    # Draw 'line'
    else:
        counter = 0  # node index
        rem = n % 10  # remainder to appear on last row
        rows = n // 10  # number of rows (not counting last row)
        lr = True  # left to right

        for i in range(rows):  # note that rows doesn't include last row
            y = -i
            for j in range(10):
                if lr:
                    x = j
                else:
                    x = 9 - j
                pos_dict[counter] = (x, y)
                counter += 1
            if lr:
                lr = False
            else:
                lr = True
        y = -rows
        for j in range(rem):  # last row
            if lr:
                x = j
            else:
                x = 9 - j
            pos_dict[counter] = (x, y)
            counter += 1
        G.set_pos(pos_dict)

    return G


def StarGraph(n, immutable=False):
    r"""
    Return a star graph with `n + 1` nodes.

    A Star graph is a basic structure where one node is connected to all other
    nodes.

    INPUT:

    - ``n`` -- a nonnegative integer; number of nodes is `n + 1`

    - ``immutable`` -- boolean (default: ``False``); whether to return an
      immutable or a mutable graph

    PLOTTING: Upon construction, the position dictionary is filled to override
    the spring-layout algorithm. By convention, each star graph will be
    displayed with the first (0) node in the center, the second node (1) at the
    top, with the rest following in a counterclockwise manner. (0) is the node
    connected to all other nodes.

    The star graph is a good opportunity to compare efficiency of filling a
    position dictionary vs. using the spring-layout algorithm for plotting. As
    far as display, the spring-layout should push all other nodes away from the
    (0) node, and thus look very similar to this constructor's positioning.

    EXAMPLES::

        sage: import networkx                                                           # needs networkx

    Compare the plots::

        sage: # needs networkx sage.plot
        sage: n = networkx.star_graph(23)
        sage: spring23 = Graph(n)
        sage: posdict23 = graphs.StarGraph(23)
        sage: spring23.show()                   # long time
        sage: posdict23.show()  # long time

    View many star graphs as a Sage Graphics Array

    With this constructor (i.e., the position dictionary filled)

    ::

        sage: # needs sage.plot
        sage: g = []
        sage: j = []
        sage: for i in range(9):
        ....:     k = graphs.StarGraph(i+3)
        ....:     g.append(k)
        sage: for i in range(3):
        ....:     n = []
        ....:     for m in range(3):
        ....:         n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ....:     j.append(n)
        sage: G = graphics_array(j)
        sage: G.show()                          # long time

    Compared to plotting with the spring-layout algorithm

    ::

        sage: # needs networkx sage.plot
        sage: g = []
        sage: j = []
        sage: for i in range(9):
        ....:     spr = networkx.star_graph(i+3)
        ....:     k = Graph(spr)
        ....:     g.append(k)
        sage: for i in range(3):
        ....:     n = []
        ....:     for m in range(3):
        ....:         n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ....:     j.append(n)
        sage: G = graphics_array(j)
        sage: G.show()                          # long time

    TESTS:

    Check the behavior of parameter ``immutable``::

        sage: graphs.StarGraph(4, immutable=True).is_immutable()
        True
    """
    G = Graph({0: list(range(1, n + 1))}, format='dict_of_lists',
              immutable=immutable, name="Star graph")
    G.set_pos({0: (0, 0)})
    G._circle_embedding(list(range(1, n + 1)), angle=pi/2)
    return G
