# -*- coding: utf-8 -*-
r"""
Various families of graphs

The methods defined here appear in :mod:`sage.graphs.graph_generators`.
"""

# ****************************************************************************
#       Copyright (C) 2006 Robert L. Miller <rlmillster@gmail.com>
#                          Emily A. Kirkman
#                     2009 Michael C. Yurko <myurko@gmail.com>
#                     2016 Rowan Schrecker <rowan.schrecker@hertford.ox.ac.uk>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from copy import copy
from math import sin, cos, pi
from sage.graphs.graph import Graph
from itertools import combinations
import subprocess


def JohnsonGraph(n, k):
    r"""
    Returns the Johnson graph with parameters `n, k`.

    Johnson graphs are a special class of undirected graphs defined from systems
    of sets. The vertices of the Johnson graph `J(n,k)` are the `k`-element
    subsets of an `n`-element set; two vertices are adjacent when they meet in a
    `(k-1)`-element set. See the :wikipedia:`Johnson_graph` for more
    information.

    EXAMPLES:

    The Johnson graph is a Hamiltonian graph::

        sage: g = graphs.JohnsonGraph(7, 3)
        sage: g.is_hamiltonian()                                                        # needs sage.numerical.mip
        True

    Every Johnson graph is vertex transitive::

        sage: g = graphs.JohnsonGraph(6, 4)
        sage: g.is_vertex_transitive()                                                  # needs sage.groups
        True

    The complement of the Johnson graph `J(n,2)` is isomorphic to the Kneser
    Graph `K(n,2)`.  In particular the complement of `J(5,2)` is isomorphic to
    the Petersen graph.::

        sage: g = graphs.JohnsonGraph(5,2)
        sage: g.complement().is_isomorphic(graphs.PetersenGraph())
        True
    """

    g = Graph(name="Johnson graph with parameters "+str(n)+","+str(k))
    from sage.combinat.subset import Set, Subsets

    S = Set(range(n))
    g.add_vertices(Subsets(S, k))

    for sub in Subsets(S, k-1):
        elem_left = S - sub
        for i in elem_left:
            for j in elem_left:
                if j <= i:
                    continue
                g.add_edge(sub + Set([i]), sub + Set([j]))

    return g


def KneserGraph(n, k):
    r"""
    Returns the Kneser Graph with parameters `n, k`.

    The Kneser Graph with parameters `n,k` is the graph
    whose vertices are the `k`-subsets of `[0,1,\dots,n-1]`, and such
    that two vertices are adjacent if their corresponding sets
    are disjoint.

    For example, the Petersen Graph can be defined
    as the Kneser Graph with parameters `5,2`.

    EXAMPLES::

        sage: KG = graphs.KneserGraph(5,2)
        sage: sorted(KG.vertex_iterator(), key=str)
        [{1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 3}, {2, 4}, {2, 5},
         {3, 4}, {3, 5}, {4, 5}]
        sage: P = graphs.PetersenGraph()
        sage: P.is_isomorphic(KG)
        True

    TESTS::

        sage: KG = graphs.KneserGraph(0,0)
        Traceback (most recent call last):
        ...
        ValueError: Parameter n should be a strictly positive integer
        sage: KG = graphs.KneserGraph(5,6)
        Traceback (most recent call last):
        ...
        ValueError: Parameter k should be a strictly positive integer inferior to n
    """

    if n <= 0:
        raise ValueError("Parameter n should be a strictly positive integer")
    if k <= 0 or k > n:
        raise ValueError("Parameter k should be a strictly positive integer inferior to n")

    g = Graph(name="Kneser graph with parameters {},{}".format(n, k))

    from sage.combinat.subset import Subsets
    S = Subsets(n, k)
    if 2 * k > n:
        g.add_vertices(S)

    s0 = S.underlying_set()    # {1,2,...,n}
    for s in S:
        for t in Subsets(s0.difference(s), k):
            g.add_edge(s, t)

    return g


def FurerGadget(k, prefix=None):
    r"""
    Return a Furer gadget of order ``k`` and their coloring.

    Construct the Furer gadget described in [CFI1992]_,
    a graph composed by a middle layer of `2^(k-1)` nodes
    and two sets of nodes `(a_0, ... , a_{k-1})` and
    `(b_0, ... , b_{k-1})`.
    Each node in the middle is connected to either `a_i` or `b_i`,
    for each i in [0,k[.
    To read about the complete construction, see [CFI1992]_.
    The returned coloring colors the middle section with one color, and
    then each pair `(a_i, b_i)` with another color.
    Since this method is mainly used to create Furer gadgets for the
    Cai-Furer-Immerman construction, returning gadgets that don't
    always have the same vertex labels is important, that's why there is
    a parameter to manually set a prefix to be appended to each vertex label.

    INPUT:

    - ``k``      -- The order of the returned Furer gadget, greater than 0.

    - ``prefix`` -- Prefix of to be appended to each vertex label,
                    so as to individualise the returned Furer gadget.
                    Must be comparable for equality and hashable.

    OUTPUT:

    - ``G``        -- The Furer gadget of order ``k``

    - ``coloring`` -- A list of list of vertices, representing the
                      partition induced by the coloring of ``G``'s
                      vertices

    EXAMPLES:

    Furer gadget of order 3, without any prefix. ::

        sage: G, p = graphs.FurerGadget(3)
        sage: sorted(G, key=str)
        [(), (0, 'a'), (0, 'b'), (0, 1), (0, 2),
         (1, 'a'), (1, 'b'), (1, 2), (2, 'a'), (2, 'b')]
        sage: sorted(G.edge_iterator(), key=str)
        [((), (0, 'b'), None), ((), (1, 'b'), None),
         ((), (2, 'b'), None), ((0, 'b'), (1, 2), None),
         ((0, 1), (0, 'a'), None), ((0, 1), (1, 'a'), None),
         ((0, 1), (2, 'b'), None), ((0, 2), (0, 'a'), None),
         ((0, 2), (1, 'b'), None), ((0, 2), (2, 'a'), None),
         ((1, 2), (1, 'a'), None), ((1, 2), (2, 'a'), None)]

    Furer gadget of order 3, with a prefix. ::

        sage: G, p = graphs.FurerGadget(3, 'Prefix')
        sage: sorted(G, key=str)
        [('Prefix', ()), ('Prefix', (0, 'a')), ('Prefix', (0, 'b')),
         ('Prefix', (0, 1)), ('Prefix', (0, 2)), ('Prefix', (1, 'a')),
         ('Prefix', (1, 'b')), ('Prefix', (1, 2)), ('Prefix', (2, 'a')),
         ('Prefix', (2, 'b'))]
        sage: sorted(G.edge_iterator(), key=str)
        [(('Prefix', ()), ('Prefix', (0, 'b')), None),
         (('Prefix', ()), ('Prefix', (1, 'b')), None),
         (('Prefix', ()), ('Prefix', (2, 'b')), None),
         (('Prefix', (0, 'b')), ('Prefix', (1, 2)), None),
         (('Prefix', (0, 1)), ('Prefix', (0, 'a')), None),
         (('Prefix', (0, 1)), ('Prefix', (1, 'a')), None),
         (('Prefix', (0, 1)), ('Prefix', (2, 'b')), None),
         (('Prefix', (0, 2)), ('Prefix', (0, 'a')), None),
         (('Prefix', (0, 2)), ('Prefix', (1, 'b')), None),
         (('Prefix', (0, 2)), ('Prefix', (2, 'a')), None),
         (('Prefix', (1, 2)), ('Prefix', (1, 'a')), None),
         (('Prefix', (1, 2)), ('Prefix', (2, 'a')), None)]
    """
    from itertools import repeat as rep, chain
    if k <= 0:
        raise ValueError("The order of the Furer gadget must be greater than zero")
    G = Graph()
    V_a = list(enumerate(rep('a', k)))
    V_b = list(enumerate(rep('b', k)))
    if prefix is not None:
        V_a = list(zip(rep(prefix, k), V_a))
        V_b = list(zip(rep(prefix, k), V_b))
    G.add_vertices(V_a)
    G.add_vertices(V_b)
    powerset = list(chain.from_iterable(combinations(range(k), r) for r in range(0, k + 1, 2)))
    if prefix is not None:
        G.add_edges(chain.from_iterable([((prefix, s), (prefix, (i, 'a'))) for i in s] for s in powerset))
        G.add_edges(chain.from_iterable([((prefix, s), (prefix, (i, 'b'))) for i in range(k) if i not in s] for s in powerset))
    else:
        G.add_edges(chain.from_iterable([(s, (i, 'a')) for i in s] for s in powerset))
        G.add_edges(chain.from_iterable([(s, (i, 'b')) for i in range(k) if i not in s] for s in powerset))
    partition = []
    for i in range(k):
        partition.append([V_a[i], V_b[i]])
    if prefix is not None:
        powerset = [(prefix, s) for s in powerset]
    partition.append(powerset)
    return G, partition


def CaiFurerImmermanGraph(G, twisted=False):
    r"""
    Return the a Cai-Furer-Immerman graph from `G`, possibly a twisted
    one, and a partition of its nodes.

    A Cai-Furer-Immerman graph from/on `G` is a graph created by
    applying the transformation described in [CFI1992]_ on a graph
    `G`, that is substituting every vertex v in `G` with a
    Furer gadget `F(v)` of order d equal to the degree of the vertex,
    and then substituting every edge `(v,u)` in `G`
    with a pair of edges, one connecting the two "a" nodes of
    `F(v)` and `F(u)` and the other their two "b" nodes.
    The returned coloring of the vertices is made by the union of the
    colorings of each single Furer gadget, individualised for each
    vertex of `G`.
    To understand better what these "a" and "b" nodes are, see the
    documentation on  Furer gadgets.

    Furthermore, this method can apply what is described in the paper
    mentioned above as a "twist" on an edge, that is taking only one of
    the pairs of edges introduced in the new graph and swap two of their
    extremes, making each edge go from an "a" node to a "b" node.
    This is only doable if the original graph G is connected.

    A CaiFurerImmerman graph on a graph with no balanced vertex
    separators smaller than s and its twisted version
    cannot be distinguished by k-WL for any k < s.

    INPUT:

    - ``G``       -- An undirected graph on which to construct the
                     Cai-Furer-Immerman graph

    - ``twisted`` -- A boolean indicating if the version to construct
                     is a twisted one or not

    OUTPUT:

    - ``H``        -- The Cai-Furer-Immerman graph on ``G``

    - ``coloring`` -- A list of list of vertices, representing the
                      partition induced by the coloring on ``H``

    EXAMPLES:

    CaiFurerImmerman graph with no balanced vertex separator smaller
    than 2 ::

        sage: G = graphs.CycleGraph(4)
        sage: CFI, p = graphs.CaiFurerImmermanGraph(G)
        sage: sorted(CFI, key=str)
        [(0, ()), (0, (0, 'a')), (0, (0, 'b')), (0, (0, 1)), (0, (1, 'a')),
         (0, (1, 'b')), (1, ()), (1, (0, 'a')), (1, (0, 'b')), (1, (0, 1)),
         (1, (1, 'a')), (1, (1, 'b')), (2, ()), (2, (0, 'a')), (2, (0, 'b')),
         (2, (0, 1)), (2, (1, 'a')), (2, (1, 'b')), (3, ()), (3, (0, 'a')),
         (3, (0, 'b')), (3, (0, 1)), (3, (1, 'a')), (3, (1, 'b'))]
        sage: sorted(CFI.edge_iterator(), key=str)
        [((0, ()), (0, (0, 'b')), None),
         ((0, ()), (0, (1, 'b')), None),
         ((0, (0, 'a')), (1, (0, 'a')), None),
         ((0, (0, 'b')), (1, (0, 'b')), None),
         ((0, (0, 1)), (0, (0, 'a')), None),
         ((0, (0, 1)), (0, (1, 'a')), None),
         ((0, (1, 'a')), (3, (0, 'a')), None),
         ((0, (1, 'b')), (3, (0, 'b')), None),
         ((1, ()), (1, (0, 'b')), None),
         ((1, ()), (1, (1, 'b')), None),
         ((1, (0, 1)), (1, (0, 'a')), None),
         ((1, (0, 1)), (1, (1, 'a')), None),
         ((1, (1, 'a')), (2, (0, 'a')), None),
         ((1, (1, 'b')), (2, (0, 'b')), None),
         ((2, ()), (2, (0, 'b')), None),
         ((2, ()), (2, (1, 'b')), None),
         ((2, (0, 1)), (2, (0, 'a')), None),
         ((2, (0, 1)), (2, (1, 'a')), None),
         ((2, (1, 'a')), (3, (1, 'a')), None),
         ((2, (1, 'b')), (3, (1, 'b')), None),
         ((3, ()), (3, (0, 'b')), None),
         ((3, ()), (3, (1, 'b')), None),
         ((3, (0, 1)), (3, (0, 'a')), None),
         ((3, (0, 1)), (3, (1, 'a')), None)]
    """
    isConnected = G.is_connected()
    newG = Graph()
    total_partition = []
    edge_index = {}
    for v in G:
        Fk, p = FurerGadget(G.degree(v), v)
        total_partition += p
        newG = newG.union(Fk)
        edge_index[v] = 0
    for v, u in G.edge_iterator(labels=False):
        i = edge_index[v]
        edge_index[v] += 1
        j = edge_index[u]
        edge_index[u] += 1
        edge_va = (v, (i, 'a'))
        edge_vb = (v, (i, 'b'))
        edge_ua = (u, (j, 'a'))
        edge_ub = (u, (j, 'b'))
        if isConnected and twisted:
            temp = edge_ua
            edge_ua = edge_ub
            edge_ub = temp
            isConnected = False
        newG.add_edge(edge_va, edge_ua)
        newG.add_edge(edge_vb, edge_ub)
    if twisted and G.is_connected():
        s = " twisted"
    else:
        s = ""
    newG.name("CaiFurerImmerman" + s + " graph constructed from a " + G.name())
    return newG, total_partition


def EgawaGraph(p, s):
    r"""
    Return the Egawa graph with parameters `p`, `s`.

    Egawa graphs are a peculiar family of graphs devised by Yoshimi
    Egawa in  [Ega1981]_ .
    The Shrikhande graph is a special case of this family of graphs,
    with parameters `(1,0)`.
    All the graphs in this family are not recognizable by 1-WL
    (Weisfeiler Lehamn algorithm of the first order) and 2-WL, that is
    their orbits are not correctly returned by k-WL for k lower than 3.

    Furthermore, all the graphs in this family are distance-regular, but
    they are not distance-transitive if `p \neq 0`.

    The Egawa graph with parameters `(0, s)` is isomorphic to the
    Hamming graph with parameters `(s, 4)`, when the underlying
    set of the Hamming graph is `[0,1,2,3]`

    INPUT:

    - ``p`` -- power to which the graph named `Y` in the reference
               provided above will be raised

    - ``s`` -- power to which the graph named `X` in the reference
               provided above will be raised

    OUTPUT:

    - ``G`` -- The Egawa graph with parameters (p,s)

    EXAMPLES:

    Every Egawa graph is distance regular.  ::

        sage: g = graphs.EgawaGraph(1, 2)
        sage: g.is_distance_regular()
        True

    An Egawa graph with parameters (0,s) is isomorphic to the Hamming
    graph with parameters (s, 4).  ::

        sage: g = graphs.EgawaGraph(0, 4)
        sage: g.is_isomorphic(graphs.HammingGraph(4,4))
        True
    """
    from sage.graphs.generators.basic import CompleteGraph
    from itertools import product, chain, repeat
    g = Graph(name="Egawa Graph with parameters " + str(p) + "," + str(s), multiedges=False)
    X = CompleteGraph(4)
    Y = Graph('O?Wse@UgqqT_LUebWkbT_')
    g.add_vertices(product(*chain(repeat(Y, p), repeat(X, s))))
    for v in g:
        for i in range(p):
            prefix = v[:i]
            suffix = v[i+1:]
            for el in Y.neighbor_iterator(v[i]):
                u = prefix + (el,) + suffix
                g.add_edge(v, u)
        for i in range(p, s + p):
            prefix = v[:i]
            suffix = v[i+1:]
            for el in X:
                if el == v[i]:
                    continue
                u = prefix + (el,) + suffix
                g.add_edge(v, u)
    return g


def HammingGraph(n, q, X=None):
    r"""
    Returns the Hamming graph with parameters ``n``, ``q`` over ``X``.

    Hamming graphs are graphs over the cartesian product of n copies
    of ``X``, where `q = |X|`, where the vertices, labelled with the
    corresponding tuple in `X^n`, are connected if the Hamming distance
    between their labels is 1. All Hamming graphs are regular,
    vertex-transitive and distance-regular.

    Hamming graphs with parameters `(1,q)` represent the complete graph
    with q vertices over the set ``X``.

    INPUT:

    - ``n`` -- power to which ``X`` will be raised to provide vertices
               for the Hamming graph

    - ``q`` -- cardinality of ``X``

    - ``X`` -- list of labels representing the vertices of the
                underlying graph the Hamming graph will be based on; if
                ``None`` (or left unused), the list `[0, ... , q-1]`
                will be used

    OUTPUT:

    - ``G`` -- The Hamming graph with parameters `(n,q,X)`

    EXAMPLES:

    Every Hamming graph is distance-regular, regular and
    vertex-transitive.  ::

        sage: g = graphs.HammingGraph(3, 7)
        sage: g.is_distance_regular()
        True
        sage: g.is_regular()
        True
        sage: g.is_vertex_transitive()                                                  # needs sage.groups
        True

    A Hamming graph with parameters (1,q) is isomorphic to the
    Complete graph with parameter q.  ::

        sage: g = graphs.HammingGraph(1, 23)
        sage: g.is_isomorphic(graphs.CompleteGraph(23))
        True

    If a parameter ``q`` is provided which is not equal to ``X``'s
    cardinality, an exception is raised. ::

        sage: X = ['a','b','c','d','e']
        sage: g = graphs.HammingGraph(2, 3, X)
        Traceback (most recent call last):
        ...
        ValueError: q must be the cardinality of X

    REFERENCES:

    For a more accurate description, see the following wikipedia page:
    :wikipedia:`Hamming_graph`
    """
    from itertools import product, repeat
    if not X:
        X = list(range(q))
    if q != len(X):
        raise ValueError("q must be the cardinality of X")
    g = Graph(name="Hamming Graph with parameters " + str(n) + "," + str(q), multiedges=False)
    g.add_vertices(product(*repeat(X, n)))
    for v in g:
        for i in range(n):
            prefix = v[:i]
            suffix = v[i+1:]
            for el in X:
                if el == v[i]:
                    continue
                u = prefix + (el,) + suffix
                g.add_edge(v, u)
    return g


def BalancedTree(r, h):
    r"""
    Returns the perfectly balanced tree of height `h \geq 1`,
    whose root has degree `r \geq 2`.

    The number of vertices of this graph is
    `1 + r + r^2 + \cdots + r^h`, that is,
    `\frac{r^{h+1} - 1}{r - 1}`. The number of edges is one
    less than the number of vertices.

    INPUT:

    - ``r`` -- positive integer `\geq 2`. The degree of the root node.

    - ``h`` -- positive integer `\geq 1`. The height of the balanced tree.

    OUTPUT:

    The perfectly balanced tree of height `h \geq 1` and whose root has
    degree `r \geq 2`. A ``NetworkXError`` is returned if `r < 2` or
    `h < 1`.

    ALGORITHM:

    Uses `NetworkX <http://networkx.lanl.gov>`_.

    EXAMPLES:

    A balanced tree whose root node has degree `r = 2`, and of height
    `h = 1`, has order 3 and size 2::

        sage: G = graphs.BalancedTree(2, 1); G                                          # needs networkx
        Balanced tree: Graph on 3 vertices
        sage: G.order(); G.size()                                                       # needs networkx
        3
        2
        sage: r = 2; h = 1
        sage: v = 1 + r
        sage: v; v - 1
        3
        2

    Plot a balanced tree of height 5, whose root node has degree `r = 3`::

        sage: G = graphs.BalancedTree(3, 5)                                             # needs networkx
        sage: G.show()                          # long time                             # needs networkx sage.plot

    A tree is bipartite. If its vertex set is finite, then it is planar. ::

        sage: # needs networkx
        sage: r = randint(2, 5); h = randint(1, 7)
        sage: T = graphs.BalancedTree(r, h)
        sage: T.is_bipartite()
        True
        sage: T.is_planar()
        True
        sage: v = (r^(h + 1) - 1) / (r - 1)
        sage: T.order() == v
        True
        sage: T.size() == v - 1
        True

    TESTS:

    Normally we would only consider balanced trees whose root node
    has degree `r \geq 2`, but the construction degenerates
    gracefully::

        sage: graphs.BalancedTree(1, 10)                                                # needs networkx
        Balanced tree: Graph on 11 vertices

    Similarly, we usually want the tree must have height `h \geq 1`
    but the algorithm also degenerates gracefully here::

        sage: graphs.BalancedTree(3, 0)                                                 # needs networkx
        Balanced tree: Graph on 1 vertex
    """
    import networkx
    return Graph(networkx.balanced_tree(r, h), name="Balanced tree")


def BarbellGraph(n1, n2):
    r"""
    Returns a barbell graph with ``2*n1 + n2`` nodes. The argument ``n1``
    must be greater than or equal to 2.

    A barbell graph is a basic structure that consists of a path graph
    of order ``n2`` connecting two complete graphs of order ``n1`` each.

    INPUT:

    - ``n1`` -- integer `\geq 2`. The order of each of the two
      complete graphs.

    - ``n2`` -- nonnegative integer. The order of the path graph
      connecting the two complete graphs.

    OUTPUT:

    A barbell graph of order ``2*n1 + n2``. A ``ValueError`` is
    returned if ``n1 < 2`` or ``n2 < 0``.

    PLOTTING:

    Upon construction, the position dictionary is filled to
    override the spring-layout algorithm. By convention, each barbell
    graph will be displayed with the two complete graphs in the
    lower-left and upper-right corners, with the path graph connecting
    diagonally between the two. Thus the ``n1``-th node will be drawn at a
    45 degree angle from the horizontal right center of the first
    complete graph, and the ``n1 + n2 + 1``-th node will be drawn 45
    degrees below the left horizontal center of the second complete graph.

    EXAMPLES:

    Construct and show a barbell graph ``Bar = 4``, ``Bells = 9``::

        sage: g = graphs.BarbellGraph(9, 4); g
        Barbell graph: Graph on 22 vertices
        sage: g.show()  # long time

    An ``n1 >= 2``, ``n2 >= 0`` barbell graph has order ``2*n1 + n2``. It
    has the complete graph on ``n1`` vertices as a subgraph. It also has
    the path graph on ``n2`` vertices as a subgraph. ::

        sage: n1 = randint(2, 2*10^2)
        sage: n2 = randint(0, 2*10^2)
        sage: g = graphs.BarbellGraph(n1, n2)
        sage: v = 2*n1 + n2
        sage: g.order() == v
        True
        sage: K_n1 = graphs.CompleteGraph(n1)
        sage: P_n2 = graphs.PathGraph(n2)

        sage: # needs sage.modules
        sage: s_K = g.subgraph_search(K_n1, induced=True)
        sage: s_P = g.subgraph_search(P_n2, induced=True)
        sage: K_n1.is_isomorphic(s_K)
        True
        sage: P_n2.is_isomorphic(s_P)
        True

    TESTS::

        sage: n1, n2 = randint(3, 10), randint(0, 10)
        sage: g = graphs.BarbellGraph(n1, n2)
        sage: g.num_verts() == 2 * n1 + n2
        True
        sage: g.num_edges() == 2 * binomial(n1, 2) + n2 + 1                             # needs sage.symbolic
        True
        sage: g.is_connected()
        True
        sage: g.girth() == 3
        True

    The input ``n1`` must be `\geq 2`::

        sage: graphs.BarbellGraph(1, randint(0, 10^6))
        Traceback (most recent call last):
        ...
        ValueError: invalid graph description, n1 should be >= 2
        sage: graphs.BarbellGraph(randint(-10^6, 1), randint(0, 10^6))
        Traceback (most recent call last):
        ...
        ValueError: invalid graph description, n1 should be >= 2

    The input ``n2`` must be `\geq 0`::

        sage: graphs.BarbellGraph(randint(2, 10^6), -1)
        Traceback (most recent call last):
        ...
        ValueError: invalid graph description, n2 should be >= 0
        sage: graphs.BarbellGraph(randint(2, 10^6), randint(-10^6, -1))
        Traceback (most recent call last):
        ...
        ValueError: invalid graph description, n2 should be >= 0
        sage: graphs.BarbellGraph(randint(-10^6, 1), randint(-10^6, -1))
        Traceback (most recent call last):
        ...
        ValueError: invalid graph description, n1 should be >= 2
    """
    # sanity checks
    if n1 < 2:
        raise ValueError("invalid graph description, n1 should be >= 2")
    if n2 < 0:
        raise ValueError("invalid graph description, n2 should be >= 0")

    G = Graph(name="Barbell graph")
    G.add_clique(list(range(n1)))
    G.add_path(list(range(n1 - 1, n1 + n2 + 1)))
    G.add_clique(list(range(n1 + n2, n1 + n2 + n1)))

    G._circle_embedding(list(range(n1)), shift=1, angle=pi/4)
    G._line_embedding(list(range(n1, n1 + n2)), first=(2, 2), last=(n2 + 1, n2 + 1))
    G._circle_embedding(list(range(n1 + n2, n1 + n2 + n1)), center=(n2 + 3, n2 + 3), angle=5*pi/4)
    return G


def LollipopGraph(n1, n2):
    r"""
    Returns a lollipop graph with n1+n2 nodes.

    A lollipop graph is a path graph (order n2) connected to a complete
    graph (order n1). (A barbell graph minus one of the bells).

    PLOTTING: Upon construction, the position dictionary is filled to
    override the spring-layout algorithm. By convention, the complete
    graph will be drawn in the lower-left corner with the (n1)th node
    at a 45 degree angle above the right horizontal center of the
    complete graph, leading directly into the path graph.

    EXAMPLES:

    Construct and show a lollipop graph Candy = 13, Stick = 4::

        sage: g = graphs.LollipopGraph(13,4); g
        Lollipop graph: Graph on 17 vertices
        sage: g.show()  # long time

    TESTS::

        sage: n1, n2 = randint(3, 10), randint(0, 10)
        sage: g = graphs.LollipopGraph(n1, n2)
        sage: g.num_verts() == n1 + n2
        True
        sage: g.num_edges() == binomial(n1, 2) + n2                                     # needs sage.symbolic
        True
        sage: g.is_connected()
        True
        sage: g.girth() == 3
        True
        sage: graphs.LollipopGraph(n1, 0).is_isomorphic(graphs.CompleteGraph(n1))
        True
        sage: graphs.LollipopGraph(0, n2).is_isomorphic(graphs.PathGraph(n2))
        True
        sage: graphs.LollipopGraph(0, 0).is_isomorphic(graphs.EmptyGraph())
        True

        The input ``n1`` must be `\geq 0`::

        sage: graphs.LollipopGraph(-1, randint(0, 10^6))
        Traceback (most recent call last):
        ...
        ValueError: invalid graph description, n1 should be >= 0

    The input ``n2`` must be `\geq 0`::

        sage: graphs.LollipopGraph(randint(2, 10^6), -1)
        Traceback (most recent call last):
        ...
        ValueError: invalid graph description, n2 should be >= 0
    """
    # sanity checks
    if n1 < 0:
        raise ValueError("invalid graph description, n1 should be >= 0")
    if n2 < 0:
        raise ValueError("invalid graph description, n2 should be >= 0")

    G = Graph(n1 + n2, name="Lollipop graph")
    G.add_clique(list(range(n1)))
    G.add_path(list(range(n1, n1 + n2)))
    if n1 * n2 > 0:
        G.add_edge(n1 - 1, n1)
    if n1 == 1:
        G.set_pos({0: (0, 0)})
    else:
        G._circle_embedding(list(range(n1)), shift=1, angle=pi/4)
    G._line_embedding(list(range(n1, n1 + n2)), first=(2, 2), last=(n2 + 1, n2 + 1))
    return G


def TadpoleGraph(n1, n2):
    r"""
    Return a tadpole graph with n1+n2 nodes.

    A tadpole graph is a path graph (order n2) connected to a cycle graph
    (order n1).

    PLOTTING: Upon construction, the position dictionary is filled to override
    the spring-layout algorithm. By convention, the cycle graph will be drawn
    in the lower-left corner with the (n1)th node at a 45 degree angle above
    the right horizontal center of the cycle graph, leading directly into the
    path graph.

    EXAMPLES:

    Construct and show a tadpole graph Cycle = 13, Stick = 4::

        sage: g = graphs.TadpoleGraph(13, 4); g
        Tadpole graph: Graph on 17 vertices
        sage: g.show()  # long time

    TESTS::

        sage: n1, n2 = randint(3, 10), randint(0, 10)
        sage: g = graphs.TadpoleGraph(n1, n2)
        sage: g.num_verts() == n1 + n2
        True
        sage: g.num_edges() == n1 + n2
        True
        sage: g.girth() == n1
        True
        sage: graphs.TadpoleGraph(n1, 0).is_isomorphic(graphs.CycleGraph(n1))
        True

    The input ``n1`` must be `\geq 3`::

        sage: graphs.TadpoleGraph(2, randint(0, 10^6))
        Traceback (most recent call last):
        ...
        ValueError: invalid graph description, n1 should be >= 3

    The input ``n2`` must be `\geq 0`::

        sage: graphs.TadpoleGraph(randint(2, 10^6), -1)
        Traceback (most recent call last):
        ...
        ValueError: invalid graph description, n2 should be >= 0
    """
    # sanity checks
    if n1 < 3:
        raise ValueError("invalid graph description, n1 should be >= 3")
    if n2 < 0:
        raise ValueError("invalid graph description, n2 should be >= 0")

    G = Graph(n1 + n2, name="Tadpole graph")
    G.add_cycle(list(range(n1)))
    G.add_path(list(range(n1, n1 + n2)))
    if n1 * n2 > 0:
        G.add_edge(n1 - 1, n1)
    G._circle_embedding(list(range(n1)), shift=1, angle=pi/4)
    G._line_embedding(list(range(n1, n1 + n2)), first=(2, 2), last=(n2 + 1, n2 + 1))
    return G


def AztecDiamondGraph(n):
    """
    Return the Aztec Diamond graph of order ``n``.

    See the :wikipedia:`Aztec_diamond` for more information.

    EXAMPLES::

        sage: graphs.AztecDiamondGraph(2)
        Aztec Diamond graph of order 2

        sage: [graphs.AztecDiamondGraph(i).num_verts() for i in range(8)]
        [0, 4, 12, 24, 40, 60, 84, 112]

        sage: [graphs.AztecDiamondGraph(i).num_edges() for i in range(8)]
        [0, 4, 16, 36, 64, 100, 144, 196]

        sage: G = graphs.AztecDiamondGraph(3)
        sage: sum(1 for p in G.perfect_matchings())
        64
    """
    from sage.graphs.generators.basic import Grid2dGraph
    if n:
        N = 2 * n
        G = Grid2dGraph(N, N)
        H = G.subgraph([(i, j) for i in range(N) for j in range(N)
                        if i - n <= j <= n + i and
                        n - 1 - i <= j <= 3 * n - i - 1])
    else:
        H = Graph()
    H.rename('Aztec Diamond graph of order {}'.format(n))
    return H


def DipoleGraph(n):
    r"""
    Returns a dipole graph with n edges.

    A dipole graph is a multigraph consisting of 2 vertices connected with n
    parallel edges.

    EXAMPLES:

    Construct and show a dipole graph with 13 edges::

        sage: g = graphs.DipoleGraph(13); g
        Dipole graph: Multi-graph on 2 vertices
        sage: g.show()  # long time

    TESTS::

        sage: n = randint(0, 10)
        sage: g = graphs.DipoleGraph(n)
        sage: g.num_verts() == 2
        True
        sage: g.num_edges() == n
        True
        sage: g.is_connected() == (n > 0)
        True
        sage: g.diameter() == (1 if n > 0 else infinity)
        True

    The input ``n`` must be `\geq 0`::

        sage: graphs.DipoleGraph(-randint(1, 10))
        Traceback (most recent call last):
        ...
        ValueError: invalid graph description, n should be >= 0
    """
    # sanity checks
    if n < 0:
        raise ValueError("invalid graph description, n should be >= 0")

    return Graph([[0, 1], [(0, 1)]*n], name="Dipole graph", multiedges=True)


def BubbleSortGraph(n):
    r"""
    Returns the bubble sort graph `B(n)`.

    The vertices of the bubble sort graph are the set of permutations
    on `n` symbols. Two vertices are adjacent if one can be obtained
    from the other by swapping the labels in the `i`-th and `(i+1)`-th
    positions for `1 \leq i \leq n-1`. In total, `B(n)` has order
    `n!`. Swapping two labels as described previously corresponds to
    multiplying on the right the permutation corresponding to the node
    by an elementary transposition in the
    :class:`~sage.groups.perm_gps.permgroup_named.SymmetricGroup`.

    The bubble sort graph is the underlying graph of the
    :meth:`~sage.geometry.polyhedron.library.Polytopes.permutahedron`.

    INPUT:

    - ``n`` -- positive integer. The number of symbols to permute.

    OUTPUT:

    The bubble sort graph `B(n)` on `n` symbols. If `n < 1`, a
    ``ValueError`` is returned.

    EXAMPLES::

        sage: g = graphs.BubbleSortGraph(4); g
        Bubble sort: Graph on 24 vertices
        sage: g.plot()  # long time
        Graphics object consisting of 61 graphics primitives

    The bubble sort graph on `n = 1` symbol is the trivial graph `K_1`::

        sage: graphs.BubbleSortGraph(1)
        Bubble sort: Graph on 1 vertex

    If `n \geq 1`, then the order of `B(n)` is `n!`::

        sage: n = randint(1, 8)
        sage: g = graphs.BubbleSortGraph(n)
        sage: g.order() == factorial(n)
        True

    .. SEEALSO::

        * :meth:`~sage.geometry.polyhedron.library.Polytopes.permutahedron`

    TESTS:

    Input ``n`` must be positive::

        sage: graphs.BubbleSortGraph(0)
        Traceback (most recent call last):
        ...
        ValueError: Invalid number of symbols to permute, n should be >= 1
        sage: graphs.BubbleSortGraph(randint(-10^6, 0))
        Traceback (most recent call last):
        ...
        ValueError: Invalid number of symbols to permute, n should be >= 1

    AUTHORS:

    - Michael Yurko (2009-09-01)
    """
    # sanity checks
    if n < 1:
        raise ValueError(
            "Invalid number of symbols to permute, n should be >= 1")
    if n == 1:
        from sage.graphs.generators.basic import CompleteGraph
        return Graph(CompleteGraph(n), name="Bubble sort")
    from sage.combinat.permutation import Permutations
    # create set from which to permute
    label_set = [str(i) for i in range(1, n + 1)]
    d = {}
    # iterate through all vertices
    for v in Permutations(label_set):
        v = list(v)  # So we can easily mutate it
        tmp_dict = {}
        # add all adjacencies
        for i in range(n - 1):
            # swap entries
            v[i], v[i + 1] = v[i + 1], v[i]
            # add new vertex
            new_vert = ''.join(v)
            tmp_dict[new_vert] = None
            # swap back
            v[i], v[i + 1] = v[i + 1], v[i]
        # add adjacency dict
        d[''.join(v)] = tmp_dict
    return Graph(d, name="Bubble sort")


def chang_graphs():
    r"""
    Return the three Chang graphs.

    Three of the four strongly regular graphs of parameters `(28,12,6,4)` are
    called the Chang graphs. The fourth is the line graph of `K_8`. For more
    information about the Chang graphs, see the :wikipedia:`Chang_graphs` or
    https://www.win.tue.nl/~aeb/graphs/Chang.html.

    EXAMPLES: check that we get 4 non-isomorphic s.r.g.'s with the
    same parameters::

        sage: chang_graphs = graphs.chang_graphs()
        sage: K8 = graphs.CompleteGraph(8)
        sage: T8 = K8.line_graph()
        sage: four_srg = chang_graphs + [T8]
        sage: for g in four_srg:
        ....:     print(g.is_strongly_regular(parameters=True))
        (28, 12, 6, 4)
        (28, 12, 6, 4)
        (28, 12, 6, 4)
        (28, 12, 6, 4)
        sage: from itertools import combinations
        sage: for g1,g2 in combinations(four_srg,2):
        ....:     assert not g1.is_isomorphic(g2)

    Construct the Chang graphs by Seidel switching::

        sage: c3c5 = graphs.CycleGraph(3).disjoint_union(graphs.CycleGraph(5))
        sage: c8 = graphs.CycleGraph(8)
        sage: s = [K8.subgraph_search(c8).edges(sort=False),                            # needs sage.modules
        ....:      [(0,1,None),(2,3,None),(4,5,None),(6,7,None)],
        ....:      K8.subgraph_search(c3c5).edges(sort=False)]
        sage: [T8.seidel_switching(x, inplace=False).is_isomorphic(G)                   # needs sage.modules
        ....:  for x, G in zip(s, chang_graphs)]
        [True, True, True]

    """
    g1 = Graph("[}~~EebhkrRb_~SoLOIiAZ?LBBxDb?bQcggjHKEwoZFAaiZ?Yf[?dxb@@tdWGkwn",
               loops=False, multiedges=False)
    g2 = Graph("[~z^UipkkZPr_~Y_LOIiATOLBBxPR@`acoojBBSoWXTaabN?Yts?Yji_QyioClXZ",
               loops=False, multiedges=False)
    g3 = Graph(r"[~~vVMWdKFpV`^UGIaIERQ`\DBxpA@g`CbGRI`AxICNaFM[?fM\?Ytj@CxrGGlYt",
               loops=False, multiedges=False)
    return [g1, g2, g3]


def CirculantGraph(n, adjacency):
    r"""
    Returns a circulant graph with n nodes.

    A circulant graph has the property that the vertex `i` is connected
    with the vertices `i+j` and `i-j` for each j in ``adjacency``.

    INPUT:


    -  ``n`` - number of vertices in the graph

    -  ``adjacency`` - the list of j values


    PLOTTING: Upon construction, the position dictionary is filled to
    override the spring-layout algorithm. By convention, each circulant
    graph will be displayed with the first (0) node at the top, with
    the rest following in a counterclockwise manner.

    Filling the position dictionary in advance adds O(n) to the
    constructor.

    .. SEEALSO::

        * :meth:`sage.graphs.generic_graph.GenericGraph.is_circulant`
          -- checks whether a (di)graph is circulant, and/or returns
          all possible sets of parameters.

    EXAMPLES: Compare plotting using the predefined layout and
    networkx::

        sage: # needs networkx
        sage: import networkx
        sage: n = networkx.cycle_graph(23)
        sage: spring23 = Graph(n)
        sage: posdict23 = graphs.CirculantGraph(23,2)
        sage: spring23.show()                   # long time
        sage: posdict23.show()  # long time

    We next view many cycle graphs as a Sage graphics array. First we
    use the ``CirculantGraph`` constructor, which fills in
    the position dictionary::

        sage: # needs sage.plot
        sage: g = []
        sage: j = []
        sage: for i in range(9):
        ....:     k = graphs.CirculantGraph(i+4, i+1)
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
        ....:  n = []
        ....:  for m in range(3):
        ....:      n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ....:  j.append(n)
        sage: G = graphics_array(j)
        sage: G.show()                          # long time

    Passing a 1 into adjacency should give the cycle.

    ::

        sage: graphs.CirculantGraph(6,1) == graphs.CycleGraph(6)
        True
        sage: graphs.CirculantGraph(7,[1,3]).edges(sort=True, labels=false)
        [(0, 1),
        (0, 3),
        (0, 4),
        (0, 6),
        (1, 2),
        (1, 4),
        (1, 5),
        (2, 3),
        (2, 5),
        (2, 6),
        (3, 4),
        (3, 6),
        (4, 5),
        (5, 6)]
    """
    if not isinstance(adjacency, list):
        adjacency = [adjacency]

    G = Graph(n, name="Circulant graph (" + str(adjacency) + ")")
    G._circle_embedding(list(range(n)))

    for v in G:
        G.add_edges([(v, (v + j) % n) for j in adjacency])

    return G


def CubeGraph(n, embedding=1):
    r"""
    Return the `n`-cube graph, also called the hypercube in `n` dimensions.

    The hypercube in `n` dimension is build upon the binary strings on `n` bits,
    two of them being adjacent if they differ in exactly one bit. Hence, the
    distance between two vertices in the hypercube is the Hamming distance.

    INPUT:

    - ``n`` -- integer; the dimension of the cube graph

    - ``embedding`` -- integer (default: ``1``); two embeddings of the `n`-cube
      are available:

      - ``1``: the `n`-cube is projected inside a regular `2n`-gonal polygon by
        a skew orthogonal projection. See the :wikipedia:`Hypercube` for more
        details.

      - ``2``: orthogonal projection of the `n`-cube. This orientation shows
        columns of independent vertices such that the neighbors of a vertex are
        located in the columns on the left and on the right. The number of
        vertices in each column represents rows in Pascal's triangle. See for
        instance the :wikipedia:`10-cube` for more details.

      - ``None`` or ``O``: no embedding is provided

    EXAMPLES:

    The distance between `0100110` and `1011010` is `5`, as expected::

        sage: g = graphs.CubeGraph(7)
        sage: g.distance('0100110','1011010')
        5

    Plot several `n`-cubes in a Sage Graphics Array::

        sage: # needs sage.plot
        sage: g = []
        sage: j = []
        sage: for i in range(6):
        ....:  k = graphs.CubeGraph(i+1)
        ....:  g.append(k)
        ...
        sage: for i in range(2):
        ....:  n = []
        ....:  for m in range(3):
        ....:      n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ....:  j.append(n)
        ...
        sage: G = graphics_array(j)
        sage: G.show(figsize=[6,4])             # long time

    Use the plot options to display larger `n`-cubes::

        sage: g = graphs.CubeGraph(9, embedding=1)
        sage: g.show(figsize=[12,12],vertex_labels=False, vertex_size=20)       # long time, needs sage.plot
        sage: g = graphs.CubeGraph(9, embedding=2)
        sage: g.show(figsize=[12,12],vertex_labels=False, vertex_size=20)       # long time, needs sage.plot

    AUTHORS:

    - Robert Miller
    - David Coudert
    """
    if embedding == 1:
        # construct recursively the adjacency dict and the embedding
        theta = float(pi/n)
        d = {'': []}
        dn = {}
        p = {'': (float(0), float(0))}
        pn = {}

        for i in range(n):
            ci = float(cos(i*theta))
            si = float(sin(i*theta))
            for v, e in d.items():
                v0 = v + '0'
                v1 = v + '1'
                l0 = [v1]
                l1 = [v0]
                for m in e:
                    l0.append(m + '0')
                    l1.append(m + '1')
                dn[v0] = l0
                dn[v1] = l1
                x, y = p[v]
                pn[v0] = (x, y)
                pn[v1] = (x + ci, y + si)
            d, dn = dn, {}
            p, pn = pn, {}

        # construct the graph
        G = Graph(d, format='dict_of_lists', pos=p, name="%d-Cube" % n)

    else:
        # construct recursively the adjacency dict
        d = {'': []}
        dn = {}

        for i in range(n):
            for v, e in d.items():
                v0 = v + '0'
                v1 = v + '1'
                l0 = [v1]
                l1 = [v0]
                for m in e:
                    l0.append(m + '0')
                    l1.append(m + '1')
                dn[v0] = l0
                dn[v1] = l1
            d, dn = dn, {}

        # construct the graph
        G = Graph(d, name="%d-Cube" % n, format='dict_of_lists')

        if embedding == 2:
            # Orthogonal projection
            s = '0'*n
            L = [[] for _ in range(n + 1)]
            for u, d in G.breadth_first_search(s, report_distance=True):
                L[d].append(u)

            p = G._circle_embedding(list(range(2*n)), radius=(n + 1)//2, angle=pi, return_dict=True)
            for i in range(n + 1):
                y = p[i][1] / 1.5
                G._line_embedding(L[i], first=(i, y), last=(i, -y), return_dict=False)

    return G


def GoethalsSeidelGraph(k, r):
    r"""
    Returns the graph `\text{Goethals-Seidel}(k,r)`.

    The graph `\text{Goethals-Seidel}(k,r)` comes from a construction presented
    in Theorem 2.4 of [GS1970]_. It relies on a :func:`(v,k)-BIBD
    <sage.combinat.designs.bibd.balanced_incomplete_block_design>` with `r`
    blocks and a
    :func:`~sage.combinat.matrices.hadamard_matrix.hadamard_matrix` of order
    `r+1`. The result is a
    :func:`sage.graphs.strongly_regular_db.strongly_regular_graph` on `v(r+1)`
    vertices with degree `k=(n+r-1)/2`.

    It appears under this name in Andries Brouwer's `database of strongly
    regular graphs <https://www.win.tue.nl/~aeb/graphs/srg/srgtab.html>`__.

    INPUT:

    - ``k,r`` -- integers

    .. SEEALSO::

        - :func:`~sage.graphs.strongly_regular_db.is_goethals_seidel`

    EXAMPLES::

        sage: graphs.GoethalsSeidelGraph(3,3)                                           # needs sage.modules
        Graph on 28 vertices
        sage: graphs.GoethalsSeidelGraph(3,3).is_strongly_regular(parameters=True)      # needs sage.modules
        (28, 15, 6, 10)
    """
    from sage.combinat.designs.bibd import balanced_incomplete_block_design
    from sage.combinat.matrices.hadamard_matrix import hadamard_matrix
    from sage.matrix.constructor import Matrix
    from sage.matrix.constructor import block_matrix

    v = (k-1)*r + 1
    n = v*(r + 1)

    # N is the (v times b) incidence matrix of a bibd
    N = balanced_incomplete_block_design(v, k).incidence_matrix()

    # L is a (r+1 times r) matrix, where r is the row sum of N
    L = hadamard_matrix(r + 1).submatrix(0, 1)
    L = [Matrix(C).transpose() for C in L.columns()]
    zero = Matrix(r + 1, 1, [0]*(r + 1))

    # For every row of N, we replace the 0s with a column of zeros, and we
    # replace the ith 1 with the ith column of L. The result is P.
    P = []
    for row in N:
        Ltmp = L[:]
        P.append([Ltmp.pop(0) if i else zero
                  for i in row])

    P = block_matrix(P)

    # The final graph
    PP = P*P.transpose()
    for i in range(n):
        PP[i, i] = 0

    G = Graph(PP, format="seidel_adjacency_matrix")
    return G


def DorogovtsevGoltsevMendesGraph(n):
    """
    Construct the n-th generation of the Dorogovtsev-Goltsev-Mendes
    graph.

    EXAMPLES::

        sage: G = graphs.DorogovtsevGoltsevMendesGraph(8)                               # needs networkx
        sage: G.size()                                                                  # needs networkx
        6561

    REFERENCE:

    - [1] Dorogovtsev, S. N., Goltsev, A. V., and Mendes, J.
      F. F., Pseudofractal scale-free web, Phys. Rev. E 066122
      (2002).
    """
    import networkx
    return Graph(networkx.dorogovtsev_goltsev_mendes_graph(n),
                 name="Dorogovtsev-Goltsev-Mendes Graph, %d-th generation" % n)


def FoldedCubeGraph(n):
    r"""
    Returns the folded cube graph of order `2^{n-1}`.

    The folded cube graph on `2^{n-1}` vertices can be obtained from a cube
    graph on `2^n` vertices by merging together opposed
    vertices. Alternatively, it can be obtained from a cube graph on
    `2^{n-1}` vertices by adding an edge between opposed vertices. This
    second construction is the one produced by this method.

    See the :wikipedia:`Folded_cube_graph` for more information.

    EXAMPLES:

    The folded cube graph of order five is the Clebsch graph::

        sage: fc = graphs.FoldedCubeGraph(5)
        sage: clebsch = graphs.ClebschGraph()
        sage: fc.is_isomorphic(clebsch)
        True
    """
    if n < 1:
        raise ValueError("The value of n must be at least 2")

    g = CubeGraph(n - 1)
    g.name("Folded Cube Graph")

    # Complementing the binary word
    def complement(x):
        x = x.replace('0', 'a')
        x = x.replace('1', '0')
        x = x.replace('a', '1')
        return x

    for x in g:
        if x[0] == '0':
            g.add_edge(x, complement(x))

    return g


def FriendshipGraph(n):
    r"""
    Return the friendship graph `F_n`.

    The friendship graph is also known as the Dutch windmill graph. Let
    `C_3` be the cycle graph on 3 vertices. Then `F_n` is constructed by
    joining `n \geq 1` copies of `C_3` at a common vertex. If `n = 1`,
    then `F_1` is isomorphic to `C_3` (the triangle graph). If `n = 2`,
    then `F_2` is the butterfly graph, otherwise known as the bowtie
    graph. For more information, see the :wikipedia:`Friendship_graph`.

    INPUT:

    - ``n`` -- positive integer; the number of copies of `C_3` to use in
      constructing `F_n`.

    OUTPUT:

    - The friendship graph `F_n` obtained from `n` copies of the cycle
      graph `C_3`.

    .. SEEALSO::

        - :meth:`GraphGenerators.ButterflyGraph`

    EXAMPLES:

    The first few friendship graphs. ::

        sage: # needs sage.plot
        sage: A = []; B = []
        sage: for i in range(9):
        ....:     g = graphs.FriendshipGraph(i + 1)
        ....:     A.append(g)
        sage: for i in range(3):
        ....:     n = []
        ....:     for j in range(3):
        ....:         n.append(A[3*i + j].plot(vertex_size=20, vertex_labels=False))
        ....:     B.append(n)
        sage: G = graphics_array(B)
        sage: G.show()                          # long time

    For `n = 1`, the friendship graph `F_1` is isomorphic to the cycle
    graph `C_3`, whose visual representation is a triangle. ::

        sage: G = graphs.FriendshipGraph(1); G
        Friendship graph: Graph on 3 vertices
        sage: G.show()                          # long time                             # needs sage.plot
        sage: G.is_isomorphic(graphs.CycleGraph(3))
        True

    For `n = 2`, the friendship graph `F_2` is isomorphic to the
    butterfly graph, otherwise known as the bowtie graph. ::

        sage: G = graphs.FriendshipGraph(2); G
        Friendship graph: Graph on 5 vertices
        sage: G.is_isomorphic(graphs.ButterflyGraph())
        True

    If `n \geq 2`, then the friendship graph `F_n` has `2n + 1` vertices
    and `3n` edges. It has radius 1, diameter 2, girth 3, and
    chromatic number 3. Furthermore, `F_n` is planar and Eulerian. ::

        sage: n = randint(2, 10^3)
        sage: G = graphs.FriendshipGraph(n)
        sage: G.order() == 2*n + 1
        True
        sage: G.size() == 3*n
        True
        sage: G.radius()
        1
        sage: G.diameter()
        2
        sage: G.girth()
        3
        sage: G.chromatic_number()
        3
        sage: G.is_planar()
        True
        sage: G.is_eulerian()
        True

    TESTS:

    The input ``n`` must be a positive integer. ::

        sage: graphs.FriendshipGraph(randint(-10^5, 0))
        Traceback (most recent call last):
        ...
        ValueError: n must be a positive integer
    """
    # sanity checks
    if n < 1:
        raise ValueError("n must be a positive integer")
    # construct the friendship graph
    if n == 1:
        from sage.graphs.generators.basic import CycleGraph
        G = CycleGraph(3)
        G.name("Friendship graph")
        return G
    # build the edges and position dictionaries
    N = 2 * n + 1           # order of F_n
    center = 2 * n
    G = Graph(N, name="Friendship graph")
    for i in range(0, N - 1, 2):
        G.add_cycle([center, i, i + 1])
    G.set_pos({center: (0, 0)})
    G._circle_embedding(list(range(N - 1)), radius=1)
    return G


def FuzzyBallGraph(partition, q):
    r"""
    Construct a Fuzzy Ball graph with the integer partition
    ``partition`` and ``q`` extra vertices.

    Let `q` be an integer and let `m_1,m_2,...,m_k` be a set of positive
    integers.  Let `n=q+m_1+...+m_k`.  The Fuzzy Ball graph with partition
    `m_1,m_2,...,m_k` and `q` extra vertices is the graph constructed from the
    graph `G=K_n` by attaching, for each `i=1,2,...,k`, a new vertex `a_i` to
    `m_i` distinct vertices of `G`.

    For given positive integers `k` and `m` and nonnegative
    integer `q`, the set of graphs ``FuzzyBallGraph(p, q)`` for
    all partitions `p` of `m` with `k` parts are cospectral with
    respect to the normalized Laplacian.

    EXAMPLES::

        sage: F = graphs.FuzzyBallGraph([3,1],2)
        sage: F.adjacency_matrix(vertices=list(F))                                      # needs sage.modules
        [0 0 1 1 1 0 0 0]
        [0 0 0 0 0 1 0 0]
        [1 0 0 1 1 1 1 1]
        [1 0 1 0 1 1 1 1]
        [1 0 1 1 0 1 1 1]
        [0 1 1 1 1 0 1 1]
        [0 0 1 1 1 1 0 1]
        [0 0 1 1 1 1 1 0]

    Pick positive integers `m` and `k` and a nonnegative integer `q`.
    All the FuzzyBallGraphs constructed from partitions of `m` with
    `k` parts should be cospectral with respect to the normalized
    Laplacian::

        sage: m = 4; q = 2; k = 2
        sage: g_list = [graphs.FuzzyBallGraph(p,q)                                      # needs sage.combinat sage.modules
        ....:           for p in Partitions(m, length=k)]
        sage: set(g.laplacian_matrix(normalized=True,   # long time (7s on sage.math, 2011), needs sage.combinat sage.modules
        ....:                        vertices=list(g)).charpoly()
        ....:     for g in g_list)
        {x^8 - 8*x^7 + 4079/150*x^6 - 68689/1350*x^5 + 610783/10800*x^4
          - 120877/3240*x^3 + 1351/100*x^2 - 931/450*x}
    """
    from sage.graphs.generators.basic import CompleteGraph
    if len(partition) < 1:
        raise ValueError("partition must be a nonempty list of positive integers")
    n = q + sum(partition)
    g = CompleteGraph(n)
    curr_vertex = 0
    for e, p in enumerate(partition):
        g.add_edges([(curr_vertex + i, 'a{0}'.format(e + 1)) for i in range(p)])
        curr_vertex += p
    return g


def FibonacciTree(n):
    r"""
    Return the graph of the Fibonacci Tree `F_{i}` of order `n`.

    The Fibonacci tree `F_{i}` is recursively defined as the tree
    with a root vertex and two attached child trees `F_{i-1}` and
    `F_{i-2}`, where `F_{1}` is just one vertex and `F_{0}` is empty.

    INPUT:

    - ``n`` - the recursion depth of the Fibonacci Tree

    EXAMPLES::

        sage: g = graphs.FibonacciTree(3)                                               # needs sage.libs.pari
        sage: g.is_tree()                                                               # needs sage.libs.pari
        True

    ::

        sage: l1 = [ len(graphs.FibonacciTree(_)) + 1 for _ in range(6) ]               # needs sage.libs.pari
        sage: l2 = list(fibonacci_sequence(2,8))                                        # needs sage.libs.pari
        sage: l1 == l2                                                                  # needs sage.libs.pari
        True

    AUTHORS:

    - Harald Schilly and Yann Laigle-Chapuy (2010-03-25)
    """
    T = Graph(name="Fibonacci-Tree-%d" % n)
    if n == 1:
        T.add_vertex(0)
    if n < 2:
        return T

    from sage.combinat.combinat import fibonacci_sequence
    F = list(fibonacci_sequence(n + 2))
    s = 1.618 ** (n / 1.618 - 1.618)
    pos = {}

    def fib(level, node, y):
        pos[node] = (node, y)
        if level < 2:
            return
        level -= 1
        y -= s
        diff = F[level]
        T.add_edge(node, node - diff)
        if level == 1:  # only one child
            pos[node - diff] = (node, y)
            return
        T.add_edge(node, node + diff)
        fib(level, node - diff, y)
        fib(level - 1, node + diff, y)

    T.add_vertices(range(sum(F[:-1])))
    fib(n, F[n + 1] - 1, 0)
    T.set_pos(pos)

    return T


def GeneralizedPetersenGraph(n, k):
    r"""
    Returns a generalized Petersen graph with `2n` nodes. The variables
    `n`, `k` are integers such that `n>2` and `0<k\leq\lfloor(n-1)`/`2\rfloor`

    For `k=1` the result is a graph isomorphic to the circular ladder graph
    with the same `n`. The regular Petersen Graph has `n=5` and `k=2`.
    Other named graphs that can be described using this notation include
    the Desargues graph and the Möbius-Kantor graph.

    INPUT:

    - ``n`` - the number of nodes is `2*n`.

    - ``k`` - integer `0<k\leq\lfloor(n-1)`/`2\rfloor`. Decides
      how inner vertices are connected.

    PLOTTING: Upon construction, the position dictionary is filled to
    override the spring-layout algorithm. By convention, the generalized
    Petersen graphs are displayed as an inner and outer cycle pair, with
    the first n nodes drawn on the outer circle. The first (0) node is
    drawn at the top of the outer-circle, moving counterclockwise after that.
    The inner circle is drawn with the (n)th node at the top, then
    counterclockwise as well.

    EXAMPLES: For `k=1` the resulting graph will be isomorphic to a circular
    ladder graph. ::

        sage: g = graphs.GeneralizedPetersenGraph(13,1)
        sage: g2 = graphs.CircularLadderGraph(13)
        sage: g.is_isomorphic(g2)
        True

    The Desargues graph::

        sage: g = graphs.GeneralizedPetersenGraph(10,3)
        sage: g.girth()
        6
        sage: g.is_bipartite()
        True

    AUTHORS:

    - Anders Jonsson (2009-10-15)
    """
    if n < 3:
        raise ValueError("n must be larger than 2")
    if k < 1 or k > (n - 1) // 2:
        raise ValueError("k must be in 1<= k <=floor((n-1)/2)")
    G = Graph(2 * n, name="Generalized Petersen graph (n="+str(n)+",k="+str(k)+")")
    for i in range(n):
        G.add_edge(i, (i+1) % n)
        G.add_edge(i, i+n)
        G.add_edge(i+n, n + (i+k) % n)
    G._circle_embedding(list(range(n)), radius=1, angle=pi/2)
    G._circle_embedding(list(range(n, 2*n)), radius=.5, angle=pi/2)
    return G


def IGraph(n, j, k):
    r"""
    Return an I-graph with `2n` nodes.

    The I-Graph family as been proposed in [BCMS1988]_ as a generalization of
    the generalized Petersen graphs.  The variables `n`, `j`, `k` are integers
    such that `n > 2` and `0 < j, k \leq \lfloor (n - 1) / 2 \rfloor`.
    When `j = 1` the resulting graph is isomorphic to the generalized Petersen
    graph with the same `n` and `k`.

    INPUT:

    - ``n`` -- the number of nodes is `2 * n`

    - ``j`` -- integer such that `0 < j \leq \lfloor (n-1) / 2 \rfloor`
      determining how outer vertices are connected

    - ``k`` -- integer such that `0 < k \leq \lfloor (n-1) / 2 \rfloor`
      determining how inner vertices are connected

    PLOTTING: Upon construction, the position dictionary is filled to override
    the spring-layout algorithm. By convention, the I-graphs are displayed as an
    inner and outer cycle pair, with the first n nodes drawn on the outer
    circle.  The first (0) node is drawn at the top of the outer-circle, moving
    counterclockwise after that. The inner circle is drawn with the (n)th node
    at the top, then counterclockwise as well.

    EXAMPLES:

    When `j = 1` the resulting graph will be isomorphic to a generalized
    Petersen graph::

        sage: g = graphs.IGraph(7,1,2)
        sage: g2 = graphs.GeneralizedPetersenGraph(7,2)
        sage: g.is_isomorphic(g2)
        True

    The IGraph with parameters `(n, j, k)` is isomorphic to the IGraph with
    parameters `(n, k, j)`::

        sage: g = graphs.IGraph(7, 2, 3)
        sage: h = graphs.IGraph(7, 3, 2)
        sage: g.is_isomorphic(h)
        True

    TESTS::

        sage: graphs.IGraph(1, 1, 1)
        Traceback (most recent call last):
        ...
        ValueError: n must be larger than 2
        sage: graphs.IGraph(3, 0, 1)
        Traceback (most recent call last):
        ...
        ValueError: j must be in 1 <= j <= floor((n - 1) / 2)
        sage: graphs.IGraph(3, 33, 1)
        Traceback (most recent call last):
        ...
        ValueError: j must be in 1 <= j <= floor((n - 1) / 2)
        sage: graphs.IGraph(3, 1, 0)
        Traceback (most recent call last):
        ...
        ValueError: k must be in 1 <= k <= floor((n - 1) / 2)
        sage: graphs.IGraph(3, 1, 3)
        Traceback (most recent call last):
        ...
        ValueError: k must be in 1 <= k <= floor((n - 1) / 2)
    """
    if n < 3:
        raise ValueError("n must be larger than 2")
    if j < 1 or j > (n - 1) // 2:
        raise ValueError("j must be in 1 <= j <= floor((n - 1) / 2)")
    if k < 1 or k > (n - 1) // 2:
        raise ValueError("k must be in 1 <= k <= floor((n - 1) / 2)")

    G = Graph(2 * n, name="I-graph (n={}, j={}, k={})".format(n, j, k))
    for i in range(n):
        G.add_edge(i, (i + j) % n)
        G.add_edge(i, i + n)
        G.add_edge(i + n, n + (i + k) % n)
    G._circle_embedding(list(range(n)), radius=1, angle=pi/2)
    G._circle_embedding(list(range(n, 2 * n)), radius=.5, angle=pi/2)
    return G


def DoubleGeneralizedPetersenGraph(n, k):
    r"""
    Return a double generalized Petersen graph with `4n` nodes.

    The double generalized Petersen graphs is a family of graphs proposed in
    [ZF2012]_ as a variant of generalized Petersen graphs.  The variables `n`,
    `k` are integers such that `n > 2` and `0 < k \leq \lfloor (n-1) / 2
    \rfloor`.

    INPUT:

    - ``n`` -- the number of nodes is `4 * n`

    - ``k`` -- integer such that `0 < k \leq \lfloor (n-1) / 2 \rfloor`
      determining how vertices on second and third inner rims are connected

    PLOTTING: Upon construction, the position dictionary is filled to override
    the spring-layout algorithm. By convention, the double generalized Petersen
    graphs are displayed as 4 cocentric cycles, with the first n nodes drawn on
    the outer circle.  The first (0) node is drawn at the top of the
    outer-circle, moving counterclockwise after that. The second circle is drawn
    with the (n)th node at the top, then counterclockwise as well. The tird
    cycle is drawn with the (2n)th node at the top, then counterclockwise.  And
    the fourth cycle is drawn with the (3n)th node at the top, then again
    counterclockwise.

    EXAMPLES:

    When `n` is even the resulting graph will be isomorphic to a double
    generalized Petersen graph with `k' = n / 2 - k`::

        sage: g = graphs.DoubleGeneralizedPetersenGraph(10, 2)
        sage: g2 = graphs.DoubleGeneralizedPetersenGraph(10, 3)
        sage: g.is_isomorphic(g2)
        True

    TESTS::

        sage: graphs.DoubleGeneralizedPetersenGraph(1, 1)
        Traceback (most recent call last):
        ...
        ValueError: n must be larger than 2
        sage: graphs.DoubleGeneralizedPetersenGraph(3, 0)
        Traceback (most recent call last):
        ...
        ValueError: k must be in 1 <= k <= floor((n - 1) / 2)
        sage: graphs.DoubleGeneralizedPetersenGraph(3, 3)
        Traceback (most recent call last):
        ...
        ValueError: k must be in 1 <= k <= floor((n - 1) / 2)
    """
    if n < 3:
        raise ValueError("n must be larger than 2")
    if k < 1 or k > (n - 1) // 2:
        raise ValueError("k must be in 1 <= k <= floor((n - 1) / 2)")

    G = Graph(4 * n, name="Double generalized Petersen graph (n={}, k={})".format(n, k))
    for i in range(n):
        G.add_edge(i, (i + 1) % n)
        G.add_edge(i + 3 * n, (i + 1) % n + 3 * n)
        G.add_edge(i, i + n)
        G.add_edge(i + 2 * n, i + 3 * n)
        G.add_edge(i + n, (i + k) % n + 2 * n)
        G.add_edge(i + 2 * n, (i + k) % n + n)
    G._circle_embedding(list(range(n)), radius=3, angle=pi/2)
    G._circle_embedding(list(range(n, 2 * n)), radius=2, angle=pi/2)
    G._circle_embedding(list(range(2 * n, 3 * n)), radius=1.5, angle=pi/2)
    G._circle_embedding(list(range(3 * n, 4 * n)), radius=0.5, angle=pi/2)
    return G


def RoseWindowGraph(n, a, r):
    r"""
    Return a rose window graph with `2n` nodes.

    The rose window graphs is a family of tetravalant graphs introduced in
    [Wilson2008]_. The parameters `n`, `a` and `r` are integers such that
    `n > 2`, `1 \leq a, r < n`, and `r \neq n / 2`.

    INPUT:

    - ``n`` -- the number of nodes is `2 * n`

    - ``a`` -- integer such that `1 \leq a < n` determining a-spoke edges

    - ``r`` -- integer such that `1 \leq r < n` and `r \neq n / 2` determining
      how inner vertices are connected

    PLOTTING: Upon construction, the position dictionary is filled to override
    the spring-layout algorithm. By convention, the rose window graphs are
    displayed as an inner and outer cycle pair, with the first n nodes drawn on
    the outer circle.  The first (0) node is drawn at the top of the
    outer-circle, moving counterclockwise after that. The inner circle is drawn
    with the (n)th node at the top, then counterclockwise as well.  Vertices in
    the outer circle are connected in the circular manner, vertices in the inner
    circle are connected when their label have difference `r` (mod n).  Vertices
    on the outer rim are connected with the vertices on the inner rim when they
    are at the same position and when they are `a` apart.

    EXAMPLES:

    The vertices of a rose window graph have all degree 4::

        sage: G = graphs.RoseWindowGraph(5, 1, 2)
        sage: all(G.degree(u) == 4 for u in G)
        True

    The smallest rose window graph as parameters `(3, 2, 1)`::

        sage: G = graphs.RoseWindowGraph(3, 2, 1)
        sage: all(G.degree(u) == 4 for u in G)
        True

    TESTS::

        sage: graphs.RoseWindowGraph(1, 1, 1)
        Traceback (most recent call last):
        ...
        ValueError: n must be larger than 2
        sage: graphs.RoseWindowGraph(6, 0, 2)
        Traceback (most recent call last):
        ...
        ValueError: a must be an integer such that 1 <= a < n
        sage: graphs.RoseWindowGraph(6, 6, 2)
        Traceback (most recent call last):
        ...
        ValueError: a must be an integer such that 1 <= a < n
        sage: graphs.RoseWindowGraph(6, 3, 0)
        Traceback (most recent call last):
        ...
        ValueError: r must be an integer such that 1 <= r < n
        sage: graphs.RoseWindowGraph(6, 3, 6)
        Traceback (most recent call last):
        ...
        ValueError: r must be an integer such that 1 <= r < n
        sage: graphs.RoseWindowGraph(6, 3, 3)
        Traceback (most recent call last):
        ...
        ValueError: r must be different than n / 2
    """
    if n < 3:
        raise ValueError("n must be larger than 2")
    if a < 1 or a >= n:
        raise ValueError("a must be an integer such that 1 <= a < n")
    if r < 1 or r >= n:
        raise ValueError("r must be an integer such that 1 <= r < n")
    if r == n / 2:
        raise ValueError("r must be different than n / 2")

    G = Graph(2 * n, name="rose window graph (n={}, a={}, r={})".format(n, a, r))
    for i in range(n):
        G.add_edge(i, (i + 1) % n)
        G.add_edge(i, i + n)
        G.add_edge((i + a) % n, i + n)
        G.add_edge(i + n, (i + r) % n + n)
    G._circle_embedding(list(range(n)), radius=1, angle=pi/2)
    G._circle_embedding(list(range(n, 2 * n)), radius=0.5, angle=pi/2)
    return G


def TabacjnGraph(n, a, b, r):
    r"""
    Return a Tabačjn graph with `2n` nodes.

    The Tabačjn graphs is a family of pentavalent bicirculants graphs proposed
    in [AHKOS2014]_ as a generalization of generalized Petersen graphs. The
    parameters `n`, `a`, `b`, `r` are integers such that `n \geq 3`, `1 \leq a,
    b, r \leq n - 1`, with `a \neq b` and `r \neq n / 2`.

    INPUT:

    - ``n`` -- the number of nodes is `2 * n`

    - ``a`` -- integer such that `0 < a < n` and `a \neq b`, that determines
      a-spoke edges

    - ``b`` -- integer such that `0 < b < n` and `b \neq a`, that determines
      b-spoke edges

    - ``r`` -- integer such that `0 < r < n` and `r \neq n/2` determining how
      inner vertices are connected

    PLOTTING: Upon construction, the position dictionary is filled to override
    the spring-layout algorithm. By convention, the rose window graphs are
    displayed as an inner and outer cycle pair, with the first n nodes drawn on
    the outer circle.  The first (0) node is drawn at the top of the
    outer-circle, moving counterclockwise after that. The inner circle is drawn
    with the (n)th node at the top, then counterclockwise as well. Vertices in
    the outer circle are connected in the circular manner, vertices in the inner
    circle are connected when their label have difference `r` (mod n). Vertices
    on the outer rim are connected with the vertices on the inner rim when they
    are at the same position and when they are `a` and `b` apart.

    EXAMPLES::

        sage: G = graphs.TabacjnGraph(3, 1, 2, 1)
        sage: G.degree()
        [5, 5, 5, 5, 5, 5]
        sage: G.is_isomorphic(graphs.CompleteGraph(6))
        True
        sage: G = graphs.TabacjnGraph(6, 1, 5, 2)
        sage: I = graphs.IcosahedralGraph()
        sage: G.is_isomorphic(I)
        True

    TESTS::

        sage: graphs.TabacjnGraph(1, 1, 1, 1)
        Traceback (most recent call last):
        ...
        ValueError: n must be larger than 2
        sage: graphs.TabacjnGraph(3, 0, 1, 1)
        Traceback (most recent call last):
        ...
        ValueError: a must be an integer such that 1 <= a < n
        sage: graphs.TabacjnGraph(3, 3, 1, 1)
        Traceback (most recent call last):
        ...
        ValueError: a must be an integer such that 1 <= a < n
        sage: graphs.TabacjnGraph(3, 1, 0, 1)
        Traceback (most recent call last):
        ...
        ValueError: b must be an integer such that 1 <= b < n
        sage: graphs.TabacjnGraph(3, 1, 3, 1)
        Traceback (most recent call last):
        ...
        ValueError: b must be an integer such that 1 <= b < n
        sage: graphs.TabacjnGraph(3, 1, 1, 1)
        Traceback (most recent call last):
        ...
        ValueError: a must be different than b
        sage: graphs.TabacjnGraph(3, 1, 2, 0)
        Traceback (most recent call last):
        ...
        ValueError: r must be an integer such that 1 <= r < n
        sage: graphs.TabacjnGraph(3, 1, 2, 3)
        Traceback (most recent call last):
        ...
        ValueError: r must be an integer such that 1 <= r < n
        sage: graphs.TabacjnGraph(4, 1, 2, 2)
        Traceback (most recent call last):
        ...
        ValueError: r must be different than n / 2
    """
    if n < 3:
        raise ValueError("n must be larger than 2")
    if a < 1 or a >= n:
        raise ValueError("a must be an integer such that 1 <= a < n")
    if b < 1 or b >= n:
        raise ValueError("b must be an integer such that 1 <= b < n")
    if a == b:
        raise ValueError("a must be different than b")
    if r < 1 or r >= n:
        raise ValueError("r must be an integer such that 1 <= r < n")
    if r == n/2:
        raise ValueError("r must be different than n / 2")

    G = Graph(2 * n, name="Tabačjn graph (n={}, a={}, b={}, r={})".format(n, a, b, r))
    for i in range(n):
        G.add_edge(i, (i + 1) % n)
        G.add_edge(i, i + n)
        G.add_edge(i + n, n + (i + r) % n)
        G.add_edge(i, (i + a) % n + n)
        G.add_edge(i, (i + b) % n + n)
    G._circle_embedding(list(range(n)), radius=1, angle=pi/2)
    G._circle_embedding(list(range(n, 2 * n)), radius=0.5, angle=pi/2)
    return G


def HararyGraph(k, n):
    r"""
    Returns the Harary graph on `n` vertices and connectivity `k`, where
    `2 \leq k < n`.

    A `k`-connected graph `G` on `n` vertices requires the minimum degree
    `\delta(G)\geq k`, so the minimum number of edges `G` should have is
    `\lceil kn/2\rceil`. Harary graphs achieve this lower bound, that is,
    Harary graphs are minimal `k`-connected graphs on `n` vertices.

    The construction provided uses the method CirculantGraph.  For more
    details, see the book D. B. West, Introduction to Graph Theory, 2nd
    Edition, Prentice Hall, 2001, p. 150--151; or the `MathWorld article on
    Harary graphs <http://mathworld.wolfram.com/HararyGraph.html>`_.

    EXAMPLES:

    Harary graphs `H_{k,n}`::

        sage: h = graphs.HararyGraph(5,9); h
        Harary graph 5, 9: Graph on 9 vertices
        sage: h.order()
        9
        sage: h.size()
        23
        sage: h.vertex_connectivity()                                                   # needs sage.numerical.mip
        5

    TESTS:

    Connectivity of some Harary graphs::

        sage: n = 10
        sage: for k in range(2,n):                                                      # needs sage.numerical.mip
        ....:     g = graphs.HararyGraph(k,n)
        ....:     if k != g.vertex_connectivity():
        ....:        print("Connectivity of Harary graphs not satisfied.")
    """
    if k < 2:
        raise ValueError("Connectivity parameter k should be at least 2.")
    if k >= n:
        raise ValueError("Number of vertices n should be greater than k.")

    if k % 2 == 0:
        G = CirculantGraph(n, list(range(1, k//2 + 1)))
    else:
        if n % 2 == 0:
            G = CirculantGraph(n, list(range(1, (k - 1)//2 + 1)))
            for i in range(n):
                G.add_edge(i, (i + n//2) % n)
        else:
            G = HararyGraph(k - 1, n)
            for i in range((n - 1)//2 + 1):
                G.add_edge(i, (i + (n - 1)//2) % n)
    G.name('Harary graph {0}, {1}'.format(k, n))
    return G


def HyperStarGraph(n, k):
    r"""
    Return the hyper-star graph `HS(n, k)`.

    The vertices of the hyper-star graph are the set of binary strings of length
    `n` which contain `k` 1s. Two vertices, `u` and `v`, are adjacent only if
    `u` can be obtained from `v` by swapping the first bit with a different
    symbol in another position. For instance, vertex ``'011100'`` of `HS(6, 3)`
    is adjacent to vertices ``'101100'``, ``'110100'`` and ``'111000'``.
    See [LKOL2002]_ for more details.

    INPUT:

    - ``n`` -- non-negative integer; length of the binary strings

    - ``k`` -- non-negative integer; number of 1s per binary string

    EXAMPLES::

        sage: g = graphs.HyperStarGraph(6,3)
        sage: sorted(g.neighbors('011100'))
        ['101100', '110100', '111000']
        sage: g.plot()  # long time
        Graphics object consisting of 51 graphics primitives

    TESTS::

        sage: graphs.HyperStarGraph(-1, 1)
        Traceback (most recent call last):
        ...
        ValueError: parameters n and k must be non-negative integers satisfying n >= k >= 0
        sage: graphs.HyperStarGraph(1, -1)
        Traceback (most recent call last):
        ...
        ValueError: parameters n and k must be non-negative integers satisfying n >= k >= 0
        sage: graphs.HyperStarGraph(1, 2)
        Traceback (most recent call last):
        ...
        ValueError: parameters n and k must be non-negative integers satisfying n >= k >= 0

    AUTHORS:

    - Michael Yurko (2009-09-01)
    """
    if n < 0 or k < 0 or k > n:
        raise ValueError("parameters n and k must be non-negative integers "
                         "satisfying n >= k >= 0")
    if not n:
        adj = {}
    elif not k:
        adj = {'0'*n: []}
    elif k == n:
        adj = {'1'*n: []}
    else:
        from sage.data_structures.bitset import Bitset
        adj = dict()
        # We consider the strings of n bits with k 1s and starting with a 0
        for c in combinations(range(1, n), k):
            u = str(Bitset(c, capacity=n))
            L = []
            c = list(c)
            # The neighbors of u are all the strings obtained by swapping a 1
            # with the first bit (0)
            for i in range(k):
                one = c[i]
                c[i] = 0
                L.append(str(Bitset(c, capacity=n)))
                c[i] = one
            adj[u] = L

    return Graph(adj, format='dict_of_lists', name="HS(%d,%d)" % (n, k))


def LCFGraph(n, shift_list, repeats):
    r"""
    Return the cubic graph specified in LCF notation.

    LCF (Lederberg-Coxeter-Fruchte) notation is a concise way of
    describing cubic Hamiltonian graphs. The way a graph is constructed
    is as follows. Since there is a Hamiltonian cycle, we first create
    a cycle on n nodes. The variable shift_list = [s_0, s_1, ...,
    s_k-1] describes edges to be created by the following scheme: for
    each i, connect vertex i to vertex (i + s_i). Then, repeats
    specifies the number of times to repeat this process, where on the
    jth repeat we connect vertex (i + j\*len(shift_list)) to vertex (
    i + j\*len(shift_list) + s_i).

    INPUT:


    -  ``n`` - the number of nodes.

    -  ``shift_list`` - a list of integer shifts mod n.

    -  ``repeats`` - the number of times to repeat the
       process.


    EXAMPLES::

        sage: G = graphs.LCFGraph(4, [2,-2], 2)                                         # needs networkx
        sage: G.is_isomorphic(graphs.TetrahedralGraph())                                # needs networkx
        True

    ::

        sage: G = graphs.LCFGraph(20, [10,7,4,-4,-7,10,-4,7,-7,4], 2)                   # needs networkx
        sage: G.is_isomorphic(graphs.DodecahedralGraph())                               # needs networkx
        True

    ::

        sage: G = graphs.LCFGraph(14, [5,-5], 7)                                        # needs networkx
        sage: G.is_isomorphic(graphs.HeawoodGraph())                                    # needs networkx
        True

    The largest cubic nonplanar graph of diameter three::

        sage: # needs networkx
        sage: G = graphs.LCFGraph(20, [-10,-7,-5,4,7,-10,-7,-4,5,7,
        ....:                          -10,-7,6,-5,7,-10,-7,5,-6,7], 1)
        sage: G.degree()
        [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
        sage: G.diameter()
        3
        sage: G.show()                          # long time                             # needs sage.plot

    PLOTTING: LCF Graphs are plotted as an n-cycle with edges in the
    middle, as described above.

    REFERENCES:

    - [1] Frucht, R. "A Canonical Representation of Trivalent
      Hamiltonian Graphs." J. Graph Th. 1, 45-60, 1976.

    - [2] Grunbaum, B.  Convex Polytope es. New York: Wiley,
      pp. 362-364, 1967.

    - [3] Lederberg, J. 'DENDRAL-64: A System for Computer
      Construction, Enumeration and Notation of Organic Molecules
      as Tree Structures and Cyclic Graphs. Part II. Topology of
      Cyclic Graphs.' Interim Report to the National Aeronautics
      and Space Administration. Grant NsG 81-60. December 15,
      1965.  http://profiles.nlm.nih.gov/BB/A/B/I/U/_/bbabiu.pdf.
    """
    import networkx
    G = Graph(networkx.LCF_graph(n, shift_list, repeats), name="LCF Graph")
    G._circle_embedding(list(range(n)), radius=1, angle=pi/2)
    return G


def MycielskiGraph(k=1, relabel=True):
    r"""
    Returns the `k`-th Mycielski Graph.

    The graph `M_k` is triangle-free and has chromatic number
    equal to `k`. These graphs show, constructively, that there
    are triangle-free graphs with arbitrarily high chromatic
    number.

    The Mycielski graphs are built recursively starting with
    `M_0`, an empty graph; `M_1`, a single vertex graph; and `M_2`
    is the graph `K_2`.  `M_{k+1}` is then built from `M_k`
    as follows:

    If the vertices of `M_k` are `v_1,\ldots,v_n`, then the
    vertices of `M_{k+1}` are
    `v_1,\ldots,v_n,w_1,\ldots,w_n,z`. Vertices `v_1,\ldots,v_n`
    induce a copy of `M_k`. Vertices `w_1,\ldots,w_n` are an
    independent set. Vertex `z` is adjacent to all the
    `w_i`-vertices. Finally, vertex `w_i` is adjacent to vertex
    `v_j` iff `v_i` is adjacent to `v_j`.

    INPUT:

    - ``k`` Number of steps in the construction process.

    - ``relabel`` Relabel the vertices so their names are the integers
      ``range(n)`` where ``n`` is the number of vertices in the graph.

    EXAMPLES:

    The Mycielski graph `M_k` is triangle-free and has chromatic
    number equal to `k`. ::

        sage: g = graphs.MycielskiGraph(5)
        sage: g.is_triangle_free()
        True
        sage: g.chromatic_number()
        5

    The graphs `M_4` is (isomorphic to) the Grotzsch graph. ::

        sage: g = graphs.MycielskiGraph(4)
        sage: g.is_isomorphic(graphs.GrotzschGraph())
        True

    REFERENCES:

    -  [1] Weisstein, Eric W. "Mycielski Graph."
       From MathWorld--A Wolfram Web Resource.
       http://mathworld.wolfram.com/MycielskiGraph.html

    """
    g = Graph()
    g.name("Mycielski Graph " + str(k))

    if k < 0:
        raise ValueError("parameter k must be a nonnegative integer")

    if k == 0:
        return g

    if k == 1:
        g.add_vertex(0)
        return g

    if k == 2:
        g.add_edge(0, 1)
        return g

    g0 = MycielskiGraph(k - 1)
    g = MycielskiStep(g0)
    g.name("Mycielski Graph " + str(k))
    if relabel:
        g.relabel()

    return g


def MycielskiStep(g):
    r"""
    Perform one iteration of the Mycielski construction.

    See the documentation for ``MycielskiGraph`` which uses this
    method. We expose it to all users in case they may find it
    useful.

    EXAMPLE. One iteration of the Mycielski step applied to the
    5-cycle yields a graph isomorphic to the Grotzsch graph ::

        sage: g = graphs.CycleGraph(5)
        sage: h = graphs.MycielskiStep(g)
        sage: h.is_isomorphic(graphs.GrotzschGraph())
        True
    """
    # Make a copy of the input graph g
    gg = copy(g)

    # rename a vertex v of gg as (1,v)
    renamer = {v: (1, v) for v in g}
    gg.relabel(renamer)

    # add the w vertices to gg as (2,v)
    wlist = [(2, v) for v in g]
    gg.add_vertices(wlist)

    # add the z vertex as (0,0)
    gg.add_vertex((0, 0))

    # add the edges from z to w_i
    gg.add_edges([((0, 0), (2, v)) for v in g])

    # make the v_i w_j edges
    for v in g:
        gg.add_edges([((1, v), (2, vv)) for vv in g.neighbors(v)])

    return gg


def NKStarGraph(n, k):
    r"""
    Returns the (n,k)-star graph.

    The vertices of the (n,k)-star graph are the set of all arrangements of
    n symbols into labels of length k. There are two adjacency rules for
    the (n,k)-star graph. First, two vertices are adjacent if one can be
    obtained from the other by swapping the first symbol with another
    symbol. Second, two vertices are adjacent if one can be obtained from
    the other by swapping the first symbol with an external symbol (a
    symbol not used in the original label).

    INPUT:

    -  ``n``

    -  ``k``

    EXAMPLES::

        sage: g = graphs.NKStarGraph(4,2)
        sage: g.plot()  # long time
        Graphics object consisting of 31 graphics primitives

    REFERENCES:

    - Wei-Kuo, Chiang, and Chen Rong-Jaye. "The (n, k)-star graph: A
      generalized star graph." Information Processing Letters 56,
      no. 5 (December 8, 1995): 259-264.

    AUTHORS:

    - Michael Yurko (2009-09-01)
    """
    from sage.combinat.permutation import Arrangements
    # set from which to permute
    set = [str(i) for i in range(1, n + 1)]
    # create dict
    d = {}
    for v in Arrangements(set, k):
        v = list(v)  # So we can easily mutate it
        tmp_dict = {}
        # add edges of dimension i
        for i in range(1, k):
            # swap 0th and ith element
            v[0], v[i] = v[i], v[0]
            # convert to str and add to list
            vert = "".join(v)
            tmp_dict[vert] = None
            # swap back
            v[0], v[i] = v[i], v[0]
        # add other edges
        tmp_bit = v[0]
        for i in set:
            # check if external
            if i not in v:
                v[0] = i
                # add edge
                vert = "".join(v)
                tmp_dict[vert] = None
            v[0] = tmp_bit
        d["".join(v)] = tmp_dict
    return Graph(d, name="(%d,%d)-star" % (n, k))


def NStarGraph(n):
    r"""
    Returns the n-star graph.

    The vertices of the n-star graph are the set of permutations on n
    symbols. There is an edge between two vertices if their labels differ
    only in the first and one other position.

    INPUT:

    -  ``n``

    EXAMPLES::

        sage: g = graphs.NStarGraph(4)
        sage: g.plot()  # long time
        Graphics object consisting of 61 graphics primitives

    REFERENCES:

    - S.B. Akers, D. Horel and B. Krishnamurthy, The star graph: An
      attractive alternative to the previous n-cube. In: Proc. Internat.
      Conf. on Parallel Processing (1987), pp. 393--400.

    AUTHORS:

    - Michael Yurko (2009-09-01)
    """
    from sage.combinat.permutation import Permutations
    # set from which to permute
    set = [str(i) for i in range(1, n + 1)]
    # create dictionary of lists
    # vertices are adjacent if the first element is swapped with the ith element
    d = {}
    for v in Permutations(set):
        v = list(v)  # So we can easily mutate it
        tmp_dict = {}
        for i in range(1, n):
            if v[0] != v[i]:
                # swap 0th and ith element
                v[0], v[i] = v[i], v[0]
                # convert to str and add to list
                vert = "".join(v)
                tmp_dict[vert] = None
                # swap back
                v[0], v[i] = v[i], v[0]
        d["".join(v)] = tmp_dict
    return Graph(d, name="%d-star" % n)


def OddGraph(n):
    r"""
    Returns the Odd Graph with parameter `n`.

    The Odd Graph with parameter `n` is defined as the
    Kneser Graph with parameters `2n-1,n-1`.
    Equivalently, the Odd Graph is the graph whose vertices
    are the `n-1`-subsets of `[0,1,\dots,2(n-1)]`, and such
    that two vertices are adjacent if their corresponding sets
    are disjoint.

    For example, the Petersen Graph can be defined
    as the Odd Graph with parameter `3`.

    EXAMPLES::

        sage: OG = graphs.OddGraph(3)
        sage: sorted(OG.vertex_iterator(), key=str)
        [{1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 3}, {2, 4}, {2, 5},
         {3, 4}, {3, 5}, {4, 5}]
        sage: P = graphs.PetersenGraph()
        sage: P.is_isomorphic(OG)
        True

    TESTS::

        sage: KG = graphs.OddGraph(1)
        Traceback (most recent call last):
        ...
        ValueError: Parameter n should be an integer strictly greater than 1
    """
    if n <= 1:
        raise ValueError("Parameter n should be an integer strictly greater than 1")
    g = KneserGraph(2*n - 1, n - 1)
    g.name("Odd Graph with parameter %s" % n)
    return g


def PaleyGraph(q):
    r"""
    Paley graph with `q` vertices

    Parameter `q` must be the power of a prime number and congruent
    to 1 mod 4.

    EXAMPLES::

        sage: G = graphs.PaleyGraph(9); G                                               # needs sage.rings.finite_rings
        Paley graph with parameter 9: Graph on 9 vertices
        sage: G.is_regular()                                                            # needs sage.rings.finite_rings
        True

    A Paley graph is always self-complementary::

        sage: G.is_self_complementary()                                                 # needs sage.rings.finite_rings
        True

    TESTS:

    Wrong parameter::

        sage: graphs.PaleyGraph(6)
        Traceback (most recent call last):
        ...
        ValueError: parameter q must be a prime power
        sage: graphs.PaleyGraph(3)
        Traceback (most recent call last):
        ...
        ValueError: parameter q must be congruent to 1 mod 4
    """
    from sage.rings.finite_rings.integer_mod import mod
    from sage.rings.finite_rings.finite_field_constructor import FiniteField
    from sage.arith.misc import is_prime_power
    if not is_prime_power(q):
        raise ValueError("parameter q must be a prime power")
    if not mod(q, 4) == 1:
        raise ValueError("parameter q must be congruent to 1 mod 4")
    g = Graph([FiniteField(q, 'a'), lambda i, j: (i - j).is_square()],
              loops=False, name="Paley graph with parameter {}".format(q))
    return g


def PasechnikGraph(n):
    r"""
    Pasechnik strongly regular graph on `(4n-1)^2` vertices

    A strongly regular graph with parameters of the orthogonal array graph
    :func:`~sage.graphs.graph_generators.GraphGenerators.OrthogonalArrayBlockGraph`,
    also known as pseudo Latin squares graph `L_{2n-1}(4n-1)`, constructed from
    a skew Hadamard matrix of order `4n` following [Pas1992]_.

    .. SEEALSO::

        - :func:`~sage.graphs.strongly_regular_db.is_orthogonal_array_block_graph`

    EXAMPLES::

        sage: graphs.PasechnikGraph(4).is_strongly_regular(parameters=True)             # needs sage.combinat sage.modules
        (225, 98, 43, 42)
        sage: graphs.PasechnikGraph(5).is_strongly_regular(parameters=True)     # long time, needs sage.combinat sage.modules
        (361, 162, 73, 72)
        sage: graphs.PasechnikGraph(9).is_strongly_regular(parameters=True)  # not tested
        (1225, 578, 273, 272)

    TESTS::

        sage: graphs.PasechnikGraph(0)
        Traceback (most recent call last):
        ...
        ValueError: parameter n must be >= 1
    """
    if n < 1:
        raise ValueError("parameter n must be >= 1")
    from sage.combinat.matrices.hadamard_matrix import skew_hadamard_matrix
    from sage.matrix.constructor import identity_matrix
    H = skew_hadamard_matrix(4 * n)
    M = H[1:].T[1:] - identity_matrix(4 * n - 1)
    G = Graph(M.tensor_product(M.T), format='seidel_adjacency_matrix')
    G.relabel()
    G.name("Pasechnik Graph_{}".format(n))
    return G


def SquaredSkewHadamardMatrixGraph(n):
    r"""
    Pseudo-`OA(2n,4n-1)`-graph from a skew Hadamard matrix of order `4n`

    A strongly regular graph with parameters of the orthogonal array graph
    :func:`~sage.graphs.graph_generators.GraphGenerators.OrthogonalArrayBlockGraph`,
    also known as pseudo Latin squares graph `L_{2n}(4n-1)`, constructed from a
    skew Hadamard matrix of order `4n`, due to Goethals and Seidel, see
    [BL1984]_.

    .. SEEALSO::

        - :func:`~sage.graphs.strongly_regular_db.is_orthogonal_array_block_graph`

    EXAMPLES::

        sage: # needs sage.modules
        sage: G = graphs.SquaredSkewHadamardMatrixGraph(4)
        sage: G.is_strongly_regular(parameters=True)
        (225, 112, 55, 56)
        sage: G = graphs.SquaredSkewHadamardMatrixGraph(5)
        sage: G.is_strongly_regular(parameters=True)    # long time
        (361, 180, 89, 90)
        sage: G = graphs.SquaredSkewHadamardMatrixGraph(9)
        sage: G.is_strongly_regular(parameters=True)    # not tested
        (1225, 612, 305, 306)

    TESTS::

        sage: graphs.SquaredSkewHadamardMatrixGraph(0)
        Traceback (most recent call last):
        ...
        ValueError: parameter n must be >= 1
    """
    if n < 1:
        raise ValueError("parameter n must be >= 1")
    from sage.combinat.matrices.hadamard_matrix import skew_hadamard_matrix
    from sage.matrix.constructor import identity_matrix, matrix
    idm = identity_matrix(4 * n - 1)
    e = matrix([1] * (4 * n - 1))
    H = skew_hadamard_matrix(4 * n)
    M = H[1:].T[1:] - idm
    s = M.tensor_product(M.T) - idm.tensor_product(e.T * e - idm)
    G = Graph(s, format='seidel_adjacency_matrix')
    G.relabel()
    G.name("skewhad^2_{}".format(n))
    return G


def SwitchedSquaredSkewHadamardMatrixGraph(n):
    r"""
    A strongly regular graph in Seidel switching class of
    `SquaredSkewHadamardMatrixGraph`

    A strongly regular graph in the :meth:`Seidel switching
    <Graph.seidel_switching>` class of the disjoint union of a 1-vertex graph
    and the one produced by :func:`Pseudo-L_{2n}(4n-1)
    <sage.graphs.graph_generators.GraphGenerators.SquaredSkewHadamardMatrixGraph>`

    In this case, the other possible parameter set of a strongly regular graph
    in the Seidel switching class of the latter graph (see [BH2012]_) coincides
    with the set of parameters of the complement of the graph returned by this
    function.

    .. SEEALSO::

        - :func:`~sage.graphs.strongly_regular_db.is_switch_skewhad`

    EXAMPLES::

        sage: g = graphs.SwitchedSquaredSkewHadamardMatrixGraph(4)                      # needs sage.modules
        sage: g.is_strongly_regular(parameters=True)                                    # needs sage.modules
        (226, 105, 48, 49)
        sage: from sage.combinat.designs.twographs import twograph_descendant
        sage: twograph_descendant(g, 0).is_strongly_regular(parameters=True)            # needs sage.modules
        (225, 112, 55, 56)
        sage: gc = g.complement()                                                       # needs sage.modules
        sage: twograph_descendant(gc, 0).is_strongly_regular(parameters=True)           # needs sage.modules
        (225, 112, 55, 56)

    TESTS::

        sage: graphs.SwitchedSquaredSkewHadamardMatrixGraph(0)
        Traceback (most recent call last):
        ...
        ValueError: parameter n must be >= 1
    """
    G = SquaredSkewHadamardMatrixGraph(n).complement()
    G.add_vertex((4 * n - 1)**2)
    G.seidel_switching(list(range((4 * n - 1) * (2 * n - 1))))
    G.name("switch skewhad^2+*_" + str((n)))
    return G


def HanoiTowerGraph(pegs, disks, labels=True, positions=True):
    r"""
    Returns the graph whose vertices are the states of the
    Tower of Hanoi puzzle, with edges representing legal moves between states.

    INPUT:

    - ``pegs`` - the number of pegs in the puzzle, 2 or greater
    - ``disks`` - the number of disks in the puzzle, 1 or greater
    - ``labels`` - default: ``True``, if ``True`` the graph contains
      more meaningful labels, see explanation below.  For large instances,
      turn off labels for much faster creation of the graph.
    - ``positions`` - default: ``True``, if ``True`` the graph contains
      layout information.  This creates a planar layout for the case
      of three pegs.  For large instances, turn off layout information
      for much faster creation of the graph.

    OUTPUT:

    The Tower of Hanoi puzzle has a certain number of identical pegs
    and a certain number of disks, each of a different radius.
    Initially the disks are all on a single peg, arranged
    in order of their radii, with the largest on the bottom.

    The goal of the puzzle is to move the disks to any other peg,
    arranged in the same order.  The one constraint is that the
    disks resident on any one peg must always be arranged with larger
    radii lower down.

    The vertices of this graph represent all the possible states
    of this puzzle.  Each state of the puzzle is a tuple with length
    equal to the number of disks, ordered by largest disk first.
    The entry of the tuple is the peg where that disk resides.
    Since disks on a given peg must go down in size as we go
    up the peg, this totally describes the state of the puzzle.

    For example ``(2,0,0)`` means the large disk is on peg 2, the
    medium disk is on peg 0, and the small disk is on peg 0
    (and we know the small disk must be above the medium disk).
    We encode these tuples as integers with a base equal to
    the number of pegs, and low-order digits to the right.

    Two vertices are adjacent if we can change the puzzle from
    one state to the other by moving a single disk.  For example,
    ``(2,0,0)`` is adjacent to ``(2,0,1)`` since we can move
    the small disk off peg 0 and onto (the empty) peg 1.
    So the solution to a 3-disk puzzle (with at least
    two pegs) can be expressed by the shortest path between
    ``(0,0,0)`` and ``(1,1,1)``.  For more on this representation
    of the graph, or its properties, see [AD2010]_.

    For greatest speed we create graphs with integer vertices,
    where we encode the tuples as integers with a base equal
    to the number of pegs, and low-order digits to the right.
    So for example, in a 3-peg puzzle with 5 disks, the
    state ``(1,2,0,1,1)`` is encoded as
    `1\ast 3^4 + 2\ast 3^3 + 0\ast 3^2 + 1\ast 3^1 + 1\ast 3^0 = 139`.

    For smaller graphs, the labels that are the tuples are informative,
    but slow down creation of the graph.  Likewise computing layout
    information also incurs a significant speed penalty. For maximum
    speed, turn off labels and layout and decode the
    vertices explicitly as needed.  The
    :meth:`sage.rings.integer.Integer.digits`
    with the ``padsto`` option is a quick way to do this, though you
    may want to reverse the list that is output.

    .. SEEALSO::

        - :meth:`~sage.graphs.generators.families.GeneralizedSierpinskiGraph`

    PLOTTING:

    The layout computed when ``positions = True`` will
    look especially good for the three-peg case, when the graph is known
    to be planar.  Except for two small cases on 4 pegs, the graph is
    otherwise not planar, and likely there is a better way to layout
    the vertices.

    EXAMPLES:

    A classic puzzle uses 3 pegs.  We solve the 5 disk puzzle using
    integer labels and report the minimum number of moves required.
    Note that `3^5-1` is the state where all 5 disks
    are on peg 2. ::

        sage: H = graphs.HanoiTowerGraph(3, 5, labels=False, positions=False)
        sage: H.distance(0, 3^5-1)
        31

    A slightly larger instance. ::

        sage: H = graphs.HanoiTowerGraph(4, 6, labels=False, positions=False)
        sage: H.num_verts()
        4096
        sage: H.distance(0, 4^6-1)
        17

    For a small graph, labels and layout information can be useful.
    Here we explicitly list a solution as a list of states. ::

        sage: H = graphs.HanoiTowerGraph(3, 3, labels=True, positions=True)
        sage: H.shortest_path((0,0,0), (1,1,1))
        [(0, 0, 0), (0, 0, 1), (0, 2, 1), (0, 2, 2), (1, 2, 2), (1, 2, 0), (1, 1, 0), (1, 1, 1)]

    Some facts about this graph with `p` pegs and `d` disks:

    - only automorphisms are the "obvious" ones - renumber the pegs.
    - chromatic number is less than or equal to `p`
    - independence number is `p^{d-1}`

    ::

        sage: H = graphs.HanoiTowerGraph(3, 4, labels=False, positions=False)
        sage: H.automorphism_group().is_isomorphic(SymmetricGroup(3))                   # needs sage.groups
        True
        sage: H.chromatic_number()
        3
        sage: len(H.independent_set()) == 3^(4-1)
        True

    TESTS:

    It is an error to have just one peg (or less). ::

        sage: graphs.HanoiTowerGraph(1, 5)
        Traceback (most recent call last):
        ...
        ValueError: Pegs for Tower of Hanoi graph should be two or greater (not 1)

    It is an error to have zero disks (or less). ::

        sage: graphs.HanoiTowerGraph(2, 0)
        Traceback (most recent call last):
        ...
        ValueError: Disks for Tower of Hanoi graph should be one or greater (not 0)

    AUTHOR:

    - Rob Beezer, (2009-12-26), with assistance from Su Doree

    """
    # sanitize input
    from sage.rings.integer import Integer
    pegs = Integer(pegs)
    if pegs < 2:
        raise ValueError("Pegs for Tower of Hanoi graph should be two or greater (not %d)" % pegs)
    disks = Integer(disks)
    if disks < 1:
        raise ValueError("Disks for Tower of Hanoi graph should be one or greater (not %d)" % disks)

    # Each state of the puzzle is a tuple with length
    # equal to the number of disks, ordered by largest disk first
    # The entry of the tuple is the peg where that disk resides
    # Since disks on a given peg must go down in size as we go
    # up the peg, this totally describes the puzzle
    # We encode these tuples as integers with a base equal to
    # the number of pegs, and low-order digits to the right

    # complete graph on number of pegs when just a single disk
    edges = [[i, j] for i in range(pegs) for j in range(i + 1, pegs)]

    nverts = 1
    for d in range(2, disks+1):
        prevedges = edges      # remember subgraph to build from
        nverts = pegs*nverts   # pegs^(d-1)
        edges = []

        # Take an edge, change its two states in the same way by adding
        # a large disk to the bottom of the same peg in each state
        # This is accomplished by adding a multiple of pegs^(d-1)
        for p in range(pegs):
            largedisk = p*nverts
            for anedge in prevedges:
                edges.append([anedge[0] + largedisk, anedge[1] + largedisk])

        # Two new states may only differ in the large disk
        # being the only disk on two different pegs, thus
        # otherwise being a common state with one less disk
        # We construct all such pairs of new states and add as edges
        from sage.combinat.subset import Subsets
        for state in range(nverts):
            emptypegs = list(range(pegs))
            reduced_state = state
            for i in range(d-1):
                apeg = reduced_state % pegs
                if apeg in emptypegs:
                    emptypegs.remove(apeg)
                reduced_state = reduced_state//pegs
            for freea, freeb in Subsets(emptypegs, 2):
                edges.append([freea*nverts + state, freeb*nverts + state])

    H = Graph({}, loops=False, multiedges=False)
    H.add_edges(edges)

    # Making labels and/or computing positions can take a long time,
    # relative to just constructing the edges on integer vertices.
    # We try to minimize coercion overhead, but need Sage
    # Integers in order to use digits() for labels.
    # Getting the digits with custom code was no faster.
    # Layouts are circular (symmetric on the number of pegs)
    # radiating outward to the number of disks (radius)
    # Algorithm uses some combination of alternate
    # clockwise/counterclockwise placements, which
    # works well for three pegs (planar layout)
    #
    if labels or positions:
        mapping = {}
        pos = {}
        a = Integer(-1)
        one = Integer(1)
        if positions:
            radius_multiplier = 1 + 1/sin(pi/pegs)
            sine = []
            cosine = []
            for i in range(pegs):
                angle = 2*i*pi/float(pegs)
                sine.append(sin(angle))
                cosine.append(cos(angle))
        for i in range(pegs**disks):
            a += one
            state = a.digits(base=pegs, padto=disks)
            if labels:
                state.reverse()
                mapping[i] = tuple(state)
                state.reverse()
            if positions:
                locx = 0.0
                locy = 0.0
                radius = 1.0
                parity = -1.0
                for index in range(disks):
                    p = state[index]
                    radius *= radius_multiplier
                    parity *= -1.0
                    locx_temp = cosine[p]*locx - parity*sine[p]*locy + radius*cosine[p]
                    locy_temp = parity*sine[p]*locx + cosine[p]*locy - radius*parity*sine[p]
                    locx = locx_temp
                    locy = locy_temp
                pos[i] = (locx, locy)
        # set positions, then relabel (not vice versa)
        if positions:
            H.set_pos(pos)
        if labels:
            H.relabel(mapping)

    return H


def line_graph_forbidden_subgraphs():
    r"""
    Returns the 9 forbidden subgraphs of a line graph.

    See the :wikipedia:`Line_graph` for more information.

    The graphs are returned in the ordering given by the Wikipedia
    drawing, read from left to right and from top to bottom.

    EXAMPLES::

        sage: graphs.line_graph_forbidden_subgraphs()
        [Claw graph: Graph on 4 vertices,
        Graph on 6 vertices,
        Graph on 6 vertices,
        Graph on 5 vertices,
        Graph on 6 vertices,
        Graph on 6 vertices,
        Graph on 6 vertices,
        Graph on 6 vertices,
        Graph on 5 vertices]

    """
    from sage.graphs.graph import Graph
    from sage.graphs.generators.basic import ClawGraph
    graphs = [ClawGraph()]

    graphs.append(Graph({
                0: [1, 2, 3],
                1: [2, 3],
                4: [2],
                5: [3]
                }))

    graphs.append(Graph({
                0: [1, 2, 3, 4],
                1: [2, 3, 4],
                3: [4],
                2: [5]
                }))

    graphs.append(Graph({
                0: [1, 2, 3],
                1: [2, 3],
                4: [2, 3]
                }))

    graphs.append(Graph({
                0: [1, 2, 3],
                1: [2, 3],
                4: [2],
                5: [3, 4]
                }))

    graphs.append(Graph({
                0: [1, 2, 3, 4],
                1: [2, 3, 4],
                3: [4],
                5: [2, 0, 1]
                }))

    graphs.append(Graph({
                5: [0, 1, 2, 3, 4],
                0: [1, 4],
                2: [1, 3],
                3: [4]
                }))

    graphs.append(Graph({
                1: [0, 2, 3, 4],
                3: [0, 4],
                2: [4, 5],
                4: [5]
                }))

    graphs.append(Graph({
                0: [1, 2, 3],
                1: [2, 3, 4],
                2: [3, 4],
                3: [4]
                }))

    return graphs


def petersen_family(generate=False):
    r"""
    Returns the Petersen family

    The Petersen family is a collection of 7 graphs which are the forbidden
    minors of the linklessly embeddable graphs. For more information see the
    :wikipedia:`Petersen_family`.

    INPUT:

    - ``generate`` (boolean) -- whether to generate the family from the
      `\Delta-Y` transformations. When set to ``False`` (default) a hardcoded
      version of the graphs (with a prettier layout) is returned.

    EXAMPLES::

        sage: graphs.petersen_family()
        [Petersen graph: Graph on 10 vertices,
         Complete graph: Graph on 6 vertices,
         Multipartite Graph with set sizes [3, 3, 1]: Graph on 7 vertices,
         Graph on 8 vertices,
         Graph on 9 vertices,
         Graph on 7 vertices,
         Graph on 8 vertices]

    The two different inputs generate the same graphs::

        sage: F1 = graphs.petersen_family(generate=False)
        sage: F2 = graphs.petersen_family(generate=True)                                # needs sage.modules
        sage: F1 = [g.canonical_label().graph6_string() for g in F1]
        sage: F2 = [g.canonical_label().graph6_string() for g in F2]                    # needs sage.modules
        sage: set(F1) == set(F2)                                                        # needs sage.modules
        True
    """
    from sage.graphs.generators.smallgraphs import PetersenGraph
    if not generate:
        from sage.graphs.generators.basic import CompleteGraph, \
             CompleteBipartiteGraph, CompleteMultipartiteGraph
        l = [PetersenGraph(), CompleteGraph(6),
             CompleteMultipartiteGraph([3, 3, 1])]
        g = CompleteBipartiteGraph(4, 4)
        g.delete_edge(0, 4)
        g.name("")
        l.append(g)
        g = Graph('HKN?Yeb')
        g._circle_embedding([1, 2, 4, 3, 0, 5])
        g._circle_embedding([6, 7, 8], radius=.6, shift=1.25)
        l.append(g)
        g = Graph('Fs\\zw')
        g._circle_embedding([1, 2, 3])
        g._circle_embedding([4, 5, 6], radius=.7)
        g._pos[0] = (0, 0)
        l.append(g)
        g = Graph('GYQ[p{')
        g._circle_embedding([1, 4, 6, 0, 5, 7, 3], shift=0.25)
        g._pos[2] = (0, 0)
        l.append(g)
        return l

    def DeltaYTrans(G, triangle):
        """
        Apply a Delta-Y transformation to a given triangle of G.
        """
        a, b, c = triangle
        G = G.copy()
        G.delete_edges([(a, b), (b, c), (c, a)])
        v = G.order()
        G.add_edges([(a, v), (b, v), (c, v)])
        return G.canonical_label()

    def YDeltaTrans(G, v):
        """
        Apply a Y-Delta transformation to a given vertex v of G.
        """
        G = G.copy()
        a, b, c = G.neighbors(v)
        G.delete_vertex(v)
        G.add_cycle([a, b, c])
        return G.canonical_label()

    # We start from the Petersen Graph, and apply Y-Delta transform
    # for as long as we generate new graphs.
    P = PetersenGraph()

    l = set([])
    l_new = [P.canonical_label().graph6_string()]

    while l_new:
        g = l_new.pop(0)
        if g in l:
            continue
        l.add(g)
        g = Graph(g)
        # All possible Delta-Y transforms
        for t in g.subgraph_search_iterator(Graph({1: [2, 3], 2: [3]}), return_graphs=False):
            l_new.append(DeltaYTrans(g, t).graph6_string())
        # All possible Y-Delta transforms
        for v in g:
            if g.degree(v) == 3:
                l_new.append(YDeltaTrans(g, v).graph6_string())

    return [Graph(x) for x in l]


def SierpinskiGasketGraph(n):
    """
    Return the Sierpinski Gasket graph of generation `n`.

    All vertices but 3 have valence 4.

    INPUT:

    - `n` -- an integer

    OUTPUT:

    a graph `S_n` with `3 (3^{n-1}+1)/2` vertices and
    `3^n` edges, closely related to the famous Sierpinski triangle
    fractal.

    All these graphs have a triangular shape, and three special
    vertices at top, bottom left and bottom right. These are the only
    vertices of valence 2, all the other ones having valence 4.

    The graph `S_1` (generation `1`) is a triangle.

    The graph `S_{n+1}` is obtained from the disjoint union of
    three copies A,B,C of `S_n` by identifying pairs of vertices:
    the top vertex of A with the bottom left vertex of B,
    the bottom right vertex of B with the top vertex of C,
    and the bottom left vertex of C with the bottom right vertex of A.

    .. PLOT::

        sphinx_plot(graphs.SierpinskiGasketGraph(4).plot(vertex_labels=False))

    .. SEEALSO::

        - :meth:`~sage.graphs.generators.families.HanoiTowerGraph`. There is
          another family of graphs called Sierpinski graphs, where all vertices
          but 3 have valence 3. They are available using
          ``graphs.HanoiTowerGraph(3, n)``.
        - :meth:`~sage.graphs.generators.families.GeneralizedSierpinskiGraph`

    EXAMPLES::

        sage: # needs sage.modules
        sage: s4 = graphs.SierpinskiGasketGraph(4); s4
        Graph on 42 vertices
        sage: s4.size()
        81
        sage: s4.degree_histogram()
        [0, 0, 3, 0, 39]
        sage: s4.is_hamiltonian()
        True

    REFERENCES:

    [LLWC2011]_
    """
    from sage.modules.free_module_element import vector
    from sage.rings.rational_field import QQ

    if n <= 0:
        raise ValueError('n should be at least 1')

    def next_step(triangle_list):
        # compute the next subdivision
        resu = []
        for a, b, c in triangle_list:
            ab = (a + b) / 2
            bc = (b + c) / 2
            ac = (a + c) / 2
            resu += [(a, ab, ac), (ab, b, bc), (ac, bc, c)]
        return resu

    tri_list = [list(vector(QQ, u) for u in [(0, 0), (0, 1), (1, 0)])]
    for k in range(n - 1):
        tri_list = next_step(tri_list)
    dg = Graph()
    dg.add_edges([(tuple(a), tuple(b)) for a, b, c in tri_list])
    dg.add_edges([(tuple(b), tuple(c)) for a, b, c in tri_list])
    dg.add_edges([(tuple(c), tuple(a)) for a, b, c in tri_list])
    dg.set_pos({(x, y): (x + y / 2, y * 3 / 4)
                for (x, y) in dg})
    dg.relabel()
    return dg


def GeneralizedSierpinskiGraph(G, k, stretch=None):
    r"""
    Return the generalized Sierpinski graph of `G` of dimension `k`.

    Generalized Sierpinski graphs have been introduced in [GKP2011]_ to
    generalize the notion of Sierpinski graphs [KM1997]_.

    Given a graph `G = (V, E)` of order `n` and a parameter `k`, the generalized
    Sierpinski graph of `G` of dimension `k`, denoted by `S(G, k)`, can be
    constructed recursively from `G` as follows. `S(G, 1)` is isomorphic to
    `G`. To construct `S(G, k)` for `k > 1`, copy `n` times `S(G, k - 1)`, once
    per vertex `u \in V`, and add `u` at the beginning of the labels of each
    vertex in the copy of `S(G, k - 1)` corresponding to vertex `u`. Then for
    any edge `\{u, v\} \in E`, add an edge between vertex `(u, v, \ldots, v)`
    and vertex `(v, u, \ldots, u)`.

    INPUT:

    - ``G`` -- a sage Graph

    - ``k`` -- integer; the dimension

    - ``stretch`` -- integer (default: ``None``); stretching factor used to
      determine the positions of the vertices of the output graph. By default
      (``None``), this value is set to twice the maximum Euclidian distance
      between the vertices of `G`. This parameter is used only when the vertices
      of `G` have positions.

    .. SEEALSO::

        - :meth:`~sage.graphs.generators.families.SierpinskiGasketGraph`
        - :meth:`~sage.graphs.generators.families.HanoiTowerGraph`

    EXAMPLES:

    The generalized Sierpinski graph of dimension 1 of any graph `G`
    is isomorphic to `G`::

        sage: G = graphs.RandomGNP(10, .5)
        sage: S = graphs.GeneralizedSierpinskiGraph(G, 1)
        sage: S.is_isomorphic(G)
        True

    When `G` is a clique of order 3, the generalized Sierpinski graphs
    of `G` are isomorphic to Hanoi Tower graphs::

        sage: k = randint(1, 5)
        sage: S = graphs.GeneralizedSierpinskiGraph(graphs.CompleteGraph(3), k)         # needs sage.modules
        sage: H = graphs.HanoiTowerGraph(3, k)
        sage: S.is_isomorphic(H)                                                        # needs sage.modules
        True

    The generalized Sierpinski graph of dimension `k` of any graph `G` with `n`
    vertices and `m` edges has `n^k` vertices and `m\sum_{i=0}^{k-1}n^i` edges::

        sage: # needs sage.modules
        sage: n = randint(2, 6)
        sage: k = randint(1, 5)
        sage: G = graphs.RandomGNP(n, .5)
        sage: m = G.size()
        sage: S = graphs.GeneralizedSierpinskiGraph(G, k)
        sage: S.order() == n**k
        True
        sage: S.size() == m*sum([n**i for i in range(k)])
        True
        sage: G = graphs.CompleteGraph(n)
        sage: S = graphs.GeneralizedSierpinskiGraph(G, k)
        sage: S.order() == n**k
        True
        sage: S.size() == (n*(n - 1)/2)*sum([n**i for i in range(k)])
        True

    The positions of the vertices of the output graph are determined from the
    positions of the vertices of `G`, if any::

        sage: G = graphs.HouseGraph()
        sage: G.get_pos() is not None
        True
        sage: H = graphs.GeneralizedSierpinskiGraph(G, 2)                               # needs sage.symbolic
        sage: H.get_pos() is not None                                                   # needs sage.symbolic
        True
        sage: G = Graph([(0, 1)])
        sage: G.get_pos() is not None
        False
        sage: H = graphs.GeneralizedSierpinskiGraph(G, 2)                               # needs sage.symbolic
        sage: H.get_pos() is not None                                                   # needs sage.symbolic
        False

    .. PLOT::

        sphinx_plot(graphs.GeneralizedSierpinskiGraph(graphs.HouseGraph(), 2).plot(vertex_labels=False))

    TESTS::

        sage: # needs sage.modules
        sage: graphs.GeneralizedSierpinskiGraph(Graph(), 3)
        Generalized Sierpinski Graph of Graph on 0 vertices of dimension 3: Graph on 0 vertices
        sage: graphs.GeneralizedSierpinskiGraph(Graph(1), 3).vertices(sort=False)
        [(0, 0, 0)]
        sage: G = graphs.GeneralizedSierpinskiGraph(Graph(2), 3)
        sage: G.order(), G.size()
        (8, 0)
        sage: graphs.GeneralizedSierpinskiGraph("foo", 1)
        Traceback (most recent call last):
        ...
        ValueError: parameter G must be a Graph
        sage: graphs.GeneralizedSierpinskiGraph(Graph(), 0)
        Traceback (most recent call last):
        ...
        ValueError: parameter k must be >= 1
    """
    if not isinstance(G, Graph):
        raise ValueError("parameter G must be a Graph")
    if k < 1:
        raise ValueError("parameter k must be >= 1")
    loops = G.allows_loops()
    multiedges = G.allows_multiple_edges()

    def rec(H, kk):
        if kk == 1:
            return H
        I = Graph(loops=loops, multiedges=multiedges)
        # add one copy of H per vertex of G
        for i in G:
            J = H.relabel(perm={u: (i,) + u for u in H}, inplace=False)
            I.add_vertices(J)
            I.add_edges(J.edge_iterator(labels=False, sort_vertices=False))
        # For each edge {u, v} of G, add edge {(u, v, ..., v), (v, u, ..., u)}
        l = len(next(H.vertex_iterator()))
        for u, v in G.edges(sort=True, labels=False):
            I.add_edge((u,) + (v,)*l, (v,) + (u,)*l)
        return rec(I, kk - 1)

    H = G.relabel(perm={u: (u,) for u in G}, inplace=False)
    if H and k > 1:
        H = rec(H, k)
    H.name("Generalized Sierpinski Graph of {} of dimension {}".format(G, k))

    # If the vertices of G have positions, we set the positions of vertices of H
    pos = G.get_pos()
    if pos:
        if stretch is None:
            # Find the geometric diameter
            from sage.modules.free_module_element import vector
            L = [vector(p) for p in pos.values()]
            stretch = 2 * max((u - v).norm() for u, v in combinations(L, 2))

        H.set_pos({u: (sum(pos[x][0]*stretch**(k-i) for i, x in enumerate(u)),
                       sum(pos[y][1]*stretch**(k-i) for i, y in enumerate(u)))
                   for u in H})
    return H


def WheelGraph(n):
    """
    Returns a Wheel graph with n nodes.

    A Wheel graph is a basic structure where one node is connected to all other
    nodes and those (outer) nodes are connected cyclically.

    PLOTTING: Upon construction, the position dictionary is filled to override
    the spring-layout algorithm. By convention, each wheel graph will be
    displayed with the first (0) node in the center, the second node at the top,
    and the rest following in a counterclockwise manner.

    With the wheel graph, we see that it doesn't take a very large n at all for
    the spring-layout to give a counter-intuitive display. (See Graphics Array
    examples below).

    EXAMPLES:

    We view many wheel graphs with a Sage Graphics Array, first with this
    constructor (i.e., the position dictionary filled)::

        sage: # needs sage.plot
        sage: g = []
        sage: j = []
        sage: for i in range(9):
        ....:  k = graphs.WheelGraph(i+3)
        ....:  g.append(k)
        ...
        sage: for i in range(3):
        ....:  n = []
        ....:  for m in range(3):
        ....:      n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ....:  j.append(n)
        ...
        sage: G = graphics_array(j)
        sage: G.show()                          # long time

    Next, using the spring-layout algorithm::

        sage: # needs networkx sage.plot
        sage: import networkx
        sage: g = []
        sage: j = []
        sage: for i in range(9):
        ....:  spr = networkx.wheel_graph(i+3)
        ....:  k = Graph(spr)
        ....:  g.append(k)
        ...
        sage: for i in range(3):
        ....:  n = []
        ....:  for m in range(3):
        ....:      n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ....:  j.append(n)
        ...
        sage: G = graphics_array(j)
        sage: G.show()                          # long time

    Compare the plotting::

        sage: # needs networkx sage.plot
        sage: n = networkx.wheel_graph(23)
        sage: spring23 = Graph(n)
        sage: posdict23 = graphs.WheelGraph(23)
        sage: spring23.show()                   # long time
        sage: posdict23.show()  # long time
    """
    from sage.graphs.generators.basic import CycleGraph
    if n < 4:
        G = CycleGraph(n)
    else:
        G = CycleGraph(n-1)
        G.relabel(perm=list(range(1, n)), inplace=True)
        G.add_edges([(0, i) for i in range(1, n)])
        G._pos[0] = (0, 0)
    G.name("Wheel graph")
    return G


def WindmillGraph(k, n):
    r"""
    Return the Windmill graph `Wd(k, n)`.

    The windmill graph `Wd(k, n)` is an undirected graph constructed for `k \geq
    2` and `n \geq 2` by joining `n` copies of the complete graph `K_k` at a
    shared vertex. It has `(k-1)n+1` vertices and `nk(k-1)/2` edges, girth 3 (if
    `k > 2`), radius 1 and diameter 2. It has vertex connectivity 1 because its
    central vertex is an articulation point; however, like the complete graphs
    from which it is formed, it is `(k-1)`-edge-connected. It is trivially
    perfect and a block graph.

    .. SEEALSO::

        - :wikipedia:`Windmill_graph`
        - :meth:`GraphGenerators.StarGraph`
        - :meth:`GraphGenerators.FriendshipGraph`

    EXAMPLES:

    The Windmill graph `Wd(2, n)` is a star graph::

        sage: n = 5
        sage: W = graphs.WindmillGraph(2, n)
        sage: W.is_isomorphic( graphs.StarGraph(n) )
        True

    The Windmill graph `Wd(3, n)` is the Friendship graph `F_n`::

        sage: n = 5
        sage: W = graphs.WindmillGraph(3, n)
        sage: W.is_isomorphic( graphs.FriendshipGraph(n) )
        True

    The Windmill graph `Wd(3, 2)` is the Butterfly graph::

        sage: W = graphs.WindmillGraph(3, 2)
        sage: W.is_isomorphic( graphs.ButterflyGraph() )
        True

    The Windmill graph `Wd(k, n)` has chromatic number `k`::

        sage: n,k = 5,6
        sage: W = graphs.WindmillGraph(k, n)
        sage: W.chromatic_number() == k
        True

    TESTS:

    Giving too small parameters::

        sage: graphs.WindmillGraph(1, 2)
        Traceback (most recent call last):
        ...
        ValueError: parameters k and n must be >= 2
        sage: graphs.WindmillGraph(2, 1)
        Traceback (most recent call last):
        ...
        ValueError: parameters k and n must be >= 2
    """
    if k < 2 or n < 2:
        raise ValueError('parameters k and n must be >= 2')

    if k == 2:
        from sage.graphs.generators.basic import StarGraph
        G = StarGraph(n)
    else:
        sector = 2*pi/n
        slide = 1/sin(sector/4)

        pos_dict = {}
        for i in range(0, k):
            x = float(cos(i*pi/(k-2)))
            y = float(sin(i*pi/(k-2))) + slide
            pos_dict[i] = (x, y)

        G = Graph()
        pos = {0: [0, 0]}
        for i in range(n):
            V = list(range(i*(k - 1) + 1, (i + 1)*(k - 1) + 1))
            G.add_clique([0]+V)
            for j, v in enumerate(V):
                x, y = pos_dict[j]
                xv = x*cos(i*sector) - y*sin(i*sector)
                yv = x*sin(i*sector) + y*cos(i*sector)
                pos[v] = [xv, yv]

        G.set_pos(pos)

    G.name("Windmill graph Wd({}, {})".format(k, n))
    return G


def trees(vertices):
    r"""
    Returns a generator of the distinct trees on a fixed number of vertices.

    INPUT:

    -  ``vertices`` - the size of the trees created.

    OUTPUT:

    A generator which creates an exhaustive, duplicate-free listing
    of the connected free (unlabeled) trees with ``vertices`` number
    of vertices.  A tree is a graph with no cycles.

    ALGORITHM:

    Uses an algorithm that generates each new tree
    in constant time.  See the documentation for, and implementation
    of, the :mod:`sage.graphs.trees` module, including a citation.

    EXAMPLES:

    We create an iterator, then loop over its elements. ::

        sage: tree_iterator = graphs.trees(7)
        sage: for T in tree_iterator:
        ....:     print(T.degree_sequence())
        [2, 2, 2, 2, 2, 1, 1]
        [3, 2, 2, 2, 1, 1, 1]
        [3, 2, 2, 2, 1, 1, 1]
        [4, 2, 2, 1, 1, 1, 1]
        [3, 3, 2, 1, 1, 1, 1]
        [3, 3, 2, 1, 1, 1, 1]
        [4, 3, 1, 1, 1, 1, 1]
        [3, 2, 2, 2, 1, 1, 1]
        [4, 2, 2, 1, 1, 1, 1]
        [5, 2, 1, 1, 1, 1, 1]
        [6, 1, 1, 1, 1, 1, 1]

    The number of trees on the first few vertex counts.
    This is sequence A000055 in Sloane's OEIS. ::

        sage: [len(list(graphs.trees(i))) for i in range(0, 15)]
        [1, 1, 1, 1, 2, 3, 6, 11, 23, 47, 106, 235, 551, 1301, 3159]
    """
    from sage.graphs.trees import TreeIterator
    return iter(TreeIterator(vertices))


def nauty_gentreeg(options="", debug=False):
    r"""
    Return a generator which creates non-isomorphic trees from nauty's gentreeg
    program.

    INPUT:

    - ``options`` -- string (default: ``""``); a string passed to ``gentreeg``
      as if it was run at a system command line. At a minimum, you *must* pass
      the number of vertices you desire. Sage expects the graphs to be in
      nauty's "sparse6" format, do not set an option to change this default or
      results will be unpredictable.

    - ``debug`` -- boolean (default: ``False``); if ``True`` the first line of
      ``gentreeg``'s output to standard error is captured and the first call to
      the generator's ``next()`` function will return this line as a string. A
      line leading with ">A" indicates a successful initiation of the program
      with some information on the arguments, while a line beginning with ">E"
      indicates an error with the input.

    The possible options, obtained as output of ``gentreeg -help``::

           n            : the number of vertices. Must be in range 1..128
        res/mod         : only generate subset res out of subsets 0..mod-1
          -D<int>       : an upper bound for the maximum degree
          -Z<int>:<int> : bounds on the diameter
          -q            : suppress auxiliary output

    Options which cause ``gentreeg`` to use an output format different than the
    sparse6 format are not listed above (-p, -l, -u) as they will confuse the
    creation of a Sage graph. The res/mod option can be useful when using the
    output in a routine run several times in parallel.

    OUTPUT:

    A generator which will produce the graphs as Sage graphs. These will be
    simple graphs: no loops, no multiple edges, no directed edges.

    .. SEEALSO::

        :meth:`trees` -- another generator of trees

    EXAMPLES:

    The generator can be used to construct trees for testing, one at a time
    (usually inside a loop). Or it can be used to create an entire list all at
    once if there is sufficient memory to contain it::

        sage: gen = graphs.nauty_gentreeg("4")
        sage: next(gen)
        Graph on 4 vertices
        sage: next(gen)
        Graph on 4 vertices
        sage: next(gen)
        Traceback (most recent call last):
        ...
        StopIteration

    The number of trees on the first few vertex counts. This agrees with
    :oeis:`A000055`::

        sage: [len(list(graphs.nauty_gentreeg(str(i)))) for i in range(1, 15)]
        [1, 1, 1, 2, 3, 6, 11, 23, 47, 106, 235, 551, 1301, 3159]

    The ``debug`` switch can be used to examine ``gentreeg``'s reaction to the
    input in the ``options`` string.  We illustrate success. (A failure will be
    a string beginning with ">E".)  Passing the "-q" switch to ``gentreeg`` will
    suppress the indicator of a successful initiation, and so the first returned
    value might be an empty string if ``debug`` is ``True``::

        sage: gen = graphs.nauty_gentreeg("4", debug=True)
        sage: print(next(gen))
        >A ...gentreeg ...
        sage: gen = graphs.nauty_gentreeg("4 -q", debug=True)
        sage: next(gen)
        ''

    TESTS:

    The number `n` of vertices must be in range 1..128::

        sage: list(graphs.nauty_gentreeg("0", debug=False))
        Traceback (most recent call last):
        ...
        ValueError: wrong format of parameter options
        sage: list(graphs.nauty_gentreeg("0", debug=True))
        ['>E gentreeg: n must be in the range 1..128\n']
        sage: list(graphs.nauty_gentreeg("200", debug=True))
        ['>E gentreeg: n must be in the range 1..128\n']

    Wrong input::

        sage: list(graphs.nauty_gentreeg("3 -x", debug=False))
        Traceback (most recent call last):
        ...
        ValueError: wrong format of parameter options
        sage: list(graphs.nauty_gentreeg("3 -x", debug=True))
        ['>E Usage: ...gentreeg [-D#] [-Z#:#] [-ulps] [-q] n [res/mod] ...
        sage: list(graphs.nauty_gentreeg("3", debug=True))
        ['>A ...gentreeg ...\n', Graph on 3 vertices]
    """
    import shlex
    from sage.features.nauty import NautyExecutable
    gen_path = NautyExecutable("gentreeg").absolute_filename()
    sp = subprocess.Popen(shlex.quote(gen_path) + " {0}".format(options), shell=True,
                          stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, close_fds=True,
                          encoding='latin-1')
    msg = sp.stderr.readline()
    if debug:
        yield msg
    elif msg.startswith('>E'):
        raise ValueError('wrong format of parameter options')
    gen = sp.stdout
    while True:
        try:
            s = next(gen)
        except StopIteration:
            # Exhausted list of graphs from nauty geng
            return
        G = Graph(s[:-1], format='sparse6', loops=False, multiedges=False)
        yield G


def RingedTree(k, vertex_labels=True):
    r"""
    Return the ringed tree on k-levels.

    A ringed tree of level `k` is a binary tree with `k` levels (counting
    the root as a level), in which all vertices at the same level are connected
    by a ring.

    More precisely, in each layer of the binary tree (i.e. a layer is the set of
    vertices `[2^i...2^{i+1}-1]`) two vertices `u,v` are adjacent if `u=v+1` or
    if `u=2^i` and `v=`2^{i+1}-1`.

    Ringed trees are defined in [CFHM2013]_.

    INPUT:

    - ``k`` -- the number of levels of the ringed tree.

    - ``vertex_labels`` (boolean) -- whether to label vertices as binary words
      (default) or as integers.

    EXAMPLES::

        sage: # needs networkx
        sage: G = graphs.RingedTree(5)
        sage: P = G.plot(vertex_labels=False, vertex_size=10)                           # needs sage.plot
        sage: P.show()                          # long time                             # needs sage.plot
        sage: G.vertices(sort=True)
        ['', '0', '00', '000', '0000', '0001', '001', '0010', '0011', '01',
         '010', '0100', '0101', '011', '0110', '0111', '1', '10', '100',
         '1000', '1001', '101', '1010', '1011', '11', '110', '1100', '1101',
         '111', '1110', '1111']

    TESTS::

        sage: G = graphs.RingedTree(-1)
        Traceback (most recent call last):
        ...
        ValueError: The number of levels must be >= 1.
        sage: G = graphs.RingedTree(5, vertex_labels=False)                             # needs networkx
        sage: G.vertices(sort=True)                                                     # needs networkx
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
        18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]
    """
    if k < 1:
        raise ValueError('The number of levels must be >= 1.')

    # Creating the Balanced tree, which contains most edges already
    g = BalancedTree(2, k - 1)
    g.name('Ringed Tree on ' + str(k) + ' levels')

    # We consider edges layer by layer
    for i in range(1, k):
        vertices = list(range(2**(i) - 1, 2**(i + 1) - 1))

        # Add the missing edges
        g.add_cycle(vertices)

        # And set the vertices' positions
        radius = i if i <= 1 else 1.5**i
        shift = -2**(i - 2) + .5 if i > 1 else 0
        g._circle_embedding(vertices, radius=radius, shift=shift)

    # Specific position for the central vertex
    g._pos[0] = (0, 0.2)

    # Relabel vertices as binary words
    if not vertex_labels:
        return g

    vertices = ['']
    for i in range(k - 1):
        for j in range(2**(i) - 1, 2**(i + 1) - 1):
            v = vertices[j]
            vertices.append(v + '0')
            vertices.append(v + '1')

    g.relabel(vertices)

    return g


def MathonPseudocyclicMergingGraph(M, t):
    r"""
    Mathon's merging of classes in a pseudo-cyclic 3-class association scheme

    Construct strongly regular graphs from p.97 of [BL1984]_.

    INPUT:

    - ``M`` -- the list of matrices in a pseudo-cyclic 3-class association scheme.
      The identity matrix must be the first entry.

    - ``t`` (integer) -- the number of the graph, from 0 to 2.

    .. SEEALSO::

        - :func:`~sage.graphs.strongly_regular_db.is_muzychuk_S6`

    TESTS::

        sage: from sage.graphs.generators.families import MathonPseudocyclicMergingGraph as mer
        sage: from sage.graphs.generators.smallgraphs import _EllipticLinesProjectivePlaneScheme as ES
        sage: G = mer(ES(3), 0)  # long time
        sage: G.is_strongly_regular(parameters=True)    # long time
        (784, 243, 82, 72)
        sage: G = mer(ES(3), 1)  # long time
        sage: G.is_strongly_regular(parameters=True)    # long time
        (784, 270, 98, 90)
        sage: G = mer(ES(3), 2)  # long time
        sage: G.is_strongly_regular(parameters=True)    # long time
        (784, 297, 116, 110)
        sage: G = mer(ES(2), 2)                                                         # needs sage.libs.gap
        Traceback (most recent call last):
        ...
        AssertionError...
        sage: M = ES(3)                                                                 # needs sage.libs.gap
        sage: M = [M[1],M[0],M[2],M[3]]                                                 # needs sage.libs.gap
        sage: G = mer(M, 2)                                                             # needs sage.libs.gap
        Traceback (most recent call last):
        ...
        AssertionError...
    """
    from sage.graphs.graph import Graph
    from sage.matrix.constructor import identity_matrix
    assert len(M) == 4
    assert M[0] == identity_matrix(M[0].nrows())
    A = sum(x.tensor_product(x) for x in M[1:])
    if t > 0:
        A += sum(x.tensor_product(M[0]) for x in M[1:])
    if t > 1:
        A += sum(M[0].tensor_product(x) for x in M[1:])
    return Graph(A)


def MathonPseudocyclicStronglyRegularGraph(t, G=None, L=None):
    r"""
    Return a strongly regular graph on `(4t+1)(4t-1)^2` vertices from
    [Mat1978]_.

    Let `4t-1` be a prime power, and `4t+1` be such that there exists
    a strongly regular graph `G` with parameters `(4t+1,2t,t-1,t)`. In
    particular, `4t+1` must be a sum of two squares [Mat1978]_. With
    this input, Mathon [Mat1978]_ gives a construction of a strongly regular
    graph with parameters `(4 \mu + 1, 2 \mu, \mu-1, \mu)`, where
    `\mu =  t(4t(4t-1)-1)`. The construction is optionally parametrised by an
    a skew-symmetric Latin square of order `4t+1`, with entries in
    `-2t,...,-1,0,1,...,2t`.

    Our implementation follows a description given in [ST1981]_.

    INPUT:

    - ``t`` -- a positive integer

    - ``G`` -- if ``None`` (default), try to construct the necessary graph
      with parameters `(4t+1,2t,t-1,t)`, otherwise use the user-supplied one,
      with vertices labelled from `0` to `4t`.

    - ``L`` -- if ``None`` (default), construct a necessary skew Latin square,
      otherwise use the user-supplied one. Here non-isomorphic Latin squares
      -- one constructed from `Z/9Z`, and the other from `(Z/3Z)^2` --
      lead to non-isomorphic graphs.

    .. SEEALSO::

        - :func:`~sage.graphs.strongly_regular_db.is_mathon_PC_srg`

    EXAMPLES:

    Using default ``G`` and ``L``. ::

        sage: from sage.graphs.generators.families import MathonPseudocyclicStronglyRegularGraph
        sage: G = MathonPseudocyclicStronglyRegularGraph(1); G                          # needs sage.modules sage.rings.finite_rings
        Mathon's PC SRG on 45 vertices: Graph on 45 vertices
        sage: G.is_strongly_regular(parameters=True)                                    # needs sage.modules sage.rings.finite_rings
        (45, 22, 10, 11)

    Supplying ``G`` and ``L`` (constructed from the automorphism group of ``G``). ::

        sage: # needs sage.groups sage.libs.gap sage.rings.finite_rings
        sage: G = graphs.PaleyGraph(9)
        sage: a = G.automorphism_group(partition=[sorted(G)])
        sage: it = (x for x in a.normal_subgroups() if x.order() == 9)
        sage: subg = next(iter(it))
        sage: r = [matrix(libgap.PermutationMat(libgap(z), 9).sage())
        ....:      for z in subg]
        sage: ff = list(map(lambda y: (y[0]-1,y[1]-1),
        ....:          Permutation(map(lambda x: 1+r.index(x^-1), r)).cycle_tuples()[1:]))
        sage: L = sum(i*(r[a]-r[b]) for i,(a,b) in zip(range(1,len(ff)+1), ff)); L
        [ 0  1 -1 -3 -2 -4  3  4  2]
        [-1  0  1 -4 -3 -2  2  3  4]
        [ 1 -1  0 -2 -4 -3  4  2  3]
        [ 3  4  2  0  1 -1 -3 -2 -4]
        [ 2  3  4 -1  0  1 -4 -3 -2]
        [ 4  2  3  1 -1  0 -2 -4 -3]
        [-3 -2 -4  3  4  2  0  1 -1]
        [-4 -3 -2  2  3  4 -1  0  1]
        [-2 -4 -3  4  2  3  1 -1  0]

        sage: # needs sage.modules sage.rings.finite_rings sage.groups sage.libs.gap
        sage: G.relabel(range(9))
        sage: G3x3 = graphs.MathonPseudocyclicStronglyRegularGraph(2, G=G, L=L)
        sage: G3x3.is_strongly_regular(parameters=True)
        (441, 220, 109, 110)
        sage: G3x3.automorphism_group(algorithm="bliss").order()  # optional - bliss
        27
        sage: G9 = graphs.MathonPseudocyclicStronglyRegularGraph(2)
        sage: G9.is_strongly_regular(parameters=True)
        (441, 220, 109, 110)
        sage: G9.automorphism_group(algorithm="bliss").order()  # optional - bliss
        9

    TESTS::

        sage: graphs.MathonPseudocyclicStronglyRegularGraph(5)                          # needs sage.modules
        Traceback (most recent call last):
        ...
        ValueError: 21  must be a sum of two squares!...
    """
    from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
    from sage.rings.integer_ring import ZZ
    from sage.matrix.constructor import matrix, block_matrix, \
        ones_matrix, identity_matrix
    from sage.arith.misc import two_squares
    p = 4*t + 1
    try:
        x = two_squares(p)
    except ValueError:
        raise ValueError(str(p)+" must be a sum of two squares!")
    if G is None:
        from sage.graphs.strongly_regular_db import strongly_regular_graph as SRG
        G = SRG(p, 2*t, t - 1)
        G.relabel(range(p))
    if L is None:
        from sage.matrix.constructor import circulant
        L = circulant(list(range(2 * t + 1)) + list(range(-2 * t, 0)))
    q = 4*t - 1
    K = GF(q, prefix='x')
    K_pairs = set(frozenset([x, -x]) for x in K)
    K_pairs.discard(frozenset([0]))
    a = [None]*(q-1)    # order the non-0 elements of K as required
    for i, (x, y) in enumerate(K_pairs):
        a[i] = x
        a[-i-1] = y
    a.append(K(0))      # and append the 0 of K at the end
    P = [matrix(ZZ, q, q, lambda i, j: 1 if a[j] == a[i] + b else 0)
         for b in a]
    g = K.primitive_element()
    F = sum(P[a.index(g**(2*i))] for i in range(1, 2*t))
    E = matrix(ZZ, q, q, lambda i, j: 0 if (a[j] - a[0]).is_square() else 1)

    def B(m):
        I = identity_matrix(q)
        J = ones_matrix(q)
        if m == 0:
            def f(i, j):
                if i == j:
                    return 0 * I
                elif (a[j] - a[i]).is_square():
                    return I + F
                else:
                    return J - F
        elif m < 2*t:
            def f(i, j):
                return F * P[a.index(g**(2*m) * (a[i] + a[j]))]
        elif m == 2*t:
            def f(i, j):
                return E * P[i]
        return block_matrix(q, q, [f(i, j) for i in range(q) for j in range(q)])

    def Acon(i, j):
        J = ones_matrix(q**2)
        if i == j:
            return B(0)
        if L[i, j] > 0:
            if G.has_edge(i, j):
                return B(L[i, j])
            return J - B(L[i, j])
        if G.has_edge(i, j):
            return B(-L[i, j]).T
        return J - B(-L[i, j]).T

    A = Graph(block_matrix(p, p, [Acon(i, j) for i in range(p) for j in range(p)]))
    A.name("Mathon's PC SRG on " + str(p*q**2) + " vertices")
    A.relabel()
    return A


def TuranGraph(n, r):
    r"""
    Returns the Turan graph with parameters `n, r`.

    Turan graphs are complete multipartite graphs with `n` vertices and `r`
    subsets, denoted `T(n,r)`, with the property that the sizes of the subsets
    are as close to equal as possible. The graph `T(n,r)` will have `n \pmod r`
    subsets of size `\lfloor n/r \rfloor` and `r - (n \pmod r)` subsets of size
    `\lceil n/r \rceil`. See the :wikipedia:`Turan_graph` for more information.

    INPUT:

    - ``n`` -- integer; the number of vertices in the graph

    - ``r`` -- integer; the number of partitions of the graph

    EXAMPLES:

    The Turan graph is a complete multipartite graph::

        sage: g = graphs.TuranGraph(13, 4)
        sage: k = graphs.CompleteMultipartiteGraph([3,3,3,4])
        sage: g.is_isomorphic(k)
        True

    The Turan graph `T(n,r)` has `\frac{(r-1)(n^2-s^2)}{2r} + \frac{s(s-1)}{2}`
    edges, where `s = n \mod r` (:trac:`34249`)::

        sage: n = 12
        sage: r = 8
        sage: g = graphs.TuranGraph(n, r)
        sage: def count(n, r):
        ....:     s = n % r
        ....:     return (r - 1) * (n**2 - s**2) / (2*r) + s*(s - 1)/2
        sage: g.size() == count(n, r)
        True
        sage: n = randint(3, 100)
        sage: r = randint(2, n - 1)
        sage: g = graphs.TuranGraph(n, r)
        sage: g.size() == count(n, r)
        True

    TESTS::

        sage: g = graphs.TuranGraph(3,6)
        Traceback (most recent call last):
        ...
        ValueError: input parameters must satisfy "1 < r < n"
    """
    if n < 1 or n < r or r < 1:
        raise ValueError('input parameters must satisfy "1 < r < n"')

    from sage.graphs.generators.basic import CompleteMultipartiteGraph

    p = n // r
    s = n % r
    vertex_sets = [p]*(r - s) + [p + 1]*s

    g = CompleteMultipartiteGraph(vertex_sets)
    g.name('Turan Graph with n: {}, r: {}'.format(n, r))

    return g


def MuzychukS6Graph(n, d, Phi='fixed', Sigma='fixed', verbose=False):
    r"""
    Return a strongly regular graph of S6 type from [Muz2007]_ on
    `n^d((n^d-1)/(n-1)+1)` vertices.

    The construction depends upon a number of parameters, two of them, `n` and
    `d`, mandatory, and `\Phi` and `\Sigma` mappings defined in [Muz2007]_.
    These graphs have parameters `(mn^d, n^{d-1}(m-1) - 1,\mu - 2,\mu)`, where
    `\mu=\frac{n^{d-1}-1}{n-1}n^{d-1}` and `m:=\frac{n^d-1}{n-1}+1`.

    Some details on `\Phi` and `\Sigma` are as follows.  Let `L` be the
    complete graph on `M:=\{0,..., m-1\}` with the matching
    `\{(2i,2i+1) | i=0,...,m/2\}` removed.
    Then one arbitrarily chooses injections `\Phi_i`
    from the edges of `L` on `i \in M` into sets of parallel classes of affine
    `d`-dimensional designs; our implementation uses the designs of hyperplanes
    in `d`-dimensional affine geometries over `GF(n)`. Finally, for each edge
    `ij` of `L` one arbitrarily chooses bijections `\Sigma_{ij}` between
    `\Phi_i` and `\Phi_j`. More details, in particular how these choices lead
    to non-isomorphic graphs, are in [Muz2007]_.

    INPUT:

    - ``n`` (integer)-- a prime power

    - ``d`` (integer)-- must be odd if `n` is odd

    - ``Phi`` is an optional parameter of the construction; it must be either

        - 'fixed'-- this will generate fixed default `\Phi_i`, for `i \in M`, or

        - 'random'-- `\Phi_i` are generated at random, or

        - A dictionary describing the functions `\Phi_i`; for `i \in M`,
          Phi[(i, T)] in `M`, for each edge T of `L` on `i`.
          Also, each `\Phi_i` must be injective.

    - ``Sigma`` is an optional parameter of the construction; it must be either

        - 'fixed'-- this will generate a fixed default `\Sigma`, or

        - 'random'-- `\Sigma` is generated at random.

    - ``verbose`` (Boolean)-- default is False. If True, print progress information

    .. SEEALSO::

        - :func:`~sage.graphs.strongly_regular_db.is_muzychuk_S6`

    .. TODO::

        Implement the possibility to explicitly supply the parameter `\Sigma`
        of the construction.

    EXAMPLES::

        sage: graphs.MuzychukS6Graph(3, 3).is_strongly_regular(parameters=True)         # needs sage.modules sage.rings.finite_rings
        (378, 116, 34, 36)
        sage: phi = {(2,(0,2)):0, (1,(1,3)):1, (0,(0,3)):1, (2,(1,2)):1,                # needs sage.modules
        ....:        (1,(1,2)):0, (0,(0,2)):0, (3,(0,3)):0, (3,(1,3)):1}
        sage: graphs.MuzychukS6Graph(2, 2,                                              # needs sage.modules sage.rings.finite_rings
        ....:                        Phi=phi).is_strongly_regular(parameters=True)
        (16, 5, 0, 2)

    TESTS::

        sage: # needs sage.modules
        sage: graphs.MuzychukS6Graph(2,2,Phi='random',Sigma='random').is_strongly_regular(parameters=True)              # needs sage.rings.finite_rings
        (16, 5, 0, 2)
        sage: graphs.MuzychukS6Graph(3,3,Phi='random',Sigma='random').is_strongly_regular(parameters=True)              # needs sage.rings.finite_rings
        (378, 116, 34, 36)
        sage: graphs.MuzychukS6Graph(3,2)
        Traceback (most recent call last):
        ...
        AssertionError: n must be even or d must be odd
        sage: graphs.MuzychukS6Graph(6,2)
        Traceback (most recent call last):
        ...
        AssertionError: n must be a prime power
        sage: graphs.MuzychukS6Graph(3,1)
        Traceback (most recent call last):
        ...
        AssertionError: d must be at least 2
        sage: graphs.MuzychukS6Graph(3,3,Phi=42)                                        # needs sage.rings.finite_rings
        Traceback (most recent call last):
        ...
        AssertionError: Phi must be a dictionary or 'random' or 'fixed'
        sage: graphs.MuzychukS6Graph(3,3,Sigma=42)                                      # needs sage.rings.finite_rings
        Traceback (most recent call last):
        ...
        ValueError: Sigma must be 'random' or 'fixed'
    """
    # TO DO: optimise
    #        add option to return phi, sigma? generate phi, sigma from seed? (int say?)

    from sage.combinat.designs.block_design import ProjectiveGeometryDesign
    from sage.misc.prandom import randrange
    from sage.misc.functional import is_even
    from sage.arith.misc import is_prime_power
    from sage.graphs.generators.basic import CompleteGraph
    from sage.rings.finite_rings.finite_field_constructor import GF
    from sage.matrix.special import ones_matrix
    from sage.matrix.constructor import matrix
    from sage.rings.rational_field import QQ
    from sage.rings.integer_ring import ZZ
    from time import time

    assert d > 1,              'd must be at least 2'
    assert is_even(n * (d-1)), 'n must be even or d must be odd'
    assert is_prime_power(n),  'n must be a prime power'
    t = time()

    # build L, L_i and the design
    m = int((n**d - 1)/(n - 1) + 1)  # from m = p + 1, p = (n^d-1) / (n-1)
    L = CompleteGraph(m)
    L.delete_edges([(2 * x, 2 * x + 1) for x in range(m // 2)])
    L_i = [L.edges_incident(x, labels=False) for x in range(m)]
    Design = ProjectiveGeometryDesign(d, d-1, GF(n, 'a'), point_coordinates=False)
    projBlocks = Design.blocks()
    atInf = projBlocks[-1]
    Blocks = [[x for x in block if x not in atInf] for block in projBlocks[:-1]]
    if verbose:
        print('finished preamble at %f (+%f)' % (time() - t, time() - t))
    t1 = time()

    # sort the hyperplanes into parallel classes
    ParClasses = [Blocks]
    while ParClasses[0]:
        nextHyp = ParClasses[0].pop()
        for C in ParClasses[1:]:
            listC = sum(C, [])
            for x in nextHyp:
                if x in listC:
                    break
            else:
                C.append(nextHyp)
                break
        else:
            ParClasses.append([nextHyp])
    del ParClasses[0]
    if verbose:
        print('finished ParClasses at %f (+%f)' % (time() - t, time() - t1))
    t1 = time()

    # build E^C_j
    E = {}
    v = ZZ(n**d)
    k = ZZ(n**(d-1))
    ones = ones_matrix(v)
    ones_v = ones/v
    for C in ParClasses:
        EC = matrix(QQ, v)
        for line in C:
            for i, j in combinations(line, 2):
                EC[i, j] = EC[j, i] = 1/k
        EC -= ones_v
        E[tuple(C[0])] = EC
    if verbose:
        print('finished E at %f (+%f)' % (time() - t, time() - t1))
    t1 = time()

    # handle Phi
    if Phi == 'random':
        Phi = {}
        for x in range(m):
            temp = list(range(len(ParClasses)))
            for line in L_i[x]:
                rand = randrange(0, len(temp))
                Phi[(x, line)] = temp.pop(rand)
    elif Phi == 'fixed':
        Phi = {(x, line): val for x in range(m) for val, line in enumerate(L_i[x])}
    else:
        assert isinstance(Phi, dict), \
               "Phi must be a dictionary or 'random' or 'fixed'"
        assert set(Phi.keys()) == set([(x, line) for x in range(m) for line in L_i[x]]), \
               'each Phi_i must have domain L_i'
        for x in range(m):
            assert m - 2 == len(set([val for (key, val) in Phi.items() if key[0] == x])), \
                   'each phi_i must be injective'
        for val in Phi.values():
            assert val in range(m - 1), \
                   'codomain should be {0,..., (n^d - 1)/(n - 1) - 1}'
    phi = {(x, line): ParClasses[Phi[(x, line)]] for x in range(m) for line in L_i[x]}
    if verbose:
        print('finished phi at %f (+%f)' % (time() - t, time() - t1))
    t1 = time()

    # handle sigma
    sigma = {}
    if Sigma == 'random':
        for x in range(m):
            for line in L_i[x]:
                [i, j] = line
                temp = phi[(j, line)][:]
                for hyp in phi[(i, line)]:
                    rand = randrange(0, len(temp))
                    sigma[(i, j, tuple(hyp))] = temp[rand]
                    sigma[(j, i, tuple(temp[rand]))] = hyp
                    del temp[rand]
    elif Sigma == 'fixed':
        for x in range(m):
            for line in L_i[x]:
                [i, j] = line
                temp = phi[(j, line)][:]
                for hyp in phi[(i, line)]:
                    val = temp.pop()
                    sigma[(i, j, tuple(hyp))] = val
                    sigma[(j, i, tuple(val))] = hyp
    else:
        raise ValueError("Sigma must be 'random' or 'fixed'")
    if verbose:
        print('finished sigma at %f (+%f)' % (time() - t, time() - t1))
    t1 = time()

    # build V
    edges = []  # how many? *m^2*n^2
    for (i, j) in L.edges(sort=True, labels=False):
        for hyp in phi[(i, (i, j))]:
            for x in hyp:
                newEdges = [((i, x), (j, y))
                            for y in sigma[(i, j, tuple(hyp))]]
                edges.extend(newEdges)
    if verbose:
        print('finished edges at %f (+%f)' % (time() - t, time() - t1))
    t1 = time()
    V = Graph(edges)
    if verbose:
        print('finished V at %f (+%f)' % (time() - t, time() - t1))
    t1 = time()

    # build D_i, F_i and A_i
    D_i = [0]*m
    for x in range(m):
        D_i[x] = sum([E[tuple(phi[x, line][0])] for line in L_i[x]])
    F_i = [1 - D_i[x] - ones_v for x in range(m)]
    # as the sum of (1/v)*J_\Omega_i, D_i, F_i is identity
    A_i = [(v-k)*ones_v - k*F_i[x] for x in range(m)]
    #   we know A_i = k''*(1/v)*J_\Omega_i + r''*D_i + s''*F_i,
    #   and (k'', s'', r'') = (v - k, 0, -k)
    if verbose:
        print('finished D, F and A at %f (+%f)' % (time() - t, time() - t1))
    t1 = time()

    # add the edges of the graph of B to V
    for i in range(m):
        V.add_edges([((i, x), (i, y)) for x in range(v)
                     for y in range(v) if not A_i[i][(x, y)]])

    V.name('Muzychuk S6 graph with parameters ('+str(n)+','+str(d)+')')
    if verbose:
        print('finished at %f (+%f)' % ((time() - t), time() - t1))
    return V


def CubeConnectedCycle(d):
    r"""
    Return the cube-connected cycle of dimension `d`.

    The cube-connected cycle of order `d` is the `d`-dimensional hypercube
    with each of its vertices replaced by a cycle of length `d`. This graph has
    order `d \times 2^d`.
    The construction is as follows:
    Construct vertex `(x,y)` for `0 \leq x < 2^d`, `0 \leq y < d`.
    For each vertex, `(x,y)`, add an edge between it and `(x, (y-1) \mod d))`,
    `(x,(y+1) \mod d)`, and `(x \oplus 2^y, y)`, where `\oplus` is the bitwise
    xor operator.

    For `d=1` and `2`, the cube-connected cycle graph contains self-loops or
    multiple edges between a pair of vertices, but for all other `d`, it is
    simple.

    INPUT:

    - ``d`` -- The dimension of the desired hypercube as well as the length
      of the cycle to be placed at each vertex of the `d`-dimensional
      hypercube. `d` must be a positive integer.

    EXAMPLES:

    The order of the graph is `d \times 2^d` ::

        sage: d = 3
        sage: g = graphs.CubeConnectedCycle(d)
        sage: len(g) == d*2**d
        True

    The diameter of cube-connected cycles for `d > 3` is
    `2d + \lfloor \frac{d}{2} \rfloor - 2` ::

        sage: d = 4
        sage: g = graphs.CubeConnectedCycle(d)
        sage: g.diameter() == 2*d+d//2-2
        True

    All vertices have degree `3` when `d > 1` ::

        sage: g = graphs.CubeConnectedCycle(5)
        sage: all(g.degree(v) == 3 for v in g)
        True

    TESTS::

        sage: g = graphs.CubeConnectedCycle(0)
        Traceback (most recent call last):
        ...
        ValueError: the dimension d must be greater than 0
    """
    if d < 1:
        raise ValueError('the dimension d must be greater than 0')

    G = Graph(name="Cube-Connected Cycle of dimension {}".format(d))

    if d == 1:
        G.allow_loops(True)
        # only d = 1 requires loops
        G.add_edges([((0, 0), (0, 1)), ((0, 0), (0, 0)), ((0, 1), (0, 1))])
        return G

    if d == 2:
        # only d = 2 require multiple edges
        G.allow_multiple_edges(True)
        G.add_edges([((0, 0), (0, 1)), ((0, 0), (0, 1)), ((0, 0), (1, 0)),
                     ((0, 1), (2, 1)), ((1, 0), (1, 1)), ((1, 0), (1, 1)),
                     ((1, 1), (3, 1)), ((2, 0), (2, 1)), ((2, 0), (2, 1)),
                     ((2, 0), (3, 0)), ((3, 0), (3, 1)), ((3, 0), (3, 1))])
        return G

    for x in range(1 << d):
        G.add_cycle([(x, y) for y in range(d)])

    for x, y in G:
        G.add_edge((x, y), (x ^ (1 << y), y))

    return G
