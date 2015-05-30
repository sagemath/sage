# -*- coding: utf-8 -*-
r"""
Random Graphs

The methods defined here appear in :mod:`sage.graphs.graph_generators`.
"""
###########################################################################
#
#           Copyright (C) 2006 Robert L. Miller <rlmillster@gmail.com>
#                              and Emily A. Kirkman
#           Copyright (C) 2009 Michael C. Yurko <myurko@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
###########################################################################

# import from Sage library
from sage.graphs.graph import Graph
from sage.misc.randstate import current_randstate

def RandomGNP(n, p, seed=None, fast=True, method='Sage'):
    r"""
    Returns a random graph on `n` nodes. Each edge is inserted independently
    with probability `p`.

    INPUTS:

    - ``n`` -- number of nodes of the graph

    - ``p`` -- probability of an edge

    - ``seed`` -- integer seed for random number generator (default=None).

    - ``fast`` -- boolean set to True (default) to use the algorithm with
      time complexity in `O(n+m)` proposed in [BatBra2005]_. It is designed
      for generating large sparse graphs. It is faster than other methods for
      *LARGE* instances (try it to know whether it is useful for you).

    - ``method`` -- By default (```method='Sage'``), this function uses the
      method implemented in ```sage.graphs.graph_generators_pyx.pyx``. When
      ``method='networkx'``, this function calls the NetworkX function
      ``fast_gnp_random_graph``, unless ``fast=False``, then
      ``gnp_random_graph``. Try them to know which method is the best for
      you. The ``fast`` parameter is not taken into account by the 'Sage'
      method so far.

    REFERENCES:

    .. [ErdRen1959] P. Erdos and A. Renyi. On Random Graphs, Publ.
       Math. 6, 290 (1959).

    .. [Gilbert1959] E. N. Gilbert. Random Graphs, Ann. Math. Stat.,
       30, 1141 (1959).

    .. [BatBra2005] V. Batagelj and U. Brandes. Efficient generation of
       large random networks. Phys. Rev. E, 71, 036113, 2005.

    PLOTTING: When plotting, this graph will use the default spring-layout
    algorithm, unless a position dictionary is specified.

    EXAMPLES: We show the edge list of a random graph on 6 nodes with
    probability `p = .4`::

        sage: set_random_seed(0)
        sage: graphs.RandomGNP(6, .4).edges(labels=False)
        [(0, 1), (0, 5), (1, 2), (2, 4), (3, 4), (3, 5), (4, 5)]

    We plot a random graph on 12 nodes with probability `p = .71`::

        sage: gnp = graphs.RandomGNP(12,.71)
        sage: gnp.show() # long time

    We view many random graphs using a graphics array::

        sage: g = []
        sage: j = []
        sage: for i in range(9):
        ....:     k = graphs.RandomGNP(i+3,.43)
        ....:     g.append(k)
        sage: for i in range(3):
        ....:     n = []
        ....:     for m in range(3):
        ....:         n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ....:     j.append(n)
        sage: G = sage.plot.graphics.GraphicsArray(j)
        sage: G.show() # long time
        sage: graphs.RandomGNP(4,1)
        Complete graph: Graph on 4 vertices

    TESTS::

        sage: graphs.RandomGNP(50,.2,method=50)
        Traceback (most recent call last):
        ...
        ValueError: 'method' must be equal to 'networkx' or to 'Sage'.
        sage: set_random_seed(0)
        sage: graphs.RandomGNP(50,.2, method="Sage").size()
        243
        sage: graphs.RandomGNP(50,.2, method="networkx").size()
        258
    """
    if n < 0:
        raise ValueError("The number of nodes must be positive or null.")
    if 0.0 > p or 1.0 < p:
        raise ValueError("The probability p must be in [0..1].")

    if seed is None:
        seed = current_randstate().long_seed()
    if p == 1:
        from sage.graphs.generators.basic import CompleteGraph
        return CompleteGraph(n)

    if method == 'networkx':
        import networkx
        if fast:
            G = networkx.fast_gnp_random_graph(n, p, seed=seed)
        else:
            G = networkx.gnp_random_graph(n, p, seed=seed)
        return Graph(G)
    elif method in ['Sage', 'sage']:
        # We use the Sage generator
        from sage.graphs.graph_generators_pyx import RandomGNP as sageGNP
        return sageGNP(n, p)
    else:
        raise ValueError("'method' must be equal to 'networkx' or to 'Sage'.")

def RandomBarabasiAlbert(n, m, seed=None):
    u"""
    Return a random graph created using the Barabasi-Albert preferential
    attachment model.

    A graph with m vertices and no edges is initialized, and a graph of n
    vertices is grown by attaching new vertices each with m edges that are
    attached to existing vertices, preferentially with high degree.

    INPUT:

    - ``n`` - number of vertices in the graph

    - ``m`` - number of edges to attach from each new node

    - ``seed`` - for random number generator

    EXAMPLES:

    We show the edge list of a random graph on 6 nodes with m = 2.

    ::

        sage: graphs.RandomBarabasiAlbert(6,2).edges(labels=False)
        [(0, 2), (0, 3), (0, 4), (1, 2), (2, 3), (2, 4), (2, 5), (3, 5)]

    We plot a random graph on 12 nodes with m = 3.

    ::

        sage: ba = graphs.RandomBarabasiAlbert(12,3)
        sage: ba.show()  # long time

    We view many random graphs using a graphics array::

        sage: g = []
        sage: j = []
        sage: for i in range(1,10):
        ....:     k = graphs.RandomBarabasiAlbert(i+3, 3)
        ....:     g.append(k)
        sage: for i in range(3):
        ....:     n = []
        ....:     for m in range(3):
        ....:         n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ....:     j.append(n)
        sage: G = sage.plot.graphics.GraphicsArray(j)
        sage: G.show()  # long time

    """
    if seed is None:
        seed = current_randstate().long_seed()
    import networkx
    return Graph(networkx.barabasi_albert_graph(n,m,seed=seed))

def RandomBipartite(n1,n2, p):
    r"""
    Returns a bipartite graph with `n1+n2` vertices
    such that any edge from `[n1]` to `[n2]` exists
    with probability `p`.

    INPUT:

        - ``n1,n2`` : Cardinalities of the two sets
        - ``p``   : Probability for an edge to exist


    EXAMPLE::

        sage: g=graphs.RandomBipartite(5,2,0.5)
        sage: g.vertices()
        [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (1, 0), (1, 1)]

    TESTS::

        sage: g=graphs.RandomBipartite(5,-3,0.5)
        Traceback (most recent call last):
        ...
        ValueError: n1 and n2 should be integers strictly greater than 0
        sage: g=graphs.RandomBipartite(5,3,1.5)
        Traceback (most recent call last):
        ...
        ValueError: Parameter p is a probability, and so should be a real value between 0 and 1

    Trac ticket #12155::

        sage: graphs.RandomBipartite(5,6,.2).complement()
        complement(Random bipartite graph of size 5+6 with edge probability 0.200000000000000): Graph on 11 vertices
    """
    if not (p>=0 and p<=1):
        raise ValueError("Parameter p is a probability, and so should be a real value between 0 and 1")
    if not (n1>0 and n2>0):
        raise ValueError("n1 and n2 should be integers strictly greater than 0")

    from numpy.random import uniform

    g=Graph(name="Random bipartite graph of size "+str(n1) +"+"+str(n2)+" with edge probability "+str(p))

    S1=[(0,i) for i in range(n1)]
    S2=[(1,i) for i in range(n2)]
    g.add_vertices(S1)
    g.add_vertices(S2)

    for w in range(n2):
        for v in range(n1):
            if uniform()<=p :
                g.add_edge((0,v),(1,w))

    pos = {}
    for i in range(n1):
        pos[(0,i)] = (0, i/(n1-1.0))
    for i in range(n2):
        pos[(1,i)] = (1, i/(n2-1.0))

    g.set_pos(pos)

    return g

def RandomBoundedToleranceGraph(n):
    r"""
    Returns a random bounded tolerance graph.

    The random tolerance graph is built from a random bounded
    tolerance representation by using the function
    `ToleranceGraph`. This representation is a list
    `((l_0,r_0,t_0), (l_1,r_1,t_1), ..., (l_k,r_k,t_k))` where
    `k = n-1` and `I_i = (l_i,r_i)` denotes a random interval and
    `t_i` a random positive value less then or equal to the length
    of the interval `I_i`. The width of the representation is
    limited to n**2 * 2**n.

    .. NOTE::

        The tolerance representation used to create the graph can
        be recovered using ``get_vertex()`` or ``get_vertices()``.

    INPUT:

    - ``n`` -- number of vertices of the random graph.

    EXAMPLE:

    Every (bounded) tolerance graph is perfect. Hence, the
    chromatic number is equal to the clique number ::

        sage: g = graphs.RandomBoundedToleranceGraph(8)
        sage: g.clique_number() == g.chromatic_number()
        True
    """
    from sage.misc.prandom import randint
    from sage.graphs.generators.intersection import ToleranceGraph

    W = n ** 2 * 2 ** n

    tolrep = map(lambda l_r: (l_r[0], l_r[1], randint(0, l_r[1] - l_r[0])),
                 [sorted((randint(0, W), randint(0, W))) for i in range(n)])

    return ToleranceGraph(tolrep)

def RandomGNM(n, m, dense=False, seed=None):
    """
    Returns a graph randomly picked out of all graphs on n vertices
    with m edges.

    INPUT:

    -  ``n`` - number of vertices.

    -  ``m`` - number of edges.

    -  ``dense`` - whether to use NetworkX's
       dense_gnm_random_graph or gnm_random_graph


    EXAMPLES: We show the edge list of a random graph on 5 nodes with
    10 edges.

    ::

        sage: graphs.RandomGNM(5, 10).edges(labels=False)
        [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]

    We plot a random graph on 12 nodes with m = 12.

    ::

        sage: gnm = graphs.RandomGNM(12, 12)
        sage: gnm.show()  # long time

    We view many random graphs using a graphics array::

        sage: g = []
        sage: j = []
        sage: for i in range(9):
        ....:     k = graphs.RandomGNM(i+3, i^2-i)
        ....:     g.append(k)
        sage: for i in range(3):
        ....:     n = []
        ....:     for m in range(3):
        ....:         n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ....:     j.append(n)
        sage: G = sage.plot.graphics.GraphicsArray(j)
        sage: G.show()  # long time
    """
    if seed is None:
        seed = current_randstate().long_seed()
    import networkx
    if dense:
        return Graph(networkx.dense_gnm_random_graph(n, m, seed=seed))
    else:
        return Graph(networkx.gnm_random_graph(n, m, seed=seed))

def RandomNewmanWattsStrogatz(n, k, p, seed=None):
    """
    Returns a Newman-Watts-Strogatz small world random graph on n
    vertices.

    From the NetworkX documentation: First create a ring over n nodes.
    Then each node in the ring is connected with its k nearest
    neighbors. Then shortcuts are created by adding new edges as
    follows: for each edge u-v in the underlying "n-ring with k nearest
    neighbors"; with probability p add a new edge u-w with
    randomly-chosen existing node w. In contrast with
    watts_strogatz_graph(), no edges are removed.

    INPUT:

    -  ``n`` - number of vertices.

    -  ``k`` - each vertex is connected to its k nearest
       neighbors

    -  ``p`` - the probability of adding a new edge for
       each edge

    -  ``seed`` - for the random number generator


    EXAMPLE: We show the edge list of a random graph on 7 nodes with 2
    "nearest neighbors" and probability `p = 0.2`::

        sage: graphs.RandomNewmanWattsStrogatz(7, 2, 0.2).edges(labels=False)
        [(0, 1), (0, 2), (0, 3), (0, 6), (1, 2), (2, 3), (2, 4), (3, 4), (3, 6), (4, 5), (5, 6)]

    ::

        sage: G = graphs.RandomNewmanWattsStrogatz(12, 2, .3)
        sage: G.show()  # long time

    REFERENCE:

    .. [NWS99] Newman, M.E.J., Watts, D.J. and Strogatz, S.H.  Random
      graph models of social networks. Proc. Nat. Acad. Sci. USA
      99, 2566-2572.
    """
    if seed is None:
        seed = current_randstate().long_seed()
    import networkx
    return Graph(networkx.newman_watts_strogatz_graph(n, k, p, seed=seed))

def RandomHolmeKim(n, m, p, seed=None):
    """
    Returns a random graph generated by the Holme and Kim algorithm for
    graphs with power law degree distribution and approximate average
    clustering.

    INPUT:

    -  ``n`` - number of vertices.

    -  ``m`` - number of random edges to add for each new
       node.

    -  ``p`` - probability of adding a triangle after
       adding a random edge.

    -  ``seed`` - for the random number generator.


    From the NetworkX documentation: The average clustering has a hard
    time getting above a certain cutoff that depends on m. This cutoff
    is often quite low. Note that the transitivity (fraction of
    triangles to possible triangles) seems to go down with network
    size. It is essentially the Barabasi-Albert growth model with an
    extra step that each random edge is followed by a chance of making
    an edge to one of its neighbors too (and thus a triangle). This
    algorithm improves on B-A in the sense that it enables a higher
    average clustering to be attained if desired. It seems possible to
    have a disconnected graph with this algorithm since the initial m
    nodes may not be all linked to a new node on the first iteration
    like the BA model.

    EXAMPLE: We show the edge list of a random graph on 8 nodes with 2
    random edges per node and a probability `p = 0.5` of
    forming triangles.

    ::

        sage: graphs.RandomHolmeKim(8, 2, 0.5).edges(labels=False)
        [(0, 2), (0, 5), (1, 2), (1, 3), (2, 3), (2, 4), (2, 6), (2, 7),
         (3, 4), (3, 6), (3, 7), (4, 5)]

    ::

        sage: G = graphs.RandomHolmeKim(12, 3, .3)
        sage: G.show()  # long time

    REFERENCE:

    .. [HolmeKim2002] Holme, P. and Kim, B.J. Growing scale-free networks
      with tunable clustering, Phys. Rev. E (2002). vol 65, no 2, 026107.
    """
    if seed is None:
        seed = current_randstate().long_seed()
    import networkx
    return Graph(networkx.powerlaw_cluster_graph(n, m, p, seed=seed))

def RandomIntervalGraph(n):
    """
    Returns a random interval graph.

    An interval graph is built from a list `(a_i,b_i)_{1\leq i \leq n}`
    of intervals : to each interval of the list is associated one
    vertex, two vertices being adjacent if the two corresponding
    intervals intersect.

    A random interval graph of order `n` is generated by picking
    random values for the `(a_i,b_j)`, each of the two coordinates
    being generated from the uniform distribution on the interval
    `[0,1]`.

    This definitions follows [boucheron2001]_.

    .. NOTE::

        The vertices are named 0, 1, 2, and so on. The intervals
        used to create the graph are saved with the graph and can
        be recovered using ``get_vertex()`` or ``get_vertices()``.

    INPUT:

    - ``n`` (integer) -- the number of vertices in the random
      graph.

    EXAMPLE:

    As for any interval graph, the chromatic number is equal to
    the clique number ::

        sage: g = graphs.RandomIntervalGraph(8)
        sage: g.clique_number() == g.chromatic_number()
        True

    REFERENCE:

    .. [boucheron2001] Boucheron, S. and FERNANDEZ de la VEGA, W.,
       On the Independence Number of Random Interval Graphs,
       Combinatorics, Probability and Computing v10, issue 05,
       Pages 385--396,
       Cambridge Univ Press, 2001
    """

    from sage.misc.prandom import random
    from sage.graphs.generators.intersection import IntervalGraph

    intervals = [tuple(sorted((random(), random()))) for i in range(n)]
    return IntervalGraph(intervals,True)

def RandomLobster(n, p, q, seed=None):
    """
    Returns a random lobster.

    A lobster is a tree that reduces to a caterpillar when pruning all
    leaf vertices. A caterpillar is a tree that reduces to a path when
    pruning all leaf vertices (q=0).

    INPUT:

    -  ``n`` - expected number of vertices in the backbone

    -  ``p`` - probability of adding an edge to the
       backbone

    -  ``q`` - probability of adding an edge (claw) to the
       arms

    -  ``seed`` - for the random number generator


    EXAMPLE: We show the edge list of a random graph with 3 backbone
    nodes and probabilities `p = 0.7` and `q = 0.3`::

        sage: graphs.RandomLobster(3, 0.7, 0.3).edges(labels=False)
        [(0, 1), (1, 2)]

    ::

        sage: G = graphs.RandomLobster(9, .6, .3)
        sage: G.show()  # long time
    """
    if seed is None:
        seed = current_randstate().long_seed()
    import networkx
    return Graph(networkx.random_lobster(n, p, q, seed=seed))

def RandomTree(n):
    """
    Returns a random tree on `n` nodes numbered `0` through `n-1`.

    By Cayley's theorem, there are `n^{n-2}` trees with vertex
    set `\{0,1,...,n-1\}`. This constructor chooses one of these uniformly
    at random.

    ALGORITHM:

    The algoritm works by generating an `(n-2)`-long
    random sequence of numbers chosen independently and uniformly
    from `\{0,1,\ldots,n-1\}` and then applies an inverse
    Prufer transformation.

    INPUT:

    -  ``n`` - number of vertices in the tree

    EXAMPLE::

        sage: G = graphs.RandomTree(10)
        sage: G.is_tree()
        True
        sage: G.show() # long

    TESTS:

    Ensuring that we encounter no unexpected surprise ::

        sage: all( graphs.RandomTree(10).is_tree()
        ....:      for i in range(100) )
        True

    """
    from sage.misc.prandom import randint
    g = Graph()

    # create random Prufer code
    code = [ randint(0,n-1) for i in xrange(n-2) ]

    # We count the number of symbols of each type.
    # count[k] is the no. of times k appears in code
    #
    # (count[k] is set to -1 when the corresponding vertex is not
    # available anymore)
    count = [ 0 for i in xrange(n) ]
    for k in code:
        count[k] += 1

    g.add_vertices(range(n))

    for s in code:
        for x in range(n):
            if count[x] == 0:
                break

        count[x] = -1
        g.add_edge(x,s)
        count[s] -= 1

    # Adding as an edge the last two available vertices
    last_edge = [ v for v in range(n) if count[v] != -1 ]
    g.add_edge(last_edge)

    return g

def RandomTreePowerlaw(n, gamma=3, tries=100, seed=None):
    """
    Returns a tree with a power law degree distribution. Returns False
    on failure.

    From the NetworkX documentation: A trial power law degree sequence
    is chosen and then elements are swapped with new elements from a
    power law distribution until the sequence makes a tree (size = order
    - 1).

    INPUT:

    -  ``n`` - number of vertices

    -  ``gamma`` - exponent of power law

    -  ``tries`` - number of attempts to adjust sequence to
       make a tree

    -  ``seed`` - for the random number generator


    EXAMPLE: We show the edge list of a random graph with 10 nodes and
    a power law exponent of 2.

    ::

        sage: graphs.RandomTreePowerlaw(10, 2).edges(labels=False)
        [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (6, 8), (6, 9)]

    ::

        sage: G = graphs.RandomTreePowerlaw(15, 2)
        sage: if G:
        ....:     G.show()  # random output, long time
    """
    if seed is None:
        seed = current_randstate().long_seed()
    import networkx
    try:
        return Graph(networkx.random_powerlaw_tree(n, gamma, seed=seed, tries=tries))
    except networkx.NetworkXError:
        return False

def RandomRegular(d, n, seed=None):
    """
    Returns a random d-regular graph on n vertices, or returns False on
    failure.

    Since every edge is incident to two vertices, n\*d must be even.

    INPUT:

    -  ``n`` - number of vertices

    -  ``d`` - degree

    -  ``seed`` - for the random number generator


    EXAMPLE: We show the edge list of a random graph with 8 nodes each
    of degree 3.

    ::

        sage: graphs.RandomRegular(3, 8).edges(labels=False)
        [(0, 1), (0, 4), (0, 7), (1, 5), (1, 7), (2, 3), (2, 5), (2, 6), (3, 4), (3, 6), (4, 5), (6, 7)]

    ::

        sage: G = graphs.RandomRegular(3, 20)
        sage: if G:
        ....:     G.show()  # random output, long time

    REFERENCES:

    .. [KimVu2003] Kim, Jeong Han and Vu, Van H. Generating random regular
      graphs. Proc. 35th ACM Symp. on Thy. of Comp. 2003, pp
      213-222. ACM Press, San Diego, CA, USA.
      http://doi.acm.org/10.1145/780542.780576

    .. [StegerWormald1999] Steger, A. and Wormald, N. Generating random
      regular graphs quickly. Prob. and Comp. 8 (1999), pp 377-396.
    """
    if seed is None:
        seed = current_randstate().long_seed()
    import networkx
    try:
        N = networkx.random_regular_graph(d, n, seed=seed)
        if N is False: return False
        return Graph(N, sparse=True)
    except Exception:
        return False

def RandomShell(constructor, seed=None):
    """
    Returns a random shell graph for the constructor given.

    INPUT:

    -  ``constructor`` - a list of 3-tuples (n,m,d), each
       representing a shell

    -  ``n`` - the number of vertices in the shell

    -  ``m`` - the number of edges in the shell

    -  ``d`` - the ratio of inter (next) shell edges to
       intra shell edges

    -  ``seed`` - for the random number generator


    EXAMPLE::

        sage: G = graphs.RandomShell([(10,20,0.8),(20,40,0.8)])
        sage: G.edges(labels=False)
        [(0, 3), (0, 7), (0, 8), (1, 2), (1, 5), (1, 8), (1, 9), (3, 6), (3, 11), (4, 6), (4, 7), (4, 8), (4, 21), (5, 8), (5, 9), (6, 9), (6, 10), (7, 8), (7, 9), (8, 18), (10, 11), (10, 13), (10, 19), (10, 22), (10, 26), (11, 18), (11, 26), (11, 28), (12, 13), (12, 14), (12, 28), (12, 29), (13, 16), (13, 21), (13, 29), (14, 18), (16, 20), (17, 18), (17, 26), (17, 28), (18, 19), (18, 22), (18, 27), (18, 28), (19, 23), (19, 25), (19, 28), (20, 22), (24, 26), (24, 27), (25, 27), (25, 29)]
        sage: G.show()  # long time
    """
    if seed is None:
        seed = current_randstate().long_seed()
    import networkx
    return Graph(networkx.random_shell_graph(constructor, seed=seed))

def RandomToleranceGraph(n):
    r"""
    Returns a random tolerance graph.

    The random tolerance graph is built from a random tolerance representation
    by using the function `ToleranceGraph`. This representation is a list
    `((l_0,r_0,t_0), (l_1,r_1,t_1), ..., (l_k,r_k,t_k))` where `k = n-1` and
    `I_i = (l_i,r_i)` denotes a random interval and `t_i` a random positive
    value. The width of the representation is limited to n**2 * 2**n.

    .. NOTE::

        The vertices are named 0, 1, ..., n-1. The tolerance representation used
        to create the graph is saved with the graph and can be recovered using
        ``get_vertex()`` or ``get_vertices()``.

    INPUT:

    - ``n`` -- number of vertices of the random graph.

    EXAMPLE:

    Every tolerance graph is perfect. Hence, the chromatic number is equal to
    the clique number ::

        sage: g = graphs.RandomToleranceGraph(8)
        sage: g.clique_number() == g.chromatic_number()
        True

    TEST:

        sage: g = graphs.RandomToleranceGraph(-2)
        Traceback (most recent call last):
        ...
        ValueError: The number `n` of vertices must be >= 0.
    """
    from sage.misc.prandom import randint
    from sage.graphs.generators.intersection import ToleranceGraph

    if n<0:
        raise ValueError('The number `n` of vertices must be >= 0.')

    W = n**2 * 2**n

    tolrep = [tuple(sorted((randint(0,W), randint(0,W)))) + (randint(0,W),) for i in range(n)]

    return ToleranceGraph(tolrep)

