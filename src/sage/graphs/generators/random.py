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
from six.moves import range
import sys
# import from Sage library
from sage.graphs.graph import Graph
from sage.misc.randstate import current_randstate
from sage.misc.prandom import randint

def RandomGNP(n, p, seed=None, fast=True, algorithm='Sage'):
    r"""
    Returns a random graph on `n` nodes. Each edge is inserted independently
    with probability `p`.

    INPUT:

    - ``n`` -- number of nodes of the graph

    - ``p`` -- probability of an edge

    - ``seed`` - a ``random.Random`` seed or a Python ``int`` for the random
      number generator (default: ``None``).

    - ``fast`` -- boolean set to True (default) to use the algorithm with
      time complexity in `O(n+m)` proposed in [BatBra2005]_. It is designed
      for generating large sparse graphs. It is faster than other algorithms for
      *LARGE* instances (try it to know whether it is useful for you).

    - ``algorithm`` -- By default (```algorithm='Sage'``), this function uses the
      algorithm implemented in ```sage.graphs.graph_generators_pyx.pyx``. When
      ``algorithm='networkx'``, this function calls the NetworkX function
      ``fast_gnp_random_graph``, unless ``fast=False``, then
      ``gnp_random_graph``. Try them to know which algorithm is the best for
      you. The ``fast`` parameter is not taken into account by the 'Sage'
      algorithm so far.

    REFERENCES:

    .. [ErdRen1959] \P. Erdos and A. Renyi. On Random Graphs, Publ.
       Math. 6, 290 (1959).

    .. [Gilbert1959] \E. N. Gilbert. Random Graphs, Ann. Math. Stat.,
       30, 1141 (1959).

    .. [BatBra2005] \V. Batagelj and U. Brandes. Efficient generation of
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

        sage: graphs.RandomGNP(50,.2,algorithm=50)
        Traceback (most recent call last):
        ...
        ValueError: 'algorithm' must be equal to 'networkx' or to 'Sage'.
        sage: set_random_seed(0)
        sage: graphs.RandomGNP(50,.2, algorithm="Sage").size()
        243
        sage: graphs.RandomGNP(50,.2, algorithm="networkx").size()
        260     # 32-bit 
        245     # 64-bit
    """
    if n < 0:
        raise ValueError("The number of nodes must be positive or null.")
    if 0.0 > p or 1.0 < p:
        raise ValueError("The probability p must be in [0..1].")

    if seed is None:
        seed = int(current_randstate().long_seed() % sys.maxsize)
    if p == 1:
        from sage.graphs.generators.basic import CompleteGraph
        return CompleteGraph(n)

    if algorithm == 'networkx':
        import networkx
        if fast:
            G = networkx.fast_gnp_random_graph(n, p, seed=seed)
        else:
            G = networkx.gnp_random_graph(n, p, seed=seed)
        return Graph(G)
    elif algorithm in ['Sage', 'sage']:
        # We use the Sage generator
        from sage.graphs.graph_generators_pyx import RandomGNP as sageGNP
        return sageGNP(n, p)
    else:
        raise ValueError("'algorithm' must be equal to 'networkx' or to 'Sage'.")

def RandomBarabasiAlbert(n, m, seed=None):
    r"""
    Return a random graph created using the Barabasi-Albert preferential
    attachment model.

    A graph with `m` vertices and no edges is initialized, and a graph of `n`
    vertices is grown by attaching new vertices each with `m` edges that are
    attached to existing vertices, preferentially with high degree.

    INPUT:

    - ``n`` -- number of vertices in the graph

    - ``m`` -- number of edges to attach from each new node

    - ``seed`` -- a ``random.Random`` seed or a Python ``int`` for the random
      number generator (default: ``None``)

    EXAMPLES:

    We show the edge list of a random graph on 6 nodes with `m = 2`::

        sage: G = graphs.RandomBarabasiAlbert(6,2)
        sage: G.order(), G.size()
        (6, 8)
        sage: G.degree_sequence()  # random
        [4, 3, 3, 2, 2, 2]

    We plot a random graph on 12 nodes with `m = 3`::

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

    When `m = 1`, the generated graph is a tree::

        sage: graphs.RandomBarabasiAlbert(6, 1).is_tree()
        True
    """
    if seed is None:
        seed = int(current_randstate().long_seed() % sys.maxsize)
    import networkx
    return Graph(networkx.barabasi_albert_graph(n, m, seed=seed))

def RandomBipartite(n1, n2, p, set_position=False):
    r"""
    Returns a bipartite graph with `n1+n2` vertices such that any edge
    from `[n1]` to `[n2]` exists with probability `p`.

    INPUT:

    - ``n1, n2`` -- Cardinalities of the two sets

    - ``p`` -- Probability for an edge to exist

    - ``set_position`` -- boolean (default ``False``); if set to ``True``, we
      assign positions to the vertices so that the set of cardinality `n1` is
      on the line `y=1` and the set of cardinality `n2` is on the line `y=0`.

    EXAMPLES::

        sage: g = graphs.RandomBipartite(5, 2, 0.5)
        sage: g.vertices()
        [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (1, 0), (1, 1)]

    TESTS::

        sage: g = graphs.RandomBipartite(5, -3, 0.5)
        Traceback (most recent call last):
        ...
        ValueError: n1 and n2 should be integers strictly greater than 0
        sage: g = graphs.RandomBipartite(5, 3, 1.5)
        Traceback (most recent call last):
        ...
        ValueError: parameter p is a probability, and so should be a real value between 0 and 1

    :trac:`12155`::

        sage: graphs.RandomBipartite(5, 6, .2).complement()
        complement(Random bipartite graph of order 5+6 with edge probability 0.200000000000000): Graph on 11 vertices

    Test assigned positions::

        sage: graphs.RandomBipartite(1, 2, .1, set_position=True).get_pos()
        {(0, 0): (1, 1.0), (1, 0): (0, 0), (1, 1): (2.0, 0.0)}
        sage: graphs.RandomBipartite(2, 1, .1, set_position=True).get_pos()
        {(0, 0): (0, 1), (0, 1): (2.0, 1.0), (1, 0): (1, 0.0)}
        sage: graphs.RandomBipartite(2, 2, .1, set_position=True).get_pos()
        {(0, 0): (0, 1), (0, 1): (2.0, 1.0), (1, 0): (0, 0), (1, 1): (2.0, 0.0)}
        sage: graphs.RandomBipartite(2, 2, .1, set_position=False).get_pos()

    """
    if not (p >= 0 and p <= 1):
        raise ValueError("parameter p is a probability, and so should be a real value between 0 and 1")
    if not (n1 > 0 and n2 > 0):
        raise ValueError("n1 and n2 should be integers strictly greater than 0")

    from numpy.random import uniform

    g=Graph(name="Random bipartite graph of order "+str(n1) +"+"+str(n2)+" with edge probability "+str(p))

    S1 = [(0,i) for i in range(n1)]
    S2 = [(1,i) for i in range(n2)]
    g.add_vertices(S1)
    g.add_vertices(S2)

    for w in range(n2):
        for v in range(n1):
            if uniform() <= p :
                g.add_edge((0, v), (1, w))

    # We now assign positions to vertices:
    # - vertices in S1 are placed on the line from (0, 1) to (max(n1, n2), 1)
    # - vertices in S2 are placed on the line from (0, 0) to (max(n1, n2), 0)
    # If S1 or S2 has a single vertex, it is centered in the line.
    if set_position:
        nmax = max(n1, n2)
        g._line_embedding(S1, first=(0, 1), last=(nmax, 1))
        g._line_embedding(S2, first=(0, 0), last=(nmax, 0))

    return g

def RandomRegularBipartite(n1, n2, d1, set_position=False):
    r"""
    Return a random regular bipartite graph on `n1 + n2` vertices.

    The bipartite graph has `n1 * d1` edges. Hence, `n2` must divide `n1 * d1`.
    Each vertex of the set of cardinality `n1` has degree `d1` (which can be at
    most `n2`) and each vertex in the set of cardinality `n2` has degree 
    `(n1 * d1) / n2`. The bipartite graph has no multiple edges.

    This generator implements an algorithm inspired by that of [MW1990]_ for 
    the uniform generation of random regular bipartite graphs. It performs well
    when `d1 = o(n2^{1/3})` or (`n2 - d1 = o(n2^{1/3})`). In other cases, the
    running time can be huge. Note that the currently implemented algorithm
    does not generate uniformly random graphs.

    INPUT:

    - ``n1, n2`` -- number of vertices in each side

    - ``d1`` -- degree of the vertices in the set of cardinality `n1`.

    - ``set_position`` -- boolean (default ``False``); if set to ``True``, we
      assign positions to the vertices so that the set of cardinality `n1` is 
      on the line `y=1` and the set of cardinality `n2` is on the line `y=0`.

    EXAMPLES::

        sage: g = graphs.RandomRegularBipartite(4, 6, 3)
        sage: g.order(), g.size()
        (10, 12)
        sage: set(g.degree())
        {2, 3}

        sage: graphs.RandomRegularBipartite(1, 2, 2, set_position=True).get_pos()
        {0: (1, 1.0), 1: (0, 0), 2: (2.0, 0.0)}
        sage: graphs.RandomRegularBipartite(2, 1, 1, set_position=True).get_pos()
        {0: (0, 1), 1: (2.0, 1.0), 2: (1, 0.0)}
        sage: graphs.RandomRegularBipartite(2, 3, 3, set_position=True).get_pos()
        {0: (0, 1), 1: (3.0, 1.0), 2: (0, 0), 3: (1.5, 0.0), 4: (3.0, 0.0)}
        sage: graphs.RandomRegularBipartite(2, 3, 3, set_position=False).get_pos()

    TESTS:

    Giving invalid parameters::

        sage: graphs.RandomRegularBipartite(0, 2, 1)
        Traceback (most recent call last):
        ...
        ValueError: n1 and n2 must be integers greater than 0
        sage: graphs.RandomRegularBipartite(2, 3, 2)
        Traceback (most recent call last):
        ...
        ValueError: the product n1 * d1 must be a multiple of n2
        sage: graphs.RandomRegularBipartite(1, 1, 2)
        Traceback (most recent call last):
        ...
        ValueError: d1 must be less than or equal to n2
    """
    if n1 < 1 or n2 < 1:
        raise ValueError("n1 and n2 must be integers greater than 0")
    if d1 > n2:
        raise ValueError("d1 must be less than or equal to n2")
    d2 = (n1 * d1) // n2
    if n1 * d1 != n2 * d2:
        raise ValueError("the product n1 * d1 must be a multiple of n2")

    complement = False
    if d1 > n2/2 or d2 > n1/2:
        # We build the complement graph instead
        complement = True
        d1 = n2 - d1
        d2 = n1 - d2

    E = set()
    F = set()

    if d1:
        from sage.misc.prandom import shuffle, choice

        M1 = n1 * d1 * (d1 - 1)
        M2 = n2 * d2 * (d2 - 1)
        M = n1 * d1 + n2 * d2
        UB_parallel = (M1 * M2) / M**2

        # We create a set of n1 * d1 random edges with possible repetitions. We
        # require that the number of repeated edges is bounded and that an edge
        # can be repeated only once.
        L = [u for u in range(n1) for i in range(d1)]
        R = [u for u in range(n1, n1 + n2) for i in range(d2)]
        restart = True
        while restart:
            restart = False
            shuffle(R)
            E = set()
            F = set()
            for e in zip(L, R):
                if e in E:
                    if e in F:
                        # We have more than 2 times e => restart
                        restart = True
                        break
                    else:
                        F.add(e)
                    if len(F) >= UB_parallel:
                        # We have too many parallel edges
                        restart = True
                        break
                else:
                    E.add(e)

    # We remove multiple edges by applying random forward d-switching. That is,
    # given edge e that is repeated twice, we select single edges f and g with
    # no common end points, and then create 4 new edges. We forbid creating new
    # multiple edges.
    while F:
        # random forward d-switching
        e = F.pop()
        E.discard(e)
        TE = tuple(E.difference(F))
        # We select 2 vertex disjoint edges
        while True:
            f = choice(TE)
            if e[0] == f[0] or e[1] == f[1]:
                continue
            g = choice(TE)
            if e[0] != g[0] and e[1] != g[1] and f[0] != g[0] and f[1] != g[1]:
                new_edges = [(f[0], e[1]), (e[0], f[1]), (e[0], g[1]), (g[0], e[1])]
                if not E.intersection(new_edges):
                    # We are not creating new parallel edges.
                    # To generate uniformly random graphs we would have to
                    # implement a probabilistic restart of the whole algorithm
                    # here, see [MW1990].
                    break
        E.discard(f)
        E.discard(g)
        E.update(new_edges)

    if complement:
        from sage.graphs.generators.basic import CompleteBipartiteGraph
        E = E.symmetric_difference(CompleteBipartiteGraph(n1, n2).edges(labels=False))
        d1, d2 = n2 - d1, n1 - d2

    name = "Random regular bipartite graph of order {}+{} and degrees {} and {}".format(n1, n2, d1, d2)
    G = Graph(list(E), name=name)

    # We now assign positions to vertices:
    # - vertices 0,..,n1-1 are placed on the line (0, 1) to (max(n1, n2), 1)
    # - vertices n1,..,n1+n2-1 are placed on the line (0, 0) to (max(n1, n2), 0)
    # If n1 (or n2) is 1, the vertex is centered in the line.
    if set_position:
        nmax = max(n1, n2)
        G._line_embedding(list(range(n1)), first=(0, 1), last=(nmax, 1))
        G._line_embedding(list(range(n1, n1+n2)), first=(0, 0), last=(nmax, 0))

    return G


def RandomBlockGraph(m, k, kmax=None, incidence_structure=False):
    r"""
    Return a Random Block Graph.

    A block graph is a connected graph in which every biconnected component
    (block) is a clique.

    .. SEEALSO::

        - :wikipedia:`Block_graph` for more details on these graphs
        - :meth:`~sage.graphs.graph.Graph.is_block_graph` -- test if a graph is a block graph
        - :meth:`~sage.graphs.generic_graph.GenericGraph.blocks_and_cut_vertices`
        - :meth:`~sage.graphs.generic_graph.GenericGraph.blocks_and_cuts_tree`
        - :meth:`~sage.combinat.designs.incidence_structures.IncidenceStructure` 

    INPUT:

    - ``m`` -- integer; number of blocks (at least one).

    - ``k`` -- integer; minimum number of vertices of a block (at least two).

    - ``kmax`` -- integer (default: ``None``) By default, each block has `k`
      vertices. When the parameter `kmax` is specified (with `kmax \geq k`), the
      number of vertices of each block is randomly chosen between `k` and
      `kmax`.

    - ``incidence_structure`` -- boolean (default: ``False``) when set to
      ``True``, the incidence structure of the graphs is returned instead of the
      graph itself, that is the list of the lists of vertices in each
      block. This is useful for the creation of some hypergraphs.

    OUTPUT:

    A Graph when ``incidence_structure==False`` (default), and otherwise an
    incidence structure.

    EXAMPLES:

    A block graph with a single block is a clique::

        sage: B = graphs.RandomBlockGraph(1, 4)
        sage: B.is_clique()
        True

    A block graph with blocks of order 2 is a tree::

        sage: B = graphs.RandomBlockGraph(10, 2)
        sage: B.is_tree()
        True

    Every biconnected component of a block graph is a clique::

        sage: B = graphs.RandomBlockGraph(5, 3, kmax=6)
        sage: blocks,cuts = B.blocks_and_cut_vertices()
        sage: all(B.is_clique(block) for block in blocks)
        True

    A block graph with blocks of order `k` has `m*(k-1)+1` vertices::

        sage: m, k = 6, 4
        sage: B = graphs.RandomBlockGraph(m, k)
        sage: B.order() == m*(k-1)+1
        True

    Test recognition methods::

        sage: B = graphs.RandomBlockGraph(6, 2, kmax=6)
        sage: B.is_block_graph()
        True
        sage: B in graph_classes.Block
        True

    Asking for the incidence structure::

        sage: m, k = 6, 4
        sage: IS = graphs.RandomBlockGraph(m, k, incidence_structure=True)
        sage: from sage.combinat.designs.incidence_structures import IncidenceStructure
        sage: IncidenceStructure(IS)
        Incidence structure with 19 points and 6 blocks
        sage: m*(k-1)+1
        19

    TESTS:

    A block graph has at least one block, so `m\geq 1`::

        sage: B = graphs.RandomBlockGraph(0, 1)
        Traceback (most recent call last):
        ...
        ValueError: the number `m` of blocks must be >= 1

    A block has at least 2 vertices, so `k\geq 2`::

        sage: B = graphs.RandomBlockGraph(1, 1)
        Traceback (most recent call last):
        ...
        ValueError: the minimum number `k` of vertices in a block must be >= 2

    The maximum size of a block is at least its minimum size, so `k\leq kmax`::

        sage: B = graphs.RandomBlockGraph(1, 3, kmax=2)
        Traceback (most recent call last):
        ...
        ValueError: the maximum number `kmax` of vertices in a block must be >= `k`
    """
    from sage.misc.prandom import choice
    from sage.sets.disjoint_set import DisjointSet

    if m < 1:
        raise ValueError("the number `m` of blocks must be >= 1")
    if k < 2:
        raise ValueError("the minimum number `k` of vertices in a block must be >= 2")
    if kmax is None:
        kmax = k
    elif kmax < k:
        raise ValueError("the maximum number `kmax` of vertices in a block must be >= `k`")

    if m == 1:
        # A block graph with a single block is a clique
        IS = [ list(range(randint(k, kmax))) ]
        
    elif kmax == 2:
        # A block graph with blocks of order 2 is a tree
        IS = [ list(e) for e in RandomTree(m+1).edges(labels=False) ]

    else:
        # We start with a random tree of order m
        T = RandomTree(m)

        # We create a block of order in range [k,kmax] per vertex of the tree
        B = {u:[(u,i) for i in range(randint(k, kmax))] for u in T}

        # For each edge of the tree, we choose 1 vertex in each of the
        # corresponding blocks and we merge them. We use a disjoint set data
        # structure to keep a unique identifier per merged vertices
        DS = DisjointSet([i for u in B for i in B[u]])
        for u,v in T.edges(labels=0):
            DS.union(choice(B[u]), choice(B[v]))

        # We relabel vertices in the range [0, m*(k-1)] and build the incidence
        # structure
        new_label = {root:i for i,root in enumerate(DS.root_to_elements_dict())}
        IS = [ [new_label[DS.find(v)] for v in B[u]] for u in B ]

    if incidence_structure:
        return IS
    
    # We finally build the block graph
    if k == kmax:
        BG = Graph(name = "Random Block Graph with {} blocks of order {}".format(m, k))
    else:
        BG = Graph(name = "Random Block Graph with {} blocks of order {} to {}".format(m, k, kmax))
    for block in IS:
        BG.add_clique( block )
    return BG


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

    EXAMPLES:

    Every (bounded) tolerance graph is perfect. Hence, the
    chromatic number is equal to the clique number ::

        sage: g = graphs.RandomBoundedToleranceGraph(8)
        sage: g.clique_number() == g.chromatic_number()
        True
    """
    from sage.misc.prandom import randint
    from sage.graphs.generators.intersection import ToleranceGraph

    W = n ** 2 * 2 ** n

    tolrep = [(l_r[0], l_r[1], randint(0, l_r[1] - l_r[0])) for l_r in [sorted((randint(0, W), randint(0, W))) for i in range(n)]]

    return ToleranceGraph(tolrep)

def RandomGNM(n, m, dense=False, seed=None):
    """
    Returns a graph randomly picked out of all graphs on n vertices
    with m edges.

    INPUT:

    - ``n`` - number of vertices.

    - ``m`` - number of edges.

    - ``dense`` - whether to use NetworkX's
      dense_gnm_random_graph or gnm_random_graph

    - ``seed`` - a ``random.Random`` seed or a Python ``int`` for the random
      number generator (default: ``None``).


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
        seed = int(current_randstate().long_seed() % sys.maxsize)
    import networkx
    if dense:
        return Graph(networkx.dense_gnm_random_graph(n, m, seed=seed))
    else:
        return Graph(networkx.gnm_random_graph(n, m, seed=seed))

def RandomNewmanWattsStrogatz(n, k, p, seed=None):
    r"""
    Return a Newman-Watts-Strogatz small world random graph on `n` vertices.

    From the NetworkX documentation: first create a ring over `n` nodes.  Then
    each node in the ring is connected with its `k` nearest neighbors. Then
    shortcuts are created by adding new edges as follows: for each edge `u-v` in
    the underlying "`n`-ring with `k` nearest neighbors"; with probability `p`
    add a new edge `u-w` with randomly-chosen existing node `w`. In contrast
    with ``networkx.watts_strogatz_graph()``, no edges are removed.

    INPUT:

    - ``n`` -- number of vertices

    - ``k`` -- each vertex is connected to its `k` nearest neighbors

    - ``p`` -- the probability of adding a new edge for each edge

    - ``seed`` -- a ``random.Random`` seed or a Python ``int`` for the random
      number generator (default: ``None``)

    EXAMPLES:

    We check that the generated graph contains a cycle of order `n`::

        sage: G = graphs.RandomNewmanWattsStrogatz(7, 2, 0.2)
        sage: G.order(), G.size()
        (7, 9)
        sage: C7 = graphs.CycleGraph(7)
        sage: G.subgraph_search(C7)
        Subgraph of (): Graph on 7 vertices
        sage: G.diameter() <= C7.diameter()
        True

    ::

        sage: G = graphs.RandomNewmanWattsStrogatz(12, 2, .3)
        sage: G.show()  # long time

    TESTS:

    We check that when `k = 2` and `p = 0`, the generated graph is a cycle::

        sage: G = graphs.RandomNewmanWattsStrogatz(7, 2, 0)
        sage: G.is_cycle()
        True

    We check that when `k = 4` and `p = 0`, the generated graph is a circulant
    graph of parameters ``[1, 2]``::

        sage: G = graphs.RandomNewmanWattsStrogatz(7, 4, 0)
        sage: G.is_isomorphic(graphs.CirculantGraph(7, [1, 2]))
        True

    REFERENCE:

    .. [NWS99] Newman, M.E.J., Watts, D.J. and Strogatz, S.H.  Random
      graph models of social networks. Proc. Nat. Acad. Sci. USA
      99, 2566-2572.
    """
    if seed is None:
        seed = int(current_randstate().long_seed() % sys.maxsize)
    import networkx
    return Graph(networkx.newman_watts_strogatz_graph(n, k, p, seed=seed))

def RandomHolmeKim(n, m, p, seed=None):
    r"""
    Return a random graph generated by the Holme and Kim algorithm for
    graphs with power law degree distribution and approximate average
    clustering.

    INPUT:

    - ``n`` -- number of vertices

    - ``m`` -- number of random edges to add for each new node

    - ``p`` -- probability of adding a triangle after adding a random edge

    - ``seed`` -- a ``random.Random`` seed or a Python ``int`` for the random
      number generator (default: ``None``)

    From the NetworkX documentation: the average clustering has a hard time
    getting above a certain cutoff that depends on `m`. This cutoff is often
    quite low. Note that the transitivity (fraction of triangles to possible
    triangles) seems to go down with network size. It is essentially the
    Barabasi-Albert growth model with an extra step that each random edge is
    followed by a chance of making an edge to one of its neighbors too (and thus
    a triangle). This algorithm improves on B-A in the sense that it enables a
    higher average clustering to be attained if desired. It seems possible to
    have a disconnected graph with this algorithm since the initial `m` nodes
    may not be all linked to a new node on the first iteration like the BA
    model.

    EXAMPLES:

    We check that a random graph on 8 nodes with 2 random edges per node and a
    probability `p = 0.5` of forming triangles contains a triangle::

        sage: G = graphs.RandomHolmeKim(8, 2, 0.5)
        sage: G.order(), G.size()
        (8, 12)
        sage: C3 = graphs.CycleGraph(3)
        sage: G.subgraph_search(C3)
        Subgraph of (): Graph on 3 vertices

    ::

        sage: G = graphs.RandomHolmeKim(12, 3, .3)
        sage: G.show()  # long time

    REFERENCE:

    .. [HolmeKim2002] Holme, P. and Kim, B.J. Growing scale-free networks
      with tunable clustering, Phys. Rev. E (2002). vol 65, no 2, 026107.
    """
    if seed is None:
        seed = int(current_randstate().long_seed() % sys.maxsize)
    import networkx
    return Graph(networkx.powerlaw_cluster_graph(n, m, p, seed=seed))


def RandomIntervalGraph(n):
    r"""
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

    EXAMPLES:

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

# Random Chordal Graphs

def growing_subtrees(T, k):
    r"""
    Return a list of the vertex sets of ``n`` randomly chosen subtrees of ``T``.

    For a tree of order `n`, the collection contains `n` subtrees with maximum
    order `k` and average order `\frac{k + 1}{2}`.

    This method is part of
    :meth:`~sage.graphs.generators.random.RandomChordalGraph`.

    ALGORITHM:

    For each subtree `T_i`, the algorithm picks a size `k_i` randomly from
    `[1,k]`. Then a random node of `T` is chosen as the first node of `T_i`. In
    each of the subsequent `k_i - 1` iterations, it picks a random node in the
    neighborhood of `T_i` and adds it to `T_i`.

    See [SHET2018]_ for more details.

    INPUT:

    - ``T`` -- a tree

    - ``k`` -- a strictly positive integer; maximum size of a subtree

    EXAMPLES::

        sage: from sage.graphs.generators.random import growing_subtrees
        sage: T = graphs.RandomTree(10)
        sage: S = growing_subtrees(T, 5)
        sage: len(S)
        10
    """
    from sage.misc.prandom import sample
    n = T.order()
    S = []
    for _ in range(n):
        ki = randint(1, k)
        if ki == n:
            Vi = frozenset(T)
        else:
            x = T.random_vertex()
            Ti = set([x])
            neighbors = set(T.neighbor_iterator(x))
            for j in range(ki - 1):
                # Select a random neighbor z outside of Ti and add it to Ti
                z = sample(neighbors, 1)[0]
                Ti.add(z)
                neighbors.update(y for y in T.neighbor_iterator(z) if y not in Ti)
            Vi = frozenset(Ti)
        S.append(Vi)

    return S

def connecting_nodes(T, l):
    r"""
    Return a list of the vertex sets of ``n`` randomly chosen subtrees of ``T``.

    This method is part of
    :meth:`~sage.graphs.generators.random.RandomChordalGraph`.

    ALGORITHM:

    For each subtree `T_i`, we first select `k_i` nodes of `T`, where `k_i` is a
    random integer from a Poisson distribution with mean `l`. `T_i` is then
    generated to be the minimal subtree that contains the selected `k_i`
    nodes. This implies that a subtree will most likely have many more nodes
    than those selected initially, and this must be taken into consideration
    when choosing `l`.

    See [SHET2018]_ for more details.

    INPUT:

    - ``T`` -- a tree

    - ``l`` -- a strictly positive real number; mean of a Poisson distribution

    EXAMPLES::

        sage: from sage.graphs.generators.random import connecting_nodes
        sage: T = graphs.RandomTree(10)
        sage: S = connecting_nodes(T, 5)
        sage: len(S)
        10
    """
    from sage.combinat.permutation import Permutations
    from sage.data_structures.bitset import Bitset
    from numpy.random import poisson

    n = T.order()
    V = list(T)
    P = Permutations(V)
    active = Bitset(capacity=n)

    # Choose a root
    root = T.random_vertex()

    # Perform BFS from root and identify parent in root to leaf orientation
    parent = {root: root}
    dist = {root: 0}
    bfs = [root]
    i = 0
    while i < n:
        u = bfs[i]
        d = dist[u]
        for v in T.neighbor_iterator(u):
            if v not in parent:
                parent[v] = u
                dist[v] = d + 1
                bfs.append(v)
        i += 1

    S = []
    for _ in range(n):
        ki = poisson(l)
        if not ki:
            ki = 1
        elif ki >= n:
            Ti = frozenset(V)

        if ki < n:
            # Select ki vertices at random
            Vi = set(P.random_element()[:ki])
            # Arrange them by distance to root and mark them as active
            d = max(dist[u] for u in Vi)
            Li = [set() for _ in range(d + 1)]
            active.clear()
            for u in Vi:
                Li[dist[u]].add(u)
                active.add(u)
            # Add to Vi the vertices of a minimal subtree containing Vi.
            # To do so, add the parents of the vertices at distance d to Vi,
            # mark them as active and add them to the set of vertices at
            # distance d - 1. Then mark the vertices at distance d as
            # inactive. Repeat the same procedure for the vertices at distance
            # d - 1, d - 2, etc. This procedure ends when at most one active
            # vertex remains.
            while len(active) > 1:
                for u in Li[d]:
                    p = parent[u]
                    Vi.add(p)
                    Li[d - 1].add(p)
                    active.add(p)
                    active.discard(u)
                d -= 1
            Ti = frozenset(Vi)

        S.append(Ti)

    return S

def pruned_tree(T, f, s):
    r"""
    Return a list of the vertex sets of ``n`` randomly chosen subtrees of ``T``.

    This method is part of
    :meth:`~sage.graphs.generators.random.RandomChordalGraph`.

    ALGORITHM:

    For each subtree `T_i`, it randomly selects a fraction `f` of the edges on
    the tree and removes them. The number of edges to delete, say `l`, is
    calculated as `\lfloor((n - 1)f\rfloor`, which will leave `l + 1` subtrees
    in total. Then, it determines the sizes of the `l + 1` subtrees and stores
    the distinct values. Finally, it picks a random size `k_i` from the set of
    largest `100(1-s)\%` of distinct values, and randomly chooses a subtree with
    size `k_i`.

    See [SHET2018]_ for more details.

    INPUT:

    - ``T`` -- a tree

    - ``f`` -- a rational number; the edge deletion fraction. This value must be
      choosen in `[0..1]`.

    - ``s`` -- a real number between 0 and 1; selection barrier for the size of
      trees

    EXAMPLES::

        sage: from sage.graphs.generators.random import pruned_tree
        sage: T = graphs.RandomTree(11)
        sage: S = pruned_tree(T, 1/10, 0.5)
        sage: len(S)
        11
    """
    n = T.order()
    ke = int((n - 1) * f)
    if not ke:
        # No removed edge. Only one possible subtree
        return [tuple(T)] * n
    elif ke == n - 1:
        # All edges are removed. Only n possible subtrees
        return [(u,) for u in T]

    random_edge_iterator = T.random_edge_iterator(labels=False)
    TT = T.copy()
    S = []

    for _ in range(n):
        # Choose ke = (n - 1) * f edges and remove them from TT
        E = set()
        while len(E) < ke:
            E.add(next(random_edge_iterator))
        TT.delete_edges(E)

        # Compute the connected components of TT and arrange them by sizes
        CC = {}
        for c in TT.connected_components(sort=False):
            l = len(c)
            if l in CC:
                CC[l].append(c)
            else:
                CC[l] = [c]

        # Randomly select a subtree size ki from the highest 100(1 - s) %
        # subtree sizes
        sizes = sorted(set(CC.keys()), reverse=True)
        ki = sizes[randint(0, int(len(sizes) * (1 - s)))]

        # Randomly select a subtree of size ki
        Ti = frozenset(CC[ki][randint(0, len(CC[ki]) - 1)])

        S.append(Ti)

        TT.add_edges(E)

    return S

def RandomChordalGraph(n, algorithm="growing", k=None, l=None, f=None, s=None):
    r"""
    Return a random chordal graph of order ``n``.

    A Graph `G` is said to be chordal if it contains no induced hole (a cycle of
    length at least 4). Equivalently, `G` is chordal if it has a perfect
    elimination orderings, if each minimal separator is a clique, or if it is
    the intersection graphs of subtrees of a tree. See the
    :wikipedia:`Chordal_graph`.

    This generator implements the algorithms proposed in [SHET2018]_ for
    generating random chordal graphs as the intersection graph of `n` subtrees
    of a tree of order `n`.

    The returned graph is not necessarily connected.

    INPUT:

    - ``n`` -- integer; the number of nodes of the graph

    - ``algorithm`` -- string (default: ``"growing"``); the choice of the
      algorithm for randomly selecting `n` subtrees of a random tree of order
      `n`. Possible choices are:

      - ``"growing"`` -- for each subtree `T_i`, the algorithm picks a size
        `k_i` randomly from `[1,k]`. Then a random node of `T` is chosen as the
        first node of `T_i`. In each of the subsequent `k_i - 1` iterations, it
        picks a random node in the neighborhood of `T_i` and adds it to `T_i`.

      - ``"connecting"`` -- for each subtree `T_i`, it first selects `k_i` nodes
        of `T`, where `k_i` is a random integer from a Poisson distribution with
        mean `l`. `T_i` is then generated to be the minimal subtree containing
        the selected `k_i` nodes. This implies that a subtree will most likely
        have many more nodes than those selected initially, and this must be
        taken into consideration when choosing `l`.

      - ``"pruned"`` -- for each subtree `T_i`, it randomly selects a fraction
        `f` of the edges on the tree and removes them. The number of edges to
        delete, say `l`, is calculated as `\lfloor (n - 1) f \rfloor`, which will
        leave `l + 1` subtrees in total. Then, it determines the sizes of the `l
        + 1` subtrees and stores the distinct values. Finally, it picks a random
        size `k_i` from the set of largest `100(1-s)\%` of distinct values, and
        randomly chooses a subtree with size `k_i`.

    - ``k`` -- integer (default: ``None``); maximum size of a subtree. If not
      specified (``None``), the maximum size is set to `\sqrt{n}`.
      This parameter is used only when ``algorithm="growing"``. See
      :meth:`~sage.graphs.generators.random.growing_subtrees` for more details.

    - ``l`` -- a strictly positive real number (default: ``None``); mean of a
      Poisson distribution. If not specified, the mean in set to `\log_2{n}`.
      This parameter is used only when ``algorithm="connecting"``. See
      :meth:`~sage.graphs.generators.random.connecting_nodes` for more details.

    - ``f`` -- a rational number (default: ``None``); the edge deletion
      fraction. This value must be choosen in `[0..1]`. If not specified, this
      parameter is set to `\frac{1}{n-1}`.
      This parameter is used only when ``algorithm="pruned"``.
      See :meth:`~sage.graphs.generators.random.pruned_tree` for more details.

    - ``s`` -- a real number between 0 and 1 (default: ``None``); selection
      barrier for the size of trees. If not specified, this parameter is set to
      `0.5`. This parameter is used only when ``algorithm="pruned"``.
      See :meth:`~sage.graphs.generators.random.pruned_tree` for more details.

    EXAMPLES::

        sage: from sage.graphs.generators.random import RandomChordalGraph
        sage: T = RandomChordalGraph(20, algorithm="growing", k=5)
        sage: T.is_chordal()
        True
        sage: T = RandomChordalGraph(20, algorithm="connecting", l=3)
        sage: T.is_chordal()
        True
        sage: T = RandomChordalGraph(20, algorithm="pruned", f=1/3, s=.5)
        sage: T.is_chordal()
        True

    TESTS::

        sage: from sage.graphs.generators.random import RandomChordalGraph
        sage: all(RandomChordalGraph(i).is_chordal() for i in range(4))
        True
        sage: RandomChordalGraph(3, algorithm="Carmen Cru")
        Traceback (most recent call last):
        ...
        NotImplementedError: unknown algorithm 'Carmen Cru'
        sage: RandomChordalGraph(3, algorithm="growing", k=0)
        Traceback (most recent call last):
        ...
        ValueError: parameter k must be >= 1
        sage: RandomChordalGraph(3, algorithm="connecting", l=0)
        Traceback (most recent call last):
        ...
        ValueError: parameter l must be > 0
        sage: RandomChordalGraph(3, algorithm="pruned", f=2)
        Traceback (most recent call last):
        ...
        ValueError: parameter f must be 0 <= f <= 1
        sage: RandomChordalGraph(3, algorithm="pruned", s=1)
        Traceback (most recent call last):
        ...
        ValueError: parameter s must be 0 < s < 1

    .. SEEALSO::

        - :meth:`~sage.graphs.graph_generators.growing_subtrees`
        - :meth:`~sage.graphs.graph_generators.connecting_nodes`
        - :meth:`~sage.graphs.graph_generators.pruned_tree`
        - :wikipedia:`Chordal_graph`
        - :meth:`~sage.graphs.generic_graph.GenericGraph.is_chordal`
        - :meth:`~sage.graphs.graph_generators.GraphGenerators.IntersectionGraph`
    """
    if n < 2:
        return Graph(n, name="Random Chordal Graph")

    # 1. Generate a random tree of order n
    T = RandomTree(n)

    # 2. Generate n non-empty subtrees of T: {T1,...,Tn}
    if algorithm == "growing":
        if k is None:
            from sage.rings.integer import Integer
            k = int(Integer(n).sqrt())
        elif k < 1:
            raise ValueError("parameter k must be >= 1")

        S = growing_subtrees(T, k)

    elif algorithm == "connecting":
        if l is None:
            from sage.rings.integer import Integer
            l = Integer(n).log(2)
        elif l <= 0:
            raise ValueError("parameter l must be > 0")

        S = connecting_nodes(T, l)

    elif algorithm == "pruned":
        if f is None:
            from sage.rings.rational import Rational
            f = 1 / Rational(n - 1)
        elif f < 0 or f > 1:
            raise ValueError("parameter f must be 0 <= f <= 1")
        if s is None:
            s = .5
        elif s <= 0 or s >= 1:
            raise ValueError("parameter s must be 0 < s < 1")

        S = pruned_tree(T, f, s)

    else:
        raise NotImplementedError("unknown algorithm '{}'".format(algorithm))

    # 3. Build the intersection graph of {V(T1),...,V(Tn)}
    vertex_to_subtrees = [[] for _ in range(n)]
    for i,s in enumerate(S):
        for x in s:
            vertex_to_subtrees[x].append(i)
    G = Graph(n, name="Random Chordal Graph")
    for X in vertex_to_subtrees:
        G.add_clique(X)

    return G


def RandomLobster(n, p, q, seed=None):
    """
    Returns a random lobster.

    A lobster is a tree that reduces to a caterpillar when pruning all
    leaf vertices. A caterpillar is a tree that reduces to a path when
    pruning all leaf vertices (q=0).

    INPUT:

    - ``n`` - expected number of vertices in the backbone

    - ``p`` - probability of adding an edge to the
      backbone

    - ``q`` - probability of adding an edge (claw) to the
      arms

    - ``seed`` - a ``random.Random`` seed or a Python ``int`` for the random
      number generator (default: ``None``).


    EXAMPLES: We show the edge list of a random graph with 3 backbone
    nodes and probabilities `p = 0.7` and `q = 0.3`::

        sage: graphs.RandomLobster(3, 0.7, 0.3).edges(labels=False)
        []                                                                  # 32-bit
        [(0, 1), (0, 5), (1, 2), (1, 6), (2, 3), (2, 7), (3, 4), (3, 8)]    # 64-bit

    ::

        sage: G = graphs.RandomLobster(9, .6, .3)
        sage: G.show()  # long time
    """
    if seed is None:
        seed = int(current_randstate().long_seed() % sys.maxsize)
    import networkx
    return Graph(networkx.random_lobster(n, p, q, seed=seed))


def RandomTree(n):
    r"""
    Returns a random tree on `n` nodes numbered `0` through `n-1`.

    By Cayley's theorem, there are `n^{n-2}` trees with vertex
    set `\{0,1,...,n-1\}`. This constructor chooses one of these uniformly
    at random.

    ALGORITHM:

    The algorithm works by generating an `(n-2)`-long
    random sequence of numbers chosen independently and uniformly
    from `\{0,1,\ldots,n-1\}` and then applies an inverse
    Prufer transformation.

    INPUT:

    -  ``n`` - number of vertices in the tree

    EXAMPLES::

        sage: G = graphs.RandomTree(10)
        sage: G.is_tree()
        True
        sage: G.show() # long time

    TESTS:

    Ensuring that we encounter no unexpected surprise ::

        sage: all( graphs.RandomTree(10).is_tree()
        ....:      for i in range(100) )
        True

    """
    from sage.misc.prandom import randint
    g = Graph()

    # create random Prufer code
    code = [ randint(0,n-1) for i in range(n-2) ]

    # We count the number of symbols of each type.
    # count[k] is the no. of times k appears in code
    #
    # (count[k] is set to -1 when the corresponding vertex is not
    # available anymore)
    count = [0] * n
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


def RandomTreePowerlaw(n, gamma=3, tries=1000, seed=None):
    """
    Return a tree with a power law degree distribution, or ``False`` on failure.

    From the NetworkX documentation: a trial power law degree sequence is chosen
    and then elements are swapped with new elements from a power law
    distribution until the sequence makes a tree (size = order - 1).

    INPUT:

    - ``n`` -- number of vertices

    - ``gamma`` -- exponent of power law distribution

    - ``tries`` -- number of attempts to adjust sequence to make a tree

    - ``seed`` -- a ``random.Random`` seed or a Python ``int`` for the random
      number generator (default: ``None``)


    EXAMPLES:

    We check that the generated graph is a tree::

        sage: G = graphs.RandomTreePowerlaw(10, 3)
        sage: G.is_tree()
        True
        sage: G.order(), G.size()
        (10, 9)

    ::

        sage: G = graphs.RandomTreePowerlaw(15, 2)
        sage: if G:
        ....:     G.show()  # random output, long time
    """
    if seed is None:
        seed = int(current_randstate().long_seed() % sys.maxsize)
    import networkx
    try:
        return Graph(networkx.random_powerlaw_tree(n, gamma, seed=seed, tries=tries))
    except networkx.NetworkXError:
        return False


def RandomRegular(d, n, seed=None):
    r"""
    Return a random `d`-regular graph on `n` vertices, or ``False`` on failure.

    Since every edge is incident to two vertices, `n\times d` must be even.

    INPUT:

    - ``d`` -- degree

    - ``n`` -- number of vertices

    - ``seed`` -- a ``random.Random`` seed or a Python ``int`` for the random
      number generator (default: ``None``)

    EXAMPLES:

    We check that a random graph with 8 nodes each of degree 3 is 3-regular::

        sage: G = graphs.RandomRegular(3, 8)
        sage: G.is_regular(k=3)
        True
        sage: G.degree_histogram()
        [0, 0, 0, 8]

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
        seed = int(current_randstate().long_seed() % sys.maxsize)
    import networkx
    try:
        N = networkx.random_regular_graph(d, n, seed=seed)
        if N is False:
            return False
        return Graph(N, sparse=True)
    except Exception:
        return False

def RandomShell(constructor, seed=None):
    """
    Return a random shell graph for the constructor given.

    INPUT:

    - ``constructor`` -- a list of 3-tuples `(n, m, d)`, each representing a
      shell, where:

      - ``n`` -- the number of vertices in the shell

      - ``m`` -- the number of edges in the shell

      - ``d`` -- the ratio of inter (next) shell edges to intra shell edges

    - ``seed`` -- a ``random.Random`` seed or a Python ``int`` for the random
      number generator (default: ``None``)

    EXAMPLES::

        sage: G = graphs.RandomShell([(10,20,0.8),(20,40,0.8)])
        sage: G.order(), G.size()
        (30, 52)
        sage: G.show()  # long time
    """
    if seed is None:
        seed = int(current_randstate().long_seed() % sys.maxsize)
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

    EXAMPLES:

    Every tolerance graph is perfect. Hence, the chromatic number is equal to
    the clique number ::

        sage: g = graphs.RandomToleranceGraph(8)
        sage: g.clique_number() == g.chromatic_number()
        True

    TESTS::

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


# uniform random triangulation using Schaeffer-Poulalhon algorithm


def _auxiliary_random_word(n):
    r"""
    Return a random word used to generate random triangulations.

    INPUT:

    n -- an integer

    OUTPUT:

    A binary sequence `w` of length `4n-2` with `n-1` ones, such that any proper
    prefix `u` of `w` satisfies `3|u|_1 - |u|_0 > -2` (where `|u|_1` and `|u|_0`
    are respectively the number of 1s and 0s in `u`). Those words are the
    expected input of :func:`_contour_and_graph_from_word`.

    ALGORITHM:

    A random word with these numbers of `0` and `1` is chosen. This
    word is then rotated in order to give an admissible code for a
    tree (as explained in Proposition 4.2, [PS2006]_). There are
    exactly two such rotations, one of which is chosen at random.

    Let us consider a word `w` satisfying the expected conditions. By
    drawing a step (1,3) for each 1 and a step (1,-1) for each 0 in
    `w`, one gets a path starting at height 0, ending at height -2 and
    staying above (or on) the horizontal line of height -1 except at the
    end point. By cutting the word at the first position of height -1,
    let us write `w=uv`. One can then see that `v` can only touch the line
    of height -1 at its initial point and just before its end point
    (these two points may be the same).

    Now consider a word `w'` obtained from `w` by any
    rotation. Because `vu` is another word satisfying the expected
    conditions, one can assume that `w'` is obtained from `w` by
    starting at some point in `u`. The algorithm must then recognize
    the end of `u` and the end of `v` inside `w'`. The end of `v` is
    the unique point of minimal height `h`. The end of `u` is the first
    point reaching the height `h+1`.

    EXAMPLES::

        sage: from sage.graphs.generators.random import _auxiliary_random_word
        sage: _auxiliary_random_word(4)  # random
        [1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0]

        sage: def check(w):
        ....:     steps = {1: 3, 0: -1}
        ....:     return all(sum(steps[u] for u in w[:i]) >= -1 for i in range(len(w)))

        sage: for n in range(1, 10):
        ....:     w = _auxiliary_random_word(n)
        ....:     assert len(w) == 4 * n - 2
        ....:     assert w.count(0) == 3 * n - 1
        ....:     assert check(w)
    """
    from sage.misc.prandom import shuffle
    w = [0] * (3 * n - 1) + [1] * (n - 1)
    shuffle(w)

    # Finding the two admissible shifts.
    # The 'if height' is true at least once.
    # If it is true just once, then the word is admissible
    # and cuts = [0, first position of -1] (ok)
    # Otherwise, cuts will always contain
    # [first position of hmin, first position of hmin - 1] (ok)
    cuts = [0, 0]
    height = 0
    height_min = 0
    for i in range(4 * n - 3):
        if w[i] == 1:
            height += 3
        else:
            height -= 1
            if height < height_min:
                height_min = height
                cuts = cuts[1], i + 1

    # random choice of one of the two possible cuts
    idx = cuts[randint(0, 1)]
    return w[idx:] + w[:idx]


def _contour_and_graph_from_word(w):
    r"""
    Return the contour word and the graph of inner vertices of the tree
    associated with the word `w`.

    INPUT:

    - `w` -- a word in `0` and `1` as given by :func:`_auxiliary_random_word`

    This word must satisfy the conditions described in Proposition 4.2 of
    [PS2006]_ (see :func:`_auxiliary_random_word`).

    OUTPUT:

    a pair ``(seq, G)`` where:

    - ``seq`` is a sequence of pairs (label, integer) representing the
      contour walk along the tree associated with `w`

    - ``G`` is the tree obtained by restriction to the set of inner vertices

    The underlying bijection from words to trees is given by lemma 4.1
    in [PS2006]_. It maps the admissible words to planar trees where
    every inner vertex has two leaves.

    In the word `w`, the letter `1` means going away from the root ("up") from
    an inner vertex to another inner vertex. The letter `0` denotes all other
    steps of the discovery, i.e. either discovering a leaf vertex or going
    toward the root ("down"). Thus, the length of `w` is twice the number of
    edges between inner vertices, plus the number of leaves.

    Inner vertices are tagged with 'in' and leaves are tagged with
    'lf'. Inner vertices are moreover labelled by integers, and leaves
    by the label of the neighbor inner vertex.

    EXAMPLES::

        sage: from sage.graphs.generators.random import _contour_and_graph_from_word
        sage: seq, G = _contour_and_graph_from_word([1,0,0,0,0,0])
        sage: seq
        [('in', 0),
         ('in', 1),
         ('lf', 1),
         ('in', 1),
         ('lf', 1),
         ('in', 1),
         ('in', 0),
         ('lf', 0),
         ('in', 0),
         ('lf', 0)]
        sage: G
        Graph on 2 vertices

        sage: from sage.graphs.generators.random import _auxiliary_random_word
        sage: seq, G = _contour_and_graph_from_word(_auxiliary_random_word(20))
        sage: G.is_tree()
        True

        sage: longw = [1,1,0,1,0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0]
        sage: seq, G = _contour_and_graph_from_word(longw)
        sage: G.get_embedding()
        {0: [1], 1: [0, 2], 2: [1, 3, 4], 3: [2], 4: [2, 5, 6], 5: [4], 6: [4]}
    """
    index       = 0          # numbering of inner vertices
    word        = [('in', 0)]  # initial vertex is inner
    leaf_stack  = [0, 0]     # stack of leaves still to be created
    inner_stack = [0]        # stack of active inner nodes
    edges = []
    embedding = {0: []}  # records the planar embedding of the tree
    for x in w:
        if x == 1:  # going up to a new inner vertex
            index += 1
            embedding[index] = inner_stack[-1:]
            embedding[inner_stack[-1]].append(index)
            leaf_stack.extend([index, index])
            inner_stack.append(index)
            edges.append(inner_stack[-2:])
            word.append(('in', index))
        else:
            if leaf_stack and inner_stack[-1] == leaf_stack[-1]:  # up and down to a new leaf
                leaf_stack.pop()
                word.extend([('lf', inner_stack[-1]), ('in', inner_stack[-1])])
            else:  # going down to a known inner vertex
                inner_stack.pop()
                word.append(('in', inner_stack[-1]))
    G = Graph(edges, format='list_of_edges')
    G.set_embedding(embedding)
    return word[:-1], G


def RandomTriangulation(n, set_position=False):
    r"""
    Return a random triangulation on `n` vertices.

    A triangulation is a planar graph all of whose faces are
    triangles (3-cycles).

    INPUT:

    - `n` -- an integer

    - ``set_position`` -- boolean (default ``False``) if set to ``True``, this
      will compute coordinates for a planar drawing of the graph.

    OUTPUT:

    A random triangulation chosen uniformly among the *rooted* triangulations on
    `n` vertices. This is a planar graph and comes with a combinatorial
    embedding.

    Because some triangulations have nontrivial automorphism
    groups, this may not be equal to the uniform distribution among unrooted
    triangulations.

    ALGORITHM:

    The algorithm is taken from [PS2006]_, section 2.1.

    Starting from a planar tree (represented by its contour as a
    sequence of vertices), one first performs local closures, until no
    one is possible. A local closure amounts to replace in the cyclic
    contour word a sequence ``in1,in2,in3,lf,in3`` by
    ``in1,in3``. After all local closures are done, one has reached
    the partial closure, as in [PS2006]_, figure 5 (a).

    Then one has to perform complete closure by adding two more
    vertices, in order to reach the situation of [PS2006]_, figure 5
    (b). For this, it is necessary to find inside the final contour
    one of the two subsequences ``lf,in,lf``.

    At every step of the algorithm, newly created edges are recorded
    in a graph, which will be returned at the end.

    The combinatorial embedding is also computed and recorded in the
    output graph.

    .. SEEALSO::

        :meth:`~sage.graphs.graph_generators.GraphGenerators.triangulations`,
        :func:`~sage.homology.examples.RandomTwoSphere`.

    EXAMPLES::

        sage: G = graphs.RandomTriangulation(6, True); G
        Graph on 6 vertices
        sage: G.is_planar()
        True
        sage: G.girth()
        3
        sage: G.plot(vertex_size=0, vertex_labels=False)
        Graphics object consisting of 13 graphics primitives

    TESTS::

        sage: G.get_embedding() is not None
        True
        sage: for i in range(10):
        ....:     g = graphs.RandomTriangulation(30)
        ....:     assert g.is_planar()
        sage: for i in range(10):
        ....:     g = graphs.RandomTriangulation(10)
        ....:     assert g.is_planar(on_embedding=g.get_embedding())

    REFERENCES:

    .. [PS2006] Dominique Poulalhon and Gilles Schaeffer,
       *Optimal coding and sampling of triangulations*,
       Algorithmica 46 (2006), no. 3-4, 505-527,
       http://www.lix.polytechnique.fr/~poulalho/Articles/PoSc_Algorithmica06.pdf

    """
    if n < 3:
        raise ValueError('only defined for n >= 3')
    w = _auxiliary_random_word(n - 2)
    word, graph = _contour_and_graph_from_word(w)
    edges = []

    embedding = graph.get_embedding()

    # 'partial closures' described in 2.1 of [PS2006]_.
    pattern = ['in', 'in', 'in', 'lf', 'in']

    def rotate_word_to_next_occurrence(word):
        """
        Rotate ``word`` so that the given pattern occurs at the beginning.

        If the given pattern is not found, return the empty list.
        """
        N = len(word)
        for i in range(N):
            if all(word[(i + j) % N][0] == pattern[j] for j in range(5)):
                return word[i:] + word[:i]
        return []

    # We greedily perform the replacements 'in1,in2,in3,lf,in3'->'in1,in3'.
    while True:
        word2 = rotate_word_to_next_occurrence(word)
        if len(word2) >= 5:
            word = [word2[0]] + word2[4:]
            in1, in2, in3 = [u[1] for u in word2[:3]]
            edges.append([in1, in3])  # edge 'in1,in3'
            idx = embedding[in1].index(in2)
            embedding[in1].insert(idx, in3)
            idx = embedding[in3].index(in2)
            embedding[in3].insert(idx + 1, in1)
        else:
            break

    graph.add_edges(edges)
    # This is the end of partial closure.

    # There remains to add two new vertices a and b.
    a = -1
    b = -2
    graph.add_edge((a, b))

    # Every remaining 'lf' vertex is linked either to a or to b.
    # Switching a/b happens when one meets the sequence 'lf','in','lf'.
    a_or_b = a
    embedding[a] = []
    embedding[b] = []
    last_lf_occurrence = -42
    change = {}
    for x in word:
        last_lf_occurrence -= 1
        if x[0] == 'lf':
            if last_lf_occurrence == -2:
                change[a_or_b] = x[1]
                a_or_b = b if a_or_b == a else a
            graph.add_edge((a_or_b, x[1]))
            embedding[a_or_b].insert(0, x[1])
            last_lf_occurrence = 0

    # conjugates the embeddings of a and b
    # in a way that helps to complete the embedding
    for a_or_b in [a, b]:
        emba = embedding[a_or_b]
        idx = emba.index(change[a_or_b])
        embedding[a_or_b] = emba[idx:] + emba[:idx]
    embedding[a].append(b)
    embedding[b].append(a)

    # completes the embedding by inserting missing half-edges
    for a_or_b in [a, b]:
        emb = embedding[a_or_b]
        for i, v in enumerate(emb[:-1]):
            if i == 0:
                embedding[v].insert(embedding[v].index(emb[1]) + 1, a_or_b)
            else:
                embedding[v].insert(embedding[v].index(emb[i - 1]), a_or_b)

    assert graph.num_edges() == 3 * (n - 2)
    assert graph.num_verts() == n

    graph.set_embedding(embedding)

    if set_position:
        graph.layout(layout="planar", save_pos=True)

    return graph


def blossoming_contour(t, shift=0):
    """
    Return a random blossoming of a binary tree `t`, as a contour word.

    This is doing several things simultaneously:

    - complete the binary tree, by adding leaves labelled ``xb``,
    - add a vertex labelled ``n`` at the middle of every inner
      edge, with a leaf labelled ``x`` either on the left or on the
      right (at random),
    - number all vertices (but not leaves) by integers starting from `shift`,
    - compute the counter-clockwise contour word of the result.

    Initial vertices receive the label ``i``.

    This is an auxiliary function, used for the generation of random
    planar bicubic maps.

    INPUT:

    - `t` -- a binary tree (non-empty)

    - ``shift`` -- an integer (default `0`), used as a starting index

    OUTPUT:

    contour word of a random blossoming of `t`

    EXAMPLES::

        sage: from sage.graphs.generators.random import blossoming_contour
        sage: print(blossoming_contour(BinaryTrees(1).an_element()))
        [('i', 0), ('xb',), ('i', 0), ('xb',), ('i', 0)]

        sage: t = BinaryTrees(2).random_element()
        sage: print(blossoming_contour(t))  # random
        [('i', 0), ('xb',), ('i', 0), ('n', 2), ('i', 1), ('xb',), ('i', 1),
        ('xb',), ('i', 1), ('n', 2), ('x',), ('n', 2), ('i', 0)]

        sage: w = blossoming_contour(BinaryTrees(3).random_element()); len(w)
        21
        sage: w.count(('xb',))
        4
        sage: w.count(('x',))
        2

    TESTS::

        sage: from sage.graphs.generators.random import blossoming_contour
        sage: blossoming_contour(BinaryTrees(0).an_element())
        Traceback (most recent call last):
        ...
        ValueError: tree must be non-empty
    """
    if not t:
        raise ValueError('tree must be non-empty')
    t1, t2 = t
    leaf_xb = ('xb',)
    leaf_x = ('x',)
    n1 = t1.node_number()
    n = t.node_number()

    # adding buds on edges in t1
    if not t1:
        tt1 = [leaf_xb]
    elif randint(0, 1):
        label1 = ('n', shift)
        tt1 = [label1, leaf_x, label1] + blossoming_contour(t1, shift + 1)
        tt1 += [label1]
    else:
        label1 = ('n', shift + 2 * n1 - 1)
        tt1 = [label1] + blossoming_contour(t1, shift)
        tt1 += [label1, leaf_x, label1]

    # adding buds on edges in t2
    if not t2:
        tt2 = [leaf_xb]
    elif randint(0, 1):
        label2 = ('n', shift + 2 * n1 + 1)
        tt2 = [label2, leaf_x, label2]
        tt2 += blossoming_contour(t2, shift + 2 * n1 + 2) + [label2]
    else:
        label2 = ('n', shift + 2 * n - 2)
        tt2 = [label2] + blossoming_contour(t2, shift + 2 * n1 + 1)
        tt2 += [label2, leaf_x, label2]

    label = [('i', shift + 2 * n1)]
    return label + tt1 + label + tt2 + label


def RandomBicubicPlanar(n):
    """
    Return the graph of a random bipartite cubic map with `3 n` edges.

    INPUT:

    `n` -- an integer (at least `1`)

    OUTPUT:

    a graph with multiple edges (no embedding is provided)

    The algorithm used is described in [Schaeffer99]_. This samples
    a random rooted bipartite cubic map, chosen uniformly at random.

    First one creates a random binary tree with `n` vertices. Next one
    turns this into a blossoming tree (at random) and reads the
    contour word of this blossoming tree.

    Then one performs a rotation on this word so that this becomes a
    balanced word. There are three ways to do that, one is picked at
    random. Then a graph is build from the balanced word by iterated
    closure (adding edges).

    In the returned graph, the three edges incident to any given
    vertex are colored by the integers 0, 1 and 2.

    .. SEEALSO:: the auxiliary method :func:`blossoming_contour`

    EXAMPLES::

        sage: n = randint(200, 300)
        sage: G = graphs.RandomBicubicPlanar(n)
        sage: G.order() == 2*n
        True
        sage: G.size() == 3*n
        True
        sage: G.is_bipartite() and G.is_planar() and G.is_regular(3)
        True
        sage: dic = {'red':[v for v in G.vertices() if v[0] == 'n'],
        ....:        'blue': [v for v in G.vertices() if v[0] != 'n']}
        sage: G.plot(vertex_labels=False,vertex_size=20,vertex_colors=dic)
        Graphics object consisting of ... graphics primitives

    .. PLOT::
        :width: 300 px

        G = graphs.RandomBicubicPlanar(200)
        V0 = [v for v in G.vertices() if v[0] == 'n']
        V1 = [v for v in G.vertices() if v[0] != 'n']
        dic = {'red': V0, 'blue': V1}
        sphinx_plot(G.plot(vertex_labels=False,vertex_colors=dic))

    REFERENCES:

    .. [Schaeffer99] Gilles Schaeffer, *Random Sampling of Large Planar Maps and Convex Polyhedra*,
       Annual ACM Symposium on Theory of Computing (Atlanta, GA, 1999)
    """
    from sage.combinat.binary_tree import BinaryTrees
    from sage.rings.finite_rings.integer_mod_ring import Zmod
    if not n:
        raise ValueError("n must be at least 1")
    # first pick a random binary tree
    t = BinaryTrees(n).random_element()

    # next pick a random blossoming of this tree, compute its contour
    contour = blossoming_contour(t) + [('xb',)]   # adding the final xb

    # first step : rotate the contour word to one of 3 balanced
    N = len(contour)
    double_contour = contour + contour
    pile = []
    not_touched = [i for i in range(N) if contour[i][0] in ['x', 'xb']]
    for i, w in enumerate(double_contour):
        if w[0] == 'x' and i < N:
            pile.append(i)
        elif w[0] == 'xb' and (i % N) in not_touched:
            if pile:
                j = pile.pop()
                not_touched.remove(i % N)
                not_touched.remove(j)

    # random choice among 3 possibilities for a balanced word
    idx = not_touched[randint(0, 2)]
    w = contour[idx + 1:] + contour[:idx + 1]

    # second step : create the graph by closure from the balanced word
    G = Graph(multiedges=True)

    pile = []
    Z3 = Zmod(3)
    colour = Z3.zero()
    not_touched = [i for i, v in enumerate(w) if v[0] in ['x', 'xb']]
    for i, v in enumerate(w):
        # internal edges
        if v[0] == 'i':
            colour += 1
            if w[i + 1][0] == 'n':
                G.add_edge((w[i], w[i + 1], colour))
        elif v[0] == 'n':
            colour += 2
        elif v[0] == 'x':
            pile.append(i)
        elif v[0] == 'xb' and i in not_touched:
            if pile:
                j = pile.pop()
                G.add_edge((w[i + 1], w[j - 1], colour))
                not_touched.remove(i)
                not_touched.remove(j)

    # there remains to add three edges to elements of "not_touched"
    # from a new vertex labelled "n"
    for i in not_touched:
        taken_colours = [edge[2] for edge in G.edges_incident(w[i - 1])]
        colour = [u for u in Z3 if u not in taken_colours][0]
        G.add_edge((('n', -1), w[i - 1], colour))

    return G
