r"""
Random graphs

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

import sys
# import from Sage library
from sage.graphs.graph import Graph
from sage.misc.randstate import current_randstate
from sage.misc.randstate import set_random_seed
from sage.misc.prandom import random
from sage.misc.prandom import randint


def RandomGNP(n, p, seed=None, fast=True, algorithm='Sage'):
    r"""
    Returns a random graph on `n` nodes. Each edge is inserted independently
    with probability `p`.

    INPUT:

    - ``n`` -- number of nodes of the graph

    - ``p`` -- probability of an edge

    - ``seed`` -- a ``random.Random`` seed or a Python ``int`` for the random
      number generator (default: ``None``)

    - ``fast`` -- boolean set to True (default) to use the algorithm with
      time complexity in `O(n+m)` proposed in [BB2005a]_. It is designed
      for generating large sparse graphs. It is faster than other algorithms for
      *LARGE* instances (try it to know whether it is useful for you).

    - ``algorithm`` -- By default (``algorithm='Sage'``), this function uses the
      algorithm implemented in ```sage.graphs.graph_generators_pyx.pyx``. When
      ``algorithm='networkx'``, this function calls the NetworkX function
      ``fast_gnp_random_graph``, unless ``fast=False``, then
      ``gnp_random_graph``. Try them to know which algorithm is the best for
      you. The ``fast`` parameter is not taken into account by the 'Sage'
      algorithm so far.

    REFERENCES:

    - [ER1959]_

    - [Gil1959]_

    PLOTTING: When plotting, this graph will use the default spring-layout
    algorithm, unless a position dictionary is specified.

    EXAMPLES: We show the edge list of a random graph on 6 nodes with
    probability `p = .4`::

        sage: set_random_seed(0)
        sage: graphs.RandomGNP(6, .4).edges(sort=true, labels=False)
        [(0, 3), (1, 2), (2, 3), (2, 4)]

    We plot a random graph on 12 nodes with probability `p = .71`::

        sage: gnp = graphs.RandomGNP(12,.71)
        sage: gnp.show()                        # long time                             # needs sage.plot

    We view many random graphs using a graphics array::

        sage: g = []
        sage: j = []
        sage: for i in range(9):
        ....:     k = graphs.RandomGNP(i+3,.43)
        ....:     g.append(k)
        sage: for i in range(3):                                                        # needs sage.plot
        ....:     n = []
        ....:     for m in range(3):
        ....:         n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ....:     j.append(n)
        sage: G = graphics_array(j)                                                     # needs sage.plot
        sage: G.show()                          # long time                             # needs sage.plot
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
        sage: graphs.RandomGNP(50,.2, algorithm="networkx").size()                      # needs networkx
        279     # 32-bit
        209     # 64-bit
    """
    if n < 0:
        raise ValueError("The number of nodes must be positive or null.")
    if 0.0 > p or 1.0 < p:
        raise ValueError("The probability p must be in [0..1].")

    if p == 1:
        from sage.graphs.generators.basic import CompleteGraph
        return CompleteGraph(n)

    if algorithm == 'networkx':
        if seed is None:
            seed = int(current_randstate().long_seed() % sys.maxsize)
        import networkx
        if fast:
            G = networkx.fast_gnp_random_graph(n, p, seed=seed)
        else:
            G = networkx.gnp_random_graph(n, p, seed=seed)
        return Graph(G)
    elif algorithm in ['Sage', 'sage']:
        # We use the Sage generator
        from sage.graphs.graph_generators_pyx import RandomGNP as sageGNP
        return sageGNP(n, p, seed=seed)
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

        sage: G = graphs.RandomBarabasiAlbert(6,2)                                      # needs networkx
        sage: G.order(), G.size()                                                       # needs networkx
        (6, 8)
        sage: G.degree_sequence()  # random                                             # needs networkx
        [4, 3, 3, 2, 2, 2]

    We plot a random graph on 12 nodes with `m = 3`::

        sage: ba = graphs.RandomBarabasiAlbert(12,3)                                    # needs networkx
        sage: ba.show()                         # long time                             # needs networkx sage.plot

    We view many random graphs using a graphics array::

        sage: # needs networkx sage.plot
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
        sage: G = graphics_array(j)
        sage: G.show()                          # long time

    When `m = 1`, the generated graph is a tree::

        sage: graphs.RandomBarabasiAlbert(6, 1).is_tree()                               # needs networkx
        True
    """
    if seed is None:
        seed = int(current_randstate().long_seed() % sys.maxsize)
    import networkx
    return Graph(networkx.barabasi_albert_graph(int(n), int(m), seed=seed))


def RandomBipartite(n1, n2, p, set_position=False, seed=None):
    r"""
    Returns a bipartite graph with `n1+n2` vertices such that any edge
    from `[n1]` to `[n2]` exists with probability `p`.

    INPUT:

    - ``n1, n2`` -- Cardinalities of the two sets

    - ``p`` -- Probability for an edge to exist

    - ``set_position`` -- boolean (default ``False``); if set to ``True``, we
      assign positions to the vertices so that the set of cardinality `n1` is
      on the line `y=1` and the set of cardinality `n2` is on the line `y=0`.

    - ``seed`` -- a ``random.Random`` seed or a Python ``int`` for the random
      number generator (default: ``None``)

    EXAMPLES::

        sage: g = graphs.RandomBipartite(5, 2, 0.5)                                     # needs numpy
        sage: g.vertices(sort=True)                                                     # needs numpy
        [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (1, 0), (1, 1)]

    TESTS::

        sage: g = graphs.RandomBipartite(5, -3, 0.5)                                    # needs numpy
        Traceback (most recent call last):
        ...
        ValueError: n1 and n2 should be integers strictly greater than 0
        sage: g = graphs.RandomBipartite(5, 3, 1.5)                                     # needs numpy
        Traceback (most recent call last):
        ...
        ValueError: parameter p is a probability, and so should be a real value between 0 and 1

    :trac:`12155`::

        sage: graphs.RandomBipartite(5, 6, .2).complement()                             # needs numpy
        complement(Random bipartite graph of order 5+6 with edge probability 0.200000000000000): Graph on 11 vertices

    Test assigned positions::

        sage: # needs numpy
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
    if seed is not None:
        set_random_seed(seed)

    from numpy.random import uniform

    g = Graph(name=f"Random bipartite graph of order {n1}+{n2} with edge probability {p}")

    S1 = [(0, i) for i in range(n1)]
    S2 = [(1, i) for i in range(n2)]
    g.add_vertices(S1)
    g.add_vertices(S2)

    for w in range(n2):
        for v in range(n1):
            if uniform() <= p:
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


def RandomRegularBipartite(n1, n2, d1, set_position=False, seed=None):
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

    - ``seed`` -- a ``random.Random`` seed or a Python ``int`` for the random
      number generator (default: ``None``)

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
    if seed is not None:
        set_random_seed(seed)

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
        E = E.symmetric_difference(CompleteBipartiteGraph(n1, n2).edges(sort=False, labels=False))
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
        G._line_embedding(list(range(n1, n1 + n2)), first=(0, 0), last=(nmax, 0))

    return G


def RandomBlockGraph(m, k, kmax=None, incidence_structure=False, seed=None):
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

    - ``seed`` -- a ``random.Random`` seed or a Python ``int`` for the random
      number generator (default: ``None``)

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
        sage: IncidenceStructure(IS)                                                    # needs sage.modules
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
    if seed is not None:
        set_random_seed(seed)

    if m == 1:
        # A block graph with a single block is a clique
        IS = [list(range(randint(k, kmax)))]

    elif kmax == 2:
        # A block graph with blocks of order 2 is a tree
        IS = [list(e) for e in RandomTree(m + 1).edges(sort=True, labels=False)]

    else:
        # We start with a random tree of order m
        T = RandomTree(m)

        # We create a block of order in range [k,kmax] per vertex of the tree
        B = {u: [(u, i) for i in range(randint(k, kmax))] for u in T}

        # For each edge of the tree, we choose 1 vertex in each of the
        # corresponding blocks and we merge them. We use a disjoint set data
        # structure to keep a unique identifier per merged vertices
        DS = DisjointSet([i for u in B for i in B[u]])
        for u, v in T.edges(sort=True, labels=0):
            DS.union(choice(B[u]), choice(B[v]))

        # We relabel vertices in the range [0, m*(k-1)] and build the incidence
        # structure
        new_label = {root: i for i, root in enumerate(DS.root_to_elements_dict())}
        IS = [[new_label[DS.find(v)] for v in B[u]] for u in B]

    if incidence_structure:
        return IS

    # We finally build the block graph
    if k == kmax:
        BG = Graph(name="Random Block Graph with {} blocks of order {}".format(m, k))
    else:
        BG = Graph(name="Random Block Graph with {} blocks of order {} to {}".format(m, k, kmax))
    for block in IS:
        BG.add_clique(block)
    return BG


def RandomBoundedToleranceGraph(n, seed=None):
    r"""
    Return a random bounded tolerance graph.

    The random tolerance graph is built from a random bounded tolerance
    representation by using the function `ToleranceGraph`. This representation
    is a list `((l_0,r_0,t_0), (l_1,r_1,t_1), ..., (l_k,r_k,t_k))` where `k =
    n-1` and `I_i = (l_i,r_i)` denotes a random interval and `t_i` a random
    positive value less than or equal to the length of the interval `I_i`. The
    width of the representation is limited to `n^2 * 2^n`.

    .. NOTE::

        The tolerance representation used to create the graph can
        be recovered using ``get_vertex()`` or ``get_vertices()``.

    INPUT:

    - ``n`` -- number of vertices of the random graph.

    - ``seed`` -- a ``random.Random`` seed or a Python ``int`` for the random
      number generator (default: ``None``)

    EXAMPLES:

    Every (bounded) tolerance graph is perfect. Hence, the
    chromatic number is equal to the clique number ::

        sage: g = graphs.RandomBoundedToleranceGraph(8)
        sage: g.clique_number() == g.chromatic_number()
        True

    TESTS:

    Check that :trac:`32186` is fixed::

        sage: for _ in range(100): _ = graphs.RandomBoundedToleranceGraph(1)

    Check input parameter::

        sage: g = graphs.RandomToleranceGraph(-2)
        Traceback (most recent call last):
        ...
        ValueError: the number `n` of vertices must be >= 0
    """
    if n < 0:
        raise ValueError('the number `n` of vertices must be >= 0')
    if seed is not None:
        set_random_seed(seed)

    from sage.graphs.generators.intersection import ToleranceGraph

    W = n ** 2 * 2 ** n
    tolrep = []
    for _ in range(n):
        left = randint(0, W - 1)
        right = randint(0, W)
        if left >= right:
            left, right = right, left + 1
        tolrep.append((left, right, randint(1, right - left)))

    return ToleranceGraph(tolrep)


def RandomGNM(n, m, dense=False, seed=None):
    """
    Returns a graph randomly picked out of all graphs on n vertices
    with m edges.

    INPUT:

    - ``n`` -- number of vertices.

    - ``m`` -- number of edges.

    - ``dense`` -- whether to use NetworkX's
      :func:`dense_gnm_random_graph` or :func:`gnm_random_graph`

    - ``seed`` -- a ``random.Random`` seed or a Python ``int`` for the random
      number generator (default: ``None``)

    EXAMPLES:

    We show the edge list of a random graph on 5 nodes with 10 edges::

        sage: graphs.RandomGNM(5, 10).edges(sort=True, labels=False)                    # needs networkx
        [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]

    We plot a random graph on 12 nodes with m = 12::

        sage: gnm = graphs.RandomGNM(12, 12)                                            # needs networkx
        sage: gnm.show()                        # long time                             # needs networkx sage.plot

    We view many random graphs using a graphics array::

        sage: # needs networkx sage.plot
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
        sage: G = graphics_array(j)
        sage: G.show()                          # long time
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

        sage: # needs networkx
        sage: G = graphs.RandomNewmanWattsStrogatz(7, 2, 0.2)
        sage: G.order()
        7
        sage: C7 = graphs.CycleGraph(7)
        sage: G.subgraph_search(C7)
        Subgraph of (): Graph on 7 vertices
        sage: G.diameter() <= C7.diameter()
        True

    ::

        sage: G = graphs.RandomNewmanWattsStrogatz(12, 2, .3)                           # needs networkx
        sage: G.show()                          # long time                             # needs networkx sage.plot

    TESTS:

    We check that when `k = 2` and `p = 0`, the generated graph is a cycle::

        sage: G = graphs.RandomNewmanWattsStrogatz(7, 2, 0)                             # needs networkx
        sage: G.is_cycle()                                                              # needs networkx
        True

    We check that when `k = 4` and `p = 0`, the generated graph is a circulant
    graph of parameters ``[1, 2]``::

        sage: G = graphs.RandomNewmanWattsStrogatz(7, 4, 0)                             # needs networkx
        sage: G.is_isomorphic(graphs.CirculantGraph(7, [1, 2]))                         # needs networkx
        True

    REFERENCE:

    [NWS2002]_
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

    EXAMPLES::

        sage: G = graphs.RandomHolmeKim(12, 3, .3)                                      # needs networkx
        sage: G.show()                          # long time                             # needs networkx sage.plot

    REFERENCE:

    [HK2002a]_
    """
    if seed is None:
        seed = int(current_randstate().long_seed() % sys.maxsize)
    import networkx
    return Graph(networkx.powerlaw_cluster_graph(n, m, p, seed=seed))


def RandomIntervalGraph(n, seed=None):
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

    This definitions follows [BF2001]_.

    .. NOTE::

        The vertices are named 0, 1, 2, and so on. The intervals
        used to create the graph are saved with the graph and can
        be recovered using ``get_vertex()`` or ``get_vertices()``.

    INPUT:

    - ``n`` -- integer; the number of vertices in the random graph

    - ``seed`` -- a ``random.Random`` seed or a Python ``int`` for the random
      number generator (default: ``None``)

    EXAMPLES:

    As for any interval graph, the chromatic number is equal to
    the clique number ::

        sage: g = graphs.RandomIntervalGraph(8)
        sage: g.clique_number() == g.chromatic_number()
        True
    """
    if seed is not None:
        set_random_seed(seed)
    from sage.misc.prandom import random
    from sage.graphs.generators.intersection import IntervalGraph

    intervals = [tuple(sorted((random(), random()))) for i in range(n)]
    return IntervalGraph(intervals, True)


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
    from sage.misc.prandom import choice
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
                z = choice(tuple(neighbors))
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
        sage: S = connecting_nodes(T, 5)                                                # needs numpy
        sage: len(S)                                                                    # needs numpy
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
      chosen in `[0..1]`.

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


def RandomChordalGraph(n, algorithm="growing", k=None, l=None, f=None, s=None, seed=None):
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
      fraction. This value must be chosen in `[0..1]`. If not specified, this
      parameter is set to `\frac{1}{n-1}`.
      This parameter is used only when ``algorithm="pruned"``.
      See :meth:`~sage.graphs.generators.random.pruned_tree` for more details.

    - ``s`` -- a real number between 0 and 1 (default: ``None``); selection
      barrier for the size of trees. If not specified, this parameter is set to
      `0.5`. This parameter is used only when ``algorithm="pruned"``.
      See :meth:`~sage.graphs.generators.random.pruned_tree` for more details.

    - ``seed`` -- a ``random.Random`` seed or a Python ``int`` for the random
      number generator (default: ``None``)

    EXAMPLES::

        sage: from sage.graphs.generators.random import RandomChordalGraph
        sage: T = RandomChordalGraph(20, algorithm="growing", k=5)
        sage: T.is_chordal()
        True
        sage: T = RandomChordalGraph(20, algorithm="connecting", l=3)                   # needs numpy
        sage: T.is_chordal()                                                            # needs numpy
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

    if seed is not None:
        set_random_seed(seed)

    # 1. Generate a random tree of order n
    T = RandomTree(n)

    # 2. Generate n non-empty subtrees of T: {T1,...,Tn}
    if algorithm == "growing":
        if k is None:
            from sage.misc.functional import isqrt
            k = isqrt(n)
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
    for i, s in enumerate(S):
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


    EXAMPLES:

    We check a random graph with 12 backbone
    nodes and probabilities `p = 0.7` and `q = 0.3`::

        sage: # needs networkx
        sage: G = graphs.RandomLobster(12, 0.7, 0.3)
        sage: leaves = [v for v in G.vertices(sort=False) if G.degree(v) == 1]
        sage: G.delete_vertices(leaves)                                 # caterpillar
        sage: leaves = [v for v in G.vertices(sort=False) if G.degree(v) == 1]
        sage: G.delete_vertices(leaves)                                 # path
        sage: s = G.degree_sequence()
        sage: if G:
        ....:     if G.num_verts() == 1:
        ....:         assert s == [0]
        ....:     else:
        ....:         assert s[-2:] == [1, 1]
        ....:     assert all(d == 2 for d in s[:-2])

    ::

        sage: G = graphs.RandomLobster(9, .6, .3)                                       # needs networkx
        sage: G.show()                          # long time                             # needs networkx sage.plot
    """
    if seed is None:
        seed = int(current_randstate().long_seed() % sys.maxsize)
    import networkx
    return Graph(networkx.random_lobster(n, p, q, seed=seed))


def RandomTree(n, seed=None):
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

    -  ``n`` -- number of vertices in the tree

    - ``seed`` -- a ``random.Random`` seed or a Python ``int`` for the random
      number generator (default: ``None``)

    EXAMPLES::

        sage: G = graphs.RandomTree(10)
        sage: G.is_tree()
        True
        sage: G.show()                          # long time                             # needs sage.plot

    TESTS:

    Ensuring that we encounter no unexpected surprise ::

        sage: all( graphs.RandomTree(10).is_tree()
        ....:      for i in range(100) )
        True

    Random tree with one and zero vertices::

        sage: graphs.RandomTree(0)
        Graph on 0 vertices
        sage: graphs.RandomTree(1)
        Graph on 1 vertex
    """
    g = Graph(n)
    if n <= 1:
        return g

    if seed is not None:
        set_random_seed(seed)

    # create random Prufer code
    code = [randint(0, n - 1) for i in range(n - 2)]

    # We count the number of symbols of each type.
    # count[k] is the number of times k appears in code
    #
    # (count[k] is set to -1 when the corresponding vertex is not
    # available anymore)
    count = [0] * n
    for k in code:
        count[k] += 1

    # We use a heap to store vertices for which count[k] == 0 and get the vertex
    # with smallest index
    from heapq import heapify, heappop, heappush
    zeros = [x for x in range(n) if not count[x]]
    heapify(zeros)

    for s in code:
        x = heappop(zeros)
        g.add_edge(x, s)
        count[x] = -1
        count[s] -= 1
        if not count[s]:
            heappush(zeros, s)

    # Adding as an edge the last two available vertices
    g.add_edge(zeros)

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

        sage: G = graphs.RandomTreePowerlaw(10, 3)                                      # needs networkx
        sage: G.is_tree()                                                               # needs networkx
        True
        sage: G.order(), G.size()                                                       # needs networkx
        (10, 9)

    ::

        sage: G = graphs.RandomTreePowerlaw(15, 2)                                      # needs networkx
        sage: if G:                             # random output         # long time, needs networkx sage.plot
        ....:     G.show()
    """
    if seed is None:
        seed = int(current_randstate().long_seed() % sys.maxsize)
    import networkx
    try:
        return Graph(networkx.random_powerlaw_tree(n, gamma, seed=seed, tries=tries))
    except networkx.NetworkXError:
        return False


def RandomKTree(n, k, seed=None):
    r"""
    Return a random `k`-tree on `n` nodes numbered `0` through `n-1`.

    ALGORITHM:

    The algorithm first generates a complete graph on `k + 1` vertices.
    Vertices are subsequently generated by randomly choosing one of the
    existing cliques in the graph, and creating a new clique by replacing
    one of the vertices in the selected clique with a newly created one.

    INPUT:

    - ``n`` -- number of vertices in the `k`-tree

    - ``k`` -- within a clique each vertex is connected to `k` vertices. `k`
      also corresponds to the treewidth of the `k`-tree

    - ``seed`` -- a ``random.Random`` seed or a Python ``int`` for the random
      number generator (default: ``None``)

    TESTS::

        sage: g=graphs.RandomKTree(50,5)
        sage: g.size()
        235
        sage: g.order()
        50
        sage: g.treewidth()
        5
    
    EXAMPLES::

    A (random) `k`-tree has treewidth `k`::

        sage: G = graphs.RandomKTree(50, 5)
        sage: G.treewidth()
        5
        sage: G.show()  # not tested 
    """    
    if n < k + 1:
        raise ValueError("n must be greater than k + 1")

    if seed is not None:
        set_random_seed(seed)

    g = Graph(name=f"Random {k}-tree of order {n}")
    g.add_clique(list(range(k + 1)))
    
    cliques = [list(range(k+1))]

    # Randomly choose a row, and copy 1 of the cliques
    # One of those vertices is then replaced with a new vertex 
    for newVertex in range(k + 1, n):
        copiedClique = cliques[randint(0, len(cliques)-1)].copy()
        copiedClique[randint(0, k)] = newVertex 
        cliques.append(copiedClique)
        for u in copiedClique:
            if u != newVertex:
                g.add_edge(u, newVertex)
    return g 


def RandomPartialKTree(n, k, x, seed=None):
    r"""
    Return a random partial `k`-tree on `n` nodes. 

    A partial `k`-tree is defined as a subgraph of a `k`-tree. This can also be
    described as a graph with treewidth at most `k`.

    INPUT:

    - ``n`` -- number of vertices in the `k`-tree

    - ``k`` -- within a clique each vertex is connected to `k` vertices. `k`
      also corresponds to the treewidth of the `k`-tree

    - ``x`` -- how many edges are deleted from the `k`-tree

    - ``seed`` -- a ``random.Random`` seed or a Python ``int`` for the random
      number generator (default: ``None``)

    TESTS::

        sage: g=graphs.RandomPartialKTree(50,5,5)
        sage: g.order()
        50
        sage: g.size()
        230
        sage: g.treewidth()
        5

    EXAMPLES::

        sage: G = graphs.RandomPartialKTree(50,5,5)
        sage: G.treewidth()
        5
        sage: G.show()  # not tested
    """
    if seed is not None:
        set_random_seed(seed)

    g = RandomKTree(n, k, seed)

    # Check that x doesn't delete too many edges
    # This formula calculates how many edges are in `k`-tree with `n` nodes
    if x > (k^2 + k) / 2 + (n - k - 1) * k:
        raise ValueError("x must be less than the number of edges in the `k`-tree with `n` nodes")

    from sage.misc.prandom import shuffle

    edges = list(g.edges())
    # Deletes x random edges from the graph
    shuffle(edges)
    g.delete_edges(edges[:x])
    return g


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

        sage: G = graphs.RandomRegular(3, 8)                                            # needs networkx
        sage: G.is_regular(k=3)                                                         # needs networkx
        True
        sage: G.degree_histogram()                                                      # needs networkx
        [0, 0, 0, 8]

    ::

        sage: G = graphs.RandomRegular(3, 20)                                           # needs networkx
        sage: if G:                             # random output         # long time, needs networkx sage.plot
        ....:     G.show()

    REFERENCES:

    - [KV2003]_

    - [SW1999]_
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

        sage: G = graphs.RandomShell([(10,20,0.8),(20,40,0.8)])                         # needs networkx
        sage: G.order(), G.size()                                                       # needs networkx
        (30, 52)
        sage: G.show()                          # long time                             # needs networkx sage.plot
    """
    if seed is None:
        seed = int(current_randstate().long_seed() % sys.maxsize)
    import networkx
    return Graph(networkx.random_shell_graph(constructor, seed=seed))


def RandomToleranceGraph(n, seed=None):
    r"""
    Return a random tolerance graph.

    The random tolerance graph is built from a random tolerance representation
    by using the function `ToleranceGraph`. This representation is a list
    `((l_0,r_0,t_0), (l_1,r_1,t_1), ..., (l_k,r_k,t_k))` where `k = n-1` and
    `I_i = (l_i,r_i)` denotes a random interval and `t_i` a random positive
    value. The width of the representation is limited to `n^2 * 2^n`.

    .. NOTE::

        The vertices are named 0, 1, ..., n-1. The tolerance representation used
        to create the graph is saved with the graph and can be recovered using
        ``get_vertex()`` or ``get_vertices()``.

    INPUT:

    - ``n`` -- number of vertices of the random graph

    - ``seed`` -- a ``random.Random`` seed or a Python ``int`` for the random
      number generator (default: ``None``)

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
        ValueError: the number `n` of vertices must be >= 0
    """
    from sage.graphs.generators.intersection import ToleranceGraph

    if n < 0:
        raise ValueError('the number `n` of vertices must be >= 0')
    if seed is not None:
        set_random_seed(seed)

    W = n**2 * 2**n

    tolrep = []
    for _ in range(n):
        left = randint(0, W)
        right = randint(0, W)
        if left > right:
            left, right = right, left
        # The tolerance value must be > 0
        tolrep.append((left, right, randint(1, W)))

    g = ToleranceGraph(tolrep)
    g.name("Random tolerance graph")
    return g


# uniform random triangulation using Schaeffer-Poulalhon algorithm

def _auxiliary_random_forest_word(n, k):
    r"""
    Return a random word used to generate random triangulations.

    INPUT:

    - ``n`` -- an integer

    - ``k`` -- an integer

    OUTPUT:

    A binary sequence `w` of length `4n+2k-4` with `n` ones, such that any
    proper prefix `u` of `w` satisfies `3|u|_1 - |u|_0 \geq -2k+4` (where
    `|u|_1` and `|u|_0` are respectively the number of 1s and 0s in `u`). Those
    words are the expected input of :func:`_contour_and_graph_from_words`.

    ALGORITHM:

    A random word with these numbers of `0` and `1` plus one additional `0` is
    chosen. This word is then rotated such the prefix property is fulfilled for
    each proper prefix and only violated by the final `0` (which is deleted
    afterwards). There is exactly one such rotation (compare Section 4.3 in
    [PS2006]_).

    Let us consider a word `w` satisfying the expected conditions. By
    drawing a step `(1,3)` for each `1` and a step `(1,-1)` for each `0` in
    `w`, one gets a path starting at height `0`, ending at height `-2k+3`
    (before removing the final `0`) and staying above (or on) the horizontal
    line of height `-2k+4` except at the end point.

    Now consider an arbitrary word `w` with `n` ones and `3n+2k-3` zeros. By
    cutting the word at the first position of minimum height, let us write
    `w=uv`. One can then see that the word `vu` touches the line of height
    `-2k+3` only after the last step. Further one can see that this is the only
    rotation of the word `w` with this property.

    EXAMPLES::

        sage: from sage.graphs.generators.random import _auxiliary_random_forest_word
        sage: with(seed(94364165)):
        ....:     _auxiliary_random_forest_word(4, 3)
        ....:     _auxiliary_random_forest_word(3, 5)
        [1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0]
        [1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    TESTS::

        sage: def partial_sums(w):
        ....:     steps = {1: 3, 0: -1}
        ....:     curr_sum = 0
        ....:     for x in w:
        ....:         curr_sum += steps[x]
        ....:         yield curr_sum

        sage: for k in range(3,6):
        ....:     for n in range(k, 10):
        ....:         w = _auxiliary_random_forest_word(n, k)
        ....:         assert len(w) == 4*n + 2*k - 4
        ....:         assert w.count(1) == n
        ....:         for partial_sum in partial_sums(w):
        ....:             assert partial_sum >= -2*k + 4
    """
    from sage.misc.prandom import shuffle
    w = [0] * (3*n + 2*k - 3) + [1] * n
    shuffle(w)

    # Finding the admissible shift
    partial_sum = 0
    min_value = 0
    min_pos = 0
    for i, x in enumerate(w):
        if x:
            partial_sum += 3
        else:
            partial_sum -= 1
        if partial_sum < min_value:
            min_value = partial_sum
            min_pos = i
    return w[min_pos+1:] + w[:min_pos]


def _contour_and_graph_from_words(pendant_word, forest_word):
    r"""
    Return the contour word and the graph of inner vertices of the `k`-gonal
    forest associated with the words ``pendant_word`` and ``forest_word``.

    INPUT:

    - ``pendant_word`` -- a word with `k-1` zeros and `k-3` ones

    - ``forest_word`` -- a word in `0` and `1` as given by
      :func:`_auxiliary_random_word` with the parameter ``k`` set to the number
      of zeros in ``pendant_word`` plus `1`

    ``forest_word`` must satisfy the conditions hinted in Proposition 5.4 of
    [PS2006]_ (see :func:`_auxiliary_random_forest_word`).

    OUTPUT:

    a pair ``(seq, G)`` where:

    - ``seq`` is a sequence of pairs (label, integer) representing the
      contour walk along the `k`-gonal forest associated with the words
      ``pendant_word`` and ``forest_word``.

    - ``G`` is the `k`-gonal forest associated with the words ``pendant_word``
      and ``forest_word``.

    The underlying bijection from words to `k`-gonal forests is described in
    Section 5.1 of [PS2006]_. The ``pendant_word`` corresponds to the factor
    `\binom{2k-4}{k-3}` in the counting formula of Proposition 5.4 and the
    ``forest_word`` corresponds to the factor `\frac{2k-3}{3m+2k-3}
    \binom{4m+2k-4}{m}`.

    In the ``forest_word``, the letter `1` means going away from the root ("up")
    from an inner vertex to another inner vertex. The letter `0` denotes all
    other steps of the discovery, i.e. either discovering a leaf vertex or going
    toward the root ("down").

    Inner vertices are tagged with 'in' and leaves are tagged with
    'lf'. Inner vertices are moreover labelled by integers, and leaves
    by the label of the neighbor inner vertex.

    EXAMPLES::

        sage: from sage.graphs.generators.random import _contour_and_graph_from_words
        sage: seq, G = _contour_and_graph_from_words([0, 0], [1, 0, 0, 0, 0, 0])
        sage: seq
        [('in', 0),
         ('in', 3),
         ('lf', 3),
         ('in', 3),
         ('lf', 3),
         ('in', 3),
         ('in', 0),
         ('in', 1),
         ('in', 2)]
        sage: G
        Graph on 4 vertices

        sage: from sage.graphs.generators.random import _auxiliary_random_forest_word
        sage: _, G = _contour_and_graph_from_words([0, 1, 0, 0, 1, 0], _auxiliary_random_forest_word(20, 5)) # random
        sage: len(G.faces())
        2

        sage: longw = [1,1,0,1,0,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0]
        sage: _, G = _contour_and_graph_from_words([0, 0], longw)
        sage: G.get_embedding()
        {0: [1, 2, 3],
         1: [2, 0],
         2: [0, 1],
         3: [0, 4],
         4: [3, 5, 6],
         5: [4],
         6: [4, 7, 8],
         7: [6],
         8: [6]}
    """
    k = (len(pendant_word) + 4) // 2

    index = 0  # numbering of inner vertices
    word = [('in', 0)]  # the word representing the contour walk

    # start with the outer face, a cycle of length k
    edges = [[i, (i + 1) % k] for i in range(k)]
    embedding = {i: [(i + 1) % k, (i - 1 + k) % k] for i in range(k)}

    # add the pendant edges
    for x in pendant_word:
        if x:
            word.extend([('lf', index), ('in', index)])
        else:
            index += 1
            word.append(('in', index))

    # add trees
    curr_word_pos = 0
    curr_forest_word_pos = 0
    while curr_forest_word_pos < len(forest_word):
        x = forest_word[curr_forest_word_pos]
        # insert a tree at current position
        if x:
            index += 1
            embedding[index] = [word[curr_word_pos][1]]
            embedding[word[curr_word_pos][1]].append(index)
            edges.append([word[curr_word_pos][1], index])
            # stack of leaves still to be created
            leaf_stack = [index, index]
            # stack of active inner nodes
            inner_stack = [word[curr_word_pos][1], index]
            word.insert(curr_word_pos+1, ('in', index))
            curr_word_pos += 1
            while len(inner_stack) > 1:
                curr_forest_word_pos += 1
                x = forest_word[curr_forest_word_pos]
                if x:
                    index += 1
                    embedding[index] = inner_stack[-1:]
                    embedding[inner_stack[-1]].append(index)
                    leaf_stack.extend([index, index])
                    inner_stack.append(index)
                    edges.append(inner_stack[-2:])
                    word.insert(curr_word_pos+1, ('in', index))
                    curr_word_pos += 1
                else:
                    # up and down to a new leaf
                    if leaf_stack and inner_stack[-1] == leaf_stack[-1]:
                        leaf_stack.pop()
                        word.insert(curr_word_pos+1, ('lf', inner_stack[-1]))
                        word.insert(curr_word_pos+2, ('in', inner_stack[-1]))
                        curr_word_pos += 2
                    # going down to a known inner vertex
                    else:
                        inner_stack.pop()
                        word.insert(curr_word_pos+1, ('in', inner_stack[-1]))
                        curr_word_pos += 1
        # go to next insertion position
        else:
            curr_word_pos += 1
            if word[curr_word_pos][0] == 'lf':
                curr_word_pos += 1
        curr_forest_word_pos += 1

    G = Graph(edges, format='list_of_edges')
    G.set_embedding(embedding)
    return word, G


def RandomTriangulation(n, set_position=False, k=3, seed=None):
    r"""
    Return a random inner triangulation of an outer face of degree ``k`` with
    ``n`` vertices in total.

    An inner triangulation is a plane graph all of whose faces (except the
    outer/unbounded face) are triangles (3-cycles).

    INPUT:

    - ``n`` -- the number of vertices of the graph

    - ``k`` -- the size of the outer face

    - ``set_position`` -- boolean (default ``False``); if set to ``True``, this
      will compute coordinates for a planar drawing of the graph.

    - ``seed`` -- a ``random.Random`` seed or a Python ``int`` for the random
      number generator (default: ``None``)

    OUTPUT:

    A random graph chosen uniformly among the inner triangulations of a *rooted*
    `k`-gon with `n` vertices (including the `k` vertices from the outer face).
    This is a planar graph and comes with a combinatorial embedding. The
    vertices of the root edge are labelled ``-1`` and ``-2`` and the outer face
    is the face returned by :meth:`Graph.faces` in which ``-1`` and ``-2`` are
    consecutive vertices in this order.

    Because some triangulations have nontrivial automorphism
    groups, this may not be equal to the uniform distribution among inner
    triangulations of unrooted `k`-gons.

    ALGORITHM:

    The algorithm is taken from [PS2006]_, Section 5.

    Starting from a planar `k`-gonal forest (represented by its contour as a
    sequence of vertices), one performs local closures, until no
    one is possible. A local closure amounts to replace in the cyclic
    contour word a sequence ``in1, in2, in3, lf, in3`` by
    ``in1, in3``.

    At every step of the algorithm, newly created edges are recorded
    in a graph, which will be returned at the end.
    The combinatorial embedding is also computed and recorded in the
    output graph.

    .. SEEALSO::

        :meth:`~sage.graphs.graph_generators.GraphGenerators.triangulations`,
        :func:`~sage.topology.simplicial_complex_examples.RandomTwoSphere`.

    EXAMPLES::

        sage: G = graphs.RandomTriangulation(6, True); G
        Graph on 6 vertices
        sage: G.is_planar()
        True
        sage: G.girth()
        3
        sage: G.plot(vertex_size=0, vertex_labels=False)                                # needs sage.plot
        Graphics object consisting of 13 graphics primitives

        sage: H = graphs.RandomTriangulation(7, k=5)
        sage: sorted(len(f) for f in H.faces())
        [3, 3, 3, 3, 3, 3, 3, 5]

    TESTS::

        sage: G.get_embedding() is not None
        True

        sage: graphs.RandomTriangulation(3, k=4)
        Traceback (most recent call last):
        ...
        ValueError: The number 'n' of vertices must be at least the size 'k' of the outer face.
        sage: graphs.RandomTriangulation(3, k=2)
        Traceback (most recent call last):
        ...
        ValueError: The size 'k' of the outer face must be at least 3.

        sage: for i in range(10):
        ....:     g = graphs.RandomTriangulation(30) # random
        ....:     assert g.is_planar()
        sage: for k in range(3, 10):
        ....:     g = graphs.RandomTriangulation(10, k=k) # random
        ....:     assert g.is_planar(on_embedding=g.get_embedding())
    """
    if k < 3:
        raise ValueError("The size 'k' of the outer face must be at least 3.")
    if n < k:
        raise ValueError("The number 'n' of vertices must be at least the size "
                         "'k' of the outer face.")
    if seed is not None:
        set_random_seed(seed)

    from sage.misc.prandom import shuffle
    pendant_word = [0] * (k-1) + [1] * (k-3)
    shuffle(pendant_word)
    forest_word = _auxiliary_random_forest_word(n-k, k)
    word, graph = _contour_and_graph_from_words(pendant_word, forest_word)
    edges = []
    embedding = graph.get_embedding()

    pattern = ['in', 'in', 'in', 'lf', 'in']  # 'partial closures'

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
            in1, in2, in3 = (u[1] for u in word2[:3])
            edges.append([in1, in3])  # edge 'in1,in3'
            idx = embedding[in1].index(in2)
            embedding[in1].insert(idx, in3)
            idx = embedding[in3].index(in2)
            embedding[in3].insert(idx + 1, in1)
        else:
            break

    graph.add_edges(edges)
    graph.set_embedding(embedding)
    graph.relabel({0: -2, 1: -1})
    assert graph.num_edges() == 3*n - 3 - k
    assert graph.num_verts() == n
    if set_position:
        graph.layout(layout="planar", save_pos=True)
    return graph


def blossoming_contour(t, shift=0, seed=None):
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

        sage: t = BinaryTrees(2).random_element()                                       # needs sage.combinat
        sage: print(blossoming_contour(t))  # random                                    # needs sage.combinat
        [('i', 0), ('xb',), ('i', 0), ('n', 2), ('i', 1), ('xb',), ('i', 1),
        ('xb',), ('i', 1), ('n', 2), ('x',), ('n', 2), ('i', 0)]

        sage: w = blossoming_contour(BinaryTrees(3).random_element()); len(w)           # needs sage.combinat
        21
        sage: w.count(('xb',))                                                          # needs sage.combinat
        4
        sage: w.count(('x',))                                                           # needs sage.combinat
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
    if seed is not None:
        set_random_seed(seed)

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


def RandomBicubicPlanar(n, seed=None):
    """
    Return the graph of a random bipartite cubic map with `3 n` edges.

    INPUT:

    `n` -- an integer (at least `1`)

    - ``seed`` -- a ``random.Random`` seed or a Python ``int`` for the random
      number generator (default: ``None``)

    OUTPUT:

    a graph with multiple edges (no embedding is provided)

    The algorithm used is described in [Sch1999]_. This samples
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

        sage: # needs sage.combinat
        sage: n = randint(200, 300)
        sage: G = graphs.RandomBicubicPlanar(n)
        sage: G.order() == 2*n
        True
        sage: G.size() == 3*n
        True
        sage: G.is_bipartite() and G.is_planar() and G.is_regular(3)
        True
        sage: dic = {'red': [v for v in G.vertices(sort=False) if v[0] == 'n'],
        ....:        'blue': [v for v in G.vertices(sort=False) if v[0] != 'n']}
        sage: G.plot(vertex_labels=False, vertex_size=20, vertex_colors=dic)            # needs sage.plot
        Graphics object consisting of ... graphics primitives

    .. PLOT::
        :width: 300 px

        G = graphs.RandomBicubicPlanar(200)
        V0 = [v for v in G.vertices(sort=False) if v[0] == 'n']
        V1 = [v for v in G.vertices(sort=False) if v[0] != 'n']
        dic = {'red': V0, 'blue': V1}
        sphinx_plot(G.plot(vertex_labels=False,vertex_colors=dic))
    """
    from sage.combinat.binary_tree import BinaryTrees
    from sage.rings.finite_rings.integer_mod_ring import Zmod
    if not n:
        raise ValueError("n must be at least 1")
    if seed is not None:
        set_random_seed(seed)

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


def RandomUnitDiskGraph(n, radius=.1, side=1, seed=None):
    r"""
    Return a random unit disk graph of order `n`.

    A unit disk graph is the intersection graph of a family of unit disks in the
    Euclidean plane. That is a graph with one vertex per disk of the family and
    an edge between two vertices whenever they lie within a unit distance of
    each other. See the :wikipedia:`Unit_disk_graph` for more details.

    INPUT:

    - ``n`` -- number of nodes

    - ``radius`` -- float (default: ``0.1``); two vertices at distance less than
      ``radius`` are connected by an edge

    - ``side`` -- float (default: ``1``); indicate the side of the area in which
      the points are drawn

    - ``seed`` -- seed of the random number generator

    EXAMPLES:

    When using twice the same seed, the vertices get the same positions::

        sage: # needs scipy
        sage: from sage.misc.randstate import current_randstate
        sage: seed = current_randstate().seed()
        sage: G = graphs.RandomUnitDiskGraph(20, radius=.5, side=1, seed=seed)
        sage: H = graphs.RandomUnitDiskGraph(20, radius=.2, side=1, seed=seed)
        sage: H.is_subgraph(G, induced=False)
        True
        sage: H.size() <= G.size()
        True
        sage: Gpos = G.get_pos()
        sage: Hpos = H.get_pos()
        sage: all(Gpos[u] == Hpos[u] for u in G)
        True

    When the radius is more than `\sqrt{2 \text{side}}`, the graph is a clique::

        sage: G = graphs.RandomUnitDiskGraph(10, radius=2, side=1)                      # needs scipy
        sage: G.is_clique()                                                             # needs scipy
        True
    """
    if seed is not None:
        set_random_seed(seed)
    from scipy.spatial import KDTree
    points = [(side*random(), side*random()) for i in range(n)]
    T = KDTree(points)
    adj = {i: [u for u in T.query_ball_point([points[i]], radius).item() if u != i]
           for i in range(n)}
    return Graph(adj, format='dict_of_lists',
                 pos={i: points[i] for i in range(n)},
                 name="Random unit disk graph")
