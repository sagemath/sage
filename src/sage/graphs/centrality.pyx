# sage_setup: distribution = sagemath-graphs
r"""
Centrality

This module is meant for all functions related to centrality in networks.

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :func:`centrality_betweenness` | Return the centrality betweenness of `G`
    :func:`centrality_closeness_top_k` | Return the k most closeness central vertices of `G`
    :func:`centrality_closeness_random_k` | Return an estimation of the closeness centrality of `G`.

Functions
---------
"""

from libc.string cimport memset
from libc.stdint cimport uint32_t
from cysignals.memory cimport check_allocarray, sig_free
from cysignals.signals cimport sig_check
from memory_allocator cimport MemoryAllocator

from sage.data_structures.bitset_base cimport *
from sage.graphs.base.static_sparse_graph cimport *
from sage.libs.gmp.mpq cimport *
from sage.rings.rational cimport Rational
from sage.graphs.base.boost_graph import shortest_paths as boost_shortest_paths
import random

ctypedef fused numerical_type:
    mpq_t
    double

cimport cython


def centrality_betweenness(G, bint exact=False, bint normalize=True):
    r"""
    Return the centrality betweenness of `G`.

    The centrality betweenness of a vertex `v\in G` is defined by:

    .. MATH::

        c(v) = \sum_{s\neq v \neq t} \frac{\#\{\text{shortest } st-\text{paths containing v}\}}
                                          {\#\{\text{shortest } st-\text{paths}\}}

    For more information, see the :wikipedia:`Betweenness_centrality`.

    INPUT:

    - ``G`` -- a (di)graph

    - ``exact`` -- boolean (default: ``False``); whether to compute over
      rationals or on ``double`` C variables

    - ``normalize`` -- boolean (default: ``True``); whether to renormalize the
      values by dividing them by `\binom {n-1} 2` (for graphs) or `2\binom {n-1}
      2` (for digraphs).

    ALGORITHM:

    To compute `c(v)`, we fix `s` and define `c_s(v)` as the centrality of `v`
    *due to* `s`, obtained from the formula above by running the sum over `t`
    only. We obtain `c(v)=\sum_{s\neq v} c_s(v)`.

    For every vertex `s`, we compute the value of `c_s(v)` for all `v`, using
    the following remark (see [Brandes01]_):

        Let `v_1,\ldots,v_k` be the out-neighbors of `v` such that
        `dist(s,v_i) = dist(s,v) + 1`. Then

        .. MATH::

            c_s(v) = \sum_{1\leq i \leq k} c_s(v_i)
                         \frac{\#\{\text{shortest } sv_i-\text{paths}\}}
                              {\#\{\text{shortest } sv  -\text{paths}\}}

    The number of shortest paths between `s` and every other vertex can be
    computed with a slightly modified BFS. While running this BFS we can also
    store the list of the vertices `v_1,\ldots,v_k` associated with each `v`.

    EXAMPLES::

        sage: from sage.graphs.centrality import centrality_betweenness
        sage: centrality_betweenness(digraphs.Circuit(6)) # abs tol 1e-10
        {0: 0.5, 1: 0.5, 2: 0.5, 3: 0.5, 4: 0.5, 5: 0.5}
        sage: centrality_betweenness(graphs.CycleGraph(6)) # abs tol 1e-10
        {0: 0.2, 1: 0.2, 2: 0.2, 3: 0.2, 4: 0.2, 5: 0.2}

    Exact computations::

        sage: graphs.PetersenGraph().centrality_betweenness(exact=True)
        {0: 1/12, 1: 1/12, 2: 1/12, 3: 1/12, 4: 1/12, 5: 1/12, 6: 1/12, 7: 1/12, 8: 1/12, 9: 1/12}

    TESTS:

    Compare with NetworkX::

        sage: # needs networkx
        sage: import networkx
        sage: g = graphs.RandomGNP(100, .2)
        sage: nw = networkx.betweenness_centrality(g.networkx_graph())
        sage: sg = centrality_betweenness(g)
        sage: max(abs(nw[x] - sg[x]) for x in g)  # abs tol 1e-10
        0

    Stupid cases::

        sage: centrality_betweenness(Graph())
        {}
        sage: centrality_betweenness(Graph(2))
        {0: 0.0, 1: 0.0}
        sage: centrality_betweenness(Graph(2), exact=1)
        {0: 0, 1: 0}
    """
    if exact:
        return centrality_betweenness_C(G, <mpq_t> 0, normalize=normalize)
    return centrality_betweenness_C(G, <double>0, normalize=normalize)


@cython.cdivision(True)
cdef dict centrality_betweenness_C(G, numerical_type _, bint normalize=True):
    r"""
    Return the centrality betweenness of G (C implementation).

    INPUT:

    - ``G`` -- a graph

    - ``_`` -- this variable is ignored, only its type matters. If it is of type
      `mpq_t` then computations are made on `Q`, if it is ``double`` the
      computations are made on ``double``.

    - ``normalize`` -- boolean (default: ``True``); whether to renormalize the
      values by dividing them by `\binom {n-1} 2` (for graphs) or `2\binom {n-1}
      2` (for digraphs).

    For more information, see the documentation of ``centrality_betweenness``.
    """
    # Trivial case
    if G.order() <= 2:
        zero = 0. if numerical_type is double else Rational(0)
        return {v: zero for v in G}

    # A copy of G, for faster neighbor enumeration
    cdef short_digraph g

    # A second copy, to remember the edges used during the BFS (see doc)
    cdef short_digraph bfs_dag

    cdef list int_to_vertex = list(G)

    cdef int n = G.order()

    cdef bitset_t seen  # Vertices whose neighbors have been explored
    cdef bitset_t next_layer  # Unexplored neighbors of vertices in 'seen'

    cdef uint32_t* queue = NULL  # BFS queue
    cdef uint32_t* degrees = NULL  # degree[v] = nb of vertices which discovered v

    cdef numerical_type* n_paths_from_source = NULL
    cdef numerical_type* betweenness_source = NULL
    cdef numerical_type* betweenness = NULL
    cdef numerical_type mpq_tmp

    cdef int layer_current_beginning
    cdef int layer_current_end
    cdef int layer_next_end

    cdef int source, i, j, u, v
    cdef uint32_t* p_tmp
    cdef list betweenness_list

    if numerical_type is mpq_t:
        mpq_init(mpq_tmp)

    try:
        init_short_digraph(g, G, edge_labelled=False, vertex_list=int_to_vertex)
        init_reverse(bfs_dag, g)

        queue = <uint32_t*> check_allocarray(n, sizeof(uint32_t))
        degrees = <uint32_t*> check_allocarray(n, sizeof(uint32_t))
        n_paths_from_source = <numerical_type*> check_allocarray(n, sizeof(numerical_type))
        betweenness_source = <numerical_type*> check_allocarray(n, sizeof(numerical_type))
        betweenness = <numerical_type*> check_allocarray(n, sizeof(numerical_type))

        bitset_init(seen, n)
        bitset_init(next_layer, n)

        if numerical_type is double:
            memset(betweenness, 0, n * sizeof(double))
        else:
            for i in range(n):
                mpq_init(betweenness[i])
                mpq_set_ui(betweenness[i], 0, 1)
                mpq_init(betweenness_source[i])
                mpq_init(n_paths_from_source[i])

        for source in range(n):

            if numerical_type is double:
                memset(betweenness_source, 0, n * sizeof(double))
                memset(n_paths_from_source, 0, n * sizeof(double))
                n_paths_from_source[source] = 1
            else:
                for i in range(n):
                    mpq_set_ui(betweenness_source[i], 0, 1)
                    mpq_set_ui(n_paths_from_source[i], 0, 1)
                mpq_set_ui(n_paths_from_source[source], 1, 1)

            # initialize data
            bitset_set_first_n(seen, 0)
            bitset_add(seen, source)
            bitset_set_first_n(next_layer, 0)

            memset(degrees, 0, n * sizeof(uint32_t))

            queue[0] = source
            layer_current_beginning = 0
            layer_current_end = 1
            layer_next_end = 1

            # The number of shortest paths from 'source' to every other vertex.
            #
            # It is a BFS. The graph is explored layer by layer.
            while layer_current_beginning<layer_current_end:

                # Looking for all non-discovered neighbors of some vertex of the
                # current layer.
                for j in range(layer_current_beginning, layer_current_end):
                    u = queue[j]

                    # List the neighbors of u
                    p_tmp = g.neighbors[u]
                    while p_tmp < g.neighbors[u + 1]:
                        v = p_tmp[0]
                        p_tmp += 1

                        # Is it a new vertex ?
                        if bitset_in(seen, v):
                            continue

                        # Is it the first time we see it ?
                        elif not bitset_in(next_layer, v):
                            bitset_add(next_layer, v)
                            queue[layer_next_end] = v
                            layer_next_end += 1

                        # update the count of paths and the BFS dag.
                        bfs_dag.neighbors[v][degrees[v]] = u
                        degrees[v] += 1
                        if numerical_type is double:
                            n_paths_from_source[v] += n_paths_from_source[u]
                        else:
                            mpq_add(n_paths_from_source[v], n_paths_from_source[v], n_paths_from_source[u])

                # 'next_layer' becomes 'current_layer'
                for j in range(layer_current_end, layer_next_end):
                    bitset_add(seen, queue[j])

                layer_current_beginning = layer_current_end
                layer_current_end = layer_next_end

            # Compute the betweenness from the number of paths
            #
            # We enumerate vertices in reverse order of discovery.
            for i in range(layer_current_end - 1, -1, -1):
                u = queue[i]
                for j in range(<int>degrees[u]):
                    v = bfs_dag.neighbors[u][j]
                    if v != source:  # better to not 'if' but set it to 0 afterwards?
                        if numerical_type is double:
                            betweenness_source[v] += (betweenness_source[u] + 1) * (n_paths_from_source[v] / n_paths_from_source[u])
                        else:
                            mpq_set_ui(mpq_tmp, 1, 1)
                            mpq_add(mpq_tmp, betweenness_source[u], mpq_tmp)
                            mpq_mul(mpq_tmp, mpq_tmp, n_paths_from_source[v])
                            mpq_div(mpq_tmp, mpq_tmp, n_paths_from_source[u])
                            mpq_add(betweenness_source[v], betweenness_source[v], mpq_tmp)

            # update betweenness from betweenness_source
            for i in range(n):
                if numerical_type is double:
                    betweenness[i] += betweenness_source[i]
                else:
                    mpq_add(betweenness[i], betweenness[i], betweenness_source[i])

            sig_check()  # check for KeyboardInterrupt

        if numerical_type is double:
            betweenness_list = [betweenness[i] for i in range(n)]
        else:
            betweenness_list = [Rational(None) for x in range(n)]

            for i in range(n):
                (<Rational> (betweenness_list[i])).set_from_mpq(betweenness[i])
            for i in range(n):
                mpq_clear(betweenness_source[i])
                mpq_clear(betweenness[i])
                mpq_clear(n_paths_from_source[i])
            mpq_clear(mpq_tmp)

    finally:
        free_short_digraph(g)
        free_short_digraph(bfs_dag)
        bitset_free(seen)
        bitset_free(next_layer)
        sig_free(queue)
        sig_free(n_paths_from_source)
        sig_free(degrees)
        sig_free(betweenness_source)
        sig_free(betweenness)

    if not G.is_directed():
        betweenness_list = [x / 2 for x in betweenness_list]

    if normalize:
        if G.is_directed():
            betweenness_list = [x / ((n - 1) * (n - 2)) for x in betweenness_list]
        else:
            betweenness_list = [2 * x / ((n - 1) * (n - 2)) for x in betweenness_list]

    return {vv: betweenness_list[i] for i, vv in enumerate(int_to_vertex)}


cdef void _estimate_reachable_vertices_dir(short_digraph g, int* reachL, int* reachU) noexcept:
    r"""
    For each vertex ``v``, bounds the number of vertices reachable from ``v``.

    The lower bound is stored in ``reachL[v]``, while the upper bound is stored
    in ``reachU[v]``. These two arrays must be pre-allocated and they must have
    size at least ``n``, where ``n`` is the number of nodes of ``g``.

    The estimate works as follows: first, we compute the graph of strongly
    connected components `\mathcal{G=(V,E)}`, then, for each SCC `C`, we set:

    .. MATH::

        L(C)=|C|+\max_{(C,C') \in \mathcal{E}}L(C') \\
        U(C)=|C|+\max_{(C,C') \in \mathcal{E}}L(C')

    By analyzing strongly connected components in reverse topological order, we
    are sure that, as soon as we process component `C`, all components `C'`
    appearing on the right hand side have already been processed.  A further
    improvement on these bounds is obtained by exactly computing the number of
    vertices reachable from the biggest strongly connected component, and handle
    this component separately.

    Then, for each vertex ``v``, we set ``reachL[v]=L(C)``, where `C` is the
    strongly connected component containing ``v``.

    INPUT:

    - ``g`` -- short_digraph; the input graph;

    OUTPUT:

    ``reachL``, ``reachU``: two arrays that should be allocated outside this
    function and that should have size at least ``g.n``. At the end,
    ``reachL[v]`` (resp., ``reachU[v]``) will contain the lower (resp., upper)
    bound on the number of reachable vertices from ``v``.
    """
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef int n = g.n
    cdef int* scc = <int*> mem.malloc(n * sizeof(int))
    cdef int i, v, w, maxscc = 0
    cdef int nscc = tarjan_strongly_connected_components_C(g, scc)
    cdef short_digraph sccgraph
    strongly_connected_components_digraph_C(g, nscc, scc, sccgraph)

    cdef int* scc_sizes = <int*> mem.calloc(nscc, sizeof(int))
    cdef int nreach_maxscc = 0
    cdef short* reach_max_scc = <short*> mem.calloc(nscc, sizeof(short))
    cdef int* reachL_scc = <int*> mem.calloc(nscc, sizeof(int))
    cdef uint64_t* reachU_scc = <uint64_t*> mem.calloc(nscc, sizeof(uint64_t))
    cdef uint64_t* reachU_without_maxscc = <uint64_t*> mem.calloc(nscc, sizeof(uint64_t))
    # We need uint64_t because these values may become much bigger than g.n,
    # up to g.n^2, during the computation. Only at the end, we set reachL and
    # reachU as the maximum between g.n and the computed value (so that they
    # can be converted to int without overflow).

    # Variables used in BFS from the largest strongly connected component
    cdef uint32_t startq, endq
    cdef int* q = <int*> mem.malloc(nscc * sizeof(int))
    cdef short* reached = <short*> mem.calloc(nscc, sizeof(short))
    cdef uint32_t* neigh_start
    cdef uint32_t* neigh_end

    # Compute scc_sizes
    for i in range(g.n):
        scc_sizes[scc[i]] += 1

    # Compute maxscc
    for i in range(nscc):
        if scc_sizes[maxscc] < scc_sizes[i]:
            maxscc = i
    reach_max_scc[maxscc] = 1

    # BFS to compute number of reachable vertices for the biggest SCC.
    q[0] = maxscc
    nreach_maxscc = scc_sizes[maxscc]
    reached[maxscc] = 1
    startq = 0
    endq = 1
    while startq < endq:
        v = q[startq]
        startq += 1
        neigh_start = sccgraph.neighbors[v]
        neigh_end = sccgraph.neighbors[v + 1]

        while neigh_start < neigh_end:
            w = neigh_start[0]
            if not reached[w]:
                reached[w] = 1
                nreach_maxscc += scc_sizes[w]
                q[endq] = w
                endq += 1
            neigh_start += 1

    reachL_scc[maxscc] = nreach_maxscc
    reachU_scc[maxscc] = nreach_maxscc
    reachU_without_maxscc[maxscc] = 0
    # Dynamic programming to estimate number of reachable vertices for other
    # SCCs
    for i in range(nscc):
        if i == maxscc:
            continue

        neigh_start = sccgraph.neighbors[i]
        neigh_end = sccgraph.neighbors[i + 1]

        while neigh_start < neigh_end:
            w = neigh_start[0]
            neigh_start += 1

            reachL_scc[i] = max(reachL_scc[i], reachL_scc[w])
            reachU_scc[i] += reachU_scc[w]
            # Note that this might become much bigger than g.n, up to g.n*g.n.
            # Hence we used uint64_t, and only at the end we take the minimum
            # between this value and g.n (since g.n is an upper bound on
            # the number of reachable vertices).
            if not reached[w]:
                reachU_without_maxscc[i] += reachU_without_maxscc[w]
            reach_max_scc[i] = reach_max_scc[i] or reach_max_scc[w]

        if reach_max_scc[i]:
            reachU_scc[i] = reachU_without_maxscc[i] + nreach_maxscc

        reachL_scc[i] += scc_sizes[i]
        reachU_scc[i] += scc_sizes[i]
        if not reached[i]:
            reachU_without_maxscc[i] += scc_sizes[i]

    for i in range(n):
        reachL[i] = reachL_scc[scc[i]]
        reachU[i] = min(<int>reachU_scc[scc[i]], g.n)


cdef void _compute_reachable_vertices_undir(short_digraph g, int* reachable) noexcept:
    r"""
    For each vertex ``v``, compute the number of vertices reachable from ``v``.

    The number of vertices reachable from ``v`` (which is the size of the
    connected component containing ``v``) is stored in variable
    ``reachable[v]``. The array `reachable` is assumed to be allocated outside
    this function, and it is assumed to have size at least ``g.n``.
    """
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef int i
    cdef int n = g.n
    cdef int* q = <int*> mem.malloc(n * sizeof(int))
    cdef short* reached = <short*> mem.calloc(n, sizeof(short))

    cdef int v, w
    cdef uint32_t* neigh_start
    cdef uint32_t* neigh_end
    cdef uint32_t startq, endq
    cdef list currentcc

    memset(reachable, 0, n * sizeof(int))

    for i in range(n):
        # BFS from i
        if reachable[i]:
            continue

        reached[i] = 1
        currentcc = [i]

        q[0] = i
        startq = 0
        endq = 1
        while startq < endq:
            v = q[startq]
            startq += 1
            neigh_start = g.neighbors[v]
            neigh_end = g.neighbors[v + 1]

            while neigh_start < neigh_end:
                w = neigh_start[0]
                if not reached[w]:
                    reached[w] = 1
                    currentcc.append(w)
                    q[endq] = w
                    endq += 1
                neigh_start += 1

        for v in currentcc:
            reachable[v] = len(currentcc)


cdef void _sort_vertices_degree(short_digraph g, int* sorted_verts) noexcept:
    r"""
    Sort vertices in decreasing order of degree.

    Uses counting sort, since degrees are between `0` and `n-1`: the running
    time is then `O(n)`.
    """
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef uint32_t* verts_of_degree = <uint32_t*> mem.calloc(g.n, sizeof(uint32_t))
    cdef uint32_t* next_vert_of_degree = <uint32_t*> mem.malloc(g.n * sizeof(uint32_t))
    cdef int d, v

    # Otherwise, segmentation fault
    if not g.n:
        return

    for v in range(g.n):
        verts_of_degree[out_degree(g, v)] += 1

    next_vert_of_degree[g.n - 1] = 0
    for i in range(g.n - 2, -1, -1):
        next_vert_of_degree[i] = next_vert_of_degree[i + 1] + verts_of_degree[i + 1]

    for v in range(g.n):
        d = out_degree(g, v)
        sorted_verts[next_vert_of_degree[d]] = v
        next_vert_of_degree[d] += 1


def centrality_closeness_top_k(G, int k=1, int verbose=0):
    r"""
    Compute the ``k`` vertices with largest closeness centrality.

    The algorithm is based on performing a breadth-first-search (BFS) from each
    vertex, and to use bounds in order to cut these BFSes as soon as possible.
    If ``k`` is small, it is much more efficient than computing all centralities
    with :meth:`~sage.graphs.generic_graph.GenericGraph.centrality_closeness`.
    Conversely, if ``k`` is close to the number of nodes, the running-time is
    approximately the same (it might even be a bit longer, because more
    computations are needed).

    For more information, see [BCM15]_. The algorithm does not work on
    weighted graphs.

    INPUT:

    - ``G`` -- a Sage Graph or DiGraph;

    - ``k`` -- integer (default: 1); the algorithm will return the ``k``
      vertices with largest closeness centrality. This value should be between 1
      and the number of vertices with positive (out)degree, because the
      closeness centrality is not defined for vertices with (out)degree 0. If
      ``k`` is bigger than this value, the output will contain all vertices of
      positive (out)degree.

    - ``verbose`` -- integer (default: 0); define how "verbose" the algorithm
      should be. If 0, nothing is printed, if 1, we print only the performance
      ratio at the end of the algorithm, if 2, we print partial results every
      1000 visits, if 3, we print partial results after every visit.

    OUTPUT:

    An ordered list of ``k`` pairs `(closv, v)`, where `v` is one of the ``k``
    most central vertices, and `closv` is its closeness centrality.  If `k` is
    bigger than the number of vertices with positive (out)degree, the list might
    be smaller.

    EXAMPLES::

        sage: from sage.graphs.centrality import centrality_closeness_top_k
        sage: g = graphs.PathGraph(10)
        sage: centrality_closeness_top_k(g, 4, 1)
        Final performance ratio: 0.711111111111...
        [(0.36, 5),
         (0.36, 4),
         (0.3333333333333333, 6),
         (0.3333333333333333, 3)]
        sage: g = digraphs.Path(10)
        sage: centrality_closeness_top_k(g, 5, 1)
        Final performance ratio: 0.422222222222...
        [(0.2, 0),
         (0.19753086419753085, 1),
         (0.19444444444444442, 2),
         (0.19047619047619047, 3),
         (0.18518518518518517, 4)]

    TESTS:

    If ``k`` or ``verbose`` is not an integer::

        sage: from sage.graphs.centrality import centrality_closeness_top_k
        sage: g = digraphs.Path(10)
        sage: centrality_closeness_top_k(g, 'abc', 1)
        Traceback (most recent call last):
        ...
        TypeError: an integer is required
        sage: centrality_closeness_top_k(g, 1, 'abc')
        Traceback (most recent call last):
        ...
        TypeError: an integer is required

    If ``k`` is bigger than the number of nodes::

        sage: from sage.graphs.centrality import centrality_closeness_top_k
        sage: g = graphs.PathGraph(5)
        sage: centrality_closeness_top_k(g, 10, 0)
        [(0.6666666666666666, 2),
         (0.5714285714285714, 3),
         (0.5714285714285714, 1),
         (0.4, 4),
         (0.4, 0)]

    Empty graph::

        sage: from sage.graphs.centrality import centrality_closeness_top_k
        sage: g = Graph()
        sage: centrality_closeness_top_k(g, 10, 0)
        []
        sage: g = Graph(10)
        sage: centrality_closeness_top_k(g, 10, 0)
        []

    The result is correct::

        sage: # needs networkx
        sage: from sage.graphs.centrality import centrality_closeness_top_k
        sage: import random
        sage: n = 20
        sage: m = random.randint(1, n * (n - 1) / 2)
        sage: k = random.randint(1, n)
        sage: g = graphs.RandomGNM(n, m)
        sage: topk = centrality_closeness_top_k(g, k)
        sage: centr = g.centrality_closeness(algorithm='BFS')
        sage: sorted_centr = sorted(centr.values(), reverse=True)
        sage: len(topk) == min(k, len(sorted_centr))
        True
        sage: all(abs(topk[i][0] - sorted_centr[i]) < 1e-12 for i in range(len(topk)))
        True

    Directed case::

        sage: # needs networkx
        sage: from sage.graphs.centrality import centrality_closeness_top_k
        sage: import random
        sage: n = 20
        sage: m = random.randint(1, n * (n - 1))
        sage: k = random.randint(1, n)
        sage: g = digraphs.RandomDirectedGNM(n, m)
        sage: topk = centrality_closeness_top_k(g, k)
        sage: centr = g.centrality_closeness(algorithm='BFS')
        sage: sorted_centr = sorted(centr.values(), reverse=True)
        sage: len(topk) == min(k, len(sorted_centr))
        True
        sage: all(abs(topk[i][0] - sorted_centr[i]) < 1e-12 for i in range(len(topk)))
        True
    """
    cdef list res
    if k >= G.order():
        closeness_dict = G.centrality_closeness(by_weight=False, algorithm='BFS')
        res = [(closz, z) for z, closz in closeness_dict.items()]
        try:
            res = sorted(res, reverse=True)
        except TypeError:
            # A TypeError may occur in Python 3 when vertex labels are of
            # different types. We then sort on values only.
            res = sorted(res, reverse=True, key=lambda zz: zz[0])
        return res

    if G.order() < 2:
        return []

    cdef MemoryAllocator mem = MemoryAllocator()
    cdef short_digraph sd
    # Copying the whole graph to obtain the list of neighbors quicker than by
    # calling out_neighbors. This data structure is well documented in the
    # module sage.graphs.base.static_sparse_graph
    cdef list V = list(G)
    init_short_digraph(sd, G, edge_labelled=False, vertex_list=V)
    cdef int n = sd.n
    cdef int* reachL = <int*> mem.malloc(n * sizeof(int))
    cdef int* reachU
    cdef int* pred = <int*> mem.calloc(n, sizeof(int))
    cdef double *farness = <double*> mem.malloc(n * sizeof(double))
    cdef int d, nd, x, v
    cdef long f, gamma
    cdef int* queue = <int*> mem.malloc(n * sizeof(int))
    cdef double tildefL, tildefU
    cdef bint stopped
    cdef uint32_t* p_tmp
    cdef int layer_current_beginning, layer_current_end, layer_next_end=0
    cdef long visited = 0
    cdef int nvis = 0
    cdef short* seen = <short*> mem.calloc(n, sizeof(short))
    cdef bint directed = G.is_directed()

    cdef int* topk = <int*> mem.malloc(k * sizeof(int))
    for i in range(k):
        topk[i] = -1
    for i in range(n):
        pred[i] = -1

    cdef double kth = n
    cdef int* sorted_vert = <int*> mem.malloc(n * sizeof(int))
    if directed:
        reachU = <int*> mem.malloc(n * sizeof(int))
        _estimate_reachable_vertices_dir(sd, reachL, reachU)
    else:
        _compute_reachable_vertices_undir(sd, reachL)
        reachU = reachL
    _sort_vertices_degree(sd, sorted_vert)

    for x in sorted_vert[:n]:
        sig_check()

        if not out_degree(sd, x):
            break
        # We start a BFSCut from x:

        # We reset variable seen:
        for v in queue[:layer_next_end]:
            seen[v] = 0
            pred[v] = -1

        layer_current_beginning = 0
        layer_current_end = 1
        layer_next_end = 1
        d = 0
        f = 0
        # We are at level 0, and gamma is the number of arcs exiting level 0
        # (hence, deg(x)).
        gamma = out_degree(sd, x)
        nd = 1
        queue[0] = x
        stopped = False
        seen[x] = 1
        nvis += 1

        # The graph is explored layer by layer.
        while layer_current_beginning < layer_current_end and not stopped:
            sig_check()

            # We update our estimate of the farness of v.
            # The estimate sets distance d+1 to gamma vertices (which is an
            # upper bound on the number of vertices at distance d+1 from v),
            # and distance d+2 to all other vertices reachable from x.
            tildefL = ((f - gamma + (d + 2) * (<double>(reachL[x] - nd))) * (n - 1)) / ((<double>(reachL[x] - 1)) * (reachL[x] - 1))
            tildefU = ((f - gamma + (d + 2) * (<double>(reachU[x] - nd))) * (n - 1)) / ((<double>(reachU[x] - 1)) * (reachU[x] - 1))
            d += 1
            gamma = 0

            if tildefL >= kth and tildefU >= kth:
                farness[x] = n
                stopped = True
                break
            # Looking for all non-discovered neighbors of some vertex of the
            # current layer.
            for j in range(layer_current_beginning, layer_current_end):
                sig_check()

                u = queue[j]

                # List the neighbors of u
                p_tmp = sd.neighbors[u]
                while p_tmp < sd.neighbors[u + 1] and not stopped:
                    sig_check()

                    visited += 1
                    v = p_tmp[0]
                    p_tmp += 1
                    # Is it a new vertex ?
                    if not seen[v]:
                        seen[v] = 1
                        queue[layer_next_end] = v
                        layer_next_end += 1
                        f = f + d
                        gamma += out_degree(sd, v) if directed else (out_degree(sd, v) - 1)
                        nd = nd + 1
                        pred[v] = u
                    elif directed or pred[u] != v:
                        tildefL += (n - 1) / (<double>(reachL[x] - 1) * (reachL[x] - 1))
                        tildefU += (n - 1) / (<double>(reachU[x] - 1) * (reachU[x] - 1))
                        if tildefL >= kth and tildefU >= kth:
                            farness[x] = n
                            stopped = True
                if stopped:
                    break
            # 'next_layer' becomes 'current_layer'
            layer_current_beginning = layer_current_end
            layer_current_end = layer_next_end

        if not stopped:
            farness[x] = ((<double> f) * (n - 1)) / (<double>(nd - 1) * (nd - 1))

        if farness[x] < kth:
            for i in range(k):
                if topk[i] == -1 or farness[topk[i]] == kth:
                    topk[i] = x
                    break
            kth = 0
            for i in range(k):
                if topk[i] == -1:
                    kth = n
                    break
                kth = max(kth, farness[topk[i]])
        if verbose >= 3 or (verbose == 2 and nvis % 1000 == 0):
            print("Visit {} from {}:".format(nvis, x))
            print("    Lower bound: {}".format(1 / kth))
            print("    Perf. ratio: {}".format(visited / (nvis * <double> (sd.neighbors[sd.n] - sd.edges))))

    if verbose > 0:
        print("Final performance ratio: {}".format(visited / (n * <double> (sd.neighbors[sd.n] - sd.edges))))

    res = [(1.0 / farness[v], V[v]) for v in topk[:k] if v != -1]
    try:
        res = sorted(res, reverse=True)
    except TypeError:
        # A TypeError may occur in Python 3 when vertex labels are of
        # different types. We then sort on values only.
        res = sorted(res, reverse=True, key=lambda vv: vv[0])
    return res


def centrality_closeness_random_k(G, int k=1):
    r"""
    Return an estimation of the closeness centrality of `G`.

    The algorithm first randomly selects a set `S` of `k` vertices.  Then it
    computes shortest path distances from each vertex in `S` (using Dijkstra for
    weighted graph and breadth-first-search (BFS) for unweighted graph) and uses
    this knowledge to estimate the closeness centrality of all vertices.

    For more information, see [EDI2014]_.

    INPUT:

    - ``G`` -- an undirected connected Graph

    - ``k`` -- integer (default: 1); number of random nodes to choose

    OUTPUT: a dictionary associating to each vertex its estimated closeness centrality

    EXAMPLES:

    Estimation of the closeness centrality of the Petersen Graph when `k == n`::

        sage: from sage.graphs.centrality import centrality_closeness_random_k
        sage: G = graphs.PetersenGraph()
        sage: centrality_closeness_random_k(G, 10)
        {0: 0.6,
         1: 0.6,
         2: 0.6,
         3: 0.6,
         4: 0.6,
         5: 0.6,
         6: 0.6,
         7: 0.6,
         8: 0.6,
         9: 0.6}

    TESTS::

        sage: from sage.graphs.centrality import centrality_closeness_random_k
        sage: G = Graph('Ihe\\n@GUA')
        sage: centrality_closeness_random_k(G, -1)
        Traceback (most recent call last):
        ...
        ValueError: parameter k must be a positive integer

        sage: centrality_closeness_random_k(DiGraph(2), 1)
        Traceback (most recent call last):
        ...
        ValueError: G must be an undirected Graph
    """
    G._scream_if_not_simple()
    if G.is_directed():
        raise ValueError("G must be an undirected Graph")

    cdef int n = G.order()
    if not n:
        return {}
    if k < 1:
        raise ValueError("parameter k must be a positive integer")
    if k > n:
        k = n

    if not G.is_connected():
        raise ValueError("G must be a connected Graph")

    # Initialization of some data structures
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef double* partial_farness = <double*> mem.malloc(n * sizeof(double))
    cdef uint32_t* distance
    cdef uint32_t* waiting_list
    cdef short_digraph sd
    cdef bitset_t seen
    cdef double farness
    cdef int i, j
    cdef dict closeness_centrality_array = {}
    cdef list int_to_vertex = list(G)
    cdef dict vertex_to_int = {v: i for i, v in enumerate(int_to_vertex)}

    # Initialize
    for i in range(n):
        partial_farness[i] = 0

    # Shuffle the vertices
    cdef list V = list(range(n))
    random.shuffle(V)

    if G.weighted():
        # For all random nodes take as a source then run Dijstra and
        # calculate closeness centrality for k random vertices from l.
        for i in range(k):
            farness = 0
            distances = boost_shortest_paths(G, int_to_vertex[V[i]], algorithm='Dijkstra')[0]
            for vertex in distances:
                farness += float(distances[vertex])
                partial_farness[vertex_to_int[vertex]] += float(distances[vertex])

            closeness_centrality_array[int_to_vertex[V[i]]] = (n - 1) / farness

    # G is unweighted graph
    else:

        # Copying the whole graph as a static_sparse_graph for fast shortest
        # paths computation in unweighted graph. This data structure is well
        # documented in module sage.graphs.base.static_sparse_graph
        init_short_digraph(sd, G, edge_labelled=False, vertex_list=int_to_vertex)
        distance = <uint32_t*> mem.malloc(n * sizeof(uint32_t))
        waiting_list = <uint32_t*> mem.malloc(n * sizeof(uint32_t))
        bitset_init(seen, n)

        # Run BFS for random k vertices
        for i in range(k):
            farness = 0
            simple_BFS(sd, V[i], distance, NULL, waiting_list, seen)
            for j in range(n):
                farness += distance[j]
                partial_farness[j] += distance[j]

            closeness_centrality_array[int_to_vertex[V[i]]] = (n - 1) / farness

        bitset_free(seen)
        free_short_digraph(sd)

    # Estimate the closeness centrality for remaining n-k vertices.
    for i in range(k, n):
        closeness_centrality_array[int_to_vertex[V[i]]] = k / partial_farness[V[i]]

    return closeness_centrality_array
