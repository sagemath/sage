# distutils: libraries = gmp
r"""
Centrality

This module is meant for all functions related to centrality in networks.

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :func:`centrality_betweenness` | Return the centrality betweenness of `G`

Functions
---------
"""
include "sage/data_structures/bitset.pxi"
include "sage/ext/interrupt.pxi"

from sage.graphs.base.static_sparse_graph cimport *
from libc.stdint cimport uint32_t
from sage.libs.gmp.mpq cimport *
from sage.rings.rational cimport Rational
from sage.ext.memory cimport check_malloc, check_calloc

ctypedef fused numerical_type:
    mpq_t
    double

import cython

def centrality_betweenness(G, exact=False, normalize=True):
    r"""
    Return the centrality betweenness of `G`

    The centrality betweenness of a vertex `v\in G` is defined by:

    .. MATH::

        c(v) = \sum_{s\neq v \neq t} \frac{\#\{\text{shortest } st-\text{paths containing v}\}}
                                          {\#\{\text{shortest } st-\text{paths}\}}

    For more information, see the :wikipedia:`Betweenness_centrality`.

    INPUT:

    - ``G`` -- a (di)graph

    - ``exact`` (boolean, default: ``False``) -- whether to compute over
      rationals or on ``double`` C variables.

    - ``normalize`` (boolean; default: ``True``) -- whether to renormalize the
      values by dividing them by `\binom {n-1} 2` (for graphs) or `2\binom {n-1}
      2` (for digraphs).

    ALGORITHM:

    To compute `c(v)`, we fix `s` and define `c_s(v)` as the centrality of `v`
    *due to* `s`, obtained from the formula above by running the sum over `t`
    only. We obtain `c(v)=\sum_{s\neq v} c_s(v)`.

    For every vertex `s`, we compute the value of `c_s(v)` for all `v`, using
    the following remark (see [Brandes01]_):

        Let `v_1,...,v_k` be the out-neighbors of `v` such that
        `dist(s,v_i)=dist(s,v)+1`. Then

        .. MATH::

            c_s(v) = \sum_{1\leq i \leq k} c_s(v_i)
                         \frac{\#\{\text{shortest } sv_i-\text{paths}\}}
                              {\#\{\text{shortest } sv  -\text{paths}\}}

    The number of shortest paths between `s` and every other vertex can be
    computed with a slightly modified BFS. While running this BFS we can also
    store the list of the vertices `v_1,...,v_k` associated with each `v`.

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

        sage: import networkx
        sage: g = graphs.RandomGNP(100,.2)
        sage: nw = networkx.betweenness_centrality(g.networkx_graph(copy=False))
        sage: sg = centrality_betweenness(g)
        sage: max(abs(nw[x]-sg[x]) for x in g) # abs tol 1e-10
        0

    Stupid cases::

        sage: centrality_betweenness(Graph())
        {}
        sage: centrality_betweenness(Graph(2))
        {0: 0.0, 1: 0.0}
        sage: centrality_betweenness(Graph(2),exact=1)
        {0: 0, 1: 0}

    REFERENCES:

    .. [Brandes01] Ulrik Brandes,
       A faster algorithm for betweenness centrality,
       Journal of Mathematical Sociology 25.2 (2001): 163-177,
       http://www.inf.uni-konstanz.de/algo/publications/b-fabc-01.pdf
    """
    if exact:
        return centrality_betweenness_C(G,<mpq_t> 0,normalize=normalize)
    else:
        return centrality_betweenness_C(G,<double>0,normalize=normalize)

@cython.cdivision(True)
cdef dict centrality_betweenness_C(G, numerical_type _, normalize=True):
    r"""
    Return the centrality betweenness of G (C implementation)

    INPUT:

    - ``G`` -- a graph

    - ``_`` -- this variable is ignored, only its type matters. If it is of type
      `mpq_t` then computations are made on `Q`, if it is ``double`` the
      computations are made on ``double``.

    - ``normalize`` (boolean; default: ``True``) -- whether to renormalize the
      values by dividing them by `\binom {n-1} 2` (for graphs) or `2\binom {n-1}
      2` (for digraphs).

    For more information, see the documentation of ``centrality_betweenness``.
    """
    # Trivial case
    if G.order() <= 2:
        zero = 0. if numerical_type is double else Rational(0)
        return {v:zero for v in G}

    # A copy of G, for faster neighbor enumeration
    cdef short_digraph g

    # A second copy, to remember the edges used during the BFS (see doc)
    cdef short_digraph bfs_dag

    cdef int n = G.order()

    cdef bitset_t seen # Vertices whose neighbors have been explored
    cdef bitset_t next_layer # Unexplored neighbors of vertices in 'seen'

    cdef uint32_t * queue # BFS queue
    cdef uint32_t * degrees # degree[v] = nb of vertices which discovered v

    cdef numerical_type * n_paths_from_source
    cdef numerical_type * betweenness_source
    cdef numerical_type * betweenness
    cdef numerical_type mpq_tmp

    cdef int layer_current_beginning
    cdef int layer_current_end
    cdef int layer_next_end

    cdef int source,i,j,u,v
    cdef uint32_t * p_tmp

    if numerical_type is mpq_t:
        mpq_init(mpq_tmp)

    try:
        init_short_digraph(g, G, edge_labelled = False)
        init_reverse(bfs_dag, g)

        queue               = <uint32_t *> check_malloc(n*sizeof(uint32_t))
        degrees             = <uint32_t *> check_malloc(n*sizeof(uint32_t))
        n_paths_from_source = <numerical_type *> check_malloc(n*sizeof(numerical_type))
        betweenness_source  = <numerical_type *> check_malloc(n*sizeof(numerical_type))
        betweenness         = <numerical_type *> check_malloc(n*sizeof(numerical_type))

        bitset_init(seen,n)
        bitset_init(next_layer,n)

        if numerical_type is double:
            memset(betweenness,0,n*sizeof(double))
        else:
            for i in range(n):
                mpq_init(betweenness[i])
                mpq_set_ui(betweenness[i],0,1)
                mpq_init(betweenness_source[i])
                mpq_init(n_paths_from_source[i])

        for source in range(n):

            if numerical_type is double:
                memset(betweenness_source ,0,n*sizeof(double))
                memset(n_paths_from_source,0,n*sizeof(double))
                n_paths_from_source[source]=1
            else:
                for i in range(n):
                    mpq_set_ui(betweenness_source[i] ,0,1)
                    mpq_set_ui(n_paths_from_source[i],0,1)
                mpq_set_ui(n_paths_from_source[source],1,1)

            # initialize data
            bitset_set_first_n(seen, 0)
            bitset_add(seen,source)
            bitset_set_first_n(next_layer, 0)

            memset(degrees,0,n*sizeof(uint32_t))

            queue[0] = source
            layer_current_beginning = 0
            layer_current_end       = 1
            layer_next_end          = 1

            # The number of shortest paths from 'source' to every other vertex.
            #
            # It is a BFS. The graph is explored layer by layer.
            while layer_current_beginning<layer_current_end:

                # Looking for all non-discovered neighbors of some vertex of the
                # current layer.
                for j in range(layer_current_beginning,layer_current_end):
                    u = queue[j]

                    # List the neighors of u
                    p_tmp = g.neighbors[u]
                    while p_tmp<g.neighbors[u+1]:
                        v = p_tmp[0]
                        p_tmp += 1

                        # Is it a new vertex ?
                        if bitset_in(seen,v):
                            continue

                        # Is it the first time we see it ?
                        elif not bitset_in(next_layer,v):
                            bitset_add(next_layer,v)
                            queue[layer_next_end] = v
                            layer_next_end += 1

                        # update the count of paths and the BFS dag.
                        bfs_dag.neighbors[v][degrees[v]] = u
                        degrees[v] += 1
                        if numerical_type is double:
                            n_paths_from_source[v] += n_paths_from_source[u]
                        else:
                            mpq_add(n_paths_from_source[v],n_paths_from_source[v],n_paths_from_source[u])

                # 'next_layer' becomes 'current_layer'
                for j in range(layer_current_end, layer_next_end):
                    bitset_add(seen,queue[j])

                layer_current_beginning = layer_current_end
                layer_current_end       = layer_next_end

            # Compute the betweenness from the number of paths
            #
            # We enumerate vertices in reverse order of discovery.
            for i in range(layer_current_end-1,-1,-1):
                u = queue[i]
                for j in range(degrees[u]):
                    v = bfs_dag.neighbors[u][j]
                    if v != source: # better to not 'if' but set it to 0 afterwards?
                        if numerical_type is double:
                            betweenness_source[v] += (betweenness_source[u]+1)*(n_paths_from_source[v]/n_paths_from_source[u])
                        else:
                            mpq_set_ui(mpq_tmp,1,1)
                            mpq_add(mpq_tmp,betweenness_source[u],mpq_tmp)
                            mpq_mul(mpq_tmp,mpq_tmp,n_paths_from_source[v])
                            mpq_div(mpq_tmp,mpq_tmp,n_paths_from_source[u])
                            mpq_add(betweenness_source[v],betweenness_source[v],mpq_tmp)

            # update betweenness from betweenness_source
            for i in range(n):
                if numerical_type is double:
                    betweenness[i] += betweenness_source[i]
                else:
                    mpq_add(betweenness[i],betweenness[i],betweenness_source[i])

            sig_check() # check for KeyboardInterrupt

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
        sage_free(queue)
        sage_free(n_paths_from_source)
        sage_free(degrees)
        sage_free(betweenness_source)
        sage_free(betweenness)

    if not G.is_directed():
        betweenness_list = [x/2 for x in betweenness_list]

    if normalize:
        if G.is_directed():
            betweenness_list = [  x/((n-1)*(n-2)) for x in betweenness_list]
        else:
            betweenness_list = [2*x/((n-1)*(n-2)) for x in betweenness_list]


    return {vv:betweenness_list[i] for i,vv in enumerate(G.vertices())}
