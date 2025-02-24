# cython: binding=True
r"""
Hyperbolicity

**Definition** :

    The hyperbolicity `\delta` of a graph `G` has been defined by Gromov
    [Gro1987]_ as follows (we give here the so-called 4-points condition):

      Let `a, b, c, d` be vertices of the graph, let `S_1`, `S_2` and `S_3` be
      defined by

      .. MATH::

          S_1 = dist(a, b) + dist(d, c)\\
          S_2 = dist(a, c) + dist(b, d)\\
          S_3 = dist(a, d) + dist(b, c)\\

      and let `M_1` and `M_2` be the two largest values among `S_1`, `S_2`, and
      `S_3`. We define `hyp(a, b, c, d) = M_1 - M_2`, and the hyperbolicity
      `\delta(G)` of the graph is the maximum of `hyp` over all possible
      4-tuples `(a, b, c, d)` divided by 2. That is, the graph is said
      `\delta`-hyperbolic when

      .. MATH::

          \delta(G) = \frac{1}{2}\max_{a,b,c,d\in V(G)}hyp(a, b, c, d)

      (note that `hyp(a, b, c, d)=0` whenever two elements among `a,b,c,d` are
      equal)

**Some known results** :

    - Trees and cliques are `0`-hyperbolic

    - `n\times n` grids are `n-1`-hyperbolic

    - Cycles are approximately `n/4`-hyperbolic

    - Chordal graphs are `\leq 1`-hyperbolic

    Besides, the hyperbolicity of a graph is the maximum over all its
    biconnected components.

**Algorithms and complexity** :

    The time complexity of the naive implementation (i.e. testing all 4-tuples)
    is `O( n^4 )`, and an algorithm with time complexity `O(n^{3.69})` has been
    proposed in [FIV2012]_.  This remains very long for large-scale graphs, and
    much harder to implement.

    Several improvements over the naive algorithm have been proposed and are
    implemented in the current module.

    - Another upper bound on `hyp(a, b, c, d)` has been proved in [CCL2015]_. It
      is used to design an algorithm with worse case time complexity in
      `O(n^4)` but that behaves much better in practice.

      Assume that `S_1 = dist(a, b) + dist(c, d)` is the largest sum among
      `S_1,S_2,S_3`. We have

      .. MATH::

          S_2 + S_3 =& dist(a, c) + dist(b, d) + dist(a, d) + dist(b, c)\\
          =& [ dist(a, c) + dist(b, c) ] + [ dist(a, d) + dist(b, d)]\\
          \geq &dist(a,b) + dist(a,b)\\
          \geq &2dist(a,b)\\

      Now, since `S_1` is the largest sum, we have

      .. MATH::

          hyp(a, b, c, d) =& S_1 - \max\{S_2, S_3\}\\
          \leq& S_1 - \frac{S_2+ S_3}{2}\\
          \leq& S_1 - dist(a, b)\\
          =& dist(c, d)\\

      We obtain similarly that `hyp(a, b, c, d) \leq dist(a, b)`. Consequently,
      in the implementation of the 'CCL' algorithm, we ensure that `S_1` is
      larger than `S_2` and `S_3` using an ordering of the pairs by decreasing
      lengths. Then, we use the best value `h` found so far to stop exploration
      as soon as `dist(a, b) \leq h`.

      The worst case time complexity of this algorithm is `O(n^4)`, but it
      performs very well in practice since it cuts the search space. This
      algorithm can be turned into an approximation algorithm since at any step
      of its execution we maintain an upper and a lower bound. We can thus stop
      execution as soon as a multiplicative approximation factor or an additive
      one is proven.

    - The notion of ''far-apart pairs'' has been introduced in [Sot2011]_ to
      further reduce the number of 4-tuples to consider. We say that the pair
      `(a,b)` is far-apart if for every `w` in `V\setminus\{a,b\}` we have

      .. MATH::

          dist(w,a)+dist(a,b) > dist(w,b) \text{ and }dist(w,b)+dist(a,b) > dist(w,a)

      Determining the set of far-apart pairs can be done in time `O(nm)` using
      BFS. Now, it is proved in [Sot2011]_ that there exists two far-apart pairs
      `(a,b)` and `(c,d)` satisfying `\delta(G) = hyp(a, b, c, d)/2`. For
      instance, the `n\times m`-grid has only two far-apart pairs, and so
      computing its hyperbolicity is immediate once the far-apart pairs are
      found. The 'CCL+FA' or 'CCL+' algorithm improves the 'CCL' algorithm
      since it uses far-apart pairs.

    - This algorithm was further improved in [BCCM2015]_: instead of iterating
      twice over all pairs of vertices, in the "inner" loop, we cut several
      pairs by exploiting properties of the underlying graph.

.. TODO::

    - Add exact methods for the hyperbolicity of chordal graphs

    - Add method for partitioning the graph with clique separators

**This module contains the following functions**

At Python level :

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~hyperbolicity` | Return the hyperbolicity of the graph or an approximation of this value.
    :meth:`~hyperbolicity_distribution` | Return the hyperbolicity distribution of the graph or a sampling of it.

AUTHORS:

- David Coudert (2012): initial version, exact and approximate algorithm,
  distribution, sampling
- David Coudert (2014): improved exact algorithm using far-apart pairs
- Michele Borassi (2015): cleaned the code and implemented the new algorithm
- Karan Desai (2016): fixed  minor typo in documentation


Methods
-------
"""

# ****************************************************************************
#       Copyright (C) 2012 David Coudert <david.coudert@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from libc.string cimport memset
from cysignals.memory cimport check_allocarray, sig_free
from cysignals.signals cimport sig_on, sig_off
from memory_allocator cimport MemoryAllocator

from sage.graphs.distances_all_pairs cimport c_distances_all_pairs
from sage.arith.misc import binomial
from sage.rings.integer_ring import ZZ
from sage.graphs.base.static_sparse_graph cimport short_digraph
from sage.graphs.base.static_sparse_graph cimport init_short_digraph
from sage.graphs.base.static_sparse_graph cimport free_short_digraph
from libc.stdint cimport uint16_t, uint32_t, uint64_t
from sage.data_structures.bitset_base cimport *


# Defining a pair of vertices as a C struct
ctypedef struct pair:
    uint32_t s
    uint32_t t


######################################################################
# Speedup functions
######################################################################

def _my_subgraph(G, vertices, relabel=False, return_map=False):
    r"""
    Return the subgraph containing the given vertices.

    This method considers only the connectivity. Therefore, edge labels are
    ignored as well as any other decoration of the graph (vertex position,
    etc.).

    If ``relabel`` is ``True``, the vertices of the new graph are relabeled with
    integers in the range '0\cdots \mid vertices \mid -1'. The relabeling map is
    returned if ``return_map`` is also ``True``.

    TESTS:

    Giving anything else than a Graph::

        sage: from sage.graphs.hyperbolicity import _my_subgraph as mysub
        sage: mysub([], [])
        Traceback (most recent call last):
        ...
        ValueError: the input parameter must be a Graph

    Subgraph of a PetersenGraph::

        sage: from sage.graphs.hyperbolicity import _my_subgraph as mysub
        sage: H = mysub(graphs.PetersenGraph(), [0,2,4,6])
        sage: H.edges(sort=True, labels=None)
        [(0, 4)]
        sage: H.vertices(sort=True)
        [0, 2, 4, 6]
    """
    from sage.graphs.graph import Graph
    if not isinstance(G, Graph):
        raise ValueError("the input parameter must be a Graph")
    H = Graph()
    if not vertices:
        return (H, {}) if (relabel and return_map) else H

    if relabel:
        map = dict(zip(iter(vertices), range(len(vertices))))
    else:
        map = dict(zip(iter(vertices), iter(vertices)))

    B = {}
    for v in G.vertex_iterator():
        B[v] = False
    for v in vertices:
        B[v] = True
        H.add_vertex(map[v])

    for u in vertices:
        for v in G.neighbor_iterator(u):
            if B[v]:
                H.add_edge(map[u], map[v])

    return (H, map) if (relabel and return_map) else H


######################################################################
# Building blocks
######################################################################

cdef inline int __hyp__(unsigned short** distances, int a, int b, int c, int d) noexcept:
    """
    Return the hyperbolicity of the given 4-tuple.
    """
    cdef int S1, S2, S3, h
    S1 = distances[a][b] + distances[c][d]
    S2 = distances[a][c] + distances[b][d]
    S3 = distances[a][d] + distances[b][c]
    if S1 >= S2:
        if S2 > S3:
            h = S1 - S2
        else:
            h = abs(S1 - S3)
    else:
        if S1 > S3:
            h = S2 - S1
        else:
            h = abs(S2 - S3)
    return h


######################################################################
# Basic algorithm for the hyperbolicity
######################################################################

cdef tuple hyperbolicity_basic_algorithm(int N,
                                         unsigned short** distances,
                                         verbose):
    """
    Return **twice** the hyperbolicity of a graph, and a certificate.

    This method implements the basic algorithm for computing the hyperbolicity
    of a graph which tests all 4-tuples of vertices not satisfying a cutting
    rule proposed in [Sot2011]_.

    INPUT:

    - ``N`` -- number of vertices of the graph

    - ``distances`` -- path distance matrix (see the distance_all_pairs
      module)

    - ``verbose`` -- boolean (default: ``False``); set to ``True`` to display
      some information during execution

    OUTPUT:

    This function returns a tuple ( h, certificate ), where:

    - ``h`` -- the maximum computed value over all 4-tuples, and so is twice
      the hyperbolicity of the graph. If no such 4-tuple is found, -1 is
      returned.

    - ``certificate`` -- 4-tuple of vertices maximizing the value `h`. If no
      such 4-tuple is found, the empty list [] is returned.
    """
    cdef int a, b, c, d, hh, h_LB
    cdef list certificate

    h_LB = -1

    for a in range(N - 3):
        for b in range(a + 1, N - 2):

            # We use the cutting rule proposed in [Sot2011]_
            if 2 * distances[a][b] <= h_LB:
                continue

            for c in range(b + 1, N - 1):

                # We use the cutting rule proposed in [Sot2011]_
                if 2 * distances[a][c] <= h_LB or 2 * distances[b][c] <= h_LB:
                    continue

                for d in range(c + 1, N):

                    # We compute the hyperbolicity of the 4-tuple
                    hh = __hyp__(distances, a, b, c, d)

                    # We compare the value with previously known bound
                    if hh > h_LB:
                        h_LB = hh
                        certificate = [a, b, c, d]

                        if verbose:
                            print('New lower bound:', ZZ(hh)/2)

    # Last, we return the computed value and the certificate
    if h_LB != -1:
        return (h_LB, certificate)
    return (-1, [])


######################################################################
# Greedy dominating set
######################################################################

def _greedy_dominating_set(H, verbose=False):
    r"""
    Return a greedy approximation of a dominating set.

    EXAMPLES::

        sage: from sage.graphs.hyperbolicity import _greedy_dominating_set
        sage: G = graphs.PetersenGraph()
        sage: _greedy_dominating_set(G)
        [0, 2, 6]
    """
    cdef list V = sorted([(d, u) for u, d in H.degree_iterator(labels=True)],
                         reverse=True, key=lambda x: x[0])
    cdef list DOM = []
    cdef set seen = set()
    for _, u in V:
        if u not in seen:
            seen.add(u)
            DOM.append(u)
            seen.update(H.neighbor_iterator(u))

    if verbose:
        print("Greedy dominating set: {}".format(DOM))

    return DOM


######################################################################
# Distances and far-apart pairs
######################################################################

cdef inline distances_and_far_apart_pairs(gg,
                                          unsigned short* distances,
                                          unsigned short* far_apart_pairs,
                                          list int_to_vertex):
    """
    Compute both distances between all pairs and far-apart pairs.

    See the module's documentation for the definition of far-apart pairs.

    This method assumes that:

    - The input graph gg is connected. If not, the result will be incorrect.

    - The arrays distances and far_apart_pairs have already been allocated with
      size `n^2`.
    """
    cdef uint32_t n = gg.order()
    cdef uint32_t i

    if not distances or not far_apart_pairs:
        raise ValueError("distances or far_apart_pairs is a NULL pointer")
    elif n > <unsigned short> -1:
        # Computing the distances/far_apart_pairs can only be done if we have
        # less than MAX_UNSIGNED_SHORT vertices.
        raise ValueError("The graph backend contains more than {} nodes and "
                         "we cannot compute the matrix of distances/far-apart "
                         "pairs on something"
                         "like that!".format(<unsigned short> -1))

    # The list of waiting vertices
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef uint32_t* waiting_list = <uint32_t*>mem.allocarray(n, sizeof(uint32_t))
    cdef unsigned short** c_far_apart = <unsigned short**>mem.allocarray(n, sizeof(unsigned short*))

    # The vertices which have already been visited
    cdef bitset_t seen
    bitset_init(seen, n)

    # the beginning and the end of the list stored in waiting_list
    cdef uint32_t waiting_beginning, waiting_end

    cdef uint32_t source
    cdef uint32_t v, u

    # All pairs are initially far-apart
    memset(far_apart_pairs, 1, n * n * sizeof(unsigned short))
    for i in range(n):
        c_far_apart[i] = far_apart_pairs + i * n
        c_far_apart[i][i] = 0

    # Copying the whole graph to obtain the list of neighbors quicker than by
    # calling out_neighbors. This data structure is well documented in the
    # module sage.graphs.base.static_sparse_graph
    cdef short_digraph sd
    init_short_digraph(sd, gg, edge_labelled=False, vertex_list=int_to_vertex)
    cdef uint32_t** p_vertices = sd.neighbors
    cdef uint32_t* p_tmp
    cdef uint32_t* end

    cdef unsigned short* c_distances = distances

    memset(distances, -1, n * n * sizeof(unsigned short))

    # We run n different BFS taking each vertex as a source
    for source in range(n):

        # The source is seen
        bitset_clear(seen)
        bitset_add(seen, source)
        c_distances[source] = 0

        # and added to the queue
        waiting_list[0] = source
        waiting_beginning = 0
        waiting_end = 0

        # For as long as there are vertices left to explore
        while waiting_beginning <= waiting_end:

            # We pick the first one
            v = waiting_list[waiting_beginning]

            p_tmp = p_vertices[v]
            end = p_vertices[v+1]

            # Iterating over all the outneighbors u of v
            while p_tmp < end:
                u = p_tmp[0]

                # If we notice one of these neighbors is not seen yet, we set
                # its parameters and add it to the queue to be explored later.
                if not bitset_in(seen, u):
                    c_distances[u] = c_distances[v] + 1
                    bitset_add(seen, u)
                    waiting_end += 1
                    waiting_list[waiting_end] = u

                if c_distances[u] == c_distances[v] + 1:
                    # v is on the path from source to u
                    c_far_apart[source][v] = 0
                    c_far_apart[v][source] = 0

                p_tmp += 1

            waiting_beginning += 1

        c_distances += n

    bitset_free(seen)
    free_short_digraph(sd)


cdef inline pair** sort_pairs(uint32_t N,
                              uint16_t D,
                              unsigned short** values,
                              unsigned short** to_include,
                              uint32_t* nb_p,
                              uint32_t* nb_pairs_of_length) noexcept:
    """
    Return an array of unordered pairs {i,j} in increasing order of values.

    Uses counting sort to list pairs {i,j} in increasing order of values(i,j).
    If to_include[i][j] = 0, the pair is ignored. We assume N and D to be
    correct with respect to the arrays values and to_include, that values and
    to_include are symmetric (that is, values[i][j] = values[j][i] and
    to_include[i][j] = to_include[j][i], and that nb_p, nb_pairs_of_length are
    already allocated.

    INPUT:

    - ``N`` -- the range of i and j (that is, the square root of the number
      of pairs to be sorted);

    - ``D`` -- the maximum value of an element;

    - ``values`` -- an array containing in position (i,j) the value of the
      pair (i,j);

    - ``to_include`` -- an array such that to_include[i][j] contains "1" if
      pair (i,j) should be included, "0" otherwise. If NULL, all elements are
      included;

    OUTPUT:

    - ``nb_p`` -- the number of pairs to be included;

    - ``nb_pairs_of_length`` -- an array containing in position k the number
      of pairs (i,j) that are included and such that values[i][j] = k

    - ``pairs_of_length`` -- this function returns this array, containing in
      position k a pointer to the first included pair (i,j) such that
      values[i][j] = k.
    """
    # pairs_of_length[d] is the list of pairs of vertices at distance d
    cdef pair** pairs_of_length = <pair**>check_allocarray(D + 1, sizeof(pair*))
    cdef unsigned short* p_to_include
    cdef uint32_t i, j, k
    nb_p[0] = 0

    # fills nb_pairs_of_length and nb_p
    memset(nb_pairs_of_length, 0, (D + 1) * sizeof(uint32_t))

    if not to_include:
        nb_p[0] = (N * (N - 1)) / 2
        for i in range(N):
            for j in range(i + 1, N):
                nb_pairs_of_length[values[i][j]] += 1
    else:
        for i in range(N):
            p_to_include = to_include[i]
            for j in range(i + 1, N):
                if p_to_include[j]:
                    nb_p[0] += 1
                    nb_pairs_of_length[values[i][j]] += 1

    pairs_of_length[0] = <pair*>check_allocarray(nb_p[0], sizeof(pair))

    # temporary variable used to fill pairs_of_length
    cdef uint32_t* cpt_pairs = <uint32_t*>check_calloc(D + 1, sizeof(uint32_t))

    # ==> Defines pairs_of_length[d] for all d
    for i in range(1, D + 1):
        pairs_of_length[i] = pairs_of_length[i - 1] + nb_pairs_of_length[i - 1]

    # ==> Fills pairs_of_length[d] for all d
    if not to_include:
        for i in range(N):
            for j in range(i + 1, N):
                k = values[i][j]
                if k:
                    pairs_of_length[k][cpt_pairs[k]].s = i
                    pairs_of_length[k][cpt_pairs[k]].t = j
                    cpt_pairs[k] += 1
    else:
        for i in range(N):
            p_to_include = to_include[i]
            for j in range(i + 1, N):
                if p_to_include[j]:
                    k = values[i][j]
                    pairs_of_length[k][cpt_pairs[k]].s = i
                    pairs_of_length[k][cpt_pairs[k]].t = j
                    cpt_pairs[k] += 1

    sig_free(cpt_pairs)
    return pairs_of_length


######################################################################
# Compute the hyperbolicity using the algorithm of [BCCM2015]_
######################################################################

cdef tuple hyperbolicity_BCCM(int N,
                              unsigned short** distances,
                              unsigned short** far_apart_pairs,
                              int D,
                              int h_LB,
                              float approximation_factor,
                              float additive_gap,
                              verbose=False):
    """
    Return the hyperbolicity of a graph.

    This method implements the exact and the approximate algorithms proposed in
    [BCCM2015]_. See the module's documentation for more details.

    This method assumes that the graph under consideration is connected.

    INPUT:

    - ``N`` -- number of vertices of the graph

    - ``distances`` -- path distance matrix

    - ``far_apart_pairs`` -- 0/1 matrix of far-apart pairs. Pair ``(i,j)`` is
      far-apart if ``far_apart_pairs[i][j]\neq 0``

    - ``D`` -- diameter of the graph

    - ``h_LB`` -- lower bound on the hyperbolicity

    - ``approximation_factor`` -- when the approximation factor is set to some
      value larger than 1.0, the function stop computations as soon as the
      ratio between the upper bound and the best found solution is less than
      the approximation factor. When the approximation factor is 1.0, the
      problem is solved optimally.

    - ``additive_gap`` -- when set to a positive number, the function stop
      computations as soon as the difference between the upper bound and the
       best found solution is less than additive gap. When the gap is 0.0, the
       problem is solved optimally.

    - ``verbose`` -- boolean (default: ``False``); set to ``True`` to display
      some information during execution

    OUTPUT:

    This function returns a tuple ( h, certificate, h_UB ), where:

    - ``h`` -- integer; when 4-tuples with hyperbolicity larger or equal
     to `h_LB are found, h is the maximum computed value and so twice the
     hyperbolicity of the graph. If no such 4-tuple is found, it returns -1.

    - ``certificate`` -- is a list of vertices; when 4-tuples with
      hyperbolicity larger that h_LB are found, certificate is the list of the
      4 vertices for which the maximum value (and so the hyperbolicity of the
      graph) has been computed. If no such 4-tuple is found, it returns the
      empty list [].

    - ``h_UB`` -- integer equal to the proven upper bound for `h`; when
      ``h == h_UB``, the returned solution is optimal
    """
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef int h = 0, hh  # can get negative value
    cdef int a, b, c, d, h_UB, n_val, n_acc, i, j
    cdef int hplusone
    cdef int condacc
    cdef int x, S1, S2, S3
    cdef list certificate = []
    cdef uint32_t nb_p  # The total number of pairs.
    cdef unsigned short *dist_a
    cdef unsigned short *dist_b
    cdef bint GOTO_RETURN = 0

    # Variable used to store "mates".
    cdef int **mate = <int**> mem.malloc(N * sizeof(int*))
    for i in range(N):
        mate[i] = <int*> mem.malloc(N * sizeof(int))
    cdef int *cont_mate = <int*> mem.calloc(N, sizeof(int))

    # The farness of all vertices (the farness of v is the sum of the distances
    # between v and all other vertices).
    cdef uint64_t *farness = <uint64_t*> mem.calloc(N, sizeof(uint64_t))
    cdef short *ecc = <short*> mem.calloc(N, sizeof(short))
    cdef int central = 0
    cdef int **mates_decr_order_value = <int**> mem.malloc(N * sizeof(int*))
    cdef int *value = <int*> mem.malloc(N * sizeof(int))
    cdef int *nvalues = <int*> mem.malloc((D + 1) * sizeof(int))
    cdef short *acc_bool = <short*> mem.calloc(N, sizeof(short))
    cdef int *acc = <int*> mem.malloc(N * sizeof(int))
    cdef int *val = <int*> mem.malloc(N * sizeof(int))
    cdef int *nvalues_cum = <int*> mem.malloc((D + 1) * sizeof(int))
    cdef uint64_t nq = 0

    # We compute the farness and the eccentricity of all vertices.
    # We set central as the vertex with minimum farness
    for a in range(N):
        dist_a = distances[a]
        for b in range(N):
            farness[a] += dist_a[b]
            ecc[a] = max(ecc[a], dist_a[b])
            if dist_a[b] >= N:
                raise ValueError("the input graph must be connected")
        if farness[a] < farness[central]:
            central = a
    cdef unsigned short *dist_central = distances[central]

    # We put in variable mates_decr_order_value[a] all vertices b, in
    # decreasing order of ecc[b]-distances[a][b]
    for a in range(N):
        mates_decr_order_value[a] = <int*> mem.malloc(N * sizeof(int))
        dist_a = distances[a]
        memset(nvalues, 0, (D + 1) * sizeof(int))

        for b in range(N):
            value[b] = ecc[b] - dist_a[b]
            nvalues[value[b]] += 1
        nvalues_cum[D] = 0

        for b in range(D - 1, -1, -1):
            nvalues_cum[b] = nvalues_cum[b + 1] + nvalues[b + 1]

        for b in range(N):
            mates_decr_order_value[a][nvalues_cum[value[b]]] = b
            nvalues_cum[value[b]] += 1

    # We sort pairs, in increasing order of distance
    cdef uint32_t * nb_pairs_of_length = <uint32_t *> mem.calloc(D + 1, sizeof(uint32_t))

    cdef pair ** pairs_of_length = sort_pairs(N, D, distances, far_apart_pairs,
                                              &nb_p, nb_pairs_of_length)

    if verbose:
        print("Current 2 connected component has %d vertices and diameter %d" % (N, D))
        if not far_apart_pairs:
            print("Number of pairs: %d" % (nb_p))
            print("Repartition of pairs:",
                  [(i, nb_pairs_of_length[i]) for i in range(1, D + 1) if nb_pairs_of_length[i] > 0])
        else:
            print("Number of far-apart pairs: %d\t(%d pairs in total)" % (nb_p, binomial(N, 2)))
            print("Repartition of far-apart pairs:",
                  [(i, nb_pairs_of_length[i]) for i in range(1, D + 1) if nb_pairs_of_length[i] > 0])

    cdef pair * sorted_pairs = pairs_of_length[0]

    approximation_factor = min(approximation_factor, D)
    additive_gap = min(additive_gap, D)

    # We start iterating from pairs with maximum distance.
    for x in range(nb_p - 1, -1, -1):
        a = sorted_pairs[x].s
        b = sorted_pairs[x].t

        # Without loss of generality, a has smaller farness than b.
        if farness[a] < farness[b]:
            a, b = b, a

        dist_a = distances[a]
        dist_b = distances[b]
        h_UB = distances[a][b]

        # If we cannot improve further, we stop
        if h_UB <= h:
            h_UB = h
            GOTO_RETURN = 1
            break

        # Termination if required approximation is found
        if (h_UB <= h * approximation_factor) or (h_UB - h <= additive_gap):
            GOTO_RETURN = 1
            break

        # We update variable mate, adding pair (a,b)
        mate[a][cont_mate[a]] = b
        cont_mate[a] += 1
        mate[b][cont_mate[b]] = a
        cont_mate[b] += 1

        # We compute acceptable and valuable vertices
        n_acc = 0
        n_val = 0

        hplusone = h + 1
        condacc = 3 * hplusone - 2 * h_UB

        for i in range(N):
            c = mates_decr_order_value[a][i]
            if cont_mate[c] > 0:
                if 2 * (ecc[c] - dist_a[c]) >= condacc:
                    if 2 * (ecc[c] - dist_b[c]) >= condacc:
                        if 2 * dist_a[c] >= hplusone and 2 * dist_b[c] >= hplusone:
                            if 2 * ecc[c] >= 2 * hplusone - h_UB + dist_a[c] + dist_b[c]:
                                # Vertex c is acceptable
                                acc_bool[c] = 1
                                acc[n_acc] = c
                                n_acc += 1
                                if 2 * dist_central[c] + h_UB - h > dist_a[c] + dist_b[c]:
                                    # Vertex c is valuable
                                    val[n_val] = c
                                    n_val += 1
                else:
                    break

        # For each pair (c,d) where c is valuable and d is acceptable, we
        # compute the hyperbolicity of (a,b,c,d), and we update h if necessary
        for i in range(n_val):
            c = val[i]
            for j in range(cont_mate[c]):
                d = mate[c][j]
                if (acc_bool[d]):
                    nq += 1
                    S1 = h_UB + distances[c][d]
                    S2 = dist_a[c] + dist_b[d]
                    S3 = dist_a[d] + dist_b[c]
                    if S2 > S3:
                        hh = S1 - S2
                    else:
                        hh = S1 - S3

                    if h < hh or not certificate:
                        # We update current bound on the hyperbolicity and the
                        # search space.
                        #
                        # Note that if hh==0, we first make sure that a,b,c,d are
                        # all distinct and are a valid certificate.
                        if hh > 0 or not (a == c or a == d or b == c or b == d):
                            h = hh
                            certificate = [a, b, c, d]

                            if verbose:
                                print("New lower bound:", ZZ(hh)/2)

        # We reset acc_bool
        for v in range(n_acc):
            acc_bool[acc[v]] = 0

    # Needed because sometimes h_UB is not updated, if the analysis is no cut.
    if not GOTO_RETURN:
        h_UB = h

    # We now free the memory
    sig_free(pairs_of_length[0])
    sig_free(pairs_of_length)

    if verbose:
        print("Visited 4-tuples:", nq)

    # Last, we return the computed value and the certificate
    if not certificate:
        return (-1, [], h_UB)

    # When using far-apart pairs, the loops may end before improving the
    # upper-bound
    return (h, certificate, h_UB)


######################################################################
# Compute the hyperbolicity using the algorithm of [CCL2015]_
######################################################################

cdef tuple hyperbolicity_CCL(int N,
                             unsigned short** distances,
                             unsigned short** far_apart_pairs,
                             int D,
                             int h_LB,
                             float approximation_factor,
                             float additive_gap,
                             verbose=False):
    """
    Return the hyperbolicity of a graph.

    This method implements the exact and the approximate algorithms proposed in
    [CCL2015]_. See the module's documentation for more details.

    This method assumes that the graph under consideration is connected.

    INPUT:

    - ``N`` -- number of vertices of the graph

    - ``distances`` -- path distance matrix

    - ``far_apart_pairs`` -- 0/1 matrix of far-apart pairs. Pair ``(i,j)`` is
      far-apart if ``far_apart_pairs[i][j]\neq 0``

    - ``D`` -- diameter of the graph

    - ``h_LB`` -- lower bound on the hyperbolicity

    - ``approximation_factor`` -- when the approximation factor is set to some
      value larger than 1.0, the function stop computations as soon as the
      ratio between the upper bound and the best found solution is less than
      the approximation factor. When the approximation factor is 1.0, the
      problem is solved optimally.

    - ``additive_gap`` -- when set to a positive number, the function stop
      computations as soon as the difference between the upper bound and the
       best found solution is less than additive gap. When the gap is 0.0, the
       problem is solved optimally.

    - ``verbose`` -- boolean (default: ``False``); set to ``True`` to display
      some information during execution

    OUTPUT:

    This function returns a tuple ( h, certificate, h_UB ), where:

    - ``h`` -- integer; when 4-tuples with hyperbolicity larger or equal
     to `h_LB are found, h is the maximum computed value and so twice the
     hyperbolicity of the graph. If no such 4-tuple is found, it returns -1.

    - ``certificate`` -- is a list of vertices; when 4-tuples with
      hyperbolicity larger that h_LB are found, certificate is the list of the
      4 vertices for which the maximum value (and so the hyperbolicity of the
      graph) has been computed. If no such 4-tuple is found, it returns the
      empty list [].

    - ``h_UB`` -- integer equal to the proven upper bound for `h`; when
      ``h == h_UB``, the returned solution is optimal
    """
    cdef int hh  # can get negative value
    cdef int a, b, c, d, h, h_UB
    cdef int l1, l2, S1, S2, S3
    cdef uint32_t x, y
    cdef list certificate = []
    cdef uint32_t nb_p  # The total number of pairs

    # Test if the distance matrix corresponds to a connected graph, i.e., if
    # distances from node 0 are all less or equal to N-1.
    for a in range(N):
        if distances[0][a] >= N:
            raise ValueError("the input graph must be connected")

    # nb_pairs_of_length[d] is the number of pairs of vertices at distance d
    cdef uint32_t* nb_pairs_of_length = <uint32_t*>check_allocarray(D + 1, sizeof(uint32_t))

    if not nb_pairs_of_length:
        raise MemoryError

    cdef pair** pairs_of_length = sort_pairs(N, D, distances, far_apart_pairs,
                                             &nb_p, nb_pairs_of_length)

    if verbose:
        print("Current 2 connected component has %d vertices and diameter %d" % (N, D))
        if not far_apart_pairs:
            print("Number of pairs: %d" % (nb_p))
            print("Repartition of pairs:",
                  [(i, nb_pairs_of_length[i]) for i in range(1, D + 1) if nb_pairs_of_length[i] > 0])
        else:
            print("Number of far-apart pairs: %d\t(%d pairs in total)" % (nb_p, binomial(N, 2)))
            print("Repartition of far-apart pairs:",
                  [(i, nb_pairs_of_length[i]) for i in range(1, D + 1) if nb_pairs_of_length[i] > 0])

    approximation_factor = min(approximation_factor, D)
    additive_gap = min(additive_gap, D)

    # We create the list of triples (sum,length1,length2) sorted in decreasing
    # lexicographic order: decreasing by sum, decreasing by length2, decreasing
    # length1. This is to ensure a valid ordering for S1, to avoid some tests,
    # and to ease computation of bounds.
    cdef list triples = []
    for l2 in range(D, 0, -1):
        if nb_pairs_of_length[l2]:
            for l1 in range(D, l2 - 1, -1):
                if nb_pairs_of_length[l1]:
                    triples.append((l1 + l2, l1, l2))

    # We use some short-cut variables for efficiency
    cdef pair* pairs_of_length_l1
    cdef pair* pairs_of_length_l2
    cdef uint32_t nb_pairs_of_length_l1, nb_pairs_of_length_l2
    cdef unsigned short * dist_a
    cdef unsigned short * dist_b
    h = h_LB
    h_UB = D
    cdef int GOTO_RETURN = 0

    # S1 = l1+l2
    # l1 = dist(a,b)
    # l2 = dist(c,d)
    # l1 >= l2
    for S1, l1, l2 in triples:

        if h_UB > l2:
            h_UB = l2

            if verbose:
                print("New upper bound:", ZZ(h_UB) / 2)

        # Termination if required approximation is found
        if certificate and ((h_UB <= h * approximation_factor) or (h_UB-h <= additive_gap)):
            GOTO_RETURN = 1
            break

        # If we cannot improve further, we stop
        #
        # See the module's documentation for a proof that this cut is
        # valid. Remember that the triples are sorted in a specific order.
        if h_UB <= h:
            h_UB = h
            break

        pairs_of_length_l1 = pairs_of_length[l1]
        pairs_of_length_l2 = pairs_of_length[l2]
        nb_pairs_of_length_l1 = nb_pairs_of_length[l1]
        nb_pairs_of_length_l2 = nb_pairs_of_length[l2]

        for x in range(nb_pairs_of_length_l1):
            a = pairs_of_length_l1[x].s
            b = pairs_of_length_l1[x].t
            dist_a = distances[a]
            dist_b = distances[b]

            # We do not want to test pairs of pairs twice if l1 == l2
            for y in range((x + 1) if l1 == l2 else 0, nb_pairs_of_length_l2):
                c = pairs_of_length_l2[y].s
                d = pairs_of_length_l2[y].t

                # We compute the hyperbolicity of the 4-tuple. We have S1 = l1 +
                # l2, and the order in which pairs are visited allow us to claim
                # that S1 = max( S1, S2, S3 ). Indeed, if S1 is not the maximum
                # value, the order ensures that the maximum value has previously
                # been checked.
                S2 = dist_a[c] + dist_b[d]
                S3 = dist_a[d] + dist_b[c]
                if S2 > S3:
                    hh = S1 - S2
                else:
                    hh = S1 - S3

                if h < hh or not certificate:
                    # We update current bound on the hyperbolicity and the
                    # search space.
                    #
                    # Note that if hh==0, we first make sure that a,b,c,d are
                    # all distinct and are a valid certificate.
                    if hh > 0 or not (a == c or a == d or b == c or b == d):
                        h = hh
                        certificate = [a, b, c, d]

                        if verbose:
                            print("New lower bound:", ZZ(hh) / 2)

                        # If we cannot improve further, we stop
                        if l2 <= h:
                            GOTO_RETURN = 1
                            h_UB = h
                            break

                        # Termination if required approximation is found
                        if (h_UB <= h * approximation_factor) or (h_UB - h <= additive_gap):
                            GOTO_RETURN = 1
                            break

            if GOTO_RETURN:
                break

        if GOTO_RETURN:
            break

    # We now free the memory
    sig_free(nb_pairs_of_length)
    sig_free(pairs_of_length[0])
    sig_free(pairs_of_length)

    # Last, we return the computed value and the certificate
    if not certificate:
        return (-1, [], h_UB)

    # When using far-apart pairs, the loops may end before improving the
    # upper-bound
    return (h, certificate, h_UB if GOTO_RETURN else h)


def hyperbolicity(G,
                  algorithm='BCCM',
                  approximation_factor=None,
                  additive_gap=None,
                  verbose=False):
    r"""
    Return the hyperbolicity of the graph or an approximation of this value.

    The hyperbolicity of a graph has been defined by Gromov [Gro1987]_ as
    follows: Let `a, b, c, d` be vertices of the graph, let `S_1 = dist(a, b) +
    dist(b, c)`, `S_2 = dist(a, c) + dist(b, d)`, and `S_3 = dist(a, d) +
    dist(b, c)`, and let `M_1` and `M_2` be the two largest values among `S_1`,
    `S_2`, and `S_3`. We have `hyp(a, b, c, d) = |M_1 - M_2|`, and the
    hyperbolicity of the graph is the maximum over all possible 4-tuples `(a,b,
    c,d)` divided by 2. The worst case time complexity is in `O( n^4 )`.

    See the documentation of :mod:`sage.graphs.hyperbolicity` for more
    information.

    INPUT:

    - ``G`` -- a connected Graph

    - ``algorithm`` -- (default: ``'BCCM'``) specifies the algorithm to use
      among:

          - ``'basic'`` is an exhaustive algorithm considering all possible
            4-tuples and so have time complexity in `O(n^4)`.

          - ``'CCL'`` is an exact algorithm proposed in [CCL2015]_. It considers
            the 4-tuples in an ordering allowing to cut the search space as soon
            as a new lower bound is found (see the module's documentation). This
            algorithm can be turned into a approximation algorithm.

          - ``'CCL+FA'`` or ``'CCL+'`` uses the notion of far-apart pairs as
            proposed in [Sot2011]_ to significantly reduce the overall
            computation time of the ``'CCL'`` algorithm.

          - ``'BCCM'`` is an exact algorithm proposed in [BCCM2015]_. It
            improves ``'CCL+FA'`` by cutting several 4-tuples (for more
            information, see the module's documentation).

          - ``'dom'`` is an approximation with additive constant four. It
            computes the hyperbolicity of the vertices of a dominating set of
            the graph. This is sometimes slower than ``'CCL'`` and sometimes
            faster. Try it to know if it is interesting for you.
            The ``additive_gap`` and ``approximation_factor`` parameters cannot
            be used in combination with this method and so are ignored.

    - ``approximation_factor`` -- (default: ``None``) when the approximation factor
      is set to some value (larger than 1.0), the function stop computations as
      soon as the ratio between the upper bound and the best found solution is
      less than the approximation factor. When the approximation factor is 1.0,
      the problem is solved optimally. This parameter is used only when the
      chosen algorithm is ``'CCL'``, ``'CCL+FA'``, or ``'BCCM'``.

    - ``additive_gap`` -- (default: ``None``) when set to a positive number, the
      function stop computations as soon as the difference between the upper
      bound and the best found solution is less than additive gap. When the gap
      is 0.0, the problem is solved optimally. This parameter is used only when
      the chosen algorithm is ``'CCL'`` or ``'CCL+FA'``, or ``'BCCM'``.

    - ``verbose`` -- boolean (default: ``False``); set to ``True`` to display
      some information during execution: new upper and lower bounds, etc.

    OUTPUT:

    This function returns the tuple ( delta, certificate, delta_UB ), where:

    - ``delta`` -- the hyperbolicity of the graph (half-integer value)

    - ``certificate`` -- is the list of the 4 vertices for which the maximum
      value has been computed, and so the hyperbolicity of the graph

    - ``delta_UB`` -- is an upper bound for ``delta``. When ``delta ==
      delta_UB``, the returned solution is optimal. Otherwise, the approximation
      factor if ``delta_UB/delta``.

    EXAMPLES:

    Hyperbolicity of a `3\times 3` grid::

        sage: from sage.graphs.hyperbolicity import hyperbolicity
        sage: G = graphs.Grid2dGraph(3, 3)
        sage: L,C,U = hyperbolicity(G, algorithm='BCCM'); L,sorted(C),U
        (2, [(0, 0), (0, 2), (2, 0), (2, 2)], 2)
        sage: L,C,U = hyperbolicity(G, algorithm='CCL'); L,sorted(C),U
        (2, [(0, 0), (0, 2), (2, 0), (2, 2)], 2)
        sage: L,C,U = hyperbolicity(G, algorithm='basic'); L,sorted(C),U
        (2, [(0, 0), (0, 2), (2, 0), (2, 2)], 2)

    Hyperbolicity of a PetersenGraph::

        sage: from sage.graphs.hyperbolicity import hyperbolicity
        sage: G = graphs.PetersenGraph()
        sage: L,C,U = hyperbolicity(G, algorithm='BCCM'); L,sorted(C),U
        (1/2, [6, 7, 8, 9], 1/2)
        sage: L,C,U = hyperbolicity(G, algorithm='CCL'); L,sorted(C),U
        (1/2, [0, 1, 2, 3], 1/2)
        sage: L,C,U = hyperbolicity(G, algorithm='CCL+'); L,sorted(C),U
        (1/2, [0, 1, 2, 3], 1/2)
        sage: L,C,U = hyperbolicity(G, algorithm='CCL+FA'); L,sorted(C),U
        (1/2, [0, 1, 2, 3], 1/2)
        sage: L,C,U = hyperbolicity(G, algorithm='basic'); L,sorted(C),U
        (1/2, [0, 1, 2, 3], 1/2)
        sage: L,C,U = hyperbolicity(G, algorithm='dom'); L,U
        (0, 1)
        sage: sorted(C)  # random
        [0, 1, 2, 6]

    Asking for an approximation in a grid graph::

        sage: from sage.graphs.hyperbolicity import hyperbolicity
        sage: G = graphs.Grid2dGraph(2, 10)
        sage: L,C,U = hyperbolicity(G, algorithm='CCL', approximation_factor=1.5); L,U
        (1, 3/2)
        sage: L,C,U = hyperbolicity(G, algorithm='CCL+', approximation_factor=1.5); L,U
        (1, 1)
        sage: L,C,U = hyperbolicity(G, algorithm='CCL', approximation_factor=4); L,U
        (1, 4)
        sage: L,C,U = hyperbolicity(G, algorithm='CCL', additive_gap=2); L,U
        (1, 3)
        sage: L,C,U = hyperbolicity(G, algorithm='dom'); L,U
        (1, 5)

    Asking for an approximation in a cycle graph::

        sage: from sage.graphs.hyperbolicity import hyperbolicity
        sage: G = graphs.CycleGraph(10)
        sage: L,C,U = hyperbolicity(G, algorithm='CCL', approximation_factor=1.5); L,U
        (2, 5/2)
        sage: L,C,U = hyperbolicity(G, algorithm='CCL+FA', approximation_factor=1.5); L,U
        (2, 5/2)
        sage: L,C,U = hyperbolicity(G, algorithm='CCL+FA', additive_gap=1); L,U
        (2, 5/2)

    Comparison of results::

        sage: from sage.graphs.hyperbolicity import hyperbolicity
        sage: for i in range(10):               # long time                             # needs networkx
        ....:     G = graphs.RandomBarabasiAlbert(100,2)
        ....:     d1,_,_ = hyperbolicity(G, algorithm='basic')
        ....:     d2,_,_ = hyperbolicity(G, algorithm='CCL')
        ....:     d3,_,_ = hyperbolicity(G, algorithm='CCL+')
        ....:     d4,_,_ = hyperbolicity(G, algorithm='CCL+FA')
        ....:     d5,_,_ = hyperbolicity(G, algorithm='BCCM')
        ....:     l3,_,u3 = hyperbolicity(G, approximation_factor=2)
        ....:     if (not d1==d2==d3==d4==d5) or l3>d1 or u3<d1:
        ....:        print("That's not good!")

        sage: from sage.graphs.hyperbolicity import hyperbolicity
        sage: import random
        sage: random.seed()
        sage: for i in range(10):               # long time                             # needs networkx
        ....:     n = random.randint(2, 20)
        ....:     m = random.randint(0, n*(n-1) / 2)
        ....:     G = graphs.RandomGNM(n, m)
        ....:     for cc in G.connected_components_subgraphs():
        ....:         d1,_,_ = hyperbolicity(cc, algorithm='basic')
        ....:         d2,_,_ = hyperbolicity(cc, algorithm='CCL')
        ....:         d3,_,_ = hyperbolicity(cc, algorithm='CCL+')
        ....:         d4,_,_ = hyperbolicity(cc, algorithm='CCL+FA')
        ....:         d5,_,_ = hyperbolicity(cc, algorithm='BCCM')
        ....:         l3,_,u3 = hyperbolicity(cc, approximation_factor=2)
        ....:         if (not d1==d2==d3==d4==d5) or l3>d1 or u3<d1:
        ....:             print("Error in graph ", cc.edges(sort=True))

    The hyperbolicity of a graph is the maximum value over all its biconnected
    components::

        sage: from sage.graphs.hyperbolicity import hyperbolicity
        sage: G = graphs.PetersenGraph() * 2
        sage: G.add_edge(0, 11)
        sage: L,C,U = hyperbolicity(G); L,sorted(C),U
        (1/2, [6, 7, 8, 9], 1/2)

    TESTS:

    Giving anything else than a Graph::

        sage: from sage.graphs.hyperbolicity import hyperbolicity
        sage: hyperbolicity([])
        Traceback (most recent call last):
        ...
        ValueError: the input parameter must be a Graph

    Giving a non connected graph::

        sage: from sage.graphs.hyperbolicity import hyperbolicity
        sage: G = Graph([(0,1),(2,3)])
        sage: hyperbolicity(G)
        Traceback (most recent call last):
        ...
        ValueError: the input Graph must be connected

    Giving wrong approximation factor::

        sage: from sage.graphs.hyperbolicity import hyperbolicity
        sage: G = graphs.PetersenGraph()
        sage: hyperbolicity(G, algorithm='CCL', approximation_factor=0.1)
        Traceback (most recent call last):
        ...
        ValueError: the approximation factor must be >= 1.0

    Giving negative additive gap::

        sage: from sage.graphs.hyperbolicity import hyperbolicity
        sage: G = Graph()
        sage: hyperbolicity(G, algorithm='CCL', additive_gap=-1)
        Traceback (most recent call last):
        ...
        ValueError: the additive gap must be a real positive number

    Asking for an unknown algorithm::

        sage: from sage.graphs.hyperbolicity import hyperbolicity
        sage: G = Graph()
        sage: hyperbolicity(G, algorithm='tip top')
        Traceback (most recent call last):
        ...
        ValueError: algorithm 'tip top' not yet implemented, please contribute
    """

    # Abbreviations for algorithms are expanded.
    if algorithm == "CCL+":
        algorithm = "CCL+FA"

    from sage.graphs.graph import Graph
    if not isinstance(G, Graph):
        raise ValueError("the input parameter must be a Graph")
    if algorithm not in ['basic', 'CCL', 'CCL+FA', 'BCCM', 'dom']:
        raise ValueError("algorithm '%s' not yet implemented, please contribute" % (algorithm))
    if approximation_factor is None:
        approximation_factor = 1.0
    elif approximation_factor == 1.0:
        pass
    elif algorithm in ['CCL', 'CCL+FA', 'BCCM']:
        from sage.rings.real_mpfr import RR

        if approximation_factor not in RR or approximation_factor < 1.0:
            raise ValueError("the approximation factor must be >= 1.0")
    else:
        raise ValueError("the approximation_factor is ignored when using"
                         "the '%s' algorithm" % (algorithm))
    if additive_gap is None:
        additive_gap = 0.0
    elif additive_gap == 0.0:
        pass
    elif algorithm in ['CCL', 'CCL+FA', 'BCCM']:
        from sage.rings.real_mpfr import RR

        if additive_gap not in RR or additive_gap < 0.0:
            raise ValueError("the additive gap must be a real positive number")
    else:
        raise ValueError("the additive_gap is ignored when using the '%s' algorithm." % (algorithm))

    # The hyperbolicity is defined on connected graphs
    if not G.is_connected():
        raise ValueError("the input Graph must be connected")

    # The hyperbolicity of some classes of graphs is known. If it is easy and
    # fast to test that a graph belongs to one of these classes, we do it.
    if G.num_verts() <= 3:
        # The hyperbolicity of a graph with 3 vertices is 0.
        # The certificate is the set of vertices.
        return 0, list(G), 0

    elif G.num_verts() == G.num_edges() + 1:
        # G is a tree
        # Any set of 4 vertices is a valid certificate
        return 0, list(G)[:4], 0

    elif G.is_clique():
        # Any set of 4 vertices is a valid certificate
        return 0, list(G)[:4], 0

    cdef int i, j, D
    cdef list certificate = []
    cdef list certif

    cdef int N = G.num_verts()
    hyp = 0
    hyp_UB = 0

    #
    # The hyperbolicity of a graph is the maximum over its 2-connected
    # components.
    #
    B, _ = G.blocks_and_cut_vertices()
    if len(B) > 1:

        if verbose:
            # we compute the distribution of size of the blocks
            L = [len(V) for V in B]
            print("Graph with %d blocks" % (len(B)))
            print("Blocks size distribution:", {x: L.count(x) for x in L})

        for V in B:

            # The hyperbolicity of a graph with 3 vertices is 0, and a graph
            # cannot have hyperbolicity larger than N/2. So we consider only
            # larger 2-connected subgraphs.
            if len(V) > max(3, 2 * hyp):

                hh, certif, hh_UB = hyperbolicity(_my_subgraph(G, V), algorithm=algorithm,
                                                  approximation_factor=approximation_factor,
                                                  additive_gap=additive_gap, verbose=verbose)

                # We test if the new computed value improves upon previous value.
                if hh > hyp or (hh == hyp and not certificate):
                    hyp = hh
                    certificate = certif

                # We update independently the upper bound for cases in which we
                # are asking for an approximation.
                hyp_UB = max(hyp_UB, hh_UB)

        # Last, we return the computed value and the certificate
        return hyp, certificate, hyp_UB

    #
    # Now the graph is 2-connected, has at least 4 vertices and is not a clique.
    #

    cdef unsigned short* _distances_
    cdef unsigned short** distances
    cdef unsigned short* _far_apart_pairs_
    cdef unsigned short** far_apart_pairs
    cdef list int_to_vertex = list(G)

    # We compute the distances and store the results in a 2D array
    distances = <unsigned short **>check_allocarray(N, sizeof(unsigned short *))
    if not distances:
        raise MemoryError("Unable to allocate array 'distances'.")

    if algorithm == 'CCL+FA' or algorithm == 'BCCM':
        _distances_ = <unsigned short *>check_allocarray(N * N, sizeof(unsigned short))
        _far_apart_pairs_ = <unsigned short *>check_allocarray(N * N, sizeof(unsigned short))
        far_apart_pairs = <unsigned short **>check_allocarray(N, sizeof(unsigned short *))

        distances_and_far_apart_pairs(G, _distances_, _far_apart_pairs_, int_to_vertex)

        for i in range(N):
            far_apart_pairs[i] = _far_apart_pairs_ + i*N

    else:
        _distances_ = c_distances_all_pairs(G, vertex_list=int_to_vertex)
        _far_apart_pairs_ = NULL
        far_apart_pairs = NULL

    D = 0
    for i in range(N):
        distances[i] = _distances_ + i * N
        for j in range(i + 1, N):
            if distances[i][j] > D:
                D = distances[i][j]

    # We call the cython function for computing the hyperbolicity with the
    # required parameters.
    if algorithm in ['CCL', 'CCL+FA']:
        sig_on()
        hyp, certif, hyp_UB = hyperbolicity_CCL(N, distances, far_apart_pairs, D, hyp,
                                                approximation_factor, 2 * additive_gap, verbose)
        sig_off()

    elif algorithm == 'BCCM':
        sig_on()
        hyp, certif, hyp_UB = hyperbolicity_BCCM(N, distances, far_apart_pairs,
                                                 D, hyp, approximation_factor,
                                                 2 * additive_gap, verbose)
        sig_off()

    elif algorithm == 'dom':
        # Computes a dominating set DOM of G, and computes the hyperbolicity
        # considering only vertices in DOM
        DOM = set(_greedy_dominating_set(G, verbose=verbose))
        # We need at least 4 vertices
        while len(DOM) < 4:
            DOM.add(G.random_vertex())
        # We map the dominating set to [0..N-1]
        v_to_int = {v: i for i, v in enumerate(G.vertex_iterator())}
        DOM_int = set(v_to_int[v] for v in DOM)
        # We set null distances to vertices outside DOM. This way these
        # vertices will not be considered anymore.
        for i in range(N):
            if i not in DOM_int:
                for j in range(N):
                    distances[i][j] = 0
                    distances[j][i] = 0
        sig_on()
        hyp, certif, hyp_UB = hyperbolicity_CCL(N, distances, NULL, D, hyp, 1.0, 0.0, verbose)
        sig_off()
        hyp_UB = min(hyp + 8, D)

    elif algorithm == 'basic':
        sig_on()
        hyp, certif = hyperbolicity_basic_algorithm(N, distances, verbose=verbose)
        sig_off()
        hyp_UB = hyp

    # We now release the memory
    sig_free(distances)
    sig_free(_distances_)
    sig_free(_far_apart_pairs_)
    sig_free(far_apart_pairs)

    # Map the certificate 'certif' with the corresponding vertices in the graph
    certificate = [int_to_vertex[i] for i in certif]

    # Last, we return the computed value and the certificate
    return ZZ(hyp)/2, certificate, ZZ(hyp_UB)/2


######################################################################
# Distribution of the hyperbolicity of 4-tuples
######################################################################

cdef dict __hyperbolicity_distribution__(int N, unsigned short** distances):
    """
    Return the distribution of the hyperbolicity of the 4-tuples of the graph.

    The hyperbolicity of a graph has been defined by Gromov [Gro1987]_ as
    follows: Let `a, b, c, d` be vertices of the graph, let `S_1 = dist(a, b) +
    dist(b, c)`, `S_2 = dist(a, c) + dist(b, d)`, and `S_3 = dist(a, d) +
    dist(b, c)`, and let `M_1` and `M_2` be the two largest values among `S_1`,
    `S_2`, and `S_3`. We have `hyp(a, b, c, d) = |M_1 - M_2|`, and the
    hyperbolicity of the graph is the maximum over all possible 4-tuples `(a, b,
    c, d)` divided by 2.

    The computation of the hyperbolicity of each 4-tuple, and so the
    hyperbolicity distribution, takes time in `O( n^4 )`.

    We use ``unsigned long int`` on 64 bits, so ``uint64_t``, to count the
    number of 4-tuples of given hyperbolicity. So we cannot exceed `2^64-1`.
    This value should be sufficient for most users.

    INPUT:

    - ``N`` -- number of vertices of the graph (and side of the matrix)

    - ``distances`` -- matrix of distances in the graph

    OUTPUT:

    - ``hdict`` -- dictionary such that hdict[i] is the number of 4-tuples of
      hyperbolicity i among the considered 4-tuples
    """
    # We initialize the table of hyperbolicity. We use an array of unsigned long
    # int instead of a dictionary since it is much faster.
    cdef int i

    cdef uint64_t* hdistr = <uint64_t*>check_calloc(N + 1, sizeof(uint64_t))
    if not hdistr:
        raise MemoryError

    # We now compute the hyperbolicity of each 4-tuple
    cdef int a, b, c, d
    for a in range(N - 3):
        for b in range(a + 1, N - 2):
            for c in range(b + 1, N - 1):
                for d in range(c + 1, N):
                    hdistr[__hyp__(distances, a, b, c, d)] += 1

    # We prepare the dictionary of hyperbolicity distribution to return
    Nchoose4 = binomial(N, 4)
    cdef dict hdict = {ZZ(i)/2: (ZZ(hdistr[i]) / Nchoose4) for i in range(N + 1) if hdistr[i] > 0}

    sig_free(hdistr)

    return hdict


# We use this trick since it is way faster than using the sage randint function.
cdef extern from "stdlib.h":
    long c_libc_random "random"()
    void c_libc_srandom "srandom"(unsigned int seed)


cdef dict __hyperbolicity_sampling__(int N, unsigned short** distances, uint64_t sampling_size):
    """
    Return a sampling of the hyperbolicity distribution of the graph.

    The hyperbolicity of a graph has been defined by Gromov [Gro1987]_ as
    follows: Let `a, b, c, d` be vertices of the graph, let `S_1 = dist(a, b) +
    dist(b, c)`, `S_2 = dist(a, c) + dist(b, d)`, and `S_3 = dist(a, d) +
    dist(b, c)`, and let `M_1` and `M_2` be the two largest values among `S_1`,
    `S_2`, and `S_3`. We have `hyp(a, b, c, d) = |M_1 - M_2|`, and the
    hyperbolicity of the graph is the maximum over all possible 4-tuples `(a, b,
    c, d)` divided by 2.

    We use ``unsigned long int`` on 64 bits, so ``uint64_t``, to count the
    number of 4-tuples of given hyperbolicity. So we cannot exceed `2^64-1`.
    This value should be sufficient for most users.

    INPUT:

    - ``N`` -- number of vertices of the graph (and side of the matrix)

    - ``distances`` -- matrix of distances in the graph

    - ``sampling_size`` -- number of 4-tuples considered. Default value is 1000

    OUTPUT:

    - ``hdict`` -- dictionary such that hdict[i] is the number of 4-tuples of
                hyperbolicity i among the considered 4-tuples
    """
    cdef int i, a, b, c, d
    cdef uint64_t j

    if N < 4:
        raise ValueError("N must be at least 4")

    # We initialize the table of hyperbolicity. We use an array of unsigned long
    # int instead of a dictionary since it is much faster.
    cdef uint64_t* hdistr = <uint64_t*>check_calloc(N + 1, sizeof(uint64_t))
    if not hdistr:
        raise MemoryError

    # We now compute the hyperbolicity of each quadruple
    for j in range(sampling_size):
        a = c_libc_random() % N
        b = c_libc_random() % N
        c = c_libc_random() % N
        d = c_libc_random() % N
        while a == b:
            b = c_libc_random() % N
        while a == c or b == c:
            c = c_libc_random() % N
        while a == d or b == d or c == d:
            d = c_libc_random() % N

        hdistr[__hyp__(distances, a, b, c, d)] += 1

    # We prepare the dictionary of hyperbolicity distribution from sampling
    cdef dict hdict = {ZZ(i)/2: ZZ(hdistr[i])/ZZ(sampling_size) for i in range(N + 1) if hdistr[i] > 0}

    sig_free(hdistr)

    return hdict


def hyperbolicity_distribution(G, algorithm='sampling', sampling_size=10**6):
    r"""
    Return the hyperbolicity distribution of the graph or a sampling of it.

    The hyperbolicity of a graph has been defined by Gromov [Gro1987]_ as
    follows: Let `a, b, c, d` be vertices of the graph, let `S_1 = dist(a, b) +
    dist(b, c)`, `S_2 = dist(a, c) + dist(b, d)`, and `S_3 = dist(a, d) +
    dist(b, c)`, and let `M_1` and `M_2` be the two largest values among `S_1`,
    `S_2`, and `S_3`. We have `hyp(a, b, c, d) = |M_1 - M_2|`, and the
    hyperbolicity of the graph is the maximum over all possible 4-tuples `(a, b,
    c, d)` divided by 2.

    The computation of the hyperbolicity of each 4-tuple, and so the
    hyperbolicity distribution, takes time in `O( n^4 )`.

    INPUT:

    - ``G`` -- a Graph

    - ``algorithm`` -- (default: ``'sampling'``) when algorithm is 'sampling', it
      returns the distribution of the hyperbolicity over a sample of
      ``sampling_size`` 4-tuples. When algorithm is 'exact', it computes the
      distribution of the hyperbolicity over all 4-tuples. Be aware that the
      computation time can be HUGE.

    - ``sampling_size`` -- (default: `10^6`) number of 4-tuples considered in
      the sampling. Used only when ``algorithm == 'sampling'``

    OUTPUT:

    - ``hdict`` -- dictionary such that hdict[i] is the number of 4-tuples of
      hyperbolicity i

    EXAMPLES:

    Exact hyperbolicity distribution of the Petersen Graph::

        sage: from sage.graphs.hyperbolicity import hyperbolicity_distribution
        sage: G = graphs.PetersenGraph()
        sage: hyperbolicity_distribution(G,algorithm='exact')
        {0: 3/7, 1/2: 4/7}

    Exact hyperbolicity distribution of a `3\times 3` grid::

        sage: from sage.graphs.hyperbolicity import hyperbolicity_distribution
        sage: G = graphs.GridGraph([3,3])
        sage: hyperbolicity_distribution(G,algorithm='exact')
        {0: 11/18, 1: 8/21, 2: 1/126}

    TESTS:

    Giving anything else than a Graph::

        sage: from sage.graphs.hyperbolicity import hyperbolicity_distribution
        sage: hyperbolicity_distribution([])
        Traceback (most recent call last):
        ...
        ValueError: the input parameter must be a Graph

    Giving a non connected graph::

        sage: from sage.graphs.hyperbolicity import hyperbolicity_distribution
        sage: G = Graph([(0,1),(2,3)])
        sage: hyperbolicity_distribution(G)
        Traceback (most recent call last):
        ...
        ValueError: the input Graph must be connected
    """
    from sage.graphs.graph import Graph
    if not isinstance(G, Graph):
        raise ValueError("the input parameter must be a Graph")
    # The hyperbolicity is defined on connected graphs
    if not G.is_connected():
        raise ValueError("the input Graph must be connected")

    # The hyperbolicity distribution of some classes of graphs is known. If it
    # is easy and fast to test that a graph belongs to one of these classes, we
    # do it.
    if (G.num_verts() == G.num_edges() + 1) or G.is_clique():
        return {0: sampling_size if algorithm=='sampling' else binomial(G.num_verts(), 4)}

    cdef int N = G.num_verts()
    cdef int i
    cdef unsigned short** distances
    cdef unsigned short* _distances_
    cdef dict hdict

    # We compute the all pairs shortest path and store the result in a 2D array
    # for faster access.
    _distances_ = c_distances_all_pairs(G, vertex_list=list(G))
    distances = <unsigned short**>check_allocarray(N, sizeof(unsigned short*))
    if not distances:
        sig_free(_distances_)
        raise MemoryError

    for i in range(N):
        distances[i] = _distances_ + i * N

    if algorithm == 'exact':
        hdict = __hyperbolicity_distribution__(N, distances)
    elif algorithm == 'sampling':
        hdict = __hyperbolicity_sampling__(N, distances, sampling_size)
    else:
        raise ValueError("algorithm '%s' not yet implemented, please contribute" % (algorithm))

    # We release memory
    sig_free(distances)
    sig_free(_distances_)

    return hdict
