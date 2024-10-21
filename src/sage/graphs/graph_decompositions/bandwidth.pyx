# sage_setup: distribution = sagemath-graphs
# cython: binding=True
r"""
Bandwidth of undirected graphs

Definition
----------

The bandwidth `bw(M)` of a matrix `M` is the smallest integer `k` such that all
nonzero entries of `M` are at distance `k` from the diagonal. The bandwidth
`bw(G)` of an undirected graph `G` is the minimum bandwidth of the adjacency
matrix of `G`, over all possible relabellings of its vertices.

**Path spanner:** alternatively, the bandwidth measures how tightly a path
represents the distance of a graph `G`. Indeed, if the vertices of `G` can be
ordered as `v_1,...,v_n` in such a way that `k \times d_G(v_i,v_j) \geq |i-j|` then
`bw(G)\leq k`.

    **Proof:** for all `v_i \sim v_j` (i.e. `d_G(v_i,v_j)=1`), the constraint
    ensures that `k\geq |i-j|`, meaning that adjacent vertices are at distance
    at most `k` in the path ordering. That alone is sufficient to ensure that
    `bw(G)\leq k`.

    As a byproduct, we obtain that `k \times d_G(v_i,v_j) \geq |i-j|` in
    general: let `v_{s_0},...,v_{s_i}` be the vertices of a shortest
    `(v_i,v_j)`-path. We have:

    .. MATH::

        k \times d_G(v_i,v_j) &= k\times d_G(v_i,v_{s_0}) + k\times d_G(v_{s_0},v_{s_1}) + \cdots + k\times d_G(v_{s_{i-1}},v_{s_i}) + k\times d_G(v_{s_i},v_j)\\
                              &\geq |v_i-v_{s_0}| + |v_{s_0}-v_{s_1}| + \cdots + |v_{s_{i-1}}-v_{s_i}| + |v_{s_i}-v_j|\\
                              &\geq |v_i-v_j|\\

Satisfiability of a partial assignment
--------------------------------------

Let us suppose that the first `i` vertices `v_1,...,v_i` of `G` have already
been assigned positions `p_1,...,p_i` in an ordering of `V(G)` of bandwidth
`\leq k`. Where can `v_{i+1}` appear ?

Because of the previous definition, `p_{i+1}` must be at distance at most
`k\times d_G(v_1,v_{i+1})` from `p_1`, and in general at distance at most
`k\times d_G(v_j,v_{i+1})` from `p_j`. Each range is an interval of
`\{1,...,n\}\backslash \{p_1,...,p_i\}`, and because the intersection of two
intervals is again an interval we deduce that in order to satisfy all these
constraints simultaneously `p_j` must belong to an interval defined from this
partial assignment.

Applying this rule to all non-assigned vertices, we deduce that each of them
must be assigned to a given interval of `\{1,...,n\}`. Note that this can also
be extended to the already assigned vertices, by saying that `v_j` with `j<i`
must be assigned within the interval `[p_j,p_j]`.

This problem is not always satisfiable, e.g. 5 vertices cannot all be assigned
to the elements of `[10,13]`. This is a matching problem which, because all
admissible sets are intervals, can be solved quickly.

Solving the matching problem
----------------------------

Let `n` points `v_1,...,v_n` be given, along with two functions `m,M:[n]\mapsto
[n]`. Is there an ordering `p_1,...,p_n` of them such that `m(v_i) \leq p_i \leq
M(v_i)` ? This is equivalent to Hall's bipartite matching theorem, and can in
this specific case be solved by the following algorithm:

- Consider all vertices `v` sorted increasingly according to `M(v)`

- For each of them, assign to `v` the smallest position in `[m(v),M(v)]` which
  has not been assigned yet. If there is none, the assignment problem is not
  satisfiable.

Note that the latest operation can be performed with very few bitset operations
(provided that `n<64`).

The algorithm
-------------

This section contains totally subjective choices, that may be changed in the
hope to get better performances.

- Try to find a satisfiable ordering by filling positions, one after the other
  (and not by trying to find each vertex' position)

- Fill the positions in this order: `0,n-1,1,n-2,3,n-3, ...`

.. NOTE::

    There is some symmetry to break as the reverse of a satisfiable ordering is
    also a satisfiable ordering.

This module contains the following methods
------------------------------------------

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`bandwidth` | Compute the bandwidth of an undirected graph
    :meth:`~sage.graphs.base.boost_graph.bandwidth_heuristics` | Use Boost heuristics to approximate the bandwidth of the input graph

Functions
---------
"""

# ****************************************************************************
#       Copyright (C) 2015 Nathann Cohen <nathann.cohen@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from libc.stdint cimport uint16_t
from cysignals.signals cimport sig_check
from memory_allocator cimport MemoryAllocator

from sage.graphs.distances_all_pairs cimport all_pairs_shortest_path_BFS

ctypedef uint16_t index_t

ctypedef struct range_t:
    index_t m
    index_t M


def bandwidth(G, k=None):
    r"""
    Compute the bandwidth of an undirected graph.

    For a definition of the bandwidth of a graph, see the documentation of the
    :mod:`~sage.graphs.graph_decompositions.bandwidth` module.

    INPUT:

    - ``G`` -- a graph

    - ``k`` -- integer (default: ``None``); set to an integer value to test
      whether `bw(G)\leq k`, or to ``None`` (default) to compute `bw(G)`

    OUTPUT:

    When `k` is an integer value, the function returns either ``False`` or an
    ordering of cost `\leq k`.

    When `k` is equal to ``None``, the function returns a pair ``(bw,
    ordering)``.

    .. SEEALSO::

        :meth:`sage.graphs.generic_graph.GenericGraph.adjacency_matrix` --
        return the adjacency matrix from an ordering of the vertices.

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.bandwidth import bandwidth
        sage: G = graphs.PetersenGraph()
        sage: bandwidth(G,3)
        False
        sage: bandwidth(G)
        (5, [0, 4, 5, 8, 1, 9, 3, 7, 6, 2])
        sage: G.adjacency_matrix(vertices=[0, 4, 5, 8, 1, 9, 3, 7, 6, 2])               # needs sage.modules
        [0 1 1 0 1 0 0 0 0 0]
        [1 0 0 0 0 1 1 0 0 0]
        [1 0 0 1 0 0 0 1 0 0]
        [0 0 1 0 0 0 1 0 1 0]
        [1 0 0 0 0 0 0 0 1 1]
        [0 1 0 0 0 0 0 1 1 0]
        [0 1 0 1 0 0 0 0 0 1]
        [0 0 1 0 0 1 0 0 0 1]
        [0 0 0 1 1 1 0 0 0 0]
        [0 0 0 0 1 0 1 1 0 0]
        sage: G = graphs.ChvatalGraph()
        sage: bandwidth(G)
        (6, [0, 5, 9, 4, 10, 1, 6, 11, 3, 8, 7, 2])
        sage: G.adjacency_matrix(vertices=[0, 5, 9, 4, 10, 1, 6, 11, 3, 8, 7, 2])       # needs sage.modules
        [0 0 1 1 0 1 1 0 0 0 0 0]
        [0 0 0 1 1 1 0 1 0 0 0 0]
        [1 0 0 0 1 0 0 1 1 0 0 0]
        [1 1 0 0 0 0 0 0 1 1 0 0]
        [0 1 1 0 0 0 1 0 0 1 0 0]
        [1 1 0 0 0 0 0 0 0 0 1 1]
        [1 0 0 0 1 0 0 1 0 0 0 1]
        [0 1 1 0 0 0 1 0 0 0 1 0]
        [0 0 1 1 0 0 0 0 0 0 1 1]
        [0 0 0 1 1 0 0 0 0 0 1 1]
        [0 0 0 0 0 1 0 1 1 1 0 0]
        [0 0 0 0 0 1 1 0 1 1 0 0]

    TESTS::

        sage: bandwidth(2*graphs.PetersenGraph())
        (5, [0, 4, 5, 8, 1, 9, 3, 7, 6, 2, 10, 14, 11, 19, 15, 16, 13, 12, 17, 18])
        sage: bandwidth(Graph())
        (0, [])
        sage: bandwidth(Graph(1))
        (0, [0])
        sage: bandwidth(Graph(3))
        (0, [0, 1, 2])

    Directed/weighted graphs::

        sage: bandwidth(digraphs.Circuit(5))
        Traceback (most recent call last):
        ...
        ValueError: this method only works on undirected graphs
        sage: bandwidth(Graph(graphs.PetersenGraph(), weighted=True))
        Traceback (most recent call last):
        ...
        ValueError: this method only works on unweighted graphs
    """
    if G.is_directed():
        raise ValueError("this method only works on undirected graphs")
    if G.weighted():
        raise ValueError("this method only works on unweighted graphs")
    # Trivial cases
    if G.order() <= 1:
        if k is None:
            return 0, list(G)
        else:
            return list(G)

    if not G.is_connected():
        max_k = 0 if k is None else k
        order = []
        for GG in G.connected_components_subgraphs():
            ans = bandwidth(GG, k=k)
            if not ans:
                return False
            if k is None:
                max_k = max(max_k, ans[0])
                ans = ans[1:]
            order.extend(ans[0])
        return (max_k, order)

    # All that this function does is allocate/free the memory for function
    # bandwidth_C

    cdef int n = G.order()
    # Must be the same order than in all_pairs_shortest_path_BFS
    cdef list int_to_vertex = list(G)

    cdef MemoryAllocator mem = MemoryAllocator()

    cdef unsigned short ** d = <unsigned short **> mem.allocarray(n, sizeof(unsigned short *))
    cdef unsigned short * distances = <unsigned short *> mem.allocarray(n*n, sizeof(unsigned short))
    cdef index_t * current = <index_t *> mem.allocarray(n, sizeof(index_t))
    cdef index_t * ordering = <index_t *> mem.allocarray(n, sizeof(index_t))
    cdef index_t * left_to_order = <index_t *> mem.allocarray(n, sizeof(index_t))
    cdef index_t * index_array_tmp = <index_t *> mem.allocarray(n, sizeof(index_t))
    cdef range_t * range_arrays = <range_t *> mem.allocarray(n*n, sizeof(range_t))
    cdef range_t ** ith_range_array = <range_t **> mem.allocarray(n, sizeof(range_t *))
    cdef range_t * range_array_tmp = <range_t *> mem.allocarray(n, sizeof(range_t))

    cdef int i, kk
    # compute the distance matrix
    all_pairs_shortest_path_BFS(G, NULL, distances, NULL, vertex_list=int_to_vertex)

    # fill d so that d[i][j] works
    for i in range(n):
        d[i] = distances + i * n

    # ith_range_array
    for i in range(n):
        ith_range_array[i] = range_arrays + i * n

    # initialize left_to_order
    for i in range(n):
        left_to_order[i] = i

    if k is None:
        for kk in range((n - 1) // G.diameter(), n):
            if bandwidth_C(n, kk, d, current, ordering, left_to_order, index_array_tmp, ith_range_array, range_array_tmp):
                ans = True
                break
    else:
        ans = bool(bandwidth_C(n, k, d, current, ordering, left_to_order, index_array_tmp, ith_range_array, range_array_tmp))

    if ans:
        order = [int_to_vertex[ordering[i]] for i in range(n)]

    if ans:
        ans = (kk, order) if k is None else order

    return ans


cdef bint bandwidth_C(int n, int k,
                      unsigned short ** d,
                      index_t * current,           # choice of vertex for the current position
                      index_t * ordering,          # the actual ordering of vertices
                      index_t * left_to_order,     # begins with the assigned vertices, ends with the others
                      index_t * index_array_tmp,   # tmp space
                      range_t ** ith_range_array,  # array of ranges, for every step of the algorithm
                      range_t * range_array_tmp) noexcept:  # tmp space

    cdef int i, v
    cdef int pi  # the position for which a vertex is being chosen
    cdef int vi  # the vertex being tested at position pi
    cdef int radius
    current[0] = -1

    # At first any vertex can be anywhere
    for v in range(n):
        ith_range_array[0][v].m = 0
        ith_range_array[0][v].M = n - 1

    i = 0
    while True:
        sig_check()

        # There are (n-i) choices for vertex i, as i-1 have already been
        # determined. Thus, i <= current[i] < n.
        current[i] += 1

        # All choices for this position i have been tested. We must change our
        # (i-1)th choice:
        if current[i] == n:
            if i == 0:
                return 0
            i -= 1
            left_to_order[i], left_to_order[current[i]] = left_to_order[current[i]], left_to_order[i]
            continue

        # The position of the ith vertex. p0 = 0, p1 = n-1, p2 = 1, p3 = n-2, ...
        pi = (n - 1 - i // 2) if i % 2 else (i // 2)

        # The ith vertex
        vi = left_to_order[current[i]]

        # If pi is not an admissible position for pi:
        if ith_range_array[i][vi].m > pi or ith_range_array[i][vi].M < pi:
            continue

        # As the choice is admissible, we update left_to_order so that
        # left_to_order[i] = vi.
        left_to_order[i], left_to_order[current[i]] = left_to_order[current[i]], left_to_order[i]

        # vi is at position pi in the final ordering.
        ordering[pi] = vi

        # If we found the position of the nth vertex, we are done.
        if i == n - 1:
            return 1

        # As vertex vi has been assigned position pi, we use that information to
        # update the intervals of admissible positions of all other vertices.
        #
        # \forall v, k*d[v][vi] >= |p_v-p_{vi}| (see module documentation)
        for v in range(n):
            radius = k * d[v][vi]
            ith_range_array[i + 1][v].m = max(<int> ith_range_array[i][v].m, pi - radius)
            ith_range_array[i + 1][v].M = min(<int> ith_range_array[i][v].M, pi + radius)

        # Check the feasibility of a matching with the updated intervals of
        # admissible positions (see module doc).
        #
        # If it is possible we explore deeper, otherwise we undo the changes as
        # pi is not a good position for vi after all.
        if is_matching_feasible(n, ith_range_array[i + 1], range_array_tmp, index_array_tmp):
            i += 1
            current[i] = i - 1
        else:
            # swap back
            left_to_order[i], left_to_order[current[i]] = left_to_order[current[i]], left_to_order[i]

cdef bint is_matching_feasible(int n, range_t * range_array, range_t * range_array_tmp, index_t * index_array_tmp) noexcept:
    r"""
    Test if the matching is feasible

    INPUT:

    - ``n`` -- integer; number of points

    - ``range_array`` -- associates to every point an interval in which the
      point must be given a position

    - ``range_array_tmp`` -- temporary spaces with the same characteristics as
      ``range_array``

    - ``index_array_tmp`` -- temporary space to associate an integer to every
      point

    OUTPUT:

    The function must return a boolean, and does not change the content of
    ``range_array``.
    """
    # Heuristic: check if some vertex has an empty range, that's an easy 'no'.
    cdef int v, M, m, j
    for v in range(n):
        if range_array[v].M < range_array[v].m:
            return 0
        index_array_tmp[v] = 0

    # Sort the guys according to increasing value of M in O(n).

    # Step 1: count the occurrences of each M
    for v in range(n):
        index_array_tmp[range_array[v].M] += 1

    # Step 2: sorted table
    for v in range(1, n):
        index_array_tmp[v] += index_array_tmp[v - 1]

    for v in range(n):
        M = range_array[v].M
        m = range_array[v].m
        index_array_tmp[M] -= 1
        range_array_tmp[index_array_tmp[M]].M = M
        range_array_tmp[index_array_tmp[M]].m = m

    # Satisfiability. We use index_array_tmp as a bitset, and mark every
    # assigned position.
    for v in range(n):
        index_array_tmp[v] = 0

    for v in range(n):
        sig_check()
        for j in range(range_array_tmp[v].m, range_array_tmp[v].M + 1):
            if index_array_tmp[j] == 0:
                index_array_tmp[j] = 1
                break
        else:
            return 0
    return 1
