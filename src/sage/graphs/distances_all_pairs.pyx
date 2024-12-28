# cython: binding=True
r"""
Distances/shortest paths between all pairs of vertices

This module implements a few functions that deal with the computation of
distances or shortest paths between all pairs of vertices.

**Efficiency** : Because these functions involve listing many times the
(out)-neighborhoods of (di)-graphs, it is useful in terms of efficiency to build
a temporary copy of the graph in a data structure that makes it easy to compute
quickly. These functions also work on large volume of data, typically dense
matrices of size `n^2`, and are expected to return corresponding dictionaries of
size `n^2`, where the integers corresponding to the vertices have first been
converted to the vertices' labels. Sadly, this last translating operation turns
out to be the most time-consuming, and for this reason it is also nice to have a
Cython module, and version of these functions that return C arrays, in order to
avoid these operations when they are not necessary.

**Memory cost** : The methods implemented in the current module sometimes need
large amounts of memory to return their result. Storing the distances between
all pairs of vertices in a graph on `1500` vertices as a dictionary of
dictionaries takes around 200MB, while storing the same information as a C array
requires 4MB.


The module's main function
--------------------------

The C function ``all_pairs_shortest_path_BFS`` actually does all the
computations, and all the others (except for ``Floyd_Warshall``) are just
wrapping it. This function begins with copying the graph in a data structure
that makes it fast to query the out-neighbors of a vertex, then starts one
Breadth First Search per vertex of the (di)graph.

**What can this function compute ?**

- The matrix of predecessors.

  This matrix `P` has size `n^2`, and is such that vertex `P[u,v]` is a
  predecessor of `v` on a shortest `uv`-path. Hence, this matrix efficiently
  encodes the information of a shortest `uv`-path for any `u,v\in G` : indeed,
  to go from `u` to `v` you should first find a shortest `uP[u,v]`-path, then
  jump from `P[u,v]` to `v` as it is one of its outneighbors. Apply recursively
  and find out what the whole path is !.

- The matrix of distances.

  This matrix has size `n^2` and associates to any `uv` the distance
  from `u` to `v`.

- The vector of eccentricities.

  This vector of size `n` encodes for each vertex `v` the distance to vertex
  which is furthest from `v` in the graph. In particular, the diameter of the
  graph is the maximum of these values.

**What does it take as input ?**

- ``gg`` a (Di)Graph.

- ``unsigned short * predecessors`` -- a pointer toward an array of size
  `n^2\cdot\text{sizeof(unsigned short)}`. Set to ``NULL`` if you do not want to
  compute the predecessors.

- ``unsigned short * distances`` -- a pointer toward an array of size
  `n^2\cdot\text{sizeof(unsigned short)}`. The computation of the distances is
  necessary for the algorithm, so this value can **not** be set to ``NULL``.

- ``int * eccentricity`` -- a pointer toward an array of size
  `n\cdot\text{sizeof(int)}`. Set to ``NULL`` if you do not want to compute the
  eccentricity.

**Technical details**

- The vertices are encoded as `1, ..., n` as they appear in the ordering of
  ``G.vertices(sort=True)``, unless another ordering is specified by the user.

- Because this function works on matrices whose size is quadratic compared to
  the number of vertices when computing all distances or predecessors, it uses
  short variables to store the vertices' names instead of long ones to divide by
  2 the size in memory. This means that only the diameter/eccentricities can be
  computed on a graph of more than 65536 nodes. For information, the current
  version of the algorithm on a graph with `65536=2^{16}` nodes creates in
  memory `2` tables on `2^{32}` short elements (2bytes each), for a total of
  `2^{33}` bytes or `8` gigabytes. In order to support larger sizes, we would
  have to replace shorts by 32-bits int or 64-bits int, which would then require
  respectively 16GB or 32GB.

- In the C version of these functions, infinite distances are represented with
  ``<unsigned short> -1 = 65535`` for ``unsigned short`` variables, and by
  ``INT32_MAX`` otherwise. These case happens when the input is a disconnected
  graph, or a non-strongly-connected digraph.

- A memory error is raised when data structures allocation failed. This could
  happen with large graphs on computers with low memory space.

.. WARNING::

    The function ``all_pairs_shortest_path_BFS`` has **no reason** to be called
    by the user, even though he would be writing his code in Cython and look for
    efficiency. This module contains wrappers for this function that feed it
    with the good parameters. As the function is inlined, using those wrappers
    actually saves time as it should avoid testing the parameters again and
    again in the main function's body.

AUTHOR:

- Nathann Cohen (2011)
- David Coudert (2014) -- 2sweep, multi-sweep and iFUB for diameter computation

Functions
---------
"""

# ****************************************************************************
#       Copyright (C) 2011 Nathann Cohen <nathann.cohen@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.data_structures.binary_matrix cimport *
from libc.string cimport memset
from libc.stdint cimport uint64_t, UINT64_MAX
from libc.stdint cimport uint32_t, INT32_MAX, UINT32_MAX
from cysignals.memory cimport sig_malloc, sig_calloc, sig_free
from cysignals.signals cimport sig_on, sig_off
from memory_allocator cimport MemoryAllocator

from sage.graphs.base.c_graph cimport CGraphBackend
from sage.graphs.base.c_graph cimport CGraph
from sage.graphs.base.static_sparse_backend cimport StaticSparseCGraph
from sage.graphs.base.static_sparse_backend cimport StaticSparseBackend
from sage.graphs.base.static_sparse_graph cimport (short_digraph,
                                                   init_short_digraph,
                                                   init_reverse,
                                                   free_short_digraph,
                                                   out_degree,
                                                   has_edge,
                                                   simple_BFS)


cdef inline c_all_pairs_shortest_path_BFS(short_digraph sd,
                                          unsigned short* predecessors,
                                          unsigned short* distances,
                                          uint32_t* eccentricity):
    r"""
    See the module's documentation.
    """
    cdef int n = sd.n

    # Computing the predecessors/distances can only be done if we have less than
    # MAX_UNSIGNED_SHORT vertices. No problem with the eccentricities though as
    # we store them on an integer vector.
    if (predecessors or distances) and n > <unsigned short> -1:
        raise ValueError("The graph backend contains more than " +
                         str(<unsigned short> -1) + " nodes and we cannot " +
                         "compute the matrix of distances/predecessors on " +
                         "something like that !")

    cdef int i
    cdef MemoryAllocator mem = MemoryAllocator()

    # The vertices which have already been visited
    cdef bitset_t seen
    bitset_init(seen, n)

    # The list of waiting vertices, the beginning and the end of the list
    cdef int* waiting_list = <int*> mem.allocarray(n, sizeof(int))
    cdef int waiting_beginning = 0
    cdef int waiting_end = 0

    cdef int source
    cdef int v, u
    cdef uint32_t* p_tmp
    cdef uint32_t* end

    cdef unsigned short *c_predecessors = predecessors
    cdef int* c_distances = <int*> mem.allocarray(n, sizeof(int))

    # The edges are stored in the vector p_edges. This vector contains, from
    # left to right The list of the first vertex's outneighbors, then the
    # second's, then the third's, ...
    #
    # The outneighbors of vertex i are enumerated from
    #
    # p_vertices[i] to p_vertices[i+1] - 1
    # (if p_vertices[i] is equal to p_vertices[i+1], then i has no outneighbours)
    #
    # This data structure is well documented in the module
    # sage.graphs.base.static_sparse_graph
    cdef uint32_t** p_vertices = sd.neighbors

    # We run n different BFS taking each vertex as a source
    for source in range(n):

        # The source is seen
        bitset_set_first_n(seen, 0)
        bitset_add(seen, source)

        # Its parameters can already be set
        c_distances[source] = 0

        if predecessors != NULL:
            c_predecessors[source] = source

        # and added to the queue
        waiting_list[0] = source
        waiting_beginning = 0
        waiting_end = 0

        # For as long as there are vertices left to explore
        while waiting_beginning <= waiting_end:

            # We pick the first one
            v = waiting_list[waiting_beginning]

            p_tmp = p_vertices[v]
            end = p_vertices[v + 1]

            # Iterating over all the outneighbors u of v
            while p_tmp < end:
                u = p_tmp[0]

                # If we notice one of these neighbors is not seen yet, we set
                # its parameters and add it to the queue to be explored later.
                if not bitset_in(seen, u):
                    c_distances[u] = c_distances[v] + 1
                    if predecessors:
                        c_predecessors[u] = v
                    bitset_add(seen, u)
                    waiting_end += 1
                    waiting_list[waiting_end] = u

                p_tmp += 1

            waiting_beginning += 1

        # If not all the vertices have been met
        if bitset_len(seen) < n:
            bitset_complement(seen, seen)
            v = bitset_next(seen, 0)
            while v >= 0:
                c_distances[v] = INT32_MAX
                if predecessors:
                    c_predecessors[v] = -1
                v = bitset_next(seen, v + 1)

            if eccentricity:
                eccentricity[source] = UINT32_MAX

        elif eccentricity:
            eccentricity[source] = c_distances[waiting_list[n - 1]]

        if predecessors:
            c_predecessors += n

        if distances:
            for i in range(n):
                distances[i] = <unsigned short> c_distances[i]
            distances += n

    bitset_free(seen)


cdef inline all_pairs_shortest_path_BFS(gg,
                                        unsigned short* predecessors,
                                        unsigned short* distances,
                                        uint32_t* eccentricity,
                                        vertex_list=None):
    r"""
    See the module's documentation.

    Optional parameter ``vertex_list`` is a list of `n` vertices
    specifying a mapping from `(0, \ldots, n-1)` to vertex labels in
    ``gg``. When ``vertex_list`` is ``None`` (default), the mapping is
    given by the ordering of ``gg.vertices(sort=True)``. When set,
    ``distances[i * n + j]`` is the shortest BFS distance between
    vertices ``vertex_list[i]`` and ``vertex_list[j]``.
    """
    from sage.rings.infinity import Infinity

    cdef list int_to_vertex
    if vertex_list is None:
        int_to_vertex = gg.vertices(sort=True)
    elif set(gg.vertex_iterator()) == set(vertex_list):
        int_to_vertex = vertex_list
    else:
        raise ValueError("parameter vertex_list is incorrect for this graph")

    cdef int n = gg.order()

    # Computing the predecessors/distances can only be done if we have less than
    # MAX_UNSIGNED_SHORT vertices. No problem with the eccentricities though as
    # we store them on an integer vector.
    if (predecessors or distances) and n > <unsigned short> -1:
        raise ValueError("The graph backend contains more than "+
                         str(<unsigned short> -1)+" nodes and we cannot "+
                         "compute the matrix of distances/predecessors on "+
                         "something like that !")

    # Copying the whole graph to obtain the list of neighbors quicker than by
    # calling out_neighbors
    cdef short_digraph sd
    init_short_digraph(sd, gg, edge_labelled=False, vertex_list=int_to_vertex)

    c_all_pairs_shortest_path_BFS(sd, predecessors, distances, eccentricity)

    free_short_digraph(sd)


################
# Predecessors #
################

cdef unsigned short* c_shortest_path_all_pairs(G, vertex_list=None) except NULL:
    r"""
    Return the matrix of predecessors in G.

    The matrix `P` returned has size `n^2`, and is such that vertex `P[u,v]` is
    a predecessor of `v` on a shortest `uv`-path. Hence, this matrix efficiently
    encodes the information of a shortest `uv`-path for any `u,v\in G` : indeed,
    to go from `u` to `v` you should first find a shortest `uP[u,v]`-path, then
    jump from `P[u,v]` to `v` as it is one of its outneighbors.

    Optional parameter ``vertex_list`` is a list of `n` vertices specifying a
    mapping from `(0, \ldots, n-1)` to vertex labels in `G`. When
    ``vertex_list`` is ``None`` (default), the mapping is given by the ordering
    of ``G.vertices(sort=True)``. When set, ``predecessors[i * n + j]`` is the
    predecessor of ``vertex_list[j]`` on the shortest path from
    ``vertex_list[i]`` to ``vertex_list[j]``.
    """

    cdef unsigned int n = G.order()
    cdef unsigned short* distances = <unsigned short*> sig_malloc(n * n * sizeof(unsigned short))
    if not distances:
        raise MemoryError()
    cdef unsigned short* predecessors = <unsigned short*> sig_malloc(n * n * sizeof(unsigned short))
    if not predecessors:
        sig_free(distances)
        raise MemoryError()
    all_pairs_shortest_path_BFS(G, predecessors, distances, NULL, vertex_list=vertex_list)

    sig_free(distances)

    return predecessors


def shortest_path_all_pairs(G):
    r"""
    Return the matrix of predecessors in G.

    The matrix `P` returned has size `n^2`, and is such that vertex `P[u,v]` is
    a predecessor of `v` on a shortest `uv`-path. Hence, this matrix efficiently
    encodes the information of a shortest `uv`-path for any `u,v\in G` : indeed,
    to go from `u` to `v` you should first find a shortest `uP[u,v]`-path, then
    jump from `P[u,v]` to `v` as it is one of its outneighbors.

    EXAMPLES::

        sage: from sage.graphs.distances_all_pairs import shortest_path_all_pairs
        sage: g = graphs.PetersenGraph()
        sage: shortest_path_all_pairs(g)
        {0: {0: None, 1: 0, 2: 1, 3: 4, 4: 0, 5: 0, 6: 1, 7: 5, 8: 5, 9: 4},
         1: {0: 1, 1: None, 2: 1, 3: 2, 4: 0, 5: 0, 6: 1, 7: 2, 8: 6, 9: 6},
         2: {0: 1, 1: 2, 2: None, 3: 2, 4: 3, 5: 7, 6: 1, 7: 2, 8: 3, 9: 7},
         3: {0: 4, 1: 2, 2: 3, 3: None, 4: 3, 5: 8, 6: 8, 7: 2, 8: 3, 9: 4},
         4: {0: 4, 1: 0, 2: 3, 3: 4, 4: None, 5: 0, 6: 9, 7: 9, 8: 3, 9: 4},
         5: {0: 5, 1: 0, 2: 7, 3: 8, 4: 0, 5: None, 6: 8, 7: 5, 8: 5, 9: 7},
         6: {0: 1, 1: 6, 2: 1, 3: 8, 4: 9, 5: 8, 6: None, 7: 9, 8: 6, 9: 6},
         7: {0: 5, 1: 2, 2: 7, 3: 2, 4: 9, 5: 7, 6: 9, 7: None, 8: 5, 9: 7},
         8: {0: 5, 1: 6, 2: 3, 3: 8, 4: 3, 5: 8, 6: 8, 7: 5, 8: None, 9: 6},
         9: {0: 4, 1: 6, 2: 7, 3: 4, 4: 9, 5: 7, 6: 9, 7: 9, 8: 6, 9: None}}
    """
    cdef int n = G.order()

    if not n:
        return {}

    # The order of vertices must be the same as in init_short_digraph
    cdef list int_to_vertex = list(G)
    cdef unsigned short* predecessors = c_shortest_path_all_pairs(G, vertex_list=int_to_vertex)
    cdef unsigned short* c_predecessors = predecessors

    cdef dict d = {}
    cdef dict d_tmp

    cdef int i, j

    for j in range(n):
        d_tmp = {}
        for i in range(n):
            if c_predecessors[i] == <unsigned short> -1:
                d_tmp[int_to_vertex[i]] = None
            else:
                d_tmp[int_to_vertex[i]] = int_to_vertex[c_predecessors[i]]

        d_tmp[int_to_vertex[j]] = None
        d[int_to_vertex[j]] = d_tmp

        c_predecessors += n

    sig_free(predecessors)
    return d


#############
# Distances #
#############

cdef unsigned short * c_distances_all_pairs(G, vertex_list=None) noexcept:
    r"""
    Return the matrix of distances in G.

    The matrix `M` returned is of length `n^2`, and the distance between
    vertices `u` and `v` is `M[u,v]`. The integer corresponding to a vertex is
    its index in the list ``G.vertices(sort=True)`` unless parameter
    ``vertex_list`` is set.

    Optional parameter ``vertex_list`` is a list of `n` vertices specifying a
    mapping from `(0, \ldots, n-1)` to vertex labels in `G`. When set,
    ``distances[i * n + j]`` is the shortest BFS distance between vertices
    ``vertex_list[i]`` and ``vertex_list[j]``.
    """
    cdef unsigned int n = G.order()
    cdef unsigned short* distances = <unsigned short*> sig_malloc(n * n * sizeof(unsigned short))
    if not distances:
        raise MemoryError()
    all_pairs_shortest_path_BFS(G, NULL, distances, NULL, vertex_list=vertex_list)

    return distances


def distances_all_pairs(G):
    r"""
    Return the matrix of distances in G.

    This function returns a double dictionary ``D`` of vertices, in which the
    distance between vertices ``u`` and ``v`` is ``D[u][v]``.

    EXAMPLES::

        sage: from sage.graphs.distances_all_pairs import distances_all_pairs
        sage: g = graphs.PetersenGraph()
        sage: distances_all_pairs(g)
        {0: {0: 0, 1: 1, 2: 2, 3: 2, 4: 1, 5: 1, 6: 2, 7: 2, 8: 2, 9: 2},
        1: {0: 1, 1: 0, 2: 1, 3: 2, 4: 2, 5: 2, 6: 1, 7: 2, 8: 2, 9: 2},
        2: {0: 2, 1: 1, 2: 0, 3: 1, 4: 2, 5: 2, 6: 2, 7: 1, 8: 2, 9: 2},
        3: {0: 2, 1: 2, 2: 1, 3: 0, 4: 1, 5: 2, 6: 2, 7: 2, 8: 1, 9: 2},
        4: {0: 1, 1: 2, 2: 2, 3: 1, 4: 0, 5: 2, 6: 2, 7: 2, 8: 2, 9: 1},
        5: {0: 1, 1: 2, 2: 2, 3: 2, 4: 2, 5: 0, 6: 2, 7: 1, 8: 1, 9: 2},
        6: {0: 2, 1: 1, 2: 2, 3: 2, 4: 2, 5: 2, 6: 0, 7: 2, 8: 1, 9: 1},
        7: {0: 2, 1: 2, 2: 1, 3: 2, 4: 2, 5: 1, 6: 2, 7: 0, 8: 2, 9: 1},
        8: {0: 2, 1: 2, 2: 2, 3: 1, 4: 2, 5: 1, 6: 1, 7: 2, 8: 0, 9: 2},
        9: {0: 2, 1: 2, 2: 2, 3: 2, 4: 1, 5: 2, 6: 1, 7: 1, 8: 2, 9: 0}}
    """
    from sage.rings.infinity import Infinity

    cdef int n = G.order()

    if not n:
        return {}

    # The order of vertices must be the same as in init_short_digraph
    cdef list int_to_vertex = list(G)

    cdef unsigned short* distances = c_distances_all_pairs(G, vertex_list=int_to_vertex)
    cdef unsigned short* c_distances = distances

    cdef dict d = {}
    cdef dict d_tmp

    cdef int i, j

    for j in range(n):
        d_tmp = {}
        for i in range(n):
            if c_distances[i] == <unsigned short> -1:
                d_tmp[int_to_vertex[i]] = Infinity
            else:
                d_tmp[int_to_vertex[i]] = c_distances[i]

        d[int_to_vertex[j]] = d_tmp
        c_distances += n

    sig_free(distances)
    return d


def is_distance_regular(G, parameters=False):
    r"""
    Test if the graph is distance-regular.

    A graph `G` is distance-regular if for any integers `j,k` the value of
    `|\{x:d_G(x,u)=j,x\in V(G)\} \cap \{y:d_G(y,v)=j,y\in V(G)\}|` is constant
    for any two vertices `u,v\in V(G)` at distance `i` from each other. In
    particular `G` is regular, of degree `b_0` (see below), as one can take
    `u=v`.

    Equivalently a graph is distance-regular if there exist integers `b_i,c_i`
    such that for any two vertices `u,v` at distance `i` we have

    * `b_i = |\{x:d_G(x,u)=i+1,x\in V(G)\}\cap N_G(v)\}|, \ 0\leq i\leq d-1`
    * `c_i = |\{x:d_G(x,u)=i-1,x\in V(G)\}\cap N_G(v)\}|, \ 1\leq i\leq d,`

    where `d` is the diameter of the graph.  For more information on
    distance-regular graphs, see the :wikipedia:`Distance-regular_graph`.

    INPUT:

    - ``parameters`` -- boolean (default: ``False``); if set to ``True``, the
      function returns the pair ``(b, c)`` of lists of integers instead of
      a boolean answer (see the definition above)

    .. SEEALSO::

        * :meth:`~sage.graphs.generic_graph.GenericGraph.is_regular`
        * :meth:`~Graph.is_strongly_regular`

    EXAMPLES::

        sage: g = graphs.PetersenGraph()
        sage: g.is_distance_regular()
        True
        sage: g.is_distance_regular(parameters = True)
        ([3, 2, None], [None, 1, 1])

    Cube graphs, which are not strongly regular, are a bit more interesting::

        sage: graphs.CubeGraph(4).is_distance_regular()
        True
        sage: graphs.OddGraph(5).is_distance_regular()
        True

    Disconnected graph::

        sage: (2*graphs.CubeGraph(4)).is_distance_regular()
        True

    TESTS::

        sage: graphs.PathGraph(2).is_distance_regular(parameters=True)
        ([1, None], [None, 1])
        sage: graphs.Tutte12Cage().is_distance_regular(parameters=True)                 # needs networkx
        ([3, 2, 2, 2, 2, 2, None], [None, 1, 1, 1, 1, 1, 3])
    """
    cdef int i, u, v, d, b, c, k
    cdef int n = G.order()
    cdef int infinity = <unsigned short> -1

    if n <= 1:
        return ([], []) if parameters else True

    if not G.is_regular():
        return False
    k = G.degree(next(G.vertex_iterator()))

    # Matrix of distances
    cdef unsigned short* distance_matrix = c_distances_all_pairs(G, vertex_list=list(G))

    # The diameter, i.e. the longest *finite* distance between two vertices
    cdef int diameter = 0
    for i in range(n * n):
        if distance_matrix[i] > diameter and distance_matrix[i] != infinity:
            diameter = distance_matrix[i]

    cdef bitset_t b_tmp
    bitset_init(b_tmp, n)

    # b_distance_matrix[d*n+v] is the set of vertices at distance d from v.
    cdef binary_matrix_t b_distance_matrix
    try:
        binary_matrix_init(b_distance_matrix, n * (diameter + 2), n)
    except MemoryError:
        sig_free(distance_matrix)
        bitset_free(b_tmp)
        raise

    # Fills b_distance_matrix
    for u in range(n):
        for v in range(u, n):
            d = distance_matrix[u * n + v]
            if d != infinity:
                binary_matrix_set1(b_distance_matrix, d * n + u, v)
                binary_matrix_set1(b_distance_matrix, d * n + v, u)

    cdef list bi = [-1 for i in range(diameter + 1)]
    cdef list ci = [-1 for i in range(diameter + 1)]

    # Applying the definition with b_i,c_i
    for u in range(n):
        for v in range(n):
            if u == v:
                continue

            d = distance_matrix[u * n + v]
            if d == infinity:
                continue

            # Computations of b_d and c_d for u,v. We intersect sets stored in
            # b_distance_matrix.
            bitset_and(b_tmp, b_distance_matrix.rows[(d + 1) * n + u], b_distance_matrix.rows[n + v])
            b = bitset_len(b_tmp)
            bitset_and(b_tmp, b_distance_matrix.rows[(d - 1) * n + u], b_distance_matrix.rows[n + v])
            c = bitset_len(b_tmp)

            # Consistency of b_d and c_d
            if bi[d] == -1:
                bi[d] = b
                ci[d] = c

            elif bi[d] != b or ci[d] != c:
                sig_free(distance_matrix)
                binary_matrix_free(b_distance_matrix)
                bitset_free(b_tmp)
                return False

    sig_free(distance_matrix)
    binary_matrix_free(b_distance_matrix)
    bitset_free(b_tmp)

    if parameters:
        bi[0] = k
        bi[diameter] = None
        ci[0] = None
        return bi, ci
    return True


###################################
# Both distances and predecessors #
###################################

def distances_and_predecessors_all_pairs(G):
    r"""
    Return the matrix of distances in G and the matrix of predecessors.

    Distances : the matrix `M` returned is of length `n^2`, and the distance
    between vertices `u` and `v` is `M[u,v]`. The integer corresponding to a
    vertex is its index in the list ``G.vertices(sort=True)``.

    Predecessors : the matrix `P` returned has size `n^2`, and is such that
    vertex `P[u,v]` is a predecessor of `v` on a shortest `uv`-path. Hence, this
    matrix efficiently encodes the information of a shortest `uv`-path for any
    `u,v\in G` : indeed, to go from `u` to `v` you should first find a shortest
    `uP[u,v]`-path, then jump from `P[u,v]` to `v` as it is one of its
    outneighbors.

    EXAMPLES::

        sage: from sage.graphs.distances_all_pairs import distances_and_predecessors_all_pairs
        sage: g = graphs.PetersenGraph()
        sage: distances_and_predecessors_all_pairs(g)
        ({0: {0: 0, 1: 1, 2: 2, 3: 2, 4: 1, 5: 1, 6: 2, 7: 2, 8: 2, 9: 2},
          1: {0: 1, 1: 0, 2: 1, 3: 2, 4: 2, 5: 2, 6: 1, 7: 2, 8: 2, 9: 2},
          2: {0: 2, 1: 1, 2: 0, 3: 1, 4: 2, 5: 2, 6: 2, 7: 1, 8: 2, 9: 2},
          3: {0: 2, 1: 2, 2: 1, 3: 0, 4: 1, 5: 2, 6: 2, 7: 2, 8: 1, 9: 2},
          4: {0: 1, 1: 2, 2: 2, 3: 1, 4: 0, 5: 2, 6: 2, 7: 2, 8: 2, 9: 1},
          5: {0: 1, 1: 2, 2: 2, 3: 2, 4: 2, 5: 0, 6: 2, 7: 1, 8: 1, 9: 2},
          6: {0: 2, 1: 1, 2: 2, 3: 2, 4: 2, 5: 2, 6: 0, 7: 2, 8: 1, 9: 1},
          7: {0: 2, 1: 2, 2: 1, 3: 2, 4: 2, 5: 1, 6: 2, 7: 0, 8: 2, 9: 1},
          8: {0: 2, 1: 2, 2: 2, 3: 1, 4: 2, 5: 1, 6: 1, 7: 2, 8: 0, 9: 2},
          9: {0: 2, 1: 2, 2: 2, 3: 2, 4: 1, 5: 2, 6: 1, 7: 1, 8: 2, 9: 0}},
         {0: {0: None, 1: 0, 2: 1, 3: 4, 4: 0, 5: 0, 6: 1, 7: 5, 8: 5, 9: 4},
          1: {0: 1, 1: None, 2: 1, 3: 2, 4: 0, 5: 0, 6: 1, 7: 2, 8: 6, 9: 6},
          2: {0: 1, 1: 2, 2: None, 3: 2, 4: 3, 5: 7, 6: 1, 7: 2, 8: 3, 9: 7},
          3: {0: 4, 1: 2, 2: 3, 3: None, 4: 3, 5: 8, 6: 8, 7: 2, 8: 3, 9: 4},
          4: {0: 4, 1: 0, 2: 3, 3: 4, 4: None, 5: 0, 6: 9, 7: 9, 8: 3, 9: 4},
          5: {0: 5, 1: 0, 2: 7, 3: 8, 4: 0, 5: None, 6: 8, 7: 5, 8: 5, 9: 7},
          6: {0: 1, 1: 6, 2: 1, 3: 8, 4: 9, 5: 8, 6: None, 7: 9, 8: 6, 9: 6},
          7: {0: 5, 1: 2, 2: 7, 3: 2, 4: 9, 5: 7, 6: 9, 7: None, 8: 5, 9: 7},
          8: {0: 5, 1: 6, 2: 3, 3: 8, 4: 3, 5: 8, 6: 8, 7: 5, 8: None, 9: 6},
          9: {0: 4, 1: 6, 2: 7, 3: 4, 4: 9, 5: 7, 6: 9, 7: 9, 8: 6, 9: None}})
    """
    from sage.rings.infinity import Infinity
    cdef unsigned int n = G.order()

    if not n:
        return {}, {}

    cdef MemoryAllocator mem = MemoryAllocator()
    cdef unsigned short* c_distances = <unsigned short*> mem.malloc(n * n * sizeof(unsigned short))
    cdef unsigned short* c_predecessor = <unsigned short*> mem.malloc(n * n * sizeof(unsigned short))

    # The order of vertices must be the same as in init_short_digraph
    cdef list int_to_vertex = list(G)

    all_pairs_shortest_path_BFS(G, c_predecessor, c_distances, NULL, vertex_list=int_to_vertex)

    cdef dict d_distance = {}
    cdef dict d_predecessor = {}
    cdef dict t_distance = {}
    cdef dict t_predecessor = {}

    cdef unsigned int i, j

    for j in range(n):
        t_distance = {}
        t_predecessor = {}

        for i in range(n):

            if c_distances[i] != <unsigned short> -1:
                t_distance[int_to_vertex[i]] = c_distances[i]
                t_predecessor[int_to_vertex[i]] = int_to_vertex[c_predecessor[i]]

        t_predecessor[int_to_vertex[j]] = None

        d_distance[int_to_vertex[j]] = t_distance
        d_predecessor[int_to_vertex[j]] = t_predecessor

        c_distances += n
        c_predecessor += n

    return d_distance, d_predecessor


################
# Eccentricity #
################

cdef uint32_t * c_eccentricity(G, vertex_list=None) except NULL:
    r"""
    Return the vector of eccentricities in G.

    The array returned is of length `n`, and by default its `i`-th component is
    the eccentricity of the `i`-th vertex in ``G.vertices(sort=True)``.

    Optional parameter ``vertex_list`` is a list of `n` vertices specifying a
    mapping from `(0, \ldots, n-1)` to vertex labels in `G`. When set,
    ``ecc[i]`` is the eccentricity of vertex ``vertex_list[i]``.
    """
    cdef unsigned int n = G.order()

    cdef uint32_t * ecc = <uint32_t *> sig_calloc(n, sizeof(uint32_t))
    if not ecc:
        raise MemoryError()
    all_pairs_shortest_path_BFS(G, NULL, NULL, ecc, vertex_list=vertex_list)

    return ecc


cdef uint32_t * c_eccentricity_bounding(short_digraph sd) except NULL:
    r"""
    Return the vector of eccentricities using the algorithm of [TK2013]_.

    The array returned is of length `n`, and its `i`-th component is the
    eccentricity of vertex `i` in ``sd``.

    This method assumes that ``sd`` is an undirected graph.

    The algorithm proposed in [TK2013]_ is based on the observation that for all
    nodes `v,w\in V`, we have `\max(ecc[v]-d(v,w), d(v,w))\leq ecc[w] \leq
    ecc[v] + d(v,w)`. Also the algorithms iteratively improves upper and lower
    bounds on the eccentricity of each node until no further improvements can be
    done. This algorithm offers good running time reduction on scale-free graphs.
    """
    cdef unsigned int n = sd.n

    # allocated some data structures
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef uint32_t * distances = <uint32_t *> mem.malloc(3 * n * sizeof(uint32_t))
    cdef uint32_t * LB = <uint32_t *> sig_calloc(n, sizeof(uint32_t))
    if not distances or not LB:
        sig_free(LB)
        free_short_digraph(sd)
        raise MemoryError()
    cdef uint32_t * waiting_list = distances + n
    cdef uint32_t * UB = distances + 2 * n
    memset(UB, -1, n * sizeof(uint32_t))
    cdef bitset_t seen
    bitset_init(seen, n)

    cdef uint32_t v, w, next_v, cpt = 0

    # The first vertex is the one with largest degree
    next_v = max((out_degree(sd, v), v) for v in range(n))[1]

    sig_on()
    while next_v != UINT32_MAX:

        v = next_v
        cpt += 1

        # Compute the exact eccentricity of v
        LB[v] = simple_BFS(sd, v, distances, NULL, waiting_list, seen)

        if LB[v] == UINT32_MAX:
            # The graph is not connected. We set maximum value and exit.
            for w in range(n):
                LB[w] = UINT32_MAX
            break

        # Improve the bounds on the eccentricity of other vertices and select
        # source of the next BFS
        next_v = UINT32_MAX
        for w in range(n):
            LB[w] = max(LB[w], max(LB[v] - distances[w], distances[w]))
            UB[w] = min(UB[w], LB[v] + distances[w])
            if LB[w] == UB[w]:
                continue
            elif next_v == UINT32_MAX or (not cpt % 2 and LB[w] < LB[next_v]) or (cpt % 2 and UB[w] > UB[next_v]):
                # The next vertex is either the vertex with largest upper bound
                # or smallest lower bound
                next_v = w

    sig_off()

    bitset_free(seen)

    return LB


cdef uint32_t * c_eccentricity_DHV(short_digraph sd) except NULL:
    r"""
    Return the vector of eccentricities using the algorithm of [Dragan2018]_.

    The array returned is of length `n`, and its `i`-th component is the
    eccentricity of vertex `i` in ``sd``.

    This method assumes that ``sd`` is an undirected graph.

    The algorithm proposed in [Dragan2018]_ is an improvement of the algorithm
    proposed in [TK2013]_. It is also based on the observation that for all
    nodes `v,w\in V`, we have `\max(ecc[v]-d(v,w), d(v,w))\leq ecc[w] \leq
    ecc[v] + d(v,w)`. Also the algorithms iteratively improves upper and lower
    bounds on the eccentricity of each vertex until no further improvements can
    be done. The difference with [TK2013]_ is in the order in which improvements
    are done.

    EXAMPLES::

        sage: from sage.graphs.distances_all_pairs import eccentricity
        sage: G = graphs.PathGraph(5)
        sage: eccentricity(G, algorithm='DHV')
        [4, 3, 2, 3, 4]

    TESTS:

        sage: G = graphs.RandomBarabasiAlbert(50, 2)                                    # needs networkx
        sage: eccentricity(G, algorithm='bounds') == eccentricity(G, algorithm='DHV')   # needs networkx
        True
    """
    cdef uint32_t n = sd.n
    if not n:
        return NULL

    cdef MemoryAllocator mem = MemoryAllocator()
    cdef uint32_t * distances = <uint32_t *> mem.malloc(3 * n * sizeof(uint32_t))
    # For storing upper bounds on eccentricity of nodes
    cdef uint32_t * ecc_upper_bound = <uint32_t *>sig_calloc(n, sizeof(uint32_t))
    if not distances or not ecc_upper_bound:
        sig_free(ecc_upper_bound)
        free_short_digraph(sd)
        raise MemoryError()

    cdef uint32_t * waiting_list = distances + n
    # For storing lower bounds on eccentricity of nodes
    cdef uint32_t * ecc_lower_bound = distances + 2 * n
    memset(ecc_upper_bound, <char>-1, n * sizeof(uint32_t))
    memset(ecc_lower_bound, 0, n * sizeof(uint32_t))

    cdef uint32_t u, ecc_u
    cdef uint32_t antipode, ecc_antipode
    cdef uint32_t v, tmp
    cdef size_t i, idx
    cdef bitset_t seen
    bitset_init(seen, n)

    cdef list active = list(range(n))

    # Algorithm
    while active:
        # Select vertex with minimum eccentricity in active and update
        # eccentricity upper bounds.
        # For this, we select u with minimum eccentricity lower bound in active
        # if ecc_u == ecc_lb[u], we are done. Otherwise, we update eccentricity
        # lower bounds and repeat

        tmp = UINT32_MAX
        for i, v in enumerate(active):
            if ecc_lower_bound[v] < tmp:
                tmp = ecc_lower_bound[v]
                idx = i
        active[idx], active[-1] = active[-1], active[idx]
        u = active.pop()
        ecc_u = simple_BFS(sd, u, distances, NULL, waiting_list, seen)
        ecc_upper_bound[u] = ecc_u

        if ecc_u == UINT32_MAX:  # Disconnected graph
            break

        if ecc_u == ecc_lower_bound[u]:
            # We found the good vertex.
            # Update eccentricity upper bounds and remove from active those
            # vertices for which gap is closed
            i = 0
            while i < len(active):
                v = active[i]
                ecc_upper_bound[v] = min(ecc_upper_bound[v], distances[v] + ecc_u)
                if ecc_upper_bound[v] == ecc_lower_bound[v]:
                    active[i] = active[-1]
                    active.pop()
                else:
                    i += 1

        else:
            # u was not a good choice.
            # We use its antipode to update eccentricity lower bounds.
            # Observe that this antipode might have already been seen.
            antipode = waiting_list[n-1]
            for i, v in enumerate(active):
                if v == antipode:
                    active[i] = active[-1]
                    active.pop()
                    break

            ecc_antipode = simple_BFS(sd, antipode, distances, NULL, waiting_list, seen)
            ecc_upper_bound[antipode] = ecc_antipode

            # Update eccentricity lower bounds and remove from active those
            # vertices for which the gap is closed
            i = 0
            while i < len(active):
                v = active[i]
                ecc_lower_bound[v] = max(ecc_lower_bound[v], distances[v])
                if ecc_upper_bound[v] == ecc_lower_bound[v]:
                    active[i] = active[-1]
                    active.pop()
                else:
                    i += 1

    bitset_free(seen)

    return ecc_upper_bound


def eccentricity(G, algorithm='standard', vertex_list=None):
    r"""
    Return the vector of eccentricities in G.

    The array returned is of length `n`, and its `i`-th component is the
    eccentricity of the `i`-th vertex in ``G.vertices(sort=True)``.

    INPUT:

    - ``G`` -- a Graph or a DiGraph

    - ``algorithm`` -- string (default: ``'standard'``); name of the method used
      to compute the eccentricity of the vertices

      - ``'standard'`` -- computes eccentricity by performing a BFS from each
        vertex

      - ``'bounds'`` -- computes eccentricity using the fast algorithm proposed
        in [TK2013]_ for undirected graphs

      - ``'DHV'`` -- computes all eccentricities of undirected graph using the
        algorithm proposed in [Dragan2018]_

    - ``vertex_list`` -- list (default: ``None``); a list of `n` vertices
      specifying a mapping from `(0, \ldots, n-1)` to vertex labels in `G`. When
      set, ``ecc[i]`` is the eccentricity of vertex ``vertex_list[i]``. When
      ``vertex_list`` is ``None``, ``ecc[i]`` is the eccentricity of vertex
      ``G.vertices(sort=True)[i]``.

    EXAMPLES::

        sage: from sage.graphs.distances_all_pairs import eccentricity
        sage: g = graphs.PetersenGraph()
        sage: eccentricity(g)
        [2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
        sage: g.add_edge(0, g.add_vertex())
        sage: V = list(g)
        sage: eccentricity(g, vertex_list=V)
        [2, 2, 3, 3, 2, 2, 3, 3, 3, 3, 3]
        sage: eccentricity(g, vertex_list=V[::-1])
        [3, 3, 3, 3, 3, 2, 2, 3, 3, 2, 2]

    TESTS:

    All algorithms are valid::

        sage: from sage.graphs.distances_all_pairs import eccentricity
        sage: g = graphs.RandomGNP(50, .1)
        sage: ecc = eccentricity(g, algorithm='standard')
        sage: ecc == eccentricity(g, algorithm='bounds')
        True
        sage: ecc == eccentricity(g, algorithm='DHV')
        True

    Case of not (strongly) connected (directed) graph::

        sage: from sage.graphs.distances_all_pairs import eccentricity
        sage: g = 2*graphs.PathGraph(2)
        sage: eccentricity(g, algorithm='bounds')
        [+Infinity, +Infinity, +Infinity, +Infinity]
        sage: g = digraphs.Path(3)
        sage: eccentricity(g, algorithm='standard')
        [2, +Infinity, +Infinity]

    The bounds algorithm is for Graph only::

        sage: from sage.graphs.distances_all_pairs import eccentricity
        sage: g = digraphs.Circuit(2)
        sage: eccentricity(g, algorithm='bounds')
        Traceback (most recent call last):
        ...
        ValueError: the 'bounds' algorithm only works on undirected graphs

    Asking for unknown algorithm::

        sage: from sage.graphs.distances_all_pairs import eccentricity
        sage: g = graphs.PathGraph(2)
        sage: eccentricity(g, algorithm='Nice Jazz Festival')
        Traceback (most recent call last):
        ...
        ValueError: unknown algorithm 'Nice Jazz Festival', please contribute

    Invalid value for parameter vertex_list::

        sage: from sage.graphs.distances_all_pairs import eccentricity
        sage: g = graphs.PathGraph(2)
        sage: eccentricity(g, vertex_list=[0, 1, 2])
        Traceback (most recent call last):
        ...
        ValueError: parameter vertex_list is incorrect for this graph
    """
    from sage.rings.infinity import Infinity
    cdef int n = G.order()
    cdef short_digraph sd

    # Trivial cases
    if algorithm not in ['standard', 'bounds', 'DHV']:
        raise ValueError("unknown algorithm '{}', please contribute".format(algorithm))
    if not n:
        return []
    elif G.is_directed() and algorithm in ['bounds', 'DHV']:
        raise ValueError("the 'bounds' algorithm only works on undirected graphs")
    elif not G.is_connected():
        return [Infinity] * n

    cdef list int_to_vertex
    if vertex_list is None:
        int_to_vertex = G.vertices(sort=True)
    elif len(vertex_list) == n and set(vertex_list) == set(G):
        int_to_vertex = vertex_list
    else:
        raise ValueError("parameter vertex_list is incorrect for this graph")

    cdef uint32_t* ecc

    if algorithm == "standard":
        ecc = c_eccentricity(G, vertex_list=int_to_vertex)

    else:
        init_short_digraph(sd, G, edge_labelled=False, vertex_list=vertex_list)

        if algorithm == "DHV":
            ecc = c_eccentricity_DHV(sd)
        else:  # "bounds"
            ecc = c_eccentricity_bounding(sd)

        free_short_digraph(sd)

    from sage.rings.integer import Integer
    cdef list l_ecc = [Integer(ecc[i]) if ecc[i] != UINT32_MAX else +Infinity for i in range(n)]

    sig_free(ecc)

    return l_ecc


############
# Diameter #
############

cdef uint32_t diameter_lower_bound_2sweep(short_digraph g,
                                          uint32_t source,
                                          uint32_t* distances,
                                          uint32_t* predecessors,
                                          uint32_t* waiting_list,
                                          bitset_t seen) noexcept:
    """
    Compute a lower bound on the diameter using the 2-sweep algorithm.

    This method computes a lower bound on the diameter of an unweighted
    undirected graph using 2 BFS, as proposed in [MLH2008]_.  It first selects a
    vertex `v` that is at largest distance from an initial vertex `source` using
    BFS. Then it performs a second BFS from `v`. The largest distance from `v`
    is returned as a lower bound on the diameter of `G`.  The time complexity of
    this method is linear in the size of the graph.


    INPUT:

    - ``g`` -- a short_digraph

    - ``source`` -- starting node of the BFS

    - ``distances`` -- array of size ``n`` to store BFS distances from `v`, the
      vertex at largest distance from ``source`` from which we start the second
      BFS. This method assumes that this array has already been allocated.
      However, there is no need to initialize it.

    - ``predecessors`` -- array of size ``n`` to store the first predecessor of
      each vertex during the BFS search from `v`. The predecessor of `v` is
      itself. This method assumes that this array has already been allocated.
      However, it is possible to pass a ``NULL`` pointer in which case the
      predecessors are not recorded.

    - ``waiting_list`` -- array of size ``n`` to store the order in which the
      vertices are visited during the BFS search from `v`. This method assumes
      that this array has already been allocated. However, there is no need to
      initialize it.

    - ``seen`` -- bitset of size ``n`` that must be initialized before calling
      this method (i.e., bitset_init(seen, n)). However, there is no need to
      clear it.
    """
    cdef uint32_t LB

    # We do a first BFS from source and get the eccentricity of source
    LB = simple_BFS(g, source, distances, NULL, waiting_list, seen)

    # If the eccentricity of the source is infinite (very large number), the
    # graph is not connected and so its diameter is infinite.
    if LB == UINT32_MAX:
        return UINT32_MAX

    # Then we perform a second BFS from the last visited vertex
    source = waiting_list[g.n - 1]
    LB = simple_BFS(g, source, distances, predecessors, waiting_list, seen)

    # We return the computed lower bound
    return LB


cdef tuple diameter_lower_bound_2Dsweep(short_digraph g,
                                        short_digraph rev_g,
                                        uint32_t source):
    r"""
    Lower bound on the diameter of digraph using directed version of 2-sweep.

    This method computes a lower bound on the diameter of an unweighted directed
    graph using directed version of the 2-sweep algorithm [Broder2000]_.
    In first part, it performs a forward BFS from `source` and selects a vertex
    `vf` at a maximum distance from `source` and then it calculates backward
    eccentricity of `vf` using a backward BFS from `vf`. In second part, it
    performs backward BFS from `source` and selects a vertex `vb` from which
    `source` is at maximum distance and then it calculates forward eccentricity
    of `vb` using a forward BFS from `vb`. It then calculates lower bound LB of
    diameter as the maximum of backward eccentricity of `vf` and forward
    eccentricity of `vb` and `s` as respective vertex.
    This method returns (`LB`, `s`, `m`, `d`), where `LB` is best found lower
    bound on diameter, `s` is vertex whose forward/backward eccentricity is
    `LB`, `d` is vertex at a distance `LB` from/to `s`, `m` is vertex at
    distance `LB/2` from/to both `s` and `d`.

    INPUT:

    - ``g`` -- a short_digraph

    - ``rev_g`` -- a copy of `g` with edges reversed

    - ``source`` -- starting node of the forward and backward BFS

    TESTS:

    Diameter of weakly connected digraph is infinity ::

        sage: from sage.graphs.distances_all_pairs import diameter
        sage: G = DiGraph([(0,1)])
        sage: diameter(G, algorithm='2Dsweep')
        +Infinity
    """
    cdef uint32_t LB_1, LB_2, LB, LB_m, m, s, d
    cdef uint32_t n = g.n
    cdef uint32_t source_1 = source
    cdef uint32_t source_2 = source
    cdef bitset_t seen_1, seen_2

    # Memory allocation
    cdef MemoryAllocator mem = MemoryAllocator()
    bitset_init(seen_1, n)
    bitset_init(seen_2, n)
    cdef uint32_t * distances_1 = <uint32_t *>mem.malloc(3 * n * sizeof(uint32_t))
    cdef uint32_t * distances_2 = <uint32_t *>mem.malloc(3 * n * sizeof(uint32_t))
    if not distances_1 or not distances_2:
        bitset_free(seen_1)
        bitset_free(seen_2)
        raise MemoryError()

    cdef uint32_t * predecessors_1 = distances_1 + n
    cdef uint32_t * predecessors_2 = distances_2 + n
    cdef uint32_t * waiting_list_1 = distances_1 + 2 * n
    cdef uint32_t * waiting_list_2 = distances_2 + 2 * n

    # we perform forward BFS from source and get its forward eccentricity
    LB_1 = simple_BFS(g, source_1, distances_1, NULL, waiting_list_1, seen_1)

    # if forward eccentricity of source is infinite, then graph is
    # not strongly connected and its diameter is infinite
    if LB_1 == UINT32_MAX:
        bitset_free(seen_1)
        bitset_free(seen_2)
        return (UINT32_MAX, 0, 0, 0)

    # we perform backward BFS from source and get its backward eccentricity
    LB_2 = simple_BFS(rev_g, source_2, distances_2, NULL, waiting_list_2, seen_2)

    # if backward eccentricity of source is infinite, then graph is
    # not strongly connected and its diameter is infinite
    if LB_2 == UINT32_MAX:
        bitset_free(seen_1)
        bitset_free(seen_2)
        return (UINT32_MAX, 0, 0, 0)

    # Then we perform backward BFS from the last visited vertex of forward BFS
    # from source and obtain its backward eccentricity.
    source_1 = waiting_list_1[n - 1]
    LB_1 = simple_BFS(rev_g, source_1, distances_1, predecessors_1, waiting_list_1, seen_1)

    # Then we perform forward BFS from the last visited vertex of backward BFS
    # from source and obtain its forward eccentricity.
    source_2 = waiting_list_2[n - 1]
    LB_2 = simple_BFS(g, source_2, distances_2, predecessors_2, waiting_list_2, seen_2)

    # we select best found lower bound as LB, s and d as source and destination
    # of that BFS call and m as vertex at a distance LB/2 from/to both s and d
    if LB_1 < LB_2:
        LB = LB_2
        s = waiting_list_2[0]
        d = waiting_list_2[n - 1]
        LB_m = LB_2 / 2
        m = d
        while distances_2[m] > LB_m:
            m = predecessors_2[m]
    else:
        LB = LB_1
        s = waiting_list_1[0]
        d = waiting_list_1[n - 1]
        LB_m = LB_1 / 2
        m = d
        while distances_1[m] > LB_m:
            m = predecessors_1[m]

    bitset_free(seen_1)
    bitset_free(seen_2)

    return (LB, s, m, d)


cdef tuple diameter_lower_bound_multi_sweep(short_digraph g,
                                            uint32_t source):
    """
    Lower bound on the diameter using multi-sweep.

    This method computes a lower bound on the diameter of an unweighted
    undirected graph using several iterations of the 2-sweep algorithms
    [CGHLM2013]_. Roughly, it first uses 2-sweep to identify two vertices `s`
    and `d` that are far apart. Then it selects a vertex `m` that is at same
    distance from `s` and `d`.  This vertex `m` will serve as the new source for
    another iteration of the 2-sweep algorithm that may improve the current
    lower bound on the diameter.  This process is repeated as long as the lower
    bound on the diameter is improved.

    The method returns a 4-tuple (LB, s, m, d), where LB is the best found lower
    bound on the diameter, s is a vertex of eccentricity LB, d is a vertex at
    distance LB from s, and m is a vertex at distance LB/2 from both s and d.

    INPUT:

    - ``g`` -- a short_digraph

    - ``source`` -- starting node of the BFS
    """
    # The while loop below might not be entered so we have to make sure that
    # s and d which are returned are initialized.
    cdef uint32_t LB, tmp, m
    cdef uint32_t s = 0
    cdef uint32_t d = 0
    cdef uint32_t n = g.n

    # Allocate some arrays and a bitset
    cdef bitset_t seen
    bitset_init(seen, n)
    cdef uint32_t * distances = <uint32_t *>sig_malloc(3 * n * sizeof(uint32_t))
    if not distances:
        bitset_free(seen)
        raise MemoryError()
    cdef uint32_t * predecessors = distances + n
    cdef uint32_t * waiting_list = distances + 2 * n

    # We perform a first 2sweep call from source. If the returned value is a
    # very large number, the graph is not connected and so the diameter is
    # infinite.
    tmp = diameter_lower_bound_2sweep(g, source, distances, predecessors, waiting_list, seen)
    if tmp == UINT32_MAX:
        sig_free(distances)
        bitset_free(seen)
        return (UINT32_MAX, 0, 0, 0)

    # We perform new 2sweep calls for as long as we are able to improve the
    # lower bound.
    LB = 0
    m = source
    while tmp > LB:

        LB = tmp

        # We store the vertices s, m, d of the last BFS call. For vertex m, we
        # search for a vertex of eccentricity LB/2. This vertex will serve as
        # the source for the next 2sweep call.
        s = waiting_list[0]
        d = waiting_list[n - 1]
        LB_2 = LB / 2
        m = d
        while distances[m] > LB_2:
            m = predecessors[m]

        # We perform a new 2sweep call from m
        tmp = diameter_lower_bound_2sweep(g, m, distances, predecessors, waiting_list, seen)

    sig_free(distances)
    bitset_free(seen)

    return (LB, s, m, d)


cdef uint32_t diameter_iFUB(short_digraph g,
                            uint32_t source) noexcept:
    """
    Compute the diameter of the input Graph using the ``iFUB`` algorithm.

    The ``iFUB`` (iterative Fringe Upper Bound) algorithm calculates the exact
    value of the diameter of a unweighted undirected graph [CGILM2010]_. This
    algorithms starts with a vertex found through a multi-sweep call (a
    refinement of the 4sweep method). The worst case time complexity of the iFUB
    algorithm is `O(nm)`, but it can be very fast in practice. See the code's
    documentation and [CGHLM2013]_ for more details.

    INPUT:

    - ``g`` -- a short_digraph

    - ``source`` -- starting node of the first BFS
    """
    cdef uint32_t i, LB, m
    cdef uint32_t n = g.n

    # We select a vertex m with low eccentricity using multi-sweep
    LB, _, m, _ = diameter_lower_bound_multi_sweep(g, source)

    # If the lower bound is a very large number, it means that the graph is not
    # connected and so the diameter is infinite.
    if LB == UINT32_MAX:
        return LB

    # We allocate some arrays and a bitset
    cdef bitset_t seen
    bitset_init(seen, n)
    cdef uint32_t * distances = <uint32_t *>sig_malloc(4 * n * sizeof(uint32_t))
    if not distances:
        bitset_free(seen)
        raise MemoryError()
    cdef uint32_t* waiting_list = distances + n
    cdef uint32_t* layer = distances + 2 * n
    cdef uint32_t* order = distances + 3 * n

    # We order the vertices by decreasing layers. This is the inverse order of a
    # BFS from m, and so the inverse order of array waiting_list. Distances are
    # stored in array layer.
    LB = simple_BFS(g, m, layer, NULL, waiting_list, seen)
    for i in range(n):
        order[i] = waiting_list[n - i - 1]

    # The algorithm:
    #
    # The diameter of the graph is equal to the maximum eccentricity of a
    # vertex. Let m be any vertex, and let V be partitioned into A u B where:
    #
    #    d(m,a)<=i for all a \in A
    #    d(m,b)>=i for all b \in B
    #
    # As all vertices from A are at distance <=2i from each other, a vertex a
    # from A with eccentricity ecc(a)>2i is at distance ecc(a) from some vertex
    # b\in B.
    #
    # Consequently, if we have already computed the eccentricity of all
    # vertices in B and know that the diameter is >2i, then we do not have to
    # compute the eccentricity of vertices in A.
    #
    # Now, we compute the maximum eccentricity of all vertices, ordered
    # decreasingly according to their distance to m. We stop when we know that
    # the eccentricity of the unexplored vertices is smaller than the max
    # eccentricity already found.
    i = 0
    while 2 * layer[order[i]] > LB and i < n:
        tmp = simple_BFS(g, order[i], distances, NULL, waiting_list, seen)
        i += 1

        # We update the lower bound
        if tmp > LB:
            LB = tmp

    sig_free(distances)
    bitset_free(seen)

    # We finally return the computed diameter
    return LB


cdef uint32_t diameter_DiFUB(short_digraph sd,
                             uint32_t source) noexcept:
    r"""
    Return the diameter of unweighted directed graph.

    The ``DiFUB`` (Directed iterative Fringe Upper Bound) algorithm calculates
    the exact value of the diameter of an unweighted directed graph [CGLM2012]_.

    This algorithm starts from a vertex found through a 2Dsweep call (a directed
    version of the 2sweep method). The worst case time complexity of the DiFUB
    algorithm is `O(nm)`, but it can be very fast in practice. See the code's
    documentation and [CGLM2012]_ for more details.

    If the digraph is not strongly connected, the returned value is infinity.

    INPUT:

    - ``sd`` -- a short_digraph

    - ``source`` -- starting node of the first BFS

    TESTS:

    The diameter of a weakly connected digraph is infinity ::

        sage: from sage.graphs.distances_all_pairs import diameter
        sage: G = digraphs.Path(5)
        sage: diameter(G, algorithm='DiFUB')
        +Infinity
    """
    cdef uint32_t n = sd.n

    if n <= 1:  # Trivial case
        return 0

    cdef short_digraph rev_sd  # Copy of sd with edges reversed
    init_reverse(rev_sd, sd)

    cdef uint32_t LB, m, LB_1, LB_2, UB
    cdef size_t i
    cdef bitset_t seen

    # We select a vertex with low eccentricity using 2Dsweep
    LB, _, m, _ = diameter_lower_bound_2Dsweep(sd, rev_sd, source)

    # If the lower bound is a very large number, it means that the digraph is
    # not strongly connected and so the diameter is infinite.
    if LB == UINT32_MAX:
        return LB

    # We allocate some arrays and a bitset
    bitset_init(seen, n)
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef uint32_t * distances = <uint32_t *>mem.malloc(6 * n * sizeof(uint32_t))

    if not distances:
        bitset_free(seen)
        raise MemoryError()

    cdef uint32_t * waiting_list = distances + n
    cdef uint32_t * order_1 = distances + 2 * n
    cdef uint32_t * order_2 = distances + 3 * n
    cdef uint32_t * layer_1 = distances + 4 * n
    cdef uint32_t * layer_2 = distances + 5 * n

    # We order the vertices by decreasing forward / backward layers. This is the
    # inverse order of a forward / backward BFS from m, and so the inverse order
    # of array waiting_list. Forward / Backward distances are stored in arrays
    # layer_1 / layer_2 respectively.
    LB_1 = simple_BFS(sd, m, layer_1, NULL, waiting_list, seen)
    for i in range(n):
        order_1[i] = waiting_list[n - i - 1]

    LB_2 = simple_BFS(rev_sd, m, layer_2, NULL, waiting_list, seen)
    for i in range(n):
        order_2[i] = waiting_list[n - i - 1]

    # update the lower bound
    LB = max(LB, LB_1, LB_2)

    if LB == UINT32_MAX:  # Not strongly connected case
        return LB

    # The algorithm:
    #
    # The diameter of the digraph is equal to the maximum forward or backward
    # eccentricity of a vertex. The algorithm is based on the following two
    # observations:
    # 1). All the nodes `x` above the level `i` in Backward BFS of `m` having
    # forward eccentricity greater than `2(i-1)` have a corresponding node `y`,
    # whose backward eccentricity is greater than or equal to the forward
    # eccentricity of `x`, below or on the level `i` in Forward BFS of `m`.
    #
    # 2). All the nodes `x` above the level `i` in Forward BFS of `m` having
    # backward eccentricity greater than `2(i-1)` have a corresponding node `y`,
    # whose forward eccentricity is greater than or equal to the backward
    # eccentricity of `x`, below or on the level `i` in Backward BFS of `m`
    #
    # Therefore, we calculate backward / forward eccentricity of all nodes at
    # level `i` in Forward / Backward BFS of `m` respectively. And their
    # maximum is `LB`. If `LB` is greater than `2(max distance at next level)`
    # then we are done, else we proceed further.

    i = 0
    UB = max(2 * layer_1[order_1[i]], 2 * layer_2[order_2[i]])

    while LB < UB:
        LB_1 = simple_BFS(rev_sd, order_1[i], distances, NULL, waiting_list, seen)
        LB_2 = simple_BFS(sd, order_2[i], distances, NULL, waiting_list, seen)

        # update the lower bound
        LB = max(LB, LB_1, LB_2)
        i += 1

        if LB == UINT32_MAX or i == n:
            break
        # maximum distance at next level
        UB = max(2 * layer_1[order_1[i]], 2 * layer_2[order_2[i]])

    bitset_free(seen)
    free_short_digraph(rev_sd)

    # Finally return the computed diameter
    return LB


cdef uint32_t diameter_DHV(short_digraph g) noexcept:
    r"""
    Return the diameter of unweighted graph `g`.

    This method computes the diameter of unweighted undirected graph using the
    algorithm proposed in [Dragan2018]_.

    This method returns Infinity if graph is not connected.

    INPUT:

    - ``g`` -- a short_digraph

    EXAMPLES::

        sage: from sage.graphs.distances_all_pairs import diameter
        sage: G = graphs.PathGraph(5)
        sage: diameter(G, algorithm='DHV')
        4

    TESTS:

        sage: G = graphs.RandomGNP(20,0.3)
        sage: G.diameter() == diameter(G, algorithm='DHV')
        True

        sage: G = Graph([(0, 1)])
        sage: diameter(G, algorithm='DHV')
        1
    """
    cdef uint32_t n = g.n
    if n <= 1:
        return 0

    cdef MemoryAllocator mem = MemoryAllocator()
    cdef uint32_t * distances = <uint32_t *>mem.malloc(4 * n * sizeof(uint32_t))
    if not distances:
        raise MemoryError()

    cdef uint32_t * waiting_list = distances + n

    # For storing upper and lower bounds on eccentricity of nodes
    cdef uint32_t * ecc_upper_bound = distances + 2 * n
    cdef uint32_t * ecc_lower_bound = distances + 3 * n
    memset(ecc_upper_bound, <char>-1, n * sizeof(uint32_t))
    memset(ecc_lower_bound, 0, n * sizeof(uint32_t))

    cdef uint32_t u, ecc_u
    cdef uint32_t x, ecc_x
    cdef uint32_t antipode, ecc_antipode
    cdef uint32_t LB = 0
    cdef uint32_t UB = UINT32_MAX
    cdef uint32_t v, tmp
    cdef size_t i, idx
    cdef bitset_t seen
    bitset_init(seen, n)

    cdef list active = list(range(n))

    # Algorithm
    while LB < UB and active:
        # 1. Select vertex u with maximum eccentricity upper bound
        tmp = 0
        for i, v in enumerate(active):
            if ecc_upper_bound[v] > tmp:
                tmp = ecc_upper_bound[v]
                idx = i
        active[idx], active[-1] = active[-1], active[idx]
        u = active.pop()

        # Compute the eccentricity of u and update eccentricity lower bounds
        ecc_u = simple_BFS(g, u, distances, NULL, waiting_list, seen)
        LB = max(LB, ecc_u)
        if LB == UINT32_MAX:  # Disconnected graph
            break

        for v in active:
            ecc_lower_bound[v] = max(ecc_lower_bound[v], distances[v])

        # 2. Select x such that dist(u, x) + ecc[x] == ecc[u].
        # Since we don't know ecc[x], we select x with minimum eccentricity
        # lower bound.  If ecc[x] == ecc_lb[x], we are done. Otherwise, we
        # update eccentricity lower bounds and repeat
        while active:
            # Select v with minimum eccentricity lower bound
            tmp = UINT32_MAX
            for i, v in enumerate(active):
                if ecc_lower_bound[v] < tmp:
                    tmp = ecc_lower_bound[v]
                    idx = i
            active[idx], active[-1] = active[-1], active[idx]
            x = active.pop()
            ecc_x = simple_BFS(g, x, distances, NULL, waiting_list, seen)
            LB = max(LB, ecc_x)

            if ecc_x == ecc_lower_bound[x]:
                # We found the good vertex x
                # We update eccentricity upper bounds of the remaining vertices,
                # set UB to the largest of these values and break
                UB = 0
                for v in active:
                    ecc_upper_bound[v] = min(ecc_upper_bound[v], distances[v] + ecc_x)
                    UB = max(UB, ecc_upper_bound[v])
                break

            else:
                # x was not a good choice
                # We use its antipode to update eccentricity lower bounds.
                # Observe that this antipode might have already been seen.
                antipode = waiting_list[n-1]
                for i, v in enumerate(active):
                    if v == antipode:
                        active[i] = active[-1]
                        tmp = active.pop()
                        break

                ecc_antipode = simple_BFS(g, antipode, distances, NULL, waiting_list, seen)
                LB = max(LB, ecc_antipode)
                for v in active:
                    ecc_lower_bound[v] = max(ecc_lower_bound[v], distances[v])

    bitset_free(seen)
    return LB


def diameter(G, algorithm=None, source=None):
    r"""
    Return the diameter of `G`.

    This method returns Infinity if the (di)graph is not connected. It can
    also quickly return a lower bound on the diameter using the ``2sweep``,
    ``2Dsweep`` and ``multi-sweep`` schemes.

    INPUT:

    - ``algorithm`` -- string (default: ``None``); specifies the algorithm to
      use among:

      - ``'standard'`` -- computes the diameter of the input (di)graph as the
        largest eccentricity of its vertices. This is the classical algorithm
        with time complexity in `O(nm)`.

      - ``'2sweep'`` -- computes a lower bound on the diameter of an
        unweighted undirected graph using 2 BFS, as proposed in [MLH2008]_.  It
        first selects a vertex `v` that is at largest distance from an initial
        vertex source using BFS. Then it performs a second BFS from `v`. The
        largest distance from `v` is returned as a lower bound on the diameter
        of `G`.  The time complexity of this algorithm is linear in the size of
        `G`.

      - ``'2Dsweep'`` -- computes lower bound on the diameter of an unweighted
        directed graph using directed version of ``2sweep`` as proposed in
        [Broder2000]_. If the digraph is not strongly connected, the returned
        value is infinity.

      - ``'DHV'`` -- computes diameter of unweighted undirected graph using the
        algorithm proposed in [Dragan2018]_

      - ``'multi-sweep'`` -- computes a lower bound on the diameter of an
        unweighted undirected graph using several iterations of the ``2sweep``
        algorithms [CGHLM2013]_. Roughly, it first uses ``2sweep`` to identify
        two vertices `u` and `v` that are far apart. Then it selects a vertex
        `w` that is at same distance from `u` and `v`.  This vertex `w` will
        serve as the new source for another iteration of the ``2sweep``
        algorithm that may improve the current lower bound on the diameter.
        This process is repeated as long as the lower bound on the diameter
        is improved.

      - ``'iFUB'`` -- the iFUB (iterative Fringe Upper Bound) algorithm,
        proposed in [CGILM2010]_, computes the exact value of the diameter of an
        unweighted undirected graph. It is based on the following observation:

            The diameter of the graph is equal to the maximum eccentricity of
            a vertex. Let `v` be any vertex, and let `V` be partitioned into
            `A\cup B` where:

            .. MATH::

                d(v,a) \leq i, \forall a \in A\\
                d(v,b) \geq i, \forall b \in B

            As all vertices from `A` are at distance `\leq 2i` from each
            other, a vertex `a\in A` with eccentricity `ecc(a)>2i` is at
            distance `ecc(a)` from some vertex `b\in B`.

            Consequently, if we have already computed the maximum eccentricity
            `m` of all vertices in `B` and if `m>2i`, then we do not need to
            compute the eccentricity of the vertices in `A`.

        Starting from a vertex `v` obtained through a multi-sweep computation
        (which refines the 4sweep algorithm used in [CGHLM2013]_), we compute
        the diameter by computing the eccentricity of all vertices sorted
        decreasingly according to their distance to `v`, and stop as allowed
        by the remark above. The worst case time complexity of the iFUB
        algorithm is `O(nm)`, but it can be very fast in practice.

      - ``'DiFUB'`` -- the directed version of iFUB (iterative Fringe Upper
        Bound) algorithm. See the code's documentation and [CGLM2012]_ for more
        details. If the digraph is not strongly connected, the returned value is
        infinity.

    - ``source`` -- (default: ``None``) vertex from which to start the first BFS.
      If ``source==None``, an arbitrary vertex of the graph is chosen. Raise an
      error if the initial vertex is not in `G`.  This parameter is not used
      when ``algorithm=='standard'``.

    .. NOTE::
        As the graph is first converted to a short_digraph, all complexity
        have an extra `O(m+n)` for ``SparseGraph`` and `O(n^2)` for
        ``DenseGraph``.

    EXAMPLES::

        sage: from sage.graphs.distances_all_pairs import diameter
        sage: G = graphs.PetersenGraph()
        sage: diameter(G, algorithm='iFUB')
        2
        sage: G = Graph({0: [], 1: [], 2: [1]})
        sage: diameter(G, algorithm='iFUB')
        +Infinity
        sage: G = digraphs.Circuit(6)
        sage: diameter(G, algorithm='2Dsweep')
        5
        sage: G = graphs.PathGraph(7).to_directed()
        sage: diameter(G, algorithm='DiFUB')
        6

    Although max( ) is usually defined as -Infinity, since the diameter will
    never be negative, we define it to be zero::

        sage: G = graphs.EmptyGraph()
        sage: diameter(G, algorithm='iFUB')
        0

    Comparison of exact algorithms for graphs::

        sage: # needs networkx
        sage: G = graphs.RandomBarabasiAlbert(100, 2)
        sage: d1 = diameter(G, algorithm='standard')
        sage: d2 = diameter(G, algorithm='iFUB')
        sage: d3 = diameter(G, algorithm='iFUB', source=G.random_vertex())
        sage: d4 = diameter(G, algorithm='DHV')
        sage: if d1 != d2 or d1 != d3 or d1 != d4: print("Something goes wrong!")

    Comparison of lower bound algorithms::

        sage: lb2 = diameter(G, algorithm='2sweep')                                     # needs networkx
        sage: lbm = diameter(G, algorithm='multi-sweep')                                # needs networkx
        sage: if not (lb2 <= lbm and lbm <= d3): print("Something goes wrong!")         # needs networkx

    Comparison of exact algorithms for digraphs::

        sage: # needs networkx
        sage: D = DiGraph(graphs.RandomBarabasiAlbert(50, 2))
        sage: d1 = diameter(D, algorithm='standard')
        sage: d2 = diameter(D, algorithm='DiFUB')
        sage: d3 = diameter(D, algorithm='DiFUB', source=D.random_vertex())
        sage: d1 == d2 and d1 == d3
        True

    TESTS:

    This was causing a segfault. Fixed in :issue:`17873` ::

        sage: G = graphs.PathGraph(1)
        sage: diameter(G, algorithm='iFUB')
        0
    """
    cdef uint32_t n = G.order()
    if n <= 1:
        return 0

    if G.is_directed():
        if algorithm is None:
            algorithm = 'DiFUB'
        elif algorithm not in ['2Dsweep', 'standard', 'DiFUB']:
            raise ValueError("unknown algorithm for computing the diameter of directed graph")
    else:
        if algorithm is None:
            algorithm = 'iFUB'
        elif algorithm not in ['2sweep', 'multi-sweep', 'iFUB', 'standard', 'DHV']:
            raise ValueError("unknown algorithm for computing the diameter of undirected graph")

    if algorithm == 'standard':
        return max(G.eccentricity())
    if source is None:
        source = next(G.vertex_iterator())
    elif not G.has_vertex(source):
        raise ValueError("the specified source is not a vertex of the input Graph")

    # Copying the whole graph to obtain the list of neighbors quicker than by
    # calling out_neighbors. This data structure is well documented in the
    # module sage.graphs.base.static_sparse_graph
    cdef list int_to_vertex = list(G)
    cdef short_digraph sd
    init_short_digraph(sd, G, edge_labelled=False, vertex_list=int_to_vertex)
    cdef short_digraph rev_sd  # to store copy of sd with edges reversed

    # and we map the source to an int in [0,n-1]
    cdef uint32_t isource = 0 if source is None else int_to_vertex.index(source)

    cdef bitset_t seen
    cdef uint32_t* tab
    cdef uint32_t LB

    if algorithm == '2sweep':
        # We need to allocate arrays and bitset
        bitset_init(seen, n)
        tab = <uint32_t*> sig_malloc(2* n * sizeof(uint32_t))
        if not tab:
            free_short_digraph(sd)
            bitset_free(seen)
            raise MemoryError()

        LB = diameter_lower_bound_2sweep(sd, isource, tab, NULL, tab + n, seen)

        bitset_free(seen)
        sig_free(tab)

    elif algorithm == '2Dsweep':
        init_reverse(rev_sd, sd)
        LB = diameter_lower_bound_2Dsweep(sd, rev_sd, isource)[0]
        free_short_digraph(rev_sd)

    elif algorithm == 'multi-sweep':
        LB = diameter_lower_bound_multi_sweep(sd, isource)[0]

    elif algorithm == 'DiFUB':
        LB = diameter_DiFUB(sd, isource)

    elif algorithm == 'DHV':
        LB = diameter_DHV(sd)

    else:  # algorithm == 'iFUB'
        LB = diameter_iFUB(sd, isource)

    free_short_digraph(sd)

    if LB < 0 or LB > n:
        from sage.rings.infinity import Infinity
        return +Infinity
    return int(LB)


###########
# Radius #
###########

def radius_DHV(G):
    r"""
    Return the radius of unweighted graph `G`.

    This method computes the radius of unweighted undirected graph using the
    algorithm given in [Dragan2018]_.

    This method returns Infinity if graph is not connected.

    EXAMPLES::

        sage: from sage.graphs.distances_all_pairs import radius_DHV
        sage: G = graphs.PetersenGraph()
        sage: radius_DHV(G)
        2
        sage: G = graphs.RandomGNP(20,0.3)
        sage: from sage.graphs.distances_all_pairs import eccentricity
        sage: radius_DHV(G) == min(eccentricity(G, algorithm='bounds'))
        True

    TESTS:

        sage: G = Graph()
        sage: radius_DHV(G)
        0
        sage: G = Graph(1)
        sage: radius_DHV(G)
        0
        sage: G = Graph(2)
        sage: radius_DHV(G)
        +Infinity
        sage: G = graphs.PathGraph(2)
        sage: radius_DHV(G)
        1
        sage: G = DiGraph(1)
        sage: radius_DHV(G)
        Traceback (most recent call last):
        ...
        TypeError: this method works for unweighted undirected graphs only
    """
    if G.is_directed():
        raise TypeError("this method works for unweighted undirected graphs only")

    cdef uint32_t n = G.order()
    if n <= 1:
        return 0

    cdef list int_to_vertex
    cdef StaticSparseCGraph cg
    cdef short_digraph sd
    if isinstance(G, StaticSparseBackend):
        cg = <StaticSparseCGraph> G._cg
        sd = <short_digraph> cg.g
        int_to_vertex = cg._vertex_to_labels
    else:
        int_to_vertex = list(G)
        init_short_digraph(sd, G, edge_labelled=False, vertex_list=int_to_vertex)

    cdef uint32_t source, ecc_source
    cdef uint32_t antipode, ecc_antipode
    cdef uint32_t UB = UINT32_MAX
    cdef uint32_t LB = 0

    cdef MemoryAllocator mem = MemoryAllocator()
    cdef uint32_t * distances = <uint32_t *>mem.malloc(3 * n * sizeof(uint32_t))
    if not distances:
        raise MemoryError()

    cdef uint32_t * waiting_list = distances + n

    # For storing lower bound on eccentricity of nodes
    cdef uint32_t * ecc_lower_bound = distances + 2 * n
    memset(ecc_lower_bound, 0, n * sizeof(uint32_t))

    cdef bitset_t seen
    bitset_init(seen, n)

    # Algorithm
    source = 0
    while LB < UB:
        # 1) pick vertex with minimum eccentricity lower bound
        # and compute its eccentricity
        ecc_source = simple_BFS(sd, source, distances, NULL, waiting_list, seen)

        if ecc_source == UINT32_MAX:  # Disconnected graph
            break

        UB = min(UB, ecc_source)  # minimum among exact computed eccentricities
        if ecc_source == ecc_lower_bound[source]:
            # we have found minimum eccentricity vertex and hence the radius
            break

        # 2) Take vertex at largest distance from source, called antipode (last
        # vertex visited in simple_BFS), and compute its BFS distances.
        # By definition of antipode, we have ecc_antipode >= ecc_source.
        antipode = waiting_list[n-1]
        ecc_antipode = simple_BFS(sd, antipode, distances, NULL, waiting_list, seen)

        # 3) Use distances from antipode to improve eccentricity lower bounds.
        # We also determine the next source
        LB = UINT32_MAX
        for v in range(n):
            ecc_lower_bound[v] = max(ecc_lower_bound[v], distances[v])
            if LB > ecc_lower_bound[v]:
                LB = ecc_lower_bound[v]
                source = v  # vertex with minimum eccentricity lower bound

    if not isinstance(G, StaticSparseBackend):
        free_short_digraph(sd)
    bitset_free(seen)
    if UB == UINT32_MAX:
        from sage.rings.infinity import Infinity
        return +Infinity

    return UB


################
# Wiener index #
################

def wiener_index(G):
    r"""
    Return the Wiener index of the graph.

    The Wiener index of an undirected graph `G` is defined as
    `W(G) = \frac{1}{2} \sum_{u,v\in G} d(u,v)` where `d(u,v)` denotes the
    distance between vertices `u` and `v` (see [KRG1996]_).

    The Wiener index of a directed graph `G` is defined as the sum of the
    distances between each pairs of vertices, `W(G) = \sum_{u,v\in G} d(u,v)`.

    EXAMPLES:

    From [GYLL1993]_, cited in [KRG1996]_::

        sage: g=graphs.PathGraph(10)
        sage: w=lambda x: (x*(x*x -1)/6)
        sage: g.wiener_index()==w(10)
        True

    Wiener index of complete (di)graphs::

        sage: n = 5
        sage: g = graphs.CompleteGraph(n)
        sage: g.wiener_index() == (n * (n - 1)) / 2
        True
        sage: g = digraphs.Complete(n)
        sage: g.wiener_index() == n * (n - 1)
        True

    Wiener index of a graph of order 1::

        sage: Graph(1).wiener_index()
        0

    The Wiener index is not defined on the empty graph::

        sage: Graph().wiener_index()
        Traceback (most recent call last):
        ...
        ValueError: Wiener index is not defined for the empty graph
    """
    if not G:
        raise ValueError("Wiener index is not defined for the empty graph")

    cdef unsigned int n = G.order()
    if n == 1:
        return 0

    # Copying the whole graph to obtain the list of neighbors quicker than by
    # calling out_neighbors.  This data structure is well documented in the
    # module sage.graphs.base.static_sparse_graph
    cdef StaticSparseCGraph cg
    cdef short_digraph sd
    if isinstance(G, StaticSparseBackend):
        cg = <StaticSparseCGraph> G._cg
        sd = <short_digraph> cg.g
    else:
        init_short_digraph(sd, G, edge_labelled=False, vertex_list=list(G))

    # allocated some data structures
    cdef bitset_t seen
    bitset_init(seen, n)
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef uint32_t * distances = <uint32_t *> mem.allocarray(2 * n, sizeof(uint32_t))
    cdef uint32_t * waiting_list = distances + n

    cdef uint64_t s = 0
    cdef uint32_t u, v
    cdef uint32_t ecc
    for u in range(n):
        ecc = simple_BFS(sd, u, distances, NULL, waiting_list, seen)
        if ecc == UINT32_MAX:
            # the graph is not  connected
            s = UINT64_MAX
            break
        for v in range(0 if G.is_directed() else (u + 1), n):
            s += distances[v]

    if not isinstance(G, StaticSparseBackend):
        free_short_digraph(sd)
    bitset_free(seen)

    if s == UINT64_MAX:
        from sage.rings.infinity import Infinity
        return +Infinity
    return s


################
# Szeged index #
################

cdef uint64_t c_szeged_index_low_memory(short_digraph sd) noexcept:
    r"""
    Return the Szeged index of the graph.

    Let `G = (V, E)` be a connected simple graph, and for any `uv\in E`, let
    `N_u(uv) = \{w\in V:d(u,w)<d(v,w)\}` and `n_u(uv)=|N_u(uv)|`. The Szeged
    index of `G` is then defined as [KRG1996]_ as `\sum_{uv \in
    E(G)}n_u(uv)\times n_v(uv)`.

    To determine `N_u(uv)`, this method performs a breadth first search (BFS)
    from each vertex `s \in V`. Then, each time an edge `uv` visited by the BFS
    is such that `d(s, u) < d(s, v)`, it adds 1 to `N_u(uv)`. Since this method
    assumes that the graph is undirected, the graph `sd` has both arcs `uv` and
    `vu`. Using one counter per arc, the counter for arc `uv` records the number
    of vertices that are closer to the side `u` of edge `uv`, and the counter
    for arc `vu` records the number of vertices that are closer to the side `v`
    of edge `uv`.

    This method assumes that the input graph has no loops or multiple edges.

    EXAMPLES::

        sage: graphs.CycleGraph(4).szeged_index(algorithm='low')
        16
    """
    cdef size_t n = sd.n
    if n <= 1:
        return 0
    if n == 2:
        return 1

    cdef MemoryAllocator mem = MemoryAllocator()
    cdef uint32_t * current_layer = <uint32_t *> mem.malloc(n * sizeof(uint32_t))
    cdef uint32_t n_current
    cdef uint32_t * next_layer = <uint32_t *> mem.malloc(n * sizeof(uint32_t))
    cdef uint32_t n_next
    cdef uint32_t * seen = <uint32_t *> mem.calloc(n, sizeof(uint32_t))
    cdef uint32_t seen_value = 0

    # For each edge e = uv, we have 2 arcs. Let p_uv be a pointer to arc uv and
    # p_vu a pointer to arc vu. The index of arc uv is p_uv - sd.edges. During a
    # BFS from source, we add 1 to counter[p_uv - sd.edges] if vertex u is closer
    # from source than v and 1 to counter[p_vu - sd.edges] if vertex v is closer
    # from source than u. The Szeged index is then
    #   sum_{e=uv} counter[p_uv - sd.edges] * counter[p_vu - sd.edges]
    cdef uint32_t * counter = <uint32_t *> mem.calloc(2 * sd.m, sizeof(uint32_t))

    cdef uint32_t source, u, v, i
    cdef uint32_t* p_uv
    cdef uint32_t* p_end

    sig_on()
    for source in range(n):

        next_layer[0] = source
        n_next = 1
        seen_value += 1

        while n_next:
            # Go to next layer
            current_layer, next_layer = next_layer, current_layer
            n_current, n_next = n_next, 0

            # Mark all vertices in current layer as seen
            for i in range(n_current):
                seen[current_layer[i]] = seen_value

            for i in range(n_current):
                u = current_layer[i]

                # Visit all (out) neighbors of u
                p_uv = sd.neighbors[u]
                p_end = sd.neighbors[u + 1]
                while p_uv < p_end:
                    v = p_uv[0]
                    if seen[v] != seen_value:
                        # u is closer to the source
                        counter[p_uv - sd.edges] += 1

                        # Ensure that v is added only once for next_level
                        if seen[v] != seen_value + 1:
                            next_layer[n_next] = v
                            n_next += 1
                            seen[v] = seen_value + 1

                    p_uv += 1
    sig_off()

    cdef uint64_t s = 0
    cdef uint32_t* p_vu

    sig_on()
    for u in range(n - 1):
        p_uv = sd.neighbors[u]
        p_end = sd.neighbors[u + 1]
        while p_uv < p_end:
            v = p_uv[0]
            if u < v:
                # Get the pointer to arc vu
                p_vu = has_edge(sd, v, u)
                s += counter[p_uv - sd.edges] * counter[p_vu - sd.edges]

            p_uv += 1
    sig_off()

    return s


cdef uint64_t c_szeged_index_high_memory(short_digraph sd) noexcept:
    r"""
    Return the Szeged index of the graph.

    Let `G = (V, E)` be a connected graph, and for any `uv\in E`, let `N_u(uv) =
    \{w\in V:d(u,w)<d(v,w)\}` and `n_u(uv)=|N_u(uv)|`. The Szeged index of `G`
    is then defined as [KRG1996]_ as `\sum_{uv \in E(G)}n_u(uv)\times n_v(uv)`.

    EXAMPLES::

        sage: graphs.CycleGraph(4).szeged_index(algorithm='high')
        16
    """
    cdef int n = sd.n
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef unsigned short* distances = <unsigned short*> mem.malloc(n * n * sizeof(unsigned short))

    # Compute all pairs shortest path
    c_all_pairs_shortest_path_BFS(sd, NULL, distances, NULL)

    cdef uint32_t* p_uv
    cdef uint32_t* p_end
    cdef uint32_t u, v, w
    cdef unsigned short* du
    cdef unsigned short* dv
    cdef uint32_t n1, n2
    cdef uint64_t s = 0

    for u in range(n):
        du = distances + u * n
        p_uv = sd.neighbors[u]
        p_end = sd.neighbors[u + 1]
        while p_uv < p_end:
            v = p_uv[0]
            if u < v:
                dv = distances + v * n
                n1 = n2 = 0
                for w in range(n):
                    if du[w] < dv[w]:
                        n1 += 1
                    elif dv[w] < du[w]:
                        n2 += 1

                s += n1 * n2
            p_uv += 1

    return s


def szeged_index(G, algorithm=None):
    r"""
    Return the Szeged index of the graph `G`.

    Let `G = (V, E)` be a connected graph, and for any `uv\in E`, let `N_u(uv) =
    \{w\in V:d(u,w)<d(v,w)\}` and `n_u(uv)=|N_u(uv)|`. The Szeged index of `G`
    is then defined as [KRG1996]_

    .. MATH::

        `\sum_{uv \in E(G)}n_u(uv)\times n_v(uv)`

    See the :wikipedia:`Szeged_index` for more details.

    INPUT:

    - ``G`` -- a Sage graph

    - ``algorithm`` -- string (default: ``None``); algorithm to use among:

      - ``'low'`` -- algorithm with time complexity in `O(nm)` and space
        complexity in `O(m)`. This implementation is currently valid only for
        simple (without loops or multiple edges) connected graphs.

      - ``'high'`` -- algorithm with time complexity in `O(nm)` and space
        complexity in `O(n^2)`. It cannot be used on graphs with more than
        `65536 = 2^{16}` vertices.

      By default (``None``), the ``'low'`` algorithm is used for graphs and the
      ``'high'`` algorithm for digraphs.

    EXAMPLES:

    True for any connected graph [KRG1996]_::

        sage: from sage.graphs.distances_all_pairs import szeged_index
        sage: g = graphs.PetersenGraph()
        sage: g.wiener_index() <= szeged_index(g)
        True

    True for all trees [KRG1996]_::

        sage: g = Graph()
        sage: g.add_edges(graphs.CubeGraph(5).min_spanning_tree())
        sage: g.wiener_index() == szeged_index(g)
        True

    Check that both algorithms return same value::

        sage: # long time, needs networkx
        sage: G = graphs.RandomBarabasiAlbert(100, 2)
        sage: a = szeged_index(G, algorithm='low')
        sage: b = szeged_index(G, algorithm='high')
        sage: a == b
        True

    The Szeged index of a directed circuit of order `n` is `(n-1)^2`::

        sage: [digraphs.Circuit(n).szeged_index() for n in range(1, 8)]
        [0, 1, 4, 9, 16, 25, 36]

    TESTS:

    Not defined when the graph is not connected (:issue:`26803`)::

        sage: szeged_index(Graph({0: [1], 2: []}))
        Traceback (most recent call last):
        ...
        ValueError: the Szeged index is defined for connected graphs only

    Directed graphs must be strongly connected::

        sage: szeged_index(digraphs.Path(2))
        Traceback (most recent call last):
        ...
        ValueError: the Szeged index is defined for strongly connected digraphs only

    Wrong name of algorithm::

        sage: szeged_index(Graph(1), algorithm='wheel')
        Traceback (most recent call last):
        ...
        ValueError: unknown algorithm 'wheel'

    Algorithm `"low"` is for graphs without loops or multiple edges::

        sage: szeged_index(Graph([(0, 0)], loops=True), algorithm='low')
        Traceback (most recent call last):
        ...
        ValueError: the 'low' algorithm is for simple connected undirected graphs only
        sage: szeged_index(Graph([(0, 1), (0, 1)], multiedges=True), algorithm='low')
        Traceback (most recent call last):
        ...
        ValueError: the 'low' algorithm is for simple connected undirected graphs only
        sage: szeged_index(digraphs.Circuit(3), algorithm='low')
        Traceback (most recent call last):
        ...
        ValueError: the 'low' algorithm cannot be used on digraphs
    """
    if not G.is_connected():
        raise ValueError("the Szeged index is defined for connected graphs only")
    if G.is_directed() and not G.is_strongly_connected():
        raise ValueError("the Szeged index is defined for "
                         "strongly connected digraphs only")
    if G.is_directed() and algorithm == "low":
        raise ValueError("the 'low' algorithm cannot be used on digraphs")

    if algorithm is None:
        algorithm = "high" if G.is_directed() else "low"

    elif algorithm not in ["low", "high"]:
        raise ValueError(f"unknown algorithm '{algorithm}'")

    if algorithm == "low" and (G.has_loops() or G.has_multiple_edges()):
        raise ValueError("the 'low' algorithm is for simple connected "
                         "undirected graphs only")

    if G.order() <= 1:
        return 0

    cdef StaticSparseCGraph cg
    cdef short_digraph sd
    if isinstance(G, StaticSparseBackend):
        cg = <StaticSparseCGraph> G._cg
        sd = <short_digraph> cg.g
    else:
        init_short_digraph(sd, G, edge_labelled=False, vertex_list=list(G))

    cdef uint64_t s

    if algorithm == "low":
        s = c_szeged_index_low_memory(sd)
    else:
        s = c_szeged_index_high_memory(sd)

    if not isinstance(G, StaticSparseBackend):
        free_short_digraph(sd)
    return s


##########################
# Distances distribution #
##########################

def distances_distribution(G):
    r"""
    Return the distances distribution of the (di)graph in a dictionary.

    This method *ignores all edge labels*, so that the distance considered is
    the topological distance.

    OUTPUT:

        A dictionary ``d`` such that the number of pairs of vertices at distance
        ``k`` (if any) is equal to `d[k] \cdot |V(G)| \cdot (|V(G)|-1)`.

    .. NOTE::

        We consider that two vertices that do not belong to the same connected
        component are at infinite distance, and we do not take the trivial pairs
        of vertices `(v, v)` at distance `0` into account. Empty (di)graphs and
        (di)graphs of order 1 have no paths and so we return the empty
        dictionary ``{}``.

    EXAMPLES:

    An empty Graph::

        sage: g = Graph()
        sage: g.distances_distribution()
        {}

    A Graph of order 1::

        sage: g = Graph()
        sage: g.add_vertex(1)
        sage: g.distances_distribution()
        {}

    A Graph of order 2 without edge::

        sage: g = Graph()
        sage: g.add_vertices([1,2])
        sage: g.distances_distribution()
        {+Infinity: 1}

    The Petersen Graph::

        sage: g = graphs.PetersenGraph()
        sage: g.distances_distribution()
        {1: 1/3, 2: 2/3}

    A graph with multiple disconnected components::

        sage: g = graphs.PetersenGraph()
        sage: g.add_edge('good','wine')
        sage: g.distances_distribution()
        {1: 8/33, 2: 5/11, +Infinity: 10/33}

    The de Bruijn digraph dB(2,3)::

        sage: D = digraphs.DeBruijn(2,3)                                                # needs sage.combinat
        sage: D.distances_distribution()                                                # needs sage.combinat
        {1: 1/4, 2: 11/28, 3: 5/14}
    """
    cdef size_t n = G.order()
    if n <= 1:
        return {}

    cdef StaticSparseCGraph cg
    cdef short_digraph sd
    if isinstance(G, StaticSparseBackend):
        cg = <StaticSparseCGraph> G._cg
        sd = <short_digraph> cg.g
    else:
        init_short_digraph(sd, G, edge_labelled=False, vertex_list=list(G))

    cdef MemoryAllocator mem = MemoryAllocator()
    cdef uint32_t * distances = <uint32_t *> mem.allocarray(2 * n, sizeof(uint32_t))
    cdef uint32_t * waiting_list = distances + n
    cdef uint64_t * count = <uint64_t *> mem.calloc(n, sizeof(uint64_t))
    cdef bitset_t seen
    bitset_init(seen, n)

    # We count the number of pairs at equal distances
    cdef uint32_t u, v
    cdef uint64_t count_inf = 0
    for u in range(n):
        ecc = simple_BFS(sd, u, distances, NULL, waiting_list, seen)
        if ecc == UINT32_MAX:
            for v in range(n):
                if bitset_in(seen, v):
                    count[distances[v]] += 1
            count_inf += n - bitset_len(seen)
        else:
            for v in range(n):
                count[distances[v]] += 1

    if not isinstance(G, StaticSparseBackend):
        free_short_digraph(sd)
    bitset_free(seen)

    from sage.rings.infinity import Infinity
    from sage.rings.rational_field import QQ

    # We normalize the distribution
    cdef uint64_t NN = n * (n - 1)
    cdef dict distr = {+Infinity: QQ((count_inf, NN))} if count_inf else {}
    cdef size_t d
    for d in range(1, n):
        if count[d]:
            distr[d] = QQ((count[d], NN))

    return distr


###################
# Antipodal graph #
###################

def antipodal_graph(G):
    r"""
    Return the antipodal graph of `G`.

    The antipodal graph of a graph `G` has the same vertex set of `G` and
    two vertices are adjacent if their distance in `G` is equal to the
    diameter of `G`.

    This method first computes the eccentricity of all vertices and determines
    the diameter of the graph. Then, it for each vertex `u` with eccentricity
    the diameter, it computes BFS distances from `u` and add an edge in the
    antipodal graph for each vertex `v` at diameter distance from `u`
    (i.e., for each antipodal vertex).

    The drawback of this method is that some BFS distances may be computed
    twice, one time to determine the eccentricities and another time is the
    vertex has eccentricity equal to the diameter. However, in practive, this is
    much more efficient. See the documentation of method
    :meth:`c_eccentricity_DHV`.

    EXAMPLES:

    The antipodal graph of a grid graph has only 2 edges::

        sage: from sage.graphs.distances_all_pairs import antipodal_graph
        sage: G = graphs.Grid2dGraph(5, 5)
        sage: A = antipodal_graph(G)
        sage: A.order(), A.size()
        (25, 2)

    The antipodal graph of a disjoint union of cliques is its complement::

        sage: from sage.graphs.distances_all_pairs import antipodal_graph
        sage: G = graphs.CompleteGraph(3) * 3
        sage: A = antipodal_graph(G)
        sage: A.is_isomorphic(G.complement())
        True

    The antipodal graph can also be constructed as the
    :meth:`sage.graphs.generic_graph.distance_graph` for diameter distance::

        sage: from sage.graphs.distances_all_pairs import antipodal_graph
        sage: G = graphs.RandomGNP(10, .2)
        sage: A = antipodal_graph(G)
        sage: B = G.distance_graph(G.diameter())
        sage: A.is_isomorphic(B)
        True

    TESTS::

        sage: from sage.graphs.distances_all_pairs import antipodal_graph
        sage: antipodal_graph(Graph())
        Traceback (most recent call last):
        ...
        ValueError: the antipodal graph of the empty graph is not defined
        sage: antipodal_graph(DiGraph(1))
        Traceback (most recent call last):
        ...
        ValueError: this method is defined for undirected graphs only
        sage: antipodal_graph(Graph(1))
        Antipodal graph of Graph on 1 vertex: Looped graph on 1 vertex
        sage: antipodal_graph(Graph(2)).edges(sort=True, labels=False)
        [(0, 1)]
    """
    if not G:
        raise ValueError("the antipodal graph of the empty graph is not defined")
    if G.is_directed():
        raise ValueError("this method is defined for undirected graphs only")

    from sage.graphs.graph import Graph

    cdef uint32_t n = G.order()
    name = f"Antipodal graph of {G}"
    if n == 1:
        return Graph(list(zip(G, G)), loops=True, name=name)

    import copy
    A = Graph(name=name, pos=copy.deepcopy(G.get_pos()))

    if not G.is_connected():
        import itertools
        CC = G.connected_components(sort=False)
        for c1, c2 in itertools.combinations(CC, 2):
            A.add_edges(itertools.product(c1, c2))
        return A

    cdef list int_to_vertex
    cdef StaticSparseCGraph cg
    cdef short_digraph sd
    if isinstance(G, StaticSparseBackend):
        cg = <StaticSparseCGraph> G._cg
        sd = <short_digraph> cg.g
        int_to_vertex = cg._vertex_to_labels
    else:
        int_to_vertex = list(G)
        init_short_digraph(sd, G, edge_labelled=False, vertex_list=int_to_vertex)

    cdef MemoryAllocator mem = MemoryAllocator()
    cdef uint32_t * distances = <uint32_t *> mem.allocarray(2 * n, sizeof(uint32_t))
    cdef uint32_t * waiting_list = distances + n
    cdef bitset_t seen
    bitset_init(seen, n)

    # Get the eccentricity of all vertices
    cdef uint32_t* ecc = c_eccentricity_DHV(sd)
    cdef uint32_t i
    cdef uint32_t diam = 0
    for i in range(n):
        if ecc[i] > diam:
            diam = ecc[i]

    cdef uint32_t ui, vj, j
    for ui in range(n):
        if ecc[ui] == diam:
            _ = simple_BFS(sd, ui, distances, NULL, waiting_list, seen)
            u = int_to_vertex[ui]
            j = n - 1
            while distances[waiting_list[j]] == diam:
                vj = waiting_list[j]
                if ui < vj:  # avoid adding twice the same edge
                    A.add_edge(u, int_to_vertex[vj])
                j -= 1

    if not isinstance(G, StaticSparseBackend):
        free_short_digraph(sd)
    bitset_free(seen)

    A.add_vertices(G)
    return A


##################
# Floyd-Warshall #
##################

def floyd_warshall(gg, paths=True, distances=False):
    r"""
    Compute the shortest path/distances between all pairs of vertices.

    For more information on the Floyd-Warshall algorithm, see
    the :wikipedia:`Floyd-Warshall_algorithm`.

    INPUT:

    - ``gg`` -- the graph on which to work

    - ``paths`` -- boolean (default: ``True``); whether to return the dictionary
      of shortest paths

    - ``distances`` -- boolean (default: ``False``); whether to return the
      dictionary of distances

    OUTPUT:

    Depending on the input, this function return the dictionary of paths, the
    dictionary of distances, or a pair of dictionaries ``(distances, paths)``
    where ``distance[u][v]`` denotes the distance of a shortest path from `u` to
    `v` and ``paths[u][v]`` denotes an inneighbor `w` of `v` such that
    `dist(u,v) = 1 + dist(u,w)`.

    .. WARNING::

        Because this function works on matrices whose size is quadratic compared
        to the number of vertices, it uses short variables instead of long ones
        to divide by 2 the size in memory. This means that the current
        implementation does not run on a graph of more than 65536 nodes (this
        can be easily changed if necessary, but would require much more
        memory. It may be worth writing two versions). For information, the
        current version of the algorithm on a graph with `65536 = 2^{16}` nodes
        creates in memory `2` tables on `2^{32}` short elements (2bytes each),
        for a total of `2^{34}` bytes or `16` gigabytes. Let us also remember
        that if the memory size is quadratic, the algorithm runs in cubic time.

    .. NOTE::

        When ``paths = False`` the algorithm saves roughly half of the memory as
        it does not have to maintain the matrix of predecessors. However,
        setting ``distances=False`` produces no such effect as the algorithm can
        not run without computing them. They will not be returned, but they will
        be stored while the method is running.

    EXAMPLES:

    Shortest paths in a small grid ::

        sage: g = graphs.Grid2dGraph(2,2)
        sage: from sage.graphs.distances_all_pairs import floyd_warshall
        sage: print(floyd_warshall(g))
        {(0, 0): {(0, 0): None, (0, 1): (0, 0), (1, 0): (0, 0), (1, 1): (0, 1)},
         (0, 1): {(0, 1): None, (0, 0): (0, 1), (1, 0): (0, 0), (1, 1): (0, 1)},
         (1, 0): {(1, 0): None, (0, 0): (1, 0), (0, 1): (0, 0), (1, 1): (1, 0)},
         (1, 1): {(1, 1): None, (0, 0): (0, 1), (0, 1): (1, 1), (1, 0): (1, 1)}}

    Checking the distances are correct ::

        sage: g = graphs.Grid2dGraph(5,5)
        sage: dist,path = floyd_warshall(g, distances=True)
        sage: all(dist[u][v] == g.distance(u, v) for u in g for v in g)
        True

    Checking a random path is valid ::

        sage: u,v = g.random_vertex(), g.random_vertex()
        sage: p = [v]
        sage: while p[0] is not None:
        ....:   p.insert(0,path[u][p[0]])
        sage: len(p) == dist[u][v] + 2
        True

    Distances for all pairs of vertices in a diamond::

        sage: g = graphs.DiamondGraph()
        sage: floyd_warshall(g, paths=False, distances=True)
        {0: {0: 0, 1: 1, 2: 1, 3: 2},
         1: {0: 1, 1: 0, 2: 1, 3: 1},
         2: {0: 1, 1: 1, 2: 0, 3: 1},
         3: {0: 2, 1: 1, 2: 1, 3: 0}}

    TESTS:

    Too large graphs::

        sage: from sage.graphs.distances_all_pairs import floyd_warshall
        sage: floyd_warshall(Graph(65536))
        Traceback (most recent call last):
        ...
        ValueError: the graph backend contains more than 65535 nodes
    """
    from sage.rings.infinity import Infinity
    cdef CGraph g = <CGraph> gg._backend.c_graph()[0]

    cdef list gverts = g.verts()

    if not gverts:
        if distances and paths:
            return {}, {}
        return {}

    cdef unsigned int n = max(gverts) + 1

    if n >= <unsigned short> -1:
        raise ValueError("the graph backend contains more than " + str(<unsigned short> -1) + " nodes")

    # All this just creates two tables prec[n][n] and dist[n][n]
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef unsigned short* t_prec = NULL
    cdef unsigned short** prec = NULL
    # init dist
    cdef unsigned short* t_dist = <unsigned short*> mem.allocarray(n * n, sizeof(unsigned short))
    cdef unsigned short** dist = <unsigned short**> mem.allocarray(n, sizeof(unsigned short*))
    dist[0] = t_dist
    cdef unsigned int i
    for i in range(1, n):
        dist[i] = dist[i - 1] + n
    memset(t_dist, -1, n * n * sizeof(short))

    cdef int v_int
    cdef int u_int
    cdef int w_int

    # Copying the adjacency matrix (vertices at distance 1)
    for v_int in gverts:
        dist[v_int][v_int] = 0
        for u_int in g.out_neighbors(v_int):
            dist[v_int][u_int] = 1

    if paths:
        # init prec
        t_prec = <unsigned short*> mem.allocarray(n * n, sizeof(unsigned short))
        prec = <unsigned short**> mem.allocarray(n, sizeof(unsigned short*))
        prec[0] = t_prec
        for i in range(1, n):
            prec[i] = prec[i - 1] + n
        memset(t_prec, 0, n * n * sizeof(short))
        # Copying the adjacency matrix (vertices at distance 1)
        for v_int in gverts:
            prec[v_int][v_int] = v_int
            for u_int in g.out_neighbors(v_int):
                prec[v_int][u_int] = v_int

    # The algorithm itself.
    cdef unsigned short* dv
    cdef unsigned short* dw
    cdef int dvw
    cdef int val

    for w_int in gverts:
        dw = dist[w_int]
        for v_int in gverts:
            dv = dist[v_int]
            dvw = dv[w_int]
            for u_int in gverts:
                val = dvw + dw[u_int]
                # If it is shorter to go from u to v through w, do it
                if dv[u_int] > val:
                    dv[u_int] = val
                    if paths:
                        prec[v_int][u_int] = prec[w_int][u_int]

    # Dictionaries of distance, precedent element, and integers
    cdef dict d_prec = {}
    cdef dict d_dist = {}
    cdef dict tmp_prec
    cdef dict tmp_dist

    cdef CGraphBackend cgb = <CGraphBackend> gg._backend

    for v_int in gverts:
        v = cgb.vertex_label(v_int)
        if paths:
            tmp_prec = {v: None}
        if distances:
            tmp_dist = {v: 0}
        dv = dist[v_int]
        for u_int in gverts:
            u = cgb.vertex_label(u_int)
            if v != u and dv[u_int] != <unsigned short> -1:
                if paths:
                    tmp_prec[u] = cgb.vertex_label(prec[v_int][u_int])

                if distances:
                    tmp_dist[u] = dv[u_int]

        if paths:
            d_prec[v] = tmp_prec
        if distances:
            d_dist[v] = tmp_dist

    if distances and paths:
        return d_dist, d_prec
    if paths:
        return d_prec
    if distances:
        return d_dist
