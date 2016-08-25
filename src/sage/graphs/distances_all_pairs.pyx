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

**Memory cost** : The methods implemented in the current module sometimes need large
amounts of memory to return their result. Storing the distances between all
pairs of vertices in a graph on `1500` vertices as a dictionary of dictionaries
takes around 200MB, while storing the same information as a C array requires
4MB.


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
      encodes the information of a shortest `uv`-path for any `u,v\in G` :
      indeed, to go from `u` to `v` you should first find a shortest
      `uP[u,v]`-path, then jump from `P[u,v]` to `v` as it is one of its
      outneighbors. Apply recursively and find out what the whole path is !.

    - The matrix of distances.

      This matrix has size `n^2` and associates to any `uv` the distance from
      `u` to `v`.

    - The vector of eccentricities.

      This vector of size `n` encodes for each vertex `v` the distance to vertex
      which is furthest from `v` in the graph. In particular, the diameter of
      the graph is the maximum of these values.

**What does it take as input ?**

    - ``gg`` a (Di)Graph.

    - ``unsigned short * predecessors`` -- a pointer toward an array of size
      `n^2\cdot\text{sizeof(unsigned short)}`. Set to ``NULL`` if you do not
      want to compute the predecessors.

    - ``unsigned short * distances`` -- a pointer toward an array of size
      `n^2\cdot\text{sizeof(unsigned short)}`. The computation of the distances
      is necessary for the algorithm, so this value can **not** be set to
      ``NULL``.

    - ``int * eccentricity`` -- a pointer toward an array of size
      `n\cdot\text{sizeof(int)}`. Set to ``NULL`` if you do not want to compute
      the eccentricity.

**Technical details**

    - The vertices are encoded as `1, ..., n` as they appear in the ordering of
      ``G.vertices()``.

    - Because this function works on matrices whose size is quadratic compared
      to the number of vertices when computing all distances or predecessors, it
      uses short variables to store the vertices' names instead of long ones to
      divide by 2 the size in memory. This means that only the
      diameter/eccentricities can be computed on a graph of more than 65536
      nodes. For information, the current version of the algorithm on a graph
      with `65536=2^{16}` nodes creates in memory `2` tables on `2^{32}` short
      elements (2bytes each), for a total of `2^{33}` bytes or `8` gigabytes. In
      order to support larger sizes, we would have to replace shorts by 32-bits
      int or 64-bits int, which would then require respectively 16GB or 32GB.

    - In the C version of these functions, infinite distances are represented
      with ``<unsigned short> -1 = 65535`` for ``unsigned short`` variables, and
      by ``INT32_MAX`` otherwise. These case happens when the input is a
      disconnected graph, or a non-strongly-connected digraph.

    - A memory error is raised when data structures allocation failed. This
      could happen with large graphs on computers with low memory space.

    .. WARNING::

        The function ``all_pairs_shortest_path_BFS`` has **no reason** to be
        called by the user, even though he would be writing his code in Cython
        and look for efficiency. This module contains wrappers for this function
        that feed it with the good parameters. As the function is inlined, using
        those wrappers actually saves time as it should avoid testing the
        parameters again and again in the main function's body.

AUTHOR:

- Nathann Cohen (2011)
- David Coudert (2014) -- 2sweep, multi-sweep and iFUB for diameter computation

REFERENCE:

.. [KRG96b] \S. Klavzar, A. Rajapakse, and I. Gutman. The Szeged and the
  Wiener index of graphs. *Applied Mathematics Letters*, 9(5):45--49, 1996.

.. [GYLL93c] \I. Gutman, Y.-N. Yeh, S.-L. Lee, and Y.-L. Luo. Some recent
  results in the theory of the Wiener number. *Indian Journal of
  Chemistry*, 32A:651--661, 1993.

.. [CGH+13] \P. Crescenzi, R. Grossi, M. Habib, L. Lanzi, A. Marino. On computing
  the diameter of real-world undirected graphs. *Theor. Comput. Sci.* 514: 84-95
  (2013) http://dx.doi.org/10.1016/j.tcs.2012.09.018

.. [CGI+10] \P. Crescenzi, R. Grossi, C. Imbrenda, L. Lanzi, and A. Marino.
  Finding the Diameter in Real-World Graphs: Experimentally Turning a Lower
  Bound into an Upper Bound. Proceedings of *18th Annual European Symposium on
  Algorithms*. Lecture Notes in Computer Science, vol. 6346, 302-313. Springer
  (2010).

.. [MLH08] \C. Magnien, M. Latapy, and M. Habib. Fast computation of empirically
  tight bounds for the diameter of massive graphs. *ACM Journal of Experimental
  Algorithms* 13 (2008) http://dx.doi.org/10.1145/1412228.1455266

.. [TK13] \F. W. Takes and W. A. Kosters. Computing the eccentricity distribution
  of large graphs. *Algorithms* 6:100-118 (2013)
  http://dx.doi.org/10.3390/a6010100

Functions
---------
"""

#*****************************************************************************
#       Copyright (C) 2011 Nathann Cohen <nathann.cohen@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

include "sage/data_structures/binary_matrix.pxi"
from libc.string cimport memset
from libc.stdint cimport uint64_t, uint32_t, INT32_MAX, UINT32_MAX
from sage.graphs.base.c_graph cimport CGraphBackend
from sage.graphs.base.c_graph cimport CGraph
from sage.ext.memory_allocator cimport MemoryAllocator

from sage.graphs.base.static_sparse_graph cimport short_digraph, init_short_digraph, free_short_digraph, out_degree

from sage.misc.decorators import rename_keyword

cdef inline all_pairs_shortest_path_BFS(gg,
                                        unsigned short * predecessors,
                                        unsigned short * distances,
                                        uint32_t       * eccentricity):
    """
    See the module's documentation.
    """

    from sage.rings.infinity import Infinity

    cdef list int_to_vertex = gg.vertices()
    cdef int i
    cdef MemoryAllocator mem = MemoryAllocator()

    cdef int n = len(int_to_vertex)

    # Computing the predecessors/distances can only be done if we have less than
    # MAX_UNSIGNED_SHORT vertices. No problem with the eccentricities though as
    # we store them on an integer vector.
    if (predecessors != NULL or distances != NULL) and n > <unsigned short> -1:
        raise ValueError("The graph backend contains more than "+
                         str(<unsigned short> -1)+" nodes and we cannot "+
                         "compute the matrix of distances/predecessors on "+
                         "something like that !")

    # The vertices which have already been visited
    cdef bitset_t seen
    bitset_init(seen, n)

    # The list of waiting vertices, the beginning and the end of the list
    cdef int * waiting_list = <int *> mem.allocarray(n, sizeof(int))
    cdef int waiting_beginning = 0
    cdef int waiting_end = 0

    cdef int source
    cdef int v, u
    cdef uint32_t * p_tmp
    cdef uint32_t * end

    cdef unsigned short * c_predecessors = predecessors
    cdef int * c_distances = <int *> mem.allocarray(n, sizeof(int))

    # Copying the whole graph to obtain the list of neighbors quicker than by
    # calling out_neighbors

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

    cdef short_digraph sd
    init_short_digraph(sd, gg)
    cdef uint32_t ** p_vertices = sd.neighbors
    cdef uint32_t * p_edges = sd.edges
    cdef uint32_t * p_next = p_edges

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
            end = p_vertices[v+1]

            # Iterating over all the outneighbors u of v
            while p_tmp < end:
                u = p_tmp[0]

                # If we notice one of these neighbors is not seen yet, we set
                # its parameters and add it to the queue to be explored later.
                if not bitset_in(seen, u):
                    c_distances[u] = c_distances[v]+1
                    if predecessors != NULL:
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
                if predecessors != NULL:
                    c_predecessors[v] = -1
                v = bitset_next(seen, v+1)

            if eccentricity != NULL:
                eccentricity[source] = UINT32_MAX

        elif eccentricity != NULL:
            eccentricity[source] = c_distances[waiting_list[n-1]]


        if predecessors != NULL:
            c_predecessors += n

        if distances != NULL:
            for i in range(n):
                distances[i] = <unsigned short> c_distances[i]
            distances += n

    bitset_free(seen)
    free_short_digraph(sd)

################
# Predecessors #
################

cdef unsigned short * c_shortest_path_all_pairs(G) except NULL:
    r"""
    Returns the matrix of predecessors in G.

    The matrix `P` returned has size `n^2`, and is such that vertex `P[u,v]` is
    a predecessor of `v` on a shortest `uv`-path. Hence, this matrix efficiently
    encodes the information of a shortest `uv`-path for any `u,v\in G` : indeed,
    to go from `u` to `v` you should first find a shortest `uP[u,v]`-path, then
    jump from `P[u,v]` to `v` as it is one of its outneighbors.
    """

    cdef unsigned int n = G.order()
    cdef unsigned short * distances = <unsigned short *> sig_malloc(n*n*sizeof(unsigned short))
    if distances == NULL:
        raise MemoryError()
    cdef unsigned short * predecessors = <unsigned short *> sig_malloc(n*n*sizeof(unsigned short))
    if predecessors == NULL:
        sig_free(distances)
        raise MemoryError()
    all_pairs_shortest_path_BFS(G, predecessors, distances, NULL)

    sig_free(distances)

    return predecessors

def shortest_path_all_pairs(G):
    r"""
    Returns the matrix of predecessors in G.

    The matrix `P` returned has size `n^2`, and is such that vertex `P[u,v]` is
    a predecessor of `v` on a shortest `uv`-path. Hence, this matrix efficiently
    encodes the information of a shortest `uv`-path for any `u,v\in G` : indeed,
    to go from `u` to `v` you should first find a shortest `uP[u,v]`-path, then
    jump from `P[u,v]` to `v` as it is one of its outneighbors.

    The integer corresponding to a vertex is its index in the list
    ``G.vertices()``.

    EXAMPLE::

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

    if n == 0:
        return {}

    cdef unsigned short * predecessors = c_shortest_path_all_pairs(G)
    cdef unsigned short * c_predecessors = predecessors

    cdef dict d = {}
    cdef dict d_tmp

    cdef CGraphBackend cg = <CGraphBackend> G._backend

    cdef list int_to_vertex = G.vertices()
    cdef int i, j

    for i, l in enumerate(int_to_vertex):
        int_to_vertex[i] = cg.get_vertex(l)

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

cdef unsigned short * c_distances_all_pairs(G):
    r"""
    Returns the matrix of distances in G.

    The matrix `M` returned is of length `n^2`, and the distance between
    vertices `u` and `v` is `M[u,v]`. The integer corresponding to a vertex is
    its index in the list ``G.vertices()``.
    """

    cdef unsigned int n = G.order()
    cdef unsigned short * distances = <unsigned short *> sig_malloc(n*n*sizeof(unsigned short))
    if distances == NULL:
        raise MemoryError()
    all_pairs_shortest_path_BFS(G, NULL, distances, NULL)

    return distances

def distances_all_pairs(G):
    r"""
    Returns the matrix of distances in G.

    This function returns a double dictionary ``D`` of vertices, in which the
    distance between vertices ``u`` and ``v`` is ``D[u][v]``.

    EXAMPLE::

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

    if n == 0:
        return {}

    cdef unsigned short * distances = c_distances_all_pairs(G)
    cdef unsigned short * c_distances = distances

    cdef dict d = {}
    cdef dict d_tmp

    cdef list int_to_vertex = G.vertices()
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

def is_distance_regular(G, parameters = False):
    r"""
    Tests if the graph is distance-regular

    A graph `G` is distance-regular if for any integers `j,k` the value
    of `|\{x:d_G(x,u)=j,x\in V(G)\} \cap \{y:d_G(y,v)=j,y\in V(G)\}|` is constant
    for any two vertices `u,v\in V(G)` at distance `i` from each other.
    In particular `G` is regular, of degree `b_0` (see below), as one can take `u=v`.

    Equivalently a graph is distance-regular if there exist integers `b_i,c_i`
    such that for any two vertices `u,v` at distance `i` we have

    * `b_i = |\{x:d_G(x,u)=i+1,x\in V(G)\}\cap N_G(v)\}|, \ 0\leq i\leq d-1`
    * `c_i = |\{x:d_G(x,u)=i-1,x\in V(G)\}\cap N_G(v)\}|, \ 1\leq i\leq d,`

    where `d` is the diameter of the graph.  For more information on
    distance-regular graphs, see its associated :wikipedia:`wikipedia
    page <Distance-regular_graph>`.

    INPUT:

    - ``parameters`` (boolean) -- if set to ``True``, the function returns the
      pair ``(b,c)`` of lists of integers instead of ``True`` (see the definition
      above). Set to ``False`` by default.

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

        sage: graphs.PathGraph(2).is_distance_regular(parameters = True)
        ([1, None], [None, 1])
        sage: graphs.Tutte12Cage().is_distance_regular(parameters=True)
        ([3, 2, 2, 2, 2, 2, None], [None, 1, 1, 1, 1, 1, 3])

    """
    cdef int i,l,u,v,d,b,c,k
    cdef int n = G.order()
    cdef int infinity = <unsigned short> -1

    if n <= 1:
        return ([],[]) if parameters else True

    if not G.is_regular():
        return False
    k = G.degree(next(G.vertex_iterator()))

    # Matrix of distances
    cdef unsigned short * distance_matrix = c_distances_all_pairs(G)

    # The diameter, i.e. the longest *finite* distance between two vertices
    cdef int diameter = 0
    for i in range(n*n):
        if distance_matrix[i] > diameter and distance_matrix[i] != infinity:
            diameter = distance_matrix[i]

    cdef bitset_t b_tmp
    bitset_init(b_tmp, n)
            
    # b_distance_matrix[d*n+v] is the set of vertices at distance d from v.
    cdef binary_matrix_t b_distance_matrix
    try:
        binary_matrix_init(b_distance_matrix,n*(diameter+2),n)
    except MemoryError:
        sig_free(distance_matrix)
        bitset_free(b_tmp)
        raise

    # Fills b_distance_matrix
    for u in range(n):
        for v in range(u,n):
            d = distance_matrix[u*n+v]
            if d != infinity:
                binary_matrix_set1(b_distance_matrix, d*n+u, v)
                binary_matrix_set1(b_distance_matrix, d*n+v, u)

    cdef list bi = [-1 for i in range(diameter +1)]
    cdef list ci = [-1 for i in range(diameter +1)]

    # Applying the definition with b_i,c_i
    for u in range(n):
        for v in range(n):
            if u == v:
               continue

            d = distance_matrix[u*n+v]
            if d == infinity:
                continue

            # Computations of b_d and c_d for u,v. We intersect sets stored in
            # b_distance_matrix.
            bitset_and(b_tmp, b_distance_matrix.rows[(d+1)*n+u], b_distance_matrix.rows[n+v])
            b = bitset_len(b_tmp)
            bitset_and(b_tmp, b_distance_matrix.rows[(d-1)*n+u], b_distance_matrix.rows[n+v])
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
    else:
        return  True

###################################
# Both distances and predecessors #
###################################

def distances_and_predecessors_all_pairs(G):
    r"""
    Returns the matrix of distances in G and the matrix of predecessors.

    Distances : the matrix `M` returned is of length `n^2`, and the distance
    between vertices `u` and `v` is `M[u,v]`. The integer corresponding to a
    vertex is its index in the list ``G.vertices()``.

    Predecessors : the matrix `P` returned has size `n^2`, and is such that
    vertex `P[u,v]` is a predecessor of `v` on a shortest `uv`-path. Hence, this
    matrix efficiently encodes the information of a shortest `uv`-path for any
    `u,v\in G` : indeed, to go from `u` to `v` you should first find a shortest
    `uP[u,v]`-path, then jump from `P[u,v]` to `v` as it is one of its
    outneighbors.

    The integer corresponding to a vertex is its index in the list
    ``G.vertices()``.


    EXAMPLE::

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

    if n == 0:
        return {}, {}

    cdef unsigned short * distances = <unsigned short *> sig_malloc(n*n*sizeof(unsigned short))
    if distances == NULL:
        raise MemoryError()
    cdef unsigned short * c_distances = distances
    cdef unsigned short * predecessor = <unsigned short *> sig_malloc(n*n*sizeof(unsigned short))
    if predecessor == NULL:
        sig_free(distances)
        raise MemoryError()
    cdef unsigned short * c_predecessor = predecessor

    all_pairs_shortest_path_BFS(G, predecessor, distances, NULL)


    cdef dict d_distance = {}
    cdef dict d_predecessor = {}
    cdef dict t_distance = {}
    cdef dict t_predecessor = {}

    cdef list int_to_vertex = G.vertices()
    cdef int i, j

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

    sig_free(distances)
    sig_free(predecessor)

    return d_distance, d_predecessor

################
# Eccentricity #
################

cdef uint32_t * c_eccentricity(G) except NULL:
    r"""
    Return the vector of eccentricities in G.

    The array returned is of length n, and its ith component is the eccentricity
    of the ith vertex in ``G.vertices()``.
    """
    cdef unsigned int n = G.order()

    cdef uint32_t * ecc = <uint32_t *> sig_calloc(n, sizeof(uint32_t))
    if ecc == NULL:
        raise MemoryError()
    all_pairs_shortest_path_BFS(G, NULL, NULL, ecc)

    return ecc

cdef uint32_t * c_eccentricity_bounding(G) except NULL:
    r"""
    Return the vector of eccentricities in G using the algorithm of [TK13]_.

    The array returned is of length n, and its ith component is the eccentricity
    of the ith vertex in ``G.vertices()``.

    The algorithm proposed in [TK13]_ is based on the observation that for all
    nodes `v,w\in V`, we have `\max(ecc[v]-d(v,w), d(v,w))\leq ecc[w] \leq
    ecc[v] + d(v,w)`. Also the algorithms iteratively improves upper and lower
    bounds on the eccentricity of each node until no further improvements can be
    done. This algorithm offers good running time reduction on scale-free graphs.
    """
    if G.is_directed():
        raise ValueError("The 'bounds' algorithm only works on undirected graphs.")

    # Copying the whole graph to obtain the list of neighbors quicker than by
    # calling out_neighbors.  This data structure is well documented in the
    # module sage.graphs.base.static_sparse_graph
    cdef unsigned int n = G.order()
    cdef short_digraph sd
    init_short_digraph(sd, G)

    # allocated some data structures
    cdef bitset_t seen
    bitset_init(seen, n)
    cdef uint32_t * distances = <uint32_t *>sig_malloc(3 * n * sizeof(uint32_t))
    cdef uint32_t * LB        = <uint32_t *>sig_calloc(n, sizeof(uint32_t))
    if distances==NULL or LB==NULL:
        bitset_free(seen)
        sig_free(LB)
        sig_free(distances)
        free_short_digraph(sd)
        raise MemoryError()
    cdef uint32_t * waiting_list = distances + n
    cdef uint32_t * UB           = distances + 2 * n
    memset(UB, -1, n * sizeof(uint32_t))

    cdef uint32_t v, w, next_v, tmp, cpt = 0

    # The first vertex is the one with largest degree
    next_v = max((out_degree(sd, v), v) for v in range(n))[1]
    
    sig_on()
    while next_v!=UINT32_MAX:

        v = next_v
        cpt += 1

        # Compute the exact eccentricity of v
        LB[v] = simple_BFS(n, sd.neighbors, v, distances, NULL, waiting_list, seen)

        if LB[v]==UINT32_MAX:
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
            if LB[w]==UB[w]:
                continue
            elif next_v==UINT32_MAX or (cpt%2==0 and LB[w]<LB[next_v]) or (cpt%2==1 and UB[w]>UB[next_v]):
                # The next vertex is either the vertex with largest upper bound
                # or smallest lower bound
                next_v = w

    sig_off()

    sig_free(distances)
    bitset_free(seen)
    free_short_digraph(sd)

    return LB

@rename_keyword(deprecation=19559 , method='algorithm')
def eccentricity(G, algorithm="standard"):
    r"""
    Return the vector of eccentricities in G.

    The array returned is of length `n`, and its `i`-th component is the
    eccentricity of the ith vertex in ``G.vertices()``.

    INPUT:

    - ``G`` -- a Graph or a DiGraph.

    - ``algorithm`` -- (default: ``'standard'``) name of the method used to
      compute the eccentricity of the vertices. Available algorithms are
      ``'standard'`` which performs a BFS from each vertex and ``'bounds'``
      which uses the fast algorithm proposed in [TK13]_ for undirected graphs.

    EXAMPLE::

        sage: from sage.graphs.distances_all_pairs import eccentricity
        sage: g = graphs.PetersenGraph()
        sage: eccentricity(g)
        [2, 2, 2, 2, 2, 2, 2, 2, 2, 2]

    TEST:

    All algorithms are valid::

        sage: from sage.graphs.distances_all_pairs import eccentricity
        sage: g = graphs.RandomGNP(50,.1)
        sage: eccentricity(g, algorithm='standard')==eccentricity(g, algorithm='bounds')
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
        ValueError: The 'bounds' algorithm only works on undirected graphs.

    Asking for unknown algorithm::

        sage: from sage.graphs.distances_all_pairs import eccentricity
        sage: g = graphs.PathGraph(2)
        sage: eccentricity(g, algorithm='Nice Jazz Festival')
        Traceback (most recent call last):
        ...
        ValueError: Unknown algorithm 'Nice Jazz Festival'. Please contribute.
    """
    from sage.rings.infinity import Infinity
    cdef int n = G.order()

    # Trivial cases
    if n == 0:
        return []
    elif G.is_directed() and algorithm=='bounds':
        raise ValueError("The 'bounds' algorithm only works on undirected graphs.")
    elif not G.is_connected():
        return [Infinity]*n

    cdef uint32_t * ecc
    if algorithm=="bounds":
        ecc = c_eccentricity_bounding(G)
    elif algorithm=="standard":
        ecc = c_eccentricity(G)
    else:
        raise ValueError("Unknown algorithm '{}'. Please contribute.".format(algorithm))

    from sage.rings.integer import Integer
    cdef list l_ecc = [Integer(ecc[i]) if ecc[i]!=UINT32_MAX else +Infinity for i in range(n)]

    sig_free(ecc)

    return l_ecc


############
# Diameter #
############

cdef inline uint32_t simple_BFS(uint32_t n,
                                uint32_t ** p_vertices,
                                uint32_t source,
                                uint32_t *distances,
                                uint32_t *predecessors,
                                uint32_t *waiting_list,
                                bitset_t seen):
    """
    Perform a breadth first search (BFS) using the same method as in
    sage.graphs.distances_all_pairs.all_pairs_shortest_path_BFS

    Furthermore, the method returns the the eccentricity of the source which is
    either the last computed distance when all vertices are seen, or a very
    large number (UINT32_MAX) when the graph is not connected.

    INPUT:

    - ``n`` -- number of vertices of the graph.

    - ``p_vertices`` -- The outneighbors of vertex i are enumerated from
      p_vertices[i] to p_vertices[i+1] - 1. If p_vertices[i] is equal to
      p_vertices[i+1], then i has no outneighbours.  This data structure is well
      documented in the module sage.graphs.base.static_sparse_graph

    - ``source`` -- Starting node of the BFS.

    - ``distances`` -- array of size ``n`` to store BFS distances from
      ``source``. This method assumes that this array has already been
      allocated. However, there is no need to initialize it.

    - ``predecessors`` -- array of size ``n`` to store the first predecessor of
      each vertex during the BFS search from ``source``. The predecessor of the
      ``source`` is itself. This method assumes that this array has already
      been allocated. However, it is possible to pass a ``NULL`` pointer in
      which case the predecessors are not recorded. 

    - ``waiting_list`` -- array of size ``n`` to store the order in which the
      vertices are visited during the BFS search from ``source``. This method
      assumes that this array has already been allocated. However, there is no
      need to initialize it.

    - ``seen`` -- bitset of size ``n`` that must be initialized before calling
      this method (i.e., bitset_init(seen, n)). However, there is no need to
      clear it.

    """
    cdef uint32_t v, u
    cdef uint32_t waiting_beginning = 0
    cdef uint32_t waiting_end = 0
    cdef uint32_t * p_tmp
    cdef uint32_t * end

    # the source is seen
    bitset_clear(seen)
    bitset_add(seen, source)
    distances[source] = 0
    if predecessors!=NULL:
        predecessors[source] = source

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

        # and we iterate over all the outneighbors u of v
        while p_tmp < end:
            u = p_tmp[0]

            # If we notice one of these neighbors is not seen yet, we set its
            # parameters and add it to the queue to be explored later.
            if not bitset_in(seen, u):
                distances[u] = distances[v]+1
                bitset_add(seen, u)
                waiting_end += 1
                waiting_list[waiting_end] = u
                if predecessors!=NULL:
                    predecessors[u] = v

            p_tmp += 1

        waiting_beginning += 1

    # We return the eccentricity of the source
    return distances[waiting_list[waiting_end]] if waiting_end==n-1 else UINT32_MAX


cdef uint32_t diameter_lower_bound_2sweep(uint32_t n,
                                          uint32_t ** p_vertices,
                                          uint32_t source,
                                          uint32_t * distances,
                                          uint32_t * predecessors,
                                          uint32_t * waiting_list,
                                          bitset_t seen):
    """
    Compute a lower bound on the diameter using the 2-sweep algorithm.

    This method computes a lower bound on the diameter of an unweighted
    undirected graph using 2 BFS, as proposed in [MLH08]_.  It first selects a
    vertex `v` that is at largest distance from an initial vertex `source` using
    BFS. Then it performs a second BFS from `v`. The largest distance from `v`
    is returned as a lower bound on the diameter of `G`.  The time complexity of
    this method is linear in the size of the graph.


    INPUT:

    - ``n`` -- number of vertices of the graph.

    - ``p_vertices`` -- The outneighbors of vertex i are enumerated from
      p_vertices[i] to p_vertices[i+1] - 1. If p_vertices[i] is equal to
      p_vertices[i+1], then i has no outneighbours.  This data structure is well
      documented in the module sage.graphs.base.static_sparse_graph

    - ``source`` -- Starting node of the BFS.

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
    cdef uint32_t LB, i, k, tmp

    # We do a first BFS from source and get the eccentricity of source
    LB = simple_BFS(n, p_vertices, source, distances, NULL, waiting_list, seen)

    # If the eccentricity of the source is infinite (very large number), the
    # graph is not connected and so its diameter is infinite.
    if LB==UINT32_MAX:
        return UINT32_MAX

    # Then we perform a second BFS from the last visited vertex
    source = waiting_list[n-1]
    LB = simple_BFS(n, p_vertices, source, distances, predecessors, waiting_list, seen)

    # We return the computed lower bound
    return LB


cdef tuple diameter_lower_bound_multi_sweep(uint32_t n,
                                            uint32_t ** p_vertices,
                                            uint32_t source):
    """
    Lower bound on the diameter using multi-sweep.

    This method computes a lower bound on the diameter of an unweighted
    undirected graph using several iterations of the 2-sweep algorithms
    [CGH+13]_. Roughly, it first uses 2-sweep to identify two vertices `s` and
    `d` that are far apart. Then it selects a vertex `m` that is at same
    distance from `s` and `d`.  This vertex `m` will serve as the new source for
    another iteration of the 2-sweep algorithm that may improve the current
    lower bound on the diameter.  This process is repeated as long as the lower
    bound on the diameter is improved.

    The method returns a 4-tuple (LB, s, m, d), where LB is the best found lower
    bound on the diameter, s is a vertex of eccentricity LB, d is a vertex at
    distance LB from s, and m is a vertex at distance LB/2 from both s and d.

    INPUT:

    - ``n`` -- number of vertices of the graph.

    - ``p_vertices`` -- The outneighbors of vertex i are enumerated from
      p_vertices[i] to p_vertices[i+1] - 1. If p_vertices[i] is equal to
      p_vertices[i+1], then i has no outneighbours.  This data structure is well
      documented in the module sage.graphs.base.static_sparse_graph

    - ``source`` -- Starting node of the BFS.

    """
    # The while loop below might not be entered so we have to make sure that
    # s and d which are returned are initialized.
    cdef uint32_t LB, tmp, s = 0, m, d = 0

    # Allocate some arrays and a bitset
    cdef bitset_t seen
    bitset_init(seen, n)
    cdef uint32_t * distances = <uint32_t *>sig_malloc(3 * n * sizeof(uint32_t))
    if distances==NULL:
        bitset_free(seen)
        raise MemoryError()
    cdef uint32_t * predecessors = distances + n
    cdef uint32_t * waiting_list = distances + 2 * n


    # We perform a first 2sweep call from source. If the returned value is a
    # very large number, the graph is not connected and so the diameter is
    # infinite.
    tmp = diameter_lower_bound_2sweep(n, p_vertices, source, distances, predecessors, waiting_list, seen)
    if tmp==UINT32_MAX:
        sig_free(distances)
        bitset_free(seen)
        return (UINT32_MAX, 0, 0, 0)

    # We perform new 2sweep calls for as long as we are able to improve the
    # lower bound.
    LB = 0
    m = source
    while tmp>LB:

        LB = tmp

        # We store the vertices s, m, d of the last BFS call. For vertex m, we
        # search for a vertex of eccentricity LB/2. This vertex will serve as
        # the source for the next 2sweep call.
        s = waiting_list[0]
        d = waiting_list[n-1]
        LB_2 = LB/2
        m = d
        while distances[m]>LB_2:
            m = predecessors[m]

        # We perform a new 2sweep call from m
        tmp = diameter_lower_bound_2sweep(n, p_vertices, m, distances, predecessors, waiting_list, seen)

    sig_free(distances)
    bitset_free(seen)

    return (LB, s, m, d)


cdef uint32_t diameter_iFUB(uint32_t n,
                            uint32_t ** p_vertices,
                            uint32_t source):
    """
    Computes the diameter of the input Graph using the iFUB algorithm.

    The iFUB (iterative Fringe Upper Bound) algorithm calculates the exact value
    of the diameter of a unweighted undirected graph [CGI+10]_. This algorithms
    starts with a vertex found through a multi-sweep call (a refinement of the
    4sweep method). The worst case time complexity of the iFUB algorithm is
    `O(nm)`, but it can be very fast in practice. See the code's documentation
    and [CGH+13]_ for more details.

    INPUT:

    - ``n`` -- number of vertices of the graph.

    - ``p_vertices`` -- The outneighbors of vertex i are enumerated from
      p_vertices[i] to p_vertices[i+1] - 1. If p_vertices[i] is equal to
      p_vertices[i+1], then i has no outneighbours.  This data structure is well
      documented in the module sage.graphs.base.static_sparse_graph

    - ``source`` -- Starting node of the first BFS.

    """
    cdef uint32_t i, LB, s, m, d

    # We select a vertex m with low eccentricity using multi-sweep
    LB, s, m, d = diameter_lower_bound_multi_sweep(n, p_vertices, source)

    # If the lower bound is a very large number, it means that the graph is not
    # connected and so the diameter is infinite.
    if LB==UINT32_MAX:
        return LB


    # We allocate some arrays and a bitset
    cdef bitset_t seen
    bitset_init(seen, n)    
    cdef uint32_t * distances = <uint32_t *>sig_malloc(4 * n * sizeof(uint32_t))
    if distances==NULL:
        bitset_free(seen)
        raise MemoryError()
    cdef uint32_t * waiting_list = distances + n
    cdef uint32_t * layer        = distances + 2 * n
    cdef uint32_t * order        = distances + 3 * n


    # We order the vertices by decreasing layers. This is the inverse order of a
    # BFS from m, and so the inverse order of array waiting_list. Distances are
    # stored in array layer.
    LB = simple_BFS(n, p_vertices, m, layer, NULL, waiting_list, seen)
    for i from 0 <= i < n:
        order[i] = waiting_list[n-i-1]

    # The algorithm:
    #
    # The diameter of the graph is equal to the maximum eccentricity of a
    # vertex. Let m be any vertex, and let V be partitionned into A u B where:
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
    while (2*layer[order[i]])>LB and i<n:
        tmp = simple_BFS(n, p_vertices, order[i], distances, NULL, waiting_list, seen)
        i += 1

        # We update the lower bound
        if tmp>LB:
            LB = tmp


    sig_free(distances)
    bitset_free(seen)

    # We finally return the computed diameter
    return LB


@rename_keyword(deprecation=19559 , method='algorithm')
def diameter(G, algorithm='iFUB', source=None):
    r"""
    Returns the diameter of `G`.

    This algorithm returns Infinity if the (di)graph is not connected. It can
    also quickly return a lower bound on the diameter using the ``2sweep`` and
    ``multi-sweep`` schemes.

    INPUT:

    - ``algorithm`` -- (default: 'iFUB') specifies the algorithm to use among:

      - ``'standard'`` -- Computes the diameter of the input (di)graph as the
        largest eccentricity of its vertices. This is the classical algorithm
        with time complexity in `O(nm)`.

      - ``'2sweep'`` -- Computes a lower bound on the diameter of an
        unweighted undirected graph using 2 BFS, as proposed in [MLH08]_.  It
        first selects a vertex `v` that is at largest distance from an initial
        vertex source using BFS. Then it performs a second BFS from `v`. The
        largest distance from `v` is returned as a lower bound on the diameter
        of `G`.  The time complexity of this algorithm is linear in the size of
        `G`.

      - ``'multi-sweep'`` -- Computes a lower bound on the diameter of an
        unweighted undirected graph using several iterations of the ``2sweep``
        algorithms [CGH+13]_. Roughly, it first uses ``2sweep`` to identify
        two vertices `u` and `v` that are far apart. Then it selects a vertex
        `w` that is at same distance from `u` and `v`.  This vertex `w` will
        serve as the new source for another iteration of the ``2sweep``
        algorithm that may improve the current lower bound on the diameter.
        This process is repeated as long as the lower bound on the diameter
        is improved.

      - ``'iFUB'`` -- The iFUB (iterative Fringe Upper Bound) algorithm,
        proposed in [CGI+10]_, computes the exact value of the diameter of an
        unweighted undirected graph. It is based on the following observation:

            The diameter of the graph is equal to the maximum eccentricity of
            a vertex. Let `v` be any vertex, and let `V` be partitionned into
            `A\cup B` where:

            .. MATH::

                d(v,a)<=i, \forall a \in A\\
                d(v,b)>=i, \forall b \in B

            As all vertices from `A` are at distance `\leq 2i` from each
            other, a vertex `a\in A` with eccentricity `ecc(a)>2i` is at
            distance `ecc(a)` from some vertex `b\in B`.

            Consequently, if we have already computed the maximum eccentricity
            `m` of all vertices in `B` and if `m>2i`, then we do not need to
            compute the eccentricity of the vertices in `A`.

        Starting from a vertex `v` obtained through a multi-sweep computation
        (which refines the 4sweep algorithm used in [CGH+13]_), we compute the
        diameter by computing the eccentricity of all vertices sorted
        decreasingly according to their distance to `v`, and stop as allowed
        by the remark above. The worst case time complexity of the iFUB
        algorithm is `O(nm)`, but it can be very fast in practice.

    - ``source`` -- (default: None) vertex from which to start the first BFS.
      If ``source==None``, an arbitrary vertex of the graph is chosen. Raise an
      error if the initial vertex is not in `G`.  This parameter is not used
      when ``algorithm=='standard'``.

    EXAMPLES::

        sage: from sage.graphs.distances_all_pairs import diameter
        sage: G = graphs.PetersenGraph()
        sage: diameter(G, algorithm='iFUB')
        2
        sage: G = Graph( { 0 : [], 1 : [], 2 : [1] } )
        sage: diameter(G, algorithm='iFUB')
        +Infinity


    Although max( ) is usually defined as -Infinity, since the diameter will
    never be negative, we define it to be zero::

        sage: G = graphs.EmptyGraph()
        sage: diameter(G, algorithm='iFUB')
        0

    Comparison of exact algorithms::

        sage: G = graphs.RandomBarabasiAlbert(100, 2)
        sage: d1 = diameter(G, algorithm='standard')
        sage: d2 = diameter(G, algorithm='iFUB')
        sage: d3 = diameter(G, algorithm='iFUB', source=G.random_vertex())
        sage: if d1!=d2 or d1!=d3: print("Something goes wrong!")

    Comparison of lower bound algorithms::

        sage: lb2 = diameter(G, algorithm='2sweep')
        sage: lbm = diameter(G, algorithm='multi-sweep')
        sage: if not (lb2<=lbm and lbm<=d3): print("Something goes wrong!")

    TEST:

    This was causing a segfault. Fixed in :trac:`17873` ::

        sage: G = graphs.PathGraph(1)
        sage: diameter(G, algorithm='iFUB')
        0
    """
    cdef int n = G.order()
    if n==0:
        return 0

    if algorithm=='standard' or G.is_directed():
        return max(G.eccentricity())
    elif algorithm is None:
        algorithm = 'iFUB'
    elif not algorithm in ['2sweep', 'multi-sweep', 'iFUB']:
        raise ValueError("Unknown algorithm for computing the diameter.")

    if source is None:
        source = next(G.vertex_iterator())
    elif not G.has_vertex(source):
        raise ValueError("The specified source is not a vertex of the input Graph.")


    # Copying the whole graph to obtain the list of neighbors quicker than by
    # calling out_neighbors. This data structure is well documented in the
    # module sage.graphs.base.static_sparse_graph
    cdef short_digraph sd
    init_short_digraph(sd, G)

    # and we map the source to an int in [0,n-1] 
    cdef uint32_t isource = 0 if source is None else G.vertices().index(source)

    cdef bitset_t seen
    cdef uint32_t * tab
    cdef int LB

    if algorithm=='2sweep':
        # We need to allocate arrays and bitset
        bitset_init(seen, n)
        tab = <uint32_t *> sig_malloc(2* n * sizeof(uint32_t))
        if tab == NULL:
            free_short_digraph(sd)
            bitset_free(seen)
            raise MemoryError()
        
        LB = diameter_lower_bound_2sweep(n, sd.neighbors, isource, tab, NULL, tab+n, seen)

        bitset_free(seen)
        sig_free(tab)

    elif algorithm=='multi-sweep':
        LB = diameter_lower_bound_multi_sweep(n, sd.neighbors, isource)[0]

    else: # algorithm=='iFUB'
        LB = diameter_iFUB(n, sd.neighbors, isource)


    free_short_digraph(sd)

    if LB<0 or LB>G.order():
        from sage.rings.infinity import Infinity
        return +Infinity
    else:
        return int(LB)



################
# Wiener index #
################

def wiener_index(G):
    r"""
    Returns the Wiener index of the graph.

    The Wiener index of a graph `G` can be defined in two equivalent
    ways [KRG96b]_ :

    - `W(G) = \frac 1 2 \sum_{u,v\in G} d(u,v)` where `d(u,v)` denotes the
      distance between vertices `u` and `v`.

    - Let `\Omega` be a set of `\frac {n(n-1)} 2` paths in `G` such that `\Omega`
      contains exactly one shortest `u-v` path for each set `\{u,v\}` of
      vertices in `G`. Besides, `\forall e\in E(G)`, let `\Omega(e)` denote the
      paths from `\Omega` containing `e`. We then have
      `W(G) = \sum_{e\in E(G)}|\Omega(e)|`.

    EXAMPLE:

    From [GYLL93c]_, cited in [KRG96b]_::

        sage: g=graphs.PathGraph(10)
        sage: w=lambda x: (x*(x*x -1)/6)
        sage: g.wiener_index()==w(10)
        True
    """
    if not G.is_connected():
        from sage.rings.infinity import Infinity
        return +Infinity

    from sage.rings.integer import Integer
    cdef unsigned short * distances = c_distances_all_pairs(G)
    cdef unsigned int NN = G.order()*G.order()
    cdef unsigned int i
    cdef uint64_t s = 0
    for 0 <= i < NN:
        s += distances[i]
    sig_free(distances)
    return Integer(s)/2

##########################
# Distances distribution #
##########################

def distances_distribution(G):
    r"""
    Returns the distances distribution of the (di)graph in a dictionary.

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

        sage: D = digraphs.DeBruijn(2,3)
        sage: D.distances_distribution()
        {1: 1/4, 2: 11/28, 3: 5/14}
    """
    if G.order() <= 1:
        return {}

    from sage.rings.infinity import Infinity
    from sage.rings.integer import Integer

    cdef unsigned short * distances = c_distances_all_pairs(G)
    cdef unsigned int NN = G.order()*G.order()
    cdef dict count = {}
    cdef dict distr = {}
    cdef unsigned int i
    NNN = Integer(NN-G.order())

    # We count the number of pairs at equal distances
    for 0 <= i < NN:
        count[ distances[i] ] = count.get(distances[i],0) + 1

    sig_free(distances)

    # We normalize the distribution
    for j in count:
        if j == <unsigned short> -1:
            distr[ +Infinity ] = Integer(count[j])/NNN
        elif j > 0:
            distr[j] = Integer(count[j])/NNN

    return distr

##################
# Floyd-Warshall #
##################

def floyd_warshall(gg, paths = True, distances = False):
    r"""
    Computes the shortest path/distances between all pairs of vertices.

    For more information on the Floyd-Warshall algorithm, see the `Wikipedia
    article on Floyd-Warshall
    <http://en.wikipedia.org/wiki/Floyd%E2%80%93Warshall_algorithm>`_.

    INPUT:

    - ``gg`` -- the graph on which to work.

    - ``paths`` (boolean) -- whether to return the dictionary of shortest
      paths. Set to ``True`` by default.

    - ``distances`` (boolean) -- whether to return the dictionary of
      distances. Set to ``False`` by default.

    OUTPUT:

        Depending on the input, this function return the dictionary of paths,
        the dictionary of distances, or a pair of dictionaries
        ``(distances, paths)`` where ``distance[u][v]`` denotes the distance of a
        shortest path from `u` to `v` and ``paths[u][v]`` denotes an inneighbor
        `w` of `v` such that `dist(u,v)= 1 + dist(u,w)`.

    .. WARNING::

        Because this function works on matrices whose size is quadratic compared
        to the number of vertices, it uses short variables instead of long ones
        to divide by 2 the size in memory. This means that the current
        implementation does not run on a graph of more than 65536 nodes (this
        can be easily changed if necessary, but would require much more
        memory. It may be worth writing two versions). For information, the
        current version of the algorithm on a graph with `65536=2^{16}` nodes
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
        {(0, 1): {(0, 1): None, (1, 0): (0, 0), (0, 0): (0, 1), (1, 1): (0, 1)},
        (1, 0): {(0, 1): (0, 0), (1, 0): None, (0, 0): (1, 0), (1, 1): (1, 0)},
        (0, 0): {(0, 1): (0, 0), (1, 0): (0, 0), (0, 0): None, (1, 1): (0, 1)},
        (1, 1): {(0, 1): (1, 1), (1, 0): (1, 1), (0, 0): (0, 1), (1, 1): None}}

    Checking the distances are correct ::

        sage: g = graphs.Grid2dGraph(5,5)
        sage: dist,path = floyd_warshall(g, distances = True)
        sage: all( dist[u][v] == g.distance(u,v) for u in g for v in g )
        True

    Checking a random path is valid ::

        sage: u,v = g.random_vertex(), g.random_vertex()
        sage: p = [v]
        sage: while p[0] is not None:
        ...     p.insert(0,path[u][p[0]])
        sage: len(p) == dist[u][v] + 2
        True

    Distances for all pairs of vertices in a diamond::

        sage: g = graphs.DiamondGraph()
        sage: floyd_warshall(g, paths = False, distances = True)
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
        ValueError: The graph backend contains more than 65535 nodes
    """

    from sage.rings.infinity import Infinity
    cdef CGraph g = <CGraph> gg._backend.c_graph()[0]

    cdef list gverts = g.verts()

    if gverts == []:
        if distances and paths:
            return {}, {}
        else:
            return {}

    cdef unsigned int n = max(gverts) + 1

    if n >= <unsigned short> -1:
        raise ValueError("The graph backend contains more than "+str(<unsigned short> -1)+" nodes")

    # All this just creates two tables prec[n][n] and dist[n][n]
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef unsigned short * t_prec
    cdef unsigned short * t_dist
    cdef unsigned short ** prec
    cdef unsigned short ** dist

    cdef int i
    cdef int v_int
    cdef int u_int
    cdef int w_int

    # init dist
    t_dist = <unsigned short *>   mem.allocarray(n*n, sizeof(unsigned short))
    dist   = <unsigned short **>  mem.allocarray(n, sizeof(unsigned short *))
    dist[0] = t_dist
    for 1 <= i< n:
        dist[i] = dist[i-1] + n
    memset(t_dist, -1, n*n*sizeof(short))
    # Copying the adjacency matrix (vertices at distance 1)
    for v_int in gverts:
        dist[v_int][v_int] =  0
        for u_int in g.out_neighbors(v_int):
            dist[v_int][u_int] = 1

    if paths:
        # init prec
        t_prec = <unsigned short *>  mem.allocarray(n*n, sizeof(unsigned short))
        prec   = <unsigned short **> mem.allocarray(n, sizeof(unsigned short *))
        prec[0] = t_prec
        for 1 <= i< n:
            prec[i] = prec[i-1] + n
        memset(t_prec, 0, n*n*sizeof(short))
        # Copying the adjacency matrix (vertices at distance 1)
        for v_int in gverts:
            prec[v_int][v_int] = v_int
            for u_int in g.out_neighbors(v_int):
                prec[v_int][u_int] = v_int

    # The algorithm itself.
    cdef unsigned short *dv
    cdef unsigned short *dw
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
    cdef dict d_prec
    cdef dict d_dist
    cdef dict tmp_prec
    cdef dict tmp_dist

    cdef CGraphBackend cgb = <CGraphBackend> gg._backend

    if paths: d_prec = {}
    if distances: d_dist = {}
    for v_int in gverts:
        v = cgb.vertex_label(v_int)
        if paths: tmp_prec = {v:None}
        if distances: tmp_dist = {v:0}
        dv = dist[v_int]
        for u_int in gverts:
            u = cgb.vertex_label(u_int)
            if v != u and dv[u_int] !=  <unsigned short> -1:
                if paths:
                    tmp_prec[u] = cgb.vertex_label(prec[v_int][u_int])
                
                if distances:
                    tmp_dist[u] = dv[u_int]

        if paths: d_prec[v] = tmp_prec
        if distances: d_dist[v] = tmp_dist

    if distances and paths:
        return d_dist, d_prec
    if paths:
        return d_prec
    if distances:
        return d_dist
