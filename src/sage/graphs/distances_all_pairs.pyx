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

REFERENCE:

.. [KRG96b] S. Klavzar, A. Rajapakse, and I. Gutman. The Szeged and the
  Wiener index of graphs. *Applied Mathematics Letters*, 9(5):45--49, 1996.

.. [GYLL93c] I. Gutman, Y.-N. Yeh, S.-L. Lee, and Y.-L. Luo. Some recent
  results in the theory of the Wiener number. *Indian Journal of
  Chemistry*, 32A:651--661, 1993.

Functions
---------
"""

##############################################################################
#       Copyright (C) 2011 Nathann Cohen <nathann.cohen@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

include "sage/misc/bitset.pxi"
from libc.stdint cimport uint64_t, uint32_t, INT32_MAX
from sage.graphs.base.c_graph cimport CGraph
from sage.graphs.base.c_graph cimport vertex_label
from sage.graphs.base.c_graph cimport get_vertex

from sage.graphs.base.static_sparse_graph cimport short_digraph, init_short_digraph, free_short_digraph


cdef inline all_pairs_shortest_path_BFS(gg,
                                        unsigned short * predecessors,
                                        unsigned short * distances,
                                        int            * eccentricity):
    """
    See the module's documentation.
    """

    from sage.rings.infinity import Infinity

    cdef CGraph cg = <CGraph> gg._backend._cg

    cdef list int_to_vertex = gg.vertices()
    cdef int i

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
    cdef int * waiting_list = <int *> sage_malloc(n*sizeof(int))
    if waiting_list == NULL:
        raise MemoryError()
    cdef int waiting_beginning = 0
    cdef int waiting_end = 0

    cdef int source
    cdef int v, u
    cdef uint32_t * p_tmp
    cdef uint32_t * end

    cdef unsigned short * c_predecessors = predecessors
    cdef int * c_eccentricity = eccentricity
    cdef int * c_distances = <int *> sage_malloc( n * sizeof(int))
    if c_distances == NULL:
        sage_free(waiting_list)
        raise MemoryError()

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
        for v in range(n):
            if not bitset_in(seen, v):
                c_distances[v] = INT32_MAX
                if predecessors != NULL:
                    c_predecessors[v] = -1

        if predecessors != NULL:
            c_predecessors += n

        if eccentricity != NULL:
            c_eccentricity[source] = 0
            for i in range(n):
                c_eccentricity[source] = c_eccentricity[source] if c_eccentricity[source] >= c_distances[i] else c_distances[i]

        if distances != NULL:
            for i in range(n):
                distances[i] = <unsigned short> c_distances[i]
            distances += n

    bitset_free(seen)
    sage_free(waiting_list)
    free_short_digraph(sd)
    sage_free(c_distances)

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
    cdef unsigned short * distances = <unsigned short *> sage_malloc(n*n*sizeof(unsigned short))
    if distances == NULL:
        raise MemoryError()
    cdef unsigned short * predecessors = <unsigned short *> sage_malloc(n*n*sizeof(unsigned short))
    if predecessors == NULL:
        sage_free(distances)
        raise MemoryError()
    all_pairs_shortest_path_BFS(G, predecessors, distances, NULL)

    sage_free(distances)

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

    cdef CGraph cg = <CGraph> G._backend._cg

    cdef list int_to_vertex = G.vertices()
    cdef int i, j

    for i, l in enumerate(int_to_vertex):
        int_to_vertex[i] = get_vertex(l, G._backend.vertex_ints, G._backend.vertex_labels, cg)

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

    sage_free(predecessors)
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
    cdef unsigned short * distances = <unsigned short *> sage_malloc(n*n*sizeof(unsigned short))
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

    sage_free(distances)
    return d

def is_distance_regular(G, parameters = False):
    r"""
    Tests if the graph is distance-regular

    A graph `G` is distance-regular if there exist integers `d_1,...,d_n` such
    that for every vertex `v\in G` we have `d_i=\#\{u:d_G(u,v) =i\}`. Thus a
    strongly-regular graph is also distance-regular, and a distance-regular
    graph is necessarily regular too.

    For more information on distance-regular graphs, see its associated
    :wikipedia:`wikipedia page <Distance-regular_graph>`.

    INPUT:

    - ``parameters`` (boolean) -- whether to replace ``True`` answers with a
      dictionary associating `d_i` to an integer `i>0` if `d_i>0` (one can then
      obtain `d_i` by doing ``dictionary.get(i,0)``). Set to ``False`` by
      default.

    .. SEEALSO::

        * :meth:`Graph.is_regular`
        * :meth:`Graph.is_strongly_regular`

    EXAMPLES::

        sage: g = graphs.PetersenGraph()
        sage: g.is_distance_regular()
        True
        sage: g.is_distance_regular(parameters = True)
        {1: 3, 2: 6}

    Cube graphs, which are not strongly regular, are a bit more interesting::

        sage: graphs.CubeGraph(4).is_distance_regular(parameters = True)
        {1: 4, 2: 6, 3: 4, 4: 1}

    TESTS::

        sage: graphs.PathGraph(2).is_distance_regular(parameters = True)
        {1: 1}

    """
    cdef int i,l
    cdef int n = G.order()

    if n <= 1:
        return {} if parameters else True

    if not G.is_regular():
        return False

    cdef unsigned short * distance_matrix = c_distances_all_pairs(G)

    # - d_array is the vector of d_i corresponding to the first vertex
    #
    # - d_tmp is a vector that we use to check that d_array is the same for
    #   every vertex v
    cdef unsigned short * d_array = <unsigned short *> sage_calloc(2*n, sizeof(unsigned short))
    cdef unsigned short * d_tmp   = d_array + n

    if d_array == NULL:
        sage_free(distance_matrix)
        raise MemoryError()

    # Filling d_array
    cdef unsigned short * pointer = distance_matrix
    for i in range(n):
        if pointer[i] < n:
            d_array[pointer[i]] += 1
    pointer += n

    # For each of the n-1 other vertices
    for l in range(1,n):

        # We set d_tmp and fill it with the data from the l^th row
        memset(d_tmp, 0, n*sizeof(unsigned short))
        for i in range(n):
            if pointer[i] < n:
                d_tmp[pointer[i]] += 1

        # If d_tmp != d_array, we are done
        if memcmp(d_array, d_tmp, n*sizeof(unsigned short)) != 0:
            sage_free(distance_matrix)
            sage_free(d_array)
            return False

        pointer += n

    cdef dict dict_parameters
    if parameters:
        dict_parameters = {i:d_array[i] for i in range(n) if i and d_array[i] > 0}

    sage_free(distance_matrix)
    sage_free(d_array)

    if parameters:
        return dict_parameters
    else:
        return True

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

    cdef unsigned short * distances = <unsigned short *> sage_malloc(n*n*sizeof(unsigned short))
    if distances == NULL:
        raise MemoryError()
    cdef unsigned short * c_distances = distances
    cdef unsigned short * predecessor = <unsigned short *> sage_malloc(n*n*sizeof(unsigned short))
    if predecessor == NULL:
        sage_free(distances)
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

            if c_distances[i] == <unsigned short> -1:
                t_distance[int_to_vertex[i]] = Infinity
                t_predecessor[int_to_vertex[i]] = None
            else:
                t_distance[int_to_vertex[i]] = c_distances[i]
                t_predecessor[int_to_vertex[i]] = int_to_vertex[c_predecessor[i]]

        t_predecessor[int_to_vertex[j]] = None

        d_distance[int_to_vertex[j]] = t_distance
        d_predecessor[int_to_vertex[j]] = t_predecessor

        c_distances += n
        c_predecessor += n

    sage_free(distances)
    sage_free(predecessor)

    return d_distance, d_predecessor

################
# Eccentricity #
################

cdef int * c_eccentricity(G) except NULL:
    r"""
    Returns the vector of eccentricities in G.

    The array returned is of length n, and its ith component is the eccentricity
    of the ith vertex in ``G.vertices()``.
    """
    cdef unsigned int n = G.order()

    cdef int * ecc = <int *> sage_malloc(n*sizeof(int))
    if ecc == NULL:
        raise MemoryError()
    all_pairs_shortest_path_BFS(G, NULL, NULL, ecc)

    return ecc

def eccentricity(G):
    r"""
    Returns the vector of eccentricities in G.

    The array returned is of length n, and its ith component is the eccentricity
    of the ith vertex in ``G.vertices()``.

    EXAMPLE::

        sage: from sage.graphs.distances_all_pairs import eccentricity
        sage: g = graphs.PetersenGraph()
        sage: eccentricity(g)
        [2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
    """
    from sage.rings.infinity import Infinity
    cdef int n = G.order()

    if n == 0:
        return []

    cdef int * ecc = c_eccentricity(G)

    cdef list l_ecc = []
    cdef int i
    for i in range(n):
        if ecc[i] == INT32_MAX:
            l_ecc.append(Infinity)
        else:
            l_ecc.append(ecc[i])

    sage_free(ecc)

    return l_ecc

def diameter(G):
    r"""
    Returns the diameter of `G`.

    EXAMPLE::

        sage: from sage.graphs.distances_all_pairs import diameter
        sage: g = graphs.PetersenGraph()
        sage: diameter(g)
        2
    """
    return max(eccentricity(G))


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
    sage_free(distances)
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

    sage_free(distances)

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
        sage: print floyd_warshall(g)
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
        sage: while p[0] != None:
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
    cdef CGraph g = <CGraph> gg._backend._cg

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
    cdef unsigned short * t_prec
    cdef unsigned short * t_dist
    cdef unsigned short ** prec
    cdef unsigned short ** dist

    cdef int i
    cdef int v_int
    cdef int u_int
    cdef int w_int

    # init dist
    t_dist = <unsigned short *> sage_malloc(n*n*sizeof(unsigned short))
    if t_dist == NULL:
        raise MemoryError()
    dist = <unsigned short **> sage_malloc(n*sizeof(unsigned short *))
    if dist == NULL:
        sage_free(t_dist)
        raise MemoryError()
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
        t_prec = <unsigned short *> sage_malloc(n*n*sizeof(unsigned short))
        if t_prec == NULL:
            sage_free(t_dist)
            sage_free(dist)
            raise MemoryError()
        prec = <unsigned short **> sage_malloc(n*sizeof(unsigned short *))
        if prec == NULL:
            sage_free(t_dist)
            sage_free(dist)
            sage_free(t_prec)
            raise MemoryError()
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
    cdef unsigned short *dv, *dw
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

    cdef dict ggbvi = gg._backend.vertex_ints
    cdef dict ggbvl = gg._backend.vertex_labels

    if paths: d_prec = {}
    if distances: d_dist = {}
    for v_int in gverts:
        if paths: tmp_prec = {}
        if distances: tmp_dist = {}
        v = vertex_label(v_int, ggbvi, ggbvl, g)
        dv = dist[v_int]
        for u_int in gverts:
            u = vertex_label(u_int, ggbvi, ggbvl, g)
            if paths:
                tmp_prec[u] = (None if v == u
                               else vertex_label(prec[v_int][u_int], ggbvi, ggbvl, g))
            if distances:
                tmp_dist[u] = (dv[u_int] if (dv[u_int] != <unsigned short> -1)
                               else Infinity)
        if paths: d_prec[v] = tmp_prec
        if distances: d_dist[v] = tmp_dist

    if paths:
        sage_free(t_prec)
        sage_free(prec)

    sage_free(t_dist)
    sage_free(dist)

    if distances and paths:
        return d_dist, d_prec
    if paths:
        return d_prec
    if distances:
        return d_dist
