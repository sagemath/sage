r"""
Static Sparse Graphs

What is the point ?
-------------------

This class implements a Cython (di)graph structure made for efficiency. The
graphs are *static*, i.e. no add/remove vertex/edges methods are available, nor
can they easily or efficiently be implemented within this data structure.

The data structure, however, is made to save the maximum amount of computations
for graph algorithms whose main operation is to *list the out-neighbours of a
vertex* (which is precisely what BFS, DFS, distance computations and the
flow-related stuff waste their life on).

The code contained in this module is written C-style. While Sage needs a class
for static graphs (not available today, i.e. 2012-01-13) it is not what we try
to address here. The purpose is efficiency and simplicity.

Author:

- Nathann Cohen (2011)

Data structure
--------------

.. image:: ../../../media/structure.png

The data structure is actually pretty simple and compact. ``short_digraph`` has
five fields

    * ``n`` (``int``) -- the number of vertices in the graph.

    * ``m`` (``int``) -- the number of edges in the graph.

    * ``edges`` (``uint32_t *``) -- array whose length is the number of edges of
      the graph.

    * ``neighbors`` (``uint32_t **``) -- this array has size `n+1`, and
      describes how the data of ``edges`` should be read : the neighbors of
      vertex `i` are the elements of ``edges`` addressed by
      ``neighbors[i]...neighbors[i+1]-1``. The element ``neighbors[n]``, which
      corresponds to no vertex (they are numbered from `0` to `n-1`) is present
      so that it remains easy to enumerate the neighbors of vertex `n-1` : the
      last of them is the element addressed by ``neighbors[n]-1``.

    * ``edge_labels`` -- this cython list associates a label to each edge of the
      graph. If a given edge is represented by ``edges[i]``, this its associated
      label can be found at ``edge_labels[i]``. This object is usually NULL,
      unless the call to ``init_short_digraph`` explicitly requires the labels
      to be stored in the data structure.

In the example given above, vertex 0 has 2,3,5,7,8 and 9 as out-neighbors, but
not 4, which is an out-neighbour of vertex 1. Vertex `n-1` has 2, 5, 8 and 9 as
out-neighbors. `\text{neighbors[n]}` points toward the cell immediately *after*
the end of `\text{edges}`, hence *outside of the allocated memory*. It is used
to indicate the end of the outneighbors of vertex `n-1`

**Iterating over the edges**

This is *the one thing* to have in mind when working with this data structure::

    cdef list_edges(short_digraph g):
        cdef int i, j
        for i in range(g.n):
            for j in range(g.neighbors[i+1]-g.neighbors[i]):
                print "There is an edge from",str(i),"to",g.neighbors[i][j]

**Advantages**

Two great points :

    * The neighbors of a vertex are C types, and are contiguous in memory.
    * Storing such graphs is incredibly cheaper than storing Python structures.

Well, I think it would be hard to have anything more efficient than that to
enumerate out-neighbors in sparse graphs ! :-)

Technical details
-----------------

    * When creating a ``fast_digraph`` from a ``Graph`` or ``DiGraph`` named
      ``G``, the `i^{\text{th}}` vertex corresponds to ``G.vertices()[i]``

    * Some methods return ``bitset_t`` objets when lists could be
      expected. There is a very useful ``bitset_list`` function for this kind of
      problems :-)

    * When the edges are labelled, most of the space taken by this graph is
      taken by edge labels. If no edge is labelled then this space is not
      allocated, but if *any* edge has a label then a (possibly empty) label is
      stored for each edge, which can double the memory needs.

    * The data structure stores the number of edges, even though it appears that
      this number can be reconstructed with
      ``g.neighbors[n]-g.neighbors[0]``. The trick is that not all elements of
      the ``g.edges`` array are necessarily used : when an undirected graph
      contains loops, only one entry of the array of size `2m` is used to store
      it, instead of the expected two. Storing the number of edges is the only
      way to avoid an uselessly costly computation to obtain the number of edges
      of an undirected, looped, AND labelled graph (think of several loops on
      the same vertex with different labels).

    * The codes of this module are well documented, and many answers can be
      found directly in the code.

Cython functions
----------------

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    ``init_short_digraph(short_digraph g, G)`` | Initializes ``short_digraph g`` from a Sage (Di)Graph.
    ``int n_edges(short_digraph g)`` | Returns the number of edges in ``g``
    ``int out_degree(short_digraph g, int i)`` | Returns the out-degree of vertex `i` in ``g``
    ``has_edge(short_digraph g, int u, int v)`` | Tests the existence of an edge.
    ``edge_label(short_digraph g, int * edge)`` | Returns the label associated with a given edge
    ``init_empty_copy(short_digraph dst, short_digraph src)`` | Allocates ``dst`` so that it can contain as many vertices and edges as ``src``.
    ``init_reverse(short_digraph dst, short_digraph src)`` | Initializes ``dst`` to a copy of ``src`` with all edges in the opposite direction.
    ``free_short_digraph(short_digraph g)`` | Free the ressources used by ``g``

**Connectivity**

``can_be_reached_from(short_digraph g, int src, bitset_t reached)``

    Assuming ``bitset_t reached`` has size at least ``g.n``, this method updates
    ``reached`` so that it represents the set of vertices that can be reached
    from ``src`` in ``g``.

``strongly_connected_component_containing_vertex(short_digraph g, short_digraph g_reversed, int v, bitset_t scc)``

    Assuming ``bitset_t reached`` has size at least ``g.n``, this method updates
    ``scc`` so that it represents the vertices of the strongly connected
    component containing ``v`` in ``g``. The variable ``g_reversed`` is assumed
    to represent the reverse of ``g``.

What is this module used for ?
------------------------------

At the moment, it is only used in the :mod:`sage.graphs.distances_all_pairs` module.

Python functions
----------------

These functions are available so that Python modules from Sage can call the
Cython routines this module implements (as they can not directly call methods
with C arguments).
"""
include "sage/misc/bitset.pxi"
cimport cpython

##############################################################################
#       Copyright (C) 2010 Nathann Cohen <nathann.cohen@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

from sage.graphs.base.c_graph cimport CGraph
from libc.stdint cimport INT32_MAX

cdef int init_short_digraph(short_digraph g, G, edge_labelled = False) except -1:
    r"""
    Initializes ``short_digraph g`` from a Sage (Di)Graph.

    If ``G`` is a ``Graph`` objet (and not a ``DiGraph``), an edge between two
    vertices `u` and `v` is replaced by two arcs in both directions.
    """
    g.edge_labels = NULL

    if G.order() >= INT32_MAX:
        raise ValueError("This structure can handle at most "+str(INT32_MAX)+" vertices !")
    else:
        g.n = G.order()

    cdef int isdigraph

    from sage.graphs.all import Graph, DiGraph

    if isinstance(G, DiGraph):
        isdigraph = 1
    elif isinstance(G, Graph):
         isdigraph = 0
    else:
        raise ValueError("The source graph must be either a DiGraph or a Graph object !")

    cdef list vertices = G.vertices()
    cdef dict v_to_id = {}
    cdef int i,j,v_id
    cdef list neighbor_label
    cdef list edge_labels

    g.m = G.size()
    cdef int n_edges = g.m if isdigraph else 2*g.m

    for i, v in enumerate(vertices):
        v_to_id[v] = i

    g.edges = <uint32_t *> sage_malloc(n_edges*sizeof(uint32_t))
    if g.edges == NULL:
        raise ValueError("Problem while allocating memory (edges)")

    g.neighbors = <uint32_t **> sage_malloc((1+<int>g.n)*sizeof(uint32_t *))
    if g.neighbors == NULL:
        raise ValueError("Problem while allocating memory (neighbors)")

    # Initializing the value of neighbors
    g.neighbors[0] = g.edges
    cdef CGraph cg = <CGraph> G._backend

    if not G.has_loops():
        # Normal case
        for i in range(1,(<int>g.n)+1):
            g.neighbors[i] = g.neighbors[i-1] + <int> (cg.out_degree(vertices[i-1]) if isdigraph else G.degree(vertices[i-1]))
    else:
        # In the presence of loops. For a funny reason, if a vertex v has a loop
        # attached to it and no other incident edge, Sage declares that it has
        # degree 2. This way, the sum of the degrees of the vertices is twice
        # the number of edges, but then the degree of a vertex is not the number
        # of its neighbors anymore. One should never try to think. It never ends
        # well.
        for i in range(1,(<int>g.n)+1):
            g.neighbors[i] = g.neighbors[i-1] + <int> len(G.edges_incident(vertices[i-1]))

    if not edge_labelled:
        for u,v in G.edge_iterator(labels = False):
            i = v_to_id[u]
            j = v_to_id[v]

            g.neighbors[i][0] = j
            g.neighbors[i] += 1

            if not isdigraph and i!=j:
                g.neighbors[j][0] = i
                g.neighbors[j] += 1

        # Reinitializing the value of neighbors
        for g.n> i >0:
            g.neighbors[i] = g.neighbors[i-1]

        g.neighbors[0] = g.edges

        # Sorting the neighbors
        for i in range(g.n):
            qsort(g.neighbors[i],g.neighbors[i+1]-g.neighbors[i],sizeof(int),compare_uint32_p)

    else:
        edge_labels = [None]*n_edges
        for v in G:
            neighbor_label = [(v_to_id[uu],l) if uu != v else (v_to_id[u],l)
                              for u,uu,l in G.edges_incident(v)]
            neighbor_label.sort()
            v_id = v_to_id[v]

            for i,(j,label) in enumerate(neighbor_label):
                g.neighbors[v_id][i] = j
                edge_labels[(g.neighbors[v_id]+i)-g.edges] = label

        g.edge_labels = <PyObject *> <void *> edge_labels
        cpython.Py_XINCREF(g.edge_labels)

cdef inline int n_edges(short_digraph g):
    # The number of edges is nothing but a difference of pointers
    return <int> (g.neighbors[g.n]-g.edges)

cdef inline int out_degree(short_digraph g, int i):
    # The out-degree is nothing but a difference of pointers
    return <int> (g.neighbors[i+1]-g.neighbors[i])

cdef int init_empty_copy(short_digraph dst, short_digraph src) except -1:
    dst.n = src.n
    dst.m = src.m
    dst.edge_labels = NULL
    cdef list edge_labels

    dst.edges = <uint32_t *> sage_malloc(n_edges(src)*sizeof(uint32_t))
    if dst.edges == NULL:
        raise ValueError("Problem while allocating memory (edges)")

    dst.neighbors = <uint32_t **> sage_malloc((src.n+1)*sizeof(uint32_t *))
    if dst.neighbors == NULL:
        raise ValueError("Problem while allocating memory (neighbors)")

    if src.edge_labels != NULL:
        edge_labels = [None]*n_edges(src)
        dst.edge_labels = <PyObject *> <void *> edge_labels
        cpython.Py_XINCREF(dst.edge_labels)

cdef int init_reverse(short_digraph dst, short_digraph src) except -1:
    cdef int i,j,v
    # Allocates memory for dst
    init_empty_copy(dst, src)

    # Avoiding a later segfault
    if dst.n == 0:
        return 0

    #### 1/3
    #
    # In a first pass, we count the in-degrees of each vertex and store it in a
    # vector. With this information, we can initialize dst.neighbors to its
    # correct value. The content of dst.edges is not touched at this level.

    cdef int * in_degree = <int *> sage_malloc(src.n*sizeof(int))
    if in_degree == NULL:
        raise ValueError("Problem while allocating memory (in_degree)")

    # Counting the degrees
    memset(in_degree, 0, src.n*sizeof(int))

    for i in range(n_edges(src)):
        in_degree[src.edges[i]] += 1

    # Updating dst.neighbors
    dst.neighbors[0] = dst.edges
    for i in range(1, src.n+1):
        dst.neighbors[i] = dst.neighbors[i-1] + in_degree[i-1]
    sage_free(in_degree)

    #### 2/3
    #
    # Second pass : we list the edges again, and add them in dst.edges. Doing
    # so, we will change the value of dst.neighbors, but that is not so bad as
    # we can fix it afterwards.
    for i in range(0, src.n):
        for j in range(out_degree(src,i)):
            v = src.neighbors[i][j]
            dst.neighbors[v][0] = i

            if dst.edge_labels != NULL:
                (<list> dst.edge_labels)[dst.neighbors[v]-dst.edges] = edge_label(src,src.neighbors[i]+j)

            dst.neighbors[v] += 1

    #### 3/3
    #
    # Final step : set the correct values of dst.neighbors again. It is easy, as
    # the correct value of dst.neighbors[i] is actually dst.neighbors[i-1]
    for src.n> i >0:
        dst.neighbors[i] = dst.neighbors[i-1]
    dst.neighbors[0] = dst.edges

    return 0

cdef int compare_uint32_p(const_void *a, const_void *b):
    return (<uint32_t *> a)[0] - (<uint32_t *> b)[0]

cdef inline uint32_t * has_edge(short_digraph g, int u, int v):
    r"""
    Tests the existence of an edge.

    Assumes that the neighbors of each vertex are sorted.
    """
    return <uint32_t *> bsearch(&v, g.neighbors[u], g.neighbors[u+1]-g.neighbors[u], sizeof(uint32_t), compare_uint32_p)

cdef inline object edge_label(short_digraph g, uint32_t * edge):
    r"""
    Returns the label associated with a given edge
    """
    if g.edge_labels == NULL:
        return None
    else:
        return (<list> g.edge_labels)[edge-g.edges]

cdef int can_be_reached_from(short_digraph g, int src, bitset_t reached) except -1:
    if g.n == 0:
        return 0

    # Initializing the set of vertices reached by setting only bit src
    bitset_set_first_n(reached, 0)
    bitset_add(reached, src)

    # We will be doing a Depth-First Search. We allocate the stack we need for
    # that, and put "src" on top of it.
    cdef int * stack = <int *> sage_malloc(g.n*sizeof(int))
    if stack == NULL:
        raise ValueError("Problem while allocating memory (stack)")

    stack[0] = src
    cdef int stack_size = 1

    # What we need to iterate over the edges
    cdef int i
    cdef uint32_t * v
    cdef uint32_t * end

    # Plain old DFS ...
    #
    #If there is something left on the stack, we remove it consider each of its
    # neighbors. If we find any which has not been reached yet, we set its
    # corresponding bit in the reached bitset, and add it on top of the stack.

    while stack_size:
        stack_size -= 1
        i = stack[stack_size]

        v = g.neighbors[i]
        end = g.neighbors[i+1]

        while v < end:
            if not bitset_in(reached, v[0]):
                bitset_add(reached, v[0])
                stack[stack_size] = v[0]
                stack_size += 1

            v += 1

    sage_free(stack)

cdef strongly_connected_component_containing_vertex(short_digraph g, short_digraph g_reversed, int v, bitset_t scc):

    # Computing the set of vertices that can be reached from v in g
    can_be_reached_from(g, v, scc)
    # Computing the set of vertices that can be reached from v in g *reversed*
    cdef bitset_t scc_reversed
    bitset_init(scc_reversed, g.n)
    can_be_reached_from(g_reversed, v, scc_reversed)
    # The scc containing v is the intersection of both sets
    bitset_intersection(scc, scc, scc_reversed)

cdef void free_short_digraph(short_digraph g):
    if g.edges != NULL:
        sage_free(g.edges)

    if g.neighbors != NULL:
        sage_free(g.neighbors)

    if g.edge_labels != NULL:
        cpython.Py_XDECREF(g.edge_labels)

def strongly_connected_components(G):
    r"""
    Returns the strongly connected components of the given DiGraph.

    INPUT:

    - ``G`` -- a DiGraph.

    .. NOTE::

        This method has been written as an attempt to solve the slowness
        reported in :trac:`12235`. It is not the one used by
        :meth:`sage.graphs.digraph.DiGraph.strongly_connected_components` as
        saving some time on the computation of the strongly connected components
        is not worth copying the whole graph, but it is a nice way to test this
        module's functions. It is also tested in the doctest or
        :meth:`sage.graphs.digraph.DiGraph.strongly_connected_components`.

    EXAMPLE::

        sage: from sage.graphs.base.static_sparse_graph import strongly_connected_components
        sage: g = digraphs.ButterflyGraph(2)
        sage: strongly_connected_components(g)
        [[('00', 0)], [('00', 1)], [('00', 2)], [('01', 0)], [('01', 1)], [('01', 2)],
        [('10', 0)], [('10', 1)], [('10', 2)], [('11', 0)], [('11', 1)], [('11', 2)]]
    """

    if G.order() == 0:
        return [[]]

    # To compute the connected component containing a given vertex v, we take
    # the intersection of the set of vertices that can be reached from v in G
    # and the set of vertices that can be reached from v in G reversed.
    #
    # That's all that happens here.

    cdef list answer = []
    cdef list vertices = G.vertices()
    cdef short_digraph g, gr

    init_short_digraph(g, G)
    init_reverse(gr, g)

    cdef bitset_t seen
    bitset_init(seen, g.n)
    bitset_set_first_n(seen, 0)

    cdef bitset_t scc
    bitset_init(scc, g.n)
    bitset_set_first_n(scc, 0)

    cdef int v
    while bitset_len(seen) < g.n:
        v = bitset_first_in_complement(seen)
        strongly_connected_component_containing_vertex(g, gr, v, scc)
        answer.append([vertices[i] for i in bitset_list(scc)])
        bitset_union(seen, seen, scc)

    bitset_free(seen)
    bitset_free(scc)
    free_short_digraph(g)
    free_short_digraph(gr)
    return answer
