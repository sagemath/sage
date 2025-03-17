# cython: binding=True
# distutils: language = c++
r"""
Static sparse graphs

What is the point ?
-------------------

This class implements a Cython (di)graph structure made for efficiency. The
graphs are *static*, i.e. no add/remove vertex/edges methods are available, nor
can they easily or efficiently be implemented within this data structure.

The data structure, however, is made to save the maximum amount of computations
for graph algorithms whose main operation is to *list the out-neighbours of a
vertex* (which is precisely what BFS, DFS, distance computations and the
flow-related stuff waste their life on).

The code contained in this module is written C-style. The purpose is efficiency
and simplicity.

For an overview of graph data structures in sage, see
:mod:`~sage.graphs.base.overview`.

Author:

- Nathann Cohen (2011)

Data structure
--------------

.. image:: ../../../media/structure.png

The data structure is actually pretty simple and compact. ``short_digraph`` has
five fields

- ``n`` -- integer; the number of vertices in the graph

- ``m`` -- integer; the number of edges in the graph

- ``edges`` -- ``uint32_t *``; array whose length is the number of edges of the
  graph

- ``neighbors`` -- ``uint32_t **``; this array has size `n+1`, and describes how
  the data of ``edges`` should be read : the neighbors of vertex `i` are the
  elements of ``edges`` addressed by ``neighbors[i]...neighbors[i+1]-1``. The
  element ``neighbors[n]``, which corresponds to no vertex (they are numbered
  from `0` to `n-1`) is present so that it remains easy to enumerate the
  neighbors of vertex `n-1` : the last of them is the element addressed by
  ``neighbors[n]-1``.
  The arrays ``neighbors[i]`` are guaranteed to be sorted so the time
  complexity for deciding if ``g`` has edge `(u, v)` is `O(\log{m})` using
  binary search.

- ``edge_labels`` -- list; this cython list associates a label to each edge
  of the graph. If a given edge is represented by ``edges[i]``, this its
  associated label can be found at ``edge_labels[i]``. This object is usually
  NULL, unless the call to ``init_short_digraph`` explicitly requires the labels
  to be stored in the data structure.

In the example given above, vertex 0 has 2,3,5,7,8 and 9 as out-neighbors, but
not 4, which is an out-neighbour of vertex 1. Vertex `n-1` has 2, 5, 8 and 9 as
out-neighbors. ``neighbors[n]`` points toward the cell immediately *after* the
end of ``edges``, hence *outside of the allocated memory*. It is used to
indicate the end of the outneighbors of vertex `n-1`

**Iterating over the edges**

This is *the one thing* to have in mind when working with this data structure::

    cdef list_edges(short_digraph g):
        cdef int i, j
        for i in range(g.n):
            for j in range(g.neighbors[i+1]-g.neighbors[i]):
                print("There is an edge from {} to {}".format(i, g.neighbors[i][j]))

**Advantages**

Two great points :

- The neighbors of a vertex are C types, and are contiguous in memory.
- Storing such graphs is incredibly cheaper than storing Python structures.

Well, I think it would be hard to have anything more efficient than that to
enumerate out-neighbors in sparse graphs ! :-)

Technical details
-----------------

* When creating a ``short_digraph`` from a ``Graph`` or ``DiGraph`` named ``G``,
  the `i^{\text{th}}` vertex corresponds *by default* to ``list(G)[i]``.
  Using optional parameter ``vertex_list``, you can specify the order of the
  vertices. Then `i^{\text{th}}` vertex will corresponds to ``vertex_list[i]``.

* Some methods return ``bitset_t`` objects when lists could be expected. There
  is a very useful ``bitset_list`` function for this kind of problems :-)

* When the edges are labelled, most of the space taken by this graph is taken by
  edge labels. If no edge is labelled then this space is not allocated, but if
  *any* edge has a label then a (possibly empty) label is stored for each edge,
  which can double the memory needs.

* The data structure stores the number of edges, even though it appears that
  this number can be reconstructed with ``g.neighbors[n]-g.neighbors[0]``. The
  trick is that not all elements of the ``g.edges`` array are necessarily used :
  when an undirected graph contains loops, only one entry of the array of size
  `2m` is used to store it, instead of the expected two. Storing the number of
  edges is the only way to avoid an uselessly costly computation to obtain the
  number of edges of an undirected, looped, AND labelled graph (think of several
  loops on the same vertex with different labels).

* The codes of this module are well documented, and many answers can be found
  directly in the code.

Cython functions
----------------

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    ``init_short_digraph(short_digraph g, G, edge_labelled, vertex_list)`` | Initialize ``short_digraph g`` from a Sage (Di)Graph.
    ``int n_edges(short_digraph g)`` | Return the number of edges in ``g``
    ``int out_degree(short_digraph g, int i)`` | Return the out-degree of vertex `i` in ``g``
    ``has_edge(short_digraph g, int u, int v)`` | Test the existence of an edge.
    ``edge_label(short_digraph g, int * edge)`` | Return the label associated with a given edge
    ``init_empty_copy(short_digraph dst, short_digraph src)`` | Allocate ``dst`` so that it can contain as many vertices and edges as ``src``.
    ``init_reverse(short_digraph dst, short_digraph src)`` | Initialize ``dst`` to a copy of ``src`` with all edges in the opposite direction.
    ``free_short_digraph(short_digraph g)`` | Free the resources used by ``g``

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

``tarjan_strongly_connected_components_C(short_digraph g, int *scc)``

    Assuming ``scc`` is already allocated and has size at least ``g.n``, this
    method computes the strongly connected components of ``g``, and outputs in
    ``scc[v]`` the number of the strongly connected component containing ``v``.
    It returns the number of strongly connected components.

``strongly_connected_components_digraph_C(short_digraph g, int nscc, int *scc, short_digraph output):``

    Assuming ``nscc`` and ``scc`` are the outputs of
    ``tarjan_strongly_connected_components_C`` on ``g``, this routine
    sets ``output`` to the
    strongly connected component digraph of ``g``, that is, the vertices of
    ``output`` are the strongly connected components of ``g`` (numbers are
    provided by ``scc``), and ``output`` contains an arc ``(C1,C2)`` if ``g``
    has an arc from a vertex in ``C1`` to a vertex in ``C2``.

What is this module used for ?
------------------------------

It is for instance used in the :mod:`sage.graphs.distances_all_pairs` module,
and in the :meth:`~sage.graphs.digraph.DiGraph.strongly_connected_components`
method.

Python functions
----------------

These functions are available so that Python modules from Sage can call the
Cython routines this module implements (as they cannot directly call methods
with C arguments).
"""
# ****************************************************************************
#       Copyright (C) 2010 Nathann Cohen <nathann.cohen@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

cimport cpython
from libc.limits cimport INT_MAX
from libc.math cimport sqrt
from libcpp.vector cimport vector
from cysignals.memory cimport check_allocarray, check_calloc, sig_free
from cysignals.signals cimport sig_on, sig_off
from cython.operator cimport postincrement
from memory_allocator cimport MemoryAllocator

from sage.data_structures.bitset_base cimport *
from sage.graphs.base.static_sparse_backend cimport StaticSparseCGraph
from sage.graphs.base.static_sparse_backend cimport StaticSparseBackend


cdef extern from "fenv.h":
    int FE_TONEAREST
    int FE_UPWARD
    int FE_DOWNWARD
    int FE_TOWARDZERO
    int fegetround ()
    int fesetround (int)


cdef int init_short_digraph(short_digraph g, G, edge_labelled=False,
                            vertex_list=None) except -1:
    r"""
    Initialize ``short_digraph g`` from a Sage (Di)Graph.

    INPUT:

    - ``g`` -- a short_digraph

    - ``G`` -- a ``Graph`` or a ``DiGraph``. If ``G`` is a ``Graph`` object,
      then any edge between two vertices `u` and `v` is replaced by two arcs in
      both directions.

    - ``edge_labelled`` -- boolean (default: ``False``); whether to store the
      label of edges or not

    - ``vertex_list`` -- list (default: ``None``); list of all vertices of ``G``
      in some order. When given, it is used to map the vertices of the graph to
      consecutive integers. Otherwise, the result of ``list(G)`` is used
      instead. Beware that if ``vertex_list`` is not ``None``, it is not checked
      and this function assumes that it contains a permutation of the vertices
      of the graph ``G``.

    COMPLEXITY:

    The time complexity for initializing ``g`` is `O(n + m)` for ``SparseGraph``
    and `O(n^2)` for ``DenseGraph``.

    TESTS:

    Indirect doctests for sorted output::

        sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
        sage: G = graphs.CompleteGraph(5).relabel(list('abcde'), inplace=False)
        sage: B = StaticSparseBackend(G, sort=False)
        sage: list(B.iterator_nbrs('a'))
        ['b', 'c', 'd', 'e']
        sage: G = graphs.CompleteGraph(5).relabel(list('badec'), inplace=False)
        sage: B = StaticSparseBackend(G, sort=False)
        sage: list(B.iterator_nbrs('a'))
        ['b', 'd', 'e', 'c']

    Same with labels::

        sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
        sage: G = graphs.CompleteGraph(5).relabel(list('abcde'), inplace=False)
        sage: for i, (u, v, _) in enumerate(G.edges()):
        ....:     G.set_edge_label(u, v, f'{u}{v}')
        sage: B = StaticSparseBackend(G, sort=False)
        sage: list(B.iterator_edges('a', True))
        [('a', 'b', 'ab'), ('a', 'c', 'ac'), ('a', 'd', 'ad'), ('a', 'e', 'ae')]
        sage: G = graphs.CompleteGraph(5).relabel(list('badec'), inplace=False)
        sage: for i, (u, v, _) in enumerate(G.edges()):
        ....:     G.set_edge_label(u, v, f'{u}{v}')
        sage: B = StaticSparseBackend(G, sort=False)
        sage: list(B.iterator_edges('a', True))
        [('a', 'b', 'ab'), ('a', 'd', 'ad'), ('a', 'e', 'ae'), ('a', 'c', 'ac')]
    """
    from sage.graphs.graph import Graph
    from sage.graphs.digraph import DiGraph

    if not isinstance(G, (Graph, DiGraph)):
        raise ValueError("The source graph must be either a DiGraph or a Graph"
                         "object !")

    if G.order() >= INT_MAX:
        raise ValueError(f"short_digraph can handle at most {INT_MAX} vertices")

    g.edge_labels = NULL
    g.n = G.order()
    g.m = G.size()

    cdef int isdigraph = G.is_directed()
    cdef uint32_t i, v_id, j
    cdef list vertices = vertex_list if vertex_list is not None else list(G)
    cdef dict v_to_id = {v: i for i, v in enumerate(vertices)}
    cdef list neighbor_label
    cdef list edge_labels
    # Loops are not stored twice for undirected graphs
    cdef int n_edges = g.m if isdigraph else 2*g.m - G.number_of_loops()

    g.edges = <uint32_t *>check_allocarray(n_edges, sizeof(uint32_t))
    g.neighbors = <uint32_t **>check_allocarray(1 + g.n, sizeof(uint32_t *))

    # Initializing the value of neighbors
    g.neighbors[0] = g.edges

    if not G.has_loops():
        # Normal case
        for i, v in enumerate(vertices):
            g.neighbors[i+1] = g.neighbors[i] + <int>(G.out_degree(v) if isdigraph else G.degree(v))
    else:
        # In the presence of loops. For a funny reason, if a vertex v has a loop
        # attached to it and no other incident edge, Sage declares that it has
        # degree 2. This way, the sum of the degrees of the vertices is twice
        # the number of edges, but then the degree of a vertex is not the number
        # of its neighbors anymore. One should never try to think. It never ends
        # well.
        for i, v in enumerate(vertices):
            g.neighbors[i+1] = g.neighbors[i] + <int> len(G.edges_incident(v))

    if edge_labelled:
        edge_labels = [None] * n_edges

    # Note that neighbors[i] will be naturally sorted by increasing id,
    # because the arrays will be built by appending vertices in the same
    # order as they appear in ``vertices``
    for i, v in enumerate(vertices):
        if isdigraph:
            edge_iterator = G.incoming_edge_iterator(v,
                                                     labels=edge_labelled)
        else:
            edge_iterator = G.edge_iterator(v, labels=edge_labelled,
                                            sort_vertices=False)
        for e in edge_iterator:
            u = e[0] if v == e[1] else e[1]
            j = v_to_id[u]
            # Handle the edge u -> v of G (= the edge j -> i of g)
            g.neighbors[j][0] = i
            # Note: cannot use the dereference Cython operator here, do not
            # known why but the following line does not compile
            # dereference(g.neighbors[j]) = i
            if edge_labelled:
                edge_labels[g.neighbors[j] - g.edges] = e[2]
            postincrement(g.neighbors[j])  # increment pointer to next item

    # Reinitializing the value of neighbors
    for i in range(g.n-1, 0, -1):
        g.neighbors[i] = g.neighbors[i-1]
    g.neighbors[0] = g.edges

    if edge_labelled:
        g.edge_labels = <PyObject *> <void *> edge_labels
        cpython.Py_XINCREF(g.edge_labels)


cdef inline int n_edges(short_digraph g) noexcept:
    """
    Return the number of edges in ``g``.

    The number of edges is nothing but a difference of pointers
    """
    return <int> (g.neighbors[g.n] - g.edges)


cdef inline int out_degree(short_digraph g, int i) noexcept:
    """
    Return the out-degree of vertex `i` in ``g``.

    The out-degree is nothing but a difference of pointers
    """
    return <int> (g.neighbors[i + 1] - g.neighbors[i])


cdef int init_empty_copy(short_digraph dst, short_digraph src) except -1:
    """
    Allocate ``dst`` so that it can contain as many vertices and edges as
    ``src``.
    """
    dst.n = src.n
    dst.m = src.m
    dst.edge_labels = NULL
    cdef list edge_labels

    dst.edges = <uint32_t *>check_allocarray(n_edges(src), sizeof(uint32_t))
    dst.neighbors = <uint32_t **>check_allocarray(src.n + 1, sizeof(uint32_t *))

    if src.edge_labels != NULL:
        edge_labels = [None] * n_edges(src)
        dst.edge_labels = <PyObject *> <void *> edge_labels
        cpython.Py_XINCREF(dst.edge_labels)


cdef int init_reverse(short_digraph dst, short_digraph src) except -1:
    """
    Initialize ``dst`` to a copy of ``src`` with all edges in the opposite
    direction.

    TESTS:

    Indirect doctests for sorted output::

        sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
        sage: D = digraphs.Complete(5).relabel(list('abcde'), inplace=False)
        sage: for i, (u, v, _) in enumerate(D.edges()):
        ....:     D.set_edge_label(u, v, f'{u}{v}')
        sage: B = StaticSparseBackend(D, sort=False)
        sage: list(B.iterator_in_edges('a', True))
        [('b', 'a', 'ba'), ('c', 'a', 'ca'), ('d', 'a', 'da'), ('e', 'a', 'ea')]
        sage: D = digraphs.Complete(5).relabel(list('badec'), inplace=False)
        sage: for i, (u, v, _) in enumerate(D.edges()):
        ....:     D.set_edge_label(u, v, f'{u}{v}')
        sage: B = StaticSparseBackend(D, sort=False)
        sage: list(B.iterator_in_edges('a', True))
        [('b', 'a', 'ba'), ('d', 'a', 'da'), ('e', 'a', 'ea'), ('c', 'a', 'ca')]
    """
    cdef int i, j, v
    # Allocates memory for dst
    init_empty_copy(dst, src)

    # Avoiding a later segfault
    if not dst.n:
        return 0

    # 1/3
    #
    # In a first pass, we count the in-degrees of each vertex and store it in a
    # vector. With this information, we can initialize dst.neighbors to its
    # correct value. The content of dst.edges is not touched at this level.
    cdef int * in_degree = <int *>check_calloc(src.n, sizeof(int))

    for i in range(n_edges(src)):
        in_degree[src.edges[i]] += 1

    # Updating dst.neighbors
    dst.neighbors[0] = dst.edges
    for i in range(1, src.n + 1):
        dst.neighbors[i] = dst.neighbors[i - 1] + in_degree[i - 1]
    sig_free(in_degree)

    # 2/3
    #
    # Second pass : we list the edges again, and add them in dst.edges. Doing
    # so, we will change the value of dst.neighbors, but that is not so bad as
    # we can fix it afterwards.
    # Note that neighbors[i] will be naturally sorted by increasing id,
    # because the arrays will be built by appending vertices in the same
    # order as they appear in ``vertices``
    for i in range(0, src.n):
        for j in range(out_degree(src, i)):
            v = src.neighbors[i][j]
            dst.neighbors[v][0] = i

            if dst.edge_labels:
                (<list> dst.edge_labels)[dst.neighbors[v] - dst.edges] = edge_label(src, src.neighbors[i] + j)

            dst.neighbors[v] += 1

    # 3/3
    #
    # Third step : set the correct values of dst.neighbors again. It is easy, as
    # the correct value of dst.neighbors[i] is actually dst.neighbors[i-1]
    for i in range(src.n - 1, 0, -1):
        dst.neighbors[i] = dst.neighbors[i - 1]
    dst.neighbors[0] = dst.edges

    return 0


cdef int compare_uint32_p(const_void *a, const_void *b) noexcept:
    """
    Comparison function needed for ``bsearch``.
    """
    return (<uint32_t *> a)[0] - (<uint32_t *> b)[0]


cdef inline uint32_t * has_edge(short_digraph g, int u, int v) noexcept:
    r"""
    Test the existence of an edge.

    Return a pointer to ``v`` in the list of neighbors of ``u`` if found and
    ``NULL`` otherwise.

    .. NOTE::

    Use the fact that the array ``g.neighbors[u]`` is guaranteed to be sorted.
    """
    # The neighbors of u are sorted by increasing label. We can use binary
    # search to decide if g has edge (u, v)
    return <uint32_t *> bsearch(&v, g.neighbors[u],
                                g.neighbors[u+1] - g.neighbors[u],
                                sizeof(uint32_t), compare_uint32_p)


cdef inline object edge_label(short_digraph g, uint32_t * edge):
    r"""
    Return the label associated with a given edge
    """
    if not g.edge_labels:
        return None
    return (<list> g.edge_labels)[edge - g.edges]


cdef uint32_t simple_BFS(short_digraph g,
                         uint32_t source,
                         uint32_t *distances,
                         uint32_t *predecessors,
                         uint32_t *waiting_list,
                         bitset_t seen) noexcept:
    """
    Perform a breadth first search (BFS) using the same method as in
    sage.graphs.distances_all_pairs.all_pairs_shortest_path_BFS

    Furthermore, the method returns the eccentricity of the source which is
    either the last computed distance when all vertices are seen, or a very
    large number (UINT32_MAX) when the graph is not connected.

    INPUT:

    - ``g`` -- a short_digraph

    - ``source`` -- starting node of the BFS

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
    cdef uint32_t n = g.n
    cdef uint32_t ** p_vertices = g.neighbors

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
        end = p_vertices[v + 1]

        # and we iterate over all the outneighbors u of v
        while p_tmp < end:
            u = p_tmp[0]

            # If we notice one of these neighbors is not seen yet, we set its
            # parameters and add it to the queue to be explored later.
            if not bitset_in(seen, u):
                distances[u] = distances[v] + 1
                bitset_add(seen, u)
                waiting_end += 1
                waiting_list[waiting_end] = u
                if predecessors:
                    predecessors[u] = v

            p_tmp += 1

        waiting_beginning += 1

    # We return the eccentricity of the source
    return distances[waiting_list[waiting_end]] if waiting_end == n - 1 else UINT32_MAX


cdef int can_be_reached_from(short_digraph g, int src, bitset_t reached) except -1:
    """
    Feed ``reached`` with the vertices reachable from ``src``.
    """
    if not g.n:
        return 0

    # Initializing the set of vertices reached by setting only bit src
    bitset_set_first_n(reached, 0)
    bitset_add(reached, src)

    # We will be doing a Depth-First Search. We allocate the stack we need for
    # that, and put "src" on top of it.
    cdef int * stack = <int *>check_allocarray(g.n, sizeof(int))

    stack[0] = src
    cdef int stack_size = 1

    # What we need to iterate over the edges
    cdef int i
    cdef uint32_t * v
    cdef uint32_t * end

    # Plain old DFS ...
    #
    # If there is something left on the stack, we remove it consider each of its
    # neighbors. If we find any which has not been reached yet, we set its
    # corresponding bit in the reached bitset, and add it on top of the stack.

    while stack_size:
        stack_size -= 1
        i = stack[stack_size]

        v = g.neighbors[i]
        end = g.neighbors[i + 1]

        while v < end:
            if not bitset_in(reached, v[0]):
                bitset_add(reached, v[0])
                stack[stack_size] = v[0]
                stack_size += 1

            v += 1

    sig_free(stack)


cdef int tarjan_strongly_connected_components_C(short_digraph g, int *scc) noexcept:
    r"""
    The Tarjan algorithm to compute strongly connected components (SCCs).

    This routine returns the number of SCCs `k` and, stores in ``scc[v]`` an
    integer between `0` and `k-1`, corresponding to the SCC containing v. SCCs
    are numbered in reverse topological order, that is, if `(v,w)` is an edge
    in the graph, ``scc[v] <= scc[w]``.

    The basic idea of the algorithm is this: a depth-first search (DFS) begins
    from an arbitrary start node (and subsequent DFSes are
    conducted on any nodes that have not yet been found). As usual with DFSes,
    the search visits every node of the graph exactly once, declining to revisit
    any node that has already been explored. Thus, the collection of search
    trees is a spanning forest of the graph. The strongly connected components
    are the subtrees of this spanning forest having no edge directed outside the
    subtree.

    To recover these components, during the DFS, we keep the index of a node,
    that is, the position in the DFS tree, and the lowlink: as soon as the
    subtree rooted at `v` has been fully explored, the lowlink of `v` is the
    smallest index reachable from `v` passing from descendants of `v`. If the
    subtree rooted at `v` has been fully explored, and the index of `v` equals
    the lowlink of `v`, that whole subtree is a new SCC.
    """
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef int u, v, w, n = g.n, current_index = 0, currentscc = 0
    cdef int *index = <int *> mem.malloc(n * sizeof(int))
    cdef int *pred = <int *> mem.malloc(n * sizeof(int))
    cdef int *lowlink = <int *> mem.malloc(n * sizeof(int))
    cdef int *dfs_stack = <int *> mem.malloc((n_edges(g) + 1) * sizeof(int))
    cdef int dfs_stack_end
    cdef int *scc_stack = <int *> mem.malloc(n * sizeof(int))  # Used to keep track of which nodes are in the "current" SCC
    cdef short *in_scc_stack = <short *> mem.calloc(n, sizeof(short))
    cdef uint32_t *p_tmp
    cdef short *visited = <short *> mem.calloc(n, sizeof(short))
    # The variable visited[v] is 0 if the vertex has never been visited, 1 if
    # it is an ancestor of the current vertex, 2 otherwise.

    for u in range(n):
        if visited[u]:
            continue

        # Perform a DFS from u
        dfs_stack_end = 1
        scc_stack_end = 0
        dfs_stack[0] = u
        pred[u] = u

        while dfs_stack_end:
            v = dfs_stack[dfs_stack_end - 1]
            if not visited[v]:
                # It means that this is the first time we visit v.
                # We set the index and the lowlink to be equal: during the
                # algorithm, the lowlink may decrease.
                visited[v] = 1
                index[v] = current_index
                lowlink[v] = current_index
                current_index = current_index + 1
                # We add v to the stack of vertices in the current SCC
                scc_stack[scc_stack_end] = v
                scc_stack_end = scc_stack_end + 1
                in_scc_stack[v] = 1

                # We iterate over all neighbors of v
                p_tmp = g.neighbors[v]
                while p_tmp < g.neighbors[v + 1]:
                    w = p_tmp[0]
                    p_tmp += 1
                    if not visited[w]:
                        # Vertex w is added to the DFS stack
                        pred[w] = v
                        dfs_stack[dfs_stack_end] = w
                        dfs_stack_end += 1
                    elif in_scc_stack[w]:
                        # We update the lowlink of v (later, we will "pass"
                        # this updated value to all ancestors of v.
                        lowlink[v] = min(lowlink[v], lowlink[w])
            else:
                # The vertex v has already been visited.
                dfs_stack_end -= 1

                if visited[v] == 1:
                    # It means that we have just processed all the DFS
                    # subtree rooted at v. Hence, the lowlink of v is the
                    # final value, and we "pass" this value to the
                    # predecessor of v.
                    lowlink[pred[v]] = min(lowlink[pred[v]], lowlink[v])

                    if lowlink[v] == index[v]:
                        # The DFS subtree rooted at v is a new SCC. We
                        # recover the SCC from scc_stack.
                        w = -1
                        while w != v:
                            scc_stack_end -= 1
                            w = scc_stack[scc_stack_end]
                            in_scc_stack[w] = 0
                            scc[w] = currentscc
                        currentscc += 1
                    visited[v] = 2

    return currentscc


def tarjan_strongly_connected_components(G):
    r"""
    Return the lists of vertices in each strongly connected components (SCCs).

    This method implements the Tarjan algorithm to compute the strongly
    connected components of the digraph. It returns a list of lists of vertices,
    each list of vertices representing a strongly connected component.

    The basic idea of the algorithm is this: a depth-first search (DFS) begins
    from an arbitrary start node (and subsequent DFSes are
    conducted on any nodes that have not yet been found). As usual with DFSes,
    the search visits every node of the graph exactly once, declining to revisit
    any node that has already been explored. Thus, the collection of search
    trees is a spanning forest of the graph. The strongly connected components
    correspond to the subtrees of this spanning forest that have no edge
    directed outside the subtree.

    To recover these components, during the DFS, we keep the index of a node,
    that is, the position in the DFS tree, and the lowlink: as soon as the
    subtree rooted at `v` has been fully explored, the lowlink of `v` is the
    smallest index reachable from `v` passing from descendants of `v`. If the
    subtree rooted at `v` has been fully explored, and the index of `v` equals
    the lowlink of `v`, that whole subtree is a new SCC.

    For more information, see the
    :wikipedia:`Tarjan%27s_strongly_connected_components_algorithm`.

    EXAMPLES::

        sage: from sage.graphs.base.static_sparse_graph import tarjan_strongly_connected_components
        sage: tarjan_strongly_connected_components(digraphs.Path(3))
        [[2], [1], [0]]
        sage: D = DiGraph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
        sage: D.connected_components(sort=True)
        [[0, 1, 2, 3], [4, 5, 6]]
        sage: D = DiGraph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
        sage: D.strongly_connected_components()
        [[3], [2], [1], [0], [6], [5], [4]]
        sage: D.add_edge([2,0])
        sage: D.strongly_connected_components()
        [[3], [0, 1, 2], [6], [5], [4]]
        sage: D = DiGraph([('a','b'), ('b','c'), ('c', 'd'), ('d', 'b'), ('c', 'e')])
        sage: [sorted(scc) for scc in D.strongly_connected_components()]
        [['e'], ['b', 'c', 'd'], ['a']]

    TESTS:

    Checking that the result is correct::

        sage: from sage.graphs.base.static_sparse_graph import tarjan_strongly_connected_components
        sage: import random
        sage: for i in range(10):                          # long time
        ....:     n = random.randint(2,20)
        ....:     m = random.randint(1, n*(n-1))
        ....:     g = digraphs.RandomDirectedGNM(n,m)
        ....:     sccs = tarjan_strongly_connected_components(g)
        ....:     for scc in sccs:
        ....:         scc_check = g.strongly_connected_component_containing_vertex(scc[0])
        ....:         assert(sorted(scc) == sorted(scc_check))

    Checking against NetworkX::

        sage: import networkx                                                                       # needs networkx
        sage: for i in range(10):               # long time                             # needs networkx
        ....:      g = digraphs.RandomDirectedGNP(100,.05)
        ....:      h = g.networkx_graph()
        ....:      scc1 = g.strongly_connected_components()
        ....:      scc2 = networkx.strongly_connected_components(h)
        ....:      s1 = Set(map(Set,scc1))
        ....:      s2 = Set(map(Set,scc2))
        ....:      if s1 != s2:
        ....:          print("Ooch !")

    Immutable digraphs::

        sage: from sage.graphs.base.static_sparse_graph import tarjan_strongly_connected_components
        sage: G = digraphs.RandomDirectedGNP(10, .4)
        sage: G._backend
        <sage.graphs.base.sparse_graph.SparseGraphBackend ...>
        sage: H = DiGraph(G, immutable=True)
        sage: H._backend
        <sage.graphs.base.static_sparse_backend.StaticSparseBackend ...>
        sage: tarjan_strongly_connected_components(G) == tarjan_strongly_connected_components(H)
        True
    """
    from sage.graphs.digraph import DiGraph

    if not isinstance(G, DiGraph):
        raise ValueError("G must be a DiGraph.")

    cdef list int_to_vertex
    cdef StaticSparseCGraph cg
    cdef short_digraph g
    if isinstance(G, StaticSparseBackend):
        cg = <StaticSparseCGraph> G._cg
        g = <short_digraph> cg.g
        int_to_vertex = cg._vertex_to_labels
    else:
        int_to_vertex = list(G)
        init_short_digraph(g, G, edge_labelled=False, vertex_list=int_to_vertex)

    cdef MemoryAllocator mem = MemoryAllocator()
    cdef int * scc = <int*> mem.malloc(g.n * sizeof(int))
    sig_on()
    cdef int nscc = tarjan_strongly_connected_components_C(g, scc)
    sig_off()

    if not isinstance(G, StaticSparseBackend):
        free_short_digraph(g)

    cdef int i
    cdef list output = [[] for i in range(nscc)]

    for i, v in enumerate(int_to_vertex):
        output[scc[i]].append(v)
    return output


cdef void strongly_connected_components_digraph_C(short_digraph g, int nscc, int *scc, short_digraph output) noexcept:
    r"""
    Compute the strongly connected components (SCCs) digraph of `g`.

    The strongly connected components digraph of `g` is a graph having a vertex
    for each SCC of `g` and an arc from component `C_1` to component `C_2` if
    and only if there is an arc in `g` from a vertex in `C_1` to a vertex in
    `C_2`. The strongly connected components digraph is acyclic by definition.

    This routine inputs the graph ``g``, the number of SCCs ``nscc``, and an
    array containing in position ``v`` the SCC of vertex ``v`` (these values
    must be already computed). The output is stored in variable ``output``,
    which should be empty at the beginning.
    """
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef size_t v, w, i
    cdef size_t s_nscc = <size_t>nscc
    cdef vector[vector[int]] scc_list = vector[vector[int]](nscc, vector[int]())
    cdef vector[vector[int]] sons = vector[vector[int]](nscc + 1, vector[int]())
    cdef short *neighbors = <short *> mem.calloc(nscc, sizeof(short))
    cdef long m = 0
    cdef uint32_t *p_tmp

    for v in range(s_nscc):
        scc_list[v] = vector[int]()
        sons[v] = vector[int]()
    sons[nscc] = vector[int]()

    for i in range(<size_t>g.n):
        scc_list[scc[i]].push_back(i)

    for v in range(s_nscc):
        for i in range(scc_list[v].size()):
            p_tmp = g.neighbors[scc_list[v][i]]
            while p_tmp<g.neighbors[scc_list[v][i] + 1]:
                w = <int> scc[p_tmp[0]]
                p_tmp += 1
                if not (neighbors[w] or w == v):
                    neighbors[w] = 1
                    sons[v].push_back(w)
                    m += 1
        for w in range(sons[v].size()):
            neighbors[sons[v][w]] = 0

    output.n = nscc
    output.m = m
    output.edge_labels = NULL

    output.neighbors = <uint32_t **> check_allocarray((1+<int>output.n), sizeof(uint32_t *))

    if not m:
        output.edges = NULL
        for v in range(1, s_nscc + 1):
            output.neighbors[v] = NULL

    output.edges = <uint32_t *> check_allocarray(m, sizeof(uint32_t))
    output.neighbors[0] = output.edges

    for v in range(1, s_nscc + 1):
        output.neighbors[v] = output.neighbors[v - 1] + sons[v - 1].size()
        for i in range(sons[v].size()):
            output.neighbors[v][i] = sons[v][i]


def strongly_connected_components_digraph(G):
    r"""
    Return the digraph of the strongly connected components (SCCs).

    This routine is used to test ``strongly_connected_components_digraph_C``,
    but it is not used by the Sage digraph. It outputs a pair ``[g_scc,scc]``,
    where ``g_scc`` is the SCC digraph of g, ``scc`` is a dictionary associating
    to each vertex ``v`` the number of the SCC of ``v``, as it appears in
    ``g_scc``.

    EXAMPLES::

        sage: from sage.graphs.base.static_sparse_graph import strongly_connected_components_digraph
        sage: strongly_connected_components_digraph(digraphs.Path(3))
        (Digraph on 3 vertices, {0: 2, 1: 1, 2: 0})
        sage: strongly_connected_components_digraph(DiGraph(4))
        (Digraph on 4 vertices, {0: 0, 1: 1, 2: 2, 3: 3})

    TESTS::

        sage: from sage.graphs.base.static_sparse_graph import strongly_connected_components_digraph
        sage: import random
        sage: for i in range(100):
        ....:     n = random.randint(2,20)
        ....:     m = random.randint(1, n*(n-1))
        ....:     g = digraphs.RandomDirectedGNM(n,m)
        ....:     scc_digraph,sccs = strongly_connected_components_digraph(g)
        ....:     assert(scc_digraph.is_directed_acyclic())
        ....:     for e in g.edges(sort=False):
        ....:         assert(sccs[e[0]]==sccs[e[1]] or scc_digraph.has_edge(sccs[e[0]],sccs[e[1]]))
        ....:         assert(sccs[e[0]] >= sccs[e[1]])

    Immutable digraphs::

        sage: from sage.graphs.base.static_sparse_graph import strongly_connected_components_digraph
        sage: G = digraphs.RandomDirectedGNP(10, .4)
        sage: G._backend
        <sage.graphs.base.sparse_graph.SparseGraphBackend ...>
        sage: H = DiGraph(G, immutable=True)
        sage: H._backend
        <sage.graphs.base.static_sparse_backend.StaticSparseBackend ...>
        sage: A = strongly_connected_components_digraph(G)[0]
        sage: B = strongly_connected_components_digraph(H)[0]
        sage: A.is_isomorphic(B)
        True
    """
    from sage.graphs.digraph import DiGraph
    if not isinstance(G, DiGraph):
        raise ValueError("G must be a DiGraph.")

    cdef list int_to_vertex
    cdef StaticSparseCGraph cg
    cdef short_digraph g, scc_g
    if isinstance(G, StaticSparseBackend):
        cg = <StaticSparseCGraph> G._cg
        g = <short_digraph> cg.g
        int_to_vertex = cg._vertex_to_labels
    else:
        int_to_vertex = list(G)
        init_short_digraph(g, G, edge_labelled=False, vertex_list=int_to_vertex)

    cdef MemoryAllocator mem = MemoryAllocator()
    cdef int * scc = <int*> mem.malloc(g.n * sizeof(int))
    cdef int i, j, nscc
    cdef list edges = []

    sig_on()
    nscc = tarjan_strongly_connected_components_C(g, scc)
    strongly_connected_components_digraph_C(g, nscc, scc, scc_g)

    output = DiGraph(nscc)

    for i in range(scc_g.n):
        for j in range(scc_g.neighbors[i + 1] - scc_g.neighbors[i]):
            edges.append((i, scc_g.neighbors[i][j]))
    output.add_edges(edges)
    sig_off()

    if not isinstance(G, StaticSparseBackend):
        free_short_digraph(g)
    free_short_digraph(scc_g)

    return output, {v: scc[i] for i, v in enumerate(int_to_vertex)}


cdef strongly_connected_component_containing_vertex(short_digraph g, short_digraph g_reversed, int v, bitset_t scc):
    """
    Feed ``scc`` with the vertices in the strongly connected component of ``v``.
    """
    # Computing the set of vertices that can be reached from v in g
    can_be_reached_from(g, v, scc)
    # Computing the set of vertices that can be reached from v in g *reversed*
    cdef bitset_t scc_reversed
    bitset_init(scc_reversed, g.n)
    can_be_reached_from(g_reversed, v, scc_reversed)
    # The scc containing v is the intersection of both sets
    bitset_intersection(scc, scc, scc_reversed)


cdef void free_short_digraph(short_digraph g) noexcept:
    """
    Free the resources used by ``g``
    """
    sig_free(g.edges)
    sig_free(g.neighbors)
    if g.edge_labels != NULL:
        cpython.Py_XDECREF(g.edge_labels)


def triangles_count(G):
    r"""
    Return the number of triangles containing `v`, for every `v`.

    INPUT:

    - ``G`` -- a graph

    EXAMPLES::

        sage: from sage.graphs.base.static_sparse_graph import triangles_count
        sage: triangles_count(graphs.PetersenGraph())
        {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0}
        sage: sum(triangles_count(graphs.CompleteGraph(15)).values()) == 3*binomial(15,3)           # needs sage.symbolic
        True

    TESTS:

    Immutable graphs::

        sage: from sage.graphs.base.static_sparse_graph import triangles_count
        sage: G = graphs.RandomGNP(10, .7)
        sage: G._backend
        <sage.graphs.base.sparse_graph.SparseGraphBackend ...>
        sage: H = Graph(G, immutable=True)
        sage: H._backend
        <sage.graphs.base.static_sparse_backend.StaticSparseBackend ...>
        sage: triangles_count(G) == triangles_count(H)
        True
    """
    from sage.rings.integer import Integer
    G._scream_if_not_simple()

    # g is a copy of G. If G is internally a static sparse graph, we use it.
    cdef list int_to_vertex
    cdef StaticSparseCGraph cg
    cdef short_digraph g
    if isinstance(G, StaticSparseBackend):
        cg = <StaticSparseCGraph> G._cg
        g = <short_digraph> cg.g
        int_to_vertex = cg._vertex_to_labels
    else:
        int_to_vertex = list(G)
        init_short_digraph(g, G, edge_labelled=False, vertex_list=int_to_vertex)

    cdef uint64_t * count = <uint64_t *> check_calloc(G.order(), sizeof(uint64_t))

    cdef uint64_t tmp_count = 0
    cdef uint32_t u, v, i
    cdef uint32_t * p1
    cdef uint32_t * p2

    for u in range(<uint32_t>g.n):
        for i in range(<uint32_t>out_degree(g, u)):
            v = g.neighbors[u][i]
            if v <= u:
                continue

            # Size of [N(u) inter N(v)]. Both are sorted lists.
            p1 = g.neighbors[u]
            p2 = g.neighbors[v]
            tmp_count = 0
            while (p1 < g.neighbors[u + 1] and p2 < g.neighbors[v + 1]):
                if p1[0] == p2[0]:
                    tmp_count += 1
                    p1 += 1
                    p2 += 1
                elif p1[0] < p2[0]:
                    p1 += 1
                else:
                    p2 += 1

            count[u] += tmp_count
            count[v] += tmp_count

    ans = {w: Integer(count[i] // 2) for i, w in enumerate(int_to_vertex)}

    if not isinstance(G, StaticSparseBackend):
        free_short_digraph(g)
    sig_free(count)
    return ans


def spectral_radius(G, prec=1e-10):
    r"""
    Return an interval of floating point number that encloses the spectral
    radius of this graph

    The input graph ``G`` must be *strongly connected*.

    INPUT:

    - ``prec`` -- (default: ``1e-10``) an upper bound for the relative precision
      of the interval

    The algorithm is iterative and uses an inequality valid for nonnegative
    matrices. Namely, if `A` is a nonnegative square matrix with
    Perron-Frobenius eigenvalue `\lambda` then the following inequality is valid
    for any vector `x`

    .. MATH::

        \min_i  \frac{(Ax)_i}{x_i} \leq \lambda \leq \max_i \frac{(Ax)_i}{x_i}

    .. NOTE::

        The speed of convergence of the algorithm is governed by the spectral
        gap (the distance to the second largest modulus of other eigenvalues).
        If this gap is small, then this function might not be appropriate.

        The algorithm is not smart and not parallel! It uses basic interval
        arithmetic and native floating point arithmetic.

    EXAMPLES::

        sage: from sage.graphs.base.static_sparse_graph import spectral_radius

        sage: G = DiGraph([(0,0),(0,1),(1,0)], loops=True)
        sage: phi = (RR(1) + RR(5).sqrt() ) / 2
        sage: phi  # abs tol 1e-14
        1.618033988749895
        sage: e_min, e_max = spectral_radius(G, 1e-14)
        sage: e_min, e_max     # abs tol 1e-14
        (1.618033988749894, 1.618033988749896)
        sage: (e_max - e_min)  # abs tol 1e-14
        1e-14
        sage: e_min < phi < e_max
        True

    This function also works for graphs::

        sage: G = Graph([(0,1),(0,2),(1,2),(1,3),(2,4),(3,4)])
        sage: e_min, e_max = spectral_radius(G, 1e-14)
        sage: e = max(G.adjacency_matrix().charpoly().roots(AA, multiplicities=False))  # needs sage.modules sage.rings.number_field
        sage: e_min < e < e_max                                                         # needs sage.modules sage.rings.number_field sage.symbolic
        True

        sage: G.spectral_radius()  # abs tol 1e-9
        (2.48119430408, 2.4811943041)

    A larger example::

        sage: # needs sage.modules
        sage: G = DiGraph()
        sage: G.add_edges((i,i+1) for i in range(200))
        sage: G.add_edge(200,0)
        sage: G.add_edge(1,0)
        sage: e_min, e_max = spectral_radius(G, 0.00001)
        sage: p = G.adjacency_matrix(sparse=True).charpoly()
        sage: p
        x^201 - x^199 - 1
        sage: r = p.roots(AA, multiplicities=False)[0]                                  # needs sage.rings.number_field
        sage: e_min < r < e_max                                                         # needs sage.rings.number_field
        True

    A much larger example::

        sage: G = DiGraph(100000)
        sage: r = list(range(100000))
        sage: while not G.is_strongly_connected():
        ....:     shuffle(r)
        ....:     G.add_edges(enumerate(r), loops=False)
        sage: spectral_radius(G, 1e-10)  # random  # long time
        (1.9997956006500042, 1.9998043797692782)

    The algorithm takes care of multiple edges::

        sage: G = DiGraph(2,loops=True,multiedges=True)
        sage: G.add_edges([(0,0),(0,0),(0,1),(1,0)])
        sage: spectral_radius(G, 1e-14)  # abs tol 1e-14
        (2.414213562373094, 2.414213562373095)
        sage: max(G.adjacency_matrix().eigenvalues(AA))                                 # needs sage.modules sage.rings.number_field
        2.414213562373095?

    Some bipartite graphs::

        sage: G = Graph([(0,1),(0,3),(2,3)])
        sage: G.spectral_radius()  # abs tol 1e-10
        (1.6180339887253428, 1.6180339887592732)

        sage: G = DiGraph([(0,1),(0,3),(2,3),(3,0),(1,0),(1,2)])
        sage: G.spectral_radius() # abs tol 1e-10
        (1.5537739740270458, 1.553773974033029)

        sage: G = graphs.CompleteBipartiteGraph(1,3)
        sage: G.spectral_radius()  # abs tol 1e-10
        (1.7320508075688772, 1.7320508075688774)

    TESTS::

        sage: from sage.graphs.base.static_sparse_graph import spectral_radius

        sage: Graph(1).spectral_radius()
        (0.0, 0.0)
        sage: Graph([(0,0)], loops=True).spectral_radius()
        (1.0, 1.0)

        sage: spectral_radius(Graph([(0,1),(0,2)]), 1e-20)
        Traceback (most recent call last):
        ...
        ValueError: precision (=1.00000000000000e-20) is too small

        sage: for _ in range(100):                                                      # needs sage.modules sage.rings.number_field
        ....:     G = digraphs.RandomDirectedGNM(10,35)
        ....:     if not G.is_strongly_connected():
        ....:         continue
        ....:     e = max(G.adjacency_matrix().charpoly().roots(AA,multiplicities=False))
        ....:     e_min, e_max = G.spectral_radius(1e-13)
        ....:     assert e_min < e < e_max

        sage: spectral_radius(Graph(), 1e-10)
        Traceback (most recent call last):
        ...
        ValueError: empty graph

        sage: G = DiGraph([(0,1),(1,2),(2,0),(2,0)], multiedges=True)
        sage: G.spectral_radius()
        Traceback (most recent call last):
        ...
        ValueError: the graph must be aperiodic
    """
    if not G:
        raise ValueError("empty graph")
    if G.is_directed():
        if not G.is_strongly_connected():
            raise ValueError("G must be strongly connected")
    elif not G.is_connected():
        raise ValueError("G must be connected")

    cdef double e_min, e_max

    if G.num_verts() == 1:
        e_min = e_max = G.num_edges()
        return (e_min, e_max)

    is_bipartite, colors = G.is_bipartite(certificate=True)
    if is_bipartite:
        # NOTE: for bipartite graph there are two eigenvalues of maximum modulus
        # and the iteration is likely to reach a cycle of length 2 and hence the
        # algorithm never terminate. Here we compute the "square" reduced to
        # one component of the bipartition.
        from sage.graphs.digraph import DiGraph
        H = DiGraph(loops=True, multiedges=True)
        if G.is_directed():
            neighbors_iterator = G.neighbor_out_iterator
        else:
            neighbors_iterator = G.neighbor_iterator
        for u0 in G:
            if colors[u0] == 0:
                for u1 in neighbors_iterator(u0):
                    for u2 in neighbors_iterator(u1):
                        H.add_edge(u0, u2)

        e_min, e_max = spectral_radius(H, prec)
        return sqrt(e_min), sqrt(e_max)

    cdef double c_prec = prec
    if 1+c_prec/2 == 1:
        raise ValueError("precision (={!r}) is too small".format(prec))

    # test if the graph is aperiodic
    from sage.graphs.digraph import DiGraph
    if not DiGraph(G).is_aperiodic():
        raise ValueError("the graph must be aperiodic")

    # make a copy of G if needed to obtain a static sparse graph
    # NOTE: the following potentially copies the labels of the graph which is
    # completely useless for the computation!
    cdef short_digraph g
    G = G.copy(immutable=True)
    g[0] = (<StaticSparseCGraph> (<StaticSparseBackend> G._backend)._cg).g[0]

    cdef size_t n = <size_t>g.n
    cdef size_t m = <size_t>g.m
    cdef uint32_t ** neighbors = g.neighbors

    # v1 and v2 are two arrays of length n, allocated as one array
    # of length 2n for efficiency.
    cdef double * vmem = <double *>check_allocarray(2*n, sizeof(double))
    cdef double * v1 = vmem
    cdef double * v2 = vmem + n
    cdef double * v3

    cdef size_t i
    cdef uint32_t *p
    cdef double s

    for i in range(n):
        v1[i] = 1
    s = n

    cdef int old_rounding = fegetround()
    e_max = m
    e_min = 0
    try:
        sig_on()
        while (e_max - e_min) > e_max * c_prec:
            # renormalize
            s = n/s
            for i in range(n):
                v1[i] *= s

            # computing e_max (with upward rounding)
            e_max = 0
            fesetround(FE_UPWARD)
            p = neighbors[0]
            for i in range(n):
                v2[i] = 0
                while p < neighbors[i + 1]:
                    v2[i] += v1[p[0]]
                    p += 1
                e = v2[i] / v1[i]
                if e > e_max:
                    e_max = e

            # computing e_min (with downward rounding)
            e_min = m
            fesetround(FE_DOWNWARD)
            p = neighbors[0]
            for i in range(n):
                v2[i] = 0
                while p < neighbors[i + 1]:
                    v2[i] += v1[p[0]]
                    p += 1
                s += v2[i]
                e = v2[i] / v1[i]
                if e < e_min:
                    e_min = e

            # computing the next vector (with nearest rounding)
            fesetround(FE_TONEAREST)
            s = 0
            p = neighbors[0]
            for i in range(n):
                v2[i] = 0
                while p < neighbors[i + 1]:
                    v2[i] += v1[p[0]]
                    p += 1
                s += v2[i]
            v3 = v1
            v1 = v2
            v2 = v3

        sig_off()
    finally:
        # be sure that the rounding is back to default
        fesetround(old_rounding)

        # and that the memory is freed
        sig_free(vmem)

    return (e_min, e_max)
