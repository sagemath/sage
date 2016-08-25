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
                print("There is an edge from {} to {}".format(i, g.neighbors[i][j]))

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

At the moment, it is used in the :mod:`sage.graphs.distances_all_pairs` module,
and in the
:meth:`~sage.graphs.digraph.DiGraph.strongly_connected_components` method.

Python functions
----------------

These functions are available so that Python modules from Sage can call the
Cython routines this module implements (as they can not directly call methods
with C arguments).
"""

#*****************************************************************************
#       Copyright (C) 2010 Nathann Cohen <nathann.cohen@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

include "sage/data_structures/bitset.pxi"
cimport cpython
from libc.string cimport memset
from libc.limits cimport INT_MAX
from sage.graphs.base.c_graph cimport CGraph
from static_sparse_backend cimport StaticSparseCGraph
from static_sparse_backend cimport StaticSparseBackend
from sage.ext.memory_allocator cimport MemoryAllocator
include "cysignals/memory.pxi"
from libcpp.vector cimport vector

cdef extern from "fenv.h":
    int FE_TONEAREST
    int FE_UPWARD
    int FE_DOWNWARD
    int FE_TOWARDZERO
    int fegetround ()
    int fesetround (int)

cdef int init_short_digraph(short_digraph g, G, edge_labelled = False) except -1:
    r"""
    Initializes ``short_digraph g`` from a Sage (Di)Graph.

    If ``G`` is a ``Graph`` objet (and not a ``DiGraph``), an edge between two
    vertices `u` and `v` is replaced by two arcs in both directions.
    """
    g.edge_labels = NULL

    if G.order() >= INT_MAX:
        raise ValueError("This structure can handle at most "+str(INT_MAX)+" vertices !")
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

    g.edges = <uint32_t *> sig_malloc(n_edges*sizeof(uint32_t))
    if g.edges == NULL:
        raise ValueError("Problem while allocating memory (edges)")

    g.neighbors = <uint32_t **> sig_malloc((1+<int>g.n)*sizeof(uint32_t *))
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

    dst.edges = <uint32_t *> sig_malloc(n_edges(src)*sizeof(uint32_t))
    if dst.edges == NULL:
        raise ValueError("Problem while allocating memory (edges)")

    dst.neighbors = <uint32_t **> sig_malloc((src.n+1)*sizeof(uint32_t *))
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

    cdef int * in_degree = <int *> sig_malloc(src.n*sizeof(int))
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
    sig_free(in_degree)

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
    cdef int * stack = <int *> sig_malloc(g.n*sizeof(int))
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

    sig_free(stack)

cdef int tarjan_strongly_connected_components_C(short_digraph g, int *scc):
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
    cdef int u,v,w, n = g.n, current_index = 0, currentscc = 0
    cdef int *index = <int *> mem.malloc(n * sizeof(int))
    cdef int *pred = <int *> mem.malloc(n * sizeof(int))
    cdef int *lowlink = <int *> mem.malloc(n * sizeof(int))
    cdef int *dfs_stack = <int *> mem.malloc((n_edges(g) + 1) * sizeof(int))
    cdef int dfs_stack_end
    cdef int *scc_stack = <int *> mem.malloc(n * sizeof(int)) # Used to keep track of which nodes are in the "current" SCC
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

        while dfs_stack_end > 0:
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
                while p_tmp<g.neighbors[v+1]:
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
    :wikipedia:`Wikipedia article on Tarjan's algorithm <Tarjan's_strongly_connected_components_algorithm>`.

    EXAMPLE::

        sage: from sage.graphs.base.static_sparse_graph import tarjan_strongly_connected_components
        sage: tarjan_strongly_connected_components(digraphs.Path(3))
        [[2], [1], [0]]
        sage: D = DiGraph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
        sage: D.connected_components()
        [[0, 1, 2, 3], [4, 5, 6]]
        sage: D = DiGraph( { 0 : [1, 3], 1 : [2], 2 : [3], 4 : [5, 6], 5 : [6] } )
        sage: D.strongly_connected_components()
        [[3], [2], [1], [0], [6], [5], [4]]
        sage: D.add_edge([2,0])
        sage: D.strongly_connected_components()
        [[3], [0, 1, 2], [6], [5], [4]]
        sage: D = DiGraph([('a','b'), ('b','c'), ('c', 'd'), ('d', 'b'), ('c', 'e')])
        sage: D.strongly_connected_components()
        [['e'], ['b', 'c', 'd'], ['a']]

    TESTS:

    Checking that the result is correct::

        sage: from sage.graphs.base.static_sparse_graph import tarjan_strongly_connected_components
        sage: import random
        sage: for i in range(10):                                     # long
        ....:     n = random.randint(2,20)
        ....:     m = random.randint(1, n*(n-1))
        ....:     g = digraphs.RandomDirectedGNM(n,m)
        ....:     sccs = tarjan_strongly_connected_components(g)
        ....:     for scc in sccs:
        ....:         scc_check = g.strongly_connected_component_containing_vertex(scc[0])
        ....:         assert(sorted(scc) == sorted(scc_check))

    Checking against NetworkX::

        sage: import networkx
        sage: for i in range(10):                                     # long
        ....:      g = digraphs.RandomDirectedGNP(100,.05)
        ....:      h = g.networkx_graph()
        ....:      scc1 = g.strongly_connected_components()
        ....:      scc2 = networkx.strongly_connected_components(h)
        ....:      s1 = Set(map(Set,scc1))
        ....:      s2 = Set(map(Set,scc2))
        ....:      if s1 != s2:
        ....:          print("Ooch !")
    """
    from sage.graphs.digraph import DiGraph

    if not isinstance(G, DiGraph):
        raise ValueError("G must be a DiGraph.")

    sig_on()
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef short_digraph g
    init_short_digraph(g, G)
    cdef int * scc = <int*> mem.malloc(g.n * sizeof(int))
    cdef int nscc = tarjan_strongly_connected_components_C(g, scc)
    cdef int i
    cdef list output = list(list() for i in range(nscc)) # We cannot use [] here

    for i,v in enumerate(G.vertices()):
        output[scc[i]].append(v)
    sig_off()
    return output

cdef void strongly_connected_components_digraph_C(short_digraph g, int nscc, int *scc, short_digraph output):
    r"""
    Computes the strongly connected components (SCCs) digraph of `g`.

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
    cdef int v, w, i
    cdef int tmp = nscc + 1
    cdef vector[vector[int]] scc_list = vector[vector[int]](nscc, vector[int]())
    cdef vector[vector[int]] sons = vector[vector[int]](nscc + 1, vector[int]())
    cdef vector[int].iterator iter
    cdef short *neighbors = <short *> mem.calloc(nscc, sizeof(short))
    cdef long m = 0
    cdef uint32_t degv
    cdef uint32_t *p_tmp

    for v in range(nscc):
        scc_list[v] = vector[int]()
        sons[v] = vector[int]()
    sons[nscc] = vector[int]()

    for i in range(g.n):
        scc_list[scc[i]].push_back(i)

    for v in range(nscc):
        for i in range(scc_list[v].size()):
            p_tmp = g.neighbors[scc_list[v][i]]
            while p_tmp<g.neighbors[scc_list[v][i]+1]:
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

    output.neighbors = <uint32_t **> check_allocarray((1+<int>output.n), sizeof(uint32_t *))

    if m == 0:
        output.edges = NULL
        for v in range(1,nscc + 1):
            output.neighbors[v] = NULL

    output.edges = <uint32_t *> check_allocarray(m, sizeof(uint32_t))
    output.neighbors[0] = output.edges

    for v in range(1,nscc + 1):
        degv = sons[v].size()
        output.neighbors[v] = output.neighbors[v-1] + sons[v-1].size()
        for i in range(sons[v].size()):
            output.neighbors[v][i] = sons[v][i]

def strongly_connected_components_digraph(G):
    r"""
    Returns the digraph of the strongly connected components (SCCs).

    This routine is used to test ``strongly_connected_components_digraph_C``,
    but it is not used by the Sage digraph. It outputs a pair ``[g_scc,scc]``,
    where ``g_scc`` is the SCC digraph of g, ``scc`` is a dictionary associating
    to each vertex ``v`` the number of the SCC of ``v``, as it appears in
    ``g_scc``.

    EXAMPLE::

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
        ....:     for e in g.edges():
        ....:         assert(sccs[e[0]]==sccs[e[1]] or scc_digraph.has_edge(sccs[e[0]],sccs[e[1]]))
        ....:         assert(sccs[e[0]] >= sccs[e[1]])
    """
    from sage.graphs.digraph import DiGraph
    if not isinstance(G, DiGraph):
        raise ValueError("G must be a DiGraph.")

    cdef MemoryAllocator mem = MemoryAllocator()
    cdef short_digraph g, scc_g
    init_short_digraph(g, G)
    cdef int * scc = <int*> mem.malloc(g.n * sizeof(int))
    cdef int i, j, nscc
    cdef list edges = []

    sig_on()
    nscc = tarjan_strongly_connected_components_C(g, scc)
    strongly_connected_components_digraph_C(g, nscc, scc, scc_g)

    output = DiGraph(nscc)

    for i in range(scc_g.n):
        for j in range(scc_g.neighbors[i+1]-scc_g.neighbors[i]):
            edges.append((i, scc_g.neighbors[i][j]))
    output.add_edges(edges)
    sig_off()
    return output, {v:scc[i] for i,v in enumerate(G.vertices())}


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
        sig_free(g.edges)

    if g.neighbors != NULL:
        sig_free(g.neighbors)

    if g.edge_labels != NULL:
        cpython.Py_XDECREF(g.edge_labels)

def triangles_count(G):
    r"""
    Return the number of triangles containing `v`, for every `v`.

    INPUT:

    - `G`-- a graph

    EXAMPLE::

        sage: from sage.graphs.base.static_sparse_graph import triangles_count
        sage: triangles_count(graphs.PetersenGraph())
        {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0}
        sage: sum(triangles_count(graphs.CompleteGraph(15)).values()) == 3*binomial(15,3)
        True
    """
    from sage.rings.integer import Integer
    G._scream_if_not_simple()

    # g is a copy of G. If G is internally a static sparse graph, we use it.
    cdef short_digraph g
    G = G.copy(immutable=True)

    cdef uint64_t * count = <uint64_t *> check_calloc(G.order(), sizeof(uint64_t))
    g[0] = (<StaticSparseCGraph> (<StaticSparseBackend> G._backend)._cg).g[0]

    cdef uint64_t tmp_count = 0
    cdef uint32_t u,v,i
    cdef uint32_t * p1
    cdef uint32_t * p2

    for u in range(g.n):
        for i in range(out_degree(g,u)):
            v = g.neighbors[u][i]
            if v<=u:
                continue

            # Size of [N(u) inter N(v)]. Both are sorted lists.
            p1 = g.neighbors[u]
            p2 = g.neighbors[v]
            tmp_count = 0
            while (p1 < g.neighbors[u+1] and p2 < g.neighbors[v+1]):
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

    ans = {w:Integer(count[i]/2)
           for i,w in enumerate(G.vertices())}

    sig_free(count)
    return ans

def spectral_radius(G, prec=1e-10):
    r"""
    Return an interval of floating point number that encloses the spectral
    radius of this graph

    The input graph ``G`` must be *strongly connected*.

    INPUT:

    - ``prec`` -- (default ``1e-10``) an upper bound for the relative precision
      of the interval

    The algorithm is iterative and uses an inequality valid for non-negative
    matrices. Namely, if `A` is a non-negative square matrix with
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
        sage: e = max(G.adjacency_matrix().charpoly().roots(AA, multiplicities=False))
        sage: e_min < e < e_max
        True

        sage: G.spectral_radius()  # abs tol 1e-9
        (2.48119430408, 2.4811943041)

    A larger example::

        sage: G = DiGraph()
        sage: G.add_edges((i,i+1) for i in range(200))
        sage: G.add_edge(200,0)
        sage: G.add_edge(1,0)
        sage: e_min, e_max = spectral_radius(G, 0.00001)
        sage: p = G.adjacency_matrix(sparse=True).charpoly()
        sage: p
        x^201 - x^199 - 1
        sage: r = p.roots(AA, multiplicities=False)[0]
        sage: e_min < r < e_max
        True

    A much larger example::

        sage: G = DiGraph(100000)
        sage: r = range(100000)
        sage: while not G.is_strongly_connected():
        ....:     shuffle(r)
        ....:     G.add_edges(enumerate(r))
        sage: spectral_radius(G, 1e-10)  # random
        (1.9997956006500042, 1.9998043797692782)

    The algorithm takes care of multiple edges::

        sage: G = DiGraph(2,loops=True,multiedges=True)
        sage: G.add_edges([(0,0),(0,0),(0,1),(1,0)])
        sage: spectral_radius(G, 1e-14)  # abs tol 1e-14
        (2.414213562373094, 2.414213562373095)
        sage: max(G.adjacency_matrix().eigenvalues(AA))
        2.414213562373095?

    TESTS::

        sage: spectral_radius(G, 1e-20)
        Traceback (most recent call last):
        ...
        ValueError: precision (=1.00000000000000e-20) is too small

        sage: for _ in range(100):
        ....:     G = digraphs.RandomDirectedGNM(10,35)
        ....:     if not G.is_strongly_connected():
        ....:         continue
        ....:     e = max(G.adjacency_matrix().charpoly().roots(AA,multiplicities=False))
        ....:     e_min, e_max = G.spectral_radius(1e-13)
        ....:     assert e_min < e < e_max
    """
    if G.is_directed():
        if not G.is_strongly_connected():
            raise ValueError("G must be strongly connected")
    elif not G.is_connected():
        raise ValueError("G must be connected")

    cdef double c_prec = prec
    if 1+c_prec/2 == 1:
        raise ValueError("precision (={!r}) is too small".format(prec))

    # make a copy of G if needed to obtain a static sparse graph
    # NOTE: the following potentially copies the labels of the graph which is
    # completely useless for the computation!
    cdef short_digraph g
    G = G.copy(immutable=True)
    g[0] = (<StaticSparseCGraph> (<StaticSparseBackend> G._backend)._cg).g[0]

    cdef long n = g.n
    cdef long m = g.m
    cdef uint32_t ** neighbors = g.neighbors

    cdef double * v1 = <double *> sig_malloc(n * sizeof(double))
    cdef double * v2 = <double *> sig_malloc(n * sizeof(double))
    cdef double * v3
    if v1 == NULL or v2 == NULL:
        sig_free(v1)
        sig_free(v2)
        raise MemoryError

    cdef size_t i
    cdef uint32_t *p
    cdef double e_min, e_max
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
                while p < neighbors[i+1]:
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
                while p < neighbors[i+1]:
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
                while p < neighbors[i+1]:
                    v2[i] += v1[p[0]]
                    p += 1
                s += v2[i]
            v3 = v1; v1 = v2; v2 = v3

        sig_off()
    finally:
        # be sure that the rounding is back to default
        fesetround(old_rounding)

        # and that the memory is freed
        sig_free(v1)
        sig_free(v2)

    return (e_min, e_max)
