r"""
Fast dense graphs

For an overview of graph data structures in sage, see
:mod:`~sage.graphs.base.overview`.

Usage Introduction
------------------

::

    sage: from sage.graphs.base.dense_graph import DenseGraph

Dense graphs are initialized as follows::

    sage: D = DenseGraph(nverts = 10, extra_vertices = 10)

This example initializes a dense graph with room for twenty vertices, the first
ten of which are in the graph. In general, the first ``nverts`` are "active."
For example, see that 9 is already in the graph::

    sage: D._num_verts()
    10
    sage: D.add_vertex(9)
    9
    sage: D._num_verts()
    10

But 10 is not, until we add it::

    sage: D._num_verts()
    10
    sage: D.add_vertex(10)
    10
    sage: D._num_verts()
    11

You can begin working right away as follows::

    sage: D.add_arc(0,1)
    sage: D.add_arc(1,2)
    sage: D.add_arc(1,0)
    sage: D.has_arc(7,3)
    False
    sage: D.has_arc(0,1)
    True
    sage: D.in_neighbors(1)
    [0]
    sage: D.out_neighbors(1)
    [0, 2]
    sage: D.del_all_arcs(0,1)
    sage: D.has_arc(0,1)
    False
    sage: D.has_arc(1,2)
    True
    sage: D.del_vertex(7)
    sage: D.has_arc(7,3)
    False
    sage: D._num_verts()
    10
    sage: D._num_arcs()
    2

Dense graphs do not support multiple or labeled edges.

::

    sage: T = DenseGraph(nverts = 3, extra_vertices = 2)
    sage: T.add_arc(0,1)
    sage: T.add_arc(1,2)
    sage: T.add_arc(2,0)
    sage: T.has_arc(0,1)
    True

::

    sage: for _ in range(10): D.add_arc(5,4)
    sage: D.has_arc(5,4)
    True

Dense graphs are by their nature directed. As of this writing, you need to do
operations in pairs to treat the undirected case (or use a backend or a Sage
graph)::

    sage: T.has_arc(1,0)
    False

The curious developer is encouraged to check out the ``unsafe`` functions,
which do not check input but which run in pure C.

Underlying Data Structure
-------------------------

The class ``DenseGraph`` contains the following variables which are inherited
from ``CGraph`` (for explanation, refer to the documentation there)::

        cdef int num_verts
        cdef int num_arcs
        cdef int *in_degrees
        cdef int *out_degrees
        cdef bitset_t active_vertices

It also contains the following variables::

        cdef int radix_div_shift
        cdef int radix_mod_mask
        cdef int num_longs
        cdef unsigned long *edges

The array ``edges`` is a series of bits which are turned on or off, and due to
this, dense graphs only support graphs without edge labels and with no multiple
edges. The ints ``radix_div_shift`` and ``radix_mod_mask`` are simply for doing
efficient division by powers of two, and ``num_longs`` stores the length of the
``edges`` array. Recall that this length reflects the number of available
vertices, not the number of "actual" vertices. For more details about this,
refer to the documentation for ``CGraph``.
"""

#*******************************************************************************
#        Copyright (C) 2008-9 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*******************************************************************************

include 'sage/data_structures/bitset.pxi'

cdef class DenseGraph(CGraph):
    """
    Compiled dense graphs.

    ::

        sage: from sage.graphs.base.dense_graph import DenseGraph

    Dense graphs are initialized as follows::

        sage: D = DenseGraph(nverts = 10, extra_vertices = 10)

    INPUT:

     - ``nverts`` - non-negative integer, the number of vertices.
     - ``extra_vertices`` - non-negative integer (default: 0), how many extra
        vertices to allocate.
     - ``verts`` - optional list of vertices to add
     - ``arcs`` - optional list of arcs to add

    The first ``nverts`` are created as vertices of the graph, and the next
    ``extra_vertices`` can be freely added without reallocation. See top level
    documentation for more details. The input ``verts`` and ``arcs`` are mainly
    for use in pickling.

    """

    def __cinit__(self, int nverts, int extra_vertices = 10, verts = None, arcs = None):
        """
        Allocation and initialization happen in one place.

        Memory usage is

        O(  (nverts + extra_vertices)^2  ).

        EXAMPLE::

            sage: from sage.graphs.base.dense_graph import DenseGraph
            sage: D = DenseGraph(nverts = 10, extra_vertices = 10)
        """
        if nverts == 0 and extra_vertices == 0:
            raise RuntimeError('Dense graphs must allocate space for vertices!')
        cdef int radix = sizeof(unsigned long) << 3
        self.radix_mod_mask = radix - 1
        cdef int i = 0
        while ((<unsigned long>1)<<i) & self.radix_mod_mask:
            i += 1
        self.radix_div_shift = i
        self.num_verts = nverts
        cdef int total_verts = nverts + extra_vertices
        self.num_arcs = 0

        i = total_verts >> self.radix_div_shift
        if total_verts & self.radix_mod_mask:
            i += 1
        self.num_longs = i

        self.edges = <unsigned long *> sage_malloc(total_verts * self.num_longs * sizeof(unsigned long))
        self.in_degrees = <int *> sage_malloc(total_verts * sizeof(int))
        self.out_degrees = <int *> sage_malloc(total_verts * sizeof(int))

        if not self.edges or not self.in_degrees or not self.out_degrees:
            if self.edges: sage_free(self.edges)
            if self.in_degrees: sage_free(self.in_degrees)
            if self.out_degrees: sage_free(self.out_degrees)
            raise MemoryError
        for i from 0 <= i < self.num_longs * total_verts:
            self.edges[i] = 0
        for i from 0 <= i < total_verts:
            self.in_degrees[i] = 0
            self.out_degrees[i] = 0
        bitset_init(self.active_vertices, total_verts)
        bitset_set_first_n(self.active_vertices, self.num_verts)

        if verts is not None:
            self.add_vertices(verts)

        if arcs is not None:
            for u,v in arcs:
                self.add_arc(u,v)

    def __dealloc__(self):
        """
        New and dealloc are both tested at class level.
        """
        sage_free(self.edges)
        sage_free(self.in_degrees)
        sage_free(self.out_degrees)
        bitset_free(self.active_vertices)

    cpdef realloc(self, int total_verts):
        """
        Reallocate the number of vertices to use, without actually adding any.

        INPUT:

         - ``total`` - integer, the total size to make the array

        Returns -1 and fails if reallocation would destroy any active vertices.

        EXAMPLES::

            sage: from sage.graphs.base.dense_graph import DenseGraph
            sage: D = DenseGraph(nverts=4, extra_vertices=4)
            sage: D.current_allocation()
            8
            sage: D.add_vertex(6)
            6
            sage: D.current_allocation()
            8
            sage: D.add_vertex(10)
            10
            sage: D.current_allocation()
            16
            sage: D.add_vertex(40)
            Traceback (most recent call last):
            ...
            RuntimeError: Requested vertex is past twice the allocated range: use realloc.
            sage: D.realloc(50)
            sage: D.add_vertex(40)
            40
            sage: D.current_allocation()
            50
            sage: D.realloc(30)
            -1
            sage: D.current_allocation()
            50
            sage: D.del_vertex(40)
            sage: D.realloc(30)
            sage: D.current_allocation()
            30

        """
        cdef int i, j
        if total_verts == 0:
            raise RuntimeError('Dense graphs must allocate space for vertices!')
        cdef bitset_t bits
        cdef int min_verts, min_longs, old_longs = self.num_longs
        if total_verts < self.active_vertices.size:
            min_verts = total_verts
            min_longs = -1
            bitset_init(bits, self.active_vertices.size)
            bitset_set_first_n(bits, total_verts)
            if not bitset_issubset(self.active_vertices, bits):
                bitset_free(bits)
                return -1
            bitset_free(bits)
        else:
            min_verts = self.active_vertices.size
            min_longs = self.num_longs

        i = total_verts >> self.radix_div_shift
        if total_verts & self.radix_mod_mask:
            i += 1
        self.num_longs = i
        if min_longs == -1: min_longs = self.num_longs

        cdef unsigned long *new_edges = <unsigned long *> sage_malloc(total_verts * self.num_longs * sizeof(unsigned long))

        for i from 0 <= i < min_verts:
            for j from 0 <= j < min_longs:
                new_edges[i*self.num_longs + j] = self.edges[i*old_longs + j]
            for j from min_longs <= j < self.num_longs:
                new_edges[i*self.num_longs + j] = 0
        for i from min_verts <= i < total_verts:
            for j from 0 <= j < self.num_longs:
                new_edges[i*self.num_longs + j] = 0
        sage_free(self.edges)
        self.edges = new_edges

        self.in_degrees = <int *> sage_realloc(self.in_degrees, total_verts * sizeof(int))
        self.out_degrees = <int *> sage_realloc(self.out_degrees, total_verts * sizeof(int))

        cdef int first_limb
        cdef unsigned long zero_gate
        if total_verts > self.active_vertices.size:
            first_limb = (self.active_vertices.size >> self.radix_div_shift)
            zero_gate = (<unsigned long>1) << (self.active_vertices.size & self.radix_mod_mask)
            zero_gate -= 1
            for i from 0 <= i < total_verts:
                self.edges[first_limb] &= zero_gate
                for j from first_limb < j < self.num_longs:
                    self.edges[j] = 0

            for i from self.active_vertices.size <= i < total_verts:
                self.in_degrees[i] = 0
                self.out_degrees[i] = 0
        bitset_realloc(self.active_vertices, total_verts)

    ###################################
    # Unlabeled arc functions
    ###################################

    cdef int add_arc_unsafe(self, int u, int v) except -1:
        """
        Adds arc (u, v) to the graph.

        INPUT:
            u, v -- non-negative integers

        """
        cdef int place = (u * self.num_longs) + (v >> self.radix_div_shift)
        cdef unsigned long word = (<unsigned long>1) << (v & self.radix_mod_mask)
        if not self.edges[place] & word:
            self.in_degrees[v] += 1
            self.out_degrees[u] += 1
            self.num_arcs += 1
            self.edges[place] |= word

    cpdef add_arc(self, int u, int v):
        """
        Adds arc ``(u, v)`` to the graph.

        INPUT:

         - ``u, v`` -- non-negative integers, must be in self

        EXAMPLE::

            sage: from sage.graphs.base.dense_graph import DenseGraph
            sage: G = DenseGraph(5)
            sage: G.add_arc(0,1)
            sage: G.add_arc(4,7)
            Traceback (most recent call last):
            ...
            LookupError: Vertex (7) is not a vertex of the graph.
            sage: G.has_arc(1,0)
            False
            sage: G.has_arc(0,1)
            True

        """
        self.check_vertex(u)
        self.check_vertex(v)
        self.add_arc_unsafe(u,v)

    cdef int has_arc_unsafe(self, int u, int v) except -1:
        """
        Checks whether arc (u, v) is in the graph.

        INPUT:
            u, v -- non-negative integers, must be in self

        OUTPUT:
            0 -- False
            1 -- True

        """
        cdef int place = (u * self.num_longs) + (v >> self.radix_div_shift)
        cdef unsigned long word = (<unsigned long>1) << (v & self.radix_mod_mask)
        return (self.edges[place] & word) >> (v & self.radix_mod_mask)

    cpdef bint has_arc(self, int u, int v) except -1:
        """
        Checks whether arc ``(u, v)`` is in the graph.

        INPUT:
            u, v -- integers

        EXAMPLE::

            sage: from sage.graphs.base.dense_graph import DenseGraph
            sage: G = DenseGraph(5)
            sage: G.add_arc(0,1)
            sage: G.has_arc(1,0)
            False
            sage: G.has_arc(0,1)
            True

        """
        if u < 0 or u >= self.active_vertices.size or not bitset_in(self.active_vertices, u):
            return False
        if v < 0 or v >= self.active_vertices.size or not bitset_in(self.active_vertices, v):
            return False
        return self.has_arc_unsafe(u,v) == 1

    cdef int del_arc_unsafe(self, int u, int v) except -1:
        """
        Deletes the arc from u to v, if it exists.

        INPUT:
            u, v -- non-negative integers, must be in self

        """
        cdef int place = (u * self.num_longs) + (v >> self.radix_div_shift)
        cdef unsigned long word = (<unsigned long>1) << (v & self.radix_mod_mask)
        if self.edges[place] & word:
            self.in_degrees[v] -= 1
            self.out_degrees[u] -= 1
            self.num_arcs -= 1
            self.edges[place] &= ~word

    cpdef del_all_arcs(self, int u, int v):
        """
        Deletes the arc from ``u`` to ``v``.

        INPUT:
         - ``u, v`` - integers

        NOTE:
        The naming of this function is for consistency with ``SparseGraph``. Of
        course, there can be at most one arc for a ``DenseGraph``.

        EXAMPLE::

            sage: from sage.graphs.base.dense_graph import DenseGraph
            sage: G = DenseGraph(5)
            sage: G.add_arc(0,1)
            sage: G.has_arc(0,1)
            True
            sage: G.del_all_arcs(0,1)
            sage: G.has_arc(0,1)
            False

        """
        self.check_vertex(u)
        self.check_vertex(v)
        self.del_arc_unsafe(u,v)

    ###################################
    # Neighbor functions
    ###################################

    cdef int out_neighbors_unsafe(self, int u, int *neighbors, int size) except -2:
        """
        Gives all v such that (u, v) is an arc of the graph.

        INPUT:
            u -- non-negative integer, must be in self
            neighbors -- must be a pointer to an (allocated) integer array
            size -- the length of the array

        OUTPUT:
            nonnegative integer -- the number of v such that (u, v) is an arc
            -1 -- indicates that the array has been filled with neighbors, but
        there were more

        """
        cdef int place = (u * self.num_longs), num_nbrs = 0
        cdef int i, v = 0
        cdef unsigned long word, data
        for i from 0 <= i < self.num_longs:
            data = self.edges[place + i]
            word = 1
            while word:
                if word & data:
                    if num_nbrs == size:
                        return -1
                    neighbors[num_nbrs] = v
                    num_nbrs += 1
                word = word << 1
                v += 1
        return num_nbrs

    cpdef list out_neighbors(self, int u):
        """
        Gives all ``v`` such that ``(u, v)`` is an arc of the graph.

        INPUT:
         - ``u`` - integer

        EXAMPLES::

            sage: from sage.graphs.base.dense_graph import DenseGraph
            sage: G = DenseGraph(5)
            sage: G.add_arc(0,1)
            sage: G.add_arc(1,2)
            sage: G.add_arc(1,3)
            sage: G.out_neighbors(0)
            [1]
            sage: G.out_neighbors(1)
            [2, 3]

        """
        cdef int i, num_nbrs
        self.check_vertex(u)
        if self.out_degrees[u] == 0:
            return []
        cdef int size = self.out_degrees[u]
        cdef int *neighbors = <int *> sage_malloc(size * sizeof(int))
        if not neighbors:
            raise MemoryError
        num_nbrs = self.out_neighbors_unsafe(u, neighbors, size)
        output = [neighbors[i] for i from 0 <= i < num_nbrs]
        sage_free(neighbors)
        return output

    cdef int in_neighbors_unsafe(self, int v, int *neighbors, int size) except -2:
        """
        Gives all u such that (u, v) is an arc of the graph.

        INPUT:
            v -- non-negative integer, must be in self
            neighbors -- must be a pointer to an (allocated) integer array
            size -- the length of the array

        OUTPUT:
            nonnegative integer -- the number of u such that (u, v) is an arc
            -1 -- indicates that the array has been filled with neighbors, but
        there were more

        """
        cdef int place = v >> self.radix_div_shift
        cdef unsigned long word = (<unsigned long>1) << (v & self.radix_mod_mask)
        cdef int i, num_nbrs = 0
        for i from 0 <= i < self.active_vertices.size:
            if self.edges[place + i*self.num_longs] & word:
                if num_nbrs == size:
                    return -1
                neighbors[num_nbrs] = i
                num_nbrs += 1
        return num_nbrs

    cpdef list in_neighbors(self, int v):
        """
        Gives all ``u`` such that ``(u, v)`` is an arc of the graph.

        INPUT:
         - ``v`` - integer

        EXAMPLES::

            sage: from sage.graphs.base.dense_graph import DenseGraph
            sage: G = DenseGraph(5)
            sage: G.add_arc(0,1)
            sage: G.add_arc(3,1)
            sage: G.add_arc(1,3)
            sage: G.in_neighbors(1)
            [0, 3]
            sage: G.in_neighbors(3)
            [1]

        """
        cdef int i, num_nbrs
        self.check_vertex(v)
        if self.in_degrees[v] == 0:
            return []
        cdef int size = self.in_degrees[v]
        cdef int *neighbors = <int *> sage_malloc(size * sizeof(int))
        if not neighbors:
            raise MemoryError
        num_nbrs = self.in_neighbors_unsafe(v, neighbors, size)
        output = [neighbors[i] for i from 0 <= i < num_nbrs]
        sage_free(neighbors)
        return output

##############################
# Further tests. Unit tests for methods, functions, classes defined with cdef.
##############################

def random_stress():
    """
    Randomly search for mistakes in the code.

    DOCTEST (No output indicates that no errors were found)::

        sage: from sage.graphs.base.dense_graph import random_stress
        sage: for _ in xrange(400):
        ...    random_stress()

    """
    cdef int i, j, k, l, n
    cdef DenseGraph Gnew
    num_verts = 10
    # This code deliberately uses random instead of sage.misc.prandom,
    # so that every time it is run it does different tests, instead of
    # doing the same random stress test every time.  (Maybe it should
    # use sage.misc.random_testing?)
    from random import randint
    from sage.graphs.all import DiGraph
    from sage.misc.misc import uniq
    Gnew = DenseGraph(num_verts)
    Gold = DiGraph(num_verts, loops=True, implementation='networkx')
    for n from 0 <= n < 100:
        i = randint(0,num_verts-1)
        j = randint(0,num_verts-1)
        k = randint(0,num_verts-1)
        if k != 0:
            Gold.add_edge(i,j)
            Gnew.add_arc_unsafe(i,j)
        else:
            Gold.delete_edge(i,j)
            Gnew.del_arc_unsafe(i,j)
    if Gnew.num_arcs != Gold.size():
        raise RuntimeError( "NO" )
    for i from 0 <= i < num_verts:
        if Gnew.out_degrees[i] != Gold.out_degree(i):
            raise RuntimeError( "NO" )
        if Gnew.in_degrees[i] != Gold.in_degree(i):
            raise RuntimeError( "NO" )

def _test_adjacency_sequence_out():
    """
    Randomly test the method ``DenseGraph.adjacency_sequence_out()``. No output
    indicates that no errors were found.

    TESTS::

        sage: from sage.graphs.base.dense_graph import _test_adjacency_sequence_out
        sage: _test_adjacency_sequence_out()  # long time
    """
    from sage.graphs.digraph import DiGraph
    from sage.graphs.graph_generators import GraphGenerators
    from sage.misc.prandom import randint, random
    low = 0
    high = 500
    randg = DiGraph(GraphGenerators().RandomGNP(randint(low, high), random()))
    n = randg.order()
    cdef DenseGraph g = DenseGraph(n,
                                   verts=randg.vertices(),
                                   arcs=randg.edges(labels=False))
    assert g._num_verts() == randg.order(), (
        "Graph order mismatch: %s vs. %s" % (g._num_verts(), randg.order()))
    assert g._num_arcs() == randg.size(), (
        "Graph size mismatch: %s vs. %s" % (g._num_arcs(), randg.size()))
    M = randg.adjacency_matrix()
    cdef int *V = <int *>sage_malloc(n * sizeof(int))
    cdef int i = 0
    for v in randg.vertex_iterator():
        V[i] = v
        i += 1
    cdef int *seq = <int *> sage_malloc(n * sizeof(int))
    for 0 <= i < randint(50, 101):
        u = randint(low, n - 1)
        g.adjacency_sequence_out(n, V, u, seq)
        A = [seq[k] for k in range(n)]
        try:
            assert A == list(M[u])
        except AssertionError:
            sage_free(V)
            sage_free(seq)
            raise AssertionError("Graph adjacency mismatch")
    sage_free(seq)
    sage_free(V)

###########################################
# Dense Graph Backend
###########################################

from c_graph cimport CGraphBackend

cdef class DenseGraphBackend(CGraphBackend):
    """
    Backend for Sage graphs using DenseGraphs.

    ::

        sage: from sage.graphs.base.dense_graph import DenseGraphBackend

    This class is only intended for use by the Sage Graph and DiGraph class.
    If you are interested in using a DenseGraph, you probably want to do
    something like the following example, which creates a Sage Graph instance
    which wraps a DenseGraph object::

        sage: G = Graph(30, implementation="c_graph", sparse=False)
        sage: G.add_edges([(0,1), (0,3), (4,5), (9, 23)])
        sage: G.edges(labels=False)
        [(0, 1), (0, 3), (4, 5), (9, 23)]

    Note that Sage graphs using the backend are more flexible than DenseGraphs
    themselves. This is because DenseGraphs (by design) do not deal with Python
    objects::

        sage: G.add_vertex((0,1,2))
        sage: G.vertices()
        [0,
        ...
         29,
         (0, 1, 2)]
        sage: from sage.graphs.base.dense_graph import DenseGraph
        sage: DG = DenseGraph(30)
        sage: DG.add_vertex((0,1,2))
        Traceback (most recent call last):
        ...
        TypeError: an integer is required

    """

    def __init__(self, n, directed=True):
        """
        Initialize a dense graph with n vertices.

        EXAMPLE::

            sage: D = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: D.add_edge(0,1,None,False)
            sage: list(D.iterator_edges(range(9), True))
            [(0, 1, None)]

        """
        self._cg = DenseGraph(n)
        self._cg_rev = None
        self._directed = directed
        self.vertex_labels = {}
        self.vertex_ints = {}

    def add_edge(self, object u, object v, object l, bint directed):
        """
        Adds the edge ``(u,v)`` to self.

        INPUT:

         - ``u,v`` - the vertices of the edge
         - ``l`` - the edge label (ignored)
         - ``directed`` - if False, also add ``(v,u)``

        NOTE:
        The input ``l`` is for consistency with other backends.

        EXAMPLE::

            sage: D = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: D.add_edge(0,1,None,False)
            sage: list(D.iterator_edges(range(9), True))
            [(0, 1, None)]

        """
        if u is None: u = self.add_vertex(None)
        if v is None: v = self.add_vertex(None)

        cdef int u_int = self.check_labelled_vertex(u, 0)
        cdef int v_int = self.check_labelled_vertex(v, 0)

        if directed or u_int == v_int:
            self._cg.add_arc(u_int, v_int)
        else:
            self._cg.add_arc(u_int, v_int)
            self._cg.add_arc(v_int, u_int)

    def add_edges(self, object edges, bint directed):
        """
        Add edges from a list.

        INPUT:

         - ``edges`` - the edges to be added - can either be of the form
           ``(u,v)`` or ``(u,v,l)``
         - ``directed`` - if False, add ``(v,u)`` as well as ``(u,v)``

        EXAMPLE::

            sage: D = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: D.add_edges([(0,1), (2,3), (4,5), (5,6)], False)
            sage: list(D.iterator_edges(range(9), True))
            [(0, 1, None),
             (2, 3, None),
             (4, 5, None),
             (5, 6, None)]

        """
        for e in edges:
            u,v = e[:2]
            self.add_edge(u,v,None,directed)

    def del_edge(self, object u, object v, object l, bint directed):
        """
        Delete edge ``(u,v)``.

        INPUT:

         - ``u,v`` - the vertices of the edge
         - ``l`` - the edge label (ignored)
         - ``directed`` - if False, also delete ``(v,u,l)``

        NOTE:
        The input ``l`` is for consistency with other backends.

        EXAMPLE::

            sage: D = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: D.add_edges([(0,1), (2,3), (4,5), (5,6)], False)
            sage: list(D.iterator_edges(range(9), True))
            [(0, 1, None),
             (2, 3, None),
             (4, 5, None),
             (5, 6, None)]
            sage: D.del_edge(0,1,None,True)
            sage: list(D.iterator_out_edges(range(9), True))
            [(1, 0, None),
             (2, 3, None),
             (3, 2, None),
             (4, 5, None),
             (5, 4, None),
             (5, 6, None),
             (6, 5, None)]

        """
        if not ( self.has_vertex(u) and self.has_vertex(v) ):
            return
        cdef int u_int = self.check_labelled_vertex(u, 0)
        cdef int v_int = self.check_labelled_vertex(v, 0)
        if v is None:
            u, v = u[:2]
        if directed:
            self._cg.del_all_arcs(u_int, v_int)
        else:
            self._cg.del_all_arcs(u_int, v_int)
            self._cg.del_all_arcs(v_int, u_int)

    def get_edge_label(self, object u, object v):
        """
        Returns the edge label for ``(u,v)``. Always None, since dense graphs
        do not support edge labels.

        INPUT:

         - ``u,v`` - the vertices of the edge

        EXAMPLE::

            sage: D = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: D.add_edges([(0,1), (2,3,7), (4,5), (5,6)], False)
            sage: list(D.iterator_edges(range(9), True))
            [(0, 1, None),
             (2, 3, None),
             (4, 5, None),
             (5, 6, None)]
            sage: D.del_edge(0,1,None,True)
            sage: list(D.iterator_out_edges(range(9), True))
            [(1, 0, None),
             (2, 3, None),
             (3, 2, None),
             (4, 5, None),
             (5, 4, None),
             (5, 6, None),
             (6, 5, None)]
            sage: D.get_edge_label(2,3)
            sage: D.get_edge_label(2,4)
            Traceback (most recent call last):
            ...
            LookupError: (2, 4) is not an edge of the graph.

        """
        if not self.has_edge(u, v, None):
            raise LookupError("({0}, {1}) is not an edge of the graph.".format(repr(u), repr(v)))
        return None

    def has_edge(self, object u, object v, object l):
        """
        Returns whether this graph has edge ``(u,v)``.

        NOTE:
        The input ``l`` is for consistency with other backends.

        INPUT:

         - ``u,v`` - the vertices of the edge
         - ``l`` - the edge label (ignored)

        EXAMPLE::

            sage: D = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: D.add_edges([(0,1), (2,3), (4,5), (5,6)], False)
            sage: D.has_edge(0,1,None)
            True

        """
        if not ( self.has_vertex(u) and self.has_vertex(v) ):
            return False
        cdef int u_int = self.get_vertex(u)
        cdef int v_int = self.get_vertex(v)
        return self._cg.has_arc(u_int, v_int)

    def iterator_edges(self, object vertices, bint labels):
        """
        Iterate over the edges incident to a sequence of vertices. Edges are
        assumed to be undirected.

        INPUT:
         - ``vertices`` - a list of vertex labels
         - ``labels`` - boolean, whether to return labels as well

        EXAMPLE::

            sage: G = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: G.add_edge(1,2,None,False)
            sage: list(G.iterator_edges(range(9), False))
            [(1, 2)]
            sage: list(G.iterator_edges(range(9), True))
            [(1, 2, None)]

        """
        cdef object v
        vertices = [self.get_vertex(v) for v in vertices if self.has_vertex(v)]
        cdef int u_int, v_int
        if labels:
            return iter([tuple(sorted(
            (self.vertex_label(v_int),
             self.vertex_label(u_int)
            )))+(None,)
                for v_int in vertices
                    for u_int in self._cg.out_neighbors(v_int)
                        if u_int >= v_int or u_int not in vertices])
        else:
            return iter([tuple(sorted(
            (self.vertex_label(v_int),
             self.vertex_label(u_int)
            )))
                for v_int in vertices
                    for u_int in self._cg.out_neighbors(v_int)
                        if u_int >= v_int or u_int not in vertices])

    def iterator_in_edges(self, object vertices, bint labels):
        """
        Iterate over the incoming edges incident to a sequence of vertices.

        INPUT:
         - ``vertices`` - a list of vertex labels
         - ``labels`` - boolean, whether to return labels as well

        EXAMPLE::

            sage: G = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: G.add_edge(1,2,None,True)
            sage: list(G.iterator_in_edges([1], False))
            []
            sage: list(G.iterator_in_edges([2], False))
            [(1, 2)]
            sage: list(G.iterator_in_edges([2], True))
            [(1, 2, None)]

        """
        cdef object v
        vertices = [self.get_vertex(v) for v in vertices if self.has_vertex(v)]
        cdef int u_int, v_int
        if labels:
            return iter([
            (self.vertex_label(u_int),
             self.vertex_label(v_int),
             None)
                for v_int in vertices
                    for u_int in self._cg.in_neighbors(v_int)])
        else:
            return iter([
            (self.vertex_label(u_int),
             self.vertex_label(v_int))
                for v_int in vertices
                    for u_int in self._cg.in_neighbors(v_int)])

    def iterator_out_edges(self, object vertices, bint labels):
        """
        Iterate over the outbound edges incident to a sequence of vertices.

        INPUT:
         - ``vertices`` - a list of vertex labels
         - ``labels`` - boolean, whether to return labels as well

        EXAMPLE::

            sage: G = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: G.add_edge(1,2,None,True)
            sage: list(G.iterator_out_edges([2], False))
            []
            sage: list(G.iterator_out_edges([1], False))
            [(1, 2)]
            sage: list(G.iterator_out_edges([1], True))
            [(1, 2, None)]

        """
        cdef object u, v
        vertices = [self.get_vertex(v) for v in vertices if self.has_vertex(v)]
        cdef int u_int, v_int
        if labels:
            return iter([
            (self.vertex_label(v_int),
             self.vertex_label(u_int),
             None)
                for v_int in vertices
                    for u_int in self._cg.out_neighbors(v_int)])
        else:
            return iter([
            (self.vertex_label(v_int),
             self.vertex_label(u_int))
                for v_int in vertices
                    for u_int in self._cg.out_neighbors(v_int)])

    def multiple_edges(self, new):
        """
        Get/set whether or not ``self`` allows multiple edges.

        INPUT:

         - ``new`` - boolean (to set) or ``None`` (to get)

        EXAMPLES::

            sage: import sage.graphs.base.dense_graph
            sage: G = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: G.multiple_edges(True)
            Traceback (most recent call last):
            ...
            NotImplementedError: Dense graphs do not support multiple edges.
            sage: G.multiple_edges(None)
            False

        """
        if new is None:
            return False
        if new:
            raise NotImplementedError("Dense graphs do not support multiple edges.")

    def set_edge_label(self, object u, object v, object l, bint directed):
        """
        Label the edge ``(u,v)`` by ``l``.

        INPUT:

         - ``u,v`` - the vertices of the edge
         - ``l`` - the edge label
         - ``directed`` - if False, also set ``(v,u)`` with label ``l``

        EXAMPLE::

            sage: import sage.graphs.base.dense_graph
            sage: G = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: G.set_edge_label(1,2,'a',True)
            Traceback (most recent call last):
            ...
            NotImplementedError: Dense graphs do not support edge labels.
        """
        raise NotImplementedError("Dense graphs do not support edge labels.")
