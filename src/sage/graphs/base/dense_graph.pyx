# sage_setup: distribution = sagemath-graphs
r"""
Fast dense graphs

For an overview of graph data structures in sage, see
:mod:`~sage.graphs.base.overview`.

Usage Introduction
------------------

::

    sage: from sage.graphs.base.dense_graph import DenseGraph

Dense graphs are initialized as follows::

    sage: D = DenseGraph(nverts=10, extra_vertices=10)

This example initializes a dense graph with room for twenty vertices, the first
ten of which are in the graph. In general, the first ``nverts`` are "active."
For example, see that 9 is already in the graph::

    sage: D.verts()
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    sage: D.add_vertex(9)
    9
    sage: D.verts()
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

But 10 is not, until we add it::

    sage: D.add_vertex(10)
    10
    sage: D.verts()
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

You can begin working right away as follows::

    sage: D.add_arc(0, 1)
    sage: D.add_arc(1, 2)
    sage: D.add_arc(1, 0)
    sage: D.has_arc(7, 3)
    False
    sage: D.has_arc(0, 1)
    True
    sage: D.in_neighbors(1)
    [0]
    sage: D.out_neighbors(1)
    [0, 2]
    sage: D.del_all_arcs(0, 1)
    sage: D.has_arc(0, 1)
    False
    sage: D.has_arc(1, 2)
    True
    sage: D.del_vertex(7)
    sage: D.has_arc(7, 3)
    False

Dense graphs do not support multiple or labeled edges.

::

    sage: T = DenseGraph(nverts=3, extra_vertices=2)
    sage: T.add_arc(0, 1)
    sage: T.add_arc(1, 2)
    sage: T.add_arc(2, 0)
    sage: T.has_arc(0, 1)
    True

::

    sage: for _ in range(10): D.add_arc(5, 4)
    sage: D.has_arc(5,4 )
    True

Dense graphs are by their nature directed. As of this writing, you need to do
operations in pairs to treat the undirected case (or use a backend or a Sage
graph)::

    sage: T.has_arc(1, 0)
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

        cdef binary_matrix_t edges

.. NOTE::

    As the edges are stored as the adjacency matrix of the graph, enumerating
    the edges of the graph has complexity `O(n^2)` and enumerating the neighbors
    of a vertex has complexity `O(n)` (where `n` in the size of the bitset
    active_vertices).
    So, the class ``DenseGraph`` should be used for graphs such that the number
    of edges is close to the square of the number of vertices.

"""

# ****************************************************************************
#       Copyright (C) 2008-9 Robert L. Miller <rlmillster@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.data_structures.bitset_base cimport *

from cysignals.memory cimport sig_calloc, sig_realloc, sig_free
from sage.data_structures.binary_matrix cimport *
from libc.string cimport memcpy


cdef extern from "Python.h":
    int unlikely(int) nogil  # Defined by Cython


cdef class DenseGraph(CGraph):
    """
    Compiled dense graphs.

    ::

        sage: from sage.graphs.base.dense_graph import DenseGraph

    Dense graphs are initialized as follows::

        sage: D = DenseGraph(nverts=10, extra_vertices=10)

    INPUT:

    - ``nverts`` -- nonnegative integer; the number of vertices
    - ``extra_vertices`` -- nonnegative integer (default: 10); how many extra
      vertices to allocate
    - ``verts`` -- list (default: ``None``); optional list of vertices to add
    - ``arcs`` -- list (default: ``None``); optional list of arcs to add
    - ``directed`` -- boolean (default: ``None``); whether the graph is directed

    The first ``nverts`` are created as vertices of the graph, and the next
    ``extra_vertices`` can be freely added without reallocation. See top level
    documentation for more details. The input ``verts`` and ``arcs`` are mainly
    for use in pickling.
    """
    def __cinit__(self, int nverts, int extra_vertices=10, verts=None, arcs=None, directed=True):
        """
        Allocation and initialization happen in one place.

        Memory usage is ``O((nverts + extra_vertices)^2)``.

        EXAMPLES::

            sage: from sage.graphs.base.dense_graph import DenseGraph
            sage: D = DenseGraph(nverts=10, extra_vertices=10)
        """
        if not nverts and not extra_vertices:
            raise RuntimeError('dense graphs must allocate space for vertices')

        self.num_verts = nverts
        self.num_arcs = 0
        cdef int total_verts = nverts + extra_vertices
        self._directed = directed

        binary_matrix_init(self.edges, total_verts, total_verts)
        self.in_degrees = <int *> sig_calloc(total_verts, sizeof(int))
        self.out_degrees = <int *> sig_calloc(total_verts, sizeof(int))

        if not self.in_degrees or not self.out_degrees:
            sig_free(self.in_degrees)
            sig_free(self.out_degrees)
            binary_matrix_free(self.edges)
            raise MemoryError

        bitset_init(self.active_vertices, total_verts)
        bitset_set_first_n(self.active_vertices, self.num_verts)

        if verts is not None:
            self.add_vertices(verts)

        if arcs is not None:
            for u, v in arcs:
                self.add_arc(u, v)

    def __dealloc__(self):
        """
        New and dealloc are both tested at class level.
        """
        binary_matrix_free(self.edges)
        sig_free(self.in_degrees)
        sig_free(self.out_degrees)
        bitset_free(self.active_vertices)

    cpdef realloc(self, int total_verts):
        """
        Reallocate the number of vertices to use, without actually adding any.

        INPUT:

        - ``total`` -- integer; the total size to make the array

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
            RuntimeError: requested vertex is past twice the allocated range: use realloc
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
        cdef int i
        if not total_verts:
            raise RuntimeError('dense graphs must allocate space for vertices')

        cdef bitset_t bits
        cdef int min_verts
        if <size_t>total_verts < self.active_vertices.size:
            min_verts = total_verts
            bitset_init(bits, self.active_vertices.size)
            bitset_set_first_n(bits, total_verts)
            if not bitset_issubset(self.active_vertices, bits):
                bitset_free(bits)
                return -1
            bitset_free(bits)
        else:
            min_verts = self.active_vertices.size

        # Resize of self.edges
        binary_matrix_realloc(self.edges, total_verts, total_verts)

        self.in_degrees = <int *> sig_realloc(self.in_degrees, total_verts * sizeof(int))
        self.out_degrees = <int *> sig_realloc(self.out_degrees, total_verts * sizeof(int))

        for i in range(self.active_vertices.size, total_verts):
            self.in_degrees[i] = 0
            self.out_degrees[i] = 0

        bitset_realloc(self.active_vertices, total_verts)

    ###################################
    # Arc functions
    ###################################

    cdef inline int _add_arc_unsafe(self, int u, int v) except -1:
        """
        .. WARNING::

            This method is for internal use only. Use :meth:`add_arc_label_unsafe` instead.

        Add arc (u, v) with label l in only one direction.
        """
        if not binary_matrix_get(self.edges, u, v):
            self.in_degrees[v] += 1
            self.out_degrees[u] += 1
            self.num_arcs += 1
            binary_matrix_set1(self.edges, u, v)

    cdef int add_arc_label_unsafe(self, int u, int v, int l) except -1:
        """
        Add arc ``(u, v)`` to the graph.

        INPUT:

        - ``u``, ``v`` -- nonnegative integers
        """
        if unlikely(l):
            raise ValueError("cannot add a labeled arc to an unlabeled graph")

        self._add_arc_unsafe(u, v)
        if u != v and not self._directed:
            self._add_arc_unsafe(v, u)

    cdef int has_arc_label_unsafe(self, int u, int v, int l) except -1:
        """
        Check whether arc ``(u, v)`` is in the graph.

        INPUT:

        - ``u``, ``v`` -- nonnegative integers, must be in self

        - ``l`` -- positive integer label, or zero for no label, or ``-1`` for any label

        OUTPUT:

        - 0 -- False
        - 1 -- True
        """
        if unlikely(l > 0):
            raise ValueError("cannot locate labeled arc in unlabeled graph")
        return binary_matrix_get(self.edges, u, v)

    cdef inline int _del_arc_unsafe(self, int u, int v) except -1:
        """
        .. WARNING::

            This method is for internal use only. Use :meth:`add_arc_label_unsafe` instead.

        Remove arc (u, v) with label l in only one direction.
        """
        if binary_matrix_get(self.edges, u, v):
            self.in_degrees[v] -= 1
            self.out_degrees[u] -= 1
            self.num_arcs -= 1
            binary_matrix_set0(self.edges, u, v)

    cdef int del_arc_unsafe(self, int u, int v) except -1:
        """
        Delete the arc from ``u`` to ``v``, if it exists.

        INPUT:

        - ``u``, ``v`` -- nonnegative integers, must be in self
        """
        self._del_arc_unsafe(u, v)
        if u != v and not self._directed:
            self._del_arc_unsafe(v, u)

    cdef inline int del_arc_label_unsafe(self, int u, int v, int l) except -1:
        if unlikely(l):
            raise ValueError("cannot delete labeled arc in unlabeled graph")
        return self.del_arc_unsafe(u, v)

    cdef inline int arc_label_unsafe(self, int u, int v) except -1:
        return 0

    cdef int all_arcs_unsafe(self, int u, int v, int* labels, int size) except -1:
        """
        Give the labels of all arcs (u, v).

        INPUT:

        - ``u``, ``v`` -- integers from `0`, ..., `n-1`, where `n` is the
          number of vertices
        - ``labels`` -- must be a pointer to an (allocated) integer array
        - ``size`` -- the length of the array

        OUTPUT: integer; the number of arcs ``(u, v)``.
        ``-1`` -- indicates that the array has been filled with labels, but
        there were more.
        """
        if self.has_arc_unsafe(u, v):
            if size > 0:
                labels[0] = 0
                return 1
            else:
                return -1
        return 0

    def complement(self):
        r"""
        Replace the graph with its complement.

        .. NOTE::

            Assumes that the graph has no loop.

        EXAMPLES::

            sage: from sage.graphs.base.dense_graph import DenseGraph
            sage: G = DenseGraph(5)
            sage: G.add_arc(0, 1)
            sage: G.has_arc(0, 1)
            True
            sage: G.complement()
            sage: G.has_arc(0, 1)
            False
        """
        cdef int num_arcs_old = self.num_arcs

        cdef size_t i
        i = bitset_next(self.active_vertices, 0)
        while i != -1:
            self.add_arc_unsafe(i, i)
            bitset_xor(self.edges.rows[i], self.edges.rows[i], self.active_vertices)
            self.in_degrees[i] = self.num_verts-self.in_degrees[i]
            self.out_degrees[i] = self.num_verts-self.out_degrees[i]
            i = bitset_next(self.active_vertices, i + 1)

        self.num_arcs = self.num_verts*(self.num_verts - 1) - num_arcs_old

    ###################################
    # Neighbor functions
    ###################################

    cdef inline int next_out_neighbor_unsafe(self, int u, int v, int* l) except -2:
        """
        Return the next out-neighbor of ``u`` that is greater than ``v``.

        If ``v`` is ``-1`` return the first neighbor of ``u``.

        Return ``-1`` in case there does not exist such an out-neighbor.

        Set ``l`` to be the label of the first arc.
        """
        l[0] = 0
        return bitset_next(self.edges.rows[u], v + 1)

    cdef inline int next_in_neighbor_unsafe(self, int v, int u, int* l) except -2:
        """
        Return the next in-neighbor of ``v`` that is greater or equal to ``u``.

        If ``u`` is ``-1`` return the first neighbor of ``v``.

        Return ``-1`` in case there does not exist such a in-neighbor.

        Set ``l`` to be the label of the first arc.
        """
        l[0] = 0
        cdef size_t i
        i = bitset_next(self.active_vertices, u + 1)
        while i != -1:
            if binary_matrix_get(self.edges, i, v):
                return i
            i = bitset_next(self.active_vertices, i + 1)
        return -1


cdef int copy_dense_graph(DenseGraph dest, DenseGraph src) except -1:
    r"""
    Unsafely copy ``dest`` over ``src``.

    .. NOTE::

        ``dest.active_vertices`` and ``src.active_vertices`` must be of same size!
    """
    if unlikely(dest.active_vertices.size != src.active_vertices.size):
        raise ValueError("``dest.active_vertices`` and ``src.active_vertices`` must be of same size")
    if unlikely(dest.edges.n_rows != src.edges.n_rows or dest.edges.n_cols != src.edges.n_cols):
        raise ValueError("the edges are not of same size")
    memcpy(dest.in_degrees,  src.in_degrees,  src.active_vertices.size * sizeof(int))
    memcpy(dest.out_degrees, src.out_degrees, src.active_vertices.size * sizeof(int))
    binary_matrix_copy(dest.edges, src.edges)
    bitset_copy(dest.active_vertices, src.active_vertices)
    dest.num_verts = src.num_verts
    dest.num_arcs = src.num_arcs


##############################
# Further tests. Unit tests for methods, functions, classes defined with cdef.
##############################

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
    while n <= 1:
        randg = DiGraph(GraphGenerators().RandomGNP(randint(low, high), random()))
        n = randg.order()

    cdef DenseGraph g = DenseGraph(n,
                                   verts=randg.vertex_iterator(),
                                   arcs=randg.edge_iterator(labels=False))
    assert g.num_verts == randg.order(), (
        "Graph order mismatch: %s vs. %s" % (g.num_verts, randg.order()))
    assert g.num_arcs == randg.size(), (
        "Graph size mismatch: %s vs. %s" % (g.num_arcs, randg.size()))
    M = randg.adjacency_matrix()
    cdef int *V = <int *>sig_malloc(n * sizeof(int))
    cdef int i = 0
    for v in randg.vertex_iterator():
        V[i] = v
        i += 1
    cdef int *seq = <int *> sig_malloc(n * sizeof(int))
    for i in range(randint(50, 101)):
        u = randint(low, n - 1)
        g.adjacency_sequence_out(n, V, u, seq)
        A = [seq[k] for k in range(n)]
        try:
            assert A == list(M[u])
        except AssertionError:
            sig_free(V)
            sig_free(seq)
            raise AssertionError("graph adjacency mismatch")
    sig_free(seq)
    sig_free(V)


###########################################
# Dense Graph Backend
###########################################

cdef class DenseGraphBackend(CGraphBackend):
    """
    Backend for Sage graphs using DenseGraphs.

    ::

        sage: from sage.graphs.base.dense_graph import DenseGraphBackend

    This class is only intended for use by the Sage Graph and DiGraph class.
    If you are interested in using a DenseGraph, you probably want to do
    something like the following example, which creates a Sage Graph instance
    which wraps a ``DenseGraph`` object::

        sage: G = Graph(30, sparse=False)
        sage: G.add_edges([(0, 1), (0, 3), (4, 5), (9, 23)])
        sage: G.edges(sort=True, labels=False)
        [(0, 1), (0, 3), (4, 5), (9, 23)]

    Note that Sage graphs using the backend are more flexible than DenseGraphs
    themselves. This is because DenseGraphs (by design) do not deal with Python
    objects::

        sage: G.add_vertex((0, 1, 2))
        sage: sorted(list(G),
        ....:        key=lambda x: (isinstance(x, tuple), x))
        [0,
        ...
         29,
         (0, 1, 2)]
        sage: from sage.graphs.base.dense_graph import DenseGraph
        sage: DG = DenseGraph(30)
        sage: DG.add_vertex((0, 1, 2))
        Traceback (most recent call last):
        ...
        TypeError: an integer is required
    """

    def __init__(self, n, directed=True):
        """
        Initialize a dense graph with ``n`` vertices.

        EXAMPLES::

            sage: D = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: D.add_edge(0, 1, None, False)
            sage: list(D.iterator_edges(range(9), True))
            [(0, 1, None)]
        """
        self._cg = DenseGraph(n, directed=directed)
        self._directed = directed
        self.vertex_labels = {}
        self.vertex_ints = {}

    cdef bint _delete_edge_before_adding(self) noexcept:
        """
        Return whether we should delete edges before adding any.

        As dense graphs do not support multiple edges, this is never required.
        """
        return False

    cdef inline int new_edge_label(self, object l) except -1:
        """
        Any label is ignored.
        """
        return 0

    cdef inline int free_edge_label(self, int l_int) except -1:
        """
        Free the label corresponding to ``l_int``.

        As the backend does not support labels this raises an error
        if ``l_int`` is non-trivial.
        """
        if unlikely(l_int):
            raise ValueError("backend does not support labels")

    def add_edges(self, object edges, bint directed, bint remove_loops=False):
        """
        Add edges from a list.

        INPUT:

        - ``edges`` -- an iterable of edges to be added; each edge can either be
          of the form ``(u, v)`` or ``(u, v, l)``

        - ``directed`` -- if ``False``, adds ``(v, u)`` as well as ``(u, v)``

        - ``remove_loops`` -- if ``True``, remove loops

        EXAMPLES::

            sage: D = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: D.add_edges([(0, 1), (2, 3), (4, 5), (5, 6)], False)
            sage: list(D.iterator_edges(range(9), True))
            [(0, 1, None),
             (2, 3, None),
             (4, 5, None),
             (5, 6, None)]
        """
        for e in edges:
            u, v = e[:2]
            if unlikely(remove_loops and u == v):
                continue
            self.add_edge(u, v, None, directed)

    def get_edge_label(self, object u, object v):
        """
        Return the edge label for ``(u, v)``.

        Always ``None``, since dense graphs do not support edge labels.

        INPUT:

        - ``u``, ``v`` -- the vertices of the edge

        EXAMPLES::

            sage: D = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: D.add_edges([(0, 1), (2, 3, 7), (4, 5), (5, 6)], False)
            sage: list(D.iterator_edges(range(9), True))
            [(0, 1, None),
             (2, 3, None),
             (4, 5, None),
             (5, 6, None)]
            sage: D.del_edge(0, 1, None, True)
            sage: list(D.iterator_out_edges(range(9), True))
            [(1, 0, None),
             (2, 3, None),
             (3, 2, None),
             (4, 5, None),
             (5, 4, None),
             (5, 6, None),
             (6, 5, None)]
            sage: D.get_edge_label(2, 3)
            sage: D.get_edge_label(2, 4)
            Traceback (most recent call last):
            ...
            LookupError: (2, 4) is not an edge of the graph
        """
        if not self.has_edge(u, v, None):
            raise LookupError("({0}, {1}) is not an edge of the graph".format(repr(u), repr(v)))
        return None

    def has_edge(self, object u, object v, object l):
        """
        Check whether this graph has edge ``(u, v)``.

        .. NOTE::

            The input ``l`` is for consistency with other backends.

        INPUT:

        - ``u``, ``v`` -- the vertices of the edge

        - ``l`` -- the edge label (ignored)

        EXAMPLES::

            sage: D = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: D.add_edges([(0, 1), (2, 3), (4, 5), (5, 6)], False)
            sage: D.has_edge(0, 1, None)
            True
        """
        cdef int u_int = self.get_vertex_checked(u)
        cdef int v_int = self.get_vertex_checked(v)
        if u_int == -1 or v_int == -1:
            return False
        return self._has_labeled_edge_unsafe(u_int, v_int, None)

    cdef inline bint _has_labeled_edge_unsafe(self, int u_int, int v_int, object l) except -1:
        """
        Return whether ``self`` has an arc specified by indices of the vertices
        and an arc label.
        """
        return 1 == self.cg().has_arc_unsafe(u_int, v_int)

    def multiple_edges(self, new):
        """
        Get/set whether or not ``self`` allows multiple edges.

        INPUT:

        - ``new`` -- boolean (to set) or ``None`` (to get)

        EXAMPLES::

            sage: import sage.graphs.base.dense_graph
            sage: G = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: G.multiple_edges(True)
            Traceback (most recent call last):
            ...
            NotImplementedError: dense graphs do not support multiple edges
            sage: G.multiple_edges(None)
            False
        """
        if new is None:
            return False
        if new:
            raise NotImplementedError("dense graphs do not support multiple edges")

    def set_edge_label(self, object u, object v, object l, bint directed):
        """
        Label the edge ``(u, v)`` by ``l``.

        INPUT:

        - ``u``, ``v`` -- the vertices of the edge

        - ``l`` -- the edge label

        - ``directed`` -- if ``False``, also set ``(v, u)`` with label ``l``

        EXAMPLES::

            sage: import sage.graphs.base.dense_graph
            sage: G = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: G.set_edge_label(1, 2, 'a', True)
            Traceback (most recent call last):
            ...
            NotImplementedError: dense graphs do not support edge labels
        """
        raise NotImplementedError("dense graphs do not support edge labels")
