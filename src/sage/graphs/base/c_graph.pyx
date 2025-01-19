# distutils: language = c++
r"""
Fast compiled graphs

This is a Cython implementation of the base class for sparse and dense graphs
in Sage. It is not intended for use on its own. Specific graph types should
extend this base class and implement missing functionalities. Whenever
possible, specific methods should also be overridden with implementations that
suit the graph type under consideration.

For an overview of graph data structures in sage, see
:mod:`~sage.graphs.base.overview`.

Data structure
--------------

The class ``CGraph`` maintains the following variables:

- ``cdef int num_verts``
- ``cdef int num_arcs``
- ``cdef int *in_degrees``
- ``cdef int *out_degrees``
- ``cdef bitset_t active_vertices``

The bitset ``active_vertices`` is a list of all available vertices for use, but
only the ones which are set are considered to actually be in the graph. The
variables ``num_verts`` and ``num_arcs`` are self-explanatory. Note that
``num_verts`` is the number of bits set in ``active_vertices``, not the full
length of the bitset. The arrays ``in_degrees`` and ``out_degrees`` are of the
same length as the bitset.

For more information about active vertices, see the documentation for the
method :meth:`realloc <sage.graphs.base.c_graph.CGraph.realloc>`.
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

from cysignals.memory cimport check_allocarray, sig_free
from libcpp.pair cimport pair
from libcpp.queue cimport queue
from libcpp.stack cimport stack

from sage.arith.long cimport pyobject_to_long
from sage.data_structures.bitset cimport FrozenBitset
from sage.data_structures.bitset_base cimport *
from sage.data_structures.pairing_heap cimport PairingHeap
from sage.rings.integer cimport smallInteger

from sage.rings.integer_ring import ZZ


cdef extern from "Python.h":
    int unlikely(int) nogil  # Defined by Cython


cdef class CGraph:
    """
    Compiled sparse and dense graphs.
    """

    ###################################
    # Vertex Functions
    ###################################

    cpdef bint has_vertex(self, int n) except -1:
        """
        Determine whether the vertex `n` is in ``self``.

        This method is different from :meth:`check_vertex`. The current method
        returns a boolean to signify whether or not `n` is a vertex of this
        graph. On the other hand, :meth:`check_vertex` raises an error if
        `n` is not a vertex of this graph.

        INPUT:

        - ``n`` -- nonnegative integer representing a vertex

        OUTPUT: ``True`` if `n` is a vertex of this graph; ``False`` otherwise

        .. SEEALSO::

            - :meth:`check_vertex`
              -- raise an error if this graph does not contain a specific
              vertex.

        EXAMPLES:

        Upon initialization, a
        :class:`SparseGraph <sage.graphs.base.sparse_graph.SparseGraph>`
        or
        :class:`DenseGraph <sage.graphs.base.dense_graph.DenseGraph>`
        has the first ``nverts`` vertices::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: S = SparseGraph(nverts=10, expected_degree=3, extra_vertices=10)
            sage: S.has_vertex(6)
            True
            sage: S.has_vertex(12)
            False
            sage: S.has_vertex(24)
            False
            sage: S.has_vertex(-19)
            False

        ::

            sage: from sage.graphs.base.dense_graph import DenseGraph
            sage: D = DenseGraph(nverts=10, extra_vertices=10)
            sage: D.has_vertex(6)
            True
            sage: D.has_vertex(12)
            False
            sage: D.has_vertex(24)
            False
            sage: D.has_vertex(-19)
            False
        """
        return (n >= 0 and
                <mp_bitcnt_t>n < self.active_vertices.size and
                bitset_in(self.active_vertices, n))

    cpdef check_vertex(self, int n):
        """
        Check that `n` is a vertex of ``self``.

        This method is different from :meth:`has_vertex`. The current method
        raises an error if `n` is not a vertex of this graph. On the other
        hand, :meth:`has_vertex` returns a boolean to signify whether or not
        `n` is a vertex of this graph.

        INPUT:

        - ``n`` -- nonnegative integer representing a vertex

        OUTPUT: raise an error if `n` is not a vertex of this graph

        .. SEEALSO::

            - :meth:`has_vertex`
              -- determine whether this graph has a specific vertex

        EXAMPLES::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: S = SparseGraph(nverts=10, expected_degree=3, extra_vertices=10)
            sage: S.check_vertex(4)
            sage: S.check_vertex(12)
            Traceback (most recent call last):
            ...
            LookupError: vertex (12) is not a vertex of the graph
            sage: S.check_vertex(24)
            Traceback (most recent call last):
            ...
            LookupError: vertex (24) is not a vertex of the graph
            sage: S.check_vertex(-19)
            Traceback (most recent call last):
            ...
            LookupError: vertex (-19) is not a vertex of the graph

        ::

            sage: from sage.graphs.base.dense_graph import DenseGraph
            sage: D = DenseGraph(nverts=10, extra_vertices=10)
            sage: D.check_vertex(4)
            sage: D.check_vertex(12)
            Traceback (most recent call last):
            ...
            LookupError: vertex (12) is not a vertex of the graph
            sage: D.check_vertex(24)
            Traceback (most recent call last):
            ...
            LookupError: vertex (24) is not a vertex of the graph
            sage: D.check_vertex(-19)
            Traceback (most recent call last):
            ...
            LookupError: vertex (-19) is not a vertex of the graph
        """
        if not self.has_vertex(n):
            raise LookupError("vertex ({0}) is not a vertex of the graph".format(n))

    cdef int add_vertex_unsafe(self, int k) except -1:
        """
        Add the vertex `k` to the graph.

        INPUT:

        - ``k`` -- nonnegative integer or ``-1``; for `k >= 0`, add the vertex
          ``k`` to this graph if the vertex is not already in the graph.  If `k
          = -1`, this function will find the first available vertex that is not
          in ``self`` and add that vertex to this graph.

        OUTPUT:

        - ``-1`` -- indicates that no vertex was added because the current
          allocation is already full or the vertex is out of range

        - nonnegative integer -- this vertex is now guaranteed to be in the
          graph

        .. WARNING::

            This method is potentially unsafe. You should instead use
            :meth:`add_vertex`.
        """
        if k == -1:
            k = bitset_first_in_complement(self.active_vertices)
        elif self.active_vertices.size <= <mp_bitcnt_t>k:
            k = -1
        if k != -1:
            if not bitset_in(self.active_vertices, k):
                self.num_verts += 1
            bitset_add(self.active_vertices, k)
        return k

    def add_vertex(self, int k=-1):
        """
        Add vertex ``k`` to the graph.

        INPUT:

        - ``k`` -- nonnegative integer or ``-1`` (default: ``-1``); if `k = -1`,
          a new vertex is added and the integer used is returned.  That is, for
          `k = -1`, this function will find the first available vertex that is
          not in ``self`` and add that vertex to this graph.

        OUTPUT:

        - ``-1`` -- indicates that no vertex was added because the current
          allocation is already full or the vertex is out of range

        - nonnegative integer -- this vertex is now guaranteed to be in the
          graph

        .. SEEALSO::

            - ``add_vertex_unsafe`` -- add a vertex to a graph. This method is
              potentially unsafe. You should instead use :meth:`add_vertex`

            - ``add_vertices`` -- add a bunch of vertices to a graph

        EXAMPLES:

        Adding vertices to a sparse graph::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(3, extra_vertices=3)
            sage: G.add_vertex(3)
            3
            sage: G.add_arc(2, 5)
            Traceback (most recent call last):
            ...
            LookupError: vertex (5) is not a vertex of the graph
            sage: G.add_arc(1, 3)
            sage: G.has_arc(1, 3)
            True
            sage: G.has_arc(2, 3)
            False

        Adding vertices to a dense graph::

            sage: from sage.graphs.base.dense_graph import DenseGraph
            sage: G = DenseGraph(3, extra_vertices=3)
            sage: G.add_vertex(3)
            3
            sage: G.add_arc(2,5)
            Traceback (most recent call last):
            ...
            LookupError: vertex (5) is not a vertex of the graph
            sage: G.add_arc(1, 3)
            sage: G.has_arc(1, 3)
            True
            sage: G.has_arc(2, 3)
            False

        Repeatedly adding a vertex using `k = -1` will allocate more memory
        as required::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(3, extra_vertices=0)
            sage: G.verts()
            [0, 1, 2]
            sage: for i in range(10):
            ....:     _ = G.add_vertex(-1);
            ...
            sage: G.verts()
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]

        ::

            sage: from sage.graphs.base.dense_graph import DenseGraph
            sage: G = DenseGraph(3, extra_vertices=0)
            sage: G.verts()
            [0, 1, 2]
            sage: for i in range(12):
            ....:     _ = G.add_vertex(-1);
            ...
            sage: G.verts()
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]

        TESTS::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(3, extra_vertices=0)
            sage: G.add_vertex(6)
            Traceback (most recent call last):
            ...
            RuntimeError: requested vertex is past twice the allocated range: use realloc

        ::

            sage: from sage.graphs.base.dense_graph import DenseGraph
            sage: G = DenseGraph(3, extra_vertices=0)
            sage: G.add_vertex(6)
            Traceback (most recent call last):
            ...
            RuntimeError: requested vertex is past twice the allocated range: use realloc
        """
        if k >= (2 * <int>self.active_vertices.size):
            raise RuntimeError(
                "requested vertex is past twice the allocated range: "
                "use realloc")
        if (k >= <int>self.active_vertices.size or
                (k == -1 and self.active_vertices.size == <mp_bitcnt_t>self.num_verts)):
            self.realloc(2 * self.active_vertices.size)
        return self.add_vertex_unsafe(k)

    cpdef add_vertices(self, verts):
        """
        Add vertices from the iterable ``verts``.

        INPUT:

        - ``verts`` -- an iterable of vertices; value -1 has a special meaning
          -- for each such value an unused vertex name is found, used to create
          a new vertex and returned.

        OUTPUT:

        List of generated labels if there is any -1 in ``verts``. ``None``
        otherwise.

        .. SEEALSO::

            - :meth:`add_vertex` -- add a vertex to a graph

        EXAMPLES:

        Adding vertices for sparse graphs::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: S = SparseGraph(nverts=4, extra_vertices=4)
            sage: S.verts()
            [0, 1, 2, 3]
            sage: S.add_vertices([3, -1, 4, 9])
            [5]
            sage: S.verts()
            [0, 1, 2, 3, 4, 5, 9]
            sage: S.realloc(20)
            sage: S.verts()
            [0, 1, 2, 3, 4, 5, 9]

        Adding vertices for dense graphs::

            sage: from sage.graphs.base.dense_graph import DenseGraph
            sage: D = DenseGraph(nverts=4, extra_vertices=4)
            sage: D.verts()
            [0, 1, 2, 3]
            sage: D.add_vertices([3, -1, 4, 9])
            [5]
            sage: D.verts()
            [0, 1, 2, 3, 4, 5, 9]
            sage: D.realloc(20)
            sage: D.verts()
            [0, 1, 2, 3, 4, 5, 9]
        """
        cdef int v
        cdef int nones = 0
        for v in verts:
            if v > -1:
                self.add_vertex(v)
            else:
                nones += 1

        new_names = []
        while nones:
            new_names.append(self.add_vertex())
            nones -= 1

        return new_names if new_names else None

    cdef int del_vertex_unsafe(self, int v) except -1:
        """
        Delete the vertex ``v``, along with all edges incident to it.

        INPUT:

        - ``v`` -- nonnegative integer representing a vertex

        .. WARNING::

            This method is potentially unsafe. Use :meth:`del_vertex` instead.
        """
        cdef int size = 0
        cdef int num_nbrs
        cdef int i
        cdef int *neighbors
        if self.in_degrees[v] > size:
            size = self.in_degrees[v]
        if self.out_degrees[v] > size:
            size = self.out_degrees[v]
        if size > 0:
            neighbors = <int *> sig_malloc(size * sizeof(int))
            if not neighbors:
                raise RuntimeError("failure allocating memory")
            # delete each arc incident with v
            num_nbrs = self.in_neighbors_unsafe(v, neighbors, size)
            for i in range(num_nbrs):
                self.del_arc_unsafe(neighbors[i], v)
            num_nbrs = self.out_neighbors_unsafe(v, neighbors, size)
            for i in range(num_nbrs):
                self.del_arc_unsafe(v, neighbors[i])
            sig_free(neighbors)

        self.num_verts -= 1
        bitset_remove(self.active_vertices, v)

    cpdef del_vertex(self, int v):
        """
        Delete the vertex `v`, along with all edges incident to it.

        If `v` is not in ``self``, fails silently.

        INPUT:

        - ``v`` -- nonnegative integer representing a vertex

        .. SEEALSO::

            - ``del_vertex_unsafe`` -- delete a vertex from a graph. This method
              is potentially unsafe. Use :meth:`del_vertex` instead

        EXAMPLES:

        Deleting vertices of sparse graphs::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(3)
            sage: G.add_arc(0, 1)
            sage: G.add_arc(0, 2)
            sage: G.add_arc(1, 2)
            sage: G.add_arc(2, 0)
            sage: G.del_vertex(2)
            sage: for i in range(2):
            ....:     for j in range(2):
            ....:         if G.has_arc(i, j):
            ....:             print("{} {}".format(i,j))
            0 1
            sage: G = SparseGraph(3)
            sage: G.add_arc(0, 1)
            sage: G.add_arc(0, 2)
            sage: G.add_arc(1, 2)
            sage: G.add_arc(2, 0)
            sage: G.del_vertex(1)
            sage: for i in range(3):
            ....:     for j in range(3):
            ....:         if G.has_arc(i, j):
            ....:             print("{} {}".format(i,j))
            0 2
            2 0

        Deleting vertices of dense graphs::

            sage: from sage.graphs.base.dense_graph import DenseGraph
            sage: G = DenseGraph(4)
            sage: G.add_arc(0, 1); G.add_arc(0, 2)
            sage: G.add_arc(3, 1); G.add_arc(3, 2)
            sage: G.add_arc(1, 2)
            sage: G.verts()
            [0, 1, 2, 3]
            sage: G.del_vertex(3); G.verts()
            [0, 1, 2]
            sage: for i in range(3):
            ....:     for j in range(3):
            ....:         if G.has_arc(i, j):
            ....:             print("{} {}".format(i,j))
            0 1
            0 2
            1 2

        If the vertex to be deleted is not in this graph, then fail silently::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(3)
            sage: G.verts()
            [0, 1, 2]
            sage: G.has_vertex(3)
            False
            sage: G.del_vertex(3)
            sage: G.verts()
            [0, 1, 2]

        ::

            sage: from sage.graphs.base.dense_graph import DenseGraph
            sage: G = DenseGraph(5)
            sage: G.verts()
            [0, 1, 2, 3, 4]
            sage: G.has_vertex(6)
            False
            sage: G.del_vertex(6)
            sage: G.verts()
            [0, 1, 2, 3, 4]
        """
        if self.has_vertex(v):
            self.del_vertex_unsafe(v)

    cpdef int current_allocation(self) noexcept:
        r"""
        Report the number of vertices allocated.

        OUTPUT:

        - The number of vertices allocated. This number is usually different
          from the order of a graph. We may have allocated enough memory for a
          graph to hold `n > 0` vertices, but the order (actual number of
          vertices) of the graph could be less than `n`.

        EXAMPLES::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: S = SparseGraph(nverts=4, extra_vertices=4)
            sage: S.current_allocation()
            8
            sage: S.add_vertex(6)
            6
            sage: S.current_allocation()
            8
            sage: S.add_vertex(10)
            10
            sage: S.current_allocation()
            16
            sage: S.add_vertex(40)
            Traceback (most recent call last):
            ...
            RuntimeError: requested vertex is past twice the allocated range: use realloc
            sage: S.realloc(50)
            sage: S.add_vertex(40)
            40
            sage: S.current_allocation()
            50
            sage: S.realloc(30)
            -1
            sage: S.current_allocation()
            50
            sage: S.del_vertex(40)
            sage: S.realloc(30)
            sage: S.current_allocation()
            30

        The actual number of vertices in a graph might be less than the number
        of vertices allocated for the graph::

            sage: from sage.graphs.base.dense_graph import DenseGraph
            sage: G = DenseGraph(nverts=3, extra_vertices=2)
            sage: order = len(G.verts())
            sage: order
            3
            sage: G.current_allocation()
            5
            sage: order < G.current_allocation()
            True
        """
        return self.active_vertices.size

    cpdef list verts(self):
        """
        Return a list of the vertices in ``self``.

        EXAMPLES::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: S = SparseGraph(nverts=4, extra_vertices=4)
            sage: S.verts()
            [0, 1, 2, 3]
            sage: S.add_vertices([3,5,7,9])
            sage: S.verts()
            [0, 1, 2, 3, 5, 7, 9]
            sage: S.realloc(20)
            sage: S.verts()
            [0, 1, 2, 3, 5, 7, 9]

        ::

            sage: from sage.graphs.base.dense_graph import DenseGraph
            sage: G = DenseGraph(3, extra_vertices=2)
            sage: G.verts()
            [0, 1, 2]
            sage: G.del_vertex(0)
            sage: G.verts()
            [1, 2]
        """
        return bitset_list(self.active_vertices)

    cpdef realloc(self, int total):
        """
        Reallocate the number of vertices to use, without actually adding any.

        INPUT:

        - ``total`` -- integer; the total size to make the array of vertices

        OUTPUT:

        - Raise a :exc:`NotImplementedError`. This method is not implemented
          in this base class. A child class should provide a suitable
          implementation.

        .. SEEALSO::

            - :meth:`realloc <sage.graphs.base.sparse_graph.SparseGraph.realloc>`
              -- a ``realloc`` implementation for sparse graphs.

            - :meth:`realloc <sage.graphs.base.dense_graph.DenseGraph.realloc>`
              -- a ``realloc`` implementation for dense graphs.

        EXAMPLES:

        First, note that :meth:`realloc` is implemented for
        :class:`SparseGraph <sage.graphs.base.sparse_graph.SparseGraph>`
        and
        :class:`DenseGraph <sage.graphs.base.dense_graph.DenseGraph>`
        differently, and is not implemented at the
        :class:`CGraph` level::

            sage: from sage.graphs.base.c_graph import CGraph
            sage: G = CGraph()
            sage: G.realloc(20)
            Traceback (most recent call last):
            ...
            NotImplementedError

        The ``realloc`` implementation for sparse graphs::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: S = SparseGraph(nverts=4, extra_vertices=4)
            sage: S.current_allocation()
            8
            sage: S.add_vertex(6)
            6
            sage: S.current_allocation()
            8
            sage: S.add_vertex(10)
            10
            sage: S.current_allocation()
            16
            sage: S.add_vertex(40)
            Traceback (most recent call last):
            ...
            RuntimeError: requested vertex is past twice the allocated range: use realloc
            sage: S.realloc(50)
            sage: S.add_vertex(40)
            40
            sage: S.current_allocation()
            50
            sage: S.realloc(30)
            -1
            sage: S.current_allocation()
            50
            sage: S.del_vertex(40)
            sage: S.realloc(30)
            sage: S.current_allocation()
            30

        The ``realloc`` implementation for dense graphs::

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
        raise NotImplementedError()

    ###################################
    # Edge Functions
    ###################################

    cdef int del_arc_unsafe(self, int u, int v) except -1:
        raise NotImplementedError()

    cpdef add_arc(self, int u, int v):
        """
        Add arc ``(u, v)`` to the graph.

        INPUT:

        - ``u``, ``v`` -- nonnegative integers; must be in self

        EXAMPLES:

        On the :class:`CGraph` level, this always produces an error, as there are no vertices::

            sage: from sage.graphs.base.c_graph import CGraph
            sage: G = CGraph()
            sage: G.add_arc(0, 1)
            Traceback (most recent call last):
            ...
            LookupError: vertex (0) is not a vertex of the graph

        It works, once there are vertices and :meth:`add_arc_unsafe` is implemented::

            sage: from sage.graphs.base.dense_graph import DenseGraph
            sage: G = DenseGraph(5)
            sage: G.add_arc(0, 1)
            sage: G.add_arc(4, 7)
            Traceback (most recent call last):
            ...
            LookupError: vertex (7) is not a vertex of the graph
            sage: G.has_arc(1, 0)
            False
            sage: G.has_arc(0, 1)
            True

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(5)
            sage: G.add_arc(0,1)
            sage: G.add_arc(4,7)
            Traceback (most recent call last):
            ...
            LookupError: vertex (7) is not a vertex of the graph
            sage: G.has_arc(1,0)
            False
            sage: G.has_arc(0,1)
            True
        """
        self.check_vertex(u)
        self.check_vertex(v)
        self.add_arc_unsafe(u, v)

    cpdef bint has_arc(self, int u, int v) except -1:
        """
        Check if the arc ``(u, v)`` is in this graph.

        INPUT:

        - ``u`` -- integer; the tail of an arc

        - ``v`` -- integer; the head of an arc

        OUTPUT:

        - Print a ``Not Implemented!`` message. This method is not implemented
          at the :class:`CGraph` level. A child class should provide a suitable
          implementation.

        EXAMPLES:

        On the :class:`CGraph` this always returns ``False``::

            sage: from sage.graphs.base.c_graph import CGraph
            sage: G = CGraph()
            sage: G.has_arc(0, 1)
            False

        It works once :class:`has_arc_unsafe` is implemented::

            sage: from sage.graphs.base.dense_graph import DenseGraph
            sage: G = DenseGraph(5)
            sage: G.add_arc(0, 1)
            sage: G.has_arc(1, 0)
            False
            sage: G.has_arc(0, 1)
            True

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(5)
            sage: G.add_arc_label(0,1)
            sage: G.has_arc(1,0)
            False
            sage: G.has_arc(0,1)
            True
        """
        if u < 0 or u >= <int>self.active_vertices.size or not bitset_in(self.active_vertices, u):
            return False
        if v < 0 or v >= <int>self.active_vertices.size or not bitset_in(self.active_vertices, v):
            return False
        return self.has_arc_unsafe(u, v) == 1

    cpdef del_all_arcs(self, int u, int v):
        """
        Delete all arcs from ``u`` to ``v``.

        INPUT:

        - ``u`` -- integer; the tail of an arc

        - ``v`` -- integer; the head of an arc

        EXAMPLES:

        On the :class:`CGraph` level, this always produces an error, as there are no vertices::

            sage: from sage.graphs.base.c_graph import CGraph
            sage: G = CGraph()
            sage: G.del_all_arcs(0,1)
            Traceback (most recent call last):
            ...
            LookupError: vertex (0) is not a vertex of the graph

        It works, once there are vertices and :meth:`del_arc_unsafe` is implemented::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(5)
            sage: G.add_arc_label(0,1,0)
            sage: G.add_arc_label(0,1,1)
            sage: G.add_arc_label(0,1,2)
            sage: G.add_arc_label(0,1,3)
            sage: G.del_all_arcs(0,1)
            sage: G.has_arc(0,1)
            False
            sage: G.arc_label(0,1)
            0
            sage: G.del_all_arcs(0,1)

            sage: from sage.graphs.base.dense_graph import DenseGraph
            sage: G = DenseGraph(5)
            sage: G.add_arc(0, 1)
            sage: G.has_arc(0, 1)
            True
            sage: G.del_all_arcs(0, 1)
            sage: G.has_arc(0, 1)
            False
        """
        self.check_vertex(u)
        self.check_vertex(v)
        self.del_arc_unsafe(u, v)

    ###################################
    # Labeled Edge Functions
    ###################################

    cdef int add_arc_label_unsafe(self, int u, int v, int l) except -1:
        raise NotImplementedError()

    cdef int has_arc_label_unsafe(self, int u, int v, int l) except -1:
        raise NotImplementedError()

    cdef int del_arc_label_unsafe(self, int u, int v, int l) except -1:
        raise NotImplementedError()

    cdef int arc_label_unsafe(self, int u, int v) except -1:
        raise NotImplementedError()

    cdef int all_arcs_unsafe(self, int u, int v, int* arc_labels, int size) except -1:
        raise NotImplementedError()

    cpdef int arc_label(self, int u, int v) noexcept:
        """
        Retrieves the first label found associated with ``(u, v)``.

        INPUT:

        - ``u``, ``v`` -- nonnegative integers; must be in self

        OUTPUT: one of

        - positive integer -- indicates that there is a label on ``(u, v)``

        - ``0`` -- either the arc ``(u, v)`` is unlabeled, or there is no arc at all

        EXAMPLES::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(5)
            sage: G.add_arc_label(3,4,7)
            sage: G.arc_label(3,4)
            7

        To this function, an unlabeled arc is indistinguishable from a non-arc::

            sage: G.add_arc_label(1,0)
            sage: G.arc_label(1,0)
            0
            sage: G.arc_label(1,1)
            0

        This function only returns the *first* label it finds from ``u`` to ``v``::

            sage: G.add_arc_label(1,2,1)
            sage: G.add_arc_label(1,2,2)
            sage: G.arc_label(1,2)
            2
        """
        self.check_vertex(u)
        self.check_vertex(v)
        return self.arc_label_unsafe(u, v)

    cpdef list all_arcs(self, int u, int v):
        """
        Give the labels of all arcs ``(u, v)``. An unlabeled arc is interpreted as
        having label 0.

        EXAMPLES::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(5)
            sage: G.add_arc_label(1,2,1)
            sage: G.add_arc_label(1,2,2)
            sage: G.add_arc_label(1,2,2)
            sage: G.add_arc_label(1,2,2)
            sage: G.add_arc_label(1,2,3)
            sage: G.add_arc_label(1,2,3)
            sage: G.add_arc_label(1,2,4)
            sage: G.all_arcs(1,2)
            [4, 3, 3, 2, 2, 2, 1]
        """
        cdef int size, num_arcs, i
        cdef int *arc_labels
        cdef list output
        self.check_vertex(u)
        self.check_vertex(v)
        if unlikely(self.in_degrees is NULL or self.out_degrees is NULL):
            raise ValueError("`self.in_degree` or `self.out_degree` not allocated")
        if self.in_degrees[v] < self.out_degrees[u]:
            size = self.in_degrees[v]
        else:
            size = self.out_degrees[u]
        arc_labels = <int *>check_allocarray(size, sizeof(int))
        num_arcs = self.all_arcs_unsafe(u, v, arc_labels, size)
        if num_arcs == -1:
            sig_free(arc_labels)
            raise RuntimeError("There was an error: there seem to be more arcs than self.in_degrees or self.out_degrees indicate.")
        output = [arc_labels[i] for i in range(num_arcs)]
        sig_free(arc_labels)
        return output

    cpdef del_arc_label(self, int u, int v, int l):
        """
        Delete an arc ``(u, v)`` with label ``l``.

        INPUT:

        - ``u``, ``v`` -- nonnegative integers; must be in self

        - ``l`` -- positive integer label, or zero for no label

        EXAMPLES::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(5)
            sage: G.add_arc_label(0,1,0)
            sage: G.add_arc_label(0,1,1)
            sage: G.add_arc_label(0,1,2)
            sage: G.add_arc_label(0,1,2)
            sage: G.add_arc_label(0,1,3)
            sage: G.del_arc_label(0,1,2)
            sage: G.all_arcs(0,1)
            [0, 3, 2, 1]
            sage: G.del_arc_label(0,1,0)
            sage: G.all_arcs(0,1)
            [3, 2, 1]
        """
        self.check_vertex(u)
        self.check_vertex(v)
        if l < 0:
            raise ValueError("Label ({0}) must be a nonnegative integer.".format(l))
        self.del_arc_label_unsafe(u, v, l)

    cpdef bint has_arc_label(self, int u, int v, int l) noexcept:
        """
        Indicates whether there is an arc ``(u, v)`` with label ``l``.

        INPUT:

        - ``u``, ``v`` -- nonnegative integers; must be in self

        - ``l`` -- positive integer label, or zero for no label

        EXAMPLES::

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(5)
            sage: G.add_arc_label(0,1,0)
            sage: G.add_arc_label(0,1,1)
            sage: G.add_arc_label(0,1,2)
            sage: G.add_arc_label(0,1,2)
            sage: G.has_arc_label(0,1,1)
            True
            sage: G.has_arc_label(0,1,2)
            True
            sage: G.has_arc_label(0,1,3)
            False
        """
        self.check_vertex(u)
        self.check_vertex(v)
        if l < 0:
            raise ValueError("Label ({0}) must be a nonnegative integer.".format(l))
        return self.has_arc_label_unsafe(u, v, l) == 1

    ###################################
    # Neighbor Functions
    ###################################

    cdef int out_neighbors_unsafe(self, int u, int *neighbors, int size) except -2:
        """
        Feed array ``neighbors`` with the out-neighbors of ``u``.

        This function will put at most ``size`` out-neighbors of ``u`` in array
        ``neighbors``. If ``u`` has more than ``size`` out-neighbors, ``size``
        of them are put in array ``neighbors`` and the function returns value
        ``-1``.  Otherwise the function returns the number of out-neighbors that
        have been put in array ``neighbors``.

        INPUT:

        - ``u`` -- nonnegative integer; must be in self

        - ``neighbors`` -- pointer to an (allocated) integer array

        - ``size`` -- the length of the array

        OUTPUT:

        - nonnegative integer -- the out-degree of ``u``

        - ``-1`` -- indicates that the array has been filled with neighbors, but
          there were more
        """
        cdef int num_nbrs = 0
        cdef int l
        cdef int v = self.next_out_neighbor_unsafe(u, -1, &l)
        while v != -1:
            if num_nbrs == size:
                return -1
            neighbors[num_nbrs] = v
            num_nbrs += 1
            v = self.next_out_neighbor_unsafe(u, v, &l)

        return num_nbrs

    cdef int in_neighbors_unsafe(self, int u, int *neighbors, int size) except -2:
        """
        Feed array ``neighbors`` with the in-neighbors of ``v``.

        This function will put at most ``size`` in-neighbors of ``v`` in array
        ``neighbors``. If ``v`` has more than ``size`` in-neighbors, ``size`` of
        them are put in array ``neighbors`` and the function returns value
        ``-1``.  Otherwise the function returns the number of in-neighbors that
        have been put in array ``neighbors``.

        INPUT:

        - ``v`` -- nonnegative integer; must be in self

        - ``neighbors`` -- pointer to an (allocated) integer array

        - ``size`` -- the length of the array

        OUTPUT:

        - nonnegative integer -- the in-degree of ``v``

        - ``-1`` -- indicates that the array has been filled with neighbors, but
          there were more
        """
        cdef int num_nbrs = 0
        cdef int l
        cdef int v = self.next_in_neighbor_unsafe(u, -1, &l)
        while v != -1:
            if num_nbrs == size:
                return -1
            neighbors[num_nbrs] = v
            num_nbrs += 1
            v = self.next_in_neighbor_unsafe(u, v, &l)

        return num_nbrs

    cdef int next_out_neighbor_unsafe(self, int u, int v, int* l) except -2:
        raise NotImplementedError()

    cdef int next_in_neighbor_unsafe(self, int v, int u, int* l) except -2:
        raise NotImplementedError()

    cdef adjacency_sequence_out(self, int n, int *vertices, int v, int* sequence):
        r"""
        Return the adjacency sequence corresponding to a list of vertices and a
        vertex.

        This method fills the array ``sequence``, whose `i`-th element is set to
        `1` iff ``(v,vertices[i])`` is an edge.

        See the function ``_test_adjacency_sequence()`` of ``dense_graph.pyx``
        and ``sparse_graph.pyx`` for unit tests.

        INPUT:

        - ``n`` -- nonnegative integer; the maximum index in ``vertices`` up to
          which we want to consider. If ``n = 0``, we only want to know if ``(v,
          vertices[0])`` is an edge. If ``n = 1``, we want to know whether ``(v,
          vertices[0])`` and ``(v, vertices[1])`` are edges.  Let ``k`` be the
          length of ``vertices``. If ``0 <= n < k``, then we want to know if
          each of ``vertices[0], vertices[1], ..., vertices[n]`` is adjacent to
          ``v``. When ``n = k - 1``, then we consider all elements in the list
          ``vertices``.

        - ``vertices`` -- list of vertices

        - ``v`` -- a vertex

        - ``sequence`` -- ``int *``; the memory segment of length `>= n` that is
          to be filled

        .. SEEALSO::

            - :meth:`adjacency_sequence_in` -- similar method for
            ``(vertices[i],v)`` instead of ``(v,vertices[i])`` (the difference
            only matters for digraphs)
        """
        cdef int i
        for i in range(n):
            sequence[i] = self.has_arc_unsafe(v, vertices[i])

    cdef adjacency_sequence_in(self, int n, int *vertices, int v, int* sequence):
        r"""
        Compute the adjacency sequence corresponding to a list of vertices and a
        vertex.

        This method fills the array ``sequence``, whose `i`-th element is set to
        `1` iff ``(v, vertices[i])`` is an edge.

        See the function ``_test_adjacency_sequence()`` of ``dense_graph.pyx``
        and ``sparse_graph.pyx`` for unit tests.

        INPUT:

        - ``n`` -- nonnegative integer; the maximum index in ``vertices`` up to
          which we want to consider. If ``n = 0``, we only want to know if
          ``(vertices[0],v)`` is an edge. If ``n = 1``, we want to know whether
          ``(vertices[0],v)`` and ``(vertices[1],v)`` are edges.  Let ``k`` be
          the length of ``vertices``. If ``0 <= n < k``, then we want to know if
          ``v`` is adjacent to each of ``vertices[0], vertices[1], ...,
          vertices[n]``. When ``n = k - 1``, then we consider all elements in
          the list ``vertices``.

        - ``vertices`` -- list of vertices

        - ``v`` -- a vertex

        - ``sequence`` -- ``int *``; the memory segment of length `>= n` that is
          to be filled

        .. SEEALSO::

            - :meth:`adjacency_sequence_out` -- similar method for ``(v,
            vertices[i])`` instead of ``(vertices[i], v)`` (the difference only
            matters for digraphs)
        """
        cdef int i
        for i in range(n):
            sequence[i] = self.has_arc_unsafe(vertices[i], v)

    cpdef list out_neighbors(self, int u):
        """
        Return the list of out-neighbors of the vertex ``u``.

        INPUT:

        - ``u`` -- integer representing a vertex of this graph

        EXAMPLES:

        On the :class:`CGraph` level, this always produces an error, as there are no vertices::

            sage: from sage.graphs.base.c_graph import CGraph
            sage: G = CGraph()
            sage: G.out_neighbors(0)
            Traceback (most recent call last):
            ...
            LookupError: vertex (0) is not a vertex of the graph

        It works, once there are vertices and :meth:`out_neighbors_unsafe` is implemented::

            sage: from sage.graphs.base.dense_graph import DenseGraph
            sage: G = DenseGraph(5)
            sage: G.add_arc(0, 1)
            sage: G.add_arc(1, 2)
            sage: G.add_arc(1, 3)
            sage: G.out_neighbors(0)
            [1]
            sage: G.out_neighbors(1)
            [2, 3]

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(5)
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
        if not self.out_degrees[u]:
            return []
        cdef int size = self.out_degrees[u]
        cdef int *neighbors = <int *>check_allocarray(size, sizeof(int))
        if not neighbors:
            raise MemoryError
        num_nbrs = self.out_neighbors_unsafe(u, neighbors, size)
        output = [neighbors[i] for i in range(num_nbrs)]
        sig_free(neighbors)
        return output

    cpdef list in_neighbors(self, int v):
        """
        Return the list of in-neighbors of the vertex ``v``.

        INPUT:

        - ``v`` -- integer representing a vertex of this graph

        OUTPUT:

        - Raise :exc:`NotImplementedError`. This method is not implemented at
          the :class:`CGraph` level. A child class should provide a suitable
          implementation.

        .. NOTE::

            Due to the implementation of SparseGraph, this method is much more
            expensive than out_neighbors_unsafe for SparseGraph's.

        EXAMPLES:

        On the :class:`CGraph` level, this always produces an error, as there are no vertices::

            sage: from sage.graphs.base.c_graph import CGraph
            sage: G = CGraph()
            sage: G.in_neighbors(0)
            Traceback (most recent call last):
            ...
            LookupError: vertex (0) is not a vertex of the graph

        It works, once there are vertices and :meth:`out_neighbors_unsafe` is implemented::

            sage: from sage.graphs.base.dense_graph import DenseGraph
            sage: G = DenseGraph(5)
            sage: G.add_arc(0, 1)
            sage: G.add_arc(3, 1)
            sage: G.add_arc(1, 3)
            sage: G.in_neighbors(1)
            [0, 3]
            sage: G.in_neighbors(3)
            [1]

            sage: from sage.graphs.base.sparse_graph import SparseGraph
            sage: G = SparseGraph(5)
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
        if not self.in_degrees[v]:
            return []
        cdef int size = self.in_degrees[v]
        cdef int *neighbors = <int *> check_allocarray(size, sizeof(int))
        if not neighbors:
            raise MemoryError
        num_nbrs = self.in_neighbors_unsafe(v, neighbors, size)
        output = [neighbors[i] for i in range(num_nbrs)]
        sig_free(neighbors)
        return output


cdef class CGraphBackend(GenericGraphBackend):
    """
    Base class for sparse and dense graph backends.

    ::

        sage: from sage.graphs.base.c_graph import CGraphBackend

    This class is extended by
    :class:`SparseGraphBackend <sage.graphs.base.sparse_graph.SparseGraphBackend>`
    and
    :class:`DenseGraphBackend <sage.graphs.base.dense_graph.DenseGraphBackend>`,
    which are fully functional backends. This class is mainly just for vertex
    functions, which are the same for both. A :class:`CGraphBackend` will not
    work on its own::

        sage: from sage.graphs.base.c_graph import CGraphBackend
        sage: CGB = CGraphBackend()
        sage: CGB.degree(0, True)
        Traceback (most recent call last):
        ...
        NotImplementedError: a derived class must return ``self._cg``

    The appropriate way to use these backends is via Sage graphs::

        sage: G = Graph(30)
        sage: G.add_edges([(0,1), (0,3), (4,5), (9, 23)])
        sage: G.edges(sort=True, labels=False)
        [(0, 1), (0, 3), (4, 5), (9, 23)]

    This class handles the labels of vertices and edges. For vertices it uses
    two dictionaries ``vertex_labels`` and ``vertex_ints``. They are just
    opposite of each other: ``vertex_ints`` makes a translation from label to
    integers (that are internally used) and ``vertex_labels`` make the
    translation from internally used integers to actual labels. This class tries
    hard to avoid translation if possible. This will work only if the graph is
    built on integers from `0` to `n-1` and the vertices are basically added in
    increasing order.

    .. SEEALSO::

        - :class:`SparseGraphBackend <sage.graphs.base.sparse_graph.SparseGraphBackend>`
          -- backend for sparse graphs.

        - :class:`DenseGraphBackend <sage.graphs.base.dense_graph.DenseGraphBackend>`
          -- backend for dense graphs.
    """

    ###################################
    # Basic Access
    ###################################

    cdef CGraph cg(self):
        r"""
        Return the attribute ``_cg`` casted into ``CGraph``.
        """
        raise NotImplementedError("a derived class must return ``self._cg``")

    def c_graph(self):
        r"""
        Return the ``._cg`` and ``._cg_rev`` attributes.

        .. NOTE::

            The ``._cg_rev`` attribute has been removed and hence ``None`` is returned.

        EXAMPLES::

            sage: cg,cg_rev = graphs.PetersenGraph()._backend.c_graph()
            sage: cg
            <sage.graphs.base.sparse_graph.SparseGraph object at ...>
        """
        return (self.cg(), None)

    def loops(self, new=None):
        """
        Check whether loops are allowed in this graph.

        INPUT:

        - ``new`` -- boolean (default: ``None``); to set or ``None`` to get

        OUTPUT:

        - If ``new=None``, return ``True`` if this graph allows self-loops or
          ``False`` if self-loops are not allowed

        - If ``new`` is a boolean, set the self-loop permission of this graph
          according to the boolean value of ``new``

        EXAMPLES::

            sage: G = Graph()
            sage: G._backend.loops()
            False
            sage: G._backend.loops(True)
            sage: G._backend.loops()
            True
        """
        if new is None:
            return self._loops
        if new:
            self._loops = True
        else:
            self._loops = False

    def num_edges(self, directed):
        """
        Return the number of edges in ``self``.

        INPUT:

        - ``directed`` -- boolean; whether to count ``(u, v)`` and ``(v, u)`` as
          one or two edges

        OUTPUT:

        - If ``directed=True``, counts the number of directed edges in this
          graph. Otherwise, return the size of this graph.

        .. SEEALSO::

            - :meth:`num_verts`
              -- return the order of this graph.

        EXAMPLES::

            sage: G = Graph(graphs.PetersenGraph())
            sage: G._backend.num_edges(False)
            15

        TESTS:

        Ensure that :issue:`8395` is fixed. ::

            sage: G = Graph({1:[1]}); G
            Looped graph on 1 vertex
            sage: G.edges(sort=True, labels=False)
            [(1, 1)]
            sage: G.size()
            1
            sage: G = Graph({1:[2,2]}); G
            Multi-graph on 2 vertices
            sage: G.edges(sort=True, labels=False)
            [(1, 2), (1, 2)]
            sage: G.size()
            2
            sage: G = Graph({1:[1,1]}); G
            Looped multi-graph on 1 vertex
            sage: G.edges(sort=True, labels=False)
            [(1, 1), (1, 1)]
            sage: G.size()
            2
            sage: D = DiGraph({1:[1]}); D
            Looped digraph on 1 vertex
            sage: D.edges(sort=True, labels=False)
            [(1, 1)]
            sage: D.size()
            1
            sage: D = DiGraph({1:[2,2], 2:[1,1]}); D
            Multi-digraph on 2 vertices
            sage: D.edges(sort=True, labels=False)
            [(1, 2), (1, 2), (2, 1), (2, 1)]
            sage: D.size()
            4
            sage: D = DiGraph({1:[1,1]}); D
            Looped multi-digraph on 1 vertex
            sage: D.edges(sort=True, labels=False)
            [(1, 1), (1, 1)]
            sage: D.size()
            2
            sage: from sage.graphs.base.sparse_graph import SparseGraphBackend
            sage: S = SparseGraphBackend(7)
            sage: S.num_edges(False)
            0
            sage: S.loops(True)
            sage: S.add_edge(1, 1, None, directed=False)
            sage: S.num_edges(False)
            1
            sage: S.multiple_edges(True)
            sage: S.add_edge(1, 1, None, directed=False)
            sage: S.num_edges(False)
            2
            sage: from sage.graphs.base.dense_graph import DenseGraphBackend
            sage: D = DenseGraphBackend(7)
            sage: D.num_edges(False)
            0
            sage: D.loops(True)
            sage: D.add_edge(1, 1, None, directed=False)
            sage: D.num_edges(False)
            1
        """
        if directed:
            return self.cg().num_arcs
        else:
            i = self.cg().num_arcs
            k = 0
            if self.loops(None):
                if self.multiple_edges(None):
                    for j in self.iterator_verts():
                        if self.has_edge(j, j, None):
                            k += len(self.get_edge_label(j, j))
                else:
                    for j in self.iterator_verts():
                        if self.has_edge(j, j, None):
                            k += 1
            i = (i - k) // 2
            return i + k

    def num_verts(self):
        """
        Return the number of vertices in ``self``.

        OUTPUT: the order of this graph

        .. SEEALSO::

            - :meth:`num_edges`
              -- return the number of (directed) edges in this graph.

        EXAMPLES::

            sage: G = Graph(graphs.PetersenGraph())
            sage: G._backend.num_verts()
            10
        """
        return self.cg().num_verts

    cdef bint _delete_edge_before_adding(self) noexcept:
        """
        Return whether we should delete edges before adding any.

        This is in particular required if the backend theoretically allows
        multiple edges but the graph should not have multiple edges.
        """
        return not self._multiple_edges

    ###################################
    # Vertex Functions
    ###################################

    cdef inline int get_vertex(self, u) except ? -2:
        """
        Return an ``int`` representing the arbitrary hashable vertex ``u``
        (whether or not ``u`` is actually in the graph), or ``-1`` if a new
        association must be made for ``u`` to be a vertex.

        TESTS:

        We check that the bug described in :issue:`8406` is gone::

            sage: # needs sage.rings.finite_rings
            sage: G = Graph()
            sage: R.<a> = GF(3**3)
            sage: S.<x> = R[]
            sage: G.add_vertex(a**2)
            sage: G.add_vertex(x)
            sage: G.vertices(sort=True)
            [a^2, x]

        And that the bug described in :issue:`9610` is gone::

            sage: n = 20
            sage: k = 3
            sage: g = DiGraph()
            sage: g.add_edges((i, Mod(i + j, n)) for i in range(n) for j in range(1, k + 1))
            sage: g.vertices(sort=True)
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
            sage: g.strongly_connected_components()
            [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]]

        The bug in :issue:`14967` and :issue:`14853` is fixed::

            sage: DiGraph({0: {}, 1/2: {}})
            Digraph on 2 vertices
            sage: A = Set([RDF.random_element(min=0, max=10) for k in range(10)])
            sage: G = Graph()
            sage: G.add_vertices(A)
            sage: Set(G.vertices(sort=False)) == A
            True
        """
        cdef dict vertex_ints = self.vertex_ints
        cdef dict vertex_labels = self.vertex_labels
        cdef CGraph G = self.cg()
        cdef long u_long
        if u in vertex_ints:
            return vertex_ints[u]
        try:
            u_long = pyobject_to_long(u)
        except Exception:
            return -1
        if u_long < 0 or u_long >= <long>G.active_vertices.size or u_long in vertex_labels:
            return -1
        return u_long

    cdef inline int get_vertex_checked(self, u) except ? -2:
        """
        As :meth:`get_vertex`, but return ``-1``,
        if ``u`` is not a vertex of the graph.
        """
        cdef int u_int = self.get_vertex(u)
        if u_int != -1 and bitset_in(self.cg().active_vertices, u_int):
            return u_int
        else:
            return -1

    cdef vertex_label(self, int u_int):
        """
        Return the object represented by ``u_int``, or ``None`` if this does not
        represent a vertex.
        """
        cdef dict vertex_labels = self.vertex_labels

        if u_int in vertex_labels:
            return vertex_labels[u_int]
        elif bitset_in(self.cg().active_vertices, u_int):
            return u_int
        else:
            return None

    cdef inline int check_labelled_vertex(self, u, bint reverse) except ? -1:
        """
        Return an ``int`` representing the arbitrary hashable vertex ``u``, and
        update, if necessary, the translation dict and list. Add a vertex if the
        label is new.
        """
        cdef CGraph G = self.cg()

        cdef int u_int = self.get_vertex(u)
        if u_int != -1:
            if not bitset_in(G.active_vertices, u_int):
                bitset_add(G.active_vertices, u_int)
                G.num_verts += 1
            return u_int

        u_int = bitset_first_in_complement(G.active_vertices)
        if u_int == -1:
            G.realloc(2 * G.active_vertices.size)
            return self.check_labelled_vertex(u, reverse)

        self.vertex_labels[u_int] = u
        self.vertex_ints[u] = u_int
        G.add_vertex(u_int)
        return u_int

    def has_vertex(self, v):
        """
        Check whether ``v`` is a vertex of ``self``.

        INPUT:

        - ``v`` -- any object

        OUTPUT: ``True`` if ``v`` is a vertex of this graph; ``False`` otherwise

        EXAMPLES::

            sage: from sage.graphs.base.sparse_graph import SparseGraphBackend
            sage: B = SparseGraphBackend(7)
            sage: B.has_vertex(6)
            True
            sage: B.has_vertex(7)
            False
        """
        cdef int v_int = self.get_vertex_checked(v)
        return v_int != -1

    def add_vertex(self, name):
        """
        Add a vertex to ``self``.

        INPUT:

        - ``name`` -- the vertex to be added (must be hashable). If ``None``,
          a new name is created

        OUTPUT:

        - If ``name = None``, the new vertex name is returned. ``None``
          otherwise.

        .. SEEALSO::

            - :meth:`add_vertices` -- add a bunch of vertices of this graph

            - :meth:`has_vertex` -- returns whether or not this graph has a
              specific vertex

        EXAMPLES::

            sage: D = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: D.add_vertex(10)
            sage: D.add_vertex([])
            Traceback (most recent call last):
            ...
            TypeError: unhashable type: 'list'

        ::

            sage: S = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: S.add_vertex(10)
            sage: S.add_vertex([])
            Traceback (most recent call last):
            ...
            TypeError: unhashable type: 'list'
        """
        retval = None
        if name is None:
            name = 0
            while (name in self.vertex_ints or
                   (name not in self.vertex_labels and
                    bitset_in(self.cg().active_vertices, <mp_bitcnt_t> name))):
                name += 1
            retval = name

        self.check_labelled_vertex(name, False)  # this will add the vertex

        return retval

    def add_vertices(self, vertices):
        """
        Add vertices to ``self``.

        INPUT:

        - ``vertices`` -- iterator of vertex labels; a new name is created, used
          and returned in the output list for all ``None`` values in
          ``vertices``

        OUTPUT:

        Generated names of new vertices if there is at least one ``None`` value
        present in ``vertices``. ``None`` otherwise.

        .. SEEALSO::

            - :meth:`add_vertex` -- add a vertex to this graph

        EXAMPLES::

            sage: D = sage.graphs.base.sparse_graph.SparseGraphBackend(1)
            sage: D.add_vertices([1, 2, 3])
            sage: D.add_vertices([None] * 4)
            [4, 5, 6, 7]

        ::

            sage: G = sage.graphs.base.sparse_graph.SparseGraphBackend(0)
            sage: G.add_vertices([0, 1])
            sage: list(G.iterator_verts(None))
            [0, 1]
            sage: list(G.iterator_edges([0, 1], True))
            []

        ::

            sage: import sage.graphs.base.dense_graph
            sage: D = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: D.add_vertices([10, 11, 12])
        """
        cdef int nones = 0
        for v in vertices:
            if v is not None:
                self.add_vertex(v)
            else:
                nones += 1

        new_names = []
        while nones:
            new_names.append(self.add_vertex(None))
            nones -= 1

        return new_names if new_names else None

    def del_vertex(self, v):
        """
        Delete a vertex in ``self``, failing silently if the vertex is not
        in the graph.

        INPUT:

        - ``v`` -- vertex to be deleted

        .. SEEALSO::

            - :meth:`del_vertices` -- delete a bunch of vertices from this graph

        EXAMPLES::

            sage: D = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: D.del_vertex(0)
            sage: D.has_vertex(0)
            False

        ::

            sage: S = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: S.del_vertex(0)
            sage: S.has_vertex(0)
            False
        """
        if not self.has_vertex(v):
            return
        cdef int v_int = self.get_vertex(v)

        # delete each arc incident with v and v
        self.cg().del_vertex(v_int)

        # add v to unused vertices
        if v_int in self.vertex_labels:
            self.vertex_ints.pop(v)
            self.vertex_labels.pop(v_int)

    def del_vertices(self, vertices):
        """
        Delete vertices from an iterable container.

        INPUT:

        - ``vertices`` -- iterator of vertex labels

        OUTPUT: same as for :meth:`del_vertex`

        .. SEEALSO::

            - :meth:`del_vertex` -- delete a vertex of this graph

        EXAMPLES::

            sage: import sage.graphs.base.dense_graph
            sage: D = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: D.del_vertices([7, 8])
            sage: D.has_vertex(7)
            False
            sage: D.has_vertex(6)
            True

        ::

            sage: D = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: D.del_vertices([1, 2, 3])
            sage: D.has_vertex(1)
            False
            sage: D.has_vertex(0)
            True
        """
        for v in vertices:
            self.del_vertex(v)

    def iterator_verts(self, verts=None):
        """
        Iterate over the vertices of ``self`` intersected with
        ``verts``.

        INPUT:

        - ``verts`` -- an iterable container of objects (default: ``None``)

        OUTPUT:

        - If ``verts=None``, return an iterator over all vertices of this graph

        - If ``verts`` is a single vertex of the graph, treat it as the
          container ``[verts]``

        - If ``verts`` is a iterable container of vertices, find the
          intersection of ``verts`` with the vertex set of this graph and return
          an iterator over the resulting intersection

        .. SEEALSO::

            - :meth:`iterator_nbrs`
              -- returns an iterator over the neighbors of a vertex.

        EXAMPLES::

            sage: P = Graph(graphs.PetersenGraph())
            sage: list(P._backend.iterator_verts(P))
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            sage: list(P._backend.iterator_verts())
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            sage: list(P._backend.iterator_verts([1, 2, 3]))
            [1, 2, 3]
            sage: list(P._backend.iterator_verts([1, 2, 10]))
            [1, 2]
        """
        cdef size_t i
        if verts is None:
            for x in self.vertex_ints:
                yield x
            i = bitset_first(self.cg().active_vertices)
            while i != <size_t>-1:
                if (i not in self.vertex_labels
                        and i not in self.vertex_ints):
                    yield i
                i = bitset_next(self.cg().active_vertices, i + 1)
            return

        try:
            hash(verts)
        except Exception:
            pass
        else:
            if self.has_vertex(verts):
                yield verts
                return

        for v in verts:
            if self.has_vertex(v):
                yield v

    def relabel(self, perm, directed):
        """
        Relabel the graph according to ``perm``.

        INPUT:

        - ``perm`` -- anything which represents a permutation as
          ``v --> perm[v]``, for example a dict or a list

        - ``directed`` -- ignored (this is here for compatibility with other
          backends)

        EXAMPLES::

            sage: G = Graph(graphs.PetersenGraph())
            sage: G._backend.relabel(range(9,-1,-1), False)
            sage: G.edges(sort=True)
            [(0, 2, None),
             (0, 3, None),
             (0, 5, None),
             (1, 3, None),
             (1, 4, None),
             (1, 6, None),
             (2, 4, None),
             (2, 7, None),
             (3, 8, None),
             (4, 9, None),
             (5, 6, None),
             (5, 9, None),
             (6, 7, None),
             (7, 8, None),
             (8, 9, None)]
        """
        cdef int i
        cdef dict new_vx_ints = {}
        cdef dict new_vx_labels = {}
        for v in self.iterator_verts(None):
            i = self.get_vertex(v)
            new_vx_ints[perm[v]] = i
            new_vx_labels[i] = perm[v]
        self.vertex_ints = new_vx_ints
        self.vertex_labels = new_vx_labels

    ###################################
    # Neighbor Functions
    ###################################

    def degree(self, v, directed):
        """
        Return the degree of the vertex ``v``.

        INPUT:

        - ``v`` -- a vertex of the graph

        - ``directed`` -- boolean; whether to take into account the
          orientation of this graph in counting the degree of ``v``

        OUTPUT: the degree of vertex ``v``

        EXAMPLES::

            sage: from sage.graphs.base.sparse_graph import SparseGraphBackend
            sage: B = SparseGraphBackend(7)
            sage: B.degree(3, False)
            0

        TESTS:

        Ensure that issue :issue:`8395` is fixed. ::

            sage: def my_add_edges(G, m, n):
            ....:     for i in range(m):
            ....:         u = randint(0, n)
            ....:         v = randint(0, n)
            ....:         G.add_edge(u, v)
            sage: G = Graph({1:[1]}); G
            Looped graph on 1 vertex
            sage: G.edges(sort=True, labels=False)
            [(1, 1)]
            sage: G.degree(); G.size()
            [2]
            1
            sage: sum(G.degree()) == 2 * G.size()
            True
            sage: G = Graph({1:[1,2], 2:[2,3]}, loops=True); G
            Looped graph on 3 vertices
            sage: my_add_edges(G, 100, 50)
            sage: sum(G.degree()) == 2 * G.size()
            True
            sage: G = Graph({1:[2,2], 2:[3]}); G
            Multi-graph on 3 vertices
            sage: G.edges(sort=True, labels=False)
            [(1, 2), (1, 2), (2, 3)]
            sage: G.degree(); G.size()
            [2, 3, 1]
            3
            sage: sum(G.degree()) == 2 * G.size()
            True
            sage: G.allow_loops(True); G
            Looped multi-graph on 3 vertices
            sage: my_add_edges(G, 100, 50)
            sage: sum(G.degree()) == 2 * G.size()
            True
            sage: D = DiGraph({1:[2], 2:[1,3]}); D
            Digraph on 3 vertices
            sage: D.edges(sort=True, labels=False)
            [(1, 2), (2, 1), (2, 3)]
            sage: D.degree(); D.size()
            [2, 3, 1]
            3
            sage: sum(D.degree()) == 2 * D.size()
            True
            sage: D.allow_loops(True); D
            Looped digraph on 3 vertices
            sage: my_add_edges(D, 100, 50)
            sage: sum(D.degree()) == 2 * D.size()
            True
            sage: D.allow_multiple_edges(True)
            sage: my_add_edges(D, 200, 50)
            sage: sum(D.degree()) == 2 * D.size()
            True
            sage: G = Graph({1:[2,2,2]})
            sage: G.allow_loops(True)
            sage: G.add_edge(1,1)
            sage: G.add_edge(1,1)
            sage: G.edges(sort=True, labels=False)
            [(1, 1), (1, 1), (1, 2), (1, 2), (1, 2)]
            sage: G.degree(1)
            7
            sage: G.allow_loops(False)
            sage: G.edges(sort=True, labels=False)
            [(1, 2), (1, 2), (1, 2)]
            sage: G.degree(1)
            3
            sage: G = Graph({1:{2:['a','a','a']}})
            sage: G.allow_loops(True)
            sage: G.add_edge(1,1,'b')
            sage: G.add_edge(1,1,'b')
            sage: G.add_edge(1,1)
            sage: G.add_edge(1,1)
            sage: G.edges_incident()
            [(1, 1, None), (1, 1, None), (1, 1, 'b'), (1, 1, 'b'), (1, 2, 'a'), (1, 2, 'a'), (1, 2, 'a')]
            sage: G.degree(1)
            11
            sage: G.allow_loops(False)
            sage: G.edges(sort=True)
            [(1, 2, 'a'), (1, 2, 'a'), (1, 2, 'a')]
            sage: G.degree(1)
            3
            sage: G = Graph({1:{2:['a','a','a']}})
            sage: G.allow_loops(True)
            sage: G.add_edge(1,1,'b')
            sage: G.add_edge(1,1,'b')
            sage: G.edges(sort=True)
            [(1, 1, 'b'), (1, 1, 'b'), (1, 2, 'a'), (1, 2, 'a'), (1, 2, 'a')]
            sage: G.degree(1)
            7
            sage: G.allow_loops(False)
            sage: G.edges(sort=True)
            [(1, 2, 'a'), (1, 2, 'a'), (1, 2, 'a')]
            sage: G.degree(1)
            3

        Ensure that :issue:`13664` is fixed ::

            sage: W = WeylGroup(["A",1])                                                # needs sage.combinat sage.groups
            sage: G = W.cayley_graph()                                                  # needs sage.combinat sage.groups
            sage: Graph(G).degree()                                                     # needs sage.combinat sage.groups
            [1, 1]
            sage: h = Graph()
            sage: h.add_edge(1,2,"a")
            sage: h.add_edge(1,2,"a")
            sage: h.degree()
            [1, 1]
        """
        cdef int v_int = self.get_vertex(v)
        if directed:
            return self.cg().in_degrees[v_int] + self.cg().out_degrees[v_int]
        cdef int d = 0
        if self._loops and self.has_edge(v, v, None):
            if self._multiple_edges:
                d += len(self.get_edge_label(v, v))
            else:
                d += 1
        return self.cg().out_degrees[v_int] + d

    def out_degree(self, v):
        r"""
        Return the out-degree of ``v``.

        INPUT:

        - ``v`` -- a vertex of the graph

        EXAMPLES::


            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: D.out_degree(1)
            2
        """
        cdef int v_int = self.get_vertex(v)
        if self._directed:
            return self.cg().out_degrees[v_int]
        cdef int d = 0
        if self._loops and self.has_edge(v, v, None):
            if self._multiple_edges:
                d += len(self.get_edge_label(v, v))
            else:
                d += 1

        return self.cg().out_degrees[v_int] + d

    def in_degree(self, v):
        r"""
        Return the in-degree of ``v``.

        INPUT:

        - ``v`` -- a vertex of the graph

        EXAMPLES::


            sage: D = DiGraph( { 0: [1,2,3], 1: [0,2], 2: [3], 3: [4], 4: [0,5], 5: [1] } )
            sage: D.out_degree(1)
            2
        """
        if not self._directed:
            return self.out_degree(v)

        cdef int v_int = self.get_vertex(v)

        return self.cg().in_degrees[v_int]

    def iterator_nbrs(self, v):
        """
        Iterate over the neighbors of ``v``.

        INPUT:

        - ``v`` -- a vertex of this graph

        OUTPUT: an iterator over the neighbors the vertex ``v``

        .. SEEALSO::

            - :meth:`iterator_in_nbrs`
              -- returns an iterator over the in-neighbors of a vertex

            - :meth:`iterator_out_nbrs`
              -- returns an iterator over the out-neighbors of a vertex

            - :meth:`iterator_verts`
              -- returns an iterator over a given set of vertices

        EXAMPLES::

            sage: P = Graph(graphs.PetersenGraph())
            sage: list(P._backend.iterator_nbrs(0))
            [1, 4, 5]
            sage: Q = DiGraph(P)
            sage: list(Q._backend.iterator_nbrs(0))
            [1, 4, 5]
        """
        if not self._directed:
            yield from self.iterator_out_nbrs(v)
            return

        cdef int u_int
        cdef int v_int = self.get_vertex(v)
        if v_int == -1 or not bitset_in(self.cg().active_vertices, v_int):
            raise LookupError("vertex ({0}) is not a vertex of the graph".format(v))

        cdef set seen = set()
        for u_int in self.cg().in_neighbors(v_int):
            if u_int not in seen:
                yield self.vertex_label(u_int)
                seen.add(u_int)
        for u_int in self.cg().out_neighbors(v_int):
            if u_int not in seen:
                yield self.vertex_label(u_int)
                seen.add(u_int)
        return

    def iterator_in_nbrs(self, v):
        """
        Iterate over the incoming neighbors of ``v``.

        INPUT:

        - ``v`` -- a vertex of this graph

        OUTPUT: an iterator over the in-neighbors of the vertex ``v``

        .. SEEALSO::

            - :meth:`iterator_nbrs`
              -- returns an iterator over the neighbors of a vertex

            - :meth:`iterator_out_nbrs`
              -- returns an iterator over the out-neighbors of a vertex

        EXAMPLES::

            sage: P = DiGraph(graphs.PetersenGraph().to_directed())
            sage: list(P._backend.iterator_in_nbrs(0))
            [1, 4, 5]

        TESTS::

            sage: P = DiGraph(graphs.PetersenGraph().to_directed())
            sage: list(P._backend.iterator_in_nbrs(63))
            Traceback (most recent call last):
            ...
            LookupError: vertex (63) is not a vertex of the graph
        """

        cdef int u_int
        cdef int v_int = self.get_vertex(v)
        if v_int == -1 or not bitset_in(self.cg().active_vertices, v_int):
            raise LookupError("vertex ({0}) is not a vertex of the graph".format(v))

        for u_int in self.cg().in_neighbors(v_int):
            yield self.vertex_label(u_int)

    def iterator_out_nbrs(self, v):
        """
        Iterate over the outgoing neighbors of ``v``.

        INPUT:

        - ``v`` -- a vertex of this graph

        OUTPUT: an iterator over the out-neighbors of the vertex ``v``

        .. SEEALSO::

            - :meth:`iterator_nbrs`
              -- returns an iterator over the neighbors of a vertex

            - :meth:`iterator_in_nbrs`
              -- returns an iterator over the in-neighbors of a vertex

        EXAMPLES::

            sage: P = DiGraph(graphs.PetersenGraph().to_directed())
            sage: list(P._backend.iterator_out_nbrs(0))
            [1, 4, 5]

        TESTS::

            sage: P = DiGraph(graphs.PetersenGraph().to_directed())
            sage: list(P._backend.iterator_out_nbrs(-41))
            Traceback (most recent call last):
            ...
            LookupError: vertex (-41) is not a vertex of the graph
        """
        cdef int u_int
        cdef int v_int = self.get_vertex(v)
        if v_int == -1 or not bitset_in(self.cg().active_vertices, v_int):
            raise LookupError("vertex ({0}) is not a vertex of the graph".format(v))

        for u_int in self.cg().out_neighbors(v_int):
            yield self.vertex_label(u_int)

    ###################################
    # Edge Functions
    ###################################

    cdef int new_edge_label(self, object l) except -1:
        raise NotImplementedError()

    def add_edges(self, object edges, bint directed, bint remove_loops=False):
        """
        Add edges from a list.

        INPUT:

        - ``edges`` -- the edges to be added; can either be of the form
          ``(u,v)`` or ``(u,v,l)``

        - ``directed`` -- if ``False``, add ``(v,u)`` as well as ``(u,v)``

        - ``remove_loops`` -- if ``True``, remove loops

        EXAMPLES::

            sage: D = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: D.add_edges([(0,1), (2,3), (4,5), (5,6)], False)
            sage: list(D.iterator_edges(range(9), True))
            [(0, 1, None),
             (2, 3, None),
             (4, 5, None),
             (5, 6, None)]
        """
        cdef object u, v, l, e
        for e in edges:
            if len(e) == 3:
                u, v, l = e
            else:
                u, v = e
                l = None
            if unlikely(remove_loops and u == v):
                continue
            self.add_edge(u, v, l, directed)

    cpdef add_edge(self, object u, object v, object l, bint directed):
        """
        Add the edge ``(u,v)`` to ``self``.

        INPUT:

        - ``u``, ``v`` -- the vertices of the edge

        - ``l`` -- the edge label

        - ``directed`` -- if ``False``, also add ``(v,u)``

        .. NOTE::

            The input ``l`` is ignored if the backend
            does not support labels.

        EXAMPLES::

            sage: D = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: D.add_edge(0,1,None,False)
            sage: list(D.iterator_edges(range(9), True))
            [(0, 1, None)]

        ::

            sage: D = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: D.add_edge(0, 1, None, False)
            sage: list(D.iterator_edges(range(9), True))
            [(0, 1, None)]

        TESTS::

            sage: D = DiGraph(sparse=True)
            sage: D.add_edge(0,1,2)
            sage: D.add_edge(0,1,3)
            sage: D.edges(sort=True)
            [(0, 1, 3)]

        Check :issue:`22991` for sparse backend::

            sage: G = Graph(3, sparse=True)
            sage: G.add_edge(0,0)
            Traceback (most recent call last):
            ...
            ValueError: cannot add edge from 0 to 0 in graph without loops
            sage: G = Graph(3, sparse=True, loops=True)
            sage: G.add_edge(0,0); G.edges(sort=True)
            [(0, 0, None)]

        Check :issue:`22991` for dense backend::

            sage: G = Graph(3, sparse=False)
            sage: G.add_edge(0,0)
            Traceback (most recent call last):
            ...
            ValueError: cannot add edge from 0 to 0 in graph without loops
            sage: G = Graph(3, sparse=True, loops=True)
            sage: G.add_edge(0, 0); G.edges(sort=True)
            [(0, 0, None)]

        Remove edges correctly when multiedges are not allowed (:issue:`28077`)::

            sage: D = DiGraph(multiedges=False)
            sage: D.add_edge(1, 2, 'A')
            sage: D.add_edge(1, 2, 'B')
            sage: D.delete_edge(1, 2)
            sage: D.incoming_edges(2)
            []
            sage: D.shortest_path(1, 2)
            []
        """
        if u is None:
            u = self.add_vertex(None)
        if v is None:
            v = self.add_vertex(None)

        cdef int u_int = self.check_labelled_vertex(u, False)
        cdef int v_int = self.check_labelled_vertex(v, False)

        cdef CGraph cg = self.cg()

        cdef int l_int
        if l is None:
            l_int = 0
        else:
            l_int = self.new_edge_label(l)

        if u_int == v_int and not self._loops:
            raise ValueError(f"cannot add edge from {u!r} to {v!r} in graph without loops")

        if self._delete_edge_before_adding():
            if cg.has_arc_label(u_int, v_int, l_int):
                return
            else:
                cg.del_all_arcs(u_int, v_int)
                if not directed and self._directed and v_int != u_int:
                    cg.del_all_arcs(v_int, u_int)

        cg.add_arc_label_unsafe(u_int, v_int, l_int)
        if not directed and self._directed and v_int != u_int:
            cg.add_arc_label_unsafe(v_int, u_int, l_int)

    def del_edges(self, object edges, bint directed):
        """
        Delete edges from a list.

        INPUT:

        - ``edges`` -- the edges to be added; can either be of the form
          ``(u,v)`` or ``(u,v,l)``

        - ``directed`` -- if ``False``, remove``(v,u)`` as well as ``(u,v)``

        EXAMPLES::

            sage: D = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: D.add_edges([(0,1), (2,3), (4,5), (5,6)], False)
            sage: D.del_edges([(0,1), (2,3), (4,5), (5,6)], False)
            sage: list(D.iterator_edges(range(9), True))
            []
        """
        cdef object u, v, l, e
        for e in edges:
            if len(e) == 3:
                u, v, l = e
            else:
                u, v = e
                l = None
            self.del_edge(u, v, l, directed)

    cpdef del_edge(self, object u, object v, object l, bint directed):
        """
        Delete edge ``(u, v, l)``.

        INPUT:

        - ``u``, ``v`` -- the vertices of the edge

        - ``l`` -- the edge label

        - ``directed`` -- if ``False``, also delete ``(v, u, l)``

        .. NOTE::

            The input ``l`` is ignored if the backend
            does not support labels.

        EXAMPLES::

            sage: D = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
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

        ::

            sage: D = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: D.add_edges([(0, 1), (2, 3), (4, 5), (5, 6)], False)
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

        TESTS::

            sage: G = Graph(sparse=True)
            sage: G.add_edge(0,1,2)
            sage: G.delete_edge(0,1)
            sage: G.edges(sort=True)
            []

            sage: G = Graph(multiedges=True, sparse=True)
            sage: G.add_edge(0,1,2)
            sage: G.add_edge(0,1,None)
            sage: G.delete_edge(0,1)
            sage: G.edges(sort=True)
            [(0, 1, 2)]

        Do we remove loops correctly? (:issue:`12135`)::

            sage: g=Graph({0:[0,0,0]}, sparse=True)
            sage: g.edges(sort=True, labels=False)
            [(0, 0), (0, 0), (0, 0)]
            sage: g.delete_edge(0,0); g.edges(sort=True, labels=False)
            [(0, 0), (0, 0)]
        """
        cdef int u_int = self.get_vertex_checked(u)
        cdef int v_int = self.get_vertex_checked(v)

        if u_int == -1 or v_int == -1:
            return

        cdef CGraph cg = self.cg()

        if l is None:
            if cg.has_arc_label_unsafe(u_int, v_int, 0):
                l_int = 0
            else:
                l_int = cg.arc_label_unsafe(u_int, v_int)
        elif not self._multiple_edges:
            l_int = cg.arc_label_unsafe(u_int, v_int)
            if not l_int or not self.edge_labels[l_int] == l:
                # The requested edge does not exist.
                return
        else:
            for l_int in cg.all_arcs(u_int, v_int):
                if l_int and self.edge_labels[l_int] == l:
                    break
            else:
                return

        cg.del_arc_label(u_int, v_int, l_int)
        if not directed and self._directed and v_int != u_int:
            cg.del_arc_label(v_int, u_int, l_int)
        self.free_edge_label(l_int)

    cdef bint _has_labeled_edge_unsafe(self, int u_int, int v_int, object l) except -1:
        """
        Return whether ``self`` has an arc specified by indices of the vertices
        and an arc label.
        """
        raise NotImplementedError

    cdef int free_edge_label(self, int l_int) except -1:
        raise NotImplementedError()

    cdef list _all_edge_labels(self, int u, int v, uint32_t* edge=NULL):
        """
        Give the labels of all arcs from ``u`` to ``v``.

        ``u`` and ``v`` are the integers corresponding to vertices.

        ``edge`` may point to an edge from ``u`` to ``v``.
        """
        cdef int l_int
        return [self.edge_labels[l_int] if l_int else None for l_int in self.cg().all_arcs(u, v)]

    ###################################
    # Edge Iterators
    ###################################

    def iterator_edges(self, object vertices, bint labels):
        """
        Iterate over the edges incident to a sequence of vertices.

        Edges are assumed to be undirected.

        .. WARNING::

            This will try to sort the two ends of every edge.

        INPUT:

        - ``vertices`` -- list of vertex labels

        - ``labels`` -- boolean, whether to return labels as well

        EXAMPLES::

            sage: G = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: G.add_edge(1,2,3,False)
            sage: list(G.iterator_edges(range(9), False))
            [(1, 2)]
            sage: list(G.iterator_edges(range(9), True))
            [(1, 2, 3)]

        TESTS::

            sage: g = graphs.PetersenGraph()
            sage: g.edges_incident([0,1,2])
            [(0, 1, None),
             (0, 4, None),
             (0, 5, None),
             (1, 2, None),
             (1, 6, None),
             (2, 3, None),
             (2, 7, None)]
        """
        yield from self._iterator_edges(vertices, labels, modus=3)

    def iterator_unsorted_edges(self, object vertices, bint labels):
        """
        Iterate over the edges incident to a sequence of vertices.

        Edges are assumed to be undirected.

        This does not sort the ends of each edge.

        INPUT:

        - ``vertices`` -- list of vertex labels

        - ``labels`` -- boolean, whether to return labels as well

        EXAMPLES::

            sage: G = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: G.add_edge(1,2,3,False)
            sage: list(G.iterator_unsorted_edges(range(9), False))
            [(2, 1)]
            sage: list(G.iterator_unsorted_edges(range(9), True))
            [(2, 1, 3)]

        TESTS::

            sage: G = Graph(sparse=True)
            sage: G.add_edge((1,'a'))
            sage: list(G._backend.iterator_unsorted_edges([1, 'a'],False))
            [(1, 'a')]
        """
        yield from self._iterator_edges(vertices, labels, modus=2)

    def iterator_out_edges(self, object vertices, bint labels):
        """
        Iterate over the outbound edges incident to a sequence of vertices.

        INPUT:

        - ``vertices`` -- list of vertex labels

        - ``labels`` -- boolean, whether to return labels as well

        EXAMPLES::

            sage: G = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: G.add_edge(1,2,3,True)
            sage: list(G.iterator_out_edges([2], False))
            []
            sage: list(G.iterator_out_edges([1], False))
            [(1, 2)]
            sage: list(G.iterator_out_edges([1], True))
            [(1, 2, 3)]
        """
        yield from self._iterator_edges(vertices, labels, modus=0)

    def iterator_in_edges(self, object vertices, bint labels):
        """
        Iterate over the incoming edges incident to a sequence of vertices.

        INPUT:

        - ``vertices`` -- list of vertex labels

        - ``labels`` -- boolean, whether to return labels as well

        EXAMPLES::

            sage: G = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: G.add_edge(1,2,3,True)
            sage: list(G.iterator_in_edges([1], False))
            []
            sage: list(G.iterator_in_edges([2], False))
            [(1, 2)]
            sage: list(G.iterator_in_edges([2], True))
            [(1, 2, 3)]
        """
        yield from self._iterator_edges(vertices, labels, modus=1)

    def _iterator_edges(self, object vertices, const bint labels, const int modus=0):
        """
        Iterate over the edges incident to a sequence of vertices.

        INPUT:

        - ``vertices`` -- list of vertex labels

        - ``labels`` -- boolean, whether to return labels as well

        - ``modus`` -- integer representing the modus of the iterator:
          - ``0`` -- outgoing edges
          - ``1`` -- ingoing edges
          - ``2`` -- unsorted edges of an undirected graph
          - ``3`` -- sorted edges of an undirected graph

        EXAMPLES::

            sage: G = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: G.add_edge(1, 2, None, False)
            sage: list(G._iterator_edges(range(9), False, 3))
            [(1, 2)]
            sage: list(G._iterator_edges(range(9), True, 3))
            [(1, 2, None)]

        ::

            sage: G = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: G.add_edge(1, 2, None, True)
            sage: list(G.iterator_in_edges([1], False))
            []
            sage: list(G.iterator_in_edges([2], False))
            [(1, 2)]
            sage: list(G.iterator_in_edges([2], True))
            [(1, 2, None)]

        ::

            sage: G = sage.graphs.base.dense_graph.DenseGraphBackend(9)
            sage: G.add_edge(1, 2, None, True)
            sage: list(G.iterator_out_edges([2], False))
            []
            sage: list(G.iterator_out_edges([1], False))
            [(1, 2)]
            sage: list(G.iterator_out_edges([1], True))
            [(1, 2, None)]
        """
        cdef object u, v, l, v_copy
        cdef int u_int, v_int, l_int, foo
        cdef CGraph cg = self.cg()
        cdef list b_vertices_2, all_arc_labels
        cdef FrozenBitset b_vertices
        cdef bint out = modus == 0

        cdef int vertices_case
        cdef object it

        if not isinstance(vertices, list):
            # ALL edges
            it = self.iterator_verts(None)
            vertices_case = 0

        elif not vertices:
            return

        elif len(vertices) == 1:
            # One vertex
            vertices_case = 1
            v_int = -1

        else:
            # Several vertices (nonempty list)
            vertices_case = 2
            b_vertices_2 = [self.get_vertex_checked(v) for v in vertices]
            try:
                b_vertices = FrozenBitset(foo for foo in b_vertices_2 if foo >= 0)
            except ValueError:
                # Avoiding "Bitset must not be empty"
                # in case none of the vertices is active.
                return
            it = iter(b_vertices)

        while True:
            # Think of this as a loop through ``vertices``.
            # We pick the next vertex according to three cases.

            if vertices_case == 0:
                # ALL edges
                try:
                    v = next(it)
                    v_int = self.get_vertex(v)
                except StopIteration:
                    return

            elif vertices_case == 1:
                # One vertex
                if v_int != -1:
                    # Only visit one vertex once.
                    return
                v = vertices[0]
                v_int = self.get_vertex_checked(v)
                if v_int == -1:
                    return

            else:
                # Several vertices (nonempty list)
                try:
                    v_int = -1
                    while v_int == -1:
                        v_int = next(it)
                    v = self.vertex_label(v_int)
                except StopIteration:
                    return

            # WARNING
            # If you modify this, you must keep in mind the documentation in the
            # corresponding method in ``generic_graph.py`` in the method ``edge_iterator``.
            # E.g. code assumes that you can use an iterator to relabel or delete arcs.

            u_int = cg._next_neighbor_unsafe(v_int, -1, out, &l_int)
            while u_int != -1:
                if (modus < 2 or                                            # Do not delete duplicates.
                        vertices_case == 1 or                               # Only one vertex, so no duplicates.
                        u_int >= v_int or                                   # We visit if u_int >= v_int ...
                        (vertices_case == 2 and
                            u_int < b_vertices.capacity() and
                            not bitset_in(b_vertices._bitset, u_int))):     # ... or if u_int is not in ``vertices``.
                    u = self.vertex_label(u_int)
                    if labels:
                        l = self.edge_labels[l_int] if l_int else None

                    # Yield the arc/arcs.
                    v_copy = v
                    if _reorganize_edge(v, u, modus):
                        u, v = v, u

                    if not self._multiple_edges:
                        if labels:
                            yield (v, u, l)
                        else:
                            yield (v, u)
                    else:
                        if out:
                            all_arc_labels = cg.all_arcs(v_int, u_int)
                        else:
                            all_arc_labels = cg.all_arcs(u_int, v_int)

                        for l_int in all_arc_labels:
                            if labels:
                                l = self.edge_labels[l_int] if l_int else None
                                yield (v, u, l)
                            else:
                                yield (v, u)
                    v = v_copy

                if unlikely(not bitset_in(self.cg().active_vertices, v_int)):
                    raise IndexError("the vertices were modified while iterating the edges")

                u_int = cg._next_neighbor_unsafe(v_int, u_int, out, &l_int)

    ###################################
    # Using Edge Iterators
    ###################################

    def is_subgraph(self, CGraphBackend other, object vertices, bint ignore_labels=False):
        """
        Return whether the subgraph of ``self`` induced by ``vertices`` is a subgraph of ``other``.

        If ``vertices`` are the vertices of ``self``, return whether ``self`` is a subgraph of ``other``.

        INPUT:

            - ``other`` -- a subclass of :class:`CGraphBackend`
            - ``vertices`` -- a iterable over the vertex labels
            - ``ignore_labels`` -- boolean (default: ``False``); whether to ignore the labels

        EXAMPLES::

            sage: G = sage.graphs.base.dense_graph.DenseGraphBackend(4, directed=True)
            sage: H = sage.graphs.base.dense_graph.DenseGraphBackend(4, directed=True)
            sage: G.add_edges([[0,1],[0,2],[0,3],[1,2]], True)
            sage: H.add_edges([[0,1],[0,2],[0,3]], True)
            sage: G.is_subgraph(H, range(4))
            False
            sage: H.is_subgraph(G, range(4))
            True
            sage: G.is_subgraph(H, [0,1,3])
            True

        Ignore the labels or not::

            sage: G = sage.graphs.base.sparse_graph.SparseGraphBackend(3, directed=True)
            sage: G.multiple_edges(True)
            sage: H = sage.graphs.base.sparse_graph.SparseGraphBackend(3, directed=True)
            sage: H.multiple_edges(True)
            sage: G.add_edges([[0,1,'a'], [0,1,'b'], [0,2,'c'], [0,2,'d'], [0,2,'e']], True)
            sage: H.add_edges([[0,1,'a'], [0,1,'foo'], [0,2,'c'], [0,2,'d'], [0,2,'e'], [0,2,'e']], True)
            sage: G.is_subgraph(H, range(3))
            False
            sage: G.is_subgraph(H, range(3), ignore_labels=True)
            True

        Multiplicities of edges are considered::

            sage: G.is_subgraph(H, [0,2])
            True
            sage: H.is_subgraph(G, [0,2])
            False
        """
        if not ignore_labels:
            return 1 == self._use_edge_iterator_on_subgraph(other, vertices, 1)
        else:
            return 1 == self._use_edge_iterator_on_subgraph(other, vertices, 2)

    def subgraph_given_vertices(self, CGraphBackend other, object vertices):
        """
        Initialize ``other`` to be the subgraph of ``self`` with given vertices.

        INPUT:

        - ``other`` -- a (mutable) subclass of :class:`CGraphBackend`
        - ``vertices`` -- list of vertex labels

        .. NOTE::

            ``other`` is assumed to be the empty graph.

        EXAMPLES:

        Make a dense copy::

            sage: G = sage.graphs.base.dense_graph.DenseGraphBackend(9, directed=True)
            sage: G.loops(True)
            sage: G.add_edges([[0,1], [1,2], [2,3], [3,4], [4,5], [5,6], [7,8], [3,3]], True)
            sage: H = sage.graphs.base.dense_graph.DenseGraphBackend(0, directed=True)
            sage: H.loops(True)
            sage: G.subgraph_given_vertices(H, range(9))
            sage: list(H.iterator_out_edges(list(range(9)), False)) == list(G.iterator_out_edges(list(range(9)), False))
            True

        Make a sparse copy::

            sage: H = sage.graphs.base.sparse_graph.SparseGraphBackend(0, directed=True)
            sage: H.loops(True)
            sage: G.subgraph_given_vertices(H, range(9))
            sage: sorted(list(H.iterator_out_edges(list(range(9)), False))) == sorted(list(G.iterator_out_edges(list(range(9)), False)))
            True

        Initialize a proper subgraph::

            sage: H = sage.graphs.base.sparse_graph.SparseGraphBackend(0, directed=True)
            sage: H.loops(True)
            sage: G.subgraph_given_vertices(H, [2,3,4,5])
            sage: list(H.iterator_out_edges(list(range(9)), False))
            [(2, 3), (3, 3), (3, 4), (4, 5)]

        Loops are removed, if the other graph does not allow loops::

            sage: H = sage.graphs.base.sparse_graph.SparseGraphBackend(0, directed=True)
            sage: H.loops(False)
            sage: G.subgraph_given_vertices(H, [2,3,4,5])
            sage: list(H.iterator_out_edges(list(range(9)), False))
            [(2, 3), (3, 4), (4, 5)]

        Multiple edges and labels are copied::

            sage: G = sage.graphs.base.sparse_graph.SparseGraphBackend(4, directed=False)
            sage: G.multiple_edges(True)
            sage: G.add_edges([[0,1,'a'], [1,2,'b'], [2,3,'c'], [0,1,'d']], False)
            sage: H = sage.graphs.base.sparse_graph.SparseGraphBackend(0, directed=False)
            sage: H.multiple_edges(True)
            sage: G.subgraph_given_vertices(H, [0,1,2])
            sage: list(H.iterator_edges(list(range(4)), True))
            [(0, 1, 'a'), (0, 1, 'd'), (1, 2, 'b')]

        Multiple edges are removed, if the other graph does not allow them::

            sage: H = sage.graphs.base.sparse_graph.SparseGraphBackend(0, directed=False)
            sage: H.multiple_edges(False)
            sage: G.subgraph_given_vertices(H, [0,1,2])
            sage: list(H.iterator_edges(list(range(4)), True))
            [(0, 1, 'd'), (1, 2, 'b')]

        Labels are removed, if the other graph does not allow them::

            sage: H = sage.graphs.base.dense_graph.DenseGraphBackend(0, directed=False)
            sage: G.subgraph_given_vertices(H, [0,1,2])
            sage: list(H.iterator_edges(list(range(4)), True))
            [(0, 1, None), (1, 2, None)]

        A directed subgraph of an undirected graph is taken by initializing
        with edges in both directions::

            sage: G = sage.graphs.base.sparse_graph.SparseGraphBackend(4, directed=True)
            sage: G.loops(True)
            sage: G.multiple_edges(True)
            sage: G.add_edges([[0,1,'a'], [1,2,'b'], [2,3,'c'], [0,1,'d'], [2,2,'e']], False)
            sage: H = sage.graphs.base.sparse_graph.SparseGraphBackend(0, directed=True)
            sage: H.multiple_edges(True)
            sage: H.loops(True)
            sage: G.subgraph_given_vertices(H, [0,1,2])
            sage: list(H.iterator_out_edges(list(range(4)), True))
            [(0, 1, 'a'),
             (0, 1, 'd'),
             (1, 0, 'a'),
             (1, 0, 'd'),
             (1, 2, 'b'),
             (2, 1, 'b'),
             (2, 2, 'e')]

        An undirected subgraph of a directeed graph is not defined::

            sage: G = sage.graphs.base.sparse_graph.SparseGraphBackend(4, directed=True)
            sage: G.add_edges([[0,1,'a'], [1,2,'b'], [2,3,'c']], False)
            sage: H = sage.graphs.base.sparse_graph.SparseGraphBackend(0, directed=False)
            sage: G.subgraph_given_vertices(H, [0,1,2])
            Traceback (most recent call last):
            ...
            ValueError: cannot obtain an undirected subgraph of a directed graph

        TESTS:

        All the examples for ``self`` a static sparse graph.

        Make a dense copy::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            sage: G = Graph(loops=True)
            sage: G.add_edges([[0,1], [1,2], [2,3], [3,4], [4,5], [5,6], [7,8], [3,3]])
            sage: G = StaticSparseBackend(G)
            sage: H = sage.graphs.base.dense_graph.DenseGraphBackend(0, directed=False)
            sage: H.loops(True)
            sage: G.subgraph_given_vertices(H, range(9))
            sage: list(H.iterator_edges(list(range(9)), False)) == list(G.iterator_edges(list(range(9)), False))
            True

        Make a sparse copy::

            sage: H = sage.graphs.base.sparse_graph.SparseGraphBackend(0, directed=False)
            sage: H.loops(True)
            sage: G.subgraph_given_vertices(H, range(9))
            sage: sorted(list(H.iterator_edges(list(range(9)), False))) == sorted(list(G.iterator_edges(list(range(9)), False)))
            True

        Initialize a proper subgraph::

            sage: H = sage.graphs.base.sparse_graph.SparseGraphBackend(0, directed=False)
            sage: H.loops(True)
            sage: G.subgraph_given_vertices(H, [2,3,4,5])
            sage: list(H.iterator_edges(list(range(9)), False))
            [(2, 3), (3, 3), (3, 4), (4, 5)]

        Loops are removed, if the other graph does not allow loops::

            sage: H = sage.graphs.base.sparse_graph.SparseGraphBackend(0, directed=False)
            sage: H.loops(False)
            sage: G.subgraph_given_vertices(H, [2,3,4,5])
            sage: list(H.iterator_edges(list(range(9)), False))
            [(2, 3), (3, 4), (4, 5)]

        Multiple edges and labels are copied::

            sage: G = Graph(multiedges=True)
            sage: G.add_edges([[0,1,'a'], [1,2,'b'], [2,3,'c'], [0,1,'d']], False)
            sage: G = StaticSparseBackend(G)
            sage: H = sage.graphs.base.sparse_graph.SparseGraphBackend(0, directed=False)
            sage: H.multiple_edges(True)
            sage: G.subgraph_given_vertices(H, [0,1,2])
            sage: list(H.iterator_edges(list(range(4)), True))
            [(0, 1, 'a'), (0, 1, 'd'), (1, 2, 'b')]

        Multiple edges are removed, if the other graph does not allow them::

            sage: H = sage.graphs.base.sparse_graph.SparseGraphBackend(0, directed=False)
            sage: H.multiple_edges(False)
            sage: G.subgraph_given_vertices(H, [0,1,2])
            sage: list(H.iterator_edges(list(range(4)), True))
            [(0, 1, 'a'), (1, 2, 'b')]

        Labels are removed, if the other graph does not allow them::

            sage: H = sage.graphs.base.dense_graph.DenseGraphBackend(0, directed=False)
            sage: G.subgraph_given_vertices(H, [0,1,2])
            sage: list(H.iterator_edges(list(range(4)), True))
            [(0, 1, None), (1, 2, None)]

        A directed subgraph of an undirected graph is taken by initializing
        with edges in both directions::

            sage: G = Graph(multiedges=True, loops=True)
            sage: G.add_edges([[0,1,'a'], [1,2,'b'], [2,3,'c'], [0,1,'d'], [2,2,'e']])
            sage: G = StaticSparseBackend(G)
            sage: H = sage.graphs.base.sparse_graph.SparseGraphBackend(0, directed=True)
            sage: H.multiple_edges(True)
            sage: H.loops(True)
            sage: G.subgraph_given_vertices(H, [0,1,2])
            sage: list(H.iterator_out_edges(list(range(4)), True))
            [(0, 1, 'a'),
             (0, 1, 'd'),
             (1, 0, 'a'),
             (1, 0, 'd'),
             (1, 2, 'b'),
             (2, 1, 'b'),
             (2, 2, 'e')]

        An undirected subgraph of a directeed graph is not defined::

            sage: G = DiGraph()
            sage: G.add_edges([[0,1,'a'], [1,2,'b'], [2,3,'c']])
            sage: G = StaticSparseBackend(G)
            sage: H = sage.graphs.base.sparse_graph.SparseGraphBackend(0, directed=False)
            sage: G.subgraph_given_vertices(H, [0,1,2])
            Traceback (most recent call last):
            ...
            ValueError: cannot obtain an undirected subgraph of a directed graph
        """
        self._use_edge_iterator_on_subgraph(other, vertices, 0)

    cdef int _use_edge_iterator_on_subgraph(self, CGraphBackend other, object vertices, const int modus) except -1:
        """
        Use an edge iterator on the subgraph induced by ``vertices`` and do something according to ``modus``.

        INPUT:

        - ``other`` -- a (mutable) subclass of :class:`CGraphBackend`
        - ``vertices`` -- list of vertex labels
        - ``modus`` -- integer representing the modus:
          - ``0`` -- initialize ``other`` to be the subgraph induced by the vertices;
            see :meth:`subgraph_given_vertices`
          - ``1`` -- test whether subgraph of ``self`` induced by the vertices is a subgraph of ``other``
          - ``2`` -- as ``1`` but ignore the labels
        """
        cdef object v, l
        cdef int u_int, v_int, l_int, l_int_other, foo
        cdef CGraph cg = self.cg()
        cdef CGraph cg_other = other.cg()
        cdef list b_vertices_2, all_arc_labels, all_arc_labels_other
        cdef FrozenBitset b_vertices
        cdef int n_vertices = len(vertices)
        cdef bint loops = other.loops()
        cdef bint multiple_edges
        if modus == 0:
            multiple_edges = self.multiple_edges(None) and other.multiple_edges(None)
        elif 1 <= modus <= 2:
            multiple_edges = self.multiple_edges(None)

        if self._directed and not other._directed and modus == 0:
            raise ValueError("cannot obtain an undirected subgraph of a directed graph")

        if self._directed != other._directed and 1 <= modus <= 2:
            if self._directed:
                raise ValueError("cannot check if directed graph is a subgraph of an undirected")
            else:
                raise ValueError("cannot check if undirected graph is a subgraph of a directed")

        b_vertices_2 = [self.get_vertex_checked(v) for v in vertices]
        try:
            b_vertices = FrozenBitset(foo for foo in b_vertices_2 if foo >= 0)
        except ValueError:
            # Avoiding "Bitset must not be empty"
            # in this case there is nothing to do
            return 1

        cdef int length = len(b_vertices)
        cdef int i
        cdef int* vertices_translation = <int *> sig_malloc(b_vertices.capacity() * sizeof(int))

        try:
            # Iterate through the vertices.
            if cg_other.active_vertices.size < length:
                cg_other.realloc(length)
            for j in range(n_vertices):
                i = b_vertices_2[j]
                if i >= 0:
                    v = self.vertex_label(i)
                    if modus == 0:
                        # Add the vertex and obtain the corresponding index.
                        vertices_translation[i] = other.check_labelled_vertex(v, False)
                    elif 1 <= modus <= 2:
                        # Obtain the corresponding index if the vertex is contained in ``other``.
                        foo = other.get_vertex_checked(v)
                        if foo >= 0:
                            vertices_translation[i] = foo
                        else:
                            # Not a subgraph.
                            return 0

            # Iterate through the edges.
            for v_int in b_vertices:
                u_int = cg.next_out_neighbor_unsafe(v_int, -1, &l_int)
                while u_int != -1:
                    if (u_int < b_vertices.capacity() and bitset_in(b_vertices._bitset, u_int)
                            and (u_int >= v_int or other._directed)):
                        # If ``other`` is directed, we should add the arcs in both directions.

                        if modus == 0:
                            # We are adding each arc to ``other``.

                            if unlikely(not loops and u_int == v_int):
                                # Delete loops if ``other`` does not allow loops.
                                u_int = cg.next_out_neighbor_unsafe(v_int, u_int, &l_int)
                                continue

                            if not multiple_edges:
                                if l_int:
                                    l = self.edge_labels[l_int]

                                    # Will return ``0``, if ``other`` does not support edges labels.
                                    l_int_other = other.new_edge_label(l)
                                else:
                                    l_int_other = 0
                                cg_other.add_arc_label_unsafe(vertices_translation[v_int], vertices_translation[u_int], l_int_other)

                            else:
                                all_arc_labels = cg.all_arcs(v_int, u_int)

                                for l_int in all_arc_labels:
                                    if l_int:
                                        l = self.edge_labels[l_int]

                                        # Will return ``0``, if ``other`` does not support edges labels.
                                        l_int_other = other.new_edge_label(l)
                                    else:
                                        l_int_other = 0

                                    cg_other.add_arc_label_unsafe(vertices_translation[v_int], vertices_translation[u_int], l_int_other)

                        else:
                            # Modus is 1 or 2 and we are checking if ``self`` is a subgraph of ``other``.

                            if not multiple_edges:
                                if modus == 1:
                                    l = self.edge_labels[l_int] if l_int else None
                                    if not other._has_labeled_edge_unsafe(vertices_translation[v_int], vertices_translation[u_int], l):
                                        return 0
                                else:
                                    # Ignore the label.
                                    if not cg_other.has_arc_unsafe(vertices_translation[v_int], vertices_translation[u_int]):
                                        return 0

                            else:
                                all_arc_labels = cg.all_arcs(v_int, u_int)

                                if modus == 1:
                                    if len(all_arc_labels) == 1:
                                        l = self.edge_labels[l_int] if l_int else None
                                        if not other._has_labeled_edge_unsafe(vertices_translation[v_int], vertices_translation[u_int], l):
                                            return 0
                                    elif other.multiple_edges(None):
                                        all_arc_labels_other = other._all_edge_labels(vertices_translation[v_int], vertices_translation[u_int])
                                        all_arc_labels = [self.edge_labels[l_int] if l_int else None for l_int in all_arc_labels]
                                        for l in all_arc_labels:
                                            try:
                                                all_arc_labels_other.remove(l)
                                            except ValueError:
                                                return 0

                                    else:
                                        # ``other`` does not allow multiple edges.
                                        # As ``self`` has a multiple edges (not only allows), it cannot be a subgraph.
                                        return 0

                                else:
                                    # Ignore the labels.
                                    if len(all_arc_labels) == 1:
                                        if not cg_other.has_arc_unsafe(vertices_translation[v_int], vertices_translation[u_int]):
                                            return 0
                                    else:
                                        all_arc_labels_other = other._all_edge_labels(vertices_translation[v_int], vertices_translation[u_int])
                                        if len(all_arc_labels) > len(all_arc_labels_other):
                                            return 0

                    u_int = cg.next_out_neighbor_unsafe(v_int, u_int, &l_int)

        finally:
            sig_free(vertices_translation)

        return 1

    ###################################
    # Paths
    ###################################

    def shortest_path_special(self, x, y, exclude_vertices=None, exclude_edges=None, distance_flag=False):
        r"""
        Return the shortest path or distance from ``x`` to ``y``.

        This method is an extension of :meth:`shortest_path` method enabling to
        exclude vertices and/or edges from the search for the shortest path
        between ``x`` and ``y``.

        INPUT:

        - ``x`` -- the starting vertex in the shortest path from ``x`` to ``y``

        - ``y`` -- the end vertex in the shortest path from ``x`` to ``y``

        - ``exclude_vertices`` -- iterable container (default: ``None``);
          iterable of vertices to exclude from the graph while calculating the
          shortest path from ``x`` to ``y``

        - ``exclude_edges`` -- iterable container (default: ``None``); iterable
          of edges to exclude from the graph while calculating the shortest path
          from ``x`` to ``y``

        - ``distance_flag`` -- boolean (default: ``False``); when set to
          ``True``, the shortest path distance from ``x`` to ``y`` is returned
          instead of the path

        OUTPUT:

        - A list of vertices in the shortest path from ``x`` to ``y`` or
          distance from ``x`` to ``y`` is returned depending upon the value of
          parameter ``distance_flag``

        EXAMPLES::

            sage: G = Graph([(1, 2), (2, 3), (3, 4), (1, 5), (5, 6), (6, 7), (7, 4)])
            sage: G._backend.shortest_path_special(1, 4)
            [1, 2, 3, 4]
            sage: G._backend.shortest_path_special(1, 4, exclude_vertices=[5,7])
            [1, 2, 3, 4]
            sage: G._backend.shortest_path_special(1, 4, exclude_vertices=[2, 3])
            [1, 5, 6, 7, 4]
            sage: G._backend.shortest_path_special(1, 4, exclude_vertices=[2], exclude_edges=[(5, 6)])
            []
            sage: G._backend.shortest_path_special(1, 4, exclude_vertices=[2], exclude_edges=[(2, 3)])
            [1, 5, 6, 7, 4]
        """
        cdef bint exclude_v = exclude_vertices
        cdef bint exclude_e = exclude_edges
        cdef bint x_excluded
        cdef bint y_excluded

        if exclude_v:
            x_excluded = x in exclude_vertices
            y_excluded = y in exclude_vertices
            if x_excluded and y_excluded:
                raise LookupError("%s and %s are excluded vertices" % (x, y))
            elif x_excluded:
                raise LookupError("no path from an excluded vertex %s" % (x))
            elif y_excluded:
                raise LookupError("no path to an excluded vertex %s" % (y))
        if x == y:
            if distance_flag:
                return 0
            else:
                return [x]

        # The function being mostly symmetric in x and y, their roles are
        # reversed at the end of each loop. For this reason is defined, for
        # example, two dictionaries dist_y and dist_x containing the distances
        # to x and y, and a dictionary dist_current and dist_other, pointing
        # toward the previous two, alternatively.
        #
        # Besides, there is another difference in the fact that for directed
        # graphs we are interested in paths leaving x toward y, so we are
        # considering the out_neighbors on x's side, and in_neighbors on
        # y's side.

        cdef int x_int = self.get_vertex(x)
        cdef int y_int = self.get_vertex(y)
        cdef int u = 0
        cdef int v = 0
        cdef int w = 0

        cdef set exclude_vertices_int = None
        cdef set exclude_edges_int = None

        if exclude_v:
            exclude_vertices_int = {self.get_vertex(v1) for v1 in exclude_vertices}
        if exclude_e:
            exclude_edges_int = {(self.get_vertex(v1), self.get_vertex(v2)) for v1, v2 in exclude_edges}

        # Each vertex knows its predecessors in the search, for each side
        cdef dict pred_x = {}
        cdef dict pred_y = {}
        cdef dict pred_current = pred_x
        cdef dict pred_other = pred_y

        # Stores the distances from x and y
        cdef dict dist_x = {}
        cdef dict dist_y = {}
        cdef dict dist_current = dist_x
        cdef dict dist_other = dist_y
        dist_x[x_int] = 0
        dist_y[y_int] = 0

        # Lists of vertices whose neighbors have not been explored yet
        cdef list next_x = [x_int]
        cdef list next_y = [y_int]
        cdef list next_current = next_x
        cdef list next_other = next_y
        cdef list next_temporary = []
        cdef list neighbors

        cdef list shortest_path = []

        # We are interested in edges leaving x and entering y, so we
        # are dealing with two different "neighbors" functions
        cdef int out = 1

        # As long as the current side (x or y) is not totally explored ...
        while next_current:
            next_temporary = []

            # Take the next vertex in the list, and study all of its neighbors.
            # When a new neighbor is found, it is added into a temporary list.
            # When all the vertices in the list are tested
            # and next_current is replaced by the temporary list
            #
            # After this, current and other are reversed, and the loop restarts
            for u in next_current:
                if out == 1:
                    nbr = self.cg().out_neighbors(u)
                else:
                    nbr = self.cg().in_neighbors(u)

                if not exclude_e and not exclude_v:
                    neighbors = nbr
                else:
                    neighbors = []
                    for w in nbr:
                        if exclude_v and w in exclude_vertices_int:
                            continue
                        if (exclude_e and
                            ((out == 1 and (u, w) in exclude_edges_int) or
                             (out == -1 and (w, u) in exclude_edges_int))):
                            continue
                        neighbors.append(w)

                for v in neighbors:
                    # If the neighbor is new, updates the distances and adds
                    # to the list.
                    if v not in dist_current:
                        dist_current[v] = dist_current[u] + 1
                        if not distance_flag:
                            pred_current[v] = u
                        next_temporary.append(v)

                        # If the new neighbor is already known by the other
                        # side ...
                        if v in dist_other:
                            # build the shortest path and returns in.
                            if distance_flag:
                                return dist_other[v] + dist_current[v]
                            w = v

                            while w != x_int:
                                shortest_path.append(self.vertex_label(w))
                                w = pred_x[w]

                            shortest_path.append(x)
                            shortest_path.reverse()

                            if v == y_int:
                                return shortest_path

                            w = pred_y[v]
                            while w != y_int:
                                shortest_path.append(self.vertex_label(w))
                                w = pred_y[w]
                            shortest_path.append(y)

                            return shortest_path

            next_current = next_temporary
            pred_current, pred_other = pred_other, pred_current
            dist_current, dist_other = dist_other, dist_current
            next_current, next_other = next_other, next_current
            out = -out

        if distance_flag:
            from sage.rings.infinity import Infinity
            return Infinity
        return []

    def shortest_path(self, x, y, distance_flag=False):
        r"""
        Return the shortest path or distance from ``x`` to ``y``.

        INPUT:

        - ``x`` -- the starting vertex in the shortest path from ``x`` to ``y``

        - ``y`` -- the end vertex in the shortest path from ``x`` to ``y``

        - ``distance_flag`` -- boolean (default: ``False``); when set to
          ``True``, the shortest path distance from ``x`` to ``y`` is returned
          instead of the path

        OUTPUT:

        - A list of vertices in the shortest path from ``x`` to ``y`` or
          distance from ``x`` to ``y`` is returned depending upon the value of
          parameter ``distance_flag``

        EXAMPLES::

            sage: G = Graph(graphs.PetersenGraph())
            sage: G.shortest_path(0, 1)
            [0, 1]
            sage: G.shortest_path_length(0, 1)
            1
        """
        if x == y:
            if distance_flag:
                return 0
            else:
                return [x]

        # The function being mostly symmetric in x and y, their roles are
        # reversed at the end of each loop. For this reason is defined, for
        # example, two dictionaries dist_y and dist_x containing the distances
        # to x and y, and a dictionary dist_current and dist_other, pointing
        # toward the previous two, alternatively.
        #
        # Besides, there is another difference in the fact that for directed
        # graphs we are interested in paths leaving x toward y, so we are
        # considering the out_neighbors on x's side, and in_neighbors on
        # y's side.

        cdef int x_int = self.get_vertex(x)
        cdef int y_int = self.get_vertex(y)
        cdef int u = 0
        cdef int v = 0
        cdef int w = 0

        # Each vertex knows its predecessors in the search, for each side
        cdef dict pred_x = {}
        cdef dict pred_y = {}
        cdef dict pred_current = pred_x
        cdef dict pred_other = pred_y

        # Stores the distances from x and y
        cdef dict dist_x = {}
        cdef dict dist_y = {}
        cdef dict dist_current = dist_x
        cdef dict dist_other = dist_y
        dist_x[x_int] = 0
        dist_y[y_int] = 0

        # Lists of vertices whose neighbors have not been explored yet
        cdef list next_x = [x_int]
        cdef list next_y = [y_int]
        cdef list next_current = next_x
        cdef list next_other = next_y
        cdef list next_temporary = []
        cdef list neighbors

        cdef list shortest_path = []

        # We are interested in edges leaving x and entering y, so we
        # are dealing with two different "neighbors" functions
        cdef int out = 1

        # As long as the current side (x or y) is not totally explored ...
        while next_current:
            next_temporary = []

            # Take the next vertex in the list, and study all of its neighbors.
            # When a new neighbor is found, it is added into a temporary list.
            # When all the vertices in the list are tested
            # and next_current is replaced by the temporary list
            #
            # After this, current and other are reversed, and the loop restarts
            for u in next_current:
                if out == 1:
                    neighbors = self.cg().out_neighbors(u)
                else:
                    neighbors = self.cg().in_neighbors(u)
                for v in neighbors:
                    # If the neighbor is new, updates the distances and adds
                    # to the list.
                    if v not in dist_current:
                        dist_current[v] = dist_current[u] + 1
                        if not distance_flag:
                            pred_current[v] = u
                        next_temporary.append(v)

                        # If the new neighbor is already known by the other
                        # side ...
                        if v in dist_other:
                            # build the shortest path and returns in.
                            if distance_flag:
                                return dist_other[v] + dist_current[v]
                            w = v

                            while w != x_int:
                                shortest_path.append(self.vertex_label(w))
                                w = pred_x[w]

                            shortest_path.append(x)
                            shortest_path.reverse()

                            if v == y_int:
                                return shortest_path

                            w = pred_y[v]
                            while w != y_int:
                                shortest_path.append(self.vertex_label(w))
                                w = pred_y[w]
                            shortest_path.append(y)

                            return shortest_path

            next_current = next_temporary
            pred_current, pred_other = pred_other, pred_current
            dist_current, dist_other = dist_other, dist_current
            next_current, next_other = next_other, next_current
            out = -out

        if distance_flag:
            from sage.rings.infinity import Infinity
            return Infinity
        return []

    def bidirectional_dijkstra_special(self, x, y, weight_function=None,
                                       exclude_vertices=None, exclude_edges=None,
                                       include_vertices=None, distance_flag=False,
                                       reduced_weight=None):
        r"""
        Return the shortest path or distance from ``x`` to ``y`` using a
        bidirectional version of Dijkstra's algorithm.

        This method is an extension of :meth:`bidirectional_dijkstra` method
        enabling to exclude vertices and/or edges from the search for the
        shortest path between ``x`` and ``y``.

        This method also has ``include_vertices`` option enabling to include the
        vertices which will be used to search for the shortest path between
        ``x`` and ``y``.

        INPUT:

        - ``x`` -- the starting vertex in the shortest path from ``x`` to ``y``

        - ``y`` -- the end vertex in the shortest path from ``x`` to ``y``

        - ``exclude_vertices`` -- iterable container (default: ``None``);
          iterable of vertices to exclude from the graph while calculating the
          shortest path from ``x`` to ``y``

        - ``exclude_edges`` -- iterable container (default: ``None``); iterable
          of edges to exclude from the graph while calculating the shortest path
          from ``x`` to ``y``

        - ``include_vertices`` -- iterable container (default: ``None``);
          iterable of vertices to consider in the graph while calculating the
          shortest path from ``x`` to ``y``

        - ``weight_function`` -- function (default: ``None``); a function that
          inputs an edge ``(u, v, l)`` and outputs its weight. If ``None``, we
          use the edge label ``l`` as a weight, if ``l`` is not ``None``, else
          ``1`` as a weight.

        - ``distance_flag`` -- boolean (default: ``False``); when set to
          ``True``, the shortest path distance from ``x`` to ``y`` is returned
          instead of the path.

        - ``reduced_weight`` -- dictionary (default: ``None``); a dictionary
          that takes as input an edge ``(u, v)`` and outputs its reduced weight

        OUTPUT:

        - A list of vertices in the shortest path from ``x`` to ``y`` or
          distance from ``x`` to ``y`` is returned depending upon the value of
          parameter ``distance_flag``

        EXAMPLES::

            sage: G = Graph([(1, 2, 20), (2, 3, 10), (3, 4, 30), (1, 5, 20), (5, 6, 10), (6, 4, 50), (4, 7, 5)])
            sage: G._backend.bidirectional_dijkstra_special(1, 4, weight_function=lambda e:e[2])
            [1, 2, 3, 4]
            sage: G._backend.bidirectional_dijkstra_special(1, 4, weight_function=lambda e:e[2], exclude_vertices=[2], exclude_edges=[(3, 4)])
            [1, 5, 6, 4]
            sage: G._backend.bidirectional_dijkstra_special(1, 4, weight_function=lambda e:e[2], exclude_vertices=[2, 7])
            [1, 5, 6, 4]
            sage: G._backend.bidirectional_dijkstra_special(1, 4, weight_function=lambda e:e[2],  exclude_edges=[(5, 6)])
            [1, 2, 3, 4]
            sage: G._backend.bidirectional_dijkstra_special(1, 4, weight_function=lambda e:e[2],  include_vertices=[1, 5, 6, 4])
            [1, 5, 6, 4]
        """
        cdef bint exclude_v = exclude_vertices
        cdef bint exclude_e = exclude_edges
        cdef bint include_v = include_vertices
        cdef bint x_excluded
        cdef bint y_excluded

        if exclude_v:
            x_excluded = x in exclude_vertices
            y_excluded = y in exclude_vertices
            if x_excluded and y_excluded:
                raise LookupError("%s and %s are excluded vertices" % (x, y))
            elif x_excluded:
                raise LookupError("no path from an excluded vertex %s" % (x))
            elif y_excluded:
                raise LookupError("no path to an excluded vertex %s" % (y))
        if x == y:
            if distance_flag:
                return 0
            else:
                return [x]

        # As for shortest_path, the roles of x and y are symmetric, hence we
        # define dictionaries like pred_current and pred_other, which
        # represent alternatively pred_x or pred_y according to the side
        # studied.
        cdef int x_int = self.get_vertex(x)
        cdef int y_int = self.get_vertex(y)
        cdef int v = 0
        cdef int w = 0
        cdef int pred
        cdef int side
        cdef double distance
        cdef set exclude_vertices_int = None
        cdef set exclude_edges_int = None

        if exclude_v:
            exclude_vertices_int = {self.get_vertex(v1) for v1 in exclude_vertices}
        if exclude_e:
            exclude_edges_int = {(self.get_vertex(v1), self.get_vertex(v2)) for v1, v2 in exclude_edges}
        if include_v:
            include_vertices_int = {self.get_vertex(v1) for v1 in include_vertices}

        # Each vertex knows its predecessors in the search, for each side
        cdef dict pred_x = {}
        cdef dict pred_y = {}
        cdef dict pred_current

        # Stores the distances from x and y
        cdef dict dist_x = {}
        cdef dict dist_y = {}
        cdef dict dist_current
        cdef dict dist_other

        # We use 2 min-heap data structures (pairing heaps), one for the
        # exploration from x and the other for the reverse exploration to y.
        # Each heap associates to a vertex a pair (distance, pred).
        cdef PairingHeap[int, pair[double, int]] px = PairingHeap[int, pair[double, int]]()
        cdef PairingHeap[int, pair[double, int]] py = PairingHeap[int, pair[double, int]]()
        cdef PairingHeap[int, pair[double, int]] * ptmp
        px.push(x_int, (0, x_int))
        py.push(y_int, (0, y_int))

        # Meeting_vertex is a vertex discovered through x and through y
        # which defines the shortest path found
        # (of length shortest_path_length).
        cdef int meeting_vertex = -1
        cdef double shortest_path_length
        cdef double f_tmp

        if reduced_weight is not None:
            def weight_function(e):
                return reduced_weight[(e[0], e[1])]

        # As long as the current side (x or y) is not totally explored ...
        while not (px.empty() and py.empty()):
            if (px.empty() or
                    (not py.empty() and px.top_value().first > py.top_value().first)):
                side = -1
                ptmp = &py
            else:  # px is not empty
                side = 1
                ptmp = &px
            v, (distance, pred) = ptmp.top()
            if meeting_vertex != -1 and distance > shortest_path_length:
                break
            ptmp.pop()

            if side == 1:
                dist_current, dist_other = dist_x, dist_y
                pred_current = pred_x
                nbr_iter = self.cg().out_neighbors(v)
            else:
                dist_current, dist_other = dist_y, dist_x
                pred_current = pred_y
                nbr_iter = self.cg().in_neighbors(v)

            dist_current[v] = distance
            if not distance_flag:
                pred_current[v] = pred

            if v in dist_other:
                f_tmp = distance + dist_other[v]
                if meeting_vertex == -1 or f_tmp < shortest_path_length:
                    meeting_vertex = v
                    shortest_path_length = f_tmp

            if not exclude_e and not exclude_v:
                neighbors = (w for w in nbr_iter
                             if not include_v or w in include_vertices_int)
            else:
                neighbors = (w for w in nbr_iter
                             if ((not exclude_v or w not in exclude_vertices_int) and
                                 (not exclude_e or
                                  ((side == 1 and (v, w) not in exclude_edges_int) or
                                   (side == -1 and (w, v) not in exclude_edges_int))) and
                                 (not include_v or w in include_vertices_int)))

            for w in neighbors:
                # If w has not yet been extracted from the heap, we check if we
                # can improve its path
                if w not in dist_current:
                    v_obj = self.vertex_label(v)
                    w_obj = self.vertex_label(w)
                    if side == -1:
                        v_obj, w_obj = w_obj, v_obj
                    if self._multiple_edges:
                        edge_label = min(weight_function((v_obj, w_obj, l)) for l in self.get_edge_label(v_obj, w_obj))
                    else:
                        edge_label = weight_function((v_obj, w_obj, self.get_edge_label(v_obj, w_obj)))
                    if edge_label < 0:
                        raise ValueError("the graph contains an edge with negative weight")
                    f_tmp = distance + edge_label
                    if ptmp.contains(w):
                        if ptmp.value(w).first > f_tmp:
                            ptmp.decrease(w, (f_tmp, v))
                    else:
                        ptmp.push(w, (f_tmp, v))

        # No meeting point has been found
        if meeting_vertex == -1:
            if distance_flag:
                from sage.rings.infinity import Infinity
                return Infinity
            return []

        if distance_flag:
            if shortest_path_length in ZZ:
                return int(shortest_path_length)
            return shortest_path_length

        # build the shortest path and return it.
        cdef list shortest_path = []
        w = meeting_vertex
        while w != x_int:
            shortest_path.append(self.vertex_label(w))
            w = pred_x[w]

        shortest_path.append(x)
        shortest_path.reverse()

        if meeting_vertex == y_int:
            return shortest_path

        w = pred_y[meeting_vertex]
        while w != y_int:
            shortest_path.append(self.vertex_label(w))
            w = pred_y[w]

        shortest_path.append(y)

        return shortest_path

    def bidirectional_dijkstra(self, x, y, weight_function=None,
                               distance_flag=False):
        r"""
        Return the shortest path or distance from ``x`` to ``y`` using a
        bidirectional version of Dijkstra's algorithm.

        INPUT:

        - ``x`` -- the starting vertex in the shortest path from ``x`` to ``y``

        - ``y`` -- the end vertex in the shortest path from ``x`` to ``y``

        - ``weight_function`` -- function (default: ``None``); a function that
          inputs an edge ``(u, v, l)`` and outputs its weight. If ``None``, we
          use the edge label ``l`` as a weight, if ``l`` is not ``None``, else
          ``1`` as a weight.

        - ``distance_flag`` -- boolean (default: ``False``); when set to
          ``True``, the shortest path distance from ``x`` to ``y`` is returned
          instead of the path.

        OUTPUT:

        - A list of vertices in the shortest path from ``x`` to ``y`` or
          distance from ``x`` to ``y`` is returned depending upon the value of
          parameter ``distance_flag``

        EXAMPLES::

            sage: G = Graph(graphs.PetersenGraph())
            sage: for (u, v) in G.edges(sort=True, labels=None):
            ....:     G.set_edge_label(u, v, 1)
            sage: G.shortest_path(0, 1, by_weight=True)
            [0, 1]
            sage: G.shortest_path_length(0, 1, by_weight=True)
            1
            sage: G = DiGraph([(1, 2, {'weight':1}), (1, 3, {'weight':5}), (2, 3, {'weight':1})])
            sage: G.shortest_path(1, 3, weight_function=lambda e:e[2]['weight'])
            [1, 2, 3]
            sage: G.shortest_path_length(1, 3, weight_function=lambda e:e[2]['weight'])
            2

        TESTS:

        Bugfix from :issue:`7673` ::

            sage: G = Graph([(0, 1, 9), (0, 2, 8), (1, 2, 7)])
            sage: G.shortest_path_length(0, 1, by_weight=True)
            9

        Bugfix from :issue:`28221` ::

            sage: G = Graph([(0, 1, 9.2), (0, 2, 4.5), (1, 2, 4.6)])
            sage: G.shortest_path_length(0, 1, by_weight=True)
            9.1

        Bugfix from :issue:`27464` ::

            sage: G = DiGraph({0: [1, 2], 1: [4], 2: [3, 4], 4: [5], 5: [6]}, multiedges=True)
            sage: for u, v in list(G.edges(labels=None, sort=False)):
            ....:     G.set_edge_label(u, v, 1)
            sage: G.distance(0, 5, by_weight=true)
            3
        """
        if x == y:
            if distance_flag:
                return 0
            else:
                return [x]

        # As for shortest_path, the roles of x and y are symmetric, hence we
        # define dictionaries like pred_current and pred_other, which
        # represent alternatively pred_x or pred_y according to the side
        # studied.
        cdef int x_int = self.get_vertex(x)
        cdef int y_int = self.get_vertex(y)
        cdef int v = 0
        cdef int w = 0
        cdef int pred
        cdef int side
        cdef double distance

        # Each vertex knows its predecessors in the search, for each side
        cdef dict pred_x = {}
        cdef dict pred_y = {}
        cdef dict pred_current

        # Stores the distances from x and y
        cdef dict dist_x = {}
        cdef dict dist_y = {}
        cdef dict dist_current
        cdef dict dist_other

        # We use 2 min-heap data structures (pairing heaps), one for the
        # exploration from x and the other for the reverse exploration to y.
        # Each heap associates to a vertex a pair (distance, pred).
        cdef PairingHeap[int, pair[double, int]] px = PairingHeap[int, pair[double, int]]()
        cdef PairingHeap[int, pair[double, int]] py = PairingHeap[int, pair[double, int]]()
        cdef PairingHeap[int, pair[double, int]] * ptmp
        px.push(x_int, (0, x_int))
        py.push(y_int, (0, y_int))

        cdef list neighbors

        # Meeting_vertex is a vertex discovered through x and through y
        # which defines the shortest path found
        # (of length shortest_path_length).
        cdef int meeting_vertex = -1
        cdef double shortest_path_length
        cdef double f_tmp

        if weight_function is None:
            def weight_function(e):
                return 1 if e[2] is None else e[2]

        # As long as the current side (x or y) is not totally explored ...
        while not (px.empty() and py.empty()):
            if (px.empty() or
                    (not py.empty() and px.top_value().first > py.top_value().first)):
                side = -1
                ptmp = &py
            else:  # px is not empty
                side = 1
                ptmp = &px
            v, (distance, pred) = ptmp.top()
            if meeting_vertex != -1 and distance > shortest_path_length:
                break
            ptmp.pop()

            if side == 1:
                dist_current, dist_other = dist_x, dist_y
                pred_current = pred_x
                neighbors = self.cg().out_neighbors(v)
            else:
                dist_current, dist_other = dist_y, dist_x
                pred_current = pred_y
                neighbors = self.cg().in_neighbors(v)

            dist_current[v] = distance
            if not distance_flag:
                pred_current[v] = pred

            if v in dist_other:
                f_tmp = distance + dist_other[v]
                if meeting_vertex == -1 or f_tmp < shortest_path_length:
                    meeting_vertex = v
                    shortest_path_length = f_tmp

            for w in neighbors:
                # If w has not yet been extracted from the heap, we check if we
                # can improve its path
                if w not in dist_current:
                    v_obj = self.vertex_label(v)
                    w_obj = self.vertex_label(w)
                    if side == -1:
                        v_obj, w_obj = w_obj, v_obj
                    if self._multiple_edges:
                        edge_label = min(weight_function((v_obj, w_obj, l)) for l in self.get_edge_label(v_obj, w_obj))
                    else:
                        edge_label = weight_function((v_obj, w_obj, self.get_edge_label(v_obj, w_obj)))
                    if edge_label < 0:
                        raise ValueError("the graph contains an edge with negative weight")
                    f_tmp = distance + edge_label
                    if ptmp.contains(w):
                        if ptmp.value(w).first > f_tmp:
                            ptmp.decrease(w, (f_tmp, v))
                    else:
                        ptmp.push(w, (f_tmp, v))

        # No meeting point has been found
        if meeting_vertex == -1:
            if distance_flag:
                from sage.rings.infinity import Infinity
                return Infinity
            return []

        if distance_flag:
            if shortest_path_length in ZZ:
                return int(shortest_path_length)
            return shortest_path_length

        # build the shortest path and return it.
        cdef list shortest_path = []
        w = meeting_vertex
        while w != x_int:
            shortest_path.append(self.vertex_label(w))
            w = pred_x[w]

        shortest_path.append(x)
        shortest_path.reverse()

        if meeting_vertex == y_int:
            return shortest_path

        w = pred_y[meeting_vertex]
        while w != y_int:
            shortest_path.append(self.vertex_label(w))
            w = pred_y[w]

        shortest_path.append(y)

        return shortest_path

    def shortest_path_all_vertices(self, v, cutoff=None,
                                   distance_flag=False):
        r"""
        Return for each reachable vertex ``u`` a shortest ``v-u`` path or
        distance from ``v`` to ``u``.

        INPUT:

        - ``v`` -- a starting vertex in the shortest path

        - ``cutoff`` -- integer (default: ``None``); maximal distance of
          returned paths (longer paths will not be returned), ignored when set
          to ``None``

        - ``distance_flag`` -- boolean (default: ``False``); when set to
          ``True``, each vertex ``u`` connected to ``v`` is mapped to shortest
          path distance from ``v`` to ``u`` instead of the shortest path in the
          output dictionary.

        OUTPUT:

        - A dictionary which maps each vertex ``u`` connected to ``v`` to the
          shortest path list or distance from ``v`` to ``u`` depending upon the
          value of parameter ``distance_flag``

        .. NOTE::

            The weight of edges is not taken into account.

        ALGORITHM:

        This is just a breadth-first search.

        EXAMPLES:

        On the Petersen Graph::

            sage: g = graphs.PetersenGraph()
            sage: paths = g._backend.shortest_path_all_vertices(0)
            sage: all((not paths[v] or len(paths[v])-1 == g.distance(0,v)) for v in g)
            True
            sage: g._backend.shortest_path_all_vertices(0, distance_flag=True)
            {0: 0, 1: 1, 2: 2, 3: 2, 4: 1, 5: 1, 6: 2, 7: 2, 8: 2, 9: 2}

        On a disconnected graph::

            sage: g = 2 * graphs.RandomGNP(20, .3)
            sage: paths = g._backend.shortest_path_all_vertices(0)
            sage: all((v not in paths and g.distance(0, v) == +Infinity) or len(paths[v]) - 1 == g.distance(0, v) for v in g)
            True

        TESTS::

            sage: graphs.KrackhardtKiteGraph().eccentricity("a")
            Traceback (most recent call last):
            ...
            LookupError: vertex 'a' is not a vertex of the graph
        """
        cdef list current_layer
        cdef list next_layer
        cdef bitset_t seen
        cdef int v_int
        cdef int u_int
        cdef dict distances
        cdef int d

        distances = {}
        d = 0

        v_int = self.get_vertex(v)
        if v_int == -1:
            raise LookupError(f"vertex {v!r} is not a vertex of the graph")

        bitset_init(seen, self.cg().active_vertices.size)
        bitset_set_first_n(seen, 0)
        bitset_add(seen, v_int)

        current_layer = [(u_int, v_int)
                         for u_int in self.cg().out_neighbors(v_int)]
        next_layer = []

        distances[v] = 0 if distance_flag else [v]

        while current_layer:
            if cutoff is not None and d >= cutoff:
                break

            d += 1
            while current_layer:
                v_int, u_int = current_layer.pop()

                if bitset_not_in(seen, v_int):
                    bitset_add(seen, v_int)
                    if distance_flag:
                        distances[self.vertex_label(v_int)] = d
                    else:
                        distances[self.vertex_label(v_int)] = distances[self.vertex_label(u_int)] + [self.vertex_label(v_int)]
                    next_layer.extend([(u_int, v_int) for u_int in self.cg().out_neighbors(v_int)])

            current_layer = next_layer
            next_layer = []

        bitset_free(seen)
        return distances

    ###################################
    # Searching
    ###################################

    def depth_first_search(self, v, reverse=False, ignore_direction=False,
                           forbidden_vertices=None):
        r"""
        Return a depth-first search from vertex ``v``.

        INPUT:

        - ``v`` -- a vertex from which to start the depth-first search

        - ``reverse`` -- boolean (default: ``False``); this is only relevant to
          digraphs. If this is a digraph, consider the reversed graph in which
          the out-neighbors become the in-neighbors and vice versa.

        - ``ignore_direction`` -- boolean (default: ``False``); this is only
          relevant to digraphs. If this is a digraph, ignore all orientations
          and consider the graph as undirected.

        - ``forbidden_vertices`` -- list (default: ``None``); set of vertices to
          avoid during the search. The start vertex ``v`` cannot be in this set.

        ALGORITHM:

        Below is a general template for depth-first search.

        - **Input:** A directed or undirected graph `G = (V, E)` of order
          `n > 0`. A vertex `s` from which to start the search. The vertices
          are numbered from 1 to `n = |V|`, i.e. `V = \{1, 2, \dots, n\}`.

        - **Output:** A list `D` of distances of all vertices from `s`. A tree
          `T` rooted at `s`.

        #. `S \leftarrow [s]`  # a stack of nodes to visit
        #. `D \leftarrow [\infty, \infty, \dots, \infty]`  # `n` copies of `\infty`
        #. `D[s] \leftarrow 0`
        #. `T \leftarrow [\,]`
        #. while `\text{length}(S) > 0` do

           #. `v \leftarrow \text{pop}(S)`
           #. for each `w \in \text{adj}(v)` do  # for digraphs, use out-neighbor set `\text{oadj}(v)`

              #. if `D[w] = \infty` then

                 #. `D[w] \leftarrow D[v] + 1`
                 #. `\text{push}(S, w)`
                 #. `\text{append}(T, vw)`
        #. return `(D, T)`

        .. SEEALSO::

            - :meth:`breadth_first_search`
              -- breadth-first search for fast compiled graphs.

            - :meth:`breadth_first_search <sage.graphs.generic_graph.GenericGraph.breadth_first_search>`
              -- breadth-first search for generic graphs.

            - :meth:`depth_first_search <sage.graphs.generic_graph.GenericGraph.depth_first_search>`
              -- depth-first search for generic graphs.

        EXAMPLES:

        Traversing the Petersen graph using depth-first search::

            sage: G = graphs.PetersenGraph()
            sage: list(G.depth_first_search(0))
            [0, 5, 8, 6, 9, 7, 2, 3, 4, 1]

        Visiting German cities using depth-first search::

            sage: G = Graph({"Mannheim": ["Frankfurt","Karlsruhe"],
            ....: "Frankfurt": ["Mannheim","Wurzburg","Kassel"],
            ....: "Kassel": ["Frankfurt","Munchen"],
            ....: "Munchen": ["Kassel","Nurnberg","Augsburg"],
            ....: "Augsburg": ["Munchen","Karlsruhe"],
            ....: "Karlsruhe": ["Mannheim","Augsburg"],
            ....: "Wurzburg": ["Frankfurt","Erfurt","Nurnberg"],
            ....: "Nurnberg": ["Wurzburg","Stuttgart","Munchen"],
            ....: "Stuttgart": ["Nurnberg"], "Erfurt": ["Wurzburg"]})
            sage: list(G.depth_first_search("Stuttgart"))
            ['Stuttgart', 'Nurnberg', ...]

        Avoiding some cities:

            sage: list(G.depth_first_search("Stuttgart",
            ....:                 forbidden_vertices=["Frankfurt", "Munchen"]))
            ['Stuttgart', 'Nurnberg', 'Wurzburg', 'Erfurt']

        TESTS:

        The start vertex cannot be forbidden::

            sage: G = graphs.PetersenGraph()
            sage: list(G.depth_first_search(0, forbidden_vertices=[0, 1]))
            Traceback (most recent call last):
            ...
            ValueError: the start vertex is in the set of forbidden vertices
        """
        return Search_iterator(self,
                               v,
                               direction=-1,
                               reverse=reverse,
                               ignore_direction=ignore_direction,
                               forbidden_vertices=forbidden_vertices)

    def breadth_first_search(self, v, reverse=False, ignore_direction=False,
                             report_distance=False, edges=False,
                             forbidden_vertices=None):
        r"""
        Return a breadth-first search from vertex ``v``.

        INPUT:

        - ``v`` -- a vertex from which to start the breadth-first search

        - ``reverse`` -- boolean (default: ``False``); this is only relevant to
          digraphs. If this is a digraph, consider the reversed graph in which
          the out-neighbors become the in-neighbors and vice versa.

        - ``ignore_direction`` -- boolean (default: ``False``); this is only
          relevant to digraphs. If this is a digraph, ignore all orientations
          and consider the graph as undirected.

        - ``report_distance`` -- boolean (default: ``False``); if ``True``,
          reports pairs ``(vertex, distance)`` where ``distance`` is the
          distance from the ``start`` nodes. If ``False`` only the vertices are
          reported.

        - ``edges`` -- boolean (default: ``False``); whether to return the edges
          of the BFS tree in the order of visit or the vertices (default).
          Edges are directed in root to leaf orientation of the tree.

          Note that parameters ``edges`` and ``report_distance`` cannot be
          ``True`` simultaneously.

        - ``forbidden_vertices`` -- list (default: ``None``); set of vertices to
          avoid during the search. The start vertex ``v`` cannot be in this set.

        ALGORITHM:

        Below is a general template for breadth-first search.

        - **Input:** A directed or undirected graph `G = (V, E)` of order
          `n > 0`. A vertex `s` from which to start the search. The vertices
          are numbered from 1 to `n = |V|`, i.e. `V = \{1, 2, \dots, n\}`.

        - **Output:** A list `D` of distances of all vertices from `s`. A
          tree `T` rooted at `s`.

        #. `Q \leftarrow [s]`  # a queue of nodes to visit
        #. `D \leftarrow [\infty, \infty, \dots, \infty]`  # `n` copies of `\infty`
        #. `D[s] \leftarrow 0`
        #. `T \leftarrow [\,]`
        #. while `\text{length}(Q) > 0` do

           #. `v \leftarrow \text{dequeue}(Q)`
           #. for each `w \in \text{adj}(v)` do  # for digraphs, use out-neighbor set `\text{oadj}(v)`

              #. if `D[w] = \infty` then

                 #. `D[w] \leftarrow D[v] + 1`
                 #. `\text{enqueue}(Q, w)`
                 #. `\text{append}(T, vw)`
        #. return `(D, T)`

        .. SEEALSO::

            - :meth:`breadth_first_search <sage.graphs.generic_graph.GenericGraph.breadth_first_search>`
              -- breadth-first search for generic graphs.

            - :meth:`depth_first_search <sage.graphs.generic_graph.GenericGraph.depth_first_search>`
              -- depth-first search for generic graphs.

            - :meth:`depth_first_search`
              -- depth-first search for fast compiled graphs.

        EXAMPLES:

        Breadth-first search of the Petersen graph starting at vertex 0::

            sage: G = Graph(graphs.PetersenGraph())
            sage: list(G.breadth_first_search(0))
            [0, 1, 4, 5, 2, 6, 3, 9, 7, 8]

        Visiting European countries using breadth-first search::

            sage: G = graphs.EuropeMap(continental=True)
            sage: list(G.breadth_first_search("Portugal"))
            ['Portugal', 'Spain', ..., 'Greece']

        Avoiding some countries:

            sage: list(G.breadth_first_search("Portugal",
            ....:                      forbidden_vertices=["Germany","Italy"]))
            ['Portugal', 'Spain', ..., 'Sweden']

        TESTS:

        The start vertex cannot be forbidden::

            sage: G = graphs.PetersenGraph()
            sage: list(G.breadth_first_search(0, forbidden_vertices=[0]))
            Traceback (most recent call last):
            ...
            ValueError: the start vertex is in the set of forbidden vertices
        """
        return Search_iterator(self,
                               v,
                               direction=0,
                               reverse=reverse,
                               ignore_direction=ignore_direction,
                               report_distance=report_distance,
                               edges=edges,
                               forbidden_vertices=forbidden_vertices)

    ###################################
    # Connectedness
    ###################################

    def is_connected(self, forbidden_vertices=None):
        r"""
        Check whether the graph is connected.

        INPUT:

        - ``forbidden_vertices`` -- list (default: ``None``); set of vertices to
          avoid during the search

        EXAMPLES:

        Petersen's graph is connected::

           sage: DiGraph(graphs.PetersenGraph()).is_connected()
           True

        While the disjoint union of two of them is not::

           sage: DiGraph(2*graphs.PetersenGraph()).is_connected()
           False

        A graph with non-integer vertex labels::

            sage: Graph(graphs.CubeGraph(3)).is_connected()
            True

        A graph with forbidden vertices::

            sage: G = graphs.PathGraph(5)
            sage: G._backend.is_connected()
            True
            sage: G._backend.is_connected(forbidden_vertices=[1])
            False
            sage: G._backend.is_connected(forbidden_vertices=[0, 1])
            True

        TESTS::

            sage: P = posets.PentagonPoset()                                            # needs sage.modules
            sage: H = P._hasse_diagram                                                  # needs sage.modules
            sage: H._backend.is_connected()                                             # needs sage.modules
            True
        """
        cdef int v_int
        cdef CGraph cg = self.cg()

        if cg.num_arcs < cg.num_verts - 1:
            return False

        v_int = bitset_first(cg.active_vertices)

        if v_int == -1:
            return True
        v = self.vertex_label(v_int)
        cdef set forbidden = set(forbidden_vertices) if forbidden_vertices else set()
        while v in forbidden:
            v_int = bitset_next(cg.active_vertices, v_int + 1)
            if v_int == -1:
                # The empty graph is connected. So the graph with only forbidden
                # vertices also is
                return True
            v = self.vertex_label(v_int)

        cdef size_t n = len(forbidden)
        for _ in self.depth_first_search(v, ignore_direction=True,
                                         forbidden_vertices=forbidden):
            n += 1
        return n == cg.num_verts

    def is_strongly_connected(self):
        r"""
        Check whether the graph is strongly connected.

        EXAMPLES:

        The circuit on 3 vertices is obviously strongly connected::

            sage: g = DiGraph({0: [1], 1: [2], 2: [0]})
            sage: g.is_strongly_connected()
            True

        But a transitive triangle is not::

            sage: g = DiGraph({0: [1,2], 1: [2]})
            sage: g.is_strongly_connected()
            False
        """
        cdef int v_int = 0
        cdef CGraph cg = self.cg()

        # Pick one vertex
        v_int = bitset_first(cg.active_vertices)

        if v_int == -1:
            return True

        v = self.vertex_label(v_int)

        cdef size_t n = 0
        for _ in self.depth_first_search(v):
            n += 1
        if cg.num_verts != n:
            return False
        n = 0
        for _ in self.depth_first_search(v, reverse=True):
            n += 1
        return cg.num_verts == n

    def strongly_connected_component_containing_vertex(self, v):
        r"""
        Return the strongly connected component containing the given vertex.

        INPUT:

        - ``v`` -- a vertex

        EXAMPLES:

        The digraph obtained from the ``PetersenGraph`` has an unique strongly
        connected component::

            sage: g = DiGraph(graphs.PetersenGraph())
            sage: g.strongly_connected_component_containing_vertex(0)
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

        In the Butterfly DiGraph, each vertex is a strongly connected
        component::

            sage: g = digraphs.ButterflyGraph(3)
            sage: all([v] == g.strongly_connected_component_containing_vertex(v) for v in g)
            True
        """
        cdef set ans = set(self.depth_first_search(v))
        ans.intersection_update(self.depth_first_search(v, reverse=True))
        return list(ans)

    ###################################
    # Miscellaneous
    ###################################

    def is_directed_acyclic(self, certificate=False):
        r"""
        Check whether the graph is both directed and acyclic (possibly with a
        certificate)

        INPUT:

        - ``certificate`` -- boolean (default: ``False``); whether to return a
          certificate

        OUTPUT:

        When ``certificate=False``, returns a boolean value. When
        ``certificate=True`` :

        * If the graph is acyclic, returns a pair ``(True, ordering)`` where
          ``ordering`` is a list of the vertices such that ``u`` appears before
          ``v`` in ``ordering`` if ``u, v`` is an edge.

        * Else, returns a pair ``(False, cycle)`` where ``cycle`` is a list of
          vertices representing a circuit in the graph.

        ALGORITHM:

        We pick a vertex at random, think hard and find out that if we are
        to remove the vertex from the graph we must remove all of its
        out-neighbors in the first place. So we put all of its out-neighbours in
        a stack, and repeat the same procedure with the vertex on top of the
        stack (when a vertex on top of the stack has no out-neighbors, we remove
        it immediately). Of course, for each vertex we only add its outneighbors
        to the end of the stack once : if for some reason the previous algorithm
        leads us to do it twice, it means we have found a circuit.

        We keep track of the vertices whose out-neighborhood has been added to
        the stack once with a variable named ``tried``.

        There is no reason why the graph should be empty at the end of this
        procedure, so we run it again on the remaining vertices until none are
        left or a circuit is found.

        .. NOTE::

            The graph is assumed to be directed. An exception is raised if it is
            not.

        EXAMPLES:

        At first, the following graph is acyclic::

            sage: D = DiGraph({ 0:[1,2,3], 4:[2,5], 1:[8], 2:[7], 3:[7], 5:[6,7], 7:[8], 6:[9], 8:[10], 9:[10] })
            sage: D.plot(layout='circular').show()                                      # needs sage.plot
            sage: D.is_directed_acyclic()
            True

        Adding an edge from `9` to `7` does not change it::

            sage: D.add_edge(9,7)
            sage: D.is_directed_acyclic()
            True

        We can obtain as a proof an ordering of the vertices such that `u`
        appears before `v` if `uv` is an edge of the graph::

            sage: D.is_directed_acyclic(certificate = True)
            (True, [4, 5, 6, 9, 0, 1, 2, 3, 7, 8, 10])

        Adding an edge from 7 to 4, though, makes a difference::

            sage: D.add_edge(7,4)
            sage: D.is_directed_acyclic()
            False

        Indeed, it creates a circuit `7, 4, 5`::

            sage: D.is_directed_acyclic(certificate = True)
            (False, [7, 4, 5])

        Checking acyclic graphs are indeed acyclic ::

            sage: def random_acyclic(n, p):
            ....:  g = graphs.RandomGNP(n, p)
            ....:  h = DiGraph()
            ....:  h.add_edges([ ((u,v) if u<v else (v,u)) for u,v,_ in g.edges(sort=True) ])
            ....:  return h
            ...
            sage: all( random_acyclic(100, .2).is_directed_acyclic()    # long time
            ....:      for i in range(50))
            True

        TESTS::

            sage: m = Matrix(3,[0, 1, 1, 0, 0, 0, 0, 1, 0])                             # needs sage.modules
            sage: g = DiGraph(m)                                                        # needs sage.modules
            sage: g.is_directed_acyclic(certificate=True)                               # needs sage.modules
            (True, [0, 2, 1])
        """
        if not self._directed:
            raise ValueError("Input must be a directed graph.")

        # Activated vertices
        cdef bitset_t activated
        bitset_init(activated, self.cg().active_vertices.size)
        bitset_set_first_n(activated, self.cg().active_vertices.size)

        cdef mp_bitcnt_t uu, u, v

        # Vertices whose neighbors have already been added to the stack
        cdef bitset_t tried
        bitset_init(tried, self.cg().active_vertices.size)
        bitset_set_first_n(tried, 0)

        # Parent of a vertex in the discovery tree
        cdef dict parent = {}

        # The vertices left to be visited
        cdef list stack = []

        # Final ordering, if the graph turns out to be acyclic
        cdef list ordering = []

        # Circuit, if the graph turns out to contain one
        cdef list cycle

        # We try any vertex as the source of the exploration tree
        for v in self.cg().verts():

            # We are not interested in trying de-activated vertices
            if bitset_not_in(activated, v):
                continue

            stack = [v]

            # For as long as some vertices are to be visited
            while stack:

                # We take the last one (depth-first search)
                u = stack[-1]

                # This vertex may have been deactivated since we added it.
                if bitset_not_in(activated, u):
                    stack.pop(-1)
                    continue

                # If we tried this vertex already, it means that all of its
                # out-neighbors have been de-activated already, for we put them
                # *after* u in the stack.
                if bitset_in(tried, u):
                    ordering.append(self.vertex_label(u))
                    bitset_discard(tried, u)
                    bitset_discard(activated, u)
                    stack.pop(-1)
                    continue

                # If we never tried it, now is the time to do it. We also must
                # remember it
                bitset_add(tried, u)

                # We append its out-neighbours to the stack.
                for uu in self.cg().out_neighbors(u):

                    # If we have found a new vertex, we put it at the end of the
                    # stack. We ignored de-activated vertices.
                    if bitset_not_in(tried, uu):
                        if bitset_in(activated, uu):
                            parent[uu] = u
                            stack.append(uu)

                    # If we have already met this vertex, it means we have found
                    # a circuit !
                    else:
                        bitset_free(activated)
                        bitset_free(tried)

                        if not certificate:
                            return False

                        # We build it, then return it
                        # // answer = [u]
                        cycle = [self.vertex_label(u)]

                        while u != uu:
                            u = parent.get(u, uu)
                            cycle.append(self.vertex_label(u))

                        cycle.reverse()
                        return (False, cycle)

        # No Cycle... Good news ! Let's return it.
        bitset_free(activated)
        bitset_free(tried)

        if certificate:
            ordering.reverse()
            return (True, ordering)
        else:
            return True


cdef class Search_iterator:
    r"""
    An iterator for traversing a (di)graph.

    This class is commonly used to perform a depth-first or breadth-first
    search. The class does not build all at once in memory the whole list of
    visited vertices. The class maintains the following variables:

    - ``graph`` -- a graph whose vertices are to be iterated over

    - ``direction`` -- integer; this determines the position at which vertices
      to be visited are removed from the list. For breadth-first search (BFS),
      element removal follow a first-in first-out (FIFO) protocol, as signified
      by the value ``direction=0``. We use a queue to maintain the list of
      vertices to visit in this case. For depth-first search (DFS), element
      removal follow a last-in first-out (LIFO) protocol, as signified by the
      value ``direction=-1``. In this case, we use a stack to maintain the list
      of vertices to visit.

    - ``stack`` -- list of vertices to visit, used only when ``direction=-1``

    - ``queue`` -- a queue of vertices to visit, used only when ``direction=0``

    - ``seen`` -- list of vertices that are already visited

    - ``test_out`` -- boolean; whether we want to consider the out-neighbors
      of the graph to be traversed. For undirected graphs, we consider both
      the in- and out-neighbors. However, for digraphs we only traverse along
      out-neighbors.

    - ``test_in`` -- boolean; whether we want to consider the in-neighbors of
      the graph to be traversed. For undirected graphs, we consider both
      the in- and out-neighbors.

    EXAMPLES::

        sage: g = graphs.PetersenGraph()
        sage: list(g.breadth_first_search(0))
        [0, 1, 4, 5, 2, 6, 3, 9, 7, 8]
    """

    cdef CGraphBackend graph
    cdef int direction
    cdef stack[int] lifo
    cdef queue[int] fifo
    cdef queue[int] fifo_edges
    cdef int first_with_new_distance
    cdef int current_distance
    cdef int n
    cdef bitset_t seen
    cdef bint test_out
    cdef bint test_in
    cdef bint report_distance  # assumed to be constant after initialization
    cdef bint edges            # assumed to be constant after initialization
    cdef in_neighbors

    def __init__(self, graph, v, direction=0, reverse=False,
                 ignore_direction=False, report_distance=False, edges=False,
                 forbidden_vertices=None):
        r"""
        Initialize an iterator for traversing a (di)graph.

        INPUT:

        - ``graph`` -- a graph to be traversed

        - ``v`` -- a vertex in ``graph`` from which to start the traversal

        - ``direction`` -- integer (default: `0`); this determines the
          position at which vertices to be visited are removed from the
          list. For breadth-first search (BFS), element removal follow a
          first-in first-out (FIFO) protocol, as signified by the value
          ``direction=0``. We use a queue to maintain the list of vertices to
          visit in this case. For depth-first search (DFS), element removal
          follow a last-in first-out (LIFO) protocol, as signified by the value
          ``direction=-1``. In this case, we use a stack to maintain the list of
          vertices to visit.

        - ``reverse`` -- boolean (default: ``False``); this is only relevant to
          digraphs. If ``graph`` is a digraph, consider the reversed graph in
          which the out-neighbors become the in-neighbors and vice versa.

        - ``ignore_direction`` -- boolean (default: ``False``); this is only
          relevant to digraphs. If ``graph`` is a digraph, ignore all
          orientations and consider the graph as undirected.

        - ``report_distance`` -- boolean (default: ``False``); if ``True``,
          reports pairs ``(vertex, distance)`` where ``distance`` is the
          distance from the ``start`` nodes. If ``False`` only the vertices are
          reported.
          Only allowed for ``direction=0``, i.e. BFS.

        - ``edges`` -- boolean (default: ``False``); whether to return the edges
          of the BFS tree in the order of visit or the vertices (default).
          Edges are directed in root to leaf orientation of the tree.
          Only allowed for ``direction=0``, i.e. BFS.

          Note that parameters ``edges`` and ``report_distance`` cannot be
          ``True`` simultaneously.

        - ``forbidden_vertices`` -- list (default: ``None``); set of vertices to
          avoid during the search. The start vertex ``v`` cannot be in this set.

        EXAMPLES::

            sage: g = graphs.PetersenGraph()
            sage: list(g.breadth_first_search(0))
            [0, 1, 4, 5, 2, 6, 3, 9, 7, 8]
            sage: list(g.breadth_first_search(0, forbidden_vertices=[1, 2]))
            [0, 4, 5, 3, 9, 7, 8, 6]

        TESTS:

        A vertex which does not belong to the graph::

            sage: list(g.breadth_first_search(-9))
            Traceback (most recent call last):
            ...
            LookupError: vertex (-9) is not a vertex of the graph

        An empty graph::

            sage: list(Graph().breadth_first_search(''))
            Traceback (most recent call last):
            ...
            LookupError: vertex ('') is not a vertex of the graph

        Immutable graphs (see :issue:`16019`)::

            sage: DiGraph([(1, 2)], immutable=True).connected_components(sort=True)
            [[1, 2]]
        """
        self.graph = graph
        self.direction = direction
        if direction != 0 and report_distance:
            raise ValueError("can only report distance for breadth first search")
        self.report_distance = report_distance
        if direction != 0 and edges:
            raise ValueError("can only list edges for breadth first search")
        if report_distance and edges:
            raise ValueError("cannot report distance while returning the edges of the BFS tree")
        self.edges = edges

        self.n = self.graph.cg().active_vertices.size
        bitset_init(self.seen, self.n)
        bitset_set_first_n(self.seen, 0)

        cdef int v_id = self.graph.get_vertex(v)
        cdef int u_id

        if v_id == -1:
            raise LookupError("vertex ({0}) is not a vertex of the graph".format(repr(v)))

        if forbidden_vertices is not None:
            for u in forbidden_vertices:
                u_id = self.graph.get_vertex(u)
                if u_id != -1:
                    if u_id == v_id:
                        raise ValueError(f"the start vertex is in the set of forbidden vertices")
                    bitset_add(self.seen, u_id)

        if direction == 0:
            self.fifo.push(v_id)
            self.first_with_new_distance = -1
            self.current_distance = 0
            bitset_add(self.seen, v_id)
            if self.edges:
                self.fifo_edges.push(-1)
        else:
            self.lifo.push(v_id)

        if not self.graph._directed:
            ignore_direction = False

        self.test_out = (not reverse) or ignore_direction
        self.test_in = reverse or ignore_direction

        if self.test_in:  # How do we list in_neighbors ?
            self.in_neighbors = self.graph.cg().in_neighbors

        if self.edges:
            # The root is not the end of any edge and must therefore be ignored.
            self.next_breadth_first_search()

    def __dealloc__(self):
        r"""
        Freeing the memory
        """
        bitset_free(self.seen)

    def __iter__(self):
        r"""
        Iterate over a traversal of a graph.

        EXAMPLES::

            sage: g = graphs.PetersenGraph()
            sage: g.breadth_first_search(0)
            <generator object ...breadth_first_search at ...
        """
        return self

    cdef inline next_breadth_first_search(self):
        r"""
        Return the next vertex in a breadth first search traversal of a graph.

        EXAMPLES::

            sage: g = graphs.PetersenGraph()
            sage: g.breadth_first_search(0)
            <generator object ...breadth_first_search at ...
            sage: next(g.breadth_first_search(0))
            0
        """
        cdef int v_int
        cdef int w_int
        cdef int l
        cdef CGraph cg = self.graph.cg()

        if not self.fifo.empty():
            v_int = self.fifo.front()
            self.fifo.pop()
            if v_int == self.first_with_new_distance:
                self.current_distance += 1
                self.first_with_new_distance = -1
            value = self.graph.vertex_label(v_int)

            if self.edges:
                prev_int = self.fifo_edges.front()
                self.fifo_edges.pop()
                value_prev = self.graph.vertex_label(prev_int) if prev_int != -1 else None

            if self.test_out:
                w_int = cg.next_out_neighbor_unsafe(v_int, -1, &l)
                while w_int != -1:
                    if bitset_not_in(self.seen, w_int):
                        bitset_add(self.seen, w_int)
                        self.fifo.push(w_int)
                        if self.first_with_new_distance == -1:
                            self.first_with_new_distance = w_int
                        if self.edges:
                            self.fifo_edges.push(v_int)
                    w_int = cg.next_out_neighbor_unsafe(v_int, w_int, &l)
            if self.test_in:
                w_int = cg.next_in_neighbor_unsafe(v_int, -1, &l)
                while w_int != -1:
                    if bitset_not_in(self.seen, w_int):
                        bitset_add(self.seen, w_int)
                        self.fifo.push(w_int)
                        if self.first_with_new_distance == -1:
                            self.first_with_new_distance = w_int
                        if self.edges:
                            self.fifo_edges.push(v_int)
                    w_int = cg.next_in_neighbor_unsafe(v_int, w_int, &l)

        else:
            raise StopIteration

        if self.report_distance:
            return value, smallInteger(self.current_distance)
        elif self.edges:
            return value_prev, value
        return value

    cdef inline next_depth_first_search(self):
        r"""
        Return the next vertex in a depth first search traversal of a graph.

        EXAMPLES::

            sage: g = graphs.PetersenGraph()
            sage: g.depth_first_search(0)
            <generator object ...depth_first_search at ...
            sage: next(g.depth_first_search(0))
            0
        """
        cdef int v_int
        cdef int w_int
        cdef int l
        cdef CGraph cg = self.graph.cg()

        while not self.lifo.empty():
            v_int = self.lifo.top()
            self.lifo.pop()

            if bitset_not_in(self.seen, v_int):
                value = self.graph.vertex_label(v_int)
                bitset_add(self.seen, v_int)

                if self.test_out:
                    w_int = cg.next_out_neighbor_unsafe(v_int, -1, &l)
                    while w_int != -1:
                        self.lifo.push(w_int)
                        w_int = cg.next_out_neighbor_unsafe(v_int, w_int, &l)
                if self.test_in:
                    w_int = cg.next_in_neighbor_unsafe(v_int, -1, &l)
                    while w_int != -1:
                        self.lifo.push(w_int)
                        w_int = cg.next_in_neighbor_unsafe(v_int, w_int, &l)
                break

        else:
            raise StopIteration

        return value

    def __next__(self):
        r"""
        Return the next vertex in a breadth first search traversal of a graph.

        EXAMPLES::

            sage: g = graphs.PetersenGraph()
            sage: g.breadth_first_search(0)
            <generator object ...breadth_first_search at ...
            sage: next(g.breadth_first_search(0))
            0
        """
        if self.direction == 0:
            return self.next_breadth_first_search()
        return self.next_depth_first_search()

##############################
# Functions to simplify edge iterator.
##############################

cdef inline bint _reorganize_edge(object v, object u, const int modus) noexcept:
    """
    Return ``True`` if ``v`` and ``u`` should be exchanged according to the modus.

    INPUT:

    - ``v`` -- vertex

    - ``u`` -- vertex

    - ``modus`` -- integer representing the modus of the iterator:
      - ``0`` -- outgoing edges
      - ``1`` -- ingoing edges
      - ``2`` -- unsorted edges of an undirected graph
      - ``3`` -- sorted edges of an undirected graph

    OUTPUT: boolean according the modus:

    - ``modus == 0`` -- ``False``
    - ``modus == 1`` -- ``True``
    - ``modus == 2`` -- ``True
    - ``modus == 3`` -- ``False if v <= u else True``
    """
    if modus == 0:
        return False
    elif modus == 3:
        try:
            if v <= u:
                return False
        except TypeError:
            pass

    return True
