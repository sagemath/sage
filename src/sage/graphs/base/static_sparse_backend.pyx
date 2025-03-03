r"""
Static sparse graph backend

This module implement a immutable sparse graph backend using the data structure
from :mod:`sage.graphs.base.static_sparse_graph`. It supports both directed and
undirected graphs, as well as vertex/edge labels, loops and multiple edges. As
it uses a very compact C structure it should be very small in memory.

As it is a sparse data structure, you can expect it to be very efficient when
you need to list the graph's edge, or those incident to a vertex, but an
adjacency test can be much longer than in a dense data structure (i.e. like in
:mod:`sage.graphs.base.static_dense_graph`)

For an overview of graph data structures in sage, see
:mod:`~sage.graphs.base.overview`.

Two classes
-----------

This module implements two classes

* :class:`StaticSparseCGraph` extends :class:`~sage.graphs.base.c_graph.CGraph`
  and is a Cython class that manages the definition/deallocation of the
  ``short_digraph`` structure. It does not know anything about labels on
  vertices.

* :class:`StaticSparseBackend` extends
  :class:`~sage.graphs.base.c_graph.CGraphBackend` and is a Python class that
  does know about vertex labels and contains an instance of
  :class:`StaticSparseCGraph` as an internal variable. The input/output of its
  methods are labeled vertices, which it translates to integer id before
  forwarding them to the :class:`StaticSparseCGraph` instance.

Classes and methods
-------------------
"""

from cysignals.memory cimport check_calloc, sig_free

from sage.graphs.base.static_sparse_graph cimport (init_short_digraph,
                                                   init_reverse,
                                                   out_degree,
                                                   has_edge,
                                                   free_short_digraph,
                                                   edge_label)
from sage.graphs.base.c_graph cimport CGraphBackend
from sage.data_structures.bitset cimport FrozenBitset
from libc.stdint cimport uint32_t
from sage.data_structures.bitset_base cimport *

cdef extern from "Python.h":
    int unlikely(int) nogil  # Defined by Cython

cdef class StaticSparseCGraph(CGraph):
    """
    :mod:`CGraph <sage.graphs.base.c_graph>` class based on the sparse graph
    data structure :mod:`static sparse graphs
    <sage.graphs.base.static_sparse_graph>`.
    """

    def __cinit__(self, G, vertex_list=None):
        r"""
        Cython constructor.

        INPUT:

        - ``G`` -- a :class:`Graph` object

        - ``vertex_list`` -- (optional) list of all vertices of ``G``

        The optional argument ``vertex_list`` is assumed to be a list of all
        vertices of the graph ``G`` in some order.
        **Beware that no serious checks are made that this input is correct**.

        If ``vertex_list`` is given, it will be used to map vertices
        of the graph to consecutive integers. Otherwise, the result of
        ``G.vertices(sort=True)`` will be used instead. Because this
        only works if the vertices can be sorted, using
        ``vertex_list`` is useful when working with possibly
        non-sortable objects in Python 3.

        TESTS::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseCGraph
            sage: g = StaticSparseCGraph(graphs.PetersenGraph())

        Check that the digraph methods are working (see :issue:`20253`)::

            sage: G = DiGraph([(0, 1), (1, 0)])
            sage: G2 = G.copy(immutable=True)
            sage: G2.is_strongly_connected()
            True

        Using the ``vertex_list`` optional argument::

            sage: g = StaticSparseCGraph(DiGraph({0: [2]}), vertex_list=[2, 0])
            sage: g.has_arc(0, 1)
            False
            sage: g.has_arc(1, 0)
            True

            sage: g = StaticSparseCGraph(DiGraph({0: [2]}), vertex_list=[2, 0, 4])
            Traceback (most recent call last):
            ...
            ValueError: vertex_list has wrong length
        """
        cdef int i, j, tmp
        has_labels = any(l is not None for _, _, l in G.edge_iterator())
        self._directed = G.is_directed()

        if vertex_list is not None and len(vertex_list) != G.order():
            raise ValueError('vertex_list has wrong length')

        init_short_digraph(self.g, G, edge_labelled=has_labels,
                           vertex_list=vertex_list)
        if self._directed:
            init_reverse(self.g_rev, self.g)

        # Store the number of loops for undirected graphs
        elif not G.has_loops():
            self.number_of_loops = NULL
        else:
            try:
                self.number_of_loops = <int *>check_calloc(self.g.n, sizeof(int))
            except MemoryError:
                free_short_digraph(self.g)
                if self._directed:
                    free_short_digraph(self.g_rev)
                raise
            for i in range(self.g.n):
                for tmp in range(out_degree(self.g, i)):
                    j = self.g.neighbors[i][tmp]
                    if j == i:
                        self.number_of_loops[i] += 1
                    if j > i:
                        break

        # Defining the meaningless set of 'active' vertices. Because of CGraph.
        # As well as num_verts and num_edges
        bitset_init(self.active_vertices,  self.g.n + 1)
        bitset_set_first_n(self.active_vertices, self.g.n)

        self.num_verts = self.g.n
        self.num_arcs = self.g.m

    def __dealloc__(self):
        r"""
        Free the memory.

        TESTS::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseCGraph
            sage: g = StaticSparseCGraph(graphs.PetersenGraph())
        """
        bitset_free(self.active_vertices)
        free_short_digraph(self.g)
        sig_free(self.number_of_loops)
        if self._directed:
            free_short_digraph(self.g_rev)

    cpdef bint has_vertex(self, int v) except -1:
        r"""
        Test if a vertex belongs to the graph

        INPUT:

        - ``n`` -- integer

        TESTS::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseCGraph
            sage: g = StaticSparseCGraph(graphs.PetersenGraph())
            sage: g.has_vertex(1)
            True
            sage: g.has_vertex(10)
            False
        """
        return 0 <= v and v < self.g.n

    cdef int add_vertex_unsafe(self, int v) except -1:
        raise ValueError("graph is immutable; please change a copy instead (use function copy())")

    cdef int del_vertex_unsafe(self, int v) except -1:
        raise ValueError("graph is immutable; please change a copy instead (use function copy())")

    def add_vertex(self, int k):
        r"""
        Add a vertex to the graph. No way.

        TESTS::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseCGraph
            sage: g = StaticSparseCGraph(graphs.PetersenGraph())
            sage: g.add_vertex(45)
            Traceback (most recent call last):
            ...
            ValueError: graph is immutable; please change a copy instead (use function copy())
        """
        self.add_vertex_unsafe(k)

    cpdef del_vertex(self, int k):
        r"""
        Remove a vertex from the graph. No way.

        TESTS::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseCGraph
            sage: g = StaticSparseCGraph(graphs.PetersenGraph())
            sage: g.del_vertex(45)
            Traceback (most recent call last):
            ...
            ValueError: graph is immutable; please change a copy instead (use function copy())
        """
        self.del_vertex_unsafe(k)

    cpdef list verts(self):
        r"""
        Return the list of vertices.

        TESTS::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseCGraph
            sage: g = StaticSparseCGraph(graphs.PetersenGraph())
            sage: g.verts()
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        """
        return list(range(self.g.n))

    cdef int has_arc_label_unsafe(self, int u, int v, int l) except -1:
        """
        Label is ignored.
        """
        return ((0 <= u) and
                (0 <= v) and
                (u < self.g.n) and
                (v < self.g.n) and
                has_edge(self.g, u, v) != NULL)

    cpdef bint has_arc(self, int u, int v) except -1:
        r"""
        Test if `uv` is an edge of the graph

        INPUT:

        - ``u``, ``v`` -- integers

        TESTS::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseCGraph
            sage: g = StaticSparseCGraph(graphs.PetersenGraph())
            sage: g.has_arc(0, 1)
            True
            sage: g.has_arc(0, 7)
            False
        """
        return self.has_arc_unsafe(u, v)

    cdef inline int next_out_neighbor_unsafe(self, int u, int v, int* l) except -2:
        """
        Return the next out-neighbor of ``u`` after ``v``.

        If ``v`` is ``-1`` return the first neighbor of ``u``.

        Return ``-1`` in case there does not exist such an out-neighbor.

        .. NOTE::

            A caller may not alter ``l``.
            It is used to keep track of the current position.
        """
        cdef int degree = out_degree(self.g, u)
        if v == -1:
            l[0] = -1
        for i in range(l[0] + 1, degree):
            if self.g.neighbors[u][i] != v:
                l[0] = i
                return self.g.neighbors[u][i]
        else:
            return -1

    cdef inline int next_in_neighbor_unsafe(self, int u, int v, int* l) except -2:
        """
        Return the next in-neighbor of ``u`` after ``v``.

        If ``v`` is ``-1`` return the first neighbor of ``u``.

        Return ``-1`` in case there does not exist such an in-neighbor.

        .. NOTE::

            A caller may not alter ``l``.
            It is used to keep track of the current position.
        """
        if not self._directed:
            return self.next_out_neighbor_unsafe(u, v, l)
        cdef int degree = out_degree(self.g_rev, u)
        if v == -1:
            l[0] = -1
        for i in range(l[0] + 1, degree):
            if self.g_rev.neighbors[u][i] != v:
                l[0] = i
                return self.g_rev.neighbors[u][i]
        else:
            return -1

    cdef int out_neighbors_unsafe(self, int u, int *neighbors, int size) except -2:
        cdef int degree = self.g.neighbors[u+1] - self.g.neighbors[u]
        cdef int i
        for i in range(min(degree, size)):
            neighbors[i] = self.g.neighbors[u][i]
        return -1 if size < degree else degree

    cdef int in_neighbors_unsafe(self, int u, int *neighbors, int size) except -2:
        if not self._directed:
            return self.out_neighbors_unsafe(u, neighbors, size)

        cdef int degree = self.g_rev.neighbors[u+1] - self.g_rev.neighbors[u]
        cdef int i
        for i in range(min(degree, size)):
            neighbors[i] = self.g_rev.neighbors[u][i]
        return -1 if size < degree else degree

    cpdef list out_neighbors(self, int u):
        r"""
        List the out-neighbors of a vertex.

        INPUT:

        - ``u`` -- a vertex

        TESTS::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseCGraph
            sage: g = StaticSparseCGraph(graphs.PetersenGraph())
            sage: g.out_neighbors(0)
            [1, 4, 5]
            sage: g.out_neighbors(10)
            Traceback (most recent call last):
            ...
            LookupError: the vertex does not belong to the graph
        """
        if u < 0 or u >= self.g.n:
            raise LookupError("the vertex does not belong to the graph")

        cdef int i
        return [<int> self.g.neighbors[u][i] for i in range(out_degree(self.g, u))]

    cpdef list in_neighbors(self, int u):
        r"""
        Return the in-neighbors of a vertex.

        INPUT:

        - ``u`` -- a vertex

        TESTS::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseCGraph
            sage: g = StaticSparseCGraph(graphs.PetersenGraph())
            sage: g.in_neighbors(0)
            [1, 4, 5]
            sage: g.in_neighbors(10)
            Traceback (most recent call last):
            ...
            LookupError: the vertex does not belong to the graph
        """
        if not self._directed:
            return self.out_neighbors(u)

        if u < 0 or u >= self.g.n:
            raise LookupError("the vertex does not belong to the graph")

        cdef int i
        return [<int> self.g_rev.neighbors[u][i] for i in range(out_degree(self.g_rev, u))]

    cpdef int out_degree(self, int u) except -1:
        r"""
        Return the out-degree of a vertex

        INPUT:

        - ``u`` -- a vertex

        TESTS::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseCGraph
            sage: g = StaticSparseCGraph(graphs.PetersenGraph())
            sage: g.out_degree(0)
            3
            sage: g.out_degree(10)
            Traceback (most recent call last):
            ...
            LookupError: the vertex does not belong to the graph
        """
        if u < 0 or u >= self.g.n:
            raise LookupError("the vertex does not belong to the graph")

        return self.g.neighbors[u+1] - self.g.neighbors[u]

    cpdef int in_degree(self, int u) except -1:
        r"""
        Return the in-degree of a vertex

        INPUT:

        - ``u`` -- a vertex

        TESTS::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseCGraph
            sage: g = StaticSparseCGraph(graphs.PetersenGraph())
            sage: g.in_degree(0)
            3
            sage: g.in_degree(10)
            Traceback (most recent call last):
            ...
            LookupError: the vertex does not belong to the graph
        """
        if u < 0 or u >= self.g.n:
            raise LookupError("the vertex does not belong to the graph")

        if not self._directed:
            return self.g.neighbors[u+1] - self.g.neighbors[u]
        else:
            return self.g_rev.neighbors[u+1] - self.g_rev.neighbors[u]


cdef class StaticSparseBackend(CGraphBackend):

    def __init__(self, G, loops=False, multiedges=False, sort=True):
        """
        A graph :mod:`backend <sage.graphs.base.graph_backends>` for static
        sparse graphs.

        EXAMPLES::

            sage: D = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: D.add_edge(0, 1, None, False)
            sage: list(D.iterator_edges(range(9), True))
            [(0, 1, None)]

        ::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            sage: g = StaticSparseBackend(graphs.PetersenGraph())
            sage: list(g.iterator_edges([0], 1))
            [(0, 1, None), (0, 4, None), (0, 5, None)]

        ::

            sage: # needs sage.combinat
            sage: g = DiGraph(digraphs.DeBruijn(4, 3), data_structure='static_sparse')
            sage: gi = DiGraph(g, data_structure='static_sparse')
            sage: gi.edges(sort=True)[0]
            ('000', '000', '0')
            sage: sorted(gi.edges_incident('111'))
            [('111', '110', '0'),
            ('111', '111', '1'),
            ('111', '112', '2'),
            ('111', '113', '3')]

            sage: set(g.edges(sort=False)) == set(gi.edges(sort=False))                 # needs sage.combinat
            True

        ::

            sage: g = graphs.PetersenGraph()
            sage: gi = Graph(g, data_structure='static_sparse')
            sage: g == gi
            True
            sage: set(g.edges(sort=False)) == set(gi.edges(sort=False))
            True

        ::

            sage: gi = Graph({ 0: {1: 1}, 1: {2: 1}, 2: {3: 1}, 3: {4: 2}, 4: {0: 2}}, data_structure='static_sparse')
            sage: (0, 4, 2) in gi.edges(sort=False)
            True
            sage: gi.has_edge(0, 4)
            True

        ::

            sage: G = Graph({1:{2:28, 6:10}, 2:{3:16, 7:14}, 3:{4:12}, 4:{5:22, 7:18}, 5:{6:25, 7:24}})
            sage: GI = Graph({1:{2:28, 6:10}, 2:{3:16, 7:14}, 3:{4:12}, 4:{5:22, 7:18}, 5:{6:25, 7:24}}, data_structure='static_sparse')
            sage: G == GI
            True

        ::

            sage: G = graphs.OddGraph(4)
            sage: d = G.diameter()
            sage: H = G.distance_graph(list(range(d + 1)))
            sage: HI = Graph(H, data_structure='static_sparse')
            sage: HI.size() == len(HI.edges(sort=False))
            True

        ::

            sage: g = Graph({1: {1: [1, 2, 3]}}, data_structure='static_sparse')
            sage: g.size()
            3
            sage: g.order()
            1
            sage: g.vertices(sort=False)
            [1]
            sage: g.edges(sort=True)
            [(1, 1, 1), (1, 1, 2), (1, 1, 3)]

        :issue:`15810` is fixed::

            sage: DiGraph({1: {2: ['a', 'b'], 3: ['c']}, 2: {3: ['d']}}, immutable=True).is_directed_acyclic()
            True
        """
        vertices = list(G)
        if sort:
            try:
                vertices.sort()
            except TypeError:
                pass
        cdef StaticSparseCGraph cg = <StaticSparseCGraph> StaticSparseCGraph(G, vertices)
        self._cg = cg

        self._directed = cg._directed

        self._order = G.order()

        # Does it allow loops/multiedges ?
        self._loops = loops
        self._multiedges = multiedges

        # Dictionary translating a vertex int to a label, and the other way around.
        self._vertex_to_labels = vertices
        self._vertex_to_int = {v: i for i, v in enumerate(vertices)}

        # Needed by CGraph. The first one is just an alias, and the second is
        # useless : accessing _vertex_to_labels (which is a list) is faster than
        # vertex_labels (which is a dictionary)
        self.vertex_ints = self._vertex_to_int
        self.vertex_labels = {i: v for i, v in enumerate(vertices)}
        self._multiple_edges = self._multiedges

    def has_vertex(self, v):
        r"""
        Test if the vertex belongs to the graph.

        INPUT:

        - ``v`` -- a vertex (or not?)

        TESTS::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            sage: g = StaticSparseBackend(graphs.PetersenGraph())
            sage: g.has_vertex(0)
            True
            sage: g.has_vertex("Hey")
            False
        """
        return v in self._vertex_to_int

    def add_vertex(self, object v):
        r"""
        Add a vertex to the graph. No way

        INPUT:

        - ``v`` -- a vertex (or not?)

        TESTS::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            sage: g = StaticSparseBackend(graphs.PetersenGraph())
            sage: g.add_vertex(123)
            Traceback (most recent call last):
            ...
            ValueError: graph is immutable; please change a copy instead (use function copy())
        """
        raise ValueError("graph is immutable; please change a copy instead (use function copy())")

    cpdef add_edge(self, object u, object v, object l, bint directed):
        r"""
        Add an edge to the graph. No way.

        TESTS::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            sage: g = StaticSparseBackend(graphs.PetersenGraph())
            sage: g.add_edge(1,2,3,True)
            Traceback (most recent call last):
            ...
            ValueError: graph is immutable; please change a copy instead (use function copy())
        """
        raise ValueError("graph is immutable; please change a copy instead (use function copy())")

    def add_edges(self, edges, directed, remove_loops=False):
        r"""
        Add edges to the graph. No way.

        INPUT:

        - ``edges`` -- a list of edges (or not?)

        TESTS::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            sage: g = StaticSparseBackend(graphs.PetersenGraph())
            sage: g.add_edges([[1, 2]], True)
            Traceback (most recent call last):
            ...
            ValueError: graph is immutable; please change a copy instead (use function copy())
        """
        raise ValueError("graph is immutable; please change a copy instead (use function copy())")

    def add_vertices(self, vertices):
        r"""
        Add vertices to the graph. No way.

        INPUT:

        - ``vertices`` -- a list of vertices (or not?)

        TESTS::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            sage: g = StaticSparseBackend(graphs.PetersenGraph())
            sage: g.add_vertices([1, 2])
            Traceback (most recent call last):
            ...
            ValueError: graph is immutable; please change a copy instead (use function copy())
        """
        raise ValueError("graph is immutable; please change a copy instead (use function copy())")

    def del_vertex(self, object v):
        r"""
        Delete a vertex from the graph. No way

        INPUT:

        - ``v`` -- a vertex (or not?)

        TESTS::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            sage: g = StaticSparseBackend(graphs.PetersenGraph())
            sage: g.del_vertex(123)
            Traceback (most recent call last):
            ...
            ValueError: graph is immutable; please change a copy instead (use function copy())

        Check that :issue:`39270` is fixed::

            sage: g.del_vertex('a')
            Traceback (most recent call last):
            ...
            ValueError: graph is immutable; please change a copy instead (use function copy())
        """
        raise ValueError("graph is immutable; please change a copy instead (use function copy())")

    def del_vertices(self, vertices):
        r"""
        Delete vertices from the graph. No way

        INPUT:

        - ``vertices`` -- a list of vertices (or not?)

        TESTS::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            sage: g = StaticSparseBackend(graphs.PetersenGraph())
            sage: g.del_vertices([123, 234])
            Traceback (most recent call last):
            ...
            ValueError: graph is immutable; please change a copy instead (use function copy())
        """
        raise ValueError("graph is immutable; please change a copy instead (use function copy())")

    cpdef del_edge(self, object u, object v, object l, bint directed):
        r"""
        Delete an edge of the graph. No way.

        TESTS::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            sage: g = StaticSparseBackend(graphs.PetersenGraph())
            sage: g.del_edge(1,2,3,True)
            Traceback (most recent call last):
            ...
            ValueError: graph is immutable; please change a copy instead (use function copy())
        """
        raise ValueError("graph is immutable; please change a copy instead (use function copy())")

    def del_edges(self, edges, directed):
        r"""
        Delete edges of the graph. No way.

        TESTS::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            sage: g = StaticSparseBackend(graphs.PetersenGraph())
            sage: g.del_edges([[1,2,3]], True)
            Traceback (most recent call last):
            ...
            ValueError: graph is immutable; please change a copy instead (use function copy())
        """
        raise ValueError("graph is immutable; please change a copy instead (use function copy())")

    def set_edge_label(self, u, v, l, directed):
        r"""
        Set edge label. No way.

        TESTS::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            sage: g = StaticSparseBackend(graphs.PetersenGraph())
            sage: g.set_edge_label(1,2,3,True)
            Traceback (most recent call last):
            ...
            ValueError: graph is immutable; please change a copy instead (use function copy())
        """
        raise ValueError("graph is immutable; please change a copy instead (use function copy())")

    def relabel(self, perm, directed):
        r"""
        Relabel the graphs' vertices. No way.

        TESTS::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            sage: g = StaticSparseBackend(graphs.PetersenGraph())
            sage: g.relabel([],True)
            Traceback (most recent call last):
            ...
            ValueError: graph is immutable; please change a copy instead (use function copy())
        """
        raise ValueError("graph is immutable; please change a copy instead (use function copy())")

    def get_edge_label(self, object u, object v):
        """
        Return the edge label for ``(u, v)``.

        INPUT:

        - ``u``, ``v`` -- two vertices

        TESTS::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            sage: g = StaticSparseBackend(graphs.PetersenGraph())
            sage: print(g.get_edge_label(0, 1))
            None
            sage: print(g.get_edge_label(0, "Hey"))
            Traceback (most recent call last):
            ...
            LookupError: one of the two vertices does not belong to the graph
            sage: print(g.get_edge_label(0, 7))
            Traceback (most recent call last):
            ...
            LookupError: the edge does not exist

        ::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            sage: g = StaticSparseBackend(digraphs.DeBruijn(3, 2))                      # needs sage.combinat
            sage: g.has_edge('00', '01', '1')                                           # needs sage.combinat
            True
            sage: g.has_edge('00', '01', '0')                                           # needs sage.combinat
            False
        """
        try:
            u = self._vertex_to_int[u]
            v = self._vertex_to_int[v]
        except KeyError:
            raise LookupError("one of the two vertices does not belong to the graph")

        cdef StaticSparseCGraph cg = self._cg

        cdef uint32_t * edge = has_edge(cg.g, u, v)
        if not edge:
            raise LookupError("the edge does not exist")

        # At this level, edge points toward a edge from u to v in the graph, but
        # not necessarily to the leftmost edge. Hence, we first decrease edge to
        # make it point toward the leftmost such edge, then build the list of
        # all labels.
        if self.multiple_edges(None):
            return self._all_edge_labels(u, v, edge)
        return edge_label(cg.g, edge)

    cdef inline list _all_edge_labels(self, int u, int v, uint32_t* edge=NULL):
        """
        Give the labels of all arcs from ``u`` to ``v``.

        ``u`` and ``v`` are the integers corresponding to vertices.

        ``edge`` may point to an edge from ``u`` to ``v``.
        """
        cdef StaticSparseCGraph cg = self._cg
        if edge is NULL:
            edge = has_edge(cg.g, u, v)

        while edge > cg.g.neighbors[u] and (edge - 1)[0] == v:
            edge -= 1
        cdef list l = []
        while edge < cg.g.neighbors[u+1] and edge[0] == v:
            l.append(edge_label(cg.g, edge))
            edge += 1
        return l

    def has_edge(self, object u, object v, object l):
        """
        Return whether this graph has edge ``(u, v)`` with label ``l``.

        If ``l`` is ``None``, return whether this graph has an edge ``(u, v)``
        with any label.

        INPUT:

        - ``u``, ``v`` -- two vertices

        - ``l`` -- a label

        TESTS::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            sage: g = StaticSparseBackend(graphs.PetersenGraph())
            sage: g.has_edge(0, 1, 'e')
            False
            sage: g.has_edge(0, 4, None)
            True
        """
        try:
            u = self._vertex_to_int[u]
            v = self._vertex_to_int[v]
        except KeyError:
            raise LookupError("one of the two vertices does not belong to the graph")

        return self._has_labeled_edge_unsafe(u, v, l)

    cdef inline bint _has_labeled_edge_unsafe(self, int u, int v, object l) except -1:
        """
        Return whether ``self`` has an arc specified by indices of the vertices
        and an arc label.
        """
        cdef uint32_t * edge = NULL
        cdef StaticSparseCGraph cg = <StaticSparseCGraph> (self._cg)
        edge = has_edge(cg.g, u, v)
        if not edge:
            return False
        if l is None:
            return True

        # At this level, edge points toward a edge from u to v in the graph, but
        # not necessarily toward the right label. As there may be many uv edges
        # with different labels, we first make edge point toward the leftmost uv
        # edge, then scan them all to find the right label.
        while edge > cg.g.neighbors[u] and (edge - 1)[0] == v:
            edge -= 1

        while edge[0] == v and edge < cg.g.neighbors[u+1]:
            if edge_label(cg.g, edge) == l:
                return True
            edge += 1

        return False

    def iterator_in_edges(self, object vertices, bint labels):
        """
        Iterate over the incoming edges incident to a sequence of vertices.

        INPUT:

        - ``vertices`` -- list of vertices

        - ``labels`` -- whether to return labels too

        TESTS::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            sage: g = StaticSparseBackend(graphs.PetersenGraph())
            sage: list(g.iterator_in_edges([0], False))
            [(0, 1), (0, 4), (0, 5)]
            sage: list(g.iterator_in_edges([0], True))
            [(0, 1, None), (0, 4, None), (0, 5, None)]

        ::

            sage: DiGraph(digraphs.Path(5), immutable=False).incoming_edges([2])
            [(1, 2, None)]
            sage: DiGraph(digraphs.Path(5), immutable=True).incoming_edges([2])
            [(1, 2, None)]
        """
        cdef StaticSparseCGraph cg = self._cg
        if not cg._directed:
            for x in self.iterator_out_edges(vertices, labels):
                yield x
            return

        try:
            vertices = [self._vertex_to_int[x] for x in vertices]
        except KeyError:
            raise LookupError("one of the vertices does not belong to the graph")

        cdef int i, j
        for i in vertices:
            vi = self._vertex_to_labels[i]
            for j in range(out_degree(cg.g_rev, i)):
                if labels:
                    yield (self._vertex_to_labels[cg.g_rev.neighbors[i][j]],
                           vi,
                           edge_label(cg.g_rev, cg.g_rev.neighbors[i] + j))
                else:
                    yield self._vertex_to_labels[cg.g_rev.neighbors[i][j]], vi

    def iterator_out_edges(self, object vertices, bint labels):
        """
        Iterate over the outbound edges incident to a sequence of vertices.

        INPUT:

        - ``vertices`` -- list of vertices

        - ``labels`` -- whether to return labels too

        TESTS::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            sage: g = StaticSparseBackend(graphs.PetersenGraph())
            sage: list(g.iterator_out_edges([0], False))
            [(0, 1), (0, 4), (0, 5)]
            sage: list(g.iterator_out_edges([0], True))
            [(0, 1, None), (0, 4, None), (0, 5, None)]
        """
        try:
            vertices = [self._vertex_to_int[x] for x in vertices]
        except KeyError:
            raise LookupError("one of the vertices does not belong to the graph")

        cdef StaticSparseCGraph cg = self._cg
        cdef int i, j
        for i in vertices:
            vi = self._vertex_to_labels[i]
            for j in range(out_degree(cg.g, i)):
                if labels:
                    yield (vi,
                           self._vertex_to_labels[cg.g.neighbors[i][j]],
                           edge_label(cg.g, cg.g.neighbors[i] + j))
                else:
                    yield vi, self._vertex_to_labels[cg.g.neighbors[i][j]]

    def iterator_verts(self, vertices):
        r"""
        Iterate over the vertices.

        INPUT:

        - ``vertices`` -- list of objects; the method will only return the
          elements of the graph which are contained in ``vertices``. It's not
          very efficient. If ``vertices`` is equal to ``None``, all the vertices
          are returned.

        TESTS::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            sage: g = StaticSparseBackend(graphs.PetersenGraph())
            sage: list(g.iterator_verts(None))
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            sage: list(g.iterator_verts([1, "Hey", "I am a french fry"]))
            [1]
        """
        if vertices is None:
            for x in self._vertex_to_labels:
                yield x
            return

        cdef set V = set(vertices)
        for x in self._vertex_to_labels:
            if x in V:
                yield x
        return

    def num_verts(self):
        r"""
        Return the number of vertices.

        TESTS::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            sage: g = StaticSparseBackend(graphs.PetersenGraph())
            sage: g.num_verts()
            10
        """
        return self._order

    def allows_loops(self, value=None):
        r"""
        Return whether the graph allows loops.

        INPUT:

        - ``value`` -- only useful for compatibility with other graph backends,
          where this method can be used to define this boolean. This method
          raises an exception if ``value`` is not equal to ``None``.

        TESTS::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            sage: g = StaticSparseBackend(graphs.PetersenGraph())
            sage: g.allows_loops()
            False
            sage: g = StaticSparseBackend(graphs.PetersenGraph(), loops=True)
            sage: g.allows_loops()
            True
        """
        if value is None:
            return self._loops
        else:
            raise ValueError("the graph is immutable and cannot be changed in any way")

    def multiple_edges(self, value=None):
        r"""
        Return whether the graph allows multiple edges.

        INPUT:

        - ``value`` -- only useful for compatibility with other graph backends,
          where this method can be used to define this boolean. This method
          raises an exception if ``value`` is not equal to ``None``.

        TESTS::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            sage: g = StaticSparseBackend(graphs.PetersenGraph())
            sage: g.multiple_edges()
            False
            sage: g = StaticSparseBackend(graphs.PetersenGraph(), multiedges=True)
            sage: g.multiple_edges()
            True
        """
        if value is None:
            return self._multiedges
        else:
            raise ValueError("the graph is immutable and cannot be changed in any way")

    def num_edges(self, directed):
        r"""
        Return the number of edges.

        INPUT:

        - ``directed`` -- boolean; whether to consider the graph as directed or
          not

        TESTS::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            sage: g = StaticSparseBackend(graphs.PetersenGraph())
            sage: g.num_edges(False)
            15

        Testing the exception::

            sage: g = StaticSparseBackend(digraphs.Circuit(4))
            sage: g.num_edges(False)
            Traceback (most recent call last):
            ...
            NotImplementedError: Sorry, I have no idea what is expected in this situation. I don't think that it is well-defined either, especially for multigraphs.

        :issue:`15491`::

            sage: g = digraphs.RandomDirectedGNP(10, .3)
            sage: gi = DiGraph(g, data_structure='static_sparse')
            sage: gi.size() == len(gi.edges(sort=False))
            True
        """
        cdef StaticSparseCGraph cg = <StaticSparseCGraph> self._cg

        if directed:
            if cg._directed:
                # Returns the real number of directed arcs
                return int(cg.g.m)
            else:
                # Returns twice the number of edges, minus the number of
                # loops. This is actually equal to the index of
                # cg.g.neighbors[cg.g.n] in the array `cg.g.edges`
                return int(cg.g.neighbors[cg.g.n] - cg.g.edges)
        else:
            if cg._directed:
                raise NotImplementedError("Sorry, I have no idea what is expected "
                                          "in this situation. I don't think "
                                          "that it is well-defined either, "
                                          "especially for multigraphs.")
            else:
                # Returns the number of edges
                return int(cg.g.m)

    def iterator_edges(self, vertices, bint labels):
        r"""
        Iterate over the graph's edges.

        INPUT:

        - ``vertices`` -- list; only returns the edges incident to at least one
          vertex of ``vertices``

        - ``labels`` -- boolean; whether to return edge labels too

        TESTS::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            sage: g = StaticSparseBackend(graphs.PetersenGraph())
            sage: list(g.iterator_edges(g.iterator_verts(None), False))
            [(0, 1), (0, 4), (0, 5), (1, 2), (1, 6), (2, 3), (2, 7),
            (3, 4), (3, 8), (4, 9), (5, 7), (5, 8), (6, 8), (6, 9), (7, 9)]

        :issue:`15665`::

            sage: Graph(immutable=True).edges(sort=False)
            []
        """
        cdef FrozenBitset b_vertices

        if not vertices:
            return

        if self._directed:
            raise RuntimeError("this is not meant for directed graphs")

        try:
            vertices = [self._vertex_to_int[x] for x in vertices]
            b_vertices = FrozenBitset(vertices)
        except KeyError:
            raise LookupError("one of the vertices does not belong to the graph")

        cdef StaticSparseCGraph cg = self._cg
        cdef int i, j, tmp

        for i in vertices:
            vi = self._vertex_to_labels[i]
            for tmp in range(out_degree(cg.g, i)):
                j = cg.g.neighbors[i][tmp]
                if j < i and j in b_vertices:
                    continue
                if labels:
                    yield (vi,
                           self._vertex_to_labels[j],
                           edge_label(cg.g, cg.g.neighbors[i] + tmp))
                else:
                    yield vi, self._vertex_to_labels[j]

    iterator_unsorted_edges = iterator_edges

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
        cdef int u_int, prev_u_int, v_int, l_int_other, tmp
        cdef StaticSparseCGraph cg = self._cg
        cdef CGraph cg_other = other.cg()
        cdef list b_vertices_2
        cdef FrozenBitset b_vertices
        cdef int n_vertices = len(vertices)
        cdef bint loops = other.loops()
        cdef bint ignore_multiple_edges = modus == 0 and self.multiple_edges(None) and not other.multiple_edges(None)

        if self._directed and not other._directed and modus == 0:
            raise ValueError("cannot obtain an undirected subgraph of a directed graph")

        if self._directed != other._directed and 1 <= modus <= 2:
            if self._directed:
                raise ValueError("cannot check if directed graph is a subgraph of an undirected")
            else:
                raise ValueError("cannot check if undirected graph is a subgraph of a directed")

        try:
            b_vertices_2 = [self._vertex_to_int[x] for x in vertices]
            b_vertices = FrozenBitset(b_vertices_2)
        except KeyError:
            raise LookupError("one of the vertices does not belong to the graph")
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
                prev_u_int = -1
                for tmp in range(out_degree(cg.g, v_int)):
                    u_int = cg.g.neighbors[v_int][tmp]
                    if (u_int < b_vertices.capacity() and bitset_in(b_vertices._bitset, u_int)
                            and (u_int >= v_int or other._directed)):

                        if unlikely(ignore_multiple_edges and u_int == prev_u_int):
                            # Delete multiple edges, if ``other`` does not allow them.
                            continue

                        if modus == 0:
                            prev_u_int = u_int

                            if unlikely(not loops and u_int == v_int):
                                # Ignore loops, if ``other`` does not allow them.
                                continue

                            l = edge_label(cg.g, cg.g.neighbors[v_int] + tmp)

                            # Will return ``0``, if ``other`` does not support edge labels.
                            l_int_other = other.new_edge_label(l)

                            cg_other.add_arc_label_unsafe(vertices_translation[v_int], vertices_translation[u_int], l_int_other)

                        else:
                            # Modus is 1 or 2.

                            # Check if the arc is contained in ``other``.

                            if unlikely(u_int == prev_u_int):
                                # Check if all of the multiple edges are contained.
                                if not other.multiple_edges(None):
                                    # ``other`` does not allow multiple edges.
                                    # As ``self`` has a multiple edges (not only allows), it cannot be a subgraph.
                                    return 0

                                all_arc_labels = self._all_edge_labels(v_int, u_int)
                                all_arc_labels_other = other._all_edge_labels(vertices_translation[v_int], vertices_translation[u_int])
                                if modus == 2:
                                    # Ignore the labels.
                                    if len(all_arc_labels) > len(all_arc_labels_other):
                                        return 0
                                else:
                                    for l in all_arc_labels:
                                        try:
                                            all_arc_labels_other.remove(l)
                                        except ValueError:
                                            return 0

                                continue
                            prev_u_int = u_int

                            l = edge_label(cg.g, cg.g.neighbors[v_int] + tmp)

                            if modus == 1:
                                if not other._has_labeled_edge_unsafe(vertices_translation[v_int], vertices_translation[u_int], l):
                                    return 0
                            else:
                                # Ignore the label.
                                if not cg_other.has_arc_unsafe(vertices_translation[v_int], vertices_translation[u_int]):
                                    return 0

        finally:
            sig_free(vertices_translation)

        return 1

    def degree(self, v, directed):
        r"""
        Return the degree of a vertex.

        INPUT:

        - ``v`` -- a vertex

        - ``directed`` -- boolean; whether to take into account the orientation
          of this graph in counting the degree of ``v``

        EXAMPLES::

            sage: g = Graph(graphs.PetersenGraph(), data_structure='static_sparse')
            sage: g.degree(0)
            3

        :issue:`17225` about the degree of a vertex with a loop::

            sage: Graph({0: [0]}, immutable=True).degree(0)
            2
            sage: Graph({0: [0], 1: [0, 1, 1, 1]}, immutable=True).degree(1)
            7
        """
        try:
            v = self._vertex_to_int[v]
        except KeyError:
            raise LookupError("the vertex does not belong to the graph")

        cdef StaticSparseCGraph cg = self._cg

        if directed:
            if cg._directed:
                return cg.in_degree(v) + cg.out_degree(v)
            else:
                return 2 * cg.out_degree(v)
        else:
            if cg._directed:
                raise NotImplementedError("Sorry, I have no idea what is expected "
                                          "in this situation. I don't think "
                                          "that it is well-defined either, "
                                          "especially for multigraphs.")
            else:
                return cg.out_degree(v) + (0 if not cg.number_of_loops else cg.number_of_loops[v])

    def in_degree(self, v):
        r"""
        Return the in-degree of a vertex.

        INPUT:

        - ``v`` -- a vertex

        EXAMPLES::

            sage: g = DiGraph(graphs.PetersenGraph(), data_structure='static_sparse')
            sage: g.in_degree(0)
            3
        """
        try:
            v = self._vertex_to_int[v]
        except KeyError:
            raise LookupError("the vertex does not belong to the graph")

        cdef StaticSparseCGraph cg = self._cg

        if cg._directed:
            return cg.in_degree(v)
        else:
            return cg.out_degree(v)

    def out_degree(self, v):
        r"""
        Return the out-degree of a vertex.

        INPUT:

        - ``v`` -- a vertex

        EXAMPLES::

            sage: g = DiGraph(graphs.PetersenGraph(), data_structure='static_sparse')
            sage: g.out_degree(0)
            3
        """
        try:
            v = self._vertex_to_int[v]
        except KeyError:
            raise LookupError("the vertex does not belong to the graph")

        cdef StaticSparseCGraph cg = self._cg

        return cg.out_degree(v)

    def iterator_nbrs(self, v):
        r"""
        Iterate over the neighbors of a vertex.

        INPUT:

        - ``v`` -- a vertex

        EXAMPLES::

            sage: g = Graph(graphs.PetersenGraph(), data_structure='static_sparse')
            sage: g.neighbors(0)
            [1, 4, 5]

        TESTS:

        Issue :issue:`25550` is fixed::

            sage: g = DiGraph({0: [1]}, immutable=True)
            sage: g.neighbors(1)
            [0]
            sage: g = DiGraph({0: [0, 1, 1]}, loops=True, multiedges=True, immutable=True)
            sage: g.neighbors(0)
            [0, 1]
            sage: g = DiGraph({0: [1, 1], 1:[0, 0]}, multiedges=True, immutable=True)
            sage: g.neighbors(0)
            [1]
            sage: g.neighbors(1)
            [0]
        """
        try:
            v = self._vertex_to_int[v]
        except KeyError:
            raise LookupError("the vertex does not belong to the graph")

        cdef StaticSparseCGraph cg = self._cg
        cdef int i, u
        cdef set seen = set()

        if cg._directed:
            for i in range(out_degree(cg.g, v)):
                u = cg.g.neighbors[v][i]
                if u not in seen:
                    yield self._vertex_to_labels[u]
                    seen.add(u)
            for i in range(out_degree(cg.g_rev, v)):
                u = cg.g_rev.neighbors[v][i]
                if u not in seen:
                    yield self._vertex_to_labels[u]
                    seen.add(u)
        else:
            for i in range(out_degree(cg.g, v)):
                u = cg.g.neighbors[v][i]
                if u not in seen:
                    yield self._vertex_to_labels[cg.g.neighbors[v][i]]
                    seen.add(u)

    def iterator_out_nbrs(self, v):
        r"""
        Iterate over the out-neighbors of a vertex.

        INPUT:

        - ``v`` -- a vertex

        EXAMPLES::

            sage: g = DiGraph(graphs.PetersenGraph(), data_structure='static_sparse')
            sage: g.neighbors_out(0)
            [1, 4, 5]
        """
        try:
            v = self._vertex_to_int[v]
        except KeyError:
            raise LookupError("the vertex does not belong to the graph")

        cdef StaticSparseCGraph cg = self._cg
        cdef int i, u
        cdef set seen = set()

        for i in range(out_degree(cg.g, v)):
            u = cg.g.neighbors[v][i]
            if u not in seen:
                yield self._vertex_to_labels[u]
                seen.add(u)

    def iterator_in_nbrs(self, v):
        r"""
        Iterate over the in-neighbors of a vertex.

        INPUT:

        - ``v`` -- a vertex

        EXAMPLES::

            sage: g = DiGraph(graphs.PetersenGraph(), data_structure='static_sparse')
            sage: g.neighbors_in(0)
            [1, 4, 5]

        TESTS::

            sage: g = DiGraph({0: [1]}, immutable=True)
            sage: print(g.neighbors_in(1))
            [0]
        """
        try:
            v = self._vertex_to_int[v]
        except KeyError:
            raise LookupError("the vertex does not belong to the graph")

        cdef StaticSparseCGraph cg = self._cg
        cdef int i, u
        cdef set seen = set()

        if cg._directed:
            for i in range(out_degree(cg.g_rev, v)):
                u = cg.g_rev.neighbors[v][i]
                if u not in seen:
                    yield self._vertex_to_labels[u]
                    seen.add(u)
        else:
            for i in range(out_degree(cg.g, v)):
                u = cg.g.neighbors[v][i]
                if u not in seen:
                    yield self._vertex_to_labels[u]
                    seen.add(u)

    def add_vertex(self, v):
        r"""
        Addition of vertices is not available on an immutable graph.

        EXAMPLES::

            sage: g = DiGraph(graphs.PetersenGraph(), data_structure='static_sparse')
            sage: g.add_vertex(1)
            Traceback (most recent call last):
            ...
            ValueError: graph is immutable; please change a copy instead (use function copy())
            sage: g.add_vertices([1,2,3])
            Traceback (most recent call last):
            ...
            ValueError: graph is immutable; please change a copy instead (use function copy())
        """
        (<StaticSparseCGraph> self._cg).add_vertex(v)

    def del_vertex(self, v):
        r"""
        Removal of vertices is not available on an immutable graph.

        EXAMPLES::

            sage: g = DiGraph(graphs.PetersenGraph(), data_structure='static_sparse')
            sage: g.delete_vertex(1)
            Traceback (most recent call last):
            ...
            ValueError: graph is immutable; please change a copy instead (use function copy())
            sage: g.delete_vertices([1,2,3])
            Traceback (most recent call last):
            ...
            ValueError: graph is immutable; please change a copy instead (use function copy())
        """
        (<StaticSparseCGraph> self._cg).del_vertex(v)


def _run_it_on_static_instead(f):
    r"""
    A decorator function to force the (Di)Graph functions to compute from a
    static sparse graph

    This decorator can be used on methods from (Di)Graph. When it is applied,
    the method that was meant to compute something on a graph first converts
    this graph to a static sparse graph, then does what it had to do on this new
    graph. Of course, it makes no sense to decorate Graph.add_vertex with it as
    such a method will never work on an immutable graph. But it can help find
    new bugs, from time to time.

    EXAMPLES::

        sage: from sage.graphs.base.static_sparse_backend import _run_it_on_static_instead
        sage: @_run_it_on_static_instead
        ....: def new_graph_method(g):
        ....:    print("My backend is of type {}".format(type(g._backend)))
        sage: Graph.new_graph_method = new_graph_method
        sage: g = Graph(5)
        sage: print("My backend is of type {}".format(type(g._backend)))
        My backend is of type <class 'sage.graphs.base.sparse_graph.SparseGraphBackend'>
        sage: g.new_graph_method()
        My backend is of type <class 'sage.graphs.base.static_sparse_backend.StaticSparseBackend'>
    """
    def same_function_on_static_version(*kwd, **kwds):
        if not isinstance(kwd[0]._backend, StaticSparseBackend):
            gcopy = kwd[0].copy(data_structure='static_sparse')
            return getattr(gcopy, f.__name__)(*kwd[1:], **kwds)
        else:
            return f(*kwd, **kwds)

    return same_function_on_static_version
