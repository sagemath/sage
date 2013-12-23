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
from sage.graphs.base.static_sparse_graph cimport (init_short_digraph,
                                                   init_reverse,
                                                   out_degree,
                                                   has_edge,
                                                   free_short_digraph,
                                                   edge_label)
from c_graph import CGraphBackend
from sage.misc.bitset cimport FrozenBitset
from libc.stdint cimport uint32_t

cdef class StaticSparseCGraph(CGraph):
    """
    :mod:`CGraph <sage.graphs.base.c_graph>` class based on the sparse graph
    data structure :mod:`static sparse graphs
    <sage.graphs.base.static_sparse_graph>`.
    """

    def __cinit__(self, G):
        r"""
        Cython constructor

        INPUT:

        - ``G`` -- a :class:`Graph` object.

        TEST::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseCGraph
            sage: g = StaticSparseCGraph(graphs.PetersenGraph())
        """
        has_labels = any(not l is None for _,_,l in G.edge_iterator())
        self.directed = G.is_directed()

        init_short_digraph(self.g,G, edge_labelled = has_labels)
        if self.directed:
            init_reverse(self.g_rev,self.g)

    def __dealloc__(self):
        r"""
        Freeing the memory

        TEST::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseCGraph
            sage: g = StaticSparseCGraph(graphs.PetersenGraph())
        """
        free_short_digraph(self.g)
        if self.g_rev != NULL:
            free_short_digraph(self.g_rev)

    cpdef bint has_vertex(self, int n):
        r"""
        Tests if a vertex belongs to the graph

        INPUT:

        - ``n`` -- an integer

        TEST::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseCGraph
            sage: g = StaticSparseCGraph(graphs.PetersenGraph())
            sage: g.has_vertex(1)
            True
            sage: g.has_vertex(10)
            False
        """
        return 0 <= n and n < self.g.n

    cdef int add_vertex_unsafe(self, int k):
        raise ValueError("Thou shalt not add a vertex to an immutable graph")

    def add_vertex(self, int k):
        r"""
        Adds a vertex to the graph. No way.

        TEST::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseCGraph
            sage: g = StaticSparseCGraph(graphs.PetersenGraph())
            sage: g.add_vertex(45)
            Traceback (most recent call last):
            ...
            ValueError: Thou shalt not add a vertex to an immutable graph

        """
        raise ValueError("Thou shalt not add a vertex to an immutable graph")

    def del_vertex(self, int k):
        r"""
        Removes a vertex from the graph. No way.

        TEST::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseCGraph
            sage: g = StaticSparseCGraph(graphs.PetersenGraph())
            sage: g.del_vertex(45)
            Traceback (most recent call last):
            ...
            ValueError: Thou shalt not remove a vertex from an immutable graph

        """
        raise ValueError("Thou shalt not remove a vertex from an immutable graph")

    cdef int del_vertex_unsafe(self, int v):
        raise ValueError("Thou shalt not remove a vertex from an immutable graph")

    cpdef list verts(self):
        r"""
        Returns the list of vertices

        TEST::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseCGraph
            sage: g = StaticSparseCGraph(graphs.PetersenGraph())
            sage: g.verts()
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        """
        return range(self.g.n)

    cdef int has_arc_unsafe(self, int u, int v):
        return ((0 <= u) and
                (0 <= v) and
                (u < self.g.n) and
                (v < self.g.n) and
                has_edge(self.g, u, v) != NULL)

    cpdef bint has_arc(self, int u, int v):
        r"""
        Tests if uv is an edge of the graph

        INPUT:

        - ``u,v`` -- integers

        TEST::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseCGraph
            sage: g = StaticSparseCGraph(graphs.PetersenGraph())
            sage: g.has_arc(0,1)
            True
            sage: g.has_arc(0,7)
            False
        """
        return self.has_arc_unsafe(u, v)

    cdef int out_neighbors_unsafe(self, int u, int *neighbors, int size) except? -2:
        cdef int degree = self.g.neighbors[u+1] - self.g.neighbors[u]
        cdef int i
        for i in range(min(degree,size)):
            neighbors[i] = self.g.neighbors[u][i]
        return -1 if size < degree else degree

    cdef int in_neighbors_unsafe(self, int u, int *neighbors, int size) except? -2:
        if not self.directed:
            return self.out_neighbors_unsafe(u,neighbors,size)

        cdef int degree = self.g_rev.neighbors[u+1] - self.g_rev.neighbors[u]
        cdef int i
        for i in range(min(degree,size)):
            neighbors[i] = self.g_rev.neighbors[u][i]
        return -1 if size < degree else degree

    cpdef list out_neighbors(self, int u):
        r"""
        List the out-neighbors of a vertex

        INPUT:

        - ``u`` -- a vertex

        TEST::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseCGraph
            sage: g = StaticSparseCGraph(graphs.PetersenGraph())
            sage: g.out_neighbors(0)
            [1, 4, 5]
        """
        if u<0 or u>self.g.n:
            raise LookupError("The vertex does not belong to the graph")

        cdef int i
        return [<int> self.g.neighbors[u][i] for i in range(out_degree(self.g,u))]

    cpdef list in_neighbors(self, int u):
        r"""
        Returns the in-neighbors of a vertex

        INPUT:

        - ``u`` -- a vertex

        TEST::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseCGraph
            sage: g = StaticSparseCGraph(graphs.PetersenGraph())
            sage: g.in_neighbors(0)
            [1, 4, 5]
        """
        if not self.directed:
            return self.out_neighbors(u)

        if u<0 or u>self.g.n:
            raise LookupError("The vertex does not belong to the graph")

        cdef int i
        return [<int> self.g_rev.neighbors[u][i] for i in range(out_degree(self.g_rev,u))]

    cpdef int out_degree(self, int u):
        r"""
        Returns the out-degree of a vertex

        INPUT:

        - ``u`` -- a vertex

        TEST::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseCGraph
            sage: g = StaticSparseCGraph(graphs.PetersenGraph())
            sage: g.out_degree(0)
            3
        """
        if u<0 or u>self.g.n:
            raise LookupError("The vertex does not belong to the graph")

        return self.g.neighbors[u+1] - self.g.neighbors[u]

    cpdef int in_degree(self, int u):
        r"""
        Returns the in-degree of a vertex

        INPUT:

        - ``u`` -- a vertex

        TEST::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseCGraph
            sage: g = StaticSparseCGraph(graphs.PetersenGraph())
            sage: g.in_degree(0)
            3
        """
        if u<0 or u>self.g.n:
            raise LookupError("The vertex does not belong to the graph")

        if not self.directed:
            return self.g.neighbors[u+1] - self.g.neighbors[u]
        else:
            return self.g_rev.neighbors[u+1] - self.g_rev.neighbors[u]

class StaticSparseBackend(CGraphBackend):

    def __init__(self, G, loops = False, multiedges=False):
        """
        A graph :mod:`backend <sage.graphs.base.graph_backends>` for static
        sparse graphs.

        EXAMPLE::

            sage: D = sage.graphs.base.sparse_graph.SparseGraphBackend(9)
            sage: D.add_edge(0,1,None,False)
            sage: list(D.iterator_edges(range(9), True))
            [(0, 1, None)]

        ::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            sage: g = StaticSparseBackend(graphs.PetersenGraph())
            sage: list(g.iterator_edges([0],1))
            [(0, 1, None), (0, 4, None), (0, 5, None)]

        ::

            sage: g=DiGraph(digraphs.DeBruijn(4,3),data_structure="static_sparse")
            sage: gi=DiGraph(g,data_structure="static_sparse")
            sage: gi.edges()[0]
            ('000', '000', '0')
            sage: gi.edges_incident('111')
            [('111', '110', '0'), ('111', '111', '1'), ('111', '112', '2'), ('111', '113', '3')]
            sage: sorted(g.edges()) == sorted(gi.edges())
            True

        ::

            sage: g = graphs.PetersenGraph()
            sage: gi=Graph(g,data_structure="static_sparse")
            sage: g == gi
            True
            sage: sorted(g.edges()) == sorted(gi.edges())
            True

        ::

            sage: gi = Graph( { 0: {1: 1}, 1: {2: 1}, 2: {3: 1}, 3: {4: 2}, 4: {0: 2} }, data_structure="static_sparse")
            sage: (0,4,2) in gi.edges()
            True
            sage: gi.has_edge(0,4)
            True

        ::

            sage: G = Graph({1:{2:28, 6:10}, 2:{3:16, 7:14}, 3:{4:12}, 4:{5:22, 7:18}, 5:{6:25, 7:24}})
            sage: GI = Graph({1:{2:28, 6:10}, 2:{3:16, 7:14}, 3:{4:12}, 4:{5:22, 7:18}, 5:{6:25, 7:24}}, data_structure="static_sparse")
            sage: G == GI
            True

        ::

            sage: G = graphs.OddGraph(4)
            sage: d = G.diameter()
            sage: H = G.distance_graph(range(d+1))
            sage: HI = Graph(H,data_structure="static_sparse")
            sage: HI.size() == len(HI.edges())
            True

        ::

            sage: g = Graph({1:{1:[1,2,3]}}, data_structure="static_sparse")
            sage: g.size()
            3
            sage: g.order()
            1
            sage: g.vertices()
            [1]
            sage: g.edges()
            [(1, 1, 1), (1, 1, 2), (1, 1, 3)]
        """
        cdef StaticSparseCGraph cg = <StaticSparseCGraph> StaticSparseCGraph(G)
        self._cg = cg
        self._directed = cg.directed

        vertices = G.vertices()
        self._order = len(vertices)

        # Does it allow loops/multiedges ?
        self._loops = loops
        self._multiedges = multiedges

        # Dictionary translating a vertex int to a label, and the other way around.
        self._vertex_to_labels = vertices
        self._vertex_to_int = {v:i for i,v in enumerate(vertices)}

    def has_vertex(self, v):
        r"""
        Tests if the vertex belongs to the graph

        INPUT:

        - ``v`` -- a vertex (or not?)

        TEST::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            sage: g = StaticSparseBackend(graphs.PetersenGraph())
            sage: g.has_vertex(0)
            True
            sage: g.has_vertex("Hey")
            False
        """
        return v in self._vertex_to_int

    def get_edge_label(self, object u, object v):
        """
        Returns the edge label for ``(u,v)``.

        INPUT:

        - ``u,v`` -- two vertices

        TEST::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            sage: g = StaticSparseBackend(graphs.PetersenGraph())
            sage: print g.get_edge_label(0,1)
            None
            sage: print g.get_edge_label(0,"Hey")
            Traceback (most recent call last):
            ...
            LookupError: One of the two vertices does not belong to the graph
            sage: print g.get_edge_label(0,7)
            Traceback (most recent call last):
            ...
            LookupError: The edge does not exist

        ::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            sage: g = StaticSparseBackend(digraphs.DeBruijn(3,2))
            sage: g.has_edge('00','01','1')
            True
            sage: g.has_edge('00','01','0')
            False
        """
        try:
            u = self._vertex_to_int[u]
            v = self._vertex_to_int[v]
        except KeyError:
            raise LookupError("One of the two vertices does not belong to the graph")

        cdef StaticSparseCGraph cg = self._cg
        cdef list l

        cdef uint32_t * edge = has_edge(cg.g,u,v)
        if edge == NULL:
            raise LookupError("The edge does not exist")

        # At this level, edge points toward a edge from u to v in the graph, but
        # not necessarily to the leftmost edge. Hence, we first decrease edge to
        # make it point toward the leftmost such edge, then build the list of
        # all labels.
        if self.multiple_edges(None):
            while edge > cg.g.neighbors[u] and (edge-1)[0] == v:
                edge -= 1
            l = []
            while edge < cg.g.neighbors[u+1] and edge[0] == v:
                l.append(edge_label(cg.g,edge))
                edge += 1
            return l

        else:
            return edge_label(cg.g,edge)

    def has_edge(self, object u, object v, object l):
        """
        Returns whether this graph has edge ``(u,v)`` with label ``l``.

        If ``l`` is ``None``, return whether this graph has an edge ``(u,v)``
        with any label.

        INPUT:

        - ``u,v`` -- two vertices

        - ``l`` -- a label

        TEST::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            sage: g = StaticSparseBackend(graphs.PetersenGraph())
            sage: g.has_edge(0,1,'e')
            False
            sage: g.has_edge(0,4,None)
            True
        """
        cdef uint32_t * edge = NULL
        cdef StaticSparseCGraph cg = <StaticSparseCGraph> (self._cg)
        try:
            u = self._vertex_to_int[u]
            v = self._vertex_to_int[v]
        except KeyError:
            raise LookupError("One of the two vertices does not belong to the graph")

        edge = has_edge(cg.g,u,v)
        if edge == NULL:
            return False
        if l is None:
            return True

        # At this level, edge points toward a edge from u to v in the graph, but
        # not necessarily toward the right label. As there may be many uv edges
        # with different labels, we first make edge point toward the leftmost uv
        # edge, then scan them all to find the right label.
        while edge > cg.g.neighbors[u] and (edge-1)[0] == v :
            edge -= 1

        while edge[0] == v and edge < cg.g.neighbors[u+1]:
            if edge_label(cg.g,edge) == l:
                return True
            edge += 1

        return False

    def iterator_in_edges(self, object vertices, bint labels):
        """
        Iterate over the incoming edges incident to a sequence of vertices.

        INPUT:

        - ``vertices`` -- a list of vertices

        - ``labels`` -- whether to return labels too

        TEST::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            sage: g = StaticSparseBackend(graphs.PetersenGraph())
            sage: list(g.iterator_in_edges([0],False))
            [(0, 1), (0, 4), (0, 5)]
            sage: list(g.iterator_in_edges([0],True))
            [(0, 1, None), (0, 4, None), (0, 5, None)]
        """
        cdef StaticSparseCGraph cg = self._cg
        if not cg.directed:
            for x in self.iterator_out_edges(vertices, labels):
                yield x
            return

        try:
            vertices = [self._vertex_to_int[x] for x in vertices]
        except KeyError:
            raise LookupError("One of the vertices does not belong to the graph")

        cdef int i,j
        for i in vertices:
            vi = self._vertex_to_labels[i]
            for j in range(out_degree(cg.g_rev,i)):
                if labels:
                    yield (vi,
                           self._vertex_to_labels[cg.g_rev.neighbors[i][j]],
                           edge_label(cg.g_rev,cg.g_rev.neighbors[i]+j))
                else:
                    yield vi,self._vertex_to_labels[cg.g_rev.neighbors[i][j]]

    def iterator_out_edges(self, object vertices, bint labels):
        """
        Iterate over the outbound edges incident to a sequence of vertices.

        INPUT:

        - ``vertices`` -- a list of vertices

        - ``labels`` -- whether to return labels too

        TEST::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            sage: g = StaticSparseBackend(graphs.PetersenGraph())
            sage: list(g.iterator_out_edges([0], False))
            [(0, 1), (0, 4), (0, 5)]
            sage: list(g.iterator_out_edges([0],True))
            [(0, 1, None), (0, 4, None), (0, 5, None)]
        """
        try:
            vertices = [self._vertex_to_int[x] for x in vertices]
        except KeyError:
            raise LookupError("One of the vertices does not belong to the graph")

        cdef StaticSparseCGraph cg = self._cg
        cdef int i,j
        for i in vertices:
            vi = self._vertex_to_labels[i]
            for j in range(out_degree(cg.g,i)):
                if labels:
                    yield (vi,
                           self._vertex_to_labels[cg.g.neighbors[i][j]],
                           edge_label(cg.g,cg.g.neighbors[i]+j))
                else:
                    yield vi,self._vertex_to_labels[cg.g.neighbors[i][j]]

    def iterator_verts(self, vertices):
        r"""
        Returns an iterator over the vertices

        INPUT:

        - ``vertices`` -- a list of objects. The method will only return the
          elements of the graph which are contained in ``vertices``. It's not
          very efficient. If ``vertices`` is equal to ``None``, all the vertices
          are returned.

        TEST::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            sage: g = StaticSparseBackend(graphs.PetersenGraph())
            sage: list(g.iterator_verts(None))
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            sage: list(g.iterator_verts([1,"Hey","I am a french fry"]))
            [1]
        """
        if vertices is None:
            return iter(self._vertex_to_labels)
        else:
            return (x for x in self._vertex_to_labels if x in vertices)

    def num_verts(self):
        r"""
        Returns the number of vertices

        TEST::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            sage: g = StaticSparseBackend(graphs.PetersenGraph())
            sage: g.num_verts()
            10
        """
        return self._order

    def allows_loops(self, value=None):
        r"""
        Returns whether the graph allows loops

        INPUT:

        - ``value`` -- only useful for compatibility with other graph backends,
          where this method can be used to define this boolean. This method
          raises an exception if ``value`` is not equal to ``None``.

        TEST::

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
            raise ValueError("The graph is immutable. You cannot change it in any way !")

    def multiple_edges(self, value=None):
        r"""
        Returns whether the graph allows multiple edges

        INPUT:

        - ``value`` -- only useful for compatibility with other graph backends,
          where this method can be used to define this boolean. This method
          raises an exception if ``value`` is not equal to ``None``.

        TEST::

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
            raise ValueError("The graph is immutable. You cannot change it in any way !")

    def num_edges(self,directed):
        r"""
        Returns the number of edges

        INPUT:

        - ``directed`` (boolean) -- whether to consider the graph as directed or
          not.

        TEST::

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

        :trac:`15491`::

            sage: g=digraphs.RandomDirectedGNP(10,.3)
            sage: gi=DiGraph(g,data_structure="static_sparse")
            sage: gi.size() == len(gi.edges())
            True
        """
        cdef StaticSparseCGraph cg = <StaticSparseCGraph> self._cg

        if directed:
            if cg.directed:
                # Returns the real number of directed arcs
                return int(cg.g.m)
            else:
                # Returns twice the number of edges, minus the number of
                # loops. This is actually equal to the index of
                # cg.g.neighbors[cg.g.n] in the array `cg.g.edges`
                return int(cg.g.neighbors[cg.g.n]-cg.g.edges)
        else:
            if cg.directed:
                raise NotImplementedError("Sorry, I have no idea what is expected "
                                          "in this situation. I don't think "
                                          "that it is well-defined either, "
                                          "especially for multigraphs.")
            else:
                # Returns the number of edges
                return int(cg.g.m)

    def iterator_edges(self, vertices, bint labels):
        r"""
        Returns an iterator over the graph's edges.

        INPUT:

        - ``vertices`` -- only returns the edges incident to at least one vertex
          of ``vertices``.

        - ``labels`` -- whether to return edge labels too

        TEST::

            sage: from sage.graphs.base.static_sparse_backend import StaticSparseBackend
            sage: g = StaticSparseBackend(graphs.PetersenGraph())
            sage: list(g.iterator_edges(g.iterator_verts(None), False))
            [(0, 1), (0, 4), (0, 5), (1, 2), (1, 6), (2, 3), (2, 7),
            (3, 4), (3, 8), (4, 9), (5, 7), (5, 8), (6, 8), (6, 9), (7, 9)]
        """
        cdef FrozenBitset fb
        try:
            vertices = FrozenBitset([self._vertex_to_int[x] for x in vertices])
        except KeyError:
            raise LookupError("One of the vertices does not belong to the graph")

        cdef StaticSparseCGraph cg = self._cg
        cdef int i,j,tmp

        for i in vertices:
            vi = self._vertex_to_labels[i]
            for tmp in range(out_degree(cg.g,i)):
                j = cg.g.neighbors[i][tmp]
                if j < i and j in vertices:
                    continue
                if labels:
                    yield (vi,
                           self._vertex_to_labels[j],
                           edge_label(cg.g,cg.g.neighbors[i]+tmp))
                else:
                    yield vi,self._vertex_to_labels[j]

    def degree(self, v, directed):
        r"""
        Returns the degree of a vertex

        INPUT:

        - ``v`` -- a vertex

        - ``directed`` -- boolean; whether to take into account the
          orientation of this graph in counting the degree of ``v``.

        EXAMPLE::

            sage: g = Graph(graphs.PetersenGraph(), data_structure="static_sparse")
            sage: g.degree(0)
            3
        """
        try:
            v = self._vertex_to_int[v]
        except KeyError:
            raise LookupError("The vertex does not belong to the graph")

        cdef StaticSparseCGraph cg = self._cg

        if directed:
            if cg.directed:
                return cg.in_degree(v) + cg.out_degree(v)
            else:
                return 2*cg.out_degree(v)
        else:
            if cg.directed:
                raise NotImplementedError("Sorry, I have no idea what is expected "
                                          "in this situation. I don't think "
                                          "that it is well-defined either, "
                                          "especially for multigraphs.")
            else:
                return cg.out_degree(v)

    def in_degree(self, v):
        r"""
        Returns the in-degree of a vertex

        INPUT:

        - ``v`` -- a vertex

        EXAMPLE::

            sage: g = DiGraph(graphs.PetersenGraph(), data_structure="static_sparse")
            sage: g.in_degree(0)
            3
        """
        try:
            v = self._vertex_to_int[v]
        except KeyError:
            raise LookupError("The vertex does not belong to the graph")

        cdef StaticSparseCGraph cg = self._cg

        if cg.directed:
            return cg.in_degree(v)
        else:
            return cg.out_degree(v)

    def out_degree(self, v):
        r"""
        Returns the out-degree of a vertex

        INPUT:

        - ``v`` -- a vertex

        EXAMPLE::

            sage: g = DiGraph(graphs.PetersenGraph(), data_structure="static_sparse")
            sage: g.out_degree(0)
            3
        """
        try:
            v = self._vertex_to_int[v]
        except KeyError:
            raise LookupError("The vertex does not belong to the graph")

        cdef StaticSparseCGraph cg = self._cg

        return cg.out_degree(v)

    def iterator_nbrs(self, v):
        r"""
        Returns the neighbors of a vertex

        INPUT:

        - ``v`` -- a vertex

        EXAMPLE::

            sage: g = Graph(graphs.PetersenGraph(), data_structure="static_sparse")
            sage: g.neighbors(0)
            [1, 4, 5]
        """
        try:
            v = self._vertex_to_int[v]
        except KeyError:
            raise LookupError("The vertex does not belong to the graph")

        cdef StaticSparseCGraph cg = self._cg
        cdef int i

        for i in range(out_degree(cg.g,v)):
            yield self._vertex_to_labels[cg.g.neighbors[v][i]]

    def iterator_out_nbrs(self, v):
        r"""
        Returns the out-neighbors of a vertex

        INPUT:

        - ``v`` -- a vertex

        EXAMPLE::

            sage: g = DiGraph(graphs.PetersenGraph(), data_structure="static_sparse")
            sage: g.neighbors_out(0)
            [1, 4, 5]
        """
        try:
            v = self._vertex_to_int[v]
        except KeyError:
            raise LookupError("The vertex does not belong to the graph")

        cdef StaticSparseCGraph cg = self._cg
        cdef int i

        for i in range(out_degree(cg.g,v)):
            yield self._vertex_to_labels[cg.g.neighbors[v][i]]

    def iterator_in_nbrs(self, v):
        r"""
        Returns the out-neighbors of a vertex

        INPUT:

        - ``v`` -- a vertex

        EXAMPLE::

            sage: g = DiGraph(graphs.PetersenGraph(), data_structure="static_sparse")
            sage: g.neighbors_in(0)
            [1, 4, 5]
        """
        try:
            v = self._vertex_to_int[v]
        except KeyError:
            raise LookupError("The vertex does not belong to the graph")

        cdef StaticSparseCGraph cg = self._cg
        cdef short_digraph g

        if cg.directed:
            for i in range(out_degree(cg.g_rev,v)):
                yield self._vertex_to_labels[cg.g_rev.neighbors[v][i]]
        else:
            for i in range(out_degree(cg.g,v)):
                yield self._vertex_to_labels[cg.g.neighbors[v][i]]

def _run_it_on_static_instead(f):
    r"""
    A decorator function to force the (Di)Graph functions to compute from a
    static sparse graph3

    This decorator can be used on methods from (Di)Graph. When it is applied,
    the method that was meant to compute something on a graph first converts
    this graph to a static sparse graph, then does what it had to do on this new
    graph. Of course, it makes no sense to decorate Graph.add_vertex with it as
    such a method will never work on an immutable graph. But it can help find
    new bugs, from time to time.

    EXAMPLE::

        sage: from sage.graphs.base.static_sparse_backend import _run_it_on_static_instead
        sage: @_run_it_on_static_instead
        ....: def new_graph_method(g):
        ....:    print "My backend is of type", g._backend
        sage: Graph.new_graph_method = new_graph_method
        sage: g = Graph(5)
        sage: print "My backend is of type", g._backend
        My backend is of type <class 'sage.graphs.base.sparse_graph.SparseGraphBackend'>
        sage: g.new_graph_method()
        My backend is of type <class 'sage.graphs.base.static_sparse_backend.StaticSparseBackend'>
    """
    def same_function_on_static_version(*kwd,**kwds):
        if not isinstance(kwd[0]._backend,StaticSparseBackend):
            gcopy = kwd[0].copy(data_structure="static_sparse")
            return getattr(gcopy,f.__name__)(*kwd[1:],**kwds)
        else:
            return f(*kwd,**kwds)

    return same_function_on_static_version
