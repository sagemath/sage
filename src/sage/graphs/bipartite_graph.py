r"""
Bipartite graphs

This module implements bipartite graphs.

AUTHORS:

- Robert L. Miller (2008-01-20): initial version

- Ryan W. Hinton (2010-03-04): overrides for adding and deleting vertices
  and edges

TESTS::

    sage: B = graphs.CompleteBipartiteGraph(7, 9)
    sage: loads(dumps(B)) == B
    True

::

    sage: B = BipartiteGraph(graphs.CycleGraph(4))
    sage: B == B.copy()
    True
    sage: type(B.copy())
    <class 'sage.graphs.bipartite_graph.BipartiteGraph'>
"""

#*****************************************************************************
#         Copyright (C) 2008 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function
from __future__ import absolute_import
from six.moves import range

from .graph import Graph

class BipartiteGraph(Graph):
    r"""
    Bipartite graph.

    INPUT:

    - ``data`` -- can be any of the following:

      #. Empty or ``None`` (creates an empty graph).
      #. An arbitrary graph.
      #. A reduced adjacency matrix.
      #. A file in alist format.
      #. From a NetworkX bipartite graph.

    A reduced adjacency matrix contains only the non-redundant portion of the
    full adjacency matrix for the bipartite graph.  Specifically, for zero
    matrices of the appropriate size, for the reduced adjacency matrix ``H``,
    the full adjacency matrix is ``[[0, H'], [H, 0]]``.

    The alist file format is described at
    http://www.inference.phy.cam.ac.uk/mackay/codes/alist.html

    - ``partition`` -- (default: ``None``) a tuple defining vertices of the left and right
      partition of the graph. Partitions will be determined automatically
      if ``partition``=``None``.

    - ``check`` -- (default: ``True``) if ``True``, an invalid input partition
      raises an exception. In the other case offending edges simply won't
      be included.

    .. NOTE::

        All remaining arguments are passed to the ``Graph`` constructor

    EXAMPLES:

    1. No inputs or ``None`` for the input creates an empty graph::

        sage: B = BipartiteGraph()
        sage: type(B)
        <class 'sage.graphs.bipartite_graph.BipartiteGraph'>
        sage: B.order()
        0
        sage: B == BipartiteGraph(None)
        True

    2. From a graph: without any more information, finds a bipartition::

        sage: B = BipartiteGraph(graphs.CycleGraph(4))
        sage: B = BipartiteGraph(graphs.CycleGraph(5))
        Traceback (most recent call last):
        ...
        TypeError: Input graph is not bipartite!
        sage: G = Graph({0:[5,6], 1:[4,5], 2:[4,6], 3:[4,5,6]})
        sage: B = BipartiteGraph(G)
        sage: B == G
        True
        sage: B.left
        {0, 1, 2, 3}
        sage: B.right
        {4, 5, 6}
        sage: B = BipartiteGraph({0:[5,6], 1:[4,5], 2:[4,6], 3:[4,5,6]})
        sage: B == G
        True
        sage: B.left
        {0, 1, 2, 3}
        sage: B.right
        {4, 5, 6}

    You can specify a partition using ``partition`` argument. Note that if such graph
    is not bipartite, then Sage will raise an error. However, if one specifies
    ``check=False``, the offending edges are simply deleted (along with
    those vertices not appearing in either list).  We also lump creating
    one bipartite graph from another into this category::

        sage: P = graphs.PetersenGraph()
        sage: partition = [list(range(5)), list(range(5,10))]
        sage: B = BipartiteGraph(P, partition)
        Traceback (most recent call last):
        ...
        TypeError: Input graph is not bipartite with respect to the given partition!

        sage: B = BipartiteGraph(P, partition, check=False)
        sage: B.left
        {0, 1, 2, 3, 4}
        sage: B.show()

      ::

        sage: G = Graph({0:[5,6], 1:[4,5], 2:[4,6], 3:[4,5,6]})
        sage: B = BipartiteGraph(G)
        sage: B2 = BipartiteGraph(B)
        sage: B == B2
        True
        sage: B3 = BipartiteGraph(G, [list(range(4)), list(range(4,7))])
        sage: B3
        Bipartite graph on 7 vertices
        sage: B3 == B2
        True

      ::

        sage: G = Graph({0:[], 1:[], 2:[]})
        sage: part = (list(range(2)), [2])
        sage: B = BipartiteGraph(G, part)
        sage: B2 = BipartiteGraph(B)
        sage: B == B2
        True

    4. From a reduced adjacency matrix::

        sage: M = Matrix([(1,1,1,0,0,0,0), (1,0,0,1,1,0,0),
        ....:             (0,1,0,1,0,1,0), (1,1,0,1,0,0,1)])
        sage: M
        [1 1 1 0 0 0 0]
        [1 0 0 1 1 0 0]
        [0 1 0 1 0 1 0]
        [1 1 0 1 0 0 1]
        sage: H = BipartiteGraph(M); H
        Bipartite graph on 11 vertices
        sage: H.edges()
        [(0, 7, None),
         (0, 8, None),
         (0, 10, None),
         (1, 7, None),
         (1, 9, None),
         (1, 10, None),
         (2, 7, None),
         (3, 8, None),
         (3, 9, None),
         (3, 10, None),
         (4, 8, None),
         (5, 9, None),
         (6, 10, None)]

      ::

        sage: M = Matrix([(1, 1, 2, 0, 0), (0, 2, 1, 1, 1), (0, 1, 2, 1, 1)])
        sage: B = BipartiteGraph(M, multiedges=True, sparse=True)
        sage: B.edges()
        [(0, 5, None),
         (1, 5, None),
         (1, 6, None),
         (1, 6, None),
         (1, 7, None),
         (2, 5, None),
         (2, 5, None),
         (2, 6, None),
         (2, 7, None),
         (2, 7, None),
         (3, 6, None),
         (3, 7, None),
         (4, 6, None),
         (4, 7, None)]

      ::

         sage: F.<a> = GF(4)
         sage: MS = MatrixSpace(F, 2, 3)
         sage: M = MS.matrix([[0, 1, a+1], [a, 1, 1]])
         sage: B = BipartiteGraph(M, weighted=True, sparse=True)
         sage: B.edges()
         [(0, 4, a), (1, 3, 1), (1, 4, 1), (2, 3, a + 1), (2, 4, 1)]
         sage: B.weighted()
         True

    5. From an alist file::

         sage: file_name = os.path.join(SAGE_TMP, 'deleteme.alist.txt')
         sage: fi = open(file_name, 'w')
         sage: fi.write("7 4 \n 3 4 \n 3 3 1 3 1 1 1 \n 3 3 3 4 \n\
                         1 2 4 \n 1 3 4 \n 1 0 0 \n 2 3 4 \n\
                         2 0 0 \n 3 0 0 \n 4 0 0 \n\
                         1 2 3 0 \n 1 4 5 0 \n 2 4 6 0 \n 1 2 4 7 \n")
         sage: fi.close();
         sage: B = BipartiteGraph(file_name)
         sage: B == H
         True

    6. From a NetworkX bipartite graph::

        sage: import networkx
        sage: G = graphs.OctahedralGraph()
        sage: N = networkx.make_clique_bipartite(G.networkx_graph())
        sage: B = BipartiteGraph(N)

    TESTS:

    Make sure we can create a ``BipartiteGraph`` with keywords but no
    positional arguments (:trac:`10958`).

    ::

        sage: B = BipartiteGraph(multiedges=True)
        sage: B.allows_multiple_edges()
        True

    Ensure that we can construct a ``BipartiteGraph`` with isolated vertices
    via the reduced adjacency matrix (:trac:`10356`)::

        sage: a=BipartiteGraph(matrix(2,2,[1,0,1,0]))
        sage: a
        Bipartite graph on 4 vertices
        sage: a.vertices()
        [0, 1, 2, 3]
        sage: g = BipartiteGraph(matrix(4,4,[1]*4+[0]*12))
        sage: g.vertices()
        [0, 1, 2, 3, 4, 5, 6, 7]
        sage: sorted(g.left.union(g.right))
        [0, 1, 2, 3, 4, 5, 6, 7]


    """

    def __init__(self, data=None, partition=None, check=True, *args, **kwds):
        """
        Create a bipartite graph. See documentation ``BipartiteGraph?`` for
        detailed information.

        EXAMPLE::

            sage: P = graphs.PetersenGraph()
            sage: partition = [list(range(5)), list(range(5,10))]
            sage: B = BipartiteGraph(P, partition, check=False)
        """
        if data is None:
            if partition is not None and check:
                if partition[0] or partition[1]:
                    raise ValueError("Invalid partition.")
            Graph.__init__(self, **kwds)
            self.left = set()
            self.right = set()
            return

        # need to turn off partition checking for Graph.__init__() adding
        # vertices and edges; methods are restored ad the end of big "if"
        # statement below
        import types
        self.add_vertex = types.MethodType(Graph.add_vertex,
                                           self,
                                           BipartiteGraph)
        self.add_vertices = types.MethodType(Graph.add_vertices,
                                             self,
                                             BipartiteGraph)
        self.add_edge = types.MethodType(Graph.add_edge, self, BipartiteGraph)

        from sage.structure.element import is_Matrix
        if isinstance(data, BipartiteGraph):
            Graph.__init__(self, data, *args, **kwds)
            self.left = set(data.left)
            self.right = set(data.right)
        elif isinstance(data, str):
            Graph.__init__(self, *args, **kwds)
            # will call self.load_afile after restoring add_vertex() instance
            # methods; initialize left and right attributes
            self.left = set()
            self.right = set()
        elif is_Matrix(data):
            # sanity check for mutually exclusive keywords
            if kwds.get("multiedges", False) and kwds.get("weighted", False):
                raise TypeError(
                    "Weighted multi-edge bipartite graphs from reduced " +
                    "adjacency matrix not supported.")
            Graph.__init__(self, *args, **kwds)
            ncols = data.ncols()
            nrows = data.nrows()
            self.left = set(range(ncols))
            self.right = set(range(ncols, nrows + ncols))

            # ensure that the vertices exist even if there
            # are no associated edges (trac #10356)
            self.add_vertices(self.left)
            self.add_vertices(self.right)

            if kwds.get("multiedges", False):
                for ii in range(ncols):
                    for jj in range(nrows):
                        if data[jj][ii] != 0:
                            self.add_edges([(ii, jj + ncols)] * data[jj][ii])
            elif kwds.get("weighted", False):
                for ii in range(ncols):
                    for jj in range(nrows):
                        if data[jj][ii] != 0:
                            self.add_edge((ii, jj + ncols, data[jj][ii]))
            else:
                for ii in range(ncols):
                    for jj in range(nrows):
                        if data[jj][ii] != 0:
                            self.add_edge((ii, jj + ncols))
        elif (isinstance(data, Graph) and partition is not None):
            from copy import copy
            left, right = partition
            left = copy(left)
            right = copy(right)
            verts = set(left) | set(right)
            if set(data.vertices()) != verts:
                data = data.subgraph(list(verts))
            Graph.__init__(self, data, *args, **kwds)
            if check:
                while len(left) > 0:
                    a = left.pop(0)
                    if len(set(data.neighbors(a)) & set(left)) != 0:
                        raise TypeError(
                            "Input graph is not bipartite with " +
                            "respect to the given partition!")
                while len(right) > 0:
                    a = right.pop(0)
                    if len(set(data.neighbors(a)) & set(right)) != 0:
                        raise TypeError(
                            "Input graph is not bipartite with " +
                            "respect to the given partition!")
            else:
                while len(left) > 0:
                    a = left.pop(0)
                    a_nbrs = set(data.neighbors(a)) & set(left)
                    if len(a_nbrs) != 0:
                        self.delete_edges([(a, b) for b in a_nbrs])
                while len(right) > 0:
                    a = right.pop(0)
                    a_nbrs = set(data.neighbors(a)) & set(right)
                    if len(a_nbrs) != 0:
                        self.delete_edges([(a, b) for b in a_nbrs])
            self.left, self.right = set(partition[0]), set(partition[1])
        elif isinstance(data, Graph):
            Graph.__init__(self, data, *args, **kwds)
            try:
                self.left, self.right = self.bipartite_sets()
            except Exception:
                raise TypeError("Input graph is not bipartite!")
        else:
            import networkx
            Graph.__init__(self, data, *args, **kwds)
            if isinstance(data, (networkx.MultiGraph, networkx.Graph)):
                if hasattr(data, "node_type"):
                    # Assume the graph is bipartite
                    self.left = set()
                    self.right = set()
                    for v in data.nodes_iter():
                        if data.node_type[v] == "Bottom":
                            self.left.add(v)
                        elif data.node_type[v] == "Top":
                            self.right.add(v)
                        else:
                            raise TypeError(
                                "NetworkX node_type defies bipartite " +
                                "assumption (is not 'Top' or 'Bottom')")
            # make sure we found a bipartition
            if not (hasattr(self, "left") and hasattr(self, "right")):
                try:
                    self.left, self.right = self.bipartite_sets()
                except Exception:
                    raise TypeError("Input graph is not bipartite!")

        # restore vertex partition checking
        self.add_vertex = types.MethodType(BipartiteGraph.add_vertex,
                                           self,
                                           BipartiteGraph)
        self.add_vertices = types.MethodType(BipartiteGraph.add_vertices,
                                             self,
                                             BipartiteGraph)
        self.add_edge = types.MethodType(BipartiteGraph.add_edge,
                                         self,
                                         BipartiteGraph)

        # post-processing
        if isinstance(data, str):
            self.load_afile(data)

        return

    def __repr__(self):
        r"""
        Returns a short string representation of self.

        EXAMPLE::

            sage: B = BipartiteGraph(graphs.CycleGraph(16))
            sage: B
            Bipartite cycle graph: graph on 16 vertices
        """
        s = Graph._repr_(self).lower()
        if "bipartite" in s:
            return s.capitalize()
        else:
            return "".join(["Bipartite ", s])

    def add_vertex(self, name=None, left=False, right=False):
        """
        Creates an isolated vertex. If the vertex already exists, then
        nothing is done.

        INPUT:

        - ``name`` -- (default: ``None``) name of the new vertex.  If no name
          is specified, then the vertex will be represented by the least
          non-negative integer not already representing a vertex.  Name must
          be an immutable object and cannot be ``None``.

        - ``left`` -- (default: ``False``) if ``True``, puts the new vertex
          in the left partition.

        - ``right`` -- (default: ``False``) if ``True``, puts the new vertex
          in the right partition.

        Obviously, ``left`` and ``right`` are mutually exclusive.

        As it is implemented now, if a graph `G` has a large number
        of vertices with numeric labels, then ``G.add_vertex()`` could
        potentially be slow, if name is ``None``.

        OUTPUT:

        - If ``name``=``None``, the new vertex name is returned. ``None`` otherwise.

        EXAMPLES::

            sage: G = BipartiteGraph()
            sage: G.add_vertex(left=True)
            0
            sage: G.add_vertex(right=True)
            1
            sage: G.vertices()
            [0, 1]
            sage: G.left
            {0}
            sage: G.right
            {1}

        TESTS:

        Exactly one of ``left`` and ``right`` must be true::

            sage: G = BipartiteGraph()
            sage: G.add_vertex()
            Traceback (most recent call last):
            ...
            RuntimeError: Partition must be specified (e.g. left=True).
            sage: G.add_vertex(left=True, right=True)
            Traceback (most recent call last):
            ...
            RuntimeError: Only one partition may be specified.

        Adding the same vertex must specify the same partition::

            sage: bg = BipartiteGraph()
            sage: bg.add_vertex(0, right=True)
            sage: bg.add_vertex(0, right=True)
            sage: bg.vertices()
            [0]
            sage: bg.add_vertex(0, left=True)
            Traceback (most recent call last):
            ...
            RuntimeError: Cannot add duplicate vertex to other partition.
        """
        # sanity check on partition specifiers
        if left and right:
            raise RuntimeError("Only one partition may be specified.")
        if not (left or right):
            raise RuntimeError("Partition must be specified (e.g. left=True).")

        # do nothing if we already have this vertex (idempotent)
        if (name is not None) and (name in self):
            if (((name in self.left) and left) or
                ((name in self.right) and right)):
                return
            else:
                raise RuntimeError(
                    "Cannot add duplicate vertex to other partition.")

        # add the vertex
        retval = Graph.add_vertex(self, name)
        if retval is not None: name = retval

        # add to proper partition
        if left:
            self.left.add(name)
        else:
            self.right.add(name)

        return retval

    def add_vertices(self, vertices, left=False, right=False):
        """
        Add vertices to the bipartite graph from an iterable container of
        vertices.  Vertices that already exist in the graph will not be added
        again.

        INPUT:

        - ``vertices`` -- sequence of vertices to add.

        - ``left`` -- (default: ``False``) either ``True`` or sequence of
          same length as ``vertices`` with ``True``/``False`` elements.

        - ``right`` -- (default: ``False``) either ``True`` or sequence of
          the same length as ``vertices`` with ``True``/``False`` elements.

        Only one of ``left`` and ``right`` keywords should be provided.  See
        the examples below.

        EXAMPLES::

            sage: bg = BipartiteGraph()
            sage: bg.add_vertices([0,1,2], left=True)
            sage: bg.add_vertices([3,4,5], left=[True, False, True])
            sage: bg.add_vertices([6,7,8], right=[True, False, True])
            sage: bg.add_vertices([9,10,11], right=True)
            sage: bg.left
            {0, 1, 2, 3, 5, 7}
            sage: bg.right
            {4, 6, 8, 9, 10, 11}

        TEST::

            sage: bg = BipartiteGraph()
            sage: bg.add_vertices([0,1,2], left=True)
            sage: bg.add_vertices([0,1,2], left=[True,True,True])
            sage: bg.add_vertices([0,1,2], right=[False,False,False])
            sage: bg.add_vertices([0,1,2], right=[False,False,False])
            sage: bg.add_vertices([0,1,2])
            Traceback (most recent call last):
            ...
            RuntimeError: Partition must be specified (e.g. left=True).
            sage: bg.add_vertices([0,1,2], left=True, right=True)
            Traceback (most recent call last):
            ...
            RuntimeError: Only one partition may be specified.
            sage: bg.add_vertices([0,1,2], right=True)
            Traceback (most recent call last):
            ...
            RuntimeError: Cannot add duplicate vertex to other partition.
            sage: (bg.left, bg.right)
            ({0, 1, 2}, set())
        """
        # sanity check on partition specifiers
        if left and right:  # also triggered if both lists are specified
            raise RuntimeError("Only one partition may be specified.")
        if not (left or right):
            raise RuntimeError("Partition must be specified (e.g. left=True).")

        # handle partitions
        if left and (not hasattr(left, "__iter__")):
            new_left = set(vertices)
            new_right = set()
        elif right and (not hasattr(right, "__iter__")):
            new_left = set()
            new_right = set(vertices)
        else:
            # simplify to always work with left
            if right:
                left = [not tf for tf in right]
            new_left = set()
            new_right = set()
            for tf, vv in zip(left, vertices):
                if tf:
                    new_left.add(vv)
                else:
                    new_right.add(vv)

        # check that we're not trying to add vertices to the wrong sets
        # or that a vertex is to be placed in both
        if ((new_left & self.right) or
            (new_right & self.left) or
            (new_right & new_left)):
            raise RuntimeError(
                "Cannot add duplicate vertex to other partition.")

        # add vertices
        Graph.add_vertices(self, vertices)
        self.left.update(new_left)
        self.right.update(new_right)

        return

    def delete_vertex(self, vertex, in_order=False):
        """
        Deletes vertex, removing all incident edges. Deleting a non-existent
        vertex will raise an exception.

        INPUT:

        - ``vertex`` -- a vertex to delete.

        - ``in_order`` -- (default ``False``) if ``True``, this deletes the
          `i`-th vertex in the sorted list of vertices,
          i.e. ``G.vertices()[i]``.

        EXAMPLES::

            sage: B = BipartiteGraph(graphs.CycleGraph(4))
            sage: B
            Bipartite cycle graph: graph on 4 vertices
            sage: B.delete_vertex(0)
            sage: B
            Bipartite cycle graph: graph on 3 vertices
            sage: B.left
            {2}
            sage: B.edges()
            [(1, 2, None), (2, 3, None)]
            sage: B.delete_vertex(3)
            sage: B.right
            {1}
            sage: B.edges()
            [(1, 2, None)]
            sage: B.delete_vertex(0)
            Traceback (most recent call last):
            ...
            RuntimeError: Vertex (0) not in the graph.

        ::

            sage: g = Graph({'a':['b'], 'c':['b']})
            sage: bg = BipartiteGraph(g)  # finds bipartition
            sage: bg.vertices()
            ['a', 'b', 'c']
            sage: bg.delete_vertex('a')
            sage: bg.edges()
            [('b', 'c', None)]
            sage: bg.vertices()
            ['b', 'c']
            sage: bg2 = BipartiteGraph(g)
            sage: bg2.delete_vertex(0, in_order=True)
            sage: bg2 == bg
            True
        """
        # cache vertex lookup if requested
        if in_order:
            vertex = self.vertices()[vertex]

        # delete from the graph
        Graph.delete_vertex(self, vertex)

        # now remove from partition (exception already thrown for non-existant
        # vertex)
        try:
            self.left.remove(vertex)
        except Exception:
            try:
                self.right.remove(vertex)
            except Exception:
                raise RuntimeError(
                    "Vertex (%s) not found in partitions" % vertex)

    def delete_vertices(self, vertices):
        """
        Remove vertices from the bipartite graph taken from an iterable
        sequence of vertices. Deleting a non-existent vertex will raise an
        exception.

        INPUT:

        - ``vertices`` -- a sequence of vertices to remove.

        EXAMPLES::

            sage: B = BipartiteGraph(graphs.CycleGraph(4))
            sage: B
            Bipartite cycle graph: graph on 4 vertices
            sage: B.delete_vertices([0,3])
            sage: B
            Bipartite cycle graph: graph on 2 vertices
            sage: B.left
            {2}
            sage: B.right
            {1}
            sage: B.edges()
            [(1, 2, None)]
            sage: B.delete_vertices([0])
            Traceback (most recent call last):
            ...
            RuntimeError: Vertex (0) not in the graph.
        """
        # remove vertices from the graph
        Graph.delete_vertices(self, vertices)

        # now remove vertices from partition lists (exception already thrown
        # for non-existant vertices)
        for vertex in vertices:
            try:
                self.left.remove(vertex)
            except Exception:
                try:
                    self.right.remove(vertex)
                except Exception:
                    raise RuntimeError(
                        "Vertex (%s) not found in partitions" % vertex)

    def add_edge(self, u, v=None, label=None):
        """
        Adds an edge from ``u`` and ``v``.

        INPUT:

        - ``u`` -- the tail of an edge.

        - ``v`` -- (default: ``None``) the head of an edge. If ``v=None``, then
          attempt to understand ``u`` as a edge tuple.

        - ``label`` -- (default: ``None``) the label of the edge ``(u, v)``.

        The following forms are all accepted:

        - ``G.add_edge(1, 2)``
        - ``G.add_edge((1, 2))``
        - ``G.add_edges([(1, 2)])``
        - ``G.add_edge(1, 2, 'label')``
        - ``G.add_edge((1, 2, 'label'))``
        - ``G.add_edges([(1, 2, 'label')])``

        See ``Graph.add_edge`` for more detail.

        This method simply checks that the edge endpoints are in different
        partitions. If a new vertex is to be created, it will be added
        to the proper partition. If both vertices are created, the first
        one will be added to the left partition, the second to the right
        partition.

        TEST::

            sage: bg = BipartiteGraph()
            sage: bg.add_vertices([0,1,2], left=[True,False,True])
            sage: bg.add_edges([(0,1), (2,1)])
            sage: bg.add_edge(0,2)
            Traceback (most recent call last):
            ...
            RuntimeError: Edge vertices must lie in different partitions.
            sage: bg.add_edge(0,3); list(bg.right)
            [1, 3]
            sage: bg.add_edge(5,6); 5 in bg.left; 6 in bg.right
            True
            True
        """
        # logic for getting endpoints copied from generic_graph.py
        if label is None:
            if v is None:
                try:
                    u, v, label = u
                except Exception:
                    u, v = u
                    label = None
        else:
            if v is None:
                u, v = u

        # check for endpoints in different partitions
        if self.left.issuperset((u, v)) or self.right.issuperset((u, v)):
            raise RuntimeError(
                "Edge vertices must lie in different partitions.")

        # automatically decide partitions for the newly created vertices
        if u not in self:
            self.add_vertex(u, left=(v in self.right or v not in self), right=(v in self.left))
        if v not in self:
            self.add_vertex(v, left=(u in self.right), right=(u in self.left))

        # add the edge
        Graph.add_edge(self, u, v, label)
        return

    def complement(self):
        """
        Return a complement of this graph.

        EXAMPLES::

            sage: B = BipartiteGraph({1: [2, 4], 3: [4, 5]})
            sage: G = B.complement(); G
            Graph on 5 vertices
            sage: G.edges(labels=False)
            [(1, 3), (1, 5), (2, 3), (2, 4), (2, 5), (4, 5)]
        """
        # This is needed because complement() of generic graph
        # would return a graph of class BipartiteGraph that is
        # not bipartite. See ticket #12376.
        return Graph(self).complement()

    def to_undirected(self):
        """
        Return an undirected Graph (without bipartite constraint) of the given
        object.

        EXAMPLES::

            sage: BipartiteGraph(graphs.CycleGraph(6)).to_undirected()
            Cycle graph: Graph on 6 vertices
        """
        return Graph(self)

    def bipartition(self):
        r"""
        Returns the underlying bipartition of the bipartite graph.

        EXAMPLE::

            sage: B = BipartiteGraph(graphs.CycleGraph(4))
            sage: B.bipartition()
            ({0, 2}, {1, 3})
        """
        return (self.left, self.right)

    def project_left(self):
        r"""
        Projects ``self`` onto left vertices. Edges are 2-paths in the
        original.

        EXAMPLE::

            sage: B = BipartiteGraph(graphs.CycleGraph(20))
            sage: G = B.project_left()
            sage: G.order(), G.size()
            (10, 10)
        """
        G = Graph()
        G.add_vertices(self.left)
        for v in G:
            for u in self.neighbor_iterator(v):
                G.add_edges([(v, w) for w in self.neighbor_iterator(u)])
        return G

    def project_right(self):
        r"""
        Projects ``self`` onto right vertices. Edges are 2-paths in the
        original.

        EXAMPLE::

            sage: B = BipartiteGraph(graphs.CycleGraph(20))
            sage: G = B.project_right()
            sage: G.order(), G.size()
            (10, 10)
        """
        G = Graph()
        G.add_vertices(self.left)
        for v in G:
            for u in self.neighbor_iterator(v):
                G.add_edges([(v, w) for w in self.neighbor_iterator(u)])
        return G

    def plot(self, *args, **kwds):
        r"""
        Overrides Graph's plot function, to illustrate the bipartite nature.

        EXAMPLE::

            sage: B = BipartiteGraph(graphs.CycleGraph(20))
            sage: B.plot()
            Graphics object consisting of 41 graphics primitives
        """
        if "pos" not in kwds:
            kwds["pos"] = None
        if kwds["pos"] is None:
            pos = {}
            left = list(self.left)
            right = list(self.right)
            left.sort()
            right.sort()
            l_len = len(self.left)
            r_len = len(self.right)
            if l_len == 1:
                pos[left[0]] = [-1, 0]
            elif l_len > 1:
                i = 0
                d = 2.0 / (l_len - 1)
                for v in left:
                    pos[v] = [-1, 1 - i*d]
                    i += 1
            if r_len == 1:
                pos[right[0]] = [1, 0]
            elif r_len > 1:
                i = 0
                d = 2.0 / (r_len - 1)
                for v in right:
                    pos[v] = [1, 1 - i*d]
                    i += 1
            kwds["pos"] = pos
        return Graph.plot(self, *args, **kwds)

    def matching_polynomial(self, algorithm="Godsil", name=None):
        r"""
        Computes the matching polynomial.

        If `p(G, k)` denotes the number of `k`-matchings (matchings with `k` edges)
        in `G`, then the *matching polynomial* is defined as [Godsil93]_:

        .. MATH::

            \mu(x)=\sum_{k \geq 0} (-1)^k p(G,k) x^{n-2k}

        INPUT:

        - ``algorithm`` - a string which must be either "Godsil" (default)
          or "rook"; "rook" is usually faster for larger graphs.

        - ``name`` - optional string for the variable name in the polynomial.

        EXAMPLE::

            sage: BipartiteGraph(graphs.CubeGraph(3)).matching_polynomial()
            x^8 - 12*x^6 + 42*x^4 - 44*x^2 + 9

        ::

            sage: x = polygen(QQ)
            sage: g = BipartiteGraph(graphs.CompleteBipartiteGraph(16, 16))
            sage: bool(factorial(16)*laguerre(16,x^2) == g.matching_polynomial(algorithm='rook'))
            True

        Compute the matching polynomial of a line with `60` vertices::

            sage: from sage.functions.orthogonal_polys import chebyshev_U
            sage: g = next(graphs.trees(60))
            sage: chebyshev_U(60, x/2) == BipartiteGraph(g).matching_polynomial(algorithm='rook')
            True

        The matching polynomial of a tree graphs is equal to its characteristic
        polynomial::

            sage: g = graphs.RandomTree(20)
            sage: p = g.characteristic_polynomial()
            sage: p == BipartiteGraph(g).matching_polynomial(algorithm='rook')
            True

        TESTS::

            sage: g = BipartiteGraph(matrix.ones(4,3))
            sage: g.matching_polynomial()
            x^7 - 12*x^5 + 36*x^3 - 24*x
            sage: g.matching_polynomial(algorithm="rook")
            x^7 - 12*x^5 + 36*x^3 - 24*x
        """
        if algorithm == "Godsil":
            return Graph.matching_polynomial(self, complement=False, name=name)
        elif algorithm == "rook":
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            A = self.reduced_adjacency_matrix()
            a = A.rook_vector()
            m = A.nrows()
            n = A.ncols()
            b = [0]*(m + n + 1)
            for i in range(min(m, n) + 1):
                b[m + n - 2*i] = a[i]*(-1)**i
            if name is None:
                name = 'x'
            K = PolynomialRing(A.base_ring(), name)
            p = K(b)
            return p
        else:
            raise ValueError('algorithm must be one of "Godsil" or "rook".')

    def load_afile(self, fname):
        r"""
        Loads into the current object the bipartite graph specified in the
        given file name.  This file should follow David MacKay's alist format,
        see
        http://www.inference.phy.cam.ac.uk/mackay/codes/data.html
        for examples and definition of the format.

        EXAMPLE::

            sage: file_name = os.path.join(SAGE_TMP, 'deleteme.alist.txt')
            sage: fi = open(file_name, 'w')
            sage: fi.write("7 4 \n 3 4 \n 3 3 1 3 1 1 1 \n 3 3 3 4 \n\
                            1 2 4 \n 1 3 4 \n 1 0 0 \n 2 3 4 \n\
                            2 0 0 \n 3 0 0 \n 4 0 0 \n\
                            1 2 3 0 \n 1 4 5 0 \n 2 4 6 0 \n 1 2 4 7 \n")
            sage: fi.close();
            sage: B = BipartiteGraph()
            sage: B.load_afile(file_name)
            Bipartite graph on 11 vertices
            sage: B.edges()
            [(0, 7, None),
             (0, 8, None),
             (0, 10, None),
             (1, 7, None),
             (1, 9, None),
             (1, 10, None),
             (2, 7, None),
             (3, 8, None),
             (3, 9, None),
             (3, 10, None),
             (4, 8, None),
             (5, 9, None),
             (6, 10, None)]
             sage: B2 = BipartiteGraph(file_name)
             sage: B2 == B
             True
        """
        # open the file
        try:
            fi = open(fname, "r")
        except IOError:
            print("Unable to open file <<" + fname + ">>.")
            return None

        # read header information
        num_cols, num_rows = [int(_) for _ in fi.readline().split()]
        max_col_degree, max_row_degree = [int(_) for _ in fi.readline().split()]
        col_degrees = [int(_) for _ in fi.readline().split()]
        row_degrees = [int(_) for _ in fi.readline().split()]

        # sanity checks on header info
        if len(col_degrees) != num_cols:
            print("Invalid Alist format: ")
            print("Number of column degree entries does not match number " +
                  "of columns.")
            return None
        if len(row_degrees) != num_rows:
            print("Invalid Alist format: ")
            print("Number of row degree entries does not match number " +
                  "of rows.")
            return None

        # clear out self
        self.clear()
        self.add_vertices(range(num_cols), left=True)
        self.add_vertices(range(num_cols, num_cols + num_rows), right=True)

        # read adjacency information
        for cidx in range(num_cols):
            for ridx in map(int, fi.readline().split()):
                # A-list uses 1-based indices with 0's as place-holders
                if ridx > 0:
                    self.add_edge(cidx, num_cols + ridx - 1)

        #NOTE:: we could read in the row adjacency information as well to
        #       double-check....
        #NOTE:: we could check the actual node degrees against the reported
        #       node degrees....

        # now we have all the edges in our graph, just fill in the
        # bipartite partitioning
        self.left = set(range(num_cols))
        self.right = set(range(num_cols, num_cols + num_rows))

        # return self for chaining calls if desired
        return self

    def save_afile(self, fname):
        r"""
        Save the graph to file in alist format.

        Saves this graph to file in David MacKay's alist format, see
        http://www.inference.phy.cam.ac.uk/mackay/codes/data.html
        for examples and definition of the format.

        EXAMPLE::

            sage: M = Matrix([(1,1,1,0,0,0,0), (1,0,0,1,1,0,0),
            ....:             (0,1,0,1,0,1,0), (1,1,0,1,0,0,1)])
            sage: M
            [1 1 1 0 0 0 0]
            [1 0 0 1 1 0 0]
            [0 1 0 1 0 1 0]
            [1 1 0 1 0 0 1]
            sage: b = BipartiteGraph(M)
            sage: file_name = os.path.join(SAGE_TMP, 'deleteme.alist.txt')
            sage: b.save_afile(file_name)
            sage: b2 = BipartiteGraph(file_name)
            sage: b == b2
            True

        TESTS::

            sage: file_name = os.path.join(SAGE_TMP, 'deleteme.alist.txt')
            sage: for order in range(3, 13, 3):
            ....:     num_chks = int(order / 3)
            ....:     num_vars = order - num_chks
            ....:     partition = (list(range(num_vars)), list(range(num_vars, num_vars+num_chks)))
            ....:     for idx in range(100):
            ....:         g = graphs.RandomGNP(order, 0.5)
            ....:         try:
            ....:             b = BipartiteGraph(g, partition, check=False)
            ....:             b.save_afile(file_name)
            ....:             b2 = BipartiteGraph(file_name)
            ....:             if b != b2:
            ....:                 print("Load/save failed for code with edges:")
            ....:                 print(b.edges())
            ....:                 break
            ....:         except Exception:
            ....:             print("Exception encountered for graph of order "+ str(order))
            ....:             print("with edges: ")
            ....:             g.edges()
            ....:             raise
        """
        # open the file
        try:
            fi = open(fname, "w")
        except IOError:
            print("Unable to open file <<" + fname + ">>.")
            return

        # prep: handy lists, functions for extracting adjacent nodes
        vnodes = list(self.left)
        cnodes = list(self.right)
        vnodes.sort()
        cnodes.sort()
        max_vdeg = max(self.degree(vnodes))
        max_cdeg = max(self.degree(cnodes))
        num_vnodes = len(vnodes)
        vnbr_str = lambda idx: str(idx - num_vnodes + 1)
        cnbr_str = lambda idx: str(idx + 1)

        # write header information
        fi.write("%d %d\n" % (len(vnodes), len(cnodes)))
        fi.write("%d %d\n" % (max_vdeg, max_cdeg))
        fi.write(" ".join(map(str, self.degree(vnodes))) + "\n")
        fi.write(" ".join(map(str, self.degree(cnodes))) + "\n")
        for vidx in vnodes:
            nbrs = self.neighbors(vidx)
            fi.write(" ".join(map(vnbr_str, nbrs)))
            fi.write(" 0"*(max_vdeg - len(nbrs)) + "\n")
        for cidx in cnodes:
            nbrs = self.neighbors(cidx)
            fi.write(" ".join(map(cnbr_str, nbrs)))
            fi.write(" 0"*(max_cdeg - len(nbrs)) + "\n")

        # done
        fi.close()

        # return self for chaining calls if desired
        return

    def __edge2idx(self, v1, v2, left, right):
        r"""
        Translate an edge to its reduced adjacency matrix position.

        Returns (row index, column index) for the given pair of vertices.

        EXAMPLE::

            sage: P = graphs.PetersenGraph()
            sage: partition = [list(range(5)), list(range(5,10))]
            sage: B = BipartiteGraph(P, partition, check=False)
            sage: B._BipartiteGraph__edge2idx(2,7,list(range(5)),list(range(5,10)))
            (2, 2)
        """
        try:
            if v1 in self.left:  # note uses attribute for faster lookup
                return (right.index(v2), left.index(v1))
            else:
                return (right.index(v1), left.index(v2))
        except ValueError:
            raise ValueError(
                "Tried to map invalid edge (%d,%d) to vertex indices" %
                (v1, v2))

    def reduced_adjacency_matrix(self, sparse=True):
        r"""
        Return the reduced adjacency matrix for the given graph.

        A reduced adjacency matrix contains only the non-redundant portion of
        the full adjacency matrix for the bipartite graph.  Specifically, for
        zero matrices of the appropriate size, for the reduced adjacency
        matrix ``H``, the full adjacency matrix is ``[[0, H'], [H, 0]]``.

        This method supports the named argument 'sparse' which defaults to
        ``True``.  When enabled, the returned matrix will be sparse.

        EXAMPLES:

        Bipartite graphs that are not weighted will return a matrix over ZZ::

            sage: M = Matrix([(1,1,1,0,0,0,0), (1,0,0,1,1,0,0),
            ....:             (0,1,0,1,0,1,0), (1,1,0,1,0,0,1)])
            sage: B = BipartiteGraph(M)
            sage: N = B.reduced_adjacency_matrix()
            sage: N
            [1 1 1 0 0 0 0]
            [1 0 0 1 1 0 0]
            [0 1 0 1 0 1 0]
            [1 1 0 1 0 0 1]
            sage: N == M
            True
            sage: N[0,0].parent()
            Integer Ring

        Multi-edge graphs also return a matrix over ZZ::

            sage: M = Matrix([(1,1,2,0,0), (0,2,1,1,1), (0,1,2,1,1)])
            sage: B = BipartiteGraph(M, multiedges=True, sparse=True)
            sage: N = B.reduced_adjacency_matrix()
            sage: N == M
            True
            sage: N[0,0].parent()
            Integer Ring

        Weighted graphs will return a matrix over the ring given by their
        (first) weights::

            sage: F.<a> = GF(4)
            sage: MS = MatrixSpace(F, 2, 3)
            sage: M = MS.matrix([[0, 1, a+1], [a, 1, 1]])
            sage: B = BipartiteGraph(M, weighted=True, sparse=True)
            sage: N = B.reduced_adjacency_matrix(sparse=False)
            sage: N == M
            True
            sage: N[0,0].parent()
            Finite Field in a of size 2^2

        TESTS::

            sage: B = BipartiteGraph()
            sage: B.reduced_adjacency_matrix()
            []
            sage: M = Matrix([[0,0], [0,0]])
            sage: BipartiteGraph(M).reduced_adjacency_matrix() == M
            True
            sage: M = Matrix([[10,2/3], [0,0]])
            sage: B = BipartiteGraph(M, weighted=True, sparse=True)
            sage: M == B.reduced_adjacency_matrix()
            True

        """
        if self.multiple_edges() and self.weighted():
            raise NotImplementedError(
                "Don't know how to represent weights for a multigraph.")
        if self.is_directed():
            raise NotImplementedError(
                "Reduced adjacency matrix does not exist for directed graphs.")

        # create sorted lists of left and right edges
        left = list(self.left)
        right = list(self.right)
        left.sort()
        right.sort()

        # create dictionary of edges, values are weights for weighted graph,
        # otherwise the number of edges (always 1 for simple graphs)
        D = {}
        if self.weighted():
            for (v1, v2, weight) in self.edge_iterator():
                D[self.__edge2idx(v1, v2, left, right)] = weight
        else:
            # if we're normal or multi-edge, just create the matrix over ZZ
            for (v1, v2, name) in self.edge_iterator():
                idx = self.__edge2idx(v1, v2, left, right)
                if idx in D:
                    D[idx] = 1 + D[idx]
                else:
                    D[idx] = 1

        # now construct and return the matrix from the dictionary we created
        from sage.matrix.constructor import matrix
        return matrix(len(self.right), len(self.left), D, sparse=sparse)
