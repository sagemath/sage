# autopep8: off
r"""
Bipartite graphs

This module implements bipartite graphs.

AUTHORS:

- Robert L. Miller (2008-01-20): initial version

- Ryan W. Hinton (2010-03-04): overrides for adding and deleting vertices
  and edges

- Enjeck M. Cleopatra (2022): fixes incorrect partite sets and adds graph
  creation from graph6 string

TESTS::

    sage: B = graphs.CompleteBipartiteGraph(7, 9)
    sage: loads(dumps(B)) == B
    True

    sage: B = BipartiteGraph(graphs.CycleGraph(4))
    sage: B == B.copy()
    True
    sage: type(B.copy())
    <class 'sage.graphs.bipartite_graph.BipartiteGraph'>
"""
# ****************************************************************************
#         Copyright (C) 2008 Robert L. Miller <rlmillster@gmail.com>
#                       2018 Julian Rüth <julian.rueth@fsfe.org>
#                       2022 Enjeck M. Cleopatra <enjeckc1e0@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from collections import defaultdict
from collections.abc import Iterable
import itertools

from .generic_graph import GenericGraph
from .graph import Graph
from sage.rings.integer import Integer
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import

lazy_import('networkx', ['MultiGraph', 'Graph'], as_=['networkx_MultiGraph', 'networkx_Graph'])


class BipartiteGraph(Graph):
    r"""
    Bipartite graph.

    INPUT:

    - ``data`` -- can be any of the following:

      #. Empty or ``None`` (creates an empty graph).

      #. An arbitrary graph.

      #. A reduced adjacency matrix.

         A reduced adjacency matrix contains only the non-redundant
         portion of the full adjacency matrix for the bipartite graph.
         Specifically, for zero matrices of the appropriate size, for
         the reduced adjacency matrix ``H``, the full adjacency matrix
         is ``[[0, H'], [H, 0]]``. The columns correspond to vertices
         on the left, and the rows correspond to vertices on the
         right.

      #. A file in alist format.

         The alist file format is described at
         http://www.inference.phy.cam.ac.uk/mackay/codes/alist.html

      #. A ``graph6`` string (see documentation of :meth:`~graph6_string`).

      #. From a NetworkX bipartite graph.

    - ``partition`` -- (default: ``None``) a tuple defining vertices of the left
      and right partition of the graph. Partitions will be determined
      automatically if ``partition`` is ``None``.

    - ``check`` -- boolean (default: ``True``); if ``True``, an invalid input
      partition raises an exception. In the other case offending edges simply
      won't be included.

    - ``loops`` -- ignored; bipartite graphs cannot have loops

    - ``multiedges`` -- boolean (default: ``None``); whether to allow multiple
      edges

    - ``weighted`` -- boolean (default: ``None``); whether graph thinks of
      itself as weighted or not. See ``self.weighted()``

    - ``hash_labels`` -- boolean (default: ``None``); whether to include edge
      labels during hashing. This parameter defaults to ``True`` if the graph is
      weighted. This parameter is ignored if the graph is mutable.
      Beware that trying to hash unhashable labels will raise an error.

    .. NOTE::

        All remaining arguments are passed to the ``Graph`` constructor

    EXAMPLES:

    #. No inputs or ``None`` for the input creates an empty graph::

        sage: B = BipartiteGraph()
        sage: type(B)
        <class 'sage.graphs.bipartite_graph.BipartiteGraph'>
        sage: B.order()
        0
        sage: B == BipartiteGraph(None)
        True

    #. From a graph: without any more information, finds a bipartition::

        sage: B = BipartiteGraph(graphs.CycleGraph(4))
        sage: B = BipartiteGraph(graphs.CycleGraph(5))
        Traceback (most recent call last):
        ...
        ValueError: input graph is not bipartite
        sage: G = Graph({0: [5, 6], 1: [4, 5], 2: [4, 6], 3: [4, 5, 6]})
        sage: B = BipartiteGraph(G)
        sage: B == G
        True
        sage: B.left
        {0, 1, 2, 3}
        sage: B.right
        {4, 5, 6}
        sage: B = BipartiteGraph({0: [5, 6], 1: [4, 5], 2: [4, 6], 3: [4, 5, 6]})
        sage: B == G
        True
        sage: B.left
        {0, 1, 2, 3}
        sage: B.right
        {4, 5, 6}

    #. If a ``Graph`` or ``DiGraph`` is used as data, you can specify
       a partition using ``partition`` argument. Note that if such
       graph is not bipartite, then Sage will raise an error. However,
       if one specifies ``check=False``, the offending edges are
       simply deleted (along with those vertices not appearing in
       either list).  We also lump creating one bipartite graph from
       another into this category::

        sage: P = graphs.PetersenGraph()
        sage: partition = [list(range(5)), list(range(5, 10))]
        sage: B = BipartiteGraph(P, partition)
        Traceback (most recent call last):
        ...
        TypeError: input graph is not bipartite with respect to the given partition

        sage: B = BipartiteGraph(P, partition, check=False)
        sage: B.left
        {0, 1, 2, 3, 4}
        sage: B.show()                                                                  # needs sage.plot

       ::

        sage: G = Graph({0: [5, 6], 1: [4, 5], 2: [4, 6], 3: [4, 5, 6]})
        sage: B = BipartiteGraph(G)
        sage: B2 = BipartiteGraph(B)
        sage: B == B2
        True
        sage: B3 = BipartiteGraph(G, [list(range(4)), list(range(4, 7))])
        sage: B3
        Bipartite graph on 7 vertices
        sage: B3 == B2
        True

       ::

        sage: G = Graph({0: [], 1: [], 2: []})
        sage: part = (list(range(2)), [2])
        sage: B = BipartiteGraph(G, part)
        sage: B2 = BipartiteGraph(B)
        sage: B == B2
        True

       ::

        sage: d = DiGraph(6)
        sage: d.add_edge(0, 1)
        sage: part=[[1, 2, 3], [0, 4, 5]]
        sage: b = BipartiteGraph(d, part)
        sage: b.left
        {1, 2, 3}
        sage: b.right
        {0, 4, 5}

    #. From a reduced adjacency matrix::

        sage: # needs sage.modules
        sage: M = Matrix([(1,1,1,0,0,0,0), (1,0,0,1,1,0,0),
        ....:             (0,1,0,1,0,1,0), (1,1,0,1,0,0,1)])
        sage: M
        [1 1 1 0 0 0 0]
        [1 0 0 1 1 0 0]
        [0 1 0 1 0 1 0]
        [1 1 0 1 0 0 1]
        sage: H = BipartiteGraph(M); H
        Bipartite graph on 11 vertices
        sage: H.edges(sort=True)
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

        sage: M = Matrix([(1, 1, 2, 0, 0), (0, 2, 1, 1, 1), (0, 1, 2, 1, 1)])           # needs sage.modules
        sage: B = BipartiteGraph(M, multiedges=True, sparse=True)                       # needs sage.modules
        sage: B.edges(sort=True)                                                        # needs sage.modules
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

         sage: # needs sage.modules sage.rings.finite_rings
         sage: F.<a> = GF(4)
         sage: MS = MatrixSpace(F, 2, 3)
         sage: M = MS.matrix([[0, 1, a + 1], [a, 1, 1]])
         sage: B = BipartiteGraph(M, weighted=True, sparse=True)
         sage: B.edges(sort=True)
         [(0, 4, a), (1, 3, 1), (1, 4, 1), (2, 3, a + 1), (2, 4, 1)]
         sage: B.weighted()
         True

    #. From an alist file::

         sage: import tempfile
         sage: with tempfile.NamedTemporaryFile(mode='w+t') as f:
         ....:     _ = f.write("7 4 \n 3 4 \n 3 3 1 3 1 1 1 \n\
         ....:                  3 3 3 4 \n 1 2 4 \n 1 3 4 \n 1 0 0 \n\
         ....:                  2 3 4 \n 2 0 0 \n 3 0 0 \n 4 0 0 \n\
         ....:                  1 2 3 0 \n 1 4 5 0 \n 2 4 6 0 \n\
         ....:                  1 2 4 7 \n")
         ....:     f.flush()
         ....:     B = BipartiteGraph(f.name)
         sage: B.is_isomorphic(H)                                                       # needs sage.modules
         True

    #. From a ``graph6`` string::

         sage: B = BipartiteGraph('Bo')
         sage: B
         Bipartite graph on 3 vertices
         sage: B.left
         {0}
         sage: B.right
         {1, 2}

       ::

         sage: B = BipartiteGraph('F?^T_\n', format='graph6')
         sage: B.vertices(sort=True)
         [0, 1, 2, 3, 4, 5, 6]
         sage: B.edges(sort=True)
         [(0, 5, None), (0, 6, None), (1, 4, None), (1, 5, None), (2, 4, None),
          (2, 6, None), (3, 4, None), (3, 5, None), (3, 6, None)]
         sage: B.left
         {0, 1, 2, 3}
         sage: B.right
         {4, 5, 6}

       ::
         sage: B = BipartiteGraph('Bo', partition=[[0], [1, 2]])
         sage: B.left
         {0}
         sage: B.right
         {1, 2}

       ::

         sage: B = BipartiteGraph('F?^T_\n', partition=[[0, 1, 2], [3, 4, 5, 6]])
         Traceback (most recent call last):
         ...
         TypeError: input graph is not bipartite with respect to the given partition

         sage: B = BipartiteGraph('F?^T_\n', partition=[[0, 1, 2], [3, 4, 5, 6]], check=False)
         sage: B.left
         {0, 1, 2}
         sage: B.show()                                                                 # needs sage.plot

    #. From a NetworkX bipartite graph::

        sage: # needs networkx
        sage: import networkx
        sage: G = graphs.OctahedralGraph()
        sage: N = networkx.make_clique_bipartite(G.networkx_graph())
        sage: B = BipartiteGraph(N)

    TESTS:

    Make sure we can create a ``BipartiteGraph`` with keywords but no positional
    arguments (:issue:`10958`)::

        sage: B = BipartiteGraph(multiedges=True)
        sage: B.allows_multiple_edges()
        True

    Ensure that we can construct a ``BipartiteGraph`` with isolated vertices via
    the reduced adjacency matrix (:issue:`10356`)::

        sage: # needs sage.modules
        sage: a = BipartiteGraph(matrix(2, 2, [1, 0, 1, 0]))
        sage: a
        Bipartite graph on 4 vertices
        sage: a.vertices(sort=True)
        [0, 1, 2, 3]
        sage: g = BipartiteGraph(matrix(4, 4, [1] * 4 + [0] * 12))
        sage: g.vertices(sort=True)
        [0, 1, 2, 3, 4, 5, 6, 7]
        sage: sorted(g.left.union(g.right))
        [0, 1, 2, 3, 4, 5, 6, 7]

    Make sure that loops are not allowed (:issue:`23275`)::

        sage: B = BipartiteGraph(loops=True)
        Traceback (most recent call last):
        ...
        ValueError: loops are not allowed in bipartite graphs
        sage: B = BipartiteGraph(loops=None)
        sage: B.allows_loops()
        False
        sage: B.add_edge(0, 0)
        Traceback (most recent call last):
        ...
        ValueError: cannot add edge from 0 to 0 in graph without loops
    """

    def __init__(self, data=None, partition=None, check=True, hash_labels=None, *args, **kwds):
        """
        Create a bipartite graph.

        See documentation ``BipartiteGraph?`` for detailed information.

        EXAMPLES::

            sage: P = graphs.PetersenGraph()
            sage: partition = [list(range(5)), list(range(5, 10))]
            sage: B = BipartiteGraph(P, partition, check=False)

        TESTS:

        Check that :issue:`33249` is fixed::

            sage: G = BipartiteGraph({2:[1], 3:[1], 4:[5]}, partition=([2,3,4],[1,5]))
            sage: print(G.left, G.right)
            {2, 3, 4} {1, 5}
            sage: G = BipartiteGraph({2:[1], 3:[1]}, partition=([1,2],[3]), check=True)
            Traceback (most recent call last):
            ...
            TypeError: input graph is not bipartite with respect to the given partition
            sage: G = BipartiteGraph({2:[1], 3:[1], 4:[5]}, partition=([2,3,4],[1]))
            Traceback (most recent call last):
            ...
            ValueError: not all vertices appear in partition
            sage: G = BipartiteGraph({2:[1], 3:[1], 4:[5]}, partition=([2,3,4],[1, 2]))
            Traceback (most recent call last):
            ...
            ValueError: the parts are not disjoint
            sage: G = BipartiteGraph({2:[1], 3:[1], 4:[5]}, partition=([2, 3, 4], [1, 7]))
            Traceback (most recent call last):
            ...
            LookupError: vertex (7) is not a vertex of the graph
        """
        if kwds is None:
            kwds = {'loops': False}
        else:
            if 'loops' in kwds and kwds['loops']:
                raise ValueError('loops are not allowed in bipartite graphs')
            kwds['loops'] = False

        if data is None:
            if partition is not None and check:
                if partition[0] or partition[1]:
                    raise ValueError("invalid partition")
            Graph.__init__(self, **kwds)
            self.left = set()
            self.right = set()
            self._hash_labels = hash_labels
            return

        # need to turn off partition checking for Graph.__init__() adding
        # vertices and edges; methods are restored at the end of big "if"
        # statement below
        from types import MethodType
        self.add_vertex = MethodType(Graph.add_vertex, self)
        self.add_vertices = MethodType(Graph.add_vertices, self)
        self.add_edge = MethodType(Graph.add_edge, self)
        self.add_edges = MethodType(Graph.add_edges, self)
        alist_file = True

        from sage.structure.element import Matrix
        if isinstance(data, BipartiteGraph):
            Graph.__init__(self, data, *args, **kwds)
            self.left = set(data.left)
            self.right = set(data.right)
        elif isinstance(data, str):
            import os
            alist_file = os.path.exists(data)
            Graph.__init__(self, data=None if alist_file else data, *args, **kwds)

            # methods; initialize left and right attributes
            self.left = set()
            self.right = set()

            # determine partitions and populate self.left and self.right
            if not alist_file:
                if partition is not None:
                    left, right = set(partition[0]), set(partition[1])

                # Some error checking.
                    if left & right:
                        raise ValueError("the parts are not disjoint")
                    if len(left) + len(right) != self.num_verts():
                        raise ValueError("not all vertices appear in partition")

                    if check:
                        if (any(left.intersection(self.neighbor_iterator(a)) for a in left) or
                                any(right.intersection(self.neighbor_iterator(a)) for a in right)):
                            raise TypeError("input graph is not bipartite with "
                                            "respect to the given partition")
                    else:
                        for a in left:
                            a_nbrs = left.intersection(self.neighbor_iterator(a))
                            if a_nbrs:
                                self.delete_edges((a, b) for b in a_nbrs)
                        for a in right:
                            a_nbrs = right.intersection(self.neighbor_iterator(a))
                            if a_nbrs:
                                self.delete_edges((a, b) for b in a_nbrs)
                    self.left, self.right = left, right
                else:
                    # Automatically get partitions if not provided
                    self._upgrade_from_graph()
        elif isinstance(data, Matrix):
            # sanity check for mutually exclusive keywords
            if kwds.get("multiedges", False) and kwds.get("weighted", False):
                raise TypeError("weighted multi-edge bipartite graphs from "
                                "reduced adjacency matrix not supported")
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
                        if data[jj, ii]:
                            self.add_edges([(ii, jj + ncols)] * data[jj, ii])
            elif kwds.get("weighted", False):
                for ii in range(ncols):
                    for jj in range(nrows):
                        if data[jj, ii]:
                            self.add_edge((ii, jj + ncols, data[jj, ii]))
            else:
                for ii in range(ncols):
                    for jj in range(nrows):
                        if data[jj, ii]:
                            self.add_edge((ii, jj + ncols))
        else:
            if partition is not None:
                left, right = set(partition[0]), set(partition[1])
                if isinstance(data, GenericGraph):
                    verts = left | right
                    if set(data) != verts:
                        data = data.subgraph(verts)
            Graph.__init__(self, data, *args, **kwds)
            if partition is not None:
                # Some error checking.
                if left & right:
                    raise ValueError("the parts are not disjoint")
                if len(left) + len(right) != self.num_verts():
                    raise ValueError("not all vertices appear in partition")

            if isinstance(data, (networkx_MultiGraph, networkx_Graph)):
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
                            raise TypeError("NetworkX node_type defies bipartite "
                                            "assumption (is not 'Top' or 'Bottom')")
            elif partition:
                if check:
                    if (any(left.intersection(self.neighbor_iterator(a)) for a in left) or
                            any(right.intersection(self.neighbor_iterator(a)) for a in right)):
                        raise TypeError("input graph is not bipartite with "
                                        "respect to the given partition")
                else:
                    for a in left:
                        a_nbrs = left.intersection(data.neighbor_iterator(a))
                        if a_nbrs:
                            self.delete_edges((a, b) for b in a_nbrs)
                    for a in right:
                        a_nbrs = right.intersection(data.neighbor_iterator(a))
                        if a_nbrs:
                            self.delete_edges((a, b) for b in a_nbrs)
                self.left, self.right = left, right

            # make sure we found a bipartition
            if not (hasattr(self, "left") and hasattr(self, "right")):
                self._upgrade_from_graph()

        # restore vertex partition checking
        del self.add_vertex
        del self.add_vertices
        del self.add_edge
        del self.add_edges

        # post-processing
        if isinstance(data, str):
            if alist_file:
                self.load_afile(data)

        if hash_labels is None and hasattr(data, '_hash_labels'):
            hash_labels = data._hash_labels
        self._hash_labels = hash_labels

        return

    @cached_method
    def __hash__(self):
        """
        Compute a hash for ``self``, if ``self`` is immutable.

        EXAMPLES::

            sage: A = BipartiteGraph([(1, 2, 1)], immutable=True)
            sage: B = BipartiteGraph([(1, 2, 33)], immutable=True)
            sage: A.__hash__() == B.__hash__()
            True
            sage: A = BipartiteGraph([(1, 2, 1)], immutable=True, hash_labels=True)
            sage: B = BipartiteGraph([(1, 2, 33)], immutable=True, hash_labels=True)
            sage: A.__hash__() == B.__hash__()
            False
            sage: A = BipartiteGraph([(1, 2, 1)], immutable=True, weighted=True)
            sage: B = BipartiteGraph([(1, 2, 33)], immutable=True, weighted=True)
            sage: A.__hash__() == B.__hash__()
            False

        TESTS::

            sage: A = BipartiteGraph([(1, 2, 1)], immutable=False)
            sage: A.__hash__()
            Traceback (most recent call last):
            ...
            TypeError: This graph is mutable, and thus not hashable. Create an immutable copy by `g.copy(immutable=True)`
            sage: B = BipartiteGraph([(1, 2, {'length': 3})], immutable=True, hash_labels=True)
            sage: B.__hash__()
            Traceback (most recent call last):
            ...
            TypeError: unhashable type: 'dict'
        """
        if self.is_immutable():
            # Determine whether to hash edge labels
            use_labels = self._use_labels_for_hash()
            edge_items = self.edge_iterator(labels=use_labels)
            if self.allows_multiple_edges():
                from collections import Counter
                edge_items = Counter(edge_items).items()
            return hash((frozenset(self.left), frozenset(self.right), frozenset(edge_items)))

        raise TypeError("This graph is mutable, and thus not hashable. "
                        "Create an immutable copy by `g.copy(immutable=True)`")

    def _upgrade_from_graph(self):
        """
        Set the left and right sets of vertices from the input graph.

        TESTS::

            sage: B = BipartiteGraph(Graph(1))
            sage: B.left, B.right
            ({0}, set())
            sage: B = BipartiteGraph(Graph(2))
            sage: B.left, B.right
            ({0, 1}, set())
            sage: B = BipartiteGraph(graphs.PathGraph(2))
            sage: B.left, B.right
            ({0}, {1})
            sage: B = BipartiteGraph(graphs.PathGraph(3))
            sage: B.left, B.right
            ({0, 2}, {1})
            sage: B = BipartiteGraph(graphs.CycleGraph(3))
            Traceback (most recent call last):
            ...
            ValueError: input graph is not bipartite
        """
        ans, certif = GenericGraph.is_bipartite(self, certificate=True)
        if not ans:
            raise ValueError("input graph is not bipartite")
        cols = defaultdict(set)
        for k, v in certif.items():
            cols[v].add(k)
        self.left = cols[1]
        self.right = cols[0]

    def __repr__(self):
        r"""
        Return a short string representation of ``self``.

        EXAMPLES::

            sage: B = BipartiteGraph(graphs.CycleGraph(16))
            sage: B
            Bipartite cycle graph: graph on 16 vertices
        """
        s = Graph._repr_(self).lower()
        if "bipartite" in s:
            return s.capitalize()
        return "".join(["Bipartite ", s])

    def add_vertex(self, name=None, left=False, right=False):
        """
        Create an isolated vertex. If the vertex already exists, then
        nothing is done.

        INPUT:

        - ``name`` -- (default: ``None``) name of the new vertex. If no name is
          specified, then the vertex will be represented by the least
          nonnegative integer not already representing a vertex. Name must be
          an immutable object and cannot be ``None``.

        - ``left`` -- boolean (default: ``False``); if ``True``, puts the new
          vertex in the left partition

        - ``right`` -- boolean (default: ``False``); if ``True``, puts the new
          vertex in the right partition

        Obviously, ``left`` and ``right`` are mutually exclusive.

        As it is implemented now, if a graph `G` has a large number of vertices
        with numeric labels, then ``G.add_vertex()`` could potentially be slow,
        if name is ``None``.

        OUTPUT:

        - If ``name`` is ``None``, the new vertex name is returned. ``None``
          otherwise.

        EXAMPLES::

            sage: G = BipartiteGraph()
            sage: G.add_vertex(left=True)
            0
            sage: G.add_vertex(right=True)
            1
            sage: G.vertices(sort=True)
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
            RuntimeError: partition must be specified (e.g. left=True)
            sage: G.add_vertex(left=True, right=True)
            Traceback (most recent call last):
            ...
            RuntimeError: only one partition may be specified

        Adding the same vertex must specify the same partition::

            sage: bg = BipartiteGraph()
            sage: bg.add_vertex(0, right=True)
            sage: bg.add_vertex(0, right=True)
            sage: bg.vertices(sort=False)
            [0]
            sage: bg.add_vertex(0, left=True)
            Traceback (most recent call last):
            ...
            RuntimeError: cannot add duplicate vertex to other partition
        """
        # sanity check on partition specifiers
        if left and right:
            raise RuntimeError("only one partition may be specified")
        if not (left or right):
            raise RuntimeError("partition must be specified (e.g. left=True)")

        # do nothing if we already have this vertex (idempotent)
        if name is not None and name in self:
            if ((left and name in self.left) or
                    (right and name in self.right)):
                return
            raise RuntimeError("cannot add duplicate vertex to other partition")

        # add the vertex
        retval = Graph.add_vertex(self, name)
        if retval is not None:
            name = retval

        # add to proper partition
        if left:
            self.left.add(name)
        else:
            self.right.add(name)

        return retval

    def add_vertices(self, vertices, left=False, right=False):
        """
        Add vertices to the bipartite graph from an iterable container of
        vertices.

        Vertices that already exist in the graph will not be added again.

        INPUT:

        - ``vertices`` -- sequence of vertices to add

        - ``left`` -- (default: ``False``) either ``True`` or sequence of same
          length as ``vertices`` with boolean elements

        - ``right`` -- (default: ``False``) either ``True`` or sequence of the
          same length as ``vertices`` with boolean elements

        Only one of ``left`` and ``right`` keywords should be provided.  See
        the examples below.

        EXAMPLES::

            sage: bg = BipartiteGraph()
            sage: bg.add_vertices([0, 1, 2], left=True)
            sage: bg.add_vertices([3, 4, 5], left=[True, False, True])
            sage: bg.add_vertices([6, 7, 8], right=[True, False, True])
            sage: bg.add_vertices([9, 10, 11], right=True)
            sage: bg.left
            {0, 1, 2, 3, 5, 7}
            sage: bg.right
            {4, 6, 8, 9, 10, 11}

        TESTS::

            sage: bg = BipartiteGraph()
            sage: bg.add_vertices([0, 1, 2], left=True)
            sage: bg.add_vertices([0, 1, 2], left=[True, True, True])
            sage: bg.add_vertices([0, 1, 2], right=[False, False, False])
            sage: bg.add_vertices([0, 1, 2], right=[False, False, False])
            sage: bg.add_vertices([0, 1, 2])
            Traceback (most recent call last):
            ...
            RuntimeError: partition must be specified (e.g. left=True)
            sage: bg.add_vertices([0,1,2], left=True, right=True)
            Traceback (most recent call last):
            ...
            RuntimeError: only one partition may be specified
            sage: bg.add_vertices([0,1,2], right=True)
            Traceback (most recent call last):
            ...
            RuntimeError: cannot add duplicate vertex to other partition
            sage: (bg.left, bg.right)
            ({0, 1, 2}, set())
        """
        # sanity check on partition specifiers
        if left and right:  # also triggered if both lists are specified
            raise RuntimeError("only one partition may be specified")
        if not (left or right):
            raise RuntimeError("partition must be specified (e.g. left=True)")

        # handle partitions
        if left and (not isinstance(left, Iterable)):
            new_left = set(vertices)
            new_right = set()
        elif right and (not isinstance(right, Iterable)):
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
            raise RuntimeError("cannot add duplicate vertex to other partition")

        # add vertices
        Graph.add_vertices(self, vertices)
        self.left.update(new_left)
        self.right.update(new_right)

    def delete_vertex(self, vertex, in_order=False):
        """
        Delete vertex, removing all incident edges.

        Deleting a non-existent vertex will raise an exception.

        INPUT:

        - ``vertex`` -- a vertex to delete

        - ``in_order`` -- boolean (default: ``False``); if ``True``, deletes the
          `i`-th vertex in the sorted list of vertices,
          i.e. ``G.vertices(sort=True)[i]``

        EXAMPLES::

            sage: B = BipartiteGraph(graphs.CycleGraph(4))
            sage: B
            Bipartite cycle graph: graph on 4 vertices
            sage: B.delete_vertex(0)
            sage: B
            Bipartite cycle graph: graph on 3 vertices
            sage: B.left
            {2}
            sage: B.edges(sort=True)
            [(1, 2, None), (2, 3, None)]
            sage: B.delete_vertex(3)
            sage: B.right
            {1}
            sage: B.edges(sort=True)
            [(1, 2, None)]
            sage: B.delete_vertex(0)
            Traceback (most recent call last):
            ...
            ValueError: vertex (0) not in the graph

        ::

            sage: g = Graph({'a': ['b'], 'c': ['b']})
            sage: bg = BipartiteGraph(g)  # finds bipartition
            sage: bg.vertices(sort=True)
            ['a', 'b', 'c']
            sage: bg.delete_vertex('a')
            sage: bg.edges(sort=True)
            [('b', 'c', None)]
            sage: bg.vertices(sort=True)
            ['b', 'c']
            sage: bg2 = BipartiteGraph(g)
            sage: bg2.delete_vertex(0, in_order=True)
            sage: bg2 == bg
            True
        """
        # cache vertex lookup if requested
        if in_order:
            vertex = self.vertices(sort=True)[vertex]

        # delete from the graph
        Graph.delete_vertex(self, vertex)

        # now remove from partition (exception already thrown for non-existent
        # vertex)
        if vertex in self.left:
            self.left.remove(vertex)
        elif vertex in self.right:
            self.right.remove(vertex)
        else:
            raise RuntimeError("vertex (%s) not found in partitions" % vertex)

    def delete_vertices(self, vertices):
        """
        Remove vertices from the bipartite graph taken from an iterable
        sequence of vertices.

        Deleting a non-existent vertex will raise an exception.

        INPUT:

        - ``vertices`` -- a sequence of vertices to remove

        EXAMPLES::

            sage: B = BipartiteGraph(graphs.CycleGraph(4))
            sage: B
            Bipartite cycle graph: graph on 4 vertices
            sage: B.delete_vertices([0, 3])
            sage: B
            Bipartite cycle graph: graph on 2 vertices
            sage: B.left
            {2}
            sage: B.right
            {1}
            sage: B.edges(sort=True)
            [(1, 2, None)]
            sage: B.delete_vertices([0])
            Traceback (most recent call last):
            ...
            ValueError: vertex (0) not in the graph
        """
        # remove vertices from the graph
        Graph.delete_vertices(self, vertices)

        # now remove vertices from partition lists (exception already thrown
        # for non-existent vertices)
        for vertex in vertices:
            if vertex in self.left:
                self.left.remove(vertex)
            elif vertex in self.right:
                self.right.remove(vertex)
            else:
                raise RuntimeError("vertex (%s) not found in partitions" % vertex)

    def _flip_vertices(self, vertices):
        r"""
        Helper method to flip the sides of a list of vertices.

        INPUT:

        - ``vertices`` -- an iterable container of vertices

        TESTS::

            sage: G = BipartiteGraph()
            sage: G.add_vertices([0, 1, 2], left=[True, False, True])
            sage: G.bipartition()
            ({0, 2}, {1})
            sage: G._flip_vertices([0, 1])
            sage: G.bipartition()
            ({1, 2}, {0})
            sage: G._flip_vertices([7])
            Traceback (most recent call last):
            ...
            RuntimeError: vertex (7) is neither in left nor in right
        """
        for vertex in vertices:
            if vertex in self.left:
                self.left.remove(vertex)
                self.right.add(vertex)
            elif vertex in self.right:
                self.right.remove(vertex)
                self.left.add(vertex)
            else:
                raise RuntimeError(f"vertex ({vertex}) is neither in left nor in right")

    def add_edge(self, u, v=None, label=None):
        r"""
        Add an edge from `u` to `v`.

        INPUT:

        - ``u`` -- the tail of an edge

        - ``v`` -- (default: ``None``) the head of an edge. If ``v=None``, then
          attempt to understand ``u`` as a edge tuple

        - ``label`` -- (default: ``None``) the label of the edge ``(u, v)``

        The following forms are all accepted:

        - ``G.add_edge(1, 2)``
        - ``G.add_edge((1, 2))``
        - ``G.add_edges([(1, 2)])``
        - ``G.add_edge(1, 2, 'label')``
        - ``G.add_edge((1, 2, 'label'))``
        - ``G.add_edges([(1, 2, 'label')])``

        See :meth:`~sage.graphs.graph.Graph.add_edge` for more detail.

        This method simply checks that the edge endpoints are in different
        partitions. If a new vertex is to be created, it will be added to the
        proper partition. If both vertices are created, the first one will be
        added to the left partition, the second to the right partition. If
        both vertices are in the same partition but different connected
        components, one of the components will be "flipped", i.e. each vertex
        will be put into whichever partition it's not currently in. This will
        allow for the graph to remain bipartite, without changing the edges or
        vertices.

        TESTS::

            sage: bg = BipartiteGraph()
            sage: bg.add_vertices([0, 1, 2], left=[True, False, True])
            sage: bg.add_edges([(0, 1), (2, 1)])
            sage: bg.add_edge(0, 2)
            Traceback (most recent call last):
            ...
            RuntimeError: edge vertices must lie in different partitions
            sage: bg.add_edge(0, 3); list(bg.right)
            [1, 3]
            sage: bg.add_edge(5, 6); 5 in bg.left; 6 in bg.right
            True
            True
            sage: G = BipartiteGraph()
            sage: G.add_edges([(0, 1), (3, 2)])
            sage: G.bipartition()
            ({0, 3}, {1, 2})
            sage: G.add_edge(1,2)
            sage: G.bipartition()
            ({0, 2}, {1, 3})
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

        # if endpoints are in the same partition
        if self.left.issuperset((u, v)) or self.right.issuperset((u, v)):

            # get v's connected component
            v_connected_component = self.connected_component_containing_vertex(v, sort=False)

            # if u is in it, then the edge still cannot exist
            if u in v_connected_component:
                raise RuntimeError("edge vertices must lie in different partitions")

            # if not, we can "flip" the connected component
            # swapping which partition the vertices are in
            self._flip_vertices(v_connected_component)

        # automatically decide partitions for the newly created vertices
        if u not in self:
            self.add_vertex(u, left=(v in self.right or v not in self), right=(v in self.left))
        if v not in self:
            self.add_vertex(v, left=(u in self.right), right=(u in self.left))

        # add the edge
        Graph.add_edge(self, u, v, label)
        return

    def add_edges(self, edges, loops=True):
        """
        Add edges from an iterable container.

        INPUT:

        - ``edges`` -- an iterable of edges, given either as ``(u, v)``
          or ``(u, v, label)``

        - ``loops`` -- ignored

        See :meth:`~sage.graphs.graph.Graph.add_edges` for more detail.

        This method simply checks that the edge endpoints are in different
        partitions. If a new vertex is to be created, it will be added to the
        proper partition. If both vertices are created, the first one will be
        added to the left partition, the second to the right partition. If
        both vertices are in the same partition but different connected
        components, one of the components will be "flipped", i.e. each vertex
        will be put into whichever partition it's not currently in. This will
        allow for the graph to remain bipartite, without changing the edges or
        vertices.

        EXAMPLES::

            sage: bg = BipartiteGraph()
            sage: bg.add_vertices([0, 1, 2], left=[True, False, True])
            sage: bg.add_edges([(0, 1), (2, 1)])
            sage: bg.add_edges([[0, 2]])
            Traceback (most recent call last):
            ...
            ValueError: the specified set of edges cannot be added
            while still preserving the bipartition property
            sage: G = BipartiteGraph()
            sage: G.add_edges([(0, 1), (3, 2), (1, 2)])
            sage: G.bipartition()
            ({0, 2}, {1, 3})


        Loops will raise an error::

            sage: bg.add_edges([[0, 3], [3, 3]])
            Traceback (most recent call last):
            ...
            ValueError: the specified set of edges cannot be added
            while still preserving the bipartition property

        Adding edges is fine as long as there exists a valid bipartition.
        Otherwise an error is raised without modifyiong the graph::

            sage: G = BipartiteGraph()
            sage: G.add_edges([(0, 1), (2, 3)])
            sage: G.bipartition()
            ({0, 2}, {1, 3})
            sage: G.add_edges([(0,2), (0,3)])
            Traceback (most recent call last):
            ...
            ValueError: the specified set of edges cannot be added
            while still preserving the bipartition property
            sage: G.bipartition()
            ({0, 2}, {1, 3})
            sage: G.edges(labels=False, sort=True)
            [(0, 1), (2, 3)]
        """
        edges_to_add = []
        for edge in edges:
            try:
                if len(edge) == 3:
                    u, v, label = edge
                else:
                    u, v = edge
                    label = None
                edges_to_add.append((u, v, label))
            except Exception:
                raise TypeError("cannot interpret {!r} as graph edge".format(edge))

        # Check whether there exists a bipartition supporting the addition of
        # input edges to the current graph before adding any edge to the
        # graph. This way, if an error is raised, self is not modified
        vertex_in_left = self._check_bipartition_for_add_edges(edges_to_add)

        if vertex_in_left is False:
            raise ValueError("the specified set of edges cannot be added while "
                             "still preserving the bipartition property")

        # If we get here, then we've found a valid bipartition.
        # We update the bipartition
        self.left.clear()
        self.right.clear()
        for v, left in vertex_in_left.items():
            if left:
                self.left.add(v)
            else:
                self.right.add(v)

        # Each edge now has one endpoint in left and the other in right
        for u, v, label in edges_to_add:
            if u not in self:
                self.add_vertex(u, left=vertex_in_left[u], right=not vertex_in_left[u])
            if v not in self:
                self.add_vertex(v, left=vertex_in_left[v], right=not vertex_in_left[v])

            self._backend.add_edge(u, v, label, self._directed)

    def _check_bipartition_for_add_edges(self, edges):
        r"""
        Helper method for ``add_edges``.

        This method checks whether the input list of edges can be added to the
        graph. More precisely, it checks whether there exists a bipartition of
        the vertices supporting the addition of input edges. If so it returns it
        as a mapping associating to each vertex a side of the bipartition.
        Otherwise, it returns ``False``.

        INPUT:

        - ``edges`` -- an iterable of edges, given either as ``(u, v)``
          or ``(u, v, label)``

        TESTS::

            sage: bg = BipartiteGraph()
            sage: bg.add_vertices([0, 1, 2, 3], left=[True, False, True, False])
            sage: b = bg._check_bipartition_for_add_edges([(0, 1), (3, 2), (1, 2)])
            sage: sorted(b.items())
            [(0, True), (1, False), (2, True), (3, False)]
            sage: b = bg._check_bipartition_for_add_edges([(0, 2)])
            sage: sorted(b.items())
            [(0, True), (1, False), (2, False), (3, False)]
            sage: bg.add_edges([(0, 1), (3, 2), (1, 2)])
            sage: bg._check_bipartition_for_add_edges([[0, 2]])
            False
        """
        # Map each vertex of the graph to a side
        vertex_in_left = {v: True for v in self.left}
        for v in self.right:
            vertex_in_left[v] = False

        # Map each vertex to the connected component it belongs to
        vertex_to_component = {v: comp for comp in self.connected_components(sort=False)
                                   for v in comp}

        for e in edges:
            u, v = e[:2]
            # if we haven't encountered either/both vertices, we choose a side
            # and extend components
            if u not in vertex_in_left:
                if v in vertex_in_left:
                    vertex_in_left[u] = not vertex_in_left[v]
                else:
                    vertex_in_left[u] = True
                    vertex_in_left[v] = False
                    vertex_to_component[v] = [v]
                vertex_to_component[v].append(u)
                vertex_to_component[u] = vertex_to_component[v]

            elif v not in vertex_in_left:
                vertex_in_left[v] = not vertex_in_left[u]
                vertex_to_component[u].append(v)
                vertex_to_component[v] = vertex_to_component[u]

            elif vertex_in_left[u] == vertex_in_left[v]:
                if vertex_to_component[u] is vertex_to_component[v]:
                    # Same side and same component. We can't add that edge
                    return False

                # Otherwise, we flip the bipartition in v's component
                for w in vertex_to_component[v]:
                    vertex_in_left[w] = not vertex_in_left[w]

                # and merge the components
                comp_u = vertex_to_component[u]
                comp_u.extend(vertex_to_component[v])
                for w in vertex_to_component[v]:
                    vertex_to_component[w] = comp_u

        # Return the bipartition
        return vertex_in_left

    def allow_loops(self, new, check=True):
        """
        Change whether loops are permitted in the (di)graph.

        .. NOTE::

            This method overwrite the
            :meth:`~sage.graphs.generic_graph.GenericGraph.allow_loops` method
            to ensure that loops are forbidden in :class:`~BipartiteGraph`.

        INPUT:

        - ``new`` -- boolean

        EXAMPLES::

            sage: B = BipartiteGraph()
            sage: B.allow_loops(True)
            Traceback (most recent call last):
            ...
            ValueError: loops are not allowed in bipartite graphs
        """
        if new is True:
            raise ValueError("loops are not allowed in bipartite graphs")

    def is_bipartite(self, certificate=False):
        r"""
        Check whether the graph is bipartite.

        This method always returns ``True`` as first value, plus a certificate
        when ``certificate == True``.

        INPUT:

        - ``certificate`` -- boolean (default: ``False``); whether to return a
          certificate. If set to ``True``, the certificate returned is a proper
          2-coloring of the vertices.

        .. SEEALSO:: :meth:`~GenericGraph.is_bipartite`

        EXAMPLES::

            sage: g = BipartiteGraph(graphs.RandomBipartite(3, 3, .5))                  # needs numpy
            sage: g.is_bipartite()                                                      # needs numpy
            True
            sage: g.is_bipartite(certificate=True)  # random                            # needs numpy
            (True, {(0, 0): 0, (0, 1): 0, (0, 2): 0, (1, 0): 1, (1, 1): 1, (1, 2): 1})

        TESTS::

            sage: BipartiteGraph().is_bipartite()
            True
            sage: BipartiteGraph().is_bipartite(certificate=True)
            (True, {})
        """
        if certificate:
            color = {u: 0 for u in self.left}
            color.update({u: 1 for u in self.right})
            return True, color
        return True

    def complement(self):
        r"""
        Return a complement of this graph.

        Given a simple :class:`~sage.graphs.bipartite_graph.BipartiteGraph`
        `G = (L, R, E)` with vertex set `L\cup R` and edge set `E`, this method
        returns a :class:`~sage.graphs.graph.Graph` `H = (V, F)`, where
        `V = L\cup R` and `F` is the set of edges of a complete graph of order
        `|V|` minus the edges in `E`.

        .. WARNING::

            This method returns the complement of a bipartite graph
            `G = (V = L \cup R, E)` with respect the complete graph of order
            `|V|`. If looking for the complement with respect the complete
            bipartite graph `K = (L, R, L\times R)`, use method
            :meth:`~sage.graphs.bipartite_graph.BipartiteGraph.complement_bipartite`.

        .. SEEALSO::

            :meth:`~sage.graphs.bipartite_graph.BipartiteGraph.complement_bipartite`

        EXAMPLES::

            sage: B = BipartiteGraph({1: [2, 4], 3: [4, 5]})
            sage: G = B.complement(); G
            Graph on 5 vertices
            sage: G.edges(sort=True, labels=False)
            [(1, 3), (1, 5), (2, 3), (2, 4), (2, 5), (4, 5)]
            sage: B.size() + G.size() == graphs.CompleteGraph(B.order()).size()
            True
        """
        # This is needed because complement() of generic graph
        # would return a graph of class BipartiteGraph that is
        # not bipartite. See issue #12376.
        return Graph(self).complement()

    def complement_bipartite(self):
        r"""
        Return the bipartite complement of this bipartite graph.

        Given a simple :class:`~sage.graphs.bipartite_graph.BipartiteGraph`
        `G = (L, R, E)` with vertex set `L\cup R` and edge set `E`, this
        method returns a :class:`~sage.graphs.bipartite_graph.BipartiteGraph`
        `H = (L\cup R, F)`, where `F` is the set of edges of a complete
        bipartite graph between vertex sets `L` and `R` minus the edges in `E`.

        .. SEEALSO::

            :meth:`~sage.graphs.bipartite_graph.BipartiteGraph.complement`

        EXAMPLES:

            sage: B = BipartiteGraph({0: [1, 2, 3]})
            sage: C = B.complement_bipartite()
            sage: C
            Bipartite graph on 4 vertices
            sage: C.is_bipartite()
            True
            sage: B.left == C.left and B.right == C.right
            True
            sage: C.size() == len(B.left)*len(B.right) - B.size()
            True
            sage: G = B.complement()
            sage: G.is_bipartite()
            False
        """
        self._scream_if_not_simple()

        E = [e for e in itertools.product(self.left, self.right) if not self.has_edge(e)]
        H = BipartiteGraph([self, E], format='vertices_and_edges', partition=[self.left, self.right])

        if self.name():
            H.name("complement-bipartite({})".format(self.name()))
        if self.is_immutable():
            return H.copy(immutable=True)
        return H

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
        Return the underlying bipartition of the bipartite graph.

        EXAMPLES::

            sage: B = BipartiteGraph(graphs.CycleGraph(4))
            sage: B.bipartition()
            ({0, 2}, {1, 3})
        """
        return (self.left, self.right)

    def project_left(self):
        r"""
        Project ``self`` onto left vertices. Edges are 2-paths in the original.

        EXAMPLES::

            sage: B = BipartiteGraph(graphs.CycleGraph(20))
            sage: G = B.project_left()
            sage: G.order(), G.size()
            (10, 10)
        """
        G = Graph()
        G.add_vertices(self.left)
        for v in G:
            for u in self.neighbor_iterator(v):
                G.add_edges(((v, w) for w in self.neighbor_iterator(u)), loops=False)
        return G

    def project_right(self):
        r"""
        Project ``self`` onto right vertices. Edges are 2-paths in the original.

        EXAMPLES::

            sage: B = BipartiteGraph(graphs.CycleGraph(20))
            sage: G = B.project_right()
            sage: G.order(), G.size()
            (10, 10)

        TESTS:

        Issue :issue:`25985` is fixed::

            sage: B = BipartiteGraph(graphs.CycleGraph(6))
            sage: B.project_left().vertices(sort=True)
            [0, 2, 4]
            sage: B.project_right().vertices(sort=True)
            [1, 3, 5]
        """
        G = Graph()
        G.add_vertices(self.right)
        for v in G:
            for u in self.neighbor_iterator(v):
                G.add_edges(((v, w) for w in self.neighbor_iterator(u)), loops=False)
        return G

    def plot(self, *args, **kwds):
        r"""
        Override Graph's plot function, to illustrate the bipartite nature.

        EXAMPLES::

            sage: B = BipartiteGraph(graphs.CycleGraph(20))
            sage: B.plot()                                                              # needs sage.plot
            Graphics object consisting of 41 graphics primitives
        """
        if "pos" not in kwds:
            kwds["pos"] = None
        if kwds["pos"] is None:
            if self.left:
                y = 0 if len(self.left) == 1 else 1
                pos = self._line_embedding(sorted(self.left), first=(-1, y), last=(-1, -y), return_dict=True)
            else:
                pos = {}
            if self.right:
                y = 0 if len(self.right) == 1 else 1
                pos.update(self._line_embedding(sorted(self.right), first=(1, y), last=(1, -y), return_dict=True))
            kwds["pos"] = pos
        return Graph.plot(self, *args, **kwds)

    def matching_polynomial(self, algorithm='Godsil', name=None):
        r"""
        Compute the matching polynomial.

        The *matching polynomial* is defined as in [God1993]_, where `p(G, k)`
        denotes the number of `k`-matchings (matchings with `k` edges) in `G` :

        .. MATH::

            \mu(x)=\sum_{k \geq 0} (-1)^k p(G,k) x^{n-2k}

        INPUT:

        - ``algorithm`` -- string (default: ``'Godsil'``); either "Godsil" or
          "rook"; "rook" is usually faster for larger graphs

        - ``name`` -- string (default: ``None``); name of the variable in the
          polynomial, set to `x` when ``name`` is ``None``

        EXAMPLES::

            sage: BipartiteGraph(graphs.CubeGraph(3)).matching_polynomial()             # needs sage.libs.flint
            x^8 - 12*x^6 + 42*x^4 - 44*x^2 + 9

        ::

            sage: x = polygen(QQ)
            sage: g = BipartiteGraph(graphs.CompleteBipartiteGraph(16, 16))
            sage: bool(factorial(16) * laguerre(16, x^2)                                # needs sage.symbolic
            ....:       == g.matching_polynomial(algorithm='rook'))
            True

        Compute the matching polynomial of a line with `60` vertices::

            sage: from sage.functions.orthogonal_polys import chebyshev_U               # needs sage.symbolic
            sage: g = next(graphs.trees(60))
            sage: (chebyshev_U(60, x/2)                                                 # needs sage.symbolic
            ....:   == BipartiteGraph(g).matching_polynomial(algorithm='rook'))
            True

        The matching polynomial of a tree is equal to its characteristic
        polynomial::

            sage: g = graphs.RandomTree(20)
            sage: p = g.characteristic_polynomial()                                     # needs sage.modules
            sage: p == BipartiteGraph(g).matching_polynomial(algorithm='rook')          # needs sage.modules
            True

        TESTS::

            sage: # needs sage.modules
            sage: g = BipartiteGraph(matrix.ones(4, 3))
            sage: g.matching_polynomial()                                               # needs sage.libs.flint
            x^7 - 12*x^5 + 36*x^3 - 24*x
            sage: g.matching_polynomial(algorithm='rook')
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
            b = [0] * (m + n + 1)
            for i in range(min(m, n) + 1):
                b[m + n - 2 * i] = a[i] * (-1)**i
            if name is None:
                name = 'x'
            K = PolynomialRing(A.base_ring(), name)
            p = K(b)
            return p
        raise ValueError('algorithm must be one of "Godsil" or "rook"')

    def perfect_matchings(self, labels=False):
        r"""
        Return an iterator over all perfect matchings of the bipartite graph.

        ALGORITHM:

        Choose a vertex `v` in the right set of vertices, then recurse through
        all edges incident to `v`, removing one edge at a time whenever an edge
        is added to a matching.

        INPUT:

        - ``labels`` -- boolean (default: ``False``); when ``True``, the edges
          in each perfect matching are triples (containing the label as the
          third element), otherwise the edges are pairs.

        .. SEEALSO::

            :meth:`~sage.graphs.graph.Graph.perfect_matchings`
            :meth:`matching`

        EXAMPLES::

            sage: B = BipartiteGraph({0: [5, 7], 1: [4, 6, 7], 2: [4, 5, 8],
            ....:                     3: [4, 5, 6], 6: [9], 8: [9]})
            sage: len(list(B.perfect_matchings()))
            6
            sage: G = Graph(B.edges(sort=False))
            sage: len(list(G.perfect_matchings()))
            6

        The algorithm ensures that for any edge of a perfect matching, the first
        vertex is on the left set of vertices and the second vertex in the right
        set::

            sage: B = BipartiteGraph({0: [5, 7], 1: [4, 6, 7], 2: [4, 5, 8],
            ....:                     3: [4, 5, 6], 6: [9], 8: [9]})
            sage: m = next(B.perfect_matchings(labels=False))
            sage: B.left
            {0, 1, 2, 3, 9}
            sage: B.right
            {4, 5, 6, 7, 8}
            sage: sorted(m)
            [(0, 7), (1, 4), (2, 5), (3, 6), (9, 8)]
            sage: all((u in B.left and v in B.right) for u, v in m)
            True

        Multiple edges are taken into account::

            sage: B = BipartiteGraph({0: [5, 7], 1: [4, 6, 7], 2: [4, 5, 8],
            ....:                     3: [4, 5, 6], 6: [9], 8: [9]})
            sage: B.allow_multiple_edges(True)
            sage: B.add_edge(0, 7)
            sage: len(list(B.perfect_matchings()))
            10


        Empty graph::

            sage: list(BipartiteGraph().perfect_matchings())
            [[]]

        Bipartite graph without perfect matching::

            sage: B = BipartiteGraph(graphs.CompleteBipartiteGraph(3, 4))
            sage: list(B.perfect_matchings())
            []

        Check that the number of perfect matchings of a complete bipartite graph
        is consistent with the matching polynomial::

            sage: B = BipartiteGraph(graphs.CompleteBipartiteGraph(4, 4))
            sage: len(list(B.perfect_matchings()))
            24
            sage: B.matching_polynomial(algorithm='rook')(0)                            # needs sage.modules
            24

        TESTS::

            sage: B = BipartiteGraph(graphs.CompleteBipartiteGraph(3, 4))
            sage: B.left, B.right
            ({0, 1, 2}, {3, 4, 5, 6})
            sage: B.add_vertex(left=True)
            7
            sage: B.left, B.right
            ({0, 1, 2, 7}, {3, 4, 5, 6})
            sage: list(B.perfect_matchings())
            []
            sage: B = BipartiteGraph(graphs.CompleteBipartiteGraph(3, 3))
            sage: B.add_vertex(left=True)
            6
            sage: B.add_vertex(right=True)
            7
            sage: list(B.perfect_matchings())
            []
            sage: G = Graph(B)
            sage: list(G.perfect_matchings())
            []
        """
        if not self:
            yield []
            return
        if len(self.left) != len(self.right):
            return

        def rec(G):
            """
            Iterator over all perfect matchings of a simple bipartite graph `G`.
            """
            if not G:
                yield []
                return
            # Take an element from the right set
            v = next(iter(G.right))
            Nv = G.neighbors(v)
            # We must have at least one vertex in Nv
            if not Nv:
                return
            G.delete_vertex(v)
            for u in Nv:
                Nu = G.neighbors(u)
                G.delete_vertex(u)
                for partial_matching in rec(G):
                    partial_matching.append((u, v))
                    yield partial_matching
                G.add_vertex(u, left=True)
                G.add_edges((u, nu) for nu in Nu)
            G.add_vertex(v, right=True)
            G.add_edges((nv, v) for nv in Nv)

        # We create a mutable copy of self
        G = self.copy(immutable=False)

        # We create a mapping from frozen unlabeled edges to (labeled) edges.
        # This ease for instance the manipulation of multiedges (if any)
        edges = {}
        for e in G.edges(sort=False, labels=labels):
            f = frozenset(e[:2])
            if e[0] not in G.left:
                e = (e[1], e[0], e[2]) if labels else (e[1], e[0])
            if f in edges:
                edges[f].append(e)
            else:
                edges[f] = [e]

        # We now get rid of multiple edges, if any
        G.allow_multiple_edges(False)

        # For each unlabeled matching, we yield all its possible labelings
        for m in rec(G):
            yield from itertools.product(*[edges[frozenset(e)] for e in m])

    def load_afile(self, fname):
        r"""
        Load into the current object the bipartite graph specified in the given
        file name.

        This file should follow David MacKay's alist format, see
        http://www.inference.phy.cam.ac.uk/mackay/codes/data.html for examples
        and definition of the format.

        EXAMPLES::

            sage: import tempfile
            sage: with tempfile.NamedTemporaryFile(mode='w+t') as f:
            ....:     _ = f.write("7 4 \n 3 4 \n 3 3 1 3 1 1 1 \n\
            ....:                 3 3 3 4 \n 1 2 4 \n 1 3 4 \n\
            ....:                 1 0 0 \n 2 3 4 \n 2 0 0 \n 3 0 0 \n\
            ....:                 4 0 0 \n 1 2 3 0 \n 1 4 5 0 \n\
            ....:                 2 4 6 0 \n 1 2 4 7 \n")
            ....:     f.flush()
            ....:     B = BipartiteGraph()
            ....:     B2 = BipartiteGraph(f.name)
            ....:     B.load_afile(f.name)
            Bipartite graph on 11 vertices
            sage: B.edges(sort=True)
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
             sage: B2 == B
             True
        """
        # open the file
        try:
            fi = open(fname)
        except OSError:
            print("unable to open file <<" + fname + ">>")
            return None

        # read header information
        num_cols, num_rows = (int(_) for _ in fi.readline().split())
        # next are max_col_degree, max_row_degree, not used
        _ = [int(_) for _ in fi.readline().split()]
        col_degrees = [int(_) for _ in fi.readline().split()]
        row_degrees = [int(_) for _ in fi.readline().split()]

        # sanity checks on header info
        if len(col_degrees) != num_cols:
            print("invalid Alist format: ")
            print("number of column degree entries does not match number of columns")
            return None
        if len(row_degrees) != num_rows:
            print("invalid Alist format: ")
            print("number of row degree entries does not match number of rows")
            return None

        # clear out self
        self.clear()
        self.add_vertices(range(num_cols), left=True)
        self.add_vertices(range(num_cols, num_cols + num_rows), right=True)

        # read adjacency information
        for cidx in range(num_cols):
            for ridx in map(int, fi.readline().split()):
                # A-list uses 1-based indices with 0s as place-holders
                if ridx > 0:
                    self.add_edge(cidx, num_cols + ridx - 1)

        # NOTE:: we could read in the row adjacency information as well to
        #        double-check....
        # NOTE:: we could check the actual node degrees against the reported
        #        node degrees....

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

        EXAMPLES::

            sage: # needs sage.modules
            sage: M = Matrix([(1,1,1,0,0,0,0), (1,0,0,1,1,0,0),
            ....:             (0,1,0,1,0,1,0), (1,1,0,1,0,0,1)])
            sage: M
            [1 1 1 0 0 0 0]
            [1 0 0 1 1 0 0]
            [0 1 0 1 0 1 0]
            [1 1 0 1 0 0 1]
            sage: b = BipartiteGraph(M)
            sage: import tempfile
            sage: with tempfile.NamedTemporaryFile() as f:
            ....:     b.save_afile(f.name)
            ....:     b2 = BipartiteGraph(f.name)
            sage: b.is_isomorphic(b2)
            True

        TESTS::

            sage: import tempfile
            sage: f = tempfile.NamedTemporaryFile()
            sage: for order in range(3, 13, 3):                                         # needs sage.combinat
            ....:     num_chks = int(order / 3)
            ....:     num_vars = order - num_chks
            ....:     partition = (list(range(num_vars)), list(range(num_vars, num_vars+num_chks)))
            ....:     for idx in range(100):
            ....:         g = graphs.RandomGNP(order, 0.5)
            ....:         try:
            ....:             b = BipartiteGraph(g, partition, check=False)
            ....:             b.save_afile(f.name)
            ....:             b2 = BipartiteGraph(f.name)
            ....:             if not b.is_isomorphic(b2):
            ....:                 print("Load/save failed for code with edges:")
            ....:                 print(b.edges(sort=True))
            ....:                 break
            ....:         except Exception:
            ....:             print("Exception encountered for graph of order "+ str(order))
            ....:             print("with edges: ")
            ....:             g.edges(sort=True)
            ....:             raise
            sage: f.close()  # this removes the file
        """
        # open the file
        try:
            fi = open(fname, "w")
        except OSError:
            print("Unable to open file <<" + fname + ">>.")
            return

        # The alist file format assumes that vertices are labeled in [0..n-1]
        # with:
        # - vertices in left labeled in [0..|left|-1]
        # - vertices in right labeled in [|left|..n-1]
        # Then, for vertex `vidx`, it stores in the file the index + 1 of a
        # neighbor in the list right. We use appropriate mappings.
        # See http://www.inference.org.uk/mackay/codes/alist.html

        # prep: handy lists, functions for extracting adjacent nodes
        vnodes = list(self.left)
        cnodes = list(self.right)
        max_vdeg = max(self.degree(vnodes))
        max_cdeg = max(self.degree(cnodes))
        vnode_to_str = {v: str(i + 1) for i, v in enumerate(vnodes)}
        cnode_to_str = {v: str(i + 1) for i, v in enumerate(cnodes)}

        def vnbr_str(idx):
            return cnode_to_str[idx]

        def cnbr_str(idx):
            return vnode_to_str[idx]

        # write header information
        fi.write("%d %d\n" % (len(vnodes), len(cnodes)))
        fi.write("%d %d\n" % (max_vdeg, max_cdeg))
        fi.write(" ".join(map(str, self.degree(vnodes))) + "\n")
        fi.write(" ".join(map(str, self.degree(cnodes))) + "\n")
        for vidx in vnodes:
            nbrs = self.neighbors(vidx)
            fi.write(" ".join(map(vnbr_str, nbrs)))
            fi.write(" 0" * (max_vdeg - len(nbrs)) + "\n")
        for cidx in cnodes:
            nbrs = self.neighbors(cidx)
            fi.write(" ".join(map(cnbr_str, nbrs)))
            fi.write(" 0" * (max_cdeg - len(nbrs)) + "\n")

        # done
        fi.close()

        # return self for chaining calls if desired
        return

    def reduced_adjacency_matrix(self, sparse=True, *, base_ring=None, **kwds):
        r"""
        Return the reduced adjacency matrix for the given graph.

        A reduced adjacency matrix contains only the non-redundant portion of
        the full adjacency matrix for the bipartite graph.  Specifically, for
        zero matrices of the appropriate size, for the reduced adjacency
        matrix ``H``, the full adjacency matrix is ``[[0, H'], [H, 0]]``.

        By default, the matrix returned is over the integers.

        INPUT:

        - ``sparse`` -- boolean (default: ``True``); whether to return a sparse
          matrix

        - ``base_ring`` -- a ring (default: ``None``); the base ring of the
          matrix space to use. By default, the base ring is ``ZZ`` if the graph
          is not weighted and otherwise the same ring as the (first) weights.

        - ``**kwds`` -- other keywords to pass to
          :func:`~sage.matrix.constructor.matrix`

        EXAMPLES:

        Bipartite graphs that are not weighted will return a matrix over ZZ,
        unless a base ring is specified::

            sage: # needs sage.modules
            sage: M = Matrix([(1,1,1,0,0,0,0), (1,0,0,1,1,0,0),
            ....:             (0,1,0,1,0,1,0), (1,1,0,1,0,0,1)])
            sage: B = BipartiteGraph(M)
            sage: N = B.reduced_adjacency_matrix(); N
            [1 1 1 0 0 0 0]
            [1 0 0 1 1 0 0]
            [0 1 0 1 0 1 0]
            [1 1 0 1 0 0 1]
            sage: N == M
            True
            sage: N[0,0].parent()
            Integer Ring
            sage: N2 = B.reduced_adjacency_matrix(base_ring=RDF); N2
            [1.0 1.0 1.0 0.0 0.0 0.0 0.0]
            [1.0 0.0 0.0 1.0 1.0 0.0 0.0]
            [0.0 1.0 0.0 1.0 0.0 1.0 0.0]
            [1.0 1.0 0.0 1.0 0.0 0.0 1.0]
            sage: N2[0, 0].parent()
            Real Double Field

        Multi-edge graphs also return a matrix over ZZ,
        unless a base ring is specified::

            sage: # needs sage.modules
            sage: M = Matrix([(1,1,2,0,0), (0,2,1,1,1), (0,1,2,1,1)])
            sage: B = BipartiteGraph(M, multiedges=True, sparse=True)
            sage: N = B.reduced_adjacency_matrix()
            sage: N == M
            True
            sage: N[0,0].parent()
            Integer Ring
            sage: N2 = B.reduced_adjacency_matrix(base_ring=RDF)
            sage: N2[0, 0].parent()
            Real Double Field

        Weighted graphs will return a matrix over the ring given by their
        (first) weights, unless a base ring is specified::

            sage: # needs sage.modules sage.rings.finite_rings
            sage: F.<a> = GF(4)
            sage: MS = MatrixSpace(F, 2, 3)
            sage: M = MS.matrix([[0, 1, a+1], [a, 1, 1]])
            sage: B = BipartiteGraph(M, weighted=True, sparse=True)
            sage: N = B.reduced_adjacency_matrix(sparse=False)
            sage: N == M
            True
            sage: N[0,0].parent()
            Finite Field in a of size 2^2
            sage: N2 = B.reduced_adjacency_matrix(base_ring=F)
            sage: N2[0, 0].parent()
            Finite Field in a of size 2^2

        TESTS::

            sage: B = BipartiteGraph()
            sage: B.reduced_adjacency_matrix()                                          # needs sage.modules
            []
            sage: M = Matrix([[0,0], [0,0]])                                            # needs sage.modules
            sage: BipartiteGraph(M).reduced_adjacency_matrix() == M                     # needs sage.modules
            True
            sage: M = Matrix([[10,2/3], [0,0]])                                         # needs sage.modules
            sage: B = BipartiteGraph(M, weighted=True, sparse=True)                     # needs sage.modules
            sage: M == B.reduced_adjacency_matrix()                                     # needs sage.modules
            True

        An error is raised if the specified base ring is not compatible with the
        type of the weights of the bipartite graph::

            sage: # needs sage.modules sage.rings.finite_rings
            sage: F.<a> = GF(4)
            sage: MS = MatrixSpace(F, 2, 3)
            sage: M = MS.matrix([[0, 1, a+1], [a, 1, 1]])
            sage: B = BipartiteGraph(M, weighted=True, sparse=True)
            sage: B.reduced_adjacency_matrix(base_ring=RDF)
            Traceback (most recent call last):
            ...
            TypeError: float() argument must be a string or a ...number, not 'sage.rings.finite_rings.element_givaro.FiniteField_givaroElement'
        """
        if self.multiple_edges() and self.weighted():
            raise NotImplementedError(
                "don't know how to represent weights for a multigraph")
        if self.is_directed():
            raise NotImplementedError(
                "reduced adjacency matrix does not exist for directed graphs")

        # Create mappings of left and right vertices to integers.
        # These mappings are used to translate an edge to its reduced adjacency
        # matrix position, that is: (row index, column index)
        left = {v: i for i, v in enumerate(sorted(self.left))}
        right = {v: i for i, v in enumerate(sorted(self.right))}

        # create dictionary of edges, values are weights for weighted graph,
        # otherwise the number of edges (always 1 for simple graphs)
        D = {}
        if self.weighted():
            for v1, v2, weight in self.edge_iterator():
                if v1 in left:
                    D[right[v2], left[v1]] = weight
                else:
                    D[right[v1], left[v2]] = weight
        else:
            # if we're normal or multi-edge, just create the matrix over ZZ
            for v1, v2, name in self.edge_iterator():
                idx = (right[v2], left[v1]) if v1 in left else (right[v1], left[v2])
                if idx in D:
                    D[idx] += 1
                else:
                    D[idx] = 1

        # now construct and return the matrix from the dictionary we created
        from sage.matrix.constructor import matrix
        if base_ring is None:
            return matrix(len(self.right), len(self.left), D, sparse=sparse, **kwds)
        return matrix(base_ring, len(self.right), len(self.left), D, sparse=sparse, **kwds)

    def matching(self, value_only=False, algorithm=None,
                 use_edge_labels=False, solver=None, verbose=0,
                 *, integrality_tolerance=1e-3):
        r"""
        Return a maximum matching of the graph represented by the list of its
        edges.

        Given a graph `G` such that each edge `e` has a weight `w_e`, a maximum
        matching is a subset `S` of the edges of `G` of maximum weight such that
        no two edges of `S` are incident with each other.

        INPUT:

        - ``value_only`` -- boolean (default: ``False``); when set to ``True``,
          only the cardinal (or the weight) of the matching is returned

        - ``algorithm`` -- string (default: ``'Hopcroft-Karp'`` if
          ``use_edge_labels==False``, otherwise ``'Edmonds'``); algorithm to use
          among:

          - ``'Hopcroft-Karp'`` selects the default bipartite graph algorithm as
            implemented in NetworkX

          - ``'Eppstein'`` selects Eppstein's algorithm as implemented in
            NetworkX

          - ``'Edmonds'`` selects Edmonds' algorithm as implemented in NetworkX

          - ``'LP'`` uses a Linear Program formulation of the matching problem

        - ``use_edge_labels`` -- boolean (default: ``False``)

          - when set to ``True``, computes a weighted matching where each edge
            is weighted by its label (if an edge has no label, `1` is assumed);
            only if ``algorithm`` is ``'Edmonds'``, ``'LP'``

          - when set to ``False``, each edge has weight `1`

        - ``solver`` -- string (default: ``None``); specifies a Mixed
          Integer Linear Programming (MILP) solver to be used. If set
          to ``None``, the default one is used. For more information
          on MILP solvers and which default solver is used, see the
          method :meth:`solve
          <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the
          class :class:`MixedIntegerLinearProgram
          <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: 0); sets the level of
          verbosity. Set to 0 by default, which means quiet.

        - ``integrality_tolerance`` -- float; parameter for use with
          MILP solvers over an inexact base ring; see
          :meth:`MixedIntegerLinearProgram.get_values`.

        .. SEEALSO::

            - :wikipedia:`Matching_(graph_theory)`
            - :meth:`~Graph.matching`

        EXAMPLES:

        Maximum matching in a cycle graph::

            sage: G = BipartiteGraph(graphs.CycleGraph(10))
            sage: G.matching()                                                          # needs networkx
            [(0, 1, None), (2, 3, None), (4, 5, None), (6, 7, None), (8, 9, None)]

        The size of a maximum matching in a complete bipartite graph using
        Eppstein::

            sage: G = BipartiteGraph(graphs.CompleteBipartiteGraph(4,5))
            sage: G.matching(algorithm='Eppstein', value_only=True)                     # needs networkx
            4

        TESTS:

        If ``algorithm`` is not set to one of the supported algorithms, an
        exception is raised::

            sage: G = BipartiteGraph(graphs.CompleteBipartiteGraph(4,5))
            sage: G.matching(algorithm='somethingdifferent')
            Traceback (most recent call last):
            ...
            ValueError: algorithm must be "Hopcroft-Karp", "Eppstein", "Edmonds" or "LP"

        Maximum matching in a weighted bipartite graph::

            sage: G = graphs.CycleGraph(4)
            sage: B = BipartiteGraph([(u,v,2) for u,v in G.edges(sort=True, labels=0)])
            sage: sorted(B.matching(use_edge_labels=True))                              # needs networkx
            [(0, 3, 2), (1, 2, 2)]
            sage: B.matching(use_edge_labels=True, value_only=True)                     # needs networkx
            4
            sage: B.matching(use_edge_labels=True, value_only=True, algorithm='Edmonds')            # needs networkx
            4
            sage: B.matching(use_edge_labels=True, value_only=True, algorithm='LP')     # needs sage.numerical.mip
            4
            sage: B.matching(use_edge_labels=True, value_only=True, algorithm='Eppstein')
            Traceback (most recent call last):
            ...
            ValueError: use_edge_labels cannot be used with "Hopcroft-Karp" or "Eppstein"
            sage: B.matching(use_edge_labels=True, value_only=True, algorithm='Hopcroft-Karp')
            Traceback (most recent call last):
            ...
            ValueError: use_edge_labels cannot be used with "Hopcroft-Karp" or "Eppstein"
            sage: B.matching(use_edge_labels=False, value_only=True,                    # needs networkx
            ....:            algorithm='Hopcroft-Karp')
            2
            sage: B.matching(use_edge_labels=False, value_only=True,                    # needs networkx
            ....:            algorithm='Eppstein')
            2
            sage: B.matching(use_edge_labels=False, value_only=True, algorithm='Edmonds')           # needs networkx
            2
            sage: B.matching(use_edge_labels=False, value_only=True, algorithm='LP')    # needs sage.numerical.mip
            2

        With multiedges enabled::

            sage: G = BipartiteGraph(graphs.CubeGraph(3))
            sage: for e in G.edges(sort=True):
            ....:     G.set_edge_label(e[0], e[1], int(e[0]) + int(e[1]))
            sage: G.allow_multiple_edges(True)
            sage: G.matching(use_edge_labels=True, value_only=True)                     # needs networkx
            444

        Empty bipartite graph and bipartite graphs without edges::

            sage: B = BipartiteGraph()
            sage: algorithms = ["Hopcroft-Karp", "Eppstein", "Edmonds", "LP"]
            sage: not any(B.matching(algorithm=algo) for algo in algorithms)            # needs networkx
            True
            sage: all(B.matching(algorithm=algo, value_only=True) == 0 for algo in algorithms)      # needs networkx
            True
            sage: B.add_vertex(1, left=True)
            sage: B.add_vertex(2, left=True)
            sage: B.add_vertex(3, right=True)
            sage: not any(B.matching(algorithm=algo) for algo in algorithms)            # needs networkx
            True
            sage: all(B.matching(algorithm=algo, value_only=True) == 0 for algo in algorithms)      # needs networkx
            True
        """
        if algorithm is None:
            algorithm = "Edmonds" if use_edge_labels else "Hopcroft-Karp"

        if algorithm == "Hopcroft-Karp" or algorithm == "Eppstein":
            if use_edge_labels:
                raise ValueError('use_edge_labels cannot be used with '
                                 '"Hopcroft-Karp" or "Eppstein"')
            d = []
            if self.size():
                import networkx
                # NetworkX matching algorithms for bipartite graphs may fail
                # when the graph is not connected
                if not self.is_connected():
                    CC = [g for g in self.connected_components_subgraphs() if g.size()]
                else:
                    CC = [self]
                v2int = {v: i for i, v in enumerate(self)}
                for g in CC:
                    h = g.networkx_graph()
                    if algorithm == "Hopcroft-Karp":
                        m = networkx.bipartite.hopcroft_karp_matching(h)
                    else:
                        m = networkx.bipartite.eppstein_matching(h)
                    d.extend((u, v, g.edge_label(u, v)) for u, v in m.items()
                             if v2int[u] < v2int[v])

            if value_only:
                return Integer(len(d))
            return d

        elif algorithm == "Edmonds" or algorithm == "LP":
            return Graph.matching(self, value_only=value_only,
                                  algorithm=algorithm,
                                  use_edge_labels=use_edge_labels,
                                  solver=solver, verbose=verbose,
                                  integrality_tolerance=integrality_tolerance)
        raise ValueError('algorithm must be "Hopcroft-Karp", '
                         '"Eppstein", "Edmonds" or "LP"')

    def vertex_cover(self, algorithm='Konig', value_only=False,
                     reduction_rules=True, solver=None, verbose=0,
                     *, integrality_tolerance=1e-3):
        r"""
        Return a minimum vertex cover of ``self`` represented by a set of vertices.

        A minimum vertex cover of a graph is a set `S` of vertices such that
        each edge is incident to at least one element of `S`, and such that `S`
        is of minimum cardinality. For more information, see
        :wikipedia:`Vertex_cover`.

        Equivalently, a vertex cover is defined as the complement of an
        independent set.

        As an optimization problem, it can be expressed as follows:

        .. MATH::

            \mbox{Minimize : }&\sum_{v\in G} b_v\\
            \mbox{Such that : }&\forall (u,v) \in G.edges(sort=True), b_u+b_v\geq 1\\
            &\forall x\in G, b_x\mbox{ is a binary variable}

        INPUT:

        - ``algorithm`` -- string (default: ``'Konig'``); algorithm to use
          among:

          - ``'Konig'`` will compute a minimum vertex cover using Konig's
            algorithm (:wikipedia:`Kőnig%27s_theorem_(graph_theory)`)

          - ``'Cliquer'`` will compute a minimum vertex cover
            using the Cliquer package

          - ``'MILP'`` will compute a minimum vertex cover through a mixed
            integer linear program

          - ``'mcqd'`` will use the MCQD solver
            (`<http://www.sicmm.org/~konc/maxclique/>`_), and the MCQD
            package must be installed

        - ``value_only`` -- boolean (default: ``False``); if set to ``True``,
          only the size of a minimum vertex cover is returned. Otherwise,
          a minimum vertex cover is returned as a list of vertices.

        - ``reduction_rules`` -- (default: ``True``) specify if the reductions
          rules from kernelization must be applied as pre-processing or not.
          See [ACFLSS04]_ for more details. Note that depending on the instance,
          it might be faster to disable reduction rules.  This parameter is
          currently ignored when ``algorithm == "Konig"``.

        - ``solver`` -- string (default: ``None``); specifies a Mixed Integer
          Linear Programming (MILP) solver to be used. If set to ``None``, the
          default one is used. For more information on MILP solvers and which
          default solver is used, see the method :meth:`solve
          <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
          :class:`MixedIntegerLinearProgram
          <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: 0); sets the level of
          verbosity. Set to 0 by default, which means quiet.

        - ``integrality_tolerance`` -- float; parameter for use with MILP
          solvers over an inexact base ring; see
          :meth:`MixedIntegerLinearProgram.get_values`.

        EXAMPLES:

        On the Cycle Graph::

            sage: B = BipartiteGraph(graphs.CycleGraph(6))
            sage: len(B.vertex_cover())                                                 # needs networkx
            3
            sage: B.vertex_cover(value_only=True)                                       # needs networkx
            3

        The two algorithms should return the same result::

           sage: # needs networkx numpy
           sage: g = BipartiteGraph(graphs.RandomBipartite(10, 10, .5))
           sage: vc1 = g.vertex_cover(algorithm='Konig')
           sage: vc2 = g.vertex_cover(algorithm='Cliquer')
           sage: len(vc1) == len(vc2)
           True

        TESTS:

        Giving a non connected bipartite graph::

            sage: B = BipartiteGraph(graphs.CycleGraph(4) * 2)
            sage: len(B.vertex_cover())                                                 # needs networkx
            4

        Empty bipartite graph and bipartite graphs without edges::

            sage: B = BipartiteGraph()
            sage: algorithms = ["Konig", "Cliquer", "MILP"]
            sage: all(B.vertex_cover(algorithm=algo) == [] for algo in algorithms)
            True
            sage: all(B.vertex_cover(algorithm=algo, value_only=True) == 0 for algo in algorithms)
            True
            sage: B.add_vertex(1, left=True)
            sage: B.add_vertex(2, left=True)
            sage: B.add_vertex(3, right=True)
            sage: all(B.vertex_cover(algorithm=algo) == [] for algo in algorithms)
            True
            sage: all(B.vertex_cover(algorithm=algo, value_only=True) == 0 for algo in algorithms)
            True
        """
        if algorithm != "Konig":
            return Graph.vertex_cover(self, algorithm=algorithm,
                                      value_only=value_only,
                                      reduction_rules=reduction_rules,
                                      solver=solver,
                                      verbose=verbose,
                                      integrality_tolerance=integrality_tolerance)

        if not self.is_connected():
            VC = []
            for b in self.connected_components_subgraphs():
                if b.size():
                    VC.extend(b.vertex_cover(algorithm='Konig'))
            if value_only:
                return sum(VC)
            return VC

        M = Graph(self.matching())
        left = set(self.left)
        right = set(self.right)

        # Initialize Z with vertices in left that are not involved in the
        # matching
        Z = left.difference(M.vertex_iterator())

        # Alternate: extend Z with all vertices reachable by alternate paths
        # (match / non-match edges).
        X = set(Z)
        while X:
            # Follow non matched edges
            Y = set()
            for u in X:
                for v in self.neighbors(u):
                    if v not in Z and not M.has_edge(u, v):
                        Y.add(v)
            Z.update(Y)

            # Follow matched edges
            X = set()
            for u in Y:
                for v in M.neighbor_iterator(u):
                    if v not in Z:
                        X.add(v)
            Z.update(X)

        # The solution is (left \ Z) + (right \cap Z)
        VC = list((left.difference(Z)).union(right.intersection(Z)))

        if value_only:
            return len(VC)
        return VC

    def _subgraph_by_adding(self, vertices=None, edges=None, edge_property=None, immutable=None):
        r"""
        Return the subgraph containing the given vertices and edges.

        The edges also satisfy the edge_property, if it is not None. The
        subgraph is created by creating a new empty graph and adding the
        necessary vertices, edges, and other properties.

        INPUT:

        - ``vertices`` -- list (default: ``None``); list of vertices

        - ``edges`` -- (default: ``None``) either a single edge or an iterable
          container of edges (e.g., a list, set, file, numeric array, etc.). If
          no edges are specified, then all edges are assumed and the returned
          graph is an induced subgraph. In the case of multiple edges,
          specifying an edge as `(u, v)` means to keep all edges `(u, v)`,
          regardless of the label.

        - ``edge_property`` -- (default: ``None``) if specified, this is
          expected to be a function on edges, which is intersected with the
          edges specified, if any are

        - ``immutable`` -- boolean (default: ``None``); currently ignored for
          ``BipartiteGraph``

        EXAMPLES::

            sage: B = BipartiteGraph(graphs.CycleGraph(6))
            sage: H = B._subgraph_by_adding(vertices=B.left)
            sage: H.order(), H.size()
            (3, 0)
            sage: H = B._subgraph_by_adding(vertices=[0, 1])
            sage: H.order(), H.size()
            (2, 1)
            sage: H = B._subgraph_by_adding(vertices=[0, 1, 2], edges=[(0, 1)])
            sage: H.order(), H.size()
            (3, 1)

        Using the property arguments::

            sage: B = BipartiteGraph([(0, 1, 1), (0, 2, 0), (0, 3, 0), (3, 4, 1)])
            sage: H = B._subgraph_by_adding(vertices=B.vertices(sort=False), edge_property=(lambda e: e[2] == 1))
            sage: H.order(), H.size()
            (5, 2)
        """
        B = self.__class__(weighted=self._weighted, loops=self.allows_loops(),
                           multiedges=self.allows_multiple_edges())
        B.name("Subgraph of ({})".format(self.name()))
        for u in vertices:
            if u in self.left:
                B.add_vertex(u, left=True)
            elif u in self.right:
                B.add_vertex(u, right=True)
        if edges is None:
            edges_to_keep = self.edge_boundary(B.left, B.right, sort=False)
        else:
            edges_to_keep_labeled = set(e for e in edges if len(e) == 3)
            edges_to_keep_unlabeled = set(e for e in edges if len(e) == 2)
            edges_to_keep = []
            for u, v, l in self.edge_boundary(B.left, B.right, sort=False):
                if ((u, v, l) in edges_to_keep_labeled
                        or (v, u, l) in edges_to_keep_labeled
                        or (u, v) in edges_to_keep_unlabeled
                        or (v, u) in edges_to_keep_unlabeled):
                    edges_to_keep.append((u, v, l))

        if edge_property is not None:
            edges_to_keep = [e for e in edges_to_keep if edge_property(e)]

        B.add_edges(edges_to_keep)

        attributes_to_update = ('_pos', '_assoc')
        for attr in attributes_to_update:
            if hasattr(self, attr) and getattr(self, attr) is not None:
                d = getattr(self, attr)
                value = {v: d.get(v, None) for v in B}
                setattr(B, attr, value)

        return B

    def _subgraph_by_deleting(self, vertices=None, edges=None, inplace=False,
                              edge_property=None, immutable=None):
        r"""
        Return the subgraph containing the given vertices and edges.

        The edges also satisfy the edge_property, if it is not None. The
        subgraph is created by creating deleting things that are not needed.

        INPUT:

        - ``vertices`` -- list (default: ``None``); list of vertices

        - ``edges`` -- (default: ``None``) either a single edge or an iterable
          container of edges (e.g., a list, set, file, numeric array, etc.). If
          no edges are specified, then all edges are assumed and the returned
          graph is an induced subgraph. In the case of multiple edges,
          specifying an edge as `(u, v)` means to keep all edges `(u, v)`,
          regardless of the label.

        - ``inplace`` -- boolean (default: ``False``); when ``True``, the
          current graph is modified in place by deleting the extra vertices and
          edges. Otherwise a modified copy of the graph is returned

        - ``edge_property`` -- (default: ``None``) if specified, this is
          expected to be a function on edges, which is intersected with the
          edges specified, if any are

        - ``immutable`` -- boolean (default: ``None``); currently ignored for
          ``BipartiteGraph``

        EXAMPLES::

            sage: B = BipartiteGraph(graphs.CycleGraph(6))
            sage: H = B._subgraph_by_deleting(vertices=B.left)
            sage: H.order(), H.size()
            (3, 0)
            sage: H = B._subgraph_by_deleting(vertices=[0, 1])
            sage: H.order(), H.size()
            (2, 1)
            sage: H = B._subgraph_by_deleting(vertices=[0, 1, 2], edges=[(0, 1)])
            sage: H.order(), H.size()
            (3, 1)

        Using the property arguments::

            sage: B = BipartiteGraph([(0, 1, 1), (0, 2, 0), (0, 3, 0), (3, 4, 1)])
            sage: H = B._subgraph_by_deleting(vertices=B.vertices(sort=False), edge_property=(lambda e: e[2] == 1))
            sage: H.order(), H.size()
            (5, 2)
        """
        if inplace:
            B = self
        else:
            # We make a copy of the graph
            B = BipartiteGraph(data=self.edges(sort=True), partition=[self.left, self.right])
            attributes_to_update = ('_pos', '_assoc')
            for attr in attributes_to_update:
                if hasattr(self, attr) and getattr(self, attr) is not None:
                    d = getattr(self, attr)
                    value = {v: d.get(v, None) for v in B}
                    setattr(B, attr, value)
        B.name("Subgraph of ({})".format(self.name()))

        vertices = set(vertices)
        B.delete_vertices([v for v in B.vertex_iterator() if v not in vertices])

        edges_to_delete = []
        if edges is not None:
            edges_to_keep_labeled = set(e for e in edges if len(e) == 3)
            edges_to_keep_unlabeled = set(e for e in edges if len(e) == 2)
            for u, v, l in B.edge_iterator():
                if ((u, v, l) not in edges_to_keep_labeled
                        and (v, u, l) not in edges_to_keep_labeled
                        and (u, v) not in edges_to_keep_unlabeled
                        and (v, u) not in edges_to_keep_unlabeled):
                    edges_to_delete.append((u, v, l))
        if edge_property is not None:
            # We might get duplicate edges, but this does handle the case of
            # multiple edges.
            edges_to_delete.extend(e for e in B.edge_iterator() if not edge_property(e))

        B.delete_edges(edges_to_delete)
        return B

    def canonical_label(self, partition=None, certificate=False,
                        edge_labels=False, algorithm=None, return_graph=True):
        r"""
        Return the canonical graph.

        A canonical graph is the representative graph of an isomorphism
        class by some canonization function `c`. If `G` and `H` are graphs,
        then `G \cong c(G)`, and `c(G) == c(H)` if and only if `G \cong H`.

        See the :wikipedia:`Graph_canonization` for more information.

        INPUT:

        - ``partition`` -- if given, the canonical label with respect
          to this set partition will be computed. The default is the unit
          set partition.

        - ``certificate`` -- boolean (default: ``False``); when set to
          ``True``, a dictionary mapping from the vertices of the (di)graph
          to its canonical label will also be returned

        - ``edge_labels`` -- boolean (default: ``False``); when set to
          ``True``, allows only permutations respecting edge labels

        - ``algorithm`` -- string (default: ``None``); the algorithm to use.
          Currently available:

          * ``'bliss'``: use the optional package bliss
            (http://www.tcs.tkk.fi/Software/bliss/index.html);
          * ``'sage'``: always use Sage's implementation.
          * ``None`` (default): use bliss when available and possible

            .. NOTE::

                Make sure you always compare canonical forms obtained by the
                same algorithm.

        - ``return_graph`` -- boolean (default: ``True``); when set to
          ``False``, returns the list of edges of the canonical graph
          instead of the canonical graph. Only available when ``'bliss'``
          is explicitly set as algorithm.

        EXAMPLES::

            sage: B = BipartiteGraph( [(0, 4), (0, 5), (0, 6), (0, 8), (1, 5),
            ....:                      (1, 7), (1, 8), (2, 6), (2, 7), (2, 8),
            ....:                      (3, 4), (3, 7), (3, 8), (4, 9), (5, 9),
            ....:                      (6, 9), (7, 9)] )
            sage: C = B.canonical_label(partition=(B.left,B.right), algorithm='sage')
            sage: C
            Bipartite graph on 10 vertices
            sage: C.left
            {0, 1, 2, 3, 4}
            sage: C.right
            {5, 6, 7, 8, 9}

        ::

            sage: B = BipartiteGraph( [(0, 4), (0, 5), (0, 6), (0, 8), (1, 5),
            ....:                      (1, 7), (1, 8), (2, 6), (2, 7), (2, 8),
            ....:                      (3, 4), (3, 7), (3, 8), (4, 9), (5, 9),
            ....:                      (6, 9), (7, 9)] )
            sage: C, cert = B.canonical_label(partition=(B.left, B.right),
            ....:                             certificate=True, algorithm='sage')
            sage: C
            Bipartite graph on 10 vertices
            sage: C.left
            {0, 1, 2, 3, 4}
            sage: C.right
            {5, 6, 7, 8, 9}
            sage: cert == {0: 3, 1: 0, 2: 1, 3: 2, 4: 5, 5: 7, 6: 6, 7: 8, 8: 9, 9: 4}
            True

        ::

            sage: G = Graph({0: [5, 6], 1: [4, 5], 2: [4, 6], 3: [4, 5, 6]})
            sage: B = BipartiteGraph(G)
            sage: C = B.canonical_label(partition=(B.left, B.right),
            ....:                       edge_labels=True, algorithm='sage')
            sage: C.left
            {0, 1, 2, 3}
            sage: C.right
            {4, 5, 6}

        TESTS:

        Check that :issue:`38832` is fixed::

            sage: B = BipartiteGraph(matrix([[1, 1], [1, 1]]))
            sage: C = B.canonical_label()
            sage: C.left, C.right
            ({0, 1}, {2, 3})
            sage: C, certificate = B.canonical_label(certificate=True)
            sage: C.left, C.right
            ({0, 1}, {2, 3})
            sage: C = B.canonical_label(edge_labels=True)
            sage: C.left, C.right
            ({0, 1}, {2, 3})
            sage: B.allow_multiple_edges(True)
            sage: B.add_edges(B.edges())
            sage: C = B.canonical_label()
            sage: C.left, C.right
            ({0, 1}, {2, 3})

        .. SEEALSO::

            :meth:`~sage.graphs.generic_graph.GenericGraph.canonical_label()`
        """
        if certificate:
            C, cert = GenericGraph.canonical_label(self, partition=partition,
                                                   certificate=certificate,
                                                   edge_labels=edge_labels,
                                                   algorithm=algorithm,
                                                   return_graph=return_graph)

        else:
            from sage.groups.perm_gps.partn_ref.refinement_graphs import search_tree
            from sage.graphs.graph import Graph
            from sage.graphs.generic_graph import graph_isom_equivalent_non_edge_labeled_graph
            from itertools import chain

            cert = {}

            if edge_labels or self.has_multiple_edges():
                G, partition, relabeling = graph_isom_equivalent_non_edge_labeled_graph(self, partition, return_relabeling=True)
                G_vertices = list(chain(*partition))
                G_to = {u: i for i, u in enumerate(G_vertices)}
                H = Graph(len(G_vertices))
                HB = H._backend
                for u, v in G.edge_iterator(labels=False):
                    HB.add_edge(G_to[u], G_to[v], None, False)
                GC = HB.c_graph()[0]
                partition = [[G_to[vv] for vv in cell] for cell in partition]
                a, b, c = search_tree(GC, partition, certificate=True, dig=False)
                # c is a permutation to the canonical label of G,
                # which depends only on isomorphism class of self.
                cert = {v: c[G_to[relabeling[v]]] for v in self}

            else:
                if partition is None:
                    partition = self.bipartition()
                G_vertices = list(chain(*partition))
                G_to = {u: i for i, u in enumerate(G_vertices)}
                H = Graph(len(G_vertices))
                HB = H._backend
                for u, v in self.edge_iterator(labels=False):
                    HB.add_edge(G_to[u], G_to[v], None, False)
                GC = HB.c_graph()[0]
                partition = [[G_to[vv] for vv in cell] for cell in partition]
                a, b, c = search_tree(GC, partition, certificate=True, dig=False)
                cert = {v: c[G_to[v]] for v in G_to}

            C = self.relabel(perm=cert, inplace=False)

        C.left = {cert[v] for v in self.left}
        C.right = {cert[v] for v in self.right}

        if certificate:
            return C, cert
        return C
