# distutils: language = c++
# distutils: extra_compile_args = -std=c++11
r"""
Split-tree (Cunningham) decomposition of connected graphs.

This module implements the incremental split decomposition algorithm described
in [GPTC2014]_. The output is the *split-tree* of the input graph, i.e., the
unique reduced graph-labelled tree (GLT) whose accessibility graph is the input
graph (Cunningham's canonical split decomposition, reformulated in terms of
GLTs).

The algorithm inserts the vertices one at a time, following a LexBFS ordering.
For each inserted vertex `x`, the current split-tree is marked with respect to
the set `S` of already inserted neighbours of `x`; exactly one of the seven
cases of Theorem 4.1 of [GPTC2014]_ applies, and the tree is updated by the
local construction of the corresponding proposition.

This implementation slightly deviates from the quasi-linear time construction
of the data structure described in the paper.

EXAMPLES::

    sage: from sage.graphs.graph_decompositions.split_decomposition import split_tree
    sage: from sage.graphs.graph_decompositions.split_decomposition import GP2012
    sage: G = GP2012()
    sage: ST = split_tree(G); ST
    Split tree on 16 vertices with 6 internal nodes
    sage: ST.root()
    10
    sage: ST.nodes()
    [N0, N1, N2, N3, N4, N5]
    sage: sorted(ST.kind(u) for u in ST.nodes())
    ['CLIQUE', 'CLIQUE', 'PRIME', 'PRIME', 'PRIME', 'STAR']
    sage: ST.leaves()
    [10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16]
    sage: T = ST.tree()
    sage: T.is_tree()
    True
    sage: sorted(T, key=str)
    [1, 10, 11, 12, 13, 14, 15, 16, 2, 3, 4, 5, 6, 7, 8, 9, N0, N1, N2, N3, N4, N5]
    sage: A = ST.accessibility_graph()
    sage: A.is_isomorphic(G)
    True
    sage: A == G
    True
    sage: sorted([(ST.kind(u), ST.label_graph(u)) for u in ST.nodes()], key=str)
    [('CLIQUE', Graph on 3 vertices),
     ('CLIQUE', Graph on 4 vertices),
     ('PRIME', Graph on 5 vertices),
     ('PRIME', Graph on 5 vertices),
     ('PRIME', Graph on 5 vertices),
     ('STAR', Graph on 4 vertices)]
"""
# ******************************************************************************
#       Copyright (C) 2026 David Coudert <david.coudert@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ******************************************************************************

from cysignals.signals cimport sig_on, sig_off
from libcpp.vector cimport vector
from libcpp.string cimport string

from collections import namedtuple


# ------------------------------------------------------------------------------
# Enumeration class and types
# ------------------------------------------------------------------------------

SplitTreeNode = namedtuple('SplitTreeNode', ['index'])
SplitTreeNode.__repr__ = lambda self: f"N{self.index}"
SplitTreeNode.__str__ = SplitTreeNode.__repr__

_KIND_NAMES = ('PRIME', 'CLIQUE', 'STAR')


# ------------------------------------------------------------------------------
# split tree
# ------------------------------------------------------------------------------

cdef class SplitTree:
    r"""
    The split tree (reduced graph-labelled tree) of a connected graph.

    This class implements the incremental split decomposition algorithm
    described in [GPTC2014]_. The output is the *split-tree* of the input graph,
    i.e., the unique reduced graph-labelled tree (GLT) whose accessibility graph
    is the input graph (Cunningham's canonical split decomposition, reformulated
    in terms of GLTs).

    .. WARNING::

        The split-tree is a data structure that is not meant to be modified
        after its construction. The methods of this class are only accessors
        and queries on the split-tree.

        Users should not attempt to build a split-tree directly. It should be
        built via the
        :func:`~sage.graphs.graph_decompositions.split_decomposition.split_tree`
        function, which implements the incremental construction of the
        split-tree and returns an instance of this class.

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.split_decomposition import split_tree
        sage: from sage.graphs.graph_decompositions.split_decomposition import SplitTree
        sage: from sage.graphs.graph_decompositions.split_decomposition import GP2012
        sage: G = GP2012()
        sage: ST = split_tree(G); ST
        Split tree on 16 vertices with 6 internal nodes
        sage: isinstance(ST, SplitTree)
        True
        sage: ST.root()
        10
        sage: ST.nodes()
        [N0, N1, N2, N3, N4, N5]
        sage: sorted(ST.kind(u) for u in ST.nodes())
        ['CLIQUE', 'CLIQUE', 'PRIME', 'PRIME', 'PRIME', 'STAR']
        sage: ST.leaves()
        [10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16]
        sage: set(ST.leaves()) == set(G)
        True
        sage: T = ST.tree()
        sage: T.is_tree()
        True
        sage: sorted(T, key=str)
        [1, 10, 11, 12, 13, 14, 15, 16, 2, 3, 4, 5, 6, 7, 8, 9, N0, N1, N2, N3, N4, N5]
        sage: A = ST.accessibility_graph()
        sage: A.is_isomorphic(G)
        True
        sage: A == G
        True
        sage: sorted([(ST.kind(u), ST.label_graph(u)) for u in ST.nodes()], key=str)
        [('CLIQUE', Graph on 3 vertices),
         ('CLIQUE', Graph on 4 vertices),
         ('PRIME', Graph on 5 vertices),
         ('PRIME', Graph on 5 vertices),
         ('PRIME', Graph on 5 vertices),
         ('STAR', Graph on 4 vertices)]

    TESTS::

        sage: split_tree(Graph())
        Split tree on 0 vertices with 0 internal nodes
        sage: split_tree(Graph(1))
        Split tree on 1 vertices with 0 internal nodes
        sage: split_tree(Graph(2))
        Traceback (most recent call last):
        ...
        ValueError: the graph must be connected
        sage: split_tree(Graph([(0, 1)]))
        Split tree on 2 vertices with 0 internal nodes
        sage: split_tree(DiGraph())
        Traceback (most recent call last):
        ...
        ValueError: the graph must be undirected
        sage: split_tree(Graph([(0, 0)], loops=True))
        Traceback (most recent call last):
        ...
        ValueError: the graph must be without loops
        sage: split_tree(Graph([(0, 1), (0, 1)], multiedges=True))
        Split tree on 2 vertices with 0 internal nodes
        sage: split_tree(Graph([(0, 1), (0, 1)], multiedges=True)).accessibility_graph()
        Graph on 2 vertices
    """
    cdef object _root  # the leaf (graph vertex) at which the tree is rooted
    cdef list _labels  # ordered list of vertices of the graph
    cdef dict _leaf_node  # dictionary mapping each leaf to its node in the GLT
    cdef list _info  # list of dicts with information on each node of the GLT
    cdef list _edges  # list of edges of the GLT

    def __init__(self, list labels, dict exported):
        r"""
        Initialize data structures.

        This method build the internal representation of the graph-labelled
        tree.

        INPUT:

        - ``labels`` -- a list of vertices of the graph, in the order they are
          inserted

        - ``exported`` -- a dictionary containing the exported information
          from the C++ implementation of the split-tree construction

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.split_decomposition import split_tree
            sage: G = graphs.DiamondGraph()
            sage: ST = split_tree(G, check=True); ST
            Split tree on 4 vertices with 2 internal nodes
            sage: T = ST.tree()
            sage: T.is_tree()
            True
            sage: sorted(T, key=str)
            [0, 1, 2, 3, N0, N1]
            sage: H = ST.accessibility_graph()
            sage: H.is_isomorphic(G)
            True

        The construction is correct with and without internal validation tests::

            sage: from sage.graphs.graph_decompositions.split_decomposition import GP2012
            sage: G = GP2012()
            sage: ST = split_tree(G, check=True)
            sage: ST.accessibility_graph().is_isomorphic(G)
            True
            sage: ST = split_tree(G, check=False)
            sage: ST.accessibility_graph().is_isomorphic(G)
            True

        TESTS::

            sage: split_tree(Graph()).tree()
            Graph on 0 vertices
            sage: split_tree(Graph(1)).tree()
            Graph on 1 vertex
            sage: split_tree(Graph(2))
            Traceback (most recent call last):
            ...
            ValueError: the graph must be connected
            sage: split_tree(DiGraph())
            Traceback (most recent call last):
            ...
            ValueError: the graph must be undirected
            sage: split_tree(Graph([(0, 0)], loops=True))
            Traceback (most recent call last):
            ...
            ValueError: the graph must be without loops
            sage: split_tree(Graph([(0, 1), (0, 1)], multiedges=True)).tree()
            Graph on 2 vertices
        """
        self._labels = list(labels)
        self._leaf_node = {}
        self._info = []
        self._edges = []
        if len(labels) < 3:  # fewer than 3 vertices: no internal node
            self._root = self._labels[0] if self._labels else None
            if len(self._labels) == 2:
                self._edges.append(tuple(self._labels))
            return

        def name(c):
            return self._labels[c] if c >= 0 else SplitTreeNode(-c - 1)

        self._root = self._labels[exported['root_vertex']]
        cdef vector[int] ms = exported['node_marker_start']
        cdef vector[int] mn = exported['node_marker_neighbour']
        cdef vector[int] es = exported['node_edge_start']
        cdef vector[int] ep = exported['node_edge_pairs']
        cdef int j, c
        for j, kind in enumerate(exported['node_kind']):
            markers = [name(c) for c in mn[ms[j]:ms[j + 1]]]
            prime_edges = [(name(ep[2 * e]), name(ep[2 * e + 1]))
                           for e in range(es[j], es[j + 1])]
            centre = name(exported['node_star_centre'][j]) if kind == 2 else None
            self._info.append({'kind': _KIND_NAMES[kind],
                               'parent': name(exported['node_parent'][j]),
                               'markers': markers,
                               'centre': centre,
                               'edges': prime_edges,
                               'edge_set': set(frozenset(e) for e in prime_edges)})
            for m in markers:
                if not isinstance(m, SplitTreeNode):
                    self._leaf_node[m] = j

    # --------------------------------------------------------------------------
    # basic accessors
    # --------------------------------------------------------------------------

    def root(self):
        r"""
        Return the leaf (graph vertex) at which the tree is rooted.
        
        The root is the first vertex of the LexBFS order used to build the
        split tree.

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.split_decomposition import split_tree
            sage: G = graphs.DiamondGraph()
            sage: ST = split_tree(G, check=True)
            sage: ST.root()
            0
        """
        return self._root

    def nodes(self):
        r"""
        Return the list of internal nodes, in preorder from the root.

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.split_decomposition import split_tree
            sage: G = graphs.DiamondGraph()
            sage: ST = split_tree(G, check=True)
            sage: ST.nodes()
            [N0, N1]
        """
        return [SplitTreeNode(j) for j in range(len(self._info))]

    def leaves(self):
        r"""
        Return the leaves of the tree, i.e. the vertices of the graph.

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.split_decomposition import split_tree
            sage: G = graphs.DiamondGraph()
            sage: ST = split_tree(G, check=True)
            sage: ST.leaves()
            [0, 1, 2, 3]
        """
        return list(self._labels)

    def kind(self, node):
        r"""
        Return the kind of the node.
        
        The kind of a node is one of ``'PRIME'``, ``'CLIQUE'`` or ``'STAR'``.

        INPUT:

        - ``node`` -- a node of the split tree

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.split_decomposition import split_tree
            sage: G = graphs.DiamondGraph()
            sage: ST = split_tree(G, check=True)
            sage: sorted([ST.kind(n) for n in ST.nodes()])
            ['CLIQUE', 'STAR']
        """
        return self._info[node.index]['kind']

    def parent(self, node):
        r"""
        Return the parent of ``node``, i.e. the tree-neighbour towards the root.

        INPUT:

        - ``node`` -- a node of the split tree

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.split_decomposition import split_tree
            sage: G = graphs.DiamondGraph()
            sage: ST = split_tree(G, check=True)
            sage: [ST.parent(n) for n in ST.nodes()]
            [0, N0]
        """
        return self._info[node.index]['parent']

    def markers(self, node):
        r"""
        Return the markers of ``node``.
        
        The markers are named by the tree-neighbour they face (the parent included).
    
        INPUT:

        - ``node`` -- a node of the split tree

        EXAMPLES::
        
            sage: from sage.graphs.graph_decompositions.split_decomposition import split_tree
            sage: G = graphs.DiamondGraph()
            sage: ST = split_tree(G, check=True)
            sage: [sorted(ST.markers(n), key=str) for n in ST.nodes()]
            [[0, 3, N1], [1, 2, N0]]
        """
        return list(self._info[node.index]['markers'])

    def children(self, node):
        r"""
        Return the tree-neighbours of ``node`` other than its parent.

        INPUT:

        - ``node`` -- a node of the split tree

        EXAMPLES::
        
            sage: from sage.graphs.graph_decompositions.split_decomposition import split_tree
            sage: G = graphs.DiamondGraph()
            sage: ST = split_tree(G, check=True)
            sage: [sorted(ST.children(n), key=str) for n in ST.nodes()]
            [[3, N1], [1, 2]]
        """
        p = self._info[node.index]['parent']
        return [m for m in self._info[node.index]['markers'] if m != p]

    def star_centre(self, node):
        r"""
        Return the centre of a ``'STAR'`` node.

        For a ``'STAR'`` node, the marker (tree-neighbour) that is the centre
        of the star. Otherwise, raises an error.

        INPUT:

        - ``node`` -- a node of the split tree

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.split_decomposition import split_tree
            sage: G = graphs.DiamondGraph()
            sage: ST = split_tree(G, check=True)
            sage: [ST.star_centre(n) for n in ST.nodes() if ST.kind(n) == 'STAR']
            [N1]
        """
        if self.kind(node) != 'STAR':
            raise ValueError("node is not a star")
        return self._info[node.index]['centre']

    def prime_edges(self, node):
        r"""
        Return the edges of a ``'PRIME'`` node.

        For a ``'PRIME'`` node, this method returns the edges of its label
        graph as pairs of markers.

        INPUT:

        - ``node`` -- a node of the split tree

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.split_decomposition import split_tree
            sage: G = graphs.CycleGraph(5)
            sage: G.add_clique([4, 5, 6, 7])
            sage: ST = split_tree(G, check=True); ST
            Split tree on 8 vertices with 3 internal nodes
            sage: sorted([ST.kind(n) for n in ST.nodes()])
            ['CLIQUE', 'PRIME', 'STAR']
            sage: p = [u for u in ST.nodes() if ST.kind(u) == 'PRIME'][0]
            sage: sorted([sorted(e, key=str) for e in ST.prime_edges(p)], key=str)
            [[0, 1], [0, N1], [1, 2], [2, 3], [3, N1]]
        """
        if self.kind(node) != 'PRIME':
            raise ValueError("node is not prime")
        return list(self._info[node.index]['edges'])

    # --------------------------------------------------------------------------
    # derived graphs
    # --------------------------------------------------------------------------

    def tree(self, immutable=False):
        r"""
        Return the underlying tree as a Sage graph.
        
        The leaves are the graph vertices. The internal nodes are the split
        tree nodes.

        INPUT:

        - ``immutable`` -- boolean (default: ``False``); whether to create a
          mutable/immutable graph.

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.split_decomposition import split_tree
            sage: G = graphs.DiamondGraph()
            sage: ST = split_tree(G, check=True)
            sage: T = ST.tree()
            sage: T.is_tree()
            True
            sage: sorted(T, key=str)
            [0, 1, 2, 3, N0, N1]
        """
        from itertools import chain
        from sage.graphs.graph import Graph

        vertices = chain(self._labels, self.nodes())
        edges = chain(self._edges,
                      ((SplitTreeNode(j), m)
                       for j, info in enumerate(self._info)
                       for m in info['markers']))
        return Graph([vertices, edges], format="vertices_and_edges",
                     immutable=immutable)

    def label_graph(self, node, immutable=False):
        r"""
        Return the label graph of ``node`` as a Sage graph.
        
        The vertices of the graph are the markers.

        INPUT:

        - ``node`` -- a node of the split tree

        - ``immutable`` -- boolean (default: ``False``); whether to create a
          mutable/immutable graph.

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.split_decomposition import split_tree
            sage: G = graphs.DiamondGraph()
            sage: ST = split_tree(G, check=True)
            sage: [ST.label_graph(n) for n in ST.nodes()]
            [Graph on 3 vertices,
             Graph on 3 vertices]
        """
        from sage.graphs.graph import Graph
        info = self._info[node.index]
        if info['kind'] == 'CLIQUE':
            from itertools import combinations
            edges = combinations(info['markers'], 2)
        elif info['kind'] == 'STAR':
            c = info['centre']
            edges = ((c, m) for m in info['markers'] if m != c)
        else:
            edges = info['edges']
        return Graph([info['markers'], edges], format="vertices_and_edges",
                     immutable=immutable)

    cdef bint _label_adjacent(self, j, a, b):
        r"""
        Return whether two markers of a node are adjacent in the label graph.

        INPUT:

        - ``j`` -- index of the node
        - ``a``, ``b`` -- markers of the node

        OUTPUT:

        - a boolean
        """
        info = self._info[j]
        if a == b:
            return False
        if info['kind'] == 'CLIQUE':
            return True
        if info['kind'] == 'STAR':
            return a == info['centre'] or b == info['centre']
        return frozenset((a, b)) in info['edge_set']

    def accessibility_graph(self, immutable=False):
        r"""
        Return the accessibility graph of the split tree.
        
        In the accessibility graph, two leaves are adjacent iff the tree path
        between them alternates tree-edges and label-edges. For the tree built
        by :func:`split_tree` it equals the input graph (a good sanity check).

        INPUT:

        - ``immutable`` -- boolean (default: ``False``); whether to create a
          mutable/immutable graph.

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.split_decomposition import split_tree
            sage: G = graphs.DiamondGraph()
            sage: ST = split_tree(G, check=True)
            sage: ST.accessibility_graph()
            Graph on 4 vertices
            sage: ST.accessibility_graph().is_isomorphic(G)
            True
        """
        from sage.graphs.graph import Graph
        cdef list edges = list(self._edges)
        for x, j0 in self._leaf_node.items():
            stack = [(j0, x)]     # (node, marker through which we entered it)
            while stack:
                j, came = stack.pop()
                for to in self._info[j]['markers']:
                    if not self._label_adjacent(j, came, to):
                        continue
                    if isinstance(to, SplitTreeNode):
                        stack.append((to.index, SplitTreeNode(j)))
                    else:
                        edges.append((x, to))
        return Graph([self._labels, edges], format="vertices_and_edges",
                     immutable=immutable)

    def __repr__(self):
        r"""
        Return a string representation of the split tree.

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.split_decomposition import split_tree
            sage: G = graphs.DiamondGraph()
            sage: ST = split_tree(G, check=True)
            sage: repr(ST)
            'Split tree on 4 vertices with 2 internal nodes'
        """
        return f"Split tree on {len(self._labels)} vertices with {len(self._info)} internal nodes"

    def __str__(self):
        r"""
        Return a string representation of the split tree.

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.split_decomposition import split_tree
            sage: G = graphs.DiamondGraph()
            sage: ST = split_tree(G, check=True)
            sage: print(ST)
            root leaf: 0
            N0  STAR  (parent 0)
                markers: 0 N1* 3
            N1  CLIQUE  (parent N0)
                markers: N0 1 2
        """
        lines = [f"root leaf: {self._root}"]
        for j, info in enumerate(self._info):
            star = lambda m: f"{m}*" if info['kind'] == 'STAR' and m == info['centre'] else f"{m}"
            lines.append(f"N{j}  {info['kind']}  (parent {info['parent']})")
            lines.append(f"    markers: {' '.join(star(m) for m in info['markers'])}")
            if info['kind'] == 'PRIME':
                lines.append(f"    label edges: {' '.join(f'{e[0]}-{e[1]}' for e in info['edges'])}")
        return "\n".join(lines)


def split_tree(G, check=True, lbfs_order=None, initial_vertex=None):
    r"""
    Return the split-tree (Cunningham) decomposition of ``G``.

    This method returns the split-tree of `G` as an instance of
    :class:`~sage.graphs.graph_decompositions.split_decomposition.SplitTree`
    which encodes the result of the construction proposed in [GPTC2014]_.
    See the documentation of module
    :mod:`~sage.graphs.graph_decompositions.split_decomposition`
    for more details.

    INPUT:

    - ``G`` -- a connected graph without loops

    - ``check`` -- boolean (default: ``True``); whether to run internal tests
      and procedures for checking the consistency of the internal data
      structures and the the validity of the construction during its execution.

    - ``lbfs_order`` -- list of vertices of ``G`` (default: ``None``); if
      provided, it is used as the LexBFS order for the incremental construction
      of the split-tree. Otherwise, a LexBFS order is computed.
      The validity of the order is not checked, it is assumed to be a valid
      LexBFS order.

    - ``initial_vertex`` -- a vertex of ``G`` (default: ``None``); if
      provided, it is used as the initial vertex for the incremental construction
      of the split-tree. Otherwise, an initial vertex is chosen.

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.split_decomposition import split_tree
        sage: G = graphs.BullGraph()
        sage: G.relabel(perm='abcde')
        sage: G.vertices()
        ['a', 'b', 'c', 'd', 'e']
        sage: ST = split_tree(G, check=True); ST
        Split tree on 5 vertices with 3 internal nodes
        sage: ST.root()
        'a'
        sage: ST.nodes()
        [N0, N1, N2]
        sage: sorted(ST.kind(u) for u in ST.nodes())
        ['CLIQUE', 'STAR', 'STAR']
        sage: ST.leaves()
        ['a', 'b', 'c', 'd', 'e']
        sage: A = ST.accessibility_graph()
        sage: A.is_isomorphic(G)
        True
        sage: T = ST.tree()
        sage: T.order()
        8
        sage: T.is_tree()
        True
        sage: sorted(T, key=str)
        [N0, N1, N2, 'a', 'b', 'c', 'd', 'e']

        sage: order = G.lex_BFS()
        sage: ST = split_tree(G, check=True, lbfs_order=order)
        sage: A = ST.accessibility_graph()
        sage: A.is_isomorphic(G)
        True
        sage: print(ST)
        root leaf: a
        N0  CLIQUE  (parent a)
            markers: a N1 N2
        N1  STAR  (parent N0)
            markers: N0 b* d
        N2  STAR  (parent N0)
            markers: N0 c* e
        sage: order = G.lex_BFS(initial_vertex='c')
        sage: ST = split_tree(G, check=True, lbfs_order=order)
        sage: ST.root()
        'c'
        sage: A = ST.accessibility_graph()
        sage: A.is_isomorphic(G)
        True
        sage: print(ST)
        root leaf: c
        N0  STAR  (parent c)
            markers: c* N1 e
        N1  CLIQUE  (parent N0)
            markers: N0 a N2
        N2  STAR  (parent N1)
            markers: N1 b* d

    TESTS::

        sage: split_tree(Graph()).tree()
        Graph on 0 vertices
        sage: split_tree(Graph(1)).tree()
        Graph on 1 vertex
        sage: split_tree(Graph(2)).tree()
        Traceback (most recent call last):
        ...
        ValueError: the graph must be connected
        sage: split_tree(DiGraph()).tree()
        Traceback (most recent call last):
        ...
        ValueError: the graph must be undirected

    Graphs with multiple edges are accepted, graphs with loops are not::

        sage: G = graphs.BullGraph()
        sage: A = split_tree(G, check=True).accessibility_graph()
        sage: G.allow_multiple_edges(True)
        sage: G.add_edges(G.edges())
        sage: B = split_tree(G, check=True).accessibility_graph()
        sage: A.is_isomorphic(B)
        True
        sage: B.allows_multiple_edges()
        False
        sage: split_tree(Graph([(0, 0)], loops=True))
        Traceback (most recent call last):
        ...
        ValueError: the graph must be without loops
    """
    if G.is_directed():
        raise ValueError("the graph must be undirected")
    if G.has_loops():
        raise ValueError("the graph must be without loops")
    if not G.is_connected():
        raise ValueError("the graph must be connected")

    cdef Py_ssize_t n = G.order()
    cdef GltExport out

    if n < 3:
        return SplitTree(list(G), {})

    cdef Py_ssize_t i
    cdef list labels = list(G)
    cdef dict index = {v: i for i, v in enumerate(labels)}
    cdef vector[vector[int]] adjacency
    adjacency.resize(n)
    for i, u in enumerate(labels):
        for w in G.neighbor_iterator(u):
            if w != u:
                adjacency[i].push_back(<int> index[w])    

    cdef int start = 0 if initial_vertex is None else <int> index[initial_vertex]
    cdef vector[int] given
    if lbfs_order is not None:
        if len(lbfs_order) != n or len(set(lbfs_order)) != n:
            raise ValueError("lbfs_order must list every vertex exactly once")
        for v in lbfs_order:
            given.push_back(<int> index[v])

    cdef string error
    sig_on()
    cdef bint success = glt_build(adjacency, given, start, check, out, error)
    sig_off()
    if not success:
        raise RuntimeError(error.decode())

    exported = {'root_vertex': out.root_vertex,
                'node_kind': out.node_kind,
                'node_parent': out.node_parent,
                'node_star_centre': out.node_star_centre,
                'node_marker_start': out.node_marker_start,
                'node_marker_neighbour': out.node_marker_neighbour,
                'node_edge_start': out.node_edge_start,
                'node_edge_pairs': out.node_edge_pairs}
    return SplitTree(labels, exported)


def GP2012():
    r"""
    Return the graph used in [GP2012b]_.

    This is useful for testing the correctness of the algorithm.

    EXAMPLE::

        sage: from sage.graphs.graph_decompositions.split_decomposition import split_tree
        sage: from sage.graphs.graph_decompositions.split_decomposition import GP2012
        sage: G = GP2012()
        sage: ST = split_tree(G, check=True); ST
        Split tree on 16 vertices with 6 internal nodes
        sage: T = ST.tree()
        sage: T.is_tree()
        True
        sage: A = ST.accessibility_graph()
        sage: A.is_isomorphic(G)
        True
    """
    from sage.graphs.graph import Graph

    def add_join(G, A, B):
        G.add_edges((a, b) for a in A for b in B)

    G = Graph([range(1, 17), []], format="vertices_and_edges")
    G.add_cycle([1, 2, 3, 4])
    add_join(G, [1, 4], [5, 6, 7, 10, 11, 12, 13, 16])
    G.add_edge(5, 6)
    add_join(G, [5, 6], [7, 10])
    G.add_edge(7, 8)
    G.add_cycle([8, 9, 10])
    G.add_cycle([12, 13, 14])
    G.add_path([14, 15, 16])
    return G
