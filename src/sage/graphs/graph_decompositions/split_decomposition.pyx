r"""
Split-tree (Cunningham) decomposition of connected graphs.

This module implements the incremental split decomposition algorithm described
in [GPTC2016]_. The output is the *split-tree* of the input graph, i.e., the
unique reduced graph-labelled tree (GLT) whose accessibility graph is the input
graph (Cunningham's canonical split decomposition, reformulated in terms of
GLTs).

The algorithm inserts the vertices one at a time, following a LexBFS ordering.
For each inserted vertex `x`, the current split-tree is marked with respect to
the set `S` of already inserted neighbours of `x`; exactly one of the seven
cases of Theorem 4.1 of [GPTC2016]_ applies, and the tree is updated by the
local construction of the corresponding proposition.

This is a *simple* implementation: states are recomputed from scratch after
every modification. This is not the quasi-linear time construction of the data
structure described in the paper.

EXAMPLES::

    sage: from sage.graphs.graph_decompositions.split_decomposition import split_tree
    sage: from sage.graphs.graph_decompositions.split_decomposition import GP2012
    sage: G = GP2012()
    sage: ST = split_tree(G)
    sage: ST.is_reduced()
    True
    sage: ST  # random
    u32: CLIQUE [q29-q30, q29-q31, q29-q33, q30-q31, q30-q33, q31-q33]
    u43: CLIQUE [q40-q41, q40-q42, q41-q42]
    u48: PRIME [q17-q20, q17-q45, q20-q45, q20-q49, q45-q46, q46-q49]
    u52: STAR (centre q25) [q25-q27, q25-q38, q25-q50]
    u54: PRIME [q22-q23, q22-q51, q23-q55, q39-q51, q39-q55]
    u60: PRIME [q35-q56, q35-q58, q36-q56, q36-q58, q36-q61, q58-q61]
    16 -- q22 -- u54
    1 -- q45 -- u48
    4 -- q20 -- u48
    15 -- q23 -- u54
    5 -- q30 -- u32
    6 -- q31 -- u32
    7 -- q35 -- u60
    10 -- q36 -- u60
    11 -- q38 -- u52
    12 -- q41 -- u43
    13 -- q42 -- u43
    2 -- q46 -- u48
    3 -- q49 -- u48
    14 -- q55 -- u54
    8 -- q58 -- u60
    9 -- q61 -- u60
    sage: T = ST.as_graph()
    sage: T.is_tree()
    True
    sage: sorted(T, key=str)
    [CLIQUE: Graph on 4 vertices,
     CLIQUE: Graph on 3 vertices,
     LEAF(1): Graph on 1 vertex,
     LEAF(10): Graph on 1 vertex,
     LEAF(11): Graph on 1 vertex,
     LEAF(12): Graph on 1 vertex,
     LEAF(13): Graph on 1 vertex,
     LEAF(14): Graph on 1 vertex,
     LEAF(15): Graph on 1 vertex,
     LEAF(16): Graph on 1 vertex,
     LEAF(2): Graph on 1 vertex,
     LEAF(3): Graph on 1 vertex,
     LEAF(4): Graph on 1 vertex,
     LEAF(5): Graph on 1 vertex,
     LEAF(6): Graph on 1 vertex,
     LEAF(7): Graph on 1 vertex,
     LEAF(8): Graph on 1 vertex,
     LEAF(9): Graph on 1 vertex,
     PRIME: Graph on 5 vertices,
     PRIME: Graph on 5 vertices,
     PRIME: Graph on 5 vertices,
     STAR: Graph on 4 vertices]
    sage: A = ST.accessibility_graph()
    sage: A.is_isomorphic(G)
    True
    sage: set(A) == set(G)
    True
    sage: EA = set(frozenset(e) for e in A.edges(labels=False))
    sage: EA == set(frozenset(e) for e in G.edges(labels=False))
    True

.. TODO::

    Implement a truly quasi-linear time version of the construction of the
    graph-labelled tree data structure as described in [GPTC2016]_.
"""
# ******************************************************************************
#       Copyright (C) 2026 David Coudert <david.coudert@inria.fr>
#                     2026 Martin Rubey <martin.rubey at tuwien.ac.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ******************************************************************************

from collections import deque
from enum import IntEnum

from cysignals.signals cimport sig_on, sig_off

from sage.graphs.generic_graph_pyx cimport GenericGraph_pyx

# ------------------------------------------------------------------------------
# Enumeration class and types
# ------------------------------------------------------------------------------

class NodeType(IntEnum):
    r"""
    Enumeration class used to define the various types of nodes.

    The various node types defined are

    - ``CLIQUE`` -- indicates the node is a clique module

    - ``STAR`` -- indicates the node is a star module

    - ``PRIME`` -- indicates the node is a prime module
    """
    CLIQUE = 0
    STAR = 1
    PRIME = 2

    def __repr__(self) -> str:
        r"""
        Return a string representation of a ``NodeType`` object.

        TESTS::

            sage: from sage.graphs.graph_decompositions.split_decomposition import NodeType
            sage: repr(NodeType.CLIQUE)
            'CLIQUE'
            sage: str(NodeType.STAR)
            'STAR'
            sage: S = NodeType.STAR
            sage: S in NodeType
            True
            sage: 37 in NodeType
            False
        """
        return self.name

    __str__ = __repr__


# Enumeration type used to specify the state of markers / leaves
# _E_: empty; _M_: mixed; _P_: perfect
ctypedef enum State: _E_, _M_, _P_

# Enumeration type used in methods subdivise and insert to specify the role of
# a marker with respect a star
ctypedef enum Role: _THIS_, _OTHER_


# ------------------------------------------------------------------------------
# graph-labelled tree
# ------------------------------------------------------------------------------

cdef class GraphLabelledTree:
    r"""
    A graph-labelled tree (GLT).

    This class implements the incremental split decomposition algorithm
    described in [GPTC2016]_. The output is the *split-tree* of the input graph,
    i.e., the unique reduced graph-labelled tree (GLT) whose accessibility graph
    is the input graph (Cunningham's canonical split decomposition, reformulated
    in terms of GLTs).

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.split_decomposition import GraphLabelledTree
        sage: G = graphs.DiamondGraph()
        sage: ST = GraphLabelledTree(G, check=True)
        sage: sorted(ST.as_graph(), key=str)
        [CLIQUE: Graph on 3 vertices,
         LEAF(0): Graph on 1 vertex,
         LEAF(1): Graph on 1 vertex,
         LEAF(2): Graph on 1 vertex,
         LEAF(3): Graph on 1 vertex,
         STAR: Graph on 3 vertices]
        sage: H = ST.accessibility_graph()
        sage: H.is_isomorphic(G)
        True
    """
    cdef readonly GenericGraph_pyx _G
    # for internal work, use a graph relabelled in 0..n-1
    cdef readonly GenericGraph_pyx _G_int
    cdef list int_to_vertex  # vertex of G_int -> vertex of G

    # label maps each node of the GLT to its adjacency dict
    cdef dict label  # node -> {marker: set(markers)}
    # node_of maps each marker to the node of the GLT containing it
    cdef dict node_of  # marker -> node
    # opposite maps each marker to the marker at the other end of its tree-edge,
    # or the *leaf* (a vertex of _G_int) if the tree-edge is incident to a leaf
    cdef dict opposite  # marker -> marker | leaf
    # leaf_marker maps each leaf to the marker it is attached to
    cdef dict leaf_marker  # leaf -> marker
    # markers and nodes are labelled with integers. To avoid collisions with
    # vertices of _G_int, we assign distinct integers larger than n
    cdef size_t _cnt  # used to assign distinct labels to nodes and markers

    cdef bint _check_mode  # whether to apply internal controls

    def __init__(self, G, check=True):
        r"""
        Initialize data structures.

        This method build the internal representation of the graph-labelled
        tree.

        INPUT:

        - ``G`` -- a graph

        - ``check`` -- boolean (default: ``True``); whether to run internal
          tests and procedures for checking the consistency of the internal data
          structures and the the validity of the construction during its
          execution. This may significantly increase the running time.

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.split_decomposition import GraphLabelledTree
            sage: G = graphs.DiamondGraph()
            sage: ST = GraphLabelledTree(G, check=True)
            sage: sorted(ST.as_graph(), key=str)
            [CLIQUE: Graph on 3 vertices,
             LEAF(0): Graph on 1 vertex,
             LEAF(1): Graph on 1 vertex,
             LEAF(2): Graph on 1 vertex,
             LEAF(3): Graph on 1 vertex,
             STAR: Graph on 3 vertices]
            sage: H = ST.accessibility_graph()
            sage: H.is_isomorphic(G)
            True

        The construction is correct with and without internal validation tests::

            sage: from sage.graphs.graph_decompositions.split_decomposition import GP2012
            sage: G = GP2012()
            sage: ST = GraphLabelledTree(G, check=True)
            sage: ST.accessibility_graph().is_isomorphic(G)
            True
            sage: ST = GraphLabelledTree(G, check=False)
            sage: ST.accessibility_graph().is_isomorphic(G)
            True

        TESTS::

            sage: GraphLabelledTree(Graph()).as_graph()
            Graph on 0 vertices
            sage: GraphLabelledTree(Graph(1)).as_graph()
            Graph on 2 vertices
            sage: GraphLabelledTree(Graph(2))
            Traceback (most recent call last):
            ...
            ValueError: the graph must be connected
            sage: GraphLabelledTree(DiGraph())
            Traceback (most recent call last):
            ...
            ValueError: the input graph must be undirected
            sage: GraphLabelledTree(Graph([(0, 0)], loops=True))
            Traceback (most recent call last):
            ...
            ValueError: the input graph must be without loops
            sage: GraphLabelledTree(Graph([(0, 1), (0, 1)], multiedges=True)).as_graph()
            Graph on 3 vertices
        """
        if G.is_directed():
            raise ValueError("the input graph must be undirected")
        if G.has_loops():
            raise ValueError("the input graph must be without loops")
        if not G.is_connected():
            raise ValueError("the graph must be connected")

        # Store the graph and a copy relabeled with integers in 0..n-1.
        # Internal computations are performed on the relabeled graph.
        self._G = <GenericGraph_pyx?> G
        self.int_to_vertex = list(G)
        cdef dict vertex_to_int = {u: i for i, u in enumerate(self.int_to_vertex)}
        self._G_int = G.relabel(perm=vertex_to_int, inplace=False)

        # Initialize data structures
        self.label = {}
        self.node_of = {}
        self.opposite = {}
        self.leaf_marker = {}
        # markers are indexed with integers >= n
        self._cnt = G.order()
        self._check_mode = check

        # Compute the graph-labelled tree
        if G:
            sig_on()
            self.make_split_tree()
            sig_off()

    cdef size_t _fresh(self) noexcept:
        r"""
        Return a new marker or node label.
        """
        self._cnt += 1
        return self._cnt

    # --------------------------------------------------------------------------
    # split-tree computation
    # --------------------------------------------------------------------------

    cdef make_split_tree(self):
        r"""
        Build the internal representation of the graph-labelled tree.

        When ``self._check_mode`` is ``True``, the method verifies after every
        insertion that the accessibility graph is the induced subgraph on the
        inserted vertices.

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.split_decomposition import GraphLabelledTree
            sage: ST = GraphLabelledTree(graphs.DiamondGraph(), check=True)
        """
        cdef size_t n = self._G.order()
        cdef size_t x, x1, x2, q1, q2, u
        if n == 1:
            x1 = 0
            q1 = self._fresh()
            self.add_node({q1: set()})
            self.opposite[q1] = x1
            self.leaf_marker[x1] = q1
            return
    
        cdef list order = self._G_int.lex_BFS()
    
        # --- initial split-tree on the first two vertices ---------------
        x1 = order[0]
        x2 = order[1]
        q1 = self._fresh()
        q2 = self._fresh()
        self.add_node({q1: {q2}, q2: {q1}})
        self.opposite[q1] = x1
        self.opposite[q2] = x2
        self.leaf_marker[x1] = q1
        self.leaf_marker[x2] = q2
        cdef set inserted = set([x1, x2])
        if self._check_mode and n == 2:
            self._check(self._G_int.subgraph(inserted))

        cdef dict neighbors = {u: set(self._G_int.neighbors(u)) for u in order}

        # --- incremental insertion --------------------------------------
        cdef set S
        cdef dict nbrs = {v: neighbors[v] & inserted for v in inserted}
        for x in order[2:]:
            S = neighbors[x] & inserted
            if self._check_mode:
                assert S, "LBFS ordering does not have connected prefixes"
            self.insert(x, S, nbrs)
            # update inserted and nbrs
            inserted.add(x)
            nbrs[x] = S
            for x2 in S:
                nbrs[x2].add(x)
            if self._check_mode:
                self._check(self._G_int.subgraph(inserted))

    # --------------------------------------------------------------------------
    # Method for checking the validity of the construction
    # --------------------------------------------------------------------------

    cdef _check(self, G):
        r"""
        Check the validity of the construction for ``G``.

        The construction is incremental. At each step of the construction, we
        can check that the internal data structures are consistent and that the
        current accessibility graph is the induced subgraph on the inserted
        vertices.

        INPUT:

        - ``G`` -- the graph of the already inserted vertices, which is a
          subgraph of ``self._G_int``

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.split_decomposition import GraphLabelledTree
            sage: ST = GraphLabelledTree(graphs.DiamondGraph(), check=True)
        """
        # Check the internal consistency of the data structures
        for u, adj in self.label.items():
            V = set(adj)
            for m, s in adj.items():
                assert s <= V, ("stale neighbour in %s" % u, m)
                for m2 in s:
                    assert m in adj[m2], ("adjacency not symmetric", m, m2)
                assert self.node_of.get(m) == u, ("wrong node_of", m)
        assert set(self.node_of) == {m for adj in self.label.values() for m in adj}
        for q, t in self.opposite.items():
            if t in self.node_of:
                assert self.opposite[t] == q, ("opposite not symmetric", q)
            else:
                assert self.leaf_marker[t] == q
        for leaf, q in self.leaf_marker.items():
            assert self.opposite[q] == leaf

        # Check that internal data structures properly encode the graph G
        H = self.accessibility_graph(algorithm='A', internal=True)
        if (set(H) != set(G)
                or set(map(frozenset, H.edges(labels=False)))
                != set(map(frozenset, G.edges(labels=False)))):
            raise AssertionError("accessibility graph mismatch on %s" % (set(G),))

    # --------------------------------------------------------------------------
    # structural operations
    # --------------------------------------------------------------------------

    cdef tuple classify_label(self, dict adj):
        r"""
        Return the kind of label is represented by ``adj``.

        Classify the label graph ``adj`` (a dict ``marker -> set of markers``)
        as either a clique, a star or a prime graph. If the label is a star, the
        method also returns its centre.

        INPUT:

        - ``adj`` -- dictionary; adjacency dictionary of a label

        OUTPUT:

        Return a tuple ``(NodeType.CLIQUE, None)``, ``(NodeType.STAR, centre)``
        or ``(NodeType.PRIME, None)``.

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.split_decomposition import GraphLabelledTree
            sage: ST = GraphLabelledTree(graphs.CompleteGraph(5), check=True)
            sage: sorted(ST.as_graph(), key=str)
            [CLIQUE: Graph on 5 vertices,...]
            sage: ST = GraphLabelledTree(graphs.StarGraph(5), check=True)
            sage: sorted(ST.as_graph(), key=str)
            [..., STAR: Graph on 6 vertices]
            sage: ST = GraphLabelledTree(graphs.PetersenGraph(), check=True)
            sage: sorted(ST.as_graph(), key=str)
            [..., PRIME: Graph on 10 vertices]

        TESTS::

            sage: ST = GraphLabelledTree(Graph(1), check=True)
            sage: sorted(ST.as_graph(), key=str)
            [CLIQUE: Graph on 1 vertex, LEAF(0): Graph on 1 vertex]
            sage: ST = GraphLabelledTree(Graph([('a', 'b')]), check=True)
            sage: sorted(ST.as_graph(), key=str)
            [CLIQUE: Graph on 2 vertices,
             LEAF(a): Graph on 1 vertex,
             LEAF(b): Graph on 1 vertex]
        """
        cdef int n = len(adj)
        if n <= 2:
            # Degenerated case (K_2 or a single vertex), treat as clique
            return (NodeType.CLIQUE, None)
        cdef int degree_one = 0  # number of vertices of degree 1
        cdef int degree_n_minus_one = 0  # number of vertices of degree n - 1
        cdef int d
        centre = None  # potential center of a star
        for u, neighbors in adj.items():
            d = len(neighbors)
            if d == 1:
                degree_one += 1
            elif d == n - 1:
                degree_n_minus_one += 1
                centre = u
        if degree_n_minus_one == n:
            return (NodeType.CLIQUE, None)
        if centre is not None and degree_one == n - 1:
            return (NodeType.STAR, centre)
        return (NodeType.PRIME, None)

    cdef int add_node(self, dict adj) noexcept:
        r"""
        Add a new node with the given label adjacency and return it.

        INPUT:

        - ``adj`` -- dictionary; adjacency dictionary of a label

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.split_decomposition import GraphLabelledTree
            sage: ST = GraphLabelledTree(graphs.DiamondGraph(), check=True)
        """
        cdef int u = self._fresh()
        self.label[u] = adj
        cdef int m
        for m in adj:
            self.node_of[m] = u
        return u

    cdef int attach_leaf(self, int u, int leaf, set neighbours) noexcept:
        r"""
        Add a marker to node ``u``, label-adjacent to ``neighbours``
        (markers of ``u``), with tree-edge to ``leaf``.

        INPUT:

        - ``u`` -- integer; a node

        - ``leaf`` -- integer; a leaf

        - ``neighbors`` -- the set of neighbors of ``leaf`` in node ``u``

        OUTPUT: the new marker ``q`` associated to ``leaf``

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.split_decomposition import GraphLabelledTree
            sage: ST = GraphLabelledTree(graphs.DiamondGraph(), check=True)
            sage: sorted(ST.as_graph().edges(), key=str)
            [(LEAF(0): Graph on 1 vertex, STAR: Graph on 3 vertices, (0, 'q10', 'LEAF')),
             (LEAF(1): Graph on 1 vertex, CLIQUE: Graph on 3 vertices, (1, 'q6', 'LEAF')),
             (LEAF(2): Graph on 1 vertex, CLIQUE: Graph on 3 vertices, (2, 'q8', 'LEAF')),
             (LEAF(3): Graph on 1 vertex, STAR: Graph on 3 vertices, (3, 'q11', 'LEAF')),
             (STAR: Graph on 3 vertices, CLIQUE: Graph on 3 vertices, ('q9', 'q5', None))]
        """
        if self._check_mode:
            assert leaf not in self.leaf_marker
        cdef int q = self._fresh()  # marker associated to leaf
        cdef dict label_u = self.label[u]
        label_u[q] = set(neighbours)
        cdef int a
        for a in neighbours:
            label_u[a].add(q)
        self.node_of[q] = u
        self.opposite[q] = leaf
        self.leaf_marker[leaf] = q
        return q

    cdef tuple node_split(self, int u, set A, set B):
        r"""
        Split node ``u`` according the bipartition ``(A, B)`` of its markers.

        Assuming that the bipartition ``(A, B)`` is a split of node ``u``, this
        method creates two new labels ``uA`` and ``uB``, each containing the
        markers of ``A`` (resp., ``B``) plus a new marker ``qA`` (resp., ``qB``)
        encoding the split on the side of ``A`` (resp., ``B``).  ``qA`` and
        ``qB`` are the extremities of the new tree-edge and are universal to the
        frontiers ``N(B) ∩ A`` and ``N(A) ∩ B`` respectively.  Node ``u`` and
        its label are then removed from the data-structures.

        INPUT:

        - ``u`` -- integer; a node

        - ``A, B`` -- disjoint sets of size at least 2 representing a
          bipartition of the markers of the label of node ``u``.

        OUTPUT:
        
        A tuple ``(uA, uB, qA, qB)`` where ``uA`` is the new marker associated
        with part ``A``, ``uB`` is the new marker associated with part ``B``,
        ``qA`` and ``qB`` are the new markers encoding the split.

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.split_decomposition import GraphLabelledTree
            sage: ST = GraphLabelledTree(graphs.DiamondGraph(), check=True)
        """
        if self._check_mode:
            assert A | B == set(self.label[u]) and not (A & B)
            assert len(A) >= 2 and len(B) >= 2, "not a split of the label"
        cdef dict adj = self.label[u]
        cdef int m, v
        # Compute the frontiers on the sides of A and B
        cdef set Ap = {m for m in A if any(v in B for v in adj[m])}
        cdef set Bp = {m for m in B if any(v in A for v in adj[m])}
        # Create new markers to represent the split and the tree-edge
        cdef int qA = self._fresh()
        cdef int qB = self._fresh()
        # Compute the labels for A and B
        cdef dict adjA = {m: adj[m] & A for m in A}
        adjA[qA] = Ap
        for m in Ap:
            adjA[m].add(qA)
        cdef dict adjB = {m: adj[m] & B for m in B}
        adjB[qB] = Bp
        for m in Bp:
            adjB[m].add(qB)
        # Store the new nodes
        cdef int uA = self.add_node(adjA)
        cdef int uB = self.add_node(adjB)
        del self.label[u]
        # Record the edge
        self.opposite[qA] = qB
        self.opposite[qB] = qA
        return uA, uB, qA, qB

    cdef int node_join(self, int q) noexcept:
        r"""
        Join the labels along the internal tree-edge ``q, opposite[q]``.

        This is the reverse operation of ``node_split``.
        The two labels ``q`` and ``opposite[q]`` are merged; all label-edges
        between ``N(q)`` and ``N(opposite(q))`` are added; the data associated
        with ``q`` and ``opposite(q)`` are deleted.

        INPUT:

        - ``q`` -- integer; a side of the internal tree-edge to merge

        OUTPUT: the new node

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.split_decomposition import GraphLabelledTree
            sage: ST = GraphLabelledTree(graphs.DiamondGraph(), check=True)
        """
        cdef int r = self.opposite.pop(q)  # extract and remove q and opposite[q]
        if self._check_mode:
            assert r in self.node_of, "not an internal tree-edge"
        cdef int u = self.node_of.pop(q)  # extract and remove q and node_of[q]
        cdef int v = self.node_of.pop(r)
        # Extract the labels of u and v and merge them into the label of u
        cdef dict adjU = self.label.pop(u)  # extract and remove u and label[u]
        cdef dict adjV = self.label.pop(v)
        cdef set Nq = adjU.pop(q)  # extract and remove q and adjU[q]
        cdef set Nr = adjV.pop(r)
        cdef set s
        # remove traces of q from adjU and r from adjV
        for s in adjU.values():
            s.discard(q)
        for s in adjV.values():
            s.discard(r)
        # merge the labels into adjU and add the join
        adjU.update(adjV)
        cdef int a, b
        for a in Nq:
            for b in Nr:
                adjU[a].add(b)
                adjU[b].add(a)
        del self.opposite[r]
        return self.add_node(adjU)

    cdef int subdivide(self, int q, int leaf, tuple spec) noexcept:
        r"""
        Subdivide the tree-edge at marker ``q`` with a new ternary node carrying
        the new ``leaf`` (used in cases 4, 5, 6 of [GPTC2016]_).

        INPUT:

        - ``u`` -- integer; a node

        - ``leaf`` -- integer; a leaf

        - ``spec`` -- tuple; either ``(CLIQUE,)`` or ``(STAR, role)`` where
          ``role`` is ``_THIS_`` if the centre of the star is the new marker
          opposite ``q``, and ``_OTHER_`` if it is the new marker opposite the
          other extremity of the subdivided edge.

        OUTPUT: the new marker ``w``

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.split_decomposition import GraphLabelledTree
            sage: ST = GraphLabelledTree(graphs.DiamondGraph(), check=True)
        """
        cdef int t = self.opposite[q]
        if self._check_mode:
            assert t in self.node_of or t in self.leaf_marker
        cdef int a = self._fresh()
        cdef int b = self._fresh()
        cdef int c = self._fresh()
        cdef dict adj
        if spec[0] == NodeType.CLIQUE:
            adj = {a: {b, c}, b: {a, c}, c: {a, b}}
        elif spec[1] == _THIS_:
            adj = {a: {b, c}, b: {a}, c: {a}}
        else:
            adj = {b: {a, c}, a: {b}, c: {b}}
        cdef int w = self.add_node(adj)
        # rewire the tree edges
        self.opposite[q] = a
        self.opposite[a] = q
        self.opposite[b] = t
        if t in self.node_of:
            self.opposite[t] = b
        else:
            self.leaf_marker[t] = b
        self.opposite[c] = leaf
        self.leaf_marker[leaf] = c
        return w

    # --------------------------------------------------------------------------
    # states
    # --------------------------------------------------------------------------

    cdef dict _node_tree(self):
        r"""
        Return the adjacency of the underlying tree on (nodes | leaves).

        This is a helper method for method ``L_sets``.
        There is one entry in ``adj`` for each node and each leaf. 

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.split_decomposition import GraphLabelledTree
            sage: ST = GraphLabelledTree(graphs.DiamondGraph(), check=True)
        """
        cdef dict adj = {u: set() for u in self.label}
        cdef int u, v, m, t
        for m, t in self.opposite.items():
            u = self.node_of[m]
            v = self.node_of[t] if t in self.node_of else t
            adj.setdefault(v, set())
            adj[u].add(v)
            adj[v].add(u)
        return adj

    cdef dict L_sets(self):
        r"""
        Return reachable leaf sets for each ``q`` in ``opposite``.

        OUTPUT:
        
        Dictionary ``{marker: set of leaves beyond its tree-edge}`` with one
        entry per marker ``q`` in ``opposite``.

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.split_decomposition import GraphLabelledTree
            sage: ST = GraphLabelledTree(graphs.DiamondGraph(), check=True)
        """
        cdef dict tree_adj = self._node_tree()
        cdef dict memo = {}  # used to memoize calls to beyond

        def beyond(u, v):
            """
            Return the set of leaves in the sub-tree rooted at ``v`` after
            the removal of tree-edge ``(u, v)``.

            Perform a BFS from ``v`` avoiding ``u`` (crossing ``u`` is
            forbidden), and return the reached leaves.
            
            Thanks to memoization, the time complexity for the sum of the calls
            to beyond is linear in the number of tree-edges.
            """
            key = (u, v)
            if key not in memo:
                leaves = {v} if v in self.leaf_marker else set()
                for x in tree_adj[v]:
                    if x == u:
                        continue
                    leaves.update(beyond(v, x))
                memo[key] = leaves
            return memo[key]

        cdef dict res = {}
        cdef int u, v, q, t
        for q, t in self.opposite.items():
            u = self.node_of[q]
            v = self.node_of[t] if t in self.node_of else t
            res[q] = beyond(u, v)
        return res

    cdef set A_set(self, int q):
        r"""
        Return the set of leaves of ``L(q)`` accessible from the marker ``q``.

        The set ``A(q)`` of leaves is obtained using an alternating sequence
        from ``q`` to a leaf of ``L(q)`` that starts by crossing ``q``'s
        tree-edge. So it explores the alternating (label / tree) walk from
        ``opposite(q)`` onwards.

        INPUT:

        - ``q`` -- integer; a side of the internal tree-edge to merge

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.split_decomposition import GraphLabelledTree
            sage: ST = GraphLabelledTree(graphs.DiamondGraph(), check=True)
        """
        cdef int r = self.opposite[q]
        if r not in self.node_of:  # r is a leaf
            return {r}
        cdef set seen_zero = {r}
        cdef set seen_one = set()
        dq = deque([(r, 0)])
        cdef set leaves = set()
        cdef int k, m, m2, t
        while dq:
            m, k = dq.popleft()
            if k == 0:  # next move: label-edge
                for m2 in self.label[self.node_of[m]][m]:
                    if m2 not in seen_one:
                        seen_one.add(m2)
                        dq.append((m2, 1))
            else:  # next move: tree-edge
                t = self.opposite[m]
                if t in self.node_of:
                    if t not in seen_zero:
                        seen_zero.add(t)
                        dq.append((t, 0))
                else:
                    leaves.add(t)
        return leaves

    cdef tuple states(self, set S, dict nbrs):
        r"""
        Return the states of all markers and leaves with respect to ``S``.

        INPUT:

        - ``S`` -- set; a set of leaves

        - ``nbrs`` -- dictionary; map each leaf to its neighbours among the
          current leaves

        OUTPUT:

        A pair of dictionary ``(ms, ls)`` respectively mapping each marker and
        leaf to a ``State``: ``_E_``, ``_M_`` or ``_P_``.

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.split_decomposition import GraphLabelledTree
            sage: ST = GraphLabelledTree(graphs.DiamondGraph(), check=True)
        """
        cdef dict Ls = self.L_sets()
        cdef dict ms = {}
        cdef set inter
        cdef int q, leaf
        for q in self.opposite:
            inter = S & Ls[q]
            if not inter:
                ms[q] = _E_
            elif inter == self.A_set(q):
                ms[q] = _P_
            else:
                ms[q] = _M_
        cdef dict ls = {}
        for leaf in self.leaf_marker:
            inter = S - {leaf}
            if not inter:
                ls[leaf] = _E_
            elif inter == set(nbrs[leaf]):
                ls[leaf] = _P_
            else:
                ls[leaf] = _M_
        return ms, ls

    # --------------------------------------------------------------------------
    # vertex insertion (Theorem 4.1 and Props. 4.2 - 4.5 of [GPTC2016]_)
    # --------------------------------------------------------------------------

    cdef void insert(self, int x, set S, dict nbrs) noexcept:
        """
        Insert a new leaf ``x``.

        This is algorithm vertex-insertion as described in Theorem 4.1
        of [GPTC2016]_.

        INPUT:

        - ``x`` -- integer; a new leaf

        - ``S`` -- set; the neighbourhood among the current leaves of ``x``

        - ``nbrs`` -- dictionary; map each current leaf to its neighbours among
          the current leaves

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.split_decomposition import GraphLabelledTree
            sage: ST = GraphLabelledTree(graphs.DiamondGraph(), check=True)
        """
        cdef dict ms, ls
        ms, ls = self.states(S, nbrs)

        def state_of(t):
            return ms[t] if t in self.node_of else ls[t]

        def edges(pred):
            """
            Return one representative marker per tree-edge.
            
            The extremity state of each representative marker satisfy function
            ``pred``.
            """
            res = []
            for q, state_q in ms.items():
                t = self.opposite[q]
                if t in self.node_of and t < q:
                    continue  # internal edge, seen from the other end
                if pred(state_q, state_of(t)):
                    res.append(q)
            return res

        cdef list pp = edges(lambda a, b: a == _P_ and b == _P_)
        cdef dict cnt
        cdef set seen
        cdef int q, m, c
        cdef list us
        if pp:
            if len(pp) == 1:
                # case 5: unique PP tree-edge (Prop. 4.4, case 5): subdivide
                # with a ternary clique node
                self.subdivide(pp[0], x, (NodeType.CLIQUE,))
                return
            # case 1: unique clique node, all markers are perfect, (Prop. 4.2,
            # case 1): add a new marker incident to all PP tree-edges
            if self._check_mode:
                cnt = {}
                for q in pp:
                    for u in (self.node_of[q],
                              self.node_of.get(self.opposite[q])):
                        if u is not None:
                            cnt[u] = cnt.get(u, 0) + 1
                # search for the unique node u s.t. c >= 2
                us = [u for u, c in cnt.items() if c >= 2]
                assert len(us) == 1, "case 1: no clique node"
                u = us[0]
                assert all(self.node_of[q] == u
                           or self.node_of[self.opposite[q]] == u for q in pp)
                assert all(ms[m] == _P_ for m in self.label[u])
                self.attach_leaf(u, x, set(self.label[u]))
                return
            # search for the unique node u appearing twice
            seen = set()
            for q in pp:
                for u in (self.node_of[q], self.node_of.get(self.opposite[q])):
                    if u is not None:
                        if u in seen:
                            # we have the right node
                            self.attach_leaf(u, x, set(self.label[u]))
                            return
                        seen.add(u)
            raise ValueError("this case should never happen")

        cdef list pe = edges(lambda a, b: {a, b} == {_P_, _E_})
        cdef set P_u
        cdef Role role
        if pe:
            if len(pe) == 1:
                # case 6: unique PE tree-edge (Prop. 4.4, case 6): subdivide
                # with a ternary star node whose centre is opposite the empty
                # extremity
                q = pe[0]
                role = _THIS_ if ms[q] == _E_ else _OTHER_
                self.subdivide(q, x, (NodeType.STAR, role))
                return
            # case 2: unique star node, all markers empty except the perfect
            # centre (Prop. 4.2, case 2): add a new marker incident to all PE
            # tree-edges
            if self._check_mode:
                cnt = {}
                for q in pe:
                    for u in (self.node_of[q],
                              self.node_of.get(self.opposite[q])):
                        if u is not None:
                            cnt[u] = cnt.get(u, 0) + 1
                # search for the unique node u s.t. cnt[u] >= 2
                us = [u for u, c in cnt.items() if c >= 2]
                assert len(us) == 1, "case 2: no star node"
                u = us[0]
                assert all(self.node_of[q] == u
                           or self.node_of[self.opposite[q]] == u for q in pe)
                P_u = {m for m in self.label[u] if ms[m] == _P_}
                assert len(P_u) == 1
                kind, centre = self.classify_label(self.label[u])
                if kind == NodeType.STAR:
                    assert P_u == {centre}
                else:
                    # a K_2 label is at the same time a clique and a star
                    # K_{1,1}; the star reading required here (only possible
                    # when inserting the third vertex) has the perfect marker as
                    # its centre.
                    assert kind == NodeType.CLIQUE and len(self.label[u]) == 2
                    assert all(ms[m] == _E_ for m in self.label[u] if m not in P_u)
                self.attach_leaf(u, x, P_u)
                return
            # get the unique node u appearing twice
            seen = set()
            for q in pe:
                for u in (self.node_of[q], self.node_of.get(self.opposite[q])):
                    if u is not None:
                        if u in seen:
                            # we have the right node; we attach it
                            for m in self.label[u]:
                                if ms[m] == _P_:
                                    self.attach_leaf(u, x, {m})
                                    return
                        seen.add(u)
            raise ValueError("this case should never happen")

        cdef list mm = edges(lambda a, b: a == _M_ and b == _M_)
        if mm:
            # case 7: fully-mixed subtree (Prop. 4.5, case 7): clean the
            # degenerate nodes of the fully-mixed subtree, contract the
            # fully-mixed subtree to a single node by node-joins, and add a
            # marker adjacent to its perfect markers.
            self._clean_and_contract(x, S, nbrs, ms)
            return

        # unique hybrid node (all markers perfect or empty, all opposites mixed)
        cdef set Pstar, Estar
        cdef int uA, uB, qA, qB
        cdef list hybrid = [
            u for u, V in self.label.items()
            if (all(ms[m] in (_P_, _E_) for m in V) and
                all(state_of(self.opposite[m]) == _M_ for m in V))]
        assert len(hybrid) == 1, "expected a unique hybrid node"
        u = hybrid[0]
        P_u = {m for m in self.label[u] if ms[m] == _P_}
        if self._check_mode:
            assert P_u
        classified_label = self.classify_label(self.label[u])
        if classified_label[0] == NodeType.PRIME:
            # case 3: unique hybrid prime node (Prop. 4.2, case 3): add a marker
            # adjacent to the perfect markers
            self.attach_leaf(u, x, P_u)
            return
        # case 4: hybrid degenerate node, all markers are perfect or empty, so
        # the centre, if any, belongs to E* (Prop. 4.3, case 4): node-split on
        # (P*, E*), then subdivide the created edge with a ternary clique/star
        # node carrying the new leaf
        centre = classified_label[1]
        cdef int icentre = -1 if centre is None else centre
        Pstar = {m for m in self.label[u]
                 if ms[m] == _P_ and m != icentre}
        Estar = {m for m in self.label[u]
                 if ms[m] == _E_ or (m == icentre and ms[m] == _P_)}
        if self._check_mode:
            assert len(Pstar) >= 2 and len(Estar) >= 2, "case 4: not a split"
        uA, uB, qA, qB = self.node_split(u, Pstar, Estar)
        cdef dict ms2 = self.states(S, nbrs)[0]
        sA, sB = ms2[qA], ms2[qB]
        if self._check_mode:
            assert {sA, sB} <= {_P_, _E_} and _P_ in {sA, sB}, "case 4: bad states"
        if sA == _P_ and sB == _P_:
            self.subdivide(qA, x, (NodeType.CLIQUE,))
        else:
            role = _THIS_ if sA == _E_ else _OTHER_
            self.subdivide(qA, x, (NodeType.STAR, role))

    cdef _clean_and_contract(self, int x, set S, dict nbrs, dict ms):
        r"""
        Helper method for ``case 7: fully-mixed subtree`` of method ``insert``.

        This method cleans the degenerate nodes of the fully-mixed subtree,
        contract it to a single node by node-joins, and attach the new leaf
        ``x`` to a new marker adjacent to the perfect markers of that node.
        See [GPTC2016]_ for more details.

        INPUT:

        - ``x`` -- integer; a new leaf

        - ``S`` -- set; the neighbourhood among the current leaves of ``x``

        - ``nbrs`` -- dictionary; map each current leaf to its neighbours among
          the current leaves

        - ``ms`` -- dictionary; map markers to states (see method ``states``)

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.split_decomposition import GraphLabelledTree
            sage: ST = GraphLabelledTree(graphs.DiamondGraph(), check=True)
        """
        # --- cleaning ------------------------------------------------
        # nodes of the fully-mixed subtree = nodes incident to a fully-mixed
        # tree-edge
        cdef set fm_nodes = set()
        cdef set V, Pstar, Estar
        cdef int q, r, u, uA, uB, qA, qB
        cdef int icentre
        for q in ms:
            r = self.opposite[q]
            if r in self.node_of and ms[q] == _M_ and ms[r] == _M_:
                fm_nodes.add(self.node_of[q])
                fm_nodes.add(self.node_of[r])
        for u in fm_nodes:
            if u not in self.label:
                continue
            kind, centre = self.classify_label(self.label[u])
            if kind == NodeType.PRIME:
                continue
            icentre = -1 if centre is None else centre
            V = set(self.label[u])
            Pstar = {m for m in V if ms[m] == _P_ and m != icentre}
            Estar = {m for m in V
                     if ms[m] == _E_ or (m == icentre and ms[m] == _P_)}
            if len(Pstar) >= 2 and len(V - Pstar) >= 2:
                # split off the P* part, then the E* part from the rest
                uA, uB, qA, qB = self.node_split(u, Pstar, V - Pstar)
                if len(Estar) >= 2 and len(set(self.label[uB]) - Estar) >= 2:
                    self.node_split(uB, Estar, set(self.label[uB]) - Estar)
            elif len(Estar) >= 2 and len(V - Estar) >= 2:
                self.node_split(u, Estar, V - Estar)

        # --- contraction ----------------------------------------------
        cdef dict ms2 = self.states(S, nbrs)[0]
        cdef int contracted = -1
        while True:
            for q in ms2:
                if (ms2[q] == _M_ and self.opposite[q] in self.node_of
                    and ms2[self.opposite[q]] == _M_):
                    # We have the right marker q
                    break
            else:  # We are done
                break
            r = self.opposite[q]
            contracted = self.node_join(q)
            del ms2[q]
            del ms2[r]
        if self._check_mode:
            assert contracted != -1, "no fully-mixed edge"
        cdef set P_u = {m for m in self.label[contracted] if ms2[m] == _P_}
        if self._check_mode:
            assert P_u, "contracted node has no perfect marker"
        self.attach_leaf(contracted, x, P_u)

    # --------------------------------------------------------------------------
    # reducedness
    # --------------------------------------------------------------------------

    def is_reduced(self):
        r"""
        Check whether ``self`` is reduced.
        
        The graph-labelled tree ``self`` is reduced if all labels are PRIME or
        degenerate (i.e., CLIQUE, STAR, LEAF) and no clique-join or star-join
        can be performed.

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.split_decomposition import GraphLabelledTree
            sage: ST = GraphLabelledTree(graphs.DiamondGraph(), check=True)
            sage: ST.is_reduced()
            True
        """
        cdef dict info = {u: self.classify_label(label)
                          for u, label in self.label.items()}
        cdef int u, q, r
        for u in self.label:
            if info[u][0] not in NodeType:  # i.e., if not clique, star of prime
                return False
        cdef set seen_edges = set()
        cdef frozenset key
        for q, r in self.opposite.items():
            if r not in self.node_of:
                continue
            key = frozenset((q, r))
            if key in seen_edges:
                continue
            seen_edges.add(key)
            ku, cu = info[self.node_of[q]]
            kv, cv = info[self.node_of[r]]
            if ku == NodeType.CLIQUE and kv == NodeType.CLIQUE:
                return False                # clique-join possible
            if ku == NodeType.STAR and kv == NodeType.STAR:
                if (q is cu) != (r is cv):
                    return False            # star-join possible
        return True

    # --------------------------------------------------------------------------
    # description
    # --------------------------------------------------------------------------

    def _marker_repr(self, m):
        return f"q{m}"

    def _node_repr(self, u):
        return f"u{u}"

    def __repr__(self):
        r"""
        Return a string representation of ``self``.

        The representation starts with the nodes followed by the leaves.
        For each node, the string contains: the index of the node prefixed with
        ``'u'``; the kind of node (``CLIQUE``, ``STAR`` or ``PRIME``) along with
        the centre for a star; the list of edges of this node. Each internal
        marker is prefixed with ``'q'``.
        For each leaf, the string contains: the label of the leaf in the input
        graph; the marker to which it is attached; the node containing the
        marker.

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.split_decomposition import GraphLabelledTree
            sage: GraphLabelledTree(graphs.DiamondGraph(), check=True)
            u7: CLIQUE [q5-q6, q5-q8, q6-q8]
            u12: STAR (centre q9) [q9-q10, q9-q11]
            0 -- q10 -- u12
            1 -- q6 -- u7
            2 -- q8 -- u7
            3 -- q11 -- u12
            sage: GraphLabelledTree(graphs.HouseGraph(), check=True)
            u14: PRIME [q12-q15, q7-q11, q7-q12, q9-q11, q9-q12, q9-q15]
            0 -- q11 -- u14
            1 -- q7 -- u14
            2 -- q9 -- u14
            3 -- q12 -- u14
            4 -- q15 -- u14
            sage: G = graphs.BullGraph()
            sage: G.relabel(perm='abcde')
            sage: G.vertices()
            ['a', 'b', 'c', 'd', 'e']
            sage: GraphLabelledTree(G, check=True)
            u8: CLIQUE [q6-q7, q6-q9, q7-q9]
            u13: STAR (centre q11) [q10-q11, q11-q12]
            u17: STAR (centre q15) [q14-q15, q15-q16]
            a -- q6 -- u8
            b -- q11 -- u13
            c -- q15 -- u17
            d -- q12 -- u13
            e -- q16 -- u17
        """
        cdef list parts = []
        cdef list edges
        cdef int u, m, m2
        for u in self.label:
            kind, centre = self.classify_label(self.label[u])
            edges = []
            for m, s in self.label[u].items():
                m_str = self._marker_repr(m)
                for m2 in s:
                    if m < m2:
                        m2_str = self._marker_repr(m2)
                        edges.append(f"{m_str}-{m2_str}")
            extra = ""
            if kind == NodeType.STAR:
                extra = f" (centre {self._marker_repr(centre)})"
            parts.append(f"u{u}: {kind}{extra} [{', '.join(sorted(edges))}]")
        for leaf, q in self.leaf_marker.items():
            extra = f"{self.int_to_vertex[leaf]} -- {self._marker_repr(q)}"
            extra += f" -- {self._node_repr(self.node_of[q])}"
            parts.append(extra)
        return "\n".join(parts)

    # ------------------------------------------------------------------
    # Methods to build the accessibility graph
    # ------------------------------------------------------------------

    def as_graph(self, include_leaves=True, immutable=None):
        r"""
        Return a the graph-labelled tree as a sage graph.

        Each vertex corresponds to a subgraph of the decomposition and each edge
        corresponds to a split. The label of an edge indicates the pair of
        marker nodes.
        This tree can be used to reconstruct the original graph.
 
        INPUT:

        - ``immutable`` -- boolean (default: ``None``); whether to create a
          mutable/immutable graph. ``immutable=None`` (default) means that the
          accessibility graph behave the same way as the input graph.

        - ``include_leaves`` -- boolean (default: ``False``); whether to include
          the leaves (``True``) or to only return the tree connecting the
          internal nodes.

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.split_decomposition import GraphLabelledTree
            sage: G = graphs.BullGraph()
            sage: G.relabel(perm='abcde')
            sage: G.vertices()
            ['a', 'b', 'c', 'd', 'e']
            sage: ST = GraphLabelledTree(G, check=True)
            sage: T = ST.as_graph()
            sage: T.order()
            8
            sage: T.is_tree()
            True
            sage: sorted(T, key=str)
            [CLIQUE: Graph on 3 vertices,
             LEAF(a): Graph on 1 vertex,
             LEAF(b): Graph on 1 vertex,
             LEAF(c): Graph on 1 vertex,
             LEAF(d): Graph on 1 vertex,
             LEAF(e): Graph on 1 vertex,
             STAR: Graph on 3 vertices,
             STAR: Graph on 3 vertices]
            sage: T.edges(sort=True, key=str)  # random
            [(LEAF(a): Graph on 1 vertex, CLIQUE: Graph on 3 vertices, ('a', 'q6', 'LEAF')),
             (LEAF(b): Graph on 1 vertex, STAR: Graph on 3 vertices, ('b', 'q11', 'LEAF')),
             (LEAF(c): Graph on 1 vertex, STAR: Graph on 3 vertices, ('c', 'q15', 'LEAF')),
             (LEAF(d): Graph on 1 vertex, STAR: Graph on 3 vertices, ('d', 'q12', 'LEAF')),
             (LEAF(e): Graph on 1 vertex, STAR: Graph on 3 vertices, ('e', 'q16', 'LEAF')),
             (STAR: Graph on 3 vertices, CLIQUE: Graph on 3 vertices, ('q10', 'q7', None)),
             (STAR: Graph on 3 vertices, CLIQUE: Graph on 3 vertices, ('q14', 'q9', None))]
            sage: T = ST.as_graph(include_leaves=False)
            sage: T.order()
            3
            sage: T.is_tree()
            True
            sage: sorted(T, key=str)
            [CLIQUE: Graph on 3 vertices,
             STAR: Graph on 3 vertices,
             STAR: Graph on 3 vertices]
            sage: T.edges(sort=True, key=str)  # random
            [(STAR: Graph on 3 vertices, CLIQUE: Graph on 3 vertices, ('q10', 'q7', None)),
             (STAR: Graph on 3 vertices, CLIQUE: Graph on 3 vertices, ('q14', 'q9', None))]

        TESTS::

            sage: ST.as_graph(include_leaves=False).is_immutable()
            False
            sage: ST.as_graph(include_leaves=False, immutable=True).is_immutable()
            True
            sage: ST = GraphLabelledTree(G.copy(immutable=True), check=True)
            sage: ST.as_graph(include_leaves=False).is_immutable()
            True
            sage: ST.as_graph(include_leaves=False, immutable=False).is_immutable()
            False
        """
        from sage.graphs.graph import Graph
        if immutable is None:
            immutable = self._G.is_immutable()

        cdef int m, m1, m2, x, y
        cdef set Y
        cdef dict m_str = {m: self._marker_repr(m) for m in self.node_of}
        cdef dict G_label = {u: {m_str[x]: {m_str[y] for y in Y}
                                 for x, Y in adj.items()}
                             for u, adj in self.label.items()}

        # Map GLT nodes to tree-vertices
        cdef dict n_to_v = {u: Graph(adj, format='dict_of_lists', immutable=True,
                                     name=self.classify_label(adj)[0])
                            for u, adj in G_label.items()}
        if include_leaves:
            for leaf in self.leaf_marker:
                true_leaf = self.int_to_vertex[leaf]
                n_to_v[leaf] = Graph([[true_leaf], []], format="vertices_and_edges",
                                     immutable=True, name=f"LEAF({true_leaf})")

        def edges():
            for m1, m2 in self.opposite.items():
                x = n_to_v[self.node_of[m1]]
                if m2 in self.leaf_marker:
                    if include_leaves:
                        yield (x, n_to_v[m2], (self.int_to_vertex[m2], m_str[m1], "LEAF"))
                else:
                    y = n_to_v[self.node_of[m2]]
                    if x != y:
                        yield (x, y, (m_str[m1], m_str[m2], None))

        return Graph(edges(), format="list_of_edges", immutable=immutable)


    cdef set _accessible_markers(self, int m0):
        r"""
        Return the set of markers accessible from the marker ``m0``.

        This is a helper method to methof ``accessibility_graph``.

        INPUT:

        - ``m0`` -- integer; a marker

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.split_decomposition import GraphLabelledTree
            sage: ST = GraphLabelledTree(graphs.DiamondGraph(), check=True)
        """
        cdef set seen_zero = {m0}
        cdef set seen_one = {m0}
        dq = deque([(m0, 0), (m0, 1)])
        cdef set markers = {m0}
        cdef int m, m2, k, t
        while dq:
            m, k = dq.popleft()
            if k == 0:  # next move: label-edge
                for m2 in self.label[self.node_of[m]][m]:
                    if m2 not in seen_one:
                        seen_one.add(m2)
                        dq.append((m2, 1))
                        markers.add(m2)
            else:  # next move: tree-edge
                t = self.opposite[m]
                if t in self.node_of:
                    if t not in seen_zero:
                        seen_zero.add(t)
                        dq.append((t, 0))
                        markers.add(t)
        return markers

    def accessibility_graph(self, algorithm='A', immutable=None,
                            internal=False):
        r"""
        Return the accessibility graph of ``self``.

        The accessibility graph is the graph encoded in the graph-labelled
        tree. It has the vertices and edges as the input graph.

        INPUT:

        - ``algorithm`` -- string (default: ``'A'``); the algorithm to use:

          - ``'A'`` -- this algorithm first build the lists of markers
            accessible from each leaf marker and then determined valid edges
            from these data. It is used from checking the validity of the
            construction and is often faster than algorithm ``'B'``.

          - ``B`` -- this algorithm first builds the graph representation of the
            graph-labelled tree (using method ``as_graph``) and then extracts
            the graph. It helps checking the construction returned by
            ``as_graph``.

        - ``immutable`` -- boolean (default: ``None``); whether to create a
          mutable/immutable graph. ``immutable=None`` (default) means that the
          accessibility graph behave the same way as the input graph.

        - ``internal`` -- boolean (default: ``False``); whether to use internal
          vertex labels (``True``) or those of the input graph (``False``).
          This is useful for internal methods checking the construction.

        EXAMPLES::

            sage: from sage.graphs.graph_decompositions.split_decomposition import GraphLabelledTree
            sage: from sage.graphs.graph_decompositions.split_decomposition import GP2012
            sage: G = GP2012()
            sage: ST = GraphLabelledTree(G, check=True)
            sage: A = ST.accessibility_graph(algorithm='A')
            sage: A.is_isomorphic(G)
            True
            sage: set(A) == set(G)
            True
            sage: EA = set(frozenset(e) for e in A.edges(labels=False))
            sage: EA == set(frozenset(e) for e in G.edges(labels=False))
            True
            sage: B = ST.accessibility_graph(algorithm='B')
            sage: B.is_isomorphic(A)
            True
            sage: set(A) == set(B)
            True
            sage: EA == set(frozenset(e) for e in B.edges(labels=False))
            True

        TESTS::

            sage: ST.accessibility_graph(immutable=True).is_immutable()
            True
            sage: H = G.copy(immutable=True)
            sage: ST = GraphLabelledTree(H, check=True)
            sage: ST.accessibility_graph().is_immutable()
            True
            sage: ST.accessibility_graph(immutable=False).is_immutable()
            False
            sage: A = ST.accessibility_graph(algorithm='A', internal=True)
            sage: set(H) == set(A)
            False
            sage: A.is_isomorphic(H)
            True
            sage: B = ST.accessibility_graph(algorithm='B', internal=True)
            sage: set(H) == set(B)
            False
            sage: B.is_isomorphic(H)
            True
        """
        from sage.graphs.graph import Graph
        if immutable is None:
            immutable = self._G.is_immutable()

        if algorithm == 'A':
            from itertools import combinations
            acc = {leaf: self._accessible_markers(marker)
                   for leaf, marker in self.leaf_marker.items()}
            edges = ((l1, l2)
                     for (l1, _), (l2, m2) in combinations(self.leaf_marker.items(), 2)
                     if m2 in acc[l1])
            G = Graph([self.leaf_marker, edges], format="vertices_and_edges")
            perm = self.int_to_vertex
        elif algorithm == 'B':
            glt = self.as_graph(immutable=False)
            G = Graph([e for g in glt for e in g.edges()], format="list_of_edges")
            perm = {}
            for m1, m2, b in glt.edge_labels():
                if b == "LEAF":
                    perm[m2] = m1
                else:
                    X = G.neighbors(m1)
                    Y = G.neighbors(m2)
                    G.delete_vertices([m1, m2])
                    G.add_edges((x, y) for x in X for y in Y)
        else:
            raise ValueError("algorithm must be 'A' or 'B'")

        if internal:
            return G
        if immutable:
            return G.relabel(perm=perm, inplace=False, immutable=True)
        G.relabel(perm=perm)
        return G


def split_tree(G, check=True):
    r"""
    Return the split-tree (Cunningham) decomposition of ``G``.

    This method returns the split-tree of `G` as a
    :class:`~sage.graphs.graph_decompositions.split_decomposition.GraphLabelledTree`
    which implements the construction proposed in [GPTC2016]_.
    See the documentation of module
    :mod:`~sage.graphs.graph_decompositions.split_decomposition`
    for more details.

    INPUT:

    - ``G`` -- a connected graph without loops

    - ``check`` -- boolean (default: ``True``); whether to run internal tests
      and procedures for checking the consistency of the internal data
      structures and the the validity of the construction during its execution.

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.split_decomposition import split_tree
        sage: G = graphs.BullGraph()
        sage: G.relabel(perm='abcde')
        sage: G.vertices()
        ['a', 'b', 'c', 'd', 'e']
        sage: ST = split_tree(G, check=True)
        sage: A = ST.accessibility_graph(algorithm='A')
        sage: A.is_isomorphic(G)
        True
        sage: B = ST.accessibility_graph(algorithm='B')
        sage: B.is_isomorphic(G)
        True
        sage: T = ST.as_graph()
        sage: T.order()
        8
        sage: T.is_tree()
        True
        sage: sorted(T, key=str)
        [CLIQUE: Graph on 3 vertices,
         LEAF(a): Graph on 1 vertex,
         LEAF(b): Graph on 1 vertex,
         LEAF(c): Graph on 1 vertex,
         LEAF(d): Graph on 1 vertex,
         LEAF(e): Graph on 1 vertex,
         STAR: Graph on 3 vertices,
         STAR: Graph on 3 vertices]

    TESTS::

        sage: split_tree(Graph()).as_graph()
        Graph on 0 vertices
        sage: split_tree(Graph(1)).as_graph()
        Graph on 2 vertices
        sage: split_tree(Graph(2)).as_graph()
        Traceback (most recent call last):
        ...
        ValueError: the graph must be connected
        sage: split_tree(DiGraph()).as_graph()
        Traceback (most recent call last):
        ...
        ValueError: the input graph must be undirected

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
        ValueError: the input graph must be without loops
    """
    return GraphLabelledTree(G, check=check)


def GP2012():
    r"""
    Return the graph used in [GP2012b]_.

    This is useful for testing the correctness of the algorithm.

    EXAMPLE::

        sage: from sage.graphs.graph_decompositions.split_decomposition import split_tree
        sage: from sage.graphs.graph_decompositions.split_decomposition import GP2012
        sage: G = GP2012()
        sage: ST = split_tree(G, check=True)
        sage: A = ST.accessibility_graph(algorithm='A')
        sage: A.is_isomorphic(G)
        True
        sage: B = ST.accessibility_graph(algorithm='B')
        sage: B.is_isomorphic(G)
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
