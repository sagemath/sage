# distutils: language = c++
r"""
Edge connectivity

This module implements methods for computing the edge-connectivity of graphs and
digraphs. It also implements methods to extract `k` edge-disjoint spanning trees
from a `2k` edge-connected graph or a `k` edge-connected digraph, and to pack
`k` edge-disjoint spanning arborescences in a directed graph following the
matroid approach of [Gabow1995]_.

.. TODO::

    - Implement the tree-packing strengthening proposed in [BHKP2008]_
    - Extend to digraphs with multiple edges
    - Extend to weighted digraphs
"""
# ****************************************************************************
#       Copyright (c) 2022 David Coudert <david.coudert@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from memory_allocator cimport MemoryAllocator
from cysignals.signals cimport sig_check
from sage.graphs.generic_graph_pyx cimport GenericGraph_pyx
from libc.limits cimport INT_MAX
from libcpp.pair cimport pair
from libcpp.vector cimport vector
from libcpp.queue cimport queue


cdef class GabowEdgeConnectivity:
    r"""
    Gabow's algorithm for finding the edge connectivity of digraphs.

    This class implements the algorithm proposed in [Gabow1995]_ for finding the
    edge connectivity of a directed graph and `k` edge disjoint spanning trees
    if the digraph is `k` edge connected.

    .. WARNING::

        Multiple edges are currently not supported. The current implementation
        act as if the digraph is simple and so the return results might not be
        correct. We therefore raise an error if the digraph has multiple edges.

    INPUT:

    - ``D`` -- a :class:`~sage.graphs.digraph.DiGraph`

    EXAMPLES:

    A random `d`-regular digraph is `d`-edge-connected::

        sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
        sage: D = DiGraph(graphs.RandomRegular(6, 50))                                  # needs networkx
        sage: while not D.is_strongly_connected():                                      # needs networkx
        ....:     D = DiGraph(graphs.RandomRegular(6, 50))
        sage: GabowEdgeConnectivity(D).edge_connectivity()                              # needs networkx
        6

    A complete digraph with `n` vertices is `n-1`-edge-connected::

        sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
        sage: D = DiGraph(digraphs.Complete(10))
        sage: GabowEdgeConnectivity(D, use_rec=True).edge_connectivity()
        9

    Check that we get the same result when with and without the DFS-based
    speed-up initialization proposed in [GKLP2021]_::

        sage: # needs networkx
        sage: G = graphs.RandomBarabasiAlbert(100, 2)
        sage: D = DiGraph(G)
        sage: ec1 = GabowEdgeConnectivity(D,
        ....:                             dfs_preprocessing=False).edge_connectivity()
        sage: ec2 = GabowEdgeConnectivity(D,
        ....:                             dfs_preprocessing=True).edge_connectivity()
        sage: ec3 = GabowEdgeConnectivity(D, dfs_preprocessing=True,
        ....:                             use_rec=True).edge_connectivity()
        sage: ec1 == ec2 and ec2 == ec3
        True

    TESTS:

    :issue:`32169`::

        sage: dig6_string = r'[E_S?_hKIH@eos[BSg???Q@FShGC?hTHUGM?IPug?'
        sage: dig6_string += r'JOEYCdOzdkQGo@ADA@AAg?GAQW?'
        sage: dig6_string += r'[aIaSwHYcD@qQb@Dd?\hJTI@OHlJ_?C_OEIKoeCR@_BC?Q?'
        sage: dig6_string += r'?YBFosqITEA?IvCU_'
        sage: D = DiGraph(dig6_string)
        sage: GabowEdgeConnectivity(D).edge_connectivity()
        5
        sage: trees = GabowEdgeConnectivity(D).edge_disjoint_spanning_trees()
        sage: len(trees)
        5
        sage: all(T.num_edges() == D.order() - 1 for T in trees)
        True
        sage: all_edges = sum((T.edges(labels=False, sort=False) for T in trees), [])
        sage: len(all_edges) == len(set(all_edges))  # the trees are edge disjoint
        True

    Corner cases::

        sage: [GabowEdgeConnectivity(DiGraph(n)).edge_connectivity() for n in range(4)]
        [0, 0, 0, 0]
        sage: D = digraphs.Circuit(3) * 2
        sage: D.add_edge(0, 3)
        sage: GabowEdgeConnectivity(D).edge_connectivity()
        0
        sage: D.add_edge(3, 0)
        sage: GabowEdgeConnectivity(D).edge_connectivity()
        1

    Looped digraphs are supported but not digraphs with multiple edges::

        sage: D = digraphs.Complete(5, loops=True)
        sage: GabowEdgeConnectivity(D).edge_connectivity()
        4
        sage: D.allow_multiple_edges(True)
        sage: D.add_edges(D.edges(sort=False))
        sage: GabowEdgeConnectivity(D).edge_connectivity()
        Traceback (most recent call last):
        ...
        ValueError: This method is not known to work on graphs with multiedges. ...
    """
    cdef MemoryAllocator mem
    cdef Py_ssize_t n  # number of nodes
    cdef Py_ssize_t m  # number of arcs

    cdef int max_ec  # upper bound on the edge connectivity
    cdef int ec  # current (proven) value of edge connectivity
    cdef bint ec_checked  # whether we have well computed edge connectivity

    cdef int UNUSED
    cdef int FIRSTEDGE
    # edge_state sentinel for an edge temporarily pulled out of its tree
    # during a swap
    cdef int TEMP_REMOVED

    # The graph is stored as lists of incident edges
    cdef GenericGraph_pyx G  # the original graph
    cdef list int_to_vertex  # mapping from integers to vertex labels
    cdef vector[vector[int]] g_out
    cdef vector[vector[int]] g_in
    cdef vector[vector[int]] my_g  # either g_out or g_in
    cdef vector[vector[int]] my_g_reversed  # either g_out or g_in

    # values associated to edges
    cdef int* tail  # source of edge j
    cdef int* head  # target of edge j
    cdef int* my_from  # either tail or head
    cdef int* my_to  # either tail or head

    cdef int* labels  # label of each edge given by the labeling algorithm, UNUSED if unlabeled

    cdef int* edge_state_1  # index of forest Ti to which belongs the arc j of g_out, UNUSED if not used
    cdef int* edge_state_2  # index of forest Ti to which belongs the arc j of g_in, UNUSED if not used
    cdef int* my_edge_state  # either edge_state_1 or edge_state_2

    # values associated to trees and forests
    cdef int root_vertex  # 0 by default
    cdef int current_tree  # index of the current tree
    cdef int next_f_tree
    cdef int augmenting_root
    cdef bint* tree_flag  # indicate whether a tree Ti has been touched
    cdef int* root  # current root vertex of f_tree i
    cdef int* L_roots  # L_roots of the trees
    cdef bint* forests  # indicate whether the f_tree is active or inactive
    cdef bint** labeled  #

    cdef int** parent_1  # parent of v in tree/forest Ti
    cdef int** parent_2  # parent of v in tree/forest Ti
    cdef int** my_parent  # either parent_1 or parent_2
    cdef int** parent_edge_id_1  # edge id of parent of v in tree/forest Ti
    cdef int** parent_edge_id_2  # edge id  of parent of v in tree/forest Ti
    cdef int** my_parent_edge_id  # either parent_edge_id_1 or parent_edge_id_2
    cdef int** depth_1  # depth of v in tree/forest Ti
    cdef int** depth_2  # depth of v in tree/forest Ti
    cdef int** my_depth  # either depth_1 or depth_2

    # to store a path
    cdef vector[int] A_path
    cdef vector[int] left_traverse
    cdef vector[int] right_traverse

    cdef bint* seen  # for method re_init
    cdef int* stack  # stack of vertices for DFS in re_init
    cdef vector[vector[int]] tree_edges  # used to organise the edges of the trees
    cdef vector[vector[int]] F  # used to store a proven k-intersection (copy of tree_edges)
    cdef vector[vector[int]] tree_edges_incident  # lists of incident edges of a given tree

    cdef queue[int] my_Q  # queue of labeled edges
    cdef queue[pair[int, int]] joining_edges  # queue of tuples (edge id, edge state)
    cdef queue[int] incident_edges_Q  # queue of edges

    cdef int num_start_f_trees  # number of f-trees at the beginning of an iteration
    cdef int num_joins  # number of joined vertices from dfs
    cdef bint* T  # whether an edge is in the proven k-intersection
    cdef bint* visited  # for method find_dfs_tree
    cdef int * incident_edge_index  # used for DFS initialization
    cdef bint dfs_preprocessing  # whether to use DFS-based fast initialization
    cdef bint use_rec  # whether to use the recursive DFS initialization

    # State for arborescence packing (Gabow 1995). Kept separate from
    # ``self.F``, which mixes the in- and out-arborescences produced by
    # ``compute_edge_connectivity``.
    cdef vector[vector[int]] arborescence_F  # out-arborescence edge sets

    # Phase 2 (extraction) state -- used by find_tree to carve one spanning
    # out-arborescence at a time out of the k-intersection.
    cdef bint* arbor_edge      # edge is in the arborescence being extracted
    cdef bint* in_A            # vertex is in the reachable set A
    cdef bint* in_X            # vertex is exhausted (no usable out-edge)
    cdef bint* marked_edge     # edge processed in this find_tree call
    cdef int* A_order          # vertices added to A, in order (root first)
    cdef int A_size            # number of vertices currently in A
    cdef bint* extracted_edge  # edge belongs to an already-extracted tree
    cdef bint* in_union        # edge is part of the Phase 1 union
    cdef int packing_ntrees    # trees in the current partial intersection

    # Disjoint-set machinery for the dedicated Phase 2 packing search.
    # Groups vertices into sets as components merge during a swap, so the
    # fundamental-cycle traversal can jump across component boundaries on a
    # split tree.
    cdef int* myset            # set id of the vertex (-1 if none)
    cdef int* set_seen         # per-set: marked during a traversal
    cdef vector[vector[int]] root_set  # L_root of each tree, per set
    cdef int sets              # number of sets currently created

    def __init__(self, G, dfs_preprocessing=True, use_rec=False):
        r"""
        Initialize this object.

        INPUT:

        - ``G`` -- a :class:`~sage.graphs.digraph.DiGraph`

        - ``dfs_preprocessing`` -- boolean (default: ``True``); whether to use
          the DFS-based speed-up initialization proposed in [GKLP2021]_

        - ``use_rec`` -- boolean (default: ``False``); whether to use a
          recursive or non-recursive DFS for ``dfs_preprocessing``. The
          recursive DFS tends to be faster than the non-recursive version on
          complete digraphs and slower on other graphs. This parameter is
          ignored when ``dfs_preprocessing`` is ``False``.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: GabowEdgeConnectivity(D).edge_connectivity()
            4
        """
        self.dfs_preprocessing = dfs_preprocessing
        self.use_rec = use_rec
        self.ec_checked = False
        from sage.graphs.digraph import DiGraph
        if not isinstance(G, DiGraph):
            raise ValueError("this method is for directed graphs only")
        G._scream_if_not_simple(allow_loops=True)
        if G.size() > INT_MAX - 2:
            raise ValueError("the graph is too large for this code")

        # Trivial cases
        if not G or not G.is_strongly_connected():
            self.ec = 0
            self.ec_checked = True
            self.F.clear()
            return

        #
        # Initialize some data structures
        #
        self.G = <GenericGraph_pyx?>G
        self.n = G.order()
        self.m = G.size()
        self.mem = MemoryAllocator()

        # Build compact graph data structure with out and in adjacencies.
        # Loops are removed from the graph.
        self.build_graph_data_structure()
        # From now on, vertices are numbered in [0..n-1] and edges in [0..m-1]
        # where m is the number of edges after the removal of the loops

        # Set upper bound on the edge connectivity
        cdef int i, d
        self.max_ec = INT_MAX
        for i in range(self.n):
            d = self.g_out[i].size()
            if d < self.max_ec:
                self.max_ec = d
            d = self.g_in[i].size()
            if d < self.max_ec:
                self.max_ec = d

        self.labels = <int*>self.mem.calloc(self.m, sizeof(int))
        self.tree_flag = <bint*>self.mem.calloc(self.max_ec, sizeof(bint))
        self.forests = <bint*>self.mem.calloc(self.n, sizeof(bint))
        self.L_roots = <int*>self.mem.calloc(self.max_ec, sizeof(int))
        self.labeled = <bint**>self.mem.calloc(self.max_ec, sizeof(bint*))
        self.seen = <bint*>self.mem.calloc(self.n, sizeof(bint))
        self.root = <int*>self.mem.calloc(self.n, sizeof(int))
        self.edge_state_1 = <int*>self.mem.calloc(self.m, sizeof(int))
        self.edge_state_2 = <int*>self.mem.calloc(self.m, sizeof(int))
        self.parent_1 = <int**>self.mem.calloc(self.max_ec, sizeof(int*))
        self.parent_2 = <int**>self.mem.calloc(self.max_ec, sizeof(int*))
        self.parent_edge_id_1 = <int**>self.mem.calloc(self.max_ec, sizeof(int*))
        self.parent_edge_id_2 = <int**>self.mem.calloc(self.max_ec, sizeof(int*))
        self.depth_1 = <int**>self.mem.calloc(self.max_ec, sizeof(int*))
        self.depth_2 = <int**>self.mem.calloc(self.max_ec, sizeof(int*))
        self.stack = <int*>self.mem.calloc(self.n, sizeof(int))
        self.incident_edge_index = <int*>self.mem.calloc(self.n, sizeof(int))
        self.tree_edges.resize(self.max_ec)
        self.tree_edges_incident.resize(self.n)
        self.T = <bint*>self.mem.calloc(self.m, sizeof(bint))
        self.visited = <bint*>self.mem.calloc(self.n, sizeof(bint))

        # Phase 2 (arborescence extraction) buffers
        self.arbor_edge = <bint*>self.mem.calloc(self.m, sizeof(bint))
        self.marked_edge = <bint*>self.mem.calloc(self.m, sizeof(bint))
        self.in_A = <bint*>self.mem.calloc(self.n, sizeof(bint))
        self.in_X = <bint*>self.mem.calloc(self.n, sizeof(bint))
        self.A_order = <int*>self.mem.calloc(self.n, sizeof(int))
        self.extracted_edge = <bint*>self.mem.calloc(self.m, sizeof(bint))
        self.in_union = <bint*>self.mem.calloc(self.m, sizeof(bint))
        self.myset = <int*>self.mem.calloc(self.n + 2, sizeof(int))
        self.set_seen = <int*>self.mem.calloc(self.n + 2, sizeof(int))

        # Set some constants
        self.UNUSED = INT_MAX
        self.FIRSTEDGE = INT_MAX - 1
        self.TEMP_REMOVED = INT_MAX - 2

        for i in range(self.m):
            self.edge_state_1[i] = self.UNUSED  # edge i is unused
            self.edge_state_2[i] = self.UNUSED
            self.labels[i] = self.UNUSED  # edge i is unlabeled
            self.T[i] = False  # edge i doesn't belong to any k-intersection yet

        _ = self.compute_edge_connectivity()
        sig_check()

    cdef build_graph_data_structure(self):
        r"""
        Build graph data structures.

        We assign each arc (u, v) a unique id and store in arrays the tail/head
        of each arc. We use vector of vectors to quickly access incident edges.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: GabowEdgeConnectivity(D).edge_connectivity()
            4
        """
        cdef int i
        self.int_to_vertex = list(self.G)
        cdef dict vertex_to_int = {u: i for i, u in enumerate(self.int_to_vertex)}

        self.tail = <int*>self.mem.calloc(self.m, sizeof(int))
        self.head = <int*>self.mem.calloc(self.m, sizeof(int))
        self.g_out.resize(self.n)
        self.g_in.resize(self.n)
        for i in range(self.n):
            self.g_out[i].clear()
            self.g_in[i].clear()

        cdef int x, y
        cdef int e_id = 0
        for x, u in enumerate(self.int_to_vertex):
            for v in self.G.neighbor_out_iterator(u):
                y = vertex_to_int[v]
                if x != y:
                    self.g_out[x].push_back(e_id)
                    self.g_in[y].push_back(e_id)
                    self.tail[e_id] = x
                    self.head[e_id] = y
                    e_id += 1
        # Loops have been removed, so we update the number of edges
        self.m = e_id

    cdef bint compute_edge_connectivity(self) except -1:
        """
        Compute the edge connectivity using Round Robin algorithm.

        The method returns ``True`` if the computation ends normally. Otherwise
        an exception is raised, for instance due to a keyboard interruption.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: GabowEdgeConnectivity(D).edge_connectivity()
            4
        """
        cdef int i

        self.root_vertex = 0
        self.next_f_tree = 0

        # Search successively trees in g_in and g_out
        self.ec = 0
        for i in range(self.max_ec):
            if self.construct_trees(False, i) and self.construct_trees(True, i):
                # We found both an in-arborescence and an out-arborescence.
                # So we can increase the edge connectivity
                self.ec += 1
                # and save the current k-intersection
                self.save_current_k_intersection()
            sig_check()
        self.ec_checked = True
        return True

    cdef bint construct_trees(self, bint reverse, int tree) except -1:
        r"""
        Search for an in or out arborescence.

        INPUT:

        - ``reverse`` -- boolean; whether to search for an in-arborescence
          (``True``) or an out-arborescence (``False``)

        - ``tree`` -- integer; index of the tree

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: GabowEdgeConnectivity(D).edge_connectivity()
            4
        """
        if reverse:
            # Search for a spanning tree in g-reversed
            self.my_g = self.g_out
            self.my_g_reversed = self.g_in
            self.my_from = self.head
            self.my_to = self.tail
            self.my_parent = self.parent_2
            self.my_depth = self.depth_2
            self.my_parent_edge_id = self.parent_edge_id_2
            self.my_edge_state = self.edge_state_2
        else:
            # Search for a spanning tree in g using incoming arcs
            self.my_g = self.g_in
            self.my_g_reversed = self.g_out
            self.my_from = self.tail
            self.my_to = self.head
            self.my_parent = self.parent_1
            self.my_depth = self.depth_1
            self.my_parent_edge_id = self.parent_edge_id_1
            self.my_edge_state = self.edge_state_1

        self.current_tree = tree
        self.increase_memory_for_new_tree(tree)

        cdef int njoins = 0
        cdef int z

        self.num_start_f_trees = self.n - self.num_joins
        # There are fewer than n f-trees. We prepare to join them. If there's
        # only one f-tree, we just save the edges and advance to the next
        # iteration
        if self.dfs_preprocessing and self.num_start_f_trees < self.n - 1:
            self.re_init(tree)

        # There are n f-trees, and we try to join them
        while njoins < self.num_start_f_trees-1:
            # Get the root of an active subtree or INT_MAX if none exists
            z = self.choose_root()
            while z != INT_MAX:
                if self.search_joining(z):
                    # We have augmented the root of the corresponding f_tree
                    njoins += 1
                else:
                    # We cannot find a tree
                    return False

                z = self.choose_root()

            # Trace the paths in order to transfer the edges to the appropriate
            # tree Ti
            self.augmentation_algorithm()
            # Reinitialize data structures and make all f_trees active for next round
            self.re_init(tree)
            sig_check()

        return True

    cdef void increase_memory_for_new_tree(self, int tree) noexcept:
        """
        Allocate data structure for the new tree/forest.

        This method also initializes data structures for this tree index. Data
        structures for a given tree index are allocated only once.

        INPUT:

        - ``tree`` -- integer; index of the tree

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: GabowEdgeConnectivity(D).edge_connectivity()
            4
        """
        if not self.labeled[tree]:
            self.labeled[tree] = <bint*>self.mem.calloc(self.n, sizeof(bint))
        if not self.my_parent[tree]:
            self.my_parent[tree] = <int*>self.mem.calloc(self.n, sizeof(int))
        if not self.my_depth[tree]:
            self.my_depth[tree] = <int*>self.mem.calloc(self.n, sizeof(int))
        if not self.my_parent_edge_id[tree]:
            self.my_parent_edge_id[tree] = <int*>self.mem.calloc(self.n, sizeof(int))

        cdef int j
        for j in range(self.n):
            self.my_parent[tree][j] = 0
            self.my_parent_edge_id[tree][j] = self.UNUSED
            self.my_depth[tree][j] = 0
            self.labeled[tree][j] = False
            self.root[j] = j
            self.forests[j] = True

        self.num_joins = 0

        # Initialize T_k to be a DFS spanning forest of G \ T
        if self.dfs_preprocessing:
            self.compute_dfs_tree()

        # Set inactive the f-trees of the root vertex
        self.forests[self.root_vertex] = False

        self.L_roots[tree] = self.UNUSED
        self.tree_flag[tree] = False

    cdef void compute_dfs_tree(self) noexcept:
        r"""
        Find a DFS spanning forest of `G \backslash T`.

        This is the DFS-based speed-up initialization proposed in [GKLP2021]_.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: GabowEdgeConnectivity(D, dfs_preprocessing=True).edge_connectivity()
            4
            sage: GabowEdgeConnectivity(D, dfs_preprocessing=False).edge_connectivity()
            4
        """
        # Mark all vertices as unvisited
        cdef int i
        for i in range(self.n):
            self.visited[i] = False

        cdef int r
        for r in range(self.n):
            if not self.visited[r]:
                # Make this vertex the root of the following dfs tree
                self.root[r] = r
                # Make the f-tree rooted at this vertex active
                self.forests[r] = True
                # Find connected vertices of the f-tree rooted at r
                if self.use_rec:
                    self.find_dfs_tree_rec(r, r)
                else:
                    self.find_dfs_tree(r)
                # Each call of find_dfs_tree creates an f-tree
                self.num_start_f_trees += 1

    cdef void find_dfs_tree(self, int r) noexcept:
        r"""
        Find more vertices of the f-tree rooted at `r`.

        This is part of the DFS-based speed-up initialization proposed in
        [GKLP2021]_.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: GabowEdgeConnectivity(D, dfs_preprocessing=True, use_rec=False).edge_connectivity()
            4
        """
        cdef int u, v, e_id
        cdef int * stack = self.stack
        cdef int * edge_index = self.incident_edge_index
        # We initialize the stack and mark the root as visited
        cdef int t = 0  # index pointing to the top of the stack
        stack[0] = r
        edge_index[r] = self.my_g_reversed[r].size()
        self.visited[r] = True

        while t >= 0:
            u = stack[t]
            if edge_index[u]:
                # Visit the next incident edge of u
                edge_index[u] -= 1
                e_id = self.my_g_reversed[u][edge_index[u]]
                if not self.T[e_id]:
                    v = self.my_to[e_id]
                    if not self.visited[v] and v != self.root_vertex:
                        # Make v belong to the f-tree rooted at r
                        self.root[v] = r
                        self.forests[v] = False
                        self.my_edge_state[e_id] = self.current_tree
                        self.num_joins += 1
                        self.visited[v] = True
                        # add v to the stack
                        t += 1
                        stack[t] = v
                        edge_index[v] = self.my_g_reversed[v].size()
                        # and proceed with v
            else:
                # We are done with u. We pop.
                t -= 1

    cdef void find_dfs_tree_rec(self, int u, int r) noexcept:
        r"""
        Find more vertices of the f-tree rooted at `r`.

        This is part of the DFS-based speed-up initialization proposed in
        [GKLP2021]_.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: GabowEdgeConnectivity(D, dfs_preprocessing=True, use_rec=True).edge_connectivity()
            4
        """
        # Mark vertex u as visited to avoid visiting it multiple times
        self.visited[u] = True

        # Visit outgoing arcs of current vertex
        cdef int e_id, v
        for e_id in self.my_g_reversed[u]:
            v = self.my_to[e_id]
            # Ensure a vertex is not visited, is not a proven k-intersection edge
            # and root_vertex remains deficient
            if not self.visited[v] and not self.T[e_id] and v != self.root_vertex:
                # Make current vertex belong to the f_tree rooted at r
                self.root[v] = r
                self.forests[v] = False
                self.my_edge_state[e_id] = self.current_tree
                self.num_joins += 1
                # recursively find more vertices and grow the subtree rooted at r
                self.find_dfs_tree_rec(v, r)

    cdef int choose_root(self) noexcept:
        """
        Return the root of an active f_tree, or INT_MAX if none exists.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: GabowEdgeConnectivity(D).edge_connectivity()
            4
        """
        cdef int v
        cdef int i

        for i in range(self.next_f_tree, self.n):
            v = self.root[i]
            if self.forests[v]:
                # this forest is active
                self.next_f_tree = i + 1
                return v
        return INT_MAX

    cdef bint search_joining(self, int x) except -1:
        """
        Try to augment the f_tree rooted at x.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: GabowEdgeConnectivity(D).edge_connectivity()
            4
        """
        cdef int y
        cdef int joining_edge
        cdef int e_id, ep

        # Store the vertex that is about to be augmented
        self.augmenting_root = x

        # Consider the incoming arcs of x
        for e_id in self.my_g[x]:
            y = self.my_from[e_id]
            # find the root of the f_tree
            y = self.root[y]

            if self.my_edge_state[e_id] == self.UNUSED:
                # The edge is available
                if x != y:
                    # ... and the f_trees have different roots. We set the
                    # label of edges in the queue to UNUSED and clear the queue
                    while not self.my_Q.empty():
                        ep = self.my_Q.front()
                        self.my_Q.pop()
                        self.labels[ep] = self.UNUSED
                    # We then assign the edge to the current_tree
                    self.join(e_id)
                    return True
                else:
                    # The f_trees have the same root (cycle).
                    # We add the edge to the queue
                    self.my_Q.push(e_id)
                    # and indicate the first edge of the path
                    self.labels[e_id] = self.FIRSTEDGE

        # If we did not find a free joining edge, we check for a sequence of
        # swaps in order to free a joining edge

        # Initialize the L_i tree of every T_i with vertex x and make x labeled
        cdef int i
        for i in range(self.current_tree + 1):
            self.L_roots[i] = x
            self.labeled[i][x] = True

        # Start cycle_scanning algorithm
        joining_edge = self.next_joining_edge_step()
        sig_check()

        if joining_edge != INT_MAX:
            # We found a joining edge
            self.joining_edges.push((joining_edge, self.my_edge_state[joining_edge]))
            self.join(joining_edge)
            return True
        return False

    cdef void join(self, int e_id) noexcept:
        """
        Assign edge e_id to current tree.

        This method joins 2 f_trees and updates the root of the new f_tree.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: GabowEdgeConnectivity(D).edge_connectivity()
            4
        """
        cdef int x = self.my_from[e_id]
        cdef int y = self.my_to[e_id]
        cdef int root_x = self.root[x]
        cdef int root_y = self.root[y]

        # Add the edge to the current tree
        self.my_edge_state[e_id] = self.current_tree

        # Make the 2 joined f_trees inactive
        self.forests[root_x] = False
        self.forests[root_y] = False

        # Update the root of the joining f_tree
        if self.augmenting_root == root_y:
            self.root[root_y] = self.root[root_x]
        else:
            self.root[root_x] = self.root[root_y]

        # Empty the queue
        while not self.my_Q.empty():
            self.my_Q.pop()

    cdef int next_joining_edge_step(self) except -1:
        """
        Process edges in the queue and start labeling until the queue is empty
        or a joining edge is found.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: GabowEdgeConnectivity(D).edge_connectivity()
            4
        """
        cdef int e_id
        cdef int found_joining
        cdef int tree = 0

        while not self.my_Q.empty():
            e_id = self.my_Q.front()
            self.my_Q.pop()

            if self.my_edge_state[e_id] == tree:
                # edge e_id is in Ti
                tree += 1
                if tree > self.current_tree:
                    tree = 0

            self.tree_flag[tree] = True

            # Search for the fundamental cycle of e_id in Ti
            found_joining = self.fundamental_cycle_step(e_id, tree)
            sig_check()
            if found_joining != INT_MAX:
                return found_joining

        return INT_MAX

    cdef int fundamental_cycle_step(self, int e_id, int tree) except -1:
        """
        Traverse tree paths from the endpoints of edge e_id to build A_path

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: GabowEdgeConnectivity(D).edge_connectivity()
            4
        """
        cdef int x = self.my_to[e_id]
        cdef int y = self.my_from[e_id]
        cdef bint left_first = True

        if self.labeled[tree][x]:
            # Node x is labeled. We go to the root of Li
            x = self.L_roots[tree]
        elif not self.labeled[tree][y]:
            raise ValueError("error in labeling")
        if self.labeled[tree][y]:
            # Node y is labeled. We go to the root of Li
            y = self.L_roots[tree]
            left_first = False
        if x == y:
            # The fundamental cycle contains no unlabeled edge
            return INT_MAX

        cdef bint stop = False
        cdef int q
        cdef vector[int] left_traverse
        cdef vector[int] right_traverse
        left_traverse.clear()
        right_traverse.clear()

        # Start double traversal
        while True:

            while self.my_depth[tree][x] >= self.my_depth[tree][y]:
                self.labeled[tree][x] = True
                q = self.my_parent_edge_id[tree][x]
                if q == self.UNUSED:
                    raise ValueError("did not find the right edge")
                # We check if edge q is unlabeled
                if self.labels[q] == self.UNUSED:
                    # If so, we place it in left_traverse array
                    left_traverse.push_back(q)
                    if self.is_joining_edge(q):
                        self.labels[q] = e_id
                        return q
                    x = self.my_parent[tree][x]
                else:
                    # Otherwise, we stop
                    stop = True
                    break
                if x == y:
                    break

            while self.my_depth[tree][y] > self.my_depth[tree][x]:
                self.labeled[tree][y] = True
                q = self.my_parent_edge_id[tree][y]
                if q == self.UNUSED:
                    raise ValueError("did not find the right edge")
                # We check if edge q is unlabeled
                if self.labels[q] == self.UNUSED:
                    # If so, we place it in right_traverse array
                    right_traverse.push_back(q)
                    if self.is_joining_edge(q):
                        self.labels[q] = e_id
                        return q
                    y = self.my_parent[tree][y]
                else:
                    # Otherwise, we stop
                    stop = True
                    break

            if x == y or stop:
                break

        if x == y:
            # Update the L_root of the tree
            self.L_roots[tree] = x
            self.labeled[tree][x] = True

        # Compute A_path
        self.A_path.clear()
        if left_first:
            for x in left_traverse:
                self.A_path.push_back(x)
            for x in range(right_traverse.size() - 1, -1, -1):
                self.A_path.push_back(right_traverse[x])
        else:
            for x in right_traverse:
                self.A_path.push_back(x)
            for x in range(left_traverse.size() - 1, -1, -1):
                self.A_path.push_back(left_traverse[x])

        return self.label_A_path(e_id)

    cdef bint is_joining_edge(self, int e_id) noexcept:
        """
        Check if edge e_id is joining.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: GabowEdgeConnectivity(D).edge_connectivity()
            4
        """
        cdef int root_x = self.root[self.my_from[e_id]]
        cdef int root_y = self.root[self.my_to[e_id]]
        return (root_x != root_y) and (root_x == self.augmenting_root or root_y == self.augmenting_root)

    cdef int label_A_path(self, int e_id) noexcept:
        """
        Labels the incident unused edges as the label_A_step of the algorithm

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: GabowEdgeConnectivity(D).edge_connectivity()
            4
        """
        cdef int e, ep

        for e in self.A_path:
            # Run label step with edge e and label e_id
            if self.label_step(e, e_id):
                return e

            if self.any_unused_is_unlabeled(self.my_to[e]):
                while not self.incident_edges_Q.empty():
                    ep = self.incident_edges_Q.front()
                    self.incident_edges_Q.pop()
                    if e != ep:
                        # Label each unused and unlabeled edge ep with e
                        if self.label_step(ep, e):
                            while not self.incident_edges_Q.empty():
                                self.incident_edges_Q.pop()
                            return ep

        while not self.incident_edges_Q.empty():
            self.incident_edges_Q.pop()

        return INT_MAX

    cdef bint label_step(self, int e_id, int e_label) noexcept:
        """
        Label edge e_id with e_label and check whether edge e_id is joining.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: GabowEdgeConnectivity(D).edge_connectivity()
            4
        """
        self.labels[e_id] = e_label

        cdef int root_x = self.root[self.my_from[e_id]]
        cdef int root_y = self.root[self.my_to[e_id]]

        if root_x == root_y:
            self.my_Q.push(e_id)
            return False
        # The roots are different. Check whether one of them is on the f_tree
        return root_x == self.augmenting_root or root_y == self.augmenting_root

    cdef bint any_unused_is_unlabeled(self, int x) noexcept:
        """
        Check if each unused edge directed to x is unlabeled

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: GabowEdgeConnectivity(D).edge_connectivity()
            4
        """
        cdef int e_id
        for e_id in self.my_g[x]:
            if self.my_edge_state[e_id] == self.UNUSED:
                if self.labels[e_id] != self.UNUSED:
                    return False
                self.incident_edges_Q.push(e_id)

        return True

    cdef void augmentation_algorithm(self) noexcept:
        """
        Trace the path of the found joining edges

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: GabowEdgeConnectivity(D).edge_connectivity()
            4
        """
        cdef int e_id, e_state
        while not self.joining_edges.empty():
            e_id, e_state = self.joining_edges.front()
            self.joining_edges.pop()
            self.trace_back(e_id, e_state)

    cdef void trace_back(self, int e_id, int e_state) noexcept:
        """
        Trace the path of a joining edge and transfer the edges to the
        appropriate tree Ti.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: GabowEdgeConnectivity(D).edge_connectivity()
            4
        """
        # Previous state (tree Ti or unused) of an edge
        cdef int previous_state = self.FIRSTEDGE

        cdef int tree
        cdef int e = self.labels[e_id]
        cdef int ep = self.labels[e]

        if e_state == self.UNUSED:
            tree = self.my_edge_state[e]
            previous_state = self.my_edge_state[ep]

            # Transfer edge ep to tree Ti and remove edge e
            self.my_edge_state[ep] = tree
            self.my_edge_state[e] = self.UNUSED

            e = ep
            ep = self.labels[e]
        else:
            tree = e_state + 1
            if tree > self.current_tree:
                tree = 0
            e = e_id
            ep = self.labels[e]

        # Transfer edges to the appropriate Ti
        while ep != self.FIRSTEDGE:
            tree -= 1
            if tree < 0:
                tree = self.current_tree

            if previous_state == self.UNUSED:
                e = ep
                ep = self.labels[e]
                self.my_edge_state[e] = self.UNUSED

            previous_state = self.my_edge_state[ep]
            self.my_edge_state[ep] = tree
            e = ep
            ep = self.labels[e]

    cdef re_init(self, int tree):
        """
        Make f_trees active (except the f_tree of the root), update depths and
        parent values, and clear the labels.

        This method is called at the end of each round of method
        construct_trees, right after the call to augmentation_algorithm.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: GabowEdgeConnectivity(D).edge_connectivity()
            4
        """
        cdef int i, j
        for j in range(self.m):
            self.labels[j] = self.UNUSED

        # Arrange the edges of each tree
        for j in range(tree + 1):
            if self.dfs_preprocessing:
                for e_id in self.tree_edges[j]:
                    self.T[e_id] = False
            self.tree_edges[j].clear()
        for j in range(self.m):
            if self.my_edge_state[j] != self.UNUSED:
                self.tree_edges[self.my_edge_state[j]].push_back(j)
                if self.dfs_preprocessing:
                    self.T[j] = True

        for j in range(tree + 1):
            if not j or j == tree or self.tree_flag[j]:
                # Build adjacency lists of incident edges (ignore direction)
                for i in range(self.n):
                    self.tree_edges_incident[i].clear()
                for i in self.tree_edges[j]:
                    self.tree_edges_incident[self.my_from[i]].push_back(i)
                    self.tree_edges_incident[self.my_to[i]].push_back(i)

                self.update_parents_depths(j)

        for i in range(tree + 1):
            self.L_roots[i] = self.UNUSED  # clear the root of each Li
            self.tree_flag[i] = False

            # Unlabel all nodes from every Ti
            for j in range(self.n):
                self.labeled[i][j] = False

        self.next_f_tree = 0

        # Finally, set active the roots of subtrees
        for i in range(self.n):
            j = self.root[i]
            if j != self.root_vertex:
                self.forests[j] = True

    cdef void update_parents_depths(self, int tree) noexcept:
        """
        Update parents, depths, and, if current_tree is k, the vertex labels to
        the root of each f_tree.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: GabowEdgeConnectivity(D).edge_connectivity()
            4
        """
        cdef int i, v

        for i in range(self.n):
            self.my_parent[tree][i] = i
            self.my_depth[tree][i] = 0
            self.seen[i] = False

        self.update_parents_dfs(tree, self.root_vertex)

        if tree == self.current_tree:
            for i in range(self.n):
                v = self.root[i]
                if self.root[v] != v:
                    v = self.root[v]
                if not self.seen[v]:
                    self.update_parents_dfs(tree, v)
                self.root[i] = self.root[v]

    cdef void update_parents_dfs(self, int tree, int x) noexcept:
        """
        Helper method for ``update_parents_depths``.

        This method updates parents and depths in specified ``tree`` starting
        from vertex ``x`` in depth first search manner.

        INPUT:

        - ``tree`` -- integer; index of the tree in which to update data

        - ``x`` -- integer; vertex from which to start the DFS

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: GabowEdgeConnectivity(D).edge_connectivity()
            4
        """
        cdef int u, v, e_id
        cdef int depth
        cdef int i = 1
        self.stack[0] = x
        self.seen[x] = True

        while i > 0:
            i -= 1
            u = self.stack[i]
            depth = self.my_depth[tree][u] + 1
            for e_id in self.tree_edges_incident[u]:
                v = self.my_to[e_id]
                if v == u:
                    v = self.my_from[e_id]
                if not self.seen[v]:
                    self.stack[i] = v
                    i += 1
                    self.seen[v] = True
                    self.my_parent[tree][v] = u
                    self.my_parent_edge_id[tree][v] = e_id
                    self.my_depth[tree][v] = depth

    cdef void save_current_k_intersection(self) noexcept:
        """
        Save the current k-intersection.

        This method is called each time the upper bound on the edge connectivity
        has been increased. The k-intersection will be used to extract the
        edge-disjoint spanning trees. If asking for the edge connectivity only,
        there is no need to call this method.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: GabowEdgeConnectivity(D).edge_connectivity()
            4
        """
        cdef int i, j
        cdef int size = self.tree_edges.size()
        self.F.resize(size)
        for i in range(size):
            self.F[i].clear()
            for j in self.tree_edges[i]:
                self.F[i].push_back(j)

    def edge_connectivity(self):
        """
        Return the edge connectivity of the digraph.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: GabowEdgeConnectivity(D).edge_connectivity()
            4
        """
        if self.ec_checked:
            return self.ec
        raise ValueError("the value of the edge connectivity has not been "
                         "properly computed. This may result from an interruption")

    #
    # Packing arborescences
    #

    cdef void rebuild_pool_adjacency(self) noexcept:
        r"""
        Rebuild ``g_in`` / ``g_out`` to contain only the *pool* edges:
        edges in the Phase 1 union that have not yet been extracted into a
        previous arborescence.

        Instead of recreating a smaller graph, we rebuild the two adjacency
        vectors in place so that both
        :meth:`construct_trees` (which reads ``g_in`` / ``g_out``) and
        :meth:`find_tree` (which reads ``g_out``) only ever see pool edges.
        Edge ids, ``tail`` and ``head`` are unchanged.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: len(GabowEdgeConnectivity(D).edge_disjoint_spanning_trees(k=4))
            4
        """
        cdef int e
        self.g_in = vector[vector[int]](self.n)
        self.g_out = vector[vector[int]](self.n)
        for e in range(self.m):
            if self.in_union[e] and not self.extracted_edge[e]:
                self.g_out[self.tail[e]].push_back(e)
                self.g_in[self.head[e]].push_back(e)

    cdef bint compute_arborescence_packing(self, int k, int root_idx) except -1:
        r"""
        Build ``k`` edge-disjoint out-arborescences rooted at ``root_idx``.

        Two phases, following [Gabow1995]_ and [GKMN2022]_:

        - **Phase 1** builds the full ``k``-intersection to determine the
          *union* (the ``k``-in-regular subgraph rooted at ``root_idx``).
        - **Phase 2** decomposes the union into ``k`` valid spanning
          out-arborescences. For each of the ``k`` rounds with
          ``current_k`` trees still to extract, it builds a
          ``(current_k - 1)``-intersection on the remaining pool, then
          calls :meth:`find_tree` to carve out one arborescence using the
          leftover free edges plus swaps, and removes it from the pool.

        INPUT:

        - ``k`` -- integer; number of arborescences to pack

        - ``root_idx`` -- integer; internal index of the root vertex (i.e.
          the position of the root in ``self.int_to_vertex``)

        OUTPUT:

        ``True`` if ``k`` edge-disjoint spanning out-arborescences rooted at
        ``root_idx`` were successfully constructed; ``False`` otherwise.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: trees = GabowEdgeConnectivity(D).edge_disjoint_spanning_trees(k=4)
            sage: len(trees)
            4
        """
        if k <= 0:
            self.arborescence_F.clear()
            return True
        if k > self.max_ec:
            return False

        cdef int i, j, t, round_i, current_k

        # Reset out-direction edge state, tree edge lists, edge labels,
        # the DFS-tree edge marker and the Phase 2 markers, so this method
        # is independent of any prior run.
        for j in range(self.m):
            self.edge_state_1[j] = self.UNUSED
            self.labels[j] = self.UNUSED
            self.T[j] = False
            self.extracted_edge[j] = False
            self.in_union[j] = False
        cdef int sz = <int>self.tree_edges.size()
        for j in range(sz):
            self.tree_edges[j].clear()

        self.root_vertex = root_idx
        self.next_f_tree = 0

        # The DFS-based speed-up of [GKLP2021]_ is disabled for the packing.
        # It warm-starts the matroid-intersection build (Phase 1 and the
        # per-round rebuilds), but profiling shows that step is only ~5% of the
        # packing time: the cost is dominated (>90%) by the Phase 2 extraction
        # (find_tree and its swaps), which the DFS preprocessing does not
        # touch. Enabling it (via a root->vertex-0 relabeling) gave no
        # measurable speed-up, so it is kept off to stay simple and robust.
        cdef bint saved_dfs = self.dfs_preprocessing
        self.dfs_preprocessing = False

        # ----- Phase 1: build the k-intersection (find the union) -----
        cdef bint ok = True
        for i in range(k):
            if not self.construct_trees(False, i):
                ok = False
                break
            sig_check()

        if not ok:
            self.dfs_preprocessing = saved_dfs
            return False

        # Record the union: every edge assigned to some tree of the
        # k-intersection. This is the k-in-regular subgraph (rooted at
        # root_idx) that we now decompose into k arborescences.
        for j in range(self.m):
            self.in_union[j] = (self.edge_state_1[j] != self.UNUSED)

        # ----- Phase 2: decompose the union into k arborescences -----
        # Restrict the working graph to the union pool. We save the full
        # adjacency and restore it at the end so the instance stays
        # reusable.
        cdef vector[vector[int]] saved_g_in = self.g_in
        cdef vector[vector[int]] saved_g_out = self.g_out
        self.rebuild_pool_adjacency()

        self.arborescence_F.clear()
        self.arborescence_F.resize(k)

        for round_i in range(k):
            current_k = k - round_i

            # Release the pool: clear the assignment/labels of every
            # not-yet-extracted union edge so we can rebuild a fresh
            # (current_k - 1)-intersection on what remains.
            for j in range(self.m):
                self.labels[j] = self.UNUSED
                self.T[j] = False
                if self.in_union[j] and not self.extracted_edge[j]:
                    self.edge_state_1[j] = self.UNUSED
            self.next_f_tree = 0

            # Build a (current_k - 1)-intersection on the remaining pool.
            # When current_k == 1 this loop is empty: the pool then holds
            # exactly one in-edge per non-root vertex, i.e. the final
            # arborescence itself.
            for t in range(current_k - 1):
                if not self.construct_trees(False, t):
                    ok = False
                    break
            if not ok:
                break
            self.packing_ntrees = current_k - 1

            # Carve out one spanning arborescence (free pool edges + swaps).
            if not self.find_tree():
                ok = False
                break

            # Record it and remove it from the pool (G_minus_A).
            self.arborescence_F[round_i].clear()
            for j in range(self.m):
                if self.arbor_edge[j]:
                    self.arborescence_F[round_i].push_back(j)
                    self.extracted_edge[j] = True
            self.rebuild_pool_adjacency()
            sig_check()

        # Restore the full adjacency and the DFS flag.
        self.g_in = saved_g_in
        self.g_out = saved_g_out
        self.dfs_preprocessing = saved_dfs
        return ok

    cdef void save_arborescence_packing_intersection(self, int k) noexcept:
        r"""
        Snapshot the first ``k`` trees from ``tree_edges`` into
        ``arborescence_F``.

        This is the packing analog of :meth:`save_current_k_intersection`,
        but it preserves only the out-arborescence edges produced by
        :meth:`compute_arborescence_packing` and leaves ``self.F``
        untouched so the connectivity computation can use it independently.

        Must be called once, after all ``k`` invocations of
        ``construct_trees(False, i)`` have completed: matroid-intersection
        augmentation can rearrange edges between earlier trees during the
        construction of later ones, so only the final state of
        ``tree_edges`` is a consistent ``k``-intersection.

        INPUT:

        - ``k`` -- integer; number of trees to snapshot

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: len(GabowEdgeConnectivity(D).edge_disjoint_spanning_trees(k=1))
            1
        """
        cdef int i, e_id
        self.arborescence_F.resize(k)
        for i in range(k):
            self.arborescence_F[i].clear()
            for e_id in self.tree_edges[i]:
                self.arborescence_F[i].push_back(e_id)

    cdef void clear_augmentation_aux_state(self) noexcept:
        r"""
        Reset the auxiliary scratch state used by the matroid-intersection
        augmentation (``labels``, ``my_Q``, ``joining_edges``,
        ``incident_edges_Q``, per-tree ``labeled[i][v]``, ``L_roots[i]``,
        ``tree_flag[i]``).

        Must be called by :meth:`search_step` after each augmentation
        attempt (success or failure) so the next attempt starts from a
        clean scratchpad. Does NOT touch the intersection state itself
        (``edge_state_1``, ``parent_1``, ``depth_1``, ``root``, ``forests``,
        ``tree_edges``); those are owned by Phase 1 / by
        :meth:`search_step`'s own save/restore protocol.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: len(GabowEdgeConnectivity(D).edge_disjoint_spanning_trees(k=4))
            4
        """
        cdef int i, j

        for j in range(self.m):
            self.labels[j] = self.UNUSED

        while not self.my_Q.empty():
            self.my_Q.pop()
        while not self.joining_edges.empty():
            self.joining_edges.pop()
        while not self.incident_edges_Q.empty():
            self.incident_edges_Q.pop()

        for i in range(self.current_tree + 1):
            for j in range(self.n):
                self.labeled[i][j] = False
            self.L_roots[i] = self.UNUSED
            self.tree_flag[i] = False

    cdef void packing_init_set(self) noexcept:
        r"""
        Reset the disjoint-set machinery for a new packing search.
        One empty set (id 0) is created with no per-tree root yet, and
        every vertex starts unassigned.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: len(GabowEdgeConnectivity(D).edge_disjoint_spanning_trees(k=4))
            4
        """
        cdef int i
        self.sets = 0
        self.root_set.clear()
        self.root_set.resize(1)
        self.root_set[0].assign(self.packing_ntrees, -1)
        for i in range(self.n + 2):
            self.myset[i] = -1
            self.set_seen[i] = 0

    cdef void packing_merge_sets(self) noexcept:
        r"""
        Close the current set and open a fresh one. Every vertex whose set
        was marked ``set_seen`` during the last traversal is folded into
        the new set; the current ``L_roots`` are stored as the new set's
        per-tree roots.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: len(GabowEdgeConnectivity(D).edge_disjoint_spanning_trees(k=4))
            4
        """
        cdef int i, s
        if self.sets != 0:
            for i in range(self.n):
                s = self.myset[i]
                if s != -1 and self.set_seen[s] == 1:
                    self.myset[i] = self.sets
        for i in range(self.packing_ntrees):
            # Translate our UNUSED (INT_MAX) "empty L_root" marker to -1, the
            # sentinel the cycle-scan checks against. Storing INT_MAX here
            # would later be used as a vertex index -> out-of-bounds crash.
            if self.L_roots[i] == self.UNUSED:
                self.root_set[self.sets][i] = -1
            else:
                self.root_set[self.sets][i] = self.L_roots[i]
        for i in range(self.n):
            self.set_seen[i] = 0
        self.sets += 1
        self.root_set.push_back(vector[int]())
        self.root_set[self.sets].assign(self.packing_ntrees, -1)

    cdef bint packing_check_if_joining(self, int j) noexcept:
        r"""
        Return whether edge ``j`` is a joining edge for the current packing
        search: its endpoints lie in different f_trees, one of which is the
        augmenting root, and ``j`` is not already in the arborescence.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: len(GabowEdgeConnectivity(D).edge_disjoint_spanning_trees(k=4))
            4
        """
        cdef int z1 = self.root[self.head[j]]
        cdef int z2 = self.root[self.tail[j]]
        if z1 == self.augmenting_root or z2 == self.augmenting_root:
            if z1 != z2 and not self.arbor_edge[j]:
                return True
        return False

    cdef bint packing_label_step(self, int j, int e) noexcept:
        r"""
        Label edge ``j`` with ``e`` and report whether ``j`` is joining.
        Non-joining edges are pushed onto the work queue ``my_Q``.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: len(GabowEdgeConnectivity(D).edge_disjoint_spanning_trees(k=4))
            4
        """
        cdef int z1 = self.root[self.head[j]]
        cdef int z2 = self.root[self.tail[j]]
        self.labels[j] = e
        if z1 != z2 and not self.arbor_edge[e]:
            if z1 == self.augmenting_root or z2 == self.augmenting_root:
                return True
        else:
            self.my_Q.push(j)
        return False

    cdef bint packing_any_unused_is_unlabelled(self, int x) noexcept:
        r"""
        Check whether every unused, non-arborescence edge directed into
        ``x`` is also unlabeled; collect such edges in
        ``incident_edges_Q``.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: len(GabowEdgeConnectivity(D).edge_disjoint_spanning_trees(k=4))
            4
        """
        cdef int e
        for e in self.g_in[x]:
            if self.edge_state_1[e] == self.UNUSED and not self.arbor_edge[e]:
                if self.labels[e] != self.UNUSED:
                    return False
                self.incident_edges_Q.push(e)
        return True

    cdef int packing_label_A_step(self, int l, int k) noexcept:
        r"""
        Label the incident unused edges along the current ``A_path``.
        ``l`` is the label to apply, ``k`` the number of edges in
        ``A_path``.

        Returns the edge id of a joining edge if found, else ``INT_MAX``.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: len(GabowEdgeConnectivity(D).edge_disjoint_spanning_trees(k=4))
            4
        """
        cdef int i = 0, g, x, edge

        while i < k:
            g = self.A_path[i]
            i += 1
            if self.packing_label_step(g, l):
                return g
            x = self.head[g]
            if self.packing_any_unused_is_unlabelled(x):
                while not self.incident_edges_Q.empty():
                    edge = self.incident_edges_Q.front()
                    self.incident_edges_Q.pop()
                    if edge != g:
                        if self.packing_label_step(edge, g):
                            while not self.incident_edges_Q.empty():
                                self.incident_edges_Q.pop()
                            return edge

        while not self.incident_edges_Q.empty():
            self.incident_edges_Q.pop()
        return INT_MAX

    cdef int packing_fundamental_cycle_step(self, int e, int tree) noexcept:
        r"""
        Dedicated Phase 2 version of the fundamental-cycle traversal.
        Walks up tree ``tree`` from the two endpoints of edge ``e``,
        building ``A_path`` of unlabeled edges and returning the first
        joining edge found, or ``INT_MAX``.

        The disjoint-set redirects (``myset`` / ``root_set`` / ``set_seen``)
        let the traversal jump across component boundaries of a split tree
        -- this is what the connectivity version of
        :meth:`fundamental_cycle_step` lacks and why it cannot be reused
        here.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: len(GabowEdgeConnectivity(D).edge_disjoint_spanning_trees(k=4))
            4
        """
        cdef int x = self.head[e]
        cdef int y = self.tail[e]
        cdef bint stop = False
        cdef int q, side = 0, s, i, k1, k2
        cdef vector[int] left_traverse
        cdef vector[int] right_traverse

        if (not self.labeled[tree][x]) and (not self.labeled[tree][y]):
            # Reference treats this as an error; degrade gracefully.
            return INT_MAX

        if self.labeled[tree][x]:
            x = self.L_roots[tree]
            side = 1
        if self.labeled[tree][y]:
            y = self.L_roots[tree]
            side = 2

        # Disjoint-set redirect for x.
        s = self.myset[x]
        if s == -1:
            self.myset[x] = self.sets
        elif self.root_set[s][tree] != -1 and s != self.sets:
            x = self.root_set[s][tree]
            self.set_seen[s] = 1
        # Disjoint-set redirect for y.
        s = self.myset[y]
        if s == -1:
            self.myset[y] = self.sets
        elif self.root_set[s][tree] != -1 and s != self.sets:
            y = self.root_set[s][tree]
            self.set_seen[s] = 1

        if x == y:
            return INT_MAX

        while True:
            while self.depth_1[tree][x] >= self.depth_1[tree][y]:
                self.labeled[tree][x] = True
                q = self.parent_edge_id_1[tree][x]
                if q == self.UNUSED:
                    return INT_MAX
                if self.labels[q] == self.UNUSED:
                    left_traverse.push_back(q)
                    if self.packing_check_if_joining(q):
                        self.labels[q] = e
                        return q
                    x = self.parent_1[tree][x]
                    # NOTE: faithful to the reference, this redirect assigns
                    # to y (not x) inside the x-walk; kept bug-for-bug.
                    s = self.myset[x]
                    if s == -1:
                        self.myset[x] = self.sets
                    elif self.root_set[s][tree] != -1 and s != self.sets:
                        y = self.root_set[s][tree]
                        self.set_seen[s] = 1
                else:
                    stop = True
                    break
                if x == y:
                    break

            while self.depth_1[tree][y] > self.depth_1[tree][x]:
                self.labeled[tree][y] = True
                q = self.parent_edge_id_1[tree][y]
                if q == self.UNUSED:
                    return INT_MAX
                if self.labels[q] == self.UNUSED:
                    right_traverse.push_back(q)
                    if self.packing_check_if_joining(q):
                        self.labels[q] = e
                        return q
                    y = self.parent_1[tree][y]
                    s = self.myset[y]
                    if s == -1:
                        self.myset[y] = self.sets
                    elif self.root_set[s][tree] != -1 and s != self.sets:
                        y = self.root_set[s][tree]
                        self.set_seen[s] = 1
                else:
                    stop = True
                    break

            if x == y or stop:
                break

        if x == y:
            self.L_roots[tree] = x
            self.labeled[tree][x] = True

        # Assemble A_path in a consistent order for the next labeling.
        self.A_path.clear()
        k1 = <int>left_traverse.size()
        k2 = <int>right_traverse.size()
        if side == 1:
            for i in range(k1):
                self.A_path.push_back(left_traverse[i])
            for i in range(k2 - 1, -1, -1):
                self.A_path.push_back(right_traverse[i])
        elif side == 2:
            for i in range(k2):
                self.A_path.push_back(right_traverse[i])
            for i in range(k1 - 1, -1, -1):
                self.A_path.push_back(left_traverse[i])

        return self.packing_label_A_step(e, <int>self.A_path.size())

    cdef int packing_next_edge_step(self) noexcept:
        r"""
        Process queued edges and run :meth:`packing_fundamental_cycle_step`
        until the queue empties or a joining edge is found. Returns the
        joining edge id or ``INT_MAX``.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: len(GabowEdgeConnectivity(D).edge_disjoint_spanning_trees(k=4))
            4
        """
        cdef int myedge, found_joining
        cdef int i = 0

        while not self.my_Q.empty():
            myedge = self.my_Q.front()
            self.my_Q.pop()
            if self.edge_state_1[myedge] == i:
                i += 1
                if i > self.current_tree:
                    i = 0
            self.tree_flag[i] = True
            found_joining = self.packing_fundamental_cycle_step(myedge, i)
            if found_joining != INT_MAX:
                return found_joining

        return INT_MAX

    cdef void packing_join(self, int j, int tree) noexcept:
        r"""
        Assign edge ``j`` to ``tree`` and merge the two f_trees it connects
        into the root component.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: len(GabowEdgeConnectivity(D).edge_disjoint_spanning_trees(k=4))
            4
        """
        cdef int r1 = self.root[self.tail[j]]
        cdef int r2 = self.root[self.head[j]]
        self.edge_state_1[j] = tree
        self.root[r2] = self.root_vertex
        self.root[r1] = self.root_vertex
        while not self.my_Q.empty():
            self.my_Q.pop()

    cdef void packing_init_Lroots(self, int x) noexcept:
        r"""
        Initialize the ``L_i`` tree of every ``T_i`` to vertex ``x`` and
        mark ``x`` labeled.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: len(GabowEdgeConnectivity(D).edge_disjoint_spanning_trees(k=4))
            4
        """
        cdef int i
        for i in range(self.current_tree + 1):
            self.L_roots[i] = x
            self.labeled[i][x] = True

    cdef void packing_empty_the_queue(self) noexcept:
        r"""
        Empty ``my_Q``, clearing the label of every edge in it.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: len(GabowEdgeConnectivity(D).edge_disjoint_spanning_trees(k=4))
            4
        """
        cdef int f
        while not self.my_Q.empty():
            f = self.my_Q.front()
            self.my_Q.pop()
            self.labels[f] = self.UNUSED

    cdef bint search_packing_joining(self, int x, int tree) noexcept:
        r"""
        Dedicated Phase 2 joining search: try to reconnect the f_tree
        rooted at ``x`` into ``tree``, either directly via a free incoming
        edge from a different f_tree, or via a cascade of swaps found by
        :meth:`packing_next_edge_step`.

        Returns ``True`` if a joining edge was found and joined.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: len(GabowEdgeConnectivity(D).edge_disjoint_spanning_trees(k=4))
            4
        """
        cdef int y, joining_edge, e

        self.augmenting_root = self.root[x]

        # Direct reconnection: a free, non-arborescence incoming edge whose
        # source is in a different f_tree.
        for e in self.g_in[x]:
            if self.edge_state_1[e] == self.UNUSED and not self.arbor_edge[e]:
                y = self.root[self.tail[e]]
                if self.root[x] != y:
                    self.packing_empty_the_queue()
                    self.packing_join(e, tree)
                    self.tree_flag[tree] = True
                    return True
                else:
                    self.labels[e] = self.FIRSTEDGE
                    self.my_Q.push(e)

        # Otherwise run the cascade of swaps.
        self.packing_init_Lroots(x)
        joining_edge = self.packing_next_edge_step()
        if joining_edge != INT_MAX:
            self.joining_edges.push((joining_edge, self.edge_state_1[joining_edge]))
            self.packing_join(joining_edge, tree)
            return True
        return False

    cdef void rebuild_one_tree(self, int t) noexcept:
        r"""
        Recompute the undirected parent/depth structure of tree ``t`` from
        the current ``edge_state_1`` assignment, rooted at
        ``self.root_vertex``.

        Used to refresh a tree after the intersection has been mutated (by
        an augmentation, or by restoring a failed swap). Assumes the
        ``my_*`` aliases are already set to the out-arborescence direction.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: len(GabowEdgeConnectivity(D).edge_disjoint_spanning_trees(k=4))
            4
        """
        cdef int i, e_id

        for i in range(self.n):
            self.tree_edges_incident[i].clear()
        for e_id in range(self.m):
            if self.edge_state_1[e_id] == t:
                self.tree_edges_incident[self.tail[e_id]].push_back(e_id)
                self.tree_edges_incident[self.head[e_id]].push_back(e_id)

        for i in range(self.n):
            self.parent_1[t][i] = i
            self.parent_edge_id_1[t][i] = self.UNUSED
            self.depth_1[t][i] = 0
            self.seen[i] = False

        self.update_parents_dfs(t, self.root_vertex)

    cdef void update_mytree(self, int t, int e) noexcept:
        r"""
        Rebuild the f_tree partition of tree ``t`` after edge ``e`` has
        just been removed from it (``edge_state_1[e]`` already set to
        ``UNUSED``).

        Removing ``e`` splits the spanning tree ``t`` into two components:
        the one still containing the global root ``self.root_vertex`` and
        the orphaned one. This method recomputes ``parent_1[t]`` /
        ``depth_1[t]`` and sets ``root[]`` so that every vertex reachable
        from the root keeps ``root = root_vertex`` while the orphaned
        component is rooted at ``myv``. This is what lets
        :meth:`search_packing_joining` see the correct two f_trees to
        reconnect.

        Assumes the ``my_*`` aliases are set to the out-arborescence
        direction.

        INPUT:

        - ``t`` -- integer; tree index from which ``e`` was removed

        - ``e`` -- integer; the removed edge id

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: len(GabowEdgeConnectivity(D).edge_disjoint_spanning_trees(k=4))
            4
        """
        cdef int x = self.tail[e]
        cdef int y = self.head[e]
        cdef int myv
        cdef int i, e_id

        # Orphan-component root: the endpoint that was the child across e.
        if self.parent_1[t][x] == y:
            myv = x
        else:
            myv = y

        # Rebuild incident-edge lists of tree t (e already excluded).
        for i in range(self.n):
            self.tree_edges_incident[i].clear()
        for e_id in range(self.m):
            if self.edge_state_1[e_id] == t:
                self.tree_edges_incident[self.tail[e_id]].push_back(e_id)
                self.tree_edges_incident[self.head[e_id]].push_back(e_id)

        # Reset and DFS from the global root: seen[] marks the root component.
        for i in range(self.n):
            self.parent_1[t][i] = i
            self.parent_edge_id_1[t][i] = self.UNUSED
            self.depth_1[t][i] = 0
            self.seen[i] = False
        self.update_parents_dfs(t, self.root_vertex)

        # Partition vertices into the root component and the orphan component.
        for i in range(self.n):
            if self.seen[i]:
                self.root[i] = self.root_vertex
            else:
                self.root[i] = myv

        # DFS over the orphan component to set its parent/depth.
        self.update_parents_dfs(t, myv)

    cdef void refresh_all_intersection_trees(self) noexcept:
        r"""
        Recompute the parent/depth structure of every tree
        ``0 .. packing_ntrees - 1`` in the current intersection, after an
        augmentation has rearranged edges between trees. Assumes the
        ``my_*`` aliases are set to the out-arborescence direction.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: len(GabowEdgeConnectivity(D).edge_disjoint_spanning_trees(k=4))
            4
        """
        cdef int t
        for t in range(self.packing_ntrees):
            self.rebuild_one_tree(t)

    cdef bint search_step(self, int e) except -1:
        r"""
        Attempt to free edge ``e`` (currently in some tree ``t`` of the
        intersection) for use in the arborescence currently being
        extracted by :meth:`find_tree`, by swapping it out of ``t`` via
        the dedicated packing search.

        Protocol:

        1. Set up the ``my_*`` aliases for the out-arborescence direction
           and ``current_tree`` to the full intersection size so the
           cascade considers every tree.
        2. Temporarily remove ``e`` from tree ``t`` (mark it
           ``TEMP_REMOVED``, not ``UNUSED``, so the search does not
           re-select the very edge being freed), and rebuild tree ``t``'s
           f_tree partition with :meth:`update_mytree`.
        3. Call :meth:`search_packing_joining` with ``y = head(e)`` to
           reconnect the orphaned component.
        4. On success: run :meth:`augmentation_algorithm` to finalize the
           edge transfers, set ``e`` to ``UNUSED`` (now genuinely free),
           refresh all trees, and clear the scratch state.
        5. On failure: restore ``e`` to tree ``t`` and rebuild tree
           ``t``'s structure.

        INPUT:

        - ``e`` -- integer; candidate edge currently in some tree of the
          k-intersection

        OUTPUT:

        ``True`` on a successful swap (``e`` is now ``UNUSED`` in the
        intersection), ``False`` otherwise (state restored).

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: len(GabowEdgeConnectivity(D).edge_disjoint_spanning_trees(k=4))
            4
        """
        cdef int t = self.edge_state_1[e]
        cdef int y = self.head[e]

        # Set up my_* aliases for the out-arborescence direction so the
        # reused augmentation_algorithm / trace_back operate on the right
        # arrays (tail/head, parent_1/depth_1, edge_state_1).
        self.my_g = self.g_in
        self.my_g_reversed = self.g_out
        self.my_from = self.tail
        self.my_to = self.head
        self.my_parent = self.parent_1
        self.my_depth = self.depth_1
        self.my_parent_edge_id = self.parent_edge_id_1
        self.my_edge_state = self.edge_state_1
        # The cascade ranges over trees 0 .. current_tree, so set it to the
        # full intersection size for the search to consider every tree.
        self.current_tree = self.packing_ntrees - 1

        # Temporarily remove e from tree t, then rebuild tree t's f_tree
        # partition (root component vs orphaned component) so the dedicated
        # packing search sees the correct two f_trees to reconnect.
        # Mark e with TEMP_REMOVED (not UNUSED) so the search does NOT
        # re-select the very edge we are trying to free as a reconnection
        # edge.
        self.edge_state_1[e] = self.TEMP_REMOVED
        self.update_mytree(t, e)

        cdef bint found = self.search_packing_joining(y, t)

        if found:
            # The swap succeeded: e is now genuinely free for the
            # arborescence. Finalize edge transfers, refresh all trees (the
            # augmentation rearranged edges across trees), and clean up.
            self.augmentation_algorithm()
            self.edge_state_1[e] = self.UNUSED
            self.refresh_all_intersection_trees()
            self.clear_augmentation_aux_state()
            return True

        # Failure: restore e to tree t and rebuild tree t's structure.
        self.edge_state_1[e] = t
        self.rebuild_one_tree(t)
        self.clear_augmentation_aux_state()
        return False

    cdef int choose_edge_step(self, int v) noexcept:
        r"""
        Find the next outgoing edge from ``v`` that can extend the
        arborescence currently being grown by :meth:`find_tree`.

        Scans the outgoing edges of ``v`` for one whose head is not yet
        in the reachable set ``A`` and that has not been marked in this
        ``find_tree`` call. If such an edge is "free" (not in any tree
        of the k-intersection), it is added directly to the
        arborescence (``arbor_edge[e] = True``) and ``A`` is grown. If
        the edge belongs to a tree, it is returned as a candidate for
        :meth:`search_step`.

        INPUT:

        - ``v`` -- integer; internal vertex index to extend from

        OUTPUT:

        ``-2`` if a free edge was added directly, ``-1`` if no usable
        edge remains, otherwise the candidate edge id (``>= 0``).

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: len(GabowEdgeConnectivity(D).edge_disjoint_spanning_trees(k=4))
            4
        """
        cdef int j, e, x

        for j in range(<int>self.g_out[v].size()):
            e = self.g_out[v][j]
            x = self.head[e]
            if self.in_A[x]:
                continue
            if self.marked_edge[e]:
                continue
            if self.extracted_edge[e]:
                # Edge already taken by a previously extracted arborescence.
                continue
            if self.edge_state_1[e] == self.UNUSED:
                # Free edge: add directly to the arborescence.
                self.arbor_edge[e] = True
                self.in_A[x] = True
                self.A_order[self.A_size] = x
                self.A_size += 1
                return -2
            # Edge belongs to a tree of the intersection: candidate for swap.
            return e

        return -1

    cdef int period_step(self) noexcept:
        r"""
        Pick the next non-exhausted vertex from ``A`` to extend in
        :meth:`find_tree`.

        Iterates the ordered list ``A_order`` and returns the first
        vertex not marked exhausted (``in_X[v] == False``). Returns
        ``-1`` when the arborescence is complete (``A_size == n``) or
        when every vertex in ``A`` is exhausted.

        OUTPUT:

        Internal vertex index, or ``-1`` when no candidate remains.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: len(GabowEdgeConnectivity(D).edge_disjoint_spanning_trees(k=4))
            4
        """
        cdef int i, v

        # Fresh disjoint-set for the next vertex's exploration.
        self.packing_init_set()

        if self.A_size >= self.n:
            return -1

        for i in range(self.A_size):
            v = self.A_order[i]
            if not self.in_X[v]:
                return v

        return -1

    cdef bint find_tree(self) except -1:
        r"""
        Extract one spanning out-arborescence rooted at
        ``self.root_vertex`` from the current k-intersection, marking
        the chosen edges in ``self.arbor_edge``.

        Grows a reachable set ``A`` starting from the root. For each
        vertex ``v`` picked from ``A`` (via :meth:`period_step`),
        :meth:`choose_edge_step` finds an outgoing edge to a vertex not
        yet in ``A``: free edges are added directly, edges already
        assigned to a tree of the intersection are routed through
        :meth:`search_step` (the swap). When a vertex has no usable
        outgoing edge, it is marked exhausted and the next vertex from
        ``A`` is processed.

        An edge already used by a tree of the intersection is not simply
        skipped: :meth:`search_step` performs a swap (an augmenting-path
        search) that frees a usable edge while preserving the validity of the
        other trees. This is what makes the greedy growth of ``A`` always
        succeed on a genuine ``k``-intersection.

        OUTPUT:

        ``True`` if a spanning out-arborescence was successfully
        extracted (every vertex reached); ``False`` otherwise.

        EXAMPLES:

        This method drives the Phase 2 extraction and is exercised through
        the public :meth:`edge_disjoint_spanning_trees`::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: trees = GabowEdgeConnectivity(D).edge_disjoint_spanning_trees(k=4)
            sage: len(trees)
            4
        """
        cdef int j, v, found, e, x

        # Reset per-call state.
        for j in range(self.n):
            self.in_A[j] = False
            self.in_X[j] = False
            self.A_order[j] = 0
        for j in range(self.m):
            self.arbor_edge[j] = False
            self.marked_edge[j] = False

        # Seed with the root.
        self.in_A[self.root_vertex] = True
        self.A_order[0] = self.root_vertex
        self.A_size = 1
        v = self.root_vertex

        # Fresh disjoint-set for this extraction.
        self.packing_init_set()

        # Grow the arborescence one vertex at a time.
        while v != -1:
            sig_check()
            found = self.choose_edge_step(v)

            if found == -2:
                # ``choose_edge_step`` already added a free edge.
                v = self.period_step()
            elif found == -1:
                # ``v`` has no usable outgoing edge in this call.
                self.in_X[v] = True
                v = self.period_step()
            else:
                # ``found`` is a candidate edge in the intersection.
                e = found
                if self.search_step(e):
                    # Swap succeeded; take ``e`` and grow ``A``.
                    self.arbor_edge[e] = True
                    x = self.head[e]
                    self.in_A[x] = True
                    self.A_order[self.A_size] = x
                    self.A_size += 1
                    v = self.period_step()
                else:
                    # Swap failed; record the explored sets and mark ``e``,
                    # then retry with the same ``v``.
                    self.packing_merge_sets()
                    self.marked_edge[e] = True

        return self.A_size == self.n

    def edge_disjoint_spanning_trees(self, k=None, root=None):
        r"""
        Return ``k`` edge-disjoint spanning arborescences rooted at ``root``.

        The returned arborescences are out-arborescences (rooted spanning
        trees with all edges directed away from the root).
        Following [Gabow1995]_ and Edmonds' branching theorem, the maximum
        number of edge-disjoint spanning out-arborescences rooted at a
        vertex ``r`` equals the minimum ``r``-cut of the digraph. When the
        digraph is strongly connected, this is at least the edge
        connectivity ``ec`` of the digraph.

        INPUT:

        - ``k`` -- integer or ``None`` (default: ``None``); the number of
          edge-disjoint spanning arborescences to return. If ``None``, the
          method returns ``self.ec`` arborescences (the value computed by
          :meth:`edge_connectivity`).

        - ``root`` -- vertex of the input digraph or ``None``
          (default: ``None``); the common root of the returned
          arborescences. If ``None``, defaults to the first vertex of the
          digraph (the same convention used by the MILP backend in
          :meth:`~sage.graphs.generic_graph.GenericGraph.edge_disjoint_spanning_trees`).

        OUTPUT:

        A list of :class:`~sage.graphs.digraph.DiGraph`, each on the same
        vertex set as the input. Every returned arborescence is rooted at
        ``root``: the root has indegree 0, every other vertex has
        indegree 1, and every vertex is reachable from ``root``. The
        returned arborescences are pairwise edge-disjoint.

        ALGORITHM:

        This uses the matroid-intersection approach of [Gabow1995]_, in two
        phases. Phase 1 builds a ``k``-intersection whose union is a
        ``k``-in-regular subgraph rooted at ``root`` (every non-root vertex
        has in-degree ``k``). Phase 2 decomposes that union into ``k``
        spanning out-arborescences one at a time (:meth:`find_tree`): it grows
        the set of vertices reachable from ``root``, and whenever the next edge
        it needs is already used by another tree, an augmenting-path *swap*
        (:meth:`search_step`) frees a usable edge while keeping the other trees
        valid.

        The extraction scheme follows the practical implementation studied in
        [GKMN2022]_. Phase 1 runs in time `O(km\log(n^2/m))` [Gabow1995]_.
        Phase 2 performs `k` extraction rounds, each rebuilding an
        intersection on the remaining edge pool before carving out one
        arborescence, and its swaps dominate the running time in practice:
        the observed total grows roughly like `k^2 m` (about `n^4` on
        complete digraphs, where `k = n - 1`).

        EXAMPLES:

        On a complete digraph, the algorithm returns ``n-1`` arborescences::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: trees = GabowEdgeConnectivity(D).edge_disjoint_spanning_trees(k=4)
            sage: len(trees)
            4
            sage: all(T.num_edges() == 4 for T in trees)
            True

        By default (``k=None``) the method returns ``ec`` arborescences,
        where ``ec`` is the edge connectivity::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: G = GabowEdgeConnectivity(D)
            sage: G.edge_connectivity()
            4
            sage: len(G.edge_disjoint_spanning_trees())
            4

        The arborescences can be rooted at any vertex; the root then has
        indegree 0, every other vertex has indegree 1, and every vertex is
        reachable from the root::

            sage: D = digraphs.Complete(5)
            sage: trees = GabowEdgeConnectivity(D).edge_disjoint_spanning_trees(k=4, root=2)
            sage: all(T.in_degree(2) == 0 for T in trees)
            True
            sage: all(T.in_degree(v) == 1 for T in trees for v in D if v != 2)
            True
            sage: all(set(T.depth_first_search(2)) == set(D) for T in trees)
            True

        The returned arborescences are pairwise edge-disjoint::

            sage: D = digraphs.Complete(5)
            sage: trees = GabowEdgeConnectivity(D).edge_disjoint_spanning_trees(k=4)
            sage: all_edges = sum((T.edges(labels=False, sort=False) for T in trees), [])
            sage: len(all_edges) == len(set(all_edges))
            True

        On a directed cycle the edge connectivity is 1, so the packing has a
        single arborescence: the directed path from the root::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Circuit(4)
            sage: trees = GabowEdgeConnectivity(D).edge_disjoint_spanning_trees(root=0)
            sage: len(trees)
            1
            sage: sorted(trees[0].edges(labels=False, sort=False))
            [(0, 1), (1, 2), (2, 3)]

        Trivial case::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: GabowEdgeConnectivity(D).edge_disjoint_spanning_trees(k=0)
            []

        Asking for more arborescences than the digraph can support raises
        :class:`~sage.categories.sets_cat.EmptySetError`::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: from sage.categories.sets_cat import EmptySetError
            sage: D = digraphs.Complete(5)
            sage: try:
            ....:     _ = GabowEdgeConnectivity(D).edge_disjoint_spanning_trees(k=5)
            ....: except EmptySetError:
            ....:     print("EmptySetError")
            EmptySetError

        An unknown ``root`` raises a :class:`ValueError`::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: GabowEdgeConnectivity(D).edge_disjoint_spanning_trees(k=2, root=99)
            Traceback (most recent call last):
            ...
            ValueError: vertex 99 is not a vertex of the graph

        TESTS:

        Randomized validation: on random strongly connected digraphs, every
        packing has the expected size and each returned graph is a valid
        out-arborescence (root in-degree 0, every other vertex in-degree 1,
        all reachable from the root), the arborescences being pairwise
        edge-disjoint::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: def is_valid_packing(D, root, trees):
            ....:     V = set(D)
            ....:     if any(T.in_degree(root) != 0 for T in trees):
            ....:         return False
            ....:     if any(T.in_degree(v) != 1 for T in trees for v in V if v != root):
            ....:         return False
            ....:     if any(set(T.depth_first_search(root)) != V for T in trees):
            ....:         return False
            ....:     edges = [e for T in trees for e in T.edge_iterator(labels=False)]
            ....:     return len(edges) == len(set(edges))
            sage: for _ in range(20):
            ....:     n = randint(3, 7)
            ....:     D = digraphs.RandomDirectedGNP(n, 0.5)
            ....:     if not D.is_strongly_connected():
            ....:         continue
            ....:     ec = GabowEdgeConnectivity(D).edge_connectivity()
            ....:     if ec == 0:
            ....:         continue
            ....:     for r in D:
            ....:         trees = GabowEdgeConnectivity(D).edge_disjoint_spanning_trees(k=ec, root=r)
            ....:         assert len(trees) == ec
            ....:         assert is_valid_packing(D, r, trees)
        """
        from sage.graphs.digraph import DiGraph
        from sage.categories.sets_cat import EmptySetError

        # Default k to the edge connectivity
        if k is None:
            if not self.ec_checked:
                self.compute_edge_connectivity()
            k = self.ec

        # Trivial / empty packing
        if k <= 0:
            return []

        # Empty or non-strongly-connected digraphs are not handled by the
        # Gabow machinery: __init__ returns early without building the
        # compact data structures, leaving ``int_to_vertex`` unset. No
        # spanning arborescence packing of size k >= 1 exists in that case
        # via this backend (the public API may fall back to MILP).
        if self.int_to_vertex is None:
            raise EmptySetError(
                "this digraph does not contain {} edge-disjoint spanning "
                "arborescences".format(k))

        # Resolve the root vertex (default: first vertex of the graph)
        cdef int root_idx = 0
        if root is not None:
            if root not in self.int_to_vertex:
                raise ValueError("vertex {!r} is not a vertex of the graph".format(root))
            root_idx = self.int_to_vertex.index(root)

        # Run the packing
        if not self.compute_arborescence_packing(k, root_idx):
            raise EmptySetError(
                "this digraph does not contain {} edge-disjoint spanning "
                "arborescences rooted at {!r}".format(k, self.int_to_vertex[root_idx]))

        # Build DiGraphs from the saved arborescence edge sets
        cdef int i
        result = []
        for i in range(k):
            edges = ((self.int_to_vertex[self.tail[e_id]], self.int_to_vertex[self.head[e_id]])
                     for e_id in self.arborescence_F[i])
            result.append(DiGraph([self.int_to_vertex, edges], format='vertices_and_edges'))

        return result
