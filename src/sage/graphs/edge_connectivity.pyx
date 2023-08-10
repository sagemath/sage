# cython: binding=True
# distutils: language = c++
r"""
Edge connectivity

This module implements methods for computing the edge-connectivity of graphs and
digraphs. It also implements methods to extract `k` edge-disjoint spanning trees
from a `2k` edge-connected graph or a `k` edge-connected digraph.

.. TODO::

    - Implement the tree-packing algorithms proposed in [Gabow1995]_ and
      [BHKP2008]_
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

    :trac:`32169`::

        sage: dig6_string = r'[E_S?_hKIH@eos[BSg???Q@FShGC?hTHUGM?IPug?'
        sage: dig6_string += r'JOEYCdOzdkQGo@ADA@AAg?GAQW?'
        sage: dig6_string += r'[aIaSwHYcD@qQb@Dd?\hJTI@OHlJ_?C_OEIKoeCR@_BC?Q?'
        sage: dig6_string += r'?YBFosqITEA?IvCU_'
        sage: D = DiGraph(dig6_string)
        sage: GabowEdgeConnectivity(D).edge_connectivity()
        5
        sage: GabowEdgeConnectivity(D).edge_disjoint_spanning_trees()
        Traceback (most recent call last):
        ...
        NotImplementedError: this method has not been implemented yet

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

    # The graph is stored as lists of incident edges
    cdef readonly GenericGraph_pyx G  # the original graph
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

        # Set some constants
        self.UNUSED = INT_MAX
        self.FIRSTEDGE = INT_MAX - 1

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

    cdef void increase_memory_for_new_tree(self, int tree):
        """
        Allocate data structure for the new tree/forest.

        This method also initializes data structures for this tree index. Data
        structures for a given tree index are allocatated only once.

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

    cdef void compute_dfs_tree(self):
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

    cdef void find_dfs_tree(self, int r):
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

    cdef void find_dfs_tree_rec(self, int u, int r):
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

    cdef int choose_root(self):
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

    cdef void join(self, int e_id):
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

    cdef bint is_joining_edge(self, int e_id):
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

    cdef int label_A_path(self, int e_id):
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

    cdef bint label_step(self, int e_id, int e_label):
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

    cdef bint any_unused_is_unlabeled(self, int x):
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

    cdef void augmentation_algorithm(self):
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

    cdef void trace_back(self, int e_id, int e_state):
        """
        Trace the path of a joining edge and transfer the edges to the
        appropriate tree Ti.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: GabowEdgeConnectivity(D).edge_connectivity()
            4
        """
        # Target x and source y of joining edge e_id
        cdef int x = self.my_to[e_id]
        cdef int y = self.my_from[e_id]
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

    cdef void update_parents_depths(self, int tree):
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

    cdef void update_parents_dfs(self, int tree, int x):
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

    cdef void save_current_k_intersection(self):
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

    def edge_disjoint_spanning_trees(self):
        r"""
        Iterator over the edge disjoint spanning trees.

        EXAMPLES::

            sage: from sage.graphs.edge_connectivity import GabowEdgeConnectivity
            sage: D = digraphs.Complete(5)
            sage: GabowEdgeConnectivity(D).edge_disjoint_spanning_trees()
            Traceback (most recent call last):
            ...
            NotImplementedError: this method has not been implemented yet
        """
        raise NotImplementedError('this method has not been implemented yet')
