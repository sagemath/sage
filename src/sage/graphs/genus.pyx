"""
Genus

This file contains a moderately-optimized implementation to compute the
genus of simple connected graph.  It runs about a thousand times faster
than the previous version in Sage, not including asymptotic improvements.

The algorithm works by enumerating combinatorial embeddings of a graph,
and computing the genus of these via the Euler characteristic.  We view
a combinatorial embedding of a graph as a pair of permutations `v,e`
which act on a set `B` of `2|E(G)|` "darts".  The permutation `e` is an
involution, and its orbits correspond to edges in the graph.  Similarly,
The orbits of `v` correspond to the vertices of the graph, and those of
`f = ve` correspond to faces of the embedded graph.

The requirement that the group `<v,e>` acts transitively on `B` is
equivalent to the graph being connected.  We can compute the genus of a
graph by

    `2 - 2g = V - E + F`

where `E`, `V`, and `F` denote the number of orbits of `e`, `v`, and
`f` respectively.

We make several optimizations to the naive algorithm, which are
described throughout the file.
"""

# ****************************************************************************
#         Copyright (C) 2010 Tom Boothby <tomas.boothby@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from libc.string cimport memcpy
from memory_allocator cimport MemoryAllocator
from cysignals.signals cimport sig_on, sig_off

# cimport sage.combinat.permutation_cython

from sage.combinat.permutation_cython cimport next_swap, reset_swap

from sage.graphs.base.dense_graph cimport DenseGraph
from sage.graphs.graph import Graph


cdef inline int edge_map(int i) noexcept:
    """
    We might as well make the edge map nice, since the vertex map is so
    slippery.  This is the fastest way I could find to establish the
    correspondence `i <-> i + 1` if `i` is even.
    """
    return i - 2 * (i & 1) + 1


cdef class simple_connected_genus_backtracker:
    r"""
    A class which computes the genus of a DenseGraph through an extremely slow
    but relatively optimized algorithm.  This is "only" exponential for graphs
    of bounded degree, and feels pretty snappy for 3-regular graphs. The
    generic runtime is

      `|V(G)| \prod_{v \in V(G)} (deg(v)-1)!`

    which is `2^{|V(G)|}` for 3-regular graphs, and can achieve `n(n-1)!^{n}`
    for the complete graph on `n` vertices.  We can handily compute the genus of
    `K_6` in milliseconds on modern hardware, but `K_7` may take a few days.
    Don't bother with `K_8`, or any graph with more than one vertex of degree 10
    or worse, unless you can find an a priori lower bound on the genus and
    expect the graph to have that genus.

    .. WARNING::

        THIS MAY SEGFAULT OR HANG ON:
            * DISCONNECTED GRAPHS
            * DIRECTED GRAPHS
            * LOOPED GRAPHS
            * MULTIGRAPHS

    EXAMPLES::

        sage: import sage.graphs.genus
        sage: G = graphs.CompleteGraph(6)
        sage: G = Graph(G, sparse=False)
        sage: bt = sage.graphs.genus.simple_connected_genus_backtracker(G._backend.c_graph()[0])
        sage: bt.genus() #long time
        1
        sage: bt.genus(cutoff=1)
        1
        sage: G = graphs.PetersenGraph()
        sage: G = Graph(G, sparse=False)
        sage: bt = sage.graphs.genus.simple_connected_genus_backtracker(G._backend.c_graph()[0])
        sage: bt.genus()
        1
        sage: G = graphs.FlowerSnark()
        sage: G = Graph(G, sparse=False)
        sage: bt = sage.graphs.genus.simple_connected_genus_backtracker(G._backend.c_graph()[0])
        sage: bt.genus()
        2
    """
    cdef MemoryAllocator mem
    cdef int **vertex_darts
    cdef int *face_map
    cdef int *degree
    cdef int *visited
    cdef int *face_freeze
    cdef int **swappers
    cdef int num_darts, num_verts, num_cycles, record_genus

    def __init__(self, DenseGraph G):
        """
        Initialize the genus_backtracker object.

        TESTS::

            sage: import sage.graphs.genus
            sage: G = Graph(sparse=False)  # indirect doctest
            sage: gb = sage.graphs.genus.simple_connected_genus_backtracker(G._backend.c_graph()[0])
            sage: G = Graph(graphs.CompleteGraph(4), sparse=False)
            sage: gb = sage.graphs.genus.simple_connected_genus_backtracker(G._backend.c_graph()[0])
            sage: gb.genus()
            0
        """
        self.num_darts = G.num_arcs
        self.num_verts = G.num_verts

        if self.num_verts <= 1:
            return

        # Allocate arrays
        self.mem = MemoryAllocator()
        self.degree = <int *> self.mem.malloc(self.num_verts * sizeof(int))
        self.face_map = <int *> self.mem.malloc(self.num_darts * sizeof(int))
        self.visited = <int *> self.mem.malloc(self.num_darts * sizeof(int))
        self.face_freeze = <int *> self.mem.malloc(self.num_darts * sizeof(int))
        self.vertex_darts = <int **>self.mem.malloc(self.num_verts * sizeof(int *))
        self.swappers = <int **>self.mem.malloc(self.num_verts * sizeof(int *))
        cdef int *w = <int *> self.mem.malloc((self.num_verts + self.num_darts) * sizeof(int))
        cdef int *s = <int *> self.mem.malloc(2 * (self.num_darts - self.num_verts) * sizeof(int))

        cdef int i, j, dv, u, v

        for v in range(self.num_verts):
            if not G.has_vertex(v):
                raise ValueError("Please relabel G so vertices are 0, ..., n-1")

            dv = G.in_degrees[v]
            self.degree[v] = 0
            self.vertex_darts[v] = w
            w += dv + 1

            self.swappers[v] = s
            s += 2 * (dv - 1)

        i = 0
        for v in range(self.num_verts):
            dv = self.degree[v]

            # we use self.face_map as a temporary int array to hold
            # neighbors of v since it will be overwritten shortly.
            G.in_neighbors_unsafe(v, self.face_map, G.in_degrees[v])
            for j in range(G.in_degrees[v]):
                u = self.face_map[j]
                if u < v:
                    # edge hasn't been seen yet
                    self.vertex_darts[u][self.degree[u]] = i
                    self.vertex_darts[v][dv] = i + 1
                    self.degree[u] += 1
                    dv += 1
                    i += 2

            self.degree[v] = dv

        for v in range(self.num_verts):
            dv = self.degree[v]
            w = self.vertex_darts[v]
            w[dv] = w[0]
            for i in range(dv):
                u = w[i]
                self.face_map[edge_map(u)] = w[i + 1]

        self.freeze_face()

#   good for debugging
#    def dump(self):
#        cdef int v, j
#        print("vertex darts:", end="")
#        for v in range(self.num_verts):
#            print('(', end="")
#            for j in range(self.degree[v] + 1):
#                print(self.vertex_darts[v][j], end="")
#            print(')', end="")
#        print("\n")

#        print("face map: [", end="")
#        for v in range(self.num_darts):
#            print(self.face_map[v], end="")
#        print(']')

    cdef inline void freeze_face(self) noexcept:
        """
        Quickly store the current face_map so we can recover
        the embedding it corresponds to later.
        """

        memcpy(self.face_freeze, self.face_map, self.num_darts * sizeof(int))

    def get_embedding(self):
        """
        Return an embedding for the graph.

        If ``min_genus_backtrack`` has been called with ``record_embedding =
        True``, then this will return the first minimal embedding that we found.
        Otherwise, this returns the first embedding considered.

        EXAMPLES::

            sage: import sage.graphs.genus
            sage: G = Graph(graphs.CompleteGraph(5), sparse=False)
            sage: gb = sage.graphs.genus.simple_connected_genus_backtracker(G._backend.c_graph()[0])
            sage: gb.genus(record_embedding=True)
            1
            sage: gb.get_embedding()
            {0: [1, 2, 3, 4], 1: [0, 2, 3, 4], 2: [0, 1, 4, 3], 3: [0, 2, 1, 4], 4: [0, 3, 1, 2]}
            sage: G = Graph(sparse=False)
            sage: G.add_edge(0,1)
            sage: gb = sage.graphs.genus.simple_connected_genus_backtracker(G._backend.c_graph()[0])
            sage: gb.get_embedding()
            {0: [1], 1: [0]}
            sage: G = Graph(sparse=False)
            sage: gb = sage.graphs.genus.simple_connected_genus_backtracker(G._backend.c_graph()[0])
            sage: gb.get_embedding()
            {}
        """
        if not self.num_verts:
            return {}
        elif self.num_verts == 1:
            return {0: []}

        cdef int i, j, v
        cdef int *w
        cdef int *face_map = self.face_freeze
        cdef list darts_to_verts = [0 for i in range(self.num_darts)]
        cdef dict embedding = {}
        cdef list orbit_v

        for v in range(self.num_verts):
            w = self.vertex_darts[v]
            for i in range(self.degree[v]):
                darts_to_verts[w[i]] = v

        for v in range(self.num_verts):
            w = self.vertex_darts[v]
            i = w[0]
            orbit_v = [darts_to_verts[edge_map(i)]]

            j = face_map[edge_map(i)]
            while j != i:
                orbit_v.append(darts_to_verts[edge_map(j)])
                j = face_map[edge_map(j)]

            embedding[v] = orbit_v

        return embedding

    cdef int run_cycle(self, int i) noexcept:
        r"""
        Mark off the orbit of `i` under face_map.

        If `i` has been visited recently, bail immediately and return 0.
        Otherwise, visit each vertex in the orbit of `i` and set
        `self.visited[j] = k`, where `j` is the `k`-th element in the orbit
        of `i`.  Then, return 1.

        In this manner, we are able to quickly check if a particular element has
        been visited recently.  Moreover, we are able to distinguish what order
        three elements of a single orbit come in. This is important for
        `self.flip()`, and discussed in more detail there.
        """
        if self.visited[i]:
            return 0

        cdef int counter = 1
        cdef int j = self.face_map[i]

        self.visited[i] = 1
        counter += 1
        while i != j:
            self.visited[j] = 1 + counter
            counter += 1
            j = self.face_map[j]
        return 1

    cdef void flip(self, int v, int i) noexcept:
        r"""
        This is where the real work happens. Once cycles have been counted for
        the initial face_map, we make small local changes, and look at their
        effect on the number of cycles.

        Consider a vertex whose embedding is given by the cycle

        `self.vertex_darts[v] = [..., v0, v1, v2, ... ]`.

        which implies that the vertex map has the cycle

        `... -> v0 -> v1 -> v2 -> ... `

        and say we'd like to exchange a1 and a2.  Then, we'll change the vertex
        map to

        `... -> v0 -> v2 -> v1 -> ...`

        and when this happens, we change the face map orbit of `e0 = e(av)`,
        `e1 = e(v1)`, and `e2 = e(v2)`, where `e` denotes the edge map.

        In fact, the only orbits that can change are those of `e0`, `e1`, and
        `e2`. Thus, to determine the effect of the flip on the cycle structure,
        we need only consider these orbits.

        We find that the set of possibilities for a flip to change the number of
        orbits among these three elements is very small.  In particular,

        * If the three elements belong to distinct orbits, a flip joins them
          into a single orbit.

        * If the three elements are among exactly two orbits, a flip does not
          change that fact (though it does break the paired elements and make a
          new pair, all we care about is the number of cycles)

        * If all three elements are in the same orbit, a flip either disconnects
          them into three distinct orbits, or maintains status quo.

        To differentiate these situations, we need only to look at the order of
        `v0`, `v1`, and `v2` under the orbit. If `e0 -> ... -> e2 -> ... -> e1`
        before the flip, the cycle breaks into three.  Otherwise, the number of
        cycles stays the same.
        """
        cdef int *w = self.vertex_darts[v]
        cdef int *face_map = self.face_map

        cdef int v0, v1, v2, e0, e1, e2, f0, f1, f2, j, k

        v0 = w[i - 1]
        v1 = w[i]
        v2 = w[i + 1]

        e0 = edge_map(v0)
        e1 = edge_map(v1)
        e2 = edge_map(v2)

        f0 = face_map[e0]
        f1 = face_map[e1]
        f2 = face_map[e2]

        face_map[e0] = -1
        face_map[e1] = -2
        face_map[e2] = -3

        j = face_map[f0]
        while j >= 0:
            j = face_map[j]
        if j != -2:
            k = face_map[f1]
            while k >= 0:
                k = face_map[k]

            # Magic function follows. There are only four possibilities for j
            # and k since j != -2. We use magic to avoid branching.
            #  j | k  | MF(j,k)
            # ---+----+--------
            # -1 | -2 |   -2
            # -1 | -3 |    0
            # -3 | -1 |    2
            # -3 | -2 |    0

            self.num_cycles += (2 * k + 1 - j) % 4

        face_map[e0] = v2
        face_map[e1] = f2
        face_map[e2] = v1

        w[i] = v2
        w[i + 1] = v1

    cdef int count_cycles(self) noexcept:
        """
        Count all cycles.
        """
        cdef int i
        self.num_cycles = 0

        for i in range(self.num_darts):
            self.visited[i] = 0

        for i in range(self.num_darts):
            self.num_cycles += self.run_cycle(i)

    def genus(self, int style=1, int cutoff=0, bint record_embedding=False):
        r"""
        Compute the minimal or maximal genus of ``self``'s graph.

        Note, this is a remarkably naive algorithm for a very difficult problem.
        Most interesting cases will take millennia to finish, with the exception
        of graphs with max degree 3.

        INPUT:

        - ``style`` -- integer (default: `1`); find minimum genus if 1,
          maximum genus if 2

        - ``cutoff`` -- integer (default: `0`); stop searching if search style
          is 1 and ``genus`` `\leq` ``cutoff``, or if style is 2 and ``genus``
          `\geq` ``cutoff``.  This is useful where the genus of the graph has a
          known bound.

        - ``record_embedding`` -- boolean (default: ``False``); whether or not
          to remember the best embedding seen. This embedding can be retrieved
          with ``self.get_embedding()``.

        OUTPUT: the minimal or maximal genus for ``self``'s graph

        EXAMPLES::

            sage: import sage.graphs.genus
            sage: G = Graph(graphs.CompleteGraph(5), sparse=False)
            sage: gb = sage.graphs.genus.simple_connected_genus_backtracker(G._backend.c_graph()[0])
            sage: gb.genus(cutoff=2, record_embedding=True)
            2
            sage: E = gb.get_embedding()
            sage: gb.genus(record_embedding=False)
            1
            sage: gb.get_embedding() == E
            True
            sage: gb.genus(style=2, cutoff=5)
            3
            sage: G = Graph(sparse=False)
            sage: gb = sage.graphs.genus.simple_connected_genus_backtracker(G._backend.c_graph()[0])
            sage: gb.genus()
            0
        """
        cdef int g

        # in the original genus implementation, this case resulted in infinite
        # recursion.  oops.  Let's skip that.
        if self.num_verts <= 0:
            return 0
        sig_on()
        if style == 1:
            g = self.genus_backtrack(cutoff, record_embedding, &min_genus_check)
        elif style == 2:
            g = self.genus_backtrack(cutoff, record_embedding, &max_genus_check)
        sig_off()
        return g

    cdef void reset_swap(self, int v) noexcept:
        """
        Reset the swapper associated with vertex ``v``.
        """
        cdef int d = self.degree[v] - 1
        reset_swap(d, self.swappers[v], self.swappers[v] + d)

    cdef int next_swap(self, int v) noexcept:
        """
        Compute and return the next swap associated with the vertex ``v``.
        """
        cdef int d = self.degree[v] - 1
        return next_swap(d, self.swappers[v], self.swappers[v] + d)

    cdef int genus_backtrack(self,
                             int cutoff,
                             bint record_embedding,
                             (int (*)(simple_connected_genus_backtracker, int, bint, int) noexcept) check_embedding) noexcept:
        """
        Here's the main backtracking routine.

        We iterate over all embeddings of ``self``'s graph by considering all
        cyclic orderings of ``self.vertex_darts``.  We use the Steinhaus-
        Johnson-Trotter algorithm to enumerate these by walking over a poly-ary
        Gray code, and each time the Gray code would flip a bit, we apply the
        next adjacent transposition from S-J-T at that vertex.

        We start by counting the number of cycles for our initial embedding.
        From that point forward, we compute the amount that each flip changes
        the number of cycles.
        """
        cdef int next_swap, vertex

        for vertex in range(self.num_verts):
            self.reset_swap(vertex)

        vertex = self.num_verts - 1

        self.count_cycles()

        if check_embedding(self, cutoff, record_embedding, 1):
            return self.record_genus

        next_swap = self.next_swap(vertex)
        while True:
            while next_swap == -1:
                self.reset_swap(vertex)
                vertex -= 1
                if vertex < 0:
                    return self.record_genus
                next_swap = self.next_swap(vertex)
            self.flip(vertex, next_swap + 1)

            if check_embedding(self, cutoff, record_embedding, 0):
                return self.record_genus

            vertex = self.num_verts-1
            next_swap = self.next_swap(vertex)

cdef int min_genus_check(simple_connected_genus_backtracker self,
                         int cutoff,
                         bint record_embedding,
                         int initial) noexcept:
    """
    Search for the minimal genus.

    If we find a genus <= cutoff, return 1 to quit entirely.
    If we find a better genus than previously recorded, keep track of that, and
    if record_embedding is set, record the face map with self.freeze_face()
    """
    cdef int g = 1 - (self.num_verts - self.num_darts / 2 + self.num_cycles) / 2
    if g < self.record_genus or initial == 1:
        self.record_genus = g
        if record_embedding:
            self.freeze_face()
        if g <= cutoff:
            return 1
    return 0

cdef int max_genus_check(simple_connected_genus_backtracker self,
                         int cutoff,
                         bint record_embedding,
                         int initial) noexcept:
    """
    Same as min_genus_check, but search for a maximum.
    """
    cdef int g = 1 - (self.num_verts - self.num_darts / 2 + self.num_cycles) / 2
    if g > self.record_genus or initial == 1:
        self.record_genus = g
        if record_embedding:
            self.freeze_face()
        if g >= cutoff:
            return 1
    return 0


def simple_connected_graph_genus(G, set_embedding=False, check=True, minimal=True):
    """
    Compute the genus of a simple connected graph.

    .. WARNING::

        THIS MAY SEGFAULT OR HANG ON:
            * DISCONNECTED GRAPHS
            * DIRECTED GRAPHS
            * LOOPED GRAPHS
            * MULTIGRAPHS

        DO NOT CALL WITH ``check = False`` UNLESS YOU ARE CERTAIN.

    EXAMPLES::

        sage: import sage.graphs.genus
        sage: from sage.graphs.genus import simple_connected_graph_genus as genus
        sage: [genus(g) for g in graphs(6) if g.is_connected()].count(1)
        13
        sage: G = graphs.FlowerSnark()
        sage: genus(G)  # see [1]
        2
        sage: G = graphs.BubbleSortGraph(4)
        sage: genus(G)
        0
        sage: G = graphs.OddGraph(3)
        sage: genus(G)
        1

    REFERENCES:

    [1] http://www.springerlink.com/content/0776127h0r7548v7/
    """
    cdef int style, cutoff
    oG = G  # original graph

    if minimal and G.is_planar(set_embedding=set_embedding):
        return 0

    if check:
        if not G.is_connected():
            raise ValueError("Cannot compute the genus of a disconnected graph")

        if G.is_directed() or G.has_multiple_edges() or G.has_loops():
            G = G.to_simple()

    G, vmap = G.relabel(inplace=False, return_map=True)
    backmap = {u: v for v, u in vmap.items()}
    G = Graph(G, sparse=False)
    GG = simple_connected_genus_backtracker(G._backend.c_graph()[0])

    if minimal:
        style = 1
        cutoff = 1
    else:
        style = 2
        cutoff = 1 + (G.num_edges() - G.num_verts()) / 2  # rounding here is ok

    g = GG.genus(style=style, cutoff=cutoff, record_embedding=set_embedding)
    if set_embedding:
        oE = {}
        E = GG.get_embedding()
        for v in E:
            oE[backmap[v]] = [backmap[x] for x in E[v]]
        oG.set_embedding(oE)
    return g
