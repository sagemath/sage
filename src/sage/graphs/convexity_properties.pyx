# cython: binding=True
r"""
Convexity properties of graphs

This class gathers the algorithms related to convexity in a graph. It implements
the following methods:

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`ConvexityProperties.hull` | Return the convex hull of a set of vertices
    :meth:`ConvexityProperties.hull_number` | Compute the hull number of a graph and a corresponding generating set
    :meth:`geodetic_closure`| Return the geodetic closure of a set of vertices
    :meth:`is_geodetic` | Check whether the input (di)graph is geodetic

Some of these methods can be used through the :class:`ConvexityProperties`
object returned by :meth:`Graph.convexity_properties`.

Methods
-------
"""

# ****************************************************************************
#       Copyright (C) 2011 Nathann Cohen <nathann.cohen@gmail.com>
#                     2021 David Coudert <david.coudert@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.data_structures.binary_matrix cimport *
from sage.numerical.backends.generic_backend cimport GenericBackend
from sage.numerical.backends.generic_backend import get_solver
from sage.graphs.distances_all_pairs cimport c_distances_all_pairs
from cysignals.memory cimport sig_free
from memory_allocator cimport MemoryAllocator
from libc.stdint cimport uint32_t
from sage.graphs.base.static_sparse_graph cimport (short_digraph,
                                                   init_short_digraph,
                                                   free_short_digraph,
                                                   out_degree,
                                                   simple_BFS)


cdef class ConvexityProperties:
    r"""
    This class gathers the algorithms related to convexity in a graph.

    **Definitions**

    A set `S \subseteq V(G)` of vertices is said to be convex if for all `u,v\in
    S` the set `S` contains all the vertices located on a shortest path between
    `u` and `v`. Alternatively, a set `S` is said to be convex if the distances
    satisfy `\forall u,v\in S, \forall w\in V\backslash S : d_{G}(u,w) +
    d_{G}(w,v) > d_{G}(u,v)`.

    The convex hull `h(S)` of a set `S` of vertices is defined as the smallest
    convex set containing `S`.

    It is a closure operator, as trivially `S\subseteq h(S)` and `h(h(S)) =
    h(S)`.

    **What this class contains**

    As operations on convex sets generally involve the computation of distances
    between vertices, this class' purpose is to cache that information so that
    computing the convex hulls of several different sets of vertices does not
    imply recomputing several times the distances between the vertices.

    In order to compute the convex hull of a set `S` it is possible to write the
    following algorithm:

        For any pair `u,v` of elements in the set `S`, and for any vertex `w`
        outside of it, add `w` to `S` if `d_{G}(u,w) + d_{G}(w,v) =
        d_{G}(u,v)`. When no vertex can be added anymore, the set `S` is convex

    The distances are not actually that relevant. The same algorithm can be
    implemented by remembering for each pair `u, v` of vertices the list of
    elements `w` satisfying the condition, and this is precisely what this class
    remembers, encoded as bitsets to make storage and union operations more
    efficient.

    .. NOTE::

        * This class is useful if you compute the convex hulls of many sets in
          the same graph, or if you want to compute the hull number itself as it
          involves many calls to :meth:`hull`

        * Using this class on non-connected graphs is a waste of space and
          efficiency ! If your graph is disconnected, the best for you is to
          deal independently with each connected component, whatever you are
          doing.

    **Possible improvements**

    When computing a convex set, all the pairs of elements belonging to the set
    `S` are enumerated several times.

    * There should be a smart way to avoid enumerating pairs of vertices which
      have already been tested. The cost of each of them is not very high, so
      keeping track of those which have been tested already may be too expensive
      to gain any efficiency.

    * The ordering in which they are visited is currently purely lexicographic,
      while there is a Poset structure to exploit. In particular, when two
      vertices `u, v` are far apart and generate a set `h(\{u,v\})` of vertices,
      all the pairs of vertices `u', v'\in h(\{u,v\})` satisfy `h(\{u',v'\})
      \subseteq h(\{u,v\})`, and so it is useless to test the pair `u', v'` when
      both `u` and `v` where present.

    * The information cached is for any pair `u,v` of vertices the list of
      elements `z` with `d_{G}(u,w) + d_{G}(w,v) = d_{G}(u,v)`. This is not in
      general equal to `h(\{u,v\})` !

    Nothing says these recommendations will actually lead to any actual
    improvements. There are just some ideas remembered while writing this
    code. Trying to optimize may well lead to lost in efficiency on many
    instances.

    EXAMPLES::

        sage: from sage.graphs.convexity_properties import ConvexityProperties
        sage: g = graphs.PetersenGraph()
        sage: CP = ConvexityProperties(g)
        sage: CP.hull([1, 3])
        [1, 2, 3]
        sage: CP.hull_number()                                                          # needs sage.numerical.mip
        3

    TESTS::

        sage: ConvexityProperties(digraphs.Circuit(5))
        Traceback (most recent call last):
        ...
        NotImplementedError: this is currently implemented for Graphs only, but only minor updates are needed if you want to make it support DiGraphs too
    """

    def __init__(self, G):
        r"""
        Constructor.

        EXAMPLES::

            sage: from sage.graphs.convexity_properties import ConvexityProperties
            sage: g = graphs.PetersenGraph()
            sage: ConvexityProperties(g)
            <sage.graphs.convexity_properties.ConvexityProperties object at ...>
        """
        from sage.graphs.digraph import DiGraph
        if isinstance(G, DiGraph):
            raise NotImplementedError("this is currently implemented for Graphs only, "
                                      "but only minor updates are needed if you want "
                                      "to make it support DiGraphs too")

        # Cached number of vertices
        cdef int n = G.order()
        self._n = n

        cdef int i = 0
        cdef int j, k

        # Build mappings integer <-> vertices.
        # Must be consistent with the mappings used in c_distances_all_pairs
        self._list_integers_to_vertices = list(G)
        self._dict_vertices_to_integers = {v: i for i, v in enumerate(self._list_integers_to_vertices)}

        # Computation of distances between all pairs. Costly.
        cdef unsigned short* c_distances = c_distances_all_pairs(G, vertex_list=self._list_integers_to_vertices)
        # Temporary variables
        cdef unsigned short* d_i
        cdef unsigned short* d_j
        cdef int d_ij

        # We use a binary matrix with one row per pair of vertices,
        # so n * (n - 1) / 2 rows. Row u * n + v is a bitset whose 1 bits are
        # the vertices located on a shortest path from vertex u to v
        #
        # Note that  u < v
        binary_matrix_init(self._cache_hull_pairs, n * (n - 1) / 2, n)
        binary_matrix_fill(self._cache_hull_pairs, 0)
        cdef bitset_t * p_bitset = self._cache_hull_pairs.rows

        # Filling the cache
        #
        # The p_bitset variable iterates over the successive elements of the cache
        #
        # For any pair i, j of vertices (i < j), we built the bitset of all the
        # elements k which are on a shortest path from i to j

        for i in range(n):
            # Caching the distances from i to the other vertices
            d_i = c_distances + n * i

            for j in range(i + 1, n):
                # Caching the distances from j to the other vertices
                d_j = c_distances + n * j

                # Caching the distance between i and j
                d_ij = d_i[j]

                # Filling it
                for k in range(n):
                    if d_i[k] + d_j[k] == d_ij:
                        bitset_add(p_bitset[0], k)

                # Next bitset !
                p_bitset = p_bitset + 1

        sig_free(c_distances)

    def __dealloc__(self):
        r"""
        Destructor.

        EXAMPLES::

            sage: from sage.graphs.convexity_properties import ConvexityProperties
            sage: g = graphs.PetersenGraph()
            sage: ConvexityProperties(g)
            <sage.graphs.convexity_properties.ConvexityProperties object at ...>
        """
        binary_matrix_free(self._cache_hull_pairs)

    cdef list _vertices_to_integers(self, vertices):
        r"""
        Convert a list of vertices to a list of integers with the cached data.
        """
        return [self._dict_vertices_to_integers[v] for v in vertices]

    cdef list _integers_to_vertices(self, list integers):
        r"""
        Convert a list of integers to a list of vertices with the cached data.
        """
        cdef int i
        return [self._list_integers_to_vertices[i] for i in integers]

    cdef _bitset_convex_hull(self, bitset_t hull):
        r"""
        Compute the convex hull of a list of vertices given as a bitset.

        (this method returns nothing and modifies the input)
        """
        cdef int count
        cdef int tmp_count
        cdef int i, j

        cdef bitset_t * p_bitset

        # Current size of the set
        count = bitset_len(hull)

        while True:

            # Iterating over all the elements in the cache
            p_bitset = self._cache_hull_pairs.rows

            # For any vertex i
            for i in range(self._n - 1):

                # If i is not in the current set, we skip it !
                if not bitset_in(hull, i):
                    p_bitset = p_bitset + (self._n - 1 - i)
                    continue

                # If it is, we iterate over all the elements j
                for j in range(i + 1, self._n):

                    # If both i and j are inside, we add all the (cached)
                    # vertices on a shortest ij-path

                    if bitset_in(hull, j):
                        bitset_union(hull, hull, p_bitset[0])

                    # Next bitset !
                    p_bitset = p_bitset + 1

            tmp_count = bitset_len(hull)

            # If we added nothing new during the previous loop, our set is
            # convex !
            if tmp_count == count:
                return

            # Otherwise, update and back to the loop
            count = tmp_count

    cpdef hull(self, list vertices):
        r"""
        Return the convex hull of a set of vertices.

        INPUT:

        - ``vertices`` -- list of vertices

        EXAMPLES::

            sage: from sage.graphs.convexity_properties import ConvexityProperties
            sage: g = graphs.PetersenGraph()
            sage: CP = ConvexityProperties(g)
            sage: CP.hull([1, 3])
            [1, 2, 3]
        """
        cdef bitset_t bs
        bitset_init(bs, self._n)
        bitset_set_first_n(bs, 0)

        for v in vertices:
            bitset_add(bs, <mp_bitcnt_t> self._dict_vertices_to_integers[v])

        self._bitset_convex_hull(bs)

        cdef list answer = self._integers_to_vertices(bitset_list(bs))

        bitset_free(bs)

        return answer

    cdef _greedy_increase(self, bitset_t bs):
        r"""
        Given a bitset whose hull is not the whole set, greedily add vertices
        and stop before its hull is the whole set.

        .. NOTE::

            * Counting the bits at each turn is not the best way...
        """
        cdef bitset_t tmp
        bitset_init(tmp, self._n)

        for i in range(self._n):
            if not bitset_in(bs, i):
                bitset_copy(tmp, bs)
                bitset_add(tmp, i)
                self._bitset_convex_hull(tmp)
                if bitset_len(tmp) < self._n:
                    bitset_add(bs, i)

        bitset_free(tmp)

    cpdef hull_number(self, value_only=True, verbose=False):
        r"""
        Compute the hull number and a corresponding generating set.

        The hull number `hn(G)` of a graph `G` is the cardinality of a smallest
        set of vertices `S` such that `h(S)=V(G)`.

        INPUT:

        - ``value_only`` -- boolean (default: ``True``); whether to return only
          the hull number (default) or a minimum set whose convex hull is the
          whole graph

        - ``verbose`` -- boolean (default: ``False``); whether to display
          information on the LP

        **COMPLEXITY:**

        This problem is NP-Hard [HLT1993]_, but seems to be of the "nice" kind.
        Update this comment if you fall on hard instances `:-)`

        **ALGORITHM:**

        This is solved by linear programming.

        As the function `h(S)` associating to each set `S` its convex hull is a
        closure operator, it is clear that any set `S_G` of vertices such that
        `h(S_G)=V(G)` must satisfy `S_G \not \subseteq C` for any *proper*
        convex set `C \subsetneq V(G)`. The following formulation is hence
        correct

        .. MATH::

            \text{Minimize :}& \sum_{v\in G}b_v\\
            \text{Such that :}&\\
            &\forall C\subsetneq V(G)\text{ a proper convex set }\\
            &\sum_{v\in V(G)\backslash C} b_v \geq 1

        Of course, the number of convex sets -- and so the number of constraints
        -- can be huge, and hard to enumerate, so at first an incomplete
        formulation is solved (it is missing some constraints). If the answer
        returned by the LP solver is a set `S` generating the whole graph, then
        it is optimal and so is returned. Otherwise, the constraint
        corresponding to the set `h(S)` can be added to the LP, which makes the
        answer `S` infeasible, and another solution computed.

        This being said, simply adding the constraint corresponding to `h(S)` is
        a bit slow, as these sets can be large (and the corresponding constraint
        a bit weak). To improve it a bit, before being added, the set `h(S)` is
        "greedily enriched" to a set `S'` with vertices for as long as
        `h(S')\neq V(G)`. This way, we obtain a set `S'` with `h(S)\subseteq
        h(S')\subsetneq V(G)`, and the constraint corresponding to `h(S')` --
        which is stronger than the one corresponding to `h(S)` -- is added.

        This can actually be seen as a hitting set problem on the complement of
        convex sets.

        EXAMPLES:

        The Hull number of Petersen's graph::

            sage: from sage.graphs.convexity_properties import ConvexityProperties
            sage: g = graphs.PetersenGraph()
            sage: CP = ConvexityProperties(g)
            sage: CP.hull_number()                                                      # needs sage.numerical.mip
            3
            sage: generating_set = CP.hull_number(value_only=False)                     # needs sage.numerical.mip
            sage: CP.hull(generating_set)                                               # needs sage.numerical.mip
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        """
        cdef int i
        cdef list constraint  # temporary variable to add constraints to the LP

        if self._n <= 2:
            if value_only:
                return self._n
            else:
                return self._list_integers_to_vertices

        cdef GenericBackend p = <GenericBackend> get_solver(constraint_generation=True)

        # Minimization
        p.set_sense(False)

        # We have exactly n binary variables, all of them with a coefficient of
        # 1 in the objective function
        p.add_variables(self._n, 0, None, True, False, False, 1, None)

        # We know that at least 2 vertices are required to cover the whole graph
        p.add_linear_constraint([(i, 1) for i in range(self._n)], 2, None)

        # The set of vertices generated by the current LP solution
        cdef bitset_t current_hull
        bitset_init(current_hull, self._n)

        # Which is at first empty
        bitset_set_first_n(current_hull, 1)

        while True:

            # Greedily increase it to obtain a better constraint
            self._greedy_increase(current_hull)

            if verbose:
                print("Adding a constraint corresponding to convex set ",
                      end="")
                print(bitset_list(current_hull))

            # Building the corresponding constraint
            constraint = []
            for i in range(self._n):
                if not bitset_in(current_hull, i):
                    constraint.append((i, 1))

            p.add_linear_constraint(constraint, 1, None)

            p.solve()

            # Computing the current solution's convex hull
            bitset_set_first_n(current_hull, 0)

            for i in range(self._n):
                if p.get_variable_value(i) > .5:
                    bitset_add(current_hull, i)

            self._bitset_convex_hull(current_hull)

            # Are we done ?
            if bitset_len(current_hull) == self._n:
                break

        bitset_free(current_hull)

        if value_only:
            return <int> p.get_objective_value()

        constraint = []
        for i in range(self._n):
            if p.get_variable_value(i) > .5:
                constraint.append(i)

        return self._integers_to_vertices(constraint)


def geodetic_closure(G, S):
    r"""
    Return the geodetic closure of the set of vertices `S` in `G`.

    The geodetic closure `g(S)` of a subset of vertices `S` of a graph `G` is in
    [HLT1993]_ as the set of all vertices that lie on a shortest `u-v` path for
    any pair of vertices `u,v \in S`. We assume that `g(\emptyset) = \emptyset`
    and that `g(\{u\}) = \{u\}` for any `u` in `G`.

    .. WARNING::

        This operation is **not** a closure function. Indeed, a closure function
        must satisfy the property that `f(f(X))` should be equal to `f(X)`,
        which is not always the case here.  The term ``closure`` is used here to
        follow the terminology of the domain. See for instance [HLT1993]_.

    Here, we implement a simple algorithm to determine this set. Roughly, for
    each vertex `u \in S`, the algorithm first performs a breadth first search
    from `u` to get distances, and then identifies the vertices of `G` lying on
    a shortest path from `u` to any `v\in S` using a reversal traversal from
    vertices in `S`.  This algorithm has time complexity in `O(|S|(n + m))` for
    ``SparseGraph``, `O(|S|(n + m) + n^2)` for ``DenseGraph`` and
    space complexity in `O(n + m)`.

    INPUT:

    - ``G`` -- a Sage graph

    - ``S`` -- a subset of vertices of `G`

    EXAMPLES:

    The vertices of the Petersen graph can be obtained by a geodetic closure of
    four of its vertices::

        sage: from sage.graphs.convexity_properties import geodetic_closure
        sage: G = graphs.PetersenGraph()
        sage: geodetic_closure(G, [0, 2, 8, 9])
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

    The vertices of a 2D grid can be obtained by a geodetic closure of
    two vertices::

        sage: G = graphs.Grid2dGraph(4, 4)
        sage: c = G.geodetic_closure([(0, 0), (3, 3)])
        sage: len(c) == G.order()
        True

    If two vertices belong to different connected components of a graph, their
    geodetic closure is trivial::

        sage: G = Graph([(0, 1), (2, 3)])
        sage: geodetic_closure(G, [0, 2])
        [0, 2]

    The geodetic closure does not satisfy the closure function property that
    `f(f(X))` should be equal to `f(X)`::

        sage: G = graphs.DiamondGraph()
        sage: G.subdivide_edge((1, 2), 1)
        sage: geodetic_closure(G, [0, 3])
        [0, 1, 2, 3]
        sage: geodetic_closure(G, geodetic_closure(G, [0, 3]))
        [0, 1, 2, 3, 4]

    TESTS::

        sage: G = graphs.DiamondGraph()
        sage: geodetic_closure(G, [])
        []
        sage: geodetic_closure(G, [1])
        [1]
        sage: S = geodetic_closure(G, list(G))
        sage: all(u in G for u in S) and len(S) == G.order()
        True
        sage: S = geodetic_closure(G, G)
        sage: all(u in G for u in S) and len(S) == G.order()
        True
        sage: geodetic_closure(G, [1, 'foo'])
        Traceback (most recent call last):
        ...
        ValueError: S is not a subset of vertices of the graph
        sage: geodetic_closure(digraphs.Path(3), [0, 1])
        Traceback (most recent call last):
        ...
        NotImplementedError: the geodetic closure of digraphs has not been implemented yet
    """
    if G.is_directed():
        raise NotImplementedError("the geodetic closure of digraphs has not been implemented yet")
    S = list(S)
    if not S:
        return []
    elif any(v not in G for v in S):
        raise ValueError("S is not a subset of vertices of the graph")
    elif len(S) == 1 or len(S) == G.order():
        return S
    elif not G.is_connected():
        L = []
        for g in G.connected_components_subgraphs():
            Sg = [u for u in S if u in g]
            L.extend(geodetic_closure(g, Sg))
        return L

    cdef int n = G.order()
    cdef int nS = len(S)
    cdef list int_to_vertex = list(G)
    cdef dict vertex_to_int = {u: i for i, u in enumerate(int_to_vertex)}
    cdef list S_int = [vertex_to_int[u] for u in S]

    # Copy the graph as a short digraph
    cdef short_digraph sd
    init_short_digraph(sd, G, edge_labelled=False, vertex_list=int_to_vertex)

    # Allocate some data structures
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef uint32_t * distances = <uint32_t *> mem.malloc(n * sizeof(uint32_t))
    cdef uint32_t * waiting_list = <uint32_t *> mem.malloc(n * sizeof(uint32_t))
    if not distances or not waiting_list:
        free_short_digraph(sd)
        raise MemoryError()
    cdef bitset_t seen
    cdef bitset_t visited
    cdef bitset_t closure
    bitset_init(seen, n)
    bitset_init(visited, n)
    bitset_init(closure, n)

    cdef int ui, vi, xi, yi, i_begin, i_end, i, j, k

    # Vertices in S are in the closure
    bitset_clear(closure)
    for ui in S_int:
        bitset_add(closure, ui)

    # We now explore geodesics between vertices in S, and we avoid visiting
    # twice the geodesics between u and v
    for i in range(nS - 1):
        ui = S_int[i]

        # Compute distances from ui using BFS
        _ = simple_BFS(sd, ui, distances, NULL, waiting_list, seen)

        # Perform a reverse BFS from vertices in S, following only vertices
        # along a shortest path from ui to identify vertices of the geodetic
        # closure.
        bitset_clear(visited)
        i_begin = 0
        i_end = 0
        for j in range(i + 1, nS):
            vi = S_int[j]
            if not bitset_in(seen, vi) or bitset_in(visited, vi):
                # vi is not reachable from ui using BFS or has already been
                # considered for the geodetic closure from ui
                continue

            # We explore all vertices on a shortest path from ui to vi
            waiting_list[i_end] = vi
            i_end += 1
            bitset_add(visited, vi)

            while i_begin < i_end:
                xi = waiting_list[i_begin]
                i_begin += 1

                for k in range(out_degree(sd, xi)):
                    yi = sd.neighbors[xi][k]
                    if distances[xi] == distances[yi] + 1:
                        # yi in on a shortest path to ui
                        if not bitset_in(visited, yi):
                            waiting_list[i_end] = yi
                            i_end += 1
                            bitset_add(visited, yi)

        # Now, visited contains all vertices on geodesics from ui to any other
        # vertex in S
        bitset_union(closure, closure, visited)

    cdef list ret = [int_to_vertex[ui] for ui in range(n) if bitset_in(closure, ui)]

    bitset_free(seen)
    bitset_free(visited)
    bitset_free(closure)
    free_short_digraph(sd)

    return ret


def is_geodetic(G):
    r"""
    Check whether the input (di)graph is geodetic.

    A graph `G` is *geodetic* if there exists only one shortest path between
    every pair of its vertices. This can be checked in time `O(nm)` for
    ``SparseGraph`` and `O(nm+n^2)` for ``DenseGraph`` in unweighted (di)graphs
    with `n` nodes and `m` edges. Examples of geodetic graphs are trees, cliques
    and odd cycles. See the :wikipedia:`Geodetic_graph` for more details.

    (Di)graphs with multiple edges are not considered geodetic.

    INPUT:

    - ``G`` -- a graph or a digraph

    EXAMPLES:

    Trees, cliques and odd cycles are geodetic::

        sage: T = graphs.RandomTree(20)
        sage: T.is_geodetic()
        True
        sage: all(graphs.CompleteGraph(n).is_geodetic() for n in range(8))
        True
        sage: all(graphs.CycleGraph(n).is_geodetic() for n in range(3, 16, 2))
        True

    Even cycles of order at least 4 are not geodetic::

        sage: all(graphs.CycleGraph(n).is_geodetic() for n in range(4, 17, 2))
        False

    The Petersen graph is geodetic::

        sage: P = graphs.PetersenGraph()
        sage: P.is_geodetic()
        True

    Grid graphs are not geodetic::

        sage: G = graphs.Grid2dGraph(2, 3)
        sage: G.is_geodetic()
        False

    This method is also valid for digraphs::

        sage: G = DiGraph(graphs.PetersenGraph())
        sage: G.is_geodetic()
        True
        sage: G = digraphs.Path(5)
        sage: G.add_path([0, 'a', 'b', 'c', 4])
        sage: G.is_geodetic()
        False

    TESTS::

        sage: all(g.is_geodetic() for g in graphs(3))
        True
        sage: all((2*g).is_geodetic() for g in graphs(3))
        True
        sage: G = graphs.CycleGraph(5)
        sage: G.allow_loops(True)
        sage: G.add_edges([(u, u) for u in G])
        sage: G.is_geodetic()
        True
        sage: G.allow_multiple_edges(True)
        sage: G.is_geodetic()
        True
        sage: G.add_edge(G.random_edge())
        sage: G.is_geodetic()
        False
    """
    if G.has_multiple_edges():
        return False

    if G.order() < 4:
        return True

    # Copy the graph as a short digraph
    cdef int n = G.order()
    cdef short_digraph sd
    init_short_digraph(sd, G, edge_labelled=False, vertex_list=list(G))

    # Allocate some data structures
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef uint32_t * distances = <uint32_t *> mem.malloc(n * sizeof(uint32_t))
    cdef uint32_t * waiting_list = <uint32_t *> mem.malloc(n * sizeof(uint32_t))
    if not distances or not waiting_list:
        free_short_digraph(sd)
        raise MemoryError()
    cdef bitset_t seen
    bitset_init(seen, n)

    # We now explore geodesics between vertices in S, and we avoid visiting
    # twice the geodesics between u and v

    cdef uint32_t source, u, v
    cdef uint32_t waiting_beginning
    cdef uint32_t waiting_end
    cdef uint32_t * p_tmp
    cdef uint32_t * end
    cdef uint32_t ** p_vertices = sd.neighbors

    for source in range(n):

        # Compute distances from source using BFS
        bitset_clear(seen)
        bitset_add(seen, source)
        distances[source] = 0
        waiting_beginning = 0
        waiting_end = 0
        waiting_list[waiting_beginning] = source

        # For as long as there are vertices left to explore
        while waiting_beginning <= waiting_end:

            # We pick the first one
            v = waiting_list[waiting_beginning]
            p_tmp = p_vertices[v]
            end = p_vertices[v + 1]

            # and we iterate over all the outneighbors u of v
            while p_tmp < end:
                u = p_tmp[0]

                # If we notice one of these neighbors is not seen yet, we set
                # its parameters and add it to the queue to be explored later.
                # Otherwise, we check whether we have detected a second shortest
                # path between source and v.
                if not bitset_in(seen, u):
                    distances[u] = distances[v] + 1
                    bitset_add(seen, u)
                    waiting_end += 1
                    waiting_list[waiting_end] = u
                elif distances[u] == distances[v] + 1:
                    # G is not geodetic
                    bitset_free(seen)
                    free_short_digraph(sd)
                    return False

                p_tmp += 1

            # We go to the next vertex in the queue
            waiting_beginning += 1

    bitset_free(seen)
    free_short_digraph(sd)

    # The graph is geodetic
    return True
