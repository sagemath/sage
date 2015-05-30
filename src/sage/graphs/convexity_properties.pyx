r"""
Convexity properties of graphs

This class gathers the algorithms related to convexity in a graph. It implements
the following methods:

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`ConvexityProperties.hull` | Returns the convex hull of a set of vertices
    :meth:`ConvexityProperties.hull_number` | Computes the hull number of a graph and a corresponding generating set.

These methods can be used through the :class:`ConvexityProperties` object
returned by :meth:`Graph.convexity_properties`.

AUTHORS:

    -  Nathann Cohen

Methods
-------
"""

##############################################################################
#       Copyright (C) 2011 Nathann Cohen <nathann.cohen@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

include "sage/data_structures/bitset.pxi"
from sage.numerical.backends.generic_backend cimport GenericBackend
from sage.numerical.backends.generic_backend import get_solver

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
    following algorithm.

    *For any pair `u,v` of elements in the set `S`, and for any vertex `w`*
    *outside of it, add `w` to `S` if `d_{G}(u,w) + d_{G}(w,v) = d_{G}(u,v)`.*
    *When no vertex can be added anymore, the set `S` is convex*

    The distances are not actually that relevant. The same algorithm can be
    implemented by remembering for each pair `u, v` of vertices the list of
    elements `w` satisfying the condition, and this is precisely what this class
    remembers, encoded as bitsets to make storage and union operations more
    efficient.

    .. NOTE::

        * This class is useful if you compute the convex hulls of many sets in
          the same graph, or if you want to compute the hull number itself as it
          involves many calls to :meth:`hull`

        * Using this class on non-conected graphs is a waste of space and
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

    Nothing says these recommandations will actually lead to any actual
    improvements. There are just some ideas remembered while writing this
    code. Trying to optimize may well lead to lost in efficiency on many
    instances.

    EXAMPLE::

        sage: from sage.graphs.convexity_properties import ConvexityProperties
        sage: g = graphs.PetersenGraph()
        sage: CP = ConvexityProperties(g)
        sage: CP.hull([1,3])
        [1, 2, 3]
        sage: CP.hull_number()
        3

    TESTS::

        sage: ConvexityProperties(digraphs.Circuit(5))
        Traceback (most recent call last):
        ...
        ValueError: This is currenly implemented for Graphs only.Only minor updates are needed if you want to makeit support DiGraphs too.
    """

    def __init__(self, G):
        r"""
        Constructor

        EXAMPLE::

            sage: from sage.graphs.convexity_properties import ConvexityProperties
            sage: g = graphs.PetersenGraph()
            sage: ConvexityProperties(g)
            <sage.graphs.convexity_properties.ConvexityProperties object at ...>
        """
        from sage.graphs.digraph import DiGraph
        if isinstance(G, DiGraph):
            raise ValueError("This is currenly implemented for Graphs only."+
                             "Only minor updates are needed if you want to make"+
                             "it support DiGraphs too.")

        # Cached number of vertices
        cdef int n = G.order()
        self._n = n

        cdef int i = 0
        cdef int j,k

        # Temporary variables
        cdef dict d_i
        cdef dict d_j
        cdef int d_ij
        self._dict_vertices_to_integers = {}
        self._list_integers_to_vertices = []

        # Remembering integers instead of the labels, and building dictionaries
        # in both directions.
        for v in G:
            self._dict_vertices_to_integers[v] = i
            self._list_integers_to_vertices.append(v)
            i = i + 1


        # Computation of distances between all pairs. Costly.
        cdef dict distances = G.distance_all_pairs()

        # _cache_hull_pairs[u*n + v] is a bitset whose 1 bits are the vertices located on a shortest path from vertex u to v
        #
        # Note that  u < v
        self._cache_hull_pairs = <bitset_t *> sage_malloc(((n*(n-1))>>1)*sizeof(bitset_t))
        cdef bitset_t * p_bitset = self._cache_hull_pairs

        # Filling the cache
        #
        # The p_bitset variable iterates over the successive elements of the cache
        #
        # For any pair i,j of vertices (i<j), we built the bitset of all the
        # elements k which are on a shortest path from i to j

        for 0<= i < n-1:
            # Caching the distances from i to the other vertices
            d_i = distances[self._list_integers_to_vertices[i]]

            for i < j < n:
                # Caching the distances from j to the other vertices
                d_j = distances[self._list_integers_to_vertices[j]]

                # Caching the distance between i and j
                d_ij = d_i[self._list_integers_to_vertices[j]]

                # Initializing the new bitset
                bitset_init(p_bitset[0], n)
                bitset_set_first_n(p_bitset[0], 0)

                # Filling it
                for 0<= k < n:
                    if ((d_i[self._list_integers_to_vertices[k]]
                         + d_j[self._list_integers_to_vertices[k]])
                        == d_ij):
                        bitset_add(p_bitset[0], k)

                # Next bitset !
                p_bitset = p_bitset + 1


    def __destruct__(self):
        r"""
        Destructor

        EXAMPLE::

            sage: from sage.graphs.convexity_properties import ConvexityProperties
            sage: g = graphs.PetersenGraph()
            sage: ConvexityProperties(g)
            <sage.graphs.convexity_properties.ConvexityProperties object at ...>

        """
        cdef bitset_t * p_bitset = self._cache_hull_pairs
        cdef int i

        for 0 <= i < ((self._n*(self._n-1))>>1):
            bitset_free(p_bitset[0])
            p_bitset = p_bitset + 1

        sage_free(self._cache_hull_pairs)

    cdef list _vertices_to_integers(self, vertices):
        r"""
        Converts a list of vertices to a list of integers with the cached data.
        """
        cdef list answer = []
        for v in v:
            answer.append(self._dict_vertices_to_integers[v])
        return answer

    cdef list _integers_to_vertices(self, integers):
        r"""
        Converts a list of integers to a list of vertices with the cached data.
        """

        cdef list answer = []
        for v in integers:
            answer.append(self._list_integers_to_vertices[v])
        return answer

    cdef _bitset_convex_hull(self, bitset_t hull):
        r"""
        Computes the convex hull of a list of vertices given as a bitset.

        (this method returns nothing and modifies the input)
        """
        cdef int count
        cdef int tmp_count
        cdef int i,j

        cdef bitset_t * p_bitset

        # Current size of the set
        count = bitset_len(hull)

        while True:

            # Iterating over all the elements in the cache
            p_bitset = self._cache_hull_pairs

            # For any vertex i
            for 0 <= i < self._n-1:

                # If i is not in the current set, we skip it !
                if not bitset_in(hull, i):
                    p_bitset = p_bitset + (self._n-1-i)
                    continue

                # If it is, we iterate over all the elements j
                for i < j < self._n:

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
        Returns the convex hull of a set of vertices.

        INPUT:

        * ``vertices`` -- A list of vertices.

        EXAMPLE::

            sage: from sage.graphs.convexity_properties import ConvexityProperties
            sage: g = graphs.PetersenGraph()
            sage: CP = ConvexityProperties(g)
            sage: CP.hull([1,3])
            [1, 2, 3]
        """
        cdef bitset_t bs
        bitset_init(bs, self._n)
        bitset_set_first_n(bs, 0)

        for v in vertices:
            bitset_add(bs, self._dict_vertices_to_integers[v])

        self._bitset_convex_hull(bs)

        #cdef list answer = bitset_list(bs)
        cdef list answer = self._integers_to_vertices(bitset_list(bs))

        bitset_free(bs)

        return answer

    cdef _greedy_increase(self, bitset_t bs):
        r"""
        Given a bitset whose hull is not the whole set, greedily add vertices
        and stop before its hull is the whole set.

        NOTE:

        * Counting the bits at each turn is not the best way...
        """
        cdef bitset_t tmp
        bitset_init(tmp, self._n)


        for 0<= i < self._n:
            if not bitset_in(bs, i):
                bitset_copy(tmp, bs)
                bitset_add(tmp, i)
                self._bitset_convex_hull(tmp)
                if bitset_len(tmp) < self._n:
                    bitset_add(bs, i)


    cpdef hull_number(self, value_only = True, verbose = False):
        r"""
        Computes the hull number and a corresponding generating set.

        The hull number `hn(G)` of a graph `G` is the cardinality of a smallest
        set of vertices `S` such that `h(S)=V(G)`.

        INPUT:

        * ``value_only`` (boolean) -- whether to return only the hull number
          (default) or a minimum set whose convex hull is the whole graph.

        * ``verbose`` (boolean) -- whether to display information on the LP.

        **COMPLEXITY:**

        This problem is NP-Hard [CHZ02]_, but seems to be of the "nice"
        kind. Update this comment if you fall on hard instances `:-)`

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
        a bit slow, as these sets can be large (and the corresponding constrait
        a bit weak). To improve it a bit, before being added, the set `h(S)` is
        "greedily enriched" to a set `S'` with vertices for as long as
        `h(S')\neq V(G)`. This way, we obtain a set `S'` with `h(S)\subseteq
        h(S')\subsetneq V(G)`, and the constraint corresponding to `h(S')` --
        which is stronger than the one corresponding to `h(S)` -- is added.

        This can actually be seen as a hitting set problem on the complement of
        convex sets.

        EXAMPLE:

        The Hull number of Petersen's graph::

            sage: from sage.graphs.convexity_properties import ConvexityProperties
            sage: g = graphs.PetersenGraph()
            sage: CP = ConvexityProperties(g)
            sage: CP.hull_number()
            3
            sage: generating_set = CP.hull_number(value_only = False)
            sage: CP.hull(generating_set)
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

        REFERENCE:

        .. [CHZ02] F. Harary, E. Loukakis, C. Tsouros
          The geodetic number of a graph
          Mathematical and computer modelling
          vol. 17 n11 pp.89--95, 1993
        """

        cdef int i
        cdef list constraint # temporary variable to add constraints to the LP

        if self._n <= 2:
            if value_only:
                return self._n
            else:
                return self._list_integers_to_vertices

        cdef GenericBackend p = <GenericBackend> get_solver(constraint_generation = True)

        # Minimization
        p.set_sense(False)

        # We have exactly n binary variables, all of them with a coefficient of
        # 1 in the objective function
        p.add_variables(self._n, 0, None, True, False, False, 1, None)

        # We know that at least 2 vertices are required to cover the whole graph
        p.add_linear_constraint(zip(range(self._n), [1]*self._n), 2, None)

        # The set of vertices generated by the current LP solution
        cdef bitset_t current_hull
        bitset_init(current_hull, self._n)

        # Which is at first empty
        bitset_set_first_n(current_hull,1)

        while True:

            # Greedily increase it to obtain a better constraint
            self._greedy_increase(current_hull)

            if verbose:
                print "Adding a constraint corresponding to convex set ",
                print bitset_list(current_hull)

            # Building the corresponding constraint
            constraint = []
            for 0 <= i < self._n:
                if not bitset_in(current_hull, i):
                    constraint.append((i,1))

            p.add_linear_constraint(constraint, 1, None)

            p.solve()

            # Computing the current solution's convex hull
            bitset_set_first_n(current_hull,0)

            for 0 <= i < self._n:
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
        for 0 <= i < self._n:
            if p.get_variable_value(i) > .5:
                constraint.append(i)

        return self._integers_to_vertices(constraint)
