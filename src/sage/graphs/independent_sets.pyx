r"""
Independent sets

This module implements the :class:`IndependentSets` class which can be used to :

- List the independent sets (or cliques) of a graph
- Count them (which is obviously faster)
- Test whether a set of vertices is an independent set

It can also be restricted to focus on (inclusionwise) maximal independent
sets. See the documentation of :class:`IndependentSets` for actual examples.

Classes and methods
-------------------

"""
from __future__ import print_function

include "sage/data_structures/binary_matrix.pxi"
from sage.misc.cachefunc import cached_method
from sage.graphs.base.static_dense_graph cimport dense_graph_init


cdef inline int ismaximal(binary_matrix_t g, int n, bitset_t s):
    cdef int i
    for i in range(n):
        if (not bitset_in(s,i)) and bitset_are_disjoint(g.rows[i], s):
            return False

    return True

cdef class IndependentSets:
    r"""
    The set of independent sets of a graph.

    For more information on independent sets, see
    :wikipedia:`Independent_set_(graph_theory)`.

    INPUT:

    - ``G`` -- a graph

    - ``maximal`` (boolean) -- whether to only consider (inclusionwise) maximal
      independent sets. Set to ``False`` by default.

    - ``complement`` (boolean) -- whether to consider the graph's complement
      (i.e. cliques instead of independent sets). Set to ``False`` by default.

    ALGORITHM:

    The enumeration of independent sets is done naively : given an independent
    set, this implementation considers all ways to add a new vertex to it
    (while keeping it an independent set), and then creates new independent
    sets from all those that were created this way.

    The implementation, however, is not recursive.

    .. NOTE::

        This implementation of the enumeration of *maximal* independent sets is
        not much faster than NetworkX', which is surprising as it is written in
        Cython. This being said, the algorithm from NetworkX appears to be
        sligthly different from this one, and that would be a good thing to
        explore if one wants to improve the implementation.

        A simple generalization can also be done without too much modifications:
        iteration through independent sets with given size bounds
        (minimum and maximum number of vertices allowed).

    EXAMPLES:

    Listing all independent sets of the Claw graph::

        sage: from sage.graphs.independent_sets import IndependentSets
        sage: g = graphs.ClawGraph()
        sage: I = IndependentSets(g)
        sage: list(I)
        [[0], [1], [1, 2], [1, 2, 3], [1, 3], [2], [2, 3], [3], []]

    Count them::

        sage: I.cardinality()
        9

    List only the maximal independent sets::

        sage: Im = IndependentSets(g, maximal = True)
        sage: list(Im)
        [[0], [1, 2, 3]]

    And count them::

        sage: Im.cardinality()
        2

    One can easily count the number of independent sets of each
    cardinality::

        sage: g = graphs.PetersenGraph()
        sage: number_of = [0] * g.order()
        sage: for x in IndependentSets(g):
        ....:     number_of[len(x)] += 1
        sage: number_of
        [1, 10, 30, 30, 5, 0, 0, 0, 0, 0]

    It is also possible to define an an iterator over all independent sets of a
    given cardinality. Note, however, that Sage will generate them *all*, to
    return only those that satisfy the cardinality constraints. Getting the list
    of independent sets of size 4 in this way can thus take a very long time::

        sage: is4 = (x for x in IndependentSets(g) if len(x) == 4)
        sage: list(is4)
        [[0, 2, 8, 9], [0, 3, 6, 7], [1, 3, 5, 9], [1, 4, 7, 8], [2, 4, 5, 6]]

    Given a subset of the vertices, it is possible to test whether it is an
    independent set::

        sage: g = graphs.DurerGraph()
        sage: I = IndependentSets(g)
        sage: [0,2] in I
        True
        sage: [0,3,5] in I
        False

    If an element of the subset is not a vertex, then an error is raised::

        sage: [0, 'a', 'b', 'c'] in I
        Traceback (most recent call last):
        ...
        ValueError: a is not a vertex of the graph.
    """
    def __init__(self, G, maximal = False, complement = False):
        r"""
        Constructor for this class

        TESTS::

            sage: from sage.graphs.independent_sets import IndependentSets
            sage: IndependentSets(graphs.PetersenGraph())
            <sage.graphs.independent_sets.IndependentSets object...

        Compute the number of matchings, and check with Sage's implementation::

            sage: from sage.graphs.independent_sets import IndependentSets
            sage: from sage.graphs.matchpoly import matching_polynomial
            sage: def check_matching(G):
            ...       number_of_matchings = sum(map(abs,matching_polynomial(G).coefficients(sparse=False)))
            ...       if number_of_matchings != IndependentSets(G.line_graph()).cardinality():
            ...           print("Ooooch !")
            sage: for i in range(30):
            ...       check_matching(graphs.RandomGNP(11,.3))

        Compare the result with the output of :meth:`subgraph_search`::

            sage: from sage.sets.set import Set
            sage: def check_with_subgraph_search(G):
            ...       IS = set(map(Set,list(IndependentSets(G))))
            ...       if not all(G.subgraph(l).is_independent_set() for l in IS):
            ...          print("Gloops")
            ...       alpha = max(map(len,IS))
            ...       IS2 = [Set([x]) for x in range(G.order())] + [Set([])]
            ...       for n in range(2,alpha+1):
            ...           IS2.extend(map(Set,list(G.subgraph_search_iterator(Graph(n), induced = True))))
            ...       if len(IS) != len(set(IS2)):
            ...          print("Oops")
            ...          print(len(IS), len(set(IS2)))
            sage: for i in range(5):
            ...       check_with_subgraph_search(graphs.RandomGNP(11,.3))

        Empty graph::

            sage: IS0 = IndependentSets(graphs.EmptyGraph())
            sage: list(IS0)
            [[]]
            sage: IS0.cardinality()
            1
        """
        cdef int i

        # Map from Vertex to Integer, and from Integer to Vertex
        self.vertices = G.vertices()
        self.n = G.order()
        self.maximal = maximal
        self.vertex_to_int = dense_graph_init(self.g, G, translation = True)

        # If we must consider the graph's complement instead
        if complement:
            binary_matrix_complement(self.g)
            for i in range(self.n):
                binary_matrix_set0(self.g,i,i)

        self.count_only = 0

    def __iter__(self):
        r"""
        Returns an iterator over the independent sets of self.

        TESTS::

            sage: from sage.graphs.independent_sets import IndependentSets
            sage: I = IndependentSets(graphs.PetersenGraph())
            sage: iter1 = iter(I)
            sage: iter2 = iter(I)
            sage: next(iter1)      # indirect doctest
            [0]
            sage: next(iter2)      # indirect doctest
            [0]
            sage: next(iter2)
            [0, 2]
            sage: next(iter1)
            [0, 2]
        """
        if self.n == 0:
            yield []
            return

        cdef int i = 0

        cdef bitset_t current_set
        cdef bitset_t tmp
        bitset_init(current_set,self.n)
        bitset_set_first_n(current_set,0)
        bitset_add(current_set,0)
        bitset_init(tmp,self.n)

        cdef uint64_t count = 0
        cdef list ans
        cdef int j

        # At every moment of the algorithm current_set represents an independent
        # set, except for the ith bit. All bits >i are zero.

        while True:

            # If i is in current_set
            if bitset_in(current_set,i):

                # We have found an independent set !
                if bitset_are_disjoint(self.g.rows[i], current_set):

                    # Saving that set
                    bitset_copy(tmp, current_set)

                    # Preparing for the next set, except if we set the last bit.
                    if i < self.n-1:

                        # Adding (i+1)th bit
                        bitset_add(current_set,i+1)
                        i += 1
                    else:
                        bitset_discard(current_set,i)

                    # Returning the result if necessary ...
                    if self.maximal and not ismaximal(self.g,self.n, tmp):
                        continue

                    count += 1

                    if not self.count_only:
                        yield [self.vertices[j] for j in range(i+1) if bitset_in(tmp,j)]
                        continue

                else:
                    # Removing the ith bit
                    bitset_discard(current_set, i)

                    # Preparing for the next set !
                    if i < self.n-1:
                        bitset_add(current_set, i+1)
                        i += 1

            # Not already included in the set
            else:
                if i == 0:
                    break

                # Going backward, we explored all we could there !
                if bitset_in(current_set,i-1):
                    bitset_discard(current_set, i-1)
                    bitset_add(current_set,i)
                else:
                    i -= 1
                    if i == -1:
                        break

        if not self.maximal:
            count += 1
            if not self.count_only:
                yield []

        if self.count_only:
            yield count

        bitset_free(current_set)
        bitset_free(tmp)

    def __dealloc__(self):
        r"""
        Frees everything we ever allocated
        """
        if self.g.rows != NULL:
            binary_matrix_free(self.g)

    @cached_method
    def cardinality(self):
        r"""
        Computes and returns the number of independent sets

        TESTS::

            sage: from sage.graphs.independent_sets import IndependentSets
            sage: IndependentSets(graphs.PetersenGraph()).cardinality()
            76

        Only maximal ones::

            sage: from sage.graphs.independent_sets import IndependentSets
            sage: IndependentSets(graphs.PetersenGraph(), maximal = True).cardinality()
            15
        """
        if self.n == 0:
            return 1

        self.count_only = 1

        for i in self:
            pass

        self.count_only = 0

        from sage.rings.integer import Integer
        return Integer(i)

    def __contains__(self, S):
        r"""
        Checks whether the set is an independent set (possibly maximal)

        INPUT:

        - ``S`` -- a set of vertices to be tested.

        TESTS:

        All independent sets of PetersenGraph are... independent sets::

            sage: from sage.graphs.independent_sets import IndependentSets
            sage: G = graphs.PetersenGraph()
            sage: IS = IndependentSets(graphs.PetersenGraph())
            sage: all(s in IS for s in IS)
            True

        And only them are::

            sage: IS2 = [x for x in subsets(G.vertices()) if x in IS]
            sage: sorted(IS) == sorted(IS2)
            True

        Same with maximal independent sets::

            sage: IS = IndependentSets(graphs.PetersenGraph(), maximal = True)
            sage: S = Subsets(G.vertices())
            sage: all(s in IS for s in IS)
            True
            sage: IS2 = [x for x in subsets(G.vertices()) if x in IS]
            sage: sorted(IS) == sorted(IS2)
            True
        """
        # Set of vertices as a bitset
        cdef bitset_t s

        bitset_init(s,self.n)
        bitset_set_first_n(s,0)

        cdef int i
        for I in S:
            try:
                i = self.vertex_to_int[I]
            except KeyError:
                raise ValueError(str(I)+" is not a vertex of the graph.")

            # Adding the new vertex to s
            bitset_add(s, i)

            # Checking that the set s is independent
            if not bitset_are_disjoint(self.g.rows[i], s):
                return False

        if self.maximal and not ismaximal(self.g, self.n,s):
            return False

        return True
