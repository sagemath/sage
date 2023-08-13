r"""
Hasse diagrams of posets

{INDEX_OF_FUNCTIONS}

"""
# ****************************************************************************
#       Copyright (C) 2008 Peter Jipsen <jipsen@chapman.edu>
#       Copyright (C) 2008 Franco Saliola <saliola@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from __future__ import annotations

from collections import deque

from sage.arith.misc import binomial
from sage.combinat.posets.hasse_cython import IncreasingChains
from sage.graphs.digraph import DiGraph
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.lazy_import import lazy_import
from sage.misc.rest_index_of_methods import gen_rest_table_index
from sage.rings.integer_ring import ZZ

lazy_import('sage.combinat.posets.hasse_cython_flint',
            ['moebius_matrix_fast', 'coxeter_matrix_fast'])
lazy_import('sage.matrix.constructor', 'matrix')
lazy_import('sage.rings.finite_rings.finite_field_constructor', 'GF')


class LatticeError(ValueError):
    """
    Helper exception class to forward elements without meet or
    join to upper level, so that the user will see "No meet for
    a and b" instead of "No meet for 1 and 2".
    """

    def __init__(self, fail, x, y):
        """
        Initialize the exception.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import LatticeError
            sage: error = LatticeError('join', 3, 8)
            sage: error.x
            3
        """
        ValueError.__init__(self, None)
        self.fail = fail
        self.x = x
        self.y = y

    def __str__(self):
        """
        Return string representation of the exception.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import LatticeError
            sage: error = LatticeError('meet', 15, 18)
            sage: error.__str__()
            'no meet for 15 and 18'
        """
        return "no {} for {} and {}".format(self.fail, self.x, self.y)


class HasseDiagram(DiGraph):
    """
    The Hasse diagram of a poset. This is just a transitively-reduced,
    directed, acyclic graph without loops or multiple edges.

    .. note::

       We assume that ``range(n)`` is a linear extension of the poset.
       That is, ``range(n)`` is the vertex set and a topological sort of
       the digraph.

    This should not be called directly, use Poset instead; all type
    checking happens there.

    EXAMPLES::

        sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
        sage: H = HasseDiagram({0:[1,2],1:[3],2:[3],3:[]}); H
        Hasse diagram of a poset containing 4 elements
        sage: TestSuite(H).run()
    """
    def _repr_(self):
        r"""
        TESTS::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,2],1:[3],2:[3],3:[]})
            sage: H._repr_()
            'Hasse diagram of a poset containing 4 elements'
        """
        return "Hasse diagram of a poset containing %s elements" % self.order()

    def linear_extension(self):
        r"""
        Return a linear extension

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,2],1:[3],2:[3],3:[]})
            sage: H.linear_extension()
            [0, 1, 2, 3]
        """
        # Recall: we assume range(n) is a linear extension.
        return list(range(len(self)))

    def linear_extensions(self):
        r"""
        Return an iterator over all linear extensions.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,2],1:[3],2:[3],3:[]})
            sage: list(H.linear_extensions())                                           # optional - sage.modules
            [[0, 1, 2, 3], [0, 2, 1, 3]]
        """
        from sage.combinat.posets.linear_extension_iterator import linear_extension_iterator
        return linear_extension_iterator(self)

    def greedy_linear_extensions_iterator(self):
        r"""
        Return an iterator over greedy linear extensions of the Hasse diagram.

        A linear extension `[e_1, e_2, \ldots, e_n]` is *greedy* if for
        every `i` either `e_{i+1}` covers `e_i` or all upper covers
        of `e_i` have at least one lower cover that is not in
        `[e_1, e_2, \ldots, e_i]`.

        Informally said a linear extension is greedy if it "always
        goes up when possible" and so has no unnecessary jumps.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: N5 = HasseDiagram({0: [1, 2], 2: [3], 1: [4], 3: [4]})
            sage: for l in N5.greedy_linear_extensions_iterator():
            ....:     print(l)
            [0, 1, 2, 3, 4]
            [0, 2, 3, 1, 4]

        TESTS::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: list(HasseDiagram({}).greedy_linear_extensions_iterator())
            [[]]
            sage: H = HasseDiagram({0: []})
            sage: list(H.greedy_linear_extensions_iterator())
            [[0]]
        """
        N = self.order()

        def greedy_rec(H, linext):
            if len(linext) == N:
                yield linext

            S = []
            if linext:
                S = [x for x in H.neighbor_out_iterator(linext[-1])
                     if all(low in linext for low in H.neighbor_in_iterator(x))]
            if not S:
                S_ = set(self).difference(set(linext))
                S = [x for x in S_
                     if not any(low in S_
                                for low in self.neighbor_in_iterator(x))]

            for e in S:
                yield from greedy_rec(H, linext + [e])

        return greedy_rec(self, [])

    def supergreedy_linear_extensions_iterator(self):
        r"""
        Return an iterator over supergreedy linear extensions of the Hasse diagram.

        A linear extension `[e_1, e_2, \ldots, e_n]` is *supergreedy* if,
        for every `i` and `j` where `i > j`, `e_i` covers `e_j` if for
        every `i > k > j` at least one lower cover of `e_k` is not in
        `[e_1, e_2, \ldots, e_k]`.

        Informally said a linear extension is supergreedy if it "always
        goes as high possible, and withdraw so less as possible".
        These are also called depth-first linear extensions.

        EXAMPLES:

        We show the difference between "only greedy" and supergreedy
        extensions::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0: [1, 2], 2: [3, 4]})
            sage: G_ext = list(H.greedy_linear_extensions_iterator())
            sage: SG_ext = list(H.supergreedy_linear_extensions_iterator())
            sage: [0, 2, 3, 1, 4] in G_ext
            True
            sage: [0, 2, 3, 1, 4] in SG_ext
            False

            sage: len(SG_ext)
            4

        TESTS::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: list(HasseDiagram({}).supergreedy_linear_extensions_iterator())
            [[]]
            sage: list(HasseDiagram({0: [], 1: []}).supergreedy_linear_extensions_iterator())
            [[0, 1], [1, 0]]
        """
        N = self.order()

        def supergreedy_rec(H, linext):
            k = len(linext)

            if k == N:
                yield linext

            else:
                S = []
                while not S:
                    if not k:  # Start from new minimal element
                        S = [x for x in self.sources() if x not in linext]
                    else:
                        S = [x for x in self.neighbor_out_iterator(linext[k - 1])
                             if x not in linext and
                             all(low in linext
                                 for low in self.neighbor_in_iterator(x))]
                        k -= 1

                for e in S:
                    yield from supergreedy_rec(H, linext + [e])

        return supergreedy_rec(self, [])

    def is_linear_extension(self, lin_ext=None) -> bool:
        r"""
        Test if an ordering is a linear extension.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,2],1:[3],2:[3],3:[]})
            sage: H.is_linear_extension(list(range(4)))
            True
            sage: H.is_linear_extension([3,2,1,0])
            False
        """
        if lin_ext is None or lin_ext == list(range(len(self))):
            return all(x < y for x, y in self.cover_relations_iterator())
        else:
            return all(lin_ext.index(x) < lin_ext.index(y)
                       for x, y in self.cover_relations_iterator())

    def cover_relations_iterator(self):
        r"""
        Iterate over cover relations.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[2,3], 1:[3,4], 2:[5], 3:[5], 4:[5]})
            sage: list(H.cover_relations_iterator())
            [(0, 2), (0, 3), (1, 3), (1, 4), (2, 5), (3, 5), (4, 5)]
        """
        yield from self.edge_iterator(labels=False)

    def cover_relations(self):
        r"""
        Return the list of cover relations.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[2,3], 1:[3,4], 2:[5], 3:[5], 4:[5]})
            sage: H.cover_relations()
            [(0, 2), (0, 3), (1, 3), (1, 4), (2, 5), (3, 5), (4, 5)]
        """
        return list(self.cover_relations_iterator())

    def is_lequal(self, i, j) -> bool:
        """
        Return ``True`` if i is less than or equal to j in the poset, and
        ``False`` otherwise.

        .. note::

            If the :meth:`lequal_matrix` has been computed, then this method is
            redefined to use the cached data (see :meth:`_alternate_is_lequal`).

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: x,y,z = 0, 1, 4
            sage: H.is_lequal(x,y)
            False
            sage: H.is_lequal(y,x)
            False
            sage: H.is_lequal(x,z)
            True
            sage: H.is_lequal(y,z)
            True
            sage: H.is_lequal(z,z)
            True
        """
        return i == j or (i < j and j in self.breadth_first_search(i))

    def is_less_than(self, x, y) -> bool:
        r"""
        Return ``True`` if ``x`` is less than but not equal to ``y`` in the
        poset, and ``False`` otherwise.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: x,y,z = 0, 1, 4
            sage: H.is_less_than(x,y)
            False
            sage: H.is_less_than(y,x)
            False
            sage: H.is_less_than(x,z)
            True
            sage: H.is_less_than(y,z)
            True
            sage: H.is_less_than(z,z)
            False
        """
        if x == y:
            return False
        return self.is_lequal(x, y)

    def is_gequal(self, x, y) -> bool:
        r"""
        Return ``True`` if ``x`` is greater than or equal to ``y``, and
        ``False`` otherwise.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: Q = HasseDiagram({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: x,y,z = 0,1,4
            sage: Q.is_gequal(x,y)
            False
            sage: Q.is_gequal(y,x)
            False
            sage: Q.is_gequal(x,z)
            False
            sage: Q.is_gequal(z,x)
            True
            sage: Q.is_gequal(z,y)
            True
            sage: Q.is_gequal(z,z)
            True
        """
        return self.is_lequal(y, x)

    def is_greater_than(self, x, y) -> bool:
        """
        Return ``True`` if ``x`` is greater than but not equal to
        ``y``, and ``False`` otherwise.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: Q = HasseDiagram({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: x,y,z = 0,1,4
            sage: Q.is_greater_than(x,y)
            False
            sage: Q.is_greater_than(y,x)
            False
            sage: Q.is_greater_than(x,z)
            False
            sage: Q.is_greater_than(z,x)
            True
            sage: Q.is_greater_than(z,y)
            True
            sage: Q.is_greater_than(z,z)
            False
        """
        return self.is_less_than(y, x)

    def minimal_elements(self):
        """
        Return a list of the minimal elements of the poset.

        EXAMPLES::

            sage: P = Poset({0:[3],1:[3],2:[3],3:[4],4:[]})
            sage: P(0) in P.minimal_elements()
            True
            sage: P(1) in P.minimal_elements()
            True
            sage: P(2) in P.minimal_elements()
            True
        """
        return self.sources()

    def maximal_elements(self):
        """
        Return a list of the maximal elements of the poset.

        EXAMPLES::

            sage: P = Poset({0:[3],1:[3],2:[3],3:[4],4:[]})
            sage: P.maximal_elements()
            [4]
        """
        return self.sinks()

    @cached_method
    def bottom(self):
        """
        Return the bottom element of the poset, if it exists.

        EXAMPLES::

            sage: P = Poset({0:[3],1:[3],2:[3],3:[4],4:[]})
            sage: P.bottom() is None
            True
            sage: Q = Poset({0:[1],1:[]})
            sage: Q.bottom()
            0
        """
        min_elms = self.minimal_elements()
        if len(min_elms) == 1:
            return min_elms[0]
        return None

    def has_bottom(self) -> bool:
        """
        Return ``True`` if the poset has a unique minimal element.

        EXAMPLES::

            sage: P = Poset({0:[3],1:[3],2:[3],3:[4],4:[]})
            sage: P.has_bottom()
            False
            sage: Q = Poset({0:[1],1:[]})
            sage: Q.has_bottom()
            True
        """
        return self.bottom() is not None

    def top(self):
        """
        Return the top element of the poset, if it exists.

        EXAMPLES::

            sage: P = Poset({0:[3],1:[3],2:[3],3:[4,5],4:[],5:[]})
            sage: P.top() is None
            True
            sage: Q = Poset({0:[1],1:[]})
            sage: Q.top()
            1
        """
        max_elms = self.maximal_elements()
        if len(max_elms) == 1:
            return max_elms[0]
        return None

    def has_top(self) -> bool:
        """
        Return ``True`` if the poset contains a unique maximal element, and
        ``False`` otherwise.

        EXAMPLES::

            sage: P = Poset({0:[3],1:[3],2:[3],3:[4,5],4:[],5:[]})
            sage: P.has_top()
            False
            sage: Q = Poset({0:[1],1:[]})
            sage: Q.has_top()
            True
        """
        return self.top() is not None

    def is_bounded(self) -> bool:
        """
        Return ``True`` if the poset contains a unique maximal element and a
        unique minimal element, and ``False`` otherwise.

        EXAMPLES::

            sage: P = Poset({0:[3],1:[3],2:[3],3:[4,5],4:[],5:[]})
            sage: P.is_bounded()
            False
            sage: Q = Poset({0:[1],1:[]})
            sage: Q.is_bounded()
            True
        """
        return self.has_top() and self.has_bottom()

    def is_chain(self) -> bool:
        """
        Return ``True`` if the poset is totally ordered, and ``False`` otherwise.

        EXAMPLES::

            sage: L = Poset({0:[1],1:[2],2:[3],3:[4]})
            sage: L.is_chain()
            True
            sage: V = Poset({0:[1,2]})
            sage: V.is_chain()
            False

        TESTS:

        Check :trac:`15330`::

            sage: p = Poset(DiGraph({0:[1],2:[1]}))
            sage: p.is_chain()
            False
        """
        if self.cardinality() == 0:
            return True
        return (self.num_edges() + 1 == self.num_verts() and  # tree
                all(d <= 1 for d in self.out_degree()) and
                all(d <= 1 for d in self.in_degree()))

    def is_antichain_of_poset(self, elms) -> bool:
        """
        Return ``True`` if ``elms`` is an antichain of the Hasse
        diagram and ``False`` otherwise.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0: [1, 2, 3], 1: [4], 2: [4], 3: [4]})
            sage: H.is_antichain_of_poset([1, 2, 3])
            True
            sage: H.is_antichain_of_poset([0, 2, 3])
            False
        """
        from itertools import combinations
        elms_sorted = sorted(set(elms))
        return not any(self.is_lequal(a, b) for a, b in
                       combinations(elms_sorted, 2))

    def dual(self):
        """
        Return a poset that is dual to the given poset.

        This means that it has the same elements but opposite order.
        The elements are renumbered to ensure that ``range(n)``
        is a linear extension.

        EXAMPLES::

            sage: P = posets.IntegerPartitions(4)                                       # optional - sage.combinat
            sage: H = P._hasse_diagram; H                                               # optional - sage.combinat
            Hasse diagram of a poset containing 5 elements
            sage: H.dual()                                                              # optional - sage.combinat
            Hasse diagram of a poset containing 5 elements

        TESTS::

            sage: H = posets.IntegerPartitions(4)._hasse_diagram                        # optional - sage.combinat
            sage: H.is_isomorphic( H.dual().dual() )                                    # optional - sage.combinat
            True
            sage: H.is_isomorphic( H.dual() )                                           # optional - sage.combinat
            False
        """
        H = self.reverse(immutable=False)
        H.relabel(perm=list(range(H.num_verts() - 1, -1, -1)), inplace=True)
        return HasseDiagram(H)

    def _precompute_intervals(self):
        """
        Precompute all intervals of the poset.

        This will significantly speed up computing congruences. On the
        other hand, it will cost much more memory. Currently this is
        a hidden feature. See the example below for how to use it.

        EXAMPLES::

            sage: B4 = posets.BooleanLattice(4)
            sage: B4.is_isoform()  # Slow                                               # optional - sage.combinat
            True
            sage: B4._hasse_diagram._precompute_intervals()
            sage: B4 = posets.BooleanLattice(4)
            sage: B4.is_isoform()  # Faster now                                         # optional - sage.combinat
            True
        """
        n = self.order()
        v_up = (frozenset(self.depth_first_search(v)) for v in range(n))
        v_down = [frozenset(self.depth_first_search(v, neighbors=self.neighbors_in))
                  for v in range(n)]
        self._intervals = [[sorted(up.intersection(down)) for down in v_down]
                           for up in v_up]

    def interval(self, x, y):
        r"""
        Return a list of the elements `z` of ``self`` such that
        `x \leq z \leq y`.

        The order is that induced by the ordering in ``self.linear_extension``.

        INPUT:

        -  ``x`` -- any element of the poset

        -  ``y`` -- any element of the poset

        .. NOTE::

            The method :meth:`_precompute_intervals()` creates a cache
            which is used if available, making the function very fast.

        .. SEEALSO:: :meth:`interval_iterator`

        EXAMPLES::

            sage: uc = [[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]]
            sage: dag = DiGraph(dict(zip(range(len(uc)),uc)))
            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram(dag)
            sage: I = set([2,5,6,4,7])
            sage: I == set(H.interval(2,7))
            True
        """
        try:
            # when the intervals have been precomputed
            return self._intervals[x][y]
        except AttributeError:
            return list(self.interval_iterator(x, y))

    def interval_iterator(self, x, y):
        r"""
        Return an iterator of the elements `z` of ``self`` such that
        `x \leq z \leq y`.

        INPUT:

        -  ``x`` -- any element of the poset

        -  ``y`` -- any element of the poset

        .. SEEALSO:: :meth:`interval`

        .. NOTE::

            This becomes much faster when first calling :meth:`_leq_storage`,
            which precomputes the principal upper ideals.

        EXAMPLES::

            sage: uc = [[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]]
            sage: dag = DiGraph(dict(zip(range(len(uc)),uc)))
            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram(dag)
            sage: I = set([2,5,6,4,7])
            sage: I == set(H.interval_iterator(2,7))
            True
        """
        for z in range(x, y + 1):
            if self.is_lequal(x, z) and self.is_lequal(z, y):
                yield z

    closed_interval = interval

    def open_interval(self, x, y):
        """
        Return a list of the elements `z` of ``self`` such that `x < z < y`.

        The order is that induced by the ordering in ``self.linear_extension``.

        EXAMPLES::

            sage: uc = [[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]]
            sage: dag = DiGraph(dict(zip(range(len(uc)),uc)))
            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram(dag)
            sage: set([5,6,4]) == set(H.open_interval(2,7))
            True
            sage: H.open_interval(7,2)
            []
        """
        ci = self.interval(x, y)
        if not ci:
            return []
        else:
            return ci[1:-1]

    def rank_function(self):
        r"""
        Return the (normalized) rank function of the poset,
        if it exists.

        A *rank function* of a poset `P` is a function `r`
        that maps elements of `P` to integers and satisfies:
        `r(x) = r(y) + 1` if `x` covers `y`. The function `r`
        is normalized such that its minimum value on every
        connected component of the Hasse diagram of `P` is
        `0`. This determines the function `r` uniquely (when
        it exists).

        OUTPUT:

        - a lambda function, if the poset admits a rank function
        - ``None``, if the poset does not admit a rank function

        EXAMPLES::

            sage: P = Poset([[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]])
            sage: P.rank_function() is not None
            True
            sage: P = Poset(([1,2,3,4],[[1,4],[2,3],[3,4]]), facade = True)
            sage: P.rank_function() is not None
            True
            sage: P = Poset(([1,2,3,4,5],[[1,2],[2,3],[3,4],[1,5],[5,4]]), facade = True)
            sage: P.rank_function() is not None
            False
            sage: P = Poset(([1,2,3,4,5,6,7,8],[[1,4],[2,3],[3,4],[5,7],[6,7]]), facade = True)
            sage: f = P.rank_function(); f is not None
            True
            sage: f(5)
            0
            sage: f(2)
            0

        TESTS::

            sage: P = Poset([[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]])
            sage: r = P.rank_function()
            sage: for u,v in P.cover_relations_iterator():
            ....:     if r(v) != r(u) + 1:
            ....:         print("Bug in rank_function!")

        ::

            sage: Q = Poset([[1,2],[4],[3],[4],[]])
            sage: Q.rank_function() is None
            True

        test for issue :trac:`14006`::

            sage: H = Poset()._hasse_diagram
            sage: s = dumps(H)
            sage: f = H.rank_function()
            sage: s = dumps(H)
        """
        if self._rank is None:
            return None
        # the rank function is just the getitem of the list
        return self._rank.__getitem__

    @lazy_attribute
    def _rank(self):
        r"""
        Build the rank function of the poset, if it exists, i.e.
        an array ``d`` where ``d[object] = self.rank_function()(object)``

        A *rank function* of a poset `P` is a function `r`
        that maps elements of `P` to integers and satisfies:
        `r(x) = r(y) + 1` if `x` covers `y`. The function `r`
        is normalized such that its minimum value on every
        connected component of the Hasse diagram of `P` is
        `0`. This determines the function `r` uniquely (when
        it exists).

        EXAMPLES::

            sage: H = Poset()._hasse_diagram
            sage: H._rank
            []
            sage: H = Poset([[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]])._hasse_diagram
            sage: H._rank
            [0, 1, 1, 2, 2, 1, 2, 3]
            sage: H = Poset(([1,2,3,4,5],[[1,2],[2,3],[3,4],[1,5],[5,4]]))._hasse_diagram
            sage: H._rank is None
            True
        """
        # rank[i] is the rank of point i. It is equal to None until the rank of
        # i is computed
        rank = [None] * self.order()
        not_found = set(self.vertex_iterator())
        while not_found:
            y = not_found.pop()
            rank[y] = 0  # We set some vertex to have rank 0
            component = set([y])
            queue = set([y])
            while queue:
                # look at the neighbors of y and set the ranks;
                # then look at the neighbors of the neighbors ...
                y = queue.pop()
                for x in self.neighbor_out_iterator(y):
                    if rank[x] is None:
                        rank[x] = rank[y] + 1
                        queue.add(x)
                        component.add(x)
                for x in self.neighbor_in_iterator(y):
                    if rank[x] is None:
                        rank[x] = rank[y] - 1
                        queue.add(x)
                        component.add(x)
                    elif rank[x] != rank[y] - 1:
                        return None
            # Normalize the ranks of vertices in the connected component
            # so that smallest is 0:
            m = min(rank[j] for j in component)
            for j in component:
                rank[j] -= m
            not_found.difference_update(component)
        # now, all ranks are set.
        return rank

    def rank(self, element=None):
        r"""
        Return the rank of ``element``, or the rank of the poset if
        ``element`` is ``None``. (The rank of a poset is the length of
        the longest chain of elements of the poset.)

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,3,2],1:[4],2:[4,5,6],3:[6],4:[7],5:[7],6:[7],7:[]})
            sage: H.rank(5)
            2
            sage: H.rank()
            3
            sage: Q = HasseDiagram({0:[1,2],1:[3],2:[],3:[]})
            sage: Q.rank()
            2
            sage: Q.rank(1)
            1
        """
        if element is None:
            return len(self.level_sets()) - 1
        else:
            return self.rank_function()(element)

    def is_ranked(self) -> bool:
        r"""
        Return ``True`` if the poset is ranked, and ``False`` otherwise.

        A poset is *ranked* if it admits a rank function. For more information
        about the rank function, see :meth:`~rank_function`
        and :meth:`~is_graded`.

        EXAMPLES::

            sage: P = Poset([[1],[2],[3],[4],[]])
            sage: P.is_ranked()
            True
            sage: Q = Poset([[1,5],[2,6],[3],[4],[],[6,3],[4]])
            sage: Q.is_ranked()
            False
        """
        return bool(self.rank_function())

    def covers(self, x, y):
        """
        Return ``True`` if y covers x and ``False`` otherwise.

        EXAMPLES::

            sage: Q = Poset([[1,5],[2,6],[3],[4],[],[6,3],[4]])
            sage: Q.covers(Q(1),Q(6))
            True
            sage: Q.covers(Q(1),Q(4))
            False
        """
        return self.has_edge(x, y)

    def upper_covers_iterator(self, element):
        r"""
        Return the list of elements that cover ``element``.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,3,2],1:[4],2:[4,5,6],3:[6],4:[7],5:[7],6:[7],7:[]})
            sage: list(H.upper_covers_iterator(0))
            [1, 2, 3]
            sage: list(H.upper_covers_iterator(7))
            []
        """
        yield from self.neighbor_out_iterator(element)

    def lower_covers_iterator(self, element):
        r"""
        Return the list of elements that are covered by ``element``.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,3,2],1:[4],2:[4,5,6],3:[6],4:[7],5:[7],6:[7],7:[]})
            sage: list(H.lower_covers_iterator(0))
            []
            sage: list(H.lower_covers_iterator(4))
            [1, 2]
        """
        yield from self.neighbor_in_iterator(element)

    def cardinality(self):
        r"""
        Return the number of elements in the poset.

        EXAMPLES::

            sage: Poset([[1,2,3],[4],[4],[4],[]]).cardinality()
            5

        TESTS:

        For a time, this function was named ``size()``, which
        would override the same-named method of the underlying
        digraph. :trac:`8735` renamed this method to ``cardinality()``
        with a deprecation warning. :trac:`11214` removed the warning
        since code for graphs was raising the warning inadvertently.
        This tests that ``size()`` for a Hasse diagram returns the
        number of edges in the digraph. ::

            sage: L = posets.BooleanLattice(5)
            sage: H = L.hasse_diagram()
            sage: H.size()
            80
            sage: H.size() == H.num_edges()
            True
        """
        return self.order()

    def moebius_function(self, i, j):  # dumb algorithm
        r"""
        Return the value of the Möbius function of the poset
        on the elements ``i`` and ``j``.

        EXAMPLES::

            sage: P = Poset([[1,2,3],[4],[4],[4],[]])
            sage: H = P._hasse_diagram
            sage: H.moebius_function(0,4)
            2
            sage: for u,v in P.cover_relations_iterator():
            ....:     if P.moebius_function(u,v) != -1:
            ....:         print("Bug in moebius_function!")
        """
        try:
            return self._moebius_function_values[(i, j)]
        except AttributeError:
            self._moebius_function_values = {}
            return self.moebius_function(i, j)
        except KeyError:
            if i == j:
                self._moebius_function_values[(i, j)] = 1
            elif i > j:
                self._moebius_function_values[(i, j)] = 0
            else:
                ci = self.interval(i, j)
                if not ci:
                    self._moebius_function_values[(i, j)] = 0
                else:
                    self._moebius_function_values[(i, j)] = -sum(self.moebius_function(i, k) for k in ci[:-1])
        return self._moebius_function_values[(i, j)]

    def bottom_moebius_function(self, j):
        r"""
        Return the value of the Möbius function of the poset
        on the elements ``zero`` and ``j``, where ``zero`` is
        ``self.bottom()``, the unique minimal element of the poset.

        EXAMPLES::

            sage: P = Poset({0: [1,2]})
            sage: hasse = P._hasse_diagram
            sage: hasse.bottom_moebius_function(1)
            -1
            sage: hasse.bottom_moebius_function(2)
            -1
            sage: P = Poset({0: [1,3], 1:[2], 2:[4], 3:[4]})
            sage: hasse = P._hasse_diagram
            sage: for i in range(5):
            ....:   print(hasse.bottom_moebius_function(i))
            1
            -1
            0
            -1
            1

        TESTS::

            sage: P = Poset({0:[2], 1:[2]})
            sage: hasse = P._hasse_diagram
            sage: hasse.bottom_moebius_function(1)
            Traceback (most recent call last):
            ...
            ValueError: the poset does not have a bottom element
        """
        zero = self.bottom()
        if zero is None:
            raise ValueError("the poset does not have a bottom element")
        # if the value has already been computed, either by self.moebius_function
        # or by self.bottom_moebius_function, then just use the cached value.
        try:
            return self._moebius_function_values[(zero, j)]
        # if the dict has not been initialized, do that and try again
        except AttributeError:
            self._moebius_function_values = {}
            return self.bottom_moebius_function(j)
        # if mu(zero, j) has not already been computed, we'll get a key error.
        except KeyError:
            if zero == j:
                self._moebius_function_values[(zero, j)] = 1
            # since zero is the minimal element, we can ignore the case that zero > j,
            # and move on to computing the interval, which is exactly the order ideal.
            else:
                # do the depth_first_search over order_ideal, because we don't care
                # about sorting the elements of the order ideal.
                ci = self._backend.depth_first_search(j, reverse=True)
                next(ci)  # throw out the first element, which is j
                self._moebius_function_values[(zero, j)] = -sum(self.bottom_moebius_function(k) for k in ci)
        return self._moebius_function_values[(zero, j)]

    def moebius_function_matrix(self, algorithm='cython'):
        r"""
        Return the matrix of the Möbius function of this poset.

        This returns the matrix over `\ZZ` whose ``(x, y)`` entry
        is the value of the Möbius function of ``self`` evaluated on
        ``x`` and ``y``, and redefines :meth:`moebius_function` to use it.

        INPUT:

        - ``algorithm`` -- optional, ``'recursive'``, ``'matrix'``
          or ``'cython'`` (default)

        This uses either the recursive formula, a generic matrix inversion
        or a specific matrix inversion coded in Cython.

        OUTPUT:

        a dense matrix for the algorithm ``cython``, a sparse matrix otherwise

        .. NOTE::

            The result is cached in :meth:`_moebius_function_matrix`.

        .. SEEALSO:: :meth:`lequal_matrix`, :meth:`coxeter_transformation`

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,3,2],1:[4],2:[4,5,6],3:[6],4:[7],5:[7],6:[7],7:[]})
            sage: H.moebius_function_matrix()                                           # optional - sage.libs.flint sage.modules
            [ 1 -1 -1 -1  1  0  1  0]
            [ 0  1  0  0 -1  0  0  0]
            [ 0  0  1  0 -1 -1 -1  2]
            [ 0  0  0  1  0  0 -1  0]
            [ 0  0  0  0  1  0  0 -1]
            [ 0  0  0  0  0  1  0 -1]
            [ 0  0  0  0  0  0  1 -1]
            [ 0  0  0  0  0  0  0  1]

        TESTS::

            sage: # needs sage.modules sage.libs.flint
            sage: H.moebius_function_matrix().is_immutable()
            True
            sage: hasattr(H,'_moebius_function_matrix')
            True
            sage: H.moebius_function == H._moebius_function_from_matrix
            True
            sage: H = posets.TamariLattice(3)._hasse_diagram
            sage: M = H.moebius_function_matrix('matrix'); M
            [ 1 -1 -1  0  1]
            [ 0  1  0  0 -1]
            [ 0  0  1 -1  0]
            [ 0  0  0  1 -1]
            [ 0  0  0  0  1]
            sage: _ = H.__dict__.pop('_moebius_function_matrix')
            sage: H.moebius_function_matrix('cython') == M
            True
            sage: _ = H.__dict__.pop('_moebius_function_matrix')
            sage: H.moebius_function_matrix('recursive') == M
            True
            sage: _ = H.__dict__.pop('_moebius_function_matrix')
            sage: H.moebius_function_matrix('banana')
            Traceback (most recent call last):
            ...
            ValueError: unknown algorithm
        """
        if not hasattr(self, '_moebius_function_matrix'):
            if algorithm == 'recursive':
                n = self.cardinality()
                gt = self._leq_storage
                greater_than = [sorted(gt[i]) for i in range(n)]
                m = {}
                for i in range(n - 1, -1, -1):
                    m[(i, i)] = ZZ.one()
                    available = []
                    for k in greater_than[i]:
                        if k != i:
                            available.append(k)
                            m[(i, k)] = -ZZ.sum(m[(j, k)]
                                                for j in available
                                                if k in greater_than[j])
                M = matrix(ZZ, n, n, m, sparse=True)
            elif algorithm == "matrix":
                M = self.lequal_matrix().inverse_of_unit()
            elif algorithm == "cython":
                M = moebius_matrix_fast(self._leq_storage)
            else:
                raise ValueError("unknown algorithm")
            self._moebius_function_matrix = M
            self._moebius_function_matrix.set_immutable()
            self.moebius_function = self._moebius_function_from_matrix
        return self._moebius_function_matrix

    # Redefine self.moebius_function
    def _moebius_function_from_matrix(self, i, j):
        r"""
        Return the value of the Möbius function of the poset
        on the elements ``i`` and ``j``.

        EXAMPLES::

            sage: P = Poset([[1,2,3],[4],[4],[4],[]])
            sage: H = P._hasse_diagram
            sage: H.moebius_function(0,4) # indirect doctest
            2
            sage: for u,v in P.cover_relations_iterator():
            ....:     if P.moebius_function(u,v) != -1:
            ....:         print("Bug in moebius_function!")

        This uses ``self._moebius_function_matrix``, as computed by
        :meth:`moebius_function_matrix`.
        """
        return self._moebius_function_matrix[i, j]

    @cached_method
    def coxeter_transformation(self, algorithm='cython'):
        r"""
        Return the matrix of the Auslander-Reiten translation acting on
        the Grothendieck group of the derived category of modules on the
        poset, in the basis of simple modules.

        INPUT:

        - ``algorithm`` -- optional, ``'cython'`` (default) or ``'matrix'``

        This uses either a specific matrix code in Cython, or generic matrices.

        .. SEEALSO:: :meth:`lequal_matrix`, :meth:`moebius_function_matrix`

        EXAMPLES::

            sage: # needs sage.libs.flint sage.modules
            sage: P = posets.PentagonPoset()._hasse_diagram
            sage: M = P.coxeter_transformation(); M
            [ 0  0  0  0 -1]
            [ 0  0  0  1 -1]
            [ 0  1  0  0 -1]
            [-1  1  1  0 -1]
            [-1  1  0  1 -1]
            sage: P.__dict__['coxeter_transformation'].clear_cache()
            sage: P.coxeter_transformation(algorithm="matrix") == M
            True

        TESTS::

            sage: # needs sage.libs.flint sage.modules
            sage: P = posets.PentagonPoset()._hasse_diagram
            sage: M = P.coxeter_transformation()
            sage: M**8 == 1
            True
            sage: P.__dict__['coxeter_transformation'].clear_cache()
            sage: P.coxeter_transformation(algorithm="banana")
            Traceback (most recent call last):
            ...
            ValueError: unknown algorithm
        """
        if algorithm == 'matrix':
            return - self.lequal_matrix() * self.moebius_function_matrix().transpose()
        elif algorithm == 'cython':
            return coxeter_matrix_fast(self._leq_storage)
        else:
            raise ValueError("unknown algorithm")

    def order_filter(self, elements):
        r"""
        Return the order filter generated by a list of elements.

        `I` is an order filter if, for any `x` in `I` and `y` such that
        `y \ge x`, then `y` is in `I`.

        EXAMPLES::

            sage: H = posets.BooleanLattice(4)._hasse_diagram
            sage: H.order_filter([3,8])
            [3, 7, 8, 9, 10, 11, 12, 13, 14, 15]
        """
        return sorted(self.depth_first_search(elements))

    def principal_order_filter(self, i):
        """
        Return the order filter generated by ``i``.

        EXAMPLES::

            sage: H = posets.BooleanLattice(4)._hasse_diagram
            sage: H.principal_order_filter(2)
            [2, 3, 6, 7, 10, 11, 14, 15]
        """
        return self.order_filter([i])

    def order_ideal(self, elements):
        r"""
        Return the order ideal generated by a list of elements.

        `I` is an order ideal if, for any `x` in `I` and `y` such that
        `y \le x`, then `y` is in `I`.

        EXAMPLES::

            sage: H = posets.BooleanLattice(4)._hasse_diagram
            sage: H.order_ideal([7,10])
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 10]
        """
        return sorted(self.depth_first_search(elements,
                                              neighbors=self.neighbors_in))

    def order_ideal_cardinality(self, elements):
        r"""
        Return the cardinality of the order ideal generated by ``elements``.

        `I` is an order ideal if, for any `x` in `I` and `y` such that
        `y \le x`, then `y` is in `I`.

        EXAMPLES::

            sage: H = posets.BooleanLattice(4)._hasse_diagram
            sage: H.order_ideal_cardinality([7,10])
            10
        """
        seen = set()
        q = deque(elements)
        size = 0
        while q:
            v = q.popleft()
            if v in seen:
                continue
            size += 1
            seen.add(v)
            q.extend(self.neighbors_in(v))

        return ZZ(size)

    def principal_order_ideal(self, i):
        """
        Return the order ideal generated by `i`.

        EXAMPLES::

            sage: H = posets.BooleanLattice(4)._hasse_diagram
            sage: H.principal_order_ideal(6)
            [0, 2, 4, 6]
        """
        return self.order_ideal([i])

    @lazy_attribute
    def _leq_storage(self):
        """
        Store the comparison relation as a list of Python sets.

        The `i`-th item in the list is the set of elements greater than `i`.

        EXAMPLES::

            sage: H = posets.DiamondPoset(7)._hasse_diagram
            sage: H._leq_storage
            [{0, 1, 2, 3, 4, 5, 6}, {1, 6}, {2, 6}, {3, 6}, {4, 6}, {5, 6}, {6}]
        """
        n = self.order()
        greater_than = [set([i]) for i in range(n)]
        for i in range(n - 1, -1, -1):
            gt = greater_than[i]
            for j in self.neighbor_out_iterator(i):
                gt = gt.union(greater_than[j])
            greater_than[i] = gt

        # Redefine self.is_lequal
        self.is_lequal = self._alternate_is_lequal

        return greater_than

    @lazy_attribute
    def _leq_matrix_boolean(self):
        r"""
        Compute a boolean matrix whose ``(i,j)`` entry is 1 if ``i``
        is less than ``j`` in the poset, and 0 otherwise; and
        redefines ``__lt__`` to use this matrix.

        .. SEEALSO:: :meth:`_leq_matrix`

        EXAMPLES::

            sage: P = Poset([[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]])
            sage: H = P._hasse_diagram
            sage: M = H._leq_matrix_boolean; M                                          # optional - sage.modules sage.rings.finite_rings
            [1 1 1 1 1 1 1 1]
            [0 1 0 1 0 0 0 1]
            [0 0 1 1 1 0 1 1]
            [0 0 0 1 0 0 0 1]
            [0 0 0 0 1 0 0 1]
            [0 0 0 0 0 1 1 1]
            [0 0 0 0 0 0 1 1]
            [0 0 0 0 0 0 0 1]
            sage: M.base_ring()                                                         # optional - sage.modules sage.rings.finite_rings
            Finite Field of size 2
        """
        n = self.order()
        R = GF(2)
        one = R.one()
        greater_than = self._leq_storage
        D = {(i, j): one for i in range(n) for j in greater_than[i]}
        M = matrix(R, n, n, D, sparse=True)
        M.set_immutable()
        return M

    @lazy_attribute
    def _leq_matrix(self):
        r"""
        Compute an integer matrix whose ``(i,j)`` entry is 1 if ``i``
        is less than ``j`` in the poset, and 0 otherwise.

        .. SEEALSO:: :meth:`_leq_matrix_boolean`

        EXAMPLES::

            sage: P = Poset([[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]])
            sage: H = P._hasse_diagram
            sage: M = H._leq_matrix; M                                                  # optional - sage.modules
            [1 1 1 1 1 1 1 1]
            [0 1 0 1 0 0 0 1]
            [0 0 1 1 1 0 1 1]
            [0 0 0 1 0 0 0 1]
            [0 0 0 0 1 0 0 1]
            [0 0 0 0 0 1 1 1]
            [0 0 0 0 0 0 1 1]
            [0 0 0 0 0 0 0 1]
            sage: M.base_ring()                                                         # optional - sage.modules
            Integer Ring
        """
        n = self.order()
        greater_than = self._leq_storage
        D = {(i, j): 1 for i in range(n) for j in greater_than[i]}
        return matrix(ZZ, n, n, D, sparse=True, immutable=True)

    def lequal_matrix(self, boolean=False):
        r"""
        Return a matrix whose ``(i,j)`` entry is 1 if ``i`` is less
        than ``j`` in the poset, and 0 otherwise; and redefines
        ``__lt__`` to use the boolean version of this matrix.

        INPUT:

        - ``boolean`` -- optional flag (default ``False``) telling whether to
          return a matrix with coefficients in `\GF(2)` or in `\ZZ`

        .. SEEALSO::

            :meth:`moebius_function_matrix`, :meth:`coxeter_transformation`

        EXAMPLES::

            sage: P = Poset([[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]])
            sage: H = P._hasse_diagram
            sage: M = H.lequal_matrix(); M                                              # optional - sage.modules
            [1 1 1 1 1 1 1 1]
            [0 1 0 1 0 0 0 1]
            [0 0 1 1 1 0 1 1]
            [0 0 0 1 0 0 0 1]
            [0 0 0 0 1 0 0 1]
            [0 0 0 0 0 1 1 1]
            [0 0 0 0 0 0 1 1]
            [0 0 0 0 0 0 0 1]
            sage: M.base_ring()                                                         # optional - sage.modules
            Integer Ring

            sage: P = posets.DiamondPoset(6)
            sage: H = P._hasse_diagram
            sage: M = H.lequal_matrix(boolean=True)                                     # optional - sage.modules sage.rings.finite_rings
            sage: M.base_ring()                                                         # optional - sage.modules sage.rings.finite_rings
            Finite Field of size 2

        TESTS::

            sage: H.lequal_matrix().is_immutable()                                      # optional - sage.modules
            True
        """
        if boolean:
            return self._leq_matrix_boolean
        else:
            return self._leq_matrix

    def _alternate_is_lequal(self, i, j):
        r"""
        Return ``True`` if ``i`` is less than or equal to ``j`` in
        ``self``, and ``False`` otherwise.

        .. NOTE::

            If the :meth:`lequal_matrix` has been computed, then
            :meth:`is_lequal` is redefined to use the cached data.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: H.lequal_matrix()                                                     # optional - sage.modules
            [1 0 1 1 1]
            [0 1 1 1 1]
            [0 0 1 1 1]
            [0 0 0 1 1]
            [0 0 0 0 1]
            sage: x,y,z = 0, 1, 4
            sage: H._alternate_is_lequal(x,y)
            False
            sage: H._alternate_is_lequal(y,x)
            False
            sage: H._alternate_is_lequal(x,z)
            True
            sage: H._alternate_is_lequal(y,z)
            True
            sage: H._alternate_is_lequal(z,z)
            True
        """
        return j in self._leq_storage[i]

    def prime_elements(self):
        r"""
        Return the join-prime and meet-prime elements of the bounded poset.

        An element `x` of a poset `P` is join-prime if the subposet
        induced by `\{y \in P \mid y \not\ge x\}` has a top element.
        Meet-prime is defined dually.

        .. NOTE::

            The poset is expected to be bounded, and this is *not* checked.

        OUTPUT:

        A pair `(j, m)` where `j` is a list of join-prime elements
        and `m` is a list of meet-prime elements.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0: [1, 2], 1: [3], 2: [4], 3: [4]})
            sage: H.prime_elements()
            ([1, 2], [2, 3])
        """
        n = self.order()
        join_primes = []
        meet_primes = []

        def add_elements(e):
            upset = frozenset(self.depth_first_search(e))
            # The complement of the upper set of a join-prime must have
            # a top element. Maximal elements of the complement are those
            # covered by only elements in the upper set. If there is only
            # one maximal element, it is a meet-prime and 'e' is a
            # join-prime.
            meet_prime = None
            for u in upset:
                for m in self.neighbor_in_iterator(u):
                    if (m not in upset and
                        all(u_ in upset for u_ in
                            self.neighbor_out_iterator(m))):
                        if meet_prime is not None:
                            return
                        meet_prime = m
            join_primes.append(e)
            meet_primes.append(meet_prime)

        for e in range(n):
            # Join-primes are join-irreducibles, only check those.
            if self.in_degree(e) == 1:
                add_elements(e)

        return join_primes, meet_primes

    @lazy_attribute
    def _meet(self):
        r"""
        Return the matrix of meets of ``self``.

        The ``(x,y)``-entry of this matrix is the meet of ``x`` and
        ``y`` in ``self`` if the meet exists; and `-1` otherwise.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,3,2],1:[4],2:[4,5,6],3:[6],4:[7],5:[7],6:[7],7:[]})
            sage: H._meet                                                               # optional - sage.modules
            [0 0 0 0 0 0 0 0]
            [0 1 0 0 1 0 0 1]
            [0 0 2 0 2 2 2 2]
            [0 0 0 3 0 0 3 3]
            [0 1 2 0 4 2 2 4]
            [0 0 2 0 2 5 2 5]
            [0 0 2 3 2 2 6 6]
            [0 1 2 3 4 5 6 7]

            sage: H = HasseDiagram({0:[2,3],1:[2,3]})
            sage: H._meet                                                               # optional - sage.modules
            [ 0 -1  0  0]
            [-1  1  1  1]
            [ 0  1  2 -1]
            [ 0  1 -1  3]

            sage: H = HasseDiagram({0:[1,2],1:[3,4],2:[3,4]})
            sage: H._meet                                                               # optional - sage.modules
            [ 0  0  0  0  0]
            [ 0  1  0  1  1]
            [ 0  0  2  2  2]
            [ 0  1  2  3 -1]
            [ 0  1  2 -1  4]

        TESTS::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: L = LatticePoset({0:[1,2,3],1:[4],2:[4],3:[4]})                       # optional - sage.modules
            sage: P = L.dual()                                                          # optional - sage.modules
            sage: P.meet(2,3)                                                           # optional - sage.modules
            4
        """
        self._meet_semilattice_failure = ()
        n = self.cardinality()
        if n == 0:
            return matrix(0)
        meet = [[-1 for x in range(n)] for x in range(n)]
        lc = [self.neighbors_in(x) for x in range(n)]  # Lc = lower covers

        for x in range(n):
            meet[x][x] = x
            for y in range(x):
                T = [meet[y][z] for z in lc[x] if meet[y][z] != -1]
                if not T:
                    q = -1
                else:
                    q = max(T)
                    for z in T:
                        if meet[z][q] != z:
                            q = -1
                            break
                meet[x][y] = q
                meet[y][x] = q
                if q == -1:
                    self._meet_semilattice_failure += ((x, y),)
        return matrix(ZZ, meet)

    def meet_matrix(self):
        r"""
        Return the matrix of meets of ``self``, when ``self`` is a
        meet-semilattice; raise an error otherwise.

        The ``(x,y)``-entry of this matrix is the meet of ``x`` and
        ``y`` in ``self``.

        This algorithm is modelled after the algorithm of Freese-Jezek-Nation
        (p217). It can also be found on page 140 of [Gec81]_.

        .. NOTE::

            If ``self`` is a meet-semilattice, then the return of this method
            is the same as :meth:`_meet`. Once the matrix has been computed,
            it is stored in :meth:`_meet`. Delete this attribute if you want to
            recompute the matrix.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,3,2],1:[4],2:[4,5,6],3:[6],4:[7],5:[7],6:[7],7:[]})
            sage: H.meet_matrix()                                                       # optional - sage.modules
            [0 0 0 0 0 0 0 0]
            [0 1 0 0 1 0 0 1]
            [0 0 2 0 2 2 2 2]
            [0 0 0 3 0 0 3 3]
            [0 1 2 0 4 2 2 4]
            [0 0 2 0 2 5 2 5]
            [0 0 2 3 2 2 6 6]
            [0 1 2 3 4 5 6 7]

        REFERENCE:

        .. [Gec81] Fundamentals of Computation Theory
          Gecseg, F.
          Proceedings of the 1981 International Fct-Conference
          Szeged, Hungaria, August 24-28, vol 117
          Springer-Verlag, 1981

        TESTS::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[2,3],1:[2,3]})
            sage: H.meet_matrix()
            Traceback (most recent call last):
            ...
            ValueError: not a meet-semilattice: no bottom element

            sage: H = HasseDiagram({0:[1,2],1:[3,4],2:[3,4]})
            sage: H.meet_matrix()                                                       # optional - sage.modules
            Traceback (most recent call last):
            ...
            LatticeError: no meet for ...
        """
        n = self.cardinality()
        if (n != 0) and (not self.has_bottom()):
            raise ValueError("not a meet-semilattice: no bottom element")
        # call the attribute to build the matrix and _meet_semilattice_failure
        mt = self._meet
        if self._meet_semilattice_failure:
            x, y = self._meet_semilattice_failure[0]
            raise LatticeError('meet', x, y)
        return mt

    def is_meet_semilattice(self) -> bool:
        r"""
        Return ``True`` if ``self`` has a meet operation, and
        ``False`` otherwise.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,3,2],1:[4],2:[4,5,6],3:[6],4:[7],5:[7],6:[7],7:[]})
            sage: H.is_meet_semilattice()                                               # optional - sage.modules
            True

            sage: H = HasseDiagram({0:[1,2],1:[3],2:[3],3:[]})
            sage: H.is_meet_semilattice()                                               # optional - sage.modules
            True

            sage: H = HasseDiagram({0:[2,3],1:[2,3]})
            sage: H.is_meet_semilattice()                                               # optional - sage.modules
            False

            sage: H = HasseDiagram({0:[1,2],1:[3,4],2:[3,4]})
            sage: H.is_meet_semilattice()                                               # optional - sage.modules
            False
        """
        try:
            self.meet_matrix()
        except ValueError:
            return False
        else:
            return True

    @lazy_attribute
    def _join(self):
        r"""
        Computes a matrix whose ``(x,y)``-entry is the join of ``x``
        and ``y`` in ``self`` if the join exists; and `-1` otherwise.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,3,2],1:[4],2:[4,5,6],3:[6],4:[7],5:[7],6:[7],7:[]})
            sage: H._join                                                               # optional - sage.modules
            [0 1 2 3 4 5 6 7]
            [1 1 4 7 4 7 7 7]
            [2 4 2 6 4 5 6 7]
            [3 7 6 3 7 7 6 7]
            [4 4 4 7 4 7 7 7]
            [5 7 5 7 7 5 7 7]
            [6 7 6 6 7 7 6 7]
            [7 7 7 7 7 7 7 7]

            sage: H = HasseDiagram({0:[2,3],1:[2,3]})
            sage: H._join                                                               # optional - sage.modules
            [ 0 -1  2  3]
            [-1  1  2  3]
            [ 2  2  2 -1]
            [ 3  3 -1  3]

            sage: H = HasseDiagram({0:[2,3],1:[2,3],2:[4],3:[4]})
            sage: H._join                                                               # optional - sage.modules
            [ 0 -1  2  3  4]
            [-1  1  2  3  4]
            [ 2  2  2  4  4]
            [ 3  3  4  3  4]
            [ 4  4  4  4  4]

        TESTS::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: L = LatticePoset({0:[1,2,3],1:[4],2:[4],3:[4]})                       # optional - sage.modules
            sage: P = L.dual()                                                          # optional - sage.modules
            sage: P.join(2,3)                                                           # optional - sage.modules
            0
        """
        self._join_semilattice_failure = ()
        n = self.cardinality()
        if n == 0:
            return matrix(0)
        join = [[-1 for x in range(n)] for x in range(n)]
        uc = [self.neighbors_out(x) for x in range(n)]  # uc = upper covers

        for x in range(n - 1, -1, -1):
            join[x][x] = x
            for y in range(n - 1, x, -1):
                T = [join[y][z] for z in uc[x] if join[y][z] != -1]
                if not T:
                    q = -1
                else:
                    q = min(T)
                    for z in T:
                        if join[z][q] != z:
                            q = -1
                            break
                join[x][y] = q
                join[y][x] = q
                if q == -1:
                    self._join_semilattice_failure += ((x, y),)

        return matrix(ZZ, join)

    def join_matrix(self):
        r"""
        Return the matrix of joins of ``self``, when ``self`` is a
        join-semilattice; raise an error otherwise.

        The ``(x,y)``-entry of this matrix is the join of ``x`` and
        ``y`` in ``self``.

        This algorithm is modelled after the algorithm of Freese-Jezek-Nation
        (p217). It can also be found on page 140 of [Gec81]_.

        .. NOTE::

            If ``self`` is a join-semilattice, then the return of this method
            is the same as :meth:`_join`. Once the matrix has been computed,
            it is stored in :meth:`_join`. Delete this attribute if you want
            to recompute the matrix.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,3,2],1:[4],2:[4,5,6],3:[6],4:[7],5:[7],6:[7],7:[]})
            sage: H.join_matrix()                                                       # optional - sage.modules
            [0 1 2 3 4 5 6 7]
            [1 1 4 7 4 7 7 7]
            [2 4 2 6 4 5 6 7]
            [3 7 6 3 7 7 6 7]
            [4 4 4 7 4 7 7 7]
            [5 7 5 7 7 5 7 7]
            [6 7 6 6 7 7 6 7]
            [7 7 7 7 7 7 7 7]

        TESTS::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[2,3],1:[2,3]})
            sage: H.join_matrix()
            Traceback (most recent call last):
            ...
            ValueError: not a join-semilattice: no top element

            sage: H = HasseDiagram({0:[2,3],1:[2,3],2:[4],3:[4]})
            sage: H.join_matrix()                                                       # optional - sage.modules
            Traceback (most recent call last):
            ...
            LatticeError: no join for ...
        """
        n = self.cardinality()
        if (n != 0) and (not self.has_top()):
            raise ValueError("not a join-semilattice: no top element")
        # call the attribute to build the matrix and _join_semilattice_failure
        jn = self._join
        if self._join_semilattice_failure:
            x, y = self._join_semilattice_failure[0]
            raise LatticeError('join', x, y)
        return jn

    def is_join_semilattice(self) -> bool:
        r"""
        Return ``True`` if ``self`` has a join operation, and
        ``False`` otherwise.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,3,2],1:[4],2:[4,5,6],3:[6],4:[7],5:[7],6:[7],7:[]})
            sage: H.is_join_semilattice()                                               # optional - sage.modules
            True
            sage: H = HasseDiagram({0:[2,3],1:[2,3]})
            sage: H.is_join_semilattice()                                               # optional - sage.modules
            False
            sage: H = HasseDiagram({0:[2,3],1:[2,3],2:[4],3:[4]})
            sage: H.is_join_semilattice()                                               # optional - sage.modules
            False
        """
        try:
            self.join_matrix()
        except ValueError:
            return False
        else:
            return True

    def find_nonsemidistributive_elements(self, meet_or_join):
        r"""
        Check if the lattice is semidistributive or not.

        INPUT:

        - ``meet_or_join`` -- string ``'meet'`` or ``'join'``
          to decide if to check for join-semidistributivity or
          meet-semidistributivity

        OUTPUT:

        - ``None`` if the lattice is semidistributive OR
        - tuple ``(u, e, x, y)`` such that
          `u = e \vee x = e \vee y` but `u \neq e \vee (x \wedge y)`
          if ``meet_or_join=='join'`` and
          `u = e \wedge x = e \wedge y` but `u \neq e \wedge (x \vee y)`
          if ``meet_or_join=='meet'``

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1, 2], 1:[3, 4], 2:[4, 5], 3:[6],
            ....:                   4:[6], 5:[6]})
            sage: H.find_nonsemidistributive_elements('join') is None                   # optional - sage.modules
            False
            sage: H.find_nonsemidistributive_elements('meet') is None                   # optional - sage.modules
            True
        """
        if meet_or_join == 'join':
            M1 = self.join_matrix()
            M2 = self.meet_matrix()
        elif meet_or_join == 'meet':
            M1 = self.meet_matrix()
            M2 = self.join_matrix()
        else:
            raise ValueError("meet_or_join must be 'join' or 'meet'")

        n = self.order()

        for e in range(n):
            for x in range(n):
                u = M1[e, x]
                for y in range(x):
                    if u == M1[e, y]:
                        if u != M1[e, M2[x, y]]:
                            return (u, e, x, y)

        return None

    def vertical_decomposition(self, return_list=False):
        """
        Return vertical decomposition of the lattice.

        This is the backend function for vertical decomposition
        functions of lattices.

        The property of being vertically decomposable is defined for lattices.
        This is *not* checked, and the function works with any bounded poset.

        INPUT:

        - ``return_list``, a boolean. If ``False`` (the default), return
          an element that is not the top neither the bottom element of the
          lattice, but is comparable to all elements of the lattice, if
          the lattice is vertically decomposable and ``None`` otherwise.
          If ``True``, return list of decomposition elements.

        EXAMPLES::

            sage: H = posets.BooleanLattice(4)._hasse_diagram
            sage: H.vertical_decomposition() is None
            True
            sage: P = Poset( ([1,2,3,6,12,18,36], attrcall("divides")) )
            sage: P._hasse_diagram.vertical_decomposition()
            3
            sage: P._hasse_diagram.vertical_decomposition(return_list=True)
            [3]
        """
        n = self.cardinality()
        if n < 3:
            if return_list:
                return []
            else:
                return None
        result = []  # Never take the bottom element to list.
        m = 0
        for i in range(n - 1):
            for j in self.outgoing_edge_iterator(i):
                m = max(m, j[1])
            if m == i + 1:
                if not return_list:
                    if m < n - 1:
                        return m
                    else:
                        return None
                result.append(m)
        result.pop()  # Remove the top element.
        return result

    def is_complemented(self) -> int | None:
        """
        Return an element of the lattice that has no complement.

        If the lattice is complemented, return ``None``.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram

            sage: H = HasseDiagram({0:[1, 2], 1:[3], 2:[3], 3:[4]})
            sage: H.is_complemented()                                                   # optional - sage.modules
            1

            sage: H = HasseDiagram({0:[1, 2, 3], 1:[4], 2:[4], 3:[4]})
            sage: H.is_complemented() is None                                           # optional - sage.modules
            True
        """
        mt = self.meet_matrix()
        jn = self.join_matrix()
        top = self.cardinality() - 1
        has_complement = [False] * top

        for i in range(1, top):
            if has_complement[i]:
                continue
            for j in range(top, 0, -1):
                if jn[i, j] == top and mt[i, j] == 0:
                    has_complement[j] = True
                    break
            else:
                return i

        return None

    def pseudocomplement(self, element):
        """
        Return the pseudocomplement of ``element``, if it exists.

        The pseudocomplement is the greatest element whose
        meet with given element is the bottom element. It may
        not exist, and then the function returns ``None``.

        INPUT:

        - ``element`` -- an element of the lattice.

        OUTPUT:

        An element of the Hasse diagram, i.e. an integer, or
        ``None`` if the pseudocomplement does not exist.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0: [1, 2], 1: [3], 2: [4], 3: [4]})
            sage: H.pseudocomplement(2)                                                 # optional - sage.modules
            3

            sage: H = HasseDiagram({0: [1, 2, 3], 1: [4], 2: [4], 3: [4]})
            sage: H.pseudocomplement(2) is None                                         # optional - sage.modules
            True
        """
        e = self.order() - 1
        mt = self.meet_matrix()
        while mt[e, element] != 0:
            e -= 1
        e1 = e
        while e1 > 0:
            if mt[e1, element] == 0 and not self.is_lequal(e1, e):
                return None
            e1 -= 1
        return e

    def orthocomplementations_iterator(self):
        r"""
        Return an iterator over orthocomplementations of the lattice.

        OUTPUT:

        An iterator that gives plain list of integers.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,2], 1:[3,4], 3:[5], 4:[5], 2:[6,7],
            ....:                   6:[8], 7:[8], 5:[9], 8:[9]})
            sage: list(H.orthocomplementations_iterator())                              # optional - sage.groups
            [[9, 8, 5, 6, 7, 2, 3, 4, 1, 0], [9, 8, 5, 7, 6, 2, 4, 3, 1, 0]]

        ALGORITHM:

        As ``DiamondPoset(2*n+2)`` has `(2n)!/(n!2^n)` different
        orthocomplementations, the complexity of listing all of
        them is necessarily `O(n!)`.

        An orthocomplemented lattice is self-dual, so that for example
        orthocomplement of an atom is a coatom. This function
        basically just computes list of possible orthocomplementations
        for every element (i.e. they must be complements and "duals"),
        and then tries to fit them all.

        TESTS:

        Special and corner cases::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram()  # Empty
            sage: list(H.orthocomplementations_iterator())                              # optional - sage.groups
            [[]]
            sage: H = HasseDiagram({0:[]})  # One element
            sage: list(H.orthocomplementations_iterator())                              # optional - sage.groups
            [[0]]
            sage: H = HasseDiagram({0:[1]})  # Two elements
            sage: list(H.orthocomplementations_iterator())                              # optional - sage.groups
            [[1, 0]]

        Trivial cases: odd number of elements, not self-dual, not complemented::

            sage: H = posets.DiamondPoset(5)._hasse_diagram
            sage: list(H.orthocomplementations_iterator())                              # optional - sage.groups
            []
            sage: H = posets.ChainPoset(4)._hasse_diagram
            sage: list(H.orthocomplementations_iterator())                              # optional - sage.groups
            []
            sage: H = HasseDiagram( ([[0, 1], [0, 2], [0, 3], [1, 4], [1, 8], [4, 6], [4, 7], [6, 9], [7, 9], [2, 5], [3, 5], [5, 8], [8, 9]]) )
            sage: list(H.orthocomplementations_iterator())                              # optional - sage.groups
            []
            sage: H = HasseDiagram({0:[1, 2, 3], 1: [4], 2:[4], 3: [5], 4:[5]})
            sage: list(H.orthocomplementations_iterator())                              # optional - sage.groups
            []

        Complemented, self-dual and even number of elements, but
        not orthocomplemented::

            sage: H = HasseDiagram( ([[0, 1], [1, 2], [2, 3], [0, 4], [4, 5], [0, 6], [3, 7], [5, 7], [6, 7]]) )
            sage: list(H.orthocomplementations_iterator())                              # optional - sage.groups
            []

        Unique orthocomplementations; second is not uniquely complemented,
        but has only one orthocomplementation::

            sage: H = posets.BooleanLattice(4)._hasse_diagram  # Uniquely complemented
            sage: len(list(H.orthocomplementations_iterator()))                         # optional - sage.groups
            1
            sage: H = HasseDiagram({0:[1, 2], 1:[3], 2:[4], 3:[5], 4:[5]})
            sage: len([_ for _ in H.orthocomplementations_iterator()])                  # optional - sage.groups
            1

        "Lengthening diamond" must keep the number of orthocomplementations::

            sage: H = HasseDiagram( ([[0, 1], [0, 2], [0, 3], [0, 4], [1, 5], [2, 5], [3, 5], [4, 5]]) )
            sage: n = len([_ for _ in H.orthocomplementations_iterator()]); n           # optional - sage.groups
            3
            sage: H = HasseDiagram('M]??O?@??C??OA???OA??@?A??C?A??O??')
            sage: len([_ for _ in H.orthocomplementations_iterator()]) == n             # optional - sage.groups
            True

        This lattice has an unique "possible orthocomplement" for every
        element, but they can not be fit together; orthocomplement pairs
        would be 0-11, 1-7, 2-4, 3-10, 5-9 and 6-8, and then orthocomplements
        for chain 0-1-6-11 would be 11-7-8-0, which is not a chain::

            sage: H = HasseDiagram('KTGG_?AAC?O?o?@?@?E?@?@??')
            sage: list(H.orthocomplementations_iterator())                              # optional - sage.groups
            []
        """
        n = self.order()

        # Special cases first
        if n == 0:
            yield []
            return
        if n == 1:
            yield [0]
            return
        if n % 2:
            return

        dual_isomorphism = self.is_isomorphic(self.reverse(), certificate=True)[1]
        if dual_isomorphism is None:  # i.e. if the lattice is not self-dual.
            return

        # We compute possible orthocomplements, i.e. elements
        # with "dual position" and complement to each other.

        orbits = self.automorphism_group(return_group=False, orbits=True)

        orbit_number = [None] * n
        for ind, orbit in enumerate(orbits):
            for e in orbit:
                orbit_number[e] = ind

        comps = [None] * n
        mt = self.meet_matrix()
        jn = self.join_matrix()
        for e in range(n):
            # Fix following after issue #20727
            comps[e] = [x for x in range(n) if
                        mt[e, x] == 0 and jn[e, x] == n - 1 and
                        x in orbits[orbit_number[dual_isomorphism[e]]]]

        # Fitting is done by this recursive function:
        def recursive_fit(orthocomplements, unbinded):
            if not unbinded:
                yield orthocomplements
            else:
                next_to_fit = unbinded[0]
                possible_values = [x for x in comps[next_to_fit]
                                   if x not in orthocomplements]
                for x in self.lower_covers_iterator(next_to_fit):
                    if orthocomplements[x] is not None:
                        possible_values = [y for y in possible_values if self.has_edge(y, orthocomplements[x])]
                for x in self.upper_covers_iterator(next_to_fit):
                    if orthocomplements[x] is not None:
                        possible_values = [y for y in possible_values if self.has_edge(orthocomplements[x], y)]

                for e in possible_values:

                    new_binded = orthocomplements[:]
                    new_binded[next_to_fit] = e
                    new_binded[e] = next_to_fit

                    new_unbinded = unbinded[1:]  # Remove next_to_fit
                    new_unbinded.remove(e)

                    yield from recursive_fit(new_binded, new_unbinded)

        start = [None] * n
        # A little optimization
        for e in range(n):
            if not comps[e]:  # Not any possible orthocomplement
                return
            if len(comps[e]) == 1:  # Do not re-fit this every time
                e_ = comps[e][0]
                # Every element might have one possible orthocomplement,
                # but so that they don't fit together. Must check that.
                for lc in self.lower_covers_iterator(e):
                    if start[lc] is not None:
                        if not self.has_edge(e_, start[lc]):
                            return
                if start[e_] is None:
                    start[e] = e_
                    start[e_] = e
        start_unbinded = [e for e in range(n) if start[e] is None]

        yield from recursive_fit(start, start_unbinded)

    def find_nonsemimodular_pair(self, upper):
        """
        Return pair of elements showing the lattice is not modular.

        INPUT:

        - ``upper``, a Boolean -- if ``True``, test whether the lattice is
          upper semimodular; otherwise test whether the lattice is
          lower semimodular.

        OUTPUT:

        ``None``, if the lattice is semimodular. Pair `(a, b)` violating
        semimodularity otherwise.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1, 2], 1:[3, 4], 2:[4, 5], 3:[6], 4:[6], 5:[6]})
            sage: H.find_nonsemimodular_pair(upper=True) is None
            True
            sage: H.find_nonsemimodular_pair(upper=False)
            (5, 3)

            sage: H_ = HasseDiagram(H.reverse().relabel(lambda x: 6-x, inplace=False))
            sage: H_.find_nonsemimodular_pair(upper=True)
            (3, 1)
            sage: H_.find_nonsemimodular_pair(upper=False) is None
            True
        """
        if upper:
            neighbors = self.neighbors_out
        else:
            neighbors = self.neighbors_in

        n = self.order()
        for e in range(n):
            covers = neighbors(e)
            covers_len = len(covers)
            if covers_len < 2:
                continue
            for a_i in range(covers_len):
                a = covers[a_i]
                covers_a = neighbors(a)
                for b_i in range(a_i):
                    b = covers[b_i]
                    if not any(j in covers_a for j in neighbors(b)):
                        return (a, b)
        return None

    def antichains_iterator(self):
        r"""
        Return an iterator over the antichains of the poset.

        .. note::

            The algorithm is based on Freese-Jezek-Nation p. 226.
            It does a depth first search through the set of all
            antichains organized in a prefix tree.

        EXAMPLES::

            sage: # needs sage.modules
            sage: P = posets.PentagonPoset()
            sage: H = P._hasse_diagram
            sage: H.antichains_iterator()
            <generator object ...antichains_iterator at ...>
            sage: list(H.antichains_iterator())
            [[], [4], [3], [2], [1], [1, 3], [1, 2], [0]]
            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,2],1:[4],2:[3],3:[4]})
            sage: list(H.antichains_iterator())
            [[], [4], [3], [2], [1], [1, 3], [1, 2], [0]]
            sage: H = HasseDiagram({0:[],1:[],2:[]})
            sage: list(H.antichains_iterator())
            [[], [2], [1], [1, 2], [0], [0, 2], [0, 1], [0, 1, 2]]
            sage: H = HasseDiagram({0:[1],1:[2],2:[3],3:[4]})
            sage: list(H.antichains_iterator())
            [[], [4], [3], [2], [1], [0]]

        TESTS::

            sage: H = Poset()._hasse_diagram
            sage: list(H.antichains_iterator())                                         # optional - sage.modules
            [[]]
        """
        # NOTE: Ordering of antichains as a prefix tree is crucial for
        # congruences_iterator() to work. Change it, if you change this.

        # Complexity note:
        # antichains_queues never grows longer than self.cardinality().
        # Indeed, if a appears before b in antichains_queues, then
        # the largest element of a is strictly smaller than that of b.
        antichains_queues = [([], list(range(self.cardinality() - 1, -1, -1)))]
        leq = self.lequal_matrix()
        while antichains_queues:
            (antichain, queue) = antichains_queues.pop()
            # Invariant:
            #  - the elements of antichain are independent
            #  - the elements of queue are independent from those of antichain
            yield antichain
            while queue:
                x = queue.pop()
                new_antichain = antichain + [x]
                new_queue = [t for t in queue if not (leq[t, x] or leq[x, t])]
                antichains_queues.append((new_antichain, new_queue))

    def are_incomparable(self, i, j):
        """
        Return whether ``i`` and ``j`` are incomparable in the poset

        INPUT:

        - ``i``, ``j`` -- vertices of this Hasse diagram

        EXAMPLES::

            sage: # needs sage.modules
            sage: P = posets.PentagonPoset()
            sage: H = P._hasse_diagram
            sage: H.are_incomparable(1,2)
            True
            sage: V = H.vertices(sort=True)
            sage: [ (i,j) for i in V for j in V if H.are_incomparable(i,j)]
            [(1, 2), (1, 3), (2, 1), (3, 1)]
        """
        if i == j:
            return False
        if i > j:
            i, j = j, i
        mat = self._leq_matrix_boolean
        return not mat[i, j]

    def are_comparable(self, i, j):
        """
        Return whether ``i`` and ``j`` are comparable in the poset

        INPUT:

         - ``i``, ``j`` -- vertices of this Hasse diagram

        EXAMPLES::

            sage: # needs sage.modules
            sage: P = posets.PentagonPoset()
            sage: H = P._hasse_diagram
            sage: H.are_comparable(1,2)
            False
            sage: V = H.vertices(sort=True)
            sage: [ (i,j) for i in V for j in V if H.are_comparable(i,j)]
            [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (1, 0), (1, 1), (1, 4),
             (2, 0), (2, 2), (2, 3), (2, 4), (3, 0), (3, 2), (3, 3), (3, 4),
             (4, 0), (4, 1), (4, 2), (4, 3), (4, 4)]
        """
        if i == j:
            return True
        if i > j:
            i, j = j, i
        mat = self._leq_matrix_boolean
        return bool(mat[i, j])

    def antichains(self, element_class=list):
        """
        Return all antichains of ``self``, organized as a prefix tree

        INPUT:

        - ``element_class`` -- (default:list) an iterable type

        EXAMPLES::

            sage: # needs sage.modules
            sage: P = posets.PentagonPoset()
            sage: H = P._hasse_diagram
            sage: A = H.antichains()
            sage: list(A)
            [[], [0], [1], [1, 2], [1, 3], [2], [3], [4]]
            sage: A.cardinality()
            8
            sage: [1,3] in A
            True
            sage: [1,4] in A
            False

        TESTS::

            sage: # needs sage.modules
            sage: TestSuite(A).run()
            sage: A = Poset()._hasse_diagram.antichains()
            sage: list(A)
            [[]]
            sage: TestSuite(A).run()
        """
        from sage.combinat.subsets_pairwise import PairwiseCompatibleSubsets
        return PairwiseCompatibleSubsets(self.vertices(sort=True),
                                         self.are_incomparable,
                                         element_class=element_class)

    def chains(self, element_class=list, exclude=None, conversion=None):
        """
        Return all chains of ``self``, organized as a prefix tree.

        INPUT:

        - ``element_class`` -- (optional, default: ``list``) an iterable type

        - ``exclude`` -- elements of the poset to be excluded
          (optional, default: ``None``)

        - ``conversion`` -- (optional, default: ``None``) used to pass
           the list of elements of the poset in their fixed order

        OUTPUT:

        The enumerated set (with a forest structure given by prefix
        ordering) consisting of all chains of ``self``, each of
        which is given as an ``element_class``.

        If ``conversion`` is given, then the chains are converted
        to chain of elements of this list.

        EXAMPLES::

            sage: # needs sage.modules
            sage: P = posets.PentagonPoset()
            sage: H = P._hasse_diagram
            sage: A = H.chains()
            sage: list(A)
            [[], [0], [0, 1], [0, 1, 4], [0, 2], [0, 2, 3], [0, 2, 3, 4], [0, 2, 4],
             [0, 3], [0, 3, 4], [0, 4], [1], [1, 4], [2], [2, 3], [2, 3, 4], [2, 4],
             [3], [3, 4], [4]]
            sage: A.cardinality()
            20
            sage: [1,3] in A
            False
            sage: [1,4] in A
            True

        One can exclude some vertices::

            sage: # needs sage.modules
            sage: list(H.chains(exclude=[4, 3]))
            [[], [0], [0, 1], [0, 2], [1], [2]]

        The ``element_class`` keyword determines how the chains are
        being returned::

            sage: P = Poset({1: [2, 3], 2: [4]})
            sage: list(P._hasse_diagram.chains(element_class=tuple))
            [(), (0,), (0, 1), (0, 1, 2), (0, 2), (0, 3), (1,), (1, 2), (2,), (3,)]
            sage: list(P._hasse_diagram.chains())
            [[], [0], [0, 1], [0, 1, 2], [0, 2], [0, 3], [1], [1, 2], [2], [3]]

        (Note that taking the Hasse diagram has renamed the vertices.) ::

            sage: list(P._hasse_diagram.chains(element_class=tuple, exclude=[0]))
            [(), (1,), (1, 2), (2,), (3,)]

        .. SEEALSO:: :meth:`antichains`
        """
        return IncreasingChains(self._leq_storage, element_class, exclude, conversion)

    def is_linear_interval(self, t_min, t_max) -> bool:
        """
        Return whether the interval ``[t_min, t_max]`` is linear.

        This means that this interval is a total order.

        EXAMPLES::
            sage: # needs sage.modules
            sage: P = posets.PentagonPoset()
            sage: H = P._hasse_diagram
            sage: H.is_linear_interval(0, 4)
            False
            sage: H.is_linear_interval(0, 3)
            True
            sage: H.is_linear_interval(1, 3)
            False
            sage: H.is_linear_interval(1, 1)
            True

        TESTS::

            sage: P = posets.TamariLattice(3)
            sage: H = P._hasse_diagram
            sage: D = H._leq_storage
            sage: a, b = H.bottom(), H.top()
            sage: H.is_linear_interval(a, b)
            False
            sage: H.is_linear_interval(a, a)
            True
        """
        if '_leq_storage' in self.__dict__:
            if not self.is_lequal(t_min, t_max):  # very quick check
                return False
            t = t_max
            while t != t_min:
                found = False
                for u in self.neighbor_in_iterator(t):
                    if self.is_lequal(t_min, u):
                        if not found:
                            found = True
                            t = u
                        else:
                            return False
            return True

        # fall back to default implementation
        it = self.all_paths_iterator([t_min], [t_max],
                                     simple=True, trivial=True)
        try:
            next(it)
        except StopIteration:  # not comparable
            return False
        try:
            next(it)
        except StopIteration:  # one path
            return True
        return False

    def diamonds(self):
        r"""
        Return the list of diamonds of ``self``.

        A diamond is the following subgraph of the Hasse diagram::

                    z
                   / \
                  x   y
                   \ /
                    w

        Thus each edge represents a cover relation in the Hasse diagram.
        We represent his as the tuple `(w, x, y, z)`.

        OUTPUT:

        A tuple with

        - a list of all diamonds in the Hasse Diagram,
        - a boolean checking that every `w,x,y` that form a ``V``, there is a
          unique element `z`, which completes the diamond.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0: [1,2], 1: [3], 2: [3], 3: []})
            sage: H.diamonds()
            ([(0, 1, 2, 3)], True)

            sage: P = posets.YoungDiagramPoset(Partition([3, 2, 2]))                    # optional - sage.combinat
            sage: H = P._hasse_diagram                                                  # optional - sage.combinat
            sage: H.diamonds()                                                          # optional - sage.combinat
            ([(0, 1, 3, 4), (3, 4, 5, 6)], False)
        """
        diamonds = []
        all_diamonds_completed = True
        for w in self.vertices(sort=True):
            covers = self.neighbors_out(w)
            for i, x in enumerate(covers):
                for y in covers[i + 1:]:
                    zs = self.common_upper_covers([x, y])
                    if len(zs) != 1:
                        all_diamonds_completed = False
                    for z in zs:
                        diamonds.append((w, x, y, z))
        return (diamonds, all_diamonds_completed)

    def common_upper_covers(self, vertices):
        r"""
        Return the list of all common upper covers of ``vertices``.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0: [1,2], 1: [3], 2: [3], 3: []})
            sage: H.common_upper_covers([1, 2])
            [3]

            sage: from sage.combinat.posets.poset_examples import Posets
            sage: H = Posets.YoungDiagramPoset(Partition([3, 2, 2]))._hasse_diagram     # optional - sage.combinat
            sage: H.common_upper_covers([4, 5])                                         # optional - sage.combinat
            [6]
        """
        covers = set(self.neighbors_out(vertices.pop()))
        for v in vertices:
            covers = covers.intersection(self.neighbors_out(v))
        return list(covers)

    def common_lower_covers(self, vertices):
        r"""
        Return the list of all common lower covers of ``vertices``.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0: [1,2], 1: [3], 2: [3], 3: []})
            sage: H.common_lower_covers([1, 2])
            [0]

            sage: from sage.combinat.posets.poset_examples import Posets
            sage: H = Posets.YoungDiagramPoset(Partition([3, 2, 2]))._hasse_diagram     # optional - sage.combinat
            sage: H.common_lower_covers([4, 5])                                         # optional - sage.combinat
            [3]
        """
        covers = set(self.neighbors_in(vertices.pop()))
        for v in vertices:
            covers = covers.intersection(self.neighbors_in(v))
        return list(covers)

    def _trivial_nonregular_congruence(self):
        """
        Return a pair of elements giving "trivial" non-regular congruence.

        This returns a pair `a, b` such that `b` covers only `a` and
        `a` is covered by only `b`, and either `a` has one lower cover
        or `b` has one upper cover.  If no such pair exists, return
        ``None``.

        This pair gives a trivial non-regular congruence.

        The Hasse diagram is expected to be bounded.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0: [1, 2], 1: [4], 2: [3], 3: [4]})
            sage: H._trivial_nonregular_congruence()
            (2, 3)

        TESTS::

            sage: H = HasseDiagram({0: [1, 2], 1: [3], 2: [3]})
            sage: H._trivial_nonregular_congruence() is None
            True
            sage: H = HasseDiagram({0: [1, 2], 1: [3], 2: [3], 3: [4]})
            sage: H._trivial_nonregular_congruence()
            (3, 4)
            sage: H = HasseDiagram({0: [1], 1: [2, 3], 2: [4], 3: [4]})
            sage: H._trivial_nonregular_congruence()
            (0, 1)
            sage: H = HasseDiagram({0: [1]})
            sage: H._trivial_nonregular_congruence() is None
            True
        """
        n = self.order()
        if n == 2:
            return None
        if self.out_degree(0) == 1:
            return (0, 1)
        if self.in_degree(n - 1) == 1:
            return (n - 2, n - 1)
        for v in range(1, n - 1):
            if self.in_degree(v) == 1 and self.out_degree(v) == 1:
                v_ = next(self.neighbor_out_iterator(v))
                if self.in_degree(v_) == 1 and self.out_degree(v_) == 1:
                    return (v, v_)
        return None

    def sublattices_iterator(self, elms, min_e):
        """
        Return an iterator over sublattices of the Hasse diagram.

        INPUT:

        - ``elms`` -- elements already in sublattice; use set() at start
        - ``min_e`` -- smallest new element to add for new sublattices

        OUTPUT:

        List of sublattices as sets of integers.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0: [1, 2], 1:[3], 2:[3]})
            sage: it = H.sublattices_iterator(set(), 0); it
            <generator object ...sublattices_iterator at ...>
            sage: next(it)                                                              # optional - sage.modules
            set()
            sage: next(it)                                                              # optional - sage.modules
            {0}
        """
        yield elms
        mt = self.meet_matrix()
        jn = self.join_matrix()
        for e in range(min_e, self.cardinality()):
            if e in elms:
                continue
            current_set = set(elms)
            gens = set([e])
            while gens:
                g = gens.pop()
                if g < e and g not in elms:
                    break
                if g in current_set:
                    continue
                for x in current_set:
                    gens.add(mt[x, g])
                    gens.add(jn[x, g])
                current_set.add(g)
            else:
                yield from self.sublattices_iterator(current_set, e + 1)

    def maximal_sublattices(self):
        """
        Return maximal sublattices of the lattice.

        EXAMPLES::

            sage: L = posets.PentagonPoset()                                            # optional - sage.modules
            sage: ms = L._hasse_diagram.maximal_sublattices()                           # optional - sage.modules
            sage: sorted(ms, key=sorted)                                                # optional - sage.modules
            [{0, 1, 2, 4}, {0, 1, 3, 4}, {0, 2, 3, 4}]
        """
        jn = self.join_matrix()
        mt = self.meet_matrix()

        def sublattice(elms, e):
            """
            Helper function to get sublattice generated by list
            of elements.
            """
            gens_remaining = set([e])
            current_set = set(elms)

            while gens_remaining:
                g = gens_remaining.pop()
                if g in current_set:
                    continue
                for x in current_set:
                    gens_remaining.add(jn[x, g])
                    gens_remaining.add(mt[x, g])
                current_set.add(g)

            return current_set

        N = self.cardinality()
        elms = [0]
        sublats = [set([0])]
        result = []
        skip = -1

        while True:
            # First try to append an element
            found_element_to_append = False
            e = elms[-1]
            while e != skip:
                e += 1
                if e == N:
                    maybe_found = sublats[-1]
                    if not any(maybe_found.issubset(x) for x in result):
                        result.append(sublats[-1])
                    break
                if e in sublats[-1]:
                    continue
                # Let's try to add 'e' and see what happens.
                sl = sublattice(sublats[-1], e)
                if len(sl) < N:
                    # Skip this, if it generated a back-reference.
                    new_elms = sl.difference(sublats[-1])
                    if not any(x < e for x in new_elms):
                        found_element_to_append = True
                        break
                # Now sl is whole lattice, so we continue and try
                # appending another element.

            if found_element_to_append:
                elms.append(e)
                sublats.append(sl)
                continue

            # Can not append. Try to increment last element.
            e = elms.pop()
            sublats.pop()

            last_element_increment = True
            while True:
                e += 1
                if e == N:
                    last_element_increment = False
                    break
                if e in sublats[-1]:
                    continue
                sl = sublattice(sublats[-1], e)
                if len(sl) == N:
                    continue

                new_elms = sl.difference(set(sublats[-1]))
                if any(x < e for x in new_elms):
                    continue

                elms.append(e)
                sublats.append(sl)
                break

            if not last_element_increment:
                # Can not append nor increment. "Backtracking".
                skip = elms[-1]
                if skip == 0:
                    break

        # Special case to handle at last.
        if len(self.neighbors_out(0)) == 1:
            result.append(set(range(1, N)))

        return result

    def frattini_sublattice(self):
        """
        Return the list of elements of the Frattini sublattice of the lattice.

        EXAMPLES::

            sage: H = posets.PentagonPoset()._hasse_diagram                             # optional - sage.modules
            sage: H.frattini_sublattice()                                               # optional - sage.modules
            [0, 4]
        """
        # Just a direct computation, no optimization at all.
        n = self.cardinality()
        if n == 0 or n == 2:
            return []
        if n == 1:
            return [0]
        max_sublats = self.maximal_sublattices()
        return [e for e in range(self.cardinality()) if
                all(e in ms for ms in max_sublats)]

    def kappa_dual(self, a):
        r"""
        Return the minimum element smaller than the element covering
        ``a`` but not smaller than ``a``.

        Define `\kappa^*(a)` as the minimum element of
        `(\downarrow a_*) \setminus (\downarrow a)`, where `a_*` is the element
        covering `a`. It is always a join-irreducible element, if it exists.

        .. NOTE::

            Element ``a`` is expected to be meet-irreducible, and
            this is *not* checked.

        INPUT:

        - ``a`` -- a join-irreducible element of the lattice

        OUTPUT:

        The element `\kappa^*(a)` or ``None`` if there
        is not a unique smallest element with given constraints.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0: [1, 2], 1: [3, 4], 2: [4, 5], 3: [6], 4: [6], 5: [6]})
            sage: H.kappa_dual(3)
            2
            sage: H.kappa_dual(4) is None
            True

        TESTS::

            sage: H = HasseDiagram({0: [1]})
            sage: H.kappa_dual(0)
            1
        """
        uc = next(self.neighbor_out_iterator(a))
        if self.in_degree(uc) == 1:
            return uc
        lt_a = set(self.depth_first_search(a, neighbors=self.neighbors_in))
        tmp = list(self.depth_first_search(uc, neighbors=lambda v: [v_ for v_ in self.neighbor_in_iterator(v) if v_ not in lt_a]))
        result = None
        for e in tmp:
            if all(x not in tmp for x in self.neighbor_in_iterator(e)):
                if result:
                    return None
                result = e
        return result

    def skeleton(self):
        """
        Return the skeleton of the lattice.

        The lattice is expected to be pseudocomplemented and non-empty.

        The skeleton of the lattice is the subposet induced by
        those elements that are the pseudocomplement to at least one
        element.

        OUTPUT:

        List of elements such that the subposet induced by them is
        the skeleton of the lattice.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0: [1, 2], 1: [3, 4], 2: [4],
            ....:                   3: [5], 4: [5]})
            sage: H.skeleton()                                                          # optional - sage.modules
            [5, 2, 0, 3]
        """
        p_atoms = []
        for atom in self.neighbor_out_iterator(0):
            p_atom = self.pseudocomplement(atom)
            if p_atom is None:
                raise ValueError("lattice is not pseudocomplemented")
            p_atoms.append(p_atom)
        n = len(p_atoms)
        mt = self.meet_matrix()
        pos = [0] * n
        meets = [self.order() - 1] * n
        result = [self.order() - 1]
        i = 0

        while i >= 0:
            new_meet = mt[meets[i - 1], p_atoms[pos[i]]]
            result.append(new_meet)
            if pos[i] == n - 1:
                i -= 1
                pos[i] += 1
            else:
                meets[i] = new_meet
                pos[i + 1] = pos[i] + 1
                i += 1

        return result

    def is_convex_subset(self, S) -> bool:
        r"""
        Return ``True`` if `S` is a convex subset of the poset,
        and ``False`` otherwise.

        A subset `S` is *convex* in the poset if `b \in S` whenever
        `a, c \in S` and `a \le b \le c`.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: B3 = HasseDiagram({0: [1, 2, 4], 1: [3, 5], 2: [3, 6],
            ....:                    3: [7], 4: [5, 6], 5: [7], 6: [7]})
            sage: B3.is_convex_subset([1, 3, 5, 4])  # Also connected
            True
            sage: B3.is_convex_subset([1, 3, 4])  # Not connected
            True

            sage: B3.is_convex_subset([0, 1, 2, 3, 6])  # No, 0 < 4 < 6
            False
            sage: B3.is_convex_subset([0, 1, 2, 7])  # No, 1 < 3 < 7.
            False

        TESTS::

            sage: B3.is_convex_subset([])
            True
            sage: B3.is_convex_subset([6])
            True
        """
        if not S:  # S is empty set
            return True
        s_max = max(S)
        ok = set()  # Already checked elements not less than any element is S.

        for a in S:
            for b in self.neighbor_out_iterator(a):
                if b >= s_max or b in S:
                    continue
                # Now b not in S, b > a and a in S.

                def neighbors(v_):
                    return [v for v in self.neighbor_out_iterator(v_)
                            if v <= s_max and v not in ok]
                for c in self.depth_first_search(b, neighbors=neighbors):
                    if c in S:  # Now c in S, b not in S, a in S, a < b < c.
                        return False
                    ok.add(c)  # Do not re-check this for being our b.

        return True

    def neutral_elements(self):
        """
        Return the list of neutral elements of the lattice.

        An element `a` in a lattice is neutral if the sublattice
        generated by `a`, `x` and `y` is distributive for every
        `x`, `y` in the lattice.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0: [1, 2], 1: [4], 2: [3], 3: [4, 5],
            ....:                   4: [6], 5: [6]})
            sage: sorted(H.neutral_elements())                                          # optional - sage.modules
            [0, 4, 6]

        ALGORITHM:

        Basically we just check the distributivity against all element
        pairs `x, y` to see if element `a` is neutral or not.

        If we found that `a, x, y` is not a distributive triple, we add
        all three to list of non-neutral elements. If we found `a` to
        be neutral, we add it to list of neutral elements. When testing
        we skip already found neutral elements, as they can't be our `x`
        or `y`.

        We skip `a, x, y` as trivial if it is a chain. We do that by
        letting `x` to be a non-comparable to `a`; `y` can be any element.

        We first try to found `x` and `y` from elements not yet tested,
        so that we could get three birds with one stone.

        And last, the top and bottom elements are always neutral and
        need not be tested.
        """
        n = self.order()
        if n < 5:
            return set(range(n))

        todo = set(range(1, n - 1))
        neutrals = set([0, n - 1])
        notneutrals = set()
        all_elements = set(range(n))

        mt = self.meet_matrix()
        jn = self.join_matrix()

        def is_neutral(a):
            noncomp = all_elements.difference(self.depth_first_search(a))
            noncomp.difference_update(self.depth_first_search(a, neighbors=self.neighbors_in))

            for x in noncomp.intersection(todo):
                meet_ax = mt[a, x]
                join_ax = jn[a, x]
                for y in todo:
                    if (mt[mt[join_ax, jn[a, y]], jn[x, y]] !=
                            jn[jn[meet_ax, mt[a, y]], mt[x, y]]):
                        notneutrals.add(x)
                        notneutrals.add(y)
                        return False
                for y in notneutrals:
                    if (mt[mt[join_ax, jn[a, y]], jn[x, y]] !=
                            jn[jn[meet_ax, mt[a, y]], mt[x, y]]):
                        notneutrals.add(x)
                        return False
            for x in noncomp.difference(todo):
                meet_ax = mt[a, x]
                join_ax = jn[a, x]
                for y in todo:
                    if (mt[mt[join_ax, jn[a, y]], jn[x, y]] !=
                            jn[jn[meet_ax, mt[a, y]], mt[x, y]]):
                        notneutrals.add(y)
                        return False
                for y in notneutrals:
                    if (mt[mt[join_ax, jn[a, y]], jn[x, y]] !=
                            jn[jn[meet_ax, mt[a, y]], mt[x, y]]):
                        return False
            return True

        while todo:
            e = todo.pop()
            if is_neutral(e):
                neutrals.add(e)
            else:
                notneutrals.add(e)

        return neutrals

    def kappa(self, a):
        r"""
        Return the maximum element greater than the element covered
        by ``a`` but not greater than ``a``.

        Define `\kappa(a)` as the maximum element of
        `(\uparrow a_*) \setminus (\uparrow a)`, where `a_*` is the element
        covered by `a`. It is always a meet-irreducible element, if it exists.

        .. NOTE::

            Element ``a`` is expected to be join-irreducible, and
            this is *not* checked.

        INPUT:

        - ``a`` -- a join-irreducible element of the lattice

        OUTPUT:

        The element `\kappa(a)` or ``None`` if there
        is not a unique greatest element with given constraints.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0: [1, 2, 3], 1: [4], 2: [4, 5], 3: [5], 4: [6], 5: [6]})
            sage: H.kappa(1)
            5
            sage: H.kappa(2) is None
            True

        TESTS::

            sage: H = HasseDiagram({0: [1]})
            sage: H.kappa(1)
            0
        """
        lc = next(self.neighbor_in_iterator(a))
        if self.out_degree(lc) == 1:
            return lc
        gt_a = set(self.depth_first_search(a))
        tmp = list(self.depth_first_search(lc, neighbors=lambda v: [v_ for v_ in self.neighbor_out_iterator(v) if v_ not in gt_a]))
        result = None
        for e in tmp:
            if all(x not in tmp for x in self.neighbor_out_iterator(e)):
                if result:
                    return None
                result = e
        return result

    def atoms_of_congruence_lattice(self):
        r"""
        Return atoms of the congruence lattice.

        In other words, return "minimal non-trivial" congruences:
        A congruence is minimal if the only finer (as a partition
        of set of elements) congruence is the trivial congruence
        where every block contains only one element.

        .. SEEALSO:: :meth:`congruence`

        OUTPUT:

        List of congruences, every congruence as
        :class:`sage.combinat.set_partition.SetPartition`

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: N5 = HasseDiagram({0: [1, 2], 1: [4], 2: [3], 3:[4]})
            sage: N5.atoms_of_congruence_lattice()                                      # optional - sage.combinat
            [{{0}, {1}, {2, 3}, {4}}]
            sage: Hex = HasseDiagram({0: [1, 2], 1: [3], 2: [4], 3: [5], 4: [5]})
            sage: Hex.atoms_of_congruence_lattice()                                     # optional - sage.combinat
            [{{0}, {1}, {2, 4}, {3}, {5}}, {{0}, {1, 3}, {2}, {4}, {5}}]

        ALGORITHM:

        Every atom is a join-irreducible. Every join-irreducible of
        `\mathrm{Con}(L)` is a principal congruence generated by a
        meet-irreducible element and the only element covering it (and also
        by a join-irreducible element and the only element covered by it).
        Hence we check those principal congruences to find the minimal ones.
        """
        # Note: A lattice L if subdirectly reducible (i.e. is a sublattice
        # of a Cartesian product of two smaller lattices) iff Con(L) has
        # at least two atoms. That's were this is used for.

        from sage.combinat.set_partition import SetPartitions

        # Get smaller set, meet- or join-irreducibles
        join_irreducibles = [v for v in self if self.in_degree(v) == 1]
        meet_irreducibles = [v for v in self if self.out_degree(v) == 1]
        if len(join_irreducibles) < len(meet_irreducibles):
            irr = [(v, next(self.neighbor_in_iterator(v))) for v in join_irreducibles]
        else:
            irr = [(next(self.neighbor_out_iterator(v)), v) for v in meet_irreducibles]

        S = SetPartitions(range(self.order()))
        min_congruences = []
        already_tried = []

        while irr:
            next_pair = irr.pop()
            cong = self.congruence([next_pair], stop_pairs=already_tried)
            already_tried.append(next_pair)
            if cong is not None:
                cong = S(cong)
                min_congruences = [c for c in min_congruences if c != cong and not S.is_less_than(cong, c)]
                if not any(S.is_less_than(c, cong) for c in min_congruences):
                    min_congruences.append(cong)

        return min_congruences

    def congruence(self, parts, start=None, stop_pairs=[]):
        """
        Return the congruence ``start`` "extended" by ``parts``.

        ``start`` is assumed to be a valid congruence of the lattice,
        and this is *not* checked.

        INPUT:

        - ``parts`` -- a list of lists; congruences to add
        - ``start`` -- a disjoint set; already computed congruence (or ``None``)
        - ``stop_pairs`` -- a list of pairs; list of pairs for stopping computation

        OUTPUT:

        ``None``, if the congruence generated by ``start`` and ``parts``
        together contains a block that has elements `a, b` so that ``(a, b)``
        is in the list ``stop_pairs``. Otherwise the least congruence that
        contains a block whose subset is `p` for every `p` in ``parts`` or
        ``start``, given as :class:`sage.sets.disjoint_set.DisjointSet_class`.

        ALGORITHM:

        Use the quadrilateral argument from page 120 of [Dav1997]_.

        Basically we take one block from todo-list, search quadrilateral
        blocks up and down against the block, and then complete them to
        closed intervals and add to todo-list.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0: [1, 2], 1: [3], 2: [4], 3: [4]})
            sage: cong = H.congruence([[0, 1]]); cong                                   # optional - sage.modules
            {{0, 1, 3}, {2, 4}}
            sage: H.congruence([[0, 2]], start=cong)                                    # optional - sage.modules
            {{0, 1, 2, 3, 4}}

            sage: H.congruence([[0, 1]], stop_pairs=[(1, 3)]) is None                   # optional - sage.modules
            True

        TESTS::

            sage: H = HasseDiagram('HT@O?GO?OE?G@??')
            sage: H.congruence([[0, 1]]).number_of_subsets()                            # optional - sage.modules
            1
            sage: H = HasseDiagram('HW_oC?@@O@?O@??')
            sage: H.congruence([[0, 1]]).number_of_subsets()                            # optional - sage.modules
            1

        Check :trac:`21861`::

            sage: H = HasseDiagram({0: [1, 2], 1: [3], 2: [4], 3: [4]})
            sage: tmp = H.congruence([[1, 3]])                                          # optional - sage.modules
            sage: tmp.number_of_subsets()                                               # optional - sage.modules
            4
            sage: H.congruence([[0, 1]], start=tmp).number_of_subsets()                 # optional - sage.modules
            2
            sage: tmp.number_of_subsets()                                               # optional - sage.modules
            4
        """
        from sage.sets.disjoint_set import DisjointSet
        from copy import copy

        n = self.order()
        mt = self.meet_matrix()
        jn = self.join_matrix()

        def fill_to_interval(S):
            """
            Return the smallest interval containing elements in the set S.
            """
            m = n - 1
            for e in S:
                m = mt[m, e]
            j = 0
            for e in S:
                j = jn[j, e]
            return self.interval(m, j)

        cong = copy(start) if start else DisjointSet(n)
        t = -1

        while t != cong.number_of_subsets():
            for part in parts:
                if part:  # Skip empty parts
                    c = part[0]
                    for e in fill_to_interval(part):
                        cong.union(e, c)
            t = cong.number_of_subsets()

            # Following is needed for cases like
            # posets.BooleanLattice(3).congruence([(0,1), (0,2), (0,4)])
            for c in list(cong):
                r = c[0]
                for v in fill_to_interval(c):
                    cong.union(r, v)

        todo = set(cong.find(e) for part in parts for e in part)

        while todo:

            # First check if we should stop now.
            for a, b in stop_pairs:
                if cong.find(a) == cong.find(b):
                    return None

            # We take one block and try to find as big interval
            # as possible to unify as a new block by the quadrilateral
            # argument.
            block = sorted(cong.root_to_elements_dict()[cong.find(todo.pop())])

            b = block[-1]
            for a in block:  # Quadrilateral up
                for c in self.neighbor_out_iterator(a):
                    if c not in block:
                        d = jn[c, b]
                        if cong.find(d) != cong.find(c):
                            break
                else:
                    continue
                break

            else:  # Not found, so...
                a = block[0]
                for b in reversed(block):  # ...quadrilateral down
                    for d in self.neighbor_in_iterator(b):
                        if d not in block:
                            c = mt[d, a]
                            if cong.find(c) != cong.find(d):
                                break
                    else:
                        continue
                    break
                else:  # Nothing found
                    continue

            # Something was found, so we put this block back to todo
            # together with just found new block.
            todo.add(a)
            todo.add(c)

            # Now the interval [c, d] will be of the same block.
            # It may "crab" other blocks within, and that can be
            # recursive process. In particular it may also combine to
            # [a, b] block we just used.
            while c is not None:
                newblock = cong.find(c)
                for i in self.interval(c, d):
                    cong.union(newblock, i)
                C = cong.root_to_elements_dict()[cong.find(newblock)]
                mins = [i for i in C if all(i_ not in C for i_ in self.neighbor_in_iterator(i))]
                maxs = [i for i in C if all(i_ not in C for i_ in self.neighbor_out_iterator(i))]
                c = None  # To stop loop, if this is not changed below.
                if len(mins) > 1 or len(maxs) > 1:
                    c = n - 1
                    for m in mins:
                        c = mt[c, m]
                    d = 0
                    for m in maxs:
                        d = jn[d, m]

            # This removes duplicates from todo.
            todo = set(cong.find(x) for x in todo)

        return cong

    def find_nontrivial_congruence(self):
        r"""
        Return a pair that generates non-trivial congruence or
        ``None`` if there is not any.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0: [1, 2], 1: [5], 2: [3, 4], 3: [5], 4: [5]})
            sage: H.find_nontrivial_congruence()                                        # optional - sage.modules
            {{0, 1}, {2, 3, 4, 5}}

            sage: H = HasseDiagram({0: [1, 2, 3], 1: [4], 2: [4], 3: [4]})
            sage: H.find_nontrivial_congruence() is None                                # optional - sage.modules
            True

        ALGORITHM:

        See https://www.math.hawaii.edu/~ralph/Preprints/conlat.pdf:

        If `\Theta` is a join irreducible element of a `\mathrm{Con}(L)`,
        then there is at least one join-irreducible `j` and one
        meet-irreducible `m` such that `\Theta` is both the principal
        congruence generated by `(j^*, j)`, where `j^*` is the unique
        lower cover of `j`, and the principal congruence generated by
        `(m, m^*)`, where `m^*` is the unique upper cover of `m`.

        So, we only check join irreducibles or meet irreducibles,
        whichever is a smaller set. To optimize more we stop computation
        whenever it finds a pair that we know to generate one-element
        congruence.
        """
        join_irreducibles = [v for v in self if self.in_degree(v) == 1]
        meet_irreducibles = [v for v in self if self.out_degree(v) == 1]
        if len(join_irreducibles) < len(meet_irreducibles):
            irr = [(v, self.neighbors_in(v)[0]) for v in join_irreducibles]
        else:
            irr = [(self.neighbors_out(v)[0], v) for v in meet_irreducibles]
        irr.sort(key=lambda x: x[0] - x[1])
        tried = []
        for pair in irr:
            cong = self.congruence([pair], stop_pairs=tried)
            if cong is not None and cong.number_of_subsets() > 1:
                return cong
            tried.append(pair)
        return None

    def principal_congruences_poset(self):
        r"""
        Return the poset of join-irreducibles of the congruence lattice.

        OUTPUT:

        A pair `(P, D)` where `P` is a poset and `D` is a dictionary.

        Elements of `P` are pairs `(x, y)` such that `x` is an element
        of the lattice and `y` is an element covering it. In the poset
        `(a, b)` is less than `(c, d)` iff the principal congruence
        generated by `(a, b)` is refinement of the principal congruence
        generated by `(c, d)`.

        `D` is a dictionary from pairs `(x, y)` to the congruence
        (given as DisjointSet) generated by the pair.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: N5 = HasseDiagram({0: [1, 2], 1: [4], 2: [3], 3: [4]})
            sage: P, D = N5.principal_congruences_poset()                               # optional - sage.combinat
            sage: P                                                                     # optional - sage.combinat
            Finite poset containing 3 elements
            sage: P.bottom()                                                            # optional - sage.combinat
            (2, 3)
            sage: D[(2, 3)]                                                             # optional - sage.combinat
            {{0}, {1}, {2, 3}, {4}}
        """
        from sage.combinat.set_partition import SetPartition, SetPartitions
        from sage.combinat.posets.posets import Poset

        n = self.order()

        # Select smaller set, meet- or join-irreducibles
        if self.in_degree_sequence().count(1) > self.out_degree_sequence().count(1):
            irr = [(e, next(self.neighbor_out_iterator(e))) for e in range(n) if self.out_degree(e) == 1]
        else:
            irr = [(next(self.neighbor_in_iterator(e)), e) for e in range(n) if self.in_degree(e) == 1]

        D = {}
        P = {}
        uniq_congs = set()
        for ab in irr:
            cong = self.congruence([ab])
            cong_ = SetPartition(cong)
            if cong_ not in uniq_congs:
                uniq_congs.add(cong_)
                D[ab] = cong
                P[ab] = cong_

        # TODO: Make a function that creates the poset from a set
        # by comparison function with minimal number of comparisons.

        T = SetPartitions(n)
        P = DiGraph([D, lambda a, b: T.is_less_than(P[a], P[b])])
        return (Poset(P), D)

    def congruences_iterator(self):
        """
        Return an iterator over all congruences of the lattice.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram('GY@OQ?OW@?O?')
            sage: it = H.congruences_iterator(); it
            <generator object ...>
            sage: sorted([cong.number_of_subsets() for cong in it])                     # optional - sage.combinat
            [1, 2, 2, 2, 4, 4, 4, 8]
        """
        from sage.sets.disjoint_set import DisjointSet

        P, congs = self.principal_congruences_poset()
        for a in P.antichains_iterator():
            achain = tuple(a)
            n = len(achain)
            if n == 0:
                yield DisjointSet(self.order())
            if n == 1:
                # We have congs[(x,y)], but we want congs[((x,y))].
                congs[achain] = congs[a[0]]
                yield congs[achain[0]]
            if n > 1:
                c = congs[achain[:-1]]
                c = self.congruence([achain[-1]], start=c)
                yield c
                congs[achain] = c

    def is_congruence_normal(self) -> bool:
        """
        Return ``True`` if the lattice can be constructed from the one-element
        lattice with Day doubling constructions of convex subsets.

        Subsets to double does not need to be lower nor upper pseudo-intervals.
        On the other hand they must be convex, i.e. doubling a non-convex but
        municipal subset will give a lattice that returns ``False`` from
        this function.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram('IX?Q@?AG?OG?W?O@??')
            sage: H.is_congruence_normal()                                              # optional - sage.combinat
            True

        The 5-element diamond is the smallest non-example::

            sage: H = HasseDiagram({0: [1, 2, 3], 1: [4], 2: [4], 3: [4]})
            sage: H.is_congruence_normal()                                              # optional - sage.combinat
            False

        This is done by doubling a non-convex subset::

            sage: H = HasseDiagram('OQC?a?@CO?G_C@?GA?O??_??@?BO?A_?G??C??_?@???')
            sage: H.is_congruence_normal()                                              # optional - sage.combinat
            False

        TESTS::

            sage: HasseDiagram().is_congruence_normal()                                 # optional - sage.combinat
            True
            sage: HasseDiagram({0: []}).is_congruence_normal()                          # optional - sage.combinat
            True

        ALGORITHM:

        See http://www.math.hawaii.edu/~jb/inflation.pdf
        """
        from sage.combinat.set_partition import SetPartition

        n = self.order()
        congs_ji: dict[SetPartition, list] = {}

        for ji in range(n):
            if self.in_degree(ji) == 1:
                cong = SetPartition(self.congruence([[ji, next(self.neighbor_in_iterator(ji))]]))  # type:ignore
                if cong not in congs_ji:
                    congs_ji[cong] = []
                congs_ji[cong].append(ji)

        for mi in range(n):
            if self.out_degree(mi) == 1:
                cong = SetPartition(self.congruence([[mi, next(self.neighbor_out_iterator(mi))]]))  # type:ignore
                if any(self.is_lequal(ji, mi) for ji in congs_ji[cong]):
                    return False

        return True

    @staticmethod
    def _glue_spectra(a_spec, b_spec, orientation):
        r"""
        Return the `a`-spectrum of a poset by merging ``a_spec`` and ``b_spec``.

        ``a_spec`` and ``b_spec`` are the `a`-spectrum and `b`-spectrum of two different
        posets (see :meth:`atkinson` for the definition of `a`-spectrum).

        The orientation determines whether `a < b` or `b < a` in the combined poset.

        This is a helper method for :meth:`atkinson`.

        INPUT:

        - ``a_spec`` -- list; the `a`-spectrum of a poset `P`

        - ``b_spec`` -- list; the `b`-spectrum of a poset `Q`

        - ``orientation`` -- boolean; ``True`` if `a < b`, ``False`` otherwise

        OUTPUT:

        The `a`-spectrum (or `b`-spectrum, depending on orientation),
        returned as a list, of the poset which is a disjoint union
        of `P` and `Q`, together with the additional
        covering relation `a < b`.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: Pdata = [0, 1, 2, 0]
            sage: Qdata = [1, 1, 0]
            sage: HasseDiagram._glue_spectra(Pdata, Qdata, True)
            [0, 20, 28, 18, 0, 0, 0]

            sage: Pdata = [0, 0, 2]
            sage: Qdata = [0, 1]
            sage: HasseDiagram._glue_spectra(Pdata, Qdata, False)
            [0, 0, 0, 0, 8]
        """
        new_a_spec = []

        if not orientation:
            a_spec, b_spec = b_spec, a_spec

        p = len(a_spec)
        q = len(b_spec)

        for r in range(1, p + q + 1):
            new_a_spec.append(0)
            for i in range(max(1, r - q), min(p, r) + 1):
                k_val = binomial(r - 1, i - 1) * binomial(p + q - r, p - i)
                if orientation:
                    inner_sum = sum(b_spec[j - 1] for j in range(r - i + 1, len(b_spec) + 1))
                else:
                    inner_sum = sum(b_spec[j - 1] for j in range(1, r - i + 1))
                new_a_spec[-1] = new_a_spec[-1] + (a_spec[i - 1] * k_val * inner_sum)

        return new_a_spec

    def _split(self, a, b):
        r"""
        Return the two connected components obtained by deleting the covering
        relation `a < b` from a Hasse diagram that is a tree.

        This is a helper method for :meth:`FinitePoset.atkinson`.

        INPUT:

        - ``a`` -- an element of the poset
        - ``b`` -- an element of the poset which covers ``a``

        OUTPUT:

        A list containing two posets which are the connected components
        of this poset after deleting the covering relation `a < b`.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0: [1, 2], 1: [], 2: []})
            sage: H._split(0, 1)
            [Hasse diagram of a poset containing 2 elements, Hasse diagram of a poset containing 1 elements]

            sage: H = HasseDiagram({0: [1], 1: [2], 2: [3], 3: [4], 4: []})
            sage: H._split(1, 2)
            [Hasse diagram of a poset containing 2 elements, Hasse diagram of a poset containing 3 elements]

        TESTS::
            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0: [1,2,3], 1: [4, 5], 2: [4, 6], 3: [5, 6], 4: [7], 5: [7], 6: [7], 7: []})
            sage: H._split(0, 1)
            Traceback (most recent call last):
            ...
            ValueError: wrong number of connected components after the covering relation is deleted

            sage: H = HasseDiagram({0: [1], 1: [], 2: []})
            sage: H._split(0, 1)
            Traceback (most recent call last):
            ...
            ValueError: wrong number of connected components after the covering relation is deleted
        """
        split_hasse = self.copy(self)
        split_hasse.delete_edge(a, b)
        components = split_hasse.connected_components_subgraphs()
        if not len(components) == 2:
            raise ValueError("wrong number of connected components after the covering relation is deleted")

        c1, c2 = components
        if a in c2:
            c1, c2 = c2, c1

        return [c1, c2]

    def _spectrum_of_tree(self, a):
        r"""
        Return the `a`-spectrum of a poset whose underlying graph is a tree.

        This is a helper method for :meth:`FinitePoset.atkinson`.

        INPUT:

        - ``a`` -- an element of the poset

        OUTPUT:

        The `a`-spectrum of this poset, returned as a list.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0: [2], 1: [2], 2: [3, 4], 3: [], 4: []})
            sage: H._spectrum_of_tree(0)
            [2, 2, 0, 0, 0]

            sage: H = HasseDiagram({0: [2], 1: [2], 2: [3, 4], 3: [], 4: []})
            sage: H._spectrum_of_tree(2)
            [0, 0, 4, 0, 0]

            sage: H = HasseDiagram({0: [2], 1: [2], 2: [3, 4], 3: [], 4: []})
            sage: H._spectrum_of_tree(3)
            [0, 0, 0, 2, 2]
        """
        upper_covers = self.neighbors_out(a)
        lower_covers = self.neighbors_in(a)
        if not upper_covers and not lower_covers:
            return [1]
        if upper_covers:
            b = upper_covers[0]
            orientation = True
        else:
            (a, b) = (lower_covers[0], a)
            orientation = False
        P, Q = self._split(a, b)
        a_spec = P._spectrum_of_tree(a)
        b_spec = Q._spectrum_of_tree(b)
        return HasseDiagram._glue_spectra(a_spec, b_spec, orientation)


__doc__ = __doc__.format(INDEX_OF_FUNCTIONS=gen_rest_table_index(HasseDiagram))
