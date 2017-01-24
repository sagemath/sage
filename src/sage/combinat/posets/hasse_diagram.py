# -*- coding: utf-8 -*-
r"""
Hasse diagrams of posets

{INDEX_OF_FUNCTIONS}

"""

#*****************************************************************************
#       Copyright (C) 2008 Peter Jipsen <jipsen@chapman.edu>
#       Copyright (C) 2008 Franco Saliola <saliola@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

from six.moves import range

from sage.graphs.digraph import DiGraph
from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.cachefunc import cached_method
from sage.misc.superseded import deprecated_function_alias

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
        return "Hasse diagram of a poset containing %s elements"%self.order()

    def linear_extension(self):
        r"""
        Return a linear extension

        TESTS::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,2],1:[3],2:[3],3:[]})
            sage: H.linear_extension()
            [0, 1, 2, 3]
        """
        # Recall: we assume range(n) is a linear extension.
        return list(range(len(self)))

    def linear_extensions(self):
        r"""
        Return all linear extensions

        TESTS::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,2],1:[3],2:[3],3:[]})
            sage: H.linear_extensions()
            [[0, 1, 2, 3], [0, 2, 1, 3]]
        """
        return self.topological_sort_generator()

    def is_linear_extension(self, lin_ext=None):
        r"""
        Test if an ordering is a linear extension.

        TESTS::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,2],1:[3],2:[3],3:[]})
            sage: H.is_linear_extension(list(range(4)))
            True
            sage: H.is_linear_extension([3,2,1,0])
            False
        """
        if lin_ext is None or lin_ext == list(range(len(self))):
            for x, y in self.cover_relations_iterator():
                if not x < y:
                    return False
            return True
        else:
            for x, y in self.cover_relations_iterator():
                if not lin_ext.index(x) < lin_ext.index(y):
                    return False
            return True

    def cover_relations_iterator(self):
        r"""
        Iterate over cover relations.

        TESTS::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[2,3], 1:[3,4], 2:[5], 3:[5], 4:[5]})
            sage: list(H.cover_relations_iterator())
            [(0, 2), (0, 3), (1, 3), (1, 4), (2, 5), (3, 5), (4, 5)]
        """
        for u,v,l in self.edge_iterator():
            yield (u,v)

    def cover_relations(self):
        r"""
        Return the list of cover relations.

        TESTS::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[2,3], 1:[3,4], 2:[5], 3:[5], 4:[5]})
            sage: H.cover_relations()
            [(0, 2), (0, 3), (1, 3), (1, 4), (2, 5), (3, 5), (4, 5)]
        """
        return list(self.cover_relations_iterator())

    def is_lequal(self, i, j):
        """
        Returns True if i is less than or equal to j in the poset, and
        False otherwise.

        .. note::

            If the :meth:`lequal_matrix` has been computed, then this method is
            redefined to use the cached matrix (see :meth:`_alternate_is_lequal`).

        TESTS::

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
        return i == j or \
                (i < j and j in self.breadth_first_search(i))

    def is_less_than(self, x, y):
        r"""
        Returns True if ``x`` is less than or equal to ``y`` in the
        poset, and False otherwise.

        TESTS::

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
        else:
            return self.is_lequal(x,y)

    def is_gequal(self, x, y):
        r"""
        Returns ``True`` if ``x`` is greater than or equal to ``y``, and
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
        return self.is_lequal(y,x)

    def is_greater_than(self, x, y):
        """
        Returns ``True`` if ``x`` is greater than but not equal to
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
        return self.is_less_than(y,x)

    def minimal_elements(self):
        """
        Returns a list of the minimal elements of the poset.

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
        Returns a list of the maximal elements of the poset.

        EXAMPLES::

            sage: P = Poset({0:[3],1:[3],2:[3],3:[4],4:[]})
            sage: P.maximal_elements()
            [4]
        """
        return self.sinks()

    def bottom(self):
        """
        Returns the bottom element of the poset, if it exists.

        EXAMPLES::

            sage: P = Poset({0:[3],1:[3],2:[3],3:[4],4:[]})
            sage: P.bottom() is None
            True
            sage: Q = Poset({0:[1],1:[]})
            sage: Q.bottom()
            0
        """
        min_elms = self.minimal_elements()
        if len(min_elms) == 1: return min_elms[0]
        return None

    def has_bottom(self):
        """
        Returns True if the poset has a unique minimal element.

        EXAMPLES::

            sage: P = Poset({0:[3],1:[3],2:[3],3:[4],4:[]})
            sage: P.has_bottom()
            False
            sage: Q = Poset({0:[1],1:[]})
            sage: Q.has_bottom()
            True
        """
        if self.bottom() is not None: return True
        return False

    def top(self):
        """
        Returns the top element of the poset, if it exists.

        EXAMPLES::

            sage: P = Poset({0:[3],1:[3],2:[3],3:[4,5],4:[],5:[]})
            sage: P.top() is None
            True
            sage: Q = Poset({0:[1],1:[]})
            sage: Q.top()
            1
        """
        max_elms = self.maximal_elements()
        if len(max_elms) == 1: return max_elms[0]
        return None

    def has_top(self):
        """
        Returns ``True`` if the poset contains a unique maximal element, and
        ``False`` otherwise.

        EXAMPLES::

            sage: P = Poset({0:[3],1:[3],2:[3],3:[4,5],4:[],5:[]})
            sage: P.has_top()
            False
            sage: Q = Poset({0:[1],1:[]})
            sage: Q.has_top()
            True
        """
        if not self.top() is None: return True
        return False

    def is_bounded(self):
        """
        Returns True if the poset contains a unique maximal element and a
        unique minimal element, and False otherwise.

        EXAMPLES::

            sage: P = Poset({0:[3],1:[3],2:[3],3:[4,5],4:[],5:[]})
            sage: P.is_bounded()
            False
            sage: Q = Poset({0:[1],1:[]})
            sage: Q.is_bounded()
            True
        """
        return self.has_top() and self.has_bottom()

    def is_chain(self):
        """
        Returns True if the poset is totally ordered, and False otherwise.

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
        return (self.num_edges()+1 == self.num_verts() and # Hasse Diagram is a tree
                all(d<=1 for d in self.out_degree())   and # max outdegree is <= 1
                all(d<=1 for d in self.in_degree()))       # max  indegree is <= 1

    def is_antichain_of_poset(self, elms):
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
        from sage.misc.misc import uniq

        elms_sorted = uniq(elms)
        return not any(self.is_lequal(a, b) for a, b in
                       combinations(elms_sorted, 2))

    def dual(self):
        """
        Returns a poset that is dual to the given poset.

        EXAMPLES::

            sage: P = Posets.IntegerPartitions(4)
            sage: H = P._hasse_diagram; H
            Hasse diagram of a poset containing 5 elements
            sage: H.dual()
            Hasse diagram of a poset containing 5 elements

        TESTS::

            sage: H = Posets.IntegerPartitions(4)._hasse_diagram
            sage: H.is_isomorphic( H.dual().dual() )
            True
            sage: H.is_isomorphic( H.dual() )
            False
        """
        H = self.reverse()
        H.relabel(perm=list(range(H.num_verts()-1, -1, -1)), inplace=True)
        return HasseDiagram(H)

    def interval(self, x, y):
        """
        Return a list of the elements `z` of ``self`` such that
        `x \leq z \leq y`. The order is that induced by the
        ordering in ``self.linear_extension``.

        INPUT:

        -  ``x`` -- any element of the poset

        -  ``y`` -- any element of the poset

        EXAMPLES::

            sage: uc = [[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]]
            sage: dag = DiGraph(dict(zip(range(len(uc)),uc)))
            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram(dag)
            sage: I = set([2,5,6,4,7])
            sage: I == set(H.interval(2,7))
            True
        """
        return [z for z in range(x, y+1) if
                self.is_lequal(x, z) and self.is_lequal(z, y)]

    closed_interval = interval

    def open_interval(self, x, y):
        """
        Return a list of the elements `z` of ``self`` such that
        `x < z < y`. The order is that induced by the ordering in
        ``self.linear_extension``.

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
        ci = self.interval(x,y)
        if len(ci) == 0:
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

        test for ticket :trac:`14006`::

            sage: H = Poset()._hasse_diagram
            sage: s = dumps(H)
            sage: f = H.rank_function()
            sage: s = dumps(H)
        """
        if(self._rank is None):
            return None
        return self._rank.__getitem__ # the rank function is just the getitem of the list

    @lazy_attribute
    def _rank(self):
        r"""
        Builds the rank function of the poset, if it exists, i.e.
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
        rank = [None]*self.order()
        not_found = set(self.vertices())
        while not_found:
            y = not_found.pop()
            rank[y] = 0  # We set some vertex to have rank 0
            component = set([y])
            queue = set([y])
            while queue:  # look at the neighbors of y and set the ranks;
                          # then look at the neighbors of the neighbors ...
                y = queue.pop()
                for x in self.neighbors_out(y):
                    if rank[x] is None:
                        rank[x] = rank[y] + 1
                        queue.add(x)
                        component.add(x)
                for x in self.neighbors_in(y):
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
        #now, all ranks are set.
        return rank

    def rank(self,element=None):
        r"""
        Returns the rank of ``element``, or the rank of the poset if
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
            return len(self.level_sets())-1
        else:
            return self.rank_function()(element)

    def is_ranked(self):
        r"""
        Returns True if the poset is ranked, and False otherwise.

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

    def covers(self,x,y):
        """
        Returns True if y covers x and False otherwise.

        EXAMPLES::

            sage: Q = Poset([[1,5],[2,6],[3],[4],[],[6,3],[4]])
            sage: Q.covers(Q(1),Q(6))
            True
            sage: Q.covers(Q(1),Q(4))
            False
        """
        return self.has_edge(x,y)

    def upper_covers_iterator(self,element):
        r"""
        Returns the list of elements that cover ``element``.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,3,2],1:[4],2:[4,5,6],3:[6],4:[7],5:[7],6:[7],7:[]})
            sage: list(H.upper_covers_iterator(0))
            [1, 2, 3]
            sage: list(H.upper_covers_iterator(7))
            []
        """
        for x in self.neighbor_out_iterator(element):
            yield x

    def lower_covers_iterator(self,element):
        r"""
        Returns the list of elements that are covered by ``element``.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,3,2],1:[4],2:[4,5,6],3:[6],4:[7],5:[7],6:[7],7:[]})
            sage: list(H.lower_covers_iterator(0))
            []
            sage: list(H.lower_covers_iterator(4))
            [1, 2]
        """
        for x in self.neighbor_in_iterator(element):
            yield x

    def cardinality(self):
        r"""
        Returns the number of elements in the poset.

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

            sage: L = Posets.BooleanLattice(5)
            sage: H = L.hasse_diagram()
            sage: H.size()
            80
            sage: H.size() == H.num_edges()
            True
        """
        return self.order()

    def moebius_function(self,i,j): # dumb algorithm
        r"""
        Returns the value of the Möbius function of the poset
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
            return self._moebius_function_values[(i,j)]
        except AttributeError:
            self._moebius_function_values = {}
            return self.moebius_function(i,j)
        except KeyError:
            if i == j:
                self._moebius_function_values[(i,j)] = 1
            elif i > j:
                self._moebius_function_values[(i,j)] = 0
            else:
                ci = self.closed_interval(i,j)
                if len(ci) == 0:
                    self._moebius_function_values[(i,j)] = 0
                else:
                    self._moebius_function_values[(i,j)] = \
                     -sum([self.moebius_function(i,k) for k in ci[:-1]])
        return self._moebius_function_values[(i,j)]
    mobius_function = deprecated_function_alias(19855, moebius_function)

    def moebius_function_matrix(self):
        r"""
        Returns the matrix of the Möbius function of this poset

        This returns the sparse matrix over `\ZZ` whose ``(x, y)`` entry
        is the value of the Möbius function of ``self`` evaluated on
        ``x`` and ``y``, and redefines :meth:`moebius_function` to use
        it.

        .. NOTE::

            The result is cached in :meth:`_moebius_function_matrix`.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,3,2],1:[4],2:[4,5,6],3:[6],4:[7],5:[7],6:[7],7:[]})
            sage: H.moebius_function_matrix()
            [ 1 -1 -1 -1  1  0  1  0]
            [ 0  1  0  0 -1  0  0  0]
            [ 0  0  1  0 -1 -1 -1  2]
            [ 0  0  0  1  0  0 -1  0]
            [ 0  0  0  0  1  0  0 -1]
            [ 0  0  0  0  0  1  0 -1]
            [ 0  0  0  0  0  0  1 -1]
            [ 0  0  0  0  0  0  0  1]

        TESTS::

            sage: H.moebius_function_matrix().is_immutable()
            True
            sage: hasattr(H,'_moebius_function_matrix')
            True

            sage: H.moebius_function == H._moebius_function_from_matrix
            True
        """
        if not hasattr(self,'_moebius_function_matrix'):
            self._moebius_function_matrix = self.lequal_matrix().inverse().change_ring(ZZ)
            self._moebius_function_matrix.set_immutable()
            self.moebius_function = self._moebius_function_from_matrix
        return self._moebius_function_matrix
    mobius_function_matrix = deprecated_function_alias(19855, moebius_function_matrix)

    # Redefine self.moebius_function
    def _moebius_function_from_matrix(self, i,j):
        r"""
        Returns the value of the Möbius function of the poset
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
        return self._moebius_function_matrix[i,j]
    _mobius_function_from_matrix = deprecated_function_alias(19855, _moebius_function_from_matrix)

    @cached_method
    def coxeter_transformation(self):
        r"""
        Returns the matrix of the Auslander-Reiten translation acting on
        the Grothendieck group of the derived category of modules on the
        poset, in the basis of simple modules.

        EXAMPLES::

            sage: M = Posets.PentagonPoset()._hasse_diagram.coxeter_transformation(); M
            [ 0  0  0  0 -1]
            [ 0  0  0  1 -1]
            [ 0  1  0  0 -1]
            [-1  1  1  0 -1]
            [-1  1  0  1 -1]

        TESTS::

            sage: M = Posets.PentagonPoset()._hasse_diagram.coxeter_transformation()
            sage: M**8 == 1
            True
        """
        return - self.lequal_matrix()*self.moebius_function_matrix().transpose()

    def order_filter(self, elements):
        """
        Return the order filter generated by a list of elements.

        `I` is an order filter if, for any `x` in `I` and `y` such that
        `y \ge x`, then `y` is in `I`.

        EXAMPLES::

            sage: H = Posets.BooleanLattice(4)._hasse_diagram
            sage: H.order_filter([3,8])
            [3, 7, 8, 9, 10, 11, 12, 13, 14, 15]
        """
        return sorted(list(self.depth_first_search(elements)))

    def principal_order_filter(self, i):
        """
        Returns the order filter generated by ``i``.

        EXAMPLES::

            sage: H = Posets.BooleanLattice(4)._hasse_diagram
            sage: H.principal_order_filter(2)
            [2, 3, 6, 7, 10, 11, 14, 15]
        """
        return self.order_filter([i])

    def order_ideal(self, elements):
        """
        Return the order ideal generated by a list of elements.

        `I` is an order ideal if, for any `x` in `I` and `y` such that
        `y \le x`, then `y` is in `I`.

        EXAMPLES::

            sage: H = Posets.BooleanLattice(4)._hasse_diagram
            sage: H.order_ideal([7,10])
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 10]
        """
        return sorted(list(
            self.depth_first_search(elements, neighbors=self.neighbors_in)))

    def principal_order_ideal(self, i):
        """
        Returns the order ideal generated by `i`.

        EXAMPLES::

            sage: H = Posets.BooleanLattice(4)._hasse_diagram
            sage: H.principal_order_ideal(6)
            [0, 2, 4, 6]
        """
        return self.order_ideal([i])

    @lazy_attribute
    def _leq_matrix(self):
        r"""
        Computes a matrix whose ``(i,j)`` entry is 1 if ``i`` is less than
        ``j`` in the poset, and 0 otherwise; and redefines ``__lt__`` to
        use this matrix.

        EXAMPLES::

            sage: P = Poset([[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]])
            sage: H = P._hasse_diagram
            sage: H._leq_matrix
            [1 1 1 1 1 1 1 1]
            [0 1 0 1 0 0 0 1]
            [0 0 1 1 1 0 1 1]
            [0 0 0 1 0 0 0 1]
            [0 0 0 0 1 0 0 1]
            [0 0 0 0 0 1 1 1]
            [0 0 0 0 0 0 1 1]
            [0 0 0 0 0 0 0 1]

        """
        # Create the matrix
        n = self.order()
        D = {}
        for i in range(n):
            for v in self.breadth_first_search(i):
                D[(i,v)] = 1
        M = matrix(ZZ, n, n, D, sparse=True)
        M.set_immutable()
        # Redefine self.is_lequal
        self.is_lequal = self._alternate_is_lequal
        # Return the matrix
        return M

    def lequal_matrix(self):
        """
        Returns the matrix whose ``(i,j)`` entry is 1 if ``i`` is less
        than ``j`` in the poset, and 0 otherwise; and redefines
        ``__lt__`` to use this matrix.

        EXAMPLES::

            sage: P = Poset([[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]])
            sage: H = P._hasse_diagram
            sage: H.lequal_matrix()
            [1 1 1 1 1 1 1 1]
            [0 1 0 1 0 0 0 1]
            [0 0 1 1 1 0 1 1]
            [0 0 0 1 0 0 0 1]
            [0 0 0 0 1 0 0 1]
            [0 0 0 0 0 1 1 1]
            [0 0 0 0 0 0 1 1]
            [0 0 0 0 0 0 0 1]

        TESTS::

            sage: H.lequal_matrix().is_immutable()
            True
        """
        return self._leq_matrix

    def _alternate_is_lequal(self,i,j):
        r"""
        Returns ``True`` if ``i`` is less than or equal to ``j`` in
        ``self``, and ``False`` otherwise.

        .. NOTE::

            If the :meth:`lequal_matrix` has been computed, then
            :meth:`is_lequal` is redefined to use the cached matrix.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: H.lequal_matrix()
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
        return bool(self._leq_matrix[i,j])

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
        Return the matrix of meets of ``self``. The ``(x,y)``-entry of
        this matrix is the meet of ``x`` and ``y`` in ``self``.

        EXAMPLES::

           sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
           sage: H = HasseDiagram({0:[1,3,2],1:[4],2:[4,5,6],3:[6],4:[7],5:[7],6:[7],7:[]})
           sage: H._meet
           [0 0 0 0 0 0 0 0]
           [0 1 0 0 1 0 0 1]
           [0 0 2 0 2 2 2 2]
           [0 0 0 3 0 0 3 3]
           [0 1 2 0 4 2 2 4]
           [0 0 2 0 2 5 2 5]
           [0 0 2 3 2 2 6 6]
           [0 1 2 3 4 5 6 7]

        TESTS::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[2,3],1:[2,3]})
            sage: H.meet_matrix()
            Traceback (most recent call last):
            ...
            ValueError: not a meet-semilattice: no bottom element

            sage: H = HasseDiagram({0:[1,2],1:[3,4],2:[3,4]})
            sage: H.meet_matrix()
            Traceback (most recent call last):
            ...
            LatticeError: no meet for ...

            sage: L = LatticePoset({0:[1,2,3],1:[4],2:[4],3:[4]})
            sage: P = L.dual()
            sage: P.meet(2,3)
            4
        """
        n = self.cardinality()
        if n == 0:
            return matrix(0)
        if not self.has_bottom():
            raise ValueError("not a meet-semilattice: no bottom element")
        meet = [[0 for x in range(n)] for x in range(n)]
        lc = [self.neighbors_in(x) for x in range(n)]  # Lc = lower covers

        for x in range(n):
            meet[x][x] = x
            for y in range(x):
                T = [meet[y][z] for z in lc[x]]

                q = max(T)
                for z in T:
                    if meet[z][q] != z:
                        raise LatticeError('meet', x, y)
                meet[x][y] = q
                meet[y][x] = q

        return matrix(ZZ, meet)

    def meet_matrix(self):
        r"""
        Returns the matrix of meets of ``self``. The ``(x,y)``-entry of
        this matrix is the meet of ``x`` and ``y`` in ``self``.

        This algorithm is modelled after the algorithm of Freese-Jezek-Nation
        (p217). It can also be found on page 140 of [Gec81]_.

        .. NOTE::

            Once the matrix has been computed, it is stored in
            :meth:`_meet_matrix`. Delete this attribute if you want to
            recompute the matrix.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,3,2],1:[4],2:[4,5,6],3:[6],4:[7],5:[7],6:[7],7:[]})
            sage: H.meet_matrix()
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
            sage: H.meet_matrix()
            Traceback (most recent call last):
            ...
            LatticeError: no meet for ...
        """
        return self._meet

    def is_meet_semilattice(self):
        r"""
        Returns ``True`` if ``self`` has a meet operation, and
        ``False`` otherwise.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,3,2],1:[4],2:[4,5,6],3:[6],4:[7],5:[7],6:[7],7:[]})
            sage: H.is_meet_semilattice()
            True

            sage: H = HasseDiagram({0:[1,2],1:[3],2:[3],3:[]})
            sage: H.is_meet_semilattice()
            True

            sage: H = HasseDiagram({0:[2,3],1:[2,3]})
            sage: H.is_meet_semilattice()
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
        and ``y`` in ``self``

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,3,2],1:[4],2:[4,5,6],3:[6],4:[7],5:[7],6:[7],7:[]})
            sage: H.join_matrix() # indirect doctest
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
            sage: H.join_matrix()
            Traceback (most recent call last):
            ...
            LatticeError: no join for ...

            sage: L = LatticePoset({0:[1,2,3],1:[4],2:[4],3:[4]})
            sage: P = L.dual()
            sage: P.join(2,3)
            0
        """
        n = self.cardinality()
        if n == 0:
            return matrix(0)
        if not self.has_top():
            raise ValueError("not a join-semilattice: no top element")
        join = [[n for x in range(n)] for x in range(n)]
        uc = [self.neighbors_out(x) for x in range(n)]  # uc = upper covers

        for x in range(n-1, -1, -1):
            join[x][x] = x
            for y in range(n-1, x, -1):
                T = [join[y][z] for z in uc[x]]

                q = min(T)
                for z in T:
                    if join[z][q] != z:
                        raise LatticeError('join', x, y)
                join[x][y] = q
                join[y][x] = q

        return matrix(ZZ, join)

    def join_matrix(self):
        r"""
        Returns the matrix of joins of ``self``. The ``(x,y)``-entry
        of this matrix is the join of ``x`` and ``y`` in ``self``.

        This algorithm is modelled after the algorithm of Freese-Jezek-Nation
        (p217). It can also be found on page 140 of [Gec81]_.

        .. note::

            Once the matrix has been computed, it is stored in
            :meth:`_join_matrix`. Delete this attribute if you want
            to recompute the matrix.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,3,2],1:[4],2:[4,5,6],3:[6],4:[7],5:[7],6:[7],7:[]})
            sage: H.join_matrix()
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
            sage: H.join_matrix()
            Traceback (most recent call last):
            ...
            LatticeError: no join for ...
        """
        return self._join

    def is_join_semilattice(self):
        r"""
        Returns ``True`` if ``self`` has a join operation, and
        ``False`` otherwise.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,3,2],1:[4],2:[4,5,6],3:[6],4:[7],5:[7],6:[7],7:[]})
            sage: H.is_join_semilattice()
            True
            sage: H = HasseDiagram({0:[2,3],1:[2,3]})
            sage: H.is_join_semilattice()
            False
            sage: H = HasseDiagram({0:[2,3],1:[2,3],2:[4],3:[4]})
            sage: H.is_join_semilattice()
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
            sage: H.find_nonsemidistributive_elements('join') is None
            False
            sage: H.find_nonsemidistributive_elements('meet') is None
            True
        """
        if meet_or_join == 'join':
            M1 = self._join
            M2 = self._meet
        elif meet_or_join == 'meet':
            M1 = self._meet
            M2 = self._join
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

    def is_distributive_lattice(self): # still a dumb algorithm...
        r"""
        Returns ``True`` if ``self`` is the Hasse diagram of a
        distributive lattice, and ``False`` otherwise.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,3,2],1:[4],2:[4,5,6],3:[6],4:[7],5:[7],6:[7],7:[]})
            sage: H.is_distributive_lattice()
            False
            sage: H = HasseDiagram({0:[1,2],1:[3],2:[3]})
            sage: H.is_distributive_lattice()
            True
            sage: H = HasseDiagram({0:[1,2,3],1:[4],2:[4],3:[4]})
            sage: H.is_distributive_lattice()
            False
        """
        try:
            jn = self.join_matrix()
            mt = self.meet_matrix()
        except ValueError:
            return False
        n = jn.ncols()
        for x in range(n):
            for y in range(n):
                for z in range(n):
                    if mt[x][jn[y][z]]!=jn[mt[x][y]][mt[x][z]]: return False
        return True

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

            sage: H = Posets.BooleanLattice(4)._hasse_diagram
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
        result = [] # Never take the bottom element to list.
        e = 0
        m = 0
        for i in range(n-1):
            for j in self.outgoing_edge_iterator(i):
                m = max(m, j[1])
            if m == i+1:
                if not return_list:
                    if m < n-1:
                        return m
                    else:
                        return None
                result.append(m)
        result.pop() # Remove the top element.
        return result

    def is_complemented(self):
        """
        Return an element of the lattice that has no complement.

        If the lattice is complemented, return ``None``.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram

            sage: H = HasseDiagram({0:[1, 2], 1:[3], 2:[3], 3:[4]})
            sage: H.is_complemented()
            1

            sage: H = HasseDiagram({0:[1, 2, 3], 1:[4], 2:[4], 3:[4]})
            sage: H.is_complemented() is None
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
            sage: H.pseudocomplement(2)
            3

            sage: H = HasseDiagram({0: [1, 2, 3], 1: [4], 2: [4], 3: [4]})
            sage: H.pseudocomplement(2) is None
            True
        """
        e = self.order() - 1
        while self._meet[e, element] != 0:
            e -= 1
        e1 = e
        while e1 > 0:
            if self._meet[e1, element] == 0 and not self.is_lequal(e1, e):
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
            sage: list(H.orthocomplementations_iterator())
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
            sage: list(H.orthocomplementations_iterator())
            [[]]
            sage: H = HasseDiagram({0:[]})  # One element
            sage: list(H.orthocomplementations_iterator())
            [[0]]
            sage: H = HasseDiagram({0:[1]})  # Two elements
            sage: list(H.orthocomplementations_iterator())
            [[1, 0]]

        Trivial cases: odd number of elements, not self-dual, not complemented::

            sage: H = Posets.DiamondPoset(5)._hasse_diagram
            sage: list(H.orthocomplementations_iterator())
            []
            sage: H = Posets.ChainPoset(4)._hasse_diagram
            sage: list(H.orthocomplementations_iterator())
            []
            sage: H = HasseDiagram( ([[0, 1], [0, 2], [0, 3], [1, 4], [1, 8], [4, 6], [4, 7], [6, 9], [7, 9], [2, 5], [3, 5], [5, 8], [8, 9]]) )
            sage: list(H.orthocomplementations_iterator())
            []
            sage: H = HasseDiagram({0:[1, 2, 3], 1: [4], 2:[4], 3: [5], 4:[5]})
            sage: list(H.orthocomplementations_iterator())
            []

        Complemented, self-dual and even number of elements, but
        not orthocomplemented::

            sage: H = HasseDiagram( ([[0, 1], [1, 2], [2, 3], [0, 4], [4, 5], [0, 6], [3, 7], [5, 7], [6, 7]]) )
            sage: list(H.orthocomplementations_iterator())
            []

        Unique orthocomplementations; second is not uniquely complemented,
        but has only one orthocomplementation.

            sage: H = Posets.BooleanLattice(4)._hasse_diagram  # Uniquely complemented
            sage: len(list(H.orthocomplementations_iterator()))
            1
            sage: H = HasseDiagram({0:[1, 2], 1:[3], 2:[4], 3:[5], 4:[5]})
            sage: len([_ for _ in H.orthocomplementations_iterator()])
            1

        "Lengthening diamond" must keep the number of orthocomplementations::

            sage: H = HasseDiagram( ([[0, 1], [0, 2], [0, 3], [0, 4], [1, 5], [2, 5], [3, 5], [4, 5]]) )
            sage: n = len([_ for _ in H.orthocomplementations_iterator()]); n
            3
            sage: H = HasseDiagram('M]??O?@??C??OA???OA??@?A??C?A??O??')
            sage: len([_ for _ in H.orthocomplementations_iterator()]) == n
            True

        This lattice has an unique "possible orthocomplement" for every
        element, but they can not be fit together; orthocomplement pairs
        would be 0-11, 1-7, 2-4, 3-10, 5-9 and 6-8, and then orthocomplements
        for chain 0-1-6-11 would be 11-7-8-0, which is not a chain::

            sage: H = HasseDiagram('KTGG_?AAC?O?o?@?@?E?@?@??')
            sage: list([_ for _ in H.orthocomplementations_iterator()])
            []
        """
        n = self.order()

        # Special cases first
        if n == 0:
            yield []
            raise(StopIteration)
        if n == 1:
            yield [0]
            raise(StopIteration)
        if n % 2 == 1:
            raise(StopIteration)

        dual_isomorphism = self.is_isomorphic(self.reverse(), certificate=True)[1]
        if dual_isomorphism is None:  # i.e. if the lattice is not self-dual.
            raise(StopIteration)

        # We compute possible orthocomplements, i.e. elements
        # with "dual position" and complement to each other.

        orbits = self.automorphism_group(return_group=False, orbits=True)

        orbit_number = [None] * n
        for ind, orbit in enumerate(orbits):
            for e in orbit:
                orbit_number[e] = ind

        comps = [None] * n
        for e in range(n):
            # Fix following after ticket #20727
            comps[e] = [x for x in range(n) if
                        self._meet[e, x] == 0 and self._join[e, x] == n-1 and
                        x in orbits[orbit_number[dual_isomorphism[e]]]]

        # Fitting is done by this recursive function:
        def recursive_fit(orthocomplements, unbinded):
            if not unbinded:
                yield orthocomplements
            else:
                next_to_fit = unbinded[0]
                possible_values = [x for x in comps[next_to_fit] if not x in orthocomplements]
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

                    for i_want_python3_yield_from in recursive_fit(new_binded, new_unbinded):
                        yield i_want_python3_yield_from

        start = [None] * n
        # A little optimization
        for e in range(n):
            if len(comps[e]) == 0:  # Not any possible orthocomplement
                raise(StopIteration)
            if len(comps[e]) == 1:  # Do not re-fit this every time
                e_ = comps[e][0]
                # Every element might have one possible orthocomplement,
                # but so that they don't fit together. Must check that.
                for lc in self.lower_covers_iterator(e):
                    if start[lc] is not None:
                        if not self.has_edge(e_, start[lc]):
                            raise(StopIteration)
                if start[e_] is None:
                    start[e] = e_
                    start[e_] = e
        start_unbinded = [e for e in range(n) if start[e] is None]

        for i_want_python3_yield_from in recursive_fit(start, start_unbinded):
            yield i_want_python3_yield_from

    def find_nonsemimodular_pair(self, upper):
        """
        Return pair of elements showing the lattice is not modular.

        INPUT:

        - upper, a Boolean -- if ``True``, test wheter the lattice is
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

            sage: P = posets.PentagonPoset()
            sage: H = P._hasse_diagram
            sage: H.antichains_iterator()
            <generator object antichains_iterator at ...>
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
            sage: list(H.antichains_iterator())
            [[]]
        """
        # Complexity note:
        # antichains_queues never grows longer than self.cardinality().
        # Indeed, if a appears before b in antichains_queues, then
        # the largest element of a is strictly smaller than that of b.
        antichains_queues = [([], list(range(self.cardinality()-1, -1, -1)))]
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
                new_queue = [t for t in queue if not (leq[t,x] or leq[x,t])]
                antichains_queues.append((new_antichain, new_queue))

    def are_incomparable(self, i, j):
        """
        Returns whether ``i`` and ``j`` are incomparable in the poset

        INPUT:

         - ``i``, ``j`` -- vertices of this Hasse diagram

        EXAMPLES::

            sage: P = posets.PentagonPoset()
            sage: H = P._hasse_diagram
            sage: H.are_incomparable(1,2)
            True
            sage: [ (i,j) for i in H.vertices() for j in H.vertices() if H.are_incomparable(i,j)]
            [(1, 2), (1, 3), (2, 1), (3, 1)]
        """
        mat = self._leq_matrix
        return not mat[i,j] and not mat[j,i]

    def are_comparable(self, i, j):
        """
        Returns whether ``i`` and ``j`` are comparable in the poset

        INPUT:

         - ``i``, ``j`` -- vertices of this Hasse diagram

        EXAMPLES::

            sage: P = posets.PentagonPoset()
            sage: H = P._hasse_diagram
            sage: H.are_comparable(1,2)
            False
            sage: [ (i,j) for i in H.vertices() for j in H.vertices() if H.are_comparable(i,j)]
            [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (1, 0), (1, 1), (1, 4), (2, 0), (2, 2), (2, 3), (2, 4), (3, 0), (3, 2), (3, 3), (3, 4), (4, 0), (4, 1), (4, 2), (4, 3), (4, 4)]
        """
        mat = self._leq_matrix
        return bool(mat[i,j]) or bool(mat[j,i])

    def antichains(self, element_class = list):
        """
        Returns all antichains of ``self``, organized as a
        prefix tree

        INPUT:

         - ``element_class`` -- (default:list) an iterable type

        EXAMPLES::

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

            sage: TestSuite(A).run(skip = "_test_pickling")

        .. note:: It's actually the pickling of the cached method
            :meth:`coxeter_transformation` that fails ...

        TESTS::

            sage: A = Poset()._hasse_diagram.antichains()
            sage: list(A)
            [[]]
            sage: TestSuite(A).run()
        """
        from sage.combinat.subsets_pairwise import PairwiseCompatibleSubsets
        return PairwiseCompatibleSubsets(self.vertices(),
                                         self.are_incomparable,
                                         element_class = element_class)

    def chains(self, element_class=list, exclude=None):
        """
        Return all chains of ``self``, organized as a prefix tree.

        INPUT:

        - ``element_class`` -- (default: ``list``) an iterable type

        - ``exclude`` -- elements of the poset to be excluded
          (default: ``None``)

        OUTPUT:

        The enumerated set (with a forest structure given by prefix
        ordering) consisting of all chains of ``self``, each of
        which is given as an ``element_class``.

        EXAMPLES::

            sage: P = posets.PentagonPoset()
            sage: H = P._hasse_diagram
            sage: A = H.chains()
            sage: list(A)
            [[], [0], [0, 1], [0, 1, 4], [0, 2], [0, 2, 3], [0, 2, 3, 4], [0, 2, 4], [0, 3], [0, 3, 4], [0, 4], [1], [1, 4], [2], [2, 3], [2, 3, 4], [2, 4], [3], [3, 4], [4]]
            sage: A.cardinality()
            20
            sage: [1,3] in A
            False
            sage: [1,4] in A
            True

        One can exclude some vertices::

            sage: list(H.chains(exclude=[4, 3]))
            [[], [0], [0, 1], [0, 2], [1], [2]]

        The ``element_class`` keyword determines how the chains are
        being returned:

            sage: P = Poset({1: [2, 3], 2: [4]})
            sage: list(P._hasse_diagram.chains(element_class=tuple))
            [(), (0,), (0, 1), (0, 1, 2), (0, 2), (0, 3), (1,), (1, 2), (2,), (3,)]
            sage: list(P._hasse_diagram.chains())
            [[], [0], [0, 1], [0, 1, 2], [0, 2], [0, 3], [1], [1, 2], [2], [3]]

        (Note that taking the Hasse diagram has renamed the vertices.)

            sage: list(P._hasse_diagram.chains(element_class=tuple, exclude=[0]))
            [(), (1,), (1, 2), (2,), (3,)]

        .. SEEALSO:: :meth:`antichains`
        """
        from sage.combinat.subsets_pairwise import PairwiseCompatibleSubsets
        if not(exclude is None):
            vertices = [u for u in self.vertices() if not u in exclude]
        else:
            vertices = self.vertices()
        return PairwiseCompatibleSubsets(vertices,
                                         self.are_comparable,
                                         element_class = element_class)

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
            <generator object sublattices_iterator at ...>
            sage: next(it)
            set()
            sage: next(it)
            {0}
        """
        # Python3-note: "yield from" would be simpler.
        yield elms
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
                    gens.add(self._meet[x, g])
                    gens.add(self._join[x, g])
                current_set.add(g)
            else:
                for x in self.sublattices_iterator(current_set, e+1):
                    yield x

    def maximal_sublattices(self):
        """
        Return maximal sublattices of the lattice.

        EXAMPLES::

            sage: L = Posets.PentagonPoset()
            sage: ms = L._hasse_diagram.maximal_sublattices()
            sage: sorted(ms, key=sorted)
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

            sage: H = Posets.PentagonPoset()._hasse_diagram
            sage: H.frattini_sublattice()
            [0, 4]
        """
        # Just a direct computation, no optimization at all.
        n = self.cardinality()
        if n == 0 or n == 2: return []
        if n == 1: return [0]
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
        tmp = list(self.depth_first_search(uc, neighbors=lambda v: [v_ for v_ in self.neighbors_in(v) if v_ not in lt_a]))
        result = None
        for e in tmp:
            if all(x not in tmp for x in self.neighbors_in(e)):
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
            sage: H.skeleton()
            [5, 2, 0, 3]
        """
        p_atoms = []
        for atom in self.neighbor_out_iterator(0):
            p_atom = self.pseudocomplement(atom)
            if p_atom is None:
                raise ValueError("lattice is not pseudocomplemented")
            p_atoms.append(p_atom)
        n = len(p_atoms)
        mt = self._meet
        pos = [0] * n
        meets = [self.order()-1] * n
        result = [self.order()-1]
        i = 0

        while i >= 0:
            new_meet = mt[meets[i-1], p_atoms[pos[i]]]
            result.append(new_meet)
            if pos[i] == n-1:
                i -= 1
                pos[i] = pos[i]+1
            else:
                meets[i] = new_meet
                pos[i+1] = pos[i]+1
                i += 1

        return result

    def is_convex_subset(self, S):
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
                neighbors = lambda v_: [v for v in self.neighbor_out_iterator(v_)
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
            sage: H= HasseDiagram({0: [1, 2], 1: [4], 2: [3], 3: [4, 5],
            ....:                  4: [6], 5:[6]})
            sage: sorted(H.neutral_elements())
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

        todo = set(range(1, n-1))
        neutrals = set([0, n-1])
        notneutrals = set()
        all_elements = set(range(n))

        mt = self._meet
        jn = self._join

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
        tmp = list(self.depth_first_search(lc, neighbors=lambda v: [v_ for v_ in self.neighbors_out(v) if v_ not in gt_a]))
        result = None
        for e in tmp:
            if all(x not in tmp for x in self.neighbors_out(e)):
                if result:
                    return None
                result = e
        return result

from sage.misc.rest_index_of_methods import gen_rest_table_index
__doc__ = __doc__.format(INDEX_OF_FUNCTIONS=gen_rest_table_index(HasseDiagram))
