r"""
Hasse diagrams of posets
"""
#*****************************************************************************
#       Copyright (C) 2008 Peter Jipsen <jipsen@chapman.edu>,
#                          Franco Saliola <saliola@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from copy import copy
from sage.graphs.digraph import DiGraph
from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ
from sage.misc.misc import uniq
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.cachefunc import cached_method

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
        TESTS::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,2],1:[3],2:[3],3:[]})
            sage: H.linear_extension()
            [0, 1, 2, 3]
        """
        # Recall: we assume range(n) is a linear extension.
        return range(len(self))

    def linear_extensions(self):
        r"""
        TESTS::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,2],1:[3],2:[3],3:[]})
            sage: H.linear_extensions()
            [[0, 1, 2, 3], [0, 2, 1, 3]]
        """
        return self.topological_sort_generator()

    def is_linear_extension(self,lin_ext=None):
        r"""
        TESTS::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,2],1:[3],2:[3],3:[]})
            sage: H.is_linear_extension(range(4))
            True
            sage: H.is_linear_extension([3,2,1,0])
            False
        """
        if lin_ext is None or lin_ext == range(len(self)):
            for x,y in self.cover_relations_iterator():
                if not x < y:
                    return False
            return True
        else:
            for x,y in self.cover_relations_iterator():
                if not lin_ext.index(x) < lin_ext.index(y):
                    return False
            return True

    def cover_relations_iterator(self):
        r"""
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
        TESTS::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[2,3], 1:[3,4], 2:[5], 3:[5], 4:[5]})
            sage: H.cover_relations()
            [(0, 2), (0, 3), (1, 3), (1, 4), (2, 5), (3, 5), (4, 5)]
        """
        return [c for c in self.cover_relations_iterator()]

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
        return (self.num_edges()+1 == self.num_verts() and # Hasse Diagram is a tree
                all(d<=1 for d in self.out_degree())   and # max outdegree is <= 1
                all(d<=1 for d in self.in_degree()))       # max  indegree is <= 1

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
        H.relabel(perm=range(H.num_verts()-1,-1,-1), inplace=True)
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
        return [z for z in range(self.order())[x:y+1] if
                self.is_lequal(x,z) and self.is_lequal(z,y)]

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
            ...    if r(v) != r(u) + 1:
            ...        print "Bug in rank_function!"

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

    def is_graded(self):
        r"""
        Deprecated, has conflicting definition of "graded" vs. "ranked"
        with posets.

        Return ``True`` if the Hasse diagram is ranked. For definition
        of ranked see :meth:`~rank_function`.
        """
        from sage.misc.superseded import deprecation
        deprecation(16998, "Use is_ranked(). Definition conflict with posets.")
        return self.is_ranked()

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

    def mobius_function(self,i,j): # dumb algorithm
        r"""
        Returns the value of the M\"obius function of the poset
        on the elements ``i`` and ``j``.

        EXAMPLES::

            sage: P = Poset([[1,2,3],[4],[4],[4],[]])
            sage: H = P._hasse_diagram
            sage: H.mobius_function(0,4)
            2
            sage: for u,v in P.cover_relations_iterator():
            ...    if P.mobius_function(u,v) != -1:
            ...        print "Bug in mobius_function!"
        """
        try:
            return self._mobius_function_values[(i,j)]
        except AttributeError:
            self._mobius_function_values = {}
            return self.mobius_function(i,j)
        except KeyError:
            if i == j:
                self._mobius_function_values[(i,j)] = 1
            elif i > j:
                self._mobius_function_values[(i,j)] = 0
            else:
                ci = self.closed_interval(i,j)
                if len(ci) == 0:
                    self._mobius_function_values[(i,j)] = 0
                else:
                    self._mobius_function_values[(i,j)] = \
                     -sum([self.mobius_function(i,k) for k in ci[:-1]])
        return self._mobius_function_values[(i,j)]

    def mobius_function_matrix(self):
        r"""
        Returns the matrix of the Mobius function of this poset

        This returns the sparse matrix over `\ZZ` whose ``(x, y)`` entry
        is the value of the M\"obius function of ``self`` evaluated on
        ``x`` and ``y``, and redefines :meth:`mobius_function` to use
        it.

        .. NOTE::

            The result is cached in :meth:`_mobius_function_matrix`.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,3,2],1:[4],2:[4,5,6],3:[6],4:[7],5:[7],6:[7],7:[]})
            sage: H.mobius_function_matrix()
            [ 1 -1 -1 -1  1  0  1  0]
            [ 0  1  0  0 -1  0  0  0]
            [ 0  0  1  0 -1 -1 -1  2]
            [ 0  0  0  1  0  0 -1  0]
            [ 0  0  0  0  1  0  0 -1]
            [ 0  0  0  0  0  1  0 -1]
            [ 0  0  0  0  0  0  1 -1]
            [ 0  0  0  0  0  0  0  1]

        TESTS::

            sage: H.mobius_function_matrix().is_immutable()
            True
            sage: hasattr(H,'_mobius_function_matrix')
            True

            sage: H.mobius_function == H._mobius_function_from_matrix
            True
        """
        if not hasattr(self,'_mobius_function_matrix'):
            self._mobius_function_matrix = self.lequal_matrix().inverse().change_ring(ZZ)
            self._mobius_function_matrix.set_immutable()
            self.mobius_function = self._mobius_function_from_matrix
        return self._mobius_function_matrix

    # Redefine self.mobius_function
    def _mobius_function_from_matrix(self, i,j):
        r"""
        Returns the value of the M\"obius function of the poset
        on the elements ``i`` and ``j``.

        EXAMPLES::

            sage: P = Poset([[1,2,3],[4],[4],[4],[]])
            sage: H = P._hasse_diagram
            sage: H.mobius_function(0,4) # indirect doctest
            2
            sage: for u,v in P.cover_relations_iterator():
            ...    if P.mobius_function(u,v) != -1:
            ...        print "Bug in mobius_function!"

        This uses ``self._mobius_function_matrix``, as computed by
        :meth:`mobius_function_matrix`.
        """
        return self._mobius_function_matrix[i,j]

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
        return - self.lequal_matrix()*self.mobius_function_matrix().transpose()

    def order_filter(self,elements):
        """
        Returns the order filter generated by a list of elements.

        `I` is an order filter if, for any `x` in `I` and `y` such that
        `y \ge x`, then `y` is in `I`.

        EXAMPLES::

            sage: H = Posets.BooleanLattice(4)._hasse_diagram
            sage: H.order_filter([3,8])
            [3, 7, 8, 9, 10, 11, 12, 13, 14, 15]
        """
        of = []
        for i in elements:
            for j in self.breadth_first_search(i):
                of.append(j)
        return uniq(of)

    def principal_order_filter(self, i):
        """
        Returns the order filter generated by ``i``.

        EXAMPLES::

            sage: H = Posets.BooleanLattice(4)._hasse_diagram
            sage: H.principal_order_filter(2)
            [2, 3, 6, 7, 10, 11, 14, 15]
        """
        return self.order_filter([i])

    def order_ideal(self,elements):
        """
        Returns the order ideal generated by a list of elements.

        `I` is an order ideal if, for any `x` in `I` and `y` such that
        `y \le x`, then `y` is in `I`.

        EXAMPLES::

            sage: H = Posets.BooleanLattice(4)._hasse_diagram
            sage: H.order_ideal([7,10])
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 10]
        """
        H = copy(self).reverse()
        oi = []
        for i in elements:
            for j in H.breadth_first_search(i):
                oi.append(j)
        return uniq(oi)

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

    @lazy_attribute
    def _meet(self):
        r"""
        Computes the matrix of meets of ``self``. The ``(x,y)``-entry of
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
            ValueError: Not a meet-semilattice: no bottom element.

            sage: H = HasseDiagram({0:[1,2],1:[3,4],2:[3,4]})
            sage: H.meet_matrix()
            Traceback (most recent call last):
            ...
            ValueError: No meet for x=...

            sage: L = LatticePoset({0:[1,2,3],1:[4],2:[4],3:[4]})
            sage: P = L.dual()
            sage: P.meet(2,3)
            4
        """
        n = self.cardinality()
        if n == 0:
            return matrix(0)
        if not self.has_bottom():
            raise ValueError("Not a meet-semilattice: no bottom element.")
        le = self._leq_matrix
        meet = [[0 for x in range(n)] for x in range(n)]
        lc = [self.neighbors_in(x) for x in range(n)]

        for x in range(n): # x=x_k
            meet[x][x] = x
            for y in range(x):
                T = [meet[y][z] for z in lc[x]] # T = {x_i \wedge z : z>-x_k}

                q = max(T)
                for z in T:
                    if not le[z,q]:
                        raise ValueError("No meet for x=%s y=%s"%(x,y))
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
            ValueError: Not a meet-semilattice: no bottom element.

            sage: H = HasseDiagram({0:[1,2],1:[3,4],2:[3,4]})
            sage: H.meet_matrix()
            Traceback (most recent call last):
            ...
            ValueError: No meet for x=...
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
            ValueError: Not a join-semilattice: no top element.

            sage: H = HasseDiagram({0:[2,3],1:[2,3],2:[4],3:[4]})
            sage: H.join_matrix()
            Traceback (most recent call last):
            ...
            ValueError: No join for x=...

            sage: L = LatticePoset({0:[1,2,3],1:[4],2:[4],3:[4]})
            sage: P = L.dual()
            sage: P.join(2,3)
            0
        """
        n = self.cardinality()
        if n == 0:
            return matrix(0)
        if not self.has_top():
            raise ValueError("Not a join-semilattice: no top element.")
        join = [[0 for x in range(n)] for x in range(n)]
        le = self.lequal_matrix()
        uc = [sorted([n-1-y for y in self.neighbors_out(x)]) for
                x in reversed(range(n))]

        for x in range(n): # x=x_k
            join[x][x] = x

            for y in range(x):
                T = [join[y][z] for z in uc[x]]

                q = max(T)
                for z in T:
                    if not le[n-1-q, n-1-z]:
                        raise ValueError("No join for x=%s y=%s"%(x,y))
                join[x][y] = q
                join[y][x] = q

        return matrix(ZZ, [[n-1-join[n-1-x][n-1-y] for y in range(n)]
                           for x in range(n)])

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
            ValueError: Not a join-semilattice: no top element.

            sage: H = HasseDiagram({0:[2,3],1:[2,3],2:[4],3:[4]})
            sage: H.join_matrix()
            Traceback (most recent call last):
            ...
            ValueError: No join for x=...
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

    def is_complemented_lattice(self):
        r"""
        Returns ``True`` if ``self`` is the Hasse diagram of a
        complemented lattice, and ``False`` otherwise.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,2,3],1:[4],2:[4],3:[4]})
            sage: H.is_complemented_lattice()
            True

            sage: H = HasseDiagram({0:[1,2],1:[3],2:[3],3:[4]})
            sage: H.is_complemented_lattice()
            False
        """
        try:
            jn = self.join_matrix()
            mt = self.meet_matrix()
        except ValueError:
            return False
        n = self.cardinality()
        c = [-1 for x in range(n)]
        for x in range(n):
            for y in range(x,n):
                if jn[x][y]==n-1 and mt[x][y]==0:
                    c[x]=y
                    c[y]=x
        return all([c[x]!=-1 for x in range(n)])

    def complements(self):
        r"""
        Deprecated.
        """
        from sage.misc.superseded import deprecation
        deprecation(17138, "This function is broken. Do not use.")
        jn = self.join_matrix()
        mt = self.meet_matrix()
        n = self.cardinality()
        c = [None for x in range(n)]
        for x in range(n):
            for y in range(x,n):
                if jn[x][y]==n-1 and mt[x][y]==0:
                    c[x]=y
                    c[y]=x
        return c

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
        antichains_queues = [([], range(self.cardinality()-1,-1,-1))]
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

        .. seealso:: :meth:`antichains`
        """
        from sage.combinat.subsets_pairwise import PairwiseCompatibleSubsets
        if not(exclude is None):
            vertices = [u for u in self.vertices() if not u in exclude]
        else:
            vertices = self.vertices()
        return PairwiseCompatibleSubsets(vertices,
                                         self.are_comparable,
                                         element_class = element_class)
