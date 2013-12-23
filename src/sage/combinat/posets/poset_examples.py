r"""
A collection of posets and lattices.
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

import random
from sage.misc.classcall_metaclass import ClasscallMetaclass
import sage.categories.posets
from sage.combinat.permutation import Permutations, Permutation
from sage.combinat.posets.posets import Poset, FinitePosets_n
from sage.combinat.posets.lattices import LatticePoset
from sage.graphs.digraph import DiGraph
from sage.rings.integer import Integer

class Posets(object):
    r"""
    A collection of posets and lattices.

    EXAMPLES::

        sage: Posets.BooleanLattice(3)
        Finite lattice containing 8 elements
        sage: Posets.ChainPoset(3)
        Finite lattice containing 3 elements
        sage: Posets.RandomPoset(17,.15)
        Finite poset containing 17 elements

    The category of all posets::

        sage: Posets()
        Category of posets

    The enumerated set of all posets on `3` vertices, up to an
    isomorphism::

        sage: Posets(3)
        Posets containing 3 vertices

    .. seealso:: :class:`~sage.categories.posets.Posets`, :class:`FinitePosets`, :func:`Poset`

    TESTS::

        sage: P = Posets
        sage: TestSuite(P).run()
    """

    __metaclass__ = ClasscallMetaclass
    @staticmethod
    def __classcall__(cls, n = None):
        r"""
        Return either the category of all posets, or the finite
        enumerated set of all finite posets on ``n`` elements up to an
        isomorphism.

        EXAMPLES::

            sage: Posets()
            Category of posets
            sage: Posets(4)
            Posets containing 4 vertices
        """
        if n is None:
            return sage.categories.posets.Posets()
        return FinitePosets_n(n)

    @staticmethod
    def BooleanLattice(n):
        """
        Returns the Boolean lattice containing `2^n` elements.

        EXAMPLES::

            sage: Posets.BooleanLattice(5)
            Finite lattice containing 32 elements
        """
        return LatticePoset([[Integer(x|(1<<y)) for y in range(0,n) if x&(1<<y)==0] for
            x in range(0,2**n)])

    @staticmethod
    def ChainPoset(n):
        """
        Returns a chain (a totally ordered poset) containing ``n`` elements.

        EXAMPLES::

            sage: C = Posets.ChainPoset(6); C
            Finite lattice containing 6 elements
            sage: C.linear_extension()
            [0, 1, 2, 3, 4, 5]
            sage: for i in range(5):
            ...       for j in range(5):
            ...           if C.covers(C(i),C(j)) and j != i+1:
            ...              print "TEST FAILED"

        TESTS:

        Check that #8422 is solved::

            sage: Posets.ChainPoset(0)
            Finite lattice containing 0 elements
            sage: C = Posets.ChainPoset(1); C
            Finite lattice containing 1 elements
            sage: C.cover_relations()
            []
            sage: C = Posets.ChainPoset(2); C
            Finite lattice containing 2 elements
            sage: C.cover_relations()
            [[0, 1]]
        """
        return LatticePoset((range(n), [[x,x+1] for x in range(n-1)]))

    @staticmethod
    def AntichainPoset(n):
        """
        Returns an antichain (a poset with no comparable elements)
        containing `n` elements.

        EXAMPLES::

            sage: A = Posets.AntichainPoset(6); A
            Finite poset containing 6 elements
            sage: for i in range(5):
            ...       for j in range(5):
            ...           if A.covers(A(i),A(j)):
            ...              print "TEST FAILED"

        TESTS:

        Check that #8422 is solved::

            sage: Posets.AntichainPoset(0)
            Finite poset containing 0 elements
            sage: C = Posets.AntichainPoset(1); C
            Finite poset containing 1 elements
            sage: C.cover_relations()
            []
            sage: C = Posets.AntichainPoset(2); C
            Finite poset containing 2 elements
            sage: C.cover_relations()
            []
        """
        return Poset((range(n), []))

    @staticmethod
    def PentagonPoset(facade = False):
        """
        Returns the "pentagon poset".

        EXAMPLES::

            sage: P = Posets.PentagonPoset(); P
            Finite lattice containing 5 elements
            sage: P.cover_relations()
            [[0, 1], [0, 2], [1, 4], [2, 3], [3, 4]]

        This lattice and the diamond poset on 5 elements are the two
        smallest lattices which are not distributive::

            sage: P.is_distributive()
            False
            sage: Posets.DiamondPoset(5).is_distributive()
            False
        """
        p = LatticePoset([[1,2],[4],[3],[4],[]], facade = facade)
        p.hasse_diagram()._pos = {0:[2,0],1:[0,2],2:[3,1],3:[3,3],4:[2,4]}
        return p

    @staticmethod
    def DiamondPoset(n, facade = False):
        """
        Returns the lattice of rank two containing ``n`` elements.

        EXAMPLES::

            sage: Posets.DiamondPoset(7)
            Finite lattice containing 7 elements
        """
        c = [[n-1] for x in range(n)]
        c[0] = [x for x in range(1,n-1)]
        c[n-1] = []
        if n > 2:
            return LatticePoset(c, facade = facade)
        else:
            return Poset(c, facade = facade)

    @staticmethod
    def IntegerCompositions(n):
        """
        Returns the poset of integer compositions of the integer ``n``.

        A composition of a positive integer `n` is a list of positive
        integers that sum to `n`. The order is reverse refinement:
        `[p_1,p_2,...,p_l] < [q_1,q_2,...,q_m]` if `q` consists
        of an integer composition of `p_1`, followed by an integer
        composition of `p_2`, and so on.

        EXAMPLES::

            sage: P = Posets.IntegerCompositions(7); P
            Finite poset containing 64 elements
            sage: len(P.cover_relations())
            192
        """
        from sage.combinat.composition import Compositions
        C = Compositions(n)
        return Poset((C, [[c,d] for c in C for d in C if d.is_finer(c)]), cover_relations=False)

    @staticmethod
    def IntegerPartitions(n):
        """
        Returns the poset of integer partitions on the integer ``n``.

        A partition of a positive integer `n` is a non-increasing list
        of positive integers that sum to `n`. If `p` and `q` are
        integer partitions of `n`, then `p` covers `q` if and only
        if `q` is obtained from `p` by joining two parts of `p`
        (and sorting, if necessary).

        EXAMPLES::

            sage: P = Posets.IntegerPartitions(7); P
            Finite poset containing 15 elements
            sage: len(P.cover_relations())
            28
        """
        def lower_covers(partition):
            r"""
            Nested function for computing the lower covers
            of elements in the poset of integer partitions.
            """
            lc = []
            for i in range(0,len(partition)-1):
                for j in range(i+1,len(partition)):
                    new_partition = partition[:]
                    del new_partition[j]
                    del new_partition[i]
                    new_partition.append(partition[i]+partition[j])
                    new_partition.sort(reverse=True)
                    tup = tuple(new_partition)
                    if tup not in lc:
                        lc.append(tup)
            return lc
        from sage.combinat.partition import Partitions
        H = DiGraph(dict([[tuple(p),lower_covers(p)] for p in Partitions(n)]))
        return Poset(H.reverse())

    @staticmethod
    def RestrictedIntegerPartitions(n):
        """
        Returns the poset of integer partitions on the integer `n`
        ordered by restricted refinement. That is, if `p` and `q`
        are integer partitions of `n`, then `p` covers `q` if and
        only if `q` is obtained from `p` by joining two distinct
        parts of `p` (and sorting, if necessary).

        EXAMPLES::

            sage: P = Posets.RestrictedIntegerPartitions(7); P
            Finite poset containing 15 elements
            sage: len(P.cover_relations())
            17
        """
        def lower_covers(partition):
            r"""
            Nested function for computing the lower covers of elements in the
            restricted poset of integer partitions.
            """
            lc = []
            for i in range(0,len(partition)-1):
                for j in range(i+1,len(partition)):
                    if partition[i] != partition[j]:
                        new_partition = partition[:]
                        del new_partition[j]
                        del new_partition[i]
                        new_partition.append(partition[i]+partition[j])
                        new_partition.sort(reverse=True)
                        tup = tuple(new_partition)
                        if tup not in lc:
                            lc.append(tup)
            return lc
        from sage.combinat.partition import Partitions
        H = DiGraph(dict([[tuple(p),lower_covers(p)] for p in Partitions(n)]))
        return Poset(H.reverse())

    @staticmethod
    def RandomPoset(n,p):
        r"""
        Generate a random poset on ``n`` vertices according to a
        probability ``p``.

        INPUT:

        - ``n`` - number of vertices, a non-negative integer

        - ``p`` - a probability, a real number between 0 and 1 (inclusive)

        OUTPUT:

        A poset on ``n`` vertices.  The construction decides to make an
        ordered pair of vertices comparable in the poset with probability
        ``p``, however a pair is not made comparable if it would violate
        the defining properties of a poset, such as transitivity.

        So in practice, once the probability exceeds a small number the
        generated posets may be very similar to a chain.  So to create
        interesting examples, keep the probability small, perhaps on the
        order of `1/n`.

        EXAMPLES::

            sage: Posets.RandomPoset(17,.15)
            Finite poset containing 17 elements

        TESTS::

            sage: Posets.RandomPoset('junk', 0.5)
            Traceback (most recent call last):
            ...
            TypeError: number of elements must be an integer, not junk

            sage: Posets.RandomPoset(-6, 0.5)
            Traceback (most recent call last):
            ...
            ValueError: number of elements must be non-negative, not -6

            sage: Posets.RandomPoset(6, 'garbage')
            Traceback (most recent call last):
            ...
            TypeError: probability must be a real number, not garbage

            sage: Posets.RandomPoset(6, -0.5)
            Traceback (most recent call last):
            ...
            ValueError: probability must be between 0 and 1, not -0.5
        """
        try:
            n = Integer(n)
        except TypeError:
            raise TypeError("number of elements must be an integer, not {0}".format(n))
        if n < 0:
            raise ValueError("number of elements must be non-negative, not {0}".format(n))
        try:
            p = float(p)
        except StandardError:
            raise TypeError("probability must be a real number, not {0}".format(p))
        if p < 0 or p> 1:
            raise ValueError("probability must be between 0 and 1, not {0}".format(p))

        D = DiGraph(loops=False,multiedges=False)
        D.add_vertices(range(n))
        for i in range(n):
            for j in range(n):
                if random.random() < p:
                    D.add_edge(i,j)
                    if not D.is_directed_acyclic():
                        D.delete_edge(i,j)
        return Poset(D,cover_relations=False)

    @staticmethod
    def SSTPoset(s,f=None):
        """
        The poset on semistandard tableaux of shape s and largest entry f that is ordered by componentwise comparison of the entries.

        INPUT:

        - ``s`` - shape of the tableaux

        - ``f`` - maximum fill number.  This is an optional argument.  If no maximal number is given, it will use the number of cells in the shape.

        NOTE: This is basic implementation and most certainly not the most efficient.

        EXAMPLES::

            sage: Posets.SSTPoset([2,1])
            Finite poset containing 8 elements

            sage: Posets.SSTPoset([2,1],4)
            Finite poset containing 20 elements

            sage: Posets.SSTPoset([2,1],2).cover_relations()
            [[[[1, 1], [2]], [[1, 2], [2]]]]

            sage: Posets.SSTPoset([3,2]).bottom()  # long time (6s on sage.math, 2012)
            [[1, 1, 1], [2, 2]]

            sage: Posets.SSTPoset([3,2],4).maximal_elements()
            [[[3, 3, 4], [4, 4]]]
        """
        from sage.combinat.tableau import SemistandardTableaux
        def tableaux_is_less_than(a,b):
            atstring = []
            btstring = []
            c=0
            for i in range(len(a)):
                atstring=atstring+a[i]
            for i in range(len(b)):
                btstring=btstring+b[i]
            for i in range(len(atstring)):
                if atstring[i] > btstring[i]:
                    c = c+1
            if c == 0:
                return True
            else:
                return False
        if f is None:
            f=0
            for i in range(len(s)):
                f = f+s[i]
        E = SemistandardTableaux(s,max_entry=f)
        return Poset((E, tableaux_is_less_than ))

    @staticmethod
    def SymmetricGroupBruhatOrderPoset(n):
        """
        The poset of permutations with respect to Bruhat order.

        EXAMPLES::

            sage: Posets.SymmetricGroupBruhatOrderPoset(4)
            Finite poset containing 24 elements
        """
        if n < 10:
            element_labels = dict([[s,"".join(map(str,s))] for s in Permutations(n)])
        return Poset(dict([[s,s.bruhat_succ()]
                for s in Permutations(n)]),element_labels)

    @staticmethod
    def SymmetricGroupBruhatIntervalPoset(start, end):
        """
        The poset of permutations with respect to Bruhat order.

        INPUT:

        - ``start`` - list permutation

        - ``end`` - list permutation (same n, of course)

        .. note::

           Must have ``start`` <= ``end``.

        EXAMPLES:

        Any interval is rank symmetric if and only if it avoids these
        permutations::

            sage: P1 = Posets.SymmetricGroupBruhatIntervalPoset([1,2,3,4], [3,4,1,2])
            sage: P2 = Posets.SymmetricGroupBruhatIntervalPoset([1,2,3,4], [4,2,3,1])
            sage: ranks1 = [P1.rank(v) for v in P1]
            sage: ranks2 = [P2.rank(v) for v in P2]
            sage: [ranks1.count(i) for i in uniq(ranks1)]
            [1, 3, 5, 4, 1]
            sage: [ranks2.count(i) for i in uniq(ranks2)]
            [1, 3, 5, 6, 4, 1]

        """
        start = Permutation(start)
        end = Permutation(end)
        if len(start) != len(end):
            raise TypeError("Start (%s) and end (%s) must have same length."%(start, end))
        if not start.bruhat_lequal(end):
            raise TypeError("Must have start (%s) <= end (%s) in Bruhat order."%(start, end))
        unseen = [start]
        nodes = {}
        while len(unseen) > 0:
            perm = unseen.pop(0)
            nodes[perm] = [succ_perm for succ_perm in perm.bruhat_succ()
                                if succ_perm.bruhat_lequal(end)]
            for succ_perm in nodes[perm]:
                if not nodes.has_key(succ_perm):
                    unseen.append(succ_perm)
        return Poset(nodes)

    @staticmethod
    def SymmetricGroupWeakOrderPoset(n, labels="permutations", side="right"):
        r"""
        The poset of permutations of `\{ 1, 2, \ldots, n \}` with respect
        to the weak order (also known as the permutohedron order, cf.
        :meth:`~sage.combinat.permutation.Permutation.permutohedron_lequal`).

        The optional variable ``labels`` (default: ``"permutations"``)
        determines the labelling of the elements if `n < 10`. The optional
        variable ``side`` (default: ``"right"``) determines whether the
        right or the left permutohedron order is to be used.

        EXAMPLES::

            sage: Posets.SymmetricGroupWeakOrderPoset(4)
            Finite poset containing 24 elements
        """
        if n < 10 and labels == "permutations":
            element_labels = dict([[s,"".join(map(str,s))] for s in Permutations(n)])
        if n < 10 and labels == "reduced_words":
            element_labels = dict([[s,"".join(map(str,s.reduced_word_lexmin()))] for s in Permutations(n)])
        if side == "left":
            def weak_covers(s):
                r"""
                Nested function for computing the covers of elements in the
                poset of left weak order for the symmetric group.
                """
                return [v for v in s.bruhat_succ() if
                    s.length() + (s.inverse().right_action_product(v)).length() == v.length()]
        else:
            def weak_covers(s):
                r"""
                Nested function for computing the covers of elements in the
                poset of right weak order for the symmetric group.
                """
                return [v for v in s.bruhat_succ() if
                    s.length() + (s.inverse().left_action_product(v)).length() == v.length()]
        return Poset(dict([[s,weak_covers(s)] for s in Permutations(n)]),element_labels)

posets = Posets
