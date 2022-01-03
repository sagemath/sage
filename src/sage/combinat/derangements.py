"""
Derangements

AUTHORS:

- Alasdair McAndrew (2010-05): Initial version
- Travis Scrimshaw (2013-03-30): Put derangements into category framework
"""

# ****************************************************************************
#       Copyright (C) 2010 Alasdair McAndrew <amca01@gmail.com>,
#                     2013 Travis Scrimshaw <tscrim@ucdavis.edu>
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
# ****************************************************************************

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.misc.misc_c import prod
from sage.misc.prandom import random, randrange
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.integer import Integer
from sage.combinat.combinat import CombinatorialElement
from sage.combinat.permutation import Permutation, Permutations


class Derangement(CombinatorialElement):
    r"""
    A derangement.

    A derangement on a set `S` is a permutation `\sigma` such that `\sigma(x)
    \neq x` for all `x \in S`, i.e. `\sigma` is a permutation of `S` with no
    fixed points.

    EXAMPLES::

        sage: D = Derangements(4)
        sage: elt = D([4,3,2,1])
        sage: TestSuite(elt).run()
    """
    def to_permutation(self):
        """
        Return the permutation corresponding to ``self``.

        EXAMPLES::

            sage: D = Derangements(4)
            sage: p = D([4,3,2,1]).to_permutation(); p
            [4, 3, 2, 1]
            sage: type(p)
            <class 'sage.combinat.permutation.StandardPermutations_all_with_category.element_class'>
            sage: D = Derangements([1, 3, 3, 4])
            sage: D[0].to_permutation()
            Traceback (most recent call last):
            ...
            ValueError: Can only convert to a permutation for derangements of [1, 2, ..., n]
        """
        if self.parent()._set != tuple(range(1, len(self) + 1)):
            raise ValueError("Can only convert to a permutation for derangements of [1, 2, ..., n]")
        return Permutation(list(self))


class Derangements(UniqueRepresentation, Parent):
    r"""
    The class of all derangements of a set or multiset.

    A derangement on a set `S` is a permutation `\sigma` such that `\sigma(x)
    \neq x` for all `x \in S`, i.e. `\sigma` is a permutation of `S` with no
    fixed points.

    For an integer, or a list or string with all elements
    distinct, the derangements are obtained by a standard result described
    in [BV2004]_. For a list or string with repeated elements, the derangements
    are formed by computing all permutations of the input and discarding all
    non-derangements.

    INPUT:

    - ``x`` -- Can be an integer which corresponds to derangements of
      `\{1, 2, 3, \ldots, x\}`, a list, or a string

    REFERENCES:

    - [BV2004]_
    - :wikipedia:`Derangement`

    EXAMPLES::

        sage: D1 = Derangements([2,3,4,5])
        sage: D1.list()
        [[3, 4, 5, 2],
         [5, 4, 2, 3],
         [3, 5, 2, 4],
         [4, 5, 3, 2],
         [4, 2, 5, 3],
         [5, 2, 3, 4],
         [5, 4, 3, 2],
         [4, 5, 2, 3],
         [3, 2, 5, 4]]
        sage: D1.cardinality()
        9
        sage: D1.random_element() # random
        [4, 2, 5, 3]
        sage: D2 = Derangements([1,2,3,1,2,3])
        sage: D2.cardinality()
        10
        sage: D2.list()
        [[2, 1, 1, 3, 3, 2],
         [2, 1, 2, 3, 3, 1],
         [2, 3, 1, 2, 3, 1],
         [2, 3, 1, 3, 1, 2],
         [2, 3, 2, 3, 1, 1],
         [3, 1, 1, 2, 3, 2],
         [3, 1, 2, 2, 3, 1],
         [3, 1, 2, 3, 1, 2],
         [3, 3, 1, 2, 1, 2],
         [3, 3, 2, 2, 1, 1]]
        sage: D2.random_element() # random
        [2, 3, 1, 3, 1, 2]
    """
    @staticmethod
    def __classcall_private__(cls, x):
        """
        Normalize ``x`` to ensure a unique representation.

        EXAMPLES::

            sage: D = Derangements(4)
            sage: D2 = Derangements([1, 2, 3, 4])
            sage: D3 = Derangements((1, 2, 3, 4))
            sage: D is D2
            True
            sage: D is D3
            True
        """
        if x in ZZ:
            x = list(range(1, x + 1))
        return super(Derangements, cls).__classcall__(cls, tuple(x))

    def __init__(self, x):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: D = Derangements(4)
            sage: TestSuite(D).run()
            sage: D = Derangements('abcd')
            sage: TestSuite(D).run()
            sage: D = Derangements([2, 2, 1, 1])
            sage: TestSuite(D).run()
        """
        Parent.__init__(self, category=FiniteEnumeratedSets())
        self._set = x
        self.__multi = len(set(x)) < len(x)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: Derangements(4)
            Derangements of the set [1, 2, 3, 4]
            sage: Derangements('abcd')
            Derangements of the set ['a', 'b', 'c', 'd']
            sage: Derangements([2,2,1,1])
            Derangements of the multiset [2, 2, 1, 1]
        """
        if self.__multi:
            return "Derangements of the multiset %s" % list(self._set)
        return "Derangements of the set %s" % list(self._set)

    def _element_constructor_(self, der):
        """
        Construct an element of ``self`` from ``der``.

        EXAMPLES::

            sage: D = Derangements(4)
            sage: elt = D([3,1,4,2]); elt
            [3, 1, 4, 2]
            sage: elt.parent() is D
            True
        """
        if isinstance(der, Derangement):
            if der.parent() is self:
                return der
            raise ValueError("Cannot convert %s to an element of %s" % (der, self))
        return self.element_class(self, der)

    Element = Derangement

    def __iter__(self):
        """
        Iterate through ``self``.

        EXAMPLES::

            sage: D = Derangements(4)
            sage: D.list() # indirect doctest
            [[2, 3, 4, 1],
             [4, 3, 1, 2],
             [2, 4, 1, 3],
             [3, 4, 2, 1],
             [3, 1, 4, 2],
             [4, 1, 2, 3],
             [4, 3, 2, 1],
             [3, 4, 1, 2],
             [2, 1, 4, 3]]
            sage: D = Derangements([1,44,918,67])
            sage: D.list()
            [[44, 918, 67, 1],
             [67, 918, 1, 44],
             [44, 67, 1, 918],
             [918, 67, 44, 1],
             [918, 1, 67, 44],
             [67, 1, 44, 918],
             [67, 918, 44, 1],
             [918, 67, 1, 44],
             [44, 1, 67, 918]]
            sage: D = Derangements(['A','AT','CAT','CATS'])
            sage: D.list()
            [['AT', 'CAT', 'CATS', 'A'],
             ['CATS', 'CAT', 'A', 'AT'],
             ['AT', 'CATS', 'A', 'CAT'],
             ['CAT', 'CATS', 'AT', 'A'],
             ['CAT', 'A', 'CATS', 'AT'],
             ['CATS', 'A', 'AT', 'CAT'],
             ['CATS', 'CAT', 'AT', 'A'],
             ['CAT', 'CATS', 'A', 'AT'],
             ['AT', 'A', 'CATS', 'CAT']]
            sage: D = Derangements('CART')
            sage: D.list()
            [['A', 'R', 'T', 'C'],
             ['T', 'R', 'C', 'A'],
             ['A', 'T', 'C', 'R'],
             ['R', 'T', 'A', 'C'],
             ['R', 'C', 'T', 'A'],
             ['T', 'C', 'A', 'R'],
             ['T', 'R', 'A', 'C'],
             ['R', 'T', 'C', 'A'],
             ['A', 'C', 'T', 'R']]
            sage: D = Derangements([1,1,2,2,3,3])
            sage: D.list()
            [[2, 2, 3, 3, 1, 1],
             [2, 3, 1, 3, 1, 2],
             [2, 3, 1, 3, 2, 1],
             [2, 3, 3, 1, 1, 2],
             [2, 3, 3, 1, 2, 1],
             [3, 2, 1, 3, 1, 2],
             [3, 2, 1, 3, 2, 1],
             [3, 2, 3, 1, 1, 2],
             [3, 2, 3, 1, 2, 1],
             [3, 3, 1, 1, 2, 2]]
            sage: D = Derangements('SATTAS')
            sage: D.list()
            [['A', 'S', 'S', 'A', 'T', 'T'],
             ['A', 'S', 'A', 'S', 'T', 'T'],
             ['A', 'T', 'S', 'S', 'T', 'A'],
             ['A', 'T', 'S', 'A', 'S', 'T'],
             ['A', 'T', 'A', 'S', 'S', 'T'],
             ['T', 'S', 'S', 'A', 'T', 'A'],
             ['T', 'S', 'A', 'S', 'T', 'A'],
             ['T', 'S', 'A', 'A', 'S', 'T'],
             ['T', 'T', 'S', 'A', 'S', 'A'],
             ['T', 'T', 'A', 'S', 'S', 'A']]
            sage: D = Derangements([1,1,2,2,2])
            sage: D.list()
            []
        """
        if self.__multi:
            for p in Permutations(self._set):
                if not self._fixed_point(p):
                    yield self.element_class(self, list(p))
        else:
            for d in self._iter_der(len(self._set)):
                yield self.element_class(self, [self._set[i - 1] for i in d])

    def _iter_der(self, n):
        r"""
        Iterate through all derangements of the list `[1, 2, 3, \ldots, n]`
        using the method given in [BV2004]_.

        EXAMPLES::

            sage: D = Derangements(4)
            sage: list(D._iter_der(4))
            [[2, 3, 4, 1],
             [4, 3, 1, 2],
             [2, 4, 1, 3],
             [3, 4, 2, 1],
             [3, 1, 4, 2],
             [4, 1, 2, 3],
             [4, 3, 2, 1],
             [3, 4, 1, 2],
             [2, 1, 4, 3]]
        """
        if n <= 1:
            return
        elif n == 2:
            yield [2, 1]
        elif n == 3:
            yield [2, 3, 1]
            yield [3, 1, 2]
        elif n >= 4:
            for d in self._iter_der(n - 1):
                for i in range(1, n):
                    s = d[:]
                    ii = d.index(i)
                    s[ii] = n
                    yield s + [i]
            for d in self._iter_der(n - 2):
                for i in range(1, n):
                    s = d[:]
                    s = [x >= i and x + 1 or x for x in s]
                    s.insert(i - 1, n)
                    yield s + [i]

    def _fixed_point(self, a):
        """
        Return ``True`` if ``a`` has a point in common with ``self._set``.

        EXAMPLES::

            sage: D = Derangements(5)
            sage: D._fixed_point([3,1,2,5,4])
            False
            sage: D._fixed_point([5,4,3,2,1])
            True
        """
        return any(x == y for (x, y) in zip(a, self._set))

    def _count_der(self, n):
        """
        Count the number of derangements of `n` using the recursion
        `D_2 = 1, D_3 = 2, D_n = (n-1) (D_{n-1} + D_{n-2})`.

        EXAMPLES::

            sage: D = Derangements(5)
            sage: D._count_der(2)
            1
            sage: D._count_der(3)
            2
            sage: D._count_der(5)
            44
        """
        if n <= 1:
            return Integer(0)
        if n == 2:
            return Integer(1)
        if n == 3:
            return Integer(2)
        # n >= 4
        last = Integer(2)
        second_last = Integer(1)
        for i in range(4, n + 1):
            current = (i - 1) * (last + second_last)
            second_last = last
            last = current
        return last

    def cardinality(self):
        r"""
        Counts the number of derangements of a positive integer, a
        list, or a string.  The list or string may contain repeated
        elements.  If an integer `n` is given, the value returned
        is the number of derangements of `[1, 2, 3, \ldots, n]`.

        For an integer, or a list or string with all elements
        distinct, the value is obtained by the standard result
        `D_2 = 1, D_3 = 2, D_n = (n-1) (D_{n-1} + D_{n-2})`.

        For a list or string with repeated elements, the number of
        derangements is computed by Macmahon's theorem. If the numbers
        of repeated elements are `a_1, a_2, \ldots, a_k` then the number
        of derangements is given by the coefficient of `x_1 x_2 \cdots
        x_k` in the expansion of `\prod_{i=0}^k (S - s_i)^{a_i}` where
        `S = x_1 + x_2 + \cdots + x_k`.

        EXAMPLES::

            sage: D = Derangements(5)
            sage: D.cardinality()
            44
            sage: D = Derangements([1,44,918,67,254])
            sage: D.cardinality()
            44
            sage: D = Derangements(['A','AT','CAT','CATS','CARTS'])
            sage: D.cardinality()
            44
            sage: D = Derangements('UNCOPYRIGHTABLE')
            sage: D.cardinality()
            481066515734
            sage: D = Derangements([1,1,2,2,3,3])
            sage: D.cardinality()
            10
            sage: D = Derangements('SATTAS')
            sage: D.cardinality()
            10
            sage: D = Derangements([1,1,2,2,2])
            sage: D.cardinality()
            0
        """
        if self.__multi:
            sL = set(self._set)
            A = [self._set.count(i) for i in sL]
            R = PolynomialRing(QQ, 'x', len(A))
            S = sum(i for i in R.gens())
            e = prod((S - x)**y for (x, y) in zip(R.gens(), A))
            return Integer(e.coefficient(dict([(x, y) for (x, y) in zip(R.gens(), A)])))
        return self._count_der(len(self._set))

    def _rand_der(self):
        r"""
        Produces a random derangement of `[1, 2, \ldots, n]`.

        This is an
        implementation of the algorithm described by Martinez et. al. in
        [MPP2008]_.

        EXAMPLES::

            sage: D = Derangements(4)
            sage: d = D._rand_der()
            sage: d in D
            True
        """
        n = len(self._set)
        A = list(range(1, n + 1))
        mark = [x < 0 for x in A]
        i, u = n, n
        while u >= 2:
            if not(mark[i - 1]):
                while True:
                    j = randrange(1, i)
                    if not(mark[j - 1]):
                        A[i - 1], A[j - 1] = A[j - 1], A[i - 1]
                        break
                p = random()
                if p < (u - 1) * self._count_der(u - 2) // self._count_der(u):
                    mark[j - 1] = True
                    u -= 1
                u -= 1
            i -= 1
        return A

    def random_element(self):
        r"""
        Produces all derangements of a positive integer, a list, or
        a string.  The list or string may contain repeated elements.
        If an integer `n` is given, then a random
        derangements of `[1, 2, 3, \ldots, n]` is returned

        For an integer, or a list or string with all elements
        distinct, the value is obtained by an algorithm described in
        [MPP2008]_. For a list or string with repeated elements the
        derangement is formed by choosing an element at random from the list of
        all possible derangements.

        OUTPUT:

        A single list or string containing a derangement, or an
        empty list if there are no derangements.

        EXAMPLES::

            sage: D = Derangements(4)
            sage: D.random_element() # random
            [2, 3, 4, 1]
            sage: D = Derangements(['A','AT','CAT','CATS','CARTS','CARETS'])
            sage: D.random_element() # random
            ['AT', 'CARTS', 'A', 'CAT', 'CARETS', 'CATS']
            sage: D = Derangements('UNCOPYRIGHTABLE')
            sage: D.random_element() # random
            ['C', 'U', 'I', 'H', 'O', 'G', 'N', 'B', 'E', 'L', 'A', 'R', 'P', 'Y', 'T']
            sage: D = Derangements([1,1,1,1,2,2,2,2,3,3,3,3])
            sage: D.random_element() # random
            [3, 2, 2, 3, 1, 3, 1, 3, 2, 1, 1, 2]
            sage: D = Derangements('ESSENCES')
            sage: D.random_element() # random
            ['N', 'E', 'E', 'C', 'S', 'S', 'S', 'E']
            sage: D = Derangements([1,1,2,2,2])
            sage: D.random_element()
            []

        TESTS:

        Check that index error discovered in :trac:`29974` is fixed::

            sage: D = Derangements([1,1,2,2])
            sage: _ = [D.random_element() for _ in range(20)]
        """
        if self.__multi:
            L = list(self)
            if len(L) == 0:
                return self.element_class(self, [])
            i = randrange(len(L))
            return L[i]
        temp = self._rand_der()
        return self.element_class(self, [self._set[ii - 1] for ii in temp])
