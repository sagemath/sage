# sage_setup: distribution = sagemath-categories
r"""
Combinations

AUTHORS:

- Mike Hansen (2007): initial implementation

- Vincent Delecroix (2011): cleaning, bug corrections, doctests

- Antoine Genitrini (2020) : new implementation of the lexicographic unranking of combinations
"""
# ****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
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
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import itertools

from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.arith.misc import binomial
from .integer_vector import IntegerVectors
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.structure.parent import Parent
from sage.misc.persist import register_unpickle_override


def Combinations(mset, k=None):
    """
    Return the combinatorial class of combinations of the multiset
    ``mset``. If ``k`` is specified, then it returns the combinatorial
    class of combinations of ``mset`` of size ``k``.

    A *combination* of a multiset `M` is an unordered selection of `k`
    objects of `M`, where every object can appear at most as many
    times as it appears in `M`.

    The combinatorial classes correctly handle the cases where ``mset`` has
    duplicate elements.

    EXAMPLES::

        sage: C = Combinations(range(4)); C
        Combinations of [0, 1, 2, 3]
        sage: C.list()
        [[],
         [0],
         [1],
         [2],
         [3],
         [0, 1],
         [0, 2],
         [0, 3],
         [1, 2],
         [1, 3],
         [2, 3],
         [0, 1, 2],
         [0, 1, 3],
         [0, 2, 3],
         [1, 2, 3],
         [0, 1, 2, 3]]
         sage: C.cardinality()
         16

    ::

        sage: C2 = Combinations(range(4),2); C2
        Combinations of [0, 1, 2, 3] of length 2
        sage: C2.list()
        [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]
        sage: C2.cardinality()
        6

    ::

        sage: Combinations([1,2,2,3]).list()
        [[],
         [1],
         [2],
         [3],
         [1, 2],
         [1, 3],
         [2, 2],
         [2, 3],
         [1, 2, 2],
         [1, 2, 3],
         [2, 2, 3],
         [1, 2, 2, 3]]

    ::

        sage: Combinations([1,2,3], 2).list()
        [[1, 2], [1, 3], [2, 3]]

    ::

        sage: mset = [1,1,2,3,4,4,5]
        sage: Combinations(mset,2).list()
        [[1, 1],
         [1, 2],
         [1, 3],
         [1, 4],
         [1, 5],
         [2, 3],
         [2, 4],
         [2, 5],
         [3, 4],
         [3, 5],
         [4, 4],
         [4, 5]]

    ::

        sage: mset = ["d","a","v","i","d"]
        sage: Combinations(mset,3).list()
        [['d', 'd', 'a'],
         ['d', 'd', 'v'],
         ['d', 'd', 'i'],
         ['d', 'a', 'v'],
         ['d', 'a', 'i'],
         ['d', 'v', 'i'],
         ['a', 'v', 'i']]

    ::

        sage: X = Combinations([1,2,3,4,5],3)
        sage: [x for x in X]
        [[1, 2, 3],
         [1, 2, 4],
         [1, 2, 5],
         [1, 3, 4],
         [1, 3, 5],
         [1, 4, 5],
         [2, 3, 4],
         [2, 3, 5],
         [2, 4, 5],
         [3, 4, 5]]

    It is possible to take combinations of Sage objects::

        sage: Combinations([vector([1,1]), vector([2,2]), vector([3,3])], 2).list()     # needs sage.modules
        [[(1, 1), (2, 2)], [(1, 1), (3, 3)], [(2, 2), (3, 3)]]

    TESTS:

    Run the test suites::

        sage: C = Combinations([2,3])
        sage: TestSuite(C).run()
        sage: C = Combinations([2,3], 1)
        sage: TestSuite(C).run()

    We check that the code works even for non mutable objects::

        sage: l = [vector((0,0)), vector((0,1))]                                        # needs sage.modules
        sage: Combinations(l).list()                                                    # needs sage.modules
        [[], [(0, 0)], [(0, 1)], [(0, 0), (0, 1)]]
    """
    # Check to see if everything in mset is unique
    is_unique = False
    if isinstance(mset, (int, Integer)):
        mset = list(range(mset))
        is_unique = True
    elif isinstance(mset, range):
        mset = list(mset)
        is_unique = True
    else:
        mset = list(mset)
        for i, e in enumerate(mset):
            if mset.index(e) != i:
                break
        else:
            is_unique = True

    if is_unique:
        if k is None:
            return Combinations_set(mset)
        else:
            return Combinations_setk(mset, k)
    else:
        if k is None:
            return Combinations_mset(mset)
        else:
            return Combinations_msetk(mset, k)


class Combinations_mset(Parent):
    def __init__(self, mset):
        """
        TESTS::

            sage: C = Combinations(range(4))
            sage: C == loads(dumps(C))
            True
        """
        self.mset = mset
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def __contains__(self, x) -> bool:
        """
        EXAMPLES::

            sage: c = Combinations(range(4))
            sage: all( i in c for i in c )
            True
            sage: [3,4] in c
            False
            sage: [0,0] in c
            False
        """
        try:
            x = list(x)
        except TypeError:
            return False

        return all(i in self.mset for i in x) and len(set(x)) == len(x)

    def __eq__(self, other) -> bool:
        """
        Test for equality.

        EXAMPLES::

            sage: c = Combinations([1,2,2,3])
            sage: c == Combinations((1,2,2,3))
            True
            sage: c == Combinations([3,4,4,6])
            False
        """
        return isinstance(other, Combinations_mset) and self.mset == other.mset

    def __ne__(self, other) -> bool:
        """
        Test for unequality.

        EXAMPLES::

            sage: c = Combinations([1,2,2])
            sage: c != Combinations([1,2,3,3])
            True
        """
        return not (self == other)

    def __repr__(self) -> str:
        """
        TESTS::

            sage: repr(Combinations(range(4)))
            'Combinations of [0, 1, 2, 3]'
        """
        return "Combinations of {}".format(self.mset)

    def __iter__(self):
        """
        TESTS::

            sage: Combinations(['a','a','b']).list() #indirect doctest
            [[], ['a'], ['b'], ['a', 'a'], ['a', 'b'], ['a', 'a', 'b']]
        """
        for k in range(len(self.mset) + 1):
            yield from Combinations_msetk(self.mset, k)

    def cardinality(self) -> Integer:
        """
        TESTS::

            sage: Combinations([1,2,3]).cardinality()
            8
            sage: Combinations(['a','a','b']).cardinality()                             # needs sage.libs.gap
            6
        """
        return ZZ.sum(Combinations_msetk(self.mset, k).cardinality()
                      for k in range(len(self.mset) + 1))


class Combinations_set(Combinations_mset):
    def __iter__(self):
        """
        EXAMPLES::

            sage: Combinations([1,2,3]).list() #indirect doctest
            [[], [1], [2], [3], [1, 2], [1, 3], [2, 3], [1, 2, 3]]
        """
        for k in range(len(self.mset) + 1):
            yield from Combinations_setk(self.mset, k)

    def unrank(self, r):
        """
        EXAMPLES::

            sage: c = Combinations([1,2,3])
            sage: c.list() == list(map(c.unrank, range(c.cardinality())))
            True
        """
        k = 0
        n = len(self.mset)
        b = binomial(n, k)
        while r >= b:
            r -= b
            k += 1
            b = binomial(n, k)

        return [self.mset[i] for i in from_rank(r, n, k)]

    def rank(self, x):
        """
        EXAMPLES::

            sage: c = Combinations([1,2,3])
            sage: list(range(c.cardinality())) == list(map(c.rank, c))
            True
        """
        x = [self.mset.index(i) for i in x]
        r = 0
        n = len(self.mset)
        for i in range(len(x)):
            r += binomial(n, i)
        r += rank(x, n)
        return r

    def cardinality(self):
        """
        Return the size of Combinations(set).

        EXAMPLES::

            sage: Combinations(range(16000)).cardinality() == 2^16000
            True
        """
        return ZZ(2)**len(self.mset)


class Combinations_msetk(Parent):
    def __init__(self, mset, k):
        """
        TESTS::

            sage: C = Combinations([1,2,3],2)
            sage: C == loads(dumps(C))
            True
        """
        self.mset = mset
        self.k = k
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def __contains__(self, x) -> bool:
        """
        EXAMPLES::

            sage: c = Combinations(range(4),2)
            sage: all( i in c for i in c )
            True
            sage: [0,1] in c
            True
            sage: [0,1,2] in c
            False
            sage: [3,4] in c
            False
            sage: [0,0] in c
            False
        """
        try:
            x = list(x)
        except TypeError:
            return False
        return x in Combinations_mset(self.mset) and len(x) == self.k

    def __eq__(self, other) -> bool:
        """
        Test for equality.

        EXAMPLES::

            sage: c = Combinations([1,2,2,3],3)
            sage: c == Combinations((1,2,2,3), 3)
            True
            sage: c == Combinations([1,2,2,3], 2)
            False
        """
        return (isinstance(other, Combinations_msetk) and
                self.mset == other.mset and self.k == other.k)

    def __ne__(self, other) -> bool:
        """
        Test for unequality.

        EXAMPLES::

            sage: c = Combinations([1,2,2,3],3)
            sage: c != Combinations((1,2,2,3), 2)
            True
        """
        return not (self == other)

    def __repr__(self) -> str:
        """
        TESTS::

            sage: repr(Combinations([1,2,2,3],2))
            'Combinations of [1, 2, 2, 3] of length 2'
        """
        return "Combinations of {} of length {}".format(self.mset, self.k)

    def __iter__(self):
        """
        EXAMPLES::

            sage: Combinations(['a','a','b'],2).list() # indirect doctest
            [['a', 'a'], ['a', 'b']]
        """
        items = [self.mset.index(x) for x in self.mset]
        indices = sorted(set(items))
        counts = [0] * len(indices)
        for i in items:
            counts[indices.index(i)] += 1
        for iv in IntegerVectors(self.k, len(indices), outer=counts):
            yield sum([[self.mset[indices[i]]] * iv[i]
                       for i in range(len(indices))], [])

    def cardinality(self) -> Integer:
        """
        Return the size of combinations(mset, k).

        IMPLEMENTATION: Wraps GAP's NrCombinations.

        EXAMPLES::

            sage: mset = [1,1,2,3,4,4,5]
            sage: Combinations(mset,2).cardinality()                                    # needs sage.libs.gap
            12
        """
        from sage.libs.gap.libgap import libgap
        items = [self.mset.index(i) for i in self.mset]
        nc = libgap.function_factory('NrCombinations')
        return ZZ(nc(items, ZZ(self.k)))


class Combinations_setk(Combinations_msetk):
    def _iterator(self, items, n):
        """
        An iterator for all the n-combinations of items.

        EXAMPLES::

            sage: it = Combinations([1,2,3,4],3)._iterator([1,2,3,4],3)
            sage: list(it)
            [[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]]
        """
        for combination in itertools.combinations(items, n):
            yield list(combination)

    def _iterator_zero(self):
        """
        An iterator which just returns the empty list.

        EXAMPLES::

            sage: it = Combinations([1,2,3,4,5],3)._iterator_zero()
            sage: list(it)
            [[]]
        """
        yield []

    def __iter__(self):
        r"""
        Uses Python's :func:`itertools.combinations` to iterate through all
        of the combinations.

        EXAMPLES::

            sage: Combinations([1,2,3,4,5],3).list() # indirect doctest
            [[1, 2, 3],
             [1, 2, 4],
             [1, 2, 5],
             [1, 3, 4],
             [1, 3, 5],
             [1, 4, 5],
             [2, 3, 4],
             [2, 3, 5],
             [2, 4, 5],
             [3, 4, 5]]
        """
        if self.k == 0:
            return self._iterator_zero()
        else:
            return self._iterator(self.mset, self.k)

    def list(self) -> list:
        """
        EXAMPLES::

            sage: Combinations([1,2,3,4,5],3).list()
            [[1, 2, 3],
             [1, 2, 4],
             [1, 2, 5],
             [1, 3, 4],
             [1, 3, 5],
             [1, 4, 5],
             [2, 3, 4],
             [2, 3, 5],
             [2, 4, 5],
             [3, 4, 5]]
        """
        return list(self)

    def unrank(self, r):
        """
        EXAMPLES::

            sage: c = Combinations([1,2,3], 2)
            sage: c.list() == list(map(c.unrank, range(c.cardinality())))
            True
        """
        return [self.mset[i] for i in from_rank(r, len(self.mset), self.k)]

    def rank(self, x):
        """
        EXAMPLES::

            sage: c = Combinations([1,2,3], 2)
            sage: list(range(c.cardinality())) == list(map(c.rank, c.list()))
            True
        """
        x = [self.mset.index(i) for i in x]
        return rank(x, len(self.mset))

    def cardinality(self) -> Integer:
        """
        Return the size of combinations(set, k).

        EXAMPLES::

            sage: Combinations(range(16000), 5).cardinality()
            8732673194560003200
        """
        return ZZ(binomial(len(self.mset), self.k))


def rank(comb, n, check=True):
    """
    Return the rank of ``comb`` in the subsets of ``range(n)`` of size ``k``
    where ``k`` is the length of ``comb``.

    The algorithm used is based on combinadics and James McCaffrey's
    MSDN article. See: :wikipedia:`Combinadic`.

    EXAMPLES::

        sage: import sage.combinat.combination as combination
        sage: combination.rank((), 3)
        0
        sage: combination.rank((0,), 3)
        0
        sage: combination.rank((1,), 3)
        1
        sage: combination.rank((2,), 3)
        2
        sage: combination.rank((0,1), 3)
        0
        sage: combination.rank((0,2), 3)
        1
        sage: combination.rank((1,2), 3)
        2
        sage: combination.rank((0,1,2), 3)
        0

        sage: combination.rank((0,1,2,3), 3)
        Traceback (most recent call last):
        ...
        ValueError: len(comb) must be <= n
        sage: combination.rank((0,0), 2)
        Traceback (most recent call last):
        ...
        ValueError: comb must be a subword of (0,1,...,n)

        sage: combination.rank([1,2], 3)
        2
        sage: combination.rank([0,1,2], 3)
        0
    """
    k = len(comb)
    if check:
        if k > n:
            raise ValueError("len(comb) must be <= n")
        comb = [int(i) for i in comb]
        for i in range(k - 1):
            if comb[i + 1] <= comb[i]:
                raise ValueError("comb must be a subword of (0,1,...,n)")

    # Generate the combinadic from the combination

    # w = [n-1-comb[i] for i in range(k)]

    # Calculate the integer that is the dual of
    # the lexicographic index of the combination
    r = k
    t = 0
    for i in range(k):
        t += binomial(n - 1 - comb[i], r)
        r -= 1

    return binomial(n, k) - t - 1


def from_rank(r, n, k):
    r"""
    Return the combination of rank ``r`` in the subsets of
    ``range(n)`` of size ``k`` when listed in lexicographic order.

    The algorithm used is based on factoradics and presented in [DGH2020]_.
    It is there compared to the other from the literature.

    EXAMPLES::

        sage: import sage.combinat.combination as combination
        sage: combination.from_rank(0,3,0)
        ()
        sage: combination.from_rank(0,3,1)
        (0,)
        sage: combination.from_rank(1,3,1)
        (1,)
        sage: combination.from_rank(2,3,1)
        (2,)
        sage: combination.from_rank(0,3,2)
        (0, 1)
        sage: combination.from_rank(1,3,2)
        (0, 2)
        sage: combination.from_rank(2,3,2)
        (1, 2)
        sage: combination.from_rank(0,3,3)
        (0, 1, 2)

    TESTS::

        sage: from sage.combinat.combination import from_rank
        sage: def _comb_largest(a,b,x):
        ....:     w = a - 1
        ....:     while binomial(w,b) > x:
        ....:         w -= 1
        ....:     return w
        sage: def from_rank_comb_largest(r, n, k):
        ....:     a = n
        ....:     b = k
        ....:     x = binomial(n, k) - 1 - r  # x is the 'dual' of m
        ....:     comb = [None] * k
        ....:     for i in range(k):
        ....:         comb[i] = _comb_largest(a, b, x)
        ....:         x = x - binomial(comb[i], b)
        ....:         a = comb[i]
        ....:         b = b - 1
        ....:     for i in range(k):
        ....:         comb[i] = (n - 1) - comb[i]
        ....:     return tuple(comb)
        sage: all(from_rank(r, n, k) == from_rank_comb_largest(r, n, k)                 # needs sage.symbolic
        ....:     for n in range(10) for k in range(n+1) for r in range(binomial(n,k)))
        True
    """
    if k < 0:
        raise ValueError("k must be > 0")
    if k > n:
        raise ValueError("k must be <= n")
    if n == 0 or k == 0:
        return ()
    if n < 0:
        raise ValueError("n must be >= 0")
    B = binomial(n, k)
    if r < 0 or r >= B:
        raise ValueError("r must satisfy  0 <= r < binomial(n, k)")
    if k == 1:
        return (r,)

    n0 = n
    D = [0] * k
    inverse = False
    if k < n0 / 2:
        inverse = True
        k = n - k
        r = B - 1 - r

    B = (B * k) // n0
    m = 0
    i = 0
    j = 0
    m2 = 0
    d = 0
    while d < k - 1:
        if B > r:
            if i < k - 2:
                if n0 - 1 - m == 0:
                    B = 1
                else:
                    B = (B * (k - 1 - i)) // (n0 - 1 - m)
            d += 1
            if inverse:
                for e in range(m2, m + i):
                    D[j] = e
                    j += 1
                m2 = m + i + 1
            else:
                D[i] = m + i
            i += 1
            n0 -= 1
        else:
            r -= B
            if n0 - 1 - m == 0:
                B = 1
            else:
                B = (B * (n0 - m - k + i)) // (n0 - 1 - m)
            m += 1
    if inverse:
        for e in range(m2, n0 + r + i - B):
            D[j] = e
            j += 1
        for e in range(n0 + r + i + 1 - B, n):
            D[j] = e
            j += 1
    else:
        D[k - 1] = n0 + r + k - 1 - B
    return tuple(D)

##########################################################
# Deprecations


class ChooseNK(Combinations_setk):
    def __setstate__(self, state):
        r"""
        For unpickling old ``ChooseNK`` objects.

        TESTS::

            sage: loads(b"x\x9ck`J.NLO\xd5K\xce\xcfM\xca\xccK,\xd1K\xce\xc8\xcf"
            ....:   b"/N\x8d\xcf\xcb\xe6r\x06\xb3\xfc\xbc\xb9\n\x195\x1b\x0b"
            ....:   b"\x99j\x0b\x995B\x99\xe2\xf3\nY :\x8a2\xf3\xd2\x8b\xf52"
            ....:   b"\xf3JR\xd3S\x8b\xb8r\x13\xb3S\xe3a\x9cB\xd6PF\xd3\xd6\xa0"
            ....:   b"B6\xa0\xfa\xecB\xf6\x0c \xd7\x08\xc8\xe5(M\xd2\x03\x00{"
            ....:   b"\x82$\xd8")
            Combinations of [0, 1, 2, 3, 4] of length 2
        """
        self.__class__ = Combinations_setk
        Combinations_setk.__init__(self, list(range(state['_n'])), state['_k'])


register_unpickle_override("sage.combinat.choose_nk", "ChooseNK", ChooseNK)
