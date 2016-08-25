r"""
Combinations

AUTHORS:

- Mike Hansen (2007): initial implementation

- Vincent Delecroix (2011): cleaning, bug corrections, doctests

"""
#*****************************************************************************
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
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import

from sage.interfaces.all import gap
from sage.rings.all import ZZ, Integer
from sage.arith.all import binomial
from .combinat import CombinatorialClass
from .integer_vector import IntegerVectors
from sage.misc.misc import uniq


def Combinations(mset, k=None):
    """
    Return the combinatorial class of combinations of the multiset
    ``mset``. If ``k`` is specified, then it returns the combinatorial
    class of combinations of ``mset`` of size ``k``.

    A *combination* of a multiset `M` is an unordered selection of `k`
    objects of `M`, where every object can appear at most as many
    times as it appears in `M`.

    The combinatorial classes correctly handle the cases where mset has
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

        sage: Combinations([vector([1,1]), vector([2,2]), vector([3,3])], 2).list()
        [[(1, 1), (2, 2)], [(1, 1), (3, 3)], [(2, 2), (3, 3)]]

    TESTS:

    We check that the code works even for non mutable objects::

        sage: l = [vector((0,0)), vector((0,1))]
        sage: Combinations(l).list()
        [[], [(0, 0)], [(0, 1)], [(0, 0), (0, 1)]]
    """



    #Check to see if everything in mset is unique
    if isinstance(mset, (int, Integer)):
        mset = range(mset)
    else:
        mset = list(mset)

    d = {}
    for i in mset:
        d[mset.index(i)] = 1

    if len(d) == len(mset):
        if k is None:
            return Combinations_set(mset)
        else:
            return Combinations_setk(mset,k)
    else:
        if k is None:
            return Combinations_mset(mset)
        else:
            return Combinations_msetk(mset,k)

class Combinations_mset(CombinatorialClass):
    def __init__(self, mset):
        """
        TESTS::

            sage: C = Combinations(range(4))
            sage: C == loads(dumps(C))
            True
        """
        self.mset = mset

    def __contains__(self, x):
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

        return all(i in self.mset for i in x) and len(uniq(x)) == len(x)


    def __repr__(self):
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
        for k in range(len(self.mset)+1):
            for comb in Combinations_msetk(self.mset, k):
                yield comb

    def cardinality(self):
        """
        TESTS::

            sage: Combinations([1,2,3]).cardinality()
            8
            sage: Combinations(['a','a','b']).cardinality()
            6
        """
        c = 0
        for k in range(len(self.mset) + 1):
            c += Combinations_msetk(self.mset, k).cardinality()
        return c

class Combinations_set(Combinations_mset):
    def __iter__(self):
        """
        EXAMPLES::

            sage: Combinations([1,2,3]).list() #indirect doctest
            [[], [1], [2], [3], [1, 2], [1, 3], [2, 3], [1, 2, 3]]
        """
        for k in range(len(self.mset) + 1):
            for comb in Combinations_setk(self.mset, k):
                yield comb


    def unrank(self, r):
        """
        EXAMPLES::

            sage: c = Combinations([1,2,3])
            sage: c.list() == map(c.unrank, range(c.cardinality()))
            True
        """
        k = 0
        n = len(self.mset)
        b = binomial(n, k)
        while r >= b:
            r -= b
            k += 1
            b = binomial(n,k)

        return [self.mset[i] for i in from_rank(r, n, k)]


    def rank(self, x):
        """
        EXAMPLES::

            sage: c = Combinations([1,2,3])
            sage: range(c.cardinality()) == map(c.rank, c)
            True
        """
        x = [self.mset.index(_) for _ in x]
        r = 0
        n = len(self.mset)
        for i in range(len(x)):
            r += binomial(n, i)
        r += rank(x, n)
        return r

class Combinations_msetk(CombinatorialClass):
    def __init__(self, mset, k):
        """
        TESTS::

            sage: C = Combinations([1,2,3],2)
            sage: C == loads(dumps(C))
            True
        """
        self.mset = mset
        self.k = k

    def __contains__(self, x):
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


    def __repr__(self):
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
        items = map(self.mset.index, self.mset)
        indices = uniq(sorted(items))
        counts = [0] * len(indices)
        for i in items:
            counts[indices.index(i)] += 1
        for iv in IntegerVectors(self.k, len(indices), outer=counts):
            yield sum([[self.mset[indices[i]]]*iv[i] for i in range(len(indices))],[])

    def cardinality(self):
        """
        Returns the size of combinations(mset,k). IMPLEMENTATION: Wraps
        GAP's NrCombinations.

        EXAMPLES::

            sage: mset = [1,1,2,3,4,4,5]
            sage: Combinations(mset,2).cardinality()
            12
        """
        items = [self.mset.index(_) for _ in self.mset]
        return ZZ(gap.eval("NrCombinations({},{})".format(items, ZZ(self.k))))



class Combinations_setk(Combinations_msetk):
    def _iterator(self, items, len_items,  n):
        """
        An iterator for all the n-combinations of items.

        EXAMPLES::

            sage: it = Combinations([1,2,3,4],3)._iterator([1,2,3,4],4,3)
            sage: list(it)
            [[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]]
        """
        for i in range(len_items):
            v = items[i:i+1]
            if n == 1:
                yield v
            else:
                rest = items[i+1:]
                for c in self._iterator(rest, len_items-i-1,  n-1):
                    yield v + c

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
        Posted by Raymond Hettinger, 2006/03/23, to the Python Cookbook:
        http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/474124

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
            return self._iterator(self.mset, len(self.mset), self.k)


    def list(self):
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
            sage: c.list() == map(c.unrank, range(c.cardinality()))
            True
        """
        return [self.mset[i] for i in from_rank(r, len(self.mset), self.k)]


    def rank(self, x):
        """
        EXAMPLES::

            sage: c = Combinations([1,2,3], 2)
            sage: range(c.cardinality()) == map(c.rank, c.list())
            True
        """
        x = [self.mset.index(_) for _ in x]
        return rank(x, len(self.mset))


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
        comb = [int(_) for _ in comb]
        for i in xrange(k - 1):
            if comb[i + 1] <= comb[i]:
                raise ValueError("comb must be a subword of (0,1,...,n)")

    #Generate the combinadic from the
    #combination

    #w = [n-1-comb[i] for i in xrange(k)]

    #Calculate the integer that is the dual of
    #the lexicographic index of the combination
    r = k
    t = 0
    for i in range(k):
        t += binomial(n - 1 - comb[i], r)
        r -= 1

    return binomial(n,k)-t-1

def _comb_largest(a,b,x):
    r"""
    Returns the largest `w < a` such that `binomial(w,b) <= x`.

    EXAMPLES::

        sage: from sage.combinat.combination import _comb_largest
        sage: _comb_largest(6,3,10)
        5
        sage: _comb_largest(6,3,5)
        4
    """
    w = a - 1

    while binomial(w,b) > x:
        w -= 1

    return w

def from_rank(r, n, k):
    r"""
    Returns the combination of rank ``r`` in the subsets of
    ``range(n)`` of size ``k`` when listed in lexicographic order.

    The algorithm used is based on combinadics and James McCaffrey's
    MSDN article. See: :wikipedia:`Combinadic`

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
    """
    if k < 0:
        raise ValueError("k must be > 0")
    if k > n:
        raise ValueError("k must be <= n")

    a = n
    b = k
    x = binomial(n, k) - 1 - r  # x is the 'dual' of m
    comb = [None] * k

    for i in xrange(k):
        comb[i] = _comb_largest(a, b, x)
        x = x - binomial(comb[i], b)
        a = comb[i]
        b = b - 1

    for i in xrange(k):
        comb[i] = (n - 1) - comb[i]

    return tuple(comb)

##########################################################
# Deprecations

class ChooseNK(Combinations_setk):
    def __setstate__(self, state):
        r"""
        For unpickling old ``ChooseNK`` objects.

        TESTS::

            sage: loads("x\x9ck`J.NLO\xd5K\xce\xcfM\xca\xccK,\xd1K\xce\xc8\xcf"
            ....:   "/N\x8d\xcf\xcb\xe6r\x06\xb3\xfc\xbc\xb9\n\x195\x1b\x0b"
            ....:   "\x99j\x0b\x995B\x99\xe2\xf3\nY :\x8a2\xf3\xd2\x8b\xf52"
            ....:   "\xf3JR\xd3S\x8b\xb8r\x13\xb3S\xe3a\x9cB\xd6PF\xd3\xd6\xa0"
            ....:   "B6\xa0\xfa\xecB\xf6\x0c \xd7\x08\xc8\xe5(M\xd2\x03\x00{"
            ....:   "\x82$\xd8")
            Combinations of [0, 1, 2, 3, 4] of length 2
        """
        self.__class__ = Combinations_setk
        Combinations_setk.__init__(self, range(state['_n']), state['_k'])

from sage.structure.sage_object import register_unpickle_override
register_unpickle_override("sage.combinat.choose_nk", "ChooseNK", ChooseNK)
