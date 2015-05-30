"""
Deprecated combinations

AUTHORS:

- Mike Hansen (2007): initial implementation

- Vincent Delecroix (2014): deprecation
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
from sage.rings.arith import binomial

def ChooseNK(n, k):
    """
    All possible choices of k elements out of range(n) without repetitions.

    The elements of the output are tuples of Python int (and not Sage Integer).

    This was deprecated in :trac:`10534` for :func:`Combinations`
    (or ``itertools.combinations`` for doing iteration).

    EXAMPLES::

        sage: from sage.combinat.choose_nk import ChooseNK
        sage: c = ChooseNK(4,2)
        doctest:...: DeprecationWarning: ChooseNk is deprecated and will be
        removed. Use Combinations instead (or combinations from the itertools
        module for iteration)
        See http://trac.sagemath.org/10534 for details.
        sage: c.first()
        [0, 1]
        sage: c.list()
        [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]
    """
    from sage.misc.superseded import deprecation
    deprecation(10534, "ChooseNk is deprecated and will be removed. Use Combinations instead (or combinations from the itertools module for iteration)")
    from sage.combinat.combination import Combinations
    return Combinations(n,k)

#TODO: the following functions are used sage.combinat.combination and
#      sage.combinat.subset. It might be good to move them somewhere else.
def rank(comb, n, check=True):
    """
    Return the rank of ``comb`` in the subsets of ``range(n)`` of size ``k``
    where ``k`` is the length of ``comb``.

    The algorithm used is based on combinadics and James McCaffrey's
    MSDN article. See: :wikipedia:`Combinadic`.

    EXAMPLES::

        sage: import sage.combinat.choose_nk as choose_nk
        sage: choose_nk.rank((), 3)
        0
        sage: choose_nk.rank((0,), 3)
        0
        sage: choose_nk.rank((1,), 3)
        1
        sage: choose_nk.rank((2,), 3)
        2
        sage: choose_nk.rank((0,1), 3)
        0
        sage: choose_nk.rank((0,2), 3)
        1
        sage: choose_nk.rank((1,2), 3)
        2
        sage: choose_nk.rank((0,1,2), 3)
        0

        sage: choose_nk.rank((0,1,2,3), 3)
        Traceback (most recent call last):
        ...
        ValueError: len(comb) must be <= n
        sage: choose_nk.rank((0,0), 2)
        Traceback (most recent call last):
        ...
        ValueError: comb must be a subword of (0,1,...,n)

        sage: choose_nk.rank([1,2], 3)
        2
        sage: choose_nk.rank([0,1,2], 3)
        0
    """
    k = len(comb)
    if check:
        if k > n:
            raise ValueError("len(comb) must be <= n")
        comb = map(int, comb)
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
    """
    Returns the largest w < a such that binomial(w,b) <= x.

    EXAMPLES::

        sage: from sage.combinat.choose_nk import _comb_largest
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
    """
    Returns the combination of rank r in the subsets of range(n) of
    size k when listed in lexicographic order.

    The algorithm used is based on combinadics and James McCaffrey's
    MSDN article. See: http://en.wikipedia.org/wiki/Combinadic

    EXAMPLES::

        sage: import sage.combinat.choose_nk as choose_nk
        sage: choose_nk.from_rank(0,3,0)
        ()
        sage: choose_nk.from_rank(0,3,1)
        (0,)
        sage: choose_nk.from_rank(1,3,1)
        (1,)
        sage: choose_nk.from_rank(2,3,1)
        (2,)
        sage: choose_nk.from_rank(0,3,2)
        (0, 1)
        sage: choose_nk.from_rank(1,3,2)
        (0, 2)
        sage: choose_nk.from_rank(2,3,2)
        (1, 2)
        sage: choose_nk.from_rank(0,3,3)
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
