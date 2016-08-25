# -*- coding: utf-8 -*-
"""
Lyndon words
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent

from sage.combinat.composition import Composition, Compositions
from sage.rings.all import Integer
from sage.arith.all import factorial, divisors, gcd, moebius
from sage.misc.all import prod
from six.moves import builtins
from . import necklace
from .integer_vector import IntegerVectors
from sage.combinat.words.words import FiniteWords


def LyndonWords(e=None, k=None):
    """
    Returns the combinatorial class of Lyndon words.

    A Lyndon word `w` is a word that is lexicographically less than all of
    its rotations.  Equivalently, whenever `w` is split into two non-empty
    substrings, `w` is lexicographically less than the right substring.

    INPUT:

    - no input at all

    or

    - ``e`` - integer, size of alphabet
    - ``k`` - integer, length of the words

    or

    - ``e`` - a composition

    OUTPUT:

    A combinatorial class of Lyndon words.

    EXAMPLES::

        sage: LyndonWords()
        Lyndon words

    If e is an integer, then e specifies the length of the
    alphabet; k must also be specified in this case::

        sage: LW = LyndonWords(3,3); LW
        Lyndon words from an alphabet of size 3 of length 3
        sage: LW.first()
        word: 112
        sage: LW.last()
        word: 233
        sage: LW.random_element() # random
        word: 112
        sage: LW.cardinality()
        8

    If e is a (weak) composition, then it returns the class of Lyndon
    words that have evaluation e::

        sage: LyndonWords([2, 0, 1]).list()
        [word: 113]
        sage: LyndonWords([2, 0, 1, 0, 1]).list()
        [word: 1135, word: 1153, word: 1315]
        sage: LyndonWords([2, 1, 1]).list()
        [word: 1123, word: 1132, word: 1213]
    """
    if e is None and k is None:
        return LyndonWords_class()
    elif isinstance(e, (int, Integer)):
        if e > 0:
            if not isinstance(k, (int, Integer)):
                raise TypeError("k must be a non-negative integer")
            if k < 0:
                raise TypeError("k must be a non-negative integer")
            return LyndonWords_nk(Integer(e), Integer(k))
    elif e in Compositions():
        return LyndonWords_evaluation(Composition(e))

    raise TypeError("e must be a positive integer or a composition")

def LyndonWord(data, check=True):
    r"""
    Construction of a Lyndon word.

    INPUT:

    - ``data`` - list
    - ``check`` - bool (optional, default: True) if True, a
      verification that the input data represent a Lyndon word.

    OUTPUT:

    A Lyndon word.

    EXAMPLES::

        sage: LyndonWord([1,2,2])
        word: 122
        sage: LyndonWord([1,2,3])
        word: 123
        sage: LyndonWord([2,1,2,3])
        Traceback (most recent call last):
        ...
        ValueError: not a Lyndon word

    If check is False, then no verification is done::

        sage: LyndonWord([2,1,2,3], check=False)
        word: 2123
    """
    return LyndonWords()(data, check=check)

class LyndonWords_class(UniqueRepresentation, Parent):
    r"""
    The set of all Lyndon words.
    """
    def __init__(self, alphabet=None):
        r"""
        INPUT:

        - ``alphabet`` -- the underlying alphabet

        TESTS::

            sage: loads(dumps(LyndonWords())) is LyndonWords()
            True
        """
        from sage.categories.sets_cat import Sets
        self._words = FiniteWords()
        Parent.__init__(self, category=Sets().Infinite(), facade=(self._words))

    def __call__(self, *args, **kwds):
        r"""
        TESTS::

            sage: L = LyndonWords()
            sage: L('aababc')
            word: aababc
            sage: L([2,0,1])
            Traceback (most recent call last):
            ...
            ValueError: not a Lyndon word
        """
        w = self._words(*args, **kwds)
        if kwds.get('check', True) and not w.is_lyndon():
            raise ValueError("not a Lyndon word")
        return w

    def __repr__(self):
        r"""
        String representation.

        EXAMPLES::

            sage: LyndonWords()
            Lyndon words
        """
        return "Lyndon words"

    def __contains__(self, w):
        """
        TESTS::

            sage: LW33 = LyndonWords(3,3)
            sage: all([lw in LyndonWords() for lw in LW33])
            True
        """
        if isinstance(w, list):
            w = self._words(w)
        return w.is_lyndon()

class LyndonWords_evaluation(UniqueRepresentation, Parent):
    r"""
    The set of Lyndon words on a fixed multiset of letters.

    EXAMPLES::

        sage: L = LyndonWords([1,2,1])
        sage: L
        Lyndon words with evaluation [1, 2, 1]
        sage: L.list()
        [word: 1223, word: 1232, word: 1322]
    """
    def __init__(self, e):
        """
        TESTS::

            sage: LW21 = LyndonWords([2,1]); LW21
            Lyndon words with evaluation [2, 1]
            sage: LW21 == loads(dumps(LW21))
            True
        """
        self._e = e
        self._words = FiniteWords(len(e))

        from sage.categories.enumerated_sets import EnumeratedSets
        Parent.__init__(self,
                        category=EnumeratedSets().Finite(),
                        facade=(self._words,)
                        )

    def __repr__(self):
        """
        TESTS::

            sage: repr(LyndonWords([2,1,1]))
            'Lyndon words with evaluation [2, 1, 1]'
        """
        return "Lyndon words with evaluation %s"%self._e

    def __call__(self, *args, **kwds):
        r"""
        TESTS::

            sage: L = LyndonWords([1,2,1])
            sage: L([1,2,2,3])
            word: 1223
            sage: L([2,1,2,3])
            Traceback (most recent call last):
            ...
            ValueError: not a Lyndon word
            sage: L([1,2])
            Traceback (most recent call last):
            ...
            ValueError: evaluation is not [1, 2, 1]
        """
        w = self._words(*args, **kwds)
        if kwds.get('check', True) and not w.is_lyndon():
            raise ValueError("not a Lyndon word")
        if kwds.get('check', True) and w.evaluation() != self._e:
            raise ValueError("evaluation is not {}".format(self._e))
        return w

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: [1,2,1,2] in LyndonWords([2,2])
            False
            sage: [1,1,2,2] in LyndonWords([2,2])
            True
            sage: all([ lw in LyndonWords([2,1,3,1]) for lw in LyndonWords([2,1,3,1])])
            True
        """
        if isinstance(x, list):
            x = self._words(x)
        return x in self._words and x.is_lyndon() and x.evaluation() == self._e

    def cardinality(self):
        """
        Returns the number of Lyndon words with the evaluation e.

        EXAMPLES::

            sage: LyndonWords([]).cardinality()
            0
            sage: LyndonWords([2,2]).cardinality()
            1
            sage: LyndonWords([2,3,2]).cardinality()
            30

        Check to make sure that the count matches up with the number of
        Lyndon words generated.

        ::

            sage: comps = [[],[2,2],[3,2,7],[4,2]]+Compositions(4).list()
            sage: lws = [ LyndonWords(comp) for comp in comps]
            sage: all( [ lw.cardinality() == len(lw.list()) for lw in lws] )
            True
        """
        evaluation = self._e
        le = builtins.list(evaluation)
        if len(evaluation) == 0:
            return 0

        n = sum(evaluation)

        return sum([moebius(j)*factorial(n/j) / prod([factorial(ni/j) for ni in evaluation]) for j in divisors(gcd(le))])/n

    def __iter__(self):
        """
        An iterator for the Lyndon words with evaluation e.

        EXAMPLES::

            sage: LyndonWords([1]).list()    #indirect doctest
            [word: 1]
            sage: LyndonWords([2]).list()    #indirect doctest
            []
            sage: LyndonWords([3]).list()    #indirect doctest
            []
            sage: LyndonWords([3,1]).list()  #indirect doctest
            [word: 1112]
            sage: LyndonWords([2,2]).list()  #indirect doctest
            [word: 1122]
            sage: LyndonWords([1,3]).list()  #indirect doctest
            [word: 1222]
            sage: LyndonWords([3,3]).list()  #indirect doctest
            [word: 111222, word: 112122, word: 112212]
            sage: LyndonWords([4,3]).list()  #indirect doctest
            [word: 1111222, word: 1112122, word: 1112212, word: 1121122, word: 1121212]

        TESTS:

        Check that :trac:`12997` is fixed::

            sage: LyndonWords([0,1]).list()
            [word: 2]
            sage: LyndonWords([0,2]).list()
            []
            sage: LyndonWords([0,0,1,0,1]).list()
            [word: 35]
        """
        if not self._e:
            return

        k = 0
        while self._e[k] == 0:
            k += 1

        for z in necklace._sfc(self._e[k:], equality=True):
            yield self._words([i+k+1 for i in z], check=False)

class LyndonWords_nk(UniqueRepresentation, Parent):
    r"""
    Lyndon words of fixed length `n` over the alphabet `{1, 2, ..., k}`.

    EXAMPLES::

        sage: L = LyndonWords(3, 4)
        sage: L.list()
        [word: 1112,
         word: 1113,
         word: 1122,
         word: 1123,
         ...
         word: 1333,
         word: 2223,
         word: 2233,
         word: 2333]
    """
    def __init__(self, n, k):
        """
        INPUT:

        - ``n`` -- the length of the words

        - ``k`` -- the size of the alphabet

        TESTS::

            sage: LW23 = LyndonWords(2,3); LW23
            Lyndon words from an alphabet of size 2 of length 3
            sage: LW23== loads(dumps(LW23))
            True
        """
        self._n = n
        self._k = k
        self._words = FiniteWords(self._n)

        from sage.categories.enumerated_sets import EnumeratedSets
        Parent.__init__(self,
                        category=EnumeratedSets().Finite(),
                        facade=(self._words,)
                        )

    def __repr__(self):
        """
        TESTS::

            sage: repr(LyndonWords(2, 3))
            'Lyndon words from an alphabet of size 2 of length 3'
        """
        return "Lyndon words from an alphabet of size %s of length %s"%(self._n, self._k)

    def __call__(self, *args, **kwds):
        r"""
        TESTS::

            sage: L = LyndonWords(3,3)
            sage: L([1,2,3])
            word: 123
            sage: L([2,3,4])
            Traceback (most recent call last):
            ...
            ValueError: 4 not in alphabet!
            sage: L([2,1,3])
            Traceback (most recent call last):
            ...
            ValueError: not a Lyndon word
            sage: L([1,2,2,3,3])
            Traceback (most recent call last):
            ...
            ValueError: length is not n=3
        """
        w = self._words(*args, **kwds)
        if kwds.get('check', True) and not w.is_lyndon():
            raise ValueError("not a Lyndon word")
        if w.length() != self._n:
            raise ValueError("length is not n={}".format(self._n))
        return w

    def __contains__(self, w):
        """
        TESTS::

            sage: LW33 = LyndonWords(3,3)
            sage: all([lw in LW33 for lw in LW33])
            True
        """
        if isinstance(w, list):
            w = self._words(w)
        return w in self._words and w.length() == self._k and len(set(w)) <= self._n

    def cardinality(self):
        """
        TESTS::

            sage: [ LyndonWords(3,i).cardinality() for i in range(1, 11) ]
            [3, 3, 8, 18, 48, 116, 312, 810, 2184, 5880]
        """
        if self._k == 0:
            return Integer(1)
        else:
            s = Integer(0)
            for d in divisors(self._k):
                s += moebius(d)*(self._n**(self._k/d))
        return s/self._k

    def __iter__(self):
        """
        TESTS::

            sage: LyndonWords(3,3).list() # indirect doctest
            [word: 112, word: 113, word: 122, word: 123, word: 132, word: 133, word: 223, word: 233]
        """
        for c in IntegerVectors(self._k, self._n):
            cf = []
            nonzero_indices = []
            for i,x in enumerate(c):
                if x:
                    nonzero_indices.append(i)
                    cf.append(x)
            for lw in LyndonWords_evaluation(Composition(cf)):
                yield self._words([nonzero_indices[x-1]+1 for x in lw], check=False)

def StandardBracketedLyndonWords(n, k):
    """
    Returns the combinatorial class of standard bracketed Lyndon words
    from [1, ..., n] of length k. These are in one to one
    correspondence with the Lyndon words and form a basis for the
    subspace of degree k of the free Lie algebra of rank n.

    EXAMPLES::

        sage: SBLW33 = StandardBracketedLyndonWords(3,3); SBLW33
        Standard bracketed Lyndon words from an alphabet of size 3 of length 3
        sage: SBLW33.first()
        [1, [1, 2]]
        sage: SBLW33.last()
        [[2, 3], 3]
        sage: SBLW33.cardinality()
        8
        sage: SBLW33.random_element()
        [1, [1, 2]]
    """
    return StandardBracketedLyndonWords_nk(n,k)

class StandardBracketedLyndonWords_nk(UniqueRepresentation, Parent):
    def __init__(self, n, k):
        """
        TESTS::

            sage: SBLW = StandardBracketedLyndonWords(3, 2)
            sage: SBLW == loads(dumps(SBLW))
            True
        """
        self._n = n
        self._k = k
        self._lyndon = LyndonWords(self._n, self._k)

        from sage.categories.enumerated_sets import EnumeratedSets
        Parent.__init__(self, category=EnumeratedSets().Finite())

    def __repr__(self):
        """
        TESTS::

            sage: repr(StandardBracketedLyndonWords(3, 3))
            'Standard bracketed Lyndon words from an alphabet of size 3 of length 3'
        """
        return "Standard bracketed Lyndon words from an alphabet of size %s of length %s"%(self._n, self._k)

    def cardinality(self):
        """
        EXAMPLES::

            sage: StandardBracketedLyndonWords(3, 3).cardinality()
            8
            sage: StandardBracketedLyndonWords(3, 4).cardinality()
            18
        """
        return self._lyndon.cardinality()

    def __call__(self, *args, **kwds):
        r"""
        EXAMPLES::

            sage: S = StandardBracketedLyndonWords(3, 3)
            sage: S([1,2,3])
            [1, [2, 3]]
        """
        return standard_bracketing(self._lyndon(*args, **kwds))

    def __iter__(self):
        """
        EXAMPLES::

            sage: StandardBracketedLyndonWords(3, 3).list()
            [[1, [1, 2]],
             [1, [1, 3]],
             [[1, 2], 2],
             [1, [2, 3]],
             [[1, 3], 2],
             [[1, 3], 3],
             [2, [2, 3]],
             [[2, 3], 3]]
        """
        from itertools import imap
        return imap(standard_bracketing, self._lyndon)

def standard_bracketing(lw):
    """
    Returns the standard bracketing of a Lyndon word lw.

    EXAMPLES::

        sage: import sage.combinat.lyndon_word as lyndon_word
        sage: map( lyndon_word.standard_bracketing, LyndonWords(3,3) )
        [[1, [1, 2]],
         [1, [1, 3]],
         [[1, 2], 2],
         [1, [2, 3]],
         [[1, 3], 2],
         [[1, 3], 3],
         [2, [2, 3]],
         [[2, 3], 3]]
    """
    if len(lw) == 1:
        return lw[0]

    for i in range(1,len(lw)):
        if lw[i:] in LyndonWords():
            return [ standard_bracketing( lw[:i] ), standard_bracketing(lw[i:]) ]
