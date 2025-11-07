"""
Lyndon words
"""
# ****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.arith.misc import divisors, gcd, moebius, multinomial
from sage.combinat.combinat_cython import lyndon_word_iterator
from sage.combinat.composition import Composition, Compositions
from sage.combinat.necklace import _sfc
from sage.combinat.words.finite_word import FiniteWord_class
from sage.combinat.words.words import FiniteWords
from sage.rings.integer import Integer
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation


def LyndonWords(e=None, k=None):
    """
    Return the combinatorial class of Lyndon words.

    A Lyndon word `w` is a word that is lexicographically less than all of
    its rotations.  Equivalently, whenever `w` is split into two non-empty
    substrings, `w` is lexicographically less than the right substring.

    See :wikipedia:`Lyndon_word`

    INPUT:

    - no input at all

    or

    - ``e`` -- integer; size of alphabet
    - ``k`` -- integer; length of the words

    or

    - ``e`` -- a composition

    OUTPUT: a combinatorial class of Lyndon words

    EXAMPLES::

        sage: LyndonWords()
        Lyndon words

    If e is an integer, then e specifies the length of the
    alphabet; k must also be specified in this case::

        sage: LW = LyndonWords(3, 4); LW
        Lyndon words from an alphabet of size 3 of length 4
        sage: LW.first()
        word: 1112
        sage: LW.last()
        word: 2333
        sage: LW.random_element()  # random                                             # needs sage.libs.pari
        word: 1232
        sage: LW.cardinality()                                                          # needs sage.libs.pari
        18

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
                raise TypeError("k must be a nonnegative integer")
            if k < 0:
                raise TypeError("k must be a nonnegative integer")
            return LyndonWords_nk(Integer(e), Integer(k))
    elif e in Compositions():
        return LyndonWords_evaluation(Composition(e))

    raise TypeError("e must be a positive integer or a composition")


def LyndonWord(data, check=True):
    r"""
    Construction of a Lyndon word.

    INPUT:

    - ``data`` -- list
    - ``check`` -- boolean (default: ``True``); if ``True``,
      check that the input data represents a Lyndon word

    OUTPUT: a Lyndon word

    EXAMPLES::

        sage: LyndonWord([1,2,2])
        word: 122
        sage: LyndonWord([1,2,3])
        word: 123
        sage: LyndonWord([2,1,2,3])
        Traceback (most recent call last):
        ...
        ValueError: not a Lyndon word

    If ``check`` is ``False``, then no verification is done::

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
            sage: all(lw in LyndonWords() for lw in LW33)
            True
        """
        if isinstance(w, list):
            w = self._words(w, check=False)
        return isinstance(w, FiniteWord_class) and w.is_lyndon()


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
        return "Lyndon words with evaluation %s" % self._e

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

    def __contains__(self, w):
        """
        EXAMPLES::

            sage: [1,2,1,2] in LyndonWords([2,2])
            False
            sage: [1,1,2,2] in LyndonWords([2,2])
            True
            sage: all(lw in LyndonWords([2,1,3,1]) for lw in LyndonWords([2,1,3,1]))
            True
        """
        if isinstance(w, list):
            w = self._words(w, check=False)
        if isinstance(w, FiniteWord_class) and all(x in self._words.alphabet() for x in w):
            ev_dict = w.evaluation_dict()
            evaluation = [ev_dict.get(x, 0) for x in self._words.alphabet()]
            return evaluation == self._e and w.is_lyndon()
        else:
            return False

    def cardinality(self):
        """
        Return the number of Lyndon words with the evaluation e.

        EXAMPLES::

            sage: LyndonWords([]).cardinality()
            0
            sage: LyndonWords([2,2]).cardinality()                                      # needs sage.libs.pari
            1
            sage: LyndonWords([2,3,2]).cardinality()                                    # needs sage.libs.pari
            30

        Check to make sure that the count matches up with the number of
        Lyndon words generated::

            sage: comps = [[],[2,2],[3,2,7],[4,2]] + Compositions(4).list()
            sage: lws = [LyndonWords(comp) for comp in comps]
            sage: all(lw.cardinality() == len(lw.list()) for lw in lws)                 # needs sage.libs.pari
            True
        """
        evaluation = self._e
        le = list(evaluation)
        if not evaluation:
            return Integer(0)
        n = sum(evaluation)
        return sum(moebius(j) * multinomial([ni // j for ni in evaluation])
                   for j in divisors(gcd(le))) // n

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

        Check that :issue:`12997` is fixed::

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
        for z in _sfc(self._e[k:], equality=True):
            yield self._words([i + k + 1 for i in z], check=False)


class LyndonWords_nk(UniqueRepresentation, Parent):
    r"""
    Lyndon words of fixed length `k` over the alphabet `\{1, 2, \ldots, n\}`.

    INPUT:

    - ``n`` -- the size of the alphabet
    - ``k`` -- the length of the words

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
        Initialize ``self``.

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
        return "Lyndon words from an alphabet of size %s of length %s" % (self._n, self._k)

    def __call__(self, *args, **kwds):
        r"""
        TESTS::

            sage: L = LyndonWords(3,3)
            sage: L([1,2,3])
            word: 123
            sage: L([2,3,4])
            Traceback (most recent call last):
            ...
            ValueError: 4 not in alphabet
            sage: L([2,1,3])
            Traceback (most recent call last):
            ...
            ValueError: not a Lyndon word
            sage: L([1,2,2,3,3])
            Traceback (most recent call last):
            ...
            ValueError: length is not k=3

        Make sure that the correct length is checked (:issue:`30186`)::

            sage: L = LyndonWords(2, 4)
            sage: _ = L(L.random_element())                                             # needs sage.libs.pari
        """
        w = self._words(*args, **kwds)
        if kwds.get('check', True) and not w.is_lyndon():
            raise ValueError("not a Lyndon word")
        if kwds.get('check', True) and w.length() != self._k:
            raise ValueError("length is not k={}".format(self._k))
        return w

    def __contains__(self, w):
        """
        TESTS::

            sage: LW33 = LyndonWords(3,3)
            sage: all(lw in LW33 for lw in LW33)
            True
        """
        if isinstance(w, list):
            w = self._words(w, check=False)
        return isinstance(w, FiniteWord_class) and w.length() == self._k \
            and all(x in self._words.alphabet() for x in w) and w.is_lyndon()

    def cardinality(self):
        """
        TESTS::

            sage: [ LyndonWords(3,i).cardinality() for i in range(1, 11) ]              # needs sage.libs.pari
            [3, 3, 8, 18, 48, 116, 312, 810, 2184, 5880]
        """
        if self._k == 0:
            return Integer(1)
        else:
            s = Integer(0)
            for d in divisors(self._k):
                s += moebius(d) * self._n**(self._k // d)
        return s // self._k

    def __iter__(self):
        """
        TESTS::

            sage: LyndonWords(3,3).list()  # indirect doctest
            [word: 112, word: 113, word: 122, word: 123, word: 132, word: 133, word: 223, word: 233]

            sage: sum(1 for lw in LyndonWords(11, 6))
            295020

            sage: sum(1 for lw in LyndonWords(1000, 1))
            1000

            sage: sum(1 for lw in LyndonWords(1, 1000))
            0

            sage: list(LyndonWords(1, 1))                                               # needs sage.libs.pari
            [word: 1]
        """
        W = self._words._element_classes['list']
        for lw in lyndon_word_iterator(self._n, self._k):
            yield W(self._words, [i + 1 for i in lw])


def StandardBracketedLyndonWords(n, k):
    """
    Return the combinatorial class of standard bracketed Lyndon words
    from [1, ..., n] of length k.

    These are in one to one correspondence with the Lyndon words and
    form a basis for the subspace of degree k of the free Lie algebra
    of rank n.

    EXAMPLES::

        sage: SBLW33 = StandardBracketedLyndonWords(3,3); SBLW33
        Standard bracketed Lyndon words from an alphabet of size 3 of length 3
        sage: SBLW33.first()
        [1, [1, 2]]
        sage: SBLW33.last()
        [[2, 3], 3]
        sage: SBLW33.cardinality()
        8
        sage: SBLW33.random_element() in SBLW33
        True
    """
    return StandardBracketedLyndonWords_nk(n, k)


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
        return "Standard bracketed Lyndon words from an alphabet of size %s of length %s" % (self._n, self._k)

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

    def __contains__(self, sblw):
        """
        EXAMPLES::

            sage: S = StandardBracketedLyndonWords(2, 3)
            sage: [[1, 2], 2] in S
            True
            sage: [1, [2, 2]] in S
            False
            sage: [1, [2, 3]] in S
            False
            sage: [1, 2] in S
            False
        """
        try:
            lw = standard_unbracketing(sblw)
        except ValueError:
            return False
        return len(lw) == self._k and all(a in self._lyndon._words.alphabet() for a in lw.parent().alphabet())

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
        for x in self._lyndon:
            yield standard_bracketing(x)


def standard_bracketing(lw):
    """
    Return the standard bracketing of a Lyndon word ``lw``.

    EXAMPLES::

        sage: import sage.combinat.words.lyndon_word as lyndon_word
        sage: [lyndon_word.standard_bracketing(u) for u in LyndonWords(3,3)]
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

    for i in range(1, len(lw)):
        if lw[i:] in LyndonWords():
            return [standard_bracketing(lw[:i]), standard_bracketing(lw[i:])]


def standard_unbracketing(sblw):
    """
    Return flattened ``sblw`` if it is a standard bracketing of a Lyndon word,
    otherwise raise an error.

    EXAMPLES::

        sage: from sage.combinat.words.lyndon_word import standard_unbracketing
        sage: standard_unbracketing([1, [2, 3]])
        word: 123
        sage: standard_unbracketing([[1, 2], 3])
        Traceback (most recent call last):
        ...
        ValueError: not a standard bracketing of a Lyndon word

    TESTS::

        sage: standard_unbracketing(1) # Letters don't use brackets.
        word: 1
        sage: standard_unbracketing([1])
        Traceback (most recent call last):
        ...
        ValueError: not a standard bracketing of a Lyndon word
    """
    # Nested helper function that not only returns (flattened) w, but also its
    # right factor in the standard Lyndon factorization.
    def standard_unbracketing_rec(w):
        if not isinstance(w, list):
            return [w], []
        if len(w) != 2:
            raise ValueError("not a standard bracketing of a Lyndon word")
        x, t = standard_unbracketing_rec(w[0])
        y, _ = standard_unbracketing_rec(w[1])
        # If x = st is a standard Lyndon factorization, and y is a Lyndon word
        # such that y <= t, then xy is standard (but not necessarily Lyndon).
        if x < y and (len(t) == 0 or y <= t):
            x += y
            return x, y
        else:
            raise ValueError("not a standard bracketing of a Lyndon word")
    lw, _ = standard_unbracketing_rec(sblw)
    return FiniteWords(list(set(lw)))(lw, datatype='list', check=False)
