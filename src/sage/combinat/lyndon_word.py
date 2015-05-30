# -*- coding: utf-8 -*-
"""
Lyndon words
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

from combinat import CombinatorialClass
from sage.combinat.composition import Composition, Compositions
from sage.rings.all import divisors, gcd, moebius, Integer
from sage.rings.arith import factorial
from sage.misc.all import prod
import __builtin__
import necklace
from integer_vector import IntegerVectors

from sage.combinat.words.word import FiniteWord_list
from sage.combinat.words.words import Words_all, FiniteWords_length_k_over_OrderedAlphabet
from sage.combinat.words.alphabet import build_alphabet

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
            return LyndonWords_nk(e, k)
    elif e in Compositions():
        return LyndonWords_evaluation(e)

    raise TypeError("e must be a positive integer or a composition")

class LyndonWord(FiniteWord_list):
    def __init__(self, data, check=True):
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
            ValueError: Not a Lyndon word

        If check is False, then no verification is done::

            sage: LyndonWord([2,1,2,3], check=False)
            word: 2123
        """
        super(LyndonWord,self).__init__(parent=LyndonWords(),data=data)
        if check and not self.is_lyndon():
            raise ValueError("Not a Lyndon word")

class LyndonWords_class(Words_all):
    def __repr__(self):
        r"""
        String representation.

        EXAMPLES::

            sage: LyndonWords()
            Lyndon words
        """
        return "Lyndon words"

    def __contains__(self, x):
        """
        TESTS::

            sage: LW33 = LyndonWords(3,3)
            sage: all([lw in LyndonWords() for lw in LW33])
            True
        """
        try:
            return LyndonWord(x)
        except ValueError:
            return False

class LyndonWords_evaluation(CombinatorialClass):
    def __init__(self, e):
        """
        TESTS::

            sage: LW21 = LyndonWords([2,1]); LW21
            Lyndon words with evaluation [2, 1]
            sage: LW21 == loads(dumps(LW21))
            True
        """
        self.e = Composition(e)

    def __repr__(self):
        """
        TESTS::

            sage: repr(LyndonWords([2,1,1]))
            'Lyndon words with evaluation [2, 1, 1]'
        """
        return "Lyndon words with evaluation %s"%self.e

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
        try:
            w = LyndonWord(x)
        except ValueError:
            return False
        ed = w.evaluation_dict()
        for (i,a) in enumerate(self.e):
            if ed[i+1] != a:
                return False
        return True

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
        evaluation = self.e
        le = __builtin__.list(evaluation)
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
        if self.e == []:
            return

        k = 0
        while self.e[k] == 0:
            k += 1

        for z in necklace._sfc(self.e[k:], equality=True):
            yield LyndonWord([i+k+1 for i in z], check=False)

class LyndonWords_nk(FiniteWords_length_k_over_OrderedAlphabet):
    def __init__(self, n, k):
        """
        TESTS::

            sage: LW23 = LyndonWords(2,3); LW23
            Lyndon words from an alphabet of size 2 of length 3
            sage: LW23== loads(dumps(LW23))
            True
        """
        self.n = Integer(n)
        self.k = Integer(k)
        alphabet = build_alphabet(range(1,self.n+1))
        super(LyndonWords_nk,self).__init__(alphabet,self.k)

    def __repr__(self):
        """
        TESTS::

            sage: repr(LyndonWords(2, 3))
            'Lyndon words from an alphabet of size 2 of length 3'
        """
        return "Lyndon words from an alphabet of size %s of length %s"%(self.n, self.k)

    def __contains__(self, x):
        """
        TESTS::

            sage: LW33 = LyndonWords(3,3)
            sage: all([lw in LW33 for lw in LW33])
            True
        """
        try:
            w = LyndonWord(x)
            return w.length() == self.k and len(set(w)) <= self.n
        except ValueError:
            return False

    def cardinality(self):
        """
        TESTS::

            sage: [ LyndonWords(3,i).cardinality() for i in range(1, 11) ]
            [3, 3, 8, 18, 48, 116, 312, 810, 2184, 5880]
        """
        if self.k == 0:
            return 1
        else:
            s = 0
            for d in divisors(self.k):
                s += moebius(d)*(self.n**(self.k/d))
        return s/self.k

    def __iter__(self):
        """
        TESTS::

            sage: LyndonWords(3,3).list() # indirect doctest
            [word: 112, word: 113, word: 122, word: 123, word: 132, word: 133, word: 223, word: 233]
        """
        for c in IntegerVectors(self.k, self.n):
            cf = [x for x in c if x != 0]
            nonzero_indices = []
            for i in range(len(c)):
                if c[i] != 0:
                    nonzero_indices.append(i)
            for lw in LyndonWords_evaluation(cf):
                yield LyndonWord(map(lambda x: nonzero_indices[x-1]+1, lw), check=False)

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

class StandardBracketedLyndonWords_nk(CombinatorialClass):
    def __init__(self, n, k):
        """
        TESTS::

            sage: SBLW = StandardBracketedLyndonWords(3, 2)
            sage: SBLW == loads(dumps(SBLW))
            True
        """
        self.n = Integer(n)
        self.k = Integer(k)

    def __repr__(self):
        """
        TESTS::

            sage: repr(StandardBracketedLyndonWords(3, 3))
            'Standard bracketed Lyndon words from an alphabet of size 3 of length 3'
        """
        return "Standard bracketed Lyndon words from an alphabet of size %s of length %s"%(self.n, self.k)

    def cardinality(self):
        """
        EXAMPLES::

            sage: StandardBracketedLyndonWords(3, 3).cardinality()
            8
            sage: StandardBracketedLyndonWords(3, 4).cardinality()
            18
        """
        return LyndonWords_nk(self.n,self.k).cardinality()

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
        for lw in LyndonWords_nk(self.n, self.k):
            yield standard_bracketing(lw)

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
