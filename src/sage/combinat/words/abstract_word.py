# sage_setup: distribution = sagemath-combinat
r"""
Abstract word (finite or infinite)

This module gathers functions that works for both finite and infinite
words.

AUTHORS:

- Sébastien Labbé
- Franco Saliola

EXAMPLES::

    sage: a = 0.618
    sage: g = words.CodingOfRotationWord(alpha=a, beta=1-a, x=a)
    sage: f = words.FibonacciWord()
    sage: p = f.longest_common_prefix(g, length='finite')
    sage: p
    word: 0100101001001010010100100101001001010010...
    sage: p.length()
    231
"""
# ****************************************************************************
#       Copyright (C) 2008-2010 Sebastien Labbe <slabqc@gmail.com>,
#                     2008-2010 Franco Saliola <saliola@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from itertools import islice, groupby

from sage.structure.sage_object import SageObject
from sage.combinat.words.word_options import word_options
from sage.rings.finite_rings.integer_mod_ring import Integers
from sage.rings.infinity import Infinity
from sage.rings.integer_ring import ZZ
from sage.structure.richcmp import richcmp_method, rich_to_bool, richcmp_item


@richcmp_method
class Word_class(SageObject):
    def parent(self):
        r"""
        Return the parent of ``self``.

        TESTS::

            sage: Word(iter([1,2,3]), length='unknown').parent()
            Finite words over Set of Python objects of class 'object'
            sage: Word(range(12)).parent()
            Finite words over Set of Python objects of class 'object'
            sage: Word(range(4), alphabet=list(range(6))).parent()
            Finite words over {0, 1, 2, 3, 4, 5}
            sage: Word(iter('abac'), alphabet='abc').parent()
            Finite words over {'a', 'b', 'c'}
        """
        return self._parent

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: Word(iter([1,2,3]), length='unknown')._repr_()
            'word: 123'
            sage: Word(range(100), length='unknown')._repr_()
            'word: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,...'
            sage: Word(lambda x:x%3)._repr_()
            'word: 0120120120120120120120120120120120120120...'
        """
        if word_options['old_repr']:
            return "Word over %s" % (str(self.parent().alphabet())[17:])
        return word_options['identifier'] + self.string_rep()

    def string_rep(self):
        r"""
        Return the (truncated) raw sequence of letters as a string.

        EXAMPLES::

            sage: Word('abbabaab').string_rep()
            'abbabaab'
            sage: Word([0, 1, 0, 0, 1]).string_rep()
            '01001'
            sage: Word([0,1,10,101]).string_rep()
            '0,1,10,101'
            sage: WordOptions(letter_separator='-')
            sage: Word([0,1,10,101]).string_rep()
            '0-1-10-101'
            sage: WordOptions(letter_separator=',')

        TESTS:

        Insertion in a str::

            sage: from itertools import count
            sage: w = Word((i % 5 for i in count()), length='unknown')
            sage: "w = %s in this string." % w
            'w = 0123401234012340123401234012340123401234... in this string.'

        Using LatexExpr::

            sage: from sage.misc.latex import LatexExpr
            sage: LatexExpr(w)
            0123401234012340123401234012340123401234...

        With the print statement::

            sage: print(w)
            0123401234012340123401234012340123401234...

        Truncation is done for possibly infinite words::

            sage: print(w)
            0123401234012340123401234012340123401234...
        """
        l = word_options['truncate_length']
        letters = list(islice(self, int(l + 1)))
        if len(letters) == l + 1:
            letters.pop()
            suffix = "..."
        else:
            suffix = ""
        if word_options['display'] == 'string':
            ls = word_options['letter_separator']
            letters = [str(a) for a in letters]
            if all(len(a) == 1 for a in letters):
                return ''.join(letters) + suffix
            elif suffix == "...":
                return ls.join(letters) + ls + suffix
            else:
                return ls.join(letters)
        elif word_options['display'] == 'list':
            if suffix == "...":
                return "[%s, %s]" % (str(letters)[1:-1], suffix)
            else:
                return str(letters)

    __str__ = string_rep

    def __iter__(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.words.word import Word_class
            sage: w = Word_class()
            sage: w.__iter__()
            Traceback (most recent call last):
            ...
            NotImplementedError: you need to define an iterator in __iter__
        """
        raise NotImplementedError("you need to define an iterator in __iter__")

    def length(self):
        r"""
        Return the length of ``self``.

        TESTS::

            sage: from sage.combinat.words.word import Word_class
            sage: w = Word(iter('abba'*100), length='unknown')
            sage: w.length() is None
            True
            sage: w = Word(iter('abba'), length='finite')
            sage: w.length()
            4
            sage: w = Word(iter([0,1,1,0,1,0,0,1]*100), length='unknown')
            sage: w.length() is None
            True
            sage: w = Word(iter([0,1,1,0,1,0,0,1]), length='finite')
            sage: w.length()
            8
        """
        return self._len

    def is_finite(self):
        r"""
        Return whether this word is known to be finite.

        .. WARNING::

            A word defined by an iterator such that its end has
            never been reached will returns False.

        EXAMPLES::

            sage: Word([]).is_finite()
            True
            sage: Word('a').is_finite()
            True
            sage: TM = words.ThueMorseWord()
            sage: TM.is_finite()
            False

        ::

            sage: w = Word(iter('a'*100))
            sage: w.is_finite()
            False
        """
        return False

    def __len__(self):
        r"""
        Return the length of ``self`` (as a Python integer).

        .. NOTE::

            For infinite words or words of unknown length,
            use `length()` method instead.

        OUTPUT: positive integer

        EXAMPLES::

            sage: len(Word(lambda n:n, length=1000))
            1000
            sage: len(Word(iter('a'*200), length='finite'))
            200

        We make sure :issue:`8574` is fixed::

            sage: s = WordMorphism('0->000,1->%s'%('1'*100))
            sage: len(s('1'))
            100

        For infinite words::

            sage: len(Word(lambda n:n))
            Traceback (most recent call last):
            ...
            TypeError: Python len method cannot return a non integer value (=+Infinity): use length method instead.
            sage: len(Word(iter('a'*200)))
            Traceback (most recent call last):
            ...
            TypeError: Python len method cannot return a non integer value (=None): use length method instead.

        For words of unknown length::

            sage: len(Word(iter('a'*200), length='unknown'))
            Traceback (most recent call last):
            ...
            TypeError: Python len method cannot return a non integer value (=None): use length method instead.
        """
        L = self.length()
        if L is None or L is Infinity:
            msg = "Python len method cannot return a non integer value (=%s): " % L
            msg += "use length method instead."
            raise TypeError(msg)
        return int(L)

    def __richcmp__(self, other, op):
        r"""
        Compare two words lexicographically according to the ordering
        defined by the parent of ``self``.

        This corresponds to Python's built-in
        ordering when no parent nor alphabet was used to defined the word.

        Provides for all normal comparison operators.

        .. NOTE::

            This function will not terminate if ``self`` and ``other``
            are equal infinite words!

        EXAMPLES::

            sage: W = Word
            sage: from itertools import count
            sage: W(range(1,10)) > W(range(10))
            True
            sage: W(range(10)) < W(range(1,10))
            True
            sage: W(range(10)) == W(range(10))
            True
            sage: W(range(10)) < W(count())
            True
            sage: W(count()) > W(range(10))
            True

        ::

            sage: W = Words(['a', 'b', 'c'])
            sage: W('a') > W([])
            True
            sage: W([]) < W('a')
            True

            sage: Word('abc') == Word(['a','b','c'])
            True
            sage: Words([0,1])([0,1,0,1]) == Words([0,1])([0,1,0,1])
            True
            sage: Words([0,1])([0,1,0,1]) == Words([0,1])([0,1,0,0])
            False

        It works even when parents are not the same::

            sage: Words([0,1,2])([0,1,0,1]) ==  Words([0,1])([0,1,0,1])
            True
            sage: Word('ababa') == Words('abcd')('ababa')
            True

        Or when one word is finite while the other is infinite::

            sage: Word(range(20)) == Word(lambda n:n)
            False
            sage: Word(lambda n:n) == Word(range(20))
            False

        Beware the following does not halt! ::

            sage: from itertools import count
            sage: Word(lambda n:n) == Word(count()) #not tested

        Examples for unequality::

            sage: w = Word(range(10))
            sage: z = Word(range(10))
            sage: w != z
            False
            sage: u = Word(range(12))
            sage: u != w
            True

        TESTS::

            sage: Word(count())[:20] == Word(range(20))
            True
            sage: Word(range(20)) == Word(count())[:20]
            True
            sage: Word(range(20)) == Word(lambda n:n)[:20]
            True
            sage: Word(range(20)) == Word(lambda n:n,length=20)
            True
        """
        if not isinstance(other, Word_class):
            return NotImplemented
        self_it, other_it = iter(self), iter(other)
        cmp_key = self._parent.sortkey_letters
        while True:
            try:
                cs = next(self_it)
            except StopIteration:
                try:
                    next(other_it)
                except StopIteration:
                    # If both self_it and other_it are exhausted then
                    # self == other. Return 0.
                    return rich_to_bool(op, 0)
                else:
                    # If self_it is exhausted, but not other_it, then
                    # self is a proper prefix of other: return -1
                    return rich_to_bool(op, -1)
            else:
                try:
                    co = next(other_it)
                except StopIteration:
                    # If self_it is not exhausted but other_it is, then
                    # other is a proper prefix of self: return 1.
                    return rich_to_bool(op, 1)
                else:
                    key_cs = cmp_key(cs)
                    key_co = cmp_key(co)
                    res = richcmp_item(key_cs, key_co, op)
                    if res is not NotImplemented:
                        return res

    def _longest_common_prefix_iterator(self, other):
        r"""
        Return an iterator of the longest common prefix of ``self`` and ``other``.

        INPUT:

        - ``other`` -- word

        OUTPUT: iterator

        EXAMPLES::

            sage: f = words.FibonacciWord()
            sage: it = f._longest_common_prefix_iterator(f)
            sage: w = Word(it, length='unknown'); w
            word: 0100101001001010010100100101001001010010...
            sage: w[:6]
            word: 010010
            sage: it = w._longest_common_prefix_iterator(w[:10])
            sage: w = Word(it, length='finite'); w
            word: 0100101001
        """
        for (b, c) in zip(self, other):
            if b == c:
                yield b
            else:
                break

    def longest_common_prefix(self, other, length='unknown'):
        r"""
        Return the longest common prefix of ``self`` and ``other``.

        INPUT:

        - ``other`` -- word

        - ``length`` -- string (default: ``'unknown'``)
          the length type of the resulting word if known. It may be one of
          the following:

          - ``'unknown'``
          - ``'finite'``
          - ``'infinite'``

        EXAMPLES::

            sage: f = lambda n : add(Integer(n).digits(2)) % 2
            sage: t = Word(f)
            sage: u = t[:10]
            sage: t.longest_common_prefix(u)
            word: 0110100110

        The longest common prefix of two equal infinite words::

            sage: t1 = Word(f)
            sage: t2 = Word(f)
            sage: t1.longest_common_prefix(t2)
            word: 0110100110010110100101100110100110010110...

        Useful to study the approximation of an infinite word::

            sage: a = 0.618
            sage: g = words.CodingOfRotationWord(alpha=a, beta=1-a, x=a)
            sage: f = words.FibonacciWord()
            sage: p = f.longest_common_prefix(g, length='finite')
            sage: p.length()
            231

        TESTS::

            sage: w = Word('12345')
            sage: y = Word('1236777')
            sage: w.longest_common_prefix(y)
            word: 123
            sage: w.longest_common_prefix(w)
            word: 12345
            sage: y.longest_common_prefix(w)
            word: 123
            sage: y.longest_common_prefix(y)
            word: 1236777
            sage: Word().longest_common_prefix(w)
            word:
            sage: w.longest_common_prefix(Word())
            word:
            sage: w.longest_common_prefix(w[:3])
            word: 123
            sage: Word("11").longest_common_prefix(Word("1"))
            word: 1
            sage: Word("1").longest_common_prefix(Word("11"))
            word: 1

        With infinite words::

            sage: t = words.ThueMorseWord('ab')
            sage: u = t[:10]
            sage: u.longest_common_prefix(t)
            word: abbabaabba
            sage: u.longest_common_prefix(u)
            word: abbabaabba

        Check length::

            sage: w1 = Word(iter('ab'*200))
            sage: w2 = Word(iter('bcd'*200))
            sage: w1.longest_common_prefix(w2, length=19)
            Traceback (most recent call last):
            ...
            ValueError: invalid argument length (=19)
        """
        it = self._longest_common_prefix_iterator(other)

        if length == "finite" or \
           (length == "unknown" and (self.is_finite() or other.is_finite())):
            parent = self._parent.factors()
        elif length == "infinite":
            parent = self._parent.shift()
        elif length == "unknown":
            from sage.combinat.words.words import FiniteOrInfiniteWords
            parent = FiniteOrInfiniteWords(self._parent.alphabet())
        else:
            raise ValueError("invalid argument length (={})".format(length))

        return parent(it)

    def _longest_periodic_prefix_iterator(self, period=1):
        r"""
        Return an iterator of the longest prefix of ``self`` having the given
        period.

        INPUT:

        - ``period`` -- positive integer (default: 1)

        OUTPUT: iterator

        EXAMPLES::

            sage: list(Word([])._longest_periodic_prefix_iterator())
            []
            sage: list(Word([1])._longest_periodic_prefix_iterator())
            [1]
            sage: list(Word([1,2])._longest_periodic_prefix_iterator())
            [1]
            sage: list(Word([1,1,2])._longest_periodic_prefix_iterator())
            [1, 1]
            sage: list(Word([1,1,1,2])._longest_periodic_prefix_iterator())
            [1, 1, 1]
            sage: Word(Word(lambda n:0)._longest_periodic_prefix_iterator())
            word: 0000000000000000000000000000000000000000...
            sage: list(Word([1,2,1,2,1,3])._longest_periodic_prefix_iterator(2))
            [1, 2, 1, 2, 1]
        """
        for i, l in enumerate(self):
            if self[i % period] == l:
                yield l
            else:
                break

    def longest_periodic_prefix(self, period=1):
        r"""
        Return the longest prefix of ``self`` having the given period.

        INPUT:

        - ``period`` -- positive integer (default: 1)

        OUTPUT: word

        EXAMPLES::

            sage: Word([]).longest_periodic_prefix()
            word:
            sage: Word([1]).longest_periodic_prefix()
            word: 1
            sage: Word([1,2]).longest_periodic_prefix()
            word: 1
            sage: Word([1,1,2]).longest_periodic_prefix()
            word: 11
            sage: Word([1,2,1,2,1,3]).longest_periodic_prefix(2)
            word: 12121
            sage: type(_)
            <class 'sage.combinat.words.word.FiniteWord_iter_with_caching'>
            sage: Word(lambda n:0).longest_periodic_prefix()
            word: 0000000000000000000000000000000000000000...
        """
        if self.is_finite():
            parent = self._parent.factors()
        else:
            from sage.combinat.words.words import FiniteOrInfiniteWords
            parent = FiniteOrInfiniteWords(self._parent.alphabet())
        return parent(self._longest_periodic_prefix_iterator(period))

    def is_empty(self):
        r"""
        Return ``True`` if the length of ``self`` is zero, and ``False`` otherwise.

        EXAMPLES::

            sage: it = iter([])
            sage: Word(it).is_empty()
            True
            sage: it = iter([1,2,3])
            sage: Word(it).is_empty()
            False
            sage: from itertools import count
            sage: Word(count()).is_empty()
            False
        """
        try:
            next(iter(self))
            return False
        except StopIteration:
            return True

    def _to_integer_iterator(self, use_parent_alphabet=False):
        r"""
        Return an iterator over the letters of an integer representation of
        ``self``.

        INPUT:

        - ``use_parent_alphabet`` -- boolean (default: ``False``); when ``True``
          and if the ``self`` parent's alphabet is finite, it uses the index of
          the letters in the alphabet. Otherwise, the first letter occurring in
          ``self`` is mapped to zero, and every letter that hasn't yet occurred in
          the word is mapped to the next available integer.

        EXAMPLES::

            sage: from itertools import count
            sage: w = Word(count())
            sage: ir = w._to_integer_iterator()
            sage: [next(ir) for _ in range(10)]
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            sage: w = Word(iter("abbacabba"))
            sage: ir = w._to_integer_iterator()
            sage: list(ir)
            [0, 1, 1, 0, 2, 0, 1, 1, 0]

        ::

            sage: w = Words('abc')('abbccc')
            sage: list(w._to_integer_iterator(True))
            [0, 1, 1, 2, 2, 2]
            sage: w = Words('acb')('abbccc')
            sage: list(w._to_integer_iterator(True))
            [0, 2, 2, 1, 1, 1]
            sage: w = Words('xabc')('abbccc')
            sage: list(w._to_integer_iterator(True))
            [1, 2, 2, 3, 3, 3]
        """
        from sage.combinat.words.words import FiniteWords, InfiniteWords
        if use_parent_alphabet and\
                isinstance(self.parent(), (FiniteWords, InfiniteWords)):
            A = self.parent().alphabet()
            for letter in self:
                yield A.rank(letter)

        else:
            mapping = {}
            next_value = 0
            for letter in self:
                if letter not in mapping:
                    mapping[letter] = next_value
                    next_value += 1
                yield mapping[letter]

    def to_integer_word(self):
        r"""
        Return a word over the integers whose letters are those output by
        self._to_integer_iterator()

        EXAMPLES::

            sage: from itertools import count
            sage: w = Word(count()); w
            word: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,...
            sage: w.to_integer_word()
            word: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,...
            sage: w = Word(iter("abbacabba"), length='finite'); w
            word: abbacabba
            sage: w.to_integer_word()
            word: 011020110
            sage: w = Word(iter("abbacabba"), length='unknown'); w
            word: abbacabba
            sage: w.to_integer_word()
            word: 011020110
        """
        length = "unknown" if self._len is None else self._len
        from sage.combinat.words.word import Word
        return Word(self._to_integer_iterator(), length=length)

    def lex_less(self, other):
        r"""
        Return ``True`` if ``self`` is lexicographically less than ``other``.

        EXAMPLES::

            sage: w = Word([1,2,3])
            sage: u = Word([1,3,2])
            sage: v = Word([3,2,1])
            sage: w.lex_less(u)
            True
            sage: v.lex_less(w)
            False
            sage: a = Word("abba")
            sage: b = Word("abbb")
            sage: a.lex_less(b)
            True
            sage: b.lex_less(a)
            False

        For infinite words::

            sage: t = words.ThueMorseWord()
            sage: t.lex_less(t[:10])
            False
            sage: t[:10].lex_less(t)
            True
        """
        return self < other

    def lex_greater(self, other):
        r"""
        Return ``True`` if ``self`` is lexicographically greater than ``other``.

        EXAMPLES::

            sage: w = Word([1,2,3])
            sage: u = Word([1,3,2])
            sage: v = Word([3,2,1])
            sage: w.lex_greater(u)
            False
            sage: v.lex_greater(w)
            True
            sage: a = Word("abba")
            sage: b = Word("abbb")
            sage: a.lex_greater(b)
            False
            sage: b.lex_greater(a)
            True

        For infinite words::

            sage: t = words.ThueMorseWord()
            sage: t[:10].lex_greater(t)
            False
            sage: t.lex_greater(t[:10])
            True
        """
        return self > other

    def apply_morphism(self, morphism):
        r"""
        Return the word obtained by applying the morphism to ``self``.

        INPUT:

        - ``morphism`` -- can be an instance of WordMorphism, or
          anything that can be used to construct one

        EXAMPLES::

            sage: w = Word("ab")
            sage: d = {'a':'ab', 'b':'ba'}
            sage: w.apply_morphism(d)
            word: abba
            sage: w.apply_morphism(WordMorphism(d))
            word: abba

        ::

            sage: w = Word('ababa')
            sage: d = dict(a='ab', b='ba')
            sage: d
            {'a': 'ab', 'b': 'ba'}
            sage: w.apply_morphism(d)
            word: abbaabbaab

        For infinite words::

            sage: t = words.ThueMorseWord([0,1]); t
            word: 0110100110010110100101100110100110010110...
            sage: t.apply_morphism({0:8,1:9})
            word: 8998988998898998988989988998988998898998...
        """
        from sage.combinat.words.morphism import WordMorphism
        if not isinstance(morphism, WordMorphism):
            morphism = WordMorphism(morphism)
        return morphism(self)

    def _delta_iterator(self):
        r"""
        Return an iterator of the image of ``self`` under the delta morphism.
        This is the word composed of the length of consecutive runs of the
        same letter in a given word.

        OUTPUT: generator object

        EXAMPLES::

            sage: W = Words('0123456789')
            sage: it=W('22112122')._delta_iterator()
            sage: Word(it)
            word: 22112
            sage: Word(W('555008')._delta_iterator())
            word: 321
            sage: Word(W()._delta_iterator())
            word:

        For infinite words::

            sage: t = words.ThueMorseWord()
            sage: it = t._delta_iterator()
            sage: Word(it)
            word: 1211222112112112221122211222112112112221...
        """
        for letter, run in groupby(self):
            yield len(list(run))

    def delta(self):
        r"""
        Return the image of ``self`` under the delta morphism.

        This is the word composed of the length of consecutive runs of
        the same letter in a given word.

        OUTPUT: word over integers

        EXAMPLES:

        For finite words::

            sage: W = Words('0123456789')
            sage: W('22112122').delta()
            word: 22112
            sage: W('555008').delta()
            word: 321
            sage: W().delta()
            word:
            sage: Word('aabbabaa').delta()
            word: 22112

        For infinite words::

            sage: t = words.ThueMorseWord()
            sage: t.delta()
            word: 1211222112112112221122211222112112112221...
        """
        from sage.combinat.words.word import Word
        from sage.rings.semirings.non_negative_integer_semiring import NN
        return Word(self._delta_iterator(), alphabet=NN)

    def _iterated_right_palindromic_closure_iterator(self, f=None):
        r"""
        Return an iterator over the iterated (`f`-)palindromic closure of ``self``.

        INPUT:

        - ``f`` -- involution on the alphabet of ``self`` (default: ``None``);
          it must be callable on letters as well as words (e.g. WordMorphism)

        OUTPUT: iterator; the iterated (`f`-)palindromic closure of ``self``

        EXAMPLES::

            sage: w = Word('abc')
            sage: it = w._iterated_right_palindromic_closure_iterator()
            sage: Word(it)
            word: abacaba

        ::

            sage: w = Word('aaa')
            sage: it = w._iterated_right_palindromic_closure_iterator()
            sage: Word(it)
            word: aaa

        ::

            sage: w = Word('abbab')
            sage: it = w._iterated_right_palindromic_closure_iterator()
            sage: Word(it)
            word: ababaabababaababa

        An infinite word::

            sage: t = words.ThueMorseWord('ab')
            sage: it = t._iterated_right_palindromic_closure_iterator()
            sage: Word(it)
            word: ababaabababaababaabababaababaabababaabab...

        TESTS:

        The empty word::

            sage: w = Word()
            sage: it = w._iterated_right_palindromic_closure_iterator()
            sage: next(it)
            Traceback (most recent call last):
            ...
            StopIteration

        REFERENCES:

        -   [1] A. de Luca, A. De Luca, Pseudopalindrome closure operators
            in free monoids, Theoret. Comput. Sci. 362 (2006) 282--300.
        """
        par = self.parent().factors()
        w = self[:0]
        for letter in self:
            length_before = w.length()
            w = (w*par([letter])).palindromic_closure(f=f)
            length_after = w.length()
            d = length_after - length_before
            yield from w[-d:]

    def _iterated_right_palindromic_closure_recursive_iterator(self, f=None):
        r"""
        Return an iterator over the iterated (`f`-)palindromic closure of ``self``.

        INPUT:

        - ``f`` -- involution (default: ``None``) on the alphabet of ``self``;
          it must be callable on letters as well as words (e.g. WordMorphism)

        OUTPUT: iterator; the iterated (`f`-)palindromic closure of ``self``

        ALGORITHM:

            For the case of palindromes only, it has been shown in [2] that
            the iterated right palindromic closure of a given word `w`,
            denoted by `IRPC(w)`, may be obtained as follows.
            Let `w` be any word and `x` be a letter. Then

            #. If `x` does not occur in `w`,
               `IRPC(wx) = IRPC(w) \cdot x \cdot IRPC(w)`
            #. Otherwise, write `w = w_1xw_2` such that `x` does not
               occur in `w_2`. Then `IRPC(wx) = IRPC(w) \cdot IRPC(w_1)^{-1}
               \cdot IRPC(w)`

            This formula is directly generalized to the case of `f`-palindromes.
            See [1] for more details.

        EXAMPLES::

            sage: w = Word('abc')
            sage: it = w._iterated_right_palindromic_closure_recursive_iterator()
            sage: Word(it)
            word: abacaba

        ::

            sage: w = Word('aaa')
            sage: it = w._iterated_right_palindromic_closure_recursive_iterator()
            sage: Word(it)
            word: aaa

        ::

            sage: w = Word('abbab')
            sage: it = w._iterated_right_palindromic_closure_recursive_iterator()
            sage: Word(it)
            word: ababaabababaababa

        An infinite word::

            sage: t = words.ThueMorseWord('ab')
            sage: it = t._iterated_right_palindromic_closure_recursive_iterator()
            sage: Word(it)
            word: ababaabababaababaabababaababaabababaabab...

        TESTS:

        The empty word::

            sage: w = Word()
            sage: it = w._iterated_right_palindromic_closure_recursive_iterator()
            sage: next(it)
            Traceback (most recent call last):
            ...
            StopIteration

        REFERENCES:

        -   [1] A. de Luca, A. De Luca, Pseudopalindrome closure operators
            in free monoids, Theoret. Comput. Sci. 362 (2006) 282--300.
        -   [2] J. Justin, Episturmian morphisms and a Galois theorem on
            continued fractions, RAIRO Theoret. Informatics Appl. 39 (2005)
            207-215.
        """
        parent = self.parent().factors()
        ipcw = self[:0]
        lengths = []
        for i, letter in enumerate(self):
            lengths.append(ipcw.length())
            w = self[:i]
            pos = w.rfind(parent([letter]))
            if pos == -1:
                to_append = parent([letter]).palindromic_closure(f=f) + ipcw
            else:
                to_append = ipcw[lengths[pos]:]
            ipcw += to_append
            yield from to_append

    def iterated_right_palindromic_closure(self, f=None, algorithm='recursive'):
        r"""
        Return the iterated (`f`-)palindromic closure of ``self``.

        INPUT:

        - ``f`` -- involution (default: ``None``) on the alphabet of ``self``;
          it must be callable on letters as well as words (e.g. WordMorphism)

        - ``algorithm`` -- string (default: ``'recursive'``); specifying which
          algorithm to be used when computing the iterated palindromic closure.
          It must be one of the two following values:

          - ``'definition'`` -- computed using the definition
          - ``'recursive'`` -- computation based on an efficient formula
            that recursively computes the iterated right palindromic closure
            without having to recompute the longest `f`-palindromic suffix
            at each iteration [2].

        OUTPUT: word; the iterated (`f`-)palindromic closure of ``self``

        EXAMPLES::

            sage: Word('123').iterated_right_palindromic_closure()
            word: 1213121

        ::

            sage: w = Word('abc')
            sage: w.iterated_right_palindromic_closure()
            word: abacaba

        ::

            sage: w = Word('aaa')
            sage: w.iterated_right_palindromic_closure()
            word: aaa

        ::

            sage: w = Word('abbab')
            sage: w.iterated_right_palindromic_closure()
            word: ababaabababaababa

        A right `f`-palindromic closure::

            sage: f = WordMorphism('a->b,b->a')
            sage: w = Word('abbab')
            sage: w.iterated_right_palindromic_closure(f=f)
            word: abbaabbaababbaabbaabbaababbaabbaab

        An infinite word::

            sage: t = words.ThueMorseWord('ab')
            sage: t.iterated_right_palindromic_closure()
            word: ababaabababaababaabababaababaabababaabab...

        There are two implementations computing the iterated right
        `f`-palindromic closure, the latter being much more efficient::

            sage: w = Word('abaab')
            sage: u = w.iterated_right_palindromic_closure(algorithm='definition')
            sage: v = w.iterated_right_palindromic_closure(algorithm='recursive')
            sage: u
            word: abaabaababaabaaba
            sage: u == v
            True
            sage: w = words.RandomWord(8)
            sage: u = w.iterated_right_palindromic_closure(algorithm='definition')
            sage: v = w.iterated_right_palindromic_closure(algorithm='recursive')
            sage: u == v
            True

        TESTS:

        The empty word::

            sage: w = Word()
            sage: w.iterated_right_palindromic_closure()
            word:

        The length-`1` word::

            sage: Word('1').iterated_right_palindromic_closure()
            word: 1

        If the word is finite, so is the result::

            sage: w = Word([0,1]*7)
            sage: c = w.iterated_right_palindromic_closure()
            sage: type(c)
            <class 'sage.combinat.words.word.FiniteWord_iter_with_caching'>

        REFERENCES:

        -   [1] A. de Luca, A. De Luca, Pseudopalindrome closure operators
            in free monoids, Theoret. Comput. Sci. 362 (2006) 282--300.
        -   [2] J. Justin, Episturmian morphisms and a Galois theorem on
            continued fractions, RAIRO Theoret. Informatics Appl. 39 (2005)
            207-215.
        """
        if algorithm == 'definition':
            it = self._iterated_right_palindromic_closure_iterator(f=f)
        elif algorithm == 'recursive':
            it = self._iterated_right_palindromic_closure_recursive_iterator(f=f)
        else:
            raise ValueError("algorithm (=%s) must be either 'definition' or 'recursive'")

        if self.is_finite():
            return self._parent(it)
        else:
            from sage.combinat.words.words import Words
            parent = Words(self._parent.alphabet())
            return parent(it)

    def prefixes_iterator(self, max_length=None):
        r"""
        Return an iterator over the prefixes of ``self``.

        INPUT:

        - ``max_length`` -- nonnegative integer or ``None`` (default); the
          maximum length of the prefixes

        OUTPUT: iterator

        EXAMPLES::

            sage: w = Word('abaaba')
            sage: for p in w.prefixes_iterator(): p
            word:
            word: a
            word: ab
            word: aba
            word: abaa
            word: abaab
            word: abaaba
            sage: for p in w.prefixes_iterator(max_length=3): p
            word:
            word: a
            word: ab
            word: aba

        You can iterate over the prefixes of an infinite word::

            sage: f = words.FibonacciWord()
            sage: for p in f.prefixes_iterator(max_length=8): p
            word:
            word: 0
            word: 01
            word: 010
            word: 0100
            word: 01001
            word: 010010
            word: 0100101
            word: 01001010

        TESTS::

            sage: list(f.prefixes_iterator(max_length=0))
            [word: ]
        """
        to_consider = self if max_length is None else self[:max_length]
        yield self[:0]
        for (i, a) in enumerate(to_consider):
            yield self[:i + 1]

    def palindrome_prefixes_iterator(self, max_length=None):
        r"""
        Return an iterator over the palindrome prefixes of ``self``.

        INPUT:

        - ``max_length`` -- nonnegative integer or ``None`` (default); the
          maximum length of the prefixes

        OUTPUT: iterator

        EXAMPLES::

            sage: w = Word('abaaba')
            sage: for pp in w.palindrome_prefixes_iterator(): pp
            word:
            word: a
            word: aba
            word: abaaba
            sage: for pp in w.palindrome_prefixes_iterator(max_length=4): pp
            word:
            word: a
            word: aba

        You can iterate over the palindrome prefixes of an infinite word::

            sage: f = words.FibonacciWord()
            sage: for pp in f.palindrome_prefixes_iterator(max_length=20): pp
            word:
            word: 0
            word: 010
            word: 010010
            word: 01001010010
            word: 0100101001001010010
        """
        for p in self.prefixes_iterator(max_length):
            if p.is_palindrome():
                yield p

    def _partial_sums_iterator(self, start, mod=None):
        r"""
        Iterator over the partial sums of the prefixes of ``self``.

        INPUT:

        - ``self`` -- a word over the integers
        - ``start`` -- integer; the first letter of the resulting word
        - ``mod`` -- (default: ``None``) it can be one of the following:
            - None or 0 : result is over the integers
            - integer : result is over the integers modulo ``mod``.

        EXAMPLES::

            sage: w = Word(range(8))
            sage: list(w._partial_sums_iterator(0, mod=10))
            [0, 0, 1, 3, 6, 0, 5, 1, 8]

        ::

            sage: w = Word([1,1,1,1,1,1,1,1,1,1,1,1])
            sage: list(w._partial_sums_iterator(0, mod=10))
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2]

        ::

            sage: w = Word([1,1,1,1,1,1,1,1,1,1,1,1])
            sage: list(w._partial_sums_iterator(0))
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
        """
        if mod in (None, 0):
            sum = start

        elif mod in ZZ:
            Zn = Integers(mod)
            sum = Zn(start)

        else:
            raise TypeError('mod(=%s) must be None or an integer' % mod)

        yield sum
        for letter in self:
            sum += letter
            yield sum

    def partial_sums(self, start, mod=None):
        r"""
        Return the word defined by the partial sums of its prefixes.

        INPUT:

        - ``self`` -- a word over the integers
        - ``start`` -- integer; the first letter of the resulting word
        - ``mod`` -- (default: ``None``) it can be one of the following:
            - None or 0 : result is over the integers
            - integer : result is over the integers modulo ``mod``.

        EXAMPLES::

            sage: w = Word(range(10))
            sage: w.partial_sums(0)
            word: 0,0,1,3,6,10,15,21,28,36,45
            sage: w.partial_sums(1)
            word: 1,1,2,4,7,11,16,22,29,37,46

        ::

            sage: w = Word([1,2,3,1,2,3,2,2,2,2])
            sage: w.partial_sums(0, mod=None)
            word: 0,1,3,6,7,9,12,14,16,18,20
            sage: w.partial_sums(0, mod=0)
            word: 0,1,3,6,7,9,12,14,16,18,20
            sage: w.partial_sums(0, mod=8)
            word: 01367146024
            sage: w.partial_sums(0, mod=4)
            word: 01323102020
            sage: w.partial_sums(0, mod=2)
            word: 01101100000
            sage: w.partial_sums(0, mod=1)
            word: 00000000000

        TESTS:

        If the word is infinite, so is the result::

            sage: w = Word(lambda n:1)
            sage: u = w.partial_sums(0)
            sage: type(u)
            <class 'sage.combinat.words.word.InfiniteWord_iter_with_caching'>
        """
        it = self._partial_sums_iterator(start=start, mod=mod)

        if mod in (None, 0):
            alphabet = None
        elif mod in ZZ:
            alphabet = Integers(mod)

        if self.is_finite():
            length = "finite"
        elif self.length() == Infinity:
            length = "infinite"
        else:
            length = "unknown"
        from sage.combinat.words.word import Word
        return Word(it, alphabet=alphabet, length=length)

    def _finite_differences_iterator(self, mod=None):
        r"""
        Iterator over the differences of consecutive letters of ``self``.

        INPUT:

        - ``self`` -- a word over the integers
        - ``mod`` -- (default: ``None``) it can be one of the following:
            - None or 0 : result is over the integers
            - integer : result is over the integers modulo ``mod``.

        EXAMPLES::

            sage: w = Word(x^2 for x in range(10))
            sage: list(w._finite_differences_iterator())
            [1, 3, 5, 7, 9, 11, 13, 15, 17]

        ::

            sage: w = Word([1,6,8,4,2,6,8,2,3])
            sage: list(w._finite_differences_iterator())
            [5, 2, -4, -2, 4, 2, -6, 1]
            sage: list(w._finite_differences_iterator(4))
            [1, 2, 0, 2, 0, 2, 2, 1]
            sage: list(w._finite_differences_iterator(5))
            [0, 2, 1, 3, 4, 2, 4, 1]

        TESTS::

            sage: w = Word([2,3,6])
            sage: list(w._finite_differences_iterator())
            [1, 3]
            sage: w = Word([2,6])
            sage: list(w._finite_differences_iterator())
            [4]
            sage: w = Word([2])
            sage: list(w._finite_differences_iterator())
            []
            sage: w = Word()
            sage: list(w._finite_differences_iterator())
            []

        ::

            sage: list(w._finite_differences_iterator('a'))
            Traceback (most recent call last):
            ...
            TypeError: mod(=a) must be None or an integer
        """
        if mod in (None, 0):
            i = iter(self)
            j = iter(self)

            try:
                next(j)
                while True:
                    yield next(j) - next(i)
            except StopIteration:
                return

        elif mod in ZZ:
            Zn = Integers(mod)
            i = iter(self)
            j = iter(self)

            try:
                next(j)
                while True:
                    yield Zn(next(j) - next(i))
            except StopIteration:
                return

        else:
            raise TypeError('mod(=%s) must be None or an integer' % mod)

    def finite_differences(self, mod=None):
        r"""
        Return the word obtained by the differences of consecutive letters
        of ``self``.

        INPUT:

        - ``self`` -- a word over the integers
        - ``mod`` -- (default: ``None``) it can be one of the following:
            - None or 0 : result is over the integers
            - integer : result is over the integers modulo ``mod``.

        EXAMPLES::

            sage: w = Word([x^2 for x in range(10)])
            sage: w.finite_differences()
            word: 1,3,5,7,9,11,13,15,17
            sage: w.finite_differences(mod=4)
            word: 131313131
            sage: w.finite_differences(mod=0)
            word: 1,3,5,7,9,11,13,15,17

        TESTS::

            sage: w = Word([2,3,6])
            sage: w.finite_differences()
            word: 13
            sage: w = Word([2,6])
            sage: w.finite_differences()
            word: 4
            sage: w = Word([2])
            sage: w.finite_differences()
            word:
            sage: w = Word()
            sage: w.finite_differences()
            word:

        If the word is infinite, so is the result::

            sage: w = Word(lambda n:n)
            sage: u = w.finite_differences()
            sage: u
            word: 1111111111111111111111111111111111111111...
            sage: type(u)
            <class 'sage.combinat.words.word.InfiniteWord_iter_with_caching'>
        """
        it = self._finite_differences_iterator(mod=mod)

        if mod in (None, 0):
            alphabet = None
        elif mod in ZZ:
            alphabet = Integers(mod)

        if self.is_finite():
            length = "finite"
        elif self.length() == Infinity:
            length = "infinite"
        else:
            length = "unknown"
        from sage.combinat.words.word import Word
        return Word(it, alphabet=alphabet, length=length)

    def sum_digits(self, base=2, mod=None):
        r"""
        Return the sequence of the sum modulo ``mod`` of the digits written
        in base ``base`` of ``self``.

        INPUT:

        - ``self`` -- word over natural numbers

        - ``base`` -- integer (default: 2) greater or equal to 2

        - ``mod`` -- modulo (default: ``None``); can take the following
          values:
          - ``integer`` -- the modulo
          - ``None`` -- the value ``base`` is considered for the modulo

        EXAMPLES:

        The Thue-Morse word::

            sage: from itertools import count
            sage: Word(count()).sum_digits()
            word: 0110100110010110100101100110100110010110...

        Sum of digits modulo 2 of the prime numbers written in base 2::

            sage: Word(primes(1000)).sum_digits()                                       # needs sage.libs.pari
            word: 1001110100111010111011001011101110011011...

        Sum of digits modulo 3 of the prime numbers written in base 3::

            sage: Word(primes(1000)).sum_digits(base=3)                                 # needs sage.libs.pari
            word: 2100002020002221222121022221022122111022...
            sage: Word(primes(1000)).sum_digits(base=3, mod=3)                          # needs sage.libs.pari
            word: 2100002020002221222121022221022122111022...

        Sum of digits modulo 2 of the prime numbers written in base 3::

            sage: Word(primes(1000)).sum_digits(base=3, mod=2)                          # needs sage.libs.pari
            word: 0111111111111111111111111111111111111111...

        Sum of digits modulo 7 of the prime numbers written in base 10::

            sage: Word(primes(1000)).sum_digits(base=10, mod=7)                         # needs sage.libs.pari
            word: 2350241354435041006132432241353546006304...

        Negative entries::

            sage: w = Word([-1,0,1,2,3,4,5])
            sage: w.sum_digits()
            Traceback (most recent call last):
            ...
            NotImplementedError: nth digit of Thue-Morse word is not implemented for negative value of n

        TESTS:

        The Thue-Morse word::

            sage: from itertools import count
            sage: w = Word(count()).sum_digits()
            sage: t = words.ThueMorseWord()
            sage: w[:100]  == t[:100]
            True

        ::

            sage: type(Word(range(10)).sum_digits())
            <class 'sage.combinat.words.word.FiniteWord_iter_with_caching'>
        """
        from functools import partial
        from sage.combinat.words.word_generators import words

        # The alphabet
        if mod is None and base >= 2:
            alphabet = list(range(base))
        elif mod in ZZ and mod >= 2:
            alphabet = list(range(mod))
        else:
            raise ValueError("base (=%s) and mod (=%s) must be integers greater or equal to 2" % (base, mod))

        # The iterator
        f = partial(words._ThueMorseWord_nth_digit, alphabet=alphabet, base=base)
        it = (f(a) for a in self)

        # The length
        if self.is_finite():
            length = "finite"
        elif self.length() == Infinity:
            length = None
        else:
            length = "unknown"

        from sage.combinat.words.word import Word
        return Word(it, alphabet=alphabet, length=length, datatype='iter')

    def first_occurrence(self, other, start=0):
        r"""
        Return the position of the first occurrence of ``other`` in ``self``.

        If ``other`` is not a factor of ``self``, it returns ``None``
        or loops forever when ``self`` is an infinite word.

        INPUT:

        - ``other`` -- a finite word
        - ``start`` -- integer (default: `0`) where the search starts

        OUTPUT: integer or ``None``

        EXAMPLES::

            sage: w = Word('01234567890123456789')
            sage: w.first_occurrence(Word('3456'))
            3
            sage: w.first_occurrence(Word('3456'), start=7)
            13

        When the factor is not present, ``None`` is returned::

            sage: w.first_occurrence(Word('3456'), start=17) is None
            True
            sage: w.first_occurrence(Word('3333')) is None
            True

        Also works for searching a finite word in an infinite word::

            sage: w = Word('0123456789')^oo
            sage: w.first_occurrence(Word('3456'))
            3
            sage: w.first_occurrence(Word('3456'), start=1000)
            1003

        But it will loop for ever if the factor is not found::

            sage: w.first_occurrence(Word('3333')) # not tested -- infinite loop

        The empty word occurs in a word::

            sage: Word('123').first_occurrence(Word(''), 0)
            0
            sage: Word('').first_occurrence(Word(''), 0)
            0
        """
        lf = other.length()
        lm = self.length()
        if lf == 0:
            return start
        elif lm == 0:
            return None
        occ = other.last_position_dict()
        suff = other.good_suffix_table()
        s = start
        while s <= lm - lf:
            for j in range(lf-1, -1, -1):
                a = self[s+j]
                if other[j] != a:
                    s += max(suff[j + 1], j - occ.get(a, -1))
                    break
            else:
                return s
        return None

    def factor_occurrences_iterator(self, fact):
        r"""
        Return an iterator over all occurrences (including overlapping ones)
        of fact in ``self`` in their order of appearance.

        INPUT:

        - ``fact`` -- a non empty finite word

        OUTPUT: iterator

        EXAMPLES::

            sage: TM = words.ThueMorseWord()
            sage: fact = Word([0,1,1,0,1])
            sage: it = TM.factor_occurrences_iterator(fact)
            sage: next(it)
            0
            sage: next(it)
            12
            sage: next(it)
            24

        ::

            sage: u = Word('121')
            sage: w = Word('121213211213')
            sage: list(w.factor_occurrences_iterator(u))
            [0, 2, 8]
        """
        if fact.is_empty():
            raise NotImplementedError("The factor must be non empty")
        if not fact.is_finite():
            raise ValueError("The factor must be finite")
        p = self.first_occurrence(fact, start=0)
        while p is not None:
            yield p
            p = self.first_occurrence(fact, start=p+1)

    def return_words_iterator(self, fact):
        r"""
        Return an iterator over all the return words of fact in self
        (without unicity).

        INPUT:

        - ``fact`` -- a non empty finite word

        OUTPUT: iterator

        EXAMPLES::

            sage: w = Word('baccabccbacbca')
            sage: b = Word('b')
            sage: list(w.return_words_iterator(b))
            [word: bacca, word: bcc, word: bac]

        ::

            sage: TM = words.ThueMorseWord()
            sage: fact = Word([0,1,1,0,1])
            sage: it = TM.return_words_iterator(fact)
            sage: next(it)
            word: 011010011001
            sage: next(it)
            word: 011010010110
            sage: next(it)
            word: 0110100110010110
            sage: next(it)
            word: 01101001
            sage: next(it)
            word: 011010011001
            sage: next(it)
            word: 011010010110
        """
        it = self.factor_occurrences_iterator(fact)
        try:
            i = next(it)
            while True:
                j = next(it)
                yield self[i:j]
                i = j
        except StopIteration:
            return

    def complete_return_words_iterator(self, fact):
        r"""
        Return an iterator over all the complete return words of fact in
        ``self`` (without unicity).

        A complete return words `u` of a factor `v`  is a factor starting
        by the given factor `v` and ending just after the next occurrence
        of this factor `v`. See for instance [1].

        INPUT:

        - ``fact`` -- a non empty finite word

        OUTPUT: iterator

        EXAMPLES::

            sage: TM = words.ThueMorseWord()
            sage: fact = Word([0,1,1,0,1])
            sage: it = TM.complete_return_words_iterator(fact)
            sage: next(it)
            word: 01101001100101101
            sage: next(it)
            word: 01101001011001101
            sage: next(it)
            word: 011010011001011001101
            sage: next(it)
            word: 0110100101101
            sage: next(it)
            word: 01101001100101101
            sage: next(it)
            word: 01101001011001101

        REFERENCES:

        -   [1] J. Justin, L. Vuillon, Return words in Sturmian and
            episturmian words, Theor. Inform. Appl. 34 (2000) 343--356.
        """
        it = self.factor_occurrences_iterator(fact)
        L = fact.length()
        try:
            i = next(it)
            while True:
                j = next(it)
                yield self[i:j+L]
                i = j
        except StopIteration:
            return
