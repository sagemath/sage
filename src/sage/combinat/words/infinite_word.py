# sage_setup: distribution = sagemath-combinat
r"""
Infinite word

AUTHORS:

- Sebastien Labbe
- Franco Saliola

EXAMPLES:

============================
Creation of an infinite word
============================

Periodic infinite words::

    sage: v = Word([0, 4, 8, 8, 3])
    sage: vv = v^Infinity
    sage: vv
    word: 0488304883048830488304883048830488304883...

Infinite words from a function `f:\mathbb{N}\rightarrow A`
over an alphabet `A`::

    sage: Word(lambda n: n%3)
    word: 0120120120120120120120120120120120120120...

::

    sage: def t(n):
    ....:     return add(Integer(n).digits(base=2)) % 2
    sage: Word(t, alphabet = [0, 1])
    word: 0110100110010110100101100110100110010110...

or as a one-liner::

    sage: Word(lambda n : add(Integer(n).digits(base=2)) % 2, alphabet = [0, 1])
    word: 0110100110010110100101100110100110010110...

Infinite words from iterators::

    sage: from itertools import count,repeat
    sage: Word( repeat(4) )
    word: 4444444444444444444444444444444444444444...
    sage: Word( count() )
    word: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,...

Infinite words from morphism

For example, let `A=\{a,b\}` and `\mu : A^* \rightarrow A^*`
be the morphism defined by `a\mapsto ab, b\mapsto ba`::

    sage: mu = WordMorphism('a->ab,b->ba'); mu
    WordMorphism: a->ab, b->ba
    sage: mu.fixed_point('a')
    word: abbabaabbaababbabaababbaabbabaabbaababba...

Infinite words in a specific combinatorial class::

    sage: W = InfiniteWords("ab"); W
    Infinite words over {'a', 'b'}
    sage: f = lambda n : 'a' if n % 2 == 1 else 'b'
    sage: W(f)
    word: babababababababababababababababababababa...
"""
# ****************************************************************************
#       Copyright (C) 2008 Sebastien Labbe <slabqc@gmail.com>,
#                          Franco Saliola <saliola@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.combinat.words.abstract_word import Word_class
from sage.combinat.words.word_options import word_options
from sage.rings.infinity import Infinity


class InfiniteWord_class(Word_class):
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
        global word_options
        if word_options['old_repr']:
            return "Infinite word over %s" % str(self.parent().alphabet())[17:]
        return word_options['identifier'] + self.string_rep()

    def length(self):
        r"""
        Return the length of ``self``.

        EXAMPLES::

            sage: f = lambda n : n % 6
            sage: w = Word(f); w
            word: 0123450123450123450123450123450123450123...
            sage: w.length()
            +Infinity
        """
        return Infinity
