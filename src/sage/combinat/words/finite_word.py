r"""
Finite word

AUTHORS:

- Arnaud Bergeron
- Amy Glen
- Sébastien Labbé
- Franco Saliola
- Julien Leroy (March 2010): reduced_rauzy_graph

EXAMPLES:

=========================
Creation of a finite word
=========================

Finite words from Python strings, lists and tuples::

    sage: Word("abbabaab")
    word: abbabaab
    sage: Word([0, 1, 1, 0, 1, 0, 0, 1])
    word: 01101001
    sage: Word( ('a', 0, 5, 7, 'b', 9, 8) )
    word: a057b98

Finite words from functions::

    sage: f = lambda n : n%3
    sage: Word(f, length=13)
    word: 0120120120120

Finite words from iterators::

    sage: from itertools import count
    sage: Word(count(), length=10)
    word: 0123456789

::

    sage: Word( iter('abbccdef') )
    word: abbccdef

Finite words from words via concatenation::

    sage: u = Word("abcccabba")
    sage: v = Word([0, 4, 8, 8, 3])
    sage: u * v
    word: abcccabba04883
    sage: v * u
    word: 04883abcccabba
    sage: u + v
    word: abcccabba04883
    sage: u^3 * v^(8/5)
    word: abcccabbaabcccabbaabcccabba04883048

Finite words from infinite words::

    sage: vv = v^Infinity
    sage: vv[10000:10015]
    word: 048830488304883

Finite words in a specific combinatorial class::

    sage: W = Words("ab")
    sage: W
    Finite and infinite words over {'a', 'b'}
    sage: W("abbabaab")
    word: abbabaab
    sage: W(["a","b","b","a","b","a","a","b"])
    word: abbabaab
    sage: W( iter('ababab') )
    word: ababab

Finite word as the image under a morphism::

    sage: m = WordMorphism({0:[4,4,5,0],5:[0,5,5],4:[4,0,0,0]})
    sage: m(0)
    word: 4450
    sage: m(0, order=2)
    word: 400040000554450
    sage: m(0, order=3)
    word: 4000445044504450400044504450445044500550...

.. NOTE::

    The following two finite words have the same string representation::

        sage: w = Word('010120')
        sage: z = Word([0, 1, 0, 1, 2, 0])
        sage: w
        word: 010120
        sage: z
        word: 010120

    but are not equal::

        sage: w == z
        False

    Indeed, w and z are defined on different alphabets::

        sage: w[2]
        '0'
        sage: z[2]
        0

========================
Functions and algorithms
========================

There are more than 100 functions defined on a finite word. Here are some
of them::

    sage: w = Word('abaabbba'); w
    word: abaabbba
    sage: w.is_palindrome()
    False
    sage: w.is_lyndon()
    False
    sage: w.number_of_factors()
    28
    sage: w.critical_exponent()
    3

::

    sage: print(w.lyndon_factorization())
    (ab, aabbb, a)
    sage: print(w.crochemore_factorization())
    (a, b, a, ab, bb, a)

::

    sage: st = w.suffix_tree()
    sage: st
    Implicit Suffix Tree of the word: abaabbba
    sage: st.show(word_labels=True)                                                     # needs sage.plot

::

    sage: T = words.FibonacciWord('ab')
    sage: T.longest_common_prefix(Word('abaabababbbbbb'))
    word: abaababa

As matrix and many other sage objects, words have a parent::

    sage: u = Word('xyxxyxyyy')
    sage: u.parent()
    Finite words over Set of Python objects of class 'object'

::

    sage: v = Word('xyxxyxyyy', alphabet='xy')
    sage: v.parent()
    Finite words over {'x', 'y'}

========================
Factors and Rauzy Graphs
========================

Enumeration of factors, the successive values returned by ``next(it)``
can appear in a different order depending on hardware. Therefore we
mark the three first results of the test ``random``. The important test
is that the iteration stops properly on the fourth call::

    sage: w = Word([4,5,6])^7
    sage: it = w.factor_iterator(4)
    sage: next(it) # random
    word: 6456
    sage: next(it) # random
    word: 5645
    sage: next(it) # random
    word: 4564
    sage: next(it)
    Traceback (most recent call last):
    ...
    StopIteration

The set of factors::

    sage: sorted(w.factor_set(3))
    [word: 456, word: 564, word: 645]
    sage: sorted(w.factor_set(4))
    [word: 4564, word: 5645, word: 6456]
    sage: w.factor_set().cardinality()
    61

Rauzy graphs::

    sage: f = words.FibonacciWord()[:30]
    sage: f.rauzy_graph(4)                                                              # needs sage.graphs
    Looped digraph on 5 vertices
    sage: f.reduced_rauzy_graph(4)                                                      # needs sage.graphs
    Looped multi-digraph on 2 vertices

Left-special and bispecial factors::

    sage: f.number_of_left_special_factors(7)
    1
    sage: f.bispecial_factors()
    [word: , word: 0, word: 010, word: 010010, word: 01001010010]
"""
# ****************************************************************************
#       Copyright (C) 2008 Arnaud Bergeron <abergeron@gmail.com>,
#                     2008 Amy Glen <amy.glen@gmail.com>,
#                     2008-2012 Sébastien Labbé <slabqc@gmail.com>,
#                     2008-2010 Franco Saliola <saliola@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from itertools import repeat
from collections import defaultdict
from itertools import islice, cycle

from sage.combinat.words.abstract_word import Word_class
from sage.combinat.words.words import Words
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import
from sage.combinat.words.word_options import word_options
from sage.rings.infinity import Infinity
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.sets.set import Set

lazy_import('sage.groups.perm_gps.permgroup_element', 'PermutationGroupElement')


class FiniteWord_class(Word_class):
    def __str__(self):
        r"""
        Return the full (not truncated) string representation of the word
        without identifier.

        TESTS::

            sage: Word('abc').__str__()
            'abc'
            sage: Word([0, 1, 0, 0, 1] * 10).__str__()
            '01001010010100101001010010100101001010010100101001'
            sage: Word([0,1,10,101]).__str__()
            '0,1,10,101'

        Insertion into a ``str``::

            sage: w = Word(range(5))
            sage: "Let's insert the word w = %s in this string." % w
            "Let's insert the word w = 01234 in this string."

        Using ``LatexExpr``::

            sage: from sage.misc.latex import LatexExpr
            sage: LatexExpr(w)
            01234

        With the ``print`` statement::

            sage: print(w)
            01234

        No truncation is done for finite words::

            sage: w = Word([i % 5 for i in range(60)])
            sage: print(w)
            012340123401234012340123401234012340123401234012340123401234
        """
        if word_options['display'] == 'string':
            ls = word_options['letter_separator']
            letters = [str(a) for a in self]
            if all(len(a) == 1 for a in letters):
                return ''.join(letters)
            else:
                return ls.join(letters)
        elif word_options['display'] == 'list':
            return str(list(self))

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: Word(range(10))._repr_()
            'word: 0123456789'
            sage: Word(range(100))._repr_()
            'word: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,...'
        """
        if word_options['old_repr']:
            if word_options['truncate'] and \
                    self.length() > word_options['truncate_length']:
                return "Finite word of length {} over {}".format(self.length(), str(self.parent().alphabet())[17:])
        return word_options['identifier'] + self.string_rep()

    def coerce(self, other):
        r"""
        Try to return a pair of words with a common parent; raise an
        exception if this is not possible.

        This function begins by checking if both words have the same
        parent. If this is the case, then no work is done and both words
        are returned as-is.

        Otherwise it will attempt to convert ``other`` to the domain of ``self``.
        If that fails, it will attempt to convert ``self`` to the domain of
        ``other``. If both attempts fail, it raises a :exc:`TypeError`
        to signal failure.

        EXAMPLES::

            sage: W1 = Words('abc'); W2 = Words('ab')
            sage: w1 = W1('abc'); w2 = W2('abba'); w3 = W1('baab')
            sage: w1.parent() is w2.parent()
            False
            sage: a, b = w1.coerce(w2)
            sage: a.parent() is b.parent()
            True
            sage: w1.parent() is w2.parent()
            False
        """
        if self.parent() != other.parent():
            try:
                other = self.parent()(other)
                other.parent()._check(other, length=None)
            except Exception:
                try:
                    self = other.parent()(self)
                    self.parent()._check(self, length=None)
                except Exception:
                    raise TypeError("no coercion rule between {!r} and {!r}".format(self.parent(), other.parent()))
        return self, other

    def __hash__(self):
        r"""
        Return the hash for this word.

        TESTS::

             sage: h = hash(Word('abc'))    # indirect test
             sage: Word('abc').__hash__() == Word('abc').__hash__()
             True
        """
        if self._hash is None:
            res = 5381
            for s in self._to_integer_iterator():
                res = ((res << 5) + res) + s
            self._hash = res
        return self._hash

    def concatenate(self, other):
        r"""
        Return the concatenation of ``self`` and ``other``.

        INPUT:

        - ``other`` -- a word over the same alphabet as ``self``

        EXAMPLES:

        Concatenation may be made using ``+`` or ``*`` operations::

            sage: w = Word('abadafd')
            sage: y = Word([5,3,5,8,7])
            sage: w * y
            word: abadafd53587
            sage: w + y
            word: abadafd53587
            sage: w.concatenate(y)
            word: abadafd53587

        Both words must be defined over the same alphabet::

            sage: z = Word('12223', alphabet = '123')
            sage: z + y
            Traceback (most recent call last):
            ...
            ValueError: 5 not in alphabet

        Eventually, it should work::

            sage: z = Word('12223', alphabet = '123')
            sage: z + y                   #todo: not implemented
            word: 1222353587

        TESTS:

        The empty word is not considered by concatenation::

            sage: type(Word([]) * Word('abcd'))
            <class 'sage.combinat.words.word.FiniteWord_str'>
            sage: type(Word('abcd') * Word())
            <class 'sage.combinat.words.word.FiniteWord_str'>
            sage: type(Word('abcd') * Word([]))
            <class 'sage.combinat.words.word.FiniteWord_str'>
            sage: type(Word('abcd') * Word(()))
            <class 'sage.combinat.words.word.FiniteWord_str'>
            sage: type(Word([1,2,3]) * Word(''))
            <class 'sage.combinat.words.word.FiniteWord_list'>

        Concatenation of finite words with infinite words works as expected::

            sage: from itertools import repeat
            sage: W = Words('ab')
            sage: w1 = W('aba')
            sage: w2 = W(repeat('b'), length='infinite')
            sage: w1*w2
            word: ababbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb...
            sage: _.parent()
            Infinite words over {'a', 'b'}
        """
        if self.is_empty():
            return other
        if isinstance(other, Word_class) and other.is_empty():
            return self
        f = CallableFromListOfWords([self, other])
        length = self.length() + other.length()
        parent = self._parent
        if length == Infinity:
            parent = parent.shift()
            return parent(f, datatype='callable', caching=True)
        else:
            return parent(f, length=length, datatype='callable', caching=True)

    __mul__ = concatenate

    __add__ = concatenate

    # TODO: This function is using domain=range(n) for Word but
    # should be a domain=slice(n) # Seb : Feb 23th : I think this is fine now!!
    def __pow__(self, exp):
        r"""
        Return the ``exp``-th power of ``self``.

        If ``exp`` is `\infty`, returns the infinite periodic word of base ``self``.
        Otherwise, `|w|\cdot exp` must be a nonnegative integer.

        INPUT:

        - ``exp`` -- integer; a rational, a float number or plus infinity

        OUTPUT: word; the ``exp``-th power of ``self``

        EXAMPLES:

        You can take nonnegative integer powers::

            sage: w = Word(range(6)); w
            word: 012345
            sage: w^2
            word: 012345012345
            sage: w^1
            word: 012345
            sage: w^0
            word:
            sage: w^(-1)
            Traceback (most recent call last):
            ...
            ValueError: Power of the word is not defined on the exponent -1: the length of the word (6) times the exponent (-1) must be a positive integer


        You can take nonnegative rational powers::

            sage: w = Word(range(6)); w
            word: 012345
            sage: w^(1/2)
            word: 012
            sage: w^(1/3)
            word: 01
            sage: (w*w)^(1/2) == w
            True
            sage: w^(5/2)
            word: 012345012345012

        ...but the length of the word times the exponent must be an integer::

            sage: w = Word(range(6))
            sage: w^(1/4)
            Traceback (most recent call last):
            ...
            ValueError: Power of the word is not defined on the exponent 1/4: the length of the word (6) times the exponent (1/4) must be a positive integer

        You can take infinite power::

            sage: w = Word(range(6)); w
            word: 012345
            sage: u = w^oo; u
            word: 0123450123450123450123450123450123450123...
            sage: u[10000000:20000000]
            word: 4501234501234501234501234501234501234501...
            sage: u[10000000:10000020]
            word: 45012345012345012345
            sage: Word()^oo
            word:
        """
        # powers of the empty word
        if self.is_empty():
            return self

        # infinite power of a non-empty word
        def fcn(n):
            return self[n % self.length()]

        if exp is Infinity:
            return self._parent.shift()(fcn)

        # If exp*|self| is not an integer
        length = exp * self.length()
        if length in ZZ and length >= 0:
            return self._parent(fcn, length=length)
        else:
            raise ValueError("Power of the word is not defined on the exponent {}:"
                    " the length of the word ({}) times the exponent ({}) must"
                    " be a positive integer".format(exp, self.length(), exp))

    def length(self):
        r"""
        Return the length of ``self``.

        TESTS::

            sage: from sage.combinat.words.word import Word_class
            sage: w = Word(iter('abba'*40), length='finite')
            sage: w._len is None
            True
            sage: w.length()
            160
            sage: w = Word(iter('abba'), length=4)
            sage: w._len
            4
            sage: w.length()
            4
            sage: def f(n):
            ....:   return list(range(2,12,2))[n]
            sage: w = Word(f, length=5)
            sage: w.length()
            5
        """
        if self._len is None:
            self._len = Integer(sum(1 for _ in self))
        return self._len

    def content(self, n=None):
        r"""
        Return content of ``self``.

        INPUT:

        - ``n`` -- (optional) an integer specifying the maximal
          letter in the alphabet

        OUTPUT: list where the `i`-th entry indicates the multiplicity
        of the `i`-th letter in the alphabet in ``self``

        EXAMPLES::

            sage: w = Word([1,2,4,3,2,2,2])
            sage: w.content()
            [1, 4, 1, 1]
            sage: w = Word([3,1])
            sage: w.content()
            [1, 1]
            sage: w.content(n=3)
            [1, 0, 1]
            sage: w = Word([2,4],alphabet=[1,2,3,4])
            sage: w.content(n=3)
            [0, 1, 0]
            sage: w.content()
            [0, 1, 0, 1]
        """
        from collections import Counter
        c = Counter(self)
        if n is not None:
            alphabet = range(1, n + 1)
        elif not self.parent().alphabet().cardinality() == +Infinity:
            alphabet = self.parent().alphabet()
        else:
            alphabet = sorted(c.keys())
        return [Integer(c[a]) for a in alphabet]

    def is_yamanouchi(self, n=None):
        r"""
        Return whether ``self`` is Yamanouchi.

        A word `w` is Yamanouchi if, when read from right to left, it
        always has weakly more `i`'s than `i+1`'s for all `i` that
        appear in `w`.

        INPUT:

        - ``n`` -- (optional) an integer specifying the maximal
          letter in the alphabet

        EXAMPLES::

            sage: w = Word([1,2,4,3,2,2,2])
            sage: w.is_yamanouchi()
            False
            sage: w = Word([2,3,4,3,1,2,1,1,2,1])
            sage: w.is_yamanouchi()
            True
            sage: w = Word([3,1])
            sage: w.is_yamanouchi(n=3)
            False
            sage: w.is_yamanouchi()
            True
            sage: w = Word([3,1],alphabet=[1,2,3])
            sage: w.is_yamanouchi()
            False
            sage: w = Word([2,1,1,2])
            sage: w.is_yamanouchi()
            False
        """
        from sage.combinat.words.word import Word
        if n is not None:
            w = Word(self, alphabet=list(range(1, n + 1)))
        elif not self.parent().alphabet().cardinality() == +Infinity:
            w = self
        else:
            w = Word(self, alphabet=sorted(self.letters()))
        l = w.length()
        for a in range(l - 1, -1, -1):
            mu = w.parent()(self[a:]).content()
            if not all(mu[i] >= mu[i+1] for i in range(len(mu)-1)):
                return False
        return True

    def schuetzenberger_involution(self, n=None):
        r"""
        Return the Schützenberger involution of the word ``self``, which is obtained
        by reverting the word and then complementing all letters within the
        underlying ordered alphabet. If ``n`` is specified, the underlying
        alphabet is assumed to be `[1,2,\ldots,n]`. If no alphabet is specified,
        `n` is the maximal letter appearing in ``self``.

        INPUT:

        - ``self`` -- a word
        - ``n`` -- integer specifying the maximal letter in the alphabet (optional)

        OUTPUT: a word, the Schützenberger involution of ``self``

        EXAMPLES::

            sage: w = Word([9,7,4,1,6,2,3])
            sage: v = w.schuetzenberger_involution(); v
            word: 7849631
            sage: v.parent()
            Finite words over Set of Python objects of class 'object'

            sage: w = Word([1,2,3],alphabet=[1,2,3,4,5])
            sage: v = w.schuetzenberger_involution();v
            word: 345
            sage: v.parent()
            Finite words over {1, 2, 3, 4, 5}

            sage: w = Word([1,2,3])
            sage: v = w.schuetzenberger_involution(n=5);v
            word: 345
            sage: v.parent()
            Finite words over Set of Python objects of class 'object'

            sage: w = Word([11,32,69,2,53,1,2,3,18,41])
            sage: w.schuetzenberger_involution()
            word: 29,52,67,68,69,17,68,1,38,59

            sage: w = Word([],alphabet=[1,2,3,4,5])
            sage: w.schuetzenberger_involution()
            word:

            sage: w = Word([])
            sage: w.schuetzenberger_involution()
            word:
        """
        if self.length() == 0:
            return self
        r = self.reversal()
        w = list(r)
        parent = self.parent()
        if n is None:
            alphsize = parent.alphabet().cardinality()
            if not alphsize == +Infinity:
                n = max(parent.alphabet())
            elif r.length() > 0:
                n = max(w)
        for k in range(r.length()):
            w[k] = n+1 - w[k]
        return parent(w, check=False)

    def foata_bijection(self):
        r"""
        Return word ``self`` under the Foata bijection.

        The Foata bijection `\phi` is a bijection on the set of words
        of given content (by a slight generalization of Section 2 in [FS1978]_).
        It can be defined by induction on the size of the word: Given a word
        `w_1 w_2 \cdots w_n`, start with `\phi(w_1) = w_1`. At the `i`-th step, if
        `\phi(w_1 w_2 \cdots w_i) = v_1 v_2 \cdots v_i`, we define
        `\phi(w_1 w_2 \cdots w_i w_{i+1})` by placing `w_{i+1}` on the end of
        the word `v_1 v_2 \cdots v_i` and breaking the word up into blocks
        as follows. If `w_{i+1} \ge v_i`, place a vertical line to the right
        of each `v_k` for which `w_{i+1} \ge v_k`. Otherwise, if
        `w_{i+1} < v_i`, place a vertical line to the right of each `v_k`
        for which `w_{i+1} < v_k`. In either case, place a vertical line at
        the start of the word as well. Now, within each block between
        vertical lines, cyclically shift the entries one place to the
        right.

        For instance, to compute `\phi([4,1,5,4,2,2,3])`, the sequence of
        words is

        * `4`,
        * `|4|1 \to 41`,
        * `|4|1|5 \to 415`,
        * `|415|4 \to 5414`,
        * `|5|4|14|2 \to 54412`,
        * `|5441|2|2 \to 154422`,
        * `|1|5442|2|3 \to 1254423`.

        So `\phi([4,1,5,4,2,2,3]) = [1,2,5,4,4,2,3]`.

        .. SEEALSO::

            :meth:`Foata bijection on Permutations <sage.combinat.permutation.Permutation.foata_bijection()>`.

        EXAMPLES::

            sage: w = Word([2,2,2,1,1,1])
            sage: w.foata_bijection()
            word: 112221
            sage: w = Word([2,2,1,2,2,2,1,1,2,1])
            sage: w.foata_bijection()
            word: 2122212211
            sage: w = Word([4,1,5,4,2,2,3])
            sage: w.foata_bijection()
            word: 1254423

        TESTS::

            sage: w = Word('121314')
            sage: w.foata_bijection()
            word: 231114
            sage: w = Word('1133a1')
            sage: w.foata_bijection()
            word: 3113a1
        """
        s = self.standard_permutation()
        ordered_alphabet = sorted(self.letters(),
                                  key=self.parent().sortkey_letters)
        eval_dict = self.evaluation_dict()
        weight = [eval_dict[a] for a in ordered_alphabet]
        return (s.foata_bijection()).destandardize(weight, ordered_alphabet=ordered_alphabet)

    def major_index(self, final_descent=False):
        r"""
        Return the major index of ``self``.

        The major index of a word `w` is the sum of the descents of `w`.

        With the ``final_descent`` option, the last position of a
        non-empty word is also considered as a descent.

        .. SEEALSO::

            :meth:`major index on Permutations <sage.combinat.permutation.Permutation.major_index()>`.

        EXAMPLES::

            sage: w = Word([2,1,3,3,2])
            sage: w.major_index()
            5
            sage: w = Word([2,1,3,3,2])
            sage: w.major_index(final_descent=True)
            10
        """
        return (self.standard_permutation()).major_index(final_descent=final_descent)

    def number_of_inversions(self):
        r"""
        Return the number of inversions in ``self``.

        An inversion of a word `w = w_1 \ldots w_n` is a pair of indices `(i, j)`
        with `i < j` and `w_i > w_j`.

        .. SEEALSO::

            :meth:`number of inversions on Permutations <sage.combinat.permutation.Permutation.number_of_inversions()>`.

        EXAMPLES::

            sage: w = Word([2,1,3,3,2])
            sage: w.number_of_inversions()
            3
        """
        return (self.standard_permutation()).number_of_inversions()

    def is_empty(self):
        r"""
        Return ``True`` if the length of ``self`` is zero,
        and ``False`` otherwise.

        EXAMPLES::

            sage: Word([]).is_empty()
            True
            sage: Word('a').is_empty()
            False
        """
        return self.length() == 0

    def is_finite(self):
        r"""
        Return ``True``.

        EXAMPLES::

            sage: Word([]).is_finite()
            True
            sage: Word('a').is_finite()
            True
        """
        return True

    def to_integer_word(self):
        r"""
        Return a word over the alphabet ``[0,1,...,self.length()-1]``
        whose letters are in the same relative order as the letters
        of ``self`` in the parent.

        EXAMPLES::

            sage: from itertools import count
            sage: w = Word('abbabaab')
            sage: w.to_integer_word()
            word: 01101001
            sage: w = Word(iter("cacao"), length='finite')
            sage: w.to_integer_word()
            word: 10102

            sage: w = Words([3,2,1])([2,3,3,1])
            sage: w.to_integer_word()
            word: 1002
        """
        from sage.combinat.words.word import Word
        return Word(self.to_integer_list())

    def to_integer_list(self):
        r"""
        Return a list of integers from ``[0,1,...,self.length()-1]`` in the
        same relative order as the letters in ``self`` in the parent.

        EXAMPLES::

            sage: from itertools import count
            sage: w = Word('abbabaab')
            sage: w.to_integer_list()
            [0, 1, 1, 0, 1, 0, 0, 1]
            sage: w = Word(iter("cacao"), length='finite')
            sage: w.to_integer_list()
            [1, 0, 1, 0, 2]
            sage: w = Words([3,2,1])([2,3,3,1])
            sage: w.to_integer_list()
            [1, 0, 0, 2]
        """
        cmp_key = self._parent.sortkey_letters
        ordered_alphabet = sorted(self.letters(), key=cmp_key)
        index = {b: a for a, b in enumerate(ordered_alphabet)}
        return [index[a] for a in self]

    def to_ordered_set_partition(self):
        r"""
        Return the ordered set partition correspond to ``self``.

        If `w` is a finite word of length `n`, then the corresponding
        ordered set partition is an ordered set partition
        `(P_1, P_2, \ldots, P_k)` of `\{1, 2, \ldots, n\}`, where
        each block `P_i` is the set of positions at which the `i`-th
        smallest letter occurring in `w` occurs in `w`.

        EXAMPLES::

            sage: w = Word('abbabaab')
            sage: w.to_ordered_set_partition()
            [{1, 4, 6, 7}, {2, 3, 5, 8}]
            sage: Word([-10, 3, -10, 2]).to_ordered_set_partition()
            [{1, 3}, {4}, {2}]
            sage: Word([]).to_ordered_set_partition()
            []
            sage: Word('aaaaa').to_ordered_set_partition()
            [{1, 2, 3, 4, 5}]
        """
        from sage.combinat.set_partition_ordered import OrderedSetPartition
        return OrderedSetPartition(word_to_ordered_set_partition(self))

    # To fix : do not slice here ! (quite expensive in copy)
    def is_suffix(self, other):
        r"""
        Return ``True`` if ``self`` is a suffix of ``other``, and ``False`` otherwise.

        EXAMPLES::

            sage: w = Word('0123456789')
            sage: y = Word('56789')
            sage: y.is_suffix(w)
            True
            sage: w.is_suffix(y)
            False
            sage: Word('579').is_suffix(w)
            False
            sage: Word().is_suffix(y)
            True
            sage: w.is_suffix(Word())
            False
            sage: Word().is_suffix(Word())
            True
        """
        return self.is_empty() or self == other[-self.length():]

    def is_proper_suffix(self, other):
        r"""
        Return ``True`` if ``self`` is a proper suffix of ``other``, and ``False`` otherwise.

        EXAMPLES::

            sage: Word('23').is_proper_suffix(Word('123'))
            True
            sage: Word('12').is_proper_suffix(Word('12'))
            False
            sage: Word().is_proper_suffix(Word('123'))
            True
            sage: Word('123').is_proper_suffix(Word('12'))
            False
        """
        return self.is_suffix(other) and self.length() < other.length()

    def has_suffix(self, other):
        """
        Test whether ``self`` has ``other`` as a suffix.

        .. NOTE::

           Some word datatype classes, like :class:`WordDatatype_str`,
           override this method.

        INPUT:

        - ``other`` -- a word, or data describing a word

        OUTPUT: boolean

        EXAMPLES::

            sage: w = Word("abbabaabababa")
            sage: u = Word("ababa")
            sage: w.has_suffix(u)
            True
            sage: u.has_suffix(w)
            False
            sage: u.has_suffix("ababa")
            True

        ::

            sage: w = Word([0,1,1,0,1,0,0,1,0,1,0,1,0])
            sage: u = Word([0,1,0,1,0])
            sage: w.has_suffix(u)
            True
            sage: u.has_suffix(w)
            False
            sage: u.has_suffix([0,1,0,1,0])
            True
        """
        from sage.combinat.words.word import Word
        w = Word(other)
        return w.is_suffix(self)

    def is_prefix(self, other):
        r"""
        Return ``True`` if ``self`` is a prefix of ``other``, and ``False`` otherwise.

        EXAMPLES::

            sage: w = Word('0123456789')
            sage: y = Word('012345')
            sage: y.is_prefix(w)
            True
            sage: w.is_prefix(y)
            False
            sage: w.is_prefix(Word())
            False
            sage: Word().is_prefix(w)
            True
            sage: Word().is_prefix(Word())
            True
        """
        return self == other[:self.length()]

    def is_proper_prefix(self, other):
        r"""
        Return ``True`` if ``self`` is a proper prefix of ``other``, and ``False`` otherwise.

        EXAMPLES::

            sage: Word('12').is_proper_prefix(Word('123'))
            True
            sage: Word('12').is_proper_prefix(Word('12'))
            False
            sage: Word().is_proper_prefix(Word('123'))
            True
            sage: Word('123').is_proper_prefix(Word('12'))
            False
            sage: Word().is_proper_prefix(Word())
            False
        """
        return self.is_prefix(other) and self.length() < other.length()

    def has_prefix(self, other):
        r"""
        Test whether ``self`` has ``other`` as a prefix.

        INPUT:

        - ``other`` -- a word, or data describing a word

        OUTPUT: boolean

        EXAMPLES::

            sage: w = Word("abbabaabababa")
            sage: u = Word("abbab")
            sage: w.has_prefix(u)
            True
            sage: u.has_prefix(w)
            False
            sage: u.has_prefix("abbab")
            True

        ::

            sage: w = Word([0,1,1,0,1,0,0,1,0,1,0,1,0])
            sage: u = Word([0,1,1,0,1])
            sage: w.has_prefix(u)
            True
            sage: u.has_prefix(w)
            False
            sage: u.has_prefix([0,1,1,0,1])
            True
        """
        from sage.combinat.words.word import Word
        w = Word(other)
        return w.is_prefix(self)

    def reversal(self):
        r"""
        Return the reversal of ``self``.

        EXAMPLES::

            sage: Word('124563').reversal()
            word: 365421
        """
        return self[::-1]

    @cached_method
    def prefix_function_table(self):
        r"""
        Return a vector containing the length of the proper prefix-suffixes
        for all the non-empty prefixes of ``self``.

        EXAMPLES::

            sage: Word('121321').prefix_function_table()
            [0, 0, 1, 0, 0, 1]
            sage: Word('1241245').prefix_function_table()
            [0, 0, 0, 1, 2, 3, 0]
            sage: Word().prefix_function_table()
            []
        """
        k = 0
        res = [0]*self.length()
        for q in range(1, self.length()):
            while k > 0 and self[k] != self[q]:
                k = res[k-1]
            if self[k] == self[q]:
                k += 1
            res[q] = k
        return res

    @cached_method
    def good_suffix_table(self):
        r"""
        Return a table of the maximum skip you can do in order not to miss
        a possible occurrence of ``self`` in a word.

        This is a part of the Boyer-Moore algorithm to find factors.
        See [BM1977]_.

        EXAMPLES::

            sage: Word('121321').good_suffix_table()
            [5, 5, 5, 5, 3, 3, 1]
            sage: Word('12412').good_suffix_table()
            [3, 3, 3, 3, 3, 1]
        """
        l = self.length()
        p = self.reversal().prefix_function_table()
        res = [l - p[-1]]*(l+1)
        for i in range(1, l+1):
            j = l - p[i - 1]
            if res[j] > (i - p[i-1]):
                res[j] = i - p[i-1]
        return res

    @cached_method
    def suffix_trie(self):
        r"""
        Return the suffix trie of ``self``.

        The *suffix trie* of a finite word `w` is a data structure
        representing the factors of `w`. It is a tree whose edges are
        labelled with letters of `w`, and whose leafs correspond to
        suffixes of `w`.

        Type ``sage.combinat.words.suffix_trees.SuffixTrie?`` for more information.

        EXAMPLES::

            sage: w = Word("cacao")
            sage: w.suffix_trie()
            Suffix Trie of the word: cacao

        ::

            sage: w = Word([0,1,0,1,1])
            sage: w.suffix_trie()
            Suffix Trie of the word: 01011
        """
        from sage.combinat.words.suffix_trees import SuffixTrie
        return SuffixTrie(self)

    def implicit_suffix_tree(self):
        r"""
        Return the implicit suffix tree of ``self``.

        The *suffix tree* of a word `w` is a compactification of the
        suffix trie for `w`. The compactification removes all nodes that have
        exactly one incoming edge and exactly one outgoing edge. It consists of
        two components: a tree and a word. Thus, instead of labelling the edges
        by factors of `w`, we can label them by indices of the occurrence of
        the factors in `w`.

        Type ``sage.combinat.words.suffix_trees.ImplicitSuffixTree?`` for more information.

        EXAMPLES::

            sage: w = Word("cacao")
            sage: w.implicit_suffix_tree()
            Implicit Suffix Tree of the word: cacao

        ::

            sage: w = Word([0,1,0,1,1])
            sage: w.implicit_suffix_tree()
            Implicit Suffix Tree of the word: 01011
        """
        from sage.combinat.words.suffix_trees import ImplicitSuffixTree
        return ImplicitSuffixTree(self)

    @cached_method
    def suffix_tree(self):
        r"""
        Alias for ``implicit_suffix_tree()``.

        EXAMPLES::

            sage: Word('abbabaab').suffix_tree()
            Implicit Suffix Tree of the word: abbabaab
        """
        return self.implicit_suffix_tree()

    def number_of_factors(self, n=None, algorithm='suffix tree'):
        r"""
        Count the number of distinct factors of ``self``.

        INPUT:

        - ``n`` -- integer or ``None``
        - ``algorithm`` -- string (default: ``'suffix tree'``); takes the
          following values:

          - ``'suffix tree'`` -- construct and use the suffix tree of the word
          - ``'naive'`` -- algorithm uses a sliding window

        OUTPUT:

        If ``n`` is an integer, returns the number of distinct factors
        of length ``n``. If ``n`` is ``None``, returns the total number of
        distinct factors.

        EXAMPLES::

            sage: w = Word([1,2,1,2,3])
            sage: w.number_of_factors()
            13
            sage: [w.number_of_factors(i) for i in range(6)]
            [1, 3, 3, 3, 2, 1]

        ::

            sage: w = words.ThueMorseWord()[:100]
            sage: [w.number_of_factors(i) for i in range(10)]
            [1, 2, 4, 6, 10, 12, 16, 20, 22, 24]

        ::

            sage: Word('1213121').number_of_factors()
            22
            sage: Word('1213121').number_of_factors(1)
            3

        ::

            sage: Word('a'*100).number_of_factors()
            101
            sage: Word('a'*100).number_of_factors(77)
            1

        ::

            sage: Word().number_of_factors()
            1
            sage: Word().number_of_factors(17)
            0

        ::

            sage: blueberry = Word("blueberry")
            sage: blueberry.number_of_factors()
            43
            sage: [blueberry.number_of_factors(i) for i in range(10)]
            [1, 6, 8, 7, 6, 5, 4, 3, 2, 1]
        """
        if algorithm == 'suffix tree':
            return ZZ(self.suffix_tree().number_of_factors(n))
        elif algorithm == 'naive':
            return ZZ(len(self.factor_set(n, algorithm='naive')))
        else:
            raise ValueError(f'Unknown algorithm (={algorithm})')

    def factor_iterator(self, n=None):
        r"""
        Generate distinct factors of ``self``.

        INPUT:

        - ``n`` -- integer or ``None``

        OUTPUT:

        If ``n`` is an integer, returns an iterator over all distinct
        factors of length ``n``. If ``n`` is ``None``, returns an iterator
        generating all distinct factors.

        EXAMPLES::

            sage: w = Word('1213121')
            sage: sorted( w.factor_iterator(0) )
            [word: ]
            sage: sorted( w.factor_iterator(10) )
            []
            sage: sorted( w.factor_iterator(1) )
            [word: 1, word: 2, word: 3]
            sage: sorted( w.factor_iterator(4) )
            [word: 1213, word: 1312, word: 2131, word: 3121]
            sage: sorted( w.factor_iterator() )
            [word: , word: 1, word: 12, word: 121, word: 1213, word: 12131, word: 121312, word: 1213121, word: 13, word: 131, word: 1312, word: 13121, word: 2, word: 21, word: 213, word: 2131, word: 21312, word: 213121, word: 3, word: 31, word: 312, word: 3121]

        ::

            sage: u = Word([1,2,1,2,3])
            sage: sorted( u.factor_iterator(0) )
            [word: ]
            sage: sorted( u.factor_iterator(10) )
            []
            sage: sorted( u.factor_iterator(1) )
            [word: 1, word: 2, word: 3]
            sage: sorted( u.factor_iterator(5) )
            [word: 12123]
            sage: sorted( u.factor_iterator() )
            [word: , word: 1, word: 12, word: 121, word: 1212, word: 12123, word: 123, word: 2, word: 21, word: 212, word: 2123, word: 23, word: 3]

        ::

            sage: xxx = Word("xxx")
            sage: sorted( xxx.factor_iterator(0) )
            [word: ]
            sage: sorted( xxx.factor_iterator(4) )
            []
            sage: sorted( xxx.factor_iterator(2) )
            [word: xx]
            sage: sorted( xxx.factor_iterator() )
            [word: , word: x, word: xx, word: xxx]

        ::

            sage: e = Word()
            sage: sorted( e.factor_iterator(0) )
            [word: ]
            sage: sorted( e.factor_iterator(17) )
            []
            sage: sorted( e.factor_iterator() )
            [word: ]

        TESTS::

            sage: type( Word('cacao').factor_iterator() )
            <class 'generator'>
        """
        return self.suffix_tree().factor_iterator(n)

    def factor_complexity(self, n):
        r"""
        Return the number of distinct factors of length ``n`` of ``self``.

        INPUT:

        - ``n`` -- the length of the factors

        EXAMPLES::

            sage: w = words.FibonacciWord()[:100]
            sage: [w.factor_complexity(i) for i in range(20)]
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]

        ::

            sage: w = words.ThueMorseWord()[:1000]
            sage: [w.factor_complexity(i) for i in range(20)]
            [1, 2, 4, 6, 10, 12, 16, 20, 22, 24, 28, 32, 36, 40, 42, 44, 46, 48, 52, 56]
        """
        return len(list(self.factor_iterator(n)))

    def factor_set(self, n=None, algorithm='suffix tree'):
        r"""
        Return the set of factors (of length ``n``) of ``self``.

        INPUT:

        - ``n`` -- integer or ``None`` (default: ``None``)
        - ``algorithm`` -- string (default: ``'suffix tree'``), takes the
          following values:

          - ``'suffix tree'`` -- construct and use the suffix tree of the word
          - ``'naive'`` -- algorithm uses a sliding window

        OUTPUT:

        If ``n`` is an integer, returns the set of all distinct
        factors of length ``n``. If ``n`` is ``None``, returns the set
        of all distinct factors.

        EXAMPLES::

            sage: w = Word('121')
            sage: sorted(w.factor_set())
            [word: , word: 1, word: 12, word: 121, word: 2, word: 21]
            sage: sorted(w.factor_set(algorithm='naive'))
            [word: , word: 1, word: 12, word: 121, word: 2, word: 21]

        ::

            sage: w = Word('1213121')
            sage: for i in range(w.length()): sorted(w.factor_set(i))
            [word: ]
            [word: 1, word: 2, word: 3]
            [word: 12, word: 13, word: 21, word: 31]
            [word: 121, word: 131, word: 213, word: 312]
            [word: 1213, word: 1312, word: 2131, word: 3121]
            [word: 12131, word: 13121, word: 21312]
            [word: 121312, word: 213121]

        ::

            sage: w = Word([1,2,1,2,3])
            sage: s = w.factor_set()
            sage: sorted(s)
            [word: , word: 1, word: 12, word: 121, word: 1212, word: 12123, word: 123, word: 2, word: 21, word: 212, word: 2123, word: 23, word: 3]

        TESTS::

            sage: w = Word("xx")
            sage: s = w.factor_set()
            sage: sorted(s)
            [word: , word: x, word: xx]

        ::

            sage: Set(Word().factor_set())
            {word: }

        ::

            sage: w = Word(range(10), alphabet=range(10))
            sage: S1 = w.factor_set(3, algorithm='suffix tree')
            sage: S2 = w.factor_set(3, algorithm='naive')
            sage: S1 == S2
            True
        """
        if algorithm == 'suffix tree':
            return Set(self.factor_iterator(n))
        elif algorithm == 'naive':
            if n is None:
                S = {self[0:0]}
                for n in range(1, self.length()+1):
                    for i in range(self.length()-n+1):
                        S.add(self[i:i+n])
                return Set(S)
            else:
                S = set()
                for i in range(self.length()-n+1):
                    S.add(self[i:i+n])
                return Set(S)
        else:
            raise ValueError(f'Unknown algorithm (={algorithm})')

    def topological_entropy(self, n):
        r"""
        Return the topological entropy for the factors of length ``n``.

        The topological entropy of a sequence `u` is defined as the
        exponential growth rate of the complexity of `u` as the length
        increases: `H_{top}(u)=\lim_{n\to\infty}\frac{\log_d(p_u(n))}{n}`
        where `d` denotes the cardinality of the alphabet and `p_u(n)` is
        the complexity function, i.e. the number of factors of length `n`
        in the sequence `u` [Fog2002]_.

        INPUT:

        - ``self`` -- a word defined over a finite alphabet
        - ``n`` -- positive integer

        OUTPUT: real number (a symbolic expression)

        EXAMPLES::

            sage: W = Words([0, 1])
            sage: w = W([0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1])
            sage: t = w.topological_entropy(3); t                                       # needs sage.symbolic
            1/3*log(7)/log(2)
            sage: n(t)                                                                  # needs sage.symbolic
            0.935784974019201

        ::

            sage: w = words.ThueMorseWord()[:100]
            sage: topo = w.topological_entropy
            sage: for i in range(0, 41, 5):                                             # needs sage.symbolic
            ....:     print("{} {}".format(i, n(topo(i), digits=5)))
            0 1.0000
            5 0.71699
            10 0.48074
            15 0.36396
            20 0.28774
            25 0.23628
            30 0.20075
            35 0.17270
            40 0.14827

        If no alphabet is specified, an error is raised::

            sage: w = Word(range(20))
            sage: w.topological_entropy(3)
            Traceback (most recent call last):
            ...
            TypeError: The word must be defined over a finite alphabet

        The following is ok::

            sage: W = Words(range(20))
            sage: w = W(range(20))
            sage: w.topological_entropy(3)                                              # needs sage.symbolic
            1/3*log(18)/log(20)
        """
        d = self.parent().alphabet().cardinality()
        if d is Infinity:
            raise TypeError("The word must be defined over a finite alphabet")
        if n == 0:
            return 1
        pn = self.number_of_factors(n)
        from sage.functions.log import log
        return log(pn, base=d)/n

    def rauzy_graph(self, n):
        r"""
        Return the Rauzy graph of the factors of length ``n`` of ``self``.

        The vertices are the factors of length `n` and there is an edge from
        `u` to `v` if `ua = bv` is a factor of length `n+1` for some letters
        `a` and `b`.

        INPUT:

        - ``n`` -- integer

        EXAMPLES::

            sage: w = Word(range(10)); w
            word: 0123456789
            sage: g = w.rauzy_graph(3); g                                               # needs sage.graphs
            Looped digraph on 8 vertices
            sage: WordOptions(identifier='')
            sage: g.vertices(sort=True)                                                 # needs sage.graphs
            [012, 123, 234, 345, 456, 567, 678, 789]
            sage: g.edges(sort=True)                                                    # needs sage.graphs
            [(012, 123, 3),
             (123, 234, 4),
             (234, 345, 5),
             (345, 456, 6),
             (456, 567, 7),
             (567, 678, 8),
             (678, 789, 9)]
            sage: WordOptions(identifier='word: ')

        ::

            sage: f = words.FibonacciWord()[:100]
            sage: f.rauzy_graph(8)                                                      # needs sage.graphs
            Looped digraph on 9 vertices

        ::

            sage: w = Word('1111111')
            sage: g = w.rauzy_graph(3)                                                  # needs sage.graphs
            sage: g.edges(sort=True)                                                    # needs sage.graphs
            [(word: 111, word: 111, word: 1)]

        ::

            sage: w = Word('111')
            sage: for i in range(5): w.rauzy_graph(i)                                   # needs sage.graphs
            Looped multi-digraph on 1 vertex
            Looped digraph on 1 vertex
            Looped digraph on 1 vertex
            Looped digraph on 1 vertex
            Looped digraph on 0 vertices

        Multi-edges are allowed for the empty word::

            sage: W = Words('abcde')
            sage: w = W('abc')
            sage: w.rauzy_graph(0)                                                      # needs sage.graphs
            Looped multi-digraph on 1 vertex
            sage: _.edges(sort=True)                                                    # needs sage.graphs
            [(word: , word: , word: a),
             (word: , word: , word: b),
             (word: , word: , word: c)]
        """
        from sage.graphs.digraph import DiGraph
        multiedges = n == 0
        g = DiGraph(loops=True, multiedges=multiedges)
        if n == self.length():
            g.add_vertex(self)
        else:
            for w in self.factor_iterator(n+1):
                u = w[:-1]
                v = w[1:]
                a = w[-1:]
                g.add_edge(u, v, a)
        return g

    def reduced_rauzy_graph(self, n):
        r"""
        Return the reduced Rauzy graph of order ``n`` of ``self``.

        INPUT:

        - ``n`` -- nonnegative integer; every vertex of a reduced
          Rauzy graph of order ``n`` is a factor of length ``n`` of ``self``

        OUTPUT: a looped multi-digraph

        DEFINITION:

        For infinite periodic words (resp. for finite words of type `u^i
        u[0:j]`), the reduced Rauzy graph of order `n` (resp. for `n`
        smaller or equal to `(i-1)|u|+j`) is the directed graph whose
        unique vertex is the prefix `p` of length `n` of ``self`` and which has
        an only edge which is a loop on `p` labelled by `w[n+1:|w|] p`
        where `w` is the unique return word to `p`.

        In other cases, it is the directed graph defined as followed.  Let
        `G_n` be the Rauzy graph of order `n` of ``self``. The vertices are the
        vertices of `G_n` that are either special or not prolongable to the
        right or to the left. For each couple (`u`, `v`) of such vertices
        and each directed path in `G_n` from `u` to `v` that contains no
        other vertices that are special, there is an edge from `u` to `v`
        in the reduced Rauzy graph of order `n` whose label is the label of
        the path in `G_n`.

        .. NOTE::

            In the case of infinite recurrent non-periodic words, this
            definition corresponds to the following one that can be
            found in [BDLGZ2009]_ and [BPS2008]_ where a simple path is a
            path that begins with a special factor, ends with a
            special factor and contains no other vertices that are
            special:

            The reduced Rauzy graph of factors of length `n` is obtained
            from `G_n` by replacing each simple path `P=v_1 v_2 ...
            v_{\ell}` with an edge `v_1 v_{\ell}` whose label is the
            concatenation of the labels of the edges of `P`.

        EXAMPLES::

            sage: w = Word(range(10)); w
            word: 0123456789
            sage: g = w.reduced_rauzy_graph(3); g                                       # needs sage.graphs
            Looped multi-digraph on 2 vertices
            sage: g.vertices(sort=True)                                                 # needs sage.graphs
            [word: 012, word: 789]
            sage: g.edges(sort=True)                                                    # needs sage.graphs
            [(word: 012, word: 789, word: 3456789)]

        For the Fibonacci word::

            sage: f = words.FibonacciWord()[:100]
            sage: g = f.reduced_rauzy_graph(8);g                                        # needs sage.graphs
            Looped multi-digraph on 2 vertices
            sage: g.vertices(sort=True)                                                 # needs sage.graphs
            [word: 01001010, word: 01010010]
            sage: g.edges(sort=True)                                                    # needs sage.graphs
            [(word: 01001010, word: 01010010, word: 010),
             (word: 01010010, word: 01001010, word: 01010),
             (word: 01010010, word: 01001010, word: 10)]

        For periodic words::

            sage: from itertools import cycle
            sage: w = Word(cycle('abcd'))[:100]
            sage: g = w.reduced_rauzy_graph(3)                                          # needs sage.graphs
            sage: g.edges(sort=True)                                                    # needs sage.graphs
            [(word: abc, word: abc, word: dabc)]

        ::

            sage: w = Word('111')
            sage: for i in range(5): w.reduced_rauzy_graph(i)                           # needs sage.graphs
            Looped digraph on 1 vertex
            Looped digraph on 1 vertex
            Looped digraph on 1 vertex
            Looped multi-digraph on 1 vertex
            Looped multi-digraph on 0 vertices

        For ultimately periodic words::

            sage: sigma = WordMorphism('a->abcd,b->cd,c->cd,d->cd')
            sage: w = sigma.fixed_point('a')[:100]; w                                   # needs sage.modules
            word: abcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd...
            sage: g = w.reduced_rauzy_graph(5)                                          # needs sage.graphs
            sage: g.vertices(sort=True)                                                 # needs sage.graphs
            [word: abcdc, word: cdcdc]
            sage: g.edges(sort=True)                                                    # needs sage.graphs
            [(word: abcdc, word: cdcdc, word: dc), (word: cdcdc, word: cdcdc, word: dc)]

        AUTHOR:

        Julien Leroy (March 2010): initial version
        """
        from sage.graphs.digraph import DiGraph
        from copy import copy
        g = copy(self.rauzy_graph(n))
        # Otherwise it changes the rauzy_graph function.
        l = [v for v in g if g.in_degree(v) == 1 == g.out_degree(v)]
        if g.num_verts() != 0 and len(l) == g.num_verts():
            # In this case, the Rauzy graph is simply a cycle.
            g = DiGraph()
            g.allow_loops(True)
            g.add_vertex(self[:n])
            g.add_edge(self[:n], self[:n], self[n:n + len(l)])
        else:
            g.allow_loops(True)
            g.allow_multiple_edges(True)
            for v in l:
                [i] = g.neighbors_in(v)
                [o] = g.neighbors_out(v)
                g.add_edge(i, o, g.edge_label(i, v)[0]*g.edge_label(v, o)[0])
                g.delete_vertex(v)
        return g

    def left_special_factors_iterator(self, n=None):
        r"""
        Return an iterator over the left special factors (of length ``n``).

        A factor `u` of a word `w` is *left special* if there are
        two distinct letters `a` and `b` such that `au` and `bu`
        are factors of `w`.

        INPUT:

        - ``n`` -- integer (default: ``None``); if ``None``, it returns
          an iterator over all left special factors

        EXAMPLES::

            sage: alpha, beta, x = 0.54, 0.294, 0.1415
            sage: w = words.CodingOfRotationWord(alpha, beta, x)[:40]
            sage: sorted(w.left_special_factors_iterator(3))
            [word: 000, word: 010]
            sage: sorted(w.left_special_factors_iterator(4))
            [word: 0000, word: 0101]
            sage: sorted(w.left_special_factors_iterator(5))
            [word: 00000, word: 01010]
        """
        if n is None:
            for i in range(self.length()):
                yield from self.left_special_factors_iterator(i)
        else:
            left_extensions = defaultdict(set)
            for w in self.factor_iterator(n+1):
                v = w[1:]
                left_extensions[v].add(w[0])
            for v in left_extensions:
                if len(left_extensions[v]) > 1:
                    yield v

    def left_special_factors(self, n=None):
        r"""
        Return the left special factors (of length ``n``).

        A factor `u` of a word `w` is *left special* if there are
        two distinct letters `a` and `b` such that `au` and `bu`
        are factors of `w`.

        INPUT:

        - ``n`` -- integer (default: ``None``); if ``None``, it
          returns all left special factors

        OUTPUT: list of words

        EXAMPLES::

            sage: alpha, beta, x = 0.54, 0.294, 0.1415
            sage: w = words.CodingOfRotationWord(alpha, beta, x)[:40]
            sage: for i in range(5):
            ....:     print("{} {}".format(i, sorted(w.left_special_factors(i))))
            0 [word: ]
            1 [word: 0]
            2 [word: 00, word: 01]
            3 [word: 000, word: 010]
            4 [word: 0000, word: 0101]
        """
        return list(self.left_special_factors_iterator(n))

    def right_special_factors_iterator(self, n=None):
        r"""
        Return an iterator over the right special factors (of length ``n``).

        A factor `u` of a word `w` is *right special* if there are
        two distinct letters `a` and `b` such that `ua` and `ub`
        are factors of `w`.

        INPUT:

        - ``n`` -- integer (default: ``None``); if ``None``, it returns
          an iterator over all right special factors

        EXAMPLES::

            sage: alpha, beta, x = 0.61, 0.54, 0.3
            sage: w = words.CodingOfRotationWord(alpha, beta, x)[:40]
            sage: sorted(w.right_special_factors_iterator(3))
            [word: 010, word: 101]
            sage: sorted(w.right_special_factors_iterator(4))
            [word: 0101, word: 1010]
            sage: sorted(w.right_special_factors_iterator(5))
            [word: 00101, word: 11010]
        """
        if n is None:
            for i in range(self.length()):
                yield from self.right_special_factors_iterator(i)
        else:
            right_extensions = defaultdict(set)
            for w in self.factor_iterator(n+1):
                v = w[:-1]
                right_extensions[v].add(w[-1])
            for v in right_extensions:
                if len(right_extensions[v]) > 1:
                    yield v

    def right_special_factors(self, n=None):
        r"""
        Return the right special factors (of length ``n``).

        A factor `u` of a word `w` is *right special* if there are
        two distinct letters `a` and `b` such that `ua` and `ub`
        are factors of `w`.

        INPUT:

        - ``n`` -- integer (default: ``None``); if ``None``, it returns
          all right special factors

        OUTPUT: list of words

        EXAMPLES::

            sage: w = words.ThueMorseWord()[:30]
            sage: for i in range(5):
            ....:     print("{} {}".format(i, sorted(w.right_special_factors(i))))
            0 [word: ]
            1 [word: 0, word: 1]
            2 [word: 01, word: 10]
            3 [word: 001, word: 010, word: 101, word: 110]
            4 [word: 0110, word: 1001]
        """
        return list(self.right_special_factors_iterator(n))

    def bispecial_factors_iterator(self, n=None):
        r"""
        Return an iterator over the bispecial factors (of length ``n``).

        A factor `u` of a word `w` is *bispecial* if it is right special
        and left special.

        INPUT:

        - ``n`` -- integer (default: ``None``); if ``None``, it returns
          an iterator over all bispecial factors

        EXAMPLES::

            sage: w = words.ThueMorseWord()[:30]
            sage: for i in range(10):
            ....:     for u in sorted(w.bispecial_factors_iterator(i)):
            ....:         print("{} {}".format(i,u))
            0
            1 0
            1 1
            2 01
            2 10
            3 010
            3 101
            4 0110
            4 1001
            6 011001
            6 100110
            8 10010110

        ::

            sage: key = lambda u : (len(u), u)
            sage: for u in sorted(w.bispecial_factors_iterator(), key=key): u
            word:
            word: 0
            word: 1
            word: 01
            word: 10
            word: 010
            word: 101
            word: 0110
            word: 1001
            word: 011001
            word: 100110
            word: 10010110
        """
        if n is None:
            for i in range(self.length()):
                yield from self.bispecial_factors_iterator(i)
        else:
            left_extensions = defaultdict(set)
            right_extensions = defaultdict(set)
            for w in self.factor_iterator(n + 2):
                v = w[1:-1]
                left_extensions[v].add(w[0])
                right_extensions[v].add(w[-1])
            for v in left_extensions:
                if (len(left_extensions[v]) > 1 and
                        len(right_extensions[v]) > 1):
                    yield v

    def bispecial_factors(self, n=None):
        r"""
        Return the bispecial factors (of length ``n``).

        A factor `u` of a word `w` is *bispecial* if it is right special
        and left special.

        INPUT:

        - ``n`` -- integer (default: ``None``); if ``None``, it returns
          all bispecial factors

        OUTPUT: list of words

        EXAMPLES::

            sage: w = words.FibonacciWord()[:30]
            sage: w.bispecial_factors()
            [word: , word: 0, word: 010, word: 010010, word: 01001010010]

        ::

            sage: w = words.ThueMorseWord()[:30]
            sage: for i in range(10):
            ....:     print("{} {}".format(i, sorted(w.bispecial_factors(i))))
            0 [word: ]
            1 [word: 0, word: 1]
            2 [word: 01, word: 10]
            3 [word: 010, word: 101]
            4 [word: 0110, word: 1001]
            5 []
            6 [word: 011001, word: 100110]
            7 []
            8 [word: 10010110]
            9 []
        """
        return list(self.bispecial_factors_iterator(n))

    def number_of_left_special_factors(self, n):
        r"""
        Return the number of left special factors of length `n`.

        A factor `u` of a word `w` is *left special* if there are
        two distinct letters `a` and `b` such that `au` and `bu`
        are factors of `w`.

        INPUT:

        - ``n`` -- integer

        OUTPUT: nonnegative integer

        EXAMPLES::

            sage: w = words.FibonacciWord()[:100]
            sage: [w.number_of_left_special_factors(i) for i in range(10)]
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

        ::

            sage: w = words.ThueMorseWord()[:100]
            sage: [w.number_of_left_special_factors(i) for i in range(10)]
            [1, 2, 2, 4, 2, 4, 4, 2, 2, 4]
        """
        it = self.left_special_factors_iterator(n)
        return sum(1 for _ in it)

    def number_of_right_special_factors(self, n):
        r"""
        Return the number of right special factors of length ``n``.

        A factor `u` of a word `w` is *right special* if there are
        two distinct letters `a` and `b` such that `ua` and `ub`
        are factors of `w`.

        INPUT:

        - ``n`` -- integer

        OUTPUT: nonnegative integer

        EXAMPLES::

            sage: w = words.FibonacciWord()[:100]
            sage: [w.number_of_right_special_factors(i) for i in range(10)]
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

        ::

            sage: w = words.ThueMorseWord()[:100]
            sage: [w.number_of_right_special_factors(i) for i in range(10)]
            [1, 2, 2, 4, 2, 4, 4, 2, 2, 4]
        """
        it = self.right_special_factors_iterator(n)
        return sum(1 for _ in it)

    def commutes_with(self, other):
        r"""
        Return ``True`` if ``self`` commutes with ``other``, and ``False`` otherwise.

        EXAMPLES::

            sage: Word('12').commutes_with(Word('12'))
            True
            sage: Word('12').commutes_with(Word('11'))
            False
            sage: Word().commutes_with(Word('21'))
            True
        """
        return (self * other) == (other * self)

    def conjugate(self, pos):
        r"""
        Return the conjugate at ``pos`` of ``self``.

        ``pos`` can be any integer, the distance used is the modulo by the length
        of ``self``.

        EXAMPLES::

            sage: Word('12112').conjugate(1)
            word: 21121
            sage: Word().conjugate(2)
            word:
            sage: Word('12112').conjugate(8)
            word: 12121
            sage: Word('12112').conjugate(-1)
            word: 21211
        """
        if self.is_empty():
            return self
        pos_mod = pos % self.length()
        return self[pos_mod:] * self[:pos_mod]

    def _conjugates_list(self):
        r"""
        Return the list of conjugates of ``self``, ordered from the `0`-th to the
        `(L-1)`-st conjugate, where `L` is the length of ``self``.

        TESTS::

            sage: Word('cbbca')._conjugates_list()
            [word: cbbca, word: bbcac, word: bcacb, word: cacbb, word: acbbc]
            sage: Word('abcabc')._conjugates_list()
            [word: abcabc,
             word: bcabca,
             word: cabcab,
             word: abcabc,
             word: bcabca,
             word: cabcab]
            sage: Word()._conjugates_list()
            [word: ]
            sage: Word('a')._conjugates_list()
            [word: a]
        """
        S = [self]
        S.extend(self.conjugate(i) for i in range(1, self.length()))
        return S

    def conjugates_iterator(self):
        r"""
        Return an iterator over the conjugates of ``self``.

        EXAMPLES::

            sage: it = Word(range(4)).conjugates_iterator()
            sage: for w in it: w
            word: 0123
            word: 1230
            word: 2301
            word: 3012
        """
        yield self
        for i in range(1, self.primitive_length()):
            yield self.conjugate(i)

    def conjugates(self):
        r"""
        Return the list of unique conjugates of ``self``.

        EXAMPLES::

            sage: Word(range(6)).conjugates()
            [word: 012345,
             word: 123450,
             word: 234501,
             word: 345012,
             word: 450123,
             word: 501234]
            sage: Word('cbbca').conjugates()
            [word: cbbca, word: bbcac, word: bcacb, word: cacbb, word: acbbc]

        The result contains each conjugate only once::

            sage: Word('abcabc').conjugates()
            [word: abcabc, word: bcabca, word: cabcab]

        TESTS::

            sage: Word().conjugates()
            [word: ]
            sage: Word('a').conjugates()
            [word: a]
        """
        return list(self.conjugates_iterator())

    def conjugate_position(self, other):
        r"""
        Return the position where ``self`` is conjugate with ``other``.
        Return ``None`` if there is no such position.

        EXAMPLES::

            sage: Word('12113').conjugate_position(Word('31211'))
            1
            sage: Word('12131').conjugate_position(Word('12113')) is None
            True
            sage: Word().conjugate_position(Word('123')) is None
            True

        TESTS:

        We check that :issue:`11128` is fixed::

            sage: w = Word([0,0,1,0,2,1])
            sage: [w.conjugate(i).conjugate_position(w) for i in range(w.length())]
            [0, 1, 2, 3, 4, 5]
        """
        if self.length() != other.length():
            return None
        other_square = other * other
        pos = other_square.find(self)
        return pos if pos != -1 else None

    def is_conjugate_with(self, other):
        r"""
        Return ``True`` if ``self`` is a conjugate of ``other``, and ``False`` otherwise.

        INPUT:

        - ``other`` -- a finite word

        OUTPUT: boolean

        EXAMPLES::

            sage: w = Word([0..20])
            sage: z = Word([7..20] + [0..6])
            sage: w
            word: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20
            sage: z
            word: 7,8,9,10,11,12,13,14,15,16,17,18,19,20,0,1,2,3,4,5,6
            sage: w.is_conjugate_with(z)
            True
            sage: z.is_conjugate_with(w)
            True
            sage: u = Word([4]*21)
            sage: u.is_conjugate_with(w)
            False
            sage: u.is_conjugate_with(z)
            False

        Both words must be finite::

            sage: w = Word(iter([2]*100),length='unknown')
            sage: z = Word([2]*100)
            sage: z.is_conjugate_with(w) #TODO: Not implemented for word of unknown length
            True
            sage: wf = Word(iter([2]*100),length='finite')
            sage: z.is_conjugate_with(wf)
            True
            sage: wf.is_conjugate_with(z)
            True

        TESTS::

            sage: Word('11213').is_conjugate_with(Word('31121'))
            True
            sage: Word().is_conjugate_with(Word('123'))
            False
            sage: Word('112131').is_conjugate_with(Word('11213'))
            False
            sage: Word('12131').is_conjugate_with(Word('11213'))
            True

        We make sure that :issue:`11128` is fixed::

            sage: Word('abaa').is_conjugate_with(Word('aaba'))
            True
            sage: Word('aaba').is_conjugate_with(Word('abaa'))
            True
        """
        return self.length() == other.length() and self.is_factor(other * other)

    def is_cadence(self, seq):
        r"""
        Return ``True`` if ``seq`` is a cadence of ``self``, and ``False`` otherwise.

        A *cadence* is an increasing sequence of indexes that all map to
        the same letter.

        EXAMPLES::

            sage: Word('121132123').is_cadence([0, 2, 6])
            True
            sage: Word('121132123').is_cadence([0, 1, 2])
            False
            sage: Word('121132123').is_cadence([])
            True
        """
        if not seq:
            return True
        try:
            it = iter(self)
            s = next(islice(it, int(seq[0]), None))
            for i in range(1, len(seq)):
                steps = seq[i] - seq[i - 1]
                for n in range(steps - 1):
                    next(it)
                if next(it) != s:
                    return False
        except StopIteration:
            return False
        return True

    def longest_forward_extension(self, x, y):
        r"""
        Compute the length of the longest factor of ``self`` that
        starts at ``x`` and that matches a factor that starts at ``y``.

        INPUT:

        - ``x``, ``y`` -- positions in ``self``

        EXAMPLES::

            sage: w = Word('0011001')
            sage: w.longest_forward_extension(0, 4)
            3
            sage: w.longest_forward_extension(0, 2)
            0

        The method also accepts negative positions indicating the distance from
        the end of the word (in order to be consist with how negative indices
        work with lists). For instance, for a word of length `7`, using
        positions `-3` and `2` is the same as using positions `4` and `2`::

            sage: w.longest_forward_extension(1, -2)
            2
            sage: w.longest_forward_extension(4, -3)
            3

        TESTS::

            sage: w = Word('0011001')
            sage: w.longest_forward_extension(-10, 2)
            Traceback (most recent call last):
            ...
            ValueError: x and y must be valid positions in self
        """
        length = self.length()
        if not (-length <= x < length and -length <= y < length):
            raise ValueError("x and y must be valid positions in self")
        if x < 0:
            x = x + length
        if y < 0:
            y = y + length
        l = 0
        while x < length and y < length and self[x] == self[y]:
            l += 1
            x += 1
            y += 1
        return l

    def longest_backward_extension(self, x, y):
        r"""
        Compute the length of the longest factor of ``self`` that
        ends at ``x`` and that matches a factor that ends at ``y``.

        INPUT:

        - ``x``, ``y`` -- positions in ``self``

        EXAMPLES::

            sage: w = Word('0011001')
            sage: w.longest_backward_extension(6, 2)
            3
            sage: w.longest_backward_extension(1, 4)
            1
            sage: w.longest_backward_extension(1, 3)
            0

        The method also accepts negative positions indicating the distance from
        the end of the word (in order to be consist with how negative indices
        work with lists). For instance, for a word of length `7`, using
        positions `6` and `-5` is the same as using positions `6` and `2`::

            sage: w.longest_backward_extension(6, -5)
            3
            sage: w.longest_backward_extension(-6, 4)
            1

        TESTS::

            sage: w = Word('0011001')
            sage: w.longest_backward_extension(4, 23)
            Traceback (most recent call last):
            ...
            ValueError: x and y must be valid positions in self
            sage: w.longest_backward_extension(-9, 4)
            Traceback (most recent call last):
            ...
            ValueError: x and y must be valid positions in self
        """
        length = self.length()
        if not (-length <= x < length and -length <= y < length):
            raise ValueError("x and y must be valid positions in self")
        if x < 0:
            x = x + length
        if y < 0:
            y = y + length
        l = 0
        while x >= 0 and y >= 0 and self[x] == self[y]:
            l += 1
            x -= 1
            y -= 1
        return l

    def longest_common_suffix(self, other):
        r"""
        Return the longest common suffix of ``self`` and ``other``.

        EXAMPLES::

            sage: w = Word('112345678')
            sage: u = Word('1115678')
            sage: w.longest_common_suffix(u)
            word: 5678
            sage: u.longest_common_suffix(u)
            word: 1115678
            sage: u.longest_common_suffix(w)
            word: 5678
            sage: w.longest_common_suffix(w)
            word: 112345678
            sage: y = Word('549332345')
            sage: w.longest_common_suffix(y)
            word:

        TESTS:

        With the empty word::

            sage: w.longest_common_suffix(Word())
            word:
            sage: Word().longest_common_suffix(w)
            word:
            sage: Word().longest_common_suffix(Word())
            word:

        With an infinite word::

            sage: t = words.ThueMorseWord('ab')
            sage: w.longest_common_suffix(t)
            Traceback (most recent call last):
            ...
            TypeError: other must be a finite word
        """
        if not isinstance(other, FiniteWord_class):
            raise TypeError("other must be a finite word")

        if self.is_empty():
            return self
        if other.is_empty():
            return other

        iter = enumerate(zip(reversed(self), reversed(other)))
        i, (b, c) = next(iter)
        if b != c:
            # In this case, return the empty word
            return self[:0]

        for i, (b, c) in iter:
            if b != c:
                return self[-i:]
        else:
            return self[-i-1:]

    def is_palindrome(self, f=None):
        r"""
        Return ``True`` if ``self`` is a palindrome (or a ``f``-palindrome), and
        ``False`` otherwise.

        Let `f : \Sigma \rightarrow \Sigma` be an involution that extends
        to a morphism on `\Sigma^*`. We say that `w\in\Sigma^*` is a
        *`f`-palindrome* if `w=f(\tilde{w})` [Lab2008]_. Also called
        *`f`-pseudo-palindrome* [AZZ2005]_.

        INPUT:

        - ``f`` -- involution (default: ``None``) on the alphabet of ``self``.
          It must be callable on letters as well as words (e.g.
          :class:`~sage.combinat.words.morphism.WordMorphism`). The
          default value corresponds to usual palindromes, i.e., ``f``
          equal to the identity.

        EXAMPLES::

            sage: Word('esope reste ici et se repose').is_palindrome()
            False
            sage: Word('esoperesteicietserepose').is_palindrome()
            True
            sage: Word('I saw I was I').is_palindrome()
            True
            sage: Word('abbcbba').is_palindrome()
            True
            sage: Word('abcbdba').is_palindrome()
            False

        Some `f`-palindromes::

            sage: f = WordMorphism('a->b,b->a')
            sage: Word('aababb').is_palindrome(f)
            True

        ::

            sage: f = WordMorphism('a->b,b->a,c->c')
            sage: Word('abacbacbab').is_palindrome(f)
            True

        ::

            sage: f = WordMorphism({'a':'b','b':'a'})
            sage: Word('aababb').is_palindrome(f)
            True

        ::

            sage: f = WordMorphism({0:[1],1:[0]})
            sage: w = words.ThueMorseWord()[:8]; w
            word: 01101001
            sage: w.is_palindrome(f)
            True

        The word must be in the domain of the involution::

            sage: f = WordMorphism('a->a')
            sage: Word('aababb').is_palindrome(f)
            Traceback (most recent call last):
            ...
            KeyError: 'b'

        TESTS:

        If the given involution is not an involution::

            sage: f = WordMorphism('a->b,b->b')
            sage: Word('abab').is_palindrome(f)
            Traceback (most recent call last):
            ...
            TypeError: self (=a->b, b->b) is not an endomorphism

        ::

            sage: Y = Word
            sage: Y().is_palindrome()
            True
            sage: Y('a').is_palindrome()
            True
            sage: Y('ab').is_palindrome()
            False
            sage: Y('aba').is_palindrome()
            True
            sage: Y('aa').is_palindrome()
            True
            sage: E = WordMorphism('a->b,b->a')
            sage: Y().is_palindrome(E)
            True
            sage: Y('a').is_palindrome(E)
            False
            sage: Y('ab').is_palindrome(E)
            True
            sage: Y('aa').is_palindrome(E)
            False
            sage: Y('aba').is_palindrome(E)
            False
            sage: Y('abab').is_palindrome(E)
            True
        """
        l = self.length()
        if f is None:
            return self[:l//2] == self[l//2 + l % 2:].reversal()
        else:
            from sage.combinat.words.morphism import WordMorphism
            if not isinstance(f, WordMorphism):
                f = WordMorphism(f)
            if not f.is_involution():
                raise ValueError("f must be an involution")
            return self[:l//2 + l % 2] == f(self[l//2:].reversal())

    def lps(self, f=None, l=None):
        r"""
        Return the longest palindromic (or ``f``-palindromic) suffix of ``self``.

        INPUT:

        - ``f`` -- involution (default: ``None``) on the alphabet of ``self``.
          It must be callable on letters as well as words (e.g. ``WordMorphism``)
        - ``l`` -- integer (default: ``None``); the length of the longest
          palindrome suffix of ``self[:-1]``, if known

        OUTPUT: word; if ``f`` is ``None``, the longest palindromic suffix of
        ``self``. Otherwise, the longest ``f``-palindromic suffix of ``self``.

        EXAMPLES::

            sage: Word('0111').lps()
            word: 111
            sage: Word('011101').lps()
            word: 101
            sage: Word('6667').lps()
            word: 7
            sage: Word('abbabaab').lps()
            word: baab
            sage: Word().lps()
            word:
            sage: f = WordMorphism('a->b,b->a')
            sage: Word('abbabaab').lps(f=f)
            word: abbabaab
            sage: w = Word('33412321')
            sage: w.lps(l=3)
            word: 12321
            sage: Y = Word
            sage: w = Y('01101001')
            sage: w.lps(l=2)
            word: 1001
            sage: w.lps()
            word: 1001
            sage: w.lps(l=None)
            word: 1001
            sage: Y().lps(l=2)
            Traceback (most recent call last):
            ...
            IndexError: list index out of range
            sage: v = Word('abbabaab')
            sage: pal = v[:0]
            sage: for i in range(1, v.length()+1):
            ....:   pal = v[:i].lps(l=pal.length())
            ....:   pal
            word: a
            word: b
            word: bb
            word: abba
            word: bab
            word: aba
            word: aa
            word: baab
            sage: f = WordMorphism('a->b,b->a')
            sage: v = Word('abbabaab')
            sage: pal = v[:0]
            sage: for i in range(1, v.length()+1):
            ....:   pal = v[:i].lps(f=f, l=pal.length())
            ....:   pal
            word:
            word: ab
            word:
            word: ba
            word: ab
            word: baba
            word: bbabaa
            word: abbabaab
        """
        # If the length of the lps of self[:-1] is not known:
        if l is None:
            l = self.lps_lengths(f)[-1]
            return self[len(self)-l:]

        # If l == w[:-1].length(), there is no shortcut
        if self.length() == l + 1:
            return self.lps(f=f)

        # Obtain the letter to the left (g) and to the right (d) of the
        # precedent lps of self
        g = self[-l-2]
        d = self[-1]

        # If the word g*d is a `f`-palindrome, the result follows
        if f is None:
            if g == d:
                return self[-l-2:]
            else:
                # Otherwise, the length of the lps of self is smallest than l+2
                return self[-l-1:].lps()
        else:
            from sage.combinat.words.morphism import WordMorphism
            f = WordMorphism(f)
            if f(g)[0] == d:
                return self[-l-2:]
            else:
                return self[-l-1:].lps(f=f)

    @cached_method
    def palindromic_lacunas_study(self, f=None):
        r"""
        Return interesting statistics about longest (``f``-)palindromic suffixes
        and lacunas of ``self`` (see [BMBL2008]_ and [BMBFLR2008]_).

        Note that a word `w` has at most `|w| + 1` different palindromic factors
        (see [DJP2001]_). For `f`-palindromes (or pseudopalindromes or theta-palindromes),
        the maximum number of `f`-palindromic factors is `|w|+1-g_f(w)`, where
        `g_f(w)` is the number of pairs `\{a, f(a)\}` such that `a` is a letter,
        `a` is not equal to `f(a)`, and `a` or `f(a)` occurs in `w`, see [Star2011]_.

        INPUT:

        - ``f`` -- involution (default: ``None``) on the alphabet of ``self``.
          It must be callable on letters as well as words (e.g. ``WordMorphism``).
          The default value corresponds to usual palindromes, i.e.,
          ``f`` equal to the identity.

        OUTPUT:

        - ``list`` -- list of the length of the longest palindromic
          suffix (lps) for each non-empty prefix of ``self``
        - ``list`` -- list of all the lacunas, i.e. positions where there is no
          unioccurrent lps
        - ``set`` -- set of palindromic factors of ``self``

        EXAMPLES::

            sage: a,b,c = Word('abbabaabbaab').palindromic_lacunas_study()
            sage: a
            [1, 1, 2, 4, 3, 3, 2, 4, 2, 4, 6, 8]
            sage: b
            [8, 9]
            sage: c          # random order
            set([word: , word: b, word: bab, word: abba, word: bb, word: aa, word: baabbaab, word: baab, word: aba, word: aabbaa, word: a])

        ::

            sage: f = WordMorphism('a->b,b->a')
            sage: a,b,c = Word('abbabaab').palindromic_lacunas_study(f=f)
            sage: a
            [0, 2, 0, 2, 2, 4, 6, 8]
            sage: b
            [0, 2, 4]
            sage: c           # random order
            set([word: , word: ba, word: baba, word: ab, word: bbabaa, word: abbabaab])
            sage: c == set([Word(), Word('ba'), Word('baba'), Word('ab'), Word('bbabaa'), Word('abbabaab')])
            True
        """
        # Initialize the results of computations
        palindromes = set()
        lengths_lps = [None] * self.length()
        lacunas = []

        # Initialize the first lps
        pal = self[:0]
        palindromes.add(pal)

        # For all the non-empty prefixes of self,
        for i in range(self.length()):

            # Compute its longest `f`-palindromic suffix using the preceding lps (pal)
            pal = self[:i+1].lps(l=pal.length(), f=f)

            lengths_lps[i] = pal.length()

            if pal in palindromes:
                lacunas.append(i)
            else:
                palindromes.add(pal)

        return lengths_lps, lacunas, palindromes

    def lacunas(self, f=None):
        r"""
        Return the list of all the lacunas of ``self``.

        A *lacuna* is a position in a word where the longest (`f`-)palindromic
        suffix is not unioccurrent (see [BMBL2008]_).

        INPUT:

        - ``f`` -- involution (default: ``None``) on the alphabet of ``self``. It must
          be callable on letters as well as words (e.g. ``WordMorphism``). The
          default value corresponds to usual palindromes, i.e., ``f`` equal to
          the identity.

        OUTPUT: list of all the lacunas of self

        EXAMPLES::

            sage: w = Word([0,1,1,2,3,4,5,1,13,3])
            sage: w.lacunas()
            [7, 9]
            sage: words.ThueMorseWord()[:100].lacunas()
            [8, 9, 24, 25, 32, 33, 34, 35, 36, 37, 38, 39, 96, 97, 98, 99]
            sage: f = WordMorphism({0:[1],1:[0]})
            sage: words.ThueMorseWord()[:50].lacunas(f)
            [0, 2, 4, 12, 16, 17, 18, 19, 48, 49]
        """
        return self.palindromic_lacunas_study(f=f)[1]

    def lengths_unioccurrent_lps(self, f=None):
        r"""
        Return the list of the lengths of the unioccurrent longest
        (``f``)-palindromic suffixes (lps) for each non-empty prefix of ``self.`` No
        unioccurrent lps are indicated by ``None``.

        It corresponds to the function `H_w` defined in [BMBL2008]_ and [BMBFLR2008]_.

        INPUT:

        - ``f`` -- involution (default: ``None``) on the alphabet of ``self``. It must
          be callable on letters as well as words (e.g. ``WordMorphism``). The
          default value corresponds to usual palindromes, i.e., ``f`` equal to
          the identity.

        OUTPUT:

        a list -- list of the length of the unioccurrent longest palindromic
        suffix (lps) for each non-empty prefix of ``self``.
        No unioccurrent lps are indicated by ``None``.

        EXAMPLES::

            sage: w = Word([0,1,1,2,3,4,5,1,13,3])
            sage: w.lengths_unioccurrent_lps()
            [1, 1, 2, 1, 1, 1, 1, None, 1, None]
            sage: f = words.FibonacciWord()[:20]
            sage: f.lengths_unioccurrent_lps() == f.lps_lengths()[1:]
            True
            sage: t = words.ThueMorseWord()
            sage: t[:20].lengths_unioccurrent_lps()
            [1, 1, 2, 4, 3, 3, 2, 4, None, None, 6, 8, 10, 12, 14, 16, 6, 8, 10, 12]
            sage: f = WordMorphism({1:[0],0:[1]})
            sage: t[:15].lengths_unioccurrent_lps(f)
            [None, 2, None, 2, None, 4, 6, 8, 4, 6, 4, 6, None, 4, 6]
        """
        l = self.lps_lengths(f=f)[1:]
        for i in self.lacunas(f=f):
            l[i] = None
        return l

    def length_maximal_palindrome(self, j, m=None, f=None):
        r"""
        Return the length of the longest palindrome centered at position ``j``.

        INPUT:

        - ``j`` -- rational; position of the symmetry axis of the palindrome.
          Must return an integer when doubled. It is an integer when the
          center of the palindrome is a letter.

        - ``m`` -- integer (default: ``None``); minimal length of palindrome, if known.
          The parity of ``m`` can't be the same as the parity of ``2j``.

        - ``f`` -- involution (default: ``None``) on the alphabet; it must be
          callable on letters as well as words (e.g. ``WordMorphism``)

        OUTPUT: length of the longest ``f``-palindrome centered at position ``j``

        EXAMPLES::

            sage: Word('01001010').length_maximal_palindrome(3/2)
            0
            sage: Word('01101001').length_maximal_palindrome(3/2)
            4
            sage: Word('01010').length_maximal_palindrome(j=3, f='0->1,1->0')
            0
            sage: Word('01010').length_maximal_palindrome(j=2.5, f='0->1,1->0')
            4
            sage: Word('0222220').length_maximal_palindrome(3, f='0->1,1->0,2->2')
            5

        ::

            sage: w = Word('abcdcbaxyzzyx')
            sage: w.length_maximal_palindrome(3)
            7
            sage: w.length_maximal_palindrome(3, 3)
            7
            sage: w.length_maximal_palindrome(3.5)
            0
            sage: w.length_maximal_palindrome(9.5)
            6
            sage: w.length_maximal_palindrome(9.5, 2)
            6

        TESTS:

        These are wrong inputs::

            sage: w.length_maximal_palindrome(9.6)
            Traceback (most recent call last):
            ...
            ValueError: j must be positive, inferior to length of self
            sage: w.length_maximal_palindrome(3, 2)
            Traceback (most recent call last):
            ...
            ValueError: (2*j-m-1)/2(=3/2) must be an integer, i.e., 2*j(=6) and
            m(=2) can't have the same parity
            sage: w.length_maximal_palindrome(9.5, 3)
            Traceback (most recent call last):
            ...
            ValueError: (2*j-m-1)/2(=15/2) must be an integer, i.e., 2*j(=19) and
            m(=3) can't have the same parity
        """
        # Ensure `f` is an involutory word morphism
        if f is not None:
            from sage.combinat.words.morphism import WordMorphism
            if not isinstance(f, WordMorphism):
                f = WordMorphism(f)
            if not f.is_involution():
                raise ValueError("f must be an involution")

        # Ensure j is a valid entry
        jj = 2*j
        if not jj.is_integer() or j < 0 or j >= len(self):
            raise ValueError("j must be positive, inferior to length of self")
        jj = Integer(jj)

        # Initialize length of the known palindrome
        if m is None:
            m = 0 if jj % 2 else -1

        # Initialize the next (left) position to check
        i = (jj - m - 1) / 2
        if not i.is_integer():
            raise ValueError(f"(2*j-m-1)/2(={i}) must be an integer, i.e., "
                             f"2*j(={jj}) and m(={m}) can't "
                             "have the same parity")
        i = Integer(i)

        # Compute
        if f is None:
            while i >= 0 and jj-i < len(self) and self[i] == self[jj-i]:
                i -= 1
        else:
            while i >= 0 and jj-i < len(self) and self[i] == f(self[jj-i])[0]:
                i -= 1
        if jj == 2 * i:
            return 0
        return jj - 2 * i - 1

    def lengths_maximal_palindromes(self, f=None):
        r"""
        Return the length of maximal palindromes centered at each position.

        INPUT:

        - ``f`` -- involution (default: ``None``) on the alphabet of ``self``;
          it must be callable on letters as well as words (e.g.
          ``WordMorphism``)

        OUTPUT: list; the length of the maximal palindrome (or ``f``-palindrome)
        with a given symmetry axis (letter or space between two letters)

        EXAMPLES::

            sage: Word('01101001').lengths_maximal_palindromes()
            [0, 1, 0, 1, 4, 1, 0, 3, 0, 3, 0, 1, 4, 1, 0, 1, 0]
            sage: Word('00000').lengths_maximal_palindromes()
            [0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0]
            sage: Word('0').lengths_maximal_palindromes()
            [0, 1, 0]
            sage: Word('').lengths_maximal_palindromes()
            [0]
            sage: Word().lengths_maximal_palindromes()
            [0]
            sage: f = WordMorphism('a->b,b->a')
            sage: Word('abbabaab').lengths_maximal_palindromes(f)
            [0, 0, 2, 0, 0, 0, 2, 0, 8, 0, 2, 0, 0, 0, 2, 0, 0]
        """
        if f is not None:
            from sage.combinat.words.morphism import WordMorphism
            if not isinstance(f, WordMorphism):
                f = WordMorphism(f)
            if not f.is_involution():
                raise ValueError("f must be an involution")

        LPC = []  # lengths of the maximal palindromes centered at a position
        LPC.append(0)

        k = 0  # index, center of rightmost-ending `f`-palindrome encountered

        for j in range(1, 2 * len(self) + 1):
            if j >= k + LPC[k]:
                p = self.length_maximal_palindrome((j - 1)*0.5, -(j % 2), f)
                LPC.append(p)
                if j + p > k + LPC[k]:
                    k = j

            # If the center is included in an encountered `f`-palindrome
            else:
                # If the `f`-palindrome centered at position j is not the
                # longest proper `f`-palindromic suffix of the maximal
                # `f`-palindrome centered at k
                if LPC[k] + k - j != LPC[2*k - j]:
                    LPC.append(min(LPC[k] + k - j, LPC[2*k - j]))

                else:
                    mp = LPC[k] + k - j
                    p = self.length_maximal_palindrome((j-1)*0.5, mp, f)
                    LPC.append(p)
                    k = j
        return LPC

    def lps_lengths(self, f=None):
        r"""
        Return the length of the longest palindromic suffix of each prefix.

        INPUT:

        - ``f`` -- involution (default: ``None``) on the alphabet of ``self``.
          It must be callable on letters as well as words (e.g.
          ``WordMorphism``).

        OUTPUT: list; the length of the longest palindromic (or
        ``f``-palindromic) suffix of each prefix of ``self``

        EXAMPLES::

            sage: Word('01101001').lps_lengths()
            [0, 1, 1, 2, 4, 3, 3, 2, 4]
            sage: Word('00000').lps_lengths()
            [0, 1, 2, 3, 4, 5]
            sage: Word('0').lps_lengths()
            [0, 1]
            sage: Word('').lps_lengths()
            [0]
            sage: Word().lps_lengths()
            [0]
            sage: f = WordMorphism('a->b,b->a')
            sage: Word('abbabaab').lps_lengths(f)
            [0, 0, 2, 0, 2, 2, 4, 6, 8]
        """
        LPC = self.lengths_maximal_palindromes(f)
        Nk = LPC[0]
        LPS = [0]  # lengths of the longest palindromic suffix of prefixes

        for j in range(1, 2 * len(self) + 1):
            Nj = j + LPC[j]
            if Nj > Nk:
                LPS.extend(i - j for i in range(Nk + 2 - (Nk % 2), Nj + 1, 2))
                Nk = Nj
        return LPS

    def palindromes(self, f=None):
        r"""
        Return the set of all palindromic (or ``f``-palindromic) factors of ``self``.

        INPUT:

        - ``f`` -- involution (default: ``None``) on the alphabet of ``self``;
          it must be callable on letters as well as words (e.g. ``WordMorphism``).

        OUTPUT: a set -- If ``f`` is ``None``, the set of all palindromic
        factors of ``self``; otherwise, the set of all ``f``-palindromic
        factors of ``self``

        EXAMPLES::

            sage: sorted(Word('01101001').palindromes())
            [word: , word: 0, word: 00, word: 010, word: 0110, word: 1, word: 1001, word: 101, word: 11]
            sage: sorted(Word('00000').palindromes())
            [word: , word: 0, word: 00, word: 000, word: 0000, word: 00000]
            sage: sorted(Word('0').palindromes())
            [word: , word: 0]
            sage: sorted(Word('').palindromes())
            [word: ]
            sage: sorted(Word().palindromes())
            [word: ]
            sage: f = WordMorphism('a->b,b->a')
            sage: sorted(Word('abbabaab').palindromes(f))
            [word: , word: ab, word: abbabaab, word: ba, word: baba, word: bbabaa]
        """
        LPS = self.lps_lengths(f)
        return {self[i - LPS[i]: i] for i in range(len(self) + 1)}

    def palindromic_complexity(self, n):
        r"""
        Return the number of distinct palindromic factors of length ``n`` of ``self``.

        INPUT:

        - ``n`` -- the length of the factors

        EXAMPLES::

            sage: w = words.FibonacciWord()[:100]
            sage: [w.palindromic_complexity(i) for i in range(20)]
            [1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2]

        ::

            sage: w = words.ThueMorseWord()[:1000]
            sage: [w.palindromic_complexity(i) for i in range(20)]
            [1, 2, 2, 2, 2, 0, 4, 0, 4, 0, 4, 0, 4, 0, 2, 0, 2, 0, 4, 0]
        """
        return len([1 for x in self.palindromes() if len(x) == n])

    def palindrome_prefixes(self):
        r"""
        Return a list of all palindrome prefixes of ``self``.

        EXAMPLES::

            sage: w = Word('abaaba')
            sage: w.palindrome_prefixes()
            [word: , word: a, word: aba, word: abaaba]
            sage: w = Word('abbbbbbbbbb')
            sage: w.palindrome_prefixes()
            [word: , word: a]
        """
        return list(self.palindrome_prefixes_iterator())

    def defect(self, f=None):
        r"""
        Return the defect of ``self``.

        The *defect* of a finite word `w` is given by the difference between
        the maximum number of possible palindromic factors in a word of length
        `|w|` and the actual number of palindromic factors contained in `w`.
        It is well known that the maximum number of palindromic factors in `w`
        is `|w|+1` (see [DJP2001]_).

        An optional involution on letters ``f`` can be given. In that case, the
        *f-palindromic defect* (or *pseudopalindromic defect*, or
        *theta-palindromic defect*) of `w` is returned. It is a
        generalization of defect to f-palindromes. More precisely, the defect is
        `D(w)=|w|+1-g_f(w)-|PAL_f(w)|`, where `PAL_f(w)` denotes the set of
        f-palindromic factors of `w` (including the empty word) and `g_f(w)` is
        the number of pairs `\{a, f(a)\}` such that `a` is a letter, `a` is not
        equal to `f(a)`, and `a` or `f(a)` occurs in `w`. In the case of usual
        palindromes (i.e., for ``f`` not given or equal to the identity),
        `g_f(w) = 0` for all `w`. See [BHNR2004]_ for usual palindromes and [Star2011]_
        for f-palindromes.

        INPUT:

        - ``f`` -- involution (default: ``None``) on the alphabet of ``self``. It must
          be callable on letters as well as words (e.g. ``WordMorphism``). The
          default value corresponds to usual palindromes, i.e., ``f`` equal to
          the identity.

        OUTPUT:

        an integer -- If ``f`` is ``None``, the palindromic defect of ``self``;
        otherwise, the ``f``-palindromic defect of ``self``.

        EXAMPLES::

            sage: Word('ara').defect()
            0
            sage: Word('abcacba').defect()
            1

        It is known that Sturmian words (see [DJP2001]_) have zero defect::

            sage: words.FibonacciWord()[:100].defect()
            0

            sage: sa = WordMorphism('a->ab,b->b')
            sage: sb = WordMorphism('a->a,b->ba')
            sage: w = (sa*sb*sb*sa*sa*sa*sb).fixed_point('a')
            sage: w[:30].defect()                                                       # needs sage.modules
            0
            sage: w[110:140].defect()                                                   # needs sage.modules
            0

        It is even conjectured that the defect of an aperiodic word which is
        a fixed point of a primitive morphism is either `0` or infinite
        (see [BBGL2008]_)::

            sage: w = words.ThueMorseWord()
            sage: w[:50].defect()                                                       # needs sage.modules
            12
            sage: w[:100].defect()                                                      # needs sage.modules
            16
            sage: w[:300].defect()                                                      # needs sage.modules
            52

        For generalized defect with an involution different from the identity,
        there is always a letter which is not a palindrome! This is the reason
        for the modification of the definition::

            sage: f = WordMorphism('a->b,b->a')
            sage: Word('a').defect(f)
            0
            sage: Word('ab').defect(f)
            0
            sage: Word('aa').defect(f)
            1
            sage: Word('abbabaabbaababba').defect(f)
            3

        ::

            sage: f = WordMorphism('a->b,b->a,c->c')
            sage: Word('cabc').defect(f)
            0
            sage: Word('abcaab').defect(f)
            2

        Other examples::

            sage: Word('000000000000').defect()
            0
            sage: Word('011010011001').defect()
            2
            sage: Word('0101001010001').defect()
            0
            sage: Word().defect()
            0
            sage: Word('abbabaabbaababba').defect()
            2
        """
        g_w = 0
        if f is not None:
            from sage.combinat.words.morphism import WordMorphism
            if not isinstance(f, WordMorphism):
                f = WordMorphism(f)
            if not f.is_involution():
                raise ValueError("f must be an involution")
            D = f.domain()
            A = set(map(D, self.letters()))
            while A:
                x = A.pop()
                if f(x) != x:  # count only non f-palindromic letters
                    if f(x) in A:
                        A.remove(f(x))
                    g_w += 1

        return self.length()+1-g_w-len(self.palindromes(f=f))

    def is_full(self, f=None):
        r"""
        Return ``True`` if ``self`` has defect `0`, and ``False`` otherwise.

        A word is *full* (or *rich*) if its defect is zero (see [BHNR2004]_).

        If ``f`` is given, then the ``f``-palindromic defect is used (see [PeSt2011]_).

        INPUT:

        - ``f`` -- involution (default: ``None``) on the alphabet of ``self``;
          it must be callable on letters as well as words (e.g. ``WordMorphism``)

        OUTPUT:

        boolean -- If ``f`` is ``None``, whether ``self`` is full;
        otherwise, whether ``self`` is full of ``f``-palindromes.

        EXAMPLES::

            sage: words.ThueMorseWord()[:100].is_full()
            False
            sage: words.FibonacciWord()[:100].is_full()
            True
            sage: Word('000000000000000').is_full()
            True
            sage: Word('011010011001').is_full()
            False
            sage: Word('2194').is_full()
            True
            sage: Word().is_full()
            True

        ::

            sage: f = WordMorphism('a->b,b->a')
            sage: Word().is_full(f)
            True
            sage: w = Word('ab')
            sage: w.is_full()
            True
            sage: w.is_full(f)
            True

        ::

            sage: f = WordMorphism('a->b,b->a')
            sage: Word('abab').is_full(f)
            True
            sage: Word('abba').is_full(f)
            False

        A simple example of an infinite word full of f-palindromes::

            sage: p = WordMorphism({0:'abc',1:'ab'})
            sage: f = WordMorphism('a->b,b->a,c->c')
            sage: p(words.FibonacciWord()[:50]).is_full(f)
            True
            sage: p(words.FibonacciWord()[:150]).is_full(f)
            True
        """
        return self.defect(f=f) == 0

    is_rich = is_full

    def palindromic_closure(self, side='right', f=None):
        r"""
        Return the shortest palindrome having ``self`` as a prefix
        (or as a suffix if ``side`` is ``'left'``).

        See [DeLuca2006]_.

        INPUT:

        - ``side`` -- ``'right'`` or ``'left'`` (default: ``'right'``) the
          direction of the  closure

        - ``f`` -- involution (default: ``None``) on the alphabet of ``self``;
          it must be callable on letters as well as words (e.g. ``WordMorphism``)

        OUTPUT:

        a word -- If ``f`` is ``None``, the right palindromic closure of ``self``;
        otherwise, the right ``f``-palindromic closure of ``self``.
        If ``side`` is ``'left'``, the left palindromic closure.

        EXAMPLES::

            sage: Word('1233').palindromic_closure()
            word: 123321
            sage: Word('12332').palindromic_closure()
            word: 123321
            sage: Word('0110343').palindromic_closure()
            word: 01103430110
            sage: Word('0110343').palindromic_closure(side='left')
            word: 3430110343
            sage: Word('01105678').palindromic_closure(side='left')
            word: 876501105678
            sage: w = Word('abbaba')
            sage: w.palindromic_closure()
            word: abbababba

        ::

            sage: f = WordMorphism('a->b,b->a')
            sage: w.palindromic_closure(f=f)
            word: abbabaab
            sage: w.palindromic_closure(f=f, side='left')
            word: babaabbaba

        TESTS::

            sage: f = WordMorphism('a->c,c->a')
            sage: w.palindromic_closure(f=f, side='left')
            Traceback (most recent call last):
            ...
            ValueError: b not in alphabet
        """
        if f is None:
            if side == 'right':
                l = self.lps().length()
                # return self * self[-(l+1)::-1]
                return self * self[:self.length() - l].reversal()
            elif side == 'left':
                l = self.reversal().lps().length()
                return self[:l-1:-1] * self
            else:
                raise ValueError("side must be either 'left' or 'right' (not %s) " % side)
        else:
            from sage.combinat.words.morphism import WordMorphism
            f = WordMorphism(f)
            if not f.is_involution():
                raise ValueError("f must be an involution")
            if side == 'right':
                l = self.lps(f=f).length()
                return self * f(self[-(l+1)::-1])
            elif side == 'left':
                l = self.reversal().lps(f=f).length()
                return f(self[:l-1:-1]) * self
            else:
                raise ValueError("side must be either 'left' or 'right' (not %s) " % side)

    def is_symmetric(self, f=None):
        r"""
        Return ``True`` if ``self`` is symmetric (or ``f``-symmetric), and
        ``False`` otherwise.

        A word is *symmetric* (resp. `f`-*symmetric*) if it is the
        product of two palindromes (resp. `f`-palindromes).
        See [BHNR2004]_ and [DeLuca2006]_.

        INPUT:

        - ``f`` -- involution (default: ``None``) on the alphabet of ``self``;
          it must be callable on letters as well as words (e.g. ``WordMorphism``)

        EXAMPLES::

            sage: Word('abbabab').is_symmetric()
            True
            sage: Word('ababa').is_symmetric()
            True
            sage: Word('aababaabba').is_symmetric()
            False
            sage: Word('aabbbaababba').is_symmetric()
            False
            sage: f = WordMorphism('a->b,b->a')
            sage: Word('aabbbaababba').is_symmetric(f)
            True
        """
        square = self * self
        return square.lps_lengths(f)[-1] >= self.length()

    def length_border(self):
        r"""
        Return the length of the border of ``self``.

        The *border* of a word is the longest word that is both a proper
        prefix and a proper suffix of ``self``.

        EXAMPLES::

            sage: Word('121').length_border()
            1
            sage: Word('1').length_border()
            0
            sage: Word('1212').length_border()
            2
            sage: Word('111').length_border()
            2
            sage: Word().length_border() is None
            True
        """
        if self.is_empty():
            return None
        return self.prefix_function_table()[-1]

    def border(self):
        r"""
        Return the longest word that is both a proper prefix and a proper
        suffix of ``self``.

        EXAMPLES::

            sage: Word('121212').border()
            word: 1212
            sage: Word('12321').border()
            word: 1
            sage: Word().border() is None
            True
        """
        if self.is_empty():
            return None
        return self[:self.length_border()]

    def minimal_period(self):
        r"""
        Return the period of ``self``.

        Let `A` be an alphabet. An integer `p\geq 1` is a *period* of a
        word `w=a_1a_2\cdots a_n` where `a_i\in A` if `a_i=a_{i+p}` for
        `i=1,\ldots,n-p`. The smallest period of `w` is called *the*
        period of `w`. See Chapter 1 of [Lot2002]_.

        EXAMPLES::

            sage: Word('aba').minimal_period()
            2
            sage: Word('abab').minimal_period()
            2
            sage: Word('ababa').minimal_period()
            2
            sage: Word('ababaa').minimal_period()
            5
            sage: Word('ababac').minimal_period()
            6
            sage: Word('aaaaaa').minimal_period()
            1
            sage: Word('a').minimal_period()
            1
            sage: Word().minimal_period()
            1
        """
        if self.is_empty():
            return 1
        return self.length() - self.length_border()

    def order(self):
        r"""
        Return the order of ``self``.

        Let `p(w)` be the period of a word `w`. The positive rational number
        `|w|/p(w)` is the *order* of `w`. See Chapter 8 of [Lot2002]_.

        OUTPUT: rational; the order

        EXAMPLES::

            sage: Word('abaaba').order()
            2
            sage: Word('ababaaba').order()
            8/5
            sage: Word('a').order()
            1
            sage: Word('aa').order()
            2
            sage: Word().order()
            0
        """
        from sage.rings.rational import Rational
        return Rational((self.length(), self.minimal_period()))

    def critical_exponent(self):
        r"""
        Return the critical exponent of ``self``.

        The *critical exponent* of a word is the supremum of the order of
        all its (finite) factors. See [Dej1972]_.

        .. NOTE::

            The implementation here uses the suffix tree to enumerate all the
            factors. It should be improved (especially when the critical
            exponent is larger than 2).

        EXAMPLES::

            sage: Word('aaba').critical_exponent()
            2
            sage: Word('aabaa').critical_exponent()
            2
            sage: Word('aabaaba').critical_exponent()
            7/3
            sage: Word('ab').critical_exponent()
            1
            sage: Word('aba').critical_exponent()
            3/2
            sage: words.ThueMorseWord()[:20].critical_exponent()
            2

        For the Fibonacci word, the critical exponent is known to be
        `(5+\sqrt(5))/2`. With a prefix of length 500, we obtain a lower bound::

            sage: words.FibonacciWord()[:500].critical_exponent()
            320/89

        It is an error to compute the critical exponent of the empty word::

            sage: Word('').critical_exponent()
            Traceback (most recent call last):
            ...
            ValueError: no critical exponent for empty word
        """
        if not self:
            raise ValueError("no critical exponent for empty word")
        else:
            st = self.suffix_tree()
            pft = [0] * self.length()  # the prefix function table
            queue = [(0, 0, -1, 0)]    # suffix tree vertices to visit for Depth First Search
            best_exp = 1               # best exponent so far
            while queue:
                v, i, j, l = queue.pop()
                for k in range(i, j+1):
                    if l-j+k-1 != 0:
                        m = pft[l-j+k-2]
                        while m > 0 and self[j-l+m] != self[k-1]:
                            m = pft[m-1]
                        if self[j-l+m] == self[k-1]:
                            m += 1
                    else:
                        m = 0
                    current_pos = k-j+l-1
                    pft[current_pos] = m
                    current_exp = QQ((current_pos+1, current_pos+1-m))
                    if current_exp > best_exp:
                        best_exp = current_exp
                for ((i, j), u) in st._transition_function[v].items():
                    if j is None:
                        j = self.length()
                    queue.append((u, i, j, l+j-i+1))
            return best_exp

    def is_overlap(self):
        r"""
        Return ``True`` if ``self`` is an overlap, and ``False`` otherwise.

        EXAMPLES::

            sage: Word('12121').is_overlap()
            True
            sage: Word('123').is_overlap()
            False
            sage: Word('1231').is_overlap()
            False
            sage: Word('123123').is_overlap()
            False
            sage: Word('1231231').is_overlap()
            True
            sage: Word().is_overlap()
            False
        """
        if self.length() == 0:
            return False
        return self.length_border() > self.length()//2

    def primitive_length(self):
        r"""
        Return the length of the primitive of ``self``.

        EXAMPLES::

            sage: Word('1231').primitive_length()
            4
            sage: Word('121212').primitive_length()
            2
        """
        l = self.length()
        if l == 0:
            return 0
        p = self.minimal_period()
        return p if l % p == 0 else l

    def is_primitive(self):
        r"""
        Return ``True`` if ``self`` is primitive, and ``False`` otherwise.

        A finite word `w` is *primitive* if it is not a positive integer
        power of a shorter word.

        EXAMPLES::

            sage: Word('1231').is_primitive()
            True
            sage: Word('111').is_primitive()
            False
        """
        return self.length() == self.primitive_length()

    def primitive(self):
        r"""
        Return the primitive of ``self``.

        EXAMPLES::

            sage: Word('12312').primitive()
            word: 12312
            sage: Word('121212').primitive()
            word: 12
        """
        return self[:self.primitive_length()]

    def exponent(self):
        r"""
        Return the exponent of ``self``.

        OUTPUT: integer; the exponent

        EXAMPLES::

            sage: Word('1231').exponent()
            1
            sage: Word('121212').exponent()
            3
            sage: Word().exponent()
            0
        """
        if self.length() == 0:
            return 0
        return self.length() // self.primitive_length()

    def has_period(self, p):
        r"""
        Return ``True`` if ``self`` has the period `p`,
        ``False`` otherwise.

        .. NOTE::

            By convention, integers greater than the length
            of ``self`` are periods of ``self``.

        INPUT:

        - ``p`` -- integer to check if it is a period
          of ``self``

        EXAMPLES::

            sage: w = Word('ababa')
            sage: w.has_period(2)
            True
            sage: w.has_period(3)
            False
            sage: w.has_period(4)
            True
            sage: w.has_period(-1)
            False
            sage: w.has_period(5)
            True
            sage: w.has_period(6)
            True
        """
        if p < 0:
            return False
        elif p >= len(self):
            return True
        else:
            for i in range(len(self) - p):
                if self[i] != self[i + p]:
                    return False
            return True

    def periods(self, divide_length=False):
        r"""
        Return a list containing the periods of ``self``
        between `1` and `n - 1`, where `n` is the length
        of ``self``.

        INPUT:

        - ``divide_length`` -- boolean (default: ``False``);
          when set to ``True``, then only periods that divide
          the length of ``self`` are considered

        OUTPUT: list of positive integers

        EXAMPLES::

            sage: w = Word('ababab')
            sage: w.periods()
            [2, 4]
            sage: w.periods(divide_length=True)
            [2]
            sage: w = Word('ababa')
            sage: w.periods()
            [2, 4]
            sage: w.periods(divide_length=True)
            []
        """
        n = len(self)
        if divide_length:
            possible = (i for i in range(1, n) if not n % i)
        else:
            possible = range(1, n)
        return [x for x in possible if self.has_period(x)]

    def longest_common_subword(self, other):
        r"""
        Return a longest subword of ``self`` and ``other``.

        A subword of a word is a subset of the word's letters, read in the
        order in which they appear in the word.

        For more information, see
        :wikipedia:`Longest_common_subsequence_problem`.

        INPUT:

        - ``other`` -- a word

        ALGORITHM:

        For any indices `i,j`, we compute the longest common subword ``lcs[i,j]`` of
        ``self[:i]`` and ``other[:j]``. This can be easily obtained as the longest
        of

        - ``lcs[i-1,j]``

        - ``lcs[i,j-1]``

        - ``lcs[i-1,j-1]+self[i]`` if ``self[i]==other[j]``

        EXAMPLES::

            sage: v1 = Word("abc")
            sage: v2 = Word("ace")
            sage: v1.longest_common_subword(v2)
            word: ac

            sage: w1 = Word("1010101010101010101010101010101010101010")
            sage: w2 = Word("0011001100110011001100110011001100110011")
            sage: w1.longest_common_subword(w2)
            word: 00110011001100110011010101010

        TESTS::

            sage: Word().longest_common_subword(Word())
            word:

        .. SEEALSO::

            :meth:`is_subword_of`
        """
        from sage.combinat.words.word import Word
        if len(self) == 0 or len(other) == 0:
            return Word()

        w2 = list(other)

        # In order to avoid storing lcs[i,j] for each pair i,j of indices, we
        # only store the lcs[i,j] for two consecutive values of i. At any step
        # of the algorithm, lcs[i,j] is stored at lcs[0][j] and lcs[i-1,j] is
        # stored at lcs[1][j]

        # The weird +1 that follows exists to make sure that lcs[i,-1] returns
        # the empty word.
        lcs = [[[] for _ in repeat(None, len(w2) + 1)] for j in range(2)]

        for i, l1 in enumerate(self):
            for j, l2 in enumerate(other):
                lcs[0][j] = max(lcs[0][j-1], lcs[1][j],
                                lcs[1][j-1] + ([l1] if l1 == l2 else []), key=len)

            # Maintaining the meaning of lcs for the next loop
            lcs.pop(1)
            lcs.insert(0, [[] for _ in repeat(None, len(w2) + 1)])

        return Word(lcs[1][-2])

    def is_subword_of(self, other):
        r"""
        Return ``True`` if ``self`` is a subword of ``other``, and ``False`` otherwise.

        A finite word `u` is a *subword* of a finite word `v` if `u` is a
        subsequence of `v`. See Chapter 6 on Subwords in [Lot1997]_.

        Some references define subword as a consecutive subsequence. Use
        :meth:`is_factor` if this is what you need.

        INPUT:

        - ``other`` -- a finite word

        EXAMPLES::

            sage: Word('bb').is_subword_of(Word('ababa'))
            True
            sage: Word('bbb').is_subword_of(Word('ababa'))
            False

        ::

            sage: Word().is_subword_of(Word('123'))
            True
            sage: Word('123').is_subword_of(Word('3211333213233321'))
            True
            sage: Word('321').is_subword_of(Word('11122212112122133111222332'))
            False

        .. SEEALSO::

            :meth:`longest_common_subword`
            :meth:`number_of_subword_occurrences`
            :meth:`is_factor`
        """
        its = iter(self)
        try:
            s = next(its)
            for e in other:
                if s == e:
                    s = next(its)
            else:
                return False
        except StopIteration:
            return True

    def subword_complementaries(self, other):
        """
        Return the possible complementaries ``other`` minus ``self`` if
        ``self`` is a subword of ``other`` (empty list otherwise).
        The complementary is made of all the letters that are in ``other`` once
        we removed the letters of ``self``.
        There can be more than one.

        To check whether ``self`` is a subword of ``other`` (without knowing its
        complementaries), use ``self.is_subword_of(other)``, and to count the
        number of occurrences of ``self`` in ``other``, use
        ``other.number_of_subword_occurrences(self)``.

        INPUT:

        - ``other`` -- finite word

        OUTPUT: list of all the complementary subwords of ``self`` in ``other``

        EXAMPLES::

            sage: Word('tamtam').subword_complementaries(Word('ta'))
            []

            sage: Word('mta').subword_complementaries(Word('tamtam'))
            [word: tam]

            sage: Word('ta').subword_complementaries(Word('tamtam'))
            [word: mtam, word: amtm, word: tamm]

            sage: Word('a').subword_complementaries(Word('a'))
            [word: ]
        """

        ls = self.length()
        lo = other.length()

        # Create a matrix to declare when letters in ``self`` and ``other`` are
        # equal or not
        Eq = [[self[i] == other[j] for j in range(lo)] for i in range(ls)]

        # Create a matrix that tells the positions of subwords of the suffixes
        Mpos = [[[] for _ in repeat(None, lo)] for i in range(ls)]
        for j in range(lo):
            if Eq[ls-1][j]:
                Mpos[ls-1][j] = [[j]]
        for i in range(ls-2, -1, -1):
            for j in range(lo):
                if Eq[i][j]:
                    temp = []
                    for k in range(j+1, lo):
                        if Eq[i+1][k]:
                            m = Mpos[i+1][k]
                            if len(m) == 1:
                                temp.append([j]+m[0])
                            if len(m) > 1:
                                temp.extend([j] + sw for sw in m)
                    Mpos[i][j] = temp

        # Create the list of positions for occurrences of `self` as a subword
        selfpos = []
        for j in range(lo):
            selfpos.extend(Mpos[0][j])

        # Create the list of the complementaries of `self`
        from sage.combinat.words.word import Word
        comp_words = []
        for sp in selfpos:  # list with positions of one occurrence of `self`
            comp_pos = (i for i in range(lo) if i not in set(sp))
            comp_words.append(Word([other[i] for i in comp_pos]))
        return comp_words

    def is_lyndon(self) -> bool:
        r"""
        Return ``True`` if ``self`` is a Lyndon word, and ``False``
        otherwise.

        A *Lyndon word* is a non-empty word that is lexicographically
        smaller than each of its proper suffixes (for the given order
        on its alphabet). That is, `w` is a Lyndon word if `w` is non-empty
        and for each factorization `w = uv` (with `u`, `v` both non-empty),
        we have `w < v`.

        Equivalently, `w` is a Lyndon word iff `w` is a non-empty word that is
        lexicographically smaller than each of its proper conjugates for the
        given order on its alphabet.

        See for instance [Lot1983]_.

        EXAMPLES::

            sage: Word('123132133').is_lyndon()
            True
            sage: Word().is_lyndon()
            False
            sage: Word('122112').is_lyndon()
            False

        TESTS:

        A sanity check: ``LyndonWords`` generates Lyndon words, so we
        filter all words of length `n<10` on the alphabet [1,2,3] for
        Lyndon words, and compare with the ``LyndonWords`` generator::

            sage: for n in range(1,10):
            ....:     lw1 = [w for w in Words([1,2,3], n) if w.is_lyndon()]
            ....:     lw2 = LyndonWords(3,n)
            ....:     if set(lw1) != set(lw2): print(False)

        Filter all words of length 8 on the alphabet [c,a,b] for Lyndon
        words, and compare with the :class:`LyndonWords` generator after
        mapping [a,b,c] to [2,3,1]::

            sage: lw = [w for w in Words('cab', 8) if w.is_lyndon()]
            sage: phi = WordMorphism({'a':2,'b':3,'c':1})
            sage: set(map(phi, lw)) == set(LyndonWords(3,8))
            True
        """
        if self.is_empty():
            return False
        key = self.parent().sortkey_letters
        n = self.length()
        i, j = 0, 1
        while j < n:
            ki = key(self[i])
            kj = key(self[j])
            if ki == kj:
                # increment i and j
                i += 1
                j += 1
            elif ki < kj:
                # reset i, increment j
                i = 0
                j += 1
            else:
                # we found the first word in the lyndon factorization;
                return False
        else:
            return i == 0

    def lyndon_factorization(self):
        r"""
        Return the Lyndon factorization of ``self``.

        The *Lyndon factorization* of a finite word `w` is the unique
        factorization of `w` as a non-increasing product of Lyndon words,
        i.e., `w = l_1\cdots l_n` where each `l_i` is a Lyndon word and
        `l_1\geq \cdots \geq l_n`. See for instance [Duv1983]_.

        OUTPUT: the list `[l_1, \ldots, l_n]` of factors obtained

        EXAMPLES::

            sage: Word('010010010001000').lyndon_factorization()
            (01, 001, 001, 0001, 0, 0, 0)
            sage: Words('10')('010010010001000').lyndon_factorization()
            (0, 10010010001000)
            sage: Word('abbababbaababba').lyndon_factorization()
            (abb, ababb, aababb, a)
            sage: Words('ba')('abbababbaababba').lyndon_factorization()
            (a, bbababbaaba, bba)
            sage: Word([1,2,1,3,1,2,1]).lyndon_factorization()
            (1213, 12, 1)

        TESTS::

            sage: Words('01')('').lyndon_factorization()
            ()
            sage: Word('01').lyndon_factorization()
            (01)
            sage: Words('10')('01').lyndon_factorization()
            (0, 1)
            sage: lynfac = Word('abbababbaababba').lyndon_factorization()
            sage: [x.is_lyndon() for x in lynfac]
            [True, True, True, True]
            sage: lynfac = Words('ba')('abbababbaababba').lyndon_factorization()
            sage: [x.is_lyndon() for x in lynfac]
            [True, True, True]
            sage: w = words.ThueMorseWord()[:1000]
            sage: w == prod(w.lyndon_factorization())
            True

        See [Me1997]_.
        """
        key = self.parent().sortkey_letters
        # We compute the indexes of the factorization.
        n = self.length()
        k = -1
        F = [0]
        while k < n-1:
            i = k+1
            j = k+2
            while j < n:
                ki = key(self[i])
                kj = key(self[j])
                if ki < kj:
                    i = k+1
                    j += 1
                elif ki == kj:
                    i += 1
                    j += 1
                else:
                    break
            while k < i:
                F.append(k + j - i + 1)
                k = k + j - i
        return Factorization([self[F[l]:F[l+1]] for l in range(len(F)-1)])

    def inversions(self):
        r"""
        Return a list of the inversions of ``self``. An inversion is a pair
        `(i,j)` of nonnegative integers `i < j` such that ``self[i] > self[j]``.

        EXAMPLES::

            sage: Word([1,2,3,2,2,1]).inversions()
            [[1, 5], [2, 3], [2, 4], [2, 5], [3, 5], [4, 5]]
            sage: Words([3,2,1])([1,2,3,2,2,1]).inversions()
            [[0, 1], [0, 2], [0, 3], [0, 4], [1, 2]]
            sage: Word('abbaba').inversions()
            [[1, 3], [1, 5], [2, 3], [2, 5], [4, 5]]
            sage: Words('ba')('abbaba').inversions()
            [[0, 1], [0, 2], [0, 4], [3, 4]]
        """
        inversion_list = []
        cmp_key = self._parent.sortkey_letters
        for i1, letter1 in enumerate(self):
            k1 = cmp_key(letter1)
            for i2, letter2 in enumerate(self[i1 + 1:]):
                k2 = cmp_key(letter2)
                if k1 > k2:
                    inversion_list.append([i1, i1 + i2 + 1])
        return inversion_list

    # TODO: This function should be defined for words of integers, but it
    # naturally is defined over an alphabet with a rank function....
    def degree(self, weights=None):
        r"""
        Return the weighted degree of ``self``, where the weighted degree of
        each letter in the ordered alphabet is given by ``weights``, which
        defaults to ``[1, 2, 3, ...]``.

        INPUT:

        - ``weights`` -- list or tuple, or dictionary keyed by the
          letters occurring in ``self``

        EXAMPLES::

            sage: Word([1,2,3]).degree()
            6
            sage: Word([3,2,1]).degree()
            6
            sage: Words("ab")("abba").degree()
            6
            sage: Words("ab")("abba").degree([0,2])
            4
            sage: Words("ab")("abba").degree([-1,-1])
            -4
            sage: Words("ab")("aabba").degree([1,1])
            5
            sage: Words([1,2,4])([1,2,4]).degree()
            6
            sage: Word([1,2,4]).degree()
            7
            sage: Word("aabba").degree({'a':1,'b':2})
            7
            sage: Word([0,1,0]).degree({0:17,1:0})
            34
        """
        if isinstance(weights, dict):
            deg = 0
            for a in self:
                deg += weights[a]
            return deg

        if hasattr(self._parent._alphabet, "rank"):
            rank_fcn = self._parent._alphabet.rank
            deg = 0
            if weights is None:
                rank = {}
                for a in self:
                    if a not in rank:
                        rank[a] = rank_fcn(a)
                    deg += rank[a] + 1
            elif isinstance(weights, (list, tuple)):
                rank = {}
                for a in self:
                    if a not in rank:
                        rank[a] = rank_fcn(a)
                    deg += weights[rank[a]]
            return deg

        if all(x in ZZ for x in self):
            return sum(self)

        raise TypeError("degree is not defined for your word")

    def deg_lex_less(self, other, weights=None):
        r"""
        Return ``True`` if ``self`` is degree lexicographically less than ``other``,
        and ``False`` otherwise. The weight of each letter in the ordered
        alphabet is given by ``weights``, which defaults to ``[1, 2, 3, ...]``.

        EXAMPLES::

            sage: Word([1,2,3]).deg_lex_less(Word([1,3,2]))
            True
            sage: Word([3,2,1]).deg_lex_less(Word([1,2,3]))
            False
            sage: W = Words(range(5))
            sage: W([1,2,4]).deg_lex_less(W([1,3,2]))
            False
            sage: Word("abba").deg_lex_less(Word("abbb"), dict(a=1,b=2))
            True
            sage: Word("abba").deg_lex_less(Word("baba"), dict(a=1,b=2))
            True
            sage: Word("abba").deg_lex_less(Word("aaba"), dict(a=1,b=2))
            False
            sage: Word("abba").deg_lex_less(Word("aaba"), dict(a=1,b=0))
            True
        """
        deg_self = self.degree(weights)
        deg_other = other.degree(weights)
        if deg_self != deg_other:
            return deg_self < deg_other
        return self.lex_less(other)

    def inv_lex_less(self, other):
        r"""
        Return ``True`` if ``self`` is inverse lexicographically less than ``other``.

        EXAMPLES::

            sage: Word([1,2,4]).inv_lex_less(Word([1,3,2]))
            False
            sage: Word([3,2,1]).inv_lex_less(Word([1,2,3]))
            True
        """
        if self.length() != len(other):
            return self.length() < len(other)
        return self.reversal() < other.reversal()

    def deg_inv_lex_less(self, other, weights=None):
        r"""
        Return ``True`` if the word ``self`` is degree inverse lexicographically
        less than ``other``.

        EXAMPLES::

            sage: Word([1,2,4]).deg_inv_lex_less(Word([1,3,2]))
            False
            sage: Word([3,2,1]).deg_inv_lex_less(Word([1,2,3]))
            True
        """
        d1 = self.degree(weights)
        d2 = other.degree(weights)
        if d1 != d2:
            return d1 < d2
        return self.inv_lex_less(other)

    def rev_lex_less(self, other):
        r"""
        Return ``True`` if the word ``self`` is reverse
        lexicographically less than ``other``.

        EXAMPLES::

            sage: Word([1,2,4]).rev_lex_less(Word([1,3,2]))
            True
            sage: Word([3,2,1]).rev_lex_less(Word([1,2,3]))
            False
        """
        if self.length() != len(other):
            return self.length() > len(other)
        return self.reversal() > other.reversal()

    def deg_rev_lex_less(self, other, weights=None):
        r"""
        Return ``True`` if ``self`` is degree reverse
        lexicographically less than ``other``.

        EXAMPLES::

            sage: Word([3,2,1]).deg_rev_lex_less(Word([1,2,3]))
            False
            sage: Word([1,2,4]).deg_rev_lex_less(Word([1,3,2]))
            False
            sage: Word([1,2,3]).deg_rev_lex_less(Word([1,2,4]))
            True
        """
        d1 = self.degree(weights)
        d2 = other.degree(weights)
        if d1 != d2:
            return d1 < d2
        return self.rev_lex_less(other)

    @cached_method
    def last_position_dict(self):
        r"""
        Return a dictionary that contains the last position of each letter
        in ``self``.

        EXAMPLES::

            sage: Word('1231232').last_position_dict()
            {'1': 3, '2': 6, '3': 5}
        """
        d = {}
        for i, letter in enumerate(self):
            d[letter] = i
        return d

    def _pos_in(self, other, p):
        r"""
        Return the position of the first occurrence of ``self`` starting at
        position ``p`` in ``other``.

        .. WARNING::

            This method is deprecated since 2020 and will be removed in a
            later version of SageMath.
            Use :meth:`first_occurrence` instead.

        EXAMPLES::

            sage: Word('12')._pos_in(Word('131231'), 2)
            doctest:warning
            ...
            DeprecationWarning: f._pos_in(w, start) is deprecated.
            Use w.first_occurrence(f, start) instead.
            See https://github.com/sagemath/sage/issues/30187 for details.
            2
            sage: Word('12')._pos_in(Word('131231'), 3) is None
            True
            sage: Word('32')._pos_in(Word('131231'), 0) is None
            True

        The empty word occurs in a word::

            sage: Word('')._pos_in(Word('123'), 0)
            0
            sage: Word('')._pos_in(Word(''), 0)
            0
        """
        from sage.misc.superseded import deprecation
        deprecation(30187, 'f._pos_in(w, start) is deprecated.'
                ' Use w.first_occurrence(f, start) instead.')
        return other.first_occurrence(self, p)

    def first_pos_in(self, other):
        r"""
        Return the position of the first occurrence of ``self`` in ``other``,
        or ``None`` if ``self`` is not a factor of ``other``.

        .. WARNING::

            This method is deprecated since 2020 and will be removed in a
            later version of SageMath.
            Use :meth:`first_occurrence` instead.

        EXAMPLES::

            sage: Word('12').first_pos_in(Word('131231'))
            doctest:warning
            ...
            DeprecationWarning: f.first_pos_in(w) is deprecated.
            Use w.first_occurrence(f) instead.
            See https://github.com/sagemath/sage/issues/30187 for details.
            2
            sage: Word('32').first_pos_in(Word('131231')) is None
            True
        """
        from sage.misc.superseded import deprecation
        deprecation(30187, 'f.first_pos_in(w) is deprecated.'
                ' Use w.first_occurrence(f) instead.')
        return other.first_occurrence(self)

    def find(self, sub, start=0, end=None):
        r"""
        Return the index of the first occurrence of ``sub`` in ``self``,
        such that ``sub`` is contained within ``self[start:end]``.
        Return `-1` on failure.

        INPUT:

        - ``sub`` -- string, list, tuple or word to search for

        - ``start`` -- nonnegative integer (default: `0`) specifying
          the position from which to start the search

        - ``end`` -- nonnegative integer (default: ``None``); specifying
          the position at which the search must stop. If ``None``, then
          the search is performed up to the end of the string.

        OUTPUT: nonnegative integer or `-1`

        EXAMPLES::

            sage: w = Word([0,1,0,0,1])
            sage: w.find(Word([1,0]))
            1

        The ``sub`` argument can also be a tuple or a list::

            sage: w.find([1,0])
            1
            sage: w.find((1,0))
            1

        Examples using ``start`` and ``end``::

            sage: w.find(Word([0,1]), start=1)
            3
            sage: w.find(Word([0,1]), start=1, end=5)
            3
            sage: w.find(Word([0,1]), start=1, end=4) == -1
            True
            sage: w.find(Word([1,1])) == -1
            True
            sage: w.find("aa")
            -1

        Instances of ``Word_str`` handle string inputs as well::

            sage: w = Word('abac')
            sage: w.find('a')
            0
            sage: w.find('ba')
            1

        TESTS:

        Check that :issue:`12804` is fixed::

            sage: w = Word(iter("ababab"), length='finite')
            sage: w.find("ab")
            0
            sage: w.find("ab", start=1)
            2
            sage: w.find("aa")
            -1
            sage: w.find("abc")
            -1
            sage: w = Words('ab')(tuple('babaabaaab'))
            sage: w.find('abc')
            -1
        """
        if not isinstance(sub, FiniteWord_class):
            try:
                sub = self.parent()(sub)
            except (ValueError, TypeError):
                return -1
        p = self[start:end].first_occurrence(sub)
        return -1 if p is None else p+start

    def rfind(self, sub, start=0, end=None):
        r"""
        Return the index of the last occurrence of ``sub`` in ``self``,
        such that ``sub`` is contained within ``self[start:end]``.
        Return ``-1`` on failure.

        INPUT:

        - ``sub`` -- string, list, tuple or word to search for

        - ``start`` -- nonnegative integer (default: `0`); specifying
          the position at which the search must stop

        - ``end`` -- nonnegative integer (default: ``None``); specifying
          the position from which to start the search. If ``None``, then
          the search is performed up to the end of the string.

        OUTPUT: nonnegative integer or `-1`

        EXAMPLES::

            sage: w = Word([0,1,0,0,1])
            sage: w.rfind(Word([0,1]))
            3

        The ``sub`` parameter can also be a list or a tuple::

            sage: w.rfind([0,1])
            3
            sage: w.rfind((0,1))
            3

        Examples using the argument ``start`` and ``end``::

            sage: w.rfind(Word([0,1]), end=4)
            0
            sage: w.rfind(Word([0,1]), end=5)
            3
            sage: w.rfind(Word([0,0]), start=2, end=5)
            2
            sage: w.rfind(Word([0,0]), start=3, end=5)
            -1

        Instances of ``Word_str`` handle string inputs as well::

            sage: w = Word('abac')
            sage: w.rfind('a')
            2
            sage: w.rfind(Word('a'))
            2
            sage: w.rfind([0,1])
            -1

        TESTS:

        Check that :issue:`12804` is fixed::

            sage: w = Word(iter("abab"), length='finite')
            sage: w.rfind("ab")
            2
            sage: w.rfind("ab", end=3)
            0
            sage: w.rfind("aa")
            -1
            sage: w.rfind([0,0,0])
            -1
        """
        if not isinstance(sub, FiniteWord_class):
            try:
                sub = self.parent()(sub)
            except (ValueError, TypeError):
                return -1
        L = len(sub)
        start = max(0, int(start))
        if end is None:
            i = len(self) - L
        else:
            i = min(end, len(self)) - L
        while i >= start:
            if self[i:i + L] == sub:
                return i
            i -= 1
        return -1

    def is_factor(self, other):
        r"""
        Return ``True`` if ``self`` is a factor of ``other``, and ``False`` otherwise.

        A finite word `u\in A^*` is a *factor* of a finite word `v\in A^*`
        if there exists `p,s\in A^*` such that `v=pus`.

        EXAMPLES::

            sage: u = Word('2113')
            sage: w = Word('123121332131233121132123')
            sage: u.is_factor(w)
            True
            sage: u = Word('321')
            sage: w = Word('1231241231312312312')
            sage: u.is_factor(w)
            False

        The empty word is factor of another word::

            sage: Word().is_factor(Word())
            True
            sage: Word().is_factor(Word('a'))
            True
            sage: Word().is_factor(Word([1,2,3]))
            True
            sage: Word().is_factor(Word(lambda n:n, length=5))
            True
        """
        return other.first_occurrence(self) is not None

    def factor_occurrences_in(self, other):
        r"""
        Return an iterator over all occurrences (including overlapping ones)
        of ``self`` in ``other`` in their order of appearance.

        .. WARNING::

            This method is deprecated since 2020 and will be removed in a
            later version of SageMath.
            Use :meth:`factor_occurrences_iterator` instead.

        EXAMPLES::

            sage: u = Word('121')
            sage: w = Word('121213211213')
            sage: list(u.factor_occurrences_in(w))
            doctest:warning
            ...
            DeprecationWarning: f.factor_occurrences_in(w) is deprecated.
            Use w.factor_occurrences_iterator(f) instead.
            See https://github.com/sagemath/sage/issues/30187 for details.
            [0, 2, 8]
        """
        from sage.misc.superseded import deprecation
        deprecation(30187, 'f.factor_occurrences_in(w) is deprecated.'
                ' Use w.factor_occurrences_iterator(f) instead.')
        return other.factor_occurrences_iterator(self)

    def nb_factor_occurrences_in(self, other):
        r"""
        Return the number of times ``self`` appears as a factor
        in ``other``.

        .. WARNING::

            This method is deprecated since 2020 and will be removed in a
            later version of SageMath.
            Use :meth:`number_of_factor_occurrences` instead.

        EXAMPLES::

            sage: Word('123').nb_factor_occurrences_in(Word('112332312313112332121123'))
            doctest:warning
            ...
            DeprecationWarning: f.nb_factor_occurrences_in(w) is deprecated.
            Use w.number_of_factor_occurrences(f) instead.
            See https://github.com/sagemath/sage/issues/30187 for details.
            4
            sage: Word('321').nb_factor_occurrences_in(Word('11233231231311233221123'))
            0

        An error is raised for the empty word::

            sage: Word().nb_factor_occurrences_in(Word('123'))
            Traceback (most recent call last):
            ...
            NotImplementedError: The factor must be non empty
        """
        from sage.misc.superseded import deprecation
        deprecation(30187, 'f.nb_factor_occurrences_in(w) is deprecated.'
                ' Use w.number_of_factor_occurrences(f) instead.')
        return other.number_of_factor_occurrences(self)

    def nb_subword_occurrences_in(self, other):
        r"""
        Return the number of times ``self`` appears in ``other`` as a subword.

        This corresponds to the notion of `binomial coefficient` of two
        finite words whose properties are presented in the chapter of
        Lothaire's book written by Sakarovitch and Simon [Lot1997]_.

        .. WARNING::

            This method is deprecated since 2020 and will be removed in a
            later version of SageMath.
            Use :meth:`number_of_subword_occurrences` instead.

        INPUT:

        - ``other`` -- finite word

        EXAMPLES::

            sage: tm = words.ThueMorseWord()

            sage: u = Word([0,1,0,1])
            sage: u.nb_subword_occurrences_in(tm[:1000])
            doctest:warning
            ...
            DeprecationWarning: f.nb_subword_occurrences_in(w) is deprecated.
            Use w.number_of_subword_occurrences(f) instead.
            See https://github.com/sagemath/sage/issues/30187 for details.
            2604124996

            sage: u = Word([0,1,0,1,1,0])
            sage: u.nb_subword_occurrences_in(tm[:100])
            20370432

        .. NOTE::

            This code, based on [MSSY2001]_, actually compute the number of
            occurrences of all prefixes of ``self`` as subwords in all
            prefixes of ``other``.  In particular, its complexity is
            bounded by ``len(self) * len(other)``.

        TESTS::

            sage: Word('').nb_subword_occurrences_in(Word(''))
            1
            sage: parent(_)
            Integer Ring
            sage: v,u = Word(), Word('123')
            sage: v.nb_subword_occurrences_in(u)
            1
            sage: v,u = Word('123'), Word('1133432311132311112')
            sage: v.nb_subword_occurrences_in(u)
            11
            sage: v,u = Word('4321'), Word('1132231112233212342231112')
            sage: v.nb_subword_occurrences_in(u)
            0
            sage: v,u = Word('3'), Word('122332112321213')
            sage: v.nb_subword_occurrences_in(u)
            4
            sage: v,u = Word([]), words.ThueMorseWord()[:1000]
            sage: v.nb_subword_occurrences_in(u)
            1
        """
        from sage.misc.superseded import deprecation
        deprecation(30187, 'f.nb_subword_occurrences_in(w) is deprecated.'
                ' Use w.number_of_subword_occurrences(f) instead.')
        return other.number_of_subword_occurrences(self)

    def number_of_factor_occurrences(self, other):
        r"""
        Return the number of times ``other`` appears as a factor
        in ``self``.

        INPUT:

        - ``other`` -- a non empty word

        EXAMPLES::

            sage: w = Word('112332312313112332121123')
            sage: w.number_of_factor_occurrences(Word('123'))
            4
            sage: w = Word('11233231231311233221123')
            sage: w.number_of_factor_occurrences(Word('321'))
            0

        ::

            sage: Word().number_of_factor_occurrences(Word('123'))
            0

        An error is raised for the empty word::

            sage: Word('123').number_of_factor_occurrences(Word())
            Traceback (most recent call last):
            ...
            NotImplementedError: The factor must be non empty
        """
        return sum(1 for _ in self.factor_occurrences_iterator(other))

    def number_of_subword_occurrences(self, other):
        r"""
        Return the number of times ``other`` appears in ``self`` as a subword.

        This corresponds to the notion of `binomial coefficient` of two
        finite words whose properties are presented in the chapter of
        Lothaire's book written by Sakarovitch and Simon [Lot1997]_.

        INPUT:

        - ``other`` -- finite word

        EXAMPLES::

            sage: tm = words.ThueMorseWord()
            sage: u = Word([0,1,0,1])
            sage: tm[:1000].number_of_subword_occurrences(u)
            2604124996

            sage: u = Word([0,1,0,1,1,0])
            sage: tm[:100].number_of_subword_occurrences(u)
            20370432

        .. NOTE::

            This code, based on [MSSY2001]_, actually compute the number of
            occurrences of all prefixes of ``self`` as subwords in all
            prefixes of ``other``.  In particular, its complexity is
            bounded by ``len(self) * len(other)``.

        TESTS::

            sage: Word('').number_of_subword_occurrences(Word(''))
            1
            sage: parent(_)
            Integer Ring
            sage: v,u = Word(), Word('123')
            sage: u.number_of_subword_occurrences(v)
            1
            sage: v,u = Word('123'), Word('1133432311132311112')
            sage: u.number_of_subword_occurrences(v)
            11
            sage: v,u = Word('4321'), Word('1132231112233212342231112')
            sage: u.number_of_subword_occurrences(v)
            0
            sage: v,u = Word('3'), Word('122332112321213')
            sage: u.number_of_subword_occurrences(v)
            4
            sage: v,u = Word([]), words.ThueMorseWord()[:1000]
            sage: u.number_of_subword_occurrences(v)
            1
        """
        # record the position of letters in other
        pos = defaultdict(list)
        for i, a in enumerate(other):
            pos[a].append(i)
        for a in pos:
            pos[a].reverse()

        # compute the occurrences of all prefixes of other as subwords in self
        occ = [ZZ.zero()] * (len(other)+1)
        occ[0] = ZZ.one()
        for a in self:
            for i in pos[a]:
                occ[i+1] += occ[i]

        # return only the number of occurrences of other
        return occ[-1]

    def number_of_letter_occurrences(self, letter):
        r"""
        Return the number of occurrences of ``letter`` in ``self``.

        INPUT:

        - ``letter`` -- a letter

        OUTPUT: integer

        EXAMPLES::

            sage: w = Word('abbabaab')
            sage: w.number_of_letter_occurrences('a')
            4
            sage: w.number_of_letter_occurrences('ab')
            0

        This methods is equivalent to ``list(w).count(letter)`` and
        ``tuple(w).count(letter)``, thus ``count`` is an alias for the method
        ``number_of_letter_occurrences``::

            sage: list(w).count('a')
            4
            sage: w.count('a')
            4

        But notice that if ``s`` and ``w`` are strings,
        ``Word(s).count(w)`` counts the number occurrences of ``w`` as a
        letter in ``Word(s)`` which is not the same as ``s.count(w)`` which
        counts the number of occurrences of the string ``w`` inside ``s``::

            sage: s = 'abbabaab'
            sage: s.count('ab')
            3
            sage: Word(s).count('ab')
            0

        .. SEEALSO::

            :meth:`sage.combinat.words.finite_word.FiniteWord_class.number_of_factor_occurrences`
        """
        return Integer(sum(1 for a in self if a == letter))
    count = number_of_letter_occurrences

    def _return_words_list(self, fact):
        r"""
        Return the return words as a list in the order they appear in the word.

        INPUT:

        - ``fact`` -- a non-empty finite word

        OUTPUT: a Python list of finite words

        TESTS::

            sage: Word('baccabccbacbca')._return_words_list(Word('b'))
            [word: bacca, word: bcc, word: bac]
        """
        return list(self.return_words_iterator(fact))

    def return_words(self, fact):
        r"""
        Return the set of return words of ``fact`` in ``self``.

        This is the set of all factors starting by the given factor and ending
        just before the next occurrence of this factor.
        See [Dur1998]_ and [HZ1999]_.

        INPUT:

        - ``fact`` -- a non-empty finite word

        OUTPUT: a Python set of finite words

        EXAMPLES::

            sage: Word('21331233213231').return_words(Word('2'))
            {word: 213, word: 21331, word: 233}
            sage: Word().return_words(Word('213'))
            set()
            sage: Word('121212').return_words(Word('1212'))
            {word: 12}

        ::

            sage: TM = words.ThueMorseWord()[:1000]
            sage: sorted(TM.return_words(Word([0])))
            [word: 0, word: 01, word: 011]
        """
        return set(self.return_words_iterator(fact))

    def complete_return_words(self, fact):
        r"""
        Return the set of complete return words of ``fact`` in ``self``.

        This is the set of all factors starting by the given factor and ending
        just after the next occurrence of this factor.
        See for instance [JV2000]_.

        INPUT:

        - ``fact`` -- a non-empty finite word

        OUTPUT: a Python set of finite words

        EXAMPLES::

            sage: s = Word('21331233213231').complete_return_words(Word('2'))
            sage: sorted(s)
            [word: 2132, word: 213312, word: 2332]
            sage: Word('').complete_return_words(Word('213'))
            set()
            sage: Word('121212').complete_return_words(Word('1212'))
            {word: 121212}
        """
        return set(self.complete_return_words_iterator(fact))

    def return_words_derivate(self, fact):
        r"""
        Return the word generated by mapping a letter to each occurrence of
        the return words for the given factor dropping any dangling prefix and
        suffix. See for instance [Dur1998]_.

        EXAMPLES::

            sage: Word('12131221312313122').return_words_derivate(Word('1'))
            word: 123242
        """
        tab = {}
        ret = [tab.setdefault(w, len(tab)) + 1 for w in self._return_words_list(fact)]
        from sage.combinat.words.word import Word
        return Word(ret)

    def is_quasiperiodic(self):
        r"""
        Return ``True`` if ``self`` is quasiperiodic, and ``False`` otherwise.

        A finite or infinite word `w` is *quasiperiodic* if it can be
        constructed by concatenations and superpositions of one of its proper
        factors `u`, which is called a *quasiperiod* of `w`.
        See for instance [AE1993]_, [Mar2004]_, and [GLR2008]_.

        EXAMPLES::

            sage: Word('abaababaabaababaaba').is_quasiperiodic()
            True
            sage: Word('abacaba').is_quasiperiodic()
            False
            sage: Word('a').is_quasiperiodic()
            False
            sage: Word().is_quasiperiodic()
            False
            sage: Word('abaaba').is_quasiperiodic()
            True
        """
        l = self.length()
        if l <= 1:
            return False
        for i in range(1, l - 1):
            return_lengths = [x.length() for x in self.return_words(self[:i])]
            if return_lengths:
                if max(return_lengths) <= i and self[l - i:l] == self[:i]:
                    return True
        return False

    def quasiperiods(self):
        r"""
        Return the quasiperiods of ``self`` as a list ordered from shortest to
        longest.

        Let `w` be a finite or infinite word. A *quasiperiod* of `w` is a
        proper factor `u` of `w` such that the occurrences of `u` in `w`
        entirely cover `w`, i.e., every position of `w` falls within some
        occurrence of `u` in `w`. See for instance [AE1993]_, [Mar2004]_,
        and [GLR2008]_.

        EXAMPLES::

            sage: Word('abaababaabaababaaba').quasiperiods()
            [word: aba, word: abaaba, word: abaababaaba]
            sage: Word('abaaba').quasiperiods()
            [word: aba]
            sage: Word('abacaba').quasiperiods()
            []
        """
        l = self.length()
        if l <= 1:
            return []
        Q = []
        for i in range(1, l - 1):
            return_lengths = [x.length() for x in self.return_words(self[:i])]
            if return_lengths:
                if max(return_lengths) <= i and self[l - i:l] == self[:i]:
                    Q.append(self[:i])
        return Q

    def crochemore_factorization(self):
        r"""
        Return the Crochemore factorization of ``self`` as an ordered list of
        factors.

        The *Crochemore factorization* or the *Lempel-Ziv decomposition* of a
        finite word `w` is the unique factorization: `(x_1, x_2, \ldots, x_n)`
        of `w` with each `x_i` satisfying either: C1. `x_i` is a letter that
        does not appear in `u = x_1\ldots x_{i-1}`; C2. `x_i` is the longest
        prefix of `v = x_i\ldots x_n` that also has an occurrence beginning
        within `u = x_1\ldots x_{i-1}`. See [Cro1983]_.

        EXAMPLES::

            sage: x = Word('abababb')
            sage: x.crochemore_factorization()
            (a, b, abab, b)
            sage: mul(x.crochemore_factorization()) == x
            True
            sage: y = Word('abaababacabba')
            sage: y.crochemore_factorization()
            (a, b, a, aba, ba, c, ab, ba)
            sage: mul(y.crochemore_factorization()) == y
            True
            sage: x = Word([0,1,0,1,0,1,1])
            sage: x.crochemore_factorization()
            (0, 1, 0101, 1)
            sage: mul(x.crochemore_factorization()) == x
            True
        """
        T = self.implicit_suffix_tree()
        cuts = T.LZ_decomposition()
        c = Factorization([self[cuts[i]:cuts[i+1]] for i in range(len(cuts)-1)])
        return c

    LZ_decomposition = crochemore_factorization

    def evaluation_dict(self):
        r"""
        Return a dictionary keyed by the letters occurring in ``self`` with
        values the number of occurrences of the letter.

        EXAMPLES::

            sage: Word([2,1,4,2,3,4,2]).evaluation_dict()
            {1: 1, 2: 3, 3: 1, 4: 2}
            sage: Word('badbcdb').evaluation_dict()
            {'a': 1, 'b': 3, 'c': 1, 'd': 2}
            sage: Word().evaluation_dict()
            {}

        ::

            sage: f = Word('1213121').evaluation_dict() # keys appear in random order
            {'1': 4, '2': 2, '3': 1}

        TESTS::

            sage: f = Word('1213121').evaluation_dict()
            sage: f['1'] == 4
            True
            sage: f['2'] == 2
            True
            sage: f['3'] == 1
            True
        """
        return evaluation_dict(self)

    def evaluation_sparse(self):
        r"""
        Return a list representing the evaluation of ``self``. The entries of
        the list are two-element lists ``[a, n]``, where ``a`` is a letter
        occurring in ``self`` and ``n`` is the number of occurrences of ``a`` in ``self``.

        EXAMPLES::

            sage: sorted(Word([4,4,2,5,2,1,4,1]).evaluation_sparse())
            [(1, 2), (2, 2), (4, 3), (5, 1)]
            sage: sorted(Word("abcaccab").evaluation_sparse())
            [('a', 3), ('b', 2), ('c', 3)]
        """
        return list(self.evaluation_dict().items())

    def evaluation_partition(self):
        r"""
        Return the evaluation of the word w as a partition.

        EXAMPLES::

            sage: Word("acdabda").evaluation_partition()
            [3, 2, 1, 1]
            sage: Word([2,1,4,2,3,4,2]).evaluation_partition()
            [3, 2, 1, 1]
        """
        p = sorted(self.evaluation_dict().values(), reverse=True)
        from sage.combinat.partition import Partition
        if 0 in p:
            return Partition(p[:p.index(0)])
        else:
            return Partition(p)

    def overlap_partition(self, other, delay=0, p=None, involution=None):
        r"""
        Return the partition of the alphabet induced by the overlap of
        ``self`` and ``other`` with the given ``delay``.

        The partition of the alphabet is given by the equivalence
        relation obtained from the symmetric, reflexive and transitive
        closure of the set of pairs of letters
        `R_{u,v,d} = \{ (u_k, v_{k-d}) : 0 \leq k < n, 0\leq k-d < m \}`
        where `u = u_0 u_1 \cdots u_{n-1}`, `v = v_0v_1\cdots v_{m-1}` are
        two words on the alphabet `A` and `d` is an integer.

        The equivalence relation defined by `R` is inspired from [Lab2008]_.

        INPUT:

        - ``other`` -- word on the same alphabet as ``self``
        - ``delay`` -- integer (default: `0`)
        - ``p`` -- disjoint sets data structure (default: ``None``),
          a partition of the alphabet into disjoint sets to start with.
          If ``None``, each letter start in distinct equivalence classes.
        - ``involution`` -- callable (default: ``None``); an
          involution on the alphabet. If ``involution`` is not ``None``, the relation
          `R_{u,v,d} \cup R_{involution(u),involution(v),d}` is considered.

        OUTPUT: a disjoint set data structure

        EXAMPLES::

            sage: W = Words(list('abc012345'))
            sage: u = W('abc')
            sage: v = W('01234')
            sage: u.overlap_partition(v)
            {{'0', 'a'}, {'1', 'b'}, {'2', 'c'}, {'3'}, {'4'}, {'5'}}
            sage: u.overlap_partition(v, 2)
            {{'0', 'c'}, {'1'}, {'2'}, {'3'}, {'4'}, {'5'}, {'a'}, {'b'}}
            sage: u.overlap_partition(v, -1)
            {{'0'}, {'1', 'a'}, {'2', 'b'}, {'3', 'c'}, {'4'}, {'5'}}

        You can re-use the same disjoint set and do more than one overlap::

            sage: p = u.overlap_partition(v, 2)
            sage: p
            {{'0', 'c'}, {'1'}, {'2'}, {'3'}, {'4'}, {'5'}, {'a'}, {'b'}}
            sage: u.overlap_partition(v, 1, p)
            {{'0', '1', 'b', 'c'}, {'2'}, {'3'}, {'4'}, {'5'}, {'a'}}

        The function  ``overlap_partition`` can be used to study equations
        on words. For example, if a word `w` overlaps itself with delay `d`, then
        `d` is a period of `w`::

            sage: W = Words(range(20))
            sage: w = W(range(14)); w
            word: 0,1,2,3,4,5,6,7,8,9,10,11,12,13
            sage: d = 5
            sage: p = w.overlap_partition(w, d)
            sage: m = WordMorphism(p.element_to_root_dict())
            sage: w2 = m(w); w2
            word: 56789567895678
            sage: w2.minimal_period() == d
            True

        If a word is equal to its reversal, then it is a palindrome::

            sage: W = Words(range(20))
            sage: w = W(range(17)); w
            word: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16
            sage: p = w.overlap_partition(w.reversal(), 0)
            sage: m = WordMorphism(p.element_to_root_dict())
            sage: w2 = m(w); w2
            word: 01234567876543210
            sage: w2.parent()
            Finite words over {0, 1, 2, 3, 4, 5, 6, 7, 8, 17, 18, 19}
            sage: w2.is_palindrome()
            True

        If the reversal of a word `w` is factor of its square `w^2`, then
        `w` is symmetric, i.e. the product of two palindromes::

            sage: W = Words(range(10))
            sage: w = W(range(10)); w
            word: 0123456789
            sage: p = (w*w).overlap_partition(w.reversal(), 4)
            sage: m = WordMorphism(p.element_to_root_dict())
            sage: w2 = m(w); w2
            word: 0110456654
            sage: w2.is_symmetric()
            True

        If the image of the reversal of a word `w` under an involution `f`
        is factor of its square `w^2`, then `w` is `f`-symmetric::

            sage: W = Words([-11,-9,..,11])
            sage: w = W([1,3,..,11])
            sage: w
            word: 1,3,5,7,9,11
            sage: inv = lambda x:-x
            sage: f = WordMorphism(dict( (a, inv(a)) for a in W.alphabet()))
            sage: p = (w*w).overlap_partition(f(w).reversal(), 2, involution=f)
            sage: m = WordMorphism(p.element_to_root_dict())
            sage: m(w)
            word: 1,-1,5,7,-7,-5
            sage: m(w).is_symmetric(f)
            True

        TESTS::

            sage: W = Words('abcdef')
            sage: w = W('abc')
            sage: y = W('def')
            sage: w.overlap_partition(y, -3)
            {{'a'}, {'b'}, {'c'}, {'d'}, {'e'}, {'f'}}
            sage: w.overlap_partition(y, -2)
            {{'a', 'f'}, {'b'}, {'c'}, {'d'}, {'e'}}
            sage: w.overlap_partition(y, -1)
            {{'a', 'e'}, {'b', 'f'}, {'c'}, {'d'}}
            sage: w.overlap_partition(y, 0)
            {{'a', 'd'}, {'b', 'e'}, {'c', 'f'}}
            sage: w.overlap_partition(y, 1)
            {{'a'}, {'b', 'd'}, {'c', 'e'}, {'f'}}
            sage: w.overlap_partition(y, 2)
            {{'a'}, {'b'}, {'c', 'd'}, {'e'}, {'f'}}
            sage: w.overlap_partition(y, 3)
            {{'a'}, {'b'}, {'c'}, {'d'}, {'e'}, {'f'}}
            sage: w.overlap_partition(y, 4)
            {{'a'}, {'b'}, {'c'}, {'d'}, {'e'}, {'f'}}

        ::

            sage: W = Words(range(2))
            sage: w = W([0,1,0,1,0,1]); w
            word: 010101
            sage: w.overlap_partition(w, 0)
            {{0}, {1}}
            sage: w.overlap_partition(w, 1)
            {{0, 1}}

        ::

            sage: empty = Word()
            sage: empty.overlap_partition(empty, 'yo')
            Traceback (most recent call last):
            ...
            TypeError: delay (=yo) must be an integer
            sage: empty.overlap_partition(empty,2,'yo')
            Traceback (most recent call last):
            ...
            TypeError: p(=yo) is not a DisjointSet

        The ``involution`` input can be any callable::

            sage: w = Words([-5,..,5])([-5..5])
            sage: inv = lambda x:-x
            sage: w.overlap_partition(w, 2, involution=inv)
            {{-4, -2, 0, 2, 4}, {-5, -3, -1, 1, 3, 5}}
        """
        if not isinstance(delay, (int, Integer)):
            raise TypeError("delay (=%s) must be an integer" % delay)
        elif delay < 0:
            return other.overlap_partition(self, -delay, p)

        from sage.sets.disjoint_set import DisjointSet_class
        if p is None:
            if self.parent().alphabet().cardinality() is Infinity:
                raise ValueError("The alphabet of the parent must be finite")
            from sage.sets.disjoint_set import DisjointSet
            p = DisjointSet(self.parent().alphabet())
        elif not isinstance(p, DisjointSet_class):
            raise TypeError("p(=%s) is not a DisjointSet" % p)

        # Join the classes of each pair of letters that are one above the other
        from sage.combinat.words.morphism import WordMorphism
        S = zip(islice(self, int(delay), None), other)
        if involution is None:
            for a, b in S:
                p.union(a, b)
        elif isinstance(involution, WordMorphism):
            for a, b in S:
                p.union(a, b)
                # take the first letter of the word
                p.union(involution(a)[0], involution(b)[0])
        elif callable(involution):
            for a, b in S:
                p.union(a, b)
                p.union(involution(a), involution(b))
        else:
            raise TypeError("involution (=%s) must be callable" % involution)
        return p

    # TODO: requires a parent with a sortkey_letters method
    def standard_permutation(self):
        r"""
        Return the standard permutation of the word
        ``self`` on the ordered alphabet. It is defined as
        the permutation with exactly the same inversions as
        ``self``. Equivalently, it is the permutation
        of minimal length whose inverse sorts ``self``.

        EXAMPLES::

            sage: w = Word([1,2,3,2,2,1]); w
            word: 123221
            sage: p = w.standard_permutation(); p
            [1, 3, 6, 4, 5, 2]
            sage: v = Word(p.inverse().action(w)); v
            word: 112223
            sage: [q for q in Permutations(w.length())
            ....:      if q.length() <= p.length() and
            ....:      q.inverse().action(w) == list(v)]
            [[1, 3, 6, 4, 5, 2]]

        ::

            sage: w = Words([1,2,3])([1,2,3,2,2,1,2,1]); w
            word: 12322121
            sage: p = w.standard_permutation(); p
            [1, 4, 8, 5, 6, 2, 7, 3]
            sage: Word(p.inverse().action(w))
            word: 11122223

        ::

            sage: w = Words([3,2,1])([1,2,3,2,2,1,2,1]); w
            word: 12322121
            sage: p = w.standard_permutation(); p
            [6, 2, 1, 3, 4, 7, 5, 8]
            sage: Word(p.inverse().action(w))
            word: 32222111

        ::

            sage: w = Words('ab')('abbaba'); w
            word: abbaba
            sage: p = w.standard_permutation(); p
            [1, 4, 5, 2, 6, 3]
            sage: Word(p.inverse().action(w))
            word: aaabbb

        ::

            sage: w = Words('ba')('abbaba'); w
            word: abbaba
            sage: p = w.standard_permutation(); p
            [4, 1, 2, 5, 3, 6]
            sage: Word(p.inverse().action(w))
            word: bbbaaa
        """
        from sage.combinat.permutation import to_standard
        return to_standard(self, key=self.parent().sortkey_letters)

    def _s(self, i):
        r"""
        Implement Lascoux and Schützenberger `s_i` operator, swap the
        number of ``i`` and ``i+1`` in a word.

        EXAMPLES::

            sage: w = Word([1,1,2,1,2,2,1,3,1,2])
            sage: w._s(1)
            word: 1221221312

        TESTS::

            sage: w = Word([])
            sage: w._s(1)
            word:
            sage: w = Words(3)([2,1,2])
            sage: w._s(1).parent()
            Finite words over {1, 2, 3}
        """
        unpaired_i = []  # positions of unpaired is
        unpaired_ip = []  # positions of unpaired i+1s
        for p, x in enumerate(self):
            if x == i:
                if unpaired_ip:
                    unpaired_ip.pop()
                else:
                    unpaired_i.append(p)
            elif x == i + 1:
                unpaired_ip.append(p)

        unpaired = unpaired_i + unpaired_ip

        out = list(self)

        # replace the unpaired subword i^a (i+1)^b
        # with i^b (i+1)^a

        for j, p in enumerate(unpaired):
            if j < len(unpaired_ip):
                out[p] = i
            else:
                out[p] = i + 1
        return self.parent()(out, check=False)

    def _to_partition_content(self):
        r"""
        Return the conversion of ``self`` to a word with partition content using
        the `s_i` operators of Lascoux and Schützenberger.

        EXAMPLES::

            sage: w = Word([1,3,2,1,2,3,4,6,4,2,3,2])
            sage: w._to_partition_content()
            word: 132112454132
            sage: Word([])._to_partition_content()
            word:
        """
        if self.length() == 0:
            return self

        from sage.combinat.words.word import Word
        n = max(self)
        ev = Word(Words(n)(self).evaluation())
        sig = ev.reversal().standard_permutation().reduced_word()

        # sig is now the reverse complement of a reduced word for a minimal
        # length permutation which would sort the letters of ev into a
        # partition

        out = self
        for i in reversed(sig):
            out = out._s(n-i)
        return out

    def cocharge(self):
        r"""
        Return the cocharge of ``self``.  For a word `w`, this can be defined as
        `n_{ev} - ch(w)`, where `ch(w)` is the charge of `w` and `ev` is the
        evaluation of `w`, and `n_{ev}` is `\sum_{i<j} min(ev_i, ev_j)`.

        EXAMPLES::

            sage: Word([1,2,3]).cocharge()
            0
            sage: Word([3,2,1]).cocharge()
            3
            sage: Word([1,1,2]).cocharge()
            0
            sage: Word([2,1,2]).cocharge()
            1

        TESTS::

            sage: Word([]).cocharge()
            0
        """
        return self.evaluation_partition().weighted_size() - self.charge()

    def charge(self, check=True):
        r"""
        Return the charge of ``self``.  This is defined as follows.

        If `w` is a permutation of length `n`, (in other words, the evaluation
        of `w` is `(1, 1, \dots, 1)`), the statistic charge(`w`) is given by
        `\sum_{i=1}^n c_i(w)` where `c_1(w) = 0` and `c_i(w)` is defined
        recursively by setting `p_i` equal to `1` if `i` appears to the right
        of `i-1` in `w` and `0` otherwise.  Then we set `c_i(w) = c_{i-1}(w) +
        p_i`.

        EXAMPLES::

            sage: Word([1, 2, 3]).charge()
            3
            sage: Word([3, 5, 1, 4, 2]).charge() == 0 + 1 + 1 + 2 + 2
            True

        If `w` is not a permutation, but the evaluation of `w` is a partition,
        the charge of `w` is defined to be the sum of its charge subwords
        (each of which will be a permutation).  The first charge subword is
        found by starting at the end of `w` and moving left until the first
        `1` is found.  This is marked, and we continue to move to the left
        until the first `2` is found, wrapping around from the beginning of
        the word back to the end, if necessary.  We mark this `2`, and
        continue on until we have marked the largest letter in `w`.  The
        marked letters, with relative order preserved, form the first charge
        subword of `w`.  This subword is removed, and the next charge subword
        is found in the same manner from the remaining letters.  In the
        following example, `w1, w2, w3` are the charge subwords of `w`.

        EXAMPLES::

            sage: w = Word([5,2,3,4,4,1,1,1,2,2,3])
            sage: w1 = Word([5, 2, 4, 1, 3])
            sage: w2 = Word([3, 4, 1, 2])
            sage: w3 = Word([1, 2])
            sage: w.charge() == w1.charge() + w2.charge() + w3.charge()
            True

        Finally, if `w` does not have partition content, we apply the
        Lascoux-Schützenberger standardization operators `s_i` in such a
        manner as to obtain a word with partition content. (The word we obtain
        is independent of the choice of operators.)  The charge is then
        defined to be the charge of this word::

            sage: Word([3,3,2,1,1]).charge()
            0
            sage: Word([1,2,3,1,2]).charge()
            2

        Note that this differs from the definition of charge given in
        Macdonald's book.  The difference amounts to a choice of
        reading a word from left-to-right or right-to-left.  The choice in
        Sage was made to agree with the definition of a reading word of a
        tableau in Sage, and seems to be the more common convention in the
        literature.

        See [Mac1995]_, [LLM2003]_, and [LLT]_.

        TESTS::

            sage: Word([1,1,2,2,3]).charge()
            4
            sage: Word([3,1,1,2,2]).charge()
            3
            sage: Word([2,1,1,2,3]).charge()
            2
            sage: Word([2,1,1,3,2]).charge()
            2
            sage: Word([3,2,1,1,2]).charge()
            1
            sage: Word([2,2,1,1,3]).charge()
            1
            sage: Word([3,2,2,1,1]).charge()
            0
            sage: Word([]).charge()
            0
        """
        if check:
            ev_dict = self.evaluation_dict()
            ordered_alphabet = sorted(ev_dict,
                                      key=self.parent().sortkey_letters)
            evaluation = [ev_dict[a] for a in ordered_alphabet]
            from sage.combinat.partition import Partitions
            if evaluation not in Partitions():
                return self._to_partition_content().charge()
        res = 0
        w = self.to_integer_list()
        while w:
            i = len(w) - 1
            l = min(w)
            index = 0
            while w and l <= max(w):
                while w[i] != l:
                    i -= 1
                    if i < 0:
                        i = len(w) - 1
                        index += 1
                res += index
                l += 1
                w.pop(i)
                i -= 1
                if i < 0:
                    i = len(w) - 1
                    index += 1
        return res

    def BWT(self):
        r"""
        Return the Burrows-Wheeler Transform (BWT) of ``self``.

        The *Burrows-Wheeler transform* of a finite word `w` is obtained
        from `w` by first listing the conjugates of `w` in lexicographic order
        and then concatenating the final letters of the conjugates in this
        order. See [BW1994]_.

        EXAMPLES::

            sage: Word('abaccaaba').BWT()
            word: cbaabaaca
            sage: Word('abaab').BWT()
            word: bbaaa
            sage: Word('bbabbaca').BWT()
            word: cbbbbaaa
            sage: Word('aabaab').BWT()
            word: bbaaaa
            sage: Word().BWT()
            word:
            sage: Word('a').BWT()
            word: a
        """
        if self.is_empty():
            return self
        conjugates = sorted(self._conjugates_list())
        return self.parent()([x[x.length() - 1] for x in conjugates],
                             check=False)

    def iterated_left_palindromic_closure(self, f=None):
        r"""
        Return the iterated left (``f``-)palindromic closure of ``self``.

        INPUT:

        - ``f`` -- involution (default: ``None``) on the alphabet of ``self``;
          it must be callable on letters as well as words (e.g. ``WordMorphism``)

        OUTPUT: word; the left iterated ``f``-palindromic closure of ``self``

        EXAMPLES::

            sage: Word('123').iterated_left_palindromic_closure()
            word: 3231323
            sage: f = WordMorphism('a->b,b->a')
            sage: Word('ab').iterated_left_palindromic_closure(f=f)
            word: abbaab
            sage: Word('aab').iterated_left_palindromic_closure(f=f)
            word: abbaabbaab

        TESTS:

        If ``f`` is not an involution::

            sage: f = WordMorphism('a->b,b->b')
            sage: Word('aab').iterated_left_palindromic_closure(f=f)
            Traceback (most recent call last):
            ...
            TypeError: self (=a->b, b->b) is not an endomorphism

        See [DeLuca2006]_.
        """
        if f is None:
            return self.reversal().iterated_right_palindromic_closure(f=f)
        else:
            from sage.combinat.words.morphism import WordMorphism
            f = WordMorphism(f)
            return f(self).reversal().iterated_right_palindromic_closure(f=f)

    def balance(self):
        r"""
        Return the balance of ``self``.

        The balance of a word is the smallest number `q` such that ``self`` is
        `q`-balanced [FV2002]_.

        A finite or infinite word `w` is said to be `q`-*balanced* if for
        any two factors `u`, `v` of `w` of the same length, the difference
        between the number of `x`'s in each of `u` and `v` is at most `q`
        for all letters `x` in the alphabet of `w`. A `1`-balanced word is
        simply said to be balanced. See Chapter 2 of [Lot2002]_.

        OUTPUT: integer

        EXAMPLES::

            sage: Word('1111111').balance()
            0
            sage: Word('001010101011').balance()
            2
            sage: Word('0101010101').balance()
            1

        ::

            sage: w = Word('11112222')
            sage: w.is_balanced(2)
            False
            sage: w.is_balanced(3)
            False
            sage: w.is_balanced(4)
            True
            sage: w.is_balanced(5)
            True
            sage: w.balance()
            4

        TESTS::

            sage: Word('1111122222').balance()
            5
            sage: Word('').balance()
            0
            sage: Word('1').balance()
            0
            sage: Word('12').balance()
            1
            sage: Word('1112').balance()
            1
        """
        alphabet = self.letters()
        best = 0
        for i in range(1, self.length()):
            start = iter(self)
            end = iter(self)
            abelian = dict(zip(alphabet, [0] * len(alphabet)))
            for _ in range(i):
                abelian[next(end)] += 1
            abel_max = abelian.copy()
            abel_min = abelian.copy()
            for _ in range(self.length() - i):
                lost = next(start)
                gain = next(end)
                abelian[gain] += 1
                abelian[lost] -= 1
                abel_max[gain] = max(abel_max[gain], abelian[gain])
                abel_min[lost] = min(abel_min[lost], abelian[lost])
            best = max(best, max(abel_max[a] - abel_min[a] for a in alphabet))
        return best

    def is_balanced(self, q=1):
        r"""
        Return ``True`` if ``self`` is ``q``-balanced, and ``False`` otherwise.

        A finite or infinite word `w` is said to be `q`-*balanced* if for
        any two factors `u`, `v` of `w` of the same length, the difference
        between the number of `x`'s in each of `u` and `v` is at most `q`
        for all letters `x` in the alphabet of `w`. A `1`-balanced word is
        simply said to be balanced. See for instance [CFZ2000]_ and Chapter
        2 of [Lot2002]_.

        INPUT:

        - ``q`` -- integer (default: `1`); the balance level

        EXAMPLES::

            sage: Word('1213121').is_balanced()
            True
            sage: Word('1122').is_balanced()
            False
            sage: Word('121333121').is_balanced()
            False
            sage: Word('121333121').is_balanced(2)
            False
            sage: Word('121333121').is_balanced(3)
            True
            sage: Word('121122121').is_balanced()
            False
            sage: Word('121122121').is_balanced(2)
            True

        TESTS::

            sage: Word('121122121').is_balanced(-1)
            Traceback (most recent call last):
            ...
            TypeError: the balance level must be a positive integer
            sage: Word('121122121').is_balanced(0)
            Traceback (most recent call last):
            ...
            TypeError: the balance level must be a positive integer
            sage: Word('121122121').is_balanced('a')
            Traceback (most recent call last):
            ...
            TypeError: the balance level must be a positive integer
        """
        if not isinstance(q, (int, Integer)) or q <= 0:
            raise TypeError("the balance level must be a positive integer")
        alphabet = self.letters()
        for i in range(2, self.length()):
            empty_sets = [set() for _ in range(len(alphabet))]
            tab = dict(zip(alphabet, empty_sets))
            for fact in self.factor_iterator(i):
                evaluation_dict = fact.evaluation_dict()
                for a in alphabet:
                    tab[a].add(evaluation_dict.get(a, 0))
            for t in tab.values():
                if len(t) > q+1:
                    return False
        return True

    def abelian_vectors(self, n):
        r"""
        Return the abelian vectors of factors of length ``n`` of ``self``.

        The vectors are defined w.r.t. the order of the alphabet of the
        parent.

        OUTPUT: a set of tuples

        EXAMPLES::

            sage: W = Words([0,1,2])
            sage: w = W([0,1,1,0,1,2,0,2,0,2])
            sage: w.abelian_vectors(3)
            {(1, 0, 2), (1, 1, 1), (1, 2, 0), (2, 0, 1)}
            sage: w[:5].abelian_vectors(3)
            {(1, 2, 0)}
            sage: w[5:].abelian_vectors(3)
            {(1, 0, 2), (2, 0, 1)}

        ::

            sage: w = words.FibonacciWord()[:100]
            sage: sorted(w.abelian_vectors(0))
            [(0, 0)]
            sage: sorted(w.abelian_vectors(1))
            [(0, 1), (1, 0)]
            sage: sorted(w.abelian_vectors(7))
            [(4, 3), (5, 2)]

        The word must be defined with a parent on a finite alphabet::

            sage: from itertools import count
            sage: w = Word(count(), alphabet=NN)
            sage: w[:2].abelian_vectors(2)
            Traceback (most recent call last):
            ...
            TypeError: The alphabet of the parent is infinite; define the
            word with a parent on a finite alphabet

        TESTS::

            sage: W = Words([0, 1])
            sage: w = W([0,0,0])
            sage: sorted(w.abelian_vectors(3))
            [(3, 0)]
            sage: w = W([0,0,0,1])
            sage: sorted(w.abelian_vectors(3))
            [(2, 1), (3, 0)]
            sage: w = W([0,0,0,1,1])
            sage: sorted(w.abelian_vectors(3))
            [(1, 2), (2, 1), (3, 0)]
            sage: w = W([0,0,0,1,1,1])
            sage: sorted(w.abelian_vectors(3))
            [(0, 3), (1, 2), (2, 1), (3, 0)]

        ::

            sage: w = Word([0,1,0], alphabet=[0,1])
            sage: w.abelian_complexity(3)
            1
            sage: w.abelian_complexity(4)
            0
        """
        alphabet = self.parent().alphabet()
        size = alphabet.cardinality()
        if size == float('inf'):
            raise TypeError("The alphabet of the parent is infinite; define"
                   " the word with a parent on a finite alphabet")
        S = set()
        if n > self.length():
            return S
        rank = {letter: i for i, letter in enumerate(alphabet)}
        start = iter(self)
        end = iter(self)
        abelian = [0] * size
        for _ in range(n):
            abelian[rank[next(end)]] += 1
        S.add(tuple(abelian))
        for letter in end:
            abelian[rank[letter]] += 1
            abelian[rank[next(start)]] -= 1
            S.add(tuple(abelian))
        return S

    def abelian_complexity(self, n):
        r"""
        Return the number of abelian vectors of factors of length ``n`` of ``self``.

        EXAMPLES::

            sage: w = words.FibonacciWord()[:100]
            sage: [w.abelian_complexity(i) for i in range(20)]
            [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]

        ::

            sage: w = words.ThueMorseWord()[:100]
            sage: [w.abelian_complexity(i) for i in range(20)]
            [1, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2]
        """
        return len(self.abelian_vectors(n))

    def sturmian_desubstitute_as_possible(self):
        r"""
        Sturmian-desubstitute the word ``self`` as much as possible.

        The finite word ``self`` must be defined on a two-letter
        alphabet or use at most two letters.

        It can be Sturmian desubstituted if one letter appears
        isolated: the Sturmian desubstitution consists in removing one
        letter per run of the non-isolated letter. The accelerated
        Sturmian desubstitution consists in removing a run equal to
        the length of the shortest inner run from any run of the
        non-isolated letter (including possible leading and trailing
        runs even if they have shorter length). The (accelerated)
        Sturmian desubstitution is done as much as possible. A word is
        a factor of a Sturmian word if, and only if, the result is the
        empty word.

        OUTPUT: a finite word defined on a two-letter alphabet

        EXAMPLES::

            sage: u = Word('10111101101110111',alphabet='01') ; u
            word: 10111101101110111
            sage: v = u.sturmian_desubstitute_as_possible() ; v
            word: 01100101
            sage: v == v.sturmian_desubstitute_as_possible()
            True

            sage: Word('azaazaaazaaazaazaaaz', alphabet='az').sturmian_desubstitute_as_possible()
            word:


        TESTS::

            sage: w = Word('azazaza', alphabet='aze')
            sage: w.sturmian_desubstitute_as_possible()
            word:
            sage: Word('aze').sturmian_desubstitute_as_possible()
            Traceback (most recent call last):
            ...
            TypeError: your word must be defined on a binary alphabet or use at most two different letters
            sage: Word('azaaazaazaazaaazaaza', alphabet='az').sturmian_desubstitute_as_possible()
            word:
            sage: Word('azaaazaazaazaaazaaaza', alphabet='az').sturmian_desubstitute_as_possible()
            word: azzaa

        Boundary effects::

            sage: Word('', alphabet='az').sturmian_desubstitute_as_possible()
            word:
            sage: Word('azzzzz', alphabet='az').sturmian_desubstitute_as_possible()
            word:
            sage: Word('zzzzza', alphabet='az').sturmian_desubstitute_as_possible()
            word:
            sage: Word('aaaaazaaaaaaaaa', alphabet='az').sturmian_desubstitute_as_possible()
            word:
            sage: Word('aaaaaaaaaaaaaa', alphabet='az').sturmian_desubstitute_as_possible()
            word:

        Boundary effects without alphabet::

            sage: Word('').sturmian_desubstitute_as_possible()
            word:
            sage: Word('azzzzz').sturmian_desubstitute_as_possible()
            word:
            sage: Word('zzzzza').sturmian_desubstitute_as_possible()
            word:
            sage: Word('aaaaazaaaaaaaaa').sturmian_desubstitute_as_possible()
            word:
            sage: Word('aaaaaaaaaaaaaa').sturmian_desubstitute_as_possible()
            word:

        Idempotence::

            sage: r = words.RandomWord(randint(1,15)).sturmian_desubstitute_as_possible() ; r == r.sturmian_desubstitute_as_possible()
            True

        AUTHOR:

        -   Thierry Monteil
        """
        if self.is_empty():
            return self
        W = self.parent()
        if W.alphabet().cardinality() == 2:
            alphabet = W.alphabet()
        else:
            alphabet = self.letters()
            if len(alphabet) > 2:
                raise TypeError('your word must be defined on a binary alphabet or use at most two different letters')
            elif len(alphabet) < 2:
                return W()
        word_from_letter = {l: W([l], datatype='list', check=False) for l in alphabet}
        is_prefix = True
        current_run_length = 0
        prefix_length = 0
        prefix_letter = self[0]
        is_isolated = {alphabet[0]: True, alphabet[1]: True}
        minimal_run = {alphabet[0]: Infinity, alphabet[1]: Infinity}
        maximal_run = {alphabet[0]: 0, alphabet[1]: 0}
        runs = {alphabet[0]: [], alphabet[1]: []}
        for i in self:
            if is_prefix:
                if i == prefix_letter:
                    prefix_length += 1
                    if prefix_length > 1:
                        is_isolated[i] = False
                else:
                    is_prefix = False
                    current_run_length = 1
                    previous_letter = i
            else:
                if i == previous_letter:
                    current_run_length += 1
                    if current_run_length > 1:
                        is_isolated[i] = False
                else:
                    runs[previous_letter].append(current_run_length)
                    minimal_run[previous_letter] = min(minimal_run[previous_letter], current_run_length)
                    maximal_run[previous_letter] = max(maximal_run[previous_letter], current_run_length)
                    current_run_length = 1
                    previous_letter = i
        # at this point, previous_letter is the suffix letter and current_run_length is the suffix length
        if not (is_isolated[alphabet[0]] or is_isolated[alphabet[1]]):
            return self
        elif is_isolated[alphabet[0]] and is_isolated[alphabet[1]]:
            return W()
        else:
            if is_isolated[alphabet[0]]:
                l_isolated = alphabet[0]  # the isolated letter
                l_running = alphabet[1]  # the running letter (non-isolated)
            else:
                l_isolated = alphabet[1]
                l_running = alphabet[0]
            w_isolated = word_from_letter[l_isolated]  # the word associated to the isolated letter
            w_running = word_from_letter[l_running]  # the word associated to the running letter
            min_run = minimal_run[l_running]
            if prefix_letter == l_isolated or prefix_length <= min_run:
                desubstitued_word = W()
            else:
                desubstitued_word = w_running ** (prefix_length - min_run)
            for i in runs[l_running]:
                desubstitued_word = desubstitued_word + w_isolated + w_running ** (i - min_run)
            if current_run_length > 0:
                desubstitued_word = desubstitued_word + w_isolated
                if previous_letter == l_running and current_run_length > min_run:
                    desubstitued_word = desubstitued_word + w_running ** (current_run_length - min_run)
            return desubstitued_word.sturmian_desubstitute_as_possible()

    def is_sturmian_factor(self):
        r"""
        Tell whether ``self`` is a factor of a Sturmian word.

        The finite word ``self`` must be defined on a two-letter alphabet.

        Equivalently, tells whether ``self`` is balanced. The
        advantage over the ``is_balanced`` method is that this one runs in
        linear time whereas ``is_balanced`` runs in quadratic time.

        OUTPUT: boolean

        EXAMPLES::

            sage: w = Word('0111011011011101101',alphabet='01')
            sage: w.is_sturmian_factor()
            True

        ::

            sage: words.LowerMechanicalWord(random(),alphabet='01')[:100].is_sturmian_factor()
            True
            sage: words.CharacteristicSturmianWord(random())[:100].is_sturmian_factor()             # needs sage.rings.real_mpfr
            True

        ::

            sage: w = Word('aabb',alphabet='ab')
            sage: w.is_sturmian_factor()
            False

            sage: s1 = WordMorphism('a->ab,b->b')
            sage: s2 = WordMorphism('a->ba,b->b')
            sage: s3 = WordMorphism('a->a,b->ba')
            sage: s4 = WordMorphism('a->a,b->ab')
            sage: W = Words('ab')
            sage: w = W('ab')
            sage: for i in range(8): w = choice([s1,s2,s3,s4])(w)
            sage: w.is_sturmian_factor()
            True

        Famous words::

            sage: words.FibonacciWord()[:100].is_sturmian_factor()
            True
            sage: words.ThueMorseWord()[:1000].is_sturmian_factor()
            False
            sage: words.KolakoskiWord()[:1000].is_sturmian_factor()
            False

        See [Arn2002]_, [Ser1985]_, and [SU2009]_.

        AUTHOR:

        -   Thierry Monteil
        """
        return self.sturmian_desubstitute_as_possible().is_empty()

    def is_tangent(self):
        r"""
        Tell whether ``self`` is a tangent word.

        The finite word ``self`` must be defined on a two-letter alphabet.

        A binary word is said to be *tangent* if it can appear in
        infinitely many cutting sequences of a smooth curve, where each
        cutting sequence is observed on a progressively smaller grid.

        This class of words strictly contains the class of `1`-balanced
        words, and is strictly contained in the class of `2`-balanced words.

        This method runs in linear time.

        OUTPUT: boolean

        EXAMPLES::

            sage: w = Word('01110110110111011101',alphabet='01')
            sage: w.is_tangent()
            True

        Some tangent words may not be balanced::

            sage: Word('aabb',alphabet='ab').is_balanced()
            False
            sage: Word('aabb',alphabet='ab').is_tangent()
            True

        Some `2`-balanced words may not be tangent::

            sage: Word('aaabb',alphabet='ab').is_tangent()
            False
            sage: Word('aaabb',alphabet='ab').is_balanced(2)
            True

        Famous words::

            sage: words.FibonacciWord()[:100].is_tangent()
            True
            sage: words.ThueMorseWord()[:1000].is_tangent()
            True
            sage: words.KolakoskiWord()[:1000].is_tangent()
            False

        See [Mon2010]_.

        AUTHOR:

        -   Thierry Monteil
        """
        if self.parent().alphabet().cardinality() != 2:
            raise TypeError('your word must be defined on a binary alphabet')
        a, b = self.parent().alphabet()
        mini = 0
        maxi = 0
        height = 0
        for i in self.sturmian_desubstitute_as_possible():
            if i == a:
                height = height + 1
                maxi = max(maxi, height)
            if i == b:
                height = height - 1
                mini = min(mini, height)
        return maxi - mini <= 2

    # TODO.
    # 1. Those three swap functions should use the cmp of python
    # 2. The actual code should then be copied as is in the Word_over_Alphabet
    # and continue to use the parent cmp
    # 3. Once Word can define Words over alphabet, the examples
    # should be updated appropriately.

    def swap(self, i, j=None):
        r"""
        Return the word `w` with entries at positions ``i`` and
        ``j`` swapped. By default, ``j = i+1``.

        EXAMPLES::

            sage: Word([1,2,3]).swap(0,2)
            word: 321
            sage: Word([1,2,3]).swap(1)
            word: 132
            sage: Word("abba").swap(1,-1)
            word: aabb
        """
        if j is None:
            j = i+1
        new = list(self)
        (new[i], new[j]) = (new[j], new[i])
        from sage.combinat.words.word import Word
        return Word(new)

    def swap_increase(self, i):
        r"""
        Return the word with positions ``i`` and ``i+1`` exchanged
        if ``self[i] > self[i+1]``. Otherwise, it returns ``self``.

        EXAMPLES::

            sage: w = Word([1,3,2])
            sage: w.swap_increase(1)
            word: 123
            sage: w.swap_increase(0)
            word: 132
            sage: w.swap_increase(0) is w
            True
            sage: Words("ab")("abba").swap_increase(0)
            word: abba
            sage: Words("ba")("abba").swap_increase(0)
            word: baba
        """
        key = self._parent.sortkey_letters
        if key(self[i]) > key(self[i + 1]):
            return self.swap(i)
        else:
            return self

    def swap_decrease(self, i):
        r"""
        Return the word with positions ``i`` and ``i+1`` exchanged
        if ``self[i] < self[i+1]``. Otherwise, it returns ``self``.

        EXAMPLES::

            sage: w = Word([1,3,2])
            sage: w.swap_decrease(0)
            word: 312
            sage: w.swap_decrease(1)
            word: 132
            sage: w.swap_decrease(1) is w
            True
            sage: Words("ab")("abba").swap_decrease(0)
            word: baba
            sage: Words("ba")("abba").swap_decrease(0)
            word: abba
        """
        key = self._parent.sortkey_letters
        if key(self[i]) < key(self[i + 1]):
            return self.swap(i)
        else:
            return self

    def abelian_vector(self):
        r"""
        Return the abelian vector of ``self`` counting the occurrences of each letter.

        The vector is defined w.r.t. the order of the alphabet of the
        parent. See also :meth:`evaluation_dict`.

        INPUT:

        - ``self`` -- word having a parent on a finite alphabet

        OUTPUT: list

        EXAMPLES::

            sage: W = Words('ab')
            sage: W('aaabbbbb').abelian_vector()
            [3, 5]
            sage: W('a').abelian_vector()
            [1, 0]
            sage: W().abelian_vector()
            [0, 0]

        The result depends on the alphabet of the parent::

            sage: W = Words('abc')
            sage: W('aabaa').abelian_vector()
            [4, 1, 0]

        TESTS::

            sage: W = Words()
            sage: W('aabaa').abelian_vector()
            Traceback (most recent call last):
            ...
            TypeError: The alphabet of the parent is infinite; define the
            word with a parent on a finite alphabet or use
            evaluation_dict() instead
        """
        alphabet = self.parent().alphabet()
        if alphabet.cardinality() is Infinity:
            raise TypeError("The alphabet of the parent is infinite; define "
                    "the word with a parent on a finite alphabet or use "
                    "evaluation_dict() instead")
        ev_dict = self.evaluation_dict()
        return [ev_dict.get(a, 0) for a in alphabet]

    evaluation = abelian_vector

    def robinson_schensted(self):
        """
        Return the semistandard tableau and standard tableau pair
        obtained by running the Robinson-Schensted algorithm on ``self``.

        This can also be done by running
        :func:`~sage.combinat.rsk.RSK` on ``self``.

        EXAMPLES::

            sage: Word([1,1,3,1,2,3,1]).robinson_schensted()
            [[[1, 1, 1, 1, 3], [2], [3]], [[1, 2, 3, 5, 6], [4], [7]]]
        """
        from sage.combinat.rsk import RSK
        return RSK(self)

    def _rsk_iter(self):
        r"""
        Return an iterator for :func:`~sage.combinat.rsk.RSK`.

        Yields pairs `(i, w_i)` for a word `w = w_1 w_2 \cdots w_k`.

        EXAMPLES::

            sage: for x in Word([1,1,3,1,2,3,1])._rsk_iter(): x
            (1, 1)
            (2, 1)
            (3, 3)
            (4, 1)
            (5, 2)
            (6, 3)
            (7, 1)
        """
        return zip(range(1, len(self) + 1), self)

    def shuffle(self, other, overlap=0):
        r"""
        Return the combinatorial class representing the shuffle product
        between words ``self`` and ``other``. This consists of all words of length
        ``self.length()+other.length()`` that have both ``self`` and ``other`` as
        subwords.

        If ``overlap`` is nonzero, then the combinatorial class representing
        the shuffle product with overlaps is returned. The calculation of
        the shift in each overlap is done relative to the order of the
        alphabet. For example, `a` shifted by `a` is `b` in the alphabet
        `[a, b, c]` and `0` shifted by `1` in `[0, 1, 2, 3]` is `2`.

        INPUT:

        - ``other`` -- finite word
        - ``overlap`` -- (default: ``0``) integer or ``True``

        OUTPUT: combinatorial class of shuffle product of ``self`` and ``other``

        EXAMPLES::

            sage: ab = Word("ab")
            sage: cd = Word("cd")
            sage: sp = ab.shuffle(cd); sp
            Shuffle product of word: ab and word: cd
            sage: sp.cardinality()
            6
            sage: sp.list()
            [word: abcd, word: acbd, word: acdb, word: cabd, word: cadb, word: cdab]
            sage: w = Word([0,1])
            sage: u = Word([2,3])
            sage: w.shuffle(w)
            Shuffle product of word: 01 and word: 01
            sage: u.shuffle(u)
            Shuffle product of word: 23 and word: 23
            sage: w.shuffle(u)
            Shuffle product of word: 01 and word: 23
            sage: sp2 = w.shuffle(u,2); sp2
            Overlapping shuffle product of word: 01 and word: 23 with 2 overlaps
            sage: list(sp2)
            [word: 24]
        """
        if overlap == 0:
            from sage.combinat.words.shuffle_product import ShuffleProduct_w1w2
            return ShuffleProduct_w1w2(self, other)
        if any(a not in ZZ for a in self) or any(a not in ZZ for a in other):
            raise ValueError("for a nonzero overlap, words must contain integers as letters")
        if overlap is True:
            from sage.combinat.shuffle import ShuffleProduct_overlapping
            return ShuffleProduct_overlapping(self, other, self.parent())
        elif isinstance(overlap, (int, Integer)):
            from sage.combinat.shuffle import ShuffleProduct_overlapping_r
            return ShuffleProduct_overlapping_r(self, other, overlap, self.parent())
        raise ValueError('overlapping must be True or an integer')

    def shifted_shuffle(self, other, shift=None):
        r"""
        Return the combinatorial class representing the shifted shuffle
        product between words ``self`` and ``other``. This is the same as the
        shuffle product of ``self`` with the word obtained from ``other`` by
        incrementing its values (i.e. its letters) by the given ``shift``.

        INPUT:

        - ``other`` -- finite word over the integers
        - ``shift`` -- integer or ``None`` (default: ``None``); added to each letter of
          ``other``. When ``shift`` is ``None``, it is replaced by ``self.length()``

        OUTPUT: combinatorial class of shifted shuffle products of ``self`` and ``other``

        EXAMPLES::

            sage: w = Word([0,1,1])
            sage: sp = w.shifted_shuffle(w); sp
            Shuffle product of word: 011 and word: 344
            sage: sp = w.shifted_shuffle(w, 2); sp
            Shuffle product of word: 011 and word: 233
            sage: sp.cardinality()
            20
            sage: WordOptions(identifier='')
            sage: sp.list()
            [011233, 012133, 012313, 012331, 021133, 021313, 021331, 023113, 023131, 023311, 201133, 201313, 201331, 203113, 203131, 203311, 230113, 230131, 230311, 233011]
            sage: WordOptions(identifier='word: ')
            sage: y = Word('aba')
            sage: y.shifted_shuffle(w,2)
            Traceback (most recent call last):
            ...
            ValueError: for shifted shuffle, words must only contain integers as letters
        """
        if any(a not in ZZ for a in self) or any(a not in ZZ for a in other):
            raise ValueError("for shifted shuffle, words must only contain integers as letters")
        if shift is None:
            from sage.combinat.words.shuffle_product import ShuffleProduct_shifted
            return ShuffleProduct_shifted(self, other)
        else:
            return self.shuffle(self._parent([x + shift for x in other], check=False))

    ######################################################################
    # XXX TODO: The description is inconsistent w.r.t. "W" and "alphabet".
    ######################################################################
    def delta_inv(self, W=None, s=None):
        r"""
        Lift ``self`` via the delta operator to obtain a word containing the
        letters in alphabet (default: ``[0, 1]``). The letters used in the
        construction start with ``s`` (default: ``alphabet[0]``) and cycle
        through alphabet.

        INPUT:

        - ``alphabet`` -- an iterable
        - ``s`` -- an object in the iterable

        EXAMPLES::

            sage: W = Words([1, 2])
            sage: W([2, 2, 1, 1]).delta_inv()
            word: 112212
            sage: W([1, 1, 1, 1]).delta_inv(Words('123'))
            word: 1231
            sage: W([2, 2, 1, 1, 2]).delta_inv(s=2)
            word: 22112122
        """
        alphabet = [1, 2] if W is None else W.alphabet()
        cycle_alphabet = cycle(alphabet)
        if self.is_empty():
            return Words(alphabet)()
        if s is None:
            s = next(cycle_alphabet)
        else:
            if s not in alphabet:
                raise ValueError("starting letter not in alphabet")
            t = next(cycle_alphabet)
            while t != s:
                t = next(cycle_alphabet)
        w = []
        for i in self:
            w.extend([s] * i)
            s = next(cycle_alphabet)
        return Words(alphabet)(w)

    def delta(self):
        r"""
        Return the image of ``self`` under the delta morphism.

        The delta morphism, also known as the run-length encoding,
        is the word composed of the length of consecutive runs of
        the same letter in a given word.

        EXAMPLES::

            sage: W = Words('0123456789')
            sage: W('22112122').delta()
            word: 22112
            sage: W('555008').delta()
            word: 321
            sage: W().delta()
            word:
            sage: Word('aabbabaa').delta()
            word: 22112
        """
        if self.is_empty():
            return Words()([])
        ss = self[0]
        c = 0
        v = []
        max_c = 0
        for s in self:
            if s == ss:
                c += 1
                if c > max_c:
                    max_c = c
            else:
                v.append(c)
                ss = s
                c = 1
        v.append(c)
        return Words(list(range(1, 1 + max_c)))(v)

    # TODO. Decide whether delta_derivate* really need W.alphabet().last()....
    # RENAME: Should "derivate" be derivative?!

    def delta_derivate(self, W=None):
        r"""
        Return the derivative under delta for ``self``.

        EXAMPLES::

            sage: W = Words('12')
            sage: W('12211').delta_derivate()
            word: 22
            sage: W('1').delta_derivate(Words([1]))
            word: 1
            sage: W('2112').delta_derivate()
            word: 2
            sage: W('2211').delta_derivate()
            word: 22
            sage: W('112').delta_derivate()
            word: 2
            sage: W('11222').delta_derivate(Words([1, 2, 3]))
            word: 3
        """
        d = self.delta()
        if len(d) == 0:
            return d
        if W is None:
            W = d.parent()
        if d[0] != W.alphabet().last():
            d = d[1:]
        if d[-1] != W.alphabet().last():
            d = d[:-1]
        return d

    def delta_derivate_left(self, W=None):
        r"""
        Return the derivative under delta for ``self``.

        EXAMPLES::

            sage: W = Words('12')
            sage: W('12211').delta_derivate_left()
            word: 22
            sage: W('1').delta_derivate_left(Words([1]))
            word: 1
            sage: W('2112').delta_derivate_left()
            word: 21
            sage: W('2211').delta_derivate_left()
            word: 22
            sage: W('112').delta_derivate_left()
            word: 21
            sage: W('11222').delta_derivate_left(Words([1, 2, 3]))
            word: 3
        """
        d = self.delta()
        if len(d) == 0:
            return d
        if W is None:
            W = d.parent()
        if d[0] != W.alphabet().last():
            d = d[1:]
        return d

    def delta_derivate_right(self, W=None):
        r"""
        Return the right derivative under delta for ``self``.

        EXAMPLES::

            sage: W = Words('12')
            sage: W('12211').delta_derivate_right()
            word: 122
            sage: W('1').delta_derivate_right(Words([1]))
            word: 1
            sage: W('2112').delta_derivate_right()
            word: 12
            sage: W('2211').delta_derivate_right()
            word: 22
            sage: W('112').delta_derivate_right()
            word: 2
            sage: W('11222').delta_derivate_right(Words([1, 2, 3]))
            word: 23
        """
        d = self.delta()
        if len(d) == 0:
            return d
        if W is None:
            W = d.parent()
        if d[-1] != W.alphabet().last():
            d = d[:-1]
        return d

    def phi(self):
        r"""
        Apply the phi function to ``self`` and return the result. This is
        the word obtained by taking the first letter of the words obtained
        by iterating delta on ``self``.

        OUTPUT: a word -- the result of the phi function

        EXAMPLES::

            sage: W = Words([1, 2])
            sage: W([2,2,1,1,2,1,2,2,1,2,2,1,1,2]).phi()
            word: 222222
            sage: W([2,1,2,2,1,2,2,1,2,1]).phi()
            word: 212113
            sage: W().phi()
            word:
            sage: Word([2,1,2,2,1,2,2,1,2,1]).phi()
            word: 212113
            sage: Word([2,3,1,1,2,1,2,3,1,2,2,3,1,2]).phi()
            word: 21215
            sage: Word("aabbabaabaabba").phi()
            word: a22222
            sage: w = Word([2,3,1,1,2,1,2,3,1,2,2,3,1,2])

        See [BL2003]_ and [BDLV2006]_.
        """
        if self.is_empty():
            return self
        v = [self[0]]
        m = self.delta()
        while m.length() > 1:
            v.append(m[0])
            m = m.delta()
        v.append(m[0])
        return Words()(v)

    def phi_inv(self, W=None):
        r"""
        Apply the inverse of the phi function to ``self``.

        INPUT:

        - ``self`` -- a word over the integers
        - ``W`` -- a parent object of words defined over integers

        OUTPUT: a word -- the inverse of the phi function

        EXAMPLES::

            sage: W = Words([1, 2])
            sage: W([2, 2, 2, 2, 1, 2]).phi_inv()
            word: 22112122
            sage: W([2, 2, 2]).phi_inv(Words([2, 3]))
            word: 2233
        """
        if W is None:
            W = self.parent()
        if self.is_empty():
            return W()
        v = self.parent()((self[-1],), check=False)
        for i in range(self.length()-2, -1, -1):
            v = v.delta_inv(W, self[i])
        return v

    def _phi_inv_tab(self, tab):
        r"""
        Specialized version of ``phi_inv()`` for long or incremental words.

        TESTS::

            sage: Words([1, 2])([1, 1, 2, 2])._phi_inv_tab([2])
            word: 12211
        """
        res = self.delta_inv(s=tab[0])
        res = res[1:]
        for i in range(1, len(tab)):
            res = res.delta_inv(s=tab[i])
        return res

    def is_smooth_prefix(self):
        r"""
        Return ``True`` if ``self`` is the prefix of a smooth word, and ``False``
        otherwise.

        Let `A_k = \{1, \ldots ,k\}`, `k \geq 2`. An infinite word `w` in
        `A_k^\omega` is said to be *smooth* if and only if for all positive
        integers `m`, `\Delta^m(w)` is in `A_k^\omega`, where `\Delta(w)` is
        the word obtained from `w` by composing the length of consecutive
        runs of the same letter in `w`. See for instance [BL2003]_ and [BDLV2006]_.

        INPUT:

        - ``self`` -- must be a word over the integers to get something other
          than ``False``

        OUTPUT: boolean; whether ``self`` is a smooth prefix or not

        EXAMPLES::

            sage: W = Words([1, 2])
            sage: W([1, 1, 2, 2, 1, 2, 1, 1]).is_smooth_prefix()
            True
            sage: W([1, 2, 1, 2, 1, 2]).is_smooth_prefix()
            False
        """
        m = self
        W = self.parent()
        while m.length() > 1:
            m = m.delta_derivate_right()
            if not all(a in W.alphabet() for a in m.letters()):
                return False
        return True

    def letters(self):
        r"""
        Return the list of letters that appear in this word, listed in the
        order of first appearance.

        EXAMPLES::

            sage: Word([0,1,1,0,1,0,0,1]).letters()
            [0, 1]
            sage: Word("cacao").letters()
            ['c', 'a', 'o']

        TESTS::

            sage: Word().letters()
            []
        """
        seen = set()
        res = []
        for x in self:
            if x not in seen:
                res.append(x)
                seen.add(x)
        return res

    def standard_factorization(self):
        r"""
        Return the standard factorization of ``self``.

        The *standard factorization* of a word `w` of length greater than
        `1` is the factorization `w = uv` where `v` is the longest proper
        suffix of `w` that is a Lyndon word.

        Note that if `w` is a Lyndon word of length greater than `1` with
        standard factorization `w = uv`, then `u` and `v` are also Lyndon
        words and `u < v`.

        See for instance [CFL1958]_, [Duv1983]_ and [Lot2002]_.

        INPUT:

        - ``self`` -- finite word of length greater than `1`

        OUTPUT: `2`-tuple `(u, v)`

        EXAMPLES::

            sage: Words('01')('0010110011').standard_factorization()
            (word: 001011, word: 0011)
            sage: Words('123')('1223312').standard_factorization()
            (word: 12233, word: 12)
            sage: Word([3,2,1]).standard_factorization()
            (word: 32, word: 1)

        ::

            sage: w = Word('0010110011',alphabet='01')
            sage: w.standard_factorization()
            (word: 001011, word: 0011)
            sage: w = Word('0010110011',alphabet='10')
            sage: w.standard_factorization()
            (word: 001011001, word: 1)
            sage: w = Word('1223312',alphabet='123')
            sage: w.standard_factorization()
            (word: 12233, word: 12)

        TESTS::

            sage: w = Word()
            sage: w.standard_factorization()
            Traceback (most recent call last):
            ...
            ValueError: Standard factorization not defined on words of
            length less than 2
            sage: w = Word('a')
            sage: w.standard_factorization()
            Traceback (most recent call last):
            ...
            ValueError: Standard factorization not defined on words of
            length less than 2
        """
        selflen = self.length()
        if selflen < 2:
            raise ValueError("Standard factorization not defined on"
                             " words of length less than 2")
        for l in range(1, selflen):
            suff = self[l:]
            if suff.is_lyndon():
                return self[:l], suff

    def apply_permutation_to_positions(self, permutation):
        r"""
        Return the word obtained by permuting the positions of the letters
        in ``self`` according to the permutation ``permutation``.

        EXAMPLES::

            sage: w = Words('abcd')('abcd')
            sage: w.apply_permutation_to_positions([2,1,4,3])
            word: badc
            sage: u = Words('dabc')('abcd')
            sage: u.apply_permutation_to_positions([2,1,4,3])
            word: badc
            sage: w.apply_permutation_to_positions(Permutation([2,1,4,3]))
            word: badc
            sage: w.apply_permutation_to_positions(PermutationGroupElement([2,1,4,3]))  # needs sage.groups
            word: badc
            sage: Word([1,2,3,4]).apply_permutation_to_positions([3,4,2,1])
            word: 3421
        """
        from sage.combinat.permutation import Permutation
        if not isinstance(permutation, Permutation):
            if isinstance(permutation, PermutationGroupElement):
                permutation = Permutation(permutation.domain())
            else:
                permutation = Permutation(permutation)
        return self.parent()(permutation.action(self), check=False)

    def apply_permutation_to_letters(self, permutation):
        r"""
        Return the word obtained by applying the permutation
        ``permutation`` of the alphabet of ``self`` to each letter of
        ``self``.

        EXAMPLES::

            sage: w = Words('abcd')('abcd')
            sage: p = [2,1,4,3]
            sage: w.apply_permutation_to_letters(p)
            word: badc
            sage: u = Words('dabc')('abcd')
            sage: u.apply_permutation_to_letters(p)
            word: dcba
            sage: w.apply_permutation_to_letters(Permutation(p))
            word: badc
            sage: w.apply_permutation_to_letters(PermutationGroupElement(p))            # needs sage.groups
            word: badc
        """
        from sage.combinat.permutation import Permutation
        if not isinstance(permutation, Permutation):
            if isinstance(permutation, PermutationGroupElement):
                permutation = Permutation(permutation.domain())
            else:
                permutation = Permutation(permutation)
        alphabet = self.parent().alphabet()
        morphism = dict(zip(alphabet, permutation.action(alphabet)))
        return self.apply_morphism(morphism)

    def colored_vector(self, x=0, y=0, width='default', height=1, cmap='hsv', thickness=1, label=None):
        r"""
        Return a vector (Graphics object) illustrating ``self``. Each letter
        is represented by a coloured rectangle.

        If the parent of ``self`` is a class of words over a finite alphabet,
        then each letter in the alphabet is assigned a unique colour, and
        this colour will be the same every time this method is called. This
        is especially useful when plotting and comparing words defined on
        the same alphabet.

        If the alphabet is infinite, then the letters appearing in the word
        are used as the alphabet.

        INPUT:

        - ``x`` -- (default: ``0``) bottom left x-coordinate of the vector
        - ``y`` -- (default: ``0``) bottom left y-coordinate of the vector
        - ``width`` -- (default: ``'default'``) width of the vector. By default,
          the width is the length of ``self``.
        - ``height`` -- (default: ``1``) height of the vector
        - ``thickness`` -- (default: ``1``) thickness of the contour
        - ``cmap`` -- (default: ``'hsv'``) color map; for available color map names
          type: ``import matplotlib.cm; list(matplotlib.cm.datad)``
        - ``label`` -- string (default: ``None``); a label to add on the colored vector

        OUTPUT: Graphics

        EXAMPLES::

            sage: # needs sage.plot
            sage: Word(range(20)).colored_vector()
            Graphics object consisting of 21 graphics primitives
            sage: Word(range(100)).colored_vector(0,0,10,1)
            Graphics object consisting of 101 graphics primitives
            sage: Words(range(100))(range(10)).colored_vector()
            Graphics object consisting of 11 graphics primitives
            sage: w = Word('abbabaab')
            sage: w.colored_vector()
            Graphics object consisting of 9 graphics primitives
            sage: w.colored_vector(cmap='autumn')
            Graphics object consisting of 9 graphics primitives
            sage: Word(range(20)).colored_vector(label='Rainbow')
            Graphics object consisting of 23 graphics primitives

        When two words are defined under the same parent, same letters are
        mapped to same colors::

            sage: W = Words(range(20))
            sage: w = W(range(20))
            sage: y = W(range(10,20))
            sage: y.colored_vector(y=1, x=10) + w.colored_vector()                      # needs sage.plot
            Graphics object consisting of 32 graphics primitives

        TESTS:

        The empty word::

            sage: Word().colored_vector()                                               # needs sage.plot
            Graphics object consisting of 1 graphics primitive
            sage: Word().colored_vector(label='empty')                                  # needs sage.plot
            Graphics object consisting of 3 graphics primitives

        Unknown cmap::

            sage: Word(range(100)).colored_vector(cmap='jolies')                        # needs sage.plot
            Traceback (most recent call last):
            ...
            RuntimeError: Color map jolies not known
            sage: Word(range(100)).colored_vector(cmap='__doc__')                       # needs sage.plot
            Traceback (most recent call last):
            ...
            RuntimeError: Color map __doc__ not known
        """
        # Recognize the color map
        import matplotlib.cm as cm
        from matplotlib.colors import LinearSegmentedColormap as C
        key_error = False
        try:
            mpl_cmap = cm.__dict__[cmap]
        except KeyError:
            key_error = True

        if key_error or not isinstance(mpl_cmap, C):
            possibilities = ', '.join(str(x) for x, val in cm.__dict__.items()
                                      if isinstance(val, C))
            import sage.misc.verbose
            sage.misc.verbose.verbose("The possible color maps include: %s" % possibilities, level=0)
            raise RuntimeError("Color map %s not known" % cmap)

        # Drawing the colored vector...
        from sage.plot.line import line
        from sage.plot.polygon import polygon
        from sage.plot.text import text

        # The default width of the vector
        if width == 'default':
            width = self.length()

        # The black frame of the vector
        ymax = y + height
        L = [(x, y), (x+width, y), (x+width, ymax), (x, ymax), (x, y)]
        rep = line(L, rgbcolor=(0, 0, 0), thickness=thickness)

        # The label
        if label is not None:
            hl = height/2.0  # height of the label rectangle
            ymax2 = ymax + hl
            rep += text(str(label), (x+width/2.0, ymax + hl/2.0), rgbcolor=(1, 0, 0))
            L = [(x, ymax), (x+width, ymax), (x+width, ymax2), (x, ymax2), (x, ymax)]
            rep += line(L, rgbcolor=(0, 0, 0), thickness=thickness)

        # base : the width of each rectangle
        base = width / float(self.length()) if not self.is_empty() else None

        # A colored rectangle for each letter
        dim = self.parent().alphabet().cardinality()
        if dim is Infinity:
            ordered_alphabet = sorted(self.letters(),
                                      key=self.parent().sortkey_letters)
            dim = float(len(ordered_alphabet))
        else:
            ordered_alphabet = self.parent().alphabet()
            dim = float(self.parent().alphabet().cardinality())
        letter_to_integer_dict = {a: i for i, a in
                enumerate(ordered_alphabet)}
        xp = x
        for a in self:
            i = letter_to_integer_dict[a]
            xq = xp + base
            L = [(xp, y), (xq, y), (xq, ymax), (xp, ymax)]
            rgbcolor = mpl_cmap(i / dim)[:3]
            rep += polygon(L, rgbcolor=rgbcolor)
            xp = xq
        rep.axes(False)
        return rep

    def is_square(self):
        r"""
        Return ``True`` if ``self`` is a square, and ``False`` otherwise.

        EXAMPLES::

            sage: Word([1,0,0,1]).is_square()
            False
            sage: Word('1212').is_square()
            True
            sage: Word('1213').is_square()
            False
            sage: Word('12123').is_square()
            False
            sage: Word().is_square()
            True
        """
        if self.length() % 2:
            return False
        else:
            l = self.length() // 2
            return self[:l] == self[l:]

    def is_square_free(self):
        r"""
        Return ``True`` if ``self`` does not contain squares, and ``False``
        otherwise.

        EXAMPLES::

            sage: Word('12312').is_square_free()
            True
            sage: Word('31212').is_square_free()
            False
            sage: Word().is_square_free()
            True

        TESTS:

        We make sure that :issue:`8490` is fixed::

            sage: Word('11').is_square_free()
            False
            sage: Word('211').is_square_free()
            False
            sage: Word('3211').is_square_free()
            False
        """
        from sage.combinat.words.suffix_trees import DecoratedSuffixTree
        T = DecoratedSuffixTree(self)
        return T.square_vocabulary() == [(0, 0)]

    def squares(self):
        r"""
        Return a set of all distinct squares of ``self``.

        EXAMPLES::

            sage: sorted(Word('cacao').squares())
            [word: , word: caca]
            sage: sorted(Word('1111').squares())
            [word: , word: 11, word: 1111]
            sage: w = Word('00110011010')
            sage: sorted(w.squares())
            [word: , word: 00, word: 00110011, word: 01100110, word: 1010, word: 11]
        """
        from sage.combinat.words.suffix_trees import DecoratedSuffixTree
        T = DecoratedSuffixTree(self)
        return set(T.square_vocabulary(output='word'))

    def is_cube(self):
        r"""
        Return ``True`` if ``self`` is a cube, and ``False`` otherwise.

        EXAMPLES::

            sage: Word('012012012').is_cube()
            True
            sage: Word('01010101').is_cube()
            False
            sage: Word().is_cube()
            True
            sage: Word('012012').is_cube()
            False
        """
        if self.length() % 3 != 0:
            return False
        l = self.length() // 3
        return self[:l] == self[l:2*l] == self[2*l:]

    def is_cube_free(self):
        r"""
        Return ``True`` if ``self`` does not contain cubes, and ``False`` otherwise.

        EXAMPLES::

            sage: Word('12312').is_cube_free()
            True
            sage: Word('32221').is_cube_free()
            False
            sage: Word().is_cube_free()
            True

        TESTS:

        We make sure that :issue:`8490` is fixed::

            sage: Word('111').is_cube_free()
            False
            sage: Word('2111').is_cube_free()
            False
            sage: Word('32111').is_cube_free()
            False
        """
        L = self.length()
        if L < 3:
            return True
        for start in range(L - 2):
            for end in range(start + 3, L + 1, 3):
                if self[start:end].is_cube():
                    return False
        return True

    def to_monoid_element(self):
        """
        Return ``self`` as an element of the free monoid with the same alphabet
        as ``self``.

        EXAMPLES::

            sage: w = Word('aabb')
            sage: w.to_monoid_element()
            a^2*b^2
            sage: W = Words('abc')
            sage: w = W(w)
            sage: w.to_monoid_element()
            a^2*b^2

        TESTS:

        Check that ``w == w.to_monoid_element().to_word()``::

            sage: all(w.to_monoid_element().to_word() == w for i in range(6) for w in Words('abc', i))
            True
        """
        from sage.monoids.free_monoid import FreeMonoid
        try:
            l = list(self.parent().alphabet())
        except AttributeError:
            l = self.letters()
        M = FreeMonoid(len(l), l)
        return M(self)

    def is_christoffel(self):
        r"""
        Return ``True`` if ``self`` is a Christoffel word, and ``False`` otherwise.

        The *Christoffel word* of slope `p/q` is obtained from the Cayley
        graph of `\ZZ/(p+q)\ZZ` with generator `q` as follows. If `u
        \rightarrow v` is an edge in the Cayley graph, then, `v = u + p
        \mod{p+q}`. Let `a`,`b` be the alphabet of `w`. Label the edge
        `u \rightarrow v` by `a` if `u < v` and `b` otherwise. The Christoffel
        word is the word obtained by reading the edge labels along the cycle
        beginning from `0`.

        Equivalently, `w` is a Christoffel word iff `w` is a symmetric
        non-empty word and `w[1:n-1]` is a palindrome.

        See for instance [Ber2007]_ and [BLRS2009]_.

        INPUT:

        - ``self`` -- word

        OUTPUT: boolean; ``True`` if ``self`` is a Christoffel word,
        ``False`` otherwise

        EXAMPLES::

            sage: Word('00100101').is_christoffel()
            True
            sage: Word('aab').is_christoffel()
            True
            sage: Word().is_christoffel()
            False
            sage: Word('123123123').is_christoffel()
            False
            sage: Word('00100').is_christoffel()
            False
            sage: Word('0').is_christoffel()
            True

        TESTS::

            sage: words.LowerChristoffelWord(5,4).is_christoffel()
            True
            sage: words.UpperChristoffelWord(5,4).is_christoffel()
            True
            sage: Word('aaaaaaaaa').is_christoffel()
            False
        """
        if len(self) == 0 or len(self.letters()) > 2 or (self.is_palindrome() and len(self) > 1):
            return False
        elif self.is_symmetric() and self[1:len(self) - 1].is_palindrome():
            return True
        else:
            return False

    def minimal_conjugate(self):
        r"""
        Return the lexicographically minimal conjugate of this word (see
        :wikipedia:`Lexicographically_minimal_string_rotation`).

        EXAMPLES::

            sage: Word('213').minimal_conjugate()
            word: 132
            sage: Word('11').minimal_conjugate()
            word: 11
            sage: Word('12112').minimal_conjugate()
            word: 11212
            sage: Word('211').minimal_conjugate()
            word: 112
            sage: Word('211211211').minimal_conjugate()
            word: 112112112

        TESTS::

            sage: Word().minimal_conjugate()
            word:
        """
        if not self:
            return self
        p = self.primitive()
        q = self.length() // p.length()
        end = 0
        for factor in (p ** 2).lyndon_factorization():
            end += factor.length()
            if end >= p.length():
                return factor ** q


class CallableFromListOfWords(tuple):
    r"""
    A class to create a callable from a list of words. The concatenation of
    a list of words is obtained by creating a word from this callable.
    """
    def __new__(cls, words):
        r"""
        TESTS::

            sage: from sage.combinat.words.finite_word import CallableFromListOfWords
            sage: w,u,x = Word([1,2,3]),Word([4,5]),Word([6,7,8])
            sage: f = CallableFromListOfWords([w,u,x]); f
            (word: 123, word: 45, word: 678)
            sage: f == loads(dumps(f))
            True
        """
        l = []
        for w in words:
            from .word_infinite_datatypes import WordDatatype_callable
            if isinstance(w, WordDatatype_callable) and \
                    isinstance(w._func, CallableFromListOfWords):
                l.extend(w._func)
            else:
                l.append(w)
        return tuple.__new__(cls, l)

    def __call__(self, i):
        r"""
        Return the character at position ``i``.

        TESTS::

            sage: from sage.combinat.words.finite_word import CallableFromListOfWords
            sage: w,u,x = Word([1,2,3]),Word([4,5]),Word([6,7,8])
            sage: f = CallableFromListOfWords([w,u,x])
            sage: [f(i) for i in range(8)]
            [1, 2, 3, 4, 5, 6, 7, 8]
        """
        j = i
        for c in self:
            if j < c.length():
                return c[j]
            j -= c.length()
        raise IndexError(f"index (={i}) out of range")


class Factorization(list):
    r"""
    A list subclass having a nicer representation for factorization of words.

    TESTS::

        sage: f = sage.combinat.words.finite_word.Factorization()
        sage: f == loads(dumps(f))
        True
    """
    def __repr__(self):
        r"""
        Return a string representation of the object.

        TESTS::

            sage: sage.combinat.words.finite_word.Factorization()
            ()
            sage: sage.combinat.words.finite_word.Factorization([Word('ab'), Word('ba')])
            (ab, ba)
        """
        return '(%s)' % ', '.join(w.string_rep() for w in self)


#######################################################################

def evaluation_dict(w):
    r"""
    Return a dictionary keyed by the letters occurring in ``w`` with
    values the number of occurrences of the letter.

    INPUT:

    - ``w`` -- a word

    TESTS::

        sage: from sage.combinat.words.finite_word import evaluation_dict
        sage: evaluation_dict([2,1,4,2,3,4,2])
        {1: 1, 2: 3, 3: 1, 4: 2}
        sage: evaluation_dict('badbcdb')
        {'a': 1, 'b': 3, 'c': 1, 'd': 2}
        sage: evaluation_dict([])
        {}

    ::

        sage: evaluation_dict('1213121') # keys appear in random order
        {'1': 4, '2': 2, '3': 1}
    """
    d = defaultdict(int)
    for a in w:
        d[a] += 1
    return dict(d)


def word_to_ordered_set_partition(w):
    r"""
    Return the ordered set partition corresponding to a finite
    word `w`.

    If `w` is a finite word of length `n`, then the corresponding
    ordered set partition is an ordered set partition
    `(P_1, P_2, \ldots, P_k)` of `\{1, 2, \ldots, n\}`, where
    each block `P_i` is the set of positions at which the `i`-th
    smallest letter occurring in `w` occurs in `w`.
    (Positions are `1`-based.)

    This is the same functionality that
    :meth:`~sage.combinat.words.finite_word.FiniteWord_class.to_ordered_set_partition`
    provides, but without the wrapping: The input `w` can be given as
    a list or tuple, not necessarily as a word; and the output is
    returned as a list of lists (which are the blocks of the ordered
    set partition in increasing order), not as an ordered set partition.

    EXAMPLES::

        sage: from sage.combinat.words.finite_word import word_to_ordered_set_partition
        sage: word_to_ordered_set_partition([3, 6, 3, 1])
        [[4], [1, 3], [2]]
        sage: word_to_ordered_set_partition((1, 3, 3, 7))
        [[1], [2, 3], [4]]
        sage: word_to_ordered_set_partition("noob")
        [[4], [1], [2, 3]]
        sage: word_to_ordered_set_partition(Word("hell"))
        [[2], [1], [3, 4]]
        sage: word_to_ordered_set_partition([1])
        [[1]]
        sage: word_to_ordered_set_partition([])
        []
    """
    vals = sorted(set(w))
    dc = {val: i for i, val in enumerate(vals)}
    P = [[] for _ in vals]
    for i, val in enumerate(w):
        P[dc[val]].append(i + 1)
    return P
