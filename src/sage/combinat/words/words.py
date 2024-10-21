# sage_setup: distribution = sagemath-combinat
r"""
Set of words

To define a new class of words, please refer to the documentation file:
sage/combinat/words/notes/word_inheritance_howto.rst

AUTHORS:

    - Franco Saliola (2008-12-17): merged into sage
    - Sebastien Labbe (2008-12-17): merged into sage
    - Arnaud Bergeron (2008-12-17): merged into sage
    - Sebastien Labbe (2009-07-21): Improved morphism iterator (:issue:`6571`).
    - Vincent Delecroix (2015): classes simplifications (:issue:`19619`)

EXAMPLES::

    sage: Words()
    Finite and infinite words over Set of Python objects of class 'object'
    sage: Words(4)
    Finite and infinite words over {1, 2, 3, 4}
    sage: Words(4,5)
    Words of length 5 over {1, 2, 3, 4}

    sage: FiniteWords('ab')
    Finite words over {'a', 'b'}
    sage: InfiniteWords('natural numbers')
    Infinite words over Non negative integers
"""
# ****************************************************************************
#       Copyright (C) 2008 Arnaud Bergeron <abergeron@gmail.com>,
#                          Sébastien Labbé <slabqc@gmail.com>,
#                          Franco Saliola <saliola@gmail.com>
#                     2015 Vincent Delecroix <20100.delecroix@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import itertools
from collections.abc import Iterable

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.persist import register_unpickle_override

from sage.structure.parent import Parent

from sage.categories.sets_cat import Sets

from sage.combinat.combinat import CombinatorialObject
from sage.structure.list_clone import ClonableElement
from sage.combinat.words.alphabet import build_alphabet

from sage.rings.infinity import Infinity
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ


def Words(alphabet=None, length=None, finite=True, infinite=True):
    """
    Return the combinatorial class of words of length k over an alphabet.

    EXAMPLES::

        sage: Words()
        Finite and infinite words over Set of Python objects of class 'object'
        sage: Words(length=7)
        Words of length 7 over Set of Python objects of class 'object'
        sage: Words(5)
        Finite and infinite words over {1, 2, 3, 4, 5}
        sage: Words(5, 3)
        Words of length 3 over {1, 2, 3, 4, 5}
        sage: Words(5, infinite=False)
        Finite words over {1, 2, 3, 4, 5}
        sage: Words(5, finite=False)
        Infinite words over {1, 2, 3, 4, 5}
        sage: Words('ab')
        Finite and infinite words over {'a', 'b'}
        sage: Words('ab', 2)
        Words of length 2 over {'a', 'b'}
        sage: Words('ab', infinite=False)
        Finite words over {'a', 'b'}
        sage: Words('ab', finite=False)
        Infinite words over {'a', 'b'}
        sage: Words('positive integers', finite=False)
        Infinite words over Positive integers
        sage: Words('natural numbers')
        Finite and infinite words over Non negative integers
    """
    if isinstance(alphabet, FiniteWords) or \
       isinstance(alphabet, InfiniteWords) or \
       isinstance(alphabet, FiniteOrInfiniteWords) or \
       isinstance(alphabet, Words_n):
        return alphabet

    if length is None:
        if finite and infinite:
            return FiniteOrInfiniteWords(alphabet)
        elif finite:
            return FiniteWords(alphabet)
        elif infinite:
            return InfiniteWords(alphabet)
        else:
            raise ValueError("either finite or infinite must be True")

    elif isinstance(length, (int, Integer)):
        return Words_n(FiniteWords(alphabet), length)

    raise ValueError("do not know how to make a combinatorial class of words from your input")


class AbstractLanguage(Parent):
    r"""
    Abstract base class.

    This is *not* to be used by any means. This class gather previous features
    of set of words (prior to :issue:`19619`). In the future that class might
    simply disappear or become a common base class for all languages. In the
    latter case, its name would possibly change to ``Language``.
    """
    def __init__(self, alphabet=None, category=None):
        r"""
        INPUT:

        - ``alphabet`` -- the underlying alphabet

        TESTS::

            sage: loads(dumps(FiniteWords('ab'))) == FiniteWords('ab')
            True
            sage: loads(dumps(InfiniteWords('ab'))) == InfiniteWords('ab')
            True

            sage: Words('abc').sortkey_letters
            <bound method AbstractLanguage._sortkey_trivial of ...>
            sage: Words('bac').sortkey_letters
            <bound method AbstractLanguage._sortkey_letters of ...>
        """
        if isinstance(alphabet, (int, Integer)):
            from sage.sets.integer_range import IntegerRange
            alphabet = IntegerRange(1, alphabet + 1)
        elif (alphabet == "integers" or
              alphabet == "positive integers" or
              alphabet == "natural numbers"):
            alphabet = build_alphabet(name=alphabet)
        else:
            alphabet = build_alphabet(alphabet)

        self._alphabet = alphabet

        # Default sorting key: use rank()
        self.sortkey_letters = self._sortkey_letters

        # Check if we should use the trivial sorting key
        N = alphabet.cardinality()
        if N == Infinity:
            self.sortkey_letters = self._sortkey_trivial
        elif N < 36:
            try:
                if all(alphabet.unrank(i) > alphabet.unrank(j)
                       for i in range(N) for j in range(i)):
                    self.sortkey_letters = self._sortkey_trivial
            except TypeError:
                pass

        if category is None:
            category = Sets()

        Parent.__init__(self, category=category)

    def alphabet(self):
        r"""
        EXAMPLES::

            sage: Words(NN).alphabet()
            Non negative integer semiring

            sage: InfiniteWords([1,2,3]).alphabet()
            {1, 2, 3}
            sage: InfiniteWords('ab').alphabet()
            {'a', 'b'}

            sage: FiniteWords([1,2,3]).alphabet()
            {1, 2, 3}
            sage: FiniteWords().alphabet()
            Set of Python objects of class 'object'
        """
        return self._alphabet

    def identity_morphism(self):
        r"""
        Return the identity morphism from ``self`` to itself.

        EXAMPLES::

            sage: W = Words('ab')
            sage: W.identity_morphism()
            WordMorphism: a->a, b->b

        ::

            sage: W = Words(range(3))
            sage: W.identity_morphism()
            WordMorphism: 0->0, 1->1, 2->2

        There is no support yet for infinite alphabet::

            sage: W = Words(alphabet=Alphabet(name='NN'))
            sage: W
            Finite and infinite words over Non negative integers
            sage: W.identity_morphism()
            Traceback (most recent call last):
            ...
            NotImplementedError: size of alphabet must be finite
        """
        if self.alphabet().cardinality() not in ZZ:
            raise NotImplementedError('size of alphabet must be finite')
        from sage.combinat.words.morphism import WordMorphism
        return WordMorphism({a: a for a in self.alphabet()})

    def _check(self, w, length=40):
        r"""
        Check that the first length elements are actually in the alphabet.

        INPUT:

        - ``w`` -- word

        - ``length`` -- integer (default: `40`)

        EXAMPLES::

            sage: W = FiniteWords(['a','b','c'])
            sage: W._check('abcabc') is None
            True
            sage: W._check('abcabcd')
            Traceback (most recent call last):
            ...
            ValueError: d not in alphabet
            sage: W._check('abcabc'*10+'z') is None
            True
            sage: W._check('abcabc'*10+'z', length=80)
            Traceback (most recent call last):
            ...
            ValueError: z not in alphabet
        """
        stop = None if length is None else int(length)
        for a in itertools.islice(w, stop):
            if a not in self.alphabet():
                raise ValueError("%s not in alphabet" % a)

    def _sortkey_trivial(self, letter1):
        """
        Trivial function, used to sort the letters by their names.

        INPUT:

        - ``letter1`` -- a letter in the alphabet

        EXAMPLES::

            sage: W = FiniteWords('ade')
            sage: W.sortkey_letters('d')  # indirect doctest
            'd'
        """
        return letter1

    def _sortkey_letters(self, letter1):
        r"""
        Return the default value used to sort the letters.

        INPUT:

        - ``letter1`` -- a letter in the alphabet

        EXAMPLES::

            sage: W = FiniteWords('woa')
            sage: W.sortkey_letters('w')  # indirect doctest
            0
            sage: W.sortkey_letters('o')  # indirect doctest
            1
            sage: W.sortkey_letters('a')  # indirect doctest
            2

        TESTS::

            sage: assert W.sortkey_letters == W._sortkey_letters
        """
        rk = self.alphabet().rank
        return rk(letter1)

    def __eq__(self, other) -> bool:
        r"""
        TESTS::

            sage: FiniteWords() == FiniteWords()
            True
            sage: FiniteWords() == InfiniteWords()
            False
            sage: InfiniteWords() == Words()
            False
            sage: FiniteWords([0,1]) == FiniteWords([0,1,2,3])
            False
        """
        return self is other or (type(self) is type(other) and
                                 self.alphabet() == other.alphabet())

    def __ne__(self, other) -> bool:
        r"""
        TESTS::

            sage: InfiniteWords() != InfiniteWords()
            False
            sage: FiniteWords() != Words()
            True
            sage: Words('ab') != Words('ab')
            False
        """
        return not (self == other)


class FiniteWords(AbstractLanguage):
    r"""
    The set of finite words over a fixed alphabet.

    EXAMPLES::

        sage: W = FiniteWords('ab')
        sage: W
        Finite words over {'a', 'b'}
    """
    def cardinality(self):
        r"""
        Return the cardinality of this set.

        EXAMPLES::

            sage: FiniteWords('').cardinality()
            1
            sage: FiniteWords('a').cardinality()
            +Infinity
        """
        if not self.alphabet():
            return ZZ.one()
        return Infinity

    def __hash__(self):
        r"""
        TESTS::

            sage: hash(FiniteWords('ab')) # random
            12
        """
        return hash(self.alphabet()) ^ hash('finite words')

    @cached_method
    def shift(self):
        r"""
        Return the set of infinite words on the same alphabet.

        EXAMPLES::

            sage: FiniteWords('ab').shift()
            Infinite words over {'a', 'b'}
        """
        return InfiniteWords(self.alphabet())

    def factors(self):
        r"""
        Return itself.

        EXAMPLES::

            sage: FiniteWords('ab').factors()
            Finite words over {'a', 'b'}
        """
        return self

    @lazy_attribute
    def _element_classes(self):
        r"""
        Return a dictionary that gives the class of the element of ``self``.

        The word may be finite, infinite or of unknown length.
        Its data may be string, list, tuple, a callable or an iterable.
        For callable and iterable, the data may be cached.

        EXAMPLES:

        Once you get the class, it can be used to create a new word::

            sage: W = FiniteWords([0,1,2])
            sage: L = [0,1,0] * 100
            sage: cls = W._element_classes['list']
            sage: w = cls(W, L)
            sage: type(w)
            <class 'sage.combinat.words.word.FiniteWord_list'>
            sage: w
            word: 0100100100100100100100100100100100100100...
            sage: w.parent()
            Finite words over {0, 1, 2}

        TESTS::

            sage: d = FiniteWords()._element_classes
            sage: type(d)
            <class 'dict'>
            sage: len(d)
            7
            sage: e = FiniteWords('abcdefg')._element_classes
            sage: d == e
            True
        """
        import sage.combinat.words.word as word
        classes = {
            'list': word.FiniteWord_list,
            'str': word.FiniteWord_str,
            'tuple': word.FiniteWord_tuple,
            'callable_with_caching': word.FiniteWord_callable_with_caching,
            'callable': word.FiniteWord_callable,
            'iter_with_caching': word.FiniteWord_iter_with_caching,
            'iter': word.FiniteWord_iter}

        # test whether or not we can use the class Finiteword_char
        if (self.alphabet().cardinality() <= 256 and
                all(isinstance(i, (int, Integer)) and
                    0 <= i < 256 for i in self.alphabet())):
            L = self.alphabet().list()
            key = self.sortkey_letters
            if (all(L[i] < L[i + 1] for i in range(len(L) - 1)) and
                    all(key(L[i]) < key(L[i + 1]) for i in range(len(L) - 1))):
                classes['char'] = word.FiniteWord_char

        return classes

    def _word_from_word(self, data):
        r"""
        Return a word from a word.

        The data is assumed to be ok, no check is performed.

        INPUT:

        - ``data`` -- word

        EXAMPLES::

            sage: W = FiniteWords([0,1,2])
            sage: w = W([0,1,2,0,1,2])
            sage: z = W._word_from_word(w)
            sage: z
            word: 012012
            sage: w is z
            True
        """
        ####################
        # If `data` is already a word and if its parent is self, then
        # return `data`.
        ###########################
        if data.parent() is self or data.parent() == self:
            return data

        ###########################
        # Otherwise, if self is not the parent of `data`, then we try to
        # recover the data, the length and the datatype of the input `data`
        # To minimize the impact of the import, we do it only at the time there
        # are needed
        ###########################
        from sage.combinat.words.word_char import WordDatatype_char
        if isinstance(data, WordDatatype_char):
            data = list(data)
            if 'char' in self._element_classes:
                return self._element_classes['char'](self, data)
            else:
                return self._element_classes['list'](self, data)

        from sage.combinat.words.word_datatypes import (WordDatatype_str,
                                                        WordDatatype_list,
                                                        WordDatatype_tuple)
        if isinstance(data, WordDatatype_str):
            return self._element_classes['str'](self, data._data)
        if isinstance(data, WordDatatype_tuple):
            return self._element_classes['tuple'](self, data._data)
        if isinstance(data, WordDatatype_list):
            return self._element_classes['list'](self, data._data)

        from sage.combinat.words.word_infinite_datatypes import \
            (WordDatatype_callable, WordDatatype_iter)
        if isinstance(data, WordDatatype_callable):
            length = data.length()
            data = data._func
            return self._word_from_callable(data, length, caching=False)
        if isinstance(data, WordDatatype_iter):
            length = data.length()
            data = iter(data)
            return self._word_from_iter(data, length, caching=False)

        raise TypeError("any instance of Word_class must be an instance of WordDatatype")

    def _word_from_callable(self, data, length, caching=True):
        r"""
        Return a word represented by a callable.

        The data is assumed to be ok, no check is performed.

        INPUT:

        - ``data`` -- callable
        - ``length`` -- integer or ``None`` or ``'infinite'`` or ``Infinity``
        - ``caching`` -- boolean (default: ``True``); whether to keep a cache
          of the letters computed by the callable

        EXAMPLES::

            sage: W = FiniteWords([0,1,2])
            sage: f = lambda n : n % 3
            sage: W._word_from_callable(f, 100)
            word: 0120120120120120120120120120120120120120...
        """
        wc = '_with_caching' if caching else ""
        if length not in ZZ or length < 0:
            raise ValueError("not a correct value for length (%s)" % length)
        return self._element_classes['callable' + wc](self, data, length)

    def _word_from_iter(self, data, length=None, caching=True):
        r"""
        Return a word represented by an iterator.

        The data is assumed to be ok, no check is performed.

        INPUT:

        - ``data`` -- iterable

        - ``length`` -- (optional) integer

        - ``caching`` -- boolean (default: ``True``); whether to keep a cache
          of the letters computed by the iterator

        EXAMPLES::

            sage: W = FiniteWords([0,1,2])
            sage: W._word_from_iter(iter([1]*10))
            word: 1111111111
            sage: W._word_from_iter(iter([1]*10), 5)
            word: 11111
        """
        wc = '_with_caching' if caching else ""
        if length is None or length == "finite":
            length = "finite"
        elif length not in ZZ or length < 0:
            raise ValueError("not a correct value for length (%s)" % length)
        return self._element_classes['iter' + wc](self, data, length)

    def __call__(self, data=None, length=None, datatype=None, caching=True, check=True):
        r"""
        Construct a new word object with parent ``self``.

        INPUT:

        - ``data`` -- (default: ``None``) list, string, tuple, iterator, ``None``
          (shorthand for []), or a callable defined on [0,1,...,length]

        - ``length`` -- integer (default: ``None``); only used if the data is
          an iterator or a callable. It determines the length of the word

        - ``datatype`` -- (default: ``None``) ``None``, "char", "list", "str",
          "tuple", "iter", "callable" or "pickled_function"; if ``None``, then
          the function tries to guess this from the data

        - ``caching`` -- boolean (default: ``True``); whether to keep a cache
          of the letters computed by an iterator or callable

        - ``check`` -- boolean (default: ``True``); whether to check if
          the 40 first letters are in the parent alphabet. This is a
          check done to test for small programming errors. Since we also
          support infinite words, we cannot really implement a more
          accurate check.

        .. NOTE::

           The check makes this method about 10 times slower (20µs instead
           of 2µs), so make sure to set it to False if you know the
           alphabet is OK. Fast creation (about 1µs) of a word can be
           done using the class directly (see :meth:`_element_classes`).

        .. WARNING::

           Be careful when defining words using callables and iterators. It
           appears that islice does not pickle correctly causing various errors
           when reloading. Also, most iterators do not support copying and
           should not support pickling by extension.

        EXAMPLES::

            sage: W = FiniteWords()

        Empty word::

            sage: W()
            word:

        Word with string::

            sage: W("abbabaab")
            word: abbabaab

        Word with string constructed from other types::

            sage: W([0,1,1,0,1,0,0,1], datatype='str')
            word: 01101001
            sage: W((0,1,1,0,1,0,0,1), datatype='str')
            word: 01101001

        Word with list::

            sage: W([0,1,1,0,1,0,0,1])
            word: 01101001

        Word with list constructed from other types::

            sage: W("01101001", datatype='list')
            word: 01101001
            sage: W((0,1,1,0,1,0,0,1), datatype='list')
            word: 01101001

        Word with tuple::

            sage: W((0,1,1,0,1,0,0,1))
            word: 01101001

        Word with tuple constructed from other types::

            sage: W([0,1,1,0,1,0,0,1], datatype='tuple')
            word: 01101001
            sage: W("01101001", datatype='str')
            word: 01101001

        Word with iterator::

            sage: from itertools import count
            sage: W(count(), length=10)
            word: 0123456789
            sage: W(iter("abbabaab"))
            word: abbabaab

        Word with function (a 'callable')::

            sage: f = lambda n : add(Integer(n).digits(2)) % 2
            sage: W(f, length=12)
            word: 011010011001
            sage: FiniteWords([0,1,2])(f, length=12)
            word: 011010011001

        Word over a string with a parent::

            sage: w = FiniteWords('abc')("abbabaab"); w
            word: abbabaab
            sage: w.parent()
            Finite words over {'a', 'b', 'c'}

        The fourty first letters of the word are checked if they are in the
        parent alphabet::

            sage: FiniteWords("ab")("abca")
            Traceback (most recent call last):
            ...
            ValueError: c not in alphabet
            sage: FiniteWords("ab")("abca", check=False)
            word: abca

        The default parent is the combinatorial class of all words::

            sage: w = Word("abbabaab"); w
            word: abbabaab
            sage: w.parent()
            Finite words over Set of Python objects of class 'object'

        Creation of a word from a word::

            sage: FiniteWords([0,1,2,3])(FiniteWords([2,3])([2,2,2,3,3,2]))
            word: 222332
            sage: _.parent()
            Finite words over {0, 1, 2, 3}

        ::

            sage: FiniteWords([3,2,1])(FiniteWords([2,3])([2,2,2,3,3,2]))
            word: 222332
            sage: _.parent()
            Finite words over {3, 2, 1}

        Construction of a word from a word when the parents are the same::

            sage: W = FiniteWords()
            sage: w = W(range(8))
            sage: z = W(w)
            sage: w is z
            True

        Construction of a word path from a finite word::

            sage: W = FiniteWords('abcd')
            sage: P = WordPaths('abcd')                                                 # needs sage.modules
            sage: w = W('aaab')
            sage: P(w)                                                                  # needs sage.modules
            Path: aaab

        Construction of a word path from a Christoffel word::

            sage: w = words.ChristoffelWord(5,8)
            sage: w
            word: 0010010100101
            sage: P = WordPaths([0,1,2,3])                                              # needs sage.modules
            sage: P(w)                                                                  # needs sage.modules
            Path: 0010010100101

        Construction of a word represented by a list from a word
        represented by a str ::

            sage: w = W('ababbbabab')
            sage: type(w)
            <class 'sage.combinat.words.word.FiniteWord_str'>
            sage: z = Word(w, datatype='list')
            sage: type(z)
            <class 'sage.combinat.words.word.FiniteWord_list'>
            sage: y = Word(w, alphabet='abc', datatype='list')
            sage: type(y)
            <class 'sage.combinat.words.word.FiniteWord_list'>

        Creation of a word from a concatenation of words::

            sage: W = FiniteWords()
            sage: w = W() * W('a')
            sage: Z = Words('ab')
            sage: Z(w)
            word: a

        Creation of a word path from a FiniteWord_iter::

            sage: w = words.FibonacciWord()
            sage: f = w[:100]
            sage: P = WordPaths([0,1,2,3])                                              # needs sage.modules
            sage: p = P(f); p                                                           # needs sage.modules
            Path: 0100101001001010010100100101001001010010...
            sage: p.length()                                                            # needs sage.modules
            100

        Creation of a word path from a FiniteWord_callable::

            sage: g = W(lambda n:n%2, length = 100)
            sage: P = WordPaths([0,1,2,3])                                              # needs sage.modules
            sage: p = P(g); p                                                           # needs sage.modules
            Path: 0101010101010101010101010101010101010101...
            sage: p.length()                                                            # needs sage.modules
            100

        Creation of a word from a pickled function::

            sage: f = lambda n : n % 10
            sage: from sage.misc.fpickle import pickle_function
            sage: s = pickle_function(f)
            sage: W(s, length=10, datatype='pickled_function')
            word: 0123456789

        If the alphabet is a subset of [0, 255], then it uses char as datatype::

            sage: type(Word([0,1,1,2,0], alphabet=list(range(256))))
            <class 'sage.combinat.words.word.FiniteWord_char'>

        If the alphabet is a subset of [0, 255], then the letters must
        convert to an unsigned char. Otherwise an error is raised before
        the check is done::

            sage: type(Word([0,1,1,2,0,257], alphabet=list(range(256))))
            Traceback (most recent call last):
            ...
            OverflowError: value too large to convert to unsigned char
            sage: type(Word([0,1,1,2,0,258], alphabet=list(range(257))))
            Traceback (most recent call last):
            ...
            ValueError: 258 not in alphabet
            sage: type(Word([0,1,1,2,0,103], alphabet=list(range(100))))
            Traceback (most recent call last):
            ...
            ValueError: 103 not in alphabet
        """
        if datatype is not None:
            if datatype == 'list':
                w = self._element_classes['list'](self, data)
            elif datatype == 'char':
                w = self._element_classes['char'](self, data)
            elif datatype == 'tuple':
                w = self._element_classes['tuple'](self, data)
            elif datatype == 'str':
                w = self._element_classes['str'](self, data)
            elif datatype == 'callable':
                w = self._word_from_callable(data, length, caching)
            elif datatype == 'iter':
                w = self._word_from_iter(data, length, caching)
            elif datatype == 'pickled_function':
                from sage.misc.fpickle import unpickle_function
                data = unpickle_function(data)
                w = self._word_from_callable(data, length, caching)
            else:
                raise ValueError("unknown datatype (={})".format(datatype))

        elif 'char' in self._element_classes:
            if data is None:
                data = []
            elif callable(data):
                data = [data(i) for i in range(length)]
            elif not isinstance(data, (tuple, list)):
                data = list(data)
            w = self._element_classes['char'](self, data)

        elif isinstance(data, list):
            w = self._element_classes['list'](self, data)

        elif data is None:
            w = self._element_classes['list'](self, [])

        elif isinstance(data, str):
            w = self._element_classes['str'](self, data)

        elif isinstance(data, tuple):
            w = self._element_classes['tuple'](self, data)

        elif isinstance(data, (CombinatorialObject, ClonableElement)):
            w = self._element_classes['list'](self, list(data))

        elif callable(data):
            w = self._word_from_callable(data, length, caching)

        elif isinstance(data, Iterable):
            from sage.combinat.words.abstract_word import Word_class
            if isinstance(data, Word_class):
                w = self._word_from_word(data)
            else:
                w = self._word_from_iter(data, length, caching)

        else:
            raise ValueError("cannot guess a datatype from data (=%s); please specify one" % data)

        if check:
            self._check(w)
        return w

    def _repr_(self) -> str:
        """
        EXAMPLES::

            sage: FiniteWords() # indirect doctest
            Finite words over Set of Python objects of class 'object'
        """
        return 'Finite words over {!r}'.format(self.alphabet())

    def _an_element_(self):
        r"""
        Return an element of ``self``.

        EXAMPLES::

            sage: FiniteWords(4).an_element() # indirect doctest
            word: 212
            sage: FiniteWords([5, 1, 9]).an_element() # indirect doctest
            word: 151
            sage: FiniteWords([1]).an_element() # indirect doctest
            word: 111
            sage: FiniteWords(NN).an_element() # indirect doctest
            word: 101
        """
        try:
            some_letters = list(self.alphabet().some_elements())
        except (TypeError, ValueError, AttributeError, NotImplementedError):
            return self([])

        if len(some_letters) == 1:
            return self([some_letters[0]] * 3)

        a, b = some_letters[:2]
        return self([b, a, b])

    def iterate_by_length(self, l=1):
        r"""
        Return an iterator over all the words of ``self`` of length `l`.

        INPUT:

        - ``l`` -- integer (default: 1); the length of the desired words

        EXAMPLES::

            sage: W = FiniteWords('ab')
            sage: list(W.iterate_by_length(1))
            [word: a, word: b]
            sage: list(W.iterate_by_length(2))
            [word: aa, word: ab, word: ba, word: bb]
            sage: list(W.iterate_by_length(3))
            [word: aaa,
             word: aab,
             word: aba,
             word: abb,
             word: baa,
             word: bab,
             word: bba,
             word: bbb]
            sage: list(W.iterate_by_length('a'))
            Traceback (most recent call last):
            ...
            TypeError: the parameter l (='a') must be an integer

        TESTS::

            sage: W = FiniteWords(NN)
            sage: list(W.iterate_by_length(1))
            Traceback (most recent call last):
            ...
            NotImplementedError: cannot iterate over words for infinite alphabets
        """
        if not isinstance(l, (int, Integer)):
            raise TypeError("the parameter l (=%r) must be an integer" % l)
        cls = self._element_classes['tuple']
        if not self.alphabet().is_finite():
            raise NotImplementedError("cannot iterate over words for infinite alphabets")
        for w in itertools.product(self.alphabet(), repeat=l):
            yield cls(self, w)

    def __iter__(self):
        r"""
        Return an iterator over all the words of ``self``.

        The iterator outputs the words in shortlex order (see
        :wikipedia:`Shortlex_order`), i.e. first by increasing length and then
        lexicographically.

        EXAMPLES::

            sage: W = Words([4,5], infinite=False)
            sage: for w in W:
            ....:   if len(w)>3:
            ....:       break
            ....:   else:
            ....:       w
            word:
            word: 4
            word: 5
            word: 44
            word: 45
            ...
            word: 554
            word: 555
            sage: W = Words([5,4], infinite=False)
            sage: for w in W:
            ....:   if len(w)>3:
            ....:       break
            ....:   else:
            ....:       w
            word:
            word: 5
            word: 4
            word: 55
            word: 54
            ...
            word: 445
            word: 444
        """
        for l in itertools.count():
            yield from self.iterate_by_length(l)

    def __contains__(self, x) -> bool:
        """
        Test whether ``self`` contains ``x``.

        OUTPUT:

            This method returns ``True`` if ``x`` is a word of the appropriate
            length and the alphabets of the parents match. It returns ``False``
            otherwise.

        EXAMPLES::

            sage: W = FiniteWords('ab')
            sage: W('ab') in W
            True
            sage: W('aa') in FiniteWords('aa')
            False
            sage: FiniteWords('a')('aa') in FiniteWords('ab')
            False
            sage: 2 in FiniteWords([1,2,3])
            False
            sage: [2] in FiniteWords([1,2,3])
            False
            sage: [1, 'a'] in FiniteWords([1,2,3])
            False
        """
        from sage.combinat.words.finite_word import FiniteWord_class
        return isinstance(x, FiniteWord_class) and x.parent().alphabet() == self.alphabet()

    def random_element(self, length=None, *args, **kwds):
        r"""
        Return a random finite word on the given alphabet.

        INPUT:

        - ``length`` -- (optional) the length of the word; if not set, will use
          a uniformly random number between 0 and 10

        - all other argument are transmitted to the random generator of the
          alphabet

        EXAMPLES::

            sage: W = FiniteWords(5)
            sage: W.random_element() # random
            word: 5114325445423521544531411434451152142155...

            sage: W = FiniteWords(ZZ)
            sage: W.random_element() # random
            word: 5,-1,-1,-1,0,0,0,0,-3,-11
            sage: W.random_element(length=4, x=0, y=4) # random
            word: 1003

        TESTS::

            sage: _ = FiniteWords(GF(5)).random_element()                               # needs sage.rings.finite_rings
        """
        if length is None:
            length = ZZ.random_element(0, 10)
        return self([self.alphabet().random_element(*args, **kwds)
                     for x in range(length)])

    def iter_morphisms(self, arg=None, codomain=None, min_length=1):
        r"""
        Iterate over all morphisms with domain ``self`` and the given
        codomain.

        INPUT:

        - ``arg`` -- (default: ``None``) it can be one of the following:

          - ``None`` -- then the method iterates through all morphisms

          - tuple `(a, b)` of two integers -- it specifies the range
            ``range(a, b)`` of values to consider for the sum of the length
            of the image of each letter in the alphabet

          - list of nonnegative integers -- the length of the list must be
            equal to the size of the alphabet, and the `i`-th integer of
            ``arg`` determines the length of the word mapped to by the `i`-th
            letter of the (ordered) alphabet

        - ``codomain`` -- (default: ``None``) a combinatorial class of words;
          by default, ``codomain`` is ``self``

        - ``min_length`` -- nonnegative integer (default: 1); if ``arg`` is
          not specified, then iterate through all the morphisms where the
          length of the images of each letter in the alphabet is at least
          ``min_length``. This is ignored if ``arg`` is a list.

        OUTPUT: iterator

        EXAMPLES:

        Iterator over all non-erasing morphisms::

            sage: W = FiniteWords('ab')
            sage: it = W.iter_morphisms()
            sage: for _ in range(7): next(it)
            WordMorphism: a->a, b->a
            WordMorphism: a->a, b->b
            WordMorphism: a->b, b->a
            WordMorphism: a->b, b->b
            WordMorphism: a->aa, b->a
            WordMorphism: a->aa, b->b
            WordMorphism: a->ab, b->a

        Iterator over all morphisms including erasing morphisms::

            sage: W = FiniteWords('ab')
            sage: it = W.iter_morphisms(min_length=0)
            sage: for _ in range(7): next(it)
            WordMorphism: a->, b->
            WordMorphism: a->a, b->
            WordMorphism: a->b, b->
            WordMorphism: a->, b->a
            WordMorphism: a->, b->b
            WordMorphism: a->aa, b->
            WordMorphism: a->ab, b->

        Iterator over morphisms where the sum of the lengths of the images
        of the letters is in a specific range::

            sage: for m in W.iter_morphisms((0, 3), min_length=0): m
            WordMorphism: a->aa, b->
            WordMorphism: a->ab, b->
            WordMorphism: a->ba, b->
            WordMorphism: a->bb, b->
            WordMorphism: a->a, b->a
            WordMorphism: a->a, b->b
            WordMorphism: a->b, b->a
            WordMorphism: a->b, b->b
            WordMorphism: a->a, b->
            WordMorphism: a->b, b->
            WordMorphism: a->, b->aa
            WordMorphism: a->, b->ab
            WordMorphism: a->, b->ba
            WordMorphism: a->, b->bb
            WordMorphism: a->, b->a
            WordMorphism: a->, b->b
            WordMorphism: a->, b->

        ::

            sage: for m in W.iter_morphisms( (2, 4) ): m
            WordMorphism: a->aa, b->a
            WordMorphism: a->aa, b->b
            WordMorphism: a->ab, b->a
            WordMorphism: a->ab, b->b
            WordMorphism: a->ba, b->a
            WordMorphism: a->ba, b->b
            WordMorphism: a->bb, b->a
            WordMorphism: a->bb, b->b
            WordMorphism: a->a, b->aa
            WordMorphism: a->a, b->ab
            WordMorphism: a->a, b->ba
            WordMorphism: a->a, b->bb
            WordMorphism: a->b, b->aa
            WordMorphism: a->b, b->ab
            WordMorphism: a->b, b->ba
            WordMorphism: a->b, b->bb
            WordMorphism: a->a, b->a
            WordMorphism: a->a, b->b
            WordMorphism: a->b, b->a
            WordMorphism: a->b, b->b

        Iterator over morphisms with specific image lengths::

            sage: for m in W.iter_morphisms([0, 0]): m
            WordMorphism: a->, b->
            sage: for m in W.iter_morphisms([0, 1]): m
            WordMorphism: a->, b->a
            WordMorphism: a->, b->b
            sage: for m in W.iter_morphisms([2, 1]): m
            WordMorphism: a->aa, b->a
            WordMorphism: a->aa, b->b
            WordMorphism: a->ab, b->a
            WordMorphism: a->ab, b->b
            WordMorphism: a->ba, b->a
            WordMorphism: a->ba, b->b
            WordMorphism: a->bb, b->a
            WordMorphism: a->bb, b->b
            sage: for m in W.iter_morphisms([2, 2]): m
            WordMorphism: a->aa, b->aa
            WordMorphism: a->aa, b->ab
            WordMorphism: a->aa, b->ba
            WordMorphism: a->aa, b->bb
            WordMorphism: a->ab, b->aa
            WordMorphism: a->ab, b->ab
            WordMorphism: a->ab, b->ba
            WordMorphism: a->ab, b->bb
            WordMorphism: a->ba, b->aa
            WordMorphism: a->ba, b->ab
            WordMorphism: a->ba, b->ba
            WordMorphism: a->ba, b->bb
            WordMorphism: a->bb, b->aa
            WordMorphism: a->bb, b->ab
            WordMorphism: a->bb, b->ba
            WordMorphism: a->bb, b->bb

        The codomain may be specified as well::

            sage: Y = FiniteWords('xyz')
            sage: for m in W.iter_morphisms([0, 2], codomain=Y): m
            WordMorphism: a->, b->xx
            WordMorphism: a->, b->xy
            WordMorphism: a->, b->xz
            WordMorphism: a->, b->yx
            WordMorphism: a->, b->yy
            WordMorphism: a->, b->yz
            WordMorphism: a->, b->zx
            WordMorphism: a->, b->zy
            WordMorphism: a->, b->zz
            sage: for m in Y.iter_morphisms([0,2,1], codomain=W): m
            WordMorphism: x->, y->aa, z->a
            WordMorphism: x->, y->aa, z->b
            WordMorphism: x->, y->ab, z->a
            WordMorphism: x->, y->ab, z->b
            WordMorphism: x->, y->ba, z->a
            WordMorphism: x->, y->ba, z->b
            WordMorphism: x->, y->bb, z->a
            WordMorphism: x->, y->bb, z->b
            sage: it = W.iter_morphisms(codomain=Y)
            sage: for _ in range(10): next(it)
            WordMorphism: a->x, b->x
            WordMorphism: a->x, b->y
            WordMorphism: a->x, b->z
            WordMorphism: a->y, b->x
            WordMorphism: a->y, b->y
            WordMorphism: a->y, b->z
            WordMorphism: a->z, b->x
            WordMorphism: a->z, b->y
            WordMorphism: a->z, b->z
            WordMorphism: a->xx, b->x

        TESTS::

            sage: list(W.iter_morphisms([1,0]))
            [WordMorphism: a->a, b->, WordMorphism: a->b, b->]
            sage: list(W.iter_morphisms([0,0], codomain=Y))
            [WordMorphism: a->, b->]
            sage: list(W.iter_morphisms([0, 1, 2]))
            Traceback (most recent call last):
            ...
            TypeError: arg (=[0, 1, 2]) must be an iterable of 2 integers
            sage: list(W.iter_morphisms([0, 'a']))
            Traceback (most recent call last):
            ...
            TypeError: arg (=[0, 'a']) must be an iterable of 2 integers
            sage: list(W.iter_morphisms([0, 1], codomain='a'))
            Traceback (most recent call last):
            ...
            TypeError: codomain (=a) must be an instance of FiniteWords
        """
        n = self.alphabet().cardinality()
        min_length = max(min_length, 0)
        # create an iterable of compositions (all "compositions" if arg is
        # None, or [arg] otherwise)
        if arg is None:
            from sage.combinat.integer_lists.nn import IntegerListsNN
            compositions = IntegerListsNN(length=n, min_part=min_length)
        elif isinstance(arg, tuple):
            from sage.combinat.integer_lists import IntegerListsLex
            a, b = arg
            compositions = IntegerListsLex(min_sum=a, max_sum=b - 1,
                                           length=n, min_part=min_length)
        else:
            arg = list(arg)
            if (not len(arg) == n or not
                    all(isinstance(a, (int, Integer)) for a in arg)):
                raise TypeError(
                    "arg (=%s) must be an iterable of %s integers" % (arg, n))
            compositions = [arg]

        # set the codomain
        if codomain is None:
            codomain = self
        elif isinstance(codomain, FiniteOrInfiniteWords):
            codomain = codomain.finite_words()
        elif not isinstance(codomain, FiniteWords):
            raise TypeError("codomain (=%s) must be an instance of FiniteWords" % codomain)

        # iterate through the morphisms
        from sage.combinat.words.morphism import WordMorphism
        for composition in compositions:
            cuts = [0] + list(composition)
            for i in range(1, len(cuts)):
                cuts[i] += cuts[i - 1]
            s = cuts[-1]  # same but better than s = sum(composition)
            for big_word in codomain.iterate_by_length(s):
                d = {}
                i = 0
                for a in self.alphabet():
                    d[a] = big_word[cuts[i]:cuts[i + 1]]
                    i += 1
                yield WordMorphism(d, codomain=codomain)


class InfiniteWords(AbstractLanguage):
    def cardinality(self):
        r"""
        Return the cardinality of this set.

        EXAMPLES::

            sage: InfiniteWords('ab').cardinality()
            +Infinity
            sage: InfiniteWords('a').cardinality()
            1
            sage: InfiniteWords('').cardinality()
            0
        """
        if not self.alphabet().cardinality():
            return ZZ.zero()
        elif self.alphabet().cardinality().is_one():
            return ZZ.one()
        else:
            return Infinity

    def __hash__(self):
        r"""
        TESTS::

            sage: hash(InfiniteWords('ab')) # random
            12
        """
        return hash(self.alphabet()) ^ hash('infinite words')

    @cached_method
    def factors(self):
        r"""
        Return the set of finite words on the same alphabet.

        EXAMPLES::

            sage: InfiniteWords('ab').factors()
            Finite words over {'a', 'b'}
        """
        return FiniteWords(self.alphabet())

    def shift(self):
        r"""
        Return itself.

        EXAMPLES::

            sage: InfiniteWords('ab').shift()
            Infinite words over {'a', 'b'}
        """
        return self

    @lazy_attribute
    def _element_classes(self):
        r"""
        Return a dictionary that gives the class of the element of ``self``.

        The word may be finite, infinite or of unknown length.
        Its data may be string, list, tuple, a callable or an iterable.
        For callable and iterable, the data may be cached.

        EXAMPLES:

        Once you get the class, it can be used to create a new word::

            sage: W = InfiniteWords([0,1,2])
            sage: cls = W._element_classes['iter_with_caching']
            sage: from itertools import count
            sage: w = cls(W, (i%3 for i in count()))
            sage: type(w)
            <class 'sage.combinat.words.word.InfiniteWord_iter_with_caching'>
            sage: w
            word: 0120120120120120120120120120120120120120...
            sage: w.parent()
            Infinite words over {0, 1, 2}

        TESTS::

            sage: d = InfiniteWords()._element_classes
            sage: type(d)
            <class 'dict'>
            sage: len(d)
            4
            sage: e = InfiniteWords('abcdefg')._element_classes
            sage: d == e
            True
        """
        import sage.combinat.words.word as word
        return {
            'callable_with_caching': word.InfiniteWord_callable_with_caching,
            'callable': word.InfiniteWord_callable,
            'iter_with_caching': word.InfiniteWord_iter_with_caching,
            'iter': word.InfiniteWord_iter}

    def random_element(self, *args, **kwds):
        r"""
        Return a random infinite word.

        EXAMPLES::

            sage: W = InfiniteWords('ab')
            sage: W.random_element() # random
            word: abbbabbaabbbabbabbaabaabbabbbbbbbbaabbbb...

            sage: W = InfiniteWords(ZZ)
            sage: W.random_element(x=2,y=4) # random
            word: 3333223322232233333223323223222233233233...
        """
        rd = self.alphabet().random_element
        from itertools import count
        return self._word_from_iter(rd(*args, **kwds) for i in count())

    def _word_from_word(self, data):
        r"""
        Return a word from a word.

        The data is assumed to be ok, no check is performed.

        INPUT:

        - ``data`` -- word

        EXAMPLES::

            sage: W = InfiniteWords([0,1,2])
            sage: w = W(words.FibonacciWord())
            sage: w
            word: 0100101001001010010100100101001001010010...
            sage: w.parent() is W
            True
            sage: z = W._word_from_word(w)
            sage: w is z
            True
        """
        ####################
        # If `data` is already a word and if its parent is self, then
        # return `data` (no matter what the parameter length, datatype)
        ###########################
        if data.parent() is self or data.parent() == self:
            return data
        elif data.length() != Infinity:
            raise ValueError("cannot build an infinite word from a finite one")

        ###########################
        # Otherwise, if self is not the parent of `data`, then we try to
        # recover the data, the length and the datatype of the input `data`
        ###########################
        from sage.combinat.words.word_infinite_datatypes import (WordDatatype_callable,
                                                                 WordDatatype_iter)
        if isinstance(data, WordDatatype_callable):
            data = data._func
            return self._word_from_callable(data, caching=False)
        elif isinstance(data, WordDatatype_iter):
            data = iter(data)
            return self._word_from_iter(data, caching=False)
        else:
            raise TypeError("any instance of Word_class must be an instance of WordDatatype")

    def _word_from_callable(self, data, caching=True):
        r"""
        Return a word represented by a callable.

        The data is assumed to be ok, no check is performed.

        INPUT:

        - ``data`` -- callable

        - ``caching`` -- boolean (default: ``True``); whether to keep a cache
          of the letters computed by the callable

        EXAMPLES::

            sage: W = InfiniteWords([0,1,2])
            sage: f = lambda n : n % 3
            sage: W._word_from_callable(f)
            word: 0120120120120120120120120120120120120120...
        """
        wc = '_with_caching' if caching else ""
        return self._element_classes['callable' + wc](self, data, Infinity)

    def _word_from_iter(self, data, caching=True):
        r"""
        Return a word represented by an iterator.

        The data is assumed to be ok, no check is performed.

        INPUT:

        - ``data`` -- iterable

        - ``caching`` -- boolean (default: ``True``); whether to keep a cache
          of the letters computed by the iterator

        EXAMPLES::

            sage: W = InfiniteWords([0,1,2])
            sage: from itertools import count
            sage: W._word_from_iter((i % 3 for i in count()))
            word: 0120120120120120120120120120120120120120...
        """
        wc = '_with_caching' if caching else ""
        return self._element_classes['iter' + wc](self, data, Infinity)

    def __call__(self, data=None, datatype=None, caching=True, check=True):
        r"""
        Construct a new word object with parent ``self``.

        INPUT:

        - ``data`` -- iterator or a callable

        - ``datatype`` -- (default: ``None``) ``None``, "iter", "callable" or
          "pickled_function"; if ``None``, then the function tries to guess
          this from the data

        - ``caching`` -- boolean (default: ``True``); whether to keep a
          cache of the letters computed by an iterator or callable

        - ``check`` -- boolean (default: ``True``); whether to check if
          the 40 first letters are in the parent alphabet. This is a
          check done to test for small programming errors. Since we also
          support infinite words, we cannot really implement a more
          accurate check.

        .. NOTE::

           The check makes this method about 10 times slower (20µs instead
           of 2µs), so make sure to set it to False if you know the
           alphabet is OK. Fast creation (about 1µs) of a word can be
           done using the class directly (see :meth:`_element_classes`).

        .. WARNING::

           Be careful when defining words using callables and iterators. It
           appears that islice does not pickle correctly causing various errors
           when reloading. Also, most iterators do not support copying and
           should not support pickling by extension.

        EXAMPLES:

        Word with iterator::

            sage: from itertools import count
            sage: InfiniteWords()(count())
            word: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,...

        Word with function (a 'callable')::

            sage: f = lambda n : add(Integer(n).digits(2)) % 2
            sage: InfiniteWords()(f)
            word: 0110100110010110100101100110100110010110...

        The fourty first letters of the word are checked if they are in the
        parent alphabet::

            sage: from itertools import count
            sage: InfiniteWords("ab")(("c" if i == 0 else "a" for i in count()))
            Traceback (most recent call last):
            ...
            ValueError: c not in alphabet

        Creation of a word from a word::

            sage: w = InfiniteWords([0,1,2,3])(words.FibonacciWord())
            sage: w
            word: 0100101001001010010100100101001001010010...
            sage: w.parent()
            Infinite words over {0, 1, 2, 3}
            sage: InfiniteWords([0,1,2,3])(w) is w
            True

        Creation of a word from a pickled function::

            sage: f = lambda n : n % 10
            sage: from sage.misc.fpickle import pickle_function
            sage: s = pickle_function(f)
            sage: Word(s, datatype='pickled_function')
            word: 0123456789012345678901234567890123456789...
        """
        if datatype is not None:
            if datatype == 'callable':
                w = self._word_from_callable(data, caching)
            elif datatype == 'iter':
                w = self._word_from_iter(data, caching)
            elif datatype == 'pickled_function':
                from sage.misc.fpickle import unpickle_function
                data = unpickle_function(data)
                w = self._word_from_callable(data, caching)
            else:
                raise ValueError("unknown datatype (={})".format(datatype))

        elif callable(data):
            w = self._word_from_callable(data, caching)

        elif isinstance(data, Iterable):
            from sage.combinat.words.abstract_word import Word_class
            if isinstance(data, Word_class):
                w = self._word_from_word(data)
            else:
                w = self._word_from_iter(data, caching)

        else:
            raise ValueError("cannot guess a datatype from data (=%s); please specify one" % data)

        if check:
            self._check(w)
        return w

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: Words('ab', finite=False)  # indirect doctest
            Infinite words over {'a', 'b'}
        """
        return "Infinite words over {!r}".format(self.alphabet())

    def _an_element_(self):
        r"""
        Return an element of ``self``.

        EXAMPLES::

            sage: W = Words('ac', finite=False); W
            Infinite words over {'a', 'c'}
            sage: W.an_element()
            word: accacaaccaacaccacaacaccaaccacaaccaacacca...

            sage: W = Words(NN, finite=False); W
            Infinite words over Non negative integer semiring
            sage: W.an_element()
            word: 0110100110010110100101100110100110010110...

            sage: W = Words('z', finite=False); W
            Infinite words over {'z'}
            sage: W.an_element()
            word: zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz...
        """
        some_letters = list(self.alphabet().some_elements())
        if len(some_letters) > 1:
            from sage.combinat.words.word_generators import words
            letters = some_letters[:2]
            return self(words.ThueMorseWord(alphabet=letters))
        else:
            letter = some_letters[0]
            return self(lambda n: letter)


class FiniteOrInfiniteWords(AbstractLanguage):
    def __init__(self, alphabet):
        r"""
        INPUT:

        - ``alphabet`` -- the underlying alphabet

        TESTS::

            sage: loads(dumps(Words())) == Words()
            True
        """
        AbstractLanguage.__init__(self, alphabet)

    def __setstate__(self, state):
        r"""
        TESTS::

            sage: import os, tempfile
            sage: W = Words('ab')
            sage: with tempfile.TemporaryDirectory() as d:
            ....:     filename = os.path.join(d, 'test.sobj')
            ....:     W.save(filename)
            ....:     load(filename)
            Finite and infinite words over {'a', 'b'}
        """
        # add a default to support old pickles from #19619
        self._alphabet = state.get('_alphabet', build_alphabet())

    def cardinality(self):
        r"""
        Return the cardinality of this set of words.

        EXAMPLES::

            sage: Words('abcd').cardinality()
            +Infinity
            sage: Words('a').cardinality()
            +Infinity
            sage: Words('').cardinality()
            1
        """
        return self.finite_words().cardinality()

    @lazy_attribute
    def _element_classes(self):
        r"""
        Return the element classes corresponding to words of unknown length.

        EXAMPLES::

            sage: Words('ab')._element_classes
            {'iter': <class 'sage.combinat.words.word.Word_iter'>,
             'iter_with_caching': <class 'sage.combinat.words.word.Word_iter_with_caching'>}
        """
        import sage.combinat.words.word as word
        return {'iter_with_caching': word.Word_iter_with_caching,
                'iter': word.Word_iter}

    def __hash__(self):
        r"""
        TESTS::

            sage: hash(Words('ab')) # random
            12
        """
        return hash(self.alphabet()) ^ hash('words')

    @cached_method
    def finite_words(self):
        r"""
        Return the set of finite words.

        EXAMPLES::

            sage: Words('ab').finite_words()
            Finite words over {'a', 'b'}
        """
        return FiniteWords(self.alphabet())

    factors = finite_words

    @cached_method
    def infinite_words(self):
        r"""
        Return the set of infinite words.

        EXAMPLES::

            sage: Words('ab').infinite_words()
            Infinite words over {'a', 'b'}
        """
        return InfiniteWords(self.alphabet())

    shift = infinite_words

    def iterate_by_length(self, length):
        r"""
        Return an iterator over the words of given length.

        EXAMPLES::

            sage: [w.string_rep() for w in Words('ab').iterate_by_length(3)]
            ['aaa', 'aab', 'aba', 'abb', 'baa', 'bab', 'bba', 'bbb']
        """
        return self.finite_words().iterate_by_length(length)

    def _word_from_word(self, data):
        r"""
        TESTS::

            sage: W = Words('ab')
            sage: w = FiniteWords('abc')('abba')
            sage: W._word_from_word(w)
            word: abba
            sage: _.parent()
            Finite words over {'a', 'b'}

            sage: w = InfiniteWords('abc')(lambda i: 'a')
            sage: W._word_from_word(w)
            word: aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa...
            sage: _.parent()
            Infinite words over {'a', 'b'}
        """
        P = data.parent()
        if P is self or P is self.finite_words() or P is self.infinite_words() or \
           P == self or P == self.finite_words() or P == self.infinite_words():
            return data
        elif data.is_finite():
            return self.finite_words()._word_from_word(data)
        else:
            return self.infinite_words()._word_from_word(data)

    def _word_from_iter(self, data, caching=True):
        r"""
        TESTS::

            sage: W = Words([0,1,2])
            sage: u = Word(iter("abcabc"*100))
            sage: type(u)
            <class 'sage.combinat.words.word.Word_iter_with_caching'>
            sage: u.length() is None
            True

            sage: u = Word(iter("abcabc"))
            sage: type(u)
            <class 'sage.combinat.words.word.FiniteWord_iter_with_caching'>
            sage: u.length()
            6
        """
        wc = '_with_caching' if caching else ''
        cls = self._element_classes['iter' + wc]
        return cls(self, data, None)

    def __call__(self, data=None, length=None, datatype=None, caching=True, check=True):
        r"""
        Construct a new word object with parent ``self``.

        INPUT:

        - ``data`` -- (default: ``None``) list, string, tuple, iterator, ``None``
          (shorthand for []), or a callable defined on [0,1,...,length]

        - ``length`` -- (default: ``None``) this is dependent on the type of data.
          It is ignored for words defined by lists, strings, tuples,
          etc., because they have a naturally defined length.
          For callables, this defines the domain of definition,
          which is assumed to be [0, 1, 2, ..., length-1].
          For iterators: Infinity if you know the iterator will not
          terminate (default); ``'unknown'`` if you do not know whether the
          iterator terminates; ``'finite'`` if you know that the iterator
          terminates, but do not know the length.

        - ``datatype`` -- (default: ``None``) ``None``, "char", "list", "str",
          "tuple", "iter", "callable" or "pickled_function"; if ``None``, then
          the function tries to guess this from the data.

        - ``caching`` -- boolean (default: ``True``); whether to keep a cache
          of the letters computed by an iterator or callable

        - ``check`` -- boolean (default: ``True``); whether to check if
          the 40 first letters are in the parent alphabet. This is a
          check done to test for small programming errors. Since we also
          support infinite words, we cannot really implement a more
          accurate check.

        .. NOTE::

           The check makes this method about 10 times slower (20µs instead
           of 2µs), so make sure to set it to False if you know the
           alphabet is OK. Fast creation (about 1µs) of a word can be
           done using the class directly (see :meth:`_element_classes`).

        .. WARNING::

           Be careful when defining words using callables and iterators. It
           appears that islice does not pickle correctly causing various errors
           when reloading. Also, most iterators do not support copying and
           should not support pickling by extension.

        EXAMPLES:

        Empty word::

            sage: Words()()
            word:

        Word with string::

            sage: Words()("abbabaab")
            word: abbabaab

        Word with string constructed from other types::

            sage: Words()([0,1,1,0,1,0,0,1], datatype='str')
            word: 01101001
            sage: Words()((0,1,1,0,1,0,0,1), datatype='str')
            word: 01101001

        Word with list::

            sage: Words()([0,1,1,0,1,0,0,1])
            word: 01101001

        Word with list constructed from other types::

            sage: Words()("01101001", datatype='list')
            word: 01101001
            sage: Words()((0,1,1,0,1,0,0,1), datatype='list')
            word: 01101001

        Word with tuple::

            sage: Words()((0,1,1,0,1,0,0,1))
            word: 01101001

        Word with tuple constructed from other types::

            sage: Words()([0,1,1,0,1,0,0,1], datatype='tuple')
            word: 01101001
            sage: Words()("01101001", datatype='str')
            word: 01101001

        Word with iterator::

            sage: from itertools import count
            sage: Words()(count())
            word: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,...
            sage: Words()(iter("abbabaab")) # iterators default to infinite words
            word: abbabaab
            sage: Words()(iter("abbabaab"), length='unknown')
            word: abbabaab
            sage: Words()(iter("abbabaab"), length='finite')
            word: abbabaab

        Word with function (a 'callable')::

            sage: f = lambda n : add(Integer(n).digits(2)) % 2
            sage: Words()(f)
            word: 0110100110010110100101100110100110010110...
            sage: Words()(f, length=8)
            word: 01101001

        Word over a string with a parent::

            sage: w = Words('abc')("abbabaab"); w
            word: abbabaab
            sage: w.parent()
            Finite words over {'a', 'b', 'c'}

        The fourty first letters of the word are checked if they are in the
        parent alphabet::

            sage: Words("ab")("abca")
            Traceback (most recent call last):
            ...
            ValueError: c not in alphabet
            sage: Words("ab")("abca", check=False)
            word: abca

        The default parent is the combinatorial class of all words::

            sage: w = Words()("abbabaab"); w
            word: abbabaab
            sage: w.parent()
            Finite words over Set of Python objects of class 'object'

        Creation of a word from a word::

            sage: Words([0,1,2,3])(Words([2,3])([2,2,2,3,3,2]))
            word: 222332
            sage: _.parent()
            Finite words over {0, 1, 2, 3}

        ::

            sage: Words([3,2,1])(Words([2,3])([2,2,2,3,3,2]))
            word: 222332
            sage: _.parent()
            Finite words over {3, 2, 1}

        Construction of a word from a word when the parents are the same::

            sage: W = Words()
            sage: w = W(range(8))
            sage: z = W(w)
            sage: w is z
            True

        Construction of a word path from a finite word::

            sage: W = Words('abcd')
            sage: P = WordPaths('abcd')                                                 # needs sage.modules
            sage: w = W('aaab')
            sage: P(w)                                                                  # needs sage.modules
            Path: aaab

        Construction of a word path from a Christoffel word::

            sage: w = words.ChristoffelWord(5,8)
            sage: w
            word: 0010010100101
            sage: P = WordPaths([0,1,2,3])                                              # needs sage.modules
            sage: P(w)                                                                  # needs sage.modules
            Path: 0010010100101

        Construction of a word represented by a list from a word
        represented by a str ::

            sage: w = Word('ababbbabab')
            sage: type(w)
            <class 'sage.combinat.words.word.FiniteWord_str'>
            sage: z = Word(w, datatype='list')
            sage: type(z)
            <class 'sage.combinat.words.word.FiniteWord_list'>
            sage: y = Word(w, alphabet='abc', datatype='list')
            sage: type(y)
            <class 'sage.combinat.words.word.FiniteWord_list'>

        Creation of a word from a concatenation of words::

            sage: W = Words()
            sage: w = W() * W('a')
            sage: Z = Words('ab')
            sage: Z(w)
            word: a

        Creation of a word path from a FiniteWord_iter::

            sage: w = words.FibonacciWord()
            sage: f = w[:100]
            sage: P = WordPaths([0,1,2,3])                                              # needs sage.modules
            sage: p = P(f); p                                                           # needs sage.modules
            Path: 0100101001001010010100100101001001010010...
            sage: p.length()                                                            # needs sage.modules
            100

        Creation of a word path from a :class:`FiniteWord_callable`::

            sage: g = Word(lambda n: n%2, length=100)
            sage: P = WordPaths([0,1,2,3])                                              # needs sage.modules
            sage: p = P(g); p                                                           # needs sage.modules
            Path: 0101010101010101010101010101010101010101...
            sage: p.length()                                                            # needs sage.modules
            100

        Creation of a word from a pickled function::

            sage: f = lambda n: n % 10
            sage: from sage.misc.fpickle import pickle_function
            sage: s = pickle_function(f)
            sage: Word(s, datatype='pickled_function')
            word: 0123456789012345678901234567890123456789...

        If the alphabet is a subset of [0, 255], then it uses char as datatype::

            sage: type(Word([0,1,1,2,0], alphabet=list(range(256))))
            <class 'sage.combinat.words.word.FiniteWord_char'>

        If the alphabet is a subset of [0, 255], then the letters must
        convert to an unsigned char. Otherwise an error is raised before
        the check is done::

            sage: type(Word([0,1,1,2,0,257], alphabet=list(range(256))))
            Traceback (most recent call last):
            ...
            OverflowError: value too large to convert to unsigned char
            sage: type(Word([0,1,1,2,0,258], alphabet=list(range(257))))
            Traceback (most recent call last):
            ...
            ValueError: 258 not in alphabet
            sage: type(Word([0,1,1,2,0,103], alphabet=list(range(100))))
            Traceback (most recent call last):
            ...
            ValueError: 103 not in alphabet

        Check that the type is rightly guessed for parking functions which are
        callable::

            sage: p = ParkingFunction([2,2,1])
            sage: Word(p).parent()
            Finite words over Set of Python objects of class 'object'
        """
        # try to guess `length` from the `datatype` or `data` if not given
        if length is None or length == 'unknown':
            if data is None:
                length = 'finite'
            elif datatype in ('callable', 'pickled_function'):
                length = 'infinite'
            elif datatype in ('list', 'char', 'str', 'tuple'):
                length = 'finite'
            elif datatype is None:
                try:
                    length = len(data)
                except TypeError:
                    if callable(data):
                        length = 'infinite'

        # now build finite/infinite or unknown length words
        if length == 'finite' or length in ZZ:
            return self.finite_words()(data, datatype=datatype, length=length, caching=caching, check=check)

        elif length == 'infinite' or length == Infinity:
            return self.infinite_words()(data, datatype=datatype, check=check, caching=caching)

        elif length == 'unknown' or length is None:
            from sage.combinat.words.abstract_word import Word_class
            if isinstance(data, Word_class):
                w = self._word_from_word(data)
            elif isinstance(data, Iterable):
                w = self._word_from_iter(data, caching)
            else:
                raise ValueError("cannot guess a datatype from data (={!r}); please specify one".format(data))

            if check:
                w.parent()._check(w)
            return w

        else:
            raise ValueError("invalid argument length (={!r})".format(length))

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: Words('ab', finite=False)._repr_()
            "Infinite words over {'a', 'b'}"
        """
        return "Finite and infinite words over {!r}".format(self.alphabet())


class Words_n(Parent):
    r"""
    The set of words of fixed length on a given alphabet.
    """
    def __init__(self, words, n):
        r"""
        INPUT:

        - ``words`` -- set of finite words

        - ``n`` -- nonnegative integer

        TESTS::

            sage: Words([0,1], length=-42)
            Traceback (most recent call last):
            ...
            ValueError: n = -42 must be nonnegative
        """
        n = ZZ(n)
        if n < 0:
            raise ValueError("n = {} must be nonnegative".format(n))
        self._words = words
        self._n = n

        Parent.__init__(self, category=Sets(), facade=(words,))

    def __setstate__(self, state):
        r"""
        TESTS::

            sage: import os, tempfile
            sage: W = Words('ab', 10)
            sage: with tempfile.TemporaryDirectory() as d:
            ....:     filename = os.path.join(d, 'test.sobj')
            ....:     W.save(filename)
            ....:     load(filename)
            Words of length 10 over {'a', 'b'}
        """
        # add a default to support old pickles from #19619
        self._n = state.get('_n')
        self._words = state.get('_words', FiniteWords())

    def alphabet(self):
        r"""
        Return the underlying alphabet.

        EXAMPLES::

            sage: Words([0,1], 4).alphabet()
            {0, 1}
        """
        return self._words.alphabet()

    def __call__(self, data, *args, **kwds):
        r"""
        INPUT:

        - all arguments are sent directly to the underlying set of finite words.
          See the documentation there for the actual input.

        TESTS::

            sage: Words(5,3)([1,2,3])
            word: 123
            sage: Words(5,3)([1,2,3,1])
            Traceback (most recent call last):
            ...
            ValueError: wrong length
        """
        if 'length' in kwds:
            if kwds['length'] != self._n:
                raise ValueError("wrong length")
        else:
            kwds['length'] = self._n
        w = self._words(data, *args, **kwds)

        if kwds.get('check', True):
            if w.length() != self._n:
                raise ValueError("wrong length")
        return w

    def list(self):
        r"""
        Return a list of all the words contained in ``self``.

        EXAMPLES::

            sage: Words(0,0).list()
            [word: ]
            sage: Words(5,0).list()
            [word: ]
            sage: Words(['a','b','c'],0).list()
            [word: ]
            sage: Words(5,1).list()
            [word: 1, word: 2, word: 3, word: 4, word: 5]
            sage: Words(['a','b','c'],2).list()
            [word: aa, word: ab, word: ac, word: ba, word: bb, word: bc, word: ca, word: cb, word: cc]
        """
        return list(self)

    def _an_element_(self):
        r"""
        Return an element of ``self``.

        EXAMPLES::

            sage: W = Words(2, 3); W
            Words of length 3 over {1, 2}
            sage: W.an_element()
            word: 121

            sage: W = Words("bac", 7); W
            Words of length 7 over {'b', 'a', 'c'}
            sage: W.an_element()
            word: bacbacb

            sage: W = Words("baczxy", 5); W
            Words of length 5 over {'b', 'a', 'c', 'z', 'x', 'y'}
            sage: W.an_element()
            word: baczx
        """
        letters = list(self.alphabet().some_elements())
        r = self._n % len(letters)
        q = (self._n - r) / len(letters)
        return self(letters * int(q) + letters[:r])

    def random_element(self, *args, **kwds):
        r"""
        Return a random word in this set.

        EXAMPLES::

            sage: W = Words('ab', 4)
            sage: W.random_element()  # random
            word: bbab
            sage: W.random_element() in W
            True

            sage: W = Words(ZZ, 5)
            sage: W.random_element()  # random
            word: 1,2,2,-1,12
            sage: W.random_element() in W
            True

        TESTS::

            sage: _ = Words(GF(5),4).random_element()                                   # needs sage.rings.finite_rings

        Check that :issue:`18283` is fixed::

            sage: w = Words('abc', 5).random_element()
            sage: w.length()
            5
        """
        return self._words.random_element(length=self._n, *args, **kwds)

    def _repr_(self):
        """
        EXAMPLES::

            sage: Words(3,5) # indirect doctest
            Words of length 5 over {1, 2, 3}
        """
        from sage.combinat.words.word_options import word_options
        if word_options['old_repr']:
            return "Words over {} of length {}".format(self.alphabet(), self._n)
        return "Words of length {} over {}".format(self._n, self.alphabet())

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: W = Words(3,5)
            sage: W.an_element() in W
            True

            sage: 2 in Words(length=3)
            False
            sage: [1,'a',3] in Words(length=3)
            False
            sage: [1,2] in Words(length=3)
            False
            sage: "abc" in Words(length=3)
            False
            sage: Words("abc")("ababc") in Words(length=3)
            False
            sage: Words([0,1])([1,0,1]) in Words([0,1], length=3)
            True
        """
        return x in self._words and x.length() == self._n

    def cardinality(self):
        r"""
        Return the number of words of length `n` from alphabet.

        EXAMPLES::

            sage: Words(['a','b','c'], 4).cardinality()
            81
            sage: Words(3, 4).cardinality()
            81
            sage: Words(0,0).cardinality()
            1
            sage: Words(5,0).cardinality()
            1
            sage: Words(['a','b','c'],0).cardinality()
            1
            sage: Words(0,1).cardinality()
            0
            sage: Words(5,1).cardinality()
            5
            sage: Words(['a','b','c'],1).cardinality()
            3
            sage: Words(7,13).cardinality()
            96889010407
            sage: Words(['a','b','c','d','e','f','g'],13).cardinality()
            96889010407
        """
        return self.alphabet().cardinality() ** self._n

    __len__ = cardinality

    def __iter__(self):
        r"""
        TESTS::

            sage: [w for w in Words(['a', 'b'], 2)]
            [word: aa, word: ab, word: ba, word: bb]
            sage: [w for w in Words(['b', 'a'], 2)]
            [word: bb, word: ba, word: ab, word: aa]
            sage: [w for w in Words(['a', 'b'], 0)]
            [word: ]
            sage: [w for w in Words([], 3)]
            []
        """
        return self._words.iterate_by_length(self._n)

    def iterate_by_length(self, length):
        r"""
        All words in this class are of the same length, so use iterator
        instead.

        TESTS::

            sage: W = Words(['a', 'b'], 2)
            sage: list(W.iterate_by_length(2))
            [word: aa, word: ab, word: ba, word: bb]
            sage: list(W.iterate_by_length(1))
            []
        """
        if length == self._n:
            return iter(self)
        else:
            return iter([])


###############
# old pickles #
###############


register_unpickle_override("sage.combinat.words.words", "Words_over_OrderedAlphabet", FiniteOrInfiniteWords)
register_unpickle_override("sage.combinat.words.words", "Words_over_Alphabet", FiniteOrInfiniteWords)
register_unpickle_override("sage.combinat.words.words", "FiniteWords_length_k_over_OrderedAlphabet", Words_n)
register_unpickle_override("sage.combinat.words.words", "FiniteWords_over_OrderedAlphabet", FiniteWords)
register_unpickle_override("sage.combinat.words.words", "InfiniteWords_over_OrderedAlphabet", InfiniteWords)
