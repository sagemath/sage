# coding=utf-8
r"""
Combinatorial classes of words.

To define a new class of words, please refer to the documentation file:
sage/combinat/words/notes/word_inheritance_howto.txt

AUTHORS:

    - Franco Saliola (2008-12-17): merged into sage
    - Sebastien Labbe (2008-12-17): merged into sage
    - Arnaud Bergeron (2008-12-17): merged into sage
    - Sebastien Labbe (2009-07-21): Improved morphism iterator (#6571).

EXAMPLES::

    sage: Words()
    Words
    sage: Words(4)
    Words over {1, 2, 3, 4}
    sage: Words('ab')
    Words over {'a', 'b'}
    sage: Words('natural numbers')
    Words over Non negative integers
"""
#*****************************************************************************
#       Copyright (C) 2008 Arnaud Bergeron <abergeron@gmail.com>,
#                          Sébastien Labbé <slabqc@gmail.com>,
#                          Franco Saliola <saliola@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.combinat.combinat import InfiniteAbstractCombinatorialClass
from sage.combinat.combinat import CombinatorialObject
from sage.combinat.words.alphabet import build_alphabet
from sage.misc.lazy_attribute import lazy_attribute
from sage.plot.misc import rename_keyword
from sage.misc.mrange import xmrange
from sage.rings.all import Infinity
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
import itertools
from sage.structure.parent import Set_PythonType

def Words(alphabet=None, length=None, finite=True, infinite=True):
    """
    Returns the combinatorial class of words of length k over an alphabet.

    EXAMPLES::

        sage: Words()
        Words
        sage: Words(length=7)
        Words of length 7
        sage: Words(5)
        Words over {1, 2, 3, 4, 5}
        sage: Words(5, 3)
        Words of length 3 over {1, 2, 3, 4, 5}
        sage: Words(5, infinite=False)
        Finite Words over {1, 2, 3, 4, 5}
        sage: Words(5, finite=False)
        Infinite Words over {1, 2, 3, 4, 5}
        sage: Words('ab')
        Words over {'a', 'b'}
        sage: Words('ab', 2)
        Words of length 2 over {'a', 'b'}
        sage: Words('ab', infinite=False)
        Finite Words over {'a', 'b'}
        sage: Words('ab', finite=False)
        Infinite Words over {'a', 'b'}
        sage: Words('positive integers', finite=False)
        Infinite Words over Positive integers
        sage: Words('natural numbers')
        Words over Non negative integers
    """
    if isinstance(alphabet, Words_all):
        return alphabet
    if alphabet is None:
        if length is None:
            if finite and infinite:
                return Words_all()
            elif finite:
                raise NotImplementedError
            else:
                raise NotImplementedError
        elif isinstance(length, (int,Integer)) and finite:
            return Words_n(length)
    else:
        if isinstance(alphabet, (int,Integer)):
            from sage.sets.integer_range import IntegerRange
            alphabet = IntegerRange(1,alphabet+1)
        elif alphabet == "integers" \
                or alphabet == "positive integers" \
                or alphabet == "natural numbers":
            alphabet = build_alphabet(name=alphabet)
        else:
            alphabet = build_alphabet(data=alphabet)
        if length is None:
            if finite and infinite:
                return Words_over_OrderedAlphabet(alphabet)
            elif finite:
                return FiniteWords_over_OrderedAlphabet(alphabet)
            else:
                return InfiniteWords_over_OrderedAlphabet(alphabet)
        elif isinstance(length, (int,Integer)):
                return FiniteWords_length_k_over_OrderedAlphabet(alphabet, length)
    raise ValueError("do not know how to make a combinatorial class of words from your input")

from sage.structure.unique_representation import UniqueRepresentation
class Words_all(InfiniteAbstractCombinatorialClass):
    r"""
    TESTS::

        sage: from sage.combinat.words.words import Words_all
        sage: list(Words_all())
        Traceback (most recent call last):
        ...
        NotImplementedError
        sage: Words_all().list()
        Traceback (most recent call last):
        ...
        NotImplementedError: infinite list
        sage: Words_all().cardinality()
        +Infinity

        sage: isinstance(Words('ab'), Words_all)
        True
        sage: isinstance(33, Words_all)
        False

    We would like the instance of this class to be unique::

        sage: Words() is Words()   # todo: not implemented
        True

    .. WARNING::

       The design of these classes is not particularly robust so extra care must
       be taken when extending this class in order to prevent unintended
       side-effects. This is particularly evident in the equality test
       :meth:`__eq__` for words.
    """
    @lazy_attribute
    def _element_classes(self):
        r"""
        Returns a dictionary that gives the class of the element of self.

        The word may be finite, infinite or of unknown length.
        Its data may be str, list, tuple, a callable or an iterable.
        For callable and iterable, the data may be cached.

        EXAMPLES:

        Once you get the class, it can be used to create a new word::

            sage: W = Words([0,1,2])
            sage: L = [0,1,0] * 100
            sage: cls = W._element_classes['FiniteWord_list']
            sage: w = cls(W, L)
            sage: type(w)
            <class 'sage.combinat.words.word.FiniteWord_list'>
            sage: w
            word: 0100100100100100100100100100100100100100...
            sage: w.parent()
            Words over {0, 1, 2}

        TESTS::

            sage: d = Words()._element_classes
            sage: type(d)
            <type 'dict'>
            sage: len(d)
            13
            sage: e = Words('abcdefg')._element_classes
            sage: d == e
            True
        """
        import sage.combinat.words.word as word
        classes = {
            'FiniteWord_list': word.FiniteWord_list,
            'FiniteWord_str': word.FiniteWord_str,
            'FiniteWord_tuple': word.FiniteWord_tuple,
            'FiniteWord_callable_with_caching': word.FiniteWord_callable_with_caching,
            'FiniteWord_callable': word.FiniteWord_callable,
            'FiniteWord_iter_with_caching': word.FiniteWord_iter_with_caching,
            'FiniteWord_iter': word.FiniteWord_iter,
            'InfiniteWord_callable_with_caching': word.InfiniteWord_callable_with_caching,
            'InfiniteWord_callable': word.InfiniteWord_callable,
            'InfiniteWord_iter_with_caching': word.InfiniteWord_iter_with_caching,
            'InfiniteWord_iter': word.InfiniteWord_iter,
            'Word_iter_with_caching': word.Word_iter_with_caching,
            'Word_iter': word.Word_iter
            }

        # test whether or not we can use the class Finiteword_char
        if (self.alphabet().cardinality() <= 256 and
                all(isinstance(i, (int,Integer)) and 0 <= i < 256 for i in self.alphabet())):
            L = self.alphabet().list()
            if (all(L[i] < L[i+1] for i in range(len(L)-1)) and
                    all(self.cmp_letters(L[i],L[i+1]) == -1 for i in range(len(L)-1))):
                classes['FiniteWord_char'] = word.FiniteWord_char

        return classes

    def _an_element_(self):
        r"""
        Return an element of self.

        EXAMPLES::

            sage: W = Words(4)
            sage: W.an_element()
            word:

        TESTS:

        Check that :trac:`16125` is fixed::

            sage: W = Words(4)
            sage: w = W([3,1,2,1])
            sage: w in ZZ
            False
        """
        return self([])

    def __call__(self, data=None, length=None, datatype=None, caching=True, check=True):
        r"""
        Construct a new word object with parent self.

        INPUT:

        -  ``data`` - (default: None) list, string, tuple, iterator, None
           (shorthand for []), or a callable defined on [0,1,...,length].

        -  ``length`` - (default: None) This is dependent on the type of data.
           It is ignored for words defined by lists, strings, tuples,
           etc., because they have a naturally defined length.
           For callables, this defines the domain of definition,
           which is assumed to be [0, 1, 2, ..., length-1].
           For iterators: Infinity if you know the iterator will not
           terminate (default); "unknown" if you do not know whether the
           iterator terminates; "finite" if you know that the iterator
           terminates, but do know know the length.

        -  ``datatype`` - (default: None) None, "char", "list", "str",
           "tuple", "iter", "callable" or "pickled_function". If None, then
           the function tries to guess this from the data.

        -  ``caching`` - (default: True) True or False. Whether to keep a cache
           of the letters computed by an iterator or callable.

        -  ``check`` - (default: True) True or False. Whether to check if
           the 40 first letters are in the parent alphabet. This is a
           check done to test for small programming errors. Since we also
           support infinite words, we cannot really implement a more
           accurate check.

        .. NOTE::

           The check makes this method about 10 times slower (20µs instead
           of 2µs), so make sure to set it to False if you know the
           alphabet is OK. Fast creation (about 1µs) of a word can be
           done using the class directly (see :meth:``_element_classes``).

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

            sage: Words()([0,1,1,0,1,0,0,1], datatype="str")
            word: 01101001
            sage: Words()((0,1,1,0,1,0,0,1), datatype="str")
            word: 01101001

        Word with list::

            sage: Words()([0,1,1,0,1,0,0,1])
            word: 01101001

        Word with list constructed from other types::

            sage: Words()("01101001", datatype="list")
            word: 01101001
            sage: Words()((0,1,1,0,1,0,0,1), datatype="list")
            word: 01101001

        Word with tuple::

            sage: Words()((0,1,1,0,1,0,0,1))
            word: 01101001

        Word with tuple constructed from other types::

            sage: Words()([0,1,1,0,1,0,0,1], datatype="tuple")
            word: 01101001
            sage: Words()("01101001", datatype="str")
            word: 01101001

        Word with iterator::

            sage: from itertools import count
            sage: Words()(count())
            word: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,...
            sage: Words()(iter("abbabaab")) # iterators default to infinite words
            word: abbabaab
            sage: Words()(iter("abbabaab"), length="unknown")
            word: abbabaab
            sage: Words()(iter("abbabaab"), length="finite")
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
            Words over {'a', 'b', 'c'}

        The fourty first letters of the word are checked if they are in the
        parent alphbet::

            sage: Words("ab")("abca")
            Traceback (most recent call last):
            ...
            ValueError: c not in alphabet!
            sage: Words("ab")("abca", check=False)
            word: abca

        The default parent is the combinatorial class of all words::

            sage: w = Words()("abbabaab"); w
            word: abbabaab
            sage: w.parent()
            Words

        Creation of a word from a word::

            sage: Words([0,1,2,3])(Words([2,3])([2,2,2,3,3,2]))
            word: 222332
            sage: _.parent()
            Words over {0, 1, 2, 3}

        ::

            sage: Words([3,2,1])(Words([2,3])([2,2,2,3,3,2]))
            word: 222332
            sage: _.parent()
            Words over {3, 2, 1}

        Construction of a word from a word when the parents are the same::

            sage: W = Words()
            sage: w = W(range(8))
            sage: z = W(w)
            sage: w is z
            True

        Construction of a word path from a finite word::

            sage: W = Words('abcd')
            sage: P = WordPaths('abcd')
            sage: w = W('aaab')
            sage: P(w)
            Path: aaab

        Construction of a word path from a Christoffel word::

            sage: w = words.ChristoffelWord(5,8)
            sage: w
            word: 0010010100101
            sage: P = WordPaths([0,1,2,3])
            sage: P(w)
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
            sage: P = WordPaths([0,1,2,3])
            sage: p = P(f); p
            Path: 0100101001001010010100100101001001010010...
            sage: p.length()
            100

        Creation of a word path from a FiniteWord_callable::

            sage: g = Word(lambda n:n%2, length = 100)
            sage: P = WordPaths([0,1,2,3])
            sage: p = P(g); p
            Path: 0101010101010101010101010101010101010101...
            sage: p.length()
            100

        Creation of a word from a pickled function::

            sage: f = lambda n : n % 10
            sage: from sage.misc.fpickle import pickle_function
            sage: s = pickle_function(f)
            sage: Word(s, datatype='pickled_function')
            word: 0123456789012345678901234567890123456789...

        If the alphabet is a subset of [0, 255], then it uses char as datatype::

            sage: type(Word([0,1,1,2,0], alphabet=range(256)))
            <class 'sage.combinat.words.word.FiniteWord_char'>

        If the alphabet is a subset of [0, 255], then the letters must
        convert to an unsigned char. Otherwise an error is raised before
        the check is done::

            sage: type(Word([0,1,1,2,0,257], alphabet=range(256)))
            Traceback (most recent call last):
            ...
            OverflowError: value too large to convert to unsigned char
            sage: type(Word([0,1,1,2,0,258], alphabet=range(257)))
            Traceback (most recent call last):
            ...
            ValueError: 258 not in alphabet!
            sage: type(Word([0,1,1,2,0,103], alphabet=range(100)))
            Traceback (most recent call last):
            ...
            ValueError: 103 not in alphabet!

        """
        if datatype is not None:
            if datatype == 'list':
                w = self._element_classes['FiniteWord_list'](self, data)
            elif datatype == 'char':
                w = self._element_classes['FiniteWord_char'](self, data)
            elif datatype == 'tuple':
                w = self._element_classes['FiniteWord_tuple'](self, data)
            elif datatype == 'str':
                w = self._element_classes['FiniteWord_str'](self, data)
            elif datatype == 'callable':
                w = self._word_from_callable(data, length, caching)
            elif datatype == 'iter':
                w = self._word_from_iter(data, length, caching)
            elif datatype == 'pickled_function':
                from sage.misc.fpickle import unpickle_function
                data = unpickle_function(data)
                w = self._word_from_callable(data, length, caching)
            else:
                raise ValueError("Unknown datatype (={})".format(datatype))

        elif isinstance(data, list):
            if 'FiniteWord_char' in self._element_classes:
                w = self._element_classes['FiniteWord_char'](self, data)
            else:
                w = self._element_classes['FiniteWord_list'](self, data)

        elif data is None:
            if 'FiniteWord_char' in self._element_classes:
                w = self._element_classes['FiniteWord_char'](self, [])
            else:
                w = self._element_classes['FiniteWord_list'](self, [])

        elif isinstance(data, str):
            w = self._element_classes['FiniteWord_str'](self, data)

        elif isinstance(data, tuple):
            if 'FiniteWord_char' in self._element_classes:
                w = self._element_classes['FiniteWord_char'](self, data)
            else:
                w = self._element_classes['FiniteWord_tuple'](self, data)

        elif isinstance(data, CombinatorialObject):
            w = self._element_classes['FiniteWord_list'](self, list(data))

        elif callable(data):
            w = self._word_from_callable(data, length, caching)

        elif hasattr(data, "__iter__"):
            from sage.combinat.words.abstract_word import Word_class
            if isinstance(data, Word_class):
                w = self._word_from_word(data)
            else:
                w = self._word_from_iter(data, length, caching)

        else:
            raise ValueError("Cannot guess a datatype from data (=%s); please specify one" % data)

        if check:
            self._check(w)
        return w

    def _check(self, w, length=40):
        r"""
        Check that the first length elements are actually in the alphabet.

        INPUT:

        - ``w`` -- word

        - ``length`` -- integer (default: ``40``)

        EXAMPLES::

            sage: from sage.combinat.words.words import Words_over_Alphabet
            sage: W = Words_over_Alphabet(['a','b','c'])
            sage: W._check('abcabc') is None
            True
            sage: W._check('abcabcd')
            Traceback (most recent call last):
            ...
            ValueError: d not in alphabet!
            sage: W._check('abcabc'*10+'z') is None
            True
            sage: W._check('abcabc'*10+'z', length=80)
            Traceback (most recent call last):
            ...
            ValueError: z not in alphabet!
        """
        for a in itertools.islice(w, length):
            if a not in self._alphabet:
                raise ValueError("%s not in alphabet!" % a)

    def _word_from_word(self, data):
        r"""
        Return a word from a word.

        The data is assumed to be ok, no check is performed.

        INPUT:

        -  ``data`` - word

        EXAMPLES::

            sage: W = Words([0,1,2])
            sage: w = W([0,1,2,0,1,2])
            sage: z = W._word_from_word(w)
            sage: z
            word: 012012
            sage: w is z
            True
        """
        ####################
        # If `data` is already a word and if its parent is self, then
        # return `data` (no matter what the parameter length, datatype and
        # length are).
        ###########################
        if data.parent() is self:
            return data
        ###########################
        # Otherwise, if self is not the parent of `data`, then we try to
        # recover the data, the length and the datatype of the input `data`
        ###########################
        from sage.combinat.words.word_char import WordDatatype_char
        from sage.combinat.words.word_infinite_datatypes import (WordDatatype_callable,
                                                                 WordDatatype_iter)
        from sage.combinat.words.word_datatypes import (WordDatatype_str,
                          WordDatatype_list, WordDatatype_tuple)
        if isinstance(data, WordDatatype_callable):
            from sage.combinat.words.finite_word import CallableFromListOfWords
            if isinstance(data._func, CallableFromListOfWords):
                # The following line is important because, in this case,
                # data._func is also a tuple (indeed
                # CallableFromListOfWords inherits from tuple)
                datatype = "callable"
            length = data.length()
            data = data._func
            return self._word_from_callable(data, length, caching=False)
        elif isinstance(data, WordDatatype_iter):
            length = data.length()
            data = iter(data)
            return self._word_from_iter(data, length, caching=False)
        elif isinstance(data, WordDatatype_char):
            data = list(data)
            if 'FiniteWord_char' in self._element_classes:
                return self._element_classes['FiniteWord_char'](self, data)
            else:
                return self._element_classes['FiniteWord_list'](self, data)
        elif isinstance(data, WordDatatype_str):
            return self._element_classes['FiniteWord_str'](self, data._data)
        elif isinstance(data, WordDatatype_tuple):
            return self._element_classes['FiniteWord_tuple'](self, data._data)
        elif isinstance(data, WordDatatype_list):
            return self._element_classes['FiniteWord_list'](self, data._data)
        else:
            raise TypeError("Any instance of Word_class must be an instance of WordDatatype.")

    def _word_from_callable(self, data, length, caching=True):
        r"""
        Return a word represented by a callable.

        The data is assumed to be ok, no check is performed.

        INPUT:

        -  ``data`` - callable
        -  ``length`` - integer or ``None`` or "infinite" or ``Infinity``
        -  ``caching`` - (default: True) True or False. Whether to keep a cache
           of the letters computed by the callable.

        EXAMPLES::

            sage: W = Words([0,1,2])
            sage: f = lambda n : n % 3
            sage: W._word_from_callable(f, 100)
            word: 0120120120120120120120120120120120120120...
        """
        wc = '_with_caching' if caching else ""
        if length in (None, Infinity, 'infinite'):
            return self._element_classes['InfiniteWord_callable'+wc](self, data, length)
        elif length in ZZ and length >= 0:
            return self._element_classes['FiniteWord_callable'+wc](self, data, length)
        else:
            raise ValueError("not a correct value for length (%s)" % length)

    def _word_from_iter(self, data, length, caching=True):
        r"""
        Return a word represented by an iterator.

        The data is assumed to be ok, no check is performed.

        INPUT:

        -  ``data`` - iterable

        -  ``length`` - integer, ``"unknown"``, ``"infinite"``,
           ``Infinity``, ``None`` or ``"finite"``

        -  ``caching`` - (default: True) True or False. Whether to keep a cache
           of the letters computed by the iterator.

        EXAMPLES::

            sage: W = Words([0,1,2])
            sage: W._word_from_iter(iter([1]*10), 10)
            word: 1111111111
        """
        wc = '_with_caching' if caching else ""
        if length is None:
            length = getattr(self,'_length',None)
        if length in (None, Infinity, 'infinite'):
            return self._element_classes['InfiniteWord_iter'+wc](self, data, length)
        elif (length == 'finite') or (length in ZZ and length >= 0):
            return self._element_classes['FiniteWord_iter'+wc](self, data, length)
        elif length == 'unknown':
            return self._element_classes['Word_iter'+wc](self, data, length)
        else:
            raise ValueError("not a correct value for length (%s)" % length)

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.combinat.words.words import Words_all
            sage: Words_all()._repr_()
            'Words'
        """
        return 'Words'

    def __contains__(self, x):
        """
        Returns True if x is contained in self.

        EXAMPLES::

            sage: from sage.combinat.words.words import Words_all
            sage: 2 in Words_all()
            False
            sage: [1,2] in Words_all()
            False
            sage: Words('ab')('abba') in Words_all()
            True
        """
        from sage.combinat.words.abstract_word import Word_class
        return isinstance(x, Word_class)

    def __eq__(self, other):
        r"""
        Returns True if self is equal to other and False otherwise.

        EXAMPLES::

            sage: Words('ab') == Words()
            False
            sage: Words() == Words('ab')
            False
            sage: Words('ab') == Words('ab')
            True
            sage: Words('ab') == Words('ba')
            False
            sage: Words('ab') == Words('abc')
            False
            sage: Words('abc') == Words('ab')
            False

        ::

            sage: WordPaths('abcd') == Words('abcd')
            True
            sage: Words('abcd') == WordPaths('abcd')
            True
            sage: Words('bacd') == WordPaths('abcd')
            False
            sage: WordPaths('bacd') == WordPaths('abcd')
            False

        TESTS:

        :trac:`15480`::

            sage: Words(3, 10) == Words(3,900)
            False
            sage: Words(2, finite=False) == Words(2)
            False
            sage: Words(2) == Words(2,30)
            False
            sage: Words(10,0) == Words(20,0)
            True
            sage: WordPaths('abcd') == Words("abcd",3)
            False
            sage: Words(3) == Words(3,finite=False)
            False
            sage: Words(3) == Words(3,infinite=False)
            False
        """

        from paths import WordPaths_all
        # Specific case of Words_over_Alphabet and WordPath. See #15480
        # i.e. when self,other in Words_over_Alphabet, WordPath and one of them at least is a wordpath
        if ((isinstance(self,WordPaths_all) and isinstance(other,WordPaths_all)) or
            (type(self) is Words_over_OrderedAlphabet and isinstance(other,WordPaths_all)) or
            (type(other) is Words_over_OrderedAlphabet and isinstance(self,WordPaths_all))):
            return self.alphabet() == other.alphabet()

        if not (type(self) is type(other)):
            return False

        cardinality = self.cardinality()

        if cardinality != other.cardinality():
            return False
        if cardinality == 1:
            return True
        if self.alphabet() != other.alphabet():
            return False

        # This method's code cannot be trusted. It's the only way I see to fix
        # the wrong results reported in #15480. But really, this kind of
        # code should not be trusted. It is likely to return wrong results if
        # whenever new classes extending Words_all are added.
        return True

    def __ne__(self, other):
        r"""
        Returns True if self is not equal to other and False otherwise.

        TESTS::

            sage: Words('ab') != Words('ab')
            False
            sage: Words('ab') != Words('abc')
            True
            sage: Words('abc') != Words('ab')
            True

        ::

            sage: WordPaths('abcd') != Words('abcd')
            False
            sage: Words('abcd') != WordPaths('abcd')
            False

        ::

            Words('ab') != 2
            True
        """
        if isinstance(other, Words_all):
            return not self.__eq__(other)
        else:
            return NotImplemented

    _alphabet = Set_PythonType(object)

    def alphabet(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.words.words import Words_over_Alphabet
            sage: W = Words_over_Alphabet([1,2,3])
            sage: W.alphabet()
            [1, 2, 3]
            sage: from sage.combinat.words.alphabet import build_alphabet
            sage: W = Words_over_Alphabet(build_alphabet('ab'))
            sage: W.alphabet()
            {'a', 'b'}

            sage: w = Word('abaccefa')
            sage: w.parent().alphabet()
            Set of Python objects of type 'object'
            sage: Words('456').alphabet()
            {'4', '5', '6'}
        """
        return self._alphabet

    def size_of_alphabet(self):
        r"""
        Returns the size of the alphabet.

        EXAMPLES::

            sage: Words().size_of_alphabet()
            +Infinity
            sage: Word('abaccefa').parent().size_of_alphabet()
            +Infinity
        """
        return Infinity

    cmp_letters = cmp

    def has_letter(self, letter):
        r"""
        Returns True if the alphabet of self contains the given letter.

        INPUT:

        -  ``letter`` - a letter

        EXAMPLES::

            sage: W = Words()
            sage: W.has_letter('a')
            True
            sage: W.has_letter(1)
            True
            sage: W.has_letter({})
            True
            sage: W.has_letter([])
            True
            sage: W.has_letter(range(5))
            True
            sage: W.has_letter(Permutation([]))
            True

            sage: from sage.combinat.words.words import Words_over_Alphabet
            sage: W = Words_over_Alphabet(['a','b','c'])
            sage: W.has_letter('a')
            True
            sage: W.has_letter('d')
            False
            sage: W.has_letter(8)
            False
        """
        return letter in self._alphabet

class Words_over_Alphabet(Words_all):
    def __init__(self, alphabet):
        """
        Words over Alphabet.

        INPUT:

        -  ``alphabet`` - assumed to be an instance of Alphabet, but no
           type checking is done here.

        EXAMPLES::

            sage: from sage.combinat.words.words import Words_over_Alphabet
            sage: W = Words_over_Alphabet([1,2,3])
            sage: W == loads(dumps(W))
            True

        The input alphabet must be an instance of Alphabet::

            sage: W = Words_over_Alphabet(Alphabet([1,2,3]))
            sage: W([1,2,2,3])
            word: 1223
        """
        self._alphabet = alphabet

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.combinat.words.words import Words_over_Alphabet
            sage: Words_over_Alphabet([1,2,3])._repr_()
            'Words over [1, 2, 3]'
        """
        return "Words over %s"%self._alphabet

    def __contains__(self, x):
        """
        Tests whether self contains x.

        OUTPUT:
            This method returns True if x is a word of the appropriate
            length and the alphabets of the parents match. Returns False
            otherwise.

        EXAMPLES::

            sage: from sage.combinat.words.words import Words_over_Alphabet
            sage: from sage.combinat.words.alphabet import build_alphabet
            sage: A = build_alphabet('ab')
            sage: Words(A)('abba') in Words_over_Alphabet(A)
            True
            sage: Words(A)('aa') in Words_over_Alphabet(A)
            True
            sage: Words('a')('aa') in Words_over_Alphabet(A)
            False
            sage: 2 in Words_over_Alphabet([1,2,3])
            False
            sage: [2] in Words_over_Alphabet([1,2,3])
            False
            sage: [1, 'a'] in Words_over_Alphabet([1,2,3])
            False
        """
        from sage.combinat.words.abstract_word import Word_class
        return isinstance(x, Word_class) and x.parent().alphabet() == self.alphabet()

    def size_of_alphabet(self):
        r"""
        Returns the size of the alphabet.

        EXAMPLES::

            sage: Words('abcdef').size_of_alphabet()
            6
            sage: Words('').size_of_alphabet()
            0
            sage: Words('456').size_of_alphabet()
            3
        """
        return self.alphabet().cardinality()

    def identity_morphism(self):
        r"""
        Returns the identity morphism from self to itself.

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
            Words over Non negative integers
            sage: W.identity_morphism()
            Traceback (most recent call last):
            ...
            NotImplementedError: size of alphabet must be finite
        """
        if self.size_of_alphabet() not in ZZ:
            raise NotImplementedError('size of alphabet must be finite')
        from sage.combinat.words.morphism import WordMorphism
        return WordMorphism(dict((a,a) for a in self.alphabet()))

class Words_n(Words_all):
    def __init__(self, n):
        """
        EXAMPLES::

            sage: from sage.combinat.words.words import Words_n
            sage: w = Words_n(3)
            sage: w == loads(dumps(w))
            True
        """
        self._n = n

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.combinat.words.words import Words_n
            sage: Words_n(3).__repr__()
            'Words of length 3'
        """
        return "Words of length %s"%self._n

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: from sage.combinat.words.words import Words_n
            sage: 2 in Words_n(3)
            False
            sage: [1,'a',3] in Words_n(3)
            False
            sage: [1,2] in Words_n(3)
            False
            sage: "abc" in Words_n(3)
            False
            sage: Words("abc")("ababc") in Words_n(3)
            False
            sage: Words([0,1])([1,0,1]) in Words_n(3)
            True
        """
        from sage.combinat.words.finite_word import FiniteWord_class
        return isinstance(x, FiniteWord_class) and x.length() == self._n

class Words_over_OrderedAlphabet(Words_over_Alphabet):
    def __init__(self, alphabet):
        r"""
        EXAMPLES::

            sage: from sage.combinat.words.words import Words_over_OrderedAlphabet
            sage: from sage.combinat.words.alphabet import build_alphabet
            sage: A = build_alphabet("abc")
            sage: W = Words_over_OrderedAlphabet(A)
            sage: W == loads(dumps(W))
            True

        TESTS:

        Impossible to iterate over the uncountable set of all words on a given
        alphabet::

            sage: Words(4)[1]
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        super(Words_over_OrderedAlphabet, self).__init__(alphabet)

    def iterate_by_length(self, l=1):
        r"""
        Returns an iterator over all the words of self of length l.

        INPUT:

        - ``l`` - integer (default: 1), the length of the desired words

        EXAMPLES::

            sage: W = Words('ab')
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
        """
        if not isinstance(l, (int,Integer)):
            raise TypeError("the parameter l (=%r) must be an integer"%l)
        A = self.alphabet()
        cls = self._element_classes['FiniteWord_tuple']
        for w in itertools.product(self.alphabet(), repeat=l):
            yield cls(self, w)

    def random_element(self):
        r"""
        Returns a random (infinite) word on the given alphabet.

        The word returned has infinite length, and is built by randomly picking
        a letter from the alphabet, one after the other.

        EXAMPLE::

            sage: W = Words(5)
            sage: W.random_element() # random
            word: 5114325445423521544531411434451152142155...

            sage: W = Words(ZZ)
            sage: W.random_element()
            Traceback (most recent call last):
            ...
            ValueError: How can I pick a random word with an infinite aphabet?

        TESTS::

            sage: _ = Words(GF(5)).random_element()
        """
        if self.alphabet().cardinality() == Infinity:
            raise ValueError("How can I pick a random word with an infinite aphabet?")
        else:
            return self(self.alphabet().random_element() for x in itertools.count())

    @rename_keyword(deprecation=10134, l='arg')
    def iter_morphisms(self, arg=None, codomain=None, min_length=1):
        r"""
        Iterate over all morphisms with domain ``self`` and the given
        codomain.

        INPUT:

        - ``arg`` - (optional, default: None) It can be one of the following :

          - ``None`` - then the method iterates through all morphisms.

          - tuple `(a, b)` of two integers  - It specifies the range
            ``range(a, b)`` of values to consider for the sum of the length
            of the image of each letter in the alphabet.

          - list of nonnegative integers - The length of the list must be
            equal to the size of the alphabet, and the i-th integer of
            ``arg`` determines the length of the word mapped to by the i-th
            letter of the (ordered) alphabet.

        - ``codomain`` - (default: None) a combinatorial class of words.
          By default, ``codomain`` is ``self``.

        - ``min_length`` - (default: 1) nonnegative integer. If ``arg`` is
          not specified, then iterate through all the morphisms where the
          length of the images of each letter in the alphabet is at least
          ``min_length``. This is ignored if ``arg`` is a list.

        OUTPUT:

        iterator

        EXAMPLES:

        Iterator over all non-erasing morphisms::

            sage: W = Words('ab')
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

            sage: W = Words('ab')
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

            sage: Y = Words('xyz')
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
            TypeError: codomain (=a) must be an instance of Words_over_OrderedAlphabet

        """
        n = self.size_of_alphabet()
        if min_length < 0:
            min_length = 0
        # create an iterable of compositions (all "compositions" if arg is
        # None, or [arg] otherwise)
        if arg is None:
            # TODO in #17927: use IntegerVectors(length=n, min_part=min_length)
            from sage.combinat.integer_list import IntegerListsNN
            compositions = IntegerListsNN(length=n, min_part=min_length)
        elif isinstance(arg, tuple):
            from sage.combinat.integer_list import IntegerListsLex
            a, b = arg
            compositions = IntegerListsLex(min_sum=a, max_sum=b-1,
                    length=n, min_part=min_length)
        else:
            arg = list(arg)
            if (not len(arg) == n or not
                    all(isinstance(a, (int,Integer)) for a in arg)):
                raise TypeError(
                    "arg (=%s) must be an iterable of %s integers" %(arg, n))
            compositions = [arg]

        # set the codomain
        if codomain is None:
            codomain = self
        elif not isinstance(codomain, Words_over_OrderedAlphabet):
            raise TypeError("codomain (=%s) must be an instance of Words_over_OrderedAlphabet"%codomain)

        # iterate through the morphisms
        from sage.combinat.words.morphism import WordMorphism
        for composition in compositions:
            cuts = [0] + list(composition)
            for i in range(1,len(cuts)):
                cuts[i] += cuts[i-1]
            s = cuts[-1] # same but better than s = sum(composition)
            for big_word in codomain.iterate_by_length(s):
                d = {}
                i = 0
                for a in self.alphabet():
                    d[a] = big_word[cuts[i]:cuts[i+1]]
                    i += 1
                yield WordMorphism(d, codomain=codomain)

    def cmp_letters(self, letter1, letter2):
        r"""
        Returns a negative number, zero or a positive number if
        ``letter1`` < ``letter2``, ``letter1`` == ``letter2`` or
        ``letter1`` > ``letter2`` respectively.

        INPUT:

        - ``letter1`` - a letter in the alphabet
        - ``letter2`` - a letter in the alphabet

        EXAMPLES::

            sage: from sage.combinat.words.words import Words_over_OrderedAlphabet
            sage: from sage.combinat.words.alphabet import build_alphabet
            sage: A = build_alphabet('woa')
            sage: W = Words_over_OrderedAlphabet(A)
            sage: W.cmp_letters('w','a')
            -2
            sage: W.cmp_letters('w','o')
            -1
            sage: W.cmp_letters('w','w')
            0
        """
        return int(self._alphabet.rank(letter1) - self._alphabet.rank(letter2))

class InfiniteWords_over_OrderedAlphabet(Words_over_OrderedAlphabet):
    def __init__(self, alphabet):
        r"""
        EXAMPLES::

            sage: from sage.combinat.words.words import InfiniteWords_over_OrderedAlphabet
            sage: from sage.combinat.words.alphabet import build_alphabet
            sage: A = build_alphabet("abc")
            sage: W = InfiniteWords_over_OrderedAlphabet(A)
            sage: W == loads(dumps(W))
            True
        """
        super(InfiniteWords_over_OrderedAlphabet, self).__init__(alphabet)

    def _repr_(self):
        r"""
        Returns a string representation of self.

        EXAMPLES::

            sage: Words('ab', finite=False)._repr_()
            "Infinite Words over {'a', 'b'}"
        """
        return "Infinite Words over %s" % self.alphabet()

class FiniteWords_over_OrderedAlphabet(Words_over_OrderedAlphabet):
    def __init__(self, alphabet):
        r"""
        EXAMPLES::

            sage: from sage.combinat.words.words import FiniteWords_over_OrderedAlphabet
            sage: from sage.combinat.words.alphabet import build_alphabet
            sage: A = build_alphabet("abc")
            sage: W = FiniteWords_over_OrderedAlphabet(A)
            sage: W == loads(dumps(W))
            True
        """
        super(FiniteWords_over_OrderedAlphabet, self).__init__(alphabet)

    def _repr_(self):
        r"""
        Returns a string representation of self.

        EXAMPLES::

            sage: Words('ab', infinite=False)._repr_()
            "Finite Words over {'a', 'b'}"
        """
        return "Finite Words over %s" % self.alphabet()

    def random_element(self):
        r"""
        Returns a random (infinite) word on the given alphabet.

        EXAMPLE::

            sage: W = Words(5, infinite=False)
            sage: W.random_element()
            Traceback (most recent call last):
            ...
            ValueError: What does it mean to pick a random finite word over a given alphabet?
        """
        raise ValueError("What does it mean to pick a random finite word over a given alphabet?")

    def __iter__(self):
        r"""
        Returns an iterator over all the words of self.

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
            ....:
            word:
            word: 4
            word: 5
            word: 44
            word: 45
            word: 54
            word: 55
            word: 444
            word: 445
            word: 454
            word: 455
            word: 544
            word: 545
            word: 554
            word: 555
            sage: W = Words([5,4], infinite=False)
            sage: for w in W:
            ....:   if len(w)>3:
            ....:       break
            ....:   else:
            ....:       w
            ....:
            word:
            word: 5
            word: 4
            word: 55
            word: 54
            word: 45
            word: 44
            word: 555
            word: 554
            word: 545
            word: 544
            word: 455
            word: 454
            word: 445
            word: 444
        """
        for l in itertools.count():
            for w in self.iterate_by_length(l):
                yield w

class FiniteWords_length_k_over_OrderedAlphabet(FiniteWords_over_OrderedAlphabet):
    def __init__(self, alphabet, length):
        """
        TESTS::

            sage: from sage.combinat.words.words import FiniteWords_length_k_over_OrderedAlphabet
            sage: A = sage.combinat.words.alphabet.build_alphabet([0,1])
            sage: W = FiniteWords_length_k_over_OrderedAlphabet(A, 3)
            sage: W == loads(dumps(W))
            True
        """
        super(FiniteWords_length_k_over_OrderedAlphabet, \
                self).__init__(alphabet)
        self._length = length

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: from sage.combinat.words.words import FiniteWords_length_k_over_OrderedAlphabet
            sage: A = sage.combinat.words.alphabet.build_alphabet([0,1])
            sage: W = FiniteWords_length_k_over_OrderedAlphabet(A, 3)
            sage: [1,2,3] in W
            False
            sage: [1,2] in W
            False
            sage: Words([0,1])([1,0,1]) in W
            True
            sage: Words([1,0])([1,0,1]) in W
            False
            sage: W([1,0,1]) in W
            True
            sage: Word([2,0]) in W
            False
        """
        if super(FiniteWords_length_k_over_OrderedAlphabet, \
                self).__contains__(x) and x.length() == self._length:
            return True
        else:
            return False

    def _repr_(self):
        """
        TESTS::

            sage: from sage.combinat.words.words import FiniteWords_length_k_over_OrderedAlphabet
            sage: A = sage.combinat.words.alphabet.build_alphabet([1,0])
            sage: FiniteWords_length_k_over_OrderedAlphabet(A,3)._repr_()
            'Words of length 3 over {1, 0}'
        """
        from sage.combinat.words.word_options import word_options
        if word_options['old_repr']:
            return "Words over %s of length %s"%(self.alphabet(), self._length)
        return "Words of length %s over %s"%(self._length, self.alphabet())

    def cardinality(self):
        r"""
        Returns the number of words of length `n` from alphabet.

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
        n = self.size_of_alphabet()
        return n**self._length

    def random_element(self):
        r"""
        Returns a random word uniformly.

        The word returned has infinite length, and is built by randomly picking
        a letter from the alphabet, one after the other.

        EXAMPLE::

            sage: W = Words(2,20)
            sage: W.random_element() # random
            word: 11212222211111221112

            sage: W = Words(ZZ,20)
            sage: W.random_element()
            Traceback (most recent call last):
            ...
            ValueError: How can I pick a random word with an infinite aphabet?

        TESTS::

            sage: _ = Words(GF(5),4).random_element()

        Check that :trac:`18283` is fixed::

            sage: w = Words('abc', 5).random_element()
            sage: w.length()
            5
        """
        if self.alphabet().cardinality() == Infinity:
            raise ValueError("How can I pick a random word with an infinite aphabet?")
        else:
            return self((self.alphabet().random_element() for _ in range(self._length)))

    def list(self):
        r"""
        Returns a list of all the words contained in self.

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

    def __iter__(self):
        """
        Returns an iterator for all of the words of length k from
        ``self.alphabet()``. The iterator outputs the words in lexicographic
        order, with respect to the ordering of the alphabet.

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
        return super(FiniteWords_length_k_over_OrderedAlphabet, \
                self).iterate_by_length(self._length)

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
        if length == self._length:
            return iter(self)
        else:
            return iter([])
