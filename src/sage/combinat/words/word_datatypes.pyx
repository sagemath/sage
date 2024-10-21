# sage_setup: distribution = sagemath-combinat
r"""
Datatypes for finite words
"""
# ****************************************************************************
#       Copyright (C) 2009 Franco Saliola <saliola@gmail.com>
#                          Vincent Delecroix <20100.delecroix@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cpython.object cimport Py_EQ, Py_NE
from itertools import islice


cdef class WordDatatype():
    r"""
    The generic WordDatatype class.

    Any word datatype must contain two attributes (at least):

    - ``_parent``
    - ``_hash``

    They are automatically defined here and it's not necessary (and forbidden)
    to define them anywhere else.

    TESTS::

        sage: w = Word([0,1,1,0,0,1])
        sage: isinstance(w, sage.combinat.words.word_datatypes.WordDatatype)
        True
    """
    def __reduce__(self):
        r"""
        Default pickle support.

        TESTS::

            sage: w = Word([0,1,1,0,0,1])
            sage: w.__reduce__()
            (Finite words over Set of Python objects of class 'object', ([0, 1, 1, 0, 0, 1],))
        """
        return self._parent, (list(self),)

    def __hash__(self):
        r"""
        Return the hash for this word.

        TESTS::

             sage: h = hash(Word('abc'))    # indirect test
             sage: Word('abc').__hash__() == Word('abc').__hash__()
             True

             sage: tm = words.ThueMorseWord()
             sage: hash(tm)
             -973965563
        """
        cdef int res
        if self._hash is None:
            res = 5381
            for s in islice(self, 1024):
                res = ((res << 5) + res) + hash(s)
            self._hash = res
        return self._hash


cdef class WordDatatype_list(WordDatatype):
    r"""
    Datatype class for words defined by lists.
    """
    cdef public list _data

    def __init__(self, parent=None, data=None):
        r"""
        Construct a word with a given parent.

        .. NOTE::

           It is slower than WordDatatype_str and WordDatatype_tuple.

        INPUT:

        - ``parent`` -- an instance of :class:`Words_all`
        - ``data`` -- an iterable

        EXAMPLES::

            sage: w = Word([0,1,1,0])
            sage: isinstance(w, sage.combinat.words.word_datatypes.WordDatatype_list)
            True
        """
        self._parent = parent
        if isinstance(data, list):
            self._data = data
        else:
            self._data = list(data)
        self._hash = None

    __hash__ = WordDatatype.__hash__

    def __contains__(self, a):
        r"""
        Test whether ``a`` is a letter of ``self``.

        INPUT:

        - ``a`` -- anything

        OUTPUT: boolean

        EXAMPLES::

            sage: w = Word([0,1,1,0])
            sage: 0 in w
            True
            sage: 3 in w
            False
        """
        return a in self._data

    def __iter__(self):
        r"""
        Return an iterator that iterates through the letters of ``self``.

        EXAMPLES::

            sage: w = Word([0,1,1,0])
            sage: list(iter(w))
            [0, 1, 1, 0]
        """
        return iter(self._data)

    def __richcmp__(self, other, int op):
        r"""
        Equality test for ``self`` and ``other`` if other is an instance of
        ``WordDatype_list``.

        INPUT:

        - ``other`` -- a word
        - ``op`` -- integer; 0, 1, 2, 3, 4 or 5

        OUTPUT: boolean or NotImplemented

        EXAMPLES::

            sage: w = Word(range(10))
            sage: w == w
            True
            sage: z = Word(range(20))
            sage: w == z
            False
            sage: z == w
            False

        It works even if the parents are not the same::

            sage: Words([0,1])([0,1,1]) == Words([0,1,2])([0,1,1])
            True

        REFERENCES:

        http://docs.cython.org/docs/special_methods.html
        """
        if isinstance(other, WordDatatype_list):
            if op == Py_EQ:
                return self._data == other._data
            elif op == Py_NE:
                return self._data != other._data

        # Otherwise, force FiniteWord_class.__richcmp__ to do it
        from sage.combinat.words.word import FiniteWord_class
        return FiniteWord_class.__richcmp__(self, other, op)

    def __len__(self):
        r"""
        Return the length of the word.

        .. NOTE::

           This function will be deprecated in a future version
           of Sage. Use ``self.length()`` instead.

        EXAMPLES::

            sage: w = Word([0,1,1,0])
            sage: len(w)
            4
        """
        return len(self._data)

    def length(self):
        r"""
        Return the length of the word.

        EXAMPLES::

            sage: w = Word([0,1,1,0])
            sage: w.length()
            4
        """
        return len(self._data)

    def __getitem__(self, key):
        r"""
        Implement :meth:`__getitem__` for words stored as lists.

        INPUT:

        - ``key`` -- integer

        EXAMPLES::

            sage: L = list(range(100))
            sage: w = Word(L)
            sage: w[4]
            4
            sage: w[-1]
            99
            sage: w[3:10:2]
            word: 3579
        """
        if isinstance(key, slice):
            return self._parent(self._data[key])
        else:
            return self._data[key]

    def __mul__(self, other):
        r"""
        Return the concatenation of ``self`` and ``other``.

        INPUT:

        - ``other`` -- word represented by a list

        OUTPUT: word

        EXAMPLES::

            sage: w = Word(list(range(10)))
            sage: w * w
            word: 01234567890123456789

        The type of the concatenation is preserved::

            sage: type(w)
            <class 'sage.combinat.words.word.FiniteWord_list'>
            sage: type(w * w)
            <class 'sage.combinat.words.word.FiniteWord_list'>
        """
        if isinstance(other, WordDatatype_list):
            return self._parent(self._data + other._data)
        else:
            return super().__mul__(other)

    __add__ = __mul__

    def number_of_letter_occurrences(self, a):
        r"""
        Return the number of occurrences of the letter ``a`` in the word
        ``self``.

        INPUT:

        - ``a`` -- a letter

        OUTPUT: integer

        EXAMPLES::

            sage: w = Word([0,1,1,0,1])
            sage: w.number_of_letter_occurrences(0)
            2
            sage: w.number_of_letter_occurrences(1)
            3
            sage: w.number_of_letter_occurrences(2)
            0

        .. SEEALSO::

            :meth:`sage.combinat.words.finite_word.FiniteWord_class.number_of_factor_occurrences`
        """
        return self._data.count(a)

cdef class WordDatatype_str(WordDatatype):
    """
    Datatype for words defined by strings.
    """
    cdef public str _data

    # TODO : allow initialization from non string data
    def __init__(self, parent=None, data=None):
        r"""
        Construct a word with parent ``parent`` from the string ``data``.

        INPUT:

        - ``parent`` -- instance of :class:`Words_all`
        - ``data`` -- string

        EXAMPLES::

            sage: w = Word("abba")
            sage: isinstance(w, sage.combinat.words.word_datatypes.WordDatatype_str)
            True
        """
        self._parent = parent
        if isinstance(data, str):
            self._data = data
        else:
            self._data = "".join(str(u) for u in data)
        self._hash = None

    __hash__ = WordDatatype.__hash__

    def __iter__(self):
        r"""
        Return an iterator that iterates through the letters of ``self``.

        EXAMPLES::

            sage: w = Word('abba')
            sage: list(iter(w))
            ['a', 'b', 'b', 'a']
        """
        return iter(self._data)

    def __richcmp__(self, other, int op):
        r"""
        Equality test for ``self`` and ``other`` if other is an instance of
        ``WordDatype_str``.

        INPUT:

        - ``other`` -- a word
        - ``op`` -- integer; 0, 1, 2, 3, 4 or 5

        OUTPUT: boolean or NotImplemented

        EXAMPLES::

            sage: w = Word('abcde')
            sage: w == w
            True
            sage: z = Word('epoisudfafgh')
            sage: w == z
            False
            sage: z == w
            False

        It works even if the parents are not the same::

            sage: Words('ab')('ababa') == Words('abcd')('ababa')
            True
            sage: Words('ab')('ababa') == Word('ababa')
            True

        REFERENCES:

        http://docs.cython.org/docs/special_methods.html
        """
        if isinstance(other, WordDatatype_str):
            if op == Py_EQ:
                return self._data == other._data
            elif op == Py_NE:
                return self._data != other._data

        # Otherwise, force FiniteWord_class.__richcmp__ to do it
        from sage.combinat.words.word import FiniteWord_class
        return FiniteWord_class.__richcmp__(self, other, op)

    def __contains__(self, a):
        r"""
        Test whether ``a`` is a letter of ``self``.

        INPUT:

        - ``a`` -- anything

        EXAMPLES::

            sage: w = Word('abba')
            sage: 'a' in w
            True
            sage: 'c' in w
            False
        """
        # we need to override the non standard behaviour of
        # the __contains__ of python str
        if not isinstance(a, str):
            return False
        if len(a) != 1:
            return False
        else:
            return a in self._data

    cpdef _has_factor_naive(self, w):
        r"""
        A naive test for testing whether the word contains ``w`` as a factor.

        .. NOTE::

           This just wraps Python's builtin :meth:`__contains__` for :class:`str`.

        INPUT:

        - ``w`` -- a word, or something that behaves like one (``list``,
          ``tuple``, ``str``, ...)

        OUTPUT: boolean

        EXAMPLES::

            sage: w = Word('abba')
            sage: w._has_factor_naive('ba')
            True
            sage: w._has_factor_naive('bab')
            False
        """
        if isinstance(w, WordDatatype_str):
            return w._data in self._data
        elif isinstance(w, str):
            return w in self._data
        raise ValueError

    cpdef find(self, sub, start=0, end=None):
        r"""
        Return the index of the first occurrence of sub in self,
        such that sub is contained within self[start:end].
        Returns -1 on failure.

        INPUT:

        - ``sub`` -- string or word to search for
        - ``start`` -- nonnegative integer (default: 0) specifying
          the position from which to start the search.
        - ``end`` -- nonnegative integer (default: ``None``); specifying
          the position at which the search must stop. If ``None``, then
          the search is performed up to the end of the string.

        OUTPUT: nonnegative integer or `-1`

        EXAMPLES::

            sage: w = Word("abbabaabababa")
            sage: w.find("a")
            0
            sage: w.find("a", 4)
            5
            sage: w.find("a", 4, 5)
            -1
        """
        if end is None:
            end = len(self._data)
        if isinstance(sub, WordDatatype_str):
            return self._data.find(sub._data, start, end)
        elif isinstance(sub, str):
            return self._data.find(sub, start, end)
        else:
            return super().find(sub, start, end)

    def rfind(self, sub, start=0, end=None):
        r"""
        Return the index of the last occurrence of sub in self,
        such that sub is contained within self[start:end].
        Returns -1 on failure.

        INPUT:

        - ``sub`` -- string or word to search for
        - ``start`` -- nonnegative integer (default: 0) specifying
          the position at which the search must stop.
        - ``end`` -- nonnegative integer (default: ``None``); specifying
          the position from which to start the search. If ``None``, then
          the search is performed up to the end of the string.

        OUTPUT: nonnegative integer or `-1`

        EXAMPLES::

            sage: w = Word("abbabaabababa")
            sage: w.rfind("a")
            12
            sage: w.rfind("a", 4, 8)
            6
            sage: w.rfind("a", 4, 5)
            -1
        """
        if end is None:
            end = len(self._data)
        if isinstance(sub, WordDatatype_str):
            return self._data.rfind(sub._data, start, end)
        elif isinstance(sub, str):
            return self._data.rfind(sub, start, end)
        else:
            return super().rfind(sub, start, end)

    def __len__(self):
        r"""
        Return the length of the word.

        .. NOTE::

           This function will be deprecated in a future version
           of Sage. Use ``self.length()`` instead.

        EXAMPLES::

            sage: w = Word("abbabaabababa")
            sage: len(w)
            13
        """
        return len(self._data)

    def length(self):
        r"""
        Return the length of the word.

        EXAMPLES::

            sage: w = Word("abbabaabababa")
            sage: w.length()
            13
        """
        return len(self._data)

    def __getitem__(self, key):
        r"""
        Implement the :meth:`__getitem__`.

        TESTS::

            sage: alphabet = [chr(i) for i in range(97, 123)]
            sage: w = Word(alphabet)
            sage: w[4]
            'e'
            sage: w[-1]
            'z'
            sage: w[3:10:2]
            word: dfhj
            sage: all(chr(i+97) == w[i] for i in range(w.length()))
            True
        """
        if isinstance(key, slice):
            return self._parent(self._data[key])
        return self._data[key]

    def __mul__(self, other):
        r"""
        Return the concatenation of ``self`` and ``other``.

        INPUT:

        - ``other`` -- word represented by a string

        OUTPUT: word

        EXAMPLES::

            sage: w = Word('abcdef')
            sage: w * w
            word: abcdefabcdef

        The type of the concatenation is preserved::

            sage: type(w)
            <class 'sage.combinat.words.word.FiniteWord_str'>
            sage: type(w * w)
            <class 'sage.combinat.words.word.FiniteWord_str'>
        """
        if isinstance(other, WordDatatype_str):
            return self._parent(self._data + other._data)
        else:
            return super().__mul__(other)

    __add__ = __mul__

    def number_of_letter_occurrences(self, letter):
        r"""
        Count the number of occurrences of ``letter``.

        INPUT:

        - ``letter`` -- a letter

        OUTPUT: integer

        EXAMPLES::

            sage: w = Word("abbabaabababa")
            sage: w.number_of_letter_occurrences('a')
            7
            sage: w.number_of_letter_occurrences('b')
            6
            sage: w.number_of_letter_occurrences('c')
            0

        ::

            sage: w.number_of_letter_occurrences('abb')
            0

        .. SEEALSO::

            :meth:`sage.combinat.words.finite_word.FiniteWord_class.number_of_factor_occurrences`
        """
        if len(letter) == 1:
            return self._data.count(letter)
        else:
            return 0

    def split(self, sep=None, maxsplit=None):
        r"""
        Return a list of words, using sep as a delimiter string.
        If maxsplit is given, at most maxsplit splits are done.

        See also the partition method.

        .. NOTE::

           This just wraps Python's builtin :meth:`str::split` for
           :class:`str`.

        INPUT:

        - ``sep`` -- string or word (default: ``None``)

        - ``maxsplit`` -- positive integer (default: ``None``)

        OUTPUT: list of words

        EXAMPLES:

        You can split along white space to find words in a sentence::

            sage: w = Word("My tailor is poor")
            sage: w.split(" ")
            [word: My, word: tailor, word: is, word: poor]

        The python behavior is kept when no argument is given::

            sage: w.split()
            [word: My, word: tailor, word: is, word: poor]

        You can split in two words letters to get the length of blocks in the
        other letter::

            sage: w = Word("ababbabaaba")
            sage: w.split('a')
            [word: , word: b, word: bb, word: b, word: , word: b, word: ]
            sage: w.split('b')
            [word: a, word: a, word: , word: a, word: aa, word: a]

        You can split along words::

            sage: w = Word("3230301030323212323032321")
            sage: w.split("32")
            [word: , word: 30301030, word: , word: 12, word: 30, word: , word: 1]

        If the separator is not a string a :exc:`ValueError` is raised::

            sage: w = Word("le papa du papa du papa etait un petit pioupiou")
            sage: w.split(Word(['p','a','p','a']))
            Traceback (most recent call last):
            ...
            ValueError: the separator must be a string
        """
        if sep is None or isinstance(sep, str):
            pass
        elif isinstance(sep, WordDatatype_str):
            sep = sep._data
        else:
            raise ValueError("the separator must be a string")
        if maxsplit is None:
            maxsplit = -1
        return [self._parent(z) for z in self._data.split(sep=sep,
                                                          maxsplit=maxsplit)]

    def partition(self, sep):
        r"""
        Search for the separator sep in S, and return the part before it,
        the separator itself, and the part after it. The concatenation of
        the terms in the list gives back the initial word.

        See also the split method.

        .. NOTE::

           This just wraps Python's builtin :meth:`str::partition` for
           :class:`str`.

        INPUT:

        - ``sep`` -- string or word

        EXAMPLES::

            sage: w = Word("MyTailorIsPoor")
            sage: w.partition("Tailor")
            [word: My, word: Tailor, word: IsPoor]

        ::

            sage: w = Word("3230301030323212323032321210121232121010")
            sage: l = w.partition("323")
            sage: print(l)
            [word: , word: 323, word: 0301030323212323032321210121232121010]
            sage: sum(l, Word('')) == w
            True

        If the separator is not a string an error is raised::

            sage: w = Word("le papa du papa du papa etait un petit pioupiou")
            sage: w.partition(Word(['p','a','p','a']))
            Traceback (most recent call last):
            ...
            ValueError: the separator must be a string
        """
        if isinstance(sep, str):
            return [self._parent(z) for z in self._data.partition(sep)]
        elif isinstance(sep, WordDatatype_str):
            return [self._parent(z) for z in self._data.partition(sep._data)]
        raise ValueError("the separator must be a string")

    def is_suffix(self, other):
        r"""
        Test whether ``self`` is a suffix of ``other``.

        INPUT:

        - ``other`` -- a word (an instance of :class:`Word_class`) or a
          :class:`str`

        OUTPUT: boolean

        EXAMPLES::

            sage: w = Word("abbabaabababa")
            sage: u = Word("ababa")
            sage: w.is_suffix(u)
            False
            sage: u.is_suffix(w)
            True
            sage: u.is_suffix("abbabaabababa")
            True

        TESTS::

            sage: w = Word("abbabaabababa")
            sage: u = Word(['a','b','a','b','a'])
            sage: w.is_suffix(u)
            False
            sage: u.is_suffix(w)
            True
        """
        if isinstance(other, WordDatatype_str):
            return other._data.endswith(self._data)
        elif isinstance(other, str):
            return other.endswith(self._data)
        else:
            return super().is_suffix(other)

    def has_suffix(self, other):
        """
        Test whether ``self`` has ``other`` as a suffix.

        INPUT:

        - ``other`` -- a word (an instance of :class:`Word_class`) or a
          :class:`str`

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
        """
        if isinstance(other, WordDatatype_str):
            return self._data.endswith(other._data)
        elif isinstance(other, str):
            return self._data.endswith(other)
        else:
            return super().has_suffix(other)

    def is_prefix(self, other):
        r"""
        Test whether ``self`` is a prefix of ``other``.

        INPUT:

        - ``other`` -- a word (an instance of :class:`Word_class`) or a
          :class:`str`

        OUTPUT: boolean

        EXAMPLES::

            sage: w = Word("abbabaabababa")
            sage: u = Word("abbab")
            sage: w.is_prefix(u)
            False
            sage: u.is_prefix(w)
            True
            sage: u.is_prefix("abbabaabababa")
            True

        TESTS::

            sage: ab = Word('ab')
            sage: abba = Word(['a','b','b','a'])
            sage: ab.is_prefix(abba)
            True
            sage: abba.is_prefix(ab)
            False
        """
        if isinstance(other, WordDatatype_str):
            return other._data.startswith(self._data)
        if isinstance(other, str):
            return other.startswith(self._data)
        return super().is_prefix(other)

    def has_prefix(self, other):
        r"""
        Test whether ``self`` has ``other`` as a prefix.

        INPUT:

        - ``other`` -- a word (an instance of :class:`Word_class`) or a
          :class:`str`

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

        TESTS::

            sage: ab = Word('ab')
            sage: abba = Word(['a','b','b','a'])
            sage: ab.has_prefix(abba)
            False
            sage: abba.has_prefix(ab)
            True
        """
        if isinstance(other, WordDatatype_str):
            return self._data.startswith(other._data)
        if isinstance(other, str):
            return self._data.startswith(other)
        else:
            return super().has_prefix(other)

cdef class WordDatatype_tuple(WordDatatype):
    r"""
    Datatype class for words defined by tuples.
    """
    cdef public tuple _data

    def __init__(self, parent=None, data=None):
        r"""
        Construct a word with parent ``parent`` from an iterable ``data``.

        INPUT:

        - ``parent`` -- instance of :class:`Words_all`
        - ``data`` -- iterable

        EXAMPLES::

            sage: w = Word((0,1,1,0))
            sage: isinstance(w, sage.combinat.words.word_datatypes.WordDatatype_tuple)
            True
            sage: u = Word([0,1,1,0], datatype='tuple')
            sage: isinstance(u, sage.combinat.words.word_datatypes.WordDatatype_tuple)
            True
        """
        self._parent = parent
        if isinstance(data, tuple):
            self._data = data
        else:
            self._data = tuple(data)
        self._hash = None

    __hash__ = WordDatatype.__hash__

    def __iter__(self):
        r"""
        Return an iterator that iterates through the letters of ``self``.

        EXAMPLES::

            sage: w = Word((0,1,1,0))
            sage: list(iter(w))
            [0, 1, 1, 0]
        """
        return iter(self._data)

    def __richcmp__(self, other, int op):
        r"""
        Equality test for ``self`` and ``other`` if other is an instance of
        ``WordDatype_tuple``.

        INPUT:

        - ``other`` -- a word
        - ``op`` -- integer; 0, 1, 2, 3, 4 or 5

        OUTPUT: boolean or NotImplemented

        EXAMPLES::

            sage: Word((1,2,3)) == Word((1,2,3))
            True
            sage: Word((1,2,3)) == Word(())
            False
            sage: Word((1,2,3)) == Word((1,2,3,4))
            False
            sage: Word((1,2,3)) == Word((1,2,3,'a'))
            False

        It works even if the parents are not the same::

            sage: Words([1,2])((1,1,1,2)) == Words([1,2,3])((1,1,1,2))
            True
            sage: Words([1,2])((1,1,1,2)) == Word((1,1,1,2))
            True

        REFERENCES:

        http://docs.cython.org/docs/special_methods.html
        """
        if isinstance(other, WordDatatype_tuple):
            if op == Py_EQ:
                return self._data == other._data
            elif op == Py_NE:
                return self._data != other._data

        # Otherwise, force FiniteWord_class.__richcmp__ to do it
        from sage.combinat.words.word import FiniteWord_class
        return FiniteWord_class.__richcmp__(self, other, op)

    def __len__(self):
        r"""
        Return the length of the word.

        .. NOTE::

           This function will be deprecated in a future version
           of Sage. Use ``self.length()`` instead.

        EXAMPLES::

            sage: w = Word((0,1,1,0))
            sage: len(w)
            4
        """
        return len(self._data)

    def length(self):
        r"""
        Return the length of the word.

        EXAMPLES::

            sage: w = Word((0,1,1,0))
            sage: w.length()
            4
        """
        return len(self._data)

    def __contains__(self, a):
        r"""
        Test whether ``a`` is a letter of ``self``.

        INPUT:

        - ``a`` -- anything

        EXAMPLES::

            sage: w = Word((0,1,1,0))
            sage: 0 in w
            True
            sage: 3 in w
            False
        """
        return a in self._data

    def __getitem__(self, key):
        r"""
        Implement ``__getitem__`` for words stored as tuples.

        INPUT:

        - ``key`` -- integer

        OUTPUT:

        - can be anything (an object contained in the word)

        EXAMPLES::

            sage: w = Word(tuple(range(100)))
            sage: w[4]
            4
            sage: w[-1]
            99
            sage: w[3:10:2]
            word: 3579
            sage: all(w[i] == i for i in range(100))
            True
        """
        if isinstance(key, slice):
            return self._parent(self._data[key])
        return self._data[key]

    def __mul__(self, other):
        r"""
        Return the concatenation of ``self`` and ``other``.

        INPUT:

        - ``other`` -- word represented by a tuple

        OUTPUT: word

        EXAMPLES::

            sage: w = Word((1,2,3,4))
            sage: w * w
            word: 12341234

        The type of the concatenation is preserved::

            sage: type(w)
            <class 'sage.combinat.words.word.FiniteWord_tuple'>
            sage: type(w * w)
            <class 'sage.combinat.words.word.FiniteWord_tuple'>
            sage: type(w + w)
            <class 'sage.combinat.words.word.FiniteWord_tuple'>
        """
        if isinstance(other, WordDatatype_tuple):
            return self._parent(self._data + other._data)
        else:
            return super().__mul__(other)

    __add__ = __mul__
