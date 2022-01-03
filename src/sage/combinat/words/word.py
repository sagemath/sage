# -*- coding: utf-8 -*-
r"""
Word classes

AUTHORS:

- Arnaud Bergeron
- Amy Glen
- Sébastien Labbé
- Franco Saliola


"""
#*****************************************************************************
#       Copyright (C) 2008 Arnaud Bergeron <abergeron@gmail.com>,
#                          Amy Glen <amy.glen@gmail.com>,
#                          Sébastien Labbé <slabqc@gmail.com>,
#                          Franco Saliola <saliola@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.combinat.words.word_char import WordDatatype_char
from sage.combinat.words.abstract_word import Word_class
from sage.combinat.words.finite_word import FiniteWord_class
from sage.combinat.words.infinite_word import InfiniteWord_class
from .word_datatypes import (WordDatatype_str,
                            WordDatatype_list,
                            WordDatatype_tuple)
from .word_infinite_datatypes import (
                            WordDatatype_iter_with_caching,
                            WordDatatype_iter,
                            WordDatatype_callable_with_caching,
                            WordDatatype_callable)
from .morphic import WordDatatype_morphic
from sage.monoids.free_monoid_element import FreeMonoidElement

# TODO. Word needs to be replaced by Word. Consider renaming
# Word_class to Word and imbedding Word as its __call__ method.

def Word(data=None, alphabet=None, length=None, datatype=None, caching=True, RSK_data=None):
    r"""
    Construct a word.

    INPUT:

    -  ``data`` -- (default: ``None``) list, string, tuple, iterator, free
       monoid element, ``None`` (shorthand for ``[]``), or a callable defined
       on ``[0,1,...,length]``.

    -  ``alphabet`` -- any argument accepted by Words

    -  ``length`` -- (default: ``None``) This is dependent on the type of data.
       It is ignored for words defined by lists, strings, tuples,
       etc., because they have a naturally defined length.
       For callables, this defines the domain of definition,
       which is assumed to be ``[0, 1, 2, ..., length-1]``.
       For iterators: Infinity if you know the iterator will not
       terminate (default); ``"unknown"`` if you do not know whether the
       iterator terminates; ``"finite"`` if you know that the iterator
       terminates, but do not know the length.

    -  ``datatype`` -- (default: ``None``) ``None``, ``"list"``, ``"str"``,
       ``"tuple"``, ``"iter"``, ``"callable"``. If ``None``, then the function
       tries to guess this from the data.

    -  ``caching`` -- (default: ``True``) ``True`` or ``False``. Whether to
       keep a cache of the letters computed by an iterator or callable.

    -  ``RSK_data`` -- (Optional. Default: ``None``) A semistandard and a
       standard Young tableau to run the inverse RSK bijection on.

    .. NOTE::

       Be careful when defining words using callables and iterators. It
       appears that islice does not pickle correctly causing various errors
       when reloading. Also, most iterators do not support copying and
       should not support pickling by extension.

    EXAMPLES:

    Empty word::

        sage: Word()
        word:

    Word with string::

        sage: Word("abbabaab")
        word: abbabaab

    Word with string constructed from other types::

        sage: Word([0,1,1,0,1,0,0,1], datatype="str")
        word: 01101001
        sage: Word((0,1,1,0,1,0,0,1), datatype="str")
        word: 01101001

    Word with list::

        sage: Word([0,1,1,0,1,0,0,1])
        word: 01101001

    Word with list constructed from other types::

        sage: Word("01101001", datatype="list")
        word: 01101001
        sage: Word((0,1,1,0,1,0,0,1), datatype="list")
        word: 01101001

    Word with tuple::

        sage: Word((0,1,1,0,1,0,0,1))
        word: 01101001

    Word with tuple constructed from other types::

        sage: Word([0,1,1,0,1,0,0,1], datatype="tuple")
        word: 01101001
        sage: Word("01101001", datatype="str")
        word: 01101001

    Word with iterator::

        sage: from itertools import count
        sage: Word(count())
        word: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,...
        sage: Word(iter("abbabaab")) # iterators default to infinite words
        word: abbabaab
        sage: Word(iter("abbabaab"), length="unknown")
        word: abbabaab
        sage: Word(iter("abbabaab"), length="finite")
        word: abbabaab

    Word with function (a 'callable')::

        sage: f = lambda n : add(Integer(n).digits(2)) % 2
        sage: Word(f)
        word: 0110100110010110100101100110100110010110...
        sage: Word(f, length=8)
        word: 01101001

    Word over a string with a parent::

        sage: w = Word("abbabaab", alphabet="abc"); w
        word: abbabaab
        sage: w.parent()
        Finite words over {'a', 'b', 'c'}

    Word from a free monoid element::

        sage: M.<x,y,z> = FreeMonoid(3)
        sage: Word(x^3*y*x*z^2*x)
        word: xxxyxzzx

    The default parent is the combinatorial class of all words::

        sage: w = Word("abbabaab"); w
        word: abbabaab
        sage: w.parent()
        Finite words over Set of Python objects of class 'object'

    We can also input a semistandard tableau and a standard tableau to
    obtain a word from the inverse RSK algorithm using the
    ``RSK_data`` option::

        sage: p = Tableau([[1,2,2],[3]]); q = Tableau([[1,2,4],[3]])
        sage: Word(RSK_data=[p, q])
        word: 1322

    TESTS::

        sage: Word(5)
        Traceback (most recent call last):
        ...
        ValueError: Cannot guess a datatype from data (=5); please specify one

    ::

        sage: W = Words()
        sage: w  = W('abc')
        sage: w is W(w)
        True
        sage: w is Word(w, alphabet='abc')
        False
    """
    if isinstance(data, FreeMonoidElement):
        return data.to_word(alphabet)

    if RSK_data is not None:
        #if a list of a semistandard and a standard tableau or a pair of lists
        from sage.combinat.tableau import Tableau
        if isinstance(RSK_data, (tuple, list)) and len(RSK_data) == 2 and \
            all((isinstance(x, Tableau) for x in RSK_data)):
            from sage.combinat.rsk import RSK_inverse
            return RSK_inverse(*RSK_data, output='word')
        elif isinstance(RSK_data, (tuple, list)) and len(RSK_data) == 2 and \
            all((isinstance(x, (list, tuple)) for x in RSK_data)):
            from sage.combinat.rsk import RSK_inverse
            P,Q = map(Tableau, RSK_data)
            return RSK_inverse(P, Q, 'word')
        raise ValueError("Invalid input. Must be a pair of tableaux")

    # Create the parent object
    from .words import Words
    parent = Words(alphabet)

    return parent(data=data, length=length, datatype=datatype, caching=caching)

#######################################################################
#                                                                     #
#                    Concrete word classes                            #
#                                                                     #
#######################################################################

##### Finite Words #####

class FiniteWord_char(WordDatatype_char, FiniteWord_class):
    r"""
    Finite word represented by an array ``unsigned char *`` (i.e. integers
    between 0 and 255).

    For any word ``w``, type ``w.<TAB>`` to see the functions that can be applied
    to ``w``.

    EXAMPLES::

        sage: W = Words(range(20))

        sage: w = W(list(range(1, 10)) * 2)
        sage: type(w)
        <class 'sage.combinat.words.word.FiniteWord_char'>
        sage: w
        word: 123456789123456789

        sage: w.is_palindrome()
        False
        sage: (w*w[::-1]).is_palindrome()
        True
        sage: (w[:-1:]*w[::-1]).is_palindrome()
        True

        sage: w.is_lyndon()
        False
        sage: W(list(range(10)) + [10, 10]).is_lyndon()
        True

        sage: w.is_square_free()
        False
        sage: w[:-1].is_square_free()
        True

        sage: u = W([randint(0,10) for i in range(10)])
        sage: (u*u).is_square()
        True
        sage: (u*u*u).is_cube()
        True

        sage: len(w.factor_set())
        127
        sage: w.rauzy_graph(5)
        Looped digraph on 9 vertices

        sage: u = W([1,2,3])
        sage: w.first_occurrence(u)
        0
        sage: w.first_occurrence(u, start=1)
        9

    TESTS::

        sage: W = Words([0,1,2])
        sage: w = W([0,1,1,0])
        sage: w == loads(dumps(w))
        True
    """
    pass

class FiniteWord_list(WordDatatype_list, FiniteWord_class):
    r"""
    Finite word represented by a Python list.

    For any word `w`, type ``w.`` and hit TAB key to see the list of
    functions defined on `w`.

    EXAMPLES::

        sage: w = Word(range(10))
        sage: w.iterated_right_palindromic_closure()
        word: 0102010301020104010201030102010501020103...

    TESTS::

        sage: w = Word([0,1,1,0])
        sage: w == loads(dumps(w))
        True
    """
    pass

class FiniteWord_str(WordDatatype_str, FiniteWord_class):
    r"""
    Finite word represented by a Python str.

    For such word `w`, type ``w.`` and hit TAB key to see the list of
    functions defined on `w`.

    EXAMPLES::

        sage: w = Word('abcdef')
        sage: w.is_square()
        False

    TESTS::

        sage: w = Word('abba')
        sage: w == loads(dumps(w))
        True
    """
    pass

class FiniteWord_tuple(WordDatatype_tuple, FiniteWord_class):
    r"""
    Finite word represented by a Python tuple.

    For such word `w`, type ``w.`` and hit TAB key to see the list of
    functions defined on `w`.

    EXAMPLES::

        sage: w = Word(())
        sage: w.is_empty()
        True

    TESTS::

        sage: w = Word((0,1,1,0))
        sage: w == loads(dumps(w))
        True
    """
    pass

class FiniteWord_iter_with_caching(WordDatatype_iter_with_caching, FiniteWord_class):
    r"""
    Finite word represented by an iterator (with caching).

    For such word `w`, type ``w.`` and hit TAB key to see the list of
    functions defined on `w`.

    EXAMPLES::

        sage: w = Word(iter('abcdef'))
        sage: w.conjugate(2)
        word: cdefab

    TESTS::

        sage: w = Word(iter(range(10)))
        sage: type(w)
        <class 'sage.combinat.words.word.FiniteWord_iter_with_caching'>
        sage: z = loads(dumps(w))
        sage: w == z
        True
        sage: type(z)
        <class 'sage.combinat.words.word.FiniteWord_list'>
    """
    pass

class FiniteWord_iter(WordDatatype_iter, FiniteWord_class):
    r"""
    Finite word represented by an iterator.

    For such word `w`, type  ``w.`` and hit TAB key to see the list of
    functions defined on `w`.

    EXAMPLES::

        sage: w = Word(iter(range(10)), caching=False)
        sage: w
        word: 0123456789
        sage: w.finite_differences()
        word: 111111111

    TESTS::

        sage: w = Word(iter(range(10)), caching=False)
        sage: type(w)
        <class 'sage.combinat.words.word.FiniteWord_iter'>
        sage: z = loads(dumps(w))
        sage: w == z
        True
        sage: type(z)
        <class 'sage.combinat.words.word.FiniteWord_list'>
    """
    pass

class FiniteWord_callable_with_caching(WordDatatype_callable_with_caching, FiniteWord_class):
    r"""
    Finite word represented by a callable (with caching).

    For such word `w`, type ``w.`` and hit TAB key to see the list of
    functions defined on `w`.

    EXAMPLES::

        sage: f = lambda n : n % 3
        sage: w = Word(f, length=32)
        sage: w
        word: 01201201201201201201201201201201
        sage: w.border()
        word: 01201201201201201201201201201

    TESTS::

        sage: w = Word(lambda n:n, length=10)
        sage: type(w)
        <class 'sage.combinat.words.word.FiniteWord_callable_with_caching'>
        sage: z = loads(dumps(w))
        sage: w == z
        True
        sage: type(z)
        <class 'sage.combinat.words.word.FiniteWord_callable_with_caching'>

    Pickle also works for concatenation of words::

        sage: w = Word(range(10)) * Word('abcdef')
        sage: type(w)
        <class 'sage.combinat.words.word.FiniteWord_callable_with_caching'>
        sage: z = loads(dumps(w))
        sage: w == z
        True
        sage: type(z)
        <class 'sage.combinat.words.word.FiniteWord_list'>

    Pickle also works for power of words::

        sage: w = Word(range(10)) ^ 2
        sage: type(w)
        <class 'sage.combinat.words.word.FiniteWord_callable_with_caching'>
        sage: z = loads(dumps(w))
        sage: w == z
        True
        sage: type(z)
        <class 'sage.combinat.words.word.FiniteWord_list'>
    """
    pass

class FiniteWord_callable(WordDatatype_callable, FiniteWord_class):
    r"""
    Finite word represented by a callable.

    For such word `w`, type ``w.`` and hit TAB key to see the list of
    functions defined on `w`.

    EXAMPLES::

        sage: f = lambda n : 3 if n > 8 else 6
        sage: w = Word(f, length=30, caching=False)
        sage: w
        word: 666666666333333333333333333333
        sage: w.is_symmetric()
        True

    TESTS::

        sage: w = Word(lambda n:n, length=10, caching=False)
        sage: type(w)
        <class 'sage.combinat.words.word.FiniteWord_callable'>
        sage: z = loads(dumps(w))
        sage: w == z
        True
        sage: type(z)
        <class 'sage.combinat.words.word.FiniteWord_callable'>
    """
    pass

##### Infinite Words #####

class InfiniteWord_iter_with_caching(WordDatatype_iter_with_caching, InfiniteWord_class):
    r"""
    Infinite word represented by an iterable (with caching).

    For such word `w`, type ``w.`` and hit TAB key to see the list of
    functions defined on `w`.

    Infinite words behave like a Python list : they can be sliced using
    square braquets to define for example a prefix or a factor.

    EXAMPLES::

        sage: from itertools import cycle
        sage: w = Word(cycle([9,8,4]))
        sage: w
        word: 9849849849849849849849849849849849849849...
        sage: prefix = w[:23]
        sage: prefix
        word: 98498498498498498498498
        sage: prefix.minimal_period()
        3

    TESTS::

        sage: from itertools import count
        sage: w = Word(count())
        sage: type(w)
        <class 'sage.combinat.words.word.Word_iter_with_caching'>

    Pickle is not supported for infinite word defined by an iterator::

        sage: dumps(w)
        Traceback (most recent call last):
        ...
        TypeError: can...t...pickle...generator...object...
    """
    pass

class InfiniteWord_iter(WordDatatype_iter, InfiniteWord_class):
    r"""
    Infinite word represented by an iterable.

    For such word `w`, type ``w.`` and hit TAB key to see the list of
    functions defined on `w`.

    Infinite words behave like a Python list : they can be sliced using
    square braquets to define for example a prefix or a factor.

    EXAMPLES::

        sage: from itertools import chain, cycle
        sage: w = Word(chain('letsgo', cycle('forever')), caching=False)
        sage: w
        word: letsgoforeverforeverforeverforeverforeve...
        sage: prefix = w[:100]
        sage: prefix
        word: letsgoforeverforeverforeverforeverforeve...
        sage: prefix.is_lyndon()
        False

    TESTS::

        sage: from itertools import count
        sage: w = Word(count(), caching=False)
        sage: type(w)
        <class 'sage.combinat.words.word.Word_iter'>

    Pickle is not supported for infinite word defined by an iterator::

        sage: dumps(w)
        Traceback (most recent call last):
        ...
        TypeError: can...t...pickle...generator...object...
    """
    pass

class InfiniteWord_callable_with_caching(WordDatatype_callable_with_caching, InfiniteWord_class):
    r"""
    Infinite word represented by a callable (with caching).

    For such word `w`, type ``w.`` and hit TAB key to see the list of
    functions defined on `w`.

    Infinite words behave like a Python list : they can be sliced using
    square braquets to define for example a prefix or a factor.

    EXAMPLES::

        sage: w = Word(lambda n:n)
        sage: factor = w[4:13]
        sage: factor
        word: 4,5,6,7,8,9,10,11,12

    TESTS::

        sage: w = Word(lambda n:n)
        sage: type(w)
        <class 'sage.combinat.words.word.InfiniteWord_callable_with_caching'>
        sage: z = loads(dumps(w))
        sage: z
        word: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,...
        sage: type(z)
        <class 'sage.combinat.words.word.InfiniteWord_callable_with_caching'>
    """
    pass

class InfiniteWord_callable(WordDatatype_callable, InfiniteWord_class):
    r"""
    Infinite word represented by a callable.

    For such word `w`, type ``w.`` and hit TAB key to see the list of
    functions defined on `w`.

    Infinite words behave like a Python list : they can be sliced using
    square braquets to define for example a prefix or a factor.

    EXAMPLES::

        sage: w = Word(lambda n:n, caching=False)
        sage: w
        word: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,...
        sage: w.iterated_right_palindromic_closure()
        word: 0102010301020104010201030102010501020103...

    TESTS::

        sage: w = Word(lambda n:n, caching=False)
        sage: type(w)
        <class 'sage.combinat.words.word.InfiniteWord_callable'>
        sage: z = loads(dumps(w))
        sage: z
        word: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,...
        sage: type(z)
        <class 'sage.combinat.words.word.InfiniteWord_callable'>
    """
    pass

##### Words of unknown length #####

class Word_iter_with_caching(WordDatatype_iter_with_caching, Word_class):
    r"""
    Word of unknown length (finite or infinite) represented by an
    iterable (with caching).

    For such word `w`, type ``w.`` and hit TAB key to see the list of
    functions defined on `w`.

    Words behave like a Python list : they can be sliced using
    square braquets to define for example a prefix or a factor.

    EXAMPLES::

        sage: w = Word(iter([1,2,3]*1000), length='unknown')
        sage: w
        word: 1231231231231231231231231231231231231231...
        sage: w.finite_differences(mod=2)
        word: 1101101101101101101101101101101101101101...

    TESTS::

        sage: w = Word(iter('abcd'*100), length='unknown')
        sage: type(w)
        <class 'sage.combinat.words.word.Word_iter_with_caching'>
        sage: w
        word: abcdabcdabcdabcdabcdabcdabcdabcdabcdabcd...

    Pickle is not supported for word of unknown length defined by an iterator::

        sage: dumps(w)
        Traceback (most recent call last):
        ...
        TypeError: can...t...pickle...generator...object...
    """
    pass

class Word_iter(WordDatatype_iter, Word_class):
    r"""
    Word of unknown length (finite or infinite) represented by an
    iterable.

    For such word `w`, type ``w.`` and hit TAB key to see the list of
    functions defined on `w`.

    Words behave like a Python list : they can be sliced using
    square braquets to define for example a prefix or a factor.

    EXAMPLES::

        sage: w = Word(iter([1,1,4,9]*1000), length='unknown', caching=False)
        sage: w
        word: 1149114911491149114911491149114911491149...
        sage: w.delta()
        word: 2112112112112112112112112112112112112112...

    TESTS::

        sage: w = Word(iter('abcd'*100), length='unknown', caching=False)
        sage: type(w)
        <class 'sage.combinat.words.word.Word_iter'>
        sage: w
        word: abcdabcdabcdabcdabcdabcdabcdabcdabcdabcd...

    Pickle is not supported for word of unknown length defined by an iterator::

        sage: dumps(w)
        Traceback (most recent call last):
        ...
        TypeError: can...t...pickle...generator...object...
    """
    pass

##### Morphic Words #####
class FiniteWord_morphic(WordDatatype_morphic, FiniteWord_class):
    r"""
    Finite morphic word.

    For such word `w`, type  ``w.`` and hit TAB key to see the list of
    functions defined on `w`.

    EXAMPLES::

        sage: m = WordMorphism("a->ab,b->")
        sage: w = m.fixed_point("a")
        sage: w
        word: ab


    TESTS::

        sage: m = WordMorphism("a->ab,b->")
        sage: w = m.fixed_point("a")
        sage: type(w)
        <class 'sage.combinat.words.word.FiniteWord_morphic'>
        sage: loads(dumps(w))
        word: ab

    """
    pass

class InfiniteWord_morphic(WordDatatype_morphic, InfiniteWord_class):
    r"""
    Morphic word of infinite length.

    For such word `w`, type ``w.`` and hit TAB key to see the list of
    functions defined on `w`.

    Infinite words behave like a Python list : they can be sliced using
    square braquets to define for example a prefix or a factor.

    EXAMPLES::

        sage: m = WordMorphism('a->ab,b->a')
        sage: w = m.fixed_point('a')
        sage: w
        word: abaababaabaababaababaabaababaabaababaaba...

    TESTS:

    Pickle is supported::

        sage: loads(dumps(w))
        word: abaababaabaababaababaabaababaabaababaaba...

    """
    pass
