# -*- coding: utf-8 -*-
r"""
Alphabet

AUTHORS:

- Franco Saliola (2008-12-17) : merged into sage
- Vincent Delecroix and Stepan Starosta (2012): remove classes for alphabet and
  use other Sage classes otherwise (TotallyOrderedFiniteSet,
  FiniteEnumeratedSet, ...). More shortcut to standard alphabets.

EXAMPLES::

    sage: build_alphabet("ab")
    {'a', 'b'}
    sage: build_alphabet([0,1,2])
    {0, 1, 2}
    sage: build_alphabet(name="PP")
    Positive integers
    sage: build_alphabet(name="NN")
    Non negative integers
    sage: build_alphabet(name="lower")
    {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z'}
"""
# ****************************************************************************
#       Copyright (C) 2008 Franco Saliola <saliola@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.sets_cat import Sets

from sage.sets.totally_ordered_finite_set import TotallyOrderedFiniteSet
from sage.sets.family import Family

from sage.rings.integer import Integer
from sage.rings.infinity import Infinity

from sage.sets.non_negative_integers import NonNegativeIntegers


set_of_letters = {
    'lower'       : "abcdefghijklmnopqrstuvwxyz",
    'upper'       : "ABCDEFGHIJKLMNOPQRSTUVWXYZ",
    'space'       : " ",
    'underscore'  : "_",
    'punctuation' : " ,.;:!?",
    'printable'   : "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~",
    'binary'      : "01",
    'octal'       : "01234567",
    'decimal'     : "0123456789",
    'hexadecimal' : "0123456789abcdef",
    'radix64'     : "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"}


def build_alphabet(data=None, names=None, name=None):
    r"""
    Return an object representing an ordered alphabet.

    INPUT:

    - ``data`` -- the letters of the alphabet; it can be:

      * a list/tuple/iterable of letters; the iterable may be infinite
      * an integer `n` to represent `\{1, \ldots, n\}`, or infinity to
        represent `\NN`

    - ``names`` -- (optional) a list for the letters (i.e. variable names) or
      a string for prefix for all letters; if given a list, it must have the
      same cardinality as the set represented by ``data``

    - ``name`` -- (optional) if given, then return a named set and can be
      equal to : ``'lower', 'upper', 'space',
      'underscore', 'punctuation', 'printable', 'binary', 'octal', 'decimal',
      'hexadecimal', 'radix64'``.

      You can use many of them at once, separated by spaces : ``'lower
      punctuation'`` represents the union of the two alphabets ``'lower'`` and
      ``'punctuation'``.

      Alternatively, ``name`` can be set to ``"positive integers"`` (or
      ``"PP"``) or ``"natural numbers"`` (or ``"NN"``).

      ``name`` cannot be combined with ``data``.

    EXAMPLES:

    If the argument is a Set, it just returns it::

        sage: build_alphabet(ZZ) is ZZ
        True
        sage: F = FiniteEnumeratedSet('abc')
        sage: build_alphabet(F) is F
        True

    If a list, tuple or string is provided, then it builds a proper Sage class
    (:class:`~sage.sets.totally_ordered_finite_set.TotallyOrderedFiniteSet`)::

        sage: build_alphabet([0,1,2])
        {0, 1, 2}
        sage: F = build_alphabet('abc'); F
        {'a', 'b', 'c'}
        sage: print(type(F).__name__)
        TotallyOrderedFiniteSet_with_category

    If an integer and a set is given, then it constructs a
    :class:`~sage.sets.totally_ordered_finite_set.TotallyOrderedFiniteSet`::

        sage: build_alphabet(3, ['a','b','c'])
        {'a', 'b', 'c'}

    If an integer and a string is given, then it considers that string as a
    prefix::

        sage: build_alphabet(3, 'x')
        {'x0', 'x1', 'x2'}

    If no data is provided, ``name`` may be a string which describe an alphabet.
    The available names decompose into two families. The first one are 'positive
    integers', 'PP', 'natural numbers' or 'NN' which refer to standard set of
    numbers::

        sage: build_alphabet(name="positive integers")
        Positive integers
        sage: build_alphabet(name="PP")
        Positive integers
        sage: build_alphabet(name="natural numbers")
        Non negative integers
        sage: build_alphabet(name="NN")
        Non negative integers

    The other families for the option ``name`` are among 'lower', 'upper',
    'space', 'underscore', 'punctuation', 'printable', 'binary', 'octal',
    'decimal', 'hexadecimal', 'radix64' which refer to standard set of
    characters. Theses names may be combined by separating them by a space::

        sage: build_alphabet(name="lower")
        {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z'}
        sage: build_alphabet(name="hexadecimal")
        {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd', 'e', 'f'}
        sage: build_alphabet(name="decimal punctuation")
        {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', ' ', ',', '.', ';', ':', '!', '?'}

    In the case the alphabet is built from a list or a tuple, the order on the
    alphabet is given by the elements themselves::

        sage: A = build_alphabet([0,2,1])
        sage: A(0) < A(2)
        True
        sage: A(2) < A(1)
        False

    If a different order is needed, you may use
    :class:`~sage.sets.totally_ordered_finite_set.TotallyOrderedFiniteSet` and
    set the option ``facade`` to ``False``. That way, the comparison fits the
    order of the input::

        sage: A = TotallyOrderedFiniteSet([4,2,6,1], facade=False)
        sage: A(4) < A(2)
        True
        sage: A(1) < A(6)
        False

    Be careful, the element of the set in the last example are no more
    integers and do not compare equal with integers::

        sage: type(A.an_element())
        <class 'sage.sets.totally_ordered_finite_set.TotallyOrderedFiniteSet_with_category.element_class'>
        sage: A(1) == 1
        False
        sage: 1 == A(1)
        False

    We give an example of an infinite alphabet indexed by the positive
    integers and the prime numbers::

        sage: build_alphabet(oo, 'x')
        Lazy family (x(i))_{i in Non negative integers}
        sage: build_alphabet(Primes(), 'y')
        Lazy family (y(i))_{i in Set of all prime numbers: 2, 3, 5, 7, ...}

    TESTS::

        sage: Alphabet(3, name="punctuation")
        Traceback (most recent call last):
        ...
        ValueError: name cannot be specified with any other argument
        sage: Alphabet(8, ['e']*10)
        Traceback (most recent call last):
        ...
        ValueError: invalid value for names
        sage: Alphabet(8, x)
        Traceback (most recent call last):
        ...
        ValueError: invalid value for names
        sage: Alphabet(name=x, names="punctuation")
        Traceback (most recent call last):
        ...
        ValueError: name cannot be specified with any other argument
        sage: Alphabet(x)
        Traceback (most recent call last):
        ...
        ValueError: unable to construct an alphabet from the given parameters
    """
    # If both 'names' and 'data' are defined
    if name is not None and (data is not None or names is not None):
        raise ValueError("name cannot be specified with any other argument")

    # Swap arguments if we need to try and make sure we have "good" user input
    if isinstance(names, (int, Integer)) or names == Infinity \
            or (data is None and names is not None):
        data, names = names, data

    # data is an integer
    if isinstance(data, (int, Integer)):
        if names is None:
            from sage.sets.integer_range import IntegerRange
            return IntegerRange(Integer(data))
        if isinstance(names, str):
            return TotallyOrderedFiniteSet([names + '%d' % i for i in range(data)])
        if len(names) == data:
            return TotallyOrderedFiniteSet(names)
        raise ValueError("invalid value for names")

    if data == Infinity:
        data = NonNegativeIntegers()

    # data is an iterable
    if isinstance(data, (tuple, list, str, range)) or data in Sets():
        if names is not None:
            if not isinstance(names, str):
                raise TypeError("names must be a string when data is a set")
            return Family(data, lambda i: names + str(i), name=names)
        if data in Sets():
            return data
        return TotallyOrderedFiniteSet(data)

    # Alphabet defined from a name
    if name is not None:
        if not isinstance(name, str):
            raise TypeError("name must be a string")
        if name == "positive integers" or name == "PP":
            from sage.sets.positive_integers import PositiveIntegers
            return PositiveIntegers()
        if name == "natural numbers" or name == "NN":
            return NonNegativeIntegers()

        data = []
        for alpha_name in name.split(' '):
            try:
                data.extend(list(set_of_letters[alpha_name]))
            except KeyError:
                raise TypeError("name is not recognized")
        return TotallyOrderedFiniteSet(data)

    # Alphabet(**nothing**)
    if data is None:  # name is also None
        from sage.sets.pythonclass import Set_PythonType
        return Set_PythonType(object)

    raise ValueError("unable to construct an alphabet from the given parameters")


# TODO: should it be deprecated as it is no more a class ?
Alphabet = build_alphabet
