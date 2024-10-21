# sage_setup: distribution = sagemath-categories
r"""
Rankers
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#                          Nicolas M. Thiery <nthiery at users.sf.net>
#  Ported from MuPAD-Combinat (combinat::rankers)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from collections.abc import Iterable, Sequence
from sage.misc.cachefunc import cached_function
from sage.misc.callable_dict import CallableDict
from sage.structure.parent import Parent
from sage.categories.enumerated_sets import EnumeratedSets


def from_list(l):
    """
    Return a ranker from the list l.

    INPUT:

    - ``l`` -- list

    OUTPUT: ``[rank, unrank]`` -- functions

    EXAMPLES::

        sage: import sage.combinat.ranker as ranker
        sage: l = [1,2,3]
        sage: r,u = ranker.from_list(l)
        sage: r(1)
        0
        sage: r(3)
        2
        sage: u(2)
        3
        sage: u(0)
        1
    """
    return [rank_from_list(l), unrank_from_list(l)]


def rank_from_list(l):
    r"""
    Return a rank function for the elements of ``l``.

    INPUT:

    - ``l`` -- a duplicate free list (or iterable) of hashable objects

    OUTPUT:

    - a function from the elements of ``l`` to ``0,...,len(l)``

    EXAMPLES::

        sage: import sage.combinat.ranker as ranker
        sage: l = ['a', 'b', 'c']
        sage: r = ranker.rank_from_list(l)
        sage: r('a')
        0
        sage: r('c')
        2

    For non elements a :exc:`ValueError` is raised, as with the usual
    ``index`` method of lists::

        sage: r('blah')
        Traceback (most recent call last):
        ...
        ValueError: 'blah' is not in dict

    Currently, the rank function is a
    :class:`~sage.misc.callable_dict.CallableDict`; but this is an
    implementation detail::

        sage: type(r)
        <class 'sage.misc.callable_dict.CallableDict'>
        sage: r
        {'a': 0, 'b': 1, 'c': 2}

    With the current implementation, no error is issued in case of
    duplicate value in ``l``. Instead, the rank function returns the
    position of some of the duplicates::

        sage: r = ranker.rank_from_list(['a', 'b', 'a', 'c'])
        sage: r('a')
        2

    Constructing the rank function itself is of complexity
    ``O(len(l))``. Then, each call to the rank function consists of an
    essentially constant time dictionary lookup.

    TESTS::

        sage: TestSuite(r).run()
    """
    return CallableDict((x,i) for i,x in enumerate(l))


def unrank_from_list(l):
    """
    Return an unrank function from a list.

    EXAMPLES::

        sage: import sage.combinat.ranker as ranker
        sage: l = [1,2,3]
        sage: u = ranker.unrank_from_list(l)
        sage: u(2)
        3
        sage: u(0)
        1
    """
    unrank = lambda j: l[j]
    return unrank


def on_fly():
    """
    Return a pair of enumeration functions rank / unrank.

    rank assigns on the fly an integer, starting from 0, to any object
    passed as argument. The object should be hashable. unrank is the
    inverse function; it returns None for indices that have not yet
    been assigned.

    EXAMPLES::

        sage: [rank, unrank] = sage.combinat.ranker.on_fly()
        sage: rank('a')
        0
        sage: rank('b')
        1
        sage: rank('c')
        2
        sage: rank('a')
        0
        sage: unrank(2)
        'c'
        sage: unrank(3)
        sage: rank('d')
        3
        sage: unrank(3)
        'd'

    .. TODO:: add tests as in combinat::rankers
    """
    def count():
        i = 0
        while True:
            yield i
            i += 1

    counter = count()

    @cached_function
    def rank(x):
        i = next(counter)
        unrank.set_cache(x, i)
        return i

    @cached_function
    def unrank(i):
        return None

    return [rank, unrank]


def unrank(L, i):
    r"""
    Return the `i`-th element of `L`.

    INPUT:

    - ``L`` -- list, tuple, finite enumerated set, etc.
    - ``i`` -- integer

    The purpose of this utility is to give a uniform idiom to recover
    the `i`-th element of an object ``L``, whether ``L`` is a list,
    tuple (or more generally a :class:`collections.abc.Sequence`), an
    enumerated set, some old parent of Sage still implementing
    unranking in the method ``__getitem__``, or an iterable (see
    :class:`collections.abc.Iterable`). See :issue:`15919`.

    EXAMPLES:

    Lists, tuples, and other :class:`sequences <collections.abc.Sequence>`::

        sage: from sage.combinat.ranker import unrank
        sage: unrank(['a','b','c'], 2)
        'c'
        sage: unrank(('a','b','c'), 1)
        'b'
        sage: unrank(range(3,13,2), 1)
        5

    Enumerated sets::

        sage: unrank(GF(7), 2)
        2
        sage: unrank(IntegerModRing(29), 10)
        10

    An iterable::

        sage: unrank(NN,4)
        4

    An iterator::

        sage: unrank(('a{}'.format(i) for i in range(20)), 0)
        'a0'
        sage: unrank(('a{}'.format(i) for i in range(20)), 2)
        'a2'

    .. WARNING::

        When unranking an iterator, it returns the ``i``-th element
        beyond where it is currently at::

            sage: from sage.combinat.ranker import unrank
            sage: it = iter(range(20))
            sage: unrank(it, 2)
            2
            sage: unrank(it, 2)
            5

    TESTS::

        sage: from sage.combinat.ranker import unrank
        sage: unrank(list(range(3)), 10)
        Traceback (most recent call last):
        ...
        IndexError: list index out of range

        sage: unrank(('a{}'.format(i) for i in range(20)), 22)
        Traceback (most recent call last):
        ...
        IndexError: index out of range
    """
    if L in EnumeratedSets():
        return L.unrank(i)
    if isinstance(L, Sequence):
        return L[i]
    if isinstance(L, Parent):
        # handle parents still implementing unranking in __getitem__
        try:
            return L[i]
        except (AttributeError, TypeError, ValueError):
            pass
    if isinstance(L, Iterable):
        try:
            it = iter(L)
            for _ in range(i):
                next(it)
            return next(it)
        except StopIteration:
            raise IndexError("index out of range")
    raise ValueError("do not know how to unrank on {}".format(L))
