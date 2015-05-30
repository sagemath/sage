"""
Iterators

Miscellaneous functions which should eventually be moved upstream into
Python's standard itertools module.

TODO: Cython this file as soon as Cython supports ``yield``.
"""
#*****************************************************************************
#  Copyright (C) 2010     Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

import itertools, heapq

def unique_merge(*lists):
    """
    INPUT:

     - ``lists``: sorted lists (or iterables)

    Return an iterator over the elements of each list in ``lists``, in
    sorted order, with duplicates removed.

        sage: from sage.misc.sage_itertools import unique_merge
        sage: list(unique_merge([1,2,2,3,4,7,9], [0,2,4], [2,5]))
        [0, 1, 2, 3, 4, 5, 7, 9]

    Inspired from: http://rosettacode.org/wiki/Create_a_Sequence_of_unique_elements#Python
    """
    return (k for k,g in itertools.groupby(heapq.merge(*lists)))

def min_cmp(L, cmp=None):
    """
    Returns the smallest item of a list (or iterable) with respect to
    a comparison function.

    INPUT:
        ``L``   -- an iterable
        ``cmp`` -- an optional comparison function.

    ``cmp(x, y)`` should return a negative value if `x < y`, `0` if
    `x == y`, and a positive value if `x > y`.

    OUTPUT: the smallest item of ``L`` with respect to ``cmp``.

    EXAMPLES::

        sage: from sage.misc.sage_itertools import min_cmp
        sage: L = [1,-1,3,-1,3,2]
        sage: min_cmp(L)
        -1
        sage: def mycmp(x,y): return y - x
        sage: min_cmp(L, mycmp)
        3

    The input can be any iterable::

        sage: min_cmp( (x^2 for x in L) )
        1
        sage: min_cmp( (x^2 for x in L), mycmp)
        9

    Computing the min of an empty iterable raises and error::

        sage: min_cmp([])
        Traceback (most recent call last):
        ...
        ValueError: min() arg is an empty sequence
        sage: min_cmp([], mycmp)
        Traceback (most recent call last):
        ...
        ValueError: min_cmp() arg is an empty sequence
    """
    if cmp is None:
        return min(L) # Resort to Python's standard min

    iterator = iter(L)
    try:
        m = next(iterator)
    except StopIteration:
        raise ValueError("min_cmp() arg is an empty sequence")
    for item in iterator:
        if cmp(item, m) < 0:
            m = item
    return m

def max_cmp(L, cmp=None):
    """
    Returns the largest item of a list (or iterable) with respect to a
    comparison function.

    INPUT:
        ``L``   -- an iterable
        ``cmp`` -- an optional comparison function.

    ``cmp(x, y)`` should return a negative value if `x < y`, `0` if
    `x == y`, and a positive value if `x > y`.

    OUTPUT: the largest item of ``L`` with respect to ``cmp``.

    EXAMPLES::

        sage: from sage.misc.sage_itertools import max_cmp
        sage: L = [1,-1,3,-1,3,2]
        sage: max_cmp(L)
        3
        sage: def mycmp(x,y): return y - x
        sage: max_cmp(L, mycmp)
        -1

    The input can be any iterable::

        sage: max_cmp( (x^2 for x in L) )
        9
        sage: max_cmp( (x^2 for x in L), mycmp)
        1

    Computing the max of an empty iterable raises and error::

        sage: max_cmp([])
        Traceback (most recent call last):
        ...
        ValueError: max() arg is an empty sequence
        sage: max_cmp([], mycmp)
        Traceback (most recent call last):
        ...
        ValueError: max_cmp() arg is an empty sequence
    """
    if cmp is None:
        return max(L) # Resort to Python's standard max

    iterator = iter(L)
    try:
        m = next(iterator)
    except StopIteration:
        raise ValueError("max_cmp() arg is an empty sequence")
    for item in iterator:
        if cmp(item, m) > 0:
            m = item
    return m

def imap_and_filter_none(function, iterable):
    r"""
    Returns an iterator over the elements ``function(x)``, where ``x``
    iterates through ``iterable``, such that ``function(x)`` is not ``None``.

    EXAMPLES::

        sage: from sage.misc.sage_itertools import imap_and_filter_none
        sage: p = imap_and_filter_none(lambda x: x if is_prime(x) else None, range(15))
        sage: [next(p), next(p), next(p), next(p), next(p), next(p)]
        [2, 3, 5, 7, 11, 13]
        sage: p = imap_and_filter_none(lambda x: x+x, ['a','b','c','d','e'])
        sage: [next(p), next(p), next(p), next(p), next(p)]
        ['aa', 'bb', 'cc', 'dd', 'ee']
    """
    for x in iterable:
        x = function(x)
        if x is not None:
            yield x
