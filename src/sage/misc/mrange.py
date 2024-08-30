"""
Multidimensional enumeration

AUTHORS:

- Joel B. Mohler (2006-10-12)

- William Stein (2006-07-19)

- Jon Hanke
"""

# *****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************
from sage.misc.misc_c import prod


def _len(L):
    """
    Determines the length of ``L``.

    Uses either ``cardinality`` or ``__len__`` as appropriate.

    EXAMPLES::

        sage: from sage.misc.mrange import _len
        sage: _len(ZZ)
        +Infinity
        sage: _len(range(4))
        4
        sage: _len(4)
        Traceback (most recent call last):
        ...
        TypeError: object of type 'sage.rings.integer.Integer' has no len()
    """
    try:
        return L.cardinality()
    except AttributeError:
        return len(L)


def _is_finite(L, fallback=True):
    """
    Determines whether ``L`` is finite.

    If ``L`` implements none of ``is_finite``, ``cardinality`` or
    ``__len__``, we assume by default that it is finite for speed
    reasons.  This choice can be overridden by passing in
    ``fallback``.

    EXAMPLES::

        sage: from sage.misc.mrange import _is_finite
        sage: _is_finite(ZZ)
        False
        sage: _is_finite(range(4))
        True
        sage: _is_finite([])
        True
        sage: _is_finite(range(10^8))
        True
        sage: from itertools import product
        sage: _is_finite(product([1],[1])) # does not provide is_finite() or __len__()
        True
    """
    try:
        return L.is_finite()
    except AttributeError:
        pass
    except ValueError:
        # L doesn't know how to determine whether it's finite
        return fallback

    try:
        n = _len(L)
    except (TypeError, AttributeError, NotImplementedError):
        # We usually assume L is finite for speed reasons
        return fallback

    from sage.rings.infinity import infinity
    if n is infinity:
        return False
    return True


def _xmrange_iter(iter_list, typ=list):
    """
    This implements the logic for :func:`mrange_iter` and :class:`xmrange_iter`.

    Note that with typ==list, we will be returning a new copy each
    iteration. This makes it OK to modified the returned list. This
    functionality is relied on in the polynomial iterators. Here's a
    doc-test to prove this::

        sage: iter = sage.misc.mrange._xmrange_iter( [[1,2],[1,3]] )
        sage: l1 = next(iter)
        sage: l2 = next(iter)
        sage: l1 is l2
        False

    However, if you would like to re-use the list object::

        sage: iter = sage.misc.mrange._xmrange_iter( [[1,2],[1,3]], lambda x: x )
        sage: l1 = next(iter)
        sage: l2 = next(iter)
        sage: l1 is l2  # eeek, this is freaky!
        True

    We check that :issue:`14285` has been resolved::

        sage: iter = sage.misc.mrange._xmrange_iter([ZZ,[]])
        sage: next(iter)
        Traceback (most recent call last):
        ...
        StopIteration

    We check that :issue:`28521` is fixed::

        sage: next(sage.misc.mrange._xmrange_iter([[], [1]]))
        Traceback (most recent call last):
        ...
        StopIteration
        sage: next(sage.misc.mrange._xmrange_iter([[1], []]))
        Traceback (most recent call last):
        ...
        StopIteration
        sage: next(sage.misc.mrange._xmrange_iter([[1], [], [2]]))
        Traceback (most recent call last):
        ...
        StopIteration
    """
    if len(iter_list) == 0:
        yield typ()
        return
    # If any iterator in the list is infinite we need to be more careful
    if any(not _is_finite(L) for L in iter_list):
        for L in iter_list:
            try:
                n = _len(L)
            except TypeError:
                continue
        if n == 0:
            return
    curr_iters = [iter(i) for i in iter_list]
    try:
        curr_elt = [next(i) for i in curr_iters[:-1]]
    except StopIteration:
        return
    curr_elt.append(None)
    place = len(iter_list) - 1
    while True:
        try:
            while True:
                curr_elt[place] = next(curr_iters[place])
                if place < len(iter_list) - 1:
                    place += 1
                    curr_iters[place] = iter(iter_list[place])
                    continue
                else:
                    yield typ(curr_elt)
        except StopIteration:
            place -= 1
            if place == -1:
                return


def mrange_iter(iter_list, typ=list):
    """
    Return the multirange list derived from the given list of iterators.

    This is the list version of :func:`xmrange_iter`. Use :class:`xmrange_iter`
    for the iterator.

    More precisely, return the iterator over all objects of type ``typ`` of
    n-tuples of Python ints with entries between 0 and the integers in
    the sizes list. The iterator is empty if sizes is empty or contains
    any nonpositive integer.

    INPUT:

    - ``iter_list`` -- a finite iterable of finite iterables

    - ``typ`` -- (default: list) a type or class; more
      generally, something that can be called with a list as input

    OUTPUT: list

    EXAMPLES::

        sage: mrange_iter([range(3),[0,2]])
        [[0, 0], [0, 2], [1, 0], [1, 2], [2, 0], [2, 2]]
        sage: mrange_iter([['Monty','Flying'],['Python','Circus']], tuple)
        [('Monty', 'Python'), ('Monty', 'Circus'), ('Flying', 'Python'), ('Flying', 'Circus')]
        sage: mrange_iter([[2,3,5,7],[1,2]], sum)
        [3, 4, 4, 5, 6, 7, 8, 9]

    Examples that illustrate empty multi-ranges::

        sage: mrange_iter([range(5),range(3),range(0)])
        []
        sage: mrange_iter([range(5), range(3), range(-2)])
        []

    This example is not empty, and should not be. See :issue:`6561`.

    ::

        sage: mrange_iter([])
        [[]]

    AUTHORS:

    - Joel B. Mohler
    """
    return list(_xmrange_iter(iter_list, typ))


class xmrange_iter:
    """
    Return the multirange iterate derived from the given iterators and
    type.

    .. NOTE::

       This basically gives you the Cartesian product of sets.

    More precisely, return the iterator over all objects of type typ of
    n-tuples of Python ints with entries between 0 and the integers in
    the sizes list. The iterator is empty if sizes is empty or contains
    any nonpositive integer.

    Use :func:`mrange_iter` for the non-iterator form.

    INPUT:

    - ``iter_list`` -- list of objects usable as iterators (possibly
      lists)

    - ``typ`` -- (default: list) a type or class; more generally,
      something that can be called with a list as input

    OUTPUT: a generator

    EXAMPLES: We create multi-range iterators, print them and also
    iterate through a tuple version. ::

        sage: z = xmrange_iter([list(range(3)),list(range(2))], tuple);z
        xmrange_iter([[0, 1, 2], [0, 1]], <... 'tuple'>)
        sage: for a in z:
        ....:     print(a)
        (0, 0)
        (0, 1)
        (1, 0)
        (1, 1)
        (2, 0)
        (2, 1)

    We illustrate a few more iterations.

    ::

        sage: list(xmrange_iter([range(3),range(2)]))
        [[0, 0], [0, 1], [1, 0], [1, 1], [2, 0], [2, 1]]
        sage: list(xmrange_iter([range(3),range(2)], tuple))
        [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)]

    Here we compute the sum of each element of the multi-range
    iterator::

        sage: list(xmrange_iter([range(3),range(2)], sum))
        [0, 1, 1, 2, 2, 3]

    Next we compute the product::

        sage: list(xmrange_iter([range(3),range(2)], prod))
        [0, 0, 0, 1, 0, 2]

    Examples that illustrate empty multi-ranges.

    ::

        sage: list(xmrange_iter([range(5),range(3),range(-2)]))
        []
        sage: list(xmrange_iter([range(5),range(3),range(0)]))
        []

    This example is not empty, and should not be. See :issue:`6561`.

    ::

        sage: list(xmrange_iter([]))
        [[]]

    We use a multi-range iterator to iterate through the Cartesian
    product of sets.

    ::

        sage: X = ['red', 'apple', 389]
        sage: Y = ['orange', 'horse']
        sage: for i,j in xmrange_iter([X, Y], tuple):
        ....:     print((i, j))
        ('red', 'orange')
        ('red', 'horse')
        ('apple', 'orange')
        ('apple', 'horse')
        (389, 'orange')
        (389, 'horse')

    AUTHORS:

    - Joel B. Mohler
    """
    def __init__(self, iter_list, typ=list):
        self.iter_list = iter_list
        self.typ = typ

    def __repr__(self):
        if self.typ == list:
            return 'xmrange_iter(%s)' % self.iter_list
        else:
            return 'xmrange_iter(%s, %s)' % (self.iter_list, self.typ)

    def __iter__(self):
        return _xmrange_iter(self.iter_list, self.typ)

    def __len__(self):
        """
        Return the cardinality of this iterator as an int.

        This raises a :exc:`TypeError` if the cardinality does not fit
        into a Python int.

        EXAMPLES::

            sage: C = cartesian_product_iterator([range(3),range(4)])
            sage: len(C)
            12
            sage: len(cartesian_product_iterator([]))
            1
            sage: len(cartesian_product_iterator([ZZ,[]]))
            0
            sage: len(cartesian_product_iterator([ZZ,ZZ]))
            Traceback (most recent call last):
            ...
            TypeError: cardinality does not fit into a Python int
        """
        n = self.cardinality()
        try:
            n = int(n)
        except TypeError:
            raise TypeError("cardinality does not fit into a Python int")
        return n

    def cardinality(self):
        """
        Return the cardinality of this iterator.

        EXAMPLES::

            sage: C = cartesian_product_iterator([range(3),range(4)])
            sage: C.cardinality()
            12
            sage: C = cartesian_product_iterator([ZZ,QQ])
            sage: C.cardinality()
            +Infinity
            sage: C = cartesian_product_iterator([ZZ,[]])
            sage: C.cardinality()
            0
        """
        from sage.rings.integer import Integer
        from sage.rings.infinity import infinity
        ans = Integer(1)
        found_infinity = False
        for L in self.iter_list:
            try:
                n = L.cardinality()
            except AttributeError:
                n = Integer(len(L))
            if n == 0:
                return Integer(0)
            elif n is infinity:
                found_infinity = True
            elif not found_infinity:
                ans *= n
        if found_infinity:
            return infinity
        else:
            return ans


def _xmrange(sizes, typ=list):
    n = len(sizes)
    if n == 0:
        yield typ([])
        return
    for i in sizes:
        if i <= 0:
            return
    v = [0] * n    # make a list of n 0's.
    v[-1] = -1
    ptr_max = n - 1
    ptr = ptr_max
    while True:
        while True:
            if ptr != -1 and v[ptr] + 1 < sizes[ptr]:
                v[ptr] += 1
                ptr = ptr_max
                break
            elif ptr != -1:
                v[ptr] = 0
                ptr -= 1
            else:
                return
        yield typ(v)   # make a copy of v!


def mrange(sizes, typ=list):
    """
    Return the multirange list with given sizes and type.

    This is the list version of :class:`xmrange`.
    Use :class:`xmrange` for the iterator.

    More precisely, return the iterator over all objects of type typ of
    n-tuples of Python ints with entries between 0 and the integers in
    the sizes list. The iterator is empty if sizes is empty or contains
    any nonpositive integer.

    INPUT:

    - ``sizes`` -- list of nonnegative integers

    - ``typ`` -- (default: list) a type or class; more
      generally, something that can be called with a list as input

    OUTPUT: list

    EXAMPLES::

        sage: mrange([3,2])
        [[0, 0], [0, 1], [1, 0], [1, 1], [2, 0], [2, 1]]
        sage: mrange([3,2], tuple)
        [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)]
        sage: mrange([3,2], sum)
        [0, 1, 1, 2, 2, 3]

    Examples that illustrate empty multi-ranges::

        sage: mrange([5,3,-2])
        []
        sage: mrange([5,3,0])
        []

    This example is not empty, and should not be. See :issue:`6561`.

    ::

        sage: mrange([])
        [[]]

    AUTHORS:

    - Jon Hanke

    - William Stein
    """
    return list(_xmrange(sizes, typ))


class xmrange:
    """
    Return the multirange iterate with given sizes and type.

    More precisely, return the iterator over all objects of type typ of
    n-tuples of Python ints with entries between 0 and the integers in
    the sizes list. The iterator is empty if sizes is empty or contains
    any nonpositive integer.

    Use mrange for the non-iterator form.

    INPUT:

    - ``sizes`` -- list of nonnegative integers

    - ``typ`` -- (default: list) a type or class; more
      generally, something that can be called with a list as input

    OUTPUT: a generator

    EXAMPLES: We create multi-range iterators, print them and also
    iterate through a tuple version.

    ::

        sage: z = xmrange([3,2]);z
        xmrange([3, 2])
        sage: z = xmrange([3,2], tuple);z
        xmrange([3, 2], <... 'tuple'>)
        sage: for a in z:
        ....:     print(a)
        (0, 0)
        (0, 1)
        (1, 0)
        (1, 1)
        (2, 0)
        (2, 1)

    We illustrate a few more iterations.

    ::

        sage: list(xmrange([3,2]))
        [[0, 0], [0, 1], [1, 0], [1, 1], [2, 0], [2, 1]]
        sage: list(xmrange([3,2], tuple))
        [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)]

    Here we compute the sum of each element of the multi-range
    iterator::

        sage: list(xmrange([3,2], sum))
        [0, 1, 1, 2, 2, 3]

    Next we compute the product::

        sage: list(xmrange([3,2], prod))
        [0, 0, 0, 1, 0, 2]

    Examples that illustrate empty multi-ranges.

    ::

        sage: list(xmrange([5,3,-2]))
        []
        sage: list(xmrange([5,3,0]))
        []

    This example is not empty, and should not be. See :issue:`6561`.

    ::

        sage: list(xmrange([]))
        [[]]

    We use a multi-range iterator to iterate through the Cartesian
    product of sets.

    ::

        sage: X = ['red', 'apple', 389]
        sage: Y = ['orange', 'horse']
        sage: for i,j in xmrange([len(X), len(Y)]):
        ....:     print((X[i], Y[j]))
        ('red', 'orange')
        ('red', 'horse')
        ('apple', 'orange')
        ('apple', 'horse')
        (389, 'orange')
        (389, 'horse')

    AUTHORS:

    - Jon Hanke

    - William Stein
    """
    def __init__(self, sizes, typ=list):
        self.sizes = [int(x) for x in sizes]
        self.typ = typ

    def __repr__(self):
        if self.typ == list:
            return 'xmrange(%s)' % self.sizes
        else:
            return 'xmrange(%s, %s)' % (self.sizes, self.typ)

    def __len__(self):
        sizes = self.sizes
        n = len(sizes)
        if n == 0:
            return 0
        for i in sizes:
            if i <= 0:
                return 0
        return prod(sizes, 1)

    def __iter__(self):
        return _xmrange(self.sizes, self.typ)


def cartesian_product_iterator(X):
    """
    Iterate over the Cartesian product.

    INPUT:

    - ``X`` -- list or tuple of lists

    OUTPUT: iterator over the Cartesian product of the elements of X

    EXAMPLES::

        sage: list(cartesian_product_iterator([[1,2], ['a','b']]))
        [(1, 'a'), (1, 'b'), (2, 'a'), (2, 'b')]
        sage: list(cartesian_product_iterator([]))
        [()]

    TESTS:

    Check that :issue:`28521` is fixed::

        sage: list(cartesian_product_iterator([[], [1]]))
        []
        sage: list(cartesian_product_iterator([[1], []]))
        []
        sage: list(cartesian_product_iterator([[1], [], [2]]))
        []
    """
    return xmrange_iter(X, tuple)


def cantor_product(*args, **kwds):
    r"""
    Return an iterator over the product of the inputs along the diagonals a la
    :wikipedia:`Cantor pairing <Pairing_function#Cantor_pairing_function>`.

    INPUT:

    - a certain number of iterables

    - ``repeat`` -- an optional integer. If it is provided, the input is
      repeated ``repeat`` times

    Other keyword arguments are passed to
    :class:`sage.combinat.integer_lists.invlex.IntegerListsLex`.

    EXAMPLES::

        sage: from sage.misc.mrange import cantor_product
        sage: list(cantor_product([0, 1], repeat=3))
        [(0, 0, 0),
         (1, 0, 0),
         (0, 1, 0),
         (0, 0, 1),
         (1, 1, 0),
         (1, 0, 1),
         (0, 1, 1),
         (1, 1, 1)]
        sage: list(cantor_product([0, 1], [0, 1, 2, 3]))
        [(0, 0), (1, 0), (0, 1), (1, 1), (0, 2), (1, 2), (0, 3), (1, 3)]

    Infinite iterators are valid input as well::

       sage: from itertools import islice
       sage: list(islice(cantor_product(ZZ, QQ), 14r))
        [(0, 0),
         (1, 0),
         (0, 1),
         (-1, 0),
         (1, 1),
         (0, -1),
         (2, 0),
         (-1, 1),
         (1, -1),
         (0, 1/2),
         (-2, 0),
         (2, 1),
         (-1, -1),
         (1, 1/2)]

    TESTS::

        sage: C = cantor_product([0, 1], [0, 1, 2, 3], [0, 1, 2])
        sage: sum(1 for _ in C) == 2*4*3
        True

        sage: from itertools import count
        sage: list(cantor_product([], count()))
        []
        sage: list(cantor_product(count(), [], count()))
        []

        sage: list(cantor_product(count(), repeat=0))
        [()]

        sage: next(cantor_product(count(), repeat=-1))
        Traceback (most recent call last):
        ...
        ValueError: repeat argument cannot be negative
        sage: next(cantor_product(count(), toto='hey'))
        Traceback (most recent call last):
        ...
        TypeError: ...__init__() got an unexpected keyword argument 'toto'

    ::

        sage: list(cantor_product(srange(5), repeat=2, min_slope=1))
        [(0, 1), (0, 2), (1, 2), (0, 3), (1, 3),
         (0, 4), (2, 3), (1, 4), (2, 4), (3, 4)]

    Check that :issue:`24897` is fixed::

        sage: from sage.misc.mrange import cantor_product
        sage: list(cantor_product([1]))
        [(1,)]
        sage: list(cantor_product([1], repeat=2))
        [(1, 1)]
        sage: list(cantor_product([1], [1,2]))
        [(1, 1), (1, 2)]
        sage: list(cantor_product([1,2], [1]))
        [(1, 1), (2, 1)]
    """
    from itertools import count
    from sage.combinat.integer_lists import IntegerListsLex

    m = len(args)                         # numer of factors
    lengths = [None] * m                  # None or length of factors
    data = [[] for _ in range(m)]         # the initial slice of each factor
    iterators = [iter(a) for a in args]   # the iterators
    repeat = int(kwds.pop('repeat', 1))
    if repeat == 0:
        yield ()
        return
    elif repeat < 0:
        raise ValueError("repeat argument cannot be negative")
    mm = m * repeat

    for n in count(0):
        # try to add one more term to each bin
        for i, a in enumerate(iterators):
            if lengths[i] is None:
                try:
                    data[i].append(next(a))
                except StopIteration:
                    assert len(data[i]) == n
                    if n == 0:
                        return
                    lengths[i] = n

        # iterate through what we have
        ceiling = [n if lengths[i] is None else lengths[i] - 1
                   for i in range(m)] * repeat
        for v in IntegerListsLex(n, length=mm, ceiling=ceiling, **kwds):
            yield tuple(data[i % m][v[i]] for i in range(mm))

        if all(l is not None for l in lengths) and repeat * sum(l - 1 for l in lengths) <= n:
            return
