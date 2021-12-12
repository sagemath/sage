.. _tutorial-comprehensions:

==================================================
Tutorial: Comprehensions, Iterators, and Iterables
==================================================

.. MODULEAUTHOR:: Florent Hivert <florent.hivert@univ-rouen.fr> and
   Nicolas M. Thiéry <nthiery at users.sf.net>

.. linkall

List comprehensions
===================

*List comprehensions* are a very handy way to construct lists in
Python. You can use either of the following idioms:

.. CODE-BLOCK:: python

      [ <expr> for <name> in <iterable> ]
      [ <expr> for <name> in <iterable> if <condition> ]

For example, here are some lists of squares::

    sage: [ i^2 for i in [1, 3, 7] ]
    [1, 9, 49]
    sage: [ i^2 for i in range(1,10) ]
    [1, 4, 9, 16, 25, 36, 49, 64, 81]
    sage: [ i^2 for i in range(1,10) if i % 2 == 1]
    [1, 9, 25, 49, 81]

And a variant on the latter::

    sage: [i^2 if i % 2 == 1 else 2 for i in range(10)]
    [2, 1, 2, 9, 2, 25, 2, 49, 2, 81]

.. TOPIC:: Exercises

    #.  Construct the list of the squares of the prime integers between 1 and 10::

            sage: # edit here

    #.  Construct the list of the perfect squares less than 100 (hint: use
        :func:`srange` to get a list of Sage integers together with the
        method ``i.sqrtrem()``)::

            sage: # edit here

One can use more than one iterable in a list comprehension::

    sage: [ (i,j) for i in range(1,6) for j in range(1,i) ]
    [(2, 1), (3, 1), (3, 2), (4, 1), (4, 2), (4, 3), (5, 1), (5, 2), (5, 3), (5, 4)]

.. warning:: Mind the order of the nested loop in the previous expression.

If instead one wants to build a list of lists, one can use nested lists as in::

    sage: [ [ binomial(n, i) for i in range(n+1) ] for n in range(10) ]
    [[1],
    [1, 1],
    [1, 2, 1],
    [1, 3, 3, 1],
    [1, 4, 6, 4, 1],
    [1, 5, 10, 10, 5, 1],
    [1, 6, 15, 20, 15, 6, 1],
    [1, 7, 21, 35, 35, 21, 7, 1],
    [1, 8, 28, 56, 70, 56, 28, 8, 1],
    [1, 9, 36, 84, 126, 126, 84, 36, 9, 1]]

.. TOPIC:: Exercises

    #.  Compute the list of pairs `(i,j)` of non negative integers such
        that ``i`` is at most `5`, ``j`` is at most ``8``, and ``i`` and
        ``j`` are co-prime::

            sage: # edit here

    #.  Compute the same list for `i < j < 10`::

            sage: # edit here


Iterators
=========

Definition
----------

To build a comprehension, Python actually uses an *iterator*. This is
a device which runs through a bunch of objects, returning one at each
call to the ``next`` method. Iterators are built using parentheses::

    sage: it = (binomial(8, i) for i in range(9))
    sage: next(it)
    1

::

    sage: next(it)
    8
    sage: next(it)
    28
    sage: next(it)
    56

You can get the list of the results that are not yet *consumed*::

    sage: list(it)
    [70, 56, 28, 8, 1]

Asking for more elements triggers a ``StopIteration`` exception::

    sage: next(it)
    Traceback (most recent call last):
    ...
    StopIteration

An iterator can be used as argument for a function. The two following
idioms give the same results; however, the second idiom is much more
memory efficient (for large examples) as it does not expand any list
in memory::

    sage: sum([binomial(8, i) for i in range(9)])
    256
    sage: sum(binomial(8, i) for i in range(9))
    256

.. TOPIC:: Exercises

    #.  Compute the sum of `\binom{10}{i}` for all even `i`::

            sage: # edit here

    #.  Compute the sum of the products of all pairs of co-prime numbers `i, j` for `i<j<10`::

            sage: # edit here


Typical usage of iterators
--------------------------

Iterators are very handy with the functions :func:`all`, :func:`any`,
and :func:`exists`::

    sage: all([True, True, True, True])
    True
    sage: all([True, False, True, True])
    False

::

    sage: any([False, False, False, False])
    False
    sage: any([False, False, True, False])
    True

Let's check that all the prime numbers larger than 2 are odd::

    sage: all( is_odd(p) for p in range(1,100) if is_prime(p) and p>2 )
    True

It is well know that if ``2^p-1`` is prime then ``p`` is prime::

    sage: def mersenne(p): return 2^p -1
    sage: [ is_prime(p) for p in range(20) if is_prime(mersenne(p)) ]
    [True, True, True, True, True, True, True]

The converse is not true::

    sage: all( is_prime(mersenne(p)) for p in range(1000) if is_prime(p) )
    False

Using a list would be much slower here::

    sage: %time all( is_prime(mersenne(p)) for p in range(1000) if is_prime(p) )    # not tested
    CPU times: user 0.00 s, sys: 0.00 s, total: 0.00 s
    Wall time: 0.00 s
    False
    sage: %time all( [ is_prime(mersenne(p)) for p in range(1000) if is_prime(p)] ) # not tested
    CPU times: user 0.72 s, sys: 0.00 s, total: 0.73 s
    Wall time: 0.73 s
    False

You can get the counterexample using :func:`exists`. It takes two
arguments: an iterator and a function which tests the property that
should hold::

    sage: exists( (p for p in range(1000) if is_prime(p)), lambda p: not is_prime(mersenne(p)) )
    (True, 11)

An alternative way to achieve this is::

    sage: counter_examples = (p for p in range(1000) if is_prime(p) and not is_prime(mersenne(p)))
    sage: next(counter_examples)
    11

.. TOPIC:: Exercises

    #.  Build the list `\{i^3 \mid -10<i<10\}`. Can you find two of those
        cubes `u` and `v` such that `u + v = 218`?

        ::

           sage: # edit here

itertools
---------

At its name suggests :mod:`itertools` is a module which defines
several handy tools for manipulating iterators::

    sage: l = [3, 234, 12, 53, 23]
    sage: [(i, l[i]) for i in range(len(l))]
    [(0, 3), (1, 234), (2, 12), (3, 53), (4, 23)]

The same results can be obtained using :func:`enumerate`::

    sage: list(enumerate(l))
    [(0, 3), (1, 234), (2, 12), (3, 53), (4, 23)]

Here is the analogue of list slicing::

    sage: list(Permutations(3))
    [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]
    sage: list(Permutations(3))[1:4]
    [[1, 3, 2], [2, 1, 3], [2, 3, 1]]

    sage: import itertools
    sage: list(itertools.islice(Permutations(3), 1r, 4r))
    [[1, 3, 2], [2, 1, 3], [2, 3, 1]]

Note that all calls to ``islice`` must have arguments of type ``int`` and
not Sage integers.

The behaviour of the functions :func:`map` and :func:`filter` has
changed between Python 2 and Python 3. In Python 3, they return an
iterator. If you want to return a list like in Python 2 you need to explicitly
wrap them in :func:`list`::

    sage: list(map(lambda z: z.cycle_type(), Permutations(3)))
    [[1, 1, 1], [2, 1], [2, 1], [3], [3], [2, 1]]

    sage: list(filter(lambda z: z.has_pattern([1,2]), Permutations(3)))
    [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2]]

.. TOPIC:: Exercises

    #.  Define an iterator for the `i`-th prime for `5<i<10`::

           sage: # edit here

Defining new iterators
----------------------

One can very easily write new iterators using the keyword
``yield``. The following function does nothing interesting beyond
demonstrating the use of ``yield``::

    sage: def f(n):
    ....:   for i in range(n):
    ....:       yield i
    sage: [ u for u in f(5) ]
    [0, 1, 2, 3, 4]

Iterators can be recursive::

    sage: def words(alphabet,l):
    ....:    if l == 0:
    ....:        yield []
    ....:    else:
    ....:        for word in words(alphabet, l-1):
    ....:            for a in alphabet:
    ....:                yield word + [a]

    sage: [ w for w in words(['a','b','c'], 3) ]
    [['a', 'a', 'a'], ['a', 'a', 'b'], ['a', 'a', 'c'], ['a', 'b', 'a'], ['a', 'b', 'b'], ['a', 'b', 'c'], ['a', 'c', 'a'], ['a', 'c', 'b'], ['a', 'c', 'c'], ['b', 'a', 'a'], ['b', 'a', 'b'], ['b', 'a', 'c'], ['b', 'b', 'a'], ['b', 'b', 'b'], ['b', 'b', 'c'], ['b', 'c', 'a'], ['b', 'c', 'b'], ['b', 'c', 'c'], ['c', 'a', 'a'], ['c', 'a', 'b'], ['c', 'a', 'c'], ['c', 'b', 'a'], ['c', 'b', 'b'], ['c', 'b', 'c'], ['c', 'c', 'a'], ['c', 'c', 'b'], ['c', 'c', 'c']]
    sage: sum(1 for w in words(['a','b','c'], 3))
    27

Here is another recursive iterator::

    sage: def dyck_words(l):
    ....:     if l==0:
    ....:         yield ''
    ....:     else:
    ....:         for k in range(l):
    ....:             for w1 in dyck_words(k):
    ....:                 for w2 in dyck_words(l-k-1):
    ....:                     yield '('+w1+')'+w2

    sage: list(dyck_words(4))
    ['()()()()',
    '()()(())',
    '()(())()',
    '()(()())',
    '()((()))',
    '(())()()',
    '(())(())',
    '(()())()',
    '((()))()',
    '(()()())',
    '(()(()))',
    '((())())',
    '((()()))',
    '(((())))']

    sage: sum(1 for w in dyck_words(5))
    42

.. TOPIC:: Exercises

    #.  Write an iterator with two parameters `n`, `l` iterating
        through the set of nondecreasing lists of integers smaller
        than `n` of length `l`::

           sage: # edit here


Standard Iterables
==================

Finally, many standard Python and Sage objects are *iterable*; that is
one may iterate through their elements::

    sage: sum( x^len(s) for s in Subsets(8) )
    x^8 + 8*x^7 + 28*x^6 + 56*x^5 + 70*x^4 + 56*x^3 + 28*x^2 + 8*x + 1

    sage: sum( x^p.length() for p in Permutations(3) )
    x^3 + 2*x^2 + 2*x + 1

    sage: factor(sum( x^p.length() for p in Permutations(3) ))
    (x^2 + x + 1)*(x + 1)

    sage: P = Permutations(5)
    sage: all( p in P for p in P )
    True

    sage: for p in GL(2, 2): print(p); print("")
    [1 0]
    [0 1]
    <BLANKLINE>
    [0 1]
    [1 0]
    <BLANKLINE>
    [0 1]
    [1 1]
    <BLANKLINE>
    [1 1]
    [0 1]
    <BLANKLINE>
    [1 1]
    [1 0]
    <BLANKLINE>
    [1 0]
    [1 1]
    <BLANKLINE>

    sage: for p in Partitions(3): print(p)
    [3]
    [2, 1]
    [1, 1, 1]


.. skip

Beware of infinite loops::

    sage: for p in Partitions(): print(p)          # not tested

.. skip

::

    sage: for p in Primes(): print(p)              # not tested

Infinite loops can nevertheless be very useful::

    sage: exists( Primes(), lambda p: not is_prime(mersenne(p)) )
    (True, 11)


    sage: counter_examples = (p for p in Primes() if not is_prime(mersenne(p)))
    sage: next(counter_examples)
    11
