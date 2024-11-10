r"""
Lazy lists

A lazy list is an iterator that behaves like a list and possesses a cache
mechanism. A lazy list is potentially infinite and speed performances of the
cache is comparable with Python lists. One major difference with original
Python list is that lazy list are immutable. The advantage is that slices
share memory.

EXAMPLES::

    sage: from sage.misc.lazy_list import lazy_list
    sage: P = lazy_list(Primes())
    sage: P[100]                                                                        # needs sage.libs.pari
    547
    sage: P[10:34]                                                                      # needs sage.libs.pari
    lazy list [31, 37, 41, ...]
    sage: P[12:23].list()                                                               # needs sage.libs.pari
    [41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83]

    sage: f = lazy_list((i**2 - 3*i for i in range(10)))
    sage: print(" ".join(str(i) for i in f))
    0 -2 -2 0 4 10 18 28 40 54
    sage: i1 = iter(f)
    sage: i2 = iter(f)
    sage: [next(i1), next(i1)]
    [0, -2]
    sage: [next(i2), next(i2)]
    [0, -2]
    sage: [next(i1), next(i2)]
    [-2, -2]

It is possible to prepend a list to a lazy list::

    sage: from itertools import count
    sage: l = [3,7] + lazy_list(i**2 for i in count())
    sage: l
    lazy list [3, 7, 0, ...]

But, naturally, not the other way around::

    sage: lazy_list(i-1 for i in count()) + [3,2,5]
    Traceback (most recent call last):
    ...
    TypeError: can only add list to lazy_list

You can easily create your own class inheriting from :class:`lazy_list_generic`. You
should call the :class:`lazy_list_generic` constructor (optionally with some
precomputed values for the cache) and implement the method ``_new_slice`` that
returns a new chunk of data at each call. Here is an example of implementation
of the Thue--Morse word that is obtained as the fixed point of the substitution
`0 \to 01` and `1 \to 10`::

    sage: from sage.misc.lazy_list import lazy_list_generic
    sage: class MyThueMorseWord(lazy_list_generic):
    ....:     def __init__(self):
    ....:         self.i = 1
    ....:         lazy_list_generic.__init__(self, cache=[0,1])
    ....:     def _new_slice(self):
    ....:         letter = self.get(self.i)
    ....:         self.i += 1
    ....:         return [0,1] if letter == 0 else [1,0]
    sage: w = MyThueMorseWord()
    sage: w
    lazy list [0, 1, 1, ...]
    sage: all(w[i] == ZZ(i).popcount()%2 for i in range(100))
    True
    sage: w[:500].list() == w[:1000:2].list()
    True

Alternatively, you can create the lazy list from an update function::

    sage: def thue_morse_update(values):
    ....:     n = len(values)
    ....:     if n == 0:
    ....:         letter = 0
    ....:     else:
    ....:         assert n%2 == 0
    ....:         letter = values[n//2]
    ....:     values.append(letter)
    ....:     values.append(1-letter)
    sage: w2 = lazy_list(update_function=thue_morse_update)
    sage: w2
    lazy list [0, 1, 1, ...]
    sage: w2[:500].list() == w[:500].list()
    True

You can also create user-defined classes (Python) and extension types (Cython)
inheriting from :class:`lazy_list_generic`.
In that case you would better implement directly the method
``_update_cache_up_to``. See the examples in this file with the classes
:class:`lazy_list_from_iterator` and :class:`lazy_list_from_function`.

Classes and Methods
===================
"""
# ****************************************************************************
#       Copyright (C) 2015 Vincent Delecroix <20100.delecroix@gmail.com>
#       Copyright (C) 2016 Daniel Krenn <dev@danielkrenn.at>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

cdef extern from "Python.h":
    Py_ssize_t PY_SSIZE_T_MAX

# make a unique instance of empty lazy lists
cdef lazy_list_generic empty_lazy_list
empty_lazy_list = lazy_list_generic.__new__(lazy_list_generic)
empty_lazy_list.start = 0
empty_lazy_list.stop = 0
empty_lazy_list.step = 1
empty_lazy_list.cache = []


def lazy_list(data=None, initial_values=None, start=None, stop=None, step=None,
        update_function=None):
    r"""
    Return a lazy list.

    INPUT:

    - ``data`` -- data to create a lazy list from. This can be

      #. a (possibly infinite) iterable,
      #. a function (that takes as input an integer ``n`` and return
         the ``n``-th term of the list),
      #. or a standard Python container ``list`` or ``tuple``.

    - ``initial_values`` -- the beginning of the sequence that will not be computed from
      the ``data`` provided

    - ``update_function`` -- you can also construct a lazy list from a function
      that takes as input a list of precomputed values and updates it with some
      more values

    .. NOTE::

        If you want finer tuning of the constructor you can directly instantiate
        the classes associated to lazy lists that are
        :class:`lazy_list_generic`, :class:`lazy_list_from_iterator`,
        :class:`lazy_list_from_function`.

    EXAMPLES:

    The basic construction of lazy lists.
    ::

        sage: from sage.misc.lazy_list import lazy_list

    #. Iterators::

        sage: from itertools import count
        sage: lazy_list(count())
        lazy list [0, 1, 2, ...]

    #. Functions::

        sage: lazy_list(lambda n: (n**2)%17)
        lazy list [0, 1, 4, ...]

    #. Plain lists::

        sage: lazy_list([1,5,7,2])
        lazy list [1, 5, 7, ...]

    If a function is only defined for large values, you can provide the beginning
    of the sequence manually::

        sage: l = lazy_list(divisors, [None])
        sage: l
        lazy list [None, [1], [1, 2], ...]

    Lazy lists behave like lists except that they are immutable::

        sage: l[3::5]
        lazy list [[1, 3], [1, 2, 4, 8], [1, 13], ...]

    If your lazy list is finite, you can obtain the underlying list with the
    method `.list()`::

        sage: l[30:50:5].list()
        [[1, 2, 3, 5, 6, 10, 15, 30],
         [1, 5, 7, 35],
         [1, 2, 4, 5, 8, 10, 20, 40],
         [1, 3, 5, 9, 15, 45]]

    TESTS::

        sage: lazy_list()
        lazy list []
        sage: lazy_list(data='hey', update_function='hello')
        Traceback (most recent call last):
        ...
        ValueError: only one of the arguments 'data' or 'update_function'
        can be used

    Applying :func:`lazy_list` is idempotent. Thus, lazy lists created
    from each other are unmodified::

        sage: C = lazy_list(count())
        sage: C[4]
        4
        sage: D = lazy_list(C)
        sage: C is D
        True
    """
    cdef lazy_list_generic l

    if data is None and update_function is None:
        return empty_lazy_list
    elif data is not None and update_function is not None:
        raise ValueError("only one of the arguments 'data' or 'update_function' can be used")

    if initial_values is None:
        cache = []
    else:
        cache = list(initial_values)

    if update_function is not None:
        assert callable(update_function)
        return lazy_list_from_update_function(update_function, cache)

    if isinstance(data, (tuple,list)):
        data = cache + list(data)
        l = lazy_list_generic(data, start=0, stop=len(data), step=1)
    elif isinstance(data, lazy_list_generic):
        return data
    else:
        # the code below is not very clean
        # we just want to differentiate on the one hand iterable (= object with a
        # proper __iter__ method)/iterator (= object with a next method) and on the
        # other hand callable (= object with __call__)
        try:
            data = iter(data)
        except TypeError:
            pass

        from sage.misc.misc import is_iterator
        if is_iterator(data):
            l = lazy_list_from_iterator(iter(data), cache)
        elif callable(data):
            l = lazy_list_from_function(data, cache)
        else:
            raise ValueError("not able to build a lazy list from {}".format(type(data)))

    return l


def slice_unpickle(master, start, stop, step):
    r"""
    Unpickle helper.

    TESTS::

        sage: from sage.misc.lazy_list import slice_unpickle
        sage: slice_unpickle(list(range(35)), 1, 3, 7) == list(range(35))[1:3:7]
        True
    """
    return master[start:stop:step]


def lazy_list_formatter(L, name='lazy list',
                        separator=', ', more='...',
                        opening_delimiter='[', closing_delimiter=']',
                        preview=3):
    r"""
    Return a string representation of ``L``.

    INPUT:

    - ``L`` -- an iterable object

    - ``name`` -- (default: ``'lazy list'``) a string appearing
      at first position (i.e., in front of the actual values)
      in the representation

    - ``opening_delimiter`` -- (default: ``'['``) a string heading
      the shown entries

    - ``closing_delimiter`` -- (default: ``']'``) a string trailing
      the shown entries

    - ``separator`` -- (default: ``', '``) a string appearing between
      two entries

    - ``more`` -- (default: ``'...'``) a string indicating that
      not all entries of the list are shown

    - ``preview`` -- (default: ``3``) an integer specifying the number of
      elements shown in the representation string

    OUTPUT: string

    EXAMPLES::

        sage: from sage.misc.lazy_list import lazy_list_formatter
        sage: lazy_list_formatter(srange(3, 1000, 5), name='list')
        'list [3, 8, 13, ...]'

    ::

        sage: # needs sage.libs.pari
        sage: from sage.misc.lazy_list import lazy_list
        sage: L = lazy_list(Primes()); L
        lazy list [2, 3, 5, ...]
        sage: repr(L) == lazy_list_formatter(L)
        True
        sage: lazy_list_formatter(L, name='primes')
        'primes [2, 3, 5, ...]'
        sage: lazy_list_formatter(L, opening_delimiter='(', closing_delimiter=')')
        'lazy list (2, 3, 5, ...)'
        sage: lazy_list_formatter(L, opening_delimiter='', closing_delimiter='')
        'lazy list 2, 3, 5, ...'
        sage: lazy_list_formatter(L, separator='--')
        'lazy list [2--3--5--...]'
        sage: lazy_list_formatter(L, more='and more')
        'lazy list [2, 3, 5, and more]'
        sage: lazy_list_formatter(L, preview=10)
        'lazy list [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, ...]'
        sage: lazy_list_formatter(L, name='primes',
        ....:                     opening_delimiter='', closing_delimiter='',
        ....:                     separator=' ', more='->', preview=7)
        'primes 2 3 5 7 11 13 17 ->'

    TESTS::

        sage: from itertools import count
        sage: lazy_list_formatter(count(), name='iterator count')
        'iterator count [0, 1, 2, ...]'
    """
    from itertools import islice

    cdef str s = name
    if s:
        s += ' '
    s += opening_delimiter
    cdef list P = list(islice(L, int(preview + 1)))
    cdef list E = list(repr(e) for e in P[:preview])
    if len(P) > preview:
        E.append(more)
    s += separator.join(E)
    s += closing_delimiter
    return s


cdef class lazy_list_generic():
    r"""
    A lazy list.

    EXAMPLES::

        sage: from sage.misc.lazy_list import lazy_list
        sage: l = lazy_list(Primes())
        sage: l                                                                         # needs sage.libs.pari
        lazy list [2, 3, 5, ...]
        sage: l[200]                                                                    # needs sage.libs.pari
        1229
    """

    def __init__(self, cache=None, start=None, stop=None, step=None):
        r"""
        No check is performed on input and bad input can result in a Sage crash.
        You are advised to use the function :func:`lazy_list` instead. The only
        case where you might want to use directly this constructor is if you
        have a list that you want to wrap (without copy) into a lazy list.
        See in the example below.

        INPUT:

        - ``cache`` -- an optional list to be used as the cache. Be careful that
          there is no copy

        - ``start``, ``stop``, ``step`` -- for slices

        .. NOTE::

            Everywhere the constant ``PY_SSIZE_T_MAX`` plays the role of infinity

        EXAMPLES::

            sage: from sage.misc.lazy_list import lazy_list_generic
            sage: l = [0,1,2]
            sage: ll = lazy_list_generic(l, 0, 2, None)
            sage: ll
            lazy list [0, 1]

        The above code may be dangerous since the lazy list holds a reference
        to the initial list::

            sage: l[0] = 'haha'
            sage: ll
            lazy list ['haha', 1]
        """
        self.cache = [] if cache is None else cache
        self.start = 0 if start is None else start
        self.stop = PY_SSIZE_T_MAX if stop is None else stop
        self.step = 1 if step is None else step

    def list(self):
        r"""
        Return the list made of the elements of ``self``.

        .. NOTE::

            If the iterator is sufficiently large, this will build a list
            of length ``PY_SSIZE_T_MAX`` which should be beyond the capacity of
            your RAM!

        EXAMPLES::

            sage: from sage.misc.lazy_list import lazy_list
            sage: P = lazy_list(Primes())
            sage: P[2:143:5].list()                                                     # needs sage.libs.pari
            [5, 19, 41, 61, 83, 107, 137, 163, 191, 223, 241, 271, 307, 337, 367,
             397, 431, 457, 487, 521, 563, 593, 617, 647, 677, 719, 751, 787, 823]
            sage: P = lazy_list(iter([1,2,3]))
            sage: P.list()
            [1, 2, 3]
            sage: P[:100000].list()
            [1, 2, 3]
            sage: P[1:7:2].list()
            [2]

        TESTS:

        Check that the cache is immutable::

            sage: lazy = lazy_list(iter(Primes()))[:5]
            sage: l = lazy.list(); l                                                    # needs sage.libs.pari
            [2, 3, 5, 7, 11]
            sage: l[0] = -1; l                                                          # needs sage.libs.pari
            [-1, 3, 5, 7, 11]
            sage: lazy.list()                                                           # needs sage.libs.pari
            [2, 3, 5, 7, 11]
        """
        self._fit(self.stop - self.step)
        return self.cache[self.start:self.stop:self.step]

    def _info(self):
        r"""
        Print information about ``self`` on standard output.

        EXAMPLES::

            sage: from sage.misc.lazy_list import lazy_list
            sage: P = lazy_list(iter(Primes()))[10:21474838:4]
            sage: P._info()
            cache length 0
            start        10
            stop         21474838
            step         4
            sage: P[0]                                                                  # needs sage.libs.pari
            31
            sage: P._info()                                                             # needs sage.libs.pari
            cache length 11
            start        10
            stop         21474838
            step         4
        """
        print("cache length", len(self.cache))
        print("start       ", self.start)
        print("stop        ", self.stop)
        print("step        ", self.step)

    def __add__(self, other):
        r"""
        If ``self`` is a list then return the lazy_list that consists of the
        concatenation of ``self`` and ``other``.

        TESTS::

            sage: from sage.misc.lazy_list import lazy_list
            sage: from itertools import count
            sage: l = lazy_list(i**3 - i + 1 for i in count()); l
            lazy list [1, 1, 7, ...]
            sage: p = ['huit', 'douze']
            sage: ll = p + l; ll
            lazy list ['huit', 'douze', 1, ...]
            sage: l[:10].list() == ll[2:12].list()
            True
            sage: p
            ['huit', 'douze']
            sage: ([0,2] + lazy_list([0,1])).list()
            [0, 2, 0, 1]
        """
        if not isinstance(self, list):
            raise TypeError("can only add list to lazy_list")

        cdef lazy_list_from_iterator l = lazy_list_from_iterator.__new__(lazy_list_from_iterator)
        l.cache = self[:]
        l.start = 0
        l.stop = PY_SSIZE_T_MAX
        l.step = 1
        l.iterator = iter(other)
        return l

    def __repr__(self):
        r"""
        Return a string representation.

        To customize the string representation, the
        :func:`lazy_list_formatter` can be used.

        TESTS::

            sage: from sage.misc.lazy_list import lazy_list
            sage: from itertools import count
            sage: r = lazy_list(count()); r  # indirect doctest
            lazy list [0, 1, 2, ...]
            sage: r[:0]
            lazy list []
            sage: r[:1]
            lazy list [0]
            sage: r[:2]
            lazy list [0, 1]
            sage: r[:3]
            lazy list [0, 1, 2]
            sage: r[:4]
            lazy list [0, 1, 2, ...]
            sage: lazy_list([0,1])
            lazy list [0, 1]
            sage: lazy_list([0,1,2,3])
            lazy list [0, 1, 2, ...]
        """
        return lazy_list_formatter(self)

    def __reduce__(self):
        r"""
        Pickling support.

        EXAMPLES::

            sage: from itertools import count
            sage: from sage.misc.lazy_list import lazy_list
            sage: m = lazy_list(count())
            sage: x = loads(dumps(m))
            sage: y = iter(x)
            sage: print("{} {} {}".format(next(y), next(y), next(y)))
            0 1 2
            sage: m2 = m[3::2]
            sage: loads(dumps(m2))
            lazy list [3, 5, 7, ...]
        """
        if self.master is None:
            raise NotImplementedError
        return slice_unpickle, (self.master, self.start, self.stop, self.step)

    cpdef int _fit(self, Py_ssize_t n) except -1:
        r"""
        Fill the cache making the term at index ``n`` available.

        You can access the term at position ``n`` from the cache when it returns
        ``0``.

        OUTPUT:

        - ``1`` -- the lazy list is actually finite and shorter than ``n``

        - ``0`` -- you can safely access the term at position ``n`` after this call

        - ``-1`` -- to handle Python errors (you can ignore it in Python code)

        EXAMPLES::

            sage: from sage.misc.lazy_list import lazy_list
            sage: l = lazy_list(iter([0,1,2,-34,3,2,-5,12,1,4,-18,5,-12]))[2::3]
            sage: l._info()
            cache length 0
            start        2
            stop         9223372036854775807    # 64-bit
            stop         2147483647             # 32-bit
            step         3
            sage: l._fit(13)
            1
            sage: l._info()
            cache length 13
            start        2
            stop         14
            step         3

            sage: l = lazy_list([0]*12)[1::2]
            sage: l._fit(100)
            1
            sage: l._info()
            cache length 12
            start        1
            stop         13
            step         2
            sage: l._fit(100)
            1
        """
        if n > self.stop - self.step:
            return 1

        if self._update_cache_up_to(n):
            self.stop = min(self.stop, len(self.cache))
            if self.master is not None:
                self.stop = min(self.stop, self.master.stop)
            if self.stop <= self.start:
                self.start = self.stop = 0
                self.step = 1
            if (self.start - self.stop) % self.step:
                self.stop += self.step + (self.start - self.stop) % self.step
            return 1
        return 0

    cpdef get(self, Py_ssize_t i):
        r"""
        Return the element at position ``i``.

        If the index is not an integer, then raise a :exc:`TypeError`.  If the
        argument is negative then raise a :exc:`ValueError`.  Finally, if the
        argument is beyond the size of that lazy list it raises a
        :exc:`IndexError`.

        EXAMPLES::

            sage: from sage.misc.lazy_list import lazy_list
            sage: from itertools import chain, repeat
            sage: f = lazy_list(chain(iter([1,2,3]), repeat('a')))
            sage: f.get(0)
            1
            sage: f.get(3)
            'a'
            sage: f.get(0)
            1
            sage: f.get(4)
            'a'

            sage: g = f[:10]
            sage: g.get(5)
            'a'
            sage: g.get(10)
            Traceback (most recent call last):
            ...
            IndexError: lazy list index out of range
            sage: g.get(1/2)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert rational 1/2 to an integer
        """
        if i < 0:
            raise ValueError("indices must be nonnegative")

        i = self.start + i * self.step
        if self._fit(i):
            raise IndexError("lazy list index out of range")
        return self.cache[i]

    def __call__(self, i):
        r"""
        An alias for :meth:`get`.

        TESTS::

            sage: from sage.misc.lazy_list import lazy_list
            sage: from itertools import chain, repeat
            sage: f = lazy_list(chain(iter([1,2,3]), repeat('a')))
            sage: f(2)
            3
            sage: f(3)
            'a'
        """
        return self.get(i)

    def __iter__(self):
        r"""
        Return an iterator.

        TESTS::

            sage: from itertools import count
            sage: from sage.misc.lazy_list import lazy_list
            sage: iter(lazy_list(count()))
            <...generator object at 0x...>

        ::

            sage: l = lazy_list(i ** 2 for i in range(5))
            sage: list(l)
            [0, 1, 4, 9, 16]
            sage: l._info()
            cache length 5
            start        0
            stop         5
            step         1
        """
        cdef Py_ssize_t i

        i = self.start
        while i < self.stop:
            if self._fit(i):
                return
            yield self.cache[i]
            i += self.step

    def __getitem__(self, key):
        r"""
        Return a lazy list which shares the same cache.

        EXAMPLES::

            sage: from sage.misc.lazy_list import lazy_list
            sage: f = lazy_list(iter([1,2,3]))
            sage: f0 = f[0:]
            sage: [f.get(0), f.get(1), f.get(2)]
            [1, 2, 3]
            sage: f1 = f[1:]
            sage: [f1.get(0), f1.get(1)]
            [2, 3]
            sage: f2 = f[2:]
            sage: f2.get(0)
            3
            sage: f3 = f[3:]
            sage: f3.get(0)
            Traceback (most recent call last):
            ...
            IndexError: lazy list index out of range

            sage: l = lazy_list([0]*12)[1::2]
            sage: l[2::3]
            lazy list [0, 0]
            sage: l[3::2]
            lazy list [0, 0]

        A lazy list automatically adjusts the indices in order that start and
        stop are congruent modulo step::

            sage: P = lazy_list(iter(Primes()))
            sage: P[1:12:4]._info()
            cache length 0
            start        1
            stop         13
            step         4
            sage: P[1:13:4]._info()
            cache length 0
            start        1
            stop         13
            step         4
            sage: P[1:14:4]._info()
            cache length 0
            start        1
            stop         17
            step         4
            sage: Q = P[100:1042233:12]
            sage: Q._info()
            cache length 0
            start        100
            stop         1042240
            step         12
            sage: R = Q[233::3]
            sage: R._info()
            cache length 0
            start        2896
            stop         1042252
            step         36
            sage: 1042252%36 == 2896%36
            True

        We check commutation::

            sage: l = lazy_list(iter(range(10000)))
            sage: l1 = l[::2][:3001]
            sage: l2 = l[:6002][::2]
            sage: l1._info()
            cache length 0
            start        0
            stop         6002
            step         2
            sage: l2._info()
            cache length 0
            start        0
            stop         6002
            step         2
            sage: l3 = l1[13::2][:50:2]
            sage: l4 = l1[:200][13:113:4]
            sage: l3._info()
            cache length 0
            start        26
            stop         226
            step         8
            sage: l4._info()
            cache length 0
            start        26
            stop         226
            step         8

        Further tests::

            sage: l = lazy_list(iter([0]*25))
            sage: l[2::3][2::3][4::5]
            lazy list []
            sage: l[2::5][3::][1::]
            lazy list [0]
            sage: l[3:24:2][1::][1:7:3]
            lazy list [0, 0]
            sage: l[::2][2::][2::3]
            lazy list [0, 0, 0]
            sage: l[4:3][:] is l[18:2]   # *the* empty_lazy_list
            True
        """
        if not isinstance(key, slice):
            return self.get(key)

        # the following make all terms > 0
        cdef Py_ssize_t start, stop, step
        start = 0 if key.start is None else key.start
        stop = PY_SSIZE_T_MAX if key.stop is None else key.stop
        step = 1 if key.step is None else key.step

        if step == 0:
            raise TypeError("step may not be 0")
        if step < 0 or start < 0 or stop < 0:
            raise ValueError("slice indices must be nonnegative")

        step = step * self.step
        start = self.start + start * self.step
        if stop != PY_SSIZE_T_MAX:
            stop = self.start + stop * self.step
        if stop > self.stop:
            stop = self.stop
        if stop != PY_SSIZE_T_MAX and stop % step != start % step:
            stop = stop - (stop - start) % step + step

        if stop <= start:
            return empty_lazy_list

        # here we return a slice of self. That is to say, a lazy list which
        # shares the same cache of values
        cdef lazy_list_generic l = lazy_list_generic.__new__(lazy_list_generic)
        l.master = self
        l.cache = self.cache
        l.start = start
        l.stop = stop
        l.step = step

        return l

    cpdef int _update_cache_up_to(self, Py_ssize_t i) except -1:
        r"""
        Update the cache up to ``i``.

        This is the default implementation that calls ``_new_slice``.

        OUTPUT:

        - ``-1`` -- a Python error occurred

        - ``0`` -- the cache has now size larger than ``i``

        - ``1`` -- the lazy list is actually finite and shorter than ``i``

        TESTS::

            sage: from sage.misc.lazy_list import lazy_list
            sage: L = lazy_list(Primes())[2:]
            sage: L._update_cache_up_to(4)                                              # needs sage.libs.pari
            0
            sage: L._info()                                                             # needs sage.libs.pari
            cache length 5
            start        2
            stop         9223372036854775807    # 64-bit
            stop         2147483647             # 32-bit
            step         1
        """
        if self.master is not None:    # this is a slice
            return self.master._fit(i)

        cdef list l
        while len(self.cache) <= i:
            l = self._new_slice()
            if not l:
                return 1
            self.cache.extend(l)
        return 0

    cpdef list _get_cache_(self):
        r"""
        Return the internal cache.

        TESTS::

            sage: from sage.misc.lazy_list import lazy_list
            sage: L = lazy_list(Primes()); L                                            # needs sage.libs.pari
            lazy list [2, 3, 5, ...]
            sage: L._get_cache_()                                                       # needs sage.libs.pari
            [2, 3, 5, 7]
        """
        return self.cache


cdef class lazy_list_from_iterator(lazy_list_generic):
    r"""
    Lazy list built from an iterator.

    EXAMPLES::

        sage: from sage.misc.lazy_list import lazy_list
        sage: from itertools import count
        sage: m = lazy_list(count()); m
        lazy list [0, 1, 2, ...]

        sage: m2 = lazy_list(count())[8:20551:2]
        sage: m2
        lazy list [8, 10, 12, ...]

        sage: x = iter(m)
        sage: [next(x), next(x), next(x)]
        [0, 1, 2]
        sage: y = iter(m)
        sage: [next(y), next(y), next(y)]
        [0, 1, 2]
        sage: [next(x), next(y)]
        [3, 3]
        sage: loads(dumps(m))
        lazy list [0, 1, 2, ...]
    """

    def __init__(self, iterator, cache=None, stop=None):
        r"""
        INPUT:

        - ``iterator`` -- an iterator

        - ``cache`` -- an optional list to be used as the cache; be careful that
          there is no copy

        - ``stop`` -- an optional stop point

        TESTS::

            sage: from sage.misc.lazy_list import lazy_list_from_iterator
            sage: from itertools import count
            sage: lazy_list_from_iterator(count())
            lazy list [0, 1, 2, ...]
            sage: lazy_list_from_iterator(count(), ['a'], 10)
            lazy list ['a', 0, 1, ...]
            sage: _._info()
            cache length 4
            start        0
            stop         10
            step         1
        """
        self.iterator = iterator
        lazy_list_generic.__init__(self, cache, None, stop, None)

    cpdef int _update_cache_up_to(self, Py_ssize_t i) except -1:
        r"""
        Update the cache up to ``i``.

        OUTPUT:

        - ``-1`` -- a Python error occurred

        - ``0`` -- everything went fine

        - ``1`` -- the iterator stopped before ``i``

        TESTS::

            sage: from sage.misc.lazy_list import lazy_list
            sage: L = lazy_list(iter(Primes()))[2:]
            sage: L._update_cache_up_to(4)                                              # needs sage.libs.pari
            0
            sage: L._info()                                                             # needs sage.libs.pari
            cache length 5
            start        2
            stop         9223372036854775807    # 64-bit
            stop         2147483647             # 32-bit
            step         1
        """
        while len(self.cache) <= i:
            try:
                o = next(self.iterator)
            except StopIteration:
                return 1
            self.cache.append(o)
        return 0

    def __reduce__(self):
        r"""
        TESTS::

            sage: from sage.misc.lazy_list import lazy_list_from_iterator
            sage: from itertools import count
            sage: loads(dumps(lazy_list_from_iterator(count())))
            lazy list [0, 1, 2, ...]
            sage: loads(dumps(lazy_list_from_iterator(count(), ['a'])))
            lazy list ['a', 0, 1, ...]
        """
        return lazy_list_from_iterator, (self.iterator, self.cache, self.stop)


cdef class lazy_list_from_function(lazy_list_generic):

    def __init__(self, function, cache=None, stop=None):
        r"""
        INPUT:

        - ``function`` -- a function that maps ``n`` to the element
          at position ``n`` (this function only needs to be defined for length
          larger than the length of the cache)

        - ``cache`` -- an optional list to be used as the cache. Be careful that
          there is no copy

        - ``stop`` -- an optional integer to specify the length of this lazy list
          (Otherwise it is considered infinite)

        EXAMPLES::

            sage: from sage.misc.lazy_list import lazy_list_from_function
            sage: lazy_list_from_function(euler_phi)                                    # needs sage.libs.pari
            lazy list [0, 1, 1, ...]
            sage: lazy_list_from_function(divisors, [None])
            lazy list [None, [1], [1, 2], ...]

        TESTS::

            sage: def f(n):
            ....:     if n >= 5: raise StopIteration
            ....:     return 5 - n
            sage: list(lazy_list_from_function(f))
            [5, 4, 3, 2, 1]
        """
        self.callable = function
        lazy_list_generic.__init__(self, cache)

    cpdef int _update_cache_up_to(self, Py_ssize_t i) except -1:
        r"""
        Update the cache up to ``i``.

        OUTPUT:

        - ``-1`` -- a Python error occurred

        - ``0`` -- everything went fine

        - ``1`` -- the iterator stopped before ``i``

        TESTS::

            sage: from sage.misc.lazy_list import lazy_list
            sage: L = lazy_list(lambda x: 2*x)[2:]
            sage: L._update_cache_up_to(4)
            0
            sage: L._info()
            cache length 5
            start        2
            stop         9223372036854775807    # 64-bit
            stop         2147483647             # 32-bit
            step         1
        """
        while len(self.cache) <= i:
            try:
                value = self.callable(len(self.cache))
            except StopIteration:
                return 1
            self.cache.append(value)

    def __reduce__(self):
        r"""
        TESTS::

            sage: from sage.misc.lazy_list import lazy_list_from_function
            sage: loads(dumps(lazy_list_from_function(euler_phi)))                      # needs sage.libs.pari
            lazy list [0, 1, 1, ...]
            sage: loads(dumps(lazy_list_from_function(divisors, [None])))
            lazy list [None, [1], [1, 2], ...]
        """
        if self.start != 0 or self.step != 1:
            raise RuntimeError
        return lazy_list_from_function, (self.callable, self.cache, self.stop)


cdef class lazy_list_from_update_function(lazy_list_generic):

    def __init__(self, function, cache=None, stop=None):
        r"""
        INPUT:

        - ``function`` -- a function that updates a list of precomputed values
          The update function should take as input a list and make it longer
          (using either the methods ``append`` or ``extend``). If after a call
          to the update function the list of values is shorter a
          :exc:`RuntimeError` will occur. If no value is added then the lazy list
          is considered finite.

        - ``cache`` -- an optional list to be used as the cache. Be careful that
          there is no copy

        - ``stop`` -- an optional integer to specify the length of this lazy list
          (otherwise it is considered infinite)

        TESTS::

            sage: from sage.misc.lazy_list import lazy_list_from_update_function
            sage: def update_function(values):
            ....:     n = len(values)+1
            ....:     values.extend([n]*n)
            sage: l = lazy_list_from_update_function(update_function)
            sage: l[:20].list()
            [1, 2, 2, 4, 4, 4, 4, 8, 8, 8, 8, 8, 8, 8, 8, 16, 16, 16, 16, 16]
        """
        self.update_function = function
        lazy_list_generic.__init__(self, cache, None, stop, None)

    cpdef int _update_cache_up_to(self, Py_ssize_t i) except -1:
        r"""
        Update the cache up to ``i``.

        OUTPUT:

        - ``-1`` -- a Python error occurred

        - ``0`` -- everything went fine

        - ``1`` -- the iterator stopped before ``i``

        TESTS::

            sage: from sage.misc.lazy_list import lazy_list_from_update_function
            sage: def update_function(values):
            ....:     n = len(values)+1
            ....:     values.extend([n]*n)
            sage: L = lazy_list_from_update_function(update_function)[2:]
            sage: L._update_cache_up_to(4)
            0
            sage: L._info()
            cache length 7
            start        2
            stop         9223372036854775807    # 64-bit
            stop         2147483647             # 32-bit
            step         1
        """
        cdef Py_ssize_t l,ll
        l = len(self.cache)
        while l <= i:
            self.update_function(self.cache)
            ll = len(self.cache)
            if ll < l:
                raise RuntimeError("the update function made the cache shorter")
            elif l == ll:
                return 1
            l = ll
        return 0

    def __reduce__(self):
        r"""
        TESTS::

            sage: from sage.misc.lazy_list import lazy_list

            sage: def my_update_function(values): values.append(ZZ(len(values)).is_prime())
            sage: l = lazy_list(update_function=my_update_function)
            sage: l[4]
            False
            sage: loads(dumps(l))   # not tested (works in console though)
            lazy list [False, False, True, ...]

            sage: def say_hey(cache): print("hey")
            sage: l = lazy_list(update_function=say_hey, initial_values=range(10))
            sage: l._fit(10)
            hey
            1
            sage: l._info()
            cache length 10
            start        0
            stop         10
            step         1
            sage: l2 = loads(dumps(l))   # not tested
            sage: l2._info()             # not tested
            sage: l2._info()             # not tested
            cache length 10
            start        0
            stop         10
            step         1
            sage: l.list() == l2.list()  # not tested
            True
        """
        if self.start != 0 or self.step != 1:
            raise RuntimeError
        return lazy_list_from_update_function, (self.update_function, self.cache, self.stop)
