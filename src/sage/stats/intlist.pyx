# sage.doctest: optional - numpy
"""
C Int Lists

This is a class for fast basic operations with lists of C ints.  It is
similar to the double precision TimeSeries class.  It has all the
standard C int semantics, of course, including overflow.  It is also
similar to the Python list class, except all elements are C ints,
which makes some operations much, much faster.  For example,
concatenating two IntLists can be over 10 times faster than
concatenating the corresponding Python lists of ints, and taking
slices is also much faster.

AUTHOR:

- William Stein, 2010-03
"""

# ****************************************************************************
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

# Global parameter that sets the maximum number of entries of an IntList to print.
max_print = 10

cimport cython
from libc.string cimport memcpy
from cysignals.memory cimport sig_malloc, sig_free
from cysignals.signals cimport sig_on, sig_off

from sage.rings.integer import Integer
from sage.stats.time_series cimport TimeSeries
from cpython.bytes cimport PyBytes_FromStringAndSize, PyBytes_AsString
from sage.structure.richcmp cimport rich_to_bool


cdef class IntList:
    """
    A list of C int's.
    """
    def __cinit__(self):
        """
        Create new empty uninitialized IntList.

        EXAMPLES::

            sage: stats.IntList(5)     # indirect test
            [0, 0, 0, 0, 0]
        """
        self._values = NULL

    def __init__(self, values):
        """
        Create an initialized list of C ints.

        INPUT:

        - ``values`` -- int, long, Integer, list of integers, or a TimeSeries

        If the input is a time series or list of floats, then the
        integer parts of the entries are taken (not the floor).

        EXAMPLES::

            sage: stats.IntList(8)
            [0, 0, 0, 0, 0, 0, 0, 0]

            sage: stats.IntList([1,5,-39392])
            [1, 5, -39392]

        We check for overflow when creating the IntList::

            sage: stats.IntList([1, 3, 2^32])
            Traceback (most recent call last):
            ...
            OverflowError: ... too large to convert to C long  # 32-bit
            OverflowError: ... too large to convert to int     # 64-bit

        Printing omits entries::

            sage: stats.IntList(1000)
            [0, 0, 0, 0, 0 ... 0, 0, 0, 0, 0]

        Floats are truncated to their integer parts::

            sage: stats.IntList([1.1, -2.6])
            [1, -2]
            sage: stats.IntList(stats.TimeSeries([1.1, -2.6]))
            [1, -2]
        """
        cdef TimeSeries T
        if isinstance(values, (int, Integer)):
            self._length = values
            values = None
        elif isinstance(values, TimeSeries):
            T = values
            self._length = T._length
        else:
            self._length = len(values)

        self._values = <int*> sig_malloc(sizeof(int)*self._length)
        if self._values == NULL:
            raise MemoryError
        cdef Py_ssize_t i
        if values is None:
            for i in range(self._length):
                self._values[i] = 0
        elif isinstance(values, TimeSeries):
            for i in range(self._length):
                self._values[i] = <int> T._values[i]
        else:
            for i in range(self._length):
                self._values[i] = values[i]

    def __richcmp__(IntList self, other, int op):
        """
        Compare ``self`` and ``other``.  This has the same semantics
        as list comparison.

        EXAMPLES::

            sage: v = stats.IntList([1,2,3]); w = stats.IntList([1,2])
            sage: v < w
            False
            sage: w < v
            True
            sage: v == v
            True
            sage: w == w
            True
        """
        cdef IntList _other
        cdef Py_ssize_t c, i
        if not isinstance(other, IntList):
            _other = IntList(other)
        else:
            _other = <IntList>other
        for i in range(min(self._length, _other._length)):
            d = self._values[i] - _other._values[i]
            if d:
                return rich_to_bool(op, -1 if d < 0 else 1)
        c = self._length - _other._length
        if c:
            return rich_to_bool(op, -1 if c < 0 else 1)
        return rich_to_bool(op, 0)

    def __dealloc__(self):
        """
        Deallocate memory used by the IntList, if it was allocated.
        """
        sig_free(self._values)

    def __repr__(self):
        """
        Return string representation of this IntList.

        EXAMPLES::

            sage: a = stats.IntList([1..15]); a.__repr__()
            '[1, 2, 3, 4, 5 ... 11, 12, 13, 14, 15]'
            sage: sage.stats.intlist.max_print = 20
            sage: a.__repr__()
            '[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]'
            sage: sage.stats.intlist.max_print = 10
            sage: a.__repr__()
            '[1, 2, 3, 4, 5 ... 11, 12, 13, 14, 15]'
        """
        if len(self) > max_print:
            v0 = self[:max_print//2]
            v1 = self[-max_print//2:]
            return '[' + ', '.join(str(x) for x in v0) + ' ... ' + \
                         ', '.join(str(x) for x in v1) + ']'
        else:
            return str(self.list())

    def __getitem__(self, i):
        """
        Return i-th entry or slice of self, following standard Python
        semantics.  The returned slice is an intlist, and the returned
        entry is a Python int.

        INPUT:

        - ``i`` -- integer or slice

        EXAMPLES::

            sage: a = stats.IntList([0..9]); a
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            sage: a[5]
            5
            sage: a[-2]
            8
            sage: a[5:-2]
            [5, 6, 7]
            sage: type(a[5:-2])
            <... 'sage.stats.intlist.IntList'>
            sage: type(a[5])
            <... 'int'>
        """
        cdef Py_ssize_t start, stop, step, j
        cdef IntList t
        if isinstance(i, slice):
            start = 0 if (i.start is None) else i.start
            stop = self._length if (i.stop is None) else i.stop
            step = 1 if (i.step is None) else i.step
            if start < 0:
                start += self._length
                if start < 0: start = 0
            elif start >= self._length:
                start = self._length - 1
            if stop < 0:
                stop += self._length
                if stop < 0: stop = 0
            elif stop > self._length:
                stop = self._length
            if start >= stop:
                return new_int_list(0)
            if step < 0:
                step = -step
                t = new_int_list((stop-start)/step)
                for j from 0 <= j < (stop-start)/step:
                    t._values[j] = self._values[stop-1 - j*step]
            elif step > 1:
                t = new_int_list((stop-start)/step)
                for j from 0 <= j < (stop-start)/step:
                    t._values[j] = self._values[j*step+start]
            else:
                t = new_int_list(stop-start)
                memcpy(t._values, self._values + start, sizeof(int)*t._length)
            return t
        else:
            j = i
            if j < 0:
                j += self._length
                if j < 0:
                    raise IndexError("IntList index out of range")
            elif j >= self._length:
                raise IndexError("IntList index out of range")
            return self._values[j]

    def __setitem__(self, Py_ssize_t i, int x):
        """
        Set the i-th entry of self, following standard Python semantics.

        INPUT:

            - ``i`` -- integer
            - ``x`` -- integer

        EXAMPLES::

            sage: a = stats.IntList([-2,3,7,-4])
            sage: a[1] = 10393; a
            [-2, 10393, 7, -4]
            sage: a[-1] = -10; a
            [-2, 10393, 7, -10]
            sage: a[100]
            Traceback (most recent call last):
            ...
            IndexError: IntList index out of range
            sage: a[-100]
            Traceback (most recent call last):
            ...
            IndexError: IntList index out of range
        """
        if i < 0:
            i += self._length
            if i < 0:
                raise IndexError("index out of range")
        elif i >= self._length:
            raise IndexError("index out of range")
        self._values[i] = x

    def __reduce__(self):
        """
        Used in pickling int lists.

        EXAMPLES::

            sage: a = stats.IntList([-2,3,7,-4])
            sage: loads(dumps(a)) == a
            True

            sage: v = stats.IntList([1,-3])
            sage: v.__reduce__()
            (<cyfunction unpickle_intlist_v1 at ...>, (..., 2))
            sage: loads(dumps(v)) == v
            True

        Note that dumping and loading with compress ``False`` is much faster,
        though dumping with compress ``True`` can save a lot of space::

            sage: v = stats.IntList([1..10^5])
            sage: loads(dumps(v, compress=False),compress=False) == v
            True
        """
        buf = PyBytes_FromStringAndSize(<char*>self._values, self._length*sizeof(int)/sizeof(char))
        return unpickle_intlist_v1, (buf, self._length)

    def list(self):
        """
        Return Python list version of ``self`` with Python ints as entries.

        EXAMPLES::

            sage: a = stats.IntList([1..15]); a
            [1, 2, 3, 4, 5 ... 11, 12, 13, 14, 15]
            sage: a.list()
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
            sage: list(a) == a.list()
            True
            sage: type(a.list()[0])
            <... 'int'>
        """
        cdef Py_ssize_t i
        return [self._values[i] for i in range(self._length)]

    cpdef int sum(self) noexcept:
        """
        Return the sum of the entries of ``self``.

        EXAMPLES::

            sage: stats.IntList([1..100]).sum()
            5050

        Note that there can be overflow, since the entries are C ints::

            sage: a = stats.IntList([2^30,2^30]); a
            [1073741824, 1073741824]
            sage: a.sum()
            -2147483648
        """
        cdef Py_ssize_t i
        cdef int s=0
        sig_on()
        for i in range(self._length):
            s += self._values[i]
        sig_off()
        return s

    cpdef int prod(self) noexcept:
        """
        Return the product of the entries of ``self``.

        EXAMPLES::

            sage: a = stats.IntList([1..10]); a
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
            sage: a.prod()
            3628800
            sage: factorial(10)
            3628800

        Note that there can be overflow::

            sage: a = stats.IntList([2^30, 2]); a
            [1073741824, 2]
            sage: a.prod()
            -2147483648
         """
        cdef Py_ssize_t i
        cdef int s=1
        sig_on()
        for i in range(self._length):
            s *= self._values[i]
        sig_off()
        return s

    def __len__(self):
        """
        Return the number of entries in this time series.

        OUTPUT: Python integer

        EXAMPLES::

            sage: len(stats.IntList([1..15]))
            15
            sage: len(stats.IntList([]))
            0
            sage: len(stats.IntList(10^6))
            1000000
        """
        return self._length

    def __add__(left, right):
        """
        Concatenate the integer lists ``self`` and ``right``.

        EXAMPLES::

            sage: stats.IntList([-2,3,5]) + stats.IntList([1,1,17])
            [-2, 3, 5, 1, 1, 17]
        """
        if not isinstance(right, IntList):
            raise TypeError("right operand must be an int list")
        if not isinstance(left, IntList):
            raise TypeError("left operand must be an int list")
        cdef IntList R = right
        cdef IntList L = left
        cdef IntList t = new_int_list(L._length + R._length)
        memcpy(t._values, L._values, sizeof(int)*L._length)
        memcpy(t._values + L._length, R._values, sizeof(int)*R._length)
        return t

    def min(self, bint index=False):
        """
        Return the smallest value in this integer list.  If this
        series has length 0 we raise a :exc:`ValueError`.

        INPUT:

        - ``index`` -- boolean (default: ``False``); if ``True``, also return
          index of minimal entry

        OUTPUT:

        - ``float`` -- smallest value
        - ``integer`` -- index of smallest value; only returned if
          ``index=True``

        EXAMPLES::

            sage: v = stats.IntList([1,-4,3,-2,-4])
            sage: v.min()
            -4
            sage: v.min(index=True)
            (-4, 1)
        """
        if self._length == 0:
            raise ValueError("min() arg is an empty sequence")
        cdef Py_ssize_t i, j
        cdef int s = self._values[0]
        j = 0
        for i in range(1, self._length):
            if self._values[i] < s:
                s = self._values[i]
                j = i
        if index:
            return s, j
        else:
            return s

    def max(self, bint index=False):
        """
        Return the largest value in this time series. If this series
        has length 0 we raise a :exc:`ValueError`

        INPUT:

        - ``index`` -- boolean (default: ``False``); if ``True``, also return
          index of maximum entry

        OUTPUT:

        - integer -- largest value
        - integer -- index of largest value; only returned if ``index=True``

        EXAMPLES::

            sage: v = stats.IntList([1,-4,3,-2,-4,3])
            sage: v.max()
            3
            sage: v.max(index=True)
            (3, 2)
        """
        if self._length == 0:
            raise ValueError("max() arg is an empty sequence")
        cdef Py_ssize_t i, j = 0
        cdef int s = self._values[0]
        for i in range(1,self._length):
            if self._values[i] > s:
                s = self._values[i]
                j = i
        if index:
            return s, j
        else:
            return s

    def time_series(self):
        """
        Return :class:`TimeSeries` version of ``self``, which involves changing
        each entry to a double.

        EXAMPLES::

            sage: T = stats.IntList([-2,3,5]).time_series(); T
            [-2.0000, 3.0000, 5.0000]
            sage: type(T)
            <... 'sage.stats.time_series.TimeSeries'>
        """
        cdef TimeSeries T = TimeSeries.__new__(TimeSeries)
        # We just reach into the data structure underlying T, since we
        # want this function to be *very* fast.
        T._length = self._length
        T._values = <double*> sig_malloc(sizeof(double)*self._length)
        cdef Py_ssize_t i
        for i in range(self._length):
            T._values[i] = self._values[i]
        return T

    def plot(self, *args, **kwds):
        """
        Return a plot of this :class:`IntList`.

        This just constructs the
        corresponding double-precision floating point :class:`TimeSeries`
        object, passing on all arguments.

        EXAMPLES::

            sage: stats.IntList([3,7,19,-2]).plot()                                     # needs sage.plot
            Graphics object consisting of 1 graphics primitive
            sage: stats.IntList([3,7,19,-2]).plot(color='red',                          # needs sage.plot
            ....:                                 pointsize=50, points=True)
            Graphics object consisting of 1 graphics primitive
        """
        return self.time_series().plot(*args, **kwds)

    def plot_histogram(self, *args, **kwds):
        """
        Return a histogram plot of this :class:`IntList`.

        This just constructs
        the corresponding double-precision floating point :class:`TimeSeries` object,
        and plots it, passing on all arguments.

        EXAMPLES::

            sage: stats.IntList([1..15]).plot_histogram()                               # needs sage.plot
            Graphics object consisting of 50 graphics primitives
        """
        return self.time_series().plot_histogram(*args, **kwds)


cdef IntList new_int_list(Py_ssize_t length):
    """
    Function that is used internally to quickly create a new intlist
    without initializing any of the allocated memory.

    INPUT:

    - ``length`` -- nonnegative integer

    OUTPUT: an IntList
    """
    if length < 0:
        raise ValueError("length must be nonnegative")
    cdef IntList t = IntList.__new__(IntList)
    t._length = length
    t._values = <int*> sig_malloc(sizeof(int)*length)
    return t


@cython.binding(True)
def unpickle_intlist_v1(bytes v, Py_ssize_t n):
    """
    Version 1 unpickle method.

    INPUT:

    - ``v`` -- a raw char buffer

    EXAMPLES::

        sage: v = stats.IntList([1,2,3])
        sage: s = v.__reduce__()[1][0]
        sage: type(s) == type(b'')
        True
        sage: sage.stats.intlist.unpickle_intlist_v1(s, 3)
        [1, 2, 3]
        sage: sage.stats.intlist.unpickle_intlist_v1(s+s,6)
        [1, 2, 3, 1, 2, 3]
        sage: sage.stats.intlist.unpickle_intlist_v1(b'',0)
        []
    """
    cdef IntList t = new_int_list(n)
    memcpy(t._values, PyBytes_AsString(v), n*sizeof(int))
    return t
