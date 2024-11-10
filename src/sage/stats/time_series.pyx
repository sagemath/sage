# sage.doctest: optional - numpy
r"""
Time Series

This is a module for working with discrete floating point time series.
It is designed so that every operation is very fast, typically much
faster than with other generic code, e.g., Python lists of doubles or
even NumPy arrays.  The semantics of time series is more similar to
Python lists of doubles than Sage real double vectors or NumPy 1-D
arrays.   In particular, time series are not endowed with much
algebraic structure and are always mutable.

.. NOTE::

    NumPy arrays are faster at slicing, since slices return
    references, and NumPy arrays have strides.  However, this speed at
    slicing makes NumPy slower at certain other operations.

EXAMPLES::

    sage: set_random_seed(1)
    sage: t = stats.TimeSeries([random()-0.5 for _ in range(10)]); t
    [0.3294, 0.0959, -0.0706, -0.4646, 0.4311, 0.2275, -0.3840, -0.3528, -0.4119, -0.2933]
    sage: t.sums()
    [0.3294, 0.4253, 0.3547, -0.1099, 0.3212, 0.5487, 0.1647, -0.1882, -0.6001, -0.8933]
    sage: t.exponential_moving_average(0.7)
    [0.0000, 0.3294, 0.1660, 0.0003, -0.3251, 0.2042, 0.2205, -0.2027, -0.3078, -0.3807]
    sage: t.standard_deviation()
    0.33729638212891383
    sage: t.mean()
    -0.08933425506929439
    sage: t.variance()
    0.1137688493972542...

AUTHOR:

- William Stein
"""

# ****************************************************************************
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

cimport cython
from cpython.bytes cimport PyBytes_FromStringAndSize, PyBytes_AsString
from libc.math cimport exp, log, pow, sqrt
from libc.string cimport memcpy
from cysignals.memory cimport sig_malloc, sig_free
from sage.structure.richcmp cimport rich_to_bool

cimport numpy as cnumpy

from sage.misc.randstate cimport randstate, current_randstate
from sage.misc.persist import register_unpickle_override
from sage.rings.integer import Integer
from sage.rings.real_double import RDF
from sage.modules.vector_real_double_dense cimport Vector_real_double_dense

max_print = 10
digits = 4

cdef class TimeSeries:
    def __cinit__(self):
        r"""
        Create new empty uninitialized time series.

        EXAMPLES:

        This implicitly calls new::

            sage: stats.TimeSeries([1,3,-4,5])
            [1.0000, 3.0000, -4.0000, 5.0000]
        """
        self._values = NULL

    def __init__(self, values, bint initialize=True):
        r"""
        Initialize new time series.

        INPUT:

        - ``values`` -- integer (number of values) or an iterable of
          floats

        - ``initialize`` -- boolean (default: ``True``); if ``False``, do not
          bother to zero out the entries of the new time series.
          For large series that you are going to just fill in,
          this can be way faster.

        EXAMPLES:

        This implicitly calls init::

            sage: stats.TimeSeries([pi, 3, 18.2])                                       # needs sage.symbolic
            [3.1416, 3.0000, 18.2000]

        Conversion from a NumPy 1-D array, which is very fast::

            sage: v = stats.TimeSeries([1..5])
            sage: w = v.numpy()
            sage: stats.TimeSeries(w)
            [1.0000, 2.0000, 3.0000, 4.0000, 5.0000]

        Conversion from an n-dimensional NumPy array also works::

            sage: import numpy
            sage: v = numpy.array([[1,2], [3,4]], dtype=float); v
            array([[1., 2.],
                   [3., 4.]])
            sage: stats.TimeSeries(v)
            [1.0000, 2.0000, 3.0000, 4.0000]
            sage: stats.TimeSeries(v[:,0])
            [1.0000, 3.0000]
            sage: u = numpy.array([[1,2],[3,4]])
            sage: stats.TimeSeries(u)
            [1.0000, 2.0000, 3.0000, 4.0000]

        For speed purposes we don't initialize (so value is garbage)::

            sage: t = stats.TimeSeries(10, initialize=False)
        """
        cdef Vector_real_double_dense z
        cdef cnumpy.ndarray np
        cdef double *np_data
        cdef unsigned int j
        if isinstance(values, (int, Integer)):
            self._length = values
            values = None
        elif isinstance(values, (Vector_real_double_dense, cnumpy.ndarray)):
            if isinstance(values, Vector_real_double_dense):
                np  = values._vector_numpy
            else:
                np = values

            if np.ndim != 1:
                np = np.reshape([np.size])

            # Make the array be the correct type and have a C array
            # for a data structure.  If the array already is the
            # correct type and has a C array, nothing is done, so this
            # should be fast in the common case.
            np = np.astype('double')
            np = cnumpy.PyArray_GETCONTIGUOUS(np)
            np_data = <double*> cnumpy.PyArray_DATA(np)
            self._length = np.shape[0]
            self._values = <double*> sig_malloc(sizeof(double) * self._length)
            if self._values == NULL:
                raise MemoryError

            memcpy(self._values, np_data, sizeof(double)*self._length)
            return
        else:
            values = [float(x) for x in values]
            self._length = len(values)

        self._values = <double*> sig_malloc(sizeof(double) * self._length)
        if self._values == NULL:
            raise MemoryError
        if not initialize: return
        cdef Py_ssize_t i
        if values is not None:
            for i from 0 <= i < self._length:
                self._values[i] = values[i]
        else:
            for i from 0 <= i < self._length:
                self._values[i] = 0

    def __reduce__(self):
        r"""
        Used in pickling time series.

        EXAMPLES::

            sage: v = stats.TimeSeries([1,-3.5])
            sage: v.__reduce__()
            (<cyfunction unpickle_time_series_v1 at ...>, (..., 2))
            sage: loads(dumps(v)) == v
            True

        Note that dumping and loading with compress ``False`` is much faster,
        though dumping with compress ``True`` can save a lot of space::

            sage: v = stats.TimeSeries([1..10^5])
            sage: loads(dumps(v, compress=False),compress=False) == v
            True
        """
        buf = PyBytes_FromStringAndSize(<char*>self._values, self._length*sizeof(double)/sizeof(char))
        return unpickle_time_series_v1, (buf, self._length)

    def __richcmp__(TimeSeries self, other, int op):
        r"""
        Compare ``self`` and ``other``.  This has the same semantics
        as list comparison.

        EXAMPLES::

            sage: v = stats.TimeSeries([1,2,3]); w = stats.TimeSeries([1,2])
            sage: v < w
            False
            sage: w < v
            True
            sage: v == v
            True
            sage: w == w
            True
        """
        cdef Py_ssize_t i
        cdef double d
        if not isinstance(other, TimeSeries):
            return NotImplemented
        _other = <TimeSeries>other
        for i in range(min(self._length, _other._length)):
            d = self._values[i] - _other._values[i]
            if d:
                return rich_to_bool(op, -1 if d < 0 else 1)
        c = self._length - _other._length
        if c:
            return rich_to_bool(op, -1 if c < 0 else 1)
        return rich_to_bool(op, 0)

    def __dealloc__(self):
        r"""
        Free up memory used by a time series.

        EXAMPLES:

        This tests ``__dealloc__`` implicitly::

            sage: v = stats.TimeSeries([1,3,-4,5])
            sage: del v
        """
        sig_free(self._values)

    def vector(self):
        r"""
        Return real double vector whose entries are the values of this
        time series.  This is useful since vectors have standard
        algebraic structure and play well with matrices.

        OUTPUT: a real double vector

        EXAMPLES::

            sage: v = stats.TimeSeries([1..10])
            sage: v.vector()
            (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0)
        """
        V = RDF**self._length
        # A copy of the numpy array is made in the vector constructor
        cdef Vector_real_double_dense x = Vector_real_double_dense(V, self.numpy(copy=False))
        return x

    def __repr__(self):
        r"""
        Return string representation of ``self``.

        EXAMPLES::

            sage: v = stats.TimeSeries([1,3.1908439,-4,5.93932])
            sage: v.__repr__()
            '[1.0000, 3.1908, -4.0000, 5.9393]'

        By default 4 digits after the decimal point are displayed.  To
        change this, change ``sage.stats.time_series.digits``. ::

            sage: sage.stats.time_series.digits = 2
            sage: v.__repr__()
            '[1.00, 3.19, -4.00, 5.94]'
            sage: v
            [1.00, 3.19, -4.00, 5.94]
            sage: sage.stats.time_series.digits = 4
            sage: v
            [1.0000, 3.1908, -4.0000, 5.9393]
        """
        return self._repr()

    def _repr(self, prec=None):
        r"""
        Print representation of a time series.

        INPUT:

        - ``prec`` -- (default: ``None``) number of digits of precision or
          ``None``. If ``None`` use the default
          ``sage.stats.time_series.digits``.

        OUTPUT: string

        EXAMPLES::

            sage: v = stats.TimeSeries([1,3.1908439,-4,5.93932])
            sage: v._repr()
            '[1.0000, 3.1908, -4.0000, 5.9393]'
            sage: v._repr(10)
            '[1.0000000000, 3.1908439000, -4.0000000000, 5.9393200000]'
            sage: v._repr(2)
            '[1.00, 3.19, -4.00, 5.94]'
        """
        if prec is None: prec = digits
        format = '%.' + str(prec) + 'f'
        if len(self) > max_print:
            v0 = self[:max_print//2]
            v1 = self[-max_print//2:]
            return '[' + ', '.join(format % x for x in v0) + ' ... ' + \
                         ', '.join(format % x for x in v1) + ']'
        else:
            return '[' + ', '.join(format % x for x in self) + ']'

    def __len__(self):
        r"""
        Return the number of entries in this time series.

        OUTPUT: Python integer

        EXAMPLES::

            sage: v = stats.TimeSeries([1,3.1908439,-4,5.93932])
            sage: v.__len__()
            4
            sage: len(v)
            4
        """
        return self._length

    def __getitem__(self, i):
        r"""
        Return `i`-th entry or slice of ``self``.

        EXAMPLES::

            sage: v = stats.TimeSeries([1,-4,3,-2.5,-4,3])
            sage: v[2]
            3.0
            sage: v[-1]
            3.0
            sage: v[-10]
            Traceback (most recent call last):
            ...
            IndexError: TimeSeries index out of range
            sage: v[5]
            3.0
            sage: v[6]
            Traceback (most recent call last):
            ...
            IndexError: TimeSeries index out of range

        Some slice examples::

            sage: v[-3:]
            [-2.5000, -4.0000, 3.0000]
            sage: v[-3:-1]
            [-2.5000, -4.0000]
            sage: v[::2]
            [1.0000, 3.0000, -4.0000]
            sage: v[3:20]
            [-2.5000, -4.0000, 3.0000]
            sage: v[3:2]
            []

        Make a copy::

            sage: v[:]
            [1.0000, -4.0000, 3.0000, -2.5000, -4.0000, 3.0000]

        Reverse the time series::

            sage: v[::-1]
            [3.0000, -4.0000, -2.5000, 3.0000, -4.0000, 1.0000]
        """
        cdef Py_ssize_t start, stop, step, j
        cdef TimeSeries t
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
                return new_time_series(0)
            if step < 0:
                step = -step
                t = new_time_series((stop-start)/step)
                for j from 0 <= j < (stop-start)/step:
                    t._values[j] = self._values[stop-1 - j*step]
            elif step > 1:
                t = new_time_series((stop-start)/step)
                for j from 0 <= j < (stop-start)/step:
                    t._values[j] = self._values[j*step+start]
            else:
                t = new_time_series(stop-start)
                memcpy(t._values, self._values + start, sizeof(double)*t._length)
            return t
        else:
            j = i
            if j < 0:
                j += self._length
                if j < 0:
                    raise IndexError("TimeSeries index out of range")
            elif j >= self._length:
                raise IndexError("TimeSeries index out of range")
            return self._values[j]

    def __setitem__(self, Py_ssize_t i, double x):
        r"""
        Set the `i`-th entry of ``self`` to ``x``.

        INPUT:

        - ``i`` -- nonnegative integer

        - ``x`` -- a float

        EXAMPLES::

            sage: v = stats.TimeSeries([1,3,-4,5.93932]); v
            [1.0000, 3.0000, -4.0000, 5.9393]
            sage: v[0] = -5.5; v
            [-5.5000, 3.0000, -4.0000, 5.9393]
            sage: v[-1] = 3.2; v
            [-5.5000, 3.0000, -4.0000, 3.2000]
            sage: v[10]
            Traceback (most recent call last):
            ...
            IndexError: TimeSeries index out of range
            sage: v[-5]
            Traceback (most recent call last):
            ...
            IndexError: TimeSeries index out of range
        """
        if i < 0:
            i += self._length
            if i < 0:
                raise IndexError("TimeSeries index out of range")
        elif i >= self._length:
            raise IndexError("TimeSeries index out of range")
        self._values[i] = x

    def __copy__(self):
        r"""
        Return a copy of ``self``.

        EXAMPLES::

            sage: v = stats.TimeSeries([1,-4,3,-2.5,-4,3])
            sage: v.__copy__()
            [1.0000, -4.0000, 3.0000, -2.5000, -4.0000, 3.0000]
            sage: copy(v)
            [1.0000, -4.0000, 3.0000, -2.5000, -4.0000, 3.0000]
            sage: copy(v) is v
            False
        """
        cdef Py_ssize_t i
        cdef TimeSeries t = new_time_series(self._length)
        memcpy(t._values, self._values , sizeof(double)*self._length)
        return t

    def __add__(left, right):
        r"""
        Concatenate the time series ``self`` and ``right``.

        .. NOTE::

            To add a single number to the entries of a time series,
            use the ``add_scalar`` method, and to add componentwise use
            the ``add_entries`` method.

        INPUT:

        - ``right`` -- a time series

        OUTPUT: a time series

        EXAMPLES::

            sage: v = stats.TimeSeries([1,2,3]); w = stats.TimeSeries([1,2])
            sage: v + w
            [1.0000, 2.0000, 3.0000, 1.0000, 2.0000]
            sage: v = stats.TimeSeries([1,2,-5]); v
            [1.0000, 2.0000, -5.0000]

        Note that both summands must be a time series::

            sage: v + list(range(4))
            Traceback (most recent call last):
            ...
            TypeError: right operand must be a time series
            sage: [1,5] + v
            Traceback (most recent call last):
            ...
            TypeError: left operand must be a time series
        """
        if not isinstance(right, TimeSeries):
            raise TypeError("right operand must be a time series")
        if not isinstance(left, TimeSeries):
            raise TypeError("left operand must be a time series")
        cdef TimeSeries R = right
        cdef TimeSeries L = left
        cdef TimeSeries t = new_time_series(L._length + R._length)
        memcpy(t._values, L._values, sizeof(double)*L._length)
        memcpy(t._values + L._length, R._values, sizeof(double)*R._length)
        return t

    def __mul__(left, right):
        r"""
        Multiply a time series by an integer n, which (like for lists)
        results in the time series concatenated with itself n times.

        .. NOTE::

            To multiply all the entries of a time series by a single
            scalar, use the ``scale`` method.

        INPUT:

        - ``left``, ``right`` -- integer and a time series

        OUTPUT: a time series

        EXAMPLES::

            sage: v = stats.TimeSeries([1,2,-5]); v
            [1.0000, 2.0000, -5.0000]
            sage: v*3
            [1.0000, 2.0000, -5.0000, 1.0000, 2.0000, -5.0000, 1.0000, 2.0000, -5.0000]
            sage: 3*v
            [1.0000, 2.0000, -5.0000, 1.0000, 2.0000, -5.0000, 1.0000, 2.0000, -5.0000]
            sage: v*v
            Traceback (most recent call last):
            ...
            TypeError: 'sage.stats.time_series.TimeSeries' object cannot be interpreted as an integer
        """
        cdef Py_ssize_t n, i
        cdef TimeSeries T
        if isinstance(left, TimeSeries):
            T = left
            n = right
        else:
            T = right
            n = left
        # Make n copies of T concatenated together
        cdef TimeSeries v = new_time_series(T._length * n)
        for i from 0 <= i < n:
            memcpy(v._values + i*T._length, T._values, sizeof(double)*T._length)
        return v

    def autoregressive_fit(self,M):
        r"""
        This method fits the time series to an autoregressive process
        of order ``M``. That is, we assume the process is given by
        `X_t-\mu=a_1(X_{t-1}-\mu)+a_2(X_{t-1}-\mu)+\cdots+a_M(X_{t-M}-\mu)+Z_t`
        where `\mu` is the mean of the process and `Z_t` is noise.
        This method returns estimates for `a_1,\dots,a_M`.

        The method works by solving the Yule-Walker equations
        `\Gamma a =\gamma`, where `\gamma=(\gamma(1),\dots,\gamma(M))`,
        `a=(a_1,\dots,a_M)`  with `\gamma(i)` the autocovariance of lag `i`
        and `\Gamma_{ij}=\gamma(i-j)`.


        .. WARNING::

            The input sequence is assumed to be stationary, which
            means that the autocovariance `\langle y_j y_k \rangle` depends
            only on the difference `|j-k|`.

        INPUT:

        - ``M`` -- integer

        OUTPUT: a time series -- the coefficients of the autoregressive process

        EXAMPLES::

            sage: set_random_seed(0)
            sage: v = stats.TimeSeries(10^4).randomize('normal').sums()
            sage: F = v.autoregressive_fit(100)
            sage: v
            [0.6767, 0.2756, 0.6332, 0.0469, -0.8897 ... 87.6759, 87.6825, 87.4120, 87.6639, 86.3194]
            sage: v.autoregressive_forecast(F)
            86.0177285042...
            sage: F
            [1.0148, -0.0029, -0.0105, 0.0067, -0.0232 ... -0.0106, -0.0068, 0.0085, -0.0131, 0.0092]

            sage: set_random_seed(0)
            sage: t = stats.TimeSeries(2000)
            sage: z = stats.TimeSeries(2000)
            sage: z.randomize('normal',1)
            [1.6767, 0.5989, 1.3576, 0.4136, 0.0635 ... 1.0057, -1.1467, 1.2809, 1.5705, 1.1095]
            sage: t[0]=1
            sage: t[1]=2
            sage: for i in range(2,2000):
            ....:     t[i]=t[i-1]-0.5*t[i-2]+z[i]
            sage: c=t[0:-1].autoregressive_fit(2)  #recovers recurrence relation
            sage: c #should be close to [1,-0.5]
            [1.0371, -0.5199]
        """
        acvs = [self.autocovariance(i) for i in range(M+1)]
        return autoregressive_fit(acvs)

    def autoregressive_forecast(self, filter):
        r"""
        Given the autoregression coefficients as outputted by the
        :meth:`autoregressive_fit` command, compute the forecast for the next
        term in the series.

        INPUT:

        - ``filter`` -- a time series outputted by the ``autoregressive_fit``
          command

        EXAMPLES::

            sage: set_random_seed(0)
            sage: v = stats.TimeSeries(100).randomize('normal').sums()
            sage: F = v[:-1].autoregressive_fit(5); F
            [1.0019, -0.0524, -0.0643, 0.1323, -0.0539]
            sage: v.autoregressive_forecast(F)
            11.7820298611...
            sage: v
            [0.6767, 0.2756, 0.6332, 0.0469, -0.8897 ... 9.2447, 9.6709, 10.4037, 10.4836, 12.1960]
        """
        cdef TimeSeries filt
        if isinstance(filter, TimeSeries):
            filt = filter
        else:
            filt = TimeSeries(filter)

        cdef double f = 0
        cdef Py_ssize_t i
        for i from 0 <= i < min(self._length, filt._length):
            f += self._values[self._length - i - 1] * filt._values[i]
        return f

    def reversed(self):
        r"""
        Return new time series obtain from this time series by
        reversing the order of the entries in this time series.

        OUTPUT: a time series

        EXAMPLES::

            sage: v = stats.TimeSeries([1..5])
            sage: v.reversed()
            [5.0000, 4.0000, 3.0000, 2.0000, 1.0000]
        """
        cdef Py_ssize_t i, n = self._length-1
        cdef TimeSeries t = new_time_series(self._length)

        for i from 0 <= i < self._length:
            t._values[i] = self._values[n - i]
        return t

    def extend(self, right):
        r"""
        Extend this time series by appending elements from the iterable
        ``right``.

        INPUT:

        - ``right`` -- iterable that can be converted to a time series

        EXAMPLES::

            sage: v = stats.TimeSeries([1,2,-5]); v
            [1.0000, 2.0000, -5.0000]
            sage: v.extend([-3.5, 2])
            sage: v
            [1.0000, 2.0000, -5.0000, -3.5000, 2.0000]
        """
        if not isinstance(right, TimeSeries):
            right = TimeSeries(right)
        if not right:
            return
        cdef TimeSeries T = right
        cdef double* z = <double*> sig_malloc(sizeof(double)*(self._length + T._length))
        if z == NULL:
            raise MemoryError
        memcpy(z, self._values, sizeof(double)*self._length)
        memcpy(z + self._length, T._values, sizeof(double)*T._length)
        sig_free(self._values)
        self._values = z
        self._length = self._length + T._length

    def list(self):
        r"""
        Return list of elements of ``self``.

        EXAMPLES::

            sage: v = stats.TimeSeries([1,-4,3,-2.5,-4,3])
            sage: v.list()
            [1.0, -4.0, 3.0, -2.5, -4.0, 3.0]
        """
        cdef Py_ssize_t i
        return [self._values[i] for i in range(self._length)]

    def log(self):
        r"""
        Return new time series got by taking the logarithms of all the
        terms in the time series.

        OUTPUT: a new time series

        EXAMPLES:

        We exponentiate then log a time series and get back
        the original series::

            sage: v = stats.TimeSeries([1,-4,3,-2.5,-4,3]); v
            [1.0000, -4.0000, 3.0000, -2.5000, -4.0000, 3.0000]
            sage: v.exp()
            [2.7183, 0.0183, 20.0855, 0.0821, 0.0183, 20.0855]
            sage: v.exp().log()
            [1.0000, -4.0000, 3.0000, -2.5000, -4.0000, 3.0000]

        Log of 0 gives ``-inf``::

            sage: stats.TimeSeries([1,0,3]).log()[1]
            -inf
        """
        cdef Py_ssize_t i
        cdef TimeSeries t = new_time_series(self._length)
        for i from 0 <= i < self._length:
            t._values[i] = log(self._values[i])
        return t

    def exp(self):
        r"""
        Return new time series got by applying the exponential map to
        all the terms in the time series.

        OUTPUT: a new time series

        EXAMPLES::

            sage: v = stats.TimeSeries([1..5]); v
            [1.0000, 2.0000, 3.0000, 4.0000, 5.0000]
            sage: v.exp()
            [2.7183, 7.3891, 20.0855, 54.5982, 148.4132]
            sage: v.exp().log()
            [1.0000, 2.0000, 3.0000, 4.0000, 5.0000]
        """
        cdef Py_ssize_t i
        cdef TimeSeries t = new_time_series(self._length)
        for i from 0 <= i < self._length:
            t._values[i] = exp(self._values[i])
        return t

    def abs(self):
        r"""
        Return new time series got by replacing all entries
        of ``self`` by their absolute value.

        OUTPUT: a new time series

        EXAMPLES::

            sage: v = stats.TimeSeries([1,3.1908439,-4,5.93932])
            sage: v
            [1.0000, 3.1908, -4.0000, 5.9393]
            sage: v.abs()
            [1.0000, 3.1908, 4.0000, 5.9393]
        """
        cdef Py_ssize_t i
        cdef TimeSeries t = new_time_series(self._length)
        for i from 0 <= i < self._length:
            t._values[i] = self._values[i] if self._values[i] >= 0 else -self._values[i]
        return t

    def diffs(self, Py_ssize_t k = 1):
        r"""
        Return the new time series got by taking the differences of
        successive terms in the time series.  So if ``self`` is the time
        series `X_0, X_1, X_2, \dots`, then this function outputs the
        series `X_1 - X_0, X_2 - X_1, \dots`.  The output series has one
        less term than the input series.  If the optional parameter
        `k` is given, return `X_k - X_0, X_{2k} - X_k, \dots`.

        INPUT:

        - ``k`` -- positive integer (default: 1)

        OUTPUT: a new time series

        EXAMPLES::

            sage: v = stats.TimeSeries([5,4,1.3,2,8]); v
            [5.0000, 4.0000, 1.3000, 2.0000, 8.0000]
            sage: v.diffs()
            [-1.0000, -2.7000, 0.7000, 6.0000]
        """
        if k != 1:
            return self.scale_time(k).diffs()
        cdef Py_ssize_t i
        cdef TimeSeries t = new_time_series(self._length - 1)
        for i from 1 <= i < self._length:
            t._values[i-1] = self._values[i] - self._values[i-1]
        return t

    def scale_time(self, Py_ssize_t k):
        r"""
        Return the new time series at scale ``k``.  If the input
        time series is `X_0, X_1, X_2, \dots`, then this function
        returns the shorter time series `X_0, X_k, X_{2k}, \dots`.

        INPUT:

        - ``k`` -- positive integer

        OUTPUT: a new time series

        EXAMPLES::

            sage: v = stats.TimeSeries([5,4,1.3,2,8,10,3,-5]); v
            [5.0000, 4.0000, 1.3000, 2.0000, 8.0000, 10.0000, 3.0000, -5.0000]
            sage: v.scale_time(1)
            [5.0000, 4.0000, 1.3000, 2.0000, 8.0000, 10.0000, 3.0000, -5.0000]
            sage: v.scale_time(2)
            [5.0000, 1.3000, 8.0000, 3.0000]
            sage: v.scale_time(3)
            [5.0000, 2.0000]
            sage: v.scale_time(10)
            []

        A series of odd length::

            sage: v = stats.TimeSeries([1..5]); v
            [1.0000, 2.0000, 3.0000, 4.0000, 5.0000]
            sage: v.scale_time(2)
            [1.0000, 3.0000, 5.0000]

        TESTS::

            sage: v.scale_time(0)
            Traceback (most recent call last):
            ...
            ValueError: k must be positive
            sage: v.scale_time(-1)
            Traceback (most recent call last):
            ...
            ValueError: k must be positive
        """
        if k <= 0:
            raise ValueError("k must be positive")

        cdef Py_ssize_t i, n
        n = self._length/k
        if self._length % 2:
            n += 1
        cdef TimeSeries t = new_time_series(n)  # in C / is floor division.
        for i from 0 <= i < n:
            t._values[i] = self._values[i*k]
        return t

    cpdef rescale(self, double s):
        r"""
        Change ``self`` by multiplying every value in the series by ``s``.

        INPUT:

        - ``s`` -- a float

        EXAMPLES::

            sage: v = stats.TimeSeries([5,4,1.3,2,8,10,3,-5]); v
            [5.0000, 4.0000, 1.3000, 2.0000, 8.0000, 10.0000, 3.0000, -5.0000]
            sage: v.rescale(0.5)
            sage: v
            [2.5000, 2.0000, 0.6500, 1.0000, 4.0000, 5.0000, 1.5000, -2.5000]
        """
        for i from 0 <= i < self._length:
            self._values[i] = self._values[i] * s

    def scale(self, double s):
        r"""
        Return new time series obtained by multiplying every value in the
        series by ``s``.

        INPUT:

        - ``s`` -- a float

        OUTPUT: a new time series with all values multiplied by ``s``

        EXAMPLES::

            sage: v = stats.TimeSeries([5,4,1.3,2,8,10,3,-5]); v
            [5.0000, 4.0000, 1.3000, 2.0000, 8.0000, 10.0000, 3.0000, -5.0000]
            sage: v.scale(0.5)
            [2.5000, 2.0000, 0.6500, 1.0000, 4.0000, 5.0000, 1.5000, -2.5000]
        """
        cdef TimeSeries t = new_time_series(self._length)
        for i from 0 <= i < self._length:
            t._values[i] = self._values[i] * s
        return t

    def add_scalar(self, double s):
        r"""
        Return new time series obtained by adding a scalar to every
        value in the series.

        .. NOTE::

            To add componentwise, use the :meth:`add_entries` method.

        INPUT:

        - ``s`` -- a float

        OUTPUT: a new time series with ``s`` added to all values

        EXAMPLES::

            sage: v = stats.TimeSeries([5,4,1.3,2,8,10,3,-5]); v
            [5.0000, 4.0000, 1.3000, 2.0000, 8.0000, 10.0000, 3.0000, -5.0000]
            sage: v.add_scalar(0.5)
            [5.5000, 4.5000, 1.8000, 2.5000, 8.5000, 10.5000, 3.5000, -4.5000]
        """
        cdef TimeSeries t = new_time_series(self._length)
        for i from 0 <= i < self._length:
            t._values[i] = self._values[i] + s
        return t

    def add_entries(self, t):
        r"""
        Add corresponding entries of ``self`` and ``t`` together,
        extending either ``self`` or ``t`` by 0s if they do
        not have the same length.

        .. NOTE::

            To add a single number to the entries of a time series,
            use the ``add_scalar`` method.

        INPUT:

        - ``t`` -- a time series

        OUTPUT:

        A time series with length the maximum of the lengths of
        ``self`` and ``t``.

        EXAMPLES::

            sage: v = stats.TimeSeries([1,2,-5]); v
            [1.0000, 2.0000, -5.0000]
            sage: v.add_entries([3,4])
            [4.0000, 6.0000, -5.0000]
            sage: v.add_entries(v)
            [2.0000, 4.0000, -10.0000]
            sage: v.add_entries([3,4,7,18.5])
            [4.0000, 6.0000, 2.0000, 18.5000]
        """
        if not isinstance(t, TimeSeries):
            t = TimeSeries(t)
        cdef Py_ssize_t i, n
        cdef TimeSeries T = t, shorter, longer
        if self._length > T._length:
            shorter = T
            longer = self
        else:
            shorter = self
            longer = T

        n = longer._length
        cdef TimeSeries v = new_time_series(n)

        # Iterate up to the length of the shorter one
        for i from 0 <= i < shorter._length:
            v._values[i] = shorter._values[i] + longer._values[i]

        # Copy over the rest of the values
        if shorter._length != n:
            memcpy(v._values + shorter._length,
                   longer._values + shorter._length,
                   sizeof(double)*(v._length - shorter._length))

        return v

    def show(self, *args, **kwds):
        r"""
        Return a plot of this time series.

        This is an alias of :meth:`plot`.

        EXAMPLES:

        Draw a plot of a time series::

            sage: stats.TimeSeries([1..10]).show()                                      # needs sage.plot
            Graphics object consisting of 1 graphics primitive
        """
        return self.plot(*args, **kwds)

    def plot(self, Py_ssize_t plot_points=1000, points=False, **kwds):
        r"""
        Return a plot of this time series.

        The plot shows a line or points through
        `(i,T(i))`, where `i` ranges over nonnegative integers up to the
        length of ``self``.

        INPUT:

        - ``plot_points`` -- (default: 1000) 0 or positive integer. Only
          plot the given number of equally spaced points in the time series.
          If 0, plot all points.

        - ``points`` -- boolean (default: ``False``); if ``True``, return just
          the points of the time series

        - ``**kwds`` -- passed to the line or point command

        EXAMPLES::

            sage: # needs sage.plot
            sage: v = stats.TimeSeries([5,4,1.3,2,8,10,3,-5]); v
            [5.0000, 4.0000, 1.3000, 2.0000, 8.0000, 10.0000, 3.0000, -5.0000]
            sage: v.plot()
            Graphics object consisting of 1 graphics primitive
            sage: v.plot(points=True)
            Graphics object consisting of 1 graphics primitive
            sage: v.plot() + v.plot(points=True, rgbcolor='red')
            Graphics object consisting of 2 graphics primitives
            sage: v.plot() + v.plot(points=True, rgbcolor='red', pointsize=50)
            Graphics object consisting of 2 graphics primitives
        """
        from sage.plot.all import line, point
        cdef Py_ssize_t s

        if self._length < plot_points:
            plot_points = 0

        if plot_points > 0:
            s = self._length/plot_points
            if plot_points <= 0:
                raise ValueError("plot_points must be a positive integer")
            v = self.scale_time(s).list()[:plot_points]
        else:
            s = 1
            v = self.list()
        w = [(i * s, y) for i,y in enumerate(v)]
        if points:
            L = point(w, **kwds)
        else:
            L = line(w, **kwds)
        L.axes_range(ymin=min(v), ymax=max(v), xmin=0, xmax=len(v)*s)
        return L

    def simple_moving_average(self, Py_ssize_t k):
        r"""
        Return the moving average time series over the last `k` time units.

        Assumes the input time series was constant with its starting value
        for negative time.  The `t`-th step of the output is the sum of
        the previous `k - 1` steps of ``self`` and the `k`-th step
        divided by `k`. Thus `k` values are averaged at each point.

        INPUT:

        - ``k`` -- positive integer

        OUTPUT: a time series with the same number of steps as ``self``

        EXAMPLES::

            sage: v = stats.TimeSeries([1,1,1,2,3]); v
            [1.0000, 1.0000, 1.0000, 2.0000, 3.0000]
            sage: v.simple_moving_average(0)
            [1.0000, 1.0000, 1.0000, 2.0000, 3.0000]
            sage: v.simple_moving_average(1)
            [1.0000, 1.0000, 1.0000, 2.0000, 3.0000]
            sage: v.simple_moving_average(2)
            [1.0000, 1.0000, 1.0000, 1.5000, 2.5000]
            sage: v.simple_moving_average(3)
            [1.0000, 1.0000, 1.0000, 1.3333, 2.0000]
        """
        if k == 0 or k == 1:
            return self.__copy__()
        if k <= 0:
            raise ValueError("k must be positive")
        cdef Py_ssize_t i
        cdef TimeSeries t = new_time_series(self._length)
        if self._length == 0:
            return t

        cdef double s = self._values[0] * k
        for i from 0 <= i < self._length:
            if i >= k:
                s -= self._values[i-k]
            else:
                s -= self._values[0]
            # I have a slight concern about accumulated rounding error given how
            # this algorithm adds and subtracts.
            s += self._values[i]
            t._values[i] = s/k
        return t

    def exponential_moving_average(self, double alpha):
        r"""
        Return the exponential moving average time series.

        Assumes the input time series was constant with its starting
        value for negative time.  The `t`-th step of the output is the
        sum of the previous `k-1` steps of ``self`` and the `k`-th
        step divided by `k`.

        The 0-th term is formally undefined, so we define it to be 0,
        and we define the first term to be ``self[0]``.

        INPUT:

        - ``alpha`` -- float; a smoothing factor with ``0 <= alpha <= 1``

        OUTPUT: a time series with the same number of steps as ``self``

        EXAMPLES::

            sage: v = stats.TimeSeries([1,1,1,2,3]); v
            [1.0000, 1.0000, 1.0000, 2.0000, 3.0000]
            sage: v.exponential_moving_average(0)
            [0.0000, 1.0000, 1.0000, 1.0000, 1.0000]
            sage: v.exponential_moving_average(1)
            [0.0000, 1.0000, 1.0000, 1.0000, 2.0000]
            sage: v.exponential_moving_average(0.5)
            [0.0000, 1.0000, 1.0000, 1.0000, 1.5000]

        Some more examples::

            sage: v = stats.TimeSeries([1,2,3,4,5])
            sage: v.exponential_moving_average(1)
            [0.0000, 1.0000, 2.0000, 3.0000, 4.0000]
            sage: v.exponential_moving_average(0)
            [0.0000, 1.0000, 1.0000, 1.0000, 1.0000]
        """
        if alpha < 0 or alpha > 1:
            raise ValueError("alpha must be between 0 and 1")
        cdef Py_ssize_t i
        cdef TimeSeries t = new_time_series(self._length)
        if self._length == 0:
            return t
        t._values[0] = 0
        if self._length == 1:
            return t
        t._values[1] = self._values[0]
        for i from 2 <= i < self._length:
            t._values[i] = alpha * self._values[i-1] + (1-alpha) *t._values[i-1]
        return t

    def sums(self, double s=0):
        r"""
        Return the new time series got by taking the running partial
        sums of the terms of this time series.

        INPUT:

        - ``s`` -- starting value for partial sums

        OUTPUT: a time series

        EXAMPLES::

            sage: v = stats.TimeSeries([1,1,1,2,3]); v
            [1.0000, 1.0000, 1.0000, 2.0000, 3.0000]
            sage: v.sums()
            [1.0000, 2.0000, 3.0000, 5.0000, 8.0000]
        """
        cdef Py_ssize_t i
        cdef TimeSeries t = new_time_series(self._length)
        for i from 0 <= i < self._length:
            s += self._values[i]
            t._values[i] = s
        return t

    cpdef double sum(self) noexcept:
        r"""
        Return the sum of all the entries of ``self``.

        If ``self`` has length 0, returns 0.

        OUTPUT: double

        EXAMPLES::

            sage: v = stats.TimeSeries([1,1,1,2,3]); v
            [1.0000, 1.0000, 1.0000, 2.0000, 3.0000]
            sage: v.sum()
            8.0
        """
        cdef double s = 0
        cdef Py_ssize_t i
        for i from 0 <= i < self._length:
            s += self._values[i]
        return s

    def prod(self):
        r"""
        Return the product of all the entries of ``self``.

        If ``self`` has length 0, returns 1.

        OUTPUT: double

        EXAMPLES::

            sage: v = stats.TimeSeries([1,1,1,2,3]); v
            [1.0000, 1.0000, 1.0000, 2.0000, 3.0000]
            sage: v.prod()
            6.0
        """
        cdef double s = 1
        cdef Py_ssize_t i
        for i from 0 <= i < self._length:
            s *= self._values[i]
        return s

    def mean(self):
        r"""
        Return the mean (average) of the elements of ``self``.

        OUTPUT: double

        EXAMPLES::

            sage: v = stats.TimeSeries([1,1,1,2,3]); v
            [1.0000, 1.0000, 1.0000, 2.0000, 3.0000]
            sage: v.mean()
            1.6
        """
        return self.sum() / self._length

    def pow(self, double k):
        r"""
        Return a new time series with all elements of ``self`` raised to the
        `k`-th power.

        INPUT:

        - ``k`` -- a float

        OUTPUT: a time series

        EXAMPLES::

            sage: v = stats.TimeSeries([1,1,1,2,3]); v
            [1.0000, 1.0000, 1.0000, 2.0000, 3.0000]
            sage: v.pow(2)
            [1.0000, 1.0000, 1.0000, 4.0000, 9.0000]
        """
        cdef Py_ssize_t i
        cdef TimeSeries t = new_time_series(self._length)
        for i from 0 <= i < self._length:
            t._values[i] = pow(self._values[i], k)
        return t

    def moment(self, int k):
        r"""
        Return the `k`-th moment of ``self``, which is just the
        mean of the `k`-th powers of the elements of ``self``.

        INPUT:

        - ``k`` -- positive integer

        OUTPUT: double

        EXAMPLES::

            sage: v = stats.TimeSeries([1,1,1,2,3]); v
            [1.0000, 1.0000, 1.0000, 2.0000, 3.0000]
            sage: v.moment(1)
            1.6
            sage: v.moment(2)
            3.2
        """
        if k <= 0:
            raise ValueError("k must be positive")
        if k == 1:
            return self.mean()
        cdef double s = 0
        cdef Py_ssize_t i
        for i from 0 <= i < self._length:
            s += pow(self._values[i], k)
        return s / self._length

    def central_moment(self, int k):
        r"""
        Return the `k`-th central moment of ``self``, which is just the mean
        of the `k`-th powers of the differences ``self[i] - mu``, where ``mu`` is
        the mean of ``self``.

        INPUT:

        - ``k`` -- positive integer

        OUTPUT: double

        EXAMPLES::

            sage: v = stats.TimeSeries([1,2,3])
            sage: v.central_moment(2)
            0.6666666666666666

        Note that the central moment is different from the moment
        here, since the mean is not `0`::

            sage: v.moment(2)     # doesn't subtract off mean
            4.666666666666667

        We compute the central moment directly::

            sage: mu = v.mean(); mu
            2.0
            sage: ((1-mu)^2 + (2-mu)^2 + (3-mu)^2) / 3
            0.6666666666666666
        """
        if k == 1:
            return float(0)
        mu = self.mean()
        # We could make this slightly faster by doing the add scalar
        # and moment calculation together.  But that would be nasty.
        return self.add_scalar(-mu).moment(k)

    def covariance(self, TimeSeries other):
        r"""
        Return the covariance of the time series ``self`` and ``other``.

        INPUT:

        - ``self``, ``other`` -- time series

        Whichever time series has more terms is truncated.

        EXAMPLES::

            sage: v = stats.TimeSeries([1,-2,3]); w = stats.TimeSeries([4,5,-10])
            sage: v.covariance(w)
            -11.777777777777779
        """
        cdef double mu = self.mean(), mu2 = other.mean()
        cdef double s = 0
        cdef Py_ssize_t i
        cdef Py_ssize_t n = min(self._length, other._length)
        for i from 0 <= i < n:
            s += (self._values[i] - mu)*(other._values[i] - mu2)
        return s / n
        # NOTE: There is also a formula
        #       self._values[i]*other._values[i] - mu*mu2
        # but that seems less numerically stable (?), and when tested
        # was not noticeably faster.

    def autocovariance(self, Py_ssize_t k=0):
        r"""
        Return the `k`-th autocovariance function `\gamma(k)` of ``self``.

        This is the covariance of ``self`` with ``self`` shifted by `k`, i.e.,

        .. MATH::

            \left.
            \left( \sum_{t=0}^{n-k-1} (x_t - \mu)(x_{t + k} - \mu) \right)
            \right/ n,

        where `n` is the length of ``self``.

        Note the denominator of `n`, which gives a "better" sample
        estimator.

        INPUT:

        - ``k`` -- nonnegative integer (default: 0)

        OUTPUT: float

        EXAMPLES::

            sage: v = stats.TimeSeries([13,8,15,4,4,12,11,7,14,12])
            sage: v.autocovariance(0)
            14.4
            sage: mu = v.mean(); sum([(a-mu)^2 for a in v])/len(v)
            14.4
            sage: v.autocovariance(1)
            -2.7
            sage: mu = v.mean()
            sage: sum((v[i]-mu)*(v[i+1]-mu) for i in range(len(v)-1))/len(v)
            -2.7
            sage: v.autocovariance(1)
            -2.7

        We illustrate with a random sample that an independently and
        identically distributed distribution with zero mean and
        variance `\sigma^2` has autocovariance function `\gamma(h)`
        with `\gamma(0) = \sigma^2` and `\gamma(h) = 0` for `h\neq 0`. ::

            sage: set_random_seed(0)
            sage: v = stats.TimeSeries(10^6)
            sage: v.randomize('normal', 0, 5)
            [3.3835, -2.0055, 1.7882, -2.9319, -4.6827 ...
             -5.1868, 9.2613, 0.9274, -6.2282, -8.7652]
            sage: v.autocovariance(0)
            24.95410689...
            sage: v.autocovariance(1)
            -0.00508390...
            sage: v.autocovariance(2)
            0.022056325...
            sage: v.autocovariance(3)
            -0.01902000...
        """
        cdef double mu = self.mean()
        cdef double s = 0
        cdef Py_ssize_t i
        cdef Py_ssize_t n = self._length - k
        for i from 0 <= i < n:
            s += (self._values[i] - mu)*(self._values[i+k] - mu)
        return s / self._length

    def correlation(self, TimeSeries other):
        r"""
        Return the correlation of ``self`` and ``other``, which is the
        covariance of ``self`` and ``other`` divided by the product of their
        standard deviation.

        INPUT:

        - ``self``, ``other`` -- time series

        Whichever time series has more terms is truncated.

        EXAMPLES::

            sage: v = stats.TimeSeries([1,-2,3]); w = stats.TimeSeries([4,5,-10])
            sage: v.correlation(w)
            -0.558041609...
            sage: v.covariance(w)/(v.standard_deviation() * w.standard_deviation())
            -0.558041609...
        """
        return self.covariance(other) / (self.standard_deviation() * other.standard_deviation())

    def autocorrelation(self, Py_ssize_t k=1):
        r"""
        Return the `k`-th sample autocorrelation of this time series
        `x_i`.

        Let `\mu` be the sample mean.  Then the sample autocorrelation
        function is

        .. MATH::

            \frac{\sum_{t=0}^{n-k-1} (x_t - \mu)(x_{t+k} - \mu) }
                 {\sum_{t=0}^{n-1}   (x_t - \mu)^2}.

        Note that the variance must be nonzero or you will get a
        :exc:`ZeroDivisionError`.

        INPUT:

        - ``k`` -- nonnegative integer (default: 1)

        OUTPUT: a time series

        EXAMPLES::

            sage: v = stats.TimeSeries([13,8,15,4,4,12,11,7,14,12])
            sage: v.autocorrelation()
            -0.1875
            sage: v.autocorrelation(1)
            -0.1875
            sage: v.autocorrelation(0)
            1.0
            sage: v.autocorrelation(2)
            -0.20138888888888887
            sage: v.autocorrelation(3)
            0.18055555555555555

            sage: stats.TimeSeries([1..1000]).autocorrelation()
            0.997
        """
        return self.autocovariance(k) / self.variance(bias=True)

    def variance(self, bias=False):
        r"""
        Return the variance of the elements of ``self``, which is the mean
        of the squares of the differences from the mean.

        INPUT:

        - ``bias`` -- boolean (default: ``False``); if ``False``, divide by
          ``self.length() - 1`` instead of ``self.length()`` to give a less
          biased estimator for the variance.

        OUTPUT: double

        EXAMPLES::

            sage: v = stats.TimeSeries([1,1,1,2,3]); v
            [1.0000, 1.0000, 1.0000, 2.0000, 3.0000]
            sage: v.variance()
            0.8
            sage: v.variance(bias=True)
            0.64

        TESTS::

            sage: stats.TimeSeries([1]).variance()
            0.0
            sage: stats.TimeSeries([]).variance()
            0.0
        """
        if self._length <= 1:
            return float(0)
        cdef double mu = self.mean()
        cdef double s = 0
        cdef double a
        cdef Py_ssize_t i
        for i from 0 <= i < self._length:
            a = self._values[i] - mu
            s += a * a
        if bias:
            return s / self._length
        else:
            return s / (self._length - 1)

    def standard_deviation(self, bias=False):
        r"""
        Return the standard deviation of the entries of ``self``.

        INPUT:

        - ``bias`` -- boolean (default: ``False``); if ``False``, divide by
          ``self.length() - 1`` instead of ``self.length()`` to give a less
          biased estimator for the variance.

        OUTPUT: double

        EXAMPLES::

            sage: v = stats.TimeSeries([1,1,1,2,3]); v
            [1.0000, 1.0000, 1.0000, 2.0000, 3.0000]
            sage: v.standard_deviation()
            0.8944271909...
            sage: v.standard_deviation(bias=True)
            0.8

        TESTS::

            sage: stats.TimeSeries([1]).standard_deviation()
            0.0
            sage: stats.TimeSeries([]).standard_deviation()
            0.0
        """
        return sqrt(self.variance(bias=bias))

    def range_statistic(self, b=None):
        r"""
        Return the rescaled range statistic `R/S` of ``self``, which is
        defined as follows (see Hurst 1951).  If the optional
        parameter ``b`` is given, return the average of `R/S` range
        statistics of disjoint blocks of size ``b``.

        Let `\sigma` be the standard deviation of the sequence of
        differences of ``self``, and let `Y_k` be the `k`-th term of ``self``.
        Let `n` be the number of terms of ``self``, and set
        `Z_k = Y_k - ((k+1)/n) \cdot Y_n`. Then

        .. MATH::

            R/S = \big( \max(Z_k) - \min(Z_k) \big) / \sigma

        where the max and min are over all `Z_k`.
        Basically replacing `Y_k` by `Z_k` allows us to measure
        the difference from the line from the origin to `(n,Y_n)`.

        INPUT:

        - ``self`` -- a time series  (*not* the series of differences)

        - ``b`` -- integer (default: ``None``); if given instead divide the
          input time series up into ``j = floor(n/b)`` disjoint
          blocks, compute the `R/S` statistic for each block,
          and return the average of those `R/S` statistics.

        OUTPUT: float

        EXAMPLES:

        Notice that if we make a Brownian motion random walk, there
        is no difference if we change the standard deviation. ::

            sage: set_random_seed(0); stats.TimeSeries(10^6).randomize('normal').sums().range_statistic()
            1897.8392602...
            sage: set_random_seed(0); stats.TimeSeries(10^6).randomize('normal',0,100).sums().range_statistic()
            1897.8392602...
        """
        cdef Py_ssize_t j, k, n = self._length

        if b is not None:
            j = n // b
            return sum([self[k*b:(k+1)*b].add_scalar(-self[k*b]).range_statistic() for k in range(j)]) / j

        cdef double Yn = self._values[n - 1]
        cdef TimeSeries Z = self.__copy__()
        for k from 0 <= k < n:
            Z._values[k] -= (k+1)*Yn / n

        sigma = self.diffs().standard_deviation()

        return (Z.max() - Z.min()) / sigma

    def hurst_exponent(self):
        r"""
        Return an estimate of the Hurst exponent of this time series.

        We use the algorithm from pages 61 -- 63 of [Peteres, Fractal
        Market Analysis (1994); see Google Books].

        We define the Hurst exponent of a constant time series to be 1.

        EXAMPLES:

        The Hurst exponent of Brownian motion is 1/2.  We approximate
        it with some specific samples.  Note that the estimator is
        biased and slightly overestimates. ::

            sage: set_random_seed(0)
            sage: bm = stats.TimeSeries(10^5).randomize('normal').sums(); bm
            [0.6767, 0.2756, 0.6332, 0.0469, -0.8897 ...
             152.2437, 151.5327, 152.7629, 152.9169, 152.9084]
            sage: bm.hurst_exponent()
            0.527450972...
        """
        # We use disjoint blocks of size 8, 16, 32, ....
        cdef Py_ssize_t k = 8
        if self._length <= 32:   # small data, estimate of Hurst will suck but...
            k = 4
        v0 = []
        v1 = []
        while k <= self._length:
            try:
                v1.append(log(self.range_statistic(k)))
                v0.append(log(k))
            except ZeroDivisionError:   # 0 standard deviation
                pass
            k *= 2
        if not v0:
            return float(1)
        if len(v0) == 1:
            return v1[0]/v0[0]
        import numpy
        if int(numpy.version.short_version[0]) > 1:
            numpy.set_printoptions(legacy="1.25")
        coeffs = numpy.polyfit(v0,v1,1)
        return coeffs[0]

    def min(self, bint index=False):
        r"""
        Return the smallest value in this time series.

        If this series has length 0, we raise a :exc:`ValueError`.

        INPUT:

        - ``index`` -- boolean (default: ``False``); if ``True``, also return
          index of minimal entry

        OUTPUT:

        - float; smallest value

        - integer; index of smallest value (only returned if ``index=True``)

        EXAMPLES::

            sage: v = stats.TimeSeries([1,-4,3,-2.5,-4])
            sage: v.min()
            -4.0
            sage: v.min(index=True)
            (-4.0, 1)
        """
        if self._length == 0:
            raise ValueError("min() arg is an empty sequence")
        cdef Py_ssize_t i, j
        cdef double s = self._values[0]
        j = 0
        for i from 1 <= i < self._length:
            if self._values[i] < s:
                s = self._values[i]
                j = i
        if index:
            return s, j
        else:
            return s

    def max(self, bint index=False):
        r"""
        Return the largest value in this time series. If this series
        has length 0 we raise a :exc:`ValueError`.

        INPUT:

        - ``index`` -- boolean (default: ``False``); if ``True``, also return
          index of maximum entry

        OUTPUT:

        - float; largest value

        - integer; index of largest value (only returned if ``index=True``)

        EXAMPLES::

            sage: v = stats.TimeSeries([1,-4,3,-2.5,-4,3])
            sage: v.max()
            3.0
            sage: v.max(index=True)
            (3.0, 2)
        """
        if self._length == 0:
            raise ValueError("max() arg is an empty sequence")
        cdef Py_ssize_t i, j = 0
        cdef double s = self._values[0]
        for i from 1 <= i < self._length:
            if self._values[i] > s:
                s = self._values[i]
                j = i
        if index:
            return s, j
        else:
            return s

    def clip_remove(self, min=None, max=None):
        r"""
        Return new time series obtained from ``self`` by removing all
        values that are less than or equal to a certain minimum value
        or greater than or equal to a certain maximum.

        INPUT:

        - ``min`` -- (default: ``None``) ``None`` or double

        - ``max`` -- (default: ``None``) ``None`` or double

        OUTPUT: a time series

        EXAMPLES::

            sage: v = stats.TimeSeries([1..10])
            sage: v.clip_remove(3,7)
            [3.0000, 4.0000, 5.0000, 6.0000, 7.0000]
            sage: v.clip_remove(3)
            [3.0000, 4.0000, 5.0000, 6.0000, 7.0000, 8.0000, 9.0000, 10.0000]
            sage: v.clip_remove(max=7)
            [1.0000, 2.0000, 3.0000, 4.0000, 5.0000, 6.0000, 7.0000]
        """
        cdef Py_ssize_t i, j
        cdef TimeSeries t
        cdef double mn, mx
        cdef double x
        if min is None and max is None:
            return self.copy()
        # This code is ugly but I think as fast as possible.
        if min is None and max is not None:
            # Return everything <= max
            j = 0
            mx = max
            for i from 0 <= i < self._length:
                if self._values[i] <= mx:
                    j += 1

            t = TimeSeries(j)
            j = 0
            for i from 0 <= i < self._length:
                if self._values[i] <= mx:
                    t._values[j] = self._values[i]
                    j += 1
        elif max is None and min is not None:
            # Return everything >= min
            j = 0
            mn = min
            for i from 0 <= i < self._length:
                if self._values[i] >= mn:
                    j += 1
            t = TimeSeries(j)
            j = 0
            for i from 0 <= i < self._length:
                if self._values[i] >= mn:
                    t._values[j] = self._values[i]
                    j += 1
        else:
            # Return everything between min and max
            j = 0
            mn = min
            mx = max
            for i from 0 <= i < self._length:
                x = self._values[i]
                if x >= mn and x <= mx:
                    j += 1
            t = TimeSeries(j)
            j = 0
            for i from 0 <= i < self._length:
                x = self._values[i]
                if x >= mn and x <= mx:
                    t._values[j] = x
                    j += 1
        return t

    def histogram(self, Py_ssize_t bins=50, bint normalize=False):
        r"""
        Return the frequency histogram of the values in
        this time series divided up into the given
        number of bins.

        INPUT:

        - ``bins`` -- positive integer (default: 50)

        - ``normalize`` -- boolean (default: ``False``); whether to normalize so the
          total area in the bars of the histogram is 1

        OUTPUT:

        - counts -- list of counts of numbers of elements in each bin

        - endpoints -- list of 2-tuples (a,b) that give the endpoints of the bins

        EXAMPLES::

            sage: v = stats.TimeSeries([5,4,1.3,2,8,10,3,-5]); v
            [5.0000, 4.0000, 1.3000, 2.0000, 8.0000, 10.0000, 3.0000, -5.0000]
            sage: v.histogram(3)
            ([1, 4, 3], [(-5.0, 0.0), (0.0, 5.0), (5.0, 10.0)])
        """
        if bins <= 0:
            raise ValueError("bins must be positive")

        cdef double mn = self.min(), mx = self.max()
        cdef double r = mx - mn, step = r/bins
        cdef Py_ssize_t j

        if r == 0:
            raise ValueError("bins have 0 width")

        v = [(mn + j*step, mn + (j+1)*step) for j in range(bins)]
        if self._length == 0:
            return [], v

        if step == 0:
            counts = [0]*bins
            counts[0] = [self._length]
            return counts, v

        cdef Py_ssize_t i
        cdef Py_ssize_t* cnts = <Py_ssize_t*>sig_malloc(sizeof(Py_ssize_t)*bins)
        if cnts == NULL:
            raise MemoryError
        for i from 0 <= i < bins:
            cnts[i] = 0

        for i from 0 <= i < self._length:
            j = int((self._values[i] - mn)/step)
            if j >= bins:
                j = bins-1
            cnts[j] += 1

        b = 1.0/(self._length * step)
        if normalize:
            counts = [cnts[i]*b for i in range(bins)]
        else:
            counts = [cnts[i] for i in range(bins)]
        sig_free(cnts)

        return counts, v

    def plot_histogram(self, bins=50, normalize=True, **kwds):
        r"""
        Return histogram plot of this time series with given number of bins.

        INPUT:

        - ``bins`` -- positive integer (default: 50)

        - ``normalize`` -- boolean (default: ``True``); whether to normalize so the
          total area in the bars of the histogram is 1

        OUTPUT: a histogram plot

        EXAMPLES::

            sage: v = stats.TimeSeries([1..50])
            sage: v.plot_histogram(bins=10)                                             # needs sage.plot
            Graphics object consisting of 10 graphics primitives

        ::

            sage: v.plot_histogram(bins=3,normalize=False,aspect_ratio=1)               # needs sage.plot
            Graphics object consisting of 3 graphics primitives
        """
        from sage.plot.all import polygon
        counts, intervals = self.histogram(bins, normalize=normalize)
        s = 0
        kwds.setdefault('aspect_ratio','automatic')
        for i, (x0,x1) in enumerate(intervals):
            s += polygon([(x0,0), (x0,counts[i]), (x1,counts[i]), (x1,0)], **kwds)
        if intervals:
            s.axes_range(ymin=0, ymax=max(counts), xmin=intervals[0][0], xmax=intervals[-1][1])
        return s

    def plot_candlestick(self, int bins=30):
        r"""
        Return a candlestick plot of this time series with the given number
        of bins.

        A candlestick plot is a style of bar-chart used to view open, high,
        low, and close stock data. At each bin, the line represents the
        high / low range. The bar represents the open / close range. The
        interval is colored blue if the open for that bin is less than the
        close. If the close is less than the open, then that bin is colored
        red instead.

        INPUT:

        - ``bins`` -- positive integer (default: 30); the number of bins
          or candles

        OUTPUT: a candlestick plot

        EXAMPLES:

        Here we look at the candlestick plot for Brownian motion::

            sage: v = stats.TimeSeries(1000).randomize()
            sage: v.plot_candlestick(bins=20)                                           # needs sage.plot
            Graphics object consisting of 40 graphics primitives
        """
        from sage.plot.all import line, polygon, Graphics

        cdef TimeSeries t = new_time_series(self._length)
        cdef TimeSeries s
        cdef int i, j, n
        bin_size = int(t._length/bins)
        rng      = t._length - bin_size

        memcpy(t._values, self._values, sizeof(double)*self._length)
        p = Graphics()

        for i from 0 <= i < bins:
            n = i*bin_size
            s = new_time_series(bin_size)
            for j from 0 <= j < bin_size:
                s._values[j] = t._values[n + j]
            low   = s.min()
            high  = s.max()
            open  = s._values[0]
            close = s._values[bin_size-1]
            left  = n + bin_size/3
            mid   = n + bin_size/2
            right = n + 2*bin_size/3

            rgbcolor =  'blue' if open < close else 'red'

            p += line([(mid, low), (mid, high)], rgbcolor=rgbcolor)
            p += polygon([(left, open), (right, open), (right, close), (left, close)], rgbcolor=rgbcolor)

        return p

    def numpy(self, copy=True):
        r"""
        Return a NumPy version of this time series.

        .. NOTE::

            If copy is ``False``, return a NumPy 1-D array reference to
            exactly the same block of memory as this time series.  This is
            very, very fast and makes it easy to quickly use all
            NumPy/SciPy functionality on ``self``.  However, it is dangerous
            because when this time series goes out of scope and is garbage
            collected, the corresponding NumPy reference object will point
            to garbage.

        INPUT:

        - ``copy`` -- boolean (default: ``True``)

        OUTPUT: a numpy 1-D array

        EXAMPLES::

            sage: v = stats.TimeSeries([1,-3,4.5,-2])
            sage: w = v.numpy(copy=False); w
            array([ 1. , -3. ,  4.5, -2. ])
            sage: type(w)
            <... 'numpy.ndarray'>
            sage: w.shape
            (4,)

        Notice that changing ``w`` also changes ``v`` too! ::

            sage: w[0] = 20
            sage: w
            array([20. , -3. ,  4.5, -2. ])
            sage: v
            [20.0000, -3.0000, 4.5000, -2.0000]

        If you want a separate copy do not give the ``copy=False`` option. ::

            sage: z = v.numpy(); z
            array([20. , -3. ,  4.5, -2. ])
            sage: z[0] = -10
            sage: v
            [20.0000, -3.0000, 4.5000, -2.0000]
        """
        cnumpy.import_array()  # This must be called before using the numpy C/api or you will get segfault
        cdef cnumpy.npy_intp dims[1]
        dims[0] = self._length
        cdef cnumpy.ndarray n = cnumpy.PyArray_SimpleNewFromData(1, dims, cnumpy.NPY_DOUBLE, self._values)
        if copy:
            return n.copy()
        else:
#            Py_INCREF(self)
#            n.base = self
            return n

    def randomize(self, distribution='uniform', loc=0, scale=1, **kwds):
        r"""
        Randomize the entries in this time series, and return a reference
        to ``self``.  Thus this function both changes ``self`` in place, and
        returns a copy of ``self``, for convenience.

        INPUT:

        - ``distribution`` -- (default: ``'uniform'``) supported values are:

          - ``'uniform'`` -- from ``loc`` to ``loc + scale``

          - ``'normal'`` -- mean ``loc`` and standard deviation ``scale``

          - ``'semicircle'`` -- with center at ``loc`` (``scale`` is ignored)

          - ``'lognormal'`` -- mean ``loc`` and standard deviation ``scale``

        - ``loc`` -- float (default: 0)

        - ``scale`` -- float (default: 1)

        .. NOTE::

            All random numbers are generated using algorithms that
            build on the high quality GMP random number function
            gmp_urandomb_ui.  Thus this function fully respects the Sage
            ``set_random_state`` command.  It's not quite as fast as the C
            library random number generator, but is of course much better
            quality, and is platform independent.

        EXAMPLES:

        We generate 5 uniform random numbers in the interval [0,1]::

            sage: set_random_seed(0)
            sage: stats.TimeSeries(5).randomize()
            [0.8685, 0.2816, 0.0229, 0.1456, 0.7314]

        We generate 5 uniform random numbers from 5 to `5+2=7`::

            sage: set_random_seed(0)
            sage: stats.TimeSeries(10).randomize('uniform', 5, 2)
            [6.7369, 5.5632, 5.0459, 5.2911, 6.4628, 5.2412, 5.2010, 5.2761, 5.5813, 5.5439]

        We generate 5 normal random values with mean 0 and variance 1. ::

            sage: set_random_seed(0)
            sage: stats.TimeSeries(5).randomize('normal')
            [0.6767, -0.4011, 0.3576, -0.5864, -0.9365]

        We generate 10 normal random values with mean 5 and variance 2. ::

            sage: set_random_seed(0)
            sage: stats.TimeSeries(10).randomize('normal', 5, 2)
            [6.3534, 4.1978, 5.7153, 3.8273, 3.1269, 2.9598, 3.7484, 6.7472, 3.8986, 4.6271]

        We generate 5 values using the semicircle distribution. ::

            sage: set_random_seed(0)
            sage: stats.TimeSeries(5).randomize('semicircle')
            [0.7369, -0.9541, 0.4628, -0.7990, -0.4187]

        We generate 1 million normal random values and create a frequency
        histogram. ::

            sage: set_random_seed(0)
            sage: a = stats.TimeSeries(10^6).randomize('normal')
            sage: a.histogram(10)[0]
            [36, 1148, 19604, 130699, 340054, 347870, 137953, 21290, 1311, 35]

        We take the above values, and compute the proportion that lie within
        1, 2, 3, 5, and 6 standard deviations of the mean (0)::

            sage: s = a.standard_deviation()
            sage: len(a.clip_remove(-s,s))/float(len(a))
            0.683094
            sage: len(a.clip_remove(-2*s,2*s))/float(len(a))
            0.954559
            sage: len(a.clip_remove(-3*s,3*s))/float(len(a))
            0.997228
            sage: len(a.clip_remove(-5*s,5*s))/float(len(a))
            0.999998

        There were no "six sigma events"::

            sage: len(a.clip_remove(-6*s,6*s))/float(len(a))
            1.0
        """
        if distribution == 'uniform':
            self._randomize_uniform(loc, loc + scale)
        elif distribution == 'normal':
            self._randomize_normal(loc, scale)
        elif distribution == 'semicircle':
            self._randomize_semicircle(loc)
        elif distribution == 'lognormal':
            self._randomize_lognormal(loc, scale)
        else:
            raise NotImplementedError
        return self

    def _randomize_uniform(self, double left, double right):
        r"""
        Generates a uniform random distribution of doubles between ``left`` and
        ``right`` and stores values in place.

        INPUT:

        - ``left`` -- left bound on random distribution

        - ``right`` -- right bound on random distribution

        EXAMPLES:

        We generate 5 values distributed with respect to the uniform
        distribution over the interval [0,1]::

            sage: v = stats.TimeSeries(5)
            sage: set_random_seed(0)
            sage: v.randomize('uniform')
            [0.8685, 0.2816, 0.0229, 0.1456, 0.7314]

        We now test that the mean is indeed 0.5::

            sage: v = stats.TimeSeries(10^6)
            sage: set_random_seed(0)
            sage: v.randomize('uniform').mean()
            0.50069085...
        """
        if left >= right:
            raise ValueError("left must be less than right")

        cdef randstate rstate = current_randstate()
        cdef Py_ssize_t k
        cdef double d = right - left
        for k from 0 <= k < self._length:
            self._values[k] = rstate.c_rand_double() * d + left

    def _randomize_normal(self, double m, double s):
        r"""
        Generates a normal random distribution of doubles with mean ``m`` and
        standard deviation ``s`` and stores values in place. Uses the
        Box-Muller algorithm.

        INPUT:

        - ``m`` -- mean

        - ``s`` -- standard deviation

        EXAMPLES:

        We generate 5 values distributed with respect to the normal
        distribution with mean 0 and standard deviation 1::

            sage: set_random_seed(0)
            sage: v = stats.TimeSeries(5)
            sage: v.randomize('normal')
            [0.6767, -0.4011, 0.3576, -0.5864, -0.9365]

        We now test that the mean is indeed 0::

            sage: set_random_seed(0)
            sage: v = stats.TimeSeries(10^6)
            sage: v.randomize('normal').mean()
            6.2705472723...

        The same test with mean equal to 2 and standard deviation equal
        to 5::

            sage: set_random_seed(0)
            sage: v = stats.TimeSeries(10^6)
            sage: v.randomize('normal', 2, 5).mean()
            2.0003135273...
        """
        # Ported from http://users.tkk.fi/~nbeijar/soft/terrain/source_o2/boxmuller.c
        # This the box muller algorithm.
        cdef randstate rstate = current_randstate()

        cdef double x1, x2, w, y1, y2
        cdef Py_ssize_t k
        for k from 0 <= k < self._length:
            while True:
                x1 = 2*rstate.c_rand_double() - 1
                x2 = 2*rstate.c_rand_double() - 1
                w = x1*x1 + x2*x2
                if w < 1: break
            w = sqrt( (-2*log(w))/w )
            y1 = x1 * w
            y2 = x2 * w
            self._values[k] = m + y1*s
            k += 1
            if k < self._length:
                self._values[k] = m + y2*s

    def _randomize_semicircle(self, double center):
        r"""
        Generates a semicircle random distribution of doubles about center
        and stores values in place. Uses the acceptance-rejection method.

        INPUT:

        - ``center`` -- the center of the semicircle distribution

        EXAMPLES:

        We generate 5 values distributed with respect to the semicircle
        distribution located at center::

            sage: v = stats.TimeSeries(5)
            sage: set_random_seed(0)
            sage: v.randomize('semicircle')
            [0.7369, -0.9541, 0.4628, -0.7990, -0.4187]

        We now test that the mean is indeed the center::

            sage: v = stats.TimeSeries(10^6)
            sage: set_random_seed(0)
            sage: v.randomize('semicircle').mean()
            0.0007207497...

        The same test with center equal to 2::

            sage: v = stats.TimeSeries(10^6)
            sage: set_random_seed(0)
            sage: v.randomize('semicircle', 2).mean()
            2.0007207497...
        """
        cdef Py_ssize_t k
        cdef double x, y, s, d = 2, left = center - 1, z
        cdef randstate rstate = current_randstate()
        z = d*d
        s = 1.5707963267948966192  # pi/2
        for k from 0 <= k < self._length:
            while True:
                x = rstate.c_rand_double() * d - 1
                y = rstate.c_rand_double() * s
                if y*y + x*x < 1:
                    break
            self._values[k] = x + center

    def _randomize_lognormal(self, double m, double s):
        r"""
        Generates a log-normal random distribution of doubles with mean ``m``
        and standard deviation ``s``. Uses Box-Muller algorithm and the
        identity: if `Y` is a random variable with normal distribution then
        `X = \exp(Y)` is a random variable with log-normal distribution.

        INPUT:

        - ``m`` -- mean

        - ``s`` -- standard deviation

        EXAMPLES:

        We generate 5 values distributed with respect to the lognormal
        distribution with mean 0 and standard deviation 1::

            sage: set_random_seed(0)
            sage: v = stats.TimeSeries(5)
            sage: v.randomize('lognormal')
            [1.9674, 0.6696, 1.4299, 0.5563, 0.3920]

        We now test that the mean is indeed `\sqrt{e}`::

            sage: set_random_seed(0)
            sage: v = stats.TimeSeries(10^6)
            sage: v.randomize('lognormal').mean()
            1.647351973...
            sage: exp(0.5)                                                              # needs sage.symbolic
            1.648721270...

        A log-normal distribution can be simply thought of as the logarithm
        of a normally distributed data set. We test that here by generating
        5 values distributed with respect to the normal distribution with mean
        0 and standard deviation 1::

            sage: set_random_seed(0)
            sage: w = stats.TimeSeries(5)
            sage: w.randomize('normal')
            [0.6767, -0.4011, 0.3576, -0.5864, -0.9365]
            sage: exp(w)
            [1.9674, 0.6696, 1.4299, 0.5563, 0.3920]
        """
        # Ported from http://users.tkk.fi/~nbeijar/soft/terrain/source_o2/boxmuller.c
        # This the box muller algorithm.
        cdef randstate rstate = current_randstate()
        cdef double x1, x2, w, y1, y2
        cdef Py_ssize_t k
        for k from 0 <= k < self._length:
            while True:
                x1 = 2*rstate.c_rand_double() - 1
                x2 = 2*rstate.c_rand_double() - 1
                w = x1*x1 + x2*x2
                if w < 1: break
            w = sqrt( (-2*log(w))/w )
            y1 = x1 * w
            y2 = x2 * w
            self._values[k] = exp(m + y1*s)
            k += 1
            if k < self._length:
                self._values[k] = exp(m + y2*s)

    def fft(self, bint overwrite=False):
        r"""
        Return the real discrete fast Fourier transform of ``self``, as a
        real time series:

        .. MATH::

            [y(0),\Re(y(1)),\Im(y(1)),\dots,\Re(y(n/2))]  \text{ if $n$ is even}

            [y(0),\Re(y(1)),\Im(y(1)),\dots,\Re(y(n/2)),\Im(y(n/2))] \text{ if $n$ is odd}

        where

        .. MATH::

            y(j) = \sum_{k=0}^{n-1} x[k] \cdot \exp(-\sqrt{-1} \cdot jk \cdot 2\pi/n)

        for `j = 0,\dots,n-1`.  Note that `y(-j) = y(n-j)`.

        EXAMPLES::

            sage: v = stats.TimeSeries([1..9]); v
            [1.0000, 2.0000, 3.0000, 4.0000, 5.0000, 6.0000, 7.0000, 8.0000, 9.0000]
            sage: w = v.fft(); w
            [45.0000, -4.5000, 12.3636, -4.5000, 5.3629, -4.5000, 2.5981, -4.5000, 0.7935]

        We get just the series of real parts of ::

            sage: stats.TimeSeries([w[0]]) + w[1:].scale_time(2)
            [45.0000, -4.5000, -4.5000, -4.5000, -4.5000]

        An example with an even number of terms::

            sage: v = stats.TimeSeries([1..10]); v
            [1.0000, 2.0000, 3.0000, 4.0000, 5.0000, 6.0000, 7.0000, 8.0000, 9.0000, 10.0000]
            sage: w = v.fft(); w
            [55.0000, -5.0000, 15.3884, -5.0000, 6.8819, -5.0000, 3.6327, -5.0000, 1.6246, -5.0000]
            sage: v.fft().ifft()
            [1.0000, 2.0000, 3.0000, 4.0000, 5.0000, 6.0000, 7.0000, 8.0000, 9.0000, 10.0000]
        """
        import scipy.fftpack
        if overwrite:
            y = self
        else:
            y = self.__copy__()
        w = y.numpy(copy=False)
        scipy.fftpack.rfft(w, overwrite_x=1)
        return y

    def ifft(self, bint overwrite=False):
        r"""
        Return the real discrete inverse fast Fourier transform of
        ``self``, which is also a real time series.

        This is the inverse of ``fft()``.

        The returned real array contains

        .. MATH::

            [y(0),y(1),\dots,y(n-1)]

        where for `n` is even

        .. MATH::

            y(j)
            =
            1/n \left(
            \sum_{k=1}^{n/2-1}
            (x[2k-1]+\sqrt{-1} \cdot x[2k])
            \cdot \exp(\sqrt{-1} \cdot jk \cdot 2pi/n)
            + c.c. + x[0] + (-1)^j x[n-1]
            \right)

        and for `n` is odd

        .. MATH::

            y(j)
            =
            1/n \left(
            \sum_{k=1}^{(n-1)/2}
            (x[2k-1]+\sqrt{-1} \cdot x[2k])
            \cdot \exp(\sqrt{-1} \cdot jk \cdot 2pi/n)
            + c.c. + x[0]
            \right)

        where `c.c.` denotes complex conjugate of preceding expression.

        EXAMPLES::

            sage: v = stats.TimeSeries([1..10]); v
            [1.0000, 2.0000, 3.0000, 4.0000, 5.0000, 6.0000, 7.0000, 8.0000, 9.0000, 10.0000]
            sage: v.ifft()
            [5.1000, -5.6876, 1.4764, -1.0774, 0.4249, -0.1000, -0.2249, 0.6663, -1.2764, 1.6988]
            sage: v.ifft().fft()
            [1.0000, 2.0000, 3.0000, 4.0000, 5.0000, 6.0000, 7.0000, 8.0000, 9.0000, 10.0000]
        """
        import scipy.fftpack
        if overwrite:
            y = self
        else:
            y = self.__copy__()
        w = y.numpy(copy=False)
        scipy.fftpack.irfft(w, overwrite_x=1)
        return y


cdef new_time_series(Py_ssize_t length):
    r"""
    Return a new uninitialized time series of the given length.
    The entries of the time series are garbage.

    INPUT:

    - ``length`` -- integer

    OUTPUT: a time series

    EXAMPLES:

    This uses ``new_time_series`` implicitly::

        sage: v = stats.TimeSeries([1,-3,4.5,-2])
        sage: v.__copy__()
        [1.0000, -3.0000, 4.5000, -2.0000]
    """
    if length < 0:
        raise ValueError("length must be nonnegative")
    cdef TimeSeries t = TimeSeries.__new__(TimeSeries)
    t._length = length
    t._values = <double*> sig_malloc(sizeof(double)*length)
    return t


@cython.binding(True)
def unpickle_time_series_v1(bytes v, Py_ssize_t n):
    r"""
    Version 1 unpickle method.

    INPUT:

    - ``v`` -- a raw char buffer

    EXAMPLES::

        sage: v = stats.TimeSeries([1,2,3])
        sage: s = v.__reduce__()[1][0]
        sage: type(s)
        <class 'bytes'>
        sage: sage.stats.time_series.unpickle_time_series_v1(s,3)
        [1.0000, 2.0000, 3.0000]
        sage: sage.stats.time_series.unpickle_time_series_v1(s+s,6)
        [1.0000, 2.0000, 3.0000, 1.0000, 2.0000, 3.0000]
        sage: sage.stats.time_series.unpickle_time_series_v1(b'',0)
        []
    """
    cdef TimeSeries t = new_time_series(n)
    memcpy(t._values, PyBytes_AsString(v), n * sizeof(double))
    return t


register_unpickle_override('sage.finance.time_series', 'unpickle_time_series_v1', unpickle_time_series_v1)


def autoregressive_fit(acvs):
    r"""
    Given a sequence of lagged autocovariances of length `M` produce
    `a_1,\dots,a_p` so that the first `M` autocovariance coefficients
    of the autoregressive processes `X_t=a_1X_{t_1}+\cdots+a_pX_{t-p}+Z_t`
    are the same as the input sequence.

    The function works by solving the Yule-Walker equations
    `\Gamma a =\gamma`, where `\gamma=(\gamma(1),\dots,\gamma(M))`,
    `a=(a_1,\dots,a_M)`, with `\gamma(i)` the autocovariance of lag `i`
    and `\Gamma_{ij}=\gamma(i-j)`.
    """
    cdef TimeSeries c
    cdef Py_ssize_t i

    M = len(acvs)-1

    if M <= 0:
        raise ValueError("M must be positive")

    if not isinstance(acvs, TimeSeries):
        c = TimeSeries(acvs)
    else:
        c = acvs

    # Also create the numpy vector v with entries c(1), ..., c(M).
    v = c[1:].numpy()

    # 2. Create the autocovariances numpy 2d array A whose i,j entry
    # is c(|i-j|).
    import numpy
    A = numpy.array([[c[abs(j-k)] for j in range(M)] for k in range(M)])

    # 3. Solve the equation A * x = v
    x = numpy.linalg.solve(A, v)

    # 4. The entries of x give the linear filter coefficients.
    return TimeSeries(x)
