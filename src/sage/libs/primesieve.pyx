# sage.doctest: needs primesieve
r"""
Cython wrapper for the primesieve library

The `primesieve library <https://github.com/kimwalisch/primesieve>`_
generates the primes below `2^{64}` using an optimized segmented sieve
of Eratosthenes.

The two entry points of this module are :func:`prime_range`, which
returns the list of primes in an interval, and :class:`prime_iterator`,
which iterates lazily over primes. They are used behind the scenes by
:func:`sage.rings.fast_arith.prime_range`, :func:`sage.arith.misc.primes`
and :class:`sage.sets.primes.Primes`; most users should call those
functions instead.

EXAMPLES::

    sage: from sage.libs.primesieve import prime_range, prime_iterator
    sage: prime_range(20)
    [2, 3, 5, 7, 11, 13, 17, 19]
    sage: it = prime_iterator(100)
    sage: next(it), next(it), next(it)
    (101, 103, 107)

AUTHORS:

- Chenxin Zhong (2026-07): initial version
"""

# ****************************************************************************
#       Copyright (C) 2026 Chenxin Zhong <chenxin.zhong@outlook.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.signals cimport sig_on, sig_off

from libc.errno cimport errno, EDOM
from libc.stdint cimport uint64_t, UINT64_MAX

from cpython.long cimport PyLong_FromUnsignedLongLong

from sage.cpython.string cimport char_to_str
from sage.ext.stdsage cimport PY_NEW
from sage.libs.gmp.mpz cimport mpz_import
from sage.rings.integer cimport Integer


cdef extern from "primesieve.h":
    enum:
        UINT64_PRIMES
    void* primesieve_generate_primes(uint64_t start, uint64_t stop,
                                     size_t* size, int type) nogil
    void primesieve_free(void* primes)
    const char* primesieve_version()

cdef extern from "primesieve/iterator.h":
    ctypedef struct primesieve_iterator:
        size_t i
        size_t size
        uint64_t start
        uint64_t* primes
        int is_error
    void primesieve_init(primesieve_iterator* it)
    void primesieve_free_iterator(primesieve_iterator* it)
    void primesieve_jump_to(primesieve_iterator* it, uint64_t start,
                            uint64_t stop_hint)
    void primesieve_generate_next_primes(primesieve_iterator* it) nogil


# The largest prime below 2^64. Asking primesieve for primes beyond it
# makes it print an error message to stderr, so we detect that case
# ourselves in prime_iterator.__next__.
cdef uint64_t largest_uint64_prime = UINT64_MAX - 58


cdef inline Integer _new_integer(uint64_t p):
    # mpz_import is used instead of mpz_set_ui for portability to
    # platforms where unsigned long is only 32 bits (e.g. Windows).
    cdef Integer z = PY_NEW(Integer)
    mpz_import(z.value, 1, 1, sizeof(uint64_t), 0, 0, &p)
    return z


def version():
    r"""
    Return the version of the primesieve library.

    EXAMPLES::

        sage: from sage.libs.primesieve import version
        sage: version()  # random
        '12.13'
    """
    return char_to_str(primesieve_version())


cpdef list prime_range(start, stop=None, bint py_ints=False):
    r"""
    Return the list of primes in the interval ``[start, stop)``.

    If ``stop`` is omitted, return the primes up to ``start``, that is,
    the primes in the interval ``[2, start)``.

    INPUT:

    - ``start`` -- integer; lower bound

    - ``stop`` -- integer (default: ``None``); upper bound (excluded).
      Bounds beyond `2^{64}` are not supported.

    - ``py_ints`` -- boolean (default: ``False``); whether to return
      Python ints instead of Sage Integers (faster)

    EXAMPLES::

        sage: from sage.libs.primesieve import prime_range
        sage: prime_range(20)
        [2, 3, 5, 7, 11, 13, 17, 19]
        sage: prime_range(10, 30)
        [11, 13, 17, 19, 23, 29]
        sage: type(prime_range(10)[0])
        <class 'sage.rings.integer.Integer'>
        sage: type(prime_range(10, py_ints=True)[0])
        <class 'int'>

    TESTS::

        sage: prime_range(0)
        []
        sage: prime_range(2, 2)
        []
        sage: prime_range(3, 2)
        []
        sage: prime_range(-10, 3)
        [2]
        sage: prime_range(2, 3)
        [2]
        sage: prime_range(2**64, 2**64 + 1)
        Traceback (most recent call last):
        ...
        ValueError: primesieve only generates primes smaller than 2^64
        sage: prime_range(2**64 - 100, 2**64)  # long time
        [18446744073709551521, 18446744073709551533, 18446744073709551557]
    """
    if stop is None:
        stop = start
        start = 2
    elif start < 0:
        start = 0

    if stop <= start or stop <= 2:
        return []
    if stop - 1 > UINT64_MAX:
        raise ValueError("primesieve only generates primes smaller than 2^64")

    cdef uint64_t c_start = start
    cdef uint64_t c_stop = stop - 1  # primesieve bounds are inclusive
    cdef size_t size, i
    cdef uint64_t* primes

    global errno
    errno = 0
    sig_on()
    primes = <uint64_t*> primesieve_generate_primes(c_start, c_stop,
                                                    &size, UINT64_PRIMES)
    sig_off()
    if primes is NULL:
        if errno == EDOM:
            raise RuntimeError("primesieve_generate_primes failed")
        return []

    cdef list res = []
    try:
        if py_ints:
            for i in range(size):
                res.append(PyLong_FromUnsignedLongLong(primes[i]))
        else:
            for i in range(size):
                res.append(_new_integer(primes[i]))
    finally:
        primesieve_free(primes)
    return res


cdef class prime_iterator:
    r"""
    Iterator over the primes in the interval ``[start, stop)``.

    INPUT:

    - ``start`` -- integer (default: `2`); the first prime returned is
      the smallest prime `\geq` ``start``, which must be smaller
      than `2^{64}`

    - ``stop`` -- integer (default: ``None``); upper bound (excluded).
      If ``None``, iterate over all primes below `2^{64}`; when that
      limit is reached, :exc:`OverflowError` is raised.

    EXAMPLES::

        sage: from sage.libs.primesieve import prime_iterator
        sage: it = prime_iterator()
        sage: [next(it) for _ in range(6)]
        [2, 3, 5, 7, 11, 13]
        sage: list(prime_iterator(10, 30))
        [11, 13, 17, 19, 23, 29]
        sage: type(next(prime_iterator()))
        <class 'sage.rings.integer.Integer'>

    TESTS::

        sage: list(prime_iterator(10, 11))
        []
        sage: list(prime_iterator(-10, 5))
        [2, 3]
        sage: list(prime_iterator(0, 0))
        []
        sage: it = prime_iterator(2**64)
        Traceback (most recent call last):
        ...
        ValueError: start must be smaller than 2^64
        sage: it = prime_iterator(0, 2**64 + 1)
        Traceback (most recent call last):
        ...
        ValueError: stop must be at most 2^64

    The iterator stays exhausted once it has stopped::

        sage: it = prime_iterator(2, 4)
        sage: list(it), list(it)
        ([2, 3], [])

    Iterating beyond `2^{64}` is not possible; note that
    :exc:`StopIteration` is raised instead when the iterator is bounded,
    since in that case all requested primes have been returned::

        sage: # long time
        sage: it = prime_iterator(2**64 - 100)
        sage: next(it), next(it), next(it)
        (18446744073709551521, 18446744073709551533, 18446744073709551557)
        sage: next(it)
        Traceback (most recent call last):
        ...
        OverflowError: primesieve cannot generate primes beyond 2^64 - 1
        sage: next(it)
        Traceback (most recent call last):
        ...
        OverflowError: primesieve cannot generate primes beyond 2^64 - 1
        sage: list(prime_iterator(2**64 - 100, 2**64))
        [18446744073709551521, 18446744073709551533, 18446744073709551557]

    Starting beyond the largest prime below `2^{64}` is detected
    without any sieving::

        sage: list(prime_iterator(2**64 - 50, 2**64))
        []
        sage: next(prime_iterator(2**64 - 50))
        Traceback (most recent call last):
        ...
        OverflowError: primesieve cannot generate primes beyond 2^64 - 1
    """
    cdef primesieve_iterator it
    cdef uint64_t last  # largest value that may be returned
    cdef bint bounded
    cdef bint exhausted
    cdef bint overflowed

    def __cinit__(self):
        primesieve_init(&self.it)

    def __dealloc__(self):
        primesieve_free_iterator(&self.it)

    def __init__(self, start=2, stop=None):
        r"""
        TESTS::

            sage: from sage.libs.primesieve import prime_iterator
            sage: next(prime_iterator(5))
            5
        """
        if start < 0:
            start = 0
        elif start > UINT64_MAX:
            raise ValueError("start must be smaller than 2^64")
        self.exhausted = False
        self.overflowed = False
        if stop is None:
            self.bounded = False
            self.last = UINT64_MAX
        else:
            self.bounded = True
            if stop - 1 > UINT64_MAX:
                raise ValueError("stop must be at most 2^64")
            if stop <= 2 or stop <= start:
                self.exhausted = True
                self.last = 0
            else:
                self.last = stop - 1
        primesieve_jump_to(&self.it, start, self.last)

    def __iter__(self):
        r"""
        TESTS::

            sage: from sage.libs.primesieve import prime_iterator
            sage: it = prime_iterator()
            sage: iter(it) is it
            True
        """
        return self

    cdef int _exhaust(self) except -1:
        # Mark the iterator as exhausted and raise the appropriate
        # exception: all primes below 2^64 have been returned, which
        # for a bounded iterator simply means the end of the iteration.
        self.exhausted = True
        if self.bounded:
            raise StopIteration
        self.overflowed = True
        raise OverflowError("primesieve cannot generate primes "
                            "beyond 2^64 - 1")

    def __next__(self):
        r"""
        TESTS::

            sage: from sage.libs.primesieve import prime_iterator
            sage: next(prime_iterator(97))
            97
        """
        cdef uint64_t p
        if self.exhausted:
            if self.overflowed:
                raise OverflowError("primesieve cannot generate primes "
                                    "beyond 2^64 - 1")
            raise StopIteration
        # Inline part of primesieve_next_prime() so that only the
        # (rare) refill of the prime buffer pays for sig_on().
        self.it.i += 1
        if self.it.i >= self.it.size:
            if (self.it.size > 0
                    and self.it.primes[self.it.size - 1] >= largest_uint64_prime):
                # the buffer already ended with the largest prime < 2^64
                self._exhaust()
            if self.it.size == 0 and self.it.start > largest_uint64_prime:
                # fresh iterator starting beyond the largest prime < 2^64
                self._exhaust()
            # If the refill is interrupted (KeyboardInterrupt), the C
            # iterator state may be inconsistent, so leave the iterator
            # exhausted rather than risk returning wrong primes.
            self.exhausted = True
            sig_on()
            primesieve_generate_next_primes(&self.it)
            sig_off()
            self.exhausted = False
        p = self.it.primes[self.it.i]
        if self.it.is_error or p == UINT64_MAX:
            # unexpected failure inside primesieve
            self._exhaust()
        if self.bounded and p > self.last:
            self.exhausted = True
            raise StopIteration
        return _new_integer(p)
