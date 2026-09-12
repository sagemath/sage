# sage.doctest: needs sage.libs.pari
"""
Basic arithmetic with C integers

TESTS:

The integer arithmetic helper classes cannot be subclassed::

    sage: from sage.rings.fast_arith import arith_int, arith_llong
    sage: type("arith_int_subclass", (arith_int,), {})
    Traceback (most recent call last):
    ...
    TypeError: type 'sage.rings.fast_arith.arith_int' is not an acceptable base type
    sage: type("arith_llong_subclass", (arith_llong,), {})
    Traceback (most recent call last):
    ...
    TypeError: type 'sage.rings.fast_arith.arith_llong' is not an acceptable base type
"""

# ****************************************************************************
#       Copyright (C) 2004 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

###################################################################
# We define the following functions in this file, both
# for int (up to bound = 2**31 - 1) and longlong (up to 2**63 - 1).
# The function definitions are identical except for the types.
# Some of their input can be at most sqrt(bound), since
# it is necessary to multiply numbers and reduce the product
# modulo n, where n is at most bound.
#
#   * abs_int -- absolute value of integer
#   * sign_int -- sign of integer
#   * c_gcd_int -- gcd of two ints
#   * gcd_int -- python export of c_gcd_int
#   * c_xgcd_int -- extended gcd of two ints
#   * c_inverse_mod_int -- inverse of an int modulo another int
#   * c_rational_recon_int -- rational reconstruction of ints
#   * rational_recon_int -- export of rational reconstruction for ints
#
#  The long long functions are the same, except they end in _longlong.
#
###################################################################

# The int definitions

from libc.limits cimport ULONG_MAX
from libc.math cimport sqrt

from cysignals.signals cimport sig_on, sig_off
from sage.ext.stdsage cimport PY_NEW
from sage.libs.gmp.mpz cimport mpz_set_ui
from sage.libs.flint.ulong_extras cimport (
    n_primes_t, n_primes_init, n_primes_clear, n_primes_next,
    n_primes_sieve_range, n_primes_jump_after
)

from sage.rings.integer cimport Integer


cdef list _flint_prime_range(unsigned long start, unsigned long stop, bint py_ints=False):
    r"""
    Return a list of primes in ``[start, stop)`` using FLINT's segmented sieve.
    """
    cdef n_primes_t it
    cdef list res = []
    cdef unsigned long seg_start = start
    cdef unsigned long seg_stop
    # FLINT requires (odd_b - odd_a) < 65536.
    # Setting chunk_size to 65534 guarantees (seg_stop - seg_start) <= 65534
    # so odd_b - odd_a <= 65534 < 65536 for all parities of seg_start and seg_stop.
    cdef unsigned long chunk_size = 65534
    cdef unsigned long p
    cdef Integer z

    if stop <= 2 or start >= stop:
        return []
    if start < 2:
        seg_start = 2

    n_primes_init(it)
    try:
        while seg_start < stop:
            if stop - seg_start > chunk_size:
                seg_stop = seg_start + chunk_size
            else:
                seg_stop = stop
            sig_on()
            n_primes_sieve_range(it, seg_start, seg_stop)
            n_primes_jump_after(it, seg_start - 1)
            sig_off()
            while True:
                p = n_primes_next(it)
                if p >= seg_stop or p == 0:
                    break
                if py_ints:
                    res.append(p)
                else:
                    z = <Integer>PY_NEW(Integer)
                    mpz_set_ui(z.value, p)
                    res.append(z)
            seg_start = seg_stop
    finally:
        n_primes_clear(it)
    return res


cpdef prime_range(start, stop=None, step=None, algorithm=None, bint py_ints=False):
    r"""
    Return a list of all primes between ``start`` and ``stop - 1``, inclusive
    (or between ``stop`` and ``start + 1`` if ``step`` is negative).

    If the second argument is omitted, this returns the primes up to the
    first argument.

    .. SEEALSO::

        - :func:`~sage.arith.misc.primes` is an alternative that
          uses less memory (but may be slower), because it returns an iterator,
          rather than building a list of the primes.

        - :class:`~sage.sets.primes.Primes` can be used to create sets of primes
          with more complicated congruence conditions.

    INPUT:

    - ``start`` -- integer; lower bound (default: 1)

    - ``stop`` -- integer; upper bound

    - ``step`` -- integer or ``None`` (default: ``None``); if not ``None``,
      the function returns only primes that are congruent to ``start`` modulo
      ``step``. If ``step`` is negative, then the returned list will be
      decreasing.

    - ``algorithm`` -- string (default: ``None``), one of:

      - ``None``: Use algorithm ``'pari_primes'`` if ``stop`` <= 436273009
        (approximately 4.36E8). Otherwise, for ranges up to `2^{64}-1`, use
        FLINT's segmented sieve (algorithm ``'flint'``) when the interval is
        sufficiently large, or ``'pari_isprime'`` for small intervals at large
        numerical offsets. If ``stop`` exceeds `2^{64}-1`, use algorithm
        ``'pari_isprime'``.

      - ``'pari_primes'``: Use PARI's :pari:`primes` function to generate all
        primes from 2 to stop. This is fast but may crash if there is
        insufficient memory. Raises an error if ``stop`` > 436273009.

      - ``'flint'``: Use FLINT's segmented sieve (:c:func:`n_primes_sieve_range`).
        Generates primes up to `2^{64}-1` using a segmented sieve, which is
        substantially faster than individual primality tests.

      - ``'pari_isprime'``: Wrapper for ``list(primes(start, stop))``. Each (odd)
        integer in the specified range is tested for primality by applying PARI's
        :pari:`isprime` function. This is slower but will work for arbitrarily large input.

    - ``py_ints`` -- boolean (default: ``False``); return Python ints rather
      than Sage Integers (faster). Supported by algorithms ``'pari_primes'``
      and ``'flint'``.

    EXAMPLES::

        sage: prime_range(10)
        [2, 3, 5, 7]
        sage: prime_range(7)
        [2, 3, 5]
        sage: prime_range(2000, 2020)
        [2003, 2011, 2017]
        sage: prime_range(2, 2)
        []
        sage: prime_range(2, 3)
        [2]
        sage: prime_range(5, 10)
        [5, 7]
        sage: prime_range(11, 100, 10)
        [11, 31, 41, 61, 71]
        sage: prime_range(-100, 10, "pari_isprime")
        [2, 3, 5, 7]
        sage: prime_range(2, 2, algorithm='pari_isprime')
        []
        sage: prime_range(10**16, 10**16+100, "pari_isprime")
        [10000000000000061, 10000000000000069, 10000000000000079, 10000000000000099]
        sage: prime_range(10**30, 10**30+100, "pari_isprime")
        [1000000000000000000000000000057, 1000000000000000000000000000099]
        sage: type(prime_range(8)[0])
        <class 'sage.rings.integer.Integer'>
        sage: type(prime_range(8, algorithm='pari_isprime')[0])
        <class 'sage.rings.integer.Integer'>

    .. NOTE::

        ``start`` and ``stop`` should be integers, but real numbers will also be accepted
        as input. In this case, they will be rounded to nearby integers start\* and
        stop\*, so the output will be the primes between start\* and stop\* - 1, which may
        not be exactly the same as the primes between ``start`` and ``stop - 1``.

    TESTS::

        sage: prime_range(-1)
        []
        sage: L = prime_range(25000, 2500000)
        sage: len(L)
        180310
        sage: L[-10:]
        [2499923, 2499941, 2499943, 2499947, 2499949, 2499953, 2499967, 2499983, 2499989, 2499997]

    A non-trivial range without primes::

        sage: prime_range(4652360, 4652400)
        []

    Test for non-existing algorithm::

        sage: prime_range(55, algorithm='banana')
        Traceback (most recent call last):
        ...
        ValueError: algorithm must be "pari_primes", "pari_isprime", or "flint"

    Confirm the fixes for :issue:`28467`::

        sage: prime_range(436273009, 436273010)
        [436273009]
        sage: prime_range(436273009, 436273010, algorithm='pari_primes')
        Traceback (most recent call last):
        ...
        ValueError: algorithm "pari_primes" is limited to primes larger than 436273008

    Tests for FLINT segmented sieve (:issue:`42751`)::

        sage: prime_range(2000, 2020, algorithm='flint')
        [2003, 2011, 2017]
        sage: prime_range(10, algorithm='flint')
        [2, 3, 5, 7]
        sage: prime_range(2, 2, algorithm='flint')
        []
        sage: prime_range(2, 3, algorithm='flint')
        [2]
        sage: prime_range(3, 4, algorithm='flint')
        [3]
        sage: prime_range(-5, 5, algorithm='flint')
        [2, 3]

    Segment boundaries around 65536 (FLINT sieve chunk size)::

        sage: prime_range(65520, 65540, algorithm='flint')
        [65521, 65537, 65539]
        sage: prime_range(65536, 65538, algorithm='flint')
        [65537]
        sage: prime_range(65535, 65537, algorithm='flint')
        []
        sage: prime_range(65537, 65538, algorithm='flint')
        [65537]
        sage: prime_range(65537, 65539, algorithm='flint')
        [65537]

    Exact and multi-segment sizes::

        sage: len(prime_range(2, 65536, algorithm='flint'))
        6542
        sage: len(prime_range(2, 65536 * 2, algorithm='flint'))
        12251
        sage: prime_range(1, 100000, algorithm='flint') == prime_range(1, 100000, algorithm='pari_primes')
        True

    Large ranges and crossover beyond 436273009::

        sage: prime_range(10^12, 10^12 + 100, algorithm='flint')
        [1000000000039, 1000000000061, 1000000000063, 1000000000091]
        sage: prime_range(10^12, 10^12 + 100) == prime_range(10^12, 10^12 + 100, algorithm='pari_isprime')
        True
        sage: prime_range(10^12, 10^12 + 10000) == prime_range(10^12, 10^12 + 10000, algorithm='flint')
        True

    Step and negative step with FLINT::

        sage: prime_range(11, 100, 10, algorithm='flint')
        [11, 31, 41, 61, 71]
        sage: prime_range(20, 10, -1, algorithm='flint')
        [19, 17, 13, 11]
        sage: prime_range(65540, 65520, -1, algorithm='flint')
        [65539, 65537, 65521]
        sage: prime_range(100, 2, -7, algorithm='flint') == [p for p in range(100, 2, -7) if is_prime(p)]
        True
        sage: prime_range(2^32 - 1, 2^32 - 100, -1, algorithm='flint')
        [4294967291, 4294967279, 4294967231, 4294967197]
        sage: prime_range(2^64 - 1, 2^64 - 100, -1, algorithm='flint')  # long time, needs !32_bit
        [18446744073709551557, 18446744073709551533, 18446744073709551521]

    Test py_ints option with FLINT::

        sage: P_flint = prime_range(10, algorithm='flint', py_ints=True)
        sage: P_flint
        [2, 3, 5, 7]
        sage: type(P_flint[0])
        <class 'int'>
        sage: type(prime_range(10, algorithm='flint', py_ints=False)[0])
        <class 'sage.rings.integer.Integer'>

    Input exceeding word size raises error with FLINT::

        sage: prime_range(10^30, 10^30 + 10, algorithm='flint')
        Traceback (most recent call last):
        ...
        OverflowError: algorithm "flint" does not support primes larger than ...

    Some step tests:

        sage: prime_range(4, 15, 3)
        [7, 13]
        sage: prime_range(10, 20, -1)
        []
        sage: prime_range(20, 10, 1)
        []
        sage: prime_range(20, 10, -1)
        [19, 17, 13, 11]
        sage: prime_range(96, 19, -1)
        [89, 83, 79, 73, 71, 67, 61, 59, 53, 47, 43, 41, 37, 31, 29, 23]


    Make sure that step behaves exactly like in range::

        sage: # needs sage.libs.pari
        sage: a = randint(1, 50)
        sage: b = randint(70, 100)
        sage: step = randint(1, 5)
        sage: prime_range(a, b, step) == list(filter(is_prime, range(a, b, step)))
        True
        sage: a = randint(50, 100)
        sage: b = randint(0, 30)
        sage: step = randint(-5, -1)
        sage: prime_range(a, b, step) == list(filter(is_prime, range(a, b, step)))
        True

    AUTHORS:

    - William Stein (original version)
    - Craig Citro (rewrote for massive speedup)
    - Kevin Stueve (added primes iterator option) 2010-10-16
    - Robert Bradshaw (speedup using Pari prime table, py_ints option)
    - Vincent Macri (added step option)
    """
    if isinstance(step, str):
        # For backwards compatibility - `algorithm` used to be the third parameter.
        # We make sure that previous code still works by treating `step` as
        # `algorithm` if `step` is a string and `algorithm` is None.
        if isinstance(algorithm, bool):
            py_ints = algorithm
        elif algorithm is not None:
            raise TypeError('step must be an integer or None')
        algorithm = step
        step = None

    # input to pari.init_primes cannot be greater than 436273290 (hardcoded bound)
    DEF init_primes_max = 436273290
    DEF small_prime_max = 436273009  # a prime < init_primes_max (preferably the largest)
    DEF prime_gap_bound = 250        # upper bound for gap between primes <= small_prime_max

    # make sure that start and stop are integers
    # First try coercing them. If that does not work, then try rounding them.
    try:
        start = Integer(start)
    except TypeError as integer_error:
        try:
            start = Integer(round(float(start)))
        except (ValueError, TypeError) as real_error:
            raise TypeError(str(integer_error)
                            + "\nand argument is also not real: "
                            + str(real_error))
    if stop is not None:
        try:
            stop = Integer(stop)
        except TypeError as integer_error:
            try:
                stop = Integer(round(float(stop)))
            except (ValueError, TypeError) as real_error:
                raise ValueError(str(integer_error)
                                 + "\nand argument is also not real: "
                                 + str(real_error))

    if step is not None:
        if not isinstance(step, (Integer, int)):
            raise TypeError('step must be an integer or None')
        step = Integer(step)
    else:
        step = 1

    if algorithm is None:
        # if 'stop' is 'None', need to change it to an integer before comparing with 'start'
        if max(start, stop or 0) <= small_prime_max:
            algorithm = "pari_primes"
        elif max(start, stop or 0) <= ULONG_MAX:
            # Segmented sieve in FLINT is substantially faster than PARI isprime,
            # except for tiny intervals at large offsets where primality-testing
            # a few integers is faster than sieving base primes up to sqrt(stop).
            if stop is None or abs(stop - start) >= 10000 or max(start, stop) <= 10**11:
                algorithm = "flint"
            else:
                algorithm = "pari_isprime"
        else:
            algorithm = "pari_isprime"

    if algorithm == "pari_primes":
        from sage.libs.pari.convert_sage import pari_maxprime, pari_prime_range
        from sage.libs.pari import pari

        if max(start, stop or 0) > small_prime_max:
            raise ValueError('algorithm "pari_primes" is limited to primes '
                             f'larger than {small_prime_max - 1}')

        congruence = start % step
        if stop is None:
            # In this case, "start" is really stop
            stop = start
            start = 1
        else:
            start = start
            stop = stop

        if step < 1:
            start, stop = stop + 1, start + 1

        if stop <= start:
            return []

        if pari_maxprime() < stop:
            # Adding prime_gap_bound should be sufficient to guarantee an
            # additional prime, given that c_stop <= small_prime_max.
            pari.init_primes(min(stop + prime_gap_bound, init_primes_max))
            assert pari_maxprime() >= stop

        res = pari_prime_range(max(start, 1), stop, py_ints)
        if step < 0:
            res = res[::-1]
        if step != 1 and step != -1:
            res = [p for p in res if p % step == congruence]

    elif algorithm == "flint":
        if max(start, stop or 0) > ULONG_MAX:
            raise OverflowError('algorithm "flint" does not support primes larger than '
                                f'{ULONG_MAX}')

        congruence = start % step
        if stop is None:
            stop = start
            start = 1

        if step < 1:
            start, stop = stop + 1, start + 1

        if stop > ULONG_MAX:
            stop = ULONG_MAX

        if stop <= start or stop <= 2:
            return []

        c_start = max(int(start), 2)
        c_stop = int(stop)

        res = _flint_prime_range(c_start, c_stop, py_ints)
        if step < 0:
            res.reverse()
        if step != 1 and step != -1:
            res = [p for p in res if p % step == congruence]

    elif algorithm == "pari_isprime":
        from sage.arith.misc import primes
        res = list(primes(start, stop, step))
    else:
        raise ValueError('algorithm must be "pari_primes", "pari_isprime", or "flint"')
    return res


cdef class arith_int:
    cdef int abs_int(self, int x) except -1:
        if x < 0:
            return -x
        return x

    cdef int sign_int(self, int n) except -2:
        if n < 0:
            return -1
        return 1

    cdef int c_gcd_int(self, int a, int b) except -1:
        cdef int c
        if a == 0:
            return self.abs_int(b)
        if b == 0:
            return self.abs_int(a)
        if a < 0:
            a = -a
        if b < 0:
            b = -b
        while b:
            c = a % b
            a = b
            b = c
        return a

    def gcd_int(self, int a, int b):
        return self.c_gcd_int(a, b)

    cdef int c_xgcd_int(self, int a, int b, int* ss, int* tt) except -1:
        cdef int psign, qsign, p, q, r, s, c, quot, new_r, new_s

        if a == 0:
            ss[0] = 0
            tt[0] = self.sign_int(b)
            return self.abs_int(b)

        if b == 0:
            ss[0] = self.sign_int(a)
            tt[0] = 0
            return self.abs_int(a)

        psign = 1
        qsign = 1

        if a < 0:
            a = -a
            psign = -1
        if b < 0:
            b = -b
            qsign = -1

        p = 1
        q = 0
        r = 0
        s = 1
        while b:
            c = a % b
            quot = a / b
            a = b
            b = c
            new_r = p - quot * r
            new_s = q - quot * s
            p = r
            q = s
            r = new_r
            s = new_s

        ss[0] = p * psign
        tt[0] = q * qsign

        return a

    def xgcd_int(self, int a, int b):
        cdef int g, s, t
        g = self.c_xgcd_int(a, b, &s, &t)
        return (g, s, t)

    cdef int c_inverse_mod_int(self, int a, int m) except -1:
        if a == 1 or m <= 1:
            return a % m   # common special case
        cdef int g, s, t
        g = self.c_xgcd_int(a, m, &s, &t)
        if g != 1:
            raise ArithmeticError("The inverse of %s modulo %s is not defined." % (a, m))
        s = s % m
        if s < 0:
            s = s + m
        return s

    def inverse_mod_int(self, int a, int m):
        return self.c_inverse_mod_int(a, m)

    cdef int c_rational_recon_int(self, int a, int m, int* n, int* d) except -1:
        cdef int u, v, u0, u1, u2, v0, v1, v2, q, t0, t1, t2, x, y
        cdef float bnd

        if m > 46340:
            raise OverflowError(f"The modulus m(={m}) should be at most 46340")

        a = a % m

        if a == 0 or m == 0:
            n[0] = 0
            d[0] = 1
            return 0

        if m < 0:
            m = -m
        if a < 0:
            a = m - a
        if a == 1:
            n[0] = 1
            d[0] = 1
            return 0

        u = m
        v = a
        bnd = sqrt(m/2.0)
        u0 = 1
        u1 = 0
        u2 = u
        v0 = 0
        v1 = 1
        v2 = v
        while self.abs_int(v2) > bnd:
            q = u2 / v2   # floor is implicit
            t0 = u0 - q * v0
            t1 = u1 - q * v1
            t2 = u2 - q * v2
            u0 = v0
            u1 = v1
            u2 = v2
            v0 = t0
            v1 = t1
            v2 = t2

        x = self.abs_int(v1)
        y = v2
        if v1 < 0:
            y = -1*y
        if x <= bnd and self.c_gcd_int(x, y) == 1:
            n[0] = y
            d[0] = x
            return 0

        n[0] = 0
        d[0] = 0
        return 0

    def rational_recon_int(self, int a, int m):
        """
        Rational reconstruction of a modulo m.
        """
        cdef int n, d
        self.c_rational_recon_int(a, m, &n, &d)
        return (n, d)


# The long long versions are next.
cdef class arith_llong:

    cdef long long abs_longlong(self, long long x) except -1:
        if x < 0:
            return -x
        return x

    cdef long long sign_longlong(self, long long n) except -2:
        if n < 0:
            return -1
        return 1

    cdef long long c_gcd_longlong(self, long long a, long long b) except -1:
        cdef long long c
        if a == 0:
            return self.abs_longlong(b)
        if b == 0:
            return self.abs_longlong(a)
        if a < 0:
            a = -a
        if b < 0:
            b = -b
        while b:
            c = a % b
            a = b
            b = c
        return a

    def gcd_longlong(self, long long a, long long b):
        return self.c_gcd_longlong(a, b)

    cdef long long c_xgcd_longlong(self, long long a, long long b,
                                   long long *ss,
                                   long long *tt) except -1:
        cdef long long psign, qsign, p, q, r, s, c, quot, new_r, new_s

        if a == 0:
            ss[0] = 0
            tt[0] = self.sign_longlong(b)
            return self.abs_longlong(b)

        if b == 0:
            ss[0] = self.sign_longlong(a)
            tt[0] = 0
            return self.abs_longlong(a)

        psign = 1
        qsign = 1

        if a < 0:
            a = -a
            psign = -1
        if b < 0:
            b = -b
            qsign = -1

        p = 1
        q = 0
        r = 0
        s = 1
        while b:
            c = a % b
            quot = a / b
            a = b
            b = c
            new_r = p - quot * r
            new_s = q - quot * s
            p = r
            q = s
            r = new_r
            s = new_s

        ss[0] = p * psign
        tt[0] = q * qsign

        return a

    cdef long long c_inverse_mod_longlong(self, long long a, long long m) except -1:
        cdef long long g, s, t
        g = self.c_xgcd_longlong(a, m, &s, &t)
        if g != 1:
            raise ArithmeticError("The inverse of %s modulo %s is not defined." % (a, m))
        s = s % m
        if s < 0:
            s = s + m
        return s

    def inverse_mod_longlong(self, long long a, long long m):
        return self.c_inverse_mod_longlong(a, m)

    cdef long long c_rational_recon_longlong(self, long long a, long long m,
                                             long long *n, long long *d) except -1:
        cdef long long u, v, u0, u1, u2, v0, v1, v2, q, t0, t1, t2, x, y
        cdef float bnd

        if m > 2147483647:
            raise OverflowError(f"The modulus m(={m}) must be at most 2147483647")

        a = a % m

        if a == 0 or m == 0:
            n[0] = 0
            d[0] = 1
            return 0

        if m < 0:
            m = -m
        if a < 0:
            a = m - a
        if a == 1:
            n[0] = 1
            d[0] = 1
            return 0

        u = m
        v = a
        bnd = sqrt(m/2.0)
        u0 = 1
        u1 = 0
        u2 = u
        v0 = 0
        v1 = 1
        v2 = v
        while self.abs_longlong(v2) > bnd:
            q = u2 / v2   # floor is implicit
            t0 = u0 - q * v0
            t1 = u1 - q * v1
            t2 = u2 - q * v2
            u0 = v0
            u1 = v1
            u2 = v2
            v0 = t0
            v1 = t1
            v2 = t2

        x = self.abs_longlong(v1)
        y = v2
        if v1 < 0:
            y = -1*y
        if x <= bnd and self.c_gcd_longlong(x, y) == 1:
            n[0] = y
            d[0] = x
            return 0

        n[0] = 0
        d[0] = 0
        return 0

    def rational_recon_longlong(self, long long a, long long m):
        """
        Rational reconstruction of a modulo m.

        EXAMPLES::

            sage: from sage.rings.fast_arith import arith_llong
            sage: arith_llong().rational_recon_longlong(1234567, 2147483629)
            (-23606, 22613)
        """
        cdef long long n, d
        self.c_rational_recon_longlong(a, m, &n, &d)
        return (n, d)
