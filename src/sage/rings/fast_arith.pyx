# sage_setup: distribution = sagemath-categories
# sage.doctest: needs sage.libs.pari
"""
Basic arithmetic with C integers
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

from libc.math cimport sqrt

from sage.rings.integer cimport Integer

cpdef prime_range(start, stop=None, algorithm=None, bint py_ints=False):
    r"""
    Return a list of all primes between ``start`` and ``stop - 1``, inclusive.

    If the second argument is omitted, this returns the primes up to the
    first argument.

    The sage command :func:`~sage.arith.misc.primes` is an alternative that
    uses less memory (but may be slower), because it returns an iterator,
    rather than building a list of the primes.

    INPUT:

    - ``start`` -- integer; lower bound (default: 1)

    - ``stop`` -- integer; upper bound

    - ``algorithm`` -- string (default: ``None``), one of:

      - ``None``: Use  algorithm ``'pari_primes'`` if ``stop`` <= 436273009
        (approximately 4.36E8). Otherwise use algorithm ``'pari_isprime'``.

      - ``'pari_primes'``: Use PARI's :pari:`primes` function to generate all
        primes from 2 to stop. This is fast but may crash if there is
        insufficient memory. Raises an error if ``stop`` > 436273009.

      - ``'pari_isprime'``: Wrapper for ``list(primes(start, stop))``. Each (odd)
        integer in the specified range is tested for primality by applying PARI's
        :pari:`isprime` function. This is slower but will work for much larger input.

    - ``py_ints`` -- boolean (default: ``False``); return Python ints rather
      than Sage Integers (faster). Ignored unless algorithm ``'pari_primes'`` is being
      used.

    EXAMPLES::

        sage: prime_range(10)
        [2, 3, 5, 7]
        sage: prime_range(7)
        [2, 3, 5]
        sage: prime_range(2000,2020)
        [2003, 2011, 2017]
        sage: prime_range(2,2)
        []
        sage: prime_range(2,3)
        [2]
        sage: prime_range(5,10)
        [5, 7]
        sage: prime_range(-100,10,"pari_isprime")
        [2, 3, 5, 7]
        sage: prime_range(2,2,algorithm='pari_isprime')
        []
        sage: prime_range(10**16,10**16+100,"pari_isprime")
        [10000000000000061, 10000000000000069, 10000000000000079, 10000000000000099]
        sage: prime_range(10**30,10**30+100,"pari_isprime")
        [1000000000000000000000000000057, 1000000000000000000000000000099]
        sage: type(prime_range(8)[0])
        <class 'sage.rings.integer.Integer'>
        sage: type(prime_range(8,algorithm='pari_isprime')[0])
        <class 'sage.rings.integer.Integer'>

    .. NOTE::

        ``start`` and ``stop`` should be integers, but real numbers will also be accepted
        as input. In this case, they will be rounded to nearby integers start\* and
        stop\*, so the output will be the primes between start\* and stop\* - 1, which may
        not be exactly the same as the primes between ``start`` and ``stop - 1``.

    TESTS::

        sage: prime_range(-1)
        []
        sage: L = prime_range(25000,2500000)
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
        ValueError: algorithm must be "pari_primes" or "pari_isprime"

    Confirm the fixes for :issue:`28467`::

        sage: prime_range(436273009, 436273010)
        [436273009]
        sage: prime_range(436273009, 436273010, algorithm='pari_primes')
        Traceback (most recent call last):
        ...
        ValueError: algorithm "pari_primes" is limited to primes larger than 436273008

    AUTHORS:

    - William Stein (original version)
    - Craig Citro (rewrote for massive speedup)
    - Kevin Stueve (added primes iterator option) 2010-10-16
    - Robert Bradshaw (speedup using Pari prime table, py_ints option)
    """
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

    if algorithm is None:
        # if 'stop' is 'None', need to change it to an integer before comparing with 'start'
        if max(start, stop or 0) <= small_prime_max:
            algorithm = "pari_primes"
        else:
            algorithm = "pari_isprime"

    if algorithm == "pari_primes":
        from sage.libs.pari.convert_sage import pari_maxprime, pari_prime_range
        from sage.libs.pari import pari

        if max(start, stop or 0) > small_prime_max:
            raise ValueError('algorithm "pari_primes" is limited to primes '
                             f'larger than {small_prime_max - 1}')

        if stop is None:
            # In this case, "start" is really stop
            stop = start
            start = 1
        else:
            start = start
            stop = stop
            if start < 1:
                start = 1
        if stop <= start:
            return []

        if pari_maxprime() < stop:
            # Adding prime_gap_bound should be sufficient to guarantee an
            # additional prime, given that c_stop <= small_prime_max.
            pari.init_primes(min(stop + prime_gap_bound, init_primes_max))
            assert pari_maxprime() >= stop

        res = pari_prime_range(start, stop, py_ints)

    elif algorithm == "pari_isprime" or algorithm == "pari_primes":
        from sage.arith.misc import primes
        res = list(primes(start, stop))
    else:
        raise ValueError('algorithm must be "pari_primes" or "pari_isprime"')
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
        if x <= bnd and self.gcd_longlong(x, y) == 1:
            n[0] = y
            d[0] = x
            return 0

        n[0] = 0
        d[0] = 0
        return 0

    def rational_recon_longlong(self, long long a, long long m):
        """
        Rational reconstruction of a modulo m.
        """
        cdef long long n, d
        self.c_rational_recon_longlong(a, m, &n, &d)
        return (n, d)
