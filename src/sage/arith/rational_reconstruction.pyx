"""
Rational reconstruction

This file is a Cython implementation of rational reconstruction, using
direct MPIR calls.

AUTHORS:

- ??? (2006 or before)

- Jeroen Demeyer (2014-10-20): move this function from ``gmp.pxi``,
  simplify and fix some bugs, see :issue:`17180`
"""
# ****************************************************************************
#       Copyright (C) 2006 ???
#       Copyright (C) 2014 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.signals cimport sig_on, sig_off

from sage.libs.gmp.mpz cimport *
from sage.libs.gmp.mpq cimport *


cdef int mpq_rational_reconstruction(mpq_t answer, mpz_t a, mpz_t m) except -1:
    """
    Set ``answer`` to a rational number which is `a` modulo `m` and
    such that the numerator and denominator of the result is bounded by
    sqrt(m/2).

    If `m` is zero, raise :class`ZeroDivisionError`. If the rational
    reconstruction does not exist, raise :exc:`ValueError`.

    We assume that ``mpq_init`` has been called on ``answer``.

    For examples, see the ``rational_reconstruction`` method of
    :class:`Integer`.

    TESTS::

        sage: q = -113/355
        sage: for p in range(2*355^2, 3*355^2):  # long time
        ....:     if is_prime(p):
        ....:         assert Mod(q, p).rational_reconstruction() == q
    """
    cdef mpz_t bound
    cdef mpz_t u1
    cdef mpz_t u2
    cdef mpz_t v1
    cdef mpz_t v2
    cdef mpz_t q

    if mpz_sgn(m) == 0:
        raise ZeroDivisionError("rational reconstruction with zero modulus")

    sig_on()
    mpz_init(bound)
    mpz_init_set_ui(u1, 0)
    mpz_init_set_ui(v1, 1)
    mpz_init(u2)
    mpz_init(v2)
    mpz_init(q)

    try:
        mpz_abs(u2, m)                 # u2 = |m|
        mpz_mod(v2, a, u2)             # v2 = a % m
        mpz_fdiv_q_2exp(bound, u2, 1)
        mpz_sqrtrem(bound, q, bound)   # bound = floor(sqrt(|m|/2))

        # Initialization: u1 = 0; v1 = 1
        while True:
            if mpz_cmpabs(v2, bound) <= 0:
                break
            mpz_fdiv_q(q, u2, v2)  # q = floor(u2/v2)
            mpz_submul(u1, q, v1)  # tmp1 = u1 - q*v1   (store tmp1 in u1)
            mpz_submul(u2, q, v2)  # tmp2 = u2 - q*v2   (store tmp2 in u2)
            mpz_swap(u1, v1)       # u1 = v1; v1 = tmp1
            mpz_swap(u2, v2)       # u2 = v2; v2 = tmp2

        # The answer is v2/v1, but check the conditions first
        if mpz_cmpabs(v1, bound) <= 0:
            mpz_gcd(q, v1, v2)
            if mpz_cmp_ui(q, 1) == 0:
                # Set answer to v2/v1 with correct sign
                if mpz_sgn(v1) >= 0:
                    mpz_set(mpq_numref(answer), v2)
                    mpz_set(mpq_denref(answer), v1)
                else:
                    mpz_neg(mpq_numref(answer), v2)
                    mpz_neg(mpq_denref(answer), v1)
                sig_off()
                return 0

        sig_off()
        # Earlier versions used to put a string representation of
        # the input here. We don't do this, since some low-level code
        # needs to catch and handle this exception so the string
        # handling is just a waste of time. The actual user-facing
        # methods could catch ValueError and raise a better exception
        # instead.
        raise ValueError("rational reconstruction does not exist")
    finally:
        mpz_clear(bound)
        mpz_clear(u1)
        mpz_clear(v1)
        mpz_clear(u2)
        mpz_clear(v2)
        mpz_clear(q)
