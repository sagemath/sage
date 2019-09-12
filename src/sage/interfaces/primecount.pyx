r"""
Interface to the primecount library
"""
#*****************************************************************************
#       Copyright (C) 2018 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from libc.stdint cimport int64_t
from libcpp.string cimport string as cppstring
from cpython.int cimport PyInt_FromString

from cysignals.signals cimport sig_on, sig_off

cimport sage.libs.primecount as primecount

cdef inline int _do_sig(int64_t n):
    "threshold for sig_on/sig_off"
    return n >> 26

cpdef int64_t prime_pi(int64_t n, method=None) except -1:
    r"""
    Return the number of prime numbers smaller or equal than ``n``.

    INPUT:

    - ``n`` - an integer

    - ``method`` - ``None`` or a string that determines the primecount
      function to be called

        - ``"deleglise_rivat"``
        - ``"legendre"``
        - ``"lehmer"``
        - ``"lmo"``
        - ``"meissel"``
        - ``"primesieve"``

    EXAMPLES::

        sage: from sage.interfaces.primecount import prime_pi # optional - primecount

        sage: prime_pi(1000) == 168                     # optional - primecount
        True
        sage: prime_pi(1000, "deleglise_rivat") == 168  # optional - primecount
        True
        sage: prime_pi(1000, "legendre") == 168         # optional - primecount
        True
        sage: prime_pi(1000, "lehmer") == 168           # optional - primecount
        True
        sage: prime_pi(1000, "lmo") == 168              # optional - primecount
        True
        sage: prime_pi(1000, "meissel") == 168          # optional - primecount
        True
        sage: prime_pi(1000, "primesieve") == 168       # optional - primecount
        True
        sage: prime_pi(1000, "youpi")                   # optional - primecount
        Traceback (most recent call last):
        ...
        ValueError: unknown method 'youpi'
    """
    cdef int64_t ans
    if method is None:
        if _do_sig(n): sig_on()
        ans = primecount.pi(n)
        if _do_sig(n): sig_off()
    elif method == "deleglise_rivat":
        if _do_sig(n): sig_on()
        ans = primecount.pi_deleglise_rivat(n)
        if _do_sig(n): sig_off()
    elif method == "legendre":
        if _do_sig(n): sig_on()
        ans = primecount.pi_legendre(n)
        if _do_sig(n): sig_off()
    elif method == "lehmer":
        if _do_sig(n): sig_on()
        ans = primecount.pi_lehmer(n)
        if _do_sig(n): sig_off()
    elif method == "lmo":
        if _do_sig(n): sig_on()
        ans = primecount.pi_lmo(n)
        if _do_sig(n): sig_off()
    elif method == "meissel":
        if _do_sig(n): sig_on()
        ans = primecount.pi_meissel(n)
        if _do_sig(n): sig_off()
    elif method == "primesieve":
        if _do_sig(n): sig_on()
        ans = primecount.pi_primesieve(n)
        if _do_sig(n): sig_off()
    else:
        raise ValueError("unknown method {!r}".format(method))

    return ans

cpdef prime_pi_128(n):
    r"""
    Return the number of prime number smaller than ``n``.

    EXAMPLES::

        sage: from sage.interfaces.primecount import prime_pi_128 # optional - primecount

        sage: prime_pi_128(1000)     # optional - primecount
        168
        sage: nth_prime_128(2**65)   # not tested
        ?
    """
    cdef cppstring s = str(n).encode('ascii')
    cdef bytes ans
    sig_on()
    ans = primecount.pi(s)
    sig_off()
    return PyInt_FromString(ans, NULL, 10)

cpdef int64_t nth_prime(int64_t n) except -1:
    r"""
    Return the ``n``-th prime integer.

    EXAMPLES::

        sage: from sage.interfaces.primecount import nth_prime # optional - primecount

        sage: nth_prime(168) == 997  # optional - primecount
        True
    """
    if n <= 0:
        raise ValueError("n must be positive")

    cdef int64_t ans
    if _do_sig(n): sig_on()
    ans = primecount.nth_prime(n)
    if _do_sig(n): sig_off()
    return ans

cpdef int64_t phi(int64_t x, int64_t a):
    r"""
    Return the number of integers smaller or equal than ``x`` by any of the
    first ``a`` primes.

    This is sometimes called a "partial sieve function" or "Legendre-sum".

    EXAMPLES::

         sage: from sage.interfaces.primecount import phi # optional - primecount

         sage: phi(1000, 3) == 266 # optional - primecount
         True
         sage: phi(2**30, 100) == 95446716 # optional - primecount
         True
    """
    return primecount.phi(x, a)

