"""
Counting Primes

AUTHORS:

- \R. Andrew Ohana (2009): initial version of efficient prime_pi

- William Stein (2009): fix plot method

- \R. Andrew Ohana (2011): complete rewrite, ~5x speedup

EXAMPLES::

    sage: z = sage.functions.prime_pi.PrimePi()
    sage: loads(dumps(z))
    prime_pi
    sage: loads(dumps(z)) == z
    True
"""

# ****************************************************************************
#       Copyright (C) 2009,2011 R. Andrew Ohana <andrew.ohana@gmail.com>
#       Copyright (C) 2009 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cypari2.paridecl cimport *
from cysignals.signals cimport *
from cysignals.memory cimport sig_malloc, sig_realloc, sig_free

from libc.stdint cimport int_fast8_t, uint_fast16_t, uint8_t, uint32_t, uint64_t
from sage.rings.integer cimport Integer
from sage.libs.pari.all import pari
from sage.symbolic.function cimport BuiltinFunction
from sage.libs.gmp.mpz cimport *


cdef uint64_t arg_to_uint64(x, str s1, str s2) except -1:
    if not isinstance(x, Integer):
        from .other import floor
        x = Integer(floor(x))
    if mpz_sgn((<Integer>x).value) <= 0:
        return 0ull
    if mpz_sizeinbase((<Integer>x).value, 2) > 63:
        raise NotImplementedError("computation of " + s1 + " for x >= "
                + "2^63 is not implemented")
    if mpz_sizeinbase((<Integer>x).value, 2) > 43:
        import warnings
        warnings.warn("computation of %s for large x can take minutes, "%s1
                + "hours, or days depending on the size of %s"%s2)
        if mpz_sizeinbase((<Integer>x).value, 2) > 50:
            warnings.warn("computation of %s for x >= 2^50 has not "%s1
                    + "been as thoroughly tested as for smaller values")
    cdef uint64_t ret = mpz_get_ui((<Integer>x).value) & 0xfffffffful
    ret += (<uint64_t>mpz_get_ui((<Integer>(x>>32)).value)) << 32ull
    return ret

cdef class PrimePi(BuiltinFunction):
    def __init__(self):
        r"""
        The prime counting function, which counts the number of primes less
        than or equal to a given value.

        INPUT:

        - ``x`` - a real number
        - ``prime_bound`` - (default 0) a real number < 2^32, ``prime_pi`` will
          make sure to use all the primes up to ``prime_bound`` (although,
          possibly more) in computing ``prime_pi``, this can potentially
          speedup the time of computation, at a cost to memory usage.

        OUTPUT:

        integer -- the number of primes :math:`\leq` ``x``

        EXAMPLES:

        These examples test common inputs::

            sage: prime_pi(7)
            4
            sage: prime_pi(100)
            25
            sage: prime_pi(1000)
            168
            sage: prime_pi(100000)
            9592
            sage: prime_pi(500509)
            41581

        These examples test a variety of odd inputs::

            sage: prime_pi(3.5)
            2
            sage: prime_pi(sqrt(2357))
            15
            sage: prime_pi(mod(30957, 9750979))
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce arguments: positive characteristic not allowed in symbolic computations

        We test non-trivial ``prime_bound`` values::

            sage: prime_pi(100000, 10000)
            9592
            sage: prime_pi(500509, 50051)
            41581

        The following test is to verify that :trac:`4670` has been essentially
        resolved::

            sage: prime_pi(10^10)
            455052511

        The ``prime_pi`` function also has a special plotting method, so it
        plots quickly and perfectly as a step function::

            sage: P = plot(prime_pi, 50, 100)

        NOTES:

        Uses a recursive implementation, using the optimizations described in
        [Oha2011]_.

        AUTHOR:

        - \R. Andrew Ohana (2011)
        """
        super(PrimePi, self).__init__('prime_pi', latex_name=r"\pi",
                conversions={'mathematica':'PrimePi', 'pari':'primepi',
                    'sympy':'primepi'})

    cdef uint32_t *__primes
    cdef uint32_t __numPrimes, __maxSieve, __primeBound
    cdef int_fast8_t *__tabS
    cdef uint_fast16_t *__smallPi
    cdef byteptr __pariPrimePtr

    def __dealloc__(self):
        if self.__smallPi != NULL:
            sig_free(self.__smallPi)
            sig_free(self.__tabS)

    cdef void _init_tables(self):
        pari.init_primes(0xffffu)
        self.__pariPrimePtr = diffptr
        self.__smallPi = <uint_fast16_t *>sig_malloc(
                0x10000u * sizeof(uint_fast16_t))
        cdef uint32_t p=0u, i=0u, k=0u
        while i < 0xfff1u: # 0xfff1 is the last prime up to 0xffff
            NEXT_PRIME_VIADIFF(p, self.__pariPrimePtr)
            while i < p:
                self.__smallPi[i] = k
                i += 1u
            k += 1u
        while i <= 0xffffu:
            self.__smallPi[i] = k
            i += 1u

        self.__tabS = <int_fast8_t *>sig_malloc(2310*sizeof(int_fast8_t))
        for i in range(2310u):
            self.__tabS[i] = ((i+1u)/2u - (i+3u)/6u - (i+5u)/10u + (i+15u)/30u
                    - (i+7u)/14u + (i+21u)/42u + (i+35u)/70u - (i+105u)/210u
                    - (i+11u)/22u + (i+33u)/66u + (i+55u)/110u + (i+77u)/154u
                    - (i+165u)/330u - (i+231u)/462u - (i+385u)/770u
                    + (i+1155u)/2310u - ((i/77u)<<4u))


    def __call__(self, *args, coerce=True, hold=False):
        r"""
        EXAMPLES::

            sage: prime_pi.__call__(756)
            133
            sage: prime_pi.__call__(6574, 577)
            850
            sage: f(x) = prime_pi.__call__(x^2); f(x)
            prime_pi(x^2)
            sage: f(5)
            9
            sage: prime_pi.__call__(1, 2, 3)
            Traceback (most recent call last):
            ...
            TypeError: Symbolic function prime_pi takes 1 or 2 arguments (3 given)
        """
        if len(args) > 2:
            raise TypeError("Symbolic function %s takes 1 or 2"%self._name
                    + " arguments (%s given)"%len(args))
        else:
            self.__primeBound = 0u if len(args) < 2 else args[1]
        return super(PrimePi, self).__call__(args[0], coerce=coerce, hold=hold)

    def _eval_(self, x):
        r"""
        EXAMPLES::

            sage: prime_pi._eval_(7)
            4
            sage: prime_pi._eval_(100)
            25
            sage: prime_pi._eval_(1000)
            168
            sage: prime_pi._eval_(100000)
            9592
            sage: prime_pi._eval_(500509)
            41581
            sage: prime_pi._eval_(3.5)
            2
            sage: prime_pi._eval_(sqrt(2357))
            15
            sage: prime_pi._eval_(str(-2^100))
            0
            sage: prime_pi._eval_(mod(30957, 9750979))
            3337

        Make sure we actually compute correct results for 64-bit entries::

            sage: for i in (32..42): prime_pi(2^i) # long time (13s on sage.math, 2011)
            203280221
            393615806
            762939111
            1480206279
            2874398515
            5586502348
            10866266172
            21151907950
            41203088796
            80316571436
            156661034233

        This implementation uses unsigned 64-bit ints and does not support
        :math:`x \geq 2^63`::

            sage: prime_pi(2^63)
            Traceback (most recent call last):
            ...
            NotImplementedError: computation of prime_pi for x >= 2^63 is not implemented
        """
        cdef uint64_t z
        try:
            z = arg_to_uint64(x, 'prime_pi', 'x')
        except NotImplementedError:
            raise
        except TypeError:
            return None
        if self.__smallPi == NULL:
            self._init_tables()
        z = self._pi(z, self.__primeBound)
        self._clean_cache()
        return Integer(z)

    cdef uint64_t _pi(self, uint64_t x, uint64_t b) except -1:
        r"""
        Returns pi(x) under the assumption that 0 <= x < 2^64
        """
        if x <= 0xffffull: return self.__smallPi[x]
        if b*b < x:
            b = Integer(x).sqrtrem()[0]
        elif b > x:
            b = x
        self._init_primes(b)
        if not sig_on_no_except():
            self._clean_cache()
            cython_check_exception()
        b = self.__numPrimes
        b += self._phi(x, b)-1ull
        sig_off()
        return b

    cdef uint32_t _cached_count(self, uint32_t p):
        r"""
        For p < 65536, returns the value stored in ``self.__smallPi[p]``. For
        p <= ``self.__maxSieve``, uses a binary search on ``self.__primes`` to
        compute pi(p).
        """
        # inspired by Yann Laigle-Chapuy's suggestion
        if p <= 0xffffu: return self.__smallPi[p]
        cdef uint32_t size = (self.__numPrimes)>>1u
        # Use the expected density of primes for expected inputs to make an
        # educated guess
        if p>>3u < size:
            size = p>>3u
        # deal with edge case separately
        elif p >= self.__primes[self.__numPrimes-1u]:
            return self.__numPrimes
        cdef uint32_t pos = size
        cdef uint32_t prime
        while size:
            prime = self.__primes[pos]
            size >>= 1u
            if prime < p: pos += size
            elif prime > p: pos -= size
            else: return pos+1u
        if self.__primes[pos] <= p:
            while self.__primes[pos] <= p: pos += 1u
            return pos
        while self.__primes[pos] > p: pos -= 1u
        return pos+1u

    cdef void _clean_cache(self):
        if self.__numPrimes:
            sig_free(self.__primes)
            self.__numPrimes = 0u
            self.__maxSieve = 0u

    cdef uint64_t _init_primes(self, uint32_t b) except -1:
        """
        Populates ``self.__primes`` with all primes < b
        """
        cdef uint32_t *prime
        cdef uint32_t newNumPrimes, i
        pari.init_primes(b+1u)
        self.__pariPrimePtr = diffptr
        newNumPrimes = self._pi(b, 0ull)
        if self.__numPrimes:
            prime = <uint32_t *>sig_realloc(self.__primes,
                    newNumPrimes * sizeof(uint32_t))
        else:
            prime = <uint32_t *>sig_malloc(newNumPrimes*sizeof(uint32_t))
        if not sig_on_no_except():
            self.__numPrimes = newNumPrimes
            self._clean_cache()
            cython_check_exception()
        if prime == NULL:
            raise RuntimeError("not enough memory, maybe try with a smaller "
                    + "prime_bound?")
        self.__primes = prime
        prime += self.__numPrimes
        for i in range(self.__numPrimes, newNumPrimes):
            prime[0] = 0u if prime == self.__primes else prime[-1]
            NEXT_PRIME_VIADIFF(prime[0], self.__pariPrimePtr)
            prime += 1
        self.__numPrimes = newNumPrimes
        self.__maxSieve = b
        sig_off()

    cdef uint64_t _phi(self, uint64_t x, uint64_t i):
        r"""
        Legendre's formula: returns the number of primes :math:`\leq` ``x``
        that are not divisible by the first ``i`` primes
        """
        if not i: return x
        # explicitly compute for small i
        cdef uint64_t s = (x+1ull)>>1ull
        if i == 1ull: return s
        s -= (x+3ull)/6ull
        if i == 2ull: return s
        s -= (x+5ull)/10ull - (x+15ull)/30ull
        if i == 3ull: return s
        s -= ((x+7ull)/14ull - (x+21ull)/42ull - (x+35ull)/70ull +
                (x+105ull)/210ull)
        if i == 4ull: return s
        s -= ((x+11ull)/22ull - (x+33ull)/66ull - (x+55ull)/110ull -
                (x+77ull)/154ull + (x+165ull)/330ull + (x+231ull)/462ull +
                (x+385ull)/770ull - (x+1155ull)/2310ull)
        if i == 5ull: return s
        cdef uint64_t y=x/13ull, j=5ull
        cdef uint32_t *prime=self.__primes+5
        # switch to 32-bit as quickly as possible
        while y > 0xffffffffull:
            s -= self._phi(y, j)
            j += 1ull
            if j == i: return s
            prime += 1
            y = x/(<uint64_t>prime[0])
        # get y <= maxSieve so we can use a binary search with our table of
        # primes
        while y > (<uint64_t>self.__maxSieve):
            s -= self._phi32(y, j)
            j += 1ull
            if j == i: return s
            prime += 1
            y = x/(<uint64_t>prime[0])
        cdef uint64_t prime2 = prime[-1]
        # get p^2 > y so that we can use the identity phi(x,a)=pi(x)-a+1
        while prime2*prime2 <= y:
            s -= self._phi32(y, j)
            j += 1ull
            if j == i: return s
            prime2 = prime[0]
            prime += 1
            y = x/(<uint64_t>prime[0])
        s += j
        # use the identity phi(x,a) = pi(x)-a+1 and compute pi using a binary
        # search
        while prime2 < y:
            s -= self._cached_count(y)-j
            j += 1ull
            if j == i: return s-i
            prime2 = prime[0]
            prime += 1
            y = x/(<uint64_t>prime[0])
        return s-i

    cdef uint32_t _phi32(self, uint32_t x, uint32_t i):
        """
        Same as _phi except specialized for 32-bit ints
        """
        # table method for explicit computation was suggested by Yann
        # Laigle-Chapuy
        if i == 5u: return ((x/77u)<<4u) + self.__tabS[x%2310u]
        cdef uint32_t s = ((x/77u)<<4u) + self.__tabS[x%2310u]
        cdef uint32_t y = x/13u, j = 5u
        cdef uint32_t *prime = self.__primes+5
        while y > self.__maxSieve:
            s -= self._phi32(y, j)
            j += 1u
            if j == i: return s
            prime += 1
            y = x/prime[0]
        cdef uint32_t prime2 = prime[-1]
        while prime2*prime2 <= y:
            s -= self._phi32(y, j)
            j += 1u
            if j == i: return s
            prime2 = prime[0]
            prime += 1
            y = x/prime[0]
        s += j
        while 0xffffu < y:
            s -= self._cached_count(y)-j
            j += 1u
            if j == i: return s-i
            prime += 1
            y = x/prime[0]
        while prime2 < y:
            s -= self.__smallPi[y]-j
            j += 1u
            if j == i: return s-i
            prime2 = prime[0]
            prime += 1
            y = x/prime[0]
        return s-i

    def plot(self, xmin=0, xmax=100, vertical_lines=True, **kwds):
        """
        Draw a plot of the prime counting function from ``xmin`` to ``xmax``.
        All additional arguments are passed on to the line command.

        WARNING: we draw the plot of ``prime_pi`` as a stairstep function with
        explicitly drawn vertical lines where the function jumps. Technically
        there should not be any vertical lines, but they make the graph look
        much better, so we include them. Use the option ``vertical_lines=False``
        to turn these off.

        EXAMPLES::

            sage: plot(prime_pi, 1, 100)
            Graphics object consisting of 1 graphics primitive
            sage: prime_pi.plot(-2, sqrt(2501), thickness=2, vertical_lines=False)
            Graphics object consisting of 16 graphics primitives
        """
        from sage.plot.step import plot_step_function
        if xmax < xmin:
            return plot_step_function([], **kwds)
        if xmax < 2:
            return plot_step_function([(xmin,0),(xmax,0)], **kwds)
        y = self(xmin)
        v = [(xmin, y)]
        from sage.rings.all import prime_range
        for p in prime_range(xmin+1, xmax+1, py_ints=True):
            y += 1
            v.append((p,y))
        v.append((xmax,y))
        return plot_step_function(v, vertical_lines=vertical_lines, **kwds)

########
prime_pi = PrimePi()

cpdef Integer legendre_phi(x, a):
    r"""
    Legendre's formula, also known as the partial sieve function, is a useful
    combinatorial function for computing the prime counting function (the
    ``prime_pi`` method in Sage). It counts the number of positive integers
    :math:`\leq` ``x`` that are not divisible by the first ``a`` primes.

    INPUT:

    - ``x`` -- a real number

    - ``a`` -- a non-negative integer

    OUTPUT:

    integer -- the number of positive integers :math:`\leq` ``x`` that are not
    divisible by the first ``a`` primes

    EXAMPLES::

        sage: legendre_phi(100, 0)
        100
        sage: legendre_phi(29375, 1)
        14688
        sage: legendre_phi(91753, 5973)
        2893
        sage: legendre_phi(7.5, 2)
        3
        sage: legendre_phi(str(-2^100), 92372)
        0
        sage: legendre_phi(4215701455, 6450023226)
        1

    NOTES:

    Uses a recursive implementation, using the optimizations described in
    [Oha2011]_.

    AUTHOR:

    - \R. Andrew Ohana (2011)
    """
    if not isinstance(a, Integer):
        a = Integer(a)
    if a < Integer(0):
        raise ValueError("a (=%s) must be non-negative"%a)
    cdef uint64_t y = arg_to_uint64(x, 'legendre_phi', 'x and a')

    # legendre_phi(x, a) = 0 when x <= 0
    if not y: return Integer(0)

    # legendre_phi(x, 0) = x
    if a == Integer(0): return Integer(y)

    # Use knowledge about the density of primes to quickly compute for many
    # cases where a is unusually large
    if a > Integer(y>>1ull): return Integer(1)
    if y > 1916ull:
        chk = Integer(y)*Integer(13271040)//Integer(86822723)
        if a > chk: return Integer(1)

    # If a > prime_pi(2^32), we compute phi(x,a) = max(pi(x)-a+1,1)
    if a > Integer(203280221):
        ret = prime_pi(x)-a+Integer(1)
        if ret < Integer(1): return Integer(1)
        return ret

    # Deal with the general case
    if (<PrimePi>prime_pi).__smallPi == NULL:
        (<PrimePi>prime_pi)._init_tables()
    cdef uint32_t z = pari.prime(a)
    if z >= y: return Integer(1)
    (<PrimePi>prime_pi)._init_primes(z)
    if not sig_on_no_except():
        (<PrimePi>prime_pi)._clean_cache()
        cython_check_exception()
    y = (<PrimePi>prime_pi)._phi(y, mpz_get_ui((<Integer>a).value))
    sig_off()
    (<PrimePi>prime_pi)._clean_cache()
    return Integer(y)

partial_sieve_function = legendre_phi
