# sage.doctest: needs sage.libs.gap sage.libs.flint sage.libs.pari sage.modules
r"""
Functions that compute some of the sequences in Sloane's tables

EXAMPLES:

Type ``sloane.[tab]`` to see a list of the sequences that are defined.

::

    sage: a = sloane.A000005; a
    The integer sequence tau(n), which is the number of divisors of n.
    sage: a(1)
    1
    sage: a(6)
    4
    sage: a(100)
    9

Type ``d._eval??`` to see how the function that
computes an individual term of the sequence is implemented.

The input must be a positive integer::

    sage: a(0)
    Traceback (most recent call last):
    ...
    ValueError: input n (=0) must be a positive integer
    sage: a(1/3)
    Traceback (most recent call last):
    ...
    TypeError: input must be an int or Integer

You can also change how a sequence prints::

    sage: a = sloane.A000005; a
    The integer sequence tau(n), which is the number of divisors of n.
    sage: a.rename('(..., tau(n), ...)')
    sage: a
    (..., tau(n), ...)
    sage: a.reset_name()
    sage: a
    The integer sequence tau(n), which is the number of divisors of n.

TESTS::

    sage: a = sloane.A000001
    sage: a == loads(dumps(a))
    True

We agree with the online database::

    sage: for t in sloane.__dir__():    # long time; optional -- internet; known bug
    ....:     online_list = list(oeis(t).first_terms())
    ....:     L = max(2, len(online_list) // 2)
    ....:     sage_list = sloane.__getattribute__(t).list(L)
    ....:     if online_list[:L] != sage_list:
    ....:         print('{} seems wrong'.format(t))

.. SEEALSO::

    - If you want to get more informations relative to a sequence (references,
      links, examples, programs, ...), you can use the On-Line Encyclopedia of
      Integer Sequences provided by the :mod:`OEIS <sage.databases.oeis>`
      module.
    - If you plan to do a lot of automatic searches for subsequences, you
      should consider installing :mod:`SloaneEncyclopedia
      <sage.databases.sloane>`, a local partial copy of the OEIS.


AUTHORS:

- William Stein: framework

- Jaap Spies: most sequences

- Nick Alexander: updated framework
"""

########################################################################
#
# To add your own new sequence here, do the following:
#
# 1. Add a new class to Section II below, which you should
#    do by copying an existing class and modifying it.
#    Make sure to at least define _eval and _repr_.
#    NOTE:   (a) define the _eval method only, which you may
#                assume has as input a *positive* Sage integer (offset > 0).
#                Each sequence in the OEIS has an offset >= 0, indicating the
#                value of the first index. The default offset = 1.
#            (b) define the list method if there is a faster
#                way to compute the terms of the sequence than
#                just calling _eval (which is the default definition
#                of list, note: the offset is counted for, it lists n numbers).
#            (c) *AVOID* using gp.method if possible!  Use pari(obj).method()
#            (d) In many cases the function that computes a given integer
#                sequence belongs elsewhere in Sage.  Put it there and make
#                your class in this file just call it.
#            (e) _eval should always return a Sage integer.
#
# 2. Add an instance of your class in Section III below.

#
# 3. Type "sage -br" to rebuild Sage, then fire up the notebook and
#    try out your new sequence.  Click the text button to get a version
#    of your session that you then include as a docstring.
#    You can check your results with the entries of the OEIS::

#       sage: seq = oeis(45) ; seq
#       A000045: Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.
#       sage: seq.first_terms()[:12]
#       (0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89)
#
# 4. Send a patch using
#      sage: hg_sage.ci()
#      sage: hg_sage.send('patchname')
#    (Email it to sage-dev@groups.google.com or post it online.)
#
########################################################################

########################################################################
# I. Define the generic Sloane sequence class.
########################################################################

# just used for handy .load, .save, etc.
import sys
import inspect

from sage.structure.sage_object import SageObject
from sage.arith.srange import srange
from sage.rings.integer_ring import ZZ
from sage.misc.lazy_import import lazy_import
from sage.rings.integer import Integer as Integer_class
# You may have to import more here when defining new sequences
import sage.arith.all as arith
from sage.rings.rational_field import QQ
from sage.combinat import combinat
from sage.misc.misc_c import prod

lazy_import("sage.functions.all", "prime_pi")
lazy_import('sage.matrix.matrix_space', 'MatrixSpace')


class SloaneSequence(SageObject):
    r"""
    Base class for a Sloane integer sequence.
    """

    def __init__(self, offset=1):
        r"""
        A sequence starting at offset (=1 by default).

        EXAMPLES::

            sage: from sage.combinat.sloane_functions import SloaneSequence
            sage: SloaneSequence().offset
            1
            sage: SloaneSequence(4).offset
            4
        """
        self.offset = ZZ(offset)

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.combinat.sloane_functions import SloaneSequence
            sage: SloaneSequence(4)._repr_()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def __eq__(self, other):
        """
        EXAMPLES::

            sage: sloane.A000007 == sloane.A000045
            False
            sage: sloane.A000007 == sloane.A000007
            True
        """
        if not isinstance(other, SloaneSequence):
            return False
        return repr(self) == repr(other)

    def __ne__(self, other):
        """
        EXAMPLES::

            sage: sloane.A000007 != sloane.A000045
            True
            sage: sloane.A000007 != sloane.A000007
            False
        """
        return not self.__eq__(other)

    def _sage_src_(self):
        """
        Return the source code for the class of ``self``.

        EXAMPLES::

            sage: sloane.A000045._sage_src_()
            'class A000045(...'
        """
        from sage.misc.sageinspect import sage_getsource
        return sage_getsource(self.__class__)

    def __call__(self, n):
        """
        EXAMPLES::

            sage: sloane.A000007(2)
            0
            sage: sloane.A000007('a')
            Traceback (most recent call last):
            ...
            TypeError: input must be an int or Integer
            sage: sloane.A000007(-1)
            Traceback (most recent call last):
            ...
            ValueError: input n (=-1) must be an integer >= 0
            sage: sloane.A000001(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
        """
        if not isinstance(n, (int, Integer_class)):
            raise TypeError("input must be an int or Integer")
        m = ZZ(n)
        if m < self.offset:
            if self.offset == 1:
                raise ValueError("input n (=%s) must be a positive integer" % (n))
            else:
                raise ValueError("input n (=%s) must be an integer >= %s" % (n, self.offset))
        return self._eval(m)

    def _eval(self, n):
        """
        EXAMPLES::

            sage: from sage.combinat.sloane_functions import SloaneSequence
            sage: SloaneSequence(0)._eval(4)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        # this is what you implement in the derived class
        # the input n is assumed to be a *Sage* integer >= offset
        raise NotImplementedError

    def list(self, n):
        r"""
        Return ``n`` terms of the sequence::

            sequence[offset], sequence[offset+1], ..., sequence[offset+n-1].

        EXAMPLES::

            sage: sloane.A000012.list(4)
            [1, 1, 1, 1]
        """
        return [self._eval(i) for i in srange(self.offset, n+self.offset)]

    # The Python default tries repeated __getitem__ calls, which will succeed,
    # but is probably not what is wanted.
    # This prevents list(sequence) from wandering off.
    def __iter__(self):
        """
        EXAMPLES::

            sage: iter(sloane.A000012)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def __getitem__(self, n):
        r"""
        Return the ``n``-th item of ``self``.

        We interpret slices as best we can, but our sequences are infinite
        so we want to prevent some mis-incantations.

        Therefore, we arbitrarily cap slices to be at most LENGTH=100000
        elements long. Since many Sloane sequences are costly to compute,
        this is probably not an unreasonable decision, but just in case,
        list does not cap length.

        EXAMPLES::

            sage: sloane.A000012[3]
            1
            sage: sloane.A000012[:4]
            [1, 1, 1, 1]
            sage: sloane.A000012[:10]
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
            sage: sloane.A000012[4:10]
            [1, 1, 1, 1, 1, 1]
            sage: sloane.A000012[0:1000000000]
            Traceback (most recent call last):
            ...
            IndexError: slice (=slice(0, 1000000000, None)) too long

            sage: sloane.A000001[:8]
            [1, 1, 1, 2, 1, 2, 1]
            sage: sloane.A000001[0:8:2]
            [1, 2, 2]
            sage: sloane.A000001[1:8:2]
            [1, 1, 1, 1]
        """
        if not isinstance(n, slice):
            return self(n)

        LENGTH = 100000
        (start, stop, step) = n.indices(2 * LENGTH)
        if abs(stop - start) > LENGTH:
            raise IndexError("slice (=%s) too long" % n)
        return [self(i) for i in range(start, stop, step) if i >= self.offset]


########################################################################
# II. Actual implementations of Sloane sequences.
########################################################################


# This one should be here!
class A000001(SloaneSequence):
    def __init__(self):
        r"""
        Number of groups of order `n`.

        INPUT:

        - ``n`` -- positive integer

        OUTPUT: integer

        EXAMPLES::

            sage: a = sloane.A000001;a
            Number of groups of order n.
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(1)
            1
            sage: a(2)
            1
            sage: a(9)
            2
            sage: a.list(16)
            [1, 1, 1, 2, 1, 2, 1, 5, 2, 2, 1, 5, 1, 2, 1, 14]
            sage: a(60)
            13

        AUTHORS:

        - Jaap Spies (2007-02-04)
        """
        self._small = [1, 1, 1, 2, 1, 2, 1, 5, 2, 2, 1, 5, 1, 2, 1, 14, 1, 5, 1, 5, 2, 2, 1, 15, 2, 2, 5, 4, 1, 4, 1, 51, 1, 2, 1, 14, 1, 2, 2, 14, 1, 6, 1, 4, 2, 2, 1, 52, 2, 5]
        SloaneSequence.__init__(self, offset=1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000001._repr_()
            'Number of groups of order n.'
        """
        return "Number of groups of order n."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: sloane.A000001._eval(4)
            2
            sage: sloane.A000001._eval(51)
            1
            sage: sloane.A000001._eval(5000)
            Traceback (most recent call last):
            ...
            GAPError: Error, the library of groups of size 5000 is not available
        """
        if n <= 50:
            return self._small[n - 1]
        from sage.libs.gap.libgap import libgap
        return ZZ(libgap.NumberSmallGroups(n))


class A000027(SloaneSequence):
    def __init__(self):
        r"""
        The natural numbers. Also called the whole numbers, the counting
        numbers or the positive integers.

        The following examples are tests of SloaneSequence more than
        A000027.

        EXAMPLES::

            sage: s = sloane.A000027; s
            The natural numbers.
            sage: s(10)
            10

        Index n is interpreted as _eval(n)::

            sage: s[10]
            10

        Slices are interpreted with absolute offsets, so the following
        returns the terms of the sequence up to but not including the third
        term::

            sage: s[:3]
            [1, 2]
            sage: s[3:6]
            [3, 4, 5]
            sage: s.list(5)
            [1, 2, 3, 4, 5]
        """
        SloaneSequence.__init__(self, offset=1)

# is this a good idea to have a link for all sequences? Jaap
    link = "http://oeis.org/classic/A000027"

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000027._repr_()
            'The natural numbers.'
        """
        return "The natural numbers."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: sloane.A000027._eval(5)
            5
        """
        return n


class A000004(SloaneSequence):
    def __init__(self):
        r"""
        The zero sequence.

        INPUT:

        - ``n`` -- nonnegative integer

        EXAMPLES::

            sage: a = sloane.A000004; a
            The zero sequence.
            sage: a(1)
            0
            sage: a(2007)
            0
            sage: a.list(12)
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

        AUTHORS:

        - Jaap Spies (2006-12-10)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000004._repr_()
            'The zero sequence.'
        """
        return "The zero sequence."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: sloane.A000004._eval(5)
            0
        """
        return ZZ.zero()


class A000005(SloaneSequence):
    def __init__(self):
        r"""
        The sequence `tau(n)`, which is the number of divisors of `n`.

        This sequence is also denoted `d(n)` (also called
        `\tau(n)` or `\sigma_0(n)`), the number of
        divisors of `n`.

        INPUT:

        - ``n`` -- positive integer

        EXAMPLES::

            sage: d = sloane.A000005; d
            The integer sequence tau(n), which is the number of divisors of n.
            sage: d(1)
            1
            sage: d(6)
            4
            sage: d(51)
            4
            sage: d(100)
            9
            sage: d(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: d.list(10)
            [1, 2, 2, 3, 2, 4, 2, 4, 3, 4]

        AUTHORS:

        - Jaap Spies (2006-12-10)

        - William Stein (2007-01-08)
        """
        SloaneSequence.__init__(self, offset=1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000005._repr_()
            'The integer sequence tau(n), which is the number of divisors of n.'
        """
        return "The integer sequence tau(n), which is the number of divisors of n."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: sloane.A000005._eval(5)
            2
        """
        return arith.number_of_divisors(n)


class A000008(SloaneSequence):
    def __init__(self):
        r"""
        Number of ways of making change for n cents using coins
        of 1, 2, 5, 10 cents.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000008;a
            Number of ways of making change for n cents using coins of 1, 2, 5, 10 cents.
            sage: a(0)
            1
            sage: a(1)
            1
            sage: a(13)
            16
            sage: a.list(14)
            [1, 1, 2, 2, 3, 4, 5, 6, 7, 8, 11, 12, 15, 16]

        AUTHOR:

        - J. Gaski (2009-05-29)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000008._repr_()
            'Number of ways of making change for n cents using coins of 1, 2, 5, 10 cents.'
        """
        return "Number of ways of making change for n cents using coins of 1, 2, 5, 10 cents."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000008._eval(n) for n in range(14)]
            [1, 1, 2, 2, 3, 4, 5, 6, 7, 8, 11, 12, 15, 16]
        """
        from sage.combinat.partition import Partitions
        return Partitions(n, parts_in=[1, 2, 5, 10]).cardinality()


class A000009(SloaneSequence):
    def __init__(self):
        r"""
        Number of partitions of `n` into odd parts.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000009;a
            Number of partitions of n into odd parts.
            sage: a(0)
            1
            sage: a(1)
            1
            sage: a(13)
            18
            sage: a.list(14)
            [1, 1, 1, 2, 2, 3, 4, 5, 6, 8, 10, 12, 15, 18]

        AUTHOR:

        - Jaap Spies (2007-01-30)
        """
        SloaneSequence.__init__(self, offset=0)
        self._b = []
        self._precompute(2)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000009._repr_()
            'Number of partitions of n into odd parts.'
        """
        return "Number of partitions of n into odd parts."

    def cf(self):
        """
        EXAMPLES::

            sage: it = sloane.A000009.cf()
            sage: [next(it) for i in range(14)]
            [1, 1, 1, 2, 2, 3, 4, 5, 6, 8, 10, 12, 15, 18]
        """
        _, x = QQ['x'].objgen()
        k = 0
        yield ZZ.one()
        p = 1
        while True:
            k += 1
            p *= (1+x**k)
            yield ZZ(p.coefficients(sparse=False)[k])

    def _precompute(self, how_many=50):
        """
        EXAMPLES::

            sage: initial = len(sloane.A000009._b)
            sage: sloane.A000009._precompute(10)
            sage: len(sloane.A000009._b) - initial == 10
            True
        """
        try:
            f = self._f
        except AttributeError:
            self._f = self.cf()
            f = self._f
        self._b += [next(f) for i in range(how_many)]

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000009._eval(i) for i in range(14)]
            [1, 1, 1, 2, 2, 3, 4, 5, 6, 8, 10, 12, 15, 18]
        """
        if len(self._b) <= n:
            self._precompute(n - len(self._b) + 1)
        return self._b[n]

    def list(self, n):
        """
        EXAMPLES::

            sage: sloane.A000009.list(14)
            [1, 1, 1, 2, 2, 3, 4, 5, 6, 8, 10, 12, 15, 18]
        """
        self._eval(n)   # force computation
        return self._b[:n]


class A000796(SloaneSequence):
    def __init__(self):
        r"""
        Decimal expansion of `\pi`.

        INPUT:

        - ``n`` -- positive integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000796;a
            Decimal expansion of Pi.
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(1)
            3
            sage: a(13)
            9
            sage: a.list(14)
            [3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9, 7]
            sage: a(100)
            7

        AUTHOR:

        - Jaap Spies (2007-01-30)
        """
        SloaneSequence.__init__(self, offset=1)
        self._b = []

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000796._repr_()
            'Decimal expansion of Pi.'
        """
        return "Decimal expansion of Pi."

    def pi(self):
        """
        Based on an algorithm of Lambert Meertens The ABC-programming
        language!!!

        EXAMPLES::

            sage: it = sloane.A000796.pi()
            sage: [next(it) for i in range(10)]
            [3, 1, 4, 1, 5, 9, 2, 6, 5, 3]
        """
        k, a, b, a1, b1 = ZZ(2), ZZ(4), ZZ.one(), ZZ(12), ZZ(4)
        while True:
            p, q, k = k*k, 2*k+1, k+1
            a, b, a1, b1 = a1, b1, p*a+q*a1, p*b+q*b1
            d, d1 = a//b, a1//b1
            while d == d1:
                yield d
                a, a1 = 10*(a % b), 10*(a1 % b1)
                d, d1 = a//b, a1//b1

    def _precompute(self, how_many=1000):
        """
        EXAMPLES::

            sage: initial = len(sloane.A000796._b)
            sage: sloane.A000796._precompute(10)
            sage: len(sloane.A000796._b) - initial
            10
        """
        try:
            f = self._f
        except AttributeError:
            self._f = self.pi()
            f = self._f
        self._b += [next(f) for i in range(how_many)]

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000796._eval(n) for n in range(1,11)]
            [3, 1, 4, 1, 5, 9, 2, 6, 5, 3]
        """
        while len(self._b) <= n:
            self._precompute()
        return self._b[n-1]

    def list(self, n):
        """
        EXAMPLES::

            sage: sloane.A000796.list(10)
            [3, 1, 4, 1, 5, 9, 2, 6, 5, 3]
        """
        self._eval(n)   # force computation
        return self._b[:n]


class A003418(SloaneSequence):
    def __init__(self):
        r"""
        Least common multiple (or lcm) of `\{1, 2, \ldots, n\}`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A003418;a
            Least common multiple (or lcm) of {1, 2, ..., n}.
            sage: a(0)
            1
            sage: a(1)
            1
            sage: a(13)
            360360
            sage: a.list(14)
            [1, 1, 2, 6, 12, 60, 60, 420, 840, 2520, 2520, 27720, 27720, 360360]
            sage: a(20.0)
            Traceback (most recent call last):
            ...
            TypeError: input must be an int or Integer

        AUTHOR:

        - Jaap Spies (2007-01-31)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A003418._repr_()
            'Least common multiple (or lcm) of {1, 2, ..., n}.'
        """
        return "Least common multiple (or lcm) of {1, 2, ..., n}."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A003418._eval(n) for n in range(1,11)]
            [1, 2, 6, 12, 60, 60, 420, 840, 2520, 2520]
        """
        return arith.lcm(range(1, n + 1))


class A007318(SloaneSequence):
    def __init__(self):
        r"""
        Pascal's triangle read by rows:
        `C(n,k) = \binom{n}{k} = \frac {n!} {(k!(n-k)!)}`,
        `0 \le k \le n`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A007318
            sage: a(0)
            1
            sage: a(1)
            1
            sage: a(13)
            4
            sage: a.list(15)
            [1, 1, 1, 1, 2, 1, 1, 3, 3, 1, 1, 4, 6, 4, 1]
            sage: a(100)
            715

        AUTHORS:

        - Jaap Spies (2007-01-31)
        """
        SloaneSequence.__init__(self, offset=0)

    keyword = ["nonn", "tabl", "nice", "easy", "core", "triangle"]

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A007318._repr_()
            "Pascal's triangle read by rows: C(n,k) = binomial(n,k) = n!/(k!*(n-k)!), 0<=k<=n."
        """
        return "Pascal's triangle read by rows: C(n,k) = binomial(n,k) = n!/(k!*(n-k)!), 0<=k<=n."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A007318._eval(n) for n in range(10)]
            [1, 1, 1, 1, 2, 1, 1, 3, 3, 1]
        """
        m = 0
        while m*(m+1)//2 <= n:
            m += 1
        m -= 1
        k = n - m*(m+1)//2
        return arith.binomial(m, k)


class A008275(SloaneSequence):
    def __init__(self):
        r"""
        Triangle of Stirling numbers of first kind, `s(n,k)`,
        `n \ge 1`, `1 \le k \le n`.

        The unsigned numbers are also called Stirling cycle numbers:

        `|s(n,k)|` = number of permutations of `n` objects
        with exactly `k` cycles.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A008275;a
            Triangle of Stirling numbers of first kind, s(n,k), n >= 1, 1<=k<=n.
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(1)
            1
            sage: a(2)
            -1
            sage: a(3)
            1
            sage: a(11)
            24
            sage: a.list(12)
            [1, -1, 1, 2, -3, 1, -6, 11, -6, 1, 24, -50]

        AUTHORS:

        - Jaap Spies (2007-02-02)
        """
        SloaneSequence.__init__(self, offset=1)

    keyword = ["sign", "tabl", "nice", "core", "triangle"]

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A008275._repr_()
            'Triangle of Stirling numbers of first kind, s(n,k), n >= 1, 1<=k<=n.'
        """
        return "Triangle of Stirling numbers of first kind, s(n,k), n >= 1, 1<=k<=n."

    def s(self, n, k):
        """
        EXAMPLES::

            sage: sloane.A008275.s(4,2)
            11
            sage: sloane.A008275.s(5,2)
            -50
            sage: sloane.A008275.s(5,3)
            35
        """
        return (-1)**(n-k) * combinat.stirling_number1(n, k)

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A008275._eval(n) for n in range(1, 11)]
            [1, -1, 1, 2, -3, 1, -6, 11, -6, 1]
        """
        m = 0
        while m*(m+1)//2 < n:
            m += 1
        k = n - m*(m-1)//2
        return self.s(m, k)  # (-1)**(m-k) * combinat.stirling_number1(m,k)


class A008277(SloaneSequence):
    def __init__(self):
        r"""
        Triangle of Stirling numbers of 2nd kind, `S2(n,k)`,
        `n \ge 1`, `1 \le k \le n`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A008277;a
            Triangle of Stirling numbers of 2nd kind, S2(n,k), n >= 1, 1<=k<=n.
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(1)
            1
            sage: a(2)
            1
            sage: a(3)
            1
            sage: a(4.0)
            Traceback (most recent call last):
            ...
            TypeError: input must be an int or Integer
            sage: a.list(15)
            [1, 1, 1, 1, 3, 1, 1, 7, 6, 1, 1, 15, 25, 10, 1]

        AUTHORS:

        - Jaap Spies (2007-01-31)
        """
        SloaneSequence.__init__(self, offset=1)

    keyword = ["nonn", "tabl", "nice", "core", "triangle"]

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A008277._repr_()
            'Triangle of Stirling numbers of 2nd kind, S2(n,k), n >= 1, 1<=k<=n.'
        """
        return "Triangle of Stirling numbers of 2nd kind, S2(n,k), n >= 1, 1<=k<=n."

    def s2(self, n, k):
        """
        Return the Stirling number S2(n,k) of the 2nd kind.

        EXAMPLES::

            sage: sloane.A008277.s2(4,2)
            7
        """
        return combinat.stirling_number2(n, k)

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A008277._eval(n) for n in range(1,11)]
            [1, 1, 1, 1, 3, 1, 1, 7, 6, 1]
        """
        m = 0
        while m*(m+1)//2 < n:
            m += 1
        k = n - m*(m-1)//2
        return self.s2(m, k)  # combinat.stirling_number2(m,k)


class A049310(SloaneSequence):
    def __init__(self):
        r"""
        Triangle of coefficients of Chebyshev's `S(n,x)`:
        `U(n, \frac x 2)` polynomials (exponents in increasing
        order).

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A049310;a
            Triangle of coefficients of Chebyshev's S(n,x) := U(n,x/2) polynomials (exponents in increasing order).
            sage: a(0)
            1
            sage: a(1)
            0
            sage: a(13)
            0
            sage: a.list(15)
            [1, 0, 1, -1, 0, 1, 0, -2, 0, 1, 1, 0, -3, 0, 1]
            sage: a(200)
            0
            sage: a.keyword
            ['sign', 'tabl', 'nice', 'easy', 'core', 'triangle']

        AUTHORS:

        - Jaap Spies (2007-01-31)
        """
        SloaneSequence.__init__(self, offset=0)

    keyword = ["sign", "tabl", "nice", "easy", "core", "triangle"]

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A049310._repr_()
            "Triangle of coefficients of Chebyshev's S(n,x) := U(n,x/2) polynomials (exponents in increasing order)."
        """
        return "Triangle of coefficients of Chebyshev's S(n,x) := U(n,x/2) polynomials (exponents in increasing order)."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A049310._eval(n) for n in range(10)]
            [1, 0, 1, -1, 0, 1, 0, -2, 0, 1]
        """
        m = 0
        while m*(m+1)//2 <= n:
            m += 1
        m -= 1
        k = n - m*(m+1)//2
        if (m+k) % 2:
            return ZZ(0)
        sign = (-1)**((m+k)//2 + k)
        return sign * arith.binomial((m+k)//2, k)


class A000010(SloaneSequence):
    def __init__(self):
        r"""
        The integer sequence A000010 is Euler's totient function.

        Number of positive integers `i < n` that are relative prime
        to `n`. Number of totatives of `n`.

        Euler totient function `\phi(n)`: count numbers `n`
        and prime to `n`. ``euler_phi`` is a standard Sage function
        implemented in PARI.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000010; a
            Euler's totient function
            sage: a(1)
            1
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(11)
            10
            sage: a.list(12)
            [1, 1, 2, 2, 4, 2, 6, 4, 6, 4, 10, 4]
            sage: a(1/3)
            Traceback (most recent call last):
            ...
            TypeError: input must be an int or Integer

        AUTHORS:

        - Jaap Spies (2007-01-12)
        """
        SloaneSequence.__init__(self, offset=1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000010._repr_()
            "Euler's totient function"
        """
        return "Euler's totient function"

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000010._eval(n) for n in range(1,11)]
            [1, 1, 2, 2, 4, 2, 6, 4, 6, 4]
        """
        return arith.euler_phi(n)

# Theme: simple functions


class A000007(SloaneSequence):
    def __init__(self):
        r"""
        The characteristic function of 0: `a(n) = 0^n`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000007;a
            The characteristic function of 0: a(n) = 0^n.
            sage: a(0)
            1
            sage: a(2)
            0
            sage: a(12)
            0
            sage: a.list(12)
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

        AUTHORS:

        - Jaap Spies (2007-01-12)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000007._repr_()
            'The characteristic function of 0: a(n) = 0^n.'
        """
        return "The characteristic function of 0: a(n) = 0^n."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000007._eval(n) for n in range(10)]
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        """
        return ZZ(0**n)


class A005843(SloaneSequence):
    def __init__(self):
        r"""
        The even numbers: `a(n) = 2n`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A005843;a
            The even numbers: a(n) = 2n.
            sage: a(0.0)
            Traceback (most recent call last):
            ...
            TypeError: input must be an int or Integer
            sage: a(1)
            2
            sage: a(2)
            4
            sage: a(9)
            18
            sage: a.list(10)
            [0, 2, 4, 6, 8, 10, 12, 14, 16, 18]

        AUTHORS:

        - Jaap Spies (2007-02-03)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A005843._repr_()
            'The even numbers: a(n) = 2n.'
        """
        return "The even numbers: a(n) = 2n."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A005843._eval(n) for n in range(10)]
            [0, 2, 4, 6, 8, 10, 12, 14, 16, 18]
        """
        return ZZ(2*n)


class A000035(SloaneSequence):
    def __init__(self):
        r"""
        A simple periodic sequence.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000035;a
            A simple periodic sequence.
            sage: a(0.0)
            Traceback (most recent call last):
            ...
            TypeError: input must be an int or Integer
            sage: a(1)
            1
            sage: a(2)
            0
            sage: a(9)
            1
            sage: a.list(10)
            [0, 1, 0, 1, 0, 1, 0, 1, 0, 1]

        AUTHORS:

        - Jaap Spies (2007-02-02)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000035._repr_()
            'A simple periodic sequence.'
        """
        return "A simple periodic sequence."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000035._eval(n) for n in range(10)]
            [0, 1, 0, 1, 0, 1, 0, 1, 0, 1]
        """
        return ZZ(n % 2)


class A000169(SloaneSequence):
    def __init__(self):
        r"""
        Number of labeled rooted trees with `n` nodes:
        `n^{(n-1)}`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000169;a
            Number of labeled rooted trees with n nodes: n^(n-1).
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(1)
            1
            sage: a(2)
            2
            sage: a(10)
            1000000000
            sage: a.list(11)
            [1, 2, 9, 64, 625, 7776, 117649, 2097152, 43046721, 1000000000, 25937424601]

        AUTHORS:

        - Jaap Spies (2007-01-26)
        """
        SloaneSequence.__init__(self, offset=1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000169._repr_()
            'Number of labeled rooted trees with n nodes: n^(n-1).'
        """
        return "Number of labeled rooted trees with n nodes: n^(n-1)."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000169._eval(n) for n in range(1,11)]
            [1, 2, 9, 64, 625, 7776, 117649, 2097152, 43046721, 1000000000]
        """
        return ZZ(n**(n-1))


class A000272(SloaneSequence):
    def __init__(self):
        r"""
        Number of labeled rooted trees on `n` nodes: `n^{(n-2)}`.

        INPUT:

        - ``n`` -- integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000272;a
            Number of labeled rooted trees with n nodes: n^(n-2).
            sage: a(0)
            1
            sage: a(1)
            1
            sage: a(2)
            1
            sage: a(10)
            100000000
            sage: a.list(12)
            [1, 1, 1, 3, 16, 125, 1296, 16807, 262144, 4782969, 100000000, 2357947691]

        AUTHORS:

        - Jaap Spies (2007-01-26)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000272._repr_()
            'Number of labeled rooted trees with n nodes: n^(n-2).'
        """
        return "Number of labeled rooted trees with n nodes: n^(n-2)."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000272._eval(n) for n in range(1,11)]
            [1, 1, 3, 16, 125, 1296, 16807, 262144, 4782969, 100000000]
        """
        if n == 0:
            return ZZ.one()
        return ZZ(n) ** (ZZ(n) - 2)


class A000312(SloaneSequence):
    def __init__(self):
        r"""
        Number of labeled mappings from `n` points to themselves
        (endofunctions): `n^n`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000312;a
            Number of labeled mappings from n points to themselves (endofunctions): n^n.
            sage: a(-1)
            Traceback (most recent call last):
            ...
            ValueError: input n (=-1) must be an integer >= 0
            sage: a(0)
            1
            sage: a(1)
            1
            sage: a(9)
            387420489
            sage: a.list(11)
            [1, 1, 4, 27, 256, 3125, 46656, 823543, 16777216, 387420489, 10000000000]

        AUTHORS:

        - Jaap Spies (2007-01-26)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000312._repr_()
            'Number of labeled mappings from n points to themselves (endofunctions): n^n.'
        """
        return "Number of labeled mappings from n points to themselves (endofunctions): n^n."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000312._eval(n) for n in range(10)]
            [1, 1, 4, 27, 256, 3125, 46656, 823543, 16777216, 387420489]
        """
        if n == 0:
            return ZZ.one()
        else:
            return ZZ(n**n)


class A001477(SloaneSequence):
    def __init__(self):
        r"""
        The nonnegative integers.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A001477;a
            The nonnegative integers.
            sage: a(-1)
            Traceback (most recent call last):
            ...
            ValueError: input n (=-1) must be an integer >= 0
            sage: a(0)
            0
            sage: a(3382789)
            3382789
            sage: a(11)
            11
            sage: a.list(12)
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

        AUTHORS:

        - Jaap Spies (2007-01-26)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A001477._repr_()
            'The nonnegative integers.'
        """
        return "The nonnegative integers."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A001477._eval(n) for n in range(10)]
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        """
        return ZZ(n)


class A004526(SloaneSequence):
    def __init__(self):
        r"""
        The nonnegative integers repeated.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A004526;a
            The nonnegative integers repeated.
            sage: a(0)
            0
            sage: a(1)
            0
            sage: a(2)
            1
            sage: a(10)
            5
            sage: a.list(12)
            [0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5]

        AUTHORS:

        - Jaap Spies (2007-01-26)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A004526._repr_()
            'The nonnegative integers repeated.'
        """
        return "The nonnegative integers repeated."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A004526._eval(n) for n in range(10)]
            [0, 0, 1, 1, 2, 2, 3, 3, 4, 4]
        """
        return ZZ(n // 2)


class A000326(SloaneSequence):
    def __init__(self):
        r"""
        Pentagonal numbers: `n(3n-1)/2`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000326;a
            Pentagonal numbers: n(3n-1)/2.
            sage: a(0)
            0
            sage: a(1)
            1
            sage: a(2)
            5
            sage: a(10)
            145
            sage: a.list(12)
            [0, 1, 5, 12, 22, 35, 51, 70, 92, 117, 145, 176]
            sage: a(1/3)
            Traceback (most recent call last):
            ...
            TypeError: input must be an int or Integer

        AUTHORS:

        - Jaap Spies (2007-01-26)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000326._repr_()
            'Pentagonal numbers: n(3n-1)/2.'
        """
        return "Pentagonal numbers: n(3n-1)/2."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000326._eval(n) for n in range(10)]
            [0, 1, 5, 12, 22, 35, 51, 70, 92, 117]
        """
        return ZZ(n * (3*n-1) // 2)


class A002378(SloaneSequence):
    def __init__(self):
        r"""
        Oblong (or pronic, or heteromecic) numbers: `n(n+1)`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A002378;a
            Oblong (or pronic, or heteromecic) numbers: n(n+1).
            sage: a(-1)
            Traceback (most recent call last):
            ...
            ValueError: input n (=-1) must be an integer >= 0
            sage: a(0)
            0
            sage: a(1)
            2
            sage: a(11)
            132
            sage: a.list(12)
            [0, 2, 6, 12, 20, 30, 42, 56, 72, 90, 110, 132]

        AUTHORS:

        - Jaap Spies (2007-01-26)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A002378._repr_()
            'Oblong (or pronic, or heteromecic) numbers: n(n+1).'
        """
        return "Oblong (or pronic, or heteromecic) numbers: n(n+1)."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A002378._eval(n) for n in range(10)]
            [0, 2, 6, 12, 20, 30, 42, 56, 72, 90]
        """
        return ZZ(n * (n+1))


class A002620(SloaneSequence):
    def __init__(self):
        r"""
        Quarter-squares: floor(n/2)\*ceiling(n/2). Equivalently,
        `\lfloor n^2/4 \rfloor`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A002620;a
            Quarter-squares: floor(n/2)*ceiling(n/2). Equivalently, floor(n^2/4).
            sage: a(0)
            0
            sage: a(1)
            0
            sage: a(2)
            1
            sage: a(10)
            25
            sage: a.list(12)
            [0, 0, 1, 2, 4, 6, 9, 12, 16, 20, 25, 30]

        AUTHORS:

        - Jaap Spies (2007-01-26)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A002620._repr_()
            'Quarter-squares: floor(n/2)*ceiling(n/2). Equivalently, floor(n^2/4).'
        """
        return "Quarter-squares: floor(n/2)*ceiling(n/2). Equivalently, floor(n^2/4)."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A002620._eval(n) for n in range(10)]
            [0, 0, 1, 2, 4, 6, 9, 12, 16, 20]
        """
        return ZZ(n**2 // 4)


class A005408(SloaneSequence):
    def __init__(self):
        r"""
        The odd numbers a(n) = 2n + 1.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A005408;a
            The odd numbers a(n) = 2n + 1.
            sage: a(-1)
            Traceback (most recent call last):
            ...
            ValueError: input n (=-1) must be an integer >= 0
            sage: a(0)
            1
            sage: a(4)
            9
            sage: a(11)
            23
            sage: a.list(12)
            [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23]

        AUTHORS:

        - Jaap Spies (2007-01-26)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A005408._repr_()
            'The odd numbers a(n) = 2n + 1.'
        """
        return "The odd numbers a(n) = 2n + 1."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A005408._eval(n) for n in range(10)]
            [1, 3, 5, 7, 9, 11, 13, 15, 17, 19]
        """
        return ZZ(2*n + 1)


class A000012(SloaneSequence):
    def __init__(self):
        r"""
        The all 1s sequence.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000012; a
            The all 1s sequence.
            sage: a(1)
            1
            sage: a(2007)
            1
            sage: a.list(12)
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

        AUTHORS:

        - Jaap Spies (2007-01-12)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000012._repr_()
            'The all 1s sequence.'
        """
        return 'The all 1s sequence.'

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000012._eval(n) for n in range(10)]
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        """
        return ZZ.one()


class A000120(SloaneSequence):
    def __init__(self):
        r"""
        1s-counting sequence: number of 1s in binary expansion of `n`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000120;a
            1s-counting sequence: number of 1s in binary expansion of n.
            sage: a(0)
            0
            sage: a(2)
            1
            sage: a(12)
            2
            sage: a.list(12)
            [0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3]

        AUTHORS:

        - Jaap Spies (2007-01-26)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000120._repr_()
            '1s-counting sequence: number of 1s in binary expansion of n.'
        """
        return '1s-counting sequence: number of 1s in binary expansion of n.'

    def f(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000120.f(n) for n in range(10)]
            [0, 1, 1, 2, 1, 2, 2, 3, 1, 2]
        """
        if n <= 1:
            return ZZ(n)
        return self.f(n // 2) + n % 2

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000120._eval(n) for n in range(10)]
            [0, 1, 1, 2, 1, 2, 2, 3, 1, 2]
        """
        return self.f(n)


class A010060(SloaneSequence):
    def __init__(self):
        r"""
        Thue-Morse sequence.

        Let `A_k` denote the first `2^k` terms; then
        `A_0 = 0`, and for `k \ge 0`,
        `A_{k+1} = A_k B_k`, where `B_k` is obtained
        from `A_k` by interchanging 0s and 1s.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A010060;a
            Thue-Morse sequence.
            sage: a(0)
            0
            sage: a(1)
            1
            sage: a(2)
            1
            sage: a(12)
            0
            sage: a.list(13)
            [0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0]

        AUTHORS:

        - Jaap Spies (2007-02-02)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A010060._repr_()
            'Thue-Morse sequence.'
        """
        return "Thue-Morse sequence."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A010060._eval(n) for n in range(10)]
            [0, 1, 1, 0, 1, 0, 0, 1, 1, 0]
        """
        return sloane.A000120(n) % 2


class A000069(SloaneSequence):
    def __init__(self):
        r"""
        Odious numbers: odd number of 1s in binary expansion.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000069; a
            Odious numbers: odd number of 1s in binary expansion.
            sage: a(0)
            1
            sage: a(2)
            4
            sage: a.list(9)
            [1, 2, 4, 7, 8, 11, 13, 14, 16]

        AUTHORS:

        - Jaap Spies (2007-02-02)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000069._repr_()
            'Odious numbers: odd number of 1s in binary expansion.'
        """
        return 'Odious numbers: odd number of 1s in binary expansion.'

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000069._eval(n) for n in range(10)]
            [1, 2, 4, 7, 8, 11, 13, 14, 16, 19]
        """
        return ZZ(2*n + 1) - sloane.A010060(n)


class A001969(SloaneSequence):
    def __init__(self):
        r"""
        Evil numbers: even number of 1s in binary expansion.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A001969;a
            Evil numbers: even number of 1s in binary expansion.
            sage: a(0)
            0
            sage: a(1)
            3
            sage: a(2)
            5
            sage: a(12)
            24
            sage: a.list(13)
            [0, 3, 5, 6, 9, 10, 12, 15, 17, 18, 20, 23, 24]

        AUTHORS:

        - Jaap Spies (2007-02-02)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A001969._repr_()
            'Evil numbers: even number of 1s in binary expansion.'
        """
        return 'Evil numbers: even number of 1s in binary expansion.'

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A001969._eval(n) for n in range(10)]
            [0, 3, 5, 6, 9, 10, 12, 15, 17, 18]
        """
        return ZZ(2*n) + sloane.A010060(n)


class A000290(SloaneSequence):
    def __init__(self):
        r"""
        The squares: `a(n) = n^2`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000290;a
            The squares: a(n) = n^2.
            sage: a(0)
            0
            sage: a(-1)
            Traceback (most recent call last):
            ...
            ValueError: input n (=-1) must be an integer >= 0
            sage: a(16)
            256
            sage: a.list(17)
            [0, 1, 4, 9, 16, 25, 36, 49, 64, 81, 100, 121, 144, 169, 196, 225, 256]

        AUTHORS:

        - Jaap Spies (2007-01-25)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000290._repr_()
            'The squares: a(n) = n^2.'
        """
        return "The squares: a(n) = n^2."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000290._eval(n) for n in range(10)]
            [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]
        """
        return ZZ(n ** 2)


class A000225(SloaneSequence):
    def __init__(self):
        r"""
        `2^n - 1`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000225;a
            2^n - 1.
            sage: a(0)
            0
            sage: a(-1)
            Traceback (most recent call last):
            ...
            ValueError: input n (=-1) must be an integer >= 0
            sage: a(12)
            4095
            sage: a.list(12)
            [0, 1, 3, 7, 15, 31, 63, 127, 255, 511, 1023, 2047]

        AUTHORS:

        - Jaap Spies (2007-01-25)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000225._repr_()
            '2^n - 1.'
        """
        return "2^n - 1."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000225._eval(n) for n in range(10)]
            [0, 1, 3, 7, 15, 31, 63, 127, 255, 511]
        """
        return ZZ(2**n - 1)


class A000015(SloaneSequence):
    def __init__(self):
        r"""
        Smallest prime power `\geq n` (where `1` is considered a prime
        power).

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000015; a
            Smallest prime power >= n.
            sage: a(1)
            1
            sage: a(8)
            8
            sage: a(305)
            307
            sage: a(-4)
            Traceback (most recent call last):
            ...
            ValueError: input n (=-4) must be a positive integer
            sage: a.list(12)
            [1, 2, 3, 4, 5, 7, 7, 8, 9, 11, 11, 13]
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer

        AUTHORS:

        - Jaap Spies (2007-01-18)
        """
        SloaneSequence.__init__(self, offset=1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000015._repr_()
            'Smallest prime power >= n.'
        """
        return "Smallest prime power >= n."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000015._eval(n) for n in range(1,11)]
            [1, 2, 3, 4, 5, 7, 7, 8, 9, 11]
        """
        if n == 1 or arith.is_prime_power(n):
            return n
        else:
            return arith.next_prime_power(n)


class A000016(SloaneSequence):
    def __init__(self):
        r"""
        Sloane's A000016.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000016; a
            Sloane's A000016.
            sage: a(1)
            1
            sage: a(0)
            1
            sage: a(8)
            16
            sage: a(75)
            251859545753048193000
            sage: a(-4)
            Traceback (most recent call last):
            ...
            ValueError: input n (=-4) must be an integer >= 0
            sage: a.list(12)
            [1, 1, 1, 2, 2, 4, 6, 10, 16, 30, 52, 94]

        AUTHORS:

        - Jaap Spies (2007-01-18)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000016._repr_()
            "Sloane's A000016."
        """
        return "Sloane's A000016."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000016._eval(n) for n in range(10)]
            [1, 1, 1, 2, 2, 4, 6, 10, 16, 30]
        """
        if n == 0:
            return ZZ.one()
        return ZZ(sum((i % 2) * arith.euler_phi(i) * 2**(n//i) / (2*n)
                      for i in arith.divisors(n)))


class A000032(SloaneSequence):
    def __init__(self):
        r"""
        Lucas numbers (beginning at 2): `L(n) = L(n-1) + L(n-2)`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000032; a
            Lucas numbers (beginning at 2): L(n) = L(n-1) + L(n-2).
            sage: a(0)
            2
            sage: a(1)
            1
            sage: a(8)
            47
            sage: a(200)
            627376215338105766356982006981782561278127
            sage: a(-4)
            Traceback (most recent call last):
            ...
            ValueError: input n (=-4) must be an integer >= 0
            sage: a.list(12)
            [2, 1, 3, 4, 7, 11, 18, 29, 47, 76, 123, 199]

        AUTHORS:

        - Jaap Spies (2007-01-18)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000032._repr_()
            'Lucas numbers (beginning at 2): L(n) = L(n-1) + L(n-2).'
        """
        return "Lucas numbers (beginning at 2): L(n) = L(n-1) + L(n-2)."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000032._eval(n) for n in range(10)]
            [2, 1, 3, 4, 7, 11, 18, 29, 47, 76]
        """
        if n == 0:
            return ZZ(2)
        elif n == 1:
            return ZZ.one()
        else:
            return sloane.A000045(n+1) + sloane.A000045(n-1)


# Theme numbers as strings of digits


class A004086(SloaneSequence):
    def __init__(self):
        r"""
        Read n backwards (referred to as `R(n)` in many
        sequences).

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A004086;a
            Read n backwards (referred to as R(n) in many sequences).
            sage: a(0)
            0
            sage: a(1)
            1
            sage: a(2)
            2
            sage: a(3333)
            3333
            sage: a(12345)
            54321
            sage: a.list(13)
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 11, 21]

        AUTHORS:

        - Jaap Spies (2007-02-02)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A004086._repr_()
            'Read n backwards (referred to as R(n) in many sequences).'
        """
        return "Read n backwards (referred to as R(n) in many sequences)."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A004086._eval(n) for n in range(10)]
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        """
        a = ''.join(reversed(str(n)))
        return ZZ(a)


class A002113(SloaneSequence):
    def __init__(self):
        r"""
        Palindromes in base 10.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A002113;a
            Palindromes in base 10.
            sage: a(0)
            0
            sage: a(1)
            1
            sage: a(2)
            2
            sage: a(12)
            33
            sage: a.list(13)
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 22, 33]

        AUTHORS:

        - Jaap Spies (2007-02-02)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A002113._repr_()
            'Palindromes in base 10.'
        """
        return "Palindromes in base 10."

    def _precompute(self, how_many=150):
        """
        EXAMPLES::

            sage: initial = len(sloane.A002113._b)
            sage: sloane.A002113._precompute()
            sage: len(sloane.A002113._b) - initial > 0
            True
        """
        try:
            self._b
            self._n
        except AttributeError:
            self._b = []
            self._n = self.offset
        self._b += [i for i in range(self._n, self._n + how_many)
                    if sloane.A004086(i) == i]
        self._n += how_many

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A002113._eval(n) for n in range(10)]
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        """
        try:
            return self._b[n]
        except (AttributeError, IndexError):
            self._precompute()
            # try again
            return self._eval(n)

    def list(self, n):
        """
        EXAMPLES::

            sage: sloane.A002113.list(15)
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 22, 33, 44, 55]
        """
        try:
            if len(self._b) <= n:
                raise IndexError
            else:
                return self._b[:n]
        except (AttributeError, IndexError):
            self._precompute()
            # try again
            return self.list(n)


class A000030(SloaneSequence):
    def __init__(self):
        r"""
        Initial digit of `n`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000030; a
            Initial digit of n
            sage: a(0)
            0
            sage: a(1)
            1
            sage: a(8)
            8
            sage: a(454)
            4
            sage: a(-4)
            Traceback (most recent call last):
            ...
            ValueError: input n (=-4) must be an integer >= 0
            sage: a.list(12)
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 1]

        AUTHORS:

        - Jaap Spies (2007-01-18)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000030._repr_()
            'Initial digit of n'
        """
        return "Initial digit of n"

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000030._eval(n) for n in range(10)]
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        """
        if n < 10:
            return n
        else:
            return self(n//10)


# Theme: primes and factoring
class A000040(SloaneSequence):
    def __init__(self):
        r"""
        The prime numbers.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000040; a
            The prime numbers.
            sage: a(1)
            2
            sage: a(8)
            19
            sage: a(305)
            2011
            sage: a.list(12)
            [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer

        AUTHORS:

        - Jaap Spies (2007-01-17)
        """
        SloaneSequence.__init__(self, offset=1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000040._repr_()
            'The prime numbers.'
        """
        return "The prime numbers."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000040._eval(n) for n in range(1,11)]
            [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
        """
        return arith.nth_prime(n)


class A002808(SloaneSequence):
    def __init__(self):
        r"""
        The composite numbers: numbers `n` of the form `xy`
        for `x > 1` and `y > 1`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A002808;a
            The composite numbers: numbers n of the form x*y for x > 1 and y > 1.
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(2)
            6
            sage: a(11)
            20
            sage: a.list(12)
            [4, 6, 8, 9, 10, 12, 14, 15, 16, 18, 20, 21]

        AUTHORS:

        - Jaap Spies (2007-01-26)
        """
        SloaneSequence.__init__(self, offset=1)
        self._b = [4]
        self._n = 5

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A002808._repr_()
            'The composite numbers: numbers n of the form x*y for x > 1 and y > 1.'
        """
        return "The composite numbers: numbers n of the form x*y for x > 1 and y > 1."

    def _precompute(self, how_many=150):
        """
        EXAMPLES::

            sage: initial = len(sloane.A002808._b)
            sage: sloane.A002808._precompute()
            sage: len(sloane.A002808._b) - initial > 0
            True
        """
        self._b += [i for i in range(self._n, self._n+how_many) if not arith.is_prime(i)]
        self._n += how_many

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A002808._eval(n) for n in range(1,11)]
            [4, 6, 8, 9, 10, 12, 14, 15, 16, 18]
        """
        try:
            return self._b[n-1]
        except (AttributeError, IndexError):
            self._precompute()
            # try again
            return self._eval(n)

    def list(self, n):
        """
        EXAMPLES::

            sage: sloane.A002808.list(10)
            [4, 6, 8, 9, 10, 12, 14, 15, 16, 18]
        """
        try:
            if len(self._b) <= n:
                raise IndexError
            else:
                return self._b[:n]
        except (AttributeError, IndexError):
            self._precompute()
            # try again
            return self.list(n)


class A018252(SloaneSequence):
    def __init__(self):
        r"""
        The nonprime numbers, starting with 1.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A018252;a
            The nonprime numbers.
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(1)
            1
            sage: a(2)
            4
            sage: a(9)
            15
            sage: a.list(10)
            [1, 4, 6, 8, 9, 10, 12, 14, 15, 16]

        AUTHORS:

        - Jaap Spies (2007-02-04)
        """
        SloaneSequence.__init__(self, offset=1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A018252._repr_()
            'The nonprime numbers.'
        """
        return "The nonprime numbers."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A018252._eval(n) for n in range(1,11)]
            [1, 4, 6, 8, 9, 10, 12, 14, 15, 16]
        """
        if n == 1:
            return ZZ.one()
        return sloane.A002808(n - 1)


class A000043(SloaneSequence):
    def __init__(self):
        r"""
        Primes `p` such that `2^p - 1` is prime.
        `2^p - 1` is then called a Mersenne prime.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000043;a
            Primes p such that 2^p - 1 is prime. 2^p - 1 is then called a Mersenne prime.
            sage: a(1)
            2
            sage: a(2)
            3
            sage: a(39)
            13466917
            sage: a(40)
            Traceback (most recent call last):
            ...
            IndexError: list index out of range
            sage: a.list(12)
            [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127]

        AUTHORS:

        - Jaap Spies (2007-01-26)
        """
        SloaneSequence.__init__(self, offset=1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000043._repr_()
            'Primes p such that 2^p - 1 is prime. 2^p - 1 is then called a Mersenne prime.'
        """
        return "Primes p such that 2^p - 1 is prime. 2^p - 1 is then called a Mersenne prime."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000043._eval(n) for n in range(1,11)]
            [2, 3, 5, 7, 13, 17, 19, 31, 61, 89]
        """
        try:
            return ZZ(self._b[n-1])
        except (AttributeError, IndexError):
            self._b = [2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,
                       4423,9689,9941,11213,19937,21701,23209,44497,86243,110503,132049,
                       216091,756839,859433,1257787,1398269,2976221,3021377,6972593,13466917]
            return ZZ(self._b[n-1])


class A000668(SloaneSequence):
    def __init__(self):
        r"""
        Mersenne primes (of form `2^p - 1` where `p` is a
        prime).

        (See A000043 for the values of `p`.)

        Warning: a(39) has 4,053,946 digits!

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000668;a
            Mersenne primes (of form 2^p - 1 where p is a prime). (See A000043 for the values of p.)
            sage: a(1)
            3
            sage: a(2)
            7
            sage: a(12)
            170141183460469231731687303715884105727

        Warning: a(39) has 4,053,946 digits!

        ::

            sage: a(40)
            Traceback (most recent call last):
            ...
            IndexError: list index out of range
            sage: a.list(8)
            [3, 7, 31, 127, 8191, 131071, 524287, 2147483647]

        AUTHORS:

        - Jaap Spies (2007-01-25)
        """
        SloaneSequence.__init__(self, offset=1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000668._repr_()
            'Mersenne primes (of form 2^p - 1 where p is a prime). (See A000043 for the values of p.)'
        """
        return "Mersenne primes (of form 2^p - 1 where p is a prime). (See A000043 for the values of p.)"

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000668._eval(n) for n in range(1,11)]
            [3,
             7,
             31,
             127,
             8191,
             131071,
             524287,
             2147483647,
             2305843009213693951,
             618970019642690137449562111]
        """
        return ZZ(2**sloane.A000043(n) - 1)


class A000396(SloaneSequence):
    def __init__(self):
        r"""
        Perfect numbers: equal to sum of proper divisors.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000396;a
            Perfect numbers: equal to sum of proper divisors.
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(1)
            6
            sage: a(2)
            28
            sage: a(7)
            137438691328
            sage: a.list(7)
            [6, 28, 496, 8128, 33550336, 8589869056, 137438691328]

        AUTHORS:

        - Jaap Spies (2007-01-25)
        """
        SloaneSequence.__init__(self, offset=1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000396._repr_()
            'Perfect numbers: equal to sum of proper divisors.'
        """
        return "Perfect numbers: equal to sum of proper divisors."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000396._eval(n) for n in range(1,6)]
            [6, 28, 496, 8128, 33550336]
        """
        p = sloane.A000043(n)
        return ZZ(2**(p-1) * (2**p - 1))


class A005100(SloaneSequence):
    def __init__(self):
        r"""
        Deficient numbers: `\sigma(n) < 2n`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A005100;a
            Deficient numbers: sigma(n) < 2n
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(1)
            1
            sage: a(2)
            2
            sage: a(12)
            14
            sage: a.list(12)
            [1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 13, 14]

        AUTHORS:

        - Jaap Spies (2007-01-26)
        """
        SloaneSequence.__init__(self, offset=1)
        self._b = [1]
        self._n = 2

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A005100._repr_()
            'Deficient numbers: sigma(n) < 2n'
        """
        return "Deficient numbers: sigma(n) < 2n"

    def _precompute(self, how_many=150):
        """
        EXAMPLES::

            sage: initial = len(sloane.A005100._b)
            sage: sloane.A005100._precompute()
            sage: len(sloane.A005100._b) - initial > 0
            True
        """
        self._b += [i for i in range(self._n, self._n+how_many) if arith.sigma(i) < 2*i]
        self._n += how_many

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A005100._eval(n) for n in range(1,10)]
            [1, 2, 3, 4, 5, 7, 8, 9, 10]
        """
        try:
            return self._b[n-1]
        except (AttributeError, IndexError):
            self._precompute()
            # try again
            return self._eval(n)

    def list(self, n):
        """
        EXAMPLES::

            sage: sloane.A005100.list(10)
            [1, 2, 3, 4, 5, 7, 8, 9, 10, 11]
        """
        try:
            if len(self._b) <= n:
                raise IndexError
            else:
                return self._b[:n]
        except (AttributeError, IndexError):
            self._precompute()
            # try again
            return self.list(n)


class A005101(SloaneSequence):
    def __init__(self):
        r"""
        Abundant numbers (sum of divisors of `n` exceeds
        `2n`).

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A005101;a
            Abundant numbers (sum of divisors of n exceeds 2n).
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(1)
            12
            sage: a(2)
            18
            sage: a(12)
            60
            sage: a.list(12)
            [12, 18, 20, 24, 30, 36, 40, 42, 48, 54, 56, 60]

        AUTHORS:

        - Jaap Spies (2007-01-26)
        """
        SloaneSequence.__init__(self, offset=1)
        self._b = [12]
        self._n = 18

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A005101._repr_()
            'Abundant numbers (sum of divisors of n exceeds 2n).'
        """
        return "Abundant numbers (sum of divisors of n exceeds 2n)."

    def _precompute(self, how_many=150):
        """
        EXAMPLES::

            sage: initial = len(sloane.A005101._b)
            sage: sloane.A005101._precompute()
            sage: len(sloane.A005101._b) - initial > 0
            True
        """
        self._b += [i for i in range(self._n, self._n+how_many) if arith.sigma(i) > 2*i]
        self._n += how_many

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A005101._eval(n) for n in range(1,11)]
            [12, 18, 20, 24, 30, 36, 40, 42, 48, 54]
        """
        try:
            return self._b[n-1]
        except (AttributeError, IndexError):
            self._precompute()
            # try again
            return self._eval(n)

    def list(self, n):
        """
        EXAMPLES::

            sage: sloane.A005101.list(10)
            [12, 18, 20, 24, 30, 36, 40, 42, 48, 54]
        """
        try:
            if len(self._b) <= n:
                raise IndexError
            else:
                return self._b[:n]
        except (AttributeError, IndexError):
            self._precompute()
            # try again
            return self.list(n)


class A002110(SloaneSequence):
    def __init__(self):
        r"""
        Primorial numbers (first definition): product of first `n`
        primes. Sometimes written `p\#`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A002110;a
            Primorial numbers (first definition): product of first n primes. Sometimes written p#.
            sage: a(0)
            1
            sage: a(2)
            6
            sage: a(8)
            9699690
            sage: a(17)
            1922760350154212639070
            sage: a.list(9)
            [1, 2, 6, 30, 210, 2310, 30030, 510510, 9699690]

        AUTHORS:

        - Jaap Spies (2007-01-25)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A002110._repr_()
            'Primorial numbers (first definition): product of first n primes. Sometimes written p#.'
        """
        return "Primorial numbers (first definition): product of first n primes. Sometimes written p#."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A002110._eval(n) for n in range(10)]
            [1, 2, 6, 30, 210, 2310, 30030, 510510, 9699690, 223092870]
        """
        return prod([sloane.A000040(i) for i in range(1, n + 1)])  # n-th prime = A000040(n)


class A000720(SloaneSequence):
    def __init__(self):
        r"""
        `pi(n)`, the number of primes `\le n`. Sometimes
        called `PrimePi(n)`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000720;a
            pi(n), the number of primes <= n. Sometimes called PrimePi(n)
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(2)
            1
            sage: a(8)
            4
            sage: a(1000)
            168
            sage: a.list(12)
            [0, 1, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5]

        AUTHORS:

        - Jaap Spies (2007-01-25)
        """
        SloaneSequence.__init__(self, offset=1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000720._repr_()
            'pi(n), the number of primes <= n. Sometimes called PrimePi(n)'
        """
        return "pi(n), the number of primes <= n. Sometimes called PrimePi(n)"

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000720._eval(n) for n in range(1,11)]
            [0, 1, 2, 2, 3, 3, 4, 4, 4, 4]
        """
        return prime_pi(n)


class A064553(SloaneSequence):
    def __init__(self):
        r"""
        `a(1) = 1`, `a(prime(i)) = i + 1` for
        `i > 0` and `a(u \cdot v) = a(u) \cdot a(v)` for
        `u, v > 0`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A064553;a
            a(1) = 1, a(prime(i)) = i+1 for i > 0 and a(u*v) = a(u)*a(v) for u,v > 0
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(1)
            1
            sage: a(2)
            2
            sage: a(9)
            9
            sage: a.list(16)
            [1, 2, 3, 4, 4, 6, 5, 8, 9, 8, 6, 12, 7, 10, 12, 16]

        AUTHORS:

        - Jaap Spies (2007-02-04)
        """
        SloaneSequence.__init__(self, offset=1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A064553._repr_()
            'a(1) = 1, a(prime(i)) = i+1 for i > 0 and a(u*v) = a(u)*a(v) for u,v > 0'
        """
        return "a(1) = 1, a(prime(i)) = i+1 for i > 0 and a(u*v) = a(u)*a(v) for u,v > 0"

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A064553._eval(n) for n in range(1,11)]
            [1, 2, 3, 4, 4, 6, 5, 8, 9, 8]
        """
        return prod([(prime_pi(p)+1)**e for p, e in arith.factor(n)])


class A001055(SloaneSequence):
    def __init__(self):
        r"""
        Number of ways of factoring `n` with all factors 1.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A001055;a
            Number of ways of factoring n with all factors >1.
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(1)
            1
            sage: a(2)
            1
            sage: a(9)
            2
            sage: a.list(16)
            [1, 1, 1, 2, 1, 2, 1, 3, 2, 2, 1, 4, 1, 2, 2, 5]

        AUTHORS:

        - Jaap Spies (2007-02-04)
        """
        SloaneSequence.__init__(self, offset=1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A001055._repr_()
            'Number of ways of factoring n with all factors >1.'
        """
        return "Number of ways of factoring n with all factors >1."

    def nwf(self, n, m):
        """
        EXAMPLES::

            sage: sloane.A001055.nwf(4,1)
            0
            sage: sloane.A001055.nwf(4,2)
            1
            sage: sloane.A001055.nwf(4,3)
            1
            sage: sloane.A001055.nwf(4,4)
            2
        """
        if n == 1:
            return ZZ.one()
        if arith.is_prime(n):
            if m < n:
                return ZZ.zero()
            else:
                return ZZ.one()
        s = ZZ(0)
        for d in arith.divisors(n):
            if d > 1 and d <= m and d < n:
                s += self.nwf(n // d, d)
        if n <= m:
            s += 1
        return s

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A001055._eval(n) for n in range(1,11)]
            [1, 1, 1, 2, 1, 2, 1, 3, 2, 2]
        """
        return self.nwf(n, n)


class A006530(SloaneSequence):
    def __init__(self):
        r"""
        Largest prime dividing `n` (with `a(1)=1`).

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A006530;a
            Largest prime dividing n (with a(1)=1).
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(1)
            1
            sage: a(2)
            2
            sage: a(8)
            2
            sage: a(11)
            11
            sage: a.list(15)
            [1, 2, 3, 2, 5, 3, 7, 2, 3, 5, 11, 3, 13, 7, 5]

        AUTHORS:

        - Jaap Spies (2007-01-25)
        """
        SloaneSequence.__init__(self, offset=1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A006530._repr_()
            'Largest prime dividing n (with a(1)=1).'
        """
        return "Largest prime dividing n (with a(1)=1)."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A006530._eval(n) for n in range(1,11)]
            [1, 2, 3, 2, 5, 3, 7, 2, 3, 5]
        """
        if n == 1:
            return ZZ.one()
        return max(p for p,_ in arith.factor(n))


class A000961(SloaneSequence):
    def __init__(self):
        r"""
        Prime powers.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000961;a
            Prime powers.
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(2)
            2
            sage: a(12)
            17
            sage: a.list(12)
            [1, 2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17]

        AUTHORS:

        - Jaap Spies (2007-01-25)
        """
        SloaneSequence.__init__(self, offset=1)
        self._b = [1]
        self._n = 2

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000961._repr_()
            'Prime powers.'
        """
        return "Prime powers."

    def _precompute(self, how_many=150):
        """
        EXAMPLES::

            sage: initial = len(sloane.A000961._b)
            sage: sloane.A000961._precompute()
            sage: len(sloane.A000961._b) - initial > 0
            True
        """
        self._b += [i for i in range(self._n, self._n+how_many) if len([p for p,_ in arith.factor(i)]) == 1]
        self._n += how_many

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000961._eval(n) for n in range(1,11)]
            [1, 2, 3, 4, 5, 7, 8, 9, 11, 13]
        """
        try:
            return self._b[n-1]
        except (AttributeError, IndexError):
            self._precompute()
            # try again
            return self._eval(n)

    def list(self, n):
        """
        EXAMPLES::

            sage: sloane.A000961.list(10)
            [1, 2, 3, 4, 5, 7, 8, 9, 11, 13]
        """
        try:
            if len(self._b) <= n:
                raise IndexError
            else:
                return self._b[:n]
        except (AttributeError, IndexError):
            self._precompute()
            # try again
            return self.list(n)


class A005117(SloaneSequence):
    def __init__(self):
        r"""
        Square-free numbers.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A005117;a
            Square-free numbers.
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(2)
            2
            sage: a(12)
            17
            sage: a.list(12)
            [1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17]

        AUTHORS:

        - Jaap Spies (2007-01-25)
        """
        SloaneSequence.__init__(self, offset=1)
        self._b = [1]
        self._n = 2

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A005117._repr_()
            'Square-free numbers.'
        """
        return "Square-free numbers."

    def _precompute(self, how_many=150):
        """
        EXAMPLES::

            sage: initial = len(sloane.A005117._b)
            sage: sloane.A005117._precompute()
            sage: len(sloane.A005117._b) - initial > 0
            True
        """
        self._b += [i for i in range(self._n, self._n+how_many)
                    if max(e for _, e in arith.factor(i)) <= 1]
        self._n += how_many

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A005117._eval(n) for n in range(1,11)]
            [1, 2, 3, 5, 6, 7, 10, 11, 13, 14]
        """
        try:
            return self._b[n-1]
        except (AttributeError, IndexError):
            self._precompute()
            # try again
            return self._eval(n)

    def list(self, n):
        """
        EXAMPLES::

            sage: sloane.A005117.list(10)
            [1, 2, 3, 5, 6, 7, 10, 11, 13, 14]
        """
        try:
            if len(self._b) <= n:
                raise IndexError
            else:
                return self._b[:n]
        except (AttributeError, IndexError):
            self._precompute()
            # try again
            return self.list(n)


class A020639(SloaneSequence):
    def __init__(self):
        r"""
        Least prime dividing `n` with `a(1)=1`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A020639;a
            Least prime dividing n (a(1)=1).
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(1)
            1
            sage: a(13)
            13
            sage: a.list(14)
            [1, 2, 3, 2, 5, 2, 7, 2, 3, 2, 11, 2, 13, 2]

        AUTHORS:

        - Jaap Spies (2007-01-25)
        """
        SloaneSequence.__init__(self, offset=1)
        self._b = [1]
        self._n = 2

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A020639._repr_()
            'Least prime dividing n (a(1)=1).'
        """
        return "Least prime dividing n (a(1)=1)."

    def _precompute(self, how_many=50):
        """
        EXAMPLES::

            sage: initial = len(sloane.A020639._b)
            sage: sloane.A020639._precompute(10)
            sage: len(sloane.A020639._b) - initial == 10
            True
        """
        self._b += [min(p for p,_ in arith.factor(i)) for i in range(self._n, self._n+how_many)]
        self._n += how_many

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A020639._eval(n) for n in range(1,11)]
            [1, 2, 3, 2, 5, 2, 7, 2, 3, 2]
        """
        try:
            return self._b[n-self.offset]
        except (AttributeError, IndexError):
            self._precompute()
            # try again
            return self._eval(n)

    def list(self, n):
        """
        EXAMPLES::

            sage: sloane.A020639.list(10)
            [1, 2, 3, 2, 5, 2, 7, 2, 3, 2]
        """
        try:
            if len(self._b) <= n:
                raise IndexError
            else:
                return self._b[:n]
        except (AttributeError, IndexError):
            self._precompute()
            # try again
            return self.list(n)


class A000041(SloaneSequence):
    def __init__(self):
        r"""
        `a(n)` = number of partitions of `n` (the partition
        numbers).

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000041;a
            a(n) = number of partitions of n (the partition numbers).
            sage: a(0)
            1
            sage: a(2)
            2
            sage: a(8)
            22
            sage: a(200)
            3972999029388
            sage: a.list(9)
            [1, 1, 2, 3, 5, 7, 11, 15, 22]

        AUTHORS:

        - Jaap Spies (2007-01-18)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000041._repr_()
            'a(n) = number of partitions of n (the partition numbers).'
        """
        return "a(n) = number of partitions of n (the partition numbers)."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000041._eval(n) for n in range(1,11)]
            [1, 2, 3, 5, 7, 11, 15, 22, 30, 42]
        """
        from sage.combinat.partition import Partitions
        return Partitions(n).cardinality()


class A000045(SloaneSequence):
    def __init__(self):
        r"""
        Sequence of Fibonacci numbers, offset 0,4.

        REFERENCES:

        - S. Plouffe, Project Gutenberg, The First 1001 Fibonacci
          Numbers,
          http://ibiblio.org/pub/docs/books/gutenberg/etext01/fbncc10.txt

        We have one more. Our first Fibonacci number is 0.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000045; a
            Fibonacci numbers with index n >= 0
            sage: a(0)
            0
            sage: a(1)
            1
            sage: a.list(12)
            [0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89]
            sage: a(1/3)
            Traceback (most recent call last):
            ...
            TypeError: input must be an int or Integer

        AUTHORS:

        - Jaap Spies (2007-01-13)
        """
        SloaneSequence.__init__(self, offset=0)
        self._b = []
        self._precompute()  # force precomputation, e.g. a(0) will fail when asked first

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000045._repr_()
            'Fibonacci numbers with index n >= 0'
        """
        return "Fibonacci numbers with index n >= 0"

    def _precompute(self, how_many=500):
        """
        EXAMPLES::

            sage: initial = len(sloane.A000045._b)
            sage: sloane.A000045._precompute(10)
            sage: len(sloane.A000045._b) - initial > 0
            True
        """
        try:
            f = self._f
        except AttributeError:
            self._f = self.fib()
            f = self._f
        self._b += [next(f) for i in range(how_many)]

    def fib(self):
        """
        Return a generator over all Fibonacci numbers, starting with 0.

        EXAMPLES::

            sage: it = sloane.A000045.fib()
            sage: [next(it) for i in range(10)]
            [0, 1, 1, 2, 3, 5, 8, 13, 21, 34]
        """
        x, y = ZZ.zero(), ZZ.one()
        yield x
        while True:
            x, y = y, x+y
            yield x

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000045._eval(n) for n in range(1,11)]
            [1, 1, 2, 3, 5, 8, 13, 21, 34, 55]
        """
        if len(self._b) < n:
            self._precompute(n - len(self._b) + 1)
        return self._b[n]

    def list(self, n):
        """
        EXAMPLES::

            sage: sloane.A000045.list(10)
            [0, 1, 1, 2, 3, 5, 8, 13, 21, 34]
        """
        self._eval(n)   # force computation
        return self._b[:n]


class A000108(SloaneSequence):
    def __init__(self):
        r"""
        Catalan numbers:
        `C_n = \frac{\binom{2n}{n}}{n+1} = \frac{(2n)!}{n!(n+1)!}`.

        Also called Segner numbers.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000108;a
            Catalan numbers: C(n) = binomial(2n,n)/(n+1) = (2n)!/(n!(n+1)!). Also called Segner numbers.
            sage: a(0)
            1
            sage: a.offset
            0
            sage: a(8)
            1430
            sage: a(40)
            2622127042276492108820
            sage: a.list(9)
            [1, 1, 2, 5, 14, 42, 132, 429, 1430]

        AUTHORS:

        - Jaap Spies (2007-01-12)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000108._repr_()
            'Catalan numbers: C(n) = binomial(2n,n)/(n+1) = (2n)!/(n!(n+1)!). Also called Segner numbers.'
        """
        return "Catalan numbers: C(n) = binomial(2n,n)/(n+1) = (2n)!/(n!(n+1)!). Also called Segner numbers."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000108._eval(n) for n in range(10)]
            [1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862]
        """
        return combinat.catalan_number(n)


class A001006(SloaneSequence):
    def __init__(self):
        r"""
        Motzkin numbers: number of ways of drawing any number of
        nonintersecting chords among `n` points on a circle.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A001006;a
            Motzkin numbers: number of ways of drawing any number of nonintersecting chords among n points on a circle.
            sage: a(0)
            1
            sage: a(1)
            1
            sage: a(2)
            2
            sage: a(12)
            15511
            sage: a.list(13)
            [1, 1, 2, 4, 9, 21, 51, 127, 323, 835, 2188, 5798, 15511]

        AUTHORS:

        - Jaap Spies (2007-02-02)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A001006._repr_()
            'Motzkin numbers: number of ways of drawing any number of nonintersecting chords among n points on a circle.'
        """
        return "Motzkin numbers: number of ways of drawing any number of nonintersecting chords among n points on a circle."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A001006._eval(n) for n in range(10)]
            [1, 1, 2, 4, 9, 21, 51, 127, 323, 835]
        """
        return sum((-1)**(n-k)*arith.binomial(n, k)*sloane.A000108(k+1) for k in range(n+1))


class A000079(SloaneSequence):
    def __init__(self):
        r"""
        Powers of 2: `a(n) = 2^n`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000079;a
            Powers of 2: a(n) = 2^n.
            sage: a(0)
            1
            sage: a(2)
            4
            sage: a(8)
            256
            sage: a(100)
            1267650600228229401496703205376
            sage: a.list(9)
            [1, 2, 4, 8, 16, 32, 64, 128, 256]

        AUTHORS:

        - Jaap Spies (2007-01-18)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000079._repr_()
            'Powers of 2: a(n) = 2^n.'
        """
        return "Powers of 2: a(n) = 2^n."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000079._eval(n) for n in range(10)]
            [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
        """
        return ZZ(2 ** n)


class A000578(SloaneSequence):
    def __init__(self):
        r"""
        The cubes: `a(n) = n^3`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000578;a
            The cubes: n^3
            sage: a(-1)
            Traceback (most recent call last):
            ...
            ValueError: input n (=-1) must be an integer >= 0
            sage: a(0)
            0
            sage: a(3)
            27
            sage: a(11)
            1331
            sage: a.list(12)
            [0, 1, 8, 27, 64, 125, 216, 343, 512, 729, 1000, 1331]

        AUTHORS:

        - Jaap Spies (2007-01-26)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000578._repr_()
            'The cubes: n^3'
        """
        return "The cubes: n^3"

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000578._eval(n) for n in range(10)]
            [0, 1, 8, 27, 64, 125, 216, 343, 512, 729]
        """
        return ZZ(n ** 3)


class A000244(SloaneSequence):
    def __init__(self):
        r"""
        Powers of 3: `a(n) = 3^n`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000244;a
            Powers of 3: a(n) = 3^n.
            sage: a(-1)
            Traceback (most recent call last):
            ...
            ValueError: input n (=-1) must be an integer >= 0
            sage: a(0)
            1
            sage: a(3)
            27
            sage: a(11)
            177147
            sage: a.list(12)
            [1, 3, 9, 27, 81, 243, 729, 2187, 6561, 19683, 59049, 177147]

        AUTHORS:

        - Jaap Spies (2007-01-26)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000244._repr_()
            'Powers of 3: a(n) = 3^n.'
        """
        return "Powers of 3: a(n) = 3^n."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000244._eval(n) for n in range(10)]
            [1, 3, 9, 27, 81, 243, 729, 2187, 6561, 19683]
        """
        return ZZ(3**n)


class A000302(SloaneSequence):
    def __init__(self):
        r"""
        Powers of 4: `a(n) = 4^n`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000302;a
            Powers of 4: a(n) = 4^n.
            sage: a(0)
            1
            sage: a(1)
            4
            sage: a(2)
            16
            sage: a(10)
            1048576
            sage: a.list(12)
            [1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576, 4194304]

        AUTHORS:

        - Jaap Spies (2007-01-26)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000302._repr_()
            'Powers of 4: a(n) = 4^n.'
        """
        return "Powers of 4: a(n) = 4^n."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000302._eval(n) for n in range(10)]
            [1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144]
        """
        return ZZ(4 ** n)


class A000583(SloaneSequence):
    def __init__(self):
        r"""
        Fourth powers: `a(n) = n^4`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000583;a
            Fourth powers: n^4.
            sage: a(0.0)
            Traceback (most recent call last):
            ...
            TypeError: input must be an int or Integer
            sage: a(1)
            1
            sage: a(2)
            16
            sage: a(9)
            6561
            sage: a.list(10)
            [0, 1, 16, 81, 256, 625, 1296, 2401, 4096, 6561]

        AUTHORS:

        - Jaap Spies (2007-02-04)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000583._repr_()
            'Fourth powers: n^4.'
        """
        return "Fourth powers: n^4."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000583._eval(n) for n in range(10)]
            [0, 1, 16, 81, 256, 625, 1296, 2401, 4096, 6561]
        """
        return ZZ(n ** 4)


class A000142(SloaneSequence):
    def __init__(self):
        r"""
        Factorial numbers: `n! = 1 \cdot 2 \cdot 3 \cdots n`.

        Order of symmetric group `S_n`, number of permutations of
        `n` letters.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000142;a
            Factorial numbers: n! = 1*2*3*4*...*n (order of symmetric group S_n, number of permutations of n letters).
            sage: a(0)
            1
            sage: a(8)
            40320
            sage: a(40)
            815915283247897734345611269596115894272000000000
            sage: a.list(9)
            [1, 1, 2, 6, 24, 120, 720, 5040, 40320]

        AUTHORS:

        - Jaap Spies (2007-01-12)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000142._repr_()
            'Factorial numbers: n! = 1*2*3*4*...*n (order of symmetric group S_n, number of permutations of n letters).'
        """
        return "Factorial numbers: n! = 1*2*3*4*...*n (order of symmetric group S_n, number of permutations of n letters)."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000142._eval(n) for n in range(10)]
            [1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880]
        """
        return arith.factorial(n)


class A000085(SloaneSequence):
    def __init__(self):
        r"""
        Number of self-inverse permutations on `n` letters, also
        known as involutions; number of Young tableaux with `n`
        cells.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000085;a
            Number of self-inverse permutations on n letters.
            sage: a(0)
            1
            sage: a(1)
            1
            sage: a(2)
            2
            sage: a(12)
            140152
            sage: a.list(13)
            [1, 1, 2, 4, 10, 26, 76, 232, 764, 2620, 9496, 35696, 140152]

        AUTHORS:

        - Jaap Spies (2007-02-03)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000085._repr_()
            'Number of self-inverse permutations on n letters.'
        """
        return "Number of self-inverse permutations on n letters."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000085._eval(n) for n in range(10)]
            [1, 1, 2, 4, 10, 26, 76, 232, 764, 2620]
        """
        return sum(arith.factorial(n) // (arith.factorial(n-2*k) * (2**k) * arith.factorial(k))
                   for k in range(n//2+1))


class A001189(SloaneSequence):
    def __init__(self):
        r"""
        Number of degree-n permutations of order exactly 2.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A001189;a
            Number of degree-n permutations of order exactly 2.
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(1)
            0
            sage: a(2)
            1
            sage: a(12)
            140151
            sage: a.list(13)
            [0, 1, 3, 9, 25, 75, 231, 763, 2619, 9495, 35695, 140151, 568503]

        AUTHORS:

        - Jaap Spies (2007-02-03)
        """
        SloaneSequence.__init__(self, offset=1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A001189._repr_()
            'Number of degree-n permutations of order exactly 2.'
        """
        return "Number of degree-n permutations of order exactly 2."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A001189._eval(n) for n in range(1,11)]
            [0, 1, 3, 9, 25, 75, 231, 763, 2619, 9495]
        """
        return sloane.A000085(n) - 1


class A000670(SloaneSequence):
    def __init__(self):
        r"""
        Number of preferential arrangements of `n` labeled
        elements; or number of weak orders on `n` labeled
        elements.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000670;a
            Number of preferential arrangements of n labeled elements.
            sage: a(0)
            1
            sage: a(1)
            1
            sage: a(2)
            3
            sage: a(9)
            7087261
            sage: a.list(10)
            [1, 1, 3, 13, 75, 541, 4683, 47293, 545835, 7087261]

        AUTHORS:

        - Jaap Spies (2007-02-03)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000670._repr_()
            'Number of preferential arrangements of n labeled elements.'
        """
        return "Number of preferential arrangements of n labeled elements."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000670._eval(n) for n in range(1,10)]
            [1, 3, 13, 75, 541, 4683, 47293, 545835, 7087261]
        """
        # a(n) = Sum from k=1 to n of k! StirlingS2(n, k)
        if n == 0:
            return ZZ.one()
        return sum(arith.factorial(k) * combinat.stirling_number2(n, k)
                   for k in range(1, n+1))


class A006318(SloaneSequence):
    def __init__(self):
        r"""
        Large Schroeder numbers.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A006318;a
            Large Schroeder numbers.
            sage: a(0)
            1
            sage: a(1)
            2
            sage: a(2)
            6
            sage: a(9)
            206098
            sage: a.list(10)
            [1, 2, 6, 22, 90, 394, 1806, 8558, 41586, 206098]

        AUTHORS:

        - Jaap Spies (2007-02-03)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A006318._repr_()
            'Large Schroeder numbers.'
        """
        return "Large Schroeder numbers."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A006318._eval(n) for n in range(10)]
            [1, 2, 6, 22, 90, 394, 1806, 8558, 41586, 206098]
        """
        if n == 0:
            return ZZ.one()
#  (PARI) a(n)=if(n<1, 1, sum(k=0, n, 2^k*binomial(n, k)*binomial(n, k-1))/n)
        return ZZ(sum(2**k * arith.binomial(n, k) * arith.binomial(n, k-1) for k in range(n+1)) // n)


class A000165(SloaneSequence):
    def __init__(self):
        r"""
        Double factorial numbers: `(2n)!! = 2^n*n!`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000165;a
            Double factorial numbers: (2n)!! = 2^n*n!.
            sage: a(0)
            1
            sage: a.offset
            0
            sage: a(8)
            10321920
            sage: a(20)
            2551082656125828464640000
            sage: a.list(9)
            [1, 2, 8, 48, 384, 3840, 46080, 645120, 10321920]

        AUTHORS:

        - Jaap Spies (2007-01-24)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000165._repr_()
            'Double factorial numbers: (2n)!! = 2^n*n!.'
        """
        return "Double factorial numbers: (2n)!! = 2^n*n!."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000165._eval(n) for n in range(10)]
            [1, 2, 8, 48, 384, 3840, 46080, 645120, 10321920, 185794560]
        """
        return (2**n) * arith.factorial(n)


class A001147(SloaneSequence):
    def __init__(self):
        r"""
        Double factorial numbers:
        `(2n-1)!! = 1 \cdot 3 \cdot 5 \cdots (2n-1)`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A001147;a
            Double factorial numbers: (2n-1)!! = 1.3.5....(2n-1).
            sage: a(0)
            1
            sage: a.offset
            0
            sage: a(8)
            2027025
            sage: a(20)
            319830986772877770815625
            sage: a.list(9)
            [1, 1, 3, 15, 105, 945, 10395, 135135, 2027025]

        AUTHORS:

        - Jaap Spies (2007-01-24)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A001147._repr_()
            'Double factorial numbers: (2n-1)!! = 1.3.5....(2n-1).'
        """
        return "Double factorial numbers: (2n-1)!! = 1.3.5....(2n-1)."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A001147._eval(n) for n in range(10)]
            [1, 1, 3, 15, 105, 945, 10395, 135135, 2027025, 34459425]
        """
        return arith.factorial(2*n) / (arith.factorial(n)*2**n)


class A006882(SloaneSequence):
    def __init__(self):
        r"""
        Double factorials `n!!`: `a(n)=n \cdot a(n-2)`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A006882;a
            Double factorials n!!: a(n)=n*a(n-2).
            sage: a(0)
            1
            sage: a(2)
            2
            sage: a(8)
            384
            sage: a(20)
            3715891200
            sage: a.list(9)
            [1, 1, 2, 3, 8, 15, 48, 105, 384]

        AUTHORS:

        - Jaap Spies (2007-01-24)
        """
        SloaneSequence.__init__(self, offset=0)
        self._b = []
        self._precompute(2)  # force precomputation, e.g. a(0) will fail when asked first

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A006882._repr_()
            'Double factorials n!!: a(n)=n*a(n-2).'
        """
        return "Double factorials n!!: a(n)=n*a(n-2)."

    def _precompute(self, how_many=10):
        """
        EXAMPLES::

            sage: initial = len(sloane.A006882._b)
            sage: sloane.A006882._precompute(10)
            sage: len(sloane.A006882._b) - initial == 10
            True
        """
        try:
            f = self._f
        except AttributeError:
            self._f = self.df()
            f = self._f
        self._b += [next(f) for i in range(how_many)]

    def df(self):
        r"""
        Double factorials n!!: a(n)=n\*a(n-2).

        EXAMPLES::

            sage: it = sloane.A006882.df()
            sage: [next(it) for i in range(10)]
            [1, 1, 2, 3, 8, 15, 48, 105, 384, 945]
        """
        x = ZZ.one()
        k = 1
        y = x
        yield x
        while True:
            k = k+1
            x, y = y, k*x
            yield x

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A006882._eval(n) for n in range(10)]
            [1, 1, 2, 3, 8, 15, 48, 105, 384, 945]
        """
        if len(self._b) <= n:
            self._precompute(n - len(self._b) + 1)
        return self._b[n]

    def list(self, n):
        """
        EXAMPLES::

            sage: sloane.A006882.list(10)
            [1, 1, 2, 3, 8, 15, 48, 105, 384, 945]
        """
        self._eval(n)   # force computation
        return self._b[:n]


class A000984(SloaneSequence):
    def __init__(self):
        r"""
        Central binomial coefficients:
        `\binom{2n}{n} = \frac {(2n)!} {(n!)^2}`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000984;a
            Central binomial coefficients: C(2n,n) = (2n)!/(n!)^2
            sage: a(0)
            1
            sage: a(2)
            6
            sage: a(8)
            12870
            sage: a.list(9)
            [1, 2, 6, 20, 70, 252, 924, 3432, 12870]

        AUTHORS:

        - Jaap Spies (2007-01-26)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000984._repr_()
            'Central binomial coefficients: C(2n,n) = (2n)!/(n!)^2'
        """
        return "Central binomial coefficients: C(2n,n) = (2n)!/(n!)^2"

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000984._eval(n) for n in range(10)]
            [1, 2, 6, 20, 70, 252, 924, 3432, 12870, 48620]
        """
        return arith.binomial(2 * n, n)


class A001405(SloaneSequence):
    def __init__(self):
        r"""
        Central binomial coefficients:
        `\binom{n}{\lfloor \frac {n}{ 2} \rfloor}`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A001405;a
            Central binomial coefficients: C(n,floor(n/2)).
            sage: a(0)
            1
            sage: a(2)
            2
            sage: a(12)
            924
            sage: a.list(12)
            [1, 1, 2, 3, 6, 10, 20, 35, 70, 126, 252, 462]

        AUTHORS:

        - Jaap Spies (2007-01-26)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A001405._repr_()
            'Central binomial coefficients: C(n,floor(n/2)).'
        """
        return "Central binomial coefficients: C(n,floor(n/2))."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A001405._eval(n) for n in range(10)]
            [1, 1, 2, 3, 6, 10, 20, 35, 70, 126]
        """
        return arith.binomial(n, n//2)


class A000292(SloaneSequence):
    def __init__(self):
        r"""
        Tetrahedral (or pyramidal) numbers:
        `\binom{n+2}{3} = n(n+1)(n+2)/6`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000292;a
            Tetrahedral (or pyramidal) numbers: C(n+2,3) = n(n+1)(n+2)/6.
            sage: a(0)
            0
            sage: a(2)
            4
            sage: a(11)
            286
            sage: a.list(12)
            [0, 1, 4, 10, 20, 35, 56, 84, 120, 165, 220, 286]

        AUTHORS:

        - Jaap Spies (2007-01-26)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000292._repr_()
            'Tetrahedral (or pyramidal) numbers: C(n+2,3) = n(n+1)(n+2)/6.'
        """
        return "Tetrahedral (or pyramidal) numbers: C(n+2,3) = n(n+1)(n+2)/6."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000292._eval(n) for n in range(10)]
            [0, 1, 4, 10, 20, 35, 56, 84, 120, 165]
        """
        return ZZ(n * (n+1) * (n+2) // 6)  # or arith.binomial(n+2,3))


class A000330(SloaneSequence):
    def __init__(self):
        r"""
        Square pyramidal numbers"
        `0^2 + 1^2 \cdots n^2 = n(n+1)(2n+1)/6`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000330;a
            Square pyramidal numbers: 0^2+1^2+2^2+...+n^2 = n(n+1)(2n+1)/6.
            sage: a(-1)
            Traceback (most recent call last):
            ...
            ValueError: input n (=-1) must be an integer >= 0
            sage: a(0)
            0
            sage: a(3)
            14
            sage: a(11)
            506
            sage: a.list(12)
            [0, 1, 5, 14, 30, 55, 91, 140, 204, 285, 385, 506]

        AUTHORS:

        - Jaap Spies (2007-01-26)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000330._repr_()
            'Square pyramidal numbers: 0^2+1^2+2^2+...+n^2 = n(n+1)(2n+1)/6.'
        """
        return "Square pyramidal numbers: 0^2+1^2+2^2+...+n^2 = n(n+1)(2n+1)/6."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000330._eval(n) for n in range(10)]
            [0, 1, 5, 14, 30, 55, 91, 140, 204, 285]
        """
        return ZZ(n * (n+1) * (2*n+1) // 6)


# Theme:  maximal permanent of an m x n (0,1)- matrix:
# Seok-Zun Song et al.  Extremes of permanents of (0,1)-matrices, p. 201-202.

class ExtremesOfPermanentsSequence(SloaneSequence):
    def _precompute(self, how_many=20):
        """
        EXAMPLES::

            sage: sloane.A000153._precompute()
            sage: v1 = len(sloane.A000153._b)
            sage: sloane.A000153._precompute(10)
            sage: len(sloane.A000153._b) - v1
            10
        """
        try:
            f = self._f
        except AttributeError:
            self._f = self.gen(*self._a0a1d)
            f = self._f
        self._b += [next(f) for i in range(how_many)]

    def gen(self, a0, a1, d):
        """
        EXAMPLES::

            sage: it = sloane.A000153.gen(0,1,2)
            sage: [next(it) for i in range(5)]
            [0, 1, 2, 7, 32]
        """
        x, y = ZZ(a0), ZZ(a1)
        k = self._k
        yield x
        while True:
            k = k+1
            x, y = y, (k)*y+(k-d)*x
            yield x

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000153._eval(n) for n in range(8)]
            [0, 1, 2, 7, 32, 181, 1214, 9403]
        """
        if len(self._b) <= n:
            self._precompute(n - len(self._b) + 1)
        return self._b[n - self.offset]

    def list(self, n):
        """
        EXAMPLES::

            sage: sloane.A000153.list(8)
            [0, 1, 2, 7, 32, 181, 1214, 9403]
        """
        self._eval(n)   # force computation
        return self._b[:n]
    _k = 1


class A000153(ExtremesOfPermanentsSequence):
    def __init__(self):
        r"""
        `a(n) = n*a(n-1) + (n-2)*a(n-2)`, with `a(0) = 0`,
        `a(1) = 1`.

        With offset 1, permanent of (0,1)-matrix of size
        `n \times (n+d)` with `d=2` and `n` zeros
        not on a line. This is a special case of Theorem 2.3 of Seok-Zun
        Song et al. Extremes of permanents of (0,1)-matrices, p. 201-202.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000153; a
            a(n) = n*a(n-1) + (n-2)*a(n-2), with a(0) = 0, a(1) = 1.
            sage: a(0)
            0
            sage: a(1)
            1
            sage: a(8)
            82508
            sage: a(20)
            10315043624498196944
            sage: a.list(8)
            [0, 1, 2, 7, 32, 181, 1214, 9403]

        AUTHORS:

        - Jaap Spies (2007-01-13)
        """
        SloaneSequence.__init__(self, offset=0)
        self._b = []
        self._a0a1d = (0, 1, 2)
        self._precompute(2)  # force precomputation, e.g. a(0) will fail when asked first

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000153._repr_()
            'a(n) = n*a(n-1) + (n-2)*a(n-2), with a(0) = 0, a(1) = 1.'
        """
        return "a(n) = n*a(n-1) + (n-2)*a(n-2), with a(0) = 0, a(1) = 1."


class A000255(ExtremesOfPermanentsSequence):
    def __init__(self):
        r"""
        `a(n) = n*a(n-1) + (n-1)*a(n-2)`, with `a(0) = 1`,
        `a(1) = 1`.

        With offset 1, permanent of (0,1)-matrix of size
        `n \times (n+d)` with `d=1` and `n` zeros
        not on a line. This is a special case of Theorem 2.3 of Seok-Zun
        Song et al. Extremes of permanents of (0,1)-matrices, p. 201-202.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000255;a
            a(n) = n*a(n-1) + (n-1)*a(n-2), a(0) = 1, a(1) = 1.
            sage: a(0)
            1
            sage: a(1)
            1
            sage: a.offset
            0
            sage: a(8)
            148329
            sage: a(22)
            9923922230666898717143
            sage: a.list(9)
            [1, 1, 3, 11, 53, 309, 2119, 16687, 148329]

        AUTHORS:

        - Jaap Spies (2007-01-13)
        """
        SloaneSequence.__init__(self, offset=0)
        self._b = []
        self._a0a1d = (1, 1, 1)
        self._precompute(2)  # force precomputation, e.g. a(0) will fail when asked first

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000255._repr_()
            'a(n) = n*a(n-1) + (n-1)*a(n-2), a(0) = 1, a(1) = 1.'
        """
        return "a(n) = n*a(n-1) + (n-1)*a(n-2), a(0) = 1, a(1) = 1."


class A000261(ExtremesOfPermanentsSequence):
    def __init__(self):
        r"""
        `a(n) = n*a(n-1) + (n-3)*a(n-2)`, with `a(1) = 1`,
        `a(2) = 1`.

        With offset 1, permanent of (0,1)-matrix of size
        `n \times (n+d)` with `d=3` and `n` zeros
        not on a line. This is a special case of Theorem 2.3 of Seok-Zun
        Song et al. Extremes of permanents of (0,1)-matrices, p. 201-202.

        Seok-Zun Song et al., Extremes of permanents of (0,1)-matrices,
        Lin. Algebra and its Applic. 373 (2003), p. 197-210.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000261;a
            a(n) = n*a(n-1) + (n-3)*a(n-2), a(1) = 0, a(2) = 1.
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(1)
            0
            sage: a.offset
            1
            sage: a(8)
            30637
            sage: a(22)
            1801366114380914335441
            sage: a.list(9)
            [0, 1, 3, 13, 71, 465, 3539, 30637, 296967]

        AUTHORS:

        - Jaap Spies (2007-01-23)
        """
        SloaneSequence.__init__(self, offset=1)
        self._a0a1d = (0, 1, 3)
        self._b = []
        self._k = self.offset + 1

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000261._repr_()
            'a(n) = n*a(n-1) + (n-3)*a(n-2), a(1) = 0, a(2) = 1.'
        """
        return "a(n) = n*a(n-1) + (n-3)*a(n-2), a(1) = 0, a(2) = 1."


class A001909(ExtremesOfPermanentsSequence):
    def __init__(self):
        r"""
        `a(n) = n*a(n-1) + (n-4)*a(n-2)`, with `a(2) = 0`,
        `a(3) = 1`.

        With offset 1, permanent of (0,1)-matrix of size
        `n \times (n+d)` with `d=4` and `n` zeros
        not on a line. This is a special case of Theorem 2.3 of Seok-Zun
        Song et al. Extremes of permanents of (0,1)-matrices, p. 201-202.

        Seok-Zun Song et al., Extremes of permanents of (0,1)-matrices,
        Lin. Algebra and its Applic. 373 (2003), p. 197-210.

        INPUT:

        - ``n`` -- positive integer >= 2

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A001909;a
            a(n) = n*a(n-1) + (n-4)*a(n-2), a(2) = 0, a(3) = 1.
            sage: a(1)
            Traceback (most recent call last):
            ...
            ValueError: input n (=1) must be an integer >= 2
            sage: a.offset
            2
            sage: a(2)
            0
            sage: a(8)
            8544
            sage: a(22)
            470033715095287415734
            sage: a.list(9)
            [0, 1, 4, 21, 134, 1001, 8544, 81901, 870274]

        AUTHORS:

        - Jaap Spies (2007-01-13)
        """
        SloaneSequence.__init__(self, offset=2)
        self._a0a1d = (0, 1, 4)
        self._b = []
        self._k = self.offset + 1

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A001909._repr_()
            'a(n) = n*a(n-1) + (n-4)*a(n-2), a(2) = 0, a(3) = 1.'
        """
        return "a(n) = n*a(n-1) + (n-4)*a(n-2), a(2) = 0, a(3) = 1."


class A001910(ExtremesOfPermanentsSequence):
    def __init__(self):
        r"""
        `a(n) = n*a(n-1) + (n-5)*a(n-2)`, with `a(3) = 0`,
        `a(4) = 1`.

        With offset 1, permanent of (0,1)-matrix of size
        `n \times (n+d)` with `d=5` and `n` zeros
        not on a line. This is a special case of Theorem 2.3 of Seok-Zun
        Song et al. Extremes of permanents of (0,1)-matrices, p. 201-202.

        Seok-Zun Song et al., Extremes of permanents of (0,1)-matrices,
        Lin. Algebra and its Applic. 373 (2003), p. 197-210.

        INPUT:

        - ``n`` -- positive integer `\geq 3`

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A001910;a
            a(n) = n*a(n-1) + (n-5)*a(n-2), a(3) = 0, a(4) = 1.
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be an integer >= 3
            sage: a(3)
            0
            sage: a.offset
            3
            sage: a(8)
            1909
            sage: a(22)
            98125321641110663023
            sage: a.list(9)
            [0, 1, 5, 31, 227, 1909, 18089, 190435, 2203319]

        AUTHORS:

        - Jaap Spies (2007-01-13)
        """
        SloaneSequence.__init__(self, offset=3)
        self._a0a1d = (0, 1, 5)
        self._b = []
        self._k = self.offset + 1

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A001910._repr_()
            'a(n) = n*a(n-1) + (n-5)*a(n-2), a(3) = 0, a(4) = 1.'
        """
        return "a(n) = n*a(n-1) + (n-5)*a(n-2), a(3) = 0, a(4) = 1."


class ExtremesOfPermanentsSequence2(ExtremesOfPermanentsSequence):
    def gen(self, a0, a1, d):
        """
        EXAMPLES::

            sage: from sage.combinat.sloane_functions import ExtremesOfPermanentsSequence2
            sage: e = ExtremesOfPermanentsSequence2()
            sage: it = e.gen(6,43,6)
            sage: [next(it) for i in range(5)]
            [6, 43, 307, 2542, 23799]
        """
        x, y = ZZ(a0), ZZ(a1)
        k = self._k
        yield x
        while True:
            k = k+1
            x, y = y, (k-self._k1)*x+(k+d-self._k2)*y
            yield x
    _k1 = 1
    _k2 = 1


class A090010(ExtremesOfPermanentsSequence2):
    def __init__(self):
        r"""
        Permanent of (0,1)-matrix of size `n \times (n+d)` with
        `d=6` and `n` zeros not on a line.

        ` a(n) = (n+5)*a(n-1) + (n-1)*a(n-2), a(1)=6, a(2)=43`.

        This is a special case of Theorem 2.3 of Seok-Zun Song et al.
        Extremes of permanents of (0,1)-matrices, p. 201-202.

        REFERENCES:

        - Seok-Zun Song et al., Extremes of permanents of
          (0,1)-matrices, Lin. Algebra and its Applic. 373 (2003), p.
          197-210.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A090010;a
            Permanent of (0,1)-matrix of size n X (n+d) with d=6 and n zeros not on a line.
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(1)
            6
            sage: a(2)
            43
            sage: a.offset
            1
            sage: a(8)
            67741129
            sage: a(22)
            192416593029158989003270143
            sage: a.list(9)
            [6, 43, 356, 3333, 34754, 398959, 4996032, 67741129, 988344062]

        AUTHORS:

        - Jaap Spies (2007-01-23)
        """
        SloaneSequence.__init__(self, offset=1)
        self._a0a1d = (6, 43, 6)
        self._k = self.offset + 1
        self._b = []

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A090010._repr_()
            'Permanent of (0,1)-matrix of size n X (n+d) with d=6 and n zeros not on a line.'
        """
        return "Permanent of (0,1)-matrix of size n X (n+d) with d=6 and n zeros not on a line."


class A055790(ExtremesOfPermanentsSequence2):
    def __init__(self):
        r"""
        `a(n) = n*a(n-1) + (n-2)*a(n-2) [a(0) = 0, a(1) = 2]`.

        With offset 1, permanent of (0,1)-matrix of size n X (n+d) with d=1
        and n-1 zeros not on a line. This is a special case of Theorem 2.3
        of Seok-Zun Song et al. Extremes of permanents of (0,1)-matrices,
        p. 201-202.

        REFERENCES:

        - Seok-Zun Song et al., Extremes of permanents of
          (0,1)-matrices, Lin. Algebra and its Applic. 373 (2003), p.
          197-210.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A055790;a
            a(n) = n*a(n-1) + (n-2)*a(n-2) [a(0) = 0, a(1) = 2].
            sage: a(0)
            0
            sage: a(1)
            2
            sage: a(2)
            4
            sage: a.offset
            0
            sage: a(8)
            165016
            sage: a(22)
            10356214297533070441564
            sage: a.list(9)
            [0, 2, 4, 14, 64, 362, 2428, 18806, 165016]

        AUTHORS:

        - Jaap Spies (2007-01-23)
        """
        SloaneSequence.__init__(self, offset=0)
        self._a0a1d = (0, 2, 1)
        self._b = []
        self._precompute(2)
        self._k1 = 2

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A055790._repr_()
            'a(n) = n*a(n-1) + (n-2)*a(n-2) [a(0) = 0, a(1) = 2].'
        """
        return "a(n) = n*a(n-1) + (n-2)*a(n-2) [a(0) = 0, a(1) = 2]."


class A090012(SloaneSequence):
    def __init__(self):
        r"""
        Permanent of (0,1)-matrix of size `n \times (n+d)` with
        `d=2` and `n-1` zeros not on a line.

        `a(n) = (n+1)*a(n-1) + (n-2)*a(n-2)`, `a(1)=3` and
        `a(2)=9`

        This is a special case of Theorem 2.3 of Seok-Zun Song et al.
        Extremes of permanents of (0,1)-matrices, p. 201-202.

        REFERENCES:

        - Seok-Zun Song et al., Extremes of permanents of
          (0,1)-matrices, Lin. Algebra and its Applic. 373 (2003), p.
          197-210.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A090012;a
            Permanent of (0,1)-matrix of size n X (n+d) with d=2 and n-1 zeros not on a line.
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(1)
            3
            sage: a(2)
            9
            sage: a.offset
            1
            sage: a(8)
            890901
            sage: a(22)
            129020386652297208795129
            sage: a.list(9)
            [3, 9, 39, 213, 1395, 10617, 91911, 890901, 9552387]

        AUTHORS:

        - Jaap Spies (2007-01-23)
        """
        SloaneSequence.__init__(self, offset=1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A090012._repr_()
            'Permanent of (0,1)-matrix of size n X (n+d) with d=2 and n-1 zeros not on a line.'
        """
        return "Permanent of (0,1)-matrix of size n X (n+d) with d=2 and n-1 zeros not on a line."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A090012._eval(n) for n in range(1,11)]
            [3, 9, 39, 213, 1395, 10617, 91911, 890901, 9552387, 112203465]
        """
        if n == 1:
            return ZZ(3)
        return sloane.A000153(n+1) + sloane.A000153(n)


class A090013(SloaneSequence):
    def __init__(self):
        r"""
        Permanent of (0,1)-matrix of size `n \times (n+d)` with
        `d=3` and `n-1` zeros not on a line.

        `a(n) = (n+1)*a(n-1) + (n-2)*a(n-2) [a(1)=4, a(2)=16]`

        This is a special case of Theorem 2.3 of Seok-Zun Song et al.
        Extremes of permanents of (0,1)-matrices, p. 201-202.

        REFERENCES:

        - Seok-Zun Song et al., Extremes of permanents of
          (0,1)-matrices, Lin. Algebra and its Applic. 373 (2003), p.
          197-210.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A090013;a
            Permanent of (0,1)-matrix of size n X (n+d) with d=3 and n-1 zeros not on a line.
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(1)
            4
            sage: a(2)
            16
            sage: a.offset
            1
            sage: a(8)
            3481096
            sage: a(22)
            1112998577171142607670336
            sage: a.list(9)
            [4, 16, 84, 536, 4004, 34176, 327604, 3481096, 40585284]

        AUTHORS:

        - Jaap Spies (2007-01-23)
        """
        SloaneSequence.__init__(self, offset=1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A090013._repr_()
            'Permanent of (0,1)-matrix of size n X (n+d) with d=3 and n-1 zeros not on a line.'
        """
        return "Permanent of (0,1)-matrix of size n X (n+d) with d=3 and n-1 zeros not on a line."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A090013._eval(n) for n in range(1,11)]
            [4, 16, 84, 536, 4004, 34176, 327604, 3481096, 40585284, 514872176]
        """
        if n == 1:
            return ZZ(4)
        return sloane.A000261(n+2) + sloane.A000261(n+1)


class A090014(SloaneSequence):
    def __init__(self):
        r"""
        Permanent of (0,1)-matrix of size `n \times (n+d)` with
        `d=4` and `n-1` zeros not on a line.

        `a(n) = (n+1)*a(n-1) + (n-2)*a(n-2) [a(1)=5, a(2)=25]`

        This is a special case of Theorem 2.3 of Seok-Zun Song et al.
        Extremes of permanents of (0,1)-matrices, p. 201-202.

        REFERENCES:

        - Seok-Zun Song et al., Extremes of permanents of
          (0,1)-matrices, Lin. Algebra and its Applic. 373 (2003), p.
          197-210.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A090014;a
            Permanent of (0,1)-matrix of size n X (n+d) with d=4 and n-1 zeros not on a line.
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(1)
            5
            sage: a(2)
            25
            sage: a.offset
            1
            sage: a(8)
            11016595
            sage: a(22)
            7469733600354446865509725
            sage: a.list(9)
            [5, 25, 155, 1135, 9545, 90445, 952175, 11016595, 138864365]

        AUTHORS:

        - Jaap Spies (2007-01-23)
        """
        SloaneSequence.__init__(self, offset=1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A090014._repr_()
            'Permanent of (0,1)-matrix of size n X (n+d) with d=4 and n-1 zeros not on a line.'
        """
        return "Permanent of (0,1)-matrix of size n X (n+d) with d=4 and n-1 zeros not on a line."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A090014._eval(n) for n in range(1,11)]
            [5, 25, 155, 1135, 9545, 90445, 952175, 11016595, 138864365, 1893369505]
        """
        if n == 1:
            return ZZ(5)
        return sloane.A001909(n+3) + sloane.A001909(n+2)


class A090015(SloaneSequence):
    def __init__(self):
        r"""
        Permanent of (0,1)-matrix of size `n \times (n+d)` with
        `d=5` and `n-1` zeros not on a line.

        `a(n) = (n+1)*a(n-1) + (n-2)*a(n-2) [a(1)=6, a(2)=36]`

        This is a special case of Theorem 2.3 of Seok-Zun Song et al.
        Extremes of permanents of (0,1)-matrices, p. 201-202.

        REFERENCES:

        - Seok-Zun Song et al., Extremes of permanents of
          (0,1)-matrices, Lin. Algebra and its Applic. 373 (2003), p.
          197-210.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A090015;a
            Permanent of (0,1)-matrix of size n X (n+d) with d=3 and n-1 zeros not on a line.
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(1)
            6
            sage: a(2)
            36
            sage: a.offset
            1
            sage: a(8)
            29976192
            sage: a(22)
            41552258517692116794936876
            sage: a.list(9)
            [6, 36, 258, 2136, 19998, 208524, 2393754, 29976192, 406446774]

        AUTHORS:

        - Jaap Spies (2007-01-23)
        """
        SloaneSequence.__init__(self, offset=1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A090015._repr_()
            'Permanent of (0,1)-matrix of size n X (n+d) with d=3 and n-1 zeros not on a line.'
        """
        return "Permanent of (0,1)-matrix of size n X (n+d) with d=3 and n-1 zeros not on a line."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A090015._eval(n) for n in range(1,10)]
            [6, 36, 258, 2136, 19998, 208524, 2393754, 29976192, 406446774]
        """
        if n == 1:
            return ZZ(6)
        return sloane.A001910(n+4) + sloane.A001910(n+3)


class A090016(SloaneSequence):
    def __init__(self):
        r"""
        Permanent of (0,1)-matrix of size `n \times (n+d)` with
        `d=6` and `n-1` zeros not on a line.

        `a(n) = (n+1)*a(n-1) + (n-2)*a(n-2) [a(1)=7, a(2)=49]`

        `A090016 a(n) = A090010(n-1) + A090010(n), a(1)=7`

        This is a special case of Theorem 2.3 of Seok-Zun Song et al.
        Extremes of permanents of (0,1)-matrices, p. 201-202.

        REFERENCES:

        - Seok-Zun Song et al., Extremes of permanents of
          (0,1)-matrices, Lin. Algebra and its Applic. 373 (2003), p.
          197-210.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A090016;a
            Permanent of (0,1)-matrix of size n X (n+d) with d=6 and n-1 zeros not on a line.
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(1)
            7
            sage: a(2)
            49
            sage: a.offset
            1
            sage: a(8)
            72737161
            sage: a(22)
            199341969448774341802426289
            sage: a.list(9)
            [7, 49, 399, 3689, 38087, 433713, 5394991, 72737161, 1056085191]

        AUTHORS:

        - Jaap Spies (2007-01-23)
        """
        SloaneSequence.__init__(self, offset=1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A090016._repr_()
            'Permanent of (0,1)-matrix of size n X (n+d) with d=6 and n-1 zeros not on a line.'
        """
        return "Permanent of (0,1)-matrix of size n X (n+d) with d=6 and n-1 zeros not on a line."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A090016._eval(n) for n in range(1,10)]
            [7, 49, 399, 3689, 38087, 433713, 5394991, 72737161, 1056085191]
        """
        if n == 1:
            return ZZ(7)
        return sloane.A090010(n-1) + sloane.A090010(n)


class A000166(SloaneSequence):
    def __init__(self):
        r"""
        Subfactorial or rencontres numbers, or derangements: number of
        permutations of `n` elements with no fixed points.

        With offset 1 also the permanent of a (0,1)-matrix of order
        `n` with `n` 0s not on a line.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000166;a
            Subfactorial or rencontres numbers, or derangements: number of permutations of $n$ elements with no fixed points.
            sage: a(0)
            1
            sage: a(1)
            0
            sage: a(2)
            1
            sage: a.offset
            0
            sage: a(8)
            14833
            sage: a(20)
            895014631192902121
            sage: a.list(9)
            [1, 0, 1, 2, 9, 44, 265, 1854, 14833]

        AUTHORS:

        - Jaap Spies (2007-01-13)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000166._repr_()
            'Subfactorial or rencontres numbers, or derangements: number of permutations of $n$ elements with no fixed points.'
        """
        return "Subfactorial or rencontres numbers, or derangements: number of permutations of $n$ elements with no fixed points."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000166._eval(n) for n in range(9)]
            [1, 0, 1, 2, 9, 44, 265, 1854, 14833]
        """
        return arith.subfactorial(n)


class A000203(SloaneSequence):
    def __init__(self):
        r"""
        The sequence `\sigma(n)`, where `\sigma(n)` is the
        sum of the divisors of `n`. Also called
        `\sigma_1(n)`.

        The function ``sigma(n, k)`` implements
        `\sigma_k(n)` in Sage.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000203; a
            sigma(n) = sum of divisors of n. Also called sigma_1(n).
            sage: a(1)
            1
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(256)
            511
            sage: a.list(12)
            [1, 3, 4, 7, 6, 12, 8, 15, 13, 18, 12, 28]
            sage: a(1/3)
            Traceback (most recent call last):
            ...
            TypeError: input must be an int or Integer

        AUTHORS:

        - Jaap Spies (2007-01-13)
        """
        SloaneSequence.__init__(self, offset=1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000203._repr_()
            'sigma(n) = sum of divisors of n. Also called sigma_1(n).'
        """
        return "sigma(n) = sum of divisors of n. Also called sigma_1(n)."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000203._eval(n) for n in range(1,11)]
            [1, 3, 4, 7, 6, 12, 8, 15, 13, 18]
        """
        return sum(arith.divisors(n))  # alternative: return arith.sigma(n)


class A001157(SloaneSequence):
    def __init__(self):
        r"""
        The sequence `\sigma_2(n)`, sum of squares of divisors of
        `n`.

        The function sigma(n, k) implements `\sigma_k*` in Sage.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A001157;a
            sigma_2(n): sum of squares of divisors of n
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(2)
            5
            sage: a(8)
            85
            sage: a.list(9)
            [1, 5, 10, 21, 26, 50, 50, 85, 91]

        AUTHORS:

        - Jaap Spies (2007-01-13)
        """
        SloaneSequence.__init__(self, offset=1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A001157._repr_()
            'sigma_2(n): sum of squares of divisors of n'
        """
        return "sigma_2(n): sum of squares of divisors of n"

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A001157._eval(n) for n in range(1,11)]
            [1, 5, 10, 21, 26, 50, 50, 85, 91, 130]
        """
        return arith.sigma(n, 2)


class A008683(SloaneSequence):
    def __init__(self):
        r"""
        Möbius function `\mu(n)`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A008683;a
            Moebius function mu(n).
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(2)
            -1
            sage: a(12)
            0
            sage: a.list(12)
            [1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0]

        AUTHORS:

        - Jaap Spies (2007-01-13)
        """
        SloaneSequence.__init__(self, offset=1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A008683._repr_()
            'Moebius function mu(n).'
        """
        return "Moebius function mu(n)."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A008683._eval(n) for n in range(1,11)]
            [1, -1, -1, 0, -1, 1, -1, 0, 0, 1]
        """
        return arith.moebius(n)


class A000204(SloaneSequence):
    def __init__(self):
        r"""
        Lucas numbers (beginning with 1): `L(n) = L(n-1) + L(n-2)`
        with `L(1) = 1`, `L(2) = 3`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000204; a
            Lucas numbers (beginning at 1): L(n) = L(n-1) + L(n-2), L(2) = 3.
            sage: a(1)
            1
            sage: a(8)
            47
            sage: a(200)
            627376215338105766356982006981782561278127
            sage: a(-4)
            Traceback (most recent call last):
            ...
            ValueError: input n (=-4) must be a positive integer
            sage: a.list(12)
            [1, 3, 4, 7, 11, 18, 29, 47, 76, 123, 199, 322]
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer

        AUTHORS:

        - Jaap Spies (2007-01-18)
        """
        SloaneSequence.__init__(self, offset=1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000204._repr_()
            'Lucas numbers (beginning at 1): L(n) = L(n-1) + L(n-2), L(2) = 3.'
        """
        return "Lucas numbers (beginning at 1): L(n) = L(n-1) + L(n-2), L(2) = 3."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000204._eval(n) for n in range(1,11)]
            [1, 3, 4, 7, 11, 18, 29, 47, 76, 123]
        """
        if n == 1:
            return ZZ.one()
        elif n == 2:
            return 3
        return sloane.A000045(n+1) + sloane.A000045(n-1)


class A000217(SloaneSequence):
    def __init__(self):
        r"""
        Triangular numbers: `a(n) = \binom{n+1}{2} = n(n+1)/2`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000217;a
            Triangular numbers: a(n) = C(n+1,2) = n(n+1)/2 = 0+1+2+...+n.
            sage: a(0)
            0
            sage: a(2)
            3
            sage: a(8)
            36
            sage: a(2000)
            2001000
            sage: a.list(9)
            [0, 1, 3, 6, 10, 15, 21, 28, 36]

        AUTHORS:

        - Jaap Spies (2007-01-25)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000217._repr_()
            'Triangular numbers: a(n) = C(n+1,2) = n(n+1)/2 = 0+1+2+...+n.'
        """
        return "Triangular numbers: a(n) = C(n+1,2) = n(n+1)/2 = 0+1+2+...+n."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000217._eval(n) for n in range(10)]
            [0, 1, 3, 6, 10, 15, 21, 28, 36, 45]
        """
        return ZZ(n*(n+1)//2)


class A000124(SloaneSequence):
    def __init__(self):
        r"""
        Central polygonal numbers (the Lazy Caterer's sequence):
        `n(n+1)/2 + 1`.

        Or, maximal number of pieces formed when slicing a pancake with
        `n` cuts.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000124;a
            Central polygonal numbers (the Lazy Caterer's sequence): n(n+1)/2 + 1.
            sage: a(0)
            1
            sage: a(1)
            2
            sage: a(2)
            4
            sage: a(9)
            46
            sage: a.list(10)
            [1, 2, 4, 7, 11, 16, 22, 29, 37, 46]

        AUTHORS:

        - Jaap Spies (2007-01-25)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000124._repr_()
            "Central polygonal numbers (the Lazy Caterer's sequence): n(n+1)/2 + 1."
        """
        return "Central polygonal numbers (the Lazy Caterer's sequence): n(n+1)/2 + 1."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000124._eval(n) for n in range(10)]
            [1, 2, 4, 7, 11, 16, 22, 29, 37, 46]
        """
        return ZZ(n*(n+1)//2 + 1)


class A002275(SloaneSequence):
    def __init__(self):
        r"""
        Repunits: `\frac {(10^n - 1)}{9}`. Often denoted by
        `R_n`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A002275;a
            Repunits: (10^n - 1)/9. Often denoted by R_n.
            sage: a(0)
            0
            sage: a(2)
            11
            sage: a(8)
            11111111
            sage: a(20)
            11111111111111111111
            sage: a.list(9)
            [0, 1, 11, 111, 1111, 11111, 111111, 1111111, 11111111]

        AUTHORS:

        - Jaap Spies (2007-01-25)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A002275._repr_()
            'Repunits: (10^n - 1)/9. Often denoted by R_n.'
        """
        return "Repunits: (10^n - 1)/9. Often denoted by R_n."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A002275._eval(n) for n in range(10)]
            [0, 1, 11, 111, 1111, 11111, 111111, 1111111, 11111111, 111111111]
        """
        return ZZ(10**n-1)//9


# inhomogeneous second order recurrences
def recur_gen2b(a0, a1, a2, a3, b):
    r"""
    Inhomogeneous second-order linear recurrence generator with fixed
    coefficients and `b = f(n)`

    `a(0) = a0`, `a(1) = a1`,
    `a(n) = a2*a(n-1) + a3*a(n-2) +f(n)`.

    EXAMPLES::

        sage: from sage.combinat.sloane_functions import recur_gen2b
        sage: it = recur_gen2b(1,1,1,1, lambda n: 0)
        sage: [next(it) for i in range(10)]
        [1, 1, 2, 3, 5, 8, 13, 21, 34, 55]
    """
    x, y = ZZ(a0), ZZ(a1)
    n = 1
    yield x
    while True:
        n = n+1
        x, y = y, a3*x+a2*y + b(n)
        yield x


class RecurrenceSequence(SloaneSequence):
    def _precompute(self, how_many=20):
        """
        EXAMPLES::

            sage: initial = len(sloane.A001110._b)
            sage: sloane.A001110._precompute(10)
            sage: len(sloane.A001110._b) - initial == 10
            True
        """
        try:
            f = self._f
        except AttributeError:
            self._f = recur_gen2b(*self._params)
            f = self._f
        self._b += [next(f) for i in range(how_many)]

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A001110._eval(n) for n in range(5)]
            [0, 1, 36, 1225, 41616]
        """
        if len(self._b) < n:
            self._precompute(n - len(self._b) + 1)
        return self._b[n]

    def list(self, n):
        """
        EXAMPLES::

            sage: sloane.A001110.list(8)
            [0, 1, 36, 1225, 41616, 1413721, 48024900, 1631432881]
        """
        self._eval(n)   # force computation
        return self._b[:n]


class A001110(RecurrenceSequence):
    def __init__(self):
        r"""
        Numbers that are both triangular and square:
        `a(n) = 34a(n-1) - a(n-2) + 2`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A001110; a
            Numbers that are both triangular and square: a(n) = 34a(n-1) - a(n-2) + 2.
            sage: a(0)
            0
            sage: a(1)
            1
            sage: a(8)
            55420693056
            sage: a(21)
            4446390382511295358038307980025
            sage: a.list(8)
            [0, 1, 36, 1225, 41616, 1413721, 48024900, 1631432881]

        AUTHORS:

        - Jaap Spies (2007-01-19)
        """
        SloaneSequence.__init__(self, offset=0)
        self._params = (0, 1, 34, -1, self.g)
        self._b = []
        self._precompute()

    link = "http://oeis.org/classic/A001110"

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A001110._repr_()
            'Numbers that are both triangular and square: a(n) = 34a(n-1) - a(n-2) + 2.'
        """
        return "Numbers that are both triangular and square: a(n) = 34a(n-1) - a(n-2) + 2."

    def g(self, k):
        """
        EXAMPLES::

            sage: sloane.A001110.g(2)
            2
            sage: sloane.A001110.g(1)
            0
        """
        if k > 1:
            return 2
        else:
            return ZZ.zero()


class A051959(RecurrenceSequence):
    def __init__(self):
        r"""
        Linear second order recurrence. A051959.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A051959; a
            Linear second order recurrence. A051959.
            sage: a(0)
            1
            sage: a(1)
            10
            sage: a(8)
            9969
            sage: a(41)
            42834431872413650
            sage: a.list(12)
            [1, 10, 36, 104, 273, 686, 1688, 4112, 9969, 24114, 58268, 140728]

        AUTHORS:

        - Jaap Spies (2007-01-19)
        """
        SloaneSequence.__init__(self, offset=0)
        self._params = (1,10,2,1,self.g)
        self._b = []
        self._precompute(2)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A051959._repr_()
            'Linear second order recurrence. A051959.'
        """
        return "Linear second order recurrence. A051959."

    def g(self, k):
        """
        EXAMPLES::

            sage: sloane.A051959.g(2)
            15
            sage: sloane.A051959.g(1)
            0
        """
        if k > 1:
            return 7*k+1
        else:
            return ZZ.zero()


class A001221(SloaneSequence):
    def __init__(self):
        r"""
        Number of different prime divisors of `n`.

        Also called omega(n) or `\omega(n)`. Maximal number of
        terms in any factorization of `n`. Number of prime powers
        that divide `n`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A001221; a
            Number of distinct primes dividing n (also called omega(n)).
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(1)
            0
            sage: a(8)
            1
            sage: a(41)
            1
            sage: a(84792)
            3
            sage: a.list(12)
            [0, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 2]

        AUTHORS:

        - Jaap Spies (2007-01-19)
        """
        SloaneSequence.__init__(self, offset=1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A001221._repr_()
            'Number of distinct primes dividing n (also called omega(n)).'
        """
        return "Number of distinct primes dividing n (also called omega(n))."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A001221._eval(n) for n in range(1,10)]
            [0, 1, 1, 1, 1, 2, 1, 1, 1]
        """
        return len(arith.prime_divisors(n))  # there is a PARI function omega


class A001222(SloaneSequence):
    def __init__(self):
        r"""
        Number of prime divisors of `n` (counted with
        multiplicity).

        Also called bigomega(n) or `\Omega(n)`. Maximal number of
        terms in any factorization of `n`. Number of prime powers
        that divide `n`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A001222; a
            Number of prime divisors of n (counted with multiplicity).
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(1)
            0
            sage: a(8)
            3
            sage: a(41)
            1
            sage: a(84792)
            5
            sage: a.list(12)
            [0, 1, 1, 2, 1, 2, 1, 3, 2, 2, 1, 3]

        AUTHORS:

        - Jaap Spies (2007-01-19)
        """
        SloaneSequence.__init__(self, offset=1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A001222._repr_()
            'Number of prime divisors of n (counted with multiplicity).'
        """
        return "Number of prime divisors of n (counted with multiplicity)."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A001222._eval(n) for n in range(1,10)]
            [0, 1, 1, 2, 1, 2, 1, 3, 2]
        """
        return sum(e for i,e in arith.factor(n))


# A046660() = A001222(n) - A001221(n)
class A046660(SloaneSequence):
    r"""
    Excess of `n` = number of prime divisors (with
    multiplicity) - number of prime divisors (without multiplicity).

    `\Omega(n) - \omega(n)`.

    INPUT:

    - ``n`` -- positive integer

    OUTPUT: integer; function value

    EXAMPLES::

        sage: a = sloane.A046660; a
        Excess of n = Bigomega (with multiplicity) - omega (without multiplicity).
        sage: a(0)
        Traceback (most recent call last):
        ...
        ValueError: input n (=0) must be a positive integer
        sage: a(1)
        0
        sage: a(8)
        2
        sage: a(41)
        0
        sage: a(84792)
        2
        sage: a.list(12)
        [0, 0, 0, 1, 0, 0, 0, 2, 1, 0, 0, 1]

    AUTHORS:

        - Jaap Spies (2007-01-19)
    """

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A046660._repr_()
            'Excess of n = Bigomega (with multiplicity) - omega (without multiplicity).'
        """
        return "Excess of n = Bigomega (with multiplicity) - omega (without multiplicity)."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A046660._eval(n) for n in range(1,10)]
            [0, 0, 0, 1, 0, 0, 0, 2, 1]
        """
        return sloane.A001222(n) - sloane.A001221(n)


class A001227(SloaneSequence):
    def __init__(self):
        r"""
        Number of odd divisors of `n`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A001227; a
            Number of odd divisors of n
            sage: a.offset
            1
            sage: a(1)
            1
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(100)
            3
            sage: a(256)
            1
            sage: a(29)
            2
            sage: a.list(20)
            [1, 1, 2, 1, 2, 2, 2, 1, 3, 2, 2, 2, 2, 2, 4, 1, 2, 3, 2, 2]
            sage: a(-1)
            Traceback (most recent call last):
            ...
            ValueError: input n (=-1) must be a positive integer

        AUTHORS:

        - Jaap Spies (2007-01-14)
        """
        SloaneSequence.__init__(self, offset=1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A001227._repr_()
            'Number of odd divisors of n'
        """
        return "Number of odd divisors of n"

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A001227._eval(n) for n in range(1,10)]
            [1, 1, 2, 1, 2, 2, 2, 1, 3]
        """
        return sum(i % 2 for i in arith.divisors(n))


class A001358(SloaneSequence):
    def __init__(self):
        r"""
        Products of two primes.

        These numbers have been called semiprimes (or semi-primes),
        biprimes or 2-almost primes.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A001358;a
            Products of two primes.
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(2)
            6
            sage: a(8)
            22
            sage: a(200)
            669
            sage: a.list(9)
            [4, 6, 9, 10, 14, 15, 21, 22, 25]

        AUTHORS:

        - Jaap Spies (2007-01-25)
        """
        SloaneSequence.__init__(self, offset=1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A001358._repr_()
            'Products of two primes.'
        """
        return "Products of two primes."

    def _precompute(self, how_many=150):
        """
        EXAMPLES::

            sage: initial = len(sloane.A001358._b)
            sage: sloane.A001358._precompute()
            sage: len(sloane.A001358._b) - initial > 0
            True
        """
        try:
            self._b
            self._n
        except AttributeError:
            self._b = []
            self._n = 1
        self._b += [i for i in range(self._n, self._n + how_many)
                    if sum(e for _, e in arith.factor(i)) == 2]
        self._n += how_many

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A001358._eval(n) for n in range(1,10)]
            [4, 6, 9, 10, 14, 15, 21, 22, 25]
        """
        try:
            return self._b[n-1]
        except (AttributeError, IndexError):
            self._precompute()
            # try again
            return self._eval(n)

    def list(self, n):
        """
        EXAMPLES::

            sage: sloane.A001358.list(9)
            [4, 6, 9, 10, 14, 15, 21, 22, 25]
        """
        try:
            if len(self._b) < n:
                raise IndexError
            else:
                return self._b[:n]
        except (AttributeError, IndexError):
            self._precompute()
            # try again
            return self.list(n)


class A001694(SloaneSequence):
    def __init__(self):
        r"""
        This function returns the `n`-th Powerful Number:

        A positive integer `n` is powerful if for every prime
        `p` dividing `n`, `p^2` also divides
        `n`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A001694; a
            Powerful Numbers (also called squarefull, square-full or 2-full numbers).
            sage: a.offset
            1
            sage: a(1)
            1
            sage: a(4)
            9
            sage: a(100)
            3136
            sage: a(156)
            7225
            sage: a.list(19)
            [1, 4, 8, 9, 16, 25, 27, 32, 36, 49, 64, 72, 81, 100, 108, 121, 125, 128, 144]
            sage: a(-1)
            Traceback (most recent call last):
            ...
            ValueError: input n (=-1) must be a positive integer

        AUTHORS:

        - Jaap Spies (2007-01-14)
        """
        SloaneSequence.__init__(self, offset=1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A001694._repr_()
            'Powerful Numbers (also called squarefull, square-full or 2-full numbers).'
        """
        return "Powerful Numbers (also called squarefull, square-full or 2-full numbers)."

    def _precompute(self, how_many=10000):
        """
        EXAMPLES::

            sage: initial = len(sloane.A001694._b)
            sage: sloane.A001694._precompute()
            sage: len(sloane.A001694._b) - initial > 0
            True
        """
        try:
            self._b
            self._n
        except AttributeError:
            self._b = [1]
            self._n = 1
        self._b += self._powerful_numbers_in_range(self._n, self._n + how_many)
        self._n += how_many

    def _powerful_numbers_in_range(self, n, m):
        """
        EXAMPLES::

            sage: sloane.A001694._powerful_numbers_in_range(0,50)
            [4, 8, 9, 16, 25, 27, 32, 36, 49]
        """
        # This is naive -- too slow; too much overhead
        #  return [i for i in range(self._n, self._n+how_many) if self.is_powerful(i)]

        n = max(n, 4)
        # Use PARI directly -- much faster.
        from sage.libs.pari import pari
        L = pari('v=listcreate(); for(i=%s,%s,if(vecmin(factor(i)[,2])>1,listput(v,i))); v' % (n, m))
        return [ZZ(x) for x in L]  # not very many, so not much overhead

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A001694._eval(n) for n in range(1,10)]
            [1, 4, 8, 9, 16, 25, 27, 32, 36]
        """
        try:
            return self._b[n-1]
        except AttributeError:
            self._b = [1]
        except IndexError:
            pass
        while len(self._b) < n:
            self._precompute(10000)
        # try again, but we could also return self._b[n-1]
        return self._eval(n)

    def list(self, n):
        """
        EXAMPLES::

            sage: sloane.A001694.list(9)
            [1, 4, 8, 9, 16, 25, 27, 32, 36]
        """
        try:
            if len(self._b) < n:
                raise IndexError
            else:
                return self._b[:n]
        except AttributeError:
            self._b = [1]
        except IndexError:
            pass
        while len(self._b) < n:
            self._precompute(10000)
        return self._b[:n]

    def is_powerful(self, n) -> bool:
        r"""
        Return ``True`` if and only if `n` is a powerful number.

        A positive integer `n` is powerful if for every prime `p`
        dividing `n`, `p^2` also divides `n`.

        See :oeis:`A001694`.

        INPUT:

        - ``n`` -- integer

        OUTPUT: ``True`` if `n` is a powerful number, else ``False``

        EXAMPLES::

            sage: a = sloane.A001694
            sage: a.is_powerful(2500)
            True
            sage: a.is_powerful(20)
            False

        AUTHORS:

        - Jaap Spies (2006-12-07)
        """
        if n <= 1:
            return True
        return all(e >= 2 for _, e in arith.factor(n))


class A001836(SloaneSequence):
    def __init__(self):
        r"""
        Numbers `n` such that `\phi(2n-1) < \phi(2n)`,
        where `\phi` is Euler's totient function.

        Euler's totient function is also known as ``euler_phi``; ``euler_phi``
        is a standard Sage function.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A001836; a
            Numbers n such that phi(2n-1) < phi(2n), where phi is Euler's totient function A000010.
            sage: a.offset
            1
            sage: a(1)
            53
            sage: a(8)
            683
            sage: a(300)
            17798
            sage: a.list(12)
            [53, 83, 158, 263, 293, 368, 578, 683, 743, 788, 878, 893]
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer

        Compare: Searching Sloane's online database... Numbers n such that
        phi(2n-1) phi(2n), where phi is Euler's totient function A000010.
        [53, 83, 158, 263, 293, 368, 578, 683, 743, 788, 878, 893]

        AUTHORS:

        - Jaap Spies (2007-01-17)
        """
        SloaneSequence.__init__(self, offset=1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A001836._repr_()
            "Numbers n such that phi(2n-1) < phi(2n), where phi is Euler's totient function A000010."
        """
        return "Numbers n such that phi(2n-1) < phi(2n), where phi is Euler's totient function A000010."

    def _precompute(self, how_many=150):
        """
        EXAMPLES::

            sage: initial = len(sloane.A001836._b)
            sage: sloane.A001836._precompute()
            sage: len(sloane.A001836._b) - initial > 0
            True
        """
        try:
            self._b
            self._n
        except AttributeError:
            self._b = []
            self._n = self.offset
        self._b += [i for i in range(self._n, self._n + how_many)
                    if arith.euler_phi(2 * i - 1) < arith.euler_phi(2 * i)]
        self._n += how_many

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A001836._eval(n) for n in range(1,10)]
            [53, 83, 158, 263, 293, 368, 578, 683, 743]
        """
        try:
            return self._b[n-1]
        except (AttributeError, IndexError):
            self._precompute()
            # try again
            return self._eval(n)

    def list(self, n):
        """
        EXAMPLES::

            sage: sloane.A001836.list(9)
            [53, 83, 158, 263, 293, 368, 578, 683, 743]
        """
        try:
            if len(self._b) < n:
                raise IndexError
            else:
                return self._b[:n]
        except (AttributeError, IndexError):
            self._precompute()
            # try again
            return self.list(n)


# a group of sequences uses this function:
def recur_gen2(a0, a1, a2, a3):
    r"""
    Homogeneous general second-order linear recurrence generator with
    fixed coefficients.

    a(0) = a0, a(1) = a1, a(n) = a2\*a(n-1) + a3\*a(n-2)

    EXAMPLES::

        sage: from sage.combinat.sloane_functions import recur_gen2
        sage: it = recur_gen2(1,1,1,1)
        sage: [next(it) for i in range(10)]
        [1, 1, 2, 3, 5, 8, 13, 21, 34, 55]
    """
    x, y = ZZ(a0), ZZ(a1)
    n = 0
    yield x
    while True:
        n = n+1
        x, y = y, a3*x+a2*y
        yield x


# A001906 = recur_gen2(0,1,3,-1)
# This can be done much more simple: return sloane.A000045(2*n).
# but this is a proof of technology!
class RecurrenceSequence2(SloaneSequence):
    def _precompute(self, how_many=150):
        """
        EXAMPLES::

            sage: initial = len(sloane.A001906._b)
            sage: sloane.A001906._precompute(10)
            sage: len(sloane.A001906._b) - initial == 10
            True
        """
        try:
            f = self._f
        except AttributeError:
            self._f = recur_gen2(*self._params)
            f = self._f
        self._b += [next(f) for i in range(how_many)]

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A001906._eval(n) for n in range(10)]
            [0, 1, 3, 8, 21, 55, 144, 377, 987, 2584]
        """
        if len(self._b) <= n:
            self._precompute(n - len(self._b) + 1)
        return self._b[n]

    def list(self, n):
        """
        EXAMPLES::

            sage: sloane.A001906.list(10)
            [0, 1, 3, 8, 21, 55, 144, 377, 987, 2584]
        """
        self._eval(n)   # force computation
        return self._b[:n]


class A001906(RecurrenceSequence2):
    def __init__(self):
        r"""
        `F(2n) =` bisection of Fibonacci sequence:
        `a(n)=3a(n-1)-a(n-2)`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A001906; a
            F(2n) = bisection of Fibonacci sequence: a(n)=3a(n-1)-a(n-2).
            sage: a(0)
            0
            sage: a(1)
            1
            sage: a(8)
            987
            sage: a(22)
            701408733
            sage: a.list(12)
            [0, 1, 3, 8, 21, 55, 144, 377, 987, 2584, 6765, 17711]

        AUTHORS:

        - Jaap Spies (2007-01-19)
        """
        SloaneSequence.__init__(self, offset=0)
        self._params = (0,1,3,-1)
        self._b = []
        self._precompute(2)  # force precomputation

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A001906._repr_()
            'F(2n) = bisection of Fibonacci sequence: a(n)=3a(n-1)-a(n-2).'
        """
        return "F(2n) = bisection of Fibonacci sequence: a(n)=3a(n-1)-a(n-2)."


class A001333(RecurrenceSequence2):
    def __init__(self):
        r"""
        Numerators of continued fraction convergents to `\sqrt 2`.

        See also A000129

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A001333;a
            Numerators of continued fraction convergents to sqrt(2).
            sage: a(0)
            1
            sage: a(1)
            1
            sage: a(2)
            3
            sage: a(3)
            7
            sage: a(11)
            8119
            sage: a.list(12)
            [1, 1, 3, 7, 17, 41, 99, 239, 577, 1393, 3363, 8119]

        AUTHORS:

        - Jaap Spies (2007-02-01)
        """
        SloaneSequence.__init__(self, offset=0)
        self._b = []
        self._params = (1,1,2,1)
        self._precompute(2)  # force precomputation

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A001333._repr_()
            'Numerators of continued fraction convergents to sqrt(2).'
        """
        return "Numerators of continued fraction convergents to sqrt(2)."


class A001045(RecurrenceSequence2):
    def __init__(self):
        r"""
        Jacobsthal sequence: `a(n) = a(n-1) + 2a(n-2)`,
        `a(0) = 0` and `a(1) = 1`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A001045;a
            Jacobsthal sequence: a(n) = a(n-1) + 2a(n-2).
            sage: a(0)
            0
            sage: a(1)
            1
            sage: a(2)
            1
            sage: a(11)
            683
            sage: a.list(12)
            [0, 1, 1, 3, 5, 11, 21, 43, 85, 171, 341, 683]

        AUTHORS:

        - Jaap Spies (2007-01-26)
        """
        SloaneSequence.__init__(self, offset=0)
        self._params = (0,1,1,2)
        self._b = []
        self._precompute(2)  # force precomputation

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A001045._repr_()
            'Jacobsthal sequence: a(n) = a(n-1) + 2a(n-2).'
        """
        return "Jacobsthal sequence: a(n) = a(n-1) + 2a(n-2)."


class A000129(RecurrenceSequence2):
    def __init__(self):
        r"""
        Pell numbers: `a(0) = 0`, `a(1) = 1`; for
        `n > 1`, `a(n) = 2a(n-1) + a(n-2)`.

        Denominators of continued fraction convergents to
        `\sqrt 2`.

        See also A001333

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000129;a
            Pell numbers: a(0) = 0, a(1) = 1; for n > 1, a(n) = 2*a(n-1) + a(n-2).
            sage: a(0)
            0
            sage: a(2)
            2
            sage: a(12)
            13860
            sage: a.list(12)
            [0, 1, 2, 5, 12, 29, 70, 169, 408, 985, 2378, 5741]

        AUTHORS:

        - Jaap Spies (2007-01-25)
        """
        SloaneSequence.__init__(self, offset=0)
        self._b = []
        self._params = (0,1,2,1)
        self._precompute(2)  # force precomputation

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000129._repr_()
            'Pell numbers: a(0) = 0, a(1) = 1; for n > 1, a(n) = 2*a(n-1) + a(n-2).'
        """
        return "Pell numbers: a(0) = 0, a(1) = 1; for n > 1, a(n) = 2*a(n-1) + a(n-2)."


class A001109(RecurrenceSequence2):
    def __init__(self):
        r"""
        `a(n)^2` is a triangular number:
        `a(n) = 6*a(n-1) - a(n-2)` with `a(0)=0`,
        `a(1)=1`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A001109;a
            a(n)^2 is a triangular number: a(n) = 6*a(n-1) - a(n-2) with a(0)=0, a(1)=1
            sage: a(0)
            0
            sage: a(1)
            1
            sage: a(2)
            6
            sage: a.offset
            0
            sage: a(8)
            235416
            sage: a(60)
            1515330104844857898115857393785728383101709300
            sage: a.list(9)
            [0, 1, 6, 35, 204, 1189, 6930, 40391, 235416]

        AUTHORS:

        - Jaap Spies (2007-01-24)
        """
        SloaneSequence.__init__(self, offset=0)
        self._params = (0,1,6,-1)
        self._b = []
        self._precompute(2)  # force precomputation

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A001109._repr_()
            'a(n)^2 is a triangular number: a(n) = 6*a(n-1) - a(n-2) with a(0)=0, a(1)=1'
        """
        return "a(n)^2 is a triangular number: a(n) = 6*a(n-1) - a(n-2) with a(0)=0, a(1)=1"


class A015521(RecurrenceSequence2):
    def __init__(self):
        r"""
        Linear 2nd order recurrence, `a(0)=0`, `a(1)=1` and
        `a(n) = 3 a(n-1) + 4 a(n-2)`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A015521; a
            Linear 2nd order recurrence, a(n) = 3 a(n-1) + 4 a(n-2).
            sage: a(0)
            0
            sage: a(1)
            1
            sage: a(8)
            13107
            sage: a(41)
            967140655691703339764941
            sage: a.list(12)
            [0, 1, 3, 13, 51, 205, 819, 3277, 13107, 52429, 209715, 838861]

        AUTHORS:

        - Jaap Spies (2007-01-19)
        """
        SloaneSequence.__init__(self, offset=0)
        self._params = (0,1,3,4)
        self._b = []
        self._precompute(2)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A015521._repr_()
            'Linear 2nd order recurrence, a(n) = 3 a(n-1) + 4 a(n-2).'
        """
        return "Linear 2nd order recurrence, a(n) = 3 a(n-1) + 4 a(n-2)."


class A015523(RecurrenceSequence2):
    def __init__(self):
        r"""
        Linear 2nd order recurrence, `a(0)=0`, `a(1)=1` and
        `a(n) = 3 a(n-1) + 5 a(n-2)`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A015523; a
            Linear 2nd order recurrence, a(n) = 3 a(n-1) + 5 a(n-2).
            sage: a(0)
            0
            sage: a(1)
            1
            sage: a(8)
            17727
            sage: a(41)
            6173719566474529739091481
            sage: a.list(12)
            [0, 1, 3, 14, 57, 241, 1008, 4229, 17727, 74326, 311613, 1306469]

        AUTHORS:

        - Jaap Spies (2007-01-19)
        """
        SloaneSequence.__init__(self, offset=0)
        self._params = (0,1,3,5)
        self._b = []
        self._precompute(2)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A015523._repr_()
            'Linear 2nd order recurrence, a(n) = 3 a(n-1) + 5 a(n-2).'
        """
        return "Linear 2nd order recurrence, a(n) = 3 a(n-1) + 5 a(n-2)."


class A015530(RecurrenceSequence2):
    def __init__(self):
        r"""
        Linear 2nd order recurrence, `a(0)=0`, `a(1)=1` and
        `a(n) = 4 a(n-1) + 3 a(n-2)`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A015530;a
            Linear 2nd order recurrence, a(n) = 4 a(n-1) + 3 a(n-2).
            sage: a(0)
            0
            sage: a(1)
            1
            sage: a(2)
            4
            sage: a.offset
            0
            sage: a(8)
            41008
            sage: a.list(9)
            [0, 1, 4, 19, 88, 409, 1900, 8827, 41008]

        AUTHORS:

        - Jaap Spies (2007-01-19)
        """
        SloaneSequence.__init__(self, offset=0)
        self._b = []
        self._params = (0,1,4,3)
        self._precompute(2)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A015530._repr_()
            'Linear 2nd order recurrence, a(n) = 4 a(n-1) + 3 a(n-2).'
        """
        return "Linear 2nd order recurrence, a(n) = 4 a(n-1) + 3 a(n-2)."


class A015531(RecurrenceSequence2):
    def __init__(self):
        r"""
        Linear 2nd order recurrence, `a(0)=0`, `a(1)=1` and
        `a(n) = 4 a(n-1) + 5 a(n-2)`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A015531;a
            Linear 2nd order recurrence, a(n) = 4 a(n-1) + 5 a(n-2).
            sage: a(0)
            0
            sage: a(1)
            1
            sage: a(2)
            4
            sage: a.offset
            0
            sage: a(8)
            65104
            sage: a(60)
            144560289664733924534327040115992228190104
            sage: a.list(9)
            [0, 1, 4, 21, 104, 521, 2604, 13021, 65104]

        AUTHORS:

        - Jaap Spies (2007-01-19)
        """
        SloaneSequence.__init__(self, offset=0)
        self._params = (0,1,4,5)
        self._b = []
        self._precompute(2)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A015531._repr_()
            'Linear 2nd order recurrence, a(n) = 4 a(n-1) + 5 a(n-2).'
        """
        return "Linear 2nd order recurrence, a(n) = 4 a(n-1) + 5 a(n-2)."


class A015551(RecurrenceSequence2):
    def __init__(self):
        r"""
        Linear 2nd order recurrence, `a(0)=0`, `a(1)=1` and
        `a(n) = 6 a(n-1) + 5 a(n-2)`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A015551;a
            Linear 2nd order recurrence, a(n) = 6 a(n-1) + 5 a(n-2).
            sage: a(0)
            0
            sage: a(1)
            1
            sage: a(2)
            6
            sage: a.offset
            0
            sage: a(8)
            570216
            sage: a(60)
            7110606606530059736761484557155863822531970573036
            sage: a.list(9)
            [0, 1, 6, 41, 276, 1861, 12546, 84581, 570216]

        AUTHORS:

        - Jaap Spies (2007-01-19)
        """
        SloaneSequence.__init__(self, offset=0)
        self._b = []
        self._params = (0,1,6,5)
        self._precompute(2)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A015551._repr_()
            'Linear 2nd order recurrence, a(n) = 6 a(n-1) + 5 a(n-2).'
        """
        return "Linear 2nd order recurrence, a(n) = 6 a(n-1) + 5 a(n-2)."


# todo jsp
#
#
# A015552 = recur_gen2(0,1,6,7)
# A015553 = recur_gen2(0,1,6,11)
# A015555 = recur_gen2(0,1,7,2)
#
# A015565 = recur_gen2(0,1,7,8)
#
# A015585 = recur_gen2(0,1,9,10)
#
# A053404 = recur_gen2(1,1,1,12)
#
# A053428 = recur_gen2(1,1,1,20)
#
# A053430 = recur_gen2(1,1,1,30)
#
# A065874 = recur_gen2(1,1,1,42)
#
# A083858 = recur_gen2(0,1,3,6)
# and more!

# Wilf_A083216 = recur_gen2(20615674205555510, 3794765361567513,1,1)  family


class A082411(RecurrenceSequence2):
    def __init__(self):
        r"""
        Second-order linear recurrence sequence with
        `a(n) = a(n-1) + a(n-2)`.

        `a(0) = 407389224418`, `a(1) = 76343678551`. This
        is the second-order linear recurrence sequence with `a(0)`
        and `a(1)` co-prime, that R. L. Graham in 1964 stated did
        not contain any primes.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A082411;a
            Second-order linear recurrence sequence with a(n) = a(n-1) + a(n-2).
            sage: a(1)
            76343678551
            sage: a(2)
            483732902969
            sage: a(3)
            560076581520
            sage: a(20)
            2219759332689173
            sage: a.list(4)
            [407389224418, 76343678551, 483732902969, 560076581520]

        AUTHORS:

        - Jaap Spies (2007-01-23)
        """
        SloaneSequence.__init__(self, offset=0)
        self._b = []
        self._params = (407389224418,76343678551,1,1)
        self._precompute(2)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A082411._repr_()
            'Second-order linear recurrence sequence with a(n) = a(n-1) + a(n-2).'
        """
        return "Second-order linear recurrence sequence with a(n) = a(n-1) + a(n-2)."


class A083103(RecurrenceSequence2):
    def __init__(self):
        r"""
        Second-order linear recurrence sequence with
        `a(n) = a(n-1) + a(n-2)`.

        `a(0) = 1786772701928802632268715130455793`,
        `a(1) = 1059683225053915111058165141686995`. This is the
        second-order linear recurrence sequence with `a(0)` and
        `a(1)` co- prime, that R. L. Graham in 1964 stated did not
        contain any primes. It has not been verified. Graham made a mistake
        in the calculation that was corrected by D. E. Knuth in 1990.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A083103;a
            Second-order linear recurrence sequence with a(n) = a(n-1) + a(n-2).
            sage: a(1)
            1059683225053915111058165141686995
            sage: a(2)
            2846455926982717743326880272142788
            sage: a(3)
            3906139152036632854385045413829783
            sage: a.offset
            0
            sage: a(8)
            45481392851206651551714764671352204
            sage: a(20)
            14639253684254059531823985143948191708
            sage: a.list(4)
            [1786772701928802632268715130455793, 1059683225053915111058165141686995, 2846455926982717743326880272142788, 3906139152036632854385045413829783]

        AUTHORS:

        - Jaap Spies (2007-01-23)
        """
        SloaneSequence.__init__(self, offset=0)
        self._b = []
        self._params = (1786772701928802632268715130455793,1059683225053915111058165141686995,1,1)
        self._precompute(2)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A083103._repr_()
            'Second-order linear recurrence sequence with a(n) = a(n-1) + a(n-2).'
        """
        return "Second-order linear recurrence sequence with a(n) = a(n-1) + a(n-2)."


class A083104(RecurrenceSequence2):
    def __init__(self):
        r"""
        Second-order linear recurrence sequence with
        `a(n) = a(n-1) + a(n-2)`.

        `a(0) = 331635635998274737472200656430763`,
        `a(1) = 1510028911088401971189590305498785`. This is the
        second-order linear recurrence sequence with `a(0)` and
        `a(1)` co-prime. It was found by Ronald Graham in 1990.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A083104;a
            Second-order linear recurrence sequence with a(n) = a(n-1) + a(n-2).
            sage: a(3)
            3351693458175078679851381267428333
            sage: a.offset
            0
            sage: a(8)
            36021870400834012982120004949074404
            sage: a(20)
            11601914177621826012468849361236300628

        AUTHORS:

        - Jaap Spies (2007-01-23)
        """
        SloaneSequence.__init__(self, offset=0)
        self._b = []
        self._params = (331635635998274737472200656430763,1510028911088401971189590305498785,1,1)
        self._precompute(2)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A083104._repr_()
            'Second-order linear recurrence sequence with a(n) = a(n-1) + a(n-2).'
        """
        return "Second-order linear recurrence sequence with a(n) = a(n-1) + a(n-2)."


class A083105(RecurrenceSequence2):
    def __init__(self):
        r"""
        Second-order linear recurrence sequence with
        `a(n) = a(n-1) + a(n-2)`.

        `a(0) = 62638280004239857`,
        `a(1) = 49463435743205655`. This is the second-order linear
        recurrence sequence with `a(0)` and `a(1)`
        co-prime. It was found by Donald Knuth in 1990.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A083105;a
            Second-order linear recurrence sequence with a(n) = a(n-1) + a(n-2).
            sage: a(1)
            49463435743205655
            sage: a(2)
            112101715747445512
            sage: a(3)
            161565151490651167
            sage: a.offset
            0
            sage: a(8)
            1853029790662436896
            sage: a(20)
            596510791500513098192
            sage: a.list(4)
            [62638280004239857, 49463435743205655, 112101715747445512, 161565151490651167]

        AUTHORS:

        - Jaap Spies (2007-01-23)
        """
        SloaneSequence.__init__(self, offset=0)
        self._b = []
        self._params = (62638280004239857,49463435743205655,1,1)
        self._precompute(2)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A083105._repr_()
            'Second-order linear recurrence sequence with a(n) = a(n-1) + a(n-2).'
        """
        return "Second-order linear recurrence sequence with a(n) = a(n-1) + a(n-2)."


class A083216(RecurrenceSequence2):
    def __init__(self):
        r"""
        Second-order linear recurrence sequence with
        `a(n) = a(n-1) + a(n-2)`.

        `a(0) = 20615674205555510`,
        `a(1) = 3794765361567513`. This is a second-order linear
        recurrence sequence with `a(0)` and `a(1)` co-prime
        that does not contain any primes. It was found by Herbert Wilf in
        1990.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A083216; a
            Second-order linear recurrence sequence with a(n) = a(n-1) + a(n-2).
            sage: a(0)
            20615674205555510
            sage: a(1)
            3794765361567513
            sage: a(8)
            347693837265139403
            sage: a(41)
            2738025383211084205003383
            sage: a.list(4)
            [20615674205555510, 3794765361567513, 24410439567123023, 28205204928690536]

        AUTHORS:

        - Jaap Spies (2007-01-19)
        """
        SloaneSequence.__init__(self, offset=0)
        self._b = []
        self._params = (20615674205555510, 3794765361567513,1,1)
        self._precompute(2)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A083216._repr_()
            'Second-order linear recurrence sequence with a(n) = a(n-1) + a(n-2).'
        """
        return "Second-order linear recurrence sequence with a(n) = a(n-1) + a(n-2)."


class A061084(SloaneSequence):
    def __init__(self):
        r"""
        Fibonacci-type sequence based on subtraction: `a(0) = 1`,
        `a(1) = 2` and `a(n) = a(n-2)-a(n-1)`.

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A061084; a
            Fibonacci-type sequence based on subtraction: a(0) = 1, a(1) = 2 and a(n) = a(n-2)-a(n-1).
            sage: a(0)
            1
            sage: a(1)
            2
            sage: a(8)
            -29
            sage: a(22)
            -24476
            sage: a.list(12)
            [1, 2, -1, 3, -4, 7, -11, 18, -29, 47, -76, 123]
            sage: a.keyword
            ['sign', 'easy', 'nice']

        AUTHORS:

        - Jaap Spies (2007-01-18)
        """
        SloaneSequence.__init__(self, offset=0)

    keyword = ["sign", "easy","nice"]

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A061084._repr_()
            'Fibonacci-type sequence based on subtraction: a(0) = 1, a(1) = 2 and a(n) = a(n-2)-a(n-1).'
        """
        return "Fibonacci-type sequence based on subtraction: a(0) = 1, a(1) = 2 and a(n) = a(n-2)-a(n-1)."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A061084._eval(n) for n in range(10)]
            [1, 2, -1, 3, -4, 7, -11, 18, -29, 47]
        """
        if n == 0:
            return ZZ.one()
        elif n == 1:
            return 2
        else:
            return (-1)**(n-1)*sloane.A000204(n-1)


# a group of sequences uses this function:
def recur_gen3(a0, a1, a2, a3, a4, a5):
    r"""
    Homogeneous general third-order linear recurrence generator with
    fixed coefficients

    a(0) = a0, a(1) = a1, a(2) = a2, a(n) = a3\*a(n-1) + a4\*a(n-2) +
    a5\*a(n-3)

    EXAMPLES::

        sage: from sage.combinat.sloane_functions import recur_gen3
        sage: it = recur_gen3(1,1,1,1,1,1)
        sage: [next(it) for i in range(10)]
        [1, 1, 1, 3, 5, 9, 17, 31, 57, 105]
    """
    x, y, z = ZZ(a0), ZZ(a1), ZZ(a2)
    yield x
    while True:
        x, y, z = y, z, a5*x+a4*y+a3*z
        yield x


class A000213(SloaneSequence):
    def __init__(self):
        r"""
        Tribonacci numbers: a(n) = a(n-1) + a(n-2) + a(n-3). Starting with
        1, 1, 1, ...

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000213;a
            Tribonacci numbers: a(n) = a(n-1) + a(n-2) + a(n-3).
            sage: a(0)
            1
            sage: a(1)
            1
            sage: a(2)
            1
            sage: a(11)
            355
            sage: a.list(12)
            [1, 1, 1, 3, 5, 9, 17, 31, 57, 105, 193, 355]

        AUTHORS:

        - Jaap Spies (2007-01-19)
        """
        SloaneSequence.__init__(self, offset=0)
        self._b = []
        self._precompute()

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000213._repr_()
            'Tribonacci numbers: a(n) = a(n-1) + a(n-2) + a(n-3).'
        """
        return "Tribonacci numbers: a(n) = a(n-1) + a(n-2) + a(n-3)."

    def _precompute(self, how_many=20):
        """
        EXAMPLES::

            sage: initial = len(sloane.A000213._b)
            sage: sloane.A000213._precompute(10)
            sage: len(sloane.A000213._b) - initial == 10
            True
        """
        try:
            f = self._f
        except AttributeError:
            self._f = recur_gen3(1,1,1,1,1,1)
            f = self._f
        self._b += [next(f) for i in range(how_many)]

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000213._eval(n) for n in range(10)]
            [1, 1, 1, 3, 5, 9, 17, 31, 57, 105]
        """
        if len(self._b) <= n:
            self._precompute(n - len(self._b) + 1)
        return self._b[n]

    def list(self, n):
        """
        EXAMPLES::

            sage: sloane.A000213.list(10)
            [1, 1, 1, 3, 5, 9, 17, 31, 57, 105]
        """
        self._eval(n)   # force computation
        return self._b[:n]


class A000073(SloaneSequence):
    def __init__(self):
        r"""
        Tribonacci numbers: a(n) = a(n-1) + a(n-2) + a(n-3). Starting with
        0, 0, 1, ...

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000073;a
            Tribonacci numbers: a(n) = a(n-1) + a(n-2) + a(n-3).
            sage: a(0)
            0
            sage: a(1)
            0
            sage: a(2)
            1
            sage: a(11)
            149
            sage: a.list(12)
            [0, 0, 1, 1, 2, 4, 7, 13, 24, 44, 81, 149]

        AUTHORS:

        - Jaap Spies (2007-01-19)
        """
        SloaneSequence.__init__(self, offset=0)
        self._b = []
        self._precompute()

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000073._repr_()
            'Tribonacci numbers: a(n) = a(n-1) + a(n-2) + a(n-3).'
        """
        return "Tribonacci numbers: a(n) = a(n-1) + a(n-2) + a(n-3)."

    def _precompute(self, how_many=20):
        """
        EXAMPLES::

            sage: initial = len(sloane.A000073._b)
            sage: sloane.A000073._precompute(10)
            sage: len(sloane.A000073._b) - initial == 10
            True
        """
        try:
            f = self._f
        except AttributeError:
            self._f = recur_gen3(0,0,1,1,1,1)
            f = self._f
        self._b += [next(f) for i in range(how_many)]

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000073._eval(n) for n in range(10)]
            [0, 0, 1, 1, 2, 4, 7, 13, 24, 44]
        """
        if len(self._b) <= n:
            self._precompute(n - len(self._b) + 1)
        return self._b[n]

    def list(self, n):
        """
        EXAMPLES::

            sage: sloane.A000073.list(10)
            [0, 0, 1, 1, 2, 4, 7, 13, 24, 44]
        """
        self._eval(n)   # force computation
        return self._b[:n]


def perm_mh(m, h):
    r"""
    This functions calculates `f(g,h)` from Sloane's sequences
    A079908-A079928

    INPUT:

    - ``m`` -- positive integer
    - ``h`` -- nonnegative integer

    OUTPUT: permanent of the `m \times (m+h)` matrix, etc.

    EXAMPLES::

        sage: from sage.combinat.sloane_functions import perm_mh
        sage: perm_mh(3,3)
        36
        sage: perm_mh(3,4)
        76

    AUTHORS:

    - Jaap Spies (2006)
    """
    n = m + h
    M = MatrixSpace(ZZ, m, n)
    A = M(0)
    for i in range(m):
        for j in range(n):
            if i <= j and j <= i + h:
                A[i, j] = 1
    return A.permanent()


class A079922(SloaneSequence):
    r"""
    Function returns solutions to the Dancing School problem with
    `n` girls and `n+3` boys.

    The value is `per(B)`, the permanent of the (0,1)-matrix
    `B` of size `n \times n+3` with `b(i,j)=1`
    if and only if `i \le j \le i+n`.

    REFERENCES:

    - Jaap Spies, Nieuw Archief voor Wiskunde, 5/7 nr 4, December 2006

    INPUT:

    - ``n`` -- positive integer

    OUTPUT: integer; function value

    EXAMPLES::

        sage: a = sloane.A079922; a
        Solutions to the Dancing School problem with n girls and n+3 boys
        sage: a.offset
        1
        sage: a(1)
        4
        sage: a(8)
        2227
        sage: a.list(8)
        [4, 13, 36, 90, 212, 478, 1044, 2227]

    Compare: Searching Sloane's online database... Solution to the
    Dancing School Problem with n girls and n+3 boys: f(n,3). [4, 13,
    36, 90, 212, 478, 1044, 2227]

    ::

        sage: a(-1)
        Traceback (most recent call last):
        ...
        ValueError: input n (=-1) must be a positive integer

    AUTHORS:

        - Jaap Spies (2007-01-14)
    """

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A079922._repr_()
            'Solutions to the Dancing School problem with n girls and n+3 boys'
        """
        return "Solutions to the Dancing School problem with n girls and n+3 boys"

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A079922._eval(n) for n in range(1,5)]
            [4, 13, 36, 90]
        """
        return perm_mh(n, 3)


class A079923(SloaneSequence):
    r"""
    Function returns solutions to the Dancing School problem with
    `n` girls and `n+4` boys.

    The value is `per(B)`, the permanent of the (0,1)-matrix
    `B` of size `n \times n+3` with `b(i,j)=1`
    if and only if `i \le j \le i+n`.

    REFERENCES:

    - Jaap Spies, Nieuw Archief voor Wiskunde, 5/7 nr 4,
      December 2006

    INPUT:

    - ``n`` -- positive integer

    OUTPUT: integer; function value

    EXAMPLES::

        sage: a = sloane.A079923; a
        Solutions to the Dancing School problem with n girls and n+4 boys
        sage: a.offset
        1
        sage: a(1)
        5
        sage: a(8)
        15458
        sage: a.list(8)
        [5, 21, 76, 246, 738, 2108, 5794, 15458]

    Compare: Searching Sloane's online database... Solution to the
    Dancing School Problem with n girls and n+4 boys: f(n,4). [5, 21,
    76, 246, 738, 2108, 5794, 15458]

    ::

        sage: a(0)
        Traceback (most recent call last):
        ...
        ValueError: input n (=0) must be a positive integer

    AUTHORS:

        - Jaap Spies (2007-01-17)
    """

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A079923._repr_()
            'Solutions to the Dancing School problem with n girls and n+4 boys'
        """
        return "Solutions to the Dancing School problem with n girls and n+4 boys"

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A079923._eval(n) for n in range(1,11)]
            [5, 21, 76, 246, 738, 2108, 5794, 15458, 40296, 103129]
        """
        return perm_mh(n, 4)


class A109814(SloaneSequence):
    r"""
    The `n` th term of the sequence `a(n)` is the
    largest `k` such that `n` can be written as sum of
    `k` consecutive integers.

    By definition, `n` is the sum of at most `a(n)` consecutive
    positive integers. Suppose `n` is to be written as sum of `k`
    consecutive integers starting with `m`, then `2n = k(2m + k -
    1)`. Only one of the factors is odd. For each odd divisor `d`
    of `n` there is a unique corresponding `k =
    min(d,2n/d)`. `a(n)` can be alternatively defined as the
    largest among those `k` .

    .. SEEALSO::

        * :wikipedia:`Polite_number`

        * `An exercise sheet (with answers) about sums of
          consecutive integers
          <http://www.jaapspies.nl/mathfiles/problem2005-2C.pdf>`_.

    INPUT:

    - ``n`` -- nonnegative integer

    OUTPUT: integer; function value

    EXAMPLES::

        sage: a = sloane.A109814; a
        a(n) is the largest k such that n can be written as sum of k consecutive positive integers.
        sage: a(0)
        Traceback (most recent call last):
        ...
        ValueError: input n (=0) must be a positive integer
        sage: a(2)
        1
        sage: a.list(9)
        [1, 1, 2, 1, 2, 3, 2, 1, 3]

    AUTHORS:

    - Jaap Spies (2007-01-13)
    """

    def __init__(self):
        r"""
        EXAMPLES::

            sage: a = sloane.A109814; a
            a(n) is the largest k such that n can be written as sum of k consecutive positive integers.
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(2)
            1
            sage: a.list(9)
            [1, 1, 2, 1, 2, 3, 2, 1, 3]
        """
        SloaneSequence.__init__(self, offset=1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A109814._repr_()
            'a(n) is the largest k such that n can be written as sum of k consecutive positive integers.'
        """
        return "a(n) is the largest k such that n can be written as sum of k consecutive positive integers."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A109814._eval(n) for n in range(1, 10)]
            [1, 1, 2, 1, 2, 3, 2, 1, 3]
        """
        if n == 1:
            return ZZ.one()
        m = 0
        for d in arith.divisors(n):
            if not d % 2:
                continue
            # d is odd divisor
            k = min(d, 2 * n // d)
            m = max(k, m)
        return ZZ(m)


class A111774(SloaneSequence):
    def __init__(self):
        r"""
        Sequence of numbers of the third kind, i.e., numbers that can be
        written as a sum of at least three consecutive positive integers.

        Odd primes can only be written as a sum of two consecutive
        integers. Powers of 2 do not have a representation as a sum of
        `k` consecutive integers (other than the trivial
        `n = n` for `k = 1`).

        See: http://www.jaapspies.nl/mathfiles/problem2005-2C.pdf

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A111774; a
            Numbers that can be written as a sum of at least three consecutive positive integers.
            sage: a(1)
            6
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(100)
            141
            sage: a(156)
            209
            sage: a(302)
            386
            sage: a.list(12)
            [6, 9, 10, 12, 14, 15, 18, 20, 21, 22, 24, 25]
            sage: a(1/3)
            Traceback (most recent call last):
            ...
            TypeError: input must be an int or Integer

        AUTHORS:

        - Jaap Spies (2007-01-13)
        """
        SloaneSequence.__init__(self, offset=1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A111774._repr_()
            'Numbers that can be written as a sum of at least three consecutive positive integers.'
        """
        return "Numbers that can be written as a sum of at least three consecutive positive integers."

    def _precompute(self, how_many=150):
        """
        EXAMPLES::

            sage: initial = len(sloane.A111774._b)
            sage: sloane.A111774._precompute()
            sage: len(sloane.A111774._b) - initial > 0
            True
        """
        try:
            self._b
            self._n
        except AttributeError:
            self._b = []
            self._n = 1
        self._b += [i for i in range(self._n, self._n + how_many)
                    if self.is_number_of_the_third_kind(i)]
        self._n += how_many

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A111774._eval(n) for n in range(1,11)]
            [6, 9, 10, 12, 14, 15, 18, 20, 21, 22]
        """
        try:
            return self._b[n-1]
        except (AttributeError, IndexError):
            self._precompute()
            # try again
            return self._eval(n)

    def list(self, n):
        """
        EXAMPLES::

            sage: sloane.A111774.list(12)
            [6, 9, 10, 12, 14, 15, 18, 20, 21, 22, 24, 25]
        """
        try:
            if len(self._b) < n:
                raise IndexError
            else:
                return self._b[:n]
        except (AttributeError, IndexError):
            self._precompute()
            # try again
            return self.list(n)

    def is_number_of_the_third_kind(self, n) -> bool:
        r"""
        Return ``True`` if and only if `n` is a number of the third kind.

        A number is of the third kind if it can be written as a sum of at
        least three consecutive positive integers. Odd primes can only be
        written as a sum of two consecutive integers. Powers of 2 do not
        have a representation as a sum of `k` consecutive integers
        (other than the trivial `n = n` for `k = 1`).

        See: http://www.jaapspies.nl/mathfiles/problem2005-2C.pdf

        INPUT:

        - ``n`` -- positive integer

        OUTPUT: ``True`` if `n` is not prime and not a power of 2

        EXAMPLES::

            sage: a = sloane.A111774
            sage: a.is_number_of_the_third_kind(6)
            True
            sage: a.is_number_of_the_third_kind(100)
            True
            sage: a.is_number_of_the_third_kind(16)
            False
            sage: a.is_number_of_the_third_kind(97)
            False

        AUTHORS:

        - Jaap Spies (2006-12-09)
        """
        return not arith.is_prime(n) and not arith.is_power_of_two(n)


class A111775(SloaneSequence):
    def __init__(self):
        r"""
        Number of ways `n` can be written as a sum of at least
        three consecutive integers.

        Powers of 2 and (odd) primes can not be written as a sum of at
        least three consecutive integers. `a(n)` strongly depends
        on the number of odd divisors of `n` (A001227): Suppose
        `n` is to be written as sum of `k` consecutive
        integers starting with `m`, then
        `2n = k(2m + k - 1)`. Only one of the factors is odd. For
        each odd divisor of `n` there is a unique corresponding
        `k`, `k=1` and `k=2` must be excluded.

        See: http://www.jaapspies.nl/mathfiles/problem2005-2C.pdf

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A111775; a
            Number of ways n can be written as a sum of at least three consecutive integers.

        ::

            sage: a(1)
            0
            sage: a(0)
            0

        We have a(15)=2 because 15 = 4+5+6 and 15 = 1+2+3+4+5. The number
        of odd divisors of 15 is 4.

        ::

            sage: a(15)
            2

        ::

            sage: a(100)
            2
            sage: a(256)
            0
            sage: a(29)
            0
            sage: a.list(20)
            [0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 2, 0, 0, 2, 0]
            sage: a(1/3)
            Traceback (most recent call last):
            ...
            TypeError: input must be an int or Integer

        AUTHORS:

        - Jaap Spies (2006-12-09)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A111775._repr_()
            'Number of ways n can be written as a sum of at least three consecutive integers.'
        """
        return "Number of ways n can be written as a sum of at least three consecutive integers."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A111775._eval(n) for n in range(10)]
            [0, 0, 0, 0, 0, 0, 1, 0, 0, 1]
        """
        if n == 1 or n == 0:
            return ZZ.zero()
        k = sum(i % 2 for i in arith.divisors(n))  # A001227, the number of odd divisors
        if n % 2 == 0:
            return k - 1
        else:
            return k - 2


class A111787(SloaneSequence):
    def __init__(self):
        r"""
        This function returns the `n`-th number of Sloane's
        sequence A111787

        `a(n)=0` if `n` is an odd prime or a power of 2.
        For numbers of the third kind (see A111774) we proceed as follows:
        suppose `n` is to be written as sum of `k`
        consecutive integers starting with `m`, then
        `2n = k(2m + k - 1)`. Let `p` be the smallest odd
        prime divisor of `n` then `a(n) = min(p,2n/p)`.

        See: http://www.jaapspies.nl/mathfiles/problem2005-2C.pdf

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A111787; a
            a(n) is the least k >= 3 such that n can be written as sum of k consecutive integers. a(n)=0 if such a k does not exist.
            sage: a.offset
            1
            sage: a(1)
            0
            sage: a(0)
            Traceback (most recent call last):
            ...
            ValueError: input n (=0) must be a positive integer
            sage: a(100)
            5
            sage: a(256)
            0
            sage: a(29)
            0
            sage: a.list(20)
            [0, 0, 0, 0, 0, 3, 0, 0, 3, 4, 0, 3, 0, 4, 3, 0, 0, 3, 0, 5]
            sage: a(-1)
            Traceback (most recent call last):
            ...
            ValueError: input n (=-1) must be a positive integer

        AUTHORS:

        - Jaap Spies (2007-01-14)
        """
        SloaneSequence.__init__(self, offset=1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A111787._repr_()
            'a(n) is the least k >= 3 such that n can be written as sum of k consecutive integers. a(n)=0 if such a k does not exist.'
        """
        return "a(n) is the least k >= 3 such that n can be written as sum of k consecutive integers. a(n)=0 if such a k does not exist."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A111787._eval(n) for n in range(1,11)]
            [0, 0, 0, 0, 0, 3, 0, 0, 3, 4]
        """
        if arith.is_prime(n) or arith.is_power_of_two(n):
            return ZZ.zero()

        for d in srange(3, n, 2):
            if n % d == 0:
                return min(d, 2*n//d)


class ExponentialNumbers(SloaneSequence):
    def __init__(self, a):
        r"""
        A sequence of Exponential numbers.

        EXAMPLES::

            sage: from sage.combinat.sloane_functions import ExponentialNumbers
            sage: ExponentialNumbers(0)
            Sequence of Exponential numbers around 0
        """
        SloaneSequence.__init__(self, offset=0)
        self.a = a

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.combinat.sloane_functions import ExponentialNumbers
            sage: ExponentialNumbers(4)._repr_()
            'Sequence of Exponential numbers around 4'
        """
        return "Sequence of Exponential numbers around %s" % self.a

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000110._eval(n) for n in range(10)]
            [1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147]
        """
        if hasattr(self, '__n'):
            if n < self.__n:
                return self.__data[n]
        from sage.combinat.expnums import expnums
        self.__data = expnums(n + 1, self.a)
        self.__n = n + 1
        return self.__data[n]


class A000110(ExponentialNumbers):
    def __init__(self):
        r"""
        The sequence of Bell numbers.

        The Bell number `B_n` counts the number of ways to put
        `n` distinguishable things into indistinguishable boxes
        such that no box is empty.

        Let `S(n, k)` denote the Stirling number of the second
        kind. Then

        .. MATH::

            B_n = \sum{k=0}^{n} S(n, k).

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; `B_n`

        EXAMPLES::

            sage: a = sloane.A000110; a
            Sequence of Bell numbers
            sage: a.offset
            0
            sage: a(0)
            1
            sage: a(100)
            47585391276764833658790768841387207826363669686825611466616334637559114497892442622672724044217756306953557882560751
            sage: a.list(10)
            [1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147]

        AUTHORS:

        - Nick Alexander
        """
        ExponentialNumbers.__init__(self, a=1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000110._repr_()
            'Sequence of Bell numbers'
        """
        return "Sequence of Bell numbers"


class A000587(ExponentialNumbers):
    def __init__(self):
        r"""
        The sequence of Uppuluri-Carpenter numbers.

        The Uppuluri-Carpenter number `C_n` counts the imbalance
        in the number of ways to put `n` distinguishable things
        into an even number of indistinguishable boxes versus into an odd
        number of indistinguishable boxes, such that no box is empty.

        Let `S(n, k)` denote the Stirling number of the second
        kind. Then

        .. MATH::

            C_n = \sum{k=0}^{n} (-1)^k S(n, k) .

        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; `C_n`

        EXAMPLES::

            sage: a = sloane.A000587; a
            Sequence of Uppuluri-Carpenter numbers
            sage: a.offset
            0
            sage: a(0)
            1
            sage: a(100)
            397577026456518507969762382254187048845620355238545130875069912944235105204434466095862371032124545552161
            sage: a.list(10)
            [1, -1, 0, 1, 1, -2, -9, -9, 50, 267]

        AUTHORS:

        - Nick Alexander
        """
        ExponentialNumbers.__init__(self, a=-1)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000587._repr_()
            'Sequence of Uppuluri-Carpenter numbers'
        """
        return "Sequence of Uppuluri-Carpenter numbers"


# A000100  a(n) = number of compositions of n in which the maximum part size is 3. Milestone!
#  a(n+3) = Sum[k=0..n, F(k)*T(n-k) ], F(i)=A000045(i+1), T(i)=A000073(i+2).
#  0, 0, 0, 1, 2, 5, 11, 23, 47, 94, 185, 360, 694, 1328, 2526, 4781, 9012, 16929, 31709, 59247


class A000100(SloaneSequence):
    def __init__(self):
        r"""
        INPUT:

        - ``n`` -- nonnegative integer

        OUTPUT: integer; function value

        EXAMPLES::

            sage: a = sloane.A000100;a
            Number of compositions of n in which the maximum part size is 3.
            sage: a(0)
            0
            sage: a(1)
            0
            sage: a(2)
            0
            sage: a(3)
            1
            sage: a(11)
            360
            sage: a.list(12)
            [0, 0, 0, 1, 2, 5, 11, 23, 47, 94, 185, 360]

        AUTHORS:

        - Jaap Spies (2007-01-26)
        """
        SloaneSequence.__init__(self, offset=0)

    def _repr_(self):
        """
        EXAMPLES::

            sage: sloane.A000100._repr_()
            'Number of compositions of n in which the maximum part size is 3.'
        """
        return "Number of compositions of n in which the maximum part size is 3."

    def _eval(self, n):
        """
        EXAMPLES::

            sage: [sloane.A000100._eval(n) for n in range(10)]
            [0, 0, 0, 1, 2, 5, 11, 23, 47, 94]
        """
        if n <= 2:
            return ZZ.zero()
        return sum(sloane.A000045(i + 1) * sloane.A000073(n - i - 1)
                   for i in range(n - 2))


#############################################################
# III. Create the Sloane object, off which all the sequence
#      objects are members.
#############################################################

class Sloane(SageObject):
    r"""
    A collection of Sloane generating functions.

    This class inspects sage.combinat.sloane_functions, accumulating
    all the SloaneSequence classes starting with 'A'. These are listed
    for tab completion, but not instantiated until requested.

    EXAMPLES:

    Ensure we have lots of entries::

        sage: len(sloane.__dir__()) > 100
        True

    Ensure none are being incorrectly returned::

        sage: [ None for n in sloane.__dir__() if not n.startswith('A') ]
        []

    Ensure we can access dynamic constructions and cache correctly::

        sage: s = sloane.A000587
        sage: s is sloane.A000587
        True

    Ensure that we can access other functions in parent classes::

        sage: sloane.__class__
        <class 'sage.combinat.sloane_functions.Sloane'>

    AUTHORS:

    - Nick Alexander
    """

    def __dir__(self):
        r"""
        List Sloane generating functions for tab-completion.

        The member classes are inspected from module
        :mod:`sage.combinat.sloane_functions`.

        They must be sub classes of SloaneSequence and must start with 'A'.
        These restrictions are only to prevent typos, incorrect inspecting,
        etc.

        EXAMPLES::

            sage: type(sloane.__dir__())
            <class 'list'>
        """
        try:
            return self.__stored_dir
        except AttributeError:
            xs = inspect.getmembers(sys.modules[__name__], inspect.isclass)
            self.__stored_dir = [n for n, c in xs
                                 if n.startswith('A') and issubclass(c, SloaneSequence)]
            return self.__stored_dir

    def __getattribute__(self, name):
        r"""
        Construct and cache unique instances of Sloane generating
        function objects.

        EXAMPLES::

            sage: sloane.__getattribute__('A000001')
            Number of groups of order n.
            sage: sloane.__getattribute__('dog')
            Traceback (most recent call last):
            ...
            AttributeError: dog

        ::

            sage: sloane.__repr__
            <built-in method __repr__ of Sloane object at 0x...>
            sage: sloane.__name__
            Traceback (most recent call last):
            ...
            AttributeError: __name__
        """
        try:
            return SageObject.__getattribute__(self, name)
        except AttributeError:
            try:
                f = getattr(sys.modules[__name__], name)
                seq = f()
                setattr(self, name, seq)
                return seq
            except (AttributeError, TypeError):
                raise AttributeError(name)


sloane = Sloane()
