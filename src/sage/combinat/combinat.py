r"""
Combinatorial functions

This module implements some combinatorial functions, as listed
below. For a more detailed description, see the relevant
docstrings.

**Numbers:**

-  Bell numbers, :func:`bell_number`

-  Catalan numbers, :func:`catalan_number`

-  Narayana numbers, :func:`narayana_number`

-  Euler numbers, :func:`euler_number`

-  Eulerian numbers, :func:`eulerian_number`

-  Fibonacci numbers, :func:`fibonacci`

-  Lucas numbers, :func:`lucas_number1`, :func:`lucas_number2`.

-  Stirling numbers, :func:`stirling_number1`,
   :func:`stirling_number2`.

-  Polygonal numbers, :func:`polygonal_number`

**Polynomials**

-  Eulerian polynomial, :func:`eulerian_polynomial`

-  Bernoulli polynomials, :func:`bernoulli_polynomial`

**Sets:**

-  Tuples of a multiset, :func:`tuples` and :func:`number_of_tuples`. An
   ordered tuple of length `k` of a set `S` is a ordered selection with
   repetitions of `S` and is represented by a sorted list of length `k`
   containing elements from `S`.

-  Unordered tuples of a set, :func:`unordered_tuples` and
   :func:`number_of_unordered_tuples`. An unordered tuple of length `k` of a
   set `S` is an unordered selection with repetitions of `S` and is represented
   by a sorted list of length `k` containing elements from `S`.

**Combinatorial functions from other modules:**

-  :func:`sage.arith.misc.binomial` binomial coefficient

-  :func:`sage.arith.misc.factorial` factorial

-  :func:`sage.arith.misc.falling_factorial` falling power

-  :func:`sage.arith.misc.rising_factorial` rising power

-  :func:`sage.combinat.partition.number_of_partitions` number of partitions

-  :func:`sage.combinat.q_analogues.gaussian_binomial` Gaussian binomial coefficient

.. TODO::

    Add GUAVA commands:
        * VandermondeMat
        * GrayMat returns a list of all different vectors of length n over
          the field F, using Gray ordering.
    Add (not in GAP):
        * Rencontres numbers (:wikipedia:`Rencontres_number`)

REFERENCES:

- :wikipedia:`Twelvefold_way`

AUTHORS:

- David Joyner (2006-07): initial implementation
- William Stein (2006-07): editing of docs and code; many
  optimizations, refinements, and bug fixes in corner cases
- Jon Hanke (2006-08): added ``tuples``
- David Joyner (2006-09): bug fix for combinations, added
  permutations_iterator, combinations_iterator from Python Cookbook,
  edited docs.
- David Joyner (2007-11): changed permutations, added hadamard_matrix
- Blair Sutton (2009-01-26): added ``bell_polynomial``
- Florent Hivert (2009-02): combinatorial class cleanup
- Bobby Moretti (20009-02): added ``fibonacci_sequence``
- Bobby Moretti (2009-02): added ``fibonacci_xrange``
- Fredrik Johansson (2010-07): fast implementation of ``stirling_number2``
- Robert Gerbicz (2010-10): added Bell numbers
- Punarbasu Purkayastha (2012-12): deprecate arrangements, combinations,
  combinations_iterator, and clean up very old deprecated methods.
- Jeroen Demeyer (2014-10): improved implementation of Dobinski formula for Bell numbers
  with more accurate error estimates (:issue:`17157`)
- Thierry Monteil (2015-09-29): the result of ``bell_polynomial`` must
  always be a polynomial
- Kei Beauduin (2024-04-06): when univariate, the Bell polynomial is in
  variable ``x0``; extended to complete exponential, partial ordinary and
  complete ordinary Bell polynomials.
- Kwankyu Lee (2025-01): added Lah numbers
"""

# ****************************************************************************
#       Copyright (C) 2006      David Joyner <wdjoyner@gmail.com>
#                     2007-2009 Mike Hansen <mhansen@gmail.com>
#                     2006      William Stein <wstein@gmail.com>
#                     2009      Blair Sutton
#                     2009      Craig Citro
#                     2009-2010 Florent Hivert
#                     2010      Francis Clarke
#                     2010      Fredrik Johansson
#                     2010      Robert Gerbicz
#                     2010-2013 Nathann Cohen
#                     2012      Christian Stump
#                     2013-2015 Travis Scrimshaw
#                     2014      Volker Braun
#                     2014-2015 Darij Grinberg
#                     2014-2015 Jeroen Demeyer
#                     2014-2021 Frédéric Chapoton
#                     2015      Thierry Monteil
#                     2019      Christian Nassau
#                     2019-2020 Alex Shearer
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from __future__ import annotations
from collections.abc import Iterator

from sage.arith.misc import bernoulli, factorial
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.integer import Integer
from sage.rings.infinity import infinity
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.misc.misc_c import prod
from sage.misc.cachefunc import cached_function
from sage.structure.sage_object import SageObject
from sage.structure.parent import Parent
from sage.misc.lazy_import import lazy_import
from sage.misc.lazy_attribute import lazy_attribute
from .combinat_cython import _stirling_number2
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.structure.element import Element

lazy_import('sage.interfaces.maxima_lib', 'maxima')
lazy_import('sage.libs.pari.all', 'pari')
lazy_import('sage.misc.prandom', 'randint')


def bell_number(n, algorithm='flint', **options) -> Integer:
    r"""
    Return the `n`-th Bell number.

    This is the number of ways to partition a set
    of `n` elements into pairwise disjoint nonempty subsets.

    INPUT:

    - ``n`` -- positive integer

    - ``algorithm`` -- (default: ``'flint'``) any one of the following:

      - ``'dobinski'`` -- use Dobinski's formula

      - ``'flint'`` -- wrap FLINT's ``arith_bell_number``

      - ``'gap'`` -- wrap GAP's ``Bell``

      - ``'mpmath'`` -- wrap mpmath's ``bell``

    .. WARNING::

        When using the mpmath algorithm to compute Bell numbers and you specify
        ``prec``, it can return incorrect results due to low precision. See
        the examples section.


    EXAMPLES::

        sage: # needs sage.libs.flint
        sage: bell_number(10)
        115975
        sage: bell_number(2)
        2
        sage: bell_number(-10)
        Traceback (most recent call last):
        ...
        ArithmeticError: Bell numbers not defined for negative indices
        sage: bell_number(1)
        1
        sage: bell_number(1/3)
        Traceback (most recent call last):
        ...
        TypeError: no conversion of this rational to integer

    When using the mpmath algorithm, we are required have mpmath's precision
    set to at least `\log_2(B_n)` bits. If upon computing the Bell number the
    first time, we deem the precision too low, we use our guess to
    (temporarily) raise mpmath's precision and the Bell number is recomputed. ::

        sage: k = bell_number(30, 'mpmath'); k                                          # needs mpmath
        846749014511809332450147
        sage: k == bell_number(30)                                                      # needs mpmath sage.libs.flint
        True

    If you knows what precision is necessary before computing the Bell number,
    you can use the ``prec`` option::

        sage: k2 = bell_number(30, 'mpmath', prec=30); k2                               # needs mpmath
        846749014511809332450147
        sage: k == k2                                                                   # needs mpmath
        True

    .. WARNING::

            Running mpmath with the precision set too low can result in
            incorrect results::

                sage: k = bell_number(30, 'mpmath', prec=15); k                         # needs mpmath
                846749014511809388871680
                sage: k == bell_number(30)                                              # needs mpmath sage.libs.flint
                False

    TESTS::

        sage: all(bell_number(n) == bell_number(n,'dobinski') for n in range(100))      # needs sage.libs.flint
        True
        sage: all(bell_number(n) == bell_number(n,'gap') for n in range(100))           # needs sage.libs.flint sage.libs.gap
        True
        sage: all(bell_number(n) == bell_number(n,'mpmath', prec=500)                   # needs mpmath sage.libs.flint
        ....:     for n in range(200, 220))
        True

    .. NOTE::

        Let `B_n` denote the `n`-th Bell number. Dobinski's formula is:

        .. MATH::

            B_n = e^{-1} \sum_{k=0}^{\infty} \frac{k^n}{k!}.

        To show our implementation of Dobinski's method works, suppose that `n \geq 5`
        and let `k_0` be the smallest positive integer such that `\frac{k_0^n}{k_0!} < 1`.
        Note that `k_0 > n` and `k_0 \leq 2n` because we can prove that
        `\frac{(2n)^n}{(2n)!} < 1` by Stirling.

        If `k > k_0`, then we have `\frac{k^n}{k!} < \frac{1}{2^{k-k_0}}`.
        We show this by induction:
        let `c_k = \frac{k^n}{k!}`, if `k > n` then

        .. MATH::

            \frac{c_{k+1}}{c_k} = \frac{(1+k^{-1})^n}{k+1} < \frac{(1+n^{-1})^n}{n}
            < \frac{1}{2}.

        The last inequality can easily be checked numerically for `n \geq 5`.

        Using this, we can see that `\frac{c_k}{c_{k_0}} < \frac{1}{2^{k-k_0}}`
        for `k > k_0 > n`. So summing this it gives that `\sum_{k=k_0+1}^{\infty}
        \frac{k^n}{k!} < 1`, and hence

        .. MATH::

            B_n = e^{-1} \left( \sum_{k=0}^{k_0} \frac{k^n}{k!} + E_1 \right)
            = e^{-1} \sum_{k=0}^{k_0} \frac{k^n}{k!} + E_2,

        where `0 < E_1 < 1` and `0 < E_2 < e^{-1}`. Next we have for any `q > 0`

        .. MATH::

            \sum_{k=0}^{k_0} \frac{k^n}{k!} = \frac{1}{q} \sum_{k=0}^{k_0} \left\lfloor
            \frac{q k^n}{k!} \right\rfloor + \frac{E_3}{q}

        where `0 \leq E_3 \leq k_0 + 1 \leq 2n + 1`. Let `E_4 = \frac{E_3}{q}`
        and let `q = 2n + 1`. We find `0 \leq E_4 \leq 1`. These two bounds give:

        .. MATH::

            \begin{aligned}
            B_n & = \frac{e^{-1}}{q} \sum_{k=0}^{k_0} \left\lfloor
            \frac{q k^n}{k!} \right\rfloor + e^{-1} E_4 + E_2 \\
            & = \frac{e^{-1}}{q} \sum_{k=0}^{k_0} \left\lfloor \frac{q k^n}{k!}
            \right\rfloor + E_5
            \end{aligned}

        where

        .. MATH::

            0 < E_5 = e^{-1} E_4 + E_2 \leq e^{-1} + e^{-1} < \frac{3}{4}.

        It follows that

        .. MATH::

            B_n = \left\lceil \frac{e^{-1}}{q} \sum_{k=0}^{k_0} \left\lfloor
            \frac{q k^n}{k!} \right\rfloor \right\rceil.

        Now define

        .. MATH::

            b = \sum_{k=0}^{k_0} \left\lfloor \frac{q k^n}{k!} \right\rfloor.

        This `b` can be computed exactly using integer arithmetic.
        To avoid the costly integer division by `k!`, we collect
        more terms and do only one division, for example with 3 terms:

        .. MATH::

            \frac{k^n}{k!} + \frac{(k+1)^n}{(k+1)!} + \frac{(k+2)^n}{(k+2)!}
            = \frac{k^n (k+1)(k+2) + (k+1)^n (k+2) + (k+2)^n}{(k+2)!}

        In the implementation, we collect `\sqrt{n}/2` terms.

        To actually compute `B_n` from `b`,
        we let `p = \lfloor \log_2(b) \rfloor + 1` such that `b < 2^p` and
        we compute with `p` bits of precision.
        This implies that `b` (and `q < b`) can be represented exactly.

        We compute `\frac{e^{-1}}{q} b`, rounding down, and we must have an
        absolute error of at most `1/4` (given that `E_5 < 3/4`).
        This means that we need a relative error of at most

        .. MATH::

            \frac{e q}{4 b} > \frac{(e q)/4}{2^p} > \frac{7}{2^p}

        (assuming `n \geq 5`).
        With a precision of `p` bits and rounding down, every rounding
        has a relative error of at most `2^{1-p} = 2/2^p`.
        Since we do 3 roundings (`b` and `q` do not require rounding),
        we get a relative error of at most `6/2^p`.
        All this implies that the precision of `p` bits is sufficient.

    REFERENCES:

    - :wikipedia:`Bell_number`
    - http://fredrik-j.blogspot.com/2009/03/computing-generalized-bell-numbers.html
    - http://mathworld.wolfram.com/DobinskisFormula.html
    """
    n = ZZ(n)
    if n < 0:
        raise ArithmeticError('Bell numbers not defined for negative indices')
    if algorithm == 'mpmath':
        from sage.libs.mpmath.all import bell, mp, mag
        old_prec = mp.dps
        if 'prec' in options:
            mp.dps = options['prec']
            ret = ZZ(int(bell(n)))
            mp.dps = old_prec
            return ret
        ret_mp = bell(n)
        p = mag(ret_mp) + 10
        if p > mp.dps:
            mp.dps = p
            ret = ZZ(int(bell(n)))
            mp.dps = old_prec
            return ret
        return ZZ(int(ret_mp))

    elif algorithm == 'flint':
        import sage.libs.flint.arith_sage
        return sage.libs.flint.arith_sage.bell_number(n)

    elif algorithm == 'gap':
        from sage.libs.gap.libgap import libgap
        return libgap.Bell(n).sage()

    elif algorithm == 'dobinski':
        # Hardcode small cases. We only proved the algorithm below
        # for n >= 5, but it turns out that n = 4 also works.
        if n < 4:
            return Integer((1, 1, 2, 5)[n])
        b = ZZ.zero()
        fact = k = ZZ.one()
        q = 2 * n + 1
        si = Integer(n).sqrtrem()[0] // 2
        while True:
            partfact = ZZ.one()
            v = ZZ.zero()
            for i in range(si - 1, -1, -1):
                v += partfact * (k + i)**n
                partfact *= k + i
            fact *= partfact
            v = (q * v) // fact
            if not v:
                break
            b += v
            k += si
        from sage.rings.real_mpfr import RealField
        R = RealField(b.exact_log(2) + 1, rnd='RNDD')
        return ((R(-1).exp() / q) * b).ceil()

    raise ValueError("unknown algorithm %r" % algorithm)


def catalan_number(n):
    r"""
    Return the `n`-th Catalan number.

    The `n`-th Catalan number is given
    directly in terms of binomial coefficients by

    .. MATH::

        C_n = \frac{1}{n+1}\binom{2n}{n} = \frac{(2n)!}{(n+1)!\,n!}
        \qquad\mbox{ for }\quad n\ge 0.

    Consider the set `S = \{ 1, ..., n \}`. A noncrossing
    partition of `S` is a partition in which no two blocks
    "cross" each other, i.e., if `a` and `b` belong to one block and
    `x` and `y` to another, they are not arranged in the order `axby`.
    `C_n` is the number of noncrossing partitions of the set
    `S`. There are many other interpretations (see REFERENCES).

    When `n=-1`, this function returns the limit value `-1/2`. For
    other `n<0` it returns `0`.

    INPUT:

    - ``n`` -- integer

    OUTPUT: integer

    EXAMPLES::

        sage: [catalan_number(i) for i in range(7)]
        [1, 1, 2, 5, 14, 42, 132]
        sage: x = (QQ[['x']].0).O(8)
        sage: (-1/2)*sqrt(1 - 4*x)
        -1/2 + x + x^2 + 2*x^3 + 5*x^4 + 14*x^5 + 42*x^6 + 132*x^7 + O(x^8)
        sage: [catalan_number(i) for i in range(-7,7)]
        [0, 0, 0, 0, 0, 0, -1/2, 1, 1, 2, 5, 14, 42, 132]
        sage: [catalan_number(n).mod(2) for n in range(16)]
        [1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1]

    REFERENCES:

    -  :wikipedia:`Catalan_number`

    -  http://www-history.mcs.st-andrews.ac.uk/~history/Miscellaneous/CatalanNumbers/catalan.html
    """
    n = ZZ(n)
    if n < -1:
        return ZZ.zero()
    if n == -1:
        return QQ((-1, 2))
    return (2 * n).binomial(n).divide_knowing_divisible_by(n + 1)


def narayana_number(n: Integer, k: Integer) -> Integer:
    r"""
    Return the Narayana number of index ``(n, k)``.

    For every integer `n \geq 1`, the sum of Narayana numbers `\sum_k N_{n,k}`
    is the Catalan number `C_n`.

    INPUT:

    - ``n`` -- integer

    - ``k`` -- integer between ``0`` and ``n - 1``

    OUTPUT: integer

    EXAMPLES::

       sage: from sage.combinat.combinat import narayana_number
       sage: [narayana_number(3, i) for i in range(3)]
       [1, 3, 1]
       sage: sum(narayana_number(7,i) for i in range(7)) == catalan_number(7)
       True

    REFERENCES:

    - :wikipedia:`Narayana_number`
    """
    n = ZZ(n)
    if n <= 0:
        return ZZ.zero()
    return (n.binomial(k + 1) * n.binomial(k)).divide_knowing_divisible_by(n)


def euler_number(n, algorithm='flint') -> Integer:
    """
    Return the `n`-th Euler number.

    INPUT:

    - ``n`` -- positive integer

    - ``algorithm`` -- (default: ``'flint'``) any one of the following:

      - ``'maxima'`` -- wraps Maxima's ``euler``

      - ``'flint'`` -- wrap FLINT's ``arith_euler_number``

    EXAMPLES::

        sage: [euler_number(i) for i in range(10)]                                      # needs sage.libs.flint
        [1, 0, -1, 0, 5, 0, -61, 0, 1385, 0]
        sage: x = PowerSeriesRing(QQ, 'x').gen().O(10)
        sage: 2/(exp(x)+exp(-x))
        1 - 1/2*x^2 + 5/24*x^4 - 61/720*x^6 + 277/8064*x^8 + O(x^10)
        sage: [euler_number(i)/factorial(i) for i in range(11)]                         # needs sage.libs.flint
        [1, 0, -1/2, 0, 5/24, 0, -61/720, 0, 277/8064, 0, -50521/3628800]
        sage: euler_number(-1)
        Traceback (most recent call last):
        ...
        ValueError: n (=-1) must be a nonnegative integer

    TESTS::

        sage: euler_number(6, 'maxima')                                                 # needs sage.symbolic
        -61

    REFERENCES:

    - :wikipedia:`Euler_number`
    """
    n = ZZ(n)
    if n < 0:
        raise ValueError("n (=%s) must be a nonnegative integer" % n)
    if algorithm == 'maxima':
        return ZZ(maxima.euler(n))  # type:ignore
    elif algorithm == 'flint':
        import sage.libs.flint.arith_sage
        return sage.libs.flint.arith_sage.euler_number(n)
    else:
        raise ValueError("algorithm must be 'flint' or 'maxima'")


@cached_function(key=lambda n, k, a: (n, k))
def eulerian_number(n, k, algorithm='recursive') -> Integer:
    """
    Return the Eulerian number of index ``(n, k)``.

    This is the coefficient of `t^k` in the Eulerian polynomial `A_n(t)`.

    INPUT:

    - ``n`` -- integer
    - ``k`` -- integer between ``0`` and ``n - 1``
    - ``algorithm`` -- ``'recursive'`` (default) or ``'formula'``

    OUTPUT: integer

    .. SEEALSO:: :func:`eulerian_polynomial`

    EXAMPLES::

        sage: from sage.combinat.combinat import eulerian_number
        sage: [eulerian_number(5,i) for i in range(5)]
        [1, 26, 66, 26, 1]

    TESTS::

        sage: [eulerian_number(6,i,"formula") for i in range(6)]
        [1, 57, 302, 302, 57, 1]
        sage: [eulerian_number(3,i) for i in range(-1, 4)]
        [0, 1, 4, 1, 0]
    """
    n = ZZ(n)
    if k < 0 or k > n - 1:
        return ZZ.zero()
    if k == 0 or k == n - 1:
        return ZZ.one()
    if algorithm == "recursive":
        s = (n - k) * eulerian_number(n - 1, k - 1, algorithm=algorithm)
        s += (k + 1) * eulerian_number(n - 1, k, algorithm=algorithm)
        return s
    return sum((-1)**m * (n + 1).binomial(m) * (k + 1 - m)**n
               for m in range(k + 1))


@cached_function(key=lambda n, a: n)
def eulerian_polynomial(n, algorithm='derivative'):
    """
    Return the Eulerian polynomial of index ``n``.

    This is the generating polynomial counting permutations in the
    symmetric group `S_n` according to their number of descents.

    INPUT:

    - ``n`` -- integer

    - ``algorithm`` -- ``'derivative'`` (default) or ``'coeffs'``

    OUTPUT: polynomial in one variable ``t``

    .. SEEALSO:: :func:`eulerian_number`

    EXAMPLES::

        sage: from sage.combinat.combinat import eulerian_polynomial
        sage: eulerian_polynomial(5)
        t^4 + 26*t^3 + 66*t^2 + 26*t + 1

    TESTS::

        sage: eulerian_polynomial(7)(1) == factorial(7)
        True

        sage: eulerian_polynomial(6, algorithm='coeffs')
        t^5 + 57*t^4 + 302*t^3 + 302*t^2 + 57*t + 1

    REFERENCES:

    - :wikipedia:`Eulerian_number`
    """
    n = ZZ(n)
    R = PolynomialRing(ZZ, 't')
    if n < 0:
        return R.zero()
    if n == 1:
        return R.one()
    t = R.gen()
    if algorithm == 'derivative':
        A = eulerian_polynomial(n - 1, algorithm=algorithm)
        return t * (1 - t) * A.derivative() + (1 + (n - 1) * t) * A
    elif algorithm == 'coeffs':
        return R([eulerian_number(n, k, "formula") for k in range(n)])


def fibonacci(n, algorithm='pari') -> Integer:
    """
    Return the `n`-th Fibonacci number.

    The Fibonacci sequence `F_n` is defined by the initial
    conditions `F_1 = F_2 = 1` and the recurrence relation
    `F_{n+2} = F_{n+1} + F_n`. For negative `n` we
    define `F_n = (-1)^{n+1}F_{-n}`, which is consistent with
    the recurrence relation.

    INPUT:

    - ``algorithm`` -- string; one of

      * ``'pari'`` -- (default) use the PARI C library's
        :pari:`fibo` function

      * ``'gap'`` -- use GAP's Fibonacci function

    .. NOTE::

        PARI is tens to hundreds of times faster than GAP here.
        Moreover, PARI works for every large input whereas GAP does not.

    EXAMPLES::

        sage: fibonacci(10)                                                             # needs sage.libs.pari
        55
        sage: fibonacci(10, algorithm='gap')                                            # needs sage.libs.gap
        55

    ::

        sage: fibonacci(-100)                                                           # needs sage.libs.pari
        -354224848179261915075
        sage: fibonacci(100)                                                            # needs sage.libs.pari
        354224848179261915075

    ::

        sage: fibonacci(0)                                                              # needs sage.libs.pari
        0
        sage: fibonacci(1/2)
        Traceback (most recent call last):
        ...
        TypeError: no conversion of this rational to integer
    """
    n = ZZ(n)
    if algorithm == 'pari':
        return ZZ(pari(n).fibonacci())
    elif algorithm == 'gap':
        from sage.libs.gap.libgap import libgap
        return libgap.Fibonacci(n).sage()
    else:
        raise ValueError("no algorithm {}".format(algorithm))


def lucas_number1(n, P, Q):
    r"""
    Return the `n`-th Lucas number "of the first kind" (this is not
    standard terminology). The Lucas sequence `L^{(1)}_n` is
    defined by the initial conditions `L^{(1)}_1 = 0`,
    `L^{(1)}_2 = 1` and the recurrence relation
    `L^{(1)}_{n+2} = P \cdot L^{(1)}_{n+1} - Q \cdot L^{(1)}_n`.

    Wraps GAP's ``Lucas(...)[1]``.

    `P=1`, `Q=-1` gives the Fibonacci sequence.

    INPUT:

    - ``n`` -- integer

    - ``P``, ``Q`` -- integer or rational numbers

    OUTPUT: integer or rational number

    EXAMPLES::

        sage: # needs sage.libs.gap
        sage: lucas_number1(5,1,-1)
        5
        sage: lucas_number1(6,1,-1)
        8
        sage: lucas_number1(7,1,-1)
        13
        sage: lucas_number1(7,1,-2)
        43
        sage: lucas_number1(5,2,3/5)
        229/25
        sage: lucas_number1(5,2,1.5)
        1/4

    There was a conjecture that the sequence `L_n` defined by
    `L_{n+2} = L_{n+1} + L_n`, `L_1=1`,
    `L_2=3`, has the property that `n` prime implies
    that `L_n` is prime. ::

        sage: def lucas(n):
        ....:     return Integer((5/2)*lucas_number1(n,1,-1) + (1/2)*lucas_number2(n,1,-1))
        sage: [[lucas(n), is_prime(lucas(n)), n+1, is_prime(n+1)] for n in range(15)]   # needs sage.libs.gap
        [[1, False, 1, False],
         [3, True, 2, True],
         [4, False, 3, True],
         [7, True, 4, False],
         [11, True, 5, True],
         [18, False, 6, False],
         [29, True, 7, True],
         [47, True, 8, False],
         [76, False, 9, False],
         [123, False, 10, False],
         [199, True, 11, True],
         [322, False, 12, False],
         [521, True, 13, True],
         [843, False, 14, False],
         [1364, False, 15, False]]

    Can you use Sage to find a counterexample to the conjecture?
    """
    n = ZZ(n)
    P = QQ(P)
    Q = QQ(Q)
    from sage.libs.gap.libgap import libgap
    return libgap.Lucas(P, Q, n)[0].sage()


def lucas_number2(n, P, Q):
    r"""
    Return the `n`-th Lucas number "of the second kind" (this is not
    standard terminology). The Lucas sequence `L^{(2)}_n` is
    defined by the initial conditions `L^{(2)}_1 = 2`,
    `L^{(2)}_2 = P` and the recurrence relation
    `L^{(2)}_{n+2} = P \cdot L^{(2)}_{n+1} - Q \cdot L^{(2)}_n`.

    Wraps GAP's Lucas(...)[2].

    INPUT:

    - ``n`` -- integer

    - ``P``, ``Q`` -- integer or rational numbers

    OUTPUT: integer or rational number

    EXAMPLES::

        sage: [lucas_number2(i,1,-1) for i in range(10)]                                # needs sage.libs.gap
        [2, 1, 3, 4, 7, 11, 18, 29, 47, 76]
        sage: [fibonacci(i-1)+fibonacci(i+1) for i in range(10)]                        # needs sage.libs.pari
        [2, 1, 3, 4, 7, 11, 18, 29, 47, 76]

    ::

        sage: # needs sage.libs.gap
        sage: n = lucas_number2(5,2,3); n
        2
        sage: type(n)
        <class 'sage.rings.integer.Integer'>
        sage: n = lucas_number2(5,2,-3/9); n
        418/9
        sage: type(n)
        <class 'sage.rings.rational.Rational'>

    The case `P=1`, `Q=-1` is the Lucas sequence in Brualdi's Introductory
    Combinatorics, 4th ed., Prentice-Hall, 2004::

        sage: [lucas_number2(n,1,-1) for n in range(10)]                                # needs sage.libs.gap
        [2, 1, 3, 4, 7, 11, 18, 29, 47, 76]
    """
    n = ZZ(n)
    P = QQ(P)
    Q = QQ(Q)
    from sage.libs.gap.libgap import libgap
    return libgap.Lucas(P, Q, n)[1].sage()


def stirling_number1(n, k, algorithm='gap') -> Integer:
    r"""
    Return the `n`-th Stirling number `S_1(n,k)` of the first kind.

    This is the number of permutations of `n` points with `k` cycles.

    See :wikipedia:`Stirling_numbers_of_the_first_kind`.

    INPUT:

    - ``n`` -- nonnegative machine-size integer
    - ``k`` -- nonnegative machine-size integer
    - ``algorithm``:

      * ``'gap'`` -- default; use GAP's ``Stirling1`` function
      * ``'flint'`` -- use flint's ``arith_stirling_number_1u`` function

    EXAMPLES::

        sage: # needs sage.libs.gap
        sage: stirling_number1(3,2)
        3
        sage: stirling_number1(5,2)
        50
        sage: 9*stirling_number1(9,5) + stirling_number1(9,4)
        269325
        sage: stirling_number1(10,5)
        269325

    Indeed, `S_1(n,k) = S_1(n-1,k-1) + (n-1)S_1(n-1,k)`.

    TESTS::

        sage: stirling_number1(10,5, algorithm='flint')                                 # needs sage.libs.flint
        269325

        sage: s_sage = stirling_number1(50,3, algorithm='mutta')
        Traceback (most recent call last):
        ...
        ValueError: unknown algorithm: mutta
    """
    n = ZZ(n)
    k = ZZ(k)
    if k > n:
        return ZZ.zero()
    if k == 0:
        return ZZ.zero() if n else ZZ.one()
    if algorithm == 'gap':
        from sage.libs.gap.libgap import libgap
        return libgap.Stirling1(n, k).sage()
    if algorithm == 'flint':
        import sage.libs.flint.arith_sage
        return sage.libs.flint.arith_sage.stirling_number_1(n, k)
    raise ValueError("unknown algorithm: %s" % algorithm)


def stirling_number2(n, k, algorithm=None) -> Integer:
    r"""
    Return the `n`-th Stirling number `S_2(n,k)` of the second kind.

    This is the number of ways to partition a set of `n` elements into `k`
    pairwise disjoint nonempty subsets. The `n`-th Bell number is the
    sum of the `S_2(n,k)`'s, `k=0,...,n`.

    See :wikipedia:`Stirling_numbers_of_the_second_kind`.

    INPUT:

    - ``n`` -- nonnegative machine-size integer
    - ``k`` -- nonnegative machine-size integer
    - ``algorithm``:

      * ``None`` -- default; use native implementation
      * ``'flint'`` -- use flint's ``arith_stirling_number_2`` function
      * ``'gap'`` -- use GAP's ``Stirling2`` function
      * ``'maxima'`` -- use Maxima's ``stirling2`` function

    EXAMPLES:

    Print a table of the first several Stirling numbers of the second kind::

        sage: for n in range(10):
        ....:     for k in range(10):
        ....:         print(str(stirling_number2(n,k)).rjust(k and 6))
        1      0      0      0      0      0      0      0      0      0
        0      1      0      0      0      0      0      0      0      0
        0      1      1      0      0      0      0      0      0      0
        0      1      3      1      0      0      0      0      0      0
        0      1      7      6      1      0      0      0      0      0
        0      1     15     25     10      1      0      0      0      0
        0      1     31     90     65     15      1      0      0      0
        0      1     63    301    350    140     21      1      0      0
        0      1    127    966   1701   1050    266     28      1      0
        0      1    255   3025   7770   6951   2646    462     36      1

    Stirling numbers satisfy `S_2(n,k) = S_2(n-1,k-1) + kS_2(n-1,k)`::

        sage: 5*stirling_number2(9,5) + stirling_number2(9,4)
        42525
        sage: stirling_number2(10,5)
        42525

    TESTS::

        sage: stirling_number2(500,501)
        0
        sage: stirling_number2(500,500)
        1
        sage: stirling_number2(500,499)
        124750
        sage: stirling_number2(500,498)
        7739801875
        sage: stirling_number2(500,497)
        318420320812125
        sage: stirling_number2(500,0)
        0
        sage: stirling_number2(500,1)
        1
        sage: stirling_number2(500,2)
        1636695303948070935006594848413799576108321023021532394741645684048066898202337277441635046162952078575443342063780035504608628272942696526664263794687
        sage: stirling_number2(500,3)
        6060048632644989473730877846590553186337230837666937173391005972096766698597315914033083073801260849147094943827552228825899880265145822824770663507076289563105426204030498939974727520682393424986701281896187487826395121635163301632473646
        sage: stirling_number2(500,30)
        13707767141249454929449108424328432845001327479099713037876832759323918134840537229737624018908470350134593241314462032607787062188356702932169472820344473069479621239187226765307960899083230982112046605340713218483809366970996051181537181362810003701997334445181840924364501502386001705718466534614548056445414149016614254231944272872440803657763210998284198037504154374028831561296154209804833852506425742041757849726214683321363035774104866182331315066421119788248419742922490386531970053376982090046434022248364782970506521655684518998083846899028416459701847828711541840099891244700173707021989771147674432503879702222276268661726508226951587152781439224383339847027542755222936463527771486827849728880
        sage: stirling_number2(500,31)
        5832088795102666690960147007601603328246123996896731854823915012140005028360632199516298102446004084519955789799364757997824296415814582277055514048635928623579397278336292312275467402957402880590492241647229295113001728653772550743446401631832152281610081188041624848850056657889275564834450136561842528589000245319433225808712628826136700651842562516991245851618481622296716433577650218003181535097954294609857923077238362717189185577756446945178490324413383417876364657995818830270448350765700419876347023578011403646501685001538551891100379932684279287699677429566813471166558163301352211170677774072447414719380996777162087158124939742564291760392354506347716119002497998082844612434332155632097581510486912
        sage: n = stirling_number2(20,11); n
        1900842429486
        sage: type(n)
        <class 'sage.rings.integer.Integer'>
        sage: n_gap = stirling_number2(20, 11, algorithm='gap'); n_gap                  # needs sage.libs.gap
        1900842429486
        sage: type(n_gap)                                                               # needs sage.libs.gap
        <class 'sage.rings.integer.Integer'>
        sage: n_flint = stirling_number2(20, 11, algorithm='flint'); n_flint            # needs sage.libs.flint
        1900842429486
        sage: type(n_flint)                                                             # needs sage.libs.flint
        <class 'sage.rings.integer.Integer'>

    Sage's implementation splitting the computation of the Stirling
    numbers of the second kind in two cases according to `n`, let us
    check the result it gives agree with both flint and gap.

    For `n<200`::

        sage: for n in Subsets(range(100,200), 5).random_element():                     # needs sage.libs.flint sage.libs.gap
        ....:     for k in Subsets(range(n), 5).random_element():
        ....:         s_sage = stirling_number2(n,k)
        ....:         s_flint = stirling_number2(n,k, algorithm = "flint")
        ....:         s_gap = stirling_number2(n,k, algorithm = "gap")
        ....:         if not (s_sage == s_flint and s_sage == s_gap):
        ....:             print("Error with n<200")

    For `n\geq 200`::

        sage: for n in Subsets(range(200,300), 5).random_element():                     # needs sage.libs.flint sage.libs.gap
        ....:     for k in Subsets(range(n), 5).random_element():
        ....:         s_sage = stirling_number2(n,k)
        ....:         s_flint = stirling_number2(n,k, algorithm = "flint")
        ....:         s_gap = stirling_number2(n,k, algorithm = "gap")
        ....:         if not (s_sage == s_flint and s_sage == s_gap):
        ....:             print("Error with n<200")

        sage: stirling_number2(20, 3, algorithm='maxima')                               # needs sage.symbolic
        580606446

        sage: s_sage = stirling_number2(5, 3, algorithm='namba')
        Traceback (most recent call last):
        ...
        ValueError: unknown algorithm: namba
    """
    n = ZZ(n)
    k = ZZ(k)
    if k > n:
        return ZZ.zero()
    if k == 0:
        return ZZ.zero() if n else ZZ.one()
    if algorithm is None:
        return _stirling_number2(n, k)
    if algorithm == 'gap':
        from sage.libs.gap.libgap import libgap
        return libgap.Stirling2(n, k).sage()
    if algorithm == 'flint':
        import sage.libs.flint.arith_sage
        return sage.libs.flint.arith_sage.stirling_number_2(n, k)
    if algorithm == 'maxima':
        return ZZ(maxima.stirling2(n, k))  # type:ignore
    raise ValueError("unknown algorithm: %s" % algorithm)


def lah_number(n, k) -> Integer:
    r"""
    Return the Lah number `L(n,k)`

    This is the number of ways to partition a set of `n` elements into `k`
    pairwise disjoint nonempty linearly-ordered subsets.

    This is also called the Stirling number of the third kind.

    See :wikipedia:`Lah_number`.

    INPUT:

    - ``n`` -- nonnegative integer

    - ``k`` -- nonnegative integer

    EXAMPLES::

        sage: from sage.combinat.combinat import lah_number
        sage: lah_number(50, 30)
        3242322638238907670866645288893161825894400000

    We verify a well-known identity::

        sage: S1 = stirling_number1; S2 = stirling_number2
        sage: all(lah_number(n, k) == sum(S1(n, j) * S2(j, k) for j in [k..n])
        ....:     for n in range(10) for k in range(10))
        True

    TESTS:

    Verify the usual convention for the degenerate case `k = 0`::

        sage: lah_number(0, 0) == 1 and lah_number(1, 0) == 0
        True
    """
    n = ZZ(n)
    k = ZZ(k)
    if k > n:
        return ZZ.zero()
    if k == 0:
        return ZZ.zero() if n else ZZ.one()
    a = n.factorial() // k.factorial()
    return a**2 * k // (n * (n - k).factorial())


def polygonal_number(s, n):
    r"""
    Return the `n`-th `s`-gonal number.

    Polygonal sequences are represented by dots forming a regular polygon.
    Two famous sequences are the triangular numbers (3rd column of Pascal's
    Triangle) and the square numbers. The `n`-th term in a polygonal sequence
    is defined by

    .. MATH::

        P(s, n) = \frac{n^2(s-2) - n(s-4)}{2},

    where `s` is the number of sides of the polygon.

    INPUT:

    - ``s`` -- integer greater than 1; the number of sides of the polygon

    - ``n`` -- integer; the index of the returned `s`-gonal number

    OUTPUT: integer

    EXAMPLES:

    The triangular numbers::

        sage: [polygonal_number(3, n) for n in range(10)]
        [0, 1, 3, 6, 10, 15, 21, 28, 36, 45]

        sage: [polygonal_number(3, n) for n in range(-10, 0)]
        [45, 36, 28, 21, 15, 10, 6, 3, 1, 0]

    The square numbers::

        sage: [polygonal_number(4, n) for n in range(10)]
        [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]

    The pentagonal numbers::

        sage: [polygonal_number(5, n) for n in range(10)]
        [0, 1, 5, 12, 22, 35, 51, 70, 92, 117]

    The hexagonal numbers::

        sage: [polygonal_number(6, n) for n in range(10)]
        [0, 1, 6, 15, 28, 45, 66, 91, 120, 153]

    The input is converted into an integer::

        sage: polygonal_number(3.0, 2.0)
        3

    A non-integer input returns an error::

        sage: polygonal_number(3.5, 1)                                                  # needs sage.rings.real_mpfr
        Traceback (most recent call last):
        ...
        TypeError: Attempt to coerce non-integral RealNumber to Integer

    `s` must be greater than 1::

        sage: polygonal_number(1, 4)
        Traceback (most recent call last):
        ...
        ValueError: s (=1) must be greater than 1

    REFERENCES:

    - :wikipedia:`Polygonal_number`
    """
    s = ZZ(s)
    n = ZZ(n)
    if s < 2:
        raise ValueError("s (=%s) must be greater than 1" % s)
    return (((n**2) * (s - 2)) - (n * (s - 4))) // 2


class CombinatorialObject(SageObject):
    def __init__(self, l, copy=True):
        """
        CombinatorialObject provides a thin wrapper around a list. The main
        differences are that __setitem__ is disabled so that
        CombinatorialObjects are shallowly immutable, and the intention is
        that they are semantically immutable.

        Because of this, CombinatorialObjects provide a __hash__
        function which computes the hash of the string representation of a
        list and the hash of its parent's class. Thus, each
        CombinatorialObject should have a unique string representation.

        .. SEEALSO::

            :class:`CombinatorialElement` if you want a combinatorial
            object which is an element of a parent.

        .. WARNING::

            This class is slowly being deprecated. Use
            :class:`~sage.structure.list_clone.ClonableList` instead.

        INPUT:

        - ``l`` -- list or any object that can be converted to a
          list by calling ``list()``

        - ``copy`` -- boolean (default: ``True``); if ``False``, then
          ``l`` must be a ``list``, which is assigned to ``self._list``
          without copying

        EXAMPLES::

            sage: c = CombinatorialObject([1,2,3])
            sage: c == loads(dumps(c))
            True
            sage: c._list
            [1, 2, 3]
            sage: c._hash is None
            True

        For efficiency, you can specify ``copy=False`` if you know what
        you are doing::

            sage: from sage.combinat.combinat import CombinatorialObject
            sage: x = [3, 2, 1]
            sage: C = CombinatorialObject(x, copy=False)
            sage: C
            [3, 2, 1]
            sage: x[0] = 5
            sage: C
            [5, 2, 1]

        TESTS:

        Test indirectly that we copy the input (see :issue:`18184`)::

            sage: # needs sage.combinat
            sage: L = IntegerListsLex(element_class=Partition)
            sage: x = [3, 2, 1]
            sage: P = L(x)
            sage: x[0] = 5
            sage: list(P)
            [3, 2, 1]
        """
        if copy:
            self._list = list(l)
        else:
            self._list = l
        self._hash = None

    def __str__(self):
        """
        EXAMPLES::

            sage: c = CombinatorialObject([1,2,3])
            sage: str(c)
            '[1, 2, 3]'
        """
        return str(self._list)

    def _repr_(self):
        """
        EXAMPLES::

            sage: c = CombinatorialObject([1,2,3])
            sage: c.__repr__()
            '[1, 2, 3]'
        """
        return repr(self._list)

    def __eq__(self, other):
        """
        Test equality of ``self`` and ``other``.

        EXAMPLES::

            sage: c = CombinatorialObject([1,2,3])
            sage: d = CombinatorialObject([2,3,4])
            sage: c == [1,2,3]
            True
            sage: c == [2,3,4]
            False
            sage: c == d
            False
            sage: c == c
            True

        .. WARNING::

            :class:`CombinatorialObject` must come **before** :class:`Element`
            for this to work because :class:`Element` is ahead of
            :class:`CombinatorialObject` in the MRO (method resolution
            order)::

                sage: from sage.structure.element import Element
                sage: class Bar(Element, CombinatorialObject):
                ....:     def __init__(self, l):
                ....:         CombinatorialObject.__init__(self, l)
                sage: L = [Bar([4-i]) for i in range(4)]
                sage: sorted(L)
                Traceback (most recent call last):
                ...
                TypeError: '<' not supported between instances of 'Bar' and 'Bar'
        """
        if isinstance(other, CombinatorialObject):
            return self._list == other._list
        else:
            return self._list == other

    def __lt__(self, other):
        """
        EXAMPLES::

            sage: c = CombinatorialObject([1,2,3])
            sage: d = CombinatorialObject([2,3,4])
            sage: c < d
            True
            sage: c < [2,3,4]
            True
            sage: c < c
            False

        Check that :issue:`14065` is fixed::

            sage: from sage.structure.element import Element
            sage: class Foo(CombinatorialObject, Element): pass
            sage: L = [Foo([4-i]) for i in range(4)]; L
            [[4], [3], [2], [1]]
            sage: sorted(L)
            [[1], [2], [3], [4]]
            sage: f = Foo([4])
            sage: f is None
            False
            sage: f is not None
            True
        """
        if isinstance(other, CombinatorialObject):
            return self._list < other._list
        else:
            return self._list < other

    def __le__(self, other):
        """
        EXAMPLES::

            sage: c = CombinatorialObject([1,2,3])
            sage: d = CombinatorialObject([2,3,4])
            sage: c <= c
            True
            sage: c <= d
            True
            sage: c <= [1,2,3]
            True
        """
        if isinstance(other, CombinatorialObject):
            return self._list <= other._list
        else:
            return self._list <= other

    def __gt__(self, other):
        """
        EXAMPLES::

            sage: c = CombinatorialObject([1,2,3])
            sage: d = CombinatorialObject([2,3,4])
            sage: c > c
            False
            sage: c > d
            False
            sage: c > [1,2,3]
            False
        """
        if isinstance(other, CombinatorialObject):
            return self._list > other._list
        else:
            return self._list > other

    def __ge__(self, other):
        """
        EXAMPLES::

            sage: c = CombinatorialObject([1,2,3])
            sage: d = CombinatorialObject([2,3,4])
            sage: c >= c
            True
            sage: c >= d
            False
            sage: c >= [1,2,3]
            True
        """
        if isinstance(other, CombinatorialObject):
            return self._list >= other._list
        else:
            return self._list >= other

    def __ne__(self, other):
        """
        EXAMPLES::

            sage: c = CombinatorialObject([1,2,3])
            sage: d = CombinatorialObject([2,3,4])
            sage: c != c
            False
            sage: c != d
            True
            sage: c != [1,2,3]
            False
        """
        if isinstance(other, CombinatorialObject):
            return self._list != other._list
        else:
            return self._list != other

    def __add__(self, other):
        """
        EXAMPLES::

            sage: c = CombinatorialObject([1,2,3])
            sage: c + [4]
            [1, 2, 3, 4]
            sage: type(_)
            <class 'list'>
        """
        return self._list + other

    def __hash__(self):
        """
        Compute the hash of ``self`` by computing the hash of the string
        representation of self._list. The hash is cached and stored in
        self._hash.

        EXAMPLES::

            sage: c = CombinatorialObject([1,2,3])
            sage: c._hash is None
            True
            sage: hash(c) #random
            1335416675971793195
            sage: c._hash #random
            1335416675971793195
        """
        if self._hash is None:
            self._hash = hash(str(self._list))
        return self._hash

    def __bool__(self) -> bool:
        """
        Return ``True`` if ``self`` is nonzero.

        We consider a list to be zero if it has length zero.

        TESTS::

            sage: c = CombinatorialObject([1,2,3])
            sage: not c
            False
            sage: c = CombinatorialObject([])
            sage: not c
            True

        Check that :issue:`14065` is fixed::

            sage: from sage.structure.element import Element
            sage: class Foo(CombinatorialObject, Element): pass
            ...
            sage: f = Foo([4])
            sage: not f
            False
            sage: f = Foo([])
            sage: not f
            True

        .. WARNING::

            :class:`CombinatorialObject` must come **before** :class:`Element`
            for this to work because :class:`Element` is ahead of
            :class:`CombinatorialObject` in the MRO (method resolution
            order)::

                sage: from sage.structure.element import Element
                sage: class Bar(Element, CombinatorialObject):
                ....:     def __init__(self, l):
                ....:         CombinatorialObject.__init__(self, l)
                sage: b = Bar([4])
                sage: not b
                False
        """
        return bool(self._list)

    def __len__(self) -> Integer:
        """
        EXAMPLES::

            sage: c = CombinatorialObject([1,2,3])
            sage: len(c)
            3
            sage: c.__len__()
            3
        """
        return len(self._list)

    def __getitem__(self, key):
        """
        EXAMPLES::

            sage: c = CombinatorialObject([1,2,3])
            sage: c[0]
            1
            sage: c[1:]
            [2, 3]
            sage: type(_)
            <class 'list'>
        """
        return self._list[key]

    def __iter__(self):
        """
        EXAMPLES::

            sage: c = CombinatorialObject([1,2,3])
            sage: list(iter(c))
            [1, 2, 3]
        """
        return iter(self._list)

    def __contains__(self, item):
        """
        EXAMPLES::

            sage: c = CombinatorialObject([1,2,3])
            sage: 1 in c
            True
            sage: 5 in c
            False
        """
        return item in self._list

    def index(self, key):
        """
        EXAMPLES::

            sage: c = CombinatorialObject([1,2,3])
            sage: c.index(1)
            0
            sage: c.index(3)
            2
        """
        return self._list.index(key)


class CombinatorialElement(CombinatorialObject, Element,
        metaclass=InheritComparisonClasscallMetaclass):
    """
    ``CombinatorialElement`` is both a :class:`CombinatorialObject`
    and an :class:`Element`. So it represents a list which is an
    element of some parent.

    A ``CombinatorialElement`` subclass also automatically supports
    the ``__classcall__`` mechanism.

    .. WARNING::

        This class is slowly being deprecated. Use
        :class:`~sage.structure.list_clone.ClonableList` instead.

    INPUT:

    - ``parent`` -- the :class:`Parent` class for this element

    - ``lst`` -- list or any object that can be converted to a
      list by calling ``list()``

    EXAMPLES::

        sage: # needs sage.combinat
        sage: from sage.combinat.combinat import CombinatorialElement
        sage: e = CombinatorialElement(Partitions(6), [3,2,1])
        sage: e == loads(dumps(e))
        True
        sage: parent(e)
        Partitions of the integer 6
        sage: list(e)
        [3, 2, 1]

    Check classcalls::

        sage: class Foo(CombinatorialElement):                                          # needs sage.combinat
        ....:     @staticmethod
        ....:     def __classcall__(cls, x):
        ....:         return x
        sage: Foo(17)                                                                   # needs sage.combinat
        17
    """

    def __init__(self, parent, *args, **kwds):
        """
        Initialize this ``CombinatorialElement`` with a parent and a
        list.

        EXAMPLES::

            sage: from sage.combinat.combinat import CombinatorialElement
            sage: e = CombinatorialElement(ZZ, list=(3,2,1))
            sage: e._list
            [3, 2, 1]
            sage: e.parent()
            Integer Ring

        TESTS::

            sage: CombinatorialElement(ZZ)
            Traceback (most recent call last):
            ...
            TypeError: ...__init__() takes exactly 2 arguments (1 given)
            sage: CombinatorialElement(ZZ, 1, 2)
            Traceback (most recent call last):
            ...
            TypeError: ...__init__() takes exactly 2 arguments (3 given)
            sage: CombinatorialElement(ZZ, 1, list=2)
            Traceback (most recent call last):
            ...
            TypeError: ...__init__() takes exactly 2 arguments (3 given)
            sage: CombinatorialElement(ZZ, a=1, b=2)
            Traceback (most recent call last):
            ...
            TypeError: ...__init__() takes exactly 2 arguments (3 given)
        """
        # There should be one "list" argument, which can be given as
        # positional or keyword argument (in the latter case, the name
        # doesn't matter).
        if len(args) == 1 and not kwds:
            L = args[0]
        elif len(kwds) == 1 and not args:
            L, = kwds.values()
        else:
            raise TypeError("__init__() takes exactly 2 arguments ({} given)".format(1 + len(args) + len(kwds)))
        super().__init__(L)
        super(CombinatorialObject, self).__init__(parent)

#####################################################
# combinatorial sets/lists


def tuples(S, k, algorithm='itertools'):
    r"""
    Return a list of all `k`-tuples of elements of a given set ``S``.

    This function accepts the set ``S`` in the form of any iterable
    (list, tuple or iterator), and returns a list of `k`-tuples.
    If ``S`` contains duplicate entries, then you should expect the
    method to return tuples multiple times!

    Recall that `k`-tuples are ordered (in the sense that two `k`-tuples
    differing in the order of their entries count as different) and
    can have repeated entries (even if ``S`` is a list with no
    repetition).

    INPUT:

    - ``S`` -- the base set
    - ``k`` -- the length of the tuples
    - ``algorithm`` -- can be one of the following:

      * ``'itertools'`` -- (default) use python's itertools
      * ``'native'`` -- use a native Sage implementation

    .. NOTE::

        The ordering of the list of tuples depends on the algorithm.

    EXAMPLES::

        sage: S = [1,2]
        sage: tuples(S,3)
        [(1, 1, 1), (1, 1, 2), (1, 2, 1), (1, 2, 2),
         (2, 1, 1), (2, 1, 2), (2, 2, 1), (2, 2, 2)]
        sage: mset = ["s","t","e","i","n"]
        sage: tuples(mset, 2)
        [('s', 's'), ('s', 't'), ('s', 'e'), ('s', 'i'), ('s', 'n'),
         ('t', 's'), ('t', 't'), ('t', 'e'), ('t', 'i'), ('t', 'n'),
         ('e', 's'), ('e', 't'), ('e', 'e'), ('e', 'i'), ('e', 'n'),
         ('i', 's'), ('i', 't'), ('i', 'e'), ('i', 'i'), ('i', 'n'),
         ('n', 's'), ('n', 't'), ('n', 'e'), ('n', 'i'), ('n', 'n')]

    ::

        sage: K.<a> = GF(4, 'a')                                                        # needs sage.rings.finite_rings
        sage: mset = [x for x in K if x != 0]                                           # needs sage.rings.finite_rings
        sage: tuples(mset, 2)                                                           # needs sage.rings.finite_rings
        [(a, a), (a, a + 1), (a, 1), (a + 1, a), (a + 1, a + 1),
         (a + 1, 1), (1, a), (1, a + 1), (1, 1)]

    We check that the implementations agree (up to ordering)::

        sage: tuples(S, 3, 'native')
        [(1, 1, 1), (2, 1, 1), (1, 2, 1), (2, 2, 1),
         (1, 1, 2), (2, 1, 2), (1, 2, 2), (2, 2, 2)]

    Lastly we check on a multiset::

        sage: S = [1,1,2]
        sage: sorted(tuples(S, 3)) == sorted(tuples(S, 3, 'native'))
        True
    """
    if algorithm == 'itertools':
        import itertools
        return list(itertools.product(S, repeat=k))
    if algorithm == 'native':
        return _tuples_native(S, k)
    raise ValueError('invalid algorithm')


def _tuples_native(S, k):
    """
    Return a list of all `k`-tuples of elements of a given set ``S``.

    This is a helper method used in :meth:`tuples`. It returns the
    same as ``tuples(S, k, algorithm='native')``.

    EXAMPLES::

        sage: S = [1,2,2]
        sage: from sage.combinat.combinat import _tuples_native
        sage: _tuples_native(S,2)
        [(1, 1), (2, 1), (2, 1), (1, 2), (2, 2), (2, 2),
         (1, 2), (2, 2), (2, 2)]
    """
    if k <= 0:
        return [()]
    if k == 1:
        return [(x,) for x in S]
    ans = []
    for s in S:
        for x in _tuples_native(S, k - 1):
            y = list(x)
            y.append(s)
            ans.append(tuple(y))
    return ans


def number_of_tuples(S, k, algorithm='naive') -> Integer:
    """
    Return the size of ``tuples(S, k)`` when `S` is a set. More
    generally, return the size of ``tuples(set(S), k)``. (So,
    unlike :meth:`tuples`, this method removes redundant entries from
    `S`.)

    INPUT:

    - ``S`` -- the base set
    - ``k`` -- the length of the tuples
    - ``algorithm`` -- can be one of the following:

      * ``'naive'`` -- (default) use the naive counting `|S|^k`
      * ``'gap'`` -- wraps GAP's ``NrTuples``

    .. WARNING::

        When using ``algorithm='gap'``, ``S`` must be a list of objects
        that have string representations that can be interpreted by the GAP
        interpreter. If ``S`` consists of at all complicated Sage
        objects, this function might *not* do what you expect.

    EXAMPLES::

        sage: S = [1,2,3,4,5]
        sage: number_of_tuples(S,2)
        25
        sage: number_of_tuples(S,2, algorithm='gap')                                    # needs sage.libs.gap
        25
        sage: S = [1,1,2,3,4,5]
        sage: number_of_tuples(S,2)
        25
        sage: number_of_tuples(S,2, algorithm='gap')                                    # needs sage.libs.gap
        25
        sage: number_of_tuples(S,0)
        1
        sage: number_of_tuples(S,0, algorithm='gap')                                    # needs sage.libs.gap
        1
    """
    if algorithm == 'naive':
        return ZZ(len(set(S)))**k  # The set is there to avoid duplicates
    if algorithm == 'gap':
        k = ZZ(k)
        from sage.libs.gap.libgap import libgap
        S = libgap.eval(str(S))
        return libgap.NrTuples(S, k).sage()
    raise ValueError('invalid algorithm')


def unordered_tuples(S, k, algorithm='itertools'):
    r"""
    Return a list of all unordered tuples of length ``k`` of the set ``S``.

    An unordered tuple of length `k` of set `S` is a unordered selection
    with repetitions of `S` and is represented by a sorted list of length
    `k` containing elements from `S`.

    Unlike :meth:`tuples`, the result of this method does not depend on
    how often an element appears in `S`; only the *set* `S` is being
    used. For example, ``unordered_tuples([1, 1, 1], 2)`` will return
    ``[(1, 1)]``. If you want it to return
    ``[(1, 1), (1, 1), (1, 1)]``, use Python's
    ``itertools.combinations_with_replacement`` instead.

    INPUT:

    - ``S`` -- the base set
    - ``k`` -- the length of the tuples
    - ``algorithm`` -- can be one of the following:

      * ``'itertools'`` -- (default) use python's itertools
      * ``'gap'`` -- wraps GAP's ``UnorderedTuples``

    .. WARNING::

        When using ``algorithm='gap'``, ``S`` must be a list of objects
        that have string representations that can be interpreted by the GAP
        interpreter. If ``S`` consists of at all complicated Sage
        objects, this function might *not* do what you expect.

    EXAMPLES::

        sage: S = [1,2]
        sage: unordered_tuples(S, 3)
        [(1, 1, 1), (1, 1, 2), (1, 2, 2), (2, 2, 2)]

    We check that this agrees with GAP::

        sage: unordered_tuples(S, 3, algorithm='gap')                                   # needs sage.libs.gap
        [(1, 1, 1), (1, 1, 2), (1, 2, 2), (2, 2, 2)]

    We check the result on strings::

        sage: S = ["a","b","c"]
        sage: unordered_tuples(S, 2)
        [('a', 'a'), ('a', 'b'), ('a', 'c'), ('b', 'b'), ('b', 'c'), ('c', 'c')]
        sage: unordered_tuples(S, 2, algorithm='gap')                                   # needs sage.libs.gap
        [('a', 'a'), ('a', 'b'), ('a', 'c'), ('b', 'b'), ('b', 'c'), ('c', 'c')]

    Lastly we check on a multiset::

        sage: S = [1,1,2]
        sage: unordered_tuples(S, 3) == unordered_tuples(S, 3, 'gap')                   # needs sage.libs.gap
        True
        sage: unordered_tuples(S, 3)
        [(1, 1, 1), (1, 1, 2), (1, 2, 2), (2, 2, 2)]
    """
    if algorithm == 'itertools':
        import itertools
        return list(itertools.combinations_with_replacement(sorted(set(S)), k))
    if algorithm == 'gap':
        k = ZZ(k)
        from sage.libs.gap.libgap import libgap
        S = libgap.eval(str(S))
        return [tuple(x) for x in libgap.UnorderedTuples(S, k).sage()]
    raise ValueError('invalid algorithm')


def number_of_unordered_tuples(S, k, algorithm='naive') -> Integer:
    r"""
    Return the size of ``unordered_tuples(S, k)`` when `S` is a set.

    INPUT:

    - ``S`` -- the base set
    - ``k`` -- the length of the tuples
    - ``algorithm`` -- can be one of the following:

      * ``'naive'`` -- (default) use the naive counting `\binom{|S|+k-1}{k}`
      * ``'gap'`` -- wraps GAP's ``NrUnorderedTuples``

    .. WARNING::

        When using ``algorithm='gap'``, ``S`` must be a list of objects
        that have string representations that can be interpreted by the GAP
        interpreter. If ``S`` consists of at all complicated Sage
        objects, this function might *not* do what you expect.

    EXAMPLES::

        sage: S = [1,2,3,4,5]
        sage: number_of_unordered_tuples(S,2)
        15
        sage: number_of_unordered_tuples(S,2, algorithm='gap')                          # needs sage.libs.gap
        15
        sage: S = [1,1,2,3,4,5]
        sage: number_of_unordered_tuples(S,2)
        15
        sage: number_of_unordered_tuples(S,2, algorithm='gap')                          # needs sage.libs.gap
        15
        sage: number_of_unordered_tuples(S,0)
        1
        sage: number_of_unordered_tuples(S,0, algorithm='gap')                          # needs sage.libs.gap
        1
    """
    if algorithm == 'naive':
        return ZZ(len(set(S)) + k - 1).binomial(k)  # The set is there to avoid duplicates
    if algorithm == 'gap':
        k = ZZ(k)
        from sage.libs.gap.libgap import libgap
        S = libgap.eval(str(S))
        return libgap.NrUnorderedTuples(S, k).sage()
    raise ValueError('invalid algorithm')


def unshuffle_iterator(a, one=1) -> Iterator:
    r"""
    Iterate over the unshuffles of a list (or tuple) ``a``, also
    yielding the signs of the respective permutations.

    If `n` and `k` are integers satisfying `0 \leq k \leq n`, then
    a `(k, n-k)`-*unshuffle* means a permutation `\pi \in S_n` such
    that `\pi(1) < \pi(2) < \cdots < \pi(k)` and
    `\pi(k+1) < \pi(k+2) < \cdots < \pi(n)`. This method provides,
    for a list `a = (a_1, a_2, \ldots, a_n)` of length `n`, an iterator
    yielding all pairs:

    .. MATH::

        \Bigl( \bigl( (a_{\pi(1)}, a_{\pi(2)}, \ldots, a_{\pi(k)}),
        (a_{\pi(k+1)}, a_{\pi(k+2)}, \ldots, a_{\pi(n)}) \bigl),
        (-1)^{\pi} \Bigr)

    for all `k \in \{0, 1, \ldots, n\}` and all `(k, n-k)`-unshuffles
    `\pi`. The optional variable ``one`` can be set to a different
    value which results in the `(-1)^{\pi}` component being multiplied
    by said value.

    The iterator does not yield these in order of increasing `k`.

    EXAMPLES::

        sage: from sage.combinat.combinat import unshuffle_iterator
        sage: list(unshuffle_iterator([1, 3, 4]))
        [(((), (1, 3, 4)), 1), (((1,), (3, 4)), 1), (((3,), (1, 4)), -1),
         (((1, 3), (4,)), 1), (((4,), (1, 3)), 1), (((1, 4), (3,)), -1),
         (((3, 4), (1,)), 1), (((1, 3, 4), ()), 1)]
        sage: list(unshuffle_iterator([3, 1]))
        [(((), (3, 1)), 1), (((3,), (1,)), 1), (((1,), (3,)), -1),
         (((3, 1), ()), 1)]
        sage: list(unshuffle_iterator([8]))
        [(((), (8,)), 1), (((8,), ()), 1)]
        sage: list(unshuffle_iterator([]))
        [(((), ()), 1)]
        sage: list(unshuffle_iterator([3, 1], 3/2))
        [(((), (3, 1)), 3/2), (((3,), (1,)), 3/2), (((1,), (3,)), -3/2),
         (((3, 1), ()), 3/2)]
    """
    from sage.combinat.subset import powerset
    n = len(a)
    for I in powerset(range(n)):
        sorted_I = tuple(sorted(I))
        nonI = list(range(n))
        for j in reversed(sorted_I):  # probably optimizable
            nonI.pop(j)
        sorted_nonI = tuple(nonI)
        sign = True
        for i in sorted_I:
            if i % 2:  # aka i % 2 == 1
                sign = not sign
        if len(sorted_I) % 4 > 1:
            sign = not sign
        yield ((tuple([a[i] for i in sorted_I]),
                tuple([a[i] for i in sorted_nonI])),
               (one if sign else - one))


def bell_polynomial(n: Integer, k=None, ordinary=False):
    r"""
    Return the partial (or complete) exponential (or ordinary) Bell polynomial.

    The partial exponential *Bell polynomial* is defined by the formula

    .. MATH::

        B_{n,k}(x_0, x_1, \ldots, x_{n-k}) =
            \sum_{\substack{j_0 + \ldots + j_{n-k} = k \\ 1 j_0 + \ldots + (n-k+1) j_{n-k} = n}}
            \frac{n!}{j_0!j_1!\cdots j_{n-k}!}
            \left(\frac{x_0}{(0+1)!}\right)^{j_0}
            \left(\frac{x_1}{(1+1)!}\right)^{j_1} \cdots
            \left(\frac{x_{n-k}}{(n-k+1)!}\right)^{j_{n-k}}.

    The complete exponential Bell Polynomial is defined as

    .. MATH::

        B_n(x_0, x_1, \ldots, x_{n-k}) =
            \sum_{k=0}^n B_{n,k}(x_0, x_1, \ldots, x_{n-k}).

    The ordinary variant of the partial Bell polynomial is defined by

    .. MATH::

        \hat B_{n,k}(x_0, x_1, \ldots, x_{n-k}) =
            \sum_{\substack{j_0 + \ldots + j_{n-k} = k \\ 1 j_0 + \ldots + (n-k+1) j_{n-k} = n}}
            \binom{k}{j_0, j_1, \ldots, j_{n-k}}
            x_0^{j_0} x_1^{j_1} \cdots x_{n-k}^{j_{n-k}},

    where we have used the multinomial coefficient. The complete version has
    the same definition as its exponential counterpart.

    If we define `f(z) = \sum_{n=1}^\infty x_{n-1} z^n/n!`
    then these are alternative definitions for exponential Bell polynomials

    .. MATH::

        \begin{aligned}
        \exp(f(z)) & = \sum_{n=0}^\infty B_n(x_0, \ldots, x_{n-1}) \frac{z^n}{n!}, \\
        \frac{f(z)^k}{k!} & = \sum_{n=k}^\infty B_{n, k}(x_0, \ldots, x_{n-k}) \frac{z^n}{n!}.
        \end{aligned}

    Defining `g(z) = \sum_{n=1}^\infty x_{n-1} z^n`,
    we have the analoguous alternative definitions

    .. MATH::

        \begin{aligned}
        \frac1{1-f(z)} & = \sum_{n=0}^\infty \hat B_n(x_0, \ldots, x_{n-1}) z^n, \\
        f(z)^k & = \sum_{n=k}^\infty \hat B_{n, k}(x_0, \ldots, x_{n-k}) z^n.
        \end{aligned}

    INPUT:

    - ``k`` -- (optional) if specified, returns the partial Bell
      polynomial, otherwise returns the complete Bell polynomial
    - ``ordinary`` -- boolean (default: ``False``); if ``True``, returns the
      (partial) ordinary Bell polynomial, otherwise returns
      the (partial) exponential Bell polynomial

    EXAMPLES:

    The complete and partial Bell polynomials::

        sage: # needs sage.combinat
        sage: bell_polynomial(3)
        x0^3 + 3*x0*x1 + x2
        sage: bell_polynomial(4)
        x0^4 + 6*x0^2*x1 + 3*x1^2 + 4*x0*x2 + x3
        sage: bell_polynomial(6, 3)
        15*x1^3 + 60*x0*x1*x2 + 15*x0^2*x3
        sage: bell_polynomial(6, 6)
        x0^6

    The ordinary variants are::

        sage: # needs sage.combinat sage.arith
        sage: bell_polynomial(3, ordinary=True)
        x0^3 + 2*x0*x1 + x2
        sage: bell_polynomial(4, ordinary=True)
        x0^4 + 3*x0^2*x1 + x1^2 + 2*x0*x2 + x3
        sage: bell_polynomial(6, 3, True)
        x1^3 + 6*x0*x1*x2 + 3*x0^2*x3
        sage: bell_polynomial(6, 6, True)
        x0^6

    We verify the alternative definition of the different Bell polynomials
    using the functions `f` and `g` given above::

        sage: # needs sage.combinat sage.arith
        sage: n = 6 # positive integer
        sage: k = 4 # positive integer
        sage: R.<x> = InfinitePolynomialRing(QQ)
        sage: PR = PolynomialRing(QQ, 'x', n)
        sage: d = {x[i]: PR.gen(i) for i in range(n)} #substitution dictionnary
        sage: L.<z> = LazyPowerSeriesRing(R)
        sage: f = L(lambda i: x[i-1]/factorial(i), valuation=1)
        sage: all(exp(f)[i].subs(d) * factorial(i) == bell_polynomial(i) for i in range(n+1))
        True
        sage: all((f^k/factorial(k))[i].subs(d) * factorial(i) == bell_polynomial(i, k) for i in range(k, n+k))
        True
        sage: g = L(lambda i: x[i-1], valuation=1)
        sage: all((1/(1-g))[i].subs(d) == bell_polynomial(i, ordinary=True) for i in range(n+1))
        True
        sage: all((g^k)[i].subs(d) == bell_polynomial(i, k, True) for i in range(k, n+k))
        True

    TESTS:

    Check that :issue:`18338` is fixed::

        sage: bell_polynomial(0, 0).parent()                                            # needs sage.combinat
        Univariate Polynomial Ring in x0 over Integer Ring

        sage: for n in (0..4):                                                          # needs sage.combinat
        ....:     print([bell_polynomial(n,k).coefficients() for k in (0..n)])
        [[1]]
        [[], [1]]
        [[], [1], [1]]
        [[], [1], [3], [1]]
        [[], [1], [3, 4], [6], [1]]

    Further checks for :issue:`37727`::

        sage: # needs sage.combinat sage.arith
        sage: bell_polynomial(0, 0)
        1
        sage: bell_polynomial(0, 0, True)
        1
        sage: bell_polynomial(1, 1)
        x0
        sage: bell_polynomial(2, 2, True)
        x0^2
        sage: bell_polynomial(5)
        x0^5 + 10*x0^3*x1 + 15*x0*x1^2 + 10*x0^2*x2 + 10*x1*x2 + 5*x0*x3 + x4
        sage: sum(bell_polynomial(5, k) for k in range(6))
        x0^5 + 10*x0^3*x1 + 15*x0*x1^2 + 10*x0^2*x2 + 10*x1*x2 + 5*x0*x3 + x4
        sage: bell_polynomial(5, None, True)
        x0^5 + 4*x0^3*x1 + 3*x0*x1^2 + 3*x0^2*x2 + 2*x1*x2 + 2*x0*x3 + x4
        sage: sum(bell_polynomial(5, k, True) for k in range(6))
        x0^5 + 4*x0^3*x1 + 3*x0*x1^2 + 3*x0^2*x2 + 2*x1*x2 + 2*x0*x3 + x4
        sage: bell_polynomial(0).parent()
        Univariate Polynomial Ring in x0 over Integer Ring

    REFERENCES:

    - [Bel1927]_
    - [Com1974]_
    """
    from sage.combinat.partition import Partitions
    from sage.arith.misc import multinomial
    if k is None:
        partitions = Partitions(n)
        # We set k = 1 to use the correct ring
        # It is not used in the computation otherwise
        k = 1
    else:
        partitions = Partitions(n, length=k)
    if n <= k:
        R = PolynomialRing(ZZ, 'x0')
    else:
        R = PolynomialRing(ZZ, 'x', n - k + 1)
    vars = R.gens()
    result = R.zero()
    for p in partitions:
        if ordinary:
            coefficient = multinomial(p.to_exp())
        else:
            factorial_product = 1
            for part, count in p.to_exp_dict().items():
                factorial_product *= factorial(count) * factorial(part)**count
            coefficient = factorial(n) // factorial_product
        result += coefficient * prod(vars[i - 1] for i in p)
    return result


def fibonacci_sequence(start, stop=None, algorithm=None) -> Iterator:
    r"""
    Return an iterator over the Fibonacci sequence, for all fibonacci
    numbers `f_n` from ``n = start`` up to (but
    not including) ``n = stop``

    INPUT:

    - ``start`` -- starting value

    - ``stop`` -- stopping value

    - ``algorithm`` -- (default: ``None``) passed on to fibonacci function (or
      not passed on if ``None``, i.e., use the default)

    EXAMPLES::

        sage: fibs = [i for i in fibonacci_sequence(10, 20)]; fibs                      # needs sage.libs.pari
        [55, 89, 144, 233, 377, 610, 987, 1597, 2584, 4181]

    ::

        sage: sum([i for i in fibonacci_sequence(100, 110)])                            # needs sage.libs.pari
        69919376923075308730013

    .. SEEALSO::

       :func:`fibonacci_xrange`
    """
    if stop is None:
        stop = ZZ(start)
        start = ZZ(0)
    else:
        start = ZZ(start)
        stop = ZZ(stop)

    if algorithm:
        for n in range(start, stop):
            yield fibonacci(n, algorithm=algorithm)
    else:
        for n in range(start, stop):
            yield fibonacci(n)


def fibonacci_xrange(start, stop=None, algorithm='pari') -> Iterator:
    r"""
    Return an iterator over all of the Fibonacci numbers in the given
    range, including ``f_n = start`` up to, but not
    including, ``f_n = stop``.

    EXAMPLES::

        sage: fibs_in_some_range = [i for i in fibonacci_xrange(10^7, 10^8)]            # needs sage.libs.pari
        sage: len(fibs_in_some_range)                                                   # needs sage.libs.pari
        4
        sage: fibs_in_some_range                                                        # needs sage.libs.pari
        [14930352, 24157817, 39088169, 63245986]

    ::

        sage: fibs = [i for i in fibonacci_xrange(10, 100)]; fibs                       # needs sage.libs.pari
        [13, 21, 34, 55, 89]

    ::

        sage: list(fibonacci_xrange(13, 34))                                            # needs sage.libs.pari
        [13, 21]

    A solution to the second Project Euler problem::

        sage: sum([i for i in fibonacci_xrange(10^6) if is_even(i)])                    # needs sage.libs.pari
        1089154

    .. SEEALSO::

       :func:`fibonacci_sequence`
    """
    if stop is None:
        stop = ZZ(start)
        start = ZZ(0)
    else:
        start = ZZ(start)
        stop = ZZ(stop)

    # iterate until we've gotten high enough
    fn = 0
    n = 0
    while fn < start:
        n += 1
        fn = fibonacci(n)

    while True:
        fn = fibonacci(n)
        n += 1
        if fn < stop:
            yield fn
        else:
            return


def bernoulli_polynomial(x, n: Integer):
    r"""
    Return the ``n``-th Bernoulli polynomial evaluated at ``x``.

    The generating function for the Bernoulli polynomials is

    .. MATH::

       \frac{t e^{xt}}{e^t-1}= \sum_{n=0}^\infty B_n(x) \frac{t^n}{n!},

    and they are given directly by

    .. MATH::

       B_n(x) = \sum_{i=0}^n \binom{n}{i}B_{n-i}x^i.

    One has `B_n(x) = - n\zeta(1 - n,x)`, where `\zeta(s,x)` is the Hurwitz
    zeta function. Thus, in a certain sense, the Hurwitz zeta function
    generalizes the Bernoulli polynomials to non-integer values of `n`.

    EXAMPLES::

        sage: # needs sage.libs.flint
        sage: y = QQ['y'].0
        sage: bernoulli_polynomial(y, 5)
        y^5 - 5/2*y^4 + 5/3*y^3 - 1/6*y
        sage: bernoulli_polynomial(y, 5)(12)
        199870
        sage: bernoulli_polynomial(12, 5)
        199870
        sage: bernoulli_polynomial(y^2 + 1, 5)
        y^10 + 5/2*y^8 + 5/3*y^6 - 1/6*y^2
        sage: P.<t> = ZZ[]
        sage: p = bernoulli_polynomial(t, 6)
        sage: p.parent()
        Univariate Polynomial Ring in t over Rational Field

    We verify an instance of the formula which is the origin of
    the Bernoulli polynomials (and numbers)::

        sage: power_sum = sum(k^4 for k in range(10))
        sage: 5*power_sum == bernoulli_polynomial(10, 5) - bernoulli(5)                 # needs sage.libs.flint
        True

    TESTS::

        sage: x = polygen(QQ, 'x')
        sage: bernoulli_polynomial(x, 0).parent()
        Univariate Polynomial Ring in x over Rational Field

    REFERENCES:

    - :wikipedia:`Bernoulli_polynomials`
    """
    try:
        n = ZZ(n)
        if n < 0:
            raise TypeError
    except TypeError:
        raise ValueError("the second argument must be a nonnegative integer")

    if n == 0:
        return x**0   # result should be in the parent of x

    if n == 1:
        return x - ZZ.one() / 2

    k = n.mod(2)
    coeffs = [0] * k + sum(([n.binomial(i) * bernoulli(n - i), 0]
                            for i in range(k, n + 1, 2)), [])
    coeffs[-3] = -n / 2

    if isinstance(x, Polynomial):
        try:
            return x.parent()(coeffs)(x)
        except TypeError:
            pass

    x2 = x * x
    xi = x**k
    s = 0
    for i in range(k, n - 1, 2):
        s += coeffs[i] * xi
        t = xi
        xi *= x2
    s += xi - t * x * n / 2
    return s
