r"""
Combinatorial Functions

This module implements some combinatorial functions, as listed
below. For a more detailed description, see the relevant
docstrings.

**Sequences:**


-  Bell numbers, :func:`bell_number`

-  Catalan numbers, :func:`catalan_number` (not to be
   confused with the Catalan constant)

-  Narayana numbers, :func:`narayana_number`

-  Euler numbers, :func:`euler_number` (Maxima)

-  Eulerian numbers, :func:`eulerian_number`

-  Eulerian polynomial, :func:`eulerian_polynomial`

-  Fibonacci numbers, :func:`fibonacci` (PARI) and
   :func:`fibonacci_number` (GAP) The PARI version is
   better.

-  Lucas numbers, :func:`lucas_number1`, :func:`lucas_number2`.

-  Stirling numbers, :func:`stirling_number1`,
   :func:`stirling_number2`.

-  Polygonal numbers, :func:`polygonal_number`

**Set-theoretic constructions:**

-  Derangements of a multiset, :func:`derangements` and
   :func:`number_of_derangements`.

-  Tuples of a multiset, :func:`tuples` and
   :func:`number_of_tuples`. An ordered tuple of length k of
   set S is a ordered selection with repetitions of S and is
   represented by a sorted list of length k containing elements from
   S.

-  Unordered tuples of a set, :func:`unordered_tuples` and
   :func:`number_of_unordered_tuples`. An unordered tuple
   of length k of set S is an unordered selection with repetitions of S
   and is represented by a sorted list of length k containing elements
   from S.

.. WARNING::

   The following function is deprecated and will soon be removed.

    - Permutations of a multiset, :func:`permutations`,
      :func:`permutations_iterator`, :func:`number_of_permutations`. A
      permutation is a list that contains exactly the same elements but possibly
      in different order.

**Related functions:**

-  Bernoulli polynomials, :func:`bernoulli_polynomial`

**Implemented in other modules (listed for completeness):**

The ``sage.arith.all`` module contains the following
combinatorial functions:

-  binomial the binomial coefficient (wrapped from PARI)

-  factorial (wrapped from PARI)

-  partition (from the Python Cookbook) Generator of the list of
   all the partitions of the integer `n`.

-  :func:`number_of_partitions` (wrapped from PARI) the
   *number* of partitions:

-  :func:`falling_factorial` Definition: for integer
   `a \ge 0` we have `x(x-1) \cdots (x-a+1)`. In all
   other cases we use the GAMMA-function:
   `\frac {\Gamma(x+1)} {\Gamma(x-a+1)}`.

-  :func:`rising_factorial` Definition: for integer
   `a \ge 0` we have `x(x+1) \cdots (x+a-1)`. In all
   other cases we use the GAMMA-function:
   `\frac {\Gamma(x+a)} {\Gamma(x)}`.

-  gaussian_binomial the gaussian binomial

.. MATH::

    \binom{n}{k}_q = \frac{(1-q^m)(1-q^{m-1})\cdots (1-q^{m-r+1})}{(1-q)(1-q^2)\cdots (1-q^r)}.

The ``sage.groups.perm_gps.permgroup_elements``
contains the following combinatorial functions:


-  matrix method of PermutationGroupElement yielding the
   permutation matrix of the group element.

.. TODO::

    GUAVA commands:
        * VandermondeMat
        * GrayMat returns a list of all different vectors of length n over
          the field F, using Gray ordering.
    Not in GAP:
        * Rencontres numbers (:wikipedia:`Rencontres_number`)

REFERENCES:

- :wikipedia:`Twelvefold_way` (general reference)

AUTHORS:

- David Joyner (2006-07): initial implementation.

- William Stein (2006-07): editing of docs and code; many
  optimizations, refinements, and bug fixes in corner cases

- David Joyner (2006-09): bug fix for combinations, added
  permutations_iterator, combinations_iterator from Python Cookbook,
  edited docs.

- David Joyner (2007-11): changed permutations, added hadamard_matrix

- Florent Hivert (2009-02): combinatorial class cleanup

- Fredrik Johansson (2010-07): fast implementation of ``stirling_number2``

- Punarbasu Purkayastha (2012-12): deprecate arrangements, combinations,
  combinations_iterator, and clean up very old deprecated methods.

Functions and classes
---------------------
"""

# ****************************************************************************
#       Copyright (C) 2006 David Joyner <wdjoyner@gmail.com>,
#                     2007 Mike Hansen <mhansen@gmail.com>,
#                     2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from __future__ import annotations
from typing import Iterator

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.integer import Integer
from sage.rings.infinity import infinity
from sage.arith.all import bernoulli, factorial
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.libs.all import pari
from sage.misc.prandom import randint
from sage.misc.misc_c import prod
from sage.misc.cachefunc import cached_function
from sage.structure.sage_object import SageObject
from sage.structure.parent import Parent
from sage.misc.lazy_import import lazy_import
from sage.misc.lazy_attribute import lazy_attribute
from .combinat_cython import _stirling_number2
from sage.categories.enumerated_sets import EnumeratedSets
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.structure.element import Element
lazy_import('sage.interfaces.maxima_lib', 'maxima')


def bell_number(n, algorithm='flint', **options) -> Integer:
    r"""
    Return the `n`-th Bell number.

    This is the number of ways to partition a set
    of `n` elements into pairwise disjoint nonempty subsets.

    INPUT:

    - ``n`` -- a positive integer

    - ``algorithm`` -- (Default: ``'flint'``) any one of the following:

      - ``'dobinski'`` -- Use Dobinski's formula implemented in Sage

      - ``'flint'`` -- Wrap FLINT's ``arith_bell_number``

      - ``'gap'`` -- Wrap GAP's ``Bell``

      - ``'mpmath'`` -- Wrap mpmath's ``bell``

    .. WARNING::

        When using the mpmath algorithm to compute Bell numbers and you specify
        ``prec``, it can return incorrect results due to low precision. See
        the examples section.

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

    EXAMPLES::

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

        sage: k = bell_number(30, 'mpmath'); k
        846749014511809332450147
        sage: k == bell_number(30)
        True

    If you knows what precision is necessary before computing the Bell number,
    you can use the ``prec`` option::

        sage: k2 = bell_number(30, 'mpmath', prec=30); k2
        846749014511809332450147
        sage: k == k2
        True

    .. WARNING::

            Running mpmath with the precision set too low can result in
            incorrect results::

                sage: k = bell_number(30, 'mpmath', prec=15); k
                846749014511809388871680
                sage: k == bell_number(30)
                False

    TESTS::

        sage: all(bell_number(n) == bell_number(n,'dobinski') for n in range(200))
        True
        sage: all(bell_number(n) == bell_number(n,'gap') for n in range(200))
        True
        sage: all(bell_number(n) == bell_number(n,'mpmath', prec=500) for n in range(200, 220))
        True

    AUTHORS:

    - Robert Gerbicz

    - Jeroen Demeyer: improved implementation of Dobinski formula with
      more accurate error estimates (:trac:`17157`)

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
        import sage.libs.flint.arith
        return sage.libs.flint.arith.bell_number(n)

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
        from sage.rings.all import RealField
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

    OUTPUT:

    integer

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

    - ``n`` -- an integer

    - ``k`` -- an integer between ``0`` and ``n - 1``

    OUTPUT:

    an integer

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

    - ``n`` -- a positive integer

    - ``algorithm`` -- (Default: ``'flint'``) any one of the following:

      - ``'maxima'`` -- Wraps Maxima's ``euler``.

      - ``'flint'`` -- Wrap FLINT's ``arith_euler_number``

    EXAMPLES::

        sage: [euler_number(i) for i in range(10)]
        [1, 0, -1, 0, 5, 0, -61, 0, 1385, 0]
        sage: x = PowerSeriesRing(QQ, 'x').gen().O(10)
        sage: 2/(exp(x)+exp(-x))
        1 - 1/2*x^2 + 5/24*x^4 - 61/720*x^6 + 277/8064*x^8 + O(x^10)
        sage: [euler_number(i)/factorial(i) for i in range(11)]
        [1, 0, -1/2, 0, 5/24, 0, -61/720, 0, 277/8064, 0, -50521/3628800]
        sage: euler_number(-1)
        Traceback (most recent call last):
        ...
        ValueError: n (=-1) must be a nonnegative integer

    TESTS::

        sage: euler_number(6, 'maxima')
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
        import sage.libs.flint.arith
        return sage.libs.flint.arith.euler_number(n)
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

    - ``algorithm`` -- ``"recursive"`` (default) or ``"formula"``

    OUTPUT:

    an integer

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

    - ``n`` -- an integer

    - ``algorithm`` -- ``"derivative"`` (default) or ``"coeffs"``

    OUTPUT:

    polynomial in one variable ``t``

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


def fibonacci(n, algorithm="pari") -> Integer:
    """
    Return the `n`-th Fibonacci number.

    The Fibonacci sequence `F_n` is defined by the initial
    conditions `F_1 = F_2 = 1` and the recurrence relation
    `F_{n+2} = F_{n+1} + F_n`. For negative `n` we
    define `F_n = (-1)^{n+1}F_{-n}`, which is consistent with
    the recurrence relation.

    INPUT:

    - ``algorithm`` -- a string:

      * ``"pari"`` - (default) use the PARI C library's
        :pari:`fibo` function

      * ``"gap"`` - use GAP's Fibonacci function

    .. NOTE::

        PARI is tens to hundreds of times faster than GAP here.
        Moreover, PARI works for every large input whereas GAP does not.

    EXAMPLES::

        sage: fibonacci(10)
        55
        sage: fibonacci(10, algorithm='gap')
        55

    ::

        sage: fibonacci(-100)
        -354224848179261915075
        sage: fibonacci(100)
        354224848179261915075

    ::

        sage: fibonacci(0)
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

    -  ``n`` -- integer

    -  ``P, Q`` -- integer or rational numbers

    OUTPUT: integer or rational number

    EXAMPLES::

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

        sage: lucas = lambda n : Integer((5/2)*lucas_number1(n,1,-1)+(1/2)*lucas_number2(n,1,-1))
        sage: [[lucas(n),is_prime(lucas(n)),n+1,is_prime(n+1)] for n in range(15)]
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


    -  ``n`` - integer

    -  ``P, Q`` - integer or rational numbers


    OUTPUT: integer or rational number

    EXAMPLES::

        sage: [lucas_number2(i,1,-1) for i in range(10)]
        [2, 1, 3, 4, 7, 11, 18, 29, 47, 76]
        sage: [fibonacci(i-1)+fibonacci(i+1) for i in range(10)]
        [2, 1, 3, 4, 7, 11, 18, 29, 47, 76]

    ::

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

        sage: [lucas_number2(n,1,-1) for n in range(10)]
        [2, 1, 3, 4, 7, 11, 18, 29, 47, 76]
    """
    n = ZZ(n)
    P = QQ(P)
    Q = QQ(Q)
    from sage.libs.gap.libgap import libgap
    return libgap.Lucas(P, Q, n)[1].sage()


def stirling_number1(n, k) -> Integer:
    r"""
    Return the `n`-th Stirling number `S_1(n,k)` of the first kind.

    This is the number of permutations of `n` points with `k` cycles.

    This wraps GAP's ``Stirling1``.

    See :wikipedia:`Stirling_numbers_of_the_first_kind`.

    EXAMPLES::

        sage: stirling_number1(3,2)
        3
        sage: stirling_number1(5,2)
        50
        sage: 9*stirling_number1(9,5)+stirling_number1(9,4)
        269325
        sage: stirling_number1(10,5)
        269325

    Indeed, `S_1(n,k) = S_1(n-1,k-1) + (n-1)S_1(n-1,k)`.
    """
    n = ZZ(n)
    k = ZZ(k)
    from sage.libs.gap.libgap import libgap
    return libgap.Stirling1(n, k).sage()


def stirling_number2(n, k, algorithm=None) -> Integer:
    r"""
    Return the `n`-th Stirling number `S_2(n,k)` of the second
    kind.

    This is the number of ways to partition a set of `n` elements into `k`
    pairwise disjoint nonempty subsets. The `n`-th Bell number is the
    sum of the `S_2(n,k)`'s, `k=0,...,n`.

    INPUT:

       *  ``n`` - nonnegative machine-size integer
       *  ``k`` - nonnegative machine-size integer
       * ``algorithm``:

         * None (default) - use native implementation
         * ``"maxima"`` - use Maxima's stirling2 function
         * ``"gap"`` - use GAP's Stirling2 function

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
        sage: n = stirling_number2(20,11)
        sage: n
        1900842429486
        sage: type(n)
        <class 'sage.rings.integer.Integer'>
        sage: n = stirling_number2(20,11,algorithm='gap')
        sage: n
        1900842429486
        sage: type(n)
        <class 'sage.rings.integer.Integer'>
        sage: n = stirling_number2(20,11,algorithm='maxima')
        sage: n
        1900842429486
        sage: type(n)
        <class 'sage.rings.integer.Integer'>

     Sage's implementation splitting the computation of the Stirling
     numbers of the second kind in two cases according to `n`, let us
     check the result it gives agree with both maxima and gap.

     For `n<200`::

         sage: for n in Subsets(range(100,200), 5).random_element():
         ....:     for k in Subsets(range(n), 5).random_element():
         ....:         s_sage = stirling_number2(n,k)
         ....:         s_maxima = stirling_number2(n,k, algorithm = "maxima")
         ....:         s_gap = stirling_number2(n,k, algorithm = "gap")
         ....:         if not (s_sage == s_maxima and s_sage == s_gap):
         ....:             print("Error with n<200")

     For `n\geq 200`::

         sage: for n in Subsets(range(200,300), 5).random_element():
         ....:     for k in Subsets(range(n), 5).random_element():
         ....:         s_sage = stirling_number2(n,k)
         ....:         s_maxima = stirling_number2(n,k, algorithm = "maxima")
         ....:         s_gap = stirling_number2(n,k, algorithm = "gap")
         ....:         if not (s_sage == s_maxima and s_sage == s_gap):
         ....:             print("Error with n<200")


     TESTS:

     Checking an exception is raised whenever a wrong value is given
     for ``algorithm``::

         sage: s_sage = stirling_number2(50,3, algorithm = "CloudReading")
         Traceback (most recent call last):
         ...
         ValueError: unknown algorithm: CloudReading
    """
    n = ZZ(n)
    k = ZZ(k)
    if algorithm is None:
        return _stirling_number2(n, k)
    elif algorithm == 'gap':
        from sage.libs.gap.libgap import libgap
        return libgap.Stirling2(n, k).sage()
    elif algorithm == 'maxima':
        return ZZ(maxima.stirling2(n, k))  # type:ignore
    else:
        raise ValueError("unknown algorithm: %s" % algorithm)


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

    OUTPUT: an integer

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

        sage: polygonal_number(3.5, 1)
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

        -  ``l`` -- a list or any object that can be converted to a
           list by calling ``list()``.

        - ``copy`` -- (boolean, default ``True``) if ``False``, then
          ``l`` must be a ``list``, which is assigned to ``self._list``
          without copying.

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

        Test indirectly that we copy the input (see :trac:`18184`)::

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

        Check that :trac:`14065` is fixed::

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
        Return ``True`` if ``self`` is non-zero.

        We consider a list to be zero if it has length zero.

        TESTS::

            sage: c = CombinatorialObject([1,2,3])
            sage: not c
            False
            sage: c = CombinatorialObject([])
            sage: not c
            True

        Check that :trac:`14065` is fixed::

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

    __nonzero__ = __bool__

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

    -  ``parent`` -- the :class:`Parent` class for this element.

    -  ``lst`` -- a list or any object that can be converted to a
       list by calling ``list()``.

    EXAMPLES::

        sage: from sage.combinat.combinat import CombinatorialElement
        sage: e = CombinatorialElement(Partitions(6), [3,2,1])
        sage: e == loads(dumps(e))
        True
        sage: parent(e)
        Partitions of the integer 6
        sage: list(e)
        [3, 2, 1]

    Check classcalls::

        sage: class Foo(CombinatorialElement):
        ....:     @staticmethod
        ....:     def __classcall__(cls, x):
        ....:         return x
        sage: Foo(17)
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
        super(CombinatorialElement, self).__init__(L)
        super(CombinatorialObject, self).__init__(parent)


class CombinatorialClass(Parent, metaclass=ClasscallMetaclass):
    """
    This class is deprecated, and will disappear as soon as all derived
    classes in Sage's library will have been fixed. Please derive
    directly from Parent and use the category :class:`EnumeratedSets`,
    :class:`FiniteEnumeratedSets`, or :class:`InfiniteEnumeratedSets`, as
    appropriate.

    For examples, see::

        sage: FiniteEnumeratedSets().example()
        An example of a finite enumerated set: {1,2,3}
        sage: InfiniteEnumeratedSets().example()
        An example of an infinite enumerated set: the non negative integers
    """
    def __init__(self, category=None):
        """
        TESTS::

            sage: C = sage.combinat.combinat.CombinatorialClass()
            sage: C.category()
            Category of enumerated sets
            sage: C.__class__
            <class 'sage.combinat.combinat.CombinatorialClass_with_category'>
            sage: isinstance(C, Parent)
            True
            sage: C = sage.combinat.combinat.CombinatorialClass(category = FiniteEnumeratedSets())
            sage: C.category()
            Category of finite enumerated sets
        """
        Parent.__init__(self, category=EnumeratedSets().or_subcategory(category))

    def is_finite(self) -> bool:
        """
        Return whether ``self`` is finite or not.

        EXAMPLES::

            sage: Partitions(5).is_finite()
            True
            sage: Permutations().is_finite()
            False
        """
        return self.cardinality() != infinity

    def __getitem__(self, i):
        """
        Return the combinatorial object of rank i.

        EXAMPLES::

            sage: class C(CombinatorialClass):
            ....:     def __iter__(self):
            ....:         return iter([1,2,3])
            sage: c = C()
            sage: c[0]
            1
            sage: c[2]
            3
            sage: c[4]
            Traceback (most recent call last):
            ...
            ValueError: the value must be between 0 and 2 inclusive
        """
        return self.unrank(i)

    def __str__(self) -> str:
        """
        Return a string representation of self.

        EXAMPLES::

            sage: str(Partitions(5))
            'Partitions of the integer 5'
        """
        return repr(self)

    def _repr_(self) -> str:
        """
        EXAMPLES::

            sage: repr(Partitions(5))   # indirect doctest
            'Partitions of the integer 5'
        """
        if hasattr(self, '_name') and self._name:
            return self._name
        else:
            return "Combinatorial Class -- REDEFINE ME!"

    def __contains__(self, x) -> bool:
        """
        Test whether or not the combinatorial class contains the object x.
        This raises a NotImplementedError as a default since _all_
        subclasses of CombinatorialClass should override this.

        Note that we could replace this with a default implementation that
        just iterates through the elements of the combinatorial class and
        checks for equality. However, since we use __contains__ for
        type checking, this operation should be cheap and should be
        implemented manually for each combinatorial class.

        EXAMPLES::

            sage: C = CombinatorialClass()
            sage: x in C
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def __eq__(self, other):
        """
        Compare two different combinatorial classes.

        For now, the comparison is done just on their repr's.

        EXAMPLES::

            sage: p5 = Partitions(5)
            sage: p6 = Partitions(6)
            sage: repr(p5) == repr(p6)
            False
            sage: p5 == p6
            False
        """
        return repr(self) == repr(other)

    def __ne__(self, other):
        """
        Test unequality of ``self`` and ``other``.

        EXAMPLES::

            sage: p5 = Partitions(5)
            sage: p6 = Partitions(6)
            sage: p5 != p6
            True
        """
        return not (self == other)

    def __hash__(self):
        """
        Create a hash value. This is based on the string representation.

        Note that in Python 3 objects that define __eq__ do not inherit their __hash__
        function. Without an explicit __hash__ they are no longer hashable.

        TESTS::

            sage: C = CombinatorialClass()
            sage: hash(C) == hash(repr(C))
            True
        """
        return hash(repr(self))

    def __cardinality_from_iterator(self) -> Integer | infinity:
        """
        Default implementation of cardinality which just goes through the iterator
        of the combinatorial class to count the number of objects.

        EXAMPLES::

            sage: class C(CombinatorialClass):
            ....:     def __iter__(self):
            ....:         return iter([1,2,3])
            sage: C().cardinality() #indirect doctest
            3
        """
        c = Integer(0)
        one = Integer(1)
        for _ in self:
            c += one
        return c

    cardinality = __cardinality_from_iterator

    # __call__, element_class, and _element_constructor_ are poor
    # man's versions of those from Parent. This is for transition,
    # until all combinatorial classes are proper parents (in Parent)
    # and use coercion, etcc

    def __call__(self, x):
        """
        Return x as an element of the combinatorial class's object class.

        EXAMPLES::

            sage: p5 = Partitions(5)
            sage: a = [2,2,1]
            sage: type(a)
            <class 'list'>
            sage: a = p5(a)
            sage: type(a)
            <class 'sage.combinat.partition.Partitions_n_with_category.element_class'>
            sage: p5([2,1])
            Traceback (most recent call last):
            ...
            ValueError: [2, 1] is not an element of Partitions of the integer 5
        """
        if x in self:
            return self._element_constructor_(x)
        else:
            raise ValueError("%s not in %s" % (x, self))

    Element = CombinatorialObject  # mostly for backward compatibility

    @lazy_attribute
    def element_class(self):
        """
        This function is a temporary helper so that a CombinatorialClass
        behaves as a parent for creating elements. This will disappear when
        combinatorial classes will be turned into actual parents (in the
        category EnumeratedSets).

        TESTS::

            sage: P5 = Partitions(5)
            sage: P5.element_class
            <class 'sage.combinat.partition.Partitions_n_with_category.element_class'>
        """
        # assert not isinstance(self, Parent) # Raises an alert if we override the proper definition from Parent
        return self.Element

    def _element_constructor_(self, x):
        """
        This function is a temporary helper so that a CombinatorialClass
        behaves as a parent for creating elements. This will disappear when
        combinatorial classes will be turned into actual parents (in the
        category EnumeratedSets).

        TESTS::

            sage: P5 = Partitions(5)
            sage: p = P5([3,2])      # indirect doctest
            sage: type(p)
            <class 'sage.combinat.partition.Partitions_n_with_category.element_class'>
        """
        # assert not isinstance(self, Parent) # Raises an alert if we override the proper definition from Parent
        return self.element_class(x)

    def __list_from_iterator(self):
        """
        The default implementation of list which builds the list from the
        iterator.

        EXAMPLES::

            sage: class C(CombinatorialClass):
            ....:     def __iter__(self):
            ....:         return iter([1,2,3])
            sage: C().list() #indirect doctest
            [1, 2, 3]
        """
        return [x for x in self]

    # Set list to the default implementation
    list = __list_from_iterator

    # Set the default object class to be CombinatorialObject
    Element = CombinatorialObject

    def __iterator_from_next(self) -> Iterator:
        """
        An iterator to use when the .first() and .next(x) methods are provided.

        EXAMPLES::

            sage: C = CombinatorialClass()
            sage: C.first = lambda: 0
            sage: C.next  = lambda c: c+1
            sage: it = iter(C) # indirect doctest
            sage: [next(it) for _ in range(4)]
            [0, 1, 2, 3]
        """
        f = self.first()
        yield f
        while True:
            try:
                f = self.next(f)
            except (TypeError, ValueError):
                break

            if f is None or f is False:
                break
            else:
                yield f

    def __iterator_from_previous(self):
        """
        An iterator to use when .last() and .previous() are provided. Note
        that this requires the combinatorial class to be finite. It is not
        recommended to implement combinatorial classes using last and
        previous.

        EXAMPLES::

            sage: C = CombinatorialClass()
            sage: C.last = lambda: 4
            sage: def prev(c):
            ....:     if c <= 1:
            ....:         return None
            ....:     else:
            ....:         return c-1
            sage: C.previous  = prev
            sage: it = iter(C) # indirect doctest
            sage: [next(it) for _ in range(4)]
            [1, 2, 3, 4]
        """
        l = self.last()
        li = [l]
        while True:
            try:
                l = self.previous(l)
            except (TypeError, ValueError):
                break

            if l is None:
                break
            else:
                li.append(l)
        return reversed(li)

    def __iterator_from_unrank(self) -> Iterator:
        """
        An iterator to use when .unrank() is provided.

        EXAMPLES::

            sage: C = CombinatorialClass()
            sage: l = [1,2,3]
            sage: C.unrank = lambda c: l[c]
            sage: list(C) # indirect doctest
            [1, 2, 3]
        """
        r = 0
        u = self.unrank(r)
        yield u
        while True:
            r += 1
            try:
                u = self.unrank(r)
            except (TypeError, ValueError, IndexError):
                break

            if u is None:
                break
            else:
                yield u

    def __iterator_from_list(self) -> Iterator:
        """
        An iterator to use when .list() is provided()

        EXAMPLES::

            sage: C = CombinatorialClass()
            sage: C.list = lambda: [1, 2, 3]
            sage: list(C) # indirect doctest
            [1, 2, 3]
        """
        for x in self.list():
            yield x

    def __iter__(self):
        """
        Allows the combinatorial class to be treated as an iterator. Default
        implementation.

        EXAMPLES::

            sage: p5 = Partitions(5)
            sage: [i for i in p5]
            [[5], [4, 1], [3, 2], [3, 1, 1], [2, 2, 1], [2, 1, 1, 1], [1, 1, 1, 1, 1]]
            sage: C = CombinatorialClass()
            sage: iter(C)
            Traceback (most recent call last):
            ...
            NotImplementedError: iterator called but not implemented
        """
        # Check whether .first() and .next(x) are overridden in the subclass
        if (self.first != self.__first_from_iterator and
                self.next != self.__next_from_iterator):
            return self.__iterator_from_next()
        # Check whether .last() and .previous() are overridden in the subclass
        elif (self.last != self.__last_from_iterator and
              self.previous != self.__previous_from_iterator):
            return self.__iterator_from_previous()
        # Check whether .unrank() is overridden in the subclass
        elif self.unrank != self.__unrank_from_iterator:
            return self.__iterator_from_unrank()
        # Check whether .list() is overridden in the subclass
        elif self.list != self.__list_from_iterator:
            return self.__iterator_from_list()
        else:
            raise NotImplementedError("iterator called but not implemented")

    def __unrank_from_iterator(self, r):
        """
        Default implementation of unrank which goes through the iterator.

        EXAMPLES::

            sage: C = CombinatorialClass()
            sage: C.list = lambda: [1,2,3]
            sage: C.unrank(1) # indirect doctest
            2
        """
        counter = 0
        for u in self:
            if counter == r:
                return u
            counter += 1
        raise ValueError("the value must be between %s and %s inclusive" % (0, counter - 1))

    # Set the default implementation of unrank
    unrank = __unrank_from_iterator

    def __random_element_from_unrank(self):
        """
        Default implementation of random which uses unrank.

        EXAMPLES::

            sage: C = CombinatorialClass()
            sage: C.list = lambda: [1,2,3]
            sage: C.random_element()       # random  # indirect doctest
            1
        """
        c = self.cardinality()
        r = randint(0, c - 1)
        return self.unrank(r)

    # Set the default implementation of random
    random_element = __random_element_from_unrank

    def __rank_from_iterator(self, obj):
        """
        Default implementation of rank which uses iterator.

        EXAMPLES::

            sage: C = CombinatorialClass()
            sage: C.list = lambda: [1,2,3]
            sage: C.rank(3) # indirect doctest
            2
        """
        r = 0
        for i in self:
            if i == obj:
                return r
            r += 1
        raise ValueError

    rank = __rank_from_iterator

    def __first_from_iterator(self):
        """
        Default implementation for first which uses iterator.

        EXAMPLES::

            sage: C = CombinatorialClass()
            sage: C.list = lambda: [1,2,3]
            sage: C.first() # indirect doctest
            1
        """
        for i in self:
            return i

    first = __first_from_iterator

    def __last_from_iterator(self):
        """
        Default implementation for first which uses iterator.

        EXAMPLES::

            sage: C = CombinatorialClass()
            sage: C.list = lambda: [1,2,3]
            sage: C.last() # indirect doctest
            3
        """
        for i in self:
            pass
        return i

    last = __last_from_iterator

    def __next_from_iterator(self, obj):
        """
        Default implementation for next which uses iterator.

        EXAMPLES::

            sage: C = CombinatorialClass()
            sage: C.list = lambda: [1,2,3]
            sage: C.next(2) # indirect doctest
            3
        """
        found = False
        for i in self:
            if found:
                return i
            if i == obj:
                found = True
        return None

    next = __next_from_iterator

    def __previous_from_iterator(self, obj):
        """
        Default implementation for next which uses iterator.

        EXAMPLES::

            sage: C = CombinatorialClass()
            sage: C.list = lambda: [1,2,3]
            sage: C.previous(2) # indirect doctest
            1
        """
        prev = None
        for i in self:
            if i == obj:
                break
            prev = i
        return prev

    previous = __previous_from_iterator

    def filter(self, f, name=None):
        """
        Return the combinatorial subclass of f which consists of the
        elements x of ``self`` such that f(x) is ``True``.

        EXAMPLES::

            sage: from sage.combinat.combinat import Permutations_CC
            sage: P = Permutations_CC(3).filter(lambda x: x.avoids([1,2]))
            sage: P.list()
            [[3, 2, 1]]
        """
        return FilteredCombinatorialClass(self, f, name=name)

    def union(self, right_cc, name=None):
        """
        Return the combinatorial class representing the union of ``self`` and
        ``right_cc``.

        EXAMPLES::

            sage: from sage.combinat.combinat import Permutations_CC
            sage: P = Permutations_CC(2).union(Permutations_CC(1))
            sage: P.list()
            [[1, 2], [2, 1], [1]]
        """
        if not isinstance(right_cc, CombinatorialClass):
            raise TypeError("right_cc must be a CombinatorialClass")
        return UnionCombinatorialClass(self, right_cc, name=name)

    def map(self, f, name=None):
        r"""
        Return the image `\{f(x) | x \in \text{self}\}` of this combinatorial
        class by `f`, as a combinatorial class.

        `f` is supposed to be injective.

        EXAMPLES::

            sage: R = Permutations(3).map(attrcall('reduced_word')); R
            Image of Standard permutations of 3 by *.reduced_word()
            sage: R.cardinality()
            6
            sage: R.list()
            [[], [2], [1], [1, 2], [2, 1], [2, 1, 2]]
            sage: [ r for r in R]
            [[], [2], [1], [1, 2], [2, 1], [2, 1, 2]]

        If the function is not injective, then there may be repeated elements::

            sage: P = Partitions(4)
            sage: P.list()
            [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]
            sage: P.map(len).list()
            [1, 2, 2, 3, 4]

        TESTS::

            sage: R = Permutations(3).map(attrcall('reduced_word'))
            sage: R == loads(dumps(R))
            True
        """
        return MapCombinatorialClass(self, f, name)


class FilteredCombinatorialClass(CombinatorialClass):
    def __init__(self, combinatorial_class, f, name=None):
        """
        A filtered combinatorial class F is a subset of another
        combinatorial class C specified by a function f that takes in an
        element c of C and returns True if and only if c is in F.

        TESTS::

            sage: from sage.combinat.combinat import Permutations_CC
            sage: Permutations_CC(3).filter(lambda x: x.avoids([1,2]))
            Filtered subclass of Standard permutations of 3
        """
        self.f = f
        self.combinatorial_class = combinatorial_class
        self._name = name

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.combinat.combinat import Permutations_CC
            sage: P = Permutations_CC(3).filter(lambda x: x.avoids([1,2]))
            sage: P.__repr__()
            'Filtered subclass of Standard permutations of 3'
            sage: P._name = 'Permutations avoiding [1, 2]'
            sage: P.__repr__()
            'Permutations avoiding [1, 2]'
        """
        if self._name:
            return self._name
        else:
            return "Filtered subclass of " + repr(self.combinatorial_class)

    def __contains__(self, x) -> bool:
        """
        EXAMPLES::

            sage: from sage.combinat.combinat import Permutations_CC
            sage: P = Permutations_CC(3).filter(lambda x: x.avoids([1,2]))
            sage: 'cat' in P
            False
            sage: [4,3,2,1] in P
            False
            sage: Permutation([1,2,3]) in P
            False
            sage: Permutation([3,2,1]) in P
            True
        """
        return x in self.combinatorial_class and self.f(x)

    def cardinality(self) -> Integer:
        """
        EXAMPLES::

            sage: from sage.combinat.combinat import Permutations_CC
            sage: P = Permutations_CC(3).filter(lambda x: x.avoids([1,2]))
            sage: P.cardinality()
            1
        """
        c = 0
        for _ in self:
            c += 1
        return c

    def __iter__(self) -> Iterator:
        """
        EXAMPLES::

            sage: from sage.combinat.combinat import Permutations_CC
            sage: P = Permutations_CC(3).filter(lambda x: x.avoids([1,2]))
            sage: list(P)
            [[3, 2, 1]]
        """
        for x in self.combinatorial_class:
            if self.f(x):
                yield x


class UnionCombinatorialClass(CombinatorialClass):
    def __init__(self, left_cc, right_cc, name=None):
        """
        A UnionCombinatorialClass is a union of two other combinatorial
        classes.

        TESTS::

            sage: from sage.combinat.combinat import Permutations_CC
            sage: P = Permutations_CC(3).union(Permutations_CC(2))
            sage: P == loads(dumps(P))
            True
        """
        self.left_cc = left_cc
        self.right_cc = right_cc
        self._name = name

    def __repr__(self) -> str:
        """
        TESTS::

            sage: from sage.combinat.combinat import Permutations_CC
            sage: print(repr(Permutations_CC(3).union(Permutations_CC(2))))
            Union combinatorial class of
                Standard permutations of 3
            and
                Standard permutations of 2
        """
        if self._name:
            return self._name
        else:
            return "Union combinatorial class of \n    %s\nand\n    %s" % (self.left_cc, self.right_cc)

    def __contains__(self, x) -> bool:
        """
        EXAMPLES::

            sage: from sage.combinat.combinat import Permutations_CC
            sage: P = Permutations_CC(3).union(Permutations_CC(2))
            sage: [1,2] in P
            True
            sage: [3,2,1] in P
            True
            sage: [1,2,3,4] in P
            False
        """
        return x in self.left_cc or x in self.right_cc

    def cardinality(self) -> Integer | infinity:
        """
        EXAMPLES::

            sage: from sage.combinat.combinat import Permutations_CC
            sage: P = Permutations_CC(3).union(Permutations_CC(2))
            sage: P.cardinality()
            8
        """
        return self.left_cc.cardinality() + self.right_cc.cardinality()

    def list(self):
        """
        EXAMPLES::

            sage: from sage.combinat.combinat import Permutations_CC
            sage: P = Permutations_CC(3).union(Permutations_CC(2))
            sage: P.list()
            [[1, 2, 3],
             [1, 3, 2],
             [2, 1, 3],
             [2, 3, 1],
             [3, 1, 2],
             [3, 2, 1],
             [1, 2],
             [2, 1]]
        """
        return self.left_cc.list() + self.right_cc.list()

    def __iter__(self) -> Iterator:
        """
        EXAMPLES::

            sage: from sage.combinat.combinat import Permutations_CC
            sage: P = Permutations_CC(3).union(Permutations_CC(2))
            sage: list(P)
            [[1, 2, 3],
             [1, 3, 2],
             [2, 1, 3],
             [2, 3, 1],
             [3, 1, 2],
             [3, 2, 1],
             [1, 2],
             [2, 1]]
        """
        for x in self.left_cc:
            yield x
        for x in self.right_cc:
            yield x

    def first(self):
        """
        EXAMPLES::

            sage: from sage.combinat.combinat import Permutations_CC
            sage: P = Permutations_CC(3).union(Permutations_CC(2))
            sage: P.first()
            [1, 2, 3]
        """
        return self.left_cc.first()

    def last(self):
        """
        EXAMPLES::

            sage: from sage.combinat.combinat import Permutations_CC
            sage: P = Permutations_CC(3).union(Permutations_CC(2))
            sage: P.last()
            [2, 1]
        """
        return self.right_cc.last()

    def rank(self, x):
        """
        EXAMPLES::

            sage: from sage.combinat.combinat import Permutations_CC
            sage: P = Permutations_CC(3).union(Permutations_CC(2))
            sage: P.rank(Permutation([2,1]))
            7
            sage: P.rank(Permutation([1,2,3]))
            0
        """
        try:
            return self.left_cc.rank(x)
        except (TypeError, ValueError):
            return self.left_cc.cardinality() + self.right_cc.rank(x)

    def unrank(self, x):
        """
        EXAMPLES::

            sage: from sage.combinat.combinat import Permutations_CC
            sage: P = Permutations_CC(3).union(Permutations_CC(2))
            sage: P.unrank(7)
            [2, 1]
            sage: P.unrank(0)
            [1, 2, 3]
        """
        try:
            return self.left_cc.unrank(x)
        except (TypeError, ValueError):
            return self.right_cc.unrank(x - self.left_cc.cardinality())


class Permutations_CC(CombinatorialClass):
    """
    A testing class for :class:`CombinatorialClass` since :class:`Permutations`
    no longer inherits from :class:`CombinatorialClass` in :trac:`14772`.
    """
    def __init__(self, n):
        """
        EXAMPLES::

            sage: from sage.combinat.combinat import Permutations_CC
            sage: P = Permutations_CC(4)
            sage: loads(dumps(P)) == P
            True
        """
        from sage.combinat.permutation import StandardPermutations_n
        self._permutations = StandardPermutations_n(n)

    def __repr__(self) -> str:
        """
        EXAMPLES::

            sage: from sage.combinat.combinat import Permutations_CC
            sage: Permutations_CC(3)
            Standard permutations of 3
        """
        return repr(self._permutations)

    def __contains__(self, x) -> bool:
        """
        EXAMPLES::

            sage: from sage.combinat.combinat import Permutations_CC
            sage: P = Permutations_CC(3)
            sage: [1, 3, 2] in P
            True
        """
        return x in self._permutations

    def __iter__(self):
        """
        EXAMPLES::

            sage: from sage.combinat.combinat import Permutations_CC
            sage: P = Permutations_CC(3)
            sage: P.list()
            [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]
        """
        return iter(self._permutations)


##############################################################################
class MapCombinatorialClass(CombinatorialClass):
    r"""
    A MapCombinatorialClass models the image of a combinatorial
    class through a function which is assumed to be injective

    See CombinatorialClass.map for examples
    """
    def __init__(self, cc, f, name=None):
        """
        TESTS::

            sage: Partitions(3).map(attrcall('conjugate'))
            Image of Partitions of the integer 3 by *.conjugate()
        """
        self.cc = cc
        self.f = f
        self._name = name

    def __repr__(self) -> str:
        """
        TESTS::

            sage: Partitions(3).map(attrcall('conjugate'))
            Image of Partitions of the integer 3 by *.conjugate()

        """
        if self._name:
            return self._name
        else:
            return "Image of %s by %s" % (self.cc, self.f)

    def cardinality(self) -> Integer | infinity:
        """
        Return the cardinality of this combinatorial class

        EXAMPLES::

            sage: R = Permutations(10).map(attrcall('reduced_word'))
            sage: R.cardinality()
            3628800
        """
        return self.cc.cardinality()

    def __iter__(self) -> Iterator:
        """
        Return an iterator over the elements of this combinatorial class

        EXAMPLES::

            sage: R = Permutations(10).map(attrcall('reduced_word'))
            sage: R.cardinality()
            3628800
        """
        for x in self.cc:
            yield self.f(x)

    def an_element(self):
        """
        Return an element of this combinatorial class

        EXAMPLES::

            sage: R = SymmetricGroup(10).map(attrcall('reduced_word'))
            sage: R.an_element()
            [9, 8, 7, 6, 5, 4, 3, 2]
        """
        return self.f(self.cc.an_element())


##############################################################################
class InfiniteAbstractCombinatorialClass(CombinatorialClass):
    r"""
    This is an internal class that should not be used directly.  A class which
    inherits from InfiniteAbstractCombinatorialClass inherits the standard
    methods list and count.

    If self._infinite_cclass_slice exists then self.__iter__ returns an
    iterator for self, otherwise raise NotImplementedError. The method
    self._infinite_cclass_slice is supposed to accept any integer as an
    argument and return something which is iterable.
    """
    def cardinality(self) -> Integer | infinity:
        """
        Count the elements of the combinatorial class.

        EXAMPLES::

            sage: R = InfiniteAbstractCombinatorialClass()
            doctest:warning...
            DeprecationWarning: this class is deprecated, do not use
            See https://trac.sagemath.org/31545 for details.

            sage: R.cardinality()
            +Infinity
        """
        return infinity

    def list(self):
        """
        Return an error since ``self`` is an infinite combinatorial class.

        EXAMPLES::

            sage: R = InfiniteAbstractCombinatorialClass()
            sage: R.list()
            Traceback (most recent call last):
            ...
            NotImplementedError: infinite list
        """
        raise NotImplementedError("infinite list")

    def __iter__(self) -> Iterator:
        """
        Return an iterator for the infinite combinatorial class ``self`` if
        possible or raise a NotImplementedError.

        EXAMPLES::

            sage: R = InfiniteAbstractCombinatorialClass()
            sage: next(iter(R))
            Traceback (most recent call last):
            ...
            NotImplementedError

            sage: c = iter(Compositions()) # indirect doctest
            sage: next(c), next(c), next(c), next(c), next(c), next(c)
            ([], [1], [1, 1], [2], [1, 1, 1], [1, 2])
            sage: next(c), next(c), next(c), next(c), next(c), next(c)
            ([2, 1], [3], [1, 1, 1, 1], [1, 1, 2], [1, 2, 1], [1, 3])
        """
        try:
            finite = self._infinite_cclass_slice
        except AttributeError:
            raise NotImplementedError
        i = 0
        while True:
            for c in finite(i):
                yield c
            i += 1


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

      * ``'itertools'`` - (default) use python's itertools
      * ``'native'`` - use a native Sage implementation

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

        sage: K.<a> = GF(4, 'a')
        sage: mset = [x for x in K if x != 0]
        sage: tuples(mset, 2)
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

    AUTHORS:

    - Jon Hanke (2006-08)
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
    same as ``tuples(S, k, algorithm="native")``.

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

      * ``'naive'`` - (default) use the naive counting `|S|^k`
      * ``'gap'`` - wraps GAP's ``NrTuples``

    .. WARNING::

        When using ``algorithm='gap'``, ``S`` must be a list of objects
        that have string representations that can be interpreted by the GAP
        interpreter. If ``S`` consists of at all complicated Sage
        objects, this function might *not* do what you expect.

    EXAMPLES::

        sage: S = [1,2,3,4,5]
        sage: number_of_tuples(S,2)
        25
        sage: number_of_tuples(S,2, algorithm="gap")
        25
        sage: S = [1,1,2,3,4,5]
        sage: number_of_tuples(S,2)
        25
        sage: number_of_tuples(S,2, algorithm="gap")
        25
        sage: number_of_tuples(S,0)
        1
        sage: number_of_tuples(S,0, algorithm="gap")
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

      * ``'itertools'`` - (default) use python's itertools
      * ``'gap'`` - wraps GAP's ``UnorderedTuples``

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

        sage: unordered_tuples(S, 3, algorithm='gap')
        [(1, 1, 1), (1, 1, 2), (1, 2, 2), (2, 2, 2)]

    We check the result on strings::

        sage: S = ["a","b","c"]
        sage: unordered_tuples(S, 2)
        [('a', 'a'), ('a', 'b'), ('a', 'c'), ('b', 'b'), ('b', 'c'), ('c', 'c')]
        sage: unordered_tuples(S, 2, algorithm='gap')
        [('a', 'a'), ('a', 'b'), ('a', 'c'), ('b', 'b'), ('b', 'c'), ('c', 'c')]

    Lastly we check on a multiset::

        sage: S = [1,1,2]
        sage: unordered_tuples(S, 3) == unordered_tuples(S, 3, 'gap')
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

      * ``'naive'`` - (default) use the naive counting `\binom{|S|+k-1}{k}`
      * ``'gap'`` - wraps GAP's ``NrUnorderedTuples``

    .. WARNING::

        When using ``algorithm='gap'``, ``S`` must be a list of objects
        that have string representations that can be interpreted by the GAP
        interpreter. If ``S`` consists of at all complicated Sage
        objects, this function might *not* do what you expect.

    EXAMPLES::

        sage: S = [1,2,3,4,5]
        sage: number_of_unordered_tuples(S,2)
        15
        sage: number_of_unordered_tuples(S,2, algorithm="gap")
        15
        sage: S = [1,1,2,3,4,5]
        sage: number_of_unordered_tuples(S,2)
        15
        sage: number_of_unordered_tuples(S,2, algorithm="gap")
        15
        sage: number_of_unordered_tuples(S,0)
        1
        sage: number_of_unordered_tuples(S,0, algorithm="gap")
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
    from sage.misc.misc import powerset
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


def bell_polynomial(n: Integer, k: Integer):
    r"""
    Return the Bell Polynomial

    .. MATH::

       B_{n,k}(x_0, x_1, \ldots, x_{n-k}) =
            \sum_{\sum{j_i}=k, \sum{(i+1) j_i}=n}
            \frac{n!}{j_0!j_1!\cdots j_{n-k}!}
            \left(\frac{x_0}{(0+1)!}\right)^{j_0}
            \left(\frac{x_1}{(1+1)!}\right)^{j_1} \cdots
            \left(\frac{x_{n-k}}{(n-k+1)!}\right)^{j_{n-k}}.

    INPUT:

    - ``n`` -- integer

    - ``k`` -- integer

    OUTPUT:

    - a polynomial in `n-k+1` variables over `\ZZ`

    EXAMPLES::

        sage: bell_polynomial(6,2)
        10*x2^2 + 15*x1*x3 + 6*x0*x4
        sage: bell_polynomial(6,3)
        15*x1^3 + 60*x0*x1*x2 + 15*x0^2*x3

    TESTS:

    Check that :trac:`18338` is fixed::

        sage: bell_polynomial(0,0).parent()
        Multivariate Polynomial Ring in x over Integer Ring

        sage: for n in (0..4):
        ....:     print([bell_polynomial(n,k).coefficients() for k in (0..n)])
        [[1]]
        [[], [1]]
        [[], [1], [1]]
        [[], [1], [3], [1]]
        [[], [1], [3, 4], [6], [1]]


    REFERENCES:

    - [Bel1927]_

    AUTHORS:

    - Blair Sutton (2009-01-26)
    - Thierry Monteil (2015-09-29): the result must always be a polynomial.
    """
    from sage.combinat.partition import Partitions
    R = PolynomialRing(ZZ, 'x', n - k + 1)
    vars = R.gens()
    result = R.zero()
    for p in Partitions(n, length=k):  # type:ignore
        factorial_product = 1
        power_factorial_product = 1
        for part, count in p.to_exp_dict().items():
            factorial_product *= factorial(count)
            power_factorial_product *= factorial(part)**count
        coefficient = factorial(n) // (factorial_product * power_factorial_product)
        result += coefficient * prod([vars[i - 1] for i in p])
    return result


def fibonacci_sequence(start, stop=None, algorithm=None) -> Iterator:
    r"""
    Return an iterator over the Fibonacci sequence, for all fibonacci
    numbers `f_n` from ``n = start`` up to (but
    not including) ``n = stop``

    INPUT:

    -  ``start`` -- starting value

    -  ``stop`` -- stopping value

    -  ``algorithm`` -- (default: ``None``) passed on to
       fibonacci function (or not passed on if None, i.e., use the
       default)

    EXAMPLES::

        sage: fibs = [i for i in fibonacci_sequence(10, 20)]
        sage: fibs
        [55, 89, 144, 233, 377, 610, 987, 1597, 2584, 4181]

    ::

        sage: sum([i for i in fibonacci_sequence(100, 110)])
        69919376923075308730013

    .. SEEALSO::

       :func:`fibonacci_xrange`

    AUTHORS:

    - Bobby Moretti
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

        sage: fibs_in_some_range =  [i for i in fibonacci_xrange(10^7, 10^8)]
        sage: len(fibs_in_some_range)
        4
        sage: fibs_in_some_range
        [14930352, 24157817, 39088169, 63245986]

    ::

        sage: fibs = [i for i in fibonacci_xrange(10, 100)]
        sage: fibs
        [13, 21, 34, 55, 89]

    ::

        sage: list(fibonacci_xrange(13, 34))
        [13, 21]

    A solution to the second Project Euler problem::

        sage: sum([i for i in fibonacci_xrange(10^6) if is_even(i)])
        1089154

    .. SEEALSO::

       :func:`fibonacci_sequence`

    AUTHORS:

    - Bobby Moretti
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

    One has `B_n(x) = - n\zeta(1 - n,x)`, where
    `\zeta(s,x)` is the Hurwitz zeta function. Thus, in a
    certain sense, the Hurwitz zeta function generalizes the
    Bernoulli polynomials to non-integer values of n.

    EXAMPLES::

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
        sage: 5*power_sum == bernoulli_polynomial(10, 5) - bernoulli(5)
        True

    REFERENCES:

    - :wikipedia:`Bernoulli_polynomials`
    """
    try:
        n = ZZ(n)
        if n < 0:
            raise TypeError
    except TypeError:
        raise ValueError("The second argument must be a non-negative integer")

    if n == 0:
        return ZZ.one()

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
