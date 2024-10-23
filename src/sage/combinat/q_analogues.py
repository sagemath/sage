r"""
`q`-Analogues
"""

# ****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.combinat.dyck_word import DyckWords
from sage.combinat.partition import _Partitions
from sage.misc.cachefunc import cached_function
from sage.misc.misc_c import prod
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.rings.polynomial.polynomial_ring import polygen
from sage.structure.element import parent


def q_int(n, q=None):
    r"""
    Return the `q`-analogue of the integer `n`.

    The `q`-analogue of the integer `n` is given by

    .. MATH::

        [n]_q = \begin{cases}
        1 + q + \cdots + q^{n-1},  & \text{if } n \geq 0, \\
        -q^{-n} [-n]_q,            & \text{if } n \leq 0.
        \end{cases}

    Consequently, if `q = 1` then `[n]_1 = n` and if `q \neq 1` then
    `[n]_q = (q^n-1)/(q-1)`.

    If the argument `q` is not specified then it defaults to the generator `q`
    of the univariate polynomial ring over the integers.

    EXAMPLES::

        sage: from sage.combinat.q_analogues import q_int
        sage: q_int(3)
        q^2 + q + 1
        sage: q_int(-3)
        -q^-3 - q^-2 - q^-1
        sage: p = ZZ['p'].0
        sage: q_int(3,p)
        p^2 + p + 1
        sage: q_int(3/2)
        Traceback (most recent call last):
        ...
        ValueError: 3/2 must be an integer

    TESTS:

    We check that :issue:`15805` is fixed::

        sage: q_int(0).parent()
        Univariate Polynomial Ring in q over Integer Ring

    We check that :issue:`25715` is fixed::

        sage: q_int(0, 3r)
        0
    """
    if n not in ZZ:
        raise ValueError(f'{n} must be an integer')
    if q is None:
        if n >= 0:
            q = polygen(ZZ, 'q')
        else:
            q = LaurentPolynomialRing(ZZ, 'q').gen()
    if n == 0:
        return parent(q)(0)
    if n > 0:
        return sum(q**i for i in range(n))
    return -q**n * sum(q**i for i in range(-n))


def q_factorial(n, q=None):
    r"""
    Return the `q`-analogue of the factorial `n!`.

    This is the product

    .. MATH::

        [1]_q [2]_q \cdots [n]_q
        = 1 \cdot (1+q) \cdot (1+q+q^2) \cdots (1+q+q^2+\cdots+q^{n-1}) .

    If `q` is unspecified, then this function defaults to
    using the generator `q` for a univariate polynomial
    ring over the integers.

    EXAMPLES::

        sage: from sage.combinat.q_analogues import q_factorial
        sage: q_factorial(3)
        q^3 + 2*q^2 + 2*q + 1
        sage: p = ZZ['p'].0
        sage: q_factorial(3, p)
        p^3 + 2*p^2 + 2*p + 1

    The `q`-analogue of `n!` is only defined for `n` a nonnegative
    integer (:issue:`11411`)::

        sage: q_factorial(-2)
        Traceback (most recent call last):
        ...
        ValueError: argument (-2) must be a nonnegative integer

    TESTS::

        sage: q_factorial(0).parent()
        Univariate Polynomial Ring in q over Integer Ring
    """
    if n in ZZ:
        if n == 0:
            return q_int(1, q)
        elif n >= 1:
            return prod(q_int(i, q) for i in range(1, n + 1))
    raise ValueError("argument (%s) must be a nonnegative integer" % n)


def q_binomial(n, k, q=None, algorithm='auto'):
    r"""
    Return the `q`-binomial coefficient.

    This is also known as the Gaussian binomial coefficient, and is defined by

    .. MATH::

        \binom{n}{k}_q = \frac{(1-q^n)(1-q^{n-1}) \cdots (1-q^{n-k+1})}
        {(1-q)(1-q^2)\cdots (1-q^k)}.

    See :wikipedia:`Gaussian_binomial_coefficient`.

    If `q` is unspecified, then the variable is the generator `q` for
    a univariate polynomial ring over the integers.

    INPUT:

    - ``n``, ``k`` -- the values `n` and `k` defined above

    - ``q`` -- (default: ``None``) the variable `q`; if ``None``, then use a
      default variable in `\ZZ[q]`

    - ``algorithm`` -- (default: ``'auto'``) the algorithm to use and can be
      one of the following:

      - ``'auto'`` -- automatically choose the algorithm; see the algorithm
        section below
      - ``'naive'`` -- use the naive algorithm
      - ``'cyclotomic'`` -- use cyclotomic algorithm

    ALGORITHM:

    The naive algorithm uses the product formula. The cyclotomic
    algorithm uses a product of cyclotomic polynomials
    (cf. [CH2006]_).

    When the algorithm is set to ``'auto'``, we choose according to
    the following rules:

    - If ``q`` is a polynomial:

      When ``n`` is small or ``k`` is small with respect to ``n``, one
      uses the naive algorithm. When both ``n`` and ``k`` are big, one
      uses the cyclotomic algorithm.

    - If ``q`` is in the symbolic ring (or a symbolic subring), one uses
      the cyclotomic algorithm.

    - Otherwise one uses the naive algorithm, unless ``q`` is a root of
      unity, then one uses the cyclotomic algorithm.

    EXAMPLES:

    By default, the variable is the generator of `\ZZ[q]`::

        sage: from sage.combinat.q_analogues import q_binomial
        sage: g = q_binomial(5,1) ; g
        q^4 + q^3 + q^2 + q + 1
        sage: g.parent()
        Univariate Polynomial Ring in q over Integer Ring

    For `n \geq 0`, the `q`-binomial coefficient vanishes unless `0 \leq k \leq n`::

        sage: q_binomial(4,5)
        0
        sage: q_binomial(5,-1)
        0

    For `k \geq 0`, the `q`-binomial coefficient is extended as a polynomial in `n`::

        sage: q_binomial(-4,1)
        -q^-4 - q^-3 - q^-2 - q^-1
        sage: q_binomial(-2,3)
        -q^-9 - q^-8 - q^-7 - q^-6

    Other variables can be used, given as third parameter::

        sage: p = ZZ['p'].gen()
        sage: q_binomial(4,2,p)
        p^4 + p^3 + 2*p^2 + p + 1

    The third parameter can also be arbitrary values::

        sage: q_binomial(5,1,2) == g.subs(q=2)
        True
        sage: q_binomial(5,1,1)
        5
        sage: q_binomial(4,2,-1)
        2
        sage: q_binomial(4,2,3.14)
        152.030056160000
        sage: R = GF((5, 2), 't')
        sage: t = R.gen(0)
        sage: q_binomial(6, 3, t)
        2*t + 3

    We can also do this for more complicated objects such as matrices or
    symmetric functions::

        sage: q_binomial(4,2,matrix([[2,1],[-1,3]]))
        [ -6  84]
        [-84  78]
        sage: Sym = SymmetricFunctions(QQ)
        sage: s = Sym.schur()
        sage: q_binomial(4,1, s[2]+s[1])
        s[] + s[1] + s[1, 1] + s[1, 1, 1] + 2*s[2] + 4*s[2, 1] + 3*s[2, 1, 1]
         + 4*s[2, 2] + 3*s[2, 2, 1] + s[2, 2, 2] + 3*s[3] + 7*s[3, 1] + 3*s[3, 1, 1]
         + 6*s[3, 2] + 2*s[3, 2, 1] + s[3, 3] + 4*s[4] + 6*s[4, 1] + s[4, 1, 1]
         + 3*s[4, 2] + 3*s[5] + 2*s[5, 1] + s[6]

    TESTS:

    One checks that the first two arguments are integers::

        sage: q_binomial(1/2,1)
        Traceback (most recent call last):
        ...
        TypeError: no conversion of this rational to integer

    One checks that either `k` or `n` is nonnegative::

        sage: q_binomial(-4,-1)
        Traceback (most recent call last):
        ...
        ValueError: either k or n must be nonnegative

    This also works for variables in the symbolic ring::

        sage: z = var('z')                                                              # needs sage.symbolic
        sage: factor(q_binomial(4, 2, z))                                               # needs sage.symbolic
        (z^2 + z + 1)*(z^2 + 1)

    This also works for complex roots of unity::

        sage: q_binomial(10, 4, QQbar(I))                                               # needs sage.rings.number_field
        2

    Note that the symbolic computation works (see :issue:`14982`)::

        sage: q_binomial(10, 4, I)                                                      # needs sage.rings.number_field
        2

    Check that the algorithm does not matter::

        sage: q_binomial(6, 3, algorithm='naive') == q_binomial(6, 3, algorithm='cyclotomic')
        True

    One more test::

        sage: q_binomial(4, 2, Zmod(6)(2), algorithm='naive')
        5

    Check that it works with Python integers::

        sage: r = q_binomial(3r, 2r, 1r); r
        3
        sage: type(r)
        <class 'int'>

    Check that arbitrary polynomials work::

        sage: R.<x> = ZZ[]
        sage: q_binomial(2, 1, x^2 - 1, algorithm='naive')
        x^2
        sage: q_binomial(2, 1, x^2 - 1, algorithm='cyclotomic')
        x^2

    Check that the parent is always the parent of ``q``::

        sage: R.<q> = CyclotomicField(3)
        sage: for algo in ["naive", "cyclotomic"]:
        ....:     for n in range(4):
        ....:         for k in range(4):
        ....:             a = q_binomial(n, k, q, algorithm=algo)
        ....:             assert a.parent() is R

    ::

        sage: q_binomial(2, 1, x^2 - 1, algorithm='quantum')
        Traceback (most recent call last):
        ...
        ValueError: unknown algorithm 'quantum'

    REFERENCES:

    .. [CH2006] William Y.C. Chen and Qing-Hu Hou, *Factors of the Gaussian
       coefficients*, Discrete Mathematics 306 (2006), 1446-1449.
       :doi:`10.1016/j.disc.2006.03.031`

    AUTHORS:

    - Frédéric Chapoton, David Joyner and William Stein
    """
    # sanity checks
    n = ZZ(n)
    k = ZZ(k)
    if k < 0 and n < 0:
        raise ValueError('either k or n must be nonnegative')

    # polynomiality test
    if q is None:
        is_polynomial = True
        if n >= 0:
            q = polygen(ZZ, 'q')
        else:
            q = LaurentPolynomialRing(ZZ, 'q').gen()
    else:
        from sage.rings.polynomial.polynomial_element import Polynomial
        is_polynomial = isinstance(q, Polynomial)

    if n < 0:
        return (-1)**k * q**(k * n - (k * k - k) // 2) * q_binomial(-n + k - 1, k, q=q)

    k = min(n - k, k)  # Pick the smallest k

    # We support non-Sage Elements too, where parent(q) is really
    # type(q). The calls R(0) and R(1) should work in all cases to
    # generate the correct 0 and 1 elements.
    R = parent(q)
    zero = R(0)
    one = R(1)

    if k <= 0:
        return one if k == 0 else zero

    # heuristic choice of the fastest algorithm
    if algorithm == 'auto':
        if n <= 70 or k <= n // 4:
            algorithm = 'naive'
        elif is_polynomial:
            algorithm = 'cyclotomic'
        else:
            import sage.rings.abc
            if isinstance(R, sage.rings.abc.SymbolicRing):
                algorithm = 'cyclotomic'
            else:
                algorithm = 'naive'

    # the algorithms
    while algorithm == 'naive':
        denom = prod(one - q**i for i in range(1, k + 1))
        if not denom:  # q is a root of unity, use the cyclotomic algorithm
            algorithm = 'cyclotomic'
            break
        else:
            num = prod(one - q**i for i in range(n - k + 1, n + 1))
            try:
                try:
                    return num // denom
                except TypeError:
                    return num / denom
            except (TypeError, ZeroDivisionError):
                # use substitution instead
                return q_binomial(n, k)(q)
    if algorithm == 'cyclotomic':
        from sage.rings.polynomial.cyclotomic import cyclotomic_value
        return prod(cyclotomic_value(d, q)
                    for d in range(2, n + 1)
                    if (n//d) != (k//d) + ((n-k)//d))
    else:
        raise ValueError("unknown algorithm {!r}".format(algorithm))


def gaussian_binomial(n, k, q=None, algorithm='auto'):
    r"""
    This is an alias of :func:`q_binomial`.

    See :func:`q_binomial` for the full documentation.

    EXAMPLES::

        sage: gaussian_binomial(4,2)
        q^4 + q^3 + 2*q^2 + q + 1
    """
    return q_binomial(n, k, q, algorithm)


def q_multinomial(seq, q=None, binomial_algorithm='auto'):
    r"""
    Return the `q`-multinomial coefficient.

    This is also known as the Gaussian multinomial coefficient, and is
    defined by

    .. MATH::

        \binom{n}{k_1, k_2, \ldots, k_m}_q = \frac{[n]_q!}
        {[k_1]_q! [k_2]_q! \cdots [k_m]_q!}

    where `n = k_1 + k_2 + \cdots + k_m`.

    If `q` is unspecified, then the variable is the generator `q` for
    a univariate polynomial ring over the integers.

    INPUT:

    - ``seq`` -- an iterable of the values `k_1` to `k_m` defined above

    - ``q`` -- (default: ``None``) the variable `q`; if ``None``, then use a
      default variable in `\ZZ[q]`

    - ``binomial_algorithm`` -- (default: ``'auto'``) the algorithm to use
      in :meth:`~sage.combinat.q_analogues.q_binomial`; see possible values
      there

    ALGORITHM:

    We use the equivalent formula

    .. MATH::

        \binom{k_1 + \cdots + k_m}{k_1, \ldots, k_m}_q
        = \prod_{i=1}^m \binom{\sum_{j=1}^i k_j}{k_i}_q.

    EXAMPLES::

        sage: from sage.combinat.q_analogues import q_multinomial
        sage: q_multinomial([1,2,1])
        q^5 + 2*q^4 + 3*q^3 + 3*q^2 + 2*q + 1
        sage: q_multinomial([1,2,1], q=1) == multinomial([1,2,1])
        True
        sage: q_multinomial((3,2)) == q_binomial(5,3)
        True
        sage: q_multinomial([])
        1
    """
    binomials = []
    partial_sum = 0
    for elem in seq:
        partial_sum += elem
        binomials.append(q_binomial(partial_sum, elem, q=q, algorithm=binomial_algorithm))
    return prod(binomials)


gaussian_multinomial = q_multinomial


def q_catalan_number(n, q=None, m=1):
    """
    Return the `q`-Catalan number of index `n`.

    INPUT:

    - ``q`` -- (optional) variable

    - ``m`` -- (optional) integer; to get instead the ``m``-Fuss-Catalan numbers

    If `q` is unspecified, then it defaults to using the generator `q` for
    a univariate polynomial ring over the integers.

    There are several `q`-Catalan numbers. This procedure
    returns the one which can be written using the `q`-binomial coefficients.

    EXAMPLES::

        sage: from sage.combinat.q_analogues import q_catalan_number
        sage: q_catalan_number(4)
        q^12 + q^10 + q^9 + 2*q^8 + q^7 + 2*q^6 + q^5 + 2*q^4 + q^3 + q^2 + 1

        sage: p = ZZ['p'].0
        sage: q_catalan_number(4, p)
        p^12 + p^10 + p^9 + 2*p^8 + p^7 + 2*p^6 + p^5 + 2*p^4 + p^3 + p^2 + 1

        sage: q_catalan_number(3, m=2)
        q^12 + q^10 + q^9 + q^8 + q^7 + 2*q^6 + q^5 + q^4 + q^3 + q^2 + 1

    TESTS:

    The `q`-Catalan number of index `n` is only defined for `n` a
    nonnegative integer (:issue:`11411`)::

        sage: q_catalan_number(-2)
        Traceback (most recent call last):
        ...
        ValueError: argument (-2) must be a nonnegative integer

        sage: q_catalan_number(3).parent()
        Univariate Polynomial Ring in q over Integer Ring
        sage: q_catalan_number(0).parent()
        Univariate Polynomial Ring in q over Integer Ring
    """
    if n in ZZ:
        if n in {0, 1}:
            return q_int(1, q)
        if n >= 2:
            return (prod(q_int(j, q)
                         for j in range(m * n + 2, (m + 1) * n + 1)) //
                    prod(q_int(j, q)
                         for j in range(2, n + 1)))
    raise ValueError(f"argument ({n}) must be a nonnegative integer")


def qt_catalan_number(n):
    """
    Return the `q,t`-Catalan number of index `n`.

    EXAMPLES::

        sage: from sage.combinat.q_analogues import qt_catalan_number
        sage: qt_catalan_number(1)
        1
        sage: qt_catalan_number(2)
        q + t
        sage: qt_catalan_number(3)
        q^3 + q^2*t + q*t^2 + t^3 + q*t
        sage: qt_catalan_number(4)
        q^6 + q^5*t + q^4*t^2 + q^3*t^3 + q^2*t^4 + q*t^5 + t^6 + q^4*t + q^3*t^2 + q^2*t^3 + q*t^4 + q^3*t + q^2*t^2 + q*t^3

    The `q,t`-Catalan number of index `n` is only defined for `n` a
    nonnegative integer (:issue:`11411`)::

        sage: qt_catalan_number(-2)
        Traceback (most recent call last):
        ...
        ValueError: argument (-2) must be a nonnegative integer
    """
    if n in ZZ and n >= 0:
        ZZqt = ZZ['q', 't']
        d = {}
        for dw in DyckWords(n):
            tup = (dw.area(), dw.bounce())
            d[tup] = d.get(tup, 0) + 1
        return ZZqt(d)
    else:
        raise ValueError("argument (%s) must be a nonnegative integer" % n)


def q_pochhammer(n, a, q=None):
    r"""
    Return the `q`-Pochhammer `(a; q)_n`.

    The `q`-Pochhammer symbol is defined by

    .. MATH::

        (a; q)_n = \prod_{k=0}^{n-1} (1 - aq^k)

    with `(a; q)_0 = 1` for all `a, q` and `n \in \NN`.
    By using the identity

    .. MATH::

        (a; q)_n = \frac{(a; q)_{\infty}}{(aq^n; q)_{\infty}},

    we can extend the definition to `n < 0` by

    .. MATH::

        (a; q)_n = \frac{1}{(aq^n; q)_{-n}}
        = \prod_{k=1}^{-n} \frac{1}{1 - a/q^k}.

    EXAMPLES::

        sage: from sage.combinat.q_analogues import q_pochhammer
        sage: q_pochhammer(3, 1/7)
        6/343*q^3 - 6/49*q^2 - 6/49*q + 6/7
        sage: q_pochhammer(3, 3)
        -18*q^3 + 6*q^2 + 6*q - 2
        sage: q_pochhammer(3, 1)
        0

        sage: R.<q> = ZZ[]
        sage: q_pochhammer(4, q)
        q^10 - q^9 - q^8 + 2*q^5 - q^2 - q + 1
        sage: q_pochhammer(4, q^2)
        q^14 - q^12 - q^11 - q^10 + q^8 + 2*q^7 + q^6 - q^4 - q^3 - q^2 + 1
        sage: q_pochhammer(-3, q)
        1/(-q^9 + q^7 + q^6 + q^5 - q^4 - q^3 - q^2 + 1)

    TESTS::

        sage: q_pochhammer(0, 2)
        1
        sage: q_pochhammer(0, 1)
        1
        sage: q_pochhammer(0, var('a'))                                                 # needs sage.symbolic
        1

    We check that :issue:`25715` is fixed::

        sage: q_pochhammer(0, 3r)
        1

    REFERENCES:

    - :wikipedia:`Q-Pochhammer_symbol`
    """
    if q is None:
        q = polygen(ZZ, 'q')
    if n not in ZZ:
        raise ValueError("{} must be an integer".format(n))
    R = parent(q)
    one = R(1)
    if n < 0:
        return R.prod(one / (one - a/q**-k) for k in range(1, -n+1))
    return R.prod((one - a*q**k) for k in range(n))


@cached_function(key=lambda t, q: (_Partitions(t), q))
def q_jordan(t, q=None):
    r"""
    Return the `q`-Jordan number of `t`.

    If `q` is the power of a prime number, the output is the number of
    complete flags in `\GF{q}^N` (where `N` is the size of `t`) stable
    under a linear nilpotent endomorphism `f_t` whose Jordan type is
    given by `t`, i.e. such that for all `i`:

    .. MATH::

        \dim (\ker f_t^i) = t[0] + \cdots + t[i-1]

    If `q` is unspecified, then it defaults to using the generator `q` for
    a univariate polynomial ring over the integers.

    The result is cached.

    INPUT:

    - ``t`` -- integer partition, or an argument accepted by :class:`Partition`

    - ``q`` -- (default: ``None``) the variable `q`; if ``None``, then use a
      default variable in `\ZZ[q]`

    EXAMPLES::

        sage: from sage.combinat.q_analogues import q_jordan
        sage: [q_jordan(mu, 2) for mu in Partitions(5)]
        [9765, 1029, 213, 93, 29, 9, 1]
        sage: [q_jordan(mu, 2) for mu in Partitions(6)]
        [615195, 40635, 5643, 2331, 1491, 515, 147, 87, 47, 11, 1]
        sage: q_jordan([3,2,1])
        16*q^4 + 24*q^3 + 14*q^2 + 5*q + 1
        sage: q_jordan([2,1], x)                                                        # needs sage.symbolic
        2*x + 1

    If the partition is trivial (i.e. has only one part), we get
    the `q`-factorial (in this case, the nilpotent endomorphism is
    necessarily `0`)::

        sage: from sage.combinat.q_analogues import q_factorial
        sage: q_jordan([5]) == q_factorial(5)
        True
        sage: q_jordan([11], 5) == q_factorial(11, 5)
        True

    TESTS::

        sage: all(multinomial(mu.conjugate()) == q_jordan(mu, 1) for mu in Partitions(6))
        True

    AUTHOR:

    - Xavier Caruso (2012-06-29)
    """
    if q is None:
        q = polygen(ZZ, 'q')

    if all(part == 0 for part in t):
        return parent(q)(1)
    tj = 0
    res = parent(q)(0)
    for i in range(len(t) - 1, -1, -1):
        ti = t[i]
        if ti > tj:
            tp = list(t)
            tp[i] -= 1
            res += q_jordan(tp, q) * q**tj * q_int(ti - tj, q)
            tj = ti
    return res


def q_subgroups_of_abelian_group(la, mu, q=None, algorithm='birkhoff'):
    r"""
    Return the `q`-number of subgroups of type ``mu`` in a finite abelian
    group of type ``la``.

    INPUT:

    - ``la`` -- type of the ambient group as a :class:`Partition`
    - ``mu`` -- type of the subgroup as a :class:`Partition`
    - ``q`` -- (default: ``None``) an indeterminate or a prime number; if
      ``None``, this defaults to `q \in \ZZ[q]`
    - ``algorithm`` -- (default: ``'birkhoff'``) the algorithm to use can be
      one of the following:

      - ``'birkhoff`` -- use the Birkhoff formula from [Bu87]_
      - ``'delsarte'`` -- use the formula from [Delsarte48]_

    OUTPUT:

    The number of subgroups of type ``mu`` in a group of type ``la`` as a
    polynomial in ``q``.

    ALGORITHM:

    Let `q` be a prime number and `\lambda = (\lambda_1, \ldots, \lambda_l)`
    be a partition. A finite abelian `q`-group is of type `\lambda` if it
    is isomorphic to

    .. MATH::

        \ZZ / q^{\lambda_1} \ZZ \times \cdots \times \ZZ / q^{\lambda_l} \ZZ.

    The formula from [Bu87]_ works as follows:
    Let `\lambda` and `\mu` be partitions. Let `\lambda^{\prime}` and
    `\mu^{\prime}` denote the conjugate partitions to `\lambda` and `\mu`,
    respectively. The number of subgroups of type `\mu` in a group of type
    `\lambda` is given by

    .. MATH::

        \prod_{i=1}^{\mu_1} q^{\mu^{\prime}_{i+1}
        (\lambda^{\prime}_i - \mu^{\prime}_i)}
        \binom{\lambda^{\prime}_i - \mu^{\prime}_{i+1}}
        {\mu^{\prime}_i - \mu^{\prime}_{i+1}}_q

    The formula from [Delsarte48]_ works as follows:
    Let `\lambda` and `\mu` be partitions. Let `(s_1, s_2, \ldots, s_l)`
    and `(r_1, r_2, \ldots, r_k)` denote the parts of the partitions
    conjugate to `\lambda` and `\mu` respectively. Let


    .. MATH::

        \mathfrak{F}(\xi_1, \ldots, \xi_k) = \xi_1^{r_2} \xi_2^{r_3} \cdots
        \xi_{k-1}^{r_k} \prod_{i_1=r_2}^{r_1-1} (\xi_1-q^{i_1})
        \prod_{i_2=r_3}^{r_2-1} (\xi_2-q^{i_2}) \cdots
        \prod_{i_k=0}^{r_k-1} (\xi_k-q^{-i_k}).

    Then the number of subgroups of type `\mu` in a group of type `\lambda`
    is given by

    .. MATH::

        \frac{\mathfrak{F}(q^{s_1}, q^{s_2}, \ldots, q^{s_k})}{\mathfrak{F}
        (q^{r_1}, q^{r_2}, \ldots, q^{r_k})}.

    EXAMPLES::

        sage: from sage.combinat.q_analogues import q_subgroups_of_abelian_group
        sage: q_subgroups_of_abelian_group([1,1], [1])
        q + 1
        sage: q_subgroups_of_abelian_group([3,3,2,1], [2,1])
        q^6 + 2*q^5 + 3*q^4 + 2*q^3 + q^2
        sage: R.<t> = QQ[]
        sage: q_subgroups_of_abelian_group([5,3,1], [3,1], t)
        t^4 + 2*t^3 + t^2
        sage: q_subgroups_of_abelian_group([5,3,1], [3,1], 3)
        144
        sage: q_subgroups_of_abelian_group([1,1,1], [1]) == q_subgroups_of_abelian_group([1,1,1], [1,1])
        True
        sage: q_subgroups_of_abelian_group([5], [3])
        1
        sage: q_subgroups_of_abelian_group([1], [2])
        0
        sage: q_subgroups_of_abelian_group([2], [1,1])
        0

    TESTS:

    Check the same examples with ``algorithm='delsarte'``::

        sage: q_subgroups_of_abelian_group([1,1], [1], algorithm='delsarte')
        q + 1
        sage: q_subgroups_of_abelian_group([3,3,2,1], [2,1], algorithm='delsarte')
        q^6 + 2*q^5 + 3*q^4 + 2*q^3 + q^2
        sage: q_subgroups_of_abelian_group([5,3,1], [3,1], t, algorithm='delsarte')
        t^4 + 2*t^3 + t^2
        sage: q_subgroups_of_abelian_group([5,3,1], [3,1], 3, algorithm='delsarte')
        144
        sage: q_subgroups_of_abelian_group([1,1,1], [1], algorithm='delsarte') == q_subgroups_of_abelian_group([1,1,1], [1,1])
        True
        sage: q_subgroups_of_abelian_group([5], [3], algorithm='delsarte')
        1
        sage: q_subgroups_of_abelian_group([1], [2], algorithm='delsarte')
        0
        sage: q_subgroups_of_abelian_group([2], [1,1], algorithm='delsarte')
        0

    Check that :issue:`25715` is fixed::

        sage: parent(q_subgroups_of_abelian_group([2], [1], algorithm='delsarte'))
        Univariate Polynomial Ring in q over Integer Ring
        sage: q_subgroups_of_abelian_group([7,7,1], [])
        1
        sage: q_subgroups_of_abelian_group([7,7,1], [0,0])
        1

    REFERENCES:

    .. [Bu87] Butler, Lynne M. *A unimodality result in the enumeration
       of subgroups of a finite abelian group.* Proceedings of the American
       Mathematical Society 101, no. 4 (1987): 771-775.
       :doi:`10.1090/S0002-9939-1987-0911049-8`

    .. [Delsarte48] \S. Delsarte, *Fonctions de Möbius Sur Les Groupes Abéliens
       Finis*, Annals of Mathematics, second series, Vol. 45, No. 3, (Jul 1948),
       pp. 600-609. http://www.jstor.org/stable/1969047

    AUTHORS:

    - Amritanshu Prasad (2013-06-07): Implemented the Delsarte algorithm
    - Tomer Bauer (2013, 2018): Implemented the Birkhoff algorithm and refactoring
    """
    if q is None:
        q = polygen(ZZ, 'q')
    la_c = _Partitions(la).conjugate()
    mu_c = _Partitions(mu).conjugate()
    k = mu_c.length()
    if not mu_c:
        # There is only one trivial subgroup
        return parent(q)(1)
    if not la_c.contains(mu_c):
        return parent(q)(0)

    if algorithm == 'delsarte':
        def F(args):
            prd = lambda j: prod(args[j]-q**i for i in range(mu_c[j+1], mu_c[j]))
            F1 = prod(args[i]**mu_c[i+1] * prd(i) for i in range(k-1))
            return F1 * prod(args[k-1]-q**i for i in range(mu_c[k-1]))

        return F([q**ss for ss in la_c[:k]])//F([q**rr for rr in mu_c])

    if algorithm == 'birkhoff':
        fac1 = q**(sum(mu_c[i+1] * (la_c[i]-mu_c[i]) for i in range(k-1)))
        fac2 = prod(q_binomial(la_c[i]-mu_c[i+1], mu_c[i]-mu_c[i+1], q=q) for i in range(k-1))
        fac3 = q_binomial(la_c[k-1], mu_c[k-1], q=q)

        return prod([fac1, fac2, fac3])

    raise ValueError("invalid algorithm choice")


@cached_function
def q_stirling_number1(n, k, q=None):
    r"""
    Return the (unsigned) `q`-Stirling number of the first kind.

    This is a `q`-analogue of :func:`sage.combinat.combinat.stirling_number1` .

    INPUT:

    - ``n``, ``k`` -- integers with `1 \leq k \leq n`

    - ``q`` -- variable (default: `q`)

    OUTPUT: a polynomial in the variable `q`

    These polynomials satisfy the recurrence

    .. MATH::

         s_{n,k} = s_{n-1,k-1} + [n-1]_q s_{n-1, k}.

    EXAMPLES::

        sage: from sage.combinat.q_analogues import q_stirling_number1
        sage: q_stirling_number1(4,2)
        q^3 + 3*q^2 + 4*q + 3

        sage: all(stirling_number1(6,k) == q_stirling_number1(6,k)(1)                   # needs sage.libs.gap
        ....:     for k in range(1,6))
        True

        sage: x = polygen(QQ['q'],'x')
        sage: S = sum(q_stirling_number1(5,k)*x**k for k in range(1, 6))
        sage: factor(S)                                                                 # needs sage.libs.singular
        x * (x + 1) * (x + q + 1) * (x + q^2 + q + 1) * (x + q^3 + q^2 + q + 1)

    TESTS::

        sage: q_stirling_number1(-1,2)
        Traceback (most recent call last):
        ...
        ValueError: q-Stirling numbers are not defined for n < 0

    We check that :issue:`25715` is fixed::

        sage: q_stirling_number1(2,1,1r)
        1

    REFERENCES:

    - [Ca1948]_

    - [Ca1954]_
    """
    if q is None:
        q = polygen(ZZ, 'q')
    if n < 0:
        raise ValueError('q-Stirling numbers are not defined for n < 0')
    if n == 0 == k:
        return parent(q)(1)
    if k > n or k < 1:
        return parent(q)(0)
    return (q_stirling_number1(n - 1, k - 1, q=q) +
            q_int(n - 1, q=q) * q_stirling_number1(n - 1, k, q=q))


@cached_function
def q_stirling_number2(n, k, q=None):
    r"""
    Return the (unsigned) `q`-Stirling number of the second kind.

    This is a `q`-analogue of :func:`sage.combinat.combinat.stirling_number2`.

    INPUT:

    - ``n``, ``k`` -- integers with `1 \leq k \leq n`

    - ``q`` -- variable (default: `q`)

    OUTPUT: a polynomial in the variable `q`

    These polynomials satisfy the recurrence

    .. MATH::

         S_{n,k} = q^{k-1} S_{n-1,k-1} + [k]_q s_{n-1, k}.

    EXAMPLES::

        sage: from sage.combinat.q_analogues import q_stirling_number2
        sage: q_stirling_number2(4,2)
        q^3 + 3*q^2 + 3*q

        sage: all(stirling_number2(6,k) == q_stirling_number2(6,k)(1)
        ....:     for k in range(7))
        True


    TESTS::

        sage: q_stirling_number2(-1,2)
        Traceback (most recent call last):
        ...
        ValueError: q-Stirling numbers are not defined for n < 0

    We check that :issue:`25715` is fixed::

        sage: q_stirling_number2(1,0).parent()
        Univariate Polynomial Ring in q over Integer Ring
        sage: q_stirling_number2(2,1,3r)
        1

    REFERENCES:

    - [Mil1978]_
    """
    if q is None:
        q = polygen(ZZ, 'q')
    if n < 0:
        raise ValueError('q-Stirling numbers are not defined for n < 0')
    if n == 0 == k:
        return parent(q)(1)
    if k > n or k <= 0:
        return parent(q)(0)
    return (q**(k-1)*q_stirling_number2(n - 1, k - 1, q=q) +
            q_int(k, q=q) * q_stirling_number2(n - 1, k, q=q))


def number_of_irreducible_polynomials(n, q=None, m=1):
    r"""
    Return the number of monic irreducible polynomials of degree ``n``
    in ``m`` variables over the finite field with ``q`` elements.

    If ``q`` is not given, the result is returned as an integer-valued
    polynomial in `\QQ[q]`.

    INPUT:

    - ``n`` -- positive integer
    - ``q`` -- ``None`` (default) or a prime power
    - ``m`` -- positive integer (default: `1`)

    OUTPUT: integer or integer-valued polynomial over `\QQ`

    EXAMPLES::

        sage: number_of_irreducible_polynomials(8, q=2)
        30
        sage: number_of_irreducible_polynomials(9, q=9)
        43046640
        sage: number_of_irreducible_polynomials(5, q=11, m=3)
        2079650567184059145647246367401741345157369643207055703168

    ::

        sage: poly = number_of_irreducible_polynomials(12); poly
        1/12*q^12 - 1/12*q^6 - 1/12*q^4 + 1/12*q^2
        sage: poly(5) == number_of_irreducible_polynomials(12, q=5)
        True
        sage: poly = number_of_irreducible_polynomials(5, m=3); poly
        q^55 + q^54 + q^53 + q^52 + q^51 + q^50 + ... + 1/5*q^5 - 1/5*q^3 - 1/5*q^2 - 1/5*q
        sage: poly(11) == number_of_irreducible_polynomials(5, q=11, m=3)
        True

    This function is *much* faster than enumerating the polynomials::

        sage: num = number_of_irreducible_polynomials(99, q=101)
        sage: num.bit_length()
        653

    ALGORITHM:

    In the univariate case, classical formula
    `\frac1n \sum_{d\mid n} \mu(n/d) q^d`
    using the Möbius function `\mu`;
    see :func:`moebius`.

    In the multivariate case, formula from [Bodin2007]_,
    independently [Alekseyev2006]_.
    """
    n = ZZ(n)
    if n <= 0:
        raise ValueError('n must be positive')
    if m <= 0:
        raise ValueError('m must be positive')

    if q is None:
        from sage.rings.rational_field import QQ
        q = QQ['q'].gen()  # we produce an integer-valued polynomial in q, but it does not necessarily have integer coefficients

    if m == 1:
        from sage.arith.misc import moebius
        r = sum((moebius(n//d) * q**d for d in n.divisors()), parent(q).zero())
        return r // n

    from sage.functions.other import binomial
    from sage.combinat.partition import Partitions

    def monic_reducible(irreducible, d):
        """
        Compute the number of monic reducible polynomials of degree `d`
        given the numbers of irreducible polynomials up to degree `d-1`.
        """
        res = 0
        for p in Partitions(d+1, max_part=d):
            res += prod(binomial(r+t-1, t) for r, t in zip(irreducible, p.to_exp(d)))
        return res

    r = []
    for d in range(n):
        monic = (q**binomial(d + m, m - 1) - 1) * q**binomial(d + m, m) // (q - 1)
        reducible = monic_reducible(r, d)
        r.append(monic - reducible)

    return r[-1]
