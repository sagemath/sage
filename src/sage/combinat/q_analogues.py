# -*- coding: utf-8 -*-
r"""
`q`-Analogues
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
# python3
from __future__ import division

from sage.misc.cachefunc import cached_function
from sage.misc.all import prod
from sage.structure.element import parent
from sage.rings.all import ZZ
from sage.combinat.dyck_word import DyckWords
from sage.combinat.partition import Partition


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
        (-q^2 - q - 1)/q^3
        sage: p = ZZ['p'].0
        sage: q_int(3,p)
        p^2 + p + 1
        sage: q_int(3/2)
        Traceback (most recent call last):
        ...
        ValueError: 3/2 must be an integer

    TESTS:

    We check that :trac:`15805` is fixed::

        sage: from sage.combinat.q_analogues import q_int
        sage: q_int(0).parent()
        Univariate Polynomial Ring in q over Integer Ring
    """
    if not n in ZZ:
        raise ValueError('%s must be an integer' % n)

    if q is None:
        q = ZZ['q'].gen()
    if n == 0:  # Special case
        return q.parent().zero()
    if n >= 0:
        return sum(q**i for i in range(n))
    return -q**n*sum(q**i for i in range(-n))

def q_factorial(n, q=None):
    """
    Returns the `q`-analogue of the factorial `n!`.

    If `q` is unspecified, then it defaults to using the generator `q` for
    a univariate polynomial ring over the integers.

    EXAMPLES::

        sage: from sage.combinat.q_analogues import q_factorial
        sage: q_factorial(3)
        q^3 + 2*q^2 + 2*q + 1
        sage: p = ZZ['p'].0
        sage: q_factorial(3, p)
        p^3 + 2*p^2 + 2*p + 1

    The `q`-analogue of `n!` is only defined for `n` a non-negative
    integer (:trac:`11411`)::

        sage: q_factorial(-2)
        Traceback (most recent call last):
        ...
        ValueError: Argument (-2) must be a nonnegative integer.
    """
    if n in ZZ and n >= 0:
        return prod([q_int(i, q) for i in range(1, n+1)])
    else:
        raise ValueError("Argument (%s) must be a nonnegative integer." %n)

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

    - ``n, k`` -- the values `n` and `k` defined above

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

    - If ``q`` is in the symbolic ring, one uses the cyclotomic algorithm.

    - Otherwise one uses the naive algorithm, unless ``q`` is a root of
      unity, then one uses the cyclotomic algorithm.

    EXAMPLES:

    By default, the variable is the generator of `\ZZ[q]`::

        sage: from sage.combinat.q_analogues import q_binomial
        sage: g = q_binomial(5,1) ; g
        q^4 + q^3 + q^2 + q + 1
        sage: g.parent()
        Univariate Polynomial Ring in q over Integer Ring

    The `q`-binomial coefficient vanishes unless `0 \leq k \leq n`::

        sage: q_binomial(4,5)
        0
        sage: q_binomial(5,-1)
        0

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
        sage: R = GF(25, 't')
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

    One checks that `n` is nonnegative::

        sage: q_binomial(-4,1)
        Traceback (most recent call last):
        ...
        ValueError: n must be nonnegative

    This also works for variables in the symbolic ring::

        sage: z = var('z')
        sage: factor(q_binomial(4,2,z))
        (z^2 + z + 1)*(z^2 + 1)

    This also works for complex roots of unity::

        sage: q_binomial(6, 1, QQbar(I))
        I + 1

    Note that the symbolic computation works (see :trac:`14982`)::

        sage: q_binomial(6, 1, I)
        I + 1

    Check that the algorithm does not matter::

        sage: q_binomial(6,3, algorithm='naive') == q_binomial(6,3, algorithm='cyclotomic')
        True

    One more test::

        sage: q_binomial(4, 2, Zmod(6)(2), algorithm='naive')
        5

    REFERENCES:

    .. [CH2006] William Y.C. Chen and Qing-Hu Hou, *Factors of the Gaussian
       coefficients*, Discrete Mathematics 306 (2006), 1446-1449.
       :doi:`10.1016/j.disc.2006.03.031`

    AUTHORS:

    - Frederic Chapoton, David Joyner and William Stein
    """
    # sanity checks
    n = ZZ(n)
    k = ZZ(k)
    if n < 0:
        raise ValueError('n must be nonnegative')

    # polynomiality test
    if q is None:
        from sage.rings.polynomial.polynomial_ring import polygen
        q = polygen(ZZ, name='q')
        is_polynomial = True
    else:
        from sage.rings.polynomial.polynomial_element import Polynomial
        is_polynomial = isinstance(q, Polynomial)

    R = parent(q)
    zero = R.zero()
    one = R.one()

    if not(0 <= k and k <= n):
        return zero

    k = min(n-k,k) # Pick the smallest k

    # heuristic choice of the fastest algorithm
    if algorithm == 'auto':
        from sage.symbolic.ring import SR
        if is_polynomial:
            if n <= 70 or k <= n // 4:
                algorithm = 'naive'
            else:
                algorithm = 'cyclo_polynomial'
        elif q in SR:
            algorithm = 'cyclo_generic'
        else:
            algorithm = 'naive'
    elif algorithm == 'cyclotomic':
        if is_polynomial:
            algorithm = 'cyclo_polynomial'
        else:
            algorithm = 'cyclo_generic'
    elif algorithm != 'naive':
        raise ValueError("invalid algorithm choice")

    # the algorithms
    if algorithm == 'naive':
        denom = prod(one - q**i for i in range(1, k+1))
        if not denom: # q is a root of unity, use the cyclotomic algorithm
            return cyclotomic_value(n, k, q, algorithm='cyclotomic')
        else:
            num = prod(one - q**i for i in range(n-k+1, n+1))
            try:
                try:
                    return num // denom
                except TypeError:
                    return num / denom
            except (TypeError, ZeroDivisionError):
                # use substitution instead
                return q_binomial(n,k)(q)
    elif algorithm == 'cyclo_generic':
        from sage.rings.polynomial.cyclotomic import cyclotomic_value
        return prod(cyclotomic_value(d,q)
                    for d in range(2,n+1)
                    if (n//d) != (k//d) + ((n-k)//d))
    elif algorithm == 'cyclo_polynomial':
        return prod(R.cyclotomic_polynomial(d)
                    for d in range(2,n+1)
                    if (n//d) != (k//d) + ((n-k)//d))


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


def q_catalan_number(n, q=None):
    """
    Return the `q`-Catalan number of index `n`.

    If `q` is unspecified, then it defaults to using the generator `q` for
    a univariate polynomial ring over the integers.

    There are several `q`-Catalan numbers. This procedure
    returns the one which can be written using the `q`-binomial coefficients.

    EXAMPLES::

        sage: from sage.combinat.q_analogues import q_catalan_number
        sage: q_catalan_number(4)
        q^12 + q^10 + q^9 + 2*q^8 + q^7 + 2*q^6 + q^5 + 2*q^4 + q^3 + q^2 + 1
        sage: p = ZZ['p'].0
        sage: q_catalan_number(4,p)
        p^12 + p^10 + p^9 + 2*p^8 + p^7 + 2*p^6 + p^5 + 2*p^4 + p^3 + p^2 + 1

    The `q`-Catalan number of index `n` is only defined for `n` a
    nonnegative integer (:trac:`11411`)::

        sage: q_catalan_number(-2)
        Traceback (most recent call last):
        ...
        ValueError: Argument (-2) must be a nonnegative integer.

    TESTS::

        sage: q_catalan_number(3).parent()
        Univariate Polynomial Ring in q over Integer Ring
    """
    if n in ZZ and n >= 0:
        return (prod(q_int(j, q) for j in range(n + 2, 2 * n + 1)) //
                prod(q_int(j, q) for j in range(2, n + 1)))
    else:
        raise ValueError("Argument (%s) must be a nonnegative integer." % n)


def qt_catalan_number(n):
    """
    Returns the `q,t`-Catalan number of index `n`.

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
    nonnegative integer (:trac:`11411`)::

        sage: qt_catalan_number(-2)
        Traceback (most recent call last):
        ...
        ValueError: Argument (-2) must be a nonnegative integer.
    """
    if n in ZZ and n >= 0:
        ZZqt = ZZ['q','t']
        d = {}
        for dw in DyckWords(n):
            tup = (dw.area(),dw.bounce())
            d[tup] = d.get(tup,0)+1
        return ZZqt(d)
    else:
        raise ValueError("Argument (%s) must be a nonnegative integer." %n)

@cached_function
def q_jordan(t, q):
    r"""
    INPUT:

    -  `t` -- a partition of an integer

    -  `q` -- an integer or an indeterminate

    OUTPUT:

    If `q` is the power of a prime number, the output is the number of
    complete flags in `F_q^N` (where `N` is the size of `t`) stable
    under a linear nilpotent endomorphism `f` whose Jordan type is
    given by `t`, i.e. such that for all `i`:

    .. MATH::

        \dim (\ker f^i) = t[0] + \cdots + t[i-1]

    If `q` is an indeterminate, the output is a polynomial whose
    values at powers of prime numbers are the previous numbers.

    The result is cached.

    EXAMPLES::

        sage: from sage.combinat.q_analogues import q_jordan
        sage: [q_jordan(mu,2) for mu in Partitions(5)]
        [9765, 1029, 213, 93, 29, 9, 1]
        sage: [q_jordan(mu,2) for mu in Partitions(6)]
        [615195, 40635, 5643, 2331, 1491, 515, 147, 87, 47, 11, 1]

        sage: q=PolynomialRing(ZZ,'q').gen()
        sage: q_jordan(Partition([3,2,1]),q)
        16*q^4 + 24*q^3 + 14*q^2 + 5*q + 1

    If the partition is trivial (i.e. has only one part), we get
    the `q`-factorial (in this case, the nilpotent endomorphism is
    necessarily `0`)::

        sage: from sage.combinat.q_analogues import q_factorial
        sage: q_jordan(Partition([5]),3) == q_factorial(5,3)
        True
        sage: q_jordan(Partition([11]),5) == q_factorial(11,5)
        True

    TESTS::

        sage: q_jordan(Partition([4,3,1]),1)
        Traceback (most recent call last):
        ...
        ValueError: q must not be equal to 1

    AUTHOR:

    - Xavier Caruso (2012-06-29)
    """

    if q == 1:
        raise ValueError("q must not be equal to 1")

    if len(t) == 0:
        return 1
    tj = 0
    res = 0
    for i in range(len(t)-1,-1,-1):
        ti = t[i]
        if ti > tj:
            tp = t.to_list()
            tp[i] -= 1
            res += q_jordan(Partition(tp),q) * ((q**ti - q**tj) // (q-1))
            tj = ti
    return res

def q_subgroups_of_abelian_group(la, mu, q=None, algorithm='birkhoff'):
    r"""
    Return the `q`-number of subgroups of type ``mu`` in a finite abelian
    group of type ``la``.

    INPUT:

    - ``la`` -- type of the ambient group as a :class:`Partition`
    - ``mu`` -- type of the subgroup as a :class:`Partition`
    - ``q`` -- (default: ``None``) an indeterminat or a prime number; if
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

    REFERENCES:

    .. [Bu87] Butler, Lynne M. *A unimodality result in the enumeration
       of subgroups of a finite abelian group.* Proceedings of the American
       Mathematical Society 101, no. 4 (1987): 771-775.
       :doi:`10.1090/S0002-9939-1987-0911049-8`

    .. [Delsarte48] \S. Delsarte, *Fonctions de Möbius Sur Les Groupes Abeliens
       Finis*, Annals of Mathematics, second series, Vol. 45, No. 3, (Jul 1948),
       pp. 600-609. http://www.jstor.org/stable/1969047

    AUTHORS:

    - Amritanshu Prasad (2013-06-07): Implemented the Delsarte algorithm
    - Tomer Bauer (2013-09-26): Implemented the Birkhoff algorithm
    """
    if q is None:
        q = ZZ['q'].gens()[0]
    la_c = Partition(la).conjugate()
    mu_c = Partition(mu).conjugate()
    k = mu_c.length()
    if not la_c.contains(mu_c):
        return q.parent().zero()

    if algorithm == 'delsarte':
        def F(args):
            prd = lambda j: prod(args[j]-q**i for i in range(mu_c[j+1],mu_c[j]))
            F1 = prod(args[i]**mu_c[i+1] * prd(i) for i in range(k-1))
            return F1 * prod(args[k-1]-q**i for i in range(mu_c[k-1]))

        return F([q**ss for ss in la_c[:k]])/F([q**rr for rr in mu_c])

    if algorithm == 'birkhoff':
        fac1 = q**(sum(mu_c[i+1] * (la_c[i]-mu_c[i]) for i in range(k-1)))
        fac2 = prod(q_binomial(la_c[i]-mu_c[i+1], mu_c[i]-mu_c[i+1], q=q) for i in range(k-1))
        fac3 = q_binomial(la_c[k-1], mu_c[k-1], q=q)

        return prod([fac1, fac2, fac3])

    raise ValueError("invalid algorithm choice")

