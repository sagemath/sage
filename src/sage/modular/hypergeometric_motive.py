# -*- coding: utf-8 -*-
"""
Hypergeometric motives

This is largely a port of the corresponding package in Magma. One
important conventional difference: the motivic parameter `t` has been replaced
with `1/t` to match the classical literature on hypergeometric series.
(E.g., see [BeukersHeckman]_)

The computation of Euler factors is currently only supported for primes `p`
of good reduction. That is, it is required that `v_p(t) = v_p(t-1) = 0`.

AUTHORS:

- Frédéric Chapoton
- Kiran S. Kedlaya

EXAMPLES::

    sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
    sage: H = Hyp(cyclotomic=([30], [1,2,3,5]))
    sage: H.alpha_beta()
    ([1/30, 7/30, 11/30, 13/30, 17/30, 19/30, 23/30, 29/30],
    [0, 1/5, 1/3, 2/5, 1/2, 3/5, 2/3, 4/5])
    sage: H.M_value() == 30**30 / (15**15 * 10**10 * 6**6)
    True
    sage: H.euler_factor(2, 7)
    T^8 + T^5 + T^3 + 1

REFERENCES:

- [BeukersHeckman]_

- [Benasque2009]_

- [Kat1991]_

- [MagmaHGM]_

- [Fedorov2015]_

- [Roberts2017]_

- [Roberts2015]_

- [BeCoMe]_

- [Watkins]_

"""
# ****************************************************************************
#       Copyright (C) 2017     Frédéric Chapoton
#                              Kiran S. Kedlaya <kskedl@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from collections import defaultdict
from itertools import combinations
from sage.arith.misc import divisors, gcd, euler_phi, moebius, is_prime
from sage.arith.misc import gauss_sum, kronecker_symbol
from sage.combinat.integer_vector_weighted import WeightedIntegerVectors
from sage.functions.generalized import sgn
from sage.functions.log import log
from sage.functions.other import floor, ceil
from sage.misc.cachefunc import cached_method
from sage.misc.functional import cyclotomic_polynomial
from sage.misc.misc_c import prod
from sage.modular.hypergeometric_misc import hgm_coeffs
from sage.rings.fraction_field import FractionField
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.rings.integer_ring import ZZ
from sage.rings.padics.padic_generic_element import gauss_table
from sage.rings.polynomial.polynomial_ring import polygen
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.rational_field import QQ
from sage.schemes.generic.spec import Spec
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.universal_cyclotomic_field import UniversalCyclotomicField

def characteristic_polynomial_from_traces(traces, d, q, i, sign):
    r"""
    Given a sequence of traces `t_1, \dots, t_k`, return the
    corresponding characteristic polynomial with Weil numbers as roots.

    The characteristic polynomial is defined by the generating series

    .. MATH::

        P(T) = \exp\left(- \sum_{k\geq 1} t_k \frac{T^k}{k}\right)

    and should have the property that reciprocals of all roots have
    absolute value `q^{i/2}`.

    INPUT:

    - ``traces`` -- a list of integers `t_1, \dots, t_k`

    - ``d`` -- the degree of the characteristic polynomial

    - ``q`` -- power of a prime number

    - ``i`` -- integer, the weight in the motivic sense

    - ``sign`` -- integer, the sign

    OUTPUT:

    a polynomial

    EXAMPLES::

        sage: from sage.modular.hypergeometric_motive import characteristic_polynomial_from_traces
        sage: characteristic_polynomial_from_traces([1, 1], 1, 3, 0, -1)
        -T + 1
        sage: characteristic_polynomial_from_traces([25], 1, 5, 4, -1)
        -25*T + 1

        sage: characteristic_polynomial_from_traces([3], 2, 5, 1, 1)
        5*T^2 - 3*T + 1
        sage: characteristic_polynomial_from_traces([1], 2, 7, 1, 1)
        7*T^2 - T + 1

        sage: characteristic_polynomial_from_traces([20], 3, 29, 2, 1)
        24389*T^3 - 580*T^2 - 20*T + 1
        sage: characteristic_polynomial_from_traces([12], 3, 13, 2, -1)
        -2197*T^3 + 156*T^2 - 12*T + 1

        sage: characteristic_polynomial_from_traces([36,7620], 4, 17, 3, 1)
        24137569*T^4 - 176868*T^3 - 3162*T^2 - 36*T + 1
        sage: characteristic_polynomial_from_traces([-4,276], 4, 5, 3, 1)
        15625*T^4 + 500*T^3 - 130*T^2 + 4*T + 1
        sage: characteristic_polynomial_from_traces([4,-276], 4, 5, 3, 1)
        15625*T^4 - 500*T^3 + 146*T^2 - 4*T + 1
        sage: characteristic_polynomial_from_traces([22, 484], 4, 31, 2, -1)
        -923521*T^4 + 21142*T^3 - 22*T + 1

    TESTS::

        sage: characteristic_polynomial_from_traces([-36], 4, 17, 3, 1)
        Traceback (most recent call last):
        ...
        ValueError: not enough traces were given
    """
    if len(traces) < d // 2:
        raise ValueError('not enough traces were given')
    if i % 2 and d % 2:
        raise ValueError('i and d may not both be odd')
    t = PowerSeriesRing(QQ, 't').gen()
    ring = PolynomialRing(ZZ, 'T')

    series = sum(- api * t**(i + 1) / (i + 1) for i, api in enumerate(traces))
    series = series.O(d // 2 + 1).exp()
    coeffs = list(series)
    coeffs += [0] * max(0, d // 2 + 1 - len(coeffs))

    data = [0 for _ in range(d + 1)]
    for k in range(d // 2 + 1):
        data[k] = coeffs[k]
    for k in range(d // 2 + 1, d + 1):
        data[k] = sign * coeffs[d - k] * q**(i * (k - d / 2))
    return ring(data)


def enumerate_hypergeometric_data(d, weight=None):
    r"""
    Return an iterator over parameters of hypergeometric motives (up to swapping).

    INPUT:

    - ``d`` -- the degree

    - ``weight`` -- optional integer, to specify the motivic weight

    EXAMPLES::

        sage: from sage.modular.hypergeometric_motive import enumerate_hypergeometric_data as enum
        sage: l = [H for H in enum(6, weight=2) if H.hodge_numbers()[0] == 1]
        sage: len(l)
        112
    """
    bound = 2 * d * d  # to make sure that phi(n) <= d
    possible = [(i, euler_phi(i)) for i in range(1, bound + 1)
                if euler_phi(i) <= d]
    poids = [z[1] for z in possible]
    N = len(poids)
    vectors = WeightedIntegerVectors(d, poids)

    def formule(u):
        return [possible[j][0] for j in range(N) for _ in range(u[j])]

    for a, b in combinations(vectors, 2):
        if not any(a[j] and b[j] for j in range(N)):
            H = HypergeometricData(cyclotomic=(formule(a), formule(b)))
            if weight is None or H.weight() == weight:
                yield H


def possible_hypergeometric_data(d, weight=None):
    """
    Return the list of possible parameters of hypergeometric motives (up to swapping).

    INPUT:

    - ``d`` -- the degree

    - ``weight`` -- optional integer, to specify the motivic weight

    EXAMPLES::

        sage: from sage.modular.hypergeometric_motive import possible_hypergeometric_data as P
        sage: [len(P(i,weight=2)) for i in range(1, 7)]
        [0, 0, 10, 30, 93, 234]
    """
    return list(enumerate_hypergeometric_data(d, weight))


def cyclotomic_to_alpha(cyclo):
    """
    Convert a list of indices of cyclotomic polynomials
    to a list of rational numbers.

    The input represents a product of cyclotomic polynomials.

    The output is the list of arguments of the roots of the
    given product of cyclotomic polynomials.

    This is the inverse of :func:`alpha_to_cyclotomic`.

    EXAMPLES::

        sage: from sage.modular.hypergeometric_motive import cyclotomic_to_alpha
        sage: cyclotomic_to_alpha([1])
        [0]
        sage: cyclotomic_to_alpha([2])
        [1/2]
        sage: cyclotomic_to_alpha([5])
        [1/5, 2/5, 3/5, 4/5]
        sage: cyclotomic_to_alpha([1,2,3,6])
        [0, 1/6, 1/3, 1/2, 2/3, 5/6]
        sage: cyclotomic_to_alpha([2,3])
        [1/3, 1/2, 2/3]
    """
    alpha = []
    for d in cyclo:
        for k in ZZ(d).coprime_integers(d):
            alpha.append(QQ((k, d)))
    return sorted(alpha)


def alpha_to_cyclotomic(alpha):
    """
    Convert from a list of rationals arguments to a list of integers.

    The input represents arguments of some roots of unity.

    The output represent a product of cyclotomic polynomials with exactly
    the given roots. Note that the multiplicity of `r/s` in the list
    must be independent of `r`; otherwise, a ``ValueError`` will be raised.

    This is the inverse of :func:`cyclotomic_to_alpha`.

    EXAMPLES::

        sage: from sage.modular.hypergeometric_motive import alpha_to_cyclotomic
        sage: alpha_to_cyclotomic([0])
        [1]
        sage: alpha_to_cyclotomic([1/2])
        [2]
        sage: alpha_to_cyclotomic([1/5,2/5,3/5,4/5])
        [5]
        sage: alpha_to_cyclotomic([0, 1/6, 1/3, 1/2, 2/3, 5/6])
        [1, 2, 3, 6]
        sage: alpha_to_cyclotomic([1/3,2/3,1/2])
        [2, 3]
    """
    cyclo = []
    Alpha = list(alpha)
    while Alpha:
        q = QQ(Alpha.pop())
        n = q.numerator()
        d = q.denominator()
        for k in d.coprime_integers(d):
            if k != n:
                try:
                    Alpha.remove(QQ((k, d)))
                except ValueError:
                    raise ValueError("multiplicities not balanced")
        cyclo.append(d)
    return sorted(cyclo)


def capital_M(n):
    """
    Auxiliary function, used to describe the canonical scheme.

    INPUT:

    - ``n`` -- an integer

    OUTPUT:

    a rational

    EXAMPLES::

        sage: from sage.modular.hypergeometric_motive import capital_M
        sage: [capital_M(i) for i in range(1,8)]
        [1, 4, 27, 64, 3125, 432, 823543]
    """
    n = ZZ(n)
    return QQ.prod(d ** (d * moebius(n / d)) for d in divisors(n))


def cyclotomic_to_gamma(cyclo_up, cyclo_down):
    """
    Convert a quotient of products of cyclotomic polynomials
    to a quotient of products of polynomials `x^n - 1`.

    INPUT:

    - ``cyclo_up`` -- list of indices of cyclotomic polynomials in the numerator
    - ``cyclo_down`` -- list of indices of cyclotomic polynomials in the denominator

    OUTPUT:

    a dictionary mapping an integer `n` to the power of `x^n - 1` that
    appears in the given product

    EXAMPLES::

        sage: from sage.modular.hypergeometric_motive import cyclotomic_to_gamma
        sage: cyclotomic_to_gamma([6], [1])
        {2: -1, 3: -1, 6: 1}
    """
    dico = defaultdict(int)
    for d in cyclo_up:
        dico[d] += 1
    for d in cyclo_down:
        dico[d] -= 1

    resu = defaultdict(int)
    for n in dico:
        for d in divisors(n):
            resu[d] += moebius(n / d) * dico[n]

    return {d: resu[d] for d in resu if resu[d]}


def gamma_list_to_cyclotomic(galist):
    r"""
    Convert a quotient of products of polynomials `x^n - 1`
    to a quotient of products of cyclotomic polynomials.

    INPUT:

    - ``galist`` -- a list of integers, where an integer `n` represents
      the power `(x^{|n|} - 1)^{\operatorname{sgn}(n)}`

    OUTPUT:

    a pair of list of integers, where `k` represents the cyclotomic
    polynomial `\Phi_k`

    EXAMPLES::

        sage: from sage.modular.hypergeometric_motive import gamma_list_to_cyclotomic
        sage: gamma_list_to_cyclotomic([-1, -1, 2])
        ([2], [1])

        sage: gamma_list_to_cyclotomic([-1, -1, -1, -3, 6])
        ([2, 6], [1, 1, 1])

        sage: gamma_list_to_cyclotomic([-1, 2, 3, -4])
        ([3], [4])

        sage: gamma_list_to_cyclotomic([8,2,2,2,-6,-4,-3,-1])
        ([2, 2, 8], [3, 3, 6])
    """
    resu = defaultdict(int)
    for n in galist:
        eps = sgn(n)
        for d in divisors(abs(n)):
            resu[d] += eps

    return (sorted(d for d in resu for k in range(resu[d])),
            sorted(d for d in resu for k in range(-resu[d])))


class HypergeometricData(object):
    _gauss_table = {}

    def __init__(self, cyclotomic=None, alpha_beta=None, gamma_list=None):
        r"""
        Creation of hypergeometric motives.

        INPUT:

        three possibilities are offered, each describing a quotient
        of products of cyclotomic polynomials.

        - ``cyclotomic`` -- a pair of lists of nonnegative integers,
          each integer `k` represents a cyclotomic polynomial `\Phi_k`

        - ``alpha_beta`` -- a pair of lists of rationals,
          each rational represents a root of unity

        - ``gamma_list`` -- a pair of lists of nonnegative integers,
          each integer `n` represents a polynomial `x^n - 1`

        In the last case, it is also allowed to send just one list of signed
        integers where signs indicate to which part the integer belongs to.

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
            sage: Hyp(cyclotomic=([2],[1]))
            Hypergeometric data for [1/2] and [0]

            sage: Hyp(alpha_beta=([1/2],[0]))
            Hypergeometric data for [1/2] and [0]
            sage: Hyp(alpha_beta=([1/5,2/5,3/5,4/5],[0,0,0,0]))
            Hypergeometric data for [1/5, 2/5, 3/5, 4/5] and [0, 0, 0, 0]

            sage: Hyp(gamma_list=([5],[1,1,1,1,1]))
            Hypergeometric data for [1/5, 2/5, 3/5, 4/5] and [0, 0, 0, 0]
            sage: Hyp(gamma_list=([5,-1,-1,-1,-1,-1]))
            Hypergeometric data for [1/5, 2/5, 3/5, 4/5] and [0, 0, 0, 0]
        """
        if gamma_list is not None:
            if isinstance(gamma_list[0], (list, tuple)):
                pos, neg = gamma_list
                gamma_list = pos + [-u for u in neg]
            cyclotomic = gamma_list_to_cyclotomic(gamma_list)
        if cyclotomic is not None:
            cyclo_up, cyclo_down = cyclotomic
            if any(x in cyclo_up for x in cyclo_down):
                raise ValueError('overlapping parameters not allowed')
            deg = sum(euler_phi(x) for x in cyclo_down)
            up_deg = sum(euler_phi(x) for x in cyclo_up)
            if up_deg != deg:
                msg = 'not the same degree: {} != {}'.format(up_deg, deg)
                raise ValueError(msg)
            cyclo_up.sort()
            cyclo_down.sort()
            alpha = cyclotomic_to_alpha(cyclo_up)
            beta = cyclotomic_to_alpha(cyclo_down)
        elif alpha_beta is not None:
            alpha, beta = alpha_beta
            if len(alpha) != len(beta):
                raise ValueError('alpha and beta not of the same length')
            alpha = sorted(u - floor(u) for u in alpha)
            beta = sorted(u - floor(u) for u in beta)
            cyclo_up = alpha_to_cyclotomic(alpha)
            cyclo_down = alpha_to_cyclotomic(beta)
            deg = sum(euler_phi(x) for x in cyclo_down)

        self._cyclo_up = tuple(cyclo_up)
        self._cyclo_down = tuple(cyclo_down)
        self._alpha = tuple(alpha)
        self._beta = tuple(beta)
        self._deg = deg
        self._gamma_array = cyclotomic_to_gamma(cyclo_up, cyclo_down)
        self._trace_coeffs = {}
        up = QQ.prod(capital_M(d) for d in cyclo_up)
        down = QQ.prod(capital_M(d) for d in cyclo_down)
        self._M_value = up / down
        if 0 in alpha:
            self._swap = HypergeometricData(alpha_beta=(beta, alpha))
        if self.weight() % 2:
            self._sign_param = 1
        else:
            if (deg % 2) != (0 in alpha):
                self._sign_param = prod(cyclotomic_polynomial(v).disc()
                                        for v in cyclo_down)
            else:
                self._sign_param = prod(cyclotomic_polynomial(v).disc()
                                        for v in cyclo_up)

    # --- Internals ---
    def __repr__(self):
        """
        Return the string representation.

        This displays the rational arguments of the roots of unity.

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
            sage: Hyp(alpha_beta=([1/2],[0]))
            Hypergeometric data for [1/2] and [0]
        """
        txt = "Hypergeometric data for {} and {}"
        return txt.format(list(self._alpha), list(self._beta))

    def __eq__(self, other):
        """
        Return whether two data are equal.

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
            sage: H1 = Hyp(alpha_beta=([1/2],[0]))
            sage: H2 = Hyp(cyclotomic=([6,2],[1,1,1]))
            sage: H1 == H1
            True
            sage: H1 == H2
            False
        """
        return (self._alpha == other._alpha and
                self._beta == other._beta)

    def __ne__(self, other):
        """
        Return whether two data are unequal.

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
            sage: H1 = Hyp(alpha_beta=([1/2],[0]))
            sage: H2 = Hyp(cyclotomic=([6,2],[1,1,1]))
            sage: H1 != H1
            False
            sage: H1 != H2
            True
        """
        return not (self == other)

    def __hash__(self):
        """
        Return a hash for ``self``.

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
            sage: H1 = Hyp(alpha_beta=([1/2],[0]))
            sage: h = hash(H1)
        """
        return hash((self._alpha, self._beta))

    # --- Parameters and invariants ---
    def cyclotomic_data(self):
        """
        Return the pair of tuples of indices of cyclotomic polynomials.

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
            sage: Hyp(alpha_beta=([1/2],[0])).cyclotomic_data()
            ([2], [1])
        """
        return (list(self._cyclo_up), list(self._cyclo_down))

    def alpha_beta(self):
        """
        Return the pair of lists of rational arguments.

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
            sage: Hyp(alpha_beta=([1/2],[0])).alpha_beta()
            ([1/2], [0])
        """
        return (list(self._alpha), list(self._beta))

    def alpha(self):
        """
        Return the first tuple of rational arguments.

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
            sage: Hyp(alpha_beta=([1/2],[0])).alpha()
            [1/2]
        """
        return list(self._alpha)

    def beta(self):
        """
        Return the second tuple of rational arguments.

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
            sage: Hyp(alpha_beta=([1/2],[0])).beta()
            [0]
        """
        return list(self._beta)

    def defining_polynomials(self):
        """
        Return the pair of products of cyclotomic polynomials.

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
            sage: Hyp(alpha_beta=([1/4,3/4],[0,0])).defining_polynomials()
            (x^2 + 1, x^2 - 2*x + 1)
        """
        up = prod(cyclotomic_polynomial(d) for d in self._cyclo_up)
        down = prod(cyclotomic_polynomial(d) for d in self._cyclo_down)
        return (up, down)

    def gamma_array(self):
        r"""
        Return the dictionary `\{v: \gamma_v\}` for the expression

        .. MATH::

            \prod_v (T^v - 1)^{\gamma_v}

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
            sage: Hyp(alpha_beta=([1/2],[0])).gamma_array()
            {1: -2, 2: 1}
            sage: Hyp(cyclotomic=([6,2],[1,1,1])).gamma_array()
            {1: -3, 3: -1, 6: 1}
        """
        return dict(self._gamma_array)

    def gamma_list(self):
        r"""
        Return a list of integers describing the `x^n - 1` factors.

        Each integer `n` stands for `(x^{|n|} - 1)^{\operatorname{sgn}(n)}`.

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
            sage: Hyp(alpha_beta=([1/2],[0])).gamma_list()
            [-1, -1, 2]

            sage: Hyp(cyclotomic=([6,2],[1,1,1])).gamma_list()
            [-1, -1, -1, -3, 6]

            sage: Hyp(cyclotomic=([3],[4])).gamma_list()
            [-1, 2, 3, -4]
        """
        gamma = self.gamma_array()
        resu = []
        for v, n in sorted(gamma.items()):
            resu += [sgn(n) * v] * abs(n)
        return resu

    def wild_primes(self):
        r"""
        Return the wild primes.

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
            sage: Hyp(cyclotomic=([3],[4])).wild_primes()
            [2, 3]
            sage: Hyp(cyclotomic=([2,2,2,2,3,3,3,6,6],[1,1,4,5,9])).wild_primes()
            [2, 3, 5]
        """
        gamma = self.gamma_array()
        return sorted(set([p for n in gamma.keys() for (p, _) in n.factor()]))

    def zigzag(self, x, flip_beta=False):
        r"""
        Count ``alpha``'s at most ``x`` minus ``beta``'s at most ``x``.

        This function is used to compute the weight and the Hodge numbers.
        With `flip_beta` set to True, replace each `b` in `\beta` with `1-b`.

        .. SEEALSO::

            :meth:`weight`, :meth:`hodge_numbers`

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
            sage: H = Hyp(alpha_beta=([1/6,1/3,2/3,5/6],[1/8,3/8,5/8,7/8]))
            sage: [H.zigzag(x) for x in [0, 1/3, 1/2]]
            [0, 1, 0]
            sage: H = Hyp(cyclotomic=([5],[1,1,1,1]))
            sage: [H.zigzag(x) for x in [0,1/6,1/4,1/2,3/4,5/6]]
            [-4, -4, -3, -2, -1, 0]

        """
        alpha = self._alpha
        beta = self._beta
        if flip_beta:
            return(sum(1 for a in alpha if a <= x) -
                   sum(1 for b in beta if 1 - b <= x))
        else:
            return(sum(1 for a in alpha if a <= x) -
                   sum(1 for b in beta if b <= x))

    def weight(self):
        """
        Return the motivic weight of this motivic data.

        EXAMPLES:

        With rational inputs::

            sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
            sage: Hyp(alpha_beta=([1/2],[0])).weight()
            0
            sage: Hyp(alpha_beta=([1/4,3/4],[0,0])).weight()
            1
            sage: Hyp(alpha_beta=([1/6,1/3,2/3,5/6],[0,0,1/4,3/4])).weight()
            1
            sage: H = Hyp(alpha_beta=([1/6,1/3,2/3,5/6],[1/8,3/8,5/8,7/8]))
            sage: H.weight()
            1

        With cyclotomic inputs::

            sage: Hyp(cyclotomic=([6,2],[1,1,1])).weight()
            2
            sage: Hyp(cyclotomic=([6],[1,2])).weight()
            0
            sage: Hyp(cyclotomic=([8],[1,2,3])).weight()
            0
            sage: Hyp(cyclotomic=([5],[1,1,1,1])).weight()
            3
            sage: Hyp(cyclotomic=([5,6],[1,1,2,2,3])).weight()
            1
            sage: Hyp(cyclotomic=([3,8],[1,1,1,2,6])).weight()
            2
            sage: Hyp(cyclotomic=([3,3],[2,2,4])).weight()
            1

        With gamma list input::

            sage: Hyp(gamma_list=([8,2,2,2],[6,4,3,1])).weight()
            3
        """
        alpha = self._alpha
        beta = self._beta
        D = [self.zigzag(x) for x in alpha + beta]
        return ZZ(max(D) - min(D) - 1)

    def degree(self):
        """
        Return the degree.

        This is the sum of the Hodge numbers.

        .. SEEALSO::

            :meth:`hodge_numbers`

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
            sage: Hyp(alpha_beta=([1/2],[0])).degree()
            1
            sage: Hyp(gamma_list=([2,2,4],[8])).degree()
            4
            sage: Hyp(cyclotomic=([5,6],[1,1,2,2,3])).degree()
            6
            sage: Hyp(cyclotomic=([3,8],[1,1,1,2,6])).degree()
            6
            sage: Hyp(cyclotomic=([3,3],[2,2,4])).degree()
            4
        """
        return self._deg

    def hodge_numbers(self):
        """
        Return the Hodge numbers.

        .. SEEALSO::

            :meth:`degree`, :meth:`hodge_polynomial`, :meth:`hodge_polygon`

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
            sage: H = Hyp(cyclotomic=([3],[6]))
            sage: H.hodge_numbers()
            [1, 1]

            sage: H = Hyp(cyclotomic=([4],[1,2]))
            sage: H.hodge_numbers()
            [2]

            sage: H = Hyp(gamma_list=([8,2,2,2],[6,4,3,1]))
            sage: H.hodge_numbers()
            [1, 2, 2, 1]

            sage: H = Hyp(gamma_list=([5],[1,1,1,1,1]))
            sage: H.hodge_numbers()
            [1, 1, 1, 1]

            sage: H = Hyp(gamma_list=[6,1,-4,-3])
            sage: H.hodge_numbers()
            [1, 1]

            sage: H = Hyp(gamma_list=[-3]*4 + [1]*12)
            sage: H.hodge_numbers()
            [1, 1, 1, 1, 1, 1, 1, 1]

        REFERENCES:

        - [Fedorov2015]_
        """
        alpha = [(x, 'a') for x in self._alpha]
        beta = [(x, 'b') for x in self._beta]
        height = 0
        hodge = defaultdict(int)
        for x, letter in sorted(alpha + beta):
            if letter == 'a':
                hodge[height] += 1
                height += 1
            else:
                height -= 1
        return [hodge[i] for i in sorted(hodge)]

    def hodge_polynomial(self):
        """
        Return the Hodge polynomial.

        .. SEEALSO::

            :meth:`hodge_numbers`, :meth:`hodge_polygon_vertices`, :meth:`hodge_function`

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
            sage: H = Hyp(cyclotomic=([6,10],[3,12]))
            sage: H.hodge_polynomial()
            (T^3 + 2*T^2 + 2*T + 1)/T^2
            sage: H = Hyp(cyclotomic=([2,2,2,2,3,3,3,6,6],[1,1,4,5,9]))
            sage: H.hodge_polynomial()
            (T^5 + 3*T^4 + 3*T^3 + 3*T^2 + 3*T + 1)/T^2
        """
        alpha = self._alpha

        def z(x):
            return alpha.count(x)

        T = polygen(ZZ, 'T')
        return sum(T ** (self.zigzag(a, flip_beta=True) - z(a)) *
                   (T**z(a) - 1) // (T - 1)
                   for a in set(alpha))

    def hodge_function(self, x):
        """
        Evaluate the Hodge polygon as a function.

        .. SEEALSO::

            :meth:`hodge_numbers`, :meth:`hodge_polynomial`, :meth:`hodge_polygon_vertices`

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
            sage: H = Hyp(cyclotomic=([6,10],[3,12]))
            sage: H.hodge_function(3)
            2
            sage: H.hodge_function(4)
            4
        """
        d = self.degree()
        hn = self.hodge_numbers()
        if x < 0:
            return 0
        i = 0
        j = 0
        k = 0
        while (i < d and i < x):
            i += hn[k]
            j += k * hn[k]
            k += 1
        if i < x:
            return j
        return j - (i - x) * (k - 1)

    def hodge_polygon_vertices(self):
        """
        Return the vertices of the Hodge polygon.

        .. SEEALSO::

            :meth:`hodge_numbers`, :meth:`hodge_polynomial`, :meth:`hodge_function`

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
            sage: H = Hyp(cyclotomic=([6,10],[3,12]))
            sage: H.hodge_polygon_vertices()
            [(0, 0), (1, 0), (3, 2), (5, 6), (6, 9)]
            sage: H = Hyp(cyclotomic=([2,2,2,2,3,3,3,6,6],[1,1,4,5,9]))
            sage: H.hodge_polygon_vertices()
            [(0, 0), (1, 0), (4, 3), (7, 9), (10, 18), (13, 30), (14, 35)]
        """
        lst = [(0, 0)]
        hn = self.hodge_numbers()
        for i in range(len(hn)):
            lst.append((lst[-1][0] + hn[i], lst[-1][1] + i * hn[i]))
        return lst

    def M_value(self):
        """
        Return the `M` coefficient that appears in the trace formula.

        OUTPUT:

        a rational

        .. SEEALSO:: :meth:`canonical_scheme`

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
            sage: H = Hyp(alpha_beta=([1/6,1/3,2/3,5/6],[1/8,3/8,5/8,7/8]))
            sage: H.M_value()
            729/4096
            sage: Hyp(alpha_beta=(([1/2,1/2,1/2,1/2],[0,0,0,0]))).M_value()
            256
            sage: Hyp(cyclotomic=([5],[1,1,1,1])).M_value()
            3125
        """
        return self._M_value

    def is_primitive(self):
        """
        Return whether this data is primitive.

        .. SEEALSO::

            :meth:`primitive_index`, :meth:`primitive_data`

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
            sage: Hyp(cyclotomic=([3],[4])).is_primitive()
            True
            sage: Hyp(gamma_list=[-2, 4, 6, -8]).is_primitive()
            False
            sage: Hyp(gamma_list=[-3, 6, 9, -12]).is_primitive()
            False
        """
        return self.primitive_index() == 1

    def primitive_index(self):
        """
        Return the primitive index.

        .. SEEALSO::

            :meth:`is_primitive`, :meth:`primitive_data`

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
            sage: Hyp(cyclotomic=([3],[4])).primitive_index()
            1
            sage: Hyp(gamma_list=[-2, 4, 6, -8]).primitive_index()
            2
            sage: Hyp(gamma_list=[-3, 6, 9, -12]).primitive_index()
            3
        """
        return gcd(self.gamma_list())

    def has_symmetry_at_one(self):
        """
        If ``True``, the motive H(t=1) is a direct sum of two motives.

        Note that simultaneous exchange of (t,1/t) and (alpha,beta)
        always gives the same motive.

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
            sage: Hyp(alpha_beta=[[1/2]*16,[0]*16]).has_symmetry_at_one()
            True

        REFERENCES:

        - [Roberts2017]_
        """
        beta_twist = self.twist()._beta
        return self.degree() % 2 == 0 and self._alpha == beta_twist

    def canonical_scheme(self, t=None):
        """
        Return the canonical scheme.

        This is a scheme that contains this hypergeometric motive in its cohomology.

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
            sage: H = Hyp(cyclotomic=([3],[4]))
            sage: H.gamma_list()
            [-1, 2, 3, -4]
            sage: H.canonical_scheme()
            Spectrum of Quotient of Multivariate Polynomial Ring
            in X0, X1, Y0, Y1 over Fraction Field of Univariate Polynomial Ring
            in t over Rational Field by the ideal
            (X0 + X1 - 1, Y0 + Y1 - 1, (-t)*X0^2*X1^3 + 27/64*Y0*Y1^4)

            sage: H = Hyp(gamma_list=[-2, 3, 4, -5])
            sage: H.canonical_scheme()
            Spectrum of Quotient of Multivariate Polynomial Ring
            in X0, X1, Y0, Y1 over Fraction Field of Univariate Polynomial Ring
            in t over Rational Field by the ideal
            (X0 + X1 - 1, Y0 + Y1 - 1, (-t)*X0^3*X1^4 + 1728/3125*Y0^2*Y1^5)

        REFERENCES:

        [Kat1991]_, section 5.4
        """
        if t is None:
            t = FractionField(QQ['t']).gen()
        basering = t.parent()
        gamma_pos = [u for u in self.gamma_list() if u > 0]
        gamma_neg = [u for u in self.gamma_list() if u < 0]
        N_pos = len(gamma_pos)
        N_neg = len(gamma_neg)
        varX = ['X{}'.format(i) for i in range(N_pos)]
        varY = ['Y{}'.format(i) for i in range(N_neg)]
        ring = PolynomialRing(basering, varX + varY)
        gens = ring.gens()
        X = gens[:N_pos]
        Y = gens[N_pos:]
        eq0 = ring.sum(X) - 1
        eq1 = ring.sum(Y) - 1
        eq2_pos = ring.prod(X[i] ** gamma_pos[i] for i in range(N_pos))
        eq2_neg = ring.prod(Y[j] ** -gamma_neg[j] for j in range(N_neg))

        ideal = ring.ideal([eq0, eq1, self.M_value() * eq2_neg - t * eq2_pos])
        return Spec(ring.quotient(ideal))

    # --- Operations on data ---
    def twist(self):
        r"""
        Return the twist of this data.

        This is defined by adding `1/2` to each rational in `\alpha`
        and `\beta`.

        This is an involution.

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
            sage: H = Hyp(alpha_beta=([1/2],[0]))
            sage: H.twist()
            Hypergeometric data for [0] and [1/2]
            sage: H.twist().twist() == H
            True

            sage: Hyp(cyclotomic=([6],[1,2])).twist().cyclotomic_data()
            ([3], [1, 2])
        """
        alpha = [x + QQ((1, 2)) for x in self._alpha]
        beta = [x + QQ((1, 2)) for x in self._beta]
        return HypergeometricData(alpha_beta=(alpha, beta))

    def swap_alpha_beta(self):
        """
        Return the hypergeometric data with ``alpha`` and ``beta`` exchanged.

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
            sage: H = Hyp(alpha_beta=([1/2],[0]))
            sage: H.swap_alpha_beta()
            Hypergeometric data for [0] and [1/2]
        """
        alpha, beta = self.alpha_beta()
        return HypergeometricData(alpha_beta=(beta, alpha))

    def primitive_data(self):
        """
        Return a primitive version.

        .. SEEALSO::

            :meth:`is_primitive`, :meth:`primitive_index`

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
            sage: H = Hyp(cyclotomic=([3],[4]))
            sage: H2 = Hyp(gamma_list=[-2, 4, 6, -8])
            sage: H2.primitive_data() == H
            True
        """
        g = self.gamma_list()
        d = gcd(g)
        return HypergeometricData(gamma_list=[x / d for x in g])

### L-functions

    def gauss_table(self, p, f, prec):
        """
        Return (and cache) a table of Gauss sums used in the trace formula.

        .. SEEALSO::

            :meth:`gauss_table_full`

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
            sage: H = Hyp(cyclotomic=([3],[4]))
            sage: H.gauss_table(2, 2, 4)
            (4, [1 + 2 + 2^2 + 2^3, 1 + 2 + 2^2 + 2^3, 1 + 2 + 2^2 + 2^3])
        """
        try:
            prec1, gtab = self._gauss_table[p, f]
            if prec1 < prec:
                raise KeyError
        except KeyError:
            use_longs = (p ** prec < 2 ** 31)
            gtab = gauss_table(p, f, prec, use_longs)
            self._gauss_table[p, f] = (prec, gtab)
            prec1 = prec
        return prec1, gtab

    def gauss_table_full(self):
        """
        Return a dict of all stored tables of Gauss sums.

        The result is passed by reference, and is an attribute of the class;
        consequently, modifying the result has global side effects. Use with
        caution.

        .. SEEALSO::

            :meth:`gauss_table`

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
            sage: H = Hyp(cyclotomic=([3],[4]))
            sage: H.euler_factor(2, 7, cache_p=True)
            7*T^2 - 3*T + 1
            sage: H.gauss_table_full()[(7, 1)]
            (2, array('l', [-1, -29, -25, -48, -47, -22]))

        Clearing cached values::

            sage: H = Hyp(cyclotomic=([3],[4]))
            sage: H.euler_factor(2, 7, cache_p=True)
            7*T^2 - 3*T + 1
            sage: d = H.gauss_table_full()
            sage: d.clear() # Delete all entries of this dict
            sage: H1 = Hyp(cyclotomic=([5],[12]))
            sage: d1 = H1.gauss_table_full()
            sage: len(d1.keys()) # No cached values
            0
        """
        return self._gauss_table

    # --- L-functions ---
    @cached_method
    def padic_H_value(self, p, f, t, prec=None, cache_p=False):
        """
        Return the `p`-adic trace of Frobenius, computed using the
        Gross-Koblitz formula.

        If left unspecified, `prec` is set to the minimum `p`-adic precision
        needed to recover the Euler factor.

        If `cache_p` is True, then the function caches an intermediate
        result which depends only on `p` and `f`. This leads to a significant
        speedup when iterating over `t`.

        INPUT:

        - `p` -- a prime number

        - `f` -- an integer such that `q = p^f`

        - `t` -- a rational parameter

        - ``prec`` -- precision (optional)

        - ``cache_p`` - a boolean

        OUTPUT:

        an integer

        EXAMPLES:

        From Benasque report [Benasque2009]_, page 8::

            sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
            sage: H = Hyp(alpha_beta=([1/2]*4,[0]*4))
            sage: [H.padic_H_value(3,i,-1) for i in range(1,3)]
            [0, -12]
            sage: [H.padic_H_value(5,i,-1) for i in range(1,3)]
            [-4, 276]
            sage: [H.padic_H_value(7,i,-1) for i in range(1,3)]
            [0, -476]
            sage: [H.padic_H_value(11,i,-1) for i in range(1,3)]
            [0, -4972]

        From [Roberts2015]_ (but note conventions regarding `t`)::

            sage: H = Hyp(gamma_list=[-6,-1,4,3])
            sage: t = 189/125
            sage: H.padic_H_value(13,1,1/t)
            0

        TESTS:

        Check issue from :trac:`28404`::

            sage: H1 = Hyp(cyclotomic=([1,1,1],[6,2]))
            sage: H2 = Hyp(cyclotomic=([6,2],[1,1,1]))
            sage: [H1.padic_H_value(5,1,i) for i in range(2,5)]
            [1, -4, -4]
            sage: [H2.padic_H_value(5,1,i) for i in range(2,5)]
            [-4, 1, -4]

        Check for potential overflow::

            sage: H = Hyp(cyclotomic=[[10,6],[5,4]])
            sage: H.padic_H_value(101, 2, 2)
            -1560629

        Check issue from :trac:`29778`::

            sage: H = Hyp(alpha_beta=([1/5,2/5,3/5,4/5,1/5,2/5,3/5,4/5], [1/4,3/4,1/7,2/7,3/7,4/7,5/7,6/7]))
            sage: try:
            ....:     print(H.padic_H_value(373, 4, 2))
            ....: except ValueError as s:
            ....:     print(s)
            p^f cannot exceed 2^31

        Check error handling for wild and tame primes::

            sage: H = Hyp(alpha_beta=([1/5,2/5,3/5,4/5,1/5,2/5,3/5,4/5], [1/4,3/4,1/7,2/7,3/7,4/7,5/7,6/7]))
            sage: try:
            ....:     print(H.padic_H_value(5, 1, 2))
            ....: except NotImplementedError as s:
            ....:     print(s)
            p is wild
            sage: try:
            ....:     print(H.padic_H_value(3, 1, 3))
            ....: except NotImplementedError as s:
            ....:     print(s)
            p is tame

        REFERENCES:

        - [MagmaHGM]_
        """
        alpha = self._alpha
        beta = self._beta
        t = QQ(t)
        if not is_prime(p):
            raise ValueError('p not prime')
        if not all(x.denominator() % p for x in self._alpha + self._beta):
            raise NotImplementedError('p is wild')
        if (t.numerator()*t.denominator() % p == 0 or (t-1) % p == 0):
            raise NotImplementedError('p is tame')

        if 0 in alpha:
            return self._swap.padic_H_value(p, f, ~t, prec)
        q = p ** f
        if q > 2 ** 31:
            raise ValueError("p^f cannot exceed 2^31")

        m = defaultdict(int)
        for b in beta:
            u = b * (q - 1)
            if u.is_integer():
                m[u] += 1
        M = self.M_value()
        D = -min(self.zigzag(x, flip_beta=True) for x in alpha + beta)
        # also: D = (self.weight() + 1 - m[0]) // 2

        if prec is None:
            prec = ceil((self.weight() * f) / 2 + log(2*self.degree()+1, p))
        use_longs = (p ** prec < 2 ** 31)

        gamma = self._gamma_array
        if cache_p:
            try:
                trcoeffs = self._trace_coeffs[p, f]
            except KeyError:
                gtab_prec, gtab = self.gauss_table(p, f, prec)
                trcoeffs = hgm_coeffs(p, f, prec, gamma, m, D, gtab, gtab_prec, use_longs)
                self._trace_coeffs[p, f] = trcoeffs
        else:
            gtab = gauss_table(p, f, prec, use_longs)
            trcoeffs = hgm_coeffs(p, f, prec, gamma, m, D, gtab, prec, use_longs)
        sigma = trcoeffs[p-2]
        p_ring = sigma.parent()
        teich = p_ring.teichmuller(M/t)
        for i in range(p-3, -1, -1):
            sigma = sigma * teich + trcoeffs[i]
        resu = ZZ(-1) ** m[0] * sigma / (1 - q)
        return IntegerModRing(p**prec)(resu).lift_centered()

    trace = padic_H_value

    @cached_method
    def H_value(self, p, f, t, ring=None):
        """
        Return the trace of the Frobenius, computed in terms of Gauss sums
        using the hypergeometric trace formula.

        INPUT:

        - `p` -- a prime number

        - `f` -- an integer such that `q = p^f`

        - `t` -- a rational parameter

        - ``ring`` -- optional (default ``UniversalCyclotomicfield``)

        The ring could be also ``ComplexField(n)`` or ``QQbar``.

        OUTPUT:

        an integer

        .. WARNING::

            This is apparently working correctly as can be tested
            using ComplexField(70) as value ring.

            Using instead UniversalCyclotomicfield, this is much
            slower than the `p`-adic version :meth:`padic_H_value`.

        EXAMPLES:

        With values in the UniversalCyclotomicField (slow)::

            sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
            sage: H = Hyp(alpha_beta=([1/2]*4,[0]*4))
            sage: [H.H_value(3,i,-1) for i in range(1,3)]
            [0, -12]
            sage: [H.H_value(5,i,-1) for i in range(1,3)]
            [-4, 276]
            sage: [H.H_value(7,i,-1) for i in range(1,3)]  # not tested
            [0, -476]
            sage: [H.H_value(11,i,-1) for i in range(1,3)]  # not tested
            [0, -4972]
            sage: [H.H_value(13,i,-1) for i in range(1,3)]  # not tested
            [-84, -1420]

        With values in ComplexField::

            sage: [H.H_value(5,i,-1, ComplexField(60)) for i in range(1,3)]
            [-4, 276]

        Check issue from :trac:`28404`::

            sage: H1 = Hyp(cyclotomic=([1,1,1],[6,2]))
            sage: H2 = Hyp(cyclotomic=([6,2],[1,1,1]))
            sage: [H1.H_value(5,1,i) for i in range(2,5)]
            [1, -4, -4]
            sage: [H2.H_value(5,1,QQ(i)) for i in range(2,5)]
            [-4, 1, -4]

        TESTS:

        Check issue from :trac:`29778`::

            sage: H = Hyp(alpha_beta=([1/5,2/5,3/5,4/5,1/5,2/5,3/5,4/5], [1/4,3/4,1/7,2/7,3/7,4/7,5/7,6/7]))
            sage: try:
            ....:     print(H.padic_H_value(373, 4, 2))
            ....: except ValueError as s:
            ....:     print(s)
            p^f cannot exceed 2^31

        Check error handling for wild and tame primes::

            sage: H = Hyp(alpha_beta=([1/5,2/5,3/5,4/5,1/5,2/5,3/5,4/5], [1/4,3/4,1/7,2/7,3/7,4/7,5/7,6/7]))
            sage: try:
            ....:     print(H.padic_H_value(5, 1, 2))
            ....: except NotImplementedError as s:
            ....:     print(s)
            p is wild
            sage: try:
            ....:     print(H.padic_H_value(3, 1, 3))
            ....: except NotImplementedError as s:
            ....:     print(s)
            p is tame

        REFERENCES:

        - [BeCoMe]_ (Theorem 1.3)
        - [Benasque2009]_
        """
        alpha = self._alpha
        beta = self._beta
        t = QQ(t)
        if not is_prime(p):
            raise ValueError('p not prime')
        if not all(x.denominator() % p for x in self._alpha + self._beta):
            raise NotImplementedError('p is wild')
        if (t.numerator()*t.denominator() % p == 0 or (t-1) % p == 0):
            raise NotImplementedError('p is tame')

        if 0 in alpha:
            return self._swap.H_value(p, f, ~t, ring)
        if ring is None:
            ring = UniversalCyclotomicField()
        gamma = self.gamma_array()
        q = p**f

        m = {r: beta.count(QQ((r, q - 1))) for r in range(q - 1)}
        D = -min(self.zigzag(x, flip_beta=True) for x in alpha + beta)
        # also: D = (self.weight() + 1 - m[0]) // 2
        M = self.M_value()

        Fq = GF((p, f))
        gen = Fq.multiplicative_generator()
        zeta_q = ring.zeta(q - 1)

        tM = Fq(M / t)
        for k in range(q - 1):
            if gen ** k == tM:
                teich = zeta_q ** k
                break

        gauss_table = [gauss_sum(zeta_q ** r, Fq) for r in range(q - 1)]

        sigma = sum(q**(D + m[0] - m[r]) *
                    prod(gauss_table[(-v * r) % (q - 1)]**gv
                         for v, gv in gamma.items()) *
                    teich ** r
                    for r in range(q - 1))
        resu = ZZ(-1) ** m[0] / (1 - q) * sigma
        if not ring.is_exact():
            resu = resu.real_part().round()
        return resu

    def sign(self, t, p):
        """
        Return the sign of the functional equation for the Euler factor of the motive `H_t` at the prime `p`.

        For odd weight, the sign of the functional equation is +1. For even
        weight, the sign is computed by a recipe found in 11.1 of [Watkins]_
        (when 0 is not in alpha).

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
            sage: H = Hyp(cyclotomic=([6,2],[1,1,1]))
            sage: H.weight(), H.degree()
            (2, 3)
            sage: [H.sign(1/4,p) for p in [5,7,11,13,17,19]]
            [1, 1, -1, -1, 1, 1]

            sage: H = Hyp(alpha_beta=([1/12,5/12,7/12,11/12],[0,1/2,1/2,1/2]))
            sage: H.weight(), H.degree()
            (2, 4)
            sage: t = -5
            sage: [H.sign(1/t,p) for p in [11,13,17,19,23,29]]
            [-1, -1, -1, 1, 1, 1]

        We check that :trac:`28404` is fixed::

            sage: H = Hyp(cyclotomic=([1,1,1],[6,2]))
            sage: [H.sign(4,p) for p in [5,7,11,13,17,19]]
            [1, 1, -1, -1, 1, 1]

        """
        t = QQ(t)
        if 0 in self._alpha:
            return self._swap.sign(~t, p)
        d = self.degree()
        w = self.weight()

        if w % 2:  # sign is always +1 for odd weight
            sign = 1
        elif d % 2:
            sign = -kronecker_symbol((1 - t) * self._sign_param, p)
        else:
            sign = kronecker_symbol(t * (t - 1) * self._sign_param, p)
        return sign

    @cached_method
    def euler_factor(self, t, p, cache_p=False):
        """
        Return the Euler factor of the motive `H_t` at prime `p`.

        INPUT:

        - `t` -- rational number, not 0 or 1

        - `p` -- prime number of good reduction

        OUTPUT:

        a polynomial

        See [Benasque2009]_ for explicit examples of Euler factors.

        For odd weight, the sign of the functional equation is +1. For even
        weight, the sign is computed by a recipe found in 11.1 of [Watkins]_.

        EXAMPLES::

            sage: from sage.modular.hypergeometric_motive import HypergeometricData as Hyp
            sage: H = Hyp(alpha_beta=([1/2]*4,[0]*4))
            sage: H.euler_factor(-1, 5)
            15625*T^4 + 500*T^3 - 130*T^2 + 4*T + 1

            sage: H = Hyp(gamma_list=[-6,-1,4,3])
            sage: H.weight(), H.degree()
            (1, 2)
            sage: t = 189/125
            sage: [H.euler_factor(1/t,p) for p in [11,13,17,19,23,29]]
            [11*T^2 + 4*T + 1,
            13*T^2 + 1,
            17*T^2 + 1,
            19*T^2 + 1,
            23*T^2 + 8*T + 1,
            29*T^2 + 2*T + 1]

            sage: H = Hyp(cyclotomic=([6,2],[1,1,1]))
            sage: H.weight(), H.degree()
            (2, 3)
            sage: [H.euler_factor(1/4,p) for p in [5,7,11,13,17,19]]
            [125*T^3 + 20*T^2 + 4*T + 1,
             343*T^3 - 42*T^2 - 6*T + 1,
             -1331*T^3 - 22*T^2 + 2*T + 1,
             -2197*T^3 - 156*T^2 + 12*T + 1,
             4913*T^3 + 323*T^2 + 19*T + 1,
             6859*T^3 - 57*T^2 - 3*T + 1]

            sage: H = Hyp(alpha_beta=([1/12,5/12,7/12,11/12],[0,1/2,1/2,1/2]))
            sage: H.weight(), H.degree()
            (2, 4)
            sage: t = -5
            sage: [H.euler_factor(1/t,p) for p in [11,13,17,19,23,29]]
            [-14641*T^4 - 1210*T^3 + 10*T + 1,
             -28561*T^4 - 2704*T^3 + 16*T + 1,
             -83521*T^4 - 4046*T^3 + 14*T + 1,
             130321*T^4 + 14440*T^3 + 969*T^2 + 40*T + 1,
             279841*T^4 - 25392*T^3 + 1242*T^2 - 48*T + 1,
             707281*T^4 - 7569*T^3 + 696*T^2 - 9*T + 1]

        This is an example of higher degree::

            sage: H = Hyp(cyclotomic=([11], [7, 12]))
            sage: H.euler_factor(2, 13)
            371293*T^10 - 85683*T^9 + 26364*T^8 + 1352*T^7 - 65*T^6 + 394*T^5 - 5*T^4 + 8*T^3 + 12*T^2 - 3*T + 1
            sage: H.euler_factor(2, 19) # long time
            2476099*T^10 - 651605*T^9 + 233206*T^8 - 77254*T^7 + 20349*T^6 - 4611*T^5 + 1071*T^4 - 214*T^3 + 34*T^2 - 5*T + 1

        TESTS::

             sage: H1 = Hyp(alpha_beta=([1,1,1],[1/2,1/2,1/2]))
             sage: H2 = H1.swap_alpha_beta()
             sage: H1.euler_factor(-1, 3)
             27*T^3 + 3*T^2 + T + 1
             sage: H2.euler_factor(-1, 3)
             27*T^3 + 3*T^2 + T + 1
             sage: H = Hyp(alpha_beta=([0,0,0,1/3,2/3],[1/2,1/5,2/5,3/5,4/5]))
             sage: H.euler_factor(5,7)
             16807*T^5 - 686*T^4 - 105*T^3 - 15*T^2 - 2*T + 1

        Check for precision downsampling::

            sage: H = Hyp(cyclotomic=[[3],[4]])
            sage: H.euler_factor(2, 11, cache_p=True)
            11*T^2 - 3*T + 1
            sage: H = Hyp(cyclotomic=[[12],[1,2,6]])
            sage: H.euler_factor(2, 11, cache_p=True)
            -T^4 + T^3 - T + 1

        Check issue from :trac:`29778`::

            sage: H = Hyp(alpha_beta=([1/5,2/5,3/5,4/5,1/5,2/5,3/5,4/5], [1/4,3/4,1/7,2/7,3/7,4/7,5/7,6/7]))
            sage: try:
            ....:     print(H.euler_factor(2, 373))
            ....: except ValueError as s:
            ....:     print(s)
            p^f cannot exceed 2^31

        Check error handling for wild and tame primes::

            sage: H = Hyp(alpha_beta=([1/5,2/5,3/5,4/5,1/5,2/5,3/5,4/5], [1/4,3/4,1/7,2/7,3/7,4/7,5/7,6/7]))
            sage: try:
            ....:     print(H.euler_factor(2, 5))
            ....: except NotImplementedError as s:
            ....:     print(s)
            p is wild
            sage: try:
            ....:     print(H.euler_factor(3, 3))
            ....: except NotImplementedError as s:
            ....:     print(s)
            p is tame

        REFERENCES:

        - [Roberts2015]_
        - [Watkins]_
        """
        t = QQ(t)
        if t in [0, 1]:
            raise ValueError('invalid t')
        alpha = self._alpha
        if 0 in alpha:
            return self._swap.euler_factor(~t, p)

        if not is_prime(p):
            raise ValueError('p not prime')
        if not all(x.denominator() % p for x in self._alpha + self._beta):
            raise NotImplementedError('p is wild')
        if (t.numerator()*t.denominator() % p == 0 or (t-1) % p == 0):
            raise NotImplementedError('p is tame')
        # now p is good
        d = self.degree()
        bound = d // 2
        if p ** bound > 2 ** 31:
            raise ValueError("p^f cannot exceed 2^31")

        traces = [self.padic_H_value(p, i + 1, t, cache_p=cache_p)
                  for i in range(bound)]

        w = self.weight()
        sign = self.sign(t, p)
        return characteristic_polynomial_from_traces(traces, d, p, w, sign)
