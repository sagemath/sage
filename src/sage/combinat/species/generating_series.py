r"""
Generating Series

This file makes a number of extensions to lazy power series by
endowing them with some semantic content for how they're to be
interpreted.

This code is based on the work of Ralf Hemmecke and Martin Rubey's
Aldor-Combinat, which can be found at
http://www.risc.uni-linz.ac.at/people/hemmecke/aldor/combinat/index.html.
In particular, the relevant section for this file can be found at
http://www.risc.uni-linz.ac.at/people/hemmecke/AldorCombinat/combinatse10.html.
One notable difference is that we use power-sum symmetric functions
as the coefficients of our cycle index series.

TESTS::

    sage: # needs sage.modules
    sage: from sage.combinat.species.generating_series import CycleIndexSeriesRing
    sage: p = SymmetricFunctions(QQ).power()
    sage: CIS = CycleIndexSeriesRing(QQ)
    sage: geo1 = CIS(lambda i: p([1])^i)
    sage: geo2 = CIS(lambda i: p([2])^(i // 2) if is_even(i) else 0)
    sage: s = geo1 * geo2
    sage: s[0]
    p[]
    sage: s[1]
    p[1]
    sage: s[2]
    p[1, 1] + p[2]
    sage: s[3]
    p[1, 1, 1] + p[2, 1]
"""

# ****************************************************************************
#       Copyright (C) 2008 Mike Hansen <mhansen@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.lazy_series import LazyPowerSeries, LazySymmetricFunction
from sage.rings.lazy_series_ring import LazyPowerSeriesRing, LazySymmetricFunctions
from sage.rings.integer import Integer
from sage.rings.rational_field import QQ
from sage.arith.misc import divisors, factorial
from sage.combinat.partition import Partition, Partitions
from sage.misc.cachefunc import cached_function
from sage.misc.lazy_import import lazy_import

lazy_import('sage.combinat.sf.sf', 'SymmetricFunctions')


class OrdinaryGeneratingSeries(LazyPowerSeries):
    r"""
    A class for ordinary generating series.

    Note that it is just a :class:`LazyPowerSeries` whose elements
    have some extra methods.

    EXAMPLES::

        sage: from sage.combinat.species.generating_series import OrdinaryGeneratingSeriesRing
        sage: R = OrdinaryGeneratingSeriesRing(QQ)
        sage: f = R(lambda n: n)
        sage: f
        z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + O(z^7)
    """

    def count(self, n):
        """
        Return the number of structures on a set of size ``n``.

        INPUT:

        - ``n`` -- the size of the set

        EXAMPLES::

            sage: from sage.combinat.species.generating_series import OrdinaryGeneratingSeriesRing
            sage: R = OrdinaryGeneratingSeriesRing(QQ)
            sage: f = R(range(20))
            sage: f.count(10)
            10
        """
        return self.coefficient(n)

    def counts(self, n):
        """
        Return the number of structures on a set for size ``i`` for
        each ``i`` in ``range(n)``.

        EXAMPLES::

            sage: from sage.combinat.species.generating_series import OrdinaryGeneratingSeriesRing
            sage: R = OrdinaryGeneratingSeriesRing(QQ)
            sage: f = R(range(20))
            sage: f.counts(10)
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        """
        return [self.count(i) for i in range(n)]


class OrdinaryGeneratingSeriesRing(LazyPowerSeriesRing):
    r"""
    Return the ring of ordinary generating series over ``R``.

    Note that it is just a
    :class:`LazyPowerSeriesRing` whose elements have
    some extra methods.

    EXAMPLES::

        sage: from sage.combinat.species.generating_series import OrdinaryGeneratingSeriesRing
        sage: R = OrdinaryGeneratingSeriesRing(QQ); R
        Lazy Taylor Series Ring in z over Rational Field
        sage: [R(lambda n: 1).coefficient(i) for i in range(4)]
        [1, 1, 1, 1]
        sage: R(lambda n: 1).counts(4)
        [1, 1, 1, 1]
        sage: R == loads(dumps(R))
        True

    TESTS:

    We test to make sure that caching works::

        sage: R is OrdinaryGeneratingSeriesRing(QQ)
        True
    """

    def __init__(self, base_ring):
        """
        Initialize ``self``.

        TESTS::

            sage: from sage.combinat.species.generating_series import OrdinaryGeneratingSeriesRing
            sage: OrdinaryGeneratingSeriesRing.options.halting_precision(15)
            sage: R = OrdinaryGeneratingSeriesRing(QQ)
            sage: TestSuite(R).run()

            sage: OrdinaryGeneratingSeriesRing.options._reset()  # reset options
        """
        super().__init__(base_ring, names='z')

    Element = OrdinaryGeneratingSeries


class ExponentialGeneratingSeries(LazyPowerSeries):
    r"""
    A class for ordinary generating series.

    Note that it is just a
    :class:`LazyPowerSeries` whose elements have
    some extra methods.

    EXAMPLES::

        sage: from sage.combinat.species.generating_series import OrdinaryGeneratingSeriesRing
        sage: R = OrdinaryGeneratingSeriesRing(QQ)
        sage: f = R(lambda n: n)
        sage: f
        z + 2*z^2 + 3*z^3 + 4*z^4 + 5*z^5 + 6*z^6 + O(z^7)
    """

    def count(self, n):
        """
        Return the number of structures of size ``n``.

        EXAMPLES::

            sage: from sage.combinat.species.generating_series import ExponentialGeneratingSeriesRing
            sage: R = ExponentialGeneratingSeriesRing(QQ)
            sage: f = R(lambda n: 1)
            sage: [f.count(i) for i in range(7)]
            [1, 1, 2, 6, 24, 120, 720]
        """
        return factorial(n) * self.coefficient(n)

    def counts(self, n):
        """
        Return the number of structures on a set for size ``i`` for
        each ``i`` in ``range(n)``.

        EXAMPLES::

            sage: from sage.combinat.species.generating_series import ExponentialGeneratingSeriesRing
            sage: R = ExponentialGeneratingSeriesRing(QQ)
            sage: f = R(range(20))
            sage: f.counts(5)
            [0, 1, 4, 18, 96]
        """
        return [self.count(i) for i in range(n)]

    def functorial_composition(self, y):
        r"""
        Return the exponential generating series which is the functorial
        composition of ``self`` with ``y``.

        If `f = \sum_{n=0}^{\infty} f_n \frac{x^n}{n!}` and
        `g = \sum_{n=0}^{\infty} g_n \frac{x^n}{n!}`, then
        functorial composition `f \Box g` is defined as

        .. MATH::

             f \Box g = \sum_{n=0}^{\infty} f_{g_n} \frac{x^n}{n!}.

        REFERENCES:

        - Section 2.2 of [BLL1998]_.

        EXAMPLES::

            sage: G = species.SimpleGraphSpecies()
            sage: g = G.generating_series()
            sage: [g.coefficient(i) for i in range(10)]
            [1, 1, 1, 4/3, 8/3, 128/15, 2048/45, 131072/315, 2097152/315, 536870912/2835]

            sage: E = species.SetSpecies()
            sage: E2 = E.restricted(min=2, max=3)
            sage: WP = species.SubsetSpecies()
            sage: P2 = E2*E
            sage: g1 = WP.generating_series()
            sage: g2 = P2.generating_series()
            sage: g1.functorial_composition(g2)[:10]
            [1, 1, 1, 4/3, 8/3, 128/15, 2048/45, 131072/315, 2097152/315, 536870912/2835]
        """
        P = self.parent()
        return P(lambda n: self.count(y.count(n)) / factorial(n), 0)


class ExponentialGeneratingSeriesRing(LazyPowerSeriesRing):
    r"""
    Return the ring of exponential generating series over ``R``.

    Note that it is just a
    :class:`LazyPowerSeriesRing` whose elements have
    some extra methods.

    EXAMPLES::

        sage: from sage.combinat.species.generating_series import ExponentialGeneratingSeriesRing
        sage: R = ExponentialGeneratingSeriesRing(QQ); R
        Lazy Taylor Series Ring in z over Rational Field
        sage: [R(lambda n: 1).coefficient(i) for i in range(4)]
        [1, 1, 1, 1]
        sage: R(lambda n: 1).counts(4)
        [1, 1, 2, 6]

    TESTS:

    We test to make sure that caching works::

        sage: R is ExponentialGeneratingSeriesRing(QQ)
        True
    """

    def __init__(self, base_ring):
        """
        Initialize ``self``.

        TESTS::

            sage: from sage.combinat.species.generating_series import ExponentialGeneratingSeriesRing
            sage: ExponentialGeneratingSeriesRing.options.halting_precision(15)
            sage: R = ExponentialGeneratingSeriesRing(QQ)
            sage: TestSuite(R).run()

            sage: ExponentialGeneratingSeriesRing.options._reset()  # reset options
        """
        super().__init__(base_ring, names='z')

    Element = ExponentialGeneratingSeries


class CycleIndexSeries(LazySymmetricFunction):
    def count(self, t):
        """
        Return the number of structures corresponding to a certain cycle
        type ``t``.

        EXAMPLES::

            sage: # needs sage.modules
            sage: from sage.combinat.species.generating_series import CycleIndexSeriesRing
            sage: p = SymmetricFunctions(QQ).power()
            sage: CIS = CycleIndexSeriesRing(QQ)
            sage: f = CIS([0, p([1]), 2*p([1,1]), 3*p([2,1])])
            sage: f.count([1])
            1
            sage: f.count([1,1])
            4
            sage: f.count([2,1])
            6
        """
        t = Partition(t)
        return t.aut() * self.coefficient_cycle_type(t)

    def coefficient_cycle_type(self, t):
        """
        Return the coefficient of a cycle type ``t`` in ``self``.

        EXAMPLES::

            sage: # needs sage.modules
            sage: from sage.combinat.species.generating_series import CycleIndexSeriesRing
            sage: p = SymmetricFunctions(QQ).power()
            sage: CIS = CycleIndexSeriesRing(QQ)
            sage: f = CIS([0, p([1]), 2*p([1,1]),3*p([2,1])])
            sage: f.coefficient_cycle_type([1])
            1
            sage: f.coefficient_cycle_type([1,1])
            2
            sage: f.coefficient_cycle_type([2,1])
            3
        """
        t = Partition(t)
        p = self.coefficient(t.size())
        return p.coefficient(t)

    def isotype_generating_series(self):
        """
        Return the isotype generating series of ``self``.

        EXAMPLES::

            sage: P = species.PermutationSpecies()
            sage: cis = P.cycle_index_series()                                          # needs sage.modules
            sage: f = cis.isotype_generating_series()                                   # needs sage.modules
            sage: f[:10]                                                                # needs sage.modules
            [1, 1, 2, 3, 5, 7, 11, 15, 22, 30]
        """
        R = self.base_ring()
        OGS = OrdinaryGeneratingSeriesRing(R)
        return OGS(lambda n: self._ogs_gen(n, self._coeff_stream._approximate_order),
                   self._coeff_stream._approximate_order)

    def _ogs_gen(self, n, ao):
        """
        Return a generator for the coefficients of the ordinary generating
        series obtained from a cycle index series.

        EXAMPLES::

            sage: P = species.PermutationSpecies()
            sage: cis = P.cycle_index_series()                                          # needs sage.modules
            sage: [cis._ogs_gen(i, 0) for i in range(10)]                               # needs sage.modules
            [1, 1, 2, 3, 5, 7, 11, 15, 22, 30]
        """
        if n < ao:
            return 0
        return sum(self.coefficient(n).coefficients())

    def generating_series(self):
        """
        Return the generating series of ``self``.

        EXAMPLES::

            sage: P = species.PartitionSpecies()
            sage: cis = P.cycle_index_series()                                          # needs sage.modules
            sage: f = cis.generating_series()                                           # needs sage.modules
            sage: f[:5]                                                                 # needs sage.modules
            [1, 1, 1, 5/6, 5/8]
        """
        R = self.base_ring()
        EGS = ExponentialGeneratingSeriesRing(R)
        return EGS(lambda n: self._egs_gen(n, self._coeff_stream._approximate_order),
                   self._coeff_stream._approximate_order)

    def _egs_gen(self, n, ao):
        """
        Return a generator for the coefficients of the exponential
        generating series obtained from a cycle index series.

        EXAMPLES::

            sage: P = species.PermutationSpecies()
            sage: cis = P.cycle_index_series()                                          # needs sage.modules
            sage: [cis._egs_gen(i, 0) for i in range(10)]                               # needs sage.modules
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        """
        if n < ao:
            return 0
        return self.coefficient(n).coefficient([1]*n)

    def derivative(self, n=1):
        r"""
        Return the species-theoretic `n`-th derivative of ``self``.

        For a cycle index series `F (p_{1}, p_{2}, p_{3}, \ldots)`, its
        derivative is the cycle index series `F' = D_{p_{1}} F` (that is,
        the formal derivative of `F` with respect to the variable `p_{1}`).

        If `F` is the cycle index series of a species `S` then `F'` is the
        cycle index series of an associated species `S'` of `S`-structures
        with a "hole".

        EXAMPLES:

        The species `E` of sets satisfies the relationship `E' = E`::

            sage: E = species.SetSpecies().cycle_index_series()                         # needs sage.modules
            sage: E[:8] == E.derivative()[:8]                                           # needs sage.modules
            True

        The species `C` of cyclic orderings and the species `L` of linear
        orderings satisfy the relationship `C' = L`::

            sage: C = species.CycleSpecies().cycle_index_series()                       # needs sage.modules
            sage: L = species.LinearOrderSpecies().cycle_index_series()                 # needs sage.modules
            sage: L[:8] == C.derivative()[:8]                                           # needs sage.modules
            True
        """
        return self.derivative_with_respect_to_p1(n=n)

    def pointing(self):
        r"""
        Return the species-theoretic pointing of ``self``.

        For a cycle index `F`, its pointing is the cycle index series
        `F^{\bullet} = p_{1} \cdot F'`.

        If `F` is the cycle index series of a species `S` then `F^{\bullet}`
        is the cycle index series of an associated species `S^{\bullet}`
        of `S`-structures with a marked "root".

        EXAMPLES:

        The species `E^{\bullet}` of "pointed sets" satisfies
        `E^{\bullet} = X \cdot E`::

            sage: E = species.SetSpecies().cycle_index_series()                         # needs sage.modules
            sage: X = species.SingletonSpecies().cycle_index_series()                   # needs sage.modules
            sage: E.pointing()[:8] == (X*E)[:8]                                         # needs sage.modules
            True
        """
        X = self.parent()([1], valuation=1)
        return X * self.derivative_with_respect_to_p1()

    def exponential(self):
        r"""
        Return the species-theoretic exponential of ``self``.

        For a cycle index `Z_{F}` of a species `F`, its exponential is the
        cycle index series `Z_{E} \circ Z_{F}`, where `Z_{E}` is the
        :meth:`~sage.combinat.species.generating_series.ExponentialCycleIndexSeries`.

        The exponential `Z_{E} \circ Z_{F}` is then the cycle index series
        of the species `E \circ F` of "sets of `F`-structures".

        EXAMPLES:

        Let `BT` be the species of binary trees, `BF` the species of binary
        forests, and `E` the species of sets. Then we have `BF = E \circ BT`::

            sage: BT = species.BinaryTreeSpecies().cycle_index_series()                 # needs sage.modules
            sage: BF = species.BinaryForestSpecies().cycle_index_series()               # needs sage.modules
            sage: BT.exponential().isotype_generating_series()[:8] == BF.isotype_generating_series()[:8]                # needs sage.modules
            True
        """
        base_ring = self.parent().base_ring().base_ring()
        E = ExponentialCycleIndexSeries(base_ring)
        return E(self)

    def logarithm(self):
        r"""
        Return the combinatorial logarithm of ``self``.

        For a cycle index `Z_{F}` of a species `F`, its logarithm is the
        cycle index series `Z_{\Omega} \circ Z_{F}`, where `Z_{\Omega}` is the
        :meth:`~sage.combinat.species.generating_series.LogarithmCycleIndexSeries`.

        The logarithm `Z_{\Omega} \circ Z_{F}` is then the cycle index series
        of the (virtual) species `\Omega \circ F` of "connected `F`-structures".
        In particular, if `F = E^{+} \circ G` for `E^{+}` the species of
        nonempty sets and `G` some other species, then `\Omega \circ F = G`.

        EXAMPLES:

        Let `G` be the species of nonempty graphs and  `CG` be the species of nonempty connected
        graphs. Then `G = E^{+} \circ CG`, so `CG = \Omega \circ G`::

            sage: G = species.SimpleGraphSpecies().cycle_index_series() - 1             # needs sage.modules
            sage: from sage.combinat.species.generating_series import LogarithmCycleIndexSeries
            sage: CG = LogarithmCycleIndexSeries()(G)                                   # needs sage.modules
            sage: CG.isotype_generating_series()[0:8]                                   # needs sage.modules
            [0, 1, 1, 2, 6, 21, 112, 853]
        """
        base_ring = self.parent().base_ring().base_ring()
        Omega = LogarithmCycleIndexSeries(base_ring)
        return Omega(self)


class CycleIndexSeriesRing(LazySymmetricFunctions):
    r"""
    Return the ring of cycle index series over ``R``.

    This is the ring of formal power series `\Lambda[x]`, where
    `\Lambda` is the ring of symmetric functions over ``R`` in the
    `p`-basis. Its purpose is to house the cycle index series of
    species (in a somewhat nonstandard notation tailored to Sage):
    If `F` is a species, then the *cycle index series* of `F` is
    defined to be the formal power series

    .. MATH::

        \sum_{n \geq 0} \frac{1}{n!} (\sum_{\sigma \in S_n}
        \operatorname{fix} F[\sigma]
        \prod_{z \text{ is a cycle of } \sigma}
        p_{\text{length of } z}) x^n
        \in \Lambda_\QQ [x],

    where `\operatorname{fix} F[\sigma]` denotes the number of
    fixed points of the permutation `F[\sigma]` of `F[n]`. We
    notice that this power series is "equigraded" (meaning that
    its `x^n`-coefficient is homogeneous of degree `n`). A more
    standard convention in combinatorics would be to use
    `x_i` instead of `p_i`, and drop the `x` (that is, evaluate
    the above power series at `x = 1`); but this would be more
    difficult to implement in Sage, as it would be an element
    of a power series ring in infinitely many variables.

    Note that it is just a :class:`LazyPowerSeriesRing` (whose base
    ring is `\Lambda`) whose elements have some extra methods.

    EXAMPLES::

        sage: from sage.combinat.species.generating_series import CycleIndexSeriesRing
        sage: R = CycleIndexSeriesRing(QQ); R                                           # needs sage.modules
        Cycle Index Series Ring over Rational Field
        sage: p = SymmetricFunctions(QQ).p()                                            # needs sage.modules
        sage: R(lambda n: p[n])                                                         # needs sage.modules
        p[] + p[1] + p[2] + p[3] + p[4] + p[5] + p[6] + O^7

    TESTS:

    We test to make sure that caching works::

        sage: R is CycleIndexSeriesRing(QQ)                                             # needs sage.modules
        True
    """
    Element = CycleIndexSeries

    def __init__(self, base_ring, sparse=True):
        """
        Initialize ``self``.

        TESTS::

            sage: from sage.combinat.species.generating_series import CycleIndexSeriesRing
            sage: CycleIndexSeriesRing.options.halting_precision(12)
            sage: R = CycleIndexSeriesRing(QQ)                                          # needs sage.modules
            sage: TestSuite(R).run()                                                    # needs sage.modules

            sage: CycleIndexSeriesRing.options._reset()  # reset options
        """
        p = SymmetricFunctions(base_ring).power()
        super().__init__(p)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.species.generating_series import CycleIndexSeriesRing
            sage: CycleIndexSeriesRing(QQ)                                              # needs sage.modules
            Cycle Index Series Ring over Rational Field
        """
        return "Cycle Index Series Ring over %s" % self.base_ring()


@cached_function
def _exp_term(n, R=QQ):
    r"""
    Compute the order-``n`` term of the cycle index series of the species
    `E` of sets.

    EXAMPLES::

        sage: from sage.combinat.species.generating_series import _exp_term
        sage: [_exp_term(i) for i in range(4)]                                          # needs sage.modules
        [p[], p[1], 1/2*p[1, 1] + 1/2*p[2], 1/6*p[1, 1, 1] + 1/2*p[2, 1] + 1/3*p[3]]
    """
    p = SymmetricFunctions(R).power()
    return sum(p(part) / part.aut() for part in Partitions(n))


@cached_function
def ExponentialCycleIndexSeries(R=QQ):
    r"""
    Return the cycle index series of the species `E` of sets.

    This cycle index satisfies

    .. MATH::

        Z_{E} = \sum_{n \geq 0} \sum_{\lambda \vdash n}
        \frac{p_{\lambda}}{z_{\lambda}}.

    EXAMPLES::

        sage: from sage.combinat.species.generating_series import ExponentialCycleIndexSeries
        sage: ExponentialCycleIndexSeries()[:5]                                         # needs sage.modules
        [p[], p[1], 1/2*p[1, 1] + 1/2*p[2], 1/6*p[1, 1, 1] + 1/2*p[2, 1]
         + 1/3*p[3], 1/24*p[1, 1, 1, 1] + 1/4*p[2, 1, 1] + 1/8*p[2, 2]
         + 1/3*p[3, 1] + 1/4*p[4]]
    """
    CIS = CycleIndexSeriesRing(R)
    return CIS(_exp_term)


@cached_function
def _cl_term(n, R=QQ):
    r"""
    Compute the order-``n`` term of the cycle index series of
    the virtual species `\Omega`, the compositional inverse of
    the species `E^{+}` of nonempty sets.

    EXAMPLES::

        sage: from sage.combinat.species.generating_series import _cl_term
        sage: [_cl_term(i) for i in range(4)]                                           # needs sage.modules
        [0, p[1], -1/2*p[1, 1] - 1/2*p[2], 1/3*p[1, 1, 1] - 1/3*p[3]]
    """
    n = Integer(n)  # check that n is an integer

    p = SymmetricFunctions(R).power()

    res = p.zero()
    if n == 1:
        res = p([1])
    elif n > 1:
        res = 1/n * ((-1)**(n-1) * p([1])**n - sum(d * p([n // d]).plethysm(_cl_term(d, R)) for d in divisors(n)[:-1]))

    return res


@cached_function
def LogarithmCycleIndexSeries(R=QQ):
    r"""
    Return the cycle index series of the virtual species `\Omega`, the
    compositional inverse of the species `E^{+}` of nonempty sets.

    The notion of virtual species is treated thoroughly in [BLL1998]_.
    The specific algorithm used here to compute the cycle index of
    `\Omega` is found in [Labelle2008]_.

    EXAMPLES:

    The virtual species `\Omega` is 'properly virtual', in the sense that
    its cycle index has negative coefficients::

        sage: from sage.combinat.species.generating_series import LogarithmCycleIndexSeries
        sage: LogarithmCycleIndexSeries()[0:4]                                          # needs sage.modules
        [0, p[1], -1/2*p[1, 1] - 1/2*p[2], 1/3*p[1, 1, 1] - 1/3*p[3]]

    Its defining property is that `\Omega \circ E^{+} = E^{+} \circ \Omega = X`
    (that is, that composition with `E^{+}` in both directions yields the
    multiplicative identity `X`)::

        sage: Eplus = sage.combinat.species.set_species.SetSpecies(min=1).cycle_index_series()      # needs sage.modules
        sage: LogarithmCycleIndexSeries()(Eplus)[0:4]                                   # needs sage.modules
        [0, p[1], 0, 0]
    """
    CIS = CycleIndexSeriesRing(R)
    return CIS(_cl_term)
