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

    sage: from sage.combinat.species.stream import Stream, _integers_from
    sage: from sage.combinat.species.generating_series import CycleIndexSeriesRing
    sage: p = SymmetricFunctions(QQ).power()
    sage: CIS = CycleIndexSeriesRing(QQ)

::

    sage: geo1 = CIS((p([1])^i  for i in _integers_from(0)))
    sage: geo2 = CIS((p([2])^i  for i in _integers_from(0)))
    sage: s = geo1 * geo2
    sage: s[0]
    p[]
    sage: s[1]
    p[1] + p[2]
    sage: s[2]
    p[1, 1] + p[2, 1] + p[2, 2]
    sage: s[3]
    p[1, 1, 1] + p[2, 1, 1] + p[2, 2, 1] + p[2, 2, 2]

Whereas the coefficients of the above test are homogeneous with
respect to total degree, the following test groups with respect to
weighted degree where each variable x_i has weight i.

::

    sage: def g():
    ....:     for i in _integers_from(0):
    ....:         yield p([2])^i
    ....:         yield p(0)
    sage: geo1 = CIS((p([1])^i  for i in _integers_from(0)))
    sage: geo2 = CIS(g())
    sage: s = geo1 * geo2
    sage: s[0]
    p[]
    sage: s[1]
    p[1]
    sage: s[2]
    p[1, 1] + p[2]
    sage: s[3]
    p[1, 1, 1] + p[2, 1]
    sage: s[4]
    p[1, 1, 1, 1] + p[2, 1, 1] + p[2, 2]

REFERENCES:

.. [BLL] \F. Bergeron, G. Labelle, and P. Leroux.
   "Combinatorial species and tree-like structures".
   Encyclopedia of Mathematics and its Applications, vol. 67, Cambridge Univ. Press. 1998.
.. [BLL-Intro] Francois Bergeron, Gilbert Labelle, and Pierre Leroux.
   "Introduction to the Theory of Species of Structures", March 14, 2008.
"""

#*****************************************************************************
#       Copyright (C) 2008 Mike Hansen <mhansen@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .series import LazyPowerSeriesRing, LazyPowerSeries
from .stream import Stream, _integers_from
from sage.rings.all import Integer, RationalField
from sage.arith.all import moebius, gcd, lcm, divisors
from sage.combinat.partition import Partition, Partitions
from functools import partial
from sage.combinat.sf.sf import SymmetricFunctions
from sage.misc.cachefunc import cached_function
from sage.arith.misc import factorial


@cached_function
def OrdinaryGeneratingSeriesRing(R):
    """
    Return the ring of ordinary generating series over ``R``.

    Note that it is just a
    :class:`LazyPowerSeriesRing` whose elements have
    some extra methods.

    EXAMPLES::

        sage: from sage.combinat.species.generating_series import OrdinaryGeneratingSeriesRing
        sage: R = OrdinaryGeneratingSeriesRing(QQ); R
        Lazy Power Series Ring over Rational Field
        sage: R([1]).coefficients(4)
        [1, 1, 1, 1]
        sage: R([1]).counts(4)
        [1, 1, 1, 1]

    TESTS:

    We test to make sure that caching works.

    ::

        sage: R is OrdinaryGeneratingSeriesRing(QQ)
        True
    """
    return OrdinaryGeneratingSeriesRing_class(R)


class OrdinaryGeneratingSeriesRing_class(LazyPowerSeriesRing):
    def __init__(self, R):
        """
        EXAMPLES::

            sage: from sage.combinat.species.generating_series import OrdinaryGeneratingSeriesRing
            sage: R = OrdinaryGeneratingSeriesRing(QQ)
            sage: R == loads(dumps(R))
            True
        """
        LazyPowerSeriesRing.__init__(self, R, element_class=OrdinaryGeneratingSeries)


class OrdinaryGeneratingSeries(LazyPowerSeries):
    def count(self, n):
        """
        Return the number of structures on a set of size ``n``.

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


@cached_function
def ExponentialGeneratingSeriesRing(R):
    """
    Return the ring of exponential generating series over ``R``.

    Note that it is just a
    :class:`LazyPowerSeriesRing` whose elements have
    some extra methods.

    EXAMPLES::

        sage: from sage.combinat.species.generating_series import ExponentialGeneratingSeriesRing
        sage: R = ExponentialGeneratingSeriesRing(QQ); R
        Lazy Power Series Ring over Rational Field
        sage: R([1]).coefficients(4)
        [1, 1, 1, 1]
        sage: R([1]).counts(4)
        [1, 1, 2, 6]

    TESTS:

    We test to make sure that caching works.

    ::

        sage: R is ExponentialGeneratingSeriesRing(QQ)
        True
    """
    return ExponentialGeneratingSeriesRing_class(R)


class ExponentialGeneratingSeriesRing_class(LazyPowerSeriesRing):
    def __init__(self, R):
        """
        EXAMPLES::

            sage: from sage.combinat.species.generating_series import ExponentialGeneratingSeriesRing
            sage: R = ExponentialGeneratingSeriesRing(QQ)
            sage: R == loads(dumps(R))
            True
        """
        LazyPowerSeriesRing.__init__(self, R, element_class=ExponentialGeneratingSeries)

class ExponentialGeneratingSeries(LazyPowerSeries):
    def count(self, n):
        """
        Return the number of structures of size ``n``.

        EXAMPLES::

            sage: from sage.combinat.species.generating_series import ExponentialGeneratingSeriesRing
            sage: R = ExponentialGeneratingSeriesRing(QQ)
            sage: f = R([1])
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

             f \Box g = \sum_{n=0}^{\infty} f_{g_n} \frac{x^n}{n!}

        REFERENCES:

        - Section 2.2 of [BLL]_.

        EXAMPLES::

            sage: G = species.SimpleGraphSpecies()
            sage: g = G.generating_series()
            sage: g.coefficients(10)
            [1, 1, 1, 4/3, 8/3, 128/15, 2048/45, 131072/315, 2097152/315, 536870912/2835]
        """
        return self._new(partial(self._functorial_compose_gen, y), lambda a,b: 0, self, y)

    def _functorial_compose_gen(self, y, ao):
        """
        Returns a generator for the coefficients of the functorial
        composition of self with y.

        EXAMPLES::

            sage: E = species.SetSpecies()
            sage: E2 = E.restricted(min=2, max=3)
            sage: WP = species.SubsetSpecies()
            sage: P2 = E2*E
            sage: g1 = WP.generating_series()
            sage: g2 = P2.generating_series()
            sage: g = g1._functorial_compose_gen(g2, 0)
            sage: [next(g) for i in range(10)]
            [1, 1, 1, 4/3, 8/3, 128/15, 2048/45, 131072/315, 2097152/315, 536870912/2835]
        """
        n = 0
        while True:
            yield self.count(y.count(n)) / factorial(n)
            n += 1


def factorial_gen():
    """
    A generator for the factorials starting at 0.

    EXAMPLES::

        sage: from sage.combinat.species.generating_series import factorial_gen
        sage: g = factorial_gen()
        sage: [next(g) for i in range(5)]
        [1, 1, 2, 6, 24]
    """
    z = Integer(1)
    yield z
    yield z
    n = Integer(2)
    while True:
        z *= n
        yield z
        n += 1


@cached_function
def CycleIndexSeriesRing(R):
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
        sage: R = CycleIndexSeriesRing(QQ); R
        Cycle Index Series Ring over Symmetric Functions over Rational Field in the powersum basis
        sage: R([1]).coefficients(4) # This is not combinatorially
        ....:                        # meaningful.
        [1, 1, 1, 1]

    TESTS:

    We test to make sure that caching works.

    ::

        sage: R is CycleIndexSeriesRing(QQ)
        True
    """
    return CycleIndexSeriesRing_class(R)


class CycleIndexSeriesRing_class(LazyPowerSeriesRing):
    def __init__(self, R):
        """
        EXAMPLES::

            sage: from sage.combinat.species.generating_series import CycleIndexSeriesRing
            sage: R = CycleIndexSeriesRing(QQ); R
            Cycle Index Series Ring over Symmetric Functions over Rational Field in the powersum basis
            sage: R == loads(dumps(R))
            True
        """
        R = SymmetricFunctions(R).power()
        LazyPowerSeriesRing.__init__(self, R, element_class=CycleIndexSeries)

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.combinat.species.generating_series import CycleIndexSeriesRing
            sage: CycleIndexSeriesRing(QQ)
            Cycle Index Series Ring over Symmetric Functions over Rational Field in the powersum basis
        """
        return "Cycle Index Series Ring over %s"%self.base_ring()


class CycleIndexSeries(LazyPowerSeries):
    def count(self, t):
        """
        Return the number of structures corresponding to a certain cycle
        type ``t``.

        EXAMPLES::

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
        Returns the coefficient of a cycle type ``t`` in ``self``.

        EXAMPLES::

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


    def stretch(self, k):
        r"""
        Return the stretch of the cycle index series ``self`` by a positive
        integer `k`.

        If

        .. MATH::

           f = \sum_{n=0}^{\infty} f_n(p_1, p_2, p_3, \ldots ),

        then the stretch `g` of `f` by `k` is

        .. MATH::

           g = \sum_{n=0}^{\infty} f_n(p_k, p_{2k}, p_{3k}, \ldots ).

        EXAMPLES::

            sage: from sage.combinat.species.generating_series import CycleIndexSeriesRing
            sage: p = SymmetricFunctions(QQ).power()
            sage: CIS = CycleIndexSeriesRing(QQ)
            sage: f = CIS([p([]), p([1]), p([2]), p.zero()])
            sage: f.stretch(3).coefficients(10)
            [p[], 0, 0, p[3], 0, 0, p[6], 0, 0, 0]
        """
        return self._new(partial(self._stretch_gen, k), lambda ao: k*ao, self)

    def _stretch_gen(self, k, ao):
        """
        EXAMPLES::

            sage: from sage.combinat.species.generating_series import CycleIndexSeriesRing
            sage: p = SymmetricFunctions(QQ).power()
            sage: CIS = CycleIndexSeriesRing(QQ)
            sage: f = CIS([p([1])]) # This is the power series whose all coefficients
            ....:                   # are p[1]. Not combinatorially meaningful!
            sage: g = f._stretch_gen(2,0)
            sage: [next(g) for i in range(10)]
            [p[2], 0, p[2], 0, p[2], 0, p[2], 0, p[2], 0]
        """
        from sage.combinat.partition import Partition
        BR = self.base_ring()
        zero = BR.zero()

        stretch_k = lambda p: Partition([k*i for i in p])

        yield self.coefficient(0).map_support(stretch_k)

        n = 1
        while True:
            for i in range(k-1):
                yield zero
            yield self.coefficient(n).map_support(stretch_k)
            n += 1

    def isotype_generating_series(self):
        """
        EXAMPLES::

            sage: P = species.PermutationSpecies()
            sage: cis = P.cycle_index_series()
            sage: f = cis.isotype_generating_series()
            sage: f.coefficients(10)
            [1, 1, 2, 3, 5, 7, 11, 15, 22, 30]
        """
        R = self.base_ring().base_ring()
        OGS = OrdinaryGeneratingSeriesRing(R)()
        return OGS._new(self._ogs_gen, lambda ao: ao, self)

    def expand_as_sf(self, n, alphabet='x'):
        """
        Returns the expansion of a cycle index series as a symmetric function in
        ``n`` variables.

        Specifically, this returns a :class:`~sage.combinat.species.series.LazyPowerSeries` whose
        ith term is obtained by calling :meth:`~sage.combinat.sf.sfa.SymmetricFunctionAlgebra_generic_Element.expand`
        on the ith term of ``self``.

        This relies on the (standard) interpretation of a cycle index series as a symmetric function
        in the power sum basis.

        INPUT:

        - ``self`` -- a cycle index series

        - ``n`` -- a positive integer

        - ``alphabet`` -- a variable for the expansion (default: `x`)

        EXAMPLES::

            sage: from sage.combinat.species.set_species import SetSpecies
            sage: SetSpecies().cycle_index_series().expand_as_sf(2).coefficients(4)
            [1, x0 + x1, x0^2 + x0*x1 + x1^2, x0^3 + x0^2*x1 + x0*x1^2 + x1^3]

        """
        expanded_poly_ring = self.coefficient(0).expand(n, alphabet).parent()
        LPSR = LazyPowerSeriesRing(expanded_poly_ring)

        expander_gen = (LPSR.term(self.coefficient(i).expand(n, alphabet), i) for i in _integers_from(0))

        return LPSR.sum_generator(expander_gen)

    def _ogs_gen(self, ao):
        """
        Returns a generator for the coefficients of the ordinary generating
        series obtained from a cycle index series.

        EXAMPLES::

            sage: P = species.PermutationSpecies()
            sage: cis = P.cycle_index_series()
            sage: g = cis._ogs_gen(0)
            sage: [next(g) for i in range(10)]
            [1, 1, 2, 3, 5, 7, 11, 15, 22, 30]
        """
        for i in range(ao):
            yield 0
        for i in _integers_from(ao):
            yield sum( self.coefficient(i).coefficients() )

    def generating_series(self):
        """
        EXAMPLES::

            sage: P = species.PartitionSpecies()
            sage: cis = P.cycle_index_series()
            sage: f = cis.generating_series()
            sage: f.coefficients(5)
            [1, 1, 1, 5/6, 5/8]
        """
        R = self.base_ring().base_ring()
        EGS = ExponentialGeneratingSeriesRing(R)()
        return EGS._new(self._egs_gen, lambda ao: ao, self)

    def _egs_gen(self, ao):
        """
        Returns a generator for the coefficients of the exponential
        generating series obtained from a cycle index series.

        EXAMPLES::

            sage: P = species.PermutationSpecies()
            sage: cis = P.cycle_index_series()
            sage: g = cis._egs_gen(0)
            sage: [next(g) for i in range(10)]
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        """
        for i in range(ao):
            yield 0
        for i in _integers_from(ao):
            yield self.coefficient(i).coefficient([1]*i)

    def __invert__(self):
        r"""
        Return the multiplicative inverse of ``self``.

        This algorithm is derived from [BLL]_.

        EXAMPLES::

            sage: E = species.SetSpecies().cycle_index_series()
            sage: E.__invert__().coefficients(4)
            [p[], -p[1], 1/2*p[1, 1] - 1/2*p[2], -1/6*p[1, 1, 1] + 1/2*p[2, 1] - 1/3*p[3]]

        The defining characteristic of the multiplicative inverse `F^{-1}` of
        a cycle index series `F` is that `F \cdot F^{-1} = F^{-1} \cdot F = 1`
        (that is, both products with `F` yield the multiplicative identity `1`)::

            sage: E = species.SetSpecies().cycle_index_series()
            sage: (E * ~E).coefficients(6)
            [p[], 0, 0, 0, 0, 0]

        REFERENCES:

        - [BLL]_
        - [BLL-Intro]_
        - http://bergeron.math.uqam.ca/Site/bergeron_anglais_files/livre_combinatoire.pdf

        AUTHORS:

        - Andrew Gainer-Dewar
        """
        if self.coefficient(0) == 0:
            raise ValueError("constant term must be non-zero")

        def multinv_builder(i):
            return self.coefficient(0)**(-i-1) * (self.coefficient(0) + (-1)*self)**i

        return self.parent().sum_generator(multinv_builder(i) for i in _integers_from(0))

    def _div_(self, y):
        """
        TESTS::

            sage: E = species.SetSpecies().cycle_index_series()
            sage: (E / E).coefficients(6)
            [p[], 0, 0, 0, 0, 0]
        """
        return self*(~y)

    def functorial_composition(self, g):
        r"""
        Returns the functorial composition of ``self`` and ``g``.

        If `F` and `G` are species, their functorial composition is the species
        `F \Box G` obtained by setting `(F \Box G) [A] = F[ G[A] ]`.
        In other words, an `(F \Box G)`-structure on a set `A` of labels is an
        `F`-structure whose labels are the set of all `G`-structures on `A`.

        It can be shown (as in section 2.2 of [BLL]_) that there is a
        corresponding operation on cycle indices:

        .. MATH::

            Z_{F} \Box Z_{G} = \sum_{n \geq 0} \frac{1}{n!}
            \sum_{\sigma \in \mathfrak{S}_{n}}
            \operatorname{fix} F[ (G[\sigma])_{1}, (G[\sigma])_{2}, \ldots ]
            \, p_{1}^{\sigma_{1}} p_{2}^{\sigma_{2}} \cdots.

        This method implements that operation on cycle index series.

        EXAMPLES:

        The species `G` of simple graphs can be expressed in terms of a functorial
        composition: `G = \mathfrak{p} \Box \mathfrak{p}_{2}`, where
        `\mathfrak{p}` is the :class:`~sage.combinat.species.subset_species.SubsetSpecies`.
        This is how it is implemented in
        :meth:`~sage.combinat.species.library.SimpleGraphSpecies`::

            sage: S = species.SimpleGraphSpecies()
            sage: S.cycle_index_series().coefficients(5)
            [p[],
             p[1],
             p[1, 1] + p[2],
             4/3*p[1, 1, 1] + 2*p[2, 1] + 2/3*p[3],
             8/3*p[1, 1, 1, 1] + 4*p[2, 1, 1] + 2*p[2, 2] + 4/3*p[3, 1] + p[4]]
        """
        return self._new(partial(self._functorial_compose_gen, g), lambda a,b: 0, self, g)

    def _functorial_compose_gen(self, g, ao):
        """
        Return a generator for the coefficients of the functorial
        composition of ``self`` with ``g``.

        EXAMPLES::

            sage: E = species.SetSpecies()
            sage: E2 = species.SetSpecies(size=2)
            sage: WP = species.SubsetSpecies()
            sage: P2 = E2*E
            sage: P2_cis = P2.cycle_index_series()
            sage: WP_cis = WP.cycle_index_series()
            sage: g = WP_cis._functorial_compose_gen(P2_cis,0)
            sage: [next(g) for i in range(5)]
            [p[],
             p[1],
             p[1, 1] + p[2],
             4/3*p[1, 1, 1] + 2*p[2, 1] + 2/3*p[3],
             8/3*p[1, 1, 1, 1] + 4*p[2, 1, 1] + 2*p[2, 2] + 4/3*p[3, 1] + p[4]]
        """
        p = self.parent().base_ring()
        n = 0
        while True:
            res = p(0)
            for s in Partitions(n):
                t = g._cycle_type(s)
                q = self.count(t) / s.aut()
                res += q*p(s)
            yield res
            n += 1

    def arithmetic_product(self, g, check_input = True):
        r"""
        Return the arithmetic product of ``self`` with ``g``.

        For species `M` and `N` such that `M[\\varnothing] = N[\\varnothing] = \\varnothing`,
        their arithmetic product is the species `M \\boxdot N` of "`M`-assemblies of cloned `N`-structures".
        This operation is defined and several examples are given in [MM]_.

        The cycle index series for `M \\boxdot N` can be computed in terms of the component series `Z_M` and `Z_N`,
        as implemented in this method.

        INPUT:

        - ``g`` -- a cycle index series having the same parent as ``self``.

        - ``check_input`` -- (default: ``True``) a Boolean which, when set
          to ``False``, will cause input checks to be skipped.

        OUTPUT:

        The arithmetic product of ``self`` with ``g``. This is a cycle
        index series defined in terms of ``self`` and ``g`` such that
        if ``self`` and ``g`` are the cycle index series of two species
        `M` and `N`, their arithmetic product is the cycle index series
        of the species `M \\boxdot N`.

        EXAMPLES:

        For `C` the species of (oriented) cycles and `L_{+}` the species of nonempty linear orders, `C \\boxdot L_{+}` corresponds
        to the species of "regular octopuses"; a `(C \\boxdot L_{+})`-structure is a cycle of some length, each of whose elements
        is an ordered list of a length which is consistent for all the lists in the structure. ::

            sage: C = species.CycleSpecies().cycle_index_series()
            sage: Lplus = species.LinearOrderSpecies(min=1).cycle_index_series()
            sage: RegularOctopuses = C.arithmetic_product(Lplus)
            sage: RegOctSpeciesSeq = RegularOctopuses.generating_series().counts(8)
            sage: RegOctSpeciesSeq
            [0, 1, 3, 8, 42, 144, 1440, 5760]

        It is shown in [MM]_ that the exponential generating function for regular octopuses satisfies
        `(C \\boxdot L_{+}) (x) = \\sum_{n \geq 1} \\sigma (n) (n - 1)! \\frac{x^{n}}{n!}` (where `\\sigma (n)` is
        the sum of the divisors of `n`). ::

            sage: RegOctDirectSeq = [0] + [sum(divisors(i))*factorial(i-1) for i in range(1,8)]
            sage: RegOctDirectSeq == RegOctSpeciesSeq
            True

        AUTHORS:

        - Andrew Gainer-Dewar (2013)

        REFERENCES:

        .. [MM] \M. Maia and M. Mendez. "On the arithmetic product of combinatorial species".
           Discrete Mathematics, vol. 308, issue 23, 2008, pp. 5407-5427.
           :arxiv:`math/0503436v2`.

        """
        from itertools import product, repeat, chain

        p = self.base_ring()

        if check_input:
            assert self.coefficient(0) == p.zero()
            assert g.coefficient(0) == p.zero()

        # We first define an operation `\\boxtimes` on partitions as in Lemma 2.1 of [MM]_.
        def arith_prod_of_partitions(l1, l2):
            # Given two partitions `l_1` and `l_2`, we construct a new partition `l_1 \\boxtimes l_2` by
            # the following procedure: each pair of parts `a \\in l_1` and `b \\in l_2` contributes
            # `\\gcd (a, b)`` parts of size `\\lcm (a, b)` to `l_1 \\boxtimes l_2`. If `l_1` and `l_2`
            # are partitions of integers `n` and `m`, respectively, then `l_1 \\boxtimes l_2` is a
            # partition of `nm`.
            term_iterable = chain.from_iterable(repeat(lcm(pair), gcd(pair))
                                                for pair in product(l1, l2))
            return Partition(sorted(term_iterable, reverse=True))

        # We then extend this to an operation on symmetric functions as per eq. (52) of [MM]_.
        # (Maia and Mendez, in [MM]_, are talking about polynomials instead of symmetric
        # functions, but this boils down to the same: Their x_i corresponds to the i-th power
        # sum symmetric function.)
        def arith_prod_sf(x, y):
            ap_sf_wrapper = lambda l1, l2: p(arith_prod_of_partitions(l1, l2))
            return p._apply_multi_module_morphism(x, y, ap_sf_wrapper)

        # Sage stores cycle index series by degree.
        # Thus, to compute the arithmetic product `Z_M \\boxdot Z_N` it is useful
        # to compute all terms of a given degree `n` at once.
        def arith_prod_coeff(n):
            if n == 0:
                res = p.zero()
            else:
                index_set = ((d, n // d) for d in divisors(n))
                res = sum(arith_prod_sf(self.coefficient(i), g.coefficient(j)) for i,j in index_set)

            # Build a list which has res in the `n`th slot and 0's before and after
            # to feed to sum_generator
            res_in_seq = [p.zero()]*n + [res, p.zero()]

            return self.parent(res_in_seq)

        # Finally, we use the sum_generator method to assemble these results into a single
        # LazyPowerSeries object.
        return self.parent().sum_generator(arith_prod_coeff(n) for n in _integers_from(0))

    def _cycle_type(self, s):
        """
        EXAMPLES::

            sage: cis = species.PartitionSpecies().cycle_index_series()
            sage: [cis._cycle_type(p) for p in Partitions(3)]
            [[3, 1, 1], [2, 1, 1, 1], [1, 1, 1, 1, 1]]
            sage: cis = species.PermutationSpecies().cycle_index_series()
            sage: [cis._cycle_type(p) for p in Partitions(3)]
            [[3, 1, 1, 1], [2, 2, 1, 1], [1, 1, 1, 1, 1, 1]]
            sage: cis = species.SetSpecies().cycle_index_series()
            sage: [cis._cycle_type(p) for p in Partitions(3)]
            [[1], [1], [1]]
        """
        if s == []:
            return self._card(0)
        res = []
        for k in range(1, self._upper_bound_for_longest_cycle(s)+1):
            e = 0
            for d in divisors(k):
                m = moebius(d)
                if m == 0:
                    continue
                u = s.power(k/d)
                e += m*self.count(u)
            res.extend([k]*int(e/k))
        res.reverse()
        return Partition(res)


    def _upper_bound_for_longest_cycle(self, s):
        """
        EXAMPLES::

            sage: cis = species.PartitionSpecies().cycle_index_series()
            sage: cis._upper_bound_for_longest_cycle([4])
            4
            sage: cis._upper_bound_for_longest_cycle([3,1])
            3
            sage: cis._upper_bound_for_longest_cycle([2,2])
            2
            sage: cis._upper_bound_for_longest_cycle([2,1,1])
            2
            sage: cis._upper_bound_for_longest_cycle([1,1,1,1])
            1
        """
        if s == []:
            return 1
        return min(self._card(sum(s)), lcm(list(s)))

    def _card(self, n):
        r"""
        Return the number of structures on an underlying set of size ``n`` for
        the species associated with ``self``.

        This is just ``n!`` times the coefficient of ``p[1]n`` in ``self``.

        EXAMPLES::

            sage: cis = species.PartitionSpecies().cycle_index_series()
            sage: cis._card(4)
            15
        """
        p = self.coefficient(n)
        return factorial(n) * p.coefficient([1] * n)

    def _compose_gen(self, y, ao):
        """
        Return a generator for the coefficients of the composition of this
        cycle index series and the cycle index series ``y``. This overrides
        the method defined in ``LazyPowerSeries``.

        The notion "composition" means plethystic substitution here, as
        defined in Section 2.2 of [BLL-Intro]_.

        EXAMPLES::

            sage: E = species.SetSpecies(); C = species.CycleSpecies()
            sage: E_cis = E.cycle_index_series()
            sage: g = E_cis._compose_gen(C.cycle_index_series(),0)
            sage: [next(g) for i in range(4)]
            [p[], p[1], p[1, 1] + p[2], p[1, 1, 1] + p[2, 1] + p[3]]
        """
        assert y.coefficient(0) == 0
        y_powers = Stream(y._power_gen())

        parent = self.parent()
        res =  parent.sum_generator(self._compose_term(self.coefficient(i), y_powers)
                                    for i in _integers_from(0))

        for i in _integers_from(0):
            yield res.coefficient(i)

    def _compose_term(self, p, y_powers):
        r"""
        Return the composition of one term in ``self`` with `y`.

        INPUT:

        -  ``p`` - a term in ``self``
        -  ``y_powers`` - a stream for the powers of `y`
           starting with `y`

        EXAMPLES::

            sage: from sage.combinat.species.stream import Stream
            sage: E = species.SetSpecies(); C = species.CycleSpecies()
            sage: E_cis = E.cycle_index_series()
            sage: C_cis = C.cycle_index_series()
            sage: c_powers = Stream(C_cis._power_gen())
            sage: p2 = E_cis.coefficient(2); p2
            1/2*p[1, 1] + 1/2*p[2]
            sage: E_cis._compose_term(p2, c_powers).coefficients(4)
            [0, 0, 1/2*p[1, 1] + 1/2*p[2], 1/2*p[1, 1, 1] + 1/2*p[2, 1]]
        """
        parent = self.parent()
        if p == 0:
            return parent(0)

        res = []
        #Go through all the partition, coefficient pairs in the term p
        for m, c in p:
            res_t = parent.term(c, 0)

            for e,v in enumerate(m.to_exp()):
                if v == 0:
                    continue
                res_t = res_t * y_powers[v-1].stretch(e+1)
            res.append(res_t)

        return parent.sum(res)

    def weighted_composition(self, y_species):
        r"""
        Return the composition of this cycle index series with the cycle
        index series of the weighted species ``y_species``.

        Note that this is basically the same algorithm as composition
        except we can not use the optimization that the powering of cycle
        index series commutes with 'stretching'.

        EXAMPLES::

            sage: E = species.SetSpecies(); C = species.CycleSpecies()
            sage: E_cis = E.cycle_index_series()
            sage: E_cis.weighted_composition(C).coefficients(4)
            [p[], p[1], p[1, 1] + p[2], p[1, 1, 1] + p[2, 1] + p[3]]
            sage: E(C).cycle_index_series().coefficients(4)
            [p[], p[1], p[1, 1] + p[2], p[1, 1, 1] + p[2, 1] + p[3]]
        """
        base_ring = self.base_ring()
        y = y_species.cycle_index_series(base_ring)
        assert y.coefficient(0) == 0
        return self._new(partial(self._weighted_compose_gen, y_species), lambda a,b:a*b, self, y)


    def _weighted_compose_gen(self, y_species, ao):
        r"""
        Return an iterator for the composition of this cycle index series
        and the cycle index series of the weighted species ``y_species``.

        EXAMPLES::

            sage: E = species.SetSpecies(); C = species.CycleSpecies()
            sage: E_cis = E.cycle_index_series()
            sage: g = E_cis._weighted_compose_gen(C,0)
            sage: [next(g) for i in range(4)]
            [p[], p[1], p[1, 1] + p[2], p[1, 1, 1] + p[2, 1] + p[3]]
        """
        parent = self.parent()
        res =  parent.sum_generator(self._weighted_compose_term(self.coefficient(i), y_species)
                                    for i in _integers_from(0))

        for i in _integers_from(0):
            yield res.coefficient(i)

    def _weighted_compose_term(self, p, y_species):
        r"""
        Return the weighted composition of one term in ``self`` with ``y``.

        INPUT:

        -  ``p`` -- a term in ``self``
        -  ``y_species`` -- a species

        EXAMPLES::

            sage: E = species.SetSpecies(); C = species.CycleSpecies()
            sage: E_cis = E.cycle_index_series()
            sage: p2 = E_cis.coefficient(2); p2
            1/2*p[1, 1] + 1/2*p[2]
            sage: E_cis._weighted_compose_term(p2, C).coefficients(4)
            [0, 0, 1/2*p[1, 1] + 1/2*p[2], 1/2*p[1, 1, 1] + 1/2*p[2, 1]]
        """
        parent = self.parent()
        if p == 0:
            return parent(0)

        base_ring = self.base_ring().base_ring()

        res = []
        #Go through all the partition, coefficient pairs in the term p
        for m, c in p:
            res_t = parent.term(c, 0)

            for e,v in enumerate(m.to_exp()):
                if v == 0:
                    continue
                res_t = res_t * (y_species.weighted(y_species._weight**(e+1)).cycle_index_series(base_ring)**v).stretch(e+1)
            res.append(res_t)

        return parent.sum(res)

    def compositional_inverse(self):
        r"""
        Return the compositional inverse of ``self`` if possible.

        (Specifically, if ``self`` is of the form `0 + p_{1} + \cdots`.)

        The compositional inverse is the inverse with respect to
        plethystic substitution. This is the operation on cycle index
        series which corresponds to substitution, a.k.a. partitional
        composition, on the level of species. See Section 2.2 of
        [BLL]_ for a definition of this operation.

        EXAMPLES::

            sage: Eplus = species.SetSpecies(min=1).cycle_index_series()
            sage: Eplus(Eplus.compositional_inverse()).coefficients(8)
            [0, p[1], 0, 0, 0, 0, 0, 0]

        TESTS::

            sage: Eplus = species.SetSpecies(min=2).cycle_index_series()
            sage: Eplus.compositional_inverse()
            Traceback (most recent call last):
            ...
            ValueError: not an invertible series

        ALGORITHM:

        Let `F` be a species satisfying `F = 0 + X + F_2 + F_3 + \cdots` for
        `X` the species of singletons. (Equivalently, `\lvert F[\varnothing]
        \rvert = 0` and `\lvert F[\{1\}] \rvert = 1`.) Then there exists a
        (virtual) species `G` satisfying `F \circ G = G \circ F = X`.

        It follows that `(F - X) \circ G = F \circ G - X \circ G = X - G`.
        Rearranging, we obtain the recursive equation `G = X - (F - X) \circ G`,
        which can be solved using iterative methods.

        .. WARNING::

            This algorithm is functional but can be very slow.
            Use with caution!

        .. SEEALSO::

            The compositional inverse `\Omega` of the species `E_{+}`
            of nonempty sets can be handled much more efficiently
            using specialized methods. See
            :func:`~sage.combinat.species.generating_series.LogarithmCycleIndexSeries`

        AUTHORS:

        - Andrew Gainer-Dewar
        """
        cisr = self.parent()
        sfa = cisr._base

        X = cisr([0, sfa([1]), 0])

        if self.coefficients(2) != X.coefficients(2):
            raise ValueError('not an invertible series')

        res = cisr()
        res.define(X - (self - X).compose(res))

        return res

    def derivative(self, order=1):
        r"""
        Return the species-theoretic `n`-th derivative of ``self``,
        where `n` is ``order``.

        For a cycle index series `F (p_{1}, p_{2}, p_{3}, \ldots)`, its
        derivative is the cycle index series `F' = D_{p_{1}} F` (that is,
        the formal derivative of `F` with respect to the variable `p_{1}`).

        If `F` is the cycle index series of a species `S` then `F'` is the
        cycle index series of an associated species `S'` of `S`-structures
        with a "hole".

        EXAMPLES:

        The species `E` of sets satisfies the relationship `E' = E`::

            sage: E = species.SetSpecies().cycle_index_series()
            sage: E.coefficients(8) == E.derivative().coefficients(8)
            True

        The species `C` of cyclic orderings and the species `L` of linear
        orderings satisfy the relationship `C' = L`::

            sage: C = species.CycleSpecies().cycle_index_series()
            sage: L = species.LinearOrderSpecies().cycle_index_series()
            sage: L.coefficients(8) == C.derivative().coefficients(8)
            True
        """
        # Make sure that order is integral
        order = Integer(order)

        if order < 0:
            raise ValueError("Order must be a non-negative integer")

        elif order == 0:
            return self

        elif order == 1:
            parent = self.parent()
            derivative_term = lambda n: parent.term(self.coefficient(n+1).derivative_with_respect_to_p1(), n)
            return parent.sum_generator(derivative_term(i) for i in _integers_from(0))

        else:
            return self.derivative(order-1)

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

            sage: E = species.SetSpecies().cycle_index_series()
            sage: X = species.SingletonSpecies().cycle_index_series()
            sage: E.pointing().coefficients(8) == (X*E).coefficients(8)
            True

        """
        p1 = self.base_ring()([1])
        X = self.parent()([0, p1, 0])

        return X*self.derivative()

    def integral(self, *args):
        r"""
        Given a cycle index `G`, it is not in general possible to recover a
        single cycle index `F` such that `F' = G` (even up to addition of a
        constant term).

        More broadly, it may be the case that there are many non-isomorphic
        species `S` such that `S' = T` for a given species `T`.
        For example, the species `3 C_{3}` of 3-cycles from three distinct
        classes and the species `X^{3}` of 3-sets are not isomorphic, but
        `(3 C_{3})' = (X^{3})' = 3 X^{2}`.

        EXAMPLES::

            sage: C3 = species.CycleSpecies(size=3).cycle_index_series()
            sage: X = species.SingletonSpecies().cycle_index_series()
            sage: (3*C3).derivative().coefficients(8) == (3*X^2).coefficients(8)
            True
            sage: (X^3).derivative().coefficients(8) == (3*X^2).coefficients(8)
            True

        .. WARNING::

            This method has no implementation and exists only to prevent you from
            doing something strange. Calling it raises a ``NotImplementedError``!
        """
        raise NotImplementedError

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

            sage: BT = species.BinaryTreeSpecies().cycle_index_series()
            sage: BF = species.BinaryForestSpecies().cycle_index_series()
            sage: BT.exponential().isotype_generating_series().coefficients(8) == BF.isotype_generating_series().coefficients(8)
            True
        """
        base_ring = self.parent().base_ring().base_ring()
        E = ExponentialCycleIndexSeries(base_ring)
        return E.compose(self)

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

            sage: G = species.SimpleGraphSpecies().cycle_index_series() - 1
            sage: from sage.combinat.species.generating_series import LogarithmCycleIndexSeries
            sage: CG = LogarithmCycleIndexSeries().compose(G)
            sage: CG.isotype_generating_series().coefficients(8)
            [0, 1, 1, 2, 6, 21, 112, 853]
        """
        base_ring = self.parent().base_ring().base_ring()
        Omega = LogarithmCycleIndexSeries(base_ring)
        return Omega.compose(self)

@cached_function
def _exp_term(n, R = RationalField()):
    """
    Compute the order-n term of the cycle index series of the species `E` of sets.

    EXAMPLES::

        sage: from sage.combinat.species.generating_series import _exp_term
        sage: [_exp_term(i) for i in range(4)]
        [p[], p[1], 1/2*p[1, 1] + 1/2*p[2], 1/6*p[1, 1, 1] + 1/2*p[2, 1] + 1/3*p[3]]
    """
    p = SymmetricFunctions(R).power()
    return sum(p(part) / part.aut() for part in Partitions(n))


def _exp_gen(R = RationalField()):
    r"""
    Produce a generator which yields the terms of the cycle index
    series of the species `E` of sets.

    EXAMPLES::

        sage: from sage.combinat.species.generating_series import _exp_gen
        sage: g = _exp_gen()
        sage: [next(g) for i in range(4)]
        [p[], p[1], 1/2*p[1, 1] + 1/2*p[2], 1/6*p[1, 1, 1] + 1/2*p[2, 1] + 1/3*p[3]]
    """
    return (_exp_term(i, R) for i in _integers_from(0))

@cached_function
def ExponentialCycleIndexSeries(R = RationalField()):
    r"""
    Return the cycle index series of the species `E` of sets.

    This cycle index satisfies

    .. MATH::

        Z_{E} = \sum_{n \geq 0} \sum_{\lambda \vdash n}
        \frac{p_{\lambda}}{z_{\lambda}}.

    EXAMPLES::

        sage: from sage.combinat.species.generating_series import ExponentialCycleIndexSeries
        sage: ExponentialCycleIndexSeries().coefficients(5)
        [p[], p[1], 1/2*p[1, 1] + 1/2*p[2], 1/6*p[1, 1, 1] + 1/2*p[2, 1]
         + 1/3*p[3], 1/24*p[1, 1, 1, 1] + 1/4*p[2, 1, 1] + 1/8*p[2, 2]
         + 1/3*p[3, 1] + 1/4*p[4]]
    """
    CIS = CycleIndexSeriesRing(R)
    return CIS(_exp_gen(R))


@cached_function
def _cl_term(n, R = RationalField()):
    r"""
    Compute the order-n term of the cycle index series of the virtual species
    `\Omega`, the compositional inverse of the species `E^{+}` of nonempty sets.

    EXAMPLES::

        sage: from sage.combinat.species.generating_series import _cl_term
        sage: [_cl_term(i) for i in range(4)]
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


def _cl_gen (R = RationalField()):
    r"""
    Produce a generator which yields the terms of the cycle index series
    of the virtual species `\Omega`, the compositional inverse of the
    species `E^{+}` of nonempty sets.

    EXAMPLES::

        sage: from sage.combinat.species.generating_series import _cl_gen
        sage: g = _cl_gen()
        sage: [next(g) for i in range(4)]
        [0, p[1], -1/2*p[1, 1] - 1/2*p[2], 1/3*p[1, 1, 1] - 1/3*p[3]]
    """
    return (_cl_term(i, R) for i in _integers_from(0))


@cached_function
def LogarithmCycleIndexSeries(R = RationalField()):
    r"""
    Return the cycle index series of the virtual species `\Omega`, the
    compositional inverse of the species `E^{+}` of nonempty sets.

    The notion of virtual species is treated thoroughly in [BLL]_.
    The specific algorithm used here to compute the cycle index of
    `\Omega` is found in [Labelle2008]_.

    EXAMPLES:

    The virtual species `\Omega` is 'properly virtual', in the sense that
    its cycle index has negative coefficients::

        sage: from sage.combinat.species.generating_series import LogarithmCycleIndexSeries
        sage: LogarithmCycleIndexSeries().coefficients(4)
        [0, p[1], -1/2*p[1, 1] - 1/2*p[2], 1/3*p[1, 1, 1] - 1/3*p[3]]

    Its defining property is that `\Omega \circ E^{+} = E^{+} \circ \Omega = X`
    (that is, that composition with `E^{+}` in both directions yields the
    multiplicative identity `X`)::

        sage: Eplus = sage.combinat.species.set_species.SetSpecies(min=1).cycle_index_series()
        sage: LogarithmCycleIndexSeries().compose(Eplus).coefficients(4)
        [0, p[1], 0, 0]
    """
    CIS = CycleIndexSeriesRing(R)
    return CIS(_cl_gen(R))

