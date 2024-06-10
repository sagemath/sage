r"""
Group Cycle Indices

This file implements the group cycle indices of Henderson and Gainer-Dewar.

For a group `\Gamma` and a ring `R`, a `\Gamma`-cycle index over `R`
is an element of the free module over the ring of cycle index series
over `R` with basis `\Gamma`.

In other words, a `\Gamma`-cycle index series over `R` is a formal sum

.. MATH::

    F = \sum_{\gamma \in \Gamma} F[\gamma] \cdot \gamma,

where each coefficient `F[\gamma]` is a cycle index series over `R`.
(By functoriality, if `F` is the `\Gamma`-cycle index of a
`\Gamma`-species, it must be a class function of `\Gamma`.)

These objects are of interest because they can be used to enumerate `\Gamma`-species;
they serve the same role in that theory as ordinary cycle indices do for classical
species.

Just like classical species, `\Gamma`-species may be combined in ways
which correspond naturally to algebraic operations on their
`\Gamma`-cycle indices.  Many standard operations on
`\Gamma`-species are available, including addition, multiplication,
and differentiation.  The associated operations `+`, `\cdot`, and
`\prime` on `\Gamma`-cycle indices are defined componentwise and are
implemented by this class.

Additionally, there is a natural way to compose `\Gamma`-species,
corresponding to the construction of
:class:`~sage.combinat.species.composition_species.CompositionSpecies`.
The `\Gamma`-cycle index of such a composition can be computed using
a variant of plethysm, but this depends on having access to all the
terms of both component `\Gamma`-cycle indices; it cannot be computed
componentwise.  More details of this operation are available in the
documentation for :meth:`~sage.combinat.species.group_cycle_index_series.GroupCycleIndexSeries.composition`,
which implements it in Sage.

AUTHORS:

- Andrew Gainer-Dewar (2013): initial version

EXAMPLES::

    sage: from sage.combinat.species.group_cycle_index_series import GroupCycleIndexSeriesRing
    sage: GCISR = GroupCycleIndexSeriesRing(SymmetricGroup(4))
    sage: loads(dumps(GCISR))
    Ring of (Symmetric group of order 4! as a permutation group)-Cycle Index Series over Rational Field

.. TODO::

    Implement (optional?) optimizations using assumption that a
    `\Gamma`-cycle index is a class function of `\Gamma`.

"""
#*****************************************************************************
#       Copyright (C) 2013 Andrew Gainer-Dewar <andrew.gainer.dewar@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.rings.rational_field import RationalField
from sage.misc.cachefunc import cached_function
from sage.combinat.free_module import CombinatorialFreeModule

@cached_function
def GroupCycleIndexSeriesRing(G, R = RationalField()):
    r"""
    Return the ring of group cycle index series.

    EXAMPLES::

        sage: from sage.combinat.species.group_cycle_index_series import GroupCycleIndexSeriesRing
        sage: GCISR = GroupCycleIndexSeriesRing(SymmetricGroup(4)); GCISR
        Ring of (Symmetric group of order 4! as a permutation group)-Cycle Index Series over Rational Field

    TESTS: We test to make sure that caching works. ::

        sage: from sage.combinat.species.group_cycle_index_series import GroupCycleIndexSeriesRing
        sage: GCISR is GroupCycleIndexSeriesRing(SymmetricGroup(4))
        True
    """
    return GroupCycleIndexSeriesRing_class(G, R)

class GroupCycleIndexSeriesRing_class(CombinatorialFreeModule):
    def __init__(self, G, R = RationalField()):
        r"""
        EXAMPLES::

            sage: from sage.combinat.species.group_cycle_index_series import GroupCycleIndexSeriesRing
            sage: GCISR = GroupCycleIndexSeriesRing(SymmetricGroup(4)); GCISR
            Ring of (Symmetric group of order 4! as a permutation group)-Cycle Index Series over Rational Field
            sage: GCISR == loads(dumps(GCISR))
            True
        """
        from sage.combinat.species.generating_series import CycleIndexSeriesRing
        from sage.categories.algebras_with_basis import AlgebrasWithBasis

        self._coeff_ring = R
        CISR = CycleIndexSeriesRing(R)
        self._cisr = CISR
        self._group = G

        CombinatorialFreeModule.__init__(self, CISR, G, element_class = GroupCycleIndexSeries, category = AlgebrasWithBasis(CISR), prefix = 'G')

    def product_on_basis(self, left, right):
        r"""
        Return the product of two basis elements ``left`` and ``right`` of ``self``.

        Multiplication of `\Gamma`-cycle indices is defined componentwise.
        That is, if `F` and `G` are two `\Gamma`-cycle indices, then `(F \cdot G) [\gamma] = F [\gamma] \cdot G [\gamma]`,
        where the multiplication on the right-hand side is ordinary multiplication of cycle indices.

        This is handled in Sage by defining multiplication on the basis of monomials induced
        by elements of `\Gamma`.

        EXAMPLES::

            sage: from sage.combinat.species.group_cycle_index_series import GroupCycleIndexSeriesRing
            sage: GCISR = GroupCycleIndexSeriesRing(SymmetricGroup(4))
            sage: e = SymmetricGroup(4).identity()
            sage: t = SymmetricGroup(4)([4,3,2,1])
            sage: GCISR.product_on_basis(t,t) == GCISR(t)
            True

        """
        if left == right:
            return self.monomial(left)

        return self(0)

    def one(self):
        r"""
        Return the multiplicative identity element of this algebra.

        EXAMPLES::

            sage: from sage.combinat.species.group_cycle_index_series import GroupCycleIndexSeriesRing
            sage: GCISR = GroupCycleIndexSeriesRing(SymmetricGroup(2))
            sage: GCISR.one()
            p[]*G[()] + p[]*G[(1,2)]

        """
        basis = self.basis()
        return self.sum(basis[k] for k in basis.keys())

    def _an_element_(self):
        r"""
        Return a representative element of this algebra.

        EXAMPLES::

            sage: from sage.combinat.species.group_cycle_index_series import GroupCycleIndexSeriesRing
            sage: GCISR = GroupCycleIndexSeriesRing(SymmetricGroup(2))
            sage: GCISR.an_element()
            p[]*G[()] + p[]*G[(1,2)]
        """
        G = self.basis().keys()
        return self.monomial(G.one()) + self.monomial(G.an_element())

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.species.group_cycle_index_series import GroupCycleIndexSeriesRing
            sage: GroupCycleIndexSeriesRing(SymmetricGroup(4))
            Ring of (Symmetric group of order 4! as a permutation group)-Cycle Index Series over Rational Field
        """
        return "Ring of (%s)-Cycle Index Series over %s" %(self._group, self._coeff_ring)

class GroupCycleIndexSeries(CombinatorialFreeModule.Element):
    r"""

    EXAMPLES::

        sage: from sage.combinat.species.group_cycle_index_series import GroupCycleIndexSeriesRing
        sage: GCISR = GroupCycleIndexSeriesRing(SymmetricGroup(4))
        sage: GCISR.an_element()
        p[]*G[()] + p[]*G[(2,3,4)]

    """

    def quotient(self):
        r"""
        Return the quotient of this group cycle index.

        This is defined to be the ordinary cycle index `F / \Gamma` obtained from a
        `\Gamma`-cycle index `F` by:

        .. MATH::
            F / \Gamma = \\frac{1}{\\lvert \Gamma \\rvert} \sum_{\gamma \in \Gamma} F [\gamma].

        It is shown in [AGdiss]_ that, if `F` is the `\Gamma`-cycle
        index of a `\Gamma`-species, `F / \Gamma` is the ordinary
        cycle index of orbits of structures under the action of
        `\Gamma`.  This corresponds to the notion of quotient species
        introduced in [Bousquet]_.

        EXAMPLES::

            sage: from sage.combinat.species.group_cycle_index_series import GroupCycleIndexSeriesRing
            sage: S4 = SymmetricGroup(4)
            sage: GCISR = GroupCycleIndexSeriesRing(S4)
            sage: GCISR.an_element()
            p[]*G[()] + p[]*G[(2,3,4)]
            sage: GCISR.an_element().quotient().coefficients(4)
            [1/12*p[], 0, 0, 0]

        REFERENCES:

        .. [AGdiss] Andrew Gainer. "`\Gamma`-species, quotients, and graph enumeration". Ph.D. diss. Brandeis University, 2012.
        .. [Bousquet] Michel Bousquet. "Especes de structures et applications du denombrement de cartes et de cactus planaires". Ph.D. diss. Universite du Quebec a Montreal, 1999.
           http://lacim.uqam.ca/publications_pdf/24.pdf

        """
        return 1/self.parent()._group.cardinality() * sum(self.coefficients())

    def composition(self, y, test_constant_term_is_zero = True):
        r"""
        Return the group-cycle-index plethysm of ``self`` with ``y``.

        Plethysm of group cycle index series is defined by a sort of 'mixing' operation in [Hend]_:

        .. MATH::
            (F \circ G) [\gamma] (p_{1}, p_{2}, p_{3}, \dots) =
            F [\gamma] \left( G [\gamma] (p_{1}, p_{2}, p_{3}, \dots),
            G [\gamma^{2}] (p_{2}, p_{4}, p_{6}, \dots), \dots \\right).

        It is shown in [Hend]_ that this operation on $\Gamma$-cycle
        indices corresponds to the 'composition' operation on
        $\Gamma$-species, a natural analogue of the
        :class:`~sage.combinat.species.composition_species.CompositionSpecies`
        operation on ordinary species.

        It is required that each of the components of `y` has zero constant term.
        However, testing this may not be possible if `y` is of an exotic class.
        Set ``test_constant_term_is_zero`` to ``False`` to suppress any testing.

        EXAMPLES::

            sage: from sage.combinat.species.group_cycle_index_series import GroupCycleIndexSeriesRing
            sage: from sage.combinat.species.group_cycle_index_series_library import CyclicOrderWithReversalGroupCycleIndex
            sage: S2 = SymmetricGroup(2)
            sage: GCISR = GroupCycleIndexSeriesRing(S2)
            sage: E = sage.combinat.species.set_species.SetSpecies().cycle_index_series()
            sage: C = CyclicOrderWithReversalGroupCycleIndex()
            sage: example = (E*GCISR.one())(C)
            sage: example.quotient().generating_series().counts(8)
            [1, 1, 2, 5, 17, 73, 398, 2636]
            sage: x = oeis("A007868") # optional: oeis
            sage: x.first_terms(8) # optional: oeis
            (1, 1, 2, 5, 17, 73, 398, 2636)
            sage: x.name() # optional: oeis
            'Number of inverse pairs of elements in symmetric group S_n, or pairs of total orders on n nodes (average of A000085 and A000142).'
            sage: example[S2.identity()].coefficients(4)
            [p[], p[1], p[1, 1] + p[2], p[1, 1, 1] + p[2, 1] + p[3]]
            sage: example[S2.gen()].coefficients(4)
            [p[], p[1], p[1, 1] + p[2], 2/3*p[1, 1, 1] + 2*p[2, 1] + 1/3*p[3]]

            sage: example[S2.gen()].generating_series().counts(8)
            [1, 1, 2, 4, 10, 26, 76, 232]
            sage: x = oeis("A000085") # optional: oeis
            sage: x.first_terms(8) # optional: oeis
            (1, 1, 2, 4, 10, 26, 76, 232)
            sage: x.name() # optional: oeis
            'Number of self-inverse permutations on n letters, also known as involutions; number of Young tableaux with n cells.'

        REFERENCES:

        .. [Hend] Anthony Henderson. "Species over a finite field". J. Algebraic Combin., 21(2):147-161, 2005.

        """
        from sage.combinat.species.stream import Stream,_integers_from
        from sage.misc.misc_c import prod

        assert self.parent() == y.parent()

        parent = self.parent()
        cisr = parent._cisr
        group = parent._group

        if test_constant_term_is_zero:
            for ycis in y.coefficients():
                assert ycis.coefficient(0) == 0

        ypowers = { g: Stream(y[g]._power_gen()) for g in group }

        def monomial_composer( partition, g ):
            partmults = ( (i+1, j) for i,j in enumerate(partition.to_exp()) if j != 0 )
            ycontrib = lambda part, mult: ypowers[g**part][mult-1].stretch(part)
            res = cisr.one()*prod(ycontrib(part, mult) for part,mult in partmults)
            return res

        def term_map( term, g ):
            if test_constant_term_is_zero and term == 0:
                return cisr.zero()

            res = sum(coeff*monomial_composer(part, g) for part,coeff in term)
            return res

        def component_builder( g ):
            if test_constant_term_is_zero and self[g] == 0:
                res = cisr.zero()
            elif g == group.identity(): #Save time by using existing CIS composition code
                res = self[g](y[g])
            else:
                res = cisr.sum_generator(term_map(self[g].coefficient(i), g) for i in _integers_from(0))
            return res

        components_dict = { g: component_builder(g) for g in group }
        res = parent._from_dict(components_dict)
        return res

    __call__ = composition

    def functorial_composition(self, y):
        r"""
        Return the functorial composition of ``self`` with ``y``.

        Functorial composition of group cycle index series satisfies

        .. MATH::
          (F \\square G) [\\gamma] = \\sum_{n \\geq 0} \\frac{1}{n!} \\sum_{\\sigma \in \\mathfrak{S}_{n}}
             \\operatorname{fix} \\left(\gamma \\cdot F \\left[ \\gamma \\cdot G [\\sigma] \\right] \\right).

        This operation on `\Gamma`-cycle indices corresponds to the functorial composition
        operation on `\Gamma`-species. A formula for the permutation `\gamma \cdot G [\sigma]`
        is given in [AGDPolya]_.

        EXAMPLES:

        Consider the `\mathfrak{S}_{2}`-species `\mathcal{D}` of directed graphs, where the action
        of the nontrivial element of `\mathfrak{S}_{2}` reverses the directions of all the edges.

        Let `\mathcal{E}` be the species of sets with the trivial action of `\mathfrak{S}_{2}`,
        `\mathcal{P}` the species of subsets with the trivial action, and `\mathcal{L}_{2}` the species
        of linear `2`-orders with the order-reversing `\mathfrak{S}_{2}`-action as defined in
        :meth:`~sage.combinat.species.group_cycle_index_series_library.LinearOrderWithReversalGroupCycleIndex`.::

            sage: from sage.combinat.species.group_cycle_index_series import GroupCycleIndexSeriesRing
            sage: from sage.combinat.species.library import SetSpecies, SubsetSpecies
            sage: from sage.combinat.species.group_cycle_index_series_library import LinearOrderWithReversalGroupCycleIndex
            sage: S2 = SymmetricGroup(2)
            sage: GCISR = GroupCycleIndexSeriesRing(S2)
            sage: P = GCISR(SubsetSpecies().cycle_index_series())
            sage: E = GCISR(SetSpecies().cycle_index_series())
            sage: L2 = LinearOrderWithReversalGroupCycleIndex().restricted(min=2,max=3)
            sage: D = P.functorial_composition(L2*E)
            sage: D.quotient().isotype_generating_series().counts(5)
            [1, 1, 3, 13, 144]

        (Compare :oeis:`A054933`.)

        REFERENCES:

        .. [AGDPolya] Andrew Gainer-Dewar. "Species with an equivariant group action." Preprint, :arxiv:`1401.6202`.

        """
        from sage.combinat.species.stream import _integers_from
        from sage.combinat.partition import Partition, Partitions
        from sage.arith.misc import divisors, moebius, factorial
        from sage.rings.integer import Integer

        assert self.parent() == y.parent()

        parent = self.parent()
        group = parent._group
        cisr = parent.base_ring()
        p = cisr.base_ring()

        # Compute the number of k-cycles in g*y[sigma] if sigma has cycle type partition ctp
        def cyccount(k, g, ctp):
            # Helper function for cycle computation
            sumterm = lambda d: moebius(k / d) * y[g**d].coefficient_cycle_type(ctp.power(d)) * ctp.power(d).aut()
            result = sum(sumterm(d) for d in divisors(k)) / k
            return Integer(result)

        # Find the cycle type of g*y[sigma] if sigma has cycle type partition ctp
        def g_ctp_builder(g, ctp):
            n = ctp.size()
            scount = y[group.one()].coefficient_cycle_type([1]*n) * factorial(n) # Number of y-structures with n labels

            result = []
            k = 1

            while sum(result) < scount and k <= scount:
                result = [k]*cyccount(k, g, ctp) + result
                k += 1

            assert sum(result) == scount

            return Partition(result)

        # Compute the g*ctp term of the result
        def ctp_term(g, ctp):
            gpart = g_ctp_builder(g, ctp)
            ctp_coeff = self[g].coefficient_cycle_type(gpart) * gpart.aut() / ctp.aut()
            return ctp_coeff * p(ctp)

        # Compute the order-n term of the result as a sum of partitions
        n_term = lambda g, n: sum(ctp_term(g, ctp) for ctp in Partitions(n))

        # Build the g component of the result
        component_builder = lambda g: cisr(n_term(g, n) for n in _integers_from(0))

        # Construct the result
        components_dict = {g: component_builder(g) for g in group}
        result = parent._from_dict(components_dict)
        return result

    def derivative(self):
        r"""
        Return the cycle-index derivative of ``self``.

        Differentiation of group cycle index series is defined termwise:

        .. MATH::
            (F')[\gamma] = (F [\gamma])'

        EXAMPLES::

            sage: from sage.combinat.species.group_cycle_index_series import GroupCycleIndexSeriesRing
            sage: S4 = SymmetricGroup(4)
            sage: GCISR = GroupCycleIndexSeriesRing(S4)
            sage: E = GCISR(species.SetSpecies().cycle_index_series())
            sage: E.derivative()[S4.an_element()].coefficients(6) == E[S4.an_element()].coefficients(6)
            True

        """
        return self.map_coefficients(lambda x: x.derivative())

    def restricted(self, min=None, max=None):
        r"""
        Return the restriction of ``self`` with coefficients starting at degree
        ``min`` and going up to, but not including, degree ``max``.
        If ``min`` is not specified, it is assumed to be zero. If ``max`` is not
        specified, it assumed to be infinity.

        This method simply calls :meth:`~sage.combinat.species.series.LazyPowerSeries.restricted` on each term
        of ``self``.

        EXAMPLES::

            sage: from sage.combinat.species.group_cycle_index_series_library import LinearOrderWithReversalGroupCycleIndex
            sage: L = LinearOrderWithReversalGroupCycleIndex()
            sage: e,t = L.parent().basis().keys()
            sage: L.restricted(min=2,max=4)[t].coefficients(5)
            [0, 0, p[2], p[2, 1], 0]

        """
        return self.map_coefficients(lambda x: x.restricted(min=min, max=max))

    def define(self, x):
        r"""
        Set ``self`` equal to ``x``, possibly recursively.

        EXAMPLES:

        Consider the `\mathfrak{S}_{2}`-species `\mathcal{A}` of rooted ordered trees,
        where the action of the nontrivial element of `\mathfrak{S}_{2}` inverts all the
        orderings on subtrees.

        Then we have the recursive functional equation

        .. MATH::
            \mathcal{A} = X \cdot \mathcal{L} (\mathcal{A})

        where `X` is the trivial singleton `\mathfrak{S}_{2}`-species and `\mathcal{L}` is the
        `\mathfrak{S}_{2}`-species of linear orderings with the reversal action
        described previously.

        To enumerate this species using Sage, we first need to set up an environment::

            sage: from sage.combinat.species.group_cycle_index_series import GroupCycleIndexSeriesRing

            sage: S2 = SymmetricGroup(2)
            sage: GCISR = GroupCycleIndexSeriesRing(S2)
            sage: e,t = GCISR.basis().keys()
            sage: from sage.combinat.species.group_cycle_index_series_library import LinearOrderWithReversalGroupCycleIndex
            sage: L = LinearOrderWithReversalGroupCycleIndex()
            sage: X = GCISR(species.SingletonSpecies().cycle_index_series())

        We can then use ``define`` to set up `\mathcal{A}`::

            sage: A = GCISR(e)*0 + GCISR(t)*0
            sage: A.define(X*L(A))
            sage: A.quotient().isotype_generating_series().counts(8)
            [0, 1, 1, 2, 4, 10, 26, 76]

        (Compare :oeis:`A007123`.)

        """

        keys = self.parent().basis().keys()
        for key in keys:
            self[key].define(x[key])
