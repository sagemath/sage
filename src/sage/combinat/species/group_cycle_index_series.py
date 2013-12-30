"""
Group Cycle Indices

This file implements the group cycle indices of Henderson and Gainer-Dewar.

For a group `\Gamma` and a ring `R`, a `\Gamma`-cycle index over `R`
is an element of the free module over the ring of cycle index series
over `R` with basis `\Gamma`.

In other words, a `\Gamma`-cycle index series over `R` is a formal sum

.. MATH::

    F = \\sum_{\gamma \\in \Gamma} F[\gamma] \cdot \gamma,

where each coefficient `F[\gamma]` is a cycle index series over `R`.

These objects are of interest because they can be used to enumerate `\Gamma`-species;
they serve the same role in that theory as ordinary cycle indices do for classical
species.

AUTHORS:

- Andrew Gainer-Dewar (2013): initial version

EXAMPLES::

    sage: GCISR = gci.GroupCycleIndexSeriesRing(SymmetricGroup(4))
    sage: loads(dumps(GCISR))
    Ring of (Symmetric group of order 4! as a permutation group)-Cycle Index Series over Rational Field

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
from sage.combinat.free_module import CombinatorialFreeModule,CombinatorialFreeModuleElement

@cached_function
def GroupCycleIndexSeriesRing(G, R = RationalField()):
    """
    Return the ring of group cycle index series.

    EXAMPLES::

        sage: GCISR = gci.GroupCycleIndexSeriesRing(SymmetricGroup(4)); GCISR
        Ring of (Symmetric group of order 4! as a permutation group)-Cycle Index Series over Rational Field

    TESTS: We test to make sure that caching works. ::

        sage: GCISR is gci.GroupCycleIndexSeriesRing(SymmetricGroup(4))
        True
    """
    return GroupCycleIndexSeriesRing_class(G, R)

class GroupCycleIndexSeriesRing_class(CombinatorialFreeModule):
    def __init__(self, G, R = RationalField()):
        """
        EXAMPLES::

            sage: GCISR = gci.GroupCycleIndexSeriesRing(SymmetricGroup(4)); GCISR
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
        """
        Return the product of two basis elements ``left`` and ``right`` of ``self``.
        
        Multiplication of `\Gamma`-cycle indices is defined componentwise.
        That is, if `F` and `G` are two `\Gamma`-cycle indices, then `(F \cdot G) [\gamma] = F [\gamma] \cdot G [\gamma]`,
        where the multiplication on the right-hand side is ordinary multiplication of cycle indices.

        This is handled in Sage by defining multiplication on the basis of monomials induced
        by elements of `\Gamma`.

        EXAMPLES::

            sage: GCISR = gci.GroupCycleIndexSeriesRing(SymmetricGroup(4))
            sage: e = SymmetricGroup(4).identity()
            sage: t = SymmetricGroup(4)([4,3,2,1])
            sage: GCISR.product_on_basis(t,t) == GCISR(t)
            True
        
        """
        if left == right:
            return self.monomial(left)
        else:
            return self(0)

    def one(self):
        """
        Return the multiplicative identity element of this algebra.

        EXAMPLES::

            sage: GCISR = gci.GroupCycleIndexSeriesRing(SymmetricGroup(2))
            sage: GCISR.one()
            p[]*G[()] + p[]*G[(1,2)]

        """
        basis = self.basis()
        return self.sum(basis[k] for k in basis.keys())

    def _repr_(self):
        """
        EXAMPLES::

            sage: gci.GroupCycleIndexSeriesRing(SymmetricGroup(4))
            Ring of (Symmetric group of order 4! as a permutation group)-Cycle Index Series over Rational Field
        """
        return "Ring of (%s)-Cycle Index Series over %s" %(self._group, self._coeff_ring)

class GroupCycleIndexSeries(CombinatorialFreeModuleElement):
    """
    
    EXAMPLES::

        sage: GCISR = gci.GroupCycleIndexSeriesRing(SymmetricGroup(4))
        sage: GCISR.an_element()
        p[]*G[()] + 2*p[]*G[(3,4)] + 3*p[]*G[(2,3)] + p[]*G[(1,2,3,4)]
    
    """
    
    def quotient(self):
        """
        Return the quotient of this group cycle index.

        This is defined to be the ordinary cycle index `F / \Gamma` obtained from a
        `\Gamma`-cycle index `F` by:

        .. MATH::
            F / \Gamma = \\frac{1}{\\lvert \Gamma \\rvert} \sum_{\gamma \in \Gamma} F [\gamma].

        By [AGdiss]_, if `F` is the `\Gamma`-cycle index of a `\Gamma`-species, `F / \Gamma` is the ordinary
        cycle index of orbits of structures under the action of `\Gamma`.

        EXAMPLES::

            sage: S4 = SymmetricGroup(4)
            sage: GCISR = gci.GroupCycleIndexSeriesRing(S4)
            sage: GCISR.an_element()
            p[]*G[()] + 2*p[]*G[(3,4)] + 3*p[]*G[(2,3)] + p[]*G[(1,2,3,4)]
            sage: GCISR.an_element().quotient().coefficients(4)
            [7/24*p[], 0, 0, 0]

        REFERENCES:

        .. [AGdiss] Andrew Gainer. "`\Gamma`-species, quotients, and graph enumeration". Ph.D. diss. Brandeis University, 2012.
        """
        return 1/self.parent()._group.cardinality() * sum(self.coefficients())
    
    def functorial_composition(self, y):
        """
        Return the functorial composition of ``self`` with ``y``.

        Functorial composition of group cycle index series is defined componentwise:

        .. MATH::
          (F \square G) [\gamma] = (F [\gamma]) \square (G [\gamma])

        This operation on $\Gamma$-cycle indices corresponds to the functorial composition
        operation on $\Gamma$-species.

        EXAMPLES::

            sage: S4 = SymmetricGroup(4)
            sage: GCISR = gci.GroupCycleIndexSeriesRing(S4)
            sage: E = species.SetSpecies().cycle_index_series()
            sage: Eg = GCISR(E)
            sage: Egplus = GCISR(E-1)
            sage: test = Eg.functorial_composition(Egplus)
            sage: test[SymmetricGroup(4).an_element()].coefficients(4)
            [p[], p[1], 1/2*p[1, 1] + 1/2*p[2], 1/6*p[1, 1, 1] + 1/2*p[2, 1] + 1/3*p[3]]
        """

        parent = self.parent()
        group = parent._group

        components_dict = { g: self[g].functorial_composition(y[g]) for g in group }
        res = parent._from_dict(components_dict)
        return res

    def composition(self, y, test_constant_term_is_zero = True):
        """
        Return the group-cycle-index plethysm of ``self`` with ``y``.
        
        Plethysm of group cycle index series is defined by a sort of 'mixing' operation in [Hend]_:

        .. MATH::
            (F \circ G) [\gamma] (p_{1}, p_{2}, p_{3}, \dots) =
            F [\gamma] \left( G [\gamma] (p_{1}, p_{2}, p_{3}, \dots),
            G [\gamma^{2}] (p_{2}, p_{4}, p_{6}, \dots), \dots \\right).

        It is shown in [Hend]_ that this operation on $\Gamma$-cycle indices corresponds to the
        'composition' operation on $\Gamma$-species.

        It is required that each of the components of `y` has zero constant term.
        However, testing this may not be possible if `y` is of an exotic class.
        Set ``test_constant_term_is_zero`` to ``False`` to suppress any testing.

        EXAMPLES::

            sage: S2 = SymmetricGroup(2)
            sage: Eplus = sage.combinat.species.set_species.SetSpecies(min=1).cycle_index_series()
            sage: E = sage.combinat.species.set_species.SetSpecies().cycle_index_series()
            sage: GCISR = gci.GroupCycleIndexSeriesRing(S2)
            sage: e = S2.identity()
            sage: t = S2.gen()
            sage: GCIS = GCISR(e)*Eplus*E + GCISR(t)*Eplus
            sage: example = GCIS(GCIS)
            sage: example[e].coefficients(4)
            [0, p[1], 3*p[1, 1] + p[2], 41/6*p[1, 1, 1] + 9/2*p[2, 1] + 2/3*p[3]]
            sage: example[t].coefficients(4)
            [0, p[1], p[1, 1] + p[2], 5/6*p[1, 1, 1] + 3/2*p[2, 1] + 2/3*p[3]]

        REFERENCES:

        .. [Hend] Anthony Henderson. "Species over a finite field". J. Algebraic Combin., 21(2):147-161, 2005.
        """
        from sage.combinat.species.stream import Stream,_integers_from
        from sage.misc.misc_c import prod

        assert self.parent() == y.parent()

        parent = self.parent()
        cisr = parent._cisr
        sfa = cisr.base_ring()
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
            else:
                res = sum(coeff*monomial_composer(part, g) for part,coeff in term)
                return res

        def component_builder( g ):
            if test_constant_term_is_zero and self[g] == 0:
                res = cisr.zero()
            #elif g == group.identity(): #Save time by using existing CIS composition code
            #    res = self[g](y[g])
            else:
                res = cisr.sum_generator(term_map(self[g].coefficient(i), g) for i in _integers_from(0))
            return res

        components_dict = { g: component_builder(g) for g in group }
        res = parent._from_dict(components_dict)
        return res

    __call__ = composition

    def derivative(self):
        """
        Return the cycle-index derivative of ``self``.
        
        Differentiation of group cycle index series is defined termwise:

        .. MATH::
            (F')[\gamma] = (F [\gamma])'

        EXAMPLES::

            sage: S4 = SymmetricGroup(4)
            sage: GCISR = gci.GroupCycleIndexSeriesRing(S4)
            sage: E = GCISR(species.SetSpecies().cycle_index_series())
            sage: E.derivative()[S4.an_element()].coefficients(6) == E[S4.an_element()].coefficients(6)
            True

        """
        return self.map_coefficients(lambda x: x.derivative())

    def restricted(self, min=None, max=None):
        """
        Return the restriction of ``self`` with coefficients starting at degree
        ``min`` and going up to, but not including, degree ``max``.
        If ``min`` is not specified, it is assumed to be zero. If ``max`` is not
        specified, it assumed to be infinity.

        This method simply calls `sage.combinat.species.series.restricted` on each term
        of ``self``.

        EXAMPLES::

            sage: L = gci.LinearOrderWithReversalGroupCycleIndex()
            sage: e,t = L.parent().basis().keys()
            sage: L.restricted(min=2,max=4)[t].coefficients(5)
            [0, 0, p[2], p[2, 1], 0]

        """
        return self.map_coefficients(lambda x: x.restricted(min=min, max=max))

    def define(self, x):
        """
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
        
            sage: from sage.combinat.species.generating_series import CycleIndexSeriesRing
            sage: from sage.combinat.species.stream import _integers_from
            sage: S2 = SymmetricGroup(2)
            sage: GCISR = gci.GroupCycleIndexSeriesRing(S2)
            sage: CISR = CycleIndexSeriesRing(QQ)
            sage: e,t = GCISR.basis().keys()

        Next, we define the `\mathfrak{S}_{2}`-species `X` and `\mathcal{L}`.
        `X` is straightforward, since the `\mathfrak{S}_{2}`-action is trivial,
        but `\mathcal{L}` requires some work. ::
        
            sage: X = GCISR(species.SingletonSpecies().cycle_index_series())
            sage: def Ltgen():
            ....:     p = SymmetricFunctions(QQ).power()
            ....:     for n in _integers_from(0):
            ....:         yield p([2]*n)
            ....:         yield p([2]*n+[1])
            sage: L = GCISR(e)*species.LinearOrderSpecies().cycle_index_series() + GCISR(t)*CISR(Ltgen())

        Finally, we can use ``define`` to set up `\mathcal{A}`::
        
            sage: A = GCISR(e)*0 + GCISR(t)*0
            sage: A.define(X*L(A))
            sage: A.quotient().isotype_generating_series().counts(8)
            [0, 1, 1, 2, 4, 10, 26, 76]

        (Compare :oeis:`A007123`.)
        
        """        
        
        keys = self.parent().basis().keys()
        for key in keys:
            self[key].define(x[key])
