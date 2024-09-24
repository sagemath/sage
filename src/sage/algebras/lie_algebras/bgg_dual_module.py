r"""
BGG Category O Dual Modules

AUTHORS:

- Travis Scrimshaw (2024-01-07): Initial version
"""

#*****************************************************************************
#       Copyright (C) 2024 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.cachefunc import cached_method
from sage.categories.enumerated_sets import EnumeratedSets
from sage.categories.monoids import Monoids
from sage.structure.parent import Parent
from sage.structure.indexed_generators import IndexedGenerators
from sage.monoids.indexed_free_monoid import IndexedFreeAbelianMonoid, IndexedMonoid
from sage.combinat.free_module import CombinatorialFreeModule
from sage.sets.family import Family
from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ
from sage.data_structures.blas_dict import iaxpy
from sage.algebras.lie_algebras.verma_module import ModulePrinting


class BGGDualModule(CombinatorialFreeModule):
    r"""
    The dual module `M^{\vee}` in the BGG Category `\mathcal{O}`.

    Let `\tau` be the transpose map of a semisimple (finite dimensional)
    Lie algebra `\mathfrak{g}` over a field `R`. Let `M \in \mathcal{O}`.
    The *BGG dual module* is the `R`-module `M^{\vee} :=
    \bigoplus_{\lambda} M_{\lambda}^*` which has a `U(\mathfrak{g})`-module
    structure given by

    .. MATH::

        x \cdot \phi(v) := \phi(\tau(x) \cdot v),

    which is also a weight module with the same grading as `M`.

    The basis we chose to work with here is the natural dual basis to the
    distinguished basis `B` of `M`. That is, we define the dual function
    to `b` as `\phi_b(c) = \delta_{bc}`.

    EXAMPLES::

        sage: g = LieAlgebra(QQ, cartan_type=['A', 1])
        sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
        sage: M = g.verma_module(2*La[1])
        sage: Mc = M.dual()
        sage: B = Mc.basis()
        sage: it = iter(B)
        sage: elts = [next(it) for _ in range(7)]; elts
        [v[2*Lambda[1]]^*,
         f[-alpha[1]]*v[2*Lambda[1]]^*,
         f[-alpha[1]]^2*v[2*Lambda[1]]^*,
         f[-alpha[1]]^3*v[2*Lambda[1]]^*,
         f[-alpha[1]]^4*v[2*Lambda[1]]^*,
         f[-alpha[1]]^5*v[2*Lambda[1]]^*,
         f[-alpha[1]]^6*v[2*Lambda[1]]^*]
        sage: e, h, f = g.pbw_basis().algebra_generators()
        sage: [f * vec for vec in elts]
        [2*f[-alpha[1]]*v[2*Lambda[1]]^*,
         2*f[-alpha[1]]^2*v[2*Lambda[1]]^*,
         0,
         -4*f[-alpha[1]]^4*v[2*Lambda[1]]^*,
         -10*f[-alpha[1]]^5*v[2*Lambda[1]]^*,
         -18*f[-alpha[1]]^6*v[2*Lambda[1]]^*,
         -28*f[-alpha[1]]^7*v[2*Lambda[1]]^*]
        sage: [e * vec for vec in elts]
        [0,
         v[2*Lambda[1]]^*,
         f[-alpha[1]]*v[2*Lambda[1]]^*,
         f[-alpha[1]]^2*v[2*Lambda[1]]^*,
         f[-alpha[1]]^3*v[2*Lambda[1]]^*,
         f[-alpha[1]]^4*v[2*Lambda[1]]^*,
         f[-alpha[1]]^5*v[2*Lambda[1]]^*]
        sage: [h * vec for vec in elts]
        [2*v[2*Lambda[1]]^*,
         0,
         -2*f[-alpha[1]]^2*v[2*Lambda[1]]^*,
         -4*f[-alpha[1]]^3*v[2*Lambda[1]]^*,
         -6*f[-alpha[1]]^4*v[2*Lambda[1]]^*,
         -8*f[-alpha[1]]^5*v[2*Lambda[1]]^*,
         -10*f[-alpha[1]]^6*v[2*Lambda[1]]^*]

    REFERENCES:

    - [Humphreys08]_
    """
    def __init__(self, module):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['B', 2])
            sage: La = g.cartan_type().root_system().weight_space().fundamental_weights()
            sage: M = g.verma_module(2*La[1] + La[2])
            sage: Mc = M.dual()
            sage: TestSuite(Mc).run()

            sage: M = g.verma_module(2/3*La[1] - 3/5*La[2])
            sage: Mc = M.dual()
            sage: TestSuite(Mc).run()
        """
        self._module = module
        self._g = module.lie_algebra()
        self._pbw = self._g.pbw_basis()
        base_ring = module.base_ring()
        indices = module.indices()
        category = module.category()
        CombinatorialFreeModule.__init__(self, base_ring, indices, category=category,
                                         **module.print_options())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['A', 1])
            sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = g.verma_module(2*La[1])
            sage: M.dual()
            BGG Dual of Verma module with highest weight 2*Lambda[1] of
             Lie algebra of ['A', 1] in the Chevalley basis
        """
        return "BGG Dual of " + repr(self._module)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['A', 1])
            sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = g.verma_module(2*La[1])
            sage: Mc = M.dual()
            sage: latex(Mc)
            { M_{2 \Lambda_{1}} }^{\vee}
        """
        from sage.misc.latex import latex
        return "{" + latex(self._module) + "}^{\\vee}"

    def _repr_generator(self, m):
        r"""
        Return a string representation of the generator indexed by ``m``.

        EXAMPLES::

            sage: g = lie_algebras.sp(QQ, 4)
            sage: La = g.cartan_type().root_system().ambient_space().fundamental_weights()
            sage: Mc = g.verma_module(La[1] + 3/7*La[2]).dual()
            sage: f1, f2 = g.f()
            sage: x = g.pbw_basis()(g([f1, [f1, f2]]))
            sage: v = x * Mc.highest_weight_vector()
            sage: Mc._repr_generator(v.leading_support())
            'f[-alpha[1]]*f[-alpha[1] - alpha[2]]*v[(10/7, 3/7)]^*'
        """
        return self._module._repr_generator(m) + "^*"

    def _latex_generator(self, m):
        """
        Return a latex representation of the generator indexed by ``m``.

        EXAMPLES::

            sage: g = lie_algebras.sp(QQ, 4)
            sage: La = g.cartan_type().root_system().ambient_space().fundamental_weights()
            sage: Mc = g.verma_module(La[1] + 3/7*La[2]).dual()
            sage: f1, f2 = g.f()
            sage: x = g.pbw_basis()(g([f1, [f1, f2]]))
            sage: v = x * Mc.highest_weight_vector()
            sage: Mc._latex_generator(v.leading_support())
            { f_{-\alpha_{1}} f_{-\alpha_{1} - \alpha_{2}} v_{\frac{10}{7} e_{0} + \frac{3}{7} e_{1}} }^{\vee}
        """
        return "{" + self._module._latex_generator(m) + "}^{\\vee}"

    _repr_term = _repr_generator
    _latex_term = _latex_generator

    def degree_on_basis(self, m):
        r"""
        Return the degree of the basis element indexed by ``m``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['D', 5])
            sage: La = g.cartan_type().root_system().weight_space().fundamental_weights()
            sage: M = g.verma_module(La[1] + La[4] - 1/3*La[5])
            sage: Mc = M.dual()
            sage: elt = Mc.an_element(); elt
            f[-alpha[5]]^2*f[-alpha[4]]^2*f[-alpha[3]]^3*v[Lambda[1] + Lambda[4] - 1/3*Lambda[5]]^*
             + 2*f[-alpha[5]]*v[Lambda[1] + Lambda[4] - 1/3*Lambda[5]]^*
             + 3*f[-alpha[4]]*v[Lambda[1] + Lambda[4] - 1/3*Lambda[5]]^*
             + v[Lambda[1] + Lambda[4] - 1/3*Lambda[5]]^*
            sage: [M.degree_on_basis(m) for m in elt.support()]
            [Lambda[1] + 3*Lambda[2] - 2*Lambda[3] - 4/3*Lambda[5],
             Lambda[1] + Lambda[4] - 1/3*Lambda[5],
             Lambda[1] + Lambda[3] + Lambda[4] - 7/3*Lambda[5],
             Lambda[1] + Lambda[3] - Lambda[4] - 1/3*Lambda[5]]
        """
        return self._module.degree_on_basis(m)

    def highest_weight(self):
        r"""
        Return the highest weight of ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['E', 7])
            sage: La = g.cartan_type().root_system().weight_space().fundamental_weights()
            sage: M = g.verma_module(2*La[1] + 5/3*La[4] - 3*La[6])
            sage: Mc = M.dual()
            sage: Mc.highest_weight()
            2*Lambda[1] + 5/3*Lambda[4] - 3*Lambda[6]
        """
        return self._module.highest_weight()

    def highest_weight_vector(self):
        r"""
        Return the highest weight vector of ``self`` (assuming the
        defining module defines such a vector).

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['A', 1])
            sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = g.verma_module(2*La[1])
            sage: Mc = M.dual()
            sage: Mc.highest_weight_vector()
            v[2*Lambda[1]]^*
        """
        hwv = self._module.highest_weight_vector()
        return self.element_class(self, hwv.monomial_coefficients(copy=False))

    def lie_algebra(self):
        r"""
        Return the underlying Lie algebra of ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['B', 3])
            sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = g.verma_module(2*La[1] + La[3])
            sage: Mc = M.dual()
            sage: Mc.lie_algebra() is g
            True
        """
        return self._g

    def dual(self):
        r"""
        Return the dual module of ``self``.

        In Category `\mathcal{O}`, we have `(M^{\vee})^{\vee} \cong M`, so
        we return the defining module `M` of `M^{\vee}`.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['F', 4])
            sage: La = g.cartan_type().root_system().weight_space().fundamental_weights()
            sage: M = g.verma_module(La[1] - 5/3*La[2] + 3*La[4])
            sage: Mc = M.dual()
            sage: Mc.dual() is M
            True
        """
        return self._module

    @cached_method
    def _lie_algebra_on_basis(self, b, m):
        r"""
        Return the action of the Lie algebra basis element indexed by ``b``
        on the basis element of ``self`` indexed by ``m``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['B', 2])
            sage: La = g.cartan_type().root_system().weight_space().fundamental_weights()
            sage: M = g.verma_module(La[1])
            sage: Mc = M.dual()
            sage: it = iter(Mc.basis())
            sage: list(g.basis())
            [E[alpha[2]], E[alpha[1]], E[alpha[1] + alpha[2]], E[alpha[1] + 2*alpha[2]],
             h1, h2,
             E[-alpha[2]], E[-alpha[1]], E[-alpha[1] - alpha[2]], E[-alpha[1] - 2*alpha[2]]]
            sage: for _ in range(3):
            ....:     m = next(it).leading_support()
            ....:     print(m, [Mc._lie_algebra_on_basis(k, m) for k in g.basis().keys()])
            1 [0, 0, 0, 0, v[Lambda[1]]^*, 0, 0, f[-alpha[1]]*v[Lambda[1]]^*,
               -2*f[-alpha[2]]*f[-alpha[1]]*v[Lambda[1]]^* + 2*f[-alpha[1] - alpha[2]]*v[Lambda[1]]^*,
               -2*f[-alpha[2]]^2*f[-alpha[1]]*v[Lambda[1]]^*
                + 2*f[-alpha[2]]*f[-alpha[1] - alpha[2]]*v[Lambda[1]]^*
                + f[-alpha[1] - 2*alpha[2]]*v[Lambda[1]]^*]
            f[-alpha[2]] [v[Lambda[1]]^*, 0, 0, 0, 2*f[-alpha[2]]*v[Lambda[1]]^*,
                -2*f[-alpha[2]]*v[Lambda[1]]^*, -2*f[-alpha[2]]^2*v[Lambda[1]]^*,
                f[-alpha[2]]*f[-alpha[1]]*v[Lambda[1]]^* + f[-alpha[1] - alpha[2]]*v[Lambda[1]]^*,
                -4*f[-alpha[2]]^2*f[-alpha[1]]*v[Lambda[1]]^* - f[-alpha[1] - 2*alpha[2]]*v[Lambda[1]]^*,
                -6*f[-alpha[2]]^3*f[-alpha[1]]*v[Lambda[1]]^* + 2*f[-alpha[2]]^2*f[-alpha[1] - alpha[2]]*v[Lambda[1]]^*]
            f[-alpha[1]] [0, v[Lambda[1]]^*, 0, 0, -f[-alpha[1]]*v[Lambda[1]]^*,
                2*f[-alpha[1]]*v[Lambda[1]]^*,
                2*f[-alpha[2]]*f[-alpha[1]]*v[Lambda[1]]^* - 2*f[-alpha[1] - alpha[2]]*v[Lambda[1]]^*,
                0, 0, f[-alpha[1]]*f[-alpha[1] - 2*alpha[2]]*v[Lambda[1]]^* + 2*f[-alpha[1] - alpha[2]]^2*v[Lambda[1]]^*]
        """
        al = self._g.degree_on_basis(b)
        wt = self.degree_on_basis(m)
        if al == 0:  # b is indexing part of the Cartan subalgebra
            # We are assuming b is part of the coroot lattice.
            # FIXME: Add something at the category level to return this.
            ac = b
            return self.term(m, wt.scalar(ac))

        # TODO: Avoid calling homogeneous_component_basis() as the result is not cached
        gens = self._module.homogeneous_component_basis(wt + al)
        elt = self._g.basis()[b]
        # TODO: Determine if we can meaningfully store these results.
        # Computing gens is ~1/3 of the computation and vecs is ~2/3.
        vecs = {g.leading_support(): elt.transpose() * g for g in gens}
        return self.element_class(self, {k: c for k, v in vecs.items() if (c := v[m])})

    def _pbw_monomial_on_basis(self, p, m):
        r"""
        Return the action of the PBW monomial indexed by ``p`` on the basis
        element of ``self`` indexed by ``m``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['A', 1])
            sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: PBW = g.pbw_basis()
            sage: e, h, f = PBW.algebra_generators()
            sage: M = g.verma_module(2*La[1])
            sage: Mc = M.dual()
            sage: v = Mc.highest_weight_vector()
            sage: Mc._pbw_monomial_on_basis((e*f^2).leading_support(), v.leading_support())
            4*f[-alpha[1]]*v[2*Lambda[1]]^*
            sage: B = Mc.basis()
            sage: it = iter(B)
            sage: elts = [next(it) for _ in range(7)]; elts
            [v[2*Lambda[1]]^*,
             f[-alpha[1]]*v[2*Lambda[1]]^*,
             f[-alpha[1]]^2*v[2*Lambda[1]]^*,
             f[-alpha[1]]^3*v[2*Lambda[1]]^*,
             f[-alpha[1]]^4*v[2*Lambda[1]]^*,
             f[-alpha[1]]^5*v[2*Lambda[1]]^*,
             f[-alpha[1]]^6*v[2*Lambda[1]]^*]
        """
        ret = self.monomial(m)
        for b, exp in reversed(p._sorted_items()):
            for _ in range(exp):
                ret = self.linear_combination((self._lie_algebra_on_basis(b, m), mc)
                                              for m, mc in ret._monomial_coefficients.items())
        return ret

    class Element(CombinatorialFreeModule.Element):
        def _acted_upon_(self, scalar, self_on_left=False):
            r"""
            Return the action of ``scalar`` on ``self``.

            EXAMPLES::

                sage: g = LieAlgebra(QQ, cartan_type=['A', 1])
                sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
                sage: PBW = g.pbw_basis()
                sage: e, h, f = PBW.algebra_generators()
                sage: M = g.verma_module(2*La[1])
                sage: Mc = M.dual()
                sage: v = Mc.highest_weight_vector()
                sage: (h*e^2*f^2) * v
                8*v[2*Lambda[1]]^*
                sage: g.casimir_element(UEA=PBW) * v
                v[2*Lambda[1]]^*
                sage: 5 * v
                5*v[2*Lambda[1]]^*
            """
            P = self.parent()
            # Check for scalars first
            if scalar in P.base_ring():
                # Don't have this be a super call
                return CombinatorialFreeModule.Element._acted_upon_(self, scalar, self_on_left)

            # Check for Lie algebra elements
            try:
                scalar = P._g(scalar)
            except (ValueError, TypeError):
                pass
            if scalar.parent() is P._g:
                if self_on_left:  # only implemented as a left module
                    return None
                mc = scalar.monomial_coefficients(copy=False)
                return P.linear_combination((P._lie_algebra_on_basis(b, m), bc * mc)
                                            for b, bc in mc.items()
                                            for m, mc in self._monomial_coefficients.items())

            # Check for PBW elements
            try:
                scalar = P._pbw(scalar)
            except (ValueError, TypeError):
                # Cannot be made into a PBW element, so propagate it up
                return CombinatorialFreeModule.Element._acted_upon_(self,
                        scalar, self_on_left)

            # We only implement x * self, i.e., as a left module
            if self_on_left:
                return None

            mc = scalar.monomial_coefficients(copy=False)
            return P.linear_combination((P._pbw_monomial_on_basis(p, m), pc * mc)
                                        for p, pc in mc.items()
                                        for m, mc in self._monomial_coefficients.items())


#####################################################################
# Simple modules


# This is an abuse as the monoid is not free.
# TODO: Rewrite this (or the indexed monoid class) to use explicit vectors
#    since we only want to consider ordered elements.
# Note, such a rewrite would force the Lie algebra to be finite dimensional.
class SimpleModuleIndices(IndexedFreeAbelianMonoid):
    r"""
    The indices of the basis for a simple `U(\mathfrak{g})`-module.

    .. NOTE::

        The current implementation assumes the Lie algebra `\mathfrak{g}`
        is finite dimensional.
    """
    # This is only necessary because of the IndexedMonoid.__classcall__.
    @staticmethod
    def __classcall__(cls, simple, prefix='f', **kwds):
        r"""
        Normalize input to ensure a unique representation.

        TESTS::

            sage: g = LieAlgebra(QQ, cartan_type=['E', 6])
            sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: L = g.simple_module(La[1] + La[3])
            sage: from sage.algebras.lie_algebras.bgg_dual_module import SimpleModuleIndices
            sage: SimpleModuleIndices(L) is L._indices
            True
        """
        return super(IndexedMonoid, cls).__classcall__(cls, simple, prefix=prefix, **kwds)

    def __init__(self, simple, prefix, category=None, **kwds):
        r"""
        Initialize ``self``.

        TESTS::

            sage: g = LieAlgebra(QQ, cartan_type=['A', 2])
            sage: La = g.cartan_type().root_system().weight_space().fundamental_weights()
            sage: I = g.simple_module(2*La[1] + La[2]).indices()
            sage: TestSuite(I).run()

            sage: I = g.simple_module(2*La[1] - 1/3*La[2]).indices()
            sage: TestSuite(I).run(max_runs=150)  # long time
        """
        self._simple = simple
        self._g = simple.lie_algebra()
        self._reached_max_depth = False
        # Below was mostly copied from IndexedMonoid.__init__()
        self._indices = FiniteEnumeratedSet(self._g._negative_half_index_set())
        category = Monoids().or_subcategory(category)
        category = category & EnumeratedSets()
        category = category.FinitelyGeneratedAsMagma()
        if self._simple._dom_int:
            category = category.Finite()
        else:
            category = category.Infinite()
        Parent.__init__(self, category=category)

        # ignore the optional 'key' since it only affects CachedRepresentation
        kwds.pop('key', None)
        sorting_key = kwds.pop('sorting_key', self._simple._pbw._monoid_key)
        IndexedGenerators.__init__(self, self._indices, prefix, sorting_key=sorting_key, **kwds)

        self._sorted_supp = sorted(self._g._negative_half_index_set(), key=self._simple._pbw._basis_key,
                                   reverse=self.print_options()['sorting_reverse'])
        self._basis = {self.one(): self._simple._ambient.highest_weight_vector()}
        self._lead_supp_to_index = {self._simple._ambient.highest_weight_vector().leading_support(): self.one()}
        # This is used for iteration and keeps track of the current depth
        self._basis_by_depth = [dict(self._basis)]
        # The basis is given as a list of indices corresponding to basis vectors in self._basis
        self._weight_space_bases = {self._simple.highest_weight(): [self.one()]}

    def _an_element_(self):
        r"""
        Return an element of ``self``.

        The only element we can quickly guarantee is in ``self`` is 1,
        so we return this.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['E', 6])
            sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: I = g.simple_module(2*La[1] + La[2]).indices()
            sage: I._an_element_()
            1
        """
        return self.one()

    def _weight_max_depth(self, mu):
        r"""
        Return the maximum depth of the weight ``mu``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['E', 6])
            sage: P = g.cartan_type().root_system().weight_lattice()
            sage: La = P.fundamental_weights()
            sage: al = P.simple_roots()
            sage: wt = 2*La[1] + La[2]
            sage: I = g.simple_module(wt).indices()
            sage: I._weight_max_depth(wt)
            0
            sage: I._weight_max_depth(wt + al[2]) is None
            True
            sage: I._weight_max_depth(wt - 2*al[2] - 5*al[4] - 3*al[6])
            10

            sage: g = LieAlgebra(QQ, cartan_type=['F', 4])
            sage: P = g.cartan_type().root_system().weight_space()
            sage: La = P.fundamental_weights()
            sage: al = P.simple_roots()
            sage: wt = 2*La[1] - 3/2*La[2]
            sage: I = g.simple_module(wt).indices()
            sage: I._weight_max_depth(wt)
            0
            sage: I._weight_max_depth(wt + al[2]) is None
            True
            sage: I._weight_max_depth(wt - 2*al[2] - 3*al[4])
            5
            sage: I._weight_max_depth(wt - 2/3*al[1]) is None
            True
        """
        al = (self._simple.highest_weight() - mu)._to_root_vector()
        if any(c not in ZZ or c < 0 for c in al):
            return None
        return sum(al)

    def weight_space_basis(self, mu):
        r"""
        Return the indices of the ``mu`` weight space basis elements.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['A', 2])
            sage: P = g.cartan_type().root_system().weight_lattice()
            sage: La = P.fundamental_weights()
            sage: al = P.simple_roots()
            sage: wt = -3*La[1] + 3*La[2]
            sage: I = g.simple_module(wt).indices()
            sage: I.weight_space_basis(wt)
            [1]
            sage: I.weight_space_basis(wt - al[1])
            [f[-alpha[1]]]
            sage: I.weight_space_basis(wt - al[2])
            [f[-alpha[2]]]
            sage: I.weight_space_basis(wt - al[1] - al[2])
            [f[-alpha[1] - alpha[2]], f[-alpha[2]]*f[-alpha[1]]]
            sage: I.weight_space_basis(wt - 4*al[1])
            [f[-alpha[1]]^4]
            sage: I.weight_space_basis(wt - 4*al[2])
            []
        """
        if self._reached_max_depth:
            return self._weight_space_bases.get(mu, [])

        max_depth = self._weight_max_depth(mu)
        while max_depth >= len(self._basis_by_depth):
            if self._reached_max_depth:  # we've already reached everything
                break
            self._construct_next_level()
        return self._weight_space_bases.get(mu, [])

    def __contains__(self, m):
        r"""
        Check if ``m`` is contained in ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['G', 2])
            sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: L = g.simple_module(La[1])
            sage: I = L.indices()
            sage: I.one() in I
            True
            sage: it = iter(I)
            sage: for _ in range(3):
            ....:     elt = next(it)
            ....:     print(elt, elt in I)
            1 True
            f[-alpha[1]] True
            f[-alpha[1] - alpha[2]] True
            sage: gens = list(I.gens()); gens
            [f[-alpha[2]],
             f[-alpha[1]],
             f[-alpha[1] - alpha[2]],
             f[-2*alpha[1] - alpha[2]],
             f[-3*alpha[1] - alpha[2]],
             f[-3*alpha[1] - 2*alpha[2]]]
            sage: gens[1] in I
            True
            sage: gens[0] * gens[1] in I
            False
            sage: gens[2] in I
            True
            sage: gens[0]^10 in I
            False
            sage: gens[5]^6 * gens[2]^10 in I
            False
        """
        if not isinstance(m, self.Element) or m.parent() is not self:
            return False
        depth = m.length()
        while depth >= len(self._basis_by_depth):
            if self._reached_max_depth:  # we've already reached everything
                break
            self._construct_next_level()
        return m in self._basis

    def __iter__(self):
        r"""
        Iterate over ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['B', 2])
            sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: L = g.simple_module(La[2])
            sage: I = L.indices()
            sage: list(I)
            [1, f[-alpha[2]], f[-alpha[1] - alpha[2]], f[-alpha[1] - 2*alpha[2]]]

            sage: L = g.simple_module(La[1]-La[2])
            sage: I = L.indices()
            sage: it = iter(I)
            sage: [next(it) for _ in range(6)]
            [1, f[-alpha[2]], f[-alpha[1]], f[-alpha[1] - alpha[2]],
             f[-alpha[1] - 2*alpha[2]], f[-alpha[2]]^2]
        """
        depth = 0
        while True:
            while depth >= len(self._basis_by_depth):
                if self._reached_max_depth:  # we've already reached everything
                    return
                self._construct_next_level()
            yield from self._basis_by_depth[depth]
            depth += 1

    def _construct_next_level(self):
        r"""
        Construct the image for the next level of ``self``.

        ALGORITHM:

        For each image vector of `f_{\beta_1}^{b_1} \cdots f_{\beta_k}^{b_k}
        v_{\lambda}` at the current depth `b_1 + \cdots b_k`, consider the
        image under multiplication by every generator `f_{\alpha}` for
        all `\alpha \leq \beta_1` in the fixed PBW ordering of the (dual)
        Verma module.

        .. TODO::

            Avoid unnecessary computations by using the corresponding
            (combinatorial) crystal.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['C', 3])
            sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: L = g.simple_module(La[1])
            sage: I = L.indices()
            sage: len(I._basis)
            1
            sage: I._construct_next_level()
            sage: len(I._basis)
            6
            sage: I._reached_max_depth
            False
            sage: I._construct_next_level()  # long time
            sage: len(I._basis)  # long time
            6
            sage: I._reached_max_depth  # long time
            True
            sage: I._construct_next_level()  # long time
        """
        if self._reached_max_depth:  # we've already reached everything
            return  # so nothing more to do

        gens = self._g.basis()
        next_level = {}
        R = self._g.base_ring()
        ambient = self._simple._ambient
        pbw = ambient._pbw
        for m, vec in self._basis_by_depth[-1].items():
            # find the first support index
            ind = len(self._sorted_supp)
            for i, ls in enumerate(self._sorted_supp):
                if ls in m._monomial:
                    ind = i + 1
                    break

            for ls in self._sorted_supp[:ind]:
                mp = dict(m._monomial)  # make a (shallow) copy
                mp[ls] = mp.get(ls, 0) + 1
                key = self.element_class(self, mp)
                new_vec = gens[ls] * vec
                if not new_vec:
                    continue
                # Echelonize the corresponding weight space
                mu = ambient.degree_on_basis(key)
                if mu not in self._weight_space_bases:
                    # the only vector in the weight space
                    self._weight_space_bases[mu] = [key]
                    next_level[key] = new_vec
                    self._basis[key] = next_level[key]
                    lead_supp = next_level[key].trailing_support(key=pbw._monomial_key)
                    self._lead_supp_to_index[lead_supp] = key
                    continue

                supp = set()
                wt_basis = [self._basis[k] for k in self._weight_space_bases[mu]]
                wt_basis.append(new_vec)
                for b in wt_basis:
                    supp.update(b.support())
                supp = sorted(supp, key=pbw._monomial_key)
                mat = matrix(R, [[b[s] for s in supp] for b in wt_basis])
                mat.echelonize()
                for i, k in enumerate(self._weight_space_bases[mu]):
                    data = {supp[ind]: R(c) for ind, c in mat[i].iteritems() if c}
                    self._basis[k] = ambient.element_class(ambient, data)
                i = mat.nrows() - 1
                data = {supp[ind]: R(c) for ind, c in mat[i].iteritems() if c}
                if data:
                    next_level[key] = ambient.element_class(ambient, data)
                    self._basis[key] = next_level[key]
                    lead_supp = next_level[key].trailing_support(key=pbw._monomial_key)
                    self._lead_supp_to_index[lead_supp] = key
                    self._weight_space_bases[mu].append(key)

        if not next_level:
            self._reached_max_depth = True
            return
        self._basis_by_depth.append(next_level)

    @cached_method
    def cardinality(self):
        r"""
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['E', 6])
            sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: L = g.simple_module(La[1]+La[4])
            sage: L._indices.cardinality()
            51975
        """
        if self._simple._dom_int:
            weight = self._simple.highest_weight()
            Phi = self._g.cartan_type().root_system()
            P = Phi.weight_lattice()
            coroots = Phi.root_lattice().simple_coroots()
            la = P._from_dict({i: weight.scalar(ac) for i, ac in coroots.items()})
            from sage.combinat.crystals.monomial_crystals import CrystalOfNakajimaMonomials
            return CrystalOfNakajimaMonomials(la).cardinality()
        from sage.rings.infinity import infinity
        return infinity


class SimpleModule(ModulePrinting, CombinatorialFreeModule):
    r"""
    Return the simple module `L_{\lambda}` as the image of the natural
    morphism `\phi: M_{\lambda} \to M_{\lambda}^{\vee}`.
    """
    @staticmethod
    def __classcall_private__(cls, g, weight, *args, **kwds):
        r"""
        Normalize input to ensure a unique representation and return
        the correct type.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['E', 6])
            sage: La = g.cartan_type().root_system().weight_space().fundamental_weights()
            sage: type(g.simple_module(La[1] + La[2]))
            <class 'sage.algebras.lie_algebras.bgg_dual_module.FiniteDimensionalSimpleModule_with_category'>
            sage: type(g.simple_module(La[1] - La[2]))
            <class 'sage.algebras.lie_algebras.bgg_dual_module.SimpleModule_with_category'>
            sage: type(g.simple_module(La[1] + 3/2*La[2]))
            <class 'sage.algebras.lie_algebras.bgg_dual_module.SimpleModule_with_category'>
        """
        if weight.is_dominant_weight():
            return FiniteDimensionalSimpleModule(g, weight, *args, **kwds)
        return super().__classcall__(cls, g, weight, *args, **kwds)

    def __init__(self, g, weight, prefix='f', basis_key=None, **kwds):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['G', 2])
            sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: L = g.simple_module(La[1] + La[2])
            sage: TestSuite(L).run()

            sage: g = LieAlgebra(QQ, cartan_type=['A', 2])
            sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: L = g.simple_module(La[1] - La[2])
            sage: TestSuite(L).run()
        """
        self._g = g
        self._weight = weight
        self._dom_int = weight.is_dominant_weight()
        self._verma = g.verma_module(weight, basis_key=basis_key)
        self._ambient = self._verma.dual()
        self._pbw = self._verma.pbw_basis()
        base_ring = self._g.base_ring()
        indices = SimpleModuleIndices(self, prefix=prefix, **kwds)
        category = self._ambient.category().Subobjects()
        if self._dom_int:
            category = category.FiniteDimensional()
        ModulePrinting.__init__(self, 'u')
        CombinatorialFreeModule.__init__(self, base_ring, indices, category=category,
                                         **self._ambient.print_options())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['A', 1])
            sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: g.simple_module(2*La[1])
            Simple module with highest weight 2*Lambda[1] of
             Lie algebra of ['A', 1] in the Chevalley basis
        """
        return "Simple module with highest weight {} of {}".format(self._weight, self._g)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['A', 1])
            sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: L = g.simple_module(2*La[1])
            sage: latex(L)
            L_{2 \Lambda_{1}}
        """
        from sage.misc.latex import latex
        return "L_{{{}}}".format(latex(self._weight))

    def ambient(self):
        r"""
        Return the ambient module of ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['G', 2])
            sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: L = g.simple_module(2*La[1])
            sage: L.ambient()
            BGG Dual of Verma module with highest weight 2*Lambda[1] of
             Lie algebra of ['G', 2] in the Chevalley basis
        """
        return self._ambient

    @lazy_attribute
    def lift(self):
        r"""
        Return the lift map of ``self`` to the ambient dual Verma module.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['G', 2])
            sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: L = g.simple_module(La[1])
            sage: [L.lift(b) for b in L.basis()]  # long time
            [v[Lambda[1]]^*,
             f[-alpha[1]]*v[Lambda[1]]^*,
             f[-alpha[2]]*f[-alpha[1]]*v[Lambda[1]]^* - f[-alpha[1] - alpha[2]]*v[Lambda[1]]^*,
             f[-alpha[1]]*f[-alpha[1] - alpha[2]]*v[Lambda[1]]^*
              + f[-2*alpha[1] - alpha[2]]*v[Lambda[1]]^*,
             f[-alpha[1]]^2*f[-alpha[1] - alpha[2]]*v[Lambda[1]]^*
              + f[-alpha[1]]*f[-2*alpha[1] - alpha[2]]*v[Lambda[1]]^*
              + 1/2*f[-3*alpha[1] - alpha[2]]*v[Lambda[1]]^*,
             f[-alpha[2]]*f[-alpha[1]]^2*f[-alpha[1] - alpha[2]]*v[Lambda[1]]^*
              + f[-alpha[2]]*f[-alpha[1]]*f[-2*alpha[1] - alpha[2]]*v[Lambda[1]]^*
              + 1/2*f[-alpha[2]]*f[-3*alpha[1] - alpha[2]]*v[Lambda[1]]^*
              - f[-alpha[1] - alpha[2]]*f[-2*alpha[1] - alpha[2]]*v[Lambda[1]]^*
              + 1/2*f[-3*alpha[1] - 2*alpha[2]]*v[Lambda[1]]^*,
             f[-alpha[1]]*f[-alpha[1] - alpha[2]]*f[-2*alpha[1] - alpha[2]]*v[Lambda[1]]^*
              - 1/2*f[-alpha[1]]*f[-3*alpha[1] - 2*alpha[2]]*v[Lambda[1]]^*
              - 1/2*f[-alpha[1] - alpha[2]]*f[-3*alpha[1] - alpha[2]]*v[Lambda[1]]^*
              + f[-2*alpha[1] - alpha[2]]^2*v[Lambda[1]]^*]
        """
        return self.module_morphism(self._lift_on_basis, codomain=self._ambient, unitriangular="upper")

    def retract(self, x):
        r"""
        Return the retraction of ``x`` in ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['A', 2])
            sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: L = g.simple_module(2*La[1])
            sage: L.retract(L.lift(sum(L.basis())))
            f[-alpha[1]]^2*u[2*Lambda[1]] + f[-alpha[1]]*f[-alpha[1] - alpha[2]]*u[2*Lambda[1]]
             + f[-alpha[1] - alpha[2]]^2*u[2*Lambda[1]] + f[-alpha[1]]*u[2*Lambda[1]]
             + f[-alpha[1] - alpha[2]]*u[2*Lambda[1]] + u[2*Lambda[1]]
            sage: B = list(L.basis())
            sage: L.retract(3/2*L.lift(B[0]) - L.lift(B[2]) - 10/3*L.lift(B[3]))
            -10/3*f[-alpha[1]]^2*u[2*Lambda[1]]
             - f[-alpha[1] - alpha[2]]*u[2*Lambda[1]]
             + 3/2*u[2*Lambda[1]]
        """
        supp = sorted(x.support(), key=self._pbw._monomial_key)
        data = x.monomial_coefficients(copy=True)  # this is destructive to data
        R = self.base_ring()
        ret = {}
        for ls in supp:
            if ls not in data:
                continue
            if ls not in self._indices._lead_supp_to_index:
                mu = self._ambient.degree_on_basis(ls)
                # this will guarantee the computation is correct
                self._indices.weight_space_basis(mu)
            if ls not in self._indices._lead_supp_to_index:
                raise ValueError(f"not an element of the simple module of weight {self._weight}")
            key = self._indices._lead_supp_to_index[ls]
            vec = self._indices._basis[key]
            coeff = R(data[ls] / vec[ls])
            iaxpy(-coeff, vec._monomial_coefficients, data)
            ret[key] = coeff
        return self.element_class(self, ret)

    def _lift_on_basis(self, m):
        r"""
        Return the lift of the basis element indexed by ``m``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['A', 1])
            sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: L = g.simple_module(2*La[1])
            sage: I = L.indices()
            sage: gen = list(I.gens())[0]
            sage: L._lift_on_basis(gen^2)
            4*f[-alpha[1]]^2*v[2*Lambda[1]]^*
            sage: L._lift_on_basis(gen^3)
            Traceback (most recent call last):
            ...
            ValueError: f[-alpha[1]]^3 does not index a basis element
        """
        # This builds the result up to the necessary depth
        if m not in self._indices:
            raise ValueError(f"{m} does not index a basis element")
        return self._indices._basis[m]

    def dual(self):
        r"""
        Return the dual module of ``self``, which is ``self`` since simple
        modules are self-dual.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['B', 4])
            sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: L = g.simple_module(2*La[1] + 3*La[4])
            sage: L.dual() is L
            True
        """
        return self

    def highest_weight(self):
        r"""
        Return the highest weight of ``self``.

        EXAMPLES::

            sage: g = lie_algebras.so(QQ, 7)
            sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: L = g.simple_module(La[1] + La[2])
            sage: L.highest_weight()
            Lambda[1] + Lambda[2]
        """
        return self._weight

    @cached_method
    def highest_weight_vector(self):
        r"""
        Return the highest weight vector of ``self``.

        EXAMPLES::

            sage: g = lie_algebras.sp(QQ, 6)
            sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: L = g.simple_module(La[1] + La[2])
            sage: L.highest_weight_vector()
            u[Lambda[1] + Lambda[2]]
        """
        one = self.base_ring().one()
        return self._from_dict({self._indices.one(): one},
                               remove_zeros=False, coerce=False)

    def lie_algebra(self):
        r"""
        Return the underlying Lie algebra of ``self``.

        EXAMPLES::

            sage: g = lie_algebras.so(QQ, 9)
            sage: La = g.cartan_type().root_system().weight_space().fundamental_weights()
            sage: L = g.simple_module(La[3] - 1/2*La[1])
            sage: L.lie_algebra()
            Lie algebra of ['B', 4] in the Chevalley basis
        """
        return self._g

    def pbw_basis(self):
        r"""
        Return the PBW basis of the underlying Lie algebra
        used to define ``self``.

        EXAMPLES::

            sage: g = lie_algebras.so(QQ, 8)
            sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: L = g.simple_module(La[2] - 2*La[3])
            sage: L.pbw_basis()
            Universal enveloping algebra of Lie algebra of ['D', 4] in the Chevalley basis
             in the Poincare-Birkhoff-Witt basis
        """
        return self._pbw

    def homogeneous_component_basis(self, mu):
        r"""
        Return a basis for the ``mu`` weight space of ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['A', 2])
            sage: P = g.cartan_type().root_system().weight_lattice()
            sage: La = P.fundamental_weights()
            sage: la = La[1] + La[2]
            sage: L = g.simple_module(la)
            sage: from itertools import product
            sage: al = P.simple_roots()
            sage: for wts in product(range(4), repeat=2):
            ....:     mu = la - wts[0] * al[1] - wts[1] * al[2]
            ....:     print(mu)
            ....:     print(L.homogeneous_component_basis(mu))
            Lambda[1] + Lambda[2]
            Family (u[Lambda[1] + Lambda[2]],)
            2*Lambda[1] - Lambda[2]
            Family (f[-alpha[2]]*u[Lambda[1] + Lambda[2]],)
            3*Lambda[1] - 3*Lambda[2]
            Family ()
            4*Lambda[1] - 5*Lambda[2]
            Family ()
            -Lambda[1] + 2*Lambda[2]
            Family (f[-alpha[1]]*u[Lambda[1] + Lambda[2]],)
            0
            Family (f[-alpha[1] - alpha[2]]*u[Lambda[1] + Lambda[2]], f[-alpha[2]]*f[-alpha[1]]*u[Lambda[1] + Lambda[2]])
            Lambda[1] - 2*Lambda[2]
            Family (f[-alpha[2]]*f[-alpha[1] - alpha[2]]*u[Lambda[1] + Lambda[2]],)
            2*Lambda[1] - 4*Lambda[2]
            Family ()
            -3*Lambda[1] + 3*Lambda[2]
            Family ()
            -2*Lambda[1] + Lambda[2]
            Family (f[-alpha[1]]*f[-alpha[1] - alpha[2]]*u[Lambda[1] + Lambda[2]],)
            -Lambda[1] - Lambda[2]
            Family (f[-alpha[1] - alpha[2]]^2*u[Lambda[1] + Lambda[2]],)
            -3*Lambda[2]
            Family ()
            -5*Lambda[1] + 4*Lambda[2]
            Family ()
            -4*Lambda[1] + 2*Lambda[2]
            Family ()
            -3*Lambda[1]
            Family ()
            -2*Lambda[1] - 2*Lambda[2]
            Family ()
        """
        return Family([self.monomial(b) for b in self._indices.weight_space_basis(mu)])

    weight_space_basis = homogeneous_component_basis

    class Element(CombinatorialFreeModule.Element):
        def _acted_upon_(self, scalar, self_on_left=True):
            r"""
            Return the action of ``scalar`` on ``self``.

            EXAMPLES::

                sage: g = LieAlgebra(QQ, cartan_type=['A', 2])
                sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
                sage: L = g.simple_module(La[1] + La[2])
                sage: v = L.highest_weight_vector(); v
                u[Lambda[1] + Lambda[2]]
                sage: f1, f2 = g.pbw_basis().f()
                sage: 5 * v
                5*u[Lambda[1] + Lambda[2]]
                sage: f1 * f2 * v
                f[-alpha[2]]*f[-alpha[1]]*u[Lambda[1] + Lambda[2]]
                 + f[-alpha[1] - alpha[2]]*u[Lambda[1] + Lambda[2]]
                sage: f2 * f1 * v
                -f[-alpha[2]]*f[-alpha[1]]*u[Lambda[1] + Lambda[2]]
                 + 2*f[-alpha[1] - alpha[2]]*u[Lambda[1] + Lambda[2]]
                sage: f2 * f2 * f1 * v
                -2*f[-alpha[2]]*f[-alpha[1] - alpha[2]]*u[Lambda[1] + Lambda[2]]
                sage: f1 * f2 * f1 * v
                f[-alpha[1]]*f[-alpha[1] - alpha[2]]*u[Lambda[1] + Lambda[2]]
                sage: f2 * f1 * f2 * f1 * v
                f[-alpha[1] - alpha[2]]^2*u[Lambda[1] + Lambda[2]]
                sage: f1 * f2 * f2 * f1 * v
                2*f[-alpha[1] - alpha[2]]^2*u[Lambda[1] + Lambda[2]]
            """
            # check for scalars first
            ret = CombinatorialFreeModule.Element._acted_upon_(self, scalar, self_on_left)
            if ret is not None:
                return ret

            if self_on_left:  # this is a left module action
                return None

            P = self.parent()
            return P.retract(scalar * P.lift(self))

        _lmul_ = _acted_upon_
        _rmul_ = _acted_upon_


class FiniteDimensionalSimpleModule(SimpleModule):
    """
    A finite dimensional simple module.
    """
    def bgg_resolution(self):
        """
        Return the BGG resolution of ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['A', 2])
            sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: L = g.simple_module(La[1] + La[2])
            sage: L.bgg_resolution()
            BGG resolution of Simple module with highest weight Lambda[1] + Lambda[2]
             of Lie algebra of ['A', 2] in the Chevalley basis
        """
        from sage.algebras.lie_algebras.bgg_resolution import BGGResolution
        return BGGResolution(self)
