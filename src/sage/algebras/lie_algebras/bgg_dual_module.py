r"""
BGG Category O Dual Modules

AUTHORS:

- Travis Scrimshaw (2024-01-07): Initial version
"""

#*****************************************************************************
#       Copyright (C) 2023 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.cachefunc import cached_method
from sage.categories.modules import Modules
from sage.categories.morphism import Morphism
from sage.categories.homset import Hom, Homset
from sage.monoids.indexed_free_monoid import IndexedFreeAbelianMonoid
from sage.combinat.free_module import CombinatorialFreeModule

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
        [v[2*Lambda[1]],
         f[-alpha[1]]*v[2*Lambda[1]],
         f[-alpha[1]]^2*v[2*Lambda[1]],
         f[-alpha[1]]^3*v[2*Lambda[1]],
         f[-alpha[1]]^4*v[2*Lambda[1]],
         f[-alpha[1]]^5*v[2*Lambda[1]],
         f[-alpha[1]]^6*v[2*Lambda[1]]]
        sage: e, h, f = g.pbw_basis().algebra_generators()
        sage: [f * vec for vec in elts]
        [2*f[-alpha[1]]*v[2*Lambda[1]],
         2*f[-alpha[1]]^2*v[2*Lambda[1]],
         0,
         -4*f[-alpha[1]]^4*v[2*Lambda[1]],
         -10*f[-alpha[1]]^5*v[2*Lambda[1]],
         -18*f[-alpha[1]]^6*v[2*Lambda[1]],
         -28*f[-alpha[1]]^7*v[2*Lambda[1]]]
        sage: [e * vec for vec in elts]
        [0,
         v[2*Lambda[1]],
         f[-alpha[1]]*v[2*Lambda[1]],
         f[-alpha[1]]^2*v[2*Lambda[1]],
         f[-alpha[1]]^3*v[2*Lambda[1]],
         f[-alpha[1]]^4*v[2*Lambda[1]],
         f[-alpha[1]]^5*v[2*Lambda[1]]]
        sage: [h * vec for vec in elts]
        [2*v[2*Lambda[1]],
         0,
         -2*f[-alpha[1]]^2*v[2*Lambda[1]],
         -4*f[-alpha[1]]^3*v[2*Lambda[1]],
         -6*f[-alpha[1]]^4*v[2*Lambda[1]],
         -8*f[-alpha[1]]^5*v[2*Lambda[1]],
         -10*f[-alpha[1]]^6*v[2*Lambda[1]]]
    """
    def __init__(self, module):
        r"""
        Initialize ``self``.
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
        return self._module._repr_generator(m)

    def _latex_generator(self, m):
        return self._module._latex_generator(m)

    _repr_term = _repr_generator
    _latex_term = _latex_generator

    def degree_on_basis(self, m):
        return self._module.degree_on_basis(m)

    def highest_weight(self):
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
            sage: Mc
            BGG Dual of Verma module with highest weight 2*Lambda[1] of
             Lie algebra of ['A', 1] in the Chevalley basis
        """
        hwv = self._module.highest_weight_vector()
        return self.element_class(self, hwv.monomial_coefficients(copy=False))

    def lie_algebra(self):
        return self._g

    def dual(self):
        return self._module

    @cached_method
    def _lie_algebra_on_basis(self, b, m):
        r"""
        Return the action of the Lie algebra basis element indexed by ``b``
        on the basis element of ``self`` indexed by ``m``.
        """
        al = self._g.degree_on_basis(b)
        wt = self.degree_on_basis(m)
        if al == 0:  # b is indexing part of the Cartan subalgebra
            # We are assuming b is part of the coroot lattice.
            # FIXME: Add something at the category level to return this.
            ac = b
            return self.term(m, wt.scalar(ac))

        gens = self._module.homogeneous_component_basis(wt + al)
        elt = self._g.basis()[b]
        vecs = {g.leading_support(): elt.transpose() * g for g in gens}
        return self.element_class(self, {k: c for k, v in vecs.items() if (c := v[m])})

    def _pbw_monomial_on_basis(self, p, m):
        r"""
        Return the action of the PBW monomial indexed by ``p`` on the basis
        element of ``self`` indexed by ``m``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['A', 1])
            sage: La = g.cartan_type().root_system().weight_lattice().fundamental_weights()
            sage: M = g.verma_module(2*La[1])
            sage: Mc = M.dual()
            sage: B = Mc.basis()
            sage: it = iter(B)
            sage: elts = [next(it) for _ in range(7)]; elts
            [v[2*Lambda[1]],
             f[-alpha[1]]*v[2*Lambda[1]],
             f[-alpha[1]]^2*v[2*Lambda[1]],
             f[-alpha[1]]^3*v[2*Lambda[1]],
             f[-alpha[1]]^4*v[2*Lambda[1]],
             f[-alpha[1]]^5*v[2*Lambda[1]],
             f[-alpha[1]]^6*v[2*Lambda[1]]]
            sage: e, h, f = g.pbw_basis().algebra_generators()
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

class SimpleModule(CombinatorialFreeModule):
    r"""
    Return the simple module `L_{\lambda}` as the image of the natural
    morphism `\phi \colom M_{\lambda} \to M_{\lambda}^{\vee}`.
    """
