# sage_setup: distribution = sagemath-categories
r"""
Super Hopf algebras with basis
"""
# ****************************************************************************
#  Copyright (C) 2015 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************
from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.super_modules import SuperModulesCategory


class SuperHopfAlgebrasWithBasis(SuperModulesCategory):
    """
    The category of super Hopf algebras with a distinguished basis.

    EXAMPLES::

        sage: C = HopfAlgebras(ZZ).WithBasis().Super(); C
        Category of super Hopf algebras with basis over Integer Ring
        sage: sorted(C.super_categories(), key=str)
        [Category of super Hopf algebras over Integer Ring,
         Category of super algebras with basis over Integer Ring,
         Category of super coalgebras with basis over Integer Ring]

    TESTS::

        sage: C = HopfAlgebras(ZZ).WithBasis().Super()
        sage: TestSuite(C).run()
    """
    class ParentMethods:
        @lazy_attribute
        def antipode(self):
            """
            The antipode of this Hopf algebra.

            If :meth:`.antipode_basis` is available, this constructs the
            antipode morphism from ``self`` to ``self`` by extending it by
            linearity. Otherwise, :meth:`self.antipode_by_coercion` is used,
            if available.

            EXAMPLES::

                sage: A = SteenrodAlgebra(7)                                            # needs sage.combinat sage.modules
                sage: a = A.an_element()                                                # needs sage.combinat sage.modules
                sage: a, A.antipode(a)                                                  # needs sage.combinat sage.modules
                (6 Q_1 Q_3 P(2,1), Q_1 Q_3 P(2,1))

            TESTS::

                sage: E.<x,y> = ExteriorAlgebra(QQ)                                     # needs sage.modules
                sage: [b.antipode() for b in E.basis()]                                 # needs sage.modules
                [1, -x, -y, x*y]
            """
            if self.antipode_on_basis is not NotImplemented:
                # Should give the information that this is an anti-morphism of algebra
                return self._module_morphism(self.antipode_on_basis,
                                             codomain=self)
            elif hasattr(self, "antipode_by_coercion"):
                return self.antipode_by_coercion

        def _test_antipode(self, **options):
            r"""
            Test the antipode.

            An *antipode* `S` of a (super) Hopf algebra is a linear
            endomorphism of the Hopf algebra that satisfies the
            following conditions (see :wikipedia:`HopfAlgebra`).

            - If `\mu` and `\Delta` denote the product and coproduct of the
              Hopf algebra, respectively, then `S` satisfies

              .. MATH::

                  \mu \circ (S \tensor 1) \circ \Delta = unit \circ counit
                  \mu \circ (1 \tensor S) \circ \Delta = unit \circ counit

            - `S` is an *anti*-homomorphism:

               .. MATH::

                   S(ab) = (-1)^{\deg a \deg b} S(b) S(a)

               for homogeneous `a` and `b`.

            These properties are tested on :meth:`some_elements`.

            TESTS::

                sage: A = SteenrodAlgebra(7)                                            # needs sage.combinat sage.modules
                sage: A._test_antipode()        # long time                             # needs sage.combinat sage.modules
            """
            tester = self._tester(**options)

            S = self.antipode

            IS = lambda x: self.sum(c * self.monomial(t1) * S(self.monomial(t2))
                                for ((t1, t2), c) in x.coproduct())

            SI = lambda x: self.sum(c * S(self.monomial(t1)) * self.monomial(t2)
                                for ((t1, t2), c) in x.coproduct())

            for x in tester.some_elements():
                x_even = x.even_component()
                x_odd = x.odd_component()
                for y in tester.some_elements():
                    y_even = y.even_component()
                    y_odd = y.odd_component()

                    # The antipode is a graded anti-homomorphism.
                    tester.assertEqual(S(x_even) * S(y_even),
                                       S(y_even * x_even))
                    tester.assertEqual(S(x_even) * S(y_odd),
                                       S(y_odd * x_even))
                    tester.assertEqual(S(x_odd) * S(y_even),
                                       S(y_even * x_odd))
                    tester.assertEqual(S(x_odd) * S(y_odd),
                                       -S(y_odd * x_odd))

                # mu * (S # I) * delta == counit * unit
                tester.assertEqual(SI(x), self.counit(x) * self.one())

                # mu * (I # S) * delta == counit * unit
                tester.assertEqual(IS(x), self.counit(x) * self.one())
