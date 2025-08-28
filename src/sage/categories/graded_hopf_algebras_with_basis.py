# sage_setup: distribution = sagemath-categories
r"""
Graded Hopf algebras with basis
"""
#*****************************************************************************
#  Copyright (C) 2008      Teresa  Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.tensor import tensor
from sage.categories.graded_modules import GradedModulesCategory
from sage.categories.with_realizations import WithRealizationsCategory
from sage.misc.cachefunc import cached_method


class GradedHopfAlgebrasWithBasis(GradedModulesCategory):
    """
    The category of graded Hopf algebras with a distinguished basis.

    EXAMPLES::

        sage: C = GradedHopfAlgebrasWithBasis(ZZ); C
        Category of graded Hopf algebras with basis over Integer Ring
        sage: C.super_categories()
        [Category of filtered Hopf algebras with basis over Integer Ring,
         Category of graded algebras with basis over Integer Ring,
         Category of graded coalgebras with basis over Integer Ring]

        sage: C is HopfAlgebras(ZZ).WithBasis().Graded()
        True
        sage: C is HopfAlgebras(ZZ).Graded().WithBasis()
        False

    TESTS::

        sage: TestSuite(C).run()
    """
    def example(self):
        """
        Return an example of a graded Hopf algebra with
        a distinguished basis.

        TESTS::

            sage: GradedHopfAlgebrasWithBasis(QQ).example()                             # needs sage.modules
            An example of a graded connected Hopf algebra with basis over Rational Field
        """
        from sage.categories.examples.graded_connected_hopf_algebras_with_basis import \
            GradedConnectedCombinatorialHopfAlgebraWithPrimitiveGenerator
        return GradedConnectedCombinatorialHopfAlgebraWithPrimitiveGenerator(self.base())

    class ParentMethods:
        pass

    class ElementMethods:
        pass

    class WithRealizations(WithRealizationsCategory):
        @cached_method
        def super_categories(self):
            """
            EXAMPLES::

                sage: GradedHopfAlgebrasWithBasis(QQ).WithRealizations().super_categories()
                [Join of Category of Hopf algebras over Rational Field
                 and Category of graded algebras over Rational Field
                 and Category of graded coalgebras over Rational Field]

            TESTS::

                sage: TestSuite(GradedHopfAlgebrasWithBasis(QQ).WithRealizations()).run()
            """
            from sage.categories.graded_hopf_algebras import GradedHopfAlgebras
            R = self.base_category().base_ring()
            return [GradedHopfAlgebras(R)]

    class Connected(CategoryWithAxiom_over_base_ring):
        def example(self):
            """
            Return an example of a graded connected Hopf algebra with
            a distinguished basis.

            TESTS::

                sage: GradedHopfAlgebrasWithBasis(QQ).Connected().example()             # needs sage.modules
                An example of a graded connected Hopf algebra with basis over Rational Field
            """
            from sage.categories.examples.graded_connected_hopf_algebras_with_basis import \
                GradedConnectedCombinatorialHopfAlgebraWithPrimitiveGenerator
            return GradedConnectedCombinatorialHopfAlgebraWithPrimitiveGenerator(self.base())

        class ParentMethods:
            def counit_on_basis(self, i):
                r"""
                The default counit of a graded connected Hopf algebra.

                INPUT:

                - ``i`` -- an element of the index set

                OUTPUT: an element of the base ring

                .. MATH::

                    c(i) := \begin{cases}
                    1 & \hbox{if $i$ indexes the $1$ of the algebra}\\
                    0 & \hbox{otherwise}.
                    \end{cases}

                EXAMPLES::

                    sage: H = GradedHopfAlgebrasWithBasis(QQ).Connected().example()     # needs sage.modules
                    sage: H.monomial(4).counit()  # indirect doctest                    # needs sage.modules
                    0
                    sage: H.monomial(0).counit()  # indirect doctest                    # needs sage.modules
                    1
                """
                if i == self.one_basis():
                    return self.base_ring().one()
                return self.base_ring().zero()

            @cached_method
            def antipode_on_basis(self, index):
                r"""
                The antipode on the basis element indexed by ``index``.

                INPUT:

                - ``index`` -- an element of the index set

                For a graded connected Hopf algebra, we can define
                an antipode recursively by

                .. MATH::

                    S(x) := -\sum_{x^L \neq x} S(x^L) \times x^R

                when `|x| > 0`, and by `S(x) = x` when `|x| = 0`.

                TESTS::

                    sage: # needs sage.modules
                    sage: H = GradedHopfAlgebrasWithBasis(QQ).Connected().example()
                    sage: H.monomial(0).antipode()  # indirect doctest
                    P0
                    sage: H.monomial(1).antipode()  # indirect doctest
                    -P1
                    sage: H.monomial(2).antipode()  # indirect doctest
                    P2
                    sage: H.monomial(3).antipode()  # indirect doctest
                    -P3
                """
                if self.monomial(index) == self.one():
                    return self.one()

                S = self.antipode_on_basis
                x__S_Id = tensor([self, self]).module_morphism(
                    lambda ab: S(ab[0]) * self.monomial(ab[1]),
                    codomain=self)
                return -x__S_Id(
                    self.monomial(index).coproduct()
                    - tensor([self.monomial(index), self.one()])
                )

        class ElementMethods:
            pass
