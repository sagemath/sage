r"""
Graded Algebras
"""
#*****************************************************************************
#  Copyright (C) 2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.graded_modules import GradedModulesCategory
from sage.categories.signed_tensor import SignedTensorProductsCategory
from sage.misc.cachefunc import cached_method

class GradedAlgebras(GradedModulesCategory):
    """
    The category of graded algebras

    EXAMPLES::

        sage: GradedAlgebras(ZZ)
        Category of graded algebras over Integer Ring
        sage: GradedAlgebras(ZZ).super_categories()
        [Category of filtered algebras over Integer Ring,
         Category of graded modules over Integer Ring]

    TESTS::

        sage: TestSuite(GradedAlgebras(ZZ)).run()
    """
    class ParentMethods:
        def graded_algebra(self):
            """
            Return the associated graded algebra to ``self``.

            Since ``self`` is already graded, this just returns
            ``self``.

            EXAMPLES::

                sage: m = SymmetricFunctions(QQ).m()
                sage: m.graded_algebra() is m
                True
            """
            return self

    class ElementMethods:
        pass

    class SubcategoryMethods:
        def SignedTensorProducts(self):
            r"""
            Return the full subcategory of objects of ``self`` constructed
            as signed tensor products.

            .. SEEALSO::

                - :class:`~sage.categories.signed_tensor.SignedTensorProductsCategory`
                - :class:`~.covariant_functorial_construction.CovariantFunctorialConstruction`

            EXAMPLES::

                sage: AlgebrasWithBasis(QQ).Graded().SignedTensorProducts()
                Category of signed tensor products of graded algebras with basis
                 over Rational Field
            """
            return SignedTensorProductsCategory.category_of(self)

    class SignedTensorProducts(SignedTensorProductsCategory):
        @cached_method
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: Algebras(QQ).Graded().SignedTensorProducts().extra_super_categories()
                [Category of graded algebras over Rational Field]
                sage: Algebras(QQ).Graded().SignedTensorProducts().super_categories()
                [Category of graded algebras over Rational Field]

            Meaning: a signed tensor product of algebras is an algebra
            """
            return [self.base_category()]

