r"""
Commutative additive groups
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_types import AbelianCategory
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.algebra_functor import AlgebrasCategory
from sage.categories.cartesian_product import CartesianProductsCategory
from sage.categories.additive_groups import AdditiveGroups

class CommutativeAdditiveGroups(CategoryWithAxiom, AbelianCategory):
    """
    The category of abelian groups, i.e. additive abelian monoids
    where each element has an inverse.

    EXAMPLES::

        sage: C = CommutativeAdditiveGroups(); C
        Category of commutative additive groups
        sage: C.super_categories()
        [Category of additive groups, Category of commutative additive monoids]
        sage: sorted(C.axioms())
        ['AdditiveAssociative', 'AdditiveCommutative', 'AdditiveInverse', 'AdditiveUnital']
        sage: C is CommutativeAdditiveMonoids().AdditiveInverse()
        True
        sage: from sage.categories.additive_groups import AdditiveGroups
        sage: C is AdditiveGroups().AdditiveCommutative()
        True

    .. NOTE::

        This category is currently empty. It's left there for backward
        compatibility and because it is likely to grow in the future.

    TESTS::

        sage: TestSuite(CommutativeAdditiveGroups()).run()
        sage: sorted(CommutativeAdditiveGroups().CartesianProducts().axioms())
        ['AdditiveAssociative', 'AdditiveCommutative', 'AdditiveInverse', 'AdditiveUnital']

    The empty covariant functorial construction category classes
    ``CartesianProducts`` and ``Algebras`` are left here for the sake
    of nicer output since this is a commonly used category::

        sage: CommutativeAdditiveGroups().CartesianProducts()
        Category of Cartesian products of commutative additive groups
        sage: CommutativeAdditiveGroups().Algebras(QQ)
        Category of commutative additive group algebras over Rational Field

    Also, it's likely that some code will end up there at some point.
    """
    _base_category_class_and_axiom = (AdditiveGroups, "AdditiveCommutative")

    class CartesianProducts(CartesianProductsCategory):
        pass

    class Algebras(AlgebrasCategory):
        pass
