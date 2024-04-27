r"""
Commutative additive monoids
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.additive_monoids import AdditiveMonoids

class CommutativeAdditiveMonoids(CategoryWithAxiom):
    """
    The category of commutative additive monoids, that is abelian
    additive semigroups with a unit

    EXAMPLES::

        sage: C = CommutativeAdditiveMonoids(); C
        Category of commutative additive monoids
        sage: C.super_categories()
        [Category of additive monoids, Category of commutative additive semigroups]
        sage: sorted(C.axioms())
        ['AdditiveAssociative', 'AdditiveCommutative', 'AdditiveUnital']
        sage: C is AdditiveMagmas().AdditiveAssociative().AdditiveCommutative().AdditiveUnital()
        True

    .. NOTE::

        This category is currently empty and only serves as a place
        holder to make ``C.example()`` work.

    TESTS::

        sage: TestSuite(CommutativeAdditiveMonoids()).run()
    """
    _base_category_class_and_axiom = (AdditiveMonoids, "AdditiveCommutative")
