r"""
Right modules
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from .category_types import Category_over_base_ring
from sage.categories.commutative_additive_groups import CommutativeAdditiveGroups

##?class RightModules(Category_over_base_rng):
class RightModules(Category_over_base_ring):
    """
    The category of right modules
    right modules over an rng (ring not necessarily with unit), i.e.
    an abelian group with right multiplication by elements of the rng

    EXAMPLES::

        sage: RightModules(QQ)
        Category of right modules over Rational Field
        sage: RightModules(QQ).super_categories()
        [Category of commutative additive groups]

    TESTS::

        sage: TestSuite(RightModules(ZZ)).run()
    """

    def super_categories(self):
        """
        EXAMPLES::

            sage: RightModules(QQ).super_categories()
            [Category of commutative additive groups]
        """
        return [CommutativeAdditiveGroups()]

    class ParentMethods:
        pass

    class ElementMethods:
        ## x * r
        pass
