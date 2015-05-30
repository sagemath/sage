r"""
Partially ordered monoids
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_singleton import Category_singleton
from sage.categories.basic import Posets, Monoids

class PartiallyOrderedMonoids(Category_singleton):
    """
    The category of partially ordered monoids, that is partially ordered sets
    which are also monoids, and such that multiplication preserves the
    ordering: `x \leq y` implies `x*z < y*z` and `z*x < z*y`.

    http://en.wikipedia.org/wiki/Ordered_monoid

    EXAMPLES::

        sage: PartiallyOrderedMonoids()
        Category of partially ordered monoids
        sage: PartiallyOrderedMonoids().super_categories()
        [Category of posets, Category of monoids]

    TESTS::

        sage: TestSuite(PartiallyOrderedMonoids()).run()
    """

    def super_categories(self):
        """
        EXAMPLES::

            sage: PartiallyOrderedMonoids().super_categories()
            [Category of posets, Category of monoids]
        """
        return [Posets(), Monoids()]

    class ParentMethods:
        pass

    class ElementMethods:
        pass
