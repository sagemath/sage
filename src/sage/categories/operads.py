r"""
Operads
"""
# ****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.category_types import Category_over_base_ring
from sage.categories.set_operads import SetOperads
from sage.categories.modules import Modules
from sage.misc.cachefunc import cached_method


class Operads(Category_over_base_ring):
    """
    The category of operads.

    EXAMPLES::

        sage: Operads(ZZ)
        Category of operads over Integer Ring
        sage: Operads(ZZ).super_categories()
        [Category of modules over Integer Ring, Category of set operads]

    TESTS::

        sage: C = Operads(ZZ)
        sage: TestSuite(C).run()
    """
    @cached_method
    def super_categories(self):
        """
        Return the super-categories.

        EXAMPLES::

            sage: Operads(QQ).super_categories()
            [Category of vector spaces over Rational Field,
            Category of set operads]
        """
        return [Modules(self.base_ring()), SetOperads()]
