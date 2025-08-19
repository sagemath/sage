# sage_setup: distribution = sagemath-categories
r"""
Pointed sets
"""
#*****************************************************************************
#  Copyright (C) 2008 David Kohel <kohel@maths.usyd.edu> and
#                     William Stein <wstein@math.ucsd.edu>
#                     Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_singleton import Category_singleton
from sage.categories.sets_cat import Sets


class PointedSets(Category_singleton):
    """
    The category of pointed sets.

    EXAMPLES::

        sage: PointedSets()
        Category of pointed sets

    TESTS::

        sage: TestSuite(PointedSets()).run()
    """
    #def __call__(self, X, pt):
    #    import sage.sets.all
    #    return sage.sets.all.Set(X, pt)

    def super_categories(self):
        """
        EXAMPLES::

            sage: PointedSets().super_categories()
            [Category of sets]
        """
        return [Sets()] # ???
