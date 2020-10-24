"""
Category of chain complexes
"""

#*****************************************************************************
#  Copyright (C) 2007 Robert Bradshaw
#                2009 Mike Hansen
#                2013 Volker Braun
#                2013, 2015 Travis Scrimshaw
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.category_types import Category_module

# TODO: make this into a better category
#############################################################
# ChainComplex
#############################################################
class ChainComplexes(Category_module):
    """
    The category of all chain complexes over a base ring.

    EXAMPLES::

        sage: ChainComplexes(RationalField())
        Category of chain complexes over Rational Field

        sage: ChainComplexes(Integers(9))
        Category of chain complexes over Ring of integers modulo 9

     TESTS::

        sage: TestSuite(ChainComplexes(RationalField())).run()
    """

    def super_categories(self):
        """
        EXAMPLES::

            sage: ChainComplexes(Integers(9)).super_categories()
            [Category of modules over Ring of integers modulo 9]
        """
        from sage.categories.all import Fields, Modules, VectorSpaces
        base_ring = self.base_ring()
        if base_ring in Fields():
            return [VectorSpaces(base_ring)]
        return [Modules(base_ring)]
