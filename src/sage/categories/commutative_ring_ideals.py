# sage_setup: distribution = sagemath-categories
r"""
Commutative ring ideals
"""
#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_types import Category_ideal
from sage.categories.commutative_rings import CommutativeRings
from sage.categories.ring_ideals import RingIdeals


class CommutativeRingIdeals(Category_ideal):
    """
    The category of ideals in a fixed commutative ring.

    EXAMPLES::

        sage: C = CommutativeRingIdeals(IntegerRing())
        sage: C
        Category of commutative ring ideals in Integer Ring
    """
    def __init__(self, R):
        """
        EXAMPLES::

            sage: CommutativeRingIdeals(ZZ)
            Category of commutative ring ideals in Integer Ring
            sage: CommutativeRingIdeals(IntegerModRing(4))
            Category of commutative ring ideals in Ring of integers modulo 4

        TESTS::

            sage: CommutativeRingIdeals(Partitions(4))                                  # needs sage.combinat
            Traceback (most recent call last):
            ...
            TypeError: R (=Partitions of the integer 4) must be a commutative ring
            sage: TestSuite(CommutativeRingIdeals(ZZ)).run()
        """
        if R not in CommutativeRings():
            raise TypeError("R (=%s) must be a commutative ring" % R)
        Category_ideal.__init__(self, R)

    def super_categories(self):
        """
        EXAMPLES::

            sage: CommutativeRingIdeals(ZZ).super_categories()
            [Category of ring ideals in Integer Ring]
        """
        R = self.ring()
        return [RingIdeals(R)]
