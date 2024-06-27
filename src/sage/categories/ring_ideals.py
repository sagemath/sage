# sage_setup: distribution = sagemath-categories
r"""
Ring ideals
"""
# ****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.categories.category_types import Category_ideal
from sage.categories.modules import Modules
from sage.categories.rings import Rings

_Rings = Rings()


class RingIdeals(Category_ideal):
    """
    The category of two-sided ideals in a fixed ring.

    EXAMPLES::

        sage: Ideals(Integers(200))
        Category of ring ideals in Ring of integers modulo 200
        sage: C = Ideals(IntegerRing()); C
        Category of ring ideals in Integer Ring
        sage: I = C([8,12,18])
        sage: I
        Principal ideal (2) of Integer Ring

    See also: :class:`CommutativeRingIdeals`.

    .. TODO::

         - If useful, implement ``RingLeftIdeals`` and ``RingRightIdeals``
           of which ``RingIdeals`` would be a subcategory.

         - Make ``RingIdeals(R)``, return ``CommutativeRingIdeals(R)``
           when ``R`` is commutative.
    """
    def __init__(self, R):
        """
        EXAMPLES::

            sage: RingIdeals(ZZ)
            Category of ring ideals in Integer Ring
            sage: RingIdeals(3)
            Traceback (most recent call last):
            ...
            TypeError: R (=3) must be a ring

        TESTS::

            sage: TestSuite(RingIdeals(ZZ)).run()
        """
        if R not in _Rings:
            raise TypeError("R (=%s) must be a ring" % R)
        Category_ideal.__init__(self, R)

    def super_categories(self):
        """
        EXAMPLES::

            sage: RingIdeals(ZZ).super_categories()
            [Category of modules over Integer Ring]
            sage: RingIdeals(QQ).super_categories()
            [Category of vector spaces over Rational Field]
        """
        R = self.ring()
        return [Modules(R)]
