r"""
Commutative algebra ideals
"""
#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from .category_types import Category_ideal, Category_in_ambient
from .algebra_ideals import AlgebraIdeals

class CommutativeAlgebraIdeals(Category_ideal):
    """
    The category of ideals in a fixed commutative algebra `A`.

    EXAMPLES::

        sage: C = CommutativeAlgebraIdeals(QQ['x'])
        sage: C
        Category of commutative algebra ideals in Univariate Polynomial Ring in x over Rational Field
    """
    def __init__(self, A):
        """
        EXAMPLES::

            sage: CommutativeAlgebraIdeals(ZZ['x'])
            Category of commutative algebra ideals in Univariate Polynomial Ring in x over Integer Ring

            sage: CommutativeAlgebraIdeals(ZZ)
            Traceback (most recent call last):
            ...
            TypeError: A (=Integer Ring) must be a commutative algebra

            sage: CommutativeAlgebraIdeals(IntegerModRing(4))
            Traceback (most recent call last):
            ...
            TypeError: A (=Ring of integers modulo 4) must be a commutative algebra

            sage: CommutativeAlgebraIdeals(Partitions(4))
            Traceback (most recent call last):
            ...
            TypeError: A (=Partitions of the integer 4) must be a commutative algebra

        TESTS::

            sage: TestSuite(CommutativeAlgebraIdeals(QQ['x'])).run()
        """
        # TODO: replace by ``A in CommutativeAlgebras(*)`` once a
        # suitable mantra has been implemented for this.
        from sage.algebras.algebra import is_Algebra
        from sage.rings.ring import CommutativeRing
        if not (is_Algebra(A) and isinstance(A, CommutativeRing)):
            raise TypeError("A (=%s) must be a commutative algebra"%A)
        Category_in_ambient.__init__(self, A)

    def algebra(self):
        """
        EXAMPLES::

            sage: CommutativeAlgebraIdeals(QQ['x']).algebra()
            Univariate Polynomial Ring in x over Rational Field
        """
        return self.ambient()

    def super_categories(self):
        """
        EXAMPLES::

            sage: CommutativeAlgebraIdeals(QQ['x']).super_categories()
            [Category of algebra ideals in Univariate Polynomial Ring in x over Rational Field]
        """
        R = self.algebra()
        return [AlgebraIdeals(R)]
