# sage_setup: distribution = sagemath-categories
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

from sage.categories.algebra_ideals import AlgebraIdeals
from sage.categories.category_types import Category_ideal, Category_in_ambient
from sage.categories.commutative_algebras import CommutativeAlgebras
from sage.categories.commutative_rings import CommutativeRings


class CommutativeAlgebraIdeals(Category_ideal):
    """
    The category of ideals in a fixed commutative algebra `A`.

    EXAMPLES::

        sage: C = CommutativeAlgebraIdeals(QQ['x'])
        sage: C
        Category of commutative algebra ideals in
         Univariate Polynomial Ring in x over Rational Field
    """
    def __init__(self, A):
        """
        EXAMPLES::

            sage: CommutativeAlgebraIdeals(ZZ['x'])
            Category of commutative algebra ideals in
             Univariate Polynomial Ring in x over Integer Ring

            sage: CommutativeAlgebraIdeals(ZZ)
            Traceback (most recent call last):
            ...
            TypeError: A (=Integer Ring) must be a commutative algebra

            sage: CommutativeAlgebraIdeals(IntegerModRing(4))
            Traceback (most recent call last):
            ...
            TypeError: A (=Ring of integers modulo 4) must be a commutative algebra

            sage: CommutativeAlgebraIdeals(Partitions(4))                               # needs sage.combinat
            Traceback (most recent call last):
            ...
            TypeError: A (=Partitions of the integer 4) must be a commutative algebra

        TESTS::

            sage: TestSuite(CommutativeAlgebraIdeals(QQ['x'])).run()
        """
        # TODO: replace by ``A in CommutativeAlgebras(*)`` once a
        # suitable mantra has been implemented for this.
        try:
            base_ring = A.base_ring()
        except AttributeError:
            raise TypeError(f"A (={A}) must be a commutative algebra")
        else:
            if base_ring not in CommutativeRings() or A not in CommutativeAlgebras(base_ring.category()):
                raise TypeError(f"A (={A}) must be a commutative algebra")

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
