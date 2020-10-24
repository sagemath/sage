r"""
Commutative algebras
"""
#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.algebras import Algebras

class CommutativeAlgebras(CategoryWithAxiom_over_base_ring):
    """
    The category of commutative algebras with unit over a given base ring.

    EXAMPLES::

        sage: M = CommutativeAlgebras(GF(19))
        sage: M
        Category of commutative algebras over Finite Field of size 19
        sage: CommutativeAlgebras(QQ).super_categories()
        [Category of algebras over Rational Field, Category of commutative rings]

    This is just a shortcut for::

        sage: Algebras(QQ).Commutative()
        Category of commutative algebras over Rational Field

    TESTS::

        sage: Algebras(QQ).Commutative() is CommutativeAlgebras(QQ)
        True
        sage: TestSuite(CommutativeAlgebras(ZZ)).run()

    Todo:

     - product   ( = Cartesian product)
     - coproduct ( = tensor product over base ring)
    """

    def __contains__(self, A):
        """
        EXAMPLES::

            sage: QQ['a'] in CommutativeAlgebras(QQ)
            True
            sage: QQ['a,b'] in CommutativeAlgebras(QQ)
            True
            sage: FreeAlgebra(QQ,2,'a,b') in CommutativeAlgebras(QQ)
            False

        TODO: get rid of this method once all commutative algebras in
        Sage declare themselves in this category
        """
        return super(CommutativeAlgebras, self).__contains__(A) or \
            (A in Algebras(self.base_ring()) and hasattr(A, "is_commutative") and A.is_commutative())
