r"""
Matrix algebras
"""
#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from category_types import Category_over_base_ring
from algebras import Algebras

class MatrixAlgebras(Category_over_base_ring):
    """
    The category of matrix algebras over a field.

    EXAMPLES::

        sage: MatrixAlgebras(RationalField())
        Category of matrix algebras over Rational Field

    TESTS::

        sage: TestSuite(MatrixAlgebras(ZZ)).run()
    """

    def super_categories(self):
        """
        EXAMPLES::

            sage: MatrixAlgebras(QQ).super_categories()
            [Category of algebras over Rational Field]
        """
        R = self.base_ring()
        return [Algebras(R)]
