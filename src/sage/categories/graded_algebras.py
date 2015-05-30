r"""
Graded Algebras
"""
#*****************************************************************************
#  Copyright (C) 2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.graded_modules import GradedModulesCategory

class GradedAlgebras(GradedModulesCategory):
    """
    The category of graded algebras

    EXAMPLES::

        sage: GradedAlgebras(ZZ)
        Category of graded algebras over Integer Ring
        sage: GradedAlgebras(ZZ).super_categories()
        [Category of algebras over Integer Ring,
         Category of graded modules over Integer Ring]

    TESTS::

        sage: TestSuite(GradedAlgebras(ZZ)).run()
    """

    class ParentMethods:
        pass

    class ElementMethods:
        pass
