r"""
Graded bialgebras
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

def GradedBialgebras(base_ring):
    """
    The category of graded bialgebras

    EXAMPLES::

        sage: C = GradedBialgebras(QQ); C
        Join of Category of graded algebras over Rational Field
            and Category of bialgebras over Rational Field
        sage: C is Bialgebras(QQ).Graded()
        True

    TESTS::

        sage: TestSuite(C).run()
    """
    from sage.categories.all import Bialgebras
    return Bialgebras(base_ring).Graded()
