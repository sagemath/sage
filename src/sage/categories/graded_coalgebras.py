r"""
Graded Coalgebras
"""
#*****************************************************************************
#  Copyright (C) 2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2011-2013 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

def GradedCoalgebras(base_ring):
    """
    The category of graded coalgebras

    EXAMPLES::

        sage: C = GradedCoalgebras(QQ); C
        Join of Category of graded modules over Rational Field
            and Category of coalgebras over Rational Field
        sage: C is Coalgebras(QQ).Graded()
        True

    TESTS::

        sage: TestSuite(C).run()
    """
    from sage.categories.all import Coalgebras
    return Coalgebras(base_ring).Graded()

