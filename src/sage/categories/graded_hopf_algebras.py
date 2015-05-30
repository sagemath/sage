r"""
Graded Hopf algebras
"""
#*****************************************************************************
#  Copyright (C) 2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

def GradedHopfAlgebras(base_ring):
    """
    The category of graded Hopf algebras.

    EXAMPLES::

        sage: C = GradedHopfAlgebras(QQ); C
        Join of Category of hopf algebras over Rational Field
            and Category of graded algebras over Rational Field
        sage: C is HopfAlgebras(QQ).Graded()
        True

    TESTS::

        sage: TestSuite(C).run()
    """
    from sage.categories.all import HopfAlgebras
    return HopfAlgebras(base_ring).Graded()

