r"""
Graded coalgebras with basis
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

def GradedCoalgebrasWithBasis(base_ring):
    """
    The category of graded coalgebras with a distinguished basis

    EXAMPLES::

        sage: C = GradedCoalgebrasWithBasis(QQ); C
        Join of Category of graded modules with basis over Rational Field
            and Category of coalgebras with basis over Rational Field
        sage: C is Coalgebras(QQ).WithBasis().Graded()
        True

    TESTS::

        sage: TestSuite(C).run()
    """
    from sage.categories.all import CoalgebrasWithBasis
    return CoalgebrasWithBasis(base_ring).Graded()
