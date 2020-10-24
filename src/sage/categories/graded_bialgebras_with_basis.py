r"""
Graded bialgebras with basis
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

def GradedBialgebrasWithBasis(base_ring):
    """
    The category of graded bialgebras with a distinguished basis

    EXAMPLES::

        sage: C = GradedBialgebrasWithBasis(QQ); C
        Join of Category of ...
        sage: sorted(C.super_categories(), key=str)
        [Category of bialgebras with basis over Rational Field,
         Category of graded algebras with basis over Rational Field,
         Category of graded coalgebras with basis over Rational Field]

    TESTS::

        sage: TestSuite(C).run()
    """
    from sage.categories.all import BialgebrasWithBasis
    return BialgebrasWithBasis(base_ring).Graded()
