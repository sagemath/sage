r"""
Graded Hopf algebras
"""
#*****************************************************************************
#  Copyright (C) 2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
#******************************************************************************

def GradedHopfAlgebras(base_ring):
    r"""
    The category of graded Hopf algebras.

    EXAMPLES::

        sage: C = GradedHopfAlgebras(QQ); C
        Join of Category of hopf algebras over Rational Field
            and Category of graded algebras over Rational Field
            and Category of graded coalgebras over Rational Field
        sage: C is HopfAlgebras(QQ).Graded()
        True

    TESTS::

        sage: TestSuite(C).run()

    .. NOTE::

        This is not a graded Hopf algebra as is typically defined
        in algebraic topology as the product in the tensor square
        `(x \otimes y) (a \otimes b) = (xa) \otimes (yb)` does
        not carry an additional sign. For this, instead use
        :class:`super Hopf algebras
        <sage.categories.hopf_algebras.HopfAlgebras.Super>`.
    """
    from sage.categories.all import HopfAlgebras
    return HopfAlgebras(base_ring).Graded()

