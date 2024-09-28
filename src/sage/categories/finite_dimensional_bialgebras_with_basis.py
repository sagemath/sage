# sage_setup: distribution = sagemath-categories
r"""
Finite dimensional bialgebras with basis
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************


def FiniteDimensionalBialgebrasWithBasis(base_ring):
    """
    The category of finite dimensional bialgebras with a distinguished basis.

    EXAMPLES::

        sage: C = FiniteDimensionalBialgebrasWithBasis(QQ); C
        Category of finite dimensional bialgebras with basis over Rational Field
        sage: sorted(C.super_categories(), key=str)
        [Category of bialgebras with basis over Rational Field,
         Category of finite dimensional algebras with basis over Rational Field]
        sage: C is Bialgebras(QQ).WithBasis().FiniteDimensional()
        True

    TESTS::

        sage: TestSuite(C).run()
    """
    from sage.categories.bialgebras_with_basis import BialgebrasWithBasis
    return BialgebrasWithBasis(base_ring).FiniteDimensional()
