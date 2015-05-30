r"""
Finite dimensional coalgebras with basis
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

def FiniteDimensionalCoalgebrasWithBasis(base_ring):
    """
    The category of finite dimensional coalgebras with a distinguished basis

    EXAMPLES::

        sage: C = FiniteDimensionalCoalgebrasWithBasis(QQ); C
        Category of finite dimensional coalgebras with basis over Rational Field
        sage: sorted(C.super_categories(), key=str)
        [Category of coalgebras with basis over Rational Field,
         Category of finite dimensional modules with basis over Rational Field]
        sage: C is Coalgebras(QQ).WithBasis().FiniteDimensional()
        True

    TESTS::

        sage: TestSuite(C).run()
    """
    from sage.categories.all import CoalgebrasWithBasis
    return CoalgebrasWithBasis(base_ring).FiniteDimensional()
