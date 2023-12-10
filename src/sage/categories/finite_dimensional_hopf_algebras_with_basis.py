# sage_setup: distribution = sagemath-categories
r"""
Finite dimensional Hopf algebras with basis
"""
# ****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring


class FiniteDimensionalHopfAlgebrasWithBasis(CategoryWithAxiom_over_base_ring):
    """
    The category of finite dimensional Hopf algebras with a
    distinguished basis.

    EXAMPLES::

        sage: FiniteDimensionalHopfAlgebrasWithBasis(QQ)
        Category of finite dimensional Hopf algebras with basis over Rational Field
        sage: FiniteDimensionalHopfAlgebrasWithBasis(QQ).super_categories()
        [Category of Hopf algebras with basis over Rational Field,
         Category of finite dimensional algebras with basis over Rational Field]

    TESTS::

        sage: TestSuite(FiniteDimensionalHopfAlgebrasWithBasis(ZZ)).run()
    """

    class ParentMethods:
        pass

    class ElementMethods:
        pass
