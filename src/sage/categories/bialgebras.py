r"""
Bialgebras
"""
# ****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#  Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.categories.category_types import Category_over_base_ring
from sage.categories.all import Algebras, Coalgebras
from sage.categories.super_modules import SuperModulesCategory
from sage.misc.lazy_import import LazyImport


class Bialgebras(Category_over_base_ring):
    """
    The category of bialgebras

    EXAMPLES::

        sage: Bialgebras(ZZ)
        Category of bialgebras over Integer Ring
        sage: Bialgebras(ZZ).super_categories()
        [Category of algebras over Integer Ring, Category of coalgebras over Integer Ring]

    TESTS::

        sage: TestSuite(Bialgebras(ZZ)).run()
    """

    def super_categories(self):
        """
        EXAMPLES::

            sage: Bialgebras(QQ).super_categories()
            [Category of algebras over Rational Field, Category of coalgebras over Rational Field]
        """
        R = self.base_ring()
        return [Algebras(R), Coalgebras(R)]

    def additional_structure(self):
        r"""
        Return ``None``.

        Indeed, the category of bialgebras defines no additional
        structure: a morphism of coalgebras and of algebras between
        two bialgebras is a bialgebra morphism.

        .. SEEALSO:: :meth:`Category.additional_structure`

        .. TODO:: This category should be a :class:`CategoryWithAxiom`.

        EXAMPLES::

            sage: Bialgebras(QQ).additional_structure()
        """
        return None

    class ElementMethods:

        def is_primitive(self):
            """
            Return whether ``self`` is a primitive element.

            EXAMPLES::

                sage: s = SymmetricFunctions(QQ).schur()
                sage: s([5]).is_primitive()
                False
                sage: p = SymmetricFunctions(QQ).powersum()
                sage: p([5]).is_primitive()
                True
            """
            one = self.parent().one()
            return self.coproduct() == one.tensor(self) + self.tensor(one)

        def is_grouplike(self):
            """
            Return whether ``self`` is a grouplike element.

            EXAMPLES::

                sage: s = SymmetricFunctions(QQ).schur()
                sage: s([5]).is_grouplike()
                False
                sage: s([]).is_grouplike()
                True
            """
            return self.coproduct() == self.tensor(self)

    class Super(SuperModulesCategory):
        pass

    WithBasis = LazyImport('sage.categories.bialgebras_with_basis', 'BialgebrasWithBasis')
