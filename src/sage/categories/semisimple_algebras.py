r"""
Semisimple Algebras
"""
#*****************************************************************************
#  Copyright (C) 2011-2015 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.bindable_class import BoundClass
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import LazyImport
from .category_types import Category_over_base_ring
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from .algebras import Algebras

class SemisimpleAlgebras(Category_over_base_ring):
    """
    The category of semisimple algebras over a given base ring.

    EXAMPLES::

        sage: from sage.categories.semisimple_algebras import SemisimpleAlgebras
        sage: C = SemisimpleAlgebras(QQ); C
        Category of semisimple algebras over Rational Field

    This category is best constructed as::

        sage: D = Algebras(QQ).Semisimple(); D
        Category of semisimple algebras over Rational Field
        sage: D is C
        True

        sage: C.super_categories()
        [Category of algebras over Rational Field]

    Typically, finite group algebras are semisimple::

        sage: DihedralGroup(5).algebra(QQ) in SemisimpleAlgebras
        True

    Unless the characteristic of the field divides the order of the group::

        sage: DihedralGroup(5).algebra(IntegerModRing(5)) in SemisimpleAlgebras
        False

        sage: DihedralGroup(5).algebra(IntegerModRing(7)) in SemisimpleAlgebras
        True

    .. SEEALSO:: :wikipedia:`Semisimple_algebra`

    TESTS::

        sage: TestSuite(C).run()
    """
    @staticmethod
    def __classget__(cls, base_category, base_category_class):
        """
        Implement the shorthand ``Algebras(K).Semisimple()`` for ``SemisimpleAlgebras(K)``.

        This magic mimics the syntax of axioms for a smooth transition
        if ``Semisimple`` becomes one.

        EXAMPLES::

            sage: Algebras(QQ).Semisimple()
            Category of semisimple algebras over Rational Field
            sage: Algebras.Semisimple
            <class 'sage.categories.semisimple_algebras.SemisimpleAlgebras'>
        """
        if base_category is None:
            return cls
        return BoundClass(cls, base_category.base_ring())

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: Algebras(QQ).Semisimple().super_categories()
            [Category of algebras over Rational Field]
        """
        R = self.base_ring()
        return [Algebras(R)]

    class ParentMethods:

        def radical_basis(self, **keywords):
            r"""
            Return a basis of the Jacobson radical of this algebra.

            - ``keywords`` -- for compatibility; ignored.

            OUTPUT: the empty list since this algebra is semisimple.

            EXAMPLES::

                sage: A = SymmetricGroup(4).algebra(QQ)
                sage: A.radical_basis()
                ()

            TESTS::

                sage: A.radical_basis.__module__
                'sage.categories.finite_dimensional_semisimple_algebras_with_basis'
            """
            return ()

    class FiniteDimensional(CategoryWithAxiom_over_base_ring):

        WithBasis = LazyImport('sage.categories.finite_dimensional_semisimple_algebras_with_basis', 'FiniteDimensionalSemisimpleAlgebrasWithBasis')
