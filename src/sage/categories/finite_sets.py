r"""
Finite sets
"""
#*****************************************************************************
#  Copyright (C) 2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.algebra_functor import AlgebrasCategory
from sage.categories.subquotients import SubquotientsCategory

class FiniteSets(CategoryWithAxiom):
    r"""
    The category of finite sets.

    EXAMPLES::

        sage: C = FiniteSets(); C
        Category of finite sets
        sage: C.super_categories()
        [Category of sets]
        sage: C.all_super_categories()
        [Category of finite sets,
         Category of sets,
         Category of sets with partial maps,
         Category of objects]
        sage: C.example()
        NotImplemented

    TESTS::

        sage: TestSuite(C).run()
        sage: C is Sets().Finite()
        True
    """

    class ParentMethods:

        def is_finite(self):
            """
            Return ``True`` since ``self`` is finite.

            EXAMPLES::

                sage: C = FiniteEnumeratedSets().example()
                sage: C.is_finite()
                True
            """
            return True

    class Subquotients(SubquotientsCategory):

        def extra_super_categories(self):
            r"""
            EXAMPLES::

                sage: FiniteSets().Subquotients().extra_super_categories()
                [Category of finite sets]

            This implements the fact that a subquotient (and therefore
            a quotient or subobject) of a finite set is finite::

                sage: FiniteSets().Subquotients().is_subcategory(FiniteSets())
                True
                sage: FiniteSets().Quotients   ().is_subcategory(FiniteSets())
                True
                sage: FiniteSets().Subobjects  ().is_subcategory(FiniteSets())
                True
            """
            return [FiniteSets()]

    class Algebras(AlgebrasCategory):

        def extra_super_categories(self):
            r"""
            EXAMPLES::

                sage: FiniteSets().Algebras(QQ).extra_super_categories()
                [Category of finite dimensional vector spaces with basis over Rational Field]

            This implements the fact that the algebra of a finite set
            is finite dimensional::

                sage: FiniteMonoids().Algebras(QQ).is_subcategory(AlgebrasWithBasis(QQ).FiniteDimensional())
                True
            """
            from sage.categories.modules_with_basis import ModulesWithBasis
            return [ModulesWithBasis(self.base_ring()).FiniteDimensional()]
