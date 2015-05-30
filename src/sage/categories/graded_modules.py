r"""
Graded modules
"""
#*****************************************************************************
#  Copyright (C) 2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2013 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_class_attribute
from sage.categories.category_types import Category_over_base_ring
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.covariant_functorial_construction import RegressiveCovariantConstructionCategory

class GradedModulesCategory(RegressiveCovariantConstructionCategory, Category_over_base_ring):
    def __init__(self, base_category):
        """
        EXAMPLES::

            sage: C = GradedAlgebras(QQ)
            sage: C
            Category of graded algebras over Rational Field
            sage: C.base_category()
            Category of algebras over Rational Field
            sage: sorted(C.super_categories(), key=str)
            [Category of algebras over Rational Field,
             Category of graded modules over Rational Field]

            sage: AlgebrasWithBasis(QQ).Graded().base_ring()
            Rational Field
            sage: GradedHopfAlgebrasWithBasis(QQ).base_ring()
            Rational Field

        TESTS::

            sage: GradedModules(ZZ)
            Category of graded modules over Integer Ring
            sage: Modules(ZZ).Graded()
            Category of graded modules over Integer Ring
            sage: GradedModules(ZZ) is Modules(ZZ).Graded()
            True
        """
        super(GradedModulesCategory, self).__init__(base_category, base_category.base_ring())

    _functor_category = "Graded"

    def _repr_object_names(self):
        """
        EXAMPLES::

            sage: AlgebrasWithBasis(QQ).Graded()  # indirect doctest
            Category of graded algebras with basis over Rational Field
        """
        return "graded {}".format(self.base_category()._repr_object_names())

class GradedModules(GradedModulesCategory):
    """
    The category of graded modules.

    EXAMPLES::

        sage: GradedModules(ZZ)
        Category of graded modules over Integer Ring
        sage: GradedModules(ZZ).super_categories()
        [Category of modules over Integer Ring]

    The category of graded modules defines the graded structure which
    shall be preserved by morphisms::

        sage: Modules(ZZ).Graded().additional_structure()
        Category of graded modules over Integer Ring

    TESTS::

        sage: TestSuite(GradedModules(ZZ)).run()
    """

    def extra_super_categories(self):
        r"""
        Adds :class:`VectorSpaces` to the super categories of ``self`` if
        the base ring is a field.

        EXAMPLES::

            sage: Modules(QQ).Graded().extra_super_categories()
            [Category of vector spaces over Rational Field]
            sage: Modules(ZZ).Graded().extra_super_categories()
            []

        This makes sure that ``Modules(QQ).Graded()`` returns an
        instance of :class:`GradedModules` and not a join category of
        an instance of this class and of ``VectorSpaces(QQ)``::

            sage: type(Modules(QQ).Graded())
            <class 'sage.categories.graded_modules.GradedModules_with_category'>

        .. TODO::

            Get rid of this workaround once there is a more systematic
            approach for the alias ``Modules(QQ)`` -> ``VectorSpaces(QQ)``.
            Probably the later should be a category with axiom, and
            covariant constructions should play well with axioms.
        """
        from sage.categories.modules import Modules
        from sage.categories.fields import Fields
        base_ring = self.base_ring()
        if base_ring in Fields:
            return [Modules(base_ring)]
        else:
            return []

    class SubcategoryMethods:

        @cached_method
        def Connected(self):
            r"""
            Return the full subcategory of the connected objects of ``self``.

            EXAMPLES::

                sage: Modules(ZZ).Graded().Connected()
                Category of graded connected modules over Integer Ring
                sage: Coalgebras(QQ).Graded().Connected()
                Join of Category of graded connected modules over Rational Field
                    and Category of coalgebras over Rational Field
                sage: GradedAlgebrasWithBasis(QQ).Connected()
                Category of graded connected algebras with basis over Rational Field

            TESTS::

                sage: TestSuite(Modules(ZZ).Graded().Connected()).run()
                sage: Coalgebras(QQ).Graded().Connected.__module__
                'sage.categories.graded_modules'
            """
            return self._with_axiom("Connected")

    class Connected(CategoryWithAxiom_over_base_ring):
        pass

    class ParentMethods:
        pass

    class ElementMethods:
        pass
