# sage_setup: distribution = sagemath-categories
r"""
Filtered Modules

A *filtered module* over a ring `R` with a totally ordered
indexing set `I` (typically `I = \NN`) is an `R`-module `M` equipped
with a family `(F_i)_{i \in I}` of `R`-submodules satisfying
`F_i \subseteq F_j` for all `i,j \in I` having `i \leq j`, and
`M = \bigcup_{i \in I} F_i`. This family is called a *filtration*
of the given module `M`.

.. TODO::

    Implement a notion for decreasing filtrations: where `F_j \subseteq F_i`
    when `i \leq j`.

.. TODO::

    Implement filtrations for all concrete categories.

.. TODO::

    Implement `\operatorname{gr}` as a functor.
"""
# ****************************************************************************
#  Copyright (C) 2014 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.category_types import Category_over_base_ring
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.covariant_functorial_construction import RegressiveCovariantConstructionCategory


class FilteredModulesCategory(RegressiveCovariantConstructionCategory, Category_over_base_ring):
    def __init__(self, base_category):
        """
        EXAMPLES::

            sage: C = Algebras(QQ).Filtered()
            sage: C
            Category of filtered algebras over Rational Field
            sage: C.base_category()
            Category of algebras over Rational Field
            sage: sorted(C.super_categories(), key=str)
            [Category of algebras over Rational Field,
             Category of filtered vector spaces over Rational Field]

            sage: AlgebrasWithBasis(QQ).Filtered().base_ring()
            Rational Field
            sage: HopfAlgebrasWithBasis(QQ).Filtered().base_ring()
            Rational Field
        """
        super().__init__(base_category, base_category.base_ring())

    _functor_category = "Filtered"

    def _repr_object_names(self):
        """
        EXAMPLES::

            sage: AlgebrasWithBasis(QQ).Filtered()  # indirect doctest
            Category of filtered algebras with basis over Rational Field
        """
        return "filtered {}".format(self.base_category()._repr_object_names())

    def _make_named_class_key(self, name):
        r"""
        Return what the element/parent/... classes depend on.

        .. SEEALSO::

            - :meth:`.CategoryWithParameters._make_named_class_key`

        EXAMPLES::

            sage: Modules(ZZ).Filtered()._make_named_class_key('element_class')
            (<class 'sage.categories.modules.Modules'>,
             Join of Category of Dedekind domains and Category of euclidean domains and Category of noetherian rings and Category of infinite enumerated sets and Category of metric spaces)

        Note that we cannot simply return the base as in
        :meth:`.Category_over_base._make_named_class_key` because of the following
        (see :issue:`39154`)::

            sage: VectorSpacesQQ = VectorSpaces(QQ); VectorSpacesQQ
            Category of vector spaces over Rational Field
            sage: # ModulesQQ = Modules(QQ)  # doesn't work because...
            sage: Modules(QQ) is VectorSpacesQQ
            True
            sage: ModulesQQ = VectorSpacesQQ.super_categories()[0]; ModulesQQ
            Category of modules over Rational Field
            sage: VectorSpacesQQ.Filtered()
            Category of filtered vector spaces over Rational Field
            sage: ModulesQQ.Filtered()
            Category of filtered modules over Rational Field
            sage: VectorSpacesQQ.Filtered()._make_named_class_key('parent_class')
            (<class 'sage.categories.vector_spaces.VectorSpaces'>,
             Join of Category of number fields and Category of quotient fields and Category of metric spaces)
            sage: ModulesQQ.Filtered()._make_named_class_key('parent_class')
            (<class 'sage.categories.modules.Modules'>,
             Join of Category of number fields and Category of quotient fields and Category of metric spaces)
            sage: assert (VectorSpacesQQ.Filtered()._make_named_class_key('parent_class') !=
            ....:         ModulesQQ.Filtered()._make_named_class_key('parent_class'))
            sage: VectorSpacesQQ.Filtered().parent_class
            <class 'sage.categories.vector_spaces.VectorSpaces.Filtered.parent_class'>
            sage: ModulesQQ.Filtered().parent_class
            <class 'sage.categories.filtered_modules.FilteredModules.parent_class'>

        Nevertheless, as explained in :meth:`.Category_over_base._make_named_class_key`,
        ``Modules(QQ).Filtered()`` and ``Modules(QQ.category()).Filtered()`` must have
        the same parent class::

            sage: Modules(QQ).Filtered().parent_class == Modules(QQ.category()).Filtered().parent_class
            True
        """
        return (type(self._base_category).__base__, super()._make_named_class_key(name))


class FilteredModules(FilteredModulesCategory):
    r"""
    The category of filtered modules over a given ring `R`.

    A *filtered module* over a ring `R` with a totally ordered
    indexing set `I` (typically `I = \NN`) is an `R`-module `M` equipped
    with a family `(F_i)_{i \in I}` of `R`-submodules satisfying
    `F_i \subseteq F_j` for all `i,j \in I` having `i \leq j`, and
    `M = \bigcup_{i \in I} F_i`. This family is called a *filtration*
    of the given module `M`.

    EXAMPLES::

        sage: Modules(ZZ).Filtered()
        Category of filtered modules over Integer Ring
        sage: Modules(ZZ).Filtered().super_categories()
        [Category of modules over Integer Ring]

    TESTS::

        sage: TestSuite(Modules(ZZ).Filtered()).run()

    REFERENCES:

    - :wikipedia:`Filtration_(mathematics)`
    """
    def extra_super_categories(self):
        r"""
        Add :class:`VectorSpaces` to the super categories of ``self`` if
        the base ring is a field.

        EXAMPLES::

            sage: Modules(QQ).Filtered().is_subcategory(VectorSpaces(QQ))
            True
            sage: Modules(ZZ).Filtered().extra_super_categories()
            []

        This makes sure that ``Modules(QQ).Filtered()`` returns an
        instance of :class:`FilteredModules` and not a join category of
        an instance of this class and of ``VectorSpaces(QQ)``::

            sage: type(Modules(QQ).Filtered())
            <class 'sage.categories.vector_spaces.VectorSpaces.Filtered_with_category'>

        .. TODO::

            Get rid of this workaround once there is a more systematic
            approach for the alias ``Modules(QQ)`` -> ``VectorSpaces(QQ)``.
            Probably the latter should be a category with axiom, and
            covariant constructions should play well with axioms.
        """
        from sage.categories.modules import Modules
        from sage.categories.fields import Fields
        from sage.categories.category import Category
        base_ring = self.base_ring()
        if base_ring in Fields() or (isinstance(base_ring, Category) and base_ring.is_subcategory(Fields())):
            return [Modules(base_ring)]
        else:
            return []

    class SubcategoryMethods:

        @cached_method
        def Connected(self):
            r"""
            Return the full subcategory of the connected objects of ``self``.

            A filtered `R`-module `M` with filtration
            `(F_0, F_1, F_2, \ldots)` (indexed by `\NN`)
            is said to be *connected* if `F_0` is isomorphic
            to `R`.

            EXAMPLES::

                sage: Modules(ZZ).Filtered().Connected()
                Category of filtered connected modules over Integer Ring
                sage: Coalgebras(QQ).Filtered().Connected()
                Category of filtered connected coalgebras over Rational Field
                sage: AlgebrasWithBasis(QQ).Filtered().Connected()
                Category of filtered connected algebras with basis over Rational Field

            TESTS::

                sage: TestSuite(Modules(ZZ).Filtered().Connected()).run()
                sage: Coalgebras(QQ).Filtered().Connected.__module__
                'sage.categories.filtered_modules'
            """
            return self._with_axiom("Connected")

    class Connected(CategoryWithAxiom_over_base_ring):
        pass
