# sage_setup: distribution = sagemath-categories
r"""
Super modules
"""
#*****************************************************************************
#  Copyright (C) 2015 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_types import Category_over_base_ring
from sage.categories.covariant_functorial_construction import CovariantConstructionCategory

# Note, a commutative algebra is not a commutative super algebra,
#   therefore the following whitelist.
axiom_whitelist = frozenset(["Facade", "Finite", "Infinite",
                             "FiniteDimensional", "Connected", "WithBasis",
                             "FinitelyGeneratedAsLambdaBracketAlgebra",
                             # "Commutative", "Cocommutative",
                             "Supercommutative", "Supercocommutative",
                             "Associative", "Inverse", "Unital", "Division",
                             "AdditiveCommutative", "AdditiveAssociative",
                             "AdditiveInverse", "AdditiveUnital",
                             "NoZeroDivisors", "Distributive"])


class SuperModulesCategory(CovariantConstructionCategory, Category_over_base_ring):
    @classmethod
    def default_super_categories(cls, category, *args):
        """
        Return the default super categories of `F_{Cat}(A,B,...)` for
        `A,B,...` parents in `Cat`.

        INPUT:

        - ``cls`` -- the category class for the functor `F`
        - ``category`` -- a category `Cat`
        - ``*args`` -- further arguments for the functor

        OUTPUT: a join category

        This implements the property that subcategories constructed by
        the set of whitelisted axioms is a subcategory.

        EXAMPLES::

            sage: HopfAlgebras(ZZ).WithBasis().FiniteDimensional().Super()  # indirect doctest
            Category of finite dimensional super Hopf algebras with basis over Integer Ring
        """
        axioms = axiom_whitelist.intersection(category.axioms())
        C = super().default_super_categories(category, *args)
        return C._with_axioms(axioms)

    def __init__(self, base_category):
        """
        EXAMPLES::

            sage: C = Algebras(QQ).Super()
            sage: C
            Category of super algebras over Rational Field
            sage: C.base_category()
            Category of algebras over Rational Field
            sage: sorted(C.super_categories(), key=str)
            [Category of graded algebras over Rational Field,
             Category of super modules over Rational Field]

            sage: AlgebrasWithBasis(QQ).Super().base_ring()
            Rational Field
            sage: HopfAlgebrasWithBasis(QQ).Super().base_ring()
            Rational Field
        """
        super().__init__(base_category, base_category.base_ring())

    _functor_category = "Super"

    def _repr_object_names(self):
        """
        EXAMPLES::

            sage: AlgebrasWithBasis(QQ).Super()  # indirect doctest
            Category of super algebras with basis over Rational Field
        """
        return "super {}".format(self.base_category()._repr_object_names())


class SuperModules(SuperModulesCategory):
    r"""
    The category of super modules.

    An `R`-*super module* (where `R` is a ring) is an `R`-module `M` equipped
    with a decomposition `M = M_0 \oplus M_1` into two `R`-submodules
    `M_0` and `M_1` (called the *even part* and the *odd part* of `M`,
    respectively).

    Thus, an `R`-super module automatically becomes a `\ZZ / 2 \ZZ`-graded
    `R`-module, with `M_0` being the degree-`0` component and `M_1` being the
    degree-`1` component.

    EXAMPLES::

        sage: Modules(ZZ).Super()
        Category of super modules over Integer Ring
        sage: Modules(ZZ).Super().super_categories()
        [Category of graded modules over Integer Ring]

    The category of super modules defines the super structure which
    shall be preserved by morphisms::

        sage: Modules(ZZ).Super().additional_structure()
        Category of super modules over Integer Ring

    TESTS::

        sage: TestSuite(Modules(ZZ).Super()).run()
    """
    def super_categories(self):
        """
        EXAMPLES::

            sage: Modules(ZZ).Super().super_categories()
            [Category of graded modules over Integer Ring]

        Nota bene::

            sage: Modules(QQ).Super()
            Category of super modules over Rational Field
            sage: Modules(QQ).Super().super_categories()
            [Category of graded modules over Rational Field]
        """
        return [self.base_category().Graded()]

    def extra_super_categories(self):
        r"""
        Add :class:`VectorSpaces` to the super categories of ``self`` if
        the base ring is a field.

        EXAMPLES::

            sage: Modules(QQ).Super().extra_super_categories()
            [Category of vector spaces over Rational Field]
            sage: Modules(ZZ).Super().extra_super_categories()
            []

        This makes sure that ``Modules(QQ).Super()`` returns an
        instance of :class:`SuperModules` and not a join category of
        an instance of this class and of ``VectorSpaces(QQ)``::

            sage: type(Modules(QQ).Super())
            <class 'sage.categories.super_modules.SuperModules_with_category'>

        .. TODO::

            Get rid of this workaround once there is a more systematic
            approach for the alias ``Modules(QQ)`` -> ``VectorSpaces(QQ)``.
            Probably the latter should be a category with axiom, and
            covariant constructions should play well with axioms.
        """
        from sage.categories.modules import Modules
        from sage.categories.fields import Fields
        base_ring = self.base_ring()
        if base_ring in Fields():
            return [Modules(base_ring)]
        else:
            return []

    class ParentMethods:
        pass

    class ElementMethods:
        def is_even_odd(self):
            """
            Return ``0`` if ``self`` is an even element or ``1``
            if an odd element.

            .. NOTE::

                The default implementation assumes that the even/odd is
                determined by the parity of :meth:`degree`.

                Overwrite this method if the even/odd behavior is desired
                to be independent.

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: cat = Algebras(QQ).WithBasis().Super()
                sage: C = CombinatorialFreeModule(QQ, Partitions(), category=cat)
                sage: C.degree_on_basis = sum
                sage: C.basis()[2,2,1].is_even_odd()
                1
                sage: C.basis()[2,2].is_even_odd()
                0
            """
            return self.degree() % 2

        def is_even(self):
            """
            Return if ``self`` is an even element.

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: cat = Algebras(QQ).WithBasis().Super()
                sage: C = CombinatorialFreeModule(QQ, Partitions(), category=cat)
                sage: C.degree_on_basis = sum
                sage: C.basis()[2,2,1].is_even()
                False
                sage: C.basis()[2,2].is_even()
                True
            """
            return self.is_even_odd() == 0

        def is_odd(self):
            """
            Return if ``self`` is an odd element.

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: cat = Algebras(QQ).WithBasis().Super()
                sage: C = CombinatorialFreeModule(QQ, Partitions(), category=cat)
                sage: C.degree_on_basis = sum
                sage: C.basis()[2,2,1].is_odd()
                True
                sage: C.basis()[2,2].is_odd()
                False
            """
            return self.is_even_odd() == 1
