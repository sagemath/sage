# sage_setup: distribution = sagemath-categories
r"""
Manifolds
"""
#*****************************************************************************
#  Copyright (C) 2015 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.categories.category_types import Category_over_base_ring
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.sets_cat import Sets
from sage.categories.fields import Fields


class Manifolds(Category_over_base_ring):
    r"""
    The category of manifolds over any topological field.

    Let `k` be a topological field. A `d`-dimensional `k`-*manifold* `M`
    is a second countable Hausdorff space such that the neighborhood of
    any point `x \in M` is homeomorphic to `k^d`.

    EXAMPLES::

        sage: # needs sage.rings.real_mpfr
        sage: from sage.categories.manifolds import Manifolds
        sage: C = Manifolds(RR); C
        Category of manifolds over Real Field with 53 bits of precision
        sage: C.super_categories()
        [Category of topological spaces]

    TESTS::

        sage: TestSuite(C).run(skip='_test_category_over_bases')                        # needs sage.rings.real_mpfr
    """
    def __init__(self, base, name=None):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.categories.manifolds import Manifolds
            sage: C = Manifolds(RR)
            sage: TestSuite(C).run(skip='_test_category_over_bases')
        """
        if base not in Fields().Topological():
            raise ValueError("base must be a topological field")
        Category_over_base_ring.__init__(self, base, name)

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: from sage.categories.manifolds import Manifolds
            sage: Manifolds(RR).super_categories()
            [Category of topological spaces]
        """
        return [Sets().Topological()]

    def additional_structure(self):
        r"""
        Return ``None``.

        Indeed, the category of manifolds defines no new
        structure: a morphism of topological spaces between
        manifolds is a manifold morphism.

        .. SEEALSO:: :meth:`Category.additional_structure`

        EXAMPLES::

            sage: from sage.categories.manifolds import Manifolds
            sage: Manifolds(RR).additional_structure()
        """
        return None

    class ParentMethods:
        @abstract_method
        def dimension(self):
            """
            Return the dimension of ``self``.

            EXAMPLES::

                sage: from sage.categories.manifolds import Manifolds
                sage: M = Manifolds(RR).example()
                sage: M.dimension()
                3
            """

    class SubcategoryMethods:
        @cached_method
        def Connected(self):
            """
            Return the full subcategory of the connected objects of ``self``.

            EXAMPLES::

                sage: from sage.categories.manifolds import Manifolds
                sage: Manifolds(RR).Connected()                                         # needs sage.rings.real_mpfr
                Category of connected manifolds
                 over Real Field with 53 bits of precision

            TESTS::

                sage: Manifolds(RR).Connected.__module__
                'sage.categories.manifolds'
            """
            return self._with_axiom('Connected')

        @cached_method
        def FiniteDimensional(self):
            """
            Return the full subcategory of the finite dimensional
            objects of ``self``.

            EXAMPLES::

                sage: from sage.categories.manifolds import Manifolds
                sage: C = Manifolds(RR).Connected().FiniteDimensional(); C              # needs sage.rings.real_mpfr
                Category of finite dimensional connected manifolds
                 over Real Field with 53 bits of precision

            TESTS::

                sage: from sage.categories.manifolds import Manifolds
                sage: Manifolds(RR).Connected().FiniteDimensional.__module__
                'sage.categories.manifolds'
            """
            return self._with_axiom('FiniteDimensional')

        @cached_method
        def Differentiable(self):
            """
            Return the subcategory of the differentiable objects
            of ``self``.

            EXAMPLES::

                sage: from sage.categories.manifolds import Manifolds
                sage: Manifolds(RR).Differentiable()                                    # needs sage.rings.real_mpfr
                Category of differentiable manifolds
                 over Real Field with 53 bits of precision

            TESTS::

                sage: TestSuite(Manifolds(RR).Differentiable()).run()
                sage: Manifolds(RR).Differentiable.__module__
                'sage.categories.manifolds'
            """
            return self._with_axiom('Differentiable')

        @cached_method
        def Smooth(self):
            """
            Return the subcategory of the smooth objects of ``self``.

            EXAMPLES::

                sage: from sage.categories.manifolds import Manifolds
                sage: Manifolds(RR).Smooth()                                            # needs sage.rings.real_mpfr
                Category of smooth manifolds
                 over Real Field with 53 bits of precision

            TESTS::

                sage: TestSuite(Manifolds(RR).Smooth()).run()
                sage: Manifolds(RR).Smooth.__module__
                'sage.categories.manifolds'
            """
            return self._with_axiom('Smooth')

        @cached_method
        def Analytic(self):
            """
            Return the subcategory of the analytic objects of ``self``.

            EXAMPLES::

                sage: from sage.categories.manifolds import Manifolds
                sage: Manifolds(RR).Analytic()                                          # needs sage.rings.real_mpfr
                Category of analytic manifolds
                 over Real Field with 53 bits of precision

            TESTS::

                sage: TestSuite(Manifolds(RR).Analytic()).run()
                sage: Manifolds(RR).Analytic.__module__
                'sage.categories.manifolds'
            """
            return self._with_axiom('Analytic')

        @cached_method
        def AlmostComplex(self):
            """
            Return the subcategory of the almost complex objects
            of ``self``.

            EXAMPLES::

                sage: from sage.categories.manifolds import Manifolds
                sage: Manifolds(RR).AlmostComplex()                                     # needs sage.rings.real_mpfr
                Category of almost complex manifolds
                 over Real Field with 53 bits of precision

            TESTS::

                sage: TestSuite(Manifolds(RR).AlmostComplex()).run()
                sage: Manifolds(RR).AlmostComplex.__module__
                'sage.categories.manifolds'
            """
            return self._with_axiom('AlmostComplex')

        @cached_method
        def Complex(self):
            r"""
            Return the subcategory of manifolds over `\CC` of ``self``.

            EXAMPLES::

                sage: from sage.categories.manifolds import Manifolds
                sage: Manifolds(CC).Complex()                                           # needs sage.rings.real_mpfr
                Category of complex manifolds over
                 Complex Field with 53 bits of precision

            TESTS::

                sage: TestSuite(Manifolds(CC).Complex()).run()                          # needs sage.rings.real_mpfr
                sage: Manifolds(CC).Complex.__module__                                  # needs sage.rings.real_mpfr
                'sage.categories.manifolds'
            """
            return ComplexManifolds(self.base())._with_axioms(self.axioms())

    class Differentiable(CategoryWithAxiom_over_base_ring):
        """
        The category of differentiable manifolds.

        A differentiable manifold is a manifold with a differentiable atlas.
        """

    class Smooth(CategoryWithAxiom_over_base_ring):
        """
        The category of smooth manifolds.

        A smooth manifold is a manifold with a smooth atlas.
        """
        def extra_super_categories(self):
            """
            Return the extra super categories of ``self``.

            A smooth manifold is differentiable.

            EXAMPLES::

                sage: from sage.categories.manifolds import Manifolds
                sage: Manifolds(RR).Smooth().super_categories()  # indirect doctest     # needs sage.rings.real_mpfr
                [Category of differentiable manifolds
                 over Real Field with 53 bits of precision]
            """
            return [Manifolds(self.base()).Differentiable()]

    class Analytic(CategoryWithAxiom_over_base_ring):
        r"""
        The category of complex manifolds.

        An analytic manifold is a manifold with an analytic atlas.
        """
        def extra_super_categories(self):
            """
            Return the extra super categories of ``self``.

            An analytic manifold is smooth.

            EXAMPLES::

                sage: from sage.categories.manifolds import Manifolds
                sage: Manifolds(RR).Analytic().super_categories()  # indirect doctest   # needs sage.rings.real_mpfr
                [Category of smooth manifolds
                 over Real Field with 53 bits of precision]
            """
            return [Manifolds(self.base()).Smooth()]

    class AlmostComplex(CategoryWithAxiom_over_base_ring):
        r"""
        The category of almost complex manifolds.

        An *almost complex manifold* `M` is a manifold with a smooth tensor
        field `J` of rank `(1, 1)` such that `J^2 = -1` when regarded as a
        vector bundle isomorphism `J : TM \to TM` on the tangent bundle.
        The tensor field `J` is called the *almost complex structure* of `M`.
        """
        def extra_super_categories(self):
            """
            Return the extra super categories of ``self``.

            An almost complex manifold is smooth.

            EXAMPLES::

                sage: from sage.categories.manifolds import Manifolds
                sage: Manifolds(RR).AlmostComplex().super_categories()  # indirect doctest          # needs sage.rings.real_mpfr
                [Category of smooth manifolds
                 over Real Field with 53 bits of precision]
            """
            return [Manifolds(self.base()).Smooth()]

    class FiniteDimensional(CategoryWithAxiom_over_base_ring):
        """
        Category of finite dimensional manifolds.

        EXAMPLES::

            sage: from sage.categories.manifolds import Manifolds
            sage: C = Manifolds(RR).FiniteDimensional()
            sage: TestSuite(C).run(skip='_test_category_over_bases')
        """

    class Connected(CategoryWithAxiom_over_base_ring):
        """
        The category of connected manifolds.

        EXAMPLES::

            sage: from sage.categories.manifolds import Manifolds
            sage: C = Manifolds(RR).Connected()
            sage: TestSuite(C).run(skip='_test_category_over_bases')
        """


class ComplexManifolds(Category_over_base_ring):
    r"""
    The category of complex manifolds.

    A `d`-dimensional complex manifold is a manifold whose underlying
    vector space is `\CC^d` and has a holomorphic atlas.
    """
    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: from sage.categories.manifolds import Manifolds
            sage: Manifolds(RR).super_categories()
            [Category of topological spaces]
        """
        return [Manifolds(self.base()).Analytic()]
