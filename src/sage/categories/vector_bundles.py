# sage_setup: distribution = sagemath-categories
# sage.doctest: needs sage.symbolic
r"""
Vector Bundles
"""
#*****************************************************************************
#  Copyright (C) 2019 Michael Jung <micjung at uni-potsdam.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_types import Category_over_base_ring
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.misc.cachefunc import cached_method
from sage.categories.sets_cat import Sets
from sage.categories.fields import Fields


class VectorBundles(Category_over_base_ring):
    r"""
    The category of vector bundles over any base space and base field.

    .. SEEALSO:: :class:`~sage.manifolds.vector_bundle.TopologicalVectorBundle`

    EXAMPLES::

        sage: M = Manifold(2, 'M', structure='top')
        sage: from sage.categories.vector_bundles import VectorBundles
        sage: C = VectorBundles(M, RR); C
        Category of vector bundles over Real Field with 53 bits of precision
         with base space 2-dimensional topological manifold M
        sage: C.super_categories()
        [Category of topological spaces]

    TESTS::

        sage: TestSuite(C).run(skip='_test_category_over_bases')
    """
    def __init__(self, base_space, base_field, name=None):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: from sage.categories.vector_bundles import VectorBundles
            sage: C = VectorBundles(M, RR)
            sage: TestSuite(C).run(skip='_test_category_over_bases')
        """
        if base_field not in Fields().Topological():
            raise ValueError("base field must be a topological field")
        self._base_space = base_space
        Category_over_base_ring.__init__(self, base_field, name)

    @cached_method
    def super_categories(self):
        r"""
        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: from sage.categories.vector_bundles import VectorBundles
            sage: VectorBundles(M, RR).super_categories()
            [Category of topological spaces]
        """
        return [Sets().Topological()]

    def base_space(self):
        r"""
        Return the base space of this category.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='top')
            sage: from sage.categories.vector_bundles import VectorBundles
            sage: VectorBundles(M, RR).base_space()
            2-dimensional topological manifold M
        """
        return self._base_space

    def _repr_object_names(self):
        r"""
        Return the name of the objects of this category.

        .. SEEALSO:: :meth:`Category._repr_object_names`

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: from sage.categories.vector_bundles import VectorBundles
            sage: VectorBundles(M, RR)._repr_object_names()
            'vector bundles over Real Field with 53 bits of precision with base
             space 2-dimensional differentiable manifold M'
        """
        base_space = self._base_space
        return Category_over_base_ring._repr_object_names(self) + \
               " with base space %s" % base_space

    class SubcategoryMethods:
        @cached_method
        def Differentiable(self):
            r"""
            Return the subcategory of the differentiable objects
            of ``self``.

            EXAMPLES::

                sage: M = Manifold(2, 'M')
                sage: from sage.categories.vector_bundles import VectorBundles
                sage: VectorBundles(M, RR).Differentiable()
                Category of differentiable vector bundles over Real Field with
                 53 bits of precision with base space 2-dimensional
                 differentiable manifold M

            TESTS::

                sage: TestSuite(VectorBundles(M, RR).Differentiable()).run()
                sage: VectorBundles(M, RR).Differentiable.__module__
                'sage.categories.vector_bundles'
            """
            return self._with_axiom('Differentiable')

        @cached_method
        def Smooth(self):
            """
            Return the subcategory of the smooth objects of ``self``.

            EXAMPLES::

                sage: M = Manifold(2, 'M')
                sage: from sage.categories.vector_bundles import VectorBundles
                sage: VectorBundles(M, RR).Smooth()
                Category of smooth vector bundles over Real Field with 53 bits
                 of precision with base space 2-dimensional differentiable
                 manifold M

            TESTS::

                sage: TestSuite(VectorBundles(M, RR).Smooth()).run()
                sage: VectorBundles(M, RR).Smooth.__module__
                'sage.categories.vector_bundles'
            """
            return self._with_axiom('Smooth')

    class Differentiable(CategoryWithAxiom_over_base_ring):
        """
        The category of differentiable vector bundles.

        A differentiable vector bundle is a differentiable manifold with
        differentiable surjective projection on a differentiable base space.
        """

    class Smooth(CategoryWithAxiom_over_base_ring):
        """
        The category of smooth vector bundles.

        A smooth vector bundle is a smooth manifold with
        smooth surjective projection on a smooth base space.
        """
