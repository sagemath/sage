"""
Simplicial Complexes
"""
#*****************************************************************************
#  Copyright (C) 2015 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.categories.category_singleton import Category_singleton
from sage.categories.category_with_axiom import CategoryWithAxiom
#from sage.categories.cw_complexes import CWComplexes
from sage.categories.sets_cat import Sets

class SimplicialComplexes(Category_singleton):
    r"""
    The category of abstract simplicial complexes.

    An abstract simplicial complex `A` is a collection of sets `X`
    such that:

    - `\emptyset \in A`,
    - if `X \subset Y \in A`, then `X \in A`.

    .. TODO::

        Implement the category of simplicial complexes considered
        as :class:`CW complexes <sage.categories.cw_complexes.CWComplexes>`
        and rename this to the category of ``AbstractSimplicialComplexes``
        with appropriate functors.

    EXAMPLES::

        sage: from sage.categories.simplicial_complexes import SimplicialComplexes
        sage: C = SimplicialComplexes(); C
        Category of simplicial complexes

    TESTS::

        sage: TestSuite(C).run()
    """
    @cached_method
    def super_categories(self):
        """
        Return the super categories of ``self``.

        EXAMPLES::

            sage: from sage.categories.simplicial_complexes import SimplicialComplexes
            sage: SimplicialComplexes().super_categories()
            [Category of sets]
        """
        return [Sets()]

    class Finite(CategoryWithAxiom):
        """
        Category of finite simplicial complexes.
        """
        class ParentMethods:
            @cached_method
            def dimension(self):
                """
                Return the dimension of ``self``.

                EXAMPLES::

                    sage: S = SimplicialComplex([[1,3,4], [1,2],[2,5],[4,5]])
                    sage: S.dimension()
                    2
                """
                return max(c.dimension() for c in self.facets())

    class ParentMethods:
        @abstract_method
        def facets(self):
            """
            Return the facets of ``self``.

            EXAMPLES::

                sage: S = SimplicialComplex([[1,3,4], [1,2],[2,5],[4,5]])
                sage: sorted(S.facets())
                [(1, 2), (1, 3, 4), (2, 5), (4, 5)]
            """

        @abstract_method
        def faces(self):
            """
            Return the faces of ``self``.

            EXAMPLES::

                sage: S = SimplicialComplex([[1,3,4], [1,2],[2,5],[4,5]])
                sage: S.faces()
                {-1: {()},
                 0: {(1,), (2,), (3,), (4,), (5,)},
                 1: {(1, 2), (1, 3), (1, 4), (2, 5), (3, 4), (4, 5)},
                 2: {(1, 3, 4)}}
            """

    class SubcategoryMethods:
        @cached_method
        def Connected(self):
            """
            Return the full subcategory of the connected objects of ``self``.

            EXAMPLES::

                sage: from sage.categories.simplicial_complexes import SimplicialComplexes
                sage: SimplicialComplexes().Connected()
                Category of connected simplicial complexes

            TESTS::

                sage: SimplicialComplexes().Connected.__module__
                'sage.categories.simplicial_complexes'
            """
            return self._with_axiom('Connected')

    class Connected(CategoryWithAxiom):
        """
        The category of connected simplicial complexes.

        EXAMPLES::

            sage: from sage.categories.simplicial_complexes import SimplicialComplexes
            sage: C = SimplicialComplexes().Connected()
            sage: TestSuite(C).run()
        """

