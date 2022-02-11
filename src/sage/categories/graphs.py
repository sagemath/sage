"""
Graphs
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
from sage.categories.simplicial_complexes import SimplicialComplexes
from sage.categories.sets_cat import Sets

class Graphs(Category_singleton):
    r"""
    The category of graphs.

    EXAMPLES::

        sage: from sage.categories.graphs import Graphs
        sage: C = Graphs(); C
        Category of graphs

    TESTS::

        sage: TestSuite(C).run()
    """
    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: from sage.categories.graphs import Graphs
            sage: Graphs().super_categories()
            [Category of simplicial complexes]
        """
        return [SimplicialComplexes()]

    class ParentMethods:
        @abstract_method
        def vertices(self):
            """
            Return the vertices of ``self``.

            EXAMPLES::

                sage: from sage.categories.graphs import Graphs
                sage: C = Graphs().example()
                sage: C.vertices()
                [0, 1, 2, 3, 4]
            """

        @abstract_method
        def edges(self):
            """
            Return the edges of ``self``.

            EXAMPLES::

                sage: from sage.categories.graphs import Graphs
                sage: C = Graphs().example()
                sage: C.edges()
                [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]
            """

        def dimension(self):
            """
            Return the dimension of ``self`` as a CW complex.

            EXAMPLES::

                sage: from sage.categories.graphs import Graphs
                sage: C = Graphs().example()
                sage: C.dimension()
                1
            """
            if self.edges():
                return 1
            return 0

        def facets(self):
            """
            Return the facets of ``self``.

            EXAMPLES::

                sage: from sage.categories.graphs import Graphs
                sage: C = Graphs().example()
                sage: C.facets()
                [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]
            """
            return self.edges()

        def faces(self):
            """
            Return the faces of ``self``.

            EXAMPLES::

                sage: from sage.categories.graphs import Graphs
                sage: C = Graphs().example()
                sage: sorted(C.faces(), key=lambda x: (x.dimension(), x.value))
                [0, 1, 2, 3, 4, (0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]
            """
            return set(self.edges()).union(self.vertices())

    class Connected(CategoryWithAxiom):
        """
        The category of connected graphs.

        EXAMPLES::

            sage: from sage.categories.graphs import Graphs
            sage: C = Graphs().Connected()
            sage: TestSuite(C).run()
        """
        def extra_super_categories(self):
            """
            Return the extra super categories of ``self``.

            A connected graph is also a metric space.

            EXAMPLES::

                sage: from sage.categories.graphs import Graphs
                sage: Graphs().Connected().super_categories() # indirect doctest
                [Category of connected topological spaces,
                 Category of connected simplicial complexes,
                 Category of graphs,
                 Category of metric spaces]
            """
            return [Sets().Metric()]

