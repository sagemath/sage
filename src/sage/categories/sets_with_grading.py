r"""
Sets With a Grading
"""
# ****************************************************************************
#  Copyright (C) 2010-2012 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.abstract_method import abstract_method
from .category_types import Category
from sage.categories.sets_cat import Sets
from sage.categories.enumerated_sets import EnumeratedSets


class SetsWithGrading(Category):
    r"""
    The category of sets with a grading.

    A *set with a grading* is a set `S` equipped with a
    grading by some other set `I` (by default the set `\NN` of the
    non-negative integers):

    .. MATH::

         S = \biguplus_{i\in I} S_i

    where the *graded components* `S_i` are (usually finite)
    sets. The *grading* function maps each element `s` of
    `S` to its *grade* `i`, so that `s\in S_i`.

    From implementation point of view, if the graded set is enumerated then
    each graded component should be enumerated (there is a check in the method
    :meth:`~SetsWithGrading.ParentMethods._test_graded_components`). The
    contrary needs not be true.

    To implement this category, a parent must either implement
    :meth:`~SetsWithGrading.ParentMethods.graded_component()` or
    :meth:`~SetsWithGrading.ParentMethods.subset()`. If only
    :meth:`~SetsWithGrading.ParentMethods.subset()` is implemented, the first
    argument must be the grading for compatibility with
    :meth:`~SetsWithGrading.ParentMethods.graded_component()`. Additionally
    either the parent must implement
    :meth:`~SetsWithGrading.ParentMethods.grading()` or its elements must
    implement a method ``grade()``. See the example
    :class:`sage.categories.examples.sets_with_grading.NonNegativeIntegers`.

    Finally, if the graded set is enumerated (see
    :class:`~sage.categories.enumerated_sets.EnumeratedSets`) then each graded
    component should be enumerated. The contrary needs not be true.

    EXAMPLES:

    A typical example of a set with a grading is the set of non-negative
    integers graded by themselves::

        sage: N = SetsWithGrading().example(); N
        Non negative integers
        sage: N.category()
        Category of facade infinite sets with grading
        sage: N.grading_set()
        Non negative integers

    The *grading function* is given by ``N.grading``::

        sage: N.grading(4)
        4

    The graded component `N_i` is the set with one element `i`::

        sage: N.graded_component(grade=5)
        {5}
        sage: N.graded_component(grade=42)
        {42}

    Here are some information about this category::

        sage: SetsWithGrading()
        Category of sets with grading
        sage: SetsWithGrading().super_categories()
        [Category of sets]
        sage: SetsWithGrading().all_super_categories()
        [Category of sets with grading,
         Category of sets,
         Category of sets with partial maps,
         Category of objects]

    .. TODO::

        - This should be moved to ``Sets().WithGrading()``.
        - Should the grading set be a parameter for this category?
        - Does the enumeration need to be compatible with the grading? Be
          careful that the fact that graded components are allowed to be finite
          or infinite make the answer complicated.

    TESTS::

        sage: C = SetsWithGrading()
        sage: TestSuite(C).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: SetsWithGrading().super_categories()
            [Category of sets]
        """
        return [Sets()]

    class ParentMethods:

        def _test_graded_components(self, **options):
            r"""
            Test that some graded components of ``self`` are parent with
            initialized category and that the parent has a properly implemented
            ``grading()`` method.

            EXAMPLES::

                sage: SetsWithGrading().example()._test_graded_components()
            """
            tester = self._tester(**options)
            for grade in self.grading_set().some_elements():
                G = self.graded_component(grade)
                if self in EnumeratedSets():
                    tester.assertIn(G, EnumeratedSets())
                else:
                    tester.assertIn(G, Sets())
                for elt in G.some_elements():
                    tester.assertEqual(self.grading(elt), grade)

        def grading_set(self):
            """
            Return the set ``self`` is graded by. By default, this is
            the set of non-negative integers.

            EXAMPLES::

                sage: SetsWithGrading().example().grading_set()
                Non negative integers
            """
            from sage.sets.non_negative_integers import NonNegativeIntegers
            return NonNegativeIntegers()

        # TODO:
        #  - Should this method be in EnumeratedSets? With a default implementation
        #    a la ``filter``?
        #  - Do we want to enforce implementing subset rather than graded_component?
        @abstract_method(optional=True)
        def subset(self, *args, **options):
            """
            Return the subset of ``self`` described by the given parameters.

            .. SEEALSO::

                -:meth:`graded_component()`

            EXAMPLES::

                sage: W = WeightedIntegerVectors([3,2,1]); W
                Integer vectors weighted by [3, 2, 1]
                sage: W.subset(4)
                Integer vectors of 4 weighted by [3, 2, 1]
            """

        def graded_component(self, grade):
            """
            Return the graded component of ``self`` with grade ``grade``.

            The default implementation just calls the method :meth:`subset()`
            with the first argument ``grade``.

            EXAMPLES::

                sage: N = SetsWithGrading().example(); N
                Non negative integers
                sage: N.graded_component(3)
                {3}
            """
            return self.subset(grade)

        def grading(self, elt):
            """
            Return the grading of the element ``elt`` of ``self``.

            This default implementation calls ``elt.grade()``.

            EXAMPLES::

                sage: N = SetsWithGrading().example(); N
                Non negative integers
                sage: N.grading(4)
                4
            """
            return elt.grade()

        def generating_series(self):
            """
            Default implementation for generating series.

            OUTPUT:

            A series, indexed by the grading set.

            EXAMPLES::

                sage: N = SetsWithGrading().example(); N
                Non negative integers
                sage: N.generating_series()
                1/(-z + 1)
            """
            from sage.combinat.species.series import LazyPowerSeriesRing
            from sage.rings.integer_ring import ZZ
            R = LazyPowerSeriesRing(ZZ)
            R(self.graded_component(grade).cardinality() for grade in self.grading_set())

        # TODO:
        #   * asymptotic behavior: we need an object for asymptotic behavior and
        #   a default name for the method that should be here. Such method will
        #   have two goals (and perhaps need two implementations): give a
        #   theorem on asymptotic and be a tool to determine a strategy for
        #   algorithms.

