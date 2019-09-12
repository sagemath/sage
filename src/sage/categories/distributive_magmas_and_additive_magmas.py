r"""
Distributive Magmas and Additive Magmas
"""
# ****************************************************************************
#  Copyright (C) 2010 Nicolas Borie <nicolas.borie@math.u-psud.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.misc.lazy_import import LazyImport
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.cartesian_product import CartesianProductsCategory


class DistributiveMagmasAndAdditiveMagmas(CategoryWithAxiom):
    """
    The category of sets `(S,+,*)` with `*` distributing on `+`.

    This is similar to a ring, but `+` and `*` are only required to be
    (additive) magmas.

    EXAMPLES::

        sage: from sage.categories.distributive_magmas_and_additive_magmas import DistributiveMagmasAndAdditiveMagmas
        sage: C = DistributiveMagmasAndAdditiveMagmas(); C
        Category of distributive magmas and additive magmas
        sage: C.super_categories()
        [Category of magmas and additive magmas]

    TESTS::

        sage: from sage.categories.magmas_and_additive_magmas import MagmasAndAdditiveMagmas
        sage: C is MagmasAndAdditiveMagmas().Distributive()
        True
        sage: C is (Magmas() & AdditiveMagmas()).Distributive()
        True
        sage: TestSuite(C).run()
    """

    class AdditiveAssociative(CategoryWithAxiom):
        class AdditiveCommutative(CategoryWithAxiom):
            class AdditiveUnital(CategoryWithAxiom):
                class Associative(CategoryWithAxiom):
                    AdditiveInverse = LazyImport('sage.categories.rngs', 'Rngs', at_startup=True)
                    Unital          = LazyImport('sage.categories.semirings', 'Semirings', at_startup=True)

    class ParentMethods:

        def _test_distributivity(self, **options):
            r"""
            Test the distributivity of `*` on `+` on (not necessarily
            all) elements of this set.

            INPUT:

            - ``options`` -- any keyword arguments accepted by :meth:`_tester`

            EXAMPLES:

            By default, this method runs the tests only on the
            elements returned by ``self.some_elements()``::

                sage: NN.some_elements()
                [0, 1, 3, 42]
                sage: NN._test_distributivity()

            However, the elements tested can be customized with the
            ``elements`` keyword argument::

                sage: CC._test_distributivity(elements=[CC(0),CC(1),CC(3),CC(I)])

            See the documentation for :class:`TestSuite` for more information.
            """
            tester = self._tester(**options)
            tester.some_elements()
            from sage.misc.misc import some_tuples
            for x, y, z in some_tuples(tester.some_elements(), 3, tester._max_runs):
                # left distributivity
                tester.assertEqual(x * (y + z), (x * y) + (x * z))
                # right distributivity
                tester.assertEqual((x + y) * z, (x * z) + (y * z))

    class CartesianProducts(CartesianProductsCategory):
        def extra_super_categories(self):
            """
            Implement the fact that a Cartesian product of magmas distributing
            over additive magmas is a magma distributing over an
            additive magma.

            EXAMPLES::

                sage: C = (Magmas() & AdditiveMagmas()).Distributive().CartesianProducts()
                sage: C.extra_super_categories()
                [Category of distributive magmas and additive magmas]
                sage: C.axioms()
                frozenset({'Distributive'})
            """
            return [DistributiveMagmasAndAdditiveMagmas()]
