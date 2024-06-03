# sage_setup: distribution = sagemath-categories
r"""
Supercommutative Algebras
"""
#*****************************************************************************
#  Copyright (C) 2019 Travis Scrimshaw <tcscrims at gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.super_algebras import SuperAlgebras
from sage.categories.signed_tensor import SignedTensorProductsCategory
from sage.misc.cachefunc import cached_method


class SupercommutativeAlgebras(CategoryWithAxiom_over_base_ring):
    r"""
    The category of supercommutative algebras.

    An `R`-*supercommutative algebra* is an `R`-super algebra
    `A = A_0 \oplus A_1` endowed with an `R`-super algebra structure
    satisfying:

    .. MATH::

        x_0 x'_0 = x'_0 x_0, \qquad
        x_1 x'_1 = -x'_1 x_1, \qquad
        x_0 x_1 = x_1 x_0,

    for all `x_0, x'_0 \in A_0` and `x_1, x'_1 \in A_1`.

    EXAMPLES::

        sage: Algebras(ZZ).Supercommutative()
        Category of supercommutative algebras over Integer Ring

    TESTS::

        sage: TestSuite(Algebras(ZZ).Supercommutative()).run()
    """
    _base_category_class_and_axiom = (SuperAlgebras, "Supercommutative")

    class SignedTensorProducts(SignedTensorProductsCategory):
        @cached_method
        def extra_super_categories(self):
            """
            Return the extra super categories of ``self``.

            A signed tensor product of supercommutative algebras is a
            supercommutative algebra.

            EXAMPLES::

                sage: C = Algebras(ZZ).Supercommutative().SignedTensorProducts()
                sage: C.extra_super_categories()
                [Category of supercommutative algebras over Integer Ring]
            """
            return [self.base_category()]

    class WithBasis(CategoryWithAxiom_over_base_ring):
        class ParentMethods:
            def _test_supercommutativity(self, **options):
                r"""
                Test supercommutativity for (not necessarily all) elements
                of this supercommutative algebra.

                INPUT:

                - ``options`` -- any keyword arguments accepted by :meth:`_tester`

                EXAMPLES:

                By default, this method tests only the elements returned by
                ``self.some_elements()``::

                    sage: E.<x,y,z> = ExteriorAlgebra(QQ)                               # needs sage.modules
                    sage: E._test_supercommutativity()                                  # needs sage.modules

                However, the elements tested can be customized with the
                ``elements`` keyword argument, but the elements must be
                homogeneous::

                    sage: E._test_supercommutativity(elements=[x+y, x*y-3*y*z, x*y*z])  # needs sage.modules
                    sage: E._test_supercommutativity(elements=[x+x*y])                  # needs sage.modules
                    Traceback (most recent call last):
                    ...
                    ValueError: element is not homogeneous

                See the documentation for :class:`TestSuite` for more information.
                """
                elements = options.pop("elements", self.basis())
                tester = self._tester(**options)
                from sage.misc.misc import some_tuples
                for x,y in some_tuples(elements, 2, tester._max_runs):
                    tester.assertEqual((x * y),
                                       (-1)**(x.is_even_odd() * y.is_even_odd())
                                       * (y * x))
