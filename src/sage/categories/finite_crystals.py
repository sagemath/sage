r"""
Finite Crystals
"""
#*****************************************************************************
#  Copyright (C) 2010    Anne Schilling <anne at math.ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.tensor import TensorProductsCategory

class FiniteCrystals(CategoryWithAxiom):
    """
    The category of finite crystals.

    EXAMPLES::

        sage: C = FiniteCrystals()
        sage: C
        Category of finite crystals
        sage: C.super_categories()
        [Category of crystals, Category of finite enumerated sets]
        sage: C.example()
        Highest weight crystal of type A_3 of highest weight omega_1

    TESTS::

        sage: TestSuite(C).run()
        sage: B = FiniteCrystals().example()
        sage: TestSuite(B).run(verbose = True)
        running ._test_an_element() . . . pass
        running ._test_cardinality() . . . pass
        running ._test_category() . . . pass
        running ._test_construction() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
          running ._test_new() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          running ._test_stembridge_local_axioms() . . . pass
          pass
        running ._test_elements_eq_reflexive() . . . pass
        running ._test_elements_eq_symmetric() . . . pass
        running ._test_elements_eq_transitive() . . . pass
        running ._test_elements_neq() . . . pass
        running ._test_enumerated_set_contains() . . . pass
        running ._test_enumerated_set_iter_cardinality() . . . pass
        running ._test_enumerated_set_iter_list() . . . pass
        running ._test_eq() . . . pass
        running ._test_fast_iter() . . . pass
        running ._test_new() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_pickling() . . . pass
        running ._test_some_elements() . . . pass
        running ._test_stembridge_local_axioms() . . . pass
    """

    @cached_method
    def extra_super_categories(self):
        r"""
        EXAMPLES::

            sage: FiniteCrystals().extra_super_categories()
            [Category of finite enumerated sets]
        """
        return [FiniteEnumeratedSets()]

    def example(self, n = 3):
        """
        Returns an example of highest weight crystals, as per
        :meth:`Category.example`.

        EXAMPLES::

            sage: B = FiniteCrystals().example(); B
            Highest weight crystal of type A_3 of highest weight omega_1
        """
        from sage.categories.crystals import Crystals
        return Crystals().example(n)

    class TensorProducts(TensorProductsCategory):
        """
        The category of finite crystals constructed by tensor
        product of finite crystals.
        """
        @cached_method
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: FiniteCrystals().TensorProducts().extra_super_categories()
                [Category of finite crystals]
            """
            return [self.base_category()]

