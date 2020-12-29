r"""
Regular Supercrystals
"""

# ****************************************************************************
#       Copyright (C) 2017 Franco Saliola <saliola@gmail.com>
#                     2017 Anne Schilling <anne at math.ucdavis.edu>
#                     2017 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.category_singleton import Category_singleton
from sage.categories.supercrystals import SuperCrystals
from sage.categories.tensor import TensorProductsCategory


class RegularSuperCrystals(Category_singleton):
    r"""
    The category of crystals for super Lie algebras.

    EXAMPLES::

        sage: from sage.categories.regular_supercrystals import RegularSuperCrystals
        sage: C = RegularSuperCrystals()
        sage: C
        Category of regular super crystals
        sage: C.super_categories()
        [Category of finite super crystals]

    Parents in this category should implement the following methods:

    - either an attribute ``_cartan_type`` or a method ``cartan_type``

    - ``module_generators``: a list (or container) of distinct elements
      that generate the crystal using `f_i` and `e_i`

    Furthermore, their elements ``x`` should implement the following
    methods:

    - ``x.e(i)`` (returning `e_i(x)`)

    - ``x.f(i)`` (returning `f_i(x)`)

    - ``x.weight()`` (returning `\operatorname{wt}(x)`)

    EXAMPLES::

        sage: from sage.misc.abstract_method import abstract_methods_of_class
        sage: from sage.categories.regular_supercrystals import RegularSuperCrystals
        sage: abstract_methods_of_class(RegularSuperCrystals().element_class)
        {'optional': [], 'required': ['e', 'f', 'weight']}

    TESTS::

        sage: from sage.categories.regular_supercrystals import RegularSuperCrystals
        sage: C = RegularSuperCrystals()
        sage: TestSuite(C).run()
        sage: B = crystals.Letters(['A',[1,1]]); B
        The crystal of letters for type ['A', [1, 1]]
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
          pass
        running ._test_elements_eq_reflexive() . . . pass
        running ._test_elements_eq_symmetric() . . . pass
        running ._test_elements_eq_transitive() . . . pass
        running ._test_elements_neq() . . . pass
        running ._test_enumerated_set_contains() . . . pass
        running ._test_enumerated_set_iter_cardinality() . . . pass
        running ._test_enumerated_set_iter_list() . . . pass
        running ._test_eq() . . . pass
        running ._test_new() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_pickling() . . . pass
        running ._test_some_elements() . . . pass
    """
    def super_categories(self):
        r"""
        EXAMPLES::

            sage: from sage.categories.regular_supercrystals import RegularSuperCrystals
            sage: C = RegularSuperCrystals()
            sage: C.super_categories()
            [Category of finite super crystals]
        """
        return [SuperCrystals().Finite()]

    class ElementMethods:
        def epsilon(self, i):
            r"""
            Return `\varepsilon_i` of ``self``.

            EXAMPLES::

                sage: C = crystals.Tableaux(['A',[1,2]], shape = [2,1])
                sage: c = C.an_element(); c
                [[-2, -2], [-1]]
                sage: c.epsilon(2)
                0
                sage: c.epsilon(0)
                0
                sage: c.epsilon(-1)
                0
            """
            string_length = 0
            x = self
            while True:
                x = x.e(i)
                if x is None:
                    return string_length
                else:
                    string_length += 1

        def phi(self, i):
            r"""
            Return `\varphi_i` of ``self``.

            EXAMPLES::

                sage: C = crystals.Tableaux(['A',[1,2]], shape = [2,1])
                sage: c = C.an_element(); c
                [[-2, -2], [-1]]
                sage: c.phi(1)
                0
                sage: c.phi(2)
                0
                sage: c.phi(0)
                1
            """
            string_length = 0
            x = self
            while True:
                x = x.f(i)
                if x is None:
                    return string_length
                else:
                    string_length += 1

    class TensorProducts(TensorProductsCategory):
        """
        The category of regular crystals constructed by tensor
        product of regular crystals.
        """
        @cached_method
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: from sage.categories.regular_supercrystals import RegularSuperCrystals
                sage: RegularSuperCrystals().TensorProducts().extra_super_categories()
                [Category of regular super crystals]
            """
            return [self.base_category()]
