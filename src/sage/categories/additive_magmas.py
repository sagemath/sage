# sage_setup: distribution = sagemath-categories
r"""
Additive magmas
"""
# ****************************************************************************
#  Copyright (C) 2010-2014 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.misc.lazy_import import LazyImport
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.category_singleton import Category_singleton
from sage.categories.algebra_functor import AlgebrasCategory
from sage.categories.cartesian_product import CartesianProductsCategory
from sage.categories.homsets import HomsetsCategory
from sage.categories.with_realizations import WithRealizationsCategory
from sage.categories.sets_cat import Sets


class AdditiveMagmas(Category_singleton):
    """
    The category of additive magmas.

    An additive magma is a set endowed with a binary operation `+`.

    EXAMPLES::

        sage: AdditiveMagmas()
        Category of additive magmas
        sage: AdditiveMagmas().super_categories()
        [Category of sets]
        sage: AdditiveMagmas().all_super_categories()
        [Category of additive magmas,
         Category of sets,
         Category of sets with partial maps,
         Category of objects]

    The following axioms are defined by this category::

        sage: AdditiveMagmas().AdditiveAssociative()
        Category of additive semigroups
        sage: AdditiveMagmas().AdditiveUnital()
        Category of additive unital additive magmas
        sage: AdditiveMagmas().AdditiveCommutative()
        Category of additive commutative additive magmas
        sage: AdditiveMagmas().AdditiveUnital().AdditiveInverse()
        Category of additive inverse additive unital additive magmas
        sage: C = AdditiveMagmas().AdditiveAssociative().AdditiveCommutative(); C
        Category of commutative additive semigroups
        sage: C.AdditiveUnital()
        Category of commutative additive monoids
        sage: C.AdditiveUnital().AdditiveInverse()
        Category of commutative additive groups

    TESTS::

        sage: C = AdditiveMagmas()
        sage: TestSuite(C).run()
    """

    def super_categories(self):
        """
        EXAMPLES::

            sage: AdditiveMagmas().super_categories()
            [Category of sets]
        """
        return [Sets()]

    class SubcategoryMethods:

        @cached_method
        def AdditiveAssociative(self):
            r"""
            Return the full subcategory of the additive associative
            objects of ``self``.

            An :class:`additive magma <AdditiveMagmas>` `M` is
            *associative* if, for all `x,y,z \in M`,

            .. MATH:: x + (y + z) = (x + y) + z

            .. SEEALSO:: :wikipedia:`Associative_property`

            EXAMPLES::

                sage: AdditiveMagmas().AdditiveAssociative()
                Category of additive semigroups

            TESTS::

                sage: TestSuite(AdditiveMagmas().AdditiveAssociative()).run()
                sage: Rings().AdditiveAssociative.__module__
                'sage.categories.additive_magmas'
            """
            return self._with_axiom('AdditiveAssociative')

        @cached_method
        def AdditiveCommutative(self):
            r"""
            Return the full subcategory of the commutative objects of ``self``.

            An :class:`additive magma <AdditiveMagmas>` `M` is
            *commutative* if, for all `x,y \in M`,

            .. MATH:: x + y = y + x

            .. SEEALSO:: :wikipedia:`Commutative_property`

            EXAMPLES::

                sage: AdditiveMagmas().AdditiveCommutative()
                Category of additive commutative additive magmas
                sage: C = AdditiveMagmas().AdditiveAssociative().AdditiveUnital()
                sage: C.AdditiveCommutative()
                Category of commutative additive monoids
                sage: C.AdditiveCommutative() is CommutativeAdditiveMonoids()
                True

            TESTS::

                sage: TestSuite(AdditiveMagmas().AdditiveCommutative()).run()
                sage: Rings().AdditiveCommutative.__module__
                'sage.categories.additive_magmas'
            """
            return self._with_axiom('AdditiveCommutative')

        @cached_method
        def AdditiveUnital(self):
            r"""
            Return the subcategory of the unital objects of ``self``.

            An :class:`additive magma <AdditiveMagmas>` `M` is *unital*
            if it admits an element `0`, called *neutral element*,
            such that for all `x \in M`,

            .. MATH:: 0 + x = x + 0 = x

            This element is necessarily unique, and should be provided
            as ``M.zero()``.

            .. SEEALSO:: :wikipedia:`Unital_magma#unital`

            EXAMPLES::

                sage: AdditiveMagmas().AdditiveUnital()
                Category of additive unital additive magmas
                sage: from sage.categories.additive_semigroups import AdditiveSemigroups
                sage: AdditiveSemigroups().AdditiveUnital()
                Category of additive monoids
                sage: CommutativeAdditiveMonoids().AdditiveUnital()
                Category of commutative additive monoids

            TESTS::

                sage: TestSuite(AdditiveMagmas().AdditiveUnital()).run()
                sage: CommutativeAdditiveSemigroups().AdditiveUnital.__module__
                'sage.categories.additive_magmas'
            """
            return self._with_axiom("AdditiveUnital")

    AdditiveAssociative = LazyImport('sage.categories.additive_semigroups', 'AdditiveSemigroups', at_startup=True)

    class ParentMethods:

        def summation(self, x, y):
            r"""
            Return the sum of ``x`` and ``y``.

            The binary addition operator of this additive magma.

            INPUT:

            - ``x``, ``y`` -- elements of this additive magma

            EXAMPLES::

                sage: S = CommutativeAdditiveSemigroups().example()
                sage: (a,b,c,d) = S.additive_semigroup_generators()
                sage: S.summation(a, b)
                a + b

            A parent in ``AdditiveMagmas()`` must
            either implement :meth:`.summation` in the parent class or
            ``_add_`` in the element class. By default, the addition
            method on elements ``x._add_(y)`` calls
            ``S.summation(x,y)``, and reciprocally.

            As a bonus effect, ``S.summation`` by itself models the
            binary function from ``S`` to ``S``::

                sage: bin = S.summation
                sage: bin(a,b)
                a + b

            Here, ``S.summation`` is just a bound method. Whenever
            possible, it is recommended to enrich ``S.summation`` with
            extra mathematical structure. Lazy attributes can come
            handy for this.

            .. TODO:: Add an example.
            """
            return x + y

        summation_from_element_class_add = summation

        def __init_extra__(self):
            """
            TESTS::

                sage: S = CommutativeAdditiveSemigroups().example()
                sage: (a,b,c,d) = S.additive_semigroup_generators()
                sage: a + b # indirect doctest
                a + b
                sage: a.__class__._add_ == a.__class__._add_parent
                True
            """
            # This should instead register the summation to the coercion model
            # But this is not yet implemented in the coercion model
            if (self.summation != self.summation_from_element_class_add) and hasattr(self, "element_class") and hasattr(self.element_class, "_add_parent"):
                self.element_class._add_ = self.element_class._add_parent

        def addition_table(self, names='letters', elements=None):
            r"""
            Return a table describing the addition operation.

            .. NOTE::

                The order of the elements in the row and column
                headings is equal to the order given by the table's
                :meth:`~sage.matrix.operation_table.OperationTable.column_keys`
                method.  The association can also be retrieved with the
                :meth:`~sage.matrix.operation_table.OperationTable.translation`
                method.

            INPUT:

            - ``names`` -- the type of names used:

              * ``'letters'`` -- lowercase ASCII letters are used
                for a base 26 representation of the elements'
                positions in the list given by
                :meth:`~sage.matrix.operation_table.OperationTable.column_keys`,
                padded to a common width with leading 'a's.
              * ``'digits'`` -- base 10 representation of the
                elements' positions in the list given by
                :meth:`~sage.matrix.operation_table.OperationTable.column_keys`,
                padded to a common width with leading zeros.
              * ``'elements'`` -- the string representations
                of the elements themselves.
              * a list - a list of strings, where the length
                of the list equals the number of elements.

            - ``elements`` -- (default: ``None``)  A list of
              elements of the additive magma, in forms that
              can be coerced into the structure, eg. their
              string representations. This may be used to
              impose an alternate ordering on the elements,
              perhaps when this is used in the context of a
              particular structure. The default is to use
              whatever ordering the ``S.list`` method returns.
              Or the ``elements`` can be a subset which is
              closed under the operation. In particular,
              this can be used when the base set is infinite.

            OUTPUT:

            The addition table as an object of the class
            :class:`~sage.matrix.operation_table.OperationTable`
            which defines several methods for manipulating and
            displaying the table.  See the documentation there
            for full details to supplement the documentation
            here.

            EXAMPLES:

            All that is required is that an algebraic structure
            has an addition defined.The default is to represent
            elements as lowercase ASCII letters.  ::

                sage: R = IntegerModRing(5)
                sage: R.addition_table()                                                # needs sage.modules
                +  a b c d e
                 +----------
                a| a b c d e
                b| b c d e a
                c| c d e a b
                d| d e a b c
                e| e a b c d

            The ``names`` argument allows displaying the elements in
            different ways.  Requesting ``elements`` will use the
            representation of the elements of the set.  Requesting
            ``digits`` will include leading zeros as padding.  ::

                sage: R = IntegerModRing(11)
                sage: P = R.addition_table(names='elements'); P                         # needs sage.modules
                 +   0  1  2  3  4  5  6  7  8  9 10
                  +---------------------------------
                 0|  0  1  2  3  4  5  6  7  8  9 10
                 1|  1  2  3  4  5  6  7  8  9 10  0
                 2|  2  3  4  5  6  7  8  9 10  0  1
                 3|  3  4  5  6  7  8  9 10  0  1  2
                 4|  4  5  6  7  8  9 10  0  1  2  3
                 5|  5  6  7  8  9 10  0  1  2  3  4
                 6|  6  7  8  9 10  0  1  2  3  4  5
                 7|  7  8  9 10  0  1  2  3  4  5  6
                 8|  8  9 10  0  1  2  3  4  5  6  7
                 9|  9 10  0  1  2  3  4  5  6  7  8
                10| 10  0  1  2  3  4  5  6  7  8  9

                sage: T = R.addition_table(names='digits'); T                           # needs sage.modules
                +  00 01 02 03 04 05 06 07 08 09 10
                  +---------------------------------
                00| 00 01 02 03 04 05 06 07 08 09 10
                01| 01 02 03 04 05 06 07 08 09 10 00
                02| 02 03 04 05 06 07 08 09 10 00 01
                03| 03 04 05 06 07 08 09 10 00 01 02
                04| 04 05 06 07 08 09 10 00 01 02 03
                05| 05 06 07 08 09 10 00 01 02 03 04
                06| 06 07 08 09 10 00 01 02 03 04 05
                07| 07 08 09 10 00 01 02 03 04 05 06
                08| 08 09 10 00 01 02 03 04 05 06 07
                09| 09 10 00 01 02 03 04 05 06 07 08
                10| 10 00 01 02 03 04 05 06 07 08 09

            Specifying the elements in an alternative order can provide
            more insight into how the operation behaves.  ::

                sage: S = IntegerModRing(7)
                sage: elts = [0, 3, 6, 2, 5, 1, 4]
                sage: S.addition_table(elements=elts)                                   # needs sage.modules
                +  a b c d e f g
                 +--------------
                a| a b c d e f g
                b| b c d e f g a
                c| c d e f g a b
                d| d e f g a b c
                e| e f g a b c d
                f| f g a b c d e
                g| g a b c d e f

            The ``elements`` argument can be used to provide
            a subset of the elements of the structure.  The subset
            must be closed under the operation.  Elements need only
            be in a form that can be coerced into the set.  The
            ``names`` argument can also be used to request that
            the elements be represented with their usual string
            representation.  ::

                sage: T = IntegerModRing(12)
                sage: elts = [0, 3, 6, 9]
                sage: T.addition_table(names='elements', elements=elts)                 # needs sage.modules
                +  0 3 6 9
                 +--------
                0| 0 3 6 9
                3| 3 6 9 0
                6| 6 9 0 3
                9| 9 0 3 6

            The table returned can be manipulated in various ways.  See
            the documentation for
            :class:`~sage.matrix.operation_table.OperationTable` for more
            comprehensive documentation. ::

                sage: # needs sage.modules
                sage: R = IntegerModRing(3)
                sage: T = R.addition_table()
                sage: T.column_keys()
                (0, 1, 2)
                sage: sorted(T.translation().items())
                [('a', 0), ('b', 1), ('c', 2)]
                sage: T.change_names(['x', 'y', 'z'])
                sage: sorted(T.translation().items())
                [('x', 0), ('y', 1), ('z', 2)]
                sage: T
                +  x y z
                 +------
                x| x y z
                y| y z x
                z| z x y
            """
            from sage.matrix.operation_table import OperationTable
            import operator
            return OperationTable(self, operation=operator.add,
                                  names=names, elements=elements)

    class ElementMethods:

        @abstract_method(optional=True)
        def _add_(self, right):
            """
            Return the sum of ``self`` and ``right``.

            INPUT:

            - ``self``, ``right`` -- two elements with the same parent

            OUTPUT: an element of the same parent

            EXAMPLES::

                sage: F = CommutativeAdditiveSemigroups().example()
                sage: (a,b,c,d) = F.additive_semigroup_generators()
                sage: a._add_(b)
                a + b
            """

        def _add_parent(self, other):
            r"""
            Return the sum of the two elements, calculated using
            the ``summation`` method of the parent.

            This is the default implementation of _add_ if
            ``summation`` is implemented in the parent.

            INPUT:

            - ``other`` -- an element of the parent of ``self``

            OUTPUT: an element of the parent of ``self``

            EXAMPLES::

                sage: S = CommutativeAdditiveSemigroups().example()
                sage: (a,b,c,d) = S.additive_semigroup_generators()
                sage: a._add_parent(b)
                a + b
            """
            return self.parent().summation(self, other)

    class Homsets(HomsetsCategory):
        def extra_super_categories(self):
            """
            Implement the fact that a homset between two magmas is a magma.

            EXAMPLES::

                sage: AdditiveMagmas().Homsets().extra_super_categories()
                [Category of additive magmas]
                sage: AdditiveMagmas().Homsets().super_categories()
                [Category of additive magmas, Category of homsets]
            """
            return [AdditiveMagmas()]

    class CartesianProducts(CartesianProductsCategory):
        def extra_super_categories(self):
            """
            Implement the fact that a Cartesian product of additive magmas is
            an additive magma.

            EXAMPLES::

                sage: C = AdditiveMagmas().CartesianProducts()
                sage: C.extra_super_categories()
                [Category of additive magmas]
                sage: C.super_categories()
                [Category of additive magmas, Category of Cartesian products of sets]
                sage: C.axioms()
                frozenset()
            """
            return [AdditiveMagmas()]

        class ElementMethods:
            def _add_(self, right):
                r"""
                EXAMPLES::

                    sage: # needs sage.rings.finite_rings
                    sage: G5 = GF(5); G8 = GF(4, 'x'); GG = G5.cartesian_product(G8)
                    sage: e = GG((G5(1), G8.primitive_element())); e
                    (1, x)
                    sage: e + e
                    (2, 0)
                    sage: e = groups.misc.AdditiveCyclic(8)                             # needs sage.groups
                    sage: x = e.cartesian_product(e)((e(1), e(2)))
                    sage: x
                    (1, 2)
                    sage: 4 * x
                    (4, 0)
                """
                return self.parent()._cartesian_product_of_elements(
                    x + y for x, y in zip(self.cartesian_factors(),
                                          right.cartesian_factors()))

    class Algebras(AlgebrasCategory):

        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: AdditiveMagmas().Algebras(QQ).extra_super_categories()
                [Category of magmatic algebras with basis over Rational Field]

                sage: AdditiveMagmas().Algebras(QQ).super_categories()
                [Category of magmatic algebras with basis over Rational Field,
                 Category of set algebras over Rational Field]
            """
            from sage.categories.magmatic_algebras import MagmaticAlgebras
            return [MagmaticAlgebras(self.base_ring()).WithBasis()]

        class ParentMethods:

            @cached_method
            def algebra_generators(self):
                r"""
                The generators of this algebra, as per
                :meth:`MagmaticAlgebras.ParentMethods.algebra_generators()
                <.magmatic_algebras.MagmaticAlgebras.ParentMethods.algebra_generators>`.

                They correspond to the generators of the additive semigroup.

                EXAMPLES::

                    sage: S = CommutativeAdditiveSemigroups().example(); S
                    An example of a commutative semigroup:
                     the free commutative semigroup generated by ('a', 'b', 'c', 'd')
                    sage: A = S.algebra(QQ)                                             # needs sage.modules
                    sage: A.algebra_generators()                                        # needs sage.modules
                    Family (B[a], B[b], B[c], B[d])

                .. TODO::

                    This doctest does not actually test this method,
                    but rather the method of the same name for
                    ``AdditiveSemigroups``. Find a better doctest!
                """
                return self.basis().keys().additive_semigroup_generators().map(self.monomial)

            def product_on_basis(self, g1, g2):
                r"""
                Product, on basis elements, as per
                :meth:`MagmaticAlgebras.WithBasis.ParentMethods.product_on_basis()
                <.magmatic_algebras.MagmaticAlgebras.WithBasis.ParentMethods.product_on_basis>`.

                The product of two basis elements is induced by the
                addition of the corresponding elements of the group.

                EXAMPLES::

                    sage: S = CommutativeAdditiveSemigroups().example(); S
                    An example of a commutative semigroup:
                     the free commutative semigroup generated by ('a', 'b', 'c', 'd')
                    sage: A = S.algebra(QQ)                                             # needs sage.modules
                    sage: a, b, c, d = A.algebra_generators()                           # needs sage.modules
                    sage: a * d * b                                                     # needs sage.modules
                    B[a + b + d]

                .. TODO::

                    This doctest does not actually test this method,
                    but rather the method of the same name for
                    ``AdditiveSemigroups``. Find a better doctest!
                """
                return self.monomial(g1 + g2)

    class AdditiveCommutative(CategoryWithAxiom):
        class CartesianProducts(CartesianProductsCategory):
            def extra_super_categories(self):
                """
                Implement the fact that a Cartesian product of commutative
                additive magmas is a commutative additive magma.

                EXAMPLES::

                    sage: C = AdditiveMagmas().AdditiveCommutative().CartesianProducts()
                    sage: C.extra_super_categories()
                    [Category of additive commutative additive magmas]
                    sage: C.axioms()
                    frozenset({'AdditiveCommutative'})
                """
                return [AdditiveMagmas().AdditiveCommutative()]

        class Algebras(AlgebrasCategory):
            def extra_super_categories(self):
                """
                Implement the fact that the algebra of a commutative additive
                magmas is commutative.

                EXAMPLES::

                    sage: C = AdditiveMagmas().AdditiveCommutative().Algebras(QQ)
                    sage: C.extra_super_categories()
                    [Category of commutative magmas]

                    sage: C.super_categories()
                    [Category of additive magma algebras over Rational Field,
                     Category of commutative magmas]
                """
                from sage.categories.magmas import Magmas
                return [Magmas().Commutative()]

    class AdditiveUnital(CategoryWithAxiom):

        def additional_structure(self):
            r"""
            Return whether ``self`` is a structure category.

            .. SEEALSO:: :meth:`Category.additional_structure`

            The category of unital additive magmas defines the zero as
            additional structure, and this zero shall be preserved by
            morphisms.

            EXAMPLES::

                sage: AdditiveMagmas().AdditiveUnital().additional_structure()
                Category of additive unital additive magmas
            """
            return self

        class SubcategoryMethods:

            @cached_method
            def AdditiveInverse(self):
                r"""
                Return the full subcategory of the additive inverse objects
                of ``self``.

                An inverse :class:`additive magma <AdditiveMagmas>` is
                a :class:`unital additive magma <AdditiveMagmas.Unital>`
                such that every element admits both an additive
                inverse on the left and on the right. Such an additive
                magma is also called an *additive loop*.

                .. SEEALSO::

                    :wikipedia:`Inverse_element`, :wikipedia:`Quasigroup`

                EXAMPLES::

                    sage: AdditiveMagmas().AdditiveUnital().AdditiveInverse()
                    Category of additive inverse additive unital additive magmas
                    sage: from sage.categories.additive_monoids import AdditiveMonoids
                    sage: AdditiveMonoids().AdditiveInverse()
                    Category of additive groups

                TESTS::

                    sage: TestSuite(AdditiveMagmas().AdditiveUnital().AdditiveInverse()).run()
                    sage: CommutativeAdditiveMonoids().AdditiveInverse.__module__
                    'sage.categories.additive_magmas'
                """
                return self._with_axiom("AdditiveInverse")

        class ParentMethods:

            def _test_zero(self, **options):
                r"""
                Test that ``self.zero()`` is an element of ``self`` and
                is neutral for the addition.

                INPUT:

                - ``options`` -- any keyword arguments accepted
                  by :meth:`_tester`

                EXAMPLES:

                By default, this method tests only the elements returned by
                ``self.some_elements()``::

                    sage: S = CommutativeAdditiveMonoids().example()
                    sage: S._test_zero()

                However, the elements tested can be customized with the
                ``elements`` keyword argument::

                    sage: (a,b,c,d) = S.additive_semigroup_generators()
                    sage: S._test_zero(elements = (a, a+c))

                See the documentation for :class:`TestSuite` for
                more information.
                """
                tester = self._tester(**options)
                zero = self.zero()
                # TODO: also call is_zero once it will work
                tester.assertTrue(self.is_parent_of(zero))
                for x in tester.some_elements():
                    tester.assertEqual(x + zero, x)
                # Check that zero is immutable if it looks like we can:
                if hasattr(zero, "is_immutable"):
                    tester.assertEqual(zero.is_immutable(), True)
                if hasattr(zero, "is_mutable"):
                    tester.assertEqual(zero.is_mutable(), False)
                # Check that bool behave consistently on zero
                tester.assertFalse(bool(self.zero()))

            @cached_method
            def zero(self):
                """
                Return the zero of this additive magma, that is the unique
                neutral element for `+`.

                The default implementation is to coerce ``0`` into ``self``.

                It is recommended to override this method because the
                coercion from the integers:

                - is not always meaningful (except for `0`), and
                - often uses ``self.zero()`` otherwise.

                EXAMPLES::

                    sage: S = CommutativeAdditiveMonoids().example()
                    sage: S.zero()
                    0
                """
                # TODO: add a test that actually exercise this default implementation
                return self(0)

            def is_empty(self):
                r"""
                Return whether this set is empty.

                Since this set is an additive magma it has a zero element and
                hence is not empty. This method thus always returns ``False``.

                EXAMPLES::

                    sage: # needs sage.modules
                    sage: A = AdditiveAbelianGroup([3, 3])
                    sage: A in AdditiveMagmas()
                    True
                    sage: A.is_empty()
                    False

                    sage: B = CommutativeAdditiveMonoids().example()
                    sage: B.is_empty()
                    False

                TESTS:

                We check that the method ``is_empty`` is inherited from this
                category in both examples above::

                    sage: A.is_empty.__module__                                         # needs sage.modules
                    'sage.categories.additive_magmas'
                    sage: B.is_empty.__module__
                    'sage.categories.additive_magmas'
                """
                return False

        class ElementMethods:
            # TODO: merge with the implementation in Element which currently
            # overrides this one, and further requires self.parent()(0) to work.
            #
            # def is_zero(self):
            #     """
            #     Returns whether self is the zero of the magma
            #
            #     The default implementation, is to compare with ``self.zero()``.
            #
            #     TESTS::
            #
            #         sage: S = CommutativeAdditiveMonoids().example()
            #         sage: S.zero().is_zero()
            #         True
            #         sage: S("aa").is_zero()
            #         False
            #     """
            #     return self == self.parent().zero()

            @abstract_method
            def __bool__(self):
                """
                Return whether ``self`` is not zero.

                All parents in the category ``CommutativeAdditiveMonoids()``
                should implement this method.

                .. NOTE:: This is currently not useful because this method is
                   overridden by ``Element``.

                TESTS::

                    sage: S = CommutativeAdditiveMonoids().example()
                    sage: bool(S.zero())
                    False
                    sage: bool(S.an_element())
                    True
                """

            def _test_nonzero_equal(self, **options):
                r"""
                Test that ``.__bool__()`` behave consistently
                with `` == 0``.

                TESTS::

                    sage: S = CommutativeAdditiveMonoids().example()
                    sage: S.zero()._test_nonzero_equal()
                    sage: S.an_element()._test_nonzero_equal()
                """
                tester = self._tester(**options)
                tester.assertEqual(bool(self), self != self.parent().zero())
                tester.assertEqual(not self, self == self.parent().zero())

            def _sub_(left, right):
                r"""
                Default implementation of difference.

                EXAMPLES::

                    sage: F = CombinatorialFreeModule(QQ, ['a', 'b'])                   # needs sage.modules
                    sage: a, b = F.basis()                                              # needs sage.modules
                    sage: a - b                                                         # needs sage.modules
                    B['a'] - B['b']

                TESTS:

                Check that :issue:`18275` is fixed::

                    sage: C = GF(5).cartesian_product(GF(5))
                    sage: C.one() - C.one()
                    (0, 0)
                """
                return left + (-right)

            def __neg__(self):
                """
                Return the negation of ``self``, if it exists.

                This top-level implementation delegates the job to
                ``_neg_``, for those additive unital magmas which may
                choose to implement it instead of ``__neg__`` for
                consistency.

                EXAMPLES::

                    sage: F = CombinatorialFreeModule(QQ, ['a', 'b'])                   # needs sage.modules
                    sage: a, b = F.basis()                                              # needs sage.modules
                    sage: -b                                                            # needs sage.modules
                    -B['b']

                TESTS::

                    sage: # needs sage.modules
                    sage: F = CombinatorialFreeModule(ZZ, ['a', 'b'])
                    sage: a, b = F.gens()
                    sage: FF = cartesian_product((F, F))
                    sage: x = cartesian_product([a, 2*a-3*b]); x
                    B[(0, 'a')] + 2*B[(1, 'a')] - 3*B[(1, 'b')]
                    sage: x.parent() is FF
                    True
                    sage: -x
                    -B[(0, 'a')] - 2*B[(1, 'a')] + 3*B[(1, 'b')]
                """
                return self._neg_()

        class Homsets(HomsetsCategory):
            def extra_super_categories(self):
                """
                Implement the fact that a homset between two unital additive
                magmas is a unital additive magma.

                EXAMPLES::

                    sage: AdditiveMagmas().AdditiveUnital().Homsets().extra_super_categories()
                    [Category of additive unital additive magmas]
                    sage: AdditiveMagmas().AdditiveUnital().Homsets().super_categories()
                    [Category of additive unital additive magmas, Category of homsets]
                """
                return [AdditiveMagmas().AdditiveUnital()]

            class ParentMethods:

                @cached_method
                def zero(self):
                    """
                    EXAMPLES::

                        sage: R = QQ['x']
                        sage: H = Hom(ZZ, R, AdditiveMagmas().AdditiveUnital())
                        sage: f = H.zero()
                        sage: f
                        Generic morphism:
                          From: Integer Ring
                          To:   Univariate Polynomial Ring in x over Rational Field
                        sage: f(3)
                        0
                        sage: f(3) is R.zero()
                        True

                    TESTS:

                        sage: TestSuite(f).run()
                    """
                    from sage.misc.constant_function import ConstantFunction
                    return self(ConstantFunction(self.codomain().zero()))

        class AdditiveInverse(CategoryWithAxiom):
            class CartesianProducts(CartesianProductsCategory):
                def extra_super_categories(self):
                    """
                    Implement the fact that a Cartesian product of additive magmas
                    with inverses is an additive magma with inverse.

                    EXAMPLES::

                        sage: C = AdditiveMagmas().AdditiveUnital().AdditiveInverse().CartesianProducts()
                        sage: C.extra_super_categories()
                        [Category of additive inverse additive unital additive magmas]
                        sage: sorted(C.axioms())
                        ['AdditiveInverse', 'AdditiveUnital']
                    """
                    return [AdditiveMagmas().AdditiveUnital().AdditiveInverse()]

                class ElementMethods:
                    def _neg_(self):
                        """
                        Return the negation of ``self``.

                        EXAMPLES::

                           sage: x = cartesian_product((GF(7)(2), 17)); x
                           (2, 17)
                           sage: -x
                           (5, -17)

                        TESTS::

                           sage: C = AdditiveMagmas().AdditiveUnital().AdditiveInverse().CartesianProducts()
                           sage: x.parent() in C
                           True
                        """
                        return self.parent()._cartesian_product_of_elements(
                            [-x for x in self.cartesian_factors()])

        class CartesianProducts(CartesianProductsCategory):
            def extra_super_categories(self):
                """
                Implement the fact that a Cartesian product of unital additive
                magmas is a unital additive magma.

                EXAMPLES::

                    sage: C = AdditiveMagmas().AdditiveUnital().CartesianProducts()
                    sage: C.extra_super_categories()
                    [Category of additive unital additive magmas]
                    sage: C.axioms()
                    frozenset({'AdditiveUnital'})
                """
                return [AdditiveMagmas().AdditiveUnital()]

            class ParentMethods:
                def zero(self):
                    r"""
                    Return the zero of this group.

                    EXAMPLES::

                        sage: GF(8, 'x').cartesian_product(GF(5)).zero()                # needs sage.rings.finite_rings
                        (0, 0)
                    """
                    return self._cartesian_product_of_elements(
                        _.zero() for _ in self.cartesian_factors())

        class Algebras(AlgebrasCategory):

            def extra_super_categories(self):
                """
                EXAMPLES::

                    sage: C = AdditiveMagmas().AdditiveUnital().Algebras(QQ)
                    sage: C.extra_super_categories()
                    [Category of unital magmas]

                    sage: C.super_categories()
                    [Category of additive magma algebras over Rational Field,
                     Category of unital algebras with basis over Rational Field]
                """
                from sage.categories.magmas import Magmas
                return [Magmas().Unital()]

            class ParentMethods:

                @cached_method
                def one_basis(self):
                    """
                    Return the zero of this additive magma, which index the
                    one of this algebra, as per
                    :meth:`AlgebrasWithBasis.ParentMethods.one_basis()
                    <sage.categories.algebras_with_basis.AlgebrasWithBasis.ParentMethods.one_basis>`.

                    EXAMPLES::

                        sage: # needs sage.modules
                        sage: S = CommutativeAdditiveMonoids().example(); S
                        An example of a commutative monoid:
                         the free commutative monoid generated by ('a', 'b', 'c', 'd')
                        sage: A = S.algebra(ZZ)
                        sage: A.one_basis()
                        0
                        sage: A.one()
                        B[0]
                        sage: A(3)
                        3*B[0]
                    """
                    return self.basis().keys().zero()

        class WithRealizations(WithRealizationsCategory):

            class ParentMethods:

                def zero(self):
                    r"""
                    Return the zero of this unital additive magma.

                    This default implementation returns the zero of the
                    realization of ``self`` given by
                    :meth:`~Sets.WithRealizations.ParentMethods.a_realization`.

                    EXAMPLES::

                        sage: A = Sets().WithRealizations().example(); A                # needs sage.modules
                        The subset algebra of {1, 2, 3} over Rational Field
                        sage: A.zero.__module__                                         # needs sage.modules
                        'sage.categories.additive_magmas'
                        sage: A.zero()                                                  # needs sage.modules
                        0

                    TESTS::

                        sage: A.zero() is A.a_realization().zero()                      # needs sage.modules
                        True
                        sage: A._test_zero()                                            # needs sage.modules
                    """
                    return self.a_realization().zero()
