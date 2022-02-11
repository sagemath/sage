r"""
Algebras With Basis
"""
#*****************************************************************************
#  Copyright (C) 2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2013 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.lazy_import import LazyImport
from sage.categories.tensor import TensorProductsCategory, tensor
from sage.categories.cartesian_product import CartesianProductsCategory
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from .unital_algebras import UnitalAlgebras


class AlgebrasWithBasis(CategoryWithAxiom_over_base_ring):
    """
    The category of algebras with a distinguished basis.

    EXAMPLES::

        sage: C = AlgebrasWithBasis(QQ); C
        Category of algebras with basis over Rational Field
        sage: sorted(C.super_categories(), key=str)
        [Category of algebras over Rational Field,
         Category of unital algebras with basis over Rational Field]

    We construct a typical parent in this category, and do some
    computations with it::

        sage: A = C.example(); A
        An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field

        sage: A.category()
        Category of algebras with basis over Rational Field

        sage: A.one_basis()
        word:
        sage: A.one()
        B[word: ]

        sage: A.base_ring()
        Rational Field
        sage: A.basis().keys()
        Finite words over {'a', 'b', 'c'}

        sage: (a,b,c) = A.algebra_generators()
        sage: a^3, b^2
        (B[word: aaa], B[word: bb])
        sage: a*c*b
        B[word: acb]

        sage: A.product
        <bound method MagmaticAlgebras.WithBasis.ParentMethods._product_from_product_on_basis_multiply of
         An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field>
        sage: A.product(a*b,b)
        B[word: abb]

        sage: TestSuite(A).run(verbose=True)
        running ._test_additive_associativity() . . . pass
        running ._test_an_element() . . . pass
        running ._test_associativity() . . . pass
        running ._test_cardinality() . . . pass
        running ._test_category() . . . pass
        running ._test_characteristic() . . . pass
        running ._test_construction() . . . pass
        running ._test_distributivity() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
          running ._test_new() . . . pass
          running ._test_nonzero_equal() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          pass
        running ._test_elements_eq_reflexive() . . . pass
        running ._test_elements_eq_symmetric() . . . pass
        running ._test_elements_eq_transitive() . . . pass
        running ._test_elements_neq() . . . pass
        running ._test_eq() . . . pass
        running ._test_new() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_one() . . . pass
        running ._test_pickling() . . . pass
        running ._test_prod() . . . pass
        running ._test_some_elements() . . . pass
        running ._test_zero() . . . pass
        sage: A.__class__
        <class 'sage.categories.examples.algebras_with_basis.FreeAlgebra_with_category'>
        sage: A.element_class
        <class 'sage.categories.examples.algebras_with_basis.FreeAlgebra_with_category.element_class'>

    Please see the source code of `A` (with ``A??``) for how to
    implement other algebras with basis.

    TESTS::

        sage: TestSuite(AlgebrasWithBasis(QQ)).run()
    """

    def example(self, alphabet = ('a','b','c')):
        """
        Return an example of algebra with basis.

        EXAMPLES::

            sage: AlgebrasWithBasis(QQ).example()
            An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field

        An other set of generators can be specified as optional argument::

            sage: AlgebrasWithBasis(QQ).example((1,2,3))
            An example of an algebra with basis: the free algebra on the generators (1, 2, 3) over Rational Field
        """
        from sage.categories.examples.algebras_with_basis import Example
        return Example(self.base_ring(), alphabet)

    Filtered = LazyImport('sage.categories.filtered_algebras_with_basis', 'FilteredAlgebrasWithBasis')
    FiniteDimensional = LazyImport('sage.categories.finite_dimensional_algebras_with_basis', 'FiniteDimensionalAlgebrasWithBasis', at_startup=True)
    Graded = LazyImport('sage.categories.graded_algebras_with_basis', 'GradedAlgebrasWithBasis')
    Super = LazyImport('sage.categories.super_algebras_with_basis', 'SuperAlgebrasWithBasis')

    class ParentMethods:

        # For backward compatibility
        one = UnitalAlgebras.WithBasis.ParentMethods.one

        def hochschild_complex(self, M):
            """
            Return the Hochschild complex of ``self`` with coefficients
            in ``M``.

            .. SEEALSO::

                :class:`~sage.homology.hochschild_complex.HochschildComplex`

            EXAMPLES::

                sage: R.<x> = QQ[]
                sage: A = algebras.DifferentialWeyl(R)
                sage: H = A.hochschild_complex(A)

                sage: SGA = SymmetricGroupAlgebra(QQ, 3)
                sage: T = SGA.trivial_representation()
                sage: H = SGA.hochschild_complex(T)
            """
            from sage.homology.hochschild_complex import HochschildComplex
            return HochschildComplex(self, M)

    class ElementMethods:

        def __invert__(self):
            """
            Return the inverse of ``self`` if ``self`` is a multiple of one,
            and one is in the basis of this algebra. Otherwise throws
            an error.

            Caveat: this generic implementation is not complete; there
            may be invertible elements in the algebra that can't be
            inversed this way. It is correct though for graded
            connected algebras with basis.

            .. WARNING::

                This might produce a result which does not belong to
                the parent of ``self``, yet believes to do so. For
                instance, inverting 2 times the unity will produce 1/2
                times the unity, even if 1/2 is not in the base ring.
                Handle with care.

            EXAMPLES::

                sage: C = AlgebrasWithBasis(QQ).example()
                sage: x = C(2); x
                2*B[word: ]
                sage: ~x
                1/2*B[word: ]
                sage: a = C.algebra_generators().first(); a
                B[word: a]
                sage: ~a
                Traceback (most recent call last):
                ...
                ValueError: cannot invert self (= B[word: a])
            """
            # FIXME: make this generic
            mcs = self.monomial_coefficients(copy=False)
            one = self.parent().one_basis()
            if len(mcs) == 1 and one in mcs:
                return self.parent().term(one, ~mcs[one])
            else:
                raise ValueError("cannot invert self (= %s)" % self)

    class CartesianProducts(CartesianProductsCategory):
        """
        The category of algebras with basis, constructed as Cartesian
        products of algebras with basis.

        Note: this construction give the direct products of algebras with basis.
        See comment in :class:`Algebras.CartesianProducts
        <sage.categories.algebras.Algebras.CartesianProducts>`
        """

        def extra_super_categories(self):
            """
            A Cartesian product of algebras with basis is endowed with
            a natural algebra with basis structure.

            EXAMPLES::

                sage: AlgebrasWithBasis(QQ).CartesianProducts().extra_super_categories()
                [Category of algebras with basis over Rational Field]
                sage: AlgebrasWithBasis(QQ).CartesianProducts().super_categories()
                [Category of algebras with basis over Rational Field,
                 Category of Cartesian products of algebras over Rational Field,
                 Category of Cartesian products of vector spaces with basis over Rational Field]
            """
            return [self.base_category()]

        class ParentMethods:
            @cached_method
            def one_from_cartesian_product_of_one_basis(self):
                """
                Returns the one of this Cartesian product of algebras, as per ``Monoids.ParentMethods.one``

                It is constructed as the Cartesian product of the ones of the
                summands, using their :meth:`~AlgebrasWithBasis.ParentMethods.one_basis` methods.

                This implementation does not require multiplication by
                scalars nor calling cartesian_product. This might help keeping
                things as lazy as possible upon initialization.

                EXAMPLES::

                    sage: A = AlgebrasWithBasis(QQ).example(); A
                    An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field
                    sage: A.one_basis()
                    word:

                    sage: B = cartesian_product((A, A, A))
                    sage: B.one_from_cartesian_product_of_one_basis()
                    B[(0, word: )] + B[(1, word: )] + B[(2, word: )]
                    sage: B.one()
                    B[(0, word: )] + B[(1, word: )] + B[(2, word: )]

                    sage: cartesian_product([SymmetricGroupAlgebra(QQ, 3), SymmetricGroupAlgebra(QQ, 4)]).one()
                    B[(0, [1, 2, 3])] + B[(1, [1, 2, 3, 4])]
                """
                return self.sum_of_monomials( zip( self._sets_keys(), (set.one_basis() for set in self._sets)) )

            @lazy_attribute
            def one(self):
                """
                TESTS::

                    sage: A = AlgebrasWithBasis(QQ).example(); A
                    An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field
                    sage: B = cartesian_product((A, A, A))
                    sage: B.one()
                    B[(0, word: )] + B[(1, word: )] + B[(2, word: )]
                """
                if all(hasattr(module, "one_basis") for module in self._sets):
                    return self.one_from_cartesian_product_of_one_basis
                else:
                    return NotImplemented

            # def product_on_basis(self, t1, t2):
            # would be easy to implement, but without a special
            # version of module morphism, this would not take
            # advantage of the block structure


    class TensorProducts(TensorProductsCategory):
        """
        The category of algebras with basis constructed by tensor product of algebras with basis
        """

        @cached_method
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: AlgebrasWithBasis(QQ).TensorProducts().extra_super_categories()
                [Category of algebras with basis over Rational Field]
                sage: AlgebrasWithBasis(QQ).TensorProducts().super_categories()
                [Category of algebras with basis over Rational Field,
                 Category of tensor products of algebras over Rational Field,
                 Category of tensor products of vector spaces with basis over Rational Field]
            """
            return [self.base_category()]

        class ParentMethods:
            """
            implements operations on tensor products of algebras with basis
            """

            @cached_method
            def one_basis(self):
                """
                Returns the index of the one of this tensor product of
                algebras, as per ``AlgebrasWithBasis.ParentMethods.one_basis``

                It is the tuple whose operands are the indices of the
                ones of the operands, as returned by their
                :meth:`.one_basis` methods.

                EXAMPLES::

                    sage: A = AlgebrasWithBasis(QQ).example(); A
                    An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field
                    sage: A.one_basis()
                    word:
                    sage: B = tensor((A, A, A))
                    sage: B.one_basis()
                    (word: , word: , word: )
                    sage: B.one()
                    B[word: ] # B[word: ] # B[word: ]
                """
                # FIXME: this method should be conditionally defined,
                # so that B.one_basis returns NotImplemented if not
                # all modules provide one_basis
                if all(hasattr(module, "one_basis") for module in self._sets):
                    return tuple(module.one_basis() for module in self._sets)
                else:
                    raise NotImplementedError

            def product_on_basis(self, t1, t2):
                """
                The product of the algebra on the basis, as per
                ``AlgebrasWithBasis.ParentMethods.product_on_basis``.

                EXAMPLES::

                    sage: A = AlgebrasWithBasis(QQ).example(); A
                    An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field
                    sage: (a,b,c) = A.algebra_generators()

                    sage: x = tensor( (a, b, c) ); x
                    B[word: a] # B[word: b] # B[word: c]
                    sage: y = tensor( (c, b, a) ); y
                    B[word: c] # B[word: b] # B[word: a]
                    sage: x*y
                    B[word: ac] # B[word: bb] # B[word: ca]

                    sage: x = tensor( ((a+2*b), c) )    ; x
                    B[word: a] # B[word: c] + 2*B[word: b] # B[word: c]
                    sage: y = tensor( (c,       a) ) + 1; y
                    B[word: ] # B[word: ] + B[word: c] # B[word: a]
                    sage: x*y
                    B[word: a] # B[word: c] + B[word: ac] # B[word: ca] + 2*B[word: b] # B[word: c] + 2*B[word: bc] # B[word: ca]


                TODO: optimize this implementation!
                """
                return tensor( (module.monomial(x1)*module.monomial(x2) for (module, x1, x2) in zip(self._sets, t1, t2)) ) #.

        class ElementMethods:
            """
            Implements operations on elements of tensor products of algebras with basis
            """
            pass

