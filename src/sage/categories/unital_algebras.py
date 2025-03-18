# sage_setup: distribution = sagemath-categories
r"""
Unital algebras
"""
# ****************************************************************************
#  Copyright (C) 2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.category import Category
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.commutative_additive_groups import CommutativeAdditiveGroups
from sage.categories.magmas import Magmas
from sage.categories.morphism import SetMorphism
from sage.categories.homset import Hom
from sage.categories.rings import Rings
from sage.categories.magmatic_algebras import MagmaticAlgebras
from sage.categories.cartesian_product import CartesianProductsCategory


class UnitalAlgebras(CategoryWithAxiom_over_base_ring):
    """
    The category of non-associative algebras over a given base ring.

    A non-associative algebra over a ring `R` is a module over `R`
    which is also a unital magma.

    .. WARNING::

        Until :issue:`15043` is implemented, :class:`Algebras` is the
        category of associative unital algebras; thus, unlike the name
        suggests, :class:`UnitalAlgebras` is not a subcategory of
        :class:`Algebras` but of
        :class:`~.magmatic_algebras.MagmaticAlgebras`.

    EXAMPLES::

        sage: from sage.categories.unital_algebras import UnitalAlgebras
        sage: C = UnitalAlgebras(ZZ); C
        Category of unital algebras over Integer Ring

    TESTS::

        sage: from sage.categories.magmatic_algebras import MagmaticAlgebras
        sage: C is MagmaticAlgebras(ZZ).Unital()
        True
        sage: TestSuite(C).run()
    """
    _base_category_class_and_axiom = (MagmaticAlgebras, "Unital")

    class ParentMethods:
        def from_base_ring(self, r):
            """
            Return the canonical embedding of ``r`` into ``self``.

            INPUT:

            - ``r`` -- an element of ``self.base_ring()``

            EXAMPLES::

                sage: A = AlgebrasWithBasis(QQ).example(); A                            # needs sage.combinat sage.modules
                An example of an algebra with basis:
                 the free algebra on the generators ('a', 'b', 'c') over Rational Field
                sage: A.from_base_ring(1)                                               # needs sage.combinat sage.modules
                B[word: ]
            """
            return self.one()._lmul_(r)

        def __init_extra__(self):
            """
            Declare the canonical coercion from ``self.base_ring()``
            to ``self``, if there has been none before.

            EXAMPLES::

                sage: A = AlgebrasWithBasis(QQ).example(); A                            # needs sage.combinat sage.modules
                An example of an algebra with basis:
                 the free algebra on the generators ('a', 'b', 'c') over Rational Field
                sage: coercion_model = sage.structure.element.get_coercion_model()
                sage: coercion_model.discover_coercion(QQ, A)                           # needs sage.combinat sage.modules
                ((map internal to coercion system -- copy before use)
                 Generic morphism:
                   From: Rational Field
                   To:   An example of an algebra with basis:
                          the free algebra on the generators ('a', 'b', 'c') over Rational Field,
                 None)
                sage: A(1)          # indirect doctest                                  # needs sage.combinat sage.modules
                B[word: ]

            TESTS:

            Ensure that :issue:`28328` is fixed and that non-associative
            algebras are supported::

                sage: # needs sage.modules
                sage: class Foo(CombinatorialFreeModule):
                ....:     def one(self):
                ....:         return self.monomial(0)
                sage: from sage.categories.magmatic_algebras import MagmaticAlgebras
                sage: C = MagmaticAlgebras(QQ).WithBasis().Unital()
                sage: F = Foo(QQ, (1,), category=C)
                sage: F(0)
                0
                sage: F(3)
                3*B[0]
            """
            base_ring = self.base_ring()
            if base_ring is self:
                # There are rings that are their own base rings. No need to register that.
                return
            if self._is_coercion_cached(base_ring):
                # We will not use any generic stuff, since a (presumably) better conversion
                # has already been registered.
                return
            if base_ring is None:
                # It may happen that self.base_ring() is not initialised at this point.
                return

            mor = self._coerce_map_from_base_ring()
            if mor is not None:
                mor._make_weak_references()
                try:
                    self.register_coercion(mor)
                except AssertionError:
                    pass

        def _coerce_map_from_(self, other):
            """
            Return a coercion map from ``other`` to ``self``, or ``None``.

            TESTS:

            Check that :issue:`19225` is solved::

                sage: A = cartesian_product((QQ['z'],)); A
                The Cartesian product of (Univariate Polynomial Ring in z over Rational Field,)
                sage: A.coerce_map_from(ZZ)
                Composite map:
                  From: Integer Ring
                  To:   The Cartesian product of (Univariate Polynomial Ring in z over Rational Field,)
                  Defn:   Natural morphism:
                          From: Integer Ring
                          To:   Rational Field
                        then
                          Generic morphism:
                          From: Rational Field
                          To:   The Cartesian product of (Univariate Polynomial Ring in z over Rational Field,)
                sage: A(1)
                (1,)
            """
            if other is self.base_ring():
                return self._coerce_map_from_base_ring()
            else:
                return self._coerce_map_via([self.base_ring()], other)

        def _coerce_map_from_base_ring(self):
            """
            Return a suitable coercion map from the base ring of ``self``.

            TESTS::

                sage: A = cartesian_product((QQ['z'],)); A
                The Cartesian product of (Univariate Polynomial Ring in z over Rational Field,)
                sage: A.base_ring()
                Rational Field
                sage: A._coerce_map_from_base_ring()
                Generic morphism:
                From: Rational Field
                To:   The Cartesian product of (Univariate Polynomial Ring in z over Rational Field,)

            Check that :issue:`29312` is fixed::

                sage: F.<x,y,z> = FreeAlgebra(QQ, implementation='letterplace')         # needs sage.combinat sage.modules
                sage: F._coerce_map_from_base_ring()                                    # needs sage.combinat sage.modules
                Generic morphism:
                  From: Rational Field
                  To:   Free Associative Unital Algebra on 3 generators (x, y, z) over Rational Field
            """
            base_ring = self.base_ring()

            # Pick a homset for the morphism to live in...
            if self in Rings():
                # The algebra is associative, and thus a ring. The
                # base ring is also a ring. Everything is OK.
                H = Hom(base_ring, self, Rings())
            else:
                # If the algebra isn't associative, we would like to
                # use the category of unital magmatic algebras (which
                # are not necessarily associative) instead. But,
                # unfortunately, certain important rings like QQ
                # aren't in that category. As a result, we have to use
                # something weaker.
                cat = Magmas().Unital()
                cat = Category.join([cat, CommutativeAdditiveGroups()])
                cat = cat.Distributive()
                H = Hom(base_ring, self, cat)

            # We need to construct a coercion from the base ring to self.
            #
            # There is a generic method from_base_ring(), that just does
            # multiplication with the multiplicative unit. However, the
            # unit is constructed repeatedly, which is slow.
            # So, if the unit is available *now*, then we can create a
            # faster coercion map.
            #
            # This only applies for the generic from_base_ring() method.
            # If there is a specialised from_base_ring(), then it should
            # be used unconditionally.
            generic_from_base_ring = self.category().parent_class.from_base_ring
            from_base_ring = self.from_base_ring   # bound method
            if from_base_ring.__func__ != generic_from_base_ring:
                # Custom from_base_ring()
                use_from_base_ring = True
            elif isinstance(generic_from_base_ring, lazy_attribute):
                # If the category implements from_base_ring() as lazy
                # attribute, then we always use it.
                # This is for backwards compatibility, see Issue #25181
                use_from_base_ring = True
            else:
                try:
                    one = self.one()
                    use_from_base_ring = False
                except (NotImplementedError, AttributeError, TypeError):
                    # The unit is not available, yet. But there are cases
                    # in which it will be available later. So, we use
                    # the generic from_base_ring() after all.
                    use_from_base_ring = True

            mor = None
            if use_from_base_ring:
                mor = SetMorphism(function=from_base_ring, parent=H)
            else:
                # We have the multiplicative unit, so implement the
                # coercion from the base ring as multiplying with that.
                #
                # But first we check that it actually works. If not,
                # then the generic implementation of from_base_ring()
                # would fail as well so we don't use it.
                try:
                    if one._lmul_(base_ring.an_element()) is not None:
                        # There are cases in which lmul returns None,
                        # which means that it's not implemented.
                        # One example: Hecke algebras.
                        mor = SetMorphism(function=one._lmul_, parent=H)
                except (NotImplementedError, AttributeError, TypeError):
                    pass
            return mor

    class WithBasis(CategoryWithAxiom_over_base_ring):

        class ParentMethods:

            @abstract_method(optional=True)
            def one_basis(self):
                """
                When the one of an algebra with basis is an element of
                this basis, this optional method can return the index of
                this element. This is used to provide a default
                implementation of :meth:`.one`, and an optimized default
                implementation of :meth:`.from_base_ring`.

                EXAMPLES::

                    sage: # needs sage.combinat sage.modules
                    sage: A = AlgebrasWithBasis(QQ).example()
                    sage: A.one_basis()
                    word:
                    sage: A.one()
                    B[word: ]
                    sage: A.from_base_ring(4)
                    4*B[word: ]
                """

            @cached_method
            def one_from_one_basis(self):
                """
                Return the one of the algebra, as per
                :meth:`Monoids.ParentMethods.one()
                <sage.categories.monoids.Monoids.ParentMethods.one>`

                By default, this is implemented from
                :meth:`.one_basis`, if available.

                EXAMPLES::

                    sage: # needs sage.combinat sage.modules
                    sage: A = AlgebrasWithBasis(QQ).example()
                    sage: A.one_basis()
                    word:
                    sage: A.one_from_one_basis()
                    B[word: ]
                    sage: A.one()
                    B[word: ]

                TESTS:

                Try to check that :issue:`5843` Heisenbug is fixed::

                    sage: # needs sage.combinat sage.modules
                    sage: A = AlgebrasWithBasis(QQ).example()
                    sage: B = AlgebrasWithBasis(QQ).example(('a', 'c'))
                    sage: A == B
                    False
                    sage: Aone = A.one_from_one_basis
                    sage: Bone = B.one_from_one_basis
                    sage: Aone is Bone
                    False

                Even if called in the wrong order, they should returns their
                respective one::

                    sage: Bone().parent() is B                                          # needs sage.combinat sage.modules
                    True
                    sage: Aone().parent() is A                                          # needs sage.combinat sage.modules
                    True
                """
                return self.monomial(self.one_basis())

            @lazy_attribute
            def one(self):
                r"""
                Return the multiplicative unit element.

                EXAMPLES::

                    sage: A = AlgebrasWithBasis(QQ).example()                           # needs sage.combinat sage.modules
                    sage: A.one_basis()                                                 # needs sage.combinat sage.modules
                    word:
                    sage: A.one()                                                       # needs sage.combinat sage.modules
                    B[word: ]
                """
                if self.one_basis is NotImplemented:
                    return NotImplemented
                return self.one_from_one_basis

            @lazy_attribute
            def from_base_ring(self):
                """
                TESTS::

                    sage: A = AlgebrasWithBasis(QQ).example()                           # needs sage.combinat sage.modules
                    sage: A.from_base_ring(3)                                           # needs sage.combinat sage.modules
                    3*B[word: ]
                """
                if self.one_basis is NotImplemented:
                    return NotImplemented
                return self.from_base_ring_from_one_basis

            def from_base_ring_from_one_basis(self, r):
                """
                Implement the canonical embedding from the ground ring.

                INPUT:

                - ``r`` -- an element of the coefficient ring

                EXAMPLES::

                    sage: # needs sage.combinat sage.modules
                    sage: A = AlgebrasWithBasis(QQ).example()
                    sage: A.from_base_ring_from_one_basis(3)
                    3*B[word: ]
                    sage: A.from_base_ring(3)
                    3*B[word: ]
                    sage: A(3)
                    3*B[word: ]
                """
                return self.term(self.one_basis(), r)

    class CartesianProducts(CartesianProductsCategory):
        r"""
        The category of unital algebras constructed as Cartesian products
        of unital algebras.

        This construction gives the direct product of algebras. See
        discussion on:

         - http://groups.google.fr/group/sage-devel/browse_thread/thread/35a72b1d0a2fc77a/348f42ae77a66d16#348f42ae77a66d16
         - :wikipedia:`Direct_product`
        """
        def extra_super_categories(self):
            """
            A Cartesian product of algebras is endowed with a natural
            unital algebra structure.

            EXAMPLES::

                sage: from sage.categories.unital_algebras import UnitalAlgebras
                sage: C = UnitalAlgebras(QQ).CartesianProducts()
                sage: C.extra_super_categories()
                [Category of unital algebras over Rational Field]
                sage: sorted(C.super_categories(), key=str)
                [Category of Cartesian products of distributive magmas and additive magmas,
                 Category of Cartesian products of unital magmas,
                 Category of Cartesian products of vector spaces over Rational Field,
                 Category of unital algebras over Rational Field]
            """
            return [self.base_category()]

        class ParentMethods:
            @cached_method
            def one(self):
                r"""
                Return the multiplicative unit element.

                EXAMPLES::

                    sage: # needs sage.graphs sage.modules
                    sage: S2 = simplicial_complexes.Sphere(2)
                    sage: H = S2.cohomology_ring(QQ)
                    sage: C = cartesian_product([H, H])
                    sage: one = C.one()
                    sage: one
                    B[(0, (0, 0))] + B[(1, (0, 0))]
                    sage: one == one * one
                    True
                    sage: all(b == b * one for b in C.basis())
                    True
                """
                data = enumerate(self.cartesian_factors())
                return self._cartesian_product_of_elements([C.one() for i, C in data])
