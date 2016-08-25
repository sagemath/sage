r"""
Unital algebras
"""
#*****************************************************************************
#  Copyright (C) 2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.morphism import SetMorphism
from sage.categories.homset import Hom
from sage.categories.rings import Rings
from sage.categories.magmas import Magmas
from sage.categories.magmatic_algebras import MagmaticAlgebras

class UnitalAlgebras(CategoryWithAxiom_over_base_ring):
    """
    The category of non-associative algebras over a given base ring.

    A non-associative algebra over a ring `R` is a module over `R`
    which s also a unital magma.

    .. WARNING::

        Until :trac:`15043` is implemented, :class:`Algebras` is the
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

                sage: A = AlgebrasWithBasis(QQ).example(); A
                An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field
                sage: A.from_base_ring(1)
                B[word: ]
            """
            return self.one()._lmul_(r)

        def __init_extra__(self):
            """
            Declare the canonical coercion from ``self.base_ring()``
            to ``self``, if there has been none before.

            EXAMPLES::

                sage: A = AlgebrasWithBasis(QQ).example(); A
                An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field
                sage: coercion_model = sage.structure.element.get_coercion_model()
                sage: coercion_model.discover_coercion(QQ, A)
                (Generic morphism:
                  From: Rational Field
                  To:   An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field, None)
                sage: A(1)          # indirect doctest
                B[word: ]

            """
            # If self has an attribute _no_generic_basering_coercion
            # set to True, then this declaration is skipped.
            # This trick, introduced in #11900, is used in
            # sage.matrix.matrix_space.py and
            # sage.rings.polynomial.polynomial_ring.
            # It will hopefully be refactored into something more
            # conceptual later on.
            if getattr(self, '_no_generic_basering_coercion', False):
                return

            base_ring = self.base_ring()
            if base_ring is self:
                # There are rings that are their own base rings. No need to register that.
                return
            if self.is_coercion_cached(base_ring):
                # We will not use any generic stuff, since a (presumably) better conversion
                # has already been registered.
                return
            mor = None
            # This could be a morphism of Algebras(self.base_ring()); however, e.g., QQ is not in Algebras(QQ)
            H = Hom(base_ring, self, Rings()) # TODO: non associative ring!

            # Idea: There is a generic method "from_base_ring", that just does multiplication with
            # the multiplicative unit. However, the unit is constructed repeatedly, which is slow.
            # Hence, if the unit is available *now*, then we store it.
            #
            # However, if there is a specialised from_base_ring method, then it should be used!
            try:
                has_custom_conversion = self.category().parent_class.from_base_ring.__func__ is not self.from_base_ring.__func__
            except AttributeError:
                # Sometimes from_base_ring is a lazy attribute
                has_custom_conversion = True
            if has_custom_conversion:
                mor = SetMorphism(function = self.from_base_ring, parent = H)
                try:
                    self.register_coercion(mor)
                except AssertionError:
                    pass
                return

            try:
                one = self.one()
            except (NotImplementedError, AttributeError, TypeError):
                # The unit is not available, yet. But there are cases
                # in which it will be available later. Hence:
                mor = SetMorphism(function = self.from_base_ring, parent = H)
            # try sanity of one._lmul_
            if mor is None:
                try:
                    if one._lmul_(base_ring.an_element()) is None:
                        # There are cases in which lmul returns None, believe it or not.
                        # One example: Hecke algebras.
                        # In that case, the generic implementation of from_base_ring would
                        # fail as well. Hence, unless it is overruled, we will not use it.
                        #mor = SetMorphism(function = self.from_base_ring, parent = H)
                        return
                except (NotImplementedError, AttributeError, TypeError):
                    # it is possible that an_element or lmul are not implemented.
                    return
                    #mor = SetMorphism(function = self.from_base_ring, parent = H)
                mor = SetMorphism(function = one._lmul_, parent = H)
            try:
                self.register_coercion(mor)
            except AssertionError:
                pass

    class ElementMethods:
        """
        ``Magmas.Element.__mul__`` is preferable to ``Modules.Element.__mul__``
        since the later does not handle products of two elements of ``self``.

        TESTS::

            sage: A = AlgebrasWithBasis(QQ).example(); A
            An example of an algebra with basis: the free algebra
            on the generators ('a', 'b', 'c') over Rational Field
            sage: x = A.an_element()
            sage: x
            B[word: ] + 2*B[word: a] + 3*B[word: b] + B[word: bab]
            sage: x.__mul__(x)
            B[word: ] + 4*B[word: a] + 4*B[word: aa] + 6*B[word: ab]
            + 2*B[word: abab] + 6*B[word: b] + 6*B[word: ba]
            + 2*B[word: bab] + 2*B[word: baba] + 3*B[word: babb]
            + B[word: babbab] + 9*B[word: bb] + 3*B[word: bbab]
        """
        __mul__ = Magmas.ElementMethods.__mul__.__func__

    class WithBasis(CategoryWithAxiom_over_base_ring):

        class ParentMethods:

            @abstract_method(optional = True)
            def one_basis(self):
                """
                When the one of an algebra with basis is an element of
                this basis, this optional method can return the index of
                this element. This is used to provide a default
                implementation of :meth:`.one`, and an optimized default
                implementation of :meth:`.from_base_ring`.

                EXAMPLES::

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

                    sage: A = AlgebrasWithBasis(QQ).example()
                    sage: A.one_basis()
                    word: 
                    sage: A.one_from_one_basis()
                    B[word: ]
                    sage: A.one()
                    B[word: ]

                TESTS:

                Try to check that :trac:`5843` Heisenbug is fixed::

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

                    sage: Bone().parent() is B
                    True
                    sage: Aone().parent() is A
                    True
                """
                return self.monomial(self.one_basis()) #.

            @lazy_attribute
            def one(self):
                r"""
                Return the multiplicative unit element.

                EXAMPLES::

                    sage: A = AlgebrasWithBasis(QQ).example()
                    sage: A.one_basis()
                    word: 
                    sage: A.one()
                    B[word: ]
                """
                if self.one_basis is NotImplemented:
                    return NotImplemented
                return self.one_from_one_basis

            @lazy_attribute
            def from_base_ring(self):
                """
                TESTS::

                    sage: A = AlgebrasWithBasis(QQ).example()
                    sage: A.from_base_ring(3)
                    3*B[word: ]
                """
                if self.one_basis is NotImplemented:
                    return NotImplemented
                return self.from_base_ring_from_one_basis

            def from_base_ring_from_one_basis(self, r):
                """
                Implement the canonical embeding from the ground ring.

                INPUT:

                - ``r`` -- an element of the coefficient ring

                EXAMPLES::

                    sage: A = AlgebrasWithBasis(QQ).example()
                    sage: A.from_base_ring_from_one_basis(3)
                    3*B[word: ]
                    sage: A.from_base_ring(3)
                    3*B[word: ]
                    sage: A(3)
                    3*B[word: ]
                """
                return self.term(self.one_basis(), r) #.

