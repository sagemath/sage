r"""
Examples of parents endowed with multiple realizations
"""
#*****************************************************************************
#  Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.bindable_class import BindableClass
from sage.categories.all import Rings, Algebras, AlgebrasWithBasis
from sage.categories.realizations import Category_realization_of_parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.sets.set import Set
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.subset import Subsets

class SubsetAlgebra(UniqueRepresentation, Parent):
    r"""
    An example of parent endowed with several realizations

    We consider an algebra `A(S)` whose bases are indexed by the
    subsets `s` of a given set `S`. We consider three natural basis of
    this algebra: ``F``, ``In``, and ``Out``. In the first basis, the
    product is given by the union of the indexing sets. That is, for any
    `s, t\subset S`

    .. MATH::

        F_s F_t  = F_{s\cup t}

    The ``In`` basis and ``Out`` basis are defined respectively by:

    .. MATH::

        In_s  = \sum_{t\subset s} F_t
        \qquad\text{and}\qquad
        F_s   = \sum_{t\supset s} Out_t

    Each such basis gives a realization of `A`, where the elements are
    represented by their expansion in this basis.

    This parent, and its code, demonstrate how to implement this
    algebra and its three realizations, with coercions and mixed
    arithmetic between them.

    .. SEEALSO::

       - :func:`Sets().WithRealizations <sage.categories.with_realizations.WithRealizations>`
       - the `Implementing Algebraic Structures
         <../../../../../thematic_tutorials/tutorial-implementing-algebraic-structures>`_
         thematic tutorial.

    EXAMPLES::

        sage: A = Sets().WithRealizations().example(); A
        The subset algebra of {1, 2, 3} over Rational Field
        sage: A.base_ring()
        Rational Field

    The three bases of ``A``::

        sage: F   = A.F()  ; F
        The subset algebra of {1, 2, 3} over Rational Field in the Fundamental basis
        sage: In  = A.In() ; In
        The subset algebra of {1, 2, 3} over Rational Field in the In basis
        sage: Out = A.Out(); Out
        The subset algebra of {1, 2, 3} over Rational Field in the Out basis

    One can quickly define all the bases using the following shortcut::

        sage: A.inject_shorthands()
        Defining F as shorthand for The subset algebra of {1, 2, 3} over Rational Field in the Fundamental basis
        Defining In as shorthand for The subset algebra of {1, 2, 3} over Rational Field in the In basis
        Defining Out as shorthand for The subset algebra of {1, 2, 3} over Rational Field in the Out basis

    Accessing the basis elements is done with :meth:`basis()` method::

        sage: F.basis().list()
        [F[{}], F[{1}], F[{2}], F[{3}], F[{1, 2}], F[{1, 3}], F[{2, 3}], F[{1, 2, 3}]]

    To access a particular basis element, you can use the :meth:`from_set`
    method::

        sage: F.from_set(2,3)
        F[{2, 3}]
        sage: In.from_set(1,3)
        In[{1, 3}]

    or as a convenient shorthand, one can use the following notation::

        sage: F[2,3]
        F[{2, 3}]
        sage: In[1,3]
        In[{1, 3}]

    Some conversions::

        sage: F(In[2,3])
        F[{}] + F[{2}] + F[{3}] + F[{2, 3}]
        sage: In(F[2,3])
        In[{}] - In[{2}] - In[{3}] + In[{2, 3}]

        sage: Out(F[3])
        Out[{3}] + Out[{1, 3}] + Out[{2, 3}] + Out[{1, 2, 3}]
        sage: F(Out[3])
        F[{3}] - F[{1, 3}] - F[{2, 3}] + F[{1, 2, 3}]

        sage: Out(In[2,3])
        Out[{}] + Out[{1}] + 2*Out[{2}] + 2*Out[{3}] + 2*Out[{1, 2}] + 2*Out[{1, 3}] + 4*Out[{2, 3}] + 4*Out[{1, 2, 3}]

    We can now mix expressions::

        sage: (1 + Out[1]) * In[2,3]
        Out[{}] + 2*Out[{1}] + 2*Out[{2}] + 2*Out[{3}] + 2*Out[{1, 2}] + 2*Out[{1, 3}] + 4*Out[{2, 3}] + 4*Out[{1, 2, 3}]

    """

    def __init__(self, R, S):
        r"""
        EXAMPLES::

            sage: from sage.categories.examples.with_realizations import SubsetAlgebra
            sage: A = SubsetAlgebra(QQ, Set((1,2,3))); A
            The subset algebra of {1, 2, 3} over Rational Field
            sage: Sets().WithRealizations().example() is A
            True
            sage: TestSuite(A).run()

        TESTS::

            sage: A = Sets().WithRealizations().example(); A
            The subset algebra of {1, 2, 3} over Rational Field
            sage: F, In, Out = A.realizations()
            sage: type(F.coerce_map_from(In))
            <class 'sage.modules.with_basis.morphism.TriangularModuleMorphismByLinearity_with_category'>
            sage: type(In.coerce_map_from(F))
            <class 'sage.modules.with_basis.morphism.TriangularModuleMorphismByLinearity_with_category'>
            sage: type(F.coerce_map_from(Out))
            <class 'sage.modules.with_basis.morphism.TriangularModuleMorphismByLinearity_with_category'>
            sage: type(Out.coerce_map_from(F))
            <class 'sage.modules.with_basis.morphism.TriangularModuleMorphismByLinearity_with_category'>
            sage: In.coerce_map_from(Out)
            Composite map:
              From: The subset algebra of {1, 2, 3} over Rational Field in the Out basis
              To:   The subset algebra of {1, 2, 3} over Rational Field in the In basis
              Defn:   Generic morphism:
                      From: The subset algebra of {1, 2, 3} over Rational Field in the Out basis
                      To:   The subset algebra of {1, 2, 3} over Rational Field in the Fundamental basis
                    then
                      Generic morphism:
                      From: The subset algebra of {1, 2, 3} over Rational Field in the Fundamental basis
                      To:   The subset algebra of {1, 2, 3} over Rational Field in the In basis
            sage: Out.coerce_map_from(In)
            Composite map:
              From: The subset algebra of {1, 2, 3} over Rational Field in the In basis
              To:   The subset algebra of {1, 2, 3} over Rational Field in the Out basis
              Defn:   Generic morphism:
                      From: The subset algebra of {1, 2, 3} over Rational Field in the In basis
                      To:   The subset algebra of {1, 2, 3} over Rational Field in the Fundamental basis
                    then
                      Generic morphism:
                      From: The subset algebra of {1, 2, 3} over Rational Field in the Fundamental basis
                      To:   The subset algebra of {1, 2, 3} over Rational Field in the Out basis
        """
        assert(R in Rings())
        self._base = R # Won't be needed when CategoryObject won't override anymore base_ring
        self._S = S
        Parent.__init__(self, category = Algebras(R).Commutative().WithRealizations())

        # Initializes the bases and change of bases of ``self``

        F = self.F()
        In = self.In()
        Out = self.Out()

        category = self.Bases()
        key = self.indices_key
        In_to_F = In.module_morphism(F.sum_of_monomials * Subsets,
                                     codomain=F, category=category,
                                     triangular='upper', unitriangular=True,
                                     key=key)
        In_to_F   .register_as_coercion()
        (~In_to_F).register_as_coercion()

        F_to_Out = F.module_morphism(Out.sum_of_monomials * self.supsets,
                                     codomain=Out, category=category,
                                     triangular='lower', unitriangular=True,
                                     key=key)
        F_to_Out   .register_as_coercion()
        (~F_to_Out).register_as_coercion()

    _shorthands = ["F", "In", "Out"]

    def a_realization(self):
        r"""
        Returns the default realization of ``self``

        EXAMPLES::

            sage: A = Sets().WithRealizations().example(); A
            The subset algebra of {1, 2, 3} over Rational Field
            sage: A.a_realization()
            The subset algebra of {1, 2, 3} over Rational Field in the Fundamental basis
        """
        return self.F()

    def base_set(self):
        r"""
        EXAMPLES::

            sage: A = Sets().WithRealizations().example(); A
            The subset algebra of {1, 2, 3} over Rational Field
            sage: A.base_set()
            {1, 2, 3}
        """
        return self._S

    def indices(self):
        r"""
        The objects that index the basis elements of this algebra.

        EXAMPLES::

            sage: A = Sets().WithRealizations().example(); A
            The subset algebra of {1, 2, 3} over Rational Field
            sage: A.indices()
            Subsets of {1, 2, 3}
        """
        return Subsets(self._S)

    def indices_key(self, x):
        r"""
        A key function on a set which gives a linear extension
        of the inclusion order.

        INPUT:

        - ``x`` -- set

        EXAMPLES::

            sage: A = Sets().WithRealizations().example(); A
            The subset algebra of {1, 2, 3} over Rational Field
            sage: sorted(A.indices(), key=A.indices_key)
            [{}, {1}, {2}, {3}, {1, 2}, {1, 3}, {2, 3}, {1, 2, 3}]
        """
        return (len(x), list(x))

    def supsets(self, set):
        r"""
        Returns all the subsets of `S` containing ``set``

        INPUT:

        - ``set`` -- a subset of the base set `S` of ``self``

        EXAMPLES::

            sage: A = Sets().WithRealizations().example(); A
            The subset algebra of {1, 2, 3} over Rational Field
            sage: A.supsets(Set((2,)))
            [{1, 2, 3}, {2, 3}, {1, 2}, {2}]
        """
        S = self.base_set()
        return list(S.difference(s) for s in Subsets(S.difference(set)))

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: Sets().WithRealizations().example()   # indirect doctest
            The subset algebra of {1, 2, 3} over Rational Field
        """
        return "The subset algebra of %s over %s"%(self.base_set(), self.base_ring())

    class Bases(Category_realization_of_parent):
        r"""
        The category of the realizations of the subset algebra
        """

        def super_categories(self):
            r"""
            EXAMPLES::

                sage: A = Sets().WithRealizations().example(); A
                The subset algebra of {1, 2, 3} over Rational Field
                sage: C = A.Bases(); C
                Category of bases of The subset algebra of {1, 2, 3} over Rational Field
                sage: C.super_categories()
                [Category of realizations of The subset algebra of {1, 2, 3} over Rational Field,
                 Join of Category of algebras with basis over Rational Field and
                         Category of commutative algebras over Rational Field and
                         Category of realizations of unital magmas]
            """
            A = self.base()
            category = Algebras(A.base_ring()).Commutative()
            return [A.Realizations(),
                    category.Realizations().WithBasis()]


        class ParentMethods:

            def from_set(self, *args):
                r"""
                Construct the monomial indexed by the set containing the
                elements passed as arguments.

                EXAMPLES::

                    sage: In = Sets().WithRealizations().example().In(); In
                    The subset algebra of {1, 2, 3} over Rational Field in the In basis
                    sage: In.from_set(2,3)
                    In[{2, 3}]

                As a shorthand, one can construct elements using the following
                notation::

                    sage: In[2,3]
                    In[{2, 3}]
                """
                return self.monomial(Set(args))

            def __getitem__(self, s):
                r"""
                This method implements a convenient shorthand for constructing
                basis elements of this algebra.

                EXAMPLES::

                    sage: A = Sets().WithRealizations().example()
                    sage: In = A.In()
                    sage: In[2,3]
                    In[{2, 3}]
                    sage: F = A.F()
                    sage: F[1,3]
                    F[{1, 3}]
                """
                from sage.rings.integer import Integer
                if isinstance(s, Integer):
                    return self.from_set(*(s,))
                else:
                    return self.from_set(*s)

            # This could go in the super category VectorSpaces().Realizations()
            def _repr_(self):
                r"""
                EXAMPLES::

                    sage: Sets().WithRealizations().example().In()  # indirect doctest
                    The subset algebra of {1, 2, 3} over Rational Field in the In basis
                """
                return "%s in the %s basis" % (self.realization_of(), self._realization_name())

            # Could this go in the super category Monoids().Realizations() ?
            @cached_method
            def one(self):
                r"""
                Returns the unit of this algebra.

                This default implementation takes the unit in the
                fundamental basis, and coerces it in ``self``.

                EXAMPLES::

                    sage: A = Sets().WithRealizations().example(); A
                    The subset algebra of {1, 2, 3} over Rational Field
                    sage: In = A.In(); Out = A.Out()
                    sage: In.one()
                    In[{}]
                    sage: Out.one()
                    Out[{}] + Out[{1}] + Out[{2}] + Out[{3}] + Out[{1, 2}] + Out[{1, 3}] + Out[{2, 3}] + Out[{1, 2, 3}]
                """
                return self(self.realization_of().F().one())

    class Fundamental(CombinatorialFreeModule, BindableClass):
        r"""
        The Subset algebra, in the fundamental basis

        INPUT:

        - ``A`` -- a parent with realization in :class:`SubsetAlgebra`

        EXAMPLES::

            sage: A = Sets().WithRealizations().example()
            sage: A.F()
            The subset algebra of {1, 2, 3} over Rational Field in the Fundamental basis
            sage: A.Fundamental()
            The subset algebra of {1, 2, 3} over Rational Field in the Fundamental basis
        """

        def __init__(self, A):
            r"""
            EXAMPLES::

                sage: A = Sets().WithRealizations().example(); A
                The subset algebra of {1, 2, 3} over Rational Field
                sage: F = A.F(); F
                The subset algebra of {1, 2, 3} over Rational Field in the Fundamental basis
                sage: TestSuite(F).run()
            """
            CombinatorialFreeModule.__init__(self,
                A.base_ring(), A.indices(),
                category=A.Bases(), prefix='F', sorting_key=A.indices_key)

        def product_on_basis(self, left, right):
            r"""
            Product of basis elements, as per :meth:`AlgebrasWithBasis.ParentMethods.product_on_basis`.

            INPUT:

            - ``left``, ``right`` -- sets indexing basis elements

            EXAMPLES::

                sage: F = Sets().WithRealizations().example().F(); F
                The subset algebra of {1, 2, 3} over Rational Field in the Fundamental basis
                sage: S = F.basis().keys(); S
                Subsets of {1, 2, 3}
                sage: F.product_on_basis(S([]), S([]))
                F[{}]
                sage: F.product_on_basis(S({1}), S({3}))
                F[{1, 3}]
                sage: F.product_on_basis(S({1,2}), S({2,3}))
                F[{1, 2, 3}]
            """
            return self.monomial(left.union(right))

        def one_basis(self):
            r"""
            Returns the index of the basis element which is equal to '1'.

            EXAMPLES::

                sage: F = Sets().WithRealizations().example().F(); F
                The subset algebra of {1, 2, 3} over Rational Field in the Fundamental basis
                sage: F.one_basis()
                {}
                sage: F.one()
                F[{}]
            """
            return Set([])

        # Bypass the definition in SubsetAlgebra.Bases.ParentMethods
        # by the usual one from basis; alternatively one could just
        # define one instead of one_basis.
        one = AlgebrasWithBasis.ParentMethods.one

    F = Fundamental

    class In(CombinatorialFreeModule, BindableClass):
        r"""
        The Subset Algebra, in the ``In`` basis

        INPUT:

        - ``A`` -- a parent with realization in :class:`SubsetAlgebra`

        EXAMPLES::

            sage: A = Sets().WithRealizations().example()
            sage: A.In()
            The subset algebra of {1, 2, 3} over Rational Field in the In basis

        TESTS:

        The product in this basis is computed by converting to the fundamental
        basis, computing the product there, and then converting back::

            sage: In = Sets().WithRealizations().example().In(); In
            The subset algebra of {1, 2, 3} over Rational Field in the In basis
            sage: x = In.an_element()
            sage: y = In.an_element()
            sage: In.product
            <bound method ....product_by_coercion ...>
            sage: In.product.__module__
            'sage.categories.magmas'
            sage: In.product(x, y)
            -21*In[{}] - 2*In[{1}] + 19*In[{2}] + 53*In[{1, 2}]
        """

        def __init__(self, A):
            r"""
            EXAMPLES::

                sage: In = Sets().WithRealizations().example().In(); In
                The subset algebra of {1, 2, 3} over Rational Field in the In basis
                sage: TestSuite(In).run()
            """
            CombinatorialFreeModule.__init__(self,
                A.base_ring(), A.indices(),
                category=A.Bases(), prefix='In', sorting_key=A.indices_key)

    class Out(CombinatorialFreeModule, BindableClass):
        r"""
        The Subset Algebra, in the `Out` basis

        INPUT:

        - ``A`` -- a parent with realization in :class:`SubsetAlgebra`

        EXAMPLES::

            sage: A = Sets().WithRealizations().example()
            sage: A.Out()
            The subset algebra of {1, 2, 3} over Rational Field in the Out basis

        TESTS:

        The product in this basis is computed by converting to the fundamental
        basis, computing the product there, and then converting back::

            sage: Out = Sets().WithRealizations().example().Out(); Out
            The subset algebra of {1, 2, 3} over Rational Field in the Out basis
            sage: x = Out.an_element()
            sage: y = Out.an_element()
            sage: Out.product
            <bound method ....product_by_coercion ...>
            sage: Out.product.__module__
            'sage.categories.magmas'
            sage: Out.product(x, y)
            Out[{}] + 4*Out[{1}] + 9*Out[{2}] + Out[{1, 2}]
        """

        def __init__(self, A):
            r"""
            EXAMPLES::

                sage: Out = Sets().WithRealizations().example().Out(); Out
                The subset algebra of {1, 2, 3} over Rational Field in the Out basis
                sage: TestSuite(Out).run()
            """
            CombinatorialFreeModule.__init__(self,
                A.base_ring(), A.indices(),
                category=A.Bases(), prefix='Out', sorting_key=A.indices_key)
