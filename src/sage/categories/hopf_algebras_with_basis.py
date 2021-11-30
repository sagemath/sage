r"""
Hopf algebras with basis
"""
# ****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#  Copyright (C) 2008-2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.tensor import TensorProductsCategory
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.lazy_import import LazyImport


class HopfAlgebrasWithBasis(CategoryWithAxiom_over_base_ring):
    """
    The category of Hopf algebras with a distinguished basis

    EXAMPLES::

        sage: C = HopfAlgebrasWithBasis(QQ)
        sage: C
        Category of hopf algebras with basis over Rational Field
        sage: C.super_categories()
        [Category of hopf algebras over Rational Field,
         Category of bialgebras with basis over Rational Field]

    We now show how to use a simple Hopf algebra, namely the group algebra of the dihedral group
    (see also AlgebrasWithBasis)::

        sage: A = C.example(); A
        An example of Hopf algebra with basis: the group algebra of the Dihedral group of order 6 as a permutation group over Rational Field
        sage: A.__custom_name = "A"
        sage: A.category()
        Category of finite dimensional hopf algebras with basis over Rational Field

        sage: A.one_basis()
        ()
        sage: A.one()
        B[()]

        sage: A.base_ring()
        Rational Field
        sage: A.basis().keys()
        Dihedral group of order 6 as a permutation group

        sage: [a,b] = A.algebra_generators()
        sage: a, b
        (B[(1,2,3)], B[(1,3)])
        sage: a^3, b^2
        (B[()], B[()])
        sage: a*b
        B[(1,2)]

        sage: A.product           # todo: not quite ...
        <bound method MagmaticAlgebras.WithBasis.ParentMethods._product_from_product_on_basis_multiply of A>
        sage: A.product(b,b)
        B[()]

        sage: A.zero().coproduct()
        0
        sage: A.zero().coproduct().parent()
        A # A
        sage: a.coproduct()
        B[(1,2,3)] # B[(1,2,3)]

        sage: TestSuite(A).run(verbose=True)
        running ._test_additive_associativity() . . . pass
        running ._test_an_element() . . . pass
        running ._test_antipode() . . . pass
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
        <class 'sage.categories.examples.hopf_algebras_with_basis.MyGroupAlgebra_with_category'>
        sage: A.element_class
        <class 'sage.categories.examples.hopf_algebras_with_basis.MyGroupAlgebra_with_category.element_class'>

    Let us look at the code for implementing A::

        sage: A??                       # todo: not implemented

    TESTS::

        sage: TestSuite(A).run()
        sage: TestSuite(C).run()
    """

    def example(self, G = None):
        """
        Returns an example of algebra with basis::

            sage: HopfAlgebrasWithBasis(QQ['x']).example()
            An example of Hopf algebra with basis: the group algebra of the Dihedral group of order 6 as a permutation group over Univariate Polynomial Ring in x over Rational Field

        An other group can be specified as optional argument::

            sage: HopfAlgebrasWithBasis(QQ).example(SymmetricGroup(4))
            An example of Hopf algebra with basis: the group algebra of the Symmetric group of order 4! as a permutation group over Rational Field
        """
        from sage.categories.examples.hopf_algebras_with_basis import MyGroupAlgebra
        from sage.groups.perm_gps.permgroup_named import DihedralGroup
        if G is None:
            G = DihedralGroup(3)
        return MyGroupAlgebra(self.base_ring(), G)

#     This is only correct in the finite dimensional / graded case
#     def dual(self):
#         """
#         Returns the dual category

#         EXAMPLES:

#         The category of Hopf algebras over any field is self dual::

#             sage: C = HopfAlgebrasWithBasis(QQ)
#             sage: C.dual()
#             Category of hopf algebras with basis over Rational Field
#         """
#         return self

    FiniteDimensional = LazyImport('sage.categories.finite_dimensional_hopf_algebras_with_basis',
                                   'FiniteDimensionalHopfAlgebrasWithBasis')
    Filtered = LazyImport('sage.categories.filtered_hopf_algebras_with_basis',
                          'FilteredHopfAlgebrasWithBasis')
    Graded = LazyImport('sage.categories.graded_hopf_algebras_with_basis',
                        'GradedHopfAlgebrasWithBasis')
    Super = LazyImport('sage.categories.super_hopf_algebras_with_basis',
                       'SuperHopfAlgebrasWithBasis')

    class ParentMethods:

        @abstract_method(optional=True)
        def antipode_on_basis(self, x):
            """
            The antipode of the Hopf algebra on the basis (optional)

            INPUT:

             - ``x`` -- an index of an element of the basis of ``self``

            Returns the antipode of the basis element indexed by ``x``.

            If this method is implemented, then :meth:`antipode` is defined
            from this by linearity.

            EXAMPLES::

                sage: A = HopfAlgebrasWithBasis(QQ).example()
                sage: W = A.basis().keys(); W
                Dihedral group of order 6 as a permutation group
                sage: w = W.gen(0); w
                (1,2,3)
                sage: A.antipode_on_basis(w)
                B[(1,3,2)]
            """

        @lazy_attribute
        def antipode(self):
            """
            The antipode of this Hopf algebra.

            If :meth:`.antipode_basis` is available, this constructs the
            antipode morphism from ``self`` to ``self`` by extending it by
            linearity. Otherwise, :meth:`self.antipode_by_coercion` is used, if
            available.

            EXAMPLES::

                sage: A = HopfAlgebrasWithBasis(ZZ).example(); A
                An example of Hopf algebra with basis: the group algebra of the Dihedral group of order 6 as a permutation group over Integer Ring
                sage: A = HopfAlgebrasWithBasis(QQ).example()
                sage: [a,b] = A.algebra_generators()
                sage: a, A.antipode(a)
                (B[(1,2,3)], B[(1,3,2)])
                sage: b, A.antipode(b)
                (B[(1,3)], B[(1,3)])

            TESTS::

                sage: all(A.antipode(x) * x == A.one() for x in A.basis())
                True
            """
            if self.antipode_on_basis is not NotImplemented:
                # Should give the information that this is an anti-morphism of algebra
                return self._module_morphism(self.antipode_on_basis, codomain = self)
            elif hasattr(self, "antipode_by_coercion"):
                return self.antipode_by_coercion

        def _test_antipode(self, **options):
            r"""
            Test the antipode.

            An *antipode* `S` of a Hopf algebra is a linear endomorphism of the
            Hopf algebra that satisfies the following conditions (see
            :wikipedia:`HopfAlgebra`).

            - If `\mu` and `\Delta` denote the product and coproduct of the
              Hopf algebra, respectively, then `S` satisfies

              .. MATH::

                  \mu \circ (S \tensor 1) \circ \Delta = unit \circ counit
                  \mu \circ (1 \tensor S) \circ \Delta = unit \circ counit

            - `S` is an *anti*-homomorphism

            These properties are tested on :meth:`some_elements`.

            TESTS::

                sage: R = NonCommutativeSymmetricFunctions(QQ).ribbon()
                sage: R._test_antipode()

            ::

                sage: s = SymmetricFunctions(QQ).schur()
                sage: s._test_antipode()

            """
            tester = self._tester(**options)

            S = self.antipode

            IS = lambda x: self.sum(c * self.monomial(t1) * S(self.monomial(t2))
                                for ((t1, t2), c) in x.coproduct())

            SI = lambda x: self.sum(c * S(self.monomial(t1)) * self.monomial(t2)
                                for ((t1, t2), c) in x.coproduct())

            for x in tester.some_elements():

                # antipode is an anti-homomorphism
                for y in tester.some_elements():
                    tester.assertEqual(S(x) * S(y), S(y * x))

                # mu * (S # I) * delta == counit * unit
                tester.assertEqual(SI(x), self.counit(x) * self.one())

                # mu * (I # S) * delta == counit * unit
                tester.assertEqual(IS(x), self.counit(x) * self.one())

    class ElementMethods:
        pass

    class TensorProducts(TensorProductsCategory):
        """
        The category of hopf algebras with basis constructed by tensor product of hopf algebras with basis
        """

        @cached_method
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: C = HopfAlgebrasWithBasis(QQ).TensorProducts()
                sage: C.extra_super_categories()
                [Category of hopf algebras with basis over Rational Field]
                sage: sorted(C.super_categories(), key=str)
                [Category of hopf algebras with basis over Rational Field,
                 Category of tensor products of algebras with basis over Rational Field,
                 Category of tensor products of hopf algebras over Rational Field]
            """
            return [self.base_category()]

        class ParentMethods:
            # todo: antipode
            pass

        class ElementMethods:
            pass

