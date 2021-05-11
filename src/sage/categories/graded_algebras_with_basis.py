r"""
Graded algebras with basis
"""
#*****************************************************************************
#  Copyright (C) 2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.graded_modules import GradedModulesCategory
from sage.categories.signed_tensor import SignedTensorProductsCategory, tensor_signed
from sage.misc.cachefunc import cached_method

class GradedAlgebrasWithBasis(GradedModulesCategory):
    """
    The category of graded algebras with a distinguished basis

    EXAMPLES::

        sage: C = GradedAlgebrasWithBasis(ZZ); C
        Category of graded algebras with basis over Integer Ring
        sage: sorted(C.super_categories(), key=str)
        [Category of filtered algebras with basis over Integer Ring,
         Category of graded algebras over Integer Ring,
         Category of graded modules with basis over Integer Ring]

    TESTS::

        sage: TestSuite(C).run()
    """
    class ParentMethods:
        # This needs to be copied in GradedAlgebras because we need to have
        #   FilteredAlgebrasWithBasis as an extra super category
        def graded_algebra(self):
            """
            Return the associated graded algebra to ``self``.

            This is ``self``, because ``self`` is already graded.
            See :meth:`~sage.categories.filtered_algebras_with_basis.FilteredAlgebrasWithBasis.graded_algebra`
            for the general behavior of this method, and see
            :class:`~sage.algebras.associated_graded.AssociatedGradedAlgebra`
            for the definition and properties of associated graded
            algebras.

            EXAMPLES::

                sage: m = SymmetricFunctions(QQ).m()
                sage: m.graded_algebra() is m
                True

            TESTS:

            Let us check that the three methods
            :meth:`to_graded_conversion`, :meth:`from_graded_conversion`
            and :meth:`projection` (which form the interface of the
            associated graded algebra) work correctly here::

                sage: to_gr = m.to_graded_conversion()
                sage: from_gr = m.from_graded_conversion()
                sage: m[2] == to_gr(m[2]) == from_gr(m[2])
                True
                sage: u = 3*m[1] - (1/2)*m[3]
                sage: u == to_gr(u) == from_gr(u)
                True
                sage: m.zero() == to_gr(m.zero()) == from_gr(m.zero())
                True
                sage: p2 = m.projection(2)
                sage: p2(m[2] - 4*m[1,1] + 3*m[1] - 2*m[[]])
                -4*m[1, 1] + m[2]
                sage: p2(4*m[1])
                0
                sage: p2(m.zero()) == m.zero()
                True
            """
            return self

        # .. TODO::
        # Possibly override ``to_graded_conversion`` and
        # ``from_graded_conversion`` with identity morphisms?
        # I have to admit I don't know the right way to construct
        # identity morphisms other than using the identity matrix.
        # Also, ``projection`` could be overridden by, well, a
        # projection.

    class ElementMethods:
        pass

    class SignedTensorProducts(SignedTensorProductsCategory):
        """
        The category of algebras with basis constructed by signed tensor
        product of algebras with basis.
        """
        @cached_method
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: Cat = AlgebrasWithBasis(QQ).Graded()
                sage: Cat.SignedTensorProducts().extra_super_categories()
                [Category of graded algebras with basis over Rational Field]
                sage: Cat.SignedTensorProducts().super_categories()
                [Category of graded algebras with basis over Rational Field,
                 Category of signed tensor products of graded algebras over Rational Field]
            """
            return [self.base_category()]

        class ParentMethods:
            """
            Implements operations on tensor products of super algebras
            with basis.
            """
            @cached_method
            def one_basis(self):
                """
                Return the index of the one of this signed tensor product of
                algebras, as per ``AlgebrasWithBasis.ParentMethods.one_basis``.

                It is the tuple whose operands are the indices of the
                ones of the operands, as returned by their
                :meth:`.one_basis` methods.

                EXAMPLES::

                    sage: A.<x,y> = ExteriorAlgebra(QQ)
                    sage: A.one_basis()
                    ()
                    sage: B = tensor((A, A, A))
                    sage: B.one_basis()
                    ((), (), ())
                    sage: B.one()
                    1 # 1 # 1
                """
                # FIXME: this method should be conditionally defined,
                # so that B.one_basis returns NotImplemented if not
                # all modules provide one_basis
                if all(hasattr(module, "one_basis") for module in self._sets):
                    return tuple(module.one_basis() for module in self._sets)
                else:
                    raise NotImplementedError

            def product_on_basis(self, t0, t1):
                """
                The product of the algebra on the basis, as per
                ``AlgebrasWithBasis.ParentMethods.product_on_basis``.

                EXAMPLES:

                Test the sign in the super tensor product::

                    sage: A = SteenrodAlgebra(3)
                    sage: x = A.Q(0)
                    sage: y = x.coproduct()
                    sage: y^2
                    0

                TODO: optimize this implementation!
                """
                basic = tensor_signed((module.monomial(x0) * module.monomial(x1)
                                      for (module, x0, x1) in zip(self._sets, t0, t1)))
                n = len(self._sets)
                parity0 = [self._sets[idx].degree_on_basis(x0)
                           for (idx, x0) in enumerate(t0)]
                parity1 = [self._sets[idx].degree_on_basis(x1)
                           for (idx, x1) in enumerate(t1)]
                parity = sum(parity0[i] * parity1[j]
                             for j in range(n) for i in range(j+1,n))
                return (-1)**parity * basic

