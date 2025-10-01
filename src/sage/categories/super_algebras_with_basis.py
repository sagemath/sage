# sage_setup: distribution = sagemath-categories
r"""
Super algebras with basis
"""
# ****************************************************************************
#  Copyright (C) 2015,2019 Travis Scrimshaw <tcscrims at gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.categories.super_modules import SuperModulesCategory
from sage.categories.signed_tensor import SignedTensorProductsCategory
from sage.misc.cachefunc import cached_method


class SuperAlgebrasWithBasis(SuperModulesCategory):
    """
    The category of super algebras with a distinguished basis.

    EXAMPLES::

        sage: C = Algebras(ZZ).WithBasis().Super(); C
        Category of super algebras with basis over Integer Ring

    TESTS::

        sage: TestSuite(C).run()
    """
    def extra_super_categories(self):
        """
        EXAMPLES::

            sage: C = Algebras(ZZ).WithBasis().Super()
            sage: sorted(C.super_categories(), key=str) # indirect doctest
            [Category of graded algebras with basis over Integer Ring,
             Category of super algebras over Integer Ring,
             Category of super modules with basis over Integer Ring]
        """
        return [self.base_category().Graded()]

    class ParentMethods:
        def graded_algebra(self):
            r"""
            Return the associated graded module to ``self``.

            See :class:`~sage.algebras.associated_graded.AssociatedGradedAlgebra`
            for the definition and the properties of this.

            .. SEEALSO::

                :meth:`~sage.categories.filtered_modules_with_basis.ParentMethods.graded_algebra`

            EXAMPLES::

                sage: W.<x,y> = algebras.DifferentialWeyl(QQ)                           # needs sage.modules
                sage: W.graded_algebra()                                                # needs sage.combinat sage.modules
                Graded Algebra of Differential Weyl algebra of
                 polynomials in x, y over Rational Field
            """
            from sage.algebras.associated_graded import AssociatedGradedAlgebra
            return AssociatedGradedAlgebra(self)

    class ElementMethods:
        def supercommutator(self, x):
            r"""
            Return the supercommutator of ``self`` and ``x``.

            Let `A` be a superalgebra. The *supercommutator* of homogeneous
            elements `x, y \in A` is defined by

            .. MATH::

                [x, y\} = x y - (-1)^{|x| |y|} y x

            and extended to all elements by linearity.

            EXAMPLES::

                sage: # needs sage.modules
                sage: Q = QuadraticForm(ZZ, 3, [1,2,3,4,5,6])
                sage: Cl.<x,y,z> = CliffordAlgebra(Q)
                sage: a = x*y - z
                sage: b = x - y + y*z
                sage: a.supercommutator(b)
                -5*x*y + 8*x*z - 2*y*z - 6*x + 12*y - 5*z
                sage: a.supercommutator(Cl.one())
                0
                sage: Cl.one().supercommutator(a)
                0
                sage: Cl.zero().supercommutator(a)
                0
                sage: a.supercommutator(Cl.zero())
                0

                sage: # needs sage.modules
                sage: Q = QuadraticForm(ZZ, 2, [-1,1,-3])
                sage: Cl.<x,y> = CliffordAlgebra(Q)
                sage: [a.supercommutator(b) for a in Cl.basis() for b in Cl.basis()]
                [0, 0, 0, 0, 0, -2, 1, -x - 2*y, 0, 1,
                 -6, 6*x + y, 0, x + 2*y, -6*x - y, 0]
                sage: [a*b-b*a for a in Cl.basis() for b in Cl.basis()]
                [0, 0, 0, 0, 0, 0, 2*x*y - 1, -x - 2*y, 0,
                 -2*x*y + 1, 0, 6*x + y, 0, x + 2*y, -6*x - y, 0]

            Exterior algebras inherit from Clifford algebras, so
            supercommutators work as well. We verify the exterior algebra
            is supercommutative::

                sage: E.<x,y,z,w> = ExteriorAlgebra(QQ)                                 # needs sage.modules
                sage: all(b1.supercommutator(b2) == 0                                   # needs sage.modules
                ....:     for b1 in E.basis() for b2 in E.basis())
                True
            """
            P = self.parent()
            ret = P.zero()
            for ms, cs in self:
                term_s = P.term(ms, cs)
                sign_s = (-1)**P.degree_on_basis(ms)
                for mx, cx in x:
                    ret += term_s * P.term(mx, cx)
                    s = sign_s**P.degree_on_basis(mx)
                    ret -= s * P.term(mx, cx) * term_s
            return ret

    class SignedTensorProducts(SignedTensorProductsCategory):
        """
        The category of super algebras with basis constructed by tensor
        product of super algebras with basis.
        """
        @cached_method
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: Algebras(QQ).Super().SignedTensorProducts().extra_super_categories()
                [Category of super algebras over Rational Field]
                sage: Algebras(QQ).Super().SignedTensorProducts().super_categories()
                [Category of signed tensor products of graded algebras over Rational Field,
                 Category of super algebras over Rational Field]

            Meaning: a signed tensor product of super algebras is a super algebra
            """
            return [self.base_category()]
