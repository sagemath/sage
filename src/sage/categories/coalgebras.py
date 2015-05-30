r"""
Coalgebras
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#  Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from category_types import Category_over_base_ring
from sage.categories.all import Modules
from sage.categories.tensor import TensorProductsCategory, tensor
from sage.categories.dual import DualObjectsCategory
from sage.categories.realizations import RealizationsCategory
from sage.categories.with_realizations import WithRealizationsCategory
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import LazyImport

class Coalgebras(Category_over_base_ring):
    """
    The category of coalgebras

    EXAMPLES::

        sage: Coalgebras(QQ)
        Category of coalgebras over Rational Field
        sage: Coalgebras(QQ).super_categories()
        [Category of vector spaces over Rational Field]

    TESTS::

        sage: TestSuite(Coalgebras(ZZ)).run()
    """
    def super_categories(self):
        """
        EXAMPLES::

            sage: Coalgebras(QQ).super_categories()
            [Category of vector spaces over Rational Field]
        """
        return [Modules(self.base_ring())]

    WithBasis = LazyImport('sage.categories.coalgebras_with_basis',  'CoalgebrasWithBasis')

    class ParentMethods:
        #def __init_add__(self): # The analogue of initDomainAdd
        #    # Will declare the coproduct of self to the coercion mechanism when it exists
        #    pass

        @cached_method
        def tensor_square(self):
            """
            Returns the tensor square of ``self``

            EXAMPLES::

                sage: A = HopfAlgebrasWithBasis(QQ).example()
                sage: A.tensor_square()
                An example of Hopf algebra with basis: the group algebra of the Dihedral group of order 6 as a permutation group over Rational Field # An example of Hopf algebra with basis: the group algebra of the Dihedral group of order 6 as a permutation group over Rational Field
            """
            return tensor([self, self])

        @abstract_method
        def counit(self, x):
            """
            Returns the counit of x.

            Eventually, there will be a default implementation,
            delegating to the overloading mechanism and forcing the
            conversion back

            EXAMPLES::

                sage: A = HopfAlgebrasWithBasis(QQ).example(); A
                An example of Hopf algebra with basis: the group algebra of the Dihedral group of order 6 as a permutation group over Rational Field
                sage: [a,b] = A.algebra_generators()
                sage: a, A.counit(a)
                (B[(1,2,3)], 1)
                sage: b, A.counit(b)
                (B[(1,3)], 1)

            TODO: implement some tests of the axioms of coalgebras, bialgebras
            and Hopf algebras using the counit.
            """


        @abstract_method
        def coproduct(self, x):
            """
            Returns the coproduct of x.

            Eventually, there will be a default implementation,
            delegating to the overloading mechanism and forcing the
            conversion back

            EXAMPLES::

                sage: A = HopfAlgebrasWithBasis(QQ).example(); A
                An example of Hopf algebra with basis: the group algebra of the Dihedral group of order 6 as a permutation group over Rational Field
                sage: [a,b] = A.algebra_generators()
                sage: a, A.coproduct(a)
                (B[(1,2,3)], B[(1,2,3)] # B[(1,2,3)])
                sage: b, A.coproduct(b)
                (B[(1,3)], B[(1,3)] # B[(1,3)])
            """
            #return self.tensor_square()(overloaded_coproduct(x))

    class ElementMethods:
        def coproduct(self):
            """
            Returns the coproduct of ``self``

            EXAMPLES::

                sage: A = HopfAlgebrasWithBasis(QQ).example(); A
                An example of Hopf algebra with basis: the group algebra of the Dihedral group of order 6 as a permutation group over Rational Field
                sage: [a,b] = A.algebra_generators()
                sage: a, a.coproduct()
                (B[(1,2,3)], B[(1,2,3)] # B[(1,2,3)])
                sage: b, b.coproduct()
                (B[(1,3)], B[(1,3)] # B[(1,3)])
            """
            return self.parent().coproduct(self)

        def counit(self):
            """
            Returns the counit of ``self``

            EXAMPLES::

                sage: A = HopfAlgebrasWithBasis(QQ).example(); A
                An example of Hopf algebra with basis: the group algebra of the Dihedral group of order 6 as a permutation group over Rational Field
                sage: [a,b] = A.algebra_generators()
                sage: a, a.counit()
                (B[(1,2,3)], 1)
                sage: b, b.counit()
                (B[(1,3)], 1)
            """
            return self.parent().counit(self)

    class TensorProducts(TensorProductsCategory):

        @cached_method
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: Coalgebras(QQ).TensorProducts().extra_super_categories()
                [Category of coalgebras over Rational Field]
                sage: Coalgebras(QQ).TensorProducts().super_categories()
                [Category of tensor products of vector spaces over Rational Field,
                 Category of coalgebras over Rational Field]

            Meaning: a tensor product of coalgebras is a coalgebra
            """
            return [self.base_category()]

        class ParentMethods:
            # TODO: provide this default implementation of one if one_basis is not implemented
            #def one(self):
            #    return tensor(module.one() for module in self.modules)
            pass

        class ElementMethods:
            pass

    class DualObjects(DualObjectsCategory):

        def extra_super_categories(self):
            r"""
            Return the dual category.

            EXAMPLES:

            The category of coalgebras over the Rational Field is dual
            to the category of algebras over the same field::

                sage: C = Coalgebras(QQ)
                sage: C.dual()
                Category of duals of coalgebras over Rational Field
                sage: C.dual().super_categories() # indirect doctest
                [Category of algebras over Rational Field, Category of duals of vector spaces over Rational Field]

            .. WARNING::

                This is only correct in certain cases (finite dimension, ...).
                See :trac:`15647`.
            """
            from sage.categories.algebras import Algebras
            return [Algebras(self.base_category().base_ring())]

    class WithRealizations(WithRealizationsCategory):

        class ParentMethods:

            def coproduct(self, x):
                r"""
                Returns the coproduct of ``x``.

                EXAMPLES::

                    sage: N = NonCommutativeSymmetricFunctions(QQ)
                    sage: S = N.complete()
                    sage: N.coproduct.__module__
                    'sage.categories.coalgebras'
                    sage: N.coproduct(S[2])
                    S[] # S[2] + S[1] # S[1] + S[2] # S[]
                """
                return self.a_realization()(x).coproduct()

            def counit(self, x):
                r"""
                Returns the counit of ``x``.

                EXAMPLES::

                    sage: Sym = SymmetricFunctions(QQ)
                    sage: s = Sym.schur()
                    sage: f = s[2,1]
                    sage: f.counit.__module__
                    'sage.categories.coalgebras'
                    sage: f.counit()
                    0

                ::

                    sage: N = NonCommutativeSymmetricFunctions(QQ)
                    sage: N.counit.__module__
                    'sage.categories.coalgebras'
                    sage: N.counit(N.one())
                    1
                    sage: x = N.an_element(); x
                    2*S[] + 2*S[1] + 3*S[1, 1]
                    sage: N.counit(x)
                    2
                """
                return self.a_realization()(x).counit()

    class Realizations(RealizationsCategory):

        class ParentMethods:

            def coproduct_by_coercion(self, x):
                r"""
                Returns the coproduct by coercion if coproduct_by_basis is not implemented.

                EXAMPLES::

                    sage: Sym = SymmetricFunctions(QQ)
                    sage: m = Sym.monomial()
                    sage: f = m[2,1]
                    sage: f.coproduct.__module__
                    'sage.categories.coalgebras'
                    sage: m.coproduct_on_basis
                    NotImplemented
                    sage: m.coproduct == m.coproduct_by_coercion
                    True
                    sage: f.coproduct()
                    m[] # m[2, 1] + m[1] # m[2] + m[2] # m[1] + m[2, 1] # m[]

                ::

                    sage: N = NonCommutativeSymmetricFunctions(QQ)
                    sage: R = N.ribbon()
                    sage: R.coproduct_by_coercion.__module__
                    'sage.categories.coalgebras'
                    sage: R.coproduct_on_basis
                    NotImplemented
                    sage: R.coproduct == R.coproduct_by_coercion
                    True
                    sage: R[1].coproduct()
                    R[] # R[1] + R[1] # R[]
                """
                from sage.categories.tensor import tensor
                R = self.realization_of().a_realization()
                return self.tensor_square().sum(coeff * tensor([self(R[I]), self(R[J])])
                                                for ((I, J), coeff) in R(x).coproduct())
