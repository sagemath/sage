# sage_setup: distribution = sagemath-categories
# sage.doctest: needs sage.combinat
r"""
Coalgebras
"""
# ****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#  Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.categories.category_types import Category_over_base_ring
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.dual import DualObjectsCategory
from sage.categories.filtered_modules import FilteredModulesCategory
from sage.categories.modules import Modules
from sage.categories.realizations import RealizationsCategory
from sage.categories.super_modules import SuperModulesCategory
from sage.categories.tensor import TensorProductsCategory
from sage.categories.with_realizations import WithRealizationsCategory
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import LazyImport


class Coalgebras(Category_over_base_ring):
    """
    The category of coalgebras.

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

    WithBasis = LazyImport('sage.categories.coalgebras_with_basis', 'CoalgebrasWithBasis')
    Graded = LazyImport('sage.categories.graded_coalgebras', 'GradedCoalgebras')

    class ParentMethods:
        #def __init_add__(self): # The analogue of initDomainAdd
        #    # Will declare the coproduct of self to the coercion mechanism when it exists
        #    pass

        @abstract_method
        def counit(self, x):
            """
            Return the counit of ``x``.

            Eventually, there will be a default implementation,
            delegating to the overloading mechanism and forcing the
            conversion back

            EXAMPLES::

                sage: # needs sage.groups sage.modules
                sage: A = HopfAlgebrasWithBasis(QQ).example(); A
                An example of Hopf algebra with basis:
                 the group algebra of the Dihedral group of order 6 as a permutation group
                 over Rational Field
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
            Return the coproduct of ``x``.

            Eventually, there will be a default implementation,
            delegating to the overloading mechanism and forcing the
            conversion back

            EXAMPLES::

                sage: # needs sage.groups sage.modules
                sage: A = HopfAlgebrasWithBasis(QQ).example(); A
                An example of Hopf algebra with basis:
                 the group algebra of the Dihedral group of order 6 as a permutation group
                 over Rational Field
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
            Return the coproduct of ``self``.

            EXAMPLES::

                sage: # needs sage.groups sage.modules
                sage: A = HopfAlgebrasWithBasis(QQ).example(); A
                An example of Hopf algebra with basis:
                 the group algebra of the Dihedral group of order 6 as a permutation group
                 over Rational Field
                sage: [a,b] = A.algebra_generators()
                sage: a, a.coproduct()
                (B[(1,2,3)], B[(1,2,3)] # B[(1,2,3)])
                sage: b, b.coproduct()
                (B[(1,3)], B[(1,3)] # B[(1,3)])
            """
            return self.parent().coproduct(self)

        def counit(self):
            """
            Return the counit of ``self``.

            EXAMPLES::

                sage: # needs sage.groups sage.modules
                sage: A = HopfAlgebrasWithBasis(QQ).example(); A
                An example of Hopf algebra with basis:
                 the group algebra of the Dihedral group of order 6 as a permutation group
                 over Rational Field
                sage: [a,b] = A.algebra_generators()
                sage: a, a.counit()
                (B[(1,2,3)], 1)
                sage: b, b.counit()
                (B[(1,3)], 1)
            """
            return self.parent().counit(self)

    class SubcategoryMethods:
        @cached_method
        def Cocommutative(self):
            r"""
            Return the full subcategory of the cocommutative objects
            of ``self``.

            A coalgebra `C` is said to be *cocommutative* if

            .. MATH::

                \Delta(c) = \sum_{(c)} c_{(1)} \otimes c_{(2)}
                = \sum_{(c)} c_{(2)} \otimes c_{(1)}

            in Sweedler's notation for all `c \in C`.

            EXAMPLES::

                sage: C1 = Coalgebras(ZZ).Cocommutative().WithBasis(); C1
                Category of cocommutative coalgebras with basis over Integer Ring
                sage: C2 = Coalgebras(ZZ).WithBasis().Cocommutative()
                sage: C1 is C2
                True
                sage: BialgebrasWithBasis(QQ).Cocommutative()
                Category of cocommutative bialgebras with basis over Rational Field

            TESTS::

                sage: TestSuite(Coalgebras(ZZ).Cocommutative()).run()
            """
            return self._with_axiom("Cocommutative")

    class Cocommutative(CategoryWithAxiom_over_base_ring):
        """
        Category of cocommutative coalgebras.
        """

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
                [Category of algebras over Rational Field,
                 Category of duals of vector spaces over Rational Field]

            .. WARNING::

                This is only correct in certain cases (finite dimension, ...).
                See :issue:`15647`.
            """
            from sage.categories.algebras import Algebras
            return [Algebras(self.base_category().base_ring())]

    class Super(SuperModulesCategory):
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: Coalgebras(ZZ).Super().extra_super_categories()
                [Category of graded coalgebras over Integer Ring]
                sage: Coalgebras(ZZ).Super().super_categories()
                [Category of graded coalgebras over Integer Ring,
                 Category of super modules over Integer Ring]

            Compare this with the situation for bialgebras::

                sage: Bialgebras(ZZ).Super().extra_super_categories()
                []
                sage: Bialgebras(ZZ).Super().super_categories()
                [Category of super algebras over Integer Ring,
                 Category of super coalgebras over Integer Ring]

            The category of bialgebras does not occur in these results,
            since super bialgebras are not bialgebras.
            """
            return [self.base_category().Graded()]

        class SubcategoryMethods:
            @cached_method
            def Supercocommutative(self):
                r"""
                Return the full subcategory of the supercocommutative
                objects of ``self``.

                EXAMPLES::

                    sage: Coalgebras(ZZ).WithBasis().Super().Supercocommutative()
                    Category of supercocommutative super coalgebras with basis over Integer Ring
                    sage: BialgebrasWithBasis(QQ).Super().Supercocommutative()
                    Join of Category of super algebras with basis over Rational Field
                     and Category of super bialgebras over Rational Field
                     and Category of super coalgebras with basis over Rational Field
                     and Category of supercocommutative super coalgebras over Rational Field

                TESTS::

                    sage: TestSuite(HopfAlgebras(ZZ).Super().Supercocommutative()).run()
                """
                return self._with_axiom("Supercocommutative")

        class Supercocommutative(CategoryWithAxiom_over_base_ring):
            """
            Category of supercocommutative coalgebras.
            """

    class Filtered(FilteredModulesCategory):
        """
        Category of filtered coalgebras.
        """

    class WithRealizations(WithRealizationsCategory):

        class ParentMethods:

            def coproduct(self, x):
                r"""
                Return the coproduct of ``x``.

                EXAMPLES::

                    sage: # needs sage.modules
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
                Return the counit of ``x``.

                EXAMPLES::

                    sage: # needs sage.modules
                    sage: Sym = SymmetricFunctions(QQ)
                    sage: s = Sym.schur()
                    sage: f = s[2,1]
                    sage: f.counit.__module__
                    'sage.categories.coalgebras'
                    sage: f.counit()
                    0

                ::

                    sage: # needs sage.modules
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
                Return the coproduct by coercion if ``coproduct_by_basis``
                is not implemented.

                EXAMPLES::

                    sage: # needs sage.modules
                    sage: Sym = SymmetricFunctions(QQ)
                    sage: m = Sym.monomial()
                    sage: f = m[2,1]
                    sage: f.coproduct.__module__
                    'sage.categories.coalgebras'
                    sage: m.coproduct_on_basis
                    NotImplemented
                    sage: m.coproduct == m.coproduct_by_coercion
                    True
                    sage: f.coproduct()                                                 # needs lrcalc_python
                    m[] # m[2, 1] + m[1] # m[2] + m[2] # m[1] + m[2, 1] # m[]

                ::

                    sage: # needs sage.modules
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
                R = self.realization_of().a_realization()
                return self.tensor_square()(R(x).coproduct())

            def counit_by_coercion(self, x):
                r"""
                Return the counit of ``x`` if ``counit_by_basis`` is
                not implemented.

                EXAMPLES::

                    sage: # needs lrcalc_python sage.modules
                    sage: sp = SymmetricFunctions(QQ).sp()
                    sage: sp.an_element()
                    2*sp[] + 2*sp[1] + 3*sp[2]
                    sage: sp.counit(sp.an_element())
                    2
                    sage: o = SymmetricFunctions(QQ).o()
                    sage: o.an_element()
                    2*o[] + 2*o[1] + 3*o[2]
                    sage: o.counit(o.an_element())
                    -1
                """
                R = self.realization_of().a_realization()
                return R(x).counit()
