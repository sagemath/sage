# sage_setup: distribution = sagemath-categories
r"""
Modules
"""
# ****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.categories.bimodules import Bimodules
from sage.categories.cartesian_product import CartesianProductsCategory
from sage.categories.category import Category
from sage.categories.category_types import Category_module
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.dual import DualObjectsCategory
from sage.categories.fields import Fields
from sage.categories.homset import Hom
from sage.categories.homsets import HomsetsCategory
from sage.categories.morphism import SetMorphism
from sage.categories.sets_cat import Sets
from sage.categories.tensor import TensorProductFunctor, TensorProductsCategory, tensor
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import LazyImport

_Fields = Fields()


class Modules(Category_module):
    r"""
    The category of all modules over a base ring `R`.

    An `R`-module `M` is a left and right `R`-module over a
    commutative ring `R` such that:

    .. MATH::

        r*(x*s) = (r*x)*s \qquad  \forall r,s \in R \text{ and } x \in M

    INPUT:

    - ``base_ring`` -- a ring `R` or subcategory of ``Rings()``
    - ``dispatch`` -- boolean (for internal use; default: ``True``)

    When the base ring is a field, the category of vector spaces is
    returned instead (unless ``dispatch == False``).

    .. WARNING::

        Outside of the context of symmetric modules over a commutative
        ring, the specifications of this category are fuzzy and not
        yet set in stone (see below). The code in this category and
        its subcategories is therefore prone to bugs or arbitrary
        limitations in this case.

    EXAMPLES::

        sage: Modules(ZZ)
        Category of modules over Integer Ring
        sage: Modules(QQ)
        Category of vector spaces over Rational Field

        sage: Modules(Rings())
        Category of modules over rings
        sage: Modules(FiniteFields())
        Category of vector spaces over finite enumerated fields

        sage: Modules(Integers(9))
        Category of modules over Ring of integers modulo 9

        sage: Modules(Integers(9)).super_categories()
        [Category of bimodules over Ring of integers modulo 9 on the left
                                and Ring of integers modulo 9 on the right]

        sage: Modules(ZZ).super_categories()
        [Category of bimodules over Integer Ring on the left
                                and Integer Ring on the right]

        sage: Modules == RingModules
        True

        sage: Modules(ZZ['x']).is_abelian()   # see #6081
        True

    TESTS::

        sage: TestSuite(Modules(ZZ)).run()

    .. TODO::

        - Clarify the distinction, if any, with ``BiModules(R, R)``.
          In particular, if `R` is a commutative ring (e.g. a field),
          some pieces of the code possibly assume that `M` is a
          *symmetric `R`-`R`-bimodule*:

          .. MATH::

              r*x = x*r \qquad  \forall r \in R \text{ and } x \in M

        - Make sure that non symmetric modules are properly supported
          by all the code, and advertise it.

        - Make sure that non commutative rings are properly supported
          by all the code, and advertise it.

        - Add support for base semirings.

        - Implement a ``FreeModules(R)`` category, when so prompted by a
          concrete use case: e.g.  modeling a free module with several
          bases (using :meth:`Sets.SubcategoryMethods.Realizations`)
          or with an atlas of local maps (see e.g. :issue:`15916`).
    """

    @staticmethod
    def __classcall_private__(cls, base_ring, dispatch=True):
        r"""
        Implement the dispatching of ``Modules(field)`` to
        ``VectorSpaces(field)``.

        This feature will later be extended, probably as a covariant
        functorial construction, to support modules over various kinds
        of rings (principal ideal domains, ...), or even over semirings.

        TESTS::

            sage: C = Modules(ZZ); C
            Category of modules over Integer Ring
            sage: C is Modules(ZZ, dispatch = False)
            True
            sage: C is Modules(ZZ, dispatch = True)
            True
            sage: C._reduction
            (<class 'sage.categories.modules.Modules'>, (Integer Ring,), {'dispatch': False})
            sage: TestSuite(C).run()

            sage: Modules(QQ) is VectorSpaces(QQ)
            True
            sage: Modules(QQ, dispatch = True) is VectorSpaces(QQ)
            True

            sage: C = Modules(NonNegativeIntegers()); C   # todo: not implemented
            Category of semiring modules over Non negative integers

            sage: C = Modules(QQ, dispatch = False); C
            Category of modules over Rational Field
            sage: C._reduction
            (<class 'sage.categories.modules.Modules'>, (Rational Field,), {'dispatch': False})
            sage: TestSuite(C).run()
        """
        if dispatch:
            if base_ring in _Fields or (isinstance(base_ring, Category)
                                        and base_ring.is_subcategory(_Fields)):
                from .vector_spaces import VectorSpaces
                return VectorSpaces(base_ring, check=False)
        result = super().__classcall__(cls, base_ring)
        result._reduction[2]['dispatch'] = False
        return result

    def super_categories(self):
        """
        EXAMPLES::

            sage: Modules(ZZ).super_categories()
            [Category of bimodules over Integer Ring on the left
                                    and Integer Ring on the right]

        Nota bene::

            sage: Modules(QQ)
            Category of vector spaces over Rational Field
            sage: Modules(QQ).super_categories()
            [Category of modules over Rational Field]
        """
        R = self.base_ring()
        return [Bimodules(R, R)]

    def additional_structure(self):
        r"""
        Return ``None``.

        Indeed, the category of modules defines no additional structure:
        a bimodule morphism between two modules is a module morphism.

        .. SEEALSO:: :meth:`Category.additional_structure`

        .. TODO:: Should this category be a :class:`~sage.categories.category_with_axiom.CategoryWithAxiom`?

        EXAMPLES::

            sage: Modules(ZZ).additional_structure()
        """
        return None

    class SubcategoryMethods:

        @cached_method
        def base_ring(self):
            r"""
            Return the base ring (category) for ``self``.

            This implements a ``base_ring`` method for all
            subcategories of ``Modules(K)``.

            EXAMPLES::

                sage: C = Modules(QQ) & Semigroups(); C
                Join of Category of semigroups
                    and Category of vector spaces over Rational Field
                sage: C.base_ring()
                Rational Field
                sage: C.base_ring.__module__
                'sage.categories.modules'

                sage: C2 = Modules(Rings()) & Semigroups(); C2
                Join of Category of semigroups and Category of modules over rings
                sage: C2.base_ring()
                Category of rings
                sage: C2.base_ring.__module__
                'sage.categories.modules'

                sage: # needs sage.combinat sage.groups sage.modules
                sage: C3 = DescentAlgebra(QQ,3).B().category()
                sage: C3.base_ring.__module__
                'sage.categories.modules'
                sage: C3.base_ring()
                Rational Field

                sage: # needs sage.combinat sage.modules
                sage: C4 = QuasiSymmetricFunctions(QQ).F().category()
                sage: C4.base_ring.__module__
                'sage.categories.modules'
                sage: C4.base_ring()
                Rational Field
            """
            for C in self.super_categories():
                # Is there a better way to ask if C is a subcategory of Modules?
                if hasattr(C, "base_ring"):
                    return C.base_ring()
            assert False, "some super category of {} should be a category over base ring".format(self)

        def TensorProducts(self):
            r"""
            Return the full subcategory of objects of ``self`` constructed
            as tensor products.

            .. SEEALSO::

                - :class:`.tensor.TensorProductsCategory`
                - :class:`~.covariant_functorial_construction.RegressiveCovariantFunctorialConstruction`.

            EXAMPLES::

                sage: ModulesWithBasis(QQ).TensorProducts()
                Category of tensor products of vector spaces with basis over Rational Field
            """
            return TensorProductsCategory.category_of(self)

        @cached_method
        def DualObjects(self):
            r"""
            Return the category of spaces constructed as duals of
            spaces of ``self``.

            The *dual* of a vector space `V` is the space consisting of
            all linear functionals on `V` (see :wikipedia:`Dual_space`).
            Additional structure on `V` can endow its dual with
            additional structure; for example, if `V` is a finite
            dimensional algebra, then its dual is a coalgebra.

            This returns the category of spaces constructed as dual of
            spaces in ``self``, endowed with the appropriate
            additional structure.

            .. WARNING::

                - This semantic of ``dual`` and ``DualObject`` is
                  imposed on all subcategories, in particular to make
                  ``dual`` a covariant functorial construction.

                  A subcategory that defines a different notion of
                  dual needs to use a different name.

                - Typically, the category of graded modules should
                  define a separate ``graded_dual`` construction (see
                  :issue:`15647`). For now the two constructions are
                  not distinguished which is an oversimplified model.

            .. SEEALSO::

                - :class:`.dual.DualObjectsCategory`
                - :class:`~.covariant_functorial_construction.CovariantFunctorialConstruction`.

            EXAMPLES::

                sage: VectorSpaces(QQ).DualObjects()
                Category of duals of vector spaces over Rational Field

            The dual of a vector space is a vector space::

                sage: VectorSpaces(QQ).DualObjects().super_categories()
                [Category of vector spaces over Rational Field]

            The dual of an algebra is a coalgebra::

                sage: sorted(Algebras(QQ).DualObjects().super_categories(), key=str)
                [Category of coalgebras over Rational Field,
                 Category of duals of vector spaces over Rational Field]

            The dual of a coalgebra is an algebra::

                sage: sorted(Coalgebras(QQ).DualObjects().super_categories(), key=str)
                [Category of algebras over Rational Field,
                 Category of duals of vector spaces over Rational Field]

            As a shorthand, this category can be accessed with the
            :meth:`~Modules.SubcategoryMethods.dual` method::

                sage: VectorSpaces(QQ).dual()
                Category of duals of vector spaces over Rational Field

            TESTS::

                sage: C = VectorSpaces(QQ).DualObjects()
                sage: C.base_category()
                Category of vector spaces over Rational Field
                sage: C.super_categories()
                [Category of vector spaces over Rational Field]
                sage: latex(C)
                \mathbf{DualObjects}(\mathbf{VectorSpaces}_{\Bold{Q}})
                sage: TestSuite(C).run()
            """
            return DualObjectsCategory.category_of(self)

        dual = DualObjects

        @cached_method
        def FiniteDimensional(self):
            r"""
            Return the full subcategory of the finite dimensional objects of ``self``.

            EXAMPLES::

                sage: Modules(ZZ).FiniteDimensional()
                Category of finite dimensional modules over Integer Ring
                sage: Coalgebras(QQ).FiniteDimensional()
                Category of finite dimensional coalgebras over Rational Field
                sage: AlgebrasWithBasis(QQ).FiniteDimensional()
                Category of finite dimensional algebras with basis over Rational Field

            TESTS::

                sage: TestSuite(Modules(ZZ).FiniteDimensional()).run()
                sage: Coalgebras(QQ).FiniteDimensional.__module__
                'sage.categories.modules'
            """
            return self._with_axiom("FiniteDimensional")

        @cached_method
        def FinitelyPresented(self):
            r"""
            Return the full subcategory of the finitely presented objects of ``self``.

            EXAMPLES::

                sage: Modules(ZZ).FinitelyPresented()
                Category of finitely presented modules over Integer Ring
                sage: A = SteenrodAlgebra(2)                                            # needs sage.combinat sage.modules
                sage: from sage.modules.fp_graded.module import FPModule                # needs sage.combinat sage.modules
                sage: FPModule(A, [0, 1], [[Sq(2), Sq(1)]]).category()                  # needs sage.combinat sage.modules
                Category of finitely presented graded modules
                 over mod 2 Steenrod algebra, milnor basis

            TESTS::

                sage: TestSuite(Modules(ZZ).FinitelyPresented()).run()
            """
            return self._with_axiom("FinitelyPresented")

        @cached_method
        def Filtered(self, base_ring=None):
            r"""
            Return the subcategory of the filtered objects of ``self``.

            INPUT:

            - ``base_ring`` -- this is ignored

            EXAMPLES::

                sage: Modules(ZZ).Filtered()
                Category of filtered modules over Integer Ring

                sage: Coalgebras(QQ).Filtered()
                Category of filtered coalgebras over Rational Field

                sage: AlgebrasWithBasis(QQ).Filtered()
                Category of filtered algebras with basis over Rational Field

            .. TODO::

                - Explain why this does not commute with :meth:`WithBasis`
                - Improve the support for covariant functorial
                  constructions categories over a base ring so as to
                  get rid of the ``base_ring`` argument.

            TESTS::

                sage: Coalgebras(QQ).Graded.__module__
                'sage.categories.modules'
            """
            assert base_ring is None or base_ring is self.base_ring()
            from sage.categories.filtered_modules import FilteredModulesCategory
            return FilteredModulesCategory.category_of(self)

        @cached_method
        def Graded(self, base_ring=None):
            r"""
            Return the subcategory of the graded objects of ``self``.

            INPUT:

            - ``base_ring`` -- this is ignored

            EXAMPLES::

                sage: Modules(ZZ).Graded()
                Category of graded modules over Integer Ring

                sage: Coalgebras(QQ).Graded()
                Category of graded coalgebras over Rational Field

                sage: AlgebrasWithBasis(QQ).Graded()
                Category of graded algebras with basis over Rational Field

            .. TODO::

                - Explain why this does not commute with :meth:`WithBasis`
                - Improve the support for covariant functorial
                  constructions categories over a base ring so as to
                  get rid of the ``base_ring`` argument.

            TESTS::

                sage: Coalgebras(QQ).Graded.__module__
                'sage.categories.modules'
            """
            assert base_ring is None or base_ring is self.base_ring()
            from sage.categories.graded_modules import GradedModulesCategory
            return GradedModulesCategory.category_of(self)

        @cached_method
        def Super(self, base_ring=None):
            r"""
            Return the super-analogue category of ``self``.

            INPUT:

            - ``base_ring`` -- this is ignored

            EXAMPLES::

                sage: Modules(ZZ).Super()
                Category of super modules over Integer Ring

                sage: Coalgebras(QQ).Super()
                Category of super coalgebras over Rational Field

                sage: AlgebrasWithBasis(QQ).Super()
                Category of super algebras with basis over Rational Field

            .. TODO::

                - Explain why this does not commute with :meth:`WithBasis`
                - Improve the support for covariant functorial
                  constructions categories over a base ring so as to
                  get rid of the ``base_ring`` argument.

            TESTS::

                sage: Coalgebras(QQ).Super.__module__
                'sage.categories.modules'
            """
            assert base_ring is None or base_ring is self.base_ring()
            from sage.categories.super_modules import SuperModulesCategory
            return SuperModulesCategory.category_of(self)

        @cached_method
        def WithBasis(self):
            r"""
            Return the full subcategory of the objects of ``self`` with
            a distinguished basis.

            EXAMPLES::

                sage: Modules(ZZ).WithBasis()
                Category of modules with basis over Integer Ring
                sage: Coalgebras(QQ).WithBasis()
                Category of coalgebras with basis over Rational Field
                sage: AlgebrasWithBasis(QQ).WithBasis()
                Category of algebras with basis over Rational Field

            TESTS::

                sage: TestSuite(Modules(ZZ).WithBasis()).run()
                sage: Coalgebras(QQ).WithBasis.__module__
                'sage.categories.modules'
            """
            return self._with_axiom("WithBasis")

    class FiniteDimensional(CategoryWithAxiom_over_base_ring):

        def extra_super_categories(self):
            """
            Implement the fact that a finite dimensional module over a finite
            ring is finite.

            EXAMPLES::

                sage: Modules(IntegerModRing(4)).FiniteDimensional().extra_super_categories()
                [Category of finite sets]
                sage: Modules(ZZ).FiniteDimensional().extra_super_categories()
                []
                sage: Modules(GF(5)).FiniteDimensional().is_subcategory(Sets().Finite())
                True
                sage: Modules(ZZ).FiniteDimensional().is_subcategory(Sets().Finite())
                False

                sage: Modules(Rings().Finite()).FiniteDimensional().is_subcategory(Sets().Finite())
                True
                sage: Modules(Rings()).FiniteDimensional().is_subcategory(Sets().Finite())
                False
            """
            base_ring = self.base_ring()
            FiniteSets = Sets().Finite()
            if (isinstance(base_ring, Category) and
                    base_ring.is_subcategory(FiniteSets)) or \
                base_ring in FiniteSets:
                return [FiniteSets]
            else:
                return []

        class TensorProducts(TensorProductsCategory):

            def extra_super_categories(self):
                """
                Implement the fact that a (finite) tensor product of
                finite dimensional modules is a finite dimensional module.

                EXAMPLES::

                    sage: Modules(ZZ).FiniteDimensional().TensorProducts().extra_super_categories()
                    [Category of finite dimensional modules over Integer Ring]
                    sage: Modules(QQ).FiniteDimensional().TensorProducts().FiniteDimensional()
                    Category of tensor products of finite dimensional vector spaces
                     over Rational Field
                """
                return [self.base_category()]

    class FinitelyPresented(CategoryWithAxiom_over_base_ring):

        def extra_super_categories(self):
            """
            Implement the fact that a finitely presented module over a finite
            ring is finite.

            EXAMPLES::

                sage: Modules(IntegerModRing(4)).FiniteDimensional().extra_super_categories()
                [Category of finite sets]
                sage: Modules(ZZ).FiniteDimensional().extra_super_categories()
                []
                sage: Modules(GF(5)).FiniteDimensional().is_subcategory(Sets().Finite())
                True
                sage: Modules(ZZ).FiniteDimensional().is_subcategory(Sets().Finite())
                False

                sage: Modules(Rings().Finite()).FiniteDimensional().is_subcategory(Sets().Finite())
                True
                sage: Modules(Rings()).FiniteDimensional().is_subcategory(Sets().Finite())
                False
            """
            base_ring = self.base_ring()
            FiniteSets = Sets().Finite()
            if (isinstance(base_ring, Category) and
                    base_ring.is_subcategory(FiniteSets)) or \
                base_ring in FiniteSets:
                return [FiniteSets]
            else:
                return []

    Filtered = LazyImport('sage.categories.filtered_modules', 'FilteredModules')
    Graded = LazyImport('sage.categories.graded_modules', 'GradedModules')
    Super = LazyImport('sage.categories.super_modules', 'SuperModules')
    # at_startup currently needed for MatrixSpace, see #22955 (e.g., comment:20)
    WithBasis = LazyImport('sage.categories.modules_with_basis', 'ModulesWithBasis',
                           at_startup=True)

    class ParentMethods:

        def linear_combination(self, iter_of_elements_coeff, factor_on_left=True):
            r"""
            Return the linear combination `\lambda_1 v_1 + \cdots +
            \lambda_k v_k` (resp.  the linear combination `v_1 \lambda_1 +
            \cdots + v_k \lambda_k`) where ``iter_of_elements_coeff`` iterates
            through the sequence `((\lambda_1, v_1), ..., (\lambda_k, v_k))`.

            INPUT:

            - ``iter_of_elements_coeff`` -- iterator of pairs
              ``(element, coeff)`` with ``element`` in ``self`` and
              ``coeff`` in ``self.base_ring()``

            - ``factor_on_left`` -- (optional) if ``True``, the coefficients
              are multiplied from the left; if ``False``, the coefficients
              are multiplied from the right

            EXAMPLES::

                sage: m = matrix([[0,1], [1,1]])                                        # needs sage.modules
                sage: J.<a,b,c> = JordanAlgebra(m)                                      # needs sage.combinat sage.modules
                sage: J.linear_combination(((a+b, 1), (-2*b + c, -1)))                  # needs sage.combinat sage.modules
                1 + (3, -1)
            """
            if factor_on_left:
                return self.sum(coeff * element
                                for element, coeff in iter_of_elements_coeff)
            else:
                return self.sum(element * coeff
                                for element, coeff in iter_of_elements_coeff)

        @cached_method
        def tensor_square(self):
            """
            Return the tensor square of ``self``.

            EXAMPLES::

                sage: A = HopfAlgebrasWithBasis(QQ).example()                           # needs sage.groups sage.modules
                sage: A.tensor_square()                                                 # needs sage.groups sage.modules
                An example of Hopf algebra with basis:
                 the group algebra of the Dihedral group of order 6
                 as a permutation group over Rational Field # An example
                 of Hopf algebra with basis: the group algebra of the Dihedral
                 group of order 6 as a permutation group over Rational Field
            """
            return tensor([self, self])

        def module_morphism(self, *, function, category=None, codomain, **keywords):
            r"""
            Construct a module morphism from ``self`` to ``codomain``.

            Let ``self`` be a module `X` over a ring `R`.
            This constructs a morphism `f: X \to Y`.

            INPUT:

            - ``self`` -- a parent `X` in ``Modules(R)``

            - ``function`` -- a function `f` from `X` to `Y`

            - ``codomain`` -- the codomain `Y` of the morphism (default:
              ``f.codomain()`` if it's defined; otherwise it must be specified)

            - ``category`` -- a category or ``None`` (default: ``None``)

            EXAMPLES::

                sage: # needs sage.modules
                sage: V = FiniteRankFreeModule(QQ, 2)
                sage: e = V.basis('e'); e
                Basis (e_0,e_1) on the
                 2-dimensional vector space over the Rational Field
                sage: neg = V.module_morphism(function=operator.neg, codomain=V); neg
                Generic endomorphism of
                 2-dimensional vector space over the Rational Field
                sage: neg(e[0])
                Element -e_0 of the 2-dimensional vector space over the Rational Field
            """
            # Make sure that we only create a module morphism, even if
            # domain and codomain have more structure
            if category is None:
                category = Modules(self.base_ring())
            return SetMorphism(Hom(self, codomain, category), function)

        def quotient(self, submodule, check=True, **kwds):
            r"""
            Construct the quotient module ``self`` / ``submodule``.

            INPUT:

            - ``submodule`` -- a submodule with basis of ``self``, or
              something that can be turned into one via
              ``self.submodule(submodule)``

            - ``check``, other keyword arguments -- passed on to
              :meth:`quotient_module`.

            This method just delegates to :meth:`quotient_module`.
            Classes implementing modules should override that method.

            Parents in categories with additional structure may override
            :meth:`quotient`. For example, in algebras, :meth:`quotient` will
            be the same as :meth:`quotient_ring`.

            EXAMPLES::

                sage: C = CombinatorialFreeModule(QQ, ['a','b','c'])                    # needs sage.modules
                sage: TA = TensorAlgebra(C)                                             # needs sage.combinat sage.modules
                sage: TA.quotient                                                       # needs sage.combinat sage.modules
                <bound method Rings.ParentMethods.quotient of
                 Tensor Algebra of Free module generated by {'a', 'b', 'c'}
                 over Rational Field>
            """
            return self.quotient_module(submodule, check=check, **kwds)

    class ElementMethods:
        pass

    class Homsets(HomsetsCategory):
        r"""
        The category of homomorphism sets `\hom(X,Y)` for `X`, `Y` modules.
        """

        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: Modules(ZZ).Homsets().extra_super_categories()
                [Category of modules over Integer Ring]
            """
            return [Modules(self.base_category().base_ring())]

        def base_ring(self):
            """
            EXAMPLES::

                sage: Modules(ZZ).Homsets().base_ring()
                Integer Ring

            .. TODO::

                Generalize this so that any homset category of a full
                subcategory of modules over a base ring is a category over
                this base ring.
            """
            return self.base_category().base_ring()

        class ParentMethods:

            @cached_method
            def base_ring(self):
                """
                Return the base ring of ``self``.

                EXAMPLES::

                    sage: # needs sage.modules
                    sage: E = CombinatorialFreeModule(ZZ, [1,2,3])
                    sage: F = CombinatorialFreeModule(ZZ, [2,3,4])
                    sage: H = Hom(E, F)
                    sage: H.base_ring()
                    Integer Ring

                This ``base_ring`` method is actually overridden by
                :meth:`sage.structure.category_object.CategoryObject.base_ring`::

                    sage: H.base_ring.__module__                                        # needs sage.modules

                Here we call it directly::

                    sage: method = H.category().parent_class.base_ring                  # needs sage.modules
                    sage: method.__get__(H)()                                           # needs sage.modules
                    Integer Ring
                """
                return self.domain().base_ring()

            @cached_method
            def zero(self):
                """
                EXAMPLES::

                    sage: # needs sage.modules
                    sage: E = CombinatorialFreeModule(ZZ, [1,2,3])
                    sage: F = CombinatorialFreeModule(ZZ, [2,3,4])
                    sage: H = Hom(E, F)
                    sage: f = H.zero()
                    sage: f
                    Generic morphism:
                      From: Free module generated by {1, 2, 3} over Integer Ring
                      To:   Free module generated by {2, 3, 4} over Integer Ring
                    sage: f(E.monomial(2))
                    0
                    sage: f(E.monomial(3)) == F.zero()
                    True

                TESTS:

                We check that ``H.zero()`` is picklable::

                    sage: loads(dumps(f.parent().zero()))                               # needs sage.modules
                    Generic morphism:
                      From: Free module generated by {1, 2, 3} over Integer Ring
                      To:   Free module generated by {2, 3, 4} over Integer Ring
                """
                from sage.misc.constant_function import ConstantFunction
                return self(ConstantFunction(self.codomain().zero()))

        class Endset(CategoryWithAxiom_over_base_ring):
            """
            The category of endomorphism sets `End(X)` for `X`
            a module (this is not used yet)
            """
            def extra_super_categories(self):
                """
                Implement the fact that the endomorphism set of a module is an algebra.

                .. SEEALSO:: :meth:`CategoryWithAxiom.extra_super_categories`

                EXAMPLES::

                    sage: Modules(ZZ).Endsets().extra_super_categories()
                    [Category of magmatic algebras over Integer Ring]

                    sage: End(ZZ^3) in Algebras(ZZ)                                     # needs sage.modules
                    True
                """
                from .magmatic_algebras import MagmaticAlgebras
                return [MagmaticAlgebras(self.base_category().base_ring())]

    class CartesianProducts(CartesianProductsCategory):
        """
        The category of modules constructed as Cartesian products of modules.

        This construction gives the direct product of modules. The
        implementation is based on the following resources:

        - http://groups.google.fr/group/sage-devel/browse_thread/thread/35a72b1d0a2fc77a/348f42ae77a66d16#348f42ae77a66d16
        - :wikipedia:`Direct_product`
        """
        def extra_super_categories(self):
            """
            A Cartesian product of modules is endowed with a natural
            module structure.

            EXAMPLES::

                sage: Modules(ZZ).CartesianProducts().extra_super_categories()
                [Category of modules over Integer Ring]
                sage: Modules(ZZ).CartesianProducts().super_categories()
                [Category of Cartesian products of commutative additive groups,
                 Category of modules over Integer Ring]
            """
            return [self.base_category()]

        class ParentMethods:

            def __init_extra__(self):
                """
                Initialise the base ring of this Cartesian product.

                EXAMPLES::

                    sage: # needs sage.modules
                    sage: E = CombinatorialFreeModule(ZZ, [1,2,3])
                    sage: F = CombinatorialFreeModule(ZZ, [2,3,4])
                    sage: C = cartesian_product([E, F]); C
                    Free module generated by {1, 2, 3} over Integer Ring (+)
                    Free module generated by {2, 3, 4} over Integer Ring
                    sage: C.base_ring()
                    Integer Ring

                Check that :issue:`29225` is fixed::

                    sage: M = cartesian_product((ZZ^2, ZZ^3)); M                        # needs sage.modules
                    The Cartesian product of
                     (Ambient free module of rank 2 over the principal ideal domain Integer Ring,
                      Ambient free module of rank 3 over the principal ideal domain Integer Ring)
                    sage: M.category()                                                  # needs sage.modules
                    Category of Cartesian products of modules with basis
                     over (Dedekind domains and euclidean domains
                      and noetherian rings
                      and infinite enumerated sets and metric spaces)
                    sage: M.base_ring()                                                 # needs sage.modules
                    Integer Ring

                    sage: A = cartesian_product((QQ^2, QQ['x'])); A                     # needs sage.modules
                    The Cartesian product of
                     (Vector space of dimension 2 over Rational Field,
                      Univariate Polynomial Ring in x over Rational Field)
                    sage: A.category()                                                  # needs sage.modules
                    Category of Cartesian products of vector spaces with basis
                     over (number fields and quotient fields and metric spaces)
                    sage: A.base_ring()                                                 # needs sage.modules
                    Rational Field

                This currently only works if all factors have the same
                base ring::

                    sage: B = cartesian_product((ZZ['x'], QQ^3)); B                     # needs sage.modules
                    The Cartesian product of
                     (Univariate Polynomial Ring in x over Integer Ring,
                      Vector space of dimension 3 over Rational Field)
                    sage: B.category()                                                  # needs sage.modules
                    Category of Cartesian products of commutative additive groups
                    sage: B.base_ring()                                                 # needs sage.modules
                """
                factors = self._sets
                if factors:
                    R = factors[0].base_ring()
                    if all(A.base_ring() is R for A in factors):
                        self._base = R

        class ElementMethods:

            def _lmul_(self, x):
                """
                Return the product of `x` with ``self``.

                EXAMPLES::

                    sage: A = FreeModule(ZZ, 2)                                         # needs sage.modules
                    sage: B = cartesian_product([A, A]); B                              # needs sage.modules
                    The Cartesian product of
                     (Ambient free module of rank 2 over the principal ideal domain Integer Ring,
                      Ambient free module of rank 2 over the principal ideal domain Integer Ring)
                    sage: 5*B(([1, 2], [3, 4]))                                         # needs sage.modules
                    ((5, 10), (15, 20))
                """
                return self.parent()._cartesian_product_of_elements(
                    x * y for y in self.cartesian_factors())

    class TensorProducts(TensorProductsCategory):
        """
        The category of modules constructed by tensor product of modules.
        """
        @cached_method
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: Modules(ZZ).TensorProducts().extra_super_categories()
                [Category of modules over Integer Ring]
                sage: Modules(ZZ).TensorProducts().super_categories()
                [Category of modules over Integer Ring]
            """
            return [self.base_category()]

        class ParentMethods:
            """
            Implement operations on tensor products of modules.
            """
            def construction(self):
                """
                Return the construction of ``self``.

                EXAMPLES::

                    sage: A = algebras.Free(QQ, 2)                                      # needs sage.combinat sage.modules
                    sage: T = A.tensor(A)                                               # needs sage.combinat sage.modules
                    sage: T.construction()                                              # needs sage.combinat sage.modules
                    (The tensor functorial construction,
                     (Free Algebra on 2 generators (None0, None1) over Rational Field,
                      Free Algebra on 2 generators (None0, None1) over Rational Field))
                """
                try:
                    factors = self.tensor_factors()
                except (TypeError, NotImplementedError):
                    from sage.misc.superseded import deprecation
                    deprecation(34393, "implementations of Modules().TensorProducts() now must define the method tensor_factors")
                    return None
                return (TensorProductFunctor(),
                        factors)

            @abstract_method(optional=True)
            def tensor_factors(self):
                """
                Return the tensor factors of this tensor product.

                EXAMPLES::

                    sage: # needs sage.modules
                    sage: F = CombinatorialFreeModule(ZZ, [1,2])
                    sage: F.rename('F')
                    sage: G = CombinatorialFreeModule(ZZ, [3,4])
                    sage: G.rename('G')
                    sage: T = tensor([F, G]); T
                    F # G
                    sage: T.tensor_factors()
                    (F, G)

                TESTS::

                    sage: Cat = ModulesWithBasis(ZZ).FiniteDimensional().TensorProducts()
                    sage: M = CombinatorialFreeModule(ZZ,                               # needs sage.modules
                    ....:                             ((1, 1), (1, 2), (2, 1), (2, 2)),
                    ....:                             category=Cat)
                    sage: M.construction()                                              # needs sage.modules
                    doctest:warning...
                    DeprecationWarning: implementations of Modules().TensorProducts() now must define the method tensor_factors
                    See https://github.com/sagemath/sage/issues/34393 for details.
                    (VectorFunctor, Integer Ring)
                """
