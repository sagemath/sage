r"""
Modules
"""
from __future__ import absolute_import
#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import LazyImport
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.homsets import HomsetsCategory
from .category import Category, JoinCategory
from .category_types import Category_module, Category_over_base_ring
import sage.categories.coercion_methods
from sage.categories.tensor import TensorProductsCategory, tensor
from .dual import DualObjectsCategory
from sage.categories.cartesian_product import CartesianProductsCategory
from sage.categories.sets_cat import Sets
from sage.categories.bimodules import Bimodules
from sage.categories.fields import Fields
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
    - ``dispatch`` -- a boolean (for internal use; default: ``True``)

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
        Category of vector spaces over finite fields

        sage: Modules(Integers(9))
        Category of modules over Ring of integers modulo 9

        sage: Modules(Integers(9)).super_categories()
        [Category of bimodules over Ring of integers modulo 9 on the left and Ring of integers modulo 9 on the right]

        sage: Modules(ZZ).super_categories()
        [Category of bimodules over Integer Ring on the left and Integer Ring on the right]

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
          or with an atlas of local maps (see e.g. :trac:`15916`).
    """

    @staticmethod
    def __classcall_private__(cls, base_ring, dispatch = True):
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
        result = super(Modules, cls).__classcall__(cls, base_ring)
        result._reduction[2]['dispatch'] = False
        return result

    def super_categories(self):
        """
        EXAMPLES::

            sage: Modules(ZZ).super_categories()
            [Category of bimodules over Integer Ring on the left and Integer Ring on the right]

        Nota bene::

            sage: Modules(QQ)
            Category of vector spaces over Rational Field
            sage: Modules(QQ).super_categories()
            [Category of modules over Rational Field]
        """
        R = self.base_ring()
        return [Bimodules(R,R)]

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
                Join of Category of semigroups and Category of vector spaces over Rational Field
                sage: C.base_ring()
                Rational Field
                sage: C.base_ring.__module__
                'sage.categories.modules'

                sage: C = Modules(Rings()) & Semigroups(); C
                Join of Category of semigroups and Category of modules over rings
                sage: C.base_ring()
                Category of rings
                sage: C.base_ring.__module__
                'sage.categories.modules'

                sage: C = DescentAlgebra(QQ,3).B().category()
                sage: C.base_ring.__module__
                'sage.categories.modules'
                sage: C.base_ring()
                Rational Field

                sage: C = QuasiSymmetricFunctions(QQ).F().category()
                sage: C.base_ring.__module__
                'sage.categories.modules'
                sage: C.base_ring()
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
                  :trac:`15647`). For now the two constructions are
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
        def Filtered(self, base_ring=None):
            r"""
            Return the subcategory of the filtered objects of ``self``.

            INPUT:

            - ``base_ring`` -- this is ignored

            EXAMPLES::

                sage: Modules(ZZ).Filtered()
                Category of filtered modules over Integer Ring

                sage: Coalgebras(QQ).Filtered()
                Join of Category of filtered modules over Rational Field
                 and Category of coalgebras over Rational Field

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
                Join of Category of graded modules over Rational Field and Category of coalgebras over Rational Field

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
            """
            if self.base_ring() in Sets().Finite():
                return [Sets().Finite()]
            else:
                return []

    Filtered = LazyImport('sage.categories.filtered_modules', 'FilteredModules')
    Graded = LazyImport('sage.categories.graded_modules', 'GradedModules')
    Super = LazyImport('sage.categories.super_modules', 'SuperModules')
    WithBasis = LazyImport('sage.categories.modules_with_basis', 'ModulesWithBasis')

    class ParentMethods:
        @cached_method
        def tensor_square(self):
            """
            Returns the tensor square of ``self``

            EXAMPLES::

                sage: A = HopfAlgebrasWithBasis(QQ).example()
                sage: A.tensor_square()
                An example of Hopf algebra with basis:
                 the group algebra of the Dihedral group of order 6
                 as a permutation group over Rational Field # An example
                 of Hopf algebra with basis: the group algebra of the Dihedral
                 group of order 6 as a permutation group over Rational Field
            """
            return tensor([self, self])

    class ElementMethods:

        __mul__ = sage.categories.coercion_methods.Modules__mul__
        __rmul__ = sage.categories.coercion_methods.Modules__rmul__

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

                    sage: E = CombinatorialFreeModule(ZZ, [1,2,3])
                    sage: F = CombinatorialFreeModule(ZZ, [2,3,4])
                    sage: H = Hom(E, F)
                    sage: H.base_ring()
                    Integer Ring

                This ``base_ring`` method is actually overridden by
                :meth:`sage.structure.category_object.CategoryObject.base_ring`::

                    sage: H.base_ring.__module__

                Here we call it directly::

                    sage: method = H.category().parent_class.base_ring
                    sage: method.__get__(H)()
                    Integer Ring
                """
                return self.domain().base_ring()

            @cached_method
            def zero(self):
                """
                EXAMPLES::

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

                    sage: loads(dumps(f.parent().zero()))
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

                    sage: End(ZZ^3) in Algebras(ZZ)
                    True
                """
                from .magmatic_algebras import MagmaticAlgebras
                return [MagmaticAlgebras(self.base_category().base_ring())]

    class CartesianProducts(CartesianProductsCategory):
        """
        The category of modules constructed as Cartesian products of modules

        This construction gives the direct product of modules. The
        implementation is based on the following resources:

        - http://groups.google.fr/group/sage-devel/browse_thread/thread/35a72b1d0a2fc77a/348f42ae77a66d16#348f42ae77a66d16
        - http://en.wikipedia.org/wiki/Direct_product
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
            def base_ring(self):
                """
                Return the base ring of this Cartesian product.

                EXAMPLES::

                    sage: E = CombinatorialFreeModule(ZZ, [1,2,3])
                    sage: F = CombinatorialFreeModule(ZZ, [2,3,4])
                    sage: C = cartesian_product([E, F]); C
                    Free module generated by {1, 2, 3} over Integer Ring (+)
                    Free module generated by {2, 3, 4} over Integer Ring
                    sage: C.base_ring()
                    Integer Ring
                """
                return self._sets[0].base_ring()

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

