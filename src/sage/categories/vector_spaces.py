r"""
Vector Spaces
"""
#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************
from sage.categories.category import Category
from sage.categories.category_types import Category_module
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.cartesian_product import CartesianProductsCategory
from sage.categories.dual import DualObjectsCategory
from sage.categories.tensor import TensorProductsCategory
from sage.categories.filtered_modules import FilteredModulesCategory
from sage.categories.graded_modules import GradedModulesCategory
from sage.categories.fields import Fields
from sage.categories.modules import Modules
from sage.categories.modules_with_basis import ModulesWithBasis
_Fields = Fields()

class VectorSpaces(Category_module):
    """
    The category of (abstract) vector spaces over a given field

    ??? with an embedding in an ambient vector space ???

    EXAMPLES::

        sage: VectorSpaces(QQ)
        Category of vector spaces over Rational Field
        sage: VectorSpaces(QQ).super_categories()
        [Category of modules over Rational Field]
    """
    @staticmethod
    def __classcall_private__(cls, K, check=True):
        """
        INPUT:

        - `K` -- a field
        - ``check`` -- a boolean (default: True) whether to check that `K` is a field.

        EXAMPLES::

            sage: VectorSpaces(QQ) is VectorSpaces(QQ, check=False)
            True

        By default, it is checked that ``K`` is a field::

            sage: VectorSpaces(ZZ)
            Traceback (most recent call last):
            ...
            ValueError: base must be a field or a subcategory of Fields(); got Integer Ring

        With ``check=False``, the check is disabled, possibly enabling
        incorrect inputs::

            sage: VectorSpaces(ZZ, check=False)
            Category of vector spaces over Integer Ring
        """
        if check:
            if not (K in _Fields or
                    (isinstance(K, Category) and K.is_subcategory(_Fields))):
                raise ValueError("base must be a field or a subcategory of Fields();" +
                                 " got {}".format(K))
        return super(VectorSpaces, cls).__classcall__(cls, K)

    def __init__(self, K):
        """
        EXAMPLES::

            sage: VectorSpaces(QQ)
            Category of vector spaces over Rational Field
            sage: VectorSpaces(ZZ)
            Traceback (most recent call last):
            ...
            ValueError: base must be a field or a subcategory of Fields(); got Integer Ring

        TESTS::

            sage: C = QQ^10      # vector space
            sage: TestSuite(C).run()
            sage: TestSuite(VectorSpaces(QQ)).run()
        """
        Category_module.__init__(self, K)

    def _call_(self, x):
        """
        Try to coerce ``x`` into an object of this category

        EXAMPLES::

            sage: VectorSpaces(QQ)(ZZ^3)
            Vector space of dimension 3 over Rational Field

        TESTS:

        Check whether :trac:`30174` is fixed::

            sage: Q3 = FiniteRankFreeModule(QQ, 3)
            sage: Modules(QQ)(Q3) is Q3
            True

        """
        try:
            V = x.vector_space(self.base_field())
            if V.base_field() != self.base_field():
                V = V.change_ring(self.base_field())
        except (TypeError, AttributeError) as msg:
            raise TypeError("%s\nunable to coerce x (=%s) into %s"%(msg,x,self))
        return V

    def base_field(self):
        """
        Returns the base field over which the vector spaces of this
        category are all defined.

        EXAMPLES::

            sage: VectorSpaces(QQ).base_field()
            Rational Field
        """
        return self.base_ring()

    def super_categories(self):
        """
        EXAMPLES::

            sage: VectorSpaces(QQ).super_categories()
            [Category of modules over Rational Field]
        """
        R = self.base_field()
        return [Modules(R, dispatch = False)]

    def additional_structure(self):
        r"""
        Return  ``None``.

        Indeed, the category of vector spaces defines no additional
        structure: a bimodule morphism between two vector spaces is a
        vector space morphism.

        .. SEEALSO:: :meth:`Category.additional_structure`

        .. TODO:: Should this category be a :class:`CategoryWithAxiom`?

        EXAMPLES::

            sage: VectorSpaces(QQ).additional_structure()
        """
        return None

    class ParentMethods:

        def dimension(self):
            """
            Return the dimension of this vector space.

            EXAMPLES::

                sage: M = FreeModule(FiniteField(19), 100)
                sage: W = M.submodule([M.gen(50)])
                sage: W.dimension()
                1

                sage: M = FiniteRankFreeModule(QQ, 3)
                sage: M.dimension()
                3
                sage: M.tensor_module(1,2).dimension()
                27

            """
            return self.rank()

    class ElementMethods:
        pass

    class WithBasis(CategoryWithAxiom_over_base_ring):

        _call_ = ModulesWithBasis.__dict__["_call_"]

        def is_abelian(self):
            """
            Return whether this category is abelian.

            This is always ``True`` since the base ring is a field.

            EXAMPLES::

                sage: VectorSpaces(QQ).WithBasis().is_abelian()
                True
            """
            return True

        class CartesianProducts(CartesianProductsCategory):
            def extra_super_categories(self):
                r"""
                The category of vector spaces with basis is closed under Cartesian products::

                    sage: C = VectorSpaces(QQ).WithBasis()
                    sage: C.CartesianProducts()
                    Category of Cartesian products of vector spaces with basis over Rational Field
                    sage: C in C.CartesianProducts().super_categories()
                    True
                """
                return [self.base_category()]

        class TensorProducts(TensorProductsCategory):
            def extra_super_categories(self):
                r"""
                The category of vector spaces with basis is closed under tensor products::

                    sage: C = VectorSpaces(QQ).WithBasis()
                    sage: C.TensorProducts()
                    Category of tensor products of vector spaces with basis over Rational Field
                    sage: C in C.TensorProducts().super_categories()
                    True
                """
                return [self.base_category()]

        class Graded(GradedModulesCategory):
            """
            Category of graded vector spaces with basis.
            """
            def example(self, base_ring=None):
                """
                Return an example of a graded vector space with basis,
                as per :meth:`Category.example()
                <sage.categories.category.Category.example>`.

                EXAMPLES::

                    sage: Modules(QQ).WithBasis().Graded().example()
                    An example of a graded module with basis:
                     the free module on partitions over Rational Field
                """
                from sage.categories.examples.graded_modules_with_basis import GradedPartitionModule
                if base_ring is None:
                    base_ring = self.base_ring()
                return GradedPartitionModule(base_ring=base_ring)

        class Filtered(FilteredModulesCategory):
            """
            Category of filtered vector spaces with basis.
            """
            def example(self, base_ring=None):
                """
                Return an example of a graded vector space with basis,
                as per :meth:`Category.example()
                <sage.categories.category.Category.example>`.

                EXAMPLES::

                    sage: Modules(QQ).WithBasis().Graded().example()
                    An example of a graded module with basis:
                     the free module on partitions over Rational Field
                """
                from sage.categories.examples.filtered_modules_with_basis import FilteredPartitionModule
                if base_ring is None:
                    base_ring = self.base_ring()
                return FilteredPartitionModule(base_ring=base_ring)

    class DualObjects(DualObjectsCategory):

        def extra_super_categories(self):
            r"""
            Returns the dual category

            EXAMPLES:

            The category of algebras over the Rational Field is dual
            to the category of coalgebras over the same field::

                sage: C = VectorSpaces(QQ)
                sage: C.dual()
                Category of duals of vector spaces over Rational Field
                sage: C.dual().super_categories() # indirect doctest
                [Category of vector spaces over Rational Field]
            """
            return [self.base_category()]

    class CartesianProducts(CartesianProductsCategory):
        def extra_super_categories(self):
            r"""
            The category of vector spaces is closed under Cartesian products::

                sage: C = VectorSpaces(QQ)
                sage: C.CartesianProducts()
                Category of Cartesian products of vector spaces over Rational Field
                sage: C in C.CartesianProducts().super_categories()
                True
            """
            return [self.base_category()]

    class TensorProducts(TensorProductsCategory):
        def extra_super_categories(self):
            r"""
            The category of vector spaces is closed under tensor products::

                sage: C = VectorSpaces(QQ)
                sage: C.TensorProducts()
                Category of tensor products of vector spaces over Rational Field
                sage: C in C.TensorProducts().super_categories()
                True
            """
            return [self.base_category()]

    class Filtered(FilteredModulesCategory):
        """
        Category of filtered vector spaces.
        """

    class Graded(GradedModulesCategory):
        """
        Category of graded vector spaces.
        """

