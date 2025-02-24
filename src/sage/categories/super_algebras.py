# sage_setup: distribution = sagemath-categories
r"""
Super Algebras
"""
# ****************************************************************************
#  Copyright (C) 2015 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.categories.super_modules import SuperModulesCategory
from sage.categories.signed_tensor import SignedTensorProductsCategory, tensor_signed
from sage.categories.tensor import tensor
from sage.misc.lazy_import import LazyImport
from sage.misc.cachefunc import cached_method


class SuperAlgebras(SuperModulesCategory):
    r"""
    The category of super algebras.

    An `R`-*super algebra* is an `R`-super module `A` endowed with an
    `R`-algebra structure satisfying

    .. MATH::

        A_0 A_0 \subseteq A_0, \qquad
        A_0 A_1 \subseteq A_1, \qquad
        A_1 A_0 \subseteq A_1, \qquad
        A_1 A_1 \subseteq A_0

    and `1 \in A_0`.

    EXAMPLES::

        sage: Algebras(ZZ).Super()
        Category of super algebras over Integer Ring

    TESTS::

        sage: TestSuite(Algebras(ZZ).Super()).run()
    """
    def extra_super_categories(self):
        """
        EXAMPLES::

            sage: Algebras(ZZ).Super().super_categories() # indirect doctest
            [Category of graded algebras over Integer Ring,
             Category of super modules over Integer Ring]
        """
        return [self.base_category().Graded()]

    Supercommutative = LazyImport('sage.categories.supercommutative_algebras',
                                  'SupercommutativeAlgebras')

    class ParentMethods:
        def graded_algebra(self):
            r"""
            Return the associated graded algebra to ``self``.

            .. WARNING::

                Because a super module `M` is naturally `\ZZ / 2 \ZZ`-graded, and
                graded modules have a natural filtration induced by the grading, if
                `M` has a different filtration, then the associated graded module
                `\operatorname{gr} M \neq M`. This is most apparent with super
                algebras, such as the :class:`differential Weyl algebra
                <sage.algebras.weyl_algebra.DifferentialWeylAlgebra>`, and the
                multiplication may not coincide.
            """
            raise NotImplementedError

        def tensor(*parents, **kwargs):
            """
            Return the tensor product of the parents.

            EXAMPLES::

                sage: # needs sage.combinat sage.modules
                sage: A.<x,y,z> = ExteriorAlgebra(ZZ); A.rename('A')
                sage: T = A.tensor(A,A); T
                A # A # A
                sage: T in Algebras(ZZ).Graded().SignedTensorProducts()
                True
                sage: T in Algebras(ZZ).Graded().TensorProducts()
                False
                sage: A.rename(None)

            This also works when the other elements do not have
            a signed tensor product (:issue:`31266`)::

                sage: # needs sage.combinat sage.modules
                sage: a = SteenrodAlgebra(3).an_element()
                sage: M = CombinatorialFreeModule(GF(3), ['s', 't', 'u'])
                sage: s = M.basis()['s']
                sage: tensor([a, s])                                                    # needs sage.rings.finite_rings
                2*Q_1 Q_3 P(2,1) # B['s']
            """
            constructor = kwargs.pop('constructor', tensor_signed)
            try:
                cat = constructor.category_from_parents(parents)
            except AttributeError:
                # Fall back to the usual tensor product if the other parents
                #   do not support signed tensor products
                cat = tensor.category_from_parents(parents)
            return parents[0].__class__.Tensor(parents, category=cat)

    class SubcategoryMethods:
        @cached_method
        def Supercommutative(self):
            r"""
            Return the full subcategory of the supercommutative objects
            of ``self``.

            A super algebra `M` is *supercommutative* if, for all
            homogeneous `x,y\in M`,

            .. MATH::

                x \cdot y = (-1)^{|x||y|} y \cdot x.

            REFERENCES:

            :wikipedia:`Supercommutative_algebra`

            EXAMPLES::

                sage: Algebras(ZZ).Super().Supercommutative()
                Category of supercommutative algebras over Integer Ring
                sage: Algebras(ZZ).Super().WithBasis().Supercommutative()
                Category of supercommutative algebras with basis over Integer Ring
            """
            return self._with_axiom('Supercommutative')

    class SignedTensorProducts(SignedTensorProductsCategory):
        @cached_method
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: Coalgebras(QQ).Graded().SignedTensorProducts().extra_super_categories()
                [Category of graded coalgebras over Rational Field]
                sage: Coalgebras(QQ).Graded().SignedTensorProducts().super_categories()
                [Category of graded coalgebras over Rational Field]

            Meaning: a signed tensor product of coalgebras is a coalgebra
            """
            return [self.base_category()]
