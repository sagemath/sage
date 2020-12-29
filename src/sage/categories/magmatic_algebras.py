r"""
Non-unital non-associative algebras
"""
# ****************************************************************************
#  Copyright (C) 2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.lazy_import import LazyImport
from sage.categories.category import Category
from sage.categories.category_types import Category_over_base_ring
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.magmas import Magmas
from sage.categories.additive_magmas import AdditiveMagmas
from sage.categories.modules import Modules


class MagmaticAlgebras(Category_over_base_ring):
    """
    The category of algebras over a given base ring.

    An algebra over a ring `R` is a module over `R` endowed with a
    bilinear multiplication.

    .. WARNING::

        :class:`MagmaticAlgebras` will eventually replace the current
        :class:`Algebras` for consistency with
        e.g. :wikipedia:`Algebras` which assumes neither associativity
        nor the existence of a unit (see :trac:`15043`).

    EXAMPLES::

        sage: from sage.categories.magmatic_algebras import MagmaticAlgebras
        sage: C = MagmaticAlgebras(ZZ); C
        Category of magmatic algebras over Integer Ring
        sage: C.super_categories()
        [Category of additive commutative additive associative additive unital distributive magmas and additive magmas,
         Category of modules over Integer Ring]

    TESTS::

        sage: TestSuite(C).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: from sage.categories.magmatic_algebras import MagmaticAlgebras
            sage: MagmaticAlgebras(ZZ).super_categories()
            [Category of additive commutative additive associative additive unital distributive magmas and additive magmas, Category of modules over Integer Ring]

            sage: from sage.categories.additive_semigroups import AdditiveSemigroups
            sage: MagmaticAlgebras(ZZ).is_subcategory((AdditiveSemigroups() & Magmas()).Distributive())
            True

        """
        R = self.base_ring()
        # Note: The specifications impose `self` to be a subcategory
        # of the join of its super categories. Here the join is non
        # trivial, since some of the axioms of Modules (like the
        # commutativity of '+') are added to the left hand side.  We
        # might want the infrastructure to take this join for us.
        return Category.join([(Magmas() & AdditiveMagmas()).Distributive(), Modules(R)], as_list=True)

    def additional_structure(self):
        r"""
        Return ``None``.

        Indeed, the category of (magmatic) algebras defines no new
        structure: a morphism of modules and of magmas between two
        (magmatic) algebras is a (magmatic) algebra morphism.

        .. SEEALSO:: :meth:`Category.additional_structure`

        .. TODO::

            This category should be a
            :class:`~sage.categories.category_with_axiom.CategoryWithAxiom`,
            the axiom specifying the compatibility between the magma and
            module structure.

        EXAMPLES::

            sage: from sage.categories.magmatic_algebras import MagmaticAlgebras
            sage: MagmaticAlgebras(ZZ).additional_structure()
        """
        return None

    Associative = LazyImport('sage.categories.associative_algebras', 'AssociativeAlgebras', at_startup=True)
    Unital = LazyImport('sage.categories.unital_algebras', 'UnitalAlgebras', at_startup=True)

    class ParentMethods:

        @abstract_method(optional=True)
        def algebra_generators(self):
            """
            Return a family of generators of this algebra.

            EXAMPLES::

                sage: F = AlgebrasWithBasis(QQ).example(); F
                An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field
                sage: F.algebra_generators()
                Family (B[word: a], B[word: b], B[word: c])
            """

    class WithBasis(CategoryWithAxiom_over_base_ring):

        class ParentMethods:

            def algebra_generators(self):
                r"""
                Return generators for this algebra.

                This default implementation returns the basis of this algebra.

                OUTPUT: a family

                .. SEEALSO::

                    - :meth:`~sage.categories.modules_with_basis.ModulesWithBasis.ParentMethods.basis`
                    - :meth:`MagmaticAlgebras.ParentMethods.algebra_generators`

                EXAMPLES::

                    sage: D4 = DescentAlgebra(QQ, 4).B()
                    sage: D4.algebra_generators()
                    Lazy family (...)_{i in Compositions of 4}

                    sage: R.<x> = ZZ[]
                    sage: P = PartitionAlgebra(1, x, R)
                    sage: P.algebra_generators()
                    Lazy family (Term map from Partition diagrams of order 1 to
                     Partition Algebra of rank 1 with parameter x over Univariate Polynomial Ring in x
                     over Integer Ring(i))_{i in Partition diagrams of order 1}
                """
                return self.basis()

            @abstract_method(optional = True)
            def product_on_basis(self, i, j):
                """
                The product of the algebra on the basis (optional).

                INPUT:

                - ``i``, ``j`` -- the indices of two elements of the
                  basis of ``self``

                Return the product of the two corresponding basis elements
                indexed by ``i`` and ``j``.

                If implemented, :meth:`product` is defined from
                it by bilinearity.

                EXAMPLES::

                    sage: A = AlgebrasWithBasis(QQ).example()
                    sage: Word = A.basis().keys()
                    sage: A.product_on_basis(Word("abc"),Word("cba"))
                    B[word: abccba]
                """

            @lazy_attribute
            def product(self):
                """
                The product of the algebra, as per
                :meth:`Magmas.ParentMethods.product()
                <sage.categories.magmas.Magmas.ParentMethods.product>`

                By default, this is implemented using one of the following
                methods, in the specified order:

                - :meth:`.product_on_basis`
                - :meth:`.product_by_coercion`

                EXAMPLES::

                    sage: A = AlgebrasWithBasis(QQ).example()
                    sage: a, b, c = A.algebra_generators()
                    sage: A.product(a + 2*b, 3*c)
                    3*B[word: ac] + 6*B[word: bc]
                """
                if self.product_on_basis is not NotImplemented:
                    return self._product_from_product_on_basis_multiply
    #                return self._module_morphism(self._module_morphism(self.product_on_basis, position = 0, codomain=self),
    #                                                                                          position = 1)
                elif hasattr(self, "product_by_coercion"):
                    return self.product_by_coercion
                else:
                    return NotImplemented

            # Provides a product using the product_on_basis by calling linear_combination only once
            def _product_from_product_on_basis_multiply( self, left, right ):
                r"""
                Compute the product of two elements by extending
                bilinearly the method :meth:`product_on_basis`.

                EXAMPLES::

                    sage: A = AlgebrasWithBasis(QQ).example(); A
                    An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field
                    sage: (a,b,c) = A.algebra_generators()
                    sage: A._product_from_product_on_basis_multiply(a*b + 2*c, a - b)
                    B[word: aba] - B[word: abb] + 2*B[word: ca] - 2*B[word: cb]

                """
                return self.linear_combination((self.product_on_basis(mon_left, mon_right), coeff_left * coeff_right )
                                                for (mon_left, coeff_left) in left.monomial_coefficients().items()
                                                for (mon_right, coeff_right) in right.monomial_coefficients().items() )

        class FiniteDimensional(CategoryWithAxiom_over_base_ring):
            class ParentMethods:
                @cached_method
                def derivations_basis(self):
                    r"""
                    Return a basis for the Lie algebra of derivations
                    of ``self`` as matrices.

                    A derivation `D` of an algebra is an endomorphism of `A`
                    such that

                    .. MATH::

                        D(ab) = D(a) b + a D(b)

                    for all `a, b \in A`. The set of all derivations
                    form a Lie algebra.

                    EXAMPLES:

                    We construct the Heisenberg Lie algebra as a
                    multiplicative algebra::

                        sage: p_mult = matrix([[0,0,0],[0,0,-1],[0,0,0]])
                        sage: q_mult = matrix([[0,0,1],[0,0,0],[0,0,0]])
                        sage: A = algebras.FiniteDimensional(QQ,
                        ....:          [p_mult,q_mult,matrix(QQ,3,3)], 'p,q,z')
                        sage: A.inject_variables()
                        Defining p, q, z
                        sage: p*q
                        z
                        sage: q*p
                        -z
                        sage: A.derivations_basis()
                        (
                        [1 0 0]  [0 1 0]  [0 0 0]  [0 0 0]  [0 0 0]  [0 0 0]
                        [0 0 0]  [0 0 0]  [1 0 0]  [0 1 0]  [0 0 0]  [0 0 0]
                        [0 0 1], [0 0 0], [0 0 0], [0 0 1], [1 0 0], [0 1 0]
                        )

                    We construct another example using the exterior algebra
                    and verify we obtain a derivation::

                        sage: A = algebras.Exterior(QQ, 1)
                        sage: A.derivations_basis()
                        (
                        [0 0]
                        [0 1]
                        )
                        sage: D = A.module_morphism(matrix=A.derivations_basis()[0], codomain=A)
                        sage: one, e = A.basis()
                        sage: all(D(a*b) == D(a) * b + a * D(b)
                        ....:     for a in A.basis() for b in A.basis())
                        True

                    REFERENCES:

                    :wikipedia:`Derivation_(differential_algebra)`
                    """
                    R = self.base_ring()
                    B = self.basis()
                    keys = list(B.keys())
                    scoeffs = {(j,y,i): c for y in keys for i in keys
                               for j,c in (B[y]*B[i]).monomial_coefficients(copy=False).items()
                              }
                    zero = R.zero()
                    data = {}
                    N = len(keys)
                    for ii,i in enumerate(keys):
                        for ij,j in enumerate(keys):
                            for il,l in enumerate(keys):
                                row = ii + N * ij + N**2 * il
                                for ik,k in enumerate(keys):
                                    data[row,ik+N*il] = (data.get((row,ik+N*il), zero)
                                                         + scoeffs.get((k, i, j), zero))
                                    data[row,ii+N*ik] = (data.get((row,ii+N*ik), zero)
                                                         - scoeffs.get((l, k, j), zero))
                                    data[row,ij+N*ik] = (data.get((row,ij+N*ik), zero)
                                                         - scoeffs.get((l, i, k), zero))
                    from sage.matrix.constructor import matrix
                    mat = matrix(R, data, sparse=True)
                    return tuple([matrix(R, N, N, list(b))
                                  for b in mat.right_kernel().basis()])
