# sage_setup: distribution = sagemath-combinat
# sage.doctest: needs sage.combinat sage.modules sage.groups
r"""
Symmetric Group Algebra
"""
# ****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import itertools

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.algebras_with_basis import AlgebrasWithBasis
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.permutation import Permutation, Permutations, from_permutation_group_element
from sage.combinat.permutation_cython import (left_action_same_n, right_action_same_n)
from sage.combinat.partition import _Partitions, Partitions, Partitions_n
from sage.combinat.tableau import Tableau, StandardTableaux_size, StandardTableaux_shape, StandardTableaux
from sage.combinat.skew_tableau import SkewTableau
from sage.algebras.group_algebra import GroupAlgebra_class
from sage.algebras.cellular_basis import CellularBasis
from sage.categories.weyl_groups import WeylGroups
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.arith.misc import factorial
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.misc.lazy_import import lazy_import
from sage.misc.persist import register_unpickle_override

lazy_import('sage.groups.perm_gps.permgroup_element', 'PermutationGroupElement')


# TODO: Remove this function and replace it with the class
# TODO: Create parents for other bases (such as the seminormal basis)


def SymmetricGroupAlgebra(R, W, category=None):
    r"""
    Return the symmetric group algebra of order ``W`` over the ring ``R``.

    INPUT:

    - ``W`` -- a symmetric group; alternatively an integer `n` can be
      provided, as shorthand for ``Permutations(n)``.
    - ``R`` -- a base ring
    - ``category`` -- a category (default: the category of ``W``)

    This supports several implementations of the symmetric group. At
    this point this has been tested with ``W=Permutations(n)`` and
    ``W=SymmetricGroup(n)``.

    .. WARNING::

        Some features are failing in the latter case, in particular if
        the domain of the symmetric group is not `1,\ldots,n`.

    .. NOTE::

        The brave can also try setting ``W=WeylGroup(['A',n-1])``, but
        little support for this currently exists.

    EXAMPLES::

        sage: QS3 = SymmetricGroupAlgebra(QQ, 3); QS3
        Symmetric group algebra of order 3 over Rational Field
        sage: QS3(1)
        [1, 2, 3]
        sage: QS3(2)
        2*[1, 2, 3]
        sage: basis = [QS3(p) for p in Permutations(3)]
        sage: a = sum(basis); a
        [1, 2, 3] + [1, 3, 2] + [2, 1, 3] + [2, 3, 1] + [3, 1, 2] + [3, 2, 1]
        sage: a^2
        6*[1, 2, 3] + 6*[1, 3, 2] + 6*[2, 1, 3] + 6*[2, 3, 1] + 6*[3, 1, 2] + 6*[3, 2, 1]
        sage: a^2 == 6*a
        True
        sage: b = QS3([3, 1, 2])
        sage: b
        [3, 1, 2]
        sage: b*a
        [1, 2, 3] + [1, 3, 2] + [2, 1, 3] + [2, 3, 1] + [3, 1, 2] + [3, 2, 1]
        sage: b*a == a
        True

    We now construct the symmetric group algebra by providing
    explicitly the underlying group::

        sage: SGA = SymmetricGroupAlgebra(QQ, Permutations(4)); SGA
        Symmetric group algebra of order 4 over Rational Field
        sage: SGA.group()
        Standard permutations of 4
        sage: SGA.an_element()
        [1, 2, 3, 4] + 2*[1, 2, 4, 3] + 3*[1, 3, 2, 4] + [4, 1, 2, 3]

        sage: SGA = SymmetricGroupAlgebra(QQ, SymmetricGroup(4)); SGA
        Symmetric group algebra of order 4 over Rational Field
        sage: SGA.group()
        Symmetric group of order 4! as a permutation group
        sage: SGA.an_element()
        () + (2,3,4) + 2*(1,3)(2,4) + 3*(1,4)(2,3)

        sage: SGA = SymmetricGroupAlgebra(QQ, WeylGroup(["A",3], prefix='s')); SGA
        Symmetric group algebra of order 4 over Rational Field
        sage: SGA.group()
        Weyl Group of type ['A', 3] (as a matrix group acting on the ambient space)
        sage: SGA.an_element()
        s1*s2*s3 + 3*s2*s3*s1*s2 + 2*s3*s1 + 1

    The preferred way to construct the symmetric group algebra is to
    go through the usual ``algebra`` method::

        sage: SGA = Permutations(3).algebra(QQ); SGA
        Symmetric group algebra of order 3 over Rational Field
        sage: SGA.group()
        Standard permutations of 3

        sage: SGA = SymmetricGroup(3).algebra(QQ); SGA
        Symmetric group algebra of order 3 over Rational Field
        sage: SGA.group()
        Symmetric group of order 3! as a permutation group

    The canonical embedding from the symmetric group algebra of order
    `n` to the symmetric group algebra of order `p > n` is available as
    a coercion::

        sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
        sage: QS4 = SymmetricGroupAlgebra(QQ, 4)
        sage: QS4.coerce_map_from(QS3)
        Generic morphism:
          From: Symmetric group algebra of order 3 over Rational Field
          To:   Symmetric group algebra of order 4 over Rational Field

        sage: x3  = QS3([3,1,2]) + 2 * QS3([2,3,1]); x3
        2*[2, 3, 1] + [3, 1, 2]
        sage: QS4(x3)
        2*[2, 3, 1, 4] + [3, 1, 2, 4]

    This allows for mixed expressions::

        sage: x4  = 3 * QS4([3, 1, 4, 2])
        sage: x3 + x4
        2*[2, 3, 1, 4] + [3, 1, 2, 4] + 3*[3, 1, 4, 2]

        sage: QS0 = SymmetricGroupAlgebra(QQ, 0)
        sage: QS1 = SymmetricGroupAlgebra(QQ, 1)
        sage: x0 = QS0([])
        sage: x1 = QS1([1])
        sage: x0 * x1
        [1]
        sage: x3 - (2*x0 + x1) - x4
        -3*[1, 2, 3, 4] + 2*[2, 3, 1, 4] + [3, 1, 2, 4] - 3*[3, 1, 4, 2]

    Caveat: to achieve this, constructing ``SymmetricGroupAlgebra(QQ,
    10)`` currently triggers the construction of all symmetric group
    algebras of smaller order. Is this a feature we really want to have?

    .. WARNING::

        The semantics of multiplication in symmetric group algebras
        with index set ``Permutations(n)`` is determined by the order
        in which permutations are multiplied, which currently defaults
        to "in such a way that multiplication is associative with
        permutations acting on integers from the right", but can be
        changed to the opposite order at runtime by setting the global
        variable ``Permutations.options['mult']`` (see
        :meth:`sage.combinat.permutation.Permutations.options` ).
        On the other hand, the semantics of multiplication in symmetric
        group algebras with index set ``SymmetricGroup(n)`` does not
        depend on this global variable. (This has the awkward
        consequence that the coercions between these two sorts of
        symmetric group algebras do not respect multiplication when
        this global variable is set to ``'r2l'``.)
        In view of this, it is recommended that code not rely on the
        usual multiplication function, but rather use the methods
        :meth:`left_action_product` and :meth:`right_action_product`
        for multiplying permutations (these methods don't depend on the
        setting). See :issue:`14885` for more information.

    We conclude by constructing the algebra of the symmetric group as
    a monoid algebra::

        sage: QS3 = SymmetricGroupAlgebra(QQ, 3, category=Monoids())
        sage: QS3.category()
        Category of finite dimensional cellular monoid algebras over Rational Field
        sage: TestSuite(QS3).run(skip=['_test_construction'])


    TESTS::

        sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
        sage: TestSuite(QS3).run()

        sage: QS3.group()
        Standard permutations of 3

        sage: QS3.one_basis()
        [1, 2, 3]

        sage: p1 = Permutation([1,2,3])
        sage: p2 = Permutation([2,1,3])
        sage: QS3.product_on_basis(p1,p2)
        [2, 1, 3]

        sage: W = WeylGroup(["A",3])
        sage: SGA = SymmetricGroupAlgebra(QQ, W)
        sage: SGA.group() is W
        True
        sage: TestSuite(SGA).run(skip=["_test_cellular", "_test_construction"])
        sage: W = WeylGroup(["A",2])
        sage: SGA = SymmetricGroupAlgebra(QQ, W)
        sage: SGA._test_cellular()

        sage: SG = SymmetricGroupAlgebra(ZZ, 3)
        sage: SG.group().conjugacy_classes_representatives()
        [[1, 2, 3], [2, 1, 3], [2, 3, 1]]

        sage: SGg = SymmetricGroup(3).algebra(ZZ)
        sage: SGg.group().conjugacy_classes_representatives()
        [(), (1,2), (1,2,3)]
    """
    from sage.rings.semirings.non_negative_integer_semiring import NN
    if W in NN:
        W = Permutations(W)
    if category is None:
        category = W.category()
    return SymmetricGroupAlgebra_n(R, W, category.Algebras(R))


class SymmetricGroupAlgebra_n(GroupAlgebra_class):

    def __init__(self, R, W, category):
        """
        TESTS::

            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: TestSuite(QS3).run()

            sage: QS3 in GroupAlgebras(QQ)
            True
            sage: QS3 in FiniteDimensionalAlgebrasWithBasis(QQ)
            True

        Check that :issue:`16926` works::

            sage: S = SymmetricGroup(4)
            sage: SGA = S.algebra(QQ)
            sage: TestSuite(SGA).run(skip='_test_cellular')
            sage: SGA._test_cellular() # long time

        Checking that coercion works between equivalent indexing sets::

            sage: G = SymmetricGroup(4).algebra(QQ)
            sage: S = SymmetricGroupAlgebra(QQ,4)
            sage: S(G.an_element())
            [1, 2, 3, 4] + [1, 3, 4, 2] + 2*[3, 4, 1, 2] + 3*[4, 3, 2, 1]
            sage: G(S.an_element())
            () + 2*(3,4) + 3*(2,3) + (1,4,3,2)

        Checking the recovery of `n`:

            sage: SymmetricGroup(4).algebra(QQ).n
            4
            sage: SymmetricGroup(1).algebra(QQ).n
            1
            sage: SymmetricGroup(0).algebra(QQ).n
            0
            sage: Permutations(4).algebra(QQ).n
            4
            sage: Permutations(1).algebra(QQ).n
            1
            sage: Permutations(0).algebra(QQ).n
            0
            sage: SymmetricGroupAlgebra(QQ, WeylGroup(["A",3])).n
            4
            sage: SymmetricGroupAlgebra(QQ, WeylGroup(["A",1])).n
            2
            sage: SymmetricGroupAlgebra(QQ, WeylGroup(["A",0])).n # todo: not implemented
            1
        """
        if W not in WeylGroups or W.cartan_type().type() != 'A':
            raise ValueError("W (=%s) should be a symmetric group or a nonnegative integer")
        rank = W.cartan_type().rank()
        if rank == 0:   # Ambiguous: n=0 or n=1?
            # The following trick works for both SymmetricGroup(n) and
            # Permutations(n) and it's currently not possible to
            # construct the WeylGroup for n=0
            self.n = W.degree()
        else:
            self.n = W.cartan_type().rank() + 1
        self._idempotent_cache = {}
        category = category.Unital().FiniteDimensional().WithBasis().Cellular()
        GroupAlgebra_class.__init__(self, R, W, prefix='',
                                    latex_prefix='', category=category)

        # Mixin class for extra methods for representations
        from sage.combinat.specht_module import SymmetricGroupRepresentation
        self._representation_mixin_class = SymmetricGroupRepresentation

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: SymmetricGroupAlgebra(QQ, 3)
            Symmetric group algebra of order 3 over Rational Field
        """
        return "Symmetric group algebra of order {} over {}".format(self.n, self.base_ring())

    def _coerce_map_from_(self, S):
        """
        Return ``True`` or a morphism if there exists a coercion from ``S``
        into ``self`` or ``False`` otherwise.

        EXAMPLES:

        Symmetric group algebras::

            sage: SGA4 = SymmetricGroupAlgebra(QQ, 4)
            sage: SGA2 = SymmetricGroupAlgebra(QQ, 2)
            sage: SGA4.has_coerce_map_from(SGA2)
            True
            sage: SGA2Z = SymmetricGroupAlgebra(ZZ, 2)
            sage: SGA4.has_coerce_map_from(SGA2Z)
            True
            sage: p = Permutation([2,1])
            sage: SGA4(-3*SGA2Z.monomial(p))
            -3*[2, 1, 3, 4]

        Descent algebras::

            sage: DA = DescentAlgebra(QQ, 4)
            sage: SGA4 = SymmetricGroupAlgebra(QQ, 4)
            sage: SGA4.has_coerce_map_from(DA.D())
            True
            sage: SGA4.has_coerce_map_from(DA.B())
            True
            sage: SGA4.has_coerce_map_from(DA.I())
            True
            sage: x = DA.B()[4]
            sage: SGA4(x)
            [1, 2, 3, 4]

            sage: DAB = DescentAlgebra(ZZ,2).B()
            sage: SGA4.has_coerce_map_from(DAB)
            True
            sage: SGA4(DAB[2])
            [1, 2, 3, 4]

            sage: QSG4 = SymmetricGroup(4).algebra(ZZ)
            sage: DAB = DescentAlgebra(ZZ,4).B()
            sage: QSG4(DAB[1,2,1])
            () + (3,4) + (2,3,4) + (1,2) + (1,2)(3,4) + (1,2,3,4)
             + (1,3,2) + (1,3,4,2) + (1,3,4) + (1,4,3,2) + (1,4,2) + (1,4)
        """
        # Symmetric group algebras of smaller rank
        if (isinstance(S, SymmetricGroupAlgebra_n) and S.n <= self.n and
                self.base_ring().has_coerce_map_from(S.base_ring())):
            return S.canonical_embedding(self)

        # Descent algebras
        from sage.combinat.descent_algebra import DescentAlgebra
        # TODO: A better way to handle all of the bases
        if isinstance(S, (DescentAlgebra.D, DescentAlgebra.B, DescentAlgebra.I)):
            # Same rank and base ring, just the natural morphism
            if (S.realization_of()._n == self.n and
                    self.base_ring() == S.base_ring() and
                    self._indices == Permutations(self.n)):
                return S.to_symmetric_group_algebra
            # Otherwise compose with the canonical embedding in order to ensure
            # that the right base ring and the right index set are being used.
            # Slightly hacky!
            if (S.realization_of()._n <= self.n and
                    self.base_ring().has_coerce_map_from(S.base_ring())):
                phi = S.to_symmetric_group_algebra
                return phi.codomain().canonical_embedding(self) * phi

        return super()._coerce_map_from_(S)

    def _element_constructor_(self, x):
        """
        Construct an element of ``self``.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 4)
            sage: G = SymmetricGroup(3)
            sage: p = Permutation((1,2))
            sage: S(p)
            [2, 1, 3, 4]
            sage: S(G(p))
            [2, 1, 3, 4]
            sage: S(p) == S(G(p))
            True
        """
        if isinstance(x, Permutation):
            return self.monomial_from_smaller_permutation(x)
        if isinstance(x, PermutationGroupElement):
            return self.monomial_from_smaller_permutation(
                from_permutation_group_element(x))

        return super()._element_constructor_(x)

    def _sibling(self, n):
        r"""
        Return the sibling group algebra of order `n`.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, Permutations(3))._sibling(4); SGA
            Symmetric group algebra of order 4 over Rational Field
            sage: SGA.group()
            Standard permutations of 4

            sage: SGA = SymmetricGroupAlgebra(QQ, SymmetricGroup(3))._sibling(4); SGA
            Symmetric group algebra of order 4 over Rational Field
            sage: SGA.group()
            Symmetric group of order 4! as a permutation group

            sage: SGA = SymmetricGroupAlgebra(QQ, WeylGroup(["A",2]))._sibling(4); SGA
            Traceback (most recent call last):
            ...
            NotImplementedError: Constructing the sibling algebra of a different order
            only implemented for PermutationGroup and SymmetricGroup
        """
        try:
            W = self.basis().keys().__class__(n)
        except (AttributeError, TypeError, ValueError):
            raise NotImplementedError("Constructing the sibling algebra of a different order "
                                      "only implemented for PermutationGroup and SymmetricGroup")
        return SymmetricGroupAlgebra(self.base_ring(), W)

    # _repr_ customization: output the basis element indexed by [1,2,3] as [1,2,3]
    _repr_option_bracket = False

    def left_action_product(self, left, right):
        """
        Return the product of two elements ``left`` and ``right`` of
        ``self``, where multiplication is defined in such a way that
        for two permutations `p` and `q`, the product `pq` is the
        permutation obtained by first applying `q` and then applying
        `p`. This definition of multiplication is tailored to make
        multiplication of permutations associative with their action on
        numbers if permutations are to act on numbers from the left.

        EXAMPLES::

            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: p1 = Permutation([2, 1, 3])
            sage: p2 = Permutation([3, 1, 2])
            sage: QS3.left_action_product(QS3(p1), QS3(p2))
            [3, 2, 1]
            sage: x = QS3([1, 2, 3]) - 2*QS3([1, 3, 2])
            sage: y = 1/2 * QS3([3, 1, 2]) + 3*QS3([1, 2, 3])
            sage: QS3.left_action_product(x, y)
            3*[1, 2, 3] - 6*[1, 3, 2] - [2, 1, 3] + 1/2*[3, 1, 2]
            sage: QS3.left_action_product(0, x)
            0

        The method coerces its input into the algebra ``self``::

            sage: QS4 = SymmetricGroupAlgebra(QQ, 4)
            sage: QS4.left_action_product(QS3([1, 2, 3]), QS3([2, 1, 3]))
            [2, 1, 3, 4]
            sage: QS4.left_action_product(1, Permutation([4, 1, 2, 3]))
            [4, 1, 2, 3]

        TESTS::

            sage: QS4 = SymmetricGroup(4).algebra(QQ)
            sage: QS4.left_action_product(QS4((1,2)), QS4((2,3)))
            (1,2,3)
            sage: QS4.left_action_product(1, QS4((1,2)))
            (1,2)

        .. WARNING::

            Note that coercion presently works from permutations of ``n``
            into the ``n``-th symmetric group algebra, and also from all
            smaller symmetric group algebras into the ``n``-th symmetric
            group algebra, but not from permutations of integers smaller
            than ``n`` into the ``n``-th symmetric group algebra.
        """
        a = self(left)
        b = self(right)
        if not isinstance(self._indices, Permutations):
            return b * a
        P = Permutations(self.n)
        return self.sum_of_terms([(P(left_action_same_n(p._list, q._list)), x * y)
                                  for (p, x) in a for (q, y) in b])
        # Why did we use left_action_same_n instead of
        # left_action_product?
        # Because having cast a and b into self, we already know that
        # p and q are permutations of the same number of elements,
        # and thus we don't need to waste our time on the input
        # sanitizing of left_action_product.

    def right_action_product(self, left, right):
        """
        Return the product of two elements ``left`` and ``right`` of
        ``self``, where multiplication is defined in such a way that
        for two permutations `p` and `q`, the product `pq` is the
        permutation obtained by first applying `p` and then applying
        `q`. This definition of multiplication is tailored to make
        multiplication of permutations associative with their action on
        numbers if permutations are to act on numbers from the right.

        EXAMPLES::

            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: p1 = Permutation([2, 1, 3])
            sage: p2 = Permutation([3, 1, 2])
            sage: QS3.right_action_product(QS3(p1), QS3(p2))
            [1, 3, 2]
            sage: x = QS3([1, 2, 3]) - 2*QS3([1, 3, 2])
            sage: y = 1/2 * QS3([3, 1, 2]) + 3*QS3([1, 2, 3])
            sage: QS3.right_action_product(x, y)
            3*[1, 2, 3] - 6*[1, 3, 2] + 1/2*[3, 1, 2] - [3, 2, 1]
            sage: QS3.right_action_product(0, x)
            0

        The method coerces its input into the algebra ``self``::

            sage: QS4 = SymmetricGroupAlgebra(QQ, 4)
            sage: QS4.right_action_product(QS3([1, 2, 3]), QS3([2, 1, 3]))
            [2, 1, 3, 4]
            sage: QS4.right_action_product(1, Permutation([4, 1, 2, 3]))
            [4, 1, 2, 3]

        TESTS::

            sage: QS4 = SymmetricGroup(4).algebra(QQ)
            sage: QS4.right_action_product(QS4((1,2)), QS4((2,3)))
            (1,3,2)
            sage: QS4.right_action_product(1, QS4((1,2)))
            (1,2)

        .. WARNING::

            Note that coercion presently works from permutations of ``n``
            into the ``n``-th symmetric group algebra, and also from all
            smaller symmetric group algebras into the ``n``-th symmetric
            group algebra, but not from permutations of integers smaller
            than ``n`` into the ``n``-th symmetric group algebra.
        """
        a = self(left)
        b = self(right)
        if not isinstance(self._indices, Permutations):
            return a * b
        P = Permutations(self.n)
        return self.sum_of_terms([(P(right_action_same_n(p._list, q._list)), x * y)
                                  for (p, x) in a for (q, y) in b])
        # Why did we use right_action_same_n instead of
        # right_action_product?
        # Because having cast a and b into self, we already know that
        # p and q are permutations of the same number of elements,
        # and thus we don't need to waste our time on the input
        # sanitizing of right_action_product.

    def canonical_embedding(self, other):
        r"""
        Return the canonical coercion of ``self`` into a symmetric
        group algebra ``other``.

        INPUT:

        - ``other`` -- a symmetric group algebra with order `p`
          satisfying `p \geq n`, where `n` is the order of ``self``,
          over a ground ring into which the ground ring of ``self``
          coerces.

        EXAMPLES::

            sage: QS2 = SymmetricGroupAlgebra(QQ, 2)
            sage: QS4 = SymmetricGroupAlgebra(QQ, 4)
            sage: phi = QS2.canonical_embedding(QS4); phi
            Generic morphism:
              From: Symmetric group algebra of order 2 over Rational Field
              To:   Symmetric group algebra of order 4 over Rational Field

            sage: x = QS2([2,1]) + 2 * QS2([1,2])
            sage: phi(x)
            2*[1, 2, 3, 4] + [2, 1, 3, 4]

            sage: loads(dumps(phi))
            Generic morphism:
              From: Symmetric group algebra of order 2 over Rational Field
              To:   Symmetric group algebra of order 4 over Rational Field

            sage: ZS2 = SymmetricGroupAlgebra(ZZ, 2)
            sage: phi = ZS2.canonical_embedding(QS4); phi
            Generic morphism:
              From: Symmetric group algebra of order 2 over Integer Ring
              To:   Symmetric group algebra of order 4 over Rational Field

            sage: phi = ZS2.canonical_embedding(QS2); phi
            Generic morphism:
              From: Symmetric group algebra of order 2 over Integer Ring
              To:   Symmetric group algebra of order 2 over Rational Field

            sage: QS4.canonical_embedding(QS2)
            Traceback (most recent call last):
            ...
            ValueError: There is no canonical embedding from Symmetric group
             algebra of order 2 over Rational Field to Symmetric group
             algebra of order 4 over Rational Field

            sage: QS4g = SymmetricGroup(4).algebra(QQ)
            sage: QS4.canonical_embedding(QS4g)(QS4([1,3,2,4]))
            (2,3)
            sage: QS4g.canonical_embedding(QS4)(QS4g((2,3)))
            [1, 3, 2, 4]
            sage: ZS2.canonical_embedding(QS4g)(ZS2([2,1]))
            (1,2)
            sage: ZS2g = SymmetricGroup(2).algebra(ZZ)
            sage: ZS2g.canonical_embedding(QS4)(ZS2g((1,2)))
            [2, 1, 3, 4]
        """
        if not isinstance(other, SymmetricGroupAlgebra_n) or self.n > other.n:
            raise ValueError("There is no canonical embedding from {0} to {1}".format(other, self))
        return self.module_morphism(other.monomial_from_smaller_permutation, codomain=other)  # category = self.category() (currently broken)

    def monomial_from_smaller_permutation(self, permutation):
        """
        Convert ``permutation`` into a permutation, possibly extending it
        to the appropriate size, and return the corresponding basis
        element of ``self``.

        EXAMPLES::

            sage: QS5 = SymmetricGroupAlgebra(QQ, 5)
            sage: QS5.monomial_from_smaller_permutation([])
            [1, 2, 3, 4, 5]
            sage: QS5.monomial_from_smaller_permutation(Permutation([3,1,2]))
            [3, 1, 2, 4, 5]
            sage: QS5.monomial_from_smaller_permutation([5,3,4,1,2])
            [5, 3, 4, 1, 2]
            sage: QS5.monomial_from_smaller_permutation(SymmetricGroup(2)((1,2)))
            [2, 1, 3, 4, 5]

            sage: QS5g = SymmetricGroup(5).algebra(QQ)
            sage: QS5g.monomial_from_smaller_permutation([2,1])
            (1,2)

        TESTS::

            sage: QS5.monomial_from_smaller_permutation([5,3,4,1,2]).parent()
            Symmetric group algebra of order 5 over Rational Field
        """
        P = self.basis().keys()
        return self.monomial(P(permutation))

    def antipode(self, x):
        r"""
        Return the image of the element ``x`` of ``self`` under the
        antipode of the Hopf algebra ``self`` (where the
        comultiplication is the usual one on a group algebra).

        Explicitly, this is obtained by replacing each permutation
        `\sigma` by `\sigma^{-1}` in ``x`` while keeping all
        coefficients as they are.

        EXAMPLES::

            sage: QS4 = SymmetricGroupAlgebra(QQ, 4)
            sage: QS4.antipode(2 * QS4([1, 3, 4, 2]) - 1/2 * QS4([1, 4, 2, 3]))
            -1/2*[1, 3, 4, 2] + 2*[1, 4, 2, 3]
            sage: all( QS4.antipode(QS4(p)) == QS4(p.inverse())
            ....:      for p in Permutations(4) )
            True

            sage: ZS3 = SymmetricGroupAlgebra(ZZ, 3)
            sage: ZS3.antipode(ZS3.zero())
            0
            sage: ZS3.antipode(-ZS3(Permutation([2, 3, 1])))
            -[3, 1, 2]
        """
        return self.sum_of_terms([(p.inverse(), coeff) for
                                  (p, coeff) in self(x)],
                                 distinct=True)

    @cached_method
    def cell_poset(self):
        """
        Return the cell poset of ``self``.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 4)
            sage: S.cell_poset()
            Finite poset containing 5 elements
        """
        from sage.combinat.posets.posets import Poset
        return Poset([Partitions_n(self.n), lambda x, y: y.dominates(x)])

    def cell_module_indices(self, la):
        r"""
        Return the indices of the cell module of ``self``
        indexed by ``la`` .

        This is the finite set `M(\lambda)`.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 4)
            sage: S.cell_module_indices([3,1])
            Standard tableaux of shape [3, 1]
        """
        return StandardTableaux(la)

    def _from_cellular_index(self, x):
        r"""
        Return the image in ``self`` from the index of the
        cellular basis ``x``.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: C = S.cellular_basis()
            sage: [S._from_cellular_index(i) for i in C.basis().keys()]
            [1/6*[1, 2, 3] - 1/6*[1, 3, 2] - 1/6*[2, 1, 3] + 1/6*[2, 3, 1]
                 + 1/6*[3, 1, 2] - 1/6*[3, 2, 1],
             1/3*[1, 2, 3] + 1/6*[1, 3, 2] - 1/3*[2, 1, 3] - 1/6*[2, 3, 1]
                 - 1/6*[3, 1, 2] + 1/6*[3, 2, 1],
             1/3*[1, 3, 2] + 1/3*[2, 3, 1] - 1/3*[3, 1, 2] - 1/3*[3, 2, 1],
             1/4*[1, 3, 2] - 1/4*[2, 3, 1] + 1/4*[3, 1, 2] - 1/4*[3, 2, 1],
             1/3*[1, 2, 3] - 1/6*[1, 3, 2] + 1/3*[2, 1, 3] - 1/6*[2, 3, 1]
                 - 1/6*[3, 1, 2] - 1/6*[3, 2, 1],
             1/6*[1, 2, 3] + 1/6*[1, 3, 2] + 1/6*[2, 1, 3]
                 + 1/6*[2, 3, 1] + 1/6*[3, 1, 2] + 1/6*[3, 2, 1]]

            sage: S = SymmetricGroupAlgebra(GF(3), 3)
            sage: C = S.cellular_basis()
            sage: [S._from_cellular_index(i) for i in C.basis().keys()]
            [[1, 2, 3] + [1, 3, 2] + [2, 1, 3] + [2, 3, 1] + [3, 1, 2] + [3, 2, 1],
             [1, 2, 3] + [2, 1, 3],
             [1, 3, 2] + [3, 1, 2],
             [1, 3, 2] + [2, 3, 1],
             [1, 2, 3] + [3, 2, 1],
             [1, 2, 3]]
            sage: TestSuite(C).run()
        """
        if ~factorial(self.n) not in self.base_ring():
            func = self.murphy_basis_element
        else:
            func = self.epsilon_ik

        SGA = SymmetricGroupAlgebra(self.base_ring(), self.n)
        P = self.basis().keys()
        if SGA.basis().keys() is P:  # Indexed by permutations
            return func(x[1], x[2])
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        if P == SymmetricGroup(self.n):
            return func(x[1], x[2])
        ret = func(x[1], x[2], mult='r2l')
        try:
            return self(ret)
        except TypeError:
            P = self.basis().keys()
            return self._from_dict({P(i.to_matrix()): c for i, c in ret},
                                   remove_zeros=False)

    def cell_module(self, la, **kwds):
        """
        Return the cell module indexed by ``la``.

        EXAMPLES::

            sage: S = SymmetricGroupAlgebra(QQ, 3)
            sage: M = S.cell_module(Partition([2,1])); M
            Cell module indexed by [2, 1] of Cellular basis of
             Symmetric group algebra of order 3 over Rational Field

        We check that the input ``la`` is standardized::

            sage: N = S.cell_module([2,1])
            sage: M is N
            True
        """
        la = _Partitions(la)
        kwds['bracket'] = kwds.get('bracket', False)
        return super().cell_module(la, **kwds)

    def retract_plain(self, f, m):
        r"""
        Return the plain retract of the element `f \in R S_n`
        to `R S_m`, where `m \leq n` (and where `R S_n` is ``self``).

        If `m` is a nonnegative integer less or equal to `n`, then the
        plain retract from `S_n` to `S_m` is defined as an `R`-linear
        map `S_n \to S_m` which sends every permutation `p \in S_n`
        to

        .. MATH::

            \begin{cases} \mbox{pret}(p) &\mbox{if } \mbox{pret}(p)\mbox{ is defined;} \\
            0 & \mbox{otherwise} \end{cases}.

        Here `\mbox{pret}(p)` denotes the plain retract of the
        permutation `p` to `S_m`, which is defined in
        :meth:`~sage.combinat.permutation.Permutation.retract_plain`.

        EXAMPLES::

            sage: SGA3 = SymmetricGroupAlgebra(QQ, 3)
            sage: SGA3.retract_plain(2*SGA3([1,2,3]) - 4*SGA3([2,1,3]) + 7*SGA3([1,3,2]), 2)
            2*[1, 2] - 4*[2, 1]
            sage: SGA3.retract_plain(2*SGA3([1,3,2]) - 5*SGA3([2,3,1]), 2)
            0

            sage: SGA5 = SymmetricGroupAlgebra(QQ, 5)
            sage: SGA5.retract_plain(8*SGA5([1,4,2,5,3]) - 6*SGA5([1,3,2,5,4]) + 11*SGA5([3,2,1,4,5]), 4)
            11*[3, 2, 1, 4]
            sage: SGA5.retract_plain(8*SGA5([1,4,2,5,3]) - 6*SGA5([1,3,2,5,4]) + 11*SGA5([3,2,1,4,5]), 3)
            11*[3, 2, 1]
            sage: SGA5.retract_plain(8*SGA5([1,4,2,5,3]) - 6*SGA5([1,3,2,5,4]) + 11*SGA5([3,2,1,4,5]), 2)
            0
            sage: SGA5.retract_plain(8*SGA5([1,4,2,5,3]) - 6*SGA5([1,3,2,5,4]) + 11*SGA5([3,2,1,4,5]), 1)
            0

            sage: SGA5.retract_plain(8*SGA5([1,2,3,4,5]) - 6*SGA5([1,3,2,4,5]), 3)
            8*[1, 2, 3] - 6*[1, 3, 2]
            sage: SGA5.retract_plain(8*SGA5([1,2,3,4,5]) - 6*SGA5([1,3,2,4,5]), 1)
            8*[1]
            sage: SGA5.retract_plain(8*SGA5([1,2,3,4,5]) - 6*SGA5([1,3,2,4,5]), 0)
            8*[]

        TESTS:

        Check this works with other indexing sets::

            sage: G = SymmetricGroup(4).algebra(QQ)
            sage: G.retract_plain(G.an_element(), 3)
            ()

        .. SEEALSO::

            :meth:`retract_direct_product`, :meth:`retract_okounkov_vershik`
        """
        RSm = self._sibling(m)
        I = RSm.group()
        pairs = []
        P = Permutations(self.n)
        for (p, coeff) in f.monomial_coefficients().items():
            p_ret = P(p).retract_plain(m)
            if p_ret is not None:
                pairs.append((I(p_ret), coeff))
        return RSm.sum_of_terms(pairs, distinct=True)

    def retract_direct_product(self, f, m):
        r"""
        Return the direct-product retract of the element `f \in R S_n`
        to `R S_m`, where `m \leq n` (and where `R S_n` is ``self``).

        If `m` is a nonnegative integer less or equal to `n`, then the
        direct-product retract from `S_n` to `S_m` is defined as an
        `R`-linear map `S_n \to S_m` which sends every permutation
        `p \in S_n` to

        .. MATH::

            \begin{cases} \mbox{dret}(p) &\mbox{if } \mbox{dret}(p)\mbox{ is defined;} \\
            0 & \mbox{otherwise} \end{cases}.

        Here `\mbox{dret}(p)` denotes the direct-product retract of the
        permutation `p` to `S_m`, which is defined in
        :meth:`~sage.combinat.permutation.Permutation.retract_direct_product`.

        EXAMPLES::

            sage: SGA3 = SymmetricGroupAlgebra(QQ, 3)
            sage: SGA3.retract_direct_product(2*SGA3([1,2,3]) - 4*SGA3([2,1,3]) + 7*SGA3([1,3,2]), 2)
            2*[1, 2] - 4*[2, 1]
            sage: SGA3.retract_direct_product(2*SGA3([1,3,2]) - 5*SGA3([2,3,1]), 2)
            0

            sage: SGA5 = SymmetricGroupAlgebra(QQ, 5)
            sage: SGA5.retract_direct_product(8*SGA5([1,4,2,5,3]) - 6*SGA5([1,3,2,5,4]) + 11*SGA5([3,2,1,4,5]), 4)
            11*[3, 2, 1, 4]
            sage: SGA5.retract_direct_product(8*SGA5([1,4,2,5,3]) - 6*SGA5([1,3,2,5,4]) + 11*SGA5([3,2,1,4,5]), 3)
            -6*[1, 3, 2] + 11*[3, 2, 1]
            sage: SGA5.retract_direct_product(8*SGA5([1,4,2,5,3]) - 6*SGA5([1,3,2,5,4]) + 11*SGA5([3,2,1,4,5]), 2)
            0
            sage: SGA5.retract_direct_product(8*SGA5([1,4,2,5,3]) - 6*SGA5([1,3,2,5,4]) + 11*SGA5([3,2,1,4,5]), 1)
            2*[1]

            sage: SGA5.retract_direct_product(8*SGA5([1,2,3,4,5]) - 6*SGA5([1,3,2,4,5]), 3)
            8*[1, 2, 3] - 6*[1, 3, 2]
            sage: SGA5.retract_direct_product(8*SGA5([1,2,3,4,5]) - 6*SGA5([1,3,2,4,5]), 1)
            2*[1]
            sage: SGA5.retract_direct_product(8*SGA5([1,2,3,4,5]) - 6*SGA5([1,3,2,4,5]), 0)
            2*[]

        TESTS:

        Check this works with other indexing sets::

            sage: G = SymmetricGroup(4).algebra(QQ)
            sage: G.retract_direct_product(G.an_element(), 3)
            ()

        .. SEEALSO::

            :meth:`retract_plain`, :meth:`retract_okounkov_vershik`
        """
        RSm = self._sibling(m)
        I = RSm.group()
        dct = {}
        P = Permutations(self.n)
        for (p, coeff) in f.monomial_coefficients().items():
            p_ret = P(p).retract_direct_product(m)
            if p_ret is not None:
                p_ret = I(p_ret)
                if p_ret not in dct:
                    dct[p_ret] = coeff
                else:
                    dct[p_ret] += coeff
        return RSm._from_dict(dct)

    def retract_okounkov_vershik(self, f, m):
        r"""
        Return the Okounkov-Vershik retract of the element `f \in R S_n`
        to `R S_m`, where `m \leq n` (and where `R S_n` is ``self``).

        If `m` is a nonnegative integer less or equal to `n`, then the
        Okounkov-Vershik retract from `S_n` to `S_m` is defined as an
        `R`-linear map `S_n \to S_m` which sends every permutation
        `p \in S_n` to the Okounkov-Vershik retract of the permutation
        `p` to `S_m`, which is defined in
        :meth:`~sage.combinat.permutation.Permutation.retract_okounkov_vershik`.

        EXAMPLES::

            sage: SGA3 = SymmetricGroupAlgebra(QQ, 3)
            sage: SGA3.retract_okounkov_vershik(2*SGA3([1,2,3]) - 4*SGA3([2,1,3]) + 7*SGA3([1,3,2]), 2)
            9*[1, 2] - 4*[2, 1]
            sage: SGA3.retract_okounkov_vershik(2*SGA3([1,3,2]) - 5*SGA3([2,3,1]), 2)
            2*[1, 2] - 5*[2, 1]

            sage: SGA5 = SymmetricGroupAlgebra(QQ, 5)
            sage: SGA5.retract_okounkov_vershik(8*SGA5([1,4,2,5,3]) - 6*SGA5([1,3,2,5,4]) + 11*SGA5([3,2,1,4,5]), 4)
            -6*[1, 3, 2, 4] + 8*[1, 4, 2, 3] + 11*[3, 2, 1, 4]
            sage: SGA5.retract_okounkov_vershik(8*SGA5([1,4,2,5,3]) - 6*SGA5([1,3,2,5,4]) + 11*SGA5([3,2,1,4,5]), 3)
            2*[1, 3, 2] + 11*[3, 2, 1]
            sage: SGA5.retract_okounkov_vershik(8*SGA5([1,4,2,5,3]) - 6*SGA5([1,3,2,5,4]) + 11*SGA5([3,2,1,4,5]), 2)
            13*[1, 2]
            sage: SGA5.retract_okounkov_vershik(8*SGA5([1,4,2,5,3]) - 6*SGA5([1,3,2,5,4]) + 11*SGA5([3,2,1,4,5]), 1)
            13*[1]

            sage: SGA5.retract_okounkov_vershik(8*SGA5([1,2,3,4,5]) - 6*SGA5([1,3,2,4,5]), 3)
            8*[1, 2, 3] - 6*[1, 3, 2]
            sage: SGA5.retract_okounkov_vershik(8*SGA5([1,2,3,4,5]) - 6*SGA5([1,3,2,4,5]), 1)
            2*[1]
            sage: SGA5.retract_okounkov_vershik(8*SGA5([1,2,3,4,5]) - 6*SGA5([1,3,2,4,5]), 0)
            2*[]

        TESTS:

        Check this works with other indexing sets::

            sage: G = SymmetricGroup(4).algebra(QQ)
            sage: G.retract_okounkov_vershik(G.an_element(), 3)
            () + 4*(2,3) + 2*(1,3)

        .. SEEALSO::

            :meth:`retract_plain`, :meth:`retract_direct_product`
        """
        RSm = self._sibling(m)
        I = RSm.group()
        dct = {}
        P = Permutations(self.n)
        for (p, coeff) in f.monomial_coefficients().items():
            p_ret = I(P(p).retract_okounkov_vershik(m))
            if p_ret not in dct:
                dct[p_ret] = coeff
            else:
                dct[p_ret] += coeff
        return RSm._from_dict(dct)

    def central_orthogonal_idempotents(self):
        r"""
        Return a maximal list of central orthogonal idempotents for ``self``.

        This method does not require that ``self`` be semisimple, relying
        on Nakayama's Conjecture whenever ``self.base_ring()`` has
        positive characteristic.

        EXAMPLES::

            sage: QS3 = SymmetricGroupAlgebra(QQ,3)
            sage: a = QS3.central_orthogonal_idempotents()
            sage: a[0]  # [3]
            1/6*[1, 2, 3] + 1/6*[1, 3, 2] + 1/6*[2, 1, 3] + 1/6*[2, 3, 1]
             + 1/6*[3, 1, 2] + 1/6*[3, 2, 1]
            sage: a[1]  # [2, 1]
            2/3*[1, 2, 3] - 1/3*[2, 3, 1] - 1/3*[3, 1, 2]

        TESTS:

        Check this works with other indexing sets::

            sage: G = SymmetricGroup(3).algebra(QQ)
            sage: a = G.central_orthogonal_idempotents()
            sage: a[0]
            1/6*() + 1/6*(2,3) + 1/6*(1,2) + 1/6*(1,2,3) + 1/6*(1,3,2) + 1/6*(1,3)
            sage: a[1]
            2/3*() - 1/3*(1,2,3) - 1/3*(1,3,2)

            sage: G = SymmetricGroup(3).algebra(GF(2))
            sage: a = G.central_orthogonal_idempotents()
            sage: a[0]
            (1,2,3) + (1,3,2)
            sage: a[1]
            () + (1,2,3) + (1,3,2)

        Check this works in positive characteristic::

            sage: def test_n_with_primes(n, primes):
            ....:     Sn = {p:SymmetricGroupAlgebra(GF(p), n) for p in primes}
            ....:     for p in primes:
            ....:         idems = Sn[p].central_orthogonal_idempotents()
            ....:         tst = [sum(idems)==Sn[p].one()]
            ....:         for i in range(len(idems)-1):
            ....:             e = idems[i]
            ....:             for j in range(i, len(idems)):
            ....:                 f = idems[j]
            ....:                 if i == j:
            ....:                     tst.append(e*e == e)
            ....:                 else:
            ....:                     tst.append(e*f == 0)
            ....:         print("{0} blocks for p={1} ... {2}".format( len(idems), p, all(tst) ))
            sage: test_n_with_primes(5, [2,3,5,7])  # long time
            2 blocks for p=2 ... True
            3 blocks for p=3 ... True
            3 blocks for p=5 ... True
            7 blocks for p=7 ... True

        .. SEEALSO::

            - :meth:`central_orthogonal_idempotent`
        """
        return [self.central_orthogonal_idempotent(key)
                for key in sorted(self._blocks_dictionary, reverse=True)]

    def central_orthogonal_idempotent(self, la, block=True):
        r"""
        Return the central idempotent for the symmetric group of order `n`
        corresponding to the indecomposable block to which the partition
        ``la`` is associated.

        If ``self.base_ring()`` contains `\QQ`, this corresponds to the
        classical central idempotent corresponding to the irreducible
        representation indexed by ``la``.

        Alternatively, if ``self.base_ring()`` has characteristic `p > 0`,
        then Theorem 2.8 in [Mur1983]_ provides that ``la`` is associated
        to an idempotent `f_\mu`, where `\mu` is the `p`-core of ``la``.
        This `f_\mu` is a sum of classical idempotents,

        .. MATH::

            f_\mu = \sum_{c(\lambda)=\mu} e_\lambda,

        where the sum ranges over the partitions `\lambda` of `n` with
        `p`-core equal to `\mu`.

        INPUT:

        - ``la`` -- a partition of ``self.n`` or a
          ``self.base_ring().characteristic()``-core of such
          a partition

        - ``block`` -- boolean (default: ``True``); when ``False``,
          this returns the classical idempotent associated to ``la``
          (defined over `\QQ`)

        OUTPUT:

        If ``block=False`` and the corresponding coefficients are
        not defined over ``self.base_ring()``, then return ``None``.
        Otherwise return an element of ``self``.

        EXAMPLES:

        Asking for block idempotents in any characteristic, by
        passing a partition of ``self.n``::

            sage: S0 = SymmetricGroup(4).algebra(QQ)
            sage: S2 = SymmetricGroup(4).algebra(GF(2))
            sage: S3 = SymmetricGroup(4).algebra(GF(3))
            sage: S0.central_orthogonal_idempotent([2,1,1])
            3/8*() - 1/8*(3,4) - 1/8*(2,3) - 1/8*(2,4) - 1/8*(1,2)
             - 1/8*(1,2)(3,4) + 1/8*(1,2,3,4) + 1/8*(1,2,4,3)
             + 1/8*(1,3,4,2) - 1/8*(1,3) - 1/8*(1,3)(2,4)
             + 1/8*(1,3,2,4) + 1/8*(1,4,3,2) - 1/8*(1,4)
             + 1/8*(1,4,2,3) - 1/8*(1,4)(2,3)
            sage: S2.central_orthogonal_idempotent([2,1,1])
            ()
            sage: idem = S3.central_orthogonal_idempotent([4]); idem
             () + (1,2)(3,4) + (1,3)(2,4) + (1,4)(2,3)
            sage: idem == S3.central_orthogonal_idempotent([1,1,1,1])
            True
            sage: S3.central_orthogonal_idempotent([2,2])
            () + (1,2)(3,4) + (1,3)(2,4) + (1,4)(2,3)

        Asking for block idempotents in any characteristic, by
        passing `p`-cores::

            sage: S0.central_orthogonal_idempotent([1,1])
            Traceback (most recent call last):
            ...
            ValueError: [1, 1] is not a partition of integer 4
            sage: S2.central_orthogonal_idempotent([])
            ()
            sage: S2.central_orthogonal_idempotent([1])
            Traceback (most recent call last):
            ...
            ValueError: the 2-core of [1] is not a 2-core of a partition of 4
            sage: S3.central_orthogonal_idempotent([1])
            () + (1,2)(3,4) + (1,3)(2,4) + (1,4)(2,3)
            sage: S3.central_orthogonal_idempotent([7])
            () + (1,2)(3,4) + (1,3)(2,4) + (1,4)(2,3)

        Asking for classical idempotents::

            sage: S3.central_orthogonal_idempotent([2,2], block=False) is None
            True
            sage: S3.central_orthogonal_idempotent([2,1,1], block=False)
            (3,4) + (2,3) + (2,4) + (1,2) + (1,2)(3,4) + 2*(1,2,3,4)
             + 2*(1,2,4,3) + 2*(1,3,4,2) + (1,3) + (1,3)(2,4)
             + 2*(1,3,2,4) + 2*(1,4,3,2) + (1,4) + 2*(1,4,2,3)
             + (1,4)(2,3)

        .. SEEALSO::

            - :meth:`sage.combinat.partition.Partition.core`
        """
        la = _Partitions(la)
        R = self.base_ring()
        p = R.characteristic()

        if not block or not p:
            if la in self._idempotent_cache:
                return self._idempotent_cache[la]
            if la.size() != self.n:
                raise ValueError("{0} is not a partition of integer {1}".format(la, self.n))
        else:
            mu = la.core(p)
            if mu in self._idempotent_cache:
                return self._idempotent_cache[mu]
            if mu not in self._blocks_dictionary:
                raise ValueError("the {1}-core of {0} is not a {1}-core of a partition of {2}".format(la, p, self.n))

        from sage.libs.gap.libgap import libgap
        from sage.data_structures.blas_dict import iaxpy
        G = self._indices
        character_table = [c.sage() for c in libgap.Irr(libgap.SymmetricGroup(self.n))]
        Pn = Partitions_n(self.n)
        C = Pn.cardinality()
        # We get the indices of the partitions in the reverse lex order
        #   (i.e., reverse of the iteration order of partitions).
        indices = {lam: C - 1 - i for i, lam in enumerate(Pn)}

        if not block or not p:
            la_index = indices[la]
            big_coeff = character_table[la_index][0] / factorial(self.n)
            character_row = character_table[la_index]
            cpi = {g: big_coeff * character_row[indices[g.cycle_type()]]
                   for g in G}
        else:
            # We compute the cycle types of the permutations
            cycles = {}
            for g in G:
                ind = indices[g.cycle_type()]
                if ind in cycles:
                    cycles[ind].append(g)
                else:
                    cycles[ind] = [g]

            denom = factorial(self.n)
            cpi = {}
            for lam in self._blocks_dictionary[mu]:
                lam_index = indices[lam]
                big_coeff = character_table[lam_index][0] / denom
                character_row = character_table[lam_index]
                iaxpy(1,
                      {g: big_coeff * character_row[ind]
                       for ind in cycles for g in cycles[ind]},
                      cpi)

        if not all(R(cpi[g].denominator()) for g in cpi):
            return None

        ret = self.element_class(self, {g: R(cpi[g]) for g in cpi if R(cpi[g])})
        if not block or not p:
            self._idempotent_cache[la] = ret
        else:
            self._idempotent_cache[mu] = ret
        return ret

    @lazy_attribute
    def _blocks_dictionary(self):
        r"""
        Return the partitions of ``self.n``, themselves partitioned
        by their distinct `p`-cores, where `p` is the characteristic
        of ``self.base_ring()``

        If the characteristic is zero, we take the `p`-core operation
        to be the identity map on partitions.

        These lists of partitions, say with common `p`-core `\mu`,
        are components of the central orthogonal idempotent
        corresponding to `\mu`.

        TESTS::

            sage: B2 = SymmetricGroupAlgebra(GF(2), 4)._blocks_dictionary
            sage: [tuple(B2[key]) for key in sorted(B2)]
            [([4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1])]
            sage: B3 = SymmetricGroupAlgebra(GF(3), 4)._blocks_dictionary
            sage: [tuple(B3[key]) for key in sorted(B3)]
            [([4], [2, 2], [1, 1, 1, 1]), ([2, 1, 1],), ([3, 1],)]
            sage: B5 = SymmetricGroupAlgebra(GF(5), 4)._blocks_dictionary
            sage: [tuple(B5[key]) for key in sorted(B5)]
            [([1, 1, 1, 1],), ([2, 1, 1],), ([2, 2],), ([3, 1],), ([4],)]
            sage: B5 == SymmetricGroupAlgebra(QQ, 4)._blocks_dictionary
            True

        .. SEEALSO::

            :meth:`central_orthogonal_idempotent`
        """
        p = self.base_ring().characteristic()
        if not p:
            return {la: [la] for la in Partitions_n(self.n)}

        blocks = {}
        for la in Partitions_n(self.n):
            c = la.core(p)
            if c in blocks:
                blocks[c].append(la)
            else:
                blocks[c] = [la]
        return blocks

    def ladder_idemponent(self, la):
        r"""
        Return the ladder idempontent of ``self``.

        Let `F` be a field of characteristic `p`. The *ladder idempotent*
        of shape `\lambda` is the idempotent of `F[S_n]` defined as follows.
        Let `T` be the :meth:`ladder tableau
        <sage.combinat.partition.Partition.ladder_tableau>` of shape `\lambda`.
        Let `[T]` be the set of standard tableaux whose residue sequence
        is the same as for `T`. Let `\alpha` be the sizes of the ladders
        of `\lambda`. Then the ladder idempontent is constructed as

        .. MATH::

            \widetilde{e}_{\lambda} := \frac{1}{\alpha!}
              \left( \sum_{\sigma \in S_{\alpha}} \sigma \right)
              \left( \overline{\sum_{U \in [T]} E_U} \right),

        where `E_{UU}` is the :meth:`seminormal_basis` element over `\QQ`
        and we project the sum to `F`, `S_{\alpha}` is the Young subgroup
        corresponding to `\alpha`, and `\alpha! = \alpha_1! \cdots \alpha_k!`.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(GF(3), 4)
            sage: for la in Partitions(SGA.n):
            ....:     idem = SGA.ladder_idemponent(la)
            ....:     print(la)
            ....:     print(idem)
            ....:     assert idem^2 == idem
            [4]
            2*[1, 2, 3, 4] + 2*[1, 2, 4, 3] + 2*[2, 1, 3, 4] + 2*[2, 1, 4, 3]
             + 2*[3, 4, 1, 2] + 2*[3, 4, 2, 1] + 2*[4, 3, 1, 2] + 2*[4, 3, 2, 1]
            [3, 1]
            2*[1, 2, 3, 4] + 2*[1, 2, 4, 3] + 2*[2, 1, 3, 4] + 2*[2, 1, 4, 3]
             + [3, 4, 1, 2] + [3, 4, 2, 1] + [4, 3, 1, 2] + [4, 3, 2, 1]
            [2, 2]
            2*[1, 2, 3, 4] + 2*[1, 2, 4, 3] + 2*[2, 1, 3, 4] + 2*[2, 1, 4, 3]
             + 2*[3, 4, 1, 2] + 2*[3, 4, 2, 1] + 2*[4, 3, 1, 2] + 2*[4, 3, 2, 1]
            [2, 1, 1]
            2*[1, 2, 3, 4] + [1, 2, 4, 3] + 2*[1, 3, 2, 4] + [1, 3, 4, 2]
             + [1, 4, 2, 3] + 2*[1, 4, 3, 2] + 2*[2, 1, 3, 4] + [2, 1, 4, 3]
             + 2*[2, 3, 1, 4] + [2, 3, 4, 1] + [2, 4, 1, 3] + 2*[2, 4, 3, 1]
             + 2*[3, 1, 2, 4] + [3, 1, 4, 2] + 2*[3, 2, 1, 4] + [3, 2, 4, 1]
             + [4, 1, 2, 3] + 2*[4, 1, 3, 2] + [4, 2, 1, 3] + 2*[4, 2, 3, 1]
            [1, 1, 1, 1]
            2*[1, 2, 3, 4] + [1, 2, 4, 3] + [2, 1, 3, 4] + 2*[2, 1, 4, 3]
             + 2*[3, 4, 1, 2] + [3, 4, 2, 1] + [4, 3, 1, 2] + 2*[4, 3, 2, 1]

        When `p = 0`, these idempotents will generate all of the simple
        modules (which are the :meth:`Specht modules <specht_module>`
        and also projective modules)::

            sage: SGA = SymmetricGroupAlgebra(QQ, 5)
            sage: for la in Partitions(SGA.n):
            ....:     idem = SGA.ladder_idemponent(la)
            ....:     assert idem^2 == idem
            ....:     print(la, SGA.principal_ideal(idem).dimension())
            [5] 1
            [4, 1] 4
            [3, 2] 5
            [3, 1, 1] 6
            [2, 2, 1] 5
            [2, 1, 1, 1] 4
            [1, 1, 1, 1, 1] 1
            sage: [StandardTableaux(la).cardinality() for la in Partitions(SGA.n)]
            [1, 4, 5, 6, 5, 4, 1]

        REFERENCES:

        - [Ryom2015]_
        """
        R = self.base_ring()
        p = R.characteristic()
        n = self.n
        if not p:
            p = n + 1
        la = _Partitions(la)
        if sum(la) != n:
            raise ValueError(f"{la} is not a partition of {n}")
        Tlad, alpha = la.ladder_tableau(p, ladder_lengths=True)
        if not all(val < p for val in alpha):
            raise ValueError(f"{la} is not {p}-ladder restricted")
        Tclass = Tlad.residue_sequence(p).standard_tableaux()
        Elad = sum(epsilon_ik(T, T) for T in Tclass)
        Elad = self.element_class(self, {sigma: R(c) for sigma, c in Elad._monomial_coefficients.items()})
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        YG = SymmetricGroup(n).young_subgroup(alpha)
        coeff = ~R.prod(factorial(val) for val in alpha)
        G = self.group()
        eprod = self.element_class(self, {G(list(elt.tuple())): coeff
                                          for elt in YG})
        return Elad * eprod

    @cached_method
    def algebra_generators(self):
        r"""
        Return generators of this group algebra (as algebra) as a
        list of permutations.

        The generators used for the group algebra of `S_n` are the
        transposition `(2, 1)` and the `n`-cycle `(1, 2, \ldots, n)`,
        unless `n \leq 1` (in which case no generators are needed).

        EXAMPLES::

            sage: SymmetricGroupAlgebra(ZZ,5).algebra_generators()
            Family ([2, 1, 3, 4, 5], [2, 3, 4, 5, 1])

            sage: SymmetricGroupAlgebra(QQ,0).algebra_generators()
            Family ()

            sage: SymmetricGroupAlgebra(QQ,1).algebra_generators()
            Family ()

        TESTS:

        Check that :issue:`15309` is fixed::

            sage: S3 = SymmetricGroupAlgebra(QQ, 3)
            sage: S3.algebra_generators()
            Family ([2, 1, 3], [2, 3, 1])
            sage: C = CombinatorialFreeModule(ZZ, ZZ)
            sage: M = C.module_morphism(lambda x: S3.zero(), codomain=S3)
            sage: M.register_as_coercion()
        """
        from sage.sets.family import Family
        if self.n <= 1:
            return Family([])
        a = list(range(1, self.n + 1))
        a[0] = 2
        a[1] = 1
        b = list(range(2, self.n + 2))
        b[self.n - 1] = 1
        return Family([self.monomial(self._indices(a)), self.monomial(self._indices(b))])

    def _conjugacy_classes_representatives_underlying_group(self):
        r"""
        Return a complete list of representatives of conjugacy
        classes of the underlying symmetric group.

        .. WARNING::

            This currently is only implemented when ``self`` is built using
            the index set ``Permutations(n)``.

        EXAMPLES::

            sage: SG=SymmetricGroupAlgebra(ZZ,3)
            sage: SG._conjugacy_classes_representatives_underlying_group()
            [[2, 3, 1], [2, 1, 3], [1, 2, 3]]

            sage: SGg = SymmetricGroup(3).algebra(ZZ)
            sage: SGg._conjugacy_classes_representatives_underlying_group() # not tested
            [(1,2,3), (1,2), ()]
        """
        P = self.basis().keys()
        return [P.element_in_conjugacy_classes(nu) for nu in Partitions_n(self.n)]

    def rsw_shuffling_element(self, k):
        r"""
        Return the `k`-th Reiner-Saliola-Welker shuffling element in
        the group algebra ``self``.

        The `k`-th Reiner-Saliola-Welker shuffling element in the
        symmetric group algebra `R S_n` over a ring `R` is defined as the
        sum `\sum_{\sigma \in S_n} \mathrm{noninv}_k(\sigma) \cdot \sigma`,
        where for every permutation `\sigma`, the number
        `\mathrm{noninv}_k(\sigma)` is the number of all
        `k`-noninversions of `\sigma` (that is, the number of all
        `k`-element subsets of `\{ 1, 2, \ldots, n \}` on which
        `\sigma` restricts to a strictly increasing map). See
        :meth:`sage.combinat.permutation.number_of_noninversions` for
        the `\mathrm{noninv}` map.

        This element is more or less the operator `\nu_{k, 1^{n-k}}`
        introduced in [RSW2011]_; more precisely, `\nu_{k, 1^{n-k}}`
        is the left multiplication by this element.

        It is a nontrivial theorem (Theorem 1.1 in [RSW2011]_) that
        the operators `\nu_{k, 1^{n-k}}` (for fixed `n` and varying
        `k`) pairwise commute. It is a conjecture (Conjecture 1.2 in
        [RSW2011]_) that all their eigenvalues are integers (which, in
        light of their commutativity and easily established symmetry,
        yields that they can be simultaneously diagonalized over `\QQ`
        with only integer eigenvalues).

        EXAMPLES:

        The Reiner-Saliola-Welker shuffling elements on `\QQ S_3`::

            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: QS3.rsw_shuffling_element(0)
            [1, 2, 3] + [1, 3, 2] + [2, 1, 3] + [2, 3, 1] + [3, 1, 2] + [3, 2, 1]
            sage: QS3.rsw_shuffling_element(1)
            3*[1, 2, 3] + 3*[1, 3, 2] + 3*[2, 1, 3] + 3*[2, 3, 1] + 3*[3, 1, 2] + 3*[3, 2, 1]
            sage: QS3.rsw_shuffling_element(2)
            3*[1, 2, 3] + 2*[1, 3, 2] + 2*[2, 1, 3] + [2, 3, 1] + [3, 1, 2]
            sage: QS3.rsw_shuffling_element(3)
            [1, 2, 3]
            sage: QS3.rsw_shuffling_element(4)
            0

        Checking the commutativity of Reiner-Saliola-Welker shuffling
        elements (we leave out the ones for which it is trivial)::

            sage: def test_rsw_comm(n):
            ....:     QSn = SymmetricGroupAlgebra(QQ, n)
            ....:     rsws = [QSn.rsw_shuffling_element(k) for k in range(2, n)]
            ....:     return all(ri * rsws[j] == rsws[j] * ri
            ....:                for i, ri in enumerate(rsws) for j in range(i))
            sage: test_rsw_comm(3)
            True
            sage: test_rsw_comm(4)   # long time
            True
            sage: test_rsw_comm(5)   # not tested
            True

        .. NOTE::

            For large ``k`` (relative to ``n``), it might be faster to call
            ``QSn.left_action_product(QSn.semi_rsw_element(k), QSn.antipode(binary_unshuffle_sum(k)))``
            than ``QSn.rsw_shuffling_element(n)``.

        .. SEEALSO::

            :meth:`semi_rsw_element`, :meth:`binary_unshuffle_sum`
        """
        P = self.basis().keys()
        I = Permutations(self.n)
        return self.sum_of_terms([(p, I(p).number_of_noninversions(k)) for p in P],
                                 distinct=True)

    def semi_rsw_element(self, k):
        r"""
        Return the `k`-th semi-RSW element in the group algebra ``self``.

        The `k`-th semi-RSW element in the symmetric group algebra
        `R S_n` over a ring `R` is defined as the sum of all permutations
        `\sigma \in S_n` satisfying
        `\sigma(1) < \sigma(2) < \cdots < \sigma(k)`.

        This element has the property that, if it is denoted by `s_k`,
        then `s_k S(s_k)` is `(n-k)!` times the `k`-th
        Reiner-Saliola-Welker shuffling element of `R S_n` (see
        :meth:`rsw_shuffling_element`). Here, `S` denotes the antipode
        of the group algebra `R S_n`.

        The `k`-th semi-RSW element is the image of the complete
        non-commutative symmetric function `S^{(k, 1^{n-k})}` in the
        ring of non-commutative symmetric functions under the canonical
        projection on the symmetric group algebra (through the descent
        algebra).

        EXAMPLES:

        The semi-RSW elements on `\QQ S_3`::

            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: QS3.semi_rsw_element(0)
            [1, 2, 3] + [1, 3, 2] + [2, 1, 3] + [2, 3, 1] + [3, 1, 2] + [3, 2, 1]
            sage: QS3.semi_rsw_element(1)
            [1, 2, 3] + [1, 3, 2] + [2, 1, 3] + [2, 3, 1] + [3, 1, 2] + [3, 2, 1]
            sage: QS3.semi_rsw_element(2)
            [1, 2, 3] + [1, 3, 2] + [2, 3, 1]
            sage: QS3.semi_rsw_element(3)
            [1, 2, 3]
            sage: QS3.semi_rsw_element(4)
            0

        Let us check the relation with the `k`-th Reiner-Saliola-Welker
        shuffling element stated in the docstring::

            sage: def test_rsw(n):
            ....:     ZSn = SymmetricGroupAlgebra(ZZ, n)
            ....:     for k in range(1, n):
            ....:         a = ZSn.semi_rsw_element(k)
            ....:         b = ZSn.left_action_product(a, ZSn.antipode(a))
            ....:         if factorial(n-k) * ZSn.rsw_shuffling_element(k) != b:
            ....:             return False
            ....:     return True
            sage: test_rsw(3)
            True
            sage: test_rsw(4)
            True
            sage: test_rsw(5)  # long time
            True

        Let us also check the statement about the complete
        non-commutative symmetric function::

            sage: def test_rsw_ncsf(n):
            ....:     ZSn = SymmetricGroupAlgebra(ZZ, n)
            ....:     NSym = NonCommutativeSymmetricFunctions(ZZ)
            ....:     S = NSym.S()
            ....:     for k in range(1, n):
            ....:         a = S(Composition([k] + [1]*(n-k))).to_symmetric_group_algebra()
            ....:         if a != ZSn.semi_rsw_element(k):
            ....:             return False
            ....:     return True
            sage: test_rsw_ncsf(3)
            True
            sage: test_rsw_ncsf(4)
            True
            sage: test_rsw_ncsf(5)  # long time
            True
        """
        n = self.n
        if n < k:
            return self.zero()

        def complement(xs):
            res = list(range(1, n + 1))
            for x in xs:
                res.remove(x)
            return res
        P = Permutations(n)
        I = self._indices
        return self.sum_of_monomials([I(P(complement(q) + list(q)))
                                      for q in itertools.permutations(range(1, n + 1), int(n - k))])

    def binary_unshuffle_sum(self, k):
        r"""
        Return the `k`-th binary unshuffle sum in the group algebra
        ``self``.

        The `k`-th binary unshuffle sum in the symmetric group algebra
        `R S_n` over a ring `R` is defined as the sum of all permutations
        `\sigma \in S_n` satisfying
        `\sigma(1) < \sigma(2) < \cdots < \sigma(k)` and
        `\sigma(k+1) < \sigma(k+2) < \cdots < \sigma(n)`.

        This element has the property that, if it is denoted by `t_k`,
        and if the `k`-th semi-RSW element (see :meth:`semi_rsw_element`)
        is denoted by `s_k`, then `s_k S(t_k)` and `t_k S(s_k)` both
        equal the `k`-th Reiner-Saliola-Welker shuffling element of
        `R S_n` (see :meth:`rsw_shuffling_element`).

        The `k`-th binary unshuffle sum is the image of the complete
        non-commutative symmetric function `S^{(k, n-k)}` in the
        ring of non-commutative symmetric functions under the canonical
        projection on the symmetric group algebra (through the descent
        algebra).

        EXAMPLES:

        The binary unshuffle sums on `\QQ S_3`::

            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: QS3.binary_unshuffle_sum(0)
            [1, 2, 3]
            sage: QS3.binary_unshuffle_sum(1)
            [1, 2, 3] + [2, 1, 3] + [3, 1, 2]
            sage: QS3.binary_unshuffle_sum(2)
            [1, 2, 3] + [1, 3, 2] + [2, 3, 1]
            sage: QS3.binary_unshuffle_sum(3)
            [1, 2, 3]
            sage: QS3.binary_unshuffle_sum(4)
            0

        Let us check the relation with the `k`-th Reiner-Saliola-Welker
        shuffling element stated in the docstring::

            sage: def test_rsw(n):
            ....:     ZSn = SymmetricGroupAlgebra(ZZ, n)
            ....:     for k in range(1, n):
            ....:         a = ZSn.semi_rsw_element(k)
            ....:         b = ZSn.binary_unshuffle_sum(k)
            ....:         c = ZSn.left_action_product(a, ZSn.antipode(b))
            ....:         d = ZSn.left_action_product(b, ZSn.antipode(a))
            ....:         e = ZSn.rsw_shuffling_element(k)
            ....:         if c != e or d != e:
            ....:             return False
            ....:     return True
            sage: test_rsw(3)
            True
            sage: test_rsw(4)  # long time
            True
            sage: test_rsw(5)  # long time
            True

        Let us also check the statement about the complete
        non-commutative symmetric function::

            sage: def test_rsw_ncsf(n):
            ....:     ZSn = SymmetricGroupAlgebra(ZZ, n)
            ....:     NSym = NonCommutativeSymmetricFunctions(ZZ)
            ....:     S = NSym.S()
            ....:     for k in range(1, n):
            ....:         a = S(Composition([k, n-k])).to_symmetric_group_algebra()
            ....:         if a != ZSn.binary_unshuffle_sum(k):
            ....:             return False
            ....:     return True
            sage: test_rsw_ncsf(3)
            True
            sage: test_rsw_ncsf(4)
            True
            sage: test_rsw_ncsf(5)  # long time
            True
        """
        n = self.n
        if n < k:
            return self.zero()

        def complement(xs):
            res = list(range(1, n + 1))
            for x in xs:
                res.remove(x)
            return res
        P = Permutations(n)
        return self.sum_of_monomials([self._indices(P(list(q) + complement(q)))
                                      for q in itertools.combinations(range(1, n + 1), int(k))])

    def specht_module(self, D):
        r"""
        Return the Specht module of ``self`` indexed by the diagram ``D``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 5)
            sage: SM = SGA.specht_module(Partition([3,1,1]))
            sage: SM
            Specht module of [3, 1, 1] over Rational Field
            sage: SM.frobenius_image()
            s[3, 1, 1]

            sage: SM = SGA.specht_module([(1,1),(1,3),(2,2),(3,1),(3,2)])
            sage: SM
            Specht module of [(1, 1), (1, 3), (2, 2), (3, 1), (3, 2)] over Rational Field
            sage: SM.frobenius_image()
            s[2, 2, 1] + s[3, 1, 1] + s[3, 2]
        """
        from sage.combinat.specht_module import SpechtModule
        return SpechtModule(self, D)

    def tabloid_module(self, D):
        r"""
        Return the module of tabloids with the natural action of ``self``.

        .. SEEALSO::

            :class:`~sage.combinat.specht_module.TabloidModule`

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 5)
            sage: TM = SGA.tabloid_module(Partition([3,1,1]))
            sage: TM
            Tabloid module of [3, 1, 1] over Rational Field
            sage: s = SymmetricFunctions(QQ).s()
            sage: s(TM.frobenius_image())
            s[3, 1, 1] + s[3, 2] + 2*s[4, 1] + s[5]
        """
        from sage.combinat.specht_module import TabloidModule
        return TabloidModule(self, D)

    def specht_module_dimension(self, D):
        r"""
        Return the dimension of the Specht module of ``self`` indexed by ``D``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 5)
            sage: SGA.specht_module_dimension(Partition([3,1,1]))
            6
            sage: SGA.specht_module_dimension([(1,1),(1,3),(2,2),(3,1),(3,2)])
            16
        """
        from sage.combinat.specht_module import specht_module_spanning_set, _to_diagram
        D = _to_diagram(D)
        span_set = specht_module_spanning_set(D, self)
        return matrix(self.base_ring(), [v.to_vector() for v in span_set]).rank()

    def simple_module_parameterization(self):
        r"""
        Return a parameterization of the simple modules of ``self``.

        The symmetric group algebra of `S_n` over a field of characteristic `p`
        has its simple modules indexed by all `p`-regular partitions of `n`.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 6)
            sage: SGA.simple_module_parameterization()
            Partitions of the integer 6

            sage: SGA = SymmetricGroupAlgebra(GF(2), 6)
            sage: SGA.simple_module_parameterization()
            2-Regular Partitions of the integer 6
        """
        p = self.base_ring().characteristic()
        if p > 0:
            return Partitions(self.n, regular=p)
        return Partitions_n(self.n)

    def simple_module(self, la):
        r"""
        Return the simple module of ``self`` indexed by the partition ``la``.

        Over a field of characteristic `0`, this simply returns the Specht
        module.

        .. SEEALSO::

            :class:`sage.combinat.specht_module.SimpleModule`

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(GF(3), 5)
            sage: D = SGA.simple_module(Partition([3,1,1]))
            sage: D
            Simple module of [3, 1, 1] over Finite Field of size 3
            sage: D.brauer_character()
            (6, 0, -2, 0, 1)
        """
        from sage.combinat.specht_module import SpechtModule
        return SpechtModule(self, la).simple_module()

    def simple_module_dimension(self, la):
        r"""
        Return the dimension of the simple module of ``self`` indexed by the
        partition ``la``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(GF(5), 6)
            sage: SGA.simple_module_dimension(Partition([4,1,1]))
            10

        TESTS::

            sage: SGA = SymmetricGroupAlgebra(GF(5), 6)
            sage: SGA.simple_module_dimension(Partition([3,1,1]))
            Traceback (most recent call last):
            ...
            ValueError: [3, 1, 1] is not a partition of 6
        """
        if sum(la) != self.n:
            raise ValueError(f"{la} is not a partition of {self.n}")
        from sage.combinat.specht_module import simple_module_rank
        return simple_module_rank(la, self.base_ring())

    def garsia_procesi_module(self, la):
        r"""
        Return the :class:`Garsia-Procesi module
        <sage.combinat.symmetric_group_representations.GarsiaProcesiModule>`
        of ``self`` indexed by ``la``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(GF(2), 6)
            sage: SGA.garsia_procesi_module(Partition([2,2,1,1]))
            Garsia-Procesi module of shape [2, 2, 1, 1] over Finite Field of size 2
        """
        from sage.combinat.symmetric_group_representations import GarsiaProcesiModule
        return GarsiaProcesiModule(self, la)

    def jucys_murphy(self, k):
        r"""
        Return the Jucys-Murphy element `J_k` (also known as a
        Young-Jucys-Murphy element) for the symmetric group
        algebra ``self``.

        The Jucys-Murphy element `J_k` in the symmetric group algebra
        `R S_n` is defined for every `k \in \{ 1, 2, \ldots, n \}` by

        .. MATH::

            J_k = (1, k) + (2, k) + \cdots + (k-1, k) \in R S_n,

        where the addends are transpositions in `S_n` (regarded as
        elements of `R S_n`). We note that there is not a dependence on `n`,
        so it is often suppressed in the notation.

        EXAMPLES::

            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: QS3.jucys_murphy(1)
            0
            sage: QS3.jucys_murphy(2)
            [2, 1, 3]
            sage: QS3.jucys_murphy(3)
            [1, 3, 2] + [3, 2, 1]

            sage: QS4 = SymmetricGroupAlgebra(QQ, 4)
            sage: j3 = QS4.jucys_murphy(3); j3
            [1, 3, 2, 4] + [3, 2, 1, 4]
            sage: j4 = QS4.jucys_murphy(4); j4
            [1, 2, 4, 3] + [1, 4, 3, 2] + [4, 2, 3, 1]
            sage: j3*j4 == j4*j3
            True

            sage: QS5 = SymmetricGroupAlgebra(QQ, 5)
            sage: QS5.jucys_murphy(4)
            [1, 2, 4, 3, 5] + [1, 4, 3, 2, 5] + [4, 2, 3, 1, 5]

        TESTS::

            sage: QS3.jucys_murphy(4)
            Traceback (most recent call last):
            ...
            ValueError: k (= 4) must be between 1 and n (= 3) (inclusive)
        """
        if k < 1 or k > self.n:
            raise ValueError("k (= {k}) must be between 1 and n (= {n}) (inclusive)".format(k=k, n=self.n))

        res = self.zero()

        for i in range(1, k):
            p = list(range(1, self.n + 1))
            p[i - 1] = k
            p[k - 1] = i
            res += self.monomial(self._indices(p))
        return res

    def seminormal_basis(self, mult='l2r'):
        r"""
        Return a list of the seminormal basis elements of ``self``.

        The seminormal basis of a symmetric group algebra is defined as
        follows:

        Let `n` be a nonnegative integer. Let `R` be a `\QQ`-algebra.
        In the following, we will use the "left action" convention for
        multiplying permutations. This means that for all permutations
        `p` and `q` in `S_n`, the product `pq` is defined in such a way
        that `(pq)(i) = p(q(i))` for each `i \in \{ 1, 2, \ldots, n \}`
        (this is the same convention as in :meth:`left_action_product`,
        but not the default semantics of the `*` operator on
        permutations in Sage). Thus, for instance, `s_2 s_1` is the
        permutation obtained by first transposing `1` with `2` and
        then transposing `2` with `3` (where `s_i = (i, i+1)`).

        For every partition `\lambda` of `n`, let

        .. MATH::

            \kappa_{\lambda} = \frac{n!}{f^{\lambda}}

        where `f^{\lambda}` is the number of standard Young tableaux
        of shape `\lambda`. Note that `\kappa_{\lambda}` is an integer,
        namely the product of all hook lengths of `\lambda` (by the
        hook length formula). In Sage, this integer can be computed by
        using :func:`sage.combinat.symmetric_group_algebra.kappa()`.

        Let `T` be a standard tableau of size `n`.

        Let `a(T)` denote the formal sum (in `R S_n`) of all
        permutations in `S_n` which stabilize the rows of `T` (as
        sets), i. e., which map each entry `i` of `T` to an entry in
        the same row as `i`. (See
        :func:`sage.combinat.symmetric_group_algebra.a()` for
        an implementation of this.)

        Let `b(T)` denote the signed formal sum (in `R S_n`) of all
        permutations in `S_n` which stabilize the columns of `T` (as
        sets). Here, "signed" means that each permutation is
        multiplied with its sign. (This is implemented in
        :func:`sage.combinat.symmetric_group_algebra.b()`.)

        Define an element `e(T)` of `R S_n` to be `a(T) b(T)`. (This
        is implemented in :func:`sage.combinat.symmetric_group_algebra.e()`
        for `R = \QQ`.)

        Let `\mathrm{sh}(T)` denote the shape of `T`.
        (See :meth:`~sage.combinat.tableau.Tableau.shape`.)

        Let `\overline{T}` denote the standard tableau of size `n-1`
        obtained by removing the letter `n` (along with its cell) from
        `T` (if `n \geq 1`).

        Now, we define an element `\epsilon(T)` of `R S_n`. We define
        it by induction on the size `n` of `T`, so we set
        `\epsilon(\emptyset) = 1` and only need to define `\epsilon(T)`
        for `n \geq 1`, assuming that `\epsilon(\overline{T})` is
        already defined. We do this by setting

        .. MATH::

            \epsilon(T) = \frac{1}{\kappa_{\mathrm{sh}(T)}}
                          \epsilon(\overline{T})
                          e(T) \epsilon(\overline{T}).

        This element `\epsilon(T)` is implemented as
        :func:`sage.combinat.symmetric_group_algebra.epsilon` for
        `R = \QQ`, but it is also a particular case of the elements
        `\epsilon(T, S)` defined below.

        Now let `S` be a further tableau of the same shape as `T`
        (possibly equal to `T`). Let `\pi_{T, S}` denote the
        permutation in `S_n` such that applying this permutation to
        the entries of `T` yields the tableau `S`. Define an element
        `\epsilon(T, S)` of `R S_n` by

        .. MATH::

            \epsilon(T, S) = \frac{1}{\kappa_{\mathrm{sh}(T)}}
                             \epsilon(\overline S) \pi_{T, S}
                             e(T) \epsilon(\overline T)
                           = \frac{1}{\kappa_{\mathrm{sh}(T)}}
                             \epsilon(\overline S) a(S) \pi_{T, S}
                             b(T) \epsilon(\overline T).

        This element `\epsilon(T, S)` is called *Young's seminormal
        unit corresponding to the bitableau `(T, S)`*, and is the
        return value of :meth:`epsilon_ik` applied to ``T`` and
        ``S``. Note that `\epsilon(T, T) = \epsilon(T)`.

        If we let `\lambda` run through all partitions of `n`, and
        `(T, S)` run through all pairs of tableaux of shape
        `\lambda`, then the elements `\epsilon(T, S)` form a basis
        of `R S_n`. This basis is called *Young's seminormal basis*
        and has the properties that

        .. MATH::

            \epsilon(T, S) \epsilon(U, V) = \delta_{T, V} \epsilon(U, S)

        (where `\delta` stands for the Kronecker delta).

        .. WARNING::

            Because of our convention, we are multiplying our elements in
            reverse of those given in some papers, for example [Ram1997]_.
            Using the other convention of multiplying permutations, we would
            instead have
            `\epsilon(U, V) \epsilon(T, S) = \delta_{T, V} \epsilon(U, S)`.

        In other words, Young's seminormal basis consists of the matrix
        units in a (particular) Artin-Wedderburn decomposition of `R S_n`
        into a direct product of matrix algebras over `\QQ`.

        The output of :meth:`seminormal_basis` is a list of all
        elements of the seminormal basis of ``self``.

        INPUT:

        - ``mult`` -- string (default: ``'l2r'``); if set to ``'r2l'``,
          this causes the method to return the list of the
          antipodes (:meth:`antipode`) of all `\epsilon(T, S)`
          instead of the `\epsilon(T, S)` themselves.

        EXAMPLES::

            sage: QS3 = SymmetricGroupAlgebra(QQ,3)
            sage: QS3.seminormal_basis()
            [1/6*[1, 2, 3] + 1/6*[1, 3, 2] + 1/6*[2, 1, 3] + 1/6*[2, 3, 1] + 1/6*[3, 1, 2] + 1/6*[3, 2, 1],
            1/3*[1, 2, 3] + 1/6*[1, 3, 2] - 1/3*[2, 1, 3] - 1/6*[2, 3, 1] - 1/6*[3, 1, 2] + 1/6*[3, 2, 1],
            1/3*[1, 3, 2] + 1/3*[2, 3, 1] - 1/3*[3, 1, 2] - 1/3*[3, 2, 1],
            1/4*[1, 3, 2] - 1/4*[2, 3, 1] + 1/4*[3, 1, 2] - 1/4*[3, 2, 1],
            1/3*[1, 2, 3] - 1/6*[1, 3, 2] + 1/3*[2, 1, 3] - 1/6*[2, 3, 1] - 1/6*[3, 1, 2] - 1/6*[3, 2, 1],
            1/6*[1, 2, 3] - 1/6*[1, 3, 2] - 1/6*[2, 1, 3] + 1/6*[2, 3, 1] + 1/6*[3, 1, 2] - 1/6*[3, 2, 1]]

        TESTS::

            sage: QS3g = SymmetricGroup(3).algebra(QQ)
            sage: QS3g.seminormal_basis()
            [1/6*() + 1/6*(2,3) + 1/6*(1,2) + 1/6*(1,2,3) + 1/6*(1,3,2) + 1/6*(1,3),
             1/3*() + 1/6*(2,3) - 1/3*(1,2) - 1/6*(1,2,3) - 1/6*(1,3,2) + 1/6*(1,3),
             1/3*(2,3) + 1/3*(1,2,3) - 1/3*(1,3,2) - 1/3*(1,3),
             1/4*(2,3) - 1/4*(1,2,3) + 1/4*(1,3,2) - 1/4*(1,3),
             1/3*() - 1/6*(2,3) + 1/3*(1,2) - 1/6*(1,2,3) - 1/6*(1,3,2) - 1/6*(1,3),
             1/6*() - 1/6*(2,3) - 1/6*(1,2) + 1/6*(1,2,3) + 1/6*(1,3,2) - 1/6*(1,3)]
        """
        basis = []
        for part in Partitions_n(self.n):
            stp = StandardTableaux_shape(part)
            basis.extend(self.epsilon_ik(t1, t2, mult=mult)
                         for t1 in stp
                         for t2 in stp)
        return basis

    def dft(self, form=None, mult='l2r'):
        r"""
        Return the discrete Fourier transform for ``self``.

        See [Mur1983]_ for the construction of central primitive orthogonal idempotents.
        For each idempotent `e_i` we have a homomorphic projection `v \mapsto v e_i`.
        Choose a basis for each submodule spanned by `\{\sigma e_i | \sigma \in S_n\}`.
        The change-of-basis from the standard basis `\{\sigma\}_\sigma` is returned.

        INPUT:

        - ``mult`` -- string (default: `l2r`); if set to `r2l`,
          this causes the method to use the antipodes
          (:meth:`antipode`) of the seminormal basis instead of
          the seminormal basis.

        EXAMPLES::

            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: QS3.dft()
            [   1    1    1    1    1    1]
            [   1  1/2   -1 -1/2 -1/2  1/2]
            [   0  3/4    0  3/4 -3/4 -3/4]
            [   0    1    0   -1    1   -1]
            [   1 -1/2    1 -1/2 -1/2 -1/2]
            [   1   -1   -1    1    1   -1]

        Over fields of characteristic `p > 0` such that `p \mid n!`, we use the
        modular Fourier transform (:issue:`37751`)::

            sage: GF2S3 = SymmetricGroupAlgebra(GF(2), 3)
            sage: GF2S3.dft()
            [1 0 0 0 1 0]
            [0 1 0 0 0 1]
            [0 0 1 0 0 1]
            [0 0 0 1 1 0]
            [1 0 0 1 1 0]
            [0 1 1 0 0 1]
        """
        if form is None:
            form = "modular" if self.base_ring().characteristic().divides(self.group().cardinality()) else "seminormal"
        if form == "seminormal":
            if self.base_ring().characteristic().divides(self.group().cardinality()):
                raise ValueError("seminormal does not work when p | n!")
            return self._dft_seminormal(mult=mult)
        if form == "modular":
            return self._dft_modular()
        raise ValueError("invalid form (= %s)" % form)

    def _dft_seminormal(self, mult='l2r'):
        """
        Return the seminormal form of the discrete Fourier transform for ``self``.

        INPUT:

        - ``mult`` -- string (default: `l2r`); if set to `r2l`,
          this causes the method to use the antipodes
          (:meth:`antipode`) of the seminormal basis instead of
          the seminormal basis.

        EXAMPLES::

            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: QS3._dft_seminormal()
            [   1    1    1    1    1    1]
            [   1  1/2   -1 -1/2 -1/2  1/2]
            [   0  3/4    0  3/4 -3/4 -3/4]
            [   0    1    0   -1    1   -1]
            [   1 -1/2    1 -1/2 -1/2 -1/2]
            [   1   -1   -1    1    1   -1]

        .. SEEALSO::

            :meth:`seminormal_basis`
        """
        snb = self.seminormal_basis(mult=mult)
        return matrix([vector(b) for b in snb]).inverse().transpose()

    def _dft_modular(self):
        r"""
        Return the discrete Fourier transform when the characteristic divides the order of the group.

        EXAMPLES::

            sage: GF3S3 = SymmetricGroupAlgebra(GF(3), 3)
            sage: GF3S3._dft_modular()
            [1 0 0 0 0 0]
            [0 1 0 0 0 0]
            [0 0 1 0 0 0]
            [0 0 0 1 0 0]
            [0 0 0 0 1 0]
            [0 0 0 0 0 1]
        """
        idempotents = self.central_orthogonal_idempotents()
        # project v onto each block U_i = F_p[S_n]*e_i via \pi_i: v |--> v*e_i
        B = self.basis()
        blocks = [self.submodule([b * idem for b in B]) for idem in idempotents]
        # compute the list of basis vectors lifted to the SGA from each block
        block_decomposition_basis = [u.lift() for block in blocks for u in block.basis()]
        # construct the matrix to the standard basis in the order given by the group
        G = self.group()
        mat = [[b[g] for b in block_decomposition_basis] for g in G]
        return matrix(self.base_ring(), mat)

    def epsilon_ik(self, itab, ktab, star=0, mult='l2r'):
        r"""
        Return the seminormal basis element of ``self`` corresponding to the
        pair of tableaux ``itab`` and ``ktab`` (or restrictions of these
        tableaux, if the optional variable ``star`` is set).

        INPUT:

        - ``itab``, ``ktab`` -- two standard tableaux of size `n`

        - ``star`` -- integer (default: `0`)

        - ``mult`` -- string (default: `l2r`); if set to `r2l`,
          this causes the method to return the antipode
          (:meth:`antipode`) of `\epsilon(I, K)` instead of
          `\epsilon(I, K)` itself.

        OUTPUT:

        The element `\epsilon(I, K)`, where `I` and `K` are the tableaux
        obtained by removing all entries higher than `n - \mathrm{star}`
        from ``itab`` and ``ktab``, respectively. Here, we are using the
        notations from :meth:`seminormal_basis`.

        EXAMPLES::

            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: a = QS3.epsilon_ik([[1,2,3]], [[1,2,3]]); a
            1/6*[1, 2, 3] + 1/6*[1, 3, 2] + 1/6*[2, 1, 3] + 1/6*[2, 3, 1] + 1/6*[3, 1, 2] + 1/6*[3, 2, 1]
            sage: QS3.dft()*vector(a)
            (1, 0, 0, 0, 0, 0)
            sage: a = QS3.epsilon_ik([[1,2],[3]], [[1,2],[3]]); a
            1/3*[1, 2, 3] - 1/6*[1, 3, 2] + 1/3*[2, 1, 3] - 1/6*[2, 3, 1] - 1/6*[3, 1, 2] - 1/6*[3, 2, 1]
            sage: QS3.dft()*vector(a)
            (0, 0, 0, 0, 1, 0)

        Let us take some properties of the seminormal basis listed in
        the docstring of :meth:`seminormal_basis`, and verify them on
        the situation of `S_3`.

        First, check the formula

        .. MATH::

            \epsilon(T) = \frac{1}{\kappa_{\mathrm{sh}(T)}}
                          \epsilon(\overline{T})
                          e(T) \epsilon(\overline{T}).

        In fact::

            sage: from sage.combinat.symmetric_group_algebra import e
            sage: def test_sn1(n):
            ....:     QSn = SymmetricGroupAlgebra(QQ, n)
            ....:     QSn1 = SymmetricGroupAlgebra(QQ, n - 1)
            ....:     for T in StandardTableaux(n):
            ....:         TT = T.restrict(n-1)
            ....:         eTT = QSn1.epsilon_ik(TT, TT)
            ....:         eT = QSn.epsilon_ik(T, T)
            ....:         kT = prod(T.shape().hooks())
            ....:         if kT * eT != eTT * e(T) * eTT:
            ....:             return False
            ....:     return True
            sage: test_sn1(3)
            True
            sage: test_sn1(4)   # long time
            True

        Next, we check the identity

        .. MATH::

            \epsilon(T, S) = \frac{1}{\kappa_{\mathrm{sh}(T)}}
                             \epsilon(\overline S) \pi_{T, S}
                             e(T) \epsilon(\overline T)

        which we used to define `\epsilon(T, S)`. In fact::

            sage: from sage.combinat.symmetric_group_algebra import e
            sage: def test_sn2(n):
            ....:     QSn = SymmetricGroupAlgebra(QQ, n)
            ....:     mul = QSn.left_action_product
            ....:     QSn1 = SymmetricGroupAlgebra(QQ, n - 1)
            ....:     for lam in Partitions(n):
            ....:         k = prod(lam.hooks())
            ....:         for T in StandardTableaux(lam):
            ....:             for S in StandardTableaux(lam):
            ....:                 TT = T.restrict(n-1)
            ....:                 SS = S.restrict(n-1)
            ....:                 eTT = QSn1.epsilon_ik(TT, TT)
            ....:                 eSS = QSn1.epsilon_ik(SS, SS)
            ....:                 eTS = QSn.epsilon_ik(T, S)
            ....:                 piTS = [0] * n
            ....:                 for (i, j) in T.cells():
            ....:                     piTS[T[i][j] - 1] = S[i][j]
            ....:                 piTS = QSn(Permutation(piTS))
            ....:                 if k * eTS != mul(mul(eSS, piTS), mul(e(T), eTT)):
            ....:                     return False
            ....:     return True
            sage: test_sn2(3)
            True
            sage: test_sn2(4)   # long time
            True

        Let us finally check the identity

        .. MATH::

            \epsilon(T, S) \epsilon(U, V) = \delta_{T, V} \epsilon(U, S)

        In fact::

            sage: def test_sn3(lam):
            ....:     n = lam.size()
            ....:     QSn = SymmetricGroupAlgebra(QQ, n)
            ....:     mul = QSn.left_action_product
            ....:     for T in StandardTableaux(lam):
            ....:         for S in StandardTableaux(lam):
            ....:             for U in StandardTableaux(lam):
            ....:                 for V in StandardTableaux(lam):
            ....:                     lhs = mul(QSn.epsilon_ik(T, S), QSn.epsilon_ik(U, V))
            ....:                     if T == V:
            ....:                         rhs = QSn.epsilon_ik(U, S)
            ....:                     else:
            ....:                         rhs = QSn.zero()
            ....:                     if rhs != lhs:
            ....:                         return False
            ....:     return True
            sage: all( test_sn3(lam) for lam in Partitions(3) )
            True
            sage: all( test_sn3(lam) for lam in Partitions(4) )   # long time
            True
        """
        it = Tableau(itab)
        kt = Tableau(ktab)

        stn = StandardTableaux_size(self.n)

        if it not in stn:
            raise TypeError("it must be a standard tableau of size %s" % self.n)

        if kt not in stn:
            raise TypeError("kt must be a standard tableau of size %s" % self.n)

        if it.shape() != kt.shape():
            raise ValueError("it and kt must be of the same shape")

        BR = self.base_ring()
        I = self._indices
        z_elts = {}
        epik = epsilon_ik(it, kt, star=star)
        for m, c in epik._monomial_coefficients.items():
            z_elts[I(m)] = BR(c)
        z = self._from_dict(z_elts)

        if mult == 'l2r':
            return z
        else:
            return z.map_support(lambda x: x.inverse())

    def murphy_basis(self):
        r"""
        Return the :class:`Murphy basis
        <sage.combinat.symmetric_group_algebra.MurphyBasis>` of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 3)
            sage: M = SGA.murphy_basis()
            sage: M(SGA.an_element())
            -C([1, 1, 1], [[1], [2], [3]], [[1], [2], [3]])
             + C([2, 1], [[1, 2], [3]], [[1, 2], [3]])
             + C([2, 1], [[1, 2], [3]], [[1, 3], [2]])
             + 2*C([2, 1], [[1, 3], [2]], [[1, 2], [3]])
             + 4*C([2, 1], [[1, 3], [2]], [[1, 3], [2]])
             - 3*C([3], [[1, 2, 3]], [[1, 2, 3]])
        """
        return MurphyBasis(self)

    @cached_method(key=lambda s, X, Y: (StandardTableaux()(X), StandardTableaux()(Y)))
    def murphy_basis_element(self, S, T):
        r"""
        Return the Murphy basis element indexed by ``S`` and ``T``.

        .. SEEALSO::

            :class:`~sage.combinat.symmetric_group_algebra.MurphyBasis`

        EXAMPLES::

            sage: import itertools
            sage: SGA = SymmetricGroupAlgebra(QQ, 3)
            sage: for S, T in itertools.product(StandardTableaux([2,1]), repeat=2):
            ....:     print(S, T, SGA.murphy_basis_element(S, T))
            [[1, 3], [2]] [[1, 3], [2]] [1, 2, 3] + [2, 1, 3]
            [[1, 3], [2]] [[1, 2], [3]] [1, 3, 2] + [3, 1, 2]
            [[1, 2], [3]] [[1, 3], [2]] [1, 3, 2] + [2, 3, 1]
            [[1, 2], [3]] [[1, 2], [3]] [1, 2, 3] + [3, 2, 1]

        TESTS::

            sage: SGA = SymmetricGroupAlgebra(QQ, 3)
            sage: SGA.murphy_basis_element([[1,2,3,4]], [[1,2],[3,4]])
            Traceback (most recent call last):
            ...
            ValueError: [[1, 2, 3, 4]] is not an element of Standard tableaux of size 3
            sage: SGA.murphy_basis_element([[1,2,3]], [[1,2],[3]])
            Traceback (most recent call last):
            ...
            ValueError: S and T must have the same shape
        """
        std_tab = StandardTableaux(self.n)
        S = std_tab(S)
        T = std_tab(T)
        S = S.conjugate()
        T = T.conjugate()

        la = S.shape()
        if la != T.shape():
            raise ValueError("S and T must have the same shape")

        G = self.group()
        ds = G(list(sum((row for row in S), ())))
        dt = G(list(sum((row for row in T), ())))
        return self.monomial(~ds) * self._row_stabilizer(la) * self.monomial(dt)

    @cached_method
    def _row_stabilizer(self, la):
        """
        Return the row stabilizer element of a canonical standard tableau
        of shape ``la``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 3)
            sage: for la in Partitions(3):
            ....:     print(la, SGA._row_stabilizer(la))
            [3] [1, 2, 3] + [1, 3, 2] + [2, 1, 3] + [2, 3, 1] + [3, 1, 2] + [3, 2, 1]
            [2, 1] [1, 2, 3] + [2, 1, 3]
            [1, 1, 1] [1, 2, 3]
        """
        G = self.group()
        return self.sum_of_monomials(G(list(w.tuple())) for w in la.young_subgroup())

    @cached_method
    def _column_antistabilizer(self, la):
        """
        Return the column antistabilizer element of a canonical standard tableau
        of shape ``la``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 3)
            sage: for la in Partitions(3):
            ....:     print(la, SGA._column_antistabilizer(la))
            [3] [1, 2, 3]
            [2, 1] [1, 2, 3] - [3, 2, 1]
            [1, 1, 1] [1, 2, 3] - [1, 3, 2] - [2, 1, 3] + [2, 3, 1] + [3, 1, 2] - [3, 2, 1]
        """
        T = []
        total = 1 # make it 1-based
        for r in la:
            T.append(list(range(total, total+r)))
            total += r
        T = Tableau(T)
        G = self.group()
        R = self.base_ring()
        return self._from_dict({G(list(w.tuple())): R(w.sign()) for w in T.column_stabilizer()},
                               remove_zeros=False)

    @cached_method
    def _young_symmetrizer(self, la):
        """
        Return the Young symmetrizer of shape ``la`` of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 3)
            sage: for la in Partitions(3):
            ....:     print(la, SGA._young_symmetrizer(la))
            [3] [1, 2, 3] + [1, 3, 2] + [2, 1, 3] + [2, 3, 1] + [3, 1, 2] + [3, 2, 1]
            [2, 1] [1, 2, 3] + [2, 1, 3] - [3, 1, 2] - [3, 2, 1]
            [1, 1, 1] [1, 2, 3] - [1, 3, 2] - [2, 1, 3] + [2, 3, 1] + [3, 1, 2] - [3, 2, 1]
        """
        return self._column_antistabilizer(la) * self._row_stabilizer(la)

    def young_symmetrizer(self, la):
        """
        Return the Young symmetrizer of shape ``la`` of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, SymmetricGroup(3))
            sage: SGA.young_symmetrizer([2,1])
            () + (1,2) - (1,3,2) - (1,3)
            sage: SGA = SymmetricGroupAlgebra(QQ, 4)
            sage: SGA.young_symmetrizer([2,1,1])
            [1, 2, 3, 4] - [1, 2, 4, 3] + [2, 1, 3, 4] - [2, 1, 4, 3]
             - [3, 1, 2, 4] + [3, 1, 4, 2] - [3, 2, 1, 4] + [3, 2, 4, 1]
             + [4, 1, 2, 3] - [4, 1, 3, 2] + [4, 2, 1, 3] - [4, 2, 3, 1]
            sage: SGA.young_symmetrizer([5,1,1])
            Traceback (most recent call last):
            ...
            ValueError: the partition [5, 1, 1] is not of size 4
        """
        la = _Partitions(la)
        if la.size() != self.n:
            raise ValueError("the partition {} is not of size {}".format(la, self.n))
        return self._young_symmetrizer(la)

    def kazhdan_lusztig_cellular_basis(self):
        r"""
        Return the Kazhdan-Lusztig basis (at `q = 1`) of ``self``
        as a cellular basis.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 3)
            sage: KL = SGA.kazhdan_lusztig_cellular_basis()
            sage: KL(SGA.an_element())
            C([2, 1], [[1, 2], [3]], [[1, 2], [3]])
             + C([2, 1], [[1, 3], [2]], [[1, 2], [3]])
             + 2*C([2, 1], [[1, 3], [2]], [[1, 3], [2]])
             - 3*C([3], [[1, 2, 3]], [[1, 2, 3]])
        """
        return KLCellularBasis(self)

    @cached_method
    def kazhdan_lusztig_basis_element(self, w):
        r"""
        Return the Kazhdan-Lusztig `C'_w` basis element at `q = 1`.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 3)
            sage: for w in SGA.group():
            ....:     print(w, SGA.kazhdan_lusztig_basis_element(w))
            [1, 2, 3] [1, 2, 3]
            [1, 3, 2] [1, 2, 3] + [1, 3, 2]
            [2, 1, 3] [1, 2, 3] + [2, 1, 3]
            [2, 3, 1] [1, 2, 3] + [1, 3, 2] + [2, 1, 3] + [2, 3, 1]
            [3, 1, 2] [1, 2, 3] + [1, 3, 2] + [2, 1, 3] + [3, 1, 2]
            [3, 2, 1] [1, 2, 3] + [1, 3, 2] + [2, 1, 3] + [2, 3, 1] + [3, 1, 2] + [3, 2, 1]
        """
        from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet
        G = self.basis().keys()
        R = self.base_ring()
        one = R.one()
        # check if the KL polynomials can be computed using ``coxeter3``
        try:
            from sage.libs.coxeter3.coxeter_group import CoxeterGroup as Coxeter3Group
        except ImportError:
            # Falback to using the KL polynomial
            from sage.combinat.kazhdan_lusztig import KazhdanLusztigPolynomial
            from sage.groups.perm_gps.permgroup_named import SymmetricGroup
            q = PolynomialRing(R, 'q').gen()
            self._KLG = SymmetricGroup(self.n)
            self._cellular_KL = KazhdanLusztigPolynomial(self._KLG, q)
            polyfunc = self._cellular_KL.P
        else:
            self._cellular_KL = Coxeter3Group(['A', self.n+1])
            self._KLG = self._cellular_KL
            polyfunc = self._cellular_KL.kazhdan_lusztig_polynomial

        if w.parent() is not self._KLG:
            w = self._KLG.from_reduced_word(w.reduced_word())
        bruhat = RecursivelyEnumeratedSet([w], lambda u: u.bruhat_lower_covers(), structure='graded')
        return self.element_class(self, {G.from_reduced_word(v.reduced_word()): R(c(q=one))
                                         for v in bruhat if (c := polyfunc(v, w))})


epsilon_ik_cache = {}


def epsilon_ik(itab, ktab, star=0):
    r"""
    Return the seminormal basis element of the symmetric group
    algebra `\QQ S_n` corresponding to the pair of tableaux
    ``itab`` and ``ktab`` (or restrictions of these tableaux,
    if the optional variable ``star`` is set).

    INPUT:

    - ``itab``, ``ktab`` -- two standard tableaux of same size

    - ``star`` -- integer (default: `0`)

    OUTPUT:

    The element `\epsilon(I, K) \in \QQ S_n`, where `I` and `K`
    are the tableaux obtained by removing all entries higher
    than `n - \mathrm{star}` from ``itab`` and ``ktab``,
    respectively (where `n` is the size of ``itab`` and
    ``ktab``). Here, we are using the notations from
    :meth:`~sage.combinat.symmetric_group_algebra.SymmetricGroupAlgebra_n.seminormal_basis`.

    EXAMPLES::

        sage: from sage.combinat.symmetric_group_algebra import epsilon_ik
        sage: epsilon_ik([[1,2],[3]], [[1,3],[2]])
        1/4*[1, 3, 2] - 1/4*[2, 3, 1] + 1/4*[3, 1, 2] - 1/4*[3, 2, 1]
        sage: epsilon_ik([[1,2],[3]], [[1,3],[2]], star=1)
        Traceback (most recent call last):
        ...
        ValueError: the two tableaux must be of the same shape
    """
    it = Tableau(itab)
    kt = Tableau(ktab)
    if star:
        it = it.restrict(it.size() - star)
        kt = kt.restrict(kt.size() - star)

    if it.shape() != kt.shape():
        raise ValueError("the two tableaux must be of the same shape")

    if kt == it:
        res = epsilon(itab)
    elif (it, kt) in epsilon_ik_cache:
        res = epsilon_ik_cache[(it, kt)]
    else:
        eik = e_ik(it, kt, star)
        QSn = eik.parent()
        mul = QSn.right_action_product
        epsilon_ik_cache[(it, kt)] = mul(mul(epsilon(it, star + 1), eik),
                                         epsilon(kt, star + 1)) * (1 / kappa(it.shape()))
        res = epsilon_ik_cache[(it, kt)]

    return res


epsilon_cache = {}


def epsilon(tab, star=0):
    r"""
    The `(T, T)`-th element of the seminormal basis of the group
    algebra `\QQ[S_n]`, where `T` is the tableau ``tab`` (with its
    ``star`` highest entries removed if the optional variable
    ``star`` is set).

    See the docstring of
    :meth:`~sage.combinat.symmetric_group_algebra.SymmetricGroupAlgebra_n.seminormal_basis`
    for the notation used herein.

    EXAMPLES::

        sage: from sage.combinat.symmetric_group_algebra import epsilon
        sage: epsilon([[1,2]])
        1/2*[1, 2] + 1/2*[2, 1]
        sage: epsilon([[1],[2]])
        1/2*[1, 2] - 1/2*[2, 1]
    """
    t = Tableau(tab)

    if star:
        t = t.restrict(t.size() - star)

    if t in epsilon_cache:
        res = epsilon_cache[t]
    else:
        if t.size() == 2:
            epsilon_cache[t] = e(t) * (1 / kappa(t.shape()))
            res = epsilon_cache[t]
        elif t == Tableau([[1]]):
            epsilon_cache[t] = e(t)
            res = epsilon_cache[t]
        else:
            et = e(t)
            QSn = et.parent()
            mul = QSn.right_action_product
            epsilon_cache[t] = mul(mul(epsilon(t, 1), e(t)),
                                   epsilon(t, 1)) * (1 / kappa(t.shape()))
            res = epsilon_cache[t]

    return res


def pi_ik(itab, ktab):
    r"""
    Return the permutation `p` which sends every entry of the
    tableau ``itab`` to the respective entry of the tableau
    ``ktab``, as an element of the corresponding symmetric group
    algebra.

    This assumes that ``itab`` and ``ktab`` are tableaux (possibly
    given just as lists of lists) of the same shape. Both
    tableaux are allowed to be skew.

    EXAMPLES::

        sage: from sage.combinat.symmetric_group_algebra import pi_ik
        sage: pi_ik([[1,3],[2]], [[1,2],[3]])
        [1, 3, 2]

    The same with skew tableaux::

        sage: from sage.combinat.symmetric_group_algebra import pi_ik
        sage: pi_ik([[None,1,3],[2]], [[None,1,2],[3]])
        [1, 3, 2]
    """
    it = SkewTableau(itab)
    kt = SkewTableau(ktab)
    n = kt.size()

    p = [None] * n
    for i in range(len(kt)):
        for j in range(len(kt[i])):
            if it[i][j] is not None:
                p[it[i][j] - 1] = kt[i][j]

    QSn = SymmetricGroupAlgebra(QQ, n)
    p = Permutation(p)
    return QSn(p)


def kappa(alpha):
    r"""
    Return `\kappa_\alpha`, which is `n!` divided by the number
    of standard tableaux of shape `\alpha` (where `\alpha` is a
    partition of `n`).

    INPUT:

    - ``alpha`` -- integer partition (can be encoded as a list)

    OUTPUT:

    The factorial of the size of ``alpha``, divided by the number of
    standard tableaux of shape ``alpha``. Equivalently, the product
    of all hook lengths of ``alpha``.

    EXAMPLES::

        sage: from sage.combinat.symmetric_group_algebra import kappa
        sage: kappa(Partition([2,1]))
        3
        sage: kappa([2,1])
        3
    """
    try:
        n = alpha.size()
    except AttributeError:
        n = sum(alpha)
    return factorial(n) / StandardTableaux(alpha).cardinality()


def a(tableau, star=0, base_ring=QQ):
    r"""
    The row projection operator corresponding to the Young tableau
    ``tableau`` (which is supposed to contain every integer from
    `1` to its size precisely once, but may and may not be standard).

    This is the sum (in the group algebra of the relevant symmetric
    group over `\QQ`) of all the permutations which preserve
    the rows of ``tableau``. It is called `a_{\text{tableau}}` in
    [EGHLSVY]_, Section 4.2.

    INPUT:

    - ``tableau`` -- Young tableau which contains every integer
      from `1` to its size precisely once

    - ``star`` -- nonnegative integer (default: `0`); when this
      optional variable is set, the method computes not the row
      projection operator of ``tableau``, but the row projection
      operator of the restriction of ``tableau`` to the entries
      ``1, 2, ..., tableau.size() - star`` instead.

    - ``base_ring`` -- commutative ring (default: ``QQ``); when this
      optional variable is set, the row projection operator is
      computed over a user-determined base ring instead of `\QQ`.
      (Note that symmetric group algebras currently don't preserve
      coercion, so e. g. a symmetric group algebra over `\ZZ`
      does not coerce into the corresponding one over `\QQ`; so
      convert manually or choose your base rings wisely!)

    EXAMPLES::

        sage: from sage.combinat.symmetric_group_algebra import a
        sage: a([[1,2]])
        [1, 2] + [2, 1]
        sage: a([[1],[2]])
        [1, 2]
        sage: a([])
        []
        sage: a([[1, 5], [2, 3], [4]])
        [1, 2, 3, 4, 5] + [1, 3, 2, 4, 5] + [5, 2, 3, 4, 1] + [5, 3, 2, 4, 1]
        sage: a([[1,4], [2,3]], base_ring=ZZ)
        [1, 2, 3, 4] + [1, 3, 2, 4] + [4, 2, 3, 1] + [4, 3, 2, 1]

    The same with a skew tableau::

        sage: a([[None,1,4], [2,3]], base_ring=ZZ)
        [1, 2, 3, 4] + [1, 3, 2, 4] + [4, 2, 3, 1] + [4, 3, 2, 1]
    """
    t = SkewTableau(tableau)
    if star:
        t = t.restrict(t.size() - star)

    rs = t.row_stabilizer().list()
    n = t.size()

    sgalg = SymmetricGroupAlgebra(base_ring, n)
    one = base_ring.one()
    P = Permutation

    # Ugly hack for the case of an empty tableau, due to the
    # annoyance of Permutation(Tableau([]).row_stabilizer()[0])
    # being [1] rather than [] (which seems to have its origins in
    # permutation group code).
    # TODO: Fix this.
    if n <= 1:
        return sgalg.one()

    rd = dict((P(h), one) for h in rs)
    return sgalg._from_dict(rd)


def b(tableau, star=0, base_ring=QQ):
    r"""
    The column projection operator corresponding to the Young tableau
    ``tableau`` (which is supposed to contain every integer from
    `1` to its size precisely once, but may and may not be standard).

    This is the signed sum (in the group algebra of the relevant
    symmetric group over `\QQ`) of all the permutations which
    preserve the column of ``tableau`` (where the signs are the usual
    signs of the permutations). It is called `b_{\text{tableau}}` in
    [EGHLSVY]_, Section 4.2.

    INPUT:

    - ``tableau`` -- Young tableau which contains every integer
      from `1` to its size precisely once

    - ``star`` -- nonnegative integer (default: `0`). When this
      optional variable is set, the method computes not the column
      projection operator of ``tableau``, but the column projection
      operator of the restriction of ``tableau`` to the entries
      ``1, 2, ..., tableau.size() - star`` instead.

    - ``base_ring`` -- commutative ring (default: ``QQ``). When this
      optional variable is set, the column projection operator is
      computed over a user-determined base ring instead of `\QQ`.
      (Note that symmetric group algebras currently don't preserve
      coercion, so e. g. a symmetric group algebra over `\ZZ`
      does not coerce into the corresponding one over `\QQ`; so
      convert manually or choose your base rings wisely!)

    EXAMPLES::

        sage: from sage.combinat.symmetric_group_algebra import b
        sage: b([[1,2]])
        [1, 2]
        sage: b([[1],[2]])
        [1, 2] - [2, 1]
        sage: b([])
        []
        sage: b([[1, 2, 4], [5, 3]])
        [1, 2, 3, 4, 5] - [1, 3, 2, 4, 5] - [5, 2, 3, 4, 1] + [5, 3, 2, 4, 1]
        sage: b([[1, 4], [2, 3]], base_ring=ZZ)
        [1, 2, 3, 4] - [1, 2, 4, 3] - [2, 1, 3, 4] + [2, 1, 4, 3]
        sage: b([[1, 4], [2, 3]], base_ring=Integers(5))
        [1, 2, 3, 4] + 4*[1, 2, 4, 3] + 4*[2, 1, 3, 4] + [2, 1, 4, 3]

    The same with a skew tableau::

        sage: b([[None, 2, 4], [1, 3], [5]])
        [1, 2, 3, 4, 5] - [1, 3, 2, 4, 5] - [5, 2, 3, 4, 1] + [5, 3, 2, 4, 1]

    With the ``l2r`` setting for multiplication, the unnormalized
    Young symmetrizer ``e(tableau)`` should be the product
    ``b(tableau) * a(tableau)`` for every ``tableau``. Let us check
    this on the standard tableaux of size 5::

        sage: from sage.combinat.symmetric_group_algebra import a, b, e
        sage: all( e(t) == b(t) * a(t) for t in StandardTableaux(5) )
        True
    """
    t = SkewTableau(tableau)
    if star:
        t = t.restrict(t.size() - star)

    cs = t.column_stabilizer().list()
    n = t.size()

    sgalg = SymmetricGroupAlgebra(base_ring, n)
    one = base_ring.one()
    P = Permutation

    # Ugly hack for the case of an empty tableau, due to the
    # annoyance of Permutation(Tableau([]).row_stabilizer()[0])
    # being [1] rather than [] (which seems to have its origins in
    # permutation group code).
    # TODO: Fix this.
    if n <= 1:
        return sgalg.one()

    cd = dict((P(v), v.sign() * one) for v in cs)
    return sgalg._from_dict(cd)


e_cache = {}


def e(tableau, star=0):
    r"""
    The unnormalized Young projection operator corresponding to
    the Young tableau ``tableau`` (which is supposed to contain
    every integer from `1` to its size precisely once, but may
    and may not be standard).

    If `n` is a nonnegative integer, and `T` is a Young tableau
    containing every integer from `1` to `n` exactly once, then
    the unnormalized Young projection operator `e(T)` is defined by

    .. MATH::

        e(T) = a(T) b(T) \in \QQ S_n,

    where `a(T) \in \QQ S_n` is the sum of all permutations in `S_n`
    which fix the rows of `T` (as sets), and `b(T) \in \QQ S_n` is the
    signed sum of all permutations in `S_n` which fix the columns of
    `T` (as sets). Here, "signed" means that each permutation is
    multiplied with its sign; and the product on the group `S_n` is
    defined in such a way that `(pq)(i) = p(q(i))` for any
    permutations `p` and `q` and any `1 \leq i \leq n`.

    Note that the definition of `e(T)` is not uniform across
    literature. Others define it as `b(T) a(T)` instead, or include
    certain scalar factors (we do not, whence "unnormalized").

    EXAMPLES::

        sage: from sage.combinat.symmetric_group_algebra import e
        sage: e([[1,2]])
        [1, 2] + [2, 1]
        sage: e([[1],[2]])
        [1, 2] - [2, 1]
        sage: e([])
        []

    There are differing conventions for the order of the symmetrizers
    and antisymmetrizers.  This example illustrates our conventions::

        sage: e([[1,2],[3]])
        [1, 2, 3] + [2, 1, 3] - [3, 1, 2] - [3, 2, 1]

    To obtain the product `b(T) a(T)`, one has to take the antipode
    of this::

        sage: QS3 = parent(e([[1,2],[3]]))
        sage: QS3.antipode(e([[1,2],[3]]))
        [1, 2, 3] + [2, 1, 3] - [2, 3, 1] - [3, 2, 1]

    And here is an example for a skew tableau::

        sage: e([[None, 2, 1], [4, 3]])
        [1, 2, 3, 4] + [1, 2, 4, 3] - [1, 3, 2, 4] - [1, 4, 2, 3]
         + [2, 1, 3, 4] + [2, 1, 4, 3] - [2, 3, 1, 4] - [2, 4, 1, 3]

    .. SEEALSO::

        :func:`e_hat`
    """
    # TODO:
    # The current method only computes the e's over QQ. There should be
    # a way to compute them over other base rings as well. Be careful
    # with the cache.

    t = SkewTableau(tableau)
    if star:
        t = t.restrict(t.size() - star)

    if t in e_cache:
        res = e_cache[t]
    else:
        rs = t.row_stabilizer().list()
        cs = t.column_stabilizer().list()
        n = t.size()

        QSn = SymmetricGroupAlgebra(QQ, n)
        one = QQ.one()
        P = Permutation

        rd = {P(h): one for h in rs}
        sym = QSn._from_dict(rd)

        cd = {P(v): QQ(v.sign()) for v in cs}
        antisym = QSn._from_dict(cd)

        res = QSn.right_action_product(antisym, sym)

        # Ugly hack for the case of an empty tableau, due to the
        # annoyance of Permutation(Tableau([]).row_stabilizer()[0])
        # being [1] rather than [] (which seems to have its origins in
        # permutation group code).
        # TODO: Fix this.
        if n <= 1:
            return QSn.one()

        e_cache[t] = res

    return res


ehat_cache = {}


def e_hat(tab, star=0):
    r"""
    The Young projection operator corresponding to the Young tableau
    ``tab`` (which is supposed to contain every integer from `1` to
    its size precisely once, but may and may not be standard). This
    is an idempotent in the rational group algebra.

    If `n` is a nonnegative integer, and `T` is a Young tableau
    containing every integer from `1` to `n` exactly once, then
    the Young projection operator `\widehat{e}(T)` is defined by

    .. MATH::

        \widehat{e}(T) = \frac{1}{\kappa_\lambda} a(T) b(T) \in \QQ S_n,

    where `\lambda` is the shape of `T`, where `\kappa_\lambda` is
    `n!` divided by the number of standard tableaux of shape
    `\lambda`, where `a(T) \in \QQ S_n` is the sum of all
    permutations in `S_n` which fix the rows of `T` (as sets), and
    where `b(T) \in \QQ S_n` is the signed sum of all permutations
    in `S_n` which fix the columns of `T` (as sets). Here, "signed"
    means that each permutation is multiplied with its sign; and
    the product on the group `S_n` is defined in such a way that
    `(pq)(i) = p(q(i))` for any permutations `p` and `q` and any
    `1 \leq i \leq n`.

    Note that the definition of `\widehat{e}(T)` is not uniform
    across literature. Others define it as
    `\frac{1}{\kappa_\lambda} b(T) a(T)` instead.

    EXAMPLES::

        sage: from sage.combinat.symmetric_group_algebra import e_hat
        sage: e_hat([[1,2,3]])
        1/6*[1, 2, 3] + 1/6*[1, 3, 2] + 1/6*[2, 1, 3] + 1/6*[2, 3, 1] + 1/6*[3, 1, 2] + 1/6*[3, 2, 1]
        sage: e_hat([[1],[2]])
        1/2*[1, 2] - 1/2*[2, 1]

    There are differing conventions for the order of the symmetrizers
    and antisymmetrizers.  This example illustrates our conventions::

        sage: e_hat([[1,2],[3]])
        1/3*[1, 2, 3] + 1/3*[2, 1, 3] - 1/3*[3, 1, 2] - 1/3*[3, 2, 1]

    .. SEEALSO::

        :func:`e`
    """
    t = SkewTableau(tab)
    # This is for consistency's sake. This method is NOT meant
    # to be applied to skew tableaux, since the meaning of
    # \kappa is unclear in that case.
    if star:
        t = t.restrict(t.size() - star)
    if t in ehat_cache:
        res = ehat_cache[t]
    else:
        res = (1 / kappa(t.shape())) * e(t)
    return res


e_ik_cache = {}


def e_ik(itab, ktab, star=0):
    """
    EXAMPLES::

        sage: from sage.combinat.symmetric_group_algebra import e_ik
        sage: e_ik([[1,2,3]], [[1,2,3]])
        [1, 2, 3] + [1, 3, 2] + [2, 1, 3] + [2, 3, 1] + [3, 1, 2] + [3, 2, 1]
        sage: e_ik([[1,2,3]], [[1,2,3]], star=1)
        [1, 2] + [2, 1]
    """
    it = SkewTableau(itab)
    kt = SkewTableau(ktab)
    if star:
        it = it.restrict(it.size() - star)
        kt = kt.restrict(kt.size() - star)

    if it.shape() != kt.shape():
        raise ValueError("the two tableaux must be of the same shape")

    if kt == it:
        return e(it)
    if (it, kt) in e_ik_cache:
        return e_ik_cache[(it, kt)]

    pi = pi_ik(it, kt)
    QSn = pi.parent()
    res = QSn.right_action_product(e(it), pi)
    e_ik_cache[(it, kt)] = res
    return res


def seminormal_test(n):
    """
    Run a variety of tests to verify that the construction of the
    seminormal basis works as desired. The numbers appearing are
    results in James and Kerber's 'Representation Theory of the
    Symmetric Group' [JK1981]_.

    EXAMPLES::

        sage: from sage.combinat.symmetric_group_algebra import seminormal_test
        sage: seminormal_test(3)
        True
    """
    for part in Partitions_n(n):
        for tab in StandardTableaux(part):
            # Theorem 3.1.10
            if not e(tab) * (1 / kappa(part)) - e_hat(tab) == 0:
                raise ValueError("3.1.10 - %s" % tab)

            # Lemma 3.2.12 (ii)
            value = e(tab) * epsilon(tab, 1) * e(tab) - e(tab) * kappa(part)
            if value != 0:
                print(value)
                raise ValueError("3.2.12.2 - %s" % tab)

            for tab2 in StandardTableaux(part):
                # 3.2.8 (i)
                if e_ik(tab, tab2) - e(tab) * pi_ik(tab, tab2) * e(tab2) * (1 / kappa(part)) != 0:
                    raise ValueError("3.2.8.1 - %s, %s" % (tab, tab2))

                # 3.2.8 (ii)
                if e(tab) * e_ik(tab, tab2) - e_ik(tab, tab2) * kappa(part) != 0:
                    raise ValueError("3.2.8.2 - %s, %s" % (tab, tab2))

                if tab == tab2:
                    continue

                if tab.last_letter_lequal(tab2):
                    # Lemma 3.1.20
                    if e(tab2) * e(tab) != 0:
                        raise ValueError("3.1.20 - %s, %s" % (tab, tab2))
                    if e_hat(tab2) * e_hat(tab) != 0:
                        raise ValueError("3.1.20 - %s, %s" % (tab, tab2))
    return True


#######################


class SGACellularBasis(CellularBasis):
    r"""
    A cellular basis of the symmetric group algebra.
    """
    def __init__(self, SGA):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(GF(3), 3)
            sage: M = SGA.murphy_basis()
            sage: TestSuite(M).run()
            sage: KL = SGA.kazhdan_lusztig_cellular_basis()
            sage: TestSuite(KL).run()
        """
        CellularBasis.__init__(self, SGA, self._to_sga)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(GF(3), 4)
            sage: SGA.murphy_basis()
            Murphy basis of Symmetric group algebra of order 4 over Finite Field of size 3
            sage: SGA.kazhdan_lusztig_cellular_basis()
            Kazhdan-Lusztig basis of Symmetric group algebra of order 4 over Finite Field of size 3
        """
        return self._name + " basis of {}".format(self._algebra)

    @cached_method
    def one_basis(self):
        r"""
        Return the index of the basis element for the multiplicative identity.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(GF(3), 4)
            sage: M = SGA.murphy_basis()
            sage: M.one_basis()
            ([4], [[1, 2, 3, 4]], [[1, 2, 3, 4]])
        """
        la = _Partitions([self._algebra.n])
        col = la.standard_tableaux()[0]
        return (la, col, col)

    @cached_method
    def one(self):
        r"""
        Return the element `1` in ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(GF(3), 4)
            sage: M = SGA.murphy_basis()
            sage: M.one()
            C([4], [[1, 2, 3, 4]], [[1, 2, 3, 4]])
        """
        return self.monomial(self.one_basis())


class MurphyBasis(SGACellularBasis):
    r"""
    The Murphy basis of a symmetric group algebra.

    Let `R` be a commutative ring, and let `A = R[S_n]` denote the group
    algebra (over `R`) of `S_n`. The *Murphy basis* is the basis of `A`
    defined as follows. Let `S, T` be standard tableaux of shape `\lambda`.
    Define `T^{\lambda}` as the standard tableau of shape `\lambda` with
    the first row filled with `1, \ldots, \lambda_1`, the second row
    `\lambda_1+1, \ldots, \lambda_1+\lambda_2`, and so on. Let `d(S)` be
    the unique permutation such that `S = T^{\lambda} d(S)` under the natural
    action. Then the Murphy basis element indexed by `S` and `T` is

    .. MATH::

        M_{S'T'} = d(S)^{-1} R_{\lambda} d(T),

    where `S'` denotes the conjugate tableau.
    The Murphy basis is a :class:`cellular basis
    <sage.algebras.cellular_basis.CellularBasis>` of `A`.

    EXAMPLES::

        sage: SGA = SymmetricGroupAlgebra(GF(3), 5)
        sage: M = SGA.murphy_basis()
        sage: for la in M.simple_module_parameterization():
        ....:     CM = M.cell_module(la)
        ....:     print(la, CM.dimension(), CM.simple_module().dimension())
        [2, 2, 1] 5 4
        [3, 1, 1] 6 6
        [3, 2] 5 1
        [4, 1] 4 4
        [5] 1 1

    REFERENCES:

    - [DJM1998]_
    - [Mathas2004]_
    """
    _name = "Murphy"

    def _to_sga(self, ind):
        r"""
        Return the element in the symmetric group algebra indexed by ``ind``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(GF(3), 3)
            sage: M = SGA.murphy_basis()
            sage: for ind in M.basis().keys():
            ....:     print(ind, M._to_sga(ind))
            ([1, 1, 1], [[1], [2], [3]], [[1], [2], [3]]) [1, 2, 3] + [1, 3, 2] + [2, 1, 3] + [2, 3, 1] + [3, 1, 2] + [3, 2, 1]
            ([2, 1], [[1, 3], [2]], [[1, 3], [2]]) [1, 2, 3] + [2, 1, 3]
            ([2, 1], [[1, 3], [2]], [[1, 2], [3]]) [1, 3, 2] + [3, 1, 2]
            ([2, 1], [[1, 2], [3]], [[1, 3], [2]]) [1, 3, 2] + [2, 3, 1]
            ([2, 1], [[1, 2], [3]], [[1, 2], [3]]) [1, 2, 3] + [3, 2, 1]
            ([3], [[1, 2, 3]], [[1, 2, 3]]) [1, 2, 3]
        """
        return self._algebra.murphy_basis_element(ind[1], ind[2])


class KLCellularBasis(SGACellularBasis):
    """
    The Kazhdan-Lusztig `C'` basis (at `q = 1`) of the symmetric group
    algebra realized as a :class:`cellular basis
    <sage.algebras.cellular_basis.CellularBasis>`

    EXAMPLES::

        sage: SGA = SymmetricGroupAlgebra(GF(3), 5)
        sage: KL = SGA.kazhdan_lusztig_cellular_basis()
        sage: for la in KL.simple_module_parameterization():
        ....:     CM = KL.cell_module(la)
        ....:     print(la, CM.dimension(), CM.simple_module().dimension())
        [2, 2, 1] 5 4
        [3, 1, 1] 6 6
        [3, 2] 5 1
        [4, 1] 4 4
        [5] 1 1
    """
    _name = "Kazhdan-Lusztig"

    def _to_sga(self, ind):
        r"""
        Return the element in the symmetric group algebra indexed by ``ind``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(GF(3), 3)
            sage: KL = SGA.kazhdan_lusztig_cellular_basis()
            sage: for ind in KL.basis().keys():
            ....:     print(ind, KL._to_sga(ind))
            ([1, 1, 1], [[1], [2], [3]], [[1], [2], [3]]) [1, 2, 3] + [1, 3, 2] + [2, 1, 3] + [2, 3, 1] + [3, 1, 2] + [3, 2, 1]
            ([2, 1], [[1, 3], [2]], [[1, 3], [2]]) [1, 2, 3] + [2, 1, 3]
            ([2, 1], [[1, 3], [2]], [[1, 2], [3]]) [1, 2, 3] + [1, 3, 2] + [2, 1, 3] + [3, 1, 2]
            ([2, 1], [[1, 2], [3]], [[1, 3], [2]]) [1, 2, 3] + [1, 3, 2] + [2, 1, 3] + [2, 3, 1]
            ([2, 1], [[1, 2], [3]], [[1, 2], [3]]) [1, 2, 3] + [1, 3, 2]
            ([3], [[1, 2, 3]], [[1, 2, 3]]) [1, 2, 3]
        """
        from sage.combinat.rsk import RSK_inverse
        S = ind[1]
        T = ind[2]
        w = RSK_inverse(T, S, output='permutation')
        return self._algebra.kazhdan_lusztig_basis_element(w)


#######################


def HeckeAlgebraSymmetricGroupT(R, n, q=None):
    r"""
    Return the Hecke algebra of the symmetric group `S_n` on the T-basis
    with quantum parameter ``q`` over the ring `R`.

    If `R` is a commutative ring and `q` is an invertible element of `R`,
    and if `n` is a nonnegative integer, then the Hecke algebra of the
    symmetric group `S_n` over `R` with quantum parameter `q` is defined
    as the algebra generated by the generators `T_1, T_2, \ldots, T_{n-1}`
    with relations

    .. MATH::

        T_i T_{i+1} T_i = T_{i+1} T_i T_{i+1}

    for all `i < n-1` ("braid relations"),

    .. MATH::

        T_i T_j = T_j T_i

    for all `i` and `j` such that `| i-j | > 1` ("locality relations"),
    and

    .. MATH::

        T_i^2 = q + (q-1) T_i

    for all `i` (the "quadratic relations", also known in the form
    `(T_i + 1) (T_i - q) = 0`). (This is only one of several existing
    definitions in literature, not all of which are fully equivalent.
    We are following the conventions of [Go1993]_.) For any permutation
    `w \in S_n`, we can define an element `T_w` of this Hecke algebra by
    setting `T_w = T_{i_1} T_{i_2} \cdots T_{i_k}`, where
    `w = s_{i_1} s_{i_2} \cdots s_{i_k}` is a reduced word for `w`
    (with `s_i` meaning the transposition `(i, i+1)`, and the product of
    permutations being evaluated by first applying `s_{i_k}`, then
    `s_{i_{k-1}}`, etc.). This element is independent of the choice of
    the reduced decomposition, and can be computed in Sage by calling
    ``H[w]`` where ``H`` is the Hecke algebra and ``w`` is the
    permutation.

    The Hecke algebra of the symmetric group `S_n` with quantum parameter
    `q` over `R` can be seen as a deformation of the group algebra
    `R S_n`; indeed, it becomes `R S_n` when `q = 1`.

    .. WARNING::

        The multiplication on the Hecke algebra of the symmetric group
        does *not* follow the global option ``mult`` of the
        :class:`Permutations` class (see
        :meth:`~sage.combinat.permutation.Permutations.options`).
        It is always as defined above. It does not match the default
        option (``mult=l2r``) of the symmetric group algebra!

    EXAMPLES::

        sage: HeckeAlgebraSymmetricGroupT(QQ, 3)
        Hecke algebra of the symmetric group of order 3 on the T basis over Univariate Polynomial Ring in q over Rational Field

    ::

        sage: HeckeAlgebraSymmetricGroupT(QQ, 3, 2)
        Hecke algebra of the symmetric group of order 3 with q=2 on the T basis over Rational Field

    The multiplication on the Hecke algebra follows a different convention
    than the one on the symmetric group algebra does by default::

        sage: H3 = HeckeAlgebraSymmetricGroupT(QQ, 3)
        sage: H3([1,3,2]) * H3([2,1,3])
        T[3, 1, 2]
        sage: S3 = SymmetricGroupAlgebra(QQ, 3)
        sage: S3([1,3,2]) * S3([2,1,3])
        [2, 3, 1]

        sage: TestSuite(H3).run()

    .. NOTE::

        :class:`~sage.algebras.iwahori_hecke_algebra.IwahoriHeckeAlgebra` gives
        a different implementation of the Iwahori-Hecke algebras of a Coxeter
        system `(W,S)`. This includes the Hecke algebras of the symmetric group
        as special case.
    """
    return HeckeAlgebraSymmetricGroup_t(R, n, q)


class HeckeAlgebraSymmetricGroup_generic(CombinatorialFreeModule):
    def __init__(self, R, n, q=None):
        """
        TESTS::

            sage: HeckeAlgebraSymmetricGroupT(QQ, 3)
            Hecke algebra of the symmetric group of order 3 on the T basis over Univariate Polynomial Ring in q over Rational Field

        ::

            sage: HeckeAlgebraSymmetricGroupT(QQ, 3, q=1)
            Hecke algebra of the symmetric group of order 3 with q=1 on the T basis over Rational Field
        """
        self.n = n
        self._indices = Permutations(n)
        self._name = "Hecke algebra of the symmetric group of order {}".format(n)
        if q is None:
            q = PolynomialRing(R, 'q').gen()
            R = q.parent()
        else:
            if q not in R:
                raise ValueError("q must be in R (= {})".format(R))
            self._name += " with q={}".format(q)

        self._q = q

        CombinatorialFreeModule.__init__(self, R, self._indices,
                                         category=AlgebrasWithBasis(R),
                                         prefix="")

    _repr_option_bracket = False

    @cached_method
    def one_basis(self):
        """
        Return the identity permutation.

        EXAMPLES::

            sage: HeckeAlgebraSymmetricGroupT(QQ, 3).one()  # indirect doctest
            T[1, 2, 3]
        """
        return self._indices.one()

    def q(self):
        """
        Return the variable or parameter `q`.

        EXAMPLES::

            sage: HeckeAlgebraSymmetricGroupT(QQ, 3).q()
            q
            sage: HeckeAlgebraSymmetricGroupT(QQ, 3, 2).q()
            2
        """
        return self._q

    def _element_constructor_(self, x):
        """
        EXAMPLES::

            sage: H3 = HeckeAlgebraSymmetricGroupT(QQ, 3)
            sage: H3([2,1])  # indirect doc test
            T[2, 1, 3]
            sage: H3( Permutations(3).an_element() )
            T[3, 1, 2]
            sage: H3( SymmetricGroup(3).an_element() )
            [1, 3, 2]
            sage: H3( [2, 1, 3, 4, 5, 6] )
            T[2, 1, 3]
        """
        ###################################################
        # Coerce permutations of size smaller that self.n #
        ###################################################
        if not x:
            return self.one()
        if x in Permutations():
            if len(x) < self.n:
                return self.monomial(self._indices(
                            list(x) + list(range(len(x) + 1, self.n + 1))
                        ))
            if all(x[i] == i+1 for i in range(self.n, len(x))):
                return self.monomial(self._indices(x[:self.n]))

        return self._indices(x)


class HeckeAlgebraSymmetricGroup_t(HeckeAlgebraSymmetricGroup_generic):

    def __init__(self, R, n, q=None):
        """
        TESTS::

            sage: H3 = HeckeAlgebraSymmetricGroupT(QQ, 3)
            sage: H3 == loads(dumps(H3))
            True
        """
        HeckeAlgebraSymmetricGroup_generic.__init__(self, R, n, q)
        self._name += " on the T basis"
        self.print_options(prefix='T')

    def t_action_on_basis(self, perm, i):
        r"""
        Return the product `T_i \cdot T_{perm}`, where ``perm`` is a
        permutation in the symmetric group `S_n`.

        EXAMPLES::

            sage: H3 = HeckeAlgebraSymmetricGroupT(QQ, 3)
            sage: H3.t_action_on_basis(Permutation([2,1,3]), 1)
            q*T[1, 2, 3] + (q-1)*T[2, 1, 3]
            sage: H3.t_action_on_basis(Permutation([1,2,3]), 1)
            T[2, 1, 3]
            sage: H3 = HeckeAlgebraSymmetricGroupT(QQ, 3, 1)
            sage: H3.t_action_on_basis(Permutation([2,1,3]), 1)
            T[1, 2, 3]
            sage: H3.t_action_on_basis(Permutation([1,3,2]), 2)
            T[1, 2, 3]
        """
        if i not in range(1, self.n):
            raise ValueError("i (= %(i)d) must be between 1 and n (= %(n)d)" % {'i': i, 'n': self.n})

        t_i = Permutation((i, i + 1))
        perm_i = t_i.right_action_product(perm)
        # This used to be perm_i = t_i * perm. I have changed it to
        # perm_i = t_i.right_action_product(perm) because it would
        # otherwise cause TestSuite(H3) to fail when
        # Permutations.options(mult) would be set to "r2l".
        # -- Darij, 19 Nov 2013

        if perm[i - 1] < perm[i]:
            return self.monomial(self._indices(perm_i))
        else:
            # Ti^2 = (q - q^(-1))*Ti - q1*q2
            q = self.q()
            z_elt = {perm_i: q, perm: q - 1}
            return self._from_dict(z_elt)

    def t_action(self, a, i):
        r"""
        Return the product `T_i \cdot a`.

        EXAMPLES::

            sage: H3 = HeckeAlgebraSymmetricGroupT(QQ, 3)
            sage: a = H3([2,1,3])+2*H3([1,2,3])
            sage: H3.t_action(a, 1)
            q*T[1, 2, 3] + (q+1)*T[2, 1, 3]
            sage: H3.t(1)*a
            q*T[1, 2, 3] + (q+1)*T[2, 1, 3]
        """
        def t_i(x):
            return self.t_action_on_basis(x, i)
        return self._apply_module_endomorphism(a, t_i)

    def product_on_basis(self, perm1, perm2):
        """
        EXAMPLES::

            sage: H3 = HeckeAlgebraSymmetricGroupT(QQ, 3, 1)
            sage: a = H3([2,1,3])+2*H3([1,2,3])-H3([3,2,1])
            sage: a^2 #indirect doctest
            6*T[1, 2, 3] + 4*T[2, 1, 3] - T[2, 3, 1] - T[3, 1, 2] - 4*T[3, 2, 1]

        ::

            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: a = QS3([2,1,3])+2*QS3([1,2,3])-QS3([3,2,1])
            sage: a^2
            6*[1, 2, 3] + 4*[2, 1, 3] - [2, 3, 1] - [3, 1, 2] - 4*[3, 2, 1]
        """
        res = self(perm1)
        for i in perm2.reduced_word():
            res = self.t_action(res, i)
        return res

    def t(self, i):
        """
        Return the element `T_i` of the Hecke algebra ``self``.

        EXAMPLES::

            sage: H3 = HeckeAlgebraSymmetricGroupT(QQ,3)
            sage: H3.t(1)
            T[2, 1, 3]
            sage: H3.t(2)
            T[1, 3, 2]
            sage: H3.t(0)
            Traceback (most recent call last):
            ...
            ValueError: i (= 0) must be between 1 and n-1 (= 2)
        """
        if i not in range(1, self.n):
            raise ValueError("i (= %(i)d) must be between 1 and n-1 (= %(nm)d)" % {'i': i, 'nm': self.n - 1})

        P = self.basis().keys()
        return self.monomial(P(list(range(1, i)) + [i + 1, i] +
                               list(range(i + 2, self.n + 1))))
        # The permutation here is simply the transposition (i, i+1).

    def algebra_generators(self):
        """
        Return the generators of the algebra.

        EXAMPLES::

            sage: HeckeAlgebraSymmetricGroupT(QQ,3).algebra_generators()
            [T[2, 1, 3], T[1, 3, 2]]
        """
        return [self.t(i) for i in range(1, self.n)]

    def jucys_murphy(self, k):
        r"""
        Return the Jucys-Murphy element `J_k` of the Hecke algebra.

        These Jucys-Murphy elements are defined by

        .. MATH::

            J_k = (T_{k-1} T_{k-2} \cdots T_1) (T_1 T_2 \cdots T_{k-1}).

        More explicitly,

        .. MATH::

            J_k = q^{k-1} + \sum_{l=1}^{k-1} (q^l - q^{l-1}) T_{(l, k)}.

        For generic `q`, the `J_k` generate a maximal commutative
        sub-algebra of the Hecke algebra.

        .. WARNING::

            The specialization `q = 1` does *not* map these elements
            `J_k` to the Young-Jucys-Murphy elements of the group
            algebra `R S_n`. (Instead, it maps the "reduced"
            Jucys-Murphy elements `(J_k - q^{k-1}) / (q - 1)` to the
            Young-Jucys-Murphy elements of `R S_n`.)

        EXAMPLES::

            sage: H3 = HeckeAlgebraSymmetricGroupT(QQ,3)
            sage: j2 = H3.jucys_murphy(2); j2
            q*T[1, 2, 3] + (q-1)*T[2, 1, 3]
            sage: j3 = H3.jucys_murphy(3); j3
            q^2*T[1, 2, 3] + (q^2-q)*T[1, 3, 2] + (q-1)*T[3, 2, 1]
            sage: j2*j3 == j3*j2
            True
            sage: j0 = H3.jucys_murphy(1); j0 == H3.one()
            True
            sage: H3.jucys_murphy(0)
            Traceback (most recent call last):
            ...
            ValueError: k (= 0) must be between 1 and n (= 3)
        """
        if k not in range(2, self.n + 1):
            if k == 1:
                return self.one()
            raise ValueError("k (= %(k)d) must be between 1 and n (= %(n)d)" % {'k': k, 'n': self.n})

        q = self.q()
        P = self._indices
        v = self.sum_of_terms(((P(list(range(1, l)) + [k] +
                                  list(range(l + 1, k)) + [l]),
                                q**l - q**(l - 1))
                               for l in range(1, k)),
                              distinct=True)
        v += q**(k - 1) * self.one()
        return v

        # old algorithm:
        # left = 1
        # right = 1
        # for j in range(1, k):
        #    left *= self.t(k-j)
        #    right *= self.t(j)
        # return left*right


# For unpickling backward compatibility (Sage <= 4.1)
register_unpickle_override('sage.combinat.symmetric_group_algebra',
                           'HeckeAlgebraSymmetricGroupElement_t',
                           CombinatorialFreeModule.Element)
register_unpickle_override('sage.combinat.symmetric_group_algebra',
                           'SymmetricGroupAlgebraElement_n',
                           CombinatorialFreeModule.Element)
